#include "coupled_navier_stokes.h"

// ********************************************************** //
//
// This file includes all the functions not specific to 
// 3D or 1D domains.
//
// setup_3d_mesh
// setup_3d_system
// prerefine_3d_mesh
// calculate_3d_boundary_values
// write_3d_solution
// solve_3d_system_iteration
// adaptively_refine
// plot_error
//
// *********************************************************** //

 

// ******************** SETUP THE 3D MESH ************************** //
// # 0 is cylindrical pipe, 1 is single bifurcating pipe, 2 is cuboid, 3 is closed cuboid
void NavierStokesCoupled::setup_3d_mesh(EquationSystems* _es,Mesh& _mesh)
{
	std::cout << std::endl;
	std::cout << "Setting up 2D/3D mesh." << std::endl;
	
	// always reread in the mesh. forget about mesh refinement for now.
	// if we are doing a 2D simulation we need to change mesh back to having dimension 2
	if(!_es->parameters.get<bool> ("threed"))
		_mesh.set_mesh_dimension(2);

	// both the cylindrical pipe and bifurcating pipe require input mesh as well as bifurcating pipe
	if(_es->parameters.get<unsigned int> ("geometry_type") == 0 || _es->parameters.get<unsigned int> ("geometry_type") == 1 
			|| _es->parameters.get<unsigned int> ("geometry_type") == 4  || _es->parameters.get<bool> ("expanding_pipe"))
	{
		std::cout << "Reading mesh." << std::endl;
		_mesh.read(_es->parameters.get<std::string> ("mesh_file"));

		
	
		//scale the mesh to make SI units
		if(!_es->parameters.get<bool>("linear_shape_functions"))
				_mesh.all_second_order();

		std::cout << "Refining mesh." << std::endl;
		mesh_refinement.uniformly_refine (_es->parameters.get<unsigned int> ("no_refinement"));

		std::cout << "Scaling mesh." << std::endl;
		double mesh_input_scaling_3d = _es->parameters.get<double>("mesh_input_scaling_3d");
		MeshTools::Modification::scale(_mesh,mesh_input_scaling_3d,mesh_input_scaling_3d,mesh_input_scaling_3d);

	}
	// the cubioids can be auto generated
	else if(es->parameters.get<unsigned int> ("geometry_type") == 2 || es->parameters.get<unsigned int> ("geometry_type") == 3
					|| _es->parameters.get<unsigned int> ("geometry_type") == 5)
	{

		if(threed){
			if(_es->parameters.get<bool>("quads"))
				MeshTools::Generation::build_cube (_mesh,_es->parameters.get<unsigned int>("cube_length_N"),_es->parameters.get<unsigned int>("cube_width_N"),_es->parameters.get<unsigned int>("cube_height_N"),0., _es->parameters.get<double>("cube_length"),0., _es->parameters.get<double>("cube_width"),0., _es->parameters.get<double>("cube_height"),HEX8);
			else
				MeshTools::Generation::build_cube (_mesh,_es->parameters.get<unsigned int>("cube_length_N"),_es->parameters.get<unsigned int>("cube_width_N"),_es->parameters.get<unsigned int>("cube_height_N"),0., _es->parameters.get<double>("cube_length"),0., _es->parameters.get<double>("cube_width"),0., _es->parameters.get<double>("cube_height"),TET4);					
		}
		else
		{
			if(_es->parameters.get<bool>("quads"))
				MeshTools::Generation::build_square (_mesh,_es->parameters.get<unsigned int>("cube_length_N"),_es->parameters.get<unsigned int>("cube_width_N"),
																							0.,_es->parameters.get<double>("cube_length"),0., _es->parameters.get<double>("cube_width"),QUAD4);
			else
				MeshTools::Generation::build_square (_mesh,_es->parameters.get<unsigned int>("cube_length_N"),_es->parameters.get<unsigned int>("cube_width_N"),
																							0.,_es->parameters.get<double>("cube_length"),0., _es->parameters.get<double>("cube_width"),TRI3);
		}



		double mesh_input_scaling_3d = _es->parameters.get<double>("mesh_input_scaling_3d");
		MeshTools::Modification::scale(_mesh,mesh_input_scaling_3d,mesh_input_scaling_3d,mesh_input_scaling_3d);

		if(!_es->parameters.get<bool>("linear_shape_functions"))
			_mesh.all_second_order();

		mesh_refinement.uniformly_refine (_es->parameters.get<unsigned int> ("no_refinement"));

	}

	//scale the mesh to make length scale

	double length_scale = _es->parameters.get<double>("length_scale");
	MeshTools::Modification::scale(_mesh,1./length_scale,1./length_scale,1./length_scale);
	

	//********* RENAME THE BOUNDARIES AND SUBDOMAINS FOR ALL DIFF PROBLEMS *****//
  //
	// we always read in the 3d mesh first so we can set all subdomain ids to 0
	// how is everything labelled
	// 3d inflow has boundary id 0
	// interfaces have boundary ids 1-n
	// 3d subdomain has id 0
	// 1d subdomains have ids 1-n corresponding to interfaces
	//

	std::cout << "Setting up boundaries." << std::endl;

	// set node ids from side ids so that surface_boundary can identify shared boundaries.
	_mesh.boundary_info->build_node_list_from_side_list();

	// the pipe geometry or the expanding_pipe geometry
	if(_es->parameters.get<unsigned int> ("geometry_type") == 0 || _es->parameters.get<bool> ("expanding_pipe"))
	{
		MeshTools::Modification::change_boundary_id(_mesh,1010,-1);		//walls

		MeshTools::Modification::change_boundary_id(_mesh,1011,0);			//inflow
		MeshTools::Modification::change_boundary_id(_mesh,1012,1);			//outflows

		MeshTools::Modification::change_subdomain_id(_mesh,1020,0);		//volume

		if(es->parameters.get<unsigned int> ("old_geometry") == 1)
		{
			MeshTools::Modification::change_boundary_id(_mesh,0,500);			//inflow
			MeshTools::Modification::change_boundary_id(_mesh,1,0);			//outflows
			MeshTools::Modification::change_boundary_id(_mesh,500,1);			//outflows
		}
		
		// minus 1 because of the wall bc
		boundary_ids.resize(_mesh.boundary_info->n_boundary_ids() - 1);	//hmmmm not sure why but hey


		for(unsigned int i=0; i < boundary_ids.size(); i++)
			boundary_ids[i] = i;
		
		//set surfaceboundary stuff - for all boundaries (we will at least some of the info at output time)
		//inflow_surface_boundary_object.init(_mesh,0);

		for(unsigned int i=0; i < boundary_ids.size(); i++)
		{

			SurfaceBoundary* surface_boundary = new SurfaceBoundary(*_es);
			surface_boundaries.push_back(surface_boundary);
			if(!es->parameters.get<bool> ("expanding_pipe"))
				surface_boundaries[i]->init(_mesh,i,false);
			else
				surface_boundaries[i]->init(_mesh,i,true);

			std::cout << "surface " << i << ":" << std::endl;
			std::cout << "\t normal is " << surface_boundaries[i]->get_normal() << std::endl;
			std::cout << "\t centroid is " << surface_boundaries[i]->get_centroid() << std::endl;
			std::cout << "\t area is " << surface_boundaries[i]->get_area() << std::endl;
			std::cout << "\t unit parabola integral is " << surface_boundaries[i]->get_unit_parabola_integral() << std::endl;
			std::cout << std::endl;
		}

		//hard coded
		// not important to preset these 3d values because they will be reset 
		// after the first iteration
		pressure_values_3d.push_back(1.0);	//inflow
		pressure_values_3d.push_back(0.0);	//outflow

		// post processed values
		flux_values_3d.push_back(1.0);	//inflow
		flux_values_3d.push_back(0.0);	//outflow


		previous_flux_values_3d = flux_values_3d;
		previous_previous_flux_values_3d = flux_values_3d;

		input_pressure_values_3d.push_back(_es->parameters.get<double>("parent_pressure_mag"));
		input_pressure_values_3d.push_back(_es->parameters.get<double>("daughter_1_pressure_mag"));
	}
	else if(es->parameters.get<unsigned int> ("geometry_type") == 1 || es->parameters.get<unsigned int> ("geometry_type") == 4)	// bifurcating pipe
	{

		std::set<boundary_id_type> bdyids = _mesh.boundary_info->get_boundary_ids();
		std::set<subdomain_id_type> subids;
		_mesh.subdomain_ids(subids);
		std::cout << "hye there" << std::endl;
		std::cout << "bdyids.size = " << bdyids.size() << std::endl;

		// we wanna do this relabelling automagically, and yes sets are ordered in a particular way
		if(true)//_es->parameters.get<std::string> ("mesh_file").compare("meshes/coarse_mesh_truncated_labelled.msh") == 0)
		{
			int count = 0;
			for (std::set<boundary_id_type>::iterator it=bdyids.begin(); it!=bdyids.end(); ++it)
			{
				if(es->parameters.get<int> ("inflow_bdy_id") > 0)
				{
					if(count == 0)
					{
						MeshTools::Modification::change_boundary_id(_mesh,*it,-1);		//walls
						std::cout << "old bdy id = " << *it << std::endl;
						std::cout << "new bdy id = " << -1 << std::endl;
						count++;

					}
					else if(*it == es->parameters.get<int> ("inflow_bdy_id"))
					{
						MeshTools::Modification::change_boundary_id(_mesh,*it,0);		//inflow
						std::cout << "old bdy id = " << *it << std::endl;
						std::cout << "new bdy id = " << 0 << std::endl;					
					}
					else
					{
						MeshTools::Modification::change_boundary_id(_mesh,*it,count);		//outlets
						std::cout << "old bdy id = " << *it << std::endl;
						std::cout << "new bdy id = " << count << std::endl;		
						count++;			
					}
				}
				else
				{
					MeshTools::Modification::change_boundary_id(_mesh,*it,count - 1);		//walls
					std::cout << "old bdy id = " << *it << std::endl;
					std::cout << "new bdy id = " << count << std::endl;
					count++;
				
				}
			}

			/*
			MeshTools::Modification::change_boundary_id(_mesh,19,0);			//inflow
			MeshTools::Modification::change_boundary_id(_mesh,20,1);			//outflow
			MeshTools::Modification::change_boundary_id(_mesh,21,2);			//outflow
			MeshTools::Modification::change_boundary_id(_mesh,22,3);			//outflow
			MeshTools::Modification::change_boundary_id(_mesh,23,4);			//outflow
			MeshTools::Modification::change_boundary_id(_mesh,24,5);			//outflow

			MeshTools::Modification::change_subdomain_id(_mesh,3,0);		//volume
			*/

			// should only be one i hope...
			MeshTools::Modification::change_subdomain_id(_mesh,*subids.begin(),0);		//volume

			if(es->parameters.get<unsigned int> ("old_geometry") == 1)
			{
				MeshTools::Modification::change_boundary_id(_mesh,0,500);			//inflow
				MeshTools::Modification::change_boundary_id(_mesh,1,0);			//outflows
				MeshTools::Modification::change_boundary_id(_mesh,500,1);			//outflows
			}
			else if(es->parameters.get<unsigned int> ("old_geometry") == 2)
			{
				MeshTools::Modification::change_boundary_id(_mesh,0,-1);			// walls
				MeshTools::Modification::change_boundary_id(_mesh,1,0);			//outflows
				MeshTools::Modification::change_boundary_id(_mesh,2,1);			//outflows
				MeshTools::Modification::change_boundary_id(_mesh,3,2);			//outflows
			}
		}
		else
		{
			MeshTools::Modification::change_boundary_id(_mesh,1010,-1);		//walls

			MeshTools::Modification::change_boundary_id(_mesh,1011,0);			//inflow
			MeshTools::Modification::change_boundary_id(_mesh,1012,1);			//outflows

			MeshTools::Modification::change_subdomain_id(_mesh,1020,0);		//volume

			if(es->parameters.get<unsigned int> ("old_geometry") == 1)
			{
				MeshTools::Modification::change_boundary_id(_mesh,0,500);			//inflow
				MeshTools::Modification::change_boundary_id(_mesh,1,0);			//outflows
				MeshTools::Modification::change_boundary_id(_mesh,500,1);			//outflows
			}
		}





		_mesh.print_info();

		// so this has the boundary ids of the inflows and outflows, not the walls
		boundary_ids.resize(_mesh.boundary_info->n_boundary_ids() - 1);
		
		// all 3d outlets not inlet
		_es->parameters.set<unsigned int> ("num_1d_trees") = boundary_ids.size() - 1;
		



		
		for(unsigned int i=0; i < boundary_ids.size(); i++)
			boundary_ids[i] = i;

		//set surfaceboundary stuff
		//inflow_surface_boundary_object.init(_mesh,0);

		for(unsigned int i=0; i < boundary_ids.size(); i++)
		{

			// don't wanna do the possible 1d boundary at 1000
			if(boundary_ids[i] < 100)
			{
				SurfaceBoundary* surface_boundary = new SurfaceBoundary(*_es);
				surface_boundaries.push_back(surface_boundary);
				surface_boundaries[i]->init(_mesh,i,0);

				std::cout << "surface " << i << ":" << std::endl;
				std::cout << "\t normal is " << surface_boundaries[i]->get_normal() << std::endl;
				std::cout << "\t centroid is " << surface_boundaries[i]->get_centroid() << std::endl;
				std::cout << "\t area is " << surface_boundaries[i]->get_area() << std::endl;
				std::cout << "\t unit parabola integral is " << surface_boundaries[i]->get_unit_parabola_integral() << std::endl;
				std::cout << "\t max radius is " << surface_boundaries[i]->get_max_radius() << std::endl;
				std::cout << std::endl;
			}
		}


		//hard coded
		// not important to preset these 3d values because they will be reset 
		// after the first iteration
		pressure_values_3d.push_back(1.0);	//inflow
		for(unsigned int i=1; i <  boundary_ids.size(); i++)
			pressure_values_3d.push_back(0.0);	//outflow

		// post processed values
		flux_values_3d.push_back(1.0);	//inflow
		for(unsigned int i=1; i <  boundary_ids.size(); i++)
			flux_values_3d.push_back(0.0);	//outflow


		previous_flux_values_3d = flux_values_3d;
		previous_previous_flux_values_3d = flux_values_3d;

		// inputting the pressure magnitude can come in later
		input_pressure_values_3d.push_back(_es->parameters.get<double>("parent_pressure_mag"));
		for(unsigned int i=1; i <  boundary_ids.size(); i++)
			input_pressure_values_3d.push_back(_es->parameters.get<double>("daughter_1_pressure_mag"));
//		input_pressure_values_3d.push_back(_es->parameters.get<double>("daughter_2_pressure_mag"));

	}
	// the cuboid geometry
	else if(_es->parameters.get<unsigned int> ("geometry_type") == 2)
	{

	
		if(!threed)
		{
			// for some odd reason, applying dirichlet conditions on one of the walls doesn't work?
			MeshTools::Modification::change_boundary_id(_mesh,0,-1);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,1,1);		//inflow	//should be x=0
			MeshTools::Modification::change_boundary_id(_mesh,2,-1);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,3,0);		//outflow	//should be x=L

		}
		else
		{
			MeshTools::Modification::change_boundary_id(_mesh,0,0);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,1,-1);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,2,-1);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,3,-1);		//walls


			MeshTools::Modification::change_boundary_id(_mesh,4,-1);			//inflow
			MeshTools::Modification::change_boundary_id(_mesh,5,1);			//outflows
		}
	
		// subdomain id 0 by default
		//MeshTools::Modification::change_subdomain_id(_mesh,1020,0);		//volume


		// hmmm how many fn boundary ids do we actually have??
		// haven't added the 
			

		boundary_ids.resize(_mesh.boundary_info->n_boundary_ids() - 1);

		for(unsigned int i=0; i < boundary_ids.size(); i++)
			boundary_ids[i] = i;
		

		//set surfaceboundary stuff
		//inflow_surface_boundary_object.init(_mesh,0);

		// unfortunately this does not work for libmesh created geoms because ids not propagated to nodes
		
		for(unsigned int i=0; i < boundary_ids.size(); i++)
		{
			SurfaceBoundary* surface_boundary = new SurfaceBoundary(*_es);
			surface_boundaries.push_back(surface_boundary);
			surface_boundaries[i]->init(_mesh,i,true);

			std::cout << "surface " << i << ":" << std::endl;
			std::cout << "\t normal is " << surface_boundaries[i]->get_normal() << std::endl;
			std::cout << "\t centroid is " << surface_boundaries[i]->get_centroid() << std::endl;
			std::cout << "\t area is " << surface_boundaries[i]->get_area() << std::endl;
			std::cout << "\t unit parabola integral is " << surface_boundaries[i]->get_unit_parabola_integral() << std::endl;
			std::cout << "\t max radius is " << surface_boundaries[i]->get_max_radius() << std::endl;
			std::cout << std::endl;
		}
		


		//hard coded
		// after the first iteration
		pressure_values_3d.push_back(1.0);	//inflow
		pressure_values_3d.push_back(0.0);	//outflow

		// post processed values
		flux_values_3d.push_back(1.0);	//inflow
		flux_values_3d.push_back(0.0);	//outflow

		previous_flux_values_3d = flux_values_3d;
		previous_previous_flux_values_3d = flux_values_3d;
	
		input_pressure_values_3d.push_back(_es->parameters.get<double>("parent_pressure_mag"));
		input_pressure_values_3d.push_back(_es->parameters.get<double>("daughter_1_pressure_mag"));

		// scale so that corresponds to the correct generation
		//MeshTools::Modification::scale(mesh,5.6);
	}
	// the axisymmetric cuboid geometry
	else if(_es->parameters.get<unsigned int> ("geometry_type") == 5)
	{

	
		if(!threed)
		{
			// for some odd reason, applying dirichlet conditions on one of the walls doesn't work?
			MeshTools::Modification::change_boundary_id(_mesh,0,-2);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,1,1);		//inflow	//should be x=0
			MeshTools::Modification::change_boundary_id(_mesh,2,-1);		//axis
			MeshTools::Modification::change_boundary_id(_mesh,3,0);		//outflow	//should be x=L

		}
		else
		{
			std::cout << "axisymmetric geometry cuboid not physical in 3D" << std::endl;
			std::exit(0);
		}
	
		// subdomain id 0 by default
		//MeshTools::Modification::change_subdomain_id(_mesh,1020,0);		//volume


		// hmmm how many fn boundary ids do we actually have??
		// haven't added the 

		// now we have two wall boundaries we need to get rid of two of them
		boundary_ids.resize(_mesh.boundary_info->n_boundary_ids() - 2);

		for(unsigned int i=0; i < boundary_ids.size(); i++)
			boundary_ids[i] = i;
		

		//set surfaceboundary stuff
		//inflow_surface_boundary_object.init(_mesh,0);

		// unfortunately this does not work for libmesh created geoms because ids not propagated to nodes
		
		for(unsigned int i=0; i < boundary_ids.size(); i++)
		{
			SurfaceBoundary* surface_boundary = new SurfaceBoundary(*_es);
			surface_boundaries.push_back(surface_boundary);
			surface_boundaries[i]->init(_mesh,i,true);

			std::cout << "surface " << i << ": normal is " << surface_boundaries[i]->get_normal() << std::endl;
			std::cout << "  and centroid is " << surface_boundaries[i]->get_centroid() << std::endl << std::endl;;
		}
		


		//hard coded
		// after the first iteration
		pressure_values_3d.push_back(1.0);	//inflow
		pressure_values_3d.push_back(0.0);	//outflow

		// post processed values
		flux_values_3d.push_back(1.0);	//inflow
		flux_values_3d.push_back(0.0);	//outflow

		previous_flux_values_3d = flux_values_3d;
		previous_previous_flux_values_3d = flux_values_3d;
	
		input_pressure_values_3d.push_back(_es->parameters.get<double>("parent_pressure_mag"));
		input_pressure_values_3d.push_back(_es->parameters.get<double>("daughter_1_pressure_mag"));

		// scale so that corresponds to the correct generation
		//MeshTools::Modification::scale(mesh,5.6);
	}

	// the cuboid geometry
	else if(_es->parameters.get<unsigned int> ("geometry_type") == 3)
	{

		if(!threed)
		{
			MeshTools::Modification::change_boundary_id(_mesh,0,-1);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,1,-1);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,2,0);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,3,-1);		//walls	
		}
		else
		{

			MeshTools::Modification::change_boundary_id(_mesh,0,-1);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,1,-1);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,2,-1);		//walls
			MeshTools::Modification::change_boundary_id(_mesh,3,0);		//walls	
			MeshTools::Modification::change_boundary_id(_mesh,4,-1);			//inflow
			MeshTools::Modification::change_boundary_id(_mesh,5,-1);			//outflows
		}
		
		// subdomain id 0 by default
		//MeshTools::Modification::change_subdomain_id(_mesh,1020,0);		//volume



		// hmmm how many fn boundary ids do we actually have??
		// haven't added the 
		boundary_ids.resize(_mesh.boundary_info->n_boundary_ids() - 1);	//hmmm not workin
		//boundary_ids.resize(_mesh.boundary_info->n_boundary_ids());	//hmmm not workin
		for(unsigned int i=0; i < boundary_ids.size(); i++)
			boundary_ids[i] = i;
		

		//set surfaceboundary stuff
		//inflow_surface_boundary_object.init(_mesh,0);

		// unfortunately this does not work for libmesh created geoms because ids not propagated to nodes
		
		for(unsigned int i=0; i < boundary_ids.size(); i++)
		{
			SurfaceBoundary* surface_boundary = new SurfaceBoundary(*_es);
			surface_boundaries.push_back(surface_boundary);
			surface_boundaries[i]->init(_mesh,i,true);

			std::cout << "surface " << i << ": normal is " << surface_boundaries[i]->get_normal() << std::endl;
			std::cout << "  and centroid is " << surface_boundaries[i]->get_centroid() << std::endl << std::endl;;
		}
		

		//hard coded
		// not important to preset these 3d values because they will be reset 
		// after the first iteration
		pressure_values_3d.push_back(1.0);	//inflow
		pressure_values_3d.push_back(0.0);	//outflow

		// post processed values
		flux_values_3d.push_back(1.0);	//inflow
		flux_values_3d.push_back(0.0);	//outflow

		previous_flux_values_3d = flux_values_3d;
		previous_previous_flux_values_3d = flux_values_3d;

	
		input_pressure_values_3d.push_back(_es->parameters.get<double>("parent_pressure_mag"));
		input_pressure_values_3d.push_back(_es->parameters.get<double>("daughter_1_pressure_mag"));

		// scale so that corresponds to the correct generation
		//MeshTools::Modification::scale(_mesh,5.6);
	}

	// *************** SET UP SUBDOMAINS IN 3D **************** //
	// can have any number of types of elements each with a different subdomainid
	// we start at 0
	// if restart then we want to do something slightly different and find out what 
	// 3D elements were used and 

	//  assign different subdomain id to different element types
	// Loop over all elements.

	std::cout << "num subdomains before = " <<	_mesh.n_subdomains () << std::endl;
	std::cout << "num subdomains before = " <<	mesh.n_subdomains () << std::endl;
	
	std::vector<libMesh::ElemType> elem_types;	// temporary variable storing the types of elements encountered
	MeshBase::element_iterator       elem_it  = mesh.elements_begin();
	const MeshBase::element_iterator elem_end = mesh.elements_end();
	
	std::vector<libMesh::ElemType>::iterator it;

	unsigned int count = 0;
	for (; elem_it != elem_end; ++elem_it)
	{
		Elem* elem = *elem_it;
		ElemType current_elem_type = elem->type();
		it = std::find(elem_types.begin(), elem_types.end(), current_elem_type);

		// element type not encountered yet
		if(it==elem_types.end())
		{
			elem_types.push_back(current_elem_type);
			elem->subdomain_id() = count;
			count++;
				
		}
		else
		{
			elem->subdomain_id() = it - elem_types.begin();
			//count++;	// not sure if this is right but makes sense
		}
	}
	

	std::cout << std::endl;
	std::cout << "Subdomain identification summary:" << std::endl;
	std::cout << " " << elem_types.size() << " element types found." << std::endl;
	std::cout << " Given subdomain ids:" << std::endl;
	for(unsigned int i=0; i<elem_types.size(); i++)
	{
		std::cout << " subdomain " << i << ": " << Utility::enum_to_string(elem_types[i]) << " elements." << std::endl;
	}
	

	//populate subdomains_3d
	for(unsigned int i=0; i<elem_types.size(); i++)
		subdomains_3d.push_back(i);



}

// ************************************************************************** //








// ************************** SETUP 3D SYSTEM ******************************* //

void NavierStokesCoupled::setup_3d_system(TransientLinearImplicitSystem* system)
{
	std::cout << "Setting up 2D/3D system." << std::endl;

	//some parameters we need
	bool linear_shape_functions = es->parameters.get<bool>("linear_shape_functions");
	unsigned int  problem_type = es->parameters.get<unsigned int> ("problem_type");
	//Real  restart_time = es->parameters.get<Real> ("restart_time");
	//const MeshBase& mesh = es->get_mesh();

	std::set<subdomain_id_type> active_subdomains;
	active_subdomains.clear(); 
	for(unsigned int i=0; i<subdomains_3d.size(); i++)
		active_subdomains.insert(subdomains_3d[i]);

	unsigned int u_var = 0, v_var = 0, w_var = 0;
	unsigned int u_adj_var = 0, v_adj_var = 0, w_adj_var = 0;
	//unsigned int g_x_var = 0, g_y_var = 0, g_z_var = 0;

	
  // Add the variables "u" & "v" to "Navier-Stokes".  They
  // will be approximated using second-order approximation.

	if(!linear_shape_functions)
	{
	  	u_var = system->add_variable ("u", SECOND,LAGRANGE,&active_subdomains);
  		v_var = system->add_variable ("v", SECOND,LAGRANGE,&active_subdomains);
			if(threed)
	  		w_var = system->add_variable ("w", SECOND,LAGRANGE,&active_subdomains);
	}
	else
	{
	  	u_var = system->add_variable ("u", FIRST,LAGRANGE,&active_subdomains);
  		v_var = system->add_variable ("v", FIRST,LAGRANGE,&active_subdomains);
			if(threed)
	  		w_var = system->add_variable ("w", FIRST,LAGRANGE,&active_subdomains);
	}

	if(!es->parameters.get<bool> ("constant_pressure"))
  	system->add_variable ("p", FIRST,LAGRANGE,&active_subdomains);
	else
		system->add_variable ("p", CONSTANT,MONOMIAL,&active_subdomains);

	if(es->parameters.get<bool>("optimisation_stabilised"))
	{
		if(!linear_shape_functions)
		{
				u_adj_var = system->add_variable ("u_adj", SECOND,LAGRANGE,&active_subdomains);
				v_adj_var = system->add_variable ("v_adj", SECOND,LAGRANGE,&active_subdomains);
				if(threed)
					w_adj_var = system->add_variable ("w_adj", SECOND,LAGRANGE,&active_subdomains);
		}
		else
		{
				u_adj_var = system->add_variable ("u_adj", FIRST,LAGRANGE,&active_subdomains);
				v_adj_var = system->add_variable ("v_adj", FIRST,LAGRANGE,&active_subdomains);
				if(threed)
					w_adj_var = system->add_variable ("w_adj", FIRST,LAGRANGE,&active_subdomains);
		}

		if(es->parameters.get<bool> ("discontinuous_linear_neumann"))
		{
				//g_x_var = system_neumann->add_variable ("g_x", FIRST,L2_LAGRANGE,&active_subdomains);
				//g_y_var = system_neumann->add_variable ("g_y", FIRST,L2_LAGRANGE,&active_subdomains);
				system_neumann->add_variable ("g_x", FIRST,L2_LAGRANGE,&active_subdomains);
				system_neumann->add_variable ("g_y", FIRST,L2_LAGRANGE,&active_subdomains);
				if(threed)
				{
					//g_z_var = system_neumann->add_variable ("g_z", FIRST,L2_LAGRANGE,&active_subdomains);
					system_neumann->add_variable ("g_z", FIRST,L2_LAGRANGE,&active_subdomains);
				}
		}

		if(!es->parameters.get<bool> ("constant_pressure"))
			system->add_variable ("p_adj", FIRST,LAGRANGE,&active_subdomains);
		else
			system->add_variable ("p_adj", CONSTANT,MONOMIAL,&active_subdomains);
	}

	// ******** WE NEED TO ADD THE PRECONDITIONER MATRIX POSSIBLY ********** //
	if(es->parameters.get<unsigned int>("preconditioner_type_3d") || es->parameters.get<unsigned int>("preconditioner_type_3d1d"))
	{
		system->add_matrix("Pressure Laplacian Matrix");
		system->add_matrix("Pressure Mass Matrix");
		system->add_matrix("Pressure Convection Diffusion Matrix");
		system->add_matrix("Velocity Mass Matrix");
		system->add_matrix("Preconditioner");
		system->add_matrix("Velocity Matrix");	// a full matrix with only the velocity block
	}

	//******* NOW WE NEED TO TAKE CARE OF THE BOUNDARY CONDITION STUFF *****//

	std::cout << "Setting up boundary conditions." << std::endl;
	// note depending on the geometry input we may use a dirichlet boundary object
	// or not... thing is sometimes we can't label the boundaries...

	//attach dirichlet boundary condition objects

	//make variables object for boundary conditions
	std::vector<unsigned int> variables;
	variables.push_back(u_var);
	variables.push_back(v_var);

	if(threed)
 		variables.push_back(w_var);

	//make variables object for axisymmetric boundary condition
	std::vector<unsigned int> variables_axi;
	variables_axi.push_back(v_var);
	

	// the inflow boundary condition

	ZeroFunction<> zf;

	// the only ones that can run an axissymetric simulation are 0, 1, 2, 5
	if(problem_type == 0)
	{

		// the wall boundary condition
		std::set<boundary_id_type> boundary_id_wall;
			boundary_id_wall.insert(-1);


		DirichletBoundary dirichlet_bc_wall(boundary_id_wall,
			                       variables,
			                       &zf);

		system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_wall);

		std::set<boundary_id_type> boundary_id_inflow;
				boundary_id_inflow.insert(0);
				boundary_id_inflow.insert(-2);


		if(threed)
		{
			std::vector<unsigned int> w_variables;
		  		w_variables.push_back(w_var);

			//parabolic
			InflowBoundaryFunction<>* p_inflow;

			// doing boundaries the old way for libmesh created meshes
			//if(es->parameters.get<unsigned int> ("geometry_type") == 2
			//		|| es->parameters.get<unsigned int> ("geometry_type") == 3)
			//	p_inflow = new InflowBoundaryFunction<>(*es,NULL);			
			//else
				p_inflow = new InflowBoundaryFunction<>(*es,&*surface_boundaries[0]);



			DirichletBoundary dirichlet_bc_inflow(boundary_id_inflow,
		                       variables,
		                       p_inflow);


			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow);

			/*
			std::vector<unsigned int> uv_variables;
		  		uv_variables.push_back(u_var);
		  		uv_variables.push_back(v_var);

				DirichletBoundary dirichlet_bc_inflow_0(boundary_id_inflow,
			                       uv_variables,
			                       &zf);

			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow_0);	
			*/
		}
		else
		{
			//want in 2d for the velocity to go in the x direction
			std::vector<unsigned int> u_variables;
		  		u_variables.push_back(u_var);

			//parabolic
			InflowBoundaryFunction<>* p_inflow;

			// doing boundaries the old way for libmesh created meshes
			//if(es->parameters.get<unsigned int> ("geometry_type") == 2
			//		|| es->parameters.get<unsigned int> ("geometry_type") == 3)
			//{
			//	p_inflow = new InflowBoundaryFunction<>(*es,NULL);
				//p_inflow = &inflow;			
			//}
			//else
				p_inflow = new InflowBoundaryFunction<>(*es,&*surface_boundaries[0]);



			DirichletBoundary dirichlet_bc_inflow(boundary_id_inflow,
		                       variables,
		                       p_inflow);


			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow);

			/*
			std::vector<unsigned int> v_variables;
		  		v_variables.push_back(v_var);

			DirichletBoundary dirichlet_bc_inflow_0(boundary_id_inflow,
		                       v_variables,
		                       &zf);

			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow_0);	

			*/
		}


		if(es->parameters.get<bool>("optimisation_stabilised"))
		{

			//make variables object for boundary conditions
			std::vector<unsigned int> variables_adj;
			variables_adj.push_back(u_adj_var);
			variables_adj.push_back(v_adj_var);
			
			if(threed)
		 		variables_adj.push_back(w_adj_var);

  		DirichletBoundary dirichlet_bc_inflow_0_adj(boundary_id_inflow,
	                         variables_adj,
	                         &zf);
			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow_0_adj);

  		DirichletBoundary dirichlet_bc_wall_adj(boundary_id_wall,
	                         variables_adj,
	                         &zf);

			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_wall_adj);

			
			//axisymmetric boundary condition
			if(es->parameters.get<unsigned int> ("geometry_type") == 5)
			{

				std::vector<unsigned int> variables_axi_adj;
				variables_axi_adj.push_back(v_adj_var);

				std::set<boundary_id_type> boundary_id_axi_adj;
					boundary_id_axi_adj.insert(-2);


				DirichletBoundary dirichlet_bc_axi_adj(boundary_id_axi_adj,
							                   variables_axi_adj,
							                   &zf);

				system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_axi_adj);
			}
			
		}

	}
	else if(problem_type == 3)
	{

		std::vector<unsigned int> w_variables;
    		w_variables.push_back(w_var);
		std::set<boundary_id_type> boundary_id_inflow;
				boundary_id_inflow.insert(0);

 		boundary_id_inflow.insert(-1);	//we also have the boundary conditions of the wall

		//parabolic
		MovingPipeBoundaryFunction<> inflow(*es);

  		DirichletBoundary dirichlet_bc_inflow(boundary_id_inflow,
	                         w_variables,
	                         &inflow);


		system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow);

		std::vector<unsigned int> uv_variables;
    		uv_variables.push_back(u_var);
    		uv_variables.push_back(v_var);

  		DirichletBoundary dirichlet_bc_inflow_0(boundary_id_inflow,
	                         uv_variables,
	                         &zf);

		system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow_0);


	}
	else if(problem_type == 1 || problem_type == 2)
	{

		// only the wall boundary condition
		std::set<boundary_id_type> boundary_id_wall;
			boundary_id_wall.insert(-1);


		DirichletBoundary dirichlet_bc_wall(boundary_id_wall,
			                       variables,
			                       &zf);

		system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_wall);

		//axisymmetric boundary condition
		if(es->parameters.get<unsigned int> ("geometry_type") == 5)
		{
			std::set<boundary_id_type> boundary_id_axi;
				boundary_id_axi.insert(-2);


			DirichletBoundary dirichlet_bc_axi(boundary_id_axi,
					                     variables_axi,
					                     &zf);

			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_axi);
		}


		if(es->parameters.get<bool>("optimisation_stabilised"))
		{

			//make variables object for boundary conditions
			std::vector<unsigned int> variables_adj;
			variables_adj.push_back(u_adj_var);
			variables_adj.push_back(v_adj_var);
			
			if(threed)
		 		variables_adj.push_back(w_adj_var);

  		DirichletBoundary dirichlet_bc_wall_adj(boundary_id_wall,
	                         variables_adj,
	                         &zf);

			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_wall_adj);
			
			//axisymmetric boundary condition
			if(es->parameters.get<unsigned int> ("geometry_type") == 5)
			{

				std::vector<unsigned int> variables_axi_adj;
				variables_axi_adj.push_back(v_adj_var);

				std::set<boundary_id_type> boundary_id_axi_adj;
					boundary_id_axi_adj.insert(-2);


				DirichletBoundary dirichlet_bc_axi_adj(boundary_id_axi_adj,
							                   variables_axi_adj,
							                   &zf);

				system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_axi_adj);
			}
		}

	}
	else if(problem_type == 4)
	{


		if(!es->parameters.get<bool>("pearson_example"))
		{

			if(es->parameters.get<bool>("leaky_lid_driven_cavity"))
			{

				// lid boundary condition and walls in u direction
				std::vector<unsigned int> u_variables;
						u_variables.push_back(u_var);
				std::set<boundary_id_type> boundary_id_inflow;
						boundary_id_inflow.insert(0);
						boundary_id_inflow.insert(-1);

				//parabolic
				LidDrivenCavityBoundaryFunction<> inflow(*es);

					DirichletBoundary dirichlet_bc_inflow(boundary_id_inflow,
					                     u_variables,
					                     &inflow);


				system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow);

				std::vector<unsigned int> vw_variables;
				vw_variables.push_back(v_var);
				if(threed) { vw_variables.push_back(w_var); }

				DirichletBoundary dirichlet_bc_inflow_0(boundary_id_inflow,
					                     vw_variables,
					                     &zf);

				system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow_0);
			}
			else
			{
				// wall boundary condition
				std::set<boundary_id_type> boundary_id_wall;
				boundary_id_wall.insert(-1);


				DirichletBoundary dirichlet_bc_wall(boundary_id_wall,
							                   variables,
							                   &zf);

				system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_wall);


				// lid boundary condition
				std::vector<unsigned int> u_variables;
						u_variables.push_back(u_var);
				std::set<boundary_id_type> boundary_id_inflow;
						boundary_id_inflow.insert(0);

				//parabolic
				LidDrivenCavityBoundaryFunction<> inflow(*es);

					DirichletBoundary dirichlet_bc_inflow(boundary_id_inflow,
					                     u_variables,
					                     &inflow);


				system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow);

				std::vector<unsigned int> vw_variables;
				vw_variables.push_back(v_var);
				if(threed) { vw_variables.push_back(w_var); }

				DirichletBoundary dirichlet_bc_inflow_0(boundary_id_inflow,
					                     vw_variables,
					                     &zf);

				system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow_0);			
			}

		}
		else
		{
			// lid boundary condition
			std::vector<unsigned int> u_variables;
		  		u_variables.push_back(u_var);
			std::set<boundary_id_type> boundary_id_inflow;
					boundary_id_inflow.insert(0);
					boundary_id_inflow.insert(-1);	//prescribe wall condition using function so gets corners correct

			//parabolic
			LidDrivenCavityBoundaryFunction<> inflow(*es);

				DirichletBoundary dirichlet_bc_inflow(boundary_id_inflow,
			                       u_variables,
			                       &inflow);


			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow);

			std::vector<unsigned int> vw_variables;
		  vw_variables.push_back(v_var);
		  if(threed) { vw_variables.push_back(w_var); }

			DirichletBoundary dirichlet_bc_inflow_0(boundary_id_inflow,
			                       vw_variables,
			                       &zf);

			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow_0);

		}

		if(es->parameters.get<bool>("optimisation_stabilised"))
		{

			//make variables object for boundary conditions
			std::vector<unsigned int> variables_adj;
			variables_adj.push_back(u_adj_var);
			variables_adj.push_back(v_adj_var);

			std::set<boundary_id_type> boundary_id_inflow;
					boundary_id_inflow.insert(0);

			std::set<boundary_id_type> boundary_id_wall;
					boundary_id_inflow.insert(0);
			
			if(threed)
		 		variables_adj.push_back(w_adj_var);

  		DirichletBoundary dirichlet_bc_inflow_0_adj(boundary_id_inflow,
	                         variables_adj,
	                         &zf);
			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_inflow_0_adj);

  		DirichletBoundary dirichlet_bc_wall_adj(boundary_id_wall,
	                         variables_adj,
	                         &zf);

			system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_wall_adj);
		}

	}

	else if(problem_type == 5)
	{



		//velocity exact solution
		ExactSolutionVelocity<> exact_velocity(*es);

		// velocity boundary condition
		std::set<boundary_id_type> boundary_id_wall;
		boundary_id_wall.insert(-1);
		boundary_id_wall.insert(0);


		DirichletBoundary dirichlet_bc_velocity(boundary_id_wall,
			                       variables,
			                       &exact_velocity);

		system->get_dof_map().add_dirichlet_boundary(dirichlet_bc_velocity);
	}


	/*
std::cout << "lakka" << std::endl;
	// okay this is just to test the exact solution crap
	ExactSolutionVelocity<> exact_velocity_bessel(*es);
	DenseVector<Number> output;
	output.resize(2);
	std::cout << "k" << std::endl;
	Point pointert(0.,0.,0.);
	std::cout << "k" << std::endl;
	exact_velocity_bessel(pointert,0.,output);
	std::cout << "bessel_result done" << std::endl;
	std::exit(0);
	*/




}

// ************************************************************************** //








// ************************** PREREFINE 3D MESH ***************************** //

void NavierStokesCoupled::prerefine_3d_mesh()
{

 	std::cout << "Prerefining the mesh." << std::endl;

	bool uniform_refinement = false;
	//kelly error estimator
	//oh well have some selection over the type of flagging to be done here.
	if(!uniform_refinement)
	{
		//error_estimator.reset(new KellyErrorEstimator);
  	picard->estimate_error (error_vector);
  	//mesh_refinement.flag_elements_by_nelem_target(error_vector);

		// my own flagging to flag the top elements for refinement
		double y_upper = 1.;
		double y_lower = 0.;
		double x_upper = 1.0;
		double x_lower = 0.75;
		//ZPositionElementFlagging element_flagger(mesh,z_upper,z_lower);

		XYZPositionElementFlagging element_flagger(mesh,x_upper,x_lower,y_upper,y_lower);

		mesh_refinement.clean_refinement_flags();
		mesh_refinement.flag_elements_by(element_flagger);

		mesh_refinement.refine_and_coarsen_elements();
	}
	else
	{
  	mesh_refinement.uniformly_refine (1);
	}

	es->reinit ();
	init_dof_variable_vectors();

	mesh.print_info();

	Real global_error = error_vector.l2_norm();

	// if we have done adaptive refinement then we always need to
	// make sure the old global solution is the correct size (and vals)

	//std::cout << "setting up old_global solution" << std::endl;
	if(sim_type == 5)
	{
	 	old_global_solution->init(*system_coupled->solution);
		//std::cout << "huh" << std::endl;
		*old_global_solution = *system_coupled->solution;	// remember we changed solution to old global before refinement
	}
	else
	{
	 	old_global_solution->init(*system_3d->solution);
		//std::cout << "huh" << std::endl;
		*old_global_solution = *system_3d->solution;	// remember we changed solution to old global before refinement	
	}

	std::cout << "Global_error = " << global_error
        				<< std::endl;
	std::cout << "Worst element error = " << error_vector.maximum()
        				<< ", mean = " << error_vector.mean() << std::endl;
	
	std::cout << "  Adaptive refinement done." << std::endl;

	


	// my own flagging to flag the top elements for refinement
	//Point ball_centre(0.2,0,-2.4);
	//double ball_radius = 0.5;
	//BallBoundaryElementFlagging element_flagger(mesh,ball_centre,ball_radius);

	//need to clean up here cause can't put it into ZPositionElementFlagging
/*

	mesh_refinement.refine_and_coarsen_elements();

	es->reinit ();

	mesh.print_info();
	std::cout << "Refinement done." << std::endl;


            	std::cout << "  Refining the mesh..." << std::endl;


	// my own flagging to flag the top elements for refinement
	z_upper = -2.0;
	z_lower = -3.0;
	ZPositionElementFlagging element_flagger_2(mesh,z_upper,z_lower);


	// my own flagging to flag the top elements for refinement
	//Point ball_centre(0.2,0,-2.4);
	//double ball_radius = 0.5;
	//BallBoundaryElementFlagging element_flagger(mesh,ball_centre,ball_radius);

	//need to clean up here cause can't put it into ZPositionElementFlagging
	mesh_refinement.clean_refinement_flags();
	mesh_refinement.flag_elements_by(element_flagger_2);

	mesh_refinement.refine_and_coarsen_elements();

	es->reinit ();

	std::cout << "Refinement done." << std::endl;
*/
	
}

// ************************************************************************** //








// ******************** CALCULATE 3D BOUNDARY VALUES ************************ //

void NavierStokesCoupled::calculate_3d_boundary_values()
{

	

	if(sim_3d)
	{

		previous_flux_values_3d = flux_values_3d;
		previous_pressure_values_3d = pressure_values_3d;

		flux_values_3d[0] = picard->calculate_flux(0);
		pressure_values_3d[0] = picard->calculate_pressure(0);
		for(unsigned int i=1; i< flux_values_3d.size(); i++)
		{
			if(es->parameters.get<bool> ("match_1d_mesh_to_3d_mesh") && sim_1d)
			{
				flux_values_3d[boundary_id_to_tree_id[i]] = picard->calculate_flux(i);
				pressure_values_3d[boundary_id_to_tree_id[i]] = picard->calculate_pressure(i);
			}
			else
			{
				flux_values_3d[i] = picard->calculate_flux(i);
				pressure_values_3d[i] = picard->calculate_pressure(i);
			}
		}
	}
	else
	{
		std::cout << "How on earth could I calculate boundary values of a 3D model"
							<< " if I'm not even running a 3D model? You, my friend, are a retard" << std::endl;
		exit_program = true;
		//exit(0);
	}
}

// ************************************************************************** //








// *********************** WRITE OUT THE 3D SOLUTION ************************ //
// write the solution at a specific time and time_step, backup=TRUE if backup is to be written

void NavierStokesCoupled::write_3d_solution(bool backup)
{

	std::cout << "Writing 2D/3D solution" << std::endl;

/*
	const char dir_path[] = "/home/james/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results";

	boost::filesystem::path dir(dir_path);
	if(boost::filesystem::create_directory(dir)) {
		std::cout << "Success" << "\n";
*/
	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(sim_type == 5)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}

	// exodus file writer
	ExodusII_IO_Extended exo = ExodusII_IO_Extended(mesh);

	exo.set_var_scalings(var_scalings_3D);

	std::vector<std::string> variables_3d;
	variables_3d.push_back("u");
	variables_3d.push_back("v");
	if(threed)
		variables_3d.push_back("w");
	variables_3d.push_back("p");

	if(es->parameters.get<bool>("output_adjoint_variables"))
	{
		variables_3d.push_back("u_adj");
		variables_3d.push_back("v_adj");
		if(threed)
			variables_3d.push_back("w_adj");
		variables_3d.push_back("p_adj");

	}

	exo.set_output_variables(variables_3d);

	std::ostringstream file_name;
	std::ostringstream file_name_soln;
	std::ostringstream file_name_es;
	std::ostringstream file_name_mesh;

//	file_name << "results/out_3D_viscosity"
//	  	<< es->parameters.get<Real> ("viscosity");

	file_name << output_folder.str() << "out_3D";
 
	// We write the file in the ExodusII format.
	file_name_soln << file_name.str();
	//file_name_soln << ".e" ;
	file_name_soln << ".e-s." ;
	file_name_soln << std::setw(4) << std::setfill('0') << t_step;

	//exo.write_equation_systems (file_name_soln.str(),
   //                             *es);


	// here we're gonna coarsen of refine the mesh for output
	int refinement_level_difference = es->parameters.get<unsigned int>("output_no_refinement") 
																		- es->parameters.get<unsigned int>("no_refinement");

	system->update();
	if(refinement_level_difference > 0 )
	{
		std::cout << "Refining mesh for output." << std::endl;
		for(unsigned int i=0; i<(unsigned int)refinement_level_difference; i++)
		{
			mesh_refinement.uniformly_refine(1);
			es->reinit();	//NEEDED: so that correct xda file is written, of course we lose the accuracy... should actually do a temp copy
		}

		//std::cout << "cannot output solution at a more refined state (need to fix dirichlet condition stuff)" << std::endl;
		//std::cout << "ABORTING" << std::endl;
		//exit(0);
	}
	else if(refinement_level_difference < 0 )
	{
		std::cout << "Coarsening mesh for output." << std::endl;
		for(unsigned int i=0; i<(unsigned int)(-refinement_level_difference); i++)
		{
			mesh_refinement.uniformly_coarsen(1);
			es->reinit();	//NEEDED: so that correct xda file is written, of course we lose the accuracy... should actually do a temp copy
		}
	}		

	//system->update();
	


	//the time step variable needs to reflect the position in the file_name
	// since we want one file per time step we just put 1 here
	std::cout << "Actually writing file." << std::endl;
	exo.write_timestep(file_name_soln.str(), *es,1,time*time_scale_factor);


	if(es->parameters.get<unsigned int>("error_estimator"))
	{
		//calculate the error
		picard->estimate_error(error_vector);

		std::cout << "Worst element error = " << error_vector.maximum()
							<< ", best element error = " << error_vector.minimum()
       				<< ", mean = " << error_vector.mean() << std::endl;
		//maybe write the error and that stuff
		plot_error(exo);
	}

	//
	std::cout << "before write pid" << std::endl;
	write_elem_pid_3d(exo);




	// We write the equation systems and mesh object in .xda format
	// for restarting... perhaps we don't wanna do this every time step
	// but every n time_steps

	std::cout << "EXODUSII output for timestep " << t_step
		<< " written to " << file_name_soln.str() << std::endl;


	if(backup)
	{
		std::cout << "Writing backup files." << std::endl;
		file_name_es << file_name.str();
		file_name_es << "_es_";
		file_name_es << std::setw(4) << std::setfill('0') << t_step;
		file_name_es << ".xda";
		es->write(file_name_es.str(), libMeshEnums::WRITE);

		file_name_mesh << file_name.str();
		file_name_mesh << "_mesh_";
		file_name_mesh << std::setw(4) << std::setfill('0') << t_step;
		file_name_mesh << ".xda";
		mesh.write(file_name_mesh.str());

		std::cout << "Backup files written to " << file_name_mesh.str()
			<< " and " << file_name_es.str() << std::endl;
	}

 
	// now do the opposite
	if(refinement_level_difference > 0 )
	{
		
		std::cout << "Coarsening mesh to return to simulation." << std::endl;
		for(unsigned int i=0; i<(unsigned int)refinement_level_difference; i++)
		{
			mesh_refinement.uniformly_coarsen(1);
			es->reinit();	//NEEDED: so that correct xda file is written, of course we lose the accuracy... should actually do a temp copy
		}
	}
	else if(refinement_level_difference < 0 )
	{

		std::cout << "Refining mesh to return to simulation." << std::endl;
		for(unsigned int i=0; i<(unsigned int)(-refinement_level_difference); i++)
		{
			mesh_refinement.uniformly_refine(1);
			es->reinit();	//NEEDED: so that correct xda file is written, of course we lose the accuracy... should actually do a temp copy
		}
	}		
	

}


// ************************************************************************** //







// ****************** SOLVE 3D SYSTEM OR COUPLED SYSTEM ********************* //

//returns true if wanna break
bool NavierStokesCoupled::solve_3d_system_iteration(TransientLinearImplicitSystem * system)
{

	PetscErrorCode ierr;

	// always use the same command line options, can differentiate between stokes and navier stokes within the command line hopefully	
	// probably set already in read_parameters
	ierr = PetscOptionsInsertString(es->parameters.get<std::string>("petsc_solver_options").c_str()); CHKERRQ(ierr);

	previous_nonlinear_residual = current_nonlinear_residual;


	//es->reinit ();	// UPDATE TIME DEPENDENT DIRICHLET BOUNDARY CONDITIONS put after

	//petsc stuff
	//if we have different solvers for the different blocks then need a prefix
	PetscLinearSolver<Number>* system_linear_solver =
			libmesh_cast_ptr<PetscLinearSolver<Number>* >
			(system->linear_solver.get());


	//****************** SOLVE THE SYSTEM AND OUTPUT SOME INFO ************************//

	*last_nonlinear_soln = *system->solution;	//don't need to worry bout sizes
	last_nonlinear_soln->close();
	//double before_norm = last_nonlinear_soln->l2_norm();

	//set the inner tolerances here - maxits is now 20
	//ierr = KSPSetTolerances(subksp[1],1e-8,1e-10,1e3,30);

	//using the previous solution as an initial guess is bad for the iterative solvers!!
	//if(es->parameters.get<bool>("adaptive_refinement")))
	// if we do not do update here then the local won't be updated and the assembly will use solution

	// must set it to zero for initial guess of solver, but don't update to local for assembly
	system->update();						//put the previous soln in local
	//system->solution->zero();		//now make initial guess zero
	//*system->solution = *last_nonlinear_soln;
		
	//unsigned int n_linear_iterations = 0;
	double final_linear_residual = 0.;

	if(sim_type != 5)
	{
		final_linear_residual = solve_and_assemble_3d_system(system);
	}
	else
	{
		final_linear_residual = solve_and_assemble_3d_system(system);

	}

	if(es->parameters.get<bool>("compute_eigenvalues"))
		compute_and_output_eigenvalues(system);



	// need to zero the preconditioner matrix myself... after we have computed the eigenvalues
	if(es->parameters.get<unsigned int>("preconditioner_type_3d") || es->parameters.get<unsigned int>("preconditioner_type_3d1d"))
	{

		// also zero the assemble bits for next time
		system->request_matrix("Preconditioner")->close();
		system->request_matrix("Preconditioner")->zero();
	}

	// need to recollect the ksp because perhaps some shit changed in restrict_solve_to
	system_linear_solver =
			libmesh_cast_ptr<PetscLinearSolver<Number>* >
			(system->linear_solver.get());
	KSP system_ksp = system_linear_solver->ksp();

	int num_outer_its = 0;	
	int num_inner_velocity_its = 0;
	int num_inner_pressure_its = 0;
	// output some data about the iterative solve
	// this is hard to do when doing restrict_solve_to, so just leave it

	ierr = KSPGetIterationNumber(system_ksp,&num_outer_its); CHKERRQ(ierr);
	if(false)//es->parameters.get<bool>("fieldsplit"))
	{
		//get the pc
		PC pc;
		ierr = KSPGetPC(system_ksp,&pc); CHKERRQ(ierr);

		//get the inner iteration ksp
		int num_splits = 0;
		KSP* subksp;	//subksps
		ierr = PCFieldSplitGetSubKSP(pc,&num_splits,&subksp); CHKERRQ(ierr);

		ierr = KSPGetIterationNumber(subksp[0],&num_inner_velocity_its); CHKERRQ(ierr);
		ierr = KSPGetIterationNumber(subksp[1],&num_inner_pressure_its); CHKERRQ(ierr);

		//std::cout << "num outer its = " << num_outer_its << std::endl;
		//std::cout << "num inner velocity its = " << num_inner_velocity_its << std::endl;
		//std::cout << "num inner pressure (schur) its = " << num_inner_pressure_its << std::endl;

		KSPConvergedReason velocity_converged_reason;
		ierr = KSPGetConvergedReason(subksp[0],&velocity_converged_reason); CHKERRQ(ierr);
		KSPConvergedReason pressure_converged_reason;
		ierr = KSPGetConvergedReason(subksp[1],&pressure_converged_reason); CHKERRQ(ierr);

		//std::cout << "velocity converged reason = " << velocity_converged_reason << std::endl;
		//std::cout << "pressure converged reason = " << pressure_converged_reason << std::endl;
	}
	PetscInt outer_maxits = 0;
	ierr = KSPGetTolerances(system_ksp,NULL,NULL,NULL,&outer_maxits); CHKERRQ(ierr);

	int converged_reason = system_linear_solver->get_petsc_converged_reason();

	//std::cout << "KSPConvergedReason = " << converged_reason << std::endl;
	
	//system->get_dof_map().enforce_constraints_exactly(*system);

	std::vector<Real> weights;  // u, v
	weights.push_back(1.0);
	weights.push_back(1.0);
	if(threed)
		weights.push_back(1.0);
	weights.push_back(0.0);            // p

	if(sim_type == 5)
	{
		weights.push_back(0.0);            // P
		weights.push_back(0.0);            // Q

	}
	

	if(es->parameters.get<bool>("optimisation_stabilised"))
	{
		weights.push_back(0.0);
		weights.push_back(0.0);
		if(threed)
			weights.push_back(0.0);
  	weights.push_back(0.0);            // p
	}
	std::vector<FEMNormType>	norms;	//if empty will automagically use discrete l2

	double norm_3d = system->calculate_norm(*last_nonlinear_soln,SystemNorm(norms, weights));

  // Compute the difference between this solution and the last
  // nonlinear iterate.
	last_nonlinear_soln->add (-1., *system->solution);
	last_nonlinear_soln->close();


	current_nonlinear_residual = system->calculate_norm(*last_nonlinear_soln,SystemNorm(norms, weights));
	// now normalise it
	current_nonlinear_residual /= system->calculate_norm(*system->solution,SystemNorm(norms, weights));


	//current_nonlinear_residual = last_nonlinear_soln->l2_norm();
	num_linear_iterations = system->n_linear_iterations();

	// Print out convergence
	std::cout << "Lin Solver converged ("
						<< converged_reason << ") at step: "
        		<< num_outer_its
        		<< ",\n\tSolver res: "
        		<< final_linear_residual
        		<< ",\n\tNonlinear res: "
        		<< current_nonlinear_residual
        		<< std::endl;
	std::cout << "Nonlinear iteration: " << nonlinear_iteration << std::endl;
	std::cout << "Time: " << time << std::endl;
	std::cout << "Time step: " << t_step << std::endl;
	std::cout << "3D norm: " << norm_3d << std::endl << std::endl;

	es->parameters.set<double> ("last_nonlinear_iterate") = current_nonlinear_residual;

	if(es->parameters.get<bool> ("output_linear_iteration_count"))
		output_linear_iteration_count();


	//*************** DECIDE WHAT TO DO ON THE NEXT ITERATION *****************//

	//if stokes we only need a single iteration unless iterating with 1d model
	// i.e. need iterations
	
	if(es->parameters.get<bool>("stokes") && sim_type != 3)
	{
		std::cout << "Stokes so only one iteration required" << std::endl;
		return true;
	}
	

	// Terminate the solution iteration if the difference between
	// this nonlinear iterate and the last is sufficiently small, AND
	// if the most recent linear system was solved to a sufficient tolerance.
  
	//doesn't make sense to compare with the first step in a time dependent problem
	if(num_outer_its == outer_maxits)
	{
		std::cout << " Nonlinear solver is not going to converge because "
          		<< "iterative solver is taking not converging" << std::endl;
		if(es->parameters.get<bool>("adaptive_time_stepping"))
		{
			*system->solution = *old_global_solution;
			reduce_dt = true;
		}
		else
		{
			std::cout << "no adaptive time stepping. exiting program." << std::endl;
			exit_program = true;			
			//exit(0);
		}
		return true;
	}
	if(es->parameters.get<bool>("increasing_residuals_exit") && 
			(current_nonlinear_residual > previous_nonlinear_residual && nonlinear_iteration > 2))	
	{
		std::cout << " Nonlinear solver is not going to converge because "
          		<< "nonlinear residuals are increasing."
          		<< std::endl;
		if(es->parameters.get<bool>("adaptive_time_stepping"))
		{
			*system->solution = *old_global_solution;
			reduce_dt = true;
		}
		else
		{
			std::cout << "no adaptive time stepping. exiting program." << std::endl;
			exit_program = true;			
			//exit(0);
		}
		return true;
	}
	else if ((current_nonlinear_residual < es->parameters.get<Real>("nonlinear_tolerance"))
       		)//&& (system_3d->final_linear_residual() < nonlinear_tolerance))
	{
		std::cout << " Nonlinear solver converged at step "
          		<< nonlinear_iteration
          		<< std::endl;

		//well l < 3 was doing very badly..

		if(nonlinear_iteration < es->parameters.get<unsigned int>("adaptive_newton_solve_limit") && 
				steps_since_last_dt_change > es->parameters.get<unsigned int>("min_steps_between_dt_increase") 
				&& es->parameters.get<bool>("adaptive_time_stepping"))	// don't wanna increase the timestep
		{
			increase_dt = true;
		}
  	return true;
 	}
	else if(nonlinear_iteration > es->parameters.get<unsigned int>("max_newton_iterations"))
	{

  			std::cout << " Nonlinear solver did not converge by step "
              		<< nonlinear_iteration
              		<< std::endl;

		if(es->parameters.get<bool>("adaptive_time_stepping"))
		{
			*system->solution = *old_global_solution;
			reduce_dt = true;
		}
		else
		{
			std::cout << "no adaptive time stepping. exiting program." << std::endl;
			exit_program = true;			
			//exit(0);
		}
		return true;
	}


	return false;

}


// ************************************************************************** //

PetscErrorCode custom_outer_monitor(KSP ksp,PetscInt n,PetscReal rnorm,void *dummy)
{
	PetscErrorCode ierr;

	//get the pc
	PC pc;
	ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);

	//get the inner iteration ksp
	int num_splits = 0;
	KSP* subksp;	//subksps
	ierr = PCFieldSplitGetSubKSP(pc,&num_splits,&subksp); CHKERRQ(ierr);

	int num_outer_its;
	int num_inner_velocity_its;
	int num_inner_pressure_its;
	ierr = KSPGetIterationNumber(ksp,&num_outer_its); CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(subksp[0],&num_inner_velocity_its); CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(subksp[1],&num_inner_pressure_its); CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"  in ksp custom monitor: num_pressure_its = %D\n",num_inner_pressure_its); CHKERRQ(ierr);

	//total_inner_pressure_its += num_inner_pressure_its;
  return(0);
}


PetscErrorCode custom_inner_pressure_monitor(KSP ksp,PetscInt n,PetscReal rnorm,void *dummy)
{
	PetscErrorCode ierr;

	//get the pc
	PC pc;
	ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);

	//get the inner iteration ksp
	int num_splits = 0;
	KSP* subksp;	//subksps
	ierr = PCFieldSplitGetSubKSP(pc,&num_splits,&subksp); CHKERRQ(ierr);

	int num_outer_its;
	int num_inner_velocity_its;
	int num_inner_pressure_its;
	ierr = KSPGetIterationNumber(ksp,&num_outer_its); CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(subksp[0],&num_inner_velocity_its); CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(subksp[1],&num_inner_pressure_its); CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"    in ksp custom pressure monitor: num_velocity_its = %D\n",num_inner_velocity_its); CHKERRQ(ierr);

	//total_inner_velocity_its += num_inner_velocity_its;
  return(0);
}


// ************************************************************************** //
// ********              SOLVE 3D LINEAR SYSTEM (replacement of libmesh)	
// ************

// NOTES:
// - no subset functionality
double NavierStokesCoupled::solve_and_assemble_3d_system(TransientLinearImplicitSystem * system)
{
	PetscErrorCode ierr;
		
	perf_log.push("assembly");
	if (system->assemble_before_solve)
	{
		// Assemble the linear system
		system->assemble ();
	}
  perf_log.pop("assembly");


	perf_log.push("solve_pre_setup");


	// ******************* SOME SETTINGS
	es->parameters.set<Real> ("linear solver tolerance") =
		 es->parameters.get<double>("outer_solver_rtol");//1.e-12;
	//set previous nonlinear residual to value it was

	//get petsc solver
	PetscLinearSolver<Number>* system_linear_solver =
			libmesh_cast_ptr<PetscLinearSolver<Number>* >
			(system->linear_solver.get());

	//let us set up to use the same preconditioner
	system_linear_solver->reuse_preconditioner(es->parameters.get<bool>("reuse_preconditioner"));
	
	// this inits the ksp so need to set prefix before this
	std::string prefix_3d;
	if(sim_type != 5)
		prefix_3d = "ns3d_";
	else
		prefix_3d = "ns3d1d_";

	system_linear_solver->set_prefix (prefix_3d);
	KSP system_ksp = system_linear_solver->ksp();

	// all tolerances set at command line
	//ierr = KSPSetTolerances(system_ksp,es->parameters.get<double>("outer_solver_rtol"),es->parameters.get<double>("outer_solver_atol"),1e38,es->parameters.get<int>("outer_solver_maxits")); CHKERRQ(ierr);

	// If the linear solver hasn't been initialized, we do so here.
	//if (libMesh::on_command_line("--solver_system_names"))
	//	system->linear_solver->init((system->name()+"_").c_str());
	//else

	// only need to do fieldsplit IS setup on first nonlinear iteration
	// hmmm, let us check whether all the fieldsplit options have been set here.
	



	// Get the user-specifiied linear solver tolerance
	const Real tol = es->parameters.get<Real>("linear solver tolerance");

	// Get the user-specified maximum # of linear solver iterations
	const unsigned int maxits = es->parameters.get<unsigned int>("linear solver maximum iterations");


	// SET THE MATRIX OPERATORS









	if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 0 && es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 0)
		system_linear_solver->solve_simple_setup (*system->request_matrix("System Matrix"), *system->solution, *system->rhs, tol, maxits);
	else
		system_linear_solver->solve_simple_setup (*system->request_matrix("System Matrix"), *system->request_matrix("Preconditioner"), *system->solution, *system->rhs, tol, maxits);


	if(nonlinear_iteration == 1 && es->parameters.get<bool>("fieldsplit"))
		system->linear_solver->init_names(*system);

//	PetscInt numsplit;
//	ierr = PCFieldSplitGetSubKSP(sub_pc,&numsplit,NULL);// CHKERRQ(ierr);
//	std::cout << "hmmm" << std::endl;
//	std::cout << "numsplit = " << numsplit << std::endl;
//	std::cout << "hmmm" << std::endl;


/*
	if(nonlinear_iteration == 1 && es->parameters.get<bool>("fieldsplit"))
		system->linear_solver->init_names(*system);

	if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 0 && es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 0)
		system_linear_solver->solve_simple_setup (*system->request_matrix("System Matrix"), *system->solution, *system->rhs, tol, maxits);
	else
		system_linear_solver->solve_simple_setup (*system->request_matrix("System Matrix"), *system->request_matrix("Preconditioner"), *system->solution, *system->rhs, tol, maxits);
*/



	// ************************************* //


	if(es->parameters.get<bool>("nonzero_initial_guess"))
	{
	  ierr = KSPSetInitialGuessNonzero (system_ksp, PETSC_TRUE); CHKERRQ(ierr);
	}
	else
	{
	  ierr = KSPSetInitialGuessNonzero (system_ksp, PETSC_FALSE); CHKERRQ(ierr);
	}


	if(false)//es->parameters.get<bool>("fieldsplit"))
	{
		PetscErrorCode (*function_ptr)(KSP, PetscInt, PetscReal, void*);
		function_ptr = &custom_outer_monitor;

		//ierr = KSPMonitorSet(system_ksp,function_ptr,NULL,NULL); CHKERRQ(ierr);
	}

	//get the inner iteration ksp
	int num_splits = 0;


/*
	std::cout << "ja" << std::endl;
	KSP* subksp;	//subksps
	PCFieldSplitGetSubKSP(pc,&num_splits,&subksp);
	std::cout << "ja" << std::endl;
	PetscErrorCode (*function_ptr_pressure)(KSP, PetscInt, PetscReal, void*);
	function_ptr_pressure = &custom_inner_pressure_monitor;
	KSPMonitorSet(subksp[1],function_ptr_pressure,NULL,NULL);
	

	std::cout << "ja" << std::endl;
*/

	//system->request_matrix("System Matrix")->print();

	//std::ofstream matrix_to_print("results/matrix_gosh.txt");
	//system->request_matrix("System Matrix")->print(matrix_to_print);
	//matrix_to_print.close();


	// ******* SETUP PRECONDITIONERS ************** //



	Mat B;

	if(es->parameters.get<unsigned int>("preconditioner_type_3d") || es->parameters.get<unsigned int>("preconditioner_type_3d1d"))
	{

		setup_preconditioners(system,B);

	}



	if(es->parameters.get<bool>("ksp_view"))
	{
		ierr = KSPView(system_ksp,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}


	//std::cout << "System Matrix:"<< std::endl;
	//system->request_matrix("System Matrix")->print();

	std::cout << std::endl;
	perf_log.pop("solve_pre_setup");
	perf_log.push("solve");

	//ierr = PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, s->myS);CHKERRQ(ierr);

	// Solve the linear system.  Several cases:
	std::pair<unsigned int, Real> rval = std::make_pair(0,0.0);
	bool _shell_matrix = false;
	// matrix always passed by reference, preconditioner can be passed by reference of pointer

	// 2.) No shell matrix, with or without user-supplied preconditioner
	
	/*
	std::cout << "system matrix" << std::endl;
	system->request_matrix("System Matrix")->print();
	std::cout << "rhs" << std::endl;
	system->rhs->print();
	std::cout << "solution" << std::endl;
	system->solution->print();
	*/

	//std::cout << "system matrix" << std::endl;
	//system->request_matrix("System Matrix")->print();

	//std::cout << "preconditioner matrix" << std::endl;
	//system->request_matrix("Preconditioner")->print();

	std::cout << "preconditioner_type_3d = " << es->parameters.get<unsigned int>("preconditioner_type_3d") << std::endl;
	std::cout << "preconditioner_type_3d1d = " << es->parameters.get<unsigned int>("preconditioner_type_3d1d") << std::endl;
	

	if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 0 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 0)
	{
		std::cout << "Using system matrix as preconditioner" << std::endl;
		rval = system_linear_solver->solve_simple (*system->request_matrix("System Matrix"), *system->solution, *system->rhs, tol, maxits);
	}
	else
	{
		std::cout << "Using different preconditioner matrix" << std::endl;
		rval = system_linear_solver->solve_simple (*system->request_matrix("System Matrix"), *system->request_matrix("Preconditioner"), *system->solution, *system->rhs, tol, maxits);
	}
	
	/*
	std::cout << "\n\n\n\nsolving for unit rhs" << std::endl;
	system->rhs->zero();
	std::cout << "k" << std::endl;
	system->rhs->set(100,1.);
	system->rhs->close();
	std::cout << "doing the solve" << std::endl;
	rval = system_linear_solver->solve_simple (*system->request_matrix("System Matrix"), *system->solution, *system->rhs, tol, maxits);

	std::cout << "DONE" << std::endl;
	*/


/*
	if(es->parameters.get<unsigned int>("preconditioner_type") == 0)
		rval = system_linear_solver->solve (*system->request_matrix("System Matrix"), *system->solution, *system->rhs, tol, maxits);
	else if(es->parameters.get<unsigned int>("preconditioner_type"))
		rval = system_linear_solver->solve (*system->request_matrix("System Matrix"), *system->request_matrix("Preconditioner"), *system->solution, *system->rhs, tol, maxits);	
*/
	perf_log.pop("solve");

	//system->request_matrix("Preconditioner")->print();
	//if(es->parameters.get<unsigned int>("preconditioner_type") == 1)
	//	int ierr = MatView(schur_complement_approx->mat(),PETSC_VIEWER_STDOUT_SELF);

	/*
	std::cout << "post solution" << std::endl;
	system->solution->print();
	*/


	int num_outer_its = 0;	
	int num_inner_velocity_its = 0;
	int num_inner_pressure_its = 0;
	// output some data about the iterative solve
	// this is hard to do when doing restrict_solve_to, so just leave it

	ierr = KSPGetIterationNumber(system_ksp,&num_outer_its); CHKERRQ(ierr);
	
	PetscReal rnorm;
	KSPGetResidualNorm(system_ksp,&rnorm);

	PetscReal divtol;
	KSPGetTolerances(system_ksp,NULL,NULL,&divtol,NULL);
	

	//double ave_pressure_its = (double)total_inner_pressure_its/(double)num_outer_its;
	//std::cout << "Ave pressure its = " << ave_pressure_its << std::endl;
	//double ave_velocity_its = (double)total_inner_pressure_its/(double)total_inner_velocity_its;
	//std::cout << "Ave velocity its = " << ave_pressure_its << std::endl;

	total_linear_iterations += num_outer_its;
	local_linear_iterations += num_outer_its;
	if(num_outer_its > total_max_iterations)
		total_max_iterations = num_outer_its;

	std::cout << std::endl;
	std::cout << "outer its this timestep = " << local_linear_iterations << std::endl;
	std::cout << "total outer its = " << total_linear_iterations << std::endl;
	std::cout << "average outer its per nonlinear iteration this time step = " << (double)local_linear_iterations/(double)(nonlinear_iteration) << std::endl;
	std::cout << "max outer its = " << total_max_iterations << std::endl;
	std::cout << std::endl;


	//_n_linear_iterations   = rval.first;
	double final_linear_residual = rval.second;

	// Update the system after the solve
	system->update();

	// delete matrices after solve

	std::cout << "k" << std::endl;
	if(es->parameters.get<unsigned int>("preconditioner_type_3d") != 0 && es->parameters.get<unsigned int>("preconditioner_type_3d") != 8)
	{
		
		delete pressure_mass_matrix;
		delete scaled_pressure_mass_matrix;
		delete pressure_laplacian_matrix;
		delete pressure_convection_diffusion_matrix;
		delete velocity_mass_matrix;

		// also zero the assemble bits for next time
		std::cout << "yas queen 1" << std::endl;
		system->request_matrix("Pressure Mass Matrix")->zero();
		std::cout << "yas queen 2" << std::endl;
		system->request_matrix("Pressure Laplacian Matrix")->zero();
		std::cout << "yas queen 3" << std::endl;
		system->request_matrix("Pressure Convection Diffusion Matrix")->zero();
		std::cout << "yas queen 4" << std::endl;
		system->request_matrix("Velocity Mass Matrix")->zero();
		std::cout << "yas queen 5" << std::endl;
		//system->request_matrix("Velocity Matrix")->close();
		//system->request_matrix("Velocity Matrix")->zero();
		//std::cout << "yas queen 6" << std::endl;

		if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 5)
		{
			ierr = MatDestroy(&B); CHKERRQ(ierr);
		}
  	//ierr = MatDestroy(&mat_ctx->pressure_laplacian_matrix);CHKERRQ(ierr);


	}

	if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 8 
		|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 9
		|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 10
		|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 11
			&& es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 12)
	{

		std::cout << "hiya" << std::endl;
		delete velocity_matrix;
		std::cout << "hmmm" << std::endl;

		// also zero the assemble bits for next time
		system->request_matrix("Velocity Matrix")->close();
		system->request_matrix("Velocity Matrix")->zero();
	}


	if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 10 
		|| es->parameters.get<unsigned int>("preconditioner_type_3d") == 11
		|| es->parameters.get<unsigned int>("preconditioner_type_3d") == 12)
	{
		std::cout << "hiya deleting SIMPLE IS's" << std::endl;
		ierr = ISDestroy(&velocity_is); CHKERRQ(ierr);
		ierr = ISDestroy(&pressure_is); CHKERRQ(ierr);
		std::cout << "hiya deleted SIMPLE IS's" << std::endl;
	}

	
	std::cout << "hi" << std::endl;

	return final_linear_residual;
}

// ***** PCSHELL preconditioner
/*
   SampleShellPCCreate - This routine creates a user-defined
   preconditioner context.

   Output Parameter:
.  shell - user-defined preconditioner context
*/
PetscErrorCode ShellPCCreate(NSShellPC **shell)
{
  NSShellPC  *newctx;
  PetscErrorCode ierr;

  //ierr         = PetscNew(NSShellPC,&newctx);CHKERRQ(ierr);
	ierr         = PetscNew(&newctx);CHKERRQ(ierr);
  newctx->pressure_mass_matrix = 0;
  newctx->pressure_laplacian_matrix = 0;
	newctx->pressure_laplacian_preconditioner = 0;
  newctx->pressure_convection_diffusion_matrix = 0;
  newctx->velocity_matrix = 0;
  newctx->S_approx = 0;
  newctx->inner_mass_ksp = 0;
  newctx->inner_lap_ksp = 0;
  newctx->inner_velocity_ksp = 0;
  newctx->temp_vec = 0;
  newctx->temp_vec_2 = 0;
  *shell       = newctx;
  return 0;
}


// ***** PCSHELL preconditioner
/*
   SampleShellPCCreate - This routine creates a user-defined
   preconditioner context.

   Output Parameter:
.  shell - user-defined preconditioner context
*/
PetscErrorCode ShellPCCreate(SIMPLEShellPC **shell)
{
  SIMPLEShellPC  *newctx;
  PetscErrorCode ierr;

  //ierr         = PetscNew(NSShellPC,&newctx);CHKERRQ(ierr);
	ierr         = PetscNew(&newctx);CHKERRQ(ierr);
  *shell       = newctx;
  return 0;
}


/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "SampleShellPCSetUp"
/*
   SampleShellPCSetUp - This routine sets up a user-defined
   preconditioner context.

   Input Parameters:
.  pc    - preconditioner object
.  pmat  - preconditioner matrix
.  x     - vector

   Output Parameter:
.  shell - fully set up user-defined preconditioner context

   Notes:
   In this example, we define the shell preconditioner to be Jacobi's
   method.  Thus, here we create a work vector for storing the reciprocal
   of the diagonal of the preconditioner matrix; this vector is then
   used within the routine SampleShellPCApply().
*/
PetscErrorCode PressureShellPCSetUp(PC pc,Mat pressure_mass_matrix, KSP schur_ksp)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);
  //Vec            diag;

	//ierr = KSPDestroy(&shell->inner_lap_ksp);CHKERRQ(ierr);
	//ierr = VecDestroy(&shell->temp_vec);CHKERRQ(ierr);
	//ierr = VecDestroy(&shell->temp_vec_2);CHKERRQ(ierr);


  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SET MASS MATRIX ******************* //

  shell->pressure_mass_matrix = pressure_mass_matrix;

	// ********* SETUP MASS MATRIX KSP **************** //

	//ierr = KSPDestroy(&shell->inner_mass_ksp);CHKERRQ(ierr);
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_mass_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);
	ierr = KSPSetOperators(shell->inner_mass_ksp,Amat,Pmat); CHKERRQ(ierr);

	// setup up the operators for the first inner solve

	PC local_mass_pc;
	ierr = KSPGetPC(shell->inner_mass_ksp,&local_mass_pc); CHKERRQ(ierr);

	// set the pc to be a ksp
	// should probably use gmres by default okay
	ierr = PCSetType(local_mass_pc,PCKSP); CHKERRQ(ierr);
	KSP local_mass_ksp;

	ierr = PCKSPGetKSP(local_mass_pc,&local_mass_ksp); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(local_mass_ksp,"ns3d_fieldsplit_1_pc_mass_"); CHKERRQ(ierr);
	// setup up the operators for the first inner solve
	ierr = KSPSetOperators(local_mass_ksp,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);
	ierr = PCSetOperators(local_mass_pc,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (local_mass_ksp); CHKERRQ(ierr);

	

	ierr = PetscLogStagePop();

//	printf ("inside shell pc setup");
  return 0;
}

PetscErrorCode PCDShellPCSetUp(PC pc,Mat pressure_mass_matrix,Mat pressure_laplacian_matrix,Mat pressure_laplacian_preconditioner,Mat pressure_convection_diffusion_matrix, KSP schur_ksp)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);




  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SETUP MATRICES ****************** //

	shell->pressure_mass_matrix = pressure_mass_matrix;
  shell->pressure_laplacian_matrix = pressure_laplacian_matrix;
	shell->pressure_laplacian_preconditioner = pressure_laplacian_preconditioner;
  shell->pressure_convection_diffusion_matrix = pressure_convection_diffusion_matrix;


	// ********* SETUP LAPLACIAN KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_lap_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);
	ierr = KSPSetOperators(shell->inner_lap_ksp,Amat,Pmat); CHKERRQ(ierr);


	PC local_lap_pc;
	ierr = KSPGetPC(shell->inner_lap_ksp,&local_lap_pc); CHKERRQ(ierr);

	// set the pc to be a ksp
	// should probably use gmres by default okay
	ierr = PCSetType(local_lap_pc,PCKSP); CHKERRQ(ierr);


	KSP local_lap_ksp;
	ierr = PCKSPGetKSP(local_lap_pc,&local_lap_ksp); CHKERRQ(ierr);


	ierr = KSPSetOptionsPrefix(local_lap_ksp,"ns3d_fieldsplit_1_pc_lap_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (local_lap_ksp); CHKERRQ(ierr);
	

	// setup up the operators for the first inner solve
	ierr = KSPSetOperators(local_lap_ksp,shell->pressure_laplacian_matrix,shell->pressure_laplacian_matrix); CHKERRQ(ierr);
	ierr = PCSetOperators(local_lap_pc,shell->pressure_laplacian_matrix,shell->pressure_laplacian_matrix); CHKERRQ(ierr);

	


	// ********* SETUP MASS KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_mass_ksp); CHKERRQ(ierr);
	ierr = KSPSetOperators(shell->inner_mass_ksp,Amat,Pmat); CHKERRQ(ierr);

	PC local_mass_pc;
	ierr = KSPGetPC(shell->inner_mass_ksp,&local_mass_pc); CHKERRQ(ierr);
	ierr = PCSetType(local_mass_pc,PCKSP); CHKERRQ(ierr);
	KSP local_mass_ksp;
	ierr = PCKSPGetKSP(local_mass_pc,&local_mass_ksp); CHKERRQ(ierr);

	
	ierr = KSPSetOptionsPrefix(local_mass_ksp,"ns3d_fieldsplit_1_pc_mass_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (local_mass_ksp); CHKERRQ(ierr);
	
	ierr = KSPSetOperators(local_mass_ksp,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);
	ierr = PCSetOperators(local_mass_pc,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);




  // ******************** SETUP VECTORS ***************************** //
	PetscInt n_col;
	PetscInt n_col_local;
	ierr = MatGetSize(Amat,NULL,&n_col); CHKERRQ(ierr);
	ierr = MatGetLocalSize(Amat,NULL,&n_col_local); CHKERRQ(ierr);
	// create some vectors that can hold temporary info
	ierr = VecCreate(PETSC_COMM_WORLD,&shell->temp_vec); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&shell->temp_vec_2); CHKERRQ(ierr);
	ierr = VecSetSizes(shell->temp_vec, n_col_local, n_col); CHKERRQ(ierr);
	ierr = VecSetSizes(shell->temp_vec_2, n_col_local, n_col); CHKERRQ(ierr);
	ierr = VecSetFromOptions(shell->temp_vec); CHKERRQ(ierr);
	ierr = VecSetFromOptions(shell->temp_vec_2); CHKERRQ(ierr);


//	printf ("inside shell pc setup");

	ierr = PetscLogStagePop();
  return 0;
}


PetscErrorCode PCD2ShellPCSetUp(PC pc,Mat pressure_mass_matrix,Mat pressure_laplacian_matrix,Mat pressure_laplacian_preconditioner,Mat pressure_convection_diffusion_matrix, KSP schur_ksp)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);

  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	//create a KSP that can do the solving we require
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_mass_ksp); CHKERRQ(ierr);
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_lap_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	//MatStructure mat_structure;
	//ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat,&mat_structure); CHKERRQ(ierr);
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);

	ierr = KSPSetOperators(shell->inner_mass_ksp,Amat,Pmat); CHKERRQ(ierr);
	ierr = KSPSetOperators(shell->inner_lap_ksp,Amat,Pmat); CHKERRQ(ierr);


  // create some vectors for temporary pc applying
	PetscInt n_col;
	PetscInt n_col_local;
	ierr = MatGetSize(Amat,NULL,&n_col); CHKERRQ(ierr);
	ierr = MatGetLocalSize(Amat,NULL,&n_col_local); CHKERRQ(ierr);
	// create some vectors that can hold temporary info
	ierr = VecCreate(PETSC_COMM_WORLD,&shell->temp_vec); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&shell->temp_vec_2); CHKERRQ(ierr);
	ierr = VecSetSizes(shell->temp_vec, n_col_local, n_col); CHKERRQ(ierr);
	ierr = VecSetSizes(shell->temp_vec_2, n_col_local, n_col); CHKERRQ(ierr);
	ierr = VecSetFromOptions(shell->temp_vec); CHKERRQ(ierr);
	ierr = VecSetFromOptions(shell->temp_vec_2); CHKERRQ(ierr);


	shell->pressure_mass_matrix = pressure_mass_matrix;
  shell->pressure_laplacian_matrix = pressure_laplacian_matrix;
	shell->pressure_laplacian_preconditioner = pressure_laplacian_preconditioner;
  shell->pressure_convection_diffusion_matrix = pressure_convection_diffusion_matrix;

//	printf ("inside shell pc setup");
	ierr = PetscLogStagePop();
  return 0;
}


/*
   MonolithicShellPCSetUp - This routine sets up a user-defined
   preconditioner context.

   Input Parameters:
.  pc    - preconditioner object
.  pmat  - preconditioner matrix
.  x     - vector

   Output Parameter:
.  shell - fully set up user-defined preconditioner context

   Notes:
   In this example, we define the shell preconditioner to be Jacobi's
   method.  Thus, here we create a work vector for storing the reciprocal
   of the diagonal of the preconditioner matrix; this vector is then
   used within the routine SampleShellPCApply().
*/
PetscErrorCode MonolithicShellPCSetUp(PC pc,Mat velocity_matrix, KSP schur_ksp)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	//ierr = PetscLogStagePush(1); // not sure why this doesn't work anymore


  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SET THE MATRIX ******************* //

  shell->velocity_matrix = velocity_matrix;

	// ********* SETUP VELOCITY MATRIX KSP **************** //

	// create
	//ierr = KSPDestroy(&shell->inner_mass_ksp);CHKERRQ(ierr);
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_velocity_ksp); CHKERRQ(ierr);
	
	// set the operators for the velocity matrix inversion
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->velocity_matrix,shell->velocity_matrix); CHKERRQ(ierr);

	// setup up the PC for the velocity matrix inversion
	PC local_velocity_pc;
	ierr = KSPGetPC(shell->inner_velocity_ksp,&local_velocity_pc); CHKERRQ(ierr);

	// set the pc to be used, lu for now
	ierr = PCSetType(local_velocity_pc,PCLU); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(shell->inner_velocity_ksp,"ns3d1d_fieldsplit_1_pc_velocity_"); CHKERRQ(ierr);
	// setup up the operators for the first inner solve
	ierr = KSPSetFromOptions (shell->inner_velocity_ksp); CHKERRQ(ierr);



	// Some of the Mats we will extract and create
	Mat Bt;
	Mat B;
	Mat A11;
	Mat S;


	Mat Bt_dense;	// dense version of BT
	Mat T;
	Mat F;	// factored velocity matrix
	IS perm, iperm;
	MatFactorInfo info;
	Mat S_approx_dense;	// approximate dense schur complement

	ierr = KSPGetOperators(schur_ksp,&S,NULL); CHKERRQ(ierr);
	ierr = MatSchurComplementGetSubMatrices(S,NULL,NULL,&Bt,&B,&A11);
	printf ("hi\n");


	ierr = MatConvert(Bt,MATDENSE,MAT_INITIAL_MATRIX,&Bt_dense);
	ierr = MatDuplicate(Bt_dense,MAT_COPY_VALUES,&T);

	// create the factorisation
	printf ("2\n");
	
  ierr = MatGetOrdering(shell->velocity_matrix,  MATORDERINGNATURAL,  &perm,  &iperm); CHKERRQ(ierr);     
  ierr = MatFactorInfoInitialize(&info); CHKERRQ(ierr);
	
	// superlu doesn't work with matmatsolve, superlu_dist does.
	ierr = MatGetFactor(shell->velocity_matrix,MATSOLVERSUPERLU_DIST,MAT_FACTOR_LU,&F); CHKERRQ(ierr); 
	printf ("3\n");
	ierr = MatLUFactorSymbolic(F,shell->velocity_matrix,perm,iperm,&info); CHKERRQ(ierr); 
	printf ("4\n");
	ierr = MatLUFactorNumeric(F,shell->velocity_matrix,&info); CHKERRQ(ierr); 
	

	printf ("5\n");
	// if using matlufactor then use shell->velocity_matrix, is using matgetfactor use F
	ierr = MatMatSolve(F,Bt_dense,T);

	printf ("6\n");
	ierr = MatMatMult(B,T,MAT_INITIAL_MATRIX,1.0,&S_approx_dense);
	printf ("7\n");

	// convert S_approx to sparse
	ierr = MatConvert(S_approx_dense,MATAIJ,MAT_INITIAL_MATRIX,&shell->S_approx);

	// S_approx = A11 - S_approx
	ierr = MatAXPY(shell->S_approx,-1.0,A11,DIFFERENT_NONZERO_PATTERN);
	ierr = MatScale(shell->S_approx,-1.0);


	// only use A11, then we don't need the matscale below
	//ierr = MatDuplicate(A11,MAT_COPY_VALUES,&shell->S_approx);
	printf ("8\n");

	// setup up the operators for the 
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->S_approx,shell->S_approx); CHKERRQ(ierr);

	printf ("yeah\n");

	// delete the matrices we have created
	ierr = MatDestroy(&Bt_dense); CHKERRQ(ierr);
	ierr = MatDestroy(&T); CHKERRQ(ierr);
	
	ierr = MatDestroy(&F); CHKERRQ(ierr);
	ierr = ISDestroy(&perm); CHKERRQ(ierr);
	ierr = ISDestroy(&iperm); CHKERRQ(ierr);
	ierr = MatDestroy(&S_approx_dense); CHKERRQ(ierr);
	

	//ierr = PetscLogStagePop();

	printf ("inside shell pc monolithic velocity setup");
  return 0;
}


/*
   Monolithic2ShellPCSetUp - This solves the matrix in the schur complement vector by vector.
*/
PetscErrorCode Monolithic2ShellPCSetUp(PC pc,Mat velocity_matrix, Vec non_zero_cols, Vec non_zero_rows, KSP schur_ksp)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	//ierr = PetscLogStagePush(1);	// not sure why this log thing makes an error now...

  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SET MASS MATRIX ******************* //

  shell->velocity_matrix = velocity_matrix;

	// ********* SETUP VELOCITY MATRIX KSP **************** //

	// create
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_velocity_ksp); CHKERRQ(ierr);
	
	// set the operators for the velocity matrix inversion
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->velocity_matrix,shell->velocity_matrix); CHKERRQ(ierr);

	// setup up the PC for the velocity matrix inversion
	PC local_velocity_pc;
	ierr = KSPGetPC(shell->inner_velocity_ksp,&local_velocity_pc); CHKERRQ(ierr);

	// set the pc to be used, lu for now
	ierr = PCSetType(local_velocity_pc,PCLU); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(shell->inner_velocity_ksp,"ns3d1d_fieldsplit_1_pc_velocity_"); CHKERRQ(ierr);
	// setup up the operators for the first inner solve
	ierr = KSPSetFromOptions (shell->inner_velocity_ksp); CHKERRQ(ierr);


	// ************* Figure out how many 


	// Get the matrices we need from the system
	Mat Bt;
	Mat B;
	Mat A11;
	Mat S;

	ierr = KSPGetOperators(schur_ksp,&S,NULL); CHKERRQ(ierr);
	ierr = MatSchurComplementGetSubMatrices(S,NULL,NULL,&Bt,&B,&A11);

	// Create the Factor matrix object and associated items
	Mat F;	// factored velocity matrix
	IS perm, iperm;
	MatFactorInfo info;

	printf ("hi\n");


	printf ("2\n");
	
	// initialise the items for the factorisation (probably unnecessary)
  ierr = MatGetOrdering(shell->velocity_matrix,  MATORDERINGNATURAL,  &perm,  &iperm); CHKERRQ(ierr);     
  ierr = MatFactorInfoInitialize(&info); CHKERRQ(ierr);

	// Do the factorisation
	// superlu doesn't work with matmatsolve, superlu_dist does.
	// can use different arguments in here for the solver
	// 1) MATSOLVERPETSC - very slow
	// 2) MATSOLVERSUPERLU_DIST - works nicely
	ierr = MatGetFactor(shell->velocity_matrix,MATSOLVERSUPERLU_DIST,MAT_FACTOR_LU,&F); CHKERRQ(ierr); 
	printf ("3\n");
	ierr = MatLUFactorSymbolic(F,shell->velocity_matrix,perm,iperm,&info); CHKERRQ(ierr); 
	printf ("4\n");
	ierr = MatLUFactorNumeric(F,shell->velocity_matrix,&info); CHKERRQ(ierr); 
	

	printf ("5\n");

	// **** MANUALLY MAKE MATRIX BF^-1BT

	// instead of matmatsolve we want to do each column and then put it into a new sparse matrix T.
	PetscInt num_coupling_points;
	ierr = VecGetSize(non_zero_cols,&num_coupling_points);
	
	PetscInt m;
	PetscInt m_local;
	PetscInt n_1d;
	PetscInt n_1d_local;
	ierr = MatGetSize(Bt,&m,&n_1d);
	ierr = MatGetLocalSize(Bt,&m_local,&n_1d_local);
	Vec Bt_col;
	VecCreate(PETSC_COMM_WORLD,&Bt_col);
	VecSetSizes(Bt_col,m_local,m);
	VecSetFromOptions(Bt_col);
	Vec T_col;
	VecCreate(PETSC_COMM_WORLD,&T_col);
	VecSetSizes(T_col,m_local,m);
	VecSetFromOptions(T_col);

	PetscInt n;
	PetscInt n_local;
	ierr = MatGetSize(B,NULL,&n);
	ierr = MatGetLocalSize(B,NULL,&n_local);

	Vec B_row;
	VecCreate(PETSC_COMM_WORLD,&B_row);
	VecSetSizes(B_row,n_local,n);
	VecSetFromOptions(B_row);

	std::cout << "n = " << n << std::endl;
	std::cout << "m = " << m << std::endl;

	// Create a matrix, T_vals of the values that will go into the BF^-1BT matrix
	// T_rows are the what rows they go into (the flux coupling dofs)
	// T_cols are the what rows they go into (the pressure coupling dofs)
	PetscScalar T_vals[num_coupling_points][num_coupling_points];
	PetscInt T_rows[num_coupling_points];
	PetscInt T_cols[num_coupling_points];

	// Create a matrix to put the result into
	// We know that there won't be more than num_coupling_points on each row,
	// so use this as an estimate.
	Mat new_T;
	MatCreateAIJ(PETSC_COMM_WORLD,n_1d_local,n_1d_local,n_1d,n_1d,
                       num_coupling_points,NULL,num_coupling_points,NULL,&new_T);

	// loop over the columns
	for(unsigned int i=0; i<num_coupling_points; i++)
	{
		PetscScalar col_number[1];
		PetscInt coupling_point[1] = {i};	// the coupling point we are interested in (must be an array...)
		ierr = VecGetValues(non_zero_cols,1,coupling_point,col_number);

		// extract the correct column as a Vec
		ierr = MatGetColumnVector(Bt,Bt_col,col_number[0]);

		// copmute F^-1BT for this particular row
		ierr = MatSolve(F,Bt_col,T_col);

		// loop over the rows for the other bit of the matrix multiplication		
		for(unsigned int j=0; j<num_coupling_points; j++)
		{
			PetscScalar row_number[1];
			coupling_point[0] = j;			// the coupling point we are interested in (must be an array...)
			ierr = VecGetValues(non_zero_rows,1,coupling_point,row_number);		// get the actual dof index in the B matrix
		
			// create some arrays for the values and indices in the rows
			const PetscScalar *B_row_array;
			const PetscInt *B_row_cols;
			PetscInt B_row_n_cols;

			// extract the correct row
			ierr = MatGetRow(B,row_number[0],&B_row_n_cols,&B_row_cols,&B_row_array);
			VecSet(B_row,0);	// reinitialise each time

			// convert the row values to a Vec for dot product
			ierr = VecSetValues(B_row,B_row_n_cols,B_row_cols,B_row_array,INSERT_VALUES);

			// do a dot product with the column that was create in the outer loop
			PetscScalar T_val;
			ierr = VecDot(B_row,T_col,&T_val);

			// save the result in the correct place and save the row and column indices
			T_vals[i][j] = T_val;
			T_cols[i] = col_number[0];
			T_rows[j] = row_number[0];

			
		}
		
	}

	// set the values of the matrix and assemble them
	MatSetValues(new_T,num_coupling_points,T_rows,num_coupling_points,T_cols,(const PetscScalar*)T_vals,INSERT_VALUES);
	MatAssemblyBegin(new_T,MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd(new_T,MAT_FINAL_ASSEMBLY );

	// put the values into the correct matrix... for no good reason..
	ierr = MatDuplicate(new_T,MAT_COPY_VALUES,&shell->S_approx);
	
	// S_approx = A11 - S_approx
	ierr = MatAXPY(shell->S_approx,-1.0,A11,DIFFERENT_NONZERO_PATTERN);
	ierr = MatScale(shell->S_approx,-1.0);

	printf ("8\n");

	// setup up the operators for the solve later
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->S_approx,shell->S_approx); CHKERRQ(ierr);

	printf ("yeah\n");

	// delete the matrices we have created	
	ierr = MatDestroy(&F); CHKERRQ(ierr);
	ierr = ISDestroy(&perm); CHKERRQ(ierr);
	ierr = ISDestroy(&iperm); CHKERRQ(ierr);
	

	//ierr = PetscLogStagePop();

	printf ("inside shell pc monolithic 2 velocity setup");
  return 0;
}






PetscErrorCode LSCShellPCSetUp(PC pc, KSP schur_ksp)
{
	printf ("inside shell pc lsc james setup");
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);




  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SETUP LAPLACIAN KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_lap_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_lap_ksp,"ns3d_fieldsplit_1_pc_lap_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_lap_ksp); CHKERRQ(ierr);
	



	
	// ********** SETUP THE INTERMEDIATE VECS ***************	//
	Mat A,Ap,B,C;
	MatSchurComplementGetSubMatrices(Amat,&A,&Ap,&B,&C,NULL);
	MatGetVecs(A,&shell->temp_vec,&shell->temp_vec_2);
	MatGetVecs(Amat,&shell->temp_vec_3,NULL);
	VecDuplicate(shell->temp_vec,&shell->lsc_scale);

	// ********** SETUP L = A10 * A01 *************************	//
	MatMatMult(C,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&shell->lsc_laplacian_matrix);		// could possibly use MAT_REUSE_MATRIX

	// ********** SETUP DIAGONAL ******************* //
	MatGetDiagonal(Ap,shell->lsc_scale); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(shell->lsc_scale);
	
	// setup up the operators for the laplacian solve
	ierr = KSPSetOperators(shell->inner_lap_ksp,shell->lsc_laplacian_matrix,shell->lsc_laplacian_matrix); CHKERRQ(ierr);


	printf ("done shell pc lsc james setup");

	ierr = PetscLogStagePop();
  return 0;
}







PetscErrorCode LSCScaledShellPCSetUp(PC pc, Mat velocity_mass_matrix, KSP schur_ksp)
{
	printf ("inside shell pc lsc james setup");
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);




  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SETUP LAPLACIAN KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_lap_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_lap_ksp,"ns3d_fieldsplit_1_pc_lap_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_lap_ksp); CHKERRQ(ierr);
	

	// ********* SET THE VELOCITY MASS MATRIX ****************	//
  	shell->velocity_matrix = velocity_mass_matrix;


	
	// ********** SETUP THE INTERMEDIATE VECS ***************	//
	Mat A,Ap,B,C,B_scaled;
	MatSchurComplementGetSubMatrices(Amat,&A,&Ap,&B,&C,NULL);
	MatGetVecs(A,&shell->temp_vec,&shell->temp_vec_2);
	MatGetVecs(Amat,&shell->temp_vec_3,NULL);
	VecDuplicate(shell->temp_vec,&shell->lsc_scale);
	MatDuplicate(B,MAT_COPY_VALUES,&B_scaled);	// create temporary place for B_scaled

	// ********** SETUP DIAGONAL ******************* //
	MatGetDiagonal(shell->velocity_matrix,shell->lsc_scale); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(shell->lsc_scale);
	
	// ********** SETUP L = A10* Q_V^-1* A01 *************************	//
	MatDiagonalScale(B_scaled,shell->lsc_scale,NULL);
	MatMatMult(C,B_scaled,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&shell->lsc_laplacian_matrix);		// could possibly use MAT_REUSE_MATRIX

	// setup up the operators for the laplacian solve
	ierr = KSPSetOperators(shell->inner_lap_ksp,shell->lsc_laplacian_matrix,shell->lsc_laplacian_matrix); CHKERRQ(ierr);

	std::cout << "Bt" << std::endl;
	MatView(B,PETSC_VIEWER_STDOUT_SELF);

	std::cout << "B" << std::endl;
	MatView(C,PETSC_VIEWER_STDOUT_SELF);


	std::cout << "lsc_laplacian after" << std::endl;
	MatView(shell->lsc_laplacian_matrix,PETSC_VIEWER_STDOUT_SELF);

	MatDestroy(&B_scaled);
	printf ("done shell pc lsc james setup");

	ierr = PetscLogStagePop();
  return 0;
}






PetscErrorCode LSCScaledStabilisedShellPCSetUp(PC pc, Mat velocity_mass_matrix, KSP schur_ksp)
{
	std::cout << "hello?" << std::endl;
	printf ("inside shell pc lsc scaled stabilised james setup");
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);



	//std::cout << "hello?" << std::endl;

  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SETUP LAPLACIAN KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_lap_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_lap_ksp,"ns3d_fieldsplit_1_pc_lap_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_lap_ksp); CHKERRQ(ierr);
	

	// ********* SET THE VELOCITY MASS MATRIX ****************	//
  	shell->velocity_matrix = velocity_mass_matrix;

	//std::cout << "hello?" << std::endl;

	
	// ********** SETUP THE INTERMEDIATE VECS ***************	//
	Mat A,Ap,Bt,B,Bt_scaled,C;
	MatSchurComplementGetSubMatrices(Amat,&A,&Ap,&Bt,&B,&C);
	MatGetVecs(A,&shell->temp_vec,&shell->temp_vec_2);
	MatGetVecs(Amat,&shell->temp_vec_3,NULL);
	VecDuplicate(shell->temp_vec,&shell->lsc_scale);
	MatDuplicate(Bt,MAT_COPY_VALUES,&Bt_scaled);	// create temporary place for B_scaled

	//std::cout << "hello?" << std::endl;
	// ********** SETUP VELOCITY MASS DIAGONAL INVERSE ******************* //
	MatGetDiagonal(shell->velocity_matrix,shell->lsc_scale); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(shell->lsc_scale);
	
	//std::cout << "hello?" << std::endl;
	// ********** CALCULATE D_inv = inv(diag(B*diag(F)^-1*B^T - C)) ************ //
	Mat D;	// D matrix
	Mat D_temp;	// D matrix
	Mat B_F_inv_Bt_D_inv;	// D matrix
	Vec F_diag_inv;
	//std::cout << "hello?" << std::endl;
	MatDuplicate(Bt,MAT_COPY_VALUES,&D_temp);	// D = B^T
	//std::cout << "hello?" << std::endl;
	MatGetVecs(A,&F_diag_inv,NULL);
	MatGetDiagonal(A,F_diag_inv); 			// get the f_diag
	//std::cout << "hello?" << std::endl;
	VecReciprocal(F_diag_inv);			// invert it
	//std::cout << "hello?" << std::endl;
	MatDiagonalScale(D_temp,F_diag_inv,NULL);	// D = diag(F)^-1 * B^T
	//std::cout << "hello?" << std::endl;
	MatMatMult(B,D_temp,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&D);	// D = B * diag(F)^-1 * B^T
	//std::cout << "hello?" << std::endl;
	MatDuplicate(D,MAT_COPY_VALUES,&B_F_inv_Bt_D_inv);	// create temporary place for B_scaled
	//std::cout << "hello?" << std::endl;
	MatAXPY(D,-1.0,C,DIFFERENT_NONZERO_PATTERN);	// D = B * diag(F)^-1 * B^T - C
	//std::cout << "hello?" << std::endl;
	MatGetVecs(D,&shell->lsc_stab_alpha_D_inv,NULL);
	MatGetDiagonal(D,shell->lsc_stab_alpha_D_inv); // D = diag(B * diag(F)^-1 * B^T - C)
	//std::cout << "hello?" << std::endl;
	VecReciprocal(shell->lsc_stab_alpha_D_inv);	// invert it

	// make B_F_inv_Bt_D_inv
	MatDiagonalScale(B_F_inv_Bt_D_inv,NULL,shell->lsc_stab_alpha_D_inv);	// B_F_inv_Bt_D_inv (the shell var doesn't have alpha yet...)
	
	
	//std::cout << "hello?" << std::endl;

	// ********** SETUP L = A10* Q_V^-1* A01 *************************	//
	Vec B_Q_inv_Bt_diag;
	Vec C_diag_inv;
	Vec D_r_vec_sqrt;
	MatDiagonalScale(Bt_scaled,shell->lsc_scale,NULL);
	MatMatMult(B,Bt_scaled,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&shell->lsc_laplacian_matrix);		// could possibly use MAT_REUSE_MATRIX

	//std::cout << "lsc_laplacian" << std::endl;
	//MatView(shell->lsc_laplacian_matrix,PETSC_VIEWER_STDOUT_SELF);
	
	
	
	//std::cout << "hello?" << std::endl;

	// ********* CALCULATE D_r_vec ***************** //
	PetscReal D_r_infty;
	//std::cout << "hi?" << std::endl;
	MatGetVecs(shell->lsc_laplacian_matrix,&B_Q_inv_Bt_diag,NULL);
	MatGetDiagonal(shell->lsc_laplacian_matrix,B_Q_inv_Bt_diag); // 
	MatGetVecs(C,&C_diag_inv,NULL);
	MatGetVecs(C,&D_r_vec_sqrt,NULL);
	MatGetDiagonal(C,C_diag_inv); // 
	VecReciprocal(C_diag_inv);	// invert it
	VecPointwiseMult(D_r_vec_sqrt,B_Q_inv_Bt_diag,C_diag_inv);
	VecNorm(D_r_vec_sqrt,NORM_INFINITY,&D_r_infty);
	VecSqrtAbs(D_r_vec_sqrt);


	//std::cout << "hello?" << std::endl;
	// we need to calculate the spectral radius, largest eigenvalue of Q_inv * F
	Mat Q_inv_F;
	MatDuplicate(A,MAT_COPY_VALUES,&Q_inv_F);	// create temporary place for B_scaled
	MatDiagonalScale(Q_inv_F,shell->lsc_scale,NULL);	// D = diag(F)^-1 * B^T

	//std::cout << "hello?" << std::endl;
	// do power series
	Vec b_k;
	Vec b_k_plus_1;
	// initialise b_k
	MatGetVecs(A,&b_k,NULL);
	MatGetVecs(A,&b_k_plus_1,NULL);
	//std::cout << "hi?" << std::endl;
	MatGetDiagonal(A,b_k); // get the f_diag  
	PetscRandom rctx;   
	PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
	PetscRandomSetFromOptions(rctx);
     	VecSetRandom(b_k,rctx);	
	//std::cout << "hi?" << std::endl;
	//std::cout << "hi?" << std::endl;
	PetscReal norm = 1.0;
	PetscReal max_eigenvalue_Q_inv_F = 1.0;
	PetscReal max_eigenvalue_Q_inv_F_prev = 1.0;
	
	//std::cout << "hello?" << std::endl;
	
	for(int i=0; i<100; i++)
	{
		MatMult(Q_inv_F,b_k,b_k_plus_1);	//Ab_k = b_k+1		

		max_eigenvalue_Q_inv_F_prev = max_eigenvalue_Q_inv_F;
		// eigenvalue estimate
		VecDot(b_k,b_k_plus_1,&max_eigenvalue_Q_inv_F); // mu = b_k*A*b_k/b_k*b_k
		//max_eigenvalue = max_eigenvalue/norm;
		//std::cout << "1/norm " << i << ": " << norm << std::endl;

		// calc new norm
		VecNorm(b_k_plus_1,NORM_2,&norm);
		VecScale(b_k_plus_1,1./norm);
		VecCopy(b_k_plus_1,b_k);
	}

	std::cout << "max eigenvalue Q inv F : " << max_eigenvalue_Q_inv_F << std::endl;
	std::cout << "error : " << fabs(max_eigenvalue_Q_inv_F - max_eigenvalue_Q_inv_F_prev)/max_eigenvalue_Q_inv_F << std::endl;


	//std::cout << "hello?" << std::endl;
	

	MatGetVecs(B_F_inv_Bt_D_inv,&b_k,NULL);
	MatGetVecs(B_F_inv_Bt_D_inv,&b_k_plus_1,NULL);
     	VecSetRandom(b_k,rctx);	

	PetscReal max_eigenvalue_B_F_inv_Bt_D_inv = 1.0;
	PetscReal max_eigenvalue_B_F_inv_Bt_D_inv_prev = 1.0;
	norm = 1.0;

	for(int i=0; i<100; i++)
	{
		MatMult(B_F_inv_Bt_D_inv,b_k,b_k_plus_1);	//Ab_k = b_k+1		

		max_eigenvalue_B_F_inv_Bt_D_inv_prev = max_eigenvalue_B_F_inv_Bt_D_inv;
		// eigenvalue estimate
		VecDot(b_k,b_k_plus_1,&max_eigenvalue_B_F_inv_Bt_D_inv); // mu = b_k*A*b_k/b_k*b_k
		//max_eigenvalue = max_eigenvalue/norm;
		//std::cout << "1/norm " << i << ": " << norm << std::endl;

		// calc new norm
		VecNorm(b_k_plus_1,NORM_2,&norm);
		VecScale(b_k_plus_1,1./norm);
		VecCopy(b_k_plus_1,b_k);
	}

	std::cout << "max eigenvalue B_F_inv_Bt_D_inv : " << max_eigenvalue_B_F_inv_Bt_D_inv << std::endl;
	std::cout << "error : " << fabs(max_eigenvalue_B_F_inv_Bt_D_inv - max_eigenvalue_B_F_inv_Bt_D_inv_prev)/max_eigenvalue_B_F_inv_Bt_D_inv << std::endl;


	
	//std::cout << "hello?" << std::endl;	
 	VecDestroy(&b_k);
	VecDestroy(&b_k_plus_1);


	//std::cout << "hello?" << std::endl;
	// ************ Set gamma and alpha ************** //
	PetscReal gamma = max_eigenvalue_Q_inv_F/3.0;
	//PetscReal gamma = 1.0e+2/3.0;
	gamma = gamma/D_r_infty;
	PetscReal alpha = 1.0/max_eigenvalue_B_F_inv_Bt_D_inv;

	//std::cout << "hello?" << std::endl;

	// ********** SETUP L = A10* Q_V^-1* A01 - gamma*D_r_sqrt* C*D_r_sqrt *************************	//
	Mat D_r_sqrt_C_D_r_sqrt;
	MatDuplicate(C,MAT_COPY_VALUES,&D_r_sqrt_C_D_r_sqrt);	// create temporary place for B_scaled
	MatDiagonalScale(D_r_sqrt_C_D_r_sqrt,D_r_vec_sqrt,D_r_vec_sqrt);	// D = diag(F)^-1 * B^T
	MatAXPY(shell->lsc_laplacian_matrix,-gamma,D_r_sqrt_C_D_r_sqrt,DIFFERENT_NONZERO_PATTERN);		// could possibly use MAT_REUSE_MATRIX
	
	//std::cout << "hello?" << std::endl;
	// ********** SETUP \alpha D_inv
	VecScale(shell->lsc_stab_alpha_D_inv,alpha);
	


	

	//std::cout << "hello?" << std::endl;
	// setup up the operators for the laplacian solve
	ierr = KSPSetOperators(shell->inner_lap_ksp,shell->lsc_laplacian_matrix,shell->lsc_laplacian_matrix); CHKERRQ(ierr);

	//std::cout << "lsc_laplacian after" << std::endl;
	//MatView(shell->lsc_laplacian_matrix,PETSC_VIEWER_STDOUT_SELF);

	// *********** DESTROY ALL THE TEMPORARY THINGS VECs and MATs ***********	//
	MatDestroy(&D);	// D matrix
	MatDestroy(&D_temp);	// D matrix
	MatDestroy(&B_F_inv_Bt_D_inv);	// D matrix
	VecDestroy(&F_diag_inv);
	VecDestroy(&B_Q_inv_Bt_diag);
	VecDestroy(&C_diag_inv);
	VecDestroy(&D_r_vec_sqrt);
	MatDestroy(&Q_inv_F);
	MatDestroy(&D_r_sqrt_C_D_r_sqrt);
	MatDestroy(&Bt_scaled);
     	PetscRandomDestroy(&rctx);
	printf ("done shell pc lsc james setup");

	//std::cout << "hello?" << std::endl;
	ierr = PetscLogStagePop();
  return 0;
}







PetscErrorCode SIMPLEShellPCSetUp(PC pc, IS velocity_is, IS pressure_is, KSP schur_ksp)
{
	printf ("inside shell pc SIMPLE james setup");
  SIMPLEShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);




  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);


	// ************* WHAT DO WE NEED ****************** //
	// we need a ksp for the velocity solve
	// we need a ksp for the schur solve


	// ************** CREATE THE SUB MATRICES WE WILL USE ********* //
	shell->velocity_is = velocity_is;
	shell->pressure_is = pressure_is;

	Mat A;
	ierr = PCGetOperators(pc,&A,NULL); CHKERRQ(ierr);

  	MatGetSubMatrix(A,shell->velocity_is,shell->velocity_is,MAT_INITIAL_MATRIX,&shell->A);
  	MatGetSubMatrix(A,shell->velocity_is,shell->pressure_is,MAT_INITIAL_MATRIX,&shell->Bt);
  	MatGetSubMatrix(A,shell->pressure_is,shell->velocity_is,MAT_INITIAL_MATRIX,&shell->B);
  	MatGetSubMatrix(A,shell->pressure_is,shell->pressure_is,MAT_INITIAL_MATRIX,&shell->C);


	// ********* SETUP VELOCITY KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_velocity_ksp); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_velocity_ksp,"ns3d_fieldsplit_1_pc_velocity_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_velocity_ksp); CHKERRQ(ierr);

	// setup up the operators for the velocity solve
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->A,shell->A); CHKERRQ(ierr);
	
	ierr = KSPSetType(shell->inner_velocity_ksp, KSPPREONLY); CHKERRQ(ierr);

	PC local_velocity_pc;
	ierr = KSPGetPC(shell->inner_velocity_ksp,&local_velocity_pc); CHKERRQ(ierr);

	ierr = PCSetType(local_velocity_pc,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverPackage(local_velocity_pc,MATSOLVERSUPERLU_DIST);

	// ************ SETUP SCHUR MATRIX ******************** //
	
	// ********** SETUP VELOCITY  DIAGONAL INVERSE ******************* //
	Vec F_diag_inv;
	MatGetVecs(shell->A,&F_diag_inv,NULL);
	MatGetDiagonal(shell->A,F_diag_inv); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(F_diag_inv);


	// ********** SETUP S_approx = C -  B* diag(F)^-1* B^T *************************	//
	Mat Bt_scaled;
	Mat B_diag_F_inv_Bt;
	MatDuplicate(shell->Bt,MAT_COPY_VALUES,&Bt_scaled);	// create temporary place for Bt_scaled
	MatDuplicate(shell->C,MAT_COPY_VALUES,&B_diag_F_inv_Bt);	// create sapce for B_diag_F_inv_Bt
	MatDuplicate(shell->C,MAT_COPY_VALUES,&shell->S_approx);	// set values of S_approx
	MatDiagonalScale(Bt_scaled,F_diag_inv,NULL);

	MatMatMult(shell->B,Bt_scaled,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&B_diag_F_inv_Bt);		// could possibly use MAT_REUSE_MATRIX


	MatAXPY(shell->S_approx,-1.0,B_diag_F_inv_Bt,DIFFERENT_NONZERO_PATTERN);		// could possibly use MAT_REUSE_MATRIX
	//MatScale(shell->S_approx,-1.0);		// could possibly use MAT_REUSE_MATRIX


	// ********** SETUP TEMPORARY VECTORS *********************** //
	// shell->temp_vec is u_star
	// shell->temp_vec_2 is delta_p
	// shell->temp_vec_3 is Bt*delta_p
	MatGetVecs(shell->A,&shell->temp_vec,&shell->temp_vec_3);
	MatGetVecs(shell->C,&shell->temp_vec_2,NULL);
	
	

	//ierr = MatDuplicate(shell->C,MAT_COPY_VALUES,&shell->S_approx);

	// ********* SETUP SCHUR KSP ****************	//
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_schur_ksp); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_schur_ksp,"ns3d_fieldsplit_1_pc_schur_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_schur_ksp); CHKERRQ(ierr);

	// setup up the operators for the schur solve
	ierr = KSPSetOperators(shell->inner_schur_ksp,shell->S_approx,shell->S_approx); CHKERRQ(ierr);
	
	//ierr = KSPSetType(shell->inner_schur_ksp, KSPPREONLY); CHKERRQ(ierr);
	ierr = KSPSetType(shell->inner_schur_ksp, KSPGMRES); CHKERRQ(ierr);

	PC local_schur_pc;
	ierr = KSPGetPC(shell->inner_schur_ksp,&local_schur_pc); CHKERRQ(ierr);

	ierr = PCSetType(local_schur_pc,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverPackage(local_schur_pc,MATSOLVERSUPERLU_DIST);
	//ierr = PCSetType(local_schur_pc,PCHYPRE); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverPackage(local_schur_pc,MATSOLVERMUMPS);



	// ************** DESTROY TEMP VECS AND MATS ********************* //
	VecDestroy(&F_diag_inv);
	MatDestroy(&Bt_scaled);
	MatDestroy(&B_diag_F_inv_Bt);
	printf ("done shell pc SIMPLE james setup");

	ierr = PetscLogStagePop();
  return 0;
}





PetscErrorCode SIMPLECShellPCSetUp(PC pc, IS velocity_is, IS pressure_is, KSP schur_ksp)
{
	printf ("inside shell pc SIMPLE james setup");
  SIMPLEShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);




  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);


	// ************* WHAT DO WE NEED ****************** //
	// we need a ksp for the velocity solve
	// we need a ksp for the schur solve


	// ************** CREATE THE SUB MATRICES WE WILL USE ********* //
	shell->velocity_is = velocity_is;
	shell->pressure_is = pressure_is;

	Mat A;
	ierr = PCGetOperators(pc,&A,NULL); CHKERRQ(ierr);

  	MatGetSubMatrix(A,shell->velocity_is,shell->velocity_is,MAT_INITIAL_MATRIX,&shell->A);
  	MatGetSubMatrix(A,shell->velocity_is,shell->pressure_is,MAT_INITIAL_MATRIX,&shell->Bt);
  	MatGetSubMatrix(A,shell->pressure_is,shell->velocity_is,MAT_INITIAL_MATRIX,&shell->B);
  	MatGetSubMatrix(A,shell->pressure_is,shell->pressure_is,MAT_INITIAL_MATRIX,&shell->C);


	// ********* SETUP VELOCITY KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_velocity_ksp); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_velocity_ksp,"ns3d_fieldsplit_1_pc_velocity_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_velocity_ksp); CHKERRQ(ierr);

	// setup up the operators for the velocity solve
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->A,shell->A); CHKERRQ(ierr);
	
	ierr = KSPSetType(shell->inner_velocity_ksp, KSPPREONLY); CHKERRQ(ierr);

	PC local_velocity_pc;
	ierr = KSPGetPC(shell->inner_velocity_ksp,&local_velocity_pc); CHKERRQ(ierr);

	ierr = PCSetType(local_velocity_pc,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverPackage(local_velocity_pc,MATSOLVERSUPERLU_DIST);

	// ************ SETUP SCHUR MATRIX ******************** //
	
	// ********** SETUP VELOCITY  DIAGONAL INVERSE ******************* //
	Vec F_diag_inv;
	MatGetVecs(shell->A,&F_diag_inv,NULL);
	//MatGetDiagonal(shell->A,F_diag_inv); /* Should be the mass matrix, but we don't have plumbing for that yet */
	//VecReciprocal(F_diag_inv);

	// multiply by transpose, extract diagonal, square root diagonal


	// instead of the diagonal we use the absolute value of the rowsum, this doesn't work so try sum of abs values
	PetscInt          rstart,rend;
	PetscInt          row,ncols,j,nrows;
	PetscReal         val;
	const PetscInt    *cols;
	const PetscScalar *vals;
	MatGetOwnershipRange(shell->A,&rstart,&rend);
	nrows = 0;
	for (row=rstart; row<rend; row++) {
		MatGetRow(shell->A,row,&ncols,&cols,&vals);
		val  = 0.0;
		for (j=0; j<ncols; j++) {
			val += PetscAbsScalar(vals[j]);
		}
		VecSetValue(F_diag_inv,row,val,ADD_VALUES);

		MatRestoreRow(A,row,&ncols,&cols,&vals);
	}
	VecAssemblyBegin(F_diag_inv);
	VecAssemblyEnd(F_diag_inv);

	//std::cout << "F_diag  sum" << std::endl;
	//VecView(F_diag_inv,PETSC_VIEWER_STDOUT_SELF);

	VecReciprocal(F_diag_inv);

	//std::cout << "F_diag_inv sum" << std::endl;
	//VecView(F_diag_inv,PETSC_VIEWER_STDOUT_SELF);


	//std::cout << "A" << std::endl;
	//MatView(shell->A,PETSC_VIEWER_STDOUT_SELF);



	// ********** SETUP S_approx = C -  B* diag(F)^-1* B^T *************************	//
	Mat Bt_scaled;
	Mat B_diag_F_inv_Bt;
	MatDuplicate(shell->Bt,MAT_COPY_VALUES,&Bt_scaled);	// create temporary place for Bt_scaled
	MatDuplicate(shell->C,MAT_COPY_VALUES,&B_diag_F_inv_Bt);	// create sapce for B_diag_F_inv_Bt
	MatDuplicate(shell->C,MAT_COPY_VALUES,&shell->S_approx);	// set values of S_approx
	MatDiagonalScale(Bt_scaled,F_diag_inv,NULL);

	MatMatMult(shell->B,Bt_scaled,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&B_diag_F_inv_Bt);		// could possibly use MAT_REUSE_MATRIX


	MatAXPY(shell->S_approx,-1.0,B_diag_F_inv_Bt,DIFFERENT_NONZERO_PATTERN);		// could possibly use MAT_REUSE_MATRIX
	//MatScale(shell->S_approx,-1.0);		// could possibly use MAT_REUSE_MATRIX


	// ********** SETUP TEMPORARY VECTORS *********************** //
	// shell->temp_vec is u_star
	// shell->temp_vec_2 is delta_p
	// shell->temp_vec_3 is Bt*delta_p
	MatGetVecs(shell->A,&shell->temp_vec,&shell->temp_vec_3);
	MatGetVecs(shell->C,&shell->temp_vec_2,NULL);
	
	

	//ierr = MatDuplicate(shell->C,MAT_COPY_VALUES,&shell->S_approx);

	// ********* SETUP SCHUR KSP ****************	//
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_schur_ksp); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_schur_ksp,"ns3d_fieldsplit_1_pc_schur_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_schur_ksp); CHKERRQ(ierr);

	// setup up the operators for the schur solve
	ierr = KSPSetOperators(shell->inner_schur_ksp,shell->S_approx,shell->S_approx); CHKERRQ(ierr);
	
	//ierr = KSPSetType(shell->inner_schur_ksp, KSPPREONLY); CHKERRQ(ierr);
	ierr = KSPSetType(shell->inner_schur_ksp, KSPGMRES); CHKERRQ(ierr);

	PC local_schur_pc;
	ierr = KSPGetPC(shell->inner_schur_ksp,&local_schur_pc); CHKERRQ(ierr);

	ierr = PCSetType(local_schur_pc,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverPackage(local_schur_pc,MATSOLVERSUPERLU_DIST);
	//ierr = PCSetType(local_schur_pc,PCHYPRE); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverPackage(local_schur_pc,MATSOLVERMUMPS);



	// ************** DESTROY TEMP VECS AND MATS ********************* //
	VecDestroy(&F_diag_inv);
	MatDestroy(&Bt_scaled);
	MatDestroy(&B_diag_F_inv_Bt);
	printf ("done shell pc SIMPLE james setup");

	ierr = PetscLogStagePop();
  return 0;
}












/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "SampleShellPCApply"
/*
   SampleShellPCApply - This routine demonstrates the use of a
   user-provided preconditioner.

   Input Parameters:
+  pc - preconditioner object
-  x - input vector

   Output Parameter:
.  y - preconditioned vector

   Notes:
   This code implements the Jacobi preconditioner, merely as an
   example of working with a PCSHELL.  Note that the Jacobi method
   is already provided within PETSc.
*/

PetscErrorCode PressureShellPCApply(PC pc,Vec x,Vec y)
{

  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);

	//printf ("in Pressure preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	//VecView(x,PETSC_VIEWER_STDOUT_SELF);

	// setup up the operators for the first inner solve

	PC local_mass_pc;
	ierr = KSPGetPC(shell->inner_mass_ksp,&local_mass_pc); CHKERRQ(ierr);
	KSP local_mass_ksp;
	ierr = PCKSPGetKSP(local_mass_pc,&local_mass_ksp); CHKERRQ(ierr);


	ierr = PCApply(local_mass_pc,x,y); CHKERRQ(ierr);
	PetscInt num_mass_its = 0;
  ierr = KSPGetIterationNumber (local_mass_ksp, &num_mass_its);
	//printf ("num_mass_its = %d\n",num_mass_its);


	//ierr = MatMult(shell->pressure_mass_matrix,x,y);

	//VecView(y,PETSC_VIEWER_STDOUT_SELF);

	ierr = PetscLogStagePop();
  return 0;
}


PetscErrorCode PCDShellPCApply(PC pc,Vec x,Vec y)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	PetscLogStage pcd_stage;
	ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	//MatView(pressure_mass_matrix,PETSC_VIEWER_STDOUT_SELF);

	// get the operators for the first inner solve

	PC local_lap_pc;
	ierr = KSPGetPC(shell->inner_lap_ksp,&local_lap_pc); CHKERRQ(ierr);
	KSP local_lap_ksp;
	ierr = PCKSPGetKSP(local_lap_pc,&local_lap_ksp); CHKERRQ(ierr);
	
	ierr = PCApply(local_lap_pc,x,shell->temp_vec);CHKERRQ(ierr);
	PetscInt num_lap_its = 0;
  ierr = KSPGetIterationNumber (local_lap_ksp, &num_lap_its); CHKERRQ(ierr);
	//std::cout << "num_lap_its = " << num_lap_its << std::endl;
	//printf ("num_lap_its = %d\n",num_lap_its);

	PetscInt rows;
	PetscInt cols;
	ierr = MatGetSize(shell->pressure_laplacian_matrix,&rows,&cols);




	// ******* Apply F_p ************ //
	//
	ierr = MatMult(shell->pressure_convection_diffusion_matrix,shell->temp_vec,shell->temp_vec_2); CHKERRQ(ierr);
	//ierr = MatMult(shell->pressure_convection_diffusion_matrix,x,shell->temp_vec_2); CHKERRQ(ierr);

	
	

	// ******* Apply M_p^-1 ********* //

	PC local_mass_pc;
	ierr = KSPGetPC(shell->inner_mass_ksp,&local_mass_pc); CHKERRQ(ierr);
	KSP local_mass_ksp;
	ierr = PCKSPGetKSP(local_mass_pc,&local_mass_ksp); CHKERRQ(ierr);

	ierr = PCApply(local_mass_pc,shell->temp_vec_2,y); CHKERRQ(ierr);
	PetscInt num_mass_its = 0;
  ierr = KSPGetIterationNumber (local_mass_ksp, &num_mass_its); CHKERRQ(ierr);
	//std::cout << "num_mass_its = " << num_mass_its << std::endl;
	//printf ("num_mass_its = %d\n",num_mass_its);

	ierr = PetscLogStagePop();

  return 0;
}

PetscErrorCode PCD2ShellPCApply(PC pc,Vec x,Vec y)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	PetscLogStage pcd2_stage;
	ierr = PetscLogStagePush(1);

	printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	//MatView(pressure_mass_matrix,PETSC_VIEWER_STDOUT_SELF);

	// setup up the operators for the first inner solve
	PC local_mass_pc;
	ierr = KSPGetPC(shell->inner_mass_ksp,&local_mass_pc); CHKERRQ(ierr);

	// set the pc to be a ksp
	// should probably use gmres by default okay
	ierr = PCSetType(local_mass_pc,PCKSP); CHKERRQ(ierr);


	KSP local_mass_ksp;
	ierr = PCKSPGetKSP(local_mass_pc,&local_mass_ksp); CHKERRQ(ierr);


	ierr = KSPSetOptionsPrefix(local_mass_ksp,"ns3d_fieldsplit_1_pc_mass_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (local_mass_ksp); CHKERRQ(ierr);
	

	ierr = KSPSetOperators(local_mass_ksp,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);
	ierr = PCSetOperators(local_mass_pc,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);

	
	ierr = PCApply(local_mass_pc,x,shell->temp_vec);CHKERRQ(ierr);
	PetscInt num_mass_its = 0;
  ierr = KSPGetIterationNumber (local_mass_ksp, &num_mass_its); CHKERRQ(ierr);
	std::cout << "num_mass_its = " << num_mass_its << std::endl;

	//ierr = MatMult(shell->pressure_laplacian_matrix,x,shell->temp_vec_2); CHKERRQ(ierr);
	







	// ******* Apply F_p ************ //
	//
	ierr = MatMult(shell->pressure_convection_diffusion_matrix,shell->temp_vec,shell->temp_vec_2); CHKERRQ(ierr);
	//ierr = MatMult(shell->pressure_convection_diffusion_matrix,x,shell->temp_vec_2); CHKERRQ(ierr);

	
	

	// ******* Apply A_p^-1 ********* //

	PC local_lap_pc;
	ierr = KSPGetPC(shell->inner_lap_ksp,&local_lap_pc); CHKERRQ(ierr);
	ierr = PCSetType(local_lap_pc,PCKSP); CHKERRQ(ierr);
	KSP local_lap_ksp;
	ierr = PCKSPGetKSP(local_lap_pc,&local_lap_ksp); CHKERRQ(ierr);

	
	ierr = KSPSetOptionsPrefix(local_lap_ksp,"ns3d_fieldsplit_1_pc_lap_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (local_lap_ksp); CHKERRQ(ierr);
	

	ierr = KSPSetOperators(local_lap_ksp,shell->pressure_laplacian_matrix,shell->pressure_laplacian_preconditioner); CHKERRQ(ierr);
	ierr = PCSetOperators(local_lap_pc,shell->pressure_laplacian_matrix,shell->pressure_laplacian_preconditioner); CHKERRQ(ierr);

	ierr = PCApply(local_lap_pc,shell->temp_vec_2,y); CHKERRQ(ierr);
	//PCApply(local_lap_pc,x,y);
	PetscInt num_lap_its = 0;
  ierr = KSPGetIterationNumber (local_lap_ksp, &num_lap_its); CHKERRQ(ierr);
	std::cout << "num_lap_its = " << num_lap_its << std::endl;
	PetscReal rnorm = 0.;
	KSPGetResidualNorm(local_lap_ksp,&rnorm);
	std::cout << "norm = " << rnorm << std::endl;
	KSPConvergedReason creason;
	KSPGetConvergedReason(local_lap_ksp,&creason);
	std::cout << "creason = " << creason << std::endl;
	
	ierr = PetscLogStagePop();

  return 0;
}


PetscErrorCode MonolithicShellPCApply(PC pc,Vec x,Vec y)
{

  NSShellPC  *shell;
  PetscErrorCode ierr;

	//ierr = PetscLogStagePush(1);

	//printf ("in Monolithic preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	//VecView(x,PETSC_VIEWER_STDOUT_SELF);
	
	// apply the preconditioner (lu)	
	PC local_velocity_pc;
	ierr = KSPGetPC(shell->inner_velocity_ksp,&local_velocity_pc); CHKERRQ(ierr);
	ierr = PCApply(local_velocity_pc,x,y); CHKERRQ(ierr);

	//PetscInt num_mass_its = 0;
  //ierr = KSPGetIterationNumber (local_mass_ksp, &num_mass_its);
	//printf ("num_mass_its = %d\n",num_mass_its);


	printf ("hi\n");
	//ierr = MatMult(shell->pressure_mass_matrix,x,y);

	//VecView(y,PETSC_VIEWER_STDOUT_SELF);

	//ierr = PetscLogStagePop();
  return 0;
}








PetscErrorCode LSCShellPCApply(PC pc,Vec x,Vec y)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	printf ("inside shell pc lsc james pcapply");

	PetscLogStage pcd_stage;
	//ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	Mat pcmat;
	ierr = PCGetOperators(pc,&pcmat,NULL); CHKERRQ(ierr);

	Mat A,B,C;
	MatSchurComplementGetSubMatrices(pcmat,&A,NULL,&B,&C,NULL);

	KSPSolve(shell->inner_lap_ksp,x,shell->temp_vec_3);

	MatMult(B,shell->temp_vec_3,shell->temp_vec);

	//VecPointwiseMult(shell->temp_vec,shell->temp_vec,shell->lsc_scale);

	MatMult(A,shell->temp_vec,shell->temp_vec_2);

	//VecPointwiseMult(shell->temp_vec_2,shell->temp_vec_2,shell->lsc_scale);

	MatMult(C,shell->temp_vec_2,shell->temp_vec_3);
	KSPSolve(shell->inner_lap_ksp,shell->temp_vec_3,y);
	return(0);



	//ierr = PetscLogStagePop();

  return 0;
}






PetscErrorCode LSCScaledShellPCApply(PC pc,Vec x,Vec y)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	printf ("inside shell pc lsc james pcapply");

	PetscLogStage pcd_stage;
	//ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	Mat pcmat;
	ierr = PCGetOperators(pc,&pcmat,NULL); CHKERRQ(ierr);

	Mat A,B,C;
	MatSchurComplementGetSubMatrices(pcmat,&A,NULL,&B,&C,NULL);

	KSPSolve(shell->inner_lap_ksp,x,shell->temp_vec_3);

	MatMult(B,shell->temp_vec_3,shell->temp_vec);

	VecPointwiseMult(shell->temp_vec,shell->temp_vec,shell->lsc_scale);

	MatMult(A,shell->temp_vec,shell->temp_vec_2);

	VecPointwiseMult(shell->temp_vec_2,shell->temp_vec_2,shell->lsc_scale);

	MatMult(C,shell->temp_vec_2,shell->temp_vec_3);
	KSPSolve(shell->inner_lap_ksp,shell->temp_vec_3,y);
	return(0);



	//ierr = PetscLogStagePop();

  return 0;
}





PetscErrorCode LSCScaledStabilisedShellPCApply(PC pc,Vec x,Vec y)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	//printf ("inside shell pc lsc scaled stabilised james pcapply");

	PetscLogStage pcd_stage;
	//ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	Mat pcmat;
	ierr = PCGetOperators(pc,&pcmat,NULL); CHKERRQ(ierr);

	Mat A,B,C;
	MatSchurComplementGetSubMatrices(pcmat,&A,NULL,&B,&C,NULL);

	KSPSolve(shell->inner_lap_ksp,x,shell->temp_vec_3);

	//std::cout << "x: " << std::endl;
	double norm;
	VecNorm(x,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	//std::cout << "shell->temp_vec_3: " << std::endl;
	VecNorm(shell->temp_vec_3,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	MatMult(B,shell->temp_vec_3,shell->temp_vec);

	//std::cout << "shell->temp_vec: " << std::endl;
	VecNorm(shell->temp_vec,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecPointwiseMult(shell->temp_vec,shell->temp_vec,shell->lsc_scale);

	//std::cout << "shell->temp_vec: " << std::endl;
	VecNorm(shell->temp_vec,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	MatMult(A,shell->temp_vec,shell->temp_vec_2);

	//std::cout << "shell->temp_vec_2: " << std::endl;
	VecNorm(shell->temp_vec_2,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecPointwiseMult(shell->temp_vec_2,shell->temp_vec_2,shell->lsc_scale);

	//std::cout << "shell->temp_vec_2: " << std::endl;
	VecNorm(shell->temp_vec_2,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	MatMult(C,shell->temp_vec_2,shell->temp_vec_3);

	//std::cout << "shell->temp_vec_3: " << std::endl;
	VecNorm(shell->temp_vec_3,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	KSPSolve(shell->inner_lap_ksp,shell->temp_vec_3,y);

	//std::cout << "y: " << std::endl;
	VecNorm(y,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

//lsc_stab_alpha_D_inv
	// + alpha D^-1
	//Vec alpha_D_inv_x;
	VecPointwiseMult(shell->temp_vec_3,shell->lsc_stab_alpha_D_inv,x);

	//std::cout << "shell->temp_vec_3: " << std::endl;
	VecNorm(shell->temp_vec_3,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecAXPY(y,1.0,shell->temp_vec_3);
	
	//std::cout << "y: " << std::endl;
	VecNorm(y,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	return(0);



	//ierr = PetscLogStagePop();

  return 0;
}




PetscErrorCode SIMPLEShellPCApply(PC pc,Vec x,Vec y)
{
  SIMPLEShellPC  *shell;
  PetscErrorCode ierr;

	printf ("inside shell pc lsc james pcapply");

	PetscLogStage pcd_stage;
	//ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	// u_star is u_star
	// delta_p is delta_p
	// F_diag_inv_Bt_delta_p is Bt*delta_p
	// r_u is r_u
	// r_p is r_p
	// B_u_star is B*u_star
	// y_u is y_u
	// y_p is y_p

	// get some matrices we will need
	// setup some temporary vectors
	Vec F_diag_inv,delta_p,F_diag_inv_Bt_delta_p,r_u,r_p,y_u,y_p,r_p_minus_B_u_star,u_star;
	MatGetVecs(shell->A,&F_diag_inv,NULL);
	MatGetVecs(shell->A,&u_star,NULL);
	MatGetVecs(shell->C,&delta_p,NULL);
	MatGetVecs(shell->A,&F_diag_inv_Bt_delta_p,NULL);
	MatGetVecs(shell->C,&r_p_minus_B_u_star,NULL);


	// don't need to allocate
	//MatGetVecs(shell->A,&r_u,NULL);
	//MatGetVecs(shell->C,&r_p,NULL);
	VecGetSubVector(x,shell->velocity_is,&r_u);
	VecGetSubVector(x,shell->pressure_is,&r_p);

	//MatGetVecs(shell->A,&y_u,NULL);
	//MatGetVecs(shell->C,&y_p,NULL);
	VecGetSubVector(y,shell->velocity_is,&y_u);
	VecGetSubVector(y,shell->pressure_is,&y_p);

	//std::cout << "r_u:" << std::endl; 
	double norm;
	VecNorm(r_u,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	//std::cout << "r_p:" << std::endl; 
	VecNorm(r_p,NORM_2,&norm);
	//VecView(r_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	MatGetDiagonal(shell->A,F_diag_inv); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(F_diag_inv);

	//std::cout << "F_diag_inv:" << std::endl; 
	VecNorm(F_diag_inv,NORM_2,&norm);
	//VecView(r_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

		

	

	Mat pcmat;
	ierr = PCGetOperators(pc,&pcmat,NULL); CHKERRQ(ierr);


	// F^-1 * u_star = r_u
	KSPSolve(shell->inner_velocity_ksp,r_u,u_star);

	//std::cout << "u_star:" << std::endl;
	VecNorm(u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// S_approx * delta_p = r_p - B * u_star
	MatMult(shell->B,u_star,r_p_minus_B_u_star);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecAXPY(r_p_minus_B_u_star,-1.0,r_p);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecScale(r_p_minus_B_u_star,-1.0);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	KSPSolve(shell->inner_schur_ksp,r_p_minus_B_u_star,delta_p);

	//std::cout << "delta_p:" << std::endl; 
	VecNorm(delta_p,NORM_2,&norm); 
	//VecView(delta_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// u = u_star - diag(F)^-1 * Bt * delta_p
	VecCopy(u_star,y_u);
	MatMult(shell->Bt,delta_p,F_diag_inv_Bt_delta_p);
	VecPointwiseMult(F_diag_inv_Bt_delta_p,F_diag_inv,F_diag_inv_Bt_delta_p);
	VecAXPY(y_u,-1.0,F_diag_inv_Bt_delta_p);

	//std::cout << "y_u:" << std::endl; 
	VecNorm(y_u,NORM_2,&norm); 
	//VecView(y_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// p = delta_p
	VecCopy(delta_p,y_p);

	//std::cout << "y_p:" << std::endl; 
	VecNorm(y_p,NORM_2,&norm); 
	//VecView(y_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// put the subvector back in its place
	VecRestoreSubVector(y,shell->velocity_is,&y_u);
	VecRestoreSubVector(y,shell->pressure_is,&y_p);


	//std::cout << "y:" << std::endl; 
	VecNorm(y,NORM_2,&norm); 
	//VecView(y,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	
	VecCopy(x,y);

	//std::cout << "y after:" << std::endl; 
	VecNorm(y,NORM_2,&norm); 
	//VecView(y,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	
	

	return(0);



	//ierr = PetscLogStagePop();

  return 0;
}





PetscErrorCode SIMPLERShellPCApply(PC pc,Vec x,Vec y)
{
  SIMPLEShellPC  *shell;
  PetscErrorCode ierr;

	printf ("inside shell pc lsc james pcapply");

	PetscLogStage pcd_stage;
	//ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	// u_star is u_star
	// delta_p is delta_p
	// F_diag_inv_Bt_delta_p is Bt*delta_p
	// r_u is r_u
	// r_p is r_p
	// B_u_star is B*u_star
	// y_u is y_u
	// y_p is y_p

	// get some matrices we will need
	// setup some temporary vectors
	Vec F_diag_inv,delta_p,F_diag_inv_Bt_delta_p,r_u,r_p,y_u,y_p,r_p_minus_B_u_star,u_star,r_u_minus_Bt_p_star,p_star,r_p_minus_B_F_diag_inv_r_u,F_diag_inv_r_u,C_p_star;
	MatGetVecs(shell->A,&F_diag_inv,NULL);
	MatGetVecs(shell->A,&u_star,NULL);
	MatGetVecs(shell->C,&delta_p,NULL);
	MatGetVecs(shell->A,&F_diag_inv_Bt_delta_p,NULL);
	MatGetVecs(shell->C,&r_p_minus_B_u_star,NULL);
	MatGetVecs(shell->A,&r_u_minus_Bt_p_star,NULL);
	MatGetVecs(shell->C,&r_p_minus_B_F_diag_inv_r_u,NULL);
	MatGetVecs(shell->A,&F_diag_inv_r_u,NULL);
	MatGetVecs(shell->C,&p_star,NULL);
	MatGetVecs(shell->C,&C_p_star,NULL);


	// don't need to allocate
	//MatGetVecs(shell->A,&r_u,NULL);
	//MatGetVecs(shell->C,&r_p,NULL);
	VecGetSubVector(x,shell->velocity_is,&r_u);
	VecGetSubVector(x,shell->pressure_is,&r_p);

	//MatGetVecs(shell->A,&y_u,NULL);
	//MatGetVecs(shell->C,&y_p,NULL);
	VecGetSubVector(y,shell->velocity_is,&y_u);
	VecGetSubVector(y,shell->pressure_is,&y_p);

	//std::cout << "r_u:" << std::endl; 
	double norm;
	VecNorm(r_u,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	//std::cout << "r_p:" << std::endl; 
	VecNorm(r_p,NORM_2,&norm);
	//VecView(r_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	MatGetDiagonal(shell->A,F_diag_inv); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(F_diag_inv);

	//std::cout << "F_diag_inv:" << std::endl; 
	VecNorm(F_diag_inv,NORM_2,&norm);
	//VecView(r_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

		

	

	Mat pcmat;
	ierr = PCGetOperators(pc,&pcmat,NULL); CHKERRQ(ierr);


	// S_approx * p_star = r_p - B * F_diag^-1* r_u
	VecCopy(r_u,F_diag_inv_r_u);
	VecPointwiseMult(F_diag_inv_r_u,F_diag_inv,F_diag_inv_r_u);
	MatMult(shell->B,F_diag_inv_r_u,r_p_minus_B_F_diag_inv_r_u);
	VecAXPY(r_p_minus_B_F_diag_inv_r_u,-1.0,r_p);
	VecScale(r_p_minus_B_F_diag_inv_r_u,-1.0);
	KSPSolve(shell->inner_schur_ksp,r_p_minus_B_F_diag_inv_r_u,p_star);
	
	//std::cout << "p_star:" << std::endl;
	VecNorm(p_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	//std::cout << "r_p_minus_B_F_diag_inv_r_u:" << std::endl;
	VecNorm(r_p_minus_B_F_diag_inv_r_u,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	

	// F^-1 * u_star = r_u - Bt * p_star
	MatMult(shell->Bt,p_star,r_u_minus_Bt_p_star);
	VecAXPY(r_u_minus_Bt_p_star,-1.0,r_u);
	VecScale(r_u_minus_Bt_p_star,-1.0);
	
	KSPSolve(shell->inner_velocity_ksp,r_u_minus_Bt_p_star,u_star);

	//std::cout << "u_star:" << std::endl;
	VecNorm(u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// S_approx * delta_p = r_p - B * u_star - C*p_star
	MatMult(shell->B,u_star,r_p_minus_B_u_star);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecAXPY(r_p_minus_B_u_star,-1.0,r_p);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecScale(r_p_minus_B_u_star,-1.0);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	MatMult(shell->C,p_star,C_p_star);
	VecAXPY(r_p_minus_B_u_star,-1.0,C_p_star);
	KSPSolve(shell->inner_schur_ksp,r_p_minus_B_u_star,delta_p);

	//std::cout << "delta_p:" << std::endl; 
	VecNorm(delta_p,NORM_2,&norm); 
	//VecView(delta_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// u = u_star - diag(F)^-1 * Bt * delta_p
	VecCopy(u_star,y_u);
	MatMult(shell->Bt,delta_p,F_diag_inv_Bt_delta_p);
	VecPointwiseMult(F_diag_inv_Bt_delta_p,F_diag_inv,F_diag_inv_Bt_delta_p);
	VecAXPY(y_u,-1.0,F_diag_inv_Bt_delta_p);

	//std::cout << "y_u:" << std::endl; 
	VecNorm(y_u,NORM_2,&norm); 
	//VecView(y_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// p = delta_p
	VecCopy(delta_p,y_p);
	VecAXPY(y_p,1.0,p_star);
	

	//std::cout << "y_p:" << std::endl; 
	VecNorm(y_p,NORM_2,&norm); 
	//VecView(y_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// put the subvector back in its place
	VecRestoreSubVector(y,shell->velocity_is,&y_u);
	VecRestoreSubVector(y,shell->pressure_is,&y_p);


	//std::cout << "y:" << std::endl; 
	VecNorm(y,NORM_2,&norm); 
	//VecView(y,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	
	//VecCopy(x,y);

	//std::cout << "y after:" << std::endl; 
	VecNorm(y,NORM_2,&norm); 
	//VecView(y,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	
	

	return(0);



	//ierr = PetscLogStagePop();

  return 0;
}


/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "SampleNSShellDestroy"
/*
   SampleNSShellDestroy - This routine destroys a user-defined
   preconditioner context.

   Input Parameter:
.  shell - user-defined preconditioner context
*/
PetscErrorCode NSShellDestroy(PC pc)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	//ierr = PetscLogStagePush(1);
  std::cout << "in shell destroy" << std::endl;

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);
	ierr = KSPDestroy(&shell->inner_mass_ksp);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl;
	ierr = KSPDestroy(&shell->inner_lap_ksp);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl;
	ierr = KSPDestroy(&shell->inner_velocity_ksp);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl;
	ierr = KSPDestroy(&shell->inner_schur_ksp);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl;
	ierr = VecDestroy(&shell->temp_vec);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl;
	ierr = VecDestroy(&shell->temp_vec_2);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl;
	ierr = VecDestroy(&shell->temp_vec_3);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl;
	ierr = VecDestroy(&shell->lsc_scale);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl;
	ierr = MatDestroy(&shell->S_approx);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl; //here, it is getting deleted earlier...
//	ierr = MatDestroy(&shell->velocity_matrix);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl;
	ierr = MatDestroy(&shell->lsc_laplacian_matrix);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl;
	ierr = VecDestroy(&shell->lsc_stab_alpha_D_inv);CHKERRQ(ierr);
  std::cout << "in shell destroy" << std::endl;


  ierr = PetscFree(shell);CHKERRQ(ierr);
	
  std::cout << "done shell destroy" << std::endl;
	//ierr = PetscLogStagePop();

  return 0;
}




/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "SampleNSShellDestroy"
/*
   SIMPLEShellDestroy - This routine destroys a user-defined
   preconditioner context.

   Input Parameter:
.  shell - user-defined preconditioner context
*/
PetscErrorCode SIMPLEShellDestroy(PC pc)
{
  SIMPLEShellPC  *shell;
  PetscErrorCode ierr;


	//ierr = PetscLogStagePush(1);
  std::cout << "in SIMPLE shell destroy" << std::endl;

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);
	ierr = KSPDestroy(&shell->inner_schur_ksp);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = KSPDestroy(&shell->inner_velocity_ksp);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = VecDestroy(&shell->temp_vec);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = VecDestroy(&shell->temp_vec_2);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = VecDestroy(&shell->temp_vec_3);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = MatDestroy(&shell->S_approx);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl; //here, it is getting deleted earlier...
//	ierr = MatDestroy(&shell->velocity_matrix);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = MatDestroy(&shell->A);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = MatDestroy(&shell->Bt);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = MatDestroy(&shell->B);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = MatDestroy(&shell->C);CHKERRQ(ierr);



  ierr = PetscFree(shell);CHKERRQ(ierr);
	
  std::cout << "done SIMPLE shell destroy" << std::endl;
	//ierr = PetscLogStagePop();

  return 0;
}



// ************ MATSHELL MULT ***************** //
PetscErrorCode MatShellMultFull(Mat A, Vec vx, Vec vy)
{
	PetscErrorCode ierr;
	PCD2ShellMatrixCtx *ctx;


	ierr = MatShellGetContext(A, (void **)&ctx); CHKERRQ(ierr);

	KSP schur_ksp = ctx->schur_ksp;

	// setup up the operators for the first inner solve
	Mat S, Pmat;
	//MatStructure mat_structure;

	//MatView(pressure_laplacian_matrix->mat(),PETSC_VIEWER_STDOUT_SELF);
	//ierr = KSPGetOperators(schur_ksp,&S,NULL,NULL); CHKERRQ(ierr);
	ierr = KSPGetOperators(schur_ksp,&S,NULL); CHKERRQ(ierr);

	Mat A00;
	Mat Bt;
	Mat B;

	//ierr = MatSchurComplementGetSubmatrices(S,&A00,NULL,&Bt,&B,NULL);
	ierr = MatSchurComplementGetSubMatrices(S,&A00,NULL,&Bt,&B,NULL);

	//MatView(S,PETSC_VIEWER_STDOUT_SELF);
	//MatView(B,PETSC_VIEWER_STDOUT_SELF);

	std::cout << "BT baby" << std::endl;
	MatView(Bt,PETSC_VIEWER_STDOUT_SELF);

	//std::cout << "not made vecs" << std::endl;
	Vec temp_vec_1;
	Vec temp_vec_2;
	Vec inv_diag;
  // create some vectors for temporary pc applying
	PetscInt n_row;
	PetscInt n_row_local;
	//std::cout << "not got mat" << std::endl;
	ierr = MatGetSize(Bt,&n_row,NULL); CHKERRQ(ierr);
	ierr = MatGetLocalSize(Bt,&n_row_local,NULL); CHKERRQ(ierr);
	//std::cout << "got mat sizes rows: " << n_row << ", cols = " << n_col << std::endl;
	// create some vectors that can hold temporary info
	ierr = VecCreate(PETSC_COMM_WORLD,&temp_vec_1); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&temp_vec_2); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&inv_diag); CHKERRQ(ierr);
	ierr = VecSetSizes(temp_vec_1, n_row_local, n_row); CHKERRQ(ierr);
	ierr = VecSetSizes(temp_vec_2, n_row_local, n_row); CHKERRQ(ierr);
	ierr = VecSetSizes(inv_diag, n_row_local, n_row); CHKERRQ(ierr);
	ierr = VecSetFromOptions(temp_vec_1); CHKERRQ(ierr);
	ierr = VecSetFromOptions(temp_vec_2); CHKERRQ(ierr);
	ierr = VecSetFromOptions(inv_diag); CHKERRQ(ierr);

	ierr = MatGetDiagonal(ctx->velocity_mass_matrix,inv_diag);CHKERRQ(ierr);
	//VecView(inv_diag,PETSC_VIEWER_STDOUT_SELF);
	ierr = VecReciprocal(inv_diag);CHKERRQ(ierr);

	PetscInt size;
	ierr = VecGetSize(temp_vec_1, &size); CHKERRQ(ierr);
	//std::cout << "size1 =  " << size << std::endl; 
	ierr = VecGetSize(temp_vec_2, &size); CHKERRQ(ierr);
	//std::cout << "size2 =  " << size << std::endl; 

	//ierr = MatGetSize(ctx->velocity_mass_matrix,&n_row,&n_col); CHKERRQ(ierr);
	//std::cout << "got mat sizes rows: " << n_row << ", cols = " << n_col << std::endl;

	//VecView(inv_diag,PETSC_VIEWER_STDOUT_SELF);
	//MatView(Bt,PETSC_VIEWER_STDOUT_SELF);
	//MatView(ctx->velocity_mass_matrix,PETSC_VIEWER_STDOUT_SELF);
	//MatView(B,PETSC_VIEWER_STDOUT_SELF);

	//std::cout << "made vecs" << std::endl;
	ierr = MatMult(Bt,vx,temp_vec_1); CHKERRQ(ierr);
	//std::cout << "1" << std::endl;
	//ierr = MatMult(A00,temp_vec_1,temp_vec_2); CHKERRQ(ierr);
	ierr = VecPointwiseMult(temp_vec_2,inv_diag,temp_vec_1);
	//std::cout << "2" << std::endl;
	ierr = MatMult(B,temp_vec_2,vy); CHKERRQ(ierr);
	//std::cout << "3" << std::endl;

	ierr = VecDestroy(&temp_vec_1);CHKERRQ(ierr);
	ierr = VecDestroy(&temp_vec_2);CHKERRQ(ierr);
	ierr = VecDestroy(&inv_diag);CHKERRQ(ierr);

	//ierr = MatMult(ctx->pressure_laplacian_matrix,vx,vy); CHKERRQ(ierr);

	return 0;
}








// ********************** SETUP PRECONDITIONERS **************************** //
// setup the preconditioner for use in the solver
// construct submatrices from full matrices that they were assembled as.
int NavierStokesCoupled::setup_preconditioners(TransientLinearImplicitSystem * system, Mat& B)
{


	PetscErrorCode ierr;

	if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 6
			|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 7)
	{
		// for the monolithic preconditioner we currently don't need to
		// do anything, the preconditioner should be built.
		// we make sure not to put the 


	  PetscErrorCode ierr;

		//get petsc solver
		PetscLinearSolver<Number>* system_linear_solver =
				libmesh_cast_ptr<PetscLinearSolver<Number>* >
				(system->linear_solver.get());
		KSP system_ksp = system_linear_solver->ksp();
		PC pc;
		ierr = KSPGetPC(system_ksp,&pc); CHKERRQ(ierr);



		std::cout << "Setting up Preconditioner." << std::endl;
		int num_splits = 0;
		ierr = KSPSetUp(system_ksp); CHKERRQ(ierr);

		//KSP* subksp;	//subksps
		KSP* subksp;	//subksps
		ierr = PCFieldSplitGetSubKSP(pc,&num_splits,&subksp); CHKERRQ(ierr);

		// set to not reuse the preconditioner
		//ierr = KSPSetReusePreconditioner(subksp[0],PETSC_TRUE); CHKERRQ(ierr);
		// set to not reuse the preconditioner
		//ierr = KSPSetReusePreconditioner(subksp[1],PETSC_TRUE); CHKERRQ(ierr);


		// 1.) We need to take the 1D matrix from the system matrix and put it in the preconditioner matrix
	}
	else if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 8
			|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 9
			|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 10
			|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 11
			|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 12)
	{
		// set up some stuff for the velocity preconditioner
		// need to extract the navier stokes block and add ones on the diagonal of the pressure bit


		std::cout << "Constructing submatrices." << std::endl;
		const unsigned int u_var = system->variable_number ("u");
		const unsigned int v_var = system->variable_number ("v");
		unsigned int w_var = 0;
		if(threed)
			w_var = system->variable_number ("w");
		const unsigned int p_var = system->variable_number ("p");

		const unsigned int Q_var = system->variable_number ("Q");
		const unsigned int P_var = system->variable_number ("P");



		std::vector<dof_id_type> NS_var_idx;
		std::vector<dof_id_type> u_var_idx;
		std::vector<dof_id_type> v_var_idx;
		std::vector<dof_id_type> w_var_idx;
		std::vector<dof_id_type> p_var_idx;
		std::vector<dof_id_type> Q_var_idx;
		std::vector<dof_id_type> P_var_idx;
		system->get_dof_map().local_variable_indices
			(u_var_idx, system->get_mesh(), u_var);
		system->get_dof_map().local_variable_indices
			(v_var_idx, system->get_mesh(), v_var);

		if(threed)
			system->get_dof_map().local_variable_indices
				(w_var_idx, system->get_mesh(), w_var);

		system->get_dof_map().local_variable_indices
			(p_var_idx, system->get_mesh(), p_var);

		system->get_dof_map().local_variable_indices
			(Q_var_idx, system->get_mesh(), Q_var);
		system->get_dof_map().local_variable_indices
			(P_var_idx, system->get_mesh(), P_var);

		// concatenate all velocity dofs
		NS_var_idx.insert(NS_var_idx.end(),u_var_idx.begin(),u_var_idx.end());
		NS_var_idx.insert(NS_var_idx.end(),v_var_idx.begin(),v_var_idx.end());
		NS_var_idx.insert(NS_var_idx.end(),w_var_idx.begin(),w_var_idx.end());
		NS_var_idx.insert(NS_var_idx.end(),p_var_idx.begin(),p_var_idx.end());

		std::sort(NS_var_idx.begin(),NS_var_idx.end());

		if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 8 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 10
			|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 11)
		{
			velocity_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
			//system->request_matrix("Velocity Matrix")->close();
			//system->request_matrix("Velocity Matrix")->add(-.5,*system->request_matrix("Velocity Matrix"));
			std::cout << "hey babe?" << std::endl;
			//std::exit(0);
			system->request_matrix("Velocity Matrix")->create_submatrix(*velocity_matrix,NS_var_idx,NS_var_idx);
			velocity_matrix->close();
		}
		else if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") != 12)
		{

			velocity_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
			system->request_matrix("System Matrix")->create_submatrix(*velocity_matrix,NS_var_idx,NS_var_idx);
			velocity_matrix->close();
		}
		// don't need a matrix for preconditioner_type 12
	
		// figure out which are the pressure dofs on the boundary. then figure out how they fit into the 
		// submatrix...
		
	



		// set up a global Vec saying which columns have entries in the B1 matrix
		
		if(es->parameters.get<bool>("multiple_column_solve"))
		{

			PetscInt n_col = Q_var_idx.size() + P_var_idx.size();

			std::vector<dof_id_type> pressure_coupling_dofs;
			std::vector<dof_id_type> flux_coupling_dofs;
			const MeshBase& mesh = es->get_mesh();
			const unsigned int p_var = system->variable_number ("P");
			const unsigned int q_var = system->variable_number ("Q");
			const DofMap & dof_map = system->get_dof_map();
			std::vector<dof_id_type> dof_indices;
			std::vector<dof_id_type> dof_indices_p;
			std::vector<dof_id_type> dof_indices_q;

			MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
			const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
						

			for ( ; el != end_el; ++el)
			{
				const Elem* elem = *el;
				//only concerned with 1d elements
				if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
				{
					for (unsigned int s=0; s<elem->n_sides(); s++)
					{

						//for some reason it is natural to have more than one boundary id per side or even node
						std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

						if(boundary_ids.size() > 0) 
						{ 

							
							// boundary conditions that are not inflow or walls or ends of airways
							if(boundary_ids[0] > 0 && boundary_ids[0] < 1000)
							{	     
								std::cout << "boundary_ids[0] = " << boundary_ids[0] << std::endl;
								// Get the degree of freedom indices for the
								// current element.
								dof_map.dof_indices (elem, dof_indices_p, p_var);
								pressure_coupling_dofs.push_back(dof_indices_p[0]);

								dof_map.dof_indices (elem, dof_indices_q, q_var);
								flux_coupling_dofs.push_back(dof_indices_q[0]);

							}
						}
					}
				}
			}

			std::cout << "pressure coupling dofs: " << std::endl;
			for(unsigned int i=0; i<pressure_coupling_dofs.size(); i++)
			{
				std::cout << "pressure_coupling_dofs[" << i <<"] = " << pressure_coupling_dofs[i] << std::endl;
			}
	
			std::vector<dof_id_type> PQ_var_idx;
			PQ_var_idx.insert(PQ_var_idx.end(), P_var_idx.begin(), P_var_idx.end() );
			PQ_var_idx.insert(PQ_var_idx.end(), Q_var_idx.begin(), Q_var_idx.end() );
			std::sort(PQ_var_idx.begin(),PQ_var_idx.end());

			ierr = VecCreate(PETSC_COMM_WORLD,&non_zero_cols); CHKERRQ(ierr);
			ierr = VecSetSizes(non_zero_cols, pressure_coupling_dofs.size(), pressure_coupling_dofs.size()); CHKERRQ(ierr);
			ierr = VecSetFromOptions(non_zero_cols); CHKERRQ(ierr);	

			ierr = VecCreate(PETSC_COMM_WORLD,&non_zero_rows); CHKERRQ(ierr);
			ierr = VecSetSizes(non_zero_rows, flux_coupling_dofs.size(), flux_coupling_dofs.size()); CHKERRQ(ierr);
			ierr = VecSetFromOptions(non_zero_rows); CHKERRQ(ierr);	

		  std::vector<dof_id_type>::iterator p;
			unsigned int pressure_counter = 0;
			unsigned int flux_counter = 0;
			for(unsigned int i=0; i<PQ_var_idx.size(); i++)
			{

				p = std::find (pressure_coupling_dofs.begin(), pressure_coupling_dofs.end(), PQ_var_idx[i]);
				if (p != pressure_coupling_dofs.end())
				{
		    	VecSetValue(non_zero_cols,pressure_counter,i,INSERT_VALUES);
					pressure_counter++;
				}


				p = std::find (flux_coupling_dofs.begin(), flux_coupling_dofs.end(), PQ_var_idx[i]);
				if (p != flux_coupling_dofs.end())
				{
		    	VecSetValue(non_zero_rows,flux_counter,i,INSERT_VALUES);
					flux_counter++;
				}
			}

			ierr = VecAssemblyBegin(non_zero_cols); CHKERRQ(ierr);
			ierr = VecAssemblyEnd(non_zero_cols); CHKERRQ(ierr);

			ierr = VecAssemblyBegin(non_zero_rows); CHKERRQ(ierr);
			ierr = VecAssemblyEnd(non_zero_rows); CHKERRQ(ierr);



			std::cout << "PQ_var_idx.size() = " << PQ_var_idx.size() << std::endl;
			//std::cout << "the non_zero_cols of Bt" << std::endl;
 			//ierr = VecView(non_zero_cols,PETSC_VIEWER_STDOUT_WORLD);
			//std::cout << "the non_zero_rows of Bt" << std::endl;
 			//ierr = VecView(non_zero_rows,PETSC_VIEWER_STDOUT_WORLD);
			

		}
		



		//PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF,PETSC_VIEWER_ASCII_DENSE);


		//MatView(velocity_matrix->mat(),PETSC_VIEWER_STDOUT_SELF);
		/*
		std::cout << "velocity matrix in full:\n\n\n" << std::endl;
		system->request_matrix("Velocity Matrix")->print();

		std::cout << "velocity matrux in petsc\n\n\n" << std::endl;
		//PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF,PETSC_VIEWER_ASCII_MATLAB);
		PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF,PETSC_VIEWER_ASCII_DENSE);


		MatView(velocity_matrix->mat(),PETSC_VIEWER_STDOUT_SELF);


		std::exit(0);
		*/

		// now set up the shell preconditioner

		//get petsc solver
		PetscLinearSolver<Number>* system_linear_solver =
				libmesh_cast_ptr<PetscLinearSolver<Number>* >
				(system->linear_solver.get());
		KSP system_ksp = system_linear_solver->ksp();
		PC pc;
		ierr = KSPGetPC(system_ksp,&pc); CHKERRQ(ierr);



		int num_splits = 0;
		if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") != 12)
		{

			std::cout << "Setting up Velocity Monolithic preconditioner." << std::endl;
			ierr = KSPSetUp(system_ksp); CHKERRQ(ierr);
			KSP* subksp;	//subksps
			ierr = PCFieldSplitGetSubKSP(pc,&num_splits,&subksp); CHKERRQ(ierr);

			// set ksp to preonly because we are defining the inverse
			ierr = KSPSetType(subksp[1], KSPPREONLY); CHKERRQ(ierr);

			PC schur_pc;
			ierr = KSPGetPC(subksp[1],&schur_pc); CHKERRQ(ierr);

			// set to not reuse the preconditioner
			//ierr = KSPSetReusePreconditioner(subksp[0],PETSC_FALSE); CHKERRQ(ierr);
			// set to not reuse the preconditioner
			//ierr = KSPSetReusePreconditioner(subksp[1],PETSC_TRUE); CHKERRQ(ierr);


			// (Required) Indicate to PETSc that we're using a "shell" preconditioner 
			ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);

			// destroy old shell before it is created
			if(mono_shell_pc_created)
			{
				std::cout << "Destroying old contents of shell pc" << std::endl;
				ierr = NSShellDestroy(schur_pc);
			}

			// (Optional) Create a context for the user-defined preconditioner; this
			//	context can be used to contain any application-specific data. 
			ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

			// (Required) Set the user-defined routine for applying the preconditioner 
			ierr = PCShellSetApply(schur_pc,MonolithicShellPCApply); CHKERRQ(ierr);
			ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

			// (Optional) Set user-defined function to free objects used by custom preconditioner 
			ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

			// (Optional) Set a name for the preconditioner, used for PCView() 
			ierr = PCShellSetName(schur_pc,"Velocity Preconditioner"); CHKERRQ(ierr);

			// (Optional) Do any setup required for the preconditioner 
			// Note: This function could be set with PCShellSetSetUp and it would be called when necessary 

			if(es->parameters.get<bool>("multiple_column_solve"))
			{
				ierr = Monolithic2ShellPCSetUp(schur_pc,velocity_matrix->mat(),non_zero_cols,non_zero_rows,subksp[1]); CHKERRQ(ierr);
			}
			else
			{
				ierr = MonolithicShellPCSetUp(schur_pc,velocity_matrix->mat(),subksp[1]); CHKERRQ(ierr);
			}

			mono_shell_pc_created = true;
		}
		else
		{
			// don't need to do anything, should be default
		}

	}


	if(es->parameters.get<unsigned int>("preconditioner_type_3d") != 0)
	{


	  PetscErrorCode ierr;

		//get petsc solver
		PetscLinearSolver<Number>* system_linear_solver =
				libmesh_cast_ptr<PetscLinearSolver<Number>* >
				(system->linear_solver.get());
		KSP system_ksp = system_linear_solver->ksp();
		PC pc;
		ierr = KSPGetPC(system_ksp,&pc); CHKERRQ(ierr);

		// ********* CONSTRUCT SUBMATRICES ********** //

		std::cout << "Constructing 3D submatrices." << std::endl;
		const unsigned int u_var = system->variable_number ("u");
		const unsigned int v_var = system->variable_number ("v");
		unsigned int w_var = 0;
		if(threed)
			w_var = system->variable_number ("w");
		const unsigned int p_var = system->variable_number ("p");

		std::vector<dof_id_type> var_idx;
		system->get_dof_map().local_variable_indices
			(var_idx, system->get_mesh(), p_var);


		std::vector<dof_id_type> U_var_idx;
		std::vector<dof_id_type> u_var_idx;
		std::vector<dof_id_type> v_var_idx;
		std::vector<dof_id_type> w_var_idx;
		system->get_dof_map().local_variable_indices
			(u_var_idx, system->get_mesh(), u_var);
		system->get_dof_map().local_variable_indices
			(v_var_idx, system->get_mesh(), v_var);

		if(threed)
			system->get_dof_map().local_variable_indices
				(w_var_idx, system->get_mesh(), w_var);

		// concatenate all velocity dofs
		U_var_idx.insert(U_var_idx.end(),u_var_idx.begin(),u_var_idx.end());
		U_var_idx.insert(U_var_idx.end(),v_var_idx.begin(),v_var_idx.end());
		U_var_idx.insert(U_var_idx.end(),w_var_idx.begin(),w_var_idx.end());

		std::sort(U_var_idx.begin(),U_var_idx.end());

	
	

	


		std::vector<dof_id_type> rows;
		rows.push_back(0);


		// pressure mass matrix
		pressure_mass_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
		system->request_matrix("Pressure Mass Matrix")->create_submatrix(*pressure_mass_matrix,var_idx,var_idx);
		//if(!es->parameters.get<unsigned int>("problem_type") == 4 && !es->parameters.get<unsigned int>("problem_type") == 5)
		//pressure_mass_matrix->zero_rows(rows,1.0);
		pressure_mass_matrix->close();

		// scaled pressure mass matrix
		scaled_pressure_mass_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
		system->request_matrix("Pressure Mass Matrix")->create_submatrix(*scaled_pressure_mass_matrix,var_idx,var_idx);
		scaled_pressure_mass_matrix->zero();
		scaled_pressure_mass_matrix->add(es->parameters.get<Real>("reynolds_number"),*pressure_mass_matrix);
		//if(!es->parameters.get<unsigned int>("problem_type") == 4 && !es->parameters.get<unsigned int>("problem_type") == 5)
		//scaled_pressure_masNus_matrix->zero_rows(rows,1.0);
		scaled_pressure_mass_matrix->close();
	

		// pressure laplacian matrix
		pressure_laplacian_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
		// need to pin laplacian matrix... at the first pressure node
		system->request_matrix("Pressure Laplacian Matrix")->create_submatrix(*pressure_laplacian_matrix,var_idx,var_idx);
		if(es->parameters.get<bool>("pin_pressure"))
		{
			pressure_laplacian_matrix->zero_rows(rows,1.0);
		}

		pressure_laplacian_matrix->close();
		//system->request_matrix("Pressure Laplacian Matrix")->print();

	

		// pressure convection diffusion matrix
		pressure_convection_diffusion_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
		system->request_matrix("Pressure Convection Diffusion Matrix")->create_submatrix(*pressure_convection_diffusion_matrix,var_idx,var_idx);
		if(es->parameters.get<bool>("pin_pressure"))
		{
			pressure_convection_diffusion_matrix->zero_rows(rows,1.0/es->parameters.get<Real>("reynolds_number"));
		}

		pressure_convection_diffusion_matrix->close();

	

		// pressure convection diffusion matrix
		velocity_mass_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
		system->request_matrix("Velocity Mass Matrix")->create_submatrix(*velocity_mass_matrix,U_var_idx,U_var_idx);
		velocity_mass_matrix->close();


		//system->request_matrix("Pressure Laplacian Matrix")->print();
		//system->request_matrix("Pressure Convection Diffusion Matrix")->print();

		//pressure_laplacian_matrix->add(-1.,*pressure_convection_diffusion_matrix);
		//MatView(pressure_laplacian_matrix->mat(),PETSC_VIEWER_STDOUT_SELF);
	

		// set pressure in system matrix and preconditioner

		std::vector<dof_id_type> rows_full;
		rows_full.push_back(var_idx[0]);
		//system->request_matrix("System Matrix")->zero_rows(rows_full,1.0);
		//system->request_matrix("Preconditioner")->zero_rows(rows_full,1.0);
		//system->request_matrix("System Matrix")->close();
		//system->request_matrix("Preconditioner")->close();

		PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF,PETSC_VIEWER_ASCII_MATLAB);
		//MatView(scaled_pressure_mass_matrix->mat(),PETSC_VIEWER_STDOUT_SELF);
		//MatView(pressure_laplacian_matrix->mat(),PETSC_VIEWER_STDOUT_SELF);
		//MatView(velocity_mass_matrix->mat(),PETSC_VIEWER_STDOUT_SELF);
		//MatView(pressure_convection_diffusion_matrix->mat(),PETSC_VIEWER_STDOUT_SELF);

		//system->request_matrix("Pressure Mass Matrix")->print();

		std::cout << "Done setting up the sub matrices." << std::endl;
		// *********************************************************************** //












		// if not SIMPLE type preconditioners
		if(es->parameters.get<unsigned int>("preconditioner_type_3d") != 10 
			&& es->parameters.get<unsigned int>("preconditioner_type_3d") != 11
			&& es->parameters.get<unsigned int>("preconditioner_type_3d") != 12)
		{
	
			int num_splits = 0;




			// ************* CONSTRUCT THE SCHUR COMPLEMENT PRECONDITIONERS ********** //
			ierr = KSPSetUp(system_ksp); CHKERRQ(ierr);
			KSP* subksp;	//subksps
			ierr = PCFieldSplitGetSubKSP(pc,&num_splits,&subksp); CHKERRQ(ierr);
			KSP* velocity_subksp;



			std::cout << "INSIDE SETUP PC = " << std::endl;
			std::cout << "ya = " << std::endl;

			if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 0)
			{
				velocity_subksp = subksp;
			}
			else
			{


				//ierr = KSPSetUp(subksp[0]);
				PC sub_pc;
				ierr = KSPGetPC(subksp[0],&sub_pc);// CHKERRQ(ierr);
				KSP* subsubksp;	//subksps
				ierr = PCFieldSplitGetSubKSP(sub_pc,NULL,&subsubksp);// CHKERRQ(ierr);
				velocity_subksp = subsubksp;
			}


			std::cout << "Setting up the Schur preconditioner Petsc settings" << std::endl;

			// now set the preconditioner stuff for the schur preconditioners
			if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 1)
			{

				//PCFieldSplitSchurPrecondition(pc,PC_FIELDSPLIT_SCHUR_PRE_USER,scaled_pressure_mass_matrix->mat());

		
				// note even though the preconditioner says it is using A11 it is not
				Mat Amat, Pmat;
				//MatStructure mat_structure;

				//MatView(scaled_pressure_mass_matrix->mat(),PETSC_VIEWER_STDOUT_SELF);
				//ierr = KSPGetOperators(velocity_subksp[1],&Amat,&Pmat,&mat_structure); CHKERRQ(ierr);
				ierr = KSPGetOperators(velocity_subksp[1],&Amat,&Pmat); CHKERRQ(ierr);
				//ierr = PetscObjectReference((PetscObject)Amat); CHKERRQ(ierr);
				//ierr = PetscObjectReference((PetscObject)mat_structure); CHKERRQ(ierr);
				//ierr = KSPSetOperators(velocity_subksp[1],scaled_pressure_mass_matrix->mat(),scaled_pressure_mass_matrix->mat(),mat_structure); CHKERRQ(ierr);
				ierr = KSPSetOperators(velocity_subksp[1],Amat,scaled_pressure_mass_matrix->mat()); CHKERRQ(ierr);
		
		
		
			}
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 2)
			{

				std::cout << "Setting up LSC preconditioner." << std::endl;

				// set ksp to preonly because we are defining the inverse
				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);



				// (Required) Indicate to PETSc that we're using a "lsc" preconditioner 
				ierr = PCSetType(schur_pc,PCLSC); CHKERRQ(ierr);

		
		
			}
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 3)	// shell pc..
			{

				std::cout << "Setting up Pressure preconditioner." << std::endl;

				// set ksp to preonly because we are defining the inverse
				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);



				// (Required) Indicate to PETSc that we're using a "shell" preconditioner 
				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);

				// destroy old shell before it is created
				if(shell_pc_created)
				{
					std::cout << "Destroying old contents of shell pc" << std::endl;
					ierr = NSShellDestroy(schur_pc);
				}

				// (Optional) Create a context for the user-defined preconditioner; this
				//	context can be used to contain any application-specific data. 
				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

				// (Required) Set the user-defined routine for applying the preconditioner 
				ierr = PCShellSetApply(schur_pc,PressureShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				// (Optional) Set user-defined function to free objects used by custom preconditioner 
				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

				// (Optional) Set a name for the preconditioner, used for PCView() 
				ierr = PCShellSetName(schur_pc,"Pressure Preconditioner"); CHKERRQ(ierr);

				// (Optional) Do any setup required for the preconditioner 
				// Note: This function could be set with PCShellSetSetUp and it would be called when necessary 
				ierr = PressureShellPCSetUp(schur_pc,scaled_pressure_mass_matrix->mat(),velocity_subksp[1]); CHKERRQ(ierr);

				shell_pc_created = true;
			}
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 4)	// shell pc..
			{

				std::cout << "Setting up PCD preconditioner." << std::endl;

				// set ksp to preonly because we are defining the inverse
				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);


				/* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);

				if(shell_pc_created)
				{
					std::cout << "Destroying old contents of shell pc" << std::endl;
					ierr = NSShellDestroy(schur_pc);
				}

				/* (Optional) Create a context for the user-defined preconditioner; this
					context can be used to contain any application-specific data. */
				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

				/* (Required) Set the user-defined routine for applying the preconditioner */
				ierr = PCShellSetApply(schur_pc,PCDShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				/* (Optional) Set user-defined function to free objects used by custom preconditioner */
				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

				/* (Optional) Set a name for the preconditioner, used for PCView() */
				ierr = PCShellSetName(schur_pc,"PCD Preconditioner"); CHKERRQ(ierr);

				/* (Optional) Do any setup required for the preconditioner */
				/* Note: This function could be set with PCShellSetSetUp and it would be called when necessary */
				ierr = PCDShellPCSetUp(schur_pc,pressure_mass_matrix->mat(),pressure_laplacian_matrix->mat(),pressure_laplacian_matrix->mat(),pressure_convection_diffusion_matrix->mat(),velocity_subksp[1]); CHKERRQ(ierr);

				shell_pc_created = true;
			}

			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 5)	// shell pc..+ shel matrix pcd2.0
			{

		


				std::cout << "hello" << std::endl;


				mat_ctx.m = var_idx.size();
				mat_ctx.n = var_idx.size();
				mat_ctx.pressure_laplacian_matrix = pressure_laplacian_matrix->mat();
				mat_ctx.velocity_mass_matrix = velocity_mass_matrix->mat();
				mat_ctx.schur_ksp = velocity_subksp[1];

				PetscInt mn = var_idx.size() * var_idx.size();

				std::cout << "m = " << mat_ctx.m  << std::endl;
				std::cout << "n = " << mat_ctx.n  << std::endl;
				std::cout << "mn = " << mn  << std::endl;
				ierr = MatCreateShell(PETSC_COMM_WORLD,mat_ctx.m,mat_ctx.n,PETSC_DETERMINE,PETSC_DETERMINE,&mat_ctx,&B);  CHKERRQ(ierr);
				ierr = MatShellSetOperation(B, MATOP_MULT, (void(*)(void))MatShellMultFull); CHKERRQ(ierr);
		
				//get operators from the outer pcshell solver
				Mat Amat, Pmat;
				//MatStructure mat_structure;
				//ierr = KSPGetOperators(velocity_subksp[1],&Amat,&Pmat,&mat_structure); CHKERRQ(ierr);
				ierr = KSPGetOperators(velocity_subksp[1],&Amat,&Pmat); CHKERRQ(ierr);

				// set ksp to preonly because we are defining the inverse
				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);


				/* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);


				if(shell_pc_created)
				{
					std::cout << "destroying old contents of shell pc" << std::endl;
					ierr = NSShellDestroy(schur_pc);
				}

				/* (Optional) Create a context for the user-defined preconditioner; this
					context can be used to contain any application-specific data. */
				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

				/* (Required) Set the user-defined routine for applying the preconditioner */
				ierr = PCShellSetApply(schur_pc,PCD2ShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				/* (Optional) Set user-defined function to free objects used by custom preconditioner */
				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

				/* (Optional) Set a name for the preconditioner, used for PCView() */
				ierr = PCShellSetName(schur_pc,"PCD2 Preconditioner"); CHKERRQ(ierr);

				/* (Optional) Do any setup required for the preconditioner */
				/* Note: This function could be set with PCShellSetSetUp and it would be called when necessary */
				//ierr = PCDShellPCSetUp(schur_pc,pressure_mass_matrix->mat(),pressure_laplacian_matrix->mat(),pressure_laplacian_matrix->mat(),pressure_convection_diffusion_matrix->mat(),velocity_subksp[1]); CHKERRQ(ierr);
				std::cout << "ya" << std::endl;
				ierr = PCD2ShellPCSetUp(schur_pc,pressure_mass_matrix->mat(),B,pressure_laplacian_matrix->mat(),pressure_convection_diffusion_matrix->mat(),velocity_subksp[1]); CHKERRQ(ierr);

				shell_pc_created = true;
			}
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 6)	// shell pc..
			{

				std::cout << "Setting up LSC James preconditioner." << std::endl;

				// set ksp to preonly because we are defining the inverse
				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);


				/* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);

				if(shell_pc_created)
				{
					std::cout << "Destroying old contents of shell pc" << std::endl;
					ierr = NSShellDestroy(schur_pc);
				}

				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

				ierr = PCShellSetApply(schur_pc,LSCShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

				ierr = PCShellSetName(schur_pc,"LSC James Preconditioner"); CHKERRQ(ierr);

				ierr = LSCShellPCSetUp(schur_pc,velocity_subksp[1]); CHKERRQ(ierr);

				shell_pc_created = true;
			}
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 7)	// shell pc..
			{

				std::cout << "Setting up LSC James scaled preconditioner." << std::endl;

				// set ksp to preonly because we are defining the inverse
				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);


				/* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);

				if(shell_pc_created)
				{
					std::cout << "Destroying old contents of shell pc" << std::endl;
					ierr = NSShellDestroy(schur_pc);
				}

				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

				ierr = PCShellSetApply(schur_pc,LSCScaledShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

				ierr = PCShellSetName(schur_pc,"LSC James Scaled Preconditioner"); CHKERRQ(ierr);

				ierr = LSCScaledShellPCSetUp(schur_pc,velocity_mass_matrix->mat(),velocity_subksp[1]); CHKERRQ(ierr);

				shell_pc_created = true;
			}
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 8)
			{
				// do nothing, we should have set this all up on the command line
			}
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 9)
			{
			
				std::cout << "Setting up LSC James scaled stabilised preconditioner." << std::endl;

				// set ksp to preonly because we are defining the inverse
				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				std::cout << "1" << std::endl;
				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);

				std::cout << "2" << std::endl;

				/* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);

				if(shell_pc_created)
				{
					std::cout << "Destroying old contents of shell pc" << std::endl;
					ierr = NSShellDestroy(schur_pc);
				}

				std::cout << "3" << std::endl;
				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

				std::cout << "4" << std::endl;
				ierr = PCShellSetApply(schur_pc,LSCScaledStabilisedShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				std::cout << "5" << std::endl;
				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

				ierr = PCShellSetName(schur_pc,"LSC James Scaled Stabilised Preconditioner"); CHKERRQ(ierr);

				std::cout << "6" << std::endl;
				       
				ierr = LSCScaledStabilisedShellPCSetUp(schur_pc,velocity_mass_matrix->mat(),velocity_subksp[1]); CHKERRQ(ierr);

				std::cout << "7" << std::endl;
				shell_pc_created = true;
			}
		}
		else
		{

			int num_splits = 0;


			// ************* CONSTRUCT THE SIMPLE-type PRECONDITIONERS ********** //
			ierr = KSPSetUp(system_ksp); CHKERRQ(ierr);

			std::cout << "INSIDE SETUP SIMPLE-type PC = " << std::endl;
			std::cout << "ya = " << std::endl;

			KSP velocity_ksp;
			PC outer_pc;
			ierr = KSPGetPC(system_ksp,&outer_pc);// CHKERRQ(ierr);

			if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 0)
			{
				velocity_ksp = system_ksp;
			}
			else
			{
				KSP* subksp;	//subksps
				ierr = PCFieldSplitGetSubKSP(outer_pc,NULL,&subksp);// CHKERRQ(ierr);
				velocity_ksp = subksp[0];
			}

			// set up the submatrices for use
			setup_is_simple(system,outer_pc);


			std::cout << "Setting up the SIMPLE-type preconditioner Petsc settings" << std::endl;
			// SIMPLE preconditioner
			if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 10)
			{

				std::cout << "Setting up SIMPLE preconditioner." << std::endl;

				// ksp should have type set to FGMRES but let us make sure
				ierr = KSPSetType(velocity_ksp, KSPFGMRES); CHKERRQ(ierr);

				std::cout << "1" << std::endl;
				PC velocity_pc;
				ierr = KSPGetPC(velocity_ksp,&velocity_pc); CHKERRQ(ierr);

				std::cout << "2" << std::endl;

				/* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
				ierr = PCSetType(velocity_pc,PCSHELL); CHKERRQ(ierr);

				if(shell_pc_created)
				{
					std::cout << "Destroying old contents of shell pc" << std::endl;
					ierr = SIMPLEShellDestroy(velocity_pc);
				}

				std::cout << "3" << std::endl;
				ierr = ShellPCCreate(&simple_shell); CHKERRQ(ierr);

				std::cout << "4" << std::endl;
				ierr = PCShellSetApply(velocity_pc,SIMPLEShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(velocity_pc,simple_shell); CHKERRQ(ierr);

				std::cout << "5" << std::endl;
				ierr = PCShellSetDestroy(velocity_pc,SIMPLEShellDestroy); CHKERRQ(ierr);

				ierr = PCShellSetName(velocity_pc,"SIMPLE Preconditioner"); CHKERRQ(ierr);

				std::cout << "6" << std::endl;
				       
				ierr = SIMPLEShellPCSetUp(velocity_pc,velocity_is,pressure_is,velocity_ksp); CHKERRQ(ierr);

				std::cout << "7" << std::endl;
				shell_pc_created = true;
			}
			// SIMPLER preconditioner
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 11)
			{

				std::cout << "Setting up SIMPLER preconditioner." << std::endl;

				// ksp should have type set to FGMRES but let us make sure
				ierr = KSPSetType(velocity_ksp, KSPFGMRES); CHKERRQ(ierr);

				std::cout << "1" << std::endl;
				PC velocity_pc;
				ierr = KSPGetPC(velocity_ksp,&velocity_pc); CHKERRQ(ierr);

				std::cout << "2" << std::endl;

				/* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
				ierr = PCSetType(velocity_pc,PCSHELL); CHKERRQ(ierr);

				if(shell_pc_created)
				{
					std::cout << "Destroying old contents of shell pc" << std::endl;
					ierr = SIMPLEShellDestroy(velocity_pc);
				}

				std::cout << "3" << std::endl;
				ierr = ShellPCCreate(&simple_shell); CHKERRQ(ierr);

				std::cout << "4" << std::endl;
				ierr = PCShellSetApply(velocity_pc,SIMPLERShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(velocity_pc,simple_shell); CHKERRQ(ierr);

				std::cout << "5" << std::endl;
				ierr = PCShellSetDestroy(velocity_pc,SIMPLEShellDestroy); CHKERRQ(ierr);

				ierr = PCShellSetName(velocity_pc,"SIMPLER Preconditioner"); CHKERRQ(ierr);

				std::cout << "6" << std::endl;
				       
				ierr = SIMPLEShellPCSetUp(velocity_pc,velocity_is,pressure_is,velocity_ksp); CHKERRQ(ierr);

				std::cout << "7" << std::endl;
				shell_pc_created = true;
			}

			// SIMPLERC preconditioner (SIMPLER with rowsum approx)
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 12)
			{

				std::cout << "Setting up SIMPLERC preconditioner." << std::endl;

				// ksp should have type set to FGMRES but let us make sure
				ierr = KSPSetType(velocity_ksp, KSPFGMRES); CHKERRQ(ierr);

				std::cout << "1" << std::endl;
				PC velocity_pc;
				ierr = KSPGetPC(velocity_ksp,&velocity_pc); CHKERRQ(ierr);

				std::cout << "2" << std::endl;

				/* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
				ierr = PCSetType(velocity_pc,PCSHELL); CHKERRQ(ierr);

				if(shell_pc_created)
				{
					std::cout << "Destroying old contents of shell pc" << std::endl;
					ierr = SIMPLEShellDestroy(velocity_pc);
				}

				std::cout << "3" << std::endl;
				ierr = ShellPCCreate(&simple_shell); CHKERRQ(ierr);

				std::cout << "4" << std::endl;
				ierr = PCShellSetApply(velocity_pc,SIMPLERShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(velocity_pc,simple_shell); CHKERRQ(ierr);

				std::cout << "5" << std::endl;
				ierr = PCShellSetDestroy(velocity_pc,SIMPLEShellDestroy); CHKERRQ(ierr);

				ierr = PCShellSetName(velocity_pc,"SIMPLERC Preconditioner"); CHKERRQ(ierr);

				std::cout << "6" << std::endl;
				       
				ierr = SIMPLECShellPCSetUp(velocity_pc,velocity_is,pressure_is,velocity_ksp); CHKERRQ(ierr);

				std::cout << "7" << std::endl;
				shell_pc_created = true;
			}
		}
	}
}




// field_number is the IS that we want to return
void NavierStokesCoupled::setup_is_simple (TransientLinearImplicitSystem * sys,
                            PC my_pc)
{


  
  for(unsigned int i=0; i<2; i++)
  {
	std::string field_name;
	if(i==0)
	{
		field_name = "0";
	}
	else
	{
		field_name = "1";

	}
	

	  // ****************** FIGURE OUT WHAT INDICES TO RETURN ****** //

	  std::cout << "Setting up the IS for fieldsplit: " << i << std::endl;

	  std::string sys_prefix = "--solver_group_";

	  if (libMesh::on_command_line("--solver_system_names"))
	    {
	      sys_prefix = sys_prefix + sys->name() + "_";
	    }

	  std::vector<dof_id_type> indices;


	  // if just a 3d simulation then we only need to check the outer,
	  // otherwise we need to 
	  if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 0)
	  {
		  if (libMesh::on_command_line("--solver_variable_names"))
		    {
		      for (unsigned int v = 0; v != sys->n_vars(); ++v)
			{
			  const std::string& var_name = sys->variable_name(v);

			  std::vector<dof_id_type> var_idx;
			  sys->get_dof_map().local_variable_indices
			    (var_idx, sys->get_mesh(), v);

			  // the outer fieldsplit grouping
			  std::string group_command = sys_prefix + var_name;

			  const std::string empty_string;

							// JAMES EDIT:
							// changed so that actually gets the command line bit
			  std::string group_name = libMesh::command_line_next
			    (group_command, empty_string);

			  if (group_name == field_name)
			    {
			      const bool prior_indices = !indices.empty();
			      indices.insert(indices.end(), var_idx.begin(),
				             var_idx.end());
			      if (prior_indices)
				std::sort(indices.begin(), indices.end());
			    }
			}
		    }
	


	  }
	  else
	  {
		// do this so that they know what preconditioner to use for the sub systems
		PetscErrorCode ierr=0;


		// now we want to do the fieldsplit of the sub fieldsplit

		std::cout << "Setting up the IS for the sub fieldsplits." << std::endl;

		PCType pc_type;
		ierr = PCGetType(my_pc,&pc_type);
		std::cout << "outer pc type = " << pc_type << std::endl;
		std::string fieldsplit_string = "fieldsplit";
		if (pc_type == fieldsplit_string)
		{
			std::cout << "yeah..." << std::endl;
			// loop over the number of dofs, cause this is the possible amount of fieldsplits
			KSP* subksp;	//subksps
			PetscInt num_splits;
			ierr = PCFieldSplitGetSubKSP(my_pc,&num_splits,&subksp);// CHKERRQ(ierr);

			// only need to check the first split
			PC sub_pc;
			ierr = KSPGetPC(subksp[0],&sub_pc);// CHKERRQ(ierr);
			PCType sub_pc_type;
			ierr = PCGetType(sub_pc,&sub_pc_type);
			std::cout << "sub_pc type = " << sub_pc_type << std::endl;
			std::cout << "hmmm" << std::endl;

				// need to extract the pc corresponding to this 
			sys_prefix = sys_prefix + static_cast<std::ostringstream*>( &(std::ostringstream() << i) )->str() + "_";

			std::map<std::string, std::vector<dof_id_type> > group_indices;

			if (libMesh::on_command_line("--solver_variable_names"))
			{
				for (unsigned int v = 0; v != sys->n_vars(); ++v)
				{
					const std::string& var_name = sys->variable_name(v);

					std::vector<dof_id_type> var_idx;
					sys->get_dof_map().local_variable_indices
					(var_idx, sys->get_mesh(), v);

					// the outer fieldsplit grouping
					std::string group_command = sys_prefix + var_name;

					const std::string empty_string;

								// JAMES EDIT:
								// changed so that actually gets the command line bit
					std::string group_name = libMesh::command_line_next
					(group_command, empty_string);


					if (group_name == field_name)
					{
						const bool prior_indices = !indices.empty();
						indices.insert(indices.end(), var_idx.begin(),
							     var_idx.end());
						if (prior_indices)
							std::sort(indices.begin(), indices.end());
					}
				}
			}
		}
	  }

		std::cout << "hi.."<< std::endl;





	  // ****************** CREATE THE IS *************************** //

		//std::cout << "giving the indices to fieldsplit" << std::endl;
	  const PetscInt *idx = PETSC_NULL;
	  if (!indices.empty())
	    idx = reinterpret_cast<const PetscInt*>(&indices[0]);

		//std::cout << "fieldname = " << field_name << std::endl;

	  if(i==0)
	  {
	    std::cout << "num dofs " << i << " = " << indices.size() << std::endl;
	    int ierr = ISCreateLibMesh(sys->comm().get(), indices.size(),
		                     idx, PETSC_COPY_VALUES, &velocity_is);
	  }
	  else
	  {
	    std::cout << "num dofs " << i << " = " << indices.size() << std::endl;
	    int ierr = ISCreateLibMesh(sys->comm().get(), indices.size(),
		                     idx, PETSC_COPY_VALUES, &pressure_is);

	  }
  }
 

}




// ********************** COMPUTE AND OUTPUT EIGENVALUES **************************** //
// compute and output the eigenvalues of the preconditioned operator to file
// can only handle a small number of unknowns like 500.
// if we are more than 500 then
int NavierStokesCoupled::compute_and_output_eigenvalues(TransientLinearImplicitSystem * system)
{

  PetscErrorCode ierr;

	//get petsc solver
	PetscLinearSolver<Number>* system_linear_solver =
			libmesh_cast_ptr<PetscLinearSolver<Number>* >
			(system->linear_solver.get());
	KSP system_ksp = system_linear_solver->ksp();

	// check the number of dofs
	dof_id_type num_dofs = system->n_dofs();

	if(num_dofs < 1000)
	{
		std::cout << "Calculating all eigenvalues for preconditioned system." << std::endl;
		PetscInt nmax = num_dofs;
		PetscReal r[num_dofs];
		PetscReal c[num_dofs];		
		
		std::cout << "Computing eigenvalues..." << std::endl;
		ierr = KSPComputeEigenvaluesExplicitly(system_ksp,nmax,r,c);

		std::cout << "Printing eigenvalues:" << std::endl;
		for(unsigned int i=0; i<num_dofs; i++)
		{
			std::cout << " eig " << i << " = " << r[i] << " + " << c[i] << "i" << std::endl;
		}
	}
	else
	{
		std::cout << "Calculating eigenvalues for large system not implemented yet." << std::endl;
		//std::cout << "Calculating extremum eigenvalues for preconditioned system." << std::endl;
		//ierr = KSPComputeEigenvalues(KSP ksp,PetscInt n,PetscReal r[],PetscReal c[],PetscInt *neig)
	}

	

}










// ************************* ADAPTIVELY REFINE 3D MESH*********************** //

// obv not working for monolithic or much really
void NavierStokesCoupled::adaptively_refine()
{
	
	std::cout << "  Adaptively refining the mesh..." << std::endl;

	//kelly error estimator
	if(es->parameters.get<unsigned int> ("error_estimator") == 1 ||
			es->parameters.get<unsigned int> ("error_estimator") == 2||
			es->parameters.get<unsigned int> ("error_estimator") == 3)
	{
  	picard->estimate_error (error_vector);
//		error_estimator->reduce_error(error_vector,mesh.comm());
		std::vector<float> grr_vector = error_vector;
		mesh.comm().sum(grr_vector);
		for(unsigned int i=0; i<error_vector.size(); i++)
			error_vector[i] = grr_vector[i];
//		error_vector = grr_vector;
	}
	else if(es->parameters.get<unsigned int> ("error_estimator") == 4)
	{

		error_estimator.reset(new KellyErrorEstimator);
  	error_estimator->estimate_error(*system_3d, error_vector);
	}

	
 	mesh_refinement.flag_elements_by_nelem_target(error_vector);

	// l2 brute force error estimator/ uniformly refines all elements and
	// calculates error.
	if(false)
	{			
		UniformRefinementEstimator *u =
                  			new UniformRefinementEstimator;	
		u->error_norm = L2;
  
                error_estimator.reset(u);
		std::vector<Real> weights(3,1.0);  // u, v
		weights.push_back(1.0);
		weights.push_back(1.0);
		if(threed)
			weights.push_back(1.0);
  	weights.push_back(0.0);            // p

		std::vector<FEMNormType>
	    		norms(1, error_estimator->error_norm.type(0));
  	error_estimator->error_norm = SystemNorm(norms, weights);

		std::cout << "hey solving system again to estimate error..." << std::endl;
  	error_estimator->estimate_error(*system_3d, error_vector);
		
		Real global_error = error_vector.l2_norm();
		std::cout << "Global_error = " << global_error
            		<< std::endl;

		std::cout << "Worst element error = " << error_vector.maximum()
							<< ", best element error = " << error_vector.minimum()
            		<< ", mean = " << error_vector.mean() << std::endl;

		//mesh_refinement.flag_elements_by_error_tolerance(error);
		// we want to limit the number of elements that can be refined
	    	//mesh_refinement.flag_elements_by_nelem_target(error);
		mesh_refinement.flag_elements_by_mean_stddev(error_vector,
			mesh_refinement.refine_fraction(),mesh_refinement.coarsen_fraction());
		
	}

	/*
	// my own flagging to flag the top elements for refinement
	if(false)
	{			

		double z_upper = 1.5;
		double z_lower = 1.0;
		ZPositionElementFlagging element_flagger(z_upper,z_lower);

	      	mesh_refinement.flag_elements_by(element_flagger);
		
	}
	*/
	mesh_refinement.refine_and_coarsen_elements();

	es->reinit ();
	init_dof_variable_vectors();

	refine_mesh = false;
	es->parameters.set<bool>("mesh_refined") = true;
	mesh.print_info();

	Real global_error = error_vector.l2_norm();

	system_3d->update();

	// if we have done adaptive refinement then we always need to
	// make sure the old global solution is the correct size (and vals)
 	old_global_solution->init(*system_3d->solution);
	*old_global_solution = *system_3d->solution;	// remember we changed solution to old global before refinement

	std::cout << "Global_error = " << global_error
        				<< std::endl;
	std::cout << "Worst element error = " << error_vector.maximum()
							<< ", best element error = " << error_vector.minimum()
        				<< ", mean = " << error_vector.mean() << std::endl;
	
	std::cout << "  Adaptive refinement done." << std::endl;


}

// ************************************************************************** //









// *********** PLOT THE ERROR ASSOCIATED WITH THE 3D SOLUTION *************** //

void NavierStokesCoupled::plot_error(ExodusII_IO_Extended& io)
{

  AutoPtr<MeshBase> meshptr = mesh.clone();
  MeshBase &temp_mesh = *meshptr;
  temp_mesh.all_first_order();
  EquationSystems temp_es (temp_mesh);
  ExplicitSystem& error_system
    = temp_es.add_system<ExplicitSystem> ("Error");
  error_system.add_variable("error", CONSTANT, MONOMIAL);
  temp_es.init();

  const DofMap& error_dof_map = error_system.get_dof_map();

  MeshBase::const_element_iterator       el     =
    temp_mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    temp_mesh.active_local_elements_end();
  std::vector<dof_id_type> dof_indices;

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    error_dof_map.dof_indices(elem, dof_indices);

    const dof_id_type elem_id = elem->id();

    //0 for the monomial basis
    const dof_id_type solution_index = dof_indices[0];

    // libMesh::out << "elem_number=" << elem_number << std::endl;
    libmesh_assert_less (elem_id, error_vector.size());

    // We may have zero error values in special circumstances
    // libmesh_assert_greater ((*this)[elem_id], 0.);
    error_system.solution->set(solution_index, error_vector[elem_id]);
  }

  // We may have to renumber if the original numbering was not
  // contiguous.  Since this is just a temporary mesh, that's probably
  // fine.
  if (temp_mesh.max_elem_id() != temp_mesh.n_elem() ||
      temp_mesh.max_node_id() != temp_mesh.n_nodes())
  {
    temp_mesh.allow_renumbering(true);
    temp_mesh.renumber_nodes_and_elements();
  }

  //ExodusII_IO io(mesh);
  //io.write(filename);
  //io.write_element_data(temp_es,time);
  io.write_element_data(temp_es);

}


// ************************************************************************** //



// *********** PLOT THE ERROR ASSOCIATED WITH THE 3D SOLUTION *************** //

void NavierStokesCoupled::write_elem_pid_3d(ExodusII_IO_Extended& io)
{

	std::cout << "hi" << std::endl;
  AutoPtr<MeshBase> meshptr = mesh.clone();
	std::cout << "hi" << std::endl;
  MeshBase &temp_mesh = *meshptr;
	temp_mesh.partitioner()->set_custom_partitioning(es->parameters.set<bool>("custom_partitioning"));
	std::cout << "hi" << std::endl;
  temp_mesh.all_first_order();
	std::cout << "hi" << std::endl;
	std::cout << "hi" << std::endl;
  EquationSystems temp_es (temp_mesh);
  ExplicitSystem& processor_system
    = temp_es.add_system<ExplicitSystem> ("Processor");
  processor_system.add_variable("proc_id", CONSTANT, MONOMIAL);
  temp_es.init();

  const DofMap& processor_dof_map = processor_system.get_dof_map();

  MeshBase::const_element_iterator       el     =
    temp_mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    temp_mesh.active_local_elements_end();
  std::vector<dof_id_type> dof_indices;

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;;
		if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
		{

		  processor_dof_map.dof_indices(elem, dof_indices);

		  //0 for the monomial basis
		  const dof_id_type solution_index = dof_indices[0];

		  // We may have zero error values in special circumstances
		  // libmesh_assert_greater ((*this)[elem_id], 0.);
		  processor_system.solution->set(solution_index, elem->processor_id());
		}
  }

  // We may have to renumber if the original numbering was not
  // contiguous.  Since this is just a temporary mesh, that's probably
  // fine.
  if (temp_mesh.max_elem_id() != temp_mesh.n_elem() ||
      temp_mesh.max_node_id() != temp_mesh.n_nodes())
  {
    temp_mesh.allow_renumbering(true);
    temp_mesh.renumber_nodes_and_elements();
  }

	std::cout << "about to write this shit" << std::endl;
  //ExodusII_IO io(mesh);
  //io.write(filename);
  //io.write_element_data(temp_es,time);
  io.write_element_data(temp_es);
	std::cout << "written" << std::endl;

}


// ************************************************************************** //



// convert the 1d part of a vector from nodal to monomial
void NavierStokesCoupled::scale_3d_solution_vector(double velocity_scaling=1.0,double pressure_scaling=1.0)
{

	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(sim_type == 5)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}

  // Numeric ids corresponding to each variable in the system
  int u_var = system->variable_number ("u");
  int v_var = system->variable_number ("v");
	int w_var = -1;
	if(threed)
  	w_var = system->variable_number ("w");
  int p_var = system->variable_number ("p");

	std::cout << "velocity_scaling = " << velocity_scaling << std::endl;
	std::cout << "pressure_scaling = " << pressure_scaling << std::endl;
	std::cout << "norm is inside before = " << system->solution->l2_norm() << std::endl;

	double value = 0.0;

	// need to make the solution vector global
  std::vector<Number> sys_soln;
  system->update_global_solution (sys_soln);

	for(unsigned int i=0; i<dof_variable_type_3d.size(); i++)
	{
		if(dof_variable_type_3d[i] == u_var
			 || dof_variable_type_3d[i] == v_var
			 || dof_variable_type_3d[i] == w_var)
		{
			value = sys_soln[i];
			system->solution->set(i,velocity_scaling * value);
		}
		else if(dof_variable_type_3d[i] == p_var)
		{
			value = sys_soln[i];
			system->solution->set(i,pressure_scaling * value);
		}
			
	}
	

	for(unsigned int i=0; i<dof_variable_type_coupled.size(); i++)
	{
		if(dof_variable_type_coupled[i] == u_var
			 || dof_variable_type_coupled[i] == v_var
			 || dof_variable_type_coupled[i] == w_var)
		{
			value = sys_soln[i];
			system->solution->set(i,velocity_scaling * value);
		}
		else if(dof_variable_type_coupled[i] == p_var)
		{
			value = sys_soln[i];
			system->solution->set(i,pressure_scaling * value);
		}
			
	}


	system->solution->close();

	std::cout << "norm is inside after = " << system->solution->l2_norm() << std::endl;
}




