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

	// add variable for proc id
	int proc_id_var = extra_3d_data_system->add_variable("proc_id", CONSTANT, MONOMIAL,&active_subdomains);


	// ******** WE NEED TO ADD THE PRECONDITIONER MATRIX POSSIBLY ********** //

	if(es->parameters.get<bool>("assemble_pressure_laplacian_matrix"))
		system->add_matrix("Pressure Laplacian Matrix");

	if(es->parameters.get<bool>("assemble_pressure_mass_matrix")
		|| es->parameters.get<bool>("assemble_scaled_pressure_mass_matrix"))
		system->add_matrix("Pressure Mass Matrix");


	if(es->parameters.get<bool>("assemble_pressure_convection_diffusion_matrix"))
		system->add_matrix("Pressure Convection Diffusion Matrix");

	if(es->parameters.get<bool>("assemble_velocity_mass_matrix"))
		system->add_matrix("Velocity Mass Matrix");

	if(es->parameters.get<unsigned int>("preconditioner_type_3d") || es->parameters.get<unsigned int>("preconditioner_type_3d1d"))
	{
		system->add_matrix("Preconditioner");
		system->add_matrix("Velocity Matrix");	// a full matrix with only the velocity block
	}


	if(es->parameters.get<bool>("moghadam_coupling"))
	{
		// need false flag to zero the vector
		system->add_vector("Moghadam Vector",false);
		system->add_vector("Moghadam Vector BC",false);
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



	// Renumbering
	if (es->parameters.get<bool> ("renumber_nodes_and_elements"))
	{
		std::cout << "\n\n\n\n\n\n\n\n\n" << std::endl;
		std::cout << "RENUMBERING" << std::endl;
		
	  mesh.allow_renumbering(true);
	  mesh.renumber_nodes_and_elements();
	}


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

	variables_3d.push_back("proc_id");

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
	perf_log.push("output_normal");
	exo.write_timestep(file_name_soln.str(), *es,1,time*time_scale_factor);
	perf_log.pop("output_normal");


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

	//std::cout << "before write pid" << std::endl;
	//write_elem_pid_3d(exo);




	// We write the equation systems and mesh object in .xda format
	// for restarting... perhaps we don't wanna do this every time step
	// but every n time_steps

	std::cout << "EXODUSII output for timestep " << t_step
		<< " written to " << file_name_soln.str() << std::endl;


	perf_log.push("output_backup");
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
	perf_log.pop("output_backup");

 
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
	// only use ksp view in the first iteration
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


	// need to zero the moghadam vectors in between each iteration
	if(es->parameters.get<bool>("moghadam_coupling"))
	{
		system->request_vector("Moghadam Vector")->close();
		system->request_vector("Moghadam Vector")->zero();
		system->request_vector("Moghadam Vector BC")->close();
		system->request_vector("Moghadam Vector BC")->zero();

	}

	// need to recollect the ksp because perhaps some shit changed in restrict_solve_to
	system_linear_solver = libmesh_cast_ptr<PetscLinearSolver<Number>* > (system->linear_solver.get());
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
	
	if(es->parameters.get<bool>("stokes") && sim_type != 3 )
	{
		std::cout << "Stokes so only one iteration required" << std::endl;
		return true;
	}

	if(es->parameters.get<unsigned int>("newton") == 3 && sim_type != 3)
	{
		std::cout << "Semi-implicit so only one iteration required" << std::endl;
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
	else if(nonlinear_iteration >= es->parameters.get<unsigned int>("max_newton_iterations"))
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

	std::cout << "before setup" << std::endl;


	if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 0 && es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 0)
		system_linear_solver->solve_simple_setup (*system->request_matrix("System Matrix"), *system->solution, *system->rhs, tol, maxits);
	else
		system_linear_solver->solve_simple_setup (*system->request_matrix("System Matrix"), *system->request_matrix("Preconditioner"), *system->solution, *system->rhs, tol, maxits);


//	PetscInt numsplit;
//	ierr = PCFieldSplitGetSubKSP(sub_pc,&numsplit,NULL);// CHKERRQ(ierr);
//	std::cout << "hmmm" << std::endl;
//	std::cout << "numsplit = " << numsplit << std::endl;
//	std::cout << "hmmm" << std::endl;

	std::cout << "after setup" << std::endl;

	//ierr = KSPSetUp(system_ksp); CHKERRQ(ierr);

	if(nonlinear_iteration == 1 && es->parameters.get<bool>("fieldsplit"))
		system->linear_solver->init_names(*system);

	std::cout << "after init names" << std::endl;

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
		//PetscErrorCode (*function_ptr)(KSP, PetscInt, PetscReal, void*);
		//function_ptr = &custom_outer_monitor;

		//ierr = KSPMonitorSet(system_ksp,function_ptr,NULL,NULL); CHKERRQ(ierr);
	}

	//get the inner iteration ksp
	//int num_splits = 0;	// unused


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

	std::cout << "after setup preconditioner" << std::endl;


	if(es->parameters.get<bool>("ksp_view_before"))
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
	//bool _shell_matrix = false;	// unused
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

	// first time step or zeroth time step (setting up IC)
	if(es->parameters.get<bool>("ksp_view") && nonlinear_iteration == 1 && (t_step - es->parameters.get<unsigned int>("restart_time_step")) <= 1)
	{
		// we need to set this directly
		ierr = KSPView(system_ksp,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
	}









	if(es->parameters.get<bool>("post_solve"))
		test_post_solve(system);




















	int num_outer_its = 0;	
	//int num_inner_velocity_its = 0;	// unused
	//int num_inner_pressure_its = 0;	// unused
	// output some data about the iterative solve
	// this is hard to do when doing restrict_solve_to, so just leave it

	ierr = KSPGetIterationNumber(system_ksp,&num_outer_its); CHKERRQ(ierr);

	if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") >= 6
		&& es->parameters.get<unsigned int>("preconditioner_type_3d") >= 2)
	{
		num_outer_its = mono_ctx->total_velocity_iterations;
	}
	
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
	std::cout << "outer its this nlin it = " << num_outer_its << std::endl;
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
	// only the pressure convection diffusion matrix and the pressure laplacian change at each step the others we don't need to delete
	if(es->parameters.get<bool>("assemble_pressure_convection_diffusion_matrix"))
	{
		delete pressure_convection_diffusion_matrix;
		system->request_matrix("Pressure Convection Diffusion Matrix")->zero();
	}

	if(es->parameters.get<bool>("assemble_pressure_laplacian_matrix"))
	{
		delete pressure_laplacian_matrix;
		system->request_matrix("Pressure Laplacian Matrix")->zero();
	}

	// delete the shell matrix from pcd 2.0
	if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 5)
	{
		ierr = MatDestroy(&B); CHKERRQ(ierr);
	}


	// the navier stokes 3d1d preconditioners need to be deleted and zeroed
	if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 8 
		|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 9)
	{

		std::cout << "hiya" << std::endl;
		delete velocity_matrix;
		std::cout << "hmmm" << std::endl;

		// also zero the assemble bits for next time
		system->request_matrix("Velocity Matrix")->close();
		system->request_matrix("Velocity Matrix")->zero();
	}


	// the navier stokes 3d1d preconditioners need to be deleted, not zeroed cause not reassembled
	if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 10 
		|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 11
		|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 12
		|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 6)
	{

		std::cout << "hiya" << std::endl;
		delete velocity_matrix;
		std::cout << "hmmm" << std::endl;

		// also zero the assemble bits for next time
		system->request_matrix("Velocity Matrix")->close();
		system->request_matrix("Velocity Matrix")->zero();
	}

	// the navier stokes 3d1d preconditioners need to be deleted and zeroed
	if((es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 10)
		&& es->parameters.get<bool>("multiple_column_solve"))
	{

		std::cout << "hiya" << std::endl;
		VecDestroy(&non_zero_rows);
		VecDestroy(&non_zero_cols);
		std::cout << "hmmm" << std::endl;

	}




	// delete the IS's we created in the SIMPLE paradigm
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

// ************************************************************************* //












// ********************** SETUP PRECONDITIONERS **************************** //
// setup the preconditioner for use in the solver
// construct submatrices from full matrices that they were assembled as.
int NavierStokesCoupled::setup_preconditioners(TransientLinearImplicitSystem * system, Mat& B)
{


	PetscErrorCode ierr;

	std::cout << "Setting up submatrices for preconditioners." << std::endl;

	//std::exit(0);

	if(es->parameters.get<unsigned int>("preconditioner_type_3d") != 0 
		|| es->parameters.get<unsigned int>("preconditioner_type_schur_stokes") != 0)
	{


	  	
		// ********* CONSTRUCT SUBMATRICES ********** //

		std::cout << "Constructing Navier Stokes Preconditioner submatrices." << std::endl;



		// set up variable dofs to set up submatrices
		const unsigned int u_var = system->variable_number ("u");
		const unsigned int v_var = system->variable_number ("v");
		unsigned int w_var = 0;
		if(threed)
			w_var = system->variable_number ("w");
		const unsigned int p_var = system->variable_number ("p");

		std::vector<dof_id_type> 
			U_var_idx, u_var_idx, v_var_idx, w_var_idx, p_var_idx;

		system->get_dof_map().local_variable_indices
			(u_var_idx, system->get_mesh(), u_var);
		system->get_dof_map().local_variable_indices
			(v_var_idx, system->get_mesh(), v_var);

		if(threed)
			system->get_dof_map().local_variable_indices
				(w_var_idx, system->get_mesh(), w_var);

		system->get_dof_map().local_variable_indices
			(p_var_idx, system->get_mesh(), p_var);

		// concatenate all velocity dofs
		U_var_idx.insert(U_var_idx.end(),u_var_idx.begin(),u_var_idx.end());
		U_var_idx.insert(U_var_idx.end(),v_var_idx.begin(),v_var_idx.end());
		U_var_idx.insert(U_var_idx.end(),w_var_idx.begin(),w_var_idx.end());

		std::sort(U_var_idx.begin(),U_var_idx.end());

	
	



		// object to enable us to pin the matrices at a pressure node if necessary
		std::vector<dof_id_type> rows;
		rows.push_back(0);



		std::cout << "hello" << std::endl << std::flush;

		// the pressure mass matrix, scaled mass matrix and pressure laplacian only need to be constructed once,
		// while the convection diffusion matrix needs to be computed at each step.

		if(!shell_pc_created)
		{

			if(es->parameters.get<bool>("assemble_pressure_mass_matrix"))
			{
				// ************ pressure mass matrix ************ //
				pressure_mass_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
				system->request_matrix("Pressure Mass Matrix")->create_submatrix(*pressure_mass_matrix,p_var_idx,p_var_idx);
				//if(!es->parameters.get<unsigned int>("problem_type") == 4 && !es->parameters.get<unsigned int>("problem_type") == 5)
				//pressure_mass_matrix->zero_rows(rows,1.0);
				pressure_mass_matrix->close();

				// now that this matrix has been set we can not assemble it again
				es->parameters.set<bool>("assemble_pressure_mass_matrix") = false;
			}


			if(es->parameters.get<bool>("assemble_scaled_pressure_mass_matrix"))
			{

				// ************ pressure mass matrix ************ //
				pressure_mass_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
				system->request_matrix("Pressure Mass Matrix")->create_submatrix(*pressure_mass_matrix,p_var_idx,p_var_idx);
				//if(!es->parameters.get<unsigned int>("problem_type") == 4 && !es->parameters.get<unsigned int>("problem_type") == 5)
				//pressure_mass_matrix->zero_rows(rows,1.0);
				pressure_mass_matrix->close();

				std::cout << "1" << std::endl;
				//************* scaled pressure mass matrix ************** //
				scaled_pressure_mass_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());

				std::cout << "1" << std::endl;
				// to initialise it
				system->request_matrix("Pressure Mass Matrix")->create_submatrix(*scaled_pressure_mass_matrix,p_var_idx,p_var_idx);

				std::cout << "1" << std::endl;
				scaled_pressure_mass_matrix->zero();
				scaled_pressure_mass_matrix->add(es->parameters.get<Real>("reynolds_number"),*pressure_mass_matrix);
				//if(!es->parameters.get<unsigned int>("problem_type") == 4 && !es->parameters.get<unsigned int>("problem_type") == 5)
				//scaled_pressure_masNus_matrix->zero_rows(rows,1.0);

				std::cout << "1" << std::endl;
				scaled_pressure_mass_matrix->close();


				// now that this matrix has been set we can not assemble it again
				es->parameters.set<bool>("assemble_scaled_pressure_mass_matrix") = false;
			}
	
			std::cout << "k?" << std::endl;


			if(es->parameters.get<bool>("assemble_velocity_mass_matrix"))
			{
				// *********** velocity mass matrix ********************************** //
				velocity_mass_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
				system->request_matrix("Velocity Mass Matrix")->create_submatrix(*velocity_mass_matrix,U_var_idx,U_var_idx);
				velocity_mass_matrix->close();

				std::cout << " yeah baby" << std::endl << std::flush;

				// now that this matrix has been set we can not assemble it again
				es->parameters.set<bool>("assemble_velocity_mass_matrix") = false;
			}
			std::cout << "k?" << std::endl << std::flush;
		}


			std::cout << "k babe?" << std::endl << std::flush;
		if(es->parameters.get<bool>("assemble_pressure_laplacian_matrix"))
		{
			// ************ pressure laplacian matrix *************** //
			pressure_laplacian_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
			// need to pin laplacian matrix... at the first pressure node
			system->request_matrix("Pressure Laplacian Matrix")->create_submatrix(*pressure_laplacian_matrix,p_var_idx,p_var_idx);
			if(es->parameters.get<bool>("pin_pressure"))
			{
				pressure_laplacian_matrix->zero_rows(rows,1.0);
			}
			pressure_laplacian_matrix->close();

		}
			std::cout << "k?" << std::endl << std::flush;

			//std::exit(0);
		if(es->parameters.get<bool>("assemble_pressure_convection_diffusion_matrix"))
		{
			// ************ pressure convection diffusion matrix ***************** //
			pressure_convection_diffusion_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
			system->request_matrix("Pressure Convection Diffusion Matrix")->create_submatrix(*pressure_convection_diffusion_matrix,p_var_idx,p_var_idx);
			if(es->parameters.get<bool>("pin_pressure"))
			{
				pressure_convection_diffusion_matrix->zero_rows(rows,1.0/es->parameters.get<Real>("reynolds_number"));
			}
			pressure_convection_diffusion_matrix->close();
		}

	
			std::cout << "k?" << std::endl << std::flush;




		//std::exit(0)



		// set pressure in system matrix and preconditioner
		std::vector<dof_id_type> rows_full;
		// rows_full.push_back(p_var_idx[0]);		//doesn't do anything, but doesn't work cause need to only use one pressure dof for the whole system


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

		std::cout << "Done setting up Navier Stokes Preconditioner sub matrices." << std::endl;


		// *********************************************************************** //


	}

	

	// give all monolithic preconditioners the mono_ctx for counting
	if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") >= 6)
	{


	  	PetscErrorCode ierr;

		//get petsc solver
		PetscLinearSolver<Number>* system_linear_solver =
				libmesh_cast_ptr<PetscLinearSolver<Number>* >
				(system->linear_solver.get());
		KSP system_ksp = system_linear_solver->ksp();
		PC pc;
		ierr = KSPGetPC(system_ksp,&pc); CHKERRQ(ierr);



		std::cout << "Setting up Diagonal Preconditioner." << std::endl;
		int num_splits = 0;
		ierr = KSPSetUp(system_ksp); CHKERRQ(ierr);

		//KSP* subksp;	//subksps
		KSP* subksp;	//subksps
		ierr = PCFieldSplitGetSubKSP(pc,&num_splits,&subksp); CHKERRQ(ierr);


		if(es->parameters.get<bool>("nonzero_initial_guess_inner_3d1d") && !es->parameters.get<bool>("monolithic_navier_stokes_preonly"))
		{
			ierr = KSPSetInitialGuessNonzero(subksp[0],PETSC_TRUE); CHKERRQ(ierr);
		}




		if(es->parameters.get<unsigned int>("preconditioner_type_3d") >= 2)
		{
			PetscErrorCode (*function_ptr)(KSP, PetscInt, PetscReal, void*);
			function_ptr = &custom_outer_monitor;

			mono_ctx->total_velocity_iterations = 0;
			mono_ctx->total_convection_diffusion_iterations = 0;

			ierr = KSPMonitorSet(system_ksp,function_ptr,mono_ctx,NULL); CHKERRQ(ierr);
		}

	}

	// note we don't need to do anything for the schur 0D (12) preconditioner as this is all specified in the command line
	if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 8
			|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 9
			|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 10
			|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 11
			|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 12
			|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 6)
		{

			// now set up the shell preconditioner

			//get petsc solver
			PetscLinearSolver<Number>* system_linear_solver =
					libmesh_cast_ptr<PetscLinearSolver<Number>* >
					(system->linear_solver.get());
			KSP system_ksp = system_linear_solver->ksp();
			PC pc;
			ierr = KSPGetPC(system_ksp,&pc); CHKERRQ(ierr);



			std::cout << "Setting up Velocity Monolithic preconditioner." << std::endl;
	

			//ierr = KSPSetUp(system_ksp); CHKERRQ(ierr);	// maybe
			KSP* subksp;	//subksps
			ierr = PCFieldSplitGetSubKSP(pc,NULL,&subksp); CHKERRQ(ierr);





			// set up some stuff for the velocity preconditioner
			// need to extract the navier stokes block and add ones on the diagonal of the pressure bit


			std::cout << "Constructing submatrix for monolithic schur complement approximation." << std::endl;
			const unsigned int u_var = system->variable_number ("u");
			const unsigned int v_var = system->variable_number ("v");
			unsigned int w_var = 0;
			if(threed)
				w_var = system->variable_number ("w");
			const unsigned int p_var = system->variable_number ("p");

			std::vector<dof_id_type> 
				NS_var_idx,	u_var_idx, v_var_idx, w_var_idx, p_var_idx;


			// get the NS dofs
			system->get_dof_map().local_variable_indices
				(u_var_idx, system->get_mesh(), u_var);
			system->get_dof_map().local_variable_indices
				(v_var_idx, system->get_mesh(), v_var);

			if(threed)
				system->get_dof_map().local_variable_indices
					(w_var_idx, system->get_mesh(), w_var);

			system->get_dof_map().local_variable_indices
				(p_var_idx, system->get_mesh(), p_var);


			// concatenate all velocity dofs and then sort
			NS_var_idx.insert(NS_var_idx.end(),u_var_idx.begin(),u_var_idx.end());
			NS_var_idx.insert(NS_var_idx.end(),v_var_idx.begin(),v_var_idx.end());
			NS_var_idx.insert(NS_var_idx.end(),w_var_idx.begin(),w_var_idx.end());
			NS_var_idx.insert(NS_var_idx.end(),p_var_idx.begin(),p_var_idx.end());

			std::sort(NS_var_idx.begin(),NS_var_idx.end());

			// extract the correct matrix
			// for schur stokes (10) we need the stokes matrix that is assembled into Velocity Matrix
			// for schur velocity (8) we need the velocity only matrix that is assembled into Velocity Matrix
			// for schur stokes velocity (11) we need the stokes velocity only matrix that is assembled into Velocity Matrix
			// for the other methods we just require the system matrix
			if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 8 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 10
				|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 11	|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 12
				|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 6)
			{
				velocity_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
				system->request_matrix("Velocity Matrix")->create_submatrix(*velocity_matrix,NS_var_idx,NS_var_idx);
				velocity_matrix->close();
			}
			else
			{

				velocity_matrix = cast_ptr<PetscMatrix<Number>*>(SparseMatrix<Number>::build(mesh.comm()).release());
				system->request_matrix("System Matrix")->create_submatrix(*velocity_matrix,NS_var_idx,NS_var_idx);
				velocity_matrix->close();
			}
	
			





			// set up a global Vec saying which columns have entries in the B1 matrix
			// remember we only need to do this once, hence the !mono_shell_pc_created

			if(es->parameters.get<bool>("multiple_column_solve") && es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 10)
			{
				std::cout << "Identifying rows and columns of submatrix involved in monolithic schur complement." << std::endl;

				std::vector<dof_id_type> pressure_coupling_dofs;
				std::vector<dof_id_type> flux_coupling_dofs;
				// set up some
				const unsigned int Q_var = system->variable_number ("Q");
				const unsigned int P_var = system->variable_number ("P");

				std::vector<dof_id_type> 
					Q_var_idx, P_var_idx;

				system->get_dof_map().local_variable_indices
					(Q_var_idx, system->get_mesh(), Q_var);
				system->get_dof_map().local_variable_indices
					(P_var_idx, system->get_mesh(), P_var);

			//PetscInt n_col = Q_var_idx.size() + P_var_idx.size();	// unused

				const MeshBase& mesh = es->get_mesh();
				const DofMap & dof_map = system->get_dof_map();

				std::vector<dof_id_type> dof_indices;
				std::vector<dof_id_type> dof_indices_p;
				std::vector<dof_id_type> dof_indices_q;



				MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
				const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
						

				// loop over the 1d elements and figure out what the P and Q dofs on the coupling boundary are
				// put these dofs into pressure_coupling_dofs and flux_coupling_dofs
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

									// Get the pressure and flow dof indices for the current element.
									dof_map.dof_indices (elem, dof_indices_p, P_var);
									pressure_coupling_dofs.push_back(dof_indices_p[0]);

									dof_map.dof_indices (elem, dof_indices_q, Q_var);
									flux_coupling_dofs.push_back(dof_indices_q[0]);

								}
							}
						}
					}
				}


	
				// now we have the global dof positions of the coupling dofs, we want to find where they occur in the local 1D block
				std::vector<dof_id_type> PQ_var_idx;
				PQ_var_idx.insert(PQ_var_idx.end(), P_var_idx.begin(), P_var_idx.end() );
				PQ_var_idx.insert(PQ_var_idx.end(), Q_var_idx.begin(), Q_var_idx.end() );
				std::sort(PQ_var_idx.begin(),PQ_var_idx.end());

				// set up the petsc Vec that is needed
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

				// put the values into the Vec
				ierr = VecAssemblyBegin(non_zero_cols); CHKERRQ(ierr);
				ierr = VecAssemblyEnd(non_zero_cols); CHKERRQ(ierr);

				ierr = VecAssemblyBegin(non_zero_rows); CHKERRQ(ierr);
				ierr = VecAssemblyEnd(non_zero_rows); CHKERRQ(ierr);


				if(nonlinear_iteration == 1 && t_step == 1)
				{
					construct_schur_stokes_matrix(system,subksp[1]);
				}


			}
		
			








			std::cout << "YES" << std::endl;






			if(es->parameters.get<bool>("nonzero_initial_guess_inner_3d1d") && !es->parameters.get<bool>("monolithic_navier_stokes_preonly"))
			{
				ierr = KSPSetInitialGuessNonzero(subksp[0],PETSC_TRUE); CHKERRQ(ierr);
			}

		
			if((es->parameters.set<double> ("last_nonlinear_iterate") < es->parameters.get<double>("monolithic_navier_stokes_preonly_switch"))
				&& nonlinear_iteration > 1)
			{
				std::cout << "changing to use preonly" << std::endl;
				ierr = KSPSetType(subksp[0], KSPPREONLY); CHKERRQ(ierr);
				ierr = KSPSetInitialGuessNonzero(subksp[0],PETSC_FALSE); CHKERRQ(ierr);
			}
			else if(nonlinear_iteration == 1 && !es->parameters.get<bool>("monolithic_navier_stokes_preonly"))
			{
				std::cout << "changing back to use gmres" << std::endl;
				ierr = KSPSetType(subksp[0], KSPGMRES); CHKERRQ(ierr);

			}


			// set ksp to preonly because we are defining the inverse
			ierr = KSPSetType(subksp[1], KSPPREONLY); CHKERRQ(ierr);

			PC schur_pc;
			ierr = KSPGetPC(subksp[1],&schur_pc); CHKERRQ(ierr);

			// think about this
			//ierr = KSPSetReusePreconditioner(subksp[0],PETSC_FALSE); CHKERRQ(ierr);
			// set to not reuse the preconditioner
			//ierr = KSPSetReusePreconditioner(subksp[1],PETSC_TRUE); CHKERRQ(ierr);


			// (Required) Indicate to PETSc that we're using a "shell" preconditioner 
			ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);

			// destroy old shell before it is created
			if(mono_shell_pc_created)
			{
				ierr = NSShellDestroy(schur_pc);
			}

			// (Optional) Create a context for the user-defined preconditioner; this
			//	context can be used to contain any application-specific data. 
			ierr = ShellPCCreate(&mono_shell); CHKERRQ(ierr);

			// (Required) Set the user-defined routine for applying the preconditioner 
			ierr = PCShellSetApply(schur_pc,MonolithicShellPCApply); CHKERRQ(ierr);
			ierr = PCShellSetContext(schur_pc,mono_shell); CHKERRQ(ierr);

			// (Optional) Set user-defined function to free objects used by custom preconditioner 
			ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

			// (Optional) Set a name for the preconditioner, used for PCView() 
			ierr = PCShellSetName(schur_pc,"Monolithic Schur Complement Preconditioner"); CHKERRQ(ierr);

			bool schur_0d = false;
			if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 12
					|| es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 6)
				schur_0d = true;

			if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 6)
				es->parameters.set<bool>("negative_mono_schur_complement") = false;	//because schur diag flips it the wrong way

			double scaling_factor = 1.0;
			if(es->parameters.get<unsigned int>("scale_mono_preconditioner"))
				scaling_factor = 1./es->parameters.get<double>("mono_preconditioner_resistance_scaling");



			std::cout << "YES" << std::endl;
			// Do setup of preconditioner
			if(es->parameters.get<bool>("multiple_column_solve"))
			{
				if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") != 10)
				{
					ierr = Monolithic3ShellPCSetUp(schur_pc,NULL,subksp[1],system_ksp,es->parameters.get<bool>("negative_mono_schur_complement"),schur_0d,scaling_factor); CHKERRQ(ierr);
				}
				else
				{
					ierr = Monolithic3ShellPCSetUp(schur_pc,schur_complement_approx,subksp[1],system_ksp,es->parameters.get<bool>("negative_mono_schur_complement"),schur_0d,scaling_factor); CHKERRQ(ierr);
				}
			}
			else if(es->parameters.get<bool>("multiple_column_solve"))
			{
				ierr = Monolithic2ShellPCSetUp(schur_pc,velocity_matrix->mat(),non_zero_cols,non_zero_rows,subksp[1]); CHKERRQ(ierr);
			}
			else
			{
				ierr = MonolithicShellPCSetUp(schur_pc,velocity_matrix->mat(),subksp[1]); CHKERRQ(ierr);
			}

			std::cout << "YES" << std::endl;


			if(es->parameters.get<unsigned int>("preconditioner_type_3d") >= 2)
			{
				PetscErrorCode (*function_ptr)(KSP, PetscInt, PetscReal, void*);
				function_ptr = &custom_outer_monitor;

				mono_ctx->total_velocity_iterations = 0;
				mono_ctx->total_convection_diffusion_iterations = 0;

				ierr = KSPMonitorSet(system_ksp,function_ptr,mono_ctx,NULL); CHKERRQ(ierr);
			}
			// let the method know we have calculated the shell preconditioner at least once before
			mono_shell_pc_created = true;

		}
	//}

	// all preconditioners except direct and straight gmres
	if(es->parameters.get<unsigned int>("preconditioner_type_3d") > 1)
	{


	  	PetscErrorCode ierr;

		//get petsc solver
		PetscLinearSolver<Number>* system_linear_solver =
				libmesh_cast_ptr<PetscLinearSolver<Number>* >
				(system->linear_solver.get());

		KSP system_ksp = system_linear_solver->ksp();
		PC pc;
		ierr = KSPGetPC(system_ksp,&pc); CHKERRQ(ierr);



		// set up variable dofs to set up submatrices
		const unsigned int u_var = system->variable_number ("u");
		const unsigned int v_var = system->variable_number ("v");
		unsigned int w_var = 0;
		if(threed)
			w_var = system->variable_number ("w");
		const unsigned int p_var = system->variable_number ("p");

		std::vector<dof_id_type> 
			U_var_idx, u_var_idx, v_var_idx, w_var_idx, p_var_idx;

		system->get_dof_map().local_variable_indices
			(u_var_idx, system->get_mesh(), u_var);
		system->get_dof_map().local_variable_indices
			(v_var_idx, system->get_mesh(), v_var);

		if(threed)
			system->get_dof_map().local_variable_indices
				(w_var_idx, system->get_mesh(), w_var);

		system->get_dof_map().local_variable_indices
			(p_var_idx, system->get_mesh(), p_var);

		// concatenate all velocity dofs
		U_var_idx.insert(U_var_idx.end(),u_var_idx.begin(),u_var_idx.end());
		U_var_idx.insert(U_var_idx.end(),v_var_idx.begin(),v_var_idx.end());
		U_var_idx.insert(U_var_idx.end(),w_var_idx.begin(),w_var_idx.end());

		std::sort(U_var_idx.begin(),U_var_idx.end());

	








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



			

			if(nonlinear_iteration == 1)
			{
				std::cout << "Setting up the convection diffusion Petsc settings" << std::endl;
				ierr = KSPSetOptionsPrefix(velocity_subksp[0],"convection_diffusion_"); CHKERRQ(ierr);

				// think about this
				//ierr = KSPSetReusePreconditioner(velocity_subksp[0],PETSC_FALSE); CHKERRQ(ierr);
				// set to not reuse the preconditioner
				if(es->parameters.get<bool>("reuse_convection_diffusion_pc"))
				{
					ierr = KSPSetReusePreconditioner(velocity_subksp[0],PETSC_TRUE); CHKERRQ(ierr);
				}
				else
				{
					ierr = KSPSetReusePreconditioner(velocity_subksp[0],PETSC_FALSE); CHKERRQ(ierr);
				}
				

				ierr = KSPSetFromOptions (velocity_subksp[0]); CHKERRQ(ierr);
			}


			std::cout << "Setting up the Schur preconditioner Petsc settings" << std::endl;


			ierr = KSPSetOptionsPrefix(velocity_subksp[1],"navier_stokes_schur_"); CHKERRQ(ierr);
			ierr = KSPSetFromOptions (velocity_subksp[1]); CHKERRQ(ierr);

			std::cout << "Old Preconditioner destroyed." << std::endl;


			std::cout << "Setting up Navier Stokes Schur complement preconditioners." << std::endl;

			// now set the preconditioner stuff for the schur preconditioners

			// scaled mass matrix with solve
			if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 1)
			{

				// note even though the preconditioner says it is using A11 it is not
				Mat Amat, Pmat;
				ierr = KSPGetOperators(velocity_subksp[1],&Amat,&Pmat); CHKERRQ(ierr);
				ierr = KSPSetOperators(velocity_subksp[1],Amat,scaled_pressure_mass_matrix->mat()); CHKERRQ(ierr);
		
		
		
			}
			// Petsc LSC preconditioner
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 2)
			{

				std::cout << "Setting up LSC preconditioner." << std::endl;

				// set ksp to preonly because we are defining the inverse
				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);
				ierr = PCSetType(schur_pc,PCLSC); CHKERRQ(ierr);

		
		
			}
			// Scaled Pressure Mass Matrix
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 3)	// shell pc..
			{

				std::cout << "Setting up Scaled Pressure Mass Matrix preconditioner." << std::endl;

				// set ksp to preonly because we are defining the inverse
				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);

				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);
				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

				ierr = PCShellSetApply(schur_pc,PressureShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);
				ierr = PCShellSetName(schur_pc,"Scaled Pressure Mass Matrix Preconditioner"); CHKERRQ(ierr);

				ierr = PressureShellPCSetUp(schur_pc,scaled_pressure_mass_matrix->mat(),velocity_subksp[1]); CHKERRQ(ierr);

				shell_pc_created = true;
			}
			// PCD
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 4)	// shell pc..
			{

				std::cout << "Setting up PCD preconditioner." << std::endl;

				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);

				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);

				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);
				ierr = PCShellSetApply(schur_pc,PCDShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);
				ierr = PCShellSetName(schur_pc,"PCD Preconditioner"); CHKERRQ(ierr);

				ierr = PCDShellPCSetUp(schur_pc,pressure_mass_matrix->mat(),pressure_laplacian_matrix->mat(),pressure_laplacian_matrix->mat(),pressure_convection_diffusion_matrix->mat(),velocity_subksp[1]); CHKERRQ(ierr);

				shell_pc_created = true;
			}
			// PCD 2.0
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 5)	// shell pc..+ shel matrix pcd2.0
			{

		
				std::cout << "Setting up PCD preconditioner." << std::endl;

				mat_ctx.m = p_var_idx.size();
				mat_ctx.n = p_var_idx.size();
				mat_ctx.pressure_laplacian_matrix = pressure_laplacian_matrix->mat();
				mat_ctx.velocity_mass_matrix = velocity_mass_matrix->mat();
				mat_ctx.schur_ksp = velocity_subksp[1];

				PetscInt mn = p_var_idx.size() * p_var_idx.size();

				std::cout << "m = " << mat_ctx.m  << std::endl;
				std::cout << "n = " << mat_ctx.n  << std::endl;
				std::cout << "mn = " << mn  << std::endl;
				ierr = MatCreateShell(PETSC_COMM_WORLD,mat_ctx.m,mat_ctx.n,PETSC_DETERMINE,PETSC_DETERMINE,&mat_ctx,&B);  CHKERRQ(ierr);
				ierr = MatShellSetOperation(B, MATOP_MULT, (void(*)(void))MatShellMultFull); CHKERRQ(ierr);
		
				//get operators from the outer pcshell solver
				Mat Amat, Pmat;
				ierr = KSPGetOperators(velocity_subksp[1],&Amat,&Pmat); CHKERRQ(ierr);

				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);

				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);


				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

				ierr = PCShellSetApply(schur_pc,PCD2ShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

				ierr = PCShellSetName(schur_pc,"PCD2 Preconditioner"); CHKERRQ(ierr);

				//ierr = PCDShellPCSetUp(schur_pc,pressure_mass_matrix->mat(),pressure_laplacian_matrix->mat(),pressure_laplacian_matrix->mat(),pressure_convection_diffusion_matrix->mat(),velocity_subksp[1]); CHKERRQ(ierr);
				ierr = PCD2ShellPCSetUp(schur_pc,pressure_mass_matrix->mat(),B,pressure_laplacian_matrix->mat(),pressure_convection_diffusion_matrix->mat(),velocity_subksp[1]); CHKERRQ(ierr);

				shell_pc_created = true;
			}
			// LSC James
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 6)	// shell pc..
			{

				std::cout << "Setting up LSC James preconditioner." << std::endl;

				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);
				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);


				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

				ierr = PCShellSetApply(schur_pc,LSCShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

				ierr = PCShellSetName(schur_pc,"LSC James Preconditioner"); CHKERRQ(ierr);

				ierr = LSCShellPCSetUp(schur_pc,velocity_subksp[1]); CHKERRQ(ierr);

				shell_pc_created = true;
			}
			// LSC James Scaled properly
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 7)	// shell pc..
			{

				std::cout << "Setting up LSC James scaled preconditioner." << std::endl;

				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);

				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);


				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

				ierr = PCShellSetApply(schur_pc,LSCScaledShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

				ierr = PCShellSetName(schur_pc,"LSC James Scaled Preconditioner"); CHKERRQ(ierr);

				ierr = LSCScaledShellPCSetUp(schur_pc,velocity_mass_matrix->mat(),velocity_subksp[1]); CHKERRQ(ierr);

				shell_pc_created = true;
			}
			// Petsc LSC command line (defunct)
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 8)
			{
				// do nothing, we should have set this all up on the command line
			}
			// LSC James scaled stabilised elements
			else if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 9)
			{
			
				std::cout << "Setting up LSC James scaled stabilised preconditioner." << std::endl;

				// set ksp to preonly because we are defining the inverse
				ierr = KSPSetType(velocity_subksp[1], KSPPREONLY); CHKERRQ(ierr);

				PC schur_pc;
				ierr = KSPGetPC(velocity_subksp[1],&schur_pc); CHKERRQ(ierr);

				ierr = PCSetType(schur_pc,PCSHELL); CHKERRQ(ierr);

				ierr = ShellPCCreate(&shell); CHKERRQ(ierr);

				ierr = PCShellSetApply(schur_pc,LSCScaledStabilisedShellPCApply); CHKERRQ(ierr);
				ierr = PCShellSetContext(schur_pc,shell); CHKERRQ(ierr);

				ierr = PCShellSetDestroy(schur_pc,NSShellDestroy); CHKERRQ(ierr);

				ierr = PCShellSetName(schur_pc,"LSC James Scaled Stabilised Preconditioner"); CHKERRQ(ierr);

				       
				ierr = LSCScaledStabilisedShellPCSetUp(schur_pc,velocity_mass_matrix->mat(),velocity_subksp[1],es->parameters.get<bool> ("negative_bfbt_schur_complement")); CHKERRQ(ierr);

				shell_pc_created = true;
			}
		}
		else
		{

			//int num_splits = 0;	// unused


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

	return 0;	
}

// ********************************************************************** //
















// ****************** CONSTRUCT THE SCHUR STOKES MATRIX ************** //

// the pressure coupling dofs and flux coupling dofs corresponding to the same 
int NavierStokesCoupled::construct_schur_stokes_matrix (TransientLinearImplicitSystem * sys, KSP schur_ksp)
{



	std::cout << "Constructing Schur Stokes Correction..." << std::endl;
	std::cout << "Could take a while..." << std::endl;


  	PetscErrorCode ierr;	// unused

	//get petsc solver
	PetscLinearSolver<Number>* system_linear_solver =
			libmesh_cast_ptr<PetscLinearSolver<Number>* >
			(sys->linear_solver.get());
	KSP system_ksp = system_linear_solver->ksp();
	PC pc;
	ierr = KSPGetPC(system_ksp,&pc); CHKERRQ(ierr);

	// some notes:
	// H^T has the pressure coupling
	// H has the flux coupling

	// setup the ksp for the stokes system
	// we have the velocity_matrix that has the stokes stuff created and restricted to the NS dofs
	// okay let us just only let it work for stabilised

	// ********* SETUP SCHUR_STOKES KSP **************** //

	// create
	KSP schur_stokes_ksp;
	ierr = KSPCreate(PETSC_COMM_WORLD,&schur_stokes_ksp);// CHKERRQ(ierr);

	PC schur_stokes_pc;
	ierr = KSPGetPC(schur_stokes_ksp,&schur_stokes_pc);// CHKERRQ(ierr);

	// set the operators for the velocity matrix inversion
	ierr = KSPSetOperators(schur_stokes_ksp,velocity_matrix->mat(),velocity_matrix->mat());// CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(schur_stokes_ksp,"schur_stokes_");// CHKERRQ(ierr);
	PetscBool initial_guess = PETSC_FALSE;
	ierr = KSPSetInitialGuessNonzero(schur_stokes_ksp,initial_guess);
	// setup up the operators for the first inner solve
	ierr = KSPSetFromOptions (schur_stokes_ksp);// CHKERRQ(ierr);

	setup_is_simple (sys,pc);

	std::string split_1 = "0";
	std::string split_2 = "1";
	PCFieldSplitSetIS(schur_stokes_pc,split_1.c_str(),velocity_is);
	PCFieldSplitSetIS(schur_stokes_pc,split_2.c_str(),pressure_is);

	ierr = KSPSetUp(schur_stokes_ksp);// CHKERRQ(ierr);	





	// Scaled Pressure Mass Matrix
	if(es->parameters.get<unsigned int>("preconditioner_type_schur_stokes") == 2)	// shell pc..
	{

		if(!es->parameters.get<bool>("stab"))
		{
			std::cout << "Quadratic elements not currently supported for fieldsplit schur stokes presolve." << std::endl;
			std::cout << "EXITING" << std::endl;
			std::exit(0);
		}

		PC schur_stokes_pc;
		ierr = KSPGetPC(schur_stokes_ksp,&schur_stokes_pc);// CHKERRQ(ierr);

		KSP* subksp;	//subksps
		ierr = PCFieldSplitGetSubKSP(schur_stokes_pc,NULL,&subksp);


		ierr = KSPSetOptionsPrefix(subksp[0],"convection_diffusion_");// CHKERRQ(ierr);
		ierr = KSPSetFromOptions (subksp[0]);// CHKERRQ(ierr);


		std::cout << "Setting up LSC James scaled stabilised preconditioner." << std::endl;

		// set ksp to preonly because we are defining the inverse
		ierr = KSPSetType(subksp[1], KSPPREONLY);// CHKERRQ(ierr);

		PC schur_pc;
		ierr = KSPGetPC(subksp[1],&schur_pc);// CHKERRQ(ierr);

		ierr = PCSetType(schur_pc,PCSHELL);// CHKERRQ(ierr);

		ierr = ShellPCCreate(&schur_stokes_shell);// CHKERRQ(ierr);

		ierr = PCShellSetApply(schur_pc,LSCScaledStabilisedShellPCApply);// CHKERRQ(ierr);
		ierr = PCShellSetContext(schur_pc,schur_stokes_shell);// CHKERRQ(ierr);

		ierr = PCShellSetDestroy(schur_pc,NSShellDestroy);// CHKERRQ(ierr);

		ierr = PCShellSetName(schur_pc,"LSC James Scaled Stabilised Preconditioner");// CHKERRQ(ierr);

		       
		ierr = LSCScaledStabilisedShellPCSetUp(schur_pc,velocity_mass_matrix->mat(),subksp[1],es->parameters.get<bool> ("negative_bfbt_schur_complement"));// CHKERRQ(ierr);








		shell_pc_created = true;


	}

	ierr = KSPView(schur_stokes_ksp,PETSC_VIEWER_STDOUT_WORLD);















	// construct the n_bc pressure coupling vectors (column vectors from H^T):
	// - it's relatively easy to retrieve the column Vecs from a Mat
	// - we won't need them again so we can create a std::vector of them
	// - we only want to do the inverse of the stokes so we need to truncate the Vecs

	// get the Bt (H^T) so that we have vecs of the correct size
	Mat S,Bt,B;
	Mat A;
	ierr = KSPGetOperators(schur_ksp,&S,NULL);// CHKERRQ(ierr);
	ierr = MatSchurComplementGetSubMatrices(S,NULL,&A,&Bt,&B,NULL);



	//std::cout << "schur stokes matrix:\n\n\n\n\n\n\n\n" << std::endl;
	//MatView(velocity_matrix->mat(),PETSC_VIEWER_STDOUT_WORLD);







	//std::cout << "stokes matrix:\n\n\n\n\n\n\n\n" << std::endl;
	//MatView(A,PETSC_VIEWER_STDOUT_WORLD);



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
	for(unsigned int i=0; i<(unsigned int)num_coupling_points; i++)
	{
		PetscScalar col_number[1];
	
		PetscInt coupling_point[1] = {(PetscInt)i};	// the coupling point we are interested in (must be an array...)
		ierr = VecGetValues(non_zero_cols,1,coupling_point,col_number);

		// extract the correct column as a Vec
		ierr = MatGetColumnVector(Bt,Bt_col,col_number[0]);

		// copmute F^-1BT for this particular row
		std::cout << "Schur Stokes Solve " << i << std::endl;
		ierr = KSPSolve(schur_stokes_ksp,Bt_col,T_col);

		//Mat schur_stokes_mat;
		//Mat schur_stokes_prec;
		//ierr = KSPGetOperators(schur_stokes_ksp,&schur_stokes_mat,&schur_stokes_prec);// CHKERRQ(ierr);

		//std::cout << "schur stokes mat:\n\n\n\n\n\n\n\n" << std::endl;
		//MatView(schur_stokes_mat,PETSC_VIEWER_STDOUT_WORLD);

		//std::cout << "schur stokes prec:\n\n\n\n\n\n\n\n" << std::endl;
		//MatView(schur_stokes_prec,PETSC_VIEWER_STDOUT_WORLD);

		//VecView(T_col,PETSC_VIEWER_STDOUT_SELF);

	
		double norm;
		VecNorm(T_col,NORM_2,&norm);
		std::cout << "Stokes result norm = " << norm << std::endl;

		// loop over the rows for the other bit of the matrix multiplication		
		for(unsigned int j=0; j<(unsigned int)num_coupling_points; j++)
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
	ierr = MatDuplicate(new_T,MAT_COPY_VALUES,&schur_complement_approx);



	// solve the equation A^-1 * f for each n_bc pressure coupling vector

	// dot product each of the results with each of the n_bc flux coupling vectors

	// put the results into a schur_complement_approx Mat in the correct positions.

	ierr = VecDestroy(&Bt_col);
	ierr = VecDestroy(&T_col);
	ierr = VecDestroy(&B_row);
	ierr = MatDestroy(&new_T);
	ierr = KSPDestroy(&schur_stokes_ksp);


	std::cout << "Done constructing Schur Stokes extras..." << std::endl;
	//std::exit(0);

	return 0;
}




// ********************************************************************** //

















// ****************** CONSTRUCT THE SCHUR STOKES MATRIX ************** //

// the pressure coupling dofs and flux coupling dofs corresponding to the same 
int NavierStokesCoupled::test_post_solve (TransientLinearImplicitSystem * sys)
{

	std::cout << "Constructing Post solve crap..." << std::endl;
	std::cout << "Could take a while..." << std::endl;
	
  	PetscErrorCode ierr;	// unused

	//get petsc solver
	PetscLinearSolver<Number>* system_linear_solver =
			libmesh_cast_ptr<PetscLinearSolver<Number>* >
			(sys->linear_solver.get());
	KSP system_ksp = system_linear_solver->ksp();
	PC pc;
	ierr = KSPGetPC(system_ksp,&pc); CHKERRQ(ierr);



	std::cout << "Setting up Velocity Monolithic preconditioner." << std::endl;


	std::cout << "1" << std::endl;
	KSP* subksp;	//subksps
	ierr = PCFieldSplitGetSubKSP(pc,NULL,&subksp);// CHKERRQ(ierr);
	std::cout << "1" << std::endl;
	KSP schur_ksp = subksp[1];
	PC ns_pc;
	KSPGetPC(subksp[0],&ns_pc);

	// get the IS for the NS dofs and the 0D dofs

	IS ns_is,rd_is;
	std::string split_1 = "0";
	std::string split_2 = "1";
	PCFieldSplitGetIS(pc,split_1.c_str(),&ns_is);
	PCFieldSplitGetIS(pc,split_2.c_str(),&rd_is);

	// some notes:
	// H^T has the pressure coupling
	// H has the flux coupling

	// setup the ksp for the stokes system
	// we have the velocity_matrix that has the stokes stuff created and restricted to the NS dofs
	// okay let us just only let it work for stabilised

	// ********* SETUP SCHUR_STOKES KSP **************** //

	std::cout << "1" << std::endl;
	// get the A so that we can use this as the matrix operator
	// get the Bt so that we can do the multiplication with the 0D vector later
	Mat S,A,Bt;
	ierr = KSPGetOperators(schur_ksp,&S,NULL);// CHKERRQ(ierr);
	ierr = MatSchurComplementGetSubMatrices(S,NULL,&A,&Bt,NULL,NULL);


	std::cout << "1" << std::endl;
	// create
	KSP post_solve_ksp;
	ierr = KSPCreate(PETSC_COMM_WORLD,&post_solve_ksp);// CHKERRQ(ierr);
	PC post_solve_pc;
	KSPGetPC(post_solve_ksp,&post_solve_pc);
	
	// set the operators for the velocity matrix inversion
	//ierr = KSPSetOperators(post_solve_ksp,velocity_matrix->mat(),velocity_matrix->mat());// CHKERRQ(ierr);

	ierr = KSPSetOperators(post_solve_ksp,A,A);// CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero (post_solve_ksp, PETSC_TRUE);// CHKERRQ(ierr);

	std::cout << "2" << std::endl;
	//ierr = KSPSetOptionsPrefix(post_solve_ksp,"post_solve_");// CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(post_solve_ksp,"post_solve_");// CHKERRQ(ierr);
	// setup up the operators for the first inner solve
	ierr = KSPSetFromOptions (post_solve_ksp);// CHKERRQ(ierr);
	std::cout << "2" << std::endl;

	setup_is_simple (sys,pc);

	PCFieldSplitSetIS(post_solve_pc,split_1.c_str(),velocity_is);
	PCFieldSplitSetIS(post_solve_pc,split_2.c_str(),pressure_is);

	//ierr = KSPView(post_solve_ksp,PETSC_VIEWER_STDOUT_WORLD);
	std::cout << "yeah?" << std::endl;
	ierr = KSPSetUp(post_solve_ksp);// CHKERRQ(ierr);
	std::cout << "yeah?" << std::endl;
	//ierr = KSPView(post_solve_ksp,PETSC_VIEWER_STDOUT_WORLD);


	std::cout << "1" << std::endl;
	// LSC James scaled stabilised elements
	if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 9
		&& !es->parameters.get<bool>("direct_post_solve"))
	{
	
		std::cout << "boo" << std::endl;
		PC post_solve_pc;
		std::cout << "boo" << std::endl;
		ierr = KSPGetPC(post_solve_ksp,&post_solve_pc);// CHKERRQ(ierr);

		std::cout << "boo" << std::endl;
		KSP* subksp;	//subksps
		std::cout << "boo" << std::endl;
		ierr = PCFieldSplitGetSubKSP(post_solve_pc,NULL,&subksp);

		std::cout << "Setting up LSC James scaled stabilised preconditioner for post solve." << std::endl;

		// set ksp to preonly because we are defining the schur inverse
		ierr = KSPSetType(subksp[1], KSPPREONLY);// CHKERRQ(ierr);

		PC schur_pc;
		ierr = KSPGetPC(subksp[1],&schur_pc);// CHKERRQ(ierr);

		ierr = PCSetType(schur_pc,PCSHELL);// CHKERRQ(ierr);

		ierr = ShellPCCreate(&post_solve_shell);// CHKERRQ(ierr);

		ierr = PCShellSetApply(schur_pc,LSCScaledStabilisedShellPCApply);// CHKERRQ(ierr);
		ierr = PCShellSetContext(schur_pc,post_solve_shell);// CHKERRQ(ierr);

		ierr = PCShellSetDestroy(schur_pc,NSShellDestroy);// CHKERRQ(ierr);

		ierr = PCShellSetName(schur_pc,"LSC James Scaled Stabilised Preconditioner");// CHKERRQ(ierr);

		       
		ierr = LSCScaledStabilisedShellPCSetUp(schur_pc,velocity_mass_matrix->mat(),subksp[1],es->parameters.get<bool> ("negative_bfbt_schur_complement"));// CHKERRQ(ierr);

		shell_pc_created = true;
	}
	else
	{
		//std::cout << "any preconditioner other than 9, not supported. EXITING..." << std::endl;
		//std::exit(0);
	}



	std::cout << "yeah?" << std::endl;
	//ierr = KSPSetUp(post_solve_ksp);// CHKERRQ(ierr);	
	std::cout << "yeah?" << std::endl;
	//ierr = KSPView(post_solve_ksp,PETSC_VIEWER_STDOUT_WORLD);




	PetscInt m;
	PetscInt m_local;
	PetscInt n;
	PetscInt n_local;
	ierr = MatGetSize(A,&m,&n);
	ierr = MatGetLocalSize(A,&m_local,&n_local);

	Vec result;
	VecCreate(PETSC_COMM_WORLD,&result);
	VecSetSizes(result,m_local,m);
	VecSetFromOptions(result);

	Vec initial_guess_vec;
	//AutoPtr<PetscVector<double> > petsc_copy;
	PetscVector<double> & petsc_copy = cast_ref<PetscVector<double>&>(const_cast<NumericVector<double>&>(*sys->solution));
	VecGetSubVector(petsc_copy.vec(),ns_is,&initial_guess_vec);

	Vec rhs_vec;
	PetscVector<double> & petsc_copy_rhs = cast_ref<PetscVector<double>&>(const_cast<NumericVector<double>&>(*sys->rhs));
	VecGetSubVector(petsc_copy_rhs.vec(),ns_is,&rhs_vec);

	Vec temp_rd_vec;
	VecGetSubVector(petsc_copy.vec(),rd_is,&temp_rd_vec);
	VecScale(temp_rd_vec,-1.0);

	MatMultAdd(Bt,temp_rd_vec,rhs_vec,rhs_vec);

	double norm;
	VecNorm(initial_guess_vec,NORM_2,&norm);
	std::cout << "before norm = " << norm << std::endl;
	std::cout << "before solve..." << std::endl;
	ierr = KSPSolve(post_solve_ksp,rhs_vec,initial_guess_vec);
	std::cout << "after solve..." << std::endl;
	VecNorm(initial_guess_vec,NORM_2,&norm);
	std::cout << "after norm = " << norm << std::endl;






	ierr = VecDestroy(&result);
	ierr = KSPDestroy(&post_solve_ksp);


	//std::cout << "exiting" << std::endl;
	//std::exit(0);

	return 0;
}
// ******************************************************************************** //



























// ****************** SETUP THE IS FOR MONOLITHIC PROBLEMS ************** //

// field_number is the IS that we want to return
int NavierStokesCoupled::setup_is_simple (TransientLinearImplicitSystem * sys,
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
		PetscErrorCode ierr=0;	// unused


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
			ierr = PCFieldSplitGetSubKSP(my_pc,&num_splits,&subksp); CHKERRQ(ierr);

			// only need to check the first split
			PC sub_pc;
			ierr = KSPGetPC(subksp[0],&sub_pc);// CHKERRQ(ierr);
			PCType sub_pc_type;
			ierr = PCGetType(sub_pc,&sub_pc_type);
			std::cout << "sub_pc type = " << sub_pc_type << std::endl;
			std::cout << "hmmm" << std::endl;

			// need to extract the pc corresponding to this 
			// we always look at the first fieldsplit cause we are interested in the navier stokes dofs
			sys_prefix = sys_prefix + static_cast<std::ostringstream*>( &(std::ostringstream() << 0) )->str() + "_";

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
					std::cout << "group_command = " << group_command << std::endl;

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
	    ISCreateLibMesh(sys->comm().get(), indices.size(),
		                     idx, PETSC_COPY_VALUES, &velocity_is);
	  }
	  else
	  {
	    std::cout << "num dofs " << i << " = " << indices.size() << std::endl;
	    ISCreateLibMesh(sys->comm().get(), indices.size(),
		                     idx, PETSC_COPY_VALUES, &pressure_is);

	  }
  }
 
	return 0;

}

// ********************************************************************************** //













// ********************** COMPUTE AND OUTPUT EIGENVALUES **************************** //
// compute and output the eigenvalues of the preconditioned operator to file
// can only handle a small number of unknowns like 500.
// if we are more than 500 then
int NavierStokesCoupled::compute_and_output_eigenvalues(TransientLinearImplicitSystem * system)
{

	PetscErrorCode ierr;	// unused

	//get petsc solver
	PetscLinearSolver<Number>* system_linear_solver =
			libmesh_cast_ptr<PetscLinearSolver<Number>* >
			(system->linear_solver.get());
	KSP system_ksp = system_linear_solver->ksp();

	// check the number of dofs
	dof_id_type num_dofs = system->n_dofs();

	if(num_dofs < 2000)
	{
		std::cout << "Calculating all eigenvalues for preconditioned system." << std::endl;
		PetscInt nmax = num_dofs;
		PetscReal r[num_dofs];
		PetscReal c[num_dofs];		
		
		std::cout << "Computing eigenvalues..." << std::endl;
		ierr = KSPComputeEigenvaluesExplicitly(system_ksp,nmax,r,c); CHKERRQ(ierr);

		std::cout << "Printing eigenvalues:" << std::endl;
		for(unsigned int i=0; i<num_dofs; i++)
		{
			std::cout << " eig " << i << " = " << r[i] << " + " << c[i] << "i" << std::endl;
		}

		if(es->parameters.get<bool> ("write_eigenvalues"))
		{
			eigenvalues_file << t_step << "\t" << nonlinear_iteration;

			for(unsigned int i=0; i<num_dofs; i++)
			{
				eigenvalues_file << "\t" << r[i] << "\t" << c[i];
			}

			eigenvalues_file << std::endl;
		}
	}
	else
	{
		std::cout << "Calculating eigenvalues for large system not implemented yet." << std::endl;
		//std::cout << "Calculating extremum eigenvalues for preconditioned system." << std::endl;
		//ierr = KSPComputeEigenvalues(KSP ksp,PetscInt n,PetscReal r[],PetscReal c[],PetscInt *neig)
	}

	
	return 0;
}

// ********************************************************************************** //















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






// *************************SET ELEM PID 3D************************** //


void NavierStokesCoupled::set_elem_proc_id_3d()
{

	//add the radius variable
	// i'm sure we only need to do this once!
	//create an equation system to hold data like the radius of the element
	std::cout << "Setting 3D proc ids." << std::endl;
	
	const DofMap& dof_map = extra_3d_data_system->get_dof_map();

	MeshBase::const_element_iterator       el     =
		mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el =
		mesh.active_local_elements_end();

	std::vector<dof_id_type> dof_indices;
	const unsigned int proc_id_var = extra_3d_data_system->variable_number ("proc_id");

	for ( ; el != end_el; ++el)
	{
		const Elem* elem = *el;
		if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
		{
			dof_map.dof_indices (elem, dof_indices,proc_id_var);

			for(unsigned int i=0; i < dof_indices.size(); i++)
				extra_3d_data_system->solution->set(dof_indices[i], elem->processor_id());
		}
	}

	//must close the vector after editing it

	extra_3d_data_system->solution->close();
	std::cout << "3D proc ids set." << std::endl;
}











// ********************** SCALE THE SOLUTION APPROPRIATELY *************************** //

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




