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

	std::cout << "WHAT THE FUCK" << std::endl;
	
	// need to scale the mesh by the length scale and/or the incoming scaling
	if(!restart)
	{
		// if we are doing a 2D simulation we need to change mesh back to having dimension 2
		if(!_es->parameters.get<bool> ("threed"))
			_mesh.set_mesh_dimension(2);

	std::cout << " _es->parameters.get<bool> (bool) = " << _es->parameters.get<bool> ("threed") << std::endl;
	std::cout << "_mesh.mesh_dimension() = " << _mesh.mesh_dimension() << std::endl;

		// both the cylindrical pipe and bifurcating pipe require input mesh as well as bifurcating pipe
		if(_es->parameters.get<unsigned int> ("geometry_type") == 0 || _es->parameters.get<unsigned int> ("geometry_type") == 1 
				|| _es->parameters.get<unsigned int> ("geometry_type") == 4  || _es->parameters.get<bool> ("expanding_pipe"))
		{
			_mesh.read(_es->parameters.get<std::string> ("mesh_file"), &mesh_data);

			//scale the mesh to make SI units

			

			if(!_es->parameters.get<bool>("linear_shape_functions"))
					_mesh.all_second_order();
			mesh_refinement.uniformly_refine (_es->parameters.get<unsigned int> ("no_refinement"));

			double mesh_input_scaling_3d = _es->parameters.get<double>("mesh_input_scaling_3d");
			MeshTools::Modification::scale(_mesh,mesh_input_scaling_3d,mesh_input_scaling_3d,mesh_input_scaling_3d);
		}
		// the cubioids can be auto generated
		else if(es->parameters.get<unsigned int> ("geometry_type") == 2 || es->parameters.get<unsigned int> ("geometry_type") == 3
						|| _es->parameters.get<unsigned int> ("geometry_type") == 5)
		{
			std::cout << "hmmm" << std::endl;

			if(threed){
				if(_es->parameters.get<bool>("quads"))
					MeshTools::Generation::build_cube (_mesh,2,2,2,0., 1.,0., 1.,0., 1.0,HEX8);
				else
					MeshTools::Generation::build_cube (_mesh,5,5,5,0., 1.,0., 1.,0., 1.0,TET4);					
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
			std::cout << "yeah" << std::endl;

			if(!_es->parameters.get<bool>("linear_shape_functions"))
				_mesh.all_second_order();

			mesh_refinement.uniformly_refine (_es->parameters.get<unsigned int> ("no_refinement"));

		}

		//scale the mesh to make length scale
	
		double length_scale = _es->parameters.get<double>("length_scale");
		MeshTools::Modification::scale(_mesh,1./length_scale,1./length_scale,1./length_scale);
		
	}
	else		// no scaling necessary
	{
		std::ostringstream file_name;
		std::ostringstream file_name_es;
		std::ostringstream file_name_mesh;

		file_name << output_folder.str() << "out_3D";
 
		if(particle_deposition == 1)
		{
			std::cout << "particle_deposition_start_time_step = " << 	_es->parameters.get<unsigned int> ("particle_deposition_start_time_step") << std::endl;

			if(_es->parameters.get<bool> ("unsteady_from_steady"))
			{
				  		// We write the file in the ExodusII format.
				file_name_es << file_name.str() << "_es_" << std::setw(4) 
					<< std::setfill('0') << 1 << ".xda";

				std::cout << "new_filename = " << file_name_es.str() << std::endl;

				file_name_mesh << file_name.str() << "_mesh_" << std::setw(4) 
					<< std::setfill('0') << 1 << ".xda";

			}
			else
			{
				  		// We write the file in the ExodusII format.
				file_name_es << file_name.str() << "_es_" << std::setw(4) 
					<< std::setfill('0') << _es->parameters.get<unsigned int> ("particle_deposition_start_time_step") << ".xda";

				std::cout << "new_filename = " << file_name_es.str() << std::endl;

				file_name_mesh << file_name.str() << "_mesh_" << std::setw(4) 
					<< std::setfill('0') << _es->parameters.get<unsigned int> ("particle_deposition_start_time_step") << ".xda";
			}

		}
		else	//normal restart
		{
			std::cout << "restart_time_step = " << 	_es->parameters.get<unsigned int> ("restart_time_step") << std::endl;

		    		// We write the file in the ExodusII format.
			file_name_es << file_name.str() << "_es_" << std::setw(4) 
				<< std::setfill('0') << _es->parameters.get<unsigned int> ("restart_time_step") << ".xda";

			std::cout << "new_filename = " << file_name_es.str() << std::endl;

			file_name_mesh << file_name.str() << "_mesh_" << std::setw(4) 
				<< std::setfill('0') << _es->parameters.get<unsigned int> ("restart_time_step") << ".xda";
		}

		_mesh.read(file_name_mesh.str());
		_es->read(file_name_es.str(), libMeshEnums::READ);

		// after reading we need to rescale, but only after dof_variable_type_3d has been set
		
		
	}

	std::cout << "WHAT THE FUCK" << std::endl;


	std::cout << " _es->parameters.get<bool> (bool) = " << _es->parameters.get<bool> ("threed") << std::endl;
	std::cout << "_mesh.mesh_dimension() = " << _mesh.mesh_dimension() << std::endl;
	//********* RENAME THE BOUNDARIES AND SUBDOMAINS FOR ALL DIFF PROBLEMS *****//
  //
	// we always read in the 3d mesh first so we can set all subdomain ids to 0
	// how is everything labelled
	// 3d inflow has boundary id 0
	// interfaces have boundary ids 1-n
	// 3d subdomain has id 0
	// 1d subdomains have ids 1-n corresponding to interfaces
	//

	
	// the pipe geometry or the expanding_pipe geometry
	if(_es->parameters.get<unsigned int> ("geometry_type") == 0 || _es->parameters.get<bool> ("expanding_pipe"))
	{
		if(!restart)
		{
			MeshTools::Modification::change_boundary_id(_mesh,1010,-1);		//walls

			MeshTools::Modification::change_boundary_id(_mesh,1011,0);			//inflow
			MeshTools::Modification::change_boundary_id(_mesh,1012,1);			//outflows

			MeshTools::Modification::change_subdomain_id(_mesh,1020,0);		//volume

			if(es->parameters.get<bool> ("old_geometry"))
			{
				MeshTools::Modification::change_boundary_id(_mesh,0,500);			//inflow
				MeshTools::Modification::change_boundary_id(_mesh,1,0);			//outflows
				MeshTools::Modification::change_boundary_id(_mesh,500,1);			//outflows
			}
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
	}
	else if(es->parameters.get<unsigned int> ("geometry_type") == 1 || es->parameters.get<unsigned int> ("geometry_type") == 4)	// bifurcating pipe
	{
		if(!restart)
		{

			std::set<boundary_id_type> bdyids = _mesh.boundary_info->get_boundary_ids();
			std::set<subdomain_id_type> subids;
			_mesh.subdomain_ids(subids);

			// we wanna do this relabelling automagically, and yes sets are ordered in a particular way
			if(true)//_es->parameters.get<std::string> ("mesh_file").compare("meshes/coarse_mesh_truncated_labelled.msh") == 0)
			{
				int count = -1;
				for (std::set<boundary_id_type>::iterator it=bdyids.begin(); it!=bdyids.end(); ++it)
				{
					MeshTools::Modification::change_boundary_id(_mesh,*it,count);		//walls
					count++;
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

				if(es->parameters.get<bool> ("old_geometry"))
				{
					MeshTools::Modification::change_boundary_id(_mesh,0,500);			//inflow
					MeshTools::Modification::change_boundary_id(_mesh,1,0);			//outflows
					MeshTools::Modification::change_boundary_id(_mesh,500,1);			//outflows
				}
			}
			else
			{
				MeshTools::Modification::change_boundary_id(_mesh,1010,-1);		//walls

				MeshTools::Modification::change_boundary_id(_mesh,1011,0);			//inflow
				MeshTools::Modification::change_boundary_id(_mesh,1012,1);			//outflows

				MeshTools::Modification::change_subdomain_id(_mesh,1020,0);		//volume

				if(es->parameters.get<bool> ("old_geometry"))
				{
					MeshTools::Modification::change_boundary_id(_mesh,0,500);			//inflow
					MeshTools::Modification::change_boundary_id(_mesh,1,0);			//outflows
					MeshTools::Modification::change_boundary_id(_mesh,500,1);			//outflows
				}
			}
		}


		_mesh.print_info();

		// so this has the boundary ids of the inflows and outflows, not the walls
		boundary_ids.resize(_mesh.boundary_info->n_boundary_ids() - 1);

		_es->parameters.set<unsigned int> ("num_1d_trees") = boundary_ids.size() - 1;


		
		for(unsigned int i=0; i < boundary_ids.size(); i++)
			boundary_ids[i] = i;

		//set surfaceboundary stuff
		//inflow_surface_boundary_object.init(_mesh,0);

		for(unsigned int i=0; i < boundary_ids.size(); i++)
		{
			SurfaceBoundary* surface_boundary = new SurfaceBoundary(*_es);
			surface_boundaries.push_back(surface_boundary);
			surface_boundaries[i]->init(_mesh,i,0);

			std::cout << "surface " << i << ": normal is " << surface_boundaries[i]->get_normal() << std::endl;
			std::cout << "  and centroid is " << surface_boundaries[i]->get_centroid() << std::endl << std::endl;;
		}

		std::cout << "voodoo" << std::endl;

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
		std::cout << "hey" << std::endl;
	}
	// the cuboid geometry
	else if(_es->parameters.get<unsigned int> ("geometry_type") == 2)
	{

		if(!restart)
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
			
		}

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
	// the axisymmetric cuboid geometry
	else if(_es->parameters.get<unsigned int> ("geometry_type") == 5)
	{

		if(!restart)
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
			
		}

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

		if(!restart)
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


		}

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

	std::cout << "WHAT THE FUCK" << std::endl;
	//populate subdomains_3d
	subdomains_3d.push_back(0);


}

// ************************************************************************** //








// ************************** SETUP 3D SYSTEM ******************************* //

void NavierStokesCoupled::setup_3d_system(TransientLinearImplicitSystem* system)
{

	//some parameters we need
	bool linear_shape_functions = es->parameters.get<bool>("linear_shape_functions");
	unsigned int  problem_type = es->parameters.get<unsigned int> ("problem_type");
	//Real  restart_time = es->parameters.get<Real> ("restart_time");
	//const MeshBase& mesh = es->get_mesh();

	std::set<subdomain_id_type> active_subdomains;
	active_subdomains.clear(); active_subdomains.insert(0);

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
	if(es->parameters.get<unsigned int>("preconditioner_type"))
	{
		system->add_matrix("Preconditioner");
	}

	//******* NOW WE NEED TO TAKE CARE OF THE BOUNDARY CONDITION STUFF *****//

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

	std::cout << "fooey" << std::endl;

		if(!es->parameters.get<bool>("pearson_example"))
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

	std::cout << "yeah we doin some xact solutions" << std::endl;


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


	if(restart)
	{
		system->update();
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

 	std::cout << "  Refining the mesh..." << std::endl;

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
			flux_values_3d[i] = picard->calculate_flux(i);
			pressure_values_3d[i] = picard->calculate_pressure(i);
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

	std::cout << "writing 3D solution" << std::endl;

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
		  &es->get_system<TransientLinearImplicitSystem> ("Coupled-Navier-Stokes");
	}
	else
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("Navier-Stokes-3D");
	}

	// exodus file writer
	ExodusII_IO_Extended exo = ExodusII_IO_Extended(mesh);

	exo.set_var_scalings(var_scalings);

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

	std::cout <<"ja" << std::endl;
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

		for(unsigned int i=0; i<(unsigned int)(-refinement_level_difference); i++)
		{
			mesh_refinement.uniformly_coarsen(1);
			es->reinit();	//NEEDED: so that correct xda file is written, of course we lose the accuracy... should actually do a temp copy
		}
	}		

	std::cout <<"ja" << std::endl;
	//system->update();
	


	//the time step variable needs to reflect the position in the file_name
	// since we want one file per time step we just put 1 here
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
	write_elem_pid_3d(exo);




	// We write the equation systems and mesh object in .xda format
	// for restarting... perhaps we don't wanna do this every time step
	// but every n time_steps

	std::cout << "EXODUSII output for timestep " << t_step
		<< " written to " << file_name_soln.str() << std::endl;


	if(backup)
	{
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
		for(unsigned int i=0; i<(unsigned int)refinement_level_difference; i++)
		{
			mesh_refinement.uniformly_coarsen(1);
			es->reinit();	//NEEDED: so that correct xda file is written, of course we lose the accuracy... should actually do a temp copy
		}
	}
	else if(refinement_level_difference < 0 )
	{

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


	// if we are doing stokes or it is the first iteration of navier stokes (which must be stokes anyway!) then use stokes otherwise use navier stokes
	// set petsc options from file, petsc command line options will overwrite
	if((es->parameters.get<unsigned int>("t_step") == 1 && nonlinear_iteration == 1) || es->parameters.get<bool>("stokes"))
		PetscOptionsInsertString(es->parameters.get<std::string>("stokes_solver_options").c_str());
	else
	{
		PetscOptionsInsertString(es->parameters.get<std::string>("solver_options").c_str());
	}

	// ******************* SOME SETTINGS
	es->parameters.set<Real> ("linear solver tolerance") =
		 es->parameters.get<double>("outer_solver_rtol");//1.e-12;
	//set previous nonlinear residual to value it was
	previous_nonlinear_residual = current_nonlinear_residual;
	//es->reinit ();	// UPDATE TIME DEPENDENT DIRICHLET BOUNDARY CONDITIONS put after

	//petsc stuff
	//if we have different solvers for the different blocks then need a prefix
	PetscLinearSolver<Number>* system_linear_solver =
			libmesh_cast_ptr<PetscLinearSolver<Number>* >
			(system->linear_solver.get());

	//let us set up to use the same preconditioner
	system_linear_solver->reuse_preconditioner(es->parameters.get<bool>("reuse_preconditioner"));
	
	KSP system_ksp = system_linear_solver->ksp();

/*
	PC             pc;
	const PetscInt ufields[] = {0,1,2},pfields[] = {3};
	KSPGetPC(system_ksp,&pc);
	PCFieldSplitSetBlockSize(pc,4);
	PCFieldSplitSetFields(pc,"0",3,ufields,ufields);
	PCFieldSplitSetFields(pc,"1",1,pfields,pfields);
*/

	//set the outer tolerances here - maxits is now 20
	KSPSetTolerances(system_ksp,
			es->parameters.get<double>("outer_solver_rtol"),
			es->parameters.get<double>("outer_solver_atol"),
			1e38,es->parameters.get<int>("outer_solver_maxits"));

	// 1e3 is the dtol need to change

	//std::cout << "hmmm" << std::endl;
	std::string prefix_3d = "ns3d_";
	system_linear_solver->set_prefix(prefix_3d);

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

	//std::cout << "fuck this" << std::endl;

	//std::cout << "prefix = " << prefix_2 <<std::endl;


	
	// need to recollect the ksp because perhaps some shit changed in restrict_solve_to
	system_linear_solver =
			libmesh_cast_ptr<PetscLinearSolver<Number>* >
			(system->linear_solver.get());
	system_ksp = system_linear_solver->ksp();

	int num_outer_its = 0;	
	int num_inner_its = 0;
	// output some data about the iterative solve
	// this is hard to do when doing restrict_solve_to, so just leave it

	KSPGetIterationNumber(system_ksp,&num_outer_its);
	if(!es->parameters.get<bool>("direct"))
	{
		//get the pc
		PC pc;
		KSPGetPC(system_ksp,&pc);

		//get the inner iteration ksp
		int num_splits = 0;
		KSP* subksp;	//subksps
		PCFieldSplitGetSubKSP(pc,&num_splits,&subksp);

		KSPGetIterationNumber(subksp[1],&num_inner_its);

		std::cout << "num outer its = " << num_outer_its << std::endl;
		std::cout << "num inner schur its = " << num_inner_its << std::endl;
	}

	int converged_reason = system_linear_solver->get_petsc_converged_reason();

	//std::cout << "KSPConvergedReason = " << converged_reason << std::endl;
	
	system->get_dof_map().enforce_constraints_exactly(*system);


  // Compute the difference between this solution and the last
  // nonlinear iterate.
	last_nonlinear_soln->add (-1., *system->solution);
	last_nonlinear_soln->close();

	std::vector<Real> weights;  // u, v
	weights.push_back(1.0);
	weights.push_back(1.0);
	if(threed)
		weights.push_back(1.0);
	weights.push_back(0.0);            // p

	if(es->parameters.get<bool>("optimisation_stabilised"))
	{
		weights.push_back(0.0);
		weights.push_back(0.0);
		if(threed)
			weights.push_back(0.0);
  	weights.push_back(0.0);            // p

	}

	std::vector<FEMNormType>	norms;	//if empty will automagically use discrete l2
	current_nonlinear_residual = system->calculate_norm(*last_nonlinear_soln,SystemNorm(norms, weights));
	// now normalise it
	current_nonlinear_residual /= system->calculate_norm(*system->solution,SystemNorm(norms, weights));


	//current_nonlinear_residual = last_nonlinear_soln->l2_norm();
	num_linear_iterations = system->n_linear_iterations();

	// Print out convergence
	std::cout << "Lin Solver converged ("
						<< converged_reason << ") at step: "
        		<< num_outer_its
        		<< ", Solver res: "
        		<< final_linear_residual
        		<< ",  Nonlinear res: "
        		<< current_nonlinear_residual
        		<< std::endl;


	total_linear_iterations += system->n_linear_iterations();

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
	

	// if iterative solver failed then adaptively refine
	// only for the steady case
	if(converged_reason < 0 && !es->parameters.get<bool>("direct"))
	{
		std::cout << " Nonlinear solver is not going to converge because "
              		<< "iterative solver is not converging."
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

	// Terminate the solution iteration if the difference between
	// this nonlinear iterate and the last is sufficiently small, AND
	// if the most recent linear system was solved to a sufficient tolerance.
  
	//doesn't make sense to compare with the first step in a time dependent problem
	if(current_nonlinear_residual > previous_nonlinear_residual && nonlinear_iteration > 2)	
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
	else if(num_outer_its == es->parameters.get<int>("outer_solver_maxits"))
	{
  			std::cout << " Nonlinear solver is not going to converge because "
              		<< "iterative solver is taking more than 20 iterations."
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


// ************************************************************************** //
// ********              SOLVE 3D LINEAR SYSTEM (replacement of libmesh)	
// ************

// NOTES:
// - no subset functionality
double NavierStokesCoupled::solve_and_assemble_3d_system(TransientLinearImplicitSystem * system)
{

		
	perf_log.push("assembly");
	if (system->assemble_before_solve)
	{
		// Assemble the linear system
		system->assemble ();
	}
  perf_log.pop("assembly");


	perf_log.push("solve");

	std::cout << "hi" << std::endl;
	//get petsc solver
	PetscLinearSolver<Number>* system_linear_solver =
			libmesh_cast_ptr<PetscLinearSolver<Number>* >
			(system->linear_solver.get());

	// If the linear solver hasn't been initialized, we do so here.
	//if (libMesh::on_command_line("--solver_system_names"))
	//	system->linear_solver->init((system->name()+"_").c_str());
	//else
	system->linear_solver->init();

	// Get the user-specifiied linear solver tolerance
	const Real tol = es->parameters.get<Real>("linear solver tolerance");

	// Get the user-specified maximum # of linear solver iterations
	const unsigned int maxits = es->parameters.get<unsigned int>("linear solver maximum iterations");

	// Solve the linear system.  Several cases:
	std::pair<unsigned int, Real> rval = std::make_pair(0,0.0);
	bool _shell_matrix = false;
	std::cout << "hi" << std::endl;
	// matrix always passed by reference, preconditioner can be passed by reference of pointer
	if(_shell_matrix)
	{
		// 1.) Shell matrix with or without user-supplied preconditioner.
		// rval = linear_solver->solve(*_shell_matrix, this->request_matrix("Preconditioner"), *solution, *rhs, tol, maxits);
	}
	else
	{
		// 2.) No shell matrix, with or without user-supplied preconditioner
		rval = system->linear_solver->solve (*system->request_matrix("System Matrix"), system->request_matrix("System Matrix"), *system->solution, *system->rhs, tol, maxits);
	}
	std::cout << "hi" << std::endl;
	perf_log.pop("solve");


	//_n_linear_iterations   = rval.first;
	double final_linear_residual = rval.second;

	// Update the system after the solve
	system->update();

	std::cout << "hi" << std::endl;
	return final_linear_residual;
}

// ************************************************************************** //






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

  AutoPtr<MeshBase> meshptr = mesh.clone();
  MeshBase &temp_mesh = *meshptr;
  temp_mesh.all_first_order();
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
		if(elem->subdomain_id() == 0)
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

  //ExodusII_IO io(mesh);
  //io.write(filename);
  //io.write_element_data(temp_es,time);
  io.write_element_data(temp_es);

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
		  &es->get_system<TransientLinearImplicitSystem> ("Coupled-Navier-Stokes");
	}
	else
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("Navier-Stokes-3D");
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




