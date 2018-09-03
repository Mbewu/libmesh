// ************* COUPLED NAVIER STOKES *********************** // 
// ******* WITH IMPLICIT EXPLICIT AND MONOLITHIC COUPLING **** //
//
// This file includes all the functions not specific to 
// 3D or 1D domains.
//
// main
// NavierStokesCoupled
// read_parameters
// set_double_parameter
// set_unsigned_int_parameter
// set_int_parameter
// set_bool_parameter
// set_string_parameterrea
// output_flux_and_pressure
// update_times
// update_time_scaling
// end_timestep_admin
// output_sim_data
// scaled_time
//
// *********************************************************** //

 
#include "coupled_navier_stokes.h"

// ************************* MAIN FUNCTION ************************ //
int main (int argc, char** argv)
{


	// Initialize libMesh.
	LibMeshInit init (argc, argv);

	//we want to get the input file from the input if possible	
	std::ostringstream input_file;
	std::ostringstream input_file_particle;

	// ******* READ IN COMMAND LINE AND CHECK IN CORRECT FORMAT **** //
	GetPot command_line (argc, argv);
	if(command_line.search("-input_file"))
	{
		input_file << command_line.next("navier.in");
		std::cout << "Using user specified input file: " << input_file.str() << std::endl;
	}
	else
	{
		input_file << "navier.in";
		std::cout << "Using default input file: " << input_file.str() << std::endl;
	}


	if(command_line.search("-particle_deposition_input_file"))
	{	
		input_file_particle << command_line.next("particle_deposition.in");
		std::cout << "Using user specified particle deposition input file: " << input_file_particle.str() << std::endl;
	}
	else
	{
		input_file_particle << "particle_deposition.in";	
		std::cout << "Using default particle deposition input file: " << input_file_particle.str() << std::endl;
	}

  //GetPot infile(input_file.str());
	NavierStokesCoupled navier_stokes_coupled(init,input_file.str(),input_file_particle.str(),command_line);
	

/*
	std::vector<double> time_steps;
	//time_steps.push_back(0.1);
	//time_steps.push_back(0.05);
	time_steps.push_back(0.01);
	//time_steps.push_back(0.005);
	time_steps.push_back(0.001);


	std::vector<double> viscosities;
	//viscosities.push_back(1.0);
	//viscosities.push_back(0.1);
	//viscosities.push_back(0.01);
	viscosities.push_back(0.001);

	input_file << "input_files/navier.in";
	for(unsigned int i=0; i<viscosities.size(); i++)
	{

		for(unsigned int j=0; j<time_steps.size();j++)
		{
			std::cout  << "TIME STEP = " <<time_steps[j] << std::endl;
			std::cout  << "viscosity = " << viscosities[i] << std::endl;
			std::cout << "EXPLICIT CALCULAION" << std::endl;
			NavierStokesCoupled navier_stokes_coupled4(init,input_file.str(),viscosities[i],time_steps[j],4);

			std::cout  << "TIME STEP = " <<time_steps[j] << std::endl;
			std::cout  << "viscosity = " << viscosities[i] << std::endl;
			std::cout << "IMPLICIT CALCULAION" << std::endl; 
			NavierStokesCoupled navier_stokes_coupled3(init,input_file.str(),viscosities[i],time_steps[j],3);

			std::cout  << "TIME STEP = " <<time_steps[j] << std::endl;
			std::cout  << "viscosity = " << viscosities[i] << std::endl;
			std::cout << "MONOLITHIC CALCULAION" << std::endl;
			NavierStokesCoupled navier_stokes_coupled5(init,input_file.str(),viscosities[i],time_steps[j],5);
			

		}
	}
*/
	

	return 0;
}









// ************************** CONSTRUCTOR ********************************** //
// set the input file in the construcor, this way we can even programmatically 
// run through the simulations in main.
NavierStokesCoupled::NavierStokesCoupled(LibMeshInit & init, std::string _input_file, std::string _input_file_particle, GetPot& command_line):
		mesh(init.comm(),3),
		mesh_3d(init.comm(),3),
		mesh_1d(init.comm(),3),
		mesh_1d_terminal(init.comm(),3),
		mesh_acinar(init.comm(),3),
		mesh_1d_acinar(init.comm(),3),
		mesh_centrelines(init.comm(),3),
		mesh_data(mesh),
		mesh_refinement (mesh),
		mesh_refinement_3d (mesh_3d),
		mesh_refinement_1d (mesh_1d),
		es(AutoPtr<EquationSystems>(new EquationSystems(mesh))),
		es_3d(AutoPtr<EquationSystems>(new EquationSystems(mesh_3d))),
		es_1d(AutoPtr<EquationSystems>(new EquationSystems(mesh_1d))),
		es_1d_terminal(AutoPtr<EquationSystems>(new EquationSystems(mesh_1d_terminal))),
		es_centrelines(AutoPtr<EquationSystems>(new EquationSystems(mesh_centrelines))),
		es_acinar(AutoPtr<EquationSystems>(new EquationSystems(mesh_acinar))),
		time(0.),
  	reduce_dt(0),
  	increase_dt(false),
  	refine_mesh(false),
		steps_since_last_dt_change(1),
		restart(false),
		perf_log("Main Program"),
		perf_log_move("move_particles"),
		sim_3d(false),
		sim_1d(false),
		input_file(_input_file),
		input_file_particle(_input_file_particle),
		comm_line(command_line),
		exit_program(false),
		threed(true),
		total_nonlinear_iterations(0),
		local_linear_iterations(0),
		total_linear_iterations(0),
		total_gmres_linear_iterations(0),
		total_cg_linear_iterations(0),
		total_max_iterations(0),
		total_max_gmres_iterations(0),
		total_max_cg_iterations(0),
		local_max_iterations(0),
		local_max_residual_iterations(0),
		particle_deposition(0),
		shell_pc_created(false),
		mono_shell_pc_created(false),
		first_3d_write(true),
		first_1d_write(true),
		init_names_done(false),
		ic_set(true),
		no_motion_end_particle_sim(false),
		total_fraction_added(0.),
		max_generations_3d(0),
		particle_fraction_0d(0),
		particle_fraction_3d(0),
		deposition_fraction_0d(0),
		deposition_fraction_3d(0),
		deposition_fraction_0d_max(0),
		deposition_fraction_3d_max(0),
		deposition_fraction_sed_0d(0),
		deposition_fraction_sed_0d_max(0),
		deposition_fraction_imp_0d(0),
		deposition_fraction_imp_0d_max(0),
		deposition_fraction_dif_0d(0),
		deposition_fraction_dif_0d_max(0),
		terminal_exit_fraction(0),
		terminal_exit_fraction_max(0),
		residual_linear_iteration(0),
		num_airways_in_full_mesh(0)
{

	perf_log.push("misc");
	perf_log.push("setup");
	PerfLog perf_log_setup("Setup");

	std::cout << "\n*** STARTING SIMULATION ***\n" << std::endl;

	//************* READ AND OUTPUT GENERAL PARAMETERS *************//
	infile = GetPot(_input_file);
	read_parameters();

	std::cout << "input file = " << _input_file << std::endl;
	std::cout << "input file particle = " << _input_file_particle << std::endl;

	//************* READ AND OUTPUT PARTICLE DEPOSITION PARAMETERS **********//
	if(es->parameters.get<unsigned int>("particle_deposition"))
	{
		infileparticle = GetPot(_input_file_particle);
		read_particle_parameters();
		output_particle_parameters();
	}

	output_parameters();
	output_command_line_options();
	print_parameters();








	std::cout << "hello" << std::endl;
	//Airway new_airway;
	//new_airway.set_generation(5);
	//std::cout << "generation of new airway = " << new_airway.get_generation() << std::endl;
	//std::exit(0);



/* some crap
	//set from input
	if(_viscosity != 0.0)
	{
		viscosity_vec[0] = _viscosity;
		es->parameters.set<Real> ("viscosity") = _viscosity;
	  es->parameters.set<Real> ("viscosity")   = viscosity_vec[0];
	}

	if(_dt != 0.0)
	{
		dt = _dt;
		es->parameters.set<Real>("dt") = dt;
	}

	// this should output to exodus every 0.1 timesteps so 80 per simulation should be okay
	//es->parameters.set<Real>("write_interval") = 0.1;

	if(_sim_type != 0.0)
	{
		sim_type = _sim_type;
	}
*/

/*
	output_folder << "results/param_sweep/";

	if(sim_type == 3)
		output_folder << "implicit/";
	if(sim_type == 4)
		output_folder << "explicit/";
	if(sim_type == 5)
		output_folder << "monolithic/";

	output_folder << "dt" << dt << "/viscosity" << es->parameters.get<Real> ("viscosity") << "/";
*/
	


	// ************ SETUP SYSTEM TYPES ************ //

	// 3d systems
	if(sim_3d && sim_type != 5)
	{
	  system_3d = &es->add_system<TransientLinearImplicitSystem>("ns3d");
	  system_neumann = &es->add_system<TransientLinearImplicitSystem>("Neumann-Variable");
		extra_3d_data_system = &es_3d->add_system<TransientExplicitSystem> ("extra_3d_data");
		system_3d_output = &es_3d->add_system<TransientExplicitSystem> ("ns3d_output");
		// setup 1d centrelines mesh
		if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))
			system_centrelines = &es_centrelines->add_system<TransientExplicitSystem> ("centrelines");
			
	}

	// 0d systems
	if(sim_1d && sim_type != 5)
	{
		system_1d = &es->add_system<TransientLinearImplicitSystem> ("ns1d");
		extra_1d_data_system = &es_1d->add_system<TransientExplicitSystem> ("extra_1d_data");
		system_1d_output = &es_1d->add_system<TransientExplicitSystem> ("ns1d_output");
		if(es->parameters.get<unsigned int>("acinar_model") == 1)
			system_acinar_output = &es_acinar->add_system<TransientExplicitSystem> ("acinar_output");

		if(particle_deposition == 2 || particle_deposition == 3 || particle_deposition == 4 || particle_deposition == 5 || particle_deposition == 6)
			particle_deposition_system_1d = &es_1d->add_system<TransientExplicitSystem> ("Particle-Deposition-1D");

	}

	// monolithic systems
	if(sim_type == 5)
	{
	
		system_coupled = &es->add_system<TransientLinearImplicitSystem> ("ns3d1d");
	  system_neumann = &es->add_system<TransientLinearImplicitSystem>("Neumann-Variable");
		extra_1d_data_system = &es_1d->add_system<TransientExplicitSystem> ("extra_1d_data");
		extra_3d_data_system = &es_3d->add_system<TransientExplicitSystem> ("extra_3d_data");
		system_3d_output = &es_3d->add_system<TransientExplicitSystem> ("ns3d_output");
		system_1d_output = &es_1d->add_system<TransientExplicitSystem> ("ns1d_output");
		if(es->parameters.get<unsigned int>("acinar_model") == 1)
			system_acinar_output = &es_acinar->add_system<TransientExplicitSystem> ("acinar_output");

		// setup 1d centrelines mesh
		if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))
			system_centrelines = &es_centrelines->add_system<TransientExplicitSystem> ("centrelines");

		// for calculating resistance etc
		system_1d = &es->add_system<TransientLinearImplicitSystem> ("ns1d");

	}
	// end setup system types









	// ************** SET UP MESHES **************** //

	mesh.partitioner()->set_custom_partitioning(es->parameters.set<unsigned int>("custom_partitioning"));
	mesh_3d.partitioner()->set_custom_partitioning(es->parameters.set<unsigned int>("custom_partitioning"));
	mesh_1d.partitioner()->set_custom_partitioning(es->parameters.set<unsigned int>("custom_partitioning"));
	mesh_1d_terminal.partitioner()->set_custom_partitioning(es->parameters.set<unsigned int>("custom_partitioning"));
	mesh_acinar.partitioner()->set_custom_partitioning(es->parameters.set<unsigned int>("custom_partitioning"));



	// ************** SET UP 3D MESH **************** //
	if(sim_3d)
	{
		perf_log_setup.push("setup_3d_mesh");
		setup_3d_mesh(&*es,mesh);

		std::cout << "Setting up output 3D mesh." << std::endl;
		setup_3d_mesh(&*es_3d,mesh_3d,true);

		perf_log_setup.pop("setup_3d_mesh");

		// in case not set already
		es->parameters.set<unsigned int>("n_initial_3d_elem") = mesh.n_elem();

		// setup 1d centrelines mesh
		if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))
			read_1d_mesh(true);


	}// setup 3d mesh





	// ******************** SETUP 1D MESH ********************* //
	if(sim_1d)
	{
		perf_log_setup.push("setup_1d_mesh");
		// the output mesh is handled within this function, 
		// so that we don't have to generate the 1d mesh and all its data twice
		setup_1d_mesh();
		// setup centrelines mesh
		perf_log_setup.pop("setup_1d_mesh");


	}// end setup 0d mesh

	// end setup meshes





	// initialise combined sparsity object
	bool augment_monolithic = false;
	bool augment_moghadam = false;
	bool augment_acinar = false;
	if(sim_type == 5 || sim_1d)
		augment_monolithic = true;
	if(es->parameters.get<unsigned int>("moghadam_coupling"))
		augment_moghadam = true;
	if(es->parameters.get<unsigned int>("acinar_model") == 1 && sim_1d)
		augment_acinar = true;


	// ************* SET UP 3D SYSTEM **************** //
	if(sim_3d)
	{
		perf_log_setup.push("setup_3d_system");
		if(!es->parameters.get<bool>("optimisation_stabilised"))
			picard = AutoPtr<Picard>(new Picard(*es,surface_boundaries,subdomains_3d,elem_to_airway,es->parameters.get<bool>("efficient_assembly")));
		else
			picard = AutoPtr<OptimisedStabilisedAssembler3D>(new OptimisedStabilisedAssembler3D(*es,surface_boundaries,subdomains_3d,elem_to_airway));

		if(sim_type == 5)
		{
			setup_3d_system(system_coupled);
			picard->set_pressure_coupled(true);
		}
		else
		{
			setup_3d_system(system_3d);
			system_3d->attach_assemble_object (*picard);
		}


		//get the extra couplings sorted for the moghadam matrix
		if(es->parameters.get<bool> ("moghadam_augment_sparsity_pattern"))
		{
			/*
			augment_sparsity_moghadam = AutoPtr<AugmentSparsityMoghadam>
										(new AugmentSparsityMoghadam(*es,subdomains_3d));
			system_3d->get_dof_map().attach_extra_sparsity_object(*augment_sparsity_moghadam);
			*/
		
			// we can attach it here because it's separate from acinar
			
			
			augment_sparsity_combined_3d = AutoPtr<AugmentSparsityCombined>(new AugmentSparsityCombined());
			augment_sparsity_combined_3d->create_moghadam_sparsity_object(*es,subdomains_3d);
			system_3d->get_dof_map().attach_extra_sparsity_object(*augment_sparsity_combined_3d);
			
		}

		perf_log_setup.pop("setup_3d_system");


		std::cout << "Setting up output 3D system." << std::endl;
		setup_3d_system(system_3d_output,true);


		// setup 1d centrelines system
		if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))
			setup_1d_system(system_centrelines,true,true);
			

	}// end setup 3d system










	// ******************* SETUP 1D SYSTEM ******************** //
	if(sim_1d)
	{
		
		perf_log_setup.push("setup_1d_system");
		if(es->parameters.get<unsigned int>("acinar_model") == 1)
		{
			ns_assembler = AutoPtr<NavierStokesAssembler>
									(new NavierStokesAssembler(*es,airway_data,subdomains_1d,elem_to_airway,subdomains_acinar,acinar_to_airway,elem_to_acinar));
		}
		else
		{
			ns_assembler = AutoPtr<NavierStokesAssembler>
									(new NavierStokesAssembler(*es,airway_data,subdomains_1d,elem_to_airway));

		}

		if(sim_type == 5)
		{
			if(es->parameters.get<unsigned int>("acinar_model") == 1)
			{
				ns_assembler_mono = AutoPtr<NavierStokesAssembler>
										(new NavierStokesAssembler(*es,airway_data,subdomains_1d,elem_to_airway,subdomains_acinar,acinar_to_airway,elem_to_acinar));
			}
			else
			{
				ns_assembler_mono = AutoPtr<NavierStokesAssembler>
										(new NavierStokesAssembler(*es,airway_data,subdomains_1d,elem_to_airway));

			}

			setup_1d_system(system_coupled);
			if(es->parameters.get<unsigned int>("acinar_model") == 1)
				setup_acinar_system(system_coupled);
			// don't give the 1d an additional prefix
			//system_1d->attach_assemble_object (*ns_assembler); //attach later

			std::cout << "cotard" << std::endl;
			// extra couplings for acinar bits
			augment_sparsity_combined_coupled = AutoPtr<AugmentSparsityCombined>(new AugmentSparsityCombined());
			if(es->parameters.get<unsigned int>("acinar_model") == 1)
			{
				/*
				std::cout << "hmm" << std::endl;
				augment_sparsity_acinar = AutoPtr<AugmentSparsityAcinar>
											(new AugmentSparsityAcinar(*es,airway_data,subdomains_1d,subdomains_acinar,elem_to_acinar,acinar_to_airway,elem_to_airway,true));
				std::cout << "hmm" << std::endl;
				system_coupled->get_dof_map().attach_extra_sparsity_object(*augment_sparsity_acinar);
				std::cout << "hmm" << std::endl;
				*/

				augment_sparsity_combined_coupled->create_acinar_sparsity_object(*es,airway_data,subdomains_1d,subdomains_acinar,elem_to_acinar,acinar_to_airway,elem_to_airway,true);

				// attach later when we setup the 3d part
				//system_coupled->get_dof_map().attach_extra_sparsity_object(*augment_sparsity_combined);
		
				
			}
			std::cout << "implared" << std::endl;
		}
		//else
		// always do this, so we have an exclusive 1d system to do calcs on
		{
			setup_1d_system(system_1d);
			
			if(es->parameters.get<unsigned int>("acinar_model") == 1)
				setup_acinar_system(system_1d);
			

			// setup output system

			system_1d->attach_assemble_object (*ns_assembler);

			//get the extra couplings sorted for the matrix, maybe don't need this?
			// yes, this is for the coupling of daughters etc
			/*
			augment_sparsity = AutoPtr<AugmentSparsityOnInterface>
										(new AugmentSparsityOnInterface(*es,airway_data,subdomains_3d,subdomains_1d,elem_to_airway));
			system_1d->get_dof_map().attach_extra_sparsity_object(*augment_sparsity);
			*/

			augment_sparsity_combined_1d = AutoPtr<AugmentSparsityCombined>(new AugmentSparsityCombined());
			augment_sparsity_combined_1d->create_monolithic_sparsity_object(*es,airway_data,subdomains_3d,subdomains_1d,elem_to_airway);

			// extra couplings for acinar bits
			if(es->parameters.get<unsigned int>("acinar_model") == 1)
			{
				/*
				augment_sparsity_acinar = AutoPtr<AugmentSparsityAcinar>
											(new AugmentSparsityAcinar(*es,airway_data,subdomains_1d,subdomains_acinar,elem_to_acinar,acinar_to_airway,elem_to_airway));
				system_1d->get_dof_map().attach_extra_sparsity_object(*augment_sparsity_acinar);
				*/

				augment_sparsity_combined_1d->create_acinar_sparsity_object(*es,airway_data,subdomains_1d,subdomains_acinar,elem_to_acinar,acinar_to_airway,elem_to_airway);
			}

			// attach after both created
			system_1d->get_dof_map().attach_extra_sparsity_object(*augment_sparsity_combined_1d);
			
		}


		setup_1d_system(system_1d_output,true);
		if(es->parameters.get<unsigned int>("acinar_model") == 1)
		{
			setup_acinar_system(system_acinar_output,true);		
		}

		if(particle_deposition == 3 ||  particle_deposition == 4 || particle_deposition == 5 ||  particle_deposition == 6)
		{		
			std::cout << "hullo" << std::endl;
			// i assume surface boundaries is just a 0 length vector by default...
			james_particle_deposition_object.setup(*es,*es_1d,airway_data,subdomains_1d,surface_boundaries,*es_1d_terminal,terminal_airway_data,terminal_airway_to_airway_map);
		}


		
		perf_log_setup.pop("setup_1d_system");
	}// end setup 0d system







	// ******************* SETUP MONOLITHIC SYSTEM ******************** //
	if(sim_type == 5)
	{
		ns_assembler_mono->set_coupled(true);
		coupled_assembler = AutoPtr<CoupledAssembler>
								(new CoupledAssembler(*picard,*ns_assembler_mono));
		system_coupled->attach_assemble_object (*coupled_assembler);

		//get the extra couplings sorted for the matrix
		/*
		augment_sparsity = AutoPtr<AugmentSparsityOnInterface>
									(new AugmentSparsityOnInterface(*es,airway_data,subdomains_3d,subdomains_1d,elem_to_airway,true));
		system_coupled->get_dof_map().attach_extra_sparsity_object(*augment_sparsity);
		*/
		augment_sparsity_combined_coupled->create_monolithic_sparsity_object(*es,airway_data,subdomains_3d,subdomains_1d,elem_to_airway,true);
		system_coupled->get_dof_map().attach_extra_sparsity_object(*augment_sparsity_combined_coupled);	// attach it now

		picard->set_subdomains_1d(subdomains_1d);
		picard->set_airway_data(airway_data);

	}




	// *********** RESET TIME STUFF ************** //
	time = es->parameters.get<Real>("restart_time");

	t_step = es->parameters.get<unsigned int> ("t_step");
	dt = es->parameters.get<Real> ("dt");

	std::cout << "dt = " << dt << std::endl;

	update_times();

	// ********






	


	// ******* INIT EQ SYSTEMS TO GET THEM READY FOR USE **************** //

	std::cout << "Init equation systems." << std::endl;
	es->init ();

	if(sim_3d)
	{
		std::cout << "initing 3d eq system for output" << std::endl;
		es_3d->init();
	}

	if(sim_1d)
	{
		std::cout << "initing 0d eq system for output" << std::endl;
		es_1d->init();

		if(es->parameters.get<unsigned int>("acinar_model") == 1)
		{
			std::cout << "initing acinar eq system for output" << std::endl;
			es_acinar->init();
		}
	}


	// setup 1d centrelines mesh
	if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))
	{
		std::cout << "initing centrelines eq system for output" << std::endl;
		es_centrelines->init();
	}

	std::cout << "done initing equation systems" << std::endl;










	// *************	INITIALISE SOME STUFF BEFORE THE LOOP ******* //

	if(sim_3d)
	{
		perf_log_setup.push("update_eq_systems");
		es->update();
		es_3d->update();
		perf_log_setup.pop("update_eq_systems");
		
		if(sim_type == 5)
		{
			last_nonlinear_soln = system_coupled->solution->clone();
			old_global_solution = system_coupled->solution->clone();
			*system_coupled->old_local_solution = *system_coupled->current_local_solution;
		}
		else
		{
			last_nonlinear_soln = system_3d->solution->clone();
			old_global_solution = system_3d->solution->clone();
			*system_3d->old_local_solution = *system_3d->current_local_solution;

			// for time dependent systems with a constant velocity inflow
			// we may want to set the old solution at t=0
			if(es->parameters.get<bool>("set_old_solution"))
			{
				// only update the constraints after the solution has been calculated
				system_3d->get_dof_map().enforce_constraints_exactly(*system_3d);	
				std::cout << "solution_norm at t=0 before solving = " << system_3d->solution->l2_norm() << std::endl;

				// for some reason this was not enough to set old_local_solution
				// probably missing an update
				*system_3d->old_local_solution = *system_3d->solution;
				std::cout << "old_solution_norm at t=0 before solving = " << system_3d->old_local_solution->l2_norm() << std::endl;
				
			}
		}
	}

	if(sim_1d)
	{
		perf_log_setup.push("update_eq_systems");
		es_1d->update();
		perf_log_setup.pop("update_eq_systems");

		// need this nonlinear vector if we are converging th 0D model
		if(true)//sim_type != 5)
		{
			old_global_solution_1d = system_1d->solution->clone();
			last_1d_nonlinear_soln = system_1d->solution->clone();
		}

		if(es->parameters.get<unsigned int>("acinar_model") == 1)
		{
			perf_log_setup.push("update_eq_systems");
			es_acinar->update();
			perf_log_setup.pop("update_eq_systems");
		}
	}


	// need to update for 0D restart (particle deposition) {not sure what this is supposed to do}
	/*
	if(false)//sim_1d && particle_deposition == 3)
	{
		perf_log_setup.push("update_eq_systems");
		es->update();
		perf_log_setup.pop("update_eq_systems");
	}
	*/
	// end initialise some stuff







	// **************************** PREFINE (DEPRECATED) *********************** //

	if(false)//es->parameters.get<unsigned int>("prerefine"))
	{
			
		if(sim_3d)
		{
			if(es->parameters.get<unsigned int>("prerefine") == 1)
				prerefine_3d_mesh();
			else if(es->parameters.get<unsigned int>("prerefine") == 2)
				adaptively_refine();
			else
			{
				std::cout << "type of prefinement " << es->parameters.get<unsigned int>("prerefine") 
					<< " does not exist. Aborting program." << std::endl;
				std::exit(0);
			}
		}
		else
			std::cout << "WE DO NOT REFINE 1D MESHES" << std::endl;
	}
	else
	{
		std::cout << "not doing prerefine anymore" << std::endl;
	}








	// **************************** SET SOME DATA WE USE *********************** //




	// set the 0d radii
	if(sim_1d)
	{
		std::cout << "Setting 1d mesh radii" << std::endl;
		set_radii(extra_1d_data_system);
		set_elem_proc_id_1d(extra_1d_data_system);


		// who cares??
		/*
		if(es->parameters.get<unsigned int>("acinar_model") == 1)
		{
			set_elem_proc_id_acinar(extra_1d_data_system);
		}
		*/
	}

	if(sim_3d)
	{
		if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))
		{
			std::cout << "Setting centreline radii" << std::endl;
			set_radii(system_centrelines,true);			
		}
		set_elem_proc_id_3d();
		if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))
		{
			set_elem_proc_id_1d(system_centrelines,true);
		}
	}











	// init dof variable vectors (first used when reading in solutions)
	init_dof_variable_vectors();

	std::cout << "after init dof var vectors" << std::endl;
	
	// after things have been setup let us set up the variable scalings for output
	if(sim_3d)
		setup_variable_scalings_3D();

	if(sim_1d)
	{
		setup_variable_scalings_1D();
		if(es->parameters.get<unsigned int>("acinar_model") == 1)
		{
			setup_variable_scalings_acinar();
		}		
	}
	std::cout << "after variable scalings" << std::endl;
	std::cout << std::endl;

	// can't remember what this is, 
	// but i think it makes calculating fluxes and pressures on boundaries more efficient
	perf_log_setup.push("boundary vectors compute");
	if(sim_3d)
		picard->setup_flow_rate_and_mean_pressure_vectors();
	perf_log_setup.pop("boundary vectors compute");










	// ********* READ IN OLD FLUID SIMULATIONS FOR RESTART AND PARTICLE DEPOSITION ************ //


	//std::cout << "hello" << std::endl;
	// **** READ IN TIMESTEPS FOR UNSTEADY RESTART ******** //
	if(particle_deposition == 1 || particle_deposition == 3 || particle_deposition == 4 || particle_deposition == 5 || particle_deposition == 6)
		setup_read_timesteps();


	// **** ACTUALLY READ IN OLD SOLUTIONS **** //
	if(restart)
	{
		std::cout << "READING IN OLD SOLUTIONS AT t=0" << std::endl;
	
		// if reading from a steady simulation then need to read timestep 1 no matter what

		unsigned int read_time_step = 0;
		double read_time = 0;

		if(es->parameters.get<bool> ("unsteady_from_steady"))
			read_time_step = 1;
		else	// normal unsteady restart of particle deposition
		{
			read_time_step = es->parameters.get<unsigned int> ("restart_time_step");
			read_time = es->parameters.get<double> ("restart_time");
		}


		// read in the solution/s at read_time_step
		read_old_solutions(read_time_step, read_time);

		// first time step, so copy the initial to the previous as well
		// should really read in two time steps, but hey
		if(sim_1d)
			copy_back_1d_solutions();
		if(sim_3d)
			copy_back_3d_solutions();

		// tell the program where these solutions came from
		read_time_1 = read_time;
		read_time_step_1 = read_time_step;
		read_time_2 = read_time;
		read_time_step_2 = read_time_step;
		
	}// end reading in old solutions

	//std::cout << "HEY BAE" << std::endl;







	// ************* READ IN BOUNDARY CONDITIONS FOR PROBLEM ****** //
	// if you want to run an uncoupled simulation from a previous coupled simulation
	// OR if you want to specify the pressure boundary conditions for example from a dirichlet problem
	if(sim_3d)
	{
		if(es->parameters.get<bool>("known_boundary_conditions")) // && sim_type == 2 
		{

			if(sim_type == 0)
			{
				// read the input boundary conditions for a pressure
				read_input_boundary_conditions_3d();
			}
			if(sim_type == 2)
			{
				// read the input boundary conditions for an uncoupled simulation
				read_input_boundary_conditions_3d0d();
			}


		}
		else if(es->parameters.get<unsigned int>("problem_type") == 1)
		{
			// if we aren't reading in the boundary conditions and we are doing a pressure type simulation
			// then we need to setup the setup pressure boundary conditions
			setup_pressure_boundary_conditions();

		}
	}





	// ************** CALCULATE THE RESISTANCE OF THE 1D TREES ************** //
	if(sim_1d)//es->parameters.get<unsigned int>("moghadam_coupling") && sim_type == 2)
		calculate_1d_pressure_deriv_values();


	// set reynolds number (could be based on the solution that was read in, but not currently)
	set_reynolds_number();

	

	double beforebefore_norm = 0.;
	if(!particle_deposition)
	{
		if(sim_3d)
		{


			if(sim_type != 5)
				beforebefore_norm = system_3d->solution->l2_norm();
			else
				beforebefore_norm = system_coupled->solution->l2_norm();					
		}
	}












	// ************* OUTPUT INFO TO SCREEN AND TO FILE ********** //

	std::cout << "\n\n\n---------- Program Mesh info ----------" << std::endl;
	mesh.print_info();
	mesh.boundary_info->print_summary();

	if(sim_3d)
	{
		std::cout << "\n\n\n---------- 3D Output Mesh info ----------" << std::endl;
		mesh_3d.print_info();
		mesh_3d.boundary_info->print_summary();

		if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))
		{
			std::cout << "\n\n\n---------- Centrelines Mesh info ----------" << std::endl;
			mesh_centrelines.print_info();
			mesh_centrelines.boundary_info->print_summary();
		}
	}

	if(sim_1d)
	{
		std::cout << "\n\n\n---------- 1D Output Mesh info ----------" << std::endl;
		mesh_1d.print_info();
		mesh_1d.boundary_info->print_summary();

		if(es->parameters.get<unsigned int> ("additional_symmetric_generations"))
		{
			std::cout << "\n\n\n---------- 1D Terminal Output Mesh info ----------" << std::endl;
			mesh_1d_terminal.print_info();
			mesh_1d_terminal.boundary_info->print_summary();
		}

		if(es->parameters.get<unsigned int>("acinar_model") == 1)
		{
			std::cout << "\n\n\n---------- 1D Acinar Output Mesh info ----------" << std::endl;
			mesh_acinar.print_info();
			mesh_acinar.boundary_info->print_summary();

		}

	}




	std::cout << "\n\n\n---------- Program Equation System ----------" << std::endl;
	es->print_info();
	if(sim_3d)
	{
		std::cout << "\n\n\n---------- 3D Output Equation System ----------" << std::endl;
		es_3d->print_info();

		if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))
		{
			std::cout << "\n\n\n---------- 3D Centrelines Equation System ----------" << std::endl;
			es_centrelines->print_info();
		}
	}

	if(sim_1d)
	{
		std::cout << "\n\n\n---------- 1D Output Equation System ----------" << std::endl;
		es_1d->print_info();

		if(es->parameters.get<unsigned int>("acinar_model") == 1)
		{
			std::cout << "\n\n\n---------- Acinar Output Equation System ----------" << std::endl;
			es_acinar->print_info();
		}

	}



	perf_log_setup.disable_logging();
	perf_log.pop("setup");






	//PerfLog perf_log_misc("Miscellaneous");

	// *************** PARTICLE DEPOSITION ******************* //
	if(particle_deposition && !(particle_deposition == 2 || particle_deposition == 3))
	{

	 	// ************ TIME AND NONLINEAR LOOPS ***************************** //

		es->parameters.set<unsigned int> ("t_step") = 
			es->parameters.get<unsigned int>("restart_time_step");


		//loop over viscositys for steady state solution num continuation
		for (unsigned int k=0; k<re_vec.size(); ++k)
		{

		 	es->parameters.set<Real> ("reynolds_number")   = re_vec[k];
		 	es->parameters.set<unsigned int> ("num_continuation_iteration")   = k;

			std::cout << "CALCULATING FOR Re = " <<
				es->parameters.get<Real> ("reynolds_number") << std::endl;



			// *************** INITIALISE PARTICLES AT t=0 **************************


			// **** init 3d particles
			// for 1 - only 3d unsteady/steady, 4 - 3d-0d steady, 6 - 3d-0d unsteady
			if(particle_deposition == 1 || particle_deposition == 4 || particle_deposition == 6)
				init_particles();

			// **** init 0d particles
			// for 0d unsteady james particle deposition, we need to initialise particles at t=0
			if(particle_deposition == 5)
			{
				// calculate the deposition efficiencies, airway flow rates, velocities and stokes number and stuff, mostly for output
				james_particle_deposition_object.calculate_airway_deposition_probability();

				// if we are doing a 0d only calculation we need to send it a fraction each time steps
				// remember this is actually time step 1, some particles could have been deposited at time step 0
				// 0 - particle deposition all at once, 1 - particle deposition at a specified rate on a surface, 2 - particle deposition all at once until all particles deposited or exited, 3 - particle deposition at a constant rate until total number have entered system, 4 -particle deposition at a constant rate until total number have entered system then wait until all have deposited or exited
				double fraction_added = 0;
				fraction_added = init_0d_particles();
				total_fraction_added += fraction_added;
				particle_fraction_added_so_far = total_fraction_added;

				// add particle fractions takes in vectors so set these up
				// empty for the first elemetn
				std::vector<double> particle_fractions(2);
				particle_fractions[1] = fraction_added;

				std::cout << "YASfraction_added = " << fraction_added << std::endl;
				james_particle_deposition_object.add_particle_fractions_to_tree(particle_fractions,airway_elem_id_starts,time);
			}
			// end initialising particles

			perf_log.push("output");
			// end setting up output files



			// *************** WRITE PARAMETERS ******************** //
			//parameters
			std::ostringstream parameter_output_data_file_name;

			parameter_output_data_file_name << output_folder.str() << "parameters" << t_step << ".dat";
			//output_data_file_name << "results/out_viscosity"	<< es->parameters.set<Real> ("viscosity") << ".dat";
			parameter_output_file.open(parameter_output_data_file_name.str().c_str());
			parameter_output_file << "# Parameters" << std::endl;
			es->parameters.print(parameter_output_file);
			parameter_output_file.close();








			perf_log.pop("output");


			// ***************** OUTPUTTING ZEROTH TIME STEP ********************* //	

			// do this 3d output for all  except 0d unsteady (5)
			std::cout << "WRITING PARTICLES AT t=0" << std::endl;


			// **** calculate 3d particle data
			perf_log.push("calculate_3d_deposition");
			if(sim_3d)
				calculate_3d_particle_data();
			perf_log.pop("calculate_3d_deposition");

			// **** calculate 0d particle data
			perf_log.push("calculate_0d_deposition");
			if(sim_1d)
				calculate_0d_particle_data();
			perf_log.pop("calculate_0d_deposition");


			// **** output 3d
			perf_log.push("output_3d_deposition");

			// if 1 - unsteady 3d, 4 - steady 3d0d, 6 - unsteady 3d0d the output
			if(particle_deposition == 1 || particle_deposition == 4 || particle_deposition == 6)
			{
				output_3d_particle_data();
				write_3d_particles();

				if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))	
					write_centrelines();
			}

			if(sim_3d && es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id") )
				output_3d_airway_deposition_data();

			perf_log.pop("output_3d_deposition");
	

			// **** output 0d
			perf_log.push("output_0d_deposition");

			// if 5 - unsteady 0d, 6 unsteady 3d0d
			if(particle_deposition == 5 || particle_deposition == 6)
				write_1d_solution();

			// output airway data
			if(sim_1d)
				output_0d_airway_deposition_data();

			perf_log.pop("output_0d_deposition");



			// **** output coupled

			//header
			output_global_deposition_metrics(true);
			// first time step
			output_global_deposition_metrics(false);


			// **** write the airway tree files
			output_airway_tree();



			// end outputting stuff






			// ************* TIME LOOP ******************** //
			// let us do this in dimensional terms
			while (time + 1e-10 < es->parameters.get<Real> ("end_time"))
			{





			
				// ************** INCREMENT TIME ************** //s

				perf_log.push("update times");
				++t_step;


				// choose dt
				/* don't do this, we want the possibility of a smaller timestep
				if(particle_deposition && !es->parameters.get<bool> ("unsteady_from_steady"))
					dt = timestep_sizes[t_step - 1];	//set timestep based on file
				*/
	
				if(time + dt > es->parameters.get<Real> ("end_time"))
				{
					dt = es->parameters.get<Real> ("end_time") - time;
				}

				time += dt;
	
				// if doing a precalculation to setup an initial condition e.g. stokes, then don't increment time
				if(es->parameters.get<bool> ("stokes_ic") && !ic_set)
				{
					--t_step;
					time -= dt;
					ic_set = true;
				}

				update_times();

				perf_log.pop("update times");

				//shell_pc_created = false;
				//mono_shell_pc_created = false;







				std::cout << "\n\n*** Solving time step " << t_step <<
				             ", time = " << time << " (" << time*time_scale_factor <<
										 "s) ***" << std::endl;


				// ******** READ IN SOLUTION AND/OR CALCULATE CURRENT FLOW RATE / VELOCITY

				if(es->parameters.get<bool>("unsteady_from_steady"))
				{
					//do nothing
				}
				else
				{
					// we need to figure out is we have passed the next time step 
					// and need to read in another solution
					
					// the last time step we read in was read_time_step_2
					double next_dt = 0.;
					if(read_time_step_2 + 1 < timestep_sizes.size())
						next_dt = timestep_sizes[read_time_step_2 + 1];
					
					// so the next time we need to read in some is
					double next_read_time = read_time_2 + next_dt;

					//std::cout << "next_dt = " << next_dt << std::endl;
					//std::cout << "next_read_time = " << next_read_time << std::endl;
					//std::cout << "read_time_2 = " << read_time_2 << std::endl;


					// if we've passed this time, we need to read in a new solution
					if(time > (read_time_2 + 1e-10))
					{
						std::cout << "READING IN OLD SOLUTION AT t = " << next_read_time << std::endl;			

						// now let's read in the solution
						// it's going to copy over only the velocity and fluxes hopefully
						perf_log.push("read_solution");

						// read in the solution/s into Q
						read_old_solutions(read_time_step_2 + 1, next_read_time);


						// copy back solutions Q_2 -> Q_1, Q -> Q_1
						if(sim_1d)
							copy_back_1d_solutions();

						if(sim_3d)
							copy_back_3d_solutions();

						// tell the program where these solutions came from
						read_time_1 = read_time_2;
						read_time_step_1++;
						read_time_2 = next_read_time;
						read_time_step_2++;
						perf_log.pop("read_solution");
					}
					else
					{
						// do nothing
					}

					perf_log.push("calculate_local_flow_solution_3d");
					if(sim_3d)
						calculate_local_particle_3d_flow_rate();
					perf_log.pop("calculate_local_flow_solution_3d");

					// only need to do if actually calculating it.
					if (time - ((int)((time + 1e-10) /es->parameters.get<Real>("particle_0d_dt")))
										*es->parameters.get<Real>("particle_0d_dt") < dt - 1e-10)// (t_step)%write_interval == 0)
					{
						perf_log.push("calculate_local_flow_solution_1d");
						// once we have read in or not, calculate the local particle flow rate
						if(sim_1d)
							calculate_local_particle_1d_flow_rate();

						perf_log.pop("calculate_local_flow_solution_1d");
					}

				}// end reading in old solutions











				// **************** DO 3D PARTICLE DEPOSITION ******************** //
				// 1 is 3d unsteady only, 4 is 3d-0d steady, 6 is coupled 3d-0d unsteady
				if(particle_deposition == 1 || particle_deposition == 4 || particle_deposition == 6)
				{

					// **** init 3d particles
					perf_log.push("init_3d_particles");
					// intialise new particles for this time step (not the first time step though) XXXX yes the first time step
					if(//t_step != 1 && 
							(es->parameters.get<unsigned int> ("particle_deposition_type") == 1
								|| es->parameters.get<unsigned int> ("particle_deposition_type") == 3
								|| es->parameters.get<unsigned int> ("particle_deposition_type") == 4
								|| es->parameters.get<unsigned int> ("particle_deposition_type") == 5))
					{
						init_particles();
					}
					perf_log.pop("init_3d_particles");




					// **** move 3d particles
					perf_log.push("move_3d");
					move_particles();
					perf_log.pop("move_3d");


					perf_log.push("calculate_3d_deposition");
					calculate_3d_particle_data();
					perf_log.pop("calculate_3d_deposition");


					// **** output 3d particles
					if (time - ((int)((time + 1e-10) /es->parameters.get<Real>("write_interval")))
										*es->parameters.get<Real>("write_interval") < dt - 1e-10)// (t_step)%write_interval == 0)
					{
						perf_log.push("output_3d_deposition");	
						write_3d_particles();
						if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))
							write_centrelines();
						perf_log.pop("output_3d_deposition");

						perf_log.push("output_3d_airway_deposition");
						if(sim_3d && es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id") )
							output_3d_airway_deposition_data();
						perf_log.pop("output_3d_airway_deposition");
					}


					
				}







				// ************ UNSTEADY JAMES PARTICLE DEPOSITION ******************* //
				if(particle_deposition == 5 || particle_deposition == 6)
				{

					std::cout << "\nADDING PARTICLE FRACTIONS TO 0D" << std::endl;
					// at every small time step we add particles from 3d to 0d or init particles
					if(particle_deposition == 6)
					{

						//std::cout << "Calculating coupling terms for coupling 3D to 0D" << std::endl;
						/*
						std::cout << "using " << fraction_exited_3d_surface.size() << " 3D exit fractions for this time step:" << std::endl;
						for(unsigned int i=0; i<fraction_exited_3d_surface.size(); i++)
							std::cout << " " << fraction_exited_3d_surface[i] << std::endl;
						*/
			
						james_particle_deposition_object.add_particle_fractions_to_tree(fraction_exited_3d_surface,airway_elem_id_starts,time);


					}
					else if(particle_deposition == 5)	// just 0D unsteady
					{
						double fraction_added = 0.;
						fraction_added = init_0d_particles();
						total_fraction_added += fraction_added;

					
						std::vector<double> input_particle_density;
						input_particle_density.push_back(0.);	// first element corresponds to the 3D
						input_particle_density.push_back(fraction_added);	// should only be one really, unless you're dumb
					
						james_particle_deposition_object.add_particle_fractions_to_tree(input_particle_density,airway_elem_id_starts,time);


					}



					// if we are at the correct time step do 0d deposition time step
					if (time - ((int)((time + 1e-10) /es->parameters.get<Real>("particle_0d_dt")))
										*es->parameters.get<Real>("particle_0d_dt") < dt - 1e-10)// (t_step)%write_interval == 0)
					{
						std::cout << "\nCALCULATING 0D UNSTEADY PARTICLE DEPOSITION" << std::endl;


						// **** setup 0d deposition
						perf_log.push("setup_0d_deposition");

						// set the airway deposition probability, flow rates, velocities etc
						james_particle_deposition_object.calculate_airway_deposition_probability();

						perf_log.pop("setup_0d_deposition");



						// **** move 3d particles
						perf_log.push("move_0d_deposition");
						// do the deposition
						james_particle_deposition_object.calculate_time_step_deposition_fraction(airway_elem_id_starts,time,es->parameters.get<double>("particle_0d_dt"));

						perf_log.pop("move_0d_deposition");

						perf_log.push("calculate_0d_deposition");
						calculate_0d_particle_data();
						perf_log.pop("calculate_0d_deposition");


						if (time - ((int)((time + 1e-10) /es->parameters.get<Real>("write_interval")))
											*es->parameters.get<Real>("write_interval") < dt - 1e-10)// (t_step)%write_interval == 0)
						{
							perf_log.push("output_0d_deposition");
							write_1d_solution();
							perf_log.pop("output_0d_deposition");

							perf_log.push("output_0d_airway_deposition");
							if(sim_1d)
								output_0d_airway_deposition_data();
							perf_log.pop("output_0d_airway_deposition");
	 	

						}
					}// end particle 0d dt
				}









				// once we have all the information from coupled sims, we can output global metrics and airway stuff (could be done before i guess) to file
				// can only do this on 0d time steps because otherwise there might be particles that are left in between
				if(particle_deposition)
				{
					if (time - ((int)((time + 1e-10) /es->parameters.get<double>("particle_write_interval")))
										*es->parameters.get<double>("particle_write_interval") < dt - 1e-10)// (t_step)%write_interval == 0)
					{
						output_global_deposition_metrics(false);						
					}
				}






				// **** CHECK THERE ARE NO PARTICLES IN MOTION ********** //
				if(particle_deposition)
				{

					// - don't actually have to run this everytime, only when we have just written output 
					// and are actually allowed to end the sim
					// - if we are doing a deposition limited simulation, 
					// - need to check if all particles have been deposited
					// - and if we have passed the time at which new particles enter the system
					if((time - ((int)((time + 1e-10) /es->parameters.get<Real>("backup_write_interval")))
										*es->parameters.get<Real>("backup_write_interval") < dt - 1e-10)
						&& (es->parameters.get<unsigned int> ("particle_deposition_type") == 2
							|| ( (es->parameters.get<unsigned int> ("particle_deposition_type") == 4 
										|| es->parameters.get<unsigned int> ("particle_deposition_type") == 5)
							&& time > es->parameters.get<double> ("particle_entry_end_time") + 1e-10)) )
					{
						std::cout << "Checking if any particles are still in motion." << std::endl;
						if(no_particles_in_motion())
							break;	// end sim after next output
					}
				}



				// ************** ADMIN AT END OF TIME STEP


				perf_log.push("end time step admin");
				end_timestep_admin();		//decide on what to do next and get vectors ready
				perf_log.pop("end time step admin");



			}// end timestep loop.








			// ****** CLOSE OUTPUT FILES ******** //

			if(particle_deposition)
			{
				particle_output_file.close();
				global_deposition_metrics_file.close();
			}

		} //end viscosity loop



	}// end particle deposition






	// ********* FLUID STUFF ************* //

	if(!es->parameters.get<bool>("compare_results") && !particle_deposition)
	{

	 	// ************ TIME AND NONLINEAR LOOPS ***************************** //

		es->parameters.set<unsigned int> ("t_step") = 
			es->parameters.get<unsigned int>("restart_time_step");


		//loop over viscositys for steady state solution num continuation
		for (unsigned int k=0; k<re_vec.size(); ++k)
		{

		 	es->parameters.set<Real> ("reynolds_number")   = re_vec[k];
		 	es->parameters.set<unsigned int> ("num_continuation_iteration")   = k;

			std::cout << "CALCULATING FOR Re = " <<
				es->parameters.get<Real> ("reynolds_number") << std::endl;










			// *************** WRITE OUTPUT ******************** //

			// *************** SETTING UP STUFF TO WRITE ******************** //
			if(!particle_deposition)
			{
				// output the parameters used to file and the header of the out.dat file
				std::ostringstream output_data_file_name;

				output_data_file_name << output_folder.str() << "out.dat";
				//output_data_file_name << "results/out_viscosity"	<< es->parameters.set<Real> ("viscosity") << ".dat";
				output_file.open(output_data_file_name.str().c_str());
			}

			if(es->parameters.get<bool> ("write_eigenvalues"))
			{
				// output the parameters used to file and the header of the out.dat file
				std::ostringstream eigenvalues_data_file_name;

				eigenvalues_data_file_name << output_folder.str() << "eigenvalues.dat";
				//output_data_file_name << "results/out_viscosity"	<< es->parameters.set<Real> ("viscosity") << ".dat";
				eigenvalues_file.open(eigenvalues_data_file_name.str().c_str());

				eigenvalues_file << "# eigenvalue data file" << std::endl;
				eigenvalues_file << "# timestep_no\tnonlinear_it\teigenvalue_real\teigenvalue_imag" << std::endl;

				// output the parameters used to file and the header of the out.dat file
				std::ostringstream eigenvalues_approx_data_file_name;

				eigenvalues_approx_data_file_name << output_folder.str() << "eigenvalues_approx.dat";
				//output_data_file_name << "results/out_viscosity"	<< es->parameters.set<Real> ("viscosity") << ".dat";
				eigenvalues_approx_file.open(eigenvalues_approx_data_file_name.str().c_str());

				eigenvalues_approx_file << "# top 10 approx eigenvalue data file" << std::endl;
				eigenvalues_approx_file << "# timestep_no\tnonlinear_it\teigenvalue_real\teigenvalue_imag" << std::endl;

			}

			if(es->parameters.get<bool> ("output_linear_iteration_count") && !particle_deposition)
			{
				std::ostringstream linear_iterations_data_file_name;

				linear_iterations_data_file_name << output_folder.str() << "linear_iterations.dat";
				//output_data_file_name << "results/out_viscosity"	<< es->parameters.set<Real> ("viscosity") << ".dat";
				linear_iterations_output_file.open(linear_iterations_data_file_name.str().c_str());
				if(!restart)
					output_linear_iteration_count(true);


			}

			if(!restart)
			{
				output_sim_data(true);
			}

			perf_log.push("output");
			// end setting up output files



			// *************** WRITE PARAMETERS ******************** //
			//parameters
			std::ostringstream parameter_output_data_file_name;

			parameter_output_data_file_name << output_folder.str() << "parameters" << t_step << ".dat";
			//output_data_file_name << "results/out_viscosity"	<< es->parameters.set<Real> ("viscosity") << ".dat";
			parameter_output_file.open(parameter_output_data_file_name.str().c_str());
			parameter_output_file << "# Parameters" << std::endl;
			es->parameters.print(parameter_output_file);
			parameter_output_file.close();








			perf_log.pop("output");

			//std::cout << "boo" << std::endl;





			// ************* ZEROING STUFF FOR THE BEGINNING IF NECESSARY ************** //

			// only do this when on the first step of num continuation
			if(k==0)
			{
				// must zero on first time step - after boundary conditions updated
				// do we really?
				if(sim_3d)
				{
					if(sim_type != 5 && !restart)
					{
						
						// don't really want to do this because we want t=0 bcs applied
						if(!es->parameters.get<bool> ("residual_formulation") || !es->parameters.get<bool> ("set_old_solution"))
						{
							system_3d->solution->zero();
							system_3d->old_local_solution->zero();
						}
						

						//if we are doing womersley we must set the inflow stuffs
						if(es->parameters.get<bool> ("prescribed_womersley"))
						{
							InflowBoundaryFunction<>* p_inflow;
							p_inflow = new InflowBoundaryFunction<>(*es);
							system_3d->project_solution(p_inflow);
						}
						system_3d->update();
					}
					else if(!restart)
					{
						
						system_coupled->solution->zero();
						system_coupled->old_local_solution->zero();
						if(es->parameters.get<unsigned int>("acinar_model") == 1 && es->parameters.get<bool>("init_acinar_volumes"))
						{
							init_acinar_volumes();
						}
						system_coupled->update();
						
					}
				}

				// must zero on first time step - might have been changed when calculating 1d resistances
				if(sim_1d)
				{
					if(sim_type != 5 && !restart)
					{
						
						system_1d->solution->zero();
						system_1d->old_local_solution->zero();
						if(es->parameters.get<unsigned int>("acinar_model") == 1 && es->parameters.get<bool>("init_acinar_volumes"))
						{
							init_acinar_volumes();
						}						
						system_1d->update();
					}

				}
			}// end zeroing stuff







			// ***************** INITING PREVIOUS TIME SOLUTION ***************** //

			if(sim_3d && sim_type != 5)
			{
				*system_3d->old_local_solution = *system_3d->current_local_solution;
				*old_global_solution = *system_3d->solution;
				//system_3d->update();
			}

			if(sim_1d && sim_type != 5)
			{
				*system_1d->old_local_solution = *system_1d->current_local_solution;
				*old_global_solution_1d = *system_1d->solution;
			}
		
			if(sim_type == 5)
			{
				*system_coupled->old_local_solution = *system_coupled->current_local_solution;
				*old_global_solution = *system_coupled->solution;
			}






			std::cout << "hi" << std::endl;
			// ***************** OUTPUTTING ZEROTH TIME STEP ********************* //	

			perf_log.push("output");
			// fluid output
			if(unsteady && !restart)
			{
				if(sim_3d) 
				{
					write_3d_solution();	// only fluid
					if(es->parameters.get<bool>("use_centreline_data") && es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))	
						write_centrelines();
				}
				if(sim_1d) 
				{
					write_1d_solution();
					if(es->parameters.get<unsigned int>("acinar_model") == 1)
					{
						write_acinar_solution();
					}
				}	// only fluid
			}









			// ************* TIME LOOP ******************** //
			// let us do this in dimensional terms
			while (time + 1e-10 < es->parameters.get<Real> ("end_time"))
			{





			
				// ************** INCREMENT TIME ************** //s

				perf_log.push("update times");
				++t_step;


				// choose dt
				/* don't do this, we want the possibility of a smaller timestep
				if(particle_deposition && !es->parameters.get<bool> ("unsteady_from_steady"))
					dt = timestep_sizes[t_step - 1];	//set timestep based on file
				*/
	
				if(time + dt > es->parameters.get<Real> ("end_time"))
				{
					dt = es->parameters.get<Real> ("end_time") - time;
				}

				time += dt;
	
				// if doing a precalculation to setup an initial condition e.g. stokes, then don't increment time
				if(es->parameters.get<bool> ("stokes_ic") && !ic_set)
				{
					--t_step;
					time -= dt;
					ic_set = true;
				}

				update_times();

				perf_log.pop("update times");

				//shell_pc_created = false;
				//mono_shell_pc_created = false;







				std::cout << "\n\n*** Solving time step " << t_step <<
				             ", time = " << time << " (" << time*time_scale_factor <<
										 "s) ***" << std::endl;







				// ************ READ IN OLD SOLUTION HERE AND CALCULATE THE FLOW RATE AT THIS TIME STEP ************ //
				// we want to read in the solution at the beginning of the time step and not 
				//  in end_timestep_admin is previously

				// ************ UPDATE TIME DEPENDENT DIRICHLET BOUNDARY CONDITIONS *********** //

				//std::cout << "read_time_2 = " << read_time_2 << std::endl;


				if(!particle_deposition)
				{
					if(sim_3d)
					{
						perf_log.push("reinit");
						update_3d_dirichlet_boundary_conditions();	// UPDATE TIME DEPENDENT DIRICHLET BOUNDARY CONDITIONS
						perf_log.pop("reinit");
						perf_log.push("calc norm");	
						if(sim_type != 5)
	 						beforebefore_norm = system_3d->solution->l2_norm();
						else
	 						beforebefore_norm = system_coupled->solution->l2_norm();
						perf_log.pop("calc norm");

					}
				}

				//std::cout << "l2_norm = " << beforebefore_norm << std::endl;
				
				//write_3d_solution(true);









				// ************** FLUID SOLVING AND OUTPUTTING ***************//

				if(!particle_deposition)
				{



					// ************** SOLVE 3D SYSTEM ********************* //
					if(sim_3d && sim_type != 5)
					{
						nonlinear_iteration = 0;
						nonlinear_iteration_1d = 0;
						es->parameters.set<unsigned int> ("nonlinear_iteration") = nonlinear_iteration;
						local_linear_iterations = 0;
						local_max_iterations = 0;
						local_max_residual_iterations = 0;

						double shift_value = 0.;
					
						//TEMP : this is to test how the different meshes compare
						//double before_norm = system_3d->solution->l2_norm();
						while(true)
						{
							nonlinear_iteration++;
							std::cout << "*********************" << std::endl;
							std::cout << "Nonlinear iteration " << nonlinear_iteration << std::endl;
							std::cout << std::endl;

							es->parameters.set<unsigned int> ("nonlinear_iteration") = nonlinear_iteration;

							perf_log.push("calculate 1d boundary values");
							if(sim_1d)
								calculate_1d_boundary_values();
							perf_log.pop("calculate 1d boundary values");

							// sim type 3 is when we tightly couple the 1d sim to the 3d sim
							if(sim_type == 0 || sim_type == 2)
							{
					
								// calculate the pressure input boundary conditions for uncoupled simulation
	
								// note: don't actually need to do this every nonlinear iteration but okay
								if(es->parameters.get<bool>("known_boundary_conditions")) // &&sim_type == 2 && )	// get input pressure from file.
								{
									std::cout << "t_step = " << t_step << std::endl;
										std::cout << "all_input_pressure_values_3d.size() = " << all_input_pressure_values_3d.size() << std::endl;
										std::cout << "all_input_pressure_values_3d[t_step - 1].size() = " << all_input_pressure_values_3d[t_step - 1].size() << std::endl;
									for(unsigned int i=0; i<input_pressure_values_3d.size(); i++)
									{
										input_pressure_values_3d[i] = all_input_pressure_values_3d[t_step - 1][i];	// t_step-1 because first timestep not output
										if(es->parameters.get<bool> ("apply_pressure_from_flux"))
											input_flux_values_3d[i] = all_input_flux_values_3d[t_step - 1][i];	// t_step-1 because first timestep not output

									}

									// when applying pressure from flux, need to calc boundary vals each time
									if(es->parameters.get<bool> ("apply_pressure_from_flux"))
										calculate_3d_boundary_values();
								}
								else if(es->parameters.get<unsigned int>("problem_type") == 1)
								{
									// calculates input_pressure_values_3d based on time scaling
									setup_pressure_boundary_conditions();
								}

								// for the dynamic pressure
								calculate_3d_boundary_values();


								// set bcs and solve 
								std::vector<double> empty_vec;
								if(es->parameters.get<unsigned int>("moghadam_coupling"))
								{
									// when applying moghadam resistance bc, need to update at each nonlinear iteration
									// solve the 0d problem to find the resistance, because we use a resistance bc, 
									// we don't know what the equivalent applied pressure would be
									if(sim_type == 2)
										calculate_1d_pressure_deriv_values();
									//picard->set_pressure_zero_flow_values_1d(pressure_zero_flow_values_1d);
									// use the pressure_deriv_values_1d and/or input_pressure_values_3d (well... we need the 1d pressures for moghadam, so yeah..)
									picard->init_bc(boundary_ids,pressure_values_1d,flux_values_1d,pressure_deriv_values_1d,
																	previous_flux_values_3d,previous_previous_flux_values_3d,
																	empty_vec,empty_vec,previous_dynamic_pressure_values_3d);
								}
								else if(es->parameters.get<bool> ("apply_pressure_from_flux"))
								{
									// for first step, we use the input pressure that we know, and send empty flow rate vec
									if(t_step == 1 && nonlinear_iteration == 1)
										picard->init_bc(boundary_ids,input_pressure_values_3d,input_flux_values_3d,pressure_deriv_values_1d,
																		flux_values_3d,previous_previous_flux_values_3d,
																		previous_pressure_values_3d,previous_previous_pressure_values_3d,
																		previous_dynamic_pressure_values_3d); 
									else
										picard->init_bc(boundary_ids,pressure_values_3d,input_flux_values_3d,pressure_deriv_values_1d,
																		flux_values_3d,previous_previous_flux_values_3d,
																		previous_pressure_values_3d,previous_previous_pressure_values_3d,
																		previous_dynamic_pressure_values_3d); 
								}
								else
									picard->init_bc(boundary_ids,input_pressure_values_3d,empty_vec,pressure_deriv_values_1d,
																	previous_flux_values_3d,previous_previous_flux_values_3d,
																	empty_vec,empty_vec,previous_dynamic_pressure_values_3d); 

								if(solve_3d_system_iteration(system_3d))
									break;


							}
							else if(sim_type == 3)
							{
								nonlinear_iteration_1d++;
								es->parameters.set<unsigned int> ("nonlinear_iteration_1d") = nonlinear_iteration_1d;


								// if we want to shift the boundary conditions, use the last element as the zero
								// from the first nonlinear iterations, well the second, after the first 3d solution has been solved.
								if(es->parameters.get<bool> ("shift_pressure_bc"))
								{
									if(nonlinear_iteration_1d == 2)
										shift_value = pressure_values_1d[pressure_values_1d.size()-1];

									for(unsigned int i=1; i< pressure_values_1d.size(); i++)
									{
										pressure_values_1d[i] = pressure_values_1d[i] - shift_value;			
									}
								}
								

								picard->init_bc(boundary_ids,pressure_values_1d,flux_values_1d,pressure_deriv_values_1d,
																previous_flux_values_3d,previous_previous_flux_values_3d,
																empty_vec,empty_vec,previous_dynamic_pressure_values_3d);

								if(solve_3d_system_iteration(system_3d))
									break;
				
								perf_log.push("calculate 3d boundary values");
								calculate_3d_boundary_values();
								perf_log.pop("calculate 3d boundary values");

								ns_assembler->init_bc(flux_values_3d);
								for(unsigned int i=0; i<flux_values_3d.size(); i++)
									std::cout << "flux_values_3d[" << i << "] = " << flux_values_3d[i] << std::endl;


								solve_1d_system(system_1d);
						
								perf_log.push("calculate 3d boundary values");
								calculate_3d_boundary_values();
								perf_log.pop("calculate 3d boundary values");
								perf_log.push("calculate 3d boundary values");
								calculate_1d_boundary_values();
								perf_log.pop("calculate 3d boundary values");
								perf_log.push("output");
								print_flux_and_pressure();
								perf_log.pop("output");
					

							}
							else if(sim_type == 4)
							{
								//in the explicitly coupled case sim_type 4 use 1d pressure vals
								picard->init_bc(boundary_ids,pressure_values_1d,flux_values_1d,pressure_deriv_values_1d,
																previous_flux_values_3d,previous_previous_flux_values_3d,
																empty_vec,empty_vec,previous_dynamic_pressure_values_3d); 
								if(solve_3d_system_iteration(system_3d))
									break;
							}

							//here we wanna output the boundary values that have been calculated and see if they are indeed converging
							if(exit_program)
								break;
						} // end nonlinear loop

						


						es->parameters.set<unsigned int> ("num_newton_steps") = nonlinear_iteration;
						total_nonlinear_iterations += nonlinear_iteration;
					}






					// ************** SOLVE 1D SYSTEM ********************* //
	
					// remember 1 is just 0D and 2 is uncoupled 0D and 3D
					// note: for moghadam we don't actually want to calculate this at the end again
					if((sim_type == 1 || sim_type == 2))// && !es->parameters.get<unsigned int>("moghadam_coupling"))
					{
						std::cout << "bye" << std::endl;

						//initialise vector to store flux input
						std::vector<double> input_flux_values(es->parameters.get<unsigned int>("num_1d_trees") + 1);

						std::cout << "input_flux_values.size() = " << input_flux_values.size() << std::endl;
						// calculate the flux input boundary conditions
						if(sim_type == 2 && es->parameters.get<bool>("known_boundary_conditions"))	// get input flux from file.
						{
							for(unsigned int i=0; i<input_flux_values.size(); i++)
							{
								input_flux_values[i] = all_input_flux_values_0d[t_step - 1][i];	//t_step-1 because first timestep not output
							}

						}
						// defunct
						else if(es->parameters.get<unsigned int>("moghadam_coupling"))	// just use the same input flux everywhere
						{
							calculate_3d_boundary_values();
							input_flux_values = flux_values_3d;
						}
						else
						{

							
							double inflow = es->parameters.get<double>("time_scaling");

							if(es->parameters.get<unsigned int>("unsteady"))
								inflow = es->parameters.get<double>("flow_mag_1d") * es->parameters.get<double>("time_scaling");//fabs(sin(2*M_PI*time/4));
							else	//put flow_mag into the nondimensional units
								inflow = es->parameters.get<double>("flow_mag_1d");	

							//input_flux_values[0] = 0.0;
						
							std::cout << "time_scaling = " << es->parameters.get<double>("time_scaling")  << std::endl;
							std::cout << "flow_mag_1d = " << es->parameters.get<double>("flow_mag_1d")  << std::endl;
							std::cout << "velocity_scale = " << es->parameters.get<double>("velocity_scale")  << std::endl;
							std::cout << "length_scale = " << es->parameters.get<double>("length_scale")  << std::endl;
							std::cout << "inflow = " << inflow  << std::endl;
							for(unsigned int i=0; i<input_flux_values.size(); i++)
							{
								input_flux_values[i] = inflow;
								std::cout << "inflow = " << inflow << std::endl;
							}
						}


						// solve the nonlinear 0D system

						ns_assembler->init_bc(input_flux_values);
						std::cout << "bye" << std::endl;
						
						solve_1d_system(system_1d);
						std::cout << "bye" << std::endl;
						
						perf_log.push("calculate 1d boundary values");
						calculate_1d_boundary_values();
						perf_log.pop("calculate 1d boundary values");
						
					}
					else if(sim_type == 4)
					{
						calculate_3d_boundary_values();
						ns_assembler->init_bc(flux_values_3d);
						solve_1d_system(system_1d);
					}










					// *********** MONOLITHIC SOLVING **************** //
					if(sim_type == 5)
					{
	
						std::cout << "hiya" << std::endl;
						local_linear_iterations = 0;
						local_max_iterations = 0;
						nonlinear_iteration = 0;
						es->parameters.set<unsigned int> ("nonlinear_iteration") = nonlinear_iteration;


						if(es->parameters.get<unsigned int>("problem_type") == 1)
						{
							// calculates input_pressure_values_3d based on time scaling
							setup_pressure_boundary_conditions();
						}


						while(true)
						{
							nonlinear_iteration++;
							es->parameters.set<unsigned int> ("nonlinear_iteration") = nonlinear_iteration;
							es->parameters.set<unsigned int> ("nonlinear_iteration_1d") = nonlinear_iteration;

							perf_log.push("calculate 1d boundary values");
							calculate_1d_boundary_values();
							perf_log.pop("calculate 1d boundary values");
							picard->init_bc(boundary_ids,input_pressure_values_3d,flux_values_1d,pressure_deriv_values_1d,
															previous_flux_values_3d,previous_previous_flux_values_3d,
															empty_vec,empty_vec,previous_dynamic_pressure_values_3d);
							std::cout << "about to calculate 3d boundary values" << std::endl;
							perf_log.push("calculate 3d boundary values");
							calculate_3d_boundary_values();
							perf_log.pop("calculate 3d boundary values");
							std::cout << "done about to calculating 3d boundary values" << std::endl;
							ns_assembler_mono->init_bc(flux_values_3d);

							if(solve_3d_system_iteration(system_coupled))
								break;
			
			
							//here we wanna output the boundary values that have been calculated and see if they are indeed converging
							if(exit_program)
								break;
						} // end nonlinear loop

						es->parameters.set<unsigned int> ("num_newton_steps") = nonlinear_iteration;
					}

					if(exit_program && !es->parameters.set<bool> ("output_then_exit"))
						break;


					// ************ WRITE OUTPUT *************** //		
			
					perf_log.push("output");

					std::cout << "yeah" << std::endl;						
					if((!reduce_dt || !unsteady) && !refine_mesh)
					{
						std::cout << "yeah" << std::endl;
						output_sim_data(false);
						std::cout << "yeah" << std::endl;

						if (!unsteady || (time - ((int)((time + 1e-10) /es->parameters.get<Real>("write_interval")))
										*es->parameters.get<Real>("write_interval") < dt - 1e-10) )// (t_step)%write_interval == 0)
						{
							if(sim_3d) 
							{
								write_3d_solution();

								if(es->parameters.get<bool>("use_centreline_data") 
										&& es->parameters.set<bool>("gmsh_diff_wall_bdy_id"))	
								{
									write_centrelines();
								}
							}
							if(sim_1d) 
							{
								write_1d_solution();
								if(es->parameters.get<unsigned int>("acinar_model") == 1)
								{
									write_acinar_solution();
								}
							}
						}
					std::cout << "yeah" << std::endl;
					}
			
					perf_log.pop("output");

					//exit(0);

					if(exit_program && es->parameters.set<bool> ("output_then_exit"))
						break;
		
				}// end fluid section







				// ************** ADMIN AT END OF TIME STEP

				std::cout << "time = " << time << std::endl;
				std::cout << "kk" << std::endl;

				perf_log.push("end time step admin");
				end_timestep_admin();		//decide on what to do next and get vectors ready
				perf_log.pop("end time step admin");

				std::cout << "AHHHHHHHHHH" << std::endl;

				//needs to be done afterwards so that solution is projected onto the correct space
				if(es->parameters.get<bool>("adaptive_refinement"))// && refine_mesh)// && !reduce_dt && sim_3d)// && refine_mesh)
				{
					std::cout << "adaptively_refine is true " << std::endl;
					adaptively_refine();
				}			

				std::cout << "time = " << time << std::endl;
				std::cout << "end_time = " << es->parameters.get<Real> ("end_time") << std::endl;

			}// end timestep loop.





			// ****** OUTPUT STUFF AT END OF FLUID SIM ************ //
			std::cout << "AVERAGE LINEAR ITERATION COUNT = " << (double)total_linear_iterations / (double)total_nonlinear_iterations;
			std::cout << "tot lin iterations = " << total_linear_iterations << std::endl;
			if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 13)
			{
				std::cout << "tot gmres lin iterations = " << total_gmres_linear_iterations << std::endl;
				std::cout << "tot cg lin iterations = " << total_cg_linear_iterations << std::endl;

				std::cout << "max gmres lin iterations = " << total_max_gmres_iterations << std::endl;
				std::cout << "max cg lin iterations = " << total_max_cg_iterations << std::endl;
			}
			std::cout << "tot nonlin iterations = " << total_nonlinear_iterations << std::endl;








			// ****** CLOSE OUTPUT FILES ******** //

			output_file.close();
			eigenvalues_file.close();
			eigenvalues_approx_file.close();
			if(particle_deposition)
			{
				particle_output_file.close();
				global_deposition_metrics_file.close();
			}

			if(es->parameters.get<bool> ("output_linear_iteration_count"))
 				linear_iterations_output_file.close();

		} //end viscosity loop




		// ********* CLEAN UP STUFF *********** //
		// petsc clean up if we have actually done some simulation, i.e. time > 0
		if(time > 1e-10 && es->parameters.get<unsigned int> ("particle_deposition") == 0)
			petsc_clean_up();

		output_logging();









		// ************* COMPARE EXACT SOLUTION (NOT TESTED FOR AGES) *************** //
		if(es->parameters.get<unsigned int> ("problem_type") == 5)
		{

			
			TransientLinearImplicitSystem * system;
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
		
			const unsigned int u_var = system->variable_number ("u");
			const unsigned int v_var = system->variable_number ("v");
			//hmmm yeah you know const an shit
			unsigned int w_var = 0;
			if(threed)
				w_var = system->variable_number ("w");
			const unsigned int p_var = system->variable_number ("p");

			double u_l2_norm,v_l2_norm,w_l2_norm,p_l2_norm;
			u_l2_norm = picard->calculate_l2_norm(u_var);
			v_l2_norm = picard->calculate_l2_norm(v_var);
			if(threed)
				w_l2_norm = picard->calculate_l2_norm(w_var);
			double vel_l2_norm = pow(u_l2_norm,2.) + pow(v_l2_norm,2.);
			if(threed)
				vel_l2_norm += pow(w_l2_norm,2.);
			vel_l2_norm = sqrt(vel_l2_norm);
			
			p_l2_norm = picard->calculate_l2_norm(p_var);

			double ref_u_l2_norm,ref_v_l2_norm,ref_w_l2_norm,ref_p_l2_norm;
			ref_u_l2_norm = picard->calculate_l2_norm(u_var,true);
			ref_v_l2_norm = picard->calculate_l2_norm(v_var,true);
			if(threed)
				ref_w_l2_norm = picard->calculate_l2_norm(w_var,true);
			double ref_vel_l2_norm = pow(ref_u_l2_norm,2.) + pow(ref_v_l2_norm,2.);
			if(threed)
				ref_vel_l2_norm += pow(ref_w_l2_norm,2.);
			ref_vel_l2_norm = sqrt(ref_vel_l2_norm);

			ref_p_l2_norm = picard->calculate_l2_norm(p_var,true);

			std::cout << "u_l2_norm = " << u_l2_norm << std::endl;

			std::cout << "ref_vel_l2_norm = " << ref_vel_l2_norm << std::endl;
			std::cout << "vel_l2_norm = " << vel_l2_norm << std::endl;
			std::cout << "vel_l2_norm_scaled = " << vel_l2_norm/ref_vel_l2_norm << std::endl;
			std::cout << "p_l2_norm = " << p_l2_norm << std::endl;
			std::cout << "ref_p_l2_norm = " << ref_p_l2_norm << std::endl;



		}// end compare exact solution

	}// end fluid sims and particle deposition that needs time steps








	// *********************** COMPARE RESULTS OF TWO PREVIOUSLY RUN SIMS (NOT TESTED IN AGES) ************************** //
	if (es->parameters.get<bool>("compare_results"))	//do the compare results stuff
	{

		/*
		std::ofstream norm_file;
		std::ostringstream norm_file_name;
		norm_file_name << output_folder.str() << "norms.dat";
		norm_file.open(norm_file_name.str().c_str());

		// get the two parameter files
		GetPot infile_1(es->parameters.get<std::string>("results_folder_1") + "navier.in");
		GetPot infile_2(es->parameters.get<std::string>("results_folder_2") + "navier.in");

		setup_results_vector(es->parameters.get<std::string>("results_folder_1"),results_vector_1);
		setup_results_vector(es->parameters.get<std::string>("results_folder_2"),results_vector_2);
		unsigned int num_results = 0;

		//do some error handling
		int n_timestep_skip = 1;
		double coarse_timestep = 0.;
		double fine_timestep = 0.;

		if(es->parameters.get<unsigned int>("unsteady"))
		{
			//get the number to skip in the fine time discretisation
			fine_timestep = infile_1("dt",1.0);
			coarse_timestep = infile_2("dt",1.0);
			n_timestep_skip = (int)(coarse_timestep/fine_timestep);
			num_results = std::min(results_vector_1.size(),results_vector_2.size());

			// set the dt of the es object to the coarse timestep
		}
		else
		{
			//use the minimum timesteps so that get a result somehow
			num_results = std::min(results_vector_1.size(),results_vector_2.size());
			n_timestep_skip = 1;	//if steady timestep skip is 1
		}


		// this may be wrong, i mean it is using the previous one
		TransientLinearImplicitSystem * system;	
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d");
	
		EquationSystems es_fine = *es;

		Mesh fine_mesh;
		es_fine.parameters.set<unsigned int>("no_refinement") = es->parameters.get<unsigned int>("fine_refinement");
		setup_3d_mesh(&es_fine,fine_mesh);

		es->reinit();
		
//			last_nonlinear_soln = system_coupled->solution->clone();
		NumericVector<Number>* exact_solution;// = *system->solution;

		NumericVector<Number>* difference;

		const unsigned int u_var = system->variable_number ("u");
		const unsigned int v_var = system->variable_number ("v");
		int w_var = -1;
		if(threed)
			w_var = system->variable_number ("w");
		const unsigned int p_var = system->variable_number ("p");

		double ref_norm_ul2 = 0.;	
		double norm_ul2 = 0.;	
		std::vector<double> norm_ul2_vec(num_results);
		double ref_norm_pl2 = 0.;	
		double norm_pl2 = 0.;	
		std::vector<double> norm_pl2_vec(num_results);
		double ref_norm_uh1_semi = 0.;	
		double norm_uh1_semi = 0.;	
		std::vector<double> norm_uh1_semi_vec(num_results);
		double norm_uh1 = 0.;	
		std::vector<double> norm_uh1_vec(num_results);
		FEMNormType l2_norm = L2;
		FEMNormType h1_norm = H1;

		norm_file << "# NORMS OVER TIME" << "\n";
		norm_file << "#time" << "\t";
		norm_file << "velocity_l2" << "\t";
		norm_file << "velocity_h1_semi" << "\t";
		norm_file << "pressure_l2" << "\t";
		norm_file << "velocity_h1" << "\n";

		//will have significant problems if don't have refined mesh!!!!!
		for(unsigned int i=0; i<num_results; i++)
		{
			std::cout << "i = " << i << std::endl; 

			// ******* CALCULATE REFERENCE NORMS *********** //
			*difference = *results_vector_1[n_timestep_skip * i];

			// velocity l2 ref norm
			ref_norm_ul2 = pow(system->calculate_norm(*difference,u_var,l2_norm),2.0);
			ref_norm_ul2 += pow(system->calculate_norm(*difference,v_var,l2_norm),2.0);
			if(threed)
				ref_norm_ul2 += pow(system->calculate_norm(*difference,w_var,l2_norm),2.0);
			ref_norm_ul2 = sqrt(ref_norm_ul2);

			// velocity h1 ref seminorm
			ref_norm_uh1_semi = pow(system->calculate_norm(*difference,u_var,h1_norm),2.0);
			ref_norm_uh1_semi += pow(system->calculate_norm(*difference,v_var,h1_norm),2.0);
			if(threed)
				ref_norm_uh1_semi += pow(system->calculate_norm(*difference,w_var,h1_norm),2.0);
			ref_norm_uh1_semi = sqrt(ref_norm_uh1_semi);

			// pressure l2 ref norm
			ref_norm_pl2 = system->calculate_norm(*difference,p_var,l2_norm);


			// ********* CALCULATE ABS NORMS ************** //

			*difference -=	*results_vector_2[i];

			// velocity l2 abs norm
			norm_ul2 = pow(system->calculate_norm(*difference,u_var,l2_norm),2.0);
			norm_ul2 += pow(system->calculate_norm(*difference,v_var,l2_norm),2.0);
			if(threed)
				norm_ul2 += pow(system->calculate_norm(*difference,w_var,l2_norm),2.0);
			norm_ul2 = sqrt(norm_ul2);

			// velocity h1 abs norm
			norm_uh1_semi = pow(system->calculate_norm(*difference,u_var,h1_norm),2.0);
			norm_uh1_semi += pow(system->calculate_norm(*difference,v_var,h1_norm),2.0);
			if(threed)
				norm_uh1_semi += pow(system->calculate_norm(*difference,w_var,h1_norm),2.0);
			norm_uh1_semi = sqrt(norm_uh1_semi);

			// pressure l2 abs norm
			norm_pl2 = system->calculate_norm(*difference,p_var,l2_norm);



			// ******* CALCULATE NORMALISED NORMS *********** //

			// normalised norms
			norm_ul2 = norm_ul2/ref_norm_ul2;
			norm_uh1_semi = norm_uh1_semi/ref_norm_uh1_semi;
			norm_pl2 = norm_pl2/ref_norm_pl2;
			norm_uh1 = norm_pl2 + norm_uh1_semi;

			// just set all of them to zero
			if(ref_norm_ul2 < 1e-10)
			{
				norm_ul2 = 0.;
				norm_uh1_semi = 0.;
				norm_pl2 = 0.;
				norm_uh1 = 0.;
			}

			norm_ul2_vec[i] = norm_ul2;
			norm_uh1_semi_vec[i] = norm_uh1_semi;
			norm_pl2_vec[i] = norm_pl2;
			norm_uh1_vec[i] = norm_uh1;

			// output the norms, i is the coarse time step number so use coarse timestep
			if(es->parameters.get<unsigned int>("unsteady"))
				norm_file << coarse_timestep*i << "\t";
			else
				norm_file << es->parameters.get<Real> ("dt")*i << "\t";
			norm_file << norm_ul2 << "\t";
			norm_file << norm_uh1_semi << "\t";
			norm_file << norm_pl2 << "\t";
			norm_file << norm_uh1 << "\n";

			
			
		}

		norm_file << "\n";
		norm_file << "# CUMULATIVE NORMS (trapezoidal)" << "\n";
		norm_file << "#velocity_l2" << "\t";
		norm_file << "velocity_h1_semi" << "\t";
		norm_file << "pressure_l2" << "\t";
		norm_file << "velocity_h1" << "\n";

		// sum by trapezoidal rule
		double cumulative_norm_ul2 = 0.;
		double cumulative_norm_uh1_semi = 0.;
		double cumulative_norm_pl2 = 0.;
		double cumulative_norm_uh1 = 0.;
		for(unsigned int i=0; i<num_results; i++)
		{
			if(i==0 || i==num_results-1)
			{
				cumulative_norm_ul2 += norm_ul2_vec[i];
				cumulative_norm_uh1_semi += norm_uh1_semi_vec[i];
				cumulative_norm_pl2 += norm_pl2_vec[i];
				cumulative_norm_uh1 += norm_uh1_vec[i];
			}
			else
			{
				cumulative_norm_ul2 += 2*norm_ul2_vec[i];
				cumulative_norm_uh1_semi += 2*norm_uh1_semi_vec[i];
				cumulative_norm_pl2 += 2*norm_pl2_vec[i];
				cumulative_norm_uh1 += 2*norm_uh1_vec[i];
			}
		}

		// assume simulation starts from zero
		std::cout << "cumulative_norm_ul2 = " << cumulative_norm_ul2 << std::endl;
		std::cout << "end time = " << es->parameters.get<double>("end_time") << std::endl;


		//cumulative_norm_ul2 *= es->parameters.get<double>("end_time") / (2 * num_results);
		//cumulative_norm_uh1_semi *= es->parameters.get<double>("end_time") / (2 * num_results);
		//cumulative_norm_pl2 *= es->parameters.get<double>("end_time") / (2 * num_results);
		//cumulative_norm_uh1 *= es->parameters.get<double>("end_time") / (2 * num_results);


		cumulative_norm_ul2 *= 1. / (2 * num_results);
		cumulative_norm_uh1_semi *= 1. / (2 * num_results);
		cumulative_norm_pl2 *= 1. / (2 * num_results);
		cumulative_norm_uh1 *= 1. / (2 * num_results);

		norm_file << cumulative_norm_ul2 << "\t";
		norm_file << cumulative_norm_uh1_semi << "\t";
		norm_file << cumulative_norm_pl2 << "\t";
		norm_file << cumulative_norm_uh1 << "\n";
		
		norm_file.close();

		*/
	}// end compare two previously run sims









	// ****************************  HOFMANN PARTICLE DEPOSITION (INCORRECT FORMULATION) ********************** //
	if (particle_deposition == 2)
	{

		HofmannParticleDeposition particle_deposition_object(*es,airway_data,num_generations[0],subdomains_1d);

		std::cout << "oh and hey babe" << std::endl;

		//intialise the flow
		particle_deposition_object.calculate_flow_rate();
		std::cout << "oh and hey again babe" << std::endl;

		//particle deposition file
		std::ostringstream deposition_output_data_file_name;

		deposition_output_data_file_name << output_folder.str() << "deposition.dat";
		//output_data_file_name << "results/out_viscosity"	<< es->parameters.set<Real> ("viscosity") << ".dat";
		deposition_output_file.open(deposition_output_data_file_name.str().c_str());

		std::cout << "yeah file = " << deposition_output_data_file_name.str().c_str() <<  std::endl;
		output_deposition_data(true);

		// ************* TIME LOOP ******************** //
		// let us do this in dimensional terms

		update_times();

		std::cout << "time = " << time << std::endl;
		std::cout << "es->parameters.get<Real> (\"end_time\") = " << es->parameters.get<Real> ("end_time") << std::endl;

		while (time + 1e-10 < 0.1)//es->parameters.get<Real> ("end_time"))
		{

			// INCREMENT TIME
			++t_step;

			if(time + dt > es->parameters.get<Real> ("end_time"))
			{
				dt = es->parameters.get<Real> ("end_time") - time;
			}
			time += dt;
			update_times();

			std::cout << "\n\n*** Solving time step " << t_step <<
			             ", time = " << time << " (" << time*time_scale_factor <<
									 "s) ***" << std::endl;

			//find efficiency of particle released at this timestep
			particle_deposition_object.calculate_deposition_efficiency(t_step);

			std::cout << "done one time step" << std::endl;
			//std::exit(0);
			
		}// end timestep loop.

		particle_deposition_object.calculate_total_efficiency_function();

		total_efficiency = particle_deposition_object.get_total_efficiency();
		alveolar_efficiency_per_generation = particle_deposition_object.get_acinar_efficiency_per_generation();
		tb_efficiency_per_generation = particle_deposition_object.get_tb_efficiency_per_generation();
		std::vector<double> non_deposited_particle_weight = particle_deposition_object.get_non_deposited_particle_weight();
		std::vector<double> output_particle_weight = particle_deposition_object.get_output_particle_weight();
		total_particles_inhaled = particle_deposition_object.get_total_particles_inhaled();

		std::cout << "non_deposited_particle_weight = " << non_deposited_particle_weight[1] << std::endl;
		std::cout << "output_particle_weight = " << output_particle_weight[1] << std::endl;

		for(unsigned int i=0; i<non_deposited_particle_weight.size(); i++)
		{
			std::cout << "non_deposited_particle_weight/tstep = " << non_deposited_particle_weight[i] << std::endl;
						
		}
			
		for(unsigned int i=0; i<output_particle_weight.size(); i++)
			std::cout << "output_particle_weight/tstep = " << output_particle_weight[i] << std::endl;


		set_efficiency();

		write_efficiency_solution();

		output_deposition_data(false);
		deposition_output_file.close();

		//output_file.close();

	}// end hofmann particle deposition












	// ************************** JAMES PARTICLE DEPOSITION STEADY ********************* //
	if(particle_deposition == 3 ||  particle_deposition == 4)
	{

		std::cout << "CALCULATING 0D STEADY PARTICLE DEPOSITION FOR Re = " <<
			es->parameters.get<Real> ("reynolds_number") << std::endl;


		//intialise the flow and associated airway deposition probability
		// don't need the particle fraction here
		james_particle_deposition_object.calculate_airway_deposition_probability();


		// setup the coupling between 3D and 0D
		if(particle_deposition == 4)
		{
			std::cout << "Setting coupling terms for coupling 3D to 0D" << std::endl;
			/*
			std::cout << "using " << total_fraction_exited_3d_surface.size() << " 3D exit fractions:" << std::endl;
			for(unsigned int i=0; i<total_fraction_exited_3d_surface.size(); i++)
				std::cout << " " << total_fraction_exited_3d_surface[i] << std::endl;
			*/


			//std::cout << "hello" << std::endl;
			// okay now we need to figure out where they go..
			// pass in the deposition fraction for the different airways
			james_particle_deposition_object.calculate_total_deposition_fraction(total_fraction_exited_3d_surface,airway_elem_id_starts);

		}
		else
		{
			// need to pass in the deposition fraction for the different airways
			james_particle_deposition_object.calculate_total_deposition_fraction();

		}







		std::cout << "you suck" << std::endl;

		double total_deposition_fraction = james_particle_deposition_object.get_total_deposition_fraction();
		double total_exit_fraction = james_particle_deposition_object.get_total_exit_fraction();

		std::cout << "total_deposition_fraction = " << total_deposition_fraction << std::endl;
		std::cout << "total_exit_fraction = " << total_exit_fraction << std::endl;
		std::cout << "total (exit+deposition) = " << total_deposition_fraction + total_exit_fraction << std::endl;


		// set time step to zero for output
		t_step = 0;
		time = 0;
		update_times();


		write_1d_solution();
		std::cout << "fuck off" << std::endl;
	}// end steady james particle deposition


	// write out parameters at the end
	print_parameters();







	if(exit_program)
		std::cout << "program aborted and exiting constructor" << std::endl;

	// barrier before end so that procs don't quit too early
	es->comm().barrier();

	perf_log.pop("misc");
}

// ******************************************************************* //






// ***************** CLASS FUNCTION DEFINITIONS ******************** //


// read in parameters from file and give them to the parameters object.
int NavierStokesCoupled::read_parameters()
{

	// - JAMES always 3D for this one
	//int dim = 3;

	// command line options from now on will not just be a letter but full desciptive 
	// word that overrides the value given in the input file.
	// -i input mesh
	// -ns navier stokes
	// -s stabilised	(0 and 1)
	// -a stability parameter
	// -v viscosity viscosity
	// -u unsteady	(0 and 1 or 2 for quadratic)
	// -r number of refinements
	// -nt number of time_steps
	// -p problem type (0 - dirichlet inflow, 1 - pressure, 2 - stress)
	// -new newton iterations
	// okay want to be able to choose solver options in the command line
	// could have mumps, superlu_dist, iterative(ml), iterative(hypre)
	// then there are other options within that we may want, but also want
	// to override at command_line - so ideally have it in a text file that
	// is read by the program at runtime and then implemented
	// -direct	using a direct solver
	// -adaptive use adaptive refinement

	std::cout << "-------------------------------------------" << std::endl;
	std::cout << "----------- SIMULATION DETAILS ------------" << std::endl;


	// ************** read in all the parameters ********************** //
	set_string_parameter(infile,"mesh_file","surface_mesh_no_holes_22.5_coarse.msh");
	sim_type = set_unsigned_int_parameter(infile,"sim_type",0);
	set_bool_parameter(infile,"stokes",true);
	set_bool_parameter(infile,"stab",true);
	set_double_parameter(infile,"alpha",5.0);
	// 0 - steady
	// 1 - sinusoidal
	// 2 - quadratic ramp up to 1 then constant for some time
	// 3 - cos shifted positive, 
	// 4 - negative sin, 
	// 5 - cos
	// 6 - hofmann breath then hold

	unsteady = set_unsigned_int_parameter(infile,"unsteady",2);
	//double viscosity = set_double_parameter(infile,"viscosity",1.0);
	set_double_parameter(infile,"viscosity",1.0);
	set_unsigned_int_parameter(infile,"no_refinement",0);
	set_unsigned_int_parameter(infile,"fine_refinement",0);
	dt = set_double_parameter(infile,"dt",0.1);
	set_double_parameter(infile,"previous_dt",0.1);
	set_double_parameter(infile,"end_time",1);
	// problem type (basically bc type)
	// 0 - dirichlet-neumann
	// 1 - pressure-pressure
	set_unsigned_int_parameter(infile,"problem_type",0);
	set_unsigned_int_parameter(infile,"newton",0);
	set_bool_parameter(infile,"direct",true);
	set_unsigned_int_parameter(infile,"prerefine",0);
	set_double_parameter(infile,"restart_time",0.0);
	bool restart_time_step = set_unsigned_int_parameter(infile,"restart_time_step",0);
	bool adaptive_refinement = set_bool_parameter(infile,"adaptive_refinement",false);
	set_bool_parameter(infile,"adaptive_time_stepping",false);
	set_double_parameter(infile,"nonlinear_tolerance",1.e-4);
	set_double_parameter(infile,"outer_solver_rtol",1.e-8);
	set_double_parameter(infile,"outer_solver_atol",1.e-8);
	set_int_parameter(infile,"outer_solver_maxits",40);
	set_double_parameter(infile,"refine_fraction",0.3);
	set_double_parameter(infile,"coarsen_fraction",0.3);
	set_double_parameter(infile,"coarsen_threshold",1.0);
	set_unsigned_int_parameter(infile,"max_h_level",2);
	set_unsigned_int_parameter(infile,"nelem_target",1000);
	set_unsigned_int_parameter(infile,"face_level_mismatch_limit",0);
	set_double_parameter(infile,"write_interval",0.1);
	set_double_parameter(infile,"backup_write_interval",es->parameters.get<double> ("write_interval")*10);
	set_unsigned_int_parameter(infile,"min_steps_between_dt_increase",0);	
	set_double_parameter(infile,"period",2.0);	
	//double density = set_double_parameter(infile,"density",1.176e-6);
	set_double_parameter(infile,"density",1.176e-6);
	//set_double_parameter(infile,"viscosity",15.23); //this is now the same
	set_double_parameter(infile,"zeta_1",-0.057);
	set_double_parameter(infile,"zeta_2",0.2096);
	set_double_parameter(infile,"zeta_3",0.00904);
	set_double_parameter(infile,"E",3.3);
	/* DEPRECATED FOR num_generations_string
	set_unsigned_int_parameter(infile,"num_generations_1",2);
	set_unsigned_int_parameter(infile,"num_generations_2",5);
	set_unsigned_int_parameter(infile,"num_generations_3",1);
	set_unsigned_int_parameter(infile,"num_generations_4",1);
	*/
	// the petsc stuff
	// different solver types. default(0), mumps_both(1), superlu_dist_both(2), mumps_1d_iter_3d(3), 
  std::string solver_options = set_string_parameter(infile,"solver_options","");
  output_folder << set_string_parameter(infile,"output_folder","results/");
  set_string_parameter(infile,"input_1d_file","");
  set_bool_parameter(infile,"_read_1d_mesh",false);
  set_unsigned_int_parameter(infile,"max_newton_iterations",6);
  set_bool_parameter(infile,"compliance_1d",false);
  set_bool_parameter(infile,"inertance_1d",false);
  set_bool_parameter(infile,"viscosity_controlled",false);
  set_bool_parameter(infile,"quadratic_ramp",false);
  set_double_parameter(infile,"ramp_duration",false);
	set_unsigned_int_parameter(infile,"geometry_type",1);	// 0 is pipe, 1 is bifurcating pipe
	set_unsigned_int_parameter(infile,"num_1d_trees",2);
	set_double_parameter(infile,"max_time_step",es->parameters.get<double>("dt")*1000.0);
	set_double_parameter(infile,"min_time_step",es->parameters.get<double>("dt")/10.0);
	set_bool_parameter(infile,"neumann_stabilised",true);
	set_bool_parameter(infile,"convective_form",false);
	set_double_parameter(infile,"backflow_stab_param",1.0);
	set_double_parameter(infile,"pressure_mag_3d",1.0);
	set_double_parameter(infile,"out_pressure_mag_3d",0.);
	set_double_parameter(infile,"time_scaling",1.0);	//this should be updated each timestep as a function of unsteady
	set_double_parameter(infile,"parent_pressure_mag",es->parameters.get<double>("pressure_mag_3d"));
	set_double_parameter(infile,"daughter_1_pressure_mag",0.0);
	set_double_parameter(infile,"daughter_2_pressure_mag",0.0);


	set_bool_parameter(infile,"neumann_stabilised_adjusted",false);
	set_bool_parameter(infile,"neumann_stabilised_adjusted_interpolated",false);
	set_bool_parameter(infile,"reynolds_number_calculation",false);
	double reynolds_number = set_double_parameter(infile,"reynolds_number",1.0);
	//double velocity_scale = set_double_parameter(infile,"velocity_scale",1.0);
	set_double_parameter(infile,"velocity_scale",1.0);
	//double length_scale = set_double_parameter(infile,"length_scale",1.0);
	set_double_parameter(infile,"length_scale",1.0);
	set_bool_parameter(infile,"constant_pressure",false);
	set_bool_parameter(infile,"compare_results",false);
	set_string_parameter(infile,"results_folder_1","results/");
	set_string_parameter(infile,"results_folder_2","results/");
	threed = set_bool_parameter(infile,"threed",true);
	set_unsigned_int_parameter(infile,"num_refinement_steps",0);
	set_unsigned_int_parameter(infile,"adaptive_newton_solve_limit",es->parameters.set<unsigned int> ("max_newton_iterations"));
	set_unsigned_int_parameter(infile,"adaptive_linear_iterations_limit",50);
	set_unsigned_int_parameter(infile,"adaptive_residual_iterations_limit",5);
	set_bool_parameter(infile,"repartition",false);
	set_bool_parameter(infile,"reuse_preconditioner",false);
	set_unsigned_int_parameter(infile,"error_estimator",0);
	set_bool_parameter(infile,"supg",false);
	set_bool_parameter(infile,"supg_newton",false);
	set_bool_parameter(infile,"supg_full_newton",false);
	set_bool_parameter(infile,"supg_picard",false);
	set_bool_parameter(infile,"pspg",false);
	set_bool_parameter(infile,"lsic",false);
	set_bool_parameter(infile,"pspg_newton",false);
	set_bool_parameter(infile,"pspg_picard",false);
	set_bool_parameter(infile,"optimisation_stabilised",false);
	set_bool_parameter(infile,"tangential_optimisation",false);
	set_double_parameter(infile,"supg_scale",1.0);
	set_bool_parameter(infile,"volume_optimisation",false);
	set_bool_parameter(infile,"volume_optimisation_h_scaled",false);
	set_double_parameter(infile,"optimisation_scaling_factor",1.0);
	set_bool_parameter(infile,"output_adjoint_variables",false);
	set_unsigned_int_parameter(infile,"inflow_profile",1);
	set_bool_parameter(infile,"reverse_inflow_profile",false);
	set_bool_parameter(infile,"quads",true);
	set_bool_parameter(infile,"fine_ref_solution",false);
	set_unsigned_int_parameter(infile,"supg_parameter",2);
	set_bool_parameter(infile,"linear_shape_functions",false);
	set_bool_parameter(infile,"supg_convection",false);
	set_bool_parameter(infile,"supg_convection_newton",false);
	set_bool_parameter(infile,"supg_test",false);
	set_bool_parameter(infile,"supg_constant_constant",false);
	set_bool_parameter(infile,"supg_laplacian",false);
	set_string_parameter(infile,"stokes_solver_options",es->parameters.get<std::string> ("solver_options"));
	set_bool_parameter(infile,"output_linear_iteration_count",false);
	set_unsigned_int_parameter(infile,"output_no_refinement",es->parameters.get<unsigned int> ("no_refinement"));	// the refinement level at which we want output

	set_double_parameter(infile,"density_1d",1.);
	set_double_parameter(infile,"viscosity_1d",1.);
	set_double_parameter(infile,"flow_mag_1d",1.);
	set_unsigned_int_parameter(infile,"resistance_type_1d",0);	//0- poiseuille, 1-pedley, 2-ertbruggen, 3- reynolds, 4 - hard coded
	set_double_parameter(infile,"length_diam_ratio",3.);
	set_unsigned_int_parameter(infile,"bc_type_1d",0);	//0- flow-pressure/acinar 1-pressure-pressure/acinar
	set_double_parameter(infile,"in_pressure_1d",0);	// pressure at inlet for pressure-pressure BCs
	set_double_parameter(infile,"out_pressure_1d",0);	// pressure at inlet for pressure-pressure BCs

	set_bool_parameter(infile,"prescribed_flow",false);

	set_unsigned_int_parameter(infile,"old_geometry",0.);
	set_bool_parameter(infile,"p2p1_stabilisation",false);
	set_bool_parameter(infile,"opt_correction_term",false);
	set_bool_parameter(infile,"discontinuous_linear_neumann",false);
	set_bool_parameter(infile,"expanding_pipe",false);
	set_double_parameter(infile,"numerical_continuation_scaling_factor",2.0);

	es->parameters.set<unsigned int> ("num_newton_steps") = 0;
	set_bool_parameter(infile,"minimise_neumann_boundary",true);

	set_bool_parameter(infile,"pressure_minimisation",false);
	set_bool_parameter(infile,"bertoglio_stabilisation",false);
	set_double_parameter(infile,"bertoglio_stab_param",1.0);
	
	std::cout << "bertoglio_stab_param = " << es->parameters.set<double> ("bertoglio_stab_param") << std::endl;
 	set_double_parameter(infile,"cube_length",1.0);
	set_double_parameter(infile,"cube_width",1.0);
	set_double_parameter(infile,"cube_height",1.0);
	set_unsigned_int_parameter(infile,"cube_length_N",5);
	set_unsigned_int_parameter(infile,"cube_width_N",5);
	set_unsigned_int_parameter(infile,"cube_height_N",5);
	set_bool_parameter(infile,"prescribed_womersley",false);
	set_string_parameter(infile,"prescribed_womersley_file","");
	set_bool_parameter(infile,"symmetric_gradient",false);
	set_double_parameter(infile,"radius",0.5);
	set_bool_parameter(infile,"neumann_stabilised_linear",false);
	set_bool_parameter(infile,"tangential_velocity_optimisation",false);
	set_bool_parameter(infile,"pearson_example",false);
	set_bool_parameter(infile,"flat_lid_driven_cavity",false);
	set_int_parameter(infile,"random_seed",0);

	set_double_parameter(infile,"mesh_input_scaling_1d",1.0);
	set_double_parameter(infile,"mesh_input_scaling_3d",1.0);
	set_double_parameter(infile,"initial_segment_length",0.04);

	set_double_parameter(infile,"flow_mag_3d",1.0);
	set_double_parameter(infile,"velocity_mag_3d",1.0);

	set_unsigned_int_parameter(infile,"preconditioner_type_3d",0);
	set_unsigned_int_parameter(infile,"preconditioner_type_3d1d",0);

	set_unsigned_int_parameter(infile,"alveolated_1d_tree",0);
	set_double_parameter(infile,"alveolar_length_diam_ratio",2.2);
	set_double_parameter(infile,"alveolar_diameter",0.0005);
	set_unsigned_int_parameter(infile,"num_alveolar_generations",0);
	set_double_parameter(infile,"hofmann_breath_hold",0.2);
	set_double_parameter(infile,"max_enhancement_factor",2.0);

	set_unsigned_int_parameter(infile,"impaction_type",0);

	set_bool_parameter(infile,"ksp_view",false);
	set_bool_parameter(infile,"ksp_view_before",false);
	set_bool_parameter(infile,"nonzero_initial_guess",false);
	set_bool_parameter(infile,"nonzero_initial_guess_inner_3d1d",false);
	set_bool_parameter(infile,"streamline_diffusion",false);
	set_bool_parameter(infile,"leaky_lid_driven_cavity",false);
	set_bool_parameter(infile,"pin_pressure",false);

	set_bool_parameter(infile,"fieldsplit",false);
	set_bool_parameter(infile,"direct",false);


	set_bool_parameter(infile,"mesh_dependent_stab_param",false);
	set_bool_parameter(infile,"gravemeier_element_length",false);
	set_unsigned_int_parameter(infile,"pcd_boundary_condition_type",0);

	set_bool_parameter(infile,"multiply_system_by_dt",false);
	set_double_parameter(infile,"numerical_continuation_starting_reynolds_number",1.0);

	set_double_parameter(infile,"element_length_scaling",1.0);

	set_bool_parameter(infile,"output_nondim",false);
	set_bool_parameter(infile,"half_initial_length",false);
	set_bool_parameter(infile,"twod_oned_tree",false);
	set_bool_parameter(infile,"initial_segment_length_from_mesh",false);
	

	es->parameters.set<double> ("last_nonlinear_iterate") = 1.0;
	
	set_int_parameter(infile,"inflow_bdy_id",-1);
	set_int_parameter(infile,"wall_bdy_id",12);	
	set_int_parameter(infile,"wall_bdy_id_2",1010);		


  set_string_parameter(infile,"input_1d_node_file_2","");
  set_string_parameter(infile,"input_1d_edge_file_2","");


  set_bool_parameter(infile,"match_1d_mesh_to_3d_mesh",false);
  set_bool_parameter(infile,"assume_symmetric_tree",false);
  set_bool_parameter(infile,"calculate_1d_info_at_coupling_nodes",false);
  set_bool_parameter(infile,"radius_on_edge",false);

  set_double_parameter(infile,"nonlinear_tolerance_1d",1e-4);
  set_bool_parameter(infile,"nonlinear_1d",false);

  set_bool_parameter(infile,"use_centreline_data",false);

	// 0 - 3D distributed everywhere, 0D distributed everywhere
	// 1 - 3D distributed on procs 0 - n-2, 0D on proc n-1
	// 2 - 3D distributed everywhere, 0D on proc 0
  set_unsigned_int_parameter(infile,"custom_partitioning",0);
  set_bool_parameter(infile,"output_resistance_for_each_tree",false);
  set_bool_parameter(infile,"output_resistance",true);

  set_bool_parameter(infile,"increasing_residuals_exit",true);
  set_unsigned_int_parameter(infile,"max_initial_picard",0);
  set_double_parameter(infile,"picard_to_newton_threshold",0.1);

  set_bool_parameter(infile,"shift_pressure_bc",0);

  set_bool_parameter(infile,"compute_eigenvalues",false);

  // ************** PETSC OPTIONS ********************************* //

	set_string_parameter(infile,"petsc_3d_fieldsplit_solver_options","");
	set_string_parameter(infile,"petsc_3d_direct_solver_options","");
	set_string_parameter(infile,"petsc_3d_gmres_solver_options","");


	set_string_parameter(infile,"petsc_convection_diffusion_solver_options","");
	set_string_parameter(infile,"petsc_navier_stokes_schur_solver_options","");

	set_string_parameter(infile,"petsc_mass_matrix_solver_options","");
	set_string_parameter(infile,"petsc_laplacian_solver_options","");

	set_string_parameter(infile,"petsc_1d_solver_options","");

	set_string_parameter(infile,"petsc_3d1d_direct_solver_options","");
	set_string_parameter(infile,"petsc_3d1d_diagonal_solver_options","");
	set_string_parameter(infile,"petsc_3d1d_schur_solver_options","");
	set_string_parameter(infile,"petsc_3d1d_gmres_solver_options","");

	set_string_parameter(infile,"petsc_3d1d_navier_stokes_direct_solver_options","");
	set_string_parameter(infile,"petsc_3d1d_navier_stokes_gmres_solver_options","");
	set_string_parameter(infile,"petsc_3d1d_navier_stokes_fieldsplit_solver_options","");
	set_string_parameter(infile,"petsc_3d1d_navier_stokes_fieldsplit_preonly_solver_options","");

	set_string_parameter(infile,"petsc_3d1d_inner_1d_solver_options","");

	set_string_parameter(infile,"petsc_3d1d_schur_stokes_direct_solver_options","");
	set_string_parameter(infile,"petsc_3d1d_schur_stokes_gmres_solver_options","");
	set_string_parameter(infile,"petsc_3d1d_schur_stokes_fieldsplit_solver_options","");

	set_string_parameter(infile,"petsc_3d1d_post_solve_fieldsplit_solver_options","");

	set_string_parameter(infile,"petsc_moghadam_velocity_gmres_solver_options","");
	set_string_parameter(infile,"petsc_moghadam_gmres_solver_options","");
	set_string_parameter(infile,"petsc_moghadam_schur_solver_options","");


	set_string_parameter(infile,"petsc_solver_options","");
	set_string_parameter(infile,"petsc_solver_options_ksp_view","");

  // *********************** //

	set_unsigned_int_parameter(infile,"random_1d_generations",0);
	set_string_parameter(infile,"num_generations_string","");
	// this one is used for all
	set_unsigned_int_parameter(infile,"num_generations",1);

	// 0 is just use the whole matrix, 1 is multiple column solve in libmesh, 2 is multiple column solve in petsc
	set_unsigned_int_parameter(infile,"multiple_column_solve",0);

	set_unsigned_int_parameter(infile,"moghadam_coupling",0);	// some files may think that this should be a bool

	set_bool_parameter(infile,"assemble_pressure_mass_matrix",false);
	set_bool_parameter(infile,"assemble_scaled_pressure_mass_matrix",false);
	set_bool_parameter(infile,"assemble_pressure_convection_diffusion_matrix",false);
	set_bool_parameter(infile,"assemble_pressure_laplacian_matrix",false);
	set_bool_parameter(infile,"assemble_velocity_mass_matrix",false);

	set_bool_parameter(infile,"schur_stokes_precompute",false);
	set_bool_parameter(infile,"monolithic_navier_stokes_preonly",false);
	set_double_parameter(infile,"monolithic_navier_stokes_preonly_switch",0.);
	set_unsigned_int_parameter(infile,"preconditioner_type_schur_stokes",0);
	set_unsigned_int_parameter(infile,"monolithic_navier_stokes_gmres_its",100);

	set_double_parameter(infile,"length_diam_ratio_1",0.);
	set_double_parameter(infile,"length_diam_ratio_2",0.);
	set_double_parameter(infile,"resistance_scale_1",0.);
	set_double_parameter(infile,"resistance_scale_2",0.);



	set_bool_parameter(infile,"post_solve",false);
	set_bool_parameter(infile,"direct_post_solve",false);
	set_bool_parameter(infile,"0D",false);

	set_bool_parameter(infile,"stokes_ic",false);
	set_double_parameter(infile,"angle_of_inflow",0.); // angle (in rad) by which inflow is rotated from normal

	set_double_parameter(infile,"length_diam_ratio_1",0.);
	set_double_parameter(infile,"length_diam_ratio_2",0.);

	set_bool_parameter(infile,"write_eigenvalues",false);

	set_bool_parameter(infile,"negative_mono_schur_complement",false);
	set_bool_parameter(infile,"negative_bfbt_schur_complement",false);

	// scale monolithic preconditioner by resistance
	// 0 - no scaling
	// 1 - scale by resistance of first pipe, max over all boundaries, min 1.0
	set_unsigned_int_parameter(infile,"scale_mono_preconditioner",0);
	set_double_parameter(infile,"mono_preconditioner_resistance_scaling",1.0);

	set_string_parameter(infile,"boundary_conditions_folder","");
	set_bool_parameter(infile,"known_boundary_conditions",false);

	set_bool_parameter(infile,"renumber_nodes_and_elements",true);
	set_double_parameter(infile,"matching_3d1d_tolerance",0.);

	set_bool_parameter(infile,"reuse_convection_diffusion_pc",true);

	set_bool_parameter(infile,"use_command_line_auto_fieldsplit_options",true);
	set_bool_parameter(infile,"output_backup_files",true);

	set_bool_parameter(infile,"efficient_assembly",true);
	set_bool_parameter(infile,"efficient_solve_setup",true);

	set_bool_parameter(infile,"create_3d_preconditioner_once",false);

	set_double_parameter(infile,"0d_bifurcation_angle",35* M_PI/180.);	// changed default to 35

	set_bool_parameter(infile,"multiple_output_files",true);

	set_bool_parameter(infile,"efficient_monolithic",true);

	set_bool_parameter(infile,"output_system_matrix",false);

	// put the characteristic length of the mesh in here (SI units)
	// - i.e. diameter of inflow surface boundary in 3D
	// - diameter of first element in 0D
	set_double_parameter(infile,"characteristic_length",1.);

	set_bool_parameter(infile,"unsteady_from_steady",false);

	set_double_parameter(infile,"read_time_step_size",0.1);
	set_bool_parameter(infile,"read_time_steps",false);
	set_double_parameter(infile,"read_end_time",0.);


	set_bool_parameter(infile,"gmsh_diff_wall_bdy_id",false);
	set_bool_parameter(infile,"bifurcation_start_1d",false);

	set_double_parameter(infile,"generation_ratio_1d",0.8);
	set_double_parameter(infile,"trachea_diameter",0.02);

	set_bool_parameter(infile,"subdomain_airways",false);

	set_bool_parameter(infile,"vary_outlet_pressure",false);

	set_bool_parameter(infile,"apply_pressure_from_flux",false);
	set_double_parameter(infile,"acceptable_flow_rate_error",0.01);
	set_bool_parameter(infile,"interpolate_pressure_from_previous_steps",false);

	// 0 - use resistance method, 1 - use newton's method
	set_unsigned_int_parameter(infile,"flow_rate_convergence_method",0);

	// number of symmetric generations to add to the geometry
	set_unsigned_int_parameter(infile,"additional_symmetric_generations",0);
	set_bool_parameter(infile,"alveolar_geometry",false);

	set_unsigned_int_parameter(infile,"nonlinear_iteration_1d",0);
	set_bool_parameter(infile,"partitioned_converge",true);

	// use the flow rate from the precious time step for the resistance
	set_bool_parameter(infile,"semi_implicit_1d",false);


	set_bool_parameter(infile,"construct_moghadam_preconditioner",false);

	set_bool_parameter(infile,"dimensional_calculation_fix",true);

	set_bool_parameter(infile,"newton_1d",false);

	// 0 - symmetric, 1 - different left and right, 2 - other random options
	set_unsigned_int_parameter(infile,"asymmetric_0d_tree",0);

	set_double_parameter(infile,"generation_ratio_1d_left",0.876);
	set_double_parameter(infile,"generation_ratio_1d_right",0.686);

	set_unsigned_int_parameter(infile,"max_nonlinear_iterations_1d",20);

	// moghadam_preconditioner will be preconditioner_3d == 13
	set_double_parameter(infile,"moghadam_outer_tol",0.4);
	set_double_parameter(infile,"moghadam_gmres_tol",0.01);
	set_double_parameter(infile,"moghadam_cg_tol",0.2);

	set_bool_parameter(infile,"moghadam_augment_sparsity_pattern",true);

	// apply left jacobi preconditioning within libmesh
	set_unsigned_int_parameter(infile,"libmesh_jacobi",0);

	// assemble residuals to form 
	// 0 is normal
	// 1 is solve normally, but allow for update of rhs based on new residual a la moghadam
	// 2 is solve using newton increments, allowing for update of rhs based on new residual
	set_bool_parameter(infile,"residual_formulation",false);
	set_bool_parameter(infile,"residual_formulation_0d",false);
	set_bool_parameter(infile,"assemble_residual_only",false);
	set_bool_parameter(infile,"residual_linear_solve",false);
	set_double_parameter(infile,"residual_linear_solver_tolerance",0.4);
	set_unsigned_int_parameter(infile,"residual_max_iterations",5);

	set_bool_parameter(infile,"output_then_exit",true);
	
	// 0 - normal moghadam preconditioner wPC, 1 - moghadam woPC, 2 - moghadam exact schur
	// only 0 needs construction of moghadam preconditioner
	set_unsigned_int_parameter(infile,"moghadam_preconditioner_type",0);

	set_bool_parameter(infile,"construct_moghadam_h_inverse",false);

	set_bool_parameter(infile,"moghadam_velocity_preconditioner",false);

	set_double_parameter(infile,"adaptive_time_step_reduce_factor",0.5);
	set_double_parameter(infile,"adaptive_time_step_increase_factor",1.1);

	set_bool_parameter(infile,"max_newton_iterations_exit",true);

	set_bool_parameter(infile,"set_old_solution",true);
	set_double_parameter(infile,"moghadam_linear_resistance_epsilon",1e-6);

	// 0 - is normal, no model, equiv to infinite compliance, 1 - simple compliance model with extra V_acin dof
	set_unsigned_int_parameter(infile,"acinar_model",0);
	set_double_parameter(infile,"total_acinar_compliance",0.2e-5);
	set_double_parameter(infile,"pl_mean",-1500.0);	// Pa
	set_double_parameter(infile,"pl_range",500.0);	// Pa
	set_double_parameter(infile,"pl_period",es->parameters.get<double>("period"));
	set_unsigned_int_parameter(infile,"pl_type",0);	// 0 - sin, 1 - cos, 2 - -cos, 3- constant, 4- cosine volume based
	set_double_parameter(infile,"pl_volume_mean",2.25e-3);	// Pa
	set_double_parameter(infile,"pl_volume_range",6e-5);	// Pa
	set_double_parameter(infile,"pl_total_resistance",2.1e3);	// Pa
	set_bool_parameter(infile,"init_acinar_volumes",true);	// set it to its steady state value -C*P_pl
	// hard coded resistances
	set_double_parameter(infile,"resistance_0",1.3e-2);
	set_double_parameter(infile,"resistance_0_a",1.3e-2);
	set_double_parameter(infile,"resistance_0_b",1.3e-2);
	set_double_parameter(infile,"resistance_1",1.);
	set_double_parameter(infile,"resistance_2",1.);
	set_double_parameter(infile,"resistance_3",1.);
	set_double_parameter(infile,"resistance_4",1.);
	set_double_parameter(infile,"resistance_5",1.);

	set_bool_parameter(infile,"nonlinear_resistance_type_1d",false);

	set_bool_parameter(infile,"increment_boundary_conditions",false);
	set_double_parameter(infile,"current_reynolds_number",0.);
	set_double_parameter(infile,"max_current_reynolds_number",0.);

	set_bool_parameter(infile,"nondimensionalised",false);

	set_bool_parameter(infile,"mono_flow_rate_penalty",false);
	set_double_parameter(infile,"mono_flow_rate_penalty_param",1e2);
	set_double_parameter(infile,"mono_flow_rate_penalty_param_2",1e4);
	es->parameters.set<double> ("mono_flow_rate_penalty_param") = es->parameters.get<double> ("mono_flow_rate_penalty_param_2");	// don't ask.................

	set_bool_parameter(infile,"twod_oned_tree_pressure",false);

	set_double_parameter(infile,"length_adjustment_factor_0d",0.99);

	set_double_parameter(infile,"start_radius_3d",0.015);
	set_unsigned_int_parameter(infile,"num_generations_3d",1); // essentially, number of bifurcations.
	set_bool_parameter(infile,"resistance_adjusted_radius",false);

	set_bool_parameter(infile,"first_gen_poiseuille",false);

	// 0 - just static pressure (as before), 
	// 1 - explicit (prev timestep) adjusted
	// 2 - semi-implicit (prev nonlin step) adjusted
	// 3 - semi-implicit (prev nonlin step) adjusted half way
	// 4 - implicit (newton)
	set_unsigned_int_parameter(infile,"total_pressure_bc",0);
	set_bool_parameter(infile,"output_total_pressure",false);

  restart_folder << set_string_parameter(infile,"restart_folder",output_folder.str());




	// ********* RESOLVE CONFLICTS IN PARAMETERS ETC ************** //




	// ****************** can only have moghadam coupling with sim_type 2 ************************** //
	if(sim_type != 2 && es->parameters.get<unsigned int> ("moghadam_coupling"))
	{
		std::cout << "- Can ony have moghadam_coupling with sim_type 2." << std::endl;
		std::cout << "Changing to use sim_type == 2." << std::endl;
		es->parameters.set<unsigned int> ("sim_type") = 2;
		sim_type = 2;
		
	}

	// ***************** can only have augment moghadam sparsity pattern if donig moghadam coupling ********** //
	if(!es->parameters.get<unsigned int> ("moghadam_coupling") && es->parameters.get<bool> ("moghadam_augment_sparsity_pattern"))
	{
		std::cout << "- No need to augment the sparsity pattern for moghadam if not using moghadam coupling" << std::endl;	
		std::cout << "Changing to moghadam_augment_sparsity_pattern == false" << std::endl;
		es->parameters.set<bool> ("moghadam_augment_sparsity_pattern") = false;
	}



	// ******************** now set some derived variables ******************** //
	// 0 - only 3d, 1 - only 1d, 2 - 3d and 1d independent, 3 - implicit coupling
	// 4 - explicit coupling (timestep), 5 - monolithic
	if(sim_type == 0 || sim_type >= 2)
		sim_3d = true;

	if(sim_type >= 1)
		sim_1d = true;

	// set the restart parameter
	// if we are doing unsteady from steady, that means restart
	if(restart_time_step > 0 || es->parameters.get<bool> ("unsteady_from_steady"))
		restart = true;

	t_step = restart_time_step;
	es->parameters.set<unsigned int> ("t_step")   = t_step;




	

	// ********************* RANDOM SEED SETUP ******************************** //
	if(es->parameters.get<int> ("random_seed") == 0)
	{
		es->parameters.set<int> ("random_seed") = std::time(NULL);
		std::cout << "- Random seed not given so using std::time(NULL)." << std::endl;
		std::cout << "\trandom_seed = " << es->parameters.get<int> ("random_seed") << std::endl;
	}
	srand( es->parameters.get<int> ("random_seed"));




	// ***************** if we aren't doing a monolithic simulation then the 3d1d preconditioner is irrelevant ********** //
	if(sim_type != 5)
	{
		es->parameters.set<unsigned int> ("preconditioner_type_3d1d") = 0;
	}	

	// ***************** if we aren't doing a 3D simulations then the 3d preconditioner is irrelevant ******************* //
	if(!sim_3d)
	{
		es->parameters.set<unsigned int> ("preconditioner_type_3d") = 0;
	}


	// ******************** SETUP PETSC OPTIONS ******************************* //
	// NOTE: need to add some options here for moghadam preconditioner

  const std::string empty_string;
	if(es->parameters.get<std::string> ("petsc_solver_options") == empty_string)
	{
		std::string petsc_solver_options = "";
		std::string petsc_solver_options_ksp_view = "";
		if(sim_3d)
		{

			// figure out if we are doing monolithic simulation or not
			if(sim_type != 5)
			{
				// we need to figure out what 3d preconditioner we are using
				if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 0)
				{
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d_direct_solver_options");
					es->parameters.set<bool> ("direct") = true;
				}
				else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 1)
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d_gmres_solver_options");
				else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13)
				{
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_moghadam_velocity_gmres_solver_options");
			
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_moghadam_schur_solver_options");

					// add the other options later
				}
				else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 14)
				{
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_moghadam_gmres_solver_options");
			
					// add the other options later
				}
				else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") > 1)
				{
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d_fieldsplit_solver_options");


					// extra navier stokes settings
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_convection_diffusion_solver_options");
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_navier_stokes_schur_solver_options");	
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_mass_matrix_solver_options");		
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_laplacian_solver_options");
					es->parameters.set<bool> ("fieldsplit") = true;
				}

				

			}
			else
			{
				// we need to figure out what 3d1d monolithic preconditioner we are using
				if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 0)
				{
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_direct_solver_options");
					es->parameters.set<bool> ("direct") = true;
				}
				else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 1)
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_gmres_solver_options");
				else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 6)
				{
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_diagonal_solver_options");
					es->parameters.set<bool> ("fieldsplit") = true;
				}
				else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") > 6)
				{
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_schur_solver_options");
					es->parameters.set<bool> ("fieldsplit") = true;
				}


				// now we need to add the navier-stokes specific preconditioner details if doing a fieldsplit monolithic solve
				if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") >= 6)
				{
					if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 0)
						petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_navier_stokes_direct_solver_options");
					else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 1)
						petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_navier_stokes_gmres_solver_options");
					else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") > 1)
					{
						if(es->parameters.get<bool> ("monolithic_navier_stokes_preonly"))
							petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_navier_stokes_fieldsplit_preonly_solver_options");
						else
							petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_navier_stokes_fieldsplit_solver_options");

						// extra navier stokes settings
						petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_convection_diffusion_solver_options");
						petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_navier_stokes_schur_solver_options");	
						petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_mass_matrix_solver_options");		
						petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_laplacian_solver_options");
						es->parameters.set<bool> ("fieldsplit") = true;
					}
				}

				// we need to add the schur stokes solver options
				if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 10)
				{
					if(es->parameters.get<unsigned int> ("preconditioner_type_schur_stokes") == 0)
						petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_schur_stokes_direct_solver_options");
					else if(es->parameters.get<unsigned int> ("preconditioner_type_schur_stokes") == 1)
						petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_schur_stokes_gmres_solver_options");
					else if(es->parameters.get<unsigned int> ("preconditioner_type_schur_stokes") == 2)
						petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_schur_stokes_fieldsplit_solver_options");
					else
					{
						std::cout << "- Schur Stokes solver type " << es->parameters.get<unsigned int> ("preconditioner_type_schur_stokes") << " not defined." << std::endl;
						std::cout << "Exiting..." << std::endl;
						std::exit(0);
					}
					
					
				}

				if(es->parameters.get<bool> ("post_solve"))
				{
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_post_solve_fieldsplit_solver_options");
				}


				
			}



		}

		if(sim_1d)
		{

			// need to figure out if we have a monolithic method
			// for various exclusive 1d computations
			petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_1d_solver_options");

			if(sim_type == 5)
			{
				// need to figure out if we have a separate solver for the 1d part of the simulation
				if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") >= 6)
					petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_3d1d_inner_1d_solver_options");
					
			}
		}


		petsc_solver_options_ksp_view = petsc_solver_options;
		if(es->parameters.get<bool> ("ksp_view"))
		{
			if(sim_3d)
			{
				es->parameters.set<bool> ("ksp_view_3d") = true;
				if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13)
				{
					petsc_solver_options_ksp_view += " -moghadam_velocity_ksp_view";
					petsc_solver_options_ksp_view += " -moghadam_schur_ksp_view";
				}
				else if(sim_type != 5)
					petsc_solver_options_ksp_view += " -ns3d_ksp_view";
				else
					petsc_solver_options_ksp_view += " -ns3d1d_ksp_view";
			}
		
			if(sim_1d)
			{
				es->parameters.set<bool> ("ksp_view_1d") = true;
				if(sim_type != 5)
					petsc_solver_options_ksp_view += " -ns1d_ksp_view";
				
			}
		}
			

		es->parameters.set<std::string> ("petsc_solver_options") = petsc_solver_options;
		es->parameters.set<std::string> ("petsc_solver_options_ksp_view") = petsc_solver_options_ksp_view;
		std::cout << "- Petsc command line options set to: " << es->parameters.get<std::string> ("petsc_solver_options") << std::endl;;
	}

	// set petsc options from file, petsc command line options will overwrite
	PetscOptionsInsertString(es->parameters.get<std::string> ("petsc_solver_options").c_str());

	if(es->parameters.get<bool> ("direct"))
	{
		if(es->parameters.get<bool> ("fieldsplit"))
		{
			std::cout << "- Using direct solver so cannot use fieldsplit. Setting to not use fieldsplit." << std::endl;	
			std::cout << "\tfieldsplit = false" << std::endl;
			es->parameters.set<bool> ("fieldsplit") = false;
		}

		if(es->parameters.get<bool> ("nonzero_initial_guess"))
		{
			std::cout << "- Using direct solver so makes no sense to use nonzero_initial_guess. Setting to not use nonzero_initial_guess." << std::endl;
			std::cout << "\tnonzero_initial_guess = false" << std::endl;
			es->parameters.set<bool> ("nonzero_initial_guess") = false;
		}

		if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") || es->parameters.get<unsigned int> ("preconditioner_type_3d"))
		{
			std::cout << "- Using direct solver so no preconditioner used. Setting to not use preconditioner." << std::endl;
			std::cout << "\tpreconditioner_type_3d = 0" << std::endl; 
			std::cout << "\tpreconditioner_type_3d1d = 0" << std::endl; 
			es->parameters.set<unsigned int> ("preconditioner_type_3d") = 0;
			es->parameters.set<unsigned int> ("preconditioner_type_3d1d") = 0;
		}
	}



	// ***************** if we aren't doing schur stokes then we can't do schur stokes precompute *********** //
	if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") != 10)
	{
		std::cout << "- No point doing schur stokes precompute if we aren't doing schur stokes." << std::endl;
		es->parameters.set<bool> ("schur_stokes_precompute") = false;
		
	}



	// **************** if we aren't doing a schur full matrix computation, no need for multiple column solve
	if(es->parameters.get<unsigned int> ("multiple_column_solve") 
		&& !(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 8
		|| es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 9
		|| es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 10
		|| es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 11
		|| es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 13))
	{
		std::cout << "- No need to do a multiple column solve if not doing preconditioner_type_3d1d 8,9,10,11." << std::endl;
		std::cout << "Changing to not use multiple_column_solve." << std::endl;
		es->parameters.set<unsigned int> ("multiple_column_solve") = 0;
	}

	if(es->parameters.get<unsigned int> ("multiple_column_solve") && es->parameters.get<unsigned int> ("custom_partitioning") != 1)
	{
		std::cout << "- Must use custom partitioning = 1 when doing multiple column solve for now." << std::endl;
		std::cout << "Changing to use custom partitioning = 1." << std::endl;
		es->parameters.set<unsigned int> ("custom_partitioning") = 1;
	}



	// ***************** figure out what sub matrices we need to assemble ********************* //
	// defaults
	es->parameters.set<bool> ("assemble_scaled_pressure_mass_matrix") = false;
	es->parameters.set<bool> ("assemble_scaled_pressure_mass_matrix") = false;
	es->parameters.set<bool> ("assemble_pressure_mass_matrix") = false;
	es->parameters.set<bool> ("assemble_pressure_convection_diffusion_matrix") = false;
	es->parameters.set<bool> ("assemble_pressure_laplacian_matrix") = false;
	es->parameters.set<bool> ("assemble_velocity_mass_matrix") = false;
	if(sim_3d)
	{

		if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 1)
			es->parameters.set<bool> ("assemble_scaled_pressure_mass_matrix") = true;

		if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 3
				|| es->parameters.get<unsigned int> ("preconditioner_type_schur_stokes") == 2)
			es->parameters.set<bool> ("assemble_scaled_pressure_mass_matrix") = true;

		if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 4 
				|| es->parameters.get<unsigned int> ("preconditioner_type_3d") == 5)
		{
			es->parameters.set<bool> ("assemble_pressure_mass_matrix") = true;
			es->parameters.set<bool> ("assemble_pressure_convection_diffusion_matrix") = true;
			es->parameters.set<bool> ("assemble_pressure_laplacian_matrix") = true;
		}
	
		if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 7
				|| es->parameters.get<unsigned int> ("preconditioner_type_3d") == 8
				|| es->parameters.get<unsigned int> ("preconditioner_type_3d") == 9)
			es->parameters.set<bool> ("assemble_velocity_mass_matrix") = true;
	}

	// if a 3d1d sim
	// defaults
	es->parameters.set<bool> ("assemble_velocity_matrix") = false;
	
	if(sim_type >= 3)
	{
		if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 8
				|| es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 9
				|| es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 10
				|| es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 11)
		{
			es->parameters.set<bool> ("assemble_velocity_matrix") = true;
			if(es->parameters.get<unsigned int> ("preconditioner_type_schur_stokes") == 2)
				es->parameters.set<bool> ("assemble_velocity_mass_matrix") = true;
				
		}

	}
	
		


	// *************** allocate the monolithic monitor context ***************** //
	PetscNew(&mono_ctx);

	PetscNew(&initial_guess_ctx);





	// **************** quadratic elements and schur stokes fieldsplit not supported ********************* //
	if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 10 && es->parameters.get<unsigned int> ("preconditioner_type_schur_stokes") == 2 && !es->parameters.get<bool> ("linear_shape_functions"))
	{
		std::cout << "- Sorry, fieldsplit and schur stokes are not supported with quadratic elements." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}


	// *************** newton 2,3 (newton (no rxn term) and semi-implicit) with optimised stabilised boundary conditions not implemented ********** //
	if((es->parameters.get<unsigned int> ("newton") == 2 || es->parameters.get<unsigned int> ("newton") == 3) 
			&& es->parameters.get<bool> ("optimisation_stabilised"))
	{
		std::cout << "- Sorry, newton (no rxn) and semi-implicit are not supported with optimised stabilised boundary conditions." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}


	// *************** can't use semi-implicit with steady **************************************** //

	if(es->parameters.get<unsigned int> ("newton") == 3 && !unsteady)
	{
		std::cout << "- Can't use a semi-implicit method with a steady simulation." << std::endl;
		std::cout << "Current setting: newton - " << es->parameters.get<unsigned int> ("newton") << std::endl;
		std::cout << "Changing to implicit Newton (with reaction terms)" << std::endl;

		es->parameters.set<unsigned int> ("newton") = 2;

		std::cout << "New setting: newton - " <<  es->parameters.get<unsigned int> ("newton") << std::endl;
		
	}


	// *************** make sure nonlinear scheme is backwards compatable ******************************* //
	
	// to be able to handle when the newton parameter was true or false.
	// i suppose if newton is true (i.e. intended to be implicit newton) it might be a large number
	// currently we have 4 nonlinear options
	// 0 - implicit picard
	// 1 - implicit newton
	// 2 - implicit newton (no rxn term)
	// 3 - semi-implicit newton (no rxn term)

	
	if(es->parameters.get<unsigned int> ("newton") > 3)
	{
		es->parameters.set<unsigned int> ("newton") = 1;
	}



	// ************** conservative form of navier stokes only implemented for newton - 0,1 ************** //
	if(!es->parameters.get<bool> ("convective_form") && es->parameters.get<unsigned int> ("newton") > 1) 
	{
		std::cout << "- Sorry, conservative Navier-Stokes formulation only supported for implicit Picard and Newton (with reaction terms)." << std::endl;
		std::cout << "Current setting: newton - " << es->parameters.get<unsigned int> ("newton") << std::endl;
		std::cout << "Changing to implicit Newton (no reaction terms)" << std::endl;

		es->parameters.set<unsigned int> ("newton") = 1;

		std::cout << "New setting: newton - " <<  es->parameters.get<unsigned int> ("newton") << std::endl;
    
	}


	// ***************************************** //
	// ********* nonlinear scheme chosen ******* //

	// ***************** if semi-implicit then we only need one nonlinear iteration ******************* //
	if(es->parameters.get<unsigned int> ("newton") == 3) 
	{
		es->parameters.set<unsigned int> ("max_newton_iterations") = 1;
	}



	// ***************** streamline diffusion not implemented for conservative form ******************* //
	if(!es->parameters.get<bool> ("convective_form") && es->parameters.get<bool> ("streamline_diffusion")) 
	{
		std::cout << "- Sorry, streamline diffusion not supported with conservative Navier-Stokes formulation." << std::endl;
		std::cout << "Changing to no streamline diffusion." << std::endl;

		es->parameters.set<bool> ("streamline_diffusion") = false;

		std::cout << "New setting: streamline_diffusion - " <<  es->parameters.get<bool> ("streamline_diffusion") << std::endl;
    
	}




	// ***************** if semi-implicit then have to do have a constant supg/pspg/lsic constant and have to have to use supg_picard with U_old ********************* //
	if(es->parameters.get<unsigned int> ("newton") == 3 && 
			(es->parameters.get<bool> ("supg") || es->parameters.get<bool> ("pspg") || es->parameters.get<bool> ("lsic")) ) 
	{

		if(!es->parameters.set<bool> ("supg_constant_constant"))
		{
			std::cout << "- Using semi-implicit method with supg/pspg/lsic. Must use have a constant supg constant." << std::endl;
			std::cout << "Changing to supg_constant_constant - true." << std::endl;

			es->parameters.set<bool> ("supg_constant_constant") = true;

			std::cout << "New setting: supg_constant_constant - " <<  es->parameters.get<bool> ("supg_constant_constant") << std::endl;
		}
	}

	// ***************** if semi-implicit then have to do supg_picard with U_old ********************* //
	if(es->parameters.get<unsigned int> ("newton") == 3 && es->parameters.get<bool> ("supg")) 
	{
		if(!es->parameters.get<bool> ("supg_picard"))
		{
			std::cout << "- Using semi-implicit method with supg. Must use have a the picard supg method (while using the velocity from the previous timestep)." << std::endl;
			std::cout << "Changing to supg_picard - true." << std::endl;

			es->parameters.set<bool> ("supg_picard") = true;
			es->parameters.set<bool> ("supg_newton") = false;
			es->parameters.set<bool> ("supg_full_newton") = false;
			es->parameters.set<bool> ("supg_convection_newton") = false;
			es->parameters.set<bool> ("supg_convection") = false;
				

			std::cout << "New setting: supg_picard - " <<  es->parameters.get<bool> ("supg_picard") << std::endl;
		}

    
	}

	// ***************** if semi-implicit then have to do pspg_picard with U_old ********************* //
	if(es->parameters.get<unsigned int> ("newton") == 3 && es->parameters.get<bool> ("pspg")) 
	{
		if(!es->parameters.get<bool> ("pspg_picard"))
		{
			std::cout << "- Using semi-implicit method with pspg. Must use have a the picard pspg method (while using the velocity from the previous timestep)." << std::endl;
			std::cout << "Changing to pspg_picard - true." << std::endl;

			es->parameters.set<bool> ("pspg_picard") = true;
			es->parameters.set<bool> ("pspg_newton") = false;
				

			std::cout << "New setting: pspg_picard - " <<  es->parameters.get<bool> ("pspg_picard") << std::endl;
		}
	}

	// ********* if semi-implicit and using neumann boundary stabilisation then need to use "linear" form ********* //
	if(es->parameters.get<unsigned int> ("newton") == 3 && es->parameters.get<bool> ("neumann_stabilised") && !es->parameters.get<bool> ("neumann_stabilised_linear")) 
	{
		std::cout << "- Using semi-implicit method with neumann stabilisation. Must linear stabilisation." << std::endl;
		std::cout << "Changing to neumann_stabilised_linear - true." << std::endl;

		es->parameters.set<bool> ("neumann_stabilised_linear") = true;
			

		std::cout << "New setting: neumann_stabilised_linear - " <<  es->parameters.get<bool> ("neumann_stabilised_linear") << std::endl;

	}



	// ***************** optimised stabilised bcs only supposed to be used with picard ******** //
	if(false)//es->parameters.get<bool> ("optimisation_stabilised"))
	{
		std::cout << "- Currently only picard iterations are supported for optimisation stabilised BCs. Setting to use picard iterations." << std::endl;
		std::cout << "\tnewton = false" << std::endl;
		es->parameters.set<unsigned int> ("newton") = 0;
	}
	else
	{
		//es->parameters.set<bool> ("output_adjoint_variables") = false;
	}
 




	if(es->parameters.get<bool> ("negative_bfbt_schur_complement") && es->parameters.get<unsigned int> ("preconditioner_type_3d") != 9)
	{
		std::cout << "- Negative BFBt schur complement currently only implemented for scaled stabilsated bfbt. Changning to positive BFBt schur complement." << std::endl;
		es->parameters.set<bool> ("negative_bfbt_schur_complement") = false;
		
	}


	// ***************** make sure that for sim_type 2 everything is correctly specified ********** //
	if(sim_type == 2 && es->parameters.get<bool> ("known_boundary_conditions"))
	{
		// check that folder has been specified
		if(es->parameters.get<std::string> ("boundary_conditions_folder").empty())
		{
			std::cout << "- Folder for getting boundary conditions for partitioned simulation not specified." << std::endl;
			std::cout << "Exiting..." << std::endl;
			std::exit(0);
		}

		// cannot have moghadam coupling and known boundary conditions
		if(es->parameters.get<unsigned int> ("moghadam_coupling"))
		{
			std::cout << "- Cannot have known_boundary_conditions and moghadam_coupling both specified." << std::endl;
			std::cout << "Exiting..." << std::endl;
			std::exit(0);
		}
	}


	// *********** adaptive refinement options ********** //
	//if not a 3D simulation then we cannot have adaptive refinement
	if(!sim_3d && adaptive_refinement)
	{
		std::cout << "- Can't have adaptive refinement if not simulating 3D. Setting to not use adaptive refinement." << std::endl;	
		std::cout << "\tadaptive_refinement = false" << std::endl;
		es->parameters.set<bool>("adaptive_refinement") = false;
		adaptive_refinement = false;
	}

	// if we have adaptive refinement then we def need an error estimator
	if(es->parameters.get<bool> ("adaptive_refinement") && !es->parameters.get<unsigned int> ("error_estimator") )
	{
		std::cout << "- Must have error estimator set when using adaptive refinement. Set to peclet number calculation." << std::endl;
		std::cout << "\terror_estimator = 3" << std::endl;
		es->parameters.set<unsigned int> ("error_estimator") = 3;
	}






	// ********** MAKE SURE STABILISATION OPTIONS ARE CONSISTENT ******************* //

	//if((!es->parameters.get<bool> ("stab") || !es->parameters.get<bool> ("threed")) && es->parameters.get<bool> ("ismail_supg_parameter"))
	if((!es->parameters.get<bool> ("linear_shape_functions")) && es->parameters.get<unsigned int> ("supg_parameter") == 0)
	{
		std::cout << "- Cannot use Ismail stabilisation parameters when not using linear elements. Set to use tezdeuyar 2 parameters." << std::endl;
		std::cout << "\tsupg_parameter = 2" << std::endl;
		es->parameters.set<unsigned int> ("supg_parameter") = 2;
	}

	if(!es->parameters.get<bool> ("pspg") && (!es->parameters.get<bool> ("stab") && es->parameters.get<bool> ("linear_shape_functions")))
	{
		std::cout << "- When not using pspg stabilisation, must use penalty stabilisation. Setting to use stabilisation with default parameter." << std::endl;
		std::cout << "\tstab = true" << std::endl;
		std::cout << "\talpha = 0.025" << std::endl;
		es->parameters.set<bool> ("stab") = true;
		es->parameters.set<double> ("alpha") = 0.025;
	}

	if((es->parameters.get<bool> ("lsic") || es->parameters.get<bool> ("pspg")) && !es->parameters.get<bool> ("linear_shape_functions"))
	{
		std::cout << "- Must use linear shape functions for lsic and/or pspg methinks. Setting to use linear shape functions." << std::endl;
		std::cout << "\tlinear_shape_functions = true" << std::endl;
		es->parameters.set<bool> ("linear_shape_functions") = true;
	}

	if(es->parameters.get<bool> ("stab") && !es->parameters.get<bool> ("linear_shape_functions"))
	{
		std::cout << "- Must use linear shape functions for penalty stabilised. Setting to use linear shape functions." << std::endl;
		std::cout << "\tlinear_shape_functions = true" << std::endl;
		es->parameters.set<bool> ("linear_shape_functions") = true;
	}

	if(es->parameters.get<bool>("stokes"))
	{
		if(es->parameters.get<bool>("supg"))
		{
			std::cout << "- Can't use SUPG when doing stokes. Setting to not use SUPG."<< std::endl;
			std::cout << "\tsupg = false" << std::endl;
			es->parameters.set<bool>("supg") = false;
		}

		if(es->parameters.get<bool>("pspg"))
		{
			std::cout << "- Can't use PSPG when doing stokes. Setting to not use PSPG."<< std::endl;
			std::cout << "\tpspg = false" << std::endl;
			es->parameters.set<bool>("pspg") = false;
		}

		if(es->parameters.get<bool>("lsic"))
		{
			std::cout << "- Can't use LSIC when doing stokes. Setting to not use LSIC."<< std::endl;
			std::cout << "\tlsic = false" << std::endl;
			es->parameters.set<bool>("lsic") = false;
		}
	}








	// ****************** SET THE REYNOLDS NUMBER **************************** //
	/*
	if(es->parameters.get<bool>("reynolds_number_calculation"))
	{
		es->parameters.set<double> ("length_scale") = 1.0;
		es->parameters.set<double> ("velocity_scale") = 1.0;
		es->parameters.set<double> ("viscosity") = 1.0/reynolds_number;
		es->parameters.set<double> ("density") = 1.0;
		es->parameters.set<double> ("in_pressure_1d") /= es->parameters.get<double> ("density") 
									* pow(es->parameters.get<double> ("velocity_scale"),2.);
	}
	else
	{
		es->parameters.set<double>("reynolds_number") = 
			es->parameters.get<double> ("density") * es->parameters.get<double> ("velocity_scale")
			 * es->parameters.set<double> ("length_scale") / es->parameters.set<double> ("viscosity");

		reynolds_number = es->parameters.get<double>("reynolds_number");
	}
	*/

	if(es->parameters.get<bool>("nondimensionalised"))
	{
		es->parameters.set<double>("reynolds_number") = 
			es->parameters.get<double> ("density") * es->parameters.get<double> ("velocity_scale")
			 * es->parameters.set<double> ("length_scale") / es->parameters.set<double> ("viscosity");
		reynolds_number = es->parameters.get<double>("reynolds_number");
	}
	else
	{
		es->parameters.set<double> ("length_scale") = 1.0;
		es->parameters.set<double> ("velocity_scale") = 1.0;
	}

	// we want all the time stuff in the program to be dimenionless
	time_scale_factor = es->parameters.get<double>("length_scale") /
							es->parameters.get<double>("velocity_scale");
	es->parameters.set<double>("time_scale_factor") = time_scale_factor;
	//std::cout << "time_scale_factor = " << time_scale_factor << std::endl;
	
	// don't change all the time stuff in the program
	/*
	es->parameters.set<double>("max_time_step") /= time_scale_factor;
	es->parameters.set<double>("min_time_step") /= time_scale_factor;
	es->parameters.set<double>("end_time") /= time_scale_factor;
	es->parameters.set<double>("restart_time") /= time_scale_factor;
	es->parameters.set<double>("period") /= time_scale_factor;
	dt /= time_scale_factor;
	es->parameters.set<double>("dt") /= time_scale_factor;
	*/









	// **************** SET PARAMETERS FOR STEADY STATE CALCULATION ************** //
	if(!unsteady)
	{
		dt = 1.0 / time_scale_factor;
		es->parameters.set<double>("dt") = dt;
		es->parameters.set<double>("end_time") = 1.0 / time_scale_factor;
	}
	es->parameters.set<double>("previous_dt") = dt;






	// ************* PIPE GEOMETRY CAN ONLY DO 1 TREE ************************* //
	if(es->parameters.get<unsigned int> ("geometry_type") == 0
			|| es->parameters.get<unsigned int> ("geometry_type") == 2)
	{
		std::cout << "- If we are doing a pipe geometry then we can only ever have 1 1D tree. Setting to use 1 1D tree." << std::endl;
		std::cout << "\tnum_1d_trees = 1" << std::endl;
		es->parameters.set<unsigned int> ("num_1d_trees") = 1;
	}


	if(es->parameters.get<unsigned int> ("geometry_type") == 5)
	{
		std::cout << "- Axisymmetric geometry is currently not working :( !!!... exiting..." << std::endl;
		std::exit(0);
	}



	// ************** ONLY DO GRAVEMEIER ELEMENT LENGTH FOR 3D **************** //
	if(es->parameters.get<bool> ("threed") && es->parameters.get<bool> ("gravemeier_element_length"))
	{
		std::cout << "- Can only use gravemeier element length in 3D. Setting to not use gravemeier element length." << std::endl;
		std::cout << "\tgravemeier_element_length = false" << std::endl;
		es->parameters.set<bool> ("gravemeier_element_length") = false;
		
	}


	



	
	// ****************** WE DON'T DO 1D ELEMENTS AT THE MOMENT *********** //
	if(es->parameters.get<bool>("0D") == false)
	{
		std::cout << "- Currently only 0D elements are supported. Setting to use 0D elements." << std::endl;
		std::cout << "\t0D = true" << std::endl;
  	es->parameters.set<bool>("0D") = true;
	}

	if(es->parameters.get<unsigned int>("resistance_type_1d") > 0 && es->parameters.get<unsigned int>("resistance_type_1d") < 4)
	{
		es->parameters.set<bool>("nonlinear_1d") = true;
	}

	if(false)//(es->parameters.get<unsigned int>("resistance_type_1d") == 0 ||
			//es->parameters.get<unsigned int>("resistance_type_1d") == 4 ||
			//es->parameters.get<unsigned int>("resistance_type_1d") == 5) && es->parameters.get<bool>("newton_1d"))
	{
		std::cout << "- No need to use newton on 0D when the resistance is linear (i.e. resistance_type_1d == 0,4,5)" << std::endl;
		std::cout << "changing to newton_1d = false" << std::endl;
		es->parameters.set<bool>("newton_1d") = false;
	}


	// **************** Stokes IC and coupled not implemented ************** //
	if(es->parameters.get<bool>("stokes_ic") && es->parameters.get<unsigned int>("sim_type") > 1)
	{
		std::cout << "- Stokes IC with coupled simulation or 0D simulation not currently supported." << std::endl;
		std::cout << "sim_type - " << es->parameters.get<unsigned int>("sim_type") << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	
	}

	// **************** if using stokes ic then we need to make sure it is done first ********* //
	if(es->parameters.get<bool>("stokes_ic"))
	{
		ic_set = false;
	}



	// **************** can only change the angle of inflow for 2D problems ******************* //
	if(es->parameters.get<bool>("threed") && fabs(es->parameters.get<double>("angle_of_inflow")) > 1e-10)
	{
		std::cout << "- Changing the angle of inflow for 3D problems not currently implemented." << std::endl;
		std::cout << "Current setting: angle_of_inflow - " << es->parameters.get<double>("angle_of_inflow") << std::endl;
		std::cout << "Changing to not alter angle of inflow." << std::endl;
		
		es->parameters.set<double>("angle_of_inflow") = 0.;

		std::cout << "New setting: angle_of_inflow - " << es->parameters.get<double>("angle_of_inflow") << std::endl;
		
	}



	// ************** IF WRITING EIGENVALUES MUST BE CALCULATING THEM ***************** //
	if(es->parameters.get<bool> ("write_eigenvalues") && !es->parameters.get<bool> ("compute_eigenvalues"))
	{
		std::cout << "- If we are writing eigenvalues to file we must calculate them too." << std::endl;
		std::cout << "Setting compute_eigenvalues to true." << std::endl;
		es->parameters.set<bool> ("compute_eigenvalues") = true;

	}

	// *************** SET DEFAULT n_initial_3d_elem ******************** //
	es->parameters.set<unsigned int>("n_initial_3d_elem") = 0;



	// *************** SET AUTO FIELDSPLIT OPTIONS ********************** //
	// doesn't work.
	//set_auto_fieldsplit_parameters();


	// *************** MUST USE CUSTOM PARTITIONING FOR MONOLITHIC ************ //
	// don't know why... lol doesn't set the boundary conditions properly
	if(false)//sim_type == 5 && es->parameters.get<unsigned int> ("custom_partitioning") < 1)
	{
		std::cout << "- Must use custom partitioning = 1 when using monolithi for some reason..." << std::endl;
		std::cout << "Setting custom_partitioning = 1." << std::endl;
		es->parameters.set<unsigned int> ("custom_partitioning") = 1;

	}


	// ************** PARAMETER TO TELL WHEN STOKES VELOCITY MATRIX HAS BEEN ALLOCATED ******** //
	es->parameters.set<bool> ("stokes_velocity_matrix_set") = false;

	// ************** PARAMETER TO TELL WHEN LINEAR ASSEMBLY HAS BEEN DONE ************ //
	es->parameters.set<bool> ("linear_assembly_done")   = false;

	// ************** EFFICIENT ASSEMBLY NOT SET UP FOR ADAPTIVE TIME-STEPPING *************** //
	// must make sure the time-step hasn't changed so no adaptive time-stepping
	// i'm not quite sure why this should be the case
	if(es->parameters.get<bool> ("adaptive_time_stepping") && es->parameters.get<bool> ("efficient_assembly"))
	{
		/*
		std::cout << "- Efficient assembly not set up for adaptive time-stepping yet." << std::endl;
		std::cout << "Changing to use regular assembly." << std::endl;
		es->parameters.set<bool> ("efficient_assembly") = false;
		*/
		std::cout << "- Efficient assembly may not set up for adaptive time-stepping yet." << std::endl;
		std::cout << "Use with caution." << std::endl;
		//es->parameters.set<bool> ("efficient_assembly") = false;
	}



	// ************** IF DOING A 1D SIMULATION THEN CAN'T GET SEGMENT LENGTH FROM MESH *************** //
	// must make sure the time-step hasn't changed so no adaptive time-stepping
	if(sim_type == 1 && es->parameters.get<bool> ("initial_segment_length_from_mesh"))
	{
		std::cout << "Can't get initial segment length from mesh if doing a 1d simulations." << std::endl;
		std::cout << "Setting initial_segment_length_from_mesh=false." << std::endl;
		es->parameters.set<bool> ("initial_segment_length_from_mesh") = false;
	}

	// ************** IF DOING A 1D SIMULATION THEN CAN'T MACTH TO MESH *************** //
	// must make sure the time-step hasn't changed so no adaptive time-stepping
	if(sim_type == 1 && es->parameters.get<bool> ("match_1d_mesh_to_3d_mesh"))
	{
		std::cout << "- Can't match 1d to 3d mesh if doing a 1d simulations." << std::endl;
		std::cout << "Setting match_1d_mesh_to_3d_mesh=false." << std::endl;
		es->parameters.set<bool> ("match_1d_mesh_to_3d_mesh") = false;
	}


	// *************** NO POINT DOING TWOD TO THREED IF DOING 3D ****************	//
	if(threed && (es->parameters.get<bool> ("twod_oned_tree") || es->parameters.get<bool> ("twod_oned_tree_pressure")))
	{
		std::cout << "- no need to do twodoned if running a 3D simulation." << std::endl;
		std::cout << "changing to false" << std::endl;
		es->parameters.set<bool> ("twod_oned_tree") = false;
		es->parameters.set<bool> ("twod_oned_tree_pressure") = false;
		std::cout << "twod_oned_tree = " << es->parameters.get<bool> ("twod_oned_tree") << std::endl;
		std::cout << "twod_oned_tree_pressure = " << es->parameters.get<bool> ("twod_oned_tree_pressure") << std::endl;
	}	

	// *************** NO POINT DOING TWOD TO THREED IF not coupled ****************	//
	// there is a point, basically if we want to simulate a 2D pipe as a 3D pipe
	if(false)//(sim_type == 0 || sim_type == 1) && (es->parameters.get<bool> ("twod_oned_tree") || es->parameters.get<bool> ("twod_oned_tree_pressure")))
	{
		std::cout << "- no need to do twodoned if not doing a coupled simulation." << std::endl;
		std::cout << "changing to false" << std::endl;
		es->parameters.set<bool> ("twod_oned_tree") = false;
		es->parameters.set<bool> ("twod_oned_tree_pressure") = false;
		std::cout << "twod_oned_tree = " << es->parameters.get<bool> ("twod_oned_tree") << std::endl;
		std::cout << "twod_oned_tree_pressure = " << es->parameters.get<bool> ("twod_oned_tree_pressure") << std::endl;
	}	

	if(es->parameters.get<bool> ("residual_formulation") && !es->parameters.get<bool> ("residual_formulation_0d") && sim_type == 5)
	{
		std::cout << "- must use 0d residual formulation when using monolithic with residual formulation." << std::endl;
		std::cout << "setting residual_formulation_0d to true." << std::endl;
		es->parameters.set<bool> ("residual_formulation_0d") = true;
	}


	// varying outlet pressure only works with sim_type 0 or sim_type 2 (no moghadam)
	if(es->parameters.get<bool> ("vary_outlet_pressure") && 
		!(sim_type == 0 || (sim_type == 2 && es->parameters.get<unsigned int> ("moghadam_coupling"))))
	{
		std::cout << "- varying outlet pressure only works with 3D and uncoupled 3D-0D calculations." << std::endl;
		std::cout << "you are currently using sim_type = " << sim_type << std::endl;
		std::cout << "and moghadam_coupling = " << es->parameters.get<unsigned int> ("moghadam_coupling") << std::endl;
		std::cout << "changing to just vary inlet pressure" << std::endl;
		es->parameters.set<bool> ("vary_outlet_pressure") = false;
	}

	// moghadam_coupling 2 not implemented yet
	if(false)//es->parameters.get<unsigned int> ("moghadam_coupling") == 2)
	{
		std::cout << "- moghadam_coupling 2 not implemented yet." << std::endl;
		std::cout << "changing to use moghadam_coupling 1." << std::endl;
		es->parameters.set<unsigned int> ("moghadam_coupling") = 1;
	}



	// only need to use moghadam_coupling 2 when we have nonlinear resistance
	if(es->parameters.get<unsigned int> ("moghadam_coupling") == 2 && 
		 (es->parameters.get<unsigned int>("resistance_type_1d") == 0 ||
			es->parameters.get<unsigned int>("resistance_type_1d") == 4 ||
			es->parameters.get<unsigned int>("resistance_type_1d") == 5))
	{
		/*
		std::cout << "- no need to use moghadam_coupling 2 when resistance_type_1d is 0 i.e. linear" << std::endl;
		std::cout << "changing to use moghadam_coupling 1" << std::endl;
		es->parameters.set<unsigned int> ("moghadam_coupling") = 1;
		*/

		std::cout << "- warning:no need to use moghadam_coupling 2 when resistance_type_1d is 0,4,5 i.e. linear" << std::endl;
		std::cout << "but okay" << std::endl;
		
	}




	// moghadam settings
	if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13)
	{
		// we want to be able to use this preconditioner formulation without moghadam too
		/*
		if(!es->parameters.get<bool> ("construct_moghadam_preconditioner"))
		{
			std::cout << "- if using moghadam preconditioner, must also construct the preconditioner." << std::endl;
			std::cout << "setting construct_moghadam_preconditioner to true." << std::endl;
			es->parameters.set<bool> ("construct_moghadam_preconditioner") = true;
		}
		*/

		if(!es->parameters.get<bool> ("residual_formulation"))
		{
			std::cout << "- must use residual formulation when using moghadam preconditioner (preconditioner 13)." << std::endl;
			std::cout << "setting residual_formulation to true." << std::endl;
			es->parameters.set<bool> ("residual_formulation") = true;
		}

		if(es->parameters.get<unsigned int> ("libmesh_jacobi") == 0)
		{
			std::cout << "- must use a jacobi preconditioner when using moghadam preconditioner (preconditioner 13)." << std::endl;
			std::cout << "setting to use jacobi symmetric preconditioning (3)." << std::endl;
			es->parameters.set<unsigned int> ("libmesh_jacobi") = 3;

		}
	}




	// actually we may want to, why not
	/*
	if(es->parameters.get<unsigned int> ("moghadam_preconditioner_type") != 0 && es->parameters.get<bool> ("construct_moghadam_preconditioner"))
	{
		std::cout << "- no need to construct moghadam preconditioner if using moghadam without pc (1) or moghadam exact (2)" << std::endl;
		std::cout << "changing to not construct moghadam preconditioner" << std::endl;
		es->parameters.set<bool> ("construct_moghadam_preconditioner") = false;

	}
	*/

	if(false)//es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13 && !es->parameters.get<unsigned int> ("moghadam_coupling"))
	{
		std::cout << "- if using moghadam preconditioner (13) then we must use moghadam coupling." << std::endl;
		std::cout << "changing to use moghadam coupling (2)" << std::endl;
		es->parameters.set<unsigned int> ("moghadam_coupling") = 2;
	}

	if(es->parameters.get<unsigned int> ("moghadam_coupling") == 0 && es->parameters.get<bool> ("construct_moghadam_preconditioner"))
	{
		std::cout << "- no need to construct moghadam preconditioner or use moghadam preconditioner is not using moghadam coupling" << std::endl;
		std::cout << "changing to not construct moghadam preconditioner" << std::endl;
		es->parameters.set<bool> ("construct_moghadam_preconditioner") = false;
	}

	if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13 && es->parameters.get<unsigned int> ("moghadam_preconditioner_type") == 0 && !es->parameters.get<bool> ("construct_moghadam_preconditioner"))
	{
		std::cout << "- need to construct moghadam preconditioner if using moghadam with pc (0)" << std::endl;
		std::cout << "changing to construct moghadam preconditioner" << std::endl;
		es->parameters.set<bool> ("construct_moghadam_preconditioner") = true;

	}

	
	if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13 && es->parameters.get<unsigned int> ("moghadam_preconditioner_type") != 0 && es->parameters.get<bool> ("construct_moghadam_preconditioner"))
	{
		std::cout << "- no need to construct moghadam preconditioner if using moghadam without pc or exact pc (1)" << std::endl;
		std::cout << "changing to not construct moghadam preconditioner" << std::endl;
		es->parameters.set<bool> ("construct_moghadam_preconditioner") = false;

	}

	


	if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13 && es->parameters.get<unsigned int> ("moghadam_preconditioner_type") == 0 && es->parameters.get<unsigned int> ("libmesh_jacobi") != 3)
	{
		std::cout << "- can only do moghdam preconditioner with pc when using symmetric jacobi." << std::endl;
		std::cout << "changing to use symmetric jacobi" << std::endl;
		es->parameters.set<unsigned int> ("libmesh_jacobi") = 3;
	}

	if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13 && !es->parameters.get<bool> ("residual_linear_solve"))
	{
		std::cout << "- can only do moghdam preconditioner (13) with residual linear solve." << std::endl;
		std::cout << "changing to use residual linear solve" << std::endl;
		es->parameters.set<bool> ("residual_linear_solve") = true;		
	}

	if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13 && es->parameters.get<bool> ("moghadam_velocity_preconditioner") && !es->parameters.get<bool> ("construct_moghadam_preconditioner"))
	{
		std::cout << "- if using moghadam velocity preconditioner then we need to construct the moghadam preconditioner." << std::endl;
		std::cout << "changing to construct moghdadam preconditioner" << std::endl;
		es->parameters.set<bool> ("construct_moghadam_preconditioner") = true;
	}



	// actually, we may want it to be possible to use preconditioner_type_3d == 13 (moghadam preconditioner) even if no moghadam coupling
	/*
	if(es->parameters.get<unsigned int> ("moghadam_coupling") == 0 && es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13)
	{
		std::cout << "- no need to construct moghadam preconditioner or use moghadam preconditioner is not using moghadam coupling" << std::endl;
		std::cout << "changing to not construct moghadam preconditioner" << std::endl;
		es->parameters.set<bool> ("construct_moghadam_preconditioner") = false;
		es->parameters.set<unsigned int> ("moghadam_preconditioner") = false;
	}
	*/

	// if linear resistance type then no point doing partitioned converge where we onverge on the correct solution
	if((es->parameters.get<unsigned int>("resistance_type_1d") == 0 ||
			es->parameters.get<unsigned int>("resistance_type_1d") == 4 ||
			es->parameters.get<unsigned int>("resistance_type_1d") == 5) && 
			es->parameters.get<bool> ("partitioned_converge"))
	{
		std::cout << "- no need to converge 0d solution if using a linear resistance law (0)" << std::endl;
		std::cout << "changing partitioned_converge to false." << std::endl;
		es->parameters.set<bool> ("partitioned_converge") = false;
	}

	// if linear resistance type then no point doing partitioned converge where we onverge on the correct solution
	if((es->parameters.get<unsigned int> ("resistance_type_1d") == 4 ||
			es->parameters.get<unsigned int> ("resistance_type_1d") == 5)
			&&  es->parameters.get<bool> ("_read_1d_mesh"))
	{
		std::cout << "- hard-coded resistances not implemented when reading in 0d mesh (doesn't make sense)." << std::endl;
		std::cout << "changing to use pedley resistance." << std::endl;
		es->parameters.set<unsigned int> ("resistance_type_1d") = 1;
	}


	// if linear resistance type then no point doing partitioned converge where we onverge on the correct solution
	if(es->parameters.get<unsigned int> ("resistance_type_1d") == 4 && es->parameters.get<unsigned int> ("num_generations") > 4)
	{
		std::cout << "- hard-coded resistances not implemented for more than 4 generations." << std::endl;
		std::cout << "changing to use pedley resistance." << std::endl;
		es->parameters.set<unsigned int> ("resistance_type_1d") = 1;
	}


	// if linear resistance type then no point doing partitioned converge where we onverge on the correct solution
	if(es->parameters.get<unsigned int> ("resistance_type_1d") == 5 && es->parameters.get<unsigned int> ("num_generations") > 1)
	{
		std::cout << "- asymmetric hard-coded resistances not implemented for more than 1 generation." << std::endl;
		std::cout << "changing to use pedley resistance." << std::endl;
		es->parameters.set<unsigned int> ("resistance_type_1d") = 1;
	}





	if(es->parameters.get<bool> ("vary_outlet_pressure") && es->parameters.get<bool> ("apply_pressure_from_flux"))
	{
		std::cout << "- varying outlet pressure and applying pressure from flux not supported." << std::endl;
		std::cout << "changing to just vary inlet pressure" << std::endl;
		es->parameters.set<bool> ("vary_outlet_pressure") = false;
	}

	if(es->parameters.get<bool> ("residual_linear_solve") && !es->parameters.get<bool> ("residual_formulation"))
	{
			std::cout << "- residual linear solve only supported with residual formulation." << std::endl;
			std::cout << "Changing to residual formulation to true" << std::endl;
			es->parameters.set<bool> ("residual_formulation") = true;
		
	}

	// ************* residual formulation restrictions
	if(es->parameters.get<bool> ("residual_formulation"))
	{
		// only implemented for efficient assembly
		if(!es->parameters.get<bool> ("efficient_assembly"))
		{
			std::cout << "- residual formulation only supported with efficient assembly." << std::endl;
			std::cout << "Exiting..." << std::endl;
			std::exit(0);
		}

		// not implemented for monolithic
		// it is now
		if(false)//sim_type == 5)
		{
			std::cout << "- residual formulation not supported with monolithic." << std::endl;
			std::cout << "changing to not use residual formulation." << std::endl;
			es->parameters.set<bool> ("residual_formulation") = false;
		}

		// only implemented dor efficient assembly
		if(es->parameters.get<bool> ("symmetric_gradient"))
		{
			std::cout << "- residual formulation not supported with symmetric gradient." << std::endl;
			std::cout << "Exiting..." << std::endl;
			std::exit(0);
		}

		if(es->parameters.get<unsigned int> ("moghadam_coupling") == 1)
		{
			std::cout << "- residual formulation not supported with moghadam_coupling 1 (strictly linear)." << std::endl;
			std::cout << "Exiting..." << std::endl;
			std::exit(0);
		}

		if(es->parameters.get<bool> ("nonzero_initial_guess"))
		{
			std::cout << "- non need to use nonzero initial guess when using residual formulation." << std::endl;
			std::cout << "setting to not use nonzero initial  guess." << std::endl;
			es->parameters.set<bool> ("nonzero_initial_guess") = false;
		}

		if(es->parameters.get<bool> ("nonzero_initial_guess_inner_3d1d"))
		{
			std::cout << "- non need to use nonzero initial guess inner 3d1d when using residual formulation." << std::endl;
			std::cout << "setting to not use nonzero initial  guess." << std::endl;
			es->parameters.set<bool> ("nonzero_initial_guess_inner_3d1d") = false;
		}
	}

	// **********
	if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 0 && es->parameters.get<unsigned int> ("libmesh_jacobi") != 0)
	{
		std::cout << "- No need to use libmesh implemented jacobi when using a direct solver." << std::endl;
		std::cout << "Setting to not use libmesh implemented jacobi." << std::endl;
		es->parameters.set<unsigned int> ("libmesh_jacobi") = 0;
	}

	if(sim_type == 5 && es->parameters.get<unsigned int> ("libmesh_jacobi") != 0)
	{
		std::cout << "- Libmesh jacobi not implemented for monolithic yet." << std::endl;
		std::cout << "Changing to not use Jacobi preconditioning." << std::endl;
		es->parameters.set<unsigned int> ("libmesh_jacobi") = 0;

	}

	if(es->parameters.get<bool> ("compute_eigenvalues") && (es->parameters.get<unsigned int>("preconditioner_type_3d") == 13 || es->parameters.get<bool>("residual_linear_solve")))
	{
		std::cout << "- Compute eigenvalues not configured for residual solve yet" << std::endl;
		std::cout << "changing to not calculate eigenvalues." << std::endl;
		es->parameters.set<bool> ("compute_eigenvalues") = 0;
		
	}

	if(es->parameters.get<unsigned int> ("acinar_model") == 1)
	{
		if(es->parameters.get<unsigned int> ("additional_symmetric_generations"))
		{
			std::cout << "- additional symmetric generations not configured for use with acinar boundary conditions." << std::endl;
			std::cout << "changing to not use additional symmetric generations." << std::endl;
			es->parameters.set<unsigned int> ("additional_symmetric_generations") = 0;
		}

		if(restart)
		{
			std::cout << "- restart not configured for use with acinar boundary conditions." << std::endl;
			std::cout << "Exiting..." << std::endl;
			std::exit(0);

		}
	}


	// a derived parameter
	if(es->parameters.get<unsigned int> ("resistance_type_1d") == 1 ||
			es->parameters.get<unsigned int> ("resistance_type_1d") == 2 ||
			es->parameters.get<unsigned int> ("resistance_type_1d") == 3)
	{
		es->parameters.set<bool> ("nonlinear_resistance_type_1d") = true;
	}

	// only use mono flow rate penalty with monolithic
	if(es->parameters.get<bool> ("mono_flow_rate_penalty") && sim_type != 5)
	{
		std::cout << "- Shouldn't use mono_flow_rate_penalty if not doing monolithic sim." << std::endl;
		std::cout << "Changing to not use mono_flow_rate_penalty." << std::endl;
		es->parameters.set<bool> ("mono_flow_rate_penalty") = false;
	}

	if(es->parameters.get<bool> ("resistance_adjusted_radius") && es->parameters.get<bool> ("threed"))
	{
		std::cout << "- no need to have resistance adjusted radius for 3D simulations." << std::endl;
		std::cout << "Changing to not use resistance_adjusted_radius." << std::endl;
		es->parameters.set<bool> ("resistance_adjusted_radius") = false;		
	}

	if(es->parameters.get<bool> ("first_gen_poiseuille") && sim_3d)
	{
		std::cout << "- no need to use poiseuille in first generation if doing 2d/3d sim." << std::endl;
		std::cout << "Changing to not use first_gen_poiseuille." << std::endl;
		es->parameters.set<bool> ("first_gen_poiseuille") = false;		
	}



	// ************** SETUP MESH REFINEMENT SETTINGS ******************* //

	if(es->parameters.get<unsigned int> ("output_no_refinement") > es->parameters.get<unsigned int> ("no_refinement") )
	{
		std::cout << "- Cannot output solution at a more refined state (need to fix dirichlet condition stuff). Setting output refinement level to same as simulation." << std::endl;
		std::cout << "\toutput_no_refinement = " << es->parameters.get<unsigned int> ("no_refinement") << std::endl;
		es->parameters.set<unsigned int> ("output_no_refinement") = es->parameters.get<unsigned int> ("no_refinement");
	}



  mesh_refinement.refine_fraction() = es->parameters.get<double>("refine_fraction");
  mesh_refinement.coarsen_fraction() = es->parameters.get<double>("coarsen_fraction");
  mesh_refinement.coarsen_threshold() = es->parameters.get<double>("coarsen_threshold");
  mesh_refinement.max_h_level() = es->parameters.get<unsigned int>("max_h_level");
	mesh_refinement.nelem_target() = es->parameters.get<unsigned int>("nelem_target");
	mesh_refinement.coarsen_by_parents() = true;
	mesh_refinement.face_level_mismatch_limit() = es->parameters.get<unsigned int>("face_level_mismatch_limit");		//without this can really get stuck...

  mesh_refinement_3d.refine_fraction() = es->parameters.get<double>("refine_fraction");
  mesh_refinement_3d.coarsen_fraction() = es->parameters.get<double>("coarsen_fraction");
  mesh_refinement_3d.coarsen_threshold() = es->parameters.get<double>("coarsen_threshold");
  mesh_refinement_3d.max_h_level() = es->parameters.get<unsigned int>("max_h_level");
	mesh_refinement_3d.nelem_target() = es->parameters.get<unsigned int>("nelem_target");
	mesh_refinement_3d.coarsen_by_parents() = true;
	mesh_refinement_3d.face_level_mismatch_limit() = es->parameters.get<unsigned int>("face_level_mismatch_limit");		//without this can really get stuck...
	//but perhaps with a larger mesh this is not the case
	es->parameters.set<bool>("mesh_refined") = false;

	/*
	if(es->parameters.get<unsigned int>("num_refinement_steps") > 0 && 
			!es->parameters.get<unsigned int>("unsteady"))
		es->parameters.set<double>("end_time") = es->parameters.get<unsigned int>("num_refinement_steps") + 2;
	*/



	// ************ SETUP STAGE NUMBERS FOR PETSC PROFILING *************** //
	PetscErrorCode ierr;	// unused
	PetscLogStage stage_num_3d_setup;
	if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 0)
		ierr = PetscLogStageRegister("Exact Solver Dummy",&stage_num_3d_setup);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 2)
		ierr = PetscLogStageRegister("Petsc LSC (BFBt) Preconditioner",&stage_num_3d_setup);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 3)
		ierr = PetscLogStageRegister("Pressure Preconditioner",&stage_num_3d_setup);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 4)
		ierr = PetscLogStageRegister("PCD Preconditioner",&stage_num_3d_setup);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 5)
		ierr = PetscLogStageRegister("PCD2 Preconditioner",&stage_num_3d_setup);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 6)
		ierr = PetscLogStageRegister("LSC James Preconditioner",&stage_num_3d_setup);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 7)
		ierr = PetscLogStageRegister("LSC James Scaled Preconditioner",&stage_num_3d_setup);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 9)
		ierr = PetscLogStageRegister("LSC James Scaled Stabilised Preconditioner Setup",&stage_num_3d_setup);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 10)
		ierr = PetscLogStageRegister("SIMPLE Preconditioner",&stage_num_3d_setup);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 11)
		ierr = PetscLogStageRegister("SIMPLER Preconditioner",&stage_num_3d_setup);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 12)
		ierr = PetscLogStageRegister("SIMPLERC Preconditioner",&stage_num_3d_setup);


	PetscLogStage stage_num_3d_apply;
	if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 9)
		ierr = PetscLogStageRegister("LSC James Scaled Stabilised Preconditioner Apply",&stage_num_3d_apply);

	

	PetscLogStage stage_num_3d1d;
	if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 6)
		ierr = PetscLogStageRegister("Mono Diagonal Preconditioner",&stage_num_3d1d);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 7)
		ierr = PetscLogStageRegister("Mono Full Preconditioner",&stage_num_3d1d);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 8)
		ierr = PetscLogStageRegister("Mono Schur Navier Stokes Velocity Preconditioner",&stage_num_3d1d);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 9)
		ierr = PetscLogStageRegister("Mono Schur Preconditioner",&stage_num_3d1d);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 10)
		ierr = PetscLogStageRegister("Mono Schur Stokes Preconditioner",&stage_num_3d1d);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 11)
		ierr = PetscLogStageRegister("Mono Schur Stokes Velocity Preconditioner",&stage_num_3d1d);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 12)
		ierr = PetscLogStageRegister("Mono Schur 0D Preconditioner",&stage_num_3d1d);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 13)
		ierr = PetscLogStageRegister("Mono Schur NS Approx Preconditioner",&stage_num_3d1d);



	particle_deposition = set_unsigned_int_parameter(infile,"particle_deposition",0);





	return 0;

}



// output parameters to screen and to file.
void NavierStokesCoupled::print_parameters()
{

	// ***************** PRINT OUT THE PARAMETERS AS THEY CURRENTLY ARE ******** //
	es->parameters.print(std::cout);

}

// output parameters to screen and to file.
void NavierStokesCoupled::output_parameters()
{

	// **************** COPY THE INPUT FILE TO THE OUTPUT FOLDER **************** //
	if(!es->parameters.get<bool>("compare_results"))
	{
		std::ifstream  input_file_src(input_file.c_str(), std::ios::binary);
		std::ofstream  input_file_dst(std::string(es->parameters.get<std::string>("output_folder") + "navier.in").c_str(),   std::ios::binary);
		input_file_dst << input_file_src.rdbuf();

		std::ifstream  mesh_file_src(es->parameters.get<std::string>("mesh_file").c_str(), std::ios::binary);
		std::ofstream  mesh_file_dst(std::string(es->parameters.get<std::string>("output_folder") + "mesh_file.msh").c_str(),   std::ios::binary);
		mesh_file_dst << mesh_file_src.rdbuf();
	}

}



// output command line options to file.
void NavierStokesCoupled::output_command_line_options()
{

	// **************** OUTPUT COMMAND LINE OPTIONS TO FILE IN OUTPUT FOLDER **************** //
	if(!es->parameters.get<bool>("compare_results"))
	{
		std::ofstream  command_line_options_file_dst(std::string(es->parameters.get<std::string>("output_folder") + "command_line_options.dat").c_str(),   std::ios::binary);
		comm_line.print(command_line_options_file_dst);
	}

}



// set autofieldsplit options
// NOTE: doesn't work...
void NavierStokesCoupled::set_auto_fieldsplit_parameters()
{
	// set to use solver_variable_names and solver_system_names 
	// (doesn't really matter if this is duplicated..)
	comm_line.set(std::string("--solver_variable_names").c_str(),"");
	comm_line.set(std::string("--solver_system_names").c_str(),"");

	// set variable numberings
	// monolithic stuff
	set_string_command_line_parameter("--solver_group_ns3d1d_u","0");
	set_string_command_line_parameter("--solver_group_ns3d1d_v","0");
	set_string_command_line_parameter("--solver_group_ns3d1d_w","0");
	set_string_command_line_parameter("--solver_group_ns3d1d_p","0");
	set_string_command_line_parameter("--solver_group_ns3d1d_Q","1");
	std::string var = set_string_command_line_parameter("--solver_group_ns3d1d_P","1");
	set_string_command_line_parameter("--solver_group_ns3d1d_0_u","0");
	set_string_command_line_parameter("--solver_group_ns3d1d_0_v","0");
	set_string_command_line_parameter("--solver_group_ns3d1d_0_w","0");
	set_string_command_line_parameter("--solver_group_ns3d1d_0_p","1");

	// non-monolithic stuff
	set_string_command_line_parameter("--solver_group_ns3d_u","0");
	set_string_command_line_parameter("--solver_group_ns3d_v","0");
	set_string_command_line_parameter("--solver_group_ns3d_w","0");
	set_string_command_line_parameter("--solver_group_ns3d_p","1");


  if (libMesh::on_command_line("--solver_system_names"))
	{
		std::cout << "solver system names working" << std::endl;
	}

  if (libMesh::on_command_line("--solver_variable_names"))
	{
		std::cout << "solver variable names working" << std::endl;
	}

	const std::string empty_string;

	std::string group_command = "--solver_group_ns3d1d_0_p";
  std::string group_name = libMesh::command_line_next
    (group_command, empty_string);

	std::cout << "group name of --solver_group_ns3d1d_0_p = " << group_name << std::endl;
}



// want functions to set values in parameters by reading from input file or overwritten by command line
// returns the value used
double NavierStokesCoupled::set_double_parameter(GetPot _infile, std::string name, double default_value)
{
	es->parameters.set<double>(name) = _infile(name,default_value);
	if(comm_line.search("-" + name))
		es->parameters.set<double>(name) = comm_line.next(default_value);
	return es->parameters.get<double>(name);
}

unsigned int NavierStokesCoupled::set_unsigned_int_parameter(GetPot _infile, std::string name, unsigned int default_value)
{
	es->parameters.set<unsigned int>(name) = _infile(name,default_value);
	if(comm_line.search("-" + name))
		es->parameters.set<unsigned int>(name) = comm_line.next(default_value);
	return es->parameters.get<unsigned int>(name);
}

int NavierStokesCoupled::set_int_parameter(GetPot _infile, std::string name, int default_value)
{
	es->parameters.set<int>(name) = _infile(name,default_value);
	if(comm_line.search("-" + name))
		es->parameters.set<int>(name) = comm_line.next(default_value);
	return es->parameters.get<int>(name);
}

bool NavierStokesCoupled::set_bool_parameter(GetPot _infile, std::string name, bool default_value)
{
	es->parameters.set<bool>(name) = _infile(name,default_value);
	if(comm_line.search("-" + name))
		es->parameters.set<bool>(name) = comm_line.next(default_value);
	return es->parameters.get<bool>(name);
}

std::string NavierStokesCoupled::set_string_parameter(GetPot _infile, std::string name, std::string default_value)
{
	es->parameters.set<std::string>(name) = _infile(name,default_value);
	if(comm_line.search("-" + name))
		es->parameters.set<std::string>(name) = comm_line.next(default_value);
	return es->parameters.get<std::string>(name);
}

// want functions to set values in command line or overwritten by command line
// returns the value used
std::string NavierStokesCoupled::set_string_command_line_parameter(std::string name, std::string default_value)
{
	// if it isn't on the command line set to default
	if(!comm_line.search(name))
	{
		comm_line.set(name.c_str(),default_value);
		std::cout << "didn't find " << name << std::endl;
	}
	else
	{
		std::cout << "found " << name << std::endl;

	}
	
	// return the value
	return comm_line(name.c_str(),"no value found.");
}


// read in timestep sizes from the out.dat file
// we are assuming that
void NavierStokesCoupled::setup_read_timesteps()
{
	if(es->parameters.get<bool>("unsteady_from_steady"))
	{
		t_step = 0;
		es->parameters.set<unsigned int> ("t_step")   = t_step;

	}
	else if(es->parameters.get<bool>("read_time_steps"))
	{

		std::ifstream timestep_file((es->parameters.get<std::string>("restart_folder") + "out.dat").c_str());
		std::cout << "input file = " << (es->parameters.get<std::string>("restart_folder") + "out.dat").c_str() << std::endl;
		std::string line;
		double begin_timestep_time = 0.;
		double end_timestep_time = 0;
		// skip the first line
		if(std::getline(timestep_file, line))
		{
			while (std::getline(timestep_file, line))
			{
				std::istringstream iss(line);
				double timestep;
				if (!(iss >> timestep >> end_timestep_time)) { break; } // error

				timestep_sizes.push_back(end_timestep_time - begin_timestep_time);
				begin_timestep_time = end_timestep_time;		
			}
			es->parameters.set<double> ("end_time")   = end_timestep_time;
		}
		else
		{
			std::cout << "ERROR: no out.dat file to read bra... EXITING." << std::endl;
			std::exit(0);
		}

		// set initial time step
		t_step = timestep_sizes[0];
		es->parameters.set<unsigned int> ("t_step")   = t_step;

		/*
		std::cout << "CHECK TIMESTEP SIZES VECTOR" << std::endl;
		for(unsigned int i=0; i< timestep_sizes.size(); i++)
			std::cout << "i = " << i << " : " << timestep_sizes[i] << std::endl;
		*/
	}
	else	// just make our own timestep vector based on a given timestep
	{
		double read_dt = es->parameters.get<double> ("read_time_step_size");
		unsigned int num_read_time_steps = (unsigned int)(es->parameters.get<double> ("read_end_time")/read_dt) + 1;
		timestep_sizes.resize(num_read_time_steps,read_dt);

		std::cout << "CHECK TIMESTEP SIZES VECTOR" << std::endl;
		for(unsigned int i=0; i< timestep_sizes.size(); i++)
			std::cout << "i = " << i << " : " << timestep_sizes[i] << std::endl;
	}


	// make sure the particle dt is smaller than the fluid time step
	if(!es->parameters.get<bool>("unsteady_from_steady"))
	{
		for(unsigned int i=0; i<timestep_sizes.size(); i++)
		{
			if(timestep_sizes[i] < dt - 1e-10)
			{
				std::cout << "particle time step size is larger than fluid time step." << std::endl;
				std::cout << "particle_dt = " << dt << std::endl;
				std::cout << "fluid_dt = " << timestep_sizes[i] << std::endl;
				std::cout << "EXITING..." << std::endl;
				std::exit(0);
			}
				
		}
	}



}


// read the input boundary conditions for an uncoupled simulation
int NavierStokesCoupled::read_input_boundary_conditions_3d0d()
{
	std::cout << "\nReading input 3D-0D boundary conditions file." << std::endl;

	// read the file 
	std::string boundary_conditions_filename = es->parameters.get<std::string>("boundary_conditions_folder") + "out.dat";

	std::cout << " file = " << boundary_conditions_filename << std::endl;

	std::ifstream input_boundary_conditions_file(boundary_conditions_filename.c_str(), std::ios::binary);
	if(!input_boundary_conditions_file.good())
	{
		std::cout << "Error reading input boundary conditions file." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}

	// put data from file into vector of vectors all_input_pressure_values_3d and all_input_flux_values_0d
	// note that we actually use the 0d pressure values as boundary conditions on the 3d (boundary condition is not applied exactly).




	// 1 - read first line and figure out if there are the correct number of boundary values in the input file

	std::string first_line;
	std::getline(input_boundary_conditions_file,first_line);
	std::istringstream iss(first_line);
	std::string temp_string;
	unsigned int count = 0;
	while (iss.good()) // while the stream is not empty
  {
  	iss >> temp_string; //get data from the stream up to next whitespace
		count++;
  }
	std::cout << "count = " << count << std::endl;

	// we need to ignore 4 strings, #, timestep, current_time and num_newton_steps
	// then, if just 3D divide by 2 for the 3d fluxes and pressures.
	// if 3D-0D divide by 4 for the 3d and 0d fluxes and pressures.
	unsigned int num_input_boundary_ids = (count - 4)/4;	

	std::cout << "number of input boundary ids = " << num_input_boundary_ids << std::endl;
	std::cout << "number of boundary ids in current simulation = " << boundary_ids.size() << std::endl;
	if(num_input_boundary_ids != boundary_ids.size())
	{
		std::cout << "Error. Number of input boundary ids not equal to number of boundary ids in current simulation." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}







	// 2 - read all the lines and put the data where it's supposed to go

	// * discard the first 2 + 2*num_input_boundary_ids (time stuff and 3d fluxes)
	// * put the next num_input_boundary_ids into all_input_flux_values_0d
	// * put the next num_input_boundary_ids into all_input_pressure_values_3d
	// * next line

	std::string current_line;
	double disc = 0.;	// temp variable for discarding unused values in input file
	while(std::getline(input_boundary_conditions_file,current_line))
	{		
		std::istringstream ciss(current_line);

		// discard timestep and current time
		ciss >> disc >> disc;

		// discard the 3d fluxes and pressures for use as boundary conditions
		for(unsigned int i=0; i<2*num_input_boundary_ids; i++)
			ciss >> disc;


		// put the next num_input_boundary_ids into a vector for flux 0d
		std::vector<double> temp_input_flux_values_0d(num_input_boundary_ids);
		for(unsigned int i=0; i<num_input_boundary_ids; i++)
			ciss >> temp_input_flux_values_0d[i];

		// put the next num_input_boundary_ids into a vector for pressure 3d
		std::vector<double> temp_input_pressure_values_3d(num_input_boundary_ids);
		for(unsigned int i=0; i<num_input_boundary_ids; i++)
			ciss >> temp_input_pressure_values_3d[i];

		// put the vectors into global vectors for the timestep
		all_input_pressure_values_3d.push_back(temp_input_pressure_values_3d);
		all_input_flux_values_0d.push_back(temp_input_flux_values_0d);
		
	}

	
	// 3 - check that the number of timesteps is the same

	unsigned int total_num_timesteps = (unsigned int)(es->parameters.get<double>("end_time")/es->parameters.get<double>("dt"));
	unsigned int input_num_timesteps = all_input_pressure_values_3d.size();

	std::cout << "total expected number of timesteps = " << total_num_timesteps << std::endl;
	std::cout << "input number of timesteps = " << input_num_timesteps << std::endl;

	if(total_num_timesteps > input_num_timesteps)
	{
		std::cout << "Error, number of input timesteps bigger than the number of expected timesteps." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}

	// 4 - set input boundary conditions for first timestep

	std::cout << "t_step = " << t_step << std::endl;
	input_pressure_values_3d = all_input_pressure_values_3d[t_step];
	input_flux_values_1d = all_input_flux_values_0d[t_step];


}





// read the input boundary conditions for an uncoupled simulation
int NavierStokesCoupled::read_input_boundary_conditions_3d()
{
	std::cout << "\nReading input 3D boundary conditions file." << std::endl;

	// read the file 
	std::string boundary_conditions_filename = es->parameters.get<std::string>("boundary_conditions_folder") + "out.dat";

	std::cout << " file = " << boundary_conditions_filename << std::endl;

	std::ifstream input_boundary_conditions_file(boundary_conditions_filename.c_str(), std::ios::binary);
	if(!input_boundary_conditions_file.good())
	{
		std::cout << "Error reading input boundary conditions file." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}

	// put data from file into vector of vectors all_input_pressure_values_3d and all_input_flux_values_0d
	// note that we actually use the 0d pressure values as boundary conditions on the 3d (boundary condition is not applied exactly).




	// 1 - read first line and figure out if there are the correct number of boundary values in the input file

	std::string first_line;
	std::getline(input_boundary_conditions_file,first_line);
	std::istringstream iss(first_line);
	std::string temp_string;
	unsigned int count = 0;
	while (iss.good()) // while the stream is not empty
  {
  	iss >> temp_string; //get data from the stream up to next whitespace
		count++;
  }
	std::cout << "count = " << count << std::endl;

	// we need to ignore 4 strings, #, timestep, current_time and num_newton_steps
	// then, if just 3D divide by 2 for the 3d fluxes and pressures.
	// if 3D-0D divide by 4 for the 3d and 0d fluxes and pressures.
	unsigned int num_input_boundary_ids = (count - 4)/2;	

	std::cout << "number of input boundary ids = " << num_input_boundary_ids << std::endl;
	std::cout << "number of boundary ids in current simulation = " << boundary_ids.size() << std::endl;
	if(num_input_boundary_ids != boundary_ids.size())
	{
		std::cout << "Error. Number of input boundary ids not equal to number of boundary ids in current simulation." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}







	// 2 - read all the lines and put the data where it's supposed to go

	// * discard the first 2 (time stuff)
	// * put the next num_input_boundary_ids into all_input_flux_values_0d
	// * put the next num_input_boundary_ids into all_input_pressure_values_3d
	// * next line

	std::string current_line;
	double disc = 0.;	// temp variable for discarding unused values in input file
	while(std::getline(input_boundary_conditions_file,current_line))
	{		
		std::istringstream ciss(current_line);

		// discard timestep and current time
		ciss >> disc >> disc;

		// put the next num_input_boundary_ids into a vector for flux 3d
		// which we won't necessarily use
		std::vector<double> temp_input_flux_values_3d(num_input_boundary_ids);
		for(unsigned int i=0; i<num_input_boundary_ids; i++)
			ciss >> temp_input_flux_values_3d[i];

		// put the next num_input_boundary_ids into a vector for pressure 3d
		std::vector<double> temp_input_pressure_values_3d(num_input_boundary_ids);
		for(unsigned int i=0; i<num_input_boundary_ids; i++)
			ciss >> temp_input_pressure_values_3d[i];

		// put the vectors into global vectors for the timestep
		all_input_pressure_values_3d.push_back(temp_input_pressure_values_3d);
		all_input_flux_values_3d.push_back(temp_input_flux_values_3d);
		
	}

	
	// 3 - check that the number of timesteps is the same

	unsigned int total_num_timesteps = (unsigned int)(es->parameters.get<double>("end_time")/es->parameters.get<double>("dt"));
	unsigned int input_num_timesteps = all_input_pressure_values_3d.size();

	std::cout << "total expected number of timesteps = " << total_num_timesteps << std::endl;
	std::cout << "input number of timesteps = " << input_num_timesteps << std::endl;

	if(total_num_timesteps > input_num_timesteps)
	{
		std::cout << "Error, number of input timesteps bigger than number of expected timesteps." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}

	// 4 - set input boundary conditions for first timestep

	std::cout << "t_step = " << t_step << std::endl;
	input_pressure_values_3d = all_input_pressure_values_3d[t_step];
	input_flux_values_3d = all_input_flux_values_3d[t_step];


}


// read the input boundary conditions for an uncoupled simulation
void NavierStokesCoupled::setup_pressure_boundary_conditions()
{
	std::cout << "\nSetting up pressure boundary conditions." << std::endl;

	input_pressure_values_3d.resize(boundary_ids.size(),0.);

	// set the parent branch
	input_pressure_values_3d[0] = es->parameters.get<double>("pressure_mag_3d") * es->parameters.get<double>("time_scaling");
	

}







void NavierStokesCoupled::print_flux_and_pressure()
{

	double flux_scaling = 1.0;
	double pressure_scaling = 1.0;
	/*
	if(!es->parameters.get<bool>("reynolds_number_calculation"))
	{
		flux_scaling = es->parameters.get<double>("velocity_scale") * 
										pow(es->parameters.get<double>("length_scale"),2.0);
		pressure_scaling = es->parameters.get<double>("density") *
													pow(es->parameters.get<double>("velocity_scale"),2.0);
	}
	*/
	if(es->parameters.get<bool>("nondimensionalised"))
	{
		flux_scaling = es->parameters.get<double>("velocity_scale") * 
										pow(es->parameters.get<double>("length_scale"),2.0);

		pressure_scaling = es->parameters.get<double>("density") *
													pow(es->parameters.get<double>("velocity_scale"),2.0);
	}



	std::cout << "pressure_scaling = " << pressure_scaling << std::endl;

	if(sim_3d)
	{
		double total_daughter_flux_3d = 0.0;
		for(unsigned int i=1; i < flux_values_3d.size(); i++)
			total_daughter_flux_3d += flux_values_3d[i];

		std::cout << "3D flux values:    " << std::endl;
		if(boundary_id_to_tree_id.size()==0)
		{
			for(unsigned int i=0; i < flux_values_3d.size(); i++)
				std::cout << " " << flux_values_3d[i] << " (" << 
									flux_values_3d[i] * flux_scaling << " m^3/s)" << std::endl;

			std::cout << "3D pressure values:    " << std::endl;
			for(unsigned int i=0; i < pressure_values_3d.size(); i++)
				std::cout << " " << pressure_values_3d[i] << " (" <<
									pressure_values_3d[i] * pressure_scaling << " Pa)" << std::endl;

			if(es->parameters.get<bool> ("output_total_pressure"))
			{
				std::cout << "3D total pressure values:    " << std::endl;
				for(unsigned int i=0; i < total_pressure_values_3d.size(); i++)
					std::cout << " " << total_pressure_values_3d[i] << " (" <<
										total_pressure_values_3d[i] * pressure_scaling << " Pa)" << std::endl;
			}
		}
		else
		{
			for(unsigned int i=0; i < flux_values_3d.size(); i++)
				std::cout << " " << flux_values_3d[tree_id_to_boundary_id[i]] << " (" << 
									flux_values_3d[tree_id_to_boundary_id[i]] * flux_scaling << " m^3/s)" << std::endl;

			std::cout << "3D pressure values:    " << std::endl;
			for(unsigned int i=0; i < pressure_values_3d.size(); i++)
				std::cout << " " << pressure_values_3d[tree_id_to_boundary_id[i]] << " (" <<
									pressure_values_3d[tree_id_to_boundary_id[i]] * pressure_scaling << " Pa)" << std::endl;
		}
		std::cout << "total daughter flux 3D = " << total_daughter_flux_3d << " (" << 
									total_daughter_flux_3d * flux_scaling << " m^3/s)" << std::endl;
	}

	if(sim_1d)
	{
		double total_daughter_flux_1d = 0.0;
		for(unsigned int i=1; i < flux_values_1d.size(); i++)
			total_daughter_flux_1d += flux_values_1d[i];

		std::cout << "1D flux values:    " << std::endl;
		for(unsigned int i=0; i < flux_values_1d.size(); i++)
			std::cout << " " << flux_values_1d[i] << " (" << 
								flux_values_1d[i] * flux_scaling << " m^3/s)" << std::endl;

		std::cout << "1D pressure values:    " << std::endl;
		for(unsigned int i=0; i < pressure_values_1d.size(); i++)
			std::cout << " " << pressure_values_1d[i] << " (" <<
								pressure_values_1d[i] * pressure_scaling << " Pa)" << std::endl;

		if(es->parameters.get<unsigned int> ("acinar_model") == 1)
		{
			double pleural_pressure = ns_assembler->calculate_pleural_pressure(es->parameters.get<double>("time"));
			std::cout << "Pleural pressure = " << pleural_pressure << " (" <<
										pleural_pressure * pressure_scaling << " Pa)" << std::endl;
		}
	}
}




void NavierStokesCoupled::update_times()
{

	if(sim_3d && sim_type != 5)
		system_3d->time = time;

	if(sim_1d && sim_type != 5)
		system_1d->time = time;

	if(sim_type == 5)
		system_coupled->time = time;

	es->parameters.set<double> ("dt") = dt;
	es->parameters.set<double> ("time") = time;
	es->parameters.set<double> ("previous_dt") = previous_dt;
	es->parameters.set<unsigned int> ("t_step") = t_step;

	update_time_scaling();

}

void NavierStokesCoupled::update_time_scaling()
{
	double time_flow = 0.0;
	double previous_time_flow = 0.0;

	//does the time flow for different values of the parameter
	// 0 - steady
	// 1 - sinusoidal
	// 2 - quadratic then steady


	const int unsteady = es->parameters.get<unsigned int>("unsteady");
	const Real period	 = es->parameters.get<Real>("period");
	//const Real end_time = es->parameters.get<Real> ("end_time");

	const double ramp_duration = es->parameters.get<double>("ramp_duration");
	const bool quadratic_ramp = es->parameters.get<bool>("quadratic_ramp");


	if(unsteady == 0)
		time_flow = 1.0;
	else if(unsteady == 1)
		time_flow = sin(2*pi*time/period);
	else if(unsteady == 2)
	{
		//quadratic (up to max of one) then constant
		if(time < ramp_duration)
		{
			if(quadratic_ramp)
				time_flow = -4.0/pow(ramp_duration*2,2.0)*time*(time-ramp_duration*2);
			else	//linear
				time_flow = time/ramp_duration;
		}
		else
		{
			time_flow = 1.0;
		}
	}
	else if(unsteady == 3)
		time_flow = 1 - cos(2*pi*time/period);
	else if(unsteady == 4)
		time_flow = -sin(2*pi*time/period);
	else if(unsteady == 5)
		time_flow = cos(2*pi*time/period);	
	else if(unsteady == 6)
	{
		//quadratic (up to max of one) then constant
		if(time < period)
		{
			time_flow = 1.0;
		}
		else if(time < period + es->parameters.get<double>("hofmann_breath_hold"))
		{
			time_flow = 0.0;
		}
		else if(time < 2*period + es->parameters.get<double>("hofmann_breath_hold"))
		{
			time_flow = -1.0;
		}
		else
		{
			time_flow = 0.0;
		}
	}
	else
	{
		//if reach here just do 1.0
		time_flow = 1.0;
	}



	es->parameters.set<double>("time_scaling") = time_flow;

	if(es->parameters.get<bool> ("increment_boundary_conditions"))
	{

	
		// previous

		double previous_time = time - dt;
		if(previous_time < 0)
			previous_time = 0;


		if(unsteady == 0)
			previous_time_flow = 1.0;
		else if(unsteady == 1)
			previous_time_flow = sin(2*pi*previous_time/period);
		else if(unsteady == 2)
		{
			//quadratic (up to max of one) then constant
			if(previous_time < ramp_duration)
			{
				if(quadratic_ramp)
					previous_time_flow = -4.0/pow(ramp_duration*2,2.0)*previous_time*(previous_time-ramp_duration*2);
				else	//linear
					previous_time_flow = previous_time/ramp_duration;
			}
			else
			{
				previous_time_flow = 1.0;
			}
		}
		else if(unsteady == 3)
			previous_time_flow = 1 - cos(2*pi*previous_time/period);
		else if(unsteady == 4)
			previous_time_flow = -sin(2*pi*previous_time/period);
		else if(unsteady == 5)
			previous_time_flow = cos(2*pi*previous_time/period);	
		else if(unsteady == 6)
		{
			//quadratic (up to max of one) then constant
			if(previous_time < period)
			{
				previous_time_flow = 1.0;
			}
			else if(previous_time < period + es->parameters.get<double>("hofmann_breath_hold"))
			{
				previous_time_flow = 0.0;
			}
			else if(previous_time < 2*period + es->parameters.get<double>("hofmann_breath_hold"))
			{
				previous_time_flow = -1.0;
			}
			else
			{
				previous_time_flow = 0.0;
			}
		}
		else
		{
			//if reach here just do 1.0
			previous_time_flow = 1.0;
		}
	}

	es->parameters.set<double>("previous_time_scaling") = previous_time_flow;

}










// need to think about where to put this and when to update vectors etc after refinement
// because we really want to refine and use the previous solution! reduces residuals (though not conv rate)
void NavierStokesCoupled::end_timestep_admin()
{


	if(particle_deposition)
	{
		if(es->parameters.get<bool>("unsteady_from_steady"))
		{
			//do nothing
		}
		else
		{
			// the reading of solutions is done at the beginning of a time step

			// move the last time step back
			if(sim_3d && sim_type != 5)
			{
				*old_global_solution = *system_3d->solution;
				system_3d->update();	
				*system_3d->old_local_solution = *system_3d->current_local_solution;
			}

			if(sim_1d && sim_type != 5)
			{
				*system_1d->old_local_solution = *system_1d->current_local_solution;
				system_1d->update();	
			}
		}


	}
	else
	{
		// when we reduce dt or refine the mesh always make sure that solution has been updated to old solution so that
		// it is correctly projected. local stuff can handle itself
		std::cout << "Performing end of time step admin." << std::endl;

		//if didn't converge or something then go back to previous time_step
		if((reduce_dt && unsteady) || refine_mesh)
		{
			std::cout << "Reducing time step size:" << std::endl;
			// check that isn't at minimum time step already
			if(fabs(dt) - 1e-10 < es->parameters.get<double>("min_time_step"))
			{
				std::cout << "At minimum time step already." << std::endl;
				std::cout << "Exiting..." << std::endl;
				exit_program = true;
			}

			time -= dt; // go back
			t_step -= 1;
		
			if(reduce_dt)
			{
				dt = es->parameters.get<double>("adaptive_time_step_reduce_factor")*dt;
				if(dt < es->parameters.get<double>("min_time_step"))
				{
					std::cout << "Cannot reduce time step passed minimum time step = " << es->parameters.get<double>("min_time_step") << std::endl;
					std::cout << "Setting to minimum time step." << std::endl;
					dt = es->parameters.get<double>("min_time_step");
				}
				else
					std::cout << " Reducing time step to dt = " << dt << std::endl;
			}

		
			update_times();

			if(sim_3d && sim_type != 5)
			{
				*system_3d->current_local_solution = *system_3d->old_local_solution;
				*system_3d->solution = *old_global_solution;
				//*old_global_solution = *system_3d->solution;	// remember we changed solution to old global before refinement
				system_3d->update();	//not sure why this reverts, because it somehow uses the solution and then distributes
				//note that it means here the solution is not correct
				// then after we're gonna need to make old_global_solution have the same size as solution
			}

			if(sim_1d && sim_type != 5)
			{
				*system_1d->current_local_solution = *system_1d->old_local_solution;
				*system_1d->solution = *old_global_solution_1d;
				system_1d->update();	// this 
			}

			if(sim_type == 5)
			{
				*system_coupled->current_local_solution = *system_coupled->old_local_solution;
				*system_coupled->solution = *old_global_solution;

			}

			steps_since_last_dt_change = 1;
			reduce_dt = 0;
		}

		//if did converge but a joke then double time step
		// (first_residual < 0.5 && system_3d->time > 30) removed this condition, made it just based on the number of iterations
		else if(increase_dt && unsteady)
		{
			std::cout << "Increasing time step size:" << std::endl;

			previous_flux_values_1d = flux_values_1d;
			previous_previous_flux_values_3d = previous_flux_values_3d;
			previous_flux_values_3d = flux_values_3d;
			previous_pressure_values_1d = pressure_values_1d;
			previous_previous_pressure_values_3d = previous_pressure_values_3d;
			previous_pressure_values_3d = pressure_values_3d;

			previous_dt = dt;
			update_times();
		

			dt = es->parameters.get<double>("adaptive_time_step_increase_factor")*dt;

			if(dt > es->parameters.get<double>("max_time_step"))
			{
					std::cout << "Cannot increase time step passed maximum time step = " << es->parameters.get<double>("min_time_step") << std::endl;
					std::cout << "Setting to maximum time step." << std::endl;
				dt = es->parameters.get<double>("max_time_step");
			}
			else
					std::cout << " Increasing time step to dt = " << dt << std::endl;

			update_times();


			// why do this??
			if(sim_3d && sim_type != 5)
			{
				*system_3d->old_local_solution = *system_3d->current_local_solution;
				*old_global_solution = *system_3d->solution;
				//system_3d->update();
			}

			if(sim_1d && sim_type != 5)
			{
				system_1d->update();	
				*old_global_solution_1d = *system_1d->solution;
				*system_1d->old_local_solution = *system_1d->current_local_solution;
				system_1d->update();	
			}

			if(sim_type == 5)
			{
				*old_global_solution = *system_coupled->solution;
				system_coupled->update();	
				*system_coupled->old_local_solution = *system_coupled->current_local_solution;
				system_coupled->update();	
			}

			steps_since_last_dt_change = 1;
			increase_dt = false;

			std::cout << " Increasing time step to " << dt << std::endl;
		}
		//if we've been worked hard enough then just 
		else
		{	

			previous_flux_values_1d = flux_values_1d;
			previous_previous_flux_values_3d = previous_flux_values_3d;
			previous_flux_values_3d = flux_values_3d;
			previous_pressure_values_1d = pressure_values_1d;
			previous_previous_pressure_values_3d = previous_pressure_values_3d;
			previous_pressure_values_3d = pressure_values_3d;

			previous_dt = dt;
			update_times();


			if(sim_3d && sim_type != 5)
			{
				*old_global_solution = *system_3d->solution;
				system_3d->update();	
				*system_3d->old_local_solution = *system_3d->current_local_solution;
			}

			if(sim_1d && sim_type != 5)
			{
				system_1d->update();	
				*old_global_solution_1d = *system_1d->solution;
				*system_1d->old_local_solution = *system_1d->current_local_solution;
				system_1d->update();	
			}

			if(sim_type == 5)
			{
				*old_global_solution = *system_coupled->solution;
				system_coupled->update();	
				*system_coupled->old_local_solution = *system_coupled->current_local_solution;
				system_coupled->update();	
			}

			steps_since_last_dt_change++;
			std::cout << " Keeping time step at " << dt << std::endl;
		}
	}
}






void NavierStokesCoupled::output_sim_data(bool header)
{

	double flux_scaling = 1.0;
	double pressure_scaling = 1.0;
	double volume_scaling = 1.0;

	// hmm, not sure whethere this should be false or true.., but works for now.
	if(es->parameters.get<bool>("nondimensionalised"))
	{
		flux_scaling = es->parameters.get<double>("velocity_scale") * 
										pow(es->parameters.get<double>("length_scale"),2.0);
		pressure_scaling = es->parameters.get<double>("density") *
													pow(es->parameters.get<double>("velocity_scale"),2.0);
		volume_scaling = pow(es->parameters.get<double>("length_scale"),3.0);
	}

	if(header)
	{
		// write geometry variables later when reading meta data file
		// write boundary conditions later with Picard class
		output_file << "# timestep\tcurrent_time";
		if(sim_3d)
		{
			for(unsigned int i=0; i < flux_values_3d.size(); i++)
				output_file <<	"\t3d_flux_value_" << i;
			for(unsigned int i=0; i < pressure_values_3d.size(); i++)
				output_file <<	"\t3d_press_value_" << i;
		}
		if(sim_1d)
		{
			for(unsigned int i=0; i < flux_values_1d.size(); i++)
				output_file <<	"\t1d_flux_value_" << i;

			for(unsigned int i=0; i < pressure_values_1d.size(); i++)
				output_file <<	"\t1d_press_value_" << i;
			
			if(es->parameters.get<unsigned int>("acinar_model") == 1)
			{
				output_file <<	"\tpleural_pressure";
				output_file <<	"\ttotal_acinar_volume";
				output_file <<	"\taverage_acinar_pressure";
			}
			
		}

		if(sim_3d)
			output_file << "\tnum_newton_steps";

		if(es->parameters.get<bool>("adaptive_time_stepping"))
			output_file << "\tadaptive_reason";

		output_file << std::endl;
	}
	else
	{

		std::cout << "k" << std::endl;

		output_file << t_step << "\t" << time;

		if(sim_3d)
		{
			std::cout << "k" << std::endl;
			calculate_3d_boundary_values();
			std::cout << "k" << std::endl;
			if(boundary_id_to_tree_id.size() == 0)
			{
				for(unsigned int i=0; i < flux_values_3d.size(); i++)
					output_file << "\t" << flux_values_3d[i] * flux_scaling;

				for(unsigned int i=0; i < pressure_values_3d.size(); i++)
					output_file << "\t" << pressure_values_3d[i] * pressure_scaling;
			}
			else
			{
				// output should be in order
				for(unsigned int i=0; i < flux_values_3d.size(); i++)
					output_file << "\t" << flux_values_3d[i] * flux_scaling;

				for(unsigned int i=0; i < pressure_values_3d.size(); i++)
					output_file << "\t" << pressure_values_3d[i] * pressure_scaling;
			}
		std::cout << "k" << std::endl;
		}
		
		if(sim_1d)
		{
			calculate_1d_boundary_values();
			// always in order of tree id
			for(unsigned int i=0; i < flux_values_1d.size(); i++)
				output_file << "\t" << flux_values_1d[i] * flux_scaling;

			for(unsigned int i=0; i < pressure_values_1d.size(); i++)
				output_file << "\t" << pressure_values_1d[i] * pressure_scaling;

			if(es->parameters.get<unsigned int>("acinar_model") == 1)
			{
				output_file <<	"\t" << ns_assembler->calculate_pleural_pressure(es->parameters.get<double>("time")) * pressure_scaling;
				output_file <<	"\t" << calculate_total_acinar_volume() * volume_scaling;
				output_file <<	"\t" << calculate_average_terminal_pressure() * pressure_scaling;
			}

		}

		if(sim_3d)
			output_file << "\t" <<	es->parameters.get<unsigned int> ("num_newton_steps");

		if(es->parameters.get<bool>("adaptive_time_stepping"))
		{
			// reduce time stepping reasons are negative
			// -1 - nonlinear, -2 linear, -3 residual
			// increase there is only one.
			int reason = 0;
			if(reduce_dt)
				reason = -1 * (int)reduce_dt;
			else if(increase_dt)
				reason = 1;

			output_file << "\t" <<	reason;
		}
		output_file << std::endl;

		std::cout << "k" << std::endl;
		print_flux_and_pressure();

			std::cout << "flux, pressure etc output for timestep " << t_step
							<< " written." << std::endl;

		std::cout << "hi" << std::endl;
		
		if(sim_1d && es->parameters.get<bool> ("output_resistance"))
		{
			output_poiseuille_resistance_per_generation();
		}
		

		std::cout << "bye" << std::endl;
	}
}


void NavierStokesCoupled::output_coupling_points()
{

	// output the parameters used to file and the header of the out.dat file
	std::ostringstream coupling_points_file_name;

	coupling_points_file_name << output_folder.str() << "coupling_points.csv";
	std::ofstream coupling_points_file;
	coupling_points_file.open(coupling_points_file_name.str().c_str());

	coupling_points_file << "x_coord y_coord z_coord coupling_point_id" << std::endl;
	for(unsigned int i=1; i < coupling_points.size(); i++)
		coupling_points_file << coupling_points[i](0) << " " << coupling_points[i](1) << " " << coupling_points[i](2) <<  " "  << i << std::endl;

	coupling_points_file.close();

}



void NavierStokesCoupled::output_linear_iteration_count(bool header)
{

	if(header)
	{
		linear_iterations_output_file << "# Output file with linear iteration counts" << std::endl;
		linear_iterations_output_file << "# Simple format with no parameters etc for ease of parsing" << std::endl;
		linear_iterations_output_file << "# Find simulation details in out.dat" << std::endl;
		// write geometry variables later when reading meta data file
		// write boundary conditions later with Picard class
		linear_iterations_output_file << "# Results" << std::endl;
		linear_iterations_output_file << "# timestep\tnonlinear_iteration\tlocal_linear_iterations\tmax_linear_iterations\ttotal_nonlinear_iterations\ttotal_linear_iterations\tnonlinear_residual";
		
		// extra ones for moghadam
		if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13)
			linear_iterations_output_file << "\ttotal_gmres_linear_iterations\ttotal_cg_linear_iterations";
		
		linear_iterations_output_file << std::endl;
	}
	else
	{
		linear_iterations_output_file << t_step << "\t" << nonlinear_iteration << 
																							"\t" << local_linear_iterations << 
																							"\t" << total_max_iterations << 
																							"\t" << total_nonlinear_iterations <<
																							"\t" << total_linear_iterations << 
																							"\t" << es->parameters.get<double> ("last_nonlinear_iterate");
		// extra ones for moghadam
		if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 13)
		{
			linear_iterations_output_file << "\t" << total_gmres_linear_iterations << 
																			 "\t" << total_cg_linear_iterations;
		}

		linear_iterations_output_file << std::endl;
	}
}



void NavierStokesCoupled::output_logging()
{

	std::ostringstream logging_data_file_name;

	logging_data_file_name << output_folder.str() << "perf_log.dat";
	std::ofstream logging_output_file(logging_data_file_name.str().c_str());

	logging_output_file << perf_log.get_log();
	logging_output_file.close();

	// do in easy format
	output_timings_to_file();
}


void NavierStokesCoupled::output_deposition_data(bool header)
{

	std::cout << "yeah" << std::endl;
	if(header)
	{
		// write geometry variables later when reading meta data file
		// write boundary conditions later with Picard class
		deposition_output_file << "# Deposition results" << std::endl;
		deposition_output_file << "# generation\t\t\ttotal_efficiency\t\t\ttb_efficiency\t\t\tacinar_efficiency\t\t\ttotal_efficiency";
	
		deposition_output_file << std::endl;
	}
	else
	{
		//calculate the generational data
		std::cout << "hmmm" << std::endl;
		std::vector<double> efficiency_per_generation(num_generations[0],0.);
		std::vector<double> efficiency_per_order(num_generations[0],0.);
		double complete_efficiency = 0.;
		for(unsigned int i=0; i<airway_data.size(); i++)
		{
			//gen
			unsigned int generation = airway_data[i].get_generation();
			if(generation+1 > efficiency_per_generation.size())
				efficiency_per_generation.resize(generation+1,0.);
			efficiency_per_generation[generation] += total_efficiency[i];

			//order
			unsigned int order = airway_data[i].get_order();
			if(order > efficiency_per_order.size())
				efficiency_per_order.resize(order+1,0.);

			efficiency_per_order[order] += total_efficiency[i];

			complete_efficiency += total_efficiency[i];
		}


		std::cout << "hmmm" << std::endl;

		std::cout << "effperg.size() = " << efficiency_per_generation.size() << std::endl;
		std::cout << "alveolar_efficiency_per_generation.size() = " << alveolar_efficiency_per_generation[1].size() << std::endl;
		std::cout << "tb_efficiency_per_generation.size() = " << tb_efficiency_per_generation[1].size() << std::endl;

		for(unsigned int i=0; i < efficiency_per_generation.size(); i++)
		{
			//calc average per unit input
			double ave_tb_per_gen = 0.;
			double ave_alv_per_gen = 0.;
			double ave_per_gen = 0.;
			for(unsigned int j=0; j<tb_efficiency_per_generation.size(); j++)
			{
				ave_tb_per_gen += tb_efficiency_per_generation[j][i]/total_particles_inhaled;
				ave_alv_per_gen += alveolar_efficiency_per_generation[j][i]/total_particles_inhaled;
			}
		
			ave_per_gen = ave_alv_per_gen + ave_tb_per_gen;

			deposition_output_file <<	i;
			deposition_output_file <<	"\t\t\t" << efficiency_per_generation[i];
			deposition_output_file <<	"\t\t\t" << ave_tb_per_gen;
			deposition_output_file <<	"\t\t\t" << ave_alv_per_gen;
			deposition_output_file <<	"\t\t\t" << ave_per_gen;
			//if(i < efficiency_per_order.size())
			//{
			//	deposition_output_file <<	"\t\t\t" << efficiency_per_order[i];
			//}
			deposition_output_file << std::endl;
		}

		std::cout << "hmmm" << std::endl;

		std::cout << "total deposition = " << complete_efficiency << std::endl;

	}
}




double NavierStokesCoupled::scaled_time()
{
	if(es->parameters.get<bool>("nondimensionalised"))
		return time*time_scale_factor;
	else
		return time;
}

void NavierStokesCoupled::init_dof_variable_vectors()
{
	std::cout << "hi" << std::endl;

	
	if(sim_type == 0 || sim_type == 2 || sim_type == 3 || sim_type == 4)
		init_dof_variable_vector(system_3d, dof_variable_type_3d);
	

	std::cout << "hi" << std::endl;
	if(sim_type == 1 || sim_type == 2 || sim_type == 3 || sim_type == 4)
	{
		init_dof_variable_vector(system_1d, dof_variable_type_1d);
	}

	std::cout << "hi" << std::endl;
	if(sim_type == 5)
	{
		init_dof_variable_vector(system_coupled, dof_variable_type_coupled);

	}

	std::cout << "hi" << std::endl;

	// dof variable vectors for output and reading output
	if(sim_3d)
	{
		init_dof_variable_vector(system_3d_output, dof_variable_type_3d_output);
		init_dof_variable_vector(extra_3d_data_system, dof_variable_type_3d_extra);
	}

	std::cout << "hi" << std::endl;
	if(sim_1d)
	{
		init_dof_variable_vector(system_1d_output, dof_variable_type_1d_output);
		if(es->parameters.get<unsigned int>("acinar_model") == 1)
		{
			init_dof_variable_vector(system_acinar_output, dof_variable_type_acinar_output);
		}

	}
	std::cout << "hi" << std::endl;
}

void NavierStokesCoupled::init_dof_variable_vector(System * system, std::vector<int>& _dof_variable_type)
{

	std::cout << "initing dof variable vector" << std::endl;	

	// DofMap things
  std::vector<dof_id_type> dof_indices_var;

	std::cout << "hi2" << std::endl;
  const DofMap & dof_map = system->get_dof_map();
	std::cout << "hi2" << std::endl;
  const int nv_sys = system->n_vars();
	std::cout << "hi2" << std::endl;
	_dof_variable_type.resize(0);
	std::cout << "hi2" << std::endl;
	_dof_variable_type.resize(system->solution->size(),0);
	std::cout << "hi2" << std::endl;

	for (int var=0; var<nv_sys; var++)
	{
		MeshBase::element_iterator       it       = system->get_mesh().active_elements_begin();
		const MeshBase::element_iterator end_elem = system->get_mesh().active_elements_end();
		for ( ; it != end_elem; ++it)
		{
			const Elem* elem = *it;

			dof_map.dof_indices (elem, dof_indices_var, var);

			for(unsigned int i=0; i<dof_indices_var.size(); i++)
			{
				_dof_variable_type[dof_indices_var[i]] = var;
			}

		}
	}
	std::cout << "hi2" << std::endl;



}

/*
// only for 3D at the moment, otherwise need separate vectors
void NavierStokesCoupled::setup_results_vector(std::string results_folder, 
				std::vector<NumericVector<Number>* >& results_vector)
{
	
	//namespace fs = boost::filesystem;
	//fs::path Path(results_folder);

	//fs::directory_iterator end_iter; // Default constructor for an iterator is the end iterator
	


	//unsigned int len = 0;//strlen(results_folder);
  DIR *Dir;
	class dirent *dp;
	std::string prefix = "out_3D_es_";
	std::vector<std::string> filenames;

	
	std::cout << "searching files in folder: " << results_folder << std::endl;

	if( (Dir = opendir(results_folder.c_str())) )
	{	
  	while ((dp = readdir(Dir)) != NULL)
		{	
			std::string filename(dp->d_name);
		 	if (!strcmp(filename.substr(0,10).c_str(), prefix.c_str())) 
			{
				filenames.push_back(filename);
				//std::cout << filename << std::endl;
		 	}
		}
  	closedir(Dir);
	}

	std::sort(filenames.begin(), filenames.end());
	for(unsigned int i=0; i< filenames.size(); i++)
	{
		std::cout << filenames[i] << std::endl;

		std::string full_path = results_folder + "/" + filenames[i];

		std::cout << full_path << std::endl;;

		AutoPtr<MeshBase> meshptr = mesh.clone();
		MeshBase &temp_mesh = *meshptr;
		EquationSystems temp_es (temp_mesh);
		// temp_es.read("/home/james/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/osc_pipe_study/nu1.0_dt0.1_v10/out_3D_es_0000.xda", libMeshEnums::READ);
		//hmmm this doesn't seem to work when have indirect path
		temp_es.read(full_path, libMeshEnums::READ);

		if(sim_3d && sim_type == 5)
		{
			//AutoPtr<NumericVector<Number> > temp_solution;
			// *temp_solution = *temp_es.get_system<TransientLinearImplicitSystem> ("ns3d1d").solution;
			// but won't this memory be deleted when out of this function?????
			NumericVector<Number>* temp_pointer = &temp_es.get_system<TransientLinearImplicitSystem> ("ns3d1d").solution;
			AutoPtr<NumericVector<Number> > temp_auto_pointer = temp_es.get_system<TransientLinearImplicitSystem> ("ns3d1d").solution->clone();
			

			results_vector.push_back(temp_pointer);
		}
		else if(sim_3d && sim_type != 5)
		{
			NumericVector<Number>* temp_pointer = &temp_es.get_system<TransientLinearImplicitSystem> ("ns3d").solution;
			results_vector.push_back(temp_pointer);
		}
	}
}

*/


// here we setup the variable scalings that we want to use. 
// the simulation is done in nondimensional variables and output in SI units. (not the particle deposition)
void NavierStokesCoupled::setup_variable_scalings_3D()
{

	std::cout << "Setting up 3D variable scalings." << std::endl;

	double velocity_scale = 1.0;
	double pressure_scale = 1.0;

	if(es->parameters.get<bool>("nondimensionalised") && !es->parameters.get<bool>("output_nondim"))
	{
		velocity_scale = es->parameters.get<double>("velocity_scale");
		pressure_scale = es->parameters.get<double>("density") *	pow(es->parameters.get<double>("velocity_scale"),2.0);

	}


	var_scalings_3D.resize(0);
	// these indices are based on the order in which they are given to exodus
	if(!threed)
	{
		var_scalings_3D.push_back(velocity_scale);	// u
		var_scalings_3D.push_back(velocity_scale);	// v
		var_scalings_3D.push_back(pressure_scale);	// p
	}
	else
	{
		var_scalings_3D.push_back(velocity_scale);	// u
		var_scalings_3D.push_back(velocity_scale);	// v
		var_scalings_3D.push_back(velocity_scale);	// w
		var_scalings_3D.push_back(pressure_scale);	// p
	}

	if(es->parameters.get<bool>("output_total_pressure"))
	{
		var_scalings_3D.push_back(pressure_scale);	// total_p
	}


	if(es->parameters.get<bool>("optimisation_stabilised"))
	{
		if(!threed)
		{
			var_scalings_3D.push_back(velocity_scale);	// u_adj
			var_scalings_3D.push_back(velocity_scale);	// v_adj
			var_scalings_3D.push_back(pressure_scale);	// p_adj
		}
		else
		{		
			var_scalings_3D.push_back(velocity_scale);	// u_adj
			var_scalings_3D.push_back(velocity_scale);	// v_adj
			var_scalings_3D.push_back(velocity_scale);	// w_adj
			var_scalings_3D.push_back(pressure_scale);	// p_adj
		}
	}


}


// here we setup the variable scalings that we want to use. 
// the simulation is done in nondimensional variables and output in SI units. (not the particle deposition)
void NavierStokesCoupled::setup_variable_scalings_1D()
{

	std::cout << "Setting up 1D variable scalings." << std::endl;

	double flow_scale = 1.0;
	double mean_pressure_scale = 1.0;

	if(es->parameters.get<bool>("nondimensionalised") && !es->parameters.get<bool>("output_nondim"))
	{
		flow_scale = es->parameters.get<double>("velocity_scale") * 
										pow(es->parameters.get<double>("length_scale"),2.0);
		mean_pressure_scale = es->parameters.get<double>("density") *
													pow(es->parameters.get<double>("velocity_scale"),2.0);
	}

	unsigned int num_airflow_vars = 7;
	var_scalings_1D.resize(7);
	var_scalings_1D[0] = mean_pressure_scale;	// P
	var_scalings_1D[1] = flow_scale;	// Q
	var_scalings_1D[2] = 1.0;	// radius ->dimensionalised already
	var_scalings_1D[3] = 1.0;	// radius vis ->dimensionalised already
	var_scalings_1D[4] = 1.0;	// poiseuille
	var_scalings_1D[5] = 1.0;	// proc_id
	var_scalings_1D[6] = 1.0;	// resistance ->dimensionalised already

	unsigned int num_part_1_vars = 0;
	if(particle_deposition == 3 ||  particle_deposition == 4 || particle_deposition == 5 ||  particle_deposition == 6)
	{
		num_part_1_vars = 7;
		var_scalings_1D.resize(num_airflow_vars + num_part_1_vars);
		var_scalings_1D[num_airflow_vars+0] = 1.0;	// p_airway
		var_scalings_1D[num_airflow_vars+1] = 1.0;	// p_total
		var_scalings_1D[num_airflow_vars+2] = 1.0;	// p_airway_sed
		var_scalings_1D[num_airflow_vars+3] = 1.0;	// p_total_sed
		var_scalings_1D[num_airflow_vars+4] = 1.0;	// p_airway_imp
		var_scalings_1D[num_airflow_vars+5] = 1.0;	// p_total_imp
		var_scalings_1D[num_airflow_vars+6] = 1.0;	// vel
	}


	unsigned int num_part_2_vars = 0;
	if(particle_deposition == 5 ||  particle_deposition == 6)
	{
		num_part_2_vars = 4;
		var_scalings_1D.resize(num_airflow_vars + num_part_1_vars + num_part_2_vars);
		var_scalings_1D[num_airflow_vars + num_part_1_vars + 0] = 1.0;	// particle_fraction
		var_scalings_1D[num_airflow_vars + num_part_1_vars + 1] = 1.0;	// terminal_exit_fraction
		var_scalings_1D[num_airflow_vars + num_part_1_vars + 2] = 1.0;	// Q_1
		var_scalings_1D[num_airflow_vars + num_part_1_vars + 3] = 1.0;	// Q_2
	}

}




// here we setup the variable scalings that we want to use. 
// the simulation is done in nondimensional variables and output in SI units. (not the particle deposition)
void NavierStokesCoupled::setup_variable_scalings_acinar()
{

	std::cout << "Setting up acinar variable scalings." << std::endl;

	double volume_scale = 1.0;

	if(es->parameters.get<bool>("nondimensionalised") && !es->parameters.get<bool>("output_nondim"))
	{
		volume_scale = pow(es->parameters.get<double>("length_scale"),3.0);
	}

	var_scalings_1D.resize(1);
	var_scalings_1D[0] = volume_scale;	// V


}

// clean up all the matrices and stuff that may not have been deleted
void NavierStokesCoupled::petsc_clean_up()
{
	std::cout << "In PETSc clean up." << std::endl;


	if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 9)
	{
		delete velocity_mass_matrix;
	}

	if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 8
			|| es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 9
			|| es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 10
			|| es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 11)
	{
		delete velocity_matrix;
		MatDestroy(&schur_complement_approx);
	
		// if we did a bfbt schur stokes, but didn't do bfbt for the rest
		if(es->parameters.get<unsigned int> ("preconditioner_type_schur_stokes") == 2
			&& es->parameters.get<unsigned int>("preconditioner_type_3d") != 9)
			delete velocity_mass_matrix;
	
	}
}


void NavierStokesCoupled::set_reynolds_number()
{

	if(!es->parameters.get<bool>("nondimensionalised"))
	{

		std::cout << "\n Setting Reynolds number" << std::endl;
		double characteristic_velocity = 0.;
		double diameter = es->parameters.get<double> ("characteristic_length");
		double radius = 0.5*diameter;
		double area = M_PI*pow(radius,2.0);

		// for 2D simulations A = D
		if(!threed)
			area = diameter;


		if(sim_type == 1)
		{
			// starting at 0D simulation
			// Re = rho*Q*D/(mu*A)

			// assume 3D
			characteristic_velocity = 8./M_PI/pow(diameter,2.0)*es->parameters.get<double> ("flow_mag_1d");
				
		}
		else
		{

			// need the velocity magnitude
			if(es->parameters.get<bool> ("prescribed_flow"))
			{
				//if 2D then flow rate = 2/3*V_max*D and not kinda assuming 2D flow
				if(!threed && !es->parameters.get<bool> ("twod_oned_tree"))
					characteristic_velocity = 3./2./diameter*es->parameters.get<double> ("flow_mag_3d");
				else
					characteristic_velocity = 8./M_PI/pow(diameter,2.0)*es->parameters.get<double> ("flow_mag_3d");
					
			}
			else
			{		
				characteristic_velocity = es->parameters.get<double> ("velocity_mag_3d");
			}

		}


		std::cout << "density = " << es->parameters.get<double> ("density") << std::endl;
		std::cout << "viscosity = " << es->parameters.get<double> ("viscosity") << std::endl;
		std::cout << "characteristic_velocity = " << characteristic_velocity << std::endl;
		std::cout << "diameter = " << diameter << std::endl;
		std::cout << "area = " << area << std::endl;


		es->parameters.set<double>("reynolds_number") = 
				es->parameters.get<double> ("density") * characteristic_velocity * diameter
				/ es->parameters.set<double> ("viscosity");


		// setup vector for reynolds numbers/viscosity if necessary

		// **************** SET REYNOLDS NUMBER FOR NUMERICAL CONTINUATION ************** //
		if(unsteady || !es->parameters.get<bool> ("viscosity_controlled"))
		{
			re_vec.push_back(es->parameters.get<double>("reynolds_number"));
		}
		else
		{
		
			double current_re = 1.0;
			current_re = es->parameters.get<double> ("numerical_continuation_starting_reynolds_number");
			std::cout << "Ramping up Reynolds number as:" << std::endl;
			std::cout << "\t";
			while(current_re < es->parameters.get<double>("reynolds_number"))
			{
				std::cout << current_re << ", ";
				re_vec.push_back(current_re);
				current_re = current_re*es->parameters.get<double> ("numerical_continuation_scaling_factor");
			}
			re_vec.push_back(es->parameters.get<double>("reynolds_number"));
	
			es->parameters.set<unsigned int> ("num_continuation_iteration")   = 0;
		}

		es->parameters.set<Real> ("reynolds_number")   = re_vec[0];


		std::cout << "Reynolds number = " << es->parameters.get<double>("reynolds_number") << std::endl;
	}
	else	//
	{
		re_vec.push_back(es->parameters.get<double>("reynolds_number"));
		std::cout << "Reynolds number = " << es->parameters.get<double>("reynolds_number") << std::endl;
	}
}





// read in the solutions at time_step and time
void NavierStokesCoupled::read_old_solutions(unsigned int read_time_step, double read_time)
{


		//std::cout << "restart_folder = " << restart_folder.str() << std::endl;

		if(sim_3d)
		{
			std::cout << "\nReading 3d solution." << std::endl;
			std::ostringstream file_name_3d;
			file_name_3d << restart_folder.str() << "out_3D";
			read_old_solution_from_file(&*es_3d,file_name_3d,read_time_step);
		}

		if(sim_1d)
		{
			std::cout << "\nReading 1d solution." << std::endl;
			std::ostringstream file_name_1d;
			file_name_1d << restart_folder.str() << "out_1D";
			read_old_solution_from_file(&*es_1d,file_name_1d,read_time_step);
		}


		// minimal time spent from here

		// convert 1d solution back to nodal
		if(sim_1d)
		{
			std::cout << "converting output solution from monomial to nodal in program." << std::endl;
			convert_1d_solution_monomial_to_nodal();
		}

		// now we have read the old solution into the output eq systems,
		// we copy it to the actual eq system that will be used
		if(sim_3d)
			copy_3d_solution_output_to_program();

		if(sim_1d)
			copy_1d_solution_output_to_program();

}






void NavierStokesCoupled::read_old_solution_from_file(EquationSystems* _es, std::ostringstream& file_name, unsigned int read_time_step)
{

	std::ostringstream file_name_es;

	std::cout << "read_time_step = " << read_time_step << std::endl;
	
	file_name_es << file_name.str() << "_es_" << std::setw(4) 
		<< std::setfill('0') << read_time_step << ".xda";

	std::cout << "Reading in data for restart from file:" << std::endl;
	std::cout << "\t" << file_name_es.str() << std::endl;

	unsigned int read_flags = (EquationSystems::READ_DATA); //(READ_HEADER | READ_DATA | READ_ADDITIONAL_DATA);
	std::cout << "hello" << std::endl;
	_es->read(file_name_es.str(), libMeshEnums::READ,read_flags);
	std::cout << "Done reading in data." << std::endl;

}





// copy the current solution to the previous solution
// copy is not necessarily the most efficient but hey
void NavierStokesCoupled::copy_back_1d_solutions()
{
	std::cout << "Copying back solution" << std::endl;


	std::cout << "1D flow" << std::endl;
	TransientExplicitSystem * system_airflow;
	// Get a reference to the airflow system object that we get the flow rate from.
	system_airflow = &es_1d->get_system<TransientExplicitSystem> ("ns1d_output");

	double flow_scale = 1.0;

	if(es->parameters.get<bool>("nondimensionalised") && !es->parameters.get<bool>("output_nondim"))
	{
		flow_scale = es->parameters.get<double>("velocity_scale") * 
										pow(es->parameters.get<double>("length_scale"),2.0);
	}

	//std::cout << "boo" << std::endl;
	// some stuff for inside the loop
	const DofMap & dof_map_airflow = system_airflow->get_dof_map();
	std::vector<dof_id_type> dof_indices_airflow;
	std::vector<dof_id_type> dof_indices_q;
	std::vector<dof_id_type> dof_indices_q_1;
	std::vector<dof_id_type> dof_indices_q_2;
	const unsigned int q_var = system_airflow->variable_number ("Q");
	const unsigned int q_1_var = system_airflow->variable_number ("Q_1");
	const unsigned int q_2_var = system_airflow->variable_number ("Q_2");


	// Get a constant reference to the mesh object.
	const MeshBase& mesh = es_1d->get_mesh();

	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	//std::cout << "boo" << std::endl;

	for ( ; el != end_el; ++el)
	{
		const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
		
			// for each element
			// - get the flow rate for the airway from solution vector
			// - calculate the deposition due to each type
			// - add them up
			// - put it in airway_data
			// - put it in a solution vector
			// note, everything is in SI units



			const int current_el_idx = elem->id();
	    dof_map_airflow.dof_indices (elem, dof_indices_airflow);
	    dof_map_airflow.dof_indices (elem, dof_indices_q, q_var);
	    dof_map_airflow.dof_indices (elem, dof_indices_q_1, q_1_var);
	    dof_map_airflow.dof_indices (elem, dof_indices_q_2, q_2_var);


			for(unsigned int i=0; i<dof_indices_q.size(); i++)
			{
				double Q = system_airflow->current_solution (dof_indices_q[i]);
				double Q_1 = system_airflow->current_solution (dof_indices_q_1[i]);
				double Q_2 = system_airflow->current_solution (dof_indices_q_2[i]);
				// Q_2 -> Q_1
				system_airflow->solution->set(dof_indices_q_1[i],Q_2 / flow_scale);
				// Q -> Q_2
				system_airflow->solution->set(dof_indices_q_2[i],Q / flow_scale);

			}

		}
	}

	system_airflow->solution->close();

}






// copy the current solution to the previous solution
// copy is not necessarily the most efficient but hey
void NavierStokesCoupled::copy_back_3d_solutions()
{
	std::cout << "Copying back solution" << std::endl;


	std::cout << "3D flow" << std::endl;
	TransientExplicitSystem * system_airflow;
	// Get a reference to the airflow system object that we get the flow rate from.
	system_airflow = &es_3d->get_system<TransientExplicitSystem> ("ns3d_output");

	double velocity_scale = 1.0;
	double pressure_scale = 1.0;

	if(es->parameters.get<bool>("nondimensionalised") && !es->parameters.get<bool>("output_nondim"))
	{
		velocity_scale = es->parameters.get<double>("velocity_scale");
		pressure_scale = es->parameters.get<double>("density") *	pow(es->parameters.get<double>("velocity_scale"),2.0);

	}

	//std::cout << "boo" << std::endl;
	// some stuff for inside the loop
	const DofMap & dof_map_airflow = system_airflow->get_dof_map();
	std::vector<dof_id_type> dof_indices_airflow;
	std::vector<dof_id_type> dof_indices_u, dof_indices_u_1, dof_indices_u_2;
	std::vector<dof_id_type> dof_indices_v, dof_indices_v_1, dof_indices_v_2;
	std::vector<dof_id_type> dof_indices_w, dof_indices_w_1, dof_indices_w_2;
	const unsigned int u_var = system_airflow->variable_number ("u");
	const unsigned int u_1_var = system_airflow->variable_number ("u_1");
	const unsigned int u_2_var = system_airflow->variable_number ("u_2");
	const unsigned int v_var = system_airflow->variable_number ("v");
	const unsigned int v_1_var = system_airflow->variable_number ("v_1");
	const unsigned int v_2_var = system_airflow->variable_number ("v_2");
	unsigned int w_var,w_1_var,w_2_var;
	if(threed)
	{
		w_var = system_airflow->variable_number ("w");
		w_1_var = system_airflow->variable_number ("w_1");
		w_2_var = system_airflow->variable_number ("w_2");
	}


	// Get a constant reference to the mesh object.
	const MeshBase& mesh = es_3d->get_mesh();

	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	//std::cout << "boo" << std::endl;

	for ( ; el != end_el; ++el)
	{
		const Elem* elem = *el;
		if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
		{
		
			// for each element
			// - get the flow rate for the airway from solution vector
			// - calculate the deposition due to each type
			// - add them up
			// - put it in airway_data
			// - put it in a solution vector
			// note, everything is in SI units



			const int current_el_idx = elem->id();
	    dof_map_airflow.dof_indices (elem, dof_indices_airflow);
	    dof_map_airflow.dof_indices (elem, dof_indices_u, u_var);
	    dof_map_airflow.dof_indices (elem, dof_indices_u_1, u_1_var);
	    dof_map_airflow.dof_indices (elem, dof_indices_u_2, u_2_var);
	    dof_map_airflow.dof_indices (elem, dof_indices_v, v_var);
	    dof_map_airflow.dof_indices (elem, dof_indices_v_1, v_1_var);
	    dof_map_airflow.dof_indices (elem, dof_indices_v_2, v_2_var);
			if(threed)
			{
			  dof_map_airflow.dof_indices (elem, dof_indices_w, w_var);
			  dof_map_airflow.dof_indices (elem, dof_indices_w_1, w_1_var);
			  dof_map_airflow.dof_indices (elem, dof_indices_w_2, w_2_var);
			}


			for(unsigned int i=0; i<dof_indices_u.size(); i++)
			{
				double u = system_airflow->current_solution (dof_indices_u[i]);
				double u_1 = system_airflow->current_solution (dof_indices_u_1[i]);
				double u_2 = system_airflow->current_solution (dof_indices_u_2[i]);
				double v = system_airflow->current_solution (dof_indices_v[i]);
				double v_1 = system_airflow->current_solution (dof_indices_v_1[i]);
				double v_2 = system_airflow->current_solution (dof_indices_v_2[i]);
				double w,w_1,w_2;
				if(threed)
				{
					w = system_airflow->current_solution (dof_indices_w[i]);
					w_1 = system_airflow->current_solution (dof_indices_w_1[i]);
					w_2 = system_airflow->current_solution (dof_indices_w_2[i]);
				}

				// u_2 -> u_1
				system_airflow->solution->set(dof_indices_u_1[i],u_2 / velocity_scale);
				system_airflow->solution->set(dof_indices_v_1[i],v_2 / velocity_scale);
				if(threed)
					system_airflow->solution->set(dof_indices_w_1[i],w_2 / velocity_scale);

				// u -> u_2
				system_airflow->solution->set(dof_indices_u_2[i],u / velocity_scale);
				system_airflow->solution->set(dof_indices_v_2[i],v / velocity_scale);
				if(threed)
					system_airflow->solution->set(dof_indices_w_2[i],w / velocity_scale);


			}

		}
	}

	system_airflow->solution->close();


}




// copy the current solution to the previous solution
// copy is not necessarily the most efficient but hey
void NavierStokesCoupled::calculate_local_particle_1d_flow_rate()
{
	std::cout << "Calculating local particle flow rate" << std::endl;


	std::cout << "1D flow" << std::endl;
	TransientExplicitSystem * system_airflow;
	// Get a reference to the airflow system object that we get the flow rate from.
	system_airflow = &es_1d->get_system<TransientExplicitSystem> ("ns1d_output");


	//std::cout << "boo" << std::endl;
	// some stuff for inside the loop
	const DofMap & dof_map_airflow = system_airflow->get_dof_map();
	std::vector<dof_id_type> dof_indices_airflow;
	std::vector<dof_id_type> dof_indices_q;
	std::vector<dof_id_type> dof_indices_q_1;
	std::vector<dof_id_type> dof_indices_q_2;
	const unsigned int q_var = system_airflow->variable_number ("Q");
	const unsigned int q_1_var = system_airflow->variable_number ("Q_1");
	const unsigned int q_2_var = system_airflow->variable_number ("Q_2");


	// Get a constant reference to the mesh object.
	const MeshBase& mesh = es_1d->get_mesh();

	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	//std::cout << "boo" << std::endl;

	for ( ; el != end_el; ++el)
	{
		const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
		
			// for each element
			// - get the flow rate for the airway from solution vector
			// - calculate the deposition due to each type
			// - add them up
			// - put it in airway_data
			// - put it in a solution vector
			// note, everything is in SI units



			const int current_el_idx = elem->id();
	    dof_map_airflow.dof_indices (elem, dof_indices_airflow);
	    dof_map_airflow.dof_indices (elem, dof_indices_q, q_var);
	    dof_map_airflow.dof_indices (elem, dof_indices_q_1, q_1_var);
	    dof_map_airflow.dof_indices (elem, dof_indices_q_2, q_2_var);

			// Q = Q_1 + (time - real_time_1)/(real_time_2 - real_time_1) * (Q_2 - Q_1)
			double fraction_of_fluid_dt = 0.;
			if(fabs(read_time_2 - read_time_1) > 1e-10)
				fraction_of_fluid_dt = (time - read_time_1) / (read_time_2 - read_time_1);
			else
				fraction_of_fluid_dt = 0.;


			for(unsigned int i=0; i<dof_indices_q.size(); i++)
			{
				double Q_1 = system_airflow->current_solution (dof_indices_q_1[i]);
				double Q_2 = system_airflow->current_solution (dof_indices_q_2[i]);

	
				// needs to be based on the new values, which won't have changed 
				double Q = Q_1 + fraction_of_fluid_dt * (Q_2 - Q_1);
				system_airflow->solution->set(dof_indices_q[i],Q);
			}

		}
	}

	system_airflow->solution->close();
	//system_airflow->update();	// don't need this




}






// copy the current solution to the previous solution
// copy is not necessarily the most efficient but hey
void NavierStokesCoupled::calculate_local_particle_3d_flow_rate()
{
	std::cout << "Calculating local particle flow rate" << std::endl;

	std::cout << "3D flow" << std::endl;

	if(true)
	{
		TransientExplicitSystem * system_airflow;
		// Get a reference to the airflow system object that we get the flow rate from.
		system_airflow = &es_3d->get_system<TransientExplicitSystem> ("ns3d_output");


		//std::cout << "boo" << std::endl;
		// some stuff for inside the loop
		const DofMap & dof_map_airflow = system_airflow->get_dof_map();
		std::vector<dof_id_type> dof_indices_airflow;
		std::vector<dof_id_type> dof_indices_u, dof_indices_u_1, dof_indices_u_2;
		std::vector<dof_id_type> dof_indices_v, dof_indices_v_1, dof_indices_v_2;
		std::vector<dof_id_type> dof_indices_w, dof_indices_w_1, dof_indices_w_2;
		const unsigned int u_var = system_airflow->variable_number ("u");
		const unsigned int u_1_var = system_airflow->variable_number ("u_1");
		const unsigned int u_2_var = system_airflow->variable_number ("u_2");
		const unsigned int v_var = system_airflow->variable_number ("v");
		const unsigned int v_1_var = system_airflow->variable_number ("v_1");
		const unsigned int v_2_var = system_airflow->variable_number ("v_2");
		unsigned int w_var, w_1_var, w_2_var;
		if(threed)
		{
			w_var = system_airflow->variable_number ("w");
			w_1_var = system_airflow->variable_number ("w_1");
			w_2_var = system_airflow->variable_number ("w_2");
		}

		// u = u_1 + (time - real_time_1)/(real_time_2 - real_time_1) * (u_2 - u_1)
		double fraction_of_fluid_dt = 0.;
		if(fabs(read_time_2 - read_time_1) > 1e-10)
			fraction_of_fluid_dt = (time - read_time_1) / (read_time_2 - read_time_1);
		else
			fraction_of_fluid_dt = 0.;

		// Get a constant reference to the mesh object.
		const MeshBase& mesh = es_3d->get_mesh();

		MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
		const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

		//std::cout << "boo" << std::endl;


		for ( ; el != end_el; ++el)
		{
			const Elem* elem = *el;
			if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
			{
		
				// for each element
				// - get the flow rate for the airway from solution vector
				// - calculate the deposition due to each type
				// - add them up
				// - put it in airway_data
				// - put it in a solution vector
				// note, everything is in SI units



				const int current_el_idx = elem->id();
			  dof_map_airflow.dof_indices (elem, dof_indices_airflow);
			  dof_map_airflow.dof_indices (elem, dof_indices_u, u_var);
			  dof_map_airflow.dof_indices (elem, dof_indices_u_1, u_1_var);
			  dof_map_airflow.dof_indices (elem, dof_indices_u_2, u_2_var);
			  dof_map_airflow.dof_indices (elem, dof_indices_v, v_var);
			  dof_map_airflow.dof_indices (elem, dof_indices_v_1, v_1_var);
			  dof_map_airflow.dof_indices (elem, dof_indices_v_2, v_2_var);
				if(threed)
				{
					dof_map_airflow.dof_indices (elem, dof_indices_w, w_var);
					dof_map_airflow.dof_indices (elem, dof_indices_w_1, w_1_var);
					dof_map_airflow.dof_indices (elem, dof_indices_w_2, w_2_var);
				}


				// u,v,w will have the same size
				for(unsigned int i=0; i<dof_indices_u.size(); i++)
				{
					double u_1 = system_airflow->current_solution (dof_indices_u_1[i]);
					double u_2 = system_airflow->current_solution (dof_indices_u_2[i]);

					double v_1 = system_airflow->current_solution (dof_indices_v_1[i]);
					double v_2 = system_airflow->current_solution (dof_indices_v_2[i]);

					double w_1,w_2;
					if(threed)
					{
						double w_1 = system_airflow->current_solution (dof_indices_w_1[i]);
						double w_2 = system_airflow->current_solution (dof_indices_w_2[i]);
					}

	
					// needs to be based on the new values, which won't have changed 
					double u = u_1 + fraction_of_fluid_dt * (u_2 - u_1);
					double v = v_1 + fraction_of_fluid_dt * (v_2 - v_1);
					double w;
					if(threed)
						w = w_1 + fraction_of_fluid_dt * (w_2 - w_1);

					system_airflow->solution->set(dof_indices_u[i],u);
					system_airflow->solution->set(dof_indices_v[i],v);
					if(threed)
						system_airflow->solution->set(dof_indices_w[i],w);
				}

			}
		}

		system_airflow->solution->close();
		//system_airflow->update();	// don't need this
	}
	else
	{
		// we want to manipulate vector variables without looping over elements

		// u = u_1 + (time - real_time_1)/(real_time_2 - real_time_1) * (u_2 - u_1)
		double fraction_of_fluid_dt = 0.;
		if(fabs(read_time_2 - read_time_1) > 1e-10)
			fraction_of_fluid_dt = (time - read_time_1) / (read_time_2 - read_time_1);
		else
			fraction_of_fluid_dt = 0.;

	}
	


}


// copy the current solution to the previous solution
// copy is not necessarily the most efficient but hey
std::string NavierStokesCoupled::get_converged_reason_string(int converged_reason)
{
	std::string converged_reason_string = "";

	// switch through all possible reasons
	switch (converged_reason)
	{
	case 1:
		{
			converged_reason_string = "KSP_CONVERGED_RTOL_NORMAL";
			break;
		}
	case 2:
		{
			converged_reason_string = "KSP_CONVERGED_RTOL";
			break;
		}
	case 3:
		{
			converged_reason_string = "KSP_CONVERGED_ATOL";
			break;
		}
	case 4:
		{
			converged_reason_string = "KSP_CONVERGED_ITS";
			break;
		}
	case -3:
		{
			converged_reason_string = "KSP_DIVERGED_ITS";
			break;
		}
	case -4:
		{
			converged_reason_string = "KSP_DIVERGED_DTOL";
			break;
		}
	case -5:
		{
			converged_reason_string = "KSP_DIVERGED_BREAKDOWN";
			break;
		}
	case -7:
		{
			converged_reason_string = "KSP_DIVERGED_NONSYMMETRIC";
			break;
		}
	case -8:
		{
			converged_reason_string = "KSP_DIVERGED_INDEFINITE_PC";
			break;
		}

	default:
	  {
	    std::cout << "ERROR: Unrecognized Petsc converged reason." << std::endl;
			break;
	  }
	}

	return converged_reason_string;
}



void NavierStokesCoupled::output_timings_to_file()
{
	std::ostringstream timings_file_name;
	std::ofstream timings_output_file;

	timings_file_name << output_folder.str() << "timings.dat";
	timings_output_file.open(timings_file_name.str().c_str());

	timings_output_file << "# Output file with timings" << std::endl;
	timings_output_file << "# Simple format for ease of parsing" << std::endl;
	timings_output_file << "# section_name\tcount\ttime\tpercent" << std::endl;

	double total_time = perf_log.get_elapsed_time();

	// vector of sections to output
	std::vector<std::string> section_labels;
	section_labels.push_back("assembly");
	section_labels.push_back("assembly_1d");
	section_labels.push_back("output");
	section_labels.push_back("setup");
	section_labels.push_back("solve");
	section_labels.push_back("solve_1d");
	section_labels.push_back("solve_setup");
	section_labels.push_back("solve_setup_1d");

	for(unsigned int i=0; i<section_labels.size(); i++)
	{
		PerfData section_data = perf_log.get_perf_data(section_labels[i]);

		// data to output
		unsigned int count = section_data.count;
		double section_time = section_data.tot_time;
		double percent = section_time/total_time;

		// output data
		timings_output_file << section_labels[i] << 
													"\t" << count << 
													"\t" << section_time << 
													"\t" << percent << std::endl;
		
		
	}

	// totals
	timings_output_file << "total" << 
												"\t" << 1 << 
												"\t" << total_time << 
												"\t" << 100 << std::endl;
	
	//close file
	timings_output_file.close();

}



