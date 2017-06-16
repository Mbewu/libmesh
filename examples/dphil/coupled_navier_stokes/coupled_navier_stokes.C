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
		mesh_data(mesh),
		mesh_refinement (mesh),
		es(AutoPtr<EquationSystems>(new EquationSystems(mesh))),
		time(0.),
	  	reduce_dt(false),
	  	increase_dt(false),
	  	refine_mesh(false),
		steps_since_last_dt_change(0),
		restart(false),
		perf_log("Main Program"),
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
		total_max_iterations(0),

		particle_deposition(false),
		shell_pc_created(false),
		mono_shell_pc_created(false),
		first_3d_write(true),
		first_1d_write(true),
		init_names_done(false),
		ic_set(true)
{

	perf_log.push("misc");
	perf_log.push("setup");
	PerfLog perf_log_setup("Setup");

	std::cout << "Starting simulation." << std::endl;

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
	// if 3D particle deposition then read in the velocity timesteps
	if(particle_deposition == 1)
		read_timesteps();

	std::cout << "hello" << std::endl;
	//Airway new_airway;
	//new_airway.set_generation(5);
	//std::cout << "generation of new airway = " << new_airway.get_generation() << std::endl;
	//std::exit(0);



/*
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
	


	// ************ SETUP SYSTEM CRAP ************ //

	if(sim_3d && sim_type != 5)
	{
	  system_3d = &es->add_system<TransientLinearImplicitSystem>("ns3d");
	  system_neumann = &es->add_system<TransientLinearImplicitSystem>("Neumann-Variable");
		extra_3d_data_system = &es->add_system<ExplicitSystem> ("extra_3d_data");
	}

	if(sim_1d && sim_type != 5)
	{
		system_1d = &es->add_system<TransientLinearImplicitSystem> ("ns1d");
		extra_1d_data_system = &es->add_system<ExplicitSystem> ("extra_1d_data");

		if(particle_deposition == 2)
			particle_deposition_system_1d = &es->add_system<ExplicitSystem> ("Particle-Deposition-1D");

	}

	if(sim_type == 5)
	{
		system_coupled = &es->add_system<TransientLinearImplicitSystem> ("ns3d1d");
	  system_neumann = &es->add_system<TransientLinearImplicitSystem>("Neumann-Variable");
		extra_1d_data_system = &es->add_system<ExplicitSystem> ("extra_1d_data");
		extra_3d_data_system = &es->add_system<ExplicitSystem> ("extra_3d_data");
	}

	mesh.partitioner()->set_custom_partitioning(es->parameters.set<unsigned int>("custom_partitioning"));

	// ************** SET UP 3D MESH **************** //
	if(sim_3d)
	{
		perf_log_setup.push("setup_3d_mesh");
		setup_3d_mesh(&*es,mesh);
		perf_log_setup.pop("setup_3d_mesh");
	// in case not set already
	es->parameters.set<unsigned int>("n_initial_3d_elem") = mesh.n_elem();

	}


	// ******************** SETUP 1D MESH ********************* //
	if(sim_1d)
	{
		perf_log_setup.push("setup_1d_mesh");
		setup_1d_mesh();
		perf_log_setup.pop("setup_1d_mesh");

	}

	// ************* SET UP 3D SYSTEM **************** //
	if(sim_3d)
	{
		perf_log_setup.push("setup_3d_system");
		if(!es->parameters.get<bool>("optimisation_stabilised"))
			picard = AutoPtr<Picard>(new Picard(*es,surface_boundaries,subdomains_3d,es->parameters.get<unsigned int>("n_initial_3d_elem"),es->parameters.get<bool>("efficient_assembly")));
		else
			picard = AutoPtr<OptimisedStabilisedAssembler3D>(new OptimisedStabilisedAssembler3D(*es,surface_boundaries,subdomains_3d,es->parameters.get<unsigned int>("n_initial_3d_elem")));

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


		perf_log_setup.pop("setup_3d_system");

	}




	// ******************* SETUP 1D SYSTEM ******************** //
	if(sim_1d)
	{
		
		perf_log_setup.push("setup_1d_system");
		ns_assembler = AutoPtr<NavierStokesAssembler>
								(new NavierStokesAssembler(*es,airway_data,subdomains_1d,es->parameters.get<unsigned int>("n_initial_3d_elem")));

		if(sim_type == 5)
		{
			setup_1d_system(system_coupled);
			// don't give the 1d an additional prefix
			//system_1d->attach_assemble_object (*ns_assembler); //attach later
		}
		else
		{
			setup_1d_system(system_1d);

			system_1d->attach_assemble_object (*ns_assembler);

			//get the extra couplings sorted for the matrix, maybe don't need this?
			augment_sparsity = AutoPtr<AugmentSparsityOnInterface>
										(new AugmentSparsityOnInterface(*es,airway_data,subdomains_3d,subdomains_1d,es->parameters.get<unsigned int>("n_initial_3d_elem")));
			system_1d->get_dof_map().attach_extra_sparsity_object(*augment_sparsity);
		}
		
		perf_log_setup.pop("setup_1d_system");
	}

	if(sim_type == 5)
	{
		ns_assembler->set_coupled(true);
		coupled_assembler = AutoPtr<CoupledAssembler>
								(new CoupledAssembler(*picard,*ns_assembler));
		system_coupled->attach_assemble_object (*coupled_assembler);

		//get the extra couplings sorted for the matrix
		augment_sparsity = AutoPtr<AugmentSparsityOnInterface>
									(new AugmentSparsityOnInterface(*es,airway_data,subdomains_3d,subdomains_1d,es->parameters.get<unsigned int>("n_initial_3d_elem"),true));
		system_coupled->get_dof_map().attach_extra_sparsity_object(*augment_sparsity);

		picard->set_subdomains_1d(subdomains_1d);
		picard->set_airway_data(airway_data);

	}

	// after things have been setup let us set up the variable scalings for output
	if(sim_3d)
		setup_variable_scalings_3D();

	if(sim_1d)
		setup_variable_scalings_1D();
	std::cout << std::endl;
	
	//init the equation systems

	std::cout << "Init equation systems." << std::endl;
	es->init ();
	std::cout << "done initing equation systems" << std::endl;

	init_dof_variable_vectors();
	
	perf_log_setup.push("boundary vectors compute");
	if(sim_3d)
		picard->setup_flow_rate_and_mean_pressure_vectors();
	perf_log_setup.pop("boundary vectors compute");

	if(restart)
	{
		std::ostringstream file_name;
		std::ostringstream file_name_es;
		std::ostringstream file_name_mesh;

		std::cout << "restart_folder = " << restart_folder.str() << std::endl;

		//file_name << output_folder.str() << "out_3D";
		file_name << restart_folder.str() << "out_3D";
 
		// if reading from a steady simulation then need to rad timestep 1 no matter what
		if(es->parameters.get<bool> ("unsteady_from_steady"))
			file_name_es << file_name.str() << "_es_" << std::setw(4) 
			<< std::setfill('0') << 1 << ".xda";
		else
			file_name_es << file_name.str() << "_es_" << std::setw(4) 
			<< std::setfill('0') << es->parameters.get<unsigned int> ("restart_time_step") << ".xda";
		

		std::cout << "Reading in data for restart from file:" << std::endl;
		std::cout << "\t" << file_name_es.str() << std::endl;

		unsigned int read_flags = (EquationSystems::READ_DATA); //(READ_HEADER | READ_DATA | READ_ADDITIONAL_DATA);
		es->read(file_name_es.str(), libMeshEnums::READ,read_flags);
		std::cout << "Done reading in data." << std::endl;

		
	}

	std::cout << "HEY BAE" << std::endl;

	// ************* READ IN BOUNDARY CONDITIONS FOR PROBLEM ****** //
	if(sim_type == 2 && es->parameters.get<bool>("known_boundary_conditions"))
	{
		std::cout << "you suck" << std::endl;
		// read the input boundary conditions for an uncoupled simulation
		read_input_boundary_conditions();

	}


	// *************	INITIALISE SOME STUFF BEFORE THE LOOP ******* //
	if(sim_3d)
	{
		perf_log_setup.push("update_eq_systems");
		es->update();
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

		}
	}


	if(es->parameters.get<unsigned int>("prerefine"))
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


	if(sim_1d)
	{
		set_radii();
	}

	// **************************** SET THE ELEM PIDs *********************** //
	if(sim_3d)
		set_elem_proc_id_3d();
	if(sim_1d)
		set_elem_proc_id_1d();


	// ************** CALCULATE THE RESISTANCE OF THE 1D TREES ************** //
	if(es->parameters.get<bool>("moghadam_coupling") && sim_type == 2)
		calculate_1d_linear_resistance_values();


	// ************* OUTPUT INFO TO SCREEN AND TO FILE ********** //
	mesh.print_info();
	mesh.boundary_info->print_summary();
	es->print_info();


	perf_log_setup.disable_logging();
	perf_log.pop("setup");

	PerfLog perf_log_misc("Miscellaneous");

 	// ************ TIME AND NONLINEAR LOOPS ***************************** //

	es->parameters.set<unsigned int> ("t_step") = 
		es->parameters.get<unsigned int>("restart_time_step");

	if(!es->parameters.get<bool>("compare_results") && !(particle_deposition == 2))
	{

		//loop over viscositys for steady state solution num continuation
		for (unsigned int k=0; k<re_vec.size(); ++k)
		{

		 	es->parameters.set<Real> ("reynolds_number")   = re_vec[k];
		 	es->parameters.set<unsigned int> ("num_continuation_iteration")   = k;

			// *************** WRITE OUTPUT ******************** //
			std::cout << "CALCULATING FOR Re = " <<
				es->parameters.get<Real> ("reynolds_number") << std::endl;


			//intialise the particles
			if(particle_deposition)
				init_particles();

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

			}

			perf_log.push("output");
			if(!restart)
			{
				output_sim_data(true);
			}


			//parameters
			std::ostringstream parameter_output_data_file_name;

			parameter_output_data_file_name << output_folder.str() << "parameters" << t_step << ".dat";
			//output_data_file_name << "results/out_viscosity"	<< es->parameters.set<Real> ("viscosity") << ".dat";
			parameter_output_file.open(parameter_output_data_file_name.str().c_str());
			parameter_output_file << "# Parameters" << std::endl;
			es->parameters.print(parameter_output_file);
			parameter_output_file.close();

			if(particle_deposition == 1)
			{
				
				std::cout << "lday" << es->parameters.get<unsigned int> ("particle_deposition_start_time_step") << std::endl;
				if(es->parameters.get<unsigned int> ("particle_deposition_start_time_step") == 0)
				{
					output_particle_data();
					print_particle_data();
				}
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

			perf_log.pop("output");

			// only do this when on the first step of num continuation
			if(k==0)
			{
				// must zero on first time step - after boundary conditions updated
				if(sim_3d)
				{
					if(sim_type != 5 && !restart)
					{
						system_3d->solution->zero();
						system_3d->old_local_solution->zero();

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
						system_1d->update();
					}

				}
			}


			if(sim_3d && sim_type != 5)
			{
				*system_3d->old_local_solution = *system_3d->current_local_solution;
				*old_global_solution = *system_3d->solution;
				//system_3d->update();
			}

			perf_log.push("output");
			if(unsteady && !restart)
			{
				if(sim_3d) {write_3d_solution();}
				if(sim_1d) {write_1d_solution();}
			}

			if(particle_deposition)
				write_particles();
			perf_log.pop("output");

			//******* reset all the time stuff
			time = es->parameters.get<Real>("restart_time");

			t_step = es->parameters.get<unsigned int> ("t_step");
			

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
			dt = es->parameters.get<Real> ("dt");


			// ************* TIME LOOP ******************** //
			// let us do this in dimensional terms
			while (time + 1e-10 < es->parameters.get<Real> ("end_time"))
			{




				perf_log_misc.push("update times");
				// INCREMENT TIME
				++t_step;


				// choose dt
				if(particle_deposition && !es->parameters.get<bool> ("unsteady_from_steady"))
					dt = timestep_sizes[t_step - 1];	//set timestep based on file

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

				perf_log_misc.pop("update times");

				//shell_pc_created = false;
				//mono_shell_pc_created = false;

				if(!particle_deposition)
				{
					if(sim_3d)
					{
						perf_log_misc.push("reinit");
						update_3d_dirichlet_boundary_conditions();	// UPDATE TIME DEPENDENT DIRICHLET BOUNDARY CONDITIONS
						perf_log_misc.pop("reinit");
						perf_log_misc.push("calc norm");	
						if(sim_type != 5)
	 						beforebefore_norm = system_3d->solution->l2_norm();
						else
	 						beforebefore_norm = system_coupled->solution->l2_norm();
						perf_log_misc.pop("calc norm");

					}
				}

				std::cout << "l2_norm = " << beforebefore_norm << std::endl;
				
				//write_3d_solution(true);

				std::cout << "\n\n*** Solving time step " << t_step <<
				             ", time = " << time << " (" << time*time_scale_factor <<
										 "s) ***" << std::endl;





				// **************** DO PARTICLE DEPOSITION ******************** //
				if(particle_deposition)
				{
					move_particles();

					perf_log.push("output");
					print_particle_data();
					
					if (time - ((int)((time + 1e-10) /es->parameters.get<Real>("write_interval")))
										*es->parameters.get<Real>("write_interval") < dt - 1e-10)// (t_step)%write_interval == 0)
					{
						output_particle_data();
						write_particles();
					}
					perf_log.pop("output");
					
					//intialise new particles for next time step
					// must be 
					if(particle_deposition && fabs(time - es->parameters.get<Real> ("end_time")) > 1e-10
						&& es->parameters.get<unsigned int> ("particle_deposition_type") == 1)
						init_particles();

				}
				else
				{



					// ************** SOLVE 3D SYSTEM ********************* //
					if(sim_3d && sim_type != 5)
					{
						nonlinear_iteration = 0;
						unsigned int nonlinear_iteration_1d = 0;
						es->parameters.set<unsigned int> ("nonlinear_iteration") = nonlinear_iteration;
						local_linear_iterations = 0;

						double shift_value = 0.;
					
						//TEMP : this is to test how the different meshes compare
						//double before_norm = system_3d->solution->l2_norm();
						while(true)
						{
							nonlinear_iteration++;
							es->parameters.set<unsigned int> ("nonlinear_iteration") = nonlinear_iteration;

							perf_log_misc.push("calculate 1d boundary values");
							if(sim_1d)
								calculate_1d_boundary_values();
							perf_log_misc.pop("calculate 1d boundary values");

							// sim type 3 is when we tightly couple the 1d sim to the 3d sim
							if(sim_type == 0 || sim_type == 2)
							{
					
								// calculate the flux input boundary conditions for uncoupled simulation
								// note: don't actually need to do this every nonlinear iteration but okay
								if(sim_type == 2 && es->parameters.get<bool>("known_boundary_conditions"))	// get input flux from file.
								{
									std::cout << "t_step = " << t_step << std::endl;
									for(unsigned int i=0; i<input_pressure_values_3d.size(); i++)
									{
										input_pressure_values_3d[i] = all_input_pressure_values_3d[t_step - 1][i];	// t_step-1 because first timestep not output
									}
								}


								// set bcs and solve 
								std::vector<double> empty_vec;
								if(es->parameters.get<bool>("moghadam_coupling"))
									picard->init_bc(boundary_ids,pressure_values_1d,flux_values_1d,linear_resistance_values_1d,previous_flux_values_3d,previous_previous_flux_values_3d);
								else
									picard->init_bc(boundary_ids,input_pressure_values_3d,empty_vec,linear_resistance_values_1d,previous_flux_values_3d,previous_previous_flux_values_3d); 

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

								picard->init_bc(boundary_ids,pressure_values_1d,flux_values_1d,linear_resistance_values_1d,previous_flux_values_3d,previous_previous_flux_values_3d);

								if(solve_3d_system_iteration(system_3d))
									break;
				
								perf_log_misc.push("calculate 3d boundary values");
								calculate_3d_boundary_values();
								perf_log_misc.pop("calculate 3d boundary values");
								ns_assembler->init_bc(flux_values_3d);
								for(unsigned int i=0; i<flux_values_3d.size(); i++)
									std::cout << "flux_values_3d[" << i << "] = " << flux_values_3d[i] << std::endl;
								solve_1d_system();

								perf_log_misc.push("calculate 3d boundary values");
								calculate_3d_boundary_values();
								perf_log_misc.pop("calculate 3d boundary values");
								perf_log_misc.push("calculate 3d boundary values");
								calculate_1d_boundary_values();
								perf_log_misc.pop("calculate 3d boundary values");
								perf_log.push("output");
								print_flux_and_pressure();
								perf_log.pop("output");
					

							}
							else if(sim_type == 4)
							{
								//in the explicitly coupled case sim_type 4 use 1d pressure vals
								picard->init_bc(boundary_ids,pressure_values_1d,flux_values_1d,linear_resistance_values_1d,previous_flux_values_3d,previous_previous_flux_values_3d); 
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
					if(sim_type == 1 || sim_type == 2)
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
						else if(es->parameters.get<bool>("moghadam_coupling"))	// just use the same input flux everywhere
						{
							calculate_3d_boundary_values();
							input_flux_values = flux_values_3d;
						}
						else
						{

							
							double inflow = es->parameters.get<double>("time_scaling");
			
							if(es->parameters.get<unsigned int>("unsteady"))
								inflow = es->parameters.get<double>("flow_mag_1d") * es->parameters.get<double>("time_scaling") / 
												(es->parameters.get<double>("velocity_scale") * pow(es->parameters.get<double>("length_scale"),2.0));//fabs(sin(2*M_PI*time/4));
							else	//put flow_mag into the nondimensional units
								inflow = es->parameters.get<double>("flow_mag_1d") / 
												(es->parameters.get<double>("velocity_scale") * pow(es->parameters.get<double>("length_scale"),2.0)) ;	

							//input_flux_values[0] = 0.0;
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

						double old_pressure = 0.;
						double new_pressure = 0.;

						unsigned int nonlinear_iteration_1d = 0;
						while(true)
						{

							nonlinear_iteration_1d++;
							es->parameters.set<unsigned int> ("nonlinear_iteration_1d") = nonlinear_iteration_1d;
							std::cout << "Nonlinear Iteration 1D = " << es->parameters.get<unsigned int> ("nonlinear_iteration_1d") << std::endl;
							old_pressure = new_pressure;
							ns_assembler->init_bc(input_flux_values);
							std::cout << "bye" << std::endl;
							solve_1d_system();
							std::cout << "bye" << std::endl;
							perf_log_misc.push("calculate 1d boundary values");
							calculate_1d_boundary_values();
							perf_log_misc.pop("calculate 1d boundary values");
							

							// stop when pressure has converged
							new_pressure = pressure_values_1d[1];
							std::cout << "resistance_type_1d = " << es->parameters.get<unsigned int> ("resistance_type_1d") << std::endl;
							std::cout << "old pressure = " << old_pressure << std::endl;
							std::cout << "new pressure = " << new_pressure << std::endl;
							std::cout << "fabs((new_pressure - old_pressure)/new_pressure) = " << fabs((new_pressure - old_pressure)/new_pressure) << std::endl;
							if(fabs((new_pressure - old_pressure)/new_pressure) < es->parameters.get<double> ("nonlinear_tolerance_1d"))
								break;

							if(!es->parameters.get<bool> ("nonlinear_1d"))
							{
								std::cout << "linear 1d so no nonlinear iterations required" << std::endl;
								break;
							}
					
						}
					}
					else if(sim_type == 4)
					{
						calculate_3d_boundary_values();
						ns_assembler->init_bc(flux_values_3d);
						solve_1d_system();
					}










					// *********** MONOLITHIC SOLVING **************** //
					if(sim_type == 5)
					{
	
						std::cout << "hiya" << std::endl;
						local_linear_iterations = 0;
						nonlinear_iteration = 0;
						es->parameters.set<unsigned int> ("nonlinear_iteration") = nonlinear_iteration;
						while(true)
						{
							nonlinear_iteration++;
							es->parameters.set<unsigned int> ("nonlinear_iteration") = nonlinear_iteration;
							es->parameters.set<unsigned int> ("nonlinear_iteration_1d") = nonlinear_iteration;

							perf_log_misc.push("calculate 1d boundary values");
							calculate_1d_boundary_values();
							perf_log_misc.pop("calculate 1d boundary values");
							picard->init_bc(boundary_ids,pressure_values_1d,flux_values_1d,linear_resistance_values_1d,previous_flux_values_3d,previous_previous_flux_values_3d);
							std::cout << "about to calculate 3d boundary values" << std::endl;
							perf_log_misc.push("calculate 3d boundary values");
							calculate_3d_boundary_values();
							perf_log_misc.pop("calculate 3d boundary values");
							std::cout << "done about to calculating 3d boundary values" << std::endl;
							ns_assembler->init_bc(flux_values_3d);

							if(solve_3d_system_iteration(system_coupled))
								break;
			
			
							//here we wanna output the boundary values that have been calculated and see if they are indeed converging
							if(exit_program)
								break;
						} // end nonlinear loop

						es->parameters.set<unsigned int> ("num_newton_steps") = nonlinear_iteration;
					}

					if(exit_program)
						break;



					// ************ WRITE OUTPUT *************** //		
			
					perf_log.push("output");

					std::cout << "yeah" << std::endl;						
					if((!reduce_dt || !unsteady) && !refine_mesh)
					{
						std::cout << "yeah" << std::endl;
						output_sim_data(false);
						std::cout << "yeah" << std::endl;

						if (time - ((int)((time + 1e-10) /es->parameters.get<Real>("write_interval")))
										*es->parameters.get<Real>("write_interval") < dt - 1e-10)// (t_step)%write_interval == 0)
						{
							if(sim_3d) {write_3d_solution();}
							if(sim_1d) {write_1d_solution();}
						}
					std::cout << "yeah" << std::endl;
					}
			
					perf_log.pop("output");

					//exit(0);
		
				}

				std::cout << "kk" << std::endl;

				perf_log_misc.push("end time step admin");
				end_timestep_admin();		//decide on what to do next and get vectors ready
				perf_log_misc.pop("end time step admin");

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

			output_file.close();
			eigenvalues_file.close();
			if(particle_deposition)
				particle_output_file.close();

			if(es->parameters.get<bool> ("output_linear_iteration_count"))
 				linear_iterations_output_file.close();

		} //end viscosity loop


		// petsc clean up if we have actually done some simulation, i.e. time > 0
		if(time > 1e-10 && es->parameters.get<unsigned int> ("particle_deposition") == 0)
			petsc_clean_up();

		output_logging();


		//now compare to exact solution if we want to
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



		}

	}
	else if (es->parameters.get<bool>("compare_results"))	//do the compare results stuff
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
	}
	else if (particle_deposition == 2)
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

	}
	else
	{
		std::cout << "error program type not specified correctly" << std::endl;
		std::cout << "not compare results, or particle deposition 2 or none of the above (IMPOSSIBLE) .. exiting" << std::endl;
		std::exit(0);
		
	}


	
	std::cout << "AVERAGE LINEAR ITERATION COUNT = " << (double)total_linear_iterations / (double)total_nonlinear_iterations;
	std::cout << "tot lin iterations = " << total_linear_iterations << std::endl;
	std::cout << "tot nonlin iterations = " << total_nonlinear_iterations << std::endl;

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
	set_unsigned_int_parameter(infile,"min_steps_between_dt_increase",3);	
	set_double_parameter(infile,"period",2.0);	
	//double density = set_double_parameter(infile,"density",1.176e-6);
	set_double_parameter(infile,"density",1.176e-6);
	//set_double_parameter(infile,"viscosity",15.23); //this is now the same
	set_double_parameter(infile,"zeta_1",-0.057);
	set_double_parameter(infile,"zeta_2",0.2096);
	set_double_parameter(infile,"zeta_3",0.00904);
	set_double_parameter(infile,"E",3.3);
	set_unsigned_int_parameter(infile,"num_generations_1",2);
	set_unsigned_int_parameter(infile,"num_generations_2",5);
	set_unsigned_int_parameter(infile,"num_generations_3",1);
	set_unsigned_int_parameter(infile,"num_generations_4",1);
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
	set_double_parameter(infile,"max_time_step",0.1);
	set_double_parameter(infile,"min_time_step",0.001);
	set_bool_parameter(infile,"neumann_stabilised",true);
	set_bool_parameter(infile,"convective_form",false);
	set_double_parameter(infile,"backflow_stab_param",1.0);
	set_double_parameter(infile,"pressure_mag",1.0);
	set_double_parameter(infile,"time_scaling",1.0);	//this should be updated each timestep as a function of unsteady
	set_double_parameter(infile,"parent_pressure_mag",1.0);
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
	set_unsigned_int_parameter(infile,"adaptive_newton_solve_limit",4);
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
	set_unsigned_int_parameter(infile,"resistance_type_1d",0);	//0- poiseuille 1-pedley
	set_double_parameter(infile,"length_diam_ratio",3.);
	set_unsigned_int_parameter(infile,"bc_type_1d",0);	//0- flow-pressure 1-pressure-pressure
	set_double_parameter(infile,"in_pressure",0);	// pressure at inlet for pressure-pressure BCs

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


  set_string_parameter(infile,"input_1d_node_file_2","");
  set_string_parameter(infile,"input_1d_edge_file_2","");


  set_bool_parameter(infile,"match_1d_mesh_to_3d_mesh",false);
  set_bool_parameter(infile,"assume_symmetric_tree",false);
  set_bool_parameter(infile,"calculate_1d_info_at_coupling_nodes",false);
  set_bool_parameter(infile,"radius_on_edge",false);

  set_double_parameter(infile,"nonlinear_tolerance_1d",1e-4);
  set_bool_parameter(infile,"nonlinear_1d",false);

  set_bool_parameter(infile,"use_centreline_data",false);

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


	set_string_parameter(infile,"petsc_solver_options","");
	set_string_parameter(infile,"petsc_solver_options_ksp_view","");

  // *********************** //

	set_unsigned_int_parameter(infile,"random_1d_generations",0);
	set_string_parameter(infile,"num_generations_string","");

	// 0 is just use the whole matrix, 1 is multiple column solve in libmesh, 2 is multiple column solve in petsc
	set_unsigned_int_parameter(infile,"multiple_column_solve",0);

	set_bool_parameter(infile,"moghadam_coupling",false);

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

	set_double_parameter(infile,"0d_bifurcation_angle",M_PI/4.0);

	set_bool_parameter(infile,"multiple_output_files",true);

	set_bool_parameter(infile,"efficient_monolithic",true);

	set_bool_parameter(infile,"output_system_matrix",false);


  restart_folder << set_string_parameter(infile,"restart_folder",output_folder.str());



	// ******************** now set some derived variables ******************** //
	// 0 - only 3d, 1 - only 1d, 2 - 3d and 1d independent, 3 - implicit coupling
	// 4 - explicit coupling (timestep), 5 - monolithic
	if(sim_type == 0 || sim_type >= 2)
		sim_3d = true;

	if(sim_type >= 1)
		sim_1d = true;

	// set the restart parameter
	if(restart_time_step > 0)
		restart = true;

	t_step = restart_time_step;
	es->parameters.set<unsigned int> ("t_step")   = t_step;




	

	// ********************* RANDOM SEED SETUP ******************************** //
	if(es->parameters.get<int> ("random_seed") == 0)
	{
		es->parameters.set<int> ("random_seed") = std::time(NULL);
		std::cout << "Random seed not given so using std::time(NULL)." << std::endl;
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
						std::cout << "Schur Stokes solver type " << es->parameters.get<unsigned int> ("preconditioner_type_schur_stokes") << " not defined." << std::endl;
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
			if(sim_type != 5)
			{
				petsc_solver_options += " " + es->parameters.get<std::string> ("petsc_1d_solver_options");

			}

			else
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
				if(sim_type != 5)
					petsc_solver_options_ksp_view += " -ns3d_ksp_view";
				else
					petsc_solver_options_ksp_view += " -ns3d1d_ksp_view";
			}
		
			if(sim_1d)
			{
				if(sim_type != 5)
					petsc_solver_options_ksp_view += " -ns1d_ksp_view";
				
			}
		}
			

		es->parameters.set<std::string> ("petsc_solver_options") = petsc_solver_options;
		es->parameters.set<std::string> ("petsc_solver_options_ksp_view") = petsc_solver_options_ksp_view;
		std::cout << "Petsc command line options set to: " << es->parameters.get<std::string> ("petsc_solver_options") << std::endl;;
	}

	// set petsc options from file, petsc command line options will overwrite
	PetscOptionsInsertString(es->parameters.get<std::string> ("petsc_solver_options").c_str());

	if(es->parameters.get<bool> ("direct"))
	{
		if(es->parameters.get<bool> ("fieldsplit"))
		{
			std::cout << "Using direct solver so cannot use fieldsplit. Setting to not use fieldsplit." << std::endl;	
			std::cout << "\tfieldsplit = false" << std::endl;
			es->parameters.set<bool> ("fieldsplit") = false;
		}

		if(es->parameters.get<bool> ("nonzero_initial_guess"))
		{
			std::cout << "Using direct solver so makes no sense to use nonzero_initial_guess. Setting to not use nonzero_initial_guess." << std::endl;
			std::cout << "\tnonzero_initial_guess = false" << std::endl;
			es->parameters.set<bool> ("nonzero_initial_guess") = false;
		}

		if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") || es->parameters.get<unsigned int> ("preconditioner_type_3d"))
		{
			std::cout << "Using direct solver so no preconditioner used. Setting to not use preconditioner." << std::endl;
			std::cout << "\tpreconditioner_type_3d = 0" << std::endl; 
			std::cout << "\tpreconditioner_type_3d1d = 0" << std::endl; 
			es->parameters.set<unsigned int> ("preconditioner_type_3d") = 0;
			es->parameters.set<unsigned int> ("preconditioner_type_3d1d") = 0;
		}
	}



	// ***************** if we aren't doing schur stokes then we can't do schur stokes precompute *********** //
	if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") != 10)
	{
		std::cout << "No point doing schur stokes precompute if we aren't doing schur stokes." << std::endl;
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
		std::cout << "No need to do a multiple column solve if not doing preconditioner_type_3d1d 8,9,10,11." << std::endl;
		std::cout << "Changing to not use multiple_column_solve." << std::endl;
		es->parameters.set<unsigned int> ("multiple_column_solve") = 0;
	}

	if(es->parameters.get<unsigned int> ("multiple_column_solve") && es->parameters.get<unsigned int> ("custom_partitioning") != 1)
	{
		std::cout << "Must use custom partitioning = 1 when doing multiple column solve for now." << std::endl;
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
		std::cout << "Sorry, fieldsplit and schur stokes are not supported with quadratic elements." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}


	// *************** newton 2,3 (newton (no rxn term) and semi-implicit) with optimised stabilised boundary conditions not implemented ********** //
	if((es->parameters.get<unsigned int> ("newton") == 2 || es->parameters.get<unsigned int> ("newton") == 3) 
			&& es->parameters.get<bool> ("optimisation_stabilised"))
	{
		std::cout << "Sorry, newton (no rxn) and semi-implicit are not supported with optimised stabilised boundary conditions." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}


	// *************** can't use semi-implicit with steady **************************************** //

	if(es->parameters.get<unsigned int> ("newton") == 3 && !unsteady)
	{
		std::cout << "Can't use a semi-implicit method with a steady simulation." << std::endl;
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
		std::cout << "Sorry, conservative Navier-Stokes formulation only supported for implicit Picard and Newton (with reaction terms)." << std::endl;
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
		std::cout << "Sorry, streamline diffusion not supported with conservative Navier-Stokes formulation." << std::endl;
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
			std::cout << "Using semi-implicit method with supg/pspg/lsic. Must use have a constant supg constant." << std::endl;
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
			std::cout << "Using semi-implicit method with supg. Must use have a the picard supg method (while using the velocity from the previous timestep)." << std::endl;
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
			std::cout << "Using semi-implicit method with pspg. Must use have a the picard pspg method (while using the velocity from the previous timestep)." << std::endl;
			std::cout << "Changing to pspg_picard - true." << std::endl;

			es->parameters.set<bool> ("pspg_picard") = true;
			es->parameters.set<bool> ("pspg_newton") = false;
				

			std::cout << "New setting: pspg_picard - " <<  es->parameters.get<bool> ("pspg_picard") << std::endl;
		}
	}

	// ********* if semi-implicit and using neumann boundary stabilisation then need to use "linear" form ********* //
	if(es->parameters.get<unsigned int> ("newton") == 3 && es->parameters.get<bool> ("neumann_stabilised") && !es->parameters.get<bool> ("neumann_stabilised_linear")) 
	{
		std::cout << "Using semi-implicit method with neumann stabilisation. Must linear stabilisation." << std::endl;
		std::cout << "Changing to neumann_stabilised_linear - true." << std::endl;

		es->parameters.set<bool> ("neumann_stabilised_linear") = true;
			

		std::cout << "New setting: neumann_stabilised_linear - " <<  es->parameters.get<bool> ("neumann_stabilised_linear") << std::endl;

	}



	// ***************** optimised stabilised bcs only supposed to be used with picard ******** //
	if(false)//es->parameters.get<bool> ("optimisation_stabilised"))
	{
		std::cout << "Currently only picard iterations are supported for optimisation stabilised BCs. Setting to use picard iterations." << std::endl;
		std::cout << "\tnewton = false" << std::endl;
		es->parameters.set<unsigned int> ("newton") = 0;
	}
	else
	{
		//es->parameters.set<bool> ("output_adjoint_variables") = false;
	}
 




	if(es->parameters.get<bool> ("negative_bfbt_schur_complement") && es->parameters.get<unsigned int> ("preconditioner_type_3d") != 9)
	{
		std::cout << "Negative BFBt schur complement currently only implemented for scaled stabilsated bfbt. Changning to positive BFBt schur complement." << std::endl;
		es->parameters.set<bool> ("negative_bfbt_schur_complement") = false;
		
	}


	// ***************** make sure that for sim_type 2 everything is correctly specified ********** //
	if(sim_type == 2 && es->parameters.get<bool> ("known_boundary_conditions"))
	{
		// check that folder has been specified
		if(es->parameters.get<std::string> ("boundary_conditions_folder").empty())
		{
			std::cout << "Folder for getting boundary conditions for partitioned simulation not specified." << std::endl;
			std::cout << "Exiting..." << std::endl;
			std::exit(0);
		}

		// cannot have moghadam coupling and known boundary conditions
		if(es->parameters.get<bool> ("moghadam_coupling"))
		{
			std::cout << "Cannot have known_boundary_conditions and moghadam_coupling both specified." << std::endl;
			std::cout << "Exiting..." << std::endl;
			std::exit(0);
		}
	}

	// ****************** can only have moghadam coupling with sim_type 2 ************************** //
	if(sim_type != 2 && es->parameters.get<bool> ("moghadam_coupling"))
	{
		std::cout << "Ca ony have moghadam_coupling with sim_type 2." << std::endl;
		std::cout << "Changing to not use moghadam coupling." << std::endl;
		es->parameters.set<bool> ("moghadam_coupling") = false;
		
	}



	// *********** adaptive refinement options ********** //
	//if not a 3D simulation then we cannot have adaptive refinement
	if(!sim_3d && adaptive_refinement)
	{
		std::cout << "Can't have adaptive refinement if not simulating 3D. Setting to not use adaptive refinement." << std::endl;	
		std::cout << "\tadaptive_refinement = false" << std::endl;
		es->parameters.set<bool>("adaptive_refinement") = false;
		adaptive_refinement = false;
	}

	// if we have adaptive refinement then we def need an error estimator
	if(es->parameters.get<bool> ("adaptive_refinement") && !es->parameters.get<unsigned int> ("error_estimator") )
	{
		std::cout << "Must have error estimator set when using adaptive refinement. Set to peclet number calculation." << std::endl;
		std::cout << "\terror_estimator = 3" << std::endl;
		es->parameters.set<unsigned int> ("error_estimator") = 3;
	}






	// ********** MAKE SURE STABILISATION OPTIONS ARE CONSISTENT ******************* //

	//if((!es->parameters.get<bool> ("stab") || !es->parameters.get<bool> ("threed")) && es->parameters.get<bool> ("ismail_supg_parameter"))
	if((!es->parameters.get<bool> ("linear_shape_functions")) && es->parameters.get<unsigned int> ("supg_parameter") == 0)
	{
		std::cout << "Cannot use Ismail stabilisation parameters when not using linear elements. Set to use tezdeuyar 2 parameters." << std::endl;
		std::cout << "\tsupg_parameter = 2" << std::endl;
		es->parameters.set<unsigned int> ("supg_parameter") = 2;
	}

	if(!es->parameters.get<bool> ("pspg") && (!es->parameters.get<bool> ("stab") && es->parameters.get<bool> ("linear_shape_functions")))
	{
		std::cout << "When not using pspg stabilisation, must use penalty stabilisation. Setting to use stabilisation with default parameter." << std::endl;
		std::cout << "\tstab = true" << std::endl;
		std::cout << "\talpha = 0.025" << std::endl;
		es->parameters.set<bool> ("stab") = true;
		es->parameters.set<double> ("alpha") = 0.025;
	}

	if((es->parameters.get<bool> ("lsic") || es->parameters.get<bool> ("pspg")) && !es->parameters.get<bool> ("linear_shape_functions"))
	{
		std::cout << "Must use linear shape functions for lsic and/or pspg methinks. Setting to use linear shape functions." << std::endl;
		std::cout << "\tlinear_shape_functions = true" << std::endl;
		es->parameters.set<bool> ("linear_shape_functions") = true;
	}

	if(es->parameters.get<bool> ("stab") && !es->parameters.get<bool> ("linear_shape_functions"))
	{
		std::cout << "Must use linear shape functions for penalty stabilised. Setting to use linear shape functions." << std::endl;
		std::cout << "\tlinear_shape_functions = true" << std::endl;
		es->parameters.set<bool> ("linear_shape_functions") = true;
	}

	if(es->parameters.get<bool>("stokes"))
	{
		if(es->parameters.get<bool>("supg"))
		{
			std::cout << "Can't use SUPG when doing stokes. Setting to not use SUPG."<< std::endl;
			std::cout << "\tsupg = false" << std::endl;
			es->parameters.set<bool>("supg") = false;
		}

		if(es->parameters.get<bool>("pspg"))
		{
			std::cout << "Can't use PSPG when doing stokes. Setting to not use PSPG."<< std::endl;
			std::cout << "\tpspg = false" << std::endl;
			es->parameters.set<bool>("pspg") = false;
		}

		if(es->parameters.get<bool>("lsic"))
		{
			std::cout << "Can't use LSIC when doing stokes. Setting to not use LSIC."<< std::endl;
			std::cout << "\tlsic = false" << std::endl;
			es->parameters.set<bool>("lsic") = false;
		}
	}








	// ****************** SET THE REYNOLDS NUMBER **************************** //
	if(es->parameters.get<bool>("reynolds_number_calculation"))
	{
		es->parameters.set<double> ("length_scale") = 1.0;
		es->parameters.set<double> ("velocity_scale") = 1.0;
		es->parameters.set<double> ("viscosity") = 1.0/reynolds_number;
		es->parameters.set<double> ("density") = 1.0;
	}
	else
	{
		es->parameters.set<double>("reynolds_number") = 
			es->parameters.get<double> ("density") * es->parameters.get<double> ("velocity_scale")
			 * es->parameters.set<double> ("length_scale") / es->parameters.set<double> ("viscosity");

		reynolds_number = es->parameters.get<double>("reynolds_number");
	}

	// we want all the time stuff in the program to be dimenionless
	time_scale_factor = es->parameters.get<double>("length_scale") /
							es->parameters.get<double>("velocity_scale");
	es->parameters.set<double>("time_scale_factor") = time_scale_factor;
	std::cout << "time_scale_factor = " << time_scale_factor << std::endl;
	
	es->parameters.set<double>("max_time_step") /= time_scale_factor;
	es->parameters.set<double>("min_time_step") /= time_scale_factor;
	es->parameters.set<double>("end_time") /= time_scale_factor;
	es->parameters.set<double>("restart_time") /= time_scale_factor;
	es->parameters.set<double>("period") /= time_scale_factor;
	dt /= time_scale_factor;
	es->parameters.set<double>("dt") /= time_scale_factor;

	es->parameters.set<double> ("in_pressure") /= es->parameters.get<double> ("density") 
								* pow(es->parameters.get<double> ("velocity_scale"),2.);





	// **************** SET REYNOLDS NUMBER FOR NUMERICAL CONTINUATION ************** //
	if(unsteady || !es->parameters.get<bool> ("viscosity_controlled"))
	{
		re_vec.push_back(reynolds_number);
	}
	else
	{
		
		double current_re = 1.0;
		current_re = es->parameters.get<double> ("numerical_continuation_starting_reynolds_number");
		std::cout << "Ramping up Reynolds number as:" << std::endl;
		std::cout << "\t";
		while(current_re < reynolds_number)
		{
			std::cout << current_re << ", ";
			re_vec.push_back(current_re);
			current_re = current_re*es->parameters.get<double> ("numerical_continuation_scaling_factor");
		}
		re_vec.push_back(reynolds_number);
		std::cout  << reynolds_number << std::endl;
	
	  es->parameters.set<unsigned int> ("num_continuation_iteration")   = 0;
	}

  es->parameters.set<Real> ("reynolds_number")   = re_vec[0];




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
		std::cout << "If we are doing a pipe geometry then we can only ever have 1 1D tree. Setting to use 1 1D tree." << std::endl;
		std::cout << "\tnum_1d_trees = 1" << std::endl;
		es->parameters.set<unsigned int> ("num_1d_trees") = 1;
	}


	if(es->parameters.get<unsigned int> ("geometry_type") == 5)
	{
		std::cout << "Axisymmetric geometry is currently not working :( !!!... exiting..." << std::endl;
		std::exit(0);
	}



	// ************** ONLY DO GRAVEMEIER ELEMENT LENGTH FOR 3D **************** //
	if(es->parameters.get<bool> ("threed") && es->parameters.get<bool> ("gravemeier_element_length"))
	{
		std::cout << "Can only use gravemeier element length in 3D. Setting to not use gravemeier element length." << std::endl;
		std::cout << "\tgravemeier_element_length = false" << std::endl;
		es->parameters.set<bool> ("gravemeier_element_length") = false;
		
	}


	



	
	// ****************** WE DON'T DO 1D ELEMENTS AT THE MOMENT *********** //
	if(es->parameters.get<bool>("0D") == false)
	{
		std::cout << "Currently only 0D elements are supported. Setting to use 0D elements." << std::endl;
		std::cout << "\t0D = true" << std::endl;
  	es->parameters.set<bool>("0D") = true;
	}

	if(es->parameters.set<unsigned int>("resistance_type_1d") > 0)
	{
		es->parameters.set<bool>("nonlinear_1d") = true;
	}


	// **************** Stokes IC and coupled not implemented ************** //
	if(es->parameters.get<bool>("stokes_ic") && es->parameters.get<unsigned int>("sim_type") > 1)
	{
		std::cout << "Stokes IC with coupled simulation or 0D simulation not currently supported." << std::endl;
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
		std::cout << "Changing the angle of inflow for 3D problems not currently implemented." << std::endl;
		std::cout << "Current setting: angle_of_inflow - " << es->parameters.get<double>("angle_of_inflow") << std::endl;
		std::cout << "Changing to not alter angle of inflow." << std::endl;
		
		es->parameters.set<double>("angle_of_inflow") = 0.;

		std::cout << "New setting: angle_of_inflow - " << es->parameters.get<double>("angle_of_inflow") << std::endl;
		
	}



	// ************** IF WRITING EIGENVALUES MUST BE CALCULATING THEM ***************** //
	if(es->parameters.get<bool> ("write_eigenvalues") && !es->parameters.get<bool> ("compute_eigenvalues"))
	{
		std::cout << "If we are writing eigenvalues to file we must calculate them too." << std::endl;
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
	if(sim_type == 5 && es->parameters.get<unsigned int> ("custom_partitioning") < 1)
	{
		std::cout << "Must use custom partitioning = 1 when using monolithi for some reason..." << std::endl;
		std::cout << "Setting custom_partitioning = 1." << std::endl;
		es->parameters.set<unsigned int> ("custom_partitioning") = 1;

	}


	// ************** PARAMETER TO TELL WHEN STOKES VELOCITY MATRIX HAS BEEN ALLOCATED ******** //
	es->parameters.set<bool> ("stokes_velocity_matrix_set") = false;

	// ************** PARAMETER TO TELL WHEN LINEAR ASSEMBLY HAS BEEN DONE ************ //
	es->parameters.set<bool> ("linear_assembly_done")   = false;

	// ************** EFFICIENT ASSEMBLY NOT SET UP FOR ADAPTIVE TIME-STEPPING *************** //
	// must make sure the time-step hasn't changed so no adaptive time-stepping
	if(es->parameters.get<bool> ("adaptive_time_stepping") && es->parameters.get<bool> ("efficient_assembly"))
	{
		std::cout << "Efficient assembly not set up for adaptive time-stepping yet." << std::endl;
		std::cout << "Changing to use regular assembly." << std::endl;
		es->parameters.set<bool> ("efficient_assembly") = false;
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
		std::cout << "Can't match 1d to 3d mesh if doing a 1d simulations." << std::endl;
		std::cout << "Setting match_1d_mesh_to_3d_mesh=false." << std::endl;
		es->parameters.set<bool> ("match_1d_mesh_to_3d_mesh") = false;
	}


	// ************** SETUP MESH REFINEMENT SETTINGS ******************* //

	if(es->parameters.get<unsigned int> ("output_no_refinement") > es->parameters.get<unsigned int> ("no_refinement") )
	{
		std::cout << "Cannot output solution at a more refined state (need to fix dirichlet condition stuff). Setting output refinement level to same as simulation." << std::endl;
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
																											//but perhaps with a larger mesh this is not the case
	es->parameters.set<bool>("mesh_refined") = false;

	/*
	if(es->parameters.get<unsigned int>("num_refinement_steps") > 0 && 
			!es->parameters.get<unsigned int>("unsteady"))
		es->parameters.set<double>("end_time") = es->parameters.get<unsigned int>("num_refinement_steps") + 2;
	*/



	// ************ SETUP STAGE NUMBERS FOR PETSC PROFILING *************** //
	PetscErrorCode ierr;	// unused
	PetscLogStage stage_num;
	if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 0)
		ierr = PetscLogStageRegister("Exact Solver Dummy",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 2)
		ierr = PetscLogStageRegister("Petsc LSC (BFBt) Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 3)
		ierr = PetscLogStageRegister("Pressure Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 4)
		ierr = PetscLogStageRegister("PCD Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 5)
		ierr = PetscLogStageRegister("PCD2 Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 6)
		ierr = PetscLogStageRegister("LSC James Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 7)
		ierr = PetscLogStageRegister("LSC James Scaled Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 9)
		ierr = PetscLogStageRegister("LSC James Scaled Stabilised Preconditioner Setup",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 10)
		ierr = PetscLogStageRegister("SIMPLE Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 11)
		ierr = PetscLogStageRegister("SIMPLER Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d") == 12)
		ierr = PetscLogStageRegister("SIMPLERC Preconditioner",&stage_num);

	

	if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 6)
		ierr = PetscLogStageRegister("Mono Diagonal Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 7)
		ierr = PetscLogStageRegister("Mono Full Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 8)
		ierr = PetscLogStageRegister("Mono Schur Navier Stokes Velocity Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 9)
		ierr = PetscLogStageRegister("Mono Schur Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 10)
		ierr = PetscLogStageRegister("Mono Schur Stokes Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 11)
		ierr = PetscLogStageRegister("Mono Schur Stokes Velocity Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 12)
		ierr = PetscLogStageRegister("Mono Schur 0D Preconditioner",&stage_num);
	else if(es->parameters.get<unsigned int> ("preconditioner_type_3d1d") == 13)
		ierr = PetscLogStageRegister("Mono Schur NS Approx Preconditioner",&stage_num);



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

// read in parameters from file and give them to the parameters object.
void NavierStokesCoupled::read_particle_parameters()
{
	std::cout << "Reading particle parameters." << std::endl;
	if(particle_deposition == 1)
		restart = true;		// 3D particle deposition is basically a restart

	set_double_parameter(infileparticle,"particle_deposition_start_time",0.0);
	t_step = set_unsigned_int_parameter(infileparticle,"particle_deposition_start_time_step",0);
	es->parameters.set<unsigned int> ("t_step")   = t_step;

	set_bool_parameter(infileparticle,"unsteady_from_steady",false);
	if(es->parameters.get<bool> ("unsteady_from_steady"))
	{
		es->parameters.set<double>("dt") = set_double_parameter(infileparticle,"particle_dt",0.1);
		dt = es->parameters.set<double>("dt");
		es->parameters.set<double>("end_time") = set_double_parameter(infileparticle,"particle_end_time",1.0);
	}


	set_unsigned_int_parameter(infileparticle,"particle_deposition_type",1);
	set_unsigned_int_parameter(infileparticle,"particle_deposition_location",0);
	set_double_parameter(infileparticle,"particle_deposition_rate",10);
	set_unsigned_int_parameter(infileparticle,"particle_deposition_surface",0);
	set_double_parameter(infileparticle,"distance_deposited_inside",1.0);
	set_unsigned_int_parameter(infileparticle,"deposition_pattern",0);
	set_double_parameter(infileparticle,"brownian_motion_magnitude",0.);
	set_unsigned_int_parameter(infileparticle,"deposit_within_radius",0.);

	//some parameters
	set_double_parameter(infileparticle,"particle_diameter",10.e-6);
	set_double_parameter(infileparticle,"particle_density",1200.);
	set_double_parameter(infileparticle,"particle_air_density",1.2);
	set_double_parameter(infileparticle,"particle_air_viscosity",2.04e-5);
	set_double_parameter(infileparticle,"mean_free_path",68.e-9);
	set_double_parameter(infileparticle,"gravity",9.81);

	set_bool_parameter(infileparticle,"particle_sedimentation",false);
	set_bool_parameter(infileparticle,"particle_impaction",false);
	set_bool_parameter(infileparticle,"particle_drag",false);

	set_unsigned_int_parameter(infileparticle,"gravity_type",0);

	set_double_parameter(infileparticle,"particle_velocity_units",1.);
	set_double_parameter(infileparticle,"initial_particle_velocity_scaling",0.99);

	set_double_parameter(infileparticle,"gravity_x",0.);
	set_double_parameter(infileparticle,"gravity_y",0.);
	set_double_parameter(infileparticle,"gravity_z",1.);

	set_bool_parameter(infileparticle,"deposition_verbose",false);


}



// output parameters to screen and to file.
void NavierStokesCoupled::output_particle_parameters()
{

	// **************** COPY THE INPUT FILE TO TH OUTPUT FOLDER **************** //
	if(!es->parameters.get<bool>("compare_results"))
	{
		std::ifstream  input_file_particle_src(input_file_particle.c_str(), std::ios::binary);
		std::ofstream  input_file_particle_dst(std::string(es->parameters.get<std::string>("output_folder") + "particle.in").c_str(),   std::ios::binary);
		input_file_particle_dst << input_file_particle_src.rdbuf();
	}

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
void NavierStokesCoupled::read_timesteps()
{
	if(es->parameters.get<bool>("unsteady_from_steady"))
	{
		t_step = 0;
		es->parameters.set<unsigned int> ("t_step")   = t_step;

	}
	else
	{

		std::ifstream timestep_file((es->parameters.get<std::string>("output_folder") + "out.dat").c_str());
		std::cout << "input file = " << (es->parameters.get<std::string>("output_folder") + "out.dat").c_str() << std::endl;
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


		std::cout << "CHECK TIMESTEP SIZES VECTOR" << std::endl;
		for(unsigned int i=0; i< timestep_sizes.size(); i++)
			std::cout << "i = " << i << " : " << timestep_sizes[i] << std::endl;
	}



}


// read the input boundary conditions for an uncoupled simulation
int NavierStokesCoupled::read_input_boundary_conditions()
{
	std::cout << "Reading input boundary conditions file." << std::endl;

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
	// then divide by 4 for the 3d and 0d fluxes and pressures.
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

	if(total_num_timesteps != input_num_timesteps)
	{
		std::cout << "Error, number of input timesteps not equal to number of expected timesteps." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}

	// 4 - set input boundary conditions for first timestep

	std::cout << "t_step = " << t_step << std::endl;
	input_pressure_values_3d = all_input_pressure_values_3d[t_step];
	input_flux_values_1d = all_input_flux_values_0d[t_step];


}


void NavierStokesCoupled::print_flux_and_pressure()
{

	double flux_scaling = 1.0;
	double pressure_scaling = 1.0;
	if(!es->parameters.get<bool>("reynolds_number_calculation"))
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
		std::cout << "total daughter flux 3D = " << total_daughter_flux_3d * flux_scaling << " m^3/s)" << std::endl;
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
	es->parameters.set<double> ("previous_dt") = previous_dt;
	es->parameters.set<unsigned int> ("t_step") = t_step;

	update_time_scaling();

}

void NavierStokesCoupled::update_time_scaling()
{
	double time_flow = 0.0;

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
			// update the timestep
			dt = timestep_sizes[t_step];		//get the next time step, time should be updated itself...
			update_times();

			// get the velocity for next timestep
			std::ostringstream file_name;
			std::ostringstream file_name_es;

			file_name << output_folder.str() << "out_3D";
			// We write the file in the ExodusII format.
			file_name_es << file_name.str() << "_es_" << std::setw(4) 
				<< std::setfill('0') << t_step << ".xda";
			es->read(file_name_es.str(), libMeshEnums::READ);

			if(sim_3d && sim_type != 5)
			{
				system_3d->update();	
			}

			if(sim_1d && sim_type != 5)
			{
				system_1d->update();	
			}

			if(sim_type == 5)
			{
				system_coupled->update();	
			}

		}


	}
	else
	{
		// when we reduce dt or refine the mesh always make sure that solution has been updated to old solution so that
		// it is correctly projected. local stuff can handle itself

		//if didn't converge or something then go back to previous time_step
		if((reduce_dt && unsteady) || refine_mesh)
		{
			time -= dt; // go back
			t_step -= 1;
		
			if(reduce_dt)
				dt = 0.5*dt;

			if(dt < es->parameters.get<double>("min_time_step"))
				dt = es->parameters.get<double>("min_time_step");
		
			update_times();

			if(sim_3d && sim_type != 5)
			{
				*system_3d->current_local_solution = *system_3d->old_local_solution;
				*old_global_solution = *system_3d->solution;	// remember we changed solution to old global before refinement
				//system_3d->update();	//not sure why this reverts, because it somehow uses the solution and then distributes
				//note that it means here the solution is not correct
				// then after we're gonna need to make old_global_solution have the same size as solution
			}

			if(sim_1d && sim_type != 5)
			{
				*system_1d->current_local_solution = *system_1d->old_local_solution;
				system_1d->update();	
			}

			if(sim_type == 5)
			{
				*system_coupled->current_local_solution = *system_coupled->old_local_solution;
				*old_global_solution = *system_coupled->solution;

			}

			steps_since_last_dt_change = 0;
			reduce_dt = false;
		}

		//if did converge but a joke then double time step
		// (first_residual < 0.5 && system_3d->time > 30) removed this condition, made it just based on the number of iterations
		else if(increase_dt && unsteady)
		{

			previous_flux_values_1d = flux_values_1d;
			previous_previous_flux_values_3d = previous_flux_values_3d;
			previous_flux_values_3d = flux_values_3d;
			previous_pressure_values_1d = pressure_values_1d;
			previous_pressure_values_3d = pressure_values_3d;

			previous_dt = dt;
			update_times();
		

			dt = 2.0*dt;

			if(dt > es->parameters.get<double>("max_time_step"))
				dt = es->parameters.get<double>("max_time_step");

			update_times();


			if(sim_3d && sim_type != 5)
			{
				*system_3d->old_local_solution = *system_3d->current_local_solution;
				*old_global_solution = *system_3d->solution;
				//system_3d->update();
			}

			if(sim_1d && sim_type != 5)
			{
				*system_1d->old_local_solution = *system_1d->current_local_solution;
				system_1d->update();	
			}

			if(sim_type == 5)
			{
				*old_global_solution = *system_coupled->solution;
				*system_coupled->old_local_solution = *system_coupled->current_local_solution;
				system_coupled->update();	
			}

			steps_since_last_dt_change = 0;
			increase_dt = false;

			std::cout << "increasing timestep" << std::endl;
		}
		//if we've been worked hard enough then just 
		else
		{	

			previous_flux_values_1d = flux_values_1d;
			previous_previous_flux_values_3d = previous_flux_values_3d;
			previous_flux_values_3d = flux_values_3d;
			previous_pressure_values_1d = pressure_values_1d;
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
		}
	}
}






void NavierStokesCoupled::output_sim_data(bool header)
{

	double flux_scaling = 1.0;
	double pressure_scaling = 1.0;

	if(!es->parameters.get<bool>("reynolds_number_calculation"))
	{
		flux_scaling = es->parameters.get<double>("velocity_scale") * 
										pow(es->parameters.get<double>("length_scale"),2.0);
		pressure_scaling = es->parameters.get<double>("density") *
													pow(es->parameters.get<double>("velocity_scale"),2.0);
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
		}

		if(sim_3d)
			output_file << "\tnum_newton_steps";

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

		}

		if(sim_3d)
			output_file << "\t" <<	es->parameters.get<unsigned int> ("num_newton_steps");

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
		linear_iterations_output_file << std::endl;
	}
	else
	{
		linear_iterations_output_file << t_step << "\t" << nonlinear_iteration << 
																							"\t" << local_linear_iterations << 
																							"\t" << total_max_iterations << 
																							"\t" << total_nonlinear_iterations <<
																							"\t" << total_linear_iterations << 
																							"\t" << es->parameters.get<double> ("last_nonlinear_iterate") << std::endl;
	}
}


// write the particle data for the timestep (old)
void NavierStokesCoupled::output_particle_data_old(bool header)
{

	std::cout << "yeah" << std::endl;
	if(header)
	{
		// write geometry variables later when reading meta data file
		// write boundary conditions later with Picard class
		particle_output_file << "# Particle results" << std::endl;
		particle_output_file << "# timestep\tcurrent_time";

		for(unsigned int i=0; i < particles_3D.size(); i++)
		{
			particle_output_file <<	"\tp" << i << "_pos_x";
			particle_output_file <<	"\tp" << i << "_pos_y";
			if(threed)
				particle_output_file <<	"\tp" << i << "_pos_z";
			particle_output_file <<	"\tp" << i << "_exit_surface";
			particle_output_file <<	"\tp" << i << "_exit_time";
			particle_output_file <<	"\tp" << i << "_entrance_time";
			particle_output_file <<	"\tp" << i << "_pid";

		}
	
		particle_output_file << std::endl;
	}
	else
	{
		particle_output_file << t_step << "\t" << time;

		for(unsigned int i=0; i < particles_3D.size(); i++)
		{
			Point particle_position = particles_3D[i].get_position();
			particle_output_file <<	"\t" << particle_position(0);
			particle_output_file <<	"\t" << particle_position(1);
			if(threed)
				particle_output_file <<	"\t" << particle_position(2);
			particle_output_file <<	"\t" << particles_3D[i].get_exit_surface();
			particle_output_file <<	"\t" << particles_3D[i].get_exit_time();
			particle_output_file <<	"\t" << particles_3D[i].get_entrance_time();
			particle_output_file <<	"\t" << particles_3D[i].get_particle_id();
		}

		particle_output_file << std::endl;
	}
}


// write the particle data for the timestep
void NavierStokesCoupled::output_particle_data()
{


	//particle deposition file
	std::ostringstream particle_output_data_file_name;

	particle_output_data_file_name << output_folder.str() << "particles";
	particle_output_data_file_name << std::setw(4) << std::setfill('0') << t_step;
	particle_output_data_file_name << ".dat";
	//output_data_file_name << "results/out_viscosity"	<< es->parameters.set<Real> ("viscosity") << ".dat";
	particle_output_file.open(particle_output_data_file_name.str().c_str());




	// write geometry variables later when reading meta data file
	// write boundary conditions later with Picard class
	particle_output_file << "# Particle results" << std::endl;
	particle_output_file << "# pid\ttimestep\tcurrent_time\tx\ty\t";
	if(threed)
		particle_output_file << "z\t";
	particle_output_file << "exit_surface\texit_time\tentrance_time\n";


	for(unsigned int i=0; i < particles_3D.size(); i++)
	{
		unsigned int pid = particles_3D[i].get_particle_id();
		Point particle_position = particles_3D[i].get_position();
		double x = particle_position(0);
		double y = particle_position(1);
		double z = 0.;
		if(threed)
			z = particle_position(2);
		int exit_surface = particles_3D[i].get_exit_surface();
		double exit_time = particles_3D[i].get_exit_time();
		double entrance_time = particles_3D[i].get_entrance_time();

		
		particle_output_file << pid << "\t";
		particle_output_file << t_step << "\t";
		particle_output_file << time << "\t";
		particle_output_file << x << "\t";
		particle_output_file << y << "\t";
		if(threed)
			particle_output_file << z << "\t";
		particle_output_file << exit_surface << "\t";
		particle_output_file << exit_time << "\t";
		particle_output_file << entrance_time << "\n";
	}

	particle_output_file.close();
}


// print some particle data for the timestep
void NavierStokesCoupled::print_particle_data()
{

	unsigned int num_valid_particles = 0;
	unsigned int num_particles_on_wall = 0;
	unsigned int num_particles_exited = 0;
	unsigned int num_lost_particles = 0;
	unsigned int num_particles_in_motion = 0;
	unsigned int num_particles_close_to_wall = 0;

	for(unsigned int i=0; i < particles_3D.size(); i++)
	{
		int exit_surface = particles_3D[i].get_exit_surface();
		double exit_time = particles_3D[i].get_exit_time();
		bool close_to_wall = particles_3D[i].was_deposited_close_to_wall();

		if(exit_time > 0.)
		{
			if(exit_surface == -1)
			{
				num_particles_on_wall++;
				num_valid_particles++;
			}
			else if(exit_surface == -2)
			{
				// none
				num_lost_particles++;
			}
			else
			{
				num_particles_exited++;
				num_valid_particles++;
			}

			if(close_to_wall)
			{
				num_particles_close_to_wall++;
			}
		}
		else
		{
			num_particles_in_motion++;
		}
		
	}

	std::cout << "num_valid_particles = " << num_valid_particles << std::endl;
	std::cout << "num_particles_on_wall = " << num_particles_on_wall << std::endl;
	std::cout << "num_particles_exited = " << num_particles_exited << std::endl;
	std::cout << "num_lost_particles = " << num_lost_particles << std::endl;
	std::cout << "num_particles_in_motion = " << num_particles_in_motion << std::endl;
	std::cout << "num_particles_close_to_wall = " << num_particles_close_to_wall << std::endl;
	std::cout << "deposition_fraction = " << (double)num_particles_on_wall/(double)num_valid_particles << std::endl;
}



void NavierStokesCoupled::output_logging()
{

	std::ostringstream logging_data_file_name;

	logging_data_file_name << output_folder.str() << "perf_log.dat";
	std::ofstream logging_output_file(logging_data_file_name.str().c_str());

	logging_output_file << perf_log.get_log();
	logging_output_file.close();
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
	if(!es->parameters.get<bool>("reynolds_number_calculation"))
		return time*time_scale_factor;
	else
		return time;
}

void NavierStokesCoupled::init_dof_variable_vectors()
{

	if(sim_type == 0 || sim_type == 2 || sim_type == 3 || sim_type == 4)
		init_dof_variable_vector(system_3d, dof_variable_type_3d);

	if(sim_type == 1 || sim_type == 2 || sim_type == 3 || sim_type == 4)
		init_dof_variable_vector(system_1d, dof_variable_type_1d);

	if(sim_type == 5)
		init_dof_variable_vector(system_coupled, dof_variable_type_coupled);
}

void NavierStokesCoupled::init_dof_variable_vector(TransientLinearImplicitSystem * system, std::vector<int>& _dof_variable_type)
{

	std::cout << "initing dof variable vector" << std::endl;	

	// DofMap things
  std::vector<dof_id_type> dof_indices_var;

  const DofMap & dof_map = system->get_dof_map();
  const int nv_sys = system->n_vars();
	_dof_variable_type.resize(0);
	_dof_variable_type.resize(system->solution->size(),0);

	for (int var=0; var<nv_sys; var++)
	{
		MeshBase::element_iterator       it       = mesh.active_elements_begin();
		const MeshBase::element_iterator end_elem = mesh.active_elements_end();
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

void NavierStokesCoupled::init_particles()
{
	// let us set gravity here
	if(es->parameters.get<unsigned int>("gravity_type") == 2)
	{
		SurfaceBoundary* deposition_boundary = surface_boundaries[es->parameters.get<unsigned int>("particle_deposition_surface")];
		Point normal = deposition_boundary->get_normal();

		if(threed)
		{
			es->parameters.set<double>("gravity_x") = -normal(0);
			es->parameters.set<double>("gravity_y") = -normal(1);
			es->parameters.set<double>("gravity_z") = -normal(2);
		}
		else
		{
			es->parameters.set<double>("gravity_x") = -normal(0);
			es->parameters.set<double>("gravity_y") = -normal(1);
		}

	}

	//check the particle is indeed in the mesh
	const PointLocatorBase & pl = mesh.point_locator();

	// figure out number of particles to deposit
	// particle_deposition_type 
	// 0 - all at once
	// 1 - at a rate
	unsigned int num_new_particles = 0;
	if(es->parameters.get<unsigned int>("particle_deposition_type") == 0) 
		num_new_particles = (unsigned int) es->parameters.get<double>("particle_deposition_rate");
	else if(es->parameters.get<unsigned int>("particle_deposition_type") == 1)
		num_new_particles = (unsigned int) es->parameters.get<double>("particle_deposition_rate") * dt;
	


	if(es->parameters.get<unsigned int>("particle_deposition_location") == 0)
	{
		if(!threed)
		{
			Point new_point(0.8,0.7,0.);

			const Elem * element_found = pl(new_point);
			if(element_found == NULL)
			{
				std::cout << "particle not in mesh EXITING" << std::endl;
				std::exit(0);
			}

			Particle new_particle(*es,new_point,element_found,particles_3D.size());
			particles_3D.push_back(new_particle);

			Point new_point_2(0.9,0.5,0.);
	
			const Elem * element_found_2 = pl(new_point_2);
			if(element_found_2 == NULL)
			{
				std::cout << "particle not in mesh EXITING" << std::endl;
				std::exit(0);
			}

			Particle new_particle_2(*es,new_point_2,element_found_2,particles_3D.size());
			particles_3D.push_back(new_particle_2);
		}
		else
		{

			Point new_point(0.2,0.1,0.5);

			const Elem * element_found = pl(new_point);
			if(element_found == NULL)
			{
				std::cout << "particle not in mesh EXITING" << std::endl;
				std::exit(0);
			}
			Particle new_particle(*es,new_point,element_found,particles_3D.size());
			particles_3D.push_back(new_particle);

			Point new_point_2(-0.3,0.2,0.5);

			const Elem * element_found_2 = pl(new_point_2);
			if(element_found_2 == NULL)
			{
				std::cout << "particle not in mesh EXITING" << std::endl;
				std::exit(0);
			}

			Particle new_particle_2(*es,new_point_2,element_found_2,particles_3D.size());
			particles_3D.push_back(new_particle_2);
		}
	}
	else if(es->parameters.get<unsigned int>("particle_deposition_location") == 1)
	{

		SurfaceBoundary* deposition_boundary = surface_boundaries[es->parameters.get<unsigned int>("particle_deposition_surface")];
		double max_radius = deposition_boundary->get_max_radius();
		Point normal = deposition_boundary->get_normal();
		Point centroid = deposition_boundary->get_centroid();
		double distance_deposited_inside = es->parameters.get<double>("distance_deposited_inside");

		std::cout << "centroid = " << centroid << std::endl;
		std::cout << "normal = " << normal << std::endl;
		std::cout << "max_radius = " << max_radius << std::endl;

		max_radius *= 1.0;


		std::cout << "depositing " << num_new_particles << " particles." << std::endl;
		
		if(threed)
		{

			//from normal get two orthogonal vectors
			Point vector_1(0.,0.,0.);
			Point vector_2(0.,0.,0.);
			if(fabs(normal(2)) > 1e-10)
			{
				vector_1(0) = 0.;
				vector_1(1) = 1.;
				vector_1(2) = -normal(1)/normal(2);
			}
			else if(fabs(normal(1)) > 1e-10)
			{
				vector_1(0) = 0.;
				vector_1(1) = -normal(2)/normal(1);
				vector_1(2) = 1.;
			}
			else if(fabs(normal(0)) > 1e-10)
			{
				vector_1(0) = -normal(1)/normal(0);
				vector_1(1) = 1.;
				vector_1(2) = 0.;
			}
			else
			{
				std::cout << "wth, fuck you" << std::endl;
				std::exit(0);
			}

			vector_1 = vector_1.unit();
			vector_2 = vector_1.cross(normal);
		
			std::cout << "vector_1 = " << vector_1 << std::endl;
			std::cout << "vector_2 = " << vector_2 << std::endl;

			unsigned int num_failed_particles = 0;

			for(unsigned int i=0; i<num_new_particles; i++)
			{
				bool particle_deposited = false;
				while(!particle_deposited)
				{
					Point new_point;				

					// square deposition pattern
					if(es->parameters.get<unsigned int>("deposition_pattern") == 0)
					{
						// random number between -1 and 1 multiplied by the radius to get -r to r
						double x = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * max_radius;
						double y = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * max_radius;

						//make parabolic distribution? i.e. more at the centre?
						x *= x;
						y *= y;

						// make negative maybe
						if(static_cast<double>(rand())/static_cast<double>(RAND_MAX) > 0.5)
							x *= -1;

						if(static_cast<double>(rand())/static_cast<double>(RAND_MAX) > 0.5)
							y *= -1;					

						new_point = centroid + x*vector_1 + y*vector_2;


					}
					// uniform circular pattern
					else if(es->parameters.get<unsigned int>("deposition_pattern") == 1)
					{
						// according to wolfram, needs to be distributed as sqrt(r)cos(theta) to be uniform
						double theta = (static_cast<double>(rand())/static_cast<double>(RAND_MAX)) * 2 * M_PI;
						double radius_scale = (static_cast<double>(rand())/static_cast<double>(RAND_MAX));

						// figure out what the radius is
						double radius_at_angle = deposition_boundary->get_radius_from_angle(theta);
						double x = sqrt(radius_scale) * radius_at_angle * cos(theta);
						double y = sqrt(radius_scale) * radius_at_angle * sin(theta);

						new_point = centroid + x*vector_1 + y*vector_2;
					}

					std::cout << "new_point = " << new_point << std::endl;
					//distance_deposited_inside = 1.e-3;

					if(deposition_boundary->is_on_surface(new_point))
					{

						const Elem * element_found = pl(new_point - distance_deposited_inside*normal);
						if(element_found != NULL)
						{
							new_point = new_point - distance_deposited_inside*normal;
							//normal is outward facing and want to deposit the particle a little bit inside
							std::cout << "hmmmm " << new_point << std::endl;
							Particle new_particle(*es,new_point,element_found,particles_3D.size());
							particles_3D.push_back(new_particle);
							std::cout << "particle deposited at " << new_point << std::endl;
							particle_deposited = true;
						}
						else
						{
							std::cout << "particle should have been found but wasn't" << std::endl;
							num_failed_particles++;
						}
					}
					else
					{
						std::cout << "particle is not on surface." << std::endl;
						num_failed_particles++;
					}
				}
			}
			std::cout << "num_failed_particles = " << num_failed_particles << std::endl;
				
		}
		else
		{

			//from normal get orthogonal vector
			Point vector_1(normal(1),-normal(0),0.);

			for(unsigned int i=0; i<num_new_particles; i++)
			{
				bool particle_deposited = false;
				while(!particle_deposited)
				{
					// random number between -1 and 1 multiplied by the radius to get -r to r
					double x = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * max_radius;
					Point new_point = centroid + x*vector_1;
					if(deposition_boundary->is_on_surface(new_point))
					{
	
						const Elem * element_found = pl(new_point - distance_deposited_inside*normal);
						if(element_found != NULL)
						{
							new_point = new_point - distance_deposited_inside*normal;
							//normal is outward facing and want to deposit the particle a little bit inside
							Particle new_particle(*es,new_point,element_found,particles_3D.size());
							particles_3D.push_back(new_particle);
							std::cout << "particle deposited at " << new_point << std::endl;
							particle_deposited = true;
						}
						else
							std::cout << "particle should have been found but wasn't..." << std::endl;
					}
				}
			}
				
		}

		
		std::cout << num_new_particles << " particles deposited." << std::endl;

	}
	

}

void NavierStokesCoupled::move_particles()
{

	// it's as easy as this
	for(unsigned int i=0; i<particles_3D.size(); i++)
	{
		particles_3D[i].try_and_move();
	}

}


// write particle data file for this time step
void NavierStokesCoupled::write_particles()
{

	// output the parameters used to file and the header of the out.dat file
	std::ostringstream particle_data_file_name;
	std::ofstream particle_vtk_output_file;

	particle_data_file_name << output_folder.str() << "particles";

	particle_data_file_name << std::setw(4) << std::setfill('0') << t_step;
	particle_data_file_name << ".vtk" ;
	//output_data_file_name << "results/out_viscosity"	<< es->parameters.set<Real> ("viscosity") << ".dat";
	particle_vtk_output_file.open(particle_data_file_name.str().c_str());

	particle_vtk_output_file << "# vtk DataFile Version 2.0" << std::endl;
	particle_vtk_output_file << "Particle file for paraview" << std::endl;
	particle_vtk_output_file << "ASCII" << std::endl;
	particle_vtk_output_file << "DATASET UNSTRUCTURED_GRID" << std::endl << std::endl;
	particle_vtk_output_file << "POINTS " << particles_3D.size() << " float" << std::endl;

	for(unsigned int i=0; i<particles_3D.size(); i++)
	{
		Point particle_position = particles_3D[i].get_position();
		particle_vtk_output_file << particle_position(0) << " ";
		particle_vtk_output_file << particle_position(1) << " ";
		if(threed)
			particle_vtk_output_file << particle_position(2) << std::endl;
		else
			particle_vtk_output_file << 0 << std::endl;
	}

	particle_vtk_output_file << std::endl;
	particle_vtk_output_file << "CELL_TYPES " << particles_3D.size() << std::endl;

	for(unsigned int i=0; i<particles_3D.size(); i++)
	{
		particle_vtk_output_file << 1 << std::endl;
	}

	particle_vtk_output_file << std::endl;

	particle_vtk_output_file << "POINT_DATA " << particles_3D.size() << std::endl;
	particle_vtk_output_file << "SCALARS status float 1" << std::endl;
	particle_vtk_output_file << "LOOKUP_TABLE default" << std::endl;

	for(unsigned int i=0; i<particles_3D.size(); i++)
	{
		particle_vtk_output_file << particles_3D[i].get_status() << " ";
	}

	particle_vtk_output_file << std::endl;
	particle_vtk_output_file << std::endl;

	particle_vtk_output_file << "SCALARS exit_surface float 1" << std::endl;
	particle_vtk_output_file << "LOOKUP_TABLE default" << std::endl;

	for(unsigned int i=0; i<particles_3D.size(); i++)
	{
		particle_vtk_output_file << particles_3D[i].get_exit_surface() << " ";
	}

	particle_vtk_output_file << std::endl;
	particle_vtk_output_file << std::endl;

	particle_vtk_output_file << "SCALARS exit_time float 1" << std::endl;
	particle_vtk_output_file << "LOOKUP_TABLE default" << std::endl;

	for(unsigned int i=0; i<particles_3D.size(); i++)
	{
		particle_vtk_output_file << particles_3D[i].get_exit_time() << " ";
	}

	particle_vtk_output_file << std::endl;
	particle_vtk_output_file << std::endl;

	particle_vtk_output_file << "SCALARS entrance_time float 1" << std::endl;
	particle_vtk_output_file << "LOOKUP_TABLE default" << std::endl;

	for(unsigned int i=0; i<particles_3D.size(); i++)
	{
		particle_vtk_output_file << particles_3D[i].get_entrance_time() << " ";
	}

	particle_vtk_output_file << std::endl;
	particle_vtk_output_file << std::endl;



	particle_vtk_output_file << "SCALARS particle_id float 1" << std::endl;
	particle_vtk_output_file << "LOOKUP_TABLE default" << std::endl;

	for(unsigned int i=0; i<particles_3D.size(); i++)
	{
		particle_vtk_output_file << particles_3D[i].get_particle_id() << " ";
	}

	particle_vtk_output_file << std::endl;
	particle_vtk_output_file << std::endl;


	particle_vtk_output_file << "VECTORS velocity float" << std::endl;
	for(unsigned int i=0; i<particles_3D.size(); i++)
	{
		NumberVectorValue particle_velocity = particles_3D[i].get_current_velocity();
		particle_vtk_output_file << particle_velocity(0) << " ";
		particle_vtk_output_file << particle_velocity(1) << " ";
		if(threed)
			particle_vtk_output_file << particle_velocity(2) << " ";
		else
			particle_vtk_output_file << 0 << "  ";
	}


	particle_vtk_output_file.close();



		

}


// here we setup the variable scalings that we want to use. 
// the simulation is done in nondimensional variables and output in SI units. (not the particle deposition)
void NavierStokesCoupled::setup_variable_scalings_3D()
{

	std::cout << "Setting up 3D variable scalings." << std::endl;

	double velocity_scale = es->parameters.get<double>("velocity_scale");
	double pressure_scale = es->parameters.get<double>("density") *	pow(es->parameters.get<double>("velocity_scale"),2.0);

	if(es->parameters.get<bool>("output_nondim"))
	{
		velocity_scale = 1.0;
		pressure_scale = 1.0;
		std::cout << "using nondim output" << std::endl;
	}

	if(es->parameters.get<bool>("reynolds_number_calculation"))
	{
		velocity_scale = 1.0;
		pressure_scale = 1.0;

	}


	// these indices are based on the order in which they are given to exodus
	if(!threed)
	{
		var_scalings_3D.resize(3);
		var_scalings_3D[0] = velocity_scale;	// u
		var_scalings_3D[1] = velocity_scale;	// v
		var_scalings_3D[2] = pressure_scale;	// p
	}
	else
	{
		var_scalings_3D.resize(4);
		var_scalings_3D[0] = velocity_scale;	// u
		var_scalings_3D[1] = velocity_scale;	// v
		var_scalings_3D[2] = velocity_scale;	// w
		var_scalings_3D[3] = pressure_scale;	// p
	}


	if(es->parameters.get<bool>("optimisation_stabilised"))
	{
		if(!threed)
		{
			var_scalings_3D.resize(6);
			var_scalings_3D[3] = velocity_scale;	// u_adj
			var_scalings_3D[4] = velocity_scale;	// v_adj
			var_scalings_3D[5] = pressure_scale;	// p_adj
		}
		else
		{		
			var_scalings_3D.resize(8);
			var_scalings_3D[4] = velocity_scale;	// u_adj
			var_scalings_3D[5] = velocity_scale;	// v_adj
			var_scalings_3D[6] = velocity_scale;	// w_adj
			var_scalings_3D[7] = pressure_scale;	// p_adj
		}
	}


}


// here we setup the variable scalings that we want to use. 
// the simulation is done in nondimensional variables and output in SI units. (not the particle deposition)
void NavierStokesCoupled::setup_variable_scalings_1D()
{

	std::cout << "Setting up 1D variable scalings." << std::endl;

	double flow_scale = es->parameters.get<double>("velocity_scale") * 
										pow(es->parameters.get<double>("length_scale"),2.0);
	double mean_pressure_scale = es->parameters.get<double>("density") *
													pow(es->parameters.get<double>("velocity_scale"),2.0);

	if(es->parameters.get<bool>("output_nondim"))
	{
		flow_scale = 1.0;
		mean_pressure_scale = 1.0;
		std::cout << "using nondim output" << std::endl;
	}

	if(es->parameters.get<bool>("reynolds_number_calculation"))
	{
		flow_scale = 1.0;
		mean_pressure_scale = 1.0;

	}

	var_scalings_1D.resize(5);
	var_scalings_1D[0] = mean_pressure_scale;	// P
	var_scalings_1D[1] = flow_scale;	// Q
	var_scalings_1D[2] = 1.0;	// radius
	var_scalings_1D[3] = 1.0;	// poiseuille
	var_scalings_1D[4] = 1.0;	// proc_id


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



