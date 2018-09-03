#include "coupled_navier_stokes.h"

// ********************************************************** //
//
// This file includes all the functions not specific to 
// particle deposition.
//
// read_particle_parameters
// output_particle_parameters
// output_particle_data_old
// output_particle_data
// print_particle_data
// no_particles_in_motion
// init_particles
// move_particles
// write_particles
// 
//
// *********************************************************** //



// read in parameters from file and give them to the parameters object.
void NavierStokesCoupled::read_particle_parameters()
{
	std::cout << "Reading particle parameters." << std::endl;
	if(particle_deposition == 1 || particle_deposition == 3 ||  particle_deposition == 4 ||  particle_deposition == 5 ||  particle_deposition == 6)
		restart = true;		// 3D particle deposition is basically a restart

	set_double_parameter(infileparticle,"particle_deposition_start_time",0.0);
	t_step = set_unsigned_int_parameter(infileparticle,"particle_deposition_start_time_step",0);
	es->parameters.set<unsigned int> ("t_step")   = t_step;

	set_bool_parameter(infileparticle,"unsteady_from_steady",false);

	// always set the time stuff from the particle file
	es->parameters.set<double>("dt") = set_double_parameter(infileparticle,"particle_dt",0.1);
	set_double_parameter(infileparticle,"particle_0d_dt",es->parameters.get<double>("particle_dt"));

	std::cout << "FUCK ME" << std::endl;

	// we must have that particle_dt is less than particle particle_0d_dt
	if(es->parameters.get<double>("particle_dt") > es->parameters.get<double>("particle_0d_dt") + 1e-10)
	{
		std::cout << "Error, 3d particle time step (particle_dt) is larger than 0d particle time step (particle_0d_dt)." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}
	es->parameters.set<double>("time") = 0.;	// make sure this is set
	dt = es->parameters.set<double>("dt");
	es->parameters.set<double>("end_time") = set_double_parameter(infileparticle,"particle_end_time",1.0);
	set_double_parameter(infileparticle,"read_end_time",0.);
	set_bool_parameter(infileparticle,"read_time_steps",false);	// in navier.in too, but doesn't matter, this overrides it
	set_double_parameter(infileparticle,"read_time_step_size",0.1);	// in navier.in too, but doesn't matter, this overrides it
	set_double_parameter(infileparticle,"particle_write_interval",0.1);	// in navier.in too, but doesn't matter, this overrides it
	set_double_parameter(infileparticle,"particle_backup_write_interval",es->parameters.get<double>("particle_write_interval"));	// in navier.in too, but doesn't matter, this overrides it

	// write interval can't be more than the particle_0d_dt
	// backup write interval can't be less than write interval
	if(es->parameters.get<double>("particle_0d_dt") > es->parameters.get<double>("particle_write_interval") + 1e-10)
	{
		std::cout << "Error, particle_0d_dt is larger than particle_write_interval." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}

	// backup write interval can't be less than write interval
	if(es->parameters.get<double>("particle_write_interval") > es->parameters.get<double>("particle_backup_write_interval") + 1e-10)
	{
		std::cout << "Error, particle_write_interval is larger than particle_backup_write_interval." << std::endl;
		std::cout << "Exiting..." << std::endl;
		std::exit(0);
	}


	// if we come here, we are running a particle simulation 
	// so set write_interval and backup_write_interval to the 
	// particle counterparts (cause 1d deposition uses the same function)
	es->parameters.set<double>("write_interval") = es->parameters.get<double>("particle_write_interval");
	es->parameters.set<double>("backup_write_interval") = es->parameters.get<double>("particle_backup_write_interval");
	


	// check that end_time is at least less than the end time
	if(!es->parameters.get<bool>("unsteady_from_steady"))
	{
		if(es->parameters.get<double>("end_time") > es->parameters.get<double>("read_end_time") + 1e-10)
		{
			std::cout << "Simulation end_time is greater than the read_end_time." << std::endl;
			std::cout << "end_time = " << es->parameters.get<double>("end_time") << std::endl;
			std::cout << "read_end_time " << es->parameters.get<double>("read_end_time") << std::endl;
			std::cout << "unsteady_from_steady = " << es->parameters.get<bool>("unsteady_from_steady") << std::endl;
			std::cout << "EXITING..." << std::endl;
			std::exit(0);
		}
	}

	set_double_parameter(infileparticle,"particle_entry_end_time",0.);

	// 0 - particle deposition all at once, 1 - particle deposition at a specified rate on a surface, 2 - particle deposition all at once until all particles deposited or exited, 3 - particle deposition at a constant rate until total number have entered system
	set_unsigned_int_parameter(infileparticle,"particle_deposition_type",1);
	set_unsigned_int_parameter(infileparticle,"particle_deposition_location",0);
	set_double_parameter(infileparticle,"particle_deposition_rate",10);
	set_unsigned_int_parameter(infileparticle,"particle_deposition_surface",0);
	set_double_parameter(infileparticle,"distance_deposited_inside",1.0);
	set_unsigned_int_parameter(infileparticle,"deposition_pattern",0);	// 0 is just a square (untested), 1 is circular uniform, 2 evenly distributed (only implemented for 2D so far), 3 is parabolic random, 4 is parabolic deterministic (2D only)
	set_double_parameter(infileparticle,"brownian_motion_magnitude",0.);
	set_unsigned_int_parameter(infileparticle,"deposit_within_radius",0.);

	//some parameters
	set_double_parameter(infileparticle,"particle_diameter",10.e-6);
	set_double_parameter(infileparticle,"particle_density",1200.);
	set_double_parameter(infileparticle,"particle_air_density",1.2);
	set_double_parameter(infileparticle,"particle_air_viscosity",1.8e-5);
	set_double_parameter(infileparticle,"mean_free_path",68.e-9);
	set_double_parameter(infileparticle,"gravity",9.81);

	set_bool_parameter(infileparticle,"particle_sedimentation",false);
	set_bool_parameter(infileparticle,"particle_impaction",false);
	set_bool_parameter(infileparticle,"particle_diffusion",false);
	set_bool_parameter(infileparticle,"particle_drag",false);

	set_unsigned_int_parameter(infileparticle,"gravity_type",0);

	set_double_parameter(infileparticle,"particle_velocity_units",1.);
	set_double_parameter(infileparticle,"initial_particle_velocity_scaling",0.99);
	set_bool_parameter(infileparticle,"particle_velocity_2d",false);

	if(threed && es->parameters.get<bool>("particle_velocity_2d"))
	{
		std::cout << "No need to use 2D particle velocity approximation in 0D pipes if 3D." << std::endl;
		std::cout << "Changing to use conventional velocity" << std::endl;	
		es->parameters.set<bool>("particle_velocity_2d") = false;
	}

	set_double_parameter(infileparticle,"gravity_x",0.);
	set_double_parameter(infileparticle,"gravity_y",0.);
	set_double_parameter(infileparticle,"gravity_z",1.);

	set_bool_parameter(infileparticle,"deposition_verbose",false);
	set_unsigned_int_parameter(infileparticle,"particle_tracking_method",2);	// default james non-verbose
	set_bool_parameter(infileparticle,"exit_on_lost_particle",false);
	set_bool_parameter(infileparticle,"output_3d_particles_readable",false);

	



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






void NavierStokesCoupled::init_particles()
{

	std::cout << "\nINITIALISING NEW 3D PARTICLES" << std::endl;

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
	const PointLocatorBase & pl = mesh_3d.point_locator();

	// figure out number of particles to deposit
	// particle_deposition_type 
	// 0 - all at once
	// 1 - at a rate
	// 2 - all at once and run simulation until all deposited or exited
	// 3 - constant rate until total number have entered system
	int num_new_particles = 0;
	if(es->parameters.get<unsigned int>("particle_deposition_type") == 0) 
		num_new_particles = (unsigned int) es->parameters.get<double>("particle_deposition_rate");
	else if(es->parameters.get<unsigned int>("particle_deposition_type") == 1)
		num_new_particles = (unsigned int) es->parameters.get<double>("particle_deposition_rate") * dt;
	else if(es->parameters.get<unsigned int>("particle_deposition_type") == 2)
		num_new_particles = (unsigned int) es->parameters.get<double>("particle_deposition_rate");
	else if(es->parameters.get<unsigned int>("particle_deposition_type") == 3
					|| es->parameters.get<unsigned int>("particle_deposition_type") == 4)
	{
		// we need to make sure we know a priori the number of particles that will enter in total
		// we just add extra particles in the last entry step, cheating, but shouldn't be too much of a difference
		
		if(es->parameters.get<double>("particle_entry_end_time") < 1e-10)
		{
			std::cout << "error, particle entry end time is zero." << std::endl;
			std::cout << "EXITING..." << std::endl;
			std::exit(0);
		}

		int num_particles_entered_so_far = (int)(particles_3D.size());
		//std::cout << "num_particles_entered_so_far = " << num_particles_entered_so_far << std::endl;
		int desired_number_of_particles_so_far = 
			(int)(es->parameters.get<double>("particle_deposition_rate") / es->parameters.get<double>("particle_entry_end_time") * time);
		//std::cout << "desired_number_of_particles_so_far = " << desired_number_of_particles_so_far << std::endl;
		num_new_particles = desired_number_of_particles_so_far - num_particles_entered_so_far;
		//std::cout << "num_new_particles = " << num_new_particles << std::endl;

		// if we are at the end time then we need to check that we are getting everything in and there aren't any rounding errors
		if(time > es->parameters.get<double>("particle_entry_end_time") - 1e-10)
		{
			num_new_particles = (int)(es->parameters.get<double>("particle_deposition_rate")) - num_particles_entered_so_far;
			//std::cout << "at end so num_new_particles = " << num_new_particles << std::endl;
		}

		// need to double check that we haven't deposited too many particles
		if(num_new_particles + num_particles_entered_so_far > (int)(es->parameters.get<double>("particle_deposition_rate")))
		{
			num_new_particles = (int)(es->parameters.get<double>("particle_deposition_rate")) - num_particles_entered_so_far;
			//std::cout << "tried to put too many so num_new_particles = " << num_new_particles << std::endl;
		}

		//std::cout << "total particles after this time step = " << num_particles_entered_so_far + num_new_particles << std::endl;
		//std::cout << "target number of particles = " << (int)(es->parameters.get<double>("particle_deposition_rate")) << std::endl;

		
	}
	else if(es->parameters.get<unsigned int>("particle_deposition_type") == 5)
	{
		// we want the particles to enter based sine wave
		// calculate the desired number of particles up until now
		// then subtract from particles entered so far to get how many we need this time step
		if(es->parameters.get<double>("particle_entry_end_time") < 1e-10)
		{
			std::cout << "error, particle entry end time is zero." << std::endl;
			std::cout << "EXITING..." << std::endl;
			std::exit(0);
		}

		int num_particles_entered_so_far = (int)(particles_3D.size());
		//std::cout << "num_particles_entered_so_far = " << num_particles_entered_so_far << std::endl;
		int desired_number_of_particles_so_far = (int)(es->parameters.get<double>("particle_deposition_rate")/2
														* ( - cos(2*M_PI/(2*es->parameters.get<double>("particle_entry_end_time")) * time) + 1));
		//std::cout << "desired_number_of_particles_so_far = " << desired_number_of_particles_so_far << std::endl;
		num_new_particles = desired_number_of_particles_so_far - num_particles_entered_so_far;
		//std::cout << "num_new_particles = " << num_new_particles << std::endl;

		// if we are at the end time then we need to check that we are getting everything in and there aren't any rounding errors
		if(time > es->parameters.get<double>("particle_entry_end_time") - 1e-10)
		{
			num_new_particles = (int)(es->parameters.get<double>("particle_deposition_rate")) - num_particles_entered_so_far;
			//std::cout << "at end so num_new_particles = " << num_new_particles << std::endl;
		}

		// need to double check that we haven't deposited too many particles
		if(num_new_particles + num_particles_entered_so_far > (int)(es->parameters.get<double>("particle_deposition_rate")))
		{
			num_new_particles = (int)(es->parameters.get<double>("particle_deposition_rate")) - num_particles_entered_so_far;
			//std::cout << "tried to put too many so num_new_particles = " << num_new_particles << std::endl;
		}

		//std::cout << "total particles after this time step = " << num_particles_entered_so_far + num_new_particles << std::endl;
		//std::cout << "target number of particles = " << (int)(es->parameters.get<double>("particle_deposition_rate")) << std::endl;
		
	}
	

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

			Particle new_particle(*es,*es_3d,new_point,element_found,particles_3D.size(),&perf_log_move);
			particles_3D.push_back(new_particle);

			Point new_point_2(0.9,0.5,0.);
	
			const Elem * element_found_2 = pl(new_point_2);
			if(element_found_2 == NULL)
			{
				std::cout << "particle not in mesh EXITING" << std::endl;
				std::exit(0);
			}

			Particle new_particle_2(*es,*es_3d,new_point_2,element_found_2,particles_3D.size(),&perf_log_move);
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
			Particle new_particle(*es,*es_3d,new_point,element_found,particles_3D.size(),&perf_log_move);
			particles_3D.push_back(new_particle);

			Point new_point_2(-0.3,0.2,0.5);

			const Elem * element_found_2 = pl(new_point_2);
			if(element_found_2 == NULL)
			{
				std::cout << "particle not in mesh EXITING" << std::endl;
				std::exit(0);
			}

			Particle new_particle_2(*es,*es_3d,new_point_2,element_found_2,particles_3D.size(),&perf_log_move);
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

		//std::cout << "centroid = " << centroid << std::endl;
		//std::cout << "normal = " << normal << std::endl;
		//std::cout << "max_radius = " << max_radius << std::endl;

		max_radius *= 1.0;


		//std::cout << "depositing " << num_new_particles << " particles." << std::endl;
		
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
		
			//std::cout << "vector_1 = " << vector_1 << std::endl;
			//std::cout << "vector_2 = " << vector_2 << std::endl;

			unsigned int num_failed_particles = 0;

			for(unsigned int i=0; i<num_new_particles; i++)
			{
				bool particle_deposited = false;
				while(!particle_deposited)
				{
					Point new_point;				

					// 3D
					// square random deposition pattern
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
					// uniform random circular pattern (to test)
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
					}// deterministic uniform
					else if(es->parameters.get<unsigned int>("deposition_pattern") == 2)
					{
						std::cout << "deterministic uniform deposition pattern not defined for 3D yet." << std::endl;
						std::exit(0);
					}// parabolic random
					else if(es->parameters.get<unsigned int>("deposition_pattern") == 3)
					{
						// according to wolfram, needs to be distributed as sqrt(r)cos(theta) to be uniform
						double theta = (static_cast<double>(rand())/static_cast<double>(RAND_MAX)) * 2 * M_PI;
						//double radius_scale = (static_cast<double>(rand())/static_cast<double>(RAND_MAX));

						// only choose radius quadratically
						// 1 - r^2
						double radius_scale = 0.;
						bool valid_particle = false;
						while(!valid_particle)
						{
							// first we get number between 0 and 1
							radius_scale = static_cast<double>(rand())/static_cast<double>(RAND_MAX); // rand 0 to 1

							// okay, so we want to accept and reject particles based on their distance from the wall, 
							// e.g. a particle at 1 is rejected 100% of the time, a particle at 0 is accepted 100% of the time
							// and a particle at 0<x<1 is accepted 1 - x^2 of the time

							double random_decision = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
							if(radius_scale < 1 - pow(random_decision,2.0))
								valid_particle = true;
						}



						// figure out what the radius is
						double radius_at_angle = deposition_boundary->get_radius_from_angle(theta);
						double x = sqrt(radius_scale) * radius_at_angle * cos(theta);
						double y = sqrt(radius_scale) * radius_at_angle * sin(theta);

						new_point = centroid + x*vector_1 + y*vector_2;

						//std::cout << "random parabolic deposition pattern not defined for 3D yet." << std::endl;
						//std::exit(0);
					}// parabolic deterministic
					else if(es->parameters.get<unsigned int>("deposition_pattern") == 4)
					{
						std::cout << "deterministic parabolic deposition pattern not defined for 3D yet (or ever)." << std::endl;
						std::exit(0);
					}


					/*
					std::cout << "new_point = " << new_point << std::endl;
					*/


					//distance_deposited_inside = 1.e-3;

					if(deposition_boundary->is_on_surface(new_point))
					{

						const Elem * element_found = pl(new_point - distance_deposited_inside*normal);
						if(element_found != NULL)
						{
							new_point = new_point - distance_deposited_inside*normal;
							//normal is outward facing and want to deposit the particle a little bit inside
							//std::cout << "hmmmm " << new_point << std::endl;
							Particle new_particle(*es,*es_3d,new_point,element_found,particles_3D.size(),&perf_log_move);
							particles_3D.push_back(new_particle);
							//std::cout << "particle deposited at " << new_point << std::endl;
							particle_deposited = true;
						}
						else
						{
							//std::cout << "particle should have been found but wasn't" << std::endl;
							num_failed_particles++;
						}
					}
					else
					{
						//std::cout << "particle is not on surface." << std::endl;
						num_failed_particles++;
					}
				}
			}
			//std::cout << "num_failed_particles = " << num_failed_particles << std::endl;
				
		}
		else	// 2D
		{

			//from normal get orthogonal vector
			Point vector_1(normal(1),-normal(0),0.);

			// calculate end points for variance free distribution
			//std::cout << "max_radius = " << max_radius << std::endl;
			Point y_0 = centroid - max_radius*vector_1;
			Point y_1 = centroid + max_radius*vector_1;
			double length = max_radius*2.0;
			double dx = length/(num_new_particles+1);

			for(unsigned int i=0; i<num_new_particles; i++)
			{
				bool particle_deposited = false;
				while(!particle_deposited)
				{
					Point new_point;

					// 2D
					// whatever, in 2D let's not do this deposition on a square crap
					if(es->parameters.get<unsigned int>("deposition_pattern") == 0 )
					{
						std::cout << "square deposition pattern not defined for 2D yet (or ever)." << std::endl;
						std::exit(0);
					}// random uniform
					else if(es->parameters.get<unsigned int>("deposition_pattern") == 1)
					{
						// random number between -1 and 1 multiplied by the radius to get -r to r
						double x = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * max_radius;
						new_point = centroid + x*vector_1;
					}
					// deposit evenly over the surface at equally spaced points
					else if(es->parameters.get<unsigned int>("deposition_pattern") == 2)
					{
						// we need the end points of the surface
						new_point = y_0 + (i+1)*dx*vector_1;

					}
					// random parabolic
					else if(es->parameters.get<unsigned int>("deposition_pattern") == 3)
					{

						double random_distance = 0.;
						bool valid_particle = false;
						while(!valid_particle)
						{
							// first we get number between 0 and 1
							random_distance = static_cast<double>(rand())/static_cast<double>(RAND_MAX); // rand 0 to 1

							// okay, so we want to accept and reject particles based on their distance from the wall, 
							// e.g. a particle at 1 is rejected 100% of the time, a particle at 0 is accepted 100% of the time
							// and a particle at 0<x<1 is accepted 1 - x^2 of the time

							double random_decision = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
							if(random_distance < 1 - pow(random_decision,2.0))
								valid_particle = true;
						}

						// randomly choose left or right
						double random_plus_minus = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
						double sign = 1.;
						if(random_plus_minus < 0.5)
							sign = -1.;

						double x = sign * random_distance * max_radius;
						new_point = centroid + x*vector_1;		

					}// deterministic uniform
					else if(es->parameters.get<unsigned int>("deposition_pattern") == 4)
					{

						std::cout << "deterministic uniform pattern not defined for 2D yet." << std::endl;
						std::exit(0);

						// need a [-1,1] mapping of distance
						// shift from 0 - 1 to -1 1
						double x = (i+1)*dx;
						std::cout << "x = " << x << std::endl;
						double x_squared = length*(2.*pow((x-length/2.)/length,2.0) + 1/2.);
						std::cout << "x_squared 1 = " << x_squared << std::endl;
						if(x < length/2. + 1e-10)
							x_squared = -2./length*x*(x-length);
						std::cout << "x_squared 2 = " << x_squared << std::endl;


						// we need the end points of the surface
						new_point = y_0 + x_squared * vector_1; //(i+1)*dx*vector_1 * x_squared;

					}
										



					/*
					std::cout << "new_point = " << new_point << std::endl;
					*/




				
					// we really shouldn't need this, but doesn't hurt to check
					if(deposition_boundary->is_on_surface(new_point))
					{

						const Elem * element_found = pl(new_point - distance_deposited_inside*normal);
						if(element_found != NULL)
						{
							new_point = new_point - distance_deposited_inside*normal;
							//normal is outward facing and want to deposit the particle a little bit inside
							Particle new_particle(*es,*es_3d,new_point,element_found,particles_3D.size(),&perf_log_move);
							particles_3D.push_back(new_particle);
							//std::cout << "particle deposited at " << new_point << std::endl;
							particle_deposited = true;
						}
						else
						{
							//std::cout << "particle should have been found but wasn't..." << std::endl;
						}
					}
				}
			}
				
		}

		
		std::cout << "number of new particles initialised = " << num_new_particles << std::endl;
		std::cout << "total particles initialised = " << particles_3D.size() << std::endl;
		

	}
	

}




double NavierStokesCoupled::init_0d_particles()
{

	std::cout << "\nINITIALISING NEW 0D PARTICLE FRACTION" << std::endl;

	// if we are doing a 0d only calculation we need to send it a fraction each time steps
	// remember this is actually time step 1, some particles could have been deposited at time step 0
	double fraction_added = 0;
	if(es->parameters.get<unsigned int>("particle_deposition_type") == 0) 
	{
		// only add particles at the beginning
		if(fabs(time) < 1e-10)
			fraction_added = 1.;
		else
			fraction_added = 0.;
	}
	else if(es->parameters.get<unsigned int>("particle_deposition_type") == 1)
	{
		// we don't want to add any more particles in thelast time step
		if(fabs(time - es->parameters.get<double>("particle_end_time")) < 1e-10)
			fraction_added = 0.;
		else
			fraction_added = dt / es->parameters.get<double>("particle_end_time");
	}
	else if(es->parameters.get<unsigned int>("particle_deposition_type") == 2)
	{
		// only add particles at the beginning
		if(fabs(time) < 1e-10)
			fraction_added = 1.;
		else
			fraction_added = 0.;
	}
	else if(es->parameters.get<unsigned int>("particle_deposition_type") == 3
					|| es->parameters.get<unsigned int>("particle_deposition_type") == 4)
	{
		// we need to make sure we know a priori the number of particles that will enter in total
		// we just add extra particles in the last entry step, cheating, but shouldn't be too much of a difference

		if(time < (es->parameters.get<double>("particle_entry_end_time") - dt + 1e-10))
			fraction_added = dt / es->parameters.get<double>("particle_entry_end_time");

	}
	else if(es->parameters.get<unsigned int>("particle_deposition_type") == 5)
	{

		// we want the particles to enter based sine wave
		// calculate the desired number of particles up until now
		// then subtract from particles entered so far to get how many we need this time step

		std::cout << "fraction added so far = " << total_fraction_added << std::endl;
		double desired_fraction_so_far = 1./2
														* ( - cos(2*M_PI/(2*es->parameters.get<double>("particle_entry_end_time")) * time) + 1);
		std::cout << "desired_fraction_so_far = " << desired_fraction_so_far << std::endl;
		fraction_added = desired_fraction_so_far - total_fraction_added;
		std::cout << "fraction to be added = " << fraction_added << std::endl;

		// if we are at the end time then we need to check that we are getting everything in and there aren't any rounding errors
		if(time > es->parameters.get<double>("particle_entry_end_time") - 1e-10)
		{
			fraction_added = 1. - total_fraction_added;
			std::cout << "at end so fraction_added = " << fraction_added << std::endl;
		}

		// need to double check that we haven't deposited too many particles
		if(fraction_added + total_fraction_added > 1. + 1e-10)
		{
			fraction_added = 1. - total_fraction_added;
			std::cout << "tried to put too many so fraction_added = " << fraction_added << std::endl;
		}

		std::cout << "total fraction after this time step = " << total_fraction_added + fraction_added << std::endl;
		std::cout << "target fraction = " << 1. << std::endl;






	}
	std::cout << "fraction added = " << fraction_added << std::endl;


	return fraction_added;

}



void NavierStokesCoupled::move_particles()
{
	std::cout << "\nMOVING 3D PARTICLES" << std::endl;

	// it's as easy as this
	for(unsigned int i=0; i<particles_3D.size(); i++)
	{		
		particles_3D[i].try_and_move();
	}

}


// write particle data file for this time step
void NavierStokesCoupled::write_3d_particles()
{

	// sometimes you may want to output the particle data in easily readable format..
	if(es->parameters.get<bool>("output_3d_particles_readable"))
	{
		output_3d_particle_data();
	}


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

	particle_vtk_output_file << "SCALARS airway_id float 1" << std::endl;
	particle_vtk_output_file << "LOOKUP_TABLE default" << std::endl;

	for(unsigned int i=0; i<particles_3D.size(); i++)
	{
		particle_vtk_output_file << particles_3D[i].get_airway_id() << " ";
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
void NavierStokesCoupled::output_3d_particle_data()
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
	particle_output_file << "exit_surface\texit_time\tentrance_time\tairway_id\n";


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
		unsigned int airway_id = particles_3D[i].get_airway_id();
		
		particle_output_file << pid << "\t";
		particle_output_file << t_step << "\t";
		particle_output_file << time << "\t";
		particle_output_file << x << "\t";
		particle_output_file << y << "\t";
		if(threed)
			particle_output_file << z << "\t";
		particle_output_file << exit_surface << "\t";
		particle_output_file << exit_time << "\t";
		particle_output_file << entrance_time << "\t";
		particle_output_file << airway_id << "\n";
	}

	particle_output_file.close();
}


// calculate and print some particle data for the timestep
void NavierStokesCoupled::calculate_3d_particle_data()
{

	unsigned int num_valid_particles = 0;
	unsigned int num_particles_on_wall = 0;
	unsigned int num_particles_exited = 0;
	unsigned int num_lost_particles = 0;
	unsigned int num_particles_in_motion = 0;
	unsigned int num_particles_close_to_wall = 0;
	unsigned int total_particles = particles_3D.size();

	// for sims with diff wall ids we want to record this
	// labelled as the negative, i.e. wall bdy -1 is in idx 1
	std::vector<unsigned int> num_particles_on_wall_vec;
	std::vector<unsigned int> num_particles_in_motion_vec;
	if(es->parameters.get<bool> ("gmsh_diff_wall_bdy_id"))
	{
		num_particles_on_wall_vec.resize(num_wall_bdy_ids+1);
		num_particles_in_motion_vec.resize(num_wall_bdy_ids+1);
	}

	// vector containing the number of particles deposited on each surface, 0 is the inflow bdy
	//std::cout << "boundary_ids.size() = " << boundary_ids.size() << std::endl;
	std::vector<unsigned int> num_particles_on_surface(boundary_ids.size());

	for(unsigned int i=0; i < particles_3D.size(); i++)
	{
		int exit_surface = particles_3D[i].get_exit_surface();
		double exit_time = particles_3D[i].get_exit_time();
		bool close_to_wall = particles_3D[i].was_deposited_close_to_wall();
		unsigned int airway_id = particles_3D[i].get_airway_id();

		if(exit_time > 0.)
		{

			// lost particles
			if(exit_surface == -1000)
			{
				// none
				num_lost_particles++;
			}
			else if(exit_surface < 0)	// on wall
			{
				num_particles_on_wall++;
				num_valid_particles++;
				if(es->parameters.get<bool> ("gmsh_diff_wall_bdy_id"))
					num_particles_on_wall_vec[(unsigned int)(-1*exit_surface)]++;
			}
			else	// on exit surface
			{
				num_particles_on_surface[(unsigned int)exit_surface]++;
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
			if(es->parameters.get<bool> ("gmsh_diff_wall_bdy_id"))
				num_particles_in_motion_vec[airway_id]++;
			num_valid_particles++;
		}
		
	}

	std::cout << "total particles created = " << total_particles << std::endl;
	std::cout << "total particles that will be created (probably) = " << (unsigned int)es->parameters.get<double>("particle_deposition_rate") << std::endl;
	std::cout << "num_valid_particles = " << num_valid_particles << std::endl;
	std::cout << "num_lost_particles = " << num_lost_particles << std::endl;
	std::cout << "num_particles_in_motion = " << num_particles_in_motion << std::endl;
	std::cout << "num_particles_exited = " << num_particles_exited << std::endl;
	std::cout << "num_particles_close_to_wall = " << num_particles_close_to_wall << std::endl;
	std::cout << "num_particles_on_wall = " << num_particles_on_wall << std::endl;

	// if diff wall bdy ids, output for each how much
	centrelines_deposition_fraction.resize(num_wall_bdy_ids+1,0.);
	centrelines_particle_fraction.resize(num_wall_bdy_ids+1,0.);
	if(es->parameters.get<bool> ("gmsh_diff_wall_bdy_id"))
	{
		std::cout << std::endl;
		for(unsigned int i=1; i<num_particles_on_wall_vec.size(); i++)
		{
			//std::cout << "num_particles_on_wall in airway " << i << " = " << num_particles_on_wall_vec[i] << std::endl;
			double deposition_fraction = (double)num_particles_on_wall_vec[i]/(double)es->parameters.get<double>("particle_deposition_rate");
			//std::cout << "fraction_on_wall in airway " << i << " = " << fraction << std::endl;	
			centrelines_deposition_fraction[i] = deposition_fraction;
			//std::cout << "and behold" << std::endl;
		}

		for(unsigned int i=1; i<num_particles_in_motion_vec.size(); i++)
		{
			//std::cout << "num_particles_on_wall in airway " << i << " = " << num_particles_on_wall_vec[i] << std::endl;
			double particle_fraction = (double)num_particles_in_motion_vec[i]/(double)es->parameters.get<double>("particle_deposition_rate");
			//std::cout << "fraction_on_wall in airway " << i << " = " << fraction << std::endl;	
			centrelines_particle_fraction[i] = particle_fraction;
			//std::cout << "and behold" << std::endl;
		}

		//std::cout << "lo" << std::endl;
		// put the data into the system object
		set_centrelines_deposition_data();
	}


	std::cout << std::endl;
	/*
	for(unsigned int i=0; i<num_particles_on_surface.size(); i++)
		std::cout << "num_particles_on_surface_" << i << " = " << num_particles_on_surface[i] << std::endl;
	*/

	std::cout << "total deposition fraction (relative to total particles that have been created):"
						 << (double)num_particles_on_wall/(double)total_particles << std::endl;

	std::cout << "total deposition fraction (relative to total particles that will be created probably):"
						 << (double)num_particles_on_wall/(double)es->parameters.get<double>("particle_deposition_rate") << std::endl;
	


	// make sure it is the correct size or the first time initialise to zero
	unsigned int num_bdy = 0.;
	if(sim_1d)
		num_bdy = airway_elem_id_starts.size();
	if(sim_3d)
		num_bdy = boundary_ids.size();


	total_fraction_exited_3d_surface.resize(num_bdy,0.);
	fraction_exited_3d_surface.resize(num_bdy,0.);
	for(unsigned int i=0; i<num_particles_on_surface.size(); i++)
	{
		double total_fraction = (double)num_particles_on_surface[i]/(double)es->parameters.get<double>("particle_deposition_rate");
		//std::cout << "total_fraction_exited_3d_surface_" << i << " = " << total_fraction << std::endl;


		double fraction = total_fraction - total_fraction_exited_3d_surface[i];	// difference from the last time step
		//std::cout << "fraction_exited_3d_surface_" << i << " = " << fraction << std::endl;


		total_fraction_exited_3d_surface[i] = total_fraction;
		fraction_exited_3d_surface[i] = fraction;
	}


	// set some of the global deposition metrics we will use
	particle_fraction_added_so_far = ((double) total_particles) / es->parameters.get<double>("particle_deposition_rate");
	particle_fraction_3d = ((double) num_particles_in_motion) / es->parameters.get<double>("particle_deposition_rate");
	deposition_fraction_3d = ((double) num_particles_on_wall) / es->parameters.get<double>("particle_deposition_rate");
	deposition_fraction_3d_max = 0;
	if(es->parameters.get<bool> ("gmsh_diff_wall_bdy_id"))
	{
		for(unsigned int i=0; i<centrelines_deposition_fraction.size(); i++)
		{
			if(centrelines_deposition_fraction[i] > deposition_fraction_3d_max)
				deposition_fraction_3d_max = centrelines_deposition_fraction[i];
		}
	}
	// looks ugly in file when not exactly zero
	if(fabs(deposition_fraction_3d_max) < 1e-10)
		deposition_fraction_3d_max = 0.;
	
	
	


	if(num_lost_particles > 0 && es->parameters.get<bool>("exit_on_lost_particle"))
	{
		std::cout << "there is a lost particle, so we are exiting now." << std::endl;
		std::cout << "GOODBYE" << std::endl;
		std::cout << "HAVEANICEDAY" << std::endl;
		std::exit(0);
	}

}


// print some particle data for the timestep
void NavierStokesCoupled::calculate_0d_particle_data()
{

	double total_deposition_fraction_this_time_step = james_particle_deposition_object.get_total_deposition_fraction_this_time_step();
	double total_exit_fraction_this_time_step = james_particle_deposition_object.get_total_exit_fraction_this_time_step();

	std::cout << "\ntotal_deposition_fraction_this_time_step = " << total_deposition_fraction_this_time_step << std::endl;
	std::cout << "total_exit_fraction_this_time_step = " << total_exit_fraction_this_time_step << std::endl;

	// calculate some max values.
	deposition_fraction_0d_max = 0;
	deposition_fraction_sed_0d_max = 0;
	deposition_fraction_imp_0d_max = 0;
	deposition_fraction_dif_0d_max = 0;
	terminal_exit_fraction_max = 0;


	for(unsigned int i=0; i<airway_data.size(); i++)
	{
		if(airway_data[i].get_deposition_fraction_total() > deposition_fraction_0d_max)
			deposition_fraction_0d_max = airway_data[i].get_deposition_fraction_total();

		if(airway_data[i].get_deposition_fraction_sed() > deposition_fraction_sed_0d_max)
			deposition_fraction_sed_0d_max = airway_data[i].get_deposition_fraction_sed();

		if(airway_data[i].get_deposition_fraction_imp() > deposition_fraction_imp_0d_max)
			deposition_fraction_imp_0d_max = airway_data[i].get_deposition_fraction_imp();

		if(airway_data[i].get_deposition_fraction_dif() > deposition_fraction_dif_0d_max)
			deposition_fraction_dif_0d_max = airway_data[i].get_deposition_fraction_dif();

		if(airway_data[i].get_terminal_exit_fraction() > terminal_exit_fraction_max)
			terminal_exit_fraction_max = airway_data[i].get_terminal_exit_fraction();
	}

	// looks bad in file when not exactly zero
	if(fabs(deposition_fraction_0d_max) < 1e-10)
		deposition_fraction_0d_max = 0.;

	if(fabs(deposition_fraction_sed_0d_max) < 1e-10)
		deposition_fraction_sed_0d_max = 0.;

	if(fabs(deposition_fraction_imp_0d_max) < 1e-10)
		deposition_fraction_imp_0d_max = 0.;

	if(fabs(deposition_fraction_dif_0d_max) < 1e-10)
		deposition_fraction_dif_0d_max = 0.;

	if(fabs(terminal_exit_fraction_max) < 1e-10)
		terminal_exit_fraction_max = 0.;

	deposition_fraction_sed_0d = james_particle_deposition_object.get_total_deposition_fraction_sed();
	deposition_fraction_imp_0d = james_particle_deposition_object.get_total_deposition_fraction_imp();
	deposition_fraction_dif_0d = james_particle_deposition_object.get_total_deposition_fraction_dif();


	deposition_fraction_0d = james_particle_deposition_object.get_total_deposition_fraction();
	terminal_exit_fraction = james_particle_deposition_object.get_total_exit_fraction();
	particle_fraction_0d = james_particle_deposition_object.get_total_particle_fraction();

	std::cout << "\ntotal_deposition_fraction_0d = " << deposition_fraction_0d << std::endl;
	std::cout << "total_terminal_exit_fraction = " << terminal_exit_fraction << std::endl;
	std::cout << "total_particle_fraction_0d = " << particle_fraction_0d << std::endl;
	std::cout << "\ntotal 0D (exit+deposition+particle) = " << deposition_fraction_0d + terminal_exit_fraction + particle_fraction_0d << std::endl << std::endl;
}


// print some particle data for the timestep
bool NavierStokesCoupled::no_particles_in_motion()
{

	for(unsigned int i=0; i < particles_3D.size(); i++)
	{
		double exit_time = particles_3D[i].get_exit_time();

		// if there is a single particle in motion, i.e. no exit time
		if(exit_time <= 0.)
		{
			return false;
		}
	}

	// check for all the deposition models that use unsteady james deposition
	if(particle_deposition >= 5)
	{
		double total_particle_fraction = james_particle_deposition_object.get_total_particle_fraction();
		std::cout << "total_particle_fraction = " << total_particle_fraction << std::endl;
		if(fabs(total_particle_fraction) > 1e-10)
			return false;
	}

	std::cout << "--- no particles in motion, ending sim after next output. soon soon. ---" << std::endl;
	return true;
	
}






void NavierStokesCoupled::set_centrelines_deposition_data()
{

	std::cout << "Setting centrelines deposition fraction" << std::endl;

	//add the radius variable
	// i'm sure we only need to do this once!
	//create an equation system to hold data like the radius of the element

	const DofMap& dof_map = system_centrelines->get_dof_map();

	MeshBase::const_element_iterator       el     =
		system_centrelines->get_mesh().active_local_elements_begin();
	const MeshBase::const_element_iterator end_el =
		system_centrelines->get_mesh().active_local_elements_end();

	std::vector<dof_id_type> dof_indices_dep_frac;
	std::vector<dof_id_type> dof_indices_part_frac;
	const unsigned int dep_frac_var = system_centrelines->variable_number ("deposition_fraction");
	const unsigned int part_frac_var = system_centrelines->variable_number ("particle_fraction");

	//std::cout << "centrelines_deposition_fraction.size() = " << centrelines_deposition_fraction.size() << std::endl;
	
	for ( ; el != end_el; ++el)
	{
		const Elem* elem = *el;
		dof_map.dof_indices (elem, dof_indices_dep_frac,dep_frac_var);
		dof_map.dof_indices (elem, dof_indices_part_frac,part_frac_var);
		const dof_id_type elem_id = elem->id();
		//std::cout << "elem_id = " << elem_id << std::endl;
		
		for(unsigned int i=0; i < dof_indices_dep_frac.size(); i++)
		{
			system_centrelines->solution->set(dof_indices_dep_frac[i], centrelines_deposition_fraction[elem_id+1]);
			system_centrelines->solution->set(dof_indices_part_frac[i], centrelines_particle_fraction[elem_id+1]);
		}

		// add the data that we can to the centrelines
		// TODO: flow_rate, velocity, stokes number; particle_fraction {need subdomain data from gmsh}
		// DONE: deposition_fraction.
		centreline_airway_data[elem_id].set_deposition_fraction_total(centrelines_deposition_fraction[elem_id+1]);
		centreline_airway_data[elem_id].set_total_particle_fraction(centrelines_particle_fraction[elem_id+1]);

	}

	//must close the vector after editing it

	system_centrelines->solution->close();
}




// writes centrelines data, possibly with deposition fraction details
void NavierStokesCoupled::write_centrelines()
{
	
	
	system_centrelines = &es_centrelines->get_system<TransientExplicitSystem> ("centrelines");
	

	ExodusII_IO_Extended exo(mesh_centrelines);

	std::vector<double> var_scalings_centrelines;
	var_scalings_centrelines.push_back(1.0);
	var_scalings_centrelines.push_back(1.0);
	if(particle_deposition)
	{
		var_scalings_centrelines.push_back(1.0);	// dep frac
		var_scalings_centrelines.push_back(1.0);	// part frac
	}
	exo.set_var_scalings(var_scalings_centrelines);

	if(!es->parameters.get<bool>("reynolds_number_calculation"))
	{
		if(!es->parameters.get<bool>("output_nondim"))
			exo.set_length_scale(es->parameters.get<double>("length_scale"));
	}

	std::vector<std::string> variables_centrelines;
	variables_centrelines.push_back("radius");
	variables_centrelines.push_back("radius_vis");

	if(particle_deposition)
	{
		variables_centrelines.push_back("deposition_fraction");
		variables_centrelines.push_back("particle_fraction");
	}

	exo.set_output_variables(variables_centrelines);

	//std::cout << "before write disc" << std::endl;

  std::ostringstream file_name;

  	// We write the file in the ExodusII format.
  	//file_name << "results/out_1D_viscosity"
	//  				<< es->parameters.get<Real> ("viscosity")

	std::ostringstream file_name_soln;

	if(!es->parameters.get<bool>("multiple_output_files"))
	{
		file_name << output_folder.str() << "out_centrelines";

		file_name_soln << file_name.str();

		file_name_soln	<< ".e";

		if(!first_1d_write)
		{
			exo.append(true);
		}
		else
			first_1d_write = false;

		// discontinuous can write the proc id
		exo.write_discontinuous_exodusII (file_name_soln.str(),
		                            *es_centrelines);

		//write_elem_pid_1d(exo);

		exo.write_time(t_step+1,time *time_scale_factor);
	}
	else
	{

		std::cout <<" hullo?" << std::endl;

		file_name << output_folder.str() << "out_centrelines";

		file_name_soln << file_name.str();

		file_name_soln	<< ".e-s."
		    << std::setw(4)
		    << std::setfill('0')
		    << std::right
		    << t_step;
		
		// discontinuous can write the proc id
		exo.write_discontinuous_exodusII (file_name_soln.str(),
		                            *es_centrelines);

		//write_elem_pid_1d(exo);

		exo.write_time(1,time *time_scale_factor);
	}


	std::cout << "EXODUSII output for timestep " << t_step
		<< " written to " << file_name_soln.str() << std::endl;

	std::ostringstream file_name_es;
	std::ostringstream file_name_mesh;

	std::cout << "hmm, backup_write_interval = " << es->parameters.get<Real>("backup_write_interval") << std::endl;
	if(es->parameters.get<bool>("output_backup_files") 
			&& ((!unsteady && !particle_deposition) || 
		(time - ((int)((time + 1e-10) /es->parameters.get<Real>("backup_write_interval")))*es->parameters.get<Real>("backup_write_interval") < dt - 1e-10)))
	{
		std::cout << "Writing backup files." << std::endl;
		file_name_es << file_name.str();
		file_name_es << "_es_";
		file_name_es << std::setw(4) << std::setfill('0') << t_step;
		file_name_es << ".xda";
		es_1d->write(file_name_es.str(), libMeshEnums::WRITE);

		file_name_mesh << file_name.str();
		file_name_mesh << "_mesh_";
		file_name_mesh << std::setw(4) << std::setfill('0') << t_step;
		file_name_mesh << ".xda";
		mesh_1d.write(file_name_mesh.str());

		std::cout << "Backup files written to " << file_name_mesh.str()
			<< " and " << file_name_es.str() << std::endl;
	}



}

// output the deposition per generation
/* DEPRECATED
void NavierStokesCoupled::output_per_generation_deposition()
{

	std::cout << "calculating per generation deposition" << std::endl;

	// need the maximum generation number
	// the max_generations should have been set for 3d and 0d
	if(sim_1d && sim_3d)
		max_generations = max_generations_3d + max_generations_0d;
	else if(sim_1d)
		max_generations = max_generations_0d;
	else if(sim_3d)
		max_generations = max_generations_3d;

	std::vector<double> deposition_fraction_per_generation(max_generations,0.);

	std::cout << "max_generations = " << max_generations << std::endl;

	// need to loop over the airways 0D and/or 2D and add to this data
	if(sim_1d)
	{
		std::cout << "Calculating 1d deposition fraction per generation." << std::endl;
		TransientExplicitSystem * system_deposition;
		// Get a reference to the deposition system object that we get the flow rate from.
		system_deposition = &es_1d->get_system<TransientExplicitSystem> ("Particle-Deposition-1D");

		const DofMap & dof_map_deposition = system_deposition->get_dof_map();
		std::vector<dof_id_type> dof_indices_p_total;
		const unsigned int p_total_var = system_deposition->variable_number ("p_total");

		// Get a constant reference to the mesh object.
		const MeshBase& mesh = es_1d->get_mesh();

		MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
		const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();


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



				const dof_id_type elem_id = elem->id();
		    dof_map_deposition.dof_indices (elem, dof_indices_p_total, p_total_var);

				// get the generation number
				unsigned int generation = airway_data[elem_id].get_generation();
				
				//std::cout << "generation = " << generation << std::endl;
				//std::cout << "dof_indices_p_total[0] = " << dof_indices_p_total[0] << std::endl;
				
	
				// should only be one really
				for(unsigned int i=0; i < dof_indices_p_total.size(); i++)
					deposition_fraction_per_generation[generation] += system_deposition->solution->el(dof_indices_p_total[i]);
	

			}
		}
	}// end 1d add to deposition vector



	// get the deposition fraction from the 2d/3d centrelines
	if(sim_3d && es->parameters.get<bool>("gmsh_diff_wall_bdy_id") && es->parameters.get<bool>("use_centreline_data"))
	{
		std::cout << "Calculating 3d deposition fraction per generation." << std::endl;

		const DofMap& dof_map = system_centrelines->get_dof_map();

		MeshBase::const_element_iterator       el     =
			system_centrelines->get_mesh().active_local_elements_begin();
		const MeshBase::const_element_iterator end_el =
			system_centrelines->get_mesh().active_local_elements_end();

		std::vector<dof_id_type> dof_indices_dep_frac;
		const unsigned int dep_frac_var = system_centrelines->variable_number ("deposition_fraction");

		
		//system_centrelines->solution->print();
		//std::cout << "centrelines_deposition_fraction.size() = " << centrelines_deposition_fraction.size() << std::endl;
		//std::cout << "1" << std::endl;
		
		for ( ; el != end_el; ++el)
		{
			const Elem* elem = *el;
			dof_map.dof_indices (elem, dof_indices_dep_frac,dep_frac_var);
			const dof_id_type elem_id = elem->id();
			
			//std::cout << "elem_id = " << elem_id << std::endl;
			//std::cout << "dep frac dof indices size = " <<dof_indices_dep_frac.size() << std::endl;
			//std::cout << "dof_indices_dep_frac[0] = " << dof_indices_dep_frac[0] << std::endl;
			//std::cout << "dep frac curr = " <<system_centrelines->solution->el(dof_indices_dep_frac[0]) << std::endl;
			

			// get the generation number
			unsigned int generation = centreline_airway_data[elem_id].get_generation();
		
			// should only be one really
			for(unsigned int i=0; i < dof_indices_dep_frac.size(); i++)
				deposition_fraction_per_generation[generation] += system_centrelines->solution->el(dof_indices_dep_frac[i]);
		}

		//system_centrelines->solution->print();
	}// end 3d add to deposition vector


	std::cout << "Outputting per generation deposition" << std::endl;





	// now output to file
	
	std::ostringstream output_data_file_name;
	std::ofstream deposition_output_file;


	output_data_file_name	<< output_folder.str()
			<< "deposition_per_generation_"
	    << std::setw(4)
	    << std::setfill('0')
	    << std::right
	    << t_step
			<< ".dat";
	
	deposition_output_file.open(output_data_file_name.str().c_str());

	

	//std::cout<< "boo" << std::endl;
	deposition_output_file << "# Deposition fraction per generation" << std::endl;
	// write geometry variables later when reading meta data file
	// write boundary conditions later with Picard class
	deposition_output_file << "# Results" << std::endl;
	deposition_output_file << "# generation";
	deposition_output_file <<	"\tdeposition_fraction";

	//unsigned int max_generations = *std::max_element(num_generations.begin(), num_generations.end());

	//std::cout << "max_generations = " << max_generations << std::endl;

	for(unsigned int i=0; i < max_generations; i++)
	{
		deposition_output_file << std::endl;
		deposition_output_file <<	i;
		//for(unsigned int j=0; j < num_1d_trees; j++)

		deposition_output_file <<	"\t" << deposition_fraction_per_generation[i];
		
				
	}
	

	deposition_output_file.close();

	std::cout << "Deposition per generation output to:\n " << output_data_file_name.str() << std::endl;

}
*/
















// output deposition data for each airway to file for post processing
void NavierStokesCoupled::output_0d_airway_deposition_data()
{
	std::cout << "Outputting 0d airway deposition data." << std::endl;





	// now output to file
	
	std::ostringstream output_data_file_name;
	std::ofstream deposition_output_file;


	output_data_file_name	<< output_folder.str()
			<< "airway_deposition_data_0d_"
	    << std::setw(4)
	    << std::setfill('0')
	    << std::right
	    << t_step
			<< ".dat";
	
	deposition_output_file.open(output_data_file_name.str().c_str());

	

	//std::cout<< "boo" << std::endl;
	deposition_output_file << "# Airway deposition data for 0d airway tree" << std::endl;
	// write geometry variables later when reading meta data file
	// write boundary conditions later with Picard class
	deposition_output_file << "# Results" << std::endl;
	// airway id data
	deposition_output_file << "# airway_id";
	deposition_output_file <<	"\tgeneration";
	// physical airway data - for this time step
	deposition_output_file <<	"\tflow_rate";
	deposition_output_file <<	"\tvelocity";
	deposition_output_file <<	"\tstokes_number";
	deposition_output_file <<	"\tdeposition_probability";
	deposition_output_file <<	"\tdeposition_probability_sed";
	deposition_output_file <<	"\tdeposition_probability_imp";
	deposition_output_file <<	"\tdeposition_probability_dif";
	// deposition data - totals at this time step
	deposition_output_file <<	"\tparticle_fraction";
	deposition_output_file <<	"\tterminal_exit_fraction";
	deposition_output_file <<	"\tdeposition_fraction";
	deposition_output_file <<	"\tdeposition_fraction_sed";
	deposition_output_file <<	"\tdeposition_fraction_imp";
	deposition_output_file <<	"\tdeposition_fraction_dif";
	deposition_output_file <<	"\tcumulative_particle_fraction";

	// okay
	//unsigned int max_generations = *std::max_element(num_generations.begin(), num_generations.end());

	//std::cout << "max_generations = " << max_generations << std::endl;

	for(unsigned int i=0; i < airway_data.size(); i++)
	{
		deposition_output_file << std::endl;
		// airway data
		deposition_output_file <<	i;
		deposition_output_file <<	"\t" << airway_data[i].get_generation();
		// physical data i.e. airflow and particle
		deposition_output_file <<	"\t" << airway_data[i].get_flow_rate();
		deposition_output_file <<	"\t" << airway_data[i].get_velocity();
		deposition_output_file <<	"\t" << airway_data[i].get_stokes_number();
		// deposition prob
		deposition_output_file <<	"\t" << airway_data[i].get_deposition_probability_total();
		deposition_output_file <<	"\t" << airway_data[i].get_deposition_probability_sed();
		deposition_output_file <<	"\t" << airway_data[i].get_deposition_probability_imp();
		deposition_output_file <<	"\t" << airway_data[i].get_deposition_probability_dif();
		// deposition data
		deposition_output_file <<	"\t" << airway_data[i].get_total_particle_fraction();
		deposition_output_file <<	"\t" << airway_data[i].get_terminal_exit_fraction();
		deposition_output_file <<	"\t" << airway_data[i].get_deposition_fraction_total();
		deposition_output_file <<	"\t" << airway_data[i].get_deposition_fraction_sed();
		deposition_output_file <<	"\t" << airway_data[i].get_deposition_fraction_imp();
		deposition_output_file <<	"\t" << airway_data[i].get_deposition_fraction_dif();
		// more data
		deposition_output_file <<	"\t" << airway_data[i].get_cumulative_particle_fraction();
		//for(unsigned int j=0; j < num_1d_trees; j++)

		
				
	}
	

	deposition_output_file.close();

	std::cout << "Deposition per generation output to:\n " << output_data_file_name.str() << std::endl;

}



// output deposition data for each airway to file for post processing
void NavierStokesCoupled::output_3d_airway_deposition_data()
{
	std::cout << "Outputting 3d airway deposition data." << std::endl;





	// now output to file
	
	std::ostringstream output_data_file_name;
	std::ofstream deposition_output_file;


	output_data_file_name	<< output_folder.str()
			<< "airway_deposition_data_3d_"
	    << std::setw(4)
	    << std::setfill('0')
	    << std::right
	    << t_step
			<< ".dat";
	
	deposition_output_file.open(output_data_file_name.str().c_str());

	

	//std::cout<< "boo" << std::endl;
	deposition_output_file << "# Airway deposition data for 3d airway tree" << std::endl;
	// write geometry variables later when reading meta data file
	// write boundary conditions later with Picard class
	deposition_output_file << "# Results" << std::endl;
	// airway id data
	deposition_output_file << "# airway_id";
	deposition_output_file <<	"\tgeneration";
	// physical airway data - for this time step
	deposition_output_file <<	"\tflow_rate";
	deposition_output_file <<	"\tvelocity";
	deposition_output_file <<	"\tstokes_number";
	deposition_output_file <<	"\tdeposition_probability";
	deposition_output_file <<	"\tdeposition_probability_sed";
	deposition_output_file <<	"\tdeposition_probability_imp";
	deposition_output_file <<	"\tdeposition_probability_dif";
	// deposition data - totals at this time step
	deposition_output_file <<	"\tparticle_fraction";
	deposition_output_file <<	"\tterminal_exit_fraction";
	deposition_output_file <<	"\tdeposition_fraction";
	deposition_output_file <<	"\tdeposition_fraction_sed";
	deposition_output_file <<	"\tdeposition_fraction_imp";
	deposition_output_file <<	"\tdeposition_fraction_dif";
	// more data
	deposition_output_file <<	"\tcumulative_particle_fraction";


	// okay
	//unsigned int max_generations = *std::max_element(num_generations.begin(), num_generations.end());

	//std::cout << "max_generations = " << max_generations << std::endl;

	for(unsigned int i=0; i < centreline_airway_data.size(); i++)
	{
		deposition_output_file << std::endl;
		// airway data
		deposition_output_file <<	i;
		deposition_output_file <<	"\t" << centreline_airway_data[i].get_generation();
		// physical data i.e. airflow and particle
		deposition_output_file <<	"\t" << 0; //centreline_airway_data[i].get_flow_rate(); //TODO
		deposition_output_file <<	"\t" << 0; //centreline_airway_data[i].get_velocity(); //TODO
		deposition_output_file <<	"\t" << 0; //centreline_airway_data[i].get_stokes_number(); //TODO
		// deposition prob
		deposition_output_file <<	"\t" << centreline_airway_data[i].get_deposition_probability_total();
		deposition_output_file <<	"\t" << 0; // //IMPOSSIBLE
		deposition_output_file <<	"\t" << 0; // //IMPOSSIBLE
		deposition_output_file <<	"\t" << 0; // //IMPOSSIBLE
		// deposition data
		deposition_output_file <<	"\t" << centreline_airway_data[i].get_total_particle_fraction(); //TODO -> done
		deposition_output_file <<	"\t" << 0; //centreline_airway_data[i].get_terminal_exit_fraction();	// UNNECESSARY
		deposition_output_file <<	"\t" << centreline_airway_data[i].get_deposition_fraction_total();
		deposition_output_file <<	"\t" << 0; // //IMPOSSIBLE
		deposition_output_file <<	"\t" << 0; // //IMPOSSIBLE
		deposition_output_file <<	"\t" << 0; // //IMPOSSIBLE	
		deposition_output_file <<	"\t" << 0; // not implemented yet	
				
	}
	

	deposition_output_file.close();

	std::cout << "Deposition per generation output to:\n " << output_data_file_name.str() << std::endl;

}


// output deposition data for each airway to file for post processing
void NavierStokesCoupled::output_global_deposition_metrics(bool header)
{
	std::cout << "Outputting global deposition metrics." << std::endl;

	if(header)
	{
		std::ostringstream global_deposition_metrics_file_name;

		global_deposition_metrics_file_name << output_folder.str() << "global_deposition_metrics.dat";
		//output_data_file_name << "results/out_viscosity"	<< es->parameters.set<Real> ("viscosity") << ".dat";
		global_deposition_metrics_file.open(global_deposition_metrics_file_name.str().c_str());

		global_deposition_metrics_file << "# Global deposition metrics" << std::endl;
		global_deposition_metrics_file << "#timestep";
		global_deposition_metrics_file << "\ttime";
		global_deposition_metrics_file << "\tpart_frac_so_far";
		global_deposition_metrics_file << "\tpart_frac_comb";
		global_deposition_metrics_file << "\tpart_frac_0d";
		global_deposition_metrics_file << "\tpart_frac_3d";
		global_deposition_metrics_file << "\tterm_exit_frac";
		global_deposition_metrics_file << "\tterm_exit_frac_max";
		global_deposition_metrics_file << "\tdep_frac_comb";
		global_deposition_metrics_file << "\tdep_frac_0d";
		global_deposition_metrics_file << "\tdep_frac_3d";
		global_deposition_metrics_file << "\tdep_frac_comb_max";
		global_deposition_metrics_file << "\tdep_frac_0d_max";
		global_deposition_metrics_file << "\tdep_frac_3d_max";
		global_deposition_metrics_file << "\tdep_frac_sed_0d";
		global_deposition_metrics_file << "\tdep_frac_sed_0d_max";
		global_deposition_metrics_file << "\tdep_frac_imp_0d";
		global_deposition_metrics_file << "\tdep_frac_imp_0d_max";
		global_deposition_metrics_file << "\tdep_frac_dif_0d";
		global_deposition_metrics_file << "\tdep_frac_dif_0d_max";
	}
	else
	{
		particle_fraction_combined = particle_fraction_0d + particle_fraction_3d;
		deposition_fraction_combined = deposition_fraction_0d + deposition_fraction_3d;
		deposition_fraction_combined_max = std::max(deposition_fraction_0d_max,deposition_fraction_3d_max);

		global_deposition_metrics_file << std::endl;
		global_deposition_metrics_file << t_step << "\t" 
																	<< time << "\t" 
																	<< particle_fraction_added_so_far << "\t"
																	<< particle_fraction_combined << "\t"
																	<< particle_fraction_0d << "\t"
																	<< particle_fraction_3d << "\t"
																	<< terminal_exit_fraction << "\t"
																	<< terminal_exit_fraction_max << "\t"
																	<< deposition_fraction_combined << "\t"
																	<< deposition_fraction_0d << "\t"
																	<< deposition_fraction_3d << "\t"
																	<< deposition_fraction_combined_max << "\t"
																	<< deposition_fraction_0d_max << "\t"
																	<< deposition_fraction_3d_max << "\t"
																	<< deposition_fraction_sed_0d << "\t"
																	<< deposition_fraction_sed_0d_max << "\t"
																	<< deposition_fraction_imp_0d << "\t"
																	<< deposition_fraction_imp_0d_max << "\t"
																	<< deposition_fraction_dif_0d << "\t"
																	<< deposition_fraction_dif_0d_max << "\t";
	}
}








