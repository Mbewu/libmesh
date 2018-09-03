
#include "james_particle_deposition.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

JamesParticleDeposition::JamesParticleDeposition ()

{

}

void JamesParticleDeposition::setup(EquationSystems& es_in,EquationSystems& es_1d_in, std::vector<Airway>& _airway_data, std::vector<unsigned int> _subdomains_1d, std::vector<SurfaceBoundary* > surface_boundaries, EquationSystems& es_1d_terminal_in, std::vector<std::vector<Airway> >& _terminal_airway_data, std::vector<unsigned int> _terminal_airway_to_airway_map)
{

	std::cout << "setting up james steady particle deposition model." << std::endl;

	es = &es_in;
	es_1d = &es_1d_in;
	es_1d_terminal = &es_1d_terminal_in;
	p_airway_data = &_airway_data;
	p_terminal_airway_data = &_terminal_airway_data;	// terminal airways
	subdomains_1d = _subdomains_1d;
	deposition_verbose = es->parameters.get<bool>("deposition_verbose");
	terminal_airway_to_airway_map = _terminal_airway_to_airway_map;

	// ***************** zero some variables ******** //

	// all airways
	total_deposition_fraction_this_time_step = 0.;
	total_exit_fraction_this_time_step = 0.;
	total_deposition_fraction = 0.;
	total_deposition_fraction_sed = 0.;
	total_deposition_fraction_imp = 0.;
	total_deposition_fraction_dif = 0.;
	total_exit_fraction = 0.;
	total_particle_fraction = 0.;

	// terminal airways
	total_terminal_deposition_fraction_this_time_step = 0.;
	total_terminal_exit_fraction_this_time_step = 0.;
	total_terminal_deposition_fraction = 0.;
	total_terminal_deposition_fraction_sed = 0.;
	total_terminal_deposition_fraction_imp = 0.;
	total_terminal_deposition_fraction_dif = 0.;
	total_terminal_exit_fraction = 0.;
	total_terminal_particle_fraction = 0.;

	// ****************** setup the particle deposition system and some vectors *********************** //
	TransientExplicitSystem * system_deposition;
  // Get a reference to the deposition system object that we get the flow rate from.
	system_deposition = &es_1d->get_system<TransientExplicitSystem> ("Particle-Deposition-1D");

	TransientExplicitSystem * system_deposition_terminal;
  // Get a reference to the deposition system object that we get the flow rate from.
	system_deposition_terminal = &es_1d_terminal->get_system<TransientExplicitSystem> ("Particle-Deposition-1D-Terminal");


	std::set<subdomain_id_type> active_subdomains;
	active_subdomains.clear(); 
	

	for(unsigned int i=0; i<_subdomains_1d.size(); i++)
		active_subdomains.insert(_subdomains_1d[i]);

	for(unsigned int i=0; i<_subdomains_1d.size(); i++)
		std::cout << "subdomains_1d = " << _subdomains_1d[i] << std::endl;


	const unsigned int p_airway_var = system_deposition->add_variable ("p_airway", CONSTANT, MONOMIAL,&active_subdomains);
	const unsigned int p_total_var = system_deposition->add_variable ("p_total", CONSTANT, MONOMIAL,&active_subdomains);
	const unsigned int p_airway_sed_var = system_deposition->add_variable ("p_airway_sed", CONSTANT, MONOMIAL,&active_subdomains);
	const unsigned int p_total_sed_var = system_deposition->add_variable ("p_total_sed", CONSTANT, MONOMIAL,&active_subdomains);
	const unsigned int p_airway_imp_var = system_deposition->add_variable ("p_airway_imp", CONSTANT, MONOMIAL,&active_subdomains);
	const unsigned int p_total_imp_var = system_deposition->add_variable ("p_total_imp", CONSTANT, MONOMIAL,&active_subdomains);
	const unsigned int p_airway_dif_var = system_deposition->add_variable ("p_airway_dif", CONSTANT, MONOMIAL,&active_subdomains);
	const unsigned int p_total_dif_var = system_deposition->add_variable ("p_total_dif", CONSTANT, MONOMIAL,&active_subdomains);
	const unsigned int particle_fraction_var = system_deposition->add_variable ("particle_fraction", CONSTANT, MONOMIAL,&active_subdomains);
	const unsigned int terminal_exit_fraction_var = system_deposition->add_variable ("terminal_exit_fraction", CONSTANT, MONOMIAL,&active_subdomains);

	// terminal, active on all subdomains
	const unsigned int p_terminal_airway_var = system_deposition_terminal->add_variable ("p_airway", CONSTANT, MONOMIAL);
	const unsigned int p_terminal_total_var = system_deposition_terminal->add_variable ("p_total", CONSTANT, MONOMIAL);
	const unsigned int p_terminal_airway_sed_var = system_deposition_terminal->add_variable ("p_airway_sed", CONSTANT, MONOMIAL);
	const unsigned int p_terminal_total_sed_var = system_deposition_terminal->add_variable ("p_total_sed", CONSTANT, MONOMIAL);
	const unsigned int p_terminal_airway_imp_var = system_deposition_terminal->add_variable ("p_airway_imp", CONSTANT, MONOMIAL);
	const unsigned int p_terminal_total_imp_var = system_deposition_terminal->add_variable ("p_total_imp", CONSTANT, MONOMIAL);
	const unsigned int p_terminal_airway_dif_var = system_deposition_terminal->add_variable ("p_airway_dif", CONSTANT, MONOMIAL);
	const unsigned int p_terminal_total_dif_var = system_deposition_terminal->add_variable ("p_total_dif", CONSTANT, MONOMIAL);
	const unsigned int terminal_particle_fraction_var = system_deposition_terminal->add_variable ("particle_fraction", CONSTANT, MONOMIAL);
	const unsigned int terminal_terminal_exit_fraction_var = system_deposition_terminal->add_variable ("terminal_exit_fraction", CONSTANT, MONOMIAL);



	airway_deposition_probability.resize(p_airway_data->size(),0.);
	airway_deposition_probability_sed.resize(p_airway_data->size(),0.);
	airway_deposition_probability_imp.resize(p_airway_data->size(),0.);
	airway_deposition_probability_dif.resize(p_airway_data->size(),0.);
	airway_flow_rate.resize(p_airway_data->size(),0.);

	// terminal
	if(es->parameters.get<unsigned int> ("additional_symmetric_generations"))
	{
		terminal_airway_deposition_probability.resize(p_terminal_airway_data->size(),std::vector<double>());
		terminal_airway_deposition_probability_sed.resize(p_terminal_airway_data->size(),std::vector<double>());
		terminal_airway_deposition_probability_imp.resize(p_terminal_airway_data->size(),std::vector<double>());
		terminal_airway_deposition_probability_dif.resize(p_terminal_airway_data->size(),std::vector<double>());
		terminal_airway_flow_rate.resize(p_terminal_airway_data->size(),std::vector<double>());

		for(unsigned int i=0; i<p_terminal_airway_data->size(); i++)
		{
			terminal_airway_deposition_probability[i].resize(p_terminal_airway_data->at(i).size(),0.);
			terminal_airway_deposition_probability_sed[i].resize(p_terminal_airway_data->at(i).size(),0.);
			terminal_airway_deposition_probability_imp[i].resize(p_terminal_airway_data->at(i).size(),0.);
			terminal_airway_deposition_probability_dif[i].resize(p_terminal_airway_data->at(i).size(),0.);
			terminal_airway_flow_rate[i].resize(p_terminal_airway_data->at(i).size(),0.);
		}
	}


	// **************** setup physical parameters ************************ //

	std::cout << "setting up global physical parameters" << std::endl;

	density = es->parameters.get<double>("density");
	particle_density = es->parameters.get<double>("particle_density");
	double lambda = es->parameters.get<double>("mean_free_path");
	particle_diameter =  es->parameters.get<double>("particle_diameter");
	cunningham_correction_factor = 1. + 2.*lambda/particle_diameter*(1.257 + 0.4*exp(-1.1*(particle_diameter/(2*lambda))));

	particle_radius = 0.5 * particle_diameter;
	
	gravity_magnitude = es->parameters.get<double>("gravity");
	viscosity = es->parameters.get<double>("viscosity");

	//from stokes-einstein relation
	double boltzmann_constant = 1.3806488e-23;
	diffusion_coefficient = boltzmann_constant * 310.15 * cunningham_correction_factor / (3*M_PI * viscosity * particle_diameter);

	double fluid_density = es->parameters.get<double>("density");
	sedimentation_velocity = 2./9.*(particle_density - fluid_density)*gravity_magnitude*pow(particle_radius,2.0)/viscosity;


	std::cout << "cunningham correction factor = " << cunningham_correction_factor << std::endl;
	std::cout << "viscosity = " << viscosity << std::endl;





	// ******************* setup branching angles and angle from gravity ****** //

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es_1d->get_mesh();


	//setup gravity direction vector
	Point gravity_direction(0.,0.,0.);
	// don't set it again if coupled, should have been set by the 3D part
	if(es->parameters.get<unsigned int>("particle_deposition") != 4
			|| es->parameters.get<unsigned int>("particle_deposition") != 6)
	{
		if(es->parameters.get<unsigned int>("gravity_type") == 0)
		{
			// gravity_direction in the z (3D) or y (2D) direction
			if(es->parameters.get<bool>("threed"))
				gravity_direction(2) = -1.;
			else
				gravity_direction(1) = -1.;	

			es->parameters.set<double>("gravity_x") = gravity_direction(0);
			es->parameters.set<double>("gravity_y") = gravity_direction(1);
			if(es->parameters.get<bool>("threed"))
				es->parameters.set<double>("gravity_z") = gravity_direction(2);


			std::cout << "gravity_direction = " << gravity_direction << std::endl;
			
		}
		else if(es->parameters.get<unsigned int>("gravity_type") == 1)
		{
			// user defined gravity_direction


			if(es->parameters.get<bool>("threed"))
			{
				std::cout << "hmmm" << std::endl;
				gravity_direction(0) = es->parameters.get<double>("gravity_x");
				gravity_direction(1) = es->parameters.get<double>("gravity_y");
				gravity_direction(2) = es->parameters.get<double>("gravity_z");
				std::cout << "gravity_direction = " << gravity_direction << std::endl;
			}
			else
			{
				std::cout << "hmmm" << std::endl;
				gravity_direction(0) = es->parameters.get<double>("gravity_x");
				gravity_direction(1) = es->parameters.get<double>("gravity_y");
				std::cout << "gravity_direction = " << gravity_direction << std::endl;	
			}
		}
		else if(es->parameters.get<unsigned int>("gravity_type") == 2)
		{
			// gravity_direction in the direction of the first element (i.e. hopefully the trachea)
			std::cout << "setting gravity direction based on first element" << std::endl;
			MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
			const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
		


			for ( ; el != end_el; ++el)
			{
				const Elem* elem = *el;
				if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
				{
					const int current_el_idx = elem->id();
					if(current_el_idx == 0)
					{
						Point p0 = *elem->get_node(0);
						Point p1 = *elem->get_node(1);
						Point direction = p1 - p0;
						direction = direction.unit();

						std::cout << "direction = " << direction << std::endl;



						if(es->parameters.get<bool>("threed"))
						{
							es->parameters.set<double>("gravity_x") = direction(0);
							es->parameters.set<double>("gravity_y") = direction(1);
							es->parameters.set<double>("gravity_z") = direction(2);
						}
						else
						{
							es->parameters.set<double>("gravity_x") = direction(0);
							es->parameters.set<double>("gravity_y") = direction(1);
						}

						gravity_direction = direction;
						std::cout << "gravity_direction = " << gravity_direction << std::endl;
						break;
					}
	
				}
			}
		}
	}
	else	// if coupled simulation, should have been set by 3D already
	{
		gravity_direction(0) = es->parameters.set<double>("gravity_x");
		gravity_direction(1) = es->parameters.set<double>("gravity_y");
		gravity_direction(2) = es->parameters.set<double>("gravity_z");
	}
	gravity_direction = gravity_direction.unit();


	// now that the gravity direction is set, calculate the branching angles and gravity angles
	
	// the gravity angle is relative to the horizontal, with no uphill or downhill differentiation so to calculate:
	// - calculate angle between pipe direction and gravity direction
	// - if < 90, gravity_angle = 90 - angle
	// - if > 90, gravity_angle = angle - 90;

	branching_angles.resize(p_airway_data->size(),0.);
	gravity_angles.resize(p_airway_data->size(),0.);

	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();


				
  for ( ; el != end_el; ++el)
  {
		const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
			

			//element data object starts numbering from 0
			//and the values referenced in it also do so need to take this into account

			const int current_el_idx = elem->id();
			Point p0 = *elem->get_node(0);
			Point p1 = *elem->get_node(1);
			Point direction = p1 - p0;
			direction = direction.unit();

			// if we have a parent calculate the branching angle
			if(p_airway_data->at(current_el_idx).has_parent())
			{
				unsigned int parent_id = p_airway_data->at(current_el_idx).get_parent();
				const Elem* parent_elem = mesh.elem(parent_id);
				Point parent_p0 = *parent_elem->get_node(0);
				Point parent_p1 = *parent_elem->get_node(1);
				Point parent_direction = parent_p1 - parent_p0;
				parent_direction = parent_direction.unit();

				double branching_angle = acos(direction*parent_direction);
				branching_angles[current_el_idx] = branching_angle;
				
			}
			else	// the first branch
			{
				// if we are running a coupled simulation, i.e. surface_boundaries in non-empty
				// we need to figure out what the angle between this and the 3D surface_boundary is
				if(surface_boundaries.size() > 0)
				{
					// the elem should have some info on which tree we are in
					// the boundary id of the 0D element is the same as the idx of the surface_boundary
					std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,0);
					if(boundary_ids.size() == 0)
					{
						std::cout << "error, no boundary id set for an element that has no parents." << std::endl;
						std::cout << "EXITING.." << std::endl;
						std::exit(0);
					}

					boundary_id_type boundary_id = boundary_ids[0];
					Point parent_direction = surface_boundaries[boundary_id]->get_normal();

					std::cout << "boundary_id = " << boundary_id << std::endl;
					std::cout << "direction = " << parent_direction << std::endl;


					double branching_angle = acos(direction*parent_direction);
					branching_angles[current_el_idx] = branching_angle;

					std::cout << "branching_angle = " << branching_angle << std::endl;
				}
				else
				{
					// just a 0D solution so first element doesn't have a branching angle
					// assuming it doesn't start with a bifurcation, unlikely and stupid
					branching_angles[current_el_idx] = 0.;
				}
			}

	
			// - calculate angle between pipe direction and gravity direction
			// - if == 90 ||, gravity_angle = 0
			// - if < 90, gravity_angle = 90 - angle
			// - if > 90, gravity_angle = angle - 90;

			// we need to check that direction*gravity_direction <= 1 and >= -1
			// this can happen when they are in the same direction and should equal 1, but numerically not
			double direction_times_gravity_direction = direction*gravity_direction;
			if(direction_times_gravity_direction > 1 - 1e-10)
				direction_times_gravity_direction = 1.0 - 1e-10;
			else if(direction_times_gravity_direction < -1 + 1e-10)
				direction_times_gravity_direction = -1.0 + 1e-10;

			double gravity_angle = acos(direction_times_gravity_direction);


			if((gravity_angle-M_PI/2.0 < 1e-10 && gravity_angle-M_PI/2.0 > -1e-10))
				gravity_angle = 0.;
			else if(gravity_angle < M_PI/2.0)
				gravity_angle = M_PI/2.0 - gravity_angle;
			else if(gravity_angle > M_PI/2.0)
				gravity_angle = gravity_angle - M_PI/2.0;

			gravity_angles[current_el_idx] = gravity_angle;			
		}
	}


	// if terminal, then branching and gravity angles are set up already




}



// calculate the deposition probability for each airway
// - also the stokes number
// - can be done indepdendently of connectivity
//   * loop over elements, based on flow rate in solution vector and airway_data calc deposition probability and save it to a solution vector and the airway_data object.
void JamesParticleDeposition::calculate_airway_deposition_probability()
{

	TransientExplicitSystem * system_airflow;
  // Get a reference to the airflow system object that we get the flow rate from.
	system_airflow = &es_1d->get_system<TransientExplicitSystem> ("ns1d_output");

	TransientExplicitSystem * system_deposition;
  // Get a reference to the deposition system object that we get the flow rate from.
	system_deposition = &es_1d->get_system<TransientExplicitSystem> ("Particle-Deposition-1D");

	// some stuff for inside the loop
  const DofMap & dof_map_airflow = system_airflow->get_dof_map();
  std::vector<dof_id_type> dof_indices_airflow;
  std::vector<dof_id_type> dof_indices_q;
  std::vector<dof_id_type> dof_indices_vel;
	const unsigned int q_var = system_airflow->variable_number ("Q");
	const unsigned int vel_var = system_airflow->variable_number ("vel");

  const DofMap & dof_map_deposition = system_deposition->get_dof_map();
  std::vector<dof_id_type> dof_indices_p_airway;
  std::vector<dof_id_type> dof_indices_p_airway_sed;
  std::vector<dof_id_type> dof_indices_p_airway_imp;
  std::vector<dof_id_type> dof_indices_p_airway_dif;
  std::vector<dof_id_type> dof_indices_p_total;
	const unsigned int p_airway_var = system_deposition->variable_number ("p_airway");
	const unsigned int p_airway_sed_var = system_deposition->variable_number ("p_airway_sed");
	const unsigned int p_airway_imp_var = system_deposition->variable_number ("p_airway_imp");
	const unsigned int p_airway_dif_var = system_deposition->variable_number ("p_airway_dif");

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



			const int current_el_idx = elem->id();
      dof_map_airflow.dof_indices (elem, dof_indices_airflow);
      dof_map_airflow.dof_indices (elem, dof_indices_q, q_var);
      dof_map_airflow.dof_indices (elem, dof_indices_vel, vel_var);

      dof_map_deposition.dof_indices (elem, dof_indices_p_airway, p_airway_var);
      dof_map_deposition.dof_indices (elem, dof_indices_p_airway_sed, p_airway_sed_var);
      dof_map_deposition.dof_indices (elem, dof_indices_p_airway_imp, p_airway_imp_var);
      dof_map_deposition.dof_indices (elem, dof_indices_p_airway_dif, p_airway_dif_var);



			// calculate common quantities

			// calculate the flow-rate and mean velocity, needed for most sims
			// the mean velocity is simply the flow rate divided by the area
			// for some reason, only the first dof has flow info
			double average_flow_rate = system_airflow->current_solution(dof_indices_q[0]);
			double pipe_radius = p_airway_data->at(current_el_idx).get_radius();
			double pipe_diameter = 2.0*pipe_radius;
			double area = M_PI*pipe_radius*pipe_radius;
			double mean_velocity = average_flow_rate/area;



			// if using 2d particle velocity calculate and set it
			if(es->parameters.get<bool>("particle_velocity_2d"))
			{
				mean_velocity = average_flow_rate/pipe_diameter;
			}



			double stokes_number = cunningham_correction_factor*particle_density*pow(particle_diameter,2.0)
												*mean_velocity/(18*viscosity*pipe_diameter);		//yeh1974

			// calculate deposition efficiencies
			double p_sed = 0.;
			double p_imp = 0.;
			double p_dif = 0.;


			if(es->parameters.get<bool>("particle_sedimentation"))
			{

				double pipe_angle_from_gravity = gravity_angles[current_el_idx];
				double pipe_length = p_airway_data->at(current_el_idx).get_length();
				// to calc time in branch, we assume that the particle is travelling with the air velocity
				double time_in_branch = pipe_length/mean_velocity;

				// calculate the sedimentation probability
				p_sed = calculate_sedimentation_probability(pipe_angle_from_gravity,pipe_radius,pipe_length,mean_velocity,time_in_branch);

			}

			if(es->parameters.get<bool>("particle_impaction"))
			{

				double branching_angle = branching_angles[current_el_idx];
				// to calc time in branch, we assume that the particle is travelling with the air velocity
				// airway reynolds number
				double reynolds_number = density * pipe_diameter * mean_velocity / viscosity;

				// calculate the sedimentation probability
				p_imp = calculate_impaction_probability(branching_angle,stokes_number,reynolds_number);

			}

			if(es->parameters.get<bool>("particle_diffusion"))
			{

				// calculate the sedimentation probability
				p_dif = 0.;

			}


			// add up the probabilities
			double p_airway = p_sed + p_imp + p_dif - p_sed*p_imp - p_sed*p_dif - p_imp*p_dif + p_imp*p_sed*p_dif;

			// calculate actual contributions of each mechanism (i.e. after subtraction of cross products)
			if(fabs(p_airway) > 1e-10)
			{
				p_sed = p_sed/(p_sed + p_imp + p_dif) * p_airway;
				p_imp = p_imp/(p_sed + p_imp + p_dif) * p_airway;
				p_dif = p_dif/(p_sed + p_imp + p_dif) * p_airway;
			}

			// put it in airway_data
			p_airway_data->at(current_el_idx).set_deposition_probability_total(p_airway);
			//airway_deposition_probability[current_el_idx] = p_airway;
			// put it in the solution vector
			for(unsigned int i=0; i < dof_indices_p_airway.size(); i++)
				system_deposition->solution->set(dof_indices_p_airway[i], p_airway);



			// sed
			if(es->parameters.get<bool>("particle_sedimentation"))
			{
				p_airway_data->at(current_el_idx).set_deposition_probability_sed(p_sed);
				//airway_deposition_probability_sed[current_el_idx] = p_sed;
				for(unsigned int i=0; i < dof_indices_p_airway_sed.size(); i++)
					system_deposition->solution->set(dof_indices_p_airway_sed[i], p_sed);
			}
			// imp
			if(es->parameters.get<bool>("particle_impaction"))
			{
				p_airway_data->at(current_el_idx).set_deposition_probability_imp(p_imp);
				//airway_deposition_probability_imp[current_el_idx] = p_imp;
				for(unsigned int i=0; i < dof_indices_p_airway_imp.size(); i++)
					system_deposition->solution->set(dof_indices_p_airway_imp[i], p_imp);
			}
			// dif
			if(es->parameters.get<bool>("particle_diffusion"))
			{
				p_airway_data->at(current_el_idx).set_deposition_probability_dif(p_dif);
				//airway_deposition_probability_dif[current_el_idx] = p_dif;
				for(unsigned int i=0; i < dof_indices_p_airway_dif.size(); i++)
					system_deposition->solution->set(dof_indices_p_airway_dif[i], p_dif);
			}





			// put the previous flow rate and velocity in their place
			double previous_flow_rate = p_airway_data->at(current_el_idx).get_flow_rate();
			double previous_velocity = p_airway_data->at(current_el_idx).get_velocity();
			p_airway_data->at(current_el_idx).set_previous_flow_rate(previous_flow_rate);
			p_airway_data->at(current_el_idx).set_previous_velocity(previous_velocity);


			// put flow rate, velocity and stokes number in airway_data for later
			//airway_flow_rate[current_el_idx] = average_flow_rate;

			for(unsigned int i=0; i<dof_indices_vel.size(); i++)
			{
				system_airflow->solution->set(dof_indices_vel[i],mean_velocity);
			}
			p_airway_data->at(current_el_idx).set_flow_rate(average_flow_rate);
			p_airway_data->at(current_el_idx).set_velocity(mean_velocity);
			p_airway_data->at(current_el_idx).set_stokes_number(stokes_number);
			
	

		}
	}

	system_deposition->solution->close();
	
}

/* DO LATER
// calculate the deposition probability for each terminal airway
// - also the stokes number
// - can be done indepdendently of connectivity
//   * loop over elements, based on flow rate in solution vector and airway_data calc deposition probability and save it to a solution vector and the airway_data object.
void JamesParticleDeposition::calculate_terminal_airway_deposition_probability()
{

	TransientExplicitSystem * system_airflow;
  // Get a reference to the airflow system object that we get the flow rate from to calculate the flow rates in the terminal branches
	system_airflow = &es_1d->get_system<TransientExplicitSystem> ("ns1d_output");

	TransientExplicitSystem * system_deposition_terminal;
  // Get a reference to the deposition system object that we get the flow rate from.
	system_deposition_terminal = &es_1d_terminal->get_system<TransientExplicitSystem> ("Particle-Deposition-1D-Terminal");

	// some stuff for inside the loop
  const DofMap & dof_map_airflow = system_airflow->get_dof_map();
  std::vector<dof_id_type> dof_indices_airflow;
  std::vector<dof_id_type> dof_indices_q;
  std::vector<dof_id_type> dof_indices_vel;
	const unsigned int q_var = system_airflow->variable_number ("Q");
	const unsigned int vel_var = system_airflow->variable_number ("vel");

  const DofMap & dof_map_deposition_terminal = system_deposition_terminal->get_dof_map();
	// finish later
  std::vector<dof_id_type> dof_indices_p_terminal_airway;
  std::vector<dof_id_type> dof_indices_p_terminal_airway_sed;
  std::vector<dof_id_type> dof_indices_p_terminal_airway_imp;
  std::vector<dof_id_type> dof_indices_p_terminal_airway_dif;
  std::vector<dof_id_type> dof_indices_p_terminal_total;
	const unsigned int p_terminal_airway_var = system_deposition_terminal->variable_number ("p_terminal_airway");
	const unsigned int p_terminal_airway_sed_var = system_deposition_terminal->variable_number ("p_terminal_airway_sed");
	const unsigned int p_terminal_airway_imp_var = system_deposition_terminal->variable_number ("p_terminal_airway_imp");
	const unsigned int p_terminal_airway_dif_var = system_deposition_terminal->variable_number ("p_terminal_airway_dif");

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



			const int current_el_idx = elem->id();

			// we need to check if it is a terminal airway
			if(!p_airway_data->at(current_el_idx).has_daughter_1())
			{


				// and we need to loop over the terminal airways corresponding to this airway

				// make sure that we have set this map up
				unsigned int::iterator terminal_airway_to_airway_map_iterator;
				terminal_airway_to_airway_map_iterator = std::find(terminal_airway_to_airway_map.begin(), terminal_airway_to_airway_map.end());
				if(terminal_airway_to_airway_map_iterator == terminal_airway_to_airway_map.end())
				{
					std::cout << "Error, airway id in terminal_airway_to_airway_map id not found..." << std::endl;
					std::cout << "Probably not set up correctly in NavierStokesCoupled::generate_1d_mesh." << std::endl;
					std::cout << "EXITING..." << std::endl;
					std::exit(0);
				}
				// otherwise just carry on

				// now we need to find the flow rate in this airway

		    dof_map_airflow.dof_indices (elem, dof_indices_airflow);
		    dof_map_airflow.dof_indices (elem, dof_indices_q, q_var);
				double parent_flow_rate = system_airflow->current_solution(dof_indices_q[0]);

				// loop over the terminal airways
				std::vector<Airway>& terminal_airway_data = p_terminal_airway_data->at(*terminal_airway_to_airway_map_iterator);
				for(unsigned int i=0; i<terminal_airway_data.size(); i++)
				{
					// calculate common quantities
					double average_flow_rate = parent_flow_rate/2.0;	//	pass it the correct flow rate and also
					double pipe_radius = terminal_airway_data[i].get_radius();
					double pipe_diameter = 2.0*pipe_radius;
					double area = M_PI*pipe_radius*pipe_radius;
					double mean_velocity = average_flow_rate/area;

					// if using 2d particle velocity calculate and set it
					if(es->parameters.get<bool>("particle_velocity_2d"))
					{
						mean_velocity = average_flow_rate/pipe_diameter;
					}

					double stokes_number = cunningham_correction_factor*particle_density*pow(particle_diameter,2.0)
														*mean_velocity/(18*viscosity*pipe_diameter);		//yeh1974

					// calculate deposition efficiencies
					double p_sed = 0.;
					double p_imp = 0.;
					double p_dif = 0.;

					if(es->parameters.get<bool>("particle_sedimentation"))
					{

						double pipe_angle_from_gravity = gravity_angles[current_el_idx];
						double pipe_length = terminal_airway_data[i].get_length();
						// to calc time in branch, we assume that the particle is travelling with the air velocity
						double time_in_branch = pipe_length/mean_velocity;

						// calculate the sedimentation probability
						p_sed = calculate_sedimentation_probability(pipe_angle_from_gravity,pipe_radius,pipe_length,mean_velocity,time_in_branch);

					}

					if(es->parameters.get<bool>("particle_impaction"))
					{

						double branching_angle = terminal_airway_data[i].get_branching_angles();
						// to calc time in branch, we assume that the particle is travelling with the air velocity
						// airway reynolds number
						double reynolds_number = density * pipe_diameter * mean_velocity / viscosity;

						// calculate the sedimentation probability
						p_imp = calculate_impaction_probability(branching_angle,stokes_number,reynolds_number);

					}

					if(es->parameters.get<bool>("particle_diffusion"))
					{

						// calculate the sedimentation probability
						p_dif = 0.;

					}


					// add up the probabilities
					double p_airway = p_sed + p_imp + p_dif - p_sed*p_imp - p_sed*p_dif - p_imp*p_dif + p_imp*p_sed*p_dif;

					// calculate actual contributions of each mechanism (i.e. after subtraction of cross products)
					if(fabs(p_airway) > 1e-10)
					{
						p_sed = p_sed/(p_sed + p_imp + p_dif) * p_airway;
						p_imp = p_imp/(p_sed + p_imp + p_dif) * p_airway;
						p_dif = p_dif/(p_sed + p_imp + p_dif) * p_airway;
					}

					// put it in airway_data
					//p_airway_data->at(current_el_idx).set_deposition_probability_total(p_airway);
					teminal_airway_data[current_el_idx].set_deposition_probability_total(p_airway);
					//airway_deposition_probability[current_el_idx] = p_airway;
					// put it in the solution vector
					for(unsigned int i=0; i < dof_indices_p_terminal_airway.size(); i++)
						system_deposition->solution->set(dof_indices_p_terminal_airway[i], p_airway);



					// sed
					if(es->parameters.get<bool>("particle_sedimentation"))
					{
						p_airway_data->at(current_el_idx).set_deposition_probability_sed(p_sed);
						//airway_deposition_probability_sed[current_el_idx] = p_sed;
						for(unsigned int i=0; i < dof_indices_p_terminal_airway_sed.size(); i++)
							system_deposition->solution->set(dof_indices_p_terminal_airway_sed[i], p_sed);
					}
					// imp
					if(es->parameters.get<bool>("particle_impaction"))
					{
						p_airway_data->at(current_el_idx).set_deposition_probability_imp(p_imp);
						//airway_deposition_probability_imp[current_el_idx] = p_imp;
						for(unsigned int i=0; i < dof_indices_p_terminal_airway_imp.size(); i++)
							system_deposition->solution->set(dof_indices_p_terminal_airway_imp[i], p_imp);
					}
					// dif
					if(es->parameters.get<bool>("particle_diffusion"))
					{
						p_airway_data->at(current_el_idx).set_deposition_probability_dif(p_dif);
						//airway_deposition_probability_dif[current_el_idx] = p_dif;
						for(unsigned int i=0; i < dof_indices_p_terminal_airway_dif.size(); i++)
							system_deposition->solution->set(dof_indices_p_terminal_airway_dif[i], p_dif);
					}





					// put the previous flow rate and velocity in their place
					double previous_flow_rate = p_airway_data->at(current_el_idx).get_flow_rate();
					double previous_velocity = p_airway_data->at(current_el_idx).get_velocity();
					p_airway_data->at(current_el_idx).set_previous_flow_rate(previous_flow_rate);
					p_airway_data->at(current_el_idx).set_previous_velocity(previous_velocity);


					// put flow rate, velocity and stokes number in airway_data for later
					//airway_flow_rate[current_el_idx] = average_flow_rate;

					for(unsigned int i=0; i<dof_indices_vel.size(); i++)
					{
						system_airflow->solution->set(dof_indices_vel[i],mean_velocity);
					}
					p_airway_data->at(current_el_idx).set_flow_rate(average_flow_rate);
					p_airway_data->at(current_el_idx).set_velocity(mean_velocity);
					p_airway_data->at(current_el_idx).set_stokes_number(stokes_number);
			
					



					parent_flow_rate = average_flow_rate;
				}





	

		}
	}

	system_deposition->solution->close();
	
}

*/




// calculate the total deposition fraction by taking into account the deposition fraction of the parent branch
// - loop over airway data
//   * get the deposition fraction from the parent branch to calc the exit fraction
//	 * divide the particles based on the flow rate or cross sectional area? flow rate gives the rate of volume coming through so makes sense
void JamesParticleDeposition::calculate_total_deposition_fraction(std::vector<double> total_fraction_exited_3d_surface, std::vector<std::vector<unsigned int> > airway_elem_id_starts)
{

	// *************** get some objects ***********

	TransientExplicitSystem * system_deposition;
  // Get a reference to the deposition system object that we get the flow rate from.
	system_deposition = &es_1d->get_system<TransientExplicitSystem> ("Particle-Deposition-1D");


	
	
	// need to tell it:
	// - what airway to do, 
	// - what particle density is coming through

	// zero the total deposition effiency and exit fraction
	total_deposition_fraction = 0.;
	total_deposition_fraction_sed = 0.;
	total_deposition_fraction_imp = 0.;
	total_deposition_fraction_dif = 0.;
	total_exit_fraction = 0.;


	// if we aren't given any info on which airway to start on, just start at 0
	if(airway_elem_id_starts.size() == 0)
	{
		std::cout << "airway_elem_id_starts vector is empty. fix it." << std::endl;
		std::cout << "EXITING..." << std::endl;
		std::exit(0);
	}
	else
	{
		// skip the first airway elem id because corresponds to the 3D inlet.
		for(unsigned int i=1; i<airway_elem_id_starts.size(); i++)
		{
			// get fraction from 3d and start at the correct airway
			double initial_particle_fraction = total_fraction_exited_3d_surface[i];

			// zero exit fraction means no need to calc anything
			if(fabs(initial_particle_fraction) > 1e-10)
			{
				if(airway_elem_id_starts[i].size() == 1)
				{
					calculate_airway_total_deposition_fraction(airway_elem_id_starts[i][0],initial_particle_fraction);
				}
				else
				{
					// calculate flow rate proportion of these branches
					std::vector<double> flow_rates;
					double total_flow_rate = 0.;
					for(unsigned int j=0; j<airway_elem_id_starts[i].size(); j++)
					{
						flow_rates.push_back(p_airway_data->at(airway_elem_id_starts[i][j]).get_flow_rate());
						total_flow_rate += flow_rates[j];
					}

					for(unsigned int j=0; j<airway_elem_id_starts[i].size(); j++)
					{
						double flow_rate_proportion = flow_rates[j]/total_flow_rate;
						calculate_airway_total_deposition_fraction(airway_elem_id_starts[i][j],flow_rate_proportion*initial_particle_fraction);
					}
				}
			}
		}
	}

	system_deposition->solution->close();
	

}







// calculate the total fraction of this airway 
// calculate the proportion of exit particles going to each daughter
// and then recursively loop to daughters
void JamesParticleDeposition::calculate_airway_total_deposition_fraction(unsigned int element_number, double entering_particle_density)
{
	
	// *************** get some objects ***********

	TransientExplicitSystem * system_deposition;
  // Get a reference to the deposition system object that we get the flow rate from.
	system_deposition = &es_1d->get_system<TransientExplicitSystem> ("Particle-Deposition-1D");

	const MeshBase& mesh = es_1d->get_mesh();



	// ******** calculate the deposition fraction in this airway
	double p_airway = p_airway_data->at(element_number).get_deposition_probability_total();//airway_deposition_probability[element_number];
	double p_airway_total = entering_particle_density * p_airway;

	// calculate the deposition due to sedimentation and the deposition due to impaction

	double p_airway_sed = p_airway_data->at(element_number).get_deposition_probability_sed(); //airway_deposition_probability_sed[element_number];
	double p_airway_total_sed = entering_particle_density * p_airway_sed;
	double p_airway_imp = p_airway_data->at(element_number).get_deposition_probability_imp(); //airway_deposition_probability_imp[element_number];
	double p_airway_total_imp = entering_particle_density * p_airway_imp;
	double p_airway_dif = p_airway_data->at(element_number).get_deposition_probability_dif(); //airway_deposition_probability_dif[element_number];
	double p_airway_total_dif = entering_particle_density * p_airway_dif;


	// ******** save the deposition in a solution vector and airway data vec
	// define some crap
  const DofMap & dof_map_deposition = system_deposition->get_dof_map();
  std::vector<dof_id_type> dof_indices_p_total;
  std::vector<dof_id_type> dof_indices_p_total_sed;
  std::vector<dof_id_type> dof_indices_p_total_imp;
  std::vector<dof_id_type> dof_indices_p_total_dif;
	const unsigned int p_total_var = system_deposition->variable_number ("p_total");
	const unsigned int p_total_sed_var = system_deposition->variable_number ("p_total_sed");
	const unsigned int p_total_imp_var = system_deposition->variable_number ("p_total_imp");
	const unsigned int p_total_dif_var = system_deposition->variable_number ("p_total_dif");

	// get the dofs for the element
	const Elem* elem = mesh.elem(element_number);
  dof_map_deposition.dof_indices (elem, dof_indices_p_total, p_total_var);
  dof_map_deposition.dof_indices (elem, dof_indices_p_total_sed, p_total_sed_var);
  dof_map_deposition.dof_indices (elem, dof_indices_p_total_imp, p_total_imp_var);
  dof_map_deposition.dof_indices (elem, dof_indices_p_total_dif, p_total_dif_var);

	// put it in the solution
	for(unsigned int i=0; i < dof_indices_p_total.size(); i++)
		system_deposition->solution->set(dof_indices_p_total[i], p_airway_total);

	for(unsigned int i=0; i < dof_indices_p_total_sed.size(); i++)
		system_deposition->solution->set(dof_indices_p_total_sed[i], p_airway_total_sed);
	for(unsigned int i=0; i < dof_indices_p_total_imp.size(); i++)
		system_deposition->solution->set(dof_indices_p_total_imp[i], p_airway_total_imp);
	for(unsigned int i=0; i < dof_indices_p_total_dif.size(); i++)
		system_deposition->solution->set(dof_indices_p_total_dif[i], p_airway_total_dif);

	// put it in the airway_data object
	p_airway_data->at(element_number).set_total_particle_fraction(entering_particle_density);	// slightly different from the defn for unsteady
	p_airway_data->at(element_number).set_deposition_fraction_total(p_airway_total);
	p_airway_data->at(element_number).set_deposition_fraction_sed(p_airway_total_sed);
	p_airway_data->at(element_number).set_deposition_fraction_imp(p_airway_total_imp);
	p_airway_data->at(element_number).set_deposition_fraction_dif(p_airway_total_dif);

	// add it to the total deposition fraction
	total_deposition_fraction += p_airway_total;
	total_deposition_fraction_sed += p_airway_total_sed;
	total_deposition_fraction_imp += p_airway_total_imp;
	total_deposition_fraction_dif += p_airway_total_dif;







	// ******** calculate the exit fraction to each daughter *********** //

	// total exit fraction
	double airway_exit_fraction = entering_particle_density - p_airway_total;

	// calculate the ratio between daughters
	// get the daughters
	std::vector<unsigned int> daughter_ids = p_airway_data->at(element_number).get_daughters();
	int daughter_1_id = -1;
	int daughter_2_id = -1;
	if(daughter_ids.size() > 0)
		daughter_1_id = daughter_ids[0];
	if(daughter_ids.size() > 1)
		daughter_2_id = daughter_ids[1];

	// has two daughter branches
	if(daughter_1_id != -1 && daughter_2_id != -1)
	{
		// calculate the flow rate
		double current_airway_flow_rate = airway_flow_rate[element_number];
		double daughter_1_flow_rate = airway_flow_rate[daughter_1_id];
		double daughter_2_flow_rate = airway_flow_rate[daughter_2_id];

		// proportion of total flow to each daughter
		double proportion_to_daughter_1 = daughter_1_flow_rate / current_airway_flow_rate;
		double proportion_to_daughter_2 = daughter_2_flow_rate / current_airway_flow_rate;

		double exit_fraction_1 = proportion_to_daughter_1 * airway_exit_fraction;
		double exit_fraction_2 = proportion_to_daughter_2 * airway_exit_fraction;


		// calculate deposition fraction in daughters
		calculate_airway_total_deposition_fraction(daughter_1_id, exit_fraction_1);
		calculate_airway_total_deposition_fraction(daughter_2_id, exit_fraction_2);
	}
	// has one daughter
	else if(daughter_1_id != -1)
	{
		// calculate deposition fraction in daughters
		calculate_airway_total_deposition_fraction(daughter_1_id, airway_exit_fraction);
	}
	else
	{
		// if no daughters then add the exit fraction to the total
		total_exit_fraction += airway_exit_fraction;

		// also add it to the terminal exit fraction of the airway
		p_airway_data->at(element_number).set_terminal_exit_fraction(airway_exit_fraction);
		
	}
	
	// recursive loop ends here if no daughters

}


// to initialise particle fraction before first time step, assumes 1 inlet and t=0
void JamesParticleDeposition::add_particle_fractions_to_tree(std::vector<double> particle_fractions, std::vector<std::vector<unsigned int> > airway_elem_id_starts, double _time)
{


	TransientExplicitSystem * system_deposition;

  // Get a reference to the deposition system object that we get the flow rate from.
	system_deposition = &es_1d->get_system<TransientExplicitSystem> ("Particle-Deposition-1D");

	const MeshBase& mesh = es_1d->get_mesh();

	// get some stuff we will use later for adding to the solution
	// define some crap
  const DofMap & dof_map_deposition = system_deposition->get_dof_map();
  std::vector<dof_id_type> dof_indices_particle_fraction;

	const unsigned int particle_fraction_var = system_deposition->variable_number ("particle_fraction");

	
	



	total_particle_fraction = 0;

	//set time and delta t
	time = _time;

	//std::cout << "airway_elem_id_starts.size() = " << airway_elem_id_starts.size() << std::endl;
	// skip the first airway elem id because corresponds to the 3D inlet.
	for(unsigned int i=1; i<airway_elem_id_starts.size(); i++)
	{
		// get fraction from 3d and start at the correct airway
		double initial_particle_fraction = particle_fractions[i];

		// add entering_particle_fraction if there is any
		// need to possibly split it up
		if(fabs(initial_particle_fraction) > 1e-10)
		{
			if(airway_elem_id_starts[i].size() == 1)
			{
				p_airway_data->at(airway_elem_id_starts[i][0]).add_particle_fraction(initial_particle_fraction,time);
				total_particle_fraction += initial_particle_fraction;

				// let's calc the total particle fraction in this airway and put it in the solution vec
				// initialise them for this element
				// get the dofs for the element
				const Elem* elem = mesh.elem(airway_elem_id_starts[i][0]);
				dof_map_deposition.dof_indices (elem, dof_indices_particle_fraction, particle_fraction_var);

				double total_airway_particle_fraction = p_airway_data->at(airway_elem_id_starts[i][0]).get_total_particle_fraction();
				// also add it to a solution vector
				for(unsigned int j=0; j < dof_indices_particle_fraction.size(); j++)
					system_deposition->solution->set(dof_indices_particle_fraction[j], total_airway_particle_fraction);
			}
			else // split between daughter branches according to flow rate
			{
				// calculate flow rate proportion of these branches
				std::vector<double> flow_rates;
				double total_flow_rate = 0.;
				for(unsigned int j=0; j<airway_elem_id_starts[i].size(); j++)
				{
					flow_rates.push_back(p_airway_data->at(airway_elem_id_starts[i][j]).get_flow_rate());
					total_flow_rate += flow_rates[j];
				}

				for(unsigned int j=0; j<airway_elem_id_starts[i].size(); j++)
				{
					double flow_rate_proportion = flow_rates[j]/total_flow_rate;
					p_airway_data->at(airway_elem_id_starts[i][j]).add_particle_fraction(flow_rate_proportion*initial_particle_fraction,time);
	
					total_particle_fraction += flow_rate_proportion * initial_particle_fraction;
					/*
					std::cout << "added " << flow_rate_proportion*initial_particle_fraction
										<< " fraction of particles, with entry time " 
										<< time << std::endl;
					std::cout << "to airway " << airway_elem_id_starts[i][j] << std::endl;
					*/


					// let's calc the total particle fraction in this airway and put it in the solution vec
					// initialise them for this element
					// get the dofs for the element
					const Elem* elem = mesh.elem(airway_elem_id_starts[i][j]);
					dof_map_deposition.dof_indices (elem, dof_indices_particle_fraction, particle_fraction_var);

					double total_airway_particle_fraction = p_airway_data->at(airway_elem_id_starts[i][j]).get_total_particle_fraction();
					// also add it to a solution vector
					for(unsigned int j=0; j < dof_indices_particle_fraction.size(); j++)
						system_deposition->solution->set(dof_indices_particle_fraction[j], total_airway_particle_fraction);


				}
				
			}
			
		}
	}

	system_deposition->solution->close();
}


// calculate the total deposition fraction by taking into account the deposition fraction of the parent branch
// - loop over airway data
//   * get the deposition fraction from the parent branch to calc the exit fraction
//	 * divide the particles based on the flow rate or cross sectional area? flow rate gives the rate of volume coming through so makes sense
// this function also adds up the total particle fraction in motion
void JamesParticleDeposition::calculate_time_step_deposition_fraction(std::vector<std::vector<unsigned int> > airway_elem_id_starts,double _time, double _dt)
{

	// *************** get some objects ***********

	TransientExplicitSystem * system_deposition;
  // Get a reference to the deposition system object that we get the flow rate from.
	system_deposition = &es_1d->get_system<TransientExplicitSystem> ("Particle-Deposition-1D");

	//set time and delta t for the loop
	time = _time;
	dt = _dt;

	
	
	// need to tell it:
	// - what airway to do, 
	// - what particle density is coming through

	// zero the total deposition effiency and exit fraction this time step
	// in the loop we calculate the total particle fraction present
	total_deposition_fraction_this_time_step = 0.;
	total_exit_fraction_this_time_step = 0.;
	total_particle_fraction = 0.;


	//std::cout << "airway_elem_id_starts.size() = " << airway_elem_id_starts.size() << std::endl;
	// skip the first airway elem id because corresponds to the 3D inlet.
	for(unsigned int i=1; i<airway_elem_id_starts.size(); i++)
	{
		//std::cout << "hello" << std::endl;
		// loop over each one
		for(unsigned int j=0; j<airway_elem_id_starts[i].size(); j++)
			calculate_airway_time_step_deposition_fraction(airway_elem_id_starts[i][j]);
	}

	system_deposition->solution->close();
}



// calculate the total fraction of this airway 
// calculate the proportion of exit particles going to each daughter
// and then recursively loop to daughters
void JamesParticleDeposition::calculate_airway_time_step_deposition_fraction(unsigned int element_number)
{
	
	// *************** get some objects ***********

	TransientExplicitSystem * system_deposition;
  // Get a reference to the deposition system object that we get the flow rate from.
	system_deposition = &es_1d->get_system<TransientExplicitSystem> ("Particle-Deposition-1D");

	const MeshBase& mesh = es_1d->get_mesh();

	// get some stuff we will use later for adding to the solution
	// define some crap
  const DofMap & dof_map_deposition = system_deposition->get_dof_map();
  std::vector<dof_id_type> dof_indices_p_total;
  std::vector<dof_id_type> dof_indices_p_total_sed;
  std::vector<dof_id_type> dof_indices_p_total_imp;
  std::vector<dof_id_type> dof_indices_p_total_dif;
  std::vector<dof_id_type> dof_indices_particle_fraction;
  std::vector<dof_id_type> dof_indices_terminal_exit_fraction;

	const unsigned int p_total_var = system_deposition->variable_number ("p_total");
	const unsigned int p_total_sed_var = system_deposition->variable_number ("p_total_sed");
	const unsigned int p_total_imp_var = system_deposition->variable_number ("p_total_imp");
	const unsigned int p_total_dif_var = system_deposition->variable_number ("p_total_dif");
	const unsigned int particle_fraction_var = system_deposition->variable_number ("particle_fraction");
	const unsigned int terminal_exit_fraction_var = system_deposition->variable_number ("terminal_exit_fraction");

	// initialise them for this element
	// get the dofs for the element
	const Elem* elem = mesh.elem(element_number);
	dof_map_deposition.dof_indices (elem, dof_indices_p_total, p_total_var);
	dof_map_deposition.dof_indices (elem, dof_indices_p_total_sed, p_total_sed_var);
	dof_map_deposition.dof_indices (elem, dof_indices_p_total_imp, p_total_imp_var);
	dof_map_deposition.dof_indices (elem, dof_indices_p_total_dif, p_total_dif_var);
	dof_map_deposition.dof_indices (elem, dof_indices_particle_fraction, particle_fraction_var);
	dof_map_deposition.dof_indices (elem, dof_indices_terminal_exit_fraction, terminal_exit_fraction_var);
	
	





	// if there are, handle their deposition, and sending particle fractions to daughters etc
	// then afterwards, we will move to the daughters

	// get the daughters
	std::vector<unsigned int> daughter_ids = p_airway_data->at(element_number).get_daughters();
	int daughter_1_id = -1;
	int daughter_2_id = -1;

	if(daughter_ids.size() > 0)
		daughter_1_id = daughter_ids[0];
	if(daughter_ids.size() > 1)
		daughter_2_id = daughter_ids[1];

	// check if there are particle fractions on this airway
	unsigned int num_particle_fractions = p_airway_data->at(element_number).num_particle_fractions();

	if(deposition_verbose)
	{
		std::cout << "\nelement_number = " << element_number << std::endl;
		std::cout << "num_particle_fractions = " << num_particle_fractions << std::endl;
	}

	if(num_particle_fractions)
	{
		// calculate the distance travelled in this time step
		double area = p_airway_data->at(element_number).get_area();
		double radius = p_airway_data->at(element_number).get_radius();
		double velocity = p_airway_data->at(element_number).get_velocity();
		double flow_rate = p_airway_data->at(element_number).get_flow_rate();
		double previous_flow_rate = p_airway_data->at(element_number).get_previous_flow_rate();
		
		// now we need to loop over the particle fractions and calculate their total distance travelled
		// we need to loop over this in a way such that removing does not affect it
		// so later, when we remove, we iterate back i--
		// if a particle arrived within the last time step, we only want to move them for that amount of time.
		for(unsigned int i=0; i<p_airway_data->at(element_number).num_particle_fractions(); i++)
		{

			double entry_time = p_airway_data->at(element_number).get_entry_time(i);
	
			// find out if particles arrived in previous timestep
			double particle_dt = std::min(time-entry_time,dt);

			if(deposition_verbose)
				std::cout << "particle_dt = " << particle_dt << std::endl;
			
			double distance_this_step = velocity * particle_dt;
			
			// get the distance so far
			double distance_so_far = p_airway_data->at(element_number).get_distance(i);

			// total distance
			double total_distance = distance_so_far + distance_this_step;

			// check if we have gone past the end of the length of the airway
			double length = p_airway_data->at(element_number).get_length();


			if(deposition_verbose)
			{
				std::cout << "radius = " << radius << std::endl;
				std::cout << "area = " << area << std::endl;
				std::cout << "flow_rate = " << flow_rate << std::endl;
				std::cout << "velocity = " << velocity << std::endl;
				std::cout << "total distance = " << total_distance << std::endl;
			}


			if(total_distance >= length)
			{
				// now we need to calculate:
				// -  how far we went past
				// - when we should have left
				// - what the mean velocity was
				// - how much we should deposit
				// - then we need to deposit the amount in the airway
				// - how much we should send to the daughters
				// - then we send it to the daughters
				// - then remove it from current airway list

				// how far we went past
				// distance_so_far + small_dt * velocity = L
				// small_dt = (L - distance_so_far) / velocity
				double small_dt = (length - distance_so_far) / velocity;
				
				// exit time
				double exit_time = time - particle_dt + small_dt;

				// mean velocity
				double mean_velocity = length/(exit_time - entry_time);

				if(deposition_verbose)
				{				
					std::cout << "small_dt = " << small_dt << std::endl;
					std::cout << "exit_time = " << exit_time << std::endl;
					std::cout << "entry_time = " << entry_time << std::endl;
					std::cout << "distance_so_far = " << distance_so_far << std::endl;
				}


				// calc dep probabilities
				double p_airway = 0;
				double p_airway_sed = 0;
				double p_airway_imp = 0;
				double p_airway_dif = 0;
				get_airway_deposition_probability(element_number, mean_velocity, p_airway, p_airway_sed, p_airway_imp, p_airway_dif);

				if(deposition_verbose)
				{				
					std::cout << "mean_velocity = " << mean_velocity << std::endl;
					std::cout << "p_airway = " << p_airway << std::endl;
				}
				// deposition and exit amount
				double particle_fraction = p_airway_data->at(element_number).get_particle_fraction(i);
				double deposition_fraction = p_airway * particle_fraction;
				double exit_fraction = particle_fraction - deposition_fraction;

				if(deposition_verbose)
					std::cout << "particle_fraction = " << particle_fraction << std::endl;

				// particle fractions for each mechanism
				double p_airway_total_sed = p_airway_sed * particle_fraction;
				double p_airway_total_imp = p_airway_imp * particle_fraction;
				double p_airway_total_dif = p_airway_dif * particle_fraction;

				// deposit the amount in the airway
				// let's do it in a solution vector and in the airway object
				// could probably do it more efficiently but hey ho
				p_airway_data->at(element_number).add_deposition_fraction_total(deposition_fraction);
				p_airway_data->at(element_number).add_deposition_fraction_sed(p_airway_total_sed);
				p_airway_data->at(element_number).add_deposition_fraction_imp(p_airway_total_imp);
				p_airway_data->at(element_number).add_deposition_fraction_dif(p_airway_total_dif);

				total_deposition_fraction_this_time_step += deposition_fraction;
				total_deposition_fraction += deposition_fraction;
				total_deposition_fraction_sed += p_airway_total_sed;
				total_deposition_fraction_imp += p_airway_total_imp;
				total_deposition_fraction_dif += p_airway_total_dif;

				if(deposition_verbose)
					std::cout << "deposited_fraction = " << deposition_fraction << std::endl;

				// put it in the solution
				for(unsigned int j=0; j < dof_indices_p_total.size(); j++)
					system_deposition->solution->add(dof_indices_p_total[j], deposition_fraction);

				for(unsigned int j=0; j < dof_indices_p_total_sed.size(); j++)
					system_deposition->solution->add(dof_indices_p_total_sed[j], p_airway_total_sed);
				for(unsigned int j=0; j < dof_indices_p_total_imp.size(); j++)
					system_deposition->solution->add(dof_indices_p_total_imp[j], p_airway_total_imp);
				for(unsigned int j=0; j < dof_indices_p_total_dif.size(); j++)
					system_deposition->solution->add(dof_indices_p_total_dif[j], p_airway_total_dif);




				// ******** remove the particle fraction from the current airway ********** //
				p_airway_data->at(element_number).remove_particle_fraction(i);
				// iterate back so that removing does not skip checking all the particle fractions
				i--;



				// ******** calculate the exit fraction to each daughter *********** //



				// has two daughter branches
				if(daughter_1_id != -1 && daughter_2_id != -1)
				{

					double daughter_1_flow_rate = p_airway_data->at(daughter_1_id).get_flow_rate();
					double daughter_1_previous_flow_rate = p_airway_data->at(daughter_1_id).get_previous_flow_rate();
					double daughter_2_flow_rate = p_airway_data->at(daughter_2_id).get_flow_rate();
					double daughter_2_previous_flow_rate = p_airway_data->at(daughter_2_id).get_previous_flow_rate();

					// calculate the flow rate at the exit time through linear interpolation
					double exit_flow_rate = previous_flow_rate + small_dt * (flow_rate - previous_flow_rate);
					double daughter_1_exit_flow_rate = daughter_1_previous_flow_rate + small_dt * (daughter_1_flow_rate - daughter_1_previous_flow_rate);
					double daughter_2_exit_flow_rate = daughter_2_previous_flow_rate + small_dt * (daughter_2_flow_rate - daughter_2_previous_flow_rate);

					// proportion of total flow to each daughter
					double proportion_to_daughter_1 = daughter_1_exit_flow_rate / exit_flow_rate;
					double proportion_to_daughter_2 = daughter_2_exit_flow_rate / exit_flow_rate;

					double exit_fraction_1 = proportion_to_daughter_1 * exit_fraction;
					double exit_fraction_2 = proportion_to_daughter_2 * exit_fraction;

					// send the particle fractions to the daughters
					if(fabs(exit_fraction_1) > 1e-10)
						p_airway_data->at(daughter_1_id).add_particle_fraction(exit_fraction_1,exit_time);

					// send the particle fractions to the daughters
					if(fabs(exit_fraction_2) > 1e-10)
						p_airway_data->at(daughter_2_id).add_particle_fraction(exit_fraction_2,exit_time);

				}
				// has one daughter, just send it all through
				else if(daughter_1_id != -1)
				{

					// send the particle fractions to the daughters
					if(fabs(exit_fraction) > 1e-10)
						p_airway_data->at(daughter_1_id).add_particle_fraction(exit_fraction,exit_time);

				}
				else
				{
					// if no daughters then add the exit fraction to the total
					total_exit_fraction += exit_fraction;
					total_exit_fraction_this_time_step += exit_fraction;

					
					if(deposition_verbose)
						std::cout << "terminal_exited_fraction = " << exit_fraction << std::endl;		

					// also add it to the exit fraction of this terminal branch
					p_airway_data->at(element_number).add_terminal_exit_fraction(exit_fraction);

					// also add it to a solution vector
					for(unsigned int j=0; j < dof_indices_terminal_exit_fraction.size(); j++)
						system_deposition->solution->add(dof_indices_terminal_exit_fraction[j], exit_fraction);
				}

			}
			else	// particles haven't left yet, add them to the 
			{
				// here, we need to update the particle fraction's distance
				p_airway_data->at(element_number).set_distance(i,total_distance);

				// we also need to add to the particle fraction in motion
				total_particle_fraction += p_airway_data->at(element_number).get_particle_fraction(i);

				if(deposition_verbose)
				{				
					std::cout << "distance = " << total_distance << std::endl;
					std::cout << "fraction added to total_particle_fraction  = " << p_airway_data->at(element_number).get_particle_fraction(i) << std::endl;
				}

			}
		}
	}

	// let's calc the total particle fraction in this airway and put it in the solution vec

	double total_airway_particle_fraction = p_airway_data->at(element_number).get_total_particle_fraction();
	// also add it to a solution vector
	for(unsigned int j=0; j < dof_indices_particle_fraction.size(); j++)
		system_deposition->solution->set(dof_indices_particle_fraction[j], total_airway_particle_fraction);


	// now we move onto the daughters, or not..
	
	// has two daughter branches
	if(daughter_1_id != -1 && daughter_2_id != -1)
	{
		// calculate deposition fraction in daughters
		calculate_airway_time_step_deposition_fraction(daughter_1_id);
		calculate_airway_time_step_deposition_fraction(daughter_2_id);
	}
	// has one daughter, just send it all through
	else if(daughter_1_id != -1)
	{
		// calculate deposition fraction in daughters
		calculate_airway_time_step_deposition_fraction(daughter_1_id);
	}
	else
	{
		// recursive loop ends here if no daughters
	}


}


// calculate the airway deposition fraction and return them in the arguments
void JamesParticleDeposition::get_airway_deposition_probability(unsigned int element_number, double mean_velocity, double& p_airway_total, double& p_airway_sed, double& p_airway_imp, double& p_airway_dif)
{


	// calculate common quantities

	// get some airway parameters
	double pipe_radius = p_airway_data->at(element_number).get_radius();
	double pipe_diameter = 2.0*pipe_radius;
	double area = M_PI*pipe_radius*pipe_radius;

	// calculate deposition efficiencies


	if(es->parameters.get<bool>("particle_sedimentation"))
	{

		double pipe_angle_from_gravity = gravity_angles[element_number];
		double pipe_length = p_airway_data->at(element_number).get_length();
		// to calc time in branch, we assume that the particle is travelling with the air velocity
		double time_in_branch = pipe_length/mean_velocity;

		// calculate the sedimentation probability
		p_airway_sed = calculate_sedimentation_probability(pipe_angle_from_gravity,pipe_radius,pipe_length,mean_velocity,time_in_branch);

	}

	if(es->parameters.get<bool>("particle_impaction"))
	{

		double branching_angle = branching_angles[element_number];
		// to calc time in branch, we assume that the particle is travelling with the air velocity
		// airway reynolds number
		double reynolds_number = density * pipe_diameter * mean_velocity / viscosity;
		double stokes_number = cunningham_correction_factor*particle_density*pow(particle_diameter,2.0)
											*mean_velocity/(18*viscosity*pipe_diameter);		//yeh1974

		// calculate the sedimentation probability
		p_airway_imp = calculate_impaction_probability(branching_angle,stokes_number,reynolds_number);

	}


	// add up the probabilities
	p_airway_total = p_airway_sed + p_airway_imp + p_airway_dif - p_airway_sed*p_airway_imp - p_airway_sed*p_airway_dif - p_airway_imp*p_airway_dif + p_airway_imp*p_airway_sed*p_airway_dif;


	// calculate the proportions due to each mechanism

	// if all are zero we have a divide by zero
	if(fabs(p_airway_sed) < 1e-10 && fabs(p_airway_imp) < 1e-10 && fabs(p_airway_dif) < 1e-10)
	{
		p_airway_sed = 0.;
		p_airway_imp = 0.;	
		p_airway_dif = 0.;	
	}
	else
	{

		// calculate the proportion of the total added without subtraction
		double sed_prop = p_airway_sed/(p_airway_sed+p_airway_imp+p_airway_dif);
		double imp_prop = p_airway_imp/(p_airway_sed+p_airway_imp+p_airway_dif);
		double dif_prop = p_airway_dif/(p_airway_sed+p_airway_imp+p_airway_dif);

		p_airway_sed = sed_prop * p_airway_total;
		p_airway_imp = imp_prop * p_airway_total;
		p_airway_dif = dif_prop * p_airway_total;
	}





}




double JamesParticleDeposition::calculate_sedimentation_probability(double pipe_angle_from_gravity,double pipe_radius,double pipe_length,double mean_pipe_velocity, double time_in_branch)
{
	if(mean_pipe_velocity > 1e-14)
	{
		return 1 - exp(-4*gravity_magnitude*cunningham_correction_factor*particle_density*pow(particle_radius,2.0)*pipe_length*cos(pipe_angle_from_gravity)
								/(9*M_PI*viscosity*pipe_radius*mean_pipe_velocity));
	}
	else if(mean_pipe_velocity < -1e-14)
	{
		mean_pipe_velocity *= -1;
		return 1 - exp(-4*gravity_magnitude*cunningham_correction_factor*particle_density*pow(particle_radius,2.0)*pipe_length*cos(pipe_angle_from_gravity)
								/(9*M_PI*viscosity*pipe_radius*mean_pipe_velocity));
	}
	else
	{
		double time_in_pipe = time_in_branch * es->parameters.get<double>("time_scale_factor");
		return 1 - exp(-4*gravity_magnitude*cunningham_correction_factor*particle_density*pow(particle_radius,2.0)*cos(pipe_angle_from_gravity) * time_in_pipe
								/(9*M_PI*viscosity*pipe_radius));
	}
}


double JamesParticleDeposition::calculate_brownian_probability(double pipe_radius, double pipe_length,double mean_pipe_velocity, double time_in_branch)
{
	mean_pipe_velocity = fabs(mean_pipe_velocity);
	double x = pipe_length*diffusion_coefficient/(2*pow(pipe_radius,2.0)*mean_pipe_velocity);

	if(mean_pipe_velocity > 1e-14)
		return 1 - 0.819*exp(-7.315*x) - 0.0976*exp(-44.61*x) - 0.0325*exp(-114.0*x) - 0.0509*exp(-79.31*pow(x,2.0/3.0));
	else
	{
		double time_in_pipe = time_in_branch * es->parameters.get<double>("time_scale_factor");
		return 1 - exp(-5.784*diffusion_coefficient*time_in_pipe/pow(pipe_radius,2.0));
	}
}

double JamesParticleDeposition::calculate_impaction_probability(double branching_angle,double stokes_number, double reynolds_number)
{

	double impaction_probability = 0.;

	if(es->parameters.get<unsigned int>("impaction_type") == 1)
	{
		if(stokes_number < 0.04)
		{
			impaction_probability = 0.000654*exp(55.7*pow(stokes_number,0.954)) * pow(reynolds_number,1./3.) * sin(branching_angle);
		}
		else
		{
			impaction_probability = (0.19-0.193*exp(-9.5*pow(stokes_number,1.565))) * pow(reynolds_number,1./3.) * sin(branching_angle);
		}
	}
	else
	{
		if(branching_angle*stokes_number < 1.0)
			impaction_probability = 1. - 2./M_PI*acos(branching_angle*stokes_number) + 1./M_PI*sin(2*acos(branching_angle*stokes_number));
		else
			impaction_probability = 1.;

	}

	// if impaction probability is super high, check we are doing stuff right
	if(impaction_probability > 1 - 1e-10 )
	{
		if(deposition_verbose)
		{		
			std::cout << "impaction_probability = " << impaction_probability << std::endl;
			std::cout << "stokes_number = " << stokes_number << std::endl;
			std::cout << "branching_angle = " << branching_angle << std::endl;
			std::cout << "reynolds_number = " << reynolds_number << std::endl;
		}

		// no probability can be > 1
		impaction_probability = 1.;
		

	}

	return impaction_probability;
}


