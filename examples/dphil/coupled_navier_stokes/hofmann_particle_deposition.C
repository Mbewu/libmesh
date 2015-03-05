
#include "hofmann_particle_deposition.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

HofmannParticleDeposition::HofmannParticleDeposition (EquationSystems& es_in, std::vector<std::vector<double> >& element_data_in): 
		es (&es_in), element_data(element_data_in), total_particles_inhaled(0), end_sim_at_end_time(true)

{

	// for the zeroth which is kinda undefined
	efficiency_function.push_back(std::vector<double> (element_data.size(),0.)); //first one is all zero.
	total_efficiency_function.resize(element_data.size(),0.); //resize total efficiency function
	particle_density = es->parameters.get<double>("particle_density");
	double lambda = es->parameters.get<double>("mean_free_path");
	double particle_diameter =  es->parameters.get<double>("particle_diameter");
	cunningham_correction_factor = 1. + 2.*lambda/particle_diameter*
																	(1.257 + 0.4*exp(-1.1*(particle_diameter/(2*lambda))));

	particle_radius = 0.5 * particle_diameter;
	
	gravity_magnitude = es->parameters.get<double>("gravity");
	viscosity = es->parameters.get<double>("viscosity");

	//from stokes-einstein relation
	double boltzmann_constant = 1.3806488e-23;
	diffusion_coefficient = boltzmann_constant * 310.15 / (3*M_PI * viscosity * particle_diameter);



	// ******************* setup branching angles and angle from gravity ****** //

	branching_angles.resize(element_data.size(),0.);
	gravity_angles.resize(element_data.size(),0.);
  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es->get_mesh();

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();


	// let us set gravity here
	if(es->parameters.get<unsigned int>("gravity_type") == 2)
	{					
		for ( ; el != end_el; ++el)
		{
			const Elem* elem = *el;
			unsigned int subdomain_id = elem->subdomain_id();
			if(subdomain_id > 0)
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
					std::cout << "gravity set" << std::endl;
					break;
				}
	
			}
		}			
	}

	//setup gravity vector
	Point gravity(0.,0.,0.);
	if(es->parameters.get<unsigned int>("gravity_type") == 0)
	{
		if(es->parameters.get<bool>("threed"))
			gravity(2) = -1.;
		else
			gravity(1) = -1.;					
	}
	else if(es->parameters.get<unsigned int>("gravity_type") == 1 || es->parameters.get<unsigned int>("gravity_type") == 2)
	{
		if(es->parameters.get<bool>("threed"))
		{
			std::cout << "hmmm" << std::endl;
			gravity(0) = es->parameters.get<double>("gravity_x");
			gravity(1) = es->parameters.get<double>("gravity_y");
			gravity(2) = es->parameters.get<double>("gravity_z");
			std::cout << "gravity = " << gravity << std::endl;
		}
		else
		{
			std::cout << "hmmm" << std::endl;
			gravity(0) = es->parameters.get<double>("gravity_x");
			gravity(1) = es->parameters.get<double>("gravity_y");
			std::cout << "gravity = " << gravity << std::endl;	
		}
	}
	gravity = gravity.unit();


  el     = mesh.active_local_elements_begin();


				
  for ( ; el != end_el; ++el)
  {
		const Elem* elem = *el;
		unsigned int subdomain_id = elem->subdomain_id();
		if(subdomain_id > 0)
		{
			

			//element data object starts numbering from 0
			//and the values referenced in it also do so need to take this into account

			const int current_el_idx = elem->id();
			Point p0 = *elem->get_node(0);
			Point p1 = *elem->get_node(1);
			Point direction = p1 - p0;
			direction = direction.unit();

			// if we have a parent calculate the branching angle
			if((int)element_data[current_el_idx][0] != -1)
			{
				unsigned int parent_id = (unsigned int)element_data[current_el_idx][0];
				const Elem* parent_elem = mesh.elem(parent_id);
				Point parent_p0 = *parent_elem->get_node(0);
				Point parent_p1 = *parent_elem->get_node(1);
				Point parent_direction = parent_p1 - parent_p0;
				parent_direction = parent_direction.unit();

				double branching_angle = acos(direction*parent_direction);
				branching_angles[current_el_idx] = branching_angle;
				
			}

			if(current_el_idx == 1)
			{
				std::cout << "direction = " << direction << std::endl;
				std::cout << "gravity = " << gravity << std::endl;
			}
			//calculate the gravity angle
			double gravity_angle = acos(direction*gravity);
			gravity_angles[current_el_idx] = gravity_angle;			
		}
	}




}


// calculate the flow rate at different times
void HofmannParticleDeposition::calculate_flow_rate()
{

	double local_time = 0.;
	unsigned int t_step = 0;
	double dt = es->parameters.get<double>("dt");

	flow_function.push_back(std::vector<double> (element_data.size(),0.)); //first one is all zero.
	times.push_back(0);
	// ************* TIME LOOP ******************** //
	// let us do this in dimensional terms
	while (local_time + 1e-10 < es->parameters.get<Real> ("end_time"))
	{

		// INCREMENT TIME
		++t_step;

		if(local_time + dt > es->parameters.get<Real> ("end_time"))
		{
			dt = es->parameters.get<Real> ("end_time") - local_time;
		}
		local_time += dt;

		double initial_flow = time_scaling(local_time) *  es->parameters.get<double>("flow_mag_1d") / 
											(es->parameters.get<double>("velocity_scale") * pow(es->parameters.get<double>("length_scale"),2.0));

		//give flow function a vector for the whole tree
		flow_function.push_back(std::vector<double> (element_data.size(),0.));
		//start at parent - should be element_data 1 - nondimensionalised
		flow_function[t_step][0] =  initial_flow;
		calculate_flow_at_daughters(t_step,0,initial_flow);
		times.push_back(local_time);
		
	}
}


void HofmannParticleDeposition::calculate_flow_at_daughters(unsigned int t_step, unsigned int parent_id, double parent_flow)
{
	int daughter_1_id = (int)element_data[parent_id][1];
	int daughter_2_id = (int)element_data[parent_id][2];
	// has two daughter branches
	if(daughter_1_id != -1 && daughter_2_id != -1)
	{

		//need to calculate the separation of flow based on area.
		double daughter_1_radius = element_data[daughter_1_id][6];
		double daughter_2_radius = element_data[daughter_2_id][6];

		double daughter_1_proportion = pow(daughter_1_radius,2.0)/(pow(daughter_1_radius,2.0) + pow(daughter_2_radius,2.0));
		double daughter_1_flow = daughter_1_proportion * parent_flow;
		flow_function[t_step][daughter_1_id] = daughter_1_flow;
		calculate_flow_at_daughters(t_step,daughter_1_id,daughter_1_flow);

		double daughter_2_proportion = pow(daughter_2_radius,2.0)/(pow(daughter_1_radius,2.0) + pow(daughter_2_radius,2.0));
		double daughter_2_flow = daughter_2_proportion * parent_flow;
		flow_function[t_step][daughter_2_id] = daughter_2_flow;
		calculate_flow_at_daughters(t_step,daughter_2_id,daughter_2_flow);
		
	}
	// has one daughter 1 branch
	else if(daughter_1_id != -1)
	{

		double daughter_1_flow = parent_flow;
		flow_function[t_step][daughter_1_id] = daughter_1_flow;
		calculate_flow_at_daughters(t_step,daughter_1_id,daughter_1_flow);
		
	}
	// has one daughter 1 branch
	else if(daughter_2_id != -1)
	{

		double daughter_2_flow = parent_flow;
		flow_function[t_step][daughter_2_id] = daughter_2_flow;
		calculate_flow_at_daughters(t_step,daughter_2_id,daughter_2_flow);
		
	}
	
}

//calculate the deposition efficiency for a particle entering the trachea at time_in
void HofmannParticleDeposition::calculate_deposition_efficiency(unsigned int t_step_in)
{

	// we only deposit particles on inflow

	double dt = es->parameters.get<double>("dt");
	double local_time = t_step_in * dt;
	double initial_flow = flow_function[t_step_in][0];

	if(initial_flow > 1e-10)
	{
		total_particles_inhaled += 1;
		efficiency_function.push_back(std::vector<double> (element_data.size(),0.)); //first one is all zero.

		//start at parent - should be element_data 0 - nondimensionalised
		efficiency_function[t_step_in][0] +=  calculate_efficiency(0,initial_flow);
		std::cout << "initial_flow = " << initial_flow << std::endl;
		std::cout << "calculate_efficiency(0,initial_flow) = " << calculate_efficiency(0,initial_flow) << std::endl;
		//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
		double parent_area = M_PI * pow(element_data[0][6],2.0);
		double parent_end_time = local_time + element_data[0][5] / (initial_flow/parent_area);
		double current_efficiency = 1 - efficiency_function[t_step_in][0];
		calculate_deposition_efficiency_at_daughters(t_step_in,0,parent_end_time,initial_flow,current_efficiency);
	}
	
	
}


//calculate the deposition efficiency for a particle entering the trachea at time_in
void HofmannParticleDeposition::calculate_deposition_efficiency_at_daughters
		(unsigned int t_step_in,  unsigned int parent_id,double parent_end_time, double parent_end_flow, double parent_end_particle_weight)
{
	int daughter_1_id = (int)element_data[parent_id][1];
	int daughter_2_id = (int)element_data[parent_id][2];
	
	// has two daughter branches
	if(daughter_1_id != -1 && daughter_2_id != -1)
	{

		//need to calculate the separation of flow based on area.
		double daughter_1_radius = element_data[daughter_1_id][6];
		double daughter_2_radius = element_data[daughter_2_id][6];

		double daughter_1_proportion = pow(daughter_1_radius,2.0)/(pow(daughter_1_radius,2.0) + pow(daughter_2_radius,2.0));
		double daughter_1_initial_weight = daughter_1_proportion * parent_end_particle_weight;
		double daughter_1_flow = 0.;
		if(calculate_flow_at_time(daughter_1_flow,parent_end_time, daughter_1_id))
		{

			//start at parent - should be element_data 0 - nondimensionalised
			efficiency_function[t_step_in][daughter_1_id] +=  calculate_efficiency(daughter_1_id,daughter_1_flow) * daughter_1_initial_weight;
			//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
			double daughter_1_area = M_PI * pow(element_data[daughter_1_id][6],2.0);
			double daughter_1_end_time = parent_end_time + element_data[daughter_1_id][5] / (daughter_1_flow/daughter_1_area);
			double current_efficiency_1 = daughter_1_initial_weight - efficiency_function[t_step_in][daughter_1_id];
			calculate_deposition_efficiency_at_daughters(t_step_in,daughter_1_id,daughter_1_end_time,daughter_1_flow,current_efficiency_1);
		}


		double daughter_2_proportion = pow(daughter_2_radius,2.0)/(pow(daughter_1_radius,2.0) + pow(daughter_2_radius,2.0));
		double daughter_2_initial_weight = daughter_2_proportion * parent_end_particle_weight;
		double daughter_2_flow = 0.;
		if(calculate_flow_at_time(daughter_2_flow,parent_end_time, daughter_2_id))
		{

			//start at parent - should be element_data 0 - nondimensionalised
			efficiency_function[t_step_in][daughter_2_id] +=  calculate_efficiency(daughter_2_id,daughter_2_flow) * daughter_2_initial_weight;
			//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
			double daughter_2_area = M_PI * pow(element_data[daughter_2_id][6],2.0);
			double daughter_2_end_time = parent_end_time + element_data[daughter_2_id][5] / (daughter_2_flow/daughter_2_area);
			double current_efficiency_2 = daughter_2_initial_weight - efficiency_function[t_step_in][daughter_2_id];
			calculate_deposition_efficiency_at_daughters(t_step_in,daughter_2_id,daughter_2_end_time,daughter_2_flow,current_efficiency_2);
		}

		
	}
	// has one daughter 1 branch
	else if(daughter_1_id != -1)
	{

		double daughter_1_initial_weight = parent_end_particle_weight;
		double daughter_1_flow = 0.;
		if(calculate_flow_at_time(daughter_1_flow,parent_end_time, daughter_1_id))
		{

			//start at parent - should be element_data 0 - nondimensionalised
			efficiency_function[t_step_in][daughter_1_id] +=  calculate_efficiency(daughter_1_id,daughter_1_flow) * daughter_1_initial_weight;
			//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
			double daughter_1_area = M_PI * pow(element_data[daughter_1_id][6],2.0);
			double daughter_1_end_time = parent_end_time + element_data[daughter_1_id][5] / (daughter_1_flow/daughter_1_area);
			double current_efficiency_1 = daughter_1_initial_weight - efficiency_function[t_step_in][daughter_1_id];
			calculate_deposition_efficiency_at_daughters(t_step_in,daughter_1_id,daughter_1_end_time,daughter_1_flow,current_efficiency_1);
		}


		
	}
	// has one daughter 1 branch
	else if(daughter_2_id != -1)
	{

		double daughter_2_initial_weight = parent_end_particle_weight;
		double daughter_2_flow = 0.;
		if(calculate_flow_at_time(daughter_2_flow,parent_end_time, daughter_2_id))
		{

			//start at parent - should be element_data 0 - nondimensionalised
			efficiency_function[t_step_in][daughter_2_id] +=  calculate_efficiency(daughter_2_id,daughter_2_flow) * daughter_2_initial_weight;
			//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
			double daughter_2_area = M_PI * pow(element_data[daughter_2_id][6],2.0);
			double daughter_2_end_time = parent_end_time + element_data[daughter_2_id][5] / (daughter_2_flow/daughter_2_area);
			double current_efficiency_2 = daughter_2_initial_weight - efficiency_function[t_step_in][daughter_2_id];
			calculate_deposition_efficiency_at_daughters(t_step_in,daughter_2_id,daughter_2_end_time,daughter_2_flow,current_efficiency_2);
		}
		
	}
	// if no daughters we need to go back up tree when velocity is negative again 
	// - here we are in the alveoli essentially, so we "wait" until velocity is negative and then go back up
	else if(daughter_1_id == -1 && daughter_2_id == -1)
	{

/*
		double daughter_2_initial_weight = parent_end_particle_weight;
		double daughter_2_flow = 0.;
		if(calculate_flow_at_time(daughter_2_flow,parent_end_time, daughter_2_id))
		{

			//start at parent - should be element_data 0 - nondimensionalised
			efficiency_function[t_step_in][daughter_2_id] +=  calculate_efficiency(daughter_2_id,daughter_2_flow) * daughter_2_initial_weight;
			//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
			double daughter_2_area = M_PI * pow(element_data[daughter_2_id][6],2.0);
			double daughter_2_end_time = parent_end_time + element_data[daughter_2_id][5] / (daughter_2_flow/daughter_2_area);
			double current_efficiency_2 = daughter_2_initial_weight - efficiency_function[t_step_in][daughter_2_id];
			calculate_deposition_efficiency_at_daughters(t_step_in,daughter_2_id,daughter_2_end_time,daughter_2_flow,current_efficiency_2);
		}
*/

	}
}

double HofmannParticleDeposition::calculate_efficiency(unsigned int branch_id,double branch_flow)
{
	// to calculate the efficiency of a branch we 
	// 1 - calculate the sedimentation efficiency of a whole pipe
	// 2 - calculate the brownian efficiency of a whole pipe
	// 3 - calculate the efficiency in a bend pipe with half the length of the pipe.

	//std::cout << "branch_flow = " << branch_flow << std::endl;

	double length_scale = es->parameters.get<double>("length_scale");
	double velocity_scale = es->parameters.get<double>("velocity_scale");
	double pipe_radius = length_scale * element_data[branch_id][6];
	double pipe_length = length_scale * element_data[branch_id][5];
	double mean_pipe_velocity = velocity_scale * branch_flow / (M_PI * pow(pipe_radius,2.0)); //flow / area
	double pipe_angle_from_gravity = gravity_angles[branch_id];
	double branching_angle = branching_angles[branch_id];
	//need to redimensionalise variables for this (particle_radius is dimensional already)
	double stokes_number = cunningham_correction_factor*particle_density*pow(particle_radius,2.0)
													*mean_pipe_velocity/(9*viscosity*pipe_radius);		//yeh1974

	//all are now dimensionalised
	double P_sed = calculate_sedimentation_efficiency(pipe_angle_from_gravity,pipe_radius,pipe_length,mean_pipe_velocity);
	double P_dif = calculate_brownian_efficiency(pipe_radius,pipe_length,mean_pipe_velocity);
	double P_imp = calculate_impaction_efficiency(branching_angle,stokes_number);

/*
	std::cout << "pipe_angle_from_gravity = " << pipe_angle_from_gravity << std::endl;
	std::cout << "P_sed = " << P_sed << std::endl;
	std::cout << "P_dif = " << P_dif << std::endl;
	std::cout << "P_imp = " << P_imp << std::endl;
*/

	double total_efficiency = 0.;
	if(mean_pipe_velocity < 1e-10)	//no impaction
		total_efficiency = P_sed + P_dif - P_sed*P_dif;
	else
		total_efficiency = P_sed + P_dif + P_imp - P_sed*P_dif - P_sed*P_imp - P_dif*P_imp + P_sed*P_dif*P_imp;

	return total_efficiency;
}

double HofmannParticleDeposition::calculate_sedimentation_efficiency(double pipe_angle_from_gravity,double pipe_radius,double pipe_length,double mean_pipe_velocity)
{
	if(mean_pipe_velocity > 1e-10)
		return 1 - exp(-4*gravity_magnitude*cunningham_correction_factor*particle_density*pow(particle_radius,2.0)*pipe_length*cos(pipe_angle_from_gravity)
								/(9*M_PI*viscosity*pipe_radius*mean_pipe_velocity));
	else if(mean_pipe_velocity < 1e-10)
	{
		mean_pipe_velocity *= -1;
		return 1 - exp(-4*gravity_magnitude*cunningham_correction_factor*particle_density*pow(particle_radius,2.0)*pipe_length*cos(pipe_angle_from_gravity)
								/(9*M_PI*viscosity*pipe_radius*mean_pipe_velocity));
	}
	else
	{
		double time_in_pipe = es->parameters.get<double>("dt") * es->parameters.get<double>("time_scale_factor");
		return 1 - exp(-4*gravity_magnitude*cunningham_correction_factor*particle_density*pow(particle_radius,2.0)*cos(pipe_angle_from_gravity) * time_in_pipe
								/(9*M_PI*viscosity*pipe_radius));
	}
}


double HofmannParticleDeposition::calculate_brownian_efficiency(double pipe_radius, double pipe_length,double mean_pipe_velocity)
{
	mean_pipe_velocity = fabs(mean_pipe_velocity);
	double x = pipe_length*diffusion_coefficient/(2*pow(pipe_radius,2.0)*mean_pipe_velocity);

	if(fabs(mean_pipe_velocity) > 1e-10)
		return 1 - 0.819*exp(-7.315*x) - 0.0976*exp(-44.61*x) - 0.0325*exp(-114.0*x) - 0.0509*exp(-79.31*pow(x,2.0/3.0));
	else
	{
		double time_in_pipe = es->parameters.get<double>("dt") * es->parameters.get<double>("time_scale_factor");
		return 1 - exp(-5.784*diffusion_coefficient*time_in_pipe/pow(pipe_radius,2.0));
	}
}

double HofmannParticleDeposition::calculate_impaction_efficiency(double branching_angle,double stokes_number)
{
	if(branching_angle*stokes_number < 1.0)
		return 1. - 2./M_PI*acos(branching_angle*stokes_number) + 1./M_PI*sin(2*acos(branching_angle*stokes_number));
	else
		return 1.;
}


//calculate the deposition efficiency for a particle entering the trachea at time_in
bool HofmannParticleDeposition::calculate_flow_at_time(double& flow, double time_in, unsigned int branch_no)
{

	for(unsigned int i=0; i<times.size(); i++)
	{
		if(time_in > times[i])
		{
			double delta_t = (time_in - times[i])/(times[i+1] - times[i]);
			double flow_i = flow_function[i][branch_no];
			double flow_i_plus_1 = flow_function[i+1][branch_no];
			flow = flow_i + delta_t * (flow_i_plus_1 - flow_i);
			return true;
		}
	}

	// if we get here it means we are at the end of the time so we can decide whether 
	if(end_sim_at_end_time)
		return false;
	else
	{
		flow = 0.;
		return true;
	}
	
		

	std::cout << "current time = " << time_in << " cannot be found. EXITING..." << std::endl;	
	std::exit(0);
}


double HofmannParticleDeposition::time_scaling(double time_in)
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
		time_flow = sin(2*pi*time_in/period);
	else if(unsteady == 2)
	{
		//quadratic (up to max of one) then constant
		if(time_in < ramp_duration)
		{
			if(quadratic_ramp)
				time_flow = -4.0/pow(ramp_duration*2,2.0)*time_in*(time_in-ramp_duration*2);
			else	//linear
				time_flow = time_in/ramp_duration;
		}
		else
		{
			time_flow = 1.0;
		}
	}
	else if(unsteady == 3)
		time_flow = 1 - cos(2*pi*time_in/period);
	else if(unsteady == 4)
		time_flow = -sin(2*pi*time_in/period);
	else if(unsteady == 5)
		time_flow = cos(2*pi*time_in/period);		
	else if(unsteady == 6)
	{
		//quadratic (up to max of one) then constant
		if(time_in < period)
		{
			time_flow = 1.0;
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

	return time_flow;

}


void HofmannParticleDeposition::calculate_total_efficiency_function()
{
	
	for(unsigned int i=0; i<efficiency_function.size(); i++)
		for(unsigned int j=0; j<total_efficiency_function.size(); j++)
			total_efficiency_function[j] += efficiency_function[i][j];
	
	for(unsigned int i=0; i<total_efficiency_function.size(); i++)
		total_efficiency_function[i] /= total_particles_inhaled;
}
