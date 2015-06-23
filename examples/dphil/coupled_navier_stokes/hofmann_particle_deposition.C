
#include "hofmann_particle_deposition.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

HofmannParticleDeposition::HofmannParticleDeposition (EquationSystems& es_in, std::vector<Airway>& _airway_data, unsigned int _num_generations,std::vector<unsigned int> _subdomains_1d): 
		es (&es_in), airway_data(_airway_data), total_particles_inhaled(0), end_sim_at_end_time(true), num_generations(_num_generations), subdomains_1d(_subdomains_1d)

{

	// for the zeroth which is kinda undefined
	efficiency_function.push_back(std::vector<double> (airway_data.size(),0.)); //first one is all zero.
	tb_efficiency_per_generation.push_back(std::vector<double> (num_generations,0.)); //first one is all zero.
	alveolar_efficiency_per_generation.push_back(std::vector<double> (num_generations,0.)); //first one is all zero.
	non_deposited_particle_weight.push_back(0.); //first one is all zero.
	output_particle_weight.push_back(0.); //first one is all zero.

	total_efficiency_function.resize(airway_data.size(),0.); //resize total efficiency function
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
	diffusion_coefficient = boltzmann_constant * 310.15 * cunningham_correction_factor / (3*M_PI * viscosity * particle_diameter);

	double fluid_density = es->parameters.get<double>("density");
	sedimentation_velocity = 2./9.*(particle_density - fluid_density)*gravity_magnitude*pow(particle_radius,2.0)/viscosity;
	alveolus_radius = 0.5*es->parameters.get<double>("alveolar_diameter");
	dt = es->parameters.get<double>("dt");


	// ******************* setup branching angles and angle from gravity ****** //

	branching_angles.resize(airway_data.size(),0.);
	gravity_angles.resize(airway_data.size(),0.);
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
			if(airway_data[current_el_idx].has_parent())
			{
				unsigned int parent_id = airway_data[current_el_idx].get_parent();
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

	flow_function.push_back(std::vector<double> (airway_data.size(),0.)); //first one is all zero.
	times.push_back(0);
	// ************* TIME LOOP ******************** //
	// let us do this in dimensional terms

	std::cout << "setting up flow crap" << std::endl;
	std::cout << "es->parameters.get<Real> (\"end_time\")" << es->parameters.get<Real> ("end_time") << std::endl;

	while (local_time + 1e-14 < es->parameters.get<Real> ("end_time"))
	{
		std::cout << "local_time = " << local_time << std::endl;
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
		flow_function.push_back(std::vector<double> (airway_data.size(),0.));
		//start at parent - should be element_data 1 - nondimensionalised
		flow_function[t_step][0] =  initial_flow;
		calculate_flow_at_daughters(t_step,0,initial_flow);
		times.push_back(local_time);
		
	}

	std::cout <<"times.size() = " << times.size() << std::endl;

}


// THIS NEEDS TO BE FIXED
void HofmannParticleDeposition::calculate_flow_at_daughters(unsigned int t_step, unsigned int parent_id, double parent_flow)
{
	std::vector<unsigned int> daughter_ids = airway_data[parent_id].get_daughters();
	int daughter_1_id = -1;
	int daughter_2_id = -1;
	if(daughter_ids.size() > 0)
		daughter_1_id = daughter_ids[0];
	if(daughter_ids.size() > 1)
		daughter_1_id = daughter_ids[1];

	// has two daughter branches
	if(daughter_1_id != -1 && daughter_2_id != -1)
	{

		//need to calculate the separation of flow based on area.
		double daughter_1_radius = airway_data[daughter_1_id].get_radius();
		double daughter_2_radius = airway_data[daughter_2_id].get_radius();

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

	std::cout << "hi" << std::endl;
	double local_time = t_step_in * dt;
	double initial_flow = flow_function[t_step_in][0];
	std::cout << "k" << std::endl;

	std::cout <<"times.size() = " << times.size() << std::endl;
	if(initial_flow > 1e-14)
	{
		total_particles_inhaled += 1;
		efficiency_function.push_back(std::vector<double> (airway_data.size(),0.)); //first one is all zero.
		tb_efficiency_per_generation.push_back(std::vector<double> (num_generations,0.)); //first one is all zero.
		alveolar_efficiency_per_generation.push_back(std::vector<double> (num_generations,0.)); //first one is all zero.
		non_deposited_particle_weight.push_back(0.); //first one is all zero.
		output_particle_weight.push_back(0.); //first one is all zero.

	std::cout << "yar = " << initial_flow << std::endl;
		//start at parent - should be element_data 0 - nondimensionalised
		double deposition_efficiency = calculate_efficiency(initial_flow,airway_data[0].get_radius(),airway_data[0].get_length(),gravity_angles[0],branching_angles[0],0.);
	std::cout << "humph" << std::endl;
		efficiency_function[t_step_in][0] +=  deposition_efficiency;
	std::cout << "yas queen" << std::endl;
		tb_efficiency_per_generation[t_step_in][airway_data[0].get_generation()] += deposition_efficiency;
		std::cout << "initial_flow = " << initial_flow << std::endl;
		std::cout << "calculate_efficiency(0,initial_flow) = " << deposition_efficiency << std::endl;
		//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
		double parent_area = M_PI * pow(airway_data[0].get_radius(),2.0);
		double parent_end_time = local_time + airway_data[0].get_length() / (initial_flow/parent_area);
		double current_efficiency = 1 - deposition_efficiency;
		std::cout << "current_efficiency" << current_efficiency << std::endl;
		calculate_deposition_efficiency_in_tb_daughters(t_step_in,0,parent_end_time,initial_flow,current_efficiency);
	}
	
	
}


//calculate the deposition efficiency for a particle entering the trachea at time_in
void HofmannParticleDeposition::calculate_deposition_efficiency_in_tb_daughters
		(unsigned int t_step_in,  unsigned int parent_id,double parent_end_time, double parent_end_flow, double parent_end_weight)
{


	std::vector<unsigned int> daughter_ids = airway_data[parent_id].get_daughters();
	int daughter_1_id = -1;
	int daughter_2_id = -1;
	if(daughter_ids.size() > 0)
		daughter_1_id = daughter_ids[0];
	if(daughter_ids.size() > 1)
		daughter_1_id = daughter_ids[1];

	int generation = (int)airway_data[parent_id].get_generation() + 1;
	double enhancement_factor = 1.;
	if(generation <= 6)
		enhancement_factor = es->parameters.get<double>("max_enhancement_factor") -  1./6. * generation;
	
	//std::cout << "hi" << std::endl;
	//std::cout <<"times.size() = " << times.size() << std::endl;
	//std::cout << "in calculate_deposition_efficiency_in_tb_daughters - gen " << element_data[parent_id][7] + 1 << std::endl;
	//std::cout << "weight = " << parent_end_weight << std::endl;
	//std::cout << "parent_end_time = " << parent_end_time << std::endl;
	//std::cout << "parent_end_flow = " << parent_end_flow << std::endl;
	//need to calculate the separation of flow based on area.
	double daughter_1_radius = 0.;
	double daughter_2_radius = 0.;
	if(daughter_1_id != -1)
		daughter_1_radius = airway_data[daughter_1_id].get_radius();

	if(daughter_2_id != -1)
		daughter_2_radius = airway_data[daughter_2_id].get_radius();

	bool end_time = false;
	//if has daughter 1 branch calculate the deposition in this branch
	if(daughter_1_id != -1)
	{
		// proportion of particles to daughter 1
		double daughter_1_proportion = pow(daughter_1_radius,2.0)/(pow(daughter_1_radius,2.0) + pow(daughter_2_radius,2.0));
		// weight of particles entering daugher 1
		double daughter_1_initial_weight = daughter_1_proportion * parent_end_weight;

		double daughter_1_flow = 0.;
		// if not at the end of time do more deposition
		if(calculate_flow_at_time(daughter_1_flow,parent_end_time, daughter_1_id))
		{

			// if we are in a breath hold phase, add up the time until flow begins again, or the end of time and do static deposition
			if(fabs(daughter_1_flow) < 1e-14)	//breath hold we simulate
			{
				// this is where we turn around...
				double current_time = parent_end_time;
				while(true)
				{
					current_time += dt;
					// should do the function and then use that value
					if(!calculate_flow_at_time(daughter_1_flow,current_time, daughter_1_id))
					{
						end_time = true;
						break;
					}					

					if(fabs(daughter_1_flow) > 1e-14)
						break;


				}
				double time_in_branch = current_time - parent_end_time;

				// static deposition for time_in_branch
				double deposition_efficiency = enhancement_factor * calculate_efficiency(0.,airway_data[daughter_1_id].get_radius(),airway_data[daughter_1_id].get_length(),gravity_angles[daughter_1_id],branching_angles[daughter_1_id],time_in_branch) * daughter_1_initial_weight;
				efficiency_function[t_step_in][daughter_1_id] +=  deposition_efficiency;
				tb_efficiency_per_generation[t_step_in][airway_data[daughter_1_id].get_generation()] += deposition_efficiency;
				//update the time coming from parent
				parent_end_time = current_time;
				// reduce the weight "coming into" the branch
				daughter_1_initial_weight = daughter_1_initial_weight - deposition_efficiency;
			}

			// potentially after a breath hold, if flow is still going down the branch calculate its deposition and 
			if(daughter_1_flow > 1e-14)
			{
				//normal deposition efficiency at this flow rate
				double deposition_efficiency = enhancement_factor * calculate_efficiency(daughter_1_flow,airway_data[daughter_1_id].get_radius(),airway_data[daughter_1_id].get_length(),gravity_angles[daughter_1_id],branching_angles[daughter_1_id],0.) * daughter_1_initial_weight;
				efficiency_function[t_step_in][daughter_1_id] +=  deposition_efficiency;
				tb_efficiency_per_generation[t_step_in][airway_data[daughter_1_id].get_generation()] += deposition_efficiency;		

				// update the time	
				// calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
				double daughter_1_area = M_PI * pow(airway_data[daughter_1_id].get_radius(),2.0);
				double daughter_1_end_time = parent_end_time + airway_data[daughter_1_id].get_length() / (daughter_1_flow/daughter_1_area);

				// update the particle weight leaving the daughter
				double daughter_1_end_weight = daughter_1_initial_weight - deposition_efficiency;
				calculate_deposition_efficiency_in_tb_daughters(t_step_in,daughter_1_id,daughter_1_end_time,daughter_1_flow,daughter_1_end_weight);
			}
			// if flow is going reversed then go to the parent, i.e. the one that called this and do no deposition
			else if(daughter_1_flow < -1e-14)	//breath hold we simulate
			{
				// should always be existing 
				calculate_deposition_efficiency_in_tb_parent(t_step_in,parent_id,parent_end_time,parent_end_flow,daughter_1_initial_weight);
			}
		}
		else
			end_time = true;

		
		// end of time save particles not deposited
		if(end_time)
		{
			non_deposited_particle_weight[t_step_in] += daughter_1_initial_weight; 
		}
	}

	end_time = false;
	// if there is a second daughter
	if(daughter_2_id != -1)
	{
		// proportion of particles to daughter 2
		double daughter_2_proportion = pow(daughter_2_radius,2.0)/(pow(daughter_1_radius,2.0) + pow(daughter_2_radius,2.0));
		// weight of particles entering daughter 2
		double daughter_2_initial_weight = daughter_2_proportion * parent_end_weight;
		double daughter_2_flow = 0.;

		// if end of time then end the simulation
		if(calculate_flow_at_time(daughter_2_flow,parent_end_time, daughter_2_id))
		{

			// breath hold
			if(fabs(daughter_2_flow) < 1e-14)
			{
				double current_time = parent_end_time;
				while(true)
				{
					current_time += dt;
					// if end of time or flow nonzero again
					if(!calculate_flow_at_time(daughter_2_flow,current_time,daughter_2_id))
					{
						end_time = true;
						break;
					}
					if(fabs(daughter_2_flow) > 1e-14)
						break;

				}
				double time_in_branch = current_time - parent_end_time;

				// static deposition in daughter 2
				double deposition_efficiency = enhancement_factor * calculate_efficiency(0.,airway_data[daughter_2_id].get_radius(),airway_data[daughter_2_id].get_length(),gravity_angles[daughter_2_id],branching_angles[daughter_2_id],time_in_branch) * daughter_2_initial_weight;
				efficiency_function[t_step_in][daughter_2_id] +=  deposition_efficiency;
				tb_efficiency_per_generation[t_step_in][airway_data[daughter_2_id].get_generation()] += deposition_efficiency;

				// update the time at which particles enter from the parent
				parent_end_time = current_time;
				daughter_2_initial_weight = daughter_2_initial_weight - deposition_efficiency;

			}


			// if still going down
			if(daughter_2_flow > 1e-14)
			{
				// normal deposition
				double deposition_efficiency = enhancement_factor * calculate_efficiency(daughter_2_flow,airway_data[daughter_2_id].get_radius(),airway_data[daughter_2_id].get_length(),gravity_angles[daughter_2_id],branching_angles[daughter_2_id],0.) * daughter_2_initial_weight;
				efficiency_function[t_step_in][daughter_2_id] +=  deposition_efficiency;
				tb_efficiency_per_generation[t_step_in][airway_data[daughter_2_id].get_generation()] += deposition_efficiency;

				// update the time
				//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
				double daughter_2_area = M_PI * pow(airway_data[daughter_2_id].get_radius(),2.0);
				double daughter_2_end_time = parent_end_time + airway_data[daughter_2_id].get_length() / (daughter_2_flow/daughter_2_area);
				double daughter_2_end_weight = daughter_2_initial_weight - deposition_efficiency;
				calculate_deposition_efficiency_in_tb_daughters(t_step_in,daughter_2_id,daughter_2_end_time,daughter_2_flow,daughter_2_end_weight);
			}
			else if(daughter_2_flow < -1e-14)	//breath hold we simulate
			{
				// should always be existing
				calculate_deposition_efficiency_in_tb_parent(t_step_in,parent_id,parent_end_time,parent_end_flow,daughter_2_initial_weight);
			}
		}
		else
			end_time = true;

		// end of time save particles not deposited
		if(end_time)
		{
			non_deposited_particle_weight[t_step_in] += daughter_2_initial_weight; 
		}

		
	}


	// if no daughters we are in the alveolar region
	// each time, we reconstruct the alveolar geometry, but it is always the same
	if(daughter_1_id == -1 && daughter_2_id == -1)
	{
		if(es->parameters.get<unsigned int>("alveolated_1d_tree"))
		{
			//generate temporary alveolar tree from
			unsigned int num_alveolar_generations = airway_data[parent_id].get_num_alveolar_generations();
			//this should construct a tree from the last parent branch
			construct_temporary_alveolated_1d_tree (parent_id);
			// there is no deposition in the first branch as this is part of the tracheobronchial tree
			// gets all the weight and the flow and knows the time
			calculate_deposition_efficiency_in_alveolar_daughters(t_step_in,0,parent_end_time,parent_end_flow,parent_end_weight);
		}
		else
		{
			std::cout << "ONLY DOING ALVEOLATED TREES FOR NOW, GOOD BYE." << std::endl;
			std::exit(0);
		}

	}
}


//calculate the deposition efficiency for a particle entering the trachea at time_in
void HofmannParticleDeposition::calculate_deposition_efficiency_in_alveolar_daughters
		(unsigned int t_step_in,  unsigned int parent_alveolar_id,double parent_end_time, double parent_end_flow, double parent_end_weight)
{

	//std::cout << "parent_end_time = " << parent_end_time << std::endl;
	//std::cout <<"times.size() = " << times.size() << std::endl;
	//std::cout << "in calculate_deposition_efficiency_in_alv_daughters - gen " << temporary_alveolar_tree[parent_alveolar_id][7] + 1 << std::endl;
	//std::cout << "weight = " << parent_end_weight << std::endl;

	int daughter_1_id = (int)temporary_alveolar_tree[parent_alveolar_id][2];
	int daughter_2_id = (int)temporary_alveolar_tree[parent_alveolar_id][3];

	bool end_time = false;
	
	// will always have two daughters
	if(daughter_1_id != -1 && daughter_2_id != -1)
	{

		//need to calculate the separation of flow based on area.
		double daughter_1_radius = temporary_alveolar_tree[daughter_1_id][6];
		double daughter_2_radius = temporary_alveolar_tree[daughter_2_id][6];

		// get the proportion of flow entering
		double daughter_1_proportion = temporary_alveolar_tree[daughter_1_id][4];	// for the flow
		double daughter_1_initial_weight = 0.5 * parent_end_weight;	//each daughter gets equal weight cause symmetry
		double daughter_1_flow = 0.;

		//std::cout << "daughter_1_initial_weight = " << daughter_1_initial_weight << std::endl;

		// if not at the end of time do some deposition
		if(calculate_flow_at_time(daughter_1_flow,parent_end_time, temporary_alveolar_tree[daughter_1_id][0]))
		{
			// flow is proportion of flow at the head of the acinar airways
			daughter_1_flow *= daughter_1_proportion;
			unsigned int alveolar_generation = (unsigned int)temporary_alveolar_tree[daughter_1_id][7];
			unsigned int num_alveolar_generations = airway_data[temporary_alveolar_tree[parent_alveolar_id][0]].get_num_alveolar_generations();

			// if breath hold then do some extra deposition, but leave the particles that were in the alveoli
			// assume that they are properly deposited
			if(fabs(daughter_1_flow) < 1e-14)
			{

				// this is where we turn around...
				double current_time = parent_end_time;
				while(true)
				{
					current_time += dt;
					// goes up during inspiration, upon exhalation should go through negative again
					if(!calculate_flow_at_time(daughter_1_flow,current_time, temporary_alveolar_tree[daughter_1_id][0]))
					{
						end_time = true;
						break;
					}				

					if(fabs(daughter_1_flow) > 1e-14)
						break;

					

				}
				double time_in_branch = current_time - parent_end_time;


				//start at parent - should be element_data 0 - nondimensionalised
				double deposition_efficiency = calculate_efficiency(0.,temporary_alveolar_tree[daughter_1_id][6],temporary_alveolar_tree[daughter_1_id][5],temporary_alveolar_tree[daughter_1_id][8],temporary_alveolar_tree[daughter_1_id][9],time_in_branch) * daughter_1_initial_weight;

				unsigned int tb_generation = airway_data[temporary_alveolar_tree[parent_alveolar_id][0]].get_generation();

				//add acinar deposition to last airway, all is deposited, none into alveoli
				efficiency_function[t_step_in][temporary_alveolar_tree[parent_alveolar_id][0]] +=  deposition_efficiency;	


				// add deposited part to alveolar
				alveolar_efficiency_per_generation[t_step_in][tb_generation + alveolar_generation] += deposition_efficiency;

				//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
				double daughter_1_area = M_PI * pow(temporary_alveolar_tree[daughter_1_id][6],2.0);
				parent_end_time = current_time;
				daughter_1_initial_weight = daughter_1_initial_weight - deposition_efficiency;

			}

			// if flow, do some deposition
			if(daughter_1_flow > 1e-14)
			{
				// check if we are at the end of the tree
				//std::cout << "alveolar_generation = " << alveolar_generation << std::endl;
				//std::cout << "num_alveolar_generations = " << num_alveolar_generations << std::endl;
				if(alveolar_generation != num_alveolar_generations)
				{
					//start at parent - should be element_data 0 - nondimensionalised
					double deposition_efficiency = calculate_efficiency(daughter_1_flow,temporary_alveolar_tree[daughter_1_id][6],temporary_alveolar_tree[daughter_1_id][5],temporary_alveolar_tree[daughter_1_id][8],temporary_alveolar_tree[daughter_1_id][9],0.) * daughter_1_initial_weight;

					//now calculate probability in acinus, generation/total generations
					unsigned int tb_generation = airway_data[temporary_alveolar_tree[parent_alveolar_id][0]].get_generation();
					double probability_in_alveoli = alveolar_generation/num_alveolar_generations;
					double deposition_in_alveoli = probability_in_alveoli * deposition_efficiency;

					// a proportion goes into the alveoli until exhalation
					temporary_alveolar_tree[daughter_1_id][13] += deposition_in_alveoli;
					if(temporary_alveolar_tree[daughter_1_id][14] < 1e-14)
						temporary_alveolar_tree[daughter_1_id][14] = parent_end_time;	// use the latest time


					//add acinar deposition to last airway, the bit that was not put in the alveoli
					efficiency_function[t_step_in][temporary_alveolar_tree[parent_alveolar_id][0]] +=  deposition_efficiency * (1-probability_in_alveoli);	

					// add deposited part to alveolar
					alveolar_efficiency_per_generation[t_step_in][tb_generation + alveolar_generation] += deposition_efficiency * (1-probability_in_alveoli);

					//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
					double daughter_1_area = M_PI * pow(temporary_alveolar_tree[daughter_1_id][6],2.0);
					double daughter_1_end_time = parent_end_time + temporary_alveolar_tree[daughter_1_id][5] / (daughter_1_flow/daughter_1_area);
					double daughter_1_end_weight = daughter_1_initial_weight - deposition_efficiency;
					calculate_deposition_efficiency_in_alveolar_daughters(t_step_in,daughter_1_id,daughter_1_end_time,daughter_1_flow,daughter_1_end_weight);
				}
				else	// if at end of tree add all particles to alveolus and iterate forward in time until same volume inhaled since start has been exhaled
				{

					//now calculate probability in acinus, generation/total generations
					unsigned int tb_generation = airway_data[temporary_alveolar_tree[parent_alveolar_id][0]].get_generation();

					// a proportion goes into the alveoli until exhalation
					temporary_alveolar_tree[daughter_1_id][13] = daughter_1_initial_weight;
					temporary_alveolar_tree[daughter_1_id][14] = parent_end_time;

					// this is where we turn around...
					double volume_inhaled = 0.;
					double current_time = parent_end_time;
					while(true)
					{
						current_time += dt;
						//std::cout << "current_time = " << current_time << std::endl;
						if(!calculate_flow_at_time(daughter_1_flow,current_time, temporary_alveolar_tree[daughter_1_id][0]))
						{
							end_time = true;
							//std::cout << "hey yah?" << std::endl;
							//std::cout << "value = " << calculate_flow_at_time(daughter_1_flow,current_time, temporary_alveolar_tree[daughter_1_id][0]) << std::endl;
							break;
						}

						volume_inhaled += daughter_1_flow * dt;
						//std::cout << "daughter_1_flow = " << daughter_1_flow << std::endl;
						//std::cout << "volume_inhaled = " << volume_inhaled << std::endl;
						// goes up during inspiration, upon exhalation should go through negative again
						if(volume_inhaled < 0)
							break;


					}

					double time_spent_in_alveoli = current_time - parent_end_time;
					double weight_in_alveoli = daughter_1_initial_weight;
					double size_of_alveoli = 1.0;

					double proportion_deposited_in_alveoli = calculate_deposition_in_alveoli(time_spent_in_alveoli);
					double deposited_in_alveoli = weight_in_alveoli* proportion_deposited_in_alveoli;

					efficiency_function[t_step_in][temporary_alveolar_tree[parent_alveolar_id][0]] +=  deposited_in_alveoli;	//add acinar deposition to last airway

					temporary_alveolar_tree[daughter_1_id][13] = 0.;
					temporary_alveolar_tree[daughter_1_id][14] = 0.;

					// add deposited part to alveolar
					alveolar_efficiency_per_generation[t_step_in][tb_generation + alveolar_generation] += deposited_in_alveoli;

					double daughter_1_end_weight = daughter_1_initial_weight - deposited_in_alveoli;

					//std::cout << "daughter_1_end_weight = " << daughter_1_end_weight << std::endl;

					calculate_deposition_efficiency_in_alveolar_parent(t_step_in,parent_alveolar_id,current_time,daughter_1_flow,daughter_1_end_weight);
				}
			}
			else if(daughter_1_flow < -1e-14)
			{
				calculate_deposition_efficiency_in_alveolar_parent(t_step_in,parent_alveolar_id,parent_end_time,daughter_1_flow,daughter_1_initial_weight);
			}
		}
		else
			end_time = true;

		if(end_time)
		{
			//if we are at the end of time then we need to calculate the amount deposited in the alveoli up until now
			// we are going down so there should be no particles stuck in the alveolus so can just go back up and calculate on the way up
			double daughter_1_end_weight = daughter_1_initial_weight;

			calculate_deposition_efficiency_in_alveolar_parent(t_step_in,parent_alveolar_id,parent_end_time,daughter_1_flow,daughter_1_end_weight);
			
		}


		end_time = false;
		double daughter_2_proportion = temporary_alveolar_tree[daughter_2_id][4];	//for the flow only
		double daughter_2_initial_weight = 0.5 * parent_end_weight;	// each daughter gets equal weight
		double daughter_2_flow = 0.;

		//std::cout << "daughter_2_initial_weight = " << daughter_2_initial_weight << std::endl;
		if(calculate_flow_at_time(daughter_2_flow,parent_end_time, temporary_alveolar_tree[daughter_2_id][0]))
		{

			// flow is proportion of flow at the head of the acinar airways
			daughter_2_flow *= daughter_2_proportion;
			unsigned int alveolar_generation = (unsigned int)temporary_alveolar_tree[daughter_2_id][7];
			unsigned int num_alveolar_generations = airway_data[temporary_alveolar_tree[parent_alveolar_id][0]].get_num_alveolar_generations();


			if(fabs(daughter_2_flow) < 1e-14)
			{

				// this is where we turn around...
				double current_time = parent_end_time;
				while(true)
				{
					current_time += dt;
					// goes up during inspiration, upon exhalation should go through negative again
					if(!calculate_flow_at_time(daughter_2_flow,current_time, temporary_alveolar_tree[daughter_2_id][0]))
					{
						end_time = true;
						break;
					}

					if(fabs(daughter_2_flow) > 1e-14)
						break;

					

				}
				double time_in_branch = current_time - parent_end_time;


				//start at parent - should be element_data 0 - nondimensionalised
				double deposition_efficiency = calculate_efficiency(0.,temporary_alveolar_tree[daughter_2_id][6],temporary_alveolar_tree[daughter_2_id][5],temporary_alveolar_tree[daughter_2_id][8],temporary_alveolar_tree[daughter_2_id][9],time_in_branch) * daughter_2_initial_weight;

				//now calculate probability in acinus, generation/total generations
				unsigned int tb_generation = airway_data[temporary_alveolar_tree[parent_alveolar_id][0]].get_generation();

				//add acinar deposition to last airway
				efficiency_function[t_step_in][temporary_alveolar_tree[parent_alveolar_id][0]] +=  deposition_efficiency;	

				// add deposited part to alveolar
				alveolar_efficiency_per_generation[t_step_in][tb_generation + alveolar_generation] += deposition_efficiency;

				//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
				double daughter_2_area = M_PI * pow(temporary_alveolar_tree[daughter_2_id][6],2.0);
				parent_end_time = current_time;
				daughter_2_initial_weight = daughter_2_initial_weight - deposition_efficiency;

			}

			if(daughter_2_flow > 1e-14)
			{
				// check if we are at the end of the tree
				if(alveolar_generation != num_alveolar_generations)
				{
					//start at parent - should be element_data 0 - nondimensionalised
					double deposition_efficiency = calculate_efficiency(daughter_2_flow,temporary_alveolar_tree[daughter_2_id][6],temporary_alveolar_tree[daughter_2_id][5],temporary_alveolar_tree[daughter_2_id][8],temporary_alveolar_tree[daughter_2_id][9],0.) * daughter_2_initial_weight;

					//now calculate probability in acinus, generation/total generations
					unsigned int tb_generation = airway_data[temporary_alveolar_tree[parent_alveolar_id][0]].get_generation();
					double probability_in_alveoli = alveolar_generation/num_alveolar_generations;
					double deposition_in_alveoli = probability_in_alveoli * deposition_efficiency;

					// a proportion goes into the alveoli until exhalation
					temporary_alveolar_tree[daughter_2_id][13] = deposition_in_alveoli;
					temporary_alveolar_tree[daughter_2_id][14] = parent_end_time;

					//add acinar deposition to last airway
					efficiency_function[t_step_in][temporary_alveolar_tree[parent_alveolar_id][0]] +=  deposition_efficiency * (1-probability_in_alveoli);	

					// add deposited part to alveolar
					alveolar_efficiency_per_generation[t_step_in][tb_generation + alveolar_generation] += deposition_efficiency * (1-probability_in_alveoli);

					//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
					double daughter_2_area = M_PI * pow(temporary_alveolar_tree[daughter_2_id][6],2.0);
					double daughter_2_end_time = parent_end_time + temporary_alveolar_tree[daughter_2_id][5] / (daughter_2_flow/daughter_2_area);
					double daughter_2_end_weight = daughter_2_initial_weight - deposition_efficiency;
					calculate_deposition_efficiency_in_alveolar_daughters(t_step_in,daughter_2_id,daughter_2_end_time,daughter_2_flow,daughter_2_end_weight);
				}
				else	// if at end of tree add all particles to alveolus and iterate forward in time until same volume inhaled since start has been exhaled
				{

					//now calculate probability in acinus, generation/total generations
					unsigned int tb_generation = airway_data[temporary_alveolar_tree[parent_alveolar_id][0]].get_generation();

					// a proportion goes into the alveoli until exhalation
					temporary_alveolar_tree[daughter_2_id][13] = daughter_2_initial_weight;
					temporary_alveolar_tree[daughter_2_id][14] = parent_end_time;

					// this is where we turn around...
					double volume_inhaled = 0.;
					double current_time = parent_end_time;
					while(true)
					{
						current_time += dt;
						if(!calculate_flow_at_time(daughter_2_flow,current_time, temporary_alveolar_tree[daughter_2_id][0]))
						{
							end_time = true;
							break;
						}

						volume_inhaled += daughter_2_flow * dt;
						//std::cout << "daughter_2_flow = " << daughter_2_flow << std::endl;
						//std::cout << "volume_inhaled = " << volume_inhaled << std::endl;
						// goes up during inspiration, upon exhalation should go through negative again
						if(volume_inhaled < 0)
							break;


					}

					double time_spent_in_alveoli = current_time - parent_end_time;
					double weight_in_alveoli = daughter_2_initial_weight;
					double size_of_alveoli = 1.0;

					double proportion_deposited_in_alveoli = calculate_deposition_in_alveoli(time_spent_in_alveoli);
					double deposited_in_alveoli = weight_in_alveoli * proportion_deposited_in_alveoli;

					temporary_alveolar_tree[daughter_2_id][13] = 0.;
					temporary_alveolar_tree[daughter_2_id][14] = 0.;

					efficiency_function[t_step_in][temporary_alveolar_tree[parent_alveolar_id][0]] +=  deposited_in_alveoli;	//add acinar deposition to last airway

					// add deposited part to alveolar
					alveolar_efficiency_per_generation[t_step_in][tb_generation + alveolar_generation] += deposited_in_alveoli;

					double daughter_2_end_weight = daughter_2_initial_weight - deposited_in_alveoli;

					//std::cout << "daughter_2_end_weight = " << daughter_2_end_weight << std::endl;

					calculate_deposition_efficiency_in_alveolar_parent(t_step_in,parent_alveolar_id,current_time,daughter_2_flow,daughter_2_end_weight);
				}
			}
			else if(daughter_2_flow < -1e-14)
			{
				calculate_deposition_efficiency_in_alveolar_parent(t_step_in,parent_alveolar_id,parent_end_time,daughter_2_flow,daughter_2_initial_weight);
			}
		}
		else
			end_time = true;

		if(end_time)
		{
			//if we are at the end of time then we need to calculate the amount deposited in the alveoli up until now
			// we are going down so there should be no particles stuck in the alveolus so can just go back up and calculate on the way up
			double daughter_2_end_weight = daughter_2_initial_weight;

			calculate_deposition_efficiency_in_alveolar_parent(t_step_in,parent_alveolar_id,parent_end_time,daughter_2_flow,daughter_2_end_weight);
			
		}
		
	}
}


//calculate the deposition efficiency for a particle entering the trachea at time_in
void HofmannParticleDeposition::calculate_deposition_efficiency_in_alveolar_parent
		(unsigned int t_step_in,  unsigned int alveolar_id,double daughter_end_time, double daughter_end_flow, double daughter_end_weight)
{

	//std::cout << "daughter_end_time = " << daughter_end_time << std::endl;
	//std::cout <<"times.size() = " << times.size() << std::endl;
	//std::cout << "in calculate_deposition_efficiency_in_alv_parent - gen " << temporary_alveolar_tree[alveolar_id][7] << std::endl;
	//std::cout << "weight = " << daughter_end_weight << std::endl;

	int parent_id = (int)temporary_alveolar_tree[alveolar_id][1];
	int super_parent_id = (int)temporary_alveolar_tree[alveolar_id][0];
	double branch_flow = 0.;
	unsigned int alveolar_generation = (unsigned int)temporary_alveolar_tree[alveolar_id][7];
	unsigned int num_alveolar_generations = airway_data[temporary_alveolar_tree[alveolar_id][0]].get_num_alveolar_generations();
	
	bool end_time = false;
	// check if we are at the end of the tree
	if(alveolar_generation != 0)
	{
		if(calculate_flow_at_time(branch_flow,daughter_end_time, temporary_alveolar_tree[alveolar_id][0]))
		{
			double flow_proportion = temporary_alveolar_tree[alveolar_id][4];
			// flow is proportion of flow at the head of the acinar airways
			branch_flow *= flow_proportion;

			// in a breath hold don't empty the alveoli
			if(fabs(branch_flow) < 1e-14)
			{

				// this is where we turn around...
				double current_time = daughter_end_time;
				while(true)
				{
					current_time += dt;
					// goes up during inspiration, upon exhalation should go through negative again
					if(!calculate_flow_at_time(branch_flow,current_time, temporary_alveolar_tree[alveolar_id][0]))
					{
						end_time = true;
						break;
					}
				
					if(fabs(branch_flow) > 1e-14)
						break;

				

				}
				double time_in_branch = current_time - daughter_end_time;

				//start at parent - should be element_data 0 - nondimensionalised
				double deposition_efficiency = calculate_efficiency(branch_flow,temporary_alveolar_tree[alveolar_id][6],temporary_alveolar_tree[alveolar_id][5],temporary_alveolar_tree[alveolar_id][8],temporary_alveolar_tree[alveolar_id][9],time_in_branch) * daughter_end_weight;

				//now calculate probability in acinus, generation/total generations
				unsigned int tb_generation = airway_data[temporary_alveolar_tree[alveolar_id][0]].get_generation();

				efficiency_function[t_step_in][temporary_alveolar_tree[alveolar_id][0]] +=  deposition_efficiency;	//add acinar deposition to last airway

				// add deposited part to alveolar
				alveolar_efficiency_per_generation[t_step_in][tb_generation + alveolar_generation] += deposition_efficiency;

				daughter_end_time = current_time;
				daughter_end_weight = daughter_end_weight - deposition_efficiency;
			}
			
			if(branch_flow < -1e-14)
			{
				//start at parent - should be element_data 0 - nondimensionalised
				double deposition_efficiency = calculate_efficiency(branch_flow,temporary_alveolar_tree[alveolar_id][6],temporary_alveolar_tree[alveolar_id][5],temporary_alveolar_tree[alveolar_id][8],temporary_alveolar_tree[alveolar_id][9],0.) * daughter_end_weight;

				//now calculate probability in acinus, generation/total generations
				unsigned int tb_generation = airway_data[temporary_alveolar_tree[alveolar_id][0]].get_generation();
				double probability_in_alveoli = alveolar_generation/num_alveolar_generations;
				double deposition_in_alveoli = probability_in_alveoli * deposition_efficiency;

				// a proportion comes from the alveoli
				double weight_in_alveoli = temporary_alveolar_tree[alveolar_id][13];
				temporary_alveolar_tree[alveolar_id][13] = 0.;
				double time_spent_in_alveoli = daughter_end_time - temporary_alveolar_tree[alveolar_id][14];
				temporary_alveolar_tree[alveolar_id][14] = 0.;
				double size_of_alveoli = 1.0;

				double proportion_deposited_in_alveoli = calculate_deposition_in_alveoli(time_spent_in_alveoli);
				double deposited_in_alveoli = weight_in_alveoli * proportion_deposited_in_alveoli;
				double released_from_alveoli = weight_in_alveoli * (1 - proportion_deposited_in_alveoli);


				efficiency_function[t_step_in][temporary_alveolar_tree[alveolar_id][0]] +=  deposition_efficiency + deposited_in_alveoli;	//add acinar deposition to last airway

				// add deposited part to alveolar
				alveolar_efficiency_per_generation[t_step_in][tb_generation + alveolar_generation] += deposition_efficiency + deposited_in_alveoli;

				//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
				double branch_area = M_PI * pow(temporary_alveolar_tree[alveolar_id][6],2.0);
				double parent_end_time = daughter_end_time + temporary_alveolar_tree[alveolar_id][5] / (fabs(branch_flow)/branch_area);
				double parent_end_weight  = daughter_end_weight - deposition_efficiency + released_from_alveoli;
				calculate_deposition_efficiency_in_alveolar_parent(t_step_in,parent_id,parent_end_time,branch_flow,parent_end_weight);
			}
			else if(branch_flow > 1e-14)
			{
				calculate_deposition_efficiency_in_alveolar_daughters(t_step_in,parent_id,daughter_end_time,branch_flow,daughter_end_weight);
			}
		}
		else
			end_time = true;

		// if we are past the end of time then we just want to deposit the alveoli and calculate the remaining particles
		if(end_time)	
		{

			//now calculate probability in acinus, generation/total generations
			unsigned int tb_generation = airway_data[temporary_alveolar_tree[alveolar_id][0]].get_generation();

			// a proportion comes from the alveoli
			double weight_in_alveoli = temporary_alveolar_tree[alveolar_id][13];
			temporary_alveolar_tree[alveolar_id][13] = 0.;
			double time_spent_in_alveoli = daughter_end_time - temporary_alveolar_tree[alveolar_id][14];
			temporary_alveolar_tree[alveolar_id][14] = 0.;
			double size_of_alveoli = 1.0;

			double proportion_deposited_in_alveoli = calculate_deposition_in_alveoli(time_spent_in_alveoli);
			double deposited_in_alveoli = weight_in_alveoli * proportion_deposited_in_alveoli;
			double released_from_alveoli = weight_in_alveoli * (1 - proportion_deposited_in_alveoli);


			efficiency_function[t_step_in][temporary_alveolar_tree[alveolar_id][0]] += deposited_in_alveoli;	//add acinar deposition to last airway

			// add deposited part to alveolar
			alveolar_efficiency_per_generation[t_step_in][tb_generation + alveolar_generation] += deposited_in_alveoli;

			//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
			double branch_area = M_PI * pow(temporary_alveolar_tree[alveolar_id][6],2.0);
			double parent_end_time = daughter_end_time + temporary_alveolar_tree[alveolar_id][5] / (fabs(branch_flow)/branch_area);
			double parent_end_weight  = daughter_end_weight + released_from_alveoli;	// only change is release from alveoli
			calculate_deposition_efficiency_in_alveolar_parent(t_step_in,parent_id,parent_end_time,branch_flow,parent_end_weight);
		}
	}
	else	// if at top of tree go straight to the tracheo bronchial tree
	{
		//std::cout << "okay, we is going back to the parent" << std::endl;
		calculate_deposition_efficiency_in_tb_parent(t_step_in,super_parent_id,daughter_end_time,branch_flow,daughter_end_weight);

	}

}



//calculate the deposition efficiency for a particle entering the trachea at time_in
void HofmannParticleDeposition::calculate_deposition_efficiency_in_tb_parent
		(unsigned int t_step_in,  unsigned int branch_id,double daughter_end_time, double daughter_end_flow, double daughter_end_weight)
{

	//std::cout << "in calculate_deposition_efficiency_in_tb_parent - gen " << element_data[branch_id][7] << std::endl;
	//std::cout << "weight = " << daughter_end_weight << std::endl;
	//std::cout << "daughter_end_time = " << daughter_end_time << std::endl;
	//std::cout << "daughter_end_flow = " << daughter_end_flow << std::endl;
	//std::cout << "parent_id = " << (int)element_data[branch_id][0] << std::endl;
	int parent_id = airway_data[branch_id].get_parent();
	
	//need to calculate the separation of flow based on area.
	double branch_radius = airway_data[branch_id].get_radius();

	double branch_flow = 0.;

	int generation = airway_data[branch_id].get_generation();
	double enhancement_factor = 1.;
	if(generation <= 6)
		enhancement_factor = es->parameters.get<double>("max_enhancement_factor") -  1./6. * generation;
	
	bool end_time = false;
	if(calculate_flow_at_time(branch_flow,daughter_end_time, branch_id))
	{

		
		if(fabs(branch_flow) < 1e-14)
		{

			// this is where we turn around...
			double current_time = daughter_end_time;
			while(true)
			{
				current_time += dt;
				// goes up during inspiration, upon exhalation should go through negative again
				if(!calculate_flow_at_time(branch_flow,current_time, branch_id))
				{
					end_time = true;
					break;
				}

				if(fabs(branch_flow) > 1e-14)
					break;


			}
			double time_in_branch = current_time - daughter_end_time;

			//start at parent - should be element_data 0 - nondimensionalised
			double deposition_efficiency = enhancement_factor * calculate_efficiency(branch_flow,airway_data[branch_id].get_radius(),airway_data[branch_id].get_length(),gravity_angles[branch_id],branching_angles[branch_id],time_in_branch) * daughter_end_weight;
			efficiency_function[t_step_in][branch_id] +=  deposition_efficiency;
			tb_efficiency_per_generation[t_step_in][airway_data[branch_id].get_generation()] += deposition_efficiency;	
			daughter_end_time = current_time;
			daughter_end_weight = daughter_end_weight - deposition_efficiency;
		}

		if(branch_flow < -1e-14)
		{
			//std::cout << "last bit" << std::endl;
			//start at parent - should be element_data 0 - nondimensionalised
			double deposition_efficiency = enhancement_factor * calculate_efficiency(branch_flow,airway_data[branch_id].get_radius(),airway_data[branch_id].get_length(),gravity_angles[branch_id],branching_angles[branch_id],0.) * daughter_end_weight;
			efficiency_function[t_step_in][branch_id] +=  deposition_efficiency;
			tb_efficiency_per_generation[t_step_in][airway_data[branch_id].get_generation()] += deposition_efficiency;			
			//calculate time particle reaches end of tube. distance / velocity =  length / (flow_rate/area)
			double branch_area = M_PI * pow(airway_data[branch_id].get_radius(),2.0);
			double branch_end_time = daughter_end_time + airway_data[branch_id].get_length() / (branch_flow/branch_area);
			double parent_end_weight = daughter_end_weight - deposition_efficiency;

			// has a parent
			if(parent_id != -1)
				calculate_deposition_efficiency_in_tb_parent(t_step_in,parent_id,branch_end_time,branch_flow,parent_end_weight);
			else	// reached the trachea and simply cancel sim
			{
				//std::cout << "end trajectory-" << parent_end_weight << std::endl;
				output_particle_weight[t_step_in] += parent_end_weight;
			}

		}
		else if(branch_flow > 1e-14)
		{
			calculate_deposition_efficiency_in_tb_daughters(t_step_in,branch_id,daughter_end_time,branch_flow,daughter_end_weight);
		}
	}
	else
		end_time = true;

	//end time then just save the amount of particles left over
	if(end_time)	
	{
		non_deposited_particle_weight[t_step_in] += daughter_end_weight; 
	}

}



double HofmannParticleDeposition::calculate_efficiency(double branch_flow,double pipe_radius, double pipe_length, double gravity_angle,double branching_angle,double time_in_branch)
{
	// to calculate the efficiency of a branch we 
	// 1 - calculate the sedimentation efficiency of a whole pipe
	// 2 - calculate the brownian efficiency of a whole pipe
	// 3 - calculate the efficiency in a bend pipe with half the length of the pipe.

	//std::cout << "branch_flow = " << branch_flow << std::endl;

	double length_scale = es->parameters.get<double>("length_scale");
	double velocity_scale = es->parameters.get<double>("velocity_scale");
	pipe_radius *= length_scale;
	pipe_length *= length_scale;
	double mean_pipe_velocity = velocity_scale * branch_flow / (M_PI * pow(pipe_radius,2.0)); //flow / area
	//need to redimensionalise variables for this (particle_radius is dimensional already)
	double stokes_number = cunningham_correction_factor*particle_density*pow(particle_radius,2.0)
													*mean_pipe_velocity/(9*viscosity*pipe_radius);		//yeh1974
	double reynolds_number = particle_density * 2.0*particle_radius * mean_pipe_velocity / viscosity;

	//all are now dimensionalised
	double P_sed = calculate_sedimentation_efficiency(gravity_angle,pipe_radius,pipe_length,mean_pipe_velocity,time_in_branch);
	double P_dif = calculate_brownian_efficiency(pipe_radius,pipe_length,mean_pipe_velocity,time_in_branch);
	double P_imp = calculate_impaction_efficiency(branching_angle,stokes_number,reynolds_number);

	// hmmmm
	//if(P_imp < 1.0)
	//	P_imp *= 1e3;

/*
	std::cout << "pipe_angle_from_gravity = " << pipe_angle_from_gravity << std::endl;
	std::cout << "P_sed = " << P_sed << std::endl;
	std::cout << "P_dif = " << P_dif << std::endl;
	std::cout << "P_imp = " << P_imp << std::endl;
*/

	double total_efficiency = 0.;
	if(mean_pipe_velocity < 1e-14)	//no impaction
		total_efficiency = P_sed + P_dif - P_sed*P_dif;
	else
		total_efficiency = P_sed + P_dif + P_imp - P_sed*P_dif - P_sed*P_imp - P_dif*P_imp + P_sed*P_dif*P_imp;

	return total_efficiency;
}

double HofmannParticleDeposition::calculate_sedimentation_efficiency(double pipe_angle_from_gravity,double pipe_radius,double pipe_length,double mean_pipe_velocity, double time_in_branch)
{
	if(mean_pipe_velocity > 1e-14)
		return 1 - exp(-4*gravity_magnitude*cunningham_correction_factor*particle_density*pow(particle_radius,2.0)*pipe_length*cos(pipe_angle_from_gravity)
								/(9*M_PI*viscosity*pipe_radius*mean_pipe_velocity));
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


double HofmannParticleDeposition::calculate_brownian_efficiency(double pipe_radius, double pipe_length,double mean_pipe_velocity, double time_in_branch)
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

double HofmannParticleDeposition::calculate_impaction_efficiency(double branching_angle,double stokes_number, double reynolds_number)
{

	if(es->parameters.get<unsigned int>("impaction_type") == 1)
	{
		if(stokes_number < 0.04)
		{
			return 0.000654*exp(55.7*pow(stokes_number,0.954)) * pow(reynolds_number,1./3.) * sin(branching_angle);
		}
		else
		{
			return (0.19-0.193*exp(-9.5*pow(stokes_number,1.565))) * pow(reynolds_number,1./3.) * sin(branching_angle);
		}
	}
	else
	{
		if(branching_angle*stokes_number < 1.0)
			return 1. - 2./M_PI*acos(branching_angle*stokes_number) + 1./M_PI*sin(2*acos(branching_angle*stokes_number));
		else
			return 1.;

	}
}


double HofmannParticleDeposition::calculate_deposition_in_alveoli(double time_spent_in_alveoli)
{
	// to calculate the efficiency of a branch we 
	// 1 - calculate the sedimentation efficiency of a whole pipe
	// 2 - calculate the brownian efficiency of a whole pipe
	// 3 - calculate the efficiency in a bend pipe with half the length of the pipe.

	//std::cout << "branch_flow = " << branch_flow << std::endl;

	time_spent_in_alveoli *= es->parameters.get<double>("time_scale_factor");

	//all are now dimensionalised
	double P_sed = calculate_sedimentation_efficiency_in_sphere(time_spent_in_alveoli);
	double P_dif = calculate_brownian_efficiency_in_sphere(time_spent_in_alveoli);

/*
	std::cout << "pipe_angle_from_gravity = " << pipe_angle_from_gravity << std::endl;
	std::cout << "P_sed = " << P_sed << std::endl;
	std::cout << "P_dif = " << P_dif << std::endl;
	std::cout << "P_imp = " << P_imp << std::endl;
*/

	double total_efficiency = 0.;
	total_efficiency = P_sed + P_dif - P_sed*P_dif;

	return total_efficiency;
}


double HofmannParticleDeposition::calculate_sedimentation_efficiency_in_sphere(double time_spent_in_alveoli)
{
	if(time_spent_in_alveoli < 2*alveolus_radius/sedimentation_velocity)
		return 1./2.*sedimentation_velocity*time_spent_in_alveoli/2./alveolus_radius*
						(3. - pow(sedimentation_velocity*time_spent_in_alveoli/(2.*alveolus_radius),2.));
	else
		return 1.;
}


double HofmannParticleDeposition::calculate_brownian_efficiency_in_sphere(double time_spent_in_alveoli)
{
	double fraction = 1.;
	unsigned int num_terms = 5;

	for(unsigned int i=1;i <=num_terms; i++)
	{
		fraction -= 6./pow(M_PI,2.)/pow(i,2.)*exp(-diffusion_coefficient*pow(i*M_PI/alveolus_radius,2.)*time_spent_in_alveoli);	
	}

	return fraction;
}

//calculate the deposition efficiency for a particle entering the trachea at time_in
bool HofmannParticleDeposition::calculate_flow_at_time(double& flow, double time_in, unsigned int branch_no)
{

	//std::cout << "time_in = " << time_in << std::endl;
	for(unsigned int i=0; i<times.size(); i++)
	{
		if(time_in < times[i])
		{
			double delta_t = ( times[i] - time_in)/(times[i] - times[i-1]);
			double flow_i = flow_function[i][branch_no];
			double flow_i_minus_1 = flow_function[i-1][branch_no];
			flow = flow_i + delta_t * (flow_i - flow_i_minus_1);
			return true;
		}
	}

	// if we get here it means we are at the end of the time so we can decide whether 
	if(end_sim_at_end_time)
	{
		return false;
	}
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
		else if(time_in < period + es->parameters.get<double>("hofmann_breath_hold"))
		{
			time_flow = 0.0;
		}
		else if(time_in < 2*period + es->parameters.get<double>("hofmann_breath_hold"))
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


void HofmannParticleDeposition::construct_temporary_alveolated_1d_tree (unsigned int parent_id)
{
	//std::cout << "constructing temporary alveolated 1d tree" << std::endl;

	//make sure alveolated tree is empty
	temporary_alveolar_tree.resize(0);

	// okay so we want to use element_data as a starting point
	// mesh uses the order of cell_vertices to add the segments with no reference
	// to element_data

	//tree params
	double bifurcation_angle = M_PI/4.0;
	double length_ratio = 1.0/1.25;	//halves the length of segments in each generation
	unsigned int num_alveolar_generations = es->parameters.get<unsigned int> ("num_alveolar_generations");
	//const double length_diam_ratio = 3.0;
	double length_diam_ratio = es->parameters.get<double> ("alveolar_length_diam_ratio");
	double left_length_ratio = length_diam_ratio;	//halves the length of segments in each generation
	double right_length_ratio = length_diam_ratio;	//halves the length of segments in each generation
	double alveolar_diameter = es->parameters.get<double> ("alveolar_diameter")/es->parameters.get<double> ("length_scale");
  const MeshBase& mesh = es->get_mesh();



	//std::cout << "we found an open daughter branch = " << parent_id << std::endl;
	num_alveolar_generations = airway_data[parent_id].get_num_alveolar_generations();

	// check the radius is not smaller than alveolar sac radius
	if(2*airway_data[parent_id].get_radius() < alveolar_diameter)
	{
		std::cout << "error, diam of terminal bronchiole (" << 2*airway_data[parent_id].get_radius()
							<< "is less than alveolar diameter(" << alveolar_diameter << ".. EXITING" << std::endl;
		std::exit(0);				
	}

	// calculate the ratios of diameters based on stepwise factor decrease
	length_ratio = pow(alveolar_diameter/(2*airway_data[parent_id].get_radius()),1./num_alveolar_generations);
	left_length_ratio = length_ratio;	//halves the length of segments in each generation
	right_length_ratio = length_ratio;	//halves the length of segments in each generation
	

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
			//std::cout << "hmmm" << std::endl;
			gravity(0) = es->parameters.get<double>("gravity_x");
			gravity(1) = es->parameters.get<double>("gravity_y");
			gravity(2) = es->parameters.get<double>("gravity_z");
			//std::cout << "gravity = " << gravity << std::endl;
		}
		else
		{
			//std::cout << "hmmm" << std::endl;
			gravity(0) = es->parameters.get<double>("gravity_x");
			gravity(1) = es->parameters.get<double>("gravity_y");
			//std::cout << "gravity = " << gravity << std::endl;	
		}
	}
	gravity = gravity.unit();

	//std::cout << "0" << std::endl;
	//helpful variables
	unsigned int parent_segment = 0;	//this is the segment we start on
	Point p0,p1;

	const Elem* elem = mesh.elem(parent_id);

	p0 = *elem->get_node(0);
	p1 = *elem->get_node(1);
	Point direction = p1 - p0;
	direction = direction.unit();	
	//calculate the gravity angle
	double gravity_angle = acos(direction*gravity);
	//gravity_angles[current_el_idx] = gravity_angle;		

	//std::cout << "1" << std::endl;
	//std::cout << "num alveolar_generations = " << num_alveolar_generations << std::endl;
	// the only difference is that when we are doing the first generation we
	// need to get data from the correct place
	// at the first generation j is only ever == 0 so this helps

	// temporary_alveolar_tree needs to have 0 - super_parent_idx 1 - parent idx, 2 - daughter_1, 3 - daughter_2, 
	// 4 - flow_fraction, 5 - length, 6 - radius, 7 - generation,
	// 8 - gravity_angle, 9 - bifurcation_angle, 10 - direction_x, 11 - direction_y, 12 - direction_z
	// 13 - weight in alveoli, 14 - time in alveoli
	
	// do parent
	temporary_alveolar_tree.push_back(std::vector<double>(15,-1));
	temporary_alveolar_tree[0][0] = parent_id;
	temporary_alveolar_tree[0][1] = -1;
	temporary_alveolar_tree[0][2] = -1;
	temporary_alveolar_tree[0][3] = -1;
	temporary_alveolar_tree[0][4] = 1;
	temporary_alveolar_tree[0][5] = airway_data[parent_id].get_length();
	temporary_alveolar_tree[0][6] = airway_data[parent_id].get_radius();
	temporary_alveolar_tree[0][7] = 0;
	temporary_alveolar_tree[0][8] = gravity_angle;
	temporary_alveolar_tree[0][9] = -1;
	temporary_alveolar_tree[0][10] = direction(0);
	temporary_alveolar_tree[0][11] = direction(1);
	temporary_alveolar_tree[0][12] = direction(2);
	temporary_alveolar_tree[0][13] = 0;
	temporary_alveolar_tree[0][14] = 0;

	parent_segment = 0;
	//std::cout << "2" << std::endl;

	// start with parent segment
	// at each generation parent_segment is the id of the first parent in the preceding generation
	// loop over each of these parents adding two branches
	// the id of the 
	for(unsigned int i=1; i<num_alveolar_generations + 1; i++)
	{

	//std::cout << "3-" << i << std::endl;
		//for each parent we generate two new points and two new segments
		for(unsigned int j=0;j<pow(2,i-1); j++)
		{
			//we need to find the direction of this segment
			direction = Point(temporary_alveolar_tree[parent_segment + j][10],
												temporary_alveolar_tree[parent_segment + j][11],
												temporary_alveolar_tree[parent_segment + j][12]);
			Point unit_direction = direction.unit();

			// rotate to new direction
			Point p0_direction;
			p0_direction(0) = cos(bifurcation_angle)*unit_direction(0) + sin(bifurcation_angle)* unit_direction(2);
			p0_direction(2) = -sin(bifurcation_angle)*unit_direction(0) + cos(bifurcation_angle)* unit_direction(2);
			p0_direction = p0_direction.unit();

			Point p1_direction;
			p1_direction(0) = cos(-bifurcation_angle)*unit_direction(0) + sin(-bifurcation_angle)* unit_direction(2);
			p1_direction(2) = -sin(-bifurcation_angle)*unit_direction(0) + cos(-bifurcation_angle)* unit_direction(2);
			p1_direction = p1_direction.unit();

			double length_1 = temporary_alveolar_tree[parent_segment + j][5]*left_length_ratio;
			double length_2 = temporary_alveolar_tree[parent_segment + j][5]*right_length_ratio;

			//update the data
			temporary_alveolar_tree.push_back(std::vector<double>(15,-1));
			temporary_alveolar_tree[parent_segment + j][2] = temporary_alveolar_tree.size() - 1;	// i.e. the latest branch
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][0] = parent_id;
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][1] = parent_segment + j;
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][4] = temporary_alveolar_tree[parent_segment + j][4] * (pow(left_length_ratio,2.0)/(pow(left_length_ratio,2.0) + pow(right_length_ratio,2.0)));
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][5] = length_1;
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][6] = temporary_alveolar_tree[parent_segment + j][5]*left_length_ratio;
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][7] = i;
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][8] = acos(p0_direction*gravity); //radius not diam
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][9] = bifurcation_angle; //generation
			//direction
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][10] = p0_direction(0);	//the flow rate in this tube for the current time step
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][11] = p0_direction(1);	//the efficiency in this tube for the current time step
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][12] = p0_direction(2);	//alveolar generation number
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][13] = 0;	//the efficiency in this tube for the current time step
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][14] = 0;	//alveolar generation number

			//update the data
			temporary_alveolar_tree.push_back(std::vector<double>(15,-1));
			temporary_alveolar_tree[parent_segment + j][3] = temporary_alveolar_tree.size() - 1;
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][0] = parent_id;
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][1] = parent_segment + j;
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][4] = temporary_alveolar_tree[parent_segment + j][4] * (pow(right_length_ratio,2.0)/(pow(left_length_ratio,2.0) + pow(right_length_ratio,2.0)));
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][5] = length_2;
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][6] = temporary_alveolar_tree[parent_segment + j][5]*right_length_ratio;
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][7] = i;
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][8] = acos(p1_direction*gravity); //radius not diam
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][9] = -bifurcation_angle; //generation
			//direction
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][10] = p1_direction(0);	//the flow rate in this tube for the current time step
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][11] = p1_direction(1);	//the efficiency in this tube for the current time step
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][12] = p1_direction(2);	//alveolar generation number
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][13] = 0;	//the efficiency in this tube for the current time step
			temporary_alveolar_tree[temporary_alveolar_tree.size() - 1][14] = 0;	//alveolar generation number

		}

		// if we are at the first we need to make the jump to the end
		parent_segment 	  += pow(2,i-1);

	}

	//std::cout << "4" << std::endl;


}
