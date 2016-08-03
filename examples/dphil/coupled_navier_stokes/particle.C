
#include "particle.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

Particle::Particle (EquationSystems& es_in,Point& p, const Elem* element, unsigned int _particle_id, double _perturbation_magnitude): 
		es (&es_in),on_wall(false),exited(false),exit_surface(0),perturbation_magnitude(_perturbation_magnitude)

{

	TransientLinearImplicitSystem * system;
	system =
  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");

	position = p;

	//const MeshBase& mesh = es->get_mesh();	// unused

	current_elem = element;	//probz illegal but oh well

	threed = es->parameters.get<bool>("threed");
	current_velocity = compute_velocity();
	current_particle_velocity = current_velocity * es->parameters.get<double>("initial_particle_velocity_scaling"); //current_velocity;
	maximum_timestep = current_elem->hmin() / current_velocity.size();

	exit_surface = 0;
	time_exited = -1;
	particle_id = _particle_id;
	entrance_time = system->time;
	broken = false;

	double lambda = es->parameters.get<double>("mean_free_path");
	double particle_diameter = es->parameters.get<double>("particle_diameter");
	cunningham_correction_factor = 1. + 2.*lambda/particle_diameter*
																	(1.257 + 0.4*exp(-1.1*(particle_diameter/(2*lambda))));


  constant_drag_force = 18. / 24. * es->parameters.get<double>("particle_air_viscosity") / 
														(es->parameters.get<double>("particle_density") * 
														pow(es->parameters.get<double>("particle_diameter"),2.0) *
														cunningham_correction_factor);


	sedimentation = es->parameters.get<bool>("particle_sedimentation");
	impaction = es->parameters.get<bool>("particle_impaction");
	drag = es->parameters.get<bool>("particle_drag");

	//std::cout << "current_velocity = " << current_velocity << std::endl;
	//std::cout << "elem number = " << current_elem->id() << std::endl;

}


double Particle::particle_reynolds_number ()
{
	// is it correct to have a reynolds nnumber for different directions
	// units converts to SI for parameters
	double Re = (current_particle_velocity - current_velocity).size() * es->parameters.get<double>("particle_velocity_units");
	//std::cout << "Re = " << Re << std::endl;
	Re *= es->parameters.get<double>("particle_air_density") * es->parameters.get<double>("particle_diameter") /
					es->parameters.get<double>("particle_air_viscosity");	//if we are not using SI unit here we are in trouble, but okay
	
	//std::cout << "particle_air_viscosity = " << es->parameters.get<double>("particle_air_viscosity") << std::endl;
	//std::cout << "Re = " << Re << std::endl;
	return Re;
	
}

double Particle::drag_coeff (double reynolds_number)
{

	//should depend on reynolds number, we use the values for reynolds number 10-100
	// see Morsi 1972, via Lambert 2011
	double a1 = 0.;
	double a2 = 0.;
	double a3 = 0.;

	if(reynolds_number < 0.1)
	{
		a1 = 0.;
		a2 = 24.0;
		a3 = 0.;
	}
	else if(reynolds_number < 1.0)
	{
		a1 = 3.69;
		a2 = 22.73;
		a3 = 0.0903;
	}
	else if(reynolds_number < 10.0)
	{
		a1 = 1.222;
		a2 = 29.1667;
		a3 = -3.8889;
	}
	else if(reynolds_number < 100.0)
	{
		a1 = 0.6167;
		a2 = 46.5;
		a3 = -116.67;
	}
	else if(reynolds_number < 1000.0)
	{
		a1 = 0.3644;
		a2 = 98.33;
		a3 = -2778;
	}
	else
	{
		a1 = 0.357;
		a2 = 148.62;
		a3 = -4.75e4;
	}

	if(fabs(reynolds_number) > 1e-10)
		return a1 + a2/reynolds_number + a3/pow(reynolds_number,2.0);
	else
		return 0.;

}

NumberVectorValue Particle::compute_velocity ()
{

  // Get a reference to the Stokes system object.
	TransientLinearImplicitSystem * system;
	system =
  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");

	// Get variable numbers
  const unsigned int u_var = system->variable_number ("u");
  const unsigned int v_var = system->variable_number ("v");
	unsigned int w_var = 0;
	if(threed)
		w_var = system->variable_number ("w");

  // Compute the velocity
  Number   u = 0., v = 0., w = 0.;
  // From the previous Newton iterate:
  u = system->point_value(u_var, position, *current_elem);
  v = system->point_value(v_var, position, *current_elem);
	if(threed)
    w = system->point_value(w_var, position, *current_elem);

  NumberVectorValue U;
	if(threed)
		U = NumberVectorValue(u, v, w);
	else
		U = NumberVectorValue(u, v);

	return U;
}

// move particles over one timestep
void Particle::try_and_move () 
{

	if(!on_wall || !exited)
	{
		std::cout << "moving particle " << particle_id << std::endl;
		move();
	}

}


// move particles over one timestep, using forward euler
NumberVectorValue Particle::compute_particle_velocity (NumberVectorValue velocity) 
{
	NumberVectorValue particle_velocity(0.,0.,0.);
	const double dt = es->parameters.get<Real>("dt");


	Point gravity(0.,0.,0.);
	if(es->parameters.get<unsigned int>("gravity_type") == 0)
	{
		if(threed)
			gravity(2) = -1.;
		else
			gravity(1) = -1.;					
	}
	else if(es->parameters.get<unsigned int>("gravity_type") == 1 || es->parameters.get<unsigned int>("gravity_type") == 2)
	{
		if(threed)
		{
			gravity(0) = es->parameters.get<double>("gravity_x");
			gravity(1) = es->parameters.get<double>("gravity_y");
			gravity(2) = es->parameters.get<double>("gravity_z");
		}
		else
		{
			gravity(0) = es->parameters.get<double>("gravity_x");
			gravity(1) = es->parameters.get<double>("gravity_y");	
		}
	}
	gravity *= es->parameters.get<double>("gravity");


	//std::cout << "gravity = " << gravity << std::endl;

	if(!sedimentation && !impaction && !drag)
	{
		particle_velocity = velocity;
	}
	else if(impaction)
	{

		//double Re_p = 0.;	// unused
		//double C_D = 0.;	// unused
		double F_D = 0.;
		
		
		if(drag)
		{
	

			double Re_p = particle_reynolds_number ();
			double C_D = drag_coeff (Re_p);
			double F_D = constant_drag_force * C_D * Re_p;

			particle_velocity = current_particle_velocity/(1. + dt*F_D);

			//std::cout << "Re_p = " << Re_p << std::endl;
			//std::cout << "constant_drag_force = " << constant_drag_force << std::endl;
			//std::cout << "C_D = " << C_D << std::endl;
			//std::cout << "F_D = " << F_D << std::endl; 
			particle_velocity += dt * F_D * (velocity)/(1. + dt*F_D);
		}
		else
		{
			particle_velocity = current_particle_velocity;
		}

		if(sedimentation)
		{

			double particle_density = es->parameters.get<double>("particle_density");
			double particle_air_density = es->parameters.get<double>("particle_air_density");

			if(drag)
			{
				particle_velocity += dt * gravity * (1. - particle_air_density/particle_density) 
																/ es->parameters.get<double>("particle_velocity_units") / (1. + dt*F_D);
			}
			else
			{
				particle_velocity += dt * gravity * (1. - particle_air_density/particle_density) 
																/ es->parameters.get<double>("particle_velocity_units");
			}

		}
	}
	else if(!impaction)
	{

		double Re_p = particle_reynolds_number ();
		double C_D = drag_coeff (Re_p);
		double F_D = constant_drag_force * C_D * Re_p;

		if(fabs(Re_p) < 1e-10)
		{
			std::cout << "can't do without impaction if the relative reynolds number is 0... EXITING" << std::endl;
			std::exit(0);
		}

		//std::cout << "Re_p = " << Re_p << std::endl;
		//std::cout << "constant_drag_force = " << constant_drag_force << std::endl;
		//std::cout << "C_D = " << C_D << std::endl;
		//std::cout << "F_D = " << F_D << std::endl; 

		if(drag)
		{
			particle_velocity = velocity;
		}
		
		if(sedimentation)
		{


			double particle_density = es->parameters.get<double>("particle_density");
			double particle_air_density = es->parameters.get<double>("particle_air_density");
			particle_velocity += gravity / F_D * (1. - particle_air_density/particle_density) 
															/ es->parameters.get<double>("particle_velocity_units");

		}
	}

	return particle_velocity;

}


// move particles over one timestep
void Particle::move () 
{

	TransientLinearImplicitSystem * system;
	system =
  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");

	//const MeshBase& mesh = es->get_mesh();	// unused

	current_velocity = compute_velocity(); //calculate velocity for the current timestep
	//std::cout << "current_velocity = " << current_velocity << std::endl;
	current_particle_velocity = compute_particle_velocity(current_velocity); //calculate velocity for the current timestep
	//std::cout << "current_particle_velocity = " << current_particle_velocity << std::endl;

	const double dt = es->parameters.get<Real>("dt");
	double local_time = system->time;
	double local_end_time = local_time + dt;
	double particle_dt = 0.;

	unsigned int count = 0;

	//std::cout << "new particle" << std::endl;

	while(local_time < local_end_time - 1e-10 && !on_wall && !exited)
	{
		count++;
		if(count > 10)
			exit(0);
		//std::cout << "local_time = " << local_time << std::endl;
		//std::cout << "dt = " << dt << std::endl;
		//std::cout << "local_end_time = " << local_end_time << std::endl;
		if(fabs(current_particle_velocity.size()) > 1e-10)
			maximum_timestep = current_elem->hmin()/current_particle_velocity.size();
		else
			maximum_timestep = dt;

		//std::cout << "current_elem->hmin()= " << current_elem->hmin() << std::endl;
		//if(maximum_timestep < dt)
		//	particle_dt = maximum_timestep;
		//else
			particle_dt = dt;
			
		// check not past end of timestep
		if(local_time + particle_dt >= local_end_time)
			particle_dt = local_end_time - local_time;

		//std::cout << "particle_dt = " << particle_dt << std::endl; 
		local_time += particle_dt;

		Point new_position = position + particle_dt*(current_particle_velocity);

		double x = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * es->parameters.get<Real>("brownian_motion_magnitude");
		double y = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * es->parameters.get<Real>("brownian_motion_magnitude");
		double z = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * es->parameters.get<Real>("brownian_motion_magnitude");
		
		Point perturbation(x,y,z);
		if(!threed)
			perturbation(2) = 0.;
		new_position = new_position + perturbation;

		//IDEALLY:
		//1 - check element
		//2 - check neighbors
		//3 - check walls of element
		//4 - find which neighbor we go through
		//5 - go back to 2

		//FOR NOW:
		//1 - check element
		//2 - check neighbors
		//3 - check all elements
		//4 - havent found anything freeze particle and say it is stuck


		//bool temporary_method = false;	// unused
		bool element_found = false;
		//now find what element we are in, search neighbours

		
		while(!element_found)
		{
			//find new position - perturbed if we are doing this again
			// note after perturbation we may be in a different element... hmmm
			//std::cout << "new_position = " << new_position << std::endl;
			//std::cout << "position = " << position << std::endl;

			std::vector<unsigned int> elements_checked;

			if(current_elem->contains_point(new_position))
			{
				//do nothing already in correct element, move point and go onto next time step.
				position = new_position;
				element_found = true;
			}
			else	//check neighbors
			{
				double old_s_param = 0.;
				check_neighbors(current_elem,new_position,old_s_param,elements_checked,element_found);				
			}
		}
	}
}


// in order to avoid going back and forth between elements we pass in the previous s_param and make sure
// that the new s_param is greater than the old one
void Particle::check_neighbors(const Elem* element, Point& new_position, double old_s_param, std::vector<unsigned int>& elements_checked,  bool& element_found)
{
	std::cout << "in_check_neighbors" << std::endl;

	TransientLinearImplicitSystem * system;
	system =
  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");

	const MeshBase& mesh = es->get_mesh();

	elements_checked.push_back(element->id());
	

	unsigned int n_neighbors =	element->n_neighbors();

	// check if particle is in one of the neighbors then we are done
	for(unsigned int i=0; i<n_neighbors; i++)
	{
		Elem* neighbor = element->neighbor(i);	//find_neighbors may have to be called
		if(neighbor != NULL)
		{
			if(neighbor->contains_point(new_position))
			{	
				
				current_elem = neighbor;
				position = new_position;
				element_found = true;
				break;
			}
		}
	}

	if(!element_found)
	{
		//std::cout << "n_neighbors = " << n_neighbors << std::endl;
		//std::cout << "elem id = " << element->id() << std::endl;
		for(unsigned int i=0; i<n_neighbors; i++)
		{
			Elem* neighbor = element->neighbor(i);	//find_neighbors may have to be called

			//std::cout << "i = " << i << std::endl;
			double s;
			// if does through side and neighbor in question has not been checked already
			if(goes_through_side(new_position,i,element,old_s_param,s) && 
				(neighbor == NULL || find(elements_checked.begin(),elements_checked.end(),neighbor->id()) == elements_checked.end()) )
			{
				//std::cout << "hey" << std::endl;
				//goes through a side of the element into the abyss.. and we are done.
				if(neighbor == NULL)
				{
					//std::cout << "grrrreat" << std::endl;
					element_found = true;
					std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(element,i);

					if(boundary_ids.size() > 0) 
					{ 
						int boundary_id = boundary_ids[0];	// should only have one
						//std::cout << "should be here b_id = " << boundary_id << std::endl;
						//std::cout << "new_position = " << new_position << std::endl;
						//std::cout << "position = " << position << std::endl;
						//std::cout << "s = " << s << std::endl;
						if(boundary_id == -1)
						{
							on_wall = true;
							exit_surface = -1;
							time_exited = system->time;
							position = position + s*(new_position - position);
							break;
						}
						else
						{
							exited = true;
							exit_surface = boundary_id;
							time_exited = system->time;
							position = position + s*(new_position - position);
							break;
						}
					}		 
				}
				else
				{
					//std::cout << "hi" << std::endl;
					check_neighbors(neighbor,new_position,s,elements_checked,element_found);
					if(element_found)
						break;
				}
			}				
		}
	}

	if(!element_found)
	{
		double x = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * perturbation_magnitude;
		double y = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * perturbation_magnitude;
		double z = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * perturbation_magnitude;
		
		Point perturbation(x,y,z);
		if(!threed)
			perturbation(2) = 0.;
		
		//std::cout << "Next element not found, probably because on edge or corner, perturbing new position by " << perturbation << std::endl;
		//new_position = new_position + perturbation;

		std::cout << "Next element not found, probably because on edge or corner, assuming on the wall." << std::endl;
		on_wall = true;
		exit_surface = -1;
		time_exited = system->time;
		element_found = true;
		broken = true;
	}
}

//returns true if goes through side side_number, also returns the parameter point at which it went through the wall
bool Particle::goes_through_side(Point new_position,int side_number,const Elem* element, double old_s_param, double& s_param)
{

	AutoPtr<Elem> side = element->build_side(side_number);

	if(threed)
	{
		//construct normal
		Point point_1 = side->point(0);
		Point point_2 = side->point(1);
		Point point_3 = side->point(2);

		Point l_1 = point_1 - point_2;
		Point l_2 = point_3 - point_2;

		Point normal(l_1(1)*l_2(2)-l_1(2)*l_2(1),
								l_1(2)*l_2(0)-l_1(0)*l_2(2),
								l_1(0)*l_2(1)-l_1(1)*l_2(0));
	
		normal = normal.unit();

		Point direction = new_position - position;	

		//now we calculate the intersection point
		s_param = normal*(point_1 - position) / (normal * (new_position - position));

		//we want to ignore if it is going along an edge
		if(fabs(normal * (new_position - position)) < 1e-10)
			s_param = -666.;
		

		//std::cout << "s_param = " << s_param << std::endl;
		//std::cout << "normal = " << normal << std::endl;
		//std::cout << "direction = " << direction << std::endl;

		// if we need to go backwards it definitely doesn't intersect, if it is the same 
		// as the previous s_param then we are just going back through the same element side,
		// unless it is approximately zero then we are not sure 
		if(s_param < -1e-10 || (old_s_param > 1e-10 && s_param < old_s_param + 1e-5))
			return false;
		else
		{
			// now we check if a point just a bit back from the intersection point is in the element
			// if this is true then we do intersect this side, note there is a tolerance on contains point that we have to beat somehow
			if(element->contains_point(position + (s_param - 1e-4)*(new_position - position)))
			{
				//std::cout << "in the element " << element->id() << std::endl;
				return true;
			}
			else
			{
				//std::cout << "not in the element " << element->id() << std::endl;
				return false;
			}
		}
	}
	else
	{
		//construct normal
		Point point_1 = side->point(0);
		Point point_2 = side->point(1);

		//get line of side
		Point l_1 = point_2 - point_1;

		//swap the coords and make one negative to get the normal
		Point normal(l_1(1),-l_1(0),0.);
		normal = normal.unit();

		Point direction = new_position - position;

		//now we calculate the intersection point
		s_param = normal*(point_1 - position) / (normal * (new_position - position));
		if(fabs(normal * (new_position - position)) < 1e-10)
			s_param = -666.;

		//std::cout << "s_param = " << s_param << std::endl;

		// if we need to go backwards it definitely doesn't intersect, if it is the same 
		// as the previous s_param then we are just going back through the same element side
		if(s_param < -1e-10 || (old_s_param > 1e-10 && s_param < old_s_param + 1e-5))
		{
			//std::cout << "hum" << std::endl;
			return false;
		}
		else
		{
				//std::cout << "kei" << std::endl;
			// now we check if a point just a bit back from the intersection point is in the element
			// if this is true then we do intersect this side
			if(element->contains_point(position + (s_param - 1e-4)*(new_position - position)))
			{
				//std::cout << "what" << std::endl;
				return true;
			}
			else
			{
				//std::cout << "do you mean" << std::endl;
				return false;
			}
		}
	}

}


