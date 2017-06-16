
#include "particle.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

Particle::Particle (EquationSystems& es_in,Point& p, const Elem* element, unsigned int _particle_id, double _perturbation_magnitude): 
		es (&es_in),on_wall(false),exited(false),exit_surface(0),perturbation_magnitude(_perturbation_magnitude),failure_perturbation_magnitude(1e-1)

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


  constant_drag_force = 18. / 24. * es->parameters.get<double>("particle_air_viscosity") / (es->parameters.get<double>("particle_density") * 
														pow(es->parameters.get<double>("particle_diameter"),2.0) *
														cunningham_correction_factor);


	sedimentation = es->parameters.get<bool>("particle_sedimentation");
	impaction = es->parameters.get<bool>("particle_impaction");
	drag = es->parameters.get<bool>("particle_drag");


	std::cout << "cunningham correction factor = " << cunningham_correction_factor << std::endl;
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

	if(!on_wall || !exited || !broken)
	{
		//std::cout << "moving particle " << particle_id << std::endl;
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

	if(!sedimentation && !drag)
	{
		particle_velocity = velocity;
	}
	else if(drag)
	{

		
	
		double Re_p = particle_reynolds_number ();
		double C_D = drag_coeff (Re_p);
		double F_D = constant_drag_force * C_D * Re_p;

		particle_velocity = current_particle_velocity/(1. + dt*F_D) + dt * F_D * (velocity)/(1. + dt*F_D);

		//std::cout << "Re_p = " << Re_p << std::endl;
		//std::cout << "constant_drag_force = " << constant_drag_force << std::endl;
		//std::cout << "C_D = " << C_D << std::endl;
		//std::cout << "F_D = " << F_D << std::endl;


		if(sedimentation)
		{

			double particle_density = es->parameters.get<double>("particle_density");
			double particle_air_density = es->parameters.get<double>("particle_air_density");

			particle_velocity += dt * gravity * (1. - particle_air_density/particle_density) / es->parameters.get<double>("particle_velocity_units") / (1. + dt*F_D);

		}
	}
	else if(!drag)
	{



		particle_velocity = current_particle_velocity;
		
		if(sedimentation)
		{


			double particle_density = es->parameters.get<double>("particle_density");
			double particle_air_density = es->parameters.get<double>("particle_air_density");

			particle_velocity += dt * gravity * (1. - particle_air_density/particle_density) / es->parameters.get<double>("particle_velocity_units");
		}
	}

	return particle_velocity;

}


// move particles over one timestep
void Particle::move () 
{

	if(es->parameters.get<bool>("deposition_verbose"))
		std::cout << "MOVING" << std::endl;

	TransientLinearImplicitSystem * system;
	system =
  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");

	const MeshBase& mesh = es->get_mesh();

	current_velocity = compute_velocity(); //calculate velocity for the current timestep
	current_particle_velocity = compute_particle_velocity(current_velocity); //calculate velocity for the current timestep

	const double dt = es->parameters.get<Real>("dt");
	double local_time = system->time;
	double local_end_time = local_time + dt;
	double particle_dt = 0.;

	unsigned int count = 0;

	// move particle in sub time steps
	while(local_time < local_end_time - 1e-10 && !on_wall && !exited && !broken)
	{
		count++;
		// if there are too many sub time steps
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

		
		if(es->parameters.get<bool>("deposition_verbose"))
		{
			std::cout << "position = " << position << std::endl;
			std::cout << "new position = " << new_position << std::endl;
		}

		
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


				bool simple = true;

				// complex check neighbors
				if(!simple)
				{
					
					//check_neighbors(current_elem,new_position,old_s_param,elements_checked,element_found);				
				}
				else
				{

					const PointLocatorBase & pl = mesh.point_locator();
				
					if(es->parameters.get<bool>("deposition_verbose"))
						std::cout << "doing point locator" << std::endl;

					const Elem * element_found_pl = pl(new_position);

					
					if(es->parameters.get<bool>("deposition_verbose"))
						std::cout << "done" << std::endl;

					if(element_found_pl == NULL)
					{
						std::cout << "element not found using point_locator..................................." << std::endl;
						check_neighbors(current_elem,new_position,old_s_param,elements_checked,element_found);	
					}
					else
					{
						current_elem = element_found_pl;
						position = new_position;
						element_found = true;
					}


					// simple check neighbors
					//check_neighbors_simple(current_elem,new_position,elements_checked,element_found);

					if(!element_found)
					{
						std::cout << "element not found using simple..................................." << std::endl;
						std::cout << elements_checked.size() << " elements checked" << std::endl;
					
						check_neighbors(current_elem,new_position,old_s_param,elements_checked,element_found);	

					}
					
				}
			}

			// check if particle is close to the wall
			if(element_found)
			{
				close_to_wall();
			}
		}
	}
}


// check if the particle is close enough to the wall to be deposited
// 1 - check if particle is in motion
// 2 - check if particle is in an element that has a wall bdy
// 3 - check how far particle is from the wall bdy
void Particle::close_to_wall()
{
	

  	// Get a reference to the Stokes system object.
	TransientLinearImplicitSystem * system;
	system =
  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");
	const MeshBase& mesh = es->get_mesh();

	// 1 - check is particle is in motion, i.e. not on wall, exited or broken
	if(!on_wall || !exited || !broken)
	{
		// 2 - check if particle is in an element that has a wall bdy

		unsigned int n_neighbors =	current_elem->n_neighbors();
		for(unsigned int i=0; i<n_neighbors; i++)
		{
			std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(current_elem,i);
			
			if(boundary_ids.size() > 0) 
			{ 
				int boundary_id = boundary_ids[0];	// should only have one
				if(boundary_id == -1)
				{
					//element has a side that is on the wall bdy
					// 3 - calculate the distance from the particle to this side

					// distance from point to a plane is v . n 
					// where n is the unit normal of the plane 
					// and v is a vector from the point to a point on the plane
					
					AutoPtr<Elem> side = current_elem->build_side(i);

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

						Point vector_to_plane = position - point_1;
				
						double distance = normal*vector_to_plane;
						
						
						// if the distance is less than the radius then it is deposited
						if(fabs(distance) < es->parameters.get<double>("particle_diameter")/2.0)
						{
							std::cout << "particle is within a radius of the wall and will be deposited." << std::endl;

							// need to calculate where it will be deposited
							// move it distance in the direction of the normal
							// note: we are not sure if the normal is outward or inward but the sign of distance will sort that out.
							position = position + distance*normal;
							on_wall = true;
							deposited_close_to_wall = true;
							exit_surface = -1;
							time_exited = system->time;
						}
					}
					else
					{
						std::cout << "close to wall not implemented for 2D sims" << std::endl;
					}	

				}
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

	// if not in a neighbor
	if(!element_found)
	{
		
		// look for where the particle went once from the current position and then again from a perturbation
		unsigned int perturbation_no = 0;
		unsigned int num_perturbations = 0;
		while(perturbation_no <= num_perturbations && !element_found)
		{

			double x = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * 1.e-4;
			double y = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * 1.e-4;
			double z = (2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1) * 1.e-4;
		
			Point perturbation(x,y,z);
			if(!threed)
				perturbation(2) = 0.;

			std::cout << "perturbation number " << perturbation_no << std::endl;
			Point centroid = element->centroid();
			position = position + 0.1*perturbation_no*(centroid - position);
			new_position = new_position + 0.1*perturbation_no*(centroid - position);
			position = centroid;
			std::cout << "looking from position " << position << std::endl;
			std::cout << "looking to new_position " << new_position << std::endl;
			//if(perturbation_no > 0)
			//{
			//	new_position = new_position - 0.5*(new_position-position);
			//	std::cout << "looking to new_position " << new_position << std::endl;
			//}
			if(perturbation_no > 0)
				std::cout << "ATTENTION PERTURBING>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

			//std::cout << "n_neighbors = " << n_neighbors << std::endl;
			//std::cout << "elem id = " << element->id() << std::endl;
			// loop over the sides (same as  neighbors in libmesh)
			for(unsigned int i=0; i<n_neighbors; i++)
			{
				Elem* neighbor = element->neighbor(i);	//find_neighbors may have to be called

				//std::cout << "i = " << i << std::endl;
				double s;
				// if goes through side and neighbor in question has not been checked already

				//std::cout << "goes through side = " << goes_through_side(new_position,i,element,old_s_param,s) << std::endl;
				if(goes_through_side(new_position,i,element,old_s_param,s))// && 
					//(neighbor == NULL || find(elements_checked.begin(),elements_checked.end(),neighbor->id()) == elements_checked.end()) )
				{
					std::cout << "hey, went through side" << std::endl;
					//goes through a side of the element into the abyss.. and we are done.
					if(neighbor == NULL)
					{
						std::cout << "grrrreat found element" << std::endl;
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
								std::cout << "particle " << particle_id << " on wall" << std::endl;
								on_wall = true;
								exit_surface = -1;
								time_exited = system->time;
								position = position + s*(new_position - position);
								break;
							}
							else
							{
								std::cout << "particle " << particle_id << " exited" << std::endl;
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
						std::cout << "hi" << std::endl;
						// check that particle is inside this new element
						std::cout << "neighbor->id() = " << neighbor->id() << std::endl;
						if(neighbor->contains_point(position))
							std::cout << "neighbor contains point" << std::endl;
						else
							std::cout << "oops neighbor doesn't contain point" << std::endl;
						check_neighbors(neighbor,new_position,s,elements_checked,element_found);
						if(element_found)
							break;
					}
				}				
			}
			perturbation_no++;
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

		std::cout << "Next element not found, probably because on edge or corner, setting broken on surface -2." << std::endl;
		//on_wall = true;
		exit_surface = -2;
		time_exited = system->time;
		position = new_position;
		element_found = true;
		broken = true;
	}
}




// in order to avoid going back and forth between elements we pass in the previous s_param and make sure
// that the new s_param is greater than the old one
void Particle::check_neighbors_simple(const Elem* element, Point& new_position, std::vector<unsigned int>& elements_checked, bool& element_found)
{
	//std::cout << "in_check_neighbors_simple" << std::endl;

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
		if(neighbor != NULL && find(elements_checked.begin(),elements_checked.end(),neighbor->id()) == elements_checked.end())
		{
			if(neighbor->contains_point(new_position))
			{	
				
				current_elem = neighbor;
				position = new_position;
				element_found = true;
				break;
			}
			// if not found in neighbors check neighbors 
			// as long as centroid is less than 3x the distance from the original position
			else
			{
				Point centroid = neighbor->centroid();
				double distance_to_centroid = (centroid - position).size();
				double particle_distance = (new_position - position).size();
				//std::cout << "not found in centroid " << centroid << std::endl;
				//std::cout << "distance_to_centroid/particle_distance " << distance_to_centroid/particle_distance << std::endl;
				if(distance_to_centroid/particle_distance < 2.0)
				{
					check_neighbors_simple(neighbor,new_position,elements_checked,element_found);
					if(element_found)
						break;
				}
			}
		}
		
	}
}

//returns true if goes through side side_number, also returns the parameter point at which it went through the wall
bool Particle::goes_through_side(Point new_position,int side_number,const Elem* element, double old_s_param, double& s_param)
{
	unsigned int n_neighbors =	element->n_neighbors();
	std::cout << "side_number = " << side_number << std::endl;
	std::cout << "n_neighbors = " << n_neighbors << std::endl;

	const MeshBase& mesh = es->get_mesh();

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
		s_param = normal*(point_2 - position) / (normal * (new_position - position));

		//we want to ignore if it is going along an edge
		if(fabs(normal * (new_position - position)) < 1e-10)
			s_param = -666.;
		

		std::cout << "s_param = " << s_param << std::endl;
		//std::cout << "normal = " << normal << std::endl;
		//std::cout << "direction = " << direction << std::endl;

		// if we need to go backwards it definitely doesn't intersect, if it is the same 
		// as the previous s_param then we are just going back through the same element side,
		// unless it is approximately zero then we are not sure 
		//if(s_param < -1e-10 || (old_s_param > 1e-10 && s_param < old_s_param + 1e-5))
		if(s_param < -1e-10 || fabs(s_param) < 1e-10)// || (old_s_param > 1e-10 && s_param < old_s_param + 1e-5))
		{
			std::cout << "goddamn" << std::endl;
			return false;
		}
		else
		{
			// now we check if a point just a bit back from the intersection point is in the element
			// if this is true then we do intersect this side, note there is a tolerance on contains point that we have to beat somehow
			Point intersection_point = position + s_param*(new_position - position);
			std::cout << "intersection point " << intersection_point << std::endl;
				
			// old tolerance of 1e-5 didn't work so well
			//if(element->contains_point(position + (s_param - 1e-8)*(new_position - position)))
			if(element->contains_point(position))
				std::cout << "contains position point" << std::endl;
			else
				std::cout << "doesn't contain position point" << std::endl;


			std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(element,side_number);

			if(boundary_ids.size() > 0) 
			{ 

			}

			// check that the intersection point is on the side
			
			// if not on a boundary check if it's in a neighbor
			if(boundary_ids.size() == 0 && element->neighbor(side_number)->contains_point(intersection_point))
			{
				std::cout << "on the bdy and in the neighboring element " << element->neighbor(side_number)->id() << std::endl;
				position = intersection_point;
				return true;
			}
			// if it is on a boundary check that it's on the current element's boundary
			else if(boundary_ids.size() > 0 && element->contains_point(intersection_point))
			{
				std::cout << "on the bdy and in the element " << element->id() << std::endl;
				position = intersection_point;
				return true;
			}
			else
			{
				std::cout << "not in the element " << element->id() << std::endl;
				std::cout << "goddamn 2" << std::endl;
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
			if(element->contains_point(position + (s_param - 1e-8)*(new_position - position)))
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


