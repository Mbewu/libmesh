#ifndef INFLOW_BOUNDARY_FUNCTION_H
#define INFLOW_BOUNDARY_FUNCTION_H

#include "libmesh/function_base.h"
#include "libmesh/equation_systems.h"
#include "surface_boundary.h"

// want this to be able to handle square boundaries as well
// assume not centred on origin and just need "radius"

namespace libMesh
{

template <typename Output=Number>
class InflowBoundaryFunction : public FunctionBase<Output>
{
private:

  /**
   * The EquationSystems object that we're using here.
   */
  EquationSystems* es;
  SurfaceBoundary* surface_boundary_object;
	std::vector<double> time_values;
	std::vector<double> position_values;
	std::vector<std::vector<double> > axial_velocity_values;

public:

  /**
   * Constructor. don't need surface boundary object for projecting the exact solution to the whole domain
   */
  InflowBoundaryFunction(EquationSystems& _es, SurfaceBoundary* _surface_boundary_object=NULL, bool read_inflow_profile=false)
			:  es(&_es), surface_boundary_object(_surface_boundary_object)	
	{
		// only if using womersley do we need to read inflow profile, but you never know
		if(es->parameters.get<unsigned int> ("inflow_profile") == 2 || read_inflow_profile)
		{
			init_inflow_profile_values();
		}
	};

  /**
   * Returns a new deep copy of the function.
   */
  virtual AutoPtr<FunctionBase<Number> > clone () const
	{
		return AutoPtr<FunctionBase<Number> > (new InflowBoundaryFunction (*es,surface_boundary_object) );
	};

  /**
   * @returns the value at point \p p and time
   * \p time, which defaults to zero.
   */
  Number operator() (const Point& p,
		     const Real t=0)
	{


		const double x = p(0);
		const double y = p(1);
		const double z = p(2);
	
		double r = sqrt(x*x + y*y);
		double radius = es->parameters.get<double> ("radius");//boundary_radius[boundary_id]
		double normalisation_constant = 1.0;
		double area = 1.0;


		// hmmm we want like  the p(2) to equal 0 for the initial condition, 
		// assuming centroid is at p(2) = 0.		
		Point p_xy(x,y,z);

		bool along_line = false;
		if(fabs(x/y - 1) < 1e-10)
		{
			along_line = true;
			//std::cout << std::endl;
			//std::cout << "x = " << x << std::endl;
		}
		//std::cout << "X/Y = " << x/y << std::endl;

		
		//double r_before_normed = r/0.5;	// unused

		if(along_line)
		{
			//std::cout << "r_before = " << r << std::endl;
			//std::cout << "r_before_normed = " << r_before_normed << std::endl;
		}
		// using new definition
		if(surface_boundary_object != NULL)
		{
			r = surface_boundary_object->get_normalised_distance_from_centroid(p_xy);
			normalisation_constant = surface_boundary_object->get_unit_parabola_integral();
			area = surface_boundary_object->get_area();
			radius = 1.0;//surface_boundary_object->get_max_radius();

			//std::cout << "r = " << r << std::endl;
			//std::cout << "radius = " << radius << std::endl;

			if(along_line)
			{
				//std::cout << "r_after = " << r << std::endl;
				//std::cout << "diff: r_before_normed - r_after = " << r_before_normed - r  << std::endl;
			}
			//std::cout << "centroid = " << surface_boundary_object->get_centroid() << std::endl;

		}
		else	//if we are not on a boundary we need to make assumptions about centreline
		{
			// centred on line x,y=0
			if(es->parameters.get<bool> ("threed"))
			{
				radius = es->parameters.get<double> ("radius");
				r = sqrt(x*x + y*y)/radius;
			}
			else	//centred on y = radius
			{
				radius = es->parameters.get<double> ("radius");
				r = fabs(y - radius)/radius;
				//r = fabs(y-radius);
				//std::cout << "radius = " << radius << std::endl;
				//std::cout << "r = " << r << std::endl;
			}
		}



		//make nondimensional
		double velocity_magnitude = 0.;
		double flow_mag = 0.;
		if(!es->parameters.get<bool> ("prescribed_flow"))
		{
			//need to scale it to be nondimensional
			velocity_magnitude = es->parameters.get<double> ("velocity_mag_3d")/es->parameters.get<double> ("velocity_scale");
		}
		else
		{
			//need parabolic integral on surface... normalisation_constant, this is done later anyway
			//need flow in dimensionless units
			flow_mag = es->parameters.get<double> ("flow_mag_3d")/(es->parameters.get<double> ("velocity_scale")*pow(es->parameters.get<double> ("length_scale"),2.0));
			
		}
			

		


		double direction = 1.0;
		if(es->parameters.get<bool> ("reverse_inflow_profile"))
			direction *= -1;
	
		

		//std::cout << "r = " << r << std::endl;
		//std::cout << "radius = " << radius << std::endl;
		//std::cout << "p = " << p << std::endl;
		// 0 - constant, 1 - parabolic (tensor)
		if(es->parameters.get<unsigned int> ("inflow_profile") == 0)
		{
			if(es->parameters.get<bool> ("prescribed_flow"))
			{
				return direction/area*es->parameters.get<double>("time_scaling")*flow_mag;
			}
			else
			{
				return direction*es->parameters.get<double>("time_scaling")*velocity_magnitude;
			}
		}
		else if(es->parameters.get<unsigned int> ("inflow_profile") == 1)
		{
			if(es->parameters.get<bool> ("prescribed_flow"))
			{
				// not right but who fuckin cares..
				// better to have max v == 1
				//std::cout << "hi - " << direction / normalisation_constant * (pow(radius,2)-pow(r,2)) * es->parameters.get<double>("time_scaling") * flow_mag << std::endl;
				//std::cout << direction << std::endl;
				//std::cout << es->parameters.get<double>("time_scaling") << std::endl;
				//std::cout << flow_mag  << std::endl;
				//std::cout << normalisation_constant  << std::endl;
				//std::cout << (pow(radius,2)-pow(r,2))  << std::endl;
				return direction / normalisation_constant * (pow(radius,2)-pow(r,2)) * es->parameters.get<double>("time_scaling") * flow_mag;//*(pow(radius,2)-r*r);
			}
			else
			{
				return direction * (pow(radius,2)-pow(r,2))/pow(radius,2)*es->parameters.get<double>("time_scaling") * velocity_magnitude;
			}
		}
		else if(es->parameters.get<unsigned int> ("inflow_profile") == 2)
		{
			for(unsigned int i=0; i<time_values.size(); i++)
			{
				if(fabs(t-time_values[i]) < 1e-10)
				{
					//std::cout << "hey baby, time_values[i] = " << time_values[i] << std::endl;
					for(unsigned int j=0;j<position_values.size()-1; j++)
					{
						// replace p(1) with r
						if(r >= position_values[j] && r <= position_values[j+1] + 1e-4)
						{
							//std::cout << "r = " << r << std::endl;
							double distance = r - position_values[j];
							double axial_velocity = axial_velocity_values[i][j] + 
																				distance*(axial_velocity_values[i][j+1] - axial_velocity_values[i][j])
																				/(position_values[j+1]-position_values[j]);
							//std::cout << "velocity = " << axial_velocity << std::endl;
							return axial_velocity;			//i think we may need negative because outward normal

						}

						if(j==position_values.size()-2)
						{
							std::cout << "FAILURE WE HAVE GONE PAST THE END OF SPACE in prescribed flow... EXITING" << std::endl;
							std::cout << "r = " << r << std::endl;
							std::cout << "x(n-1) = " << position_values[j] << ", and x(n) = " << position_values[j+1] + 1e-4 << std::endl;
							std::exit(0);
							return 0;
						}
					}
					break;
				}
			}

			std::cout << "FAILURE WE HAVE GONE PAST THE END OF TIME in prescribed flow... EXITING" << std::endl;
			std::exit(0);
			return 0;
		
		}
		else
		{
			std::cout << "FAILURE: INFLOW PROFILE " << es->parameters.get<unsigned int> ("inflow_profile")  << " IS NOT SUPPORTED... EXITING" << std::endl;
			std::exit(0);			
			return 0;
		}

	};


  /**
   * @returns the value at point \p p and time
   * \p time, which defaults to zero.
   */
  Number calculate_velocity_magnitude(const Point& p,
		     const Real t=0.,
				 const Real time_scaling = 0)
	{


		const double x = p(0);
		const double y = p(1);
		const double z = p(2);
	
		double r = sqrt(x*x + y*y);
		double radius = es->parameters.get<double> ("radius");//boundary_radius[boundary_id]
		double normalisation_constant = 1.0;
		double area = 1.0;


		// hmmm we want like  the p(2) to equal 0 for the initial condition, 
		// assuming centroid is at p(2) = 0.		
		Point p_xy(x,y,z);

		bool along_line = false;
		if(fabs(x/y - 1) < 1e-10)
		{
			along_line = true;
			//std::cout << std::endl;
			//std::cout << "x = " << x << std::endl;
		}
		//std::cout << "X/Y = " << x/y << std::endl;

		
		//double r_before_normed = r/0.5;	// unused

		if(along_line)
		{
			//std::cout << "r_before = " << r << std::endl;
			//std::cout << "r_before_normed = " << r_before_normed << std::endl;
		}
		// using new definition
		if(surface_boundary_object != NULL)
		{
			r = surface_boundary_object->get_normalised_distance_from_centroid(p_xy);
			normalisation_constant = surface_boundary_object->get_unit_parabola_integral();
			area = surface_boundary_object->get_area();
			radius = 1.0;//surface_boundary_object->get_max_radius();

			//std::cout << "r = " << r << std::endl;
			//std::cout << "radius = " << radius << std::endl;

			if(along_line)
			{
				//std::cout << "r_after = " << r << std::endl;
				//std::cout << "diff: r_before_normed - r_after = " << r_before_normed - r  << std::endl;
			}
			//std::cout << "centroid = " << surface_boundary_object->get_centroid() << std::endl;

		}
		else	//if we are not on a boundary we need to make assumptions about centreline
		{
			// centred on line x,y=0
			if(es->parameters.get<bool> ("threed"))
			{
				radius = es->parameters.get<double> ("radius");
				r = sqrt(x*x + y*y)/radius;
			}
			else	//centred on y = radius
			{
				radius = es->parameters.get<double> ("radius");
				r = fabs(y - radius)/radius;
				//r = fabs(y-radius);
				//std::cout << "radius = " << radius << std::endl;
				//std::cout << "r = " << r << std::endl;
			}
		}



		//make nondimensional
		double velocity_magnitude = 0.;
		double flow_mag = 0.;
		if(!es->parameters.get<bool> ("prescribed_flow"))
		{
			//need to scale it to be nondimensional
			velocity_magnitude = es->parameters.get<double> ("velocity_mag_3d")/es->parameters.get<double> ("velocity_scale");
		}
		else
		{
			//need parabolic integral on surface... normalisation_constant, this is done later anyway
			//need flow in dimensionless units
			flow_mag = es->parameters.get<double> ("flow_mag_3d")/(es->parameters.get<double> ("velocity_scale")*pow(es->parameters.get<double> ("length_scale"),2.0));
			
		}
			

		


		double direction = 1.0;
		if(es->parameters.get<bool> ("reverse_inflow_profile"))
			direction *= -1;
	
		

		//std::cout << "r = " << r << std::endl;
		//std::cout << "radius = " << radius << std::endl;
		//std::cout << "p = " << p << std::endl;
		// 0 - constant, 1 - parabolic (tensor)
		if(es->parameters.get<unsigned int> ("inflow_profile") == 0)
		{
			if(es->parameters.get<bool> ("prescribed_flow"))
			{
				return direction/area*time_scaling*flow_mag;
			}
			else
			{
				return direction*time_scaling*velocity_magnitude;
			}
		}
		else if(es->parameters.get<unsigned int> ("inflow_profile") == 1)
		{
			if(es->parameters.get<bool> ("prescribed_flow"))
			{
				// not right but who fuckin cares..
				// better to have max v == 1
				//std::cout << "hi - " << direction / normalisation_constant * (pow(radius,2)-pow(r,2)) * time_scaling * flow_mag << std::endl;
				//std::cout << direction << std::endl;
				//std::cout << time_scaling << std::endl;
				//std::cout << flow_mag  << std::endl;
				//std::cout << normalisation_constant  << std::endl;
				//std::cout << (pow(radius,2)-pow(r,2))  << std::endl;
				return direction / normalisation_constant * (pow(radius,2)-pow(r,2)) * time_scaling * flow_mag;//*(pow(radius,2)-r*r);
			}
			else
			{
				return direction * (pow(radius,2)-pow(r,2))/pow(radius,2)*time_scaling * velocity_magnitude;
			}
		}
		else if(es->parameters.get<unsigned int> ("inflow_profile") == 2)
		{
			for(unsigned int i=0; i<time_values.size(); i++)
			{
				if(fabs(t-time_values[i]) < 1e-10)
				{
					//std::cout << "hey baby, time_values[i] = " << time_values[i] << std::endl;
					for(unsigned int j=0;j<position_values.size()-1; j++)
					{
						// replace p(1) with r
						if(r >= position_values[j] && r <= position_values[j+1] + 1e-4)
						{
							//std::cout << "r = " << r << std::endl;
							double distance = r - position_values[j];
							double axial_velocity = axial_velocity_values[i][j] + 
																				distance*(axial_velocity_values[i][j+1] - axial_velocity_values[i][j])
																				/(position_values[j+1]-position_values[j]);
							//std::cout << "velocity = " << axial_velocity << std::endl;
							return axial_velocity;			//i think we may need negative because outward normal

						}

						if(j==position_values.size()-2)
						{
							std::cout << "FAILURE WE HAVE GONE PAST THE END OF SPACE in prescribed flow... EXITING" << std::endl;
							std::cout << "r = " << r << std::endl;
							std::cout << "x(n-1) = " << position_values[j] << ", and x(n) = " << position_values[j+1] + 1e-4 << std::endl;
							std::exit(0);
							return 0;
						}
					}
					break;
				}
			}

			std::cout << "FAILURE WE HAVE GONE PAST THE END OF TIME in prescribed flow... EXITING" << std::endl;
			std::exit(0);
			return 0;
		
		}
		else
		{
			std::cout << "FAILURE: INFLOW PROFILE " << es->parameters.get<unsigned int> ("inflow_profile")  << " IS NOT SUPPORTED... EXITING" << std::endl;
			std::exit(0);			
			return 0;
		}

	};

  /**
   * Like before, but returns the values in a
   * writable reference.
	 * note that now that we have surfaces with arb normal we need a vector function
   */
  void operator() (const Point& p,
		   const Real time,
		   DenseVector<Number>& output)
	{
		//we want the thing to go in the negative normal direction
		double time_scaling = es->parameters.get<double>("time_scaling");
		double inflow_magnitude = calculate_velocity_magnitude(p,time,time_scaling);
		//double inflow_magnitude = (*this)(p,time);
		if(es->parameters.get<bool> ("increment_boundary_conditions"))
		{
			double previous_time_scaling = es->parameters.get<double>("previous_time_scaling");
			// always start from zero on first time step, unless IC is set
			if(es->parameters.get<unsigned int> ("t_step") == 1 && !(es->parameters.set<unsigned int> ("restart_time_step") > 0 || es->parameters.set<bool> ("stokes_ic")))
				previous_time_scaling = 0.;
			double previous_inflow_magnitude = calculate_velocity_magnitude(p,time,previous_time_scaling);

			// get the difference
			inflow_magnitude -= previous_inflow_magnitude;

		}

		double angle_of_inflow = es->parameters.get<double> ("angle_of_inflow");

		// use the object if defined otherwise defaults
		Point normal;
		if(surface_boundary_object == NULL)
		{
			if(es->parameters.get<bool> ("threed"))
				normal(2) = 1.0;	
			else
				normal(0) = 1.0;
		}
		else 
			normal = surface_boundary_object->get_normal();


		// potentially rotate the flow direction
		if(fabs(angle_of_inflow) > 1e-10)
		{
			Point old_normal;
			old_normal = normal;
			normal(0) = old_normal(0)*cos(angle_of_inflow) - old_normal(1)*sin(angle_of_inflow);
			normal(1) = old_normal(0)*sin(angle_of_inflow) + old_normal(1)*cos(angle_of_inflow);
		}


		// need a special case for womersley where we want to impose over whole volume initially
		if(!es->parameters.get<bool> ("prescribed_womersley"))
		{
			if(es->parameters.get<bool> ("threed"))
			{
				output(0) = inflow_magnitude * normal(0);
				output(1) = inflow_magnitude * normal(1);
				output(2) = inflow_magnitude * normal(2);
			}
			else
			{
				output(0) = inflow_magnitude * normal(0);
				output(1) = inflow_magnitude * normal(1);
			}
		}
		else
		{
			if(es->parameters.get<bool> ("threed"))
			{
				output(0) = 0.;
				output(1) = 0.;
				output(2) = inflow_magnitude;
			}
			else
			{
				output(0) = inflow_magnitude;
				output(1) = 0.;
			}
		}
	};

	void init_inflow_profile_values()
	{

		//****** read the values from file and put into data structure **********//
    std::ifstream infile;
    infile.open (es->parameters.get<std::string> ("prescribed_womersley_file").c_str());
		std::string temp;
	
		std::cout << "holla filename = " << es->parameters.get<std::string> ("prescribed_womersley_file") << std::endl;
		
		if(std::getline(infile, temp))
		{
			//std::cout << "k" << std::endl;
			std::istringstream buffer(temp);
			std::vector<double> line((std::istream_iterator<double>(buffer)),
			                         std::istream_iterator<double>());
			line.erase (line.begin());
			position_values = line;
		}

		while (std::getline(infile, temp)) {
			std::istringstream buffer(temp);
			std::vector<double> line((std::istream_iterator<double>(buffer)),
			                         std::istream_iterator<double>());

		//std::cout << "sigh, " << std::endl;
			time_values.push_back(line[0]);
			line.erase(line.begin());
			axial_velocity_values.push_back(line);
		}

	}

};

}


#endif
