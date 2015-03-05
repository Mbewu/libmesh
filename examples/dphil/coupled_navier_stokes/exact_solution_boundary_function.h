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

public:

  /**
   * Constructor.
   */
  InflowBoundaryFunction(EquationSystems& _es, SurfaceBoundary* _surface_boundary_object)
			:  es(&_es), surface_boundary_object(_surface_boundary_object)	{};

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
		     const Real t=0.)
	{

		const double x = p(0);
		const double y = p(1);
		const double z = p(2);
	
		double r = sqrt(x*x + y*y);
		double radius = 0.5;//boundary_radius[boundary_id]
		double normalisation_constant = 1.0;

		// using new definition
		if(surface_boundary_object != NULL)
		{
			r = surface_boundary_object->get_normalised_distance_from_centroid(p);
			normalisation_constant = surface_boundary_object->get_unit_parabola_integral();
			radius = 1.0;
		}

		double direction = 1.0;
		if(es->parameters.get<bool> ("reverse_inflow_profile"))
			direction *= -1;
	
		//std::cout << "r = " << r << std::endl;
		//std::cout << "p = " << p << std::endl;
		// 0 - constant, 1 - parabolic (tensor)
		if(es->parameters.get<unsigned int> ("inflow_profile") == 0)
		{
			return direction*es->parameters.get<double>("time_scaling");
		}
		else if(es->parameters.get<unsigned int> ("inflow_profile") == 1)
		{
			bool prescribed_flow = false;
			if(prescribed_flow)
			{
				// not right but who fuckin cares..
				// better to have max v == 1
				return direction / normalisation_constant * (pow(radius,2)-pow(r,2)) * es->parameters.get<double>("time_scaling");//*(pow(radius,2)-r*r);
			}
			else
			{
				return direction * (pow(radius,2)-pow(r,2))*es->parameters.get<double>("time_scaling");
			}
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
		double inflow_magnitude = (*this)(p,time);


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
	};

};

}


#endif
