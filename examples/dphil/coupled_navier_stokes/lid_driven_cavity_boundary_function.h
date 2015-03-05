#ifndef LID_DRIVEN_CAVITY_BOUNDARY_FUNCTION_H
#define LID_DRIVEN_CAVITY_BOUNDARY_FUNCTION_H

#include "libmesh/function_base.h"
#include "libmesh/equation_systems.h"


namespace libMesh
{

template <typename Output=Number>
class LidDrivenCavityBoundaryFunction : public FunctionBase<Output>
{
private:

  /**
   * The EquationSystems object that we're using here.
   */
  EquationSystems* es; 

public:

  /**
   * Constructor.
   */
  LidDrivenCavityBoundaryFunction(EquationSystems& _es)
			:  es(&_es)	{};

  /**
   * Returns a new deep copy of the function.
   */
  virtual AutoPtr<FunctionBase<Number> > clone () const
	{
		return AutoPtr<FunctionBase<Number> > (new LidDrivenCavityBoundaryFunction (*es) );
	};

  /**
   * @returns the value at point \p p and time
   * \p time, which defaults to zero.
   */
  Number operator() (const Point& p,
		     const Real)// t=0.)
	{
		bool threed = es->parameters.get<bool>("threed");

		const double x = p(0);
		//const double y = p(1);
		double z = 0; //p(2);

		if(threed)
			z = p(2);
	
		if(threed)
		{
			if(!es->parameters.get<bool>("flat_lid_driven_cavity"))
				return (1.0 - z) * z * (1 - x) * x / (0.5*0.5*0.5*0.5)*es->parameters.get<double>("time_scaling");
			else
				return es->parameters.get<double>("time_scaling");
		}
		else
		{
			if(!es->parameters.get<bool>("flat_lid_driven_cavity"))
			{
				if(p(1) > es->parameters.get<double>("cube_width")-1e-8)
					return (es->parameters.get<double>("cube_width") - x) * x / pow(0.5*es->parameters.get<double>("cube_width"),2.0)*es->parameters.get<double>("time_scaling");
				else
					return 0.;
			}	
			else
			{
				if(p(1) > es->parameters.get<double>("cube_width")-1e-8)
					return es->parameters.get<double>("time_scaling");				
				else
					return 0.;
			}
		}

	
	};

  /**
   * Like before, but returns the values in a
   * writable reference.
   */
  void operator() (const Point& p,
		   const Real time,
		   DenseVector<Number>& output)
	{
		// only components in the z direction
		output(0) = (*this)(p,time);
	};

};

}


#endif
