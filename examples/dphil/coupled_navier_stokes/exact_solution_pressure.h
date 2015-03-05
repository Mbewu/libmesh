#ifndef EXACT_SOLUTION_PRESSURE_H
#define EXACT_SOLUTION_PRESSURE_H

#include "libmesh/function_base.h"
#include "libmesh/equation_systems.h"

// want this to be able to handle square boundaries as well
// assume not centred on origin and just need "radius"

namespace libMesh
{

template <typename Output=Number>
class ExactSolutionPressure : public FunctionBase<Output>
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
  ExactSolutionPressure(EquationSystems& _es)
			:  es(&_es)	{};

  /**
   * Returns a new deep copy of the function.
   */
  virtual AutoPtr<FunctionBase<Number> > clone () const
	{
		return AutoPtr<FunctionBase<Number> > (new ExactSolutionPressure (*es) );
	};

  /**
   * @returns the value at point \p p and time
   * \p time, which defaults to zero.
   */
  Number operator() (const Point& p,
		     const Real)		//const Real t=0.) t unused but inherited
	{

		if(es->parameters.get<bool> ("threed"))
		{
			std::cout << "exact solution does not work in 3D... exiting." << std::endl;
			std::exit(0);
			//output(0) = p(1);
			//output(1) = inflow_magnitude * normal(1);
			//output(2) = inflow_magnitude * normal(2);
		}
		else
		{
			//return 2*p(0) + 2*p(1);
			//return 0.;
			//return 1 + cos(pi*p(0)) + sin(pi*p(1));// + sin(pi*p(0)*p(1));
			return 1 + cos(pi*p(0));// + sin(pi*p(1));// + sin(pi*p(0)*p(1));
		}
	};

	/**
   * Like before, but returns the values in a
   * writable reference.
	 * note that now that we have surfaces with arb normal we need a vector function
   */
  void operator() (const Point&,	//const Point& p,
		   const Real,								//const Real time,
		   DenseVector<Number>&)			//DenseVector<Number>& output)
	{

			std::cout << "pressure exact solution is not a vector... exiting." << std::endl;
			std::exit(0);
	};


};

}


#endif
