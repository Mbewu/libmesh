#ifndef EXACT_SOLUTION_VELOCITY_H
#define EXACT_SOLUTION_VELOCITY_H

#include "libmesh/function_base.h"
#include "libmesh/equation_systems.h"

//#include <boost/math/special_functions/bessel.hpp>

// want this to be able to handle square boundaries as well
// assume not centred on origin and just need "radius"

namespace libMesh
{

template <typename Output=Number>
class ExactSolutionVelocity : public FunctionBase<Output>
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
  ExactSolutionVelocity(EquationSystems& _es)
			:  es(&_es)	{};

  /**
   * Returns a new deep copy of the function.
   */
  virtual AutoPtr<FunctionBase<Number> > clone () const
	{
		return AutoPtr<FunctionBase<Number> > (new ExactSolutionVelocity (*es) );
	};

  /**
   * @returns the value at point \p p and time
   * \p time, which defaults to zero.
   */
  Number operator() (const Point&,			//const Point& p,
		     const Real)										//const Real t=0.)
	{

		std::cout << "velocity exact solution is not a scalar... exiting." << std::endl;
		std::exit(0);
	};

  /**
   * Like before, but returns the values in a
   * writable reference.
	 * note that now that we have surfaces with arb normal we need a vector function
   */
  void operator() (const Point& p,
		   const Real,		//const Real time, //time unused but inherited
		   DenseVector<Number>& output)
	{
		std::cout << "hullo" << std::endl;

		
		if(es->parameters.get<bool> ("threed"))
		{
			std::cout << "exact solution does not work in 3D... exiting." << std::endl;
			//output(0) = p(1);
			//output(1) = inflow_magnitude * normal(1);
			//output(2) = inflow_magnitude * normal(2);
		}
		else
		{
			
			//output(0) = p(1)*p(1);
			//output(1) = p(0)*p(0);
			
			/*
			output(0) = p(1);
			output(1) = p(0);
			*/

			//output(0) = 1 + sin(pi*p(0)) + cos(pi*p(1));// + cos(pi*p(0)*p(1));
			//output(1) = 1 + cos(pi*p(0)) + sin(pi*p(1));// + cos(pi*p(0)*p(1));


			output(0) = 1 + sin(pi*p(0));// + cos(pi*p(1));// + cos(pi*p(0)*p(1));
			output(1) = 1 + cos(pi*p(0));// + sin(pi*p(1));// + cos(pi*p(0)*p(1));

			//std::cout << "yaya" << std::endl;
			//bessel function crap
			//double az = 1.;
			//double lim = 1.;
			//std::complex<float> cpxTerm = std::complex<float>(cos(az), -cos(sin(lim)));
			//std::complex<float> besselTerm = boost::math::cyl_bessel_j(0,cpxTerm);// cpxTerm); //boost can't do bessel functions of general complex variable

			//std::cout << "besselTerm = " << besselTerm << std::endl;
		}
	};

};

}


#endif
