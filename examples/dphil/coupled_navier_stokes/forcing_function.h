#ifndef FORCING_FUNCTION_H
#define FORCING_FUNCTION_H

#include "libmesh/function_base.h"
#include "libmesh/equation_systems.h"

// want this to be able to handle square boundaries as well
// assume not centred on origin and just need "radius"

namespace libMesh
{

template <typename Output=Number>
class ForcingFunction : public FunctionBase<Output>
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
  ForcingFunction(EquationSystems& _es)
			:  es(&_es)	{};

  /**
   * Returns a new deep copy of the function.
   */
  virtual AutoPtr<FunctionBase<Number> > clone () const
	{
		return AutoPtr<FunctionBase<Number> > (new ForcingFunction (*es) );
	};

	/**
   * @returns need this function for compilation
   */
  Number operator() (const Point&,// p,
		     const Real)// t=0.)
	{
		std::cout << "forcing function is not a scalar... exiting." << std::endl;
		std::exit(0);
		
	}
  /**
   * Returns value of forcing function at point p
   */
  void operator() (const Point& p,
		   const Real,// time,
		   DenseVector<Number>& output)
	{
		if(es->parameters.get<bool> ("threed"))
		{
			output(0) = 0.;
			output(1) = 0.;
			output(2) = 0.;
		}
		else
		{
			
			/*
			if(es->parameters.get<bool> ("stokes"))
			{			
				output(0) = 0.;
				output(1) = 0.;
			}
			else
			{
				output(0) = + 2. * p(1) * p(0) * p(0);
				output(1) = + 2. * p(0) * p(1) * p(1);
			}
			*/

			/*
			if(es->parameters.get<bool> ("stokes"))
			{			
				output(0) = 0.;
				output(1) = 0.;
			}
			else
			{
				output(0) = + p(1);
				output(1) = + p(0);
			}
			*/

			if(es->parameters.get<unsigned int> ("problem_type") == 5)
			{
				double x = p(0);
				double y = p(1);
				//double lap_x		= pi*pi*(-(sin(pi*x) + cos(pi*y)));// - cos(pi*x*y)*(x*x + y*y));
				//double lap_y		= pi*pi*(-(sin(pi*y) + cos(pi*x)));// - cos(pi*x*y)*(x*x + y*y));

				double lap_x		= pi*pi*(-(sin(pi*x)));// + cos(pi*y)));// - cos(pi*x*y)*(x*x + y*y));
				double lap_y		= 0;//pi*pi*(-(cos(pi*x)));// - cos(pi*x*y)*(x*x + y*y));

				double pres_x		= pi*(-sin(pi*x));// + y*cos(pi*x*y));
				double pres_y		= pi*(cos(pi*y));// + x*cos(pi*x*y));
				pres_y		= 0.;// + x*cos(pi*x*y));

				double u	= 1 + sin(pi*x) + cos(pi*y);// + cos(pi*x*y);
				double v	= 1 + cos(pi*x) + sin(pi*y);// + cos(pi*x*y);


				double conv_x		= u*pi*(cos(pi*x) - y*sin(pi*x*y)) + v*pi*(-sin(pi*y) - x*sin(pi*x*y));
				double conv_y		= u*pi*(-sin(pi*x) - y*sin(pi*x*y)) + v*pi*(cos(pi*y) - x*sin(pi*x*y));


				if(es->parameters.get<bool> ("stokes"))
				{	

					output(0) = -lap_x + pres_x;
					output(1) = -lap_y + pres_y;
				}
				else
				{
					output(0) = conv_x - lap_x + pres_x;
					output(1) = conv_y - lap_y + pres_y;
				}
			}
			else
			{
				output(0) = 0.;
				output(1) = 0.;
			}
			
		}
	};

};

}


#endif
