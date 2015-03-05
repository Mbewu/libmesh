#ifndef MOVING_PIPE_BOUNDARY_FUNCTION_H
#define MOVING_PIPE_BOUNDARY_FUNCTION_H

#include "libmesh/function_base.h"
#include "libmesh/equation_systems.h"


namespace libMesh
{

template <typename Output=Number>
class MovingPipeBoundaryFunction : public FunctionBase<Output>
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
  MovingPipeBoundaryFunction(EquationSystems& _es)
			:  es(&_es)	{};

  /**
   * Returns a new deep copy of the function.
   */
  virtual AutoPtr<FunctionBase<Number> > clone () const
	{
		return AutoPtr<FunctionBase<Number> > (new MovingPipeBoundaryFunction (*es) );
	};

  /**
   * @returns the value at point \p p and time
   * \p time, which defaults to zero.
   */
  Number operator() (const Point&,// p,
		     const Real)// t=0.)
	{

		//const double x = p(0);
		//const double y = p(1);
		//const double z = p(2);

		return -es->parameters.get<double>("time_scaling");
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
		output(2) = (*this)(p,time);
	};

};

}


#endif
