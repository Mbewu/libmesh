//class to flag elements based on the z position of the centre
// /centroid of the element. only refines for the moment.

// C++ include files that we may need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>

#include "libmesh/mesh_refinement.h"
#include "libmesh/libmesh.h"
#include "libmesh/elem.h"
#include "libmesh/mesh.h"

using namespace libMesh;

class XYZPositionElementFlagging : public MeshRefinement::ElementFlagging
{
	public:
	
		XYZPositionElementFlagging (MeshBase& mesh,double _x_upper, double _x_lower,
																double _y_upper = 0, double _y_lower = 0,double _z_upper = 0, double _z_lower = 0)
					:_mesh(mesh),x_upper(_x_upper),x_lower(_x_lower),y_upper(_y_upper),y_lower(_y_lower),z_upper(_z_upper),z_lower(_z_lower)
		{};

		virtual ~XYZPositionElementFlagging () {}

		virtual void flag_elements ()
		{
			// Need to clean up the refinement flags prior to this.

			// Let's do the element flagging
			MeshBase::element_iterator elem_it =
				_mesh.active_elements_begin();
			const MeshBase::element_iterator elem_end =
				_mesh.active_elements_end();

			for (; elem_it != elem_end; ++elem_it)
			{
				Elem* elem = *elem_it;

				// get the centroid of the element
				Point centroid = elem->centroid();
				if(elem->dim() == 2)
				{
					if (centroid(0) <= x_upper && centroid(0) >= x_lower &&
						centroid(1) <= y_upper && centroid(1) >= y_lower)
						elem->set_refinement_flag(Elem::REFINE);
				}
				else if(elem->dim() == 2)
				{
					if (centroid(0) <= x_upper && centroid(0) >= x_lower &&
						centroid(1) <= y_upper && centroid(1) >= y_lower &&
						centroid(2) <= z_upper && centroid(2) >= z_lower)
						elem->set_refinement_flag(Elem::REFINE);
				}
			}
		};

	private:

		MeshBase& _mesh;
		double x_upper;
		double x_lower;
		double y_upper;
		double y_lower;
		double z_upper;
		double z_lower;
};
