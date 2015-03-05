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

class ZPositionElementFlagging : public MeshRefinement::ElementFlagging
{
	public:
	
		ZPositionElementFlagging (MeshBase& mesh,double _z_upper, double _z_lower)
					:_mesh(mesh),z_upper(_z_upper),z_lower(_z_lower)
		{};

		virtual ~ZPositionElementFlagging () {}

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
 
				if (centroid(2) <= z_upper && centroid(2) >= z_lower)
					elem->set_refinement_flag(Elem::REFINE);
			}
		};

	private:

		MeshBase& _mesh;
		double z_upper;
		double z_lower;
};
