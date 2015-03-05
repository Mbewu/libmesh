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

class BallBoundaryElementFlagging : public MeshRefinement::ElementFlagging
{
	public:
	
		BallBoundaryElementFlagging (MeshBase& mesh,Point _ball_centre, double _ball_radius)
					:_mesh(mesh),ball_centre(_ball_centre),ball_radius(_ball_radius)
		{};

		virtual ~BallBoundaryElementFlagging () {}

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
 				Point difference = centroid - ball_centre;
				bool elem_refined = false;
				
				//check if we are within the ball
				if (difference.size() < ball_radius)
				{
					//check if neighbouring cells are on boundary
					for(unsigned int side=0; side<elem->n_sides(); side++)
					{
						//check if we are within the ball
						//check if we are on the boundary
						if(elem->neighbor(side) == NULL)
						{
							elem->set_refinement_flag(Elem::REFINE);
							break;
						}
						
						Elem* elem_neighbor = elem->neighbor(side);
						Point neighbor_centroid = elem_neighbor->centroid();
 						Point neighbor_difference = neighbor_centroid - ball_centre;

						//check if neighbor within the ball
						if (neighbor_difference.size() < ball_radius)
						{
							//check if neighbor is on the boundary
							for(unsigned int neighbor_side=0; neighbor_side<elem_neighbor->n_sides(); neighbor_side++)
							{
								if(elem_neighbor->neighbor(neighbor_side) == NULL)
								{
									elem->set_refinement_flag(Elem::REFINE);
									elem_refined = true;
									break;
								}

								Elem* elem_neighbor_neighbor = elem_neighbor->neighbor(neighbor_side);
								
								Point neighbor_neighbor_centroid = elem_neighbor_neighbor->centroid();
 								Point neighbor_neighbor_difference = neighbor_neighbor_centroid - ball_centre;
						
								//check if neighbor within the ball
								if (neighbor_neighbor_difference.size() < ball_radius)
								{
									for(unsigned int neighbor_neighbor_side=0; neighbor_neighbor_side<elem_neighbor_neighbor->n_sides(); neighbor_neighbor_side++)
									{
										if(elem_neighbor_neighbor->neighbor(neighbor_neighbor_side) == NULL)
										{
											elem->set_refinement_flag(Elem::REFINE);
											elem_refined = true;
											break;
										}
									}// end neighbor neighbor loop
								}

								//if elem refined because neighbor of boundary
								if(elem_refined)
								{
									break;
								}
								

								

							}//end neighbor side loop

							//if elem refined because neighbor of boundary
							if(elem_refined)
							{
								break;
							}
						}

					}//end side loop
				}
			}
		};

	private:

		MeshBase& _mesh;
		Point ball_centre;
		double ball_radius;
};
