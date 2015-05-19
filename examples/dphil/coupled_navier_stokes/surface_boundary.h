#ifndef SURFACE_BOUNDARY_H
#define SURFACE_BOUNDARY_H

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"


// want this to be able to handle square boundaries as well
// assume not centred on origin and just need "radius"
// well I guess this can be given the mesh and stuff and then calculate for a given boundary id

namespace libMesh
{

class SurfaceBoundary
{
private:

  /**
   * The EquationSystems object that we're using here. probz don't need.
   */
  //EquationSystems* es;

	// list of points on edge of surface sorted.
	std::vector<Point>	boundary_points;
	std::vector<double>		boundary_point_angles;
	std::vector<double>		boundary_point_radii;
	Point centroid;
	Point normal;
	Point increasing_angle_direction; // this direction in which the angle increases
	int libmesh_geometry;
	bool threed;
	double parabolic_integral;
	double area;
	EquationSystems* es;
	

public:

  /**
   * Constructor. For a given boundary id this calculates the normal centroid and map from polygon to circle
   */
  SurfaceBoundary(EquationSystems& es_in): es(&es_in)	{}

	void init(Mesh& mesh,boundary_id_type surface_boundary_id, int _libmesh_geometry)
	{


		

		libmesh_geometry = _libmesh_geometry;		
		std::cout << "libmesh_geometry = " << libmesh_geometry << std::endl;

		//libmesh_geometry = 0;
		if(mesh.mesh_dimension() == 3)
			threed = true;
		else
			threed = false;



		if(!libmesh_geometry)
		{


			//EquationSystems es = _es;

			//first we loop over all the elements in the mesh, then all the ones in the boundary, 
			//loop over all the points on this surface and check if on the edge, if it is add it to boundary points
			// -- but check it hasn't been added already (actually we can do this later when we sort the points)

			// Iterators
			MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
			const MeshBase::const_element_iterator end_el = mesh.active_elements_end();


			normal.zero();
	//  FEType fe_vel_type = system->variable_type(u_var);
			const unsigned int dim = mesh.mesh_dimension();
			AutoPtr<FEBase> fe_face (FEBase::build(dim, FEType()));
			QGauss qrule (dim-1, static_cast<Order>(1));
			const std::vector<Point>& 	face_normals = fe_face->get_normals();
			fe_face->attach_quadrature_rule (&qrule);


			std::cout << "1" << std::endl;

			for ( ; el != end_el; ++el)
			{

			  const Elem* elem = *el;

				for (unsigned int s=0; s<elem->n_sides(); s++)
				{
					//for some reason it is natural to have more than one boundary id per side or even node
					std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

					if(boundary_ids.size() > 0 && boundary_ids[0] == surface_boundary_id) 
					{
						// now we know we are on a boundary, are we on the edge of a boundary
			      AutoPtr<Elem> side (elem->build_side(s));
						fe_face->reinit(elem, s);

			      // Loop over the nodes on the side.
			      for (unsigned int ns=0; ns<side->n_nodes(); ns++)
			      {
							std::vector<boundary_id_type> side_boundary_ids = mesh.boundary_info->boundary_ids(side->get_node(ns));

							// if node has multiple boundary ids then it is on the edge and we add it to the list.
							if(side_boundary_ids.size() > 1)
								boundary_points.push_back(*side->get_node(ns));
						}


						std::vector<Point> random_point(1,side->centroid());
						normal += side->volume() * face_normals[0];
					}
				}
			}

			// rescale normal
			normal = normal.unit();

			std::cout << "2" << std::endl;
/*
			for(unsigned int i=0; i< boundary_points.size(); i++)
			{
				std::cout << "boundary_point[" << i << "] = " << boundary_points[i] << std::endl;
			}
*/			
			std::cout << "normal = " << normal << std::endl;

			// sort the boundary points vector and remove duplicate boundary points
			double tol = 1e-10;
			std::vector<Point> sorted_boundary_points;
			sorted_boundary_points.push_back(boundary_points[0]);
			boundary_points.erase(boundary_points.begin() + 0);

			for(unsigned int j=0; j < sorted_boundary_points.size(); j++)
			{

				int min_index = -1;
				double min_value = 1e10;

				for(unsigned int i=0; i<boundary_points.size(); i++)
				{

					Point diff = sorted_boundary_points[j] - boundary_points[i];
					double distance = diff.size();

					// if this is the same point then yeah delete the point and move on
					if(distance < tol)
					{
						boundary_points.erase(boundary_points.begin() + i);
						i--;
					}
					else if(distance < min_value)
					{
						min_index = i;
						min_value = distance;
					}
				}
		
				// if there was even anything to look through
				if(boundary_points.size() > 0)
				{
					sorted_boundary_points.push_back(boundary_points[min_index]);
					boundary_points.erase(boundary_points.begin() + min_index);
				}
				else if(boundary_points.size() == 0)
				{
					break;
				}
			}


			std::cout << "3" << std::endl;
			boundary_points = sorted_boundary_points;

			/*
			for(unsigned int i=0; i< boundary_points.size(); i++)
			{
				std::cout << "boundary_point[" << i << "] = " << boundary_points[i] << std::endl;
			}
			*/

			if(threed)
			{
				// calculate the centroid using the wikipedia formula
				// we need to map to a 2d surface so choose an origin
				Point origin = boundary_points[5];
				//basis in surface
				Point basis_1 = boundary_points[0] - boundary_points[1];
				Point basis_2 = boundary_points[0] - boundary_points[2];
				//find vector orthogonal to vector_in_plane_1 that is component
				basis_2 = basis_2 - basis_2*basis_1	/(basis_1*basis_1) * basis_1;
	
				basis_1 = basis_1.unit();
				basis_2 = basis_2.unit();

				std::vector<Point> boundary_points_transformed(boundary_points.size());
				for(unsigned int i=0; i< boundary_points.size(); i++)
				{
					boundary_points_transformed[i](0) = (boundary_points[i](0) - origin(0)) * basis_1(0)
																							+ (boundary_points[i](1) - origin(1)) * basis_1(1)
																							+ (boundary_points[i](2) - origin(2)) * basis_1(2);

					boundary_points_transformed[i](1) = (boundary_points[i](0) - origin(0)) * basis_2(0)
																							+ (boundary_points[i](1) - origin(1)) * basis_2(1)
																							+ (boundary_points[i](2) - origin(2)) * basis_2(2);
				}	




				// calc area of convex polygon
				double area = 0.;
				for(unsigned int i=0; i<boundary_points_transformed.size() - 1; i++)
					area += (boundary_points_transformed[i](0)*boundary_points_transformed[i+1](1) - 
										boundary_points_transformed[i+1](0)*boundary_points_transformed[i](1));

				area += (boundary_points_transformed[boundary_points_transformed.size() - 1](0)*boundary_points_transformed[0](1) - 
										boundary_points_transformed[0](0)*boundary_points_transformed[boundary_points_transformed.size() - 1](1));

				area /= 2.;

				std::cout << "area = " << area << std::endl;


				centroid = 0.;
				Point centroid_2d;


		//		std::cout << "centroid_2d = " << centroid_2d << std::endl;
				for(unsigned int i=0; i<boundary_points_transformed.size() - 1; i++)
				{
					centroid_2d(0) += 	(boundary_points_transformed[i](0) + boundary_points_transformed[i+1](0)) *
													(boundary_points_transformed[i](0)*boundary_points_transformed[i+1](1) - 
														boundary_points_transformed[i+1](0)*boundary_points_transformed[i](1));

					centroid_2d(1) +=	(boundary_points_transformed[i](1) + boundary_points_transformed[i+1](1)) *
													(boundary_points_transformed[i](0)*boundary_points_transformed[i+1](1) - 
														boundary_points_transformed[i+1](0)*boundary_points_transformed[i](1));
				}

				centroid_2d(0) += (boundary_points_transformed[boundary_points_transformed.size() - 1](0) + boundary_points_transformed[0](0)) *
												(boundary_points_transformed[boundary_points_transformed.size() - 1](0)*boundary_points_transformed[0](1) - 
												boundary_points_transformed[0](0)*boundary_points_transformed[boundary_points_transformed.size() - 1](1));

				centroid_2d(1) += (boundary_points_transformed[boundary_points_transformed.size() - 1](1) + boundary_points_transformed[0](1)) *
												(boundary_points_transformed[boundary_points_transformed.size() - 1](0)*boundary_points_transformed[0](1) - 
												boundary_points_transformed[0](0)*boundary_points_transformed[boundary_points_transformed.size() - 1](1));

				centroid_2d /= (6. * area);

				centroid(0) = origin(0) + basis_1(0) * centroid_2d(0) + basis_2(0) * centroid_2d(1);
				centroid(1) = origin(1) + basis_1(1) * centroid_2d(0) + basis_2(1) * centroid_2d(1);
				centroid(2) = origin(2) + basis_1(2) * centroid_2d(0) + basis_2(2) * centroid_2d(1);


				std::cout << "4" << std::endl;

				// want the angle to increase from zero.
				// hmmm points will be approximately in a plane, but not exactly...
				// so can do the angle thing by, noticing when angles are decreasing and simply doing 360 - angle.
				// the question is, when i get on arb point "in the plane" how do i know if it is in the top bit or bottom bit
				// idea idea is that the angle should go in the same vague direction (dot product of directions should be positive)
				// calculate the angle between the point-centroid and the ref_point-centroid (ref_point is first point and arb)
				Point ref_point = boundary_points[0];
				boundary_point_angles.resize(boundary_points.size(),0.);
				bool bottom = false;
		
				// we can also calculate the increasing angle direction from the line 
				// orthogonal to the ref_point-centroid and first_point-centroid

				Point ref_line = boundary_points[0] - centroid;
				for(unsigned int i=0; i<boundary_points.size();i++)
				{
					Point line = boundary_points[i] - centroid;

					// first find the angle between the point-centroid and the ref_point-centroid
					double angle = acos(ref_line * line / (line.size() * ref_line.size()));//dot product formula.

					if(fabs(ref_line * line / (line.size() * ref_line.size()) - 1.0) < tol) // 0 degree line
						angle = 0.;
					else if(fabs(ref_line * line / (line.size() * ref_line.size()) + 1.0) < tol)	// 180 degree line
						angle = pi;

				
					if(bottom)
					{
						boundary_point_angles[i] = 2*pi - angle;
					}
					else if(i > 0 && angle < boundary_point_angles[i-1])
					{
						bottom = true;
						boundary_point_angles[i] = 2*pi - angle;
					}
					else
					{
						boundary_point_angles[i] = angle;
					}
				
				}

				std::cout << "5" << std::endl;

				// we calculate the direction in which the angle is increasing based 
				// on the perpendicular component, relative to the ref_point
				// assume there is at least "a few points" defining the boundary lol
				Point line = boundary_points[1] - centroid;
	
				// based on dot product cos rule and projection
				increasing_angle_direction = line - ref_line*line/(ref_line*ref_line) * ref_line;

	
				// calculate the distance to the centroid of each boundary point
				boundary_point_radii.resize(boundary_points.size(),0.);
				for(unsigned int i=0; i<boundary_points.size();i++)
				{
					Point diff = boundary_points[i] - centroid;
					boundary_point_radii[i] = diff.size();
				}


				/*
				for(unsigned int i=0; i< boundary_points.size(); i++)
				{
					std::cout << "point " << i << " = " << boundary_points[i] << std::endl;
					std::cout << "angle " << i << " = " << boundary_point_angles[i] << std::endl;
					std::cout << "radius " << i << " = " << boundary_point_radii[i] << std::endl;
				}
				*/
			}
			else
			{
				if(boundary_points.size() > 2)
				{
					std::cout << "error, in 2D a side should only have two boundary points.. Exiting" << std::endl;
					exit(0);
				}

				centroid = (boundary_points[0] + boundary_points[1])/2.0;
				area = (boundary_points[1] - boundary_points[0]).size();
			}

			std::cout << "centroid = " << centroid << std::endl;

		}
		else
		{

			// ***** CALCULATE THE NORMAL ****** //
			// we assume what it is for now... 
			// because calculating centroid could be overly tricky, plus we assume stuff later too
			//hmmm this is okay for now but may want to do both boundaries... this may be the issue
			if(threed)
				normal = Point(0.,0.,1.0);
			else
				normal = Point(1.0,0.0);



			// Iterators
			MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
			const MeshBase::const_element_iterator end_el = mesh.active_elements_end();


			normal.zero();
			//  FEType fe_vel_type = system->variable_type(u_var);
			const unsigned int dim = mesh.mesh_dimension();
			AutoPtr<FEBase> fe_face (FEBase::build(dim, FEType()));
			QGauss qrule (dim-1, static_cast<Order>(1));
			const std::vector<Point>& 	face_normals = fe_face->get_normals();
			fe_face->attach_quadrature_rule (&qrule);


			for ( ; el != end_el; ++el)
			{

			  const Elem* elem = *el;

				for (unsigned int s=0; s<elem->n_sides(); s++)
				{
					//for some reason it is natural to have more than one boundary id per side or even node
					std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

					if(boundary_ids.size() > 0 && boundary_ids[0] == surface_boundary_id) 
					{
						// now we know we are on a boundary, are we on the edge of a boundary
			      AutoPtr<Elem> side (elem->build_side(s));
						fe_face->reinit(elem, s);

						normal += side->volume() * face_normals[0];
					}
				}
			}

			// rescale normal
			normal = normal.unit();

			//std::cout << "normal where i know = " << normal << std::endl;


			// ***** CALCULATE THE CENTROID ****** //

			// can calculate the centroid by the weight sum of centroids of all the centroids

			el     = mesh.active_elements_begin();
			centroid.zero();

			double total_area = 0.;

			for ( ; el != end_el; ++el)
			{

			  const Elem* elem = *el;

				for (unsigned int s=0; s<elem->n_sides(); s++)
				{
					//for some reason it is natural to have more than one boundary id per side or even node
					std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

					if(boundary_ids.size() > 0 && boundary_ids[0] == surface_boundary_id) 
					{
						// now we know we are on a boundary, are we on the edge of a boundary
			      AutoPtr<Elem> side (elem->build_side(s));
						fe_face->reinit(elem, s);

						Point elem_centroid = side->centroid();
						centroid += side->volume() * elem_centroid;
						total_area += side->volume();
					}
				}
			}

			// rescale normal
			centroid /= total_area;

			//calculate

			//if we are using an axisymmetric geometry, then centroid is at y=0
			// but we don't use this geometry anymore so yeah
			if(es->parameters.get<unsigned int> ("geometry_type") == 5)
			{
				centroid(1) = 0.;
			}

			//std::cout << "centroid where i know = " << centroid << std::endl;

		}



		//want to calculate the integral of a parabolic profile on the surface to normalise parabolic integrals with
		// Iterators
		MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
		const MeshBase::const_element_iterator end_el = mesh.active_elements_end();


		//  FEType fe_vel_type = system->variable_type(u_var);
		const unsigned int dim = mesh.mesh_dimension();
		AutoPtr<FEBase> fe_face (FEBase::build(dim, FEType()));
		QGauss qrule (dim-1, static_cast<Order>(2));
		fe_face->attach_quadrature_rule (&qrule);
		//const std::vector<Point>& 	face_normals = fe_face->get_normals();
		const std::vector<Point>&  q_point = fe_face->get_xyz();
	  const std::vector<Real>&   JxW_face = fe_face->get_JxW();

		parabolic_integral = 0.;
		area = 0.;

			std::cout << "7" << std::endl;

		for ( ; el != end_el; ++el)
		{

			const Elem* elem = *el;

			for (unsigned int s=0; s<elem->n_sides(); s++)
			{
				//for some reason it is natural to have more than one boundary id per side or even node
				std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

				//std::cout << "boundary_ids[0] = " << boundary_ids[0] << std::endl;

				if(boundary_ids.size() > 0) 
				{
					if(boundary_ids[0] == surface_boundary_id) 
					{
						// now we know we are on a boundary integrate over the dofs
			      AutoPtr<Elem> side (elem->build_side(s));
						fe_face->reinit(elem, s);
	
						
						for (unsigned int qp=0; qp<qrule.n_points(); qp++)
						{
							double r = get_normalised_distance_from_centroid (q_point[qp]);
							double radius = 1.0;
							double parabola_value = (pow(radius,2)-pow(r,2));
							parabolic_integral += parabola_value * JxW_face[qp];
							area += JxW_face[qp];
						}
					}
				}
			}

			//if axisym then need to double the parabolic integral calculated
			// but not in use currently
			if(es->parameters.get<unsigned int> ("geometry_type") == 5)
				parabolic_integral *= 2.0;
		}

			std::cout << "8" << std::endl;
		//std::cout << "parabola integral = " << parabolic_integral << std::endl;
		//std::cout << "area = " << area << std::endl;

		//if(surface_boundary_id == 0)
		//	centroid = Point(0.5,0.5,1);
		//else if(surface_boundary_id == 1)
		//	centroid = Point(0.5,0.5,0);

	};




  /**
   * @returns the normalised distance to the centroid, normalised to the distance to the edge of the polygonal surface
	 * 					if we are on a libmesh generated object square or cuboid then 
   */
  Number get_normalised_distance_from_centroid (const Point& p)
	{
		if(!libmesh_geometry)
		{

			if(threed)
			{
				Point ref_line = boundary_points[0] - centroid;
				Point line = p - centroid;

				double tol = 1e-10;

				// first find the angle between the point-centroid and the ref_point-centroid
				double angle = acos(ref_line * line / (line.size() * ref_line.size()));//dot product formula.

				if(fabs(ref_line * line / (line.size() * ref_line.size()) - 1.0) < tol)	// for 0 degrees
					angle = 0.;
				else if(fabs(ref_line * line / (line.size() * ref_line.size()) + 1.0) < tol) // for 180 degrees
					angle = pi;

				// based on dot producct cos rule and projection
				Point orthogonal_component = line - ref_line*line/(ref_line*ref_line) * ref_line;

				// if opposite direction then for the angle we need 360 - angle
				if(orthogonal_component * increasing_angle_direction < 0)
					angle = 2*pi - angle;

				// next find between which points it is, could do a nice search but whatever
				unsigned int max_index = 0;
				while(angle >= boundary_point_angles[max_index] && max_index < boundary_point_angles.size() - 1)
					max_index++;

				double large_angle = boundary_point_angles[max_index] - boundary_point_angles[max_index - 1];
				double small_angle = angle - boundary_point_angles[max_index - 1];
				double a = boundary_point_radii[max_index - 1];
				double b = boundary_point_radii[max_index];
				double c = sqrt(a*a + b*b - 2*a*b*cos(large_angle));
				double opposite_angle = 0.;

				if(a >= b)
				{
					opposite_angle = asin(b/c*sin(large_angle));
				}
				else
				{
				
					opposite_angle = pi - asin(a/c*sin(large_angle)) - large_angle;
				}

				double correct_radius = a * sin(opposite_angle)/sin(pi - small_angle - opposite_angle);

				/*
				std::cout << "in radius calc:" << std::endl;
				std::cout << "angle = " << angle << std::endl;
				std::cout << "BIG angle = " << boundary_point_angles[max_index] << std::endl;
				std::cout << "SMALL angle = " << boundary_point_angles[max_index-1] << std::endl;
				std::cout << "BIG angle radii = " << boundary_point_radii[max_index] << std::endl;
				std::cout << "SMALL angle radii = " << boundary_point_radii[max_index-1] << std::endl;
				std::cout << "correct angle radii = " << correct_radius << std::endl;
				std::cout << "opp_angle = " << opposite_angle << std::endl;
				std::cout << "small_angle = " << small_angle << std::endl;
				std::cout << "other_small_angle = " << pi - small_angle - opposite_angle << std::endl;
				*/

				return line.size() / correct_radius;
			}
			else
			{

				double distance = (p - centroid).size();

				return distance/get_max_radius();		//radius is always 0.5 dudeman
			}
		}
		else
		{

			if(!threed)
			{
				double distance = (p - centroid).size();

				//if axisymmetric, radius is cube_width if not it is cube_width/2.0
				//but not in use
				if(es->parameters.get<unsigned int> ("geometry_type") == 5)
					distance = distance/(es->parameters.get<double>("cube_width"));
				else
				{
					//distance = distance/(es->parameters.get<double>("cube_width")/2.);
					distance = distance/get_max_radius();
				}
					
				return distance;		//radius is always 0.5 dudeman

			}
			else
			{
				double distance = (p - centroid).size();
				Point transformed_p = p - centroid;
				double angle = 0.;
				double correct_radius = 0.;
				// do separately for each region
				if(fabs(transformed_p(0)) < 1e-10 || fabs(transformed_p(1)) < 1e-10)
				{
					correct_radius = 0.5;
				}
				else if(transformed_p(0) * transformed_p(1) > 0)
				{
					angle = atan(transformed_p(1)/transformed_p(0));
					if(angle < M_PI/4.)
					{
						//correct_radius = 0.5 * sin(angle);
						correct_radius = 0.5 / cos(angle);
					}
					else
					{
						//correct_radius = 0.5 * sin(M_PI - angle);
						correct_radius = 0.5 / cos(M_PI/2. - angle);												
					}
				}
				else if(transformed_p(0)*transformed_p(1) < 0)
				{
					angle = atan(-transformed_p(1)/transformed_p(0));
					if(angle < M_PI/4.)
					{
						//correct_radius = 0.5 * sin(angle);
						correct_radius = 0.5 / cos(angle);
					}
					else
					{
						//correct_radius = 0.5 * sin(M_PI - angle);
						correct_radius = 0.5 / cos(M_PI/2. - angle);												
					}
				}
				else
				{
					std::cout << "HOUSTON WE HAVE A PROBLEM" << std::endl;
				}

				return distance/correct_radius;		//radius is always 0.5 dudeman
			}

		}
			
	};

	/**
   * @returns the normal
   */
  Point get_normal()
	{

		return normal;
	}

	/**
   * @returns the centroid
   */
  Point get_centroid()
	{

		return centroid;
	}

	/**
   * @returns the integrl of transformed unit parabola
   */
  double get_unit_parabola_integral()
	{

		return parabolic_integral;
	}

	/**
   * @returns the area of the surface
   */
  double get_area()
	{

		return area;
	}

	/**
   * @returns the area of the surface
   */
  double get_max_radius()
	{

		double max_radius = 0.;
		if(threed)
		{

			if(libmesh_geometry)
			{
				std::cout << "no support for function SurfaceBoundary::get_max_radius() for 3D libmesh geometries... exiting" << std::endl;
				std::exit(0);
			}

			for(unsigned int i=0; i<boundary_point_radii.size(); i++)
				if(boundary_point_radii[i] > max_radius)
					max_radius = boundary_point_radii[i];
		}
		else
		{
			
			if(es->parameters.get<unsigned int> ("geometry_type") == 5)
			{
				std::cout << "no support for function SurfaceBoundary::get_max_radius() for expanding pipe geometries... exiting" << std::endl;
				std::exit(0);
			}
			
			if(libmesh_geometry)
				max_radius = es->parameters.get<double>("cube_width")/2.0;
			else
				max_radius = area/2.0;
		}

		

		return max_radius;
	}

	/**
   * @returns true if the point is on the surface, assumes the surface is vaguely convex
	 *          essentially you can draw a line from all boundary points to the centroid (may be too strong in future)
   */
	bool is_on_surface(Point p)
	{

		//do polygon test, calc area of all triangles formed by connecting
		if(threed)
		{

			if(libmesh_geometry)
			{
				std::cout << "no support for function SurfaceBoundary::is_on_surface(Point p) for 3D libmesh geometries... exiting" << std::endl;
				std::exit(0);
			}

			double new_area = 0.;
			for(unsigned int i=0; i<boundary_points.size(); i++)
			{
				//triangle made from points: p, boundary_points[i] and boundary_points[i+1]
				Point p1 = p;
				Point p2 = boundary_points[i];
				Point p3;
				if(i+1 < boundary_points.size())
					p3 = boundary_points[i+1];
				else
					p3 = boundary_points[0];

				double a = (p1 - p2).size();
				double b = (p2 - p3).size();
				double c = (p3 - p1).size();
	
				double s = (a + b + c)/2.;
				new_area += sqrt(s*(s - a)*(s - b)*(s - c));
				
			}

			std::cout << "area = " << area << std::endl;
			std::cout << "new_area = " << new_area << std::endl;

			// have to give a fair bit of leeway here... hmmm
			if(new_area > area*1.005)
				return false;
			else
				return true;

		}
		else
		{

			if(es->parameters.get<unsigned int> ("geometry_type") == 5)
			{
				std::cout << "no support for function SurfaceBoundary::is_on_surface(Point p) for expanding pipe geometries... exiting" << std::endl;
				std::exit(0);
			}

			if(libmesh_geometry)
			{
				if((p - centroid).size() < (es->parameters.get<double>("cube_width")/2.0))
					return true;
				else
					return false;
			}
			else
			{
				if((p - centroid).size() < (area/2.0))
					return true;
				else
					return false;
			}
		}

	}
};

}


#endif
