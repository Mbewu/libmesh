/* The libMesh Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

 // <h1>Systems Example 2 - Unsteady Nonlinear Navier-Stokes</h1>
 //
 // This example shows how a simple, unsteady, nonlinear system of equations
 // can be solved in parallel.  The system of equations are the familiar
 // Navier-Stokes equations for low-speed incompressible fluid flow.  This
 // example introduces the concept of the inner nonlinear loop for each
 // timestep, and requires a good deal of linear algebra number-crunching
 // at each step.  If you have a ExodusII viewer such as ParaView installed,
 // the script movie.sh in this directory will also take appropriate screen
 // shots of each of the solution files in the time sequence.  These rgb files
 // can then be animated with the "animate" utility of ImageMagick if it is
 // installed on your system.  On a PIII 1GHz machine in debug mode, this
 // example takes a little over a minute to run.  If you would like to see
 // a more detailed time history, or compute more timesteps, that is certainly
 // possible by changing the n_timesteps and dt variables below.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/perf_log.h"
#include "libmesh/boundary_info.h"
#include "libmesh/utility.h"

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

//for parsing command line nicely
#include "libmesh/getpot.h"

// Defines the MeshData class, which allows you to store
// data about the mesh when reading in files, etc.
#include "libmesh/mesh_data.h"

#include "libmesh/gmv_io.h"
#include "libmesh/serial_mesh.h"


//for mesh generation
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_edge4.h"


#include "navier_stokes_assembler.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;


void build_1D_mesh (Mesh& mesh, bool zero_D, std::vector<std::vector<double> >& element_data);

// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

	GetPot command_line (argc, argv);

  // Trilinos solver NaNs by default on the zero pressure block.
  // We'll skip this example for now.
  if (libMesh::default_solver_package() == TRILINOS_SOLVERS)
    {
      std::cout << "We skip example 13 when using the Trilinos solvers.\n"
                << std::endl;
      return 0;
    }

	//no command line stuff
  std::cout << "Running " << argv[0];
		
  for (int i=1; i<argc; i++)
    std::cout << " " << argv[i];
  
	std::cout << std::endl << std::endl;
	// - JAMES always 3D for this one
	//int dim = 1;

  Mesh mesh(init.comm());
	bool zero_D = true;
  std::vector<std::vector<double> > element_data;
	build_1D_mesh(mesh,zero_D,element_data);


	//MeshTools::Generation::build_line(mesh,2,0.,1.,EDGE3);
  
  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  EquationSystems empty_equation_systems (mesh);

  // Declare the system and its variables.
  // Creates a transient system named "Navier-Stokes"
  TransientLinearImplicitSystem & system =
    equation_systems.add_system<TransientLinearImplicitSystem> ("Stokes_1D");
  TransientLinearImplicitSystem & radius_system =
    equation_systems.add_system<TransientLinearImplicitSystem> ("Radius");

	//will add flow later as expicit system  
  equation_systems.parameters.set<bool>("0D") = zero_D;

	//some more parameters
  equation_systems.parameters.set<double>("density") = 1.176e-9;
  equation_systems.parameters.set<double>("viscocity") = 15.23;
  equation_systems.parameters.set<double>("zeta_1") = -0.057;
  equation_systems.parameters.set<double>("zeta_2") = 0.2096;
  equation_systems.parameters.set<double>("zeta_3") = 0.00904;
  equation_systems.parameters.set<double>("E") = 3.3;


	if(zero_D)
	{
		system.add_variable ("p", FIRST,MONOMIAL);
		system.add_variable ("q", FIRST,MONOMIAL);
	}
	else
	{
		system.add_variable ("p", FIRST);
	}

	radius_system.add_variable("radius", CONSTANT, MONOMIAL);
	

  // Give the system a pointer to the matrix assembly
  // object.
	NavierStokesAssembler ns_assembler(equation_systems,element_data);
	system.attach_assemble_object (ns_assembler);

  // Initialize the data structures for the equation system.
  equation_systems.init ();

  // Prints information about the system to the screen.
  equation_systems.print_info();

  // Create a performance-logging object for this example
  PerfLog perf_log("Stokes 1D");

  // Get a reference to the Stokes system to use later.
  TransientLinearImplicitSystem&  stokes_1d_system =
        equation_systems.get_system<TransientLinearImplicitSystem>("Stokes_1D");

  // Now we begin the timestep loop to compute the time-accurate
  // solution of the equations.
  Real dt = 0.001;
  stokes_1d_system.time     = 0.0;
  unsigned int n_timesteps = 5000;

  // Write out every nth timestep to file.
  const unsigned int write_interval = 10;
  equation_systems.parameters.set<Real>("dt") = dt;

  // We also set a standard linear solver flag in the EquationSystems object
  // which controls the maxiumum number of linear solver iterations allowed.
	// JAMES: changed from 250 for large Re
  equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 1000;

	//JAMES: add a parameter to track whether a steady or unsteady problem
	const bool unsteady = true;
  equation_systems.parameters.set<bool> ("unsteady")   = unsteady;
	//change time stepping parameters
	if(!unsteady)
	{
		n_timesteps = 1;
		dt = 1.0;
	}

	//JAMES: add a K parameter
	const Real K = 0.02;
  equation_systems.parameters.set<Real> ("K")   = K;

  // Tell the system of equations what the timestep is by using
  // the set_parameter function.  The matrix assembly routine can
  // then reference this parameter.
  equation_systems.parameters.set<Real> ("dt")   = dt;



  // The first thing to do is to get a copy of the solution at
  // the current nonlinear iteration.  This value will be used to
  // determine if we can exit the nonlinear loop.

  for (unsigned int t_step=0; t_step<n_timesteps; ++t_step)
    {


      // Incremenet the time cfaceounter, set the time step size as
      // a parameter in the EquationSystem.
      stokes_1d_system.time += dt;

      // A pretty update message
      std::cout << "\n\n*** Solving time step " << t_step <<
                   ", time = " << stokes_1d_system.time <<
                   " ***" << std::endl;

      // Now we need to update the solution vector from the
      // previous time step.  This is done directly through
      // the reference to the Stokes system.
      *stokes_1d_system.old_local_solution = *stokes_1d_system.current_local_solution;

      // At the beginning of each solve, reset the linear solver tolerance
      // to a "reasonable" starting value.
      const Real initial_linear_solver_tol = 1.e-8;
      equation_systems.parameters.set<Real> ("linear solver tolerance") = initial_linear_solver_tol;

      // Assemble & solve the linear system.
      perf_log.push("linear solve");
      equation_systems.get_system("Stokes_1D").solve();
      perf_log.pop("linear solve");

			//stokes_1d_system.solution->print();

      // How many iterations were required to solve the linear system?
      const unsigned int n_linear_iterations = stokes_1d_system.n_linear_iterations();

      // What was the final residual of the linear system?
      const Real final_linear_residual = stokes_1d_system.final_linear_residual();

      // Print out convergence information for the linear and
      // nonlinear iterations.
      std::cout << "Linear solver converged at step: "
                << n_linear_iterations
                << ", final residual: "
                << final_linear_residual
                << std::endl;

			//add the radius variable
			//create an equation system to hold data like the radius of the element
			const DofMap& dof_map = radius_system.get_dof_map();

			MeshBase::const_element_iterator       el     =
				mesh.active_local_elements_begin();
			const MeshBase::const_element_iterator end_el =
				mesh.active_local_elements_end();
			std::vector<dof_id_type> dof_indices;
			const unsigned int rad_var = radius_system.variable_number ("radius");

			for ( ; el != end_el; ++el)
			{
				const Elem* elem = *el;

				dof_map.dof_indices (elem, dof_indices);

				const dof_id_type elem_id = elem->id();
				radius_system.solution->set(dof_indices[0], element_data[elem_id][6]);
				//std::cout << "radius = " << element_data[elem_id][6] << std::endl;
			}

			//JAMES: write in exodus or vtk formats
			const bool vtk = false;
#ifdef LIBMESH_HAVE_EXODUS_API
      if ((t_step+1)%write_interval == 0)
        {
          std::ostringstream file_name;


          //char buf[14];

					  //sprintf (buf, "out.%03d.gmv", t_step);

	          //GMVIO(mesh).write_equation_systems (buf,equation_systems);
					//GMVIO gmvio(mesh);

					//gmvio.discontinuous() = true;
					//gmvio.write_discontinuous_gmv (buf,equation_systems,true);
					//gmvio.write_discontinuous_gmv (buf,temp_es,true);
					//gmvio.write_equation_systems (buf,equation_systems); //doesn't work!!


					//changed to vtk file format for paraview reading correctly hmmm
					if(vtk)
					{
          	// We write the file in the vtk format.
	          file_name << "out_"
  	                  << std::setw(2)
  	                  << std::setfill('0')
  	                  << std::right
  	                  << t_step + 1
  										<< ".vtk";
						VTKIO vtkio(mesh);

            //vtkio.write_equation_systems (file_name.str(),
	          //                                    equation_systems);
						//vtkio.write_element_data(equation_systems);

					}
					else
					{

	          // We write the file in the ExodusII format.
	          file_name << "out.e-s."
	                    << std::setw(3)
	                    << std::setfill('0')
	                    << std::right
	                    << t_step + 1;

						
						ExodusII_IO exo(mesh);
	
//	          exo.write_timestep (file_name.str(),
//                                              empty_equation_systems,1,stokes_1d_system.time);

						exo.write_discontinuous_exodusII (file_name.str(),
                                              equation_systems);
						exo.write_time(1,stokes_1d_system.time);


			      //exo.write(file_name.str());
	          //exo.write_element_data(equation_systems);
	

					}
        }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API


    } // end timestep loop.

  // All done.
  return 0;
}

void build_1D_mesh (Mesh& mesh, bool zero_D, std::vector<std::vector<double> >& element_data)
{

	/*
	MeshTools::Generation::build_cube(mesh,
               2, 0, 0,
               0., 1.,
               0., 0.,
               0., 0.,
               EDGE3);

	*/

  // Declare that we are using the indexing utility routine
  // in the "Private" part of our current namespace.  If this doesn't
  // work in GCC 2.95.3 we can either remove it or stop supporting
  // 2.95.3 altogether.
  // Changing this to import the whole namespace... just importing idx
  // causes an internal compiler error for Intel Compiler 11.0 on Linux
  // in debug mode.
  //using namespace MeshTools::Generation::Private;

  // Clear the mesh and start from scratch
  mesh.clear();

  mesh.set_mesh_dimension(1);

	//create array of nodes
	//network must start from a single segment
	// this segment must be the first one listed
	// no loops

  unsigned int num_generations = 3;
	const double bifurction_angle = M_PI/4.0;
	const double length_ratio = 1.0;	//halves the length of segments in each generation
	double left_length_ratio = 0.876;	//halves the length of segments in each generation
	double right_length_ratio = 0.686;	//halves the length of segments in each generation

	left_length_ratio = length_ratio;
	right_length_ratio = length_ratio;

	Point p0,p1;
	p0 = Point (0.0,0.0,0);
	p1 = Point (0.0,0.0,20.3);	//chosen so that wall thickness is always positive

	std::vector<Point> vertices;
	vertices.push_back(p0);
	vertices.push_back(p1);

	//create array of segments
	const unsigned int vertices_per_cell = 2;		//will always be this unless curved or something
	std::vector<std::vector<unsigned int> > cell_vertices;
	std::vector<unsigned int> segment;
	segment.push_back(0);
	segment.push_back(1);
	cell_vertices.push_back(segment);

	// for each element/segment, assuming same numbering convention
	// 0 - parent elem no
  // 1 - daughter elem no 1
  // 2 - daughter elem no 2
  // 3 - sibling elem no
  // 4 - is daughter_1 bool
  // 5 - length
  // 6 - radius
	element_data.push_back(std::vector<double>(7,-1));
	element_data[0][0] = -1;
	element_data[0][5] = 20.3;
	element_data[0][6] = element_data[0][5]/5.8;


	unsigned int parent_start_node = 1;
	unsigned int parent_segment = 0;
  for(unsigned int i=1; i<num_generations; i++)
	{
		//for each parent we generate two new points and two new segments
		for(unsigned int j=0;j<pow(2,i-1); j++)
		{
			//we need to find the direction of this segment
			Point direction = vertices[cell_vertices[parent_segment + j][1]] 
												- vertices[cell_vertices[parent_segment + j][0]];
			Point unit_direction = direction.unit();

			Point p0_direction;
			p0_direction(0) = cos(bifurction_angle)*unit_direction(0) + sin(bifurction_angle)* unit_direction(2);
			p0_direction(2) = -sin(bifurction_angle)*unit_direction(0) + cos(bifurction_angle)* unit_direction(2);

			Point p1_direction;
			p1_direction(0) = cos(-bifurction_angle)*unit_direction(0) + sin(-bifurction_angle)* unit_direction(2);
			p1_direction(2) = -sin(-bifurction_angle)*unit_direction(0) + cos(-bifurction_angle)* unit_direction(2);

			double length_1 = element_data[parent_segment + j][5]*left_length_ratio;
			double length_2 = element_data[parent_segment + j][5]*right_length_ratio;

			//p0 = vertices[parent_start_node + j] + pow(length_ratio,i)*p0_direction;
			//p1 = vertices[parent_start_node + j] + pow(length_ratio,i)*p1_direction;
			p0 = vertices[parent_start_node + j] + length_1*p0_direction;
			p1 = vertices[parent_start_node + j] + length_2*p1_direction;

			vertices.push_back(p0);
			//now let us make these two segments
			segment[0] = parent_start_node + j;
			segment[1] = vertices.size() - 1;	//should be the last vertex added
			cell_vertices.push_back(segment);

			//update the data
			element_data[parent_segment + j][1] = cell_vertices.size() - 1;
			element_data.push_back(std::vector<double>(7,-1));
			element_data[cell_vertices.size() - 1][0] = parent_segment + j;
			element_data[cell_vertices.size() - 1][3] = cell_vertices.size();
			element_data[cell_vertices.size() - 1][4] = 1;
			element_data[cell_vertices.size() - 1][5] = element_data[parent_segment + j][5]*left_length_ratio;
			element_data[cell_vertices.size() - 1][6] = element_data[cell_vertices.size() - 1][5]/5.8;

			vertices.push_back(p1);
			segment[0] = parent_start_node + j;
			segment[1] = vertices.size() - 1;
			cell_vertices.push_back(segment);

			//update the data
			element_data[parent_segment + j][2] = cell_vertices.size() - 1;
			element_data.push_back(std::vector<double>(7,-1));
			element_data[cell_vertices.size() - 1][0] = parent_segment + j;
			element_data[cell_vertices.size() - 1][3] = cell_vertices.size() - 2;
			element_data[cell_vertices.size() - 1][4] = 0;
			element_data[cell_vertices.size() - 1][5] = element_data[parent_segment + j][5]*right_length_ratio;
			element_data[cell_vertices.size() - 1][6] = element_data[cell_vertices.size() - 1][5]/5.8;
		}
		
		parent_start_node += pow(2,i-1);
	  parent_segment 	 += pow(2,i-1);

	}

	
	const unsigned int
	n_vertices = vertices.size();
	unsigned int n_segments = cell_vertices.size();
	std::cout << "num segments = " << n_segments << std::endl;

	ElemType type;
	unsigned int nx;
	if(zero_D)
	{
		type = EDGE2;		//element type
		nx = 1;//2;	//no elements
	}
	else
	{
	
		type = EDGE2;		//element type
		nx = 2;//2;	//no elements
	}
	//hmmm we need to build a map from segment vertex to node number
	std::vector<unsigned int> segment_vertex_to_node(n_vertices,-1);
	std::vector<unsigned int> vertex_count(n_vertices,0);

	//need to construct a thing that counts up how many times a vertex was used
	for(unsigned int i=0;i<n_segments;i++)
	{
		vertex_count[cell_vertices[i][0]] += 1;
		vertex_count[cell_vertices[i][1]] += 1;
	}


	//keep track of the number of nodes from when the previous element added
	unsigned int node_start = 0;
	//loop over segments adding nodes and stuff

	//keep track of last node id added
	unsigned int node_id = 0;
	for(unsigned int i=0;i<n_segments;i++)
	{

		//starting point
		Point start_point = vertices[cell_vertices[i][0]];
		Point end_point = vertices[cell_vertices[i][1]];
		Point direction = end_point - start_point;

		libmesh_assert_not_equal_to (nx, 0);

		// Reserve elements
		switch (type)
		{
		//this will stay the same because for each segment we want the number of elements simply adds
		case INVALID_ELEM:
		case EDGE2:
		case EDGE3:
		  {
		    mesh.reserve_elem (nx);
		    break;
		  }

		default:
		  {
		    libMesh::err << "ERROR: Unrecognized 1D element type." << std::endl;
		    libmesh_error();
		  }
		}

		//this won't stay the same because some nodes may be shared by different segments, 
		// if no segments = ns then no nodes = nx+1 + (ns-1)*nx
		// may have to keep track of the direction of each segment so that after first seg
		// the first node is not added, but found somehow

		// Reserve nodes
		// don't really know how many nodes cause don't know how many bif, trif and endpoints
/*
		switch (type)
		{
		case INVALID_ELEM:
		case EDGE2:
		  {
				//if first segment then add nx+1 else add nx
				if(i==0)
			    mesh.reserve_nodes(nx+1);
				else
			    mesh.reserve_nodes(nx);
			
		    break;
		  }

		case EDGE3:
		  {
				if(i==0)
			    mesh.reserve_nodes(2*nx+1);
				else
			    mesh.reserve_nodes(2*nx);
		    break;
		  }


		default:
		  {
		    libMesh::err << "ERROR: Unrecognized 1D element type." << std::endl;
		    libmesh_error();
		  }
		}
*/

			
		// Build the nodes, depends on whether we're using linears,
		// quadratics or cubics and whether using uniform grid or Gauss-Lobatto
		// we will not consider gauss lobatto 
		switch(type)
		{
		  case INVALID_ELEM:
		  case EDGE2:
		    {
					//only for the first segment do we add the first node, otherwise from you know
					if(i==0)
					{
						mesh.add_point (start_point, node_id++);

						std::cout << "after added nodes n_nodes() = " << mesh.n_nodes() << std::endl;
					}

					//ensures that right number of elements are added 
		      for (unsigned int j=1; j<=nx; j++)
		      {
		        mesh.add_point (start_point + static_cast<Real>(j)/static_cast<Real>(nx) * direction, node_id++);

						std::cout << "after added nodes n_nodes() = " << mesh.n_nodes() << std::endl;
		      }
		      break;
		    }

		  case EDGE3:
		    {
					if(i==0)
					{
						mesh.add_point (start_point, node_id++);

						std::cout << "after added nodes n_nodes() = " << mesh.n_nodes() << std::endl;
					}
		      for (unsigned int j=1; j<=2*nx; j++)
		      {
		      
		        mesh.add_point (start_point + static_cast<Real>(j)/static_cast<Real>(2*nx) * direction, node_id++);
						std::cout << "after added nodes n_nodes() = " << mesh.n_nodes() << std::endl;
		      }
		      break;
		    }


		  default:
		    {
		      libMesh::err << "ERROR: Unrecognized 1D element type." << std::endl;
		      libmesh_error();
		    }

		}

		// Build the elements of the mesh
		switch(type)
		{
		  case INVALID_ELEM:
		  case EDGE2:
		    {
					//for the first segment we add all the elements normally
					//however for the others the first element needs to be connected to the 
					std::cout << "let's go" << std::endl;
					
		      for (unsigned int j=0; j<nx; j++)
					{
						
		        Elem* elem = mesh.add_elem (new Edge2);
						//if we are starting we need to create a new node otherwise not
						if(i==0)
		        	elem->set_node(0) = mesh.node_ptr(node_start + j);
						else
		        	elem->set_node(0) = mesh.node_ptr(segment_vertex_to_node[cell_vertices[i][0]]);

		        elem->set_node(1) = mesh.node_ptr(node_start + j +1);
	
						std::cout << "yep" << std::endl;

						//inflow boundary
						if(j==0 && vertex_count[cell_vertices[i][0]] == 1)
		        	mesh.boundary_info->add_side(elem, 0, 0);

						//sorta obvious??
						segment_vertex_to_node[0] = 0;

		        if (j == (nx-1))
						{
							//hmmm need to now know a priori which elements are terminal somehow
							//if the vertex count of the current segments end is equal to one then
							// terminal and add it to the boundary information
							if(vertex_count[cell_vertices[i][1]] == 1)
		          	mesh.boundary_info->add_side(elem, 1, 1);

							//also save the node of the end of the segment
							segment_vertex_to_node[cell_vertices[i][1]] = node_start + nx;
						}
		      }
		    break;
		    }

		  case EDGE3:
		    {

					//for the first segment we add all the elements normally
					//however for the others the first element needs to be connected to the 
					
		      for (unsigned int j=0; j<nx; j++)
		      {
						
		        Elem* elem = mesh.add_elem (new Edge3);

						if(i==0)
			        elem->set_node(0) = mesh.node_ptr(node_start + 2*j);
            else
              elem->set_node(0) = mesh.node_ptr(segment_vertex_to_node[cell_vertices[i][0]]);
		        elem->set_node(2) = mesh.node_ptr(node_start + 2*j + 1);
		        elem->set_node(1) = mesh.node_ptr(node_start + 2*j + 2);

						//inflow boundary
						if(j==0 && vertex_count[cell_vertices[i][0]] == 1)
		        	mesh.boundary_info->add_side(elem, 0, 0);


		        if (j == (nx-1))
						{
							//hmmm need to now know a priori which elements are terminal somehow
							//if the vertex count of the current segments end is equal to one then
							// terminal and add it to the boundary information
							if(vertex_count[cell_vertices[i][1]] == 1)
		          	mesh.boundary_info->add_side(elem, 1, 1);

							//also save the node of the end of the segment
							segment_vertex_to_node[cell_vertices[i][1]] = node_start + 2*nx;
						}

						//sorta obvious?? - the first one
						segment_vertex_to_node[0] = 0;

					}

		    break;
		    }

		  default:
		    {
		      libMesh::err << "ERROR: Unrecognized 1D element type." << std::endl;
		      libmesh_error();
		    }

		}

		// Scale the nodal positions
		// no need for scaling in new formulation
		//for (unsigned int p=0; p<mesh.n_nodes(); p++)
		//  mesh.node(p)(1) = (mesh.node(p)(1))*(xmax-xmin) + xmin;


		//update where to start adding nodes from
		node_start = mesh.n_nodes() - 1;
	}

	// Add sideset names to boundary info
	mesh.boundary_info->sideset_name(0) = "left";
	mesh.boundary_info->sideset_name(1) = "right";

	// Add nodeset names to boundary info
	mesh.boundary_info->nodeset_name(0) = "left";
	mesh.boundary_info->nodeset_name(1) = "right";

  // Done building the mesh.  Now prepare it for use.
  mesh.prepare_for_use (/*skip_renumber =*/ false);

}


