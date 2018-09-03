
#include "ns_assembler_3d.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

// boundary ids must have the inflow boundary id followed  by the outlfow 
// boundary ids in whatever order actually, must just match up
// these values are just magnitude values
void NSAssembler3D::init_bc (std::vector<unsigned int> boundary_ids,
											std::vector<double> pressure_values,
											std::vector<double> flow_values,
											std::vector<double> pressure_deriv_values,
											std::vector<double> previous_flow_values,
											std::vector<double> previous_previous_flow_values,
											std::vector<double> previous_pressure_values,
											std::vector<double> previous_previous_pressure_values,
											std::vector<double> previous_dynamic_pressure_values) 
{

	std::cout << "\nSetting 2D/3D BCs." << std::endl;
	

	double pressure_scaling = 1.0;
	if(es->parameters.get<bool>("nondimensionalised"))
	{
		pressure_scaling = es->parameters.get<double>("density") *
													pow(es->parameters.get<double>("velocity_scale"),2.0);
	}


	//resize so that doesn't break
	pressure_values.resize(boundary_ids.size());
	flow_values.resize(boundary_ids.size());
	pressure_deriv_values.resize(boundary_ids.size());
	previous_pressure_values.resize(boundary_ids.size());
	previous_previous_pressure_values.resize(boundary_ids.size());
	previous_dynamic_pressure_values.resize(boundary_ids.size());

	// what problem type we are considering
	// 0 - dirichlet at one inflow, pressure at others
	// 1 - pressure at all pressure, 
	// 2 - some kind of stress and 
	// 3 - is pipe moving sinusoidically in a static fluid
	// 4 - lid driven cavity, i.e. dirichlet everywhere
	// 5 - exact solution comparison, i.e. dirichlet everywhere
	const int problem_type    = es->parameters.get<unsigned int>("problem_type");
	
	// ********** dirichlet at inflow
	if(problem_type == 0 || problem_type == 3)
	{
		// set the bc type:
		// dirichlet at inflow (0) and pressure at others
		bc_type[boundary_ids[0]] = "dirichlet";
		for(unsigned int i = 1; i < boundary_ids.size(); i++)
			bc_type[boundary_ids[i]] = "pressure";
		
		
		// set the bc values:
		// dirichlet (0) handled elsewhere
		// pressure values set from input values, either straight pressure or calculate moghadam resistatnce
		bc_value[boundary_ids[0]] = 0;	//dirichlet will be handled outside this class	
		for(unsigned int i = 1; i < boundary_ids.size(); i++)
		{
			// calculate moghadam resistance
			if(es->parameters.get<unsigned int>("moghadam_coupling"))
			{
				// should always have some info here
				bc_value[boundary_ids[i]] = pressure_deriv_values[i];
				if(es->parameters.get<unsigned int>("moghadam_coupling") == 2)
				{
					bc_value_pressure_1d[boundary_ids[i]] = pressure_values[i];
					if(fabs(flow_values[i]) > 1e-10)
					{
						bc_value_resistance[boundary_ids[i]] = pressure_values[i]/flow_values[i];
						bc_value_flow_rate[boundary_ids[i]] = flow_values[i];
					}
					else
					{
						bc_value_resistance[boundary_ids[i]] = pressure_deriv_values[i];
						bc_value_flow_rate[boundary_ids[i]] = 0.;
					}

				}	
				/*
				if(fabs(flow_values[i]) > 1e-10)
				{
					bc_value[boundary_ids[i]] = pressure_values[i]/flow_values[i];	//these should ideally be zero...
					std::cout << "pressure[" << i << "] = " << pressure_values[i] << std::endl;
					std::cout << "flow_values[" << i << "] = " << flow_values[i] << std::endl;
					std::cout << "bc_value[" << i << "] = " << bc_value[boundary_ids[i]] << std::endl;
				}
				else if(fabs(pressure_deriv_values[i]) > 1e-10)
				{				
					std::cout << "using linear resistance" << std::endl;
				}
				else
				{
					bc_value[boundary_ids[i]] = 0;	//these should ideally be zero...					
					std::cout << "oops" << std::endl;
				}
				*/
			}
			// use straight pressure
			else
			{
				bc_value[boundary_ids[i]] = pressure_values[i];	//these should ideally be zero...
				std::cout << "pressure[" << i << "] = " << pressure_values[i] << std::endl;
			}
		}

		// calculate the interpolated flow value
		// if steady then this is just 1.		
		if(es->parameters.get<bool>("neumann_stabilised"))
		{
			if(es->parameters.get<bool>("neumann_stabilised_adjusted"))
			{
				for(unsigned int i = 0; i < boundary_ids.size(); i++)
				{
					interp_flow_bc_value[boundary_ids[i]] = previous_flow_values[i];	//these should ideally be zero...
					if(!es->parameters.get<unsigned int>("unsteady"))
					{
						interp_flow_bc_value[boundary_ids[i]] = -2.0/3.0;
						std::cout << "heya" << std::endl;
					}
				}
			}
			else if(es->parameters.get<bool>("neumann_stabilised_adjusted_interpolated"))
			{
				for(unsigned int i = 0; i < boundary_ids.size(); i++)
				{
					interp_flow_bc_value[boundary_ids[i]] = previous_flow_values[i] + 
						(previous_flow_values[i] - previous_previous_flow_values[i])/es->parameters.get<double>("previous_dt") * es->parameters.get<double>("dt");	
					if(!es->parameters.get<unsigned int>("unsteady"))
					{
						interp_flow_bc_value[boundary_ids[i]] = -2.0/3.0;
						std::cout << "heya2" << std::endl;
					}
				}
			}
				
		}
		

	}


	// ***************** pressure at all boundaries
	else if(problem_type == 1)
	{
	
		bc_type[boundary_ids[0]] = "pressure";
		for(unsigned int i = 1; i < boundary_ids.size(); i++)
			bc_type[boundary_ids[i]] = "pressure";
		
	
		// calculate the applied pressure from the flux that is read in
		if(es->parameters.get<bool> ("apply_pressure_from_flux"))
		{

			std::cout << "previous iteration applied pressure = " << bc_value[boundary_ids[0]] << std::endl;
			std::cout << "previous iteration calc pressure = " << pressure_values[0] << std::endl;
			std::cout << "previous time step calc pressure = " << previous_pressure_values[0] << std::endl;
			std::cout << "previous previous time step calc pressure = " << previous_previous_pressure_values[0] << std::endl;
			std::cout << "previous iteration flow rate = " << previous_flow_values[0] << std::endl;
			std::cout << "desired flow rate = " << flow_values[0] << std::endl;

			// set the internal pressure anf flow rates form previous nonlinear itertions
			pressure_prev_prev_nlin = pressure_prev_nlin;
			pressure_prev_nlin = pressure_values[0];
			flow_rate_prev_prev_nlin = flow_rate_prev_nlin;
			flow_rate_prev_nlin = previous_flow_values[0];



			double new_pressure = 0.;
			// if the first time step and nonlinear iteration then just use the input pressure in pressure_values[0]
			if(es->parameters.get<unsigned int>("t_step") == 1 && es->parameters.get<unsigned int>("nonlinear_iteration") == 1)
			{
				std::cout << "using input pressure at first time step" << std::endl;
				new_pressure = pressure_values[0];
			}
			// if the first nonlinear iteration, then guess the pressure from the previous time steps
			else if(es->parameters.get<unsigned int>("nonlinear_iteration") == 1 && es->parameters.get<bool>("interpolate_pressure_from_previous_steps"))
			{
				std::cout << "interpolating pressure from previous time step" << std::endl;
				// simple linear interpolation (assuming constant time step)
				// previously applied pressure + (previous - previous_previous)
				new_pressure = bc_value[boundary_ids[0]] + (previous_pressure_values[0] - previous_previous_pressure_values[0]);
			}
			// otherwise if the flow rate is not within an acceptable range, estimate the new pressure from the previous iteration
			else
			{

				// if the error of the flux is less than acceptable_flow_rate_error then we don't change anything
				// flow_values has the read in flow rates and previous_flow_values has the previous iteration flow rates
				//double flow_rate_error = fabs(previous_flow_values[0] - flow_values[0])/fabs(flow_values[0]);

				// NOTE: using an error metric percentage doesn't seem to work so well, let's use absolute.
				double flow_rate_error = fabs(previous_flow_values[0] - flow_values[0]);
				
				std::cout << "flow_rate_error = " << flow_rate_error << std::endl;
				if(flow_rate_error > es->parameters.get<double> ("acceptable_flow_rate_error"))
				{			
					std::cout << "changing applied pressure" << std::endl;
					// use resistance if resistance method or first couple of iterations of newton
					if(es->parameters.get<unsigned int> ("flow_rate_convergence_method") == 0
						|| (es->parameters.get<unsigned int> ("flow_rate_convergence_method") == 1 && es->parameters.get<unsigned int>("nonlinear_iteration") <= 2) )
					{
						std::cout << "using resitance method" << std::endl;
						// use estimate from resistance
						if(fabs(previous_flow_values[0]) > 1e-10)
						{
							// use the previously applied pressure
							double resistance = bc_value[boundary_ids[0]]/previous_flow_values[0];
							std::cout << "resistance = " << resistance << std::endl;
							new_pressure = resistance * flow_values[0];
						}
						else
						{
							std::cout << "flow rate is zero, so no resistance" << std::endl;
							std::cout << "EXITING..." << std::endl;
							std::exit(0);
						}
					}
					else if(es->parameters.get<unsigned int> ("flow_rate_convergence_method") == 1)
					{
						std::cout << "using newton's method" << std::endl;
						// newton's method
						// we use newton's method to estimate from previous nonlinear iterations
						// the resistance is the same as the derivative assuming linearity
						// P_new = P_prev_applied - (Q_prev_calc - Q_input)/((Q_prev_calc - Q_prev_prev_calc)/(P_prev_calc - P_prev_prev_calc))
						new_pressure = bc_value[boundary_ids[0]] - (flow_rate_prev_nlin - flow_values[0])/( (flow_rate_prev_nlin-flow_rate_prev_prev_nlin)/(pressure_prev_nlin-pressure_prev_prev_nlin) );
						
					}
					else
					{
						std::cout << "flow rate convergence method " << es->parameters.get<unsigned int> ("flow_rate_convergence_method") << " not defined." << std::endl;
						std::cout << "EXITING..." << std::endl;
						std::exit(0);
					}
				}
				else
				{
					std::cout << "not changing applied pressure" << std::endl;
					new_pressure = bc_value[boundary_ids[0]]; // keep it the same
				}
			}

			// once the new pressure boundary condition has been calculated, apply it

			bc_value[boundary_ids[0]] = new_pressure;//"0.00675896";	//"0.008";

			std::cout << "new applied pressure = " << new_pressure << std::endl;

		}
		else	// normal
		{

			if(es->parameters.get<unsigned int>("moghadam_coupling"))
			{
				// set the bc values:
				// pressure at the top (0) will be set before, potentially from file or input
				// pressure values set from input values, either straight pressure or calculate moghadam resistatnce
				// NOTE: bit of a hack putting this here but hey
				bc_value[boundary_ids[0]] = es->parameters.get<double>("pressure_mag_3d") * pressure_scaling * es->parameters.get<double>("time_scaling");

				for(unsigned int i = 1; i < boundary_ids.size(); i++)
				{
					// should always have some info here
					bc_value[boundary_ids[i]] = pressure_deriv_values[i];	

					if(es->parameters.get<unsigned int>("moghadam_coupling") == 2)
					{
						bc_value_pressure_1d[boundary_ids[i]] = pressure_values[i];
						if(fabs(flow_values[i]) > 1e-10)
						{
							bc_value_resistance[boundary_ids[i]] = pressure_values[i]/flow_values[i];
							bc_value_flow_rate[boundary_ids[i]] = flow_values[i];
						}
						else
						{
							bc_value_resistance[boundary_ids[i]] = pressure_deriv_values[i];
							bc_value_flow_rate[boundary_ids[i]] = 0.;
						}
					}
					// use straight pressure
					else
					{
						bc_value[boundary_ids[i]] = pressure_values[i];	//these won't matter but let's give them a value
					}
				}
			}
			else
			{
				bc_value[boundary_ids[0]] = es->parameters.get<double>("pressure_mag_3d") * pressure_scaling * es->parameters.get<double>("time_scaling");

				for(unsigned int i = 1; i < boundary_ids.size(); i++)
				{
					//bc_value[boundary_ids[i]] = pressure_values[i];	//these won't matter but let's give them a value
					bc_value[boundary_ids[i]] = es->parameters.get<double>("out_pressure_mag_3d") * pressure_scaling * es->parameters.get<double>("time_scaling");	//these won't matter but let's give them a value
				}
			}

		}

		// this should only be applied if we aren't doing a coupled simulation
		if(es->parameters.get<bool> ("vary_outlet_pressure"))
		{
			bc_value[boundary_ids[0]] = 0.;
			for(unsigned int i = 1; i < boundary_ids.size(); i++)
			{
				bc_value[boundary_ids[i]] = -pressure_values[0];	//these should ideally be zero...
			}
		}

		

		//here we calculate the interpolated flow value
		if(es->parameters.get<bool>("neumann_stabilised"))
		{
			if(es->parameters.get<bool>("neumann_stabilised_adjusted"))
			{
				for(unsigned int i = 0; i < boundary_ids.size(); i++)
				{
					interp_flow_bc_value[boundary_ids[i]] = previous_flow_values[i];	//these should ideally be zero...
				}
			}
			else if(es->parameters.get<bool>("neumann_stabilised_adjusted_interpolated"))
			{
				for(unsigned int i = 0; i < boundary_ids.size(); i++)
				{
					interp_flow_bc_value[boundary_ids[i]] = previous_flow_values[i] + 
						(previous_flow_values[i] - previous_previous_flow_values[i])/es->parameters.get<double>("previous_dt") * es->parameters.get<double>("dt");
				}

			}
				
		}
	}
	// stress at all boundaries
	else if(problem_type == 2)
	{
	
		bc_type[boundary_ids[0]] = "stress";
		for(unsigned int i = 1; i < boundary_ids.size(); i++)
			bc_type[boundary_ids[i]] = "stress";
		

	
		bc_value[boundary_ids[0]] = pressure_values[0];	//"0.008";
		for(unsigned int i = 1; i < boundary_ids.size(); i++)	
		{
			if(es->parameters.get<unsigned int>("moghadam_coupling"))
			{
				bc_value_pressure_1d[boundary_ids[i]] = pressure_values[i];
				// if no flow information then just use zero
				if(fabs(flow_values[i]) > 1e-10)
					bc_value[boundary_ids[i]] = pressure_values[i]/flow_values[i];	//these should ideally be zero...
				else
					bc_value[boundary_ids[i]] = 0;	//these should ideally be zero...					
			}
			else
				bc_value[boundary_ids[i]] = pressure_values[i];	//these should ideally be zero...
		}

		//here we calculate the interpolated flow value
		if(es->parameters.get<bool>("neumann_stabilised"))
		{
			if(es->parameters.get<bool>("neumann_stabilised_adjusted"))
			{
				for(unsigned int i = 0; i < boundary_ids.size(); i++)
					interp_flow_bc_value[boundary_ids[i]] = previous_flow_values[i];	//these should ideally be zero...
			}
			else if(es->parameters.get<bool>("neumann_stabilised_adjusted_interpolated"))
			{
				for(unsigned int i = 0; i < boundary_ids.size(); i++)
					interp_flow_bc_value[boundary_ids[i]] = previous_flow_values[i] + 
						(previous_flow_values[i] - previous_previous_flow_values[i])/es->parameters.get<double>("previous_dt") * es->parameters.get<double>("dt");	
			}
				
		}
	
	}
	// done setting boundary values
			

	

	// not sure what this is
	if(problem_type != 4)
	{
		//note that these values should actually be read in
		Real inflow_centre[] = {0.0,0.0,1.5};
		std::vector<Real> inflow_centre_vec (inflow_centre, inflow_centre + sizeof(inflow_centre) / sizeof(Real) );
		boundary_centre[boundary_ids[0]] = inflow_centre_vec;

		// set boundary radius at inflow, for more irregular geometries will need more fancy mappings etc
		boundary_radius[boundary_ids[0]] = 0.5;
	}


	// set up monolithic boundary nodes for use in assembly
	if(pressure_coupled)
	{
		std::cout << "Setting up coupling dofs for monolithic." << std::endl;
		find_1d_boundary_nodes();
	}


	if(es->parameters.get<unsigned int>("total_pressure_bc"))
	{
		for(unsigned int i = 0; i < boundary_ids.size(); i++)	
		{
			bc_value_dynamic_pressure[boundary_ids[i]] = previous_dynamic_pressure_values[i];		
		}
	}


	std::cout << "Done setting 2D/3D BCs." << std::endl;

		
}

void NSAssembler3D::find_1d_boundary_nodes()
{
  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es->get_mesh();

	// fix this crap, tak care of the boundary ids better
	primary_pressure_boundary_nodes_1d.resize(mesh.boundary_info->n_boundary_ids() - 2);
	secondary_pressure_boundary_nodes_1d.resize(mesh.boundary_info->n_boundary_ids() - 2);
	primary_flux_boundary_nodes_1d.resize(mesh.boundary_info->n_boundary_ids() - 2);
	secondary_flux_boundary_nodes_1d.resize(mesh.boundary_info->n_boundary_ids() - 2);

	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("ns1d");
	}

  const unsigned int p_var = system->variable_number ("P");
  const unsigned int q_var = system->variable_number ("Q");

  const DofMap & dof_map = system->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_p;
  std::vector<dof_id_type> dof_indices_q;

	// could do locally but don't really wanna
  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

	/*
	for(unsigned int i=0; i<airway_data.size(); i++)
	{
		std::cout << "airway " << i << " daughter 1 status = " << airway_data[i].get_is_daughter_1() << std::endl;
	}	
	*/

  for ( ; el != end_el; ++el)
  {
		const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
			//element data object starts numbering from 0
			//and the values referenced in it also do so need to take this into account

			//const int current_el_idx = elem->id();

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_q, q_var);

			// need to figure out whether we are on daughter 1 or daughter 2			
			const int current_el_idx = elem->id();
			//unsigned int current_1d_el_idx = current_el_idx -	n_initial_3d_elem;
			unsigned int current_1d_el_idx = elem_to_airway[current_el_idx];// -	n_initial_3d_elem;
			bool is_daughter_1 = airway_data[current_1d_el_idx].get_is_daughter_1();	//this is a bool duh!
			//std::cout << "elem id " << current_1d_el_idx << " daughter 1 status = " << is_daughter_1 << std::endl;
			

      //const unsigned int n_dofs   = 4;
      //const unsigned int n_p_dofs = 2;

			std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,0);
			if(boundary_ids.size() > 0)
			{
				if(boundary_ids[0] > 0 && boundary_ids[0] < 1000 && is_daughter_1)	// shouldn't get here if looking at side 0 anyway
				{
					primary_pressure_boundary_nodes_1d[boundary_ids[0]] = dof_indices_p[0];
					secondary_pressure_boundary_nodes_1d[boundary_ids[0]] = dof_indices_p[1];
					primary_flux_boundary_nodes_1d[boundary_ids[0]] = dof_indices_q[0];
					secondary_flux_boundary_nodes_1d[boundary_ids[0]] = dof_indices_q[1];

					//std::cout << "primary_flux_boundary_nodes_1d " << boundary_ids[0] << "  = " << primary_flux_boundary_nodes_1d[boundary_ids[0]] << std::endl;
				}
			}		
		}	
	}

	

}



// calculate the flux on a single boundary
double NSAssembler3D::calculate_flux(const int boundary_id)
{

	//std::cout << "begin calculating flux on bouundary " << boundary_id << std::endl;
	
	TransientLinearImplicitSystem * system;
	// Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}




	// clear the rhs before adding things in for flux
  system->rhs->zero ();
	const MeshBase& mesh = es->get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  const unsigned int u_var = system->variable_number ("u");
  const unsigned int v_var = system->variable_number ("v");
	//hmmm yeah you know const an shit
	unsigned int w_var = 0;
	if(threed)
		w_var = system->variable_number ("w");
  const unsigned int p_var = system->variable_number ("p");

  FEType fe_vel_type = system->variable_type(u_var);
	//boundary finite element stuffs
	AutoPtr<FEBase> fe_vel_face (FEBase::build(dim, fe_vel_type));

	QGauss qface(dim-1, fe_vel_type.default_quadrature_order());
	fe_vel_face->attach_quadrature_rule (&qface);


	//variables for equation system
  const std::vector<Real>& JxW_face = fe_vel_face->get_JxW();
  const std::vector<std::vector<Real> >& phi_face = fe_vel_face->get_phi();
  const std::vector<Point>& 						 qface_normals = fe_vel_face->get_normals();


  DenseVector<Number> Fe;
  DenseSubVector<Number>  Fu(Fe), Fv(Fe), Fw(Fe), Fp(Fe);

	// some dof stuff
  const DofMap & dof_map = system->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_w;
  std::vector<dof_id_type> dof_indices_p;

  //iterators
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	//for(unsigned int i=0; i<subdomains_3d.size(); i++)
	//	std::cout << "subdomains[" << i << "] = " << subdomains_3d[i] << std::endl;

	double radius_3d = 0.;

  for ( ; el != end_el; ++el)
  {				

	  const Elem* elem = *el;
		if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
		{

		  dof_map.dof_indices (elem, dof_indices);
		  dof_map.dof_indices (elem, dof_indices_u, u_var);
		  dof_map.dof_indices (elem, dof_indices_v, v_var);
			if(threed) { dof_map.dof_indices (elem, dof_indices_w, w_var); }
		  dof_map.dof_indices (elem, dof_indices_p, p_var);

		  unsigned int n_dofs   = dof_indices.size();
		  unsigned int n_u_dofs = dof_indices_u.size();
		  unsigned int n_v_dofs = dof_indices_v.size();
			unsigned int n_w_dofs = 0;
			if(threed) { n_w_dofs = dof_indices_w.size(); }
		  unsigned int n_p_dofs = dof_indices_p.size();

		  // Compute shape functions etc of current element
		  Fe.resize (n_dofs);

		  // Reposition the submatrices
		  Fu.reposition (u_var*n_u_dofs, n_u_dofs);
		  Fv.reposition (v_var*n_u_dofs, n_v_dofs);
			if(threed) { Fw.reposition (w_var*n_u_dofs, n_w_dofs); }
		  Fp.reposition (p_var*n_u_dofs, n_p_dofs);

		
		  // The following loops over the sides of the element.
		  // If the element has no neighbor on a side then that
		  // side MUST live on a boundary of the domain.

		
		  for (unsigned int s=0; s<elem->n_sides(); s++)
			{
		
				//for some reason it is natural to have more than one boundary id per side or even node
				std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

				if(boundary_ids.size() > 0) 
				{ 
					// dirichlet boundary conditions
					if(boundary_ids[0] == boundary_id)
					{	     
						if(es->parameters.get<bool> ("twod_oned_tree"))
						{
							// this should conserve the mean 3D velocity across the boundary
							//radius_3d = (*surface_boundaries)[boundary_id]->get_max_radius();
							radius_3d = (*surface_boundaries)[0]->get_max_radius();
						}
							

						fe_vel_face->reinit(elem, s);

						for (unsigned int qp=0; qp<qface.n_points(); qp++)
						{
						  for (unsigned int i=0; i<n_u_dofs; i++)
							{
						    Fu(i) += JxW_face[qp]*(phi_face[i][qp]*qface_normals[qp](0));
						    Fv(i) += JxW_face[qp]*(phi_face[i][qp]*qface_normals[qp](1));
						    if(threed) { Fw(i) += JxW_face[qp]*(phi_face[i][qp]*qface_normals[qp](2)); }
							}

						}//end face quad loop         
					}//end if(mesh.boundary_info->boundary_id(elem,s) == boundary_id)

				}
			}//end side loop

		  system->rhs->add_vector    (Fe, dof_indices);
		}
  } // end of element loop

  system->rhs->close();
	double flux = system->rhs->dot(*system->solution);

	
	if(es->parameters.get<bool> ("twod_oned_tree"))
	{
		std::cout << "radius_3d = " << radius_3d << std::endl;
		double flow_rate_scaling = M_PI*radius_3d/2.;
		flux *= flow_rate_scaling;
	}
	

  return flux;

}

// calculate the fluxes on multiple boundaries
void NSAssembler3D::setup_flow_rate_and_mean_pressure_vectors()
{

	std::cout << "Setting up flow rate and mean pressure vectors.. " << std::endl;
	
	TransientLinearImplicitSystem * system;
	// Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}

	std::cout << "boo" << std::endl;


	// clear the rhs before adding things in for flux
	const MeshBase& mesh = es->get_mesh();

  const unsigned int dim = mesh.mesh_dimension();
  const unsigned int num_1d_trees = es->parameters.set<unsigned int> ("num_1d_trees");
	const bool twod_oned_tree_pressure = es->parameters.get<bool> ("twod_oned_tree_pressure");

  const unsigned int u_var = system->variable_number ("u");
  const unsigned int v_var = system->variable_number ("v");
	//hmmm yeah you know const an shit
	unsigned int w_var = 0;
	if(threed)
		w_var = system->variable_number ("w");
  const unsigned int p_var = system->variable_number ("p");

  FEType fe_vel_type = system->variable_type(u_var);
	//boundary finite element stuffs
	AutoPtr<FEBase> fe_vel_face (FEBase::build(dim, fe_vel_type));

	FEType fe_pres_type = system->variable_type(p_var);
	//boundary finite element stuffs
	AutoPtr<FEBase> fe_pres_face (FEBase::build(dim, fe_pres_type));

	QGauss qface(dim-1, fe_vel_type.default_quadrature_order());
	fe_vel_face->attach_quadrature_rule (&qface);
	fe_pres_face->attach_quadrature_rule (&qface);


	//variables for equation system
  const std::vector<Real>& JxW_face = fe_vel_face->get_JxW();
  const std::vector<std::vector<Real> >& phi_face = fe_vel_face->get_phi();
  const std::vector<std::vector<Real> >& psi_face = fe_pres_face->get_phi();
  const std::vector<Point>& 						 qface_normals = fe_vel_face->get_normals();


  DenseVector<Number> Fe_Q, Fe_P;
  DenseSubVector<Number>  Fu(Fe_Q), Fv(Fe_Q), Fw(Fe_Q), Fp(Fe_P);

	// some dof stuff
  const DofMap & dof_map = system->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_w;
  std::vector<dof_id_type> dof_indices_p;

  std::vector<double> area(num_1d_trees+1);

  //iterators
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	//for(unsigned int i=0; i<subdomains_3d.size(); i++)
	//	std::cout << "subdomains[" << i << "] = " << subdomains_3d[i] << std::endl;
	std::cout << "boo" << std::endl;

  for ( ; el != end_el; ++el)
  {				

  	const Elem* elem = *el;
		if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
		{

		
				for (unsigned int s=0; s<elem->n_sides(); s++)
			{
		
				//for some reason it is natural to have more than one boundary id per side or even node
				std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

				if(boundary_ids.size() > 0) 
				{ 
					// check whether boundary is in list of boundaries we want to make flux vectors for
					if(boundary_ids[0] >= 0 && boundary_ids[0] <= num_1d_trees)
					{	     
						unsigned int boundary_id = boundary_ids[0];

						dof_map.dof_indices (elem, dof_indices);
						dof_map.dof_indices (elem, dof_indices_u, u_var);
						dof_map.dof_indices (elem, dof_indices_v, v_var);
						if(threed) { dof_map.dof_indices (elem, dof_indices_w, w_var); }
						dof_map.dof_indices (elem, dof_indices_p, p_var);

						unsigned int n_dofs   = dof_indices.size();
						unsigned int n_u_dofs = dof_indices_u.size();
						unsigned int n_v_dofs = dof_indices_v.size();
						unsigned int n_w_dofs = 0;
						if(threed) { n_w_dofs = dof_indices_w.size(); }
						unsigned int n_p_dofs = dof_indices_p.size();

						// Compute shape functions etc of current element
						Fe_Q.resize (n_dofs);
						Fe_P.resize (n_dofs);

						// Reposition the submatrices
						Fu.reposition (u_var*n_u_dofs, n_u_dofs);
						Fv.reposition (v_var*n_u_dofs, n_v_dofs);
						if(threed) { Fw.reposition (w_var*n_u_dofs, n_w_dofs); }
						Fp.reposition (p_var*n_u_dofs, n_p_dofs);


						// contributions to flow rate vector
						// note we assume only one side of the element on a boundary
						fe_vel_face->reinit(elem, s);

						// scale the flow rate accordingly
						double flow_rate_scaling = 1.0;
						if(es->parameters.get<bool> ("twod_oned_tree"))
						{
							// this should conserve the mean 3D velocity across the boundary
							//double radius_3d = (*surface_boundaries)[boundary_id]->get_max_radius();
							double radius_3d = (*surface_boundaries)[0]->get_max_radius();
							flow_rate_scaling = M_PI*radius_3d/2.;
						}

						for (unsigned int qp=0; qp<qface.n_points(); qp++)
						{
								for (unsigned int i=0; i<n_u_dofs; i++)
							{
							  		Fu(i) += flow_rate_scaling * JxW_face[qp]*(phi_face[i][qp]*qface_normals[qp](0));
							  		Fv(i) += flow_rate_scaling * JxW_face[qp]*(phi_face[i][qp]*qface_normals[qp](1));
							  		if(threed) { Fw(i) += JxW_face[qp]*(phi_face[i][qp]*qface_normals[qp](2)); }
							}

						}//end face quad loop

						// contributions to mean pressure vector
						fe_pres_face->reinit(elem, s);


						// need this to calculate the total area on this boundary
						AutoPtr<Elem> side (elem->build_side(s));
						area[boundary_id] += side->volume();


						// pressure

						// scale the flow rate accordingly
						double pressure_scaling = 1.0;
	
						if(twod_oned_tree_pressure)
						{
							pressure_scaling = 12./32.;
		
						}

						// scale 2D pressure to a 3D pressure
						for (unsigned int qp=0; qp<qface.n_points(); qp++)
						{
							for (unsigned int i=0; i<n_p_dofs; i++)
							{
								Fp(i) += JxW_face[qp]/pressure_scaling*(psi_face[i][qp]);
							}

						}//end face quad loop     

						std::ostringstream number;
						number << boundary_id;
							system->request_vector("Flow Rate Vector " + number.str())->add_vector    (Fe_Q, dof_indices);
							system->request_vector("Mean Pressure Vector " + number.str())->add_vector    (Fe_P, dof_indices);

						break;   
					}//end if(mesh.boundary_info->boundary_id(elem,s) == boundary_id)

				}
			}//end side loop

		}
  } // end of element loop

	std::cout << "boo" << std::endl;

  for(unsigned int i=0; i<=num_1d_trees; i++)
  {
	std::ostringstream number;
	number << i;
	system->request_vector("Flow Rate Vector " + number.str())->close();

	system->request_vector("Mean Pressure Vector " + number.str())->close();
	es->comm().sum(area[i]);
	system->request_vector("Mean Pressure Vector " + number.str())->scale(1./area[i]);
  }


}



// calculate the fluxes on multiple boundaries
std::vector<double> NSAssembler3D::calculate_fluxes(std::vector<unsigned int> boundary_ids)
{

	//std::cout << "begin calculating flux on bouundary " << boundary_id << std::endl;
	
	TransientLinearImplicitSystem * system;
	// Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}


	
	std::vector<double> fluxes(boundary_ids.size());

	for(unsigned int i=0; i<boundary_ids.size(); i++)
	{
		std::ostringstream number;
		number << boundary_ids[i];
		fluxes[i] = system->request_vector("Flow Rate Vector " + number.str())->dot(*system->solution);
	}


  	return fluxes;

}

// calculate the fluxes on multiple boundaries
std::vector<double> NSAssembler3D::calculate_pressures(std::vector<unsigned int> boundary_ids)
{

	//std::cout << "begin calculating flux on bouundary " << boundary_id << std::endl;
	
	TransientLinearImplicitSystem * system;
	// Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}


	
	std::vector<double> pressures(boundary_ids.size());

	for(unsigned int i=0; i<boundary_ids.size(); i++)
	{
		std::ostringstream number;
		number << boundary_ids[i];
		pressures[i] = system->request_vector("Mean Pressure Vector " + number.str())->dot(*system->solution);
	}


  	return pressures;

}



double NSAssembler3D::calculate_pressure(const int boundary_id)
{

	//std::cout << "begin calculating flux on bouundary " << boundary_id << std::endl;
	
	const MeshBase& mesh = es->get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}

	// clear the rhs before adding things in for flux
  system->rhs->zero ();

	const bool twod_oned_tree_pressure = es->parameters.get<bool> ("twod_oned_tree_pressure");

  const unsigned int u_var = system->variable_number ("u");
  const unsigned int v_var = system->variable_number ("v");
	//hmmm yeah you know const an shit
	unsigned int w_var = 0;
	if(threed) { w_var = system->variable_number ("w"); }
  const unsigned int p_var = system->variable_number ("p");


	FEType fe_pres_type = system->variable_type(p_var);
	//boundary finite element stuffs
	AutoPtr<FEBase> fe_pres_face (FEBase::build(dim, fe_pres_type));

	QGauss qface(dim-1, fe_pres_type.default_quadrature_order());
	fe_pres_face->attach_quadrature_rule (&qface);

	//variables for equation system
  const std::vector<Real>& JxW_face = fe_pres_face->get_JxW();
  const std::vector<std::vector<Real> >& psi_face = fe_pres_face->get_phi();


  DenseVector<Number> Fe;
  DenseSubVector<Number>  Fu(Fe), Fv(Fe), Fw(Fe), Fp(Fe);

	// some dof stuff
  const DofMap & dof_map = system->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_w;
  std::vector<dof_id_type> dof_indices_p;

  //iterators
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	double area = 0.0;

	// do in parallel?
  for ( ; el != end_el; ++el)
  {				

	  const Elem* elem = *el;
		if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
		{
			// note each element will only have one side on the boundary
		  for (unsigned int s=0; s<elem->n_sides(); s++)
			{


				//for some reason it is natural to have more than one boundary id per side or even node
				std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

				if(boundary_ids.size() > 0) 
				{ 
					// dirichlet boundary conditions
					if(boundary_ids[0] == boundary_id)
					{	     
						dof_map.dof_indices (elem, dof_indices);
						dof_map.dof_indices (elem, dof_indices_u, u_var);
						dof_map.dof_indices (elem, dof_indices_v, v_var);
						dof_map.dof_indices (elem, dof_indices_w, w_var);
						dof_map.dof_indices (elem, dof_indices_p, p_var);

						const unsigned int n_dofs   = dof_indices.size();
						const unsigned int n_u_dofs = dof_indices_u.size();
						const unsigned int n_v_dofs = dof_indices_v.size();
						unsigned int n_w_dofs = dof_indices_w.size();
						const unsigned int n_p_dofs = dof_indices_p.size();

						// Compute shape functions etc of current element
						Fe.resize (n_dofs);

						// Reposition the submatrices
						Fu.reposition (u_var*n_u_dofs, n_u_dofs);
						Fv.reposition (v_var*n_u_dofs, n_v_dofs);
						if(threed) { Fw.reposition (w_var*n_u_dofs, n_w_dofs); }
						Fp.reposition (p_var*n_u_dofs, n_p_dofs);
				
						// dirichlet boundary conditions
						// need this to calculate the total area on this boundary
						AutoPtr<Elem> side (elem->build_side(s));
						area += side->volume();
				
						fe_pres_face->reinit(elem, s);

						for (unsigned int qp=0; qp<qface.n_points(); qp++)
						{
							for (unsigned int i=0; i<n_p_dofs; i++)
							{
							  Fp(i) += JxW_face[qp]*(psi_face[i][qp]);
							}

						}//end face quad loop     

						system->rhs->add_vector    (Fe, dof_indices);    

					}//end if(mesh.boundary_info->boundary_id(elem,s) == boundary_id)
				}
			}//end side loop

		}
  } // end of element loop

	es->comm().sum(area);


  system->rhs->close();
	double pressure = system->rhs->dot(*system->solution) / area;
	
	// scale the flow rate accordingly
	double pressure_scaling = 1.0;
	
	if(twod_oned_tree_pressure)
	{
		pressure_scaling = 12./32.;		
	}

	// take from 2D to 3D
	pressure /= pressure_scaling;

  return pressure;

}



double NSAssembler3D::calculate_previous_dynamic_pressure(const int boundary_id)
{


	std::cout << "begin calculating dyanamic pressure " << boundary_id << std::endl;
	
	const MeshBase& mesh = es->get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}

	// Dimensional calculation fixes
  const bool nondimensionalised = es->parameters.get < bool >("nondimensionalised");
  const double rho = es->parameters.get < double >("density");
	double density_scale = 1.0;
	if(es->parameters.get < bool > ("dimensional_calculation_fix"))
	{
		// constant for the density, to go on the unsteady and convection terms
		// - if reynolds number calculation, this is 1.0
		// - if dimensional calculation this is \rho
		if(!nondimensionalised)
			density_scale = rho;
	}


  const unsigned int u_var = system->variable_number ("u");
  const unsigned int v_var = system->variable_number ("v");
	//hmmm yeah you know const an shit
	unsigned int w_var = 0;
	if(threed) { w_var = system->variable_number ("w"); }


	FEType fe_vel_type = system->variable_type(u_var);
	//boundary finite element stuffs
	AutoPtr<FEBase> fe_vel_face (FEBase::build(dim, fe_vel_type));

	QGauss qface(dim-1, fe_vel_type.default_quadrature_order());
	fe_vel_face->attach_quadrature_rule (&qface);

	//variables for equation system
  const std::vector<Real>& JxW_face = fe_vel_face->get_JxW();
  const std::vector<std::vector<Real> >& phi_face = fe_vel_face->get_phi();

	// some dof stuff
  const DofMap & dof_map = system->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_w;

  //iterators
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	double area = 0.0;
	double dynamic_pressure = 0.0;

	// do in parallel?
  for ( ; el != end_el; ++el)
  {				

	  const Elem* elem = *el;
		if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
		{
			// note each element will only have one side on the boundary
		  for (unsigned int s=0; s<elem->n_sides(); s++)
			{


				//for some reason it is natural to have more than one boundary id per side or even node
				std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

				if(boundary_ids.size() > 0) 
				{ 
					// dirichlet boundary conditions
					if(boundary_ids[0] == boundary_id)
					{	     
						dof_map.dof_indices (elem, dof_indices);
						dof_map.dof_indices (elem, dof_indices_u, u_var);
						dof_map.dof_indices (elem, dof_indices_v, v_var);
						dof_map.dof_indices (elem, dof_indices_w, w_var);

						const unsigned int n_dofs   = dof_indices.size();
						const unsigned int n_u_dofs = dof_indices_u.size();

				
						// dirichlet boundary conditions
						// need this to calculate the total area on this boundary
						AutoPtr<Elem> side (elem->build_side(s));
						area += side->volume();

						fe_vel_face->reinit(elem, s);

						for (unsigned int qp=0; qp<qface.n_points(); qp++)
						{
							double u=0.,v=0.,w=0.;

							for (unsigned int l = 0; l < n_u_dofs; l++)
			 				{
								// From the previous Newton iterate:
								u += phi_face[l][qp] * system->current_solution (dof_indices_u[l]);
								v += phi_face[l][qp] * system->current_solution (dof_indices_v[l]);
								if (threed)
								{
									w += phi_face[l][qp] * system->current_solution (dof_indices_w[l]);
								}				
							}

							double velocity_mag_squared = u*u + v*v + w*w;

						  dynamic_pressure += JxW_face[qp]*0.5*density_scale*velocity_mag_squared;

						}//end face quad loop     
  

					}//end if(mesh.boundary_info->boundary_id(elem,s) == boundary_id)
				}
			}//end side loop

		}
  } // end of element loop

	es->comm().sum(area);
	es->comm().sum(dynamic_pressure);

	dynamic_pressure /= area;

  return dynamic_pressure;

}



// calculate the mass avergaed total pressure
//
// = \int {P_tot * (velocity \dot normal) dS} / \int {(velocity \dot normal)dS}
// P_tot = p + 0.5*density*velocity_mag_squared
double NSAssembler3D::calculate_mass_averaged_total_pressure(const int boundary_id)
{


	std::cout << "begin calculating total pressure " << boundary_id << std::endl;
	
	const MeshBase& mesh = es->get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}

	// Dimensional calculation fixes
	const bool twod_oned_tree_pressure = es->parameters.get<bool> ("twod_oned_tree_pressure");
  const bool nondimensionalised = es->parameters.get < bool >("nondimensionalised");
  const double rho = es->parameters.get < double >("density");
	double density_scale = 1.0;
	if(es->parameters.get < bool > ("dimensional_calculation_fix"))
	{
		// constant for the density, to go on the unsteady and convection terms
		// - if reynolds number calculation, this is 1.0
		// - if dimensional calculation this is \rho
		if(!nondimensionalised)
			density_scale = rho;
	}


  const unsigned int u_var = system->variable_number ("u");
  const unsigned int v_var = system->variable_number ("v");
	//hmmm yeah you know const an shit
	unsigned int w_var = 0;
	if(threed) { w_var = system->variable_number ("w"); }
  const unsigned int p_var = system->variable_number ("p");


	FEType fe_vel_type = system->variable_type(u_var);
	//boundary finite element stuffs
	AutoPtr<FEBase> fe_vel_face (FEBase::build(dim, fe_vel_type));

	QGauss qface(dim-1, fe_vel_type.default_quadrature_order());
	fe_vel_face->attach_quadrature_rule (&qface);


	FEType fe_pres_type = system->variable_type(p_var);
	//boundary finite element stuffs
	AutoPtr<FEBase> fe_pres_face (FEBase::build(dim, fe_pres_type));
	fe_pres_face->attach_quadrature_rule (&qface);

	//variables for equation system
  const std::vector<Real>& JxW_face = fe_vel_face->get_JxW();
  const std::vector<std::vector<Real> >& phi_face = fe_vel_face->get_phi();
  const std::vector < Point > &qface_normals = fe_vel_face->get_normals ();
  const std::vector<std::vector<Real> >& psi_face = fe_pres_face->get_phi();

	// some dof stuff
  const DofMap & dof_map = system->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_w;
  std::vector<dof_id_type> dof_indices_p;

  //iterators
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	double numerator = 0.0;
	double denominator = 0.0;

	// do in parallel?
  for ( ; el != end_el; ++el)
  {				

	  const Elem* elem = *el;
		if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
		{
			// note each element will only have one side on the boundary
		  for (unsigned int s=0; s<elem->n_sides(); s++)
			{


				//for some reason it is natural to have more than one boundary id per side or even node
				std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

				if(boundary_ids.size() > 0) 
				{ 
					// dirichlet boundary conditions
					if(boundary_ids[0] == boundary_id)
					{	     
						dof_map.dof_indices (elem, dof_indices);
						dof_map.dof_indices (elem, dof_indices_u, u_var);
						dof_map.dof_indices (elem, dof_indices_v, v_var);
						dof_map.dof_indices (elem, dof_indices_w, w_var);
						dof_map.dof_indices (elem, dof_indices_p, p_var);

						const unsigned int n_dofs   = dof_indices.size();
						const unsigned int n_u_dofs = dof_indices_u.size();
						const unsigned int n_p_dofs = dof_indices_p.size();

				
						// dirichlet boundary conditions
						// need this to calculate the total area on this boundary
						AutoPtr<Elem> side (elem->build_side(s));

						fe_vel_face->reinit(elem, s);
						fe_pres_face->reinit(elem, s);

						for (unsigned int qp=0; qp<qface.n_points(); qp++)
						{
							double u=0.,v=0.,w=0.,p=0.;

							for (unsigned int l = 0; l < n_u_dofs; l++)
			 				{
								// From the previous Newton iterate:
								u += phi_face[l][qp] * system->current_solution (dof_indices_u[l]);
								v += phi_face[l][qp] * system->current_solution (dof_indices_v[l]);
								if (threed)
								{
									w += phi_face[l][qp] * system->current_solution (dof_indices_w[l]);
								}				
							}

							NumberVectorValue U;
							if (threed)
								U = NumberVectorValue (u, v, w);
							else
								U = NumberVectorValue (u, v);


							for (unsigned int l = 0; l < n_p_dofs; l++)
			 				{
								// From the previous Newton iterate:
								p += psi_face[l][qp] * system->current_solution (dof_indices_p[l]);
							}

							double velocity_mag_squared = u*u + v*v + w*w;

						  numerator += JxW_face[qp] * (p + 0.5*density_scale*velocity_mag_squared) * (U * qface_normals[qp]);
						  denominator += JxW_face[qp] * (U * qface_normals[qp]);

						}//end face quad loop     
  

					}//end if(mesh.boundary_info->boundary_id(elem,s) == boundary_id)
				}
			}//end side loop

		}
  } // end of element loop

	es->comm().sum(numerator);
	es->comm().sum(denominator);

	double mass_averaged_total_pressure = numerator/denominator;

	double pressure_scaling = 1.0;
	if(twod_oned_tree_pressure)
	{
		pressure_scaling = 12./32.;		
	}

	std::cout << "mass averaged total pressure 2D = " << mass_averaged_total_pressure << std::endl;
	std::cout << "mass averaged total pressure 3D = " << mass_averaged_total_pressure/pressure_scaling << std::endl;

  return mass_averaged_total_pressure/pressure_scaling;

}



void NSAssembler3D::estimate_error(ErrorVector& error_vector)
{
	const MeshBase& mesh = es->get_mesh();

	// resize and set to zero
	error_vector.resize(mesh.max_elem_id());
  std::fill (error_vector.begin(), error_vector.end(), 0.);

	// tell assemble function to fill the error vector with values
	// should fill error_vector with relative sizes of convection and diffusion
	estimating_error = true;

	if(es->parameters.get<unsigned int> ("error_estimator")  == 1 ||
			es->parameters.get<unsigned int> ("error_estimator") == 2)
		assemble(error_vector);
	else if(es->parameters.get<unsigned int> ("error_estimator") == 3)
		calculate_peclet_number(error_vector);

	estimating_error = false;
}

void NSAssembler3D::calculate_peclet_number(ErrorVector& error_vector)
{

	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}

	// Problem parameters
  Real reynolds_number    = es->parameters.get<Real>("reynolds_number");


	// Mesh and system things
	const MeshBase& mesh = es->get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  const unsigned int u_var = system->variable_number ("u");
  const unsigned int v_var = system->variable_number ("v");
	unsigned int w_var = 0;
	if(threed)
		w_var = system->variable_number ("w");


	// Finite element types
  FEType fe_vel_type = system->variable_type(u_var);

  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));


	// Quadrature rules
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  fe_vel->attach_quadrature_rule (&qrule);


	// Initialise some variables
  //const std::vector<Real>& JxW = fe_vel->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();

	// DofMap things
  const DofMap & dof_map = system->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_w;


  // Iterators
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
		
  for ( ; el != end_el; ++el)
  {				
    const Elem* elem = *el;
		// only solve on the 3d subdomain
		if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
		{

		  dof_map.dof_indices (elem, dof_indices);
		  dof_map.dof_indices (elem, dof_indices_u, u_var);
		  dof_map.dof_indices (elem, dof_indices_v, v_var);
			if(threed)
				dof_map.dof_indices (elem, dof_indices_w, w_var);

		  //unsigned int n_dofs   = dof_indices.size();
		  unsigned int n_u_dofs = dof_indices_u.size();

		  // Compute shape functions etc of current element
		  fe_vel->reinit  (elem);
			Real elem_volume = elem->volume();
			//take cube root to get diameter of the element
			double h_T = 0;
			if(threed)
				h_T = pow(elem_volume,1.0/3.0);
			else
				h_T = pow(elem_volume,1.0/2.0);


	    // Values to hold the solution & its gradient at the previous timestep.
	    Number   mean_u = 0., mean_v = 0., mean_w = 0.;

		  // Build volume contribution to element matrix and rhs
		  for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		  {
		    // Values to hold the solution & its gradient at the previous timestep.
		    Number   u = 0., v = 0., w = 0.;

		    // Compute the velocity & its gradient from the previous timestep
		    // and the old Newton iterate.
		    for (unsigned int l=0; l<n_u_dofs; l++)
		    {
		      // From the previous Newton iterate:
		      u += phi[l][qp]*system->current_solution (dof_indices_u[l]);
		      v += phi[l][qp]*system->current_solution (dof_indices_v[l]);
					if(threed)
			      w += phi[l][qp]*system->current_solution (dof_indices_w[l]);
		    }

				mean_u += u;
				mean_v += v;
				mean_w += w;

		  } // end of the quadrature point qp-loop

			mean_u /= qrule.n_points();
			mean_v /= qrule.n_points();
			mean_w /= qrule.n_points();

			double u_mag = sqrt(mean_u*mean_u + mean_v*mean_v + mean_w*mean_w);

			double peclet_number = u_mag * h_T * reynolds_number;
	
			double error = peclet_number;
			error_vector[elem->id()] += static_cast<ErrorVectorReal>(error);

		}
  } // end of element loop

  return;
}


double NSAssembler3D::calculate_l2_norm(const unsigned int var_to_calc, bool reference)
{

	//std::cout << "begin calculating flux on bouundary " << boundary_id << std::endl;
	
	TransientLinearImplicitSystem * system;
	// Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}

	// clear the rhs before adding things in for flux
  system->rhs->zero ();
	const MeshBase& mesh = es->get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  const unsigned int u_var = system->variable_number ("u");
  const unsigned int v_var = system->variable_number ("v");
	//hmmm yeah you know const an shit
	unsigned int w_var = 0;
	if(threed)
		w_var = system->variable_number ("w");
  const unsigned int p_var = system->variable_number ("p");

  FEType fe_vel_type = system->variable_type(u_var);
  FEType fe_pres_type = system->variable_type(p_var);
	//boundary finite element stuffs

  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
	AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));

	//default quad order plus 2
	int default_quad_order = static_cast<int>(fe_vel_type.default_quadrature_order());
	std::cout << "default_quad_order = " << default_quad_order << std::endl; 
	std::cout << "u_var = " << u_var << std::endl; 
	//QGauss qrule(dim, static_cast<Order>(default_quad_order));

	QGauss qrule(dim, static_cast<Order>(default_quad_order + 2));
	fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
 

	

	//variables for equation system
  const std::vector<Real>& JxW = fe_vel->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();
  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<Point >& qpoint = fe_vel->get_xyz();

	std::cout << "var_to_calc = " << var_to_calc << std::endl;
	std::cout << "u_var = " << u_var << std::endl;
	if(var_to_calc == u_var)
		std::cout << "gosh, npoints = " << qrule.n_points() << std::endl;


	
  DenseVector<Number> Fe;
  DenseSubVector<Number>  Fu(Fe), Fv(Fe), Fw(Fe), Fp(Fe);

	// some dof stuff
  const DofMap & dof_map = system->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_w;
  std::vector<dof_id_type> dof_indices_p;

	ExactSolutionVelocity<> exact_solution_vel_object(*es);
	ExactSolutionPressure<> exact_solution_pres_object(*es);

	double approx_solution_mask = 1.;
	if(reference)
		approx_solution_mask = 0.;

  //iterators
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	double local_integral = 0.;

  for ( ; el != end_el; ++el)
  {				

	  const Elem* elem = *el;
		if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
		{

		  dof_map.dof_indices (elem, dof_indices);
		  dof_map.dof_indices (elem, dof_indices_u, u_var);
		  dof_map.dof_indices (elem, dof_indices_v, v_var);
			if(threed) { dof_map.dof_indices (elem, dof_indices_w, w_var); }
		  dof_map.dof_indices (elem, dof_indices_p, p_var);

		  unsigned int n_dofs   = dof_indices.size();
		  unsigned int n_u_dofs = dof_indices_u.size();
		  unsigned int n_v_dofs = dof_indices_v.size();
			unsigned int n_w_dofs = 0;
			if(threed) { n_w_dofs = dof_indices_w.size(); }
		  unsigned int n_p_dofs = dof_indices_p.size();

		  // Compute shape functions etc of current element
		  Fe.resize (n_dofs);

		  // Reposition the submatrices
		  Fu.reposition (u_var*n_u_dofs, n_u_dofs);
		  Fv.reposition (v_var*n_u_dofs, n_v_dofs);
			if(threed) { Fw.reposition (w_var*n_u_dofs, n_w_dofs); }
		  Fp.reposition (p_var*n_u_dofs, n_p_dofs);

		  fe_vel->reinit  (elem);
		  fe_pres->reinit (elem);
		
		  // The following loops over the sides of the element.
		  // If the element has no neighbor on a side then that
		  // side MUST live on a boundary of the domain.


		  // Build volume contribution to element matrix and rhs
		  for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		  {
		    // Compute the velocity & its gradient from the previous timestep
		    // and the old Newton iterate.
				DenseVector<Number> exact_solution_vel(dim);
				double exact_solution_pres = 0.;
				exact_solution_vel_object(qpoint[qp],0.,exact_solution_vel);
				exact_solution_pres = exact_solution_pres_object(qpoint[qp],0.);

				double u=0., v=0.,w=0.,p=0.;

		    // Compute the velocity & its gradient from the previous timestep
		    // and the old Newton iterate.
		    for (unsigned int l=0; l<n_u_dofs; l++)
		    {
		      // From the previous Newton iterate:
		      u += phi[l][qp]*system->current_solution (dof_indices_u[l]);
		      v += phi[l][qp]*system->current_solution (dof_indices_v[l]);
					if(threed)
			      w += phi[l][qp]*system->current_solution (dof_indices_w[l]);
				}

				for (unsigned int l=0; l<n_p_dofs; l++)
				{
			    p += psi[l][qp]*system->current_solution (dof_indices_p[l]);
				}

				if(var_to_calc == u_var)
					local_integral += pow(approx_solution_mask*u - exact_solution_vel(0),2.) * JxW[qp];

				
				if(var_to_calc == v_var)
					local_integral += pow(approx_solution_mask*v - exact_solution_vel(1),2.) * JxW[qp];

				if(threed)
					if(var_to_calc == w_var)
						local_integral += pow(approx_solution_mask*w - exact_solution_vel(2),2.) * JxW[qp];
				
				if(var_to_calc == p_var)
		     	local_integral += pow(approx_solution_mask*p - exact_solution_pres,2.) * JxW[qp];
				

		  } // end of the quadrature point qp-loop

		}
  } // end of element loop

	es->comm().sum(local_integral);

	double l2_norm = local_integral;
	l2_norm = sqrt(l2_norm);

  return l2_norm;

}


void NSAssembler3D::assemble_preconditioner ()
{

	TransientLinearImplicitSystem * system;
	// Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}

	//here we just copy the matrix
	// the stokes pressure mass matrix preconditioner
	if(es->parameters.get<unsigned int>("preconditioner_type_3d") == 1)
	{
		system->get_matrix("Preconditioner").close();
	}
}



void NSAssembler3D::assemble_efficient (ErrorVector& error)
{
	std::cout << "Efficient assembly was not implemented for you assembly class." << std::endl;
	std::cout << "Using standard assembly." << std::endl;
	this->assemble(error);
}


void NSAssembler3D::assemble_residual_rhs (ErrorVector& error)
{
	std::cout << "assemble_residual_rhs was not implemented for you assembly class." << std::endl;
	std::cout << "Exiting." << std::endl;
	std::exit(0);
}



