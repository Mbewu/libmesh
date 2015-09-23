
#include "ns_assembler_3d.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

// boundary ids must have the inflow boundary id followed  by the outlfow 
// boundary ids in whatever order actually, must just match up
// these values are just magnitude values
void NSAssembler3D::init_bc (std::vector<unsigned int> boundary_ids,
											std::vector<double> pressure_values = std::vector<double>(),
											std::vector<double> previous_flow_values = std::vector<double>(),
											std::vector<double> previous_previous_flow_values = std::vector<double>()) 
{

	std::cout << "Setting 2D/3D BCs." << std::endl;
	
	//resize so that doesn't break ;)
	pressure_values.resize(boundary_ids.size());
	const int problem_type    = es->parameters.get<unsigned int>("problem_type");
	
	// dirichlet input
	if(problem_type == 0 || problem_type == 3)
	{
		bc_type[boundary_ids[0]] = "dirichlet";
		for(unsigned int i = 1; i < boundary_ids.size(); i++)
			bc_type[boundary_ids[i]] = "pressure";
		
		
			
		bc_value[boundary_ids[0]] = 0;	//dirichlet will be handled outside this class	
		for(unsigned int i = 1; i < boundary_ids.size(); i++)
			bc_value[boundary_ids[i]] = pressure_values[i];	//these should ideally be zero...

		//here we calculate the interpolated flow value
		//if steady then this is just 1.
		
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
	else if(problem_type == 1)
	{
	
		bc_type[boundary_ids[0]] = "pressure";
		for(unsigned int i = 1; i < boundary_ids.size(); i++)
			bc_type[boundary_ids[i]] = "pressure";
		

	
		bc_value[boundary_ids[0]] = pressure_values[0];;//"0.00675896";	//"0.008";
		for(unsigned int i = 1; i < boundary_ids.size(); i++)
			bc_value[boundary_ids[i]] = pressure_values[i];//"-100.0";


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
	else if(problem_type == 2)
	{
	
		bc_type[boundary_ids[0]] = "stress";
		for(unsigned int i = 1; i < boundary_ids.size(); i++)
			bc_type[boundary_ids[i]] = "stress";
		

	
		bc_value[boundary_ids[0]] = pressure_values[0];	//"0.008";
		for(unsigned int i = 1; i < boundary_ids.size(); i++)	
			bc_value[boundary_ids[i]] = pressure_values[i];


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

			

	

	if(problem_type != 4)
	{
		//note that these values should actually be read in
		Real inflow_centre[] = {0.0,0.0,1.5};
		std::vector<Real> inflow_centre_vec (inflow_centre, inflow_centre + sizeof(inflow_centre) / sizeof(Real) );
		boundary_centre[boundary_ids[0]] = inflow_centre_vec;

		// set boundary radius at inflow, for more irregular geometries will need more fancy mappings etc
		boundary_radius[boundary_ids[0]] = 0.5;
	}

	if(pressure_coupled)
	{
		find_1d_boundary_nodes();
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
			unsigned int current_1d_el_idx = current_el_idx -	n_initial_3d_elem;
			bool is_daughter_1 = airway_data[current_1d_el_idx].get_is_daughter_1();	//this is a bool duh!

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

					//std::cout << "primary_pressure_boundary_nodes_1d" << boundary_ids[0] << "  = " << primary_pressure_boundary_nodes_1d[boundary_ids[0]] << std::endl;
				}
			}		
		}	
	}

	

}


double NSAssembler3D::calculate_flux(const int boundary_id)
{

	std::cout << "begin calculating flux on bouundary " << boundary_id << std::endl;
	
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

	for(unsigned int i=0; i<subdomains_3d.size(); i++)
		std::cout << "subdomains[" << i << "] = " << subdomains_3d[i] << std::endl;


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

  return flux;

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
	
  return pressure;

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
	if(es->parameters.get<unsigned int>("preconditioner_type") == 1)
	{
		system->get_matrix("Preconditioner").close();
	}
}

