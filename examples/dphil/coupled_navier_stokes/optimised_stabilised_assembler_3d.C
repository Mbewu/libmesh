
#include "optimised_stabilised_assembler_3d.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

void OptimisedStabilisedAssembler3D::assemble(ErrorVector&)// error_vector)
{

	if(threed)
		std::cout << "Begin 3D assembly... ";
	else
		std::cout << "Begin 2D assembly... ";

	TransientLinearImplicitSystem * system;
	TransientLinearImplicitSystem * system_neumann;
  // Get a reference to the Stokes system object.
	if(pressure_coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("Coupled-Navier-Stokes");
	}
	else
	{
		system =
	  	&es->get_system<TransientLinearImplicitSystem> ("Navier-Stokes-3D");
	}

	if(es->parameters.get<bool> ("discontinuous_linear_neumann"))
		system_neumann =
  		&es->get_system<TransientLinearImplicitSystem> ("Neumann-Variable");

	// Problem parameters
  const Real dt    = es->parameters.get<Real>("dt");
	//const Real time = system->time;
  Real Re    = es->parameters.get<Real>("reynolds_number");
  const unsigned int unsteady    = es->parameters.get<unsigned int>("unsteady");
  const bool stokes    = es->parameters.get<bool>("stokes");
  const bool stab    = es->parameters.get<bool>("stab");
  const Real alpha    = es->parameters.get<Real>("alpha");
	//const Real period	= es->parameters.get<Real>("period");
	const bool newton	= es->parameters.get<bool>("newton");
	//const double density	= es->parameters.get<double>("density");
	const bool convective_form	= es->parameters.get<bool>("convective_form");
	const double backflow_stab_param	= es->parameters.get<double>("backflow_stab_param");
	const bool supg	= es->parameters.get<bool>("supg");
	//const bool pspg	= es->parameters.get<bool>("pspg");
	//const bool lsic	= es->parameters.get<bool>("lsic");


	// Mesh and system things
	const MeshBase& mesh = es->get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  const unsigned int u_var = system->variable_number ("u");
  const unsigned int v_var = system->variable_number ("v");
	unsigned int w_var = 0;
	if(threed)
		w_var = system->variable_number ("w");
  const unsigned int p_var = system->variable_number ("p");


  const unsigned int u_adj_var = system->variable_number ("u_adj");
  const unsigned int v_adj_var = system->variable_number ("v_adj");
	unsigned int w_adj_var = 0;
	if(threed)
		w_adj_var = system->variable_number ("w_adj");

	unsigned int g_x_var=0, g_y_var=0, g_z_var=0;
	if(es->parameters.get<bool> ("discontinuous_linear_neumann"))
	{
		g_x_var = system_neumann->variable_number ("g_x");
		g_y_var = system_neumann->variable_number ("g_y");
		if(threed)
			g_z_var = system_neumann->variable_number ("g_z");
	}



	// Finite element types - adjoint variables has same fe types
  FEType fe_vel_type = system->variable_type(u_var);
	FEType fe_pres_type = system->variable_type(p_var);
  FEType fe_neumann_type;
	if(es->parameters.get<bool> ("discontinuous_linear_neumann"))
		fe_neumann_type = system_neumann->variable_type(g_x_var);
	//FEType fe_neumann_type = system->variable_type(u_var);	//temp

  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
	AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
  AutoPtr<FEBase> fe_neumann  (FEBase::build(dim, fe_neumann_type));

	// Quadrature rules
	int default_quad_order = static_cast<int>(fe_vel_type.default_quadrature_order());
	//QGauss qrule(dim, static_cast<Order>(default_quad_order));
	QGauss qrule(dim, static_cast<Order>(default_quad_order + 6));

	std::cout << "volume quad rule = " << default_quad_order << std::endl;

  //QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
  fe_neumann->attach_quadrature_rule (&qrule);


	// Boundary finite element types
	AutoPtr<FEBase> fe_vel_face (FEBase::build(dim, fe_vel_type));
	AutoPtr<FEBase> fe_pres_face (FEBase::build(dim, fe_pres_type));
	AutoPtr<FEBase> fe_neumann_face (FEBase::build(dim, fe_neumann_type));


	QGauss qface(dim-1, static_cast<Order>(default_quad_order + 6));
	std::cout << "surface quad rule = " << default_quad_order << std::endl;
	//QGauss qface(dim-1, fe_vel_type.default_quadrature_order());
	fe_vel_face->attach_quadrature_rule (&qface);
	fe_pres_face->attach_quadrature_rule (&qface);
	fe_neumann_face->attach_quadrature_rule (&qface);


	// Initialise some variables
  const std::vector<Real>& JxW = fe_vel->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();
  const std::vector<std::vector<RealTensor> >& d2phi = fe_vel->get_d2phi();
  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();
  // const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();
  const std::vector<Real>& JxW_face = fe_vel_face->get_JxW();
  const std::vector<std::vector<Real> >& phi_face = fe_vel_face->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi_face = fe_vel_face->get_dphi();
  const std::vector<Point>& 						 qface_normals = fe_vel_face->get_normals();
	const std::vector< std::vector<Point> > &		qface_tangents = fe_vel_face->get_tangents();
  const std::vector<std::vector<Real> >& xi_face = fe_neumann_face->get_phi();



  // Local element matrix
  DenseMatrix<Number> Ke;
  DenseMatrix<Number> Mb;
  DenseMatrix<Number> Nb;
  DenseMatrix<Number> Nb_T;
  DenseVector<Number> Fe;
  DenseVector<Number> Fe_coupled_p;
  DenseVector<Number> Fe_coupled_u;

	// okay we don't need to make matrices that are zero..
  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke),
    Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke),
    Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke),
    Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke),

		//some extra matrices ;)
    Mbdy_uu(Ke), Mbdy_vv(Ke), Mbdy_ww(Ke),
		Auu(Ke), Avv(Ke), Aww(Ke);

	// okay we don't need to make matrices that are zero..
  DenseSubMatrix<Number>
    Mbxx(Mb), Mbyy(Mb), Mbzz(Mb),
    Nbxx(Nb), Nbyy(Nb), Nbzz(Nb),
		Nb_Txx(Nb_T), Nb_Tyy(Nb_T), Nb_Tzz(Nb_T), Mbdy(Ke);
		//this is the usual force vector
  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe),
    Fp(Fe);

		//this is the additional one due to navier stokes correction
  DenseSubVector<Number>
    Fu_adj(Fe),
    Fv_adj(Fe),
    Fw_adj(Fe);

  DenseSubVector<Number>
    Fu_coupled_p(Fe_coupled_p),
    Fv_coupled_p(Fe_coupled_p),
    Fw_coupled_p(Fe_coupled_p),
    Fp_coupled_p(Fe_coupled_p); //unnecessary


  DenseSubVector<Number>
    Fu_coupled_u(Fe_coupled_u),
    Fv_coupled_u(Fe_coupled_u),
    Fw_coupled_u(Fe_coupled_u),
    Fp_coupled_u(Fe_coupled_u); //unnecessary


	
	


	// DofMap things
  const DofMap & dof_map = system->get_dof_map();
  const DofMap & dof_map_neumann = system_neumann->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;
  std::vector<dof_id_type> dof_indices_w;
  std::vector<dof_id_type> dof_indices_p;

  std::vector<dof_id_type> dof_indices_u_adj;
  std::vector<dof_id_type> dof_indices_v_adj;
  std::vector<dof_id_type> dof_indices_w_adj;

  std::vector<dof_id_type> dof_indices_g_x;
  std::vector<dof_id_type> dof_indices_g_y;
  std::vector<dof_id_type> dof_indices_g_z;

	// 0 corresponds to the one actually on the boundary
	// 1 corresponds to the other dof on the element
	// let us build this map before for ease of coding
	dof_id_type dof_index_1d_0 = 0;//p0
	dof_id_type dof_index_1d_1 = 0;//p1
	//dof_id_type dof_index_1d_2 = 0;//q0


  bool pin_pressure = false;
	if(es->parameters.get<unsigned int>("problem_type") == 4)
		pin_pressure = true;


  // Iterators
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	double max_u = 0.0;
	//double u_mag = 0.0;

	NumberVectorValue total_residual(0,0);
	NumberVectorValue total_residual_1(0,0);
	NumberVectorValue total_residual_2(0,0);
	NumberVectorValue total_residual_3(0,0);
	NumberVectorValue total_residual_4(0,0);

		

	// the equation system will look like
	//
	//	| A 0  Kt  Bt ||  u  |   |  0   |
	//	| 0 0  B   0  ||  p  |   |  0   |
	//	| K Bt -ML 0  ||u_adj| = |f-Nb*g|
	//	| B 0  0   0  ||v_adj|   |  0   |
	//
	//
	// ********************************

	unsigned int count = 0;

  for ( ; el != end_el; ++el)
  {				
		count += 1;
		bool boundary_element = false;
    const Elem* elem = *el;
		// only solve on the 3d subdomain
		if(elem->subdomain_id() == 0)
		{
		  dof_map.dof_indices (elem, dof_indices);
		  dof_map.dof_indices (elem, dof_indices_u, u_var);
		  dof_map.dof_indices (elem, dof_indices_v, v_var);
			if(threed)
				dof_map.dof_indices (elem, dof_indices_w, w_var);
		  dof_map.dof_indices (elem, dof_indices_p, p_var);
		  dof_map.dof_indices (elem, dof_indices_u_adj, u_adj_var);
		  dof_map.dof_indices (elem, dof_indices_v_adj, v_adj_var);
			if(threed)
				dof_map.dof_indices (elem, dof_indices_w_adj, w_adj_var);

		  dof_map_neumann.dof_indices (elem, dof_indices_g_x, g_x_var);
		  dof_map_neumann.dof_indices (elem, dof_indices_g_y, g_y_var);
			if(threed)
				dof_map_neumann.dof_indices (elem, dof_indices_g_z, g_z_var);

		  unsigned int n_dofs   = dof_indices.size();
		  unsigned int n_u_dofs = dof_indices_u.size();
		  unsigned int n_v_dofs = dof_indices_v.size();
			unsigned int n_w_dofs = 0;
			if(threed)
				n_w_dofs = dof_indices_w.size();
		  unsigned int n_p_dofs = dof_indices_p.size();

		  unsigned int n_g_x_dofs = dof_indices_g_x.size();
		  unsigned int n_g_y_dofs = dof_indices_g_y.size();
			unsigned int n_g_z_dofs = 0;
			if(threed)
				n_g_z_dofs = dof_indices_g_z.size();

		  // Compute shape functions etc of current element
		  fe_vel->reinit  (elem);
		  fe_pres->reinit (elem);
		  fe_neumann->reinit  (elem);
			Real elem_volume = elem->volume();
			//take cube root to get diameter of the element
			double h_T = 0;
			if(threed)
				h_T = pow(elem_volume,1.0/3.0);
			else
				h_T = pow(elem_volume,1.0/2.0);

			//std::cout << "n_g_x_dofs = " << n_g_x_dofs << std::endl;
			//std::cout << "n_p_dofs = " << n_p_dofs << std::endl;
			//std::cout << "n_u_dofs = " << n_u_dofs << std::endl;

		  Ke.resize (n_dofs, n_dofs);
		  Mb.resize (dim*n_g_x_dofs, dim*n_g_x_dofs);
		  Nb.resize (dim*n_u_dofs, dim*n_g_x_dofs);
		  Nb_T.resize (dim*n_g_x_dofs,dim*n_u_dofs);
		  Fe.resize (n_dofs);
		  Fe_coupled_p.resize (n_dofs);
		  Fe_coupled_u.resize (n_dofs);

		  // Reposition the submatrices
			Kuu.reposition (u_var*n_u_dofs + n_dofs/2, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
		  Kuv.reposition (u_var*n_u_dofs + n_dofs/2, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
		  Kup.reposition (u_var*n_u_dofs + n_dofs/2, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
			if(threed)
				Kuw.reposition (u_var*n_u_dofs + n_dofs/2, w_var*n_u_dofs, n_u_dofs, n_w_dofs);

		  Kvu.reposition (v_var*n_v_dofs + n_dofs/2, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
		  Kvv.reposition (v_var*n_v_dofs + n_dofs/2, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
		  Kvp.reposition (v_var*n_v_dofs + n_dofs/2, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
			if(threed)
				Kvw.reposition (v_var*n_v_dofs + n_dofs/2, w_var*n_v_dofs, n_v_dofs, n_w_dofs);

			if(threed)
			{
				Kwu.reposition (w_var*n_w_dofs + n_dofs/2, u_var*n_w_dofs, n_w_dofs, n_u_dofs);
				Kwv.reposition (w_var*n_w_dofs + n_dofs/2, v_var*n_w_dofs, n_w_dofs, n_v_dofs);
				Kww.reposition (w_var*n_w_dofs + n_dofs/2, w_var*n_w_dofs, n_w_dofs, n_w_dofs);
				Kwp.reposition (w_var*n_w_dofs + n_dofs/2, p_var*n_w_dofs, n_w_dofs, n_p_dofs);
			}
		
		  Kpu.reposition (p_var*n_u_dofs + n_dofs/2, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
		  Kpv.reposition (p_var*n_u_dofs + n_dofs/2, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
		  Kpp.reposition (p_var*n_u_dofs + n_dofs/2, p_var*n_u_dofs, n_p_dofs, n_p_dofs);
			if(threed)
				Kpw.reposition (p_var*n_u_dofs + n_dofs/2, w_var*n_u_dofs, n_p_dofs, n_w_dofs);

		  Fu.reposition (u_var*n_u_dofs + n_dofs/2, n_u_dofs);
		  Fv.reposition (v_var*n_u_dofs + n_dofs/2, n_v_dofs);
		  Fp.reposition (p_var*n_u_dofs + n_dofs/2, n_p_dofs);
			if(threed)
				Fw.reposition (w_var*n_u_dofs + n_dofs/2, n_w_dofs);

		  Fu_adj.reposition (u_var*n_u_dofs, n_u_dofs);
		  Fv_adj.reposition (v_var*n_u_dofs, n_v_dofs);
			if(threed)
				Fw_adj.reposition (w_var*n_u_dofs, n_w_dofs);

			// nothing goes in the other block of the rhs

		  // Reposition the submatrices
			Auu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
		  Avv.reposition (v_var*n_u_dofs, v_var*n_u_dofs, n_v_dofs, n_v_dofs);
			if(threed)
				Aww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);

		  // Reposition the submatrices
			Mbdy_uu.reposition (u_var*n_u_dofs + n_dofs/2, u_var*n_u_dofs + n_dofs/2, n_u_dofs, n_u_dofs);
		  Mbdy_vv.reposition (v_var*n_u_dofs + n_dofs/2, v_var*n_u_dofs + n_dofs/2, n_v_dofs, n_v_dofs);
			if(threed)
				Mbdy_ww.reposition (w_var*n_u_dofs + n_dofs/2, w_var*n_u_dofs + n_dofs/2, n_w_dofs, n_w_dofs);



		  Fu_coupled_p.reposition (u_var*n_u_dofs, n_u_dofs);
		  Fv_coupled_p.reposition (v_var*n_u_dofs, n_v_dofs);
		  Fp_coupled_p.reposition (p_var*n_u_dofs, n_p_dofs);
			if(threed)
				Fw_coupled_p.reposition (w_var*n_u_dofs, n_w_dofs);

		  Fu_coupled_u.reposition (u_var*n_u_dofs, n_u_dofs);
		  Fv_coupled_u.reposition (v_var*n_u_dofs, n_v_dofs);
		  Fp_coupled_u.reposition (p_var*n_u_dofs, n_p_dofs);
			if(threed)
				Fw_coupled_u.reposition (w_var*n_u_dofs, n_w_dofs);

			Mbxx.reposition (g_x_var*n_g_x_dofs, g_x_var*n_g_x_dofs, n_g_x_dofs, n_g_x_dofs);
			Mbyy.reposition (g_y_var*n_g_y_dofs, g_y_var*n_g_y_dofs, n_g_y_dofs, n_g_y_dofs);
			if(threed)
				Mbzz.reposition (g_z_var*n_g_z_dofs, g_z_var*n_g_z_dofs, n_g_z_dofs, n_g_z_dofs);



			Nbxx.reposition (u_var*n_u_dofs, g_x_var*n_g_x_dofs, n_u_dofs, n_g_x_dofs);
			Nbyy.reposition (v_var*n_u_dofs, g_y_var*n_g_y_dofs, n_v_dofs, n_g_y_dofs);
			if(threed)
				Nbzz.reposition (w_var*n_u_dofs, g_z_var*n_g_z_dofs, n_w_dofs, n_g_z_dofs);


			Nb_Txx.reposition (g_x_var*n_g_x_dofs,u_var*n_u_dofs,  n_g_x_dofs, n_u_dofs);
			Nb_Tyy.reposition (g_y_var*n_g_y_dofs, v_var*n_u_dofs, n_g_y_dofs, n_v_dofs);
			if(threed)
				Nb_Tzz.reposition (g_z_var*n_g_z_dofs, w_var*n_u_dofs, n_g_z_dofs, n_w_dofs);


		  // Reposition the submatrices
			Mbdy.reposition (u_var*n_u_dofs + n_dofs/2, u_var*n_u_dofs + n_dofs/2, dim*n_u_dofs, dim*n_u_dofs);

		  // Build volume contribution to element matrix and rhs
		  for (unsigned int qp=0; qp<qrule.n_points(); qp++)
		  {
		    // Values to hold the solution & its gradient at the previous timestep.
		    Number   u = 0., u_old = 0., v = 0., v_old = 0., w = 0., w_old = 0.;
		    Number   u_adj = 0., v_adj = 0., w_adj = 0.;
		    Gradient grad_u, grad_v, grad_w, grad_p;
				Number lap_u = 0., lap_v = 0., lap_w = 0.;

		    // Compute the velocity & its gradient from the previous timestep
		    // and the old Newton iterate.
		    for (unsigned int l=0; l<n_u_dofs; l++)
		    {
		      // From the old timestep:
		      u_old += phi[l][qp]*system->old_solution (dof_indices_u[l]);
		      v_old += phi[l][qp]*system->old_solution (dof_indices_v[l]);
					if(threed)
			      w_old += phi[l][qp]*system->old_solution (dof_indices_w[l]);

		      // From the previous Newton iterate:
		      u += phi[l][qp]*system->current_solution (dof_indices_u[l]);
		      v += phi[l][qp]*system->current_solution (dof_indices_v[l]);
					if(threed)
			      w += phi[l][qp]*system->current_solution (dof_indices_w[l]);

		      // Adjoint solution from the previous Newton iterate:
		      u_adj += phi[l][qp]*system->current_solution (dof_indices_u_adj[l]);
		      v_adj += phi[l][qp]*system->current_solution (dof_indices_v_adj[l]);
					if(threed)
			      w_adj += phi[l][qp]*system->current_solution (dof_indices_w_adj[l]);

		      grad_u.add_scaled (dphi[l][qp],system->current_solution (dof_indices_u[l]));
		      grad_v.add_scaled (dphi[l][qp],system->current_solution (dof_indices_v[l]));
					if(threed)
			      grad_w.add_scaled (dphi[l][qp],system->current_solution (dof_indices_w[l]));

					double laplacian_operator = 0;
					laplacian_operator += d2phi[l][qp](0,0) +	d2phi[l][qp](1,1);
					if(threed) {laplacian_operator += d2phi[l][qp](2,2); }

		      lap_u += laplacian_operator*system->current_solution (dof_indices_u[l]);
		      lap_v += laplacian_operator*system->current_solution (dof_indices_v[l]);
					if(threed)
			      lap_w += laplacian_operator*system->current_solution (dof_indices_w[l]);

		    }

				for(unsigned int l=0; l<n_p_dofs; l++)
				{
		      grad_p.add_scaled (dpsi[l][qp],system->current_solution (dof_indices_p[l]));
				}

				if(pow(pow(u,2.0) + pow(v,2.0),0.5) > max_u)
					max_u = pow(pow(u,2.0) + pow(v,2.0),0.5);

				/*
				if(threed)
					u_mag = pow(pow(u,2.0) + pow(v,2.0) + pow(w,2.0),0.5);
				else
					u_mag = pow(pow(u,2.0) + pow(v,2.0),0.5);
				*/

				//if(fabs(u_old) > max_u)
				//	max_u =fabs(u_old);

				// construct vectors we will need in assembly


		    NumberVectorValue U_old;
		    NumberVectorValue U;
		    NumberVectorValue U_adj;
		    NumberVectorValue dU_dx,dU_dy,dU_dz;

				if(threed)
				{
					U = NumberVectorValue(u, v, w);
					U_adj = NumberVectorValue(u_adj, v_adj, w_adj);
					U_old = NumberVectorValue(u_old, v_old, w_old);
					dU_dx = NumberVectorValue(grad_u(0), grad_v(0), grad_w(0));
					dU_dy = NumberVectorValue(grad_u(1), grad_v(1), grad_w(1));
					dU_dz = NumberVectorValue(grad_u(2), grad_v(2), grad_w(2));
				}
				else
				{
					U = NumberVectorValue(u, v);
					U_adj = NumberVectorValue(u_adj, v_adj);
					U_old = NumberVectorValue(u_old, v_old);
					dU_dx = NumberVectorValue(grad_u(0), grad_v(0));
					dU_dy = NumberVectorValue(grad_u(1), grad_v(1));
				}
			


				//if(fabs(w_old) > 0)
				//	std::cout << "w_old = " << w_old << std::endl;

				//	if(fabs(w) > 0)
				//		std::cout << "w = " << w << std::endl;
		    // First, an i-loop over the velocity degrees of freedom.
				// need to be in the adj block now
		    for (unsigned int i=0; i<n_u_dofs; i++)
		    {

					// Unsteady term
					if(unsteady)
					{
		        Fu(i) -= -JxW[qp]*(u_old*phi[i][qp]);
		        Fv(i) -= -JxW[qp]*(v_old*phi[i][qp]);
						if(threed) {	Fw(i) -= -JxW[qp]*(w_old*phi[i][qp]); }
					}

					// Newton term

					if(!stokes && newton)
					{
						if(convective_form)
						{
					    Fu(i) += JxW[qp]*dt*((U*grad_u)*phi[i][qp]);
					    Fv(i) += JxW[qp]*dt*((U*grad_v)*phi[i][qp]);
							if(threed) { Fw(i) += JxW[qp]*dt*((U*grad_w)*phi[i][qp]); }
						}
						else
						{
					    Fu(i) += -JxW[qp]*dt*(U*dphi[i][qp])*u;
					    Fv(i) += -JxW[qp]*dt*(U*dphi[i][qp])*v;
							if(threed) {	Fw(i) += -JxW[qp]*dt*(U*dphi[i][qp])*w; }
						}
					}		


		      // Matrix contributions for the uu and vv couplings. i.e. first equation
		      for (unsigned int j=0; j<n_u_dofs; j++)
		      {
						// Unsteady term i.e. mass matrix
						if(unsteady)
						{
		        	Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);                // mass matrix term
		        	Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);                // mass matrix term
							if(threed) {	Kww(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]); }                // mass matrix term
						}



	          Kuu(i,j) += JxW[qp]*(dt/Re*(dphi[i][qp]*dphi[j][qp]));   // diffusion
	          Kvv(i,j) += JxW[qp]*(dt/Re*(dphi[i][qp]*dphi[j][qp]));   // diffusion
						if(threed) { Kww(i,j) += JxW[qp]*(dt/Re*(dphi[i][qp]*dphi[j][qp]));}   // diffusion
						// Convection and diffusion terms
						if(!stokes)
						{
		                             
							if(convective_form)
							{
				        Kuu(i,j) += JxW[qp]*(dt*(U*dphi[j][qp])*phi[i][qp]);  // convection
				        Kvv(i,j) += JxW[qp]*(dt*(U*dphi[j][qp])*phi[i][qp]);   // convection
								if(threed){ Kww(i,j) += JxW[qp]*(dt*(U*dphi[j][qp])*phi[i][qp]);}   // convection
							}
							else
							{
				        Kuu(i,j) += -JxW[qp]*(dt*(U*dphi[i][qp])*phi[j][qp]);  // convection
				        Kvv(i,j) += -JxW[qp]*(dt*(U*dphi[i][qp])*phi[j][qp]);   // convection
				        if(threed) { Kww(i,j) += -JxW[qp]*(dt*(U*dphi[i][qp])*phi[j][qp]);}   // convection
							}

														

							if(newton)
							{
								if(convective_form)
								{
					      	Kuu(i,j) += JxW[qp]*(dt*grad_u(0)*phi[i][qp]*phi[j][qp]);
					      	Kuv(i,j) += JxW[qp]*(dt*grad_u(1)*phi[i][qp]*phi[j][qp]);                
					      	if(threed) { Kuw(i,j) += JxW[qp]*(dt*grad_u(2)*phi[i][qp]*phi[j][qp]); }          
					      	Kvu(i,j) += JxW[qp]*(dt*grad_v(0)*phi[i][qp]*phi[j][qp]);                
					      	Kvv(i,j) += JxW[qp]*(dt*grad_v(1)*phi[i][qp]*phi[j][qp]);                
					      	if(threed) { Kvw(i,j) += JxW[qp]*(dt*grad_v(2)*phi[i][qp]*phi[j][qp]); }             
									if(threed) { Kwu(i,j) += JxW[qp]*(dt*grad_w(0)*phi[i][qp]*phi[j][qp]); }             
									if(threed) { Kwv(i,j) += JxW[qp]*(dt*grad_w(1)*phi[i][qp]*phi[j][qp]); }
									if(threed) { Kww(i,j) += JxW[qp]*(dt*grad_w(2)*phi[i][qp]*phi[j][qp]); }
								}
								else
								{
					      	Kuu(i,j) += -JxW[qp]*(dt*u*dphi[i][qp](0)*phi[j][qp]);
					      	Kuv(i,j) += -JxW[qp]*(dt*u*dphi[i][qp](1)*phi[j][qp]);                
					      	if(threed) { Kuw(i,j) += -JxW[qp]*(dt*u*dphi[i][qp](2)*phi[j][qp]); }                
					      	Kvu(i,j) += -JxW[qp]*(dt*v*dphi[i][qp](0)*phi[j][qp]);                
					      	Kvv(i,j) += -JxW[qp]*(dt*v*dphi[i][qp](1)*phi[j][qp]);                
					      	if(threed) { Kvw(i,j) += -JxW[qp]*(dt*v*dphi[i][qp](2)*phi[j][qp]); }              
									if(threed) { Kwu(i,j) += -JxW[qp]*(dt*w*dphi[i][qp](0)*phi[j][qp]); }              
									if(threed) { Kwv(i,j) += -JxW[qp]*(dt*w*dphi[i][qp](1)*phi[j][qp]); }
									if(threed) { Kww(i,j) += -JxW[qp]*(dt*w*dphi[i][qp](2)*phi[j][qp]); }
								}
							}
						}

		      }

		      // Pressure gradient term
					// put the density in here
		      for (unsigned int j=0; j<n_p_dofs; j++)
		      {
		        Kup(i,j) += -JxW[qp]*(dt*psi[j][qp]*dphi[i][qp](0));
		        Kvp(i,j) += -JxW[qp]*(dt*psi[j][qp]*dphi[i][qp](1));
						if(threed) { Kwp(i,j) += -JxW[qp]*(dt*psi[j][qp]*dphi[i][qp](2)); }
		      }
		    }

		    // Continuity equation - could possibly be NOT multiplied by negative.
		    for (unsigned int i=0; i<n_p_dofs; i++)
		    {
		      for (unsigned int j=0; j<n_u_dofs; j++)
		      {
		        Kpu(i,j) += -JxW[qp]*dt*psi[i][qp]*dphi[j][qp](0);
		        Kpv(i,j) += -JxW[qp]*dt*psi[i][qp]*dphi[j][qp](1);
						if(threed) { Kpw(i,j) += -JxW[qp]*dt*psi[i][qp]*dphi[j][qp](2); }
		      }

					// Stabilisation term - should be opposite sign to continuity equation
					// stabilisation depends on if using P1-P1 or P1-P0
					if(stab || es->parameters.get<bool> ("p2p1_stabilisation") )
					{
		      	for (unsigned int j=0; j<n_p_dofs; j++)
		        {
		          //Kpp(i,j) += alpha*elem_volume*JxW[qp]*dt*dpsi[i][qp]*dpsi[j][qp];
							Kpp(i,j) += -alpha*h_T*h_T*JxW[qp]*dt*dpsi[i][qp]*dpsi[j][qp];
		        }
					}
		    }

				// extra term for making oseen problem convex i thunk
				if(!es->parameters.get<bool> ("stokes") && !es->parameters.get<bool> ("newton")
						&& es->parameters.get<bool> ("opt_correction_term"))
				{
					for (unsigned int i=0; i<n_u_dofs; i++)
		    	{
			      Fu_adj(i) += -JxW[qp]*(phi[i][qp]*dU_dx*U_adj);
			      Fv_adj(i) += -JxW[qp]*(phi[i][qp]*dU_dy*U_adj);
						if(threed) {	Fw_adj(i) += -JxW[qp]*(phi[i][qp]*dU_dz*U_adj); }
					}
				}
				
				double scaling_factor = es->parameters.get<double>("optimisation_scaling_factor");
				if(es->parameters.get<bool>("volume_optimisation_h_scaled"))
					scaling_factor /= h_T;		//want bigger so stronger and gradient is 1/h_T relative to mass matrix

				if(es->parameters.get<bool>("volume_optimisation"))
				{
					for (unsigned int i=0; i<n_u_dofs; i++)
		    	{
						for (unsigned int j=0; j<n_u_dofs; j++)
				  	{
							// term corresponding to minimising gradients, the dt is there to keep terms a consistent size
							//if(!es->parameters.get<bool>("tangential_velocity_optimisation"))

							if(!es->parameters.get<bool>("pearson_example"))
							{
								if(!es->parameters.get<bool>("tangential_velocity_optimisation"))
									Auu(i,j) += scaling_factor * JxW[qp]*(1./Re * dphi[i][qp]*dphi[j][qp]);
								Avv(i,j) += scaling_factor * JxW[qp]*(1./Re * dphi[i][qp]*dphi[j][qp]);
								if(threed) { Aww(i,j) += scaling_factor * JxW[qp]*(1./Re * dphi[i][qp]*dphi[j][qp]); }
							}
							else
							{
								Auu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
								Avv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
								if(threed) { Aww(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]); }

								Mbdy_uu(i,j) += -1./scaling_factor * JxW[qp]*(phi[i][qp]*phi[j][qp]);
								Mbdy_vv(i,j) += -1./scaling_factor * JxW[qp]*(phi[i][qp]*phi[j][qp]);
								if(threed) { Mbdy_ww(i,j) += -1./scaling_factor * JxW[qp]*(phi[i][qp]*phi[j][qp]); }
							}
						}
					}
				}

				//SUPG terms (in a picard fashion)
				if(true)//supg)
		    {
					// petrov galerkin velocity - if unsteady use previous timestep ala tedzuyar, otherwise use current (picard)
					NumberVectorValue velocity;
					if(unsteady)
						velocity = U_old;
						//velocity = NumberVectorValue(1., 0.);
					else
						velocity = U;


					velocity = U;
					//supg parameter
					double tau_supg = 1.0/sqrt(pow(2*velocity.size()/h_T,2.0) + 9.0*pow(4/(h_T*h_T*Re),2.0));
					//std::cout << "tau_supg = " << tau_supg << std::endl;					
					//tau_supg /= 2.0;
					//tau_supg = 0.0;
					
					//note in the residual we ignore the integrals over boundaries cause int by parts not done
				  NumberVectorValue residual;
					if(threed)
						residual = NumberVectorValue(0, 0, 0);
					else
						residual = NumberVectorValue(0, 0);

					if(unsteady)
						residual += (U - U_old)/dt;

					total_residual_1 += JxW[qp]*residual;
					if(!stokes)
					{
						residual(0) += U * grad_u;
						residual(1) += U * grad_v;
						if(threed) { residual(2) += U * grad_w; }

					}

					total_residual_2 += JxW[qp]*residual;
					residual += grad_p;
					total_residual_3 += velocity;
					//std::cout << "grad_p = " << sqrt(grad_p(0)*grad_p(0) + grad_p(1)*grad_p(1)) << std::endl;


				
					//if first order then this term is zero
					if(!stab)
					{
						residual(0) += -1.0/Re * lap_u;
						residual(1) += -1.0/Re * lap_v;
						if(threed) { residual(2) += -1.0/Re * lap_w;}
					}
					total_residual_4 += JxW[qp]*residual;

					total_residual += JxW[qp]*residual;


					// let us use the velocity from the previous
					// there is a dt from the actual time derivative as well as from inside the residual
					if(supg)
				  {
						if(es->parameters.get<bool>("supg_newton"))
						{

							
							for (unsigned int i=0; i<n_u_dofs; i++)
							{
		
								if(unsteady)
								{
								  Fu(i) -= -tau_supg*JxW[qp]*dt/dt*(velocity*dphi[i][qp]*u_old);
								  Fv(i) -= -tau_supg*JxW[qp]*dt/dt*(velocity*dphi[i][qp]*v_old);
									if(threed) {	Fw(i) -= -tau_supg*JxW[qp]*dt/dt*(velocity*dphi[i][qp]*w_old); }
								}

								//the newton term bro
							  Fu(i) += tau_supg*JxW[qp]*dt*((velocity*dphi[i][qp])*(U*grad_u));
							  Fv(i) += tau_supg*JxW[qp]*dt*((velocity*dphi[i][qp])*(U*grad_v));
								if(threed) {	Fw(i) += tau_supg*JxW[qp]*dt*((velocity*dphi[i][qp])*(U*grad_w)); }
								


								for (unsigned int j=0; j<n_u_dofs; j++)
								{

									double laplacian_operator = 0;
									laplacian_operator += d2phi[j][qp](0,0) +	d2phi[j][qp](1,1);
									if(threed) {laplacian_operator += d2phi[j][qp](2,2); }

							  	Kuu(i,j) += 							tau_supg*JxW[qp]*dt*(phi[j][qp]*grad_u(0) * velocity*dphi[i][qp] + (U*dphi[j][qp]) * (velocity*dphi[i][qp]) 
																												- 1.0/Re * velocity*dphi[i][qp] * laplacian_operator);
							  	Kuv(i,j) += 							tau_supg*JxW[qp]*dt*(phi[j][qp]*grad_u(1) * velocity*dphi[i][qp] + (U*dphi[j][qp]) * (velocity*dphi[i][qp])
																												- 1.0/Re * velocity*dphi[i][qp] * laplacian_operator);                
							  	if(threed) { Kuw(i,j) += tau_supg*JxW[qp]*dt*(phi[j][qp]*grad_u(2) * velocity*dphi[i][qp] + (U*dphi[j][qp]) * (velocity*dphi[i][qp])
																												- 1.0/Re * velocity*dphi[i][qp] * laplacian_operator); }          
							  	Kvu(i,j) += 							tau_supg*JxW[qp]*dt*(phi[j][qp]*grad_v(0) * velocity*dphi[i][qp] + (U*dphi[j][qp]) * (velocity*dphi[i][qp])
																												- 1.0/Re * velocity*dphi[i][qp] * laplacian_operator);                
							  	Kvv(i,j) += 							tau_supg*JxW[qp]*dt*(phi[j][qp]*grad_v(1) * velocity*dphi[i][qp] + (U*dphi[j][qp]) * (velocity*dphi[i][qp])
																												- 1.0/Re * velocity*dphi[i][qp] * laplacian_operator);                
							  	if(threed) { Kvw(i,j) += tau_supg*JxW[qp]*dt*(phi[j][qp]*grad_v(2) * velocity*dphi[i][qp] + (U*dphi[j][qp]) * (velocity*dphi[i][qp])
																												- 1.0/Re * velocity*dphi[i][qp] * laplacian_operator); }             
									if(threed) { Kwu(i,j) += tau_supg*JxW[qp]*dt*(phi[j][qp]*grad_w(0) * velocity*dphi[i][qp] + (U*dphi[j][qp]) * (velocity*dphi[i][qp])
																												- 1.0/Re * velocity*dphi[i][qp] * laplacian_operator); }             
									if(threed) { Kwv(i,j) += tau_supg*JxW[qp]*dt*(phi[j][qp]*grad_w(1) * velocity*dphi[i][qp] + (U*dphi[j][qp]) * (velocity*dphi[i][qp])
																												- 1.0/Re * velocity*dphi[i][qp] * laplacian_operator); }
									if(threed) { Kww(i,j) += tau_supg*JxW[qp]*dt*(phi[j][qp]*grad_w(2) * velocity*dphi[i][qp] + (U*dphi[j][qp]) * (velocity*dphi[i][qp])
																												- 1.0/Re * velocity*dphi[i][qp] * laplacian_operator); }

									if(unsteady)
									{
										Kuu(i,j) += tau_supg*JxW[qp]*dt*(1.0/dt*velocity*dphi[i][qp]*phi[j][qp]);
										Kuv(i,j) += tau_supg*JxW[qp]*dt*(1.0/dt*velocity*dphi[i][qp]*phi[j][qp]);                
										if(threed) { Kuw(i,j) += tau_supg*JxW[qp]*dt*(1.0/dt*velocity*dphi[i][qp]*phi[j][qp]); }          
										Kvu(i,j) += tau_supg*JxW[qp]*dt*(1.0/dt*velocity*dphi[i][qp]*phi[j][qp]);                
										Kvv(i,j) += tau_supg*JxW[qp]*dt*(1.0/dt*velocity*dphi[i][qp]*phi[j][qp]);                
										if(threed) { Kvw(i,j) += tau_supg*JxW[qp]*dt*(1.0/dt*velocity*dphi[i][qp]*phi[j][qp]); }             
										if(threed) { Kwu(i,j) += tau_supg*JxW[qp]*dt*(1.0/dt*velocity*dphi[i][qp]*phi[j][qp]); }             
										if(threed) { Kwv(i,j) += tau_supg*JxW[qp]*dt*(1.0/dt*velocity*dphi[i][qp]*phi[j][qp]); }
										if(threed) { Kww(i,j) += tau_supg*JxW[qp]*dt*(1.0/dt*velocity*dphi[i][qp]*phi[j][qp]); }

									}
								}

								if(es->parameters.get<bool>("supg_picard_pressure"))
								{

									//the newton term bro
									Fu(i) += tau_supg*JxW[qp]*dt*(velocity*dphi[i][qp]*grad_p(0));
									Fv(i) += tau_supg*JxW[qp]*dt*(velocity*dphi[i][qp]*grad_p(1));
									if(threed) {	Fw(i) += tau_supg*JxW[qp]*dt*(velocity*dphi[i][qp]*grad_p(2)); }
								}
								else
								{
									for (unsigned int j=0; j<n_p_dofs; j++)
									{
										Kup(i,j) += tau_supg*JxW[qp]*dt*(velocity*dphi[i][qp] * dpsi[j][qp](0));      
										Kvp(i,j) += tau_supg*JxW[qp]*dt*(velocity*dphi[i][qp] * dpsi[j][qp](1));
										if(threed) { Kwp(i,j) += tau_supg*JxW[qp]*dt*(velocity*dphi[i][qp] * dpsi[j][qp](2)); } 
									}
								}

							}

						}
						//else if(es->parameters.get<bool>("supg_newton"))
						//{
						//}
						else
						{
							/*
							for (unsigned int i=0; i<n_u_dofs; i++)
							{
								for (unsigned int j=0; j<n_u_dofs; j++)
								{
							  	Kuu(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](0)*residual(0)*phi[j][qp]);
							  	Kuv(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](1)*residual(0)*phi[j][qp]);                
							  	if(threed) { Kuw(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](2)*residual(0)*phi[j][qp]); }          
							  	Kvu(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](0)*residual(1)*phi[j][qp]);                
							  	Kvv(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](1)*residual(1)*phi[j][qp]);                
							  	if(threed) { Kvw(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](2)*residual(1)*phi[j][qp]); }             
									if(threed) { Kwu(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](0)*residual(2)*phi[j][qp]); }             
									if(threed) { Kwv(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](1)*residual(2)*phi[j][qp]); }
									if(threed) { Kww(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](2)*residual(2)*phi[j][qp]); }
								}
							}
							*/
						}
					}

				}

		  } // end of the quadrature point qp-loop

			

		  // Apply boundary terms/conditions
		  {
		    // The penalty value.  \f$ \frac{1}{\epsilon} \f$
		    const Real penalty = 1.e10;

		    for (unsigned int s=0; s<elem->n_sides(); s++)
				{
					//for some reason it is natural to have more than one boundary id per side or even node
					std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

					if(boundary_ids.size() > 0) 
					{ 
						int boundary_id = boundary_ids[0];	// should only have one
 
						// stress boundary conditions/terms, which are only applied if not coupled
						// otherwise we put these contributions 
						// will have to think about this wrt the stabilisation term on the boundary
						if(!pressure_coupled)
						{
				
							// mean pressure conditions/terms
							if(bc_type[boundary_id].compare("pressure") == 0 || bc_type[boundary_id].compare("stress") == 0 || bc_type[boundary_id].compare("neumann") == 0)
							{

								boundary_element = true;
								double mean_pressure = bc_value[boundary_id];

								double area = (*surface_boundaries)[boundary_id]->get_area();

								if(bc_type[boundary_id].compare("neumann") == 0)
									mean_pressure = 0;

								// if coupled then we want the inflow to possibly be timeflow controlled, but future work
								if(es->parameters.get<unsigned int>("sim_type") == 0 || es->parameters.get<unsigned int>("sim_type") == 2)
									mean_pressure *= es->parameters.get<double>("time_scaling");

								//at the moment we just to a vector in normal direction for stress
								double stress_mag = bc_value[boundary_id];
								if(es->parameters.get<unsigned int>("sim_type") == 0 || es->parameters.get<unsigned int>("sim_type") == 2)
									stress_mag *= es->parameters.get<double>("time_scaling");
					
								DenseMatrix<Number> stress;
								DenseVector<Number> normal_stress;
								unsigned int dimension = 3;
								if(!threed)
									dimension = 2;
								stress.resize (dimension,dimension);
								normal_stress.resize (dimension);

								//make stress identity times constant
								stress(0,0) = 1.0; stress(1,1) = 1.0;
								if(threed)
									stress(2,2) = 1.0;
								stress *= stress_mag;

				
								fe_vel_face->reinit(elem, s);
								fe_pres_face->reinit(elem, s);
								fe_neumann_face->reinit(elem, s);

								//std::cout << "boundary_id = " << boundary_id << std::endl;
								//std::cout << "mean_pressure = " << mean_pressure << std::endl;
								//std::cout << "area = " << area << std::endl;

								for (unsigned int qp=0; qp<qface.n_points(); qp++)
								{

									if(bc_type[boundary_id].compare("stress") == 0)
									{
										normal_stress(0) = stress(0,0) * qface_normals[qp](0)	+ stress(0,1) * qface_normals[qp](1);
										if(threed) { normal_stress(0) += stress(0,2) * qface_normals[qp](2); }

										normal_stress(1) = stress(1,0) * qface_normals[qp](0)	+ stress(1,1) * qface_normals[qp](1);
										if(threed) { normal_stress(1) += stress(1,2) * qface_normals[qp](2); }

										if(threed)
											normal_stress(2) = stress(2,0) * qface_normals[qp](0)
																				+ stress(2,1) * qface_normals[qp](1)
																				+ stress(2,2) * qface_normals[qp](2);
									}
									else if(bc_type[boundary_id].compare("pressure") == 0)
									{
										normal_stress(0) = mean_pressure*qface_normals[qp](0);
										normal_stress(1) = mean_pressure*qface_normals[qp](1);
										if(threed)
											normal_stress(2) = mean_pressure*qface_normals[qp](2);
									}


									/*	for control problem the pressuure is enforced below
									for (unsigned int i=0; i<n_u_dofs; i++)
									{
							
								    Fu(i) += -JxW_face[qp]*(dt*((normal_stress(0))*phi_face[i][qp]));
								    Fv(i) += -JxW_face[qp]*(dt*((normal_stress(1))*phi_face[i][qp]));
								    if(threed) { Fw(i) += -JxW_face[qp]*(dt*((normal_stress(2))*phi_face[i][qp])); }
									}
									*/
		
									double scaling_factor = es->parameters.get<double>("optimisation_scaling_factor");
									// terms corresponding to minising the neumann condition

									// NOTE: probably have to use FACE values for this DUH
									for (unsigned int i=0; i<n_u_dofs; i++)
									{
										// this needs to go in the adj row and does
										Fu(i) += -JxW_face[qp]*(dt*(phi_face[i][qp] * normal_stress(0)));
										Fv(i) += -JxW_face[qp]*(dt*(phi_face[i][qp] * normal_stress(1)));
										if(threed) { Fw(i) += -JxW_face[qp]*(dt*(phi_face[i][qp] * normal_stress(2))); }

										for (unsigned int j=0; j<n_u_dofs; j++)
										{
											if(!es->parameters.get<bool>("discontinuous_linear_neumann"))	
											{
												if(es->parameters.get<bool>("minimise_neumann_boundary"))// hmmm not sure bout this && !es->parameters.get<bool>("pressure_minimisation"))
												{
													if(es->parameters.get<bool>("pressure_minimisation"))
													{
														Mbdy_uu(i,j) += -JxW_face[qp]/area*(dt*(phi_face[i][qp]*qface_normals[qp](0)) * (phi_face[j][qp]*qface_normals[qp](0)));
														Mbdy_vv(i,j) += -JxW_face[qp]/area*(dt*(phi_face[i][qp]*qface_normals[qp](1)) * (phi_face[j][qp]*qface_normals[qp](1)));
														if(threed) { Mbdy_ww(i,j) += -JxW_face[qp]/area*(dt*(phi_face[i][qp]*qface_normals[qp](2)) * (phi_face[j][qp]*qface_normals[qp](2))); }
													}
													else
													{
														Mbdy_uu(i,j) += -JxW_face[qp]*(dt*(phi_face[i][qp]*phi_face[j][qp]));
														Mbdy_vv(i,j) += -JxW_face[qp]*(dt*(phi_face[i][qp]*phi_face[j][qp]));
														if(threed) { Mbdy_ww(i,j) += -JxW_face[qp]*(dt*(phi_face[i][qp]*phi_face[j][qp])); }
													}
													
												}
												
											}
							
											if(es->parameters.get<bool>("tangential_optimisation"))
											{
												unsigned int num_tangents = qface_tangents[qp].size();
												if(!threed)
													num_tangents = 1;


												for(unsigned int k=0; k<num_tangents; k++)
												{
													//std::cout << "tangent(" << k << ") = " << qface_tangents[qp][k] << std::endl;
													// term corresponding to minimising gradients, the dt is there to keep terms a consistent size
													//if(!es->parameters.get<bool>("tangential_velocity_optimisation"))
														Auu(i,j) += scaling_factor * JxW_face[qp]*(1./Re *(dphi_face[i][qp]*qface_tangents[qp][k](0)*dphi_face[j][qp]*qface_tangents[qp][k](0)));
													Avv(i,j) += scaling_factor * JxW_face[qp]*(1./Re *(dphi_face[i][qp]*qface_tangents[qp][k](1)*dphi_face[j][qp]*qface_tangents[qp][k](1)));
													if(threed) { Aww(i,j) += scaling_factor * JxW_face[qp]*(dt*(dphi_face[i][qp]*qface_tangents[qp][k](2)*dphi_face[j][qp]*qface_tangents[qp][k](2))); }
												}
											}
											else if(!es->parameters.get<bool>("volume_optimisation")) // ie normal, non tangential
											{
												// term corresponding to minimising gradients, the dt is there to keep terms a consistent size
												if(!es->parameters.get<bool>("tangential_velocity_optimisation"))
													Auu(i,j) += scaling_factor * JxW_face[qp]*(1./Re*(dphi_face[i][qp]*dphi_face[j][qp]));
												Avv(i,j) += scaling_factor * JxW_face[qp]*(1./Re*(dphi_face[i][qp]*dphi_face[j][qp]));
												if(threed) { Aww(i,j) += scaling_factor * JxW_face[qp]*(1./Re*(dphi_face[i][qp]*dphi_face[j][qp])); }
											}
									
										}
									}

									// terms corresponding to minising the neumann condition
									
									if(es->parameters.get<bool>("discontinuous_linear_neumann"))	
									{
										
										for (unsigned int i=0; i<n_g_x_dofs; i++)
										{
											for (unsigned int j=0; j<n_g_x_dofs; j++)
											{
												Mbxx(i,j) += JxW_face[qp]*(dt*(xi_face[i][qp]*xi_face[j][qp]));
												Mbyy(i,j) += JxW_face[qp]*(dt*(xi_face[i][qp]*xi_face[j][qp]));
												if(threed) { Mbzz(i,j) += JxW_face[qp]*(dt*(xi_face[i][qp]*xi_face[j][qp])); }
											}
										}
										
										/*
										for (unsigned int i=0; i<n_u_dofs; i++)
										{
											for (unsigned int j=0; j<n_u_dofs; j++)
											{
												Mbxx(i,j) += JxW_face[qp]*(dt*(phi_face[i][qp]*phi_face[j][qp]));
												Mbyy(i,j) += JxW_face[qp]*(dt*(phi_face[i][qp]*phi_face[j][qp]));
												if(threed) { Mbzz(i,j) += JxW_face[qp]*(dt*(phi_face[i][qp]*phi_face[j][qp])); }
											}
										}
										*/

										for (unsigned int i=0; i<n_u_dofs; i++)
										{
											for (unsigned int j=0; j<n_g_x_dofs; j++)
											{
												Nbxx(i,j) += JxW_face[qp]*(dt*(phi_face[i][qp]*xi_face[j][qp]));
												Nbyy(i,j) += JxW_face[qp]*(dt*(phi_face[i][qp]*xi_face[j][qp]));
												if(threed) { Nbzz(i,j) += JxW_face[qp]*(dt*(phi_face[i][qp]*xi_face[j][qp])); }
											}
										}

										for (unsigned int i=0; i<n_g_x_dofs; i++)
										{
											for (unsigned int j=0; j<n_u_dofs; j++)
											{
												Nb_Txx(i,j) += JxW_face[qp]*(dt*(xi_face[i][qp]*phi_face[j][qp]));
												Nb_Tyy(i,j) += JxW_face[qp]*(dt*(xi_face[i][qp]*phi_face[j][qp]));
												if(threed) { Nb_Tzz(i,j) += JxW_face[qp]*(dt*(xi_face[i][qp]*phi_face[j][qp])); }
											}
										}
									}

									if(!convective_form || es->parameters.get<bool>("neumann_stabilised"))
									{
										//let us calculate the inflow normal velocity
										double normal_velocity = 0.0;
										double inflow_normal_velocity_bdy = 0.0;
										double outflow_normal_velocity_bdy = 0.0;
										double previous_inflow_normal_velocity_bdy = 0.0;
										//double previous_outflow_normal_velocity_bdy = 0.0;
										Number   u = 0., v = 0., w = 0.;
										Number   previous_u = 0., previous_v = 0., previous_w = 0.;
										double previous_normal_velocity = 0.0;
													
										for (unsigned int l=0; l<n_u_dofs; l++)
										{
											// From the previous Newton iterate:
											u += phi_face[l][qp]*system->current_solution (dof_indices_u[l]);
											v += phi_face[l][qp]*system->current_solution (dof_indices_v[l]);
											if(threed) { w += phi_face[l][qp]*system->current_solution (dof_indices_w[l]); }

											// From the previous Newton iterate:
											previous_u += phi_face[l][qp]*system->old_solution (dof_indices_u[l]);
											previous_v += phi_face[l][qp]*system->old_solution (dof_indices_v[l]);
											if(threed) { previous_w += phi_face[l][qp]*system->old_solution (dof_indices_w[l]); }
										}

										normal_velocity += qface_normals[qp](0) * u;
										normal_velocity += qface_normals[qp](1) * v;
										if(threed) { normal_velocity += qface_normals[qp](2) * w; }

										previous_normal_velocity += qface_normals[qp](0) * previous_u;
										previous_normal_velocity += qface_normals[qp](1) * previous_v;
										if(threed) { previous_normal_velocity += qface_normals[qp](2) * previous_w; }
										//std::cout << "normal = (" << qface_normals[qp](0) << "," << qface_normals[qp](1) << "," << qface_normals[qp](2) << ")" << std::endl;

										inflow_normal_velocity_bdy = (normal_velocity - fabs(normal_velocity))/2.0;
										outflow_normal_velocity_bdy = (normal_velocity + fabs(normal_velocity))/2.0;


										previous_inflow_normal_velocity_bdy = (previous_normal_velocity - fabs(previous_normal_velocity))/2.0;
										//previous_outflow_normal_velocity_bdy = (previous_normal_velocity + fabs(previous_normal_velocity))/2.0;

										//says whether is on the correct boundary, inflow for convective form and outflow for conservative form
										int bdy_bool = 0;
										if(fabs(normal_velocity) > 1e-10)
										{											
											if(convective_form)
												bdy_bool = ((-inflow_normal_velocity_bdy/fabs(normal_velocity)) + 0.5);
											else
												bdy_bool = (outflow_normal_velocity_bdy/fabs(normal_velocity) + 0.5);
										}

										int previous_inflow_bdy_bool = 0;
										if(fabs(previous_normal_velocity) > 1e-10)
										{											
											previous_inflow_bdy_bool = ((-previous_inflow_normal_velocity_bdy/fabs(previous_normal_velocity)) + 0.5);
										}
										//if conservative and not neumann stabilised then we always have the terms from int by parts
										if(!convective_form && !es->parameters.get<bool>("neumann_stabilised"))
										{
											bdy_bool = 1;
										}






										// hmmm small value changes appear to be making a difference
										for (unsigned int i=0; i<n_u_dofs; i++)
										{
											if(convective_form)
											{
												// if in convective form we only add a term if doing the inflow stabilisation
												if(es->parameters.get<bool>("neumann_stabilised"))
												{													
													for (unsigned int j=0; j<n_u_dofs; j++)
													{
														// - \int u_n^{in} * w * u
														Kuu(i,j) += -backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*(phi_face[i][qp]*phi_face[j][qp]);
														Kvv(i,j) += -backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*(phi_face[i][qp]*phi_face[j][qp]);
														if(threed) { Kww(i,j) += -backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*(phi_face[i][qp]*phi_face[j][qp]); }


														if(newton)
														{
															Kuu(i,j) += -backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																					* (u*qface_normals[qp](0)*phi_face[i][qp]*phi_face[j][qp]);
															Kuv(i,j) += -backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																					* (u*qface_normals[qp](1)*phi_face[i][qp]*phi_face[j][qp]);                
															if(threed) { Kuw(i,j) += -backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																					* (u*qface_normals[qp](2)*phi_face[i][qp]*phi_face[j][qp]); }
															Kvu(i,j) += -backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																					* (v*qface_normals[qp](0)*phi_face[i][qp]*phi_face[j][qp]);
															Kvv(i,j) += -backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																					* (v*qface_normals[qp](1)*phi_face[i][qp]*phi_face[j][qp]);                
															if(threed) { Kvw(i,j) += -backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																					* (v*qface_normals[qp](2)*phi_face[i][qp]*phi_face[j][qp]); }             
															if(threed) { Kwu(i,j) += -backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																					* (w*qface_normals[qp](0)*phi_face[i][qp]*phi_face[j][qp]); }              
															if(threed) { Kwv(i,j) += -backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																					* (w*qface_normals[qp](1)*phi_face[i][qp]*phi_face[j][qp]); }
															if(threed) { Kww(i,j) += -backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																					* (w*qface_normals[qp](2)*phi_face[i][qp]*phi_face[j][qp]); }
														}

													}

													// another newton term
													if(newton)
													{
														// \int (u_old \cdot w) * (u_n^{in})
														Fu(i) += -backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*u*phi_face[i][qp];
														Fv(i) += -backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*v*phi_face[i][qp];
														if(threed) { Fw(i) += -backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*w*phi_face[i][qp]; }
													}

													//we may also want to adjust the mean pressure boundary condition based on the normal velocity
													//well, of the previous timestep to be easy and avoid issues
													if(es->parameters.get<bool>("neumann_stabilised_adjusted"))
													{
														Fu(i) += -backflow_stab_param*JxW_face[qp]*dt*previous_normal_velocity*previous_u
																		 *previous_inflow_bdy_bool*phi_face[i][qp];
														Fv(i) += -backflow_stab_param*JxW_face[qp]*dt*previous_normal_velocity*previous_v
																	 	 *previous_inflow_bdy_bool*phi_face[i][qp];
														if(threed) { Fw(i) += -backflow_stab_param*JxW_face[qp]*dt*previous_normal_velocity*previous_w
																		 *previous_inflow_bdy_bool*phi_face[i][qp]; }
													}
													//Fu(i) += -JxW_face[qp]*(dt*((phi_face[i][qp]*qface_normals[qp](0))*mean_pressure));
													//Fv(i) += -JxW_face[qp]*(dt*((phi_face[i][qp]*qface_normals[qp](1))*mean_pressure));
													//Fw(i) += -JxW_face[qp]*(dt*((phi_face[i][qp]*qface_normals[qp](2))*mean_pressure));
												}
											}
											else
											{
												//always add these terms stabilised or not
												for (unsigned int j=0; j<n_u_dofs; j++)
												{
													// + \int u_n^{out} * w * u
													Kuu(i,j) += backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*(phi_face[i][qp]*phi_face[j][qp]);
													Kvv(i,j) += backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*(phi_face[i][qp]*phi_face[j][qp]);
													if(threed) { Kww(i,j) += backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*(phi_face[i][qp]*phi_face[j][qp]); }

													if(newton)
													{
														Kuu(i,j) += backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																				* (u*qface_normals[qp](0)*phi_face[i][qp]*phi_face[j][qp]);
														Kuv(i,j) += backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																				* (u*qface_normals[qp](1)*phi_face[i][qp]*phi_face[j][qp]);                
														if(threed) { Kuw(i,j) += backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																				* (u*qface_normals[qp](2)*phi_face[i][qp]*phi_face[j][qp]); }
														Kvu(i,j) += backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																				* (v*qface_normals[qp](0)*phi_face[i][qp]*phi_face[j][qp]);                
														Kvv(i,j) += backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																				* (v*qface_normals[qp](1)*phi_face[i][qp]*phi_face[j][qp]);                
														if(threed) { Kvw(i,j) += backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																				* (v*qface_normals[qp](2)*phi_face[i][qp]*phi_face[j][qp]); }           
														if(threed) { Kwu(i,j) += backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																				* (w*qface_normals[qp](0)*phi_face[i][qp]*phi_face[j][qp]); }             
														if(threed) { Kwv(i,j) += backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																				* (w*qface_normals[qp](1)*phi_face[i][qp]*phi_face[j][qp]); }
														if(threed) { Kww(i,j) += backflow_stab_param*JxW_face[qp]*dt*bdy_bool
																				* (w*qface_normals[qp](2)*phi_face[i][qp]*phi_face[j][qp]); }
													}

												}

												// another newton term
												if(newton)
												{
													// \int (u_old \cdot w) * (u_n^{in})
													Fu(i) += backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*u*phi_face[i][qp];
													Fv(i) += backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*v*phi_face[i][qp];
													if(threed) { Fw(i) += backflow_stab_param*JxW_face[qp]*dt*normal_velocity*bdy_bool*w*phi_face[i][qp]; }
												}

												//we may also want to adjust the mean pressure boundary condition based on the normal velocity
												//well, of the previous timestep to be easy and avoid issues
												if(es->parameters.get<bool>("neumann_stabilised_adjusted"))
												{
													Fu(i) += -backflow_stab_param*JxW_face[qp]*dt*previous_normal_velocity*previous_u
																		*previous_inflow_bdy_bool*phi_face[i][qp];
													Fv(i) += -backflow_stab_param*JxW_face[qp]*dt*previous_normal_velocity*previous_v
																		*previous_inflow_bdy_bool*phi_face[i][qp];
													if(threed) { Fw(i) += -backflow_stab_param*JxW_face[qp]*dt*previous_normal_velocity*previous_w
																		*previous_inflow_bdy_bool*phi_face[i][qp]; }
												}
											}
										}
									}

								}//end face quad loop



							}//end if(bc_type[boundary_id].compare("pressure") == 0)

							// mean-flow conditions/terms
							else if(bc_type[boundary_id].compare("mean-flow") == 0)
							{

								//NOT IMPLEMENTED

							}//end if(bc_type[boundary_id].compare("mean-flow") == 0)
						}

						//construct stuff to put into 1d matrix lines
						else
						{

							// mean pressure conditions/terms
							if(bc_type[boundary_id].compare("pressure") == 0)
							{
								//set the dof_index_0 to the correct node number
								dof_index_1d_0 = primary_pressure_boundary_nodes_1d[boundary_id];
								dof_index_1d_1 = secondary_pressure_boundary_nodes_1d[boundary_id];
								//dof_index_1d_2 = primary_flux_boundary_nodes_1d[boundary_id];
				
								fe_vel_face->reinit(elem, s);
								fe_pres_face->reinit(elem, s);

								// for the pressure term in the 3d equations
								for (unsigned int qp=0; qp<qface.n_points(); qp++)
								{
									for (unsigned int i=0; i<n_u_dofs; i++)
									{
								    Fu_coupled_p(i) += JxW_face[qp]*dt*phi_face[i][qp]*qface_normals[qp](0);
								    Fv_coupled_p(i) += JxW_face[qp]*dt*phi_face[i][qp]*qface_normals[qp](1);
								    if(threed) { Fw_coupled_p(i) += JxW_face[qp]*dt*phi_face[i][qp]*qface_normals[qp](2); }
									}

								}//end face quad loop

								// for the flux condition in the 1d equations

								// for the pressure term in the 3d equations
								// EQN IS + dt*Q_1d - dt*Q_3d = 0
								for (unsigned int qp=0; qp<qface.n_points(); qp++)
								{
									for (unsigned int i=0; i<n_u_dofs; i++)
									{
								    Fu_coupled_u(i) += JxW_face[qp]*dt*phi_face[i][qp]*qface_normals[qp](0);
								    Fv_coupled_u(i) += JxW_face[qp]*dt*phi_face[i][qp]*qface_normals[qp](1);
								    if(threed) { Fw_coupled_u(i) += JxW_face[qp]*dt*phi_face[i][qp]*qface_normals[qp](2); }
									}

								}//end face quad loop
							}
						}
					}

				}//end side loop


		    // Pin the pressure to zero at global node number "pressure_node".
		    // This effectively removes the non-trivial null space of constant
		    // pressure solutions.

				// need to pin pressure for lid driven cavity
		    if (pin_pressure)
		    {
		      const unsigned int pressure_node = 0;
		      const Real p_value               = 0.0;
		      for (unsigned int c=0; c<elem->n_nodes(); c++)
					{
		        if (elem->node(c) == pressure_node)
		        {
		          Kpp(c,c) += penalty;
		          Fp(c)    += penalty*p_value;
		        }
					}
		    }

		  } // end boundary condition section


			// okay now if we are doing the multiplication of the matrices for the neumann imposition
			//std::cout << "m matrix = " << Mb << std::endl;
			if(es->parameters.get<bool>("discontinuous_linear_neumann"))
			{
				if(!threed)
				{
					//okay so the structure of the matrix is known 2x2 block then diagonal then 2x2 block then diagonal
					DenseMatrix<Number> temp(6,6);
					temp = Mb;

					//std::cout << "m matrix = " << Mb << std::endl;

					double det_upper = Mb(1,1)*Mb(2,2) - Mb(1,2)*Mb(2,1);
					Mb(1,1) = temp(2,2)/det_upper;
					Mb(1,2) = -temp(1,2)/det_upper;
					Mb(2,1) = -temp(2,1)/det_upper;
					Mb(2,2) = temp(1,1)/det_upper;


					double det_lower = Mb(5,5)*Mb(6,6) - Mb(5,6)*Mb(6,5);
					Mb(5,5) = temp(6,6)/det_lower;
					Mb(5,6) = -temp(5,6)/det_lower;
					Mb(6,5) = -temp(6,5)/det_lower;
					Mb(6,6) = temp(5,5)/det_lower;

			
					//std::cout << "m matrix = " << Mb << std::endl;
					//std::cout << "det_lower = " << det_lower << std::endl;

					//std::cout << "n matrix = " << Nb << std::endl;
					//std::cout << "det_upper = " << det_upper << std::endl;

					Mb.right_multiply(Nb_T);
					Nb.right_multiply(Mb);
					for(unsigned int i=0; i<6; i++)
						for(unsigned int j=0; j<6; j++)
							Mbdy(i,j) = -Nb(i,j);

					if(fabs(det_upper) < 1e-10)
						Mbdy.zero();

					//okay let us do M_bdy = 0;
					if(!es->parameters.get<bool>("discontinuous_linear_neumann"))
						Mbdy.zero();

				}
				else
				{
					std::cout << "we cannot do this stupid matrix inversion in 3D" << std::endl;
					std::exit(0);
				}
			}

			// if we are estimating the error then calcalate relative contributions
			if(estimating_error)
			{
				std::cout << "matrix based error estimator not supported using this model." << std::endl;
				exit(0);
			}
			else
			{

				if(boundary_element)
				{
				
				}



				//now we make the extra parts of matrices from existing matrix
				// 1 - copy Bt matrix


				unsigned int from_row_offset = dim*n_u_dofs + n_p_dofs;
				unsigned int from_col_offset = dim*n_u_dofs;
				unsigned int to_row_offset = 0;
				unsigned int to_col_offset = dim*n_u_dofs + n_p_dofs + dim*n_u_dofs;
				for(unsigned int i=0; i<(dim*n_u_dofs); i++)
				{
					for(unsigned int j=0; j<n_p_dofs; j++)
					{
						Ke(i + to_row_offset, j + to_col_offset) = Ke(i + from_row_offset, j + from_col_offset);	
					}
				}
				// 2 - copy B matrix

				from_row_offset = dim*n_u_dofs + n_p_dofs + dim*n_u_dofs;
				from_col_offset = 0;
				to_row_offset = dim*n_u_dofs;
				to_col_offset = dim*n_u_dofs + n_p_dofs;
				for(unsigned int i=0; i<n_p_dofs; i++)
				{
					for(unsigned int j=0; j<(dim*n_u_dofs); j++)
					{
						Ke(i + to_row_offset, j + to_col_offset) = Ke(i + from_row_offset, j + from_col_offset);	
					}
				}
				// 3 - copy K matrix and transpose

				from_row_offset = dim*n_u_dofs + n_p_dofs;
				from_col_offset = 0;
				to_row_offset = 0;
				to_col_offset = dim*n_u_dofs + n_p_dofs;
				for(unsigned int i=0; i<(dim*n_u_dofs); i++)
				{
					for(unsigned int j=0; j<(dim*n_u_dofs); j++)
					{
						Ke(i + to_row_offset, j + to_col_offset) = Ke(j + from_row_offset, i + from_col_offset);	//note transposed
					}
				}
				// 4 - copy C matrix

				from_row_offset = dim*n_u_dofs + n_p_dofs + dim*n_u_dofs;
				from_col_offset = dim*n_u_dofs;
				to_row_offset = dim*n_u_dofs;
				to_col_offset = dim*n_u_dofs + n_p_dofs + dim*n_u_dofs;
				for(unsigned int i=0; i<(n_p_dofs); i++)
				{
					for(unsigned int j=0; j<(n_p_dofs); j++)
					{
						Ke(i + to_row_offset, j + to_col_offset) = Ke(i + from_row_offset, j + from_col_offset);
					}
				}


				//if(count == 3)
				//{
				//	std::cout << "element matrix before" << std::endl;
				//	Ke.print();
				//	exit(0);
				//}

				// apply dirichlet conditions (e.g. wall and maybe inflow)
				//hmm really not sure bout this last argument but it works
				dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices,false);
				//dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices,false);

				if(boundary_element)
				{
					//std::cout << "element matrix after, dim = " << dim << std::endl;
					//Ke.print(std::cout);
				}
				

				system->matrix->add_matrix (Ke, dof_indices);
				system->rhs->add_vector (Fe, dof_indices);
				// if we have pressure coupled, then we need and found a boundary
				if(pressure_coupled && dof_index_1d_0 != 0)
				{

					bool compute_using_monomials = false;
					// hmmm okay so we will have two columns to put here
					// a positive one on the boundary and a negative one from the same element
					for(unsigned int i=0; i< dof_indices.size(); i++)
					{
						system->matrix->add(dof_indices[i],dof_index_1d_0,Fe_coupled_p(i));
						if(compute_using_monomials)
							system->matrix->add(dof_indices[i],dof_index_1d_1,-Fe_coupled_p(i));
					}


					//flux coupling - only in row corresponding to 1d dof is how it is defined
					// this 1d dof is now pressure 0 so that getting towards symmetry
					for(unsigned int i=0; i< dof_indices.size(); i++)
					{
						system->matrix->add(dof_index_1d_0,dof_indices[i],Fe_coupled_u(i));
					}

					dof_index_1d_0 = 0;
					dof_index_1d_1 = 0;
					//dof_index_1d_2 = 0;
				}
			}





		}
  } // end of element loop


	// an element matrix
	//std::cout << "ke element matrix" << std::endl;
	//std::cout << Ke << std::endl;

/*
	std::cout << "total residual = " << total_residual << std::endl;
	std::cout << "total residual_1 = " << total_residual_1 << std::endl;
	std::cout << "total residual_2 = " << total_residual_2 << std::endl;
	std::cout << "total residual_3 = " << total_residual_3 << std::endl;
	std::cout << "total residual_4 = " << total_residual_4 << std::endl;
*/
	if(threed)
		std::cout << "End 3D assembly.";
	else
		std::cout << "End 2D assembly.";

	std::cout << "max_u = " << max_u << std::endl;

	

  return;


}
