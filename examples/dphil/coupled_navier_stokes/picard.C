
#include "picard.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

void
Picard::assemble (ErrorVector &)	// error_vector)
{

  if (threed)
    std::cout << "Begin 3D assembly... " << std::endl;
  else
    std::cout << "Begin 2D assembly... " << std::endl;

  TransientLinearImplicitSystem *system;

  // Get a reference to the Stokes system object.
  if (pressure_coupled)
    {
      system = &es->get_system < TransientLinearImplicitSystem > ("ns3d1d");
    }
  else
    {
      system = &es->get_system < TransientLinearImplicitSystem > ("ns3d");
    }

  // Problem parameters
  const Real dt = es->parameters.get < Real > ("dt");
  const Real time = system->time;
  Real Re = es->parameters.get < Real > ("reynolds_number");
  const unsigned int unsteady =
    es->parameters.get < unsigned int >("unsteady");
  bool stokes = es->parameters.get < bool > ("stokes");
  bool stab = es->parameters.get < bool > ("stab");
  Real alpha = es->parameters.get < Real > ("alpha");
  //const Real period     = es->parameters.get<Real>("period");
  bool newton = es->parameters.get < bool > ("newton");
  //const double density  = es->parameters.get<double>("density");
  const bool convective_form =
    es->parameters.get < bool > ("convective_form");
  const double backflow_stab_param =
    es->parameters.get < double >("backflow_stab_param");
  const double bertoglio_stab_param =
    es->parameters.get < double >("bertoglio_stab_param");
  bool supg = es->parameters.get < bool > ("supg");
  bool pspg = es->parameters.get < bool > ("pspg");
  bool lsic = es->parameters.get < bool > ("lsic");
  bool supg_laplacian = es->parameters.get < bool > ("supg_laplacian");
  bool symmetric_gradient =
    es->parameters.get < bool > ("symmetric_gradient");
  bool multiply_system_by_dt =
    es->parameters.get < bool > ("multiply_system_by_dt");
  double sd_param = 0.;

  //
  std::cout << "nonlin iteration = " << es->parameters.get <
    unsigned int >("nonlinear_iteration") << std::endl;
  std::cout << "max initial picard = " << es->parameters.get <
    unsigned int >("max_initial_picard") << std::endl;
  if ((es->parameters.get < unsigned int >("nonlinear_iteration") == 1
       || (es->parameters.get < unsigned int >("nonlinear_iteration") > 1
	   && es->parameters.get < double >("last_nonlinear_iterate") > es->parameters.get < double >("picard_to_newton_threshold"))) 
		&& es->parameters.get <unsigned int >("nonlinear_iteration") <= es->parameters.get <
      unsigned int >("max_initial_picard"))
    {
      std::
	cout << "Using Picard iteration because early in nonlinear iteration."
	<< std::endl;
      newton = false;
    }

  std::cout << "newton = " << newton << std::endl;

  //for the first time step we solve stokes - so that easy
  // in any case the system is stokes by default
  if (false)			//es->parameters.get<unsigned int> ("t_step") == 1)
    {
      stokes = true;
      stab = true;
      supg = false;
      pspg = false;
      lsic = false;
    }

  if (es->parameters.get < unsigned int >("t_step") == 1
      && es->parameters.get < unsigned int >("nonlinear_iteration") == 1
      && es->parameters.get < unsigned int >("num_continuation_iteration") ==
      0)
    {
      std::cout << "using stokes for first iteration" << std::endl;
      stokes = true;
    }



  // Mesh and system things
  const MeshBase & mesh = es->get_mesh ();
  const unsigned int dim = mesh.mesh_dimension ();

  const unsigned int u_var = system->variable_number ("u");
  const unsigned int v_var = system->variable_number ("v");
  unsigned int w_var = 0;
  if (threed)
    w_var = system->variable_number ("w");
  const unsigned int p_var = system->variable_number ("p");

  ForcingFunction <> forcing_function (*es);


  // Finite element types
  FEType fe_vel_type = system->variable_type (u_var);
  FEType fe_pres_type = system->variable_type (p_var);

  AutoPtr < FEBase > fe_vel (FEBase::build (dim, fe_vel_type));
  AutoPtr < FEBase > fe_pres (FEBase::build (dim, fe_pres_type));


  // Quadrature rules
  //QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  int default_quad_order =
    static_cast < int >(fe_vel_type.default_quadrature_order ());
  std::cout << "quad_order = " << default_quad_order << std::endl;
  //QGauss qrule(dim, static_cast<Order>(default_quad_order));

  QGauss qrule (dim, static_cast < Order > (default_quad_order + 2));

  //std::cout << "default quad order = " << fe_vel_type.default_quadrature_order() << std::endl;
  //std::cout << "quad order used = " << qrule.get_order() << std::endl;

  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);


  // Boundary finite element types
  AutoPtr < FEBase > fe_vel_face (FEBase::build (dim, fe_vel_type));
  AutoPtr < FEBase > fe_pres_face (FEBase::build (dim, fe_pres_type));

  //QGauss qface(dim-1, static_cast<Order>(default_quad_order+2));
  QGauss qface (dim - 1, static_cast < Order > (default_quad_order + 2));
  //QGauss qface(dim-1, static_cast<Order>(6));
  fe_vel_face->attach_quadrature_rule (&qface);
  fe_pres_face->attach_quadrature_rule (&qface);


  // Initialise some variables
  const std::vector < Real > &JxW_elem = fe_vel->get_JxW ();
  const std::vector < std::vector < Real > >&phi = fe_vel->get_phi ();
  const std::vector < std::vector < RealGradient > >&dphi =
    fe_vel->get_dphi ();
  const std::vector < std::vector < RealTensor > >&d2phi =
    fe_vel->get_d2phi ();
  const std::vector < std::vector < Real > >&psi = fe_pres->get_phi ();
  const std::vector < std::vector < RealGradient > >&dpsi =
    fe_pres->get_dphi ();
  // const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();
  const std::vector < Real > &JxW_face_elem = fe_vel_face->get_JxW ();
  const std::vector < std::vector < Real > >&phi_face =
    fe_vel_face->get_phi ();
  const std::vector < std::vector < Real > >&psi_face =
    fe_pres_face->get_phi ();
  const std::vector < std::vector < RealGradient > >&dphi_face =
    fe_vel_face->get_dphi ();
  const std::vector < Point > &qface_normals = fe_vel_face->get_normals ();
  const std::vector < std::vector < Point > >&qface_tangents =
    fe_vel_face->get_tangents ();
  const std::vector < Point > &qpoint_face = fe_vel_face->get_xyz ();
  const std::vector < Point > &qpoint = fe_vel->get_xyz ();


  // Local element matrix
  DenseMatrix < Number > Ke;
  DenseVector < Number > Fe;
  DenseVector < Number > Fe_coupled_p;
  DenseVector < Number > Fe_coupled_u;

  DenseMatrix < Number > Ke_pre_mass;
  DenseMatrix < Number > Ke_pre_laplacian;
  DenseMatrix < Number > Ke_pre_convection_diffusion;
  DenseMatrix < Number > Ke_pre_velocity_mass;
  DenseMatrix < Number > Ke_pre_velocity;	// just the velocity block, we will copy from Ke.

  DenseSubMatrix < Number >
    Kuu (Ke), Kuv (Ke), Kuw (Ke), Kup (Ke),
    Kvu (Ke), Kvv (Ke), Kvw (Ke), Kvp (Ke),
    Kwu (Ke), Kwv (Ke), Kww (Ke), Kwp (Ke),
    Kpu (Ke), Kpv (Ke), Kpw (Ke), Kpp (Ke);

  DenseSubMatrix < Number > Kpp_pre_mass (Ke_pre_mass);

  DenseSubMatrix < Number > Kpp_pre_laplacian (Ke_pre_laplacian);

  DenseSubMatrix < Number >
    Kpp_pre_convection_diffusion (Ke_pre_convection_diffusion);

  DenseSubMatrix < Number >
    Kuu_pre_velocity_mass (Ke_pre_velocity_mass),
    Kvv_pre_velocity_mass (Ke_pre_velocity_mass),
    Kww_pre_velocity_mass (Ke_pre_velocity_mass);

  DenseSubMatrix < Number >
    Kuu_pre_velocity (Ke_pre_velocity), Kuv_pre_velocity (Ke_pre_velocity),
    Kuw_pre_velocity (Ke_pre_velocity), Kup_pre_velocity (Ke_pre_velocity),
    Kvu_pre_velocity (Ke_pre_velocity), Kvv_pre_velocity (Ke_pre_velocity),
    Kvw_pre_velocity (Ke_pre_velocity), Kvp_pre_velocity (Ke_pre_velocity),
    Kwu_pre_velocity (Ke_pre_velocity), Kwv_pre_velocity (Ke_pre_velocity),
    Kww_pre_velocity (Ke_pre_velocity), Kwp_pre_velocity (Ke_pre_velocity),
    Kpu_pre_velocity (Ke_pre_velocity), Kpv_pre_velocity (Ke_pre_velocity),
    Kpw_pre_velocity (Ke_pre_velocity), Kpp_pre_velocity (Ke_pre_velocity);

  DenseSubVector < Number > Fu (Fe), Fv (Fe), Fw (Fe), Fp (Fe);

  DenseSubVector < Number > Fu_coupled_p (Fe_coupled_p), Fv_coupled_p (Fe_coupled_p), Fw_coupled_p (Fe_coupled_p), Fp_coupled_p (Fe_coupled_p);	//unnecessary


  DenseSubVector < Number > Fu_coupled_u (Fe_coupled_u), Fv_coupled_u (Fe_coupled_u), Fw_coupled_u (Fe_coupled_u), Fp_coupled_u (Fe_coupled_u);	//unnecessary






  // DofMap things
  const DofMap & dof_map = system->get_dof_map ();

  std::vector < dof_id_type > dof_indices;
  std::vector < dof_id_type > dof_indices_u;
  std::vector < dof_id_type > dof_indices_v;
  std::vector < dof_id_type > dof_indices_w;
  std::vector < dof_id_type > dof_indices_p;

  // 0 corresponds to the one actually on the boundary
  // 1 corresponds to the other dof on the element
  // let us build this map before for ease of coding
  dof_id_type pressure_1d_index = 0;	//p0
  dof_id_type flux_1d_index = 0;	//q0


  bool pin_pressure = false;
  int pressure_dof = -1;
  const unsigned int pressure_node = 0;
  if (es->parameters.get < bool > ("pin_pressure"))
    pin_pressure = true;

  // need to add dirichlet boundary conditions on pressure nodes on the inflow dirichlet boundary
  std::vector < unsigned int >pressure_dofs_on_inflow_boundary;
  std::vector < unsigned int >all_global_pressure_dofs_on_inflow_boundary;
  double sum_fp_diag = 0.;


  // Iterators
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin ();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end ();

  double max_u = 0.0;
  double max_u_old = 0.0;
  //double u_mag = 0.0;

  NumberVectorValue total_residual (0, 0);
  NumberVectorValue total_residual_1 (0, 0);
  NumberVectorValue total_residual_2 (0, 0);
  NumberVectorValue total_residual_3 (0, 0);
  NumberVectorValue total_residual_4 (0, 0);

  double tau_sum = 0.;
  double alpha_factor_sum = 0;
  unsigned int count = 0;
  for (; el != end_el; ++el)
    {
      const Elem *elem = *el;
      // only solve on the 3d subdomain
      if (std::
	  find (subdomains_3d.begin (), subdomains_3d.end (),
		elem->subdomain_id ()) != subdomains_3d.end ())
	{
	  count++;

	  // ******************************************************************* //
	  // reinitialise and set some variables
	  dof_map.dof_indices (elem, dof_indices);
	  dof_map.dof_indices (elem, dof_indices_u, u_var);
	  dof_map.dof_indices (elem, dof_indices_v, v_var);
	  if (threed)
	    dof_map.dof_indices (elem, dof_indices_w, w_var);
	  dof_map.dof_indices (elem, dof_indices_p, p_var);

	  unsigned int n_dofs = dof_indices.size ();
	  unsigned int n_u_dofs = dof_indices_u.size ();
	  unsigned int n_v_dofs = dof_indices_v.size ();
	  unsigned int n_w_dofs = 0;
	  if (threed)
	    n_w_dofs = dof_indices_w.size ();
	  unsigned int n_p_dofs = dof_indices_p.size ();

	  // Compute shape functions etc of current element
	  fe_vel->reinit (elem);
	  fe_pres->reinit (elem);
	  Real elem_volume = elem->volume ();
	  //take cube root to get diameter of the element
	  double h_T = 0;
	  if (threed)
	    {
	      h_T = pow (elem_volume, 1.0 / 3.0);
	      if (es->parameters.get < bool > ("gravemeier_element_length"))
		h_T = pow (6 * elem_volume / M_PI, 1. / 3.) / sqrt (3.);
	    }
	  else
	    h_T = pow (elem_volume, 1.0 / 2.0);

	  h_T *= es->parameters.get < double >("element_length_scaling");

	  pressure_dof = -1;
	  pressure_dofs_on_inflow_boundary.resize (0);

	  // **************************************************************** //

	  // ************* SET UP ELEMENT MATRICES ************ //

	  Ke.resize (n_dofs, n_dofs);
	  Ke_pre_mass.resize (n_dofs, n_dofs);
	  Ke_pre_laplacian.resize (n_dofs, n_dofs);
	  Ke_pre_convection_diffusion.resize (n_dofs, n_dofs);
	  Ke_pre_velocity_mass.resize (n_dofs, n_dofs);
	  Ke_pre_velocity.resize (n_dofs, n_dofs);
	  Fe.resize (n_dofs);
	  Fe_coupled_p.resize (n_dofs);
	  Fe_coupled_u.resize (n_dofs);

	  // Reposition the submatrices
	  Kuu.reposition (u_var * n_u_dofs, u_var * n_u_dofs, n_u_dofs,
			  n_u_dofs);
	  Kuv.reposition (u_var * n_u_dofs, v_var * n_u_dofs, n_u_dofs,
			  n_v_dofs);
	  Kup.reposition (u_var * n_u_dofs, p_var * n_u_dofs, n_u_dofs,
			  n_p_dofs);
	  if (threed)
	    Kuw.reposition (u_var * n_u_dofs, w_var * n_u_dofs, n_u_dofs,
			    n_w_dofs);

	  Kvu.reposition (v_var * n_v_dofs, u_var * n_v_dofs, n_v_dofs,
			  n_u_dofs);
	  Kvv.reposition (v_var * n_v_dofs, v_var * n_v_dofs, n_v_dofs,
			  n_v_dofs);
	  Kvp.reposition (v_var * n_v_dofs, p_var * n_v_dofs, n_v_dofs,
			  n_p_dofs);
	  if (threed)
	    Kvw.reposition (v_var * n_v_dofs, w_var * n_v_dofs, n_v_dofs,
			    n_w_dofs);

	  if (threed)
	    {
	      Kwu.reposition (w_var * n_w_dofs, u_var * n_w_dofs, n_w_dofs,
			      n_u_dofs);
	      Kwv.reposition (w_var * n_w_dofs, v_var * n_w_dofs, n_w_dofs,
			      n_v_dofs);
	      Kww.reposition (w_var * n_w_dofs, w_var * n_w_dofs, n_w_dofs,
			      n_w_dofs);
	      Kwp.reposition (w_var * n_w_dofs, p_var * n_w_dofs, n_w_dofs,
			      n_p_dofs);
	    }

	  Kpu.reposition (p_var * n_u_dofs, u_var * n_u_dofs, n_p_dofs,
			  n_u_dofs);
	  Kpv.reposition (p_var * n_u_dofs, v_var * n_u_dofs, n_p_dofs,
			  n_v_dofs);
	  Kpp.reposition (p_var * n_u_dofs, p_var * n_u_dofs, n_p_dofs,
			  n_p_dofs);
	  if (threed)
	    Kpw.reposition (p_var * n_u_dofs, w_var * n_u_dofs, n_p_dofs,
			    n_w_dofs);

	  Kpp_pre_mass.reposition (p_var * n_u_dofs, p_var * n_u_dofs,
				   n_p_dofs, n_p_dofs);
	  Kpp_pre_laplacian.reposition (p_var * n_u_dofs, p_var * n_u_dofs,
					n_p_dofs, n_p_dofs);
	  Kpp_pre_convection_diffusion.reposition (p_var * n_u_dofs,
						   p_var * n_u_dofs, n_p_dofs,
						   n_p_dofs);



	  Kuu_pre_velocity_mass.reposition (u_var * n_u_dofs,
					    u_var * n_u_dofs, n_u_dofs,
					    n_u_dofs);
	  Kvv_pre_velocity_mass.reposition (v_var * n_v_dofs,
					    v_var * n_v_dofs, n_v_dofs,
					    n_v_dofs);
	  if (threed)
	    Kww_pre_velocity_mass.reposition (w_var * n_w_dofs,
					      w_var * n_w_dofs, n_w_dofs,
					      n_w_dofs);



	  // the velocity block matrix
	  Kuu_pre_velocity.reposition (u_var * n_u_dofs, u_var * n_u_dofs,
				       n_u_dofs, n_u_dofs);
	  Kuv_pre_velocity.reposition (u_var * n_u_dofs, v_var * n_u_dofs,
				       n_u_dofs, n_v_dofs);
	  Kup_pre_velocity.reposition (u_var * n_u_dofs, p_var * n_u_dofs,
				       n_u_dofs, n_p_dofs);
	  if (threed)
	    Kuw_pre_velocity.reposition (u_var * n_u_dofs, w_var * n_u_dofs,
					 n_u_dofs, n_w_dofs);

	  Kvu_pre_velocity.reposition (v_var * n_v_dofs, u_var * n_v_dofs,
				       n_v_dofs, n_u_dofs);
	  Kvv_pre_velocity.reposition (v_var * n_v_dofs, v_var * n_v_dofs,
				       n_v_dofs, n_v_dofs);
	  Kvp_pre_velocity.reposition (v_var * n_v_dofs, p_var * n_v_dofs,
				       n_v_dofs, n_p_dofs);
	  if (threed)
	    Kvw_pre_velocity.reposition (v_var * n_v_dofs, w_var * n_v_dofs,
					 n_v_dofs, n_w_dofs);

	  if (threed)
	    {
	      Kwu_pre_velocity.reposition (w_var * n_w_dofs, u_var * n_w_dofs,
					   n_w_dofs, n_u_dofs);
	      Kwv_pre_velocity.reposition (w_var * n_w_dofs, v_var * n_w_dofs,
					   n_w_dofs, n_v_dofs);
	      Kww_pre_velocity.reposition (w_var * n_w_dofs, w_var * n_w_dofs,
					   n_w_dofs, n_w_dofs);
	      Kwp_pre_velocity.reposition (w_var * n_w_dofs, p_var * n_w_dofs,
					   n_w_dofs, n_p_dofs);
	    }

	  Kpu_pre_velocity.reposition (p_var * n_u_dofs, u_var * n_u_dofs,
				       n_p_dofs, n_u_dofs);
	  Kpv_pre_velocity.reposition (p_var * n_u_dofs, v_var * n_u_dofs,
				       n_p_dofs, n_v_dofs);
	  Kpp_pre_velocity.reposition (p_var * n_u_dofs, p_var * n_u_dofs,
				       n_p_dofs, n_p_dofs);
	  if (threed)
	    Kpw_pre_velocity.reposition (p_var * n_u_dofs, w_var * n_u_dofs,
					 n_p_dofs, n_w_dofs);



	  Fu.reposition (u_var * n_u_dofs, n_u_dofs);
	  Fv.reposition (v_var * n_u_dofs, n_v_dofs);
	  Fp.reposition (p_var * n_u_dofs, n_p_dofs);
	  if (threed)
	    Fw.reposition (w_var * n_u_dofs, n_w_dofs);


	  Fu_coupled_p.reposition (u_var * n_u_dofs, n_u_dofs);
	  Fv_coupled_p.reposition (v_var * n_u_dofs, n_v_dofs);
	  Fp_coupled_p.reposition (p_var * n_u_dofs, n_p_dofs);
	  if (threed)
	    Fw_coupled_p.reposition (w_var * n_u_dofs, n_w_dofs);

	  Fu_coupled_u.reposition (u_var * n_u_dofs, n_u_dofs);
	  Fv_coupled_u.reposition (v_var * n_u_dofs, n_v_dofs);
	  Fp_coupled_u.reposition (p_var * n_u_dofs, n_p_dofs);
	  if (threed)
	    Fw_coupled_u.reposition (w_var * n_u_dofs, n_w_dofs);

	  std::vector < Real > JxW = JxW_elem;
	  //if axisym then need to multiply all integrals by "r" i.e. y = qpoint[qp](1)
	  if (es->parameters.get < unsigned int >("geometry_type") == 5)
	    {
	      for (unsigned int qp = 0; qp < qrule.n_points (); qp++)
		{
		  JxW[qp] *= qpoint[qp] (1);
		}
	    }

	  if (multiply_system_by_dt)
	    {
	      for (unsigned int qp = 0; qp < qrule.n_points (); qp++)
		{
		  JxW[qp] *= dt;
		}
	    }

	  // ****************************************************************** //
	  // ***************** CALCULATE THE STREAMLINE DIFFUSION PARAM ******* //

	  double max_u_sd = 0.;
	  sd_param = 0.;
	  if (es->parameters.get < bool > ("streamline_diffusion") && !stokes)
	    {
	      // calculate the maximum u at the gauss points
	      for (unsigned int qp = 0; qp < qrule.n_points (); qp++)
		{
		  // Values to hold the solution & its gradient at the previous timestep.
		  Number u = 0., v = 0., w = 0.;

		  // Compute the velocity & its gradient from the previous timestep
		  // and the old Newton iterate.
		  for (unsigned int l = 0; l < n_u_dofs; l++)
		    {
		      // From the previous Newton iterate:
		      u +=
			phi[l][qp] *
			system->current_solution (dof_indices_u[l]);
		      v +=
			phi[l][qp] *
			system->current_solution (dof_indices_v[l]);
		      if (threed)
			w +=
			  phi[l][qp] *
			  system->current_solution (dof_indices_w[l]);
		    }

		  if (pow (pow (u, 2.0) + pow (v, 2.0) + pow (w, 2.0), 0.5) >
		      max_u_sd)
		    {
		      max_u_sd =
			pow (pow (u, 2.0) + pow (v, 2.0) + pow (w, 2.0), 0.5);
		    }
		}

	      //max_u_sd = 1.0;
	      double Pe_k =
		Re / es->parameters.get <
		double >("velocity_scale") * max_u_sd * h_T;
	      //std::cout << "hi Pe_k = " << Pe_k << std::endl;
	      if (Pe_k > 1.)
		{
		  sd_param = 0.5 * h_T * (1. - 1. / Pe_k);
		  //std::cout << "Pe_k > 1., sd_param = " << sd_param << std::endl;
		  //std::cout << "Pe_k = " << Pe_k << std::endl;
		  //std::cout << "h_T = " << h_T << std::endl;
		}
	      else
		{
		  sd_param = 0.;
		}

	      // doesn't really help
	      //sd_param *= 100.0;

	    }

	  // ***************************************************************** //








	  // ***************************************************************** //
	  // ************* CALCULATE THE STABILISATION PARAM ***************** //

	  max_u_sd = 0.;
	  if (es->parameters.get < bool > ("mesh_dependent_stab_param")
	      && (stab
		  || (es->parameters.get < unsigned int >("t_step") == 1
		      && pspg)))
	    {

	      // ******** Build volume contribution to element matrix and rhs ************** //
	      for (unsigned int qp = 0; qp < qrule.n_points (); qp++)
		{
		  // Values to hold the solution & its gradient at the previous timestep.
		  Number u = 0., v = 0., w = 0.;

		  // Compute the velocity & its gradient from the previous timestep
		  // and the old Newton iterate.
		  for (unsigned int l = 0; l < n_u_dofs; l++)
		    {
		      // From the previous Newton iterate:
		      u +=
			phi[l][qp] *
			system->current_solution (dof_indices_u[l]);
		      v +=
			phi[l][qp] *
			system->current_solution (dof_indices_v[l]);
		      if (threed)
			w +=
			  phi[l][qp] *
			  system->current_solution (dof_indices_w[l]);
		    }

		  if (pow (pow (u, 2.0) + pow (v, 2.0) + pow (w, 2.0), 0.5) >
		      max_u_sd)
		    {
		      max_u_sd =
			pow (pow (u, 2.0) + pow (v, 2.0) + pow (w, 2.0), 0.5);
		    }
		}

	      if (stokes)
		max_u_sd = 0.;

	      double m_k = 1. / 3.;
	      double Pe_k_1 = 4. * dt / (Re * m_k * h_T * h_T);
	      double Pe_k_2 = Re * m_k * max_u_sd * h_T / 2.;
	      double xi_1 = Pe_k_1;	//will always be the max cause dt is "infinite"
	      if (unsteady)
		xi_1 = std::max (Pe_k_1, 1.);
	      double xi_2 = std::max (Pe_k_2, 1.);

	      // take out factor of h_T
	      alpha =
		1. / (1. / dt * xi_1 + 4. / m_k / Re / (h_T * h_T) * xi_2);

	      //std::cout << " using stab_param alpha = " << alpha << std::endl;                      

	    }
	  else
	    {
	      alpha = es->parameters.get < Real > ("alpha") * Re * h_T * h_T;
	    }

	  // debugging
	  alpha_factor_sum += alpha / (Re * h_T * h_T);

	  // *************************************************************************** //




	  // *************************************************************************** //
	  // ******** BUILD VOLUME CONTRIBUTION TO ELEMENT MATRIX AND RHS ************* //
	  for (unsigned int qp = 0; qp < qrule.n_points (); qp++)
	    {




	      // ************* SET UP VARIABLES FROM PREVIOUS ITERATE/TIMESTEP ************* //

	      // Values to hold the solution & its gradient at the previous timestep.
	      Number u = 0., u_old = 0., v = 0., v_old = 0., w = 0., w_old =
		0.;
	      Gradient grad_u, grad_v, grad_w, grad_p, grad_u_old, grad_v_old,
		grad_w_old;
	      Number lap_u = 0., lap_v = 0., lap_w = 0.;

	      // Compute the velocity & its gradient from the previous timestep
	      // and the old Newton iterate.
	      for (unsigned int l = 0; l < n_u_dofs; l++)
		{
		  // From the old timestep:
		  u_old +=
		    phi[l][qp] * system->old_solution (dof_indices_u[l]);
		  v_old +=
		    phi[l][qp] * system->old_solution (dof_indices_v[l]);
		  if (threed)
		    w_old +=
		      phi[l][qp] * system->old_solution (dof_indices_w[l]);


		  // From the previous Newton iterate:
		  u +=
		    phi[l][qp] * system->current_solution (dof_indices_u[l]);
		  v +=
		    phi[l][qp] * system->current_solution (dof_indices_v[l]);
		  if (threed)
		    w +=
		      phi[l][qp] *
		      system->current_solution (dof_indices_w[l]);

		  grad_u.add_scaled (dphi[l][qp],
				     system->
				     current_solution (dof_indices_u[l]));
		  grad_v.add_scaled (dphi[l][qp],
				     system->
				     current_solution (dof_indices_v[l]));
		  if (threed)
		    grad_w.add_scaled (dphi[l][qp],
				       system->
				       current_solution (dof_indices_w[l]));


		  grad_u_old.add_scaled (dphi[l][qp],
					 system->
					 old_solution (dof_indices_u[l]));
		  grad_v_old.add_scaled (dphi[l][qp],
					 system->
					 old_solution (dof_indices_v[l]));
		  if (threed)
		    grad_w_old.add_scaled (dphi[l][qp],
					   system->
					   old_solution (dof_indices_w[l]));

		  double laplacian_operator = 0;
		  laplacian_operator +=
		    d2phi[l][qp] (0, 0) + d2phi[l][qp] (1, 1);
		  if (threed)
		    {
		      laplacian_operator += d2phi[l][qp] (2, 2);
		    }

		  lap_u +=
		    laplacian_operator *
		    system->current_solution (dof_indices_u[l]);
		  lap_v +=
		    laplacian_operator *
		    system->current_solution (dof_indices_v[l]);
		  if (threed)
		    lap_w +=
		      laplacian_operator *
		      system->current_solution (dof_indices_w[l]);

		}

	      for (unsigned int l = 0; l < n_p_dofs; l++)
		{
		  grad_p.add_scaled (dpsi[l][qp],
				     system->
				     current_solution (dof_indices_p[l]));
		}

	      // vector forms of velocities
	      NumberVectorValue U_old;
	      if (threed)
		U_old = NumberVectorValue (u_old, v_old, w_old);
	      else
		U_old = NumberVectorValue (u_old, v_old);


	      NumberVectorValue U;
	      if (threed)
		U = NumberVectorValue (u, v, w);
	      else
		U = NumberVectorValue (u, v);



	      // SOME DEBUGGING

	      if (pow (pow (u, 2.0) + pow (v, 2.0), 0.5) > max_u)
		max_u = pow (pow (u, 2.0) + pow (v, 2.0), 0.5);

	      if (pow (pow (u_old, 2.0) + pow (v_old, 2.0), 0.5) > max_u_old)
		max_u_old = pow (pow (u_old, 2.0) + pow (v_old, 2.0), 0.5);

	      /*
	         if(threed)
	         u_mag = pow(pow(u,2.0) + pow(v,2.0) + pow(w,2.0),0.5);
	         else
	         u_mag = pow(pow(u,2.0) + pow(v,2.0),0.5);
	       */


	      //if(fabs(u_old) > max_u)
	      //      max_u =fabs(u_old);



	      // *************************************************************** //










	      // *************************************************************************** //
	      // ********************* REGULAR NAVIER STOKES TERMS ************************* //

	      // First, an i-loop over the velocity degrees of freedom.
	      for (unsigned int i = 0; i < n_u_dofs; i++)
		{

		  // *********** Unsteady term ************* //
		  if (unsteady)
		    {
		      Fu (i) -= -JxW[qp] / dt * (u_old * phi[i][qp]);
		      Fv (i) -= -JxW[qp] / dt * (v_old * phi[i][qp]);
		      if (threed)
			{
			  Fw (i) -= -JxW[qp] / dt * (w_old * phi[i][qp]);
			}
		    }





		  // ********** Forcing term ************** //
		  DenseVector < Number > f;
		  f.resize (dim);
		  forcing_function (qpoint[qp], time, f);
		  //std::cout << "f[" << qpoint[qp] << "] = " << f;
		  //if(unsteady)
		  //{
		  //Fu(i) += JxW[qp]*(phi[i][qp]*f(0));
		  //Fv(i) += JxW[qp]*(phi[i][qp]*f(1));
		  //if(threed) {  Fw(i) += JxW[qp]*(phi[i][qp]*f(2)); }
		  //}




		  // ********** Newton convection term ************ //
		  if (!stokes && newton)
		    {
		      if (convective_form)
			{
			  Fu (i) += JxW[qp] * ((U * grad_u) * phi[i][qp]);
			  Fv (i) += JxW[qp] * ((U * grad_v) * phi[i][qp]);
			  if (threed)
			    {
			      Fw (i) += JxW[qp] * ((U * grad_w) * phi[i][qp]);
			    }
			}
		      else
			{
			  Fu (i) += -JxW[qp] * (U * dphi[i][qp]) * u;
			  Fv (i) += -JxW[qp] * (U * dphi[i][qp]) * v;
			  if (threed)
			    {
			      Fw (i) += -JxW[qp] * (U * dphi[i][qp]) * w;
			    }
			}
		    }




		  // Matrix contributions for the uu and vv couplings. i.e. first equation
		  for (unsigned int j = 0; j < n_u_dofs; j++)
		    {



		      // *************** Unsteady term i.e. mass matrix ************ //
		      if (unsteady)
			{
			  Kuu (i, j) += JxW[qp] / dt * (phi[i][qp] * phi[j][qp]);	// mass matrix term
			  Kvv (i, j) += JxW[qp] / dt * (phi[i][qp] * phi[j][qp]);	// mass matrix term
			  if (threed)
			    {
			      Kww (i, j) +=
				JxW[qp] / dt * (phi[i][qp] * phi[j][qp]);
			    }	// mass matrix term
			}






		      // *************** Diffusion terms ************** //
		      // We can use the nonsymmetric gradient or the symmetrised gradient
		      if (!symmetric_gradient)
			{
			  Kuu (i, j) += JxW[qp] * (1. / Re * (dphi[i][qp] * dphi[j][qp]));	// diffusion
			  Kvv (i, j) += JxW[qp] * (1. / Re * (dphi[i][qp] * dphi[j][qp]));	// diffusion
			  if (threed)
			    {
			      Kww (i, j) +=
				JxW[qp] * (1. / Re *
					   (dphi[i][qp] * dphi[j][qp]));
			    }	// diffusion
			}
		      else
			{
			  Kuu (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] * dphi[j][qp]));	// diffusion
			  Kvv (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] * dphi[j][qp]));	// diffusion
			  if (threed)
			    {
			      Kww (i, j) +=
				0.5 * JxW[qp] * (1. / Re *
						 (dphi[i][qp] * dphi[j][qp]));
			    }	// diffusion

			  Kuu (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (0) * dphi[j][qp] (0)));	// diffusion
			  Kuv (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (0) * dphi[j][qp] (1)));	// diffusion
			  if (threed)
			    {
			      Kuw (i, j) +=
				0.5 * JxW[qp] * (1. / Re *
						 (dphi[i][qp] (0) *
						  dphi[j][qp] (2)));
			    }	// diffusion
			  Kvu (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (1) * dphi[j][qp] (0)));	// diffusion
			  Kvv (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (1) * dphi[j][qp] (1)));	// diffusion
			  if (threed)
			    {
			      Kvw (i, j) +=
				0.5 * JxW[qp] * (1. / Re *
						 (dphi[i][qp] (1) *
						  dphi[j][qp] (2)));
			    }	// diffusion
			  Kwu (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (2) * dphi[j][qp] (0)));	// diffusion
			  Kwv (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (2) * dphi[j][qp] (1)));	// diffusion
			  if (threed)
			    {
			      Kww (i, j) +=
				0.5 * JxW[qp] * (1. / Re *
						 (dphi[i][qp] (2) *
						  dphi[j][qp] (2)));
			    }	// diffusion

			}







		      // *************** Convection terms *************** //
		      if (!stokes)
			{

			  // use convective or conservative form 
			  if (convective_form)
			    {
			      Kuu (i, j) += JxW[qp] * ((U * dphi[j][qp]) * phi[i][qp]);	// convection
			      Kvv (i, j) += JxW[qp] * ((U * dphi[j][qp]) * phi[i][qp]);	// convection
			      if (threed)
				{
				  Kww (i, j) +=
				    JxW[qp] * ((U * dphi[j][qp]) *
					       phi[i][qp]);
				}	// convection
			    }
			  else
			    {
			      Kuu (i, j) += -JxW[qp] * ((U * dphi[i][qp]) * phi[j][qp]);	// convection
			      Kvv (i, j) += -JxW[qp] * ((U * dphi[i][qp]) * phi[j][qp]);	// convection
			      if (threed)
				{
				  Kww (i, j) +=
				    -JxW[qp] * ((U * dphi[i][qp]) *
						phi[j][qp]);
				}	// convection
			    }



			  // *************** Streamline Diffusion Terms ************** //
			  // streamline diffusion term, this is nonlinear and done in a 
			  // picard sense so may hamper convergence.
			  if (es->parameters.get < bool >
			      ("streamline_diffusion") && !stokes)
			    {
			      if (!convective_form)
				{
				  std::
				    cout <<
				    "error conservative form for streamline diffusion not implemented"
				    << std::endl;
				  std::exit (0);
				}

			      Kuu (i, j) += sd_param * JxW[qp] * ((U * dphi[i][qp]) * (U * dphi[j][qp]));	// convection
			      Kvv (i, j) += sd_param * JxW[qp] * ((U * dphi[i][qp]) * (U * dphi[j][qp]));	// convection
			      if (threed)
				{
				  Kww (i, j) +=
				    sd_param * JxW[qp] * ((U * dphi[i][qp]) *
							  (U * dphi[j][qp]));
				}	// convection

			    }





			  // **************** Newton Convection Terms ***************** //
			  // convective or conservative form.
			  if (newton)
			    {
			      if (convective_form)
				{
				  Kuu (i, j) +=
				    JxW[qp] * (grad_u (0) * phi[i][qp] *
					       phi[j][qp]);
				  Kuv (i, j) +=
				    JxW[qp] * (grad_u (1) * phi[i][qp] *
					       phi[j][qp]);
				  if (threed)
				    {
				      Kuw (i, j) +=
					JxW[qp] * (grad_u (2) * phi[i][qp] *
						   phi[j][qp]);
				    }
				  Kvu (i, j) +=
				    JxW[qp] * (grad_v (0) * phi[i][qp] *
					       phi[j][qp]);
				  Kvv (i, j) +=
				    JxW[qp] * (grad_v (1) * phi[i][qp] *
					       phi[j][qp]);
				  if (threed)
				    {
				      Kvw (i, j) +=
					JxW[qp] * (grad_v (2) * phi[i][qp] *
						   phi[j][qp]);
				    }
				  if (threed)
				    {
				      Kwu (i, j) +=
					JxW[qp] * (grad_w (0) * phi[i][qp] *
						   phi[j][qp]);
				    }
				  if (threed)
				    {
				      Kwv (i, j) +=
					JxW[qp] * (grad_w (1) * phi[i][qp] *
						   phi[j][qp]);
				    }
				  if (threed)
				    {
				      Kww (i, j) +=
					JxW[qp] * (grad_w (2) * phi[i][qp] *
						   phi[j][qp]);
				    }
				}
			      else
				{
				  Kuu (i, j) +=
				    -JxW[qp] * (u * dphi[i][qp] (0) *
						phi[j][qp]);
				  Kuv (i, j) +=
				    -JxW[qp] * (u * dphi[i][qp] (1) *
						phi[j][qp]);
				  if (threed)
				    {
				      Kuw (i, j) +=
					-JxW[qp] * (u * dphi[i][qp] (2) *
						    phi[j][qp]);
				    }
				  Kvu (i, j) +=
				    -JxW[qp] * (v * dphi[i][qp] (0) *
						phi[j][qp]);
				  Kvv (i, j) +=
				    -JxW[qp] * (v * dphi[i][qp] (1) *
						phi[j][qp]);
				  if (threed)
				    {
				      Kvw (i, j) +=
					-JxW[qp] * (v * dphi[i][qp] (2) *
						    phi[j][qp]);
				    }
				  if (threed)
				    {
				      Kwu (i, j) +=
					-JxW[qp] * (w * dphi[i][qp] (0) *
						    phi[j][qp]);
				    }
				  if (threed)
				    {
				      Kwv (i, j) +=
					-JxW[qp] * (w * dphi[i][qp] (1) *
						    phi[j][qp]);
				    }
				  if (threed)
				    {
				      Kww (i, j) +=
					-JxW[qp] * (w * dphi[i][qp] (2) *
						    phi[j][qp]);
				    }
				}
			    }
			}

		    }







		  // ************ Pressure gradient term *********************** //
		  // put the density in here
		  for (unsigned int j = 0; j < n_p_dofs; j++)
		    {
		      Kup (i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp] (0);
		      Kvp (i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp] (1);
		      if (threed)
			{
			  Kwp (i, j) +=
			    -JxW[qp] * psi[j][qp] * dphi[i][qp] (2);
			}
		    }
		}






	      for (unsigned int i = 0; i < n_p_dofs; i++)
		{
		  // ************* Continuity equation  ****************************** //
		  // could possibly be NOT multiplied by negative.
		  for (unsigned int j = 0; j < n_u_dofs; j++)
		    {
		      Kpu (i, j) += -JxW[qp] * psi[i][qp] * dphi[j][qp] (0);
		      Kpv (i, j) += -JxW[qp] * psi[i][qp] * dphi[j][qp] (1);
		      if (threed)
			{
			  Kpw (i, j) +=
			    -JxW[qp] * psi[i][qp] * dphi[j][qp] (2);
			}
		    }







		  // ************** Stabilisation term - should be opposite sign to continuity equation
		  // stabilisation depends on if using P1-P1 or P1-P0
		  // need stabilisation for first tstep of pspg because parameter is zero
		  if (stab
		      || (es->parameters.get < unsigned int >("t_step") == 1
			  && pspg))
		    {
		      for (unsigned int j = 0; j < n_p_dofs; j++)
			{
			  //Kpp(i,j) += alpha*elem_volume*JxW[qp]*dpsi[i][qp]*dpsi[j][qp];
			  Kpp (i, j) +=
			    -alpha * JxW[qp] * dpsi[i][qp] * dpsi[j][qp];
			}
		    }








		  // ***************** Preconditioner terms ************************ //
		  for (unsigned int j = 0; j < n_p_dofs; j++)
		    {

		      Kpp_pre_mass (i, j) -=
			JxW[qp] * psi[i][qp] * psi[j][qp];
		      Kpp_pre_laplacian (i, j) +=
			JxW[qp] * dpsi[i][qp] * dpsi[j][qp];
		      // I don't think the convection diffusion operator should have a 
		      // mass matrix corresponding to the unsteady terms.
		      if (!stokes)
			Kpp_pre_convection_diffusion (i, j) +=
			  JxW[qp] * (U * dpsi[j][qp]) * psi[i][qp];

		      Kpp_pre_convection_diffusion (i, j) +=
			JxW[qp] / Re * dpsi[i][qp] * dpsi[j][qp];

		      if (unsteady)
			Kpp_pre_convection_diffusion (i, j) +=
			  JxW[qp] / dt * psi[i][qp] * psi[j][qp];


		      if (es->parameters.get < bool > ("streamline_diffusion")
			  && !stokes)
			Kpp_pre_convection_diffusion (i, j) +=
			  sd_param * JxW[qp] * (U * dpsi[i][qp]) * (U *
								    dpsi[j]
								    [qp]);

		    }
		}

	      // more preconditioner terms
	      for (unsigned int i = 0; i < n_u_dofs; i++)
		{
		  for (unsigned int j = 0; j < n_u_dofs; j++)
		    {
		      Kuu_pre_velocity_mass (i, j) +=
			JxW[qp] * phi[i][qp] * phi[j][qp];
		      Kvv_pre_velocity_mass (i, j) +=
			JxW[qp] * phi[i][qp] * phi[j][qp];
		      Kww_pre_velocity_mass (i, j) +=
			JxW[qp] * phi[i][qp] * phi[j][qp];
		    }
		}





	      // temp more preconditioner terms

	      for (unsigned int i = 0; i < n_u_dofs; i++)
		{
		  for (unsigned int j = 0; j < n_u_dofs; j++)
		    {

		      // *************** Unsteady term i.e. mass matrix ************ //
		      if (unsteady)
			{
			  Kuu_pre_velocity (i, j) += JxW[qp] / dt * (phi[i][qp] * phi[j][qp]);	// mass matrix term
			  Kvv_pre_velocity (i, j) += JxW[qp] / dt * (phi[i][qp] * phi[j][qp]);	// mass matrix term
			  if (threed)
			    {
			      Kww_pre_velocity (i, j) +=
				JxW[qp] / dt * (phi[i][qp] * phi[j][qp]);
			    }	// mass matrix term
			}






		      // *************** Diffusion terms ************** //
		      // We can use the nonsymmetric gradient or the symmetrised gradient
		      if (!symmetric_gradient)
			{
			  Kuu_pre_velocity (i, j) += JxW[qp] * (1. / Re * (dphi[i][qp] * dphi[j][qp]));	// diffusion
			  Kvv_pre_velocity (i, j) += JxW[qp] * (1. / Re * (dphi[i][qp] * dphi[j][qp]));	// diffusion
			  if (threed)
			    {
			      Kww_pre_velocity (i, j) +=
				JxW[qp] * (1. / Re *
					   (dphi[i][qp] * dphi[j][qp]));
			    }	// diffusion
			}
		      else
			{
			  Kuu_pre_velocity (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] * dphi[j][qp]));	// diffusion
			  Kvv_pre_velocity (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] * dphi[j][qp]));	// diffusion
			  if (threed)
			    {
			      Kww_pre_velocity (i, j) +=
				0.5 * JxW[qp] * (1. / Re *
						 (dphi[i][qp] * dphi[j][qp]));
			    }	// diffusion

			  Kuu_pre_velocity (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (0) * dphi[j][qp] (0)));	// diffusion
			  Kuv_pre_velocity (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (0) * dphi[j][qp] (1)));	// diffusion
			  if (threed)
			    {
			      Kuw_pre_velocity (i, j) +=
				0.5 * JxW[qp] * (1. / Re *
						 (dphi[i][qp] (0) *
						  dphi[j][qp] (2)));
			    }	// diffusion
			  Kvu_pre_velocity (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (1) * dphi[j][qp] (0)));	// diffusion
			  Kvv_pre_velocity (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (1) * dphi[j][qp] (1)));	// diffusion
			  if (threed)
			    {
			      Kvw_pre_velocity (i, j) +=
				0.5 * JxW[qp] * (1. / Re *
						 (dphi[i][qp] (1) *
						  dphi[j][qp] (2)));
			    }	// diffusion
			  Kwu_pre_velocity (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (2) * dphi[j][qp] (0)));	// diffusion
			  Kwv_pre_velocity (i, j) += 0.5 * JxW[qp] * (1. / Re * (dphi[i][qp] (2) * dphi[j][qp] (1)));	// diffusion
			  if (threed)
			    {
			      Kww_pre_velocity (i, j) +=
				0.5 * JxW[qp] * (1. / Re *
						 (dphi[i][qp] (2) *
						  dphi[j][qp] (2)));
			    }	// diffusion

			}


		    }

		  // ************ Pressure gradient term *********************** //
		  // put the density in here
		  for (unsigned int j = 0; j < n_p_dofs; j++)
		    {
		      Kup_pre_velocity (i, j) +=
			-JxW[qp] * psi[j][qp] * dphi[i][qp] (0);
		      Kvp_pre_velocity (i, j) +=
			-JxW[qp] * psi[j][qp] * dphi[i][qp] (1);
		      if (threed)
			{
			  Kwp_pre_velocity (i, j) +=
			    -JxW[qp] * psi[j][qp] * dphi[i][qp] (2);
			}
		    }
		}


	      for (unsigned int i = 0; i < n_p_dofs; i++)
		{
		  // ************* Continuity equation  ****************************** //
		  // could possibly be NOT multiplied by negative.
		  for (unsigned int j = 0; j < n_u_dofs; j++)
		    {
		      Kpu_pre_velocity (i, j) +=
			-JxW[qp] * psi[i][qp] * dphi[j][qp] (0);
		      Kpv_pre_velocity (i, j) +=
			-JxW[qp] * psi[i][qp] * dphi[j][qp] (1);
		      if (threed)
			{
			  Kpw_pre_velocity (i, j) +=
			    -JxW[qp] * psi[i][qp] * dphi[j][qp] (2);
			}
		    }
		  // ************** Stabilisation term - should be opposite sign to continuity equation
		  // stabilisation depends on if using P1-P1 or P1-P0
		  // need stabilisation for first tstep of pspg because parameter is zero
		  if (stab
		      || (es->parameters.get < unsigned int >("t_step") == 1
			  && pspg))
		    {
		      for (unsigned int j = 0; j < n_p_dofs; j++)
			{
			  Kpp_pre_velocity (i, j) +=
			    -alpha * JxW[qp] * dpsi[i][qp] * dpsi[j][qp];
			}
		    }
		}





	      // *************************************************************************** //











	      // ************************** SUPG/PSPG/LSIC TERMS *************************** //
	      // we only want to apply this if there is a nonzero velocity

	      if (supg || pspg || lsic)
		{
		  double length_scale =
		    es->parameters.get < double >("length_scale");
		  double velocity_scale =
		    es->parameters.get < double >("velocity_scale");
		  double time_scale = length_scale / velocity_scale;
		  double mu = es->parameters.get < double >("viscosity");
		  double rho = es->parameters.get < double >("density");

		  // petrov galerkin velocity - if unsteady use previous timestep ala tedzuyar, otherwise use current (picard)
		  NumberVectorValue velocity;	//velocity to use for supg - could be previous timestep most likely

		  NumberVectorValue grad_V_u;
		  NumberVectorValue grad_V_v;
		  NumberVectorValue grad_V_w;
		  // not we can't use ismail parameters when steady sim
		  // obviousl can't be constant on the first time step, oh well this is to complex now
		  if (unsteady
		      && es->parameters.get < bool >
		      ("supg_constant_constant"))
		    {
		      velocity = U_old;
		      grad_V_u = grad_u_old;
		      grad_V_v = grad_v_old;
		      grad_V_w = grad_w_old;
		      //velocity = NumberVectorValue(1., 0.);
		    }
		  else
		    {
		      velocity = U;
		      grad_V_u = grad_u;
		      grad_V_v = grad_v;
		      grad_V_w = grad_w;
		    }



		  //if(es->parameters.get<bool>("supg_full_newton"))
		  //      velocity = U;

		  NumberVectorValue grad_Vx;
		  NumberVectorValue grad_Vy;
		  NumberVectorValue grad_Vz;
		  NumberVectorValue r;
		  if (threed)
		    {
		      grad_Vx =
			NumberVectorValue (grad_V_u (0), grad_V_v (0),
					   grad_V_w (0));
		      grad_Vy =
			NumberVectorValue (grad_V_u (1), grad_V_v (1),
					   grad_V_w (1));
		      grad_Vz =
			NumberVectorValue (grad_V_u (2), grad_V_v (2),
					   grad_V_w (2));
		      r =
			NumberVectorValue (grad_Vx * velocity,
					   grad_Vy * velocity,
					   grad_Vz * velocity);
		      r = r.unit ();
		    }
		  else
		    {
		      grad_Vx =
			NumberVectorValue (grad_V_u (0), grad_V_v (0));
		      grad_Vy =
			NumberVectorValue (grad_V_u (1), grad_V_v (1));
		      r =
			NumberVectorValue (grad_Vx * velocity,
					   grad_Vy * velocity);
		      r = r.unit ();
		    }

		  //velocity = U;         // always use from previous time step?

		  double velocity_mag = velocity.size ();

		  // need this to calc the local reynolds number
		  double grad_times_r = 0.;
		  for (unsigned int i = 0; i < n_u_dofs; i++)
		    {
		      grad_times_r += fabs (dphi[i][qp] * r);	//r is already unit /length_scale for the grad
		    }
		  double h_rgn = 2. / grad_times_r;

		  double reynolds_supg =
		    velocity_mag * h_rgn / (2 * mu / rho) * velocity_scale *
		    length_scale;

		  //supg parameter
		  double tau_supg =
		    1.0 / sqrt (pow (2 * velocity_mag / h_T, 2.0) +
				9.0 * pow (4 / (h_T * h_T * reynolds_supg),
					   2.0));
		  //std::cout << "tau_supg = " << tau_supg << std::endl;                                  
		  //tau_supg /= 2.0;
		  //tau_supg = h_T*h_T;

		  double tau_pspg = 0.;
		  double tau_lsic = 0.;


		  // note the supg parameter is not affected by the nondimensionalisation so need to redimensionalise
		  // including the length scales

		  if (velocity_mag < 1e-10)
		    {
		      tau_supg = 0.;
		      tau_pspg = 0.;
		      tau_lsic = 0.;
		    }
		  else
		    {
		      if (es->parameters.get < unsigned int >("supg_parameter") == 0)	//not great nonlinear conv but better
			{
			  //static characteristic element length
			  double h_p =
			    pow (6 * elem_volume / pi,
				 1.0 / 3.0) / sqrt (3.0);
			  double grad_times_velocity = 0.;
			  double grad_times_r = 0.;


			  // dimensionalised already
			  for (unsigned int i = 0; i < n_u_dofs; i++)
			    {
			      grad_times_velocity += dphi[i][qp] * velocity;	///length_scale for the grad
			    }
			  double h_u = 2. / grad_times_velocity;

			  // in this method, we definitely use the actual reynolds number we use the calc, cause comes from viscosity
			  double m_k = 1. / 3.;	//for linearly interpolated elements
			  double xi_C =
			    std::min (m_k * velocity_mag * h_p * Re / 2., 1.);
			  double tau_C = velocity_mag * h_p / 2 * xi_C;

			  double f_T = 1.0;	// depends on time integration scheme
			  double xi_M1p =
			    std::max (4 * f_T * dt / (Re * m_k * h_p * h_p),
				      1.);
			  double xi_M2p =
			    std::max (m_k * velocity_mag * h_p * Re / 2., 1.);
			  double tau_Mp =
			    1. / (1. / (f_T * dt) * xi_M1p +
				  4. / (Re * m_k * h_p * h_p) * xi_M2p);

			  double xi_M1u =
			    std::max (4 * f_T * dt / (Re * m_k * h_u * h_u),
				      1.);
			  double xi_M2u =
			    std::max (m_k * velocity_mag * h_u * Re / 2., 1.);
			  double tau_Mu =
			    1. / (1. / (f_T * dt) * xi_M1u +
				  4. / (Re * m_k * h_p * h_p) * xi_M2u);

			  //std::cout << "tau old = " << tau_supg << std::endl;
			  tau_supg = tau_Mu;
			  //std::cout << "tau new = " << tau_supg << std::endl;
			  tau_lsic = tau_C;
			  tau_pspg = tau_Mp;
			}
		      else if (es->parameters.get < unsigned int >("supg_parameter") == 1)	// slow convergence (none with changing)
			{
			  //static characteristic element length
			  //double h_p = pow(6*elem_volume/pi,1.0/3.0)/sqrt(3.0);
			  double grad_times_velocity = 0.;
			  double grad_times_r = 0.;

			  for (unsigned int i = 0; i < n_u_dofs; i++)
			    {
			      grad_times_velocity +=
				fabs (dphi[i][qp] * velocity);
			    }
			  double tau_sugn12 = 1. / grad_times_velocity;

			  double tau_sugn3 =
			    h_rgn * h_rgn * reynolds_supg / 4.;

			  tau_supg =
			    1. / sqrt (1. / (tau_sugn12 * tau_sugn12) +
				       1. / (tau_sugn3 * tau_sugn3));
			  tau_pspg = tau_supg;
			  tau_lsic = tau_supg * velocity_mag * velocity_mag;

			}
		      else if (es->parameters.get < unsigned int >("supg_parameter") == 2)	//good convecgence
			{
			  double grad_times_velocity = 0.;
			  double grad_times_r = 0.;

			  for (unsigned int i = 0; i < n_u_dofs; i++)
			    {
			      grad_times_velocity +=
				fabs (dphi[i][qp] * velocity);
			    }
			  double h_ugn =
			    2. * velocity_mag / grad_times_velocity;
			  double tau_sugn1 = h_ugn / 2. / velocity_mag;
			  double tau_sugn2 = dt / 2.;
			  double tau_sugn3 =
			    h_ugn * h_ugn * reynolds_supg / 4.;

			  double Re_ugn =
			    velocity_mag * h_ugn * reynolds_supg / 2.;
			  double z = 0.;
			  if (Re_ugn <= 3.)
			    z = Re_ugn / 3.;
			  else
			    z = 1.;

			  tau_supg =
			    1. / sqrt (1. / (tau_sugn1 * tau_sugn1) +
				       1. / (tau_sugn2 * tau_sugn2) +
				       1. / (tau_sugn3 * tau_sugn3));
			  tau_pspg = tau_supg;
			  tau_lsic = h_ugn / 2.0 * velocity_mag * z;

			}
		      else
			{
			  tau_pspg = tau_supg;

			}
		    }



		  tau_supg *= es->parameters.get < double >("supg_scale");
		  tau_pspg *= es->parameters.get < double >("supg_scale");
		  tau_lsic *= es->parameters.get < double >("supg_scale");


		  tau_sum += tau_supg;

		  /*
		     //calculate tau_supg on an element level basis
		     DenseMatrix m, c, k, k_tilde, c_tilde;
		     m.resize(n_dofs,n_u_dofs);
		     c.resize(n_dofs,n_u_dofs);
		     k.resize(n_dofs,n_u_dofs);
		     k_tilde.resize(n_dofs,n_u_dofs);
		     c_tilde.resize(n_dofs,n_u_dofs);

		     double tau_advection = 0.;
		     double tau_transient = 0.;
		     double tau_diffusion = 0.;

		     for (unsigned int i=0; i<n_u_dofs; i++)
		     {
		     for (unsigned int j=0; j<n_u_dofs; j++)
		     {
		     m(i,j) +=phi
		     }
		     }
		   */



		  // now that we have set up the parameter, we can use s different velocity for the picard/newton iterations
		  // not we can't use ismail parameters when steady sim
		  if (false)	//unsteady)
		    {
		      velocity = U_old;
		      grad_V_u = grad_u_old;
		      grad_V_v = grad_v_old;
		      grad_V_w = grad_w_old;
		      //velocity = NumberVectorValue(1., 0.);
		    }
		  else
		    {
		      velocity = U;
		      grad_V_u = grad_u;
		      grad_V_v = grad_v;
		      grad_V_w = grad_w;
		    }


		  // let us use the velocity from the previous
		  // there is a dt from the actual time derivative as well as from inside the residual
		  if (supg)
		    {
		      if (es->parameters.get < bool > ("supg_newton"))
			{


			  for (unsigned int i = 0; i < n_u_dofs; i++)
			    {

			      if (unsteady)
				{
				  Fu (i) -=
				    -tau_supg * JxW[qp] / dt * (velocity *
								dphi[i][qp] *
								u_old);
				  Fv (i) -=
				    -tau_supg * JxW[qp] / dt * (velocity *
								dphi[i][qp] *
								v_old);
				  if (threed)
				    {
				      Fw (i) -=
					-tau_supg * JxW[qp] / dt * (velocity *
								    dphi[i]
								    [qp] *
								    w_old);
				    }
				}

			      //the newton term bro
			      if (!stokes)
				{
				  Fu (i) +=
				    tau_supg * JxW[qp] *
				    ((velocity * dphi[i][qp]) * (U * grad_u));
				  Fv (i) +=
				    tau_supg * JxW[qp] *
				    ((velocity * dphi[i][qp]) * (U * grad_v));
				  if (threed)
				    {
				      Fw (i) +=
					tau_supg * JxW[qp] *
					((velocity * dphi[i][qp]) *
					 (U * grad_w));
				    }
				}



			      for (unsigned int j = 0; j < n_u_dofs; j++)
				{

				  double laplacian_operator = 0;
				  laplacian_operator +=
				    d2phi[j][qp] (0, 0) + d2phi[j][qp] (1, 1);
				  if (threed)
				    {
				      laplacian_operator +=
					d2phi[j][qp] (2, 2);
				    }

				  if (!stokes)
				    {
				      Kuu (i, j) +=
					tau_supg * JxW[qp] * (phi[j][qp] *
							      grad_u (0) *
							      velocity *
							      dphi[i][qp] +
							      (U *
							       dphi[j][qp]) *
							      (velocity *
							       dphi[i][qp]));
				      Kuv (i, j) +=
					tau_supg * JxW[qp] * (phi[j][qp] *
							      grad_u (1) *
							      velocity *
							      dphi[i][qp]);
				      if (threed)
					{
					  Kuw (i, j) +=
					    tau_supg * JxW[qp] * (phi[j][qp] *
								  grad_u (2) *
								  velocity *
								  dphi[i]
								  [qp]);
					}
				      Kvu (i, j) +=
					tau_supg * JxW[qp] * (phi[j][qp] *
							      grad_v (0) *
							      velocity *
							      dphi[i][qp]);
				      Kvv (i, j) +=
					tau_supg * JxW[qp] * (phi[j][qp] *
							      grad_v (1) *
							      velocity *
							      dphi[i][qp] +
							      (U *
							       dphi[j][qp]) *
							      (velocity *
							       dphi[i][qp]));
				      if (threed)
					{
					  Kvw (i, j) +=
					    tau_supg * JxW[qp] * (phi[j][qp] *
								  grad_v (2) *
								  velocity *
								  dphi[i]
								  [qp]);
					}
				      if (threed)
					{
					  Kwu (i, j) +=
					    tau_supg * JxW[qp] * (phi[j][qp] *
								  grad_w (0) *
								  velocity *
								  dphi[i]
								  [qp]);
					}
				      if (threed)
					{
					  Kwv (i, j) +=
					    tau_supg * JxW[qp] * (phi[j][qp] *
								  grad_w (1) *
								  velocity *
								  dphi[i]
								  [qp]);
					}
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] * (phi[j][qp] *
								  grad_w (2) *
								  velocity *
								  dphi[i][qp]
								  +
								  (U *
								   dphi[j]
								   [qp]) *
								  (velocity *
								   dphi[i]
								   [qp]));
					}
				    }

				  if (supg_laplacian)
				    {
				      Kuu (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re *
							      velocity *
							      dphi[i][qp] *
							      laplacian_operator);
				      Kuv (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re *
							      velocity *
							      dphi[i][qp] *
							      laplacian_operator);
				      if (threed)
					{
					  Kuw (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  velocity *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      Kvu (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re *
							      velocity *
							      dphi[i][qp] *
							      laplacian_operator);
				      Kvv (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re *
							      velocity *
							      dphi[i][qp] *
							      laplacian_operator);
				      if (threed)
					{
					  Kvw (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  velocity *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kwu (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  velocity *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kwv (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  velocity *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  velocity *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}

				    }

				  if (unsteady)
				    {
				      Kuu (i, j) +=
					tau_supg * JxW[qp] * (1.0 / dt *
							      velocity *
							      dphi[i][qp] *
							      phi[j][qp]);
				      Kvv (i, j) +=
					tau_supg * JxW[qp] * (1.0 / dt *
							      velocity *
							      dphi[i][qp] *
							      phi[j][qp]);
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] * (1.0 / dt *
								  velocity *
								  dphi[i][qp]
								  *
								  phi[j][qp]);
					}

				    }
				}

			      for (unsigned int j = 0; j < n_p_dofs; j++)
				{
				  Kup (i, j) +=
				    tau_supg * JxW[qp] * (velocity *
							  dphi[i][qp] *
							  dpsi[j][qp] (0));
				  Kvp (i, j) +=
				    tau_supg * JxW[qp] * (velocity *
							  dphi[i][qp] *
							  dpsi[j][qp] (1));
				  if (threed)
				    {
				      Kwp (i, j) +=
					tau_supg * JxW[qp] * (velocity *
							      dphi[i][qp] *
							      dpsi[j][qp]
							      (2));
				    }
				}

			    }

			}
		      else if (es->parameters.get < bool >
			       ("supg_full_newton"))
			{


			  for (unsigned int i = 0; i < n_u_dofs; i++)
			    {

			      if (unsteady)
				{
				  Fu (i) -=
				    -tau_supg * JxW[qp] / dt * (U *
								dphi[i][qp] *
								u);
				  Fv (i) -=
				    -tau_supg * JxW[qp] / dt * (U *
								dphi[i][qp] *
								v);
				  if (threed)
				    {
				      Fw (i) -=
					-tau_supg * JxW[qp] / dt * (U *
								    dphi[i]
								    [qp] * w);
				    }
				}

			      //the newton term bro
			      if (!stokes)
				{
				  Fu (i) +=
				    2 * tau_supg * JxW[qp] *
				    ((U * dphi[i][qp]) * (U * grad_u));
				  Fv (i) +=
				    2 * tau_supg * JxW[qp] *
				    ((U * dphi[i][qp]) * (U * grad_v));
				  if (threed)
				    {
				      Fw (i) +=
					2 * tau_supg * JxW[qp] *
					((U * dphi[i][qp]) * (U * grad_w));
				    }
				}



			      for (unsigned int j = 0; j < n_u_dofs; j++)
				{

				  double laplacian_operator = 0;
				  laplacian_operator +=
				    d2phi[j][qp] (0, 0) + d2phi[j][qp] (1, 1);
				  if (threed)
				    {
				      laplacian_operator +=
					d2phi[j][qp] (2, 2);
				    }

				  // ignore the laplacian operator because we only really want to 
				  // don't need the pressure terms for some reason
				  if (!stokes)
				    {
				      Kuu (i, j) += tau_supg * JxW[qp] * (U * grad_u * phi[j][qp] * dphi[i][qp] (0) + (U * dphi[i][qp]) * (U * dphi[j][qp]) + U * dphi[i][qp] * phi[j][qp] * grad_u (0));	// + dphi[i][qp](0) * phi[j][qp]*grad_p(0));
				      Kuv (i, j) += tau_supg * JxW[qp] * (U * grad_u * phi[j][qp] * dphi[i][qp] (1) + U * dphi[i][qp] * phi[j][qp] * grad_u (1));	// + dphi[i][qp](1) * phi[j][qp]*grad_p(0));                
				      if (threed)
					{
					  Kuw (i, j) +=
					    tau_supg * JxW[qp] * (U * grad_u *
								  phi[j][qp] *
								  dphi[i][qp]
								  (2) +
								  U *
								  dphi[i][qp]
								  *
								  phi[j][qp] *
								  grad_u (2));
					}	// + dphi[i][qp](2) * phi[j][qp]*grad_p(0)); }          
				      Kvu (i, j) += tau_supg * JxW[qp] * (U * grad_v * phi[j][qp] * dphi[i][qp] (0) + U * dphi[i][qp] * phi[j][qp] * grad_v (0));	// + dphi[i][qp](0) * phi[j][qp]*grad_p(1));                
				      Kvv (i, j) += tau_supg * JxW[qp] * (U * grad_v * phi[j][qp] * dphi[i][qp] (1) + (U * dphi[i][qp]) * (U * dphi[j][qp]) + U * dphi[i][qp] * phi[j][qp] * grad_v (1));	// + dphi[i][qp](1) * phi[j][qp]*grad_p(1));                
				      if (threed)
					{
					  Kvw (i, j) +=
					    tau_supg * JxW[qp] * (U * grad_v *
								  phi[j][qp] *
								  dphi[i][qp]
								  (2) +
								  U *
								  dphi[i][qp]
								  *
								  phi[j][qp] *
								  grad_v (2));
					}	// + dphi[i][qp](2) * phi[j][qp]*grad_p(1)); }             
				      if (threed)
					{
					  Kwu (i, j) +=
					    tau_supg * JxW[qp] * (U * grad_w *
								  phi[j][qp] *
								  dphi[i][qp]
								  (0) +
								  U *
								  dphi[i][qp]
								  *
								  phi[j][qp] *
								  grad_w (0));
					}	// + dphi[i][qp](0) * phi[j][qp]*grad_p(2)); }             
				      if (threed)
					{
					  Kwv (i, j) +=
					    tau_supg * JxW[qp] * (U * grad_w *
								  phi[j][qp] *
								  dphi[i][qp]
								  (1) +
								  U *
								  dphi[i][qp]
								  *
								  phi[j][qp] *
								  grad_w (1));
					}	// + dphi[i][qp](1) * phi[j][qp]*grad_p(2)); }
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] * (U * grad_w *
								  phi[j][qp] *
								  dphi[i][qp]
								  (2) +
								  (U *
								   dphi[i]
								   [qp]) *
								  (U *
								   dphi[j]
								   [qp]) +
								  U *
								  dphi[i][qp]
								  *
								  phi[j][qp] *
								  grad_w (2));
					}	// + dphi[i][qp](2) * phi[j][qp]*grad_p(2)); }
				    }
				  else
				    {
				      /*
				         Kuu(i,j) +=                                                  tau_supg*JxW[qp]*(dphi[i][qp](0) * phi[j][qp]*grad_p(0));
				         Kuv(i,j) +=                                                  tau_supg*JxW[qp]*(dphi[i][qp](1) * phi[j][qp]*grad_p(0));                
				         if(threed) { Kuw(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](2) * phi[j][qp]*grad_p(0)); }          
				         Kvu(i,j) +=                                                  tau_supg*JxW[qp]*(dphi[i][qp](0) * phi[j][qp]*grad_p(1));                
				         Kvv(i,j) +=                                                  tau_supg*JxW[qp]*(dphi[i][qp](1) * phi[j][qp]*grad_p(1));                
				         if(threed) { Kvw(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](2) * phi[j][qp]*grad_p(1)); }             
				         if(threed) { Kwu(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](0) * phi[j][qp]*grad_p(2)); }             
				         if(threed) { Kwv(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](1) * phi[j][qp]*grad_p(2)); }
				         if(threed) { Kww(i,j) += tau_supg*JxW[qp]*(dphi[i][qp](2) * phi[j][qp]*grad_p(2)); }
				       */
				    }

				  if (supg_laplacian)
				    {
				      Kuu (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re * U *
							      dphi[i][qp] *
							      laplacian_operator);
				      Kuv (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re * U *
							      dphi[i][qp] *
							      laplacian_operator);
				      if (threed)
					{
					  Kuw (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  U *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      Kvu (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re * U *
							      dphi[i][qp] *
							      laplacian_operator);
				      Kvv (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re * U *
							      dphi[i][qp] *
							      laplacian_operator);
				      if (threed)
					{
					  Kvw (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  U *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kwu (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  U *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kwv (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  U *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  U *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				    }


				  if (unsteady)
				    {
				      Kuu (i, j) +=
					tau_supg * JxW[qp] / dt * (U *
								   dphi[i][qp]
								   *
								   phi[j]
								   [qp]);
				      Kvv (i, j) +=
					tau_supg * JxW[qp] / dt * (U *
								   dphi[i][qp]
								   *
								   phi[j]
								   [qp]);
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] / dt * (U *
								       dphi[i]
								       [qp] *
								       phi[j]
								       [qp]);
					}

				      Kuu (i, j) +=
					tau_supg * JxW[qp] / dt * (u *
								   dphi[i][qp]
								   (0) *
								   phi[j][qp]
								   -
								   u_old *
								   dphi[i][qp]
								   (0) *
								   phi[j]
								   [qp]);
				      Kuv (i, j) +=
					tau_supg * JxW[qp] / dt * (u *
								   dphi[i][qp]
								   (1) *
								   phi[j][qp]
								   -
								   u_old *
								   dphi[i][qp]
								   (1) *
								   phi[j]
								   [qp]);
				      if (threed)
					{
					  Kuw (i, j) +=
					    tau_supg * JxW[qp] / dt * (u *
								       dphi[i]
								       [qp]
								       (2) *
								       phi[j]
								       [qp] -
								       u_old *
								       dphi[i]
								       [qp]
								       (2) *
								       phi[j]
								       [qp]);
					}
				      Kvu (i, j) +=
					tau_supg * JxW[qp] / dt * (v *
								   dphi[i][qp]
								   (0) *
								   phi[j][qp]
								   -
								   v_old *
								   dphi[i][qp]
								   (0) *
								   phi[j]
								   [qp]);
				      Kvv (i, j) +=
					tau_supg * JxW[qp] / dt * (v *
								   dphi[i][qp]
								   (1) *
								   phi[j][qp]
								   -
								   v_old *
								   dphi[i][qp]
								   (1) *
								   phi[j]
								   [qp]);
				      if (threed)
					{
					  Kvw (i, j) +=
					    tau_supg * JxW[qp] / dt * (v *
								       dphi[i]
								       [qp]
								       (2) *
								       phi[j]
								       [qp] -
								       v_old *
								       dphi[i]
								       [qp]
								       (2) *
								       phi[j]
								       [qp]);
					}
				      if (threed)
					{
					  Kwu (i, j) +=
					    tau_supg * JxW[qp] / dt * (w *
								       dphi[i]
								       [qp]
								       (0) *
								       phi[j]
								       [qp] -
								       w_old *
								       dphi[i]
								       [qp]
								       (0) *
								       phi[j]
								       [qp]);
					}
				      if (threed)
					{
					  Kwv (i, j) +=
					    tau_supg * JxW[qp] / dt * (w *
								       dphi[i]
								       [qp]
								       (1) *
								       phi[j]
								       [qp] -
								       w_old *
								       dphi[i]
								       [qp]
								       (1) *
								       phi[j]
								       [qp]);
					}
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] / dt * (w *
								       dphi[i]
								       [qp]
								       (2) *
								       phi[j]
								       [qp] -
								       w_old *
								       dphi[i]
								       [qp]
								       (2) *
								       phi[j]
								       [qp]);
					}

				    }
				}

			      for (unsigned int j = 0; j < n_p_dofs; j++)
				{
				  Kup (i, j) +=
				    tau_supg * JxW[qp] * (U * dphi[i][qp] *
							  dpsi[j][qp] (0));
				  Kvp (i, j) +=
				    tau_supg * JxW[qp] * (U * dphi[i][qp] *
							  dpsi[j][qp] (1));
				  if (threed)
				    {
				      Kwp (i, j) +=
					tau_supg * JxW[qp] * (U *
							      dphi[i][qp] *
							      dpsi[j][qp]
							      (2));
				    }
				}
			    }

			}
		      else if (es->parameters.get < bool >
			       ("supg_convection_newton"))
			{


			  for (unsigned int i = 0; i < n_u_dofs; i++)
			    {

			      //the newton term bro
			      if (!stokes)
				{
				  Fu (i) +=
				    2. * tau_supg * JxW[qp] *
				    ((U * dphi[i][qp]) * (U * grad_u));
				  Fv (i) +=
				    2. * tau_supg * JxW[qp] *
				    ((U * dphi[i][qp]) * (U * grad_v));
				  if (threed)
				    {
				      Fw (i) +=
					2 * tau_supg * JxW[qp] *
					((U * dphi[i][qp]) * (U * grad_w));
				    }
				}



			      for (unsigned int j = 0; j < n_u_dofs; j++)
				{
				  // ignore the laplacian operator because we only really want to 

				  double laplacian_operator = 0;
				  laplacian_operator +=
				    d2phi[j][qp] (0, 0) + d2phi[j][qp] (1, 1);
				  if (threed)
				    {
				      laplacian_operator +=
					d2phi[j][qp] (2, 2);
				    }

				  if (!stokes)
				    {
				      Kuu (i, j) +=
					tau_supg * JxW[qp] * (U * grad_u *
							      phi[j][qp] *
							      dphi[i][qp] (0)
							      +
							      (U *
							       dphi[i][qp]) *
							      (U *
							       dphi[j][qp]) +
							      U *
							      dphi[i][qp] *
							      phi[j][qp] *
							      grad_u (0));
				      Kuv (i, j) +=
					tau_supg * JxW[qp] * (U * grad_u *
							      phi[j][qp] *
							      dphi[i][qp] (1)
							      +
							      U *
							      dphi[i][qp] *
							      phi[j][qp] *
							      grad_u (1));
				      if (threed)
					{
					  Kuw (i, j) +=
					    tau_supg * JxW[qp] * (U * grad_u *
								  phi[j][qp] *
								  dphi[i][qp]
								  (2) +
								  U *
								  dphi[i][qp]
								  *
								  phi[j][qp] *
								  grad_u (2));
					}
				      Kvu (i, j) +=
					tau_supg * JxW[qp] * (U * grad_v *
							      phi[j][qp] *
							      dphi[i][qp] (0)
							      +
							      U *
							      dphi[i][qp] *
							      phi[j][qp] *
							      grad_v (0));
				      Kvv (i, j) +=
					tau_supg * JxW[qp] * (U * grad_v *
							      phi[j][qp] *
							      dphi[i][qp] (1)
							      +
							      (U *
							       dphi[i][qp]) *
							      (U *
							       dphi[j][qp]) +
							      U *
							      dphi[i][qp] *
							      phi[j][qp] *
							      grad_v (1));
				      if (threed)
					{
					  Kvw (i, j) +=
					    tau_supg * JxW[qp] * (U * grad_v *
								  phi[j][qp] *
								  dphi[i][qp]
								  (2) +
								  U *
								  dphi[i][qp]
								  *
								  phi[j][qp] *
								  grad_v (2));
					}
				      if (threed)
					{
					  Kwu (i, j) +=
					    tau_supg * JxW[qp] * (U * grad_w *
								  phi[j][qp] *
								  dphi[i][qp]
								  (0) +
								  U *
								  dphi[i][qp]
								  *
								  phi[j][qp] *
								  grad_w (0));
					}
				      if (threed)
					{
					  Kwv (i, j) +=
					    tau_supg * JxW[qp] * (U * grad_w *
								  phi[j][qp] *
								  dphi[i][qp]
								  (1) +
								  U *
								  dphi[i][qp]
								  *
								  phi[j][qp] *
								  grad_w (1));
					}
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] * (U * grad_w *
								  phi[j][qp] *
								  dphi[i][qp]
								  (2) +
								  (U *
								   dphi[i]
								   [qp]) *
								  (U *
								   dphi[j]
								   [qp]) +
								  U *
								  dphi[i][qp]
								  *
								  phi[j][qp] *
								  grad_w (2));
					}
				    }
				  if (supg_laplacian)
				    {
				      Kuu (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re * U *
							      dphi[i][qp] *
							      laplacian_operator);
				      Kuv (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re * U *
							      dphi[i][qp] *
							      laplacian_operator);
				      if (threed)
					{
					  Kuw (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  U *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      Kvu (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re * U *
							      dphi[i][qp] *
							      laplacian_operator);
				      Kvv (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re * U *
							      dphi[i][qp] *
							      laplacian_operator);
				      if (threed)
					{
					  Kvw (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  U *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kwu (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  U *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kwv (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  U *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  U *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				    }

				}
			    }

			}
		      else if (es->parameters.get < bool >
			       ("supg_convection"))
			{

			  if (!stokes)
			    {
			      for (unsigned int i = 0; i < n_u_dofs; i++)
				{


				  //the newton term bro
				  Fu (i) +=
				    tau_supg * JxW[qp] *
				    ((velocity * dphi[i][qp]) * (U * grad_u));
				  Fv (i) +=
				    tau_supg * JxW[qp] *
				    ((velocity * dphi[i][qp]) * (U * grad_v));
				  if (threed)
				    {
				      Fw (i) +=
					tau_supg * JxW[qp] *
					((velocity * dphi[i][qp]) *
					 (U * grad_w));
				    }



				  for (unsigned int j = 0; j < n_u_dofs; j++)
				    {

				      Kuu (i, j) +=
					tau_supg * JxW[qp] * (phi[j][qp] *
							      grad_u (0) *
							      velocity *
							      dphi[i][qp] +
							      (U *
							       dphi[j][qp]) *
							      (velocity *
							       dphi[i][qp]));
				      Kuv (i, j) +=
					tau_supg * JxW[qp] * (phi[j][qp] *
							      grad_u (1) *
							      velocity *
							      dphi[i][qp]);
				      if (threed)
					{
					  Kuw (i, j) +=
					    tau_supg * JxW[qp] * (phi[j][qp] *
								  grad_u (2) *
								  velocity *
								  dphi[i]
								  [qp]);
					}
				      Kvu (i, j) +=
					tau_supg * JxW[qp] * (phi[j][qp] *
							      grad_v (0) *
							      velocity *
							      dphi[i][qp]);
				      Kvv (i, j) +=
					tau_supg * JxW[qp] * (phi[j][qp] *
							      grad_v (1) *
							      velocity *
							      dphi[i][qp] +
							      (U *
							       dphi[j][qp]) *
							      (velocity *
							       dphi[i][qp]));
				      if (threed)
					{
					  Kvw (i, j) +=
					    tau_supg * JxW[qp] * (phi[j][qp] *
								  grad_v (2) *
								  velocity *
								  dphi[i]
								  [qp]);
					}
				      if (threed)
					{
					  Kwu (i, j) +=
					    tau_supg * JxW[qp] * (phi[j][qp] *
								  grad_w (0) *
								  velocity *
								  dphi[i]
								  [qp]);
					}
				      if (threed)
					{
					  Kwv (i, j) +=
					    tau_supg * JxW[qp] * (phi[j][qp] *
								  grad_w (1) *
								  velocity *
								  dphi[i]
								  [qp]);
					}
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] * (phi[j][qp] *
								  grad_w (2) *
								  velocity *
								  dphi[i][qp]
								  +
								  (U *
								   dphi[j]
								   [qp]) *
								  (velocity *
								   dphi[i]
								   [qp]));
					}
				    }

				}
			    }

			}
		      else if (es->parameters.get < bool > ("supg_picard"))
			{


			  for (unsigned int i = 0; i < n_u_dofs; i++)
			    {

			      if (unsteady)
				{
				  Fu (i) -=
				    -tau_supg * JxW[qp] / dt * (velocity *
								dphi[i][qp] *
								u_old);
				  Fv (i) -=
				    -tau_supg * JxW[qp] / dt * (velocity *
								dphi[i][qp] *
								v_old);
				  if (threed)
				    {
				      Fw (i) -=
					-tau_supg * JxW[qp] / dt * (velocity *
								    dphi[i]
								    [qp] *
								    w_old);
				    }
				}

			      for (unsigned int j = 0; j < n_u_dofs; j++)
				{

				  double laplacian_operator = 0;
				  laplacian_operator +=
				    d2phi[j][qp] (0, 0) + d2phi[j][qp] (1, 1);
				  if (threed)
				    {
				      laplacian_operator +=
					d2phi[j][qp] (2, 2);
				    }

				  if (!stokes)
				    {
				      Kuu (i, j) +=
					tau_supg * JxW[qp] *
					((U * dphi[j][qp]) *
					 (velocity * dphi[i][qp]));
				      Kvv (i, j) +=
					tau_supg * JxW[qp] *
					((U * dphi[j][qp]) *
					 (velocity * dphi[i][qp]));
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] *
					    ((U * dphi[j][qp]) *
					     (velocity * dphi[i][qp]));
					}
				    }
				  if (supg_laplacian)
				    {
				      Kuu (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re *
							      velocity *
							      dphi[i][qp] *
							      laplacian_operator);
				      Kuv (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re *
							      velocity *
							      dphi[i][qp] *
							      laplacian_operator);
				      if (threed)
					{
					  Kuw (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  velocity *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      Kvu (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re *
							      velocity *
							      dphi[i][qp] *
							      laplacian_operator);
				      Kvv (i, j) +=
					tau_supg * JxW[qp] * (-1.0 / Re *
							      velocity *
							      dphi[i][qp] *
							      laplacian_operator);
				      if (threed)
					{
					  Kvw (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  velocity *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kwu (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  velocity *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kwv (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  velocity *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] * (-1.0 / Re *
								  velocity *
								  dphi[i][qp]
								  *
								  laplacian_operator);
					}
				    }

				  if (unsteady)
				    {
				      Kuu (i, j) +=
					tau_supg * JxW[qp] * (1.0 / dt *
							      velocity *
							      dphi[i][qp] *
							      phi[j][qp]);
				      Kvv (i, j) +=
					tau_supg * JxW[qp] * (1.0 / dt *
							      velocity *
							      dphi[i][qp] *
							      phi[j][qp]);
				      if (threed)
					{
					  Kww (i, j) +=
					    tau_supg * JxW[qp] * (1.0 / dt *
								  velocity *
								  dphi[i][qp]
								  *
								  phi[j][qp]);
					}

				    }
				}

			      for (unsigned int j = 0; j < n_p_dofs; j++)
				{
				  Kup (i, j) +=
				    tau_supg * JxW[qp] * (velocity *
							  dphi[i][qp] *
							  dpsi[j][qp] (0));
				  Kvp (i, j) +=
				    tau_supg * JxW[qp] * (velocity *
							  dphi[i][qp] *
							  dpsi[j][qp] (1));
				  if (threed)
				    {
				      Kwp (i, j) +=
					tau_supg * JxW[qp] * (velocity *
							      dphi[i][qp] *
							      dpsi[j][qp]
							      (2));
				    }
				}

			    }

			}
		      //else if(es->parameters.get<bool>("supg_newton"))
		      //{
		      //}

		    }

		  if (lsic)
		    {
		      for (unsigned int i = 0; i < n_u_dofs; i++)
			{
			  for (unsigned int j = 0; j < n_u_dofs; j++)
			    {

			      Kuu (i, j) +=
				tau_lsic * JxW[qp] * (dphi[i][qp] (0) *
						      dphi[j][qp] (0));
			      Kuv (i, j) +=
				tau_lsic * JxW[qp] * (dphi[i][qp] (0) *
						      dphi[j][qp] (1));
			      if (threed)
				{
				  Kuw (i, j) +=
				    tau_lsic * JxW[qp] * (dphi[i][qp] (0) *
							  dphi[j][qp] (2));
				}
			      Kvu (i, j) +=
				tau_lsic * JxW[qp] * (dphi[i][qp] (1) *
						      dphi[j][qp] (0));
			      Kvv (i, j) +=
				tau_lsic * JxW[qp] * (dphi[i][qp] (1) *
						      dphi[j][qp] (1));
			      if (threed)
				{
				  Kvw (i, j) +=
				    tau_lsic * JxW[qp] * (dphi[i][qp] (1) *
							  dphi[j][qp] (2));
				}
			      if (threed)
				{
				  Kwu (i, j) +=
				    tau_lsic * JxW[qp] * (dphi[i][qp] (2) *
							  dphi[j][qp] (0));
				}
			      if (threed)
				{
				  Kwv (i, j) +=
				    tau_lsic * JxW[qp] * (dphi[i][qp] (2) *
							  dphi[j][qp] (1));
				}
			      if (threed)
				{
				  Kww (i, j) +=
				    tau_lsic * JxW[qp] * (dphi[i][qp] (2) *
							  dphi[j][qp] (2));
				}
			    }
			}
		    }

		  // note the negative dt that has been introduced to preserve symmetry
		  if (pspg)
		    {

		      if (es->parameters.get < bool > ("pspg_newton"))
			{

			  NumberVectorValue grad_Ux;
			  NumberVectorValue grad_Uy;
			  NumberVectorValue grad_Uz;
			  if (threed)
			    {
			      grad_Ux =
				NumberVectorValue (grad_u (0), grad_v (0),
						   grad_w (0));
			      grad_Uy =
				NumberVectorValue (grad_u (1), grad_v (1),
						   grad_w (1));
			      grad_Uz =
				NumberVectorValue (grad_u (2), grad_v (2),
						   grad_w (2));
			    }
			  else
			    {
			      grad_Ux =
				NumberVectorValue (grad_u (0), grad_v (0));
			      grad_Uy =
				NumberVectorValue (grad_u (1), grad_v (1));
			    }

			  for (unsigned int i = 0; i < n_p_dofs; i++)
			    {

			      if (unsteady)
				{
				  Fp (i) +=
				    -tau_pspg * JxW[qp] / dt * (dpsi[i][qp] *
								U_old);
				}

			      //the newton term bro
			      if (!stokes)
				{
				  Fp (i) +=
				    -tau_pspg * JxW[qp] * (dpsi[i][qp] (0) *
							   (U * grad_u) +
							   dpsi[i][qp] (1) *
							   (U * grad_v));
				  if (threed)
				    {
				      Fp (i) +=
					-tau_pspg * JxW[qp] *
					(dpsi[i][qp] (2) * (U * grad_w));
				    }
				}



			      for (unsigned int j = 0; j < n_u_dofs; j++)
				{

				  double laplacian_operator = 0;
				  laplacian_operator +=
				    d2phi[j][qp] (0, 0) + d2phi[j][qp] (1, 1);
				  if (threed)
				    {
				      laplacian_operator +=
					d2phi[j][qp] (2, 2);
				    }

				  if (!stokes)
				    {
				      Kpu (i, j) +=
					-tau_pspg * JxW[qp] * (phi[j][qp] *
							       grad_Ux *
							       dpsi[i][qp] +
							       (U *
								dphi[j][qp]) *
							       dpsi[i][qp]
							       (0));
				      Kpv (i, j) +=
					-tau_pspg * JxW[qp] * (phi[j][qp] *
							       grad_Uy *
							       dpsi[i][qp] +
							       (U *
								dphi[j][qp]) *
							       dpsi[i][qp]
							       (1));
				      if (threed)
					{
					  Kpw (i, j) +=
					    -tau_pspg * JxW[qp] *
					    (phi[j][qp] * grad_Uz *
					     dpsi[i][qp] +
					     (U * dphi[j][qp]) *
					     dpsi[i][qp] (2));
					}
				    }

				  if (supg_laplacian)
				    {
				      Kpu (i, j) +=
					-tau_pspg * JxW[qp] * (-1.0 / Re *
							       dpsi[i][qp] (0)
							       *
							       laplacian_operator);
				      Kpv (i, j) +=
					-tau_pspg * JxW[qp] * (-1.0 / Re *
							       dpsi[i][qp] (1)
							       *
							       laplacian_operator);
				      if (threed)
					{
					  Kpw (i, j) +=
					    -tau_pspg * JxW[qp] * (-1.0 / Re *
								   dpsi[i][qp]
								   (2) *
								   laplacian_operator);
					}
				    }

				  if (unsteady)
				    {
				      Kpu (i, j) +=
					-tau_pspg * JxW[qp] * (1.0 / dt *
							       dpsi[i][qp] (0)
							       * phi[j][qp]);
				      Kpv (i, j) +=
					-tau_pspg * JxW[qp] * (1.0 / dt *
							       dpsi[i][qp] (1)
							       * phi[j][qp]);
				      if (threed)
					{
					  Kpw (i, j) +=
					    -tau_pspg * JxW[qp] * (1.0 / dt *
								   dpsi[i][qp]
								   (2) *
								   phi[j]
								   [qp]);
					}

				    }
				}

			      // it is this term that is screwing things up grrr
			      for (unsigned int j = 0; j < n_p_dofs; j++)
				{
				  Kpp (i, j) +=
				    -tau_pspg * JxW[qp] * (dpsi[i][qp] *
							   dpsi[j][qp]);
				}

			    }

			}
		      else if (es->parameters.get < bool > ("pspg_picard"))
			{


			  for (unsigned int i = 0; i < n_p_dofs; i++)
			    {

			      if (unsteady)
				{
				  Fp (i) +=
				    -tau_pspg * JxW[qp] / dt * (dpsi[i][qp] *
								U_old);
				}

			      for (unsigned int j = 0; j < n_u_dofs; j++)
				{

				  double laplacian_operator = 0;
				  laplacian_operator +=
				    d2phi[j][qp] (0, 0) + d2phi[j][qp] (1, 1);
				  if (threed)
				    {
				      laplacian_operator +=
					d2phi[j][qp] (2, 2);
				    }

				  if (!stokes)
				    {
				      Kpu (i, j) +=
					-tau_pspg * JxW[qp] *
					((U * dphi[j][qp]) * dpsi[i][qp] (0));
				      Kpv (i, j) +=
					-tau_pspg * JxW[qp] *
					((U * dphi[j][qp]) * dpsi[i][qp] (1));
				      if (threed)
					{
					  Kpw (i, j) +=
					    -tau_pspg * JxW[qp] *
					    ((U * dphi[j][qp]) *
					     dpsi[i][qp] (2));
					}
				    }
				  if (supg_laplacian)
				    {
				      Kpu (i, j) +=
					-tau_pspg * JxW[qp] * (-1.0 / Re *
							       dpsi[i][qp] (0)
							       *
							       laplacian_operator);
				      Kpv (i, j) +=
					-tau_pspg * JxW[qp] * (-1.0 / Re *
							       dpsi[i][qp] (1)
							       *
							       laplacian_operator);
				      if (threed)
					{
					  Kpw (i, j) +=
					    -tau_pspg * JxW[qp] * (-1.0 / Re *
								   dpsi[i][qp]
								   (2) *
								   laplacian_operator);
					}
				    }


				  if (unsteady)
				    {
				      Kpu (i, j) +=
					-tau_pspg * JxW[qp] * (1.0 / dt *
							       dpsi[i][qp] (0)
							       * phi[j][qp]);
				      Kpv (i, j) +=
					-tau_pspg * JxW[qp] * (1.0 / dt *
							       dpsi[i][qp] (1)
							       * phi[j][qp]);
				      if (threed)
					{
					  Kpw (i, j) +=
					    -tau_pspg * JxW[qp] * (1.0 / dt *
								   dpsi[i][qp]
								   (2) *
								   phi[j]
								   [qp]);
					}

				    }
				}

			      for (unsigned int j = 0; j < n_p_dofs; j++)
				{
				  Kpp (i, j) +=
				    -tau_pspg * JxW[qp] * (dpsi[i][qp] *
							   dpsi[j][qp]);
				}


			    }

			}
		    }

		}
	      // ************** END SUPG TERMS ************************************//
	      // *************************************************************************** //

	    }			// end of the quadrature point qp-loop













	  // **************** BOUNDARY TERMS/CONDITIONS **************//
	  {
	    // The penalty value.  \f$ \frac{1}{\epsilon} \f$
	    const Real penalty = 1.e10;

	    for (unsigned int s = 0; s < elem->n_sides (); s++)
	      {
		//for some reason it is natural to have more than one boundary id per side or even node
		std::vector < boundary_id_type > boundary_ids =
		  mesh.boundary_info->boundary_ids (elem, s);

		// if there is a boundary id on this side
		if (boundary_ids.size () > 0)
		  {
		    int boundary_id = boundary_ids[0];	// should only have one hopefully






		    // *********** PRECONDITIONER BOUNDARY CONDITIONS ************** //
		    // want to get the pressure dofs that are on the dirichlet inflow boundary
		    // we assume that the dirichlet boundary inflow has boundary id 0
		    // but dof_indices_p includes the dofs that are in the whole element.. hmmm
		    // use FEInterface::dofs_on_side

		    if ((es->parameters.get <
			 unsigned int >("preconditioner_type_3d") == 4
			 || es->parameters.get <
			 unsigned int >("preconditioner_type_3d") ==
			 5) &&es->parameters.get <
			unsigned int >("problem_type") != 4)
		      {
			if (boundary_id == 0)
			  {

			    if (es->parameters.get <
				unsigned int >("pcd_boundary_condition_type")
				== 1
				|| es->parameters.get <
				unsigned int >("pcd_boundary_condition_type")
				== 3
				|| es->parameters.get <
				unsigned int >("pcd_boundary_condition_type")
				== 4)
			      {
				std::vector <
				  unsigned int >pressure_dofs_on_side;
				FEInterface::dofs_on_side (elem, dim,
							   fe_pres_type, s,
							   pressure_dofs_on_side);
				// presumably this returns their place in dof_indices_p
				// will get some overlap for cts elements but this doesn't matter for later on
				for (unsigned int i = 0;
				     i < pressure_dofs_on_side.size (); i++)
				  {
				    pressure_dofs_on_inflow_boundary.
				      push_back (dof_indices_p
						 [pressure_dofs_on_side[i]]);
				    if (es->parameters.get <
					unsigned int
					>("pcd_boundary_condition_type") == 3)
				      all_global_pressure_dofs_on_inflow_boundary.
					push_back (dof_indices_p
						   [pressure_dofs_on_side
						    [i]]);
				  }
			      }
			    else if (es->parameters.get <
				     unsigned int
				     >("pcd_boundary_condition_type") == 2)
			      {


				fe_vel_face->reinit (elem, s);
				fe_pres_face->reinit (elem, s);


				std::vector < Real > JxW_face = JxW_face_elem;

				if (multiply_system_by_dt)
				  {
				    for (unsigned int qp = 0;
					 qp < qface.n_points (); qp++)
				      {
					JxW_face[qp] *= dt;
				      }
				  }

				for (unsigned int qp = 0;
				     qp < qface.n_points (); qp++)
				  {
				    // calculate the velocity at gauss points on the face
				    Number u = 0., v = 0., w = 0.;

				    for (unsigned int l = 0; l < n_u_dofs;
					 l++)
				      {
					// From the previous Newton iterate:
					u +=
					  phi_face[l][qp] *
					  system->
					  current_solution (dof_indices_u[l]);
					v +=
					  phi_face[l][qp] *
					  system->
					  current_solution (dof_indices_v[l]);
					if (threed)
					  w +=
					    phi_face[l][qp] *
					    system->
					    current_solution (dof_indices_w
							      [l]);
				      }

				    NumberVectorValue U;
				    if (threed)
				      U = NumberVectorValue (u, v, w);
				    else
				      U = NumberVectorValue (u, v);


				    // add robin terms to the matrix
				    for (unsigned int i = 0; i < n_p_dofs;
					 i++)
				      {
					for (unsigned int j = 0; j < n_p_dofs;
					     j++)
					  {
					    Kpp_pre_convection_diffusion (i,
									  j)
					      -=
					      JxW_face[qp] * (U *
							      qface_normals
							      [qp]) *
					      psi_face[i][qp] *
					      psi_face[j][qp];
					  }
				      }
				  }
			      }
			  }
		      }






		    // *********** STRESS BOUNDARY CONDITIONS/TERMS **************** //
		    // stress boundary conditions/terms, which are only applied if not coupled
		    // otherwise we put these contributions in the off diagonal block.
		    // will have to think about this wrt the stabilisation term on the boundary

		    if (!pressure_coupled)
		      {
			// mean pressure conditions/terms
			if (bc_type[boundary_id].compare ("pressure") == 0
			    || bc_type[boundary_id].compare ("stress") == 0
			    || bc_type[boundary_id].compare ("neumann") == 0)
			  {

			    // mean pressure or moghadam resistance
			    double mean_pressure = bc_value[boundary_id];
			    //mean_pressure = 2.;

			    // TEST
			    //std::cout << "resistance = " << mean_pressure << std::endl;

			    if (bc_type[boundary_id].compare ("neumann") == 0)
			      mean_pressure = 0;

			    // if coupled then we want the inflow to possibly be timeflow controlled, but future work
			    if (es->parameters.get <
				unsigned int >("sim_type") == 0
				|| es->parameters.get <
				unsigned int >("sim_type") == 2)
			      mean_pressure *=
				es->parameters.get < double >("time_scaling");

			    //at the moment we just to a vector in normal direction for stress
			    double stress_mag = bc_value[boundary_id];
			    if (es->parameters.get <
				unsigned int >("sim_type") == 0
				|| es->parameters.get <
				unsigned int >("sim_type") == 2)
			      stress_mag *=
				es->parameters.get < double >("time_scaling");

			    DenseMatrix < Number > stress;
			    DenseVector < Number > normal_stress;
			    unsigned int dimension = 3;
			    if (!threed)
			      dimension = 2;
			    stress.resize (dimension, dimension);
			    normal_stress.resize (dimension);

			    //make stress identity times constant
			    stress (0, 0) = 1.0;
			    stress (1, 1) = 1.0;
			    if (threed)
			      stress (2, 2) = 1.0;
			    stress *= stress_mag;



			    fe_vel_face->reinit (elem, s);
			    fe_pres_face->reinit (elem, s);


			    std::vector < Real > JxW_face = JxW_face_elem;
			    //if axisym then need to multiply all integrals by "r" i.e. y = qpoint[qp](1)
			    if (es->parameters.get <
				unsigned int >("geometry_type") == 5)
			      {
				for (unsigned int qp = 0;
				     qp < qface.n_points (); qp++)
				  JxW_face[qp] *= qpoint_face[qp] (1);
			      }


			    if (multiply_system_by_dt)
			      {
				for (unsigned int qp = 0;
				     qp < qface.n_points (); qp++)
				  {
				    JxW_face[qp] *= dt;
				  }
			      }

			    std::vector < double >temp_vector (n_u_dofs);
			    std::vector < std::vector <
			      double > >temp_matrix_1 (n_u_dofs,
						      std::vector <
						      double >(n_u_dofs));
			    std::vector < std::vector <
			      double > >temp_matrix_2 (n_u_dofs,
						      std::vector <
						      double >(n_u_dofs));
			    double extra_terms = 0;
			    double actual_terms = 0;
			    double actual_terms_2 = 0.;
			    //std::cout << "no gps = " << qface.n_points() << std::endl;
			    //std::cout << "no u_dofs = " << n_u_dofs << std::endl;

			    for (unsigned int qp = 0; qp < qface.n_points ();
				 qp++)
			      {

				//std::cout << "face qpoint = " << qpoint_face[qp] << std::endl;

				if (bc_type[boundary_id].compare ("stress") ==
				    0)
				  {
				    normal_stress (0) =
				      stress (0,
					      0) * qface_normals[qp] (0) +
				      stress (0, 1) * qface_normals[qp] (1);
				    if (threed)
				      {
					normal_stress (0) +=
					  stress (0,
						  2) * qface_normals[qp] (2);
				      }

				    normal_stress (1) =
				      stress (1,
					      0) * qface_normals[qp] (0) +
				      stress (1, 1) * qface_normals[qp] (1);
				    if (threed)
				      {
					normal_stress (1) +=
					  stress (1,
						  2) * qface_normals[qp] (2);
				      }

				    if (threed)
				      normal_stress (2) =
					stress (2,
						0) * qface_normals[qp] (0) +
					stress (2,
						1) * qface_normals[qp] (1) +
					stress (2, 2) * qface_normals[qp] (2);
				  }
				else if (bc_type[boundary_id].
					 compare ("pressure") == 0)
				  {
				    normal_stress (0) =
				      mean_pressure * qface_normals[qp] (0);
				    normal_stress (1) =
				      mean_pressure * qface_normals[qp] (1);
				    if (threed)
				      {
					normal_stress (2) =
					  mean_pressure *
					  qface_normals[qp] (2);
				      }
				  }


				if (!es->parameters.get < bool >
				    ("moghadam_coupling"))
				  {
				    for (unsigned int i = 0; i < n_u_dofs;
					 i++)
				      {

					Fu (i) +=
					  -JxW_face[qp] *
					  (((normal_stress (0)) *
					    phi_face[i][qp]));
					Fv (i) +=
					  -JxW_face[qp] *
					  (((normal_stress (1)) *
					    phi_face[i][qp]));
					if (threed)
					  {
					    Fw (i) +=
					      -JxW_face[qp] *
					      (((normal_stress (2)) *
						phi_face[i][qp]));
					  }
				      }
				  }
				else
				  {
				    // for moghadam we need a sub loop over quad points to do a 

				    for (unsigned int i = 0; i < n_u_dofs;
					 i++)
				      {

					temp_vector[i] +=
					  JxW_face[qp] *
					  (qface_normals[qp] (0) *
					   phi_face[i][qp]);


					for (unsigned int qp2 = 0;
					     qp2 < qface.n_points (); qp2++)
					  {

					    for (unsigned int j = 0;
						 j < n_u_dofs; j++)
					      {
						temp_matrix_1[i][j] +=
						  JxW_face[qp] *
						  JxW_face[qp2] *
						  mean_pressure *
						  (qface_normals[qp] (0) *
						   phi_face[i][qp] *
						   qface_normals[qp2] (0) *
						   phi_face[j][qp2]);
						actual_terms +=
						  JxW_face[qp] *
						  JxW_face[qp2] *
						  mean_pressure *
						  (qface_normals[qp] (0) *
						   phi_face[i][qp] *
						   qface_normals[qp2] (0) *
						   phi_face[j][qp2]);

						/*
						   Kuu(i,j) += JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](0)*phi_face[i][qp] * qface_normals[qp2](0)*phi_face[j][qp2]);
						   Kuv(i,j) += JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](1)*phi_face[i][qp] * qface_normals[qp2](0)*phi_face[j][qp2]);
						   if(threed) { Kuw(i,j) += JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](2)*phi_face[i][qp] * qface_normals[qp2](0)*phi_face[j][qp2]); }
						   Kvu(i,j) += JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](0)*phi_face[i][qp] * qface_normals[qp2](1)*phi_face[j][qp2]);
						   Kvv(i,j) += JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](1)*phi_face[i][qp] * qface_normals[qp2](1)*phi_face[j][qp2]);
						   if(threed) { Kvw(i,j) += JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](2)*phi_face[i][qp] * qface_normals[qp2](1)*phi_face[j][qp2]); }
						   if(threed)
						   {                                                                                            
						   Kwu(i,j) += JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](0)*phi_face[i][qp] * qface_normals[qp2](2)*phi_face[j][qp2]);
						   Kwv(i,j) += JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](1)*phi_face[i][qp] * qface_normals[qp2](2)*phi_face[j][qp2]);
						   Kww(i,j) += JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](2)*phi_face[i][qp] * qface_normals[qp2](2)*phi_face[j][qp2]);
						   }
						   extra_terms +=JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](1)*phi_face[i][qp] * qface_normals[qp2](0)*phi_face[j][qp2]) + JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](0)*phi_face[i][qp] * qface_normals[qp2](1)*phi_face[j][qp2]) + JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](1)*phi_face[i][qp] * qface_normals[qp2](1)*phi_face[j][qp2]);
						   actual_terms += mean_pressure*(qface_normals[qp](0)*phi_face[i][qp] * qface_normals[qp2](0)*phi_face[j][qp2]);
						 */
						//std::cout << "extra terms = " <<  << std::endl;
						//std::cout << "actual term = " << JxW_face[qp]*JxW_face[qp2]*mean_pressure*(qface_normals[qp](0)*phi_face[i][qp] * qface_normals[qp2](0)*phi_face[j][qp2]) << std::endl;
						//std::cout << "mean_pressure = " << mean_pressure << std::endl;
					      }
					    //std::cout << "qp = " << qp << ", qp2 = " << qp2 << std::endl;
					  }




				      }
				  }



				if (!convective_form
				    || es->parameters.get < bool >
				    ("neumann_stabilised")
				    || es->parameters.get < bool >
				    ("bertoglio_stabilisation"))
				  {
				    std::
				      cout << "in the wrong place..." << std::
				      endl;
				    //let us calculate the inflow normal velocity
				    double normal_velocity = 0.0;
				    double previous_normal_velocity = 0.0;
				    double inflow_normal_velocity_bdy = 0.0;
				    double outflow_normal_velocity_bdy = 0.0;
				    double approx_inflow_normal_velocity_bdy = 0.0;	// either interpolation or from prev timestep
				    //double approx_outflow_normal_velocity_bdy = 0.0;
				    Number u = 0., v = 0., w = 0.;
				    Number previous_u = 0., previous_v =
				      0., previous_w = 0.;
				    Number previous_previous_u =
				      0., previous_previous_v =
				      0., previous_previous_w = 0.;
				    Number approx_u = 0., approx_v =
				      0., approx_w = 0.;
				    double approx_normal_velocity = 0.0;

				    Gradient grad_u, grad_v, grad_w;

				    for (unsigned int l = 0; l < n_u_dofs;
					 l++)
				      {
					// From the previous Newton iterate:
					u +=
					  phi_face[l][qp] *
					  system->
					  current_solution (dof_indices_u[l]);
					v +=
					  phi_face[l][qp] *
					  system->
					  current_solution (dof_indices_v[l]);
					if (threed)
					  {
					    w +=
					      phi_face[l][qp] *
					      system->
					      current_solution (dof_indices_w
								[l]);
					  }


					grad_u.add_scaled (dphi_face[l][qp],
							   system->
							   current_solution
							   (dof_indices_u
							    [l]));
					grad_v.add_scaled (dphi_face[l][qp],
							   system->
							   current_solution
							   (dof_indices_v
							    [l]));
					if (threed)
					  grad_w.add_scaled (dphi_face[l][qp],
							     system->
							     current_solution
							     (dof_indices_w
							      [l]));

					// From the previous time_step:
					previous_u +=
					  phi_face[l][qp] *
					  system->
					  old_solution (dof_indices_u[l]);
					previous_v +=
					  phi_face[l][qp] *
					  system->
					  old_solution (dof_indices_v[l]);
					if (threed)
					  {
					    previous_w +=
					      phi_face[l][qp] *
					      system->
					      old_solution (dof_indices_w[l]);
					  }

					// From the previous previous time_step:
					previous_previous_u +=
					  phi_face[l][qp] *
					  system->
					  older_solution (dof_indices_u[l]);
					previous_previous_v +=
					  phi_face[l][qp] *
					  system->
					  older_solution (dof_indices_v[l]);
					if (threed)
					  {
					    previous_previous_w +=
					      phi_face[l][qp] *
					      system->
					      older_solution (dof_indices_w
							      [l]);
					  }

				      }

				    // linearly interpolate
				    if (false)	//es->parameters.get<bool>("neumann_stabilised_adjusted_interpolated"))
				      {
					approx_u =
					  previous_u + (previous_u -
							previous_previous_u) /
					  dt * dt;
					approx_v =
					  previous_v + (previous_v -
							previous_previous_v) /
					  dt * dt;
					if (threed)
					  {
					    approx_w =
					      previous_w + (previous_w -
							    previous_previous_w)
					      / dt * dt;
					  }
				      }
				    else
				      {
					approx_u = previous_u;
					approx_v = previous_v;
					if (threed)
					  {
					    approx_w = previous_w;
					  }
				      }




				    normal_velocity +=
				      qface_normals[qp] (0) * u;
				    normal_velocity +=
				      qface_normals[qp] (1) * v;
				    if (threed)
				      {
					normal_velocity +=
					  qface_normals[qp] (2) * w;
				      }


				    previous_normal_velocity +=
				      qface_normals[qp] (0) * previous_u;
				    previous_normal_velocity +=
				      qface_normals[qp] (1) * previous_v;
				    if (threed)
				      {
					previous_normal_velocity +=
					  qface_normals[qp] (2) * previous_w;
				      }

				    approx_normal_velocity +=
				      qface_normals[qp] (0) * approx_u;
				    approx_normal_velocity +=
				      qface_normals[qp] (1) * approx_v;
				    if (threed)
				      {
					approx_normal_velocity +=
					  qface_normals[qp] (2) * approx_w;
				      }
				    //std::cout << "normal = (" << qface_normals[qp](0) << "," << qface_normals[qp](1) << "," << qface_normals[qp](2) << ")" << std::endl;

				    inflow_normal_velocity_bdy =
				      (normal_velocity -
				       fabs (normal_velocity)) / 2.0;
				    outflow_normal_velocity_bdy =
				      (normal_velocity +
				       fabs (normal_velocity)) / 2.0;


				    approx_inflow_normal_velocity_bdy =
				      (approx_normal_velocity -
				       fabs (approx_normal_velocity)) / 2.0;
				    //approx_outflow_normal_velocity_bdy = (approx_normal_velocity + fabs(approx_normal_velocity))/2.0;

				    //says whether is on the correct boundary, inflow for convective form and outflow for conservative form
				    int bdy_bool = 0;
				    if (fabs (normal_velocity) > 1e-10)
				      {
					if (convective_form)
					  bdy_bool =
					    ((-inflow_normal_velocity_bdy /
					      fabs (normal_velocity)) + 0.5);
					else
					  bdy_bool =
					    (outflow_normal_velocity_bdy /
					     fabs (normal_velocity) + 0.5);
				      }

				    //says whether is on the correct boundary, inflow for convective form and outflow for conservative form
				    int previous_bdy_bool = 0;
				    if (fabs (previous_normal_velocity) >
					1e-10)
				      {
					previous_bdy_bool =
					  ((-
					    (previous_normal_velocity -
					     fabs (previous_normal_velocity))
					    / 2.0 /
					    fabs (previous_normal_velocity)) +
					   0.5);
				      }

				    int approx_inflow_bdy_bool = 0;
				    if (fabs (approx_normal_velocity) > 1e-10)
				      {
					approx_inflow_bdy_bool =
					  ((-approx_inflow_normal_velocity_bdy
					    / fabs (approx_normal_velocity)) +
					   0.5);
				      }
				    //if conservative and not neumann stabilised then we always have the terms from int by parts
				    if (!convective_form
					&& !es->parameters.get < bool >
					("neumann_stabilised"))
				      {
					bdy_bool = 1;
				      }


				    double r =
				      (*surface_boundaries)[boundary_id]->
				      get_normalised_distance_from_centroid
				      (qpoint_face[qp]);
				    //double r = 1.0;//
				    double normalisation_constant =
				      (*surface_boundaries)[boundary_id]->
				      get_unit_parabola_integral ();
				    //double normalisation_constant = 1.0;//
				    double radius = 1.0;

				    //approx_u = 20.;
				    //approx_v = 0.;
				    //approx_w = 0.;
				    //approx_normal_velocity = 20.;

				    //std::cout << "approx normal velocity = " <<  approx_normal_velocity << std::endl;


				    // hmmm small value changes appear to be making a difference
				    for (unsigned int i = 0; i < n_u_dofs;
					 i++)
				      {
					if (convective_form)
					  {
					    // if in convective form we only add a term if doing the inflow stabilisation
					    if (es->parameters.get < bool >
						("neumann_stabilised"))
					      {
						for (unsigned int j = 0;
						     j < n_u_dofs; j++)
						  {
						    // - \int u_n^{in} * w * u
						    if (!es->parameters.get <
							bool >
							("neumann_stabilised_linear"))
						      {
							Kuu (i, j) +=
							  -backflow_stab_param
							  * JxW_face[qp] *
							  normal_velocity *
							  bdy_bool *
							  (phi_face[i][qp] *
							   phi_face[j][qp]);
							Kvv (i, j) +=
							  -backflow_stab_param
							  * JxW_face[qp] *
							  normal_velocity *
							  bdy_bool *
							  (phi_face[i][qp] *
							   phi_face[j][qp]);
							if (threed)
							  {
							    Kww (i, j) +=
							      -backflow_stab_param
							      * JxW_face[qp] *
							      normal_velocity
							      * bdy_bool *
							      (phi_face[i][qp]
							       *
							       phi_face[j]
							       [qp]);
							  }
						      }
						    else
						      {
							Kuu (i, j) +=
							  -backflow_stab_param
							  * JxW_face[qp] *
							  previous_normal_velocity
							  *
							  previous_bdy_bool *
							  (phi_face[i][qp] *
							   phi_face[j][qp]);
							Kvv (i, j) +=
							  -backflow_stab_param
							  * JxW_face[qp] *
							  previous_normal_velocity
							  *
							  previous_bdy_bool *
							  (phi_face[i][qp] *
							   phi_face[j][qp]);
							if (threed)
							  {
							    Kww (i, j) +=
							      -backflow_stab_param
							      * JxW_face[qp] *
							      previous_normal_velocity
							      *
							      previous_bdy_bool
							      *
							      (phi_face[i][qp]
							       *
							       phi_face[j]
							       [qp]);
							  }
						      }

						    if (newton
							&& !es->parameters.
							get < bool >
							("neumann_stabilised_linear"))
						      {
							Kuu (i, j) +=
							  -backflow_stab_param
							  * JxW_face[qp] *
							  bdy_bool * (u *
								      qface_normals
								      [qp] (0)
								      *
								      phi_face
								      [i][qp]
								      *
								      phi_face
								      [j]
								      [qp]);
							Kuv (i, j) +=
							  -backflow_stab_param
							  * JxW_face[qp] *
							  bdy_bool * (u *
								      qface_normals
								      [qp] (1)
								      *
								      phi_face
								      [i][qp]
								      *
								      phi_face
								      [j]
								      [qp]);
							if (threed)
							  {
							    Kuw (i, j) +=
							      -backflow_stab_param
							      * JxW_face[qp] *
							      bdy_bool * (u *
									  qface_normals
									  [qp]
									  (2)
									  *
									  phi_face
									  [i]
									  [qp]
									  *
									  phi_face
									  [j]
									  [qp]);
							  }
							Kvu (i, j) +=
							  -backflow_stab_param
							  * JxW_face[qp] *
							  bdy_bool * (v *
								      qface_normals
								      [qp] (0)
								      *
								      phi_face
								      [i][qp]
								      *
								      phi_face
								      [j]
								      [qp]);
							Kvv (i, j) +=
							  -backflow_stab_param
							  * JxW_face[qp] *
							  bdy_bool * (v *
								      qface_normals
								      [qp] (1)
								      *
								      phi_face
								      [i][qp]
								      *
								      phi_face
								      [j]
								      [qp]);
							if (threed)
							  {
							    Kvw (i, j) +=
							      -backflow_stab_param
							      * JxW_face[qp] *
							      bdy_bool * (v *
									  qface_normals
									  [qp]
									  (2)
									  *
									  phi_face
									  [i]
									  [qp]
									  *
									  phi_face
									  [j]
									  [qp]);
							  }
							if (threed)
							  {
							    Kwu (i, j) +=
							      -backflow_stab_param
							      * JxW_face[qp] *
							      bdy_bool * (w *
									  qface_normals
									  [qp]
									  (0)
									  *
									  phi_face
									  [i]
									  [qp]
									  *
									  phi_face
									  [j]
									  [qp]);
							  }
							if (threed)
							  {
							    Kwv (i, j) +=
							      -backflow_stab_param
							      * JxW_face[qp] *
							      bdy_bool * (w *
									  qface_normals
									  [qp]
									  (1)
									  *
									  phi_face
									  [i]
									  [qp]
									  *
									  phi_face
									  [j]
									  [qp]);
							  }
							if (threed)
							  {
							    Kww (i, j) +=
							      -backflow_stab_param
							      * JxW_face[qp] *
							      bdy_bool * (w *
									  qface_normals
									  [qp]
									  (2)
									  *
									  phi_face
									  [i]
									  [qp]
									  *
									  phi_face
									  [j]
									  [qp]);
							  }
						      }

						  }

						// another newton term
						if (newton
						    && !es->parameters.get <
						    bool >
						    ("neumann_stabilised_linear"))
						  {
						    // \int (u_old \cdot w) * (u_n^{in})
						    Fu (i) +=
						      -backflow_stab_param *
						      JxW_face[qp] *
						      normal_velocity *
						      bdy_bool * u *
						      phi_face[i][qp];
						    Fv (i) +=
						      -backflow_stab_param *
						      JxW_face[qp] *
						      normal_velocity *
						      bdy_bool * v *
						      phi_face[i][qp];
						    if (threed)
						      {
							Fw (i) +=
							  -backflow_stab_param
							  * JxW_face[qp] *
							  normal_velocity *
							  bdy_bool * w *
							  phi_face[i][qp];
						      }
						  }

						//we may also want to adjust the mean pressure boundary condition based on the normal velocity
						//well, of the previous timestep to be easy and avoid issues
						if ((es->parameters.get <
						     bool >
						     ("neumann_stabilised_adjusted")
						     || es->parameters.get <
						     bool >
						     ("neumann_stabilised_adjusted_interpolated"))
						    && approx_inflow_bdy_bool
						    > 0)
						  {

						    Fu (i) +=
						      -backflow_stab_param *
						      JxW_face[qp] *
						      approx_normal_velocity *
						      (pow (radius, 2) -
						       pow (r,
							    2)) /
						      normalisation_constant *
						      interp_flow_bc_value
						      [boundary_id] *
						      qface_normals[qp] (0) *
						      approx_inflow_bdy_bool *
						      phi_face[i][qp];
						    Fv (i) +=
						      -backflow_stab_param *
						      JxW_face[qp] *
						      approx_normal_velocity *
						      (pow (radius, 2) -
						       pow (r,
							    2)) /
						      normalisation_constant *
						      interp_flow_bc_value
						      [boundary_id] *
						      qface_normals[qp] (1) *
						      approx_inflow_bdy_bool *
						      phi_face[i][qp];
						    if (threed)
						      {
							Fw (i) +=
							  -backflow_stab_param
							  * JxW_face[qp] *
							  approx_normal_velocity
							  * (pow (radius, 2) -
							     pow (r,
								  2)) /
							  normalisation_constant
							  *
							  interp_flow_bc_value
							  [boundary_id] *
							  qface_normals[qp]
							  (2) *
							  approx_inflow_bdy_bool
							  * phi_face[i][qp];
						      }
						  }
						//Fu(i) += -JxW_face[qp]*(((phi_face[i][qp]*qface_normals[qp](0))*mean_pressure));
						//Fv(i) += -JxW_face[qp]*(((phi_face[i][qp]*qface_normals[qp](1))*mean_pressure));
						//Fw(i) += -JxW_face[qp]*(((phi_face[i][qp]*qface_normals[qp](2))*mean_pressure));
					      }
					    else if (es->parameters.get <
						     bool >
						     ("bertoglio_stabilisation"))
					      {
						//std::cout << "holla, num_tangents = " << qface_tangents[qp].size() << std::endl;
						//std::cout << "tan 1 = " << qface_tangents[qp][0] << std::endl;
						//std::cout << "tan 2 = " << qface_tangents[qp][1] << std::endl;
						//std::cout << "bertoglio_stab_param = " << bertoglio_stab_param << std::endl;

						// weird, in 2D there are still 2 tangents, but we only want the first one
						unsigned int num_tangents =
						  qface_tangents[qp].size ();
						if (!threed)
						  num_tangents = 1;

						for (unsigned int k = 0;
						     k < num_tangents; k++)
						  {
						    for (unsigned int j = 0;
							 j < n_u_dofs; j++)
						      {
							// - \int u_n^{in} * w * u

// gradients in tangential directions
/*
															Kuu(i,j) += -bertoglio_stab_param*JxW_face[qp]*previous_normal_velocity*previous_bdy_bool*(dphi_face[i][qp]*qface_tangents[qp][k]
																																																					*dphi_face[j][qp]*qface_tangents[qp][k]);
															Kvv(i,j) += -bertoglio_stab_param*JxW_face[qp]*previous_normal_velocity*previous_bdy_bool*(dphi_face[i][qp]*qface_tangents[qp][k]
																																																					*dphi_face[j][qp]*qface_tangents[qp][k]);
															if(threed) { Kww(i,j) += -bertoglio_stab_param*JxW_face[qp]*previous_normal_velocity*previous_bdy_bool*(dphi_face[i][qp]*qface_tangents[qp][k]
																																																					*dphi_face[j][qp]*qface_tangents[qp][k]); }
*/

// velocities in tangential directions
							Kuu (i, j) +=
							  -bertoglio_stab_param
							  * JxW_face[qp] *
							  previous_normal_velocity
							  *
							  previous_bdy_bool *
							  (dphi_face[i][qp] *
							   qface_tangents[qp]
							   [k] (0) *
							   dphi_face[j][qp] *
							   qface_tangents[qp]
							   [k] (0));
							Kvv (i, j) +=
							  -bertoglio_stab_param
							  * JxW_face[qp] *
							  previous_normal_velocity
							  *
							  previous_bdy_bool *
							  (dphi_face[i][qp] *
							   qface_tangents[qp]
							   [k] (1) *
							   dphi_face[j][qp] *
							   qface_tangents[qp]
							   [k] (1));
							if (threed)
							  {
							    Kww (i, j) +=
							      -bertoglio_stab_param
							      * JxW_face[qp] *
							      previous_normal_velocity
							      *
							      previous_bdy_bool
							      *
							      (dphi_face[i]
							       [qp] *
							       qface_tangents
							       [qp][k] (2) *
							       dphi_face[j]
							       [qp] *
							       qface_tangents
							       [qp][k] (2));
							  }



/*
															Kuu(i,j) += bertoglio_stab_param*JxW_face[qp]*normal_velocity*bdy_bool*(dphi_face[i][qp]*qface_tangents[qp][k]
																																																					*dphi_face[j][qp]*qface_tangents[qp][k]);
															Kvv(i,j) += bertoglio_stab_param*JxW_face[qp]*normal_velocity*bdy_bool*(dphi_face[i][qp]*qface_tangents[qp][k]
																																																					*dphi_face[j][qp]*qface_tangents[qp][k]);
															if(threed) { Kww(i,j) += bertoglio_stab_param*JxW_face[qp]*normal_velocity*bdy_bool*(dphi_face[i][qp]*qface_tangents[qp][k]
																																																					*dphi_face[j][qp]*qface_tangents[qp][k]); }
*/

/*
															if(newton)
															{
																Kuu(i,j) += bertoglio_stab_param*JxW_face[qp]*bdy_bool
																						* (qface_normals[qp](0)*phi_face[j][qp]* grad_u*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k]);
																Kuv(i,j) += bertoglio_stab_param*JxW_face[qp]*bdy_bool
																						* (qface_normals[qp](1)*phi_face[j][qp]* grad_u*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k]);                
																if(threed) { Kuw(i,j) += bertoglio_stab_param*JxW_face[qp]*bdy_bool
																						* (qface_normals[qp](2)*phi_face[j][qp]* grad_u*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k]); }
																Kvu(i,j) += bertoglio_stab_param*JxW_face[qp]*bdy_bool
																						* (qface_normals[qp](0)*phi_face[j][qp]* grad_v*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k]);
																Kvv(i,j) += bertoglio_stab_param*JxW_face[qp]*bdy_bool
																						* (qface_normals[qp](1)*phi_face[j][qp]* grad_v*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k]);                
																if(threed) { Kvw(i,j) += bertoglio_stab_param*JxW_face[qp]*bdy_bool
																						* (qface_normals[qp](2)*phi_face[j][qp]* grad_v*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k]); }             
																if(threed) { Kwu(i,j) += bertoglio_stab_param*JxW_face[qp]*bdy_bool
																						* (qface_normals[qp](0)*phi_face[j][qp]* grad_w*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k]); }              
																if(threed) { Kwv(i,j) += bertoglio_stab_param*JxW_face[qp]*bdy_bool
																						* (qface_normals[qp](1)*phi_face[j][qp]* grad_w*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k]); }
																if(threed) { Kww(i,j) += bertoglio_stab_param*JxW_face[qp]*bdy_bool
																						* (qface_normals[qp](2)*phi_face[j][qp]* grad_w*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k]); }
															}*/



						      }

						    // another newton term
						    if (newton)
						      {
/*
															// \int (u_old \cdot w) * (u_n^{in})
															Fu(i) += bertoglio_stab_param*JxW_face[qp]*normal_velocity*bdy_bool* grad_u*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k];
															Fv(i) += bertoglio_stab_param*JxW_face[qp]*normal_velocity*bdy_bool* grad_v*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k];
															if(threed) { Fw(i) += bertoglio_stab_param*JxW_face[qp]*normal_velocity*bdy_bool* grad_w*qface_tangents[qp][k]
																							* dphi_face[i][qp]*qface_tangents[qp][k]; }
*/
						      }


						  }
					      }
					  }
					else
					  {

					    //always add these terms stabilised or not - actually add 1.0 * stab terms if on bdy and (1.0 - stab_param) * stab terms if not
					    // fix this another time though
					    for (unsigned int j = 0;
						 j < n_u_dofs; j++)
					      {
						// + \int u_n^{out} * w * u
						Kuu (i, j) +=
						  backflow_stab_param *
						  JxW_face[qp] *
						  normal_velocity * bdy_bool *
						  (phi_face[i][qp] *
						   phi_face[j][qp]);
						Kvv (i, j) +=
						  backflow_stab_param *
						  JxW_face[qp] *
						  normal_velocity * bdy_bool *
						  (phi_face[i][qp] *
						   phi_face[j][qp]);
						if (threed)
						  {
						    Kww (i, j) +=
						      backflow_stab_param *
						      JxW_face[qp] *
						      normal_velocity *
						      bdy_bool *
						      (phi_face[i][qp] *
						       phi_face[j][qp]);
						  }

						if (newton)
						  {
						    Kuu (i, j) +=
						      backflow_stab_param *
						      JxW_face[qp] *
						      bdy_bool * (u *
								  qface_normals
								  [qp] (0) *
								  phi_face[i]
								  [qp] *
								  phi_face[j]
								  [qp]);
						    Kuv (i, j) +=
						      backflow_stab_param *
						      JxW_face[qp] *
						      bdy_bool * (u *
								  qface_normals
								  [qp] (1) *
								  phi_face[i]
								  [qp] *
								  phi_face[j]
								  [qp]);
						    if (threed)
						      {
							Kuw (i, j) +=
							  backflow_stab_param
							  * JxW_face[qp] *
							  bdy_bool * (u *
								      qface_normals
								      [qp] (2)
								      *
								      phi_face
								      [i][qp]
								      *
								      phi_face
								      [j]
								      [qp]);
						      }
						    Kvu (i, j) +=
						      backflow_stab_param *
						      JxW_face[qp] *
						      bdy_bool * (v *
								  qface_normals
								  [qp] (0) *
								  phi_face[i]
								  [qp] *
								  phi_face[j]
								  [qp]);
						    Kvv (i, j) +=
						      backflow_stab_param *
						      JxW_face[qp] *
						      bdy_bool * (v *
								  qface_normals
								  [qp] (1) *
								  phi_face[i]
								  [qp] *
								  phi_face[j]
								  [qp]);
						    if (threed)
						      {
							Kvw (i, j) +=
							  backflow_stab_param
							  * JxW_face[qp] *
							  bdy_bool * (v *
								      qface_normals
								      [qp] (2)
								      *
								      phi_face
								      [i][qp]
								      *
								      phi_face
								      [j]
								      [qp]);
						      }
						    if (threed)
						      {
							Kwu (i, j) +=
							  backflow_stab_param
							  * JxW_face[qp] *
							  bdy_bool * (w *
								      qface_normals
								      [qp] (0)
								      *
								      phi_face
								      [i][qp]
								      *
								      phi_face
								      [j]
								      [qp]);
						      }
						    if (threed)
						      {
							Kwv (i, j) +=
							  backflow_stab_param
							  * JxW_face[qp] *
							  bdy_bool * (w *
								      qface_normals
								      [qp] (1)
								      *
								      phi_face
								      [i][qp]
								      *
								      phi_face
								      [j]
								      [qp]);
						      }
						    if (threed)
						      {
							Kww (i, j) +=
							  backflow_stab_param
							  * JxW_face[qp] *
							  bdy_bool * (w *
								      qface_normals
								      [qp] (2)
								      *
								      phi_face
								      [i][qp]
								      *
								      phi_face
								      [j]
								      [qp]);
						      }
						  }

					      }

					    // another newton term
					    if (newton)
					      {
						// \int (u_old \cdot w) * (u_n^{in})
						Fu (i) +=
						  backflow_stab_param *
						  JxW_face[qp] *
						  normal_velocity * bdy_bool *
						  u * phi_face[i][qp];
						Fv (i) +=
						  backflow_stab_param *
						  JxW_face[qp] *
						  normal_velocity * bdy_bool *
						  v * phi_face[i][qp];
						if (threed)
						  {
						    Fw (i) +=
						      backflow_stab_param *
						      JxW_face[qp] *
						      normal_velocity *
						      bdy_bool * w *
						      phi_face[i][qp];
						  }
					      }

					    //we may also want to adjust the mean pressure boundary condition based on the normal velocity
					    //well, of the previous timestep to be easy and avoid issues
					    if (es->parameters.get < bool >
						("neumann_stabilised_adjusted"))
					      {
						Fu (i) +=
						  backflow_stab_param *
						  JxW_face[qp] *
						  approx_normal_velocity *
						  approx_u *
						  approx_inflow_bdy_bool *
						  phi_face[i][qp];
						Fv (i) +=
						  backflow_stab_param *
						  JxW_face[qp] *
						  approx_normal_velocity *
						  approx_v *
						  approx_inflow_bdy_bool *
						  phi_face[i][qp];
						if (threed)
						  {
						    Fw (i) +=
						      backflow_stab_param *
						      JxW_face[qp] *
						      approx_normal_velocity *
						      approx_w *
						      approx_inflow_bdy_bool *
						      phi_face[i][qp];
						  }
					      }
					  }
				      }
				  }

			      }	//end face quad loop


			    // moghadam crap
			    /*
			       for (unsigned int i=0; i<n_u_dofs; i++)
			       {
			       for (unsigned int j=0; j<n_u_dofs; j++)
			       {

			       Kuu(i,j) += mean_pressure*temp_vector[i] * temp_vector[j];
			       temp_matrix_2[i][j] += mean_pressure*temp_vector[i] * temp_vector[j];
			       actual_terms_2 += mean_pressure*temp_vector[i] * temp_vector[j];

			       }
			       }

			       std::cout << "extra_terms = " << extra_terms << std::endl;
			       std::cout << "actual_terms = " << actual_terms << std::endl;
			       std::cout << "actual_terms_2 = " << actual_terms_2 << std::endl;

			       std::cout << "temp_matrix_1:" << std::endl;
			       for(unsigned int i=0; i<n_u_dofs; i++)
			       {
			       std::cout << "row " << i;
			       for(unsigned int j=0; j<n_u_dofs; j++)
			       std::cout << " " << temp_matrix_1[i][j];
			       std::cout << std::endl;
			       }


			       std::cout << "temp_matrix_2:" << std::endl;
			       for(unsigned int i=0; i<n_u_dofs; i++)
			       {
			       std::cout << "row " << i;
			       for(unsigned int j=0; j<n_u_dofs; j++)
			       std::cout << " " << temp_matrix_2[i][j];
			       std::cout << std::endl;
			       }
			     */


			  }	//end if(bc_type[boundary_id].compare("pressure") == 0)

			// mean-flow conditions/terms
			else if (bc_type[boundary_id].compare ("mean-flow") ==
				 0)
			  {

			    //NOT IMPLEMENTED

			  }	//end if(bc_type[boundary_id].compare("mean-flow") == 0)
		      }






		    // ************** COUPLING BOUNDARY TERMS ********************* //
		    // construct stuff to put into 1d matrix lines
		    else
		      {

			// mean pressure conditions/terms
			if (bc_type[boundary_id].compare ("pressure") == 0)
			  {
			    //set the dof_index_0 to the correct node number
			    // primary pressure boundary node comes from either daughter 1 or daughter 2 (so we choose daughter 1)
			    // primary flux node comes from daughter 1
			    pressure_1d_index =
			      primary_pressure_boundary_nodes_1d[boundary_id];
			    flux_1d_index =
			      primary_flux_boundary_nodes_1d[boundary_id];

			    fe_vel_face->reinit (elem, s);
			    fe_pres_face->reinit (elem, s);

			    std::vector < Real > JxW_face = JxW_face_elem;
			    //if axisym then need to multiply all integrals by "r" i.e. y = qpoint[qp](1)
			    if (es->parameters.get <
				unsigned int >("geometry_type") == 5)
			      {
				for (unsigned int qp = 0;
				     qp < qface.n_points (); qp++)
				  JxW_face[qp] *= qpoint_face[qp] (1);
			      }


			    if (multiply_system_by_dt)
			      {
				for (unsigned int qp = 0;
				     qp < qface.n_points (); qp++)
				  {
				    JxW_face[qp] *= dt;
				  }
			      }

			    // for the pressure term in the 3d equations
			    for (unsigned int qp = 0; qp < qface.n_points ();
				 qp++)
			      {
				for (unsigned int i = 0; i < n_u_dofs; i++)
				  {
				    Fu_coupled_p (i) +=
				      JxW_face[qp] * phi_face[i][qp] *
				      qface_normals[qp] (0);
				    Fv_coupled_p (i) +=
				      JxW_face[qp] * phi_face[i][qp] *
				      qface_normals[qp] (1);
				    if (threed)
				      {
					Fw_coupled_p (i) +=
					  JxW_face[qp] * phi_face[i][qp] *
					  qface_normals[qp] (2);
				      }
				  }

			      }	//end face quad loop

			    // for the flux condition in the 1d equations

			    // for the pressure term in the 3d equations
			    // EQN IS + dt*Q_1d - dt*Q_3d = 0
			    // later on we eliminate the terms that correspond to the 
			    for (unsigned int qp = 0; qp < qface.n_points ();
				 qp++)
			      {
				for (unsigned int i = 0; i < n_u_dofs; i++)
				  {
				    Fu_coupled_u (i) +=
				      JxW_face[qp] * phi_face[i][qp] *
				      qface_normals[qp] (0);
				    Fv_coupled_u (i) +=
				      JxW_face[qp] * phi_face[i][qp] *
				      qface_normals[qp] (1);
				    if (threed)
				      {
					Fw_coupled_u (i) +=
					  JxW_face[qp] * phi_face[i][qp] *
					  qface_normals[qp] (2);
				      }

				  }
				//std::cout << "normal = " << qface_normals[qp] << std::endl;

			      }	//end face quad loop


			  }
		      }		// end coupling terms

		  }

	      }			//end side loop






	    // ******************* PIN THE PRESSURE **************************** //
	    // Pin the pressure to zero at global node number "pressure_node".
	    // This effectively removes the non-trivial null space of constant
	    // pressure solutions.

	    // need to pin pressure for lid driven cavity
	    // get dof associated with node to pin
	    // should work for lagrange linear

	    /*
	       if (pin_pressure)
	       {
	       const Real p_value               = 0.0;
	       for (unsigned int c=0; c<elem->n_nodes(); c++)
	       {
	       if (elem->node(c) == pressure_node)
	       {
	       pressure_dof = dof_indices_p[c];
	       // hmmm forgot what it was..
	       }
	       }
	       }
	     */

	    if (pin_pressure)
	      {
		const unsigned int pressure_node = 0;
		const Real p_value = 0.0;
		for (unsigned int c = 0; c < elem->n_nodes (); c++)
		  {
		    if (elem->node (c) == pressure_node)
		      {
			Kpp (c, c) += penalty;
			Fp (c) += penalty * p_value;
			if (es->parameters.get <
			    unsigned int >("preconditioner_type_3d"))
			  {
			    Kpp_pre_laplacian (c, c) += penalty * Re;
			    Kpp_pre_convection_diffusion (c, c) += penalty;
			  }
		      }
		  }
	      }



	  }			// end boundary condition section






	  // copy entries from Ke to Ke_pre_velocity
	  if ((es->parameters.get < unsigned int >("preconditioner_type_3d")
	       || es->parameters.get <
	       unsigned int >("preconditioner_type_3d1d")) &&(es->parameters.
							      get <
							      unsigned int
							      >
							      ("preconditioner_type_3d1d")
							      != 10
							      && es->
							      parameters.get <
							      unsigned int
							      >
							      ("preconditioner_type_3d1d")
							      != 11))
	    {
	      Ke_pre_velocity = Ke;

	      if (es->parameters.get <
		  unsigned int >("preconditioner_type_3d1d") != 10)
		{
		  // zero out the pressure bits
		  Kup_pre_velocity.zero ();
		  Kvp_pre_velocity.zero ();
		  if (threed)
		    Kwp_pre_velocity.zero ();

		  Kpu_pre_velocity.zero ();
		  Kpv_pre_velocity.zero ();
		  if (threed)
		    Kpw_pre_velocity.zero ();

		  Kpp_pre_velocity.zero ();
		}


	    }
	  else if (es->parameters.get <
		   unsigned int >("preconditioner_type_3d1d") == 10)
	    {
	      // don't do anything, the stokes part has been assembled
	    }
	  else if (es->parameters.get <
		   unsigned int >("preconditioner_type_3d1d") == 11)
	    {

	      // stokes part has been assembled, but need to zero out the pressure bits
	      Kup_pre_velocity.zero ();
	      Kvp_pre_velocity.zero ();
	      if (threed)
		Kwp_pre_velocity.zero ();

	      Kpu_pre_velocity.zero ();
	      Kpv_pre_velocity.zero ();
	      if (threed)
		Kpw_pre_velocity.zero ();

	      Kpp_pre_velocity.zero ();
	    }














	  // ********* ASSEMBLE AND APPLY DIRICHLET BOUNDARY CONDITIONS ******** //
	  // if we are estimating the error then calcalate relative contributions
	  if (!estimating_error)
	    {

	      // ***************** INSERT COUPLING BOUNDARY CONDITIONS ********* //
	      // if we have pressure coupled, then we need and found a boundary
	      // the only problem is that we don't want the dofs that are on the 
	      // boundary to be affected. so we zero them out and assume that there
	      // is zero contribution from them, assumming a zero boundary condition
	      // this is fairly common especially in my applications.
	      if (pressure_coupled && pressure_1d_index != 0)
		{

		  // if one of the dofs is a 
		  // constrained dof then we zero it. otherwise we add it to the matrix
		  // problem is that the pressure integral should contain stuff from the boundary
		  // maybe it wil average? who knows... seems like it, but it isn't nice.                         
		  for (unsigned int i = 0; i < dof_indices.size (); i++)
		    {
		      //if not constrained dof
		      if (!dof_map.is_constrained_dof (dof_indices[i]))
			{
			  system->matrix->add (dof_indices[i],
					       pressure_1d_index,
					       Fe_coupled_p (i));
			  if (es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") == 7
			      || es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") == 8
			      || es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") == 9
			      || es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") == 10
			      || es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") == 11
			      || es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") ==
			      12)
			    system->get_matrix ("Preconditioner").
			      add (dof_indices[i], pressure_1d_index,
				   Fe_coupled_p (i));

			}
		    }




		  //flux coupling - only in row corresponding to 1d dof is how it is defined
		  // this 1d dof is now pressure 0 so that getting towards symmetry

		  // if one of the dofs is a 
		  // constrained dof then we zero it. otherwise we add it to the matrix
		  for (unsigned int i = 0; i < dof_indices.size (); i++)
		    {
		      //if not constrained dof
		      if (!dof_map.is_constrained_dof (dof_indices[i]))
			{
			  system->matrix->add (flux_1d_index, dof_indices[i],
					       Fe_coupled_u (i));
			  if (es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") == 7
			      || es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") == 8
			      || es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") == 9
			      || es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") == 10
			      || es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") == 11
			      || es->parameters.get <
			      unsigned int >("preconditioner_type_3d1d") ==
			      12)
			    system->get_matrix ("Preconditioner").
			      add (flux_1d_index, dof_indices[i],
				   Fe_coupled_u (i));
			}
		    }

		  // reset to zero so that isn't used in next loop.
		  pressure_1d_index = 0;
		  flux_1d_index = 0;
		}






	      // ************** INSERT LOCAL MATRICES INTO BIG MATRIX ************ //
	      // apply dirichlet conditions (e.g. wall and maybe inflow)
	      //hmm really not sure bout this last argument but it works
	      //dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices,false);
	      dof_map.heterogenously_constrain_element_matrix_and_vector (Ke,
									  Fe,
									  dof_indices,
									  true);
	      dof_map.constrain_element_matrix (Ke_pre_velocity_mass,
						dof_indices, false);
	      //dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices,false);                              

	      system->matrix->add_matrix (Ke, dof_indices);
	      system->rhs->add_vector (Fe, dof_indices);
	      if (es->parameters.get <
		  unsigned int >("preconditioner_type_3d")
		  || es->parameters.get <
		  unsigned int >("preconditioner_type_3d1d"))
		{
		  system->get_matrix ("Pressure Mass Matrix").
		    add_matrix (Ke_pre_mass, dof_indices);
		  system->get_matrix ("Pressure Laplacian Matrix").
		    add_matrix (Ke_pre_laplacian, dof_indices);
		  system->get_matrix ("Pressure Convection Diffusion Matrix").
		    add_matrix (Ke_pre_convection_diffusion, dof_indices);
		  system->get_matrix ("Velocity Mass Matrix").
		    add_matrix (Ke_pre_velocity_mass, dof_indices);
		  system->get_matrix ("Preconditioner").add_matrix (Ke,
								    dof_indices);

		  // don't use any boundary conditions... for the stokes preconditioners only assemble once
		  dof_map.constrain_element_matrix (Ke_pre_velocity,
						    dof_indices, true);
		  system->get_matrix ("Velocity Matrix").
		    add_matrix (Ke_pre_velocity, dof_indices);


		}

	      // ***************************************************************** //










	      // ****************** APPLY PRECONDITIONER BOUNDARY CONDITIONS ***** //
	      // if we have a pressure dof that needs to be constrained then make all
	      // row/ column zero apart from diag
	      /*
	         if(pin_pressure)
	         {
	         if(pressure_dof >= 0)
	         {
	         std::cout << "pinning pressure" << std::endl;
	         for(unsigned int i=0; i<dof_indices.size(); i++)
	         {
	         system->matrix->set(pressure_dof,dof_indices[i],0.);
	         system->matrix->set(dof_indices[i],pressure_dof,0.);
	         if(es->parameters.get<unsigned int>("preconditioner_type"))
	         {
	         system->request_matrix("Preconditioner")->set(pressure_dof,dof_indices[i],0.);
	         system->request_matrix("Preconditioner")->set(dof_indices[i],pressure_dof,0.);
	         system->request_matrix("Pressure Mass Matrix")->set(pressure_dof,dof_indices[i],0.);
	         system->request_matrix("Pressure Mass Matrix")->set(dof_indices[i],pressure_dof,0.);
	         system->request_matrix("Pressure Laplacian Matrix")->set(pressure_dof,dof_indices[i],0.);
	         system->request_matrix("Pressure Laplacian Matrix")->set(dof_indices[i],pressure_dof,0.);
	         system->request_matrix("Pressure Convection Diffusion Matrix")->set(pressure_dof,dof_indices[i],0.);
	         system->request_matrix("Pressure Convection Diffusion Matrix")->set(dof_indices[i],pressure_dof,0.);
	         }
	         }
	         system->matrix->set(pressure_dof,pressure_dof,1.);
	         system->rhs->set(pressure_dof,0.);

	         if(es->parameters.get<unsigned int>("preconditioner_type"))
	         {
	         system->request_matrix("Preconditioner")->set(pressure_dof,pressure_dof,1.);
	         system->request_matrix("Pressure Mass Matrix")->set(pressure_dof,pressure_dof,-1.);
	         system->request_matrix("Pressure Laplacian Matrix")->set(pressure_dof,pressure_dof,1.);
	         system->request_matrix("Pressure Convection Diffusion Matrix")->set(pressure_dof,pressure_dof,1.);
	         }
	         }
	         }
	       */


	      // now we need to add the pressure dof constraints that were found to be on the inflow dirichlet boundary
	      // doesn't really matter if there is overlap.

	      if ((es->parameters.get <
		   unsigned int >("preconditioner_type_3d") == 4
		   || es->parameters.get <
		   unsigned int >("preconditioner_type_3d") ==
		   5) &&es->parameters.get < unsigned int >("problem_type") !=
		  4)
		{
		  if (es->parameters.get <
		      unsigned int >("pcd_boundary_condition_type") == 1
		      || es->parameters.get <
		      unsigned int >("pcd_boundary_condition_type") == 3
		      || es->parameters.get <
		      unsigned int >("pcd_boundary_condition_type") == 4)
		    {
		      for (unsigned int i = 0;
			   i < pressure_dofs_on_inflow_boundary.size (); i++)
			{
			  unsigned int local_pressure_dof = -1;
			  //find what local pressure dof in Ke_pre_convection_diffusion we need to take as the increment of the diagonal scaling
			  for (unsigned int j = 0; j < dof_indices.size ();
			       j++)
			    {
			      if (pressure_dofs_on_inflow_boundary[i] ==
				  dof_indices[j])
				{
				  local_pressure_dof = j;
				  break;
				}
			    }

			  // we found the correct pressure dof...
			  if (local_pressure_dof >= 0)
			    {
			      //std::cout << "hi" << std::endl;
			      double diagonal_scaling_increment =
				Ke_pre_convection_diffusion
				(local_pressure_dof, local_pressure_dof);
			      sum_fp_diag += diagonal_scaling_increment;
			      //std::cout << "pinning pressure dof " << pressure_dofs_on_inflow_boundary[i] << " on elem " << elem->id() << " on inflow bdy with scaling = " << diagonal_scaling_increment << std::endl;
			      for (unsigned int j = 0;
				   j < dof_indices.size (); j++)
				{
				  if (pressure_dofs_on_inflow_boundary[i] !=
				      dof_indices[j])
				    {
				      system->
					request_matrix
					("Pressure Laplacian Matrix")->
					set (pressure_dofs_on_inflow_boundary
					     [i], dof_indices[j], 0.);
				      system->
					request_matrix
					("Pressure Laplacian Matrix")->
					set (dof_indices[j],
					     pressure_dofs_on_inflow_boundary
					     [i], 0.);
				      system->
					request_matrix
					("Pressure Convection Diffusion Matrix")->
					set (pressure_dofs_on_inflow_boundary
					     [i], dof_indices[j], 0.);
				      system->
					request_matrix
					("Pressure Convection Diffusion Matrix")->
					set (dof_indices[j],
					     pressure_dofs_on_inflow_boundary
					     [i], 0.);
				    }
				}

			      if (es->parameters.get <
				  unsigned int
				  >("pcd_boundary_condition_type") == 1)
				{
				  // first subject the increment we in the matrix construction of this element
				  system->
				    request_matrix
				    ("Pressure Laplacian Matrix")->
				    add (pressure_dofs_on_inflow_boundary[i],
					 pressure_dofs_on_inflow_boundary[i],
					 -Ke_pre_laplacian
					 (local_pressure_dof,
					  local_pressure_dof));
				  // then add the correct increment
				  system->
				    request_matrix
				    ("Pressure Laplacian Matrix")->
				    add (pressure_dofs_on_inflow_boundary[i],
					 pressure_dofs_on_inflow_boundary[i],
					 diagonal_scaling_increment * Re);
				}
			      else if (es->parameters.get <
				       unsigned int
				       >("pcd_boundary_condition_type") == 3)
				{
				  // first subject the increment we in the matrix construction of this element, so that can add average later
				  system->
				    request_matrix
				    ("Pressure Laplacian Matrix")->
				    add (pressure_dofs_on_inflow_boundary[i],
					 pressure_dofs_on_inflow_boundary[i],
					 -Ke_pre_laplacian
					 (local_pressure_dof,
					  local_pressure_dof));

				  // first subject the increment we in the matrix construction of this element
				  system->
				    request_matrix
				    ("Pressure Convection Diffusion Matrix")->
				    add (pressure_dofs_on_inflow_boundary[i],
					 pressure_dofs_on_inflow_boundary[i],
					 -Ke_pre_convection_diffusion
					 (local_pressure_dof,
					  local_pressure_dof));


				}

			    }
			  else
			    {
			      std::
				cout <<
				"Error local pressure dof not found in Picard.C... EXITING"
				<< std::endl;
			      std::exit (0);
			    }
			}
		    }
		}

	    }
	  else
	    {
	      // estimate error TODO
	    }




	}

    }				// end of element loop









  // *********** APPLY THE AVERAGED PRECONDITIONER BOUNDARY CONDITION ****** //
  //find the average value of
  if (es->parameters.get < unsigned int >("pcd_boundary_condition_type") == 3)
    {
      // sort the dofs
      std::sort (all_global_pressure_dofs_on_inflow_boundary.begin (),
		 all_global_pressure_dofs_on_inflow_boundary.end ());
      // remove duplicate dofs
      std::vector < unsigned int >::iterator it;
      it =
	std::unique (all_global_pressure_dofs_on_inflow_boundary.begin (),
		     all_global_pressure_dofs_on_inflow_boundary.end ());
      all_global_pressure_dofs_on_inflow_boundary.
	resize (std::
		distance (all_global_pressure_dofs_on_inflow_boundary.
			  begin (), it));
      // set the parameter
      double ave_fp_diag =
	sum_fp_diag /
	(double) all_global_pressure_dofs_on_inflow_boundary.size ();

      for (unsigned int i = 0;
	   i < all_global_pressure_dofs_on_inflow_boundary.size (); i++)
	{
	  system->request_matrix ("Pressure Convection Diffusion Matrix")->
	    set (all_global_pressure_dofs_on_inflow_boundary[i],
		 all_global_pressure_dofs_on_inflow_boundary[i], ave_fp_diag);
	  system->request_matrix ("Pressure Laplacian Matrix")->
	    set (all_global_pressure_dofs_on_inflow_boundary[i],
		 all_global_pressure_dofs_on_inflow_boundary[i],
		 Re * ave_fp_diag);
	}
    }




  // put ones on the diagonal of the pressure block of the navier stokes matrix preconditioner
  // preconditioner 9 doesn't use the velocity matrix, takes it from the system matrix
  if (es->parameters.get < unsigned int >("preconditioner_type_3d1d") == 8
      || es->parameters.get < unsigned int >("preconditioner_type_3d1d") ==
      11)
    {
      const unsigned int p_var = system->variable_number ("p");
      std::vector < dof_id_type > p_var_idx;
      system->get_dof_map ().local_variable_indices
	(p_var_idx, system->get_mesh (), p_var);
      // velocity mass matrix
      for (unsigned int i = 0; i < p_var_idx.size (); i++)
	system->request_matrix ("Velocity Matrix")->set (p_var_idx[i],
							 p_var_idx[i], 1.0);
    }










  if (threed)
    std::cout << "end 3D assembly" << std::endl;
  else
    std::cout << "end 2D assembly" << std::endl;

  //std::cout << "tau_sum = " << tau_sum << std::endl;
  std::cout << "average alpha factor = " << alpha_factor_sum /
    count << std::endl;


  return;


}
