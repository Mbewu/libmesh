
#include "navier_stokes_assembler.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

// The matrix assembly function to be called at each time step to
// prepare for the linear solve. This uses a picard type linearisation
// as opposed to the newton type linearisation used in the 
// assemble_stokes function. Well now assemble_stokes is not working in 3D
void NavierStokesAssembler::assemble ()
{
	
	if(!es->parameters.get<bool> ("assemble_residual_only"))
		assemble_efficient ();
	// only do residual formulation for mono
	if(es->parameters.get<bool> ("residual_formulation_0d"))
		this->assemble_residual_rhs();
}

// The matrix assembly function to be called at each time step to
// prepare for the linear solve. This uses a picard type linearisation
// as opposed to the newton type linearisation used in the 
// assemble_stokes function. Well now assemble_stokes is not working in 3D
void NavierStokesAssembler::assemble_efficient ()
{

	std::cout << "Begin 0D assembly... ";

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es->get_mesh();

  // The dimension that we are running
  //const unsigned int dim = mesh.mesh_dimension();	 // unused

 
	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(coupled)
	{
			std::cout << "coupled" << std::endl;
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
			std::cout << "uncoupled" << std::endl;
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns1d");
	}

	// check if it is the first time assembling this
	if(es->parameters.get<bool> ("residual_formulation_0d"))
	{
		// if we have not changed timestep or not changed 3d nonlinear iteration, then we must NOT reassemble the stuff into the forcing vector
		if(!first_assembly && (current_timestep != (int)es->parameters.get<unsigned int>("t_step") || 
				current_nonlinear_iteration_3d != (int)es->parameters.get < unsigned int >("nonlinear_iteration")))
		{
			first_assembly = true;
			current_timestep = (int)es->parameters.get<unsigned int>("t_step");
			current_nonlinear_iteration_3d = (int)es->parameters.get < unsigned int >("nonlinear_iteration");
		}

		// if not coupled, we can assemble everything every time
		if(!coupled)
		{
			first_assembly = true;
			system->get_vector("Forcing Vector BC").close();
			system->get_vector("Forcing Vector BC").zero();
			system->get_vector("Forcing Vector").close();
			system->get_vector("Forcing Vector").zero();
		}
	}

  // Numeric ids corresponding to each variable in the system
  const unsigned int p_var = system->variable_number ("P");
  const unsigned int q_var = system->variable_number ("Q");
  unsigned int v_var = 0;
	if(es->parameters.get<unsigned int>("acinar_model") == 1)
	{
  	v_var = system->variable_number ("V");
	}



	//some parameters

  const double density = es->parameters.get<double>("density");
  double viscosity = es->parameters.get<double>("viscosity");
  double length_scale = es->parameters.get<double>("length_scale");
  double velocity_scale = es->parameters.get<double>("velocity_scale");
  double reynolds_number = es->parameters.get<double>("reynolds_number");
  const double zeta_1 = es->parameters.get<double>("zeta_1");
  const double zeta_2 = es->parameters.get<double>("zeta_2");
  const double zeta_3 = es->parameters.get<double>("zeta_2");
  //const double period = es->parameters.get<double>("period");	// unused
  const double E = es->parameters.get<double>("E");
  const unsigned int unsteady = es->parameters.get<unsigned int>("unsteady");
  const Real dt    = es->parameters.get<Real>("dt");
  const Real time    = es->parameters.get<Real>("time");
	//const double time = system->time;	// unused
	const bool compliance_1d = es->parameters.get<bool>("compliance_1d");
	const bool inertance_1d = es->parameters.get<bool>("inertance_1d");
  const unsigned int resistance_type_1d = es->parameters.get<unsigned int>("resistance_type_1d");
  const bool first_gen_poiseuille = es->parameters.get<bool>("first_gen_poiseuille");
  const bool semi_implicit_1d = es->parameters.get<bool>("semi_implicit_1d");
  const bool newton_1d = es->parameters.get<bool>("newton_1d");
	const unsigned int scale_mono_preconditioner = es->parameters.get<unsigned int>("scale_mono_preconditioner");
	const unsigned int bc_type_1d = es->parameters.get<unsigned int>("bc_type_1d");
	const bool residual_formulation_0d = es->parameters.get<bool>("residual_formulation_0d");

	// mono flow rate penalty param
  double mono_flow_rate_penalty_param = 1.;
	if(es->parameters.get <bool> ("mono_flow_rate_penalty"))
		mono_flow_rate_penalty_param = es->parameters.get <double> ("mono_flow_rate_penalty_param");

	// acinar stuff
	const unsigned int acinar_model = es->parameters.get<unsigned int>("acinar_model");
	const double total_acinar_compliance = es->parameters.get<double>("total_acinar_compliance");
	double total_terminal_area = 0.;
	double pl = 0.;
	double pl_diff = 0.;
	if(acinar_model == 1)
	{
		pl = calculate_pleural_pressure(time);
		double pl_previous = calculate_pleural_pressure(time - dt);
		pl_diff = pl - pl_previous;
		std::cout << "using pleural pressure = " << pl << std::endl;
		total_terminal_area = es->parameters.get<double>("total_terminal_area");
	}

	// hard-coded resistances
	const double resistance_0 = es->parameters.get<double>("resistance_0");
	const double resistance_0_a = es->parameters.get<double>("resistance_0_a");
	const double resistance_0_b = es->parameters.get<double>("resistance_0_b");
	const double resistance_1 = es->parameters.get<double>("resistance_1");
	const double resistance_2 = es->parameters.get<double>("resistance_2");
	const double resistance_3 = es->parameters.get<double>("resistance_3");
	const double resistance_4 = es->parameters.get<double>("resistance_4");

	// if we are doing a reynolds_number_calculation we need to get the viscosity from reynolds_number
	// length scale and velocity scale being 1
	if(es->parameters.get<bool> ("reynolds_number_calculation"))
		viscosity = 1./reynolds_number;

	std::cout << "resistance_type_1d = " << resistance_type_1d << std::endl;

	//van ertbruggen parameters
	std::vector<double> ertbruggen_gamma(9);
	ertbruggen_gamma[0] = 0.162;
	ertbruggen_gamma[1] = 0.239;
	ertbruggen_gamma[2] = 0.244;
	ertbruggen_gamma[3] = 0.295;
	ertbruggen_gamma[4] = 0.175;
	ertbruggen_gamma[5] = 0.303;
	ertbruggen_gamma[6] = 0.356;
	ertbruggen_gamma[7] = 0.566;
	ertbruggen_gamma[8] = 0.327;	// generation > 7


	// BCs
	const double in_pressure = es->parameters.get<double> ("in_pressure_1d")* es->parameters.get<double> ("time_scaling");
	const double out_pressure = es->parameters.get<double> ("out_pressure_1d")* es->parameters.get<double> ("time_scaling");
	const double out_pressure_diff = out_pressure - es->parameters.get<double> ("out_pressure_1d")* es->parameters.get<double> ("previous_time_scaling");

	// set if using a preconditioner matrix or not
	if((es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 6 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 7 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 8 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 9 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 10 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 11 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 12 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 13) 
		&& coupled)
		preconditioner_assemble = true;
	else
		preconditioner_assemble = false;

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap & dof_map = system->get_dof_map();

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_p;
  std::vector<dof_id_type> dof_indices_q;
  std::vector<dof_id_type> dof_indices_parent_p;
  std::vector<dof_id_type> dof_indices_parent_q;
  std::vector<dof_id_type> dof_indices_daughter_1_p;
  std::vector<std::vector<dof_id_type> > dof_indices_siblings_q;
  std::vector<std::vector<dof_id_type> > dof_indices_siblings_p;

  std::vector<dof_id_type> dof_indices_acinar;
  std::vector<dof_id_type> dof_indices_v;

	// resistance scaling value
	double mono_preconditioner_resistance_scaling = 0.;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the \p active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	std::vector<unsigned int> zero_d_dofs;

								
  for ( ; el != end_el; ++el)
  {
		const Elem* elem = *el;
		unsigned int subdomain_id = elem->subdomain_id();
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), subdomain_id) != subdomains_1d.end())
		{

			//element data object starts numbering from 0
			//and the values referenced in it also do so need to take this into account

			const int current_el_idx = elem->id();
			//unsigned int current_1d_el_idx = current_el_idx -	n_initial_3d_elem;
			unsigned int current_1d_el_idx = elem_to_airway[current_el_idx];// -	n_initial_3d_elem;
			bool has_parent = airway_data[current_1d_el_idx].has_parent();
			int parent_el_idx = 0;
			if(has_parent)
			{
				//parent_el_idx = (int)airway_data[current_1d_el_idx].get_parent() + n_initial_3d_elem;
				parent_el_idx = airway_data[(int)airway_data[current_1d_el_idx].get_parent()].get_local_elem_number();
			}
			std::vector<unsigned int> daughter_el_ids = airway_data[current_1d_el_idx].get_daughters();
			int daughter_1_el_idx = -1;
			if(daughter_el_ids.size() > 0)
			{
				//daughter_1_el_idx = daughter_el_ids[0] + n_initial_3d_elem;
				daughter_1_el_idx = airway_data[daughter_el_ids[0]].get_local_elem_number();// + n_initial_3d_elem;
			}
			
			bool is_daughter_1 = airway_data[current_1d_el_idx].get_is_daughter_1();	//this is a bool duh!
			std::vector<unsigned int> sibling_el_ids = airway_data[current_1d_el_idx].get_siblings();
			for(unsigned int i=0; i<sibling_el_ids.size(); i++)
			{
				//sibling_el_ids[i] += n_initial_3d_elem;
				sibling_el_ids[i] = airway_data[sibling_el_ids[i]].get_local_elem_number();//n_initial_3d_elem;
				//std::cout << "sibling_el_ids = " << sibling_el_ids[i] << std::endl;
			}

			int primary_sibling_el_idx = -1;
			if(!is_daughter_1)
			{
				//primary_sibling_el_idx = airway_data[current_1d_el_idx].get_primary_sibling() + n_initial_3d_elem;
				primary_sibling_el_idx = airway_data[airway_data[current_1d_el_idx].get_primary_sibling()].get_local_elem_number();// + n_initial_3d_elem;
			}
			
			const double l = airway_data[current_1d_el_idx].get_length();	//nondimensionalised length
			const double r = airway_data[current_1d_el_idx].get_radius(); //nondimensionalised radius
			int generation = airway_data[current_1d_el_idx].get_generation();

			bool has_sibling = false;
			if(sibling_el_ids.size() > 0)
				has_sibling = true;

			//very dirty hack
			if(!has_parent)
				parent_el_idx = current_el_idx;

			//if(primary_sibling_el_idx  < n_initial_3d_elem)
			if(primary_sibling_el_idx  < 0)
				primary_sibling_el_idx  = current_el_idx;

			// don't need cause we check with size of sibling_el_ids
			//if(daughter_1_el_idx  < n_initial_3d_elem)
			if(daughter_1_el_idx  < 0)
				daughter_1_el_idx  = current_el_idx;

			
			const Elem* parent_elem = mesh.elem(parent_el_idx);
			const Elem* daughter_1_elem = mesh.elem(daughter_1_el_idx);

			std::vector<const Elem*> sibling_elems;
			for(unsigned int i=0; i<sibling_el_ids.size(); i++)
				sibling_elems.push_back(mesh.elem(sibling_el_ids[i]));

			//const Elem* primary_sibling_elem = mesh.elem(primary_sibling_el_idx); // unused
			//const Elem* sibling_elem = mesh.elem(sibling_el_idx);

			/*
			if(parent_el_idx >= 0)
				parent_elem = mesh.elem(parent_el_idx);

			if(daughter_1_el_idx >= 0)
				daughter_1_elem = mesh.elem(daughter_1_el_idx);

			if(daughter_2_el_idx >= 0)
				daughter_2_elem = mesh.elem(daughter_2_el_idx);
			*/

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_q, q_var);
      dof_map.dof_indices (parent_elem, dof_indices_parent_p, p_var);
      dof_map.dof_indices (parent_elem, dof_indices_parent_q, q_var);
      dof_map.dof_indices (daughter_1_elem, dof_indices_daughter_1_p, p_var);

			dof_indices_siblings_q.resize(sibling_el_ids.size());
			dof_indices_siblings_p.resize(sibling_el_ids.size());

			for(unsigned int i=0; i< sibling_el_ids.size(); i++)
			{
		    dof_map.dof_indices (sibling_elems[i], dof_indices_siblings_q[i], q_var);
		    dof_map.dof_indices (sibling_elems[i], dof_indices_siblings_p[i], p_var);
			}


			for(unsigned int i=0; i<dof_indices_p.size(); i++)
				zero_d_dofs.push_back(dof_indices_p[i]);

			for(unsigned int i=0; i<dof_indices_q.size(); i++)
				zero_d_dofs.push_back(dof_indices_q[i]);

			// get the acinar dofs and other acinar stuff
			unsigned int acinar_v_dof_idx = 0;
			double local_acinar_compliance = 0.;
			if(acinar_model == 1)
			{
				// first check if it has an acinar i.e. terminal
				int acinar_elem_number = airway_data[current_1d_el_idx].get_acinar_elem();
				if(acinar_elem_number > -1)
				{
					// now get the correct element
					const Elem* acinar_elem = mesh.elem(acinar_elem_number);
//					std::cout << "v_var = " << v_var << std::endl;
//					std::cout << "acinar_elem_number = " << acinar_elem_number << std::endl;
      		dof_map.dof_indices (acinar_elem, dof_indices_acinar);
//					std::cout << "hey babe" << std::endl;
//					std::cout << "dof_indices.size() = " << dof_indices_acinar.size() << std::endl;
//					std::cout << "hey babe" << std::endl;
      		dof_map.dof_indices (acinar_elem, dof_indices_v, v_var);
//					std::cout << "dof_indices_v.size() = " << dof_indices_v.size() << std::endl;
					acinar_v_dof_idx = dof_indices_v[0];
					
					local_acinar_compliance = airway_data[current_1d_el_idx].get_local_acinar_compliance();
//					std::cout << "local_acinar_compliance = " << local_acinar_compliance << std::endl;
				}

			}
			







      //const unsigned int n_dofs   = 4;	// unused
      //const unsigned int n_p_dofs = 2;	// unused

			// some parameters

			// use the current solution or the previous time step's solution
			// if no compliance, flow rate is the same at the beginning and end of the pipe
			double average_flow = 0.;
			if(semi_implicit_1d)
			{
				if(compliance_1d)
					average_flow = ((*system->old_local_solution)(dof_indices_q[0]) + (*system->old_local_solution)(dof_indices_q[1]))/2.0;
				else
					average_flow = (*system->old_local_solution)(dof_indices_q[0]);
			}
			else
			{
				if(compliance_1d)
					average_flow = (system->current_solution(dof_indices_q[0]) + system->current_solution(dof_indices_q[1]))/2.0;
				else
					average_flow = system->current_solution(dof_indices_q[0]);
			}


			// construct dimensionalised parameters
			double R_poi = airway_data[current_1d_el_idx].get_poiseuille_resistance();
			double R_deriv = 0.;

			double C = 0.;
			if(compliance_1d)
			{
				C = airway_data[current_1d_el_idx].get_airway_compliance();
			}

			double I = 0.;
			if(inertance_1d)
			{
				C = airway_data[current_1d_el_idx].get_inertance();

			}


			double R = 0.;
			

/*
			std::cout << "length_scale = " << length_scale << std::endl;
			std::cout << "velocity_scale = " << velocity_scale << std::endl;
			std::cout << "R = " << R << std::endl;
			std::cout << "viscosity = " << viscosity << std::endl;
			std::cout << "density = " << density << std::endl;
			std::cout << "l = " << l << std::endl;
			std::cout << "r = " << r << std::endl;
			//std::cout << "r^4 = " << pow(r,4.0) << std::endl;
*/


			//1d resistance types
			if(resistance_type_1d == 0)
			{
				//poiseuille resistance
				R=R_poi;
				airway_data[current_1d_el_idx].set_poiseuille(true);
			}
			else if(resistance_type_1d == 1)
			{
				//pedley resistance

				double reynolds_number_0d = 2*density*velocity_scale*length_scale*fabs(average_flow)/(viscosity*M_PI*r);
				double pedley_factor = 0.327 * sqrt(reynolds_number_0d*2.*r/l);

				if(generation < 5)
				{
					//std::cout << "generation = " << generation << std::endl;
					//std::cout << "reynolds_number_0d = " << reynolds_number_0d << std::endl;
					//std::cout << "pedley factor = " << 0.327 * sqrt(reynolds_number_0d*2.*r/l) << std::endl;
					//std::cout << "length = " << l << std::endl;
					//std::cout << "radius = " << r << std::endl;
				}

				if(generation == 0 && first_gen_poiseuille)
				{
					std::cout << "using first gen poi" << std::endl;
					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
				}
				else if(pedley_factor > 1.0)
				{

					R = pedley_factor * R_poi;
					if(newton_1d)
					{
						R_deriv = 0.327 * 2*density/(viscosity*M_PI*l) / sqrt(reynolds_number_0d*2.*r/l) * R_poi;
					}

					airway_data[current_1d_el_idx].set_poiseuille(false);
				}
				else
				{

					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
					
				}
			}
			else if(resistance_type_1d == 2)
			{


				double reynolds_number_0d = 2*density*velocity_scale*length_scale*fabs(average_flow)/(viscosity*M_PI*r);

				double gamma = 0;
				if(generation < 0)		//alveolar
					gamma = ertbruggen_gamma[8];
				else if(generation < 8)
					gamma = ertbruggen_gamma[generation];
				else
					gamma = ertbruggen_gamma[8];

				// van ertbruggen resistance
				double ertbruggen_factor = gamma * sqrt(reynolds_number_0d*2.*r/l);

				if(generation == 0 && first_gen_poiseuille)
				{
					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
				}
				else if(ertbruggen_factor > 1.0)
				{

					R = ertbruggen_factor * R_poi;
					if(newton_1d)
						R_deriv =  gamma * 2*density/(viscosity*M_PI*l) / sqrt(reynolds_number_0d*2.*r/l) * R_poi;

					airway_data[current_1d_el_idx].set_poiseuille(false);
				}
				else
				{
					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
				}
			}
			else if(resistance_type_1d == 3)
			{
	

				double reynolds_number_0d = 2*density*velocity_scale*length_scale*fabs(average_flow)/(viscosity*M_PI*r);

				double constant = 0;
				//lobar bronchi defined as generation less than equal 2, but doesn't quite get it, oops
				if(generation > 2 || generation < 0)
					constant = 1.0;
				else
					constant = 3.4;

				// van ertbruggen resistance
				if(es->parameters.get<unsigned int> ("t_step") != 1)
				{
					R = (constant + 2.1e-3 * reynolds_number_0d) * R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(false);
				}
				else
				{
					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
				}
			}
			else if(resistance_type_1d == 4)
			{
				if(generation == 0)
					R = resistance_0;
				else if(generation == 1)
					R = resistance_1;
				else if(generation == 2)
					R = resistance_2;
				else if(generation == 3)
					R = resistance_3;
				else if(generation == 4)
					R = resistance_4;
				else
				{
					std::cout << "hard-coded resistances not implmeneted for meshes with more than 4 generations." << std::endl;
					std::cout << "Exiting..." << std::endl;
					std::exit(0);
				}

				std::cout << "gen " << generation << std::endl;
				std::cout << "resistance = " << R << std::endl;
			}
			else if(resistance_type_1d == 5)
			{
				if(generation == 0 && current_1d_el_idx == 0)
					R = resistance_0_a;
				else if(generation == 0 && current_1d_el_idx == 1)
					R = resistance_0_b;
				else
				{
					std::cout << "hard-coded asymmetric resistances not implmeneted for meshes with more than 1 generations." << std::endl;
					std::cout << "Exiting..." << std::endl;
					std::exit(0);
				}

				std::cout << "gen " << generation << std::endl;
				std::cout << "resistance = " << R << std::endl;
			}
			else
			{
					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
			}

			// set the resistance so that can be output
			airway_data[current_1d_el_idx].set_resistance(R);


			// set the resistance to be used by the monolithic preconditioner
			if(generation == 0)
			{
				if(R > mono_preconditioner_resistance_scaling)
				{
					mono_preconditioner_resistance_scaling = R;
				}
			}

			double old_p0 = system->old_solution(dof_indices_p[0]);
			//double old_p1 = system->old_solution(dof_indices_p[1]);	// unused
			//double old_q0 = system->old_solution(dof_indices_q[0]);	// unused
			double old_q1 = system->old_solution(dof_indices_q[1]);

			double old_p0_diff = system->older_solution(dof_indices_p[0]) - system->old_solution(dof_indices_p[0]);
			double old_q1_diff = system->older_solution(dof_indices_q[1]) - system->old_solution(dof_indices_q[1]);

			//std::cout << "R = " << R << std::endl;
			//std::cout << "R_deriv = " << R_deriv << std::endl;
		
		
			//bool compute_using_monomials = false;	// unused
			{
				// first we decide which equations go where based on the boundary conditions that exist

				// the first equation is the compliance equation (mass conservation)
				// the second equation is the resistance equation (flow = grad(p))
				// the third equation is always an equation for the inflow boundary
				// the fourth equation is always an equation for the outflow boundary 
				std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,0);

				unsigned int eqn_1_dof = 0;
				unsigned int eqn_2_dof = 0;
				unsigned int eqn_3_dof = 0;		//put influx condition in p[0], what why?
				unsigned int eqn_4_dof = 0;

				//default
				if(is_daughter_1)
				{
					eqn_1_dof = dof_indices_q[1];
					eqn_2_dof = dof_indices_p[0];
					eqn_3_dof = dof_indices_q[0];
					eqn_4_dof = dof_indices_p[1];
				}
				else
				{
					eqn_1_dof = dof_indices_q[0];
					eqn_2_dof = dof_indices_q[1];
					eqn_3_dof = dof_indices_p[0];
					eqn_4_dof = dof_indices_p[1];
				}

				//inflow bc - must be a major/daughter_1 branch
				// takes care of all inflow bc and inflow and outflow bc
				if(parent_el_idx == current_el_idx)
				{
					//flow prescribed at inflow
					if(bc_type_1d == 0)
					{
						// these actually remain the same
						// eqn_3_dof = dof_indices_q[0];		//put influx condition in eqn_3, no change
						// if we also have the pressure prescribed at the outflow then eqn_4 is also unchanged, still p_out

				
					}
					//pressure prescribed at inflow
					else if(bc_type_1d == 1)
					{
						eqn_3_dof = dof_indices_p[0];		//put pressure inflow condition in p[0]
						eqn_1_dof = dof_indices_q[0];		// need to put eqn_1 in q_in row because the pressure bc has taken it's place
						eqn_2_dof = dof_indices_q[1];		// need to put eqn_2 in q_out row because the q_in has taken it
						// if we also have the pressure prescribed at the outflow then eqn_4 is also unchanged, still p_out
					}
				}

				//only outflow bc (we only consider pressure outflow boundary conditions which leave the equations unchanged)				
						
				// ********** PUTTING THE EQUATIONS IN THE MATRIX AND RHS ********* //

				// **** equation 1 - compliance eqn
	    	add_to_matrices (system,eqn_1_dof,dof_indices_q[1],1.0);
	    	add_to_matrices (system,eqn_1_dof,dof_indices_q[0],-1.0);

				if(unsteady)
				{
	    		add_to_matrices (system,eqn_1_dof,dof_indices_p[0],C/dt);
					if(!residual_formulation_0d)
	    			system->rhs->add (eqn_1_dof,C*old_p0/dt);
					else 
					{
						if(first_assembly)
						{
	    				//system->get_vector("Forcing Vector BC").add (eqn_1_dof,C*old_p0_diff/dt);
							system->get_vector("Forcing Vector BC").add (eqn_1_dof,C*old_p0/dt);
							system->get_vector("Forcing Vector").add (eqn_1_dof,C*old_p0/dt);
						}
						else
	    				system->get_vector("Forcing Vector BC").add (eqn_1_dof,0);
					}
				}

				// **** equation 2 - resistance eqn
				if(!newton_1d)
		    	add_to_matrices (system,eqn_2_dof,dof_indices_q[1],R);
				else
				{
		    	add_to_matrices (system,eqn_2_dof,dof_indices_q[1],R);
		    	add_to_matrices (system,eqn_2_dof,dof_indices_q[1],R_deriv*fabs(average_flow));
				}

	    	add_to_matrices (system,eqn_2_dof,dof_indices_p[1],1.0);
	    	add_to_matrices (system,eqn_2_dof,dof_indices_p[0],-1.0);

				if(unsteady)
				{
	    		add_to_matrices (system,eqn_2_dof,dof_indices_q[1],I/dt);
					if(!residual_formulation_0d)
	    			system->rhs->add (eqn_2_dof,I*old_q1/dt);
					else 
					{
						if(first_assembly)
						{
	    				//system->get_vector("Forcing Vector BC").add (eqn_2_dof,I*old_q1_diff/dt);
							system->get_vector("Forcing Vector BC").add (eqn_2_dof,I*old_q1/dt);
							system->get_vector("Forcing Vector").add (eqn_2_dof,I*old_q1/dt);
						}
						else // do nothing
	    				system->get_vector("Forcing Vector BC").add (eqn_2_dof,0);
					}

				}

				if(newton_1d && !residual_formulation_0d)
    			system->rhs->add (eqn_2_dof,R_deriv*fabs(average_flow)*average_flow);


				// **** equation 3 - BC on inflow
				if(parent_el_idx == current_el_idx)
				{
					//std::cout << "parent_el_idx = " << parent_el_idx << std::endl;
					//std::cout << "current_el_idx = " << current_el_idx << std::endl;
					// flow
					if(bc_type_1d == 0)
					{
						// now there is no zero on the 3rd equation
						// equation 3 - inflow bc multiplied by dt so is like the 3d eqn - not anymore it isn't
	
						// daughter 1 handles the flux conservation to sibling and parent (1 equation)
						if(is_daughter_1)
						{
							
							add_to_matrices (system,eqn_3_dof,dof_indices_q[0],-mono_flow_rate_penalty_param);

							//std::cout << "flux dof in 0d = " << eqn_3_dof << std::endl;
							if(has_sibling)
							{
								//std::cout << "num siblings = " << sibling_el_ids.size() << std::endl;
								//std::cout << "sibling id daughter1 = " << sibling_el_ids[0] << std::endl;
								//std::cout << "current_1d_el_idx  = " << current_1d_el_idx  << std::endl;
								for(unsigned int i=0; i<sibling_el_ids.size(); i++)
									add_to_matrices (system,eqn_3_dof,dof_indices_siblings_q[i][0],-mono_flow_rate_penalty_param);
							}
			
							// boundary condition should be on side 0
							std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,0);
							if(boundary_ids.size() == 0)
								std::cout << "error, boundary does not have boundary id as expected, el idx = " << current_el_idx << std::endl;

							/*
							std::cout << "hi" << std::endl;
							std::cout << "boundary_ids.size() = " << boundary_ids.size() << std::endl;
							for(unsigned int i=0; i<boundary_ids.size(); i++)
								std::cout << boundary_ids[i] << std::endl;
							std::cout << "flux_values.size() = " << flux_values.size() << std::endl;
							for(unsigned int i=0; i<flux_values.size(); i++)
								std::cout << flux_values[i] << std::endl;
							std::cout << "flux_values[" << boundary_ids[0] <<"] = " << flux_values[boundary_ids[0]] << std::endl;
							*/

							if(!coupled)		//if coupled then the rest is put in the matrix by the 3D assembler
							{
								// if there is a bifurcation at the inflow then daughter_1 can add it
								if(is_daughter_1)
								{
									//std::cout << "double oopsie" << std::endl;
									if(!residual_formulation_0d)
										system->rhs->add (eqn_3_dof,-flux_values[boundary_ids[0]]);
									else 
									{
										if(first_assembly)
										{
											//system->get_vector("Forcing Vector BC").add (eqn_2_dof,I*old_q1_diff/dt);
											system->get_vector("Forcing Vector BC").add (eqn_3_dof,-flux_values[boundary_ids[0]]);
											system->get_vector("Forcing Vector").add (eqn_3_dof,-flux_values[boundary_ids[0]]);
										}
										else // do nothing
											system->get_vector("Forcing Vector BC").add (eqn_3_dof,0);
									}

								}
							}
						}
						else	// if we are a daughter 2 then we have a sibling and need to apply the pressure continuity condition
						{
							add_to_matrices (system,eqn_3_dof,dof_indices_p[0],1.0);
							//std::cout << "num siblings = " << sibling_el_ids.size() << std::endl;
							//std::cout << "sibling id daughter2 = " << sibling_el_ids[0] << std::endl;
//							std::cout << "actual val = " << (int)airway_data[current_1d_el_idx][3] + n_initial_3d_elem << std::endl;
							//std::cout << "n_initial_3d_elem = " << n_initial_3d_elem << std::endl;
							//std::cout << "is_daughter_1 = " << is_daughter_1 << std::endl;
							//std::cout << "current_1d_el_idx  = " << current_1d_el_idx  << std::endl;
							//if(!(bool)is_daughter_1)
							//	std::cout << "cool" << std::endl;
							//else
							//	std::cout << "couldn't get in" << std::endl;

							for(unsigned int i=0; i<sibling_el_ids.size(); i++)
								add_to_matrices (system,eqn_3_dof,dof_indices_siblings_p[i][0],-1.0);
						}
					}
					// pressure-pressure type boundary condition
					else if(bc_type_1d == 1)
					{
		    		add_to_matrices (system,eqn_3_dof,dof_indices_p[0],1.0);
						if(!residual_formulation_0d)
							system->rhs->add (eqn_3_dof,in_pressure);
						else 
						{
							if(first_assembly)
							{
								//system->get_vector("Forcing Vector BC").add (eqn_2_dof,I*old_q1_diff/dt);
								system->get_vector("Forcing Vector BC").add (eqn_3_dof,in_pressure);
								system->get_vector("Forcing Vector").add (eqn_3_dof,in_pressure);
							}
							else // do nothing
								system->get_vector("Forcing Vector BC").add (eqn_3_dof,0);
						}
					}

				}
				else
				{
					if(is_daughter_1)
					{
						// conservation of flux
						add_to_matrices (system,eqn_3_dof,dof_indices_parent_q[1],1.0);
						add_to_matrices (system,eqn_3_dof,dof_indices_q[0],-1.0);
						for(unsigned int i=0; i<sibling_el_ids.size(); i++)
							add_to_matrices (system,eqn_3_dof,dof_indices_siblings_q[i][0],-1.0);
					}
					else
					{
						// parent pressure cont
				  	add_to_matrices (system,eqn_3_dof,dof_indices_p[0],1.0);
				  	add_to_matrices (system,eqn_3_dof,dof_indices_parent_p[1],-1.0);
					}

				}


				// **** equation 4 - BC on outflow
				if(daughter_1_el_idx  == current_el_idx)
				{

					if(acinar_model == 1)
					{
//						std::cout << "gru1" << std::endl;
						// acinar compliance equation
			  		add_to_matrices (system,eqn_4_dof,dof_indices_p[1],1.0);
//						std::cout << "gru1" << std::endl;
			  		add_to_matrices (system,eqn_4_dof,acinar_v_dof_idx,-1.0/local_acinar_compliance);
//						std::cout << "gru1" << std::endl;
						if(!residual_formulation_0d)
			  			system->rhs->add (eqn_4_dof,pl);
						else 
						{
							if(first_assembly)
							{
			  				//system->get_vector("Forcing Vector BC").add (eqn_4_dof,pl_diff);
								system->get_vector("Forcing Vector BC").add (eqn_4_dof,pl);
								system->get_vector("Forcing Vector").add (eqn_4_dof,pl);
							}
							else // do nothing
			  				system->get_vector("Forcing Vector BC").add (eqn_4_dof,0);
						}
						
//						std::cout << "gru1" << std::endl;
						
					}
					else
					{
						//apply pressure
			  		add_to_matrices (system,eqn_4_dof,dof_indices_p[1],1.0);
						if(!residual_formulation_0d)
			  			system->rhs->add (eqn_4_dof,out_pressure);
						else 
						{
							if(first_assembly)
							{
			  				//system->get_vector("Forcing Vector BC").add (eqn_4_dof,out_pressure_diff);
								system->get_vector("Forcing Vector BC").add (eqn_4_dof,out_pressure);
								system->get_vector("Forcing Vector").add (eqn_4_dof,out_pressure);
							}
							else // do nothing
			  				system->get_vector("Forcing Vector BC").add (eqn_4_dof,0);
						}
					}

				}
				else
				{
					//pressure continuity to daughter_1
	    		add_to_matrices (system,eqn_4_dof,dof_indices_p[1],1.0);
	    		add_to_matrices
							(system,eqn_4_dof,dof_indices_daughter_1_p[0],-1.0);
				}
			
			}

    }// end of subdomains_1d conditional
		else if(std::find(subdomains_acinar.begin(), subdomains_acinar.end(), subdomain_id) != subdomains_acinar.end())
		{
//						std::cout << "gru2" << std::endl;
			// need to find the element that this acinar is assoc with
			const int current_el_idx = elem->id();
			const unsigned int acinar_id = elem_to_acinar[current_el_idx];
			//const unsigned int parent_el_idx = acinar_to_airway[acinar_id] + n_initial_3d_elem;
			const unsigned int parent_el_idx = airway_data[acinar_to_airway[acinar_id]].get_local_elem_number();// + n_initial_3d_elem;
			const Elem* parent_elem = mesh.elem(parent_el_idx);

//						std::cout << "gru2" << std::endl;
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (parent_elem, dof_indices_q, q_var);

//						std::cout << "gru2" << std::endl;
			// get old volume value
			double v_old = (*system->old_local_solution)(dof_indices_v[0]);
			double v_old_diff = v_old - (*system->older_local_solution)(dof_indices_v[0]);

			//std::cout << "v_old = " << v_old << std::endl;

//						std::cout << "gru2" << std::endl;
			// volume-flow equation
  		add_to_matrices (system,dof_indices_v[0],dof_indices_v[0],1.0/dt);
//						std::cout << "gru2" << std::endl;
  		add_to_matrices (system,dof_indices_v[0],dof_indices_q[1],-1.0);
//						std::cout << "gru2" << std::endl;
			if(!residual_formulation_0d)
  			system->rhs->add (dof_indices_v[0],v_old/dt);
			else 
			{
				if(first_assembly)
				{
  				//system->get_vector("Forcing Vector BC").add (dof_indices_v[0],v_old_diff/dt);
					system->get_vector("Forcing Vector BC").add (dof_indices_v[0],v_old/dt);
					system->get_vector("Forcing Vector").add (dof_indices_v[0],v_old/dt);
				}
				else // do nothing
  				system->get_vector("Forcing Vector BC").add (dof_indices_v[0],0.);
			}
//						std::cout << "gru2" << std::endl;
		}
			
	}// end of element loop





	es->parameters.set<double> ("mono_preconditioner_resistance_scaling") = mono_preconditioner_resistance_scaling;
	//std::cout << "mono_preconditioner_resistance_scaling set to be = " << es->parameters.get<double> ("mono_preconditioner_resistance_scaling") << std::endl;
	/*
	if(scale_mono_preconditioner == 0)
		std::cout << "scaling not used" << std::endl;
	*/

	//std::cout << "length_scale = " << length_scale << std::endl;
	//std::cout << "velocity_scale = " << velocity_scale << std::endl;


  // That's it.
	//system->matrix->close();

	/*
	std::cout << "0d matrix" << std::endl;
	for(unsigned int i=0; i<zero_d_dofs.size();i++)
	{
		for(unsigned int j=0;j<zero_d_dofs.size();j++)
		{
			std::cout << " " << (*(system->matrix))(zero_d_dofs[i],zero_d_dofs[j]);
		}
		std::cout << std::endl;
	}
	*/
		
	//system->matrix->print();
	//system->rhs->print();

	//system->matrix->close();
	//system->matrix->print();
	//std::cout << "true rhs" << std::endl;
	//system->rhs->print();

	std::cout << "End 0D assembly" << std::endl;
  return;
}

// This adds values to the system matrix and/or preconditioner
void NavierStokesAssembler::add_to_matrices (TransientLinearImplicitSystem * system, unsigned int& row_number, unsigned int& col_number, double value)
{
	
	system->matrix->add (row_number,col_number,value);
	if(preconditioner_assemble)
		system->request_matrix("Preconditioner")->add (row_number,col_number,value);
	
}


// assmebl the residual rhs
void NavierStokesAssembler::assemble_residual_rhs ()
{

    std::cout << "\nBeginning 0D assembly of residual rhs... " << std::endl;

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es->get_mesh();

  // The dimension that we are running
  //const unsigned int dim = mesh.mesh_dimension();	 // unused

 
	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(coupled)
	{
			std::cout << "coupled" << std::endl;
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
			std::cout << "uncoupled" << std::endl;
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns1d");
	}

	if(!coupled)
	{
		// zero this first as it may be run multiple times
		system->get_vector("Residual LHS").close();
		system->get_vector("Residual LHS").zero();
	}

  // Numeric ids corresponding to each variable in the system
  const unsigned int p_var = system->variable_number ("P");
  const unsigned int q_var = system->variable_number ("Q");
  unsigned int v_var = 0;
	if(es->parameters.get<unsigned int>("acinar_model") == 1)
	{
  	v_var = system->variable_number ("V");
	}

	// residual lhs will have been zeroed in the 

	//some parameters

  const double density = es->parameters.get<double>("density");
  double viscosity = es->parameters.get<double>("viscosity");
  double length_scale = es->parameters.get<double>("length_scale");
  double velocity_scale = es->parameters.get<double>("velocity_scale");
  double reynolds_number = es->parameters.get<double>("reynolds_number");
  const double zeta_1 = es->parameters.get<double>("zeta_1");
  const double zeta_2 = es->parameters.get<double>("zeta_2");
  const double zeta_3 = es->parameters.get<double>("zeta_2");
  //const double period = es->parameters.get<double>("period");	// unused
  const double E = es->parameters.get<double>("E");
  const unsigned int unsteady = es->parameters.get<unsigned int>("unsteady");
  const Real dt    = es->parameters.get<Real>("dt");
  const Real time    = es->parameters.get<Real>("time");
	//const double time = system->time;	// unused
	const bool compliance_1d = es->parameters.get<bool>("compliance_1d");
	const bool inertance_1d = es->parameters.get<bool>("inertance_1d");
  const unsigned int resistance_type_1d = es->parameters.get<unsigned int>("resistance_type_1d");
  const bool first_gen_poiseuille = es->parameters.get<bool>("first_gen_poiseuille");
  const bool semi_implicit_1d = es->parameters.get<bool>("semi_implicit_1d");
  const bool newton_1d = es->parameters.get<bool>("newton_1d");
	const unsigned int scale_mono_preconditioner = es->parameters.get<unsigned int>("scale_mono_preconditioner");
	const unsigned int bc_type_1d = es->parameters.get<unsigned int>("bc_type_1d");

	// mono flow rate penalty param
  double mono_flow_rate_penalty_param = 1.;
	if(es->parameters.get <bool> ("mono_flow_rate_penalty"))
		mono_flow_rate_penalty_param = es->parameters.get <double> ("mono_flow_rate_penalty_param");

	// acinar stuff
	const unsigned int acinar_model = es->parameters.get<unsigned int>("acinar_model");
	const double total_acinar_compliance = es->parameters.get<double>("total_acinar_compliance");
	double total_terminal_area = 0.;
	double pl = 0.;
	if(acinar_model == 1)
	{
		pl = calculate_pleural_pressure(time);
		std::cout << "using pleural pressure = " << pl << std::endl;
		total_terminal_area = es->parameters.get<double>("total_terminal_area");
	}

	// hard-coded resistances
	const double resistance_0 = es->parameters.get<double>("resistance_0");
	const double resistance_0_a = es->parameters.get<double>("resistance_0_a");
	const double resistance_0_b = es->parameters.get<double>("resistance_0_b");
	const double resistance_1 = es->parameters.get<double>("resistance_1");
	const double resistance_2 = es->parameters.get<double>("resistance_2");
	const double resistance_3 = es->parameters.get<double>("resistance_3");
	const double resistance_4 = es->parameters.get<double>("resistance_4");

	// if we are doing a reynolds_number_calculation we need to get the viscosity from reynolds_number
	// length scale and velocity scale being 1
	if(es->parameters.get<bool> ("reynolds_number_calculation"))
		viscosity = 1./reynolds_number;

	std::cout << "resistance_type_1d = " << resistance_type_1d << std::endl;

	//van ertbruggen parameters
	std::vector<double> ertbruggen_gamma(9);
	ertbruggen_gamma[0] = 0.162;
	ertbruggen_gamma[1] = 0.239;
	ertbruggen_gamma[2] = 0.244;
	ertbruggen_gamma[3] = 0.295;
	ertbruggen_gamma[4] = 0.175;
	ertbruggen_gamma[5] = 0.303;
	ertbruggen_gamma[6] = 0.356;
	ertbruggen_gamma[7] = 0.566;
	ertbruggen_gamma[8] = 0.327;	// generation > 7


	// BCs
	const double in_pressure = es->parameters.get<double> ("in_pressure_1d")* es->parameters.get<double> ("time_scaling");
	const double out_pressure = es->parameters.get<double> ("out_pressure_1d")* es->parameters.get<double> ("time_scaling");

	// set if using a preconditioner matrix or not
	if((es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 6 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 7 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 8 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 9 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 10 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 11 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 12 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 13) 
		&& coupled)
		preconditioner_assemble = true;
	else
		preconditioner_assemble = false;

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap & dof_map = system->get_dof_map();

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_p;
  std::vector<dof_id_type> dof_indices_q;
  std::vector<dof_id_type> dof_indices_parent_p;
  std::vector<dof_id_type> dof_indices_parent_q;
  std::vector<dof_id_type> dof_indices_daughter_1_p;
  std::vector<std::vector<dof_id_type> > dof_indices_siblings_q;
  std::vector<std::vector<dof_id_type> > dof_indices_siblings_p;

  std::vector<dof_id_type> dof_indices_acinar;
  std::vector<dof_id_type> dof_indices_v;

	// resistance scaling value
	double mono_preconditioner_resistance_scaling = 0.;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the \p active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	std::vector<unsigned int> zero_d_dofs;

								
  for ( ; el != end_el; ++el)
  {
		const Elem* elem = *el;
		unsigned int subdomain_id = elem->subdomain_id();
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), subdomain_id) != subdomains_1d.end())
		{

			//element data object starts numbering from 0
			//and the values referenced in it also do so need to take this into account

			const int current_el_idx = elem->id();
			//unsigned int current_1d_el_idx = current_el_idx -	n_initial_3d_elem;
			unsigned int current_1d_el_idx = elem_to_airway[current_el_idx];// -	n_initial_3d_elem;
			bool has_parent = airway_data[current_1d_el_idx].has_parent();
			int parent_el_idx = 0;
			if(has_parent)
			{
				//parent_el_idx = (int)airway_data[current_1d_el_idx].get_parent() + n_initial_3d_elem;
				parent_el_idx = airway_data[(int)airway_data[current_1d_el_idx].get_parent()].get_local_elem_number();
			}
			std::vector<unsigned int> daughter_el_ids = airway_data[current_1d_el_idx].get_daughters();
			int daughter_1_el_idx = -1;
			if(daughter_el_ids.size() > 0)
			{
				//daughter_1_el_idx = daughter_el_ids[0] + n_initial_3d_elem;
				daughter_1_el_idx = airway_data[daughter_el_ids[0]].get_local_elem_number();// + n_initial_3d_elem;
			}
			
			bool is_daughter_1 = airway_data[current_1d_el_idx].get_is_daughter_1();	//this is a bool duh!
			std::vector<unsigned int> sibling_el_ids = airway_data[current_1d_el_idx].get_siblings();
			for(unsigned int i=0; i<sibling_el_ids.size(); i++)
			{
				//sibling_el_ids[i] += n_initial_3d_elem;
				sibling_el_ids[i] = airway_data[sibling_el_ids[i]].get_local_elem_number();//n_initial_3d_elem;
				//std::cout << "sibling_el_ids = " << sibling_el_ids[i] << std::endl;
			}

			int primary_sibling_el_idx = -1;
			if(!is_daughter_1)
			{
				//primary_sibling_el_idx = airway_data[current_1d_el_idx].get_primary_sibling() + n_initial_3d_elem;
				primary_sibling_el_idx = airway_data[airway_data[current_1d_el_idx].get_primary_sibling()].get_local_elem_number();// + n_initial_3d_elem;
			}
			
			const double l = airway_data[current_1d_el_idx].get_length();	//nondimensionalised length
			const double r = airway_data[current_1d_el_idx].get_radius(); //nondimensionalised radius
			int generation = airway_data[current_1d_el_idx].get_generation();

			bool has_sibling = false;
			if(sibling_el_ids.size() > 0)
				has_sibling = true;

			//very dirty hack
			if(!has_parent)
				parent_el_idx = current_el_idx;

			//if(primary_sibling_el_idx  < n_initial_3d_elem)
			if(primary_sibling_el_idx  < 0)
				primary_sibling_el_idx  = current_el_idx;

			// don't need cause we check with size of sibling_el_ids
			//if(daughter_1_el_idx  < n_initial_3d_elem)
			if(daughter_1_el_idx  < 0)
				daughter_1_el_idx  = current_el_idx;

			
			const Elem* parent_elem = mesh.elem(parent_el_idx);
			const Elem* daughter_1_elem = mesh.elem(daughter_1_el_idx);

			std::vector<const Elem*> sibling_elems;
			for(unsigned int i=0; i<sibling_el_ids.size(); i++)
				sibling_elems.push_back(mesh.elem(sibling_el_ids[i]));

			//const Elem* primary_sibling_elem = mesh.elem(primary_sibling_el_idx); // unused
			//const Elem* sibling_elem = mesh.elem(sibling_el_idx);

			/*
			if(parent_el_idx >= 0)
				parent_elem = mesh.elem(parent_el_idx);

			if(daughter_1_el_idx >= 0)
				daughter_1_elem = mesh.elem(daughter_1_el_idx);

			if(daughter_2_el_idx >= 0)
				daughter_2_elem = mesh.elem(daughter_2_el_idx);
			*/

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_q, q_var);
      dof_map.dof_indices (parent_elem, dof_indices_parent_p, p_var);
      dof_map.dof_indices (parent_elem, dof_indices_parent_q, q_var);
      dof_map.dof_indices (daughter_1_elem, dof_indices_daughter_1_p, p_var);

			dof_indices_siblings_q.resize(sibling_el_ids.size());
			dof_indices_siblings_p.resize(sibling_el_ids.size());

			for(unsigned int i=0; i< sibling_el_ids.size(); i++)
			{
		    dof_map.dof_indices (sibling_elems[i], dof_indices_siblings_q[i], q_var);
		    dof_map.dof_indices (sibling_elems[i], dof_indices_siblings_p[i], p_var);
			}


			for(unsigned int i=0; i<dof_indices_p.size(); i++)
				zero_d_dofs.push_back(dof_indices_p[i]);

			for(unsigned int i=0; i<dof_indices_q.size(); i++)
				zero_d_dofs.push_back(dof_indices_q[i]);

			// get the acinar dofs and other acinar stuff
			unsigned int acinar_v_dof_idx = 0;
			double local_acinar_compliance = 0.;
			int acinar_elem_number = 0.;
			if(acinar_model == 1)
			{
				// first check if it has an acinar i.e. terminal
				acinar_elem_number = airway_data[current_1d_el_idx].get_acinar_elem();
				if(acinar_elem_number > -1)
				{
					// now get the correct element
					const Elem* acinar_elem = mesh.elem(acinar_elem_number);
//					std::cout << "v_var = " << v_var << std::endl;
//					std::cout << "acinar_elem_number = " << acinar_elem_number << std::endl;
      		dof_map.dof_indices (acinar_elem, dof_indices_acinar);
//					std::cout << "hey babe" << std::endl;
//					std::cout << "dof_indices.size() = " << dof_indices_acinar.size() << std::endl;
//					std::cout << "hey babe" << std::endl;
      		dof_map.dof_indices (acinar_elem, dof_indices_v, v_var);
//					std::cout << "dof_indices_v.size() = " << dof_indices_v.size() << std::endl;
					acinar_v_dof_idx = dof_indices_v[0];
					
					local_acinar_compliance = airway_data[current_1d_el_idx].get_local_acinar_compliance();
//					std::cout << "local_acinar_compliance = " << local_acinar_compliance << std::endl;
				}

			}
			







      //const unsigned int n_dofs   = 4;	// unused
      //const unsigned int n_p_dofs = 2;	// unused

			// some parameters

			// use the current solution or the previous time step's solution
			// if no compliance, flow rate is the same at the beginning and end of the pipe
			double average_flow = 0.;
			if(semi_implicit_1d)
			{
				if(compliance_1d)
					average_flow = ((*system->old_local_solution)(dof_indices_q[0]) + (*system->old_local_solution)(dof_indices_q[1]))/2.0;
				else
					average_flow = (*system->old_local_solution)(dof_indices_q[0]);
			}
			else
			{
				if(compliance_1d)
					average_flow = (system->current_solution(dof_indices_q[0]) + system->current_solution(dof_indices_q[1]))/2.0;
				else
					average_flow = system->current_solution(dof_indices_q[0]);
			}


			// construct dimensionalised parameters
			double R_poi = airway_data[current_1d_el_idx].get_poiseuille_resistance();
			double R_deriv = 0.;

			double C = 0.;
			if(compliance_1d)
			{
				C = airway_data[current_1d_el_idx].get_airway_compliance();
			}

			double I = 0.;
			if(inertance_1d)
			{
				C = airway_data[current_1d_el_idx].get_inertance();

			}


			double R = 0.;
			

/*
			std::cout << "length_scale = " << length_scale << std::endl;
			std::cout << "velocity_scale = " << velocity_scale << std::endl;
			std::cout << "R = " << R << std::endl;
			std::cout << "viscosity = " << viscosity << std::endl;
			std::cout << "density = " << density << std::endl;
			std::cout << "l = " << l << std::endl;
			std::cout << "r = " << r << std::endl;
			//std::cout << "r^4 = " << pow(r,4.0) << std::endl;
*/


			//1d resistance types
			if(resistance_type_1d == 0)
			{
				//poiseuille resistance
				R=R_poi;
				airway_data[current_1d_el_idx].set_poiseuille(true);
			}
			else if(resistance_type_1d == 1)
			{
				//pedley resistance

				double reynolds_number_0d = 2*density*velocity_scale*length_scale*fabs(average_flow)/(viscosity*M_PI*r);
				double pedley_factor = 0.327 * sqrt(reynolds_number_0d*2.*r/l);

				if(generation < 5)
				{
					//std::cout << "generation = " << generation << std::endl;
					//std::cout << "reynolds_number_0d = " << reynolds_number_0d << std::endl;
					//std::cout << "pedley factor = " << 0.327 * sqrt(reynolds_number_0d*2.*r/l) << std::endl;
					//std::cout << "length = " << l << std::endl;
					//std::cout << "radius = " << r << std::endl;
				}

				if(generation == 0 && first_gen_poiseuille)
				{
					std::cout << "using first gen poi" << std::endl;
					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
				}
				else if(pedley_factor > 1.0)
				{

					R = pedley_factor * R_poi;
					if(newton_1d)
					{
						R_deriv = 0.327 * 2*density/(viscosity*M_PI*l) / sqrt(reynolds_number_0d*2.*r/l) * R_poi;
					}

					airway_data[current_1d_el_idx].set_poiseuille(false);
				}
				else
				{

					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
					
				}
			}
			else if(resistance_type_1d == 2)
			{


				double reynolds_number_0d = 2*density*velocity_scale*length_scale*fabs(average_flow)/(viscosity*M_PI*r);

				double gamma = 0;
				if(generation < 0)		//alveolar
					gamma = ertbruggen_gamma[8];
				else if(generation < 8)
					gamma = ertbruggen_gamma[generation];
				else
					gamma = ertbruggen_gamma[8];

				// van ertbruggen resistance
				double ertbruggen_factor = gamma * sqrt(reynolds_number_0d*2.*r/l);
				
				if(generation == 0 && first_gen_poiseuille)
				{
					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
				}
				else if(ertbruggen_factor > 1.0)
				{

					R = ertbruggen_factor * R_poi;
					if(newton_1d)
						R_deriv =  gamma * 2*density/(viscosity*M_PI*l) / sqrt(reynolds_number_0d*2.*r/l) * R_poi;

					airway_data[current_1d_el_idx].set_poiseuille(false);
				}
				else
				{
					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
				}
			}
			else if(resistance_type_1d == 3)
			{
	

				double reynolds_number_0d = 2*density*velocity_scale*length_scale*fabs(average_flow)/(viscosity*M_PI*r);

				double constant = 0;
				//lobar bronchi defined as generation less than equal 2, but doesn't quite get it, oops
				if(generation > 2 || generation < 0)
					constant = 1.0;
				else
					constant = 3.4;

				// van ertbruggen resistance
				if(es->parameters.get<unsigned int> ("t_step") != 1)
				{
					R = (constant + 2.1e-3 * reynolds_number_0d) * R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(false);
				}
				else
				{
					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
				}
			}
			else if(resistance_type_1d == 4)
			{
				if(generation == 0)
					R = resistance_0;
				else if(generation == 1)
					R = resistance_1;
				else if(generation == 2)
					R = resistance_2;
				else if(generation == 3)
					R = resistance_3;
				else if(generation == 4)
					R = resistance_4;
				else
				{
					std::cout << "hard-coded resistances not implmeneted for meshes with more than 4 generations." << std::endl;
					std::cout << "Exiting..." << std::endl;
					std::exit(0);
				}

				std::cout << "gen " << generation << std::endl;
				std::cout << "resistance = " << R << std::endl;
			}
			else if(resistance_type_1d == 5)
			{
				if(generation == 0 && current_1d_el_idx == 0)
					R = resistance_0_a;
				else if(generation == 0 && current_1d_el_idx == 1)
					R = resistance_0_b;
				else
				{
					std::cout << "hard-coded asymmetric resistances not implmeneted for meshes with more than 1 generations." << std::endl;
					std::cout << "Exiting..." << std::endl;
					std::exit(0);
				}

				std::cout << "gen " << generation << std::endl;
				std::cout << "resistance = " << R << std::endl;
			}
			else
			{
					R=R_poi;
					airway_data[current_1d_el_idx].set_poiseuille(true);
			}

			// set the resistance so that can be output
			airway_data[current_1d_el_idx].set_resistance(R);


			// set the resistance to be used by the monolithic preconditioner
			if(generation == 0)
			{
				if(R > mono_preconditioner_resistance_scaling)
				{
					mono_preconditioner_resistance_scaling = R;
				}
			}

			double old_p0 = system->old_solution(dof_indices_p[0]);
			//double old_p1 = system->old_solution(dof_indices_p[1]);	// unused
			//double old_q0 = system->old_solution(dof_indices_q[0]);	// unused
			double old_q1 = system->old_solution(dof_indices_q[1]);
		
			double p0 = system->current_solution(dof_indices_p[0]);
			double p1 = system->current_solution(dof_indices_p[1]);
			double q0 = system->current_solution(dof_indices_q[0]);
			double q1 = system->current_solution(dof_indices_q[1]);

			//std::cout << "q0 = " << q0 << std::endl;
			//std::cout << "q1 = " << q1 << std::endl;
			//std::cout << "average_flow = " << average_flow << std::endl;

			std::vector<double> siblings_q0(sibling_el_ids.size());
			std::vector<double> siblings_p0(sibling_el_ids.size());
			for(unsigned int i=0; i<sibling_el_ids.size(); i++)
			{
				siblings_q0[i] = system->current_solution(dof_indices_siblings_q[i][0]);
				siblings_p0[i] = system->current_solution(dof_indices_siblings_p[i][0]);
			}

			double parent_p1 = 0.;
			double parent_q1 = 0.;
			if(parent_el_idx != current_el_idx)
			{
				parent_p1 = system->current_solution(dof_indices_parent_p[1]);
				parent_q1 = system->current_solution(dof_indices_parent_q[1]);
			}
		

			double v0 = 0.;
			if(acinar_model == 1 && acinar_elem_number > -1)
				v0 = system->current_solution(acinar_v_dof_idx);

			double daughter_1_p0 = 0.;
			if(daughter_1_el_idx != current_el_idx)
				daughter_1_p0 = system->current_solution(dof_indices_daughter_1_p[0]);

			//bool compute_using_monomials = false;	// unused
			{
				// first we decide which equations go where based on the boundary conditions that exist

				// the first equation is the compliance equation (mass conservation)
				// the second equation is the resistance equation (flow = grad(p))
				// the third equation is always an equation for the inflow boundary
				// the fourth equation is always an equation for the outflow boundary 
				std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,0);

				unsigned int eqn_1_dof = 0;
				unsigned int eqn_2_dof = 0;
				unsigned int eqn_3_dof = 0;		//put influx condition in p[0], what why?
				unsigned int eqn_4_dof = 0;

				//default
				if(is_daughter_1)
				{
					eqn_1_dof = dof_indices_q[1];
					eqn_2_dof = dof_indices_p[0];
					eqn_3_dof = dof_indices_q[0];
					eqn_4_dof = dof_indices_p[1];
				}
				else
				{
					eqn_1_dof = dof_indices_q[0];
					eqn_2_dof = dof_indices_q[1];
					eqn_3_dof = dof_indices_p[0];
					eqn_4_dof = dof_indices_p[1];
				}

				//inflow bc - must be a major/daughter_1 branch
				// takes care of all inflow bc and inflow and outflow bc
				if(parent_el_idx == current_el_idx)
				{
					//flow prescribed at inflow
					if(bc_type_1d == 0)
					{
						// these actually remain the same
						// eqn_3_dof = dof_indices_q[0];		//put influx condition in eqn_3, no change
						// if we also have the pressure prescribed at the outflow then eqn_4 is also unchanged, still p_out

				
					}
					//pressure prescribed at inflow
					else if(bc_type_1d == 1)
					{
						eqn_3_dof = dof_indices_p[0];		//put pressure inflow condition in p[0]
						eqn_1_dof = dof_indices_q[0];		// need to put eqn_1 in q_in row because the pressure bc has taken it's place
						eqn_2_dof = dof_indices_q[1];		// need to put eqn_2 in q_out row because the q_in has taken it
						// if we also have the pressure prescribed at the outflow then eqn_4 is also unchanged, still p_out
					}
				}

				//only outflow bc (we only consider pressure outflow boundary conditions which leave the equations unchanged)				
						
				// ********** PUTTING THE EQUATIONS IN THE MATRIX AND RHS ********* //

				// **** equation 1 - compliance eqn	
				system->get_vector("Residual LHS").add(eqn_1_dof,q1);
				system->get_vector("Residual LHS").add(eqn_1_dof,-q0);

				if(unsteady)
				{
					system->get_vector("Residual LHS").add(eqn_1_dof,C/dt*p0);
				}

				// **** equation 2 - resistance eqn
				// don't need derivative in residual retard
				system->get_vector("Residual LHS").add(eqn_2_dof,R*q1);

				system->get_vector("Residual LHS").add(eqn_2_dof,p1);
				system->get_vector("Residual LHS").add(eqn_2_dof,-p0);

				if(unsteady)
				{
					system->get_vector("Residual LHS").add(eqn_2_dof,I/dt*q1);
				}

				// **** equation 3 - BC on inflow
				if(parent_el_idx == current_el_idx)
				{
					//std::cout << "parent_el_idx = " << parent_el_idx << std::endl;
					//std::cout << "current_el_idx = " << current_el_idx << std::endl;
					// flow
					if(bc_type_1d == 0)
					{
						// now there is no zero on the 3rd equation
						// equation 3 - inflow bc multiplied by dt so is like the 3d eqn - not anymore it isn't
	
						// daughter 1 handles the flux conservation to sibling and parent (1 equation)
						if(is_daughter_1)
						{
							
							system->get_vector("Residual LHS").add(eqn_3_dof,-mono_flow_rate_penalty_param*q0);

							//std::cout << "flux dof in 0d = " << eqn_3_dof << std::endl;
							if(has_sibling)
							{
								//std::cout << "num siblings = " << sibling_el_ids.size() << std::endl;
								//std::cout << "sibling id daughter1 = " << sibling_el_ids[0] << std::endl;
								//std::cout << "current_1d_el_idx  = " << current_1d_el_idx  << std::endl;
								for(unsigned int i=0; i<sibling_el_ids.size(); i++)
									system->get_vector("Residual LHS").add(eqn_3_dof,-mono_flow_rate_penalty_param*siblings_q0[i]);
							}
			
							// boundary condition should be on side 0
							std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,0);
							if(boundary_ids.size() == 0)
								std::cout << "error, boundary does not have boundary id as expected, el idx = " << current_el_idx << std::endl;

							//std::cout << "flux_values[" << boundary_ids[0] <<"] = " << flux_values[boundary_ids[0]] << std::endl;

						}
						else	// if we are a daughter 2 then we have a sibling and need to apply the pressure continuity condition
						{
							system->get_vector("Residual LHS").add(eqn_3_dof,p0);
							for(unsigned int i=0; i<sibling_el_ids.size(); i++)
								system->get_vector("Residual LHS").add(eqn_3_dof,-siblings_p0[i]);
						}
					}
					// pressure-pressure type boundary condition
					else if(bc_type_1d == 1)
					{
						system->get_vector("Residual LHS").add(eqn_3_dof,p0);
					}

				}
				else
				{
					if(is_daughter_1)
					{
						// conservation of flux
						system->get_vector("Residual LHS").add(eqn_3_dof,parent_q1);
						system->get_vector("Residual LHS").add(eqn_3_dof,-q0);
						for(unsigned int i=0; i<sibling_el_ids.size(); i++)
							system->get_vector("Residual LHS").add(eqn_3_dof,-siblings_q0[i]);
					}
					else
					{
						// parent pressure cont
						system->get_vector("Residual LHS").add(eqn_3_dof,p0);
						system->get_vector("Residual LHS").add(eqn_3_dof,-parent_p1);
					}

				}


				// **** equation 4 - BC on outflow
				if(daughter_1_el_idx  == current_el_idx)
				{

					if(acinar_model == 1)
					{
						// acinar compliance equation
						system->get_vector("Residual LHS").add(eqn_4_dof,p1);
						system->get_vector("Residual LHS").add(eqn_4_dof,-1.0/local_acinar_compliance * v0);
						
					}
					else
					{
						//apply pressure
						system->get_vector("Residual LHS").add(eqn_4_dof,p1);
					}

				}
				else
				{
					//pressure continuity to daughter_1
					system->get_vector("Residual LHS").add(eqn_4_dof,p1);
					system->get_vector("Residual LHS").add(eqn_4_dof,-daughter_1_p0);
				}
			
			}

    }// end of subdomains_1d conditional
		else if(std::find(subdomains_acinar.begin(), subdomains_acinar.end(), subdomain_id) != subdomains_acinar.end())
		{
//						std::cout << "gru2" << std::endl;
			// need to find the element that this acinar is assoc with
			const int current_el_idx = elem->id();
			const unsigned int acinar_id = elem_to_acinar[current_el_idx];
			//const unsigned int parent_el_idx = acinar_to_airway[acinar_id] + n_initial_3d_elem;
			const unsigned int parent_el_idx = airway_data[acinar_to_airway[acinar_id]].get_local_elem_number();// + n_initial_3d_elem;
			const Elem* parent_elem = mesh.elem(parent_el_idx);

//						std::cout << "gru2" << std::endl;
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (parent_elem, dof_indices_q, q_var);

//						std::cout << "gru2" << std::endl;
			// get old volume value
			double v_old = (*system->old_local_solution)(dof_indices_v[0]);
			double q1 = system->current_solution(dof_indices_q[1]);
			double v0 = system->current_solution(dof_indices_v[0]);

			//std::cout << "v_old = " << v_old << std::endl;

//						std::cout << "gru2" << std::endl;
			// volume-flow equation
			system->get_vector("Residual LHS").add(dof_indices_v[0],1.0/dt*v0);
//						std::cout << "gru2" << std::endl;
			system->get_vector("Residual LHS").add(dof_indices_v[0],-q1);
		}
			
	}// end of element loop



	// make the rhs -Ax+f
	system->rhs->close();
	system->rhs->zero();
	system->rhs->add(-1.0,system->get_vector("Residual LHS"));

	if(es->parameters.get<bool>("increment_boundary_conditions") && es->parameters.get < unsigned int >("nonlinear_iteration") == 1)
		system->rhs->add(system->get_vector("Forcing Vector BC"));
	else
		system->rhs->add(system->get_vector("Forcing Vector"));
/*
		std::cout << "***** residual_lhs_norm = " << system->get_vector("Residual LHS").l2_norm() << std::endl;
		std::cout << "***** forcing_vector_norm = " << system->get_vector("Forcing Vector").l2_norm() << std::endl;
	*/

	std::cout << "End 0D assembly" << std::endl;
  return;
}






// total flux through boundary_id boundaries, or at the second node of an element mid mesh
double NavierStokesAssembler::calculate_flux (const int boundary_id, const int mid_mesh_element)
{
	double flux = 0.0;


  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es->get_mesh();

  // The dimension that we are running
  //const unsigned int dim = mesh.mesh_dimension();	// unused

	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns1d");
	}

  // Numeric ids corresponding to each variable in the system
  const unsigned int q_var = system->variable_number ("Q");

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples->
  const DofMap & dof_map = system->get_dof_map();

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_q;


  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
						

  for ( ; el != end_el; ++el)
  {
		const Elem* elem = *el;
		//only concerned with 1d elements
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
			if(mid_mesh_element < 0)
			{
		    for (unsigned int s=0; s<elem->n_sides(); s++)
				{

					//for some reason it is natural to have more than one boundary id per side or even node
					std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

					if(boundary_ids.size() > 0) 
					{ 
						// dirichlet boundary conditions
						if(boundary_ids[0] == boundary_id)
						{	     
					
							// Get the degree of freedom indices for the
							// current element.  These define where in the global
							// matrix and right-hand-side this element will
							// contribute to.
							dof_map.dof_indices (elem, dof_indices);
							dof_map.dof_indices (elem, dof_indices_q, q_var);

							if(s == 0)
								flux += (*system->solution) (dof_indices_q[0]);
							else
								flux += (*system->solution) (dof_indices_q[1]);

						}
					}
				}
			}
			else
			{
				
				const int current_el_idx = elem->id();
				//unsigned int current_1d_el_idx = current_el_idx -	n_initial_3d_elem;
				unsigned int current_1d_el_idx = elem_to_airway[current_el_idx];// -	n_initial_3d_elem;
				if(current_1d_el_idx == (unsigned int)mid_mesh_element)
				{
					dof_map.dof_indices (elem, dof_indices);
					dof_map.dof_indices (elem, dof_indices_q, q_var);

					flux += (*system->solution) (dof_indices_q[1]);
					
				}
			}
		}
	}

	//std::cout << "flux calculated is " << flux << std::endl;

	es->comm().sum(flux);
	//std::cout << "flux calculated is " << flux << std::endl;
			
	return flux;
}

// average pressure over boundary_id boundaries
double NavierStokesAssembler::calculate_pressure (const int boundary_id, const int mid_mesh_element)
{
	double pressure = 0.0;

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es->get_mesh();

  // The dimension that we are running
  //const unsigned int dim = mesh.mesh_dimension();	// unused

	
	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(coupled)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns1d");
	}

  // Numeric ids corresponding to each variable in the system
  const unsigned int p_var = system->variable_number ("P");

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap & dof_map = system->get_dof_map();

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_p;


  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	unsigned int total_nodes = 0;
						
  for ( ; el != end_el; ++el)
  {
		const Elem* elem = *el;
		//only concerned with 1d elements
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
			if(mid_mesh_element < 0)
			{
		    for (unsigned int s=0; s<elem->n_sides(); s++)
				{
					//for some reason it is natural to have more than one boundary id per side or even node
					std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,s);

					if(boundary_ids.size() > 0) 
					{ 
						// dirichlet boundary conditions
						if(boundary_ids[0] == boundary_id)
						{	     
							total_nodes++;
							// Get the degree of freedom indices for the
							// current element.  These define where in the global
							// matrix and right-hand-side this element will
							// contribute to.
							dof_map.dof_indices (elem, dof_indices);
							dof_map.dof_indices (elem, dof_indices_p, p_var);

							if(s == 0)
								pressure += (*system->solution) (dof_indices_p[0]);
							else
								pressure += (*system->solution) (dof_indices_p[1]);
						}
					}
				}
			}
			else
			{

				const int current_el_idx = elem->id();
				//unsigned int current_1d_el_idx = current_el_idx -	n_initial_3d_elem;
				unsigned int current_1d_el_idx = elem_to_airway[current_el_idx];// -	n_initial_3d_elem;
				if(current_1d_el_idx == (unsigned int)mid_mesh_element)
				{
					total_nodes++;
					dof_map.dof_indices (elem, dof_indices);
					dof_map.dof_indices (elem, dof_indices_p, p_var);

					pressure += (*system->solution) (dof_indices_p[1]);
					
				}
			}
		}
	}

	es->comm().sum(pressure);
	es->comm().sum(total_nodes);

	// get the average pressure of course ;) total_nodes should be 1
	pressure = pressure/total_nodes;

	//std::cout << "total_nodes = " << total_nodes << std::endl;
	//std::cout << "pressure_calculated = " << pressure << std::endl;

			
	return pressure;
}


void NavierStokesAssembler::init_bc(std::vector<double> _flux_values,std::vector<double> _pressure_values)
{
	std::cout << "inside init bc 1d" << std::endl;
	pressure_values = _pressure_values;
	flux_values = _flux_values;

	

	// cause backwards from 3D
	/*
	for(unsigned int i=0; i<flux_values.size(); i++)
		std::cout << "flux val[" << i << "] = " << flux_values[i] << std::endl;
	*/
}


double NavierStokesAssembler::calculate_pleural_pressure(double time)
{
	double pl = 0.;

	unsigned int pl_type = es->parameters.get<unsigned int>("pl_type");
	double pl_mean = es->parameters.get<double>("pl_mean");
	double pl_range = es->parameters.get<double>("pl_range");
	double pl_period = es->parameters.get<double>("pl_period");
	double pl_volume_mean = es->parameters.get<double>("pl_volume_mean");
	double pl_volume_range = es->parameters.get<double>("pl_volume_range");
	double pl_total_resistance = es->parameters.get<double>("pl_total_resistance");
	double total_acinar_compliance = es->parameters.get<double>("total_acinar_compliance");


	if(pl_type == 0)
		pl = pl_mean + pl_range * sin(2*M_PI*time/pl_period);
	else if(pl_type == 1)
		pl = pl_mean + pl_range * cos(2*M_PI*time/pl_period);
	else if(pl_type == 2)
		pl = pl_mean - pl_range * cos(2*M_PI*time/pl_period);
	else if(pl_type == 3)
		pl = pl_mean;
	else if(pl_type == 4)
	{	
		double pl_volume = pl_volume_mean - pl_volume_range * cos(2*M_PI*time/pl_period);
		double pl_volume_derivative = pl_volume_range * 2*M_PI/pl_period * sin(2*M_PI*time/pl_period);
		pl = -1./total_acinar_compliance * pl_volume - pl_total_resistance * pl_volume_derivative;
	}

	//std::cout << "calculating pleural pressure" << std::endl;
	//std::cout << "time = " << time << std::endl;
	//std::cout << "pl = " << pl << std::endl;
	//std::cout << "cos = " << cos(2*M_PI*time/pl_period) << std::endl;
	//std::cout << "cos arg = " << 2*M_PI*time/pl_period << std::endl;

	return pl;
}






