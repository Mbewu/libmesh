
#include "navier_stokes_assembler.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

// The matrix assembly function to be called at each time step to
// prepare for the linear solve. This uses a picard type linearisation
// as opposed to the newton type linearisation used in the 
// assemble_stokes function. Well now assemble_stokes is not working in 3D
void NavierStokesAssembler::assemble ()
{
	
  if(es->parameters.get<bool>("0D"))
	{
		assemble_stokes_steady_0D ();
	}
	else
		assemble_stokes_1D ();
}

// The matrix assembly function to be called at each time step to
// prepare for the linear solve. This uses a picard type linearisation
// as opposed to the newton type linearisation used in the 
// assemble_stokes function. Well now assemble_stokes is not working in 3D
void NavierStokesAssembler::assemble_stokes_steady_0D ()
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
  const unsigned int q_var = system->variable_number ("Q");

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
	//const double time = system->time;	// unused
	const bool compliance_1d = es->parameters.get<bool>("compliance_1d");
	const bool inertance_1d = es->parameters.get<bool>("inertance_1d");
  const unsigned int resistance_type_1d = es->parameters.get<unsigned int>("resistance_type_1d");

	// if we are doing a reynolds_number_calculation we need to get the viscosity from reynolds_number
	// length scale and velocity scale being 1
	if(es->parameters.get<bool> ("reynolds_number_calculation"))
		viscosity = 1./reynolds_number;

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

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the \p active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

								
  for ( ; el != end_el; ++el)
  {
		const Elem* elem = *el;
		unsigned int subdomain_id = elem->subdomain_id();
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), subdomain_id) != subdomains_1d.end())
		{
			//element data object starts numbering from 0
			//and the values referenced in it also do so need to take this into account

			const int current_el_idx = elem->id();
			unsigned int current_1d_el_idx = current_el_idx -	n_initial_3d_elem;
			bool has_parent = airway_data[current_1d_el_idx].has_parent();
			int parent_el_idx = 0;
			if(has_parent)
				parent_el_idx = (int)airway_data[current_1d_el_idx].get_parent() + n_initial_3d_elem;
			std::vector<unsigned int> daughter_el_ids = airway_data[current_1d_el_idx].get_daughters();
			int daughter_1_el_idx = -1;
			if(daughter_el_ids.size() > 0)
				daughter_1_el_idx = daughter_el_ids[0] + n_initial_3d_elem;
			
			bool is_daughter_1 = airway_data[current_1d_el_idx].get_is_daughter_1();	//this is a bool duh!
			std::vector<unsigned int> sibling_el_ids = airway_data[current_1d_el_idx].get_siblings();
			for(unsigned int i=0; i<sibling_el_ids.size(); i++)
			{
				sibling_el_ids[i] += n_initial_3d_elem;
				//std::cout << "sibling_el_ids = " << sibling_el_ids[i] << std::endl;
			}

			int primary_sibling_el_idx = -1;
			if(!is_daughter_1)
				primary_sibling_el_idx = airway_data[current_1d_el_idx].get_primary_sibling() + n_initial_3d_elem;
			
			const double l = airway_data[current_1d_el_idx].get_length();	//nondimensionalised length
			const double r = airway_data[current_1d_el_idx].get_radius(); //nondimensionalised radius
			int generation = airway_data[current_1d_el_idx].get_generation();

			bool has_sibling = false;
			if(sibling_el_ids.size() > 0)
				has_sibling = true;

			//very dirty hack
			if(!has_parent)
				parent_el_idx = current_el_idx;

			if(primary_sibling_el_idx  < n_initial_3d_elem)
				primary_sibling_el_idx  = current_el_idx;

			// don't need cause we check with size of sibling_el_ids
			if(daughter_1_el_idx  < n_initial_3d_elem)
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

      //const unsigned int n_dofs   = 4;	// unused
      //const unsigned int n_p_dofs = 2;	// unused

			// some parameters
			// these are the conventional parameters divided by density because of the pressure scaling
			double average_flow = (system->current_solution(dof_indices_q[0]) + system->current_solution(dof_indices_q[1]))/2.0;
			// average_flow /= 2.0; //huh?


			// construct dimensionalised parameters
			double reynolds_number_0d = 2*density*velocity_scale*length_scale*fabs(average_flow)/(viscosity*M_PI*r);
			const double t = zeta_1*pow(r*length_scale,2) + zeta_2*(r*length_scale) + zeta_3;	//dimensionalised thickness
			double R = 8*l*viscosity/(M_PI*pow(length_scale,3)*pow(r,4.0));	//dunno why this is pow 4??
			double C = pow(length_scale,4)*2*l*pow(r,3.0)/(E*t);
			double I = l*density/(M_PI*r*r);

			// construct nondimensionalised parameters
			R = R*pow(length_scale,2)/(density*velocity_scale);
			C = density*pow(velocity_scale,2)*C/length_scale;
			I = I*length_scale/density;

		//std::cout << "R = " << R << std::endl;
		//std::cout << "viscosity = " << viscosity << std::endl;
		//std::cout << "l = " << l << std::endl;
		//std::cout << "r = " << r << std::endl;
		//std::cout << "r^4 = " << pow(r,4.0) << std::endl;


			// use given reynolds number is reynolds number calc
			if(es->parameters.get<bool> ("reynolds_number_calculation"))
				reynolds_number = es->parameters.get<double> ("reynolds_number");


			//1d resistance types
			if(resistance_type_1d == 0)
			{
				//poiseuille resistance
				R=R;
				airway_data[current_1d_el_idx].set_poiseuille(true);
			}
			else if(resistance_type_1d == 1)
			{
				if(generation < 5)
				{
					//std::cout << "generation = " << generation << std::endl;
					//std::cout << "reynolds_number_0d = " << reynolds_number_0d << std::endl;
					//std::cout << "pedley factor = " << 0.327 * sqrt(reynolds_number_0d*2.*r/l) << std::endl;
					//std::cout << "length = " << l << std::endl;
					//std::cout << "radius = " << r << std::endl;
				}

				//pedley resistance
				// we also want to set here whether we used a pedley model
				if((es->parameters.get<unsigned int> ("t_step") != 1 || es->parameters.get<unsigned int> ("nonlinear_iteration_1d") != 1) && 0.327 * sqrt(reynolds_number_0d*2.*r/l) > 1.0)
				{
					
					R = 0.327 * sqrt(reynolds_number_0d*2.*r/l) * R;
					airway_data[current_1d_el_idx].set_poiseuille(false);
					std::cout << "gen = " << generation << ", using pedley" << std::endl;

				}
				else
				{
					R=R;
					airway_data[current_1d_el_idx].set_poiseuille(true);
				}
			}
			else if(resistance_type_1d == 2)
			{
				double gamma = 0;
				if(generation < 0)		//alveolar
					gamma = ertbruggen_gamma[8];
				else if(generation < 8)
					gamma = ertbruggen_gamma[generation];
				else
					gamma = ertbruggen_gamma[8];

				// van ertbruggen resistance
				if(es->parameters.get<unsigned int> ("t_step") != 1 && gamma * sqrt(reynolds_number_0d*2.*r/l) > 1.0)
				{
					R = gamma * sqrt(reynolds_number_0d*2.*r/l) * R;
					airway_data[current_1d_el_idx].set_poiseuille(false);
				}
				else
				{
					R=R;
					airway_data[current_1d_el_idx].set_poiseuille(true);
				}
			}
			else if(resistance_type_1d == 3)
			{
				double constant = 0;
				//lobar bronchi defined as generation less than equal 2, but doesn't quite get it, oops
				if(generation > 2 || generation < 0)
					constant = 1.0;
				else
					constant = 3.4;

				// van ertbruggen resistance
				if(es->parameters.get<unsigned int> ("t_step") != 1)
				{
					R = (constant + 2.1e-3 * reynolds_number_0d) * R;
					airway_data[current_1d_el_idx].set_poiseuille(false);
				}
				else
				{
					R=R;
					airway_data[current_1d_el_idx].set_poiseuille(true);
				}
			}
			else
			{
					R=R;
					airway_data[current_1d_el_idx].set_poiseuille(true);
			}

			// unused
			//double flow_scale = es->parameters.get<double>("velocity_scale") * 
			//									pow(es->parameters.get<double>("length_scale"),2.0);
			//double mean_pressure_scale = es->parameters.get<double>("density") *
			//												pow(es->parameters.get<double>("velocity_scale"),2.0);


			//change for correctness, to change
			//R = 8*R;


			const double in_pressure = es->parameters.get<double> ("in_pressure");
			const double out_pressure = 0.0;

			//std::cout << "resistance 1d = " << R << std::endl;

			/*
			std::cout << "I = " << I << std::endl;
			std::cout << "C = " << C << std::endl;
			std::cout << "R = " << R << std::endl;
			std::cout << "t = " << t << std::endl;
			std::cout << "r = " << r << std::endl;
			*/
			//I = 0;
			if(!compliance_1d)
				C = 0;
			if(!inertance_1d)
				I = 0;

			double old_p0 = system->old_solution(dof_indices_p[0]);
			//double old_p1 = system->old_solution(dof_indices_p[1]);	// unused
			//double old_q0 = system->old_solution(dof_indices_q[0]);	// unused
			double old_q1 = system->old_solution(dof_indices_q[1]);
		
		
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
					if(es->parameters.get<unsigned int> ("bc_type_1d") == 0)
					{
						// these actually remain the same
						// eqn_3_dof = dof_indices_q[0];		//put influx condition in eqn_3, no change
						// if we also have the pressure prescribed at the outflow then eqn_4 is also unchanged, still p_out

				
					}
					//pressure prescribed at inflow
					else if(es->parameters.get<unsigned int> ("bc_type_1d") == 1)
					{
						eqn_3_dof = dof_indices_p[0];		//put pressure inflow condition in p[0]
						eqn_1_dof = dof_indices_q[0];		// need to put eqn_1 in q_in row because the pressure bc has taken it's place
						eqn_2_dof = dof_indices_q[1];		// need to put eqn_2 in q_out row because the q_in has taken it
						// if we also have the pressure prescribed at the outflow then eqn_4 is also unchanged, still p_out
					}
				}

				//only outflow bc (we only consider pressure outflow boundary conditions which leave the equations unchanged)				
						
				// ********** PUTTING THE EQUATIONS IN THE MATRIX AND RHS ********* //

				// **** equation 1 - compliance
	    	add_to_matrices (system,eqn_1_dof,dof_indices_q[1],1.0);
	    	add_to_matrices (system,eqn_1_dof,dof_indices_q[0],-1.0);

				if(unsteady)
				{
	    		add_to_matrices (system,eqn_1_dof,dof_indices_p[0],C/dt);
	    		system->rhs->add (eqn_1_dof,C*old_p0/dt);
				}

				// **** equation 2 - resistance
	    	add_to_matrices (system,eqn_2_dof,dof_indices_q[1],R);
	    	add_to_matrices (system,eqn_2_dof,dof_indices_p[1],1.0);
	    	add_to_matrices (system,eqn_2_dof,dof_indices_p[0],-1.0);

				if(unsteady)
				{
	    		add_to_matrices (system,eqn_2_dof,dof_indices_q[1],I/dt);
	    		system->rhs->add (eqn_2_dof,I*old_q1/dt);
				}

				// **** equation 3 - BC on inflow
				if(parent_el_idx == current_el_idx)
				{
					//std::cout << "parent_el_idx = " << parent_el_idx << std::endl;
					//std::cout << "current_el_idx = " << current_el_idx << std::endl;
					// flow
					if(es->parameters.get<unsigned int> ("bc_type_1d") == 0)
					{
						// now there is no zero on the 3rd equation
						// equation 3 - inflow bc multiplied by dt so is like the 3d eqn - not anymore it isn't
	
						if(is_daughter_1)
						{
							
							add_to_matrices (system,eqn_3_dof,dof_indices_q[0],-1.0);

							std::cout << "flux dof in 0d = " << eqn_3_dof << std::endl;
							if(has_sibling)
							{
								//std::cout << "num siblings = " << sibling_el_ids.size() << std::endl;
								//std::cout << "sibling id daughter1 = " << sibling_el_ids[0] << std::endl;
								//std::cout << "current_1d_el_idx  = " << current_1d_el_idx  << std::endl;
								for(unsigned int i=0; i<sibling_el_ids.size(); i++)
									add_to_matrices (system,eqn_3_dof,dof_indices_siblings_q[i][0],-1.0);
							}
			
							// boundary condition should be on side 0
							std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,0);
							if(boundary_ids.size() == 0)
								std::cout << "error, boundary does not have boundary id as expected, el idx = " << current_el_idx << std::endl;

							//std::cout << "flux_values[" << boundary_ids[0] <<"] = " << flux_values[boundary_ids[0]] << std::endl;

							if(!coupled)		//if coupled then the rest is put in the matrix by the 3D assembler
							{
								// if there is a bifurcation at the inflow then daughter_1 can add it
								if(is_daughter_1)
								{
									std::cout << "double oopsie" << std::endl;
									system->rhs->add (eqn_3_dof,-flux_values[boundary_ids[0]]);
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
					// pressure
					else if(es->parameters.get<unsigned int> ("bc_type_1d") == 1)
					{
		    		add_to_matrices (system,eqn_3_dof,dof_indices_p[0],1.0);
		    		system->rhs->add (eqn_3_dof,in_pressure);
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
					//can only be pressure
	    		add_to_matrices (system,eqn_4_dof,dof_indices_p[1],1.0);
	    		system->rhs->add (eqn_4_dof,out_pressure);

				}
				else
				{
					//pressure continuity to daughter_1
	    		add_to_matrices (system,eqn_4_dof,dof_indices_p[1],1.0);
	    		add_to_matrices
							(system,eqn_4_dof,dof_indices_daughter_1_p[0],-1.0);
				}
			
			}

    }// end of subdomain conditional
	}// end of element loop

  // That's it.
	//system->matrix->close();
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
void NavierStokesAssembler::add_to_matrices (TransientLinearImplicitSystem * system, unsigned int row_number, unsigned int col_number, double value)
{
	
	system->matrix->add (row_number,col_number,value);
	if(es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 6 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 7 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 8 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 9 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 10 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 11 || es->parameters.get<unsigned int>("preconditioner_type_3d1d") == 12)
		system->request_matrix("Preconditioner")->add (row_number,col_number,value);
	
}


// The matrix assembly function to be called at each time step to
// prepare for the linear solve. This uses a picard type linearisation
// as opposed to the newton type linearisation used in the 
// assemble_stokes function. Well now assemble_stokes is not working in 3D
void NavierStokesAssembler::assemble_stokes_1D ()
{

	std::cout << "begin assembly" << std::endl;

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es->get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

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
  const unsigned int p_var = system->variable_number ("p");

  // Get the Finite Element type for "p".
  FEType fe_pres_type = system->variable_type(p_var);

  // Build a Finite Element object of the specified type for
  // the velocity variables.

	// element quadrature
  AutoPtr<FEBase> fe_pres  (FEBase::build(dim, fe_pres_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the \p FEType object decide what order rule is appropriate.
  QGauss qrule (dim, fe_pres_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_pres->attach_quadrature_rule (&qrule);

	//face quadrature
	AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_pres_type));
	QGauss qface(dim-1, fe_pres_type.default_quadrature_order());
	fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW = fe_pres->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  //const std::vector<std::vector<Real> >& phi = fe_pres->get_phi();

  // The element shape function gradients for the velocity
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi = fe_pres->get_dphi();
  const std::vector<std::vector<Real> >& phi = fe_pres->get_phi();

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the \p DofMap
  // in future examples.
  const DofMap & dof_map = system->get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

	//JAMES: can create the 3D ones even if working in 2D
  DenseSubMatrix<Number>
    Kpp(Ke);

  DenseSubVector<Number>
    Fp(Fe);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_p;

  // Find out what the timestep size parameter is from the system, and
  // the value of theta for the theta method.  We use implicit Euler (theta=1)
  // for this simulation even though it is only first-order accurate in time.
  // The reason for this decision is that the second-order Crank-Nicolson
  // method is notoriously oscillatory for problems with discontinuous
  // initial data such as the lid-driven cavity.  Therefore,
  // we sacrifice accuracy in time for stability, but since the solution
  // reaches steady state relatively quickly we can afford to take small
  // timesteps.  If you monitor the initial nonlinear residual for this
  // simulation, you should see that it is monotonically decreasing in time.
  //const Real dt    = es->parameters.get<Real>("dt");

	//JAMES: viscosity and unsteadyness
  const Real K    = es->parameters.get<Real>("K");

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the \p active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	// libMesh::processor_id() doesn't exist in libmesh-0.9.4, looks like it does?
	//std::cerr << " proc: " << libMesh::processor_id() << std::endl;
	std::cerr << " proc: " << libMesh::global_processor_id() << std::endl;
								
  for ( ; el != end_el; ++el)
  {				
    const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_p, p_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_p_dofs = dof_indices_p.size();

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe_pres->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      // Reposition the submatrices...  The idea is this:
			// JAMES: note this is now 3D hey
      //
      //         -           -          -  -
      //   Ke = | Kpp |;  Fe = | Fp |
      //         -           -          -  -
      //
      // The \p DenseSubMatrix.repostition () member takes the
      // (row_offset, column_offset, row_size, column_size).
      //
      // Similarly, the \p DenseSubVector.reposition () member
      // takes the (row_offset, row_size)
      Kpp.reposition (p_var*n_p_dofs, p_var*n_p_dofs, n_p_dofs, n_p_dofs);

      Fp.reposition (p_var*n_p_dofs, n_p_dofs);

      // Now we will build the element matrix and right-hand-side.
      // Constructing the RHS requires the solution and its
      // gradient from the previous timestep.  This must be
      // calculated at each quadrature point by summing the
      // solution degree-of-freedom values by the appropriate
      // weight functions.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
      {
        // First, an i-loop over the velocity degrees of freedom.
        // We know that n_p_dofs == n_v_dofs so we can compute contributions
        // for both at the same time.

        for (unsigned int i=0; i<n_p_dofs; i++)
        {
				std::cout << "phi[" << i << "][" << qp << "] = " << phi[i][qp] << std::endl;

          // Matrix contributions for the uu and vv couplings.
          for (unsigned int j=0; j<n_p_dofs; j++)
        	{
						//will need to check that 1D gradient is correct
            Kpp(i,j) += JxW[qp]*(K*(dphi[i][qp]*dphi[j][qp]));  // diffusion term
						std::cout << "gradient = " << dphi[i][qp] << std::endl;
					}

        }

      } // end of the quadrature point qp-loop

			// boundary condition
      {
        // The penalty value.  \f$ \frac{1}{\epsilon} \f$
        const Real penalty = 1.e10;

        // The following loops over the sides of the element.
        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.

				// don't need when using dirichlet on both boundaries
				
        for (unsigned int s=0; s<elem->n_sides(); s++)
				{
				/*
					// matrix and vector stuff on boundary
					const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
					const std::vector<std::vector<Real> >&  dphi_face = fe_face->get_dphi();
	        const std::vector<Real>& JxW_face = fe_face->get_JxW();
					const std::vector<Point >& qface_point = fe_face->get_xyz();
        	fe_face->reinit(elem, s);
        
  
	    		for (unsigned int qp=0; qp<qface.n_points(); qp++)
      		{
						for (unsigned int i=0; i<n_p_dofs; i++)
		        {
    		      // Matrix contributions for the uu and vv couplings.
    		      for (unsigned int j=0; j<n_p_dofs; j++)
    		    	{
            		Kpp(i,j) += -JxW_face[qp]*(K*(dphi[i][qp]*phi_face[j][qp]);  // diffusion term
							}
						}
					}
				
					*/


					// boundary condition	- note that this elem->neighbor thing don't work so well
					// use boundary info object rather

		    	if (mesh.boundary_info->boundary_id (elem, s) == 0 ||
								mesh.boundary_info->boundary_id (elem, s) == 1)
		      {
		        AutoPtr<Elem> side (elem->build_side(s));

		        // Loop over the nodes on the side.
		        for (unsigned int ns=0; ns<side->n_nodes(); ns++)
		        {
							//JAMES: want pressure =1 at one end and 0 at other

							{
								//set all equal zero
								Real p_value = 0.;

								//at inflow apply zero and outflow apply 1
								if(mesh.boundary_info->boundary_id (elem, s) == 0)
								{
									p_value = 0.0;
								}
								else if (mesh.boundary_info->boundary_id (elem, s) == 1)
								{
									p_value = 1.0;
								}

								// Find the node on the element matching this node on
		            // the side.  That defined where in the element matrix
		            // the boundary condition will be applied.
		            for (unsigned int n=0; n<elem->n_nodes(); n++)
								{
		              if (elem->node(n) == side->node(ns))
		              {
		                // Matrix contribution.
		                Kpp(n,n) += penalty;

		                // Right-hand-side contribution.
		                Fp(n) += penalty*p_value;
		              }
								}	
							}
						}
          } // end face node loop
        } // end if (elem->neighbor(side) == NULL)
      } // end boundary condition section

      // If this assembly program were to be used on an adaptive mesh,
      // we would have to apply any hanging node constraint equations
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The \p SparseMatrix::add_matrix()
      // and \p NumericVector::add_vector() members do this for us.
      system->matrix->add_matrix (Ke, dof_indices);
      system->rhs->add_vector    (Fe, dof_indices);
    } // end of element loop
	}

  // That's it.


	std::cout << "end assembly" << std::endl;
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
				unsigned int current_1d_el_idx = current_el_idx -	n_initial_3d_elem;
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
				unsigned int current_1d_el_idx = current_el_idx -	n_initial_3d_elem;
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
	pressure_values = _pressure_values;
	flux_values = _flux_values;

	

	// cause backwards from 3D
	//for(unsigned int i=0; i<flux_values.size(); i++)
	//	flux_values[i] *= -1.0;
}
