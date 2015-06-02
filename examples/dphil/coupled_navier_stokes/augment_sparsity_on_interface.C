// local includes
#include "augment_sparsity_on_interface.h"

// libMesh includes
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

void AugmentSparsityOnInterface::augment_sparsity_pattern (SparsityPattern::Graph & ,
                                                           std::vector<dof_id_type> & n_nz,
                                                           std::vector<dof_id_type> & n_oz)
{
	std::cout << "baby" << std::endl;
	
  // get a constant reference to the mesh object
  const MeshBase& mesh = es->get_mesh();
	
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


  const unsigned int p_var = system->variable_number ("P");
  const unsigned int q_var = system->variable_number ("Q");
	const bool threed = es->parameters.get<bool>("threed");

  const DofMap & dof_map = system->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_p;
  std::vector<dof_id_type> dof_indices_q;

  // TODO: This won't work with ParallelMesh.
  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  

	// these boundarynodes -2 are okay because only ever used in the case of coupling
	interface_elem_node_map.resize(mesh.boundary_info->n_boundary_ids() - 2);
	tree_elem_elem_map.resize(4);		//for all the things, parent, daughters and sibling
	boundary_nodes_1d.resize(mesh.boundary_info->n_boundary_ids() - 2, std::vector<unsigned int>(4));


  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;
    
		if(coupled)
		{
			if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
		  {
		    for (unsigned char side=0; side<elem->n_sides(); side++)
		    if (elem->neighbor(side) == NULL)
		    {
					std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,side);
					if(boundary_ids.size() > 0)
					{
						int boundary_id = boundary_ids[0];
						if(boundary_id > 0 && boundary_id < 1000)
						{
		        	interface_elem_node_map[boundary_id][elem->id()] = boundary_id;
						}
					}
		    }
		  }
		}

    if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{

			//***** first calculate the ids on the boundary
			const int current_el_idx = elem->id();

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_q, q_var);



			//****** next calculate the coupling between 1d elements
			int current_1d_el_idx = current_el_idx -	n_initial_3d_elem;
			int parent_el_idx = (int)element_data[current_1d_el_idx][0] + n_initial_3d_elem;
			int daughter_1_el_idx = (int)element_data[current_1d_el_idx][1] + n_initial_3d_elem;
			int daughter_2_el_idx = (int)element_data[current_1d_el_idx][2] + n_initial_3d_elem;
			int sibling_el_idx = (int)element_data[current_1d_el_idx][3] + n_initial_3d_elem;
			// need to figure out whether we are on daughter 1 or daughter 2			
			int is_daughter_1 = (int)element_data[current_1d_el_idx][4];	//this is a bool duh!

			//very dirty hack
			if(parent_el_idx < n_initial_3d_elem)
				parent_el_idx = current_el_idx;

			if(daughter_1_el_idx  < n_initial_3d_elem)
				daughter_1_el_idx  = current_el_idx;

			if(daughter_2_el_idx  < n_initial_3d_elem)
				daughter_2_el_idx  = current_el_idx;

			if(sibling_el_idx  < n_initial_3d_elem)
				sibling_el_idx  = current_el_idx;


			if(coupled)
			{
				std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,0);
				if(boundary_ids.size() > 0)
				{
					if(boundary_ids[0] > 0 && boundary_ids[0] < 1000 && is_daughter_1)	// shouldn't get here if looking at side 0 anyway
					{
						boundary_nodes_1d[boundary_ids[0]][0] = dof_indices_p[0];
						boundary_nodes_1d[boundary_ids[0]][1] = dof_indices_p[1];
						boundary_nodes_1d[boundary_ids[0]][2] = dof_indices_q[0];
						boundary_nodes_1d[boundary_ids[0]][3] = dof_indices_q[1];
					}
				}
			}

			tree_elem_elem_map[0][elem->id()] = parent_el_idx;
			tree_elem_elem_map[1][elem->id()] = daughter_1_el_idx;
			tree_elem_elem_map[2][elem->id()] = daughter_2_el_idx;
			tree_elem_elem_map[3][elem->id()] = sibling_el_idx;
    }

	}

  {
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




		// ******* NOW DO THE 1D TREE COUPLINGS
		// these should be fine because we will probably only ever use 1d elements, but maybe not
    const DofMap &dof_map = system->get_dof_map();
    const unsigned int sys_num = system->number();
		std::vector<int> variable_numbers_1d;
    variable_numbers_1d.push_back(system->variable_number("P"));
    variable_numbers_1d.push_back(system->variable_number("Q"));

	  std::vector<dof_id_type> dof_indices;
	  std::vector<dof_id_type> dof_indices_neighbor;

    // Create a set that will have the DoF numbers we need to augment the sparsity pattern
    std::set<dof_id_type> local_coupled_dofs, remote_coupled_dofs;

		for(unsigned int i=0; i < tree_elem_elem_map.size(); i++)
		{

		  std::map<dof_id_type, dof_id_type>::iterator it     = tree_elem_elem_map[i].begin();
		  std::map<dof_id_type, dof_id_type>::iterator it_end = tree_elem_elem_map[i].end();

		  for( ; it != it_end; ++it)
		  {
		    local_coupled_dofs.clear();
		    remote_coupled_dofs.clear();

		    // get a pointer to the one of the elements on the 3d interface
		    const Elem* this_elem = mesh.elem(it->first);
				
				// for some reason, when you try to call the generic dof_indices, or try and get the fod number from the 				
				// node it is not updated!! but the 3d ones are fine!! wtf??
		    const Elem* neighbor_elem = mesh.elem(it->second);

				// for the 4 variables p0,p1,q0,q1
				for(unsigned int vn=0; vn<variable_numbers_1d.size(); vn++)
				{	
	      	dof_map.dof_indices (neighbor_elem, dof_indices,variable_numbers_1d[vn]);
					for(unsigned int dm=0; dm<dof_indices.size(); dm++)
					{
					  const dof_id_type global_dof_number = dof_indices[dm];

					  // and finally insert it into one of the sets
					  if ((global_dof_number <  dof_map.first_dof()) ||
					      (global_dof_number >= dof_map.end_dof()))
					  {
					    remote_coupled_dofs.insert(global_dof_number);
					  }
					  else
					  {
					    local_coupled_dofs.insert(global_dof_number);
						}
					}
				}


				// dunno what this means but stays pretty much the same except the 
				// multiple variables to extend.
				// conservatively increase the preallocation for the implicit matrix contribution
				// to account for these dofs
				for(unsigned int vn=0; vn<variable_numbers_1d.size(); vn++)
				{	
	      	dof_map.dof_indices (this_elem, dof_indices,variable_numbers_1d[vn]);
					for(unsigned int dm=0; dm<dof_indices.size(); dm++)
					{
				  	const dof_id_type
				    	global_dof_number = dof_indices[dm],
				    n_local_dofs      = dof_map.n_local_dofs(),
				    n_remote_dofs     = dof_map.n_dofs() - n_local_dofs;


						// only monkey with the sparsity pattern for local dofs!
						if ((global_dof_number >= dof_map.first_dof()) &&
													(global_dof_number  < dof_map.end_dof()))
						{
							const dof_id_type
								dof_offset = global_dof_number - dof_map.first_dof();

							libmesh_assert_less (dof_offset, n_nz.size());
							libmesh_assert_less (dof_offset, n_oz.size());

							n_nz[dof_offset] += local_coupled_dofs.size();
							n_oz[dof_offset] += remote_coupled_dofs.size();

							// for crazy coarse problems on many processors we need to impose sane limits
							// since we shouldn't allow n_nz > n_local_dofs, for example
							n_nz[dof_offset] = std::min(n_nz[dof_offset], n_local_dofs);
							n_oz[dof_offset] = std::min(n_oz[dof_offset], n_remote_dofs);
						}
					}
				}
      } // end loop over nodes on target element
    }





		// ******* now do the couplings on the interface
		// needs to go both ways!
		if(coupled)
		{
		  std::vector<int> variable_numbers_3d;
		  variable_numbers_3d.push_back(system->variable_number("u"));
		  variable_numbers_3d.push_back(system->variable_number("v"));
			if(threed)
			  variable_numbers_3d.push_back(system->variable_number("w"));
		  variable_numbers_3d.push_back(system->variable_number("p")); // just to be safe...

			// vector containing the dof indices for the different variables
			std::vector<std::vector<dof_id_type> > variable_dof_indices(4);

			for(unsigned int i=1; i < interface_elem_node_map.size(); i++)
			{
				std::map<dof_id_type, dof_id_type>::iterator it     = interface_elem_node_map[i].begin();
				std::map<dof_id_type, dof_id_type>::iterator it_end = interface_elem_node_map[i].end();

				for( ; it != it_end; ++it)
				{
				  local_coupled_dofs.clear();
				  remote_coupled_dofs.clear();

				  // get a pointer to the one of the elements on the 3d interface
				  const Elem* this_elem = mesh.elem(it->first);

					for(unsigned int vn=0; vn<variable_numbers_3d.size(); vn++)
					{	
      			dof_map.dof_indices (this_elem, variable_dof_indices[vn], variable_numbers_3d[vn]);
					}

					// for the 4 variables p0,p1,q0,q1
					for(unsigned int j=0; j < 4; j++)
					{
					  const dof_id_type global_dof_number = boundary_nodes_1d[i][j];

					  // and finally insert it into one of the sets
					  if ((global_dof_number <  dof_map.first_dof()) ||
					      (global_dof_number >= dof_map.end_dof()))
					  {
					    remote_coupled_dofs.insert(global_dof_number);
					  }
					  else
					  {
					    local_coupled_dofs.insert(global_dof_number);
						}
					}

					// dunno what this means but stays pretty much the same except the 
					// multiple variables to extend. but problem is doesn't do it in a symmetric manner
					// conservatively increase the preallocation for the implicit matrix contribution
					// to account for these dofs
		
					for(unsigned int vn=0; vn<variable_numbers_3d.size(); vn++)
					{	
						for(unsigned int di=0; di<variable_dof_indices[vn].size(); di++)
						{
					  	const dof_id_type
					    	global_dof_number = variable_dof_indices[vn][di],
					    n_local_dofs      = dof_map.n_local_dofs(),
					    n_remote_dofs     = dof_map.n_dofs() - n_local_dofs;

							// only monkey with the sparsity pattern for local dofs!
							if ((global_dof_number >= dof_map.first_dof()) &&
														(global_dof_number  < dof_map.end_dof()))
							{
								const dof_id_type
									dof_offset = global_dof_number - dof_map.first_dof();

								libmesh_assert_less (dof_offset, n_nz.size());
								libmesh_assert_less (dof_offset, n_oz.size());

								n_nz[dof_offset] += local_coupled_dofs.size();
								n_oz[dof_offset] += remote_coupled_dofs.size();

								// for crazy coarse problems on many processors we need to impose sane limits
								// since we shouldn't allow n_nz > n_local_dofs, for example
								n_nz[dof_offset] = std::min(n_nz[dof_offset], n_local_dofs);
								n_oz[dof_offset] = std::min(n_oz[dof_offset], n_remote_dofs);
							}
						}
					}

					// now do the other way for the other equation!!!
				  local_coupled_dofs.clear();
				  remote_coupled_dofs.clear();

					for(unsigned int vn=0; vn<variable_numbers_3d.size(); vn++)
					{	
						for(unsigned int di=0; di<variable_dof_indices[vn].size(); di++)
						{
					  	const dof_id_type
					    	global_dof_number = variable_dof_indices[vn][di];

							// and finally insert it into one of the sets
							if ((global_dof_number <  dof_map.first_dof()) ||
							    (global_dof_number >= dof_map.end_dof()))
							{
							  remote_coupled_dofs.insert(global_dof_number);
							}
							else
							{
							  local_coupled_dofs.insert(global_dof_number);
							}
						}
					}

					// dunno what this means but stays pretty much the same except the 
					// multiple variables to extend. but problem is doesn't do it in a symmetric manner
					// conservatively increase the preallocation for the implicit matrix contribution
					// to account for these dofs
		
					// for the 4 variables p0,p1,q0,q1	- only really needs to be for the q0 equation
					for(unsigned int j=0; j < 4; j++)
					{
					  const dof_id_type global_dof_number = boundary_nodes_1d[i][j],
				    n_local_dofs      = dof_map.n_local_dofs(),
				    n_remote_dofs     = dof_map.n_dofs() - n_local_dofs;

						// only monkey with the sparsity pattern for local dofs!
						if ((global_dof_number >= dof_map.first_dof()) &&
													(global_dof_number  < dof_map.end_dof()))
						{
							const dof_id_type
								dof_offset = global_dof_number - dof_map.first_dof();

							libmesh_assert_less (dof_offset, n_nz.size());
							libmesh_assert_less (dof_offset, n_oz.size());

							n_nz[dof_offset] += local_coupled_dofs.size();
							n_oz[dof_offset] += remote_coupled_dofs.size();

							// for crazy coarse problems on many processors we need to impose sane limits
							// since we shouldn't allow n_nz > n_local_dofs, for example
							n_nz[dof_offset] = std::min(n_nz[dof_offset], n_local_dofs);
							n_oz[dof_offset] = std::min(n_oz[dof_offset], n_remote_dofs);
						}
					}

		    } // end loop over nodes on target element
		  }
		}


  }
	      
}
