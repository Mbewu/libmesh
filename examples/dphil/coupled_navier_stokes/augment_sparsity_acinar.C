// local includes
#include "augment_sparsity_acinar.h"

// libMesh includes
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

// we augment the sparsity pattern by coupling the acinar dofs to the terminal 0D dofs
void AugmentSparsityAcinar::augment_sparsity_pattern (SparsityPattern::Graph & ,
                                                           std::vector<dof_id_type> & n_nz,
                                                           std::vector<dof_id_type> & n_oz)
{
	std::cout << "Augmenting acinar sparsity pattern..." << std::endl;
	
  // get a constant reference to the mesh object
  const MeshBase& mesh = es->get_mesh();
	
	// moghadam only works with ns3d
	TransientLinearImplicitSystem * system;
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

	const bool threed = es->parameters.get<bool>("threed");

  const unsigned int dim = mesh.mesh_dimension ();
  const unsigned int p_var = system->variable_number ("P");
  const unsigned int q_var = system->variable_number ("Q");
  const unsigned int v_var = system->variable_number ("V");


  const DofMap & dof_map = system->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  //std::vector < dof_id_type > dof_indices_p;
  //std::vector < dof_id_type > dof_indices_q;
  //std::vector < dof_id_type > dof_indices_v;

  // TODO: This won't work with ParallelMesh.
  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  

	// need to resize to the number of terminals
	dofs_1d.resize(acinar_to_airway.size());
	dofs_acinar.resize(acinar_to_airway.size());

	std::cout << "before loop" << std::endl;
  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;
    
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
	  {
			const int current_el_idx = elem->id();
			unsigned int current_1d_el_idx = elem_to_airway[current_el_idx]; //current_el_idx -	n_initial_3d_elem;
			int acinar_id = airway_data[current_1d_el_idx].get_acinar_id();
			// if the element has an acinar, then add its dofs
			if(acinar_id > -1)
			{
//				std::cout << "acinar_id = " << acinar_id << std::endl;
//				std::cout << "current_1d_el_idx = " << current_1d_el_idx << std::endl;
				// get dof indices for element
				dof_map.dof_indices (elem, dof_indices);

				// add them all to the array
				for(unsigned int i=0; i<dof_indices.size(); i++)
				{
					dofs_1d[acinar_id].push_back(dof_indices[i]);
//					std::cout << " 1d dof added = " << dof_indices[i] << std::endl;
				}
			}
		}


		if(std::find(subdomains_acinar.begin(), subdomains_acinar.end(), elem->subdomain_id()) != subdomains_acinar.end())
	  {
			const int current_el_idx = elem->id();
			// get dof indices for element
			dof_map.dof_indices (elem, dof_indices);

//			std::cout << "acinar_id = " << elem_to_acinar[current_el_idx] << std::endl;

			// add them all to the array
			for(unsigned int i=0; i<dof_indices.size(); i++)
			{
				dofs_acinar[elem_to_acinar[current_el_idx]].push_back(dof_indices[i]);
//				std::cout << " acinar dof added = " << dof_indices[i] << std::endl;
			}
		}
	}


	std::cout << "after loop" << std::endl;


  // Create a set that will have the DoF numbers we need to augment the sparsity pattern
  std::set<dof_id_type> local_coupled_dofs, remote_coupled_dofs;

	// loop over each acinar
	for(unsigned int i=0; i < dofs_acinar.size(); i++)
	{

		// loop over each dof on this acinus (should be one at the moment, but you never know in future)
		for(unsigned int j=0; j<dofs_acinar[i].size(); j++)
		{

			local_coupled_dofs.clear();
			remote_coupled_dofs.clear();

			// we couple all the boundary element dofs to this element
			for(unsigned int k=0; k< dofs_1d[i].size(); k++)
			{
				local_coupled_dofs.insert(dofs_1d[i][k]);
				remote_coupled_dofs.insert(dofs_1d[i][k]);
			}

			// the dof we are coupling all the boundary dofs to
			const dof_id_type
		  	global_dof_number = dofs_acinar[i][j],
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

	std::cout << "to me" << std::endl;
	      
}


