// local includes
#include "augment_sparsity_moghadam.h"

// libMesh includes
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"

void AugmentSparsityMoghadam::augment_sparsity_pattern (SparsityPattern::Graph & ,
                                                           std::vector<dof_id_type> & n_nz,
                                                           std::vector<dof_id_type> & n_oz)
{
	std::cout << "baby don lie" << std::endl;
	
  // get a constant reference to the mesh object
  const MeshBase& mesh = es->get_mesh();
	
	// moghadam only works with ns3d
	TransientLinearImplicitSystem * system;
	system =
		&es->get_system<TransientLinearImplicitSystem> ("ns3d");

	const bool threed = es->parameters.get<bool>("threed");

  const unsigned int dim = mesh.mesh_dimension ();
  const unsigned int u_var = system->variable_number ("u");
  const unsigned int v_var = system->variable_number ("v");
  unsigned int w_var = 0;
  if (threed)
    w_var = system->variable_number ("w");


  const DofMap & dof_map = system->get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector < dof_id_type > dof_indices_u;
  std::vector < dof_id_type > dof_indices_v;
  std::vector < dof_id_type > dof_indices_w;

  // Finite element types
  FEType fe_vel_type = system->variable_type (u_var);
  AutoPtr < FEBase > fe_vel (FEBase::build (dim, fe_vel_type));


  // TODO: This won't work with ParallelMesh.
  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  

	// these boundarynodes -2 are okay because only ever used in the case of coupling
	boundary_elem_dofs.resize(mesh.boundary_info->n_boundary_ids() - 2);

	std::cout << "before loop" << std::endl;
  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;
    
		if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
	  {
			// set dof indices for element
			dof_map.dof_indices (elem, dof_indices);
			dof_map.dof_indices (elem, dof_indices_u, u_var);
			dof_map.dof_indices (elem, dof_indices_v, v_var);
			if (threed)
			  dof_map.dof_indices (elem, dof_indices_w, w_var);

	    for (unsigned char side=0; side<elem->n_sides(); side++)
			{
			  if (elem->neighbor(side) == NULL)
			  {
					std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(elem,side);
					if(boundary_ids.size() > 0)
					{

						int boundary_id = boundary_ids[0];
						if(boundary_id > 0 && boundary_id < 1000)
						{
							std::vector <unsigned int >velocity_dofs_on_side;
							FEInterface::dofs_on_side (elem, dim,fe_vel_type, side,velocity_dofs_on_side);

							// need the dofs and the element
							// first the dofs
							for (unsigned int i = 0; i < velocity_dofs_on_side.size (); i++)
							{
								boundary_elem_dofs[boundary_id].push_back (dof_indices_u[velocity_dofs_on_side[i]]);
								boundary_elem_dofs[boundary_id].push_back (dof_indices_v[velocity_dofs_on_side[i]]);
								if(threed)
									boundary_elem_dofs[boundary_id].push_back (dof_indices_w[velocity_dofs_on_side[i]]);
							}
						}
					}
			  }
			}
		}
	}

	// sort and remove duplicates of dofs
	for(unsigned int i=0; i<boundary_elem_dofs.size(); i++)
	{
		std::sort (boundary_elem_dofs[i].begin (),boundary_elem_dofs[i].end ());
		// remove duplicate dofs
		std::vector < unsigned int >::iterator it;
		it = std::unique (boundary_elem_dofs[i].begin (),boundary_elem_dofs[i].end ());
		boundary_elem_dofs[i].resize (std::distance (boundary_elem_dofs[i].begin (), it));
	}


	std::cout << "after loop" << std::endl;


  // Create a set that will have the DoF numbers we need to augment the sparsity pattern
  std::set<dof_id_type> local_coupled_dofs, remote_coupled_dofs;

	// loop over each boundary
	for(unsigned int i=1; i < boundary_elem_dofs.size(); i++)
	{

		// loop over each dof on this boundary
		for(unsigned int j=0; j<boundary_elem_dofs[i].size(); j++)
		{

			local_coupled_dofs.clear();
			remote_coupled_dofs.clear();

			// we couple all the boundary element dofs to this element
			for(unsigned int k=0; k< boundary_elem_dofs[i].size(); k++)
			{
				local_coupled_dofs.insert(boundary_elem_dofs[i][k]);
				remote_coupled_dofs.insert(boundary_elem_dofs[i][k]);
			}

			// the dof we are coupling all the boundary dofs to
			const dof_id_type
		  	global_dof_number = boundary_elem_dofs[i][j],
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


