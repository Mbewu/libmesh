// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "libmesh/petsc_auto_fieldsplit.h"

#ifdef LIBMESH_HAVE_PETSC

EXTERN_C_FOR_PETSC_BEGIN
#  include <petscksp.h>
EXTERN_C_FOR_PETSC_END

// Local includes
#include "libmesh/dof_map.h"
#include "libmesh/system.h"
#include "libmesh/parallel_implementation.h"

#if !PETSC_VERSION_LESS_THAN(3,2,0)

// C++ includes

namespace {
using namespace libMesh;

void indices_to_fieldsplit (const Parallel::Communicator& comm,
                            const std::vector<dof_id_type>& indices,
                            PC my_pc,
                            const std::string& field_name)
{

	std::cerr << "giving the indices to fieldsplit" << std::endl;
  const PetscInt *idx = PETSC_NULL;
  if (!indices.empty())
    idx = reinterpret_cast<const PetscInt*>(&indices[0]);

	std::cerr << "fieldname = " << field_name << std::endl;

  IS is;
  int ierr = ISCreateLibMesh(comm.get(), indices.size(),
                             idx, PETSC_COPY_VALUES, &is);
  CHKERRABORT(comm.get(), ierr);

  ierr = PCFieldSplitSetIS(my_pc, field_name.c_str(), is);
  CHKERRABORT(comm.get(), ierr);

	PCType type;
	ierr = PCGetType(my_pc,&type);
  CHKERRABORT(comm.get(), ierr);

	//std::cerr << "pc type = " << type << std::endl;

}

}

namespace libMesh
{

void petsc_auto_fieldsplit (PC my_pc,
                            const System &sys)
{

	std::cerr << "RANK " << sys.comm().rank() << ": Setting up the IS for fieldsplit." << std::endl;

  std::string sys_prefix = "--solver_group_";

  if (libMesh::on_command_line("--solver_system_names"))
    {
      sys_prefix = sys_prefix + sys.name() + "_";
    }

  std::map<std::string, std::vector<dof_id_type> > group_indices;

  if (libMesh::on_command_line("--solver_variable_names"))
    {
      for (unsigned int v = 0; v != sys.n_vars(); ++v)
        {
          const std::string& var_name = sys.variable_name(v);

          std::vector<dof_id_type> var_idx;
          sys.get_dof_map().local_variable_indices
            (var_idx, sys.get_mesh(), v);

	  // the outer fieldsplit grouping
          std::string group_command = sys_prefix + var_name;

          const std::string empty_string;

					// JAMES EDIT:
					// changed so that actually gets the command line bit
          std::string group_name = libMesh::command_line_next
            (group_command, empty_string);

			Mat pc_mat;
			PCGetOperators(my_pc,&pc_mat,NULL);
			PetscInt m,n;
			MatGetLocalSize(pc_mat,&m,&n);
			std::cerr << "RANK " << sys.comm().rank() << ": rows = " << m << ", cols = " << n << std::endl;
			PetscInt r,s;
			MatGetOwnershipRange(pc_mat,&r,&s);
			std::cerr << "RANK " << sys.comm().rank() << ": subrange = " << r << " to " << s << std::endl;

	std::cerr << "RANK " << sys.comm().rank() << ": yeah doing " << var_name << " figuring out now" << std::endl; 

          if (group_name != empty_string)
            {
							//std::cerr << "RANK " << sys.comm().rank() << ": yeah doing " << var_name << " u figuring out now" << std::endl; 
              std::vector<dof_id_type> &indices =
                group_indices[group_name];
              const bool prior_indices = !indices.empty();
              indices.insert(indices.end(), var_idx.begin(),
                             var_idx.end());
              if (prior_indices)
                std::sort(indices.begin(), indices.end());
            }
          else
            {
							//std::cerr << "RANK " << sys.comm().rank() << ": doing outer  now" << std::endl; 
              indices_to_fieldsplit (sys.comm(), var_idx, my_pc, var_name);
            }
        }
    }
	
	// want to do groups first
  for (std::map<std::string, std::vector<dof_id_type> >::const_iterator
         i = group_indices.begin(); i != group_indices.end(); ++i)
    {
			//std::cerr << "RANK " << sys.comm().rank() << ": k doing outer group" << std::endl;
      indices_to_fieldsplit(sys.comm(), i->second, my_pc, i->first);
    }



  if(true)
  {

	std::cerr << "RANK " << sys.comm().rank() << ": Setting up the IS for the sub fieldsplits." << std::endl;
	// do this so that they know what preconditioner to use for the sub systems
	PetscErrorCode ierr=0;
	ierr = PCSetUp(my_pc);


	// now we want to do the fieldsplit of the sub fieldsplit


	PCType pc_type;
	ierr = PCGetType(my_pc,&pc_type);
	std::cerr << "RANK " << sys.comm().rank() << ": outer pc type = " << pc_type << std::endl;
	std::string fieldsplit_string = "fieldsplit";
	if (pc_type == fieldsplit_string)
	{
		std::cerr << "RANK " << sys.comm().rank() << ": yeah..." << std::endl;
		// loop over the number of dofs, cause this is the possible amount of fieldsplits
		KSP* subksp;	//subksps
		PetscInt num_splits;
		ierr = PCFieldSplitGetSubKSP(my_pc,&num_splits,&subksp);// CHKERRQ(ierr);

		

		for (unsigned int i = 0; i < num_splits; ++i)
		{
			std::cerr << "RANK " << sys.comm().rank() << ": i = " << i << std::endl;
			PC sub_pc;
			ierr = KSPGetPC(subksp[i],&sub_pc);// CHKERRQ(ierr);

			

			PCType sub_pc_type;
			ierr = PCGetType(sub_pc,&sub_pc_type);
			std::cerr << "RANK " << sys.comm().rank() << ": sub_pc type = " << sub_pc_type << std::endl;
			Mat sub_mat;
			ierr = PCGetOperators(sub_pc,&sub_mat,NULL);
			PetscInt m,n;
			ierr = MatGetLocalSize(sub_mat,&m,&n);
			std::cerr << "RANK " << sys.comm().rank() << ": subrows = " << m << ", subcols = " << n << std::endl;
			PetscInt r,s;
			MatGetOwnershipRange(sub_mat,&r,&s);
			std::cerr << "RANK " << sys.comm().rank() << ": subrange = " << r << " to " << s << std::endl;
			std::cerr << "RANK " << sys.comm().rank() << ": hmmm" << std::endl;


			if (sub_pc_type == fieldsplit_string)
			{
				// need to extract the pc corresponding to this 
				sys_prefix = sys_prefix + static_cast<std::ostringstream*>( &(std::ostringstream() << i) )->str() + "_";

				std::map<std::string, std::vector<dof_id_type> > group_indices;

				if (libMesh::on_command_line("--solver_variable_names"))
				{
					for (unsigned int v = 0; v != sys.n_vars(); ++v)
					{
						const std::string& var_name = sys.variable_name(v);

						std::vector<dof_id_type> var_idx;
						sys.get_dof_map().local_variable_indices
						(var_idx, sys.get_mesh(), v);

						// the outer fieldsplit grouping
						std::string group_command = sys_prefix + var_name;

						const std::string empty_string;

									// JAMES EDIT:
									// changed so that actually gets the command line bit
						std::string group_name = libMesh::command_line_next
						(group_command, empty_string);


						if (group_name != empty_string)
						{
							std::cerr << "RANK " << sys.comm().rank() << ": yeah doing " << var_name << " figuring out now" << std::endl; 
							std::vector<dof_id_type> &indices =
							group_indices[group_name];
							const bool prior_indices = !indices.empty();
							indices.insert(indices.end(), var_idx.begin(),
								     var_idx.end());

							if (prior_indices)
							std::sort(indices.begin(), indices.end());
						}
					}
				  }

				// need to create a local to global vector for this fieldsplit. Need to first collect all the indices from U and p.

				std::vector<unsigned int> global_indices;
				  for (std::map<std::string, std::vector<dof_id_type> >::const_iterator
					 j = group_indices.begin(); j != group_indices.end(); ++j)
				  {

					global_indices.insert(global_indices.end(), j->second.begin(), j->second.end());
				}

					//need to put all the global indices from all the processes into a single vector so we can figure out the mapping
					sys.comm().allgather(global_indices,false);

					// need to sort because pressure dofs come after
					std::sort(global_indices.begin(), global_indices.end());

					// want to do groups first
				  for (std::map<std::string, std::vector<dof_id_type> >::const_iterator
					 j = group_indices.begin(); j != group_indices.end(); ++j)
				  {

					//need to put all the global indices from all the processes into a single vector so we can figure out the mapping
					std::vector<unsigned int> local_indices = j->second;

					// loop over the sub indices and exchange them for its local index
					// should be sorted so we can use the previous idx as a starting point
					int find_start = 0;
					for(unsigned int k=0; k<j->second.size(); k++)
					{
						int local_idx = find(global_indices.begin() + find_start,global_indices.end(),j->second[k]) - global_indices.begin();
						local_indices[k] = local_idx;
						find_start = local_idx;
					}


				      indices_to_fieldsplit(sys.comm(), local_indices, sub_pc, j->first);
				  }

				//std::exit(0);
				std::cerr << "RANK " << sys.comm().rank() << ": before ksp set up" << std::endl;
				ierr = KSPSetUp(subksp[i]);// CHKERRQ(ierr);
				std::cerr << "RANK " << sys.comm().rank() << ": before pc set up" << std::endl;
				ierr = PCSetUp(sub_pc);
				std::cerr << "RANK " << sys.comm().rank() << ": done" << std::endl;

			}


			

			
		}
	}
  }

	std::cerr << "RANK " << sys.comm().rank() << ": done fieldsplit setup" << std::endl;

	//int ierr = PCSetDiagonalScale(my_pc,sys.solution->vec());
  //CHKERRABORT(sys.comm().get(), ierr);
}

} // namespace libMesh


#else  // #PETSC_VERSION < 3.2.0

namespace libMesh
{
void petsc_auto_fieldsplit (PC /* my_pc */,
                            const System & /* sys */)
{
  if (libMesh::on_command_line("--solver_variable_names"))
    {
      libmesh_do_once(
                      libMesh::out << "WARNING: libMesh does not support setting field splits" <<
                      std::endl << "with PETSc "
                      << LIBMESH_DETECTED_PETSC_VERSION_MAJOR << '.'
                      << LIBMESH_DETECTED_PETSC_VERSION_MINOR << '.'
                      << LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR << std::endl;);
    }
}
}

#endif // #PETSC_VERSION > 3.2.0
#endif // #ifdef LIBMESH_HAVE_PETSC
