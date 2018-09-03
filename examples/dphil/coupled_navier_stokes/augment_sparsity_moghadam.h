// only works for linear 1d elements!!!!!


#ifndef AUGMENT_SPARSITY_MOGHADAM_H
#define AUGMENT_SPARSITY_MOGHADAM_H

#include "libmesh/libmesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"

// include for preconditioner getting pressure dofs on side
#include "libmesh/fe_interface.h"

using namespace libMesh;

// Convenient typedef for a map for (element id,side id) --> element neighbor id
typedef std::map< dof_id_type, dof_id_type> ElementIdMap;

class AugmentSparsityMoghadam : public DofMap::AugmentSparsityPattern
{
private:

  /**
   * The EquationSystems object that we're using here.
   */
  EquationSystems* es;

  /**
   * A map from (lower element ID, side ID) to matching upper element ID. Here "lower"
   * and "upper" refer to the lower and upper (wrt +z direction) sides of the crack in
   * our mesh.
   */
  std::vector<std::vector<unsigned int> > boundary_elem_dofs;	// boundary elements on each outlet

	std::vector<unsigned int> subdomains_3d;


public:

  /**
   * Constructor.
   */
  AugmentSparsityMoghadam(EquationSystems& _es, std::vector<unsigned int> _subdomains_3d)
			:  es(&_es), subdomains_3d(_subdomains_3d)
	{ 
		
};
  
  /**
   * User-defined function to augment the sparsity pattern.
   */
  virtual void augment_sparsity_pattern (SparsityPattern::Graph & ,
                                         std::vector<dof_id_type> & n_nz,
                                         std::vector<dof_id_type> & n_oz);

};

#endif
