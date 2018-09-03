// only works for linear 1d elements!!!!!


#ifndef AUGMENT_SPARSITY_ACINAR_H
#define AUGMENT_SPARSITY_ACINAR_H

#include "libmesh/libmesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"

// include for preconditioner getting pressure dofs on side
#include "libmesh/fe_interface.h"
#include "airway.h"

using namespace libMesh;

// Convenient typedef for a map for (element id,side id) --> element neighbor id
typedef std::map< dof_id_type, dof_id_type> ElementIdMap;

class AugmentSparsityAcinar : public DofMap::AugmentSparsityPattern
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
  std::vector<std::vector<unsigned int> > dofs_1d;	// the dofs that need to be coupled for each acinus
  std::vector<std::vector<unsigned int> > dofs_acinar;	// the dofs that need to be coupled for each acinus

	std::vector<unsigned int> subdomains_1d;
	std::vector<unsigned int> subdomains_acinar;

	std::vector<unsigned int> elem_to_airway;
//	int n_initial_3d_elem;

	bool coupled;
	std::vector<Airway> airway_data;
	std::vector<unsigned int> elem_to_acinar;
	std::vector<unsigned int> acinar_to_airway;

public:

  /**
   * Constructor.
   */
  AugmentSparsityAcinar(EquationSystems& _es, std::vector<Airway>&  _airway_data, std::vector<unsigned int> _subdomains_1d, std::vector<unsigned int> _subdomains_acinar,std::vector<unsigned int> _elem_to_acinar,std::vector<unsigned int> _acinar_to_airway, std::vector<unsigned int> _elem_to_airway, bool _coupled=false)
			:  es(&_es), subdomains_1d(_subdomains_1d), subdomains_acinar(_subdomains_acinar), elem_to_airway(_elem_to_airway), coupled(_coupled), airway_data(_airway_data), elem_to_acinar(_elem_to_acinar), acinar_to_airway(_acinar_to_airway)
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
