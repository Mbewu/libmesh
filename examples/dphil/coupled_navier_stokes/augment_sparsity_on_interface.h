// only works for linear 1d elements!!!!!


#ifndef AUGMENT_SPARSITY_ON_INTERFACE_H
#define AUGMENT_SPARSITY_ON_INTERFACE_H

#include "libmesh/libmesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "airway.h"

using namespace libMesh;

// Convenient typedef for a map for (element id,side id) --> element neighbor id
typedef std::map< dof_id_type, dof_id_type> ElementIdMap;

class AugmentSparsityOnInterface : public DofMap::AugmentSparsityPattern
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
  std::vector<ElementIdMap> interface_elem_node_map;	//maps 3d element connected to surface to boundary id
  std::vector<ElementIdMap> tree_elem_elem_map;				//vector because has for parent, daughters and sibling

	std::vector<std::vector<unsigned int> > boundary_nodes_1d;



	std::vector<unsigned int> subdomains_3d;
	std::vector<unsigned int> subdomains_1d;
 
	int n_initial_3d_elem;

	bool coupled;
	std::vector<Airway> airway_data;

public:

  /**
   * Constructor.
   */
  AugmentSparsityOnInterface(EquationSystems& _es, std::vector<Airway>&  _airway_data, std::vector<unsigned int> _subdomains_3d,std::vector<unsigned int> _subdomains_1d, unsigned int _n_initial_3d_elem=0, bool _coupled=false)
			:  es(&_es), subdomains_3d(_subdomains_3d), subdomains_1d(_subdomains_1d), n_initial_3d_elem(_n_initial_3d_elem), coupled(_coupled), airway_data(_airway_data)
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
