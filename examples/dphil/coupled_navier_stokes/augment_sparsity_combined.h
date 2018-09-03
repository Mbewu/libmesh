// only works for linear 1d elements!!!!!


#ifndef AUGMENT_SPARSITY_COMBINED_H
#define AUGMENT_SPARSITY_COMBINED_H

#include "libmesh/libmesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"

// include for preconditioner getting pressure dofs on side
#include "libmesh/fe_interface.h"
#include "airway.h"

#include "augment_sparsity_on_interface.h"
#include "augment_sparsity_moghadam.h"
#include "augment_sparsity_acinar.h"

using namespace libMesh;

// Convenient typedef for a map for (element id,side id) --> element neighbor id
typedef std::map< dof_id_type, dof_id_type> ElementIdMap;

class AugmentSparsityCombined : public DofMap::AugmentSparsityPattern
{
private:

	// it's actually monolithic and 1D
	bool monolithic;
	bool moghadam;
	bool acinar;
	
	AutoPtr<AugmentSparsityOnInterface> augment_sparsity_monolithic;
	AutoPtr<AugmentSparsityMoghadam> augment_sparsity_moghadam;
	AutoPtr<AugmentSparsityAcinar> augment_sparsity_acinar;

public:

  /**
   * Constructor.
   */
  AugmentSparsityCombined()
			:  monolithic(false), moghadam(false), acinar(false)
	{ 
		
	};
  
  /**
   * call the augment_sparsity_pattern functions of the relevant types
   */
  virtual void augment_sparsity_pattern (SparsityPattern::Graph & graph,
                                         std::vector<dof_id_type> & n_nz,
                                         std::vector<dof_id_type> & n_oz)
	{
		if(monolithic)
		{
			std::cout << "Augmenting 0D/monolithic sparsity pattern." << std::endl;
			augment_sparsity_monolithic->augment_sparsity_pattern(graph, n_nz, n_oz);
		}

		if(moghadam)
		{
			std::cout << "Augmenting moghadam sparsity pattern." << std::endl;
			augment_sparsity_moghadam->augment_sparsity_pattern(graph, n_nz, n_oz);
		}

		if(acinar)
		{
			std::cout << "Augmenting acinar sparsity pattern." << std::endl;
			augment_sparsity_acinar->augment_sparsity_pattern(graph, n_nz, n_oz);
		}
	}


  /**
   * create the monolithic sparsity pattern
   */
  void create_monolithic_sparsity_object (EquationSystems& es, 
																									std::vector<Airway>&  airway_data, 
																									std::vector<unsigned int>& subdomains_3d,
																									std::vector<unsigned int>& subdomains_1d, 
																									std::vector<unsigned int>& elem_to_airway, 
																									bool coupled=false)
	{
		augment_sparsity_monolithic = AutoPtr<AugmentSparsityOnInterface>
										(new AugmentSparsityOnInterface(es,airway_data,subdomains_3d,subdomains_1d,elem_to_airway,coupled));
		monolithic = true;
	}

  /**
   * create the moghadam sparsity pattern
   */
  void create_moghadam_sparsity_object (EquationSystems& es, 
																				std::vector<unsigned int>& subdomains_3d)
	{
		augment_sparsity_moghadam = AutoPtr<AugmentSparsityMoghadam>
										(new AugmentSparsityMoghadam(es,subdomains_3d));
		moghadam = true;
	}

  /**
   * create the moghadam sparsity pattern
   */
  void create_acinar_sparsity_object (EquationSystems& es, 
																			std::vector<Airway>&  airway_data, 
																			std::vector<unsigned int>& subdomains_1d, 
																			std::vector<unsigned int>& subdomains_acinar,
																			std::vector<unsigned int>& elem_to_acinar,
																			std::vector<unsigned int>& acinar_to_airway, 
																			std::vector<unsigned int>& elem_to_airway, 
																			bool coupled=false)
	{
		augment_sparsity_acinar = AutoPtr<AugmentSparsityAcinar>
											(new AugmentSparsityAcinar(es,airway_data,subdomains_1d,subdomains_acinar,elem_to_acinar,acinar_to_airway,elem_to_airway,coupled));
		acinar = true;
	}


};

#endif
