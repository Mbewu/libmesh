
#include "airway.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

Airway::Airway ()

{

	// set defaults
	generation = -1;
	is_daughter_1 = true;
	primary_sibling = -1;
	length = -1;
	radius = -1;
	order = -1;
	flow_rate = 0.;
	efficiency = 0.;
	num_alveolar_generations = 0.;
	parent = -1;
	node_1 = Point(0.,0.,0.);
	node_2 = Point(0.,0.,0.);
	poiseuille = true;
	local_elem_number = 0;
	tree_number = 0;
	flow = 0.;
	pressure_diff = 0.;

	previous_flow_rate = 0.;
	velocity = 0.;
	previous_velocity = 0.;

	terminal_exit_fraction = 0.;

	particle_fraction = 0.;

}

// all the get methods

int Airway::get_generation() 
{ 
	if(generation < 0)
	{
		std::cout << "error, generation has not been set. EXITING." <<  std::endl;
		std::cout << "airway data:" << std::endl;
		this->print_concise();
		std::exit(0);
	}
	else
		return generation; 
};

int Airway::get_primary_sibling() 
{ 
	if(primary_sibling < 0)
	{
		std::cout << "error, primary_sibling has not been set. EXITING."<< std::endl;
		std::cout << "airway data:" << std::endl;
		this->print_concise();
		std::exit(0);
	}
	else
		return primary_sibling; 
};

double Airway::get_length() 
{ 
	if(length < 0)
	{
		std::cout << "error, length has not been set (or has been set to be negative). EXITING." << std::endl;
		std::cout << "airway data:" << std::endl;
		this->print_concise();
		std::exit(0);
	}
	else
		return length; 
};

double Airway::get_radius() 
{ 
	if(radius < 0)
	{
		std::cout << "error, radius has not been set (or has been set to be negative). EXITING." << std::endl;
		std::cout << "airway data:" << std::endl;
		this->print_concise();
		std::exit(0);
	}
	else
		return radius; 
};

int Airway::get_order() 
{ 
	if(order < 0)
	{
		std::cout << "error, order has not been set. EXITING." << std::endl;
		std::cout << "airway data:" << std::endl;
		this->print_concise();
		std::exit(0);
	}
	else
		return order; 
};


int Airway::get_daughter_1() 
{ 
	if(daughters.size() == 0)
	{
		std::cout << "error, daughter_1 has not been set. EXITING." << std::endl;
		std::cout << "airway data:" << std::endl;
		this->print_concise();
		std::exit(0);
	}
	else
		return daughters[0]; 
};

int Airway::get_parent() 
{ 
	if(parent < 0)
	{
		std::cout << "error, parent has not been set. EXITING." << std::endl;
		std::cout << "airway data:" << std::endl;
		this->print_concise();
		std::exit(0);
	}
	else
		return parent; 
};


bool Airway::has_daughter_1() 
{ 
	if(daughters.size() > 0)
	{
		return true;
	}
	else
		return false; 
};


bool Airway::has_parent() 
{ 
	if(parent < 0)
	{
		return false;
	}
	else
		return true; 
};

bool Airway::has_primary_sibling() 
{ 
	if(primary_sibling < 0)
	{
		return false;
	}
	else
		return true; 
};


void Airway::move_element_numbers(unsigned int amount)
{
	// we need to move all of the element numbers
	for(unsigned int i=0; i<siblings.size(); i++)
		siblings[i] += amount;

	for(unsigned int i=0; i<daughters.size(); i++)
		daughters[i] += amount;

	if(primary_sibling >= 0)
		primary_sibling += amount;

	if(parent >= 0)
		parent += amount;
}

void Airway::print()
{
	std::cout << "\t local_elem_number = " << local_elem_number << std::endl;
	std::cout << "\t generation = " << generation << std::endl;
	std::cout << "\t is_daughter_1 = " << is_daughter_1 << std::endl;
	std::cout << "\t num_siblings = " << siblings.size() << std::endl;
	for(unsigned int i=0; i<siblings.size(); i++)
	{
		std::cout << "\t\t " << i << " = " << siblings[i] << std::endl;
	}
	std::cout << "\t primary_sibling = " << primary_sibling << std::endl;
	std::cout << "\t length = " << length << std::endl;
	std::cout << "\t radius = " << radius << std::endl;
	std::cout << "\t order = " << order << std::endl;
	std::cout << "\t flow_rate = " << flow_rate << std::endl;
	std::cout << "\t efficiency = " << efficiency << std::endl;
	std::cout << "\t num_alveolar_generations = " << num_alveolar_generations << std::endl;
	std::cout << "\t num_daughters = " << daughters.size() << std::endl;
	for(unsigned int i=0; i<daughters.size(); i++)
	{
		std::cout << "\t\t " << i << " = " << daughters[i] << std::endl;
	}
	std::cout << "\t parent = " << parent << std::endl;
	std::cout << "\t node1 = " << node_1 << std::endl;
	std::cout << "\t node2 = " << node_2 << std::endl;
	std::cout << "\t poiseuille = " << poiseuille << std::endl;

}


void Airway::print_concise()
{
	std::cout << "\t local_elem_number = " << local_elem_number << std::endl;
	std::cout << "\t generation = " << generation << std::endl;
	std::cout << "\t parent = " << parent << std::endl;
	std::cout << "\t is_daughter_1 = " << is_daughter_1 << std::endl;
	std::cout << "\t num_siblings = " << siblings.size() << std::endl;
	for(unsigned int i=0; i<siblings.size(); i++)
	{
		std::cout << "\t\t " << i << " = " << siblings[i] << std::endl;
	}
	std::cout << "\t primary_sibling = " << primary_sibling << std::endl;
	std::cout << "\t length = " << length << std::endl;
	std::cout << "\t radius = " << radius << std::endl;
	std::cout << "\t num_daughters = " << daughters.size() << std::endl;
	for(unsigned int i=0; i<daughters.size(); i++)
	{
		std::cout << "\t\t " << i << " = " << daughters[i] << std::endl;
	}

}



void Airway::add_particle_fraction(double amount, double entry_time)
{
	particle_fraction_amount.push_back(amount);
	particle_fraction_entry_time.push_back(entry_time);
	particle_fraction_distance.push_back(0.);
}


unsigned int Airway::num_particle_fractions()
{
	return particle_fraction_amount.size();
}

double Airway::get_total_particle_fraction()
{
	double total_particle_fraction = 0.;
	
	// we assume that if this is called during an unsteady simulation
	// Airway.particle_fraction has not been set and is zero as it should be
	if(particle_fraction_amount.size() > 0)
	{
		for(unsigned int i=0; i<particle_fraction_amount.size(); i++)
		{
			total_particle_fraction += particle_fraction_amount[i];
		}
	}
	else
	{
		total_particle_fraction = particle_fraction;
	}

	return total_particle_fraction;
}

void Airway::set_total_particle_fraction(double _particle_fraction)
{
	particle_fraction = _particle_fraction;
}


void Airway::remove_particle_fraction(unsigned int particle_fraction_number)
{
	particle_fraction_amount.erase(particle_fraction_amount.begin() + particle_fraction_number);
	particle_fraction_entry_time.erase(particle_fraction_entry_time.begin() + particle_fraction_number);
	particle_fraction_distance.erase(particle_fraction_distance.begin() + particle_fraction_number);
}
