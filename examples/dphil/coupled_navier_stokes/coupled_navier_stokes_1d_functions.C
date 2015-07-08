#include "coupled_navier_stokes.h"

// ********************************************************** //
//
// This file includes all the functions specific to 
// 1D domains.
//
// setup_1d_mesh
// read_1d_mesh
// calculate_num_alveloar_generations_for_tree
// generate_1d_mesh
// create_1d_tree
// add_1d_tree_to_mesh
// setup_1d_system
// calculate_1d_boundary_values
// write_1d_solution
// solve_1d_system
// set_radii
// convert_1d_monomial_to_nodal
// convert_1d_nodal_to_monomial
// scale_1d_solution_vector
// output_poiseuille_resistance_per_generation
// write_elem_pid_1d
// write_efficiency_solution
// set_efficiency
//
// *********************************************************** //


void NavierStokesCoupled::setup_1d_mesh ()
{

	std::cout << "ahem ahem ahem" << std::endl;
	if(es->parameters.get<bool>("_read_1d_mesh"))
		read_1d_mesh();
	else
		generate_1d_mesh();
}

void NavierStokesCoupled::read_1d_mesh ()
{
	std::cout << "setting up 1d mesh from file" << std::endl;	

	unsigned int num_1d_elements = 0;	

	pressure_values_1d.push_back(0.0);	//inflow
	flux_values_1d.push_back(0.0);	//inflow

	std::cout << "calculate_1d_info_at_coupling_nodes = " << es->parameters.get<bool> ("calculate_1d_info_at_coupling_nodes") << std::endl;
	std::cout << "num_1d_trees = " << es->parameters.get<unsigned int> ("num_1d_trees") << std::endl;


	// number of 1d meshes is all the boundaries minus the inflow boundary
	if(es->parameters.get<bool> ("match_1d_mesh_to_3d_mesh"))
	{
		es->parameters.set<unsigned int> ("num_1d_trees") = surface_boundaries.size() - 1;

		// the inflow boundary id 0 does not have a 1d subdomain
		boundary_id_to_tree_id.resize(surface_boundaries.size(),0);	// one for each boundary including inflow
		tree_id_to_boundary_id.resize(surface_boundaries.size(),0);	// one for each boundary including inflow
	}

	unsigned int tree_start = 1;
	if(es->parameters.get<bool> ("use_centreline_data"))
		tree_start = 0;

	for(unsigned int m=tree_start; m<es->parameters.get<unsigned int> ("num_1d_trees")+1; m++)
	{

		std::stringstream node_file;
		std::stringstream edge_file;

		if(es->parameters.get<bool> ("match_1d_mesh_to_3d_mesh"))
		{

			node_file << es->parameters.get<std::string> ("input_1d_file") << "_" << m << ".node";//add number to the stream
			edge_file << es->parameters.get<std::string> ("input_1d_file") << "_" << m << ".edge";//add number to the stream
		}
		else
		{

			node_file << es->parameters.get<std::string> ("input_1d_file") << ".node";//add number to the stream
			edge_file << es->parameters.get<std::string> ("input_1d_file") << ".edge";//add number to the stream

		}

		std::cout << "node_file = " << node_file.str() << std::endl;
		std::cout << "edge_file = " << edge_file.str() << std::endl;

		std::ifstream infile_node(node_file.str().c_str());
		std::ifstream infile_edge(edge_file.str().c_str());

	
		std::cout << "1";
		// okay need to parse into a vector

		//node num is idx, 0-xcoord 1-ycoord 2-zcoord 3-radius
		//edge num is idx, 0-elem_num 1-node1 2-node2 3-generation 4-order

		// parse node fileclear
		{
			std::string line;
			//read first line with node num
			{
				std::getline(infile_node,line);
				std::stringstream line_stream(line);
				line_stream >> num_nodes_1d;
			}
	
			node_1d_data.resize(0);
			while(std::getline(infile_node,line))
			{
				std::vector<double> node_data;
				std::stringstream line_stream(line);
		
				unsigned int node_num; //unused
				double x_coord = 0.;
				double y_coord = 0.;
				double z_coord = 0.;
				double radius = 0.;
				double sub_tree_number = 0;
				line_stream >> node_num;
				line_stream >> x_coord;
				line_stream >> y_coord;
				line_stream >> z_coord;
				line_stream >> radius;
				line_stream >> sub_tree_number;

				//scale the data appropriately, to SI units from the mesh
				x_coord *= es->parameters.get<double> ("mesh_input_scaling_1d");
				y_coord *= es->parameters.get<double> ("mesh_input_scaling_1d");
				z_coord *= es->parameters.get<double> ("mesh_input_scaling_1d");
				radius *= es->parameters.get<double> ("mesh_input_scaling_1d");

				//okay, now that we have the mesh in SI units, we need to scale the mesh to the length scale we want to solve on
				//now in dimensionless units
				x_coord /= es->parameters.get<double> ("length_scale");
				y_coord /= es->parameters.get<double> ("length_scale");
				z_coord /= es->parameters.get<double> ("length_scale");
				radius /= es->parameters.get<double> ("length_scale");

				node_data.push_back(x_coord);
				node_data.push_back(y_coord);
				node_data.push_back(z_coord);
				node_data.push_back(radius);
				node_data.push_back(sub_tree_number);

				node_1d_data.push_back(node_data);

			}
		}

		std::cout << "2";	

		// parse edge file
		{
			std::string line;
			//read first line with node num
			std::getline(infile_edge,line);
			std::stringstream line_stream(line);
			line_stream >> num_edges_1d;
	
			edge_1d_data.resize(0);
			while(std::getline(infile_edge,line))
			{
				std::vector<double> edge_data;
				std::stringstream line_stream(line);
		
				unsigned int edge_num = 0;
				unsigned int node_1 = 0.;
				unsigned int node_2 = 0.;
				double radius = 0.;
				unsigned int generation = 0.;
				unsigned int order = 0.;

				line_stream >> edge_num;
				line_stream >> node_1;
				line_stream >> node_2;
				line_stream >> radius;
				line_stream >> generation;
				line_stream >> order;

				if(es->parameters.set<bool> ("radius_on_edge") && radius < 1e-10)
				{
					std::cout << "error, radius < 0 on elem" << edge_num << " EXItin" << std::endl; 
					std::cout << "node_1 = " << node_1 << std::endl;
					std::cout << "node_2 = " << node_2 << std::endl;
					std::cout << "radius = " << radius << std::endl;
					std::cout << "generation = " << generation << std::endl;
					std::cout << "order = " << order << std::endl;
					std::exit(0);
				}

				radius *= es->parameters.get<double> ("mesh_input_scaling_1d");
				radius /= es->parameters.get<double> ("length_scale");

				edge_data.push_back(edge_num);
				edge_data.push_back(node_1);
				edge_data.push_back(node_2);
				edge_data.push_back(radius);
				edge_data.push_back(generation);
				edge_data.push_back(order);

				edge_1d_data.push_back(edge_data);
			}
		}

		//edge 3 to 4 and 4 to 5

		std::cout << "3";	
		//do vertex and element stuff before reordering the vector
		//std::cout << "num_nodes_1d = " << num_nodes_1d << std::endl;
		std::vector<Point> vertices(num_nodes_1d);
		std::vector<std::vector<unsigned int> > cell_vertices(num_edges_1d);


		for(unsigned int i=0; i<num_nodes_1d; i++)
		{
			Point node(node_1d_data[i][0],node_1d_data[i][1],node_1d_data[i][2]);
			vertices[i] = node;
		}

		for(unsigned int i=0; i<num_edges_1d; i++)
		{
			std::vector<unsigned int> segment;
			segment.push_back((unsigned int)edge_1d_data[i][1]);
			segment.push_back((unsigned int)edge_1d_data[i][2]);

			cell_vertices[i] = segment;
		}

		std::cout << "4";	
		// use custom comparator to sort by generation
		// make a vector of the generations


	/*
		std::vector<int> generation_vector;
		std::vector<int> index;
		for(unsigned int i=0; i<edge_1d_data.size(); i++)
		{
			generation_vector.push_back(edge_1d_data[i][3]);
			index.push_back(i);
		}

		// NOT WORKING NOW! NEW COMPILER?
		std::sort(index.begin(),index.end(), 
							boost::bind([&generation](int a, int b) { 
									return generation[a] < generation[b]; 
			}));
	*/

		/*
		std::vector<int> generation_vector;
		std::vector<std::pair<int,int> > index;
		for(unsigned int i=0; i<edge_1d_data.size(); i++)
		{
			generation_vector.push_back(edge_1d_data[i][3]);
			index.push_back(std::pair<int,int>(edge_1d_data[i][3],i));
		}

		// NOT WORKING NOW! NEW COMPILER?
		std::sort(index.begin(),index.end(), boost::bind([](std::pair<int,int>& a, std::pair<int,int>& b)
				{ return (a.first < b.first); } ) );
	*/


		// sort generation array so that can figure out how many airways in each generation etc
		// while we are doing this we can also find the minimum generation number

		unsigned int min_generation = 1000;
		unsigned int max_generation = 0;
		int generation_array[edge_1d_data.size()];
		int index[edge_1d_data.size()];
		for(unsigned int i=0; i<edge_1d_data.size(); i++)
		{
			if((unsigned int)edge_1d_data[i][4] < min_generation)
				min_generation = (unsigned int)edge_1d_data[i][4];

			if((unsigned int)edge_1d_data[i][4] > max_generation)
				max_generation = (unsigned int)edge_1d_data[i][4];
			generation_array[i] = (unsigned int)edge_1d_data[i][4];
			index[i] = i;
		}

		std::cout << "max_generation = " << max_generation << std::endl;

	
	
		PetscSortIntWithArray(edge_1d_data.size(),generation_array,index);

		//exit(0);


		std::vector<std::vector<double> > temp_edge_1d_data = edge_1d_data;
	
		std::cout << "5";	
		for(unsigned int i=0; i<edge_1d_data.size(); i++)
		{
			edge_1d_data[i] = temp_edge_1d_data[index[i]];
		}
	

		//now should be sorted by generation, this is so that it is easier computationally to find parents etc
		//now we need vector with the starting positions of each generation
		std::vector<int> generation_start;
		int current_gen = min_generation - 1;
	
		for(unsigned int i=0; i<edge_1d_data.size(); i++)
		{
			if((int)edge_1d_data[i][4] > current_gen)
			{
				generation_start.push_back(i);
				current_gen = (int)edge_1d_data[i][4];
			
			}
		}
	
		std::cout << "6";	
		//set the number of generations, num generations is NOT max generation number.
		es->parameters.set<unsigned int> ("num_generations_1") = generation_start.size();

		// this is needed for the particle deposition
		if(es->parameters.get<unsigned int> ("alveolated_1d_tree"))
			num_generations.push_back(max_generation + 1 + es->parameters.get<unsigned int> ("num_alveolar_generations"));
		else
		{
			//num_generations.push_back(generation_start.size());
			if(m>0)
				num_generations.push_back(max_generation + 1);
		}
	

		//need a dummy generation denoting the end
		generation_start.push_back(edge_1d_data.size());

	

		//now we need to create the element_data vector
		// we now add an extra bit to the element data vector with the flow and current efficiency
		// we can also get the point at which 

		Point coupling_point;
		std::vector<Airway> airway_data_temp(edge_1d_data.size(),Airway());
		for(unsigned int i=0; i<edge_1d_data.size(); i++)
		{
			unsigned int elem_number = (unsigned int)edge_1d_data[i][0];
			unsigned int generation = (unsigned int)edge_1d_data[i][4];
			unsigned int order = (unsigned int)edge_1d_data[i][5];
			unsigned int node_1 = (unsigned int)edge_1d_data[i][1];
			unsigned int node_2 = (unsigned int)edge_1d_data[i][2];
			airway_data_temp[elem_number].set_local_elem_number(elem_number);
			airway_data_temp[elem_number].set_tree_number(m);

			//okay let us find the parent number
			if(generation > min_generation)	//has a parent
			{
				for(int j=generation_start[generation-min_generation-1]; j<generation_start[generation-min_generation]; j++)
				{
					//if first node same as second node then parent
					if(node_1 == (unsigned int)edge_1d_data[j][2])
					{
						airway_data_temp[elem_number].set_parent((unsigned int)edge_1d_data[j][0]);
					}
				}
			}
			else	//has no parent
			{

				// we should check whether it has a sibling, check in it's own generation, but not itself
				// we can also set whether it is daughter 1 in a first come first serve manner
				// okay so this now needs to be able to handle any number of sibling
				std::vector<unsigned int> siblings;
				unsigned int min_sibling_idx = 99999999;	// large number
				for(int j=generation_start[generation-min_generation]; j<generation_start[generation-min_generation+1]; j++)
				{
					if(node_1 == (unsigned int)edge_1d_data[j][1] && elem_number != (unsigned int)edge_1d_data[j][0])
					{
						// add this sibling to the list
						siblings.push_back(j);
						if((unsigned int)edge_1d_data[j][0] < min_sibling_idx)
						{
							min_sibling_idx = (unsigned int)edge_1d_data[j][0];
						}
					}
				}

				// check whether is daughter 1, has an elem number lower than the rest
				// also set the sibling data
				if(elem_number < min_sibling_idx)
				{
					airway_data_temp[elem_number].set_is_daughter_1(true);	//give it the number of siblings or just 1 if no siblings	

					// if we are daughter 1 then we need to set all the siblings at the end of element_data_temp
					for(unsigned int j=0; j<siblings.size(); j++)
					{
						// set the sibling
						airway_data_temp[elem_number].add_sibling((unsigned int)edge_1d_data[siblings[j]][0]);
					}

				}
				else
				{

					airway_data_temp[elem_number].set_is_daughter_1(false);

					// if we are daughter 2 then we just need to set the sibling to the one immediately before it
					// don't need the first one
					for(unsigned int j=siblings.size()-1; j>=0; j--)
					{
						airway_data_temp[elem_number].add_sibling((unsigned int)edge_1d_data[siblings[j]][0]);
						// set the sibling to the one just before we went past it...
						if((unsigned int)edge_1d_data[siblings[j]][0] < elem_number && !airway_data_temp[elem_number].has_primary_sibling())
						{
							airway_data_temp[elem_number].set_primary_sibling((unsigned int)edge_1d_data[siblings[j]][0]);
							break;
						}
					}

				}
			
				// the coupling point because this has no parent
				coupling_point = Point(node_1d_data[node_1][0],node_1d_data[node_1][1],node_1d_data[node_1][2]);
				std::cout << "coupling_point [" << m << "] = " << coupling_point << std::endl;
			
			}


			//okay let us find the daughter numbers
			//generation size is like one more than you'd expect so that 
			if(generation-min_generation < generation_start.size() - 2)
			{

				std::vector<unsigned int> daughters;	// a full list of the daughters of this airway
				for(int j=generation_start[generation-min_generation+1]; j<generation_start[generation-min_generation+2]; j++)
				{
					//if second node same as first node then daughter. daughter 1 is the first in the list
					if(node_2 == (unsigned int)edge_1d_data[j][1])
					{
						daughters.push_back((unsigned int)edge_1d_data[j][0]);
						airway_data_temp[elem_number].add_daughter(daughters.back());
						//if given daughter 1 status let us also tell that daughter it is daughter 1 and set daughter 1
						if(daughters.size() == 1)
						{
							airway_data_temp[daughters[0]].set_is_daughter_1(true);
						}
						else
						{
							airway_data_temp[daughters.back()].set_is_daughter_1(false);
						}
					}
				
				}

				
				// once the daughters have been set (if), tell them they are siblings
				if(daughters.size() > 0)
				{
					// give daughter 1 her siblings
					for(unsigned int j=1; j<daughters.size(); j++)
					{
						airway_data_temp[daughters[0]].add_sibling(daughters[j]);
					}

					// now let us tell the rest of the daughters who their primary siblings are, should be in order already
					for(unsigned int j=1; j<daughters.size(); j++)
					{
						airway_data_temp[daughters[j]].set_primary_sibling(daughters[j-1]);
						for(unsigned int k=0; k<daughters.size(); k++)
							if(j != k)
								airway_data_temp[daughters[j]].add_sibling(daughters[k]);
								
					}
				}

			}
			else
			{
				//has no daughters - do nothing
			}


			// set the length of the segment
			Point point_1(node_1d_data[node_1][0],node_1d_data[node_1][1],node_1d_data[node_1][2]);
			Point point_2(node_1d_data[node_2][0],node_1d_data[node_2][1],node_1d_data[node_2][2]);
			Point difference = point_1 - point_2;
			double length = difference.size();

			// set the radius to the average
			double radius = 0.;
			if(es->parameters.set<bool> ("radius_on_edge"))
			{
				if(edge_1d_data[i][3] < 1e-10)
				{
					std::cout << "ERROR: radius < 0 on elem " << i << ", radius = " << edge_1d_data[i][3] <<  " EXITING" << std::endl;
					std::exit(0);
				}
				radius = edge_1d_data[i][3];
			}
			else
			{
				
				if(node_1d_data[node_1][3] + node_1d_data[node_2][3] < 1e-10)
				{
					std::cout << "ERROR: radius < 0 on nodes " << node_1 << " and " << node_2 << " EXITING" << std::endl;
					std::exit(0);
				}
				radius = 0.5 * (node_1d_data[node_1][3] + node_1d_data[node_2][3]);

			}

			airway_data_temp[elem_number].set_node_1(point_1);
			airway_data_temp[elem_number].set_node_2(point_2);

			airway_data_temp[elem_number].set_length(length);
			airway_data_temp[elem_number].set_radius(radius);

			// set the generation and order of the airway
			airway_data_temp[elem_number].set_generation(generation);
			airway_data_temp[elem_number].set_order(order);




			//check if this is the terminating edge between 3d and 1d, imaging and generated data
			if(es->parameters.set<unsigned int> ("num_1d_trees") == 1 && es->parameters.get<bool> ("calculate_1d_info_at_coupling_nodes"))
			{
				unsigned int subtree_number_1 = (unsigned int)node_1d_data[node_1][4];
				unsigned int subtree_number_2 = (unsigned int)node_1d_data[node_2][4];
				// if we are at the end of 3d and the beginning of 1d
				if(subtree_number_1 == 0 && subtree_number_2 != 0)
				{
					if(subtree_number_2 >= subtree_starting_elements.size())
					{
						subtree_starting_elements.resize(subtree_number_2+1,-1);
					}
					subtree_starting_elements[subtree_number_2] = elem_number;
					flux_values_1d.push_back(0.);
					pressure_values_1d.push_back(0.);

					// add the coupling point to the array
					if(coupling_points.size() < subtree_number_2+1)
						coupling_points.resize(subtree_number_2+1);				// resize to subtree_number_2+1 because have point for inflow boundary

					coupling_points[subtree_number_2] = point_2;



					std::cout << "coupling_point[" << subtree_number_2 << "] = " << point_2 << std::endl;
					std::cout << "subtree_starting_elements[" << subtree_number_2 << "] = " << subtree_starting_elements[subtree_number_2] << std::endl;




				}
			}





			// if this is the 3d centreline data then we need to populate the 
			// centreline_terminal_id_to_boundary_id vector
			if(m == 0 && es->parameters.get<bool> ("use_centreline_data"))
			{

				unsigned int subtree_number_2 = (unsigned int)node_1d_data[node_2][4];

				centreline_terminal_id_to_tree_id.resize(airway_data_temp.size(),-1);	
				if(subtree_number_2 > 0)
				{
					centreline_terminal_id_to_tree_id[elem_number] = subtree_number_2;
				}

				centreline_points.resize(airway_data_temp.size());
				std::vector<Point> points;
				points.push_back(point_1);
				points.push_back(point_2);
				centreline_points[elem_number] = points;
			}
			


		}

		// okay if we are doing a 3d1d then the coupling points
		if(es->parameters.get<bool> ("match_1d_mesh_to_3d_mesh"))
		{
			// add the coupling point to the array
			if(coupling_points.size() < m+1)
				coupling_points.resize(m+1);				// resize to subtree_number_2+1 because have point for inflow boundary

			coupling_points[m] = coupling_point;
		}

		//std::cout << "humm" << std::endl;



		
		// if we are reading in the data from the 3d part then we just need to put it in centreline_element_data
		if(m == 0 && es->parameters.get<bool> ("use_centreline_data"))
		{
			centreline_airway_data = airway_data_temp;

		}
		else
		{

			calculate_num_alveloar_generations_for_tree();


			std::cout << "7" << std::endl;	

			/*
			for(unsigned int i=0; i<10;i++)
			{
				for(unsigned int j=0; j<5; j++)
					std::cout << "element_data_temp[" << i << "][" << j << "] = " << element_data_temp[i][j] << std::endl;
			}
			*/


			//airway_data_temp.

			// move all the indices of element data up by num_1d_elements
			for(unsigned int i=0; i<airway_data_temp.size();i++)
			{
				airway_data_temp[i].move_element_numbers(num_1d_elements);
			}

			std::cout << "num_1d_elements = " << num_1d_elements << std::endl;


			airway_data.insert( airway_data.end(), airway_data_temp.begin(), airway_data_temp.end() );

			num_1d_elements += airway_data_temp.size();

	
			unsigned int subdomain_id = m;

			// if there is a 3D subdomain, give the 1D tree a subdomain starting after the last 3D
			if(subdomains_3d.size() > 0)
				subdomain_id += subdomains_3d.back();

			// make boundary id 0 by default, i.e. only 1d mesh
			unsigned int boundary_id = 0;
			// here we need to figure out what boundary to give it by comparing to centroid of something.
			if(es->parameters.get<bool> ("match_1d_mesh_to_3d_mesh"))
			{
				std::cout << "coupling_point is " << coupling_point << std::endl;
				Point centroid;
				double max_radius;
				// don't bother checking inflow boundary

				for(unsigned int i=1; i<surface_boundaries.size(); i++)
				{
					centroid = surface_boundaries[i]->get_centroid();
					max_radius = surface_boundaries[i]->get_max_radius();
					double distance = (centroid - coupling_point).size();

					//std::cout << "distance = " << distance << std::endl;
					//std::cout << "max_radius = " << max_radius << std::endl;
					if(distance < max_radius)
					{
						std::cout << "coupling surface found = " << i << std::endl;
						boundary_id_to_tree_id[i] = m;
						boundary_id = i;
						break;
					}

					if(i == surface_boundaries.size()-1)
					{
						std::cout << "error coupling surface not found" << std::endl;
						std::exit(0);
					}
				}
			}


			// give it boundary id 1
			add_1d_tree_to_mesh(vertices,cell_vertices, subdomain_id,boundary_id);
			mesh.prepare_for_use (/*skip_renumber =*/ false);

			subdomains_1d.push_back(subdomain_id);


			// well although the 1d stuff doesn't actually have a boundary at 0
			// we still consider it
			if(!es->parameters.get<bool> ("calculate_1d_info_at_coupling_nodes"))
			{
				pressure_values_1d.push_back(0.0);	//inflow
				flux_values_1d.push_back(0.0);	//inflow
			}
		}


	}
	//pressure_values_1d.push_back(0.0);	//outflow
	//flux_values_1d.push_back(0.0);	//outflow

	for(unsigned int i=0; i<boundary_id_to_tree_id.size(); i++)
	{
		tree_id_to_boundary_id[boundary_id_to_tree_id[i]] = i;
		std::cout << "boundary id " << i << " assoc with tree " << boundary_id_to_tree_id[i] << std::endl;
	}

	for(unsigned int i=0; i<subtree_starting_elements.size(); i++)
	{
		std::cout << "subtree_starting_elements[" << i << "] = " << subtree_starting_elements[i] << std::endl;
	}

	/*
	std::cout << "CENTRELINE DATA" << std::endl;
	for(unsigned int i=0; i<centreline_airway_data.size(); i++)
	{
		std::cout << "i = " << i << std::endl;
		for(unsigned int j=0; j<centreline_airway_data[i].size(); j++)
		{
			std::cout << " " << centreline_airway_data[i][j] << std::endl;
		}
	}
	*/

	/*
	for(unsigned int i=0; i<airway_data.size(); i++)
	{
		std::cout << "i = " << i << std::endl;
		airway_data[i].print_concise();
	}
	*/	

	if(es->parameters.get<bool> ("match_1d_mesh_to_3d_mesh") || es->parameters.get<bool> ("calculate_1d_info_at_coupling_nodes"))
	{
		std::cout << "about to output coupling points" << std::endl;
		output_coupling_points();
		std::cout << "outputted coupling points" << std::endl;
	}

	std::cout << " finished setting up 1d mesh from file" << std::endl;
}



void NavierStokesCoupled::calculate_num_alveloar_generations_for_tree()
{
	for(unsigned int i=0; i<airway_data.size(); i++)
		if(!airway_data[i].has_daughter_1())
			airway_data[i].set_num_alveolar_generations(es->parameters.get<unsigned int> ("num_alveolar_generations"));
}





void NavierStokesCoupled::generate_1d_mesh ()
{

	std::cout << "in generate_1d_mesh" << std::endl;
	std::cout << "num 1d trees = " << es->parameters.get<unsigned int> ("num_1d_trees") << std::endl;

	//so to create a 1d_tree we first create the first 1d segment and then grow the tree from that
	//these are simply parameters of the tree making thing so can be gotten from the input file
	// initial segment length, chosen so that wall thickness is always positive
	double length_scale = es->parameters.get<double> ("length_scale");
	double initial_segment_length = es->parameters.get<double> ("initial_segment_length");
	initial_segment_length /= es->parameters.get<double> ("length_scale");

	// some vectors that will temporarily be used when creating the tree, element_data
	// is owned globally.
	std::vector<Point> vertices;
	std::vector<std::vector<unsigned int> > cell_vertices;
	std::vector<unsigned int> segment;
	Point p0,p1;

	unsigned int num_generations_1 = es->parameters.get<unsigned int> ("num_generations_1");
	unsigned int num_generations_2 = es->parameters.get<unsigned int> ("num_generations_2");
	unsigned int num_generations_3 = es->parameters.get<unsigned int> ("num_generations_3");

	double length_diam_ratio = es->parameters.get<double> ("length_diam_ratio");

	//we need to make this tree in dimensionless units, assume parameters are in SI units
	
	// if we are using an arbitrary tree mesh and 
	if((es->parameters.get<unsigned int> ("geometry_type") == 4 || es->parameters.get<unsigned int> ("geometry_type") == 2) && sim_3d)
	{

		std::cout << "in generate_1d_mesh generating " << es->parameters.get<unsigned int> ("num_1d_trees") << " trees" << std::endl;

		//subdomain counting
		unsigned int subdomain_id = 0;

		// if there is a 3D subdomain, give the 1D tree a subdomain starting after the last 3D
		if(subdomains_3d.size() > 0)
			subdomain_id = subdomains_3d.back() + 1;

		// loop over all the boundaries except the inflow boundary
		// now the subdomain ids won't match up with the boundary ids
		for(unsigned int i=1; i < boundary_ids.size(); i++)
		{

			vertices.clear();
			cell_vertices.clear();
			segment.clear();

			Point normal = surface_boundaries[i]->get_normal();
			Point centroid = surface_boundaries[i]->get_centroid();
			double area = surface_boundaries[i]->get_area();


			// ********************************************************************* //
			// ********************** TREE NUMBER i ******************************** //
			// ********************************************************************* //

			// start at centroid
			p0 = centroid;
			std::cout << "centroid = " << centroid << std::endl;
			std::cout << "area = " << area << std::endl;
			//next point should be in direction of normal, with length based on diameter
			double approx_diam = 0.;
			if(es->parameters.get<bool> ("threed"))
				approx_diam = 2* sqrt(area/M_PI); // diameter approximated by assuming circular outflow
			else
				approx_diam = area; // diameter approximated by assuming circular outflow
				
			if(es->parameters.get<bool> ("initial_segment_length_from_mesh"))
			{
				initial_segment_length = length_diam_ratio * approx_diam;
				if(es->parameters.get<bool> ("half_initial_length"))
					initial_segment_length *= 0.5;

				std::cout << "initial_sgment_length = " << initial_segment_length << std::endl;
				std::cout << "approx_diam = " << approx_diam << std::endl;
			}
			p1 = centroid + normal * initial_segment_length;

			vertices.push_back(p0);
			vertices.push_back(p1);
			segment.push_back(0);
			segment.push_back(1);
			cell_vertices.push_back(segment);

			unsigned int num_generations_local = num_generations_1;
			if(i==1)
				num_generations_local = num_generations_1;
			else if(i==2)
				num_generations_local = num_generations_2;
			else if(i==3)
				num_generations_local = num_generations_3;

			// for each element/segment, assuming same numbering convention
			// 0 - parent elem no, 1 - daughter elem no 1, 2 - daughter elem no 2
			// 3 - sibling elem no, 4 - is daughter_1 bool 5 - length, 6 - radius
			element_data.push_back(std::vector<double>(12,-1));
			element_data[element_data.size() - 1][0] = -1;
			element_data[element_data.size() - 1][5] = initial_segment_length;

			double radius = element_data[element_data.size() - 1][5]/length_diam_ratio/2.0;
			if(es->parameters.get<bool> ("half_initial_length"))
			{
				radius *= 2.0;
			}

			// if we are doing a twod approx then we need to change the radius
			if(es->parameters.get<bool> ("twod_oned_tree"))
				radius = pow(16/3/M_PI*pow(radius,3.0),1./4.0);

			element_data[element_data.size() - 1][6] = radius;	//this is the radius no diameter
			element_data[element_data.size() - 1][7] = 0;	//generation
			element_data[element_data.size() - 1][8] = num_generations_local;	//generation
			//particle deposition stuff
			element_data[element_data.size() - 1][9] = 0;	//the flow rate in this tube for the current time step
			element_data[element_data.size() - 1][10] = 0;	//the efficiency in this tube for the current time step
			element_data[element_data.size() - 1][11] = 0;	//num alveolar generations


			create_1d_tree(vertices,cell_vertices,num_generations_local);
			num_generations.push_back(num_generations_local);

			// this is needed for the particle deposition in alveoli
			if(es->parameters.get<unsigned int> ("alveolated_1d_tree"))
				num_generations.push_back(num_generations_local + es->parameters.get<unsigned int> ("num_alveolar_generations"));
			else
				num_generations.push_back(num_generations_local);

	
			add_1d_tree_to_mesh(vertices,cell_vertices, subdomain_id, i);

			subdomains_1d.push_back(subdomain_id);
			subdomain_id++;
			
		}

		// only need to do this once
		if(es->parameters.get<unsigned int> ("alveolated_1d_tree"))
			calculate_num_alveloar_generations_for_tree();
		
	}
	else if(es->parameters.get<unsigned int> ("num_1d_trees") == 2)
	{

		// ********************************************************************* //
		// ********************** TREE NUMBER 1 ******************************** //
		// ********************************************************************* //

		// this will of course be read into the file somehow.. lol, calculate the end from the mesh file perhaps?
		p0 = Point(2.5e-2,0.0,-3.3e-2);		//let us say in centemetres
		p0 /= es->parameters.get<double> ("length_scale");
		//this point was before mesh was scaled so also scale it
		//p0 = 5.6*p0;

		//initial segment angle
		p1 = Point (p0(0) + initial_segment_length*sin(3*M_PI/8.0),0.0,
								p0(2) - initial_segment_length*cos(3*M_PI/8.0));	// length chosen so that wall thickness is always positive

		vertices.push_back(p0);
		vertices.push_back(p1);
		segment.push_back(0);
		segment.push_back(1);
		cell_vertices.push_back(segment);

		// for each element/segment, assuming same numbering convention
		// 0 - parent elem no, 1 - daughter elem no 1, 2 - daughter elem no 2
		// 3 - sibling elem no, 4 - is daughter_1 bool 5 - length, 6 - radius
		element_data.push_back(std::vector<double>(12,-1));
		element_data[0][0] = -1;
		element_data[0][5] = initial_segment_length;
		element_data[0][6] = element_data[0][5]/length_diam_ratio/2.0;	//this is the radius no diameter
		element_data[0][7] = 0;	//generation
		element_data[0][8] = num_generations_1;	//generation
		//particle deposition stuff
		element_data[0][9] = 0;	//the flow rate in this tube for the current time step
		element_data[0][10] = 0;	//the efficiency in this tube for the current time step
		element_data[0][11] = 0;	//num alveolar generations

	
	
		
		create_1d_tree(vertices,cell_vertices,num_generations_1);

		// this is needed for alveolar deposition
		if(es->parameters.get<unsigned int> ("alveolated_1d_tree"))
			num_generations.push_back(num_generations_1 + es->parameters.get<unsigned int> ("num_alveolar_generations"));
		else
			num_generations.push_back(num_generations_1);
	
		unsigned int subdomain_id = 0;

		// if there is a 3D subdomain, give the 1D tree a subdomain starting after the last 3D
		if(subdomains_3d.size() > 0)
			subdomain_id = subdomains_3d.back() + 1;

		// give it a boundary id of 1
		add_1d_tree_to_mesh(vertices,cell_vertices, subdomain_id,1);
		subdomains_1d.push_back(subdomain_id);
		subdomain_id++;

		// ********************************************************************* //
		// ******************* TREE NUMBER 2 ************************************** //
		// ********************************************************************* //

		// the next tree, will of course be done in a way such that can be read in and given the correct ids etc
		// clear some vectors
		vertices.clear();
		cell_vertices.clear();
		segment.clear();
	
		// this will of course be read into the file somehow.. lol, calculate the end from the mesh file perhaps?
		p0 = Point (-1.3e-2,0.0,-4.3e-2);		//let us say in centimetres
		p0 /= es->parameters.get<double> ("length_scale");
		//this point was before mesh was scaled so also scale it
		//p0 = 5.6*p0;
	
		//initial segment angle
		p1 = Point (p0(0) - initial_segment_length*sin(M_PI/8.0),0.0,
								p0(2) - initial_segment_length*cos(M_PI/8.0));	// length chosen so that wall thickness is always positive
	
		vertices.push_back(p0);
		vertices.push_back(p1);
		segment.push_back(0);
		segment.push_back(1);
		cell_vertices.push_back(segment);

		// for each element/segment, assuming same numbering convention
		// 0 - parent elem no, 1 - daughter elem no 1, 2 - daughter elem no 2
		// 3 - sibling elem no, 4 - is daughter_1 bool 5 - length, 6 - radius
		element_data.push_back(std::vector<double>(12,-1));
		element_data[element_data.size() - 1][0] = -1;
		element_data[element_data.size() - 1][5] = initial_segment_length;
		element_data[element_data.size() - 1][6] = element_data[element_data.size() - 1][5]/length_diam_ratio/2.0;	// radius not diam
		element_data[element_data.size() - 1][7] = 0;	//generation
		element_data[element_data.size() - 1][8] = num_generations_2;	//generation
		//particle deposition stuff
		element_data[element_data.size() - 1][9] = 0;	//the flow rate in this tube for the current time step
		element_data[element_data.size() - 1][10] = 0;	//the efficiency in this tube for the current time step
		element_data[element_data.size() - 1][11] = 0;	//num_alveolar_generations

		create_1d_tree(vertices,cell_vertices,num_generations_2);

		// this is needed for alveolar deposition
		if(es->parameters.get<unsigned int> ("alveolated_1d_tree"))
			num_generations.push_back(num_generations_2 + es->parameters.get<unsigned int> ("num_alveolar_generations"));
		else
			num_generations.push_back(num_generations_2);
	

		// only need to do this once
		if(es->parameters.get<unsigned int> ("alveolated_1d_tree"))
			calculate_num_alveloar_generations_for_tree();

		// subdomain has iterated already
		// give it a boundary id of 2
		add_1d_tree_to_mesh(vertices,cell_vertices,subdomain_id,2);

		subdomains_1d.push_back(subdomain_id);

	}
	else if(es->parameters.get<unsigned int> ("num_1d_trees") == 1)
	{

		std::cout << "in generate_1d_mesh generating 1 tree" << std::endl;
		// ********************************************************************* //
		// ********************** TREE NUMBER 1 ******************************** //
		// ********************************************************************* //

		// this will of course be read into the file somehow.. lol, calculate the end from the mesh file perhaps?
		p0 = Point(0.0,0.0,0.0);
		p0 /= es->parameters.get<double> ("length_scale");
		//this point was before mesh was scaled so also scale it
		//p0 = 5.6*p0;

		//initial segment angle
		p1 = Point (0.,0.0,initial_segment_length);	// length chosen so that wall thickness is always positive

		vertices.push_back(p0);
		vertices.push_back(p1);
		segment.push_back(0);
		segment.push_back(1);
		cell_vertices.push_back(segment);

		// for each element/segment, assuming same numbering convention
		// 0 - parent elem no, 1 - daughter elem no 1, 2 - daughter elem no 2
		// 3 - sibling elem no, 4 - is daughter_1 bool 5 - length, 6 - radius
		element_data.push_back(std::vector<double>(12,-1));
		element_data[0][0] = -1;
		element_data[0][5] = initial_segment_length;
		element_data[0][6] = element_data[0][5]/length_diam_ratio/2.0;	//this is the radius no diameter
		element_data[0][7] = 0;	//generation
		element_data[0][8] = num_generations_1;	//generation
		//particle deposition stuff
		element_data[0][9] = 0;	//the flow rate in this tube for the current time step
		element_data[0][10] = 0;	//the efficiency in this tube for the current time step
		element_data[0][11] = 0;	//num_alveolar_generations

		create_1d_tree(vertices,cell_vertices,num_generations_1);

		// this is needed for alveolar deposition
		if(es->parameters.get<unsigned int> ("alveolated_1d_tree"))
			num_generations.push_back(num_generations_1 + es->parameters.get<unsigned int> ("num_alveolar_generations"));
		else
			num_generations.push_back(num_generations_1);

		// only need to do this once
		if(es->parameters.get<unsigned int> ("alveolated_1d_tree"))
			calculate_num_alveloar_generations_for_tree();


		unsigned int subdomain_id = 0;

		// if there is a 3D subdomain, give the 1D tree a subdomain starting after the last 3D
		if(subdomains_3d.size() > 0)
			subdomain_id = subdomains_3d.back() + 1;

		// give it a boundary id of 1
		add_1d_tree_to_mesh(vertices,cell_vertices, subdomain_id,1);

		subdomains_1d.push_back(subdomain_id);
	}

  // Done building the mesh.  Now prepare it for use.
  mesh.prepare_for_use (/*skip_renumber =*/ false);



	// well although the 1d stuff doesn't actually have a boundary at 0
	// we still consider it
	pressure_values_1d.push_back(-1.0);	//inflow
	for(unsigned int i=1; i<es->parameters.get<unsigned int>("num_1d_trees")+1; i++)
		pressure_values_1d.push_back(0.0);	//inflow

	flux_values_1d.push_back(0.0);	//inflow
	for(unsigned int i=1; i<es->parameters.get<unsigned int>("num_1d_trees")+1; i++)
		flux_values_1d.push_back(0.0);	//inflow
	

}


// creates a 1d tree and puts the data into vertices and cell_vertices
// element_data is owened globally. add_1d_tree_to_mesh should be called immediatelly afterwards.
void NavierStokesCoupled::create_1d_tree(std::vector<Point>& vertices, 
		std::vector<std::vector<unsigned int> >& cell_vertices,unsigned int num_generations)
{
	//tree params
	double bifurcation_angle = M_PI/4.0;
	double length_ratio = 1.0/1.25;	//halves the length of segments in each generation
	double left_length_ratio = 0.876;	//halves the length of segments in each generation
	double right_length_ratio = 0.686;	//halves the length of segments in each generation
	//const double length_diam_ratio = 3.0;
	const double length_diam_ratio = es->parameters.get<double> ("length_diam_ratio");
	left_length_ratio = length_ratio;
	right_length_ratio = length_ratio;

	//helpful variables
	std::vector<unsigned int> segment(2);
	unsigned int parent_start_node = 1;
	unsigned int parent_segment = 0;
	unsigned int element_offset = element_data.size() - 1;
	Point p0,p1;

	std::cout << "hey in create_1d_tree" << std::endl;

	//unsigned int element_offset = 0;
  for(unsigned int i=1; i<num_generations; i++)
	{

		std::cout << "hmmm" << std::endl;
		//for each parent we generate two new points and two new segments
		for(unsigned int j=0;j<pow(2,i-1); j++)
		{
			//we need to find the direction of this segment
			Point direction = vertices[cell_vertices[parent_segment + j][1]] 
												- vertices[cell_vertices[parent_segment + j][0]];
			Point unit_direction = direction.unit();

			Point p0_direction;
			if(es->parameters.get<bool> ("threed"))
			{
				p0_direction(0) = cos(bifurcation_angle)*unit_direction(0) + sin(bifurcation_angle)* unit_direction(2);
				p0_direction(2) = -sin(bifurcation_angle)*unit_direction(0) + cos(bifurcation_angle)* unit_direction(2);
			}
			else
			{
				p0_direction(0) = cos(bifurcation_angle)*unit_direction(0) + sin(bifurcation_angle)* unit_direction(1);
				p0_direction(1) = -sin(bifurcation_angle)*unit_direction(0) + cos(bifurcation_angle)* unit_direction(1);
			}

			Point p1_direction;
			if(es->parameters.get<bool> ("threed"))
			{
				p1_direction(0) = cos(-bifurcation_angle)*unit_direction(0) + sin(-bifurcation_angle)* unit_direction(2);
				p1_direction(2) = -sin(-bifurcation_angle)*unit_direction(0) + cos(-bifurcation_angle)* unit_direction(2);
			}
			else
			{
				p1_direction(0) = cos(-bifurcation_angle)*unit_direction(0) + sin(-bifurcation_angle)* unit_direction(1);
				p1_direction(1) = -sin(-bifurcation_angle)*unit_direction(0) + cos(-bifurcation_angle)* unit_direction(1);
			}

			double length_1 = element_data[parent_segment + j + element_offset][5]*left_length_ratio;
			double length_2 = element_data[parent_segment + j + element_offset][5]*right_length_ratio;

			// if the first segment is a half segment then the next should double
			if(i == 1 && es->parameters.get<bool> ("half_initial_length"))
			{
				length_1 *= 2.;
				length_2 *= 2.;
			}

			//p0 = vertices[parent_start_node + j] + pow(length_ratio,i)*p0_direction;
			//p1 = vertices[parent_start_node + j] + pow(length_ratio,i)*p1_direction;
			p0 = vertices[parent_start_node + j] + length_1*p0_direction;
			p1 = vertices[parent_start_node + j] + length_2*p1_direction;

			vertices.push_back(p0);
			//now let us make these two segments
			segment[0] = parent_start_node + j;
			segment[1] = vertices.size() - 1;	//should be the last vertex added
			cell_vertices.push_back(segment);

			//update the data
			element_data[parent_segment + j + element_offset][1] = cell_vertices.size() - 1 + element_offset;
			element_data.push_back(std::vector<double>(12,-1));
			element_data[cell_vertices.size() - 1 + element_offset][0] = parent_segment + j + element_offset;
			element_data[cell_vertices.size() - 1 + element_offset][3] = cell_vertices.size() + element_offset;
			element_data[cell_vertices.size() - 1 + element_offset][4] = 1;
			element_data[cell_vertices.size() - 1 + element_offset][5] = length_1;

			double radius_1 = length_1/length_diam_ratio/2.0;
			// if we are doing a twod approx then we need to change the radius
			if(es->parameters.get<bool> ("twod_oned_tree"))
				radius_1 = pow(16/3/M_PI*pow(radius_1,3.0),1./4.0);

			element_data[cell_vertices.size() - 1 + element_offset][6] = radius_1; //radius not diam
			element_data[cell_vertices.size() - 1 + element_offset][7] = i; //generation
			element_data[cell_vertices.size() - 1 + element_offset][8] = num_generations - i; //order
			//particle deposition stuff
			element_data[cell_vertices.size() - 1 + element_offset][9] = 0;	//the flow rate in this tube for the current time step
			element_data[cell_vertices.size() - 1 + element_offset][10] = 0;	//the efficiency in this tube for the current time step
			element_data[cell_vertices.size() - 1 + element_offset][11] = 0;	//num_alveolar_generations

			vertices.push_back(p1);
			segment[0] = parent_start_node + j;
			segment[1] = vertices.size() - 1;
			cell_vertices.push_back(segment);

			//update the data
			element_data[parent_segment + j + element_offset][2] = cell_vertices.size() - 1 + element_offset;
			element_data.push_back(std::vector<double>(12,-1));
			element_data[cell_vertices.size() - 1 + element_offset][0] = parent_segment + j + element_offset;
			element_data[cell_vertices.size() - 1 + element_offset][3] = cell_vertices.size() - 2 + element_offset;
			element_data[cell_vertices.size() - 1 + element_offset][4] = 0;
			element_data[cell_vertices.size() - 1 + element_offset][5] = length_2;

			double radius_2 = length_2/length_diam_ratio/2.0;
			// if we are doing a twod approx then we need to change the radius
			if(es->parameters.get<bool> ("twod_oned_tree"))
				radius_2 = pow(16/3/M_PI*pow(radius_2,3.0),1./4.0);

			element_data[cell_vertices.size() - 1 + element_offset][6] = radius_2; //radius not diam
			element_data[cell_vertices.size() - 1 + element_offset][7] = i; //generation
			element_data[cell_vertices.size() - 1 + element_offset][8] = num_generations - i; //order
			//particle deposition stuff
			element_data[cell_vertices.size() - 1 + element_offset][9] = 0;	//the flow rate in this tube for the current time step
			element_data[cell_vertices.size() - 1 + element_offset][10] = 0;	//the efficiency in this tube for the current time step
			element_data[cell_vertices.size() - 1 + element_offset][11] = 0;	//num_alveolar_generations

		}
	
		parent_start_node += pow(2,i-1);
		parent_segment 	  += pow(2,i-1);

	}

}




void NavierStokesCoupled::add_1d_tree_to_mesh(std::vector<Point>& vertices, 
		std::vector<std::vector<unsigned int> >& cell_vertices, unsigned int subdomain_id, unsigned int boundary_id)
{



	// we made the inflow bdys the same as the subdomain_id
	// we made the outflow bdys 1000
	// okay need to figure out how we are going to distribute the mesh_1d
	// need to take into account the connectivity, but essentially we could put separate trees onto separate cores

	//loop over segments adding nodes and stuff

	unsigned int
	n_vertices = vertices.size();
	unsigned int n_segments = cell_vertices.size();

	//hmmm we need to build a map from segment vertex to node number
	std::vector<unsigned int> segment_vertex_to_node(n_vertices,-1);
	std::vector<unsigned int> vertex_count(n_vertices,0);

	//need to construct a thing that counts up how many times a vertex was used
	// if we start on a segment then vertex_count will be 1,
	// if we start on a bifurcation then vertex count will be 2,
	// if we are in the middle of a mesh the vertex count will be 3,
	// if we are at the end of a tree then the vertex count will be 1
	for(unsigned int i=0;i<n_segments;i++)
	{
		vertex_count[cell_vertices[i][0]] += 1;
		vertex_count[cell_vertices[i][1]] += 1;
	}

	//std::cout << "max vertex count = " << *std::max_element(vertex_count.begin(),vertex_count.end()) << std::endl;
	//std::cout << "min vertex count = " << *std::min_element(vertex_count.begin(),vertex_count.end()) << std::endl;

	/*
	for(unsigned int i=0; i<vertex_count.size();i++)
		if(vertex_count[i] == 2)
			std::cout << "vertex_count = 2" << std::endl;
	*/

	ElemType type;
	unsigned int nx;
	if(es->parameters.set<bool>("0D"))
	{
		type = EDGE2;		//element type
		nx = 1;//2;	//no elements
	}
	else
	{
	
		type = EDGE2;		//element type
		nx = 2;//2;	//no elements
	}

	//keep track of the number of nodes from when the previous element added
	unsigned int node_start = mesh.n_nodes(); //0;
	for(unsigned int i=0;i<n_segments;i++)
	{

		//starting point
		Point start_point = vertices[cell_vertices[i][0]];
		Point end_point = vertices[cell_vertices[i][1]];
		Point direction = end_point - start_point;

		libmesh_assert_not_equal_to (nx, 0);

		// Reserve elements
		switch (type)
		{
		//this will stay the same because for each segment we want the number of elements simply adds
		case INVALID_ELEM:
		case EDGE2:
		case EDGE3:
		  {
		    mesh.reserve_elem (nx);
		    break;
		  }

		default:
		  {
		    libMesh::err << "ERROR: Unrecognized 1D element type." << std::endl;
		    libmesh_error();
		  }
		}

		// Build the nodes, depends on whether we're using linears,
		// quadratics or cubics and whether using uniform grid or Gauss-Lobatto
		// we will not consider gauss lobatto 
		switch(type)
		{
		  case INVALID_ELEM:
		  case EDGE2:
		    {
					//only for the first segment do we add the first node, otherwise from you know
					if(i==0)
					{
						mesh.add_point (start_point);

					}

					//ensures that right number of elements are added 
		      for (unsigned int j=1; j<=nx; j++)
		      {
		        mesh.add_point (start_point + static_cast<Real>(j)/static_cast<Real>(nx) * direction);

		      }
		      break;
		    }

		  case EDGE3:
		    {
					if(i==0)
					{
						mesh.add_point (start_point);

					}
		      for (unsigned int j=1; j<=2*nx; j++)
		      {
		      
		        mesh.add_point (start_point + static_cast<Real>(j)/static_cast<Real>(2*nx) * direction);
		      }
		      break;
		    }


		  default:
		    {
		      libMesh::err << "ERROR: Unrecognized 1D element type." << std::endl;
		      libmesh_error();
		    }

		}

		// Build the elements of the mesh
		switch(type)
		{
		  case INVALID_ELEM:
		  case EDGE2:
		    {
					//for the first segment we add all the elements normally
					//however for the others the first element needs to be connected to the 
					
		      for (unsigned int j=0; j<nx; j++)
					{
		        Elem* elem = mesh.add_elem (new Edge2);
        	  elem->subdomain_id() = subdomain_id;
						//if we are starting we need to create a new node otherwise not
						//and also give it the correct boundary id
						if(i==0)
						{
		        	elem->set_node(0) = mesh.node_ptr(node_start + j);
						}
						else
		        	elem->set_node(0) = mesh.node_ptr(segment_vertex_to_node[cell_vertices[i][0]]);
		        elem->set_node(1) = mesh.node_ptr(node_start + j +1);

						//inflow boundary
						if(j==0 && vertex_count[cell_vertices[i][0]] < 3)//vertex_count[cell_vertices[i][0]] == 1)
						{
							//std::cout << "setting boundary id, vertex_count = " << vertex_count[cell_vertices[i][0]] << std::endl;
							//std::cout << "elem = " << i << std::endl;
		        	mesh.boundary_info->add_side(elem, 0, boundary_id);
						}
						//sorta obvious??
						segment_vertex_to_node[0] = node_start;
		        if (j == (nx-1))
						{
							//hmmm need to now know a priori which elements are terminal somehow
							//if the vertex count of the current segments end is equal to one then
							// terminal and add it to the boundary information
							if(vertex_count[cell_vertices[i][1]] == 1)
		          	mesh.boundary_info->add_side(elem, 1, 1000);
							//also save the node of the end of the segment
							segment_vertex_to_node[cell_vertices[i][1]] = node_start + nx;
						}
		      }
		    break;
		    }

		  case EDGE3:
		    {

					//for the first segment we add all the elements normally
					//however for the others the first element needs to be connected to the 
					
		      for (unsigned int j=0; j<nx; j++)
		      {
						
		        Elem* elem = mesh.add_elem (new Edge3);
        	  elem->subdomain_id() = subdomain_id;

						if(i==0)
			        elem->set_node(0) = mesh.node_ptr(node_start + 2*j);
            else
              elem->set_node(0) = mesh.node_ptr(segment_vertex_to_node[cell_vertices[i][0]]);
		        elem->set_node(2) = mesh.node_ptr(node_start + 2*j + 1);
		        elem->set_node(1) = mesh.node_ptr(node_start + 2*j + 2);

						//inflow boundary
						if(j==0 && vertex_count[cell_vertices[i][0]] == 1)
		        	mesh.boundary_info->add_side(elem, 0, boundary_id);


		        if (j == (nx-1))
						{
							//hmmm need to now know a priori which elements are terminal somehow
							//if the vertex count of the current segments end is equal to one then
							// terminal and add it to the boundary information
							if(vertex_count[cell_vertices[i][1]] == 1)
		          	mesh.boundary_info->add_side(elem, 1, 1000);

							//also save the node of the end of the segment
							segment_vertex_to_node[cell_vertices[i][1]] = node_start + 2*nx;
						}

						//sorta obvious?? - the first one
						segment_vertex_to_node[0] = node_start;

					}

		    break;
		    }

		  default:
		    {
		      libMesh::err << "ERROR: Unrecognized 1D element type." << std::endl;
		      libmesh_error();
		    }

		}
		

		// Scale the nodal positions
		// no need for scaling in new formulation
		//for (unsigned int p=0; p<mesh.n_nodes(); p++)
		//  mesh.node(p)(1) = (mesh.node(p)(1))*(xmax-xmin) + xmin;


		//update where to start adding nodes from
		node_start = mesh.n_nodes() - 1;
	}


}






void NavierStokesCoupled::setup_1d_system(TransientLinearImplicitSystem * system)
{

	std::set<subdomain_id_type> active_subdomains;
	active_subdomains.clear(); 
	

	for(unsigned int i=0; i<subdomains_1d.size(); i++)
		active_subdomains.insert(subdomains_1d[i]);

	for(unsigned int i=0; i<subdomains_1d.size(); i++)
		std::cout << "subdomains_1d = " << subdomains_1d[i] << std::endl;

	int P_var = 0;
	int Q_var = 0;
	//int radius_var = 0;

	std::cout << "hi 0D prarmeter = " << es->parameters.set<bool>("0D") << std::endl;

	if(es->parameters.set<bool>("0D"))
	{
		P_var = system->add_variable ("P", FIRST,MONOMIAL,&active_subdomains);
		Q_var = system->add_variable ("Q", FIRST,MONOMIAL,&active_subdomains);
	}
	else
	{
		system->add_variable ("P", FIRST,LAGRANGE,&active_subdomains);
	}

	int radius_var = extra_1d_data_system->add_variable("radius", CONSTANT, MONOMIAL,&active_subdomains);
	int poiseuille_var = extra_1d_data_system->add_variable("poiseuille", CONSTANT, MONOMIAL,&active_subdomains);

	int efficiency_var = 0;
	// variables for particle deposition system
	if(particle_deposition == 2)
	{
		efficiency_var = particle_deposition_system_1d->add_variable ("efficiency", CONSTANT, MONOMIAL,&active_subdomains);
	}
	
	std::cout << "P_var = " << P_var << std::endl;
	std::cout << "Q_var = " << Q_var << std::endl;
	std::cout << "radius_var = " << radius_var << std::endl;
	std::cout << "poiseuille_var = " << poiseuille_var << std::endl;
	if(particle_deposition == 2)
		std::cout << "efficiency_var = " << efficiency_var << std::endl;
}







void NavierStokesCoupled::calculate_1d_boundary_values()
{
	
	std::cout << "yeah" << std::endl;

	if(sim_1d)
	{

		previous_flux_values_3d = flux_values_3d;
		previous_pressure_values_3d = pressure_values_3d;

		// if we are running a 3d-1d then the flux values are the first elements
		if(sim_type == 5 || sim_type == 3)
		{

			// there is no 0 boundary flux
			// i is the tree id, but we want to calculate on the boundary
			for(unsigned int i=1; i< flux_values_1d.size(); i++)
			{
				flux_values_1d[i] = ns_assembler->calculate_flux(tree_id_to_boundary_id[i]);
				pressure_values_1d[i] = ns_assembler->calculate_pressure(tree_id_to_boundary_id[i]);			
			}
		}
		else
		{

			// hmmm why is this not working
			flux_values_1d[0] = ns_assembler->calculate_flux(0);
			pressure_values_1d[0] = ns_assembler->calculate_pressure(0);

			// i is the tree id
			if(es->parameters.get<bool> ("calculate_1d_info_at_coupling_nodes"))
			{
				for(unsigned int i=1; i<subtree_starting_elements.size();i++)
				{
					std::cout << "subtree_starting_elements[i] = " << subtree_starting_elements[i] << std::endl;
					if(subtree_starting_elements[i] > -1)
					{
						flux_values_1d[i] = ns_assembler->calculate_flux(i,subtree_starting_elements[i]);
						pressure_values_1d[i] = ns_assembler->calculate_pressure(i,subtree_starting_elements[i]);
					}
				}
			}
		}

	}
	else
	{
		std::cout << "How on earth could I calculate boundary values of a 1D model"
							<< " if I'm not even running a 1D model? You, my friend, are a retard" << std::endl;
		exit_program = true;
		//exit(0);
	}
}







void NavierStokesCoupled::write_1d_solution()
{
	


	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(sim_type == 5)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns1d");
	}
	
	
	//std::cout << "before converting nodal-mono" << std::endl;
	//before we output take solution back to the monomial space
	convert_1d_nodal_to_monomial(*system->solution);
	system->solution->close();








  std::ostringstream file_name;

  // We write the file in the ExodusII format.
  //file_name << "results/out_1D_viscosity"
	//  				<< es->parameters.get<Real> ("viscosity")

	file_name << output_folder.str() << "out_1D";

	file_name	<< ".e-s."
            << std::setw(4)
            << std::setfill('0')
            << std::right
            << t_step;

	ExodusII_IO_Extended exo(mesh);

	exo.set_var_scalings(var_scalings_1D);

	std::vector<std::string> variables_1d;
	variables_1d.push_back("P");
	variables_1d.push_back("Q");
	variables_1d.push_back("radius");
	variables_1d.push_back("poiseuille");
	exo.set_output_variables(variables_1d);

	//std::cout << "before write disc" << std::endl;

	exo.write_discontinuous_exodusII (file_name.str(),
                                    *es);

	write_elem_pid_1d(exo);

	exo.write_time(1,time *time_scale_factor);

	
	std::cout << "EXODUSII output for timestep " << t_step
		<< " written to " << file_name.str() << std::endl;




	//before we assemble take the old solution and solution back to the nodal space
	convert_1d_monomial_to_nodal(*system->solution);	
	system->solution->close();

}





void NavierStokesCoupled::solve_1d_system()
{

	PetscLinearSolver<Number>* system_linear_solver =
			libmesh_cast_ptr<PetscLinearSolver<Number>* >
			(system_1d->linear_solver.get());

	std::string prefix_1d = "ns1d_";
	system_linear_solver->set_prefix(prefix_1d);
	system_linear_solver->init();	// set the name

	// set the poiseuille data
	set_poiseuille();


	// Assemble & solve the linear system.
	perf_log.push("Linear 1D solve");
	es->get_system("ns1d").solve();
	perf_log.pop("Linear 1D solve");


	//system_1d->matrix->print(std::cout,false);
	//system_1d->rhs->print(std::cout);

	//system_1d->solution->print();

  // Print out convergence information for the linear and
  // nonlinear iterations.

	std::cout << "Solver converged at step: "
        		<< system_1d->n_linear_iterations()
        		<< ", Solver residual: "
        		<< system_1d->final_linear_residual()
        		<< std::endl;
}


void NavierStokesCoupled::set_radii()
{

	//add the radius variable
	// i'm sure we only need to do this once!
	//create an equation system to hold data like the radius of the element

	
	const DofMap& dof_map = extra_1d_data_system->get_dof_map();

	MeshBase::const_element_iterator       el     =
		mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el =
		mesh.active_local_elements_end();

	std::vector<dof_id_type> dof_indices;
	const unsigned int rad_var = extra_1d_data_system->variable_number ("radius");

	for ( ; el != end_el; ++el)
	{
		const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
			dof_map.dof_indices (elem, dof_indices,rad_var);
			const dof_id_type elem_id = elem->id() - es->parameters.get<unsigned int>("n_initial_3d_elem");


			if(airway_data[elem_id].get_generation() < 5)
				std::cout << "radius at gen " << airway_data[elem_id].get_generation() << " = " << airway_data[elem_id].get_radius() << std::endl;
			
			for(unsigned int i=0; i < dof_indices.size(); i++)
				extra_1d_data_system->solution->set(dof_indices[i], airway_data[elem_id].get_radius());
		}
	}

	//must close the vector after editing it

	extra_1d_data_system->solution->close();
}


void NavierStokesCoupled::set_poiseuille()
{

	//add the radius variable
	// i'm sure we only need to do this once!
	//create an equation system to hold data like the radius of the element

	
	const DofMap& dof_map = extra_1d_data_system->get_dof_map();

	MeshBase::const_element_iterator       el     =
		mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el =
		mesh.active_local_elements_end();

	std::vector<dof_id_type> dof_indices;
	const unsigned int rad_var = extra_1d_data_system->variable_number ("poiseuille");

	for ( ; el != end_el; ++el)
	{
		const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
			dof_map.dof_indices (elem, dof_indices,rad_var);
			const dof_id_type elem_id = elem->id() - es->parameters.get<unsigned int>("n_initial_3d_elem");

			for(unsigned int i=0; i < dof_indices.size(); i++)
			{
				if(airway_data[elem_id].get_poiseuille())
					extra_1d_data_system->solution->set(dof_indices[i], 1);
				else
					extra_1d_data_system->solution->set(dof_indices[i], 0);
			}
		}
	}

	//must close the vector after editing it

	extra_1d_data_system->solution->close();
}







// convert the 1d part of a vector from monomial to nodal
void NavierStokesCoupled::convert_1d_monomial_to_nodal(NumericVector<Number>& vector)
{

	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(sim_type == 5)
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

  const DofMap & dof_map = system->get_dof_map();
  std::vector<dof_id_type> dof_indices_p;
  std::vector<dof_id_type> dof_indices_q;

	// not sure if this vector is local or global or which one i need to edit, probably the global one
	// remember the lo
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

								
  for ( ; el != end_el; ++el)
  {
		const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
			//element data object starts numbering from 0
			//and the values referenced in it also do so need to take this into account

			//const int current_el_idx = elem->id();

      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_q, q_var);

			double p0_mono = vector(dof_indices_p[0]);
			double p1_mono = vector(dof_indices_p[1]);
			double q0_mono = vector(dof_indices_q[0]);
			double q1_mono = vector(dof_indices_q[1]);

			vector.set(dof_indices_p[0],p0_mono - p1_mono);
			vector.set(dof_indices_p[1],p0_mono + p1_mono);
			vector.set(dof_indices_q[0],q0_mono - q1_mono);
			vector.set(dof_indices_q[1],q0_mono + q1_mono);
		}
	}
}

// convert the 1d part of a vector from nodal to monomial
void NavierStokesCoupled::convert_1d_nodal_to_monomial(NumericVector<Number>& vector)
{

	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(sim_type == 5)
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

  const DofMap & dof_map = system->get_dof_map();
  std::vector<dof_id_type> dof_indices_p;
  std::vector<dof_id_type> dof_indices_q;

	// not sure if this vector is local or global or which one i need to edit, probably the global one
	// remember the lo
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

								
  for ( ; el != end_el; ++el)
  {
		const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
			//element data object starts numbering from 0
			//and the values referenced in it also do so need to take this into account

			//const int current_el_idx = elem->id();

      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_q, q_var);

			double p0_nodal = vector(dof_indices_p[0]);
			double p1_nodal = vector(dof_indices_p[1]);
			double q0_nodal = vector(dof_indices_q[0]);
			double q1_nodal = vector(dof_indices_q[1]);

			vector.set(dof_indices_p[0],(p0_nodal + p1_nodal)/2.0);
			vector.set(dof_indices_p[1],(p1_nodal - p0_nodal)/2.0);
			vector.set(dof_indices_q[0],(q0_nodal + q1_nodal)/2.0);
			vector.set(dof_indices_q[1],(q1_nodal - q0_nodal)/2.0);
		}
	}
}


// ************************************************************************** //



// convert the 1d part of a vector from nodal to monomial
void NavierStokesCoupled::scale_1d_solution_vector(double flux_scaling=1.0,double pressure_scaling=1.0)
{

	TransientLinearImplicitSystem * system;
  // Get a reference to the Stokes system object.
	if(sim_type == 5)
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
  const int P_var = system->variable_number ("P");
  const int Q_var = system->variable_number ("Q");

	double value = 0.0;

	// need to make the solution vector global
  std::vector<Number> sys_soln;
  system->update_global_solution (sys_soln);


	for(unsigned int i=0; i<dof_variable_type_1d.size(); i++)
	{
		if(dof_variable_type_1d[i] == Q_var)
		{
			value = sys_soln[i];
			system->solution->set(i,flux_scaling * value);
		}
		else if(dof_variable_type_1d[i] == P_var)
		{
			value = sys_soln[i];
			system->solution->set(i,pressure_scaling * value);
		}			
	}
	

	for(unsigned int i=0; i<dof_variable_type_coupled.size(); i++)
	{
		if(dof_variable_type_coupled[i] == Q_var)
		{
			value = sys_soln[i];
			system->solution->set(i,flux_scaling * value);
		}
		else if(dof_variable_type_coupled[i] == P_var)
		{
			value = sys_soln[i];
			system->solution->set(i,pressure_scaling * value);
		}
	}


	system->solution->close();
}


// okay so we want to output the poiseuille per generation, but also we want to 
// output the data on for each airway so that we can supplement it with 3d data
// later possibly.
void NavierStokesCoupled::output_poiseuille_resistance_per_generation()
{

	std::cout<< "hey there" << std::endl;
	
	bool per_order = false;
	
	TransientLinearImplicitSystem * system;
	TransientLinearImplicitSystem * system_threed;
  // Get a reference to the Stokes system object.
	if(sim_type == 5)
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
		system_threed =
		  &es->get_system<TransientLinearImplicitSystem> ("ns3d1d");
	}
	else
	{
		system =
		  &es->get_system<TransientLinearImplicitSystem> ("ns1d");

		if(sim_3d)
			system_threed =
			  &es->get_system<TransientLinearImplicitSystem> ("ns3d");
	}


	unsigned int n_initial_3d_elem = es->parameters.get<unsigned int>("n_initial_3d_elem");

  //const double viscosity = es->parameters.get<double>("viscosity_1d");
  const unsigned int p_var = system->variable_number ("P");
  const unsigned int q_var = system->variable_number ("Q");


  unsigned int p_var_threed = 0;
	if(sim_3d)
	  p_var_threed = system_threed->variable_number ("p");
	
  const DofMap & dof_map = system->get_dof_map();
  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_p;
  std::vector<dof_id_type> dof_indices_q;


	// trees are gonna be numbered by their subdomain_id
	unsigned int num_1d_trees = es->parameters.get<unsigned int>("num_1d_trees");

	for(unsigned int i=0; i<num_generations.size(); i++)
	{
		std::cout << "num_generations = " << num_generations[i] << std::endl;
	}


	//vector of edges per generation to calc average
	std::vector<unsigned int> total_num_edges_per_generation;
	std::vector<std::vector<unsigned int> > num_edges_per_generation(num_1d_trees);
	for(unsigned int i=0; i<num_1d_trees; i++)
		num_edges_per_generation[i] = std::vector<unsigned int> (num_generations[i],0);

	//vector of resistance per generation
	std::vector<double> total_resistance_per_generation;
	std::vector<std::vector<double> > resistance_per_generation(num_1d_trees);
	for(unsigned int i=0; i<num_1d_trees; i++)
		resistance_per_generation[i] = std::vector<double> (num_generations[i],0.);

	//vector of average_pressure_diff per generation
	std::vector<std::vector<double> > average_pressure_diff_per_generation(num_1d_trees);
	for(unsigned int i=0; i<num_1d_trees; i++)
		average_pressure_diff_per_generation[i] = std::vector<double> (num_generations[i],0.);

	std::cout << "num_generations.size() = " << num_edges_per_generation[0].size() << std::endl;

	std::cout << "num_generations = " << num_generations[0] << std::endl;

	std::cout<< "hey there" << std::endl;

	//vector of resistance per generation, not sized
	std::vector<std::vector<double> > resistance_per_order(num_1d_trees);
	std::vector<std::vector<double> > num_edges_per_order(num_1d_trees);

	double flow_scale = es->parameters.get<double>("velocity_scale") * 
										pow(es->parameters.get<double>("length_scale"),2.0);
	double mean_pressure_scale = es->parameters.get<double>("density") *
													pow(es->parameters.get<double>("velocity_scale"),2.0);

	unsigned int max_generation = 0;


	// ************* CALCULATE RESISTANCE AND PRESSURE/FLUX IN 1D TREE ******** //

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

	// note this is not a parallel loop cause i am lazy
  for ( ; el != end_el; ++el)
  {
		const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
			//element data object starts numbering from 0
			//and the values referenced in it also do so need to take this into account

			const int current_el_idx = elem->id();
			unsigned int current_1d_el_idx = current_el_idx -	n_initial_3d_elem;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_q, q_var);

			int generation = airway_data[current_1d_el_idx].get_generation();
			unsigned int order = airway_data[current_1d_el_idx].get_order();
			//if(order < 1)
			//{
				//std::cout << "order = " << order << std::endl;
			//	std::exit(0);
			//}

			if(generation > max_generation)
				max_generation = generation;

			double p0 = system->current_solution(dof_indices_p[0]);
			double p1 = system->current_solution(dof_indices_p[1]);
			double q0 = system->current_solution(dof_indices_q[0]);
			double q1 = system->current_solution(dof_indices_q[1]);
		
			// should be correct either way...
			//double pressure_diff = fabs((p0 + p1)/2 - (p0 - p1)/2);
			//double flow = fabs((q0 + q1)/2 + (q0 - q1)/2);

			double pressure_diff = fabs(p0 - p1);
			// convert pressure to SI units
			pressure_diff *= mean_pressure_scale;
			
			double flow = fabs(q0 + q1)/2.;
			flow *= flow_scale;

			double resistance = pressure_diff / flow;
			//get average pressure drop but, total flow

			airway_data[current_1d_el_idx].set_flow(flow);
			airway_data[current_1d_el_idx].set_pressure_diff(pressure_diff);

			
			//add to the correct tree
			if(!per_order)
			{
				num_edges_per_generation[elem->subdomain_id() - subdomains_1d.front()][generation] += 1;
				resistance_per_generation[elem->subdomain_id() - subdomains_1d.front()][generation] += resistance;
			}
			else
			{
				if(order > num_edges_per_order[elem->subdomain_id() - subdomains_1d.front()].size())
				{
					num_edges_per_order[elem->subdomain_id() - subdomains_1d.front()].resize(order);
					resistance_per_order[elem->subdomain_id() - subdomains_1d.front()].resize(order);
				}
				
				num_edges_per_order[elem->subdomain_id() - subdomains_1d.front()][order - 1] +=1;
				resistance_per_order[elem->subdomain_id() - subdomains_1d.front()][order - 1] += resistance;
			}
			

    }// end of subdomain conditional
	}// end of element loop


	// get average
	total_num_edges_per_generation.resize(max_generation+1,0);
	total_resistance_per_generation.resize(max_generation+1,0);

	//do a loop over the 3d section
	//we add this to the total

	// need outside cause may want to output
	std::vector<double> flux_per_3d_element(centreline_airway_data.size());

	// ************* CALCULATE RESISTANCE AND PRESSURE/FLUX IN 3D TREE ******** //
	
	if(es->parameters.get<bool>("use_centreline_data"))
	{


		// ******* MOVE POINTS SO THAT INSIDE DOMAIN ************************** //

		// percentage of airway moved to get inside 3D airway
		double tol = 1e-3;
		for(unsigned int i=0; i<centreline_airway_data.size(); i++)
		{
			// hmmm, so if we are on a boundary we may have difficulty...
			double segment_length = (centreline_points[i][1] - centreline_points[i][0]).size();

			// if we don't have a parent then move point a small distance
			if(!centreline_airway_data[i].has_parent())
			{
				unsigned int boundary_id = 0;
				Point normal = surface_boundaries[0]->get_normal();
				// move point backwards along normal
				Point adjusted_point = centreline_points[i][0] - tol*segment_length*normal;
				std::cout << "calculating pressure at adjusted point " << adjusted_point << std::endl;
				centreline_points[i][0] = adjusted_point;
			}

			if(!centreline_airway_data[i].has_daughter_1())
			{
				unsigned int boundary_id = tree_id_to_boundary_id[centreline_terminal_id_to_tree_id[i]];
				Point normal = surface_boundaries[boundary_id]->get_normal();
				// move point backwards along normal
				Point adjusted_point = centreline_points[i][1] - tol*segment_length*normal;
				std::cout << "calculating pressure at adjusted point " << adjusted_point << std::endl;
				centreline_points[i][1] = adjusted_point;
			}

		}

		// ******* FIND WHAT ELEMENT END POINTS ARE IN ************************** //

		std::vector<std::vector<unsigned int> > centreline_points_elements(centreline_points.size(),std::vector<unsigned int>(2));
		
		el     = mesh.active_elements_begin();
		const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

		unsigned int count = 0;
		// note this is not a parallel loop cause i am lazy
		for ( ; el != end_el; ++el)
		{
			const Elem* elem = *el;
			if(std::find(subdomains_3d.begin(), subdomains_3d.end(), elem->subdomain_id()) != subdomains_3d.end())
			{
				for(unsigned int i=0; i<centreline_points.size(); i++)
				{
					//check point 1
					if(elem->contains_point(centreline_points[i][0]))
					{
						std::cout << "point 0 " << i << ": " << centreline_points[i][0] << " found in elem " << elem->id() << std::endl;
						centreline_points_elements[i][0] = elem->id();
						count++;
					}

					//check point 1
					if(elem->contains_point(centreline_points[i][1]))
					{
						std::cout << "point 1 " << i << ": " << centreline_points[i][1] <<  " found in elem " << elem->id() << std::endl;
						centreline_points_elements[i][1] = elem->id();
						count++;
					}
				}
			}

		}

		std::cout << "num_points found = " << count << std::endl;


		
		// 1.) calculate the fluxes by going up from each terminal branch
		// 2.) calculate the pressures, only want to do this once per point, 
		// but it's okay only double the calculation

		for(unsigned int i=0; i<centreline_airway_data.size(); i++)
		{
			std::cout << "1" << std::endl;
			//loop over terminal elements, terminal element when has no daughter 1
			//std::cout << "centreline_airway_data[" <<  i << "].get_daughter_1() = " << centreline_airway_data[i].get_daughter_1() << std::endl;
			if(!centreline_airway_data[i].has_daughter_1())
			{
				// okay, so we need to know what 1d elements this goes into, because
				double flux_in_terminal_element = flux_values_3d[tree_id_to_boundary_id[centreline_terminal_id_to_tree_id[i]]];
				flux_in_terminal_element *= flow_scale;

			std::cout << "glux = " << flux_in_terminal_element << std::endl;
				// generation of this element
				unsigned int generation = centreline_airway_data[i].get_generation();

				//go up all the generations
				int next_element = i;
				std::cout << "generation = " << generation <<std::endl;
				for(int j=generation; j>-1; j--)
				{
					std::cout << "got here" << std::endl;
					//add the flux to the
					//flux_per_3d_element[next_element] += flux_in_terminal_element;
					centreline_airway_data[next_element].set_flow(centreline_airway_data[next_element].get_flow() + flux_in_terminal_element);
					if(centreline_airway_data[next_element].has_parent())
						next_element = centreline_airway_data[next_element].get_parent();	// make the next element the parent
					//std::cout << "next_element = " << next_element << std::endl;
				}
			}

			double start_pressure = system_threed->point_value(p_var_threed,centreline_points[i][0],*mesh.elem(centreline_points_elements[i][0]));
			double end_pressure = system_threed->point_value(p_var_threed,centreline_points[i][1],*mesh.elem(centreline_points_elements[i][1]));

			double pressure_diff = start_pressure - end_pressure;

			// convert pressure to SI units
			pressure_diff *= mean_pressure_scale;

			centreline_airway_data[i].set_pressure_diff(pressure_diff);
		}
		
		
		// calculate the resistance and add to the vector

		for(unsigned int i=0; i<centreline_airway_data.size(); i++)
		{
			// generation of this element
			unsigned int generation = centreline_airway_data[i].get_generation();

			std::cout << "1" << std::endl;
			//calculate the resistance
			std::cout << "pressure_diff_per_3d_element[i] = " << centreline_airway_data[i].get_pressure_diff() << std::endl;
			std::cout << "flux_per_3d_element[i] = " << centreline_airway_data[i].get_flow() << std::endl;
			double resistance = centreline_airway_data[i].get_pressure_diff() / centreline_airway_data[i].get_flow();

			//add to the correct tree
			if(!per_order)
			{
				total_num_edges_per_generation[generation] += 1;
				total_resistance_per_generation[generation] += resistance;
			}
			else
			{
				/*
				if(order > num_edges_per_order[elem->subdomain_id() - subdomains_1d.front()].size())
				{
					num_edges_per_order[elem->subdomain_id() - subdomains_1d.front()].resize(order);
					resistance_per_order[elem->subdomain_id() - subdomains_1d.front()].resize(order);
				}
				
				num_edges_per_order[elem->subdomain_id() - subdomains_1d.front()][order - 1] +=1;
				resistance_per_order[elem->subdomain_id() - subdomains_1d.front()][order - 1] += resistance;
				*/


			}
			
			std::cout << "2" << std::endl;
		}
		
		
	}

	std::cout<< "hey there, max_gerertaion = " << max_generation << std::endl;

		


	// ********** COMPILING RESISTANCES TOGETHER ***************************** //

	
	// loop over all the trees
	for(unsigned int i=0; i<num_edges_per_generation.size(); i++)
	{
		if(!per_order)
		{
			// loop over all the generations of each tree
			for(unsigned int j=0; j<num_edges_per_generation[i].size(); j++)
			{
				if(num_edges_per_generation[i][j] > 0)
				{
					//average_pressure_diff_per_generation[i][j] /= num_edges_per_generation[i][j];
				
					total_num_edges_per_generation[j] += num_edges_per_generation[i][j];
					total_resistance_per_generation[j] += resistance_per_generation[i][j];
					if(es->parameters.get<bool>("assume_symmetric_tree"))
					{
						// assume symmetric tree and calc parallel resistance
						resistance_per_generation[i][j] = resistance_per_generation[i][j]/pow(2.,j);
					}
					else
					{
						//find average resistance
						resistance_per_generation[i][j] = resistance_per_generation[i][j]/num_edges_per_generation[i][j];
					}
					//resistance_per_generation[i][j] = 1./resistance_per_generation[i][j];

				}
			}
		}
		else
		{
			for(unsigned int j=0; j<num_edges_per_order[i].size(); j++)
			{
				if(num_edges_per_order[i][j] > 0)
				{
					//average_pressure_diff_per_generation[i][j] /= num_edges_per_generation[i][j];
					resistance_per_order[i][j] = 1./resistance_per_order[i][j];
					//resistance_per_generation[i][j] = 1./resistance_per_generation[i][j];
				}
			}

		}
	}

	// calculate the total average resistance 
	for(unsigned int i=0; i<total_num_edges_per_generation.size(); i++)
	{
		if(!per_order)
		{
			if(es->parameters.get<bool>("assume_symmetric_tree"))
			{
				// assume symmetric tree and calc parallel resistance
				total_resistance_per_generation[i] = total_resistance_per_generation[i]/pow(2.,i);
			}
			else
			{
				//find average resistance
				total_resistance_per_generation[i] = total_resistance_per_generation[i]/total_num_edges_per_generation[i];
			}
		}
	}

	




	// ****** OUTPUTTING RESISTANCE AND PRESSURE/FLUX DATA ***************** //
	if(num_1d_trees > 0)
	{

		std::cout << "boo" << std::endl;
		std::ostringstream output_data_file_name;
		std::ofstream resistance_output_file;

		std::cout << "boo" << std::endl;
		output_data_file_name << output_folder.str() << "resistance_per_generation.dat";
		std::cout << "boo " << output_data_file_name.str() << std::endl;
		resistance_output_file.open(output_data_file_name.str().c_str());

		

		std::cout<< "boo" << std::endl;
		resistance_output_file << "# Resistance per generation per tree" << std::endl;
		// write geometry variables later when reading meta data file
		// write boundary conditions later with Picard class
		resistance_output_file << "# Results" << std::endl;
		if(!per_order)
			resistance_output_file << "# generation";
		else
			resistance_output_file << "# order";
		
		std::cout<< "boo" << std::endl;
		if(es->parameters.get<bool>("output_resistance_for_each_tree"))
		{
			for(unsigned int i=0; i < num_1d_trees; i++)
			{
				resistance_output_file <<	"\ttree_" << i+1;
				resistance_output_file <<	"\tnum_edges";
			}
		}
		//resistance_output_file <<	"\ttree_" << 1;
		resistance_output_file <<	"\ttotal";
		resistance_output_file <<	"\ttotal_num_edges";

		//unsigned int max_generations = *std::max_element(num_generations.begin(), num_generations.end());

		std::cout << "max_generations = " << max_generation << std::endl;
		std::cout << "resistance_per_generation.size() = " << resistance_per_generation[0].size() << std::endl;
		std::cout << "num_edges_per_generation.size() = " << num_edges_per_generation[0].size() << std::endl;

		if(!per_order)
		{
			for(unsigned int i=0; i < max_generation; i++)
			{
				resistance_output_file << std::endl;
				resistance_output_file <<	i;
				//for(unsigned int j=0; j < num_1d_trees; j++)

				if(es->parameters.get<bool>("output_resistance_for_each_tree"))
				{
					for(unsigned int j=0; j < num_1d_trees; j++)
					{
						if(i < num_generations[j])
						{
							resistance_output_file <<	"\t" << resistance_per_generation[j][i];
							resistance_output_file <<	"\t" << num_edges_per_generation[j][i];
						}
						else
						{
							resistance_output_file <<	"\t" << 0.;
							resistance_output_file <<	"\t" << 0.;
						}

					}
				}

				resistance_output_file <<	"\t" << total_resistance_per_generation[i];
				resistance_output_file <<	"\t" << total_num_edges_per_generation[i];
				
						
			}
		}
		else
		{
			unsigned int max_num_order = 0;

			for(unsigned int i=0; i<num_edges_per_order.size(); i++)
			{
				if(num_edges_per_order[i].size() > max_num_order)
					max_num_order = num_edges_per_order[i].size();
			}

			for(unsigned int i=0; i < max_num_order; i++)
			{
				resistance_output_file << std::endl;
				resistance_output_file <<	i+1;
				for(unsigned int j=0; j < num_1d_trees; j++)
				{
					if(num_edges_per_order[j][i] > 0)
						resistance_output_file <<	"\t" << resistance_per_order[j][i];
					else
						resistance_output_file <<	"\t" << 0.;
				}		
			}

		}
		

		resistance_output_file.close();

		// ************* THE FILE CONTAINING FLUX AND PRESSURE FOR EACH AIRWAY ******* //

		std::ostringstream flux_data_file_name;
		std::ofstream flux_output_file;
		flux_data_file_name << output_folder.str() << "flux_and_pressure_per_airway.dat";
		flux_output_file.open(flux_data_file_name.str().c_str());

	

		// write geometry variables later when reading meta data file
		// write boundary conditions later with Picard class
		flux_output_file << "# local_edge_num";
		flux_output_file << "\ttree_no";
		flux_output_file << "\tgeneration";
		flux_output_file << "\torder";
		flux_output_file << "\tpressure_diff";
		flux_output_file << "\tflux";
	
		if(es->parameters.get<bool>("use_centreline_data"))
		{

			for(unsigned int i=0; i < centreline_airway_data.size(); i++)
			{
				flux_output_file << std::endl;
				flux_output_file <<	centreline_airway_data[i].get_local_elem_number();
				flux_output_file <<	"\t";
				flux_output_file <<	centreline_airway_data[i].get_tree_number();
				flux_output_file <<	"\t";
				flux_output_file <<	centreline_airway_data[i].get_generation();
				flux_output_file <<	"\t";
				flux_output_file <<	centreline_airway_data[i].get_order();
				flux_output_file <<	"\t";				
				flux_output_file <<	centreline_airway_data[i].get_pressure_diff();
				flux_output_file <<	"\t";
				flux_output_file <<	centreline_airway_data[i].get_flow();
			}
		}

		for(unsigned int i=0; i < airway_data.size(); i++)
		{
			flux_output_file << std::endl;
			flux_output_file <<	airway_data[i].get_local_elem_number();
			flux_output_file <<	"\t";
			flux_output_file <<	airway_data[i].get_tree_number();
			flux_output_file <<	"\t";
			flux_output_file <<	airway_data[i].get_generation();
			flux_output_file <<	"\t";
			flux_output_file <<	airway_data[i].get_order();
			flux_output_file <<	"\t";				
			flux_output_file <<	airway_data[i].get_pressure_diff();
			flux_output_file <<	"\t";
			flux_output_file <<	airway_data[i].get_flow();
		}
		
		flux_output_file.close();

	}
	else
	{
		std::cout << "WARNING: no 1d trees so no need to output num trees per generation/order" << std::endl;
	}
	

	std::cout << "end outputing"<< std::endl;
}


// *********** PLOT THE ERROR ASSOCIATED WITH THE 3D SOLUTION *************** //

void NavierStokesCoupled::write_elem_pid_1d(ExodusII_IO_Extended& io)
{

  AutoPtr<MeshBase> meshptr = mesh.clone();
  MeshBase &temp_mesh = *meshptr;
	temp_mesh.partitioner()->set_custom_partitioning(es->parameters.set<bool>("custom_partitioning"));
  temp_mesh.all_first_order();
  EquationSystems temp_es (temp_mesh);
  ExplicitSystem& processor_system
    = temp_es.add_system<ExplicitSystem> ("Processor");
  processor_system.add_variable("proc_id", CONSTANT, MONOMIAL);
  temp_es.init();

  const DofMap& processor_dof_map = processor_system.get_dof_map();

  MeshBase::const_element_iterator       el     =
    temp_mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    temp_mesh.active_local_elements_end();
  std::vector<dof_id_type> dof_indices;

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{

		  processor_dof_map.dof_indices(elem, dof_indices);

		  //0 for the monomial basis
		  const dof_id_type solution_index = dof_indices[0];

		  // We may have zero error values in special circumstances
		  // libmesh_assert_greater ((*this)[elem_id], 0.);
		  processor_system.solution->set(solution_index, elem->processor_id());
		}
  }

  // We may have to renumber if the original numbering was not
  // contiguous.  Since this is just a temporary mesh, that's probably
  // fine.
  if (temp_mesh.max_elem_id() != temp_mesh.n_elem() ||
      temp_mesh.max_node_id() != temp_mesh.n_nodes())
  {
    temp_mesh.allow_renumbering(true);
    temp_mesh.renumber_nodes_and_elements();
  }

  //ExodusII_IO io(mesh);
  //io.write(filename);
  //io.write_element_data(temp_es,time);
  io.write_element_data(temp_es);

}


void NavierStokesCoupled::write_efficiency_solution()
{
	


	ExplicitSystem * system;
  // Get a reference to the Stokes system object.
	system =
	  &es->get_system<ExplicitSystem> ("Particle-Deposition-1D");
	
	//put the efficiency data into the system






  std::ostringstream file_name;

  // We write the file in the ExodusII format.
  //file_name << "results/out_1D_viscosity"
	//  				<< es->parameters.get<Real> ("viscosity")

	file_name << output_folder.str() << "out_1D_efficiency";

	file_name	<< ".e";

	ExodusII_IO_Extended exo(mesh);

	exo.set_var_scalings(var_scalings_1D);

	std::vector<std::string> variables_1d;
	variables_1d.push_back("efficiency");
	exo.set_output_variables(variables_1d);

	//std::cout << "before write disc" << std::endl;

	exo.write_discontinuous_exodusII (file_name.str(),
                                    *es);

	write_elem_pid_1d(exo);

	
	std::cout << "EXODUSII output for particle deposition Hofmann: written to " << file_name.str() << std::endl;


}


void NavierStokesCoupled::set_efficiency()
{

	//add the radius variable
	// i'm sure we only need to do this once!
	//create an equation system to hold data like the radius of the element

	//get the data from the hofmann particle deposition object
	
	const DofMap& dof_map = particle_deposition_system_1d->get_dof_map();

	MeshBase::const_element_iterator       el     =
		mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el =
		mesh.active_local_elements_end();

	std::vector<dof_id_type> dof_indices;
	const unsigned int eff_var = particle_deposition_system_1d->variable_number ("efficiency");

	for ( ; el != end_el; ++el)
	{
		const Elem* elem = *el;
		if(std::find(subdomains_1d.begin(), subdomains_1d.end(), elem->subdomain_id()) != subdomains_1d.end())
		{
			dof_map.dof_indices (elem, dof_indices,eff_var);
			const dof_id_type elem_id = elem->id() - es->parameters.get<unsigned int>("n_initial_3d_elem");
			
			for(unsigned int i=0; i < dof_indices.size(); i++)
				particle_deposition_system_1d->solution->set(dof_indices[i], total_efficiency[elem_id]);
		}
	}

	//must close the vector after editing it

	particle_deposition_system_1d->solution->close();
}




