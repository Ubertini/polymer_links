#include <fstream>
#include <string>
#include <iostream>
#include "monte_carlo_ring_links.hpp"
#include "initial_condition.hpp"


int main(){
	
	double monomer_density=0.4;
	long long int MC_steps=1e6;
	int stride = 10000000;
	std::string length= std::to_string(int(monomer_density*100));
	int polymer_length=1280;
	int number_of_polymers=160;
	int run = 1;
	std::vector<int> coordinates;
		
	std::ifstream file("ICS_RINGS/IC_RING_"+length+"_"+std::to_string(polymer_length)); 
		
	int i=0;
	double val=0;
		
	if ( file.is_open()){
		while(!file.eof()){
	     while(file >> val){
	      coordinates.push_back(int(val));
	      i+=1;
	     }
	  }
	}
	  
	else{
	  	std::cout << "unable to open file."<<std::endl;
	}

	file.close();
	  

	MC_routine_rings_links(coordinates,polymer_length,number_of_polymers,monomer_density,MC_steps,stride,run);


	
	return 0;
}
