#ifndef __initial_condition__
#define __initial_condition__

#include <fstream>
#include <string>
#include <iostream>
#include "monte_carlo_linear.hpp"

void Initial_Condition(int polymer_length,int number_of_polymers){
		
		std::cout<<polymer_length<<"\t"<<number_of_polymers<<std::endl;
		double monomer_density= double(polymer_length*number_of_polymers)/double(Lattice_size*Lattice_size*Lattice_size);
		std::cout<<monomer_density<<std::endl;
		std::string length= std::to_string(int(monomer_density*100));
		std::cout<<length<<std::endl;

		int lattice[Lattice_size][Lattice_size][Lattice_size]={0};
		
		std::vector<MC_move> mc_move;
		mc_move.push_back(MC_move(1,0,0));
		mc_move.push_back(MC_move(-1,0,0));
		mc_move.push_back(MC_move(0,1,0));
		mc_move.push_back(MC_move(0,-1,0));
		mc_move.push_back(MC_move(0,0,1));
		mc_move.push_back(MC_move(0,0,-1));
		mc_move.push_back(MC_move(1,0,1));
		mc_move.push_back(MC_move(-1,0,-1));
		mc_move.push_back(MC_move(0,1,1));
		mc_move.push_back(MC_move(0,-1,-1));
		mc_move.push_back(MC_move(1,1,1));
		mc_move.push_back(MC_move(-1,-1,-1));
		
		int number_of_created_polymer=0;
		
		std::vector<polymer> melt;

		while(number_of_created_polymer < number_of_polymers){
			
			polymer polymer_tmp{polymer_length,number_of_created_polymer};
			
			int t= rand() % Lattice_size;
			int u= rand() % Lattice_size;
			int v= rand() % Lattice_size;
			
			while(lattice[t][u][v]!=0){
				t= rand() % Lattice_size;
				u= rand() % Lattice_size;
				v= rand() % Lattice_size;
			}
			
			monomer monomer_tmp{t,u,v};
			polymer_tmp.chain.push_back(monomer_tmp);
			lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]+=1;
			int n_monomers=1;
			int occupation=0;

			for(int j=1;j<int(polymer_length/2);j++){
				
				int n_trials=1;
				auto MC_move = rand()%12; //there are 12 possible movement
				monomer monomer_tmp{t+mc_move[MC_move].dt,u+mc_move[MC_move].du,v+mc_move[MC_move].dv};
				occupation = lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image];

				while(occupation != 0 && n_trials < 12){
					auto MC_move = rand()%12; //there are 12 possible movement
					monomer monomer_trial{t+mc_move[MC_move].dt,u+mc_move[MC_move].du,v+mc_move[MC_move].dv};
					monomer_tmp = monomer_trial;
					//std::cout<<monomer_tmp.t<<"\t"<<monomer_trial.t<<"\t"<<monomer_tmp.u<<"\t"<<monomer_trial.u<<"\t";
					occupation = lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image];
					n_trials+=1;
					//std::cout<<"n_trials="<<n_trials<<"\t"<<occupation<<"\t"<<lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]<<std::endl;
				}

				if(n_trials <= 12 && occupation ==0){
					t=monomer_tmp.t;
					u=monomer_tmp.u;
					v=monomer_tmp.v;
					polymer_tmp.chain.push_back(monomer_tmp);
					//std::cout<<"occupation="<<occupation<<"  lattice = "<<lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]<<std::endl;
					lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]+=1;
					n_monomers+=1;
				}

				else{
					
					for(std::size_t i=0;i<polymer_tmp.chain.size();i++){
						lattice[polymer_tmp.chain[i].t_image][polymer_tmp.chain[i].u_image][polymer_tmp.chain[i].v_image]-=1;
					}

					//std::cout<<"break at n = "<<number_of_created_polymer<<std::endl;
					break;
				}
			
			}

			if(n_monomers==int(polymer_length/2)){
				for(int i=n_monomers-1;i>=0;i--){
					polymer_tmp.chain.push_back(polymer_tmp.chain[i]);
					n_monomers+=1;
					lattice[polymer_tmp.chain[i].t_image][polymer_tmp.chain[i].u_image][polymer_tmp.chain[i].v_image]+=1;
				}
				//std::cout<<n_monomers<<std::endl;
				melt.push_back(polymer_tmp);
				number_of_created_polymer+=1;
				//std::cout<<"number_of_polymers="<<number_of_created_polymer<<std::endl;
			}
		}
		
		bool condition=true;

		for(auto i=0;i<Lattice_size;i++){
			for(auto j=0;j<Lattice_size;j++){
				for(auto k=0;k<Lattice_size;k++){
					if(lattice[i][j][k] != 2 && lattice[i][j][k] != 0){
						condition=false;
						std::cout<<lattice[i][j][k]<<std::endl;
						std::cout<<i<<"\t"<<j<<"\t"<<k<<std::endl;
					}
				}
			}
		}

		if(condition){
			for(std::size_t k=0;k<melt.size();k++){
				for(std::size_t z=0;z<melt[k].chain.size()-1;z++){
					if(fabs(distance(melt[k].chain[z],melt[k].chain[z+1])-1)>0.1 && z != int(melt[k].chain.size()/2)-1){
						condition=false;
						std::cout<<distance(melt[k].chain[z],melt[k].chain[z+1])<<std::endl;
					}
				}
			}
		}

		if(condition){
			std::string name_file = "ICS/IC_"+length+"_"+std::to_string(polymer_length);
			std::cout<<name_file<<std::endl;
			std::ofstream traj;
			traj.open(name_file);

			for(std::size_t k=0;k<melt.size();k++){
				for(std::size_t z=0;z<melt[k].chain.size();z++){
					traj<<melt[k].chain[z].t<<"\t"<<melt[k].chain[z].u<<"\t"<<melt[k].chain[z].v<<"\t";
				}
			}
			traj.close();
		}

		else{
			std::cout<<"There's been an error somewhere during the creation of the initial condition"<<std::endl;
			std::string name_file = "IC_"+length+"_"+std::to_string(polymer_length)+"_error";
			std::ofstream traj;
			traj.open(name_file);
			traj<<"There's been an error somewhere during the creation of the initial condition"<<std::endl;
		}
}

void Initial_Condition_Ring(int polymer_length,int number_of_polymers){

		std::cout<<polymer_length<<"\t"<<number_of_polymers<<std::endl;


		double monomer_density= double(polymer_length*number_of_polymers)/double(Lattice_size*Lattice_size*Lattice_size);
		std::cout<<monomer_density<<std::endl;
		std::string length= std::to_string(int(monomer_density*100));

		int lattice[Lattice_size][Lattice_size][Lattice_size]={0};
		
		std::vector<MC_move> mc_move;
		mc_move.push_back(MC_move(1,0,0));
		mc_move.push_back(MC_move(-1,0,0));
		mc_move.push_back(MC_move(0,1,0));
		mc_move.push_back(MC_move(0,-1,0));
		mc_move.push_back(MC_move(0,0,1));
		mc_move.push_back(MC_move(0,0,-1));
		mc_move.push_back(MC_move(1,0,1));
		mc_move.push_back(MC_move(-1,0,-1));
		mc_move.push_back(MC_move(0,1,1));
		mc_move.push_back(MC_move(0,-1,-1));
		mc_move.push_back(MC_move(1,1,1));
		mc_move.push_back(MC_move(-1,-1,-1));
		
		int number_of_created_polymer=0;
		
		std::vector<polymer> melt;

		while(number_of_created_polymer < number_of_polymers){
			
			polymer polymer_tmp{polymer_length,number_of_created_polymer};
			
			int t= rand() % Lattice_size;
			int u= rand() % Lattice_size;
			int v= rand() % Lattice_size;
			
			while(lattice[t][u][v]!=0){
				t= rand() % Lattice_size;
				u= rand() % Lattice_size;
				v= rand() % Lattice_size;
			}
			
			monomer monomer_tmp{t,u,v};
			polymer_tmp.chain.push_back(monomer_tmp);
			lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]+=1;
			int n_monomers=1;
			int occupation=0;

			for(int j=1;j<int(polymer_length/2);j++){
				
				int n_trials=1;
				auto MC_move = rand()%12; //there are 12 possible movement
				monomer monomer_tmp{t+mc_move[MC_move].dt,u+mc_move[MC_move].du,v+mc_move[MC_move].dv};
				occupation = lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image];

				while(occupation != 0 && n_trials < 12){
					auto MC_move = rand()%12; //there are 12 possible movement
					monomer monomer_trial{t+mc_move[MC_move].dt,u+mc_move[MC_move].du,v+mc_move[MC_move].dv};
					monomer_tmp = monomer_trial;
					//std::cout<<monomer_tmp.t<<"\t"<<monomer_trial.t<<"\t"<<monomer_tmp.u<<"\t"<<monomer_trial.u<<"\t";
					occupation = lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image];
					n_trials+=1;
					//std::cout<<"n_trials="<<n_trials<<"\t"<<occupation<<"\t"<<lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]<<std::endl;
				}

				if(n_trials <= 12 && occupation ==0){
					t=monomer_tmp.t;
					u=monomer_tmp.u;
					v=monomer_tmp.v;
					polymer_tmp.chain.push_back(monomer_tmp);
					//std::cout<<"occupation="<<occupation<<"  lattice = "<<lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]<<std::endl;
					lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]+=1;
					n_monomers+=1;
				}

				else{
					
					for(std::size_t i=0;i<polymer_tmp.chain.size();i++){
						lattice[polymer_tmp.chain[i].t_image][polymer_tmp.chain[i].u_image][polymer_tmp.chain[i].v_image]-=1;
					}

					//std::cout<<"break at n = "<<number_of_created_polymer<<std::endl;
					break;
				}
			
			}

			if(n_monomers==int(polymer_length/2)){
				for(int i=n_monomers-1;i>=0;i--){
					polymer_tmp.chain.push_back(polymer_tmp.chain[i]);
					n_monomers+=1;
					lattice[polymer_tmp.chain[i].t_image][polymer_tmp.chain[i].u_image][polymer_tmp.chain[i].v_image]+=1;
				}
				//std::cout<<n_monomers<<std::endl;
				melt.push_back(polymer_tmp);
				number_of_created_polymer+=1;
				//std::cout<<"number_of_polymers="<<number_of_created_polymer<<std::endl;
			}
		}
		
		bool condition=true;

		for(auto i=0;i<Lattice_size;i++){
			for(auto j=0;j<Lattice_size;j++){
				for(auto k=0;k<Lattice_size;k++){
					if(lattice[i][j][k] != 2 && lattice[i][j][k] != 0){
						condition=false;
						std::cout<<lattice[i][j][k]<<std::endl;
						std::cout<<i<<"\t"<<j<<"\t"<<k<<std::endl;
					}
				}
			}
		}

		if(condition){
			for(std::size_t k=0;k<melt.size();k++){
				for(std::size_t z=0;z<melt[k].chain.size()-1;z++){
					if(fabs(distance(melt[k].chain[z],melt[k].chain[z+1])-1)>0.1 && z != int(melt[k].chain.size()/2)-1){
						condition=false;
						std::cout<<distance(melt[k].chain[z],melt[k].chain[z+1])<<std::endl;
					}
				}
			}
		}

		if(condition){
			std::string name_file = "ICS_RINGS/IC_RING_"+length+"_"+std::to_string(polymer_length);
			std::cout<<name_file<<std::endl;
			std::ofstream traj;
			traj.open(name_file);

			for(std::size_t k=0;k<melt.size();k++){
				for(std::size_t z=0;z<melt[k].chain.size();z++){
					traj<<melt[k].chain[z].t<<"\t"<<melt[k].chain[z].u<<"\t"<<melt[k].chain[z].v<<"\t";
				}
			}
			traj.close();
		}

		else{
			std::cout<<"There's been an error somewhere during the creation of the initial condition"<<std::endl;
			std::string name_file = "IC_"+length+"_"+std::to_string(polymer_length)+"_error";
			std::ofstream traj;
			traj.open(name_file);
			traj<<"There's been an error somewhere during the creation of the initial condition"<<std::endl;
		}
}

#endif

