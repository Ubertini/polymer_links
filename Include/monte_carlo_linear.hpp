#ifndef __monte_carlo__
#define __monte_carlo__

#include <string>
#include "acceptance_linear.hpp"

struct MC_move{		
	int dt;
	int du;
	int dv;
	MC_move(int t,int u,int v):dt{t},du{u},dv{v} {};
};

void MC_routine_linear(std::vector<int> coordinates, int polymer_length, int number_of_polymers,
	double monomer_density, long long int MC_steps, int stride){
	
	std::string DIR = "/net/sbp/sbpstore1/mubertin/DATA/FJC/";

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

	std::string length= std::to_string(int(monomer_density*100));
	
	for(auto i=10;i<11;i++){

		srand(time(NULL));
		std::string name_file1= DIR+"GR_real_"+length+"_"+std::to_string(polymer_length)+"_"+std::to_string(i);
		std::string name_file2= DIR+"R_EE_real_"+length+"_"+std::to_string(polymer_length)+"_"+std::to_string(i);
		
		std::ofstream R_EE{name_file2};
		std::ofstream GR{name_file1};

		std::vector<polymer> melt;
		
		int lattice[Lattice_size][Lattice_size][Lattice_size]={0};

	
		
		for(int n_pol=0;n_pol<number_of_polymers;n_pol++){
			polymer polymer_tmp{polymer_length,n_pol};
			for(int j=0;j<polymer_length;j++){
				int t = coordinates[3*polymer_length*n_pol+3*j];
				int u = coordinates[3*polymer_length*n_pol+3*j+1];
				int v = coordinates[3*polymer_length*n_pol+3*j+2];
				monomer monomer_tmp{t,u,v};
				polymer_tmp.chain.push_back(monomer_tmp);
				lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]+=1;
			}
			melt.push_back(polymer_tmp);
		}
	/*	for(int i=0;i<number_of_polymers;i++){
			int t= rand() % Lattice_size;
			int u= rand() % Lattice_size;
			int v= rand() % Lattice_size;
			
			while(lattice[t][u][v]!=0){
				t= rand() % Lattice_size;
				u= rand() % Lattice_size;
				v= rand() % Lattice_size;
			}

			lattice[t][u][v]+=polymer_length;
			polymer temp{polymer_length,i,t,u,v};
			melt.push_back(temp);
	}*/

		int n_file=0;
		int acc = 0;
		std::string name_file3 = DIR+"traj_real_"+length+"_"+std::to_string(polymer_length)+"_"+std::to_string(n_file)+"_"+std::to_string(i);
		std::string name_file_R_CM= DIR+"R_CM_real_"+length+"_"+std::to_string(polymer_length)+"_"+std::to_string(n_file)+"_"+std::to_string(i);

		std::ofstream traj;
		std::ofstream R_CM;
		traj.open(name_file3);
		R_CM.open(name_file_R_CM);

		for(long long int j=0;j<MC_steps;j++){
				


				if(j%stride==0 && j>0){
					n_file+=1;
					traj.close();
					R_CM.close();
					name_file3 = DIR+"traj_real_"+length+"_"+std::to_string(polymer_length)+"_"+std::to_string(n_file)+"_"+std::to_string(i);
					std::string name_file_R_CM= DIR+"R_CM_real_"+length+"_"+std::to_string(polymer_length)+"_"+std::to_string(n_file)+"_"+std::to_string(i);
					R_CM.open(name_file_R_CM);
					traj.open(name_file3);
				}
				
				unsigned int polymer_to_move = rand()%melt.size();
				unsigned int monomer_to_move = rand()%melt[polymer_to_move].chain.size();
				auto MC_move=rand()%12; //there are 12 possible movement
				monomer trial{melt[polymer_to_move].chain[monomer_to_move].t+mc_move[MC_move].dt,melt[polymer_to_move].chain[monomer_to_move].u+mc_move[MC_move].du,melt[polymer_to_move].chain[monomer_to_move].v+mc_move[MC_move].dv};
				
				if(acceptance_real_linear(j,lattice,trial,melt[polymer_to_move],monomer_to_move)){
					
					lattice[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]-=1;

					melt[polymer_to_move].chain[monomer_to_move].t=trial.t;
					melt[polymer_to_move].chain[monomer_to_move].u=trial.u;
					melt[polymer_to_move].chain[monomer_to_move].v=trial.v;
					melt[polymer_to_move].chain[monomer_to_move].t_image=trial.t_image;
					melt[polymer_to_move].chain[monomer_to_move].u_image=trial.u_image;
					melt[polymer_to_move].chain[monomer_to_move].v_image=trial.v_image;
					lattice[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]+=1;
					acc++;
				}

				if(j%stride==0 && j>=0){
					double GR_average = 0;
					double R_EE_average = 0;
					for(std::size_t k=0;k<melt.size();k++){
						auto r_cm=compute_R_cm(melt[k]);
						R_CM<<r_cm[0]<<"\t"<<r_cm[1]<<"\t"<<r_cm[2]<<"\t";
						GR_average += gyration_radius(melt[k]);
						R_EE_average += end_end_distance(melt[k]);
						for(std::size_t z=0;z<melt[k].chain.size();z++){
							traj<<melt[k].chain[z].t<<"\t"<<melt[k].chain[z].u<<"\t"<<melt[k].chain[z].v<<"\t";
						}
					}
					GR_average = GR_average/melt.size();
					R_EE_average = R_EE_average/melt.size();
					GR<<GR_average<<std::endl;
					R_EE<<R_EE_average<<std::endl;
					traj<<std::endl;
					R_CM<<std::endl;
				}
		}

		double temp=0;
		for(auto i=0;i<Lattice_size;i++){
			for(auto j=0;j<Lattice_size;j++){
				for(auto k=0;k<Lattice_size;k++){
					if(lattice[i][j][k]>2){
						temp+=1;
						
					}
				}
			}
		}

		std::cout<<temp/(Lattice_size*Lattice_size*Lattice_size*monomer_density)<<std::endl;
		std::cout<<double(acc)/double(MC_steps)<<std::endl;
			
	}
}















#endif
