#ifndef __monte_carlo__
#define __monte_carlo__

#include <string>
#include "acceptance_ring.hpp"

struct MC_move{		
	int dt;
	int du;
	int dv;
	MC_move(int t,int u,int v):dt{t},du{u},dv{v} {};
};

void update_configuration(int (&lattice)[Lattice_size][Lattice_size][Lattice_size],
		int (&lattice1)[Lattice_size][Lattice_size][Lattice_size],int (&lattice2)[Lattice_size][Lattice_size][Lattice_size],
		monomer& trial,unsigned int monomer_to_move, unsigned int polymer_to_move,std::vector<polymer>& melt){

		lattice[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]-=1;

		if(lattice1[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]==polymer_to_move){
			
			double distance_plus;
			double distance_minus;
		
			if(monomer_to_move>0 && monomer_to_move<melt[polymer_to_move].chain.size()-1){
				distance_plus = distance(trial,melt[polymer_to_move].chain[monomer_to_move+1]);
				distance_minus = distance(trial,melt[polymer_to_move].chain[monomer_to_move-1]);
			}

			else if(monomer_to_move==0){
				distance_minus=distance(trial,melt[polymer_to_move].chain[melt[polymer_to_move].chain.size()-1]);
				distance_plus=distance(trial,melt[polymer_to_move].chain[1]);
			}

			else{
				distance_minus = distance(trial,melt[polymer_to_move].chain[melt[polymer_to_move].chain.size()-2]);
				distance_plus = distance(trial,melt[polymer_to_move].chain[0]);
			}

			if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus-1)<0.0001)
				lattice1[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]=-1;
		}

		if(lattice2[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]==polymer_to_move){
			
			double distance_plus;
			double distance_minus;
		
			if(monomer_to_move>0 && monomer_to_move<melt[polymer_to_move].chain.size()-1){
				distance_plus = distance(trial,melt[polymer_to_move].chain[monomer_to_move+1]);
				distance_minus = distance(trial,melt[polymer_to_move].chain[monomer_to_move-1]);
			}

			else if(monomer_to_move==0){
				distance_minus=distance(trial,melt[polymer_to_move].chain[melt[polymer_to_move].chain.size()-1]);
				distance_plus=distance(trial,melt[polymer_to_move].chain[1]);
			}

			else{
				distance_minus = distance(trial,melt[polymer_to_move].chain[melt[polymer_to_move].chain.size()-2]);
				distance_plus = distance(trial,melt[polymer_to_move].chain[0]);
			}
			
			if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus-1)<0.0001)
				lattice2[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]=-1;
		}


		melt[polymer_to_move].chain[monomer_to_move].t=trial.t;
		melt[polymer_to_move].chain[monomer_to_move].u=trial.u;
		melt[polymer_to_move].chain[monomer_to_move].v=trial.v;
		melt[polymer_to_move].chain[monomer_to_move].t_image=trial.t_image;
		melt[polymer_to_move].chain[monomer_to_move].u_image=trial.u_image;
		melt[polymer_to_move].chain[monomer_to_move].v_image=trial.v_image;
		
		lattice[melt[polymer_to_move].chain[monomer_to_move].t_image][melt[polymer_to_move].chain[monomer_to_move].u_image][melt[polymer_to_move].chain[monomer_to_move].v_image]+=1;
}

void MC_routine_rings_links(std::vector<int> coordinates,int polymer_length, 
	int number_of_polymers, double monomer_density, long long int MC_steps, int stride,int run){
	std::string DIR = "/net/sbp/sbpstore1/mubertin/DATA/RINGS/";
	
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
	
	for(auto i=run;i<run+1;i++){

		srand(time(NULL));
		

		std::vector<polymer> melt;
		
		int lattice[Lattice_size][Lattice_size][Lattice_size]={0};
		int lattice1[Lattice_size][Lattice_size][Lattice_size];
		int lattice2[Lattice_size][Lattice_size][Lattice_size];

		for(auto i=0;i<Lattice_size;i++){
			for(auto j=0;j<Lattice_size;j++){
				for(auto k=0;k<Lattice_size;k++){
						lattice1[i][j][k]=-1;
						lattice2[i][j][k]=-1;
					}
				}
			}

		for(int n_pol=0;n_pol<number_of_polymers;n_pol++){
			polymer polymer_tmp{polymer_length,n_pol};
			for(int j=0;j<polymer_length;j++){
				
				int t = coordinates[3*polymer_length*n_pol+3*j];
				int u = coordinates[3*polymer_length*n_pol+3*j+1];
				int v = coordinates[3*polymer_length*n_pol+3*j+2];
				
				monomer monomer_tmp{t,u,v};
				polymer_tmp.chain.push_back(monomer_tmp);
				
				lattice[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]+=1;
				
				if(lattice1[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]==-1){
					lattice1[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]=n_pol;
				}

				else if(lattice1[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]==n_pol);

				else if(lattice1[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]!=-1 
						&& lattice1[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image] != n_pol 
						&& lattice2[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image] == -1){		
						lattice2[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]=n_pol;
				}

				else if(lattice2[monomer_tmp.t_image][monomer_tmp.u_image][monomer_tmp.v_image]==n_pol); 
				
				else{
					MC_steps=0;
					std::string name_file_error= "Error_"+length+"_"+std::to_string(polymer_length)+"_"+std::to_string(i);
					std::ofstream Error{name_file_error};
					Error<<"There are some problems with initial conditions: more than two polymers in the same lattice size."<<std::endl;
					std::cout<<"\t Error"<<std::endl;
				}			
			}
			melt.push_back(polymer_tmp);
		}

		double temp=0;
		for(auto i=0;i<Lattice_size;i++){
			for(auto j=0;j<Lattice_size;j++){
				for(auto k=0;k<Lattice_size;k++){
					if(lattice2[i][j][k]!=-1){
						temp+=1;
						std::cout<<lattice[i][j][k]<<std::endl;
						std::cout<<i<<"\t"<<j<<"\t"<<k<<std::endl;
					}
				}
			}
		}
		
		std::string name_file1= DIR+"GR_real_"+length+"_"+std::to_string(polymer_length)+"_"+std::to_string(i);
		std::ofstream GR{name_file1};
		
		int n_file=0;
		int acc = 0;
		
		std::string name_file3 = DIR+"traj_real_"+length+"_"+std::to_string(polymer_length)+"_"+std::to_string(n_file)+std::to_string(i);
		std::string name_file_R_CM= DIR+"R_CM_real_"+length+"_"+std::to_string(polymer_length)+"_"+std::to_string(n_file)+std::to_string(i);
		std::ofstream traj;
		traj.open(name_file3);
		std::ofstream R_CM;
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
				
				if(acceptance_real_ring_links(j,lattice,lattice1,lattice2,trial,melt[polymer_to_move],monomer_to_move,polymer_to_move))
				{
					update_configuration(lattice,lattice1,lattice2,trial,monomer_to_move,polymer_to_move,melt);
					acc++;
				}

				if(j%stride==0 && j>=0){
					double GR_average = 0;
					for(std::size_t k=0;k<melt.size();k++){
						auto r_cm=compute_R_cm(melt[k]);
						R_CM<<r_cm[0]<<"\t"<<r_cm[1]<<"\t"<<r_cm[2]<<"\t";
						GR_average += gyration_radius(melt[k]);
						for(std::size_t z=0;z<melt[k].chain.size();z++){
							traj<<melt[k].chain[z].t<<"\t"<<melt[k].chain[z].u<<"\t"<<melt[k].chain[z].v<<"\t";
						}
					}
					GR_average = GR_average/melt.size();
					GR<<GR_average<<std::endl;
					traj<<std::endl;
					R_CM<<std::endl;
				}
				

		}

		double temp_1=0;
		for(auto i=0;i<Lattice_size;i++){
			for(auto j=0;j<Lattice_size;j++){
				for(auto k=0;k<Lattice_size;k++){
					if(lattice[i][j][k]>2 || lattice[i][j][k]<0){
						temp_1+=1;
						std::cout<<lattice[i][j][k]<<std::endl;
						std::cout<<i<<"\t"<<j<<"\t"<<k<<std::endl;
					}
				}
			}
		}

		std::cout<<temp_1/(Lattice_size*Lattice_size*Lattice_size*monomer_density)<<std::endl;
		std::cout<<double(acc)/double(MC_steps)<<std::endl;
			
	}
}

#endif
