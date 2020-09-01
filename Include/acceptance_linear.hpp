#ifndef __acceptance__
#define __acceptance__

#include "observables.hpp"


bool distance_nn_linear(int (&lattice)[Lattice_size][Lattice_size][Lattice_size], monomer& trial, polymer& p, unsigned int index_monomer){
		
		if(index_monomer>0 && index_monomer<p.chain.size()-1){
			
			double distance_plus = distance(trial,p.chain[index_monomer+1]);
			double distance_minus = distance(trial,p.chain[index_monomer-1]); 
			//std::cout<<lattice[trial.t_image][trial.u_image][trial.v_image]<<std::endl;
			if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus-1)<0.0001 && lattice[trial.t_image][trial.u_image][trial.v_image]==0)
				return true;

			else if(fabs(distance_plus)<0.00001 && fabs(distance_minus)<0.0001)
				return true;
			
			else if(fabs(distance_plus)<0.0001 && fabs(distance_minus-1)<0.0001)
				return true;
			
			else if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus)<0.0001)
				return true;
			
			else
				return false;
		}

		else if(index_monomer==0){
			//std::cout<<lattice[trial.t_image][trial.u_image][trial.v_image]<<std::endl;
			double distance_plus=distance(trial,p.chain[index_monomer+1]);
			if(fabs(distance_plus)<0.0001)
				return true;
			
			else if(fabs(distance_plus-1)<0.0001 && lattice[trial.t_image][trial.u_image][trial.v_image]==0)
				return true;
			
			else
				return false;
		}

		else{
			double distance_minus=distance(trial,p.chain[index_monomer-1]);
			if(fabs(distance_minus)<0.0001)
				return true;
			
			else if(fabs(distance_minus-1)<0.0001 && lattice[trial.t_image][trial.u_image][trial.v_image]==0)
				return true;
			
			else
				return false;
		}
}

bool distance_nn_linear_1_store(int (&lattice)[Lattice_size][Lattice_size][Lattice_size], monomer& trial, polymer& p, unsigned int index_monomer){
		
		if(index_monomer>0 && index_monomer<p.chain.size()-1){
			
			double distance_plus = distance(trial,p.chain[index_monomer+1]);
			double distance_minus = distance(trial,p.chain[index_monomer-1]); 

			if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus-1)<0.0001 && 
			lattice[trial.t_image][trial.u_image][trial.v_image]==0)
				return true;
			
			else if(fabs(distance_minus)<0.0001 && fabs(distance_plus-1)<0.0001 &&   
			lattice[trial.t_image][trial.u_image][trial.v_image]==1)
				return true;
			
			else if(fabs(distance_plus)<0.0001 && fabs(distance_minus-1)<0.0001 &&
			lattice[trial.t_image][trial.u_image][trial.v_image]==1)
				return true;

			else 
				return false;
		}

		else if(index_monomer==0){
			//std::cout<<lattice[trial.t_image][trial.u_image][trial.v_image]<<std::endl;
			double distance_plus = distance(trial,p.chain[1]);
			if(fabs(distance_plus)<0.0001 && lattice[trial.t_image][trial.u_image][trial.v_image]==1)
				return true;
			
			else if(fabs(distance_plus-1)<0.0001 && lattice[trial.t_image][trial.u_image][trial.v_image]==0)
				return true;
			
			else
				return false;
		}

		else{
			double distance_minus=distance(trial,p.chain[p.chain.size()-2]);
			
			if(fabs(distance_minus)<0.0001 && lattice[trial.t_image][trial.u_image][trial.v_image]==1)
				return true;
			
			else if(fabs(distance_minus-1)<0.0001 && lattice[trial.t_image][trial.u_image][trial.v_image]==0)
				return true;
			
			else
				return false;
		}
}

bool distance_nn_ideal_linear(monomer& trial, polymer& p, unsigned int index_monomer){
		
		if(index_monomer>0 && index_monomer<p.chain.size()-1){
			
			double distance_plus = distance(trial,p.chain[index_monomer+1]);
			double distance_minus = distance(trial,p.chain[index_monomer-1]); 
			//std::cout<<lattice[trial.t_image][trial.u_image][trial.v_image]<<std::endl;
			if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus-1)<0.0001)
				return true;

			else if(fabs(distance_plus)<0.0001 && fabs(distance_minus)<0.0001)
				return true;
			
			else if(fabs(distance_plus)<0.0001 && fabs(distance_minus-1)<0.0001)
				return true;
			
			else if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus)<0.0001)
				return true;
			
			else
				return false;
		}

		else if(index_monomer==0){
			//std::cout<<lattice[trial.t_image][trial.u_image][trial.v_image]<<std::endl;
			double distance_plus=distance(trial,p.chain[index_monomer+1]);
			if(fabs(distance_plus)<0.0001)
				return true;
			
			else if(fabs(distance_plus-1)<0.0001)
				return true;
			
			else
				return false;
		}

		else{
			double distance_minus=distance(trial,p.chain[index_monomer-1]);
			//std::cout<<lattice[trial.t_image][trial.u_image][trial.v_image]<<std::endl;
			if(fabs(distance_minus)<0.0001)
				return true;
			
			else if(fabs(distance_minus-1)<0.0001)
				return true;
			
			else
				return false;
		}
}

bool distance_nn_ideal_linear_1_store(monomer& trial, polymer& p, unsigned int index_monomer){
		
		if(index_monomer>0 && index_monomer < p.chain.size()-1){
			
			double distance_plus = distance(trial,p.chain[index_monomer+1]);
			double distance_minus = distance(trial,p.chain[index_monomer-1]); 
			//std::cout<<lattice[trial.t_image][trial.u_image][trial.v_image]<<std::endl;
			if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus-1)<0.0001)
				return true;

			else if(fabs(distance_plus)<0.0001 && fabs(distance_minus)<0.0001)
				return false;
			
			else if(fabs(distance_minus)<0.0001 && fabs(distance_plus-1)<0.0001){
				if(index_monomer==1)
					return true;
				
				else{
					if(fabs(distance(p.chain[index_monomer-1],p.chain[index_monomer-2])-1)<0.0001)
						return true;
					else
						return false;
				}
			}
			
			else if(fabs(distance_plus)<0.0001 && fabs(distance_minus-1)<0.0001){
				if(index_monomer==p.chain.size()-2)
					return true;
				
				else{
					if(fabs(distance(p.chain[index_monomer+1],p.chain[index_monomer+2])-1)<0.0001)
						return true;
					else
						return false;
				}
			}
			
			else
				return false;
		}

		else if(index_monomer==0){
			//std::cout<<lattice[trial.t_image][trial.u_image][trial.v_image]<<std::endl;
			double distance_plus = distance(trial,p.chain[index_monomer+1]);
			if(fabs(distance_plus)<0.0001 && fabs(distance(p.chain[1],p.chain[2])-1)<0.0001)
				return true;
			
			else if(fabs(distance_plus-1)<0.0001)
				return true;
			
			else
				return false;
		}

		else{
			double distance_minus=distance(trial,p.chain[p.chain.size()-2]);
			
			if(fabs(distance_minus)<0.0001 && fabs(distance(p.chain[p.chain.size()-2],p.chain[p.chain.size()-3])-1)<0.0001)
				return true;
			
			else if(fabs(distance_minus-1)<0.0001)
				return true;
			
			else
				return false;
		}
}

bool acceptance_ideal_linear(monomer trial, polymer& a, unsigned int index_monomer){
		return distance_nn_ideal_linear_1_store(trial,a,index_monomer);
}

bool acceptance_ideal_ring(monomer trial, polymer& a, unsigned int index_monomer){
		return distance_nn_ideal_rings_1_store(trial,a,index_monomer);
}

bool acceptance_real_linear(long long int MC_step,
	int (&lattice)[Lattice_size][Lattice_size][Lattice_size],monomer trial,polymer& a,const unsigned int index_monomer){
		/*if(MC_step>2e10)
			return distance_nn_linear_1_store(lattice,trial,a,index_monomer);
		else
			return distance_nn_linear(lattice,trial,a,index_monomer);*/
		return distance_nn_linear_1_store(lattice,trial,a,index_monomer);
}

bool acceptance_real_ring(long long int MC_step,int (&lattice)[Lattice_size][Lattice_size][Lattice_size],monomer trial,polymer& a,const unsigned int index_monomer){
		/*if(MC_step>1e6)
			return distance_nn_rings_1_store(lattice,trial,a,index_monomer);
		else
			return distance_nn_rings(lattice,trial,a,index_monomer);*/
		return distance_nn_rings_1_store(lattice,trial,a,index_monomer);
}

#endif
