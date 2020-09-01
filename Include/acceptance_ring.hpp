#ifndef __acceptance__
#define __acceptance__

#include "observables.hpp"

bool distance_nn_rings(int (&lattice)[Lattice_size][Lattice_size][Lattice_size], monomer& trial, polymer& p, unsigned int index_monomer){
		
		double distance_plus;
		double distance_minus;

		if(index_monomer > 0 && index_monomer<p.chain.size()-1){
			
			distance_plus = distance(trial,p.chain[index_monomer+1]);
			distance_minus = distance(trial,p.chain[index_monomer-1]);
		}

		else if(index_monomer==0){
			distance_minus=distance(trial,p.chain[p.chain.size()-1]);
			distance_plus=distance(trial,p.chain[index_monomer+1]);
		}

		else{
			distance_minus = distance(trial,p.chain[index_monomer-1]);
			distance_plus = distance(trial,p.chain[0]);
		}
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

bool distance_nn_rings_1_store(int (&lattice)[Lattice_size][Lattice_size][Lattice_size], monomer& trial, polymer& p, unsigned int index_monomer){
		
		double distance_plus;
		double distance_minus;


		if(index_monomer>0 && index_monomer<p.chain.size()-1){
			
			distance_plus = distance(trial,p.chain[index_monomer+1]);
			distance_minus = distance(trial,p.chain[index_monomer-1]);
		}

		else if(index_monomer==0){
			distance_minus=distance(trial,p.chain[p.chain.size()-1]);
			distance_plus=distance(trial,p.chain[1]);
		}

		else{
			distance_minus = distance(trial,p.chain[p.chain.size()-2]);
			distance_plus = distance(trial,p.chain[0]);
		}
		
		if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus-1)<0.0001 && 
			lattice[trial.t_image][trial.u_image][trial.v_image]==0)
			return true;
			
		else if(fabs(distance_plus)<0.0001 && fabs(distance_minus-1)<0.0001 && 
			lattice[trial.t_image][trial.u_image][trial.v_image]==1)
			return true;
			
		else if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus)<0.0001 && 
			lattice[trial.t_image][trial.u_image][trial.v_image]==1)
			return true;
			
		else
			return false;
}

bool distance_nn_rings_1_store_links(int (&lattice)[Lattice_size][Lattice_size][Lattice_size],
		int (&lattice1)[Lattice_size][Lattice_size][Lattice_size],int (&lattice2)[Lattice_size][Lattice_size][Lattice_size],
		monomer& trial, polymer& p, const unsigned int index_monomer, const int n_pol){
		
		double distance_plus;
		double distance_minus;


		if(index_monomer>0 && index_monomer<p.chain.size()-1){
			
			distance_plus = distance(trial,p.chain[index_monomer+1]);
			distance_minus = distance(trial,p.chain[index_monomer-1]);
		}

		else if(index_monomer==0){
			distance_minus=distance(trial,p.chain[p.chain.size()-1]);
			distance_plus=distance(trial,p.chain[1]);
		}

		else{
			distance_minus = distance(trial,p.chain[p.chain.size()-2]);
			distance_plus = distance(trial,p.chain[0]);
		}
		
		if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus-1)<0.0001 && 
			lattice[trial.t_image][trial.u_image][trial.v_image]<=2){
				
				if(lattice1[trial.t_image][trial.u_image][trial.v_image]==-1){
					lattice1[trial.t_image][trial.u_image][trial.v_image]=n_pol;
					return true;
				}

				else if(lattice2[trial.t_image][trial.u_image][trial.v_image]==-1){
					lattice2[trial.t_image][trial.u_image][trial.v_image]=n_pol;
					return true;
				}

				else
					return false;
		}

		else if(fabs(distance_plus)<0.0001 && fabs(distance_minus-1)<0.0001 && 
			lattice[trial.t_image][trial.u_image][trial.v_image]<=3){
							
				double next_neighbor_plus;
				
				if(index_monomer<p.chain.size()-2){
					next_neighbor_plus=distance(p.chain[index_monomer+1],p.chain[index_monomer+2]);
				}

				else if(index_monomer==p.chain.size()-2){
					next_neighbor_plus=distance(p.chain[p.chain.size()-1],p.chain[0]);
				}

				else{
					next_neighbor_plus=distance(p.chain[0],p.chain[1]);
				}
				
				if(fabs(next_neighbor_plus-1)<0.0001)
					return true;

				else
					return false;
		}

		else if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus)<0.0001 && 
			lattice[trial.t_image][trial.u_image][trial.v_image]<=3){
			
			double next_neighbor_minus;
			
			if(index_monomer>1){
				next_neighbor_minus=distance(p.chain[index_monomer-1],p.chain[index_monomer-2]);
			}

			else if(index_monomer==1){
				next_neighbor_minus=distance(p.chain[0],p.chain[p.chain.size()-1]);
			}

			else{
				next_neighbor_minus=distance(p.chain[p.chain.size()-1],p.chain[p.chain.size()-2]);
			}

			if(fabs(next_neighbor_minus-1)<0.0001)
				return true;
				
			else
				return false;
		
		}
			
		else
			return false;
}

bool distance_nn_ideal_rings(monomer& trial, polymer& p, unsigned int index_monomer){
		
		double distance_plus;
		double distance_minus;

		if(index_monomer>0 && index_monomer<p.chain.size()-1){
			
			distance_plus = distance(trial,p.chain[index_monomer+1]);
			distance_minus = distance(trial,p.chain[index_monomer-1]);
		}

		else if(index_monomer==0){
			distance_minus=distance(trial,p.chain[p.chain.size()-1]);
			distance_plus=distance(trial,p.chain[index_monomer+1]);
		}

		else{
			distance_minus = distance(trial,p.chain[index_monomer-1]);
			distance_plus = distance(trial,p.chain[0]);
		}

		if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus-1)<0.0001)
			return true;

		/*else if(fabs(distance_plus)<0.00001 && fabs(distance_minus)<0.0001)
			return true;*/
			
		/*else if(fabs(distance_plus)<0.0001 && fabs(distance_minus-1)<0.0001)
			return true;
			
		else if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus)<0.0001)
			return true;*/
			
		else
			return false;
}

bool distance_nn_ideal_rings_1_store(monomer& trial, polymer& p, unsigned int index_monomer){
		
		double distance_plus;
		double distance_minus;
		double next_neighbor_plus;
		double next_neighbor_minus;

		if(index_monomer>1 && index_monomer<p.chain.size()-2){
			
			distance_plus = distance(trial,p.chain[index_monomer+1]);
			distance_minus = distance(trial,p.chain[index_monomer-1]);
			next_neighbor_plus=distance(p.chain[index_monomer+1],p.chain[index_monomer+2]);
			next_neighbor_minus=distance(p.chain[index_monomer-1],p.chain[index_monomer-2]);
		}

		else if(index_monomer==0){
			distance_minus=distance(trial,p.chain[p.chain.size()-1]);
			distance_plus=distance(trial,p.chain[1]);
			next_neighbor_plus=distance(p.chain[index_monomer+1],p.chain[index_monomer+2]);
			next_neighbor_minus=distance(p.chain[p.chain.size()-1],p.chain[p.chain.size()-2]);
		}

		else if(index_monomer==1){
			distance_minus=distance(trial,p.chain[0]);
			distance_plus=distance(trial,p.chain[2]);
			next_neighbor_plus=distance(p.chain[2],p.chain[3]);
			next_neighbor_minus=distance(p.chain[0],p.chain[p.chain.size()-1]);
		}

		else if(index_monomer==p.chain.size()-2){
			distance_minus=distance(trial,p.chain[p.chain.size()-3]);
			distance_plus=distance(trial,p.chain[p.chain.size()-1]);
			next_neighbor_plus=distance(p.chain[p.chain.size()-1],p.chain[0]);
			next_neighbor_minus=distance(p.chain[p.chain.size()-3],p.chain[p.chain.size()-4]);
		}

		else{
			distance_minus = distance(trial,p.chain[p.chain.size()-2]);
			distance_plus = distance(trial,p.chain[0]);
			next_neighbor_plus=distance(p.chain[0],p.chain[1]);
			next_neighbor_minus=distance(p.chain[p.chain.size()-2],p.chain[p.chain.size()-3]);
		}
			//std::cout<<lattice[trial.t_image][trial.u_image][trial.v_image]<<std::endl;
		if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus-1)<0.0001)
			return true;

		else if(fabs(distance_plus)<0.00001 && fabs(distance_minus)<0.0001)
			return false;
			
		else if(fabs(distance_plus)<0.0001 && fabs(distance_minus-1)<0.0001 && fabs(next_neighbor_plus-1)<0.0001)
			return true;
			
		else if(fabs(distance_plus-1)<0.0001 && fabs(distance_minus)<0.0001 && fabs(next_neighbor_minus-1)<0.0001)
			return true;
			
		else
			return false;
}

bool acceptance_ideal_ring(monomer trial, polymer& a, unsigned int index_monomer){
		return distance_nn_ideal_rings_1_store(trial,a,index_monomer);
}

bool acceptance_real_ring(long long int MC_step,int (&lattice)[Lattice_size][Lattice_size][Lattice_size],monomer trial,polymer& a,const unsigned int index_monomer){
		/*if(MC_step>1e6)
			return distance_nn_rings_1_store(lattice,trial,a,index_monomer);
		else
			return distance_nn_rings(lattice,trial,a,index_monomer);*/
		return distance_nn_rings_1_store(lattice,trial,a,index_monomer);
}

bool acceptance_real_ring_links(long long int MC_step,int (&lattice)[Lattice_size][Lattice_size][Lattice_size],
															int (&lattice1)[Lattice_size][Lattice_size][Lattice_size],int (&lattice2)[Lattice_size][Lattice_size][Lattice_size],
															monomer& trial, polymer& p, const unsigned int index_monomer, const int n_pol){
		
		return distance_nn_rings_1_store_links(lattice,lattice1,lattice2,trial,p,index_monomer,n_pol);

}

#endif
