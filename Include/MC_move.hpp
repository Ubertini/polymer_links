#ifndef __acceptance__
#define __acceptance__

bool acceptance_ideal_linear(monomer trial, polymer& a, unsigned int index_monomer){
		bool movement=false;	
		if(index_monomer<(a.v.size()-1) && index_monomer>0){
			if(distance_nn(trial, a.v[index_monomer-1]) && distance_nn(trial, a.v[index_monomer+1]))
				movement = true;
		}
		else if(index_monomer==a.v.size()-1){
			if(distance_nn(trial, a.v[index_monomer-1])) // && distance_nn(trial, a.v[0])) //second check for rings
				movement = true;
		}
		else{
			if(distance_nn(trial, a.v[index_monomer+1])) //&& distance_nn(trial, a.v[a.v.size()-1])) //second check for rings
				movement = true;
		}
		return movement;
}

bool acceptance_ideal_ring(monomer trial, polymer& a, unsigned int index_monomer){
		bool movement=false;	
		if(index_monomer<(a.v.size()-1) && index_monomer>0){
			if(distance_nn(trial, a.v[index_monomer-1]) && distance_nn(trial, a.v[index_monomer+1]))
				movement = true;
		}
		else if(index_monomer==a.v.size()-1){
			if(distance_nn(trial, a.v[index_monomer-1]) && distance_nn(trial, a.v[0])) //second check for rings
				movement = true;
		}
		else{
			if(distance_nn(trial, a.v[index_monomer+1]) && distance_nn(trial, a.v[a.v.size()-1])) //second check for rings
				movement = true;
		}
		return movement;
}


bool lattice_free(monomer trial, polymer& a, const unsigned int index_monomer){
		
		if(index_monomer != 0 && index_monomer != (a.v.size()-1)){
			for(unsigned int i=0;i<index_monomer-1;i++){
				//if(!(trial.t_image != a.v[i].t_image || trial.u_image != a.v[i].u_image || trial.v_image != a.v[i].v_image))
					if(fabs(distance_images(trial,a.v[i]))<0.0001)
						return false;
			}
			
			for(unsigned int i = index_monomer+2; i<a.v.size();i++){
				//if(!(trial.t_image != a.v[i].t_image || trial.u_image != a.v[i].u_image || trial.v_image != a.v[i].v_image))
				if(fabs(distance_images(trial,a.v[i]))<0.0001)
					return false;
			}
		}

		else if(index_monomer == 0){
			//to change the initial point at i=a.v.size()-1 if you want ring
			for(unsigned int i=2; i<a.v.size(); i++){
				//if(!(trial.t_image != a.v[i].t_image || trial.u_image != a.v[i].u_image || trial.v_image != a.v[i].v_image))
				if(fabs(distance_images(trial,a.v[i]))<0.0001)	
					return false;
			}																												
		}

		else{
			//to change the initial point at i=1 if you want ring
			for(unsigned int i=0;i<a.v.size()-2;i++){
				//if(!(trial.t_image != a.v[i].t_image && trial.u_image != a.v[i].u_image && trial.v_image != a.v[i].v_image))
				if(fabs(distance_images(trial,a.v[i]))<0.0001)
					return false;
			}
		}

		return true;
}

bool acceptance_real_linear(monomer trial,polymer& a,const unsigned int index_monomer){
		
		bool movement=false;	
		
		if(index_monomer<(a.v.size()-1) && index_monomer>0){
			if(distance_nn(trial, a.v[index_monomer-1]) && distance_nn(trial, a.v[index_monomer+1])){
				if(lattice_free(trial,a,index_monomer))
					movement = true;
			}
		}
		
		else if(index_monomer==a.v.size()-1){
			if(distance_nn(trial, a.v[index_monomer-1])){
				if(lattice_free(trial,a,index_monomer))
					movement = true;
			}				
		}
		
		else{
			if(distance_nn(trial, a.v[index_monomer+1])){
				if(lattice_free(trial,a,index_monomer))
				movement = true;
			}
		}
		return movement;
}

bool acceptance_real_ring(monomer trial,polymer& a,const unsigned int index_monomer){
		
		bool movement=false;	
		
		if(index_monomer<(a.v.size()-1) && index_monomer>0){
			if(distance_nn(trial, a.v[index_monomer-1]) && distance_nn(trial, a.v[index_monomer+1])){
				if(lattice_free(trial,a,index_monomer))
					movement = true;
			}
		}
		
		else if(index_monomer==a.v.size()-1){
			if(distance_nn(trial, a.v[index_monomer-1]) && distance_nn(trial, a.v[0])){ //second check for rings
				if(lattice_free(trial,a,index_monomer))
					movement = true;
			}				
		}
		
		else{
			if(distance_nn(trial, a.v[index_monomer+1])) && distance_nn(trial, a.v[a.v.size()-1]){ //second check for rings
				if(lattice_free(trial,a,index_monomer))
				movement = true;
			}
		}
		return movement;
}





#endif