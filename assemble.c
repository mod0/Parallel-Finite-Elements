#include "domain.h"
# include <math.h>
void* assemble(domain*);
void* assemble_local_K(domain*,double[3][3],double[3]);
void* set_boundary_values();



// takes domain object creates K and M matrices for it. 
void* assemble(domain* Mydomain) {
	//Todo: Initialize K and M matrices in sparse form
	int subdomain_idx = 0;
	double Ktilde[3][3] = {{1,-0.5, -0.5},{-0.5,0.5,0},{-0.5,0,0.5}};
	double h = Mydomain->cartesian_grid->h;
	double offdiag = h*h*1.0/24 ; double diag= h*h*1.0/12;
	double Mtilde[3][3] = {{diag,offdiag, offdiag},{offdiag,diag,offdiag}
	                       ,{offdiag,offdiag,diag}};
	double vol = 1.0/6 * h*h; 
	double Ftilde[3] = {vol, vol,vol};
	// Assemble K 
	assemble_local_K(Mydomain,Ktilde,Ftilde);
    
} 


void* assemble_local_K(domain* Mydomain, double[3][3] Ktilde, double[3] Ftilde){
	int subdomain_idx=0; // comes from omp
	int global_vertice_idx, global_element_idx,Nv, Nx, Ny, Nt;
	int local_i, local_j;
	Nx = Mydomain.subdomain[subdomain_idx].dimX;
	Ny = Mydomain.subdomain[subdomain_idx].dimY;
	Nv = Nx* Ny;
	Nt = (Nx-1)* (Ny -1) *2;
	double F[Nv];
	int* global_vertices;
	vector bands[4];
	int vector_sizes[4] = {Nv,Nv-1,Nv-Nx, Nv-Nx-1};
	for (int i=0; i<4; i++ ){
		bands[i] = vector_init(vector_sizes[i]);
	}

	
	for(int t_idx=1; t_idx<=Nt; t_idx++ ){
		global_element_idx = Mydomain.subdomain[subdomain_idx].elements[t_idx];
		global_vertices = _triangular_elements[global_element_idx].grid_vertex;
		for(int k=1; k<=3; k++){
			local_i =  Mydomain.subdomain_vertex_map[subdomain_idx][global_vertices[k-1]];
			for(int l=1; l<=3; l++){
				local_j =  Mydomain.subdomain_vertex_map[subdomain_idx][vertex_pointer[l-1]];
				//matrix_increment_value(subdomain_K, locali,localj,Ktilde[k-1][l-1]);
				
				if (local_j> local_i){
					switch (local_j -local_i){
						case 0:
							bands[0].elements[local_i] += Ktilde[k-1][l-1];
							break;
						case 1 :
							band[1].elements[local_i] += Ktilde[k-1][l-1];
							break;
						case Nx:
							band[2].elements[local_i] += Ktilde[k-1][l-1];
							break;
						case Nx+1:
							band[3].elements[local_i] += Ktilde[k-1][l-1];
							break;
					}
				}
			}
			F[local_i-1] = F[local_i-1] + Ftilde[k-1];
		}
	}
	// Create sparse K 
	// K=create_symmetric_sparse(bands); 
}

void* set_boundary_values(){
	int subdomain_idx, p; 
	if (subdomain_idx == 0) {


	}
	else if (subdomain_idx==p-1){

	}
	else{

	}

}