#include "domain.h"
#include <stdio.h>
#include <math.h>
#include "elements.h"
#include "error.h"
void* assemble(domain*);
void* assemble_local_K(domain* , double* , double*, int );
void* set_boundary_values(domain* D, int subdomain_idx, vector * bands, double* F );
double bc_choice(double* hat_u,double u0, int subdomain_idx, int local_coor);



// takes domain object creates K and M matrices for it. 
void* assemble(domain* D) {
	//Todo: Initialize K and M matrices in sparse form
	int subdomain_idx; int p =2; // coming from omp
	double Ktilde[3][3] = {{1,-0.5, -0.5},{-0.5,0.5,0},{-0.5,0,0.5}};
	double h = D->cartesian_grid->h;
	double offdiag = h*h*1.0/24 ; double diag= h*h*1.0/12;
	double Mtilde[3][3] = {{diag,offdiag, offdiag},{offdiag,diag,offdiag}
	                       ,{offdiag,offdiag,diag}};
	double vol = 1.0/6 * h*h; //Todo
	double Ftilde[3] = {vol, vol,vol};
	for (int subdomain_idx=0; subdomain_idx< p; subdomain_idx++){
		// Assemble K 
		assemble_local_K(D, Ktilde[0], Ftilde, subdomain_idx);
        set_boundary_values();
	}
    
} 

void* assemble_local_K(domain* D, double* Ktilde, double* Ftilde, int subdomain_idx){

	int global_vertice_idx, global_element_idx,Nv, Nx, Ny, Nt;
	int local_i, local_j;
	Nx = D->subdomains[subdomain_idx].dimX;
	Ny = D->subdomains[subdomain_idx].dimY;
	Nv = Nx* Ny;
	Nt = (Nx-1)* (Ny -1) *2;
	double F[Nv];
	vertex** triple_vertices;
	vector bands[4];
	int vector_sizes[4] = {Nv,Nv-1,Nv-Nx, Nv-Nx-1};
	for (int i=0; i<4; i++ ){
		bands[i] = vector_init(vector_sizes[i]);
	}

	
	for(int t_idx=1; t_idx<=Nt; t_idx++ ){
		global_element_idx = D->subdomains[subdomain_idx].elements[t_idx];
		triple_vertices = _triangular_elements[global_element_idx].grid_vertex;

		for(int k=1; k<=3; k++){
			local_i =  D->subdomain_vertex_map[subdomain_idx][triple_vertices[k-1]->id];
			for(int l=1; l<=3; l++){
				local_j =  D->subdomain_vertex_map[subdomain_idx][triple_vertices[l-1]->id];
				//matrix_increment_value(subdomain_K, locali,localj,Ktilde[k-1][l-1]);
				
				if (local_j>= local_i){
					int Delta = local_j -local_i;
					if(Delta == 0)
							bands[0].elements[local_i] += Ktilde[k-1][l-1];
				    else if(Delta==1)
							band[1].elements[local_i] += Ktilde[k-1][l-1];
					else if(Delta == Nx)
						  	band[2].elements[local_i] += Ktilde[k-1][l-1];
					else if (Delta == Nx+1)
							band[3].elements[local_i] += Ktilde[k-1][l-1];
					else{
						error("Error Indexing local K matrix out of band bounds !");
					}
				}
			}
			F[local_i] = F[local_i] + Ftilde[k];
		}
	}
}
// jiiouou
double bc_choice(double* hat_u,double u0, int subdomain_idx, int local_coor){
    int p=2; //get this from mpi
    if((subdomain_idx == 0) || (subdomain_idx == p-1))
        return u0;
    else{
        return(hat_u[local_coor]);
    }
}


void* set_boundary_values(domain* D, int subdomain_idx, vector * bands, double* F ) {
    int Nx, Ny, global_Nx;
    double u0 = 1; // move to domain
    double newF;
    double uhat[Nx*Ny];
    //uhat = internal_boundart_from_smoothed_sol();
    Nx = D->subdomains[subdomain_idx].dimX;
    Ny = D->subdomains[subdomain_idx].dimY;
    int local_boundary_element_idx;

    // Update K for Bottom Wall
    for (int i = 0; i < Nx; i++) {
        local_boundary_element_idx = i;
        bands[0].elements[local_boundary_element_idx] = 1;
        newF = bc_choice(hat_u,u0, subdomain_idx, local_boundary_element_idx);
        F[local_boundary_element_idx] = newF;
        for (int band_id = 1; band_id < 4; band_id++)
            if (local_boundary_element_idx < vector_sizes[band_id])
                bands[band_id].elements[local_boundary_element_idx] = 0;
    }
    // Update K for Top Wall
    for (int i = 0; i < Nx; i++) {
        local_boundary_element_idx = (Ny - 1) * Nx + i;
        bands[0].elements[local_boundary_element_idx] = 1;
        newF = bc_choice(hat_u,u0, subdomain_idx, local_boundary_element_idx);
        F[local_boundary_element_idx] = newF;
        for (int band_id = 1; band_id < 4; band_id++)
            if (local_boundary_element_idx < vector_sizes[band_id])
                bands[band_id].elements[local_boundary_element_idx] = 0;
    }
    // Update K for Right Wall
    for(int j=0; j<Ny; j++){
        local_boundary_element_idx = (j+1)*Nx -1 ;
        bands[0].elements[local_boundary_element_idx] = 1;
        newF = bc_choice(hat_u,u0, subdomain_idx, local_boundary_element_idx);
        F[local_boundary_element_idx] = newF;
        for(int band_id=1; band_id<4; band_id++)
            if (local_boundary_element_idx<vector_sizes[band_id])
                bands[band_id].elements[local_boundary_element_idx] = 0;
    }
    // Update K for Left Wall
    for (int j = 0; j < Ny; j++) {
        local_boundary_element_idx = j * Nx;
        bands[0].elements[local_boundary_element_idx] = 1;
        newF = bc_choice(hat_u,u0, subdomain_idx, local_boundary_element_idx);
        F[local_boundary_element_idx] = newF;
        for (int band_id = 1; band_id < 4; band_id++)
            if (local_boundary_element_idx < vector_sizes[band_id])
                    bands[band_id].elements[local_boundary_element_idx] = 0;
        }
    }