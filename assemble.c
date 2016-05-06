#include "domain.h"
#include <stdio.h>
#include <math.h>
#include "elements.h"
#include "error.h"
#include "vector.h"
#include "assemble.h"
#include "matrix.h"

void assemble_local_KF(sparse_matrix* K, vector* F, domain* D,
                       int subdomain_idx)
{
  int i, k, l, t_idx, Delta;
    int local_i, local_j;
  int global_element_idx,Nv, Nx, Ny, Nt;
  vector bands[7];
  double Ktilde[3][3] = {{1.0, -0.5, -0.5},
                         {-0.5, 0.5, 0.0},
                         {-0.5, 0.0, 0.5}};
  double h = D->cartesian_grid->h;
  /*double diag= h*h*1.0/12;
  double offdiag = h*h*1.0/24 ;
  double Mtilde[3][3] = {{diag,offdiag, offdiag},
                         {offdiag,diag,offdiag},
                         {offdiag,offdiag,diag}};*/
  double vol = 1.0/6 * h * h;
  double Ftilde[3] = {vol, vol, vol};
	Nx = D->subdomains[subdomain_idx].dimX;
	Ny = D->subdomains[subdomain_idx].dimY;
	Nv = Nx* Ny;
	Nt = (Nx-1)* (Ny -1) *2;
  int vector_sizes[7] = {Nv, Nv-1, Nv-Nx, Nv-Nx-1, Nv-Nx-1, Nv-Nx, Nv-1};
  int diagonal_offsets[7] = {0, 1, Nx, Nx + 1, -Nx-1, -Nx, -1};
  vector_init(F, Nv);
	vertex** triple_vertices;

  for (i = 0; i < 7; i++)
  {
		 vector_init(&bands[i],vector_sizes[i]);
	}

	for(t_idx = 0; t_idx < Nt; t_idx++)
  {
		global_element_idx = D->subdomains[subdomain_idx].elements[t_idx];
		triple_vertices = _triangular_elements[global_element_idx].grid_vertex;

		for(k = 1; k <= 3; k++)
    {
			local_i =  D->subdomain_vertex_map[subdomain_idx][triple_vertices[k-1]->id];

      for(l = 1; l <= 3; l++)
      {
				local_j =  D->subdomain_vertex_map[subdomain_idx][triple_vertices[l-1]->id];

				Delta = local_j - local_i;
				if(Delta == 0)
        {
					bands[0].elements[local_i] += Ktilde[k-1][l-1];
        }
		    else if(Delta==1)
        {
					bands[1].elements[local_i] += Ktilde[k-1][l-1];
        }
				else if(Delta == Nx)
        {
			  	bands[2].elements[local_i] += Ktilde[k-1][l-1];
        }
				else if (Delta == Nx+1)
        {
					bands[3].elements[local_i] += Ktilde[k-1][l-1];
        }
        else if(Delta== -1)
        {
          bands[6].elements[local_i - 1] += Ktilde[k-1][l-1];
        }
        else if(Delta == -Nx)
        {
          bands[5].elements[local_i - Nx] += Ktilde[k-1][l-1];
        }
        else if (Delta == -Nx-1)
        {
          bands[4].elements[local_i - Nx - 1] += Ktilde[k-1][l-1];
        }
				else
        {
					printf("Delta= %d \n",Delta);
          error("Error Indexing local K matrix out of band bounds !");
				}
			}

			F->elements[local_i] += Ftilde[k-1];
		}
	}

  // apply the boundary conditions before converting to a sparse_symmetric_banded matrix
  boundary_op_K(bands, F, D, subdomain_idx);
  sparse_matrix_banded_init(K, Nv, bands, diagonal_offsets, 7);

  for (i = 0; i < 7; i++)
  {
    vector_free(&bands[i]);
  }
}

void boundary_op_local_F(vector* F, domain* D, int subdomain_idx)
{
  int i, j;
  int Nx, Ny;
  double u0 = 0;
  double newF;
  Nx = D->subdomains[subdomain_idx].dimX;
  Ny = D->subdomains[subdomain_idx].dimY;
  int local_boundary_element_idx;

  // Update K for Bottom Wall
  for (i = 0; i < Nx; i++)
  {
      local_boundary_element_idx = i;
      newF = bc_choice(D,u0, subdomain_idx, local_boundary_element_idx, bottom);
      F->elements[local_boundary_element_idx] = newF;
  }

  // Update K for Top Wall
  for (i = 0; i < Nx; i++)
  {
      local_boundary_element_idx = (Ny - 1) * Nx + i;
      newF = bc_choice(D,u0, subdomain_idx, local_boundary_element_idx, top);
      F->elements[local_boundary_element_idx] = newF;
  }

  // Update K for Right Wall
  for(j = 0; j < Ny; j++)
  {
      local_boundary_element_idx = (j+1)*Nx -1 ;
      newF = bc_choice(D,u0, subdomain_idx, local_boundary_element_idx, right);
      F->elements[local_boundary_element_idx] = newF;
  }

  // Update K for Left Wall
  for (j = 0; j < Ny; j++)
  {
      local_boundary_element_idx = j * Nx;
      newF = bc_choice(D,u0, subdomain_idx, local_boundary_element_idx, left);
      F->elements[local_boundary_element_idx] = newF;
  }
}

double bc_choice(domain* D,double u0, int subdomain_idx,
                 int local_coor, int where){
    int p = D->subdomain_count_x;
    if (p == 1)
    {
      return u0;
    }
    if(subdomain_idx == 0)
    {
      if (where != right)
      {
        return u0;
      }
      else
      {
          return D->subdomains[subdomain_idx].subdomain_solution.elements[local_coor];
      }
    }
    else if (subdomain_idx == p-1)
    {
      if (where != left)
      {
        return u0;
      }
      else
      {
          return D->subdomains[subdomain_idx].subdomain_solution.elements[local_coor];
      }
    }
    else
    {
      if ((where == bottom) || (where == top))
			{
			  return u0;
			}
			else
      {
        return D->subdomains[subdomain_idx].subdomain_solution.elements[local_coor];
      }
    }
}


void boundary_op_K(vector* bands, vector* F, domain* D,
                    int subdomain_idx)
{
  int i, j, band_id;
  int Nx, Ny, Nv;

  Nx = D->subdomains[subdomain_idx].dimX;
  Ny = D->subdomains[subdomain_idx].dimY;
  int local_boundary_element_idx;
  Nv = Nx* Ny;
  int vector_sizes[7] = {Nv, Nv-1, Nv-Nx, Nv-Nx-1, Nv-Nx-1, Nv-Nx, Nv-1};

  // Update K for Bottom Wall
  for (i = 0; i < Nx; i++)
  {
    local_boundary_element_idx = i;               // Variable index in the matrix
    bands[0].elements[local_boundary_element_idx] = 1;
    for (band_id = 1; band_id < 4; band_id++)
    {
      if (local_boundary_element_idx < vector_sizes[band_id])
      {
        bands[band_id].elements[local_boundary_element_idx] = 0;
      }
    }
    for (band_id = 4; band_id < 7; band_id++)
    {
      if(local_boundary_element_idx >= Nv - vector_sizes[band_id])
      {
          bands[band_id].elements[local_boundary_element_idx - Nv + vector_sizes[band_id]] = 0;
      }
    }
  }

  // Update K for Top Wall
  for (i = 0; i < Nx; i++)
  {
    local_boundary_element_idx = (Ny - 1) * Nx + i;
    bands[0].elements[local_boundary_element_idx] = 1;
    for (band_id = 1; band_id < 4; band_id++)
    {
      if (local_boundary_element_idx < vector_sizes[band_id])
      {
        bands[band_id].elements[local_boundary_element_idx] = 0;
      }
    }
    for (band_id = 4; band_id < 7; band_id++)
    {
      if(local_boundary_element_idx >= Nv - vector_sizes[band_id])
      {
          bands[band_id].elements[local_boundary_element_idx - Nv + vector_sizes[band_id]] = 0;
      }
    }
  }

  // Update K for Right Wall
  for(j = 0; j < Ny; j++)
  {
    local_boundary_element_idx = (j+1) * Nx - 1 ;
    bands[0].elements[local_boundary_element_idx] = 1;
    for (band_id = 1; band_id < 4; band_id++)
    {
      if (local_boundary_element_idx < vector_sizes[band_id])
      {
        bands[band_id].elements[local_boundary_element_idx] = 0;
      }
    }
    for (band_id = 4; band_id < 7; band_id++)
    {
      if(local_boundary_element_idx >= Nv - vector_sizes[band_id])
      {
          bands[band_id].elements[local_boundary_element_idx - Nv + vector_sizes[band_id]] = 0;
      }
    }
  }

  // Update K for Left Wall
  for (j = 0; j < Ny; j++)
  {
    local_boundary_element_idx = j * Nx;
    bands[0].elements[local_boundary_element_idx] = 1;
    for (band_id = 1; band_id < 4; band_id++)
    {
      if (local_boundary_element_idx < vector_sizes[band_id])
      {
        bands[band_id].elements[local_boundary_element_idx] = 0;
      }
    }
    for (band_id = 4; band_id < 7; band_id++)
    {
      if(local_boundary_element_idx >= Nv - vector_sizes[band_id])
      {
          bands[band_id].elements[local_boundary_element_idx - Nv + vector_sizes[band_id]] = 0;
      }
    }
  }
}
