#include "domain.h"
#include <stdio.h>
#include <math.h>
#include "elements.h"
#include "error.h"
#include "vector.h"
#include "assemble.h"
#include "matrix.h"
#include "fep.h"

void assemble_local_KF(sparse_matrix* K, vector* F, domain* D,
                       int subdomain_idx)
{
  int i, k, l, t_idx, Delta;
  int local_i, local_j;
  int global_element_idx, Nv, Nx, Ny, Nt;
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

 void assemble_global_matrix(domain* D, int subdomain_idx, double** local_matrix, sparse_matrix* global_matrix
                            int vector_sizes[], int diagonal_offsets[], int diag_count)
{
  int i, j, k, l, offset, found_offset, vertex_count;
  int element_vertex_id1, element_vertex_id2, global_element_idx;
  vector* bands;
	vertex** element_vertices;

  bands = (vector*) calloc(diag_count, sizeof(vector));

  // Initialize the vector for the bands
  for (i = 0; i < diag_count; i++)
  {
		 vector_init(&bands[i], vector_sizes[i]);
	}

  #define CSD (D->subdomains)
  #define CDM (D->subdomain_vertex_map)

  // Go over each element in subdomain
	for(i = 0; i < CSD[subdomain_idx].elements_count; i++)
  {
		global_element_idx = CSD[subdomain_idx].elements[i];

    // TODO: Can we get rid of reference to triangular element
		element_vertices = _triangular_elements[global_element_idx].grid_vertex;
    vertex_count = _triangular_elements[global_element_idx].element_vertex_count;

		for(j = 1; j <= vertex_count; j++)
    {
			element_vertex_id1 =  CDM[subdomain_idx][element_vertices[j-1]->id];

      for(k = 1; k <= vertex_count; k++)
      {
				element_vertex_id2 =  CDM[subdomain_idx][element_vertices[k-1]->id];
				offset = element_vertex_id2 - element_vertex_id1;

        found_offset = 0;
        for(l = 1; l < diag_count; l++)
        {
          if(diagonal_offsets[l] == offset)
          {
            found_offset = 1;
            if(offset >= 0)
            {
              bands[l].elements[element_vertex_id1] += Ktilde[j-1][k-1];
            }
            else
            {
              bands[l].elements[element_vertex_id1 + offset] += Ktilde[j-1][k-1];
            }
          }
        }

        if(!found_offset)
        {
					printf("offset = %d \n",offset);
          error("Error indexing local matrix out of band bounds!");
				}
			}
		}
	}

  #undef CSD
  #undef CDM

  // apply the boundary conditions before converting to a sparse_symmetric_banded matrix
  boundary_op_K(bands, F, D, subdomain_idx);

  // Create sparse matrix
  // TODO: Find out what Nv is and fix it.
  sparse_matrix_banded_init(global_matrix, Nv, bands, diagonal_offsets, diag_count);

  for (i = 0; i < diag_count; i++)
  {
    vector_free(&bands[i]);
  }

  free(bands);
}


void assemble_global_load_vector(domain* D, int subdomain_idx, vector* global_load_vector)
{
  int i, j, vertex_count;
	vertex** element_vertices;

  #define CSD (D->subdomains)
  #define CDM (D->subdomain_vertex_map)
  #define LDV (*global_load_vector)

  vector = calloc(1, size(vector));

  // TODO: Need to call a function with the element_id based on which the
  // function will return the integral of v.f over the element

  // Go over each element in subdomain
	for(i = 0; i < CSD[subdomain_idx].elements_count; i++)
  {
		global_element_idx = CSD[subdomain_idx].elements[i];

    // TODO: Can we get rid of reference to triangular element
		element_vertices = _triangular_elements[global_element_idx].grid_vertex;
    vertex_count = _triangular_elements[global_element_idx].element_vertex_count;

		for(j = 1; j <= vertex_count; j++)
    {
			element_vertex_id1 =  CDM[subdomain_idx][element_vertices[j-1]->id];
			LDV.elements[element_vertex_id1] += Ftilde[j-1];
		}
	}

  #undef CSD
  #undef CDM
  #undef CDV
}


void boundary_op_local_F(vector* F, domain* D, int subdomain_idx)
{
  int i, j;
  int Nx, Ny;
  int vertex_id;

  #define CSD (D->subdomains)
  #define FE (F->elements)

  Nx = CSD[subdomain_idx].dimX;
  Ny = CSD[subdomain_idx].dimY;


  // Update F for Bottom Wall
  for (i = 0; i < Nx; i++)
  {
      vertex_id = i;
      FE[vertex_id] = get_boundary_value(D, subdomain_idx, vertex_id, bottom);
  }

  // Update F for Top Wall
  for (i = 0; i < Nx; i++)
  {
      vertex_id = (Ny - 1) * Nx + i;
      FE[vertex_id] = get_boundary_value(D, subdomain_idx, vertex_id, top);
  }

  // Update F for Right Wall
  for(j = 0; j < Ny; j++)
  {
      vertex_id = (j+1)*Nx -1 ;
      FE[vertex_id] = get_boundary_value(D, subdomain_idx, vertex_id, right);
  }

  // Update F for Left Wall
  for (j = 0; j < Ny; j++)
  {
      vertex_id = j * Nx;
      FE[vertex_id] = get_boundary_value(D, subdomain_idx, vertex_id, left);
  }

  #undef FE
  #undef CSD
}


double get_boundary_value(domain* D, int subdomain_idx, int vertex_id, int wall)
{
  vertex* grid_vertex;

  #define CSD (D->subdomains)
  #define CDM (D->subdomain_vertex_map)
  #define CSDN (D->subdomain_count_x)

  grid_vertex = CSD[subdomain_idx].subdomain_vertices[vertex_id];

  if (CSDN == 1)
  {
    return boundary_value(grid_vertex->x, grid_vertex->y, 0);
  }
  if(subdomain_idx == 0)
  {
    if (wall != right)
    {
      return boundary_value(grid_vertex->x, grid_vertex->y, 0);
    }
    else
    {
      return CSD[subdomain_idx].subdomain_solution.elements[vertex_id];
    }
  }
  else if (subdomain_idx == CSDN - 1)
  {
    if (wall != left)
    {
      return boundary_value(grid_vertex->x, grid_vertex->y, 0);
    }
    else
    {
      return CSD[subdomain_idx].subdomain_solution.elements[vertex_id];
    }
  }
  else
  {
    if ((wall == bottom) || (wall == top))
		{
		  return boundary_value(grid_vertex->x, grid_vertex->y, 0);
		}
		else
    {
      return CSD[subdomain_idx].subdomain_solution.elements[vertex_id];
    }
  }

  #undef CSD
  #undef CDM
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
