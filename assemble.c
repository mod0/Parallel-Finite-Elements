#include "domain.h"
#include <stdio.h>
#include <math.h>
#include "elements.h"
#include "error.h"
#include "vector.h"
#include "assemble.h"
#include "matrix.h"
#include "fep.h"

//TODO: Will only assemble linear lagrange basis on triangular matrices at the moment.
 void assemble_global_matrix(domain* D, int subdomain_idx, double local_matrix[][3], sparse_matrix* global_matrix,
                            int vector_sizes[], int diagonal_offsets[], int diag_count, double boundary_diagonal_value)
{
  int i, j, k, l, offset, found_offset, Nx, Ny, Nv;
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

  Nx = CSD[subdomain_idx].dimX;
  Ny = CSD[subdomain_idx].dimY;
  Nv = Nx* Ny;

  // Go over each element in subdomain
	for(i = 0; i < CSD[subdomain_idx].elements_count; i++)
  {
		global_element_idx = CSD[subdomain_idx].elements[i];

    // TODO: Can we get rid of reference to triangular element
		element_vertices = _triangular_elements[global_element_idx].grid_vertex;

		for(j = 1; j <= 3; j++)
    {
			element_vertex_id1 =  CDM[subdomain_idx][element_vertices[j-1]->id];

      for(k = 1; k <= 3; k++)
      {
				element_vertex_id2 =  CDM[subdomain_idx][element_vertices[k-1]->id];
				offset = element_vertex_id2 - element_vertex_id1;

        found_offset = 0;
        for(l = 0; l < diag_count; l++)
        {
          if(diagonal_offsets[l] == offset)
          {
            found_offset = 1;
            if(offset >= 0)
            {
              bands[l].elements[element_vertex_id1] += local_matrix[j-1][k-1];
            }
            else
            {
              bands[l].elements[element_vertex_id1 + offset] += local_matrix[j-1][k-1];
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
  apply_boundary_operator_on_matrix(D,  subdomain_idx, bands, vector_sizes,
                                    diagonal_offsets, diag_count, boundary_diagonal_value);

  // Create sparse matrix
  sparse_matrix_banded_init(global_matrix, Nv, bands, diagonal_offsets, diag_count);

  for (i = 0; i < diag_count; i++)
  {
    vector_free(&bands[i]);
  }

  free(bands);
}


void assemble_global_load_vector(domain* D, int subdomain_idx, vector* global_load_vector)
{
  int i, j, global_element_idx;
  int element_vertex_id1, Nx, Ny, Nv;
	vertex** element_vertices;

  #define CSD (D->subdomains)
  #define CDM (D->subdomain_vertex_map)
  #define LDV (*global_load_vector)

  // TODO: Need to call a function with the element_id based on which the
  // function will return the integral of v.f over the element


  Nx = CSD[subdomain_idx].dimX;
  Ny = CSD[subdomain_idx].dimY;
  Nv = Nx* Ny;

  vector_init(global_load_vector, Nv);

  // Go over each element in subdomain
	for(i = 0; i < CSD[subdomain_idx].elements_count; i++)
  {
		global_element_idx = CSD[subdomain_idx].elements[i];

    // TODO: Can we get rid of reference to triangular element
		element_vertices = _triangular_elements[global_element_idx].grid_vertex;

		for(j = 1; j <= 3; j++)
    {
			element_vertex_id1 =  CDM[subdomain_idx][element_vertices[j-1]->id];
			LDV.elements[element_vertex_id1] += triangular_element_one_point_quadrature(D, subdomain_idx, global_element_idx);
		}
	}

  #undef CSD
  #undef CDM
  #undef CDV
}


void apply_boundary_operator_on_vector(domain* D, int subdomain_idx, vector* F)
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
  #undef CSDN
}


void apply_boundary_operator_on_matrix(domain* D,  int subdomain_idx, vector* bands, int vector_sizes[],
                                      int diagonal_offsets[], int diag_count, double boundary_diagonal_value)
{
  int i, j, band_id;
  int Nx, Ny, Nv;
  int vertex_id;

  #define CSD (D->subdomains)

  Nx = CSD[subdomain_idx].dimX;
  Ny = CSD[subdomain_idx].dimY;
  Nv = Nx* Ny;

  // Update K for Bottom Wall
  for (i = 0; i < Nx; i++)
  {
    vertex_id = i;               // Variable index in the matrix

    for(band_id = 0; band_id < diag_count; band_id++)
    {
      if(diagonal_offsets[band_id] < 0 && (vertex_id >= Nv - vector_sizes[band_id]))
      {
        bands[band_id].elements[vertex_id - Nv + vector_sizes[band_id]] = 0;
      }
      else if(diagonal_offsets[band_id] > 0 && (vertex_id < vector_sizes[band_id]))
      {
        bands[band_id].elements[vertex_id] = 0;
      }
      else if(diagonal_offsets[band_id] == 0)
      {
        bands[band_id].elements[vertex_id] = boundary_diagonal_value;
      }
    }
  }

  // Update K for Top Wall
  for (i = 0; i < Nx; i++)
  {
    vertex_id = (Ny - 1) * Nx + i;

    for(band_id = 0; band_id < diag_count; band_id++)
    {
      if(diagonal_offsets[band_id] < 0 && (vertex_id >= Nv - vector_sizes[band_id]))
      {
        bands[band_id].elements[vertex_id - Nv + vector_sizes[band_id]] = 0;
      }
      else if(diagonal_offsets[band_id] > 0 && (vertex_id < vector_sizes[band_id]))
      {
        bands[band_id].elements[vertex_id] = 0;
      }
      else if(diagonal_offsets[band_id] == 0)
      {
        bands[band_id].elements[vertex_id] = boundary_diagonal_value;
      }
    }
  }

  // Update K for Right Wall
  for(j = 0; j < Ny; j++)
  {
    vertex_id = (j + 1) * Nx - 1 ;

    for(band_id = 0; band_id < diag_count; band_id++)
    {
      if(diagonal_offsets[band_id] < 0 && (vertex_id >= Nv - vector_sizes[band_id]))
      {
        bands[band_id].elements[vertex_id - Nv + vector_sizes[band_id]] = 0;
      }
      else if(diagonal_offsets[band_id] > 0 && (vertex_id < vector_sizes[band_id]))
      {
        bands[band_id].elements[vertex_id] = 0;
      }
      else if(diagonal_offsets[band_id] == 0)
      {
        bands[band_id].elements[vertex_id] = boundary_diagonal_value;
      }
    }
  }

  // Update K for Left Wall
  for (j = 0; j < Ny; j++)
  {
    vertex_id = j * Nx;

    for(band_id = 0; band_id < diag_count; band_id++)
    {
      if(diagonal_offsets[band_id] < 0 && (vertex_id >= Nv - vector_sizes[band_id]))
      {
        bands[band_id].elements[vertex_id - Nv + vector_sizes[band_id]] = 0;
      }
      else if(diagonal_offsets[band_id] > 0 && (vertex_id < vector_sizes[band_id]))
      {
        bands[band_id].elements[vertex_id] = 0;
      }
      else if(diagonal_offsets[band_id] == 0)
      {
        bands[band_id].elements[vertex_id] = boundary_diagonal_value;
      }
    }
  }

  #undef CSD
}


double triangular_element_one_point_quadrature(domain* cartesian_domain, int subdomain_idx, int element_id)
{
  double h, centroid_x, centroid_y;
  triangular_element t;

  #define CGH (cartesian_domain->cartesian_grid->h);

  h = CGH;
  t = _triangular_elements[element_id];

  centroid_x = t.grid_vertex[0]->x;
  centroid_y = t.grid_vertex[0]->y;
  centroid_x += t.grid_vertex[1]->x;
  centroid_y += t.grid_vertex[1]->y;
  centroid_x += t.grid_vertex[2]->x;
  centroid_y += t.grid_vertex[2]->y;

  centroid_x = centroid_x/3.0;
  centroid_y = centroid_y/3.0;

  #undef CGH

  return (1.0/6 * h * h * forcing_term(centroid_x, centroid_y, 0));
}
