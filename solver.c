#include <math.h>
#include <omp.h>

#include "domain.h"
#include "assemble.h"
#include "solver.h"
#include "matrix.h"
#include "output.h"
#include "error.h"
#include "fep.h"

double _triangular_local_stiffness_matrix[][3] = {{1.0, -0.5, -0.5},
                                                  {-0.5, 0.5, 0.0},
                                                  {-0.5, 0.0, 0.5}};

double _triangular_local_mass_matrix[][3] =  {{1.0/12, 1.0/24, 1.0/24},
                                              {1.0/24, 1.0/12, 1.0/24},
                                              {1.0/24, 1.0/24, 1.0/12}};

static int is_converged(domain* cartesian_domain)
{
  int i;
  int converged = 1;

#define CSDN (cartesian_domain->subdomain_count_x)
#define CSD (cartesian_domain->subdomains)
  for(i = 0; i < CSDN; i++)
    {
      converged &= (CSD[i].converged != 0);
    }
#undef CSD
#undef CSDN

  return converged;
}

static double smooth_solution_and_get_norm(domain* cartesian_domain, int idx)
{
  int i, j, count;
  double true_norm;
  double deviation_norm;
  double current_value, ghost_value, new_value;
  double weight;
#define CSDN (cartesian_domain->subdomain_count_x)
#define CSD (cartesian_domain->subdomains)

  true_norm = 0.0;
  deviation_norm = 0.0;

  if(idx < CSDN - 1)
    {
      count = 0;
      for(i = 0; i < CSD[idx].dimY; i++)
        {
          for(j = CSD[idx].dimX - CSD[idx].overlap; j < CSD[idx].dimX; j++)
            {
              weight = (CSD[idx].dimX - j - 1)/(1.0 * CSD[idx].overlap - 1);
              current_value = CSD[idx].subdomain_solution.elements[i*CSD[idx].dimX + j];
              ghost_value = CSD[idx].ghost_subdomain_right.elements[count];
              new_value = (weight * current_value) + ((1.0 - weight) * ghost_value);
              deviation_norm = deviation_norm + pow((new_value - current_value), 2);
              true_norm = true_norm + pow(new_value, 2);
              CSD[idx].subdomain_solution.elements[i*CSD[idx].dimX + j] = new_value;
              count++;
            }
        }
    }
  else
    {
      return 0.0;
    }

#undef CSD
#undef CSDN

  return sqrt(deviation_norm/true_norm);
}

int ellipticsolver(domain* cartesian_domain, elliptic_solver_parameters solver_parameters)
{
  int Nx, Ny, Nv;
  int i, itrCount, diag_count;
  vector* F_array;
  sparse_matrix* K_array;
  int ** vector_sizes_array;
  int ** diagonal_offsets_array;

#define CSDN (cartesian_domain->subdomain_count_x)
#define CSD (cartesian_domain->subdomains)
#define CG (cartesian_domain->cartesian_grid)

  // Create an array of F vectors and K arrays
  F_array = calloc(CSDN, sizeof(vector));
  K_array = calloc(CSDN, sizeof(sparse_matrix));
  vector_sizes_array = calloc(CSDN, sizeof(int*));
  diagonal_offsets_array = calloc(CSDN, sizeof(int*));

  diag_count = 7;
  for(i = 0; i < CSDN; i++)
    {
      vector_sizes_array[i] = calloc(diag_count, sizeof(int));
      diagonal_offsets_array[i] = calloc(diag_count, sizeof(int));
    }

#pragma omp parallel for private(i, Nx, Ny, Nv)
  for(i = 0 ; i < CSDN; i++)
    {
      Nx = CSD[i].dimX;
      Ny = CSD[i].dimY;
      Nv = Nx * Ny;

      vector_sizes_array[i][0] = Nv - Nx - 1;
      vector_sizes_array[i][1] = Nv - Nx;
      vector_sizes_array[i][2] = Nv - 1;
      vector_sizes_array[i][3] = Nv;
      vector_sizes_array[i][4] = Nv - 1;
      vector_sizes_array[i][5] = Nv - Nx;
      vector_sizes_array[i][6] = Nv - Nx - 1;


      diagonal_offsets_array[i][0] = - Nx - 1;
      diagonal_offsets_array[i][1] = - Nx;
      diagonal_offsets_array[i][2] = - 1;
      diagonal_offsets_array[i][3] = 0;
      diagonal_offsets_array[i][4] = 1;
      diagonal_offsets_array[i][5] = Nx;
      diagonal_offsets_array[i][6] = Nx + 1;


      // Assemble the K matrix
      assemble_global_matrix(cartesian_domain, i, _triangular_local_stiffness_matrix, &(K_array[i]),
                             vector_sizes_array[i], diagonal_offsets_array[i], diag_count, 1.0);

      // Assemble the load vector_init
      assemble_global_load_vector(cartesian_domain, i, &(F_array[i]));
    }

  itrCount = 0;

  do
    {
      printf("Elliptic Solver: Iteration count %d\n", itrCount);

#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++)
        {
          // Apply the boundary condition on K and F
          // TODO: This is redundant operation, F is fixed
          apply_boundary_operator_on_vector(cartesian_domain, i, &(F_array[i]));

          // Mark the subdomain as not converged
          CSD[i].converged = 0;

          // Solve each subdomain
          mgmres(&(K_array[i]),  &(CSD[i].subdomain_solution), &(F_array[i]), solver_parameters.mgmresParameters);

          // Send information left
          copy_overlap_to_adjacent_neighbours_ghost(cartesian_domain, i, -1);
        }

#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++) {
        // Smooth and compute the norm
        double rel_error = smooth_solution_and_get_norm(cartesian_domain, i);

        // Check for tolerance and mark as converged if converged
        if(rel_error < solver_parameters.solverRelTol)
          {
            CSD[i].converged = 1;
          }

        // Send information right
        copy_overlap_to_adjacent_neighbours_ghost(cartesian_domain, i, 1);
      }

#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++) {
        // Copy information from the left ghost cell of the current subdomain
        copy_from_my_ghost_cell(cartesian_domain, i, -1);

        // Write to file
        solver_parameters.outputProcessor(cartesian_domain, i, itrCount, write_output_for_vertex);
      }

      if(++itrCount > solver_parameters.maxItr)
        {
          warn("The elliptic solver has exceeded the number of maximum iterations");
          break;
        }
    } while(!is_converged(cartesian_domain));

  if(K_array != NULL)
    {
#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++)
        {
          sparse_matrix_free(&(K_array[i]));
        }

      free(K_array);
    }

  if(F_array != NULL)
    {
#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++)
        {
          vector_free(&(F_array[i]));
        }

      free(F_array);
    }

  if(vector_sizes_array != NULL)
    {
#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++)
        {
          free(vector_sizes_array[i]);
        }

      free(vector_sizes_array);
    }


  if(diagonal_offsets_array != NULL)
    {
#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++)
        {
          free(diagonal_offsets_array[i]);
        }

      free(diagonal_offsets_array);
    }

#undef CSD
#undef CSDN
#undef CG
  return 0;
}

int parabolicsolver(domain* cartesian_domain, parabolic_solver_parameters solver_parameters)
{
  int i,j, pseudoItrCount, timeItrCount;
  int diag_count, Nx, Ny, Nv;
  vector* F_array;
  sparse_matrix* LHS_array;
  sparse_matrix* RHS_array;
  double h, dt, x, y;
  double lhs[3][3], rhs[3][3];
  int ** vector_sizes_array;
  int ** diagonal_offsets_array;

#define CSDN (cartesian_domain->subdomain_count_x)
#define CSD (cartesian_domain->subdomains)
#define CG (cartesian_domain->cartesian_grid)

  F_array = calloc(CSDN, sizeof(vector));
  LHS_array = calloc(CSDN, sizeof(sparse_matrix));
  RHS_array = calloc(CSDN, sizeof(sparse_matrix));
  vector_sizes_array = calloc(CSDN, sizeof(int*));
  diagonal_offsets_array = calloc(CSDN, sizeof(int*));


  dt = solver_parameters.timeSpan/solver_parameters.numTimeSteps;
  h = CG->h;
  for (i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
        {
          lhs[i][j] = h*h*_triangular_local_mass_matrix[i][j] + dt/2 * _triangular_local_stiffness_matrix[i][j];
          rhs[i][j] = h*h*_triangular_local_mass_matrix[i][j] - dt/2 * _triangular_local_stiffness_matrix[i][j];
        }
    }

  diag_count = 7;
  for(i = 0; i < CSDN; i++)
    {
      vector_sizes_array[i] = calloc(diag_count, sizeof(int));
      diagonal_offsets_array[i] = calloc(diag_count, sizeof(int));
    }

#pragma omp parallel for private(i, Nx, Ny, Nv, x, y, j)
  for(i = 0 ; i < CSDN; i++)
    {
      Nx = CSD[i].dimX;
      Ny = CSD[i].dimY;
      Nv = Nx * Ny;

      for(j = 0; j < Nx * Ny; j++)
        {
          x = CSD[i].subdomain_vertices[j]->x;
          y = CSD[i].subdomain_vertices[j]->y;
          CSD[i].subdomain_solution.elements[j] = initial_value(x, y);
        }

      vector_sizes_array[i][0] = Nv - Nx - 1;
      vector_sizes_array[i][1] = Nv - Nx;
      vector_sizes_array[i][2] = Nv - 1;
      vector_sizes_array[i][3] = Nv;
      vector_sizes_array[i][4] = Nv - 1;
      vector_sizes_array[i][5] = Nv - Nx;
      vector_sizes_array[i][6] = Nv - Nx - 1;


      diagonal_offsets_array[i][0] = - Nx - 1;
      diagonal_offsets_array[i][1] = - Nx;
      diagonal_offsets_array[i][2] = - 1;
      diagonal_offsets_array[i][3] = 0;
      diagonal_offsets_array[i][4] = 1;
      diagonal_offsets_array[i][5] = Nx;
      diagonal_offsets_array[i][6] = Nx + 1;

      // Assemble the K and M matrix toghether in LHS_array and RHS_array
      // Assemble the load vector_init
      assemble_global_matrix(cartesian_domain, i, lhs, &(LHS_array[i]),
                             vector_sizes_array[i], diagonal_offsets_array[i], diag_count, 1.0);
      assemble_global_matrix(cartesian_domain, i, rhs, &(RHS_array[i]),
                             vector_sizes_array[i], diagonal_offsets_array[i], diag_count, 1.0);

      vector_init(&(F_array[i]),Nv);
      sparse_matrix_vector_multiply(&(RHS_array[i]), &(CSD[i].subdomain_solution), &(F_array[i]));
    }



  for (timeItrCount=0; timeItrCount<solver_parameters.numTimeSteps; timeItrCount++)
    {
      printf("Unsteady Elliptic Solver: Time Iteration count %d\n", timeItrCount);
      pseudoItrCount = 0;
      do
        {
          printf("Unsteady Elliptic Solver: Inner Iteration count %d\n", pseudoItrCount);
#pragma omp parallel for private(i)
          for (i = 0; i < CSDN; i++) {
            //?? Apply the boundary condition F

            // Mark the subdomain as not converged
            CSD[i].converged = 0;

            // Apply boundary op on solution
            apply_boundary_operator_on_vector(cartesian_domain, i, &(F_array[i]));

            // Solve each subdomain
            mgmres(&(LHS_array[i]), &(CSD[i].subdomain_solution), &(F_array[i]), solver_parameters.mgmresParameters);

            // Send information left
            copy_overlap_to_adjacent_neighbours_ghost(cartesian_domain, i, -1);
          }

#pragma omp parallel for private(i)
          for (i = 0; i < CSDN; i++) {
            // Smooth and compute the norm
            double rel_error = smooth_solution_and_get_norm(cartesian_domain, i);

            // Check for tolerance and mark as converged if converged
            if (rel_error < solver_parameters.solverRelTol) {
              CSD[i].converged = 1;
            }

            // Send information right
            copy_overlap_to_adjacent_neighbours_ghost(cartesian_domain, i, 1);
          }

#pragma omp parallel for private(i)
          for (i = 0; i < CSDN; i++) {
            // Copy information from the left ghost cell of the current subdomain
            copy_from_my_ghost_cell(cartesian_domain, i, -1);
          }

          if (++pseudoItrCount > solver_parameters.maxItr) {
            printf("Elliptic solver did not converge at T= %2.2f\n", timeItrCount * dt);
            break;
          }
        } while (!is_converged(cartesian_domain));

#pragma omp parallel for private(i)
      for (i = 0; i < CSDN; i++) {
        // Write to file
        solver_parameters.outputProcessor(cartesian_domain, i, timeItrCount, write_output_for_vertex);
        // vn = RHS_array*Cn;
        sparse_matrix_vector_multiply(&(RHS_array[i]), &(CSD[i].subdomain_solution), &(F_array[i]));
      }
    }

  // Deallocation

  if(LHS_array != NULL)
    {
#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++)
        {
          sparse_matrix_free(&(LHS_array[i]));
        }

      free(LHS_array);
    }

  if(RHS_array != NULL)
    {
#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++)
        {
          sparse_matrix_free(&(RHS_array[i]));
        }

      free(RHS_array);
    }

  if(F_array != NULL)
    {
#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++)
        {
          vector_free(&(F_array[i]));
        }

      free(F_array);
    }


  if(vector_sizes_array != NULL)
    {
#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++)
        {
          free(vector_sizes_array[i]);
        }

      free(vector_sizes_array);
    }


  if(diagonal_offsets_array != NULL)
    {
#pragma omp parallel for private(i)
      for(i = 0 ; i < CSDN; i++)
        {
          free(diagonal_offsets_array[i]);
        }

      free(diagonal_offsets_array);
    }


#undef CSD
#undef CSDN
#undef CG
  return 0;
}
