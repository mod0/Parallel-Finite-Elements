#include <math.h>

#include "domain.h"
#include "assemble.h"
#include "solver.h"
#include "matrix.h"
#include "output.h"
#include "error.h"

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
        weight = (CSD[idx].dimX - j)/(1.0 * CSD[idx].overlap);
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
  int i, itrCount;
  double rel_error;
  vector* F_array;
  sparse_matrix* K_array;

  #define CSDN (cartesian_domain->subdomain_count_x)
  #define CSD (cartesian_domain->subdomains)

  // Create an array of F vectors and K arrays
  F_array = calloc(CSDN, sizeof(vector));
  K_array = calloc(CSDN, sizeof(sparse_matrix));

  for(i = 0 ; i < CSDN; i++)
  {
    // Assemble the K matrix
    // Assemble the load vector_init
    assemble_local_KF(&(K_array[i]), &(F_array[i]), cartesian_domain, i);
  }

  itrCount = 0;

  do
  {
    for(i = 0 ; i < CSDN; i++)
    {
      // Apply the boundary condition on K and F
      boundary_op_local_F(&(F_array[i]), cartesian_domain, i);

      // Mark the subdomain as not converged
      CSD[i].converged = 0;

      // Solve each subdomain
      mgmres(&(K_array[i]),  &(CSD[i].subdomain_solution), &(F_array[i]), solver_parameters.mgmresParameters);

      // Send information left
      copy_overlap_to_adjacent_neighbours_ghost(cartesian_domain, i, -1);

      // Smooth and compute the norm
      rel_error = smooth_solution_and_get_norm(cartesian_domain, i);

      // Check for tolerance and mark as converged if converged
      if(rel_error < solver_parameters.solverRelTol)
      {
        CSD[i].converged = 1;
      }

      // Send information right
      copy_overlap_to_adjacent_neighbours_ghost(cartesian_domain, i, 1);

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

  #undef CSD
  #undef CSDN

  return 0;
}
