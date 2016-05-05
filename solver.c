#include "grid.h"
#include "domain.h"
#include "elements.h"
#include "fep.h"
#include "assemble.h"
#include "solver.h"
#include "matrix.h"

int ellipticsolver(domain* cartesian_domain)
{
  int i;
  vector* F_array;
  sparse_matrix* K_array;

  #define CSDN (cartesian_domain->subdomain_count_x)
  #define CSD (cartesian_domain->cartesian_subdomain)

  // Create an array of F vectors and K arrays
  F_array = calloc(CSDN, sizeof(vector));
  K_array = calloc(CSDN, sizeof(sparse_matrix));

  for(i = 0 ; i < CSDN; i++)
  {
    // Assemble the K matrix
    // Assemble the load vector_init
    assemble_local_KF(K_array[i], F_array[i], cartesian_domain, i);

    // Apply the boundary condition on K and F
  }

  do
  {
    for(i = 0 ; i < ;)
    {
      // Mark the subdomain as not converged

      // Solve each subdomain

      // Send information left

      // Smooth and compute the norm

      // Check for tolerance and mark as converged if converged

      // Send information right

      // Write to file
    }
  } while(!is_converged(cartesian_domain));

  #undef CSD
  #undef CSDN

  return 0;
}

int is_converged(domain* cartesian_domain)
{
  int i;
  int converged = 1;

  #define CSDN (cartesian_domain->subdomain_count_x)
  #define CSD (cartesian_domain->cartesian_subdomain)
  for(i = 0; i < CSDN; i++)
  {
    converged &= (CSD[i].converged != 0);
  }
  #undef CSD
  #undef CSDN

  return converged;
}
