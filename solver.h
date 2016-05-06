#ifndef __SOLVER_H__
#define __SOLVER_H__

#define RELTOL 1e-6
#define ABSTOL 1e-6

int ellipticsolver(cartesian_domain* domain);

static int is_converged(domain* cartesian_domain);

static double smooth_solution_and_get_norm(cartesian_domain* domain, int idx);

#endif
