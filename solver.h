#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "parameters.h"

int ellipticsolver(domain* cartesian_domain, elliptic_solver_parameters solver_parameters);

static int is_converged(domain* cartesian_domain);

static double smooth_solution_and_get_norm(domain* cartesian_domain, int idx);

#endif
