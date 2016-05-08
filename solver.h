#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "parameters.h"

int ellipticsolver(domain* cartesian_domain, elliptic_solver_parameters solver_parameters);
int parabolicsolver(domain* cartesian_domain, parabolic_solver_parameters solver_parameters);

#endif
