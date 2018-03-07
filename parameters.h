#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include "output.h"

typedef struct
{
  int outerItr;
  int innerItr;
  double absTol;
  double relTol;
} mgmres_parameters;

typedef struct
{
  mgmres_parameters mgmresParameters;
  output_processor outputProcessor;
  int maxItr;
  double solverRelTol;
} elliptic_solver_parameters;

typedef struct
{
  mgmres_parameters mgmresParameters;
  output_processor outputProcessor;
  int maxItr;
  double solverRelTol;
  double timeSpan;
  int numTimeSteps;
} parabolic_solver_parameters;

#endif
