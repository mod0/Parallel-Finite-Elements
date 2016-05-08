#include <math.h>
#include "fep.h"

#define PI 3.14159265358979323

// Define boundary operators
double boundary_value (double x, double y, double t)
{
  return 0.0;
}

// Define forcing term
double forcing_term (double x, double y, double t)
{
  return 0.0;
}

// Define initial value operators
double initial_value (double x, double y)
{
  return sin(PI * x) * sin(PI * y);
}
