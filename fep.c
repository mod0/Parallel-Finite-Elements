#include <math.h>
#include "fep.h"


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
  return 1.0;
}
