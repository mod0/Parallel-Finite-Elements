#ifndef __FEP_H__
#define __FEP_H__

// Define boundary operators
double boundary_value (double x, double y, double t);

// Define forcing term
double forcing_term (double x, double y, double t);

// Define initial value operators
double initial_value (double x, double y);

#endif
