#ifndef __FE_H__
#define __FE_H__

// Define boundary operators
double (boundary_value_cb*) (double x, double y);

// Define forcing term
double (forcing_term_cb*) (double x, double y);

// Define initial value operators
double (initial_value_cb*) (double x, double y);

#endif
