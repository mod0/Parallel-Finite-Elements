#ifndef __ASSEMBLE_H__
#define __ASSEMBLE_H__

void* assemble(domain* cartesian_domain);

void* assemble_local_K(domain* cartesian_domain, double* Ktilde, double* Ftilde, int subdomain_idx);

void* set_boundary_values(domain* cartesian_domain, int subdomain_idx, vector * bands, double* F );

double bc_choice(double* hat_u, double u0, int subdomain_idx, int local_coor);

#endif
