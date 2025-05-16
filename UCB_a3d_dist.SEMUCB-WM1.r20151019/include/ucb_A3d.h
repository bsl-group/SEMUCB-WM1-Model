#ifndef _UCB_A3D_H_

#define _UCB_A3D_H_

void init_ucb_ref_();
void init_ucb_a3d_();
void init_moho_();

void get_ucb_ref_(double *r, double *rho, double *vpv, double *vsv,
                  double *qk, double *qmu, double *vph, double *vsh, 
                  double *eta);

void get_a3d_perturbation_(double *theta, double *phi, double *r, char *desc, double *dp);


void get_moho_radius_(double *lat, double *lon, double *radius);

#endif
