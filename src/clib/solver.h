#include <stdlib.h>
#include <stdio.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

double qinterp(double q0, double q1, double deltat, double dt) ;
int do_network_step(double *Y, double *r, double dt, double *dtrcmd, int ispecies) ;
int integrate_network(int water_rate, double *Y, double T0, double T1, double n, double metl, double UV, double dtwant, int *nstp, code_units *my_units, int ispecies, int H2_shield, double crsHI, double k24);
//int integrate_network(double *Y, double T0, double T1, double n, double metl, double UV, double dtwant, int *nstp, code_units *my_units, int ispecies);

