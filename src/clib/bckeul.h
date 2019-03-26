#include <stdlib.h>
#include <stdio.h>

int backward_euler_step(double *r, double *Y, double dt, double *dtrcmd, int newton, double *Ylst, int ispecies);
