#include <stdlib.h>
#include <stdio.h>
//#include <parameters.h>
//#include <reactions.h>
#include <chemistry.h>

int check_solution(double *Y, double *Ylst, double *err, double *tol);

int rk45_integrate(double *Y, double *r, double *dt, double *dtrcmd);
