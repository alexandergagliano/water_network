#include <stdlib.h>
#include <stdio.h>

void clip_jacobian(double **J, int ispecies);
void get_f_vector(double *rate, double *F, double *Y, int ispecies);
void get_jacobian(double rate[], double **J, double Y[], int ispecies);
