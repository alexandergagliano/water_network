#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <err.h>
#include <bckeul.h>
#include <const.h>
//#include <icfutils.h>
//#include <parameters.h>
//#include <reactions.h>
#include <rates.h>
#include <jac_rhs.h>
#include <linalg.h>
#include <chemistry.h>
//#include <icfutils.h>


// rk4/5 cash-karp coefficients
double b21 = 1. / 5., 
       b31 = 3. / 40., 
       b32 = 9. / 40., 
       b41 = 3. / 10., 
       b42 = -9. / 10., 
       b43 = 6. / 5., 
       b51 = -11. / 54., 
       b52 = 5. * 0.5, 
       b53 = -70. / 27., 
       b54 = 35. / 27., 
       b61 = 1631. / 55296., 
       b62 = 175 / 512., 
       b63 = 575. / 13824., 
       b64 = 44275. / 110592., 
       b65 = 253. / 4096., 
       c1  = 37. / 378.,
       c1s = 2825. / 27648., 
       c2  = 0.,
       c2s = 0., 
       c3  = 250. / 621.,
       c3s = 18575. / 48384., 
       c4  = 125. / 594.,
       c4s = 13525. / 55296., 
       c5  = 0.,
       c5s = 277. / 14336., 
       c6  = 512. / 1771.,
       c6s = 1. / 4.;

double tiny = 1.e-99;

int check_solution(double *Y, double *Ylst, double *err, double *tol, int nSpecies) {
  double *X = (double *) malloc(nSpecies * sizeof(double)); // mass fractions

  double minY = Y[0];
  for (int i = 1; i < nSpecies; i++)
  {
     if (Y[i] < minY) { minY = Y[i]; }
  }
  double neg_tol,pos_tol = 1.e-27;
  double Ynegmax = fabs(minY) ? (minY < 0.) : 0.;
  //double Ynegmax = (minY < 0.) ? fabs(minY) : 0.;
  //printf("Ynegmax = %.2e\n", Ynegmax);
  int ierr;

/*
  calculate_mass_fractions(Y, X, nd, p);
  double resid = fabs(1. - sum(X,p->n));

  ierr = ALLOK;
*/
  if (Ynegmax >= neg_tol) {
    // solution converged but with "large" negative abundances, error is no
    // longer grd but the amount by which the min. abundance is below the
    // -SMALL threshold
    //*err = Ynegmax;
    ierr = NEG_X;
    goto error;
  }
  // check for negative abundances and NaNs
  double x;
  double sumM0 = 0;
  double sumM1 = 0;
  double subcycle = 0;
  for (int i = 0; i < nSpecies; i++) {
     x = Y[i];
     if (x != x || isnan(x)) {
       // abundance is NaN, reject time step and try with smaller one
       ierr = NAN_X;
       goto error;
     } 
  }
  sumM0 = calculate_mass(Ylst);
  sumM1 = calculate_mass(Y);
  if (fabs(sumM0 - sumM1) > pos_tol)
  {
    //printf("Mass Conservation violation! pChange(x) = %.5e \n", sumFrac);
    ierr = MASSC;
    goto error;
  } else {
    // check for convergence of abundances
    double eps = tiny;
    for (int i = 0; i < nSpecies; i++) {
      eps = fmax(fabs((Y[i] - Ylst[i]) / fmax(Ylst[i], SMALL)), eps);
    }

    *err = eps; *tol = pos_tol;
    ierr = (*err < *tol) ? 0 : E_TOL;
    //printf("eps = %.2e \n", eps);
  }
  error:
     free (X);
     return ierr;
}


int rk45_integrate(double *Y, double *r, double *dt, double *dtrcmd, int ispecies) {
  // MultiSpecies handler 
  int nSpecies   = (ispecies > 1) ? ((ispecies > 2) ? 26 : 23) : 21;
  int nReactions = (ispecies > 1) ? ((ispecies > 2) ? 54 : 48) : 24; 
  // RHS and initialize to zero
  double rmax = 2; double rmin = 0.2;
  double *dydt; // RHS (dY/dt)
  dydt = (double*) malloc(nSpecies * sizeof(double));
  memset(dydt, 0., nSpecies*sizeof(double*));

  // k1, k2, k3, k4, k5 and k6 or the cash-carp method - intermediate solutions
  double *k1 = (double*) malloc(nSpecies*sizeof(double));
  double *k2 = (double*) malloc(nSpecies*sizeof(double));
  double *k3 = (double*) malloc(nSpecies*sizeof(double));
  double *k4 = (double*) malloc(nSpecies*sizeof(double));
  double *k5 = (double*) malloc(nSpecies*sizeof(double));
  double *k6 = (double*) malloc(nSpecies*sizeof(double));

  // entry abundances Y0, 4th order solution Y4 and 5th order solution Y5
  double *Y0 = (double*) malloc(nSpecies*sizeof(double));
  memcpy(Y0, Y, sizeof(double)*nSpecies);
  double *Y4 = (double*) malloc(nSpecies*sizeof(double));
  double *Y5 = (double*) malloc(nSpecies*sizeof(double));

  // get RHS (dY/dt) at entry abundances Y0 and calculate k1 from it, and new
  // abundances from that
  //calculate_dydt(r, dydt, Y0, p);
  get_f_vector(r, dydt, Y0, ispecies);
  for (int i = 0; i < nSpecies; i++) {
    k1[i] = *dt * dydt[i];
    Y[i] = Y0[i] + (b21 * k1[i]);
  }

  // use k1 to get new abundances; get k2 from those
  //calculate_dydt(r, dydt, Y, p);
  get_f_vector(r, dydt, Y, ispecies);
  for (int i = 0; i < nSpecies; i++) {
    k2[i] = *dt * dydt[i];
    Y[i] = Y0[i] + (b31 * k1[i] + b32 * k2[i]);
  }

  // use k2 to get k3
  //calculate_dydt(r, dydt, Y, p);
  get_f_vector(r, dydt, Y, ispecies);
  for (int i = 0; i < nSpecies; i++) {
    k3[i] = *dt * dydt[i];
    Y[i] = Y0[i] + (b41*k1[i] + b42*k2[i] + b43*k3[i]);
  }

  // use k3 to get k4
  //calculate_dydt(r, dydt, Y, p);
  get_f_vector(r, dydt, Y, ispecies);
  for (int i = 0; i < nSpecies; i++) {
    k4[i] = *dt * dydt[i];
    Y[i] = Y0[i] + (b51*k1[i] + b52*k2[i] + b53*k3[i] + b54*k4[i]);
  }

  // use k4 to get k5
  //calculate_dydt(r, dydt, Y, p);
  get_f_vector(r, dydt, Y, ispecies);
  for (int i = 0; i < nSpecies; i++) {
    k5[i] = *dt * dydt[i];
    Y[i] = Y0[i] + (b61*k1[i] + b62*k2[i] + b63*k3[i] + b64*k4[i] + b65*k5[i]);
  }

  // get k6
  //calculate_dydt(r, dydt, Y, p);
  get_f_vector(r, dydt, Y, ispecies);
  for (int i = 0; i < nSpecies; i++) {
    k6[i] = *dt * dydt[i];
  }

  double err, tol;
  int ierr;
  // calculate 4th and 5th order solutions and set result Y to 5th order solution
  for (int i = 0; i < nSpecies; i++) {
    Y5[i] = Y0[i] + c1 *k1[i] + c2 *k2[i] + c3 *k3[i] + c4 *k4[i] + c5 *k5[i] + c6 *k6[i];
    Y4[i] = Y0[i] + c1s*k1[i] + c2s*k2[i] + c3s*k3[i] + c4s*k4[i] + c5s*k5[i] + c6s*k6[i];
    Y[i] = Y5[i];
  }

  ierr = check_solution(Y4, Y5, &err, &tol, nSpecies); 
  //ierr = 0;

  // cleanup
  free(dydt);
  free(k1); free(k2); free(k3); free(k4); free(k5); free(k6); 
  free(Y0); free(Y4); free(Y5); 
  
  return ierr;
}
