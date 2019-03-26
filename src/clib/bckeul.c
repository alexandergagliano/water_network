#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <err.h>
#include <bckeul.h>
#include <const.h>
#include <stdbool.h>

// compute a backward-euler time step. if newton is true (non-zero), this is a
// step of a newton-raphson iteration, in which case the abundances for the
// previous iteration should also be provided in the variable Ylst
int backward_euler_step(double *r, double *Y, double dt, double *dtrcmd, int newton, double *Ylst, int ispecies) {
  // MultiSpecies handler 
  int nSpecies   = (ispecies > 1) ? ((ispecies > 2) ? 26 : 23) : 21;
  int nReactions = (ispecies > 1) ? ((ispecies > 2) ? 54 : 48) : 24;

  //Add UV reactions!!
  nReactions += 8;


  double *F = (double *) malloc(nSpecies * sizeof(double)); // RHS (dY/dt)
  gsl_matrix *A;
  gsl_vector *b;
  double *pdelta = (double *) malloc(nSpecies * sizeof(double));

  double dtinv = 1. / dt;
  int i, j, ierr;

  // set up Jacobian matrix and initialize to zero
  double **J = (double **) malloc(nSpecies * sizeof(double *));
  for (i = 0; i < nSpecies; i++) {
    J[i] = (double *) malloc(nSpecies * sizeof(double));
    for (j = 0; j < nSpecies; j++) {
      J[i][j] = 0.;
    }
  }

  // get jacobian: use Y from last time step if BE or last iteration (Ylst) if NR
  (newton) ? get_jacobian(r, J, Ylst, ispecies) : get_jacobian(r, J, Y, ispecies);

  // clip abs. values to 1e-99
  clip_jacobian(J, ispecies);

  // alloc
  b = gsl_vector_alloc(nSpecies);
  A = gsl_matrix_alloc(nSpecies, nSpecies);

  // get RHS (dY/dt): use Y from last time step if BE or last iteration (Ylst) if NR
  (newton) ? get_f_vector(r, F, Ylst, ispecies) : get_f_vector(r, F, Y, ispecies);

  // set up matrix A and vector b for the Ax = b linear solve
  if (newton) {
    // NR: copy -dt*J to A and Ylst - Y_previous - dt * F to b
    for (i = 0; i < nSpecies; i++) {
      gsl_vector_set(b, i, Ylst[i] - Y[i] - dt * F[i]);
      for (j = 0; j < nSpecies; j++) {
        // if N-R multiply J by dt
        gsl_matrix_set(A, i, j, -dt *J[i][j]);
      }
    }
  }
  else {
    // BE: copy -J to A and F to b
    for (i = 0; i < nSpecies; i++) {
      gsl_vector_set(b, i, F[i]);
      for (j = 0; j < nSpecies; j++) {
        // if N-R multiply J by dt
        gsl_matrix_set(A, i, j, -J[i][j]);
      }
    }
  }

  /*printf("max/min A: %.5e, %.5e\n", gsl_matrix_max(A), gsl_matrix_min(A));*/

  // BE: A = I/dt - J; NR: A = I - dt*J
  if (newton) {
    for (i = 0; i < nSpecies; i++) {
      gsl_matrix_set(A, i, i, 1. + gsl_matrix_get(A, i, i));
    }
  }
  else {
    for (i = 0; i < nSpecies; i++) {
      gsl_matrix_set(A, i, i, dtinv + gsl_matrix_get(A, i, i));
    }
  }


  // solve system Ax = b for x (returned as b)
  ierr = solve_Ax_equals_b(A, b);
  if (ierr != 0) {
    printf("ERROR: Ax = b solve failed with error %d", ierr);
    *dtrcmd = .2 * dt;
    goto error;
  }

  // update the abundances
  if (newton) {
    for (i = 0; i < nSpecies; i++) {
      Ylst[i] -= gsl_vector_get(b, i);
    }
  }
  else {
    for (i = 0; i < nSpecies; i++) {
      Y[i] += gsl_vector_get(b, i);
    }
  }

  // check for negative abundances and NaNs
  double x;
  double sumM0 = 0;
  double sumM1 = 0;
  double subcycle = 0;
  for (i = 0; i < nSpecies; i++) {
    x = (newton) ? Ylst[i] : Y[i];

    if (x < 0.0) {
      // abundance is negative, reject time step and try with smaller one
      ierr = NEG_X;
      *dtrcmd = .2 * dt;
      goto error;
    }
    else if (x != x || isnan(x)) {
      // abundance is NaN, reject time step and try with smaller one
      ierr = NAN_X;
      *dtrcmd = .2 * dt;
      goto error;
    }
    
    // Check that no species changes by more than 20% in one iteration
    pdelta[i] = (Y[i] - Ylst[i]) / Ylst[i];
    if (pdelta[i] > 0.2) {subcycle = 1;}
  }

  if (subcycle) {
     ierr = DEL_X;
     *dtrcmd = 0.1 * dt;
     goto error;
  }

  //check for mass conservation!
  sumM0 = calculate_mass(Ylst);
  sumM1 = calculate_mass(Y);
  double sumFrac = sumM0/sumM1;
  if (fabs(sumFrac - 1.) > 1.e-20)
  {

    printf("Mass Conservation violation! pChange(x) = %.5e \n", sumFrac);
    *dtrcmd = 0.1 * dt;
    ierr = MASSC;
  }


  // free memory
  gsl_vector_free (b);
  gsl_matrix_free (A);
  for (i = 0; i < nSpecies; i++) 
  {
      free(J[i]);
  }
  free(J); free(F);
  free(pdelta);
  return 0;

error:
  gsl_vector_free (b);
  gsl_matrix_free (A);
  for (i = 0; i < nSpecies; i++) 
  { 
      free(J[i]);
  } 
  free(J); free(F);
  free(pdelta);
  return ierr;
}


// perform Newton-Raphson iterations using backward Euler to converge to solution
int backward_euler_NR(double *r, double *Y, double dt, double *dtrcmd, int ispecies) {
  int iter, ierr, i;
  // MultiSpecies handler 
  int nSpecies   = (ispecies > 1) ? ((ispecies > 2) ? 26 : 23) : 21;
  int nReactions = (ispecies > 1) ? ((ispecies > 2) ? 54 : 48) : 24;
  //Add UV reactions!!
  nReactions += 8;
  const int maxiter = 1.e3;
  double *Yold = (double *) malloc(nSpecies * sizeof(double)); // abundances at prev. time step
  double *Ycur = (double *) malloc(nSpecies * sizeof(double)); // abundances at current iteration
  double *Ylst = (double *) malloc(nSpecies * sizeof(double)); // abundances at last iteration
  const double tol = 1e-20, rmax = 2, rmin = 0.2;
  double eps = 0.;

  // initialise
  memcpy(Yold, Y   , sizeof(double)*nSpecies);
  memcpy(Ycur, Yold, sizeof(double)*nSpecies);
  memcpy(Ylst, Ycur, sizeof(double)*nSpecies);

  // perform newton iterations
  for (iter = 0; iter < maxiter; iter++) {

    // perform one backward euler step
    ierr = backward_euler_step(r, Yold, dt, dtrcmd, 1, Ycur, ispecies);
    if (ierr != 0) break;

    // check for convergence
    for (i = 0; i < nSpecies; i++) {
      eps = fmax(fabs((Ycur[i] - Ylst[i]) / fmax(Ylst[i], SMALL)), eps);
    }

    //check that all species values are positive
    bool allPos = 1;
    int check;
    for(check = 0; check < nSpecies; check++){
        if (Ycur[check] < 0.0){
            allPos = 0;
        }
    }

    if (allPos){
        
        //converged?
        if (eps < tol) 
	{
	    //printf("All vals are positive, doing okay!\n"); 
            ierr = 0; break; 
        }
    }

 
    // check whether max iterations have been reached
    if (iter == maxiter - 1) 
    {
        //printf("ERROR: HIT THE MAX NUMBER OF ITERATIONS IN NR!\n");
        ierr = MXITR; 
	break; 
    }
    
    // copy Ycur (current abundances) to Ylst (abundances at last iteration)
    memcpy(Ylst, Ycur, sizeof(double)*nSpecies);
  }

  // write solution
  memcpy(Y, Ycur, sizeof(double)*nSpecies);

  // recommended time step to try next
  if (ierr == 0 || ierr == MXITR) 
  {
    *dtrcmd = fmax(fmin(pow(( tol / eps ), HALF), rmax), rmin) * dt;
  }
  else 
  {
    *dtrcmd = 0.2 * dt;
  }

  // cleanup
  free(Yold); free(Ycur); free(Ylst);
  return ierr;
}


