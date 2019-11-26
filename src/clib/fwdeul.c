#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fwdeul.h>
#include <const.h>
#include <rates.h>
#include <jac_rhs.h>
#include <chemistry.h>
#include <stdbool.h>
#include <err.h>

// An implementation of the forward euler method. Because it's so simple, it's very fast. This also
// means it's wildly inaccurate and should probably only ever be used for a first approximation in
// a more stable or accurate scheme.
int forward_euler_integrate(double *r,double *Y,double dt,int ispecies,double UV,int water_rates){
  // Since nSpecies never changes, we will never need to reallocate memory. We'll do the dirty
  // thing of never releasing the memory and just crossing our fingers that the OS will take
  // care of it.
  static double *rhs;
  static int first = 1;
  if(first){
    rhs = (double*)malloc(nSpecies * sizeof(double));
    first = 0;
  }

  // Initialize our right-hand side.
  memset(rhs,0,nSpecies*sizeof(double));
  //get_f_vector(r, rhs, Y, ispecies, water_rates, (int) UV);
  cal_rhs(rhs, Y, my_reactions, r);


  // Apply the Euler step.
  for(int i = 0; i < nSpecies; i++){
    Y[i] = Y[i] + dt*rhs[i];
  }

  // I don't immediately know of any reason this method should fail, so it always returns success.
  // If someone else sees a glaring omission, they can add some failure case.
  return ALLOK;
}
