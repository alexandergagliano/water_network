#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <solver.h>
#include <bckeul.h>
#include <err.h>
#include <const.h>
#include <chemistry.h>
#include <calc_mass.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

int do_network_step(double *Y, double *r, double dt, double *dtrcmd, int ispecies) {
  // backward-Euler
  int ierr = backward_euler_step(r, Y, dt, dtrcmd, 1, Y, ispecies);
  //int ierr = backward_euler_NR(r, Y, dt, dtrcmd, ispecies);
  //int ierr = backward_euler_step(r, Y, dt, dtrcmd, 1, Y, ispecies);

  //printf("In do_network_step, Y[0] = %.2e\n", Y[0]);

  //int ierr = rk45_integrate(Y, r, dt, dtrcmd, ispecies);
  return ierr;
}

int integrate_network(int water_rate, double *Y, double T0, double T1, double n, double metl, double UV, double dtwant, int *nstp, code_units *my_units, int ispecies, int H2_shield, double crsHI, double k24)
  {
  // MultiSpecies handler 
  int nSpecies = (ispecies > 1) ? ((ispecies > 2) ? 26 : 23) : 21;
  int nReactions = (ispecies > 1) ? ((ispecies > 2) ? 54 : 48) : 24;

  //Add UV rates! 
  nReactions += 8;

  double *Ylst = (double *) malloc(nSpecies * sizeof(double)); // abundances at last iteration
  int s;
  memcpy(Ylst, Y, sizeof(double) * nSpecies);
  double dttry, dtdone, dtrcmd;
  const int max_substp = 5.e2; // max. number of substeps
  int j, ierr;
  double *Y_old = (double *) malloc(nSpecies * sizeof(double));

  // init
  double *r = (double *) malloc(62 * sizeof(double));

  // We try to jump the entire time interval, instead of assuming subcycling
  dttry = dtwant ;
  dtdone = 0.; 
  dtrcmd = 1.e99; 
  memcpy(Y_old, Y, sizeof(double) * nSpecies);

  // substep loop
  for (j = 1; j < max_substp; j++) {
    
    get_rates(water_rate, r, T0, n, T0, metl, UV, my_units, ispecies, Y, H2_shield, crsHI, k24);
    //printf("In solver, T0 = %.2f\n", T0);
    //get_rates(r, T0, n, T0, metl, UV, my_units);


    // perform one step integration
    ierr = do_network_step(Y, r, dttry, &dtrcmd, ispecies);
    //check success of step
    if (ierr == 0){
      // success
      dtdone += dttry;
      if (dtdone == dtwant ) {
        //cycle complete
        ierr = 0;
        break;
      } else{
        dttry = fmin(dttry, dtwant - dtdone);
      }
      memcpy(Y_old, Y, sizeof(double) * nSpecies);
    }
    else {
      //failure
      memcpy(Y, Y_old, sizeof(double) * nSpecies);
      dtrcmd = 0.1 * dttry;
      dttry = fmin(dtrcmd, dtwant - dtdone);
    } 
    //dttry = fmin(dtrcmd, dtwant - dtdone);

    // check for time step underflow
    if (dttry < SMALL) {
      dttry = SMALL;
      ierr = SMLDT;
      break;
    }
    //maximum number of iterations reached
    if ( j == max_substp - 1) {
      //Currently - DON'T reject entire timestep if not finished
      ierr = MXSTP;
      break;
    }
  }
  free(r);
  free(Ylst);
  free(Y_old);
  return ierr;
}
