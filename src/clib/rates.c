#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <chemistry.h>
#include <rates.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"


/* ---------------------------------------------------
 *
 *                   Water Network
 *               Last Modified: 02/20/19
 *                Author: Alex Gagliano
 *
 * ---------------------------------------------------
 * The following lists the chemical reactions first 
 * defined in Omukai et al. (2005) for the creation
 * of water. Some of these are already defined in 
 * Grackle; where this is the case, we take the Grackle
 * rates and zero the reactions in this file. 
 * The table is constructed as follows:
 *
 * ------- 01/02/03:    a    + b   ->   c   +  d
 *
 * Where the reaction occurs between a and b
 * and produces c and d.
 *
 * 01 gives the listing of the reaction in Omukai (2005),
 * 02 gives the reaction number in Grackle (2016), 
 * and 03 gives the reaction number in Bialy & 
 * Sternberg (2018). If an entry is missing, it means this
 * reaction is not used or given in the paper. 
 *
 * Note: reaction 8 is H2 formation on dust grains!
 *
 * --------------Hydrogen Reactions------------------
 *
 * ------- H1/02/:        HII   + e   -> HI    + phot
 * ------- H2/01/:        HI    + e   -> HII   + phot
 * ------- H3/08/:        HM    + HI  -> H2I   + e
 * ------- H4/13/:        H2I   + HI  -> HI    + 2HI
 * ------- H5/22/:        HI    + 2HI -> H2I   + HI
 * ------- H6/21/:        2HI   + H2I -> H2I   + H2I
 * ------- H7//:          H2I   + H2I -> 2HI   + H2I
 * ------- H8/h2dust/:    2HI   + gr  -> H2I
 *
 * --------------Deuterium Reactions------------------
 *
 * ------- D1/50/:    DI    + HII -> DII   + HI
 * ------- D2/51/:    DII   + HI  -> DI    + HII
 * ------- D3/54/:    DI    + H2I -> HI    + HDI
 * ------- D4/52/:    DII   + H2I -> HII   + HDI
 * ------- D5/55/:    HDI   + HI  -> H2I   + DI
 * ------- D6/53/:    HDI   + HII -> H2I   + DII
 *
 * --------------Metal Reactions------------------
 *
 * ------- Z01//:    OII   + HI   -> HII   + OI
 * ------- Z02//:    HII   + OI   -> OII   + HI
 * ------- Z03//:    OII   + H2I  -> OHII  + HI
 * ------- Z04//:    OHII  + H2I  -> H2OII + HI
 * ------- Z05//:    H2OII + H2I  -> H3OII + HI
 * ------- Z06//:    H2OII + e    -> OI    + H2I
 * ------- Z07//:    H2OII + e    -> OHI   + HI
 * ------- Z08//:    H3OII + e    -> OHI   + 2HI
 * ------- Z09//:    H3OII + e    -> H2OI  + HI
 * ------- Z10//:   OI    + HI   -> OHI   + phot
 * ------- Z11//:   OI    + H2I  -> OHI   + HI
 * ------- Z12//:   H2I   + OHI  -> H2OI  + HI
 * ------- Z13//:   OHI   + OHI  -> H2OI  + OI
 * ------- Z14//:   HII   + OHI  -> OHII  + HI 
 * ------- Z15//:   HII   + H2OI -> H2OII + HI
 * ------- Z16//:   HI    + OHI  -> H2I   + OI
 * ------- Z17//:   HI    + H2OI -> OHI   + H2I
 * ------- Z18//:   OI    + OI   -> O2I   + phot
 * ------- Z19//:   OI    + OHI  -> O2I   + HI
 * ------- Z20//:   HI    + O2I  -> OHI   + OI
 * ------- Z21//:   HII   + O2I  -> O2II  + HI
 * ------- Z22//:   O2II  + e    -> OI    + OI
 * ------- Z23//:   OI    + CHI  -> COI   + HI
 * ------- Z24//:   CI    + OHI  -> COI   + HI
 * ------- Z25//:   CI    + O2I  -> COI   + OI
 * ------- Z26//:   CII   + O2I  -> OII   + CO
 * ------- Z27//:   OHI   + COI  -> CO2I  + HI
 * ------- Z28//:   CII   + e    -> CI    + phot
 * ------- Z29//:   CII   + OHI  -> COII  + HI
 * ------- Z30//:   COII  + HI   -> HII   + COI
 * ------- Z31//:   CI    + HI   -> CHI   + phot
 * ------- Z32//:   CI    + H2I  -> CHI   + HI
 * ------- Z33//:   HI    + CHI  -> CI    + H2I
 * ------- Z34//:   H2I   + CHI  -> CH2I + HI
 * ------- Z35//:   HI    + CH2I -> CHI  + H2I
 * ------- Z36//:   H2I   + CH2I -> CH3I + HI
 * ------- Z37//:   HI    + CH3I -> CH2I + H2I
 * ------- Z38//:   H2I   + CH3I -> CH4I + HI
 * ------- Z39//:   HI    + CH4I -> H2I  + CH3I 
 * ------- Z40//:   H2I   + CI   -> CH2I + phot
 *
 * ---------------------------------------------------
 */

void get_rates(int water_rates, double *rate, double temp, double n, double t_dust, double Z, double UV, code_units *my_units, int ispecies, double Y[], int H2_shield, double crsHI, double k24) 
//void get_rates(double *rate, double temp, double n, double t_dust, double Z, double UV, code_units *my_units)
{
  
  double uxyz = my_units->length_units; 
  double utim = my_units->time_units;
  double urho = my_units->density_units;
  double uaye = my_units->a_units;
  double aye = my_units->a_value;

  double mh     = 1.67262171e-24;
  double tbase1 = utim;
  double xbase1 = uxyz/(aye*uaye);  // uxyz is [x]*a     
  double dbase1 = urho * pow((aye*uaye), 3); // urho is [dens]/a^3
  double kunit   = (pow(uaye, 3) * mh) / (dbase1 * tbase1);

//================================================================================================//

  double Tev = temp * 8.621738e-5;
  //double XRTev = XRT * 8.621738e-5;
  double logTev = log(Tev);
  double tinv = 1.0 / temp;
  double kB = 1.381e-16; // erg/K
  double pi = 3.1415;
  double G = 6.67e-8; // cm^3/s^2/g
  double mH = 1.67262171e-24;
  double Z_solar = 0.01295;
  double dom = urho*pow(aye, 3)/mh;
  /* Note: The forula for Jeans length is sqrt((15*kB*T)/(4*pi*G*mu*rho)) */
  /*       obviously mu is the mean mass per particle and rho is the density*/
  double mu = 2*mH; // mass of H2 is twice the proton mass, assuming all are neutral

  // n gets passed in not as a number density, but as a mass density in code units
  // We need to convert to physical units and then to number density  
  n = n*urho/mu;

  double tov300 = temp / 300.;
  double T2 = temp / 100.;

  double rho_H2 = n*mu; //note that we're assuming nearly everything here is molecular hydrogen (n is actually total number                         density)
  double jeans = sqrt(15.0*kB* temp/(4.0*pi*G*mu*rho_H2)); //jeans length
  double Nh2 = n*jeans; // estimating H2 column density
  //double xx = Nh2/(5.0e14); // dimensionless
  //double b5 = sqrt( 2.0 * kB*temp/mH)/100000.0; //doppler parameter, from cm/s to km/s
  // Self-shielding of molecular Hydrogen
  //double shield = 0.965/pow(1.0 + xx/b5, 1.1) + 0.035/pow(1.0+xx, 0.5) * exp(-pow(1180.0,-1)*pow(1.0 + xx, 0.5)); // eqn (26) in the textbook given below

   //these coefficients are taken from Bialy & Sternberg 2015
   
   // if shielding is high, assume we use the LW-blocked rates instead of the thin rates
// If column density of H2 is greater than 1.e22 cm^-2, then we'll have the LW-blocked rates
   
   rate[UV1] = 0.0;
   rate[UV2] = 0.0;
   rate[UV3] = 0.0;
   rate[UV4] = 0.0;
   rate[UV5] = 0.0;
   rate[UV6] = 0.0;
   rate[UV7] = 0.0;
   rate[UV8] = 0.0;

   /* H1 - H6 already treated in coll_rates_g.F and calc_rates_g.F, so are zeroed here. */
   rate[H1] = 0.;
   rate[H2] = 0.;
   rate[H3] = 0.;
   rate[H4] = 0.;
   rate[H5] = 0.;
   rate[H6] = 0.;
   rate[H7] = 0.;
   rate[H8] = 0.;

   /* H8, D1 - D6 already treated in coll_rates_g.F and calc_rates_g.F, so are zeroed here. */
   //Comment out for now! Maybe add back in later? 
   rate[D1] = 0.;
   rate[D2] = 0.;
   rate[D3] = 0.;
   rate[D4] = 0.;
   rate[D5] = 0.;
   rate[D6] = 0.;

   rate[Z1] = 0.;
   rate[Z2] = 0.;
   rate[Z3] = 0.;
   rate[Z4] = 0.;
   rate[Z5] = 0.;
   rate[Z6] = 0.;
   rate[Z7] = 0.;
   rate[Z8] = 0.;
   rate[Z9] = 0.;
   rate[Z10] = 0.;
   rate[Z11] = 0.;
   rate[Z12] = 0.;
   rate[Z13] = 0.;
   rate[Z14] = 0.;
   rate[Z15] = 0.;
   rate[Z16] = 0.;
   rate[Z17] = 0.;
   rate[Z18] = 0.;
   rate[Z19] = 0.;
   rate[Z20] = 0.;
   rate[Z21] = 0.;
   rate[Z22] = 0.;
   rate[Z23] = 0.;
   rate[Z24] = 0.;
   rate[Z5] = 0.;
   rate[Z26] = 0.;
   rate[Z27] = 0.;
   rate[Z28] = 0.;
   rate[Z29] = 0.;
   rate[Z30] = 0.;
   rate[Z31] = 0.;
   rate[Z32] = 0.;
   rate[Z33] = 0.;
   rate[Z34] = 0.;
   rate[Z35] = 0.;
   rate[Z36] = 0.;
   rate[Z37] = 0.;
   rate[Z38] = 0.;
   rate[Z39] = 0.;
   rate[Z40] = 0.;

  /* Hydrogen formation reactions */ 

  /* NOTE: H1-H6 rates are already calculated in calc_rates_g.F
           with updated values, so we provide the original 
           Omukai (2005) rates here only for the purposes of the 
           One-zone freefall test.
   */
  if (water_rates == 1)
  {
        //Deuterium reactions
    rate[D1]= 3.7e-10 * pow( temp, 0.28) * exp( -43.0 * tinv);
    rate[D2]= 3.7e-10 * pow( temp, 0.28);
    rate[D3]= 9.0e-11 * exp( -3876.0 * tinv);
    rate[D4]= 2.1e-9;
    rate[D5] = 3.2e-11 * exp( -3624.0 * tinv);
    rate[D6]= 1.0e-9 * exp( -464.0 * tinv); 

  } else  {
        if (H2_shield > 0){
        // Compute shielding factor for H
        double avgsighi = crsHI;
        double nSSh = 6.73e-3 *  pow(avgsighi/2.49e-18, -2./3.) *
        pow(temp/1.0e4, 0.17) * pow(k24/tbase1/1.0e-12, 2.0/3.0);

        // Compute the total Hydrogen number density
        double nratio = (Y[H] + Y[Hplus]);
        if (ispecies > 1){
          nratio = nratio + Y[Hmin] + Y[H2m];
          if (ispecies > 2){
             nratio = nratio + 0.5*(Y[D] + Y[Dplus]) + 2.0*Y[HD]/3.0;
          }
        }

        double nH = nratio*urho/mu;

        nratio = nratio*dom/nSSh;

/*      LW-blocked Draine */
        rate[UV1] = 2.8e-10*UV/n;
        rate[UV2] = 0.28e-10*UV/n;
        rate[UV3] = 5.5e-10*UV/n;
        rate[UV4] = 7.0e-10*UV/n;
        rate[UV5] = 8.8e-10*UV/n;
        rate[UV6] = 0.0e-10*UV/n;
        rate[UV7] = 0.0e-10*UV/n;
        rate[UV8] = 0.0e-12*UV/n; 
     }
   }
   if (water_rates == 3)
   {
       /****** BIALY RATES **********/
      //rate[H2] = 4.0e-16 * pow(tov300, 0.67) * exp(4.0/temp);
      //rate[H3] = 1.30e-9;
      //rate[H7] = 1.0e-8 * exp(-84100.0/temp);
      //rate[H8] = 3e-17 * pow(T2, 0.5) * pow(Z, 1.0) ;
      
      //Metallicity Reactions
      rate[Z1]  = 5.66e-10 * pow(tov300, 0.36) * exp(8.6 * tinv);
      rate[Z2]  = 7.31e-10 * pow(tov300, 0.23) * exp(-225.9 * tinv);
      rate[Z3]  = 1.70e-9;
      rate[Z4]  = 1.01e-9;
      rate[Z5]  = 6.40e-10;
      rate[Z6]  = 3.60e-8   * pow ( tov300, -0.50 );
      rate[Z7]  = 7.92e-8   * pow ( tov300, -0.50 );
      rate[Z8]  = 2.58e-7   * pow ( tov300, -0.50 );
      rate[Z9]  = 1.08e-7   * pow ( tov300, -0.50 );
      rate[Z10] = 9.90e-19 * pow ( tov300, -0.38 );
      rate[Z11] = 3.43e-13 * pow ( tov300, 2.70 )   * exp ( -3150.0 * tinv );
      rate[Z12] = 2.05e-12 * pow ( tov300, 1.52  )   * exp ( -1736.0 * tinv );
      rate[Z13] = 1.65e-12 * pow ( tov300, 1.14  )   * exp ( -50.0   * tinv );
      rate[Z14] = 2.10e-9;
      rate[Z15] = 6.90e-9;
      rate[Z16] = 6.99e-14 * pow ( tov300, 2.80  )   * exp ( -1950.0  * tinv );
      rate[Z17] = 1.59e-11 * pow ( tov300, 1.20  )   * exp ( -9610.0 * tinv );
      rate[Z18] = 4.90e-20 * pow ( tov300, 1.58  );
      rate[Z19] = 4.34e-11 * pow ( tov300, -0.50 )   * exp ( -30.0    * tinv );
      rate[Z20] = 2.61e-10 * exp ( -8156.0 * tinv );
      rate[Z21] = 2.00e-9;
      rate[Z22] = 1.95e-7  * pow(  tov300, -0.70);
      rate[Z23] = 6.60e-11;
      rate[Z24] = 1.10e-10 * pow(  tov300, 0.50);
      rate[Z25] = 3.30e-11;
      rate[Z26] = 6.20e-10;
      rate[Z27] = 1.17e-13 * pow( tov300, 0.95) * exp(74.0 * tinv);
      rate[Z28] = 4.67e-12 * pow(  tov300, -0.60);
      rate[Z29] = 7.70e-10;
      rate[Z30] = 7.50e-10;
      rate[Z31] = 1.00e-17;
      rate[Z32] = 6.64e-10 * exp( -11700.0 * tinv);
      rate[Z33] = 2.70e-11 * pow( tov300, 0.38);
      rate[Z34] = 5.46e-10 * exp( -1943.0 * tinv);
      rate[Z35] = 6.64e-11;
      rate[Z36] = 5.18e-11 * pow(  tov300, 0.17) * exp( -6400.0 * tinv);
      rate[Z37] = 1.00e-10* exp( -7600.0 * tinv);
      rate[Z38] = 6.86e-14 * pow(  tov300, 2.74) * exp( -4740.0 * tinv);
      rate[Z39] = 5.94e-13 * pow(  tov300, 3.00) * exp( -4045.0 * tinv);
      rate[Z40] = 1.00e-17;
    } else{

      //printf("temp = %.2f\n", temp);
      //H7 is currently in neither grackle nor Bialy, so options 2-3 all need to use Omukai's rate
      double klow = 1.18e-10 * exp(-6.95e4*tinv);
      double khigh = 8.125e-8 * pow( temp, -0.5) * exp(-5.2e4*tinv) * ( (double) 1.0 - exp(-6.0e3 *tinv));
      double ep2 = 4.845 - 1.3* log10( temp * 1.0e-4) + 1.62* pow( log10( temp * 1.0e-4) , 2.0);
      double ncrit = pow(10, ep2);
      double a = (double) 1.0 / ( (double) 1.0 + n / ncrit);

      //rate[H7]= pow( khigh, 1.0 - a) * pow( klow, a);

      //Metallicity Reactions
      rate[Z1]  = 6.80e-10;
      rate[Z2]  = 7.00e-10 * exp( -232.0 * tinv);
      rate[Z3]  = 1.70e-9;
      rate[Z4]  = 1.01e-9;
      rate[Z5]  = 8.30e-10;
      rate[Z6]  = 2.00e-7   * pow ( tov300, -0.50 );
      rate[Z7]  = 1.60e-7   * pow ( tov300, -0.50 );
      rate[Z8]  = 6.50e-7   * pow ( tov300, -0.50 );
      rate[Z9]  = 3.50e-7   * pow ( tov300, -0.50 );
      rate[Z10] = 9.90e-19 * pow ( tov300, -0.38 );
      rate[Z11] = 3.43e-13 * pow ( tov300, 2.67  )   * exp ( -3160.0 * tinv );
      rate[Z12] = 1.55e-12 * pow ( tov300, 1.60  )   * exp ( -1660.0 * tinv );
      rate[Z13] = 1.65e-12 * pow ( tov300, 1.14  )   * exp ( -50.0   * tinv );
      rate[Z14] = 2.10e-9;
      rate[Z15] = 6.90e-9;
      rate[Z16] = 7.00e-14 * pow ( tov300, 2.80  )   * exp ( -1950.0  * tinv );
      rate[Z17] = 6.83e-12 * pow ( tov300, 1.60  )   * exp ( -9720.0 * tinv );
      rate[Z18] = 4.90e-20 * pow ( tov300, 1.58  );
      rate[Z19] = 4.34e-11 * pow ( tov300, -0.50 )   * exp ( -30.0    * tinv );
      rate[Z20] = 3.30e-10 * exp ( -8460.0 * tinv );
      rate[Z21] = 2.00e-9;
      rate[Z22] = 1.95e-7  * pow(  tov300, -0.70);
      rate[Z23] = 6.60e-11;
      rate[Z24] = 1.10e-10 * pow(  tov300, 0.50);
      rate[Z25] = 3.30e-11;
      rate[Z26] = 6.20e-10;
      rate[Z27] = 1.00e-13;
      rate[Z28] = 4.40e-12 * pow(  tov300, -0.61);
      rate[Z29] = 7.70e-10;
      rate[Z30] = 7.50e-10;
      rate[Z31] = 1.00e-17;
      rate[Z32] = 6.64e-10 * exp( -11700.0 * tinv);
      rate[Z33] = 4.98e-11;
      rate[Z34] = 2.38e-10 * exp( -1760.0 * tinv);
      rate[Z35] = 2.70e-10;
      rate[Z36] = 5.18e-11 * pow(  tov300, 0.17) * exp( -6400.0 * tinv);
      rate[Z37] = 1.00e-10* exp( -7600.0 * tinv);
      rate[Z38] = 6.86e-14 * pow(  tov300, 2.74) * exp( -4740.0 * tinv);
      rate[Z39] = 5.82e-13 * pow(  tov300, 3.00) * exp( -4045.0 * tinv);
      rate[Z40] = 1.00e-17;
    }
  }
