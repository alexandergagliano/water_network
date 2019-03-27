#include <stdlib.h>
#include <stdio.h>
#include "grackle_chemistry_data.h"

/****************** Declaration of all chemical species ****************************/
// Multispecies = 1, withWater = 1
int H, Hplus, el, O, OH, H2O, O2, Oplus, OHplus, H2Oplus, H3Oplus, O2plus, Cplus, C, CH, CH2, CH3, CH4, CO, COplus, CO2; 

// Multispecies = 2, withWater = 1
int H2m, Hmin; 

// Multispecies = 3, withWater = 1
int D, Dplus, HD;

// withWater = 1, water_rates = 3 (Using Bialy's expanded network)
int SI, SII, SiI, SiII, CHII, CH2II, HeI, HeII, SiO, H3II;

/***************** Declaration of all chemical reactions ***************************/  
int H1, H2, H3, H4, H5, H6, H7, H8;
int D1, D2, D3, D4, D5, D6;
int Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z12, Z13, Z14, Z15, Z16, 
    Z17, Z18, Z19, Z20, Z21, Z22, Z23, Z24, Z25, Z26, Z27, Z28, Z29, Z30, 
    Z31, Z32, Z33, Z34, Z35, Z36, Z37, Z38, Z39, Z40, Z41, Z42, Z43, Z44, 
    Z45, Z46, Z47, Z48, Z49, Z50, Z51, Z52, Z53;
int UV1, UV2, UV3, UV4, UV5, UV6, UV7, UV8;

/*************** Initializing variables for reactions and species *****************/ 
void setup_rxns(int ispecies, int UV, int water_rates)
{
   int rxn_counter = 0;
   H1 = rxn_counter;

   Z1 = ++rxn_counter;
   Z2 = ++rxn_counter;
   Z7 = ++rxn_counter;
   Z8 = ++rxn_counter; 
   Z9 = ++rxn_counter;
   Z10 = ++rxn_counter;
   Z13 = ++rxn_counter;
   Z14 = ++rxn_counter;
   Z15 = ++rxn_counter;
   Z18 = ++rxn_counter;
   Z19 = ++rxn_counter;
   Z20 = ++rxn_counter;
   Z21 = ++rxn_counter;
   Z22 = ++rxn_counter;
   Z23 = ++rxn_counter;
   Z24 = ++rxn_counter;
   Z25 = ++rxn_counter;
   Z26 = ++rxn_counter;
   Z27 = ++rxn_counter;
   Z28 = ++rxn_counter;
   Z29 = ++rxn_counter;
   Z30 = ++rxn_counter;
   Z31 = ++rxn_counter;
   
   if (ispecies > 1){
      H2 = ++rxn_counter;
      H3 = ++rxn_counter;
      H4 = ++rxn_counter;
      H5 = ++rxn_counter;
      H6 = ++rxn_counter;
      H7 = ++rxn_counter;
      H8 = ++rxn_counter;
      Z3 = ++rxn_counter;
      Z4 = ++rxn_counter;
      Z5 = ++rxn_counter;
      Z6 = ++rxn_counter;
      Z11 = ++rxn_counter;
      Z12 = ++rxn_counter;
      Z16 = ++rxn_counter;
      Z17 = ++rxn_counter;
      Z32 = ++rxn_counter;
      Z33 = ++rxn_counter;
      Z34 = ++rxn_counter;
      Z35 = ++rxn_counter;
      Z36 = ++rxn_counter;
      Z37 = ++rxn_counter;
      Z38 = ++rxn_counter;
      Z39 = ++rxn_counter;
      Z40 = ++rxn_counter;

      if (ispecies > 2){
         D1  = ++rxn_counter;
         D2  = ++rxn_counter;
         D3  = ++rxn_counter;
         D4  = ++rxn_counter;
         D5  = ++rxn_counter;
         D6  = ++rxn_counter;
      }
   }
     
   if (UV) {
      UV1 = ++rxn_counter;
      UV2 = ++rxn_counter;
      UV3 = ++rxn_counter;
      UV4 = ++rxn_counter;
      UV5 = ++rxn_counter;
      UV6 = ++rxn_counter;
      UV7 = ++rxn_counter;
      UV8 = ++rxn_counter;
   }
    
   if (water_rates == 3){
      Z41 = ++rxn_counter;
      Z42 = ++rxn_counter;
      Z43 = ++rxn_counter;
      Z44 = ++rxn_counter;
      Z45 = ++rxn_counter;
      Z46 = ++rxn_counter;
      Z47 = ++rxn_counter;
      Z48 = ++rxn_counter;
      Z49 = ++rxn_counter;
      Z50 = ++rxn_counter;
      Z51 = ++rxn_counter;
      Z52 = ++rxn_counter;
      Z53 = ++rxn_counter;
   }
}

void setup_species(int ispecies, int withWater, int UV, int water_rates)
{
  int sp_counter = 0;
  H       = sp_counter; 

  Hplus   = ++sp_counter; 
  el      = ++sp_counter; 
  O       = ++sp_counter; 
  OH      = ++sp_counter; 
  H2O     = ++sp_counter; 
  O2      = ++sp_counter; 
  Oplus   = ++sp_counter;
  OHplus  = ++sp_counter; 
  H2Oplus = ++sp_counter; 
  H3Oplus = ++sp_counter; 
  O2plus  = ++sp_counter; 
  Cplus   = ++sp_counter; 
  C       = ++sp_counter; 
  CH      = ++sp_counter;  
  CH2     = ++sp_counter;
  CH3     = ++sp_counter; 
  CH4     = ++sp_counter; 
  CO      = ++sp_counter; 
  COplus  = ++sp_counter; 
  CO2     = ++sp_counter;

  if (ispecies > 1){
     H2m  = ++sp_counter;  
     Hmin = ++sp_counter;

     if (ispecies > 2){ 
        D     = ++sp_counter; 
        Dplus = ++sp_counter; 
        HD    = ++sp_counter;
     }
  }

/*    
  if (water_rates == 3){
     SI    = ++sp_counter;
     SII   = ++sp_counter;
     SiI   = ++sp_counter;
     SiII  = ++sp_counter;
     CHII  = ++sp_counter;
     CH2II = ++sp_counter;
     HeI   = ++sp_counter;
     HeII  = ++sp_counter;
     SiO   = ++sp_counter;
     H3II  = ++sp_counter;
  } 
*/
}
