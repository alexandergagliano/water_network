#include <stdlib.h>
#include <stdio.h>
#include "grackle_chemistry_data.h"

// Multispecies = 1, withWater = 1
int H       = 0  ; int Hplus   = 1  ; int el      = 2  ; int O       = 3  ; int OH      = 4  ; int H2O     = 5  ; int O2      = 6  ; int Oplus   = 7  ;
int OHplus  = 8  ; int H2Oplus = 9  ; int H3Oplus = 10 ; int O2plus  = 11 ; int Cplus   = 12 ; int C       = 13 ; int CH      = 14 ; int CH2     = 15 ;
int CH3     = 16 ; int CH4     = 17 ; int CO      = 18 ; int COplus  = 19 ; int CO2     = 20 ; 

// Multispecies = 2, withWater = 1
int H2m  = 21 ; int Hmin    = 22;
// Multispecies = 3, withWater = 1
int D    = 23 ; int Dplus   = 24 ; int HD     = 25; 

// Reactions for Multispecies = 1, withWater = 1
int H1  = 0  ; 

int Z1  = 1  ; int Z2  = 2  ; int Z7  = 3  ; int Z8  = 4  ; int Z9  = 5  ; int Z10 = 6 ; 
int Z13 = 7  ; int Z14 = 8  ; int Z15 = 9  ; int Z18 = 10 ; int Z19 = 11 ; int Z20 = 12 ; 
int Z21 = 13 ; int Z22 = 14 ; int Z23 = 15 ; int Z24 = 16 ; int Z25 = 17 ; int Z26 = 18 ; 
int Z27 = 19 ; int Z28 = 20 ; int Z29 = 21 ; int Z30 = 22 ; int Z31 = 23 ; 

// Reactions for Multispecies = 2, withWater = 1
int H2  = 24  ; int H3  = 25  ; int H4  = 26  ; int H5  = 27  ; int H6  = 28  ; int H7  = 29  ; 
int H8  = 30  ; 

int Z3  = 31  ; int Z4  = 32  ; int Z5  = 33  ; int Z6  = 34  ; int Z11 = 35  ;
int Z12 = 36  ; int Z16 = 37  ; int Z17 = 38  ; int Z32 = 39  ; int Z33 = 40  ; int Z34 = 41  ;
int Z35 = 42  ; int Z36 = 43  ; int Z37 = 44  ; int Z38 = 45  ; int Z39 = 46  ; int Z40 = 47  ; 

// Reactions for Multispecies = 3, withWater = 1
int D1  = 48  ; int D2  = 49  ; int D3  = 50  ; int D4  = 51  ; int D5  = 52  ; int D6  = 53  ;

// Reactions for UV = 1, withWater = 1, Multispecies = 1,2,3
int UV1 = 54 ; int UV2 = 55 ; int UV3 = 56 ;  int UV4 = 57 ;  int UV5 = 58 ;  int UV6 = 59 ;  int UV7 = 60 ;  int UV8 = 61 ; 
