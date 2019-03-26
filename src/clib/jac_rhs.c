#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <jac_rhs.h>
#include <chemistry.h>
#include <gsl/gsl_math.h>

// clip magnitude of jacobian elements to 1e-99 to avoid underflow error
void clip_jacobian(double **J, int ispecies){
  int i, j;
  double jabs, sign;
  int nSpecies = (ispecies > 1) ? ((ispecies > 2) ? 26 : 23) : 21;
  for (i = 0; i < nSpecies; i++) {
    for (j = 0; j < nSpecies; j++) {
      sign = (double) GSL_SIGN(J[i][j]);
      jabs = fabs(J[i][j]);
      if (jabs > 0.) {
        J[i][j] = sign * fmax(jabs,1e-99);
      }
      else {
        J[i][j] = 0.;
      }
    }
  }
}

void get_f_vector(double *rate, double *F, double *Y, int ispecies) {

/**************************************************/
/*              MULTISPECIES = 1                  */
/**************************************************/

/************* F[H] ******************/

  F[H] =  rate[H1]*Y[Hplus]*Y[el] 
    +     rate[Z2]*Y[Hplus]*Y[O] 
    +     rate[Z7]*Y[H2Oplus]*Y[el]
    +     rate[Z9]*Y[H3Oplus]*Y[el] 
    +     rate[Z14]*Y[Hplus]*Y[OH] 
    +     rate[Z15]*Y[Hplus]*Y[H2O] 
    +     rate[Z19]*Y[O]*Y[OH]
    +     rate[Z21]*Y[Hplus]*Y[O2] 
    +     rate[Z23]*Y[O]*Y[CH] 
    +     rate[Z24]*Y[C]*Y[OH]
    +     rate[Z27]*Y[OH]*Y[CO] 
    +     rate[Z29]*Y[Cplus]*Y[OH] 
    + 2.0*rate[Z8]*Y[H3Oplus]*Y[el] 
    +     rate[UV3]*Y[H2O] 
    -     rate[Z1]*Y[Oplus]*Y[H] 
    -     rate[Z10]*Y[O]*Y[H] 
    -     rate[Z20]*Y[O2]*Y[H] 
    -     rate[Z30]*Y[COplus]*Y[H] 
    -     rate[Z31]*Y[C]*Y[H]
    +     rate[UV5]*Y[CH]
    +     rate[UV1]*Y[OH];

/************* F[el] ******************/

  F[el] = 
    -     rate[H1]*Y[Hplus]*Y[el] 
    -     rate[Z7]*Y[H2Oplus]*Y[el] 
    -     rate[Z8]*Y[H3Oplus]*Y[el] 
    -     rate[Z9]*Y[H3Oplus]*Y[el]
    -     rate[Z22]*Y[O2plus]*Y[el] 
    -     rate[Z28]*Y[Cplus]*Y[el]
    +     rate[UV7]*Y[C]
    +     rate[UV8]*Y[O];

/************* F[Hplus] ******************/

  F[Hplus] = 
    +        rate[Z1]*Y[Oplus]*Y[H] 
    +        rate[Z30]*Y[COplus]*Y[H]
    -        rate[Z2]*Y[Hplus]*Y[O]
    -        rate[Z14]*Y[Hplus]*Y[OH] 
    -        rate[Z15]*Y[Hplus]*Y[H2O] 
    -        rate[Z21]*Y[Hplus]*Y[O2]
    -        rate[H1]*Y[Hplus]*Y[el]; 

/************* F[O] ******************/

  F[O] =  rate[Z1]*Y[Oplus]*Y[H] 
    +     rate[Z13]*Y[OH]*Y[OH] 
    +     rate[Z20]*Y[H]*Y[O2] 
    + 2.0*rate[Z22]*Y[O2plus]*Y[el] 
    +     rate[Z25]*Y[C]*Y[O2] 
    +     rate[UV6]*Y[CO] 
    +     rate[UV2]*Y[H2O]
    +     rate[UV1]*Y[OH];
    -     rate[Z2]*Y[Hplus]*Y[O] 
    -     rate[Z10]*Y[O]*Y[H] 
    - 2.0*rate[Z18]*Y[O]*Y[O] 
    -     rate[Z19]*Y[O]*Y[OH] 
    -     rate[Z23]*Y[O]*Y[CH]
    + 2.0*rate[UV4]*Y[O2]
    -     rate[UV8]*Y[O];

/************* F[OH] ******************/

  F[OH] = rate[Z7]*Y[H2Oplus]*Y[el] 
    +     rate[Z8]*Y[H3Oplus]*Y[el] 
    +     rate[Z10]*Y[O]*Y[H]
    +     rate[Z20]*Y[H]*Y[O2] 
    +     rate[UV3]*Y[H2O]
    - 2.0*rate[Z13]*Y[OH]*Y[OH] 
    -     rate[Z14]*Y[Hplus]*Y[OH]
    -     rate[Z19]*Y[O]*Y[OH] 
    -     rate[Z24]*Y[C]*Y[OH]
    -     rate[Z27]*Y[OH]*Y[CO] 
    -     rate[Z29]*Y[Cplus]*Y[OH]
    -     rate[UV1]*Y[OH];

/************* F[H2O] ******************/

  F[H2O] = rate[Z9]*Y[H3Oplus]*Y[el] 
    +      rate[Z13]*Y[OH]*Y[OH]  
    -      rate[Z15]*Y[Hplus]*Y[H2O]
    -      rate[UV3]*Y[H2O] 
    -      rate[UV2]*Y[H2O];

/************* F[O2] ******************/

  F[O2] = rate[Z18]*Y[O]*Y[O] 
    +     rate[Z19]*Y[O]*Y[OH]
    -     rate[Z20]*Y[H]*Y[O2] 
    -     rate[Z21]*Y[Hplus]*Y[O2] 
    -     rate[Z25]*Y[C]*Y[O2]
    -     rate[Z26]*Y[Cplus]*Y[O2]
    -     rate[UV4]*Y[O2];

/************* F[Oplus] ******************/

  F[Oplus] = rate[Z2]*Y[Hplus]*Y[O] 
    +        rate[Z26]*Y[Cplus]*Y[O2]
    -        rate[Z1]*Y[Oplus]*Y[H]
    +        rate[UV8]*Y[O];

/************* F[OHplus] ******************/

  F[OHplus] = rate[Z14]*Y[Hplus]*Y[OH];

/************* F[H2Oplus] ******************/

  F[H2Oplus] = rate[Z15]*Y[Hplus]*Y[H2O]
    -          rate[Z7]*Y[H2Oplus]*Y[el]; 

/************* F[H3Oplus] ******************/

  F[H3Oplus] = -rate[Z8]*Y[H3Oplus]*Y[el] 
    -           rate[Z9]*Y[H3Oplus]*Y[el]; 

/************* F[O2plus] ******************/

  F[O2plus] = rate[Z21]*Y[Hplus]*Y[O2] 
    -         rate[Z22]*Y[O2plus]*Y[el]; 

/************* F[Cplus] ******************/

  F[Cplus] = -rate[Z26]*Y[Cplus]*Y[O2] 
    -         rate[Z28]*Y[Cplus]*Y[el] 
    -         rate[Z29]*Y[Cplus]*Y[OH]
    +         rate[UV7]*Y[C];

/************* F[C] ******************/

  F[C] = rate[Z28]*Y[Cplus]*Y[el] 
    +    rate[UV6]*Y[CO]
    -    rate[Z24]*Y[C]*Y[OH] 
    -    rate[Z25]*Y[C]*Y[O2] 
    -    rate[Z31]*Y[C]*Y[H]
    +    rate[UV5]*Y[CH]
    -    rate[UV7]*Y[C];

/************* F[CH] ******************/

  F[CH] = rate[Z31]*Y[C]*Y[H]
   -      rate[Z23]*Y[O]*Y[CH]
   -      rate[UV5]*Y[CH];

/************* F[CO] ******************/

  F[CO] = rate[Z23]*Y[O]*Y[CH] 
   +      rate[Z24]*Y[C]*Y[OH] 
   +      rate[Z25]*Y[C]*Y[O2] 
   +      rate[Z26]*Y[Cplus]*Y[O2]
   +      rate[Z30]*Y[COplus]*Y[H]
   -      rate[Z27]*Y[OH]*Y[CO];
   -      rate[UV6]*Y[CO]; 

/************* F[COplus] ******************/

  F[COplus] = rate[Z29]*Y[Cplus]*Y[OH] 
            - rate[Z30]*Y[COplus]*Y[H]; 

/************* F[CO2] ******************/

  F[CO2] = rate[Z27]*Y[OH]*Y[CO]; 


/**************************************************/
/*              MULTISPECIES = 2                  */
/**************************************************/
  if (ispecies > 1)
  {
    F[H] += 3.0*rate[H4]*Y[H2m]*Y[H]
    +           rate[H5]*Y[H]*Y[H]*Y[H]
    +       2.0*rate[H7]*Y[H2m]*Y[H2m]
    +           rate[Z3]*Y[Oplus]*Y[H2m]
    +           rate[Z4]*Y[OHplus]*Y[H2m]
    +           rate[Z5]*Y[H2Oplus]*Y[H2m]
    +           rate[Z11]*Y[O]*Y[H2m]
    +           rate[Z12]*Y[H2m]*Y[OH]
    +           rate[Z32]*Y[C]*Y[H2m]
    +           rate[Z34]*Y[H2m]*Y[CH]
    +           rate[Z36]*Y[H2m]*Y[CH2]
    +           rate[Z38]*Y[H2m]*Y[CH3]
    -           rate[H2]*Y[H]*Y[el]
    -           rate[H3]*Y[Hmin]*Y[H]
    -           rate[H4]*Y[H2m]*Y[H]
    -       3.0*rate[H5]*Y[H]*Y[H]*Y[H]
    -       2.0*rate[H6]*Y[H]*Y[H]*Y[H2m]
    -       2.0*rate[H8]*Y[H]*Y[H]
    -           rate[Z16]*Y[OH]*Y[H]
    -           rate[Z17]*Y[H2O]*Y[H]
    -           rate[Z33]*Y[CH]*Y[H]
    -           rate[Z35]*Y[CH2]*Y[H]
    -           rate[Z37]*Y[CH3]*Y[H]
    -           rate[Z39]*Y[CH4]*Y[H];

    F[el] += rate[H3]*Y[Hmin]*Y[H]
    -        rate[H2]*Y[H]*Y[el]
    -     rate[Z6]*Y[H2Oplus]*Y[el];

    F[O] += rate[Z6]*Y[H2Oplus]*Y[el]
    +       rate[Z16]*Y[H]*Y[OH]
    -       rate[Z11]*Y[O]*Y[H2m];

    F[OH] += rate[Z11]*Y[O]*Y[H2m]
      +      rate[Z17]*Y[H]*Y[H2O]
      -      rate[Z12]*Y[H2m]*Y[OH]
      -      rate[Z16]*Y[H]*Y[OH];

    F[H2O] += rate[Z12]*Y[OH]*Y[H2m]
      -       rate[Z17]*Y[H]*Y[H2O];

    F[Oplus] += -rate[Z3]*Y[Oplus]*Y[H2m];

    F[OHplus] += rate[Z3]*Y[Oplus]*Y[H2m]
      -          rate[Z4]*Y[OHplus]*Y[H2m];

    F[H2Oplus] += rate[Z4]*Y[OHplus]*Y[H2m]
    -             rate[Z5]*Y[H2Oplus]*Y[H2m]
    -             rate[Z6]*Y[H2Oplus]*Y[el];

    F[H3Oplus] += rate[Z5]*Y[H2Oplus]*Y[H2m];

    F[C] += rate[Z33]*Y[H]*Y[CH]
      -     rate[Z32]*Y[C]*Y[H2m]
      -     rate[Z40]*Y[H2m]*Y[C];

    F[CH] += rate[Z32]*Y[C]*Y[H2m]
     +       rate[Z35]*Y[H]*Y[CH2]
     -       rate[Z33]*Y[H]*Y[CH]
     -       rate[Z34]*Y[H2m]*Y[CH];

    F[CH2] = rate[Z34]*Y[H2m]*Y[CH]
     +       rate[Z37]*Y[H]*Y[CH3]
     +       rate[Z40]*Y[H2m]*Y[C]
     -       rate[Z35]*Y[H]*Y[CH2]
     -       rate[Z36]*Y[H2m]*Y[CH2];

    F[CH3] = rate[Z36]*Y[H2m]*Y[CH2]
     +       rate[Z39]*Y[H]*Y[CH4]
     -       rate[Z38]*Y[H2m]*Y[CH3]
     -       rate[Z37]*Y[H]*Y[CH3];

    F[CH4] = rate[Z38]*Y[H2m]*Y[CH3]
     -       rate[Z39]*Y[H]*Y[CH4];

      F[H2m] = rate[H3]*Y[Hmin]*Y[H] 
        +  2.0*rate[H6]*Y[H]*Y[H]*Y[H2m] 
        +      rate[H7]*Y[H2m]*Y[H2m] 
        +      rate[H8]*Y[H]*Y[H] 
        +      rate[H5]*Y[H]*Y[H]*Y[H]
        +      rate[Z6]*Y[H2Oplus]*Y[el] 
        +      rate[Z16]*Y[H]*Y[OH] 
        +      rate[Z17]*Y[H]*Y[H2O]
        +      rate[Z33]*Y[H]*Y[CH] 
        +      rate[Z35]*Y[H]*Y[CH2] 
        +      rate[Z37]*Y[H]*Y[CH3] 
        +      rate[Z39]*Y[CH4]*Y[H] 
        +      rate[UV2]*Y[H2O]
        -      rate[H4]*Y[H2m]*Y[H] 
        -      rate[H6]*Y[H]*Y[H]*Y[H2m]  
        -  2.0*rate[H7]*Y[H2m]*Y[H2m] 
        -      rate[Z3]*Y[Oplus]*Y[H2m] 
        -      rate[Z4]*Y[OHplus]*Y[H2m] 
        -      rate[Z5]*Y[H2Oplus]*Y[H2m]
        -      rate[Z11]*Y[O]*Y[H2m] 
        -      rate[Z12]*Y[H2m]*Y[OH] 
        -      rate[Z32]*Y[C]*Y[H2m] 
        -      rate[Z34]*Y[H2m]*Y[CH]
        -      rate[Z36]*Y[H2m]*Y[CH2] 
        -      rate[Z38]*Y[H2m]*Y[CH3] 
        -      rate[Z40]*Y[H2m]*Y[C];

      F[Hmin] = rate[H2]*Y[H]*Y[el] 
              - rate[H3]*Y[Hmin]*Y[H];

/**************************************************/
/*              MULTISPECIES = 2                  */
/**************************************************/
      if (ispecies > 2){
          F[D] = rate[D2]*Y[Dplus]*Y[H] 
           +     rate[D5]*Y[HD]*Y[H]
           -     rate[D1]*Y[D]*Y[Hplus] 
           -     rate[D3]*Y[D]*Y[H2m];

          F[Dplus] = rate[D1]*Y[D]*Y[Hplus]  
           +         rate[D6]*Y[HD]*Y[Hplus]
           -         rate[D2]*Y[Dplus]*Y[H] 
           -         rate[D4]*Y[Dplus]*Y[H2m];

          F[HD] = rate[D3]*Y[D]*Y[H2m] 
           +      rate[D4]*Y[Dplus]*Y[H2m]
           -      rate[D5]*Y[HD]*Y[H] 
           -      rate[D6]*Y[HD]*Y[Hplus];

          F[H2m] += rate[D5]*Y[HD]*Y[H]
           +        rate[D6]*Y[HD]*Y[Hplus]
           -        rate[D3]*Y[D]*Y[H2m]
           -        rate[D4]*Y[Dplus]*Y[H2m];

          F[H] += rate[D1]*Y[D]*Y[Hplus]
           +      rate[D3]*Y[D]*Y[H2m]
           -      rate[D2]*Y[H]*Y[Dplus]
           -      rate[D5]*Y[HD]*Y[H];

          F[Hplus] += rate[D2]*Y[Dplus]*Y[H]
           +          rate[D4]*Y[Dplus]*Y[H2m]
           -          rate[D1]*Y[D]*Y[Hplus]
           -          rate[D6]*Y[HD]*Y[Hplus];
      }
  }

}

/**
 * 03/04/2018
 * Changes made by Alex Gagliano:
 * 1. in J[H2m][H], changed 2.0*rate[H6]*Y[H]*Y[H2m] to 4.0*rate[H6]*Y[H]*Y[H2m] and
 *    -rate[H6]*Y[H]*Y[H2m] to -2.0*rate[H6]*Y[H]*Y[H2m]
 * 2. in J[H2m][HD], changed rate[D5]*Y[H2m] to rate[D5]*Y[H]
 * 3. in J[Dplus][Dplus], changed -rate[D5]*Y[H2m] to -rate[D4]*Y[H2m]
**/
void get_jacobian(double rate[], double **J, double Y[], int ispecies) {
  // We are taking the derivative of the rate equations
  //   \dot{y} = f(y)

/**************************************************/
/*              MULTISPECIES = 2                  */
/**************************************************/

/*************** H *******************/

  J[H][H] =
    -           rate[Z1]*Y[Oplus]
    -           rate[Z10]*Y[O] 
    -           rate[Z20]*Y[O2]
    -           rate[Z30]*Y[COplus] 
    -           rate[Z31]*Y[C];

  J[H][Hplus] = rate[H1]*Y[el] 
    +           rate[Z2]*Y[O] 
    +           rate[Z14]*Y[OH]
    +           rate[Z15]*Y[H2O] 
    +           rate[Z21]*Y[O2] ;

  J[H][el] = rate[H1]*Y[Hplus] 
    +        rate[Z7]*Y[H2Oplus] 
    +        rate[Z9]*Y[H3Oplus]
    +    2.0*rate[Z8]*Y[H3Oplus];

  J[H][O] = rate[Z2]*Y[Hplus] 
    +       rate[Z19]*Y[OH] 
    +       rate[Z23]*Y[CH]
    -       rate[Z10]*Y[H]; 

  J[H][OH] = rate[Z14]*Y[Hplus] 
    +        rate[Z19]*Y[O] 
    +        rate[Z24]*Y[C]
    +        rate[Z27]*Y[CO] 
    +        rate[Z29]*Y[Cplus]
    +        rate[UV1];

  J[H][H2O] = rate[Z15]*Y[Hplus]
    +         rate[UV3];

  J[H][O2] =  rate[Z21]*Y[Hplus]
    -         rate[Z20]*Y[H];

  J[H][Oplus] = -rate[Z1]*Y[H];

  J[H][H2Oplus] = rate[Z7]*Y[el]; 

  J[H][H3Oplus] = rate[Z9]*Y[el] 
   +          2.0*rate[Z8]*Y[el]; 

  J[H][Cplus] = rate[Z29]*Y[OH]; 

  J[H][C] = rate[Z24]*Y[OH]
    -       rate[Z31]*Y[H]; 

  J[H][CH] = rate[Z23]*Y[O]
  +          rate[UV5];

  J[H][CO] = rate[Z27]*Y[OH]; 

  J[H][COplus] = -rate[Z30]*Y[H]; 

/************** el *******************/

  J[el][Hplus] = -rate[H1]*Y[el]; 

  J[el][el] = -rate[H1]*Y[Hplus] 
    -          rate[Z7]*Y[H2Oplus] 
    -          rate[Z8]*Y[H3Oplus]
    -          rate[Z9]*Y[H3Oplus]
    -          rate[Z22]*Y[O2plus]
    -          rate[Z28]*Y[Cplus]; 

  J[el][H2Oplus] = -rate[Z7]*Y[el]; 

  J[el][H3Oplus] = -rate[Z8]*Y[el] 
    -               rate[Z9]*Y[el]; 

  J[el][O2plus] = -rate[Z22]*Y[el]; 

  J[el][C]      =  rate[UV7];

  J[el][Cplus] = -rate[Z28]*Y[el]; 

  J[el][O]     =  rate[UV8];

/************** Hplus *******************/

  J[Hplus][H] = rate[Z1]*Y[Oplus] 
    +           rate[Z30]*Y[COplus]; 

  J[Hplus][Hplus] = -rate[Z2]*Y[O]
    -                rate[Z14]*Y[OH]
    -                rate[Z15]*Y[H2O] 
    -                rate[Z21]*Y[O2] 
    -                rate[H1]*Y[el]; 

  J[Hplus][Oplus] = rate[Z1]*Y[H]; 

  J[Hplus][COplus] = rate[Z30]*Y[H]; 

  J[Hplus][O] = -rate[Z2]*Y[Hplus]; 

  J[Hplus][OH] = -rate[Z14]*Y[Hplus]; 

  J[Hplus][H2O] = -rate[Z15]*Y[Hplus]; 

  J[Hplus][O2] = -rate[Z21]*Y[Hplus]; 

  J[Hplus][el] = -rate[H1]*Y[Hplus]; 
 
/*************** O *******************/

  J[O][H] = rate[Z1]*Y[Oplus] 
    +       rate[Z20]*Y[O2]
    -       rate[Z10]*Y[O]; 

  J[O][Hplus] = -rate[Z2]*Y[O]; 

  J[O][CO] = rate[UV6];

  J[O][Oplus] = rate[Z1]*Y[H]; 

  J[O][O2] = rate[Z20]*Y[H] 
    +        rate[Z25]*Y[C]
    +    2.0*rate[UV4]; 

  J[O][O2plus] = 2.0*rate[Z22]*Y[el]; 

  J[O][OH] = 2.0*rate[Z13]*Y[OH] 
    -            rate[Z19]*Y[O]
    +            rate[UV1];

  J[O][el] = 2.0*rate[Z22]*Y[O2plus]; 

  J[O][H2O] = rate[UV2];

  J[O][C] = rate[Z25]*Y[O2] ;

  J[O][CH] = -rate[Z23]*Y[O] ;

  J[O][O] = -rate[Z2]*Y[Hplus] 
    -        rate[Z10]*Y[H] 
    -    4.0*rate[Z18]*Y[O] 
    -        rate[Z19]*Y[OH] 
    -        rate[Z23]*Y[CH]
    -        rate[UV8];

/*************** OH *******************/

  J[OH][H] = rate[Z10]*Y[O] 
    +        rate[Z20]*Y[O2];

  J[OH][Hplus] = -rate[Z14]*Y[OH]; 

  J[OH][O2] = rate[Z20]*Y[H]; 

  J[OH][O] = rate[Z10]*Y[H] 
    -        rate[Z19]*Y[OH]; 

  J[OH][OH] =
    -      4.0*rate[Z13]*Y[OH] 
    -          rate[Z14]*Y[Hplus] 
    -          rate[Z19]*Y[O] 
    -          rate[Z24]*Y[C] 
    -          rate[Z27]*Y[CO] 
    -          rate[Z29]*Y[Cplus]
    -          rate[UV1];

  J[OH][H2O] = rate[UV3];

  J[OH][H2Oplus] = rate[Z7]*Y[el]; 

  J[OH][H3Oplus] = rate[Z8]*Y[el]; 

  J[OH][el] = rate[Z7]*Y[H2Oplus] 
   +          rate[Z8]*Y[H3Oplus]; 

  J[OH][C] = -rate[Z24]*Y[OH]; 

  J[OH][CO] = -rate[Z27]*Y[OH];

  J[OH][Cplus] = -rate[Z29]*Y[OH]; 

/************** H2O *******************/

  J[H2O][Hplus] = -rate[Z15]*Y[H2O]; 

  J[H2O][el] = rate[Z9]*Y[H3Oplus]; 

  J[H2O][OH] = 2.0*rate[Z13]*Y[OH];

  J[H2O][H2O] = -rate[Z15]*Y[Hplus]
   -             rate[UV3] 
   -             rate[UV2];

  J[H2O][H3Oplus] = rate[Z9]*Y[el]; 

/*************** O2 *******************/

  J[O2][H] = -rate[Z20]*Y[O2]; 

  J[O2][Hplus] = -rate[Z21]*Y[O2]; 

  J[O2][O] = 2.0*rate[Z18]*Y[O] 
    +            rate[Z19]*Y[OH]; 

  J[O2][OH] = rate[Z19]*Y[O]; 

  J[O2][O2] = -rate[Z20]*Y[H] 
    -          rate[Z21]*Y[Hplus]
    -          rate[Z25]*Y[C]
    -          rate[Z26]*Y[Cplus]
    -          rate[UV4];

  J[O2][C] = -rate[Z25]*Y[O2]; 

  J[O2][Cplus]= -rate[Z26]*Y[O2]; 

/*********** Oplus ****************/

  J[Oplus][H] = -rate[Z1]*Y[Oplus]; 

  J[Oplus][Hplus] = rate[Z2]*Y[O]; 

  J[Oplus][O] = rate[Z2]*Y[Hplus]
   +            rate[UV8];

  J[Oplus][O2] = rate[Z26]*Y[Cplus]; 

  J[Oplus][Oplus] = -rate[Z1]*Y[H];

  J[Oplus][Cplus] = rate[Z26]*Y[O2]; 

/************* OHplus ***************/

  J[OHplus][Hplus] = rate[Z14]*Y[OH]; 

  J[OHplus][OH] = rate[Z14]*Y[Hplus]; 

/************ H2Oplus  ***************/

  J[H2Oplus][Hplus] = rate[Z15]*Y[H2O]; 

  J[H2Oplus][H2O] = rate[Z15]*Y[Hplus]; 

  J[H2Oplus][el] = -rate[Z7]*Y[H2Oplus]; 

  J[H2Oplus][H2Oplus] = -rate[Z7]*Y[el]; 

/*************** H3Oplus **************/

  J[H3Oplus][el] = -rate[Z8]*Y[H3Oplus] 
     -              rate[Z9]*Y[H3Oplus]; 

  J[H3Oplus][H3Oplus] = -rate[Z8]*Y[el] 
    -                    rate[Z9]*Y[el]; 

/**************** O2plus ****************/

  J[O2plus][Hplus] = rate[Z21]*Y[O2]; 

  J[O2plus][O2] = rate[Z21]*Y[Hplus]; 

  J[O2plus][O2plus] = -rate[Z22]*Y[el]; 

  J[O2plus][el] = -rate[Z22]*Y[O2plus]; 

/*********** COplus *******************/

  J[Cplus][el] = -rate[Z28]*Y[Cplus]; 

  J[Cplus][O2] = -rate[Z26]*Y[Cplus]; 

  J[Cplus][OH] = -rate[Z29]*Y[Cplus]; 
 
  J[Cplus][C]  = +rate[UV7];

  J[Cplus][Cplus] = -rate[Z26]*Y[O2] 
    -                rate[Z28]*Y[el] 
    -                rate[Z29]*Y[OH]; 

/*************** CO *******************/

  J[C][H] = -rate[Z31]*Y[C]; 

  J[C][el] = rate[Z28]*Y[Cplus]; 

  J[C][Cplus] = rate[Z28]*Y[el]; 
  
  J[C][OH] = -rate[Z24]*Y[C]; 

  J[C][O2] = -rate[Z25]*Y[C]; 

  J[C][C]  = -rate[Z24]*Y[OH] 
    -         rate[Z25]*Y[O2] 
    -         rate[Z31]*Y[H]
    -         rate[UV7];

  J[C][CH] =  rate[UV5]; 

  J[C][CO] =  rate[UV6];

/************** CH *******************/

  J[CH][H] = rate[Z31]*Y[C];

  J[CH][O] = -rate[Z23]*Y[CH]; 

  J[CH][C] = rate[Z31]*Y[H];

  J[CH][CH] = -rate[Z23]*Y[O]
  -            rate[UV5];

/************** CO *******************/

  J[CO][H] = rate[Z30]*Y[COplus];

  J[CO][O2] = rate[Z25]*Y[C] 
    +         rate[Z26]*Y[Cplus]; 

  J[CO][C] = rate[Z24]*Y[OH] 
    +        rate[Z25]*Y[O2]; 


  J[CO][O] = rate[Z23]*Y[CH]; 

  J[CO][OH] = rate[Z24]*Y[C] 
    -         rate[Z27]*Y[CO]; 

  J[CO][Cplus] = rate[Z26]*Y[O2]; 

  J[CO][CH] = rate[Z23]*Y[O]; 

  J[CO][COplus] = rate[Z30]*Y[H];   

  J[CO][CO] =-rate[Z27]*Y[OH]
    -         rate[UV6];  

/*********** COplus *******************/

  J[COplus][H] = -rate[Z30]*Y[COplus]; 

  J[COplus][OH] = rate[Z29]*Y[Cplus]; 

  J[COplus][Cplus] = rate[Z29]*Y[OH]; 

  J[COplus][COplus] = -rate[Z30]*Y[H]; 

/************* CO2 *********************/

  J[CO2][OH] = rate[Z27]*Y[CO]; 

  J[CO2][CO] = rate[Z27]*Y[OH]; 

/**************************************************/
/*              MULTISPECIES = 2                  */
/**************************************************/
  if (ispecies > 1)
  {
    J[H][H] += 3.0*rate[H4]*Y[H2m]
      +        3.0*rate[H5]*Y[H]*Y[H]
      -            rate[H2]*Y[el]
      -            rate[H3]*Y[Hmin]
      -            rate[H4]*Y[H2m]
      -        9.0*rate[H5]*Y[H]*Y[H]
      -        4.0*rate[H6]*Y[H2m]*Y[H]
      -        4.0*rate[H8]*Y[H]
      -            rate[Z16]*Y[OH]
      -            rate[Z17]*Y[H2O]
      -            rate[Z33]*Y[CH]
      -            rate[Z35]*Y[CH2]
      -            rate[Z37]*Y[CH3]
      -            rate[Z39]*Y[CH4];

    J[H][el] += -rate[H2]*Y[H];

    J[H][O] += rate[Z11]*Y[H2m];

    J[H][OH] += rate[Z12]*Y[H2m]
      -         rate[Z16]*Y[H];

    J[H][H2O] += -rate[Z17]*Y[H];

    J[H][Oplus] += rate[Z3]*Y[H2m];

    J[H][OHplus] = rate[Z4]*Y[H2m];

    J[H][H2Oplus] += rate[Z5]*Y[H2m];

    J[H][C] += rate[Z32]*Y[H2m];

    J[H][CH] += rate[Z34]*Y[H2m]
      -         rate[Z33]*Y[H];

    J[H][CH2] = rate[Z36]*Y[H2m]
      -         rate[Z35]*Y[H];

    J[H][CH3] = rate[Z38]*Y[H2m]
      -         rate[Z37]*Y[H];

    J[H][CH4] = -rate[Z39]*Y[H];

    J[H][Hmin] = - rate[H3]*Y[H];

    J[H][H2m] =
      +     3.0*rate[H4]*Y[H]
      +     4.0*rate[H7]*Y[H2m]
      +         rate[Z3]*Y[Oplus]
      +         rate[Z4]*Y[OHplus]
      +         rate[Z5]*Y[H2Oplus]
      +         rate[Z11]*Y[O]
      +         rate[Z12]*Y[OH]
      +         rate[Z32]*Y[C]
      +         rate[Z34]*Y[CH]
      +         rate[Z36]*Y[CH2]
      +         rate[Z38]*Y[CH3]
      -         rate[H4]*Y[H]
      -     2.0*rate[H6]*Y[H]*Y[H];

    J[el][H] = rate[H3]*Y[Hmin]
      -        rate[H2]*Y[el];

    J[el][el] += -rate[Z6]*Y[H2Oplus]
       -          rate[H2]*Y[H];

    J[el][H2Oplus] += -rate[Z6]*Y[el];

    J[el][Hmin] = rate[H3]*Y[H];

    J[O][H] += rate[Z16]*Y[OH];

    J[O][OH] += rate[Z16]*Y[H];

    J[O][el] += rate[Z6]*Y[H2Oplus];

    J[O][H2Oplus] = rate[Z6]*Y[el];

    J[O][O] += -rate[Z11]*Y[H2m];

    J[O][H2m] = -rate[Z11]*Y[O];

    J[OH][H] += rate[Z17]*Y[H2O]
      -         rate[Z16]*Y[OH];

    J[OH][OH] += -rate[Z12]*Y[H2m]
      -           rate[Z16]*Y[H];

    J[OH][H2O] = rate[Z17]*Y[H];

    J[OH][O] += rate[Z11]*Y[H2m];

    J[OH][H2m] = rate[Z11]*Y[O]
      -          rate[Z12]*Y[OH];

    J[H2O][H] = -rate[Z17]*Y[H2O];

    J[H2O][OH] += rate[Z12]*Y[H2m];

    J[H2O][H2O] += -rate[Z17]*Y[H];

    J[H2O][H2m] = rate[Z12]*Y[OH];

    J[Oplus][Oplus] += -rate[Z3]*Y[H2m];

    J[Oplus][H2m] = -rate[Z3]*Y[Oplus];

    J[OHplus][OHplus] = -rate[Z4]*Y[H2m];

    J[OHplus][Oplus] = rate[Z3]*Y[H2m];

    J[OHplus][H2m] = rate[Z3]*Y[Oplus]
      -              rate[Z4]*Y[OHplus];

    J[H2Oplus][OHplus] = rate[Z4]*Y[H2m];

    J[H2Oplus][el] += -rate[Z6]*Y[H2Oplus];

    J[H2Oplus][H2Oplus] += -rate[Z5]*Y[H2m]
     -                   rate[Z6]*Y[el];

    J[H2Oplus][H2m] = rate[Z4]*Y[OHplus]
      -               rate[Z5]*Y[H2Oplus];

    J[H3Oplus][H2Oplus] = rate[Z5]*Y[H2m];

      J[H3Oplus][H2m] = rate[Z5]*Y[H2Oplus];

    J[C][H] += rate[Z33]*Y[CH];

    J[C][CH] = rate[Z33]*Y[H];

     J[C][C] += -rate[Z32]*Y[H2m]
       -         rate[Z40]*Y[H2m];

    J[C][H2m] = -rate[Z32]*Y[C]
      -          rate[Z40]*Y[C];

    J[CH][H] += rate[Z35]*Y[CH2]
      -         rate[Z33]*Y[CH];

    J[CH][C] += rate[Z32]*Y[H2m];

    J[CH][CH] += -rate[Z33]*Y[H]
      -           rate[Z34]*Y[H2m];

    J[CH][CH2] = rate[Z35]*Y[H];

    J[CH][H2m] = rate[Z32]*Y[C]
      -          rate[Z34]*Y[CH];

    J[CH2][H] = rate[Z37]*Y[CH3]
      -         rate[Z35]*Y[CH2];

    J[CH2][C] = rate[Z40]*Y[H2m];

    J[CH2][CH] = rate[Z34]*Y[H2m];

    J[CH2][CH2] = -rate[Z35]*Y[H]
      -            rate[Z36]*Y[H2m];

    J[CH2][CH3] = rate[Z37]*Y[H];

    J[CH2][H2m] = rate[Z34]*Y[CH]
      +           rate[Z40]*Y[C]
      -           rate[Z36]*Y[CH2];

    J[CH3][H] = rate[Z39]*Y[CH4]
      -         rate[Z37]*Y[CH3];

    J[CH3][CH2] = rate[Z36]*Y[H2m];

    J[CH3][CH3] = -rate[Z38]*Y[H2m]
      -            rate[Z37]*Y[H];

    J[CH3][CH4] = rate[Z39]*Y[H];

    J[CH3][H2m] = rate[Z36]*Y[CH2]
      -           rate[Z38]*Y[CH3];

    J[CH4][H] = -rate[Z39]*Y[CH4];

    J[CH4][CH3] = rate[Z38]*Y[H2m];

    J[CH4][CH4] = -rate[Z39]*Y[H];

    J[CH4][H2m] = rate[Z38]*Y[CH3];

    J[H2m][H] = rate[H3]*Y[Hmin]  
      +     3.0*rate[H5]*Y[H]*Y[H]  
      +     4.0*rate[H6]*Y[H]*Y[H2m]
      +     2.0*rate[H8]*Y[H] 
      +         rate[Z16]*Y[OH]
      +         rate[Z17]*Y[H2O] 
      +         rate[Z33]*Y[CH] 
      +         rate[Z35]*Y[CH2] 
      +         rate[Z37]*Y[CH3] 
      +         rate[Z39]*Y[CH4]
      -         rate[H4]*Y[H2m] 
      -     2.0*rate[H6]*Y[H]*Y[H2m];

    J[H2m][Hmin] =  rate[H3]*Y[H];

    J[H2m][H2m] = 2.0*rate[H7]*Y[H2m] 
      +           2.0*rate[H6]*Y[H]*Y[H]
      -               rate[H4]*Y[H] 
      -               rate[H6]*Y[H]*Y[H] 
      -           4.0*rate[H7]*Y[H2m] 
      -               rate[Z3]*Y[Oplus] 
      -               rate[Z4]*Y[OHplus] 
      -               rate[Z5]*Y[H2Oplus]
      -               rate[Z11]*Y[O] 
      -               rate[Z12]*Y[OH] 
      -               rate[Z32]*Y[C] 
      -               rate[Z34]*Y[CH]
      -               rate[Z36]*Y[CH2] 
      -               rate[Z38]*Y[CH3] 
      -               rate[Z40]*Y[C];
 
    J[H2m][el] = rate[Z6]*Y[H2Oplus];
    
    J[H2m][O] = -rate[Z11]*Y[H2m];

    J[H2m][OH] = rate[Z16]*Y[H]
     -           rate[Z12]*Y[H2m];

    J[H2m][H2O] = rate[Z17]*Y[H]
     +            rate[UV2];

    J[H2m][Oplus] = -rate[Z3]*Y[H2m];
  
    J[H2m][OHplus] = -rate[Z4]*Y[H2m];

    J[H2m][H2Oplus] = rate[Z6]*Y[el] 
     -                rate[Z5]*Y[H2m];

    J[H2m][C] = -rate[Z32]*Y[H2m] 
     -           rate[Z40]*Y[H2m];

    J[H2m][CH] = rate[Z33]*Y[H] 
     -           rate[Z34]*Y[H2m];

    J[H2m][CH2] = rate[Z35]*Y[H] 
     -            rate[Z36]*Y[H2m];

    J[H2m][CH3] = rate[Z37]*Y[H] 
     -            rate[Z38]*Y[H2m];

    J[H2m][CH4] = rate[Z39]*Y[H];

    J[Hmin][H] = rate[H2]*Y[el] 
     -           rate[H3]*Y[Hmin];

    J[Hmin][Hmin] = -rate[H3]*Y[H];

    J[Hmin][el] = rate[H2]*Y[H];

/**************************************************/ 
/*              MULTISPECIES = 3                  */
/**************************************************/
    if (ispecies  > 2){

      J[D][H] = rate[D2]*Y[Dplus] + rate[D5]*Y[HD];

      J[D][Hplus] = -rate[D1]*Y[D];

      J[D][Dplus] = rate[D2]*Y[H];

      J[D][D] = -rate[D1]*Y[Hplus] - rate[D3]*Y[H2m];

      J[D][HD] = rate[D5]*Y[H];

      J[D][H2m] = -rate[D3]*Y[D];

      J[Dplus][H] = -rate[D2]*Y[Dplus];

      J[Dplus][Hplus] = rate[D1]*Y[D] + rate[D6]*Y[HD];

      J[Dplus][H2m] = -rate[D4]*Y[Dplus];

      J[Dplus][HD] = rate[D6]*Y[Hplus];

      J[Dplus][D] = rate[D1]*Y[Hplus];

      J[Dplus][Dplus] = -rate[D2]*Y[H] - rate[D4]*Y[H2m];

      J[HD][H] = -rate[D5]*Y[HD];

      J[HD][Hplus] =  -rate[D6]*Y[HD];

      J[HD][H2m] = rate[D3]*Y[D] + rate[D4]*Y[Dplus];

      J[HD][D] = rate[D3]*Y[H2m];

      J[HD][Dplus] = rate[D4]*Y[H2m];

      J[HD][HD] = -rate[D5]*Y[H] - rate[D6]*Y[Hplus];

      J[H2m][H] += rate[D5]*Y[HD];

      J[H2m][Hplus] = rate[D6]*Y[HD];

      J[H2m][D] = - rate[D3]*Y[H2m];

      J[H2m][Dplus] = -rate[D4]*Y[H2m];

      J[H2m][HD] = rate[D5]*Y[H]
        +          rate[D6]*Y[Hplus];

      J[H2m][H2m] += -rate[D3]*Y[D]
        -             rate[D4]*Y[Dplus];

      J[H][H] += -rate[D2]*Y[Dplus]
        -         rate[D5]*Y[HD];

      J[H][Hplus] += rate[D1]*Y[D];

      J[H][H2m] += rate[D3]*Y[D];

      J[H][D] = rate[D1]*Y[Hplus]
        +       rate[D3]*Y[H2m];

      J[H][Dplus] = -rate[D2]*Y[H];

      J[H][HD] = -rate[D5]*Y[H];

      J[Hplus][H] += rate[D2]*Y[Dplus];

      J[Hplus][Hplus] = -rate[D1]*Y[D]
        -                rate[D6]*Y[HD];

      J[Hplus][H2m] = rate[D4]*Y[Dplus];

      J[Hplus][Dplus] = rate[D2]*Y[H]
        +               rate[D4]*Y[H2m];

      J[Hplus][D] = -rate[D1]*Y[Hplus];

      J[Hplus][HD] = -rate[D6]*Y[Hplus];

      J[D][H] = rate[D2]*Y[Dplus]
        +       rate[D5]*Y[HD];

      J[D][Hplus] = -rate[D1]*Y[D];

      J[D][Dplus] = rate[D2]*Y[H];

      J[D][D] = -rate[D1]*Y[Hplus]
        -        rate[D3]*Y[H2m];

      J[D][HD] = rate[D5]*Y[H];

      J[D][H2m] = -rate[D3]*Y[D];

      J[Dplus][H] = -rate[D2]*Y[Dplus];

      J[Dplus][Hplus] = rate[D1]*Y[D]
       +                rate[D6]*Y[HD];

      J[Dplus][H2m] = -rate[D4]*Y[Dplus];

      J[Dplus][HD] = rate[D6]*Y[Hplus];

      J[Dplus][D] = rate[D1]*Y[Hplus];

      J[Dplus][Dplus] = -rate[D2]*Y[H]
         -               rate[D4]*Y[H2m];

      J[HD][H] = -rate[D5]*Y[HD];

      J[HD][Hplus] =  -rate[D6]*Y[HD];

      J[HD][H2m] = rate[D3]*Y[D]
       +           rate[D4]*Y[Dplus];

      J[HD][D] = rate[D3]*Y[H2m];

      J[HD][Dplus] = rate[D4]*Y[H2m];

      J[HD][HD] = -rate[D5]*Y[H]
        -          rate[D6]*Y[Hplus];

      }
    }
}
