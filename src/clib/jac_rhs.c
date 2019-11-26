#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <jac_rhs.h>
#include <chemistry.h>
#include <gsl/gsl_math.h>

// clip magnitude of jacobian elements to 1e-99 to avoid underflow error
void clip_jacobian(double **J, int ispecies, int water_rates){
  int i, j;
  double jabs, sign;

  for (i = 0; i < nSpecies; i++) {
    for (j = 0; j < nSpecies; j++) {
      //sign = (double) GSL_SIGN(J[i][j]);
      jabs = fabs(J[i][j]);
      if (jabs < 1.e-50){
        J[i][j] = 0.;
      }
    }
  }
}

void get_f_vector(double *rate, double *F, double *Y, int ispecies, int water_rates, int UV) {

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

  F[el] = -rate[H1]*Y[Hplus]*Y[el] 
    -     rate[Z7]*Y[H2Oplus]*Y[el] 
    -     rate[Z8]*Y[H3Oplus]*Y[el] 
    -     rate[Z9]*Y[H3Oplus]*Y[el]
    -     rate[Z22]*Y[O2plus]*Y[el] 
    -     rate[Z28]*Y[Cplus]*Y[el]
    +     rate[UV7]*Y[C]
    +     rate[UV8]*Y[O];
    //+     rate[UV9]*Y[Hmin];

/************* F[Hplus] ******************/

  F[Hplus] = rate[Z1]*Y[Oplus]*Y[H] 
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
    +     rate[UV1]*Y[OH]
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
   -      rate[Z27]*Y[OH]*Y[CO]
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
    -           rate[Z39]*Y[CH4]*Y[H]
    +           rate[UV9]*Y[Hmin];

    F[el] += rate[H3]*Y[Hmin]*Y[H]
    +        rate[UV9]*Y[Hmin]
    -        rate[H2]*Y[H]*Y[el]
    -        rate[Z6]*Y[H2Oplus]*Y[el];

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

    F[CH3] += rate[Z36]*Y[H2m]*Y[CH2]
     +       rate[Z39]*Y[H]*Y[CH4]
     -       rate[Z38]*Y[H2m]*Y[CH3]
     -       rate[Z37]*Y[H]*Y[CH3];

    F[CH4] += rate[Z38]*Y[H2m]*Y[CH3]
     -       rate[Z39]*Y[H]*Y[CH4];

    F[H2m] += rate[H3]*Y[Hmin]*Y[H] 
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

      F[Hmin] += rate[H2]*Y[H]*Y[el] 
              - rate[H3]*Y[Hmin]*Y[H]
              - rate[UV9]*Y[Hmin];

/**************************************************/
/*              MULTISPECIES = 3                  */
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
  
  if (water_rates == 3){
     F[C] += rate[CRX2]*Y[CO] 
      +      rate[CRX3]*Y[CHplus]
      +      rate[CRX4]*Y[CH]
      -      rate[CRX1]*Y[C]
      +      rate[Z68]*Y[CHplus]*Y[el]
      +      rate[Z70]*Y[CH2plus]*Y[el]
      +      rate[Z71]*Y[CH2plus]*Y[el]
      +      rate[Z73]*Y[H]*Y[CH]
      -      rate[Z55]*Y[C]*Y[OH]
      -      rate[Z61]*Y[H3plus]*Y[C];

     F[O] +=  rate[Z45]*Y[Heplus]*Y[CO] 
      +       rate[Z46]*Y[O2]*Y[Cplus]
      +       rate[Z50]*Y[H3Oplus]*Y[el] 
      +       rate[Z51]*Y[H2Oplus]*Y[el] 
      +       rate[Z52]*Y[OHplus]*Y[el]
      +       rate[CRX2]*Y[CO]
      +       rate[CRX6]*Y[OH]
      +   2.0*rate[CRX9]*Y[O2]
      -       rate[Z43]*Y[O]*Y[CH2] 
      -       rate[Z44]*Y[O]*Y[CH2] 
      -       rate[Z48]*Y[O]*Y[H3plus]
      -       rate[Z76]*Y[O]*Y[CH];


     F[Cplus] +=  rate[CRX1]*Y[C]
      +           rate[Z45]*Y[Heplus]*Y[CO]
      -           rate[Z46]*Y[O2]*Y[Cplus]
      +           rate[Z72]*Y[H]*Y[CHplus]
      -           rate[Z56]*Y[Cplus]*Y[OH]
      -           rate[Z58]*Y[Cplus]*Y[H2O]
      -           rate[Z63]*Y[H2m]*Y[Cplus];
 
     F[H2m] += rate[Z43]*Y[O]*Y[CH2] 
      +       rate[Z48]*Y[O]*Y[H3plus] 
      +       rate[Z49]*Y[H3Oplus]*Y[el] 
      +       rate[Z50]*Y[H3Oplus]*Y[el] 
      +       rate[CRX7]*Y[CH4]
      +       rate[CRX16]*Y[CH3]
      +       rate[H11]*Y[H3plus]*Y[el]
      +       rate[Z53]*Y[H3plus]*Y[CO]
      +       rate[Z60]*Y[H3plus]*Y[H2O]
      +       rate[Z61]*Y[H3plus]*Y[C]
      +       rate[Z66]*Y[CH3plus]*Y[el]
      +       rate[Z71]*Y[CH2plus]*Y[el]
      +       rate[Z72]*Y[H]*Y[CHplus]
      +       rate[Z73]*Y[H]*Y[CH]
      -       rate[H9]*Y[H2plus]*Y[H2m]
      -       rate[H13]*Y[H2m]*Y[Heplus]
      -       rate[H14]*Y[H2m]*Y[Heplus]
      -       rate[H16]*Y[H2m]*Y[HeHplus]
      -       rate[Z62]*Y[H2m]*Y[CHplus]
      -       rate[Z63]*Y[H2m]*Y[Cplus]
      -       rate[Z64]*Y[H2m]*Y[CH2plus]
      -       rate[Z74]*Y[H2m]*Y[COplus];

     F[H] += 2.0*rate[Z44]*Y[O]*Y[CH2] 
      +      2.0*rate[Z51]*Y[H2Oplus]*Y[el] 
      +          rate[Z52]*Y[OHplus]*Y[el]
      +          rate[CRX4]*Y[CH]
      +          rate[CRX5]*Y[H2O]
      +          rate[CRX6]*Y[OH]
      +          rate[CRX14]*Y[CH2]
      +          rate[CRX15]*Y[CH3]
      +          rate[H9]*Y[H2plus]*Y[H2m]
      +      3.0*rate[H10]*Y[H3plus]*Y[el]	
      +          rate[H11]*Y[H3plus]*Y[el]
      +          rate[H14]*Y[H2m]*Y[Heplus]
      +          rate[H15]*Y[H2plus]*Y[He]
      +          rate[Z50]*Y[H3Oplus]*Y[el]
      +      2.0*rate[Z54]*Y[H3Oplus]*Y[el]
      +          rate[Z55]*Y[C]*Y[OH]
      +          rate[Z56]*Y[Cplus]*Y[OH]
      +          rate[Z57]*Y[Hplus]*Y[OH]
      +          rate[Z58]*Y[Cplus]*Y[H2O]
      +          rate[Z62]*Y[H2m]*Y[CHplus]
      +          rate[Z64]*Y[H2m]*Y[CH2plus]
      +          rate[Z65]*Y[CH3plus]*Y[el]
      +      2.0*rate[Z67]*Y[CH3plus]*Y[el]
      +          rate[Z68]*Y[CHplus]*Y[el]
      +          rate[Z69]*Y[CH2plus]*Y[el]
      +      2.0*rate[Z70]*Y[CH2plus]*Y[el]
      +          rate[Z74]*Y[H2m]*Y[COplus]
      +          rate[Z75]*Y[HCOplus]*Y[el]
      -          rate[H12]*Y[H]*Y[Heplus]
      -          rate[H17]*Y[H]*Y[HeHplus]
      -          rate[Z72]*Y[H]*Y[CHplus]
      -          rate[Z73]*Y[H]*Y[CH];
 
     F[CHplus] = -rate[CRX3]*Y[CHplus]
      +           rate[Z61]*Y[H3plus]*Y[C]
      +           rate[UV10]*Y[CH]
      -           rate[Z62]*Y[H2m]*Y[CHplus]	
      -           rate[Z68]*Y[CHplus]*Y[el]
      -           rate[Z72]*Y[H]*Y[CHplus]; 

     F[el] +=  rate[CRX1]*Y[C]
      +        rate[CRX8]*Y[O2]
      +        rate[CRX13]*Y[CH2] 
      +        rate[Z76]*Y[CH]*Y[O]
      +        rate[CRX17]*Y[He]
      +        rate[UV10]*Y[CH]
      -        rate[Z49]*Y[H3Oplus]*Y[el] 
      -        rate[Z50]*Y[H3Oplus]*Y[el]
      -        rate[Z51]*Y[H2Oplus]*Y[el]
      -        rate[Z52]*Y[OHplus]*Y[el]
      -        rate[H10]*Y[H3plus]*Y[el]
      -        rate[H11]*Y[H3plus]*Y[el]
      -        rate[Z54]*Y[H3Oplus]*Y[el]
      -        rate[Z65]*Y[CH3plus]*Y[el]
      -        rate[Z66]*Y[CH3plus]*Y[el]
      -        rate[Z67]*Y[CH3plus]*Y[el]
      -        rate[Z68]*Y[CHplus]*Y[el]
      -        rate[Z69]*Y[CH2plus]*Y[el]
      -        rate[Z70]*Y[CH2plus]*Y[el]
      -        rate[Z71]*Y[CH2plus]*Y[el]
      -        rate[Z75]*Y[HCOplus]*Y[el];

 
     F[OHplus] += rate[Z48]*Y[O]*Y[H3plus]
      +           rate[Z57]*Y[Hplus]*Y[OH];
      -           rate[Z52]*Y[OHplus]*Y[el];
  
     F[H2O] += -rate[CRX5]*Y[H2O]
      -         rate[Z58]*Y[Cplus]*Y[H2O]
      -         rate[Z59]*Y[H2O]*Y[HCOplus]
      -         rate[Z60]*Y[H3plus]*Y[H2O];
 
     F[H3Oplus] += -rate[Z49]*Y[H3Oplus]*Y[el] 
      -             rate[Z50]*Y[H3Oplus]*Y[el]
      +             rate[Z59]*Y[H2O]*Y[HCOplus]
      +             rate[Z60]*Y[H3plus]*Y[H2O]
      -             rate[Z54]*Y[H3Oplus]*Y[el];

     F[CH] += rate[CRX14]*Y[CH2]
      +       rate[CRX16]*Y[CH3]
      -       rate[CRX4]*Y[CH]
      +       rate[Z66]*Y[CH3plus]*Y[el]
      +       rate[Z67]*Y[CH3plus]*Y[el]
      +       rate[Z69]*Y[CH2plus]*Y[el]
      -       rate[Z73]*Y[H]*Y[CH]
      -       rate[Z76]*Y[CH]*Y[O]
      -       rate[UV10]*Y[CH];

     F[CH2] +=  rate[CRX7]*Y[CH4]
      +         rate[CRX15]*Y[CH3]
      +         rate[Z65]*Y[CH3plus]*Y[el]
      -         rate[Z43]*Y[O]*Y[CH2] 
      -         rate[Z44]*Y[O]*Y[CH2]
      -         rate[CRX13]*Y[CH2]
      -         rate[CRX14]*Y[CH2];
 
     F[CH3] += -rate[CRX15]*Y[CH3]
      -         rate[CRX16]*Y[CH3];

     F[CH3plus] = rate[Z64]*Y[H2m]*Y[CH2plus]
      -           rate[Z65]*Y[CH3plus]*Y[el]
      -           rate[Z66]*Y[CH3plus]*Y[el]
      -           rate[Z67]*Y[CH3plus]*Y[el];

     F[CH2plus] = rate[CRX13]*Y[CH2]
      +           rate[Z62]*Y[H2m]*Y[CHplus]
      +           rate[Z63]*Y[H2m]*Y[Cplus]
      -           rate[Z64]*Y[H2m]*Y[CH2plus]
      -           rate[Z69]*Y[CH2plus]*Y[el]
      -           rate[Z70]*Y[CH2plus]*Y[el]
      -           rate[Z71]*Y[CH2plus]*Y[el];
 
     F[CH4] += -rate[CRX7]*Y[CH4];
 
     F[He] = rate[Z45]*Y[Heplus]*Y[CO]
      +      rate[H12]*Y[H]*Y[Heplus]
      +      rate[H13]*Y[H2m]*Y[Heplus]
      +      rate[H14]*Y[H2m]*Y[Heplus]
      +      rate[H16]*Y[H2m]*Y[HeHplus]
      +      rate[H17]*Y[H]*Y[HeHplus]
      -      rate[H15]*Y[H2plus]*Y[He]
      -      rate[H18]*Y[Hplus]*Y[He];
      -      rate[CRX17]*Y[He];
 
     F[Heplus] = -rate[Z45]*Y[Heplus]*Y[CO]
      +           rate[CRX17]*Y[He]
      -           rate[H12]*Y[H]*Y[Heplus]
      -           rate[H13]*Y[H2m]*Y[Heplus]
      -           rate[H14]*Y[H2m]*Y[Heplus];	 

     F[CO] += rate[Z43]*Y[O]*Y[CH2] 
      +       rate[Z44]*Y[O]*Y[CH2] 
      -       rate[Z45]*Y[Heplus]*Y[CO]
      -       rate[CRX2]*Y[CO]
      +       rate[Z55]*Y[C]*Y[OH]
      +       rate[Z59]*Y[H2O]*Y[HCOplus]
      +       rate[Z75]*Y[HCOplus]*Y[el]
      -       rate[Z53]*Y[H3plus]*Y[CO];
 
     F[COplus] += rate[Z46]*Y[O2]*Y[Cplus]
      +           rate[Z56]*Y[Cplus]*Y[OH]
      -           rate[Z74]*Y[H2m]*Y[COplus];

     F[O2] += -rate[Z46]*Y[O2]*Y[Cplus] 
      -        rate[CRX8]*Y[O2]
      -        rate[CRX9]*Y[O2];

     F[O2plus] += rate[CRX8]*Y[O2];

     F[H2plus] =  rate[H13]*Y[H2m]*Y[Heplus]
      +           rate[H17]*Y[H]*Y[HeHplus]
      -           rate[H9]*Y[H2plus]*Y[H2m]
      -           rate[H15]*Y[H2plus]*Y[He];
 
     F[H3plus] =  -rate[Z48]*Y[O]*Y[H3plus]
      +            rate[H9]*Y[H2plus]*Y[H2m]
      +            rate[H16]*Y[H2m]*Y[HeHplus]
      -            rate[H10]*Y[H3plus]*Y[el]
      -            rate[H11]*Y[H3plus]*Y[el]
      -            rate[Z53]*Y[H3plus]*Y[CO]
      -            rate[Z60]*Y[H3plus]*Y[H2O]
      -            rate[Z61]*Y[H3plus]*Y[C];

     F[OH] += rate[Z49]*Y[H3Oplus]*Y[el]
      +       rate[CRX5]*Y[H2O]
      -       rate[CRX6]*Y[OH]
      +       rate[Z54]*Y[H3Oplus]*Y[el]
      -       rate[Z55]*Y[C]*Y[OH]
      -       rate[Z56]*Y[Cplus]*Y[OH]
      -       rate[Z57]*Y[Hplus]*Y[OH];

     F[H2Oplus] += -rate[Z51]*Y[H2Oplus]*Y[el];
 
     F[Hplus]   +=  rate[CRX3]*Y[CHplus]
      +             rate[H12]*Y[H]*Y[Heplus]
      +             rate[H14]*Y[H2m]*Y[Heplus]
      -             rate[H18]*Y[Hplus]*Y[He]
      -             rate[Z57]*Y[Hplus]*Y[OH];

     F[HeHplus] =  rate[H15]*Y[H2plus]*Y[He]
      +             rate[H18]*Y[Hplus]*Y[He]
      -             rate[H16]*Y[H2m]*Y[HeHplus]
      -             rate[H17]*Y[H]*Y[HeHplus];
 
     F[HCOplus] = rate[Z53]*Y[H3plus]*Y[CO]
      +            rate[Z58]*Y[Cplus]*Y[H2O]
      +            rate[Z74]*Y[H2m]*Y[COplus]
      +            rate[Z76]*Y[CH]*Y[O]
      -            rate[Z59]*Y[H2O]*Y[HCOplus]
      -            rate[Z75]*Y[HCOplus]*Y[el];

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
void get_jacobian(double rate[], double **J, double Y[], int ispecies, int water_rates, int UV) {
  // We are taking the derivative of the rate equations
  //   \dot{y} = f(y)

/**************************************************/
/*              MULTISPECIES = 1                  */
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

/*********** Cplus *******************/

  J[Cplus][el] = -rate[Z28]*Y[Cplus]; 

  J[Cplus][O2] = -rate[Z26]*Y[Cplus]; 

  J[Cplus][OH] = -rate[Z29]*Y[Cplus]; 
 
  J[Cplus][C]  = +rate[UV7];

  J[Cplus][Cplus] = -rate[Z26]*Y[O2] 
    -                rate[Z28]*Y[el] 
    -                rate[Z29]*Y[OH]; 

/*************** C *******************/

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
    /*********** H *********************/

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

    J[H][OHplus] += rate[Z4]*Y[H2m];

    J[H][H2Oplus] += rate[Z5]*Y[H2m];

    J[H][C] += rate[Z32]*Y[H2m];

    J[H][CH] += rate[Z34]*Y[H2m]
      -         rate[Z33]*Y[H];

    J[H][CH2] += rate[Z36]*Y[H2m]
      -         rate[Z35]*Y[H];

    J[H][CH3] += rate[Z38]*Y[H2m]
      -         rate[Z37]*Y[H];

    J[H][CH4] += -rate[Z39]*Y[H];

    J[H][Hmin] += -rate[H3]*Y[H]
                 +rate[UV9];

    J[H][H2m] +=
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

    /*********** el ********************/

    J[el][H] += rate[H3]*Y[Hmin]
      -         rate[H2]*Y[el];

    J[el][el] +=-rate[Z6]*Y[H2Oplus]
      -          rate[H2]*Y[H];

    J[el][H2Oplus] += -rate[Z6]*Y[el];

    J[el][Hmin] += rate[H3]*Y[H]
      +           rate[UV9];

    /*********** O *********************/

    J[O][H] += rate[Z16]*Y[OH];

    J[O][OH] += rate[Z16]*Y[H];

    J[O][el] += rate[Z6]*Y[H2Oplus];

    J[O][H2Oplus] += rate[Z6]*Y[el];

    J[O][O] += -rate[Z11]*Y[H2m];

    J[O][H2m] += -rate[Z11]*Y[O];

    /*********** OH ********************/

    J[OH][H] += rate[Z17]*Y[H2O]
      -         rate[Z16]*Y[OH];

    J[OH][OH] += -rate[Z12]*Y[H2m]
      -           rate[Z16]*Y[H];

    J[OH][H2O] = rate[Z17]*Y[H];

    J[OH][O] += rate[Z11]*Y[H2m];

    J[OH][H2m] = rate[Z11]*Y[O]
      -          rate[Z12]*Y[OH];

    /*********** H2O *******************/

    J[H2O][H] = -rate[Z17]*Y[H2O];

    J[H2O][OH] += rate[Z12]*Y[H2m];

    J[H2O][H2O] += -rate[Z17]*Y[H];

    J[H2O][H2m] = rate[Z12]*Y[OH];

    /*********** Oplus *****************/

    J[Oplus][Oplus] += -rate[Z3]*Y[H2m];

    J[Oplus][H2m] = -rate[Z3]*Y[Oplus];

    /*********** OHplus ****************/

    J[OHplus][OHplus] = -rate[Z4]*Y[H2m];

    J[OHplus][Oplus] = rate[Z3]*Y[H2m];

    J[OHplus][H2m] = rate[Z3]*Y[Oplus]
      -              rate[Z4]*Y[OHplus];

    /*********** H2Oplus ***************/

    J[H2Oplus][OHplus] = rate[Z4]*Y[H2m];

    J[H2Oplus][el] += -rate[Z6]*Y[H2Oplus];

    J[H2Oplus][H2Oplus] += -rate[Z5]*Y[H2m]
     -                   rate[Z6]*Y[el];

    J[H2Oplus][H2m] = rate[Z4]*Y[OHplus]
      -               rate[Z5]*Y[H2Oplus];

    /*********** H3Oplus ***************/

    J[H3Oplus][H2Oplus] = rate[Z5]*Y[H2m];

    J[H3Oplus][H2m] = rate[Z5]*Y[H2Oplus];

    /*********** C *********************/

    J[C][H] += rate[Z33]*Y[CH];

    J[C][CH] = rate[Z33]*Y[H];

     J[C][C] += -rate[Z32]*Y[H2m]
       -         rate[Z40]*Y[H2m];

    J[C][H2m] = -rate[Z32]*Y[C]
      -          rate[Z40]*Y[C];

    /*********** CH ********************/

    J[CH][H] += rate[Z35]*Y[CH2]
      -         rate[Z33]*Y[CH];

    J[CH][C] += rate[Z32]*Y[H2m];

    J[CH][CH] += -rate[Z33]*Y[H]
      -           rate[Z34]*Y[H2m];

    J[CH][CH2] = rate[Z35]*Y[H];

    J[CH][H2m] = rate[Z32]*Y[C]
      -          rate[Z34]*Y[CH];

    /*********** CH2 *******************/

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

    /*********** CH3 *******************/

    J[CH3][H] = rate[Z39]*Y[CH4]
      -         rate[Z37]*Y[CH3];

    J[CH3][CH2] = rate[Z36]*Y[H2m];

    J[CH3][CH3] = -rate[Z38]*Y[H2m]
      -            rate[Z37]*Y[H];

    J[CH3][CH4] = rate[Z39]*Y[H];

    J[CH3][H2m] = rate[Z36]*Y[CH2]
      -           rate[Z38]*Y[CH3];

    /*********** CH4 *******************/

    J[CH4][H] = -rate[Z39]*Y[CH4];

    J[CH4][CH3] = rate[Z38]*Y[H2m];

    J[CH4][CH4] = -rate[Z39]*Y[H];

    J[CH4][H2m] = rate[Z38]*Y[CH3];

    /*********** H2m *******************/

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

    J[H2m][Oplus] += -rate[Z3]*Y[H2m];
  
    J[H2m][OHplus] += -rate[Z4]*Y[H2m];

    J[H2m][H2Oplus] += rate[Z6]*Y[el] 
     -                rate[Z5]*Y[H2m];

    J[H2m][C] += -rate[Z32]*Y[H2m] 
     -           rate[Z40]*Y[H2m];

    J[H2m][CH] += rate[Z33]*Y[H] 
     -           rate[Z34]*Y[H2m];

    J[H2m][CH2] += rate[Z35]*Y[H] 
     -            rate[Z36]*Y[H2m];

    J[H2m][CH3] += rate[Z37]*Y[H] 
     -            rate[Z38]*Y[H2m];

    J[H2m][CH4] += rate[Z39]*Y[H];

    /*********** Hmin ******************/

    J[Hmin][H] += rate[H2]*Y[el] 
     -           rate[H3]*Y[Hmin];

    J[Hmin][Hmin] += -rate[H3]*Y[H]
                    -rate[UV9];

    J[Hmin][el] += rate[H2]*Y[H];

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

      }
    }

    if (water_rates == 3){
      /***************** C **************/ 

      J[C][C] += -rate[CRX1]
       -          rate[Z55]*Y[OH]
       -          rate[Z61]*Y[H3plus];


      J[C][CO] += rate[CRX2];
 
      J[C][CHplus] += rate[CRX3]
       +              rate[Z68]*Y[el];

      J[C][CH] += rate[CRX4]
       +          rate[Z73]*Y[H];

      J[C][el] += rate[Z68]*Y[CHplus]
       +          rate[Z70]*Y[CH2plus]
       +          rate[Z71]*Y[CH2plus];

      J[C][CH2plus] += rate[Z70]*Y[el]
       +               rate[Z71]*Y[el];

      J[C][H] += rate[Z73]*Y[CH];

      J[C][OH] -= rate[Z55]*Y[C];

      J[C][H3plus] -= rate[Z61]*Y[C]; 

      /***************** O **************/ 

      J[O][O]  += -rate[Z43]*Y[CH2] 
        -          rate[Z44]*Y[CH2]
        -          rate[Z48]*Y[H3plus]
        -          rate[Z76]*Y[CH];

      J[O][OH] += rate[CRX6];

      J[O][CH] += -rate[Z76]*Y[O];

      J[O][CH2] += -rate[Z43]*Y[O]  
        -            rate[Z44]*Y[O];
 
      J[O][Heplus] += rate[Z45]*Y[CO];

      J[O][CO] += rate[Z45]*Y[Heplus]
        +         rate[CRX2];

      J[O][O2] += rate[Z46]*Y[Cplus] 
        +     2.0*rate[CRX9];
 
      J[O][Cplus] += rate[Z46]*Y[O2];
     
      J[O][H3plus] += -rate[Z48]*Y[O];
 
      J[O][el] += rate[Z50]*Y[H3Oplus] 
        +          rate[Z51]*Y[H2Oplus]
        +          rate[Z52]*Y[OHplus];

      J[O][H3Oplus] += rate[Z50]*Y[el];

      J[O][H2Oplus] += 
       +               rate[Z51]*Y[el];

      J[O][OHplus] += rate[Z52]*Y[el];
 
      /***************** Cplus **********/ 

      J[Cplus][C] += rate[CRX1];

      J[Cplus][Cplus] += 
       -                  rate[Z46]*Y[O2]
       -                  rate[Z56]*Y[OH]
       -                  rate[Z58]*Y[H2O]
       -                  rate[Z63]*Y[H2m];

      J[Cplus][CO] += rate[Z45]*Y[CO];

      J[Cplus][O2] += -rate[Z46]*Y[Cplus];

      J[Cplus][H] += rate[Z72]*Y[CHplus];

      J[Cplus][Heplus] += rate[Z45]*Y[CO];

      J[Cplus][CHplus] += rate[Z72]*Y[H]; 

      J[Cplus][OH] -= rate[Z56]*Y[Cplus];

      J[Cplus][H2O] -= rate[Z58]*Y[Cplus];

      J[Cplus][H2m] -= rate[Z63]*Y[Cplus]; 

      /***************** H2m ************/ 

      J[H2m][H] += rate[Z72]*Y[CHplus]
       +           rate[Z73]*Y[CH];

      J[H2m][el] += rate[Z49]*Y[H3Oplus] 
       +            rate[Z50]*Y[H3Oplus]
       +            rate[H11]*Y[H3plus]
       +            rate[Z66]*Y[CH3plus]
       +            rate[Z71]*Y[CH2plus];

      J[H2m][H2m]-= rate[H9]*Y[H2plus]
       -            rate[H13]*Y[Heplus]
       -            rate[H14]*Y[Heplus]
       -            rate[H16]*Y[HeHplus]
       -            rate[Z62]*Y[CHplus]
       -            rate[Z63]*Y[Cplus]
       -            rate[Z64]*Y[CH2plus]
       -            rate[Z74]*Y[COplus];

      J[H2m][C] += rate[Z61]*Y[H3plus];

      J[H2m][O] += rate[Z43]*Y[CH2] 
       +          rate[Z48]*Y[H3plus];

      J[H2m][Heplus] -= rate[H13]*Y[H2m]
       -                rate[H14]*Y[H2m];

      J[H2m][Cplus] -= rate[Z63]*Y[H2m];

      J[H2m][HeHplus] -= rate[H16]*Y[H2m];

      J[H2m][COplus] -= rate[Z74]*Y[H2m];

      J[H2m][CH] += rate[Z73]*Y[H];

      J[H2m][H2plus] -= rate[H9]*Y[H2m];

      J[H2m][CHplus] += rate[Z72]*Y[H]
       -                rate[Z62]*Y[H2m];

     J[H2m][CH2] += rate[Z43]*Y[O];

     J[H2m][CH2plus] += rate[Z71]*Y[el]
      -                 rate[Z64]*Y[H2m];

     J[H2m][CH3] += rate[CRX16];

     J[H2m][CH3plus] += rate[Z66]*Y[el];
 
     J[H2m][CH4] += rate[CRX7];

     J[H2m][H2O] += rate[Z60]*Y[H3plus];

     J[H2m][H3plus] += rate[Z48]*Y[O]
      +                rate[H11]*Y[el]
      +                rate[Z53]*Y[CO]
      +                rate[Z60]*Y[H2O]
      +                rate[Z61]*Y[C];
 
     J[H2m][CO] += rate[Z53]*Y[H3plus];
     
     J[H2m][H3Oplus] += rate[Z49]*Y[el] 
      +                rate[Z50]*Y[el];
 
     /***************** H **************/ 

     J[H][H] -= rate[H12]*Y[Heplus]
      -         rate[H17]*Y[HeHplus]
      -         rate[Z72]*Y[CHplus]
      -         rate[Z73]*Y[CH];
 
     J[H][Hplus] += rate[Z57]*Y[OH];

     J[H][H2plus] += rate[H9]*Y[H2m]
      +              rate[H15]*Y[He];

     J[H][H2m] += rate[H9]*Y[H2plus]
      +           rate[H14]*Y[Heplus]
      +           rate[Z62]*Y[CHplus]
      +           rate[Z64]*Y[CH2plus]
      +           rate[Z74]*Y[COplus];
 
     J[H][O] += 2.0*rate[Z44]*Y[CH2];

     J[H][C] += rate[Z55]*Y[OH];

     J[H][Cplus] += rate[Z56]*Y[OH]
      +             rate[Z58]*Y[H2O];

     J[H][COplus] += rate[Z74]*Y[H2m];

     J[H][CH] += rate[CRX4]
      -          rate[Z73]*Y[H];

     J[H][CHplus] += rate[Z62]*Y[H2m]
      +              rate[Z68]*Y[el]
      -              rate[Z72]*Y[H];
 
     J[H][CH2] += 2.0*rate[Z44]*Y[O] 
      +               rate[CRX14];

     J[H][CH2plus] += rate[Z64]*Y[H2m]
      +               rate[Z69]*Y[el]
      +           2.0*rate[Z70]*Y[el];
  
     J[H][CH3] += rate[CRX15];

     J[H][CH3plus] += rate[Z65]*Y[el]
      +           2.0*rate[Z67]*Y[el];
  
     J[H][HCOplus] += rate[Z75]*Y[el];

     J[H][H2O] += rate[CRX5];
      +           rate[Z58]*Y[Cplus];

     J[H][H3Oplus] += rate[Z50]*Y[el]
      +           2.0*rate[Z54]*Y[el];

     J[H][H3plus] += 3.0*rate[H10]*Y[el]	
      +                  rate[H11]*Y[el];

     J[H][el] += 
      +      2.0*rate[Z51]*Y[H2Oplus]
      +          rate[Z52]*Y[OHplus]
      +      3.0*rate[H10]*Y[H3plus]	
      +          rate[H11]*Y[H3plus]
      +          rate[Z50]*Y[H3Oplus]
      +      2.0*rate[Z54]*Y[H3Oplus]
      +          rate[Z65]*Y[CH3plus]
      +      2.0*rate[Z67]*Y[CH3plus]
      +          rate[Z68]*Y[CHplus]
      +          rate[Z69]*Y[CH2plus]
      +      2.0*rate[Z70]*Y[CH2plus]
      +          rate[Z75]*Y[HCOplus];

     J[H][He] += rate[H15]*Y[H2plus];

     J[H][Heplus] += - rate[H12]*Y[H]
                     + rate[H14]*Y[H2m];
  
     J[H][HeHplus] -= rate[H17]*Y[H];

     J[H][H2Oplus] += 2.0*rate[Z51]*Y[el];

     J[H][OH] += rate[CRX6]
      +          rate[Z55]*Y[C]
      +          rate[Z56]*Y[Cplus]
      +          rate[Z57]*Y[Hplus];

     J[H][OHplus] += rate[Z52]*Y[el];

     /***************** H2plus *********/ 

     J[H2plus][H2m] += rate[H13]*Y[H3plus]
      -                rate[H9]*Y[H2plus];

     J[H2plus][Heplus] += rate[H13]*Y[H2m]; 

     J[H2plus][H] += rate[H17]*Y[HeHplus];

     J[H2plus][HeHplus] += rate[H17]*Y[H];

     J[H2plus][H2plus] -= rate[H9]*Y[H2m]
      -                   rate[H15]*Y[He];

     J[H2plus][He] -= rate[H15]*Y[H2plus];

     /***************** Hplus **********/ 

     J[Hplus][CHplus] += rate[CRX3];

     J[Hplus][H] += rate[H12]*Y[Heplus]; 

     J[Hplus][Heplus] += rate[H12]*Y[H]
      +                  rate[H14]*Y[H2m];

     J[Hplus][H2m] += rate[H14]*Y[Heplus];

     J[Hplus][Hplus] -= rate[H18]*Y[He]
      -                 rate[Z57]*Y[OH]; 

     J[Hplus][He] -= rate[H18]*Y[Hplus];

     J[Hplus][OH] -= rate[Z57]*Y[Hplus];

     /***************** HeHplus ********/ 

     J[HeHplus][H2plus] += rate[H15]*Y[He];

     J[HeHplus][He] += rate[H15]*Y[H2plus]
      +                rate[H18]*Y[Hplus]; 

     J[HeHplus][Hplus] += rate[H18]*Y[He];

     J[HeHplus][H2m] -= rate[H16]*Y[HeHplus];

     J[HeHplus][HeHplus] -= rate[H16]*Y[H2m]
      -                     rate[H17]*Y[H];

     J[HeHplus][H] -= rate[H17]*Y[HeHplus]; 

     /***************** el *************/ 

     J[el][C] += rate[CRX1];

     J[el][He] += rate[CRX17];

     J[el][O] += rate[Z76]*Y[CH];
  
     J[el][O2] += rate[CRX8];
 
     J[el][CH] += rate[Z76]*Y[O] 
      +          rate[UV10];

     J[el][CH2]  += rate[CRX13];

     J[el][CH2plus] += -rate[Z69]*Y[el]	
      -                rate[Z70]*Y[el]	
      -                rate[Z71]*Y[el];
 
     J[el][el] += -rate[Z49]*Y[H3Oplus] 
      -            rate[Z50]*Y[H3Oplus]
      -            rate[Z51]*Y[H2Oplus]
      -            rate[Z52]*Y[OHplus]
      -            rate[H10]*Y[H3plus]
      -            rate[H11]*Y[H3plus]
      -            rate[Z54]*Y[H3Oplus]
      -            rate[Z65]*Y[CH3plus]
      -            rate[Z66]*Y[CH3plus]
      -            rate[Z67]*Y[CH3plus]	
      -            rate[Z68]*Y[CHplus]	
      -            rate[Z69]*Y[CH2plus]	
      -            rate[Z70]*Y[CH2plus]	
      -            rate[Z71]*Y[CH2plus]
      -            rate[Z75]*Y[HCOplus];

     J[el][H3plus] -= rate[H10]*Y[el]
      -               rate[H11]*Y[el];

     J[el][H3Oplus] += -rate[Z49]*Y[el] 
      -                 rate[Z50]*Y[el]
      -                 rate[Z54]*Y[el];
 
     J[el][H2Oplus] += 
      -                 rate[Z51]*Y[el];
 
     J[el][OHplus] += -rate[Z52]*Y[el];

     J[el][CH3plus] -= rate[Z65]*Y[el]
      -                rate[Z66]*Y[el]
      -                rate[Z67]*Y[el];

     J[el][CHplus] -= rate[Z68]*Y[el];

     J[el][HCOplus] -= rate[Z75]*Y[el];

     /***************** OHplus *********/ 

     J[OHplus][O] += rate[Z48]*Y[H3plus];
 
     J[OHplus][H3plus] += rate[Z48]*Y[O];

     J[OHplus][OHplus] += -rate[Z52]*Y[el];

     J[OHplus][el] += -rate[Z52]*Y[OHplus];

     J[OHplus][Hplus] += rate[Z57]*Y[OH]; 

     J[OHplus][OH] += rate[Z57]*Y[Hplus];
 
     /***************** H3Oplus ********/ 

     J[H3Oplus][H3Oplus] += -rate[Z49]*Y[el] 
      -                      rate[Z50]*Y[el]
      -                      rate[Z54]*Y[el];

 
     J[H3Oplus][el] += - rate[Z49]*Y[H3Oplus]
      -                  rate[Z50]*Y[H3Oplus]
      -                  rate[Z54]*Y[H3Oplus];

     J[H3Oplus][H2O] += rate[Z59]*Y[HCOplus]
      +                 rate[Z60]*Y[H3plus];

     J[H3Oplus][HCOplus] += rate[Z59]*Y[H2O]; 

     J[H3Oplus][H3plus] += rate[Z60]*Y[H2O];
 
     /***************** CH *************/ 

     J[CH][CH2] += rate[CRX14];
 
     J[CH][CH3] += rate[CRX16];
 
     J[CH][CH] += -rate[CRX4]
      -            rate[Z73]*Y[H]
      -            rate[Z76]*Y[O]
      -            rate[UV10];


     J[CH][CH3plus] += rate[Z66]*Y[el]
      +                rate[Z67]*Y[el];

     J[CH][el] += rate[Z66]*Y[CH3plus]
      +           rate[Z67]*Y[CH3plus]
      +           rate[Z69]*Y[CH2plus];

     J[CH][CH2plus] += rate[Z69]*Y[el]; 

     J[CH][H] -= rate[Z73]*Y[H]*Y[CH];

     J[CH][O] -= rate[Z76]*Y[CH];
 
     /***************** CHplus *********/ 

     J[CHplus][CHplus] += -rate[CRX3]
      -                    rate[Z62]*Y[H2m]
      -                    rate[Z68]*Y[el]
      -                    rate[Z72]*Y[H];


     J[CHplus][H3plus] += rate[Z61]*Y[C]; 

     J[CHplus][C] += rate[Z61]*Y[H3plus];

     J[CHplus][CH] += rate[UV10];

     J[CHplus][H2m] -= rate[Z62]*Y[CHplus];

     J[CHplus][el] -= rate[Z68]*Y[CHplus];

     J[CHplus][H] -= rate[Z72]*Y[CHplus];
 
     /***************** CH2 ************/ 

     J[CH2][O] += - rate[Z43]*Y[CH2]
      -             rate[Z44]*Y[CH2];

     J[CH2][CH2] += - rate[Z43]*Y[O]
      -               rate[Z44]*Y[O]
      -               rate[CRX13] 
      -               rate[CRX14];
 
     J[CH2][CH4] += rate[CRX7];

     J[CH2][CH3] += rate[CRX15];

     J[CH2][CH3plus] += rate[Z65]*Y[el];

     J[CH2][el] += rate[Z65]*Y[CH3plus];

     /***************** CH2plus ********/ 

     J[CH2plus][CH2] += rate[CRX13];

     J[CH2plus][H2m] += rate[Z62]*Y[CHplus]
      +                 rate[Z63]*Y[Cplus]
      -                 rate[Z64]*Y[CH2plus];

     J[CH2plus][CHplus] += rate[Z62]*Y[H2m];

     J[CH2plus][Cplus] += rate[Z63]*Y[H2m];

     J[CH2plus][CH2plus] -= rate[Z64]*Y[H2m]
      -                     rate[Z69]*Y[el]
      -                     rate[Z70]*Y[el]
      -                     rate[Z71]*Y[el];

     J[CH2plus][el] -= rate[Z69]*Y[CH2plus]
	         -rate[Z70]*Y[CH2plus]
	         -rate[Z71]*Y[CH2plus];
  
     /***************** CH3 ************/ 

     J[CH3][CH3] += -rate[CRX15]
      -              rate[CRX16];

     /***************** CH3plus ********/ 

     J[CH3plus][H2m] += rate[Z64]*Y[CH2plus]; 

     J[CH3plus][CH2plus] += rate[Z64]*Y[H2m]; 

     J[CH3plus][CH3plus] -= rate[Z65]*Y[el]
      -                     rate[Z66]*Y[el]
      -                     rate[Z67]*Y[el]; 

     J[CH3plus][el] -= rate[Z65]*Y[CH3plus]
      -                rate[Z66]*Y[CH3plus]
      -                rate[Z67]*Y[CH3plus];

     /***************** CH4 ************/ 

     J[CH4][CH4] += -rate[CRX7];

     /***************** He *************/ 

     J[He][CO] += rate[Z45]*Y[Heplus];

     J[He][Heplus] += rate[Z45]*Y[CO]
      +               rate[H12]*Y[H]
      +               rate[H13]*Y[H2m]
      +               rate[H14]*Y[H2m];


     J[He][H] += rate[H12]*Y[Heplus]
      +          rate[H17]*Y[HeHplus];

     J[He][H2m] += rate[H13]*Y[Heplus]
      +            rate[H14]*Y[Heplus]
      +            rate[H16]*Y[HeHplus];

     J[He][HeHplus] += rate[H16]*Y[H2m]
      +                rate[H17]*Y[H];

     J[He][H2plus] -= rate[H15]*Y[He];

     J[He][He] -= rate[H15]*Y[H2plus]
      -           rate[H18]*Y[Hplus]
      -           rate[CRX17];

     J[He][Hplus] -= rate[H18]*Y[He];

     /***************** Heplus *********/ 

     J[Heplus][Heplus] += -rate[Z45]*Y[CO]
      -                    rate[H12]*Y[H]
      -                    rate[H13]*Y[H2m]
      -                    rate[H14]*Y[H2m];
  
     J[Heplus][CO] += -rate[Z45]*Y[Heplus];

     J[Heplus][He] += rate[CRX17]; 

     J[Heplus][H] -= rate[H12]*Y[Heplus]; 

     J[Heplus][H2m] -= rate[H13]*Y[Heplus]
	         -rate[H14]*Y[Heplus];	 

     /***************** CO *************/ 

     J[CO][O] += rate[Z43]*Y[CH2] 
      +          rate[Z44]*Y[CH2];

     J[CO][CH2] += rate[Z43]*Y[O] 
      +            rate[Z44]*Y[O]; 

     J[CO][Heplus] += -rate[Z45]*Y[CO];
 
     J[CO][CO] += -rate[Z45]*Y[Heplus]
      -            rate[CRX2]
      -            rate[Z53]*Y[H3plus];

     J[CO][C] += rate[Z55]*Y[OH];

     J[CO][OH] += rate[Z55]*Y[C];

     J[CO][H2O] += rate[Z59]*Y[HCOplus];

     J[CO][HCOplus] += rate[Z59]*Y[H2O]
      +                rate[Z75]*Y[el];

     J[CO][el] += rate[Z75]*Y[HCOplus];

     J[CO][H3plus] -rate[Z53]*Y[CO];

     /***************** COplus *********/ 

     J[COplus][O2] += rate[Z46]*Y[Cplus];
 
     J[COplus][Cplus] += rate[Z46]*Y[O2];
      +                  rate[Z56]*Y[OH];

     J[COplus][OH] += rate[Z56]*Y[Cplus];

     J[COplus][H2m] -= rate[Z74]*Y[COplus]; 

     J[COplus][COplus] -= rate[Z74]*Y[H2m];

     /***************** O2 *************/ 

     J[O2][O2] += -rate[Z46]*Y[Cplus]
      -            rate[CRX8] 
      -            rate[CRX9];

     J[O2][Cplus] += -rate[Z46]*Y[O2];

     /***************** O2plus *********/ 

     J[O2plus][O2] += rate[CRX8];
 
     /***************** H3plus *********/ 

     J[H3plus][O] += -rate[Z48]*Y[H3plus];
 
     J[H3plus][H3plus] += -rate[Z48]*Y[O]
      -                    rate[H10]*Y[el]
      -                    rate[H11]*Y[el]
      -                    rate[Z53]*Y[CO]
      -                    rate[Z60]*Y[H2O]
      -                    rate[Z61]*Y[C];

     J[H3plus][H2plus] += rate[H9]*Y[H2m];

     J[H3plus][H2m] += rate[H9]*Y[H2plus]
      +                rate[H16]*Y[HeHplus];

     J[H3plus][HeHplus] += rate[H16]*Y[H2m];

     J[H3plus][el] += -rate[H10]*Y[H3plus]
      -                rate[H11]*Y[H3plus];

     J[H3plus][CO] += -rate[Z53]*Y[H3plus];

     J[H3plus][H2O] += -rate[Z60]*Y[H3plus];

     J[H3plus][C] += -rate[Z61]*Y[H3plus];

     /***************** OH *************/ 

     J[OH][H3Oplus] += rate[Z49]*Y[el]
      +                rate[Z54]*Y[el];

     J[OH][el] += rate[Z49]*Y[H3Oplus]
      +           rate[Z54]*Y[H3Oplus];

     J[OH][H2O] += rate[CRX5];

     J[OH][OH] += -rate[CRX6]
      -            rate[Z55]*Y[C]
      -            rate[Z56]*Y[Cplus]
      -            rate[Z57]*Y[Hplus];

     J[OH][C] -= rate[Z55]*Y[OH];

     J[OH][Cplus] -= rate[Z56]*Y[OH];

     J[OH][Hplus] -= rate[Z57]*Y[OH]; 
  
     /***************** H2O ************/ 

     J[H2O][H2O] += -rate[CRX5]
      -              rate[Z58]*Y[Cplus]
      -              rate[Z59]*Y[HCOplus]
      -              rate[Z60]*Y[H3plus];


     J[H2O][Cplus] -= rate[Z58]*Y[H2O]; 

     J[H2O][HCOplus] -= rate[Z59]*Y[H2O];

     J[H2O][H3plus] -= rate[Z60]*Y[H2O];

     /***************** H2Oplus ********/ 

     J[H2Oplus][H2Oplus] +=  
      -                      rate[Z51]*Y[el];
 
     J[H2Oplus][el] +=
      -                 rate[Z51]*Y[H2Oplus];

     /***************** HCOplus ********/

     J[HCOplus][H3plus] += rate[Z53]*Y[CO];

     J[HCOplus][CO] += rate[Z53]*Y[H3plus];

     J[HCOplus][Cplus] += rate[Z58]*Y[H2O];

     J[HCOplus][H2O] += rate[Z58]*Y[Cplus]
      -                 rate[Z59]*Y[HCOplus];

     J[HCOplus][H2m] += rate[Z74]*Y[COplus];

     J[HCOplus][COplus] += rate[Z74]*Y[H2m];
 
     J[HCOplus][CH] += rate[Z76]*Y[O];

     J[HCOplus][O] += rate[Z76]*Y[CH];

     J[HCOplus][HCOplus] -= rate[Z59]*Y[H2O]
      -                     rate[Z75]*Y[el]; 

     J[HCOplus][el] -= rate[Z75]*Y[HCOplus];
  
    }
}
