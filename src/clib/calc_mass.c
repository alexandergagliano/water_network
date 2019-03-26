#include <stdlib.h>
#include <stdio.h>
#include <chemistry.h>


double
calculate_mass(double *Y)
{
    
     double weight[26];
     weight[H] = 1.00;
     weight[H2m] = 2.000;
     weight[el] = 0.000000;
     weight[Hplus] = 1.00;
     weight[Hmin] = 1.00;
     weight[D] = 2.00;
     weight[Dplus] = 2.00;
     weight[HD] = 3.00;
     weight[O] = 16.0;
     weight[OH] = 16.0 + 1.00;
     weight[OHplus] = weight[OH];
     weight[H2O] = 2.0*1.00 + 16.0;
     weight[O2] = 32.0;
     weight[Oplus] = 16.0;
     weight[H2Oplus] = weight[H2O];
     weight[H3Oplus] = 3.0 + 16.0;
     weight[O2plus] = 32.0;
     weight[Cplus] = 12.0;
     weight[C] = 12.0;
     weight[CH] = 13.0;
     weight[CH2] = 14.0;
     weight[CH3] = 15.0;
     weight[CH4] = 16.0;
     weight[CO] = 28.0;
     weight[COplus] = 28.0;
     weight[CO2] = 44.0;


     double sum = 0.0; int i;
     for (i = 0; i < 26; i++){
           sum += Y[i]*weight[i];
     }
     return sum;
}

double
calculate_metl_mass(double *Y)
{

     double weight[26];
     weight[H] = 0.00;
     weight[H2m] = 0.000;
     weight[el] = 0.000000;
     weight[Hplus] = 0.00;
     weight[Hmin] = 0.00;
     weight[D] = 0.00;
     weight[Dplus] = 0.00;
     weight[HD] = 0.00;
     weight[O] = 16.0;
     weight[OH] = 16.0 + 1.00;
     weight[OHplus] = weight[OH];
     weight[H2O] = 2.0*1.00 + 16.0;
     weight[O2] = 32.0;
     weight[Oplus] = 16.0;
     weight[H2Oplus] = weight[H2O];
     weight[H3Oplus] = 3.0 + 16.0;
     weight[O2plus] = 32.0;
     weight[Cplus] = 12.0;
     weight[C] = 12.0;
     weight[CH] = 13.0;
     weight[CH2] = 14.0;
     weight[CH3] = 15.0;
     weight[CH4] = 16.0;
     weight[CO] = 28.0;
     weight[COplus] = 28.0;
     weight[CO2] = 44.0;


     double sum = 0.0; int i;
     for (i = 0; i < 26; i++){
           sum += ((Y[i]*weight[i])/(6.022e23));
     }
     return sum;
}



void
convert_to_moles(double *B, double * n_dens)
{

     double amu = 1.660538921e-24; //in grams
     double pro = 1.007316; // in amu
     double neu = 1.008701; // in amu
     double ele = 0.000549; // in amu
     int i;

     // First get total number density
     double summ = 0.0;
        B[el]      /= amu*ele ;
        B[H]       /= amu*(1.0*pro + 1.0*ele) ;
        B[H2m]     /= amu*(2.0*pro + 2.0*ele) ;
        B[Hplus]   /= amu*(1.0*pro) ;
        B[Hmin]    /= amu*(1.0*pro + 2.0*ele) ;
        B[D]       /= amu*(1.0*pro + 1.0*neu + 1.0*ele) ;
        B[Dplus]   /= amu*(1.0*pro + 1.0*neu) ;
        B[HD]      /= amu*(2.0*pro + 1.0*neu + 2.0*ele) ;
        B[O]       /= amu*(8.0*pro + 8.0*neu + 8.0*ele) ;
        B[OH]      /= amu*(9.0*pro + 8.0*neu + 9.0*ele) ;
        B[H2O]     /= amu*(10.0*pro + 8.0*neu + 10.0*ele) ;
        B[O2]      /= amu*(16.0*pro + 16.0*neu + 16.0*ele) ;
        B[Oplus]   /= amu*(8.0*pro + 8.0*neu + 7.0*ele) ;
        B[OHplus]  /= amu*(9.0*pro + 8.0*neu + 8.0*ele) ;
        B[H2Oplus] /= amu*(10.0*pro + 8.0*neu + 9.0*ele) ;
        B[H3Oplus] /= amu*(11.0*pro + 8.0*neu + 10.0*ele) ;
        B[O2plus]  /= amu*(16.0*pro + 16.0*neu + 15.0*ele) ;
        B[Cplus]   /= amu*(6.0*pro + 6.0*neu + 5.0*ele) ;
        B[C]       /= amu*(6.0*pro + 6.0*neu + 6.0*ele) ;
        B[CH]      /= amu*(7.0*pro + 6.0*neu + 7.0*ele) ;
        B[CH2]     /= amu*(8.0*pro + 6.0*neu + 8.0*ele) ;
        B[CH3]     /= amu*(9.0*pro + 6.0*neu + 9.0*ele) ;
        B[CH4]     /= amu*(10.0*pro + 6.0*neu + 10.0*ele) ;
        B[CO]      /= amu*(14.0*pro + 14.0*neu + 14.0*ele) ;
        B[COplus]  /= amu*(14.0*pro + 14.0*neu + 13.0*ele) ;
        B[CO2]     /= amu*(22.0*pro + 22.0*neu + 22.0*ele) ;


     for (i = 0; i < 26; i++){
             summ += B[i];
     }
     *n_dens = summ;


     // Molar fraction is then B/n_dens;
     for (i = 0; i < 26; i++){
        B[i] /= summ;
     }
}


/*
void 
convert_to_dens(double *B, double n_dens)
{

     double amu = 1.660538921e-24; //in grams
     double pro = 1.007276;   // in amu
//     double pro = 1.007316; // in amu
     double neu = 1.008701; // in amu
     double ele = 0.000549; // in amu
     int i;
     
     // First get total number density
     //      double summ = 0.0;
 
        B[el]      *= amu*ele*n_dens;
        B[H]       *= amu*(1.0*pro + 1.0*ele)*n_dens;
        B[H2m]     *= amu*(2.0*pro + 2.0*ele)*n_dens;
        B[Hplus]   *= amu*(1.0*pro)*n_dens;
        B[Hmin]    *= amu*(1.0*pro + 2.0*ele)*n_dens;
        B[D]       *= amu*(1.0*pro + 1.0*neu + 1.0*ele)*n_dens;
        B[Dplus]   *= amu*(1.0*pro + 1.0*neu)*n_dens;
        B[HD]      *= amu*(2.0*pro + 1.0*neu + 2.0*ele)*n_dens;
        B[O]       *= amu*(8.0*pro + 8.0*neu + 8.0*ele)*n_dens; 
        B[OH]      *= amu*(9.0*pro + 8.0*neu + 9.0*ele)*n_dens; 
        B[H2O]     *= amu*(10.0*pro + 8.0*neu + 10.0*ele)*n_dens;
        B[O2]      *= amu*(16.0*pro + 16.0*neu + 16.0*ele)*n_dens;
        B[Oplus]   *= amu*(8.0*pro + 8.0*neu + 7.0*ele)*n_dens;
        B[OHplus]  *= amu*(9.0*pro + 8.0*neu + 8.0*ele) *n_dens;
        B[H2Oplus] *= amu*(10.0*pro + 8.0*neu + 9.0*ele)*n_dens;
        B[H3Oplus] *= amu*(11.0*pro + 8.0*neu + 10.0*ele)*n_dens;
        B[O2plus]  *= amu*(16.0*pro + 16.0*neu + 15.0*ele)*n_dens;
        B[Cplus]   *= amu*(6.0*pro + 6.0*neu + 5.0*ele)*n_dens;
        B[C]       *= amu*(6.0*pro + 6.0*neu + 6.0*ele)*n_dens;
        B[CH]      *= amu*(7.0*pro + 6.0*neu + 7.0*ele)*n_dens;
        B[CH2]     *= amu*(8.0*pro + 6.0*neu + 8.0*ele)*n_dens; 
        B[CH3]     *= amu*(9.0*pro + 6.0*neu + 9.0*ele)*n_dens; 
        B[CH4]     *= amu*(10.0*pro + 6.0*neu + 10.0*ele)*n_dens;
        B[CO]      *= amu*(14.0*pro + 14.0*neu + 14.0*ele)*n_dens;
        B[COplus]  *= amu*(14.0*pro + 14.0*neu + 13.0*ele)*n_dens;
        B[CO2]     *= amu*(22.0*pro + 22.0*neu + 22.0*ele)*n_dens;    

}
*/

void
convert_to_dens(double *B)
{

     double amu = 1.660538921e-24; //in grams
     double pro = 1.007276;   // in amu
     double neu = 1.008701; // in amu
     double ele = 0.000549; // in amu
     int i;

        B[el]      *= amu*ele;
        B[H]       *= amu*(1.0*pro + 1.0*ele);
        B[H2m]     *= amu*(2.0*pro + 2.0*ele);
        B[Hplus]   *= amu*(1.0*pro);
        B[Hmin]    *= amu*(1.0*pro + 2.0*ele);
        B[D]       *= amu*(1.0*pro + 1.0*neu + 1.0*ele);
        B[Dplus]   *= amu*(1.0*pro + 1.0*neu);
        B[HD]      *= amu*(2.0*pro + 1.0*neu + 2.0*ele);
        B[O]       *= amu*(8.0*pro + 8.0*neu + 8.0*ele);
        B[OH]      *= amu*(9.0*pro + 8.0*neu + 9.0*ele);
        B[H2O]     *= amu*(10.0*pro + 8.0*neu + 10.0*ele);
        B[O2]      *= amu*(16.0*pro + 16.0*neu + 16.0*ele);
        B[Oplus]   *= amu*(8.0*pro + 8.0*neu + 7.0*ele);
        B[OHplus]  *= amu*(9.0*pro + 8.0*neu + 8.0*ele) ;
        B[H2Oplus] *= amu*(10.0*pro + 8.0*neu + 9.0*ele);
        B[H3Oplus] *= amu*(11.0*pro + 8.0*neu + 10.0*ele);
        B[O2plus]  *= amu*(16.0*pro + 16.0*neu + 15.0*ele);
        B[Cplus]   *= amu*(6.0*pro + 6.0*neu + 5.0*ele);
        B[C]       *= amu*(6.0*pro + 6.0*neu + 6.0*ele);
        B[CH]      *= amu*(7.0*pro + 6.0*neu + 7.0*ele);
        B[CH2]     *= amu*(8.0*pro + 6.0*neu + 8.0*ele);
        B[CH3]     *= amu*(9.0*pro + 6.0*neu + 9.0*ele); 
        B[CH4]     *= amu*(10.0*pro + 6.0*neu + 10.0*ele);
        B[CO]      *= amu*(14.0*pro + 14.0*neu + 14.0*ele);
        B[COplus]  *= amu*(14.0*pro + 14.0*neu + 13.0*ele);
        B[CO2]     *= amu*(22.0*pro + 22.0*neu + 22.0*ele);

}



void
convert_from_dens(double *B)
{
     
     double amu = 1.660538921e-24; //in grams
     double pro = 1.007276;   // in amu
     double neu = 1.008701; // in amu
     double ele = 0.000549; // in amu
     int i;
     
     // First get total number density
     //      double summ = 0.0;

        B[el]      /= amu*ele;
        B[H]       /= amu*(1.0*pro + 1.0*ele);
        B[H2m]     /= amu*(2.0*pro + 2.0*ele);
        B[Hplus]   /= amu*(1.0*pro);
        B[Hmin]    /= amu*(1.0*pro + 2.0*ele);
        B[D]       /= amu*(1.0*pro + 1.0*neu + 1.0*ele);
        B[Dplus]   /= amu*(1.0*pro + 1.0*neu);
        B[HD]      /= amu*(2.0*pro + 1.0*neu + 2.0*ele);
        B[O]       /= amu*(8.0*pro + 8.0*neu + 8.0*ele);
        B[OH]      /= amu*(9.0*pro + 8.0*neu + 9.0*ele); 
        B[H2O]     /= amu*(10.0*pro + 8.0*neu + 10.0*ele);
        B[O2]      /= amu*(16.0*pro + 16.0*neu + 16.0*ele);
        B[Oplus]   /= amu*(8.0*pro + 8.0*neu + 7.0*ele);
        B[OHplus]  /= amu*(9.0*pro + 8.0*neu + 8.0*ele) ;
        B[H2Oplus] /= amu*(10.0*pro + 8.0*neu + 9.0*ele);
        B[H3Oplus] /= amu*(11.0*pro + 8.0*neu + 10.0*ele);
        B[O2plus]  /= amu*(16.0*pro + 16.0*neu + 15.0*ele);
        B[Cplus]   /= amu*(6.0*pro + 6.0*neu + 5.0*ele);
        B[C]       /= amu*(6.0*pro + 6.0*neu + 6.0*ele);
        B[CH]      /= amu*(7.0*pro + 6.0*neu + 7.0*ele);
        B[CH2]     /= amu*(8.0*pro + 6.0*neu + 8.0*ele);
        B[CH3]     /= amu*(9.0*pro + 6.0*neu + 9.0*ele); 
        B[CH4]     /= amu*(10.0*pro + 6.0*neu + 10.0*ele);
        B[CO]      /= amu*(14.0*pro + 14.0*neu + 14.0*ele);
        B[COplus]  /= amu*(14.0*pro + 14.0*neu + 13.0*ele);
        B[CO2]     /= amu*(22.0*pro + 22.0*neu + 22.0*ele);

}
