/***********************************************************************
/
/ Calculate pressure field
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;

int _calculate_pressure(chemistry_data *my_chemistry,
                        chemistry_data_storage *my_rates,
                        code_units *my_units,
                        int grid_rank, int *grid_dimension,
                        int *grid_start, int *grid_end,
                        gr_float *density, gr_float *internal_energy,
                        gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                        gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                        gr_float *H2I_density, gr_float *H2II_density,
                        gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                        gr_float *Water_density, gr_float *O_density, gr_float *OH_density,
                        gr_float *O2_density, gr_float *Oplus_density, gr_float *OHplus_density,
                        gr_float *H2Oplus_density, gr_float *H3Oplus_density, gr_float *O2plus_density,
                        gr_float *Cplus_density, gr_float *C_density, gr_float *CH_density,
                        gr_float *CH2_density, gr_float *CH3_density, gr_float *CH4_density,
                        gr_float *CO_density, gr_float *COplus_density, gr_float *CO2_density,
                        gr_float *e_density, gr_float *metal_density,
                        gr_float *pressure)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  double tiny_number = 1.e-20;
  int i, dim, size = 1;
  for (dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];

# ifdef _OPENMP
# pragma omp parallel for schedule( runtime ) private( i )
# endif
  for (i = 0; i < size; i++) {
 
    pressure[i] = (my_chemistry->Gamma - 1.0) * density[i] * internal_energy[i];
 
    if (pressure[i] < tiny_number)
      pressure[i] = tiny_number;
  }
  
  /* Correct for Gamma from H2. */
 
  if (my_chemistry->primordial_chemistry > 1) {
 
    /* Calculate temperature units. */

    double temperature_units =  mh * POW(my_units->velocity_units, 2) / kboltz;

    double number_density, nH2, GammaH2Inverse,
      GammaInverse = 1.0/(my_chemistry->Gamma-1.0), x, Gamma1, temp;
  
#   ifdef _OPENMP
#   pragma omp parallel for schedule( runtime ) \
    private( i, number_density, nH2, GammaH2Inverse, x, Gamma1, temp )
#   endif
    for (i = 0; i < size; i++) {
 
      number_density =
        0.25 * (HeI_density[i] + HeII_density[i] + HeIII_density[i]) +
        HI_density[i] + HII_density[i] + HM_density[i] +
        e_density[i];

      /* Adding water chemistry species */
      //maybe remove later??
      /*
      if (my_chemistry->withWater) {
          number_density += Water_density[i]/18. + O_density[i]/16. + OH_density[i]/17. +
          O2_density[i]/32. + Oplus_density[i]/16. + OHplus_density[i]/17. +
          H2Oplus_density[i]/18. + H3Oplus_density[i]/19. + O2plus_density[i]/32. +
          Cplus_density[i]/12. + C_density[i]/12. + CH_density[i]/13. + CH2_density[i]/14. +
          CH3_density[i]/15. + CH4_density[i]/16. + CO_density[i]/28. + COplus_density[i]/28. +
          CO2_density[i]/44.;
      }*/

      nH2 = 0.5 * (H2I_density[i] + H2II_density[i]);
 
      /* First, approximate temperature. */
 
      if (number_density == 0)
	number_density = tiny_number;
      temp = max(temperature_units * pressure[i] / (number_density + nH2), 1);
 
      /* Only do full computation if there is a reasonable amount of H2.
	 The second term in GammaH2Inverse accounts for the vibrational
	 degrees of freedom. */
 
      GammaH2Inverse = 0.5*5.0;
      if (nH2 / number_density > 1e-3) {
        x = 6100.0 / temp;
	if (x < 10.0)
	  GammaH2Inverse = 0.5*(5 + 2.0 * x*x * exp(x)/POW(exp(x)-1.0,2));
      }
 
      Gamma1 = 1.0 + (nH2 + number_density) /
	             (nH2 * GammaH2Inverse + number_density * GammaInverse);
	
      /* Correct pressure with improved Gamma. */
 
      pressure[i] *= (Gamma1 - 1.0) / (my_chemistry->Gamma - 1.0);
 
    } // end: loop over i
 
  } // end: if (my_chemistry->primordial_chemistry > 1)
 
  return SUCCESS;
}

int local_calculate_pressure(chemistry_data *my_chemistry,
                             chemistry_data_storage *my_rates,
                             code_units *my_units,
                             grackle_field_data *my_fields,
                             gr_float *pressure)
{
  if (_calculate_pressure(my_chemistry, my_rates, my_units,
                          my_fields->grid_rank, my_fields->grid_dimension,
                          my_fields->grid_start, my_fields->grid_end,
                          my_fields->density, my_fields->internal_energy,
                          my_fields->HI_density, my_fields->HII_density,
                          my_fields->HM_density,
                          my_fields->HeI_density, my_fields->HeII_density,
                          my_fields->HeIII_density,
                          my_fields->H2I_density, my_fields->H2II_density,
                          my_fields->DI_density, my_fields->DII_density,
                          my_fields->HDI_density,
                          my_fields->Water_density, my_fields->O_density, 
                          my_fields->OH_density, my_fields->O2_density, 
                          my_fields->Oplus_density, my_fields->OHplus_density,
                          my_fields->H2Oplus_density, my_fields->H3Oplus_density, 
                          my_fields->O2plus_density, my_fields->Cplus_density, 
                          my_fields->C_density, my_fields->CH_density,
                          my_fields->CH2_density, my_fields->CH3_density,  
                          my_fields->CH4_density, my_fields->CO_density, 
                          my_fields->COplus_density, my_fields->CO2_density,
                          my_fields->e_density, my_fields->metal_density,
                          pressure) == FAIL) {
    fprintf(stderr, "Error in _calculate_pressure.\n");
    return FAIL;
  }
  return SUCCESS;
}

int calculate_pressure(code_units *my_units,
                       grackle_field_data *my_fields,
                       gr_float *pressure)
{
  if (local_calculate_pressure(grackle_data, &grackle_rates, my_units,
                               my_fields, pressure) == FAIL) {
    fprintf(stderr, "Error in local_calculate_pressure.\n");
    return FAIL;
  }
  return SUCCESS;
}
