/***********************************************************************
/
/ Calculate gamma (ratio of specific heats) field
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

int _calculate_temperature(chemistry_data *my_chemistry,
                           chemistry_data_storage *my_rates,
                           code_units *my_units,
                           int grid_rank, int *grid_dimension,
                           int *grid_start, int *grid_end,
                           gr_float *density, gr_float *internal_energy,
                           gr_float *HI_density, gr_float *HII_density, 
                           gr_float *HM_density, gr_float *HeI_density, 
                           gr_float *HeII_density, gr_float *HeIII_density,
                           gr_float *H2I_density, gr_float *H2II_density,
                           gr_float *DI_density, gr_float *DII_density, 
                           gr_float *HDI_density, gr_float *Water_density, 
                           gr_float *O_density, gr_float *OH_density,
                           gr_float *O2_density, gr_float *Oplus_density,
                           gr_float *OHplus_density, gr_float *H2Oplus_density,
                           gr_float *H3Oplus_density, gr_float *O2plus_density,
                           gr_float *Cplus_density, gr_float *C_density,
                           gr_float *CH_density, gr_float *CH2_density,
                           gr_float *CH3_density, gr_float *CH4_density,
                           gr_float *CO_density, gr_float *COplus_density, 
                           gr_float *CO2_density, 
                           gr_float *CHplus_density,
                           gr_float *CH2plus_density,
                           gr_float *H3plus_density, gr_float *e_density, 
                           gr_float *metal_density, gr_float *temperature);

int _calculate_gamma(chemistry_data *my_chemistry,
                     chemistry_data_storage *my_rates,
                     code_units *my_units,
                     int grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *internal_energy,
                     gr_float *HI_density, gr_float *HII_density, 
                     gr_float *HM_density, gr_float *HeI_density, 
                     gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, 
                     gr_float *HDI_density, gr_float *Water_density, 
                     gr_float *O_density, gr_float *OH_density,
                     gr_float *O2_density, gr_float *Oplus_density,
                     gr_float *OHplus_density, gr_float *H2Oplus_density,
                     gr_float *H3Oplus_density, gr_float *O2plus_density,
                     gr_float *Cplus_density, gr_float *C_density,
                     gr_float *CH_density, gr_float *CH2_density,
                     gr_float *CH3_density, gr_float *CH4_density,
                     gr_float *CO_density, gr_float *COplus_density, 
                     gr_float *CO2_density,
                     gr_float *CHplus_density,
                     gr_float *CH2plus_density,
                     gr_float *H3plus_density, gr_float *e_density, 
                     gr_float *metal_density, gr_float *my_gamma)
{

  if (!my_chemistry->use_grackle)
    return SUCCESS;
 
  int i, dim, size = 1;
  for (dim = 0; dim < grid_rank; dim++)
    size *= grid_dimension[dim];
  
  /* If molecular hydrogen is not being used, just use monotonic.
     (this should not really be called, but provide it just in case). */
 
  for (i = 0; i < size; i++) {
    my_gamma[i] = my_chemistry->Gamma;
  }
 
  if (my_chemistry->primordial_chemistry > 1) {

    /* Compute the temperature first. */
 
    if (_calculate_temperature(my_chemistry, my_rates, my_units,
                               grid_rank, grid_dimension,
                               grid_start, grid_end,
                               density, internal_energy,
                               HI_density, HII_density, HM_density,
                               HeI_density, HeII_density, HeIII_density,
                               H2I_density, H2II_density,
                               DI_density, DII_density, HDI_density,
                               Water_density, O_density, OH_density,
                               O2_density, Oplus_density,
                               OHplus_density, H2Oplus_density,
                               H3Oplus_density, O2plus_density,
                               Cplus_density, C_density,
                               CH_density, CH2_density,
                               CH3_density, CH4_density,
                               CO_density, COplus_density, CO2_density,
                               CHplus_density,
                               CH2plus_density,
                               H3plus_density, e_density, 
                               metal_density,
                               my_gamma) == FAIL) {
      fprintf(stderr, "Error in calculate_gamma.\n");
      return FAIL;
    }

    /* Compute Gamma with molecular Hydrogen formula from Omukau \& Nishi
       astro-ph/9811308. */
 
    double x, nH2, number_density, GammaH2Inverse, 
      GammaInverse = 1 / (my_chemistry->Gamma - 1.0);

#   ifdef _OPENMP
#   pragma omp parallel for schedule( runtime ) \
    private( x, nH2, number_density, GammaH2Inverse )
#   endif
    for (i = 0; i < size; i++) {
 
      /* Compute relative number abundence of molecular hydrogen. */
 
      number_density =
        0.25 * (HeI_density[i] + HeII_density[i] + HeIII_density[i]) +
        HI_density[i] + HII_density[i] + HM_density[i] +
        e_density[i];

      /* Adding water chemistry species */ 
      nH2 = 0.5 * (H2I_density[i]  + H2II_density[i]);
 
      /* Only do full computation if there is a reasonable amount of H2.
         The second term in GammaH2Inverse accounts for the vibrational
         degrees of freedom. */
 
      GammaH2Inverse = 0.5*5.0;
      if (nH2 / number_density > 1e-3) {
	x = 6100.0 / my_gamma[i];
	if (x < 10.0)
	  GammaH2Inverse = 0.5*(5 + 2.0 * x*x * exp(x)/POW(exp(x)-1.0,2));
      }
 
      /* Add in H2. */
 
      my_gamma[i] = 1.0 + (nH2 + number_density) /
        (nH2 * GammaH2Inverse + number_density * GammaInverse);
 
    } // end: loop over i
 
  } // end: if (my_chemistry->primordial_chemistry > 1)

  return SUCCESS;
}

int local_calculate_gamma(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          code_units *my_units,
                          grackle_field_data *my_fields,
                          gr_float *my_gamma)
{
  if (_calculate_gamma(my_chemistry, my_rates, my_units,
                       my_fields->grid_rank, my_fields->grid_dimension,
                       my_fields->grid_start, my_fields->grid_end,
                       my_fields->density, my_fields->internal_energy,
                       my_fields->HI_density, my_fields->HII_density, 
                       my_fields->HM_density, my_fields->HeI_density, 
                       my_fields->HeII_density, my_fields->HeIII_density,
                       my_fields->H2I_density, my_fields->H2II_density,
                       my_fields->DI_density, my_fields->DII_density, 
                       my_fields->HDI_density, my_fields->Water_density, 
                       my_fields->O_density, my_fields->OH_density,
                       my_fields->O2_density, my_fields->Oplus_density,
                       my_fields->OHplus_density, my_fields->H2Oplus_density,
                       my_fields->H3Oplus_density, my_fields->O2plus_density,
                       my_fields->Cplus_density, my_fields->C_density,
                       my_fields->CH_density, my_fields->CH2_density,
                       my_fields->CH3_density, my_fields->CH4_density,
                       my_fields->CO_density, my_fields->COplus_density, 
                       my_fields->CO2_density, 
                       my_fields->CHplus_density,
                       my_fields->CH2plus_density, 
                       my_fields->H3plus_density, my_fields->e_density, 
                       my_fields->metal_density, my_gamma) == FAIL) {
    fprintf(stderr, "Error in _calculate_gamma.\n");
    return FAIL;
  }
  return SUCCESS;
}

int calculate_gamma(code_units *my_units,
                    grackle_field_data *my_fields,
                    gr_float *my_gamma)
{
  if (local_calculate_gamma(grackle_data, &grackle_rates, my_units,
                            my_fields, my_gamma) == FAIL) {
    fprintf(stderr, "Error in local_calculate_gamma.\n");
    return FAIL;
  }
  return SUCCESS;
}
