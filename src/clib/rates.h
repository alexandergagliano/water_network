#include <stdlib.h>
#include <stdio.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

void get_rates(int water_rate, double *rate, double temp, double n, double t_dust, double Z, double UV, code_units *my_units, int ispecies, double Y[], int H2_shield, double crsHI, double k24);
//void get_rates(double *rate, double temp, double n, double t_dust, double Z, double UV, code_units *my_units);
