/*********************************************************************************
//
//   THIS TEST IS AN ATTEMPT TO REPRODUCE THE MIDDLE PANELS OF FIGURE 6
//   IN BIALY & STERNBERG MNRAS 450, 4424-4445 (2015)
//
//   Created by Brandon Wiggins
//
//   Make sure this file is in your /grackle_waternet/src/examples folder
//   Type "make" to build and then run
//
//   python plot_bialy.py
//
//   to see the output
//
**********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

extern "C" {
#include <grackle.h>
}

#define mh     1.67262171e-24   
#define kboltz 1.3806504e-16
#define tiny_number 1.e-20

//THIS ROUTINE PRINTS OUT ALL THE ABUNDANCE CRAP
void print_abund(double time, grackle_field_data my_fields, FILE** fp){
  double ndens = (my_fields.HI_density[0] + 
                  my_fields.HII_density[0] +
                  my_fields.HM_density[0] +
              2.0*my_fields.H2I_density[0] +
              2.0*my_fields.H2II_density[0]);
  //ndens = 1.0;
  double mHI = my_fields.HI_density[0];
  double mHII = my_fields.HII_density[0];
  double mHM = my_fields.HM_density[0];
  double mH2I = my_fields.H2I_density[0]*2.0;
  double mH2II = my_fields.H2II_density[0]*2.0;

  fprintf(*fp, "%g \t %g \t %g \t %g \t %g \t %g \t %g \n",time, ndens,
          mHI/ndens,mHII/ndens,mHM/ndens,mH2I/ndens,mH2II/ndens);
}

//THIS ROUTINE INITIALIZES ALL THE CHEMISTRY FIELDS
void reset_fields(grackle_field_data  *my_fields, code_units *my_units, int field_size){

  // set temperature units
  double temperature_units = mh * pow(my_units->a_units * 
                                      my_units->length_units /
                                      my_units->time_units, 2) / kboltz;

  int i;
  for (i = 0;i < field_size;i++) {
    my_fields->density[i] = 1.0e3;
    my_fields->HI_density[i] = my_fields->density[i];
    my_fields->HII_density[i] = tiny_number * my_fields->density[i];
    my_fields->HM_density[i] = tiny_number * my_fields->density[i];
    my_fields->HeI_density[i] = tiny_number * my_fields->density[i];
    my_fields->HeII_density[i] = tiny_number* my_fields->density[i];
    my_fields->HeIII_density[i] = tiny_number* my_fields->density[i];
    my_fields->H2I_density[i] = tiny_number * my_fields->density[i];
    my_fields->H2II_density[i] = tiny_number * my_fields->density[i];
    my_fields->DII_density[i] = tiny_number * my_fields->density[i];
    my_fields->HDI_density[i] = tiny_number * my_fields->density[i];
    my_fields->DI_density[i] = tiny_number * my_fields->density[i];
    my_fields->e_density[i] =  tiny_number * my_fields->density[i];

    my_fields->O_density[i] = tiny_number*my_fields->density[i];
    my_fields->Water_density[i] = tiny_number* my_fields->density[i];
    my_fields->OH_density[i] = tiny_number * my_fields->density[i];
    my_fields->O2_density[i] = tiny_number * my_fields->density[i];
    my_fields->Oplus_density[i] = tiny_number*my_fields->density[i];
    my_fields->OHplus_density[i] = tiny_number*my_fields->density[i];
    my_fields->H2Oplus_density[i] = tiny_number*my_fields->density[i];
    my_fields->H3Oplus_density[i] = tiny_number*my_fields->density[i];
    my_fields->O2plus_density[i] = tiny_number * my_fields->density[i];
    my_fields->Cplus_density[i] = tiny_number * my_fields->density[i];
    my_fields->C_density[i] = tiny_number *my_fields->density[i];
    my_fields->CH_density[i] = tiny_number * my_fields->density[i];
    my_fields->CH2_density[i] = tiny_number * my_fields->density[i];
    my_fields->CH3_density[i] = tiny_number * my_fields->density[i];
    my_fields->CH4_density[i] = tiny_number * my_fields->density[i];

    my_fields->CO_density[i]= tiny_number * my_fields->density[i];
    my_fields->COplus_density[i]= tiny_number * my_fields->density[i];
    my_fields->CO2_density[i]= tiny_number * my_fields->density[i];
    my_fields->CHplus_density[i]= tiny_number * my_fields->density[i];
    my_fields->CH2plus_density[i]= tiny_number * my_fields->density[i];
    my_fields->H3plus_density[i] = tiny_number * my_fields->density[i];
    my_fields->HCOplus_density[i] = tiny_number * my_fields->density[i];
    my_fields->HeHplus_density[i] = tiny_number * my_fields->density[i];
    my_fields->CH3plus_density[i] = tiny_number * my_fields->density[i];
    my_fields->CH4plus_density[i] = tiny_number * my_fields->density[i];
    my_fields->CH5plus_density[i] = tiny_number * my_fields->density[i];
    my_fields->O2Hplus_density[i] = tiny_number * my_fields->density[i];

    // solar metallicity
    //my_fields->metal_density[i] = grackle_data->SolarMetalFractionByMass *
    //      my_fields->density[i];
    my_fields->metal_density[i] = tiny_number * my_fields->density[i];
    //my_fields->metal_density[i] =  0.423501;
    //my_fields->metal_density[i] = grackle_data->SolarMetalFractionByMass *
     // my_fields->density[i]*1.0e-7;

    my_fields->x_velocity[i] = 0.0;
    my_fields->y_velocity[i] = 0.0;
    my_fields->z_velocity[i] = 0.0;

    // initialize internal energy (here 100 K for no reason)
    my_fields->internal_energy[i] = 100. / temperature_units;

    my_fields->volumetric_heating_rate[i] = 0.0;
    my_fields->specific_heating_rate[i] = 0.0;

    my_fields->RT_HI_ionization_rate[i] = 0.0;
    my_fields->RT_HeI_ionization_rate[i] = 0.0;
    my_fields->RT_HeII_ionization_rate[i] = 0.0;
    my_fields->RT_H2_dissociation_rate[i] = 0.0;
    my_fields->RT_heating_rate[i] = 0.0;
  }
}


int main(int argc, char *argv[])
{ 
  /*********************************************************************
  / Initial setup of units and chemistry objects.
  / This should be done at simulation start.
  *********************************************************************/
   
  // Set initial redshift (for internal units).
  double initial_redshift = 0.;

  // Enable output
  grackle_verbose = 1;

  // First, set up the units system.
  // These are conversions from code units to cgs.
  code_units my_units;
  my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units.density_units = 1.67e-24;
  my_units.length_units = 1.0;
  my_units.time_units = 1.0e12;
  my_units.velocity_units = my_units.length_units / my_units.time_units;
  my_units.a_units = 1.0; // units for the expansion factor
  // Set expansion factor to 1 for non-cosmological simulation.
  my_units.a_value = 1. / (1. + initial_redshift) / my_units.a_units;

  // Second, create a chemistry object for parameters.  This needs to be a pointer.
  chemistry_data *my_grackle_data;
  my_grackle_data = new chemistry_data;
  if (set_default_chemistry_parameters(my_grackle_data) == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
    return EXIT_FAILURE;
  }
  // Set parameter values for chemistry.
  // Access the parameter storage with the struct you've created
  // or with the grackle_data pointer declared in grackle.h (see further below).
  grackle_data->use_grackle = 1;            // chemistry on
  grackle_data->with_radiative_cooling = 1; // cooling on
  grackle_data->primordial_chemistry = 3;   // molecular network with H, He, D
  grackle_data->metal_cooling = 1;          // metal cooling on
  grackle_data->UVbackground = 1;           // UV background on
  grackle_data->grackle_data_file = "../../input/CloudyData_UVB=FG2011_shielded.h5"; // data file
  grackle_data->withWater = 1;
  grackle_data->water_only = 1;
  grackle_data->water_rates = 3;  //USE BIALY YOU FOOL
  grackle_data->h2_on_dust=1;

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&my_units) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return EXIT_FAILURE;
  }

  // Create struct for storing grackle field data
  grackle_field_data my_fields;

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  int field_size = 1;
  my_fields.grid_rank = 3;
  my_fields.grid_dimension = new int[3];
  my_fields.grid_start = new int[3];
  my_fields.grid_end = new int[3];
  for (int i = 0;i < 3;i++) {
    my_fields.grid_dimension[i] = 1; // the active dimension not including ghost zones.
    my_fields.grid_start[i] = 0;
    my_fields.grid_end[i] = 0;
  }
  my_fields.grid_dimension[0] = field_size;
  my_fields.grid_end[0] = field_size - 1;
  my_fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

  my_fields.density         = new gr_float[field_size];
  my_fields.internal_energy = new gr_float[field_size];
  my_fields.x_velocity      = new gr_float[field_size];
  my_fields.y_velocity      = new gr_float[field_size];
  my_fields.z_velocity      = new gr_float[field_size];
  // for primordial_chemistry >= 1
  my_fields.HI_density      = new gr_float[field_size];
  my_fields.HII_density     = new gr_float[field_size];
  my_fields.HeI_density     = new gr_float[field_size];
  my_fields.HeII_density    = new gr_float[field_size];
  my_fields.HeIII_density   = new gr_float[field_size];
  my_fields.e_density       = new gr_float[field_size];
  // for primordial_chemistry >= 2
  my_fields.HM_density      = new gr_float[field_size];
  my_fields.H2I_density     = new gr_float[field_size];
  my_fields.H2II_density    = new gr_float[field_size];
  // for primordial_chemistry >= 3
  my_fields.DI_density      = new gr_float[field_size];
  my_fields.DII_density     = new gr_float[field_size];
  my_fields.HDI_density     = new gr_float[field_size];
  // for metal_cooling = 1
  my_fields.metal_density   = new gr_float[field_size];

  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  my_fields.volumetric_heating_rate = new gr_float[field_size];
  // specific heating rate (provide in units [egs s^-1 g^-1]
  my_fields.specific_heating_rate = new gr_float[field_size];

  // radiative transfer ionization / dissociation rate fields (provided in units of [1/s])
  my_fields.RT_HI_ionization_rate = new gr_float[field_size];
  my_fields.RT_HeI_ionization_rate = new gr_float[field_size];
  my_fields.RT_HeII_ionization_rate = new gr_float[field_size];
  my_fields.RT_H2_dissociation_rate = new gr_float[field_size];
  // radiative transfer heating rate (provide in units [erg s^-1 cm^-3])
  my_fields.RT_heating_rate = new gr_float[field_size];

  my_fields.O_density  = new gr_float[field_size];
  my_fields.Water_density  = new gr_float[field_size];
  my_fields.OH_density = new gr_float[field_size];
  my_fields.O2_density = new gr_float[field_size];
  my_fields.Oplus_density = new gr_float[field_size];
  my_fields.OHplus_density = new gr_float[field_size];
  my_fields.H2Oplus_density = new gr_float[field_size];
  my_fields.H3Oplus_density = new gr_float[field_size];
  my_fields.O2plus_density= new gr_float[field_size];
  my_fields.Cplus_density = new gr_float[field_size];
  my_fields.C_density = new gr_float[field_size];
  my_fields.CH_density = new gr_float[field_size];
  my_fields.CH2_density = new gr_float[field_size];
  my_fields.CH3_density = new gr_float[field_size];
  my_fields.CH4_density = new gr_float[field_size];
  my_fields.CO_density = new gr_float[field_size];
  my_fields.COplus_density = new gr_float[field_size];
  my_fields.CO2_density = new gr_float[field_size];
  my_fields.CHplus_density = new gr_float[field_size];
  my_fields.CH2plus_density = new gr_float[field_size];
  my_fields.H3plus_density = new gr_float[field_size];
  my_fields.HCOplus_density = new gr_float[field_size];
  my_fields.HeHplus_density = new gr_float[field_size];
  my_fields.CH3plus_density = new gr_float[field_size];
  my_fields.CH4plus_density =  new gr_float[field_size];
  my_fields.CH5plus_density =  new gr_float[field_size];
  my_fields.O2Hplus_density = new gr_float[field_size];

   //initialize fields
   reset_fields(&my_fields, &my_units, field_size);

   //DEFINE number density THE WAY BIALY LIKES TO
   double ndens = (my_fields.HI_density[0] + 0.5*my_fields.H2I_density[0]);

   //OPEN TEXT FILE THAT WE'LL WRITE ABUNDANCES TO
   FILE *fp;
   fp = fopen("abundances.txt", "w+");
   int j;
   //fprintf(fp, "Metal \t C \t O \t CO \t OH \t Water \t O2 \t CH \n");
   fprintf(fp, "Time \t Hydrogen \t H \t H+ \t H- \t H2 \t H2+ \n");
   fprintf(fp, "%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e\n", my_fields.C_density[0], my_fields.O_density[0], my_fields.CO_density[0], my_fields.OH_density[0], my_fields.Water_density[0], my_fields.O2_density[0], my_fields.CH_density[0]);

   //CREATING A METAL ARRAY THAT VARIES LOG-UNIFORM FROM 10.E-3 TO 1.0
   /*double metal[500]; 
   for (j = 0; j<500; j++){
       double expon = -3.0 + 3.0/500.*(double) j;
       metal[j] = pow(10.0, expon);
   }*/ 

//THIS IS THE TEST HERE: VARY METALLICITY AND FIXED TEMPERATURE, DENSITY AND IONIZATION RATES
//LOOP OVER ALL METALLICITIES (TESTING 500 SUCH METALLICITIES DISTRIBUTED IN LOG SPACE
/*for (k = 0; k<500; k++){
   reset_fields(&my_fields, &my_units); //RESET ALL ABUNDANCES TO INITIAL VALUES
   my_fields.metal_density[0] = grackle_data->SolarMetalFractionByMass *
          my_fields.density[0]*metal[k]; //LET'S FIX THE METALLICITY 
   my_fields.O_density[0] *= metal[k];   //LET'S ALSO FIX O AND C AS PER BIALY
   my_fields.C_density[0] *= metal[k];   
   
    //INTEGRATE NETWORK TO 10 MILLION YEARS
    for (j = 0; j<N_iters; j++){
    double dt = 3.15e7 * 10e7/ (float) N_iters / my_units.time_units;
       if (solve_chemistry(&my_units, &my_fields, dt) == 0) {
       fprintf(stderr, "Error in solve_chemistry.\n");
        return EXIT_FAILURE;}
}*/

   my_fields.metal_density[0] = grackle_data->SolarMetalFractionByMass * my_fields.density[0];
   int N_iters = 1.e3;  //WE'LL GIVE THE NETWORK 1000 STEPS TO GET TO EQUILIBRIUM (THE NETWORK MAY SUBCYCLE THESE UP TO 500 STEPS EACH)

   int k;
   for (k = 0; k < N_iters; k++){
      /*printf("HI: %g\n",my_fields.HI_density[0]);
      printf("HII: %g\n",my_fields.HII_density[0]);
      printf("HM: %g\n",my_fields.HM_density[0]);*/
      double dt = 3.15e7 * 10.e7/ (float) N_iters / my_units.time_units;
      print_abund(k*dt*my_units.time_units,my_fields, &fp);
      if (solve_chemistry(&my_units, &my_fields, dt) == 0) {
         fprintf(stderr, "Error in solve_chemistry.\n");
         return EXIT_FAILURE;
      }
   }

  return EXIT_SUCCESS;
}
