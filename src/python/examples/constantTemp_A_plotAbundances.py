########################################################################
#
# Free-fall example script
#
#
# Copyright (c) 2013-2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from matplotlib import pyplot
import os
import yt
import numpy as np
import pandas as pd
import sys

from pygrackle import \
    chemistry_data, \
    FluidContainer, \
    evolve_constant_density, \
    evolve_freefall, \
    evolve_freefall_metal

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    mass_electron_cgs, \
    sec_per_Myr, \
    cm_per_mpc

tiny_number = 1.e-9

amu_cgs           = 1.660538921e-24  # g
mass_oxygen_cgs = 15.9994*amu_cgs  # g
mass_carbon_cgs = 12.0107*amu_cgs  # g

if __name__=="__main__":
    current_redshift = 0.

    # Set solver parameters
    colors = ["blue", "brown", "magenta", "purple"]
    #metallicity = [1.e-6, 1.e-4, 1.e-2, 1.e-0]
    #logZ = [-6, -4, -2, 0]
    metallicity = [1.e-3]
    final_metallicity = 1.0
    for i in range(len(metallicity)):
        my_chemistry = chemistry_data()
        my_chemistry.use_grackle = 1
        my_chemistry.with_radiative_cooling = 1
        my_chemistry.primordial_chemistry = 3
        my_chemistry.UVbackground = 1
        my_chemistry.self_shielding_method = 1
        my_chemistry.H2_self_shielding = 3
        my_chemistry.Gamma = 5. / 3.
        my_chemistry.three_body_rate = 1
        my_chemistry.CaseBRecombination = 0
        my_chemistry.cie_cooling = 0
        my_chemistry.h2_optical_depth_approximation = 1
        my_chemistry.water_rates = 3
        my_chemistry.cmb_temperature_floor = 1
        # TURNING WATER OFF FOR NOW!!
        my_chemistry.withWater = 1
        #metal cooling needs to be on for this test!
        my_chemistry.metal_cooling = 1
        grackle_dir = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__)))))
        my_chemistry.grackle_data_file = os.sep.join(
            [grackle_dir, "input", "CloudyData_UVB=FG2011_shielded.h5"])
        my_chemistry.h2_on_dust = 1
        # Set units
        my_chemistry.comoving_coordinates = 0 # proper units
        my_chemistry.a_units = 1.0
        my_chemistry.a_value = 1. / (1. + current_redshift) / \
            my_chemistry.a_units
        my_chemistry.density_units  = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
        my_chemistry.length_units   = cm_per_mpc        # 1 Mpc in cm
        my_chemistry.time_units     = sec_per_Myr       # 1 Myr in s
        my_chemistry.velocity_units = my_chemistry.a_units * \
            (my_chemistry.length_units / my_chemistry.a_value) / \
            my_chemistry.time_units

        # set initial density and temperature
        initial_temperature = 100.
        initial_density     = 1.e1 * mass_hydrogen_cgs # g / cm^3
         
        rval = my_chemistry.initialize()

        fc = FluidContainer(my_chemistry, 1)
        fc["density"][:] = initial_density / my_chemistry.density_units
        fc["HII"][:] = tiny_number * fc["density"]
        fc["HI"][:] = 6.5e-1 * fc["density"]
        fc["HeI"][:] = tiny_number * fc["density"]
        fc["HeII"][:] = tiny_number * fc["density"]
        fc["HeIII"][:] = tiny_number * fc["density"]
        if my_chemistry.primordial_chemistry > 1:
            fc["H2I"][:] = 2.5e-1 * fc["density"]
            fc["H2II"][:] = tiny_number * fc["density"]
            fc["HM"][:] = tiny_number * fc["density"]
        if my_chemistry.primordial_chemistry > 2:
            fc["HDI"][:] = 1.e-9 * fc["density"]
            fc["DI"][:] = (3.e-5 - fc["HDI"][:]) * fc["density"]
            fc["DII"][:] = tiny_number * fc["density"]
        if my_chemistry.metal_cooling == 1:
            fc["metal"][:] = metallicity[i] * my_chemistry.SolarMetalFractionByMass * \
                fc["density"]
        if my_chemistry.withWater == 1:
            fc["HeI"][:] = 0.1 * fc["density"]
#            fc["Water_density"][:] = 0.0 * fc["density"]
#            fc["O_density"][:] = 4.9e-3 * metallicity[i] * fc["density"]
#            fc["OH_density"][:] = 0.0 * fc["density"]
#            fc["O2_density"][:] = 0.0 * fc["density"]
#            fc["Oplus_density"][:] = tiny_number * fc["density"]
#            fc["OHplus_density"][:] = tiny_number * fc["density"]
#            fc["H2Oplus_density"][:] = tiny_number * fc["density"]
#            fc["H3Oplus_density"][:] = tiny_number * fc["density"]
#            fc["O2plus_density"][:] = tiny_number * fc["density"]
#            fc["Cplus_density"][:] = tiny_number * fc["density"]
#            fc["C_density"][:] = 2.9e-3 * metallicity[i] * fc["density"]
#            fc["CH_density"][:] = tiny_number * fc["density"]
#            fc["CH2_density"][:] = tiny_number * fc["density"]
#            fc["CH3_density"][:] = tiny_number * fc["density"]
#            fc["CH4_density"][:] = tiny_number * fc["density"]
#            fc["CO_density"][:] = 9.e-9 * fc["density"]
#            fc["COplus_density"][:] = tiny_number  * fc["density"]
#            fc["CO2_density"][:] = tiny_number * fc["density"]

            fc["Water_density"][:] = tiny_number * fc["density"]
            fc["O_density"][:] = 4.9e-4 * metallicity[i] * fc["density"]
#            fc["OH_density"][:] = 3.e-6 * fc["density"]
#            fc["O2_density"][:] = 1.e-7 * fc["density"]
            fc["OH_density"][:] = tiny_number * fc["density"]
            fc["O2_density"][:] = tiny_number * fc["density"]
            fc["Oplus_density"][:] = tiny_number * fc["density"]
            fc["OHplus_density"][:] = tiny_number * fc["density"]
            fc["H2Oplus_density"][:] = tiny_number * fc["density"]
            fc["H3Oplus_density"][:] = tiny_number * fc["density"]
            fc["O2plus_density"][:] = tiny_number * fc["density"]
            fc["Cplus_density"][:] = tiny_number * fc["density"]
            fc["C_density"][:] = 2.9e-4 * metallicity[i] * fc["density"]
#            fc["CH_density"][:] = 1.e-11 * fc["density"]
            fc["CH_density"][:] = tiny_number * fc["density"]
            fc["CH2_density"][:] = tiny_number * fc["density"]
            fc["CH3_density"][:] = tiny_number * fc["density"]
            fc["CH4_density"][:] = tiny_number * fc["density"]
#            fc["CO_density"][:] = 2.e-9 * fc["density"]
            fc["CO_density"][:] = tiny_number * fc["density"]
            fc["COplus_density"][:] = tiny_number  * fc["density"]
            fc["CO2_density"][:] = tiny_number * fc["density"]

            #fc["de"][:] = 8.e-5 * fc["density"]
            #fc["HeI"][:] = tiny_number * fc["density"]
            #fc["HeII"][:] = tiny_number * fc["density"]
            #fc["HeIII"][:] = tiny_number * fc["density"]

        fc["energy"][:] = initial_temperature / \
            fc.chemistry_data.temperature_units / \
            fc.calculate_mean_molecular_weight() / \
        (my_chemistry.Gamma - 1.0)

        fc["x-velocity"][:] = 0.0
        fc["y-velocity"][:] = 0.0
        fc["z-velocity"][:] = 0.0
 
        #fc["n_H"][:] = fc.calculate_hydrogen_number_density()
        #abundanceSum = 0
        #for field in fc.density_fields:
        #    if (field is not "metal") and (field is not "density"):
        #        abundanceSum += (fc[field][0]/fc["density"][0])
        #print("Abundance sum = %.2e\n" % abundanceSum)
        #sys.exit()

        # timestepping safety factor
        safety_factor = 0.01
        # then begin collapse 
        # evolve density and temperature according to free-fall collapse
        data = evolve_freefall_metal(fc, final_metallicity,
                               safety_factor=safety_factor)
      
        n = data["HII"].in_units("cm**-3", equivalence="number_density") + data["HI"].in_units("cm**-3", equivalence="number_density") + 2*data["H2I"].in_units("cm**-3", equivalence="number_density") + 2*data["H2II"].in_units("cm**-3", equivalence="number_density")
        metallicity = (data["metal"] / my_chemistry.SolarMetalFractionByMass / data["density"])
 
        #pyplot.figure(dpi=500)

############################################################
#                                                          #
#   PLOT ONE: METAL MASS FRACTIONS VS METALLICITY          #
#                                                          #
############################################################
        if my_chemistry.withWater == 1:
        ## OPEN FOR COMPARISON
                with open('MassFraction_vs_Metallicity2_C.csv') as bialyData:
                        lines = bialyData.readlines()
                        bialy_Z = []
                        bialy_C = []
                        for line in lines[6:]:
                                p = line.split(',')
                                bialy_Z.append(float(p[0]))
                                bialy_C.append(float(p[1]))
                        pyplot.loglog(bialy_Z, bialy_C, label='C, Bialy', color='black', linestyle=':',lw=2)
                        bialy_Z = []

                with open('MassFraction_vs_Metallicity2_CH.csv') as bialyData:
                        lines = bialyData.readlines()
                        bialy_Z = []
                        bialy_CH = []
                        for line in lines[6:]:
                                p = line.split(',')
                                bialy_Z.append(float(p[0]))
                                bialy_CH.append(float(p[1]))
                        pyplot.loglog(bialy_Z, bialy_CH, label='CH, Bialy', color='teal', linestyle=':',lw=2)
                        bialy_Z = []

                with open('MassFraction_vs_Metallicity2_O.csv') as bialyData:
                        lines = bialyData.readlines()
                        bialy_Z = []
                        bialy_O = []
                        for line in lines[6:]:
                                p = line.split(',')
                                bialy_Z.append(float(p[0]))
                                bialy_O.append(float(p[1]))
                        pyplot.loglog(bialy_Z, bialy_O, label='O, Bialy', color='maroon', linestyle=':',lw=2)
                        bialy_Z = []

                with open('MassFraction_vs_Metallicity2_H2O.csv') as bialyData:
                        lines = bialyData.readlines()
                        bialy_Z = []
                        bialy_H2O = []
                        for line in lines[6:]:
                                p = line.split(',')
                                bialy_Z.append(float(p[0]))
                                bialy_H2O.append(float(p[1]))
                        pyplot.loglog(bialy_Z, bialy_H2O, label='H2O, Bialy', color='red', linestyle=':',lw=2)
                        bialy_Z = []

                with open('MassFraction_vs_Metallicity2_O2.csv') as bialyData:
                        lines = bialyData.readlines()
                        bialy_Z = []
                        bialy_O2 = []
                        for line in lines[6:]:
                                p = line.split(',')
                                bialy_Z.append(float(p[0]))
                                bialy_O2.append(float(p[1]))
                        pyplot.loglog(bialy_Z, bialy_O2, label='O2, Bialy', color='magenta', linestyle=':',lw=2)
                        bialy_Z = []

                with open('MassFraction_vs_Metallicity2_CO.csv') as bialyData:
                        lines = bialyData.readlines()
                        bialy_Z = []
                        bialy_CO = []
                        for line in lines[6:]:
                                p = line.split(',')
                                bialy_Z.append(float(p[0]))
                                bialy_CO.append(float(p[1]))
                        pyplot.loglog(bialy_Z, bialy_CO, label='CO, Bialy', color='blue', linestyle=':',lw=2)
                        bialy_Z = []

                with open('MassFraction_vs_Metallicity2_OH.csv') as bialyData:
                        lines = bialyData.readlines()
                        bialy_Z = []
                        bialy_OH = []
                        for line in lines[6:]:
                                p = line.split(',')
                                bialy_Z.append(float(p[0]))
                                bialy_OH.append(float(p[1]))
                        pyplot.loglog(bialy_Z, bialy_OH, label='OH, Bialy', color='goldenrod', linestyle=':',lw=2)

                pyplot.loglog(metallicity, (data["O_density"]/n), label ='O', color='maroon', lw=2)
                pyplot.loglog(metallicity, (data["C_density"]/n), label ='C', color='black', lw=2)
                pyplot.loglog(metallicity, (data["CO_density"]/n), label ='CO', color='blue', lw=2)
                pyplot.loglog(metallicity, (data["Water_density"]/n), label ='H2O', color='red', lw=2)
                pyplot.loglog(metallicity, (data["OH_density"]/n), label ='OH', color='goldenrod', lw=2)
                pyplot.loglog(metallicity, (data["O2_density"]/n), label ='O2', color='magenta', lw=2)
                pyplot.loglog(metallicity, (data["CH_density"]/n), label ='CH', color='teal', lw=2)

                pyplot.xlabel("Z'",fontsize=22)
                pyplot.ylabel("Xi",fontsize=22)

                #leg = pyplot.legend(fancybox = True, labelspacing=0.0, loc='best')
                #leg.get_frame().set_alpha(0.5)
                #pyplot.ylim(ymin=1.e-10,ymax=1.e-3)
                #pyplot.xlim(xmin=1.e-3,xmax=1.e0)
                pyplot.tight_layout()
                pyplot.savefig("metal_evolve_abundances_metals.png")

                pyplot.clf()


############################################################
#                                                          #
#             PLOT TWO: H2, HI VS METALLICITY              #
#                                                          #
############################################################

        #pyplot.figure(dpi=1000)
        ## OPEN FOR COMPARISON
        with open('MassFraction_vs_Metallicity3_H.csv') as bialyData:
                lines = bialyData.readlines()
                bialy_Z = []
                bialy_H = []
                for line in lines[6:]:
                        p = line.split(',')
                        bialy_Z.append(float(p[0]))
                        bialy_H.append(float(p[1]))
                pyplot.loglog(bialy_Z, bialy_H, label='H, Bialy', color='coral', linestyle=':', lw=2)
                bialy_Z = []

        with open('MassFraction_vs_Metallicity3_H2.csv') as bialyData:
                lines = bialyData.readlines()
                bialy_Z = []
                bialy_H2 = []
                for line in lines[6:]:
                        p = line.split(',')
                        bialy_Z.append(float(p[0]))
                        bialy_H2.append(float(p[1]))
                pyplot.loglog(bialy_Z, bialy_H2, label='H2, Bialy', color='black', linestyle=':', lw=2)
                bialy_Z = []

        H2 = (data["H2I"].in_units("cm**-3", equivalence="number_density") + data["H2II"].in_units("cm**-3", equivalence="number_density"))/n
  
        H = (data["HI"].in_units("cm**-3", equivalence="number_density") + data["HII"].in_units("cm**-3", equivalence="number_density"))/n
        pyplot.loglog(metallicity, H2, label ='H2', color='black', lw=2)
        pyplot.loglog(metallicity, H, label ='H', color='coral', lw=2)

        pyplot.xlabel("Z'",fontsize=22)
        pyplot.ylabel("Xi",fontsize=22)

        #leg = pyplot.legend(fancybox = True, labelspacing=0.0, loc='best')
        #leg.get_frame().set_alpha(0.5)
        pyplot.ylim(ymin=1.e-2,ymax=1.e0)
        pyplot.xlim(xmin=1.e-3,xmax=1.e0)
        pyplot.tight_layout()
        #pyplot.xlim(xmin=1.e-3,xmax=1.e0)
        pyplot.savefig("metal_evolve_abundances_H2.png")

        pyplot.clf()


############################################################
#                                                          #
#       PLOT THREE: IONS,ELECTRONS VS METALLICITY          #
#                                                          #
############################################################

        #pyplot.figure(dpi=1000)

        with open('MassFraction_vs_Metallicity_electrons.csv') as bialyData:
                lines = bialyData.readlines()
                bialy_Z = []
                bialy_electrons = []
                for line in lines[6:]:
                        p = line.split(',')
                        bialy_Z.append(float(p[0]))
                        bialy_electrons.append(float(p[1]))
                pyplot.loglog(bialy_Z, bialy_electrons, label='e, Bialy', color='magenta', linestyle=':', lw=2)
                bialy_Z = []
        if my_chemistry.withWater == 1:
                with open('MassFraction_vs_Metallicity_metalions.csv') as bialyData:
                        lines = bialyData.readlines()
                        bialy_Z = []
                        bialy_metal = []
                        for line in lines[6:]:
                             p = line.split(',')
                             bialy_Z.append(float(p[0]))
                             bialy_metal.append(float(p[1]))
                        pyplot.loglog(bialy_Z, bialy_metal, label='metal ions, Bialy', color='orange', linestyle=':', lw=2)
                        bialy_Z = []

        with open('MassFraction_vs_Metallicity_H+.csv') as bialyData:
                lines = bialyData.readlines()
                bialy_Z = []
                bialy_Hplus = []
                for line in lines[6:]:
                        p = line.split(',')
                        bialy_Z.append(float(p[0]))
                        bialy_Hplus.append(float(p[1]))
                pyplot.loglog(bialy_Z, bialy_Hplus, label='H+, Bialy', color='indigo', linestyle=':', lw=2)
                bialy_Z = []

        if my_chemistry.withWater == 1:
                n_metalIons = data["COplus_density"] + data["Oplus_density"] + data["OHplus_density"] + data["H2Oplus_density"] + data["H3Oplus_density"] + data["O2plus_density"] + data["Cplus_density"]

        pyplot.loglog(metallicity, (data["HII"].in_units("cm**-3", equivalence="number_density")/n), label =r'H$^+$', color='indigo', lw=2)
        pyplot.loglog(metallicity, (data["de"].in_units("cm**-3", equivalence="number_density")/n), label ='e', color='magenta', lw=2)

        if my_chemistry.withWater == 1:
                pyplot.loglog(metallicity, (n_metalIons/n), label ='metal ions', color='orange', lw=2)

        pyplot.xlabel("Z'",fontsize=22)
        pyplot.ylabel("Xi",fontsize=22)

        leg = pyplot.legend(fancybox = True, labelspacing=0.0, loc='best')
        leg.get_frame().set_alpha(0.5)
        pyplot.xlim(xmin=1.e-3)
       # pyplot.ylim(ymin=1.e-8,ymax=1.e-4)
       # pyplot.xlim(xmin=1.e-3,xmax=1.e0)
        pyplot.tight_layout()
        pyplot.savefig("metal_evolve_abundances_H+.png")




