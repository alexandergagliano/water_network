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
    evolve_freefall

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    mass_electron_cgs, \
    sec_per_Myr, \
    cm_per_mpc

tiny_number = 1.e-20

amu_cgs           = 1.660538921e-24  # g
mass_oxygen_cgs = 15.9994*amu_cgs  # g
mass_carbon_cgs = 12.0107*amu_cgs  # g

if __name__=="__main__":
    current_redshift = 0.

    # Set solver parameters
    colors = ["blue", "brown", "magenta", "purple"]
    #metallicity = [1.e-6, 1.e-4, 1.e-2, 1.e-0]
    #logZ = [-6, -4, -2, 0]
    metallicity = [1.e-6]
    logZ = [-6]
    for i in range(len(metallicity)):
        my_chemistry = chemistry_data()
        my_chemistry.use_grackle = 1
        my_chemistry.with_radiative_cooling = 1
        my_chemistry.primordial_chemistry = 3
        my_chemistry.UVbackground = 0
        my_chemistry.self_shielding_method = 0
        my_chemistry.H2_self_shielding = 0
        my_chemistry.water_rates = 1
        my_chemistry.Gamma = 5. / 3.
        my_chemistry.three_body_rate = 1
        my_chemistry.CaseBRecombination = 0
        my_chemistry.cie_cooling = 1
        my_chemistry.h2_optical_depth_approximation = 1
        my_chemistry.cmb_temperature_floor = 1
        my_chemistry.withWater = 1
        my_chemistry.metal_cooling = 1
        my_chemistry.h2_on_dust = 1
        grackle_dir = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__)))))
        my_chemistry.grackle_data_file = os.sep.join(
            [grackle_dir, "input", "cloudy_metals_2008_3D.h5"])

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
        initial_temperature = 300.
        #initial_temperature = 40000. # start the gas at this temperature
        # then begin collapse
        initial_density     = 1.e-1 * mass_hydrogen_cgs # g / cm^3
        # stopping condition
        final_density       = 1.e12 * mass_hydrogen_cgs
         
        rval = my_chemistry.initialize()

        fc = FluidContainer(my_chemistry, 1)
        fc["density"][:] = initial_density / my_chemistry.density_units
        fc["HII"][:] = 1.e-4 * fc["density"]
        fc["HeI"][:] = 8.333e-2 * fc["density"]
        fc["HeII"][:] = tiny_number * fc["density"]
        fc["HeIII"][:] = tiny_number * fc["density"]
        if my_chemistry.primordial_chemistry > 1:
            fc["H2I"][:] = 1.e-6 * fc["density"]
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
            fc["HI"][:] = (1.e0 - fc["HII"][:] - 2.e0 * fc["H2I"][:]) * fc["density"]
            fc["Water_density"][:] = tiny_number * fc["density"]
            fc["O_density"][:] = 3.468e-4 * metallicity[i] * fc["density"]
            fc["OH_density"][:] = tiny_number * fc["density"]
            fc["O2_density"][:] = tiny_number * fc["density"]
            fc["Oplus_density"][:] = tiny_number * fc["density"]
            fc["OHplus_density"][:] = tiny_number * fc["density"]
            fc["H2Oplus_density"][:] = tiny_number * fc["density"]
            fc["H3Oplus_density"][:] = tiny_number * fc["density"]
            fc["O2plus_density"][:] = tiny_number * fc["density"]
            fc["Cplus_density"][:] = 0.927e-4 * metallicity[i] * fc["density"]
            fc["C_density"][:] = tiny_number * fc["density"]
            fc["CH_density"][:] = tiny_number * fc["density"]
            fc["CH2_density"][:] = tiny_number * fc["density"]
            fc["CH3_density"][:] = tiny_number * fc["density"]
            fc["CH4_density"][:] = tiny_number * fc["density"]
            fc["CO_density"][:] = tiny_number * fc["density"]
            fc["COplus_density"][:] = tiny_number  * fc["density"]
            fc["CO2_density"][:] = tiny_number * fc["density"]
            #fc["de"][:] = (fc["HII"][:] + fc["Cplus_density"][:] + fc["OHplus_density"][:] + fc["Oplus_density"][:] + fc["O2plus_density"][:] + fc["COplus_density"][:]) * fc["density"]      
            fc["de"][:] = 1.e-4 * fc["density"]
            fc["HeI"][:] = 8.333e-2 * fc["density"]
            fc["HeII"][:] = tiny_number * fc["density"]
            fc["HeIII"][:] = tiny_number * fc["density"]

        fc["energy"][:] = initial_temperature / \
            fc.chemistry_data.temperature_units / \
            fc.calculate_mean_molecular_weight() / \
        (my_chemistry.Gamma - 1.0)

        fc["x-velocity"][:] = 0.0
        fc["y-velocity"][:] = 0.0
        fc["z-velocity"][:] = 0.0
 
        fc["n_H"][:] = fc.calculate_hydrogen_number_density()
        print(fc["n_H"][:])
        abundanceSum = 0
        for field in fc.density_fields:
            if (field is not "metal") and (field is not "density"):
                abundanceSum += (fc[field][0]/fc["density"][0])
        #sys.exit()

        # timestepping safety factor
        safety_factor = 0.1

        # let the gas cool at constant density from the starting temperature
        # down to a lower temperature to get the species fractions in a
        # reasonable state.
        #cooling_temperature = 300.

        #data0 = evolve_constant_density(
        #    fc, final_temperature=cooling_temperature,
        #    safety_factor=safety_factor)

        #sys.exit()

        # then begin collapse 
        # evolve density and temperature according to free-fall collapse
        data = evolve_freefall(fc, final_density,
                               safety_factor=safety_factor)
      
        #fH2 data!
        #logfH2_df = pd.DataFrame({"log(n_H)" : np.log10( data["n_H"] ), "log(fH2)" : np.log10( data["H2I"] / (data["HDI"] + data["HM"] + data["HI"] + data["HII"] + data["H2I"] + data["H2II"] ) )})
        #logfH2_df.to_csv("logfH2_v_logN_[Z]=%d_2005.csv" % logZ[i], index=None)
 
        #fHD data!
        #logfH2_df = pd.DataFrame({"log(n_H)" : np.log10( data["n_H"] ), "log(fHD)" : np.log10( data["HDI"] / (data["HDI"] + data["DI"] + data["DII"] ) )})
        #logfH2_df.to_csv("logfHD_v_logN_[Z]=%d_2005.csv" % logZ[i], index=None)  

#        # FIRST SET OF PLOTS: H2O, OH, O
#        pyplot.figure(dpi=1000)

        pyplot.plot(np.log10( data["n_H"]), np.log10((data["Water_density"]+data["H2Oplus_density"]) / (data["n_H"])) , label ='H2O', color=colors[0])
        pyplot.plot(np.log10( data["n_H"]), np.log10((data["OH_density"]+data["OHplus_density"]) / (data["n_H"])) , label ='OH', color=colors[1])
        pyplot.plot(np.log10( data["n_H"]), np.log10((data["O_density"] + data["Oplus_density"]) / (data["n_H"])) , label ='O', color=colors[2])


        # OPEN OMUKAI 2005 FOR COMPARISON
        with open('logy(H2O)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                pyplot.plot(Omukai_logn, Omukai_spec, label='H2O, Omukai', color=colors[0], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []

        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('logy(OH)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                pyplot.plot(Omukai_logn, Omukai_spec, label='OH, Omukai', color=colors[1], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []
#        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('logy(O)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                pyplot.plot(Omukai_logn, Omukai_spec, label='O, Omukai', color=colors[2], linestyle=':')
                Omukai_logn = []
                Omukai_spec = []
 
        pyplot.xlabel("log($n_H$ [cm$^{-3}$])")
        pyplot.ylabel("log y(O, OH, H2O)")

        pyplot.ylim(ymin = -15, ymax = -8)
        pyplot.xlim(xmin = -1, xmax = 10)
        pyplot.tight_layout()
        #leg = pyplot.legend(fancybox = True, labelspacing=0.0, loc='best')
        #leg.get_frame().set_alpha(0.5)
 
        if my_chemistry.metal_cooling:
            output = "freefall_metal"
        else:
            output = "freefall"
        pyplot.savefig("%s_Water_abundances1_%d.png" % (output, logZ[i]))
        pyplot.clf()
      
        # SECOND SET OF PLOTS: H2, e, HD 
#        pyplot.figure(dpi=1000)

        pyplot.plot(np.log10( data["n_H"]), np.log10((data["H2I"] + data["H2II"]) / (data["n_H"]*my_chemistry.density_units)) , label ='H2', color=colors[0])
        pyplot.plot(np.log10( data["n_H"]), np.log10(data["de"] / (data["n_H"]*my_chemistry.density_units)) , label ='e', color=colors[1])

        pyplot.plot(np.log10( data["n_H"]), np.log10(data["HDI"] / (data["n_H"]*my_chemistry.density_units)) , label ='HD', color=colors[2])

        
        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('logy(H2)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                pyplot.plot(Omukai_logn, Omukai_spec, label='H2, Omukai', color=colors[0], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []

        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('logy(e)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                pyplot.plot(Omukai_logn, Omukai_spec, label='de, Omukai', color=colors[1], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []

        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('logy(HD)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                pyplot.plot(Omukai_logn, Omukai_spec, label='HD, Omukai', color=colors[2], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []


        pyplot.xlabel("log($n_H$ [cm$^{-3}$])")
        pyplot.ylabel("log y(H2, e, HD)")

        pyplot.ylim(ymin = -12.5, ymax=0)
        pyplot.xlim(xmin = -1, xmax = 10)
        pyplot.tight_layout()
        #leg = pyplot.legend(fancybox = True, labelspacing=0.0, loc='best')
        #leg.get_frame().set_alpha(0.5)

        if my_chemistry.metal_cooling:
            output = "freefall_metal"
        else:
            output = "freefall"

        pyplot.savefig("%s_Water_abundances2_%d.png" % (output, logZ[i]))
        pyplot.clf()

        # THIRD SET OF PLOTS: C+, C, CO
        #pyplot.figure(dpi=1000)

        pyplot.plot(np.log10( data["n_H"]), np.log10((data["Cplus_density"]) / (data["n_H"])) , label ='C+', color=colors[0])
        pyplot.plot(np.log10( data["n_H"]), np.log10(data["C_density"] / (data["n_H"])) , label ='C', color=colors[1])
        pyplot.plot(np.log10( data["n_H"]), np.log10((data["CO_density"] + data["COplus_density"]) / (data["n_H"])) , label ='CO', color=colors[2])

        
        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('logy(C+)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                pyplot.plot(Omukai_logn, Omukai_spec, label='C+, Omukai', color=colors[0], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []

        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('logy(C)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                pyplot.plot(Omukai_logn, Omukai_spec, label='C, Omukai', color=colors[1], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []

        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('logy(CO)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                pyplot.plot(Omukai_logn, Omukai_spec, label='CO, Omukai', color=colors[2], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []


        pyplot.xlabel("log($n_H$ [cm$^{-3}$])")
        pyplot.ylabel("log y(C+, C, CO)")

       # pyplot.ylim(ymin = -20, ymax = -8)
        pyplot.ylim(ymin = -15, ymax = -8)
        pyplot.xlim(xmin = -1, xmax = 10)
        pyplot.tight_layout()
        #leg = pyplot.legend(fancybox = True, labelspacing=0.0, loc='best')
        #leg.get_frame().set_alpha(0.5)
        if my_chemistry.metal_cooling:
            output = "freefall_metal"
        else:
            output = "freefall"

        pyplot.savefig("%s_Water_abundances3_%d.png" % (output, logZ[i]))
        pyplot.clf()

