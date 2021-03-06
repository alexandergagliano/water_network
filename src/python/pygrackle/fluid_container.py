########################################################################
#
# Grackle python fluid container
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import numpy as np

from pygrackle.grackle_wrapper import \
    calculate_cooling_time, \
    calculate_gamma, \
    calculate_pressure, \
    calculate_temperature, \
    solve_chemistry

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs

_base_fluids = ["density", "metal"]
_nd_fields   = ["energy",
                "x-velocity", "y-velocity", "z-velocity",
                "temperature", "pressure",
                "gamma", "cooling_time"]
_fluid_names = {}
_fluid_names[0] = _base_fluids
_fluid_names[1] = _fluid_names[0] + \
  ["HI", "HII", "HeI", "HeII", "HeIII", "de"]
_fluid_names[2] = _fluid_names[1] + \
  ["H2I", "H2II", "HM"]
_fluid_names[3] = _fluid_names[2] + \
  ["DI", "DII", "HDI"]

_water_names = \
  ["Water_density", "O_density", "OH_density", "O2_density", "Oplus_density", "OHplus_density", "H2Oplus_density", "H3Oplus_density", "O2plus_density", "Cplus_density", "C_density", "CH_density", "CH2_density", "CH3_density", "CH4_density", "CO_density", "COplus_density", "CO2_density"]

_bialy_names = \
["CHplus_density", "CH2plus_density", "H3plus_density", "HCOplus_density", "HeHplus_density", "CH3plus_density", "CH4plus_density", "CH5plus_density", "O2Hplus_density"]

_rad_trans_names = ['RT_heating_rate', 'RT_HI_ionization_rate',
                    'RT_HeI_ionization_rate', 'RT_HeII_ionization_rate',
                    'RT_H2_dissociation_rate']

try:
    xrange
except NameError:
    xrange = range

class FluidContainer(dict):
    def __init__(self, chemistry_data, n_vals, dtype="float64",
                 itype="int64"):
        super(FluidContainer, self).__init__()
        self.dtype = dtype
        self.chemistry_data = chemistry_data
        self.n_vals = n_vals
        for fluid in _fluid_names[self.chemistry_data.primordial_chemistry] + \
            _nd_fields:
            self._setup_fluid(fluid)
        if self.chemistry_data.use_radiative_transfer:
            for fluid in _rad_trans_names:
                self._setup_fluid(fluid)
        if self.chemistry_data.withWater:
            if (self.chemistry_data.water_rates == 3):
                for fluid in (_water_names + _bialy_names):
                   self._setup_fluid(fluid)
            else: 
                for fluid in _water_names:
                   self._setup_fluid(fluid)

#Set up number density of hydrogen as well!
        self._setup_fluid("n_H")
    def _setup_fluid(self, fluid_name):
        self[fluid_name] = np.zeros(self.n_vals, self.dtype)

    @property
    def cooling_units(self):
        tbase1 = self.chemistry_data.time_units
        if self.chemistry_data.comoving_coordinates:
            xbase1 = self.chemistry_data.length_units / \
                (self.chemistry_data.a_value * self.chemistry_data.a_units)
            dbase1 = self.chemistry_data.density_units * \
                (self.chemistry_data.a_value * self.chemistry_data.a_units)**3
        else:
            xbase1 = self.chemistry_data.length_units / \
                self.chemistry_data.a_units
            dbase1 = self.chemistry_data.density_units * \
                self.chemistry_data.a_units**3

        coolunit = (self.chemistry_data.a_units**5 * xbase1**2 *
                    mass_hydrogen_cgs**2) / (tbase1**3 * dbase1)
        return coolunit

    @property
    def density_fields(self):
            return _fluid_names[self.chemistry_data.primordial_chemistry]

    @property
    def n_density_fields(self):
        if self.chemistry_data.withWater:
            if (self.chemistry_data.water_rates == 3):
               return _water_names + _bialy_names
            else:
               return _water_names
        else:
            return []

    def calculate_hydrogen_number_density(self):
        my_chemistry = self.chemistry_data
        if my_chemistry.primordial_chemistry == 0:
            return my_chemistry.HydrogenFractionByMass * \
              self["density"] * my_chemistry.density_units / mass_hydrogen_cgs
        nH = self["HI"] + self["HII"]
        if my_chemistry.primordial_chemistry > 1:
            nH += self["HM"] + self["H2I"] + self["H2II"]
        if my_chemistry.primordial_chemistry > 2:
            nH += self["HDI"] / 2.
        return nH * my_chemistry.density_units / mass_hydrogen_cgs
 
    def calculate_metallicity(self):
        my_chemistry = self.chemistry_data
        if my_chemistry.withWater == 0:
            return 0
        metal_density = 16.0*self["O_density"] + 17.0*self["OH_density"] + 18.0*self["Water_density"] + \
                32.0*self["O2_density"] + 16.0*self["Oplus_density"] + 17.0*self["OHplus_density"] + \
                18.0*self["H2Oplus_density"] + 19.0*self["H3Oplus_density"] + 32.0*self["O2plus_density"] + \
                12.0*self["Cplus_density"] + 12.0*self["C_density"] + 13.0*self["CH_density"] + \
                14.0*self["CH2_density"] + 15.0*self["CH3_density"] + 16.0*self["CH4_density"] + \
                28.0*self["CO_density"] + 28.0*self["COplus_density"] + 44.0*self["CO2_density"] 
        if my_chemistry.water_rates == 3:
                metal_density += 13.0*self["CHplus_density"] + 14.0*self["CH2plus_density"] + \
                29.0*self["HCOplus_density"] + 15.0*self["CH3plus_density"]

        return (metal_density / (self["density"] * my_chemistry.density_units))


    def calculate_mean_molecular_weight(self):
        my_chemistry = self.chemistry_data
        if (self["energy"] == 0).all():
            return np.ones(self["energy"].size)
        self.calculate_temperature()
        return (self["temperature"] / \
                (self["energy"] * (my_chemistry.Gamma - 1.) *
                 self.chemistry_data.temperature_units))

    def calculate_cooling_time(self):
        calculate_cooling_time(self)

    def calculate_gamma(self):
        calculate_gamma(self)

    def calculate_pressure(self):
        calculate_pressure(self)

    def calculate_temperature(self):
        calculate_temperature(self)

    def solve_chemistry(self, dt):
        solve_chemistry(self, dt)

class FieldNotFound(Exception):
    def __init__(self, field):
        self.field = field

    def __str__(self):
        return "Field '%s' not found!" % ((self.field, ))

class NotAGrid(Exception):
    def __str__(self):
        return "This routine needs a yt grid."


_grackle_to_yt = {
    'density': ('gas', 'density'),
    'n_H' : ('gas', 'n_H'),
    'HI': ('gas', 'H_p0_density'),
    'HII': ('gas', 'H_p1_density'),
    'HM': ('gas', 'HM_density'),
    'HeI': ('gas', 'He_p0_density'),
    'HeII': ('gas', 'He_p1_density'),
    'HeIII': ('gas', 'He_p2_density'),
    'H2I': ('gas', 'H2_p0_density'),
    'H2II': ('gas', 'H2_p1_density'),
    'DI': ('gas', 'D_p0_density'),
    'DII': ('gas', 'D_p1_density'),
    'HDI': ('gas', 'HD_p0_density'),
    'de': ('gas', 'El_density'),
    'metal': ('gas', 'metal_density'),
    'x-velocity': ('gas', 'velocity_x'),
    'y-velocity': ('gas', 'velocity_y'),
    'z-velocity': ('gas', 'velocity_z'),
    'energy': ('gas', 'thermal_energy'),
    'Water_density'    : ('gas', 'H2OI_density'),
    'O_density'        : ('gas', 'OI_density'),
    'OH_density'       : ('gas', 'OHI_density'),
    'O2_density'       : ('gas', 'O2I_density'),
    'Oplus_density'    : ('gas', 'OII_density'),
    'OHplus_density'   : ('gas', 'OHII_density'),
    'H2Oplus_density'  : ('gas', 'H2OII_density'),
    'H3Oplus_density'  : ('gas', 'H3OII_density'),
    'O2plus_density'   : ('gas', 'O2II_density'),
    'Cplus_density'    : ('gas', 'CII_density'),
    'C_density'        : ('gas', 'CI_density'),
    'CH_density'       : ('gas', 'CHI_density'),  
    'CH2_density'      : ('gas', 'CH2I_density'),
    'CH3_density'      : ('gas', 'CH3I_density'),
    'CH4_density'      : ('gas', 'CH4I_density'),
    'CO_density'       : ('gas', 'COI_density'),
    'COplus_density'   : ('gas', 'COII_density'),
    'CO2_density'      : ('gas', 'CO2I_density'),
    'CHplus_density'   : ('gas', 'CHII_density'),
    'CH2plus_density'  : ('gas', 'CH2II_density'),
    'H3plus_density'   : ('gas', 'H3II_density'),
    'HCOplus_density'  : ('gas', 'HCOplus_density'),
    'HeHplus_density'  : ('gas', 'HeHplus_density'),
    'CH3plus_density'  : ('gas', 'CH3plus_density'),
    'CH4plus_density'  : ('gas', 'CH4plus_density'),
    'CH5plus_density'  : ('gas', 'CH5plus_density'),
    'O2Hplus_density'  : ('gas', 'O2Hplus_density')

}

_skip = ("pressure", "temperature", "cooling_time", "gamma")

_yt_to_grackle = dict((b, a) for a, b in _grackle_to_yt.items())

def _units(chemistry_data, fname):
    if fname[1].endswith("density"):
        return chemistry_data.density_units
    elif fname[1].endswith("energy"):
        energy_units = (chemistry_data.velocity_units)**2.0
        return energy_units
    elif fname[1].startswith("velocity"):
        v_units = (chemistry_data.velocity_units)
        return v_units
    else:
        raise FieldNotFound(fname)

def _needed_fields(fc):
    cd = fc.chemistry_data
    for f in sorted(fc):
        if f in _skip: continue
        f2 = _grackle_to_yt[f]
        yield f, f2, _units(cd, f2)

def grid_to_grackle(chemistry_data, grid, update = True):
    if not hasattr(grid, "ActiveDimensions"):
        raise RuntimeError
    ds = grid.ds
    fields = ds.derived_field_list + ds.field_list
    ni = grid.ActiveDimensions[0]
    fc = FluidContainer(chemistry_data, ni)
    for f1, f2, conv in _needed_fields(fc):
        if f2 not in fields:
            raise FieldNotFound(f2)
    for j in range(grid.ActiveDimensions[1]):
        for k in range(grid.ActiveDimensions[2]):
            for f1, f2, conv in _needed_fields(fc):
                fc[f1][:] = grid[f2][:,j,k] / conv
            yield fc
            if not update: continue
            for f1, f2, conv in _needed_fields(fc):
                grid[f2][:,j,k] = fc[f1][:] * conv
