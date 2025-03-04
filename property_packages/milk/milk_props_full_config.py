#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Example property package for the VLE calucations for a milk_solid-water-o-Xylene
system. If using the activity coefficient models (NRTL or Wilson), the user is
expected to provide the parameters necessary for these models. Please note that
these parameters are declared as variables here to allow for use in a parameter
estimation problem if the VLE data is available.
"""

# Import Pyomo libraries
from pyomo.environ import Param, NonNegativeReals, Set, units as pyunits

# Import IDAES cores
from idaes.core import declare_process_block_class, Component
from idaes.core.util.misc import extract_data

from idaes.models.properties.activity_coeff_models.activity_coeff_prop_pack import (
    ActivityCoeffParameterData,
)
from idaes.logger import getIdaesLogger

# Some more information about this module
__author__ = "Jaffer Ghouse"
__version__ = "0.0.1"


# Set up logger
_log = getIdaesLogger(__name__)


@declare_process_block_class("MilkParameterBlock")
class MilkParameterData(ActivityCoeffParameterData):
    """Property package for mixtures of water and  milk solids"""

    def build(self):
        """
        Callable method for Block construction.
        """
        self.component_list_master = Set(initialize=["milk_solid", "water"])

        # Create component objects
        # NOTE: User needs to update this list; can be a subset or
        # equal to the master component list
        self.milk_solid = Component()
        self.water = Component()

        super(MilkParameterData, self).build()

        # List of phase equilibrium index
        self.phase_equilibrium_idx_master = Set(initialize=[1, 2])

        self.phase_equilibrium_idx = Set(initialize=[1, 2])

        self.phase_equilibrium_list_master = {
            1: ["milk_solid", ("Vap", "Liq")],
            2: ["water", ("Vap", "Liq")],
        }

        self.phase_equilibrium_list = {
            1: ["milk_solid", ("Vap", "Liq")],
            2: ["water", ("Vap", "Liq")],
        }

        # Thermodynamic reference state
        self.pressure_reference = Param(
            mutable=True,
            default=101325,
            doc="Reference pressure [Pa]",
            units=pyunits.Pa,
        )
        self.temperature_reference = Param(
            mutable=True,
            default=298.15,
            doc="Reference temperature [K]",
            units=pyunits.K,
        )

        # Critical properties
        pressure_critical_data = {
            "milk_solid": 1332.96e3,  #https://www.chemeo.com/cid/13-615-4/Oleic-Acid aprrox as Oleic acid Joback method
            "water": 22046e3, # NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI
        }

        self.pressure_critical = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=True,
            initialize=extract_data(pressure_critical_data),
            doc="Critical pressure [Pa]",
            units=pyunits.Pa,
        )

        temperature_critical_data = {
            "milk_solid": 937.21, #https://www.chemeo.com/cid/13-615-4/Oleic-Acid aprrox as Oleic acid Joback method
            "water": 647.14, # NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI
        }

        self.temperature_critical = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=True,
            initialize=extract_data(temperature_critical_data),
            doc="Critical temperature [K]",
            units=pyunits.K,
        )


        mw_comp_data = {
            "milk_solid": 232e-3, # F.Glasser et al Technical Note: Estimation of Milk Fatty Acid Yield from Milk Fat Data https://www.sciencedirect.com/science/article/pii/S0022030207717241#:~:text=The%20mean%20molecular%20weight%20of,%3D%209%20g%2Fmol).
            "water": 18.02e-3,
        }

        self.mw_comp = Param(
            self.component_list,
            mutable=True,
            initialize=extract_data(mw_comp_data),
            doc="molecular weight kg/mol",
            units=pyunits.kg / pyunits.mol,
        )

        # Constants for specific heat capacity, enthalpy valid to 70C
        #https://cameochemicals.noaa.gov/chris/OLA.pdf Oleic acid
        # Sources: For water NIST Webbook, https://webbook.nist.gov addapted for water https://webbook.nist.gov/cgi/cbook.cgi?ID=C14940637&Type=JANAFG&Plot=on#JANAFG
        #

        Cp_Liq_A_data = { 
            ("milk_solid"): 470599.7666604,
              ("water"): 92072.5211902,}
        Cp_Liq_B_data = {
            ("milk_solid"): -6423.9816705,
            ("water"):-186.9823652,
        }
        Cp_Liq_C_data = {
            ("milk_solid"): 32.8747223,
            ("water"): 0.8662807,
        }
        Cp_Liq_D_data = {("milk_solid"): -0.0747202, 
                         ("water"): -0.0019851, }

        Cp_Liq_E_data = {("milk_solid"): 0.0000637,
                          ("water"): 0.0000019
,}

        self.cp_mol_liq_comp_coeff_A = Param(
            self.component_list,
            units=pyunits.J / pyunits.kmol / pyunits.K,
            doc="Liquid phase Cp parameter A",
            initialize=extract_data(Cp_Liq_A_data),
        )

        self.cp_mol_liq_comp_coeff_B = Param(
            self.component_list,
            units=pyunits.J / pyunits.kmol / pyunits.K**2,
            doc="Liquid phase Cp parameter B",
            initialize=extract_data(Cp_Liq_B_data),
        )

        self.cp_mol_liq_comp_coeff_C = Param(
            self.component_list,
            units=pyunits.J / pyunits.kmol / pyunits.K**3,
            doc="Liquid phase Cp parameter C",
            initialize=extract_data(Cp_Liq_C_data),
        )

        self.cp_mol_liq_comp_coeff_D = Param(
            self.component_list,
            units=pyunits.J / pyunits.kmol / pyunits.K**4,
            doc="Liquid phase Cp parameter D",
            initialize=extract_data(Cp_Liq_D_data),
        )

        self.cp_mol_liq_comp_coeff_E = Param(
            self.component_list,
            units=pyunits.J / pyunits.kmol / pyunits.K**5,
            doc="Liquid phase Cp parameter E",
            initialize=extract_data(Cp_Liq_E_data),
        )

#water Cp vap from Cengel and Boles Thermodynamics an Engineering Approach
#use it stays as a liquid for milk solids
        Cp_Vap_A_data = {
            ("milk_solid"): 39.98454639,
            ("water"): 39.98454639,
        }
        Cp_Vap_B_data = {
            ("milk_solid"):-0.06627961,
            ("water"): -0.06627961,
        }
        Cp_Vap_C_data = {
            ("milk_solid"): 0.00022787,
            ("water"): 0.00022787,
        }
        Cp_Vap_D_data = {
            ("milk_solid"): -0.00000030,
            ("water"): -0.00000030,
        }
        Cp_Vap_E_data = {("milk_solid"):  0.0000001,
                          ("water"): 0, 
                          }

        self.cp_mol_vap_comp_coeff_A = Param(
            self.component_list,
            units=pyunits.J / pyunits.mol / pyunits.K,
            doc="Vapor phase Cp parameter A",
            initialize=extract_data(Cp_Vap_A_data),
        )

        self.cp_mol_vap_comp_coeff_B = Param(
            self.component_list,
            units=pyunits.J / pyunits.mol / pyunits.K**2,
            doc="Vapor phase Cp parameter B",
            initialize=extract_data(Cp_Vap_B_data),
        )

        self.cp_mol_vap_comp_coeff_C = Param(
            self.component_list,
            units=pyunits.J / pyunits.mol / pyunits.K**3,
            doc="Vapor phase Cp parameter C",
            initialize=extract_data(Cp_Vap_C_data),
        )

        self.cp_mol_vap_comp_coeff_D = Param(
            self.component_list,
            units=pyunits.J / pyunits.mol / pyunits.K**4,
            doc="Vapor phase Cp parameter D",
            initialize=extract_data(Cp_Vap_D_data),
        )

        self.cp_mol_vap_comp_coeff_E = Param(
            self.component_list,
            units=pyunits.J / pyunits.mol / pyunits.K**5,
            doc="Vapor phase Cp parameter E",
            initialize=extract_data(Cp_Vap_E_data),
        )

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid

        #This is the Wagner Equation for the calculation of the saturation pressure 
        pressure_sat_coeff_data = {
            ("milk_solid", "A"): -13.07261788, #Calcuated from the aceentic factor of Oleic acid file:///C:/Users/bjl25/Downloads/ie202379u_si_001.pdf
            ("milk_solid", "B"): 6.31986368,# Eq https://www.sciencedirect.com/science/article/pii/S0021961411000905
            ("milk_solid", "C"): -13.50398608,
            ("milk_solid", "D"): -5.45874388,
            ("water", "A"): -7.76451, #https://pubs.acs.org/doi/pdf/10.1021/i200021a023 Correlation and Prediction of the Vapor Pressures of Pure Liquids over Large Pressure Ranges
            ("water", "B"): 1.45838,
            ("water", "C"): -2.7758,
            ("water", "D"): -1.23303,
        }

        self.pressure_sat_coeff = Param(
            self.component_list,
            ["A", "B", "C", "D"],
            mutable=False,
            initialize=extract_data(pressure_sat_coeff_data),
            doc="parameters to compute P_sat",
        )

        # Standard heats of formation
        # Source: NIST Webbook, https://webbook.nist.gov
        # Retrieved 25th September 2019
        dh_form_data = {
            ("Vap", "milk_solid"): -690e3,
            ("Vap", "water"): -241.8e3 ,
            ("Liq", "milk_solid"): -764.8e3,
            ("Liq", "water"): -285.83e3,
        }

        self.dh_form = Param(
            self.phase_list,
            self.component_list,
            mutable=True,
            initialize=extract_data(dh_form_data),
            doc="Standard heats of formation [J/mol]",
            units=pyunits.J / pyunits.mol,
        )

        # Standard entropy of formation
        # Retrieved 2025 from NIST Webbook, https://webbook.nist.gov
        ds_form_data = {
            ("Vap", "milk_solid"): 1,
            ("Vap", "water"): 188.835,
            ("Liq", "milk_solid"): 1,
            ("Liq", "water"): 69.95,
        }

        self.ds_form = Param(
            self.phase_list,
            self.component_list,
            mutable=True,
            initialize=extract_data(ds_form_data),
            doc="Standard entropy of formation [J/mol.K]",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )