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
Hydrocarbon processing phase equilibrium package using Peng-Robinson EoS.

Example property package using the Generic Property Package Framework.
This example shows how to set up a property package to do hydrocarbon
processing phase equilibrium in the generic framework using Peng-Robinson
equation along with methods drawn from the pre-built IDAES property libraries.

The example includes the dictionary named configuration contains parameters
for calculating VLE phase equilibrium and properties for hydrocarbon processing.
"""

# Import Python libraries
import logging

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component, PhaseType as PT
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    LogBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity

from idaes.models.properties.modular_properties.pure import Perrys
from idaes.models.properties.modular_properties.pure import RPP4
from idaes.models.properties.modular_properties.pure import RPP5

from ..modular.builder.data.chem_sep import ChemSep

# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for a Prng Robinson alkene system

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [2] Perry's Chemical Engineers' Handbook 7th Ed.
#     Converted to J/mol.K, mol/m^3
# [3] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 15th september, 2020
# [4] The Properties of Gases and Liquids (2001)
#     5th edition, Chemical Engineering Series - Robert C. Reid

configuration = {
    # Specifying components
    
    'components': {
        '1-butene': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (36492.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (66.383, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (0.51076, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.00068154, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (2.6315e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (100270.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (86.345, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (7.7333, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (0.00096546, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (2.0281e-05, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.98, pyunits.kmol / pyunits.m**3),
                    '2': (0.25169, None),
                    '3': (419.54, pyunits.K), 
                    '4': (0.26645, None),
                    'eqn_type': 1 
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (-500000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-307400.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (56.10632, pyunits.kg / pyunits.kmol), 
                'omega': 0.194, 
                'pressure_crit': (4020000.0, pyunits.Pa), 
                'pressure_sat_comp_coeff': {
                    'A': (21.09, None), 
                    'B': (2368.5, pyunits.K),
                    'C': (-19.25, pyunits.K),
                },
                'temperature_crit': (419.5, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        '1-heptene': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (62525.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (106.62, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (1.0534, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.0014615, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (5.8492e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (58419.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (89.259, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (10.549, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (0.0039271, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (-8.6181e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.26107, pyunits.kmol / pyunits.m**3),
                    '2': (0.16952, None),
                    '3': (537.3, pyunits.K),
                    '4': (0.1874, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (-62890000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-425200.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (98.18607, pyunits.kg / pyunits.kmol),
                'omega': 0.343,
                'pressure_crit': (2920000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (20.78, None),
                    'B': (2900.8, pyunits.K),
                    'C': (-53.334, pyunits.K)
                },
                'temperature_crit': (537.3, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        '1-hexene': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (52406.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (102.7, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (0.8529, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.0011866, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (4.7459e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (120740.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (197.35, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (7.4671, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (0.012038, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (-7.6352e-06, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.42588, pyunits.kmol / pyunits.m**3),
                    '2': (0.20073, None),
                    '3': (504.0, pyunits.K),
                    '4': (0.21659, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (-41670000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-386300.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (84.15948, pyunits.kg / pyunits.kmol),
                'omega': 0.281,
                'pressure_crit': (3143000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (21.306, None),
                    'B': (2992.5, pyunits.K),
                    'C': (-30.644, pyunits.K)
                },
                'temperature_crit': (504.0, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        '1-octene': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (72710.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (110.08, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (1.2552, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.0017373, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (6.9559e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (-24253.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (87.834, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (11.35, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (0.0032002, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (-1.2467e-06, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.44684, pyunits.kmol / pyunits.m**3),
                    '2': (0.23463, None),
                    '3': (567.0, pyunits.K),
                    '4': (0.24846, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (-81940000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-463700.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (112.2126, pyunits.kg / pyunits.kmol),
                'omega': 0.393,
                'pressure_crit': (2680000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (20.716, None),
                    'B': (3032.4, pyunits.K),
                    'C': (-64.529, pyunits.K)
                },
                'temperature_crit': (567.0, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        '1-pentene': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (42229.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (99.1, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (0.65169, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.00091143, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (3.6426e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (92753.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (117.21, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (8.6537, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (0.007447, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (-2.6759e-06, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.44335, pyunits.kmol / pyunits.m**3),
                    '2': (0.18566, None),
                    '3': (473.43, pyunits.K),
                    '4': (0.23587, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (-21620000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-346200.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (70.1329, pyunits.kg / pyunits.kmol),
                'omega': 0.237,
                'pressure_crit': (3560000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (21.064, None),
                    'B': (2629.6, pyunits.K),
                    'C': (-27.471, pyunits.K)
                },
                'temperature_crit': (464.8, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        'ethane': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (34738.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-36.808, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (0.4706, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.000553, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (2.0678e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (68726.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-1953.6, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (31.772, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.10571, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (0.00019673, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (1.3749, pyunits.kmol / pyunits.m**3),
                    '2': (0.23949, None),
                    '3': (305.43, pyunits.K),
                    '4': (0.22875, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (-83820000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-229120.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (30.06904, pyunits.kg / pyunits.kmol),
                'omega': 0.099,
                'pressure_crit': (4872000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (20.962, None),
                    'B': (1665.8, pyunits.K),
                    'C': (-7.8809, pyunits.K)
                },
                'temperature_crit': (305.32, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        'ethylene': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (35095.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-73.018, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (0.48182, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.00055948, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (2.0878e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (68016.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-22414.0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (286.75, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-1.1802, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (0.0017304, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (2.3782, pyunits.kmol / pyunits.m**3),
                    '2': (0.29542, None),
                    '3': (282.36, pyunits.K),
                    '4': (0.32456, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (52510000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-219200.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (28.05316, pyunits.kg / pyunits.kmol),
                'omega': 0.087,
                'pressure_crit': (5041000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (20.89, None),
                    'B': (1502.3, pyunits.K),
                    'C': (-8.9148, pyunits.K)
                },
                'temperature_crit': (282.34, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        'hydrogen': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (23052.64, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (33.74914, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (-0.0639907, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (5.1023e-05, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (-1.37699e-08, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (14510.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-1191.1, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (156.51, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-6.1773, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (0.087907, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (2.9613, pyunits.kmol / pyunits.m**3),
                    '2': (0.25981, None),
                    '3': (33.19, pyunits.K),
                    '4': (0.19104, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (0.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-130571.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (2.01588, pyunits.kg / pyunits.kmol),
                'omega': -0.215993,
                'pressure_crit': (1313000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (18.072, None),
                    'B': (137.26, pyunits.K),
                    'C': (0.53751, pyunits.K)
                },
                'temperature_crit': (33.19, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        'isobutane': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (27862.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (148.69, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (0.45538, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.00067339, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (2.6964e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (89466.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-323.61, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (12.827, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.010476, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (2.5037e-05, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.89934, pyunits.kmol / pyunits.m**3),
                    '2': (0.25371, None),
                    '3': (407.85, pyunits.K),
                    '4': (0.25125, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (-134990000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-295500.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (58.1222, pyunits.kg / pyunits.kmol),
                'omega': 0.186,
                'pressure_crit': (3640000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (20.913, None),
                    'B': (2269.9, pyunits.K),
                    'C': (-19.458, pyunits.K)
                },
                'temperature_crit': (407.85, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        'methane': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (37981.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-74.622, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (0.3019, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.00028327, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (9.0711e-08, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (61157.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (5034.1, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (-48.913, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.22998, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (0.0022243, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (1.894, pyunits.kmol / pyunits.m**3),
                    '2': (0.23603, None),
                    '3': (191.05, pyunits.K),
                    '4': (0.21974, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (-74520000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-186270.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (16.04246, pyunits.kg / pyunits.kmol),
                'omega': 0.011,
                'pressure_crit': (4599000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (20.701, None),
                    'B': (1035.0, pyunits.K),
                    'C': (1.2704, pyunits.K)
                },
                'temperature_crit': (190.56, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        'n-butane': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (46120.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (46.029, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (0.6699, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.00087892, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (3.4372e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (115150.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-3564.7, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (41.067, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.098803, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (0.0001183, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (1.0023, pyunits.kmol / pyunits.m**3),
                    '2': (0.26457, None),
                    '3': (425.17, pyunits.K),
                    '4': (0.27138, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (-125790000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-309910.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (58.1222, pyunits.kg / pyunits.kmol),
                'omega': 0.199,
                'pressure_crit': (3796000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (20.958, None),
                    'B': (2350.4, pyunits.K),
                    'C': (-23.412, pyunits.K)
                },
                'temperature_crit': (425.12, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        'propane': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (31986.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (42.662, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (0.49978, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.00065626, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (2.56e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (87486.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-13371.0, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (156.92, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.5459, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (0.00068504, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (1.3186, pyunits.kmol / pyunits.m**3),
                    '2': (0.27005, None),
                    '3': (369.86, pyunits.K),
                    '4': (0.27852, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (-104680000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-270200.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (44.09562, pyunits.kg / pyunits.kmol),
                'omega': 0.152,
                'pressure_crit': (4248000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (21.04, None),
                    'B': (2072.9, pyunits.K),
                    'C': (-13.18, pyunits.K)
                },
                'temperature_crit': (369.83, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        'propylene': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': RPP4,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': RPP4,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (5083.56, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (225.639, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (-0.0999265, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (1.33106e-05, pyunits.J / pyunits.kmol / pyunits.K**4)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (79790.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (300.8, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (5.1342, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (0.0095615, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (1.2777e-05, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.98129, pyunits.kmol / pyunits.m**3),
                    '2': (0.22226, None),
                    '3': (365.58, pyunits.K),
                    '4': (0.24039, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (20230000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-267000.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (42.07974, pyunits.kg / pyunits.kmol),
                'omega': 0.137588,
                'pressure_crit': (4600000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (20.886, None),
                    'B': (1932.7, pyunits.K),
                    'C': (-18.939, pyunits.K)
                },
                'temperature_crit': (364.85, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        }
    },
                     
    # Specifying phases
    "phases": {
        "Liq": {
            "type": LiquidPhase,
            "equation_of_state": Cubic,
            "equation_of_state_options": {"type": CubicType.PR},
        },
        "Vap": {
            "type": VaporPhase,
            "equation_of_state": Cubic,
            "equation_of_state_options": {"type": CubicType.PR},
        },
    },

    # Set base units of measurement
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
        "temperature": (273.15, 300, 1500, pyunits.K),
        "pressure": (5e4, 1e5, 1e7, pyunits.Pa),
    },

    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": LogBubbleDew,
    "parameter_data": {
        "PR_kappa": {
            ("hydrogen", "hydrogen"): 0.000,
            ("hydrogen", "methane"): 0.000,
            ("hydrogen", "ethane"): 0.000,
            ("hydrogen", "propane"): 0.000,
            ("hydrogen", "n-butane"): 0.000,
            ("hydrogen", "isobutane"): 0.000,
            ("hydrogen", "ethylene"): 0.000,
            ("hydrogen", "propylene"): 0.000,
            ("hydrogen", "1-butene"): 0.000,
            ("hydrogen", "1-pentene"): 0.000,
            ("hydrogen", "1-hexene"): 0.000,
            ("hydrogen", "1-heptene"): 0.000,
            ("hydrogen", "1-octene"): 0.000,
            ("methane", "hydrogen"): 0.000,
            ("methane", "methane"): 0.000,
            ("methane", "ethane"): 0.000,
            ("methane", "propane"): 0.000,
            ("methane", "n-butane"): 0.000,
            ("methane", "isobutane"): 0.000,
            ("methane", "ethylene"): 0.000,
            ("methane", "propylene"): 0.000,
            ("methane", "1-butene"): 0.000,
            ("methane", "1-pentene"): 0.000,
            ("methane", "1-hexene"): 0.000,
            ("methane", "1-heptene"): 0.000,
            ("methane", "1-octene"): 0.000,
            ("ethane", "hydrogen"): 0.000,
            ("ethane", "methane"): 0.000,
            ("ethane", "ethane"): 0.000,
            ("ethane", "propane"): 0.000,
            ("ethane", "n-butane"): 0.000,
            ("ethane", "isobutane"): 0.000,
            ("ethane", "ethylene"): 0.000,
            ("ethane", "propylene"): 0.000,
            ("ethane", "1-butene"): 0.000,
            ("ethane", "1-pentene"): 0.000,
            ("ethane", "1-hexene"): 0.000,
            ("ethane", "1-heptene"): 0.000,
            ("ethane", "1-octene"): 0.000,
            ("propane", "hydrogen"): 0.000,
            ("propane", "methane"): 0.000,
            ("propane", "ethane"): 0.000,
            ("propane", "propane"): 0.000,
            ("propane", "n-butane"): 0.000,
            ("propane", "isobutane"): 0.000,
            ("propane", "ethylene"): 0.000,
            ("propane", "propylene"): 0.000,
            ("propane", "1-butene"): 0.000,
            ("propane", "1-pentene"): 0.000,
            ("propane", "1-hexene"): 0.000,
            ("propane", "1-heptene"): 0.000,
            ("propane", "1-octene"): 0.000,
            ("n-butane", "hydrogen"): 0.000,
            ("n-butane", "methane"): 0.000,
            ("n-butane", "ethane"): 0.000,
            ("n-butane", "propane"): 0.000,
            ("n-butane", "n-butane"): 0.000,
            ("n-butane", "isobutane"): 0.000,
            ("n-butane", "ethylene"): 0.000,
            ("n-butane", "propylene"): 0.000,
            ("n-butane", "1-butene"): 0.000,
            ("n-butane", "1-pentene"): 0.000,
            ("n-butane", "1-hexene"): 0.000,
            ("n-butane", "1-heptene"): 0.000,
            ("n-butane", "1-octene"): 0.000,
            ("isobutane", "hydrogen"): 0.000,
            ("isobutane", "methane"): 0.000,
            ("isobutane", "ethane"): 0.000,
            ("isobutane", "propane"): 0.000,
            ("isobutane", "n-butane"): 0.000,
            ("isobutane", "isobutane"): 0.000,
            ("isobutane", "ethylene"): 0.000,
            ("isobutane", "propylene"): 0.000,
            ("isobutane", "1-butene"): 0.000,
            ("isobutane", "1-pentene"): 0.000,
            ("isobutane", "1-hexene"): 0.000,
            ("isobutane", "1-heptene"): 0.000,
            ("isobutane", "1-octene"): 0.000,
            ("ethylene", "hydrogen"): 0.000,
            ("ethylene", "methane"): 0.000,
            ("ethylene", "ethane"): 0.000,
            ("ethylene", "propane"): 0.000,
            ("ethylene", "n-butane"): 0.000,
            ("ethylene", "isobutane"): 0.000,
            ("ethylene", "ethylene"): 0.000,
            ("ethylene", "propylene"): 0.000,
            ("ethylene", "1-butene"): 0.000,
            ("ethylene", "1-pentene"): 0.000,
            ("ethylene", "1-hexene"): 0.000,
            ("ethylene", "1-heptene"): 0.000,
            ("ethylene", "1-octene"): 0.000,
            ("propylene", "hydrogen"): 0.000,
            ("propylene", "methane"): 0.000,
            ("propylene", "ethane"): 0.000,
            ("propylene", "propane"): 0.000,
            ("propylene", "n-butane"): 0.000,
            ("propylene", "isobutane"): 0.000,
            ("propylene", "ethylene"): 0.000,
            ("propylene", "propylene"): 0.000,
            ("propylene", "1-butene"): 0.000,
            ("propylene", "1-pentene"): 0.000,
            ("propylene", "1-hexene"): 0.000,
            ("propylene", "1-heptene"): 0.000,
            ("propylene", "1-octene"): 0.000,
            ("1-butene", "hydrogen"): 0.000,
            ("1-butene", "methane"): 0.000,
            ("1-butene", "ethane"): 0.000,
            ("1-butene", "propane"): 0.000,
            ("1-butene", "n-butane"): 0.000,
            ("1-butene", "isobutane"): 0.000,
            ("1-butene", "ethylene"): 0.000,
            ("1-butene", "propylene"): 0.000,
            ("1-butene", "1-butene"): 0.000,
            ("1-butene", "1-pentene"): 0.000,
            ("1-butene", "1-hexene"): 0.000,
            ("1-butene", "1-heptene"): 0.000,
            ("1-butene", "1-octene"): 0.000,
            ("1-pentene", "hydrogen"): 0.000,
            ("1-pentene", "methane"): 0.000,
            ("1-pentene", "ethane"): 0.000,
            ("1-pentene", "propane"): 0.000,
            ("1-pentene", "n-butane"): 0.000,
            ("1-pentene", "isobutane"): 0.000,
            ("1-pentene", "ethylene"): 0.000,
            ("1-pentene", "propylene"): 0.000,
            ("1-pentene", "1-butene"): 0.000,
            ("1-pentene", "1-pentene"): 0.000,
            ("1-pentene", "1-hexene"): 0.000,
            ("1-pentene", "1-heptene"): 0.000,
            ("1-pentene", "1-octene"): 0.000,
            ("1-hexene", "hydrogen"): 0.000,
            ("1-hexene", "methane"): 0.000,
            ("1-hexene", "ethane"): 0.000,
            ("1-hexene", "propane"): 0.000,
            ("1-hexene", "n-butane"): 0.000,
            ("1-hexene", "isobutane"): 0.000,
            ("1-hexene", "ethylene"): 0.000,
            ("1-hexene", "propylene"): 0.000,
            ("1-hexene", "1-butene"): 0.000,
            ("1-hexene", "1-pentene"): 0.000,
            ("1-hexene", "1-hexene"): 0.000,
            ("1-hexene", "1-heptene"): 0.000,
            ("1-hexene", "1-octene"): 0.000,
            ("1-heptene", "hydrogen"): 0.000,
            ("1-heptene", "methane"): 0.000,
            ("1-heptene", "ethane"): 0.000,
            ("1-heptene", "propane"): 0.000,
            ("1-heptene", "n-butane"): 0.000,
            ("1-heptene", "isobutane"): 0.000,
            ("1-heptene", "ethylene"): 0.000,
            ("1-heptene", "propylene"): 0.000,
            ("1-heptene", "1-butene"): 0.000,
            ("1-heptene", "1-pentene"): 0.000,
            ("1-heptene", "1-hexene"): 0.000,
            ("1-heptene", "1-heptene"): 0.000,
            ("1-heptene", "1-octene"): 0.000,
            ("1-octene", "hydrogen"): 0.000,
            ("1-octene", "methane"): 0.000,
            ("1-octene", "ethane"): 0.000,
            ("1-octene", "propane"): 0.000,
            ("1-octene", "n-butane"): 0.000,
            ("1-octene", "isobutane"): 0.000,
            ("1-octene", "ethylene"): 0.000,
            ("1-octene", "propylene"): 0.000,
            ("1-octene", "1-butene"): 0.000,
            ("1-octene", "1-pentene"): 0.000,
            ("1-octene", "1-hexene"): 0.000,
            ("1-octene", "1-heptene"): 0.000,
            ("1-octene", "1-octene"): 0.000,
        }
    },
}
