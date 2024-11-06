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
                    'A': (36492.0, EXPR),
                    'B': (66.383, EXPR),
                    'C': (0.51076, EXPR),
                    'D': (-0.00068154, EXPR),
                    'E': (2.6315e-07, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (100270.0, EXPR),
                    'B': (86.345, EXPR),
                    'C': (7.7333, EXPR),
                    'D': (0.00096546, EXPR),
                    'E': (2.0281e-05, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.98, EXPR),
                    '2': (0.25169, None),
                    '3': (419.54, PyomoUnit),
                    '4': (0.26645, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (-500000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-307400.0, EXPR),
                'mw': (56.10632, EXPR),
                'omega': 0.194,
                'pressure_crit': (4020000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (21.09, None),
                    'B': (2368.5, PyomoUnit),
                    'C': (-19.25, PyomoUnit)
                },
                'temperature_crit': (419.5, PyomoUnit)
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
                    'A': (62525.0, EXPR),
                    'B': (106.62, EXPR),
                    'C': (1.0534, EXPR),
                    'D': (-0.0014615, EXPR),
                    'E': (5.8492e-07, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (58419.0, EXPR),
                    'B': (89.259, EXPR),
                    'C': (10.549, EXPR),
                    'D': (0.0039271, EXPR),
                    'E': (-8.6181e-07, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.26107, EXPR),
                    '2': (0.16952, None),
                    '3': (537.3, PyomoUnit),
                    '4': (0.1874, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (-62890000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-425200.0, EXPR),
                'mw': (98.18607, EXPR),
                'omega': 0.343,
                'pressure_crit': (2920000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (20.78, None),
                    'B': (2900.8, PyomoUnit),
                    'C': (-53.334, PyomoUnit)
                },
                'temperature_crit': (537.3, PyomoUnit)
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
                    'A': (52406.0, EXPR),
                    'B': (102.7, EXPR),
                    'C': (0.8529, EXPR),
                    'D': (-0.0011866, EXPR),
                    'E': (4.7459e-07, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (120740.0, EXPR),
                    'B': (197.35, EXPR),
                    'C': (7.4671, EXPR),
                    'D': (0.012038, EXPR),
                    'E': (-7.6352e-06, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.42588, EXPR),
                    '2': (0.20073, None),
                    '3': (504.0, PyomoUnit),
                    '4': (0.21659, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (-41670000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-386300.0, EXPR),
                'mw': (84.15948, EXPR),
                'omega': 0.281,
                'pressure_crit': (3143000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (21.306, None),
                    'B': (2992.5, PyomoUnit),
                    'C': (-30.644, PyomoUnit)
                },
                'temperature_crit': (504.0, PyomoUnit)
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
                    'A': (72710.0, EXPR),
                    'B': (110.08, EXPR),
                    'C': (1.2552, EXPR),
                    'D': (-0.0017373, EXPR),
                    'E': (6.9559e-07, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (-24253.0, EXPR),
                    'B': (87.834, EXPR),
                    'C': (11.35, EXPR),
                    'D': (0.0032002, EXPR),
                    'E': (-1.2467e-06, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.44684, EXPR),
                    '2': (0.23463, None),
                    '3': (567.0, PyomoUnit),
                    '4': (0.24846, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (-81940000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-463700.0, EXPR),
                'mw': (112.2126, EXPR),
                'omega': 0.393,
                'pressure_crit': (2680000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (20.716, None),
                    'B': (3032.4, PyomoUnit),
                    'C': (-64.529, PyomoUnit)
                },
                'temperature_crit': (567.0, PyomoUnit)
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
                    'A': (42229.0, EXPR),
                    'B': (99.1, EXPR),
                    'C': (0.65169, EXPR),
                    'D': (-0.00091143, EXPR),
                    'E': (3.6426e-07, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (92753.0, EXPR),
                    'B': (117.21, EXPR),
                    'C': (8.6537, EXPR),
                    'D': (0.007447, EXPR),
                    'E': (-2.6759e-06, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.44335, EXPR),
                    '2': (0.18566, None),
                    '3': (473.43, PyomoUnit),
                    '4': (0.23587, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (-21620000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-346200.0, EXPR),
                'mw': (70.1329, EXPR),
                'omega': 0.237,
                'pressure_crit': (3560000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (21.064, None),
                    'B': (2629.6, PyomoUnit),
                    'C': (-27.471, PyomoUnit)
                },
                'temperature_crit': (464.8, PyomoUnit)
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
                    'A': (34738.0, EXPR),
                    'B': (-36.808, EXPR),
                    'C': (0.4706, EXPR),
                    'D': (-0.000553, EXPR),
                    'E': (2.0678e-07, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (68726.0, EXPR),
                    'B': (-1953.6, EXPR),
                    'C': (31.772, EXPR),
                    'D': (-0.10571, EXPR),
                    'E': (0.00019673, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (1.3749, EXPR),
                    '2': (0.23949, None),
                    '3': (305.43, PyomoUnit),
                    '4': (0.22875, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (-83820000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-229120.0, EXPR),
                'mw': (30.06904, EXPR),
                'omega': 0.099,
                'pressure_crit': (4872000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (20.962, None),
                    'B': (1665.8, PyomoUnit),
                    'C': (-7.8809, PyomoUnit)
                },
                'temperature_crit': (305.32, PyomoUnit)
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
                    'A': (35095.0, EXPR),
                    'B': (-73.018, EXPR),
                    'C': (0.48182, EXPR),
                    'D': (-0.00055948, EXPR),
                    'E': (2.0878e-07, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (68016.0, EXPR),
                    'B': (-22414.0, EXPR),
                    'C': (286.75, EXPR),
                    'D': (-1.1802, EXPR),
                    'E': (0.0017304, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (2.3782, EXPR),
                    '2': (0.29542, None),
                    '3': (282.36, PyomoUnit),
                    '4': (0.32456, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (52510000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-219200.0, EXPR),
                'mw': (28.05316, EXPR),
                'omega': 0.087,
                'pressure_crit': (5041000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (20.89, None),
                    'B': (1502.3, PyomoUnit),
                    'C': (-8.9148, PyomoUnit)
                },
                'temperature_crit': (282.34, PyomoUnit)
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
                    'A': (23052.64, EXPR),
                    'B': (33.74914, EXPR),
                    'C': (-0.0639907, EXPR),
                    'D': (5.1023e-05, EXPR),
                    'E': (-1.37699e-08, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (14510.0, EXPR),
                    'B': (-1191.1, EXPR),
                    'C': (156.51, EXPR),
                    'D': (-6.1773, EXPR),
                    'E': (0.087907, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (2.9613, EXPR),
                    '2': (0.25981, None),
                    '3': (33.19, PyomoUnit),
                    '4': (0.19104, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (0.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-130571.0, EXPR),
                'mw': (2.01588, EXPR),
                'omega': -0.215993,
                'pressure_crit': (1313000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (18.072, None),
                    'B': (137.26, PyomoUnit),
                    'C': (0.53751, PyomoUnit)
                },
                'temperature_crit': (33.19, PyomoUnit)
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
                    'A': (27862.0, EXPR),
                    'B': (148.69, EXPR),
                    'C': (0.45538, EXPR),
                    'D': (-0.00067339, EXPR),
                    'E': (2.6964e-07, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (89466.0, EXPR),
                    'B': (-323.61, EXPR),
                    'C': (12.827, EXPR),
                    'D': (-0.010476, EXPR),
                    'E': (2.5037e-05, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.89934, EXPR),
                    '2': (0.25371, None),
                    '3': (407.85, PyomoUnit),
                    '4': (0.25125, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (-134990000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-295500.0, EXPR),
                'mw': (58.1222, EXPR),
                'omega': 0.186,
                'pressure_crit': (3640000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (20.913, None),
                    'B': (2269.9, PyomoUnit),
                    'C': (-19.458, PyomoUnit)
                },
                'temperature_crit': (407.85, PyomoUnit)
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
                    'A': (37981.0, EXPR),
                    'B': (-74.622, EXPR),
                    'C': (0.3019, EXPR),
                    'D': (-0.00028327, EXPR),
                    'E': (9.0711e-08, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (61157.0, EXPR),
                    'B': (5034.1, EXPR),
                    'C': (-48.913, EXPR),
                    'D': (-0.22998, EXPR),
                    'E': (0.0022243, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (1.894, EXPR),
                    '2': (0.23603, None),
                    '3': (191.05, PyomoUnit),
                    '4': (0.21974, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (-74520000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-186270.0, EXPR),
                'mw': (16.04246, EXPR),
                'omega': 0.011,
                'pressure_crit': (4599000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (20.701, None),
                    'B': (1035.0, PyomoUnit),
                    'C': (1.2704, PyomoUnit)
                },
                'temperature_crit': (190.56, PyomoUnit)
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
                    'A': (46120.0, EXPR),
                    'B': (46.029, EXPR),
                    'C': (0.6699, EXPR),
                    'D': (-0.00087892, EXPR),
                    'E': (3.4372e-07, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (115150.0, EXPR),
                    'B': (-3564.7, EXPR),
                    'C': (41.067, EXPR),
                    'D': (-0.098803, EXPR),
                    'E': (0.0001183, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (1.0023, EXPR),
                    '2': (0.26457, None),
                    '3': (425.17, PyomoUnit),
                    '4': (0.27138, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (-125790000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-309910.0, EXPR),
                'mw': (58.1222, EXPR),
                'omega': 0.199,
                'pressure_crit': (3796000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (20.958, None),
                    'B': (2350.4, PyomoUnit),
                    'C': (-23.412, PyomoUnit)
                },
                'temperature_crit': (425.12, PyomoUnit)
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
                    'A': (31986.0, EXPR),
                    'B': (42.662, EXPR),
                    'C': (0.49978, EXPR),
                    'D': (-0.00065626, EXPR),
                    'E': (2.56e-07, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (87486.0, EXPR),
                    'B': (-13371.0, EXPR),
                    'C': (156.92, EXPR),
                    'D': (-0.5459, EXPR),
                    'E': (0.00068504, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (1.3186, EXPR),
                    '2': (0.27005, None),
                    '3': (369.86, PyomoUnit),
                    '4': (0.27852, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (-104680000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-270200.0, EXPR),
                'mw': (44.09562, EXPR),
                'omega': 0.152,
                'pressure_crit': (4248000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (21.04, None),
                    'B': (2072.9, PyomoUnit),
                    'C': (-13.18, PyomoUnit)
                },
                'temperature_crit': (369.83, PyomoUnit)
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
                    'A': (5083.56, EXPR),
                    'B': (225.639, EXPR),
                    'C': (-0.0999265, EXPR),
                    'D': (1.33106e-05, EXPR)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (79790.0, EXPR),
                    'B': (300.8, EXPR),
                    'C': (5.1342, EXPR),
                    'D': (0.0095615, EXPR),
                    'E': (1.2777e-05, EXPR)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.98129, EXPR),
                    '2': (0.22226, None),
                    '3': (365.58, PyomoUnit),
                    '4': (0.24039, None),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, EXPR),
                'enth_mol_form_vap_comp_ref': (20230000.0, EXPR),
                'entr_mol_form_liq_comp_ref': (0, EXPR),
                'entr_mol_form_vap_comp_ref': (-267000.0, EXPR),
                'mw': (42.07974, EXPR),
                'omega': 0.137588,
                'pressure_crit': (4600000.0, PyomoUnit),
                'pressure_sat_comp_coeff': {
                    'A': (20.886, None),
                    'B': (1932.7, PyomoUnit),
                    'C': (-18.939, PyomoUnit)
                },
                'temperature_crit': (364.85, PyomoUnit)
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
            ("hydrogen", "nbutane"): 0.000,
            ("hydrogen", "ibutane"): 0.000,
            ("hydrogen", "ethylene"): 0.000,
            ("hydrogen", "propene"): 0.000,
            ("hydrogen", "butene"): 0.000,
            ("hydrogen", "pentene"): 0.000,
            ("hydrogen", "hexene"): 0.000,
            ("hydrogen", "heptene"): 0.000,
            ("hydrogen", "octene"): 0.000,
            ("methane", "hydrogen"): 0.000,
            ("methane", "methane"): 0.000,
            ("methane", "ethane"): 0.000,
            ("methane", "propane"): 0.000,
            ("methane", "nbutane"): 0.000,
            ("methane", "ibutane"): 0.000,
            ("methane", "ethylene"): 0.000,
            ("methane", "propene"): 0.000,
            ("methane", "butene"): 0.000,
            ("methane", "pentene"): 0.000,
            ("methane", "hexene"): 0.000,
            ("methane", "heptene"): 0.000,
            ("methane", "octene"): 0.000,
            ("ethane", "hydrogen"): 0.000,
            ("ethane", "methane"): 0.000,
            ("ethane", "ethane"): 0.000,
            ("ethane", "propane"): 0.000,
            ("ethane", "nbutane"): 0.000,
            ("ethane", "ibutane"): 0.000,
            ("ethane", "ethylene"): 0.000,
            ("ethane", "propene"): 0.000,
            ("ethane", "butene"): 0.000,
            ("ethane", "pentene"): 0.000,
            ("ethane", "hexene"): 0.000,
            ("ethane", "heptene"): 0.000,
            ("ethane", "octene"): 0.000,
            ("propane", "hydrogen"): 0.000,
            ("propane", "methane"): 0.000,
            ("propane", "ethane"): 0.000,
            ("propane", "propane"): 0.000,
            ("propane", "nbutane"): 0.000,
            ("propane", "ibutane"): 0.000,
            ("propane", "ethylene"): 0.000,
            ("propane", "propene"): 0.000,
            ("propane", "butene"): 0.000,
            ("propane", "pentene"): 0.000,
            ("propane", "hexene"): 0.000,
            ("propane", "heptene"): 0.000,
            ("propane", "octene"): 0.000,
            ("nbutane", "hydrogen"): 0.000,
            ("nbutane", "methane"): 0.000,
            ("nbutane", "ethane"): 0.000,
            ("nbutane", "propane"): 0.000,
            ("nbutane", "nbutane"): 0.000,
            ("nbutane", "ibutane"): 0.000,
            ("nbutane", "ethylene"): 0.000,
            ("nbutane", "propene"): 0.000,
            ("nbutane", "butene"): 0.000,
            ("nbutane", "pentene"): 0.000,
            ("nbutane", "hexene"): 0.000,
            ("nbutane", "heptene"): 0.000,
            ("nbutane", "octene"): 0.000,
            ("ibutane", "hydrogen"): 0.000,
            ("ibutane", "methane"): 0.000,
            ("ibutane", "ethane"): 0.000,
            ("ibutane", "propane"): 0.000,
            ("ibutane", "nbutane"): 0.000,
            ("ibutane", "ibutane"): 0.000,
            ("ibutane", "ethylene"): 0.000,
            ("ibutane", "propene"): 0.000,
            ("ibutane", "butene"): 0.000,
            ("ibutane", "pentene"): 0.000,
            ("ibutane", "hexene"): 0.000,
            ("ibutane", "heptene"): 0.000,
            ("ibutane", "octene"): 0.000,
            ("ethylene", "hydrogen"): 0.000,
            ("ethylene", "methane"): 0.000,
            ("ethylene", "ethane"): 0.000,
            ("ethylene", "propane"): 0.000,
            ("ethylene", "nbutane"): 0.000,
            ("ethylene", "ibutane"): 0.000,
            ("ethylene", "ethylene"): 0.000,
            ("ethylene", "propene"): 0.000,
            ("ethylene", "butene"): 0.000,
            ("ethylene", "pentene"): 0.000,
            ("ethylene", "hexene"): 0.000,
            ("ethylene", "heptene"): 0.000,
            ("ethylene", "octene"): 0.000,
            ("propene", "hydrogen"): 0.000,
            ("propene", "methane"): 0.000,
            ("propene", "ethane"): 0.000,
            ("propene", "propane"): 0.000,
            ("propene", "nbutane"): 0.000,
            ("propene", "ibutane"): 0.000,
            ("propene", "ethylene"): 0.000,
            ("propene", "propene"): 0.000,
            ("propene", "butene"): 0.000,
            ("propene", "pentene"): 0.000,
            ("propene", "hexene"): 0.000,
            ("propene", "heptene"): 0.000,
            ("propene", "octene"): 0.000,
            ("butene", "hydrogen"): 0.000,
            ("butene", "methane"): 0.000,
            ("butene", "ethane"): 0.000,
            ("butene", "propane"): 0.000,
            ("butene", "nbutane"): 0.000,
            ("butene", "ibutane"): 0.000,
            ("butene", "ethylene"): 0.000,
            ("butene", "propene"): 0.000,
            ("butene", "butene"): 0.000,
            ("butene", "pentene"): 0.000,
            ("butene", "hexene"): 0.000,
            ("butene", "heptene"): 0.000,
            ("butene", "octene"): 0.000,
            ("pentene", "hydrogen"): 0.000,
            ("pentene", "methane"): 0.000,
            ("pentene", "ethane"): 0.000,
            ("pentene", "propane"): 0.000,
            ("pentene", "nbutane"): 0.000,
            ("pentene", "ibutane"): 0.000,
            ("pentene", "ethylene"): 0.000,
            ("pentene", "propene"): 0.000,
            ("pentene", "butene"): 0.000,
            ("pentene", "pentene"): 0.000,
            ("pentene", "hexene"): 0.000,
            ("pentene", "heptene"): 0.000,
            ("pentene", "octene"): 0.000,
            ("hexene", "hydrogen"): 0.000,
            ("hexene", "methane"): 0.000,
            ("hexene", "ethane"): 0.000,
            ("hexene", "propane"): 0.000,
            ("hexene", "nbutane"): 0.000,
            ("hexene", "ibutane"): 0.000,
            ("hexene", "ethylene"): 0.000,
            ("hexene", "propene"): 0.000,
            ("hexene", "butene"): 0.000,
            ("hexene", "pentene"): 0.000,
            ("hexene", "hexene"): 0.000,
            ("hexene", "heptene"): 0.000,
            ("hexene", "octene"): 0.000,
            ("heptene", "hydrogen"): 0.000,
            ("heptene", "methane"): 0.000,
            ("heptene", "ethane"): 0.000,
            ("heptene", "propane"): 0.000,
            ("heptene", "nbutane"): 0.000,
            ("heptene", "ibutane"): 0.000,
            ("heptene", "ethylene"): 0.000,
            ("heptene", "propene"): 0.000,
            ("heptene", "butene"): 0.000,
            ("heptene", "pentene"): 0.000,
            ("heptene", "hexene"): 0.000,
            ("heptene", "heptene"): 0.000,
            ("heptene", "octene"): 0.000,
            ("octene", "hydrogen"): 0.000,
            ("octene", "methane"): 0.000,
            ("octene", "ethane"): 0.000,
            ("octene", "propane"): 0.000,
            ("octene", "nbutane"): 0.000,
            ("octene", "ibutane"): 0.000,
            ("octene", "ethylene"): 0.000,
            ("octene", "propene"): 0.000,
            ("octene", "butene"): 0.000,
            ("octene", "pentene"): 0.000,
            ("octene", "hexene"): 0.000,
            ("octene", "heptene"): 0.000,
            ("octene", "octene"): 0.000,
        }
    },
}
