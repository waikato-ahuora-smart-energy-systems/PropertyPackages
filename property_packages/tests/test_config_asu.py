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

config = {
    'base_units': {
        'amount': pyunits.mol,
        'length': pyunits.m,
        'mass': pyunits.kg,
        'temperature': pyunits.K,
        'time': pyunits.s
    },
    'bubble_dew_method': LogBubbleDew,
    'components': {
        'argon': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (20786.0, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (0.0, pyunits.J / (pyunits.mol * pyunits.K)),
                    'C': (0.0, pyunits.J / (pyunits.mol * pyunits.K)),
                    'D': (0.0, pyunits.J / (pyunits.mol * pyunits.K)),
                    'E': (0.0, pyunits.J / (pyunits.mol * pyunits.K))
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (46085.0, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (-1304.5, pyunits.J / (pyunits.mol * pyunits.K)),
                    'C': (21.195, pyunits.J / (pyunits.mol * pyunits.K)),
                    'D': (-0.015382, pyunits.J / (pyunits.mol * pyunits.K)),
                    'E': (3.3063e-05, pyunits.J / (pyunits.mol * pyunits.K))
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (3.803, pyunits.mol / pyunits.L),
                    '2': (0.286, pyunits.mol / pyunits.L),
                    '3': (150.86, pyunits.kg / pyunits.m**3),
                    '4': (0.2984, pyunits.mol / pyunits.L),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.mol),
                'enth_mol_form_vap_comp_ref': (0.0, pyunits.J / pyunits.mol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / (pyunits.mol * pyunits.K)),
                'entr_mol_form_vap_comp_ref': (-154732.0, pyunits.J / (pyunits.mol * pyunits.K)),
                'mw': (39.948, pyunits.g / pyunits.mol),
                'omega': -0.002,
                'pressure_crit': (4898000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (20.673, pyunits.bar),
                    'B': (804.37, pyunits.Pa),
                    'C': (0.74965, pyunits.Pa)
                },
                'temperature_crit': (150.86, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        'nitrogen': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (29425.0, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (-2.1701, pyunits.J / (pyunits.mol * pyunits.K)),
                    'C': (0.00058201, pyunits.J / (pyunits.mol * pyunits.K)),
                    'D': (1.3054e-05, pyunits.J / (pyunits.mol * pyunits.K)),
                    'E': (-8.2313e-09, pyunits.J / (pyunits.mol * pyunits.K))
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (55135.0, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (217.45, pyunits.J / (pyunits.mol * pyunits.K)),
                    'C': (-0.9071, pyunits.J / (pyunits.mol * pyunits.K)),
                    'D': (0.05327, pyunits.J / (pyunits.mol * pyunits.K)),
                    'E': (0.00024166, pyunits.J / (pyunits.mol * pyunits.K))
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (2.435, pyunits.mol / pyunits.L),
                    '2': (0.25137, pyunits.mol / pyunits.L),
                    '3': (126.27, pyunits.kg / pyunits.m**3),
                    '4': (0.249, pyunits.mol / pyunits.L),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.mol),
                'enth_mol_form_vap_comp_ref': (0.0, pyunits.J / pyunits.mol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / (pyunits.mol * pyunits.K)),
                'entr_mol_form_vap_comp_ref': (-131730.0, pyunits.J / (pyunits.mol * pyunits.K)),
                'mw': (28.014, pyunits.g / pyunits.mol),
                'omega': -0.001,
                'pressure_crit': (3390000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (21.219, pyunits.bar),
                    'B': (805.20, pyunits.Pa),
                    'C': (0.75547, pyunits.Pa)
                },
                'temperature_crit': (126.19, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        'oxygen': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (30182.0, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (-14.916, pyunits.J / (pyunits.mol * pyunits.K)),
                    'C': (0.054709, pyunits.J / (pyunits.mol * pyunits.K)),
                    'D': (-4.997e-05, pyunits.J / (pyunits.mol * pyunits.K)),
                    'E': (1.4883e-08, pyunits.J / (pyunits.mol * pyunits.K))
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (53393.0, pyunits.J / (pyunits.mol * pyunits.K)),
                    'B': (-1966.4, pyunits.J / (pyunits.mol * pyunits.K)),
                    'C': (48.21, pyunits.J / (pyunits.mol * pyunits.K)),
                    'D': (-0.31631, pyunits.J / (pyunits.mol * pyunits.K)),
                    'E': (0.0010466, pyunits.J / (pyunits.mol * pyunits.K))
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (2.6097, pyunits.mol / pyunits.L),
                    '2': (0.23614, pyunits.mol / pyunits.L),
                    '3': (154.78, pyunits.kg / pyunits.m**3),
                    '4': (0.23695, pyunits.mol / pyunits.L),
                    'eqn_type': 1
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.mol),
                'enth_mol_form_vap_comp_ref': (0.0, pyunits.J / pyunits.mol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / (pyunits.mol * pyunits.K)),
                'entr_mol_form_vap_comp_ref': (-205043.0, pyunits.J / (pyunits.mol * pyunits.K)),
                'mw': (32.0, pyunits.g / pyunits.mol),
                'omega': -0.0021,
                'pressure_crit': (5040000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (22.133, pyunits.bar),
                    'B': (788.18, pyunits.Pa),
                    'C': (0.76725, pyunits.Pa)
                },
                'temperature_crit': (154.58, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        }
    },
    'parameter_data': {
        'PR_kappa': {
            ('argon', 'argon'): 0.0,
            ('argon', 'nitrogen'): -0.01,
            ('argon', 'oxygen'): 0.0089,
            ('nitrogen', 'argon'): -0.01,
            ('nitrogen', 'nitrogen'): 0.0,
            ('nitrogen', 'oxygen'): -0.01,
            ('oxygen', 'argon'): 0.0089,
            ('oxygen', 'nitrogen'): -0.01,
            ('oxygen', 'oxygen'): 0.0
        }
    },
    'phase_equilibrium_form': {
        ('Vap', 'Liq'): log_fugacity
    },
    'pressure_sat_comp': ChemSep,
    'type': Component,
    'phase_equilibrium_state': {
        ('Vap', 'Liq'): SmoothVLE
    },
    'phases': {
        'Liq': {
            'equation_of_state': Cubic,
            'equation_of_state_options': {
                'type': CubicType.PR
            },
            'type': LiquidPhase
        },
        'Vap': {
            'equation_of_state': Cubic,
            'equation_of_state_options': {
                'type': CubicType.PR
            },
            'type': VaporPhase
        }
    },
    'phases_in_equilibrium': [('Vap', 'Liq')],
    'pressure_ref': (101325, pyunits.Pa),
    'state_bounds': {
        'flow_mol': (0, 100, 1000, pyunits.mol),
        'pressure': (50000.0, 100000.0, 1000000.0, pyunits.Pa),
        'temperature': (10, 300, 500, pyunits.K)
    },
    'state_definition': FTPx,
    'temperature_ref': (298.15, pyunits.K)
}