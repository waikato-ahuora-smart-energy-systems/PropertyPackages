


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

from property_packages.modular.builder.data.chem_sep import ChemSep

config = {
    'base_units': {
        'amount': pyunits.mol,
        'length': pyunits.m,
        'mass': pyunits.kg,
        'temperature': pyunits.K,
        'time': pyunits.s,
    },
    
    'bubble_dew_method': LogBubbleDew,

    'components': {
        'carbon dioxide': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': RPP4,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': RPP4,
            'entr_mol_liq_comp': ChemSep,
            'valid_phase_types': PT.vaporPhase,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (19800.0, pyunits.J / pyunits.mol / pyunits.K),
                    'B': (73.44, pyunits.J / pyunits.mol / pyunits.K**2),
                    'C': (-0.05602, pyunits.J / pyunits.mol / pyunits.K**3),
                    'D': (1.715e-05, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (80592.0, pyunits.J * pyunits.kmol**-1 * pyunits.K**-1),
                    'B': (108.83, pyunits.J * pyunits.kmol**-1 * pyunits.K**-2),
                    'C': (-6.9126, pyunits.J * pyunits.kmol**-1 * pyunits.K**-3),
                    'D': (0.059647, pyunits.J * pyunits.kmol**-1 * pyunits.K**-4),
                    'E': (6.9922e-06, pyunits.J * pyunits.kmol**-1 * pyunits.K**-5),
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (2.768, pyunits.kmol * pyunits.m**-3),
                    '2': (0.26212, pyunits.dimensionless),
                    '3': (304.21, pyunits.K),
                    '4': (0.2908, pyunits.dimensionless),
                    'eqn_type': 1,
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kilomol),
                'enth_mol_form_vap_comp_ref': (-393510000.0, pyunits.J / pyunits.kilomol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kilomol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-213677.0, pyunits.J / pyunits.kilomol / pyunits.K),
                'mw': (44.0095, pyunits.g / pyunits.mol),
                'omega': 0.223621,
                'pressure_crit': (7383000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (21.794, pyunits.dimensionless),
                    'B': (1725.4, pyunits.K),
                    'C': (-16.793, pyunits.K),
                },
                'temperature_crit': (304.21, pyunits.K),
            },
            'pressure_sat_comp': ChemSep,
            'type': Component,
        },
        
        'water': {
            'dens_mol_liq_comp': ChemSep,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (33444.62, pyunits.J / pyunits.kilomol / pyunits.K),
                    'B': (-5.799206, pyunits.J / pyunits.kilomol / pyunits.K ** 2),
                    'C': (0.0251681, pyunits.J / pyunits.kilomol / pyunits.K ** 3),
                    'D': (-1.43103e-05, pyunits.J / pyunits.kilomol / pyunits.K ** 4),
                    'E': (2.76249e-09, pyunits.J / pyunits.kilomol / pyunits.K ** 5),
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (75539.0, pyunits.J * pyunits.kmol**-1 * pyunits.K**-1),
                    'B': (-22297.0, pyunits.J * pyunits.kmol**-1 * pyunits.K**-2),
                    'C': (136.02, pyunits.J * pyunits.kmol**-1 * pyunits.K**-3),
                    'D': (-0.25622, pyunits.J * pyunits.kmol**-1 * pyunits.K**-4),
                    'E': (0.00018273, pyunits.J * pyunits.kmol**-1 * pyunits.K**-5),
                },
                'dens_mol_liq_comp_coeff': {
                    'A': (32.51621, pyunits.kmol * pyunits.m**-3),
                    'B': (-3.213004, pyunits.dimensionless),
                    'C': (7.92411, pyunits.dimensionless),
                    'D': (-7.359898, pyunits.dimensionless),
                    'E': (2.703522, pyunits.dimensionless),
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kilomol),
                'enth_mol_form_vap_comp_ref': (-241814000.0, pyunits.J / pyunits.kilomol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kilomol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-188724.0, pyunits.J / pyunits.kilomol / pyunits.K),
                'mw': (18.01528, pyunits.g / pyunits.mol),
                'omega': 0.344,
                'pressure_crit': (22064000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (23.401, pyunits.dimensionless),
                    'B': (3987.3, pyunits.K),
                    'C': (-37.161, pyunits.K),
                },
                'temperature_crit': (647.14, pyunits.K),
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component,
        }
    },

    'parameter_data': {
        'PR_kappa': {
            ('carbon dioxide', 'carbon dioxide'): 0.0,
            ('carbon dioxide', 'water'): -0.12155,
            ('water', 'carbon dioxide'): -0.12155,
            ('water', 'water'): 0.0,
        }
    },

    'phase_equilibrium_state': {('Vap', 'Liq'): SmoothVLE},

    'phases': {
        'Liq': {
            'equation_of_state': Cubic,
            'equation_of_state_options': {'type': CubicType.PR},
            'type': LiquidPhase,
        },
        'Vap': {
            'equation_of_state': Cubic,
            'equation_of_state_options': {'type': CubicType.PR},
            'type': VaporPhase,
        }
    },

    'phases_in_equilibrium': [('Vap', 'Liq')],
    
    'pressure_ref': (101325, pyunits.Pa),

    'state_bounds': {
        'flow_mol': (0, 100, 1000, pyunits.mol / pyunits.s),
        'pressure': (50000.0, 100000.0, 1000000.0, pyunits.Pa),
        'temperature': (216.58, 250, 400, pyunits.K),
    },

    'state_definition': FTPx,

    'temperature_ref': (298.15, pyunits.K),
}