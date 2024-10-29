from idaes.models.properties.modular_properties.phase_equil.bubble_dew import LogBubbleDew
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.pure import RPP4, Perrys
from idaes.core import LiquidPhase, VaporPhase, Component
from pyomo.environ import units as pyunits
from .chem_sep import ChemSep

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
        'benzene': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (29525.0, pyunits.J/pyunits.kilomol/pyunits.K),
                    'B': (-51.417, pyunits.J/pyunits.kilomol/pyunits.K**2),
                    'C': (1.1944, pyunits.J/pyunits.kilomol/pyunits.K**3),
                    'D': (-0.0016468, pyunits.J/pyunits.kilomol/pyunits.K**4),
                    'E': (6.8461e-07, pyunits.J/pyunits.kilomol/pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.99938, pyunits.kmol/pyunits.m**3),
                    '2': (0.26348, None),
                    '3': (562.05, pyunits.K),
                    '4': (0.27856, None),
                    'eqn_type': 1
                },
                'enth_mol_form_vap_comp_ref': (82880000.0, pyunits.J/pyunits.kilomol),
                'entr_mol_form_vap_comp_ref': (-269300.0, pyunits.J/pyunits.kilomol/pyunits.K),
                'mw': (78.11184, pyunits.kg/pyunits.kilomol),
                'omega': 0.209,
                'pressure_crit': (4895000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (88.368, None),
                    'B': (-6712.9, pyunits.K),
                    'C': (-10.022, pyunits.K)
                },
                'temperature_crit': (562.05, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        },
        'toluene': {
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': RPP4,
            'entr_mol_ig_comp': RPP4,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (-43647.49, pyunits.J/pyunits.kilomol/pyunits.K),
                    'B': (603.542, pyunits.J/pyunits.kilomol/pyunits.K**2),
                    'C': (-0.399451, pyunits.J/pyunits.kilomol/pyunits.K**3),
                    'D': (0.000104382, pyunits.J/pyunits.kilomol/pyunits.K**4)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.89799, pyunits.kmol/pyunits.m**3),
                    '2': (0.27359, None),
                    '3': (591.75, pyunits.K),
                    '4': (0.30006, None),
                    'eqn_type': 1
                },
                'enth_mol_form_vap_comp_ref': (50170000.0, pyunits.J/pyunits.kilomol),
                'entr_mol_form_vap_comp_ref': (-320990.0, pyunits.J/pyunits.kilomol/pyunits.K),
                'mw': (92.13843, pyunits.kg/pyunits.kilomol),
                'omega': 0.264,
                'pressure_crit': (4108000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (32.89891, None),
                    'B': (-5013.81, pyunits.K),
                    'C': (-1.348918, pyunits.K)
                },
                'temperature_crit': (591.75, pyunits.K)
            },
            'phase_equilibrium_form': {('Vap', 'Liq'): log_fugacity},
            'pressure_sat_comp': ChemSep,
            'type': Component
        }
    },
    'parameter_data': {
        'PR_kappa': {
            ('benzene', 'toluene'): 0,
            ('toluene', 'benzene'): 0,
            ('benzene', 'benzene'): 0,
            ('toluene', 'toluene'): 0
        }
    },
    'phase_equilibrium_state': {('Vap', 'Liq'): SmoothVLE},
    'phases': {
        'Liq': {
            'equation_of_state': Cubic,
            'equation_of_state_options': {'type': CubicType.PR},
            'type': LiquidPhase
        },
        'Vap': {
            'equation_of_state': Cubic,
            'equation_of_state_options': {'type': CubicType.PR},
            'type': VaporPhase
        }
    },
    'phases_in_equilibrium': [('Vap', 'Liq')],
    'pressure_ref': (101325, pyunits.Pa),
    'state_bounds': {
        'flow_mol': (0, 100, 1000, pyunits.mol / pyunits.s),
        'pressure': (50000.0, 100000.0, 1000000.0, pyunits.Pa),
        'temperature': (273.15, 300, 500, pyunits.K)
    },
    'state_definition': FTPx,
    'temperature_ref': (298.15, pyunits.K)
}