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

config = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    'bubble_dew_method': LogBubbleDew,
    'components': {
        'benzene': {
            'cp_mol_ig_comp': ChemSep,
            'dens_mol_liq_comp': Perrys,
            'enth_mol_ig_comp': ChemSep,
            'enth_mol_liq_comp': ChemSep,
            'entr_mol_ig_comp': ChemSep,
            'entr_mol_liq_comp': ChemSep,
            'parameter_data': {
                'cp_mol_ig_comp_coeff': {
                    'A': (34010.24, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-588.0978, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (12.81777, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.000197306, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (5.142899e-08, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'cp_mol_liq_comp_coeff': {
                    'A': (111460.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-1854.3, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (22.399, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.028936, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (2.8991e-05, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'dens_mol_liq_comp_coeff': {
                    '1': (0.99938, pyunits.kmol / pyunits.m**3),
                    '2': (0.26348, None),
                    '3': (562.05, pyunits.K),
                    '4': (0.27856, None),
                    'eqn_type': 1
                },
                'enth_entr_mol_ig_coeff': {
                    'A': (29525.0, pyunits.J / pyunits.kmol / pyunits.K),
                    'B': (-51.417, pyunits.J / pyunits.kmol / pyunits.K**2),
                    'C': (1.1944, pyunits.J / pyunits.kmol / pyunits.K**3),
                    'D': (-0.0016468, pyunits.J / pyunits.kmol / pyunits.K**4),
                    'E': (6.8461e-07, pyunits.J / pyunits.kmol / pyunits.K**5)
                },
                'enth_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol),
                'enth_mol_form_vap_comp_ref': (82880000.0, pyunits.J / pyunits.kmol),
                'entr_mol_form_liq_comp_ref': (0, pyunits.J / pyunits.kmol / pyunits.K),
                'entr_mol_form_vap_comp_ref': (-269300.0, pyunits.J / pyunits.kmol / pyunits.K),
                'mw': (78.11184, pyunits.kg / pyunits.kmol),
                'omega': 0.209,
                'pressure_crit': (4895000.0, pyunits.Pa),
                'pressure_sat_comp_coeff': {
                    'A': (21.075, None),
                    'B': (2977.3, pyunits.K),
                    'C': (-41.505, pyunits.K)
                },
                'temperature_crit': (562.05, pyunits.K)
            },
            'phase_equilibrium_form': {
                ('Vap', 'Liq'): log_fugacity
            },
            'pressure_sat_comp': ChemSep,
            'type': Component
        }
    },
    'parameter_data': {
        'PR_kappa': {('benzene', 'benzene'): 0.0}
    },
    'phase_equilibrium_state': {
        ('Liq', 'Vap'): SmoothVLE
    },
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
    'phases_in_equilibrium': [('Liq', 'Vap')],
    'pressure_ref': (101325, pyunits.Pa),
    'state_bounds': {
        'flow_mol': (0, 100, 1000, pyunits.mol/pyunits.s),
        'pressure': (50000.0, 100000.0, 1000000.0, pyunits.Pa),
        'temperature': (278.68, 300, 500, pyunits.K)
    },
    'state_definition': FTPx,
    'temperature_ref': (298.15, pyunits.K),
}