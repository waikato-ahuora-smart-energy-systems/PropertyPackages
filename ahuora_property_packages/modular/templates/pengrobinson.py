from ..builder.common_parsers import *

template = {
    'base_units': base_units_parser,
    'bubble_dew_method': bubble_dew_method_parser,
    'components': components_parser,
    'parameter_data': pr_kappa_parser,
    'phase_equilibrium_state': phase_equilibrium_state_parser,
    'phases': phases_parser,
    'phases_in_equilibrium': phases_in_equilibrium_parser,
    'pressure_ref': pressure_ref_parser,
    'state_bounds': state_bounds_parser,
    'state_definition': state_definition_parser,
    'temperature_ref': temperature_ref_parser,
}