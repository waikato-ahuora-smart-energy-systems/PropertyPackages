from typing import List, Literal
from idaes.models.properties.general_helmholtz import (
    registered_components,
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
    AmountBasis
)

# TODO: add a new directory in here for storing the parameters of more helmholtz components. 
# then use set_parameter_path to point to the directory.
# https://idaes-pse.readthedocs.io/en/2.4.0/reference_guides/model_libraries/generic/property_models/helmholtz.html#idaes.models.properties.general_helmholtz.set_parameter_path


def build_helmholtz_package(compound_list: List[str]):
    # TODO: Support additional args, based on what additional parameters can be passed (e.g phaseType etc)
    if(len(compound_list) == 0):
        raise ValueError("No compounds provided")
    elif(len(compound_list) > 1):
        raise ValueError("Helmholtz EOS only supports single component systems")
    else:
        component = compound_list[0]
        if component in registered_components():
            # Build and return a helmholtz property package for this compound
            return HelmholtzParameterBlock(pure_component=component,
                                            phase_presentation=PhaseType.MIX,
                                            state_vars=StateVars.PH,
                                            amount_basis=AmountBasis.MOLE)