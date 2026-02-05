from typing import List
from idaes.models.properties.general_helmholtz import (
    registered_components,
    PhaseType,
    StateVars,
    AmountBasis
)
from .helmholtz_extended import HelmholtzExtendedParameterBlock
from .parameters import register_compounds

# Add the "parameters" directory to the path so that the Helmholtz EOS can find the parameter files.
register_compounds()


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
            return HelmholtzExtendedParameterBlock(pure_component=component,
                                            phase_presentation=PhaseType.LG,
                                            state_vars=StateVars.PH,
                                            amount_basis=AmountBasis.MOLE)
        else:
            raise ValueError(f"Compound {component} not found in Helmholtz registered components")
