from typing import List
from idaes.models.properties.general_helmholtz import (
    registered_components,
    PhaseType,
    StateVars,
    AmountBasis
)
from ahuora_property_packages.surrogate.surrogate_extended import SurrogateExtendedParameterBlock


def build_surrogate_package(compound_list: List[str]):
    return SurrogateExtendedParameterBlock()