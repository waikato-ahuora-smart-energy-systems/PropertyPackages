from .milk_modular import milk_configuration
from ..modular.modular_extended import GenericExtendedParameterBlock

from typing import List, Literal

def build_milk_package(compound_list: List[str]):
    if len(compound_list) != 2:
        raise ValueError("Milk package only supports two components: water and milk solids")
    for compound in compound_list:
        if compound not in ["water", "milk_solid"]:
            raise ValueError(f"Compound {compound} not supported in milk package")
    return GenericExtendedParameterBlock(**milk_configuration)