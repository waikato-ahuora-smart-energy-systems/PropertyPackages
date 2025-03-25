from HumidAirSurrogate import HAirParameterBlock
from typing import List, Literal

def build_humid_air_package(compound_list: List[str]):
    if len(compound_list) != 2:
        raise ValueError("Humid Air package only supports two components: water and air")
    for compound in compound_list:
        if compound not in ["water", "air"]:
            raise ValueError(f"Compound {compound} not supported in Humid Air package")

    return HAirParameterBlock()