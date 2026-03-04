from typing import List
from .seawater_extended import SeawaterExtendedParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock


valid_components = ["H2O", "TDS"]

def build_seawater_package(compound_list: List[str]):
    component = compound_list
    if component[0] not in valid_components or component[1] not in valid_components:
        raise ValueError(f"Compound {component} not found in valid components list for seawater package. Seawater package only supports the compounds 'H2O' and 'TDS'")
    return SeawaterExtendedParameterBlock()