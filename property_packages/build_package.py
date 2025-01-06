from property_packages.modular.template_builder import build_config
from .helmholtz.helmholtz_builder import build_helmholtz_package
from property_packages.types import PackageName, States
from typing import List, Optional

def build_package(package_name: PackageName, compound_list: List[str], valid_states: List[States]=["Liq", "Vap"]): # type: ignore
    match package_name:
        case "peng-robinson":
            return build_config("peng-robinson", compound_list, valid_states)
        case "coolprop":
            return build_config("coolprop", compound_list, valid_states)
        case "helmholtz":
            return build_helmholtz_package(compound_list)
        case "nrtl":
            # TODO: Implement build_package for nrtl.
            raise ValueError("NRTL package not implemented yet.")
        case _:
            raise ValueError(f"Invalid package name {package_name}. Expected a valid property package type, e.g helmholtz")