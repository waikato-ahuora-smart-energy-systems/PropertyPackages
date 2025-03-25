from property_packages.modular.template_builder import build_config
from .helmholtz.helmholtz_builder import build_helmholtz_package
from property_packages.types import PackageName, States
from typing import List, Optional
from property_packages.milk import build_milk_package
from property_packages.humid_air import build_humid_air_package

def build_package(package_name: PackageName, compound_list: List[str], valid_states: List[States]=["Liq", "Vap"]): # type: ignore

    if valid_states is None:
        valid_states = ["Liq", "Vap"]

    match package_name:
        case "peng-robinson":
            return build_config("peng-robinson", compound_list, valid_states)
        case "helmholtz":
            return build_helmholtz_package(compound_list)
        case "nrtl":
            pass # TODO: Implement build_package for nrtl.
        case "milk":
            return build_milk_package(compound_list)
        case "humid_air":
            return build_humid_air_package(compound_list)
        case _:
            raise ValueError(f"Invalid package name {package_name}. Expected a valid property package type, e.g helmholtz")