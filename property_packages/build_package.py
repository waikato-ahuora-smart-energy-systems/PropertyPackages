from property_packages.modular.template_builder import build_config
from .helmholtz.helmholtz_builder import build_helmholtz_package
from property_packages.humid_air import build_humid_air_package
from property_packages.milk import build_milk_package
from property_packages.types import PackageName, States
from typing import List

def build_package(package_name: PackageName, compound_list: List[str], valid_states: List[States]=["Liq", "Vap"]): # type: ignore
    """ Builds a property package

    Args:
        package_name (PackageName): Name of the property package to build.
        compound_list (List[str]): List of compound names to include in the package.
        valid_states (List[States], optional): List of valid states for the compounds.
    
    Returns:
        object: IDAES ParameterBlock object.

    Raises:
        ValueError: If the package name or states are invalid
    """
    
    # Type checking states
    for state in valid_states:
        if state not in ["Liq", "Vap"]:
            raise ValueError(f"Invalid state {state}. Valid states are: Liq, Vap")

    # Type checking package
    match package_name:
        case "peng-robinson":
            return build_config("peng-robinson", compound_list, valid_states)
        case "helmholtz":
            return build_helmholtz_package(compound_list)
        case "milk":
            return build_milk_package(compound_list)
        case "humid_air":
            return build_humid_air_package(compound_list)
        case "genericML":
            raise NotImplementedError("Generic ML package is not implemented yet.")
        case _:
            raise ValueError(f"Invalid package name {package_name}. Expected a valid property package type, e.g helmholtz")