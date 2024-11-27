from pyomo.environ import units as pyunits
import xml.etree.ElementTree as ET
from pydantic import BaseModel
from typing import Dict, Any
import os


def convert_string_to_float(string: str) -> float | str:
    try:
        return float(string)
    except:
        return string

def convert_string_to_int(string: str) -> int | str:
    try:
        return int(string)
    except:
        return string


unit_dict = {
    "K": pyunits.K,
    "Pa": pyunits.Pa,
    "m3/kmol": pyunits.m**3 / pyunits.kmol,
    "kmol/m3": pyunits.kmol / pyunits.m**3,
    "_": pyunits.dimensionless,
    "kg/kmol": pyunits.kg / pyunits.kmol,
    "J/kmol": pyunits.J / pyunits.kmol,
    "J/kmol/K": pyunits.J / (pyunits.kmol * pyunits.K),
    "m": pyunits.m,
    "J0.5/m1.5": pyunits.J**0.5 / pyunits.m**1.5,
    "Coulomb.m": pyunits.C * pyunits.m, # TODO: update
    "m2/kmol": pyunits.m**2 / pyunits.kmol,
    "W/m/K": pyunits.W / (pyunits.m * pyunits.K),
    "N/m": pyunits.N / pyunits.m,
    "kg0.25.m3/s0.5/kmol": pyunits.kg**0.25 * pyunits.m**3 / (pyunits.s**0.5 * pyunits.kmol),
    "m3/kmol": pyunits.m**3 / pyunits.kmol
}


def get_unit_from_string(unit_str):
    if unit_str in unit_dict:
        return unit_dict[unit_str]
    else:
        return None


class UnitValuePair(BaseModel):
    name: str | None
    value: float | str | None
    unit: Any | None


# Helper functions
def parse_element(elem: ET.Element):
    try:
        return UnitValuePair(
            name=elem.get('name'), 
            value=convert_string_to_float(elem.get('value')), 
            unit=get_unit_from_string(elem.get('units')))
    except:
        return None

Coefficients = Dict[str, float | str | Any | None]

#
# Expected output
#
# Coefficients = {
#     "units": Any,
#     "eqno": int,
#     "A": float | None,
#     "B": float | None,
#     "C": float | None,
#     "D": float | None,
#     "E": float | None,
#     "Tmin": float,
#     "Tmax": float,
# }


def parse_coeff(element: ET.Element) -> Coefficients:
    self = {}
    for key in ['eqno', 'A', 'B', 'C', 'D', 'E', 'Tmin', 'Tmax']:
        child = element.find(key)
        if child is not None:
            value = child.get('value')
            if value is not None and key != 'eqno':
                self[key] = convert_string_to_float(value)
            elif value is not None and key == 'eqno':
                self[key] = convert_string_to_int(value)
        else:
            self[key] = 0.0
    self['units'] = get_unit_from_string(element.get('units'))
    return self


compound_template = {
    'LibraryIndex': parse_element,
    'CompoundID': parse_element,
    'StructureFormula': parse_element,
    'Family': parse_element,
    'CriticalTemperature': parse_element,
    'CriticalPressure': parse_element,
    'CriticalVolume': parse_element,
    'CriticalCompressibility': parse_element,
    'NormalBoilingPointTemperature': parse_element,
    'NormalMeltingPointTemperature': parse_element,
    'TriplePointTemperature': parse_element,
    'TriplePointPressure': parse_element,
    'MolecularWeight': parse_element,
    'LiquidVolumeAtNormalBoilingPoint': parse_element,
    'AcentricityFactor': parse_element,
    'RadiusOfGyration': parse_element,
    'SolubilityParameter': parse_element,
    'DipoleMoment': parse_element,
    'VanDerWaalsVolume': parse_element,
    'VanDerWaalsArea': parse_element,
    'HeatOfFormation': parse_element,
    'GibbsEnergyOfFormation': parse_element,
    'AbsEntropy': parse_element,
    'HeatOfFusionAtMeltingPoint': parse_element,
    'MatthiasCopemanC1': parse_element,
    'HeatOfCombustion': parse_element,
    'SolidDensity': parse_coeff,
    'LiquidDensity': parse_coeff,
    'VaporPressure': parse_coeff,
    'HeatOfVaporization': parse_coeff,
    'SolidHeatCapacityCp': parse_coeff,
    'LiquidHeatCapacityCp': parse_coeff,
    'IdealGasHeatCapacityCp': parse_coeff,
    'SecondVirialCoefficient': parse_coeff,
    'LiquidViscosity': parse_coeff,
    'VaporViscosity': parse_coeff,
    'LiquidThermalConductivity': parse_coeff,
    'VaporThermalConductivity': parse_coeff,
    'SurfaceTension': parse_coeff,
    'RPPHeatCapacityCp': parse_coeff,
    'RelativeStaticPermittivity': parse_coeff,
    'AntoineVaporPressure': parse_coeff,
    'LiquidViscosityRPS': parse_coeff,
}

Compound = Dict

def load_compound(name: str) -> Compound:
    """
    Parses the XML string and populates the object's attributes.

    Args:
        xml_string (str): The XML data as a string.
    """
    self = {}
    with open(os.path.dirname(__file__) + "/data_files/" + name.lower() + ".xml", 'r') as file:
        file = ''.join(file.readlines())

    root = ET.fromstring(file)

    # Iterating over all attributes in the template
    # and parsing the corresponding XML elements
    for attr in compound_template:
        if root.find(attr) is not None:
            self[attr] = compound_template[attr](root.find(attr))
        else:
            self[attr] = None
    return self
