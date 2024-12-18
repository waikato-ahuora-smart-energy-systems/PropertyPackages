import xml.etree.ElementTree as ET
from pydantic import BaseModel
from typing import Dict
from pyomo.environ import units
import os

# Hardcoded mapping of unit strings to unit objects
str_map = {
    "kg/mol": units.kg/units.mol,
    "m": units.m,
    "kg": units.kg,
    "mol": units.mol,
    "K": units.K,
    "m/s": units.m/units.s,
    "m^3/mol": units.m**3/units.mol,
    "Pa": units.Pa,
    "J/mol": units.J/units.mol,
    "kmol/m3": units.kmol/units.m**3,
    "J/kmol": units.J/units.kmol,
    "W/m/K": units.W/units.m/units.K,
    "Pa.s": units.Pa*units.s,
    "m3/kmol": units.m**3/units.kmol,
    "_": units.dimensionless,
    "None": None,
    "J/kmol/K": units.J/units.kmol/units.K,
    "m3/kmol": units.m**3/units.kmol,
    "m2/kmol": units.m**2/units.kmol,
    "N/m": units.N/units.m,
    "Coulomb.m": None, # Not supported in pyomo
    "kg/kmol": units.kg/units.kmol,
    "kmol/m3": units.kmol/units.m**3,
    "J0.5/m1.5": units.J**0.5/units.m**1.5,
}

def str_to_units(string: str) -> units:
    if string is not None:
        try:
            return str_map[string]
        except:
            return f"Unit {string} not found"
    else:
        return None

def str_to_float(string: str) -> float | str:
    try:
        return float(string)
    except:
        return string

# Helper functions
def parse_element(elem: ET.Element):
    try:
        return (str_to_float(elem.get('value')), str_to_units(elem.get('units')))
    except:
        return None

Coefficients = Dict[str, float | str]
def parse_coeff(element: ET.Element) -> Coefficients:
    self = {"units": str_to_units(element.get('units'))}
    for child in element:
        value = child.get('value')
        if value is not None:
            self[child.tag] = (str_to_float(value), None)
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
