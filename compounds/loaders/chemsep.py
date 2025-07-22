from compounds.loaders import loader
import xml.etree.ElementTree as ET
from pydantic import BaseModel
from typing import Dict
from compounds.CompoundDB import PropertyPackage
import os

print("Chemsep loader initialized.")

@loader("chemsep")
def load(registry):

    registry.register_package(PropertyPackage("peng-robinson"))

    print("RARRRRR")

    for file in os.listdir(os.path.dirname(__file__) + "/data/chemsep/"):
        if file.endswith(".xml"):
            compound_name = file[:-4]
            compound = load_compound(compound_name)
            print(compound_name)
            registry.register_compound(compound_name, "chemsep", compound)


def convert_string_to_float(string: str) -> float | str:
    try:
        return float(string)
    except:
        return string


class UnitValuePair(BaseModel):
    name: str | None
    value: float | str | None
    unit: str | None


# Helper functions
def parse_element(elem: ET.Element):
    try:
        return UnitValuePair(name=elem.get('name'), value=convert_string_to_float(elem.get('value')), unit=elem.get('units'))
    except:
        return None


Coefficients = Dict[str, float | str]


def parse_coeff(element: ET.Element) -> Coefficients:
    self = {}
    for child in element:
        value = child.get('value')
        if value is not None:
            self[child.tag] = convert_string_to_float(value)
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
        name (str): The XML data as a string.
    """
    self = {}
    with open(os.path.dirname(__file__) + "/data/chemsep/" + name.lower() + ".xml", 'r') as file:
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
