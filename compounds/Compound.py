from pprint import pprint
import xml.etree.ElementTree as ET
from typing import Dict
import os
import json

        
class Coefficients:
    # Allows square bracket access
    def __getitem__(self, name: str) -> any:
        return self._attributes[name]

    def __setitem__(self, name: str, value: any) -> None:
        self._attributes[name] = value
    
    def __init__(self, element: ET.Element):
        self._attributes = {}
        for child in element:
            if child.get('value') is not None:
                self[child.tag] = convert_string_to_float(child.get('value'))
        

class UnitValuePair:
    def __init__(self, name, value, unit):
        self.name = name
        self.value = convert_string_to_float(value)
        self.unit = unit

    def to_dict(self) -> Dict[str, any]:
        return {
            'name': self.name,
            'value': self.value,
            'unit': self.unit
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict())

# Helper functions
def parse_element(elem: ET.Element):
    try:
        return UnitValuePair(elem.get('name'), elem.get('value'), elem.get('units'))
    except:
        return None

def parse_coeff(element: ET.Element) -> Coefficients:
    return Coefficients(element)

def convert_string_to_float(string: str) -> float:
    try:
        return float(string)
    except:
        return string

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

class Compound:

    def __init__(self, name: str) -> None:
        """
        Initializes the Compound object.

        Args:
            name (str): Name of compound.
        """

        self._attributes = {}
        self._parse_xml(name)

    # Allows square bracket access
    def __getitem__(self, name: str) -> any:
        return self._attributes[name]

    def __setitem__(self, name: str, value: any) -> None:
        self._attributes[name] = value

    def _parse_xml(self, name: str):
        """
        Parses the XML string and populates the object's attributes.

        Args:
            xml_string (str): The XML data as a string.
        """

        with open(os.path.dirname(__file__) + "/data_files/" + name + ".xml", 'r') as file:
            file = ''.join(file.readlines())

        root = ET.fromstring(file)

        # Iterating over all attributes in the template
        # and parsing the corresponding XML elements
        for attr in compound_template:
            if root.find(attr) is not None:
                self[attr] = compound_template[attr](root.find(attr))
            else:
                self[attr] = None


