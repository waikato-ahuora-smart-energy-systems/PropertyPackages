from ahuora_compounds.loaders import loader
from ahuora_compounds.PropertyPackage import DefaultPropertyPackage

@loader("humid_air")
def load(registry):

    registry.register_package(DefaultPropertyPackage("humid_air"))
    
    registry.register_compound('air', "humid_air", {})
    registry.register_compound('water', "humid_air", {})

    registry.bind("air", "humid_air")
    registry.bind("water", "humid_air")