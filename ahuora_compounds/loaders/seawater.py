from ahuora_compounds.loaders import loader
from ahuora_compounds.PropertyPackage import DefaultPropertyPackage

@loader("seawater")
def load(registry):

    registry.register_package(DefaultPropertyPackage("seawater"))
    
    registry.register_compound('H2O', "seawater", {})
    registry.register_compound('TDS', "seawater", {})

    registry.bind("TDS", "seawater")
    registry.bind("H2O", "seawater")