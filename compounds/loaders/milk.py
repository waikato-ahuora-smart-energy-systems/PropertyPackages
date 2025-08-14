from compounds.loaders import loader
from compounds.PropertyPackage import DefaultPropertyPackage

@loader("milk")
def load(registry):

    registry.register_package(DefaultPropertyPackage("milk"))
    
    registry.register_compound('milk_solid', "milk", {})
    registry.register_compound('water', "milk", {})

    registry.bind("milk_solid", "milk")
    registry.bind("water", "milk")