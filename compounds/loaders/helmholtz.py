from compounds.loaders import loader
from compounds.PropertyPackage import DefaultPropertyPackage

@loader("helmholtz")
def load(registry):
    
    registry.register_package(DefaultPropertyPackage("helmholtz"))

    # todo: add compounds and binds