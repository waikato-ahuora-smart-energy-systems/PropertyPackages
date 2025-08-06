from compounds.loaders import loader
from compounds.PropertyPackage import DefaultPropertyPackage

@loader("biomass")
def load(registry):

    registry.register_package(DefaultPropertyPackage("biomass_and_flue"))
    registry.register_package(DefaultPropertyPackage("biomass_combustion_reaction"))
    
    registry.register_compound('water', "biomass", {})
    registry.register_compound('carbon dioxide', "biomass", {})
    registry.register_compound('oxygen', "biomass", {})
    registry.register_compound('carbon monoxide', "biomass", {})
    registry.register_compound('nitrogen', "biomass", {})
    registry.register_compound('biomass', "biomass", {})
    registry.register_compound('ash', "biomass", {})

    registry.bind("biomass", "biomass_and_flue")
    registry.bind("water", "biomass_and_flue")
    registry.bind("carbon dioxide", "biomass_and_flue")
    registry.bind("oxygen", "biomass_and_flue")
    registry.bind("carbon monoxide", "biomass_and_flue")
    registry.bind("nitrogen", "biomass_and_flue")
    registry.bind("ash", "biomass_and_flue")

    registry.bind("biomass", "biomass_combustion_reaction")
    registry.bind("water", "biomass_combustion_reaction")
    registry.bind("carbon dioxide", "biomass_combustion_reaction")
    registry.bind("oxygen", "biomass_combustion_reaction")
    registry.bind("carbon monoxide", "biomass_combustion_reaction")
    registry.bind("nitrogen", "biomass_combustion_reaction")
    registry.bind("ash", "biomass_combustion_reaction")