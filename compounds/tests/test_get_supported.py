from compounds.CompoundDB import db

def test_combustion_support_strict():
    assert db.get_supported_compounds(["biomass_combustion_reaction", "biomass_and_flue"], 
                                    strict=True) == {"biomass", "water", "carbon dioxide", "oxygen", "carbon monoxide", "nitrogen"}

def test_combustion_support_not_strict():
    assert db.get_supported_compounds(["biomass_combustion_reaction", "biomass_and_flue"], 
                                    strict=False) == {"ash", "biomass", "water", "carbon dioxide", "oxygen", "carbon monoxide", "nitrogen"}
    
def test_package_support_strict():
    assert db.get_supported_packages(["biomass", "water", "carbon dioxide", "oxygen", "carbon monoxide", "nitrogen"], 
                                    strict=True) == {"biomass_combustion_reaction", "biomass_and_flue"}

def test_package_support_not_strict():
    assert db.get_supported_packages(["biomass", "water", "carbon dioxide", "oxygen", "carbon monoxide", "nitrogen", "ash"], 
                                    strict=False) == {"biomass_combustion_reaction", "biomass_and_flue", "peng-robinson"}