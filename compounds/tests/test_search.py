from compounds.CompoundDB import db

def test_compound_search():
    assert len(db.search_compounds("ane")) == 122 # Subject to change
    assert len(db.search_compounds("ide")) == 55 # Subject to change

def test_compound_search_filter():
    assert len(db.search_compounds("ide", package_filters=["biomass_combustion_reaction"])) == 2
    assert len(db.search_compounds("milk", package_filters=["peng-robinson"])) == 0
    assert len(db.search_compounds("milk", package_filters=["milk"])) == 1

def test_compound_search_filter_strict():
    assert db.search_compounds("ide", package_filters=["biomass_combustion_reaction", "peng-robinson"], 
                               filter_strict=True) == {"carbon monoxide", "carbon dioxide"}