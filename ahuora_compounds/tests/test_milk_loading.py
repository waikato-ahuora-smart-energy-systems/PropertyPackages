
def test_milk_compounds():
    from ahuora_compounds import CompoundDB as db
    assert "milk_solid" in db.get_compound_names()
    assert "water" in db.get_compound_names()
    milk_solid = db.get_compound("milk_solid")
    assert milk_solid == None # no data related to milk_solid