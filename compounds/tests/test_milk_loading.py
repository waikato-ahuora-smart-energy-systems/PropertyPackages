
def test_milk_compounds():
    from compounds import CompoundDB as db
    assert "milk_solids" in db.get_compound_names()
    assert "water" in db.get_compound_names()
    milk_solids = db.get_compound("milk_solids")
    assert milk_solids["AbsEntropy"].value == "not for use"