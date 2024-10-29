# Simple test that loading and accessing compounds works

# FYI: How to make a singleton in python:
# https://stackoverflow.com/questions/6760685/what-is-the-best-way-of-implementing-singleton-in-python
# TLDR; don't, just use a module


def test_load_compounds():
    from compounds import CompoundDB as db
    assert "Benzene" in db.get_compound_names()
    benzene = db.get_compound("Benzene")
    assert benzene["AbsEntropy"].value == 269300
