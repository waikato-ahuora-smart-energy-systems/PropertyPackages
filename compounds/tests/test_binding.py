from compounds.CompoundDB import db


def test_binding():
  comp = db.get_compound("benzene").get_source("chemsep")
  assert comp["AbsEntropy"].value == 269300
  assert comp["AbsEntropy"].unit == "J/kmol/K"