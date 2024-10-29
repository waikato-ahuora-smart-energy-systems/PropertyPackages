from ..Compound import Compound


def test_parsing():
  comp = Compound("benzene")
  assert comp["AbsEntropy"].value == "269300"
  assert comp["AbsEntropy"].unit == "J/kmol/K"