from ..Compound import load_compound


def test_parsing():
  comp = load_compound("benzene")
  assert comp["AbsEntropy"].value == 269300
  assert comp["AbsEntropy"].unit == "J/kmol/K"