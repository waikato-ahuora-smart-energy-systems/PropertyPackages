from compounds.CompoundDB import get_compound
from property_packages.types import States
from typing import List
from pyomo.environ import units as u
from ..build_package import build_package

def test_units():

    # Build list of compound objects
    build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])