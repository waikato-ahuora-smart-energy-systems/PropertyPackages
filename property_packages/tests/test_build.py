import pytest
from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    Set,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
import pyomo.common.unittest as unittest

from idaes.core import Component
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)
from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.state_definitions import FTPx

from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from ..build_package import build_package

def test_build():
    build_package("helmholtz", ["water"])
    build_package("peng-robinson", ["benzene", "toluene"])
    build_package("helmholtz", ["water"], ["Liq"])
    build_package("peng-robinson", ["benzene", "toluene"], ["Liq"])
    build_package("helmholtz", ["water"], ["Vap", "Liq"])
    build_package("peng-robinson", ["benzene", "toluene"], ["Vap", "Liq"])