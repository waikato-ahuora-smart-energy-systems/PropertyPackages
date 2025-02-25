from property_packages.milk import build_milk_package
# Build and solve a state block.
from property_packages.build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, value

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver

from pyomo.util.check_units import assert_units_consistent

from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from pyomo.environ import  check_optimal_termination, ConcreteModel, Objective
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom
solver = get_solver(solver="ipopt")

from pyomo.environ import units

import idaes.logger as idaeslog
SOUT = idaeslog.INFO


def test_milk_1():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_milk_package()

    m.fs.state = m.fs.properties.build_state_block([1], defined_state=True)

    iscale.calculate_scaling_factors(m)

    assert degrees_of_freedom(m) == 5

    m.fs.state[1].flow_mol.fix(1*units.mol/units.s)
    m.fs.state[1].pressure.fix(101325*units.Pa)
    m.fs.state[1].temperature.fix(300*units.K)
    m.fs.state[1].mole_frac_comp["water"].fix(0.5)
    m.fs.state[1].mole_frac_comp["milk_solid"].fix(0.5)

    assert degrees_of_freedom(m) == 0

    m.fs.state.initialize()