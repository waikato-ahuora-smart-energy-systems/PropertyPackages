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

def test_milk_invalid_configuration():
    try:
        build_package("milk", ["benzene", "milk_solid"], ["Liq", "Vap"])
    except ValueError:
        assert True
        return
    assert False

def test_milk_valid_configuration():
    try:
        build_package("milk", ["water", "milk_solid"], ["Liq", "Vap"])
    except ValueError:
        assert False
    assert True

def test_milk_solve():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_package("milk", ["water", "milk_solid"], ["Liq", "Vap"])

    m.fs.state = m.fs.properties.build_state_block([1], defined_state=True)

    iscale.calculate_scaling_factors(m)

    assert degrees_of_freedom(m) == 5

    m.fs.state[1].flow_mol.fix(1*units.mol/units.s)
    m.fs.state[1].pressure.fix(101325*units.Pa)
    m.fs.state[1].temperature.fix(300*units.K)

    print(dir(m.fs.state[1].mole_frac_comp))

    m.fs.state[1].mole_frac_comp["water"].fix(0.9)
    m.fs.state[1].mole_frac_comp["milk_solid"].fix(0.1)

    assert degrees_of_freedom(m) == 0

    m.fs.state.initialize()

def test_milk_custom_props():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_package("milk", ["water", "milk_solid"], ["Liq", "Vap"])
    m.fs.state = m.fs.properties.build_state_block([1], defined_state=True)
    m.fs.state[1].flow_mol.fix(1*units.mol/units.s)
    m.fs.state[1].pressure.fix(101325*units.Pa)
    m.fs.state[1].temperature.fix(300*units.K)
    m.fs.state[1].mole_frac_comp["water"].fix(0.9)
    m.fs.state[1].mole_frac_comp["milk_solid"].fix(0.1)
    m.fs.state.initialize()

    m.fs.state[1].enth_mol # trigger enthalpy build
    m.fs.state[1].entr_mol # trigger entropy build

    solver.solve(m)

    # Values not verified

    assert value(m.fs.state[1].enth_mol_comp["water"]) == approx(-527400, abs=1e3)
    assert value(m.fs.state[1].enth_mol_comp["milk_solid"]) == approx(-1453000, abs=1e3)
    assert value(m.fs.state[1].enth_mol) == approx(-1.98e6, abs=1e4)

    assert value(m.fs.state[1].enth_mass_comp["water"]) == approx(-2.92e7, abs=1e5)
    assert value(m.fs.state[1].enth_mass_comp["milk_solid"]) == approx(-6.26e6, abs=1e4)
    assert value(m.fs.state[1].enth_mass) == approx(-3.55e7, abs=1e5)

    assert value(m.fs.state[1].flow_mol_comp["water"]) == 0.9
    assert value(m.fs.state[1].flow_mol_comp["milk_solid"]) == 0.1

    assert value(m.fs.state[1].flow_mass_comp["water"]) == approx(0.0162, abs=1e-4)
    assert value(m.fs.state[1].flow_mass_comp["milk_solid"]) == approx(0.0232, abs=1e-4)
    assert value(m.fs.state[1].flow_mass) == approx(0.0394, abs=1e-4)

    assert value(m.fs.state[1].entr_mol_comp["water"]) == approx(259, rel=1e-2)
    assert value(m.fs.state[1].entr_mol_comp["milk_solid"]) == approx(171, rel=1e-2)
    assert value(m.fs.state[1].entr_mol) == approx(430, rel=1e-2)

    assert value(m.fs.state[1].entr_mass_comp["water"]) == approx(1.44e4, abs=1e2)
    assert value(m.fs.state[1].entr_mass_comp["milk_solid"]) == approx(739, abs=1)
    assert value(m.fs.state[1].entr_mass) == approx(1.51e4, abs=1e2)

    assert value(m.fs.state[1].vapor_frac) == approx(0, abs=1e-4)
    assert value(m.fs.state[1].flow_vol) == approx(0)
    assert value(m.fs.state[1].total_energy_flow) == approx(-1400646, abs=1e3)





