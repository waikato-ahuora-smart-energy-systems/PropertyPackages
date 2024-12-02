# Build and solve a heater block.
from ..build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.heater import Heater
from idaes.core.util.tables import _get_state_from_port

def assert_approx(value, expected_value, error_margin):
    percent_error = error_margin / 100
    tolerance = abs(percent_error * expected_value)
    assert approx(value, abs=tolerance) == expected_value

def test_heater_bt():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])
    
    m.fs.heater = Heater(property_package=m.fs.properties)

    m.fs.heater.inlet.flow_mol.fix(1000/3600)
    m.fs.heater.inlet.mole_frac_comp[0, "benzene"].fix(0.4)
    m.fs.heater.inlet.mole_frac_comp[0, "toluene"].fix(0.6)
    m.fs.heater.inlet.pressure.fix(101325)
    m.fs.heater.inlet.temperature.fix(353)
    m.fs.heater.heat_duty.fix(459.10147722222354)

    m.fs.heater.initialize()

    assert degrees_of_freedom(m) == 0

    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)

    assert_approx(value(_get_state_from_port(m.fs.heater.outlet,0).temperature), 363, 0.2)
    assert value(m.fs.heater.outlet.pressure[0]) == approx(101325)
    assert value(m.fs.heater.outlet.flow_mol[0]) == approx(1000/3600)

def test_heater_asu():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("peng-robinson", ["argon", "nitrogen", "oxygen"], ["Liq", "Vap"])
    
    m.fs.heater = Heater(property_package=m.fs.properties)

    m.fs.heater.inlet.flow_mol.fix(units.convert(1 * units.kilomol/units.hour, units.mol/units.s)) # mol/s
    m.fs.heater.inlet.mole_frac_comp[0, "argon"].fix(0.33)
    m.fs.heater.inlet.mole_frac_comp[0, "nitrogen"].fix(0.33)
    m.fs.heater.inlet.mole_frac_comp[0, "oxygen"].fix(0.33)
    m.fs.heater.inlet.pressure.fix(100000)
    m.fs.heater.inlet.temperature.fix(units.convert_temp_C_to_K(25))
    m.fs.heater.outlet.temperature.fix(units.convert_temp_C_to_K(50))

    m.fs.heater.initialize()

    assert degrees_of_freedom(m) == 0

    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)

    assert_approx(value(m.fs.heater.heat_duty[0]), 183.856, 1) # 1% tolerance

    m.fs.heater2 = Heater(property_package=m.fs.properties)

    m.fs.heater2.inlet.flow_mol.fix(units.convert(1 * units.kilomol/units.hour, units.mol/units.s)) # mol/s
    m.fs.heater2.inlet.mole_frac_comp[0, "argon"].fix(0.33)
    m.fs.heater2.inlet.mole_frac_comp[0, "nitrogen"].fix(0.33)
    m.fs.heater2.inlet.mole_frac_comp[0, "oxygen"].fix(0.33)
    m.fs.heater2.inlet.pressure.fix(100000)
    m.fs.heater2.inlet.temperature.fix(units.convert_temp_C_to_K(25))
    m.fs.heater2.outlet.temperature.fix(units.convert_temp_C_to_K(-15))

    m.fs.heater2.initialize()

    assert degrees_of_freedom(m) == 0

    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)

    assert_approx(value(m.fs.heater2.heat_duty[0]), -292.5, 1) # 1% tolerance
