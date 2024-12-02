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

def test_heater():
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