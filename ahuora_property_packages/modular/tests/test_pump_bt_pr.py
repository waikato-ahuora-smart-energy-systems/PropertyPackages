# Build and solve a heater block.
from ahuora_property_packages.build_package import build_package
from pytest import approx
# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.pressure_changer import Pump

def assert_approx(value, expected_value, error_margin):
    percent_error = error_margin / 100
    tolerance = abs(percent_error * expected_value)
    assert approx(value, abs=tolerance) == expected_value

def test_pump():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])
    
    m.fs.pump = Pump(property_package=m.fs.properties)

    m.fs.pump.inlet.flow_mol[0].fix(100)
    m.fs.pump.inlet.pressure.fix(101325)
    m.fs.pump.inlet.mole_frac_comp[0, "benzene"].fix(0.4)
    m.fs.pump.inlet.mole_frac_comp[0, "toluene"].fix(0.6)
    m.fs.pump.inlet.temperature.fix(353)

    m.fs.pump.deltaP.fix(100000)
    m.fs.pump.efficiency_pump.fix(0.8)

    m.fs.pump.initialize()

    assert degrees_of_freedom(m) == 0

    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)

    assert_approx(value(m.fs.pump.outlet.temperature[0]), 353, 0.2)
    assert value(m.fs.pump.outlet.pressure[0]) == approx(201325)
    assert value(m.fs.pump.outlet.flow_mol[0]) == approx(100)