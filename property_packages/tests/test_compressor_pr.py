# Build and solve a heater block.
from ..build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.heater import Heater
from idaes.models.unit_models.pressure_changer import Pump
from idaes.core.util.tables import _get_state_from_port
from idaes.models.unit_models.pressure_changer import PressureChanger, ThermodynamicAssumption

def assert_approx(value, expected_value, error_margin):
    percent_error = error_margin / 100
    tolerance = abs(percent_error * expected_value)
    assert approx(value, abs=tolerance) == expected_value

def test_compressor_asu():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("peng-robinson", ["argon", "oxygen", "nitrogen"], ["Liq", "Vap"])
    
    m.fs.compressor = PressureChanger(
        property_package=m.fs.properties, 
        thermodynamic_assumption=ThermodynamicAssumption.isentropic,
        compressor=True)

    m.fs.compressor.inlet.flow_mol[0].fix(1*units.kilomol/units.hour)
    m.fs.compressor.inlet.pressure.fix(200000 * units.Pa) # doesn't work from 100,000
    m.fs.compressor.inlet.temperature.fix((273.15+25) * units.K)
    m.fs.compressor.inlet.mole_frac_comp[0, "argon"].fix(0.33)
    m.fs.compressor.inlet.mole_frac_comp[0, "oxygen"].fix(0.33)
    m.fs.compressor.inlet.mole_frac_comp[0, "nitrogen"].fix(0.33)

    m.fs.compressor.outlet.pressure.fix(500000*units.Pa)
    m.fs.compressor.efficiency_isentropic.fix(0.75)

    assert degrees_of_freedom(m) == 0

    m.fs.compressor.initialize()

    assert degrees_of_freedom(m) == 0

    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)

    assert_approx(value(m.fs.compressor.outlet.temperature[0]), 466.12, 0.1)

def test_expander_asu():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("peng-robinson", ["argon", "oxygen", "nitrogen"], ["Liq", "Vap"])
    
    m.fs.compressor = PressureChanger(
        property_package=m.fs.properties, 
        thermodynamic_assumption=ThermodynamicAssumption.adiabatic,
        compressor=False)

    m.fs.compressor.inlet.flow_mol[0].fix(1*units.kilomol/units.hour)
    m.fs.compressor.inlet.pressure.fix(1000000 * units.Pa)
    m.fs.compressor.inlet.temperature.fix((273.15+25) * units.K)
    m.fs.compressor.inlet.mole_frac_comp[0, "argon"].fix(0.33)
    m.fs.compressor.inlet.mole_frac_comp[0, "oxygen"].fix(0.33)
    m.fs.compressor.inlet.mole_frac_comp[0, "nitrogen"].fix(0.33)

    m.fs.compressor.outlet.pressure.fix(100000*units.Pa)

    assert degrees_of_freedom(m) == 0

    m.fs.compressor.initialize()

    assert degrees_of_freedom(m) == 0

    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)

    assert_approx(int(value(m.fs.compressor.deltaP[0])), -900000, 0)