# Build and solve a heater block.
from ahuora_property_packages.build_package import build_package
from pytest import approx
# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.heat_exchanger import HeatExchanger, delta_temperature_amtd_callback
from idaes.models.properties.iapws95 import htpx

def assert_approx(value, expected_value, error_margin):
    percent_error = error_margin / 100
    tolerance = abs(percent_error * expected_value)
    assert approx(value, abs=tolerance) == expected_value

def test_heat_exchanger_bt():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties1 = build_package("helmholtz", ["h2o"], ["Liq", "Vap"])
    m.fs.properties2 = build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])
    
    m.fs.heat_exchanger = HeatExchanger(
        delta_temperature_callback=delta_temperature_amtd_callback,
        hot_side_name="shell",
        cold_side_name="tube",
        shell={"property_package": m.fs.properties1},
        tube={"property_package": m.fs.properties2})
    
    assert degrees_of_freedom(m) == 10

    h = htpx(450*units.K, P = 101325*units.Pa) # 450 K
    m.fs.heat_exchanger.shell_inlet.flow_mol.fix(100) # mol/s
    m.fs.heat_exchanger.shell_inlet.pressure.fix(101325) # Pa
    m.fs.heat_exchanger.shell_inlet.enth_mol.fix(h) # J/mol

    assert degrees_of_freedom(m) == 7

    m.fs.heat_exchanger.tube_inlet.flow_mol.fix(250) # mol/s
    m.fs.heat_exchanger.tube_inlet.mole_frac_comp[0, "benzene"].fix(0.4)
    m.fs.heat_exchanger.tube_inlet.mole_frac_comp[0, "toluene"].fix(0.6)
    m.fs.heat_exchanger.tube_inlet.pressure.fix(101325) # Pa
    m.fs.heat_exchanger.tube_inlet.temperature[0].fix(350) # K

    assert degrees_of_freedom(m) == 2

    m.fs.heat_exchanger.area.fix(50) # m2
    m.fs.heat_exchanger.overall_heat_transfer_coefficient[0].fix(500) # W/m2/K

    assert degrees_of_freedom(m) == 0

    m.fs.heat_exchanger.initialize()

    solver = SolverFactory('ipopt')
    result = solver.solve(m)

    assert value(m.fs.heat_exchanger.shell.properties_out[0].temperature) == approx(373.162, abs=1e-2)
    assert_approx(value(m.fs.heat_exchanger.tube.properties_out[0].temperature), 369.24, 0.5)

def test_heat_exchanger_asu():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties1 = build_package("helmholtz", ["h2o"], ["Liq", "Vap"])
    m.fs.properties2 = build_package("peng-robinson", ["oxygen", "argon", "nitrogen"], ["Liq", "Vap"])
    
    m.fs.heat_exchanger = HeatExchanger(
        delta_temperature_callback=delta_temperature_amtd_callback,
        hot_side_name="shell",
        cold_side_name="tube",
        shell={"property_package": m.fs.properties1},
        tube={"property_package": m.fs.properties2})
    
    assert degrees_of_freedom(m) == 11

    h = htpx(473.15*units.K, P = 100000*units.Pa) # 200 C
    m.fs.heat_exchanger.shell_inlet.flow_mol.fix(20*units.kilomol/units.hour)
    m.fs.heat_exchanger.shell_inlet.pressure.fix(100000) # Pa
    m.fs.heat_exchanger.shell_inlet.enth_mol.fix(h) # J/mol

    assert degrees_of_freedom(m) == 8

    m.fs.heat_exchanger.tube_inlet.flow_mol.fix(1*units.kilomol/units.hour) # mol/s
    m.fs.heat_exchanger.tube_inlet.mole_frac_comp[0, "oxygen"].fix(1/3)
    m.fs.heat_exchanger.tube_inlet.mole_frac_comp[0, "argon"].fix(1/3)
    m.fs.heat_exchanger.tube_inlet.mole_frac_comp[0, "nitrogen"].fix(1/3)
    m.fs.heat_exchanger.tube_inlet.pressure.fix(100000) # Pa
    m.fs.heat_exchanger.tube_inlet.temperature[0].fix((25 + 273.15)*units.K)

    assert degrees_of_freedom(m) == 2

    m.fs.heat_exchanger.shell_outlet.pressure.fix(100000)
    m.fs.heat_exchanger.heat_duty.fix(927.7) # W

    assert degrees_of_freedom(m) == 0

    m.fs.heat_exchanger.initialize()

    solver = SolverFactory('ipopt')
    result = solver.solve(m)
    
    assert_approx(value(m.fs.heat_exchanger.shell.properties_out[0].temperature), 195.3 + 273.15, 0.1)
