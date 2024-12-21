# Build and solve a heater block.
from ..build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units, Var

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


def report_vars(blk, file):
    with open(file, "w") as f:
        for c in blk.component_data_objects(Var, descend_into=True, active=True):
            f.write(f"{c}: {c.value} {"fixed" if c.fixed else "unfixed"}\n")
    with open("file14.txt", "r") as f:
        with open(f"{file}_diff", "w") as f2:
            for c in blk.component_data_objects(Var, descend_into=True):
                line = f.readline()
                key, value, _ = line.split(" ")
                if type(c.value) not in [float, int]:
                    diff = "None - " + type(c.value).__name__
                elif value == "None":
                    diff = "var_was_none"
                else:
                    diff = abs(float(value) - c.value)
                f2.write(f"{c}: {diff} {"fixed" if c.fixed else "unfixed"}\n")

def add_report_vars_to_blk(blk):
    blk.report_vars = lambda x: report_vars(blk, x)

def test_compressor_asu():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("peng-robinson", ["argon", "oxygen", "nitrogen"], ["Liq", "Vap"])
    # m.fs.properties.default_temperature_bounds = (m.fs.properties.default_temperature_bounds[0], None)
    
    m.fs.compressor = PressureChanger(
        property_package=m.fs.properties, 
        thermodynamic_assumption=ThermodynamicAssumption.isentropic,
        compressor=True)

    sb = _get_state_from_port(m.fs.compressor.inlet,0)
    m.fs.compressor.inlet.flow_mol[0].fix(1)
    m.fs.compressor.inlet.pressure.fix(200000 * units.Pa) # doesn't work from 100,000
    # sb.temperature.fix((273.15+25) * units.K) # Temperature is not on port
    sb.enth_mol.fix(4000)
    sb.temperature.value = 5
    m.fs.compressor.inlet.mole_frac_comp[0, "argon"].fix(0.33)
    m.fs.compressor.inlet.mole_frac_comp[0, "oxygen"].fix(0.33)
    m.fs.compressor.inlet.mole_frac_comp[0, "nitrogen"].fix(0.33)
    # m.fs.compressor.outlet.enth_mol.fix(-20)

    m.fs.compressor.outlet.pressure.fix(200000*units.Pa)
    m.fs.compressor.efficiency_isentropic.fix(0.75)
    # m.fs.compressor.control_volume.properties_out[0].temperature.ub = None

    # assert degrees_of_freedom(m) == 0

    add_report_vars_to_blk(m.fs.compressor.control_volume.properties_out[0])
    add_report_vars_to_blk(m.fs.compressor.control_volume.properties_in[0])

    m.fs.compressor.control_volume.properties_in.initialize()
    m.fs.compressor.control_volume.properties_out[0].report_vars("file15.txt")
    sb.report_vars("file14.txt")
    # m.fs.compressor.outlet.flow_mol[0].value = 1*units.kilomol/units.hour
    # m.fs.compressor.control_volume.properties_out[0].temperature.value = 1357
    # m.fs.compressor.control_volume.properties_out[0].report_vars("file16.txt")
    m.fs.compressor.outlet.flow_mol[0].fix(100)
    m.fs.compressor.outlet.pressure.fix(200000 * units.Pa) # doesn't work from 100,000
    # sb.temperature.fix((273.15+25) * units.K) # Temperature is not on port
    m.fs.compressor.outlet.enth_mol.fix(5)
    m.fs.compressor.outlet.mole_frac_comp[0, "argon"].fix(0.33333)
    m.fs.compressor.outlet.mole_frac_comp[0, "oxygen"].fix(0.33333)
    m.fs.compressor.outlet.mole_frac_comp[0, "nitrogen"].fix(0.33333)

    # m.fs.compressor.control_volume.properties_out.initialize()
    return

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
    
    sb = _get_state_from_port(m.fs.compressor.inlet,0)
    m.fs.compressor.inlet.flow_mol[0].fix(1*units.kilomol/units.hour)
    m.fs.compressor.inlet.pressure.fix(1000000 * units.Pa)
    sb.temperature.fix((273.15+25) * units.K)
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