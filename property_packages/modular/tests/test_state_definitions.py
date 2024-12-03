from ..template_builder import build_config
from pytest import approx
from pyomo.environ import ConcreteModel, SolverFactory, value, units
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Compressor

def flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_config("peng-robinson",["benzene","toluene"],["Liq","Vap"])
    return m

def solve(m):
    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)

def test_state_definition2():
    m = flowsheet()
    m.fs.sb = m.fs.properties.build_state_block(m.fs.time)
    sb = m.fs.sb[0]
    sb.constrain("temperature", 280)
    sb.constrain("pressure", 101325)
    sb.constrain("flow_mass", 1)
    sb.mole_frac_comp["benzene"].fix(0.5)
    sb.mole_frac_comp["toluene"].fix(0.5)
    m.fs.sb.initialize()
    solve(m)
    assert value(sb.temperature) == approx(280)


def test_state_definition_temp():
    m = flowsheet()
    m.fs.sb = m.fs.properties.build_state_block(m.fs.time)
    sb = m.fs.sb[0]
    sb.constrain("smooth_temperature", 273.5809)
    sb.constrain("pressure", 101325)
    sb.constrain("flow_mass", 1)
    sb.mole_frac_comp["benzene"].fix(0.5)
    sb.mole_frac_comp["toluene"].fix(0.5)
    m.fs.sb.initialize()
    solve(m)
    assert value(sb.enth_mass) == approx(1878.712)


def test_initialise_compressor():
    # The purpose of this test is to make sure the compressor can initialise.
    # We have had some errors with Too Few Degrees of Freedom
    # being thrown in initialisation. This is because constraints are 
    # being fixed in the inlet state block instead of variables.
    # This asserts that the appropriate variables are unfixed.
    m = flowsheet()
    m.fs.compressor = Compressor(property_package=m.fs.properties)
    inlet = m.fs.compressor.control_volume.properties_in[0]
    outlet = m.fs.compressor.control_volume.properties_out[0]
    inlet.constrain("smooth_temperature", 273.5809)
    inlet.constrain("pressure", 101325)
    inlet.constrain("flow_mass", 1)
    inlet.mole_frac_comp["benzene"].fix(0.5)
    inlet.mole_frac_comp["toluene"].fix(0.5)
    m.fs.compressor.deltaP.fix(100000)
    m.fs.compressor.initialize()
    solve(m)
    assert value(outlet.temperature) == approx(393.5689573)

    
