from ..helmholtz_builder import build_helmholtz_package
from pytest import approx
from pyomo.environ import ConcreteModel, SolverFactory, value, assert_optimal_termination
from pyomo.core.base.constraint import Constraint, ScalarConstraint
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Compressor
from idaes.core.util.model_statistics import degrees_of_freedom

def flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_helmholtz_package(["h2o"])
    return m


def build_state(m):
    m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)
    return m.fs.state[0]


def initialize(m):
    assert degrees_of_freedom(m) == 0
    m.fs.state.initialize()
    assert degrees_of_freedom(m) == 0


def solve(m):
    solver = SolverFactory('ipopt')
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)


def test_constrain():
    m = flowsheet()
    sb = build_state(m)
    # fix flow_mol directly
    c = sb.constrain("flow_mol", 1)
    assert c == sb.flow_mol
    assert c.value == 1
    assert c.is_fixed()

    # add a constraint for flow_mass
    c = sb.constrain("flow_mass", 1)
    assert type(c) == ScalarConstraint
    assert c in sb.component_data_objects(Constraint)
    assert getattr(sb.constraints, "flow_mass") == c


def test_state_definition():
    m = flowsheet()
    sb = build_state(m)
    sb.constrain("enth_mass", 1878.87)
    sb.constrain("pressure", 101325)
    sb.constrain("flow_mass", 1)
    initialize(m)
    solve(m)
    assert value(sb.temperature) == approx(273.5809)


def test_state_definition_temp():
    m = flowsheet()
    sb = build_state(m)
    sb.constrain("temperature", 273.5809)
    sb.constrain("pressure", 101325)
    sb.constrain("flow_mass", 1)
    initialize(m)
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
    inlet.constrain("temperature", 273.5809)
    inlet.constrain("pressure", 101325)
    inlet.constrain("flow_mass", 1)
    m.fs.compressor.efficiency_isentropic.fix(1)
    m.fs.compressor.deltaP.fix(100000)
    m.fs.compressor.initialize(outlvl=1)
    solve(m)
    assert value(outlet.temperature) == approx(273.58047830928655)


def test_butane():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_helmholtz_package(["i-butane"])
    m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)
    state = m.fs.state[0]
    state.constrain("temperature",262)
    state.constrain("pressure",101325)
    state.constrain("flow_mol",100)
    m.fs.state.initialize()
    solve(m)
    # boiling point of butane at 1atm is -0.5 degC.
    assert value(state.vapor_frac) == approx(1)
    state.constraints.del_component(state.constraints.temperature)
    state.constrain("temperature",260)
    solve(m)
    assert value(state.vapor_frac) == approx(0)
    