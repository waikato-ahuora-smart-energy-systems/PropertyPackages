"""
Test initializing and solving a state block with constraints
rather than fixing the state variables directly.
"""

from pytest import approx
from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    value,
    units,
    assert_optimal_termination,
)
from pyomo.core.base.constraint import Constraint, ScalarConstraint
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from ..template_builder import build_config


def flowsheet():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_config(
        "peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"]
    )
    return m


def build_state(m):
    m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)
    return m.fs.state[0]


def initialize(m):
    assert degrees_of_freedom(m) == 0
    m.fs.state.initialize()
    assert degrees_of_freedom(m) == 0


def solve(m):
    solver = SolverFactory("ipopt")
    res = solver.solve(m, tee=True)
    assert_optimal_termination(res)


def test_constrain():
    m = flowsheet()
    sb = build_state(m)

    # fix flow_mol by string
    c = sb.constrain("flow_mol", 10)
    assert c == sb.flow_mol
    assert c.value == 10
    assert c.is_fixed()

    # constrain flow_mass by string
    c = sb.constrain("flow_mass", 10)
    assert type(c) == ScalarConstraint
    assert c in sb.component_data_objects(Constraint)
    assert getattr(sb.constraints, "flow_mass") == c


def test_constrain_component():
    m = flowsheet()
    sb = build_state(m)

    # fix flow_mol by component
    c = sb.constrain_component(sb.flow_mol, 10)
    assert c == sb.flow_mol
    assert c.value == 10
    assert c.is_fixed()

    # constrain flow_mass by component
    c = sb.constrain_component(sb.flow_mass, 10)
    assert type(c) == ScalarConstraint
    assert c in sb.component_data_objects(Constraint)
    assert getattr(sb.constraints, "flow_mass") == c


def test_state_vars():
    # default state vars: flow_mol, temperature, pressure, mole_frac_comp
    m = flowsheet()
    sb = build_state(m)
    sb.constrain("flow_mol", 1)
    sb.constrain("temperature", 290)
    sb.constrain("pressure", 100000)
    sb.mole_frac_comp["benzene"].fix(0.5)
    sb.mole_frac_comp["toluene"].fix(0.5)
    initialize(m)
    solve(m)
    assert value(sb.temperature) == approx(290)
    assert value(sb.pressure) == approx(100000)
    assert value(sb.flow_mol) == approx(1)
    assert value(sb.mole_frac_comp["benzene"]) == approx(0.5)
    assert value(sb.mole_frac_comp["toluene"]) == approx(0.5)


def test_enth_mol():
    m = flowsheet()
    sb = build_state(m)
    sb.constrain("flow_mol", 1)
    sb.constrain("enth_mol", 30611.284116732746)  # J/mol
    sb.constrain("pressure", 100000)
    sb.mole_frac_comp["benzene"].fix(0.5)
    sb.mole_frac_comp["toluene"].fix(0.5)
    initialize(m)
    solve(m)
    assert value(sb.temperature) == approx(290)
    assert value(sb.enth_mol) == approx(30611.284116732746)


def test_enth_mass():
    m = flowsheet()
    sb = build_state(m)
    sb.constrain("flow_mass", 1)
    sb.constrain("enth_mass", 359603.35471696255)
    sb.constrain("pressure", 101325)
    sb.mole_frac_comp["benzene"].fix(0.5)
    sb.mole_frac_comp["toluene"].fix(0.5)
    initialize(m)
    solve(m)
    assert value(sb.temperature) == approx(290, rel=1e-3)
    assert value(sb.enth_mass) == approx(359603.35471696255)


def test_entr_mol():
    m = flowsheet()
    sb = build_state(m)
    sb.constrain("flow_mass", 1)
    sb.constrain("entr_mol", -388.1271856851202)
    sb.constrain("pressure", 101325)
    sb.mole_frac_comp["benzene"].fix(0.5)
    sb.mole_frac_comp["toluene"].fix(0.5)
    initialize(m)
    solve(m)
    assert value(sb.temperature) == approx(290)
    assert value(sb.entr_mol) == approx(-388.1271856851202)


def test_entr_mass():
    m = flowsheet()
    sb = build_state(m)
    sb.constrain("flow_mass", 1)
    sb.constrain("entr_mass", -4559.489810913312)
    sb.constrain("pressure", 101325)
    sb.mole_frac_comp["benzene"].fix(0.5)
    sb.mole_frac_comp["toluene"].fix(0.5)
    initialize(m)
    solve(m)
    assert value(sb.temperature) == approx(290)
    assert value(sb.entr_mass) == approx(-4559.489810913312)


def test_flow_mass():
    m = flowsheet()
    sb = build_state(m)
    sb.constrain("flow_mass", 0.08512513500000002)
    sb.constrain("temperature", 290)
    sb.constrain("pressure", 101325)
    sb.mole_frac_comp["benzene"].fix(0.5)
    sb.mole_frac_comp["toluene"].fix(0.5)
    initialize(m)
    solve(m)
    assert value(sb.flow_mass) == approx(0.08512513500000002)
    assert value(sb.flow_mol) == approx(1)


def test_flow_vol():
    m = flowsheet()
    sb = build_state(m)
    sb.constrain("flow_vol", 9.660626354419042e-05)
    sb.constrain("temperature", 290)
    sb.constrain("pressure", 101325)
    sb.mole_frac_comp["benzene"].fix(0.5)
    sb.mole_frac_comp["toluene"].fix(0.5)
    initialize(m)
    solve(m)
    assert value(sb.flow_vol) == approx(9.660626354419042e-05)
    assert value(sb.flow_mol) == approx(1, rel=1e-3)


def test_enthalpy_temperature():
    # use enthalpy and temperature to define pressure
    # TODO: this test is not currently working
    # related: https://github.com/IDAES/idaes-pse/pull/1554
    m = flowsheet()
    sb = build_state(m)
    sb.constrain("flow_mass", 1)
    sb.constrain("temperature", 290)
    sb.constrain("enth_mol", 30611.284116732746)
    sb.mole_frac_comp["benzene"].fix(0.5)
    sb.mole_frac_comp["toluene"].fix(0.5)
    return  # fails to initialize
    initialize(m)
    solve(m)
    assert value(sb.pressure) == approx(100000)
