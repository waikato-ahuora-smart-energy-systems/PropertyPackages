import pytest
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Set,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.unittest import assertStructuredAlmostEqual

from idaes.core import Component
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)
from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.examples.ASU_PR import configuration
from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from ..build_package import build_package

solver = get_solver("ipopt")

def assert_approx(value, expected_value, error_margin):
    percent_error = error_margin / 100
    tolerance = abs(percent_error * expected_value)
    assert pytest.approx(value, abs=tolerance) == expected_value


def _as_quantity(x):
    unit = pyunits.get_units(x)
    if unit is None:
        unit = pyunits.dimensionless
    return value(x) * unit._get_pint_unit()


class TestParamBlock(object):
    def test_build(self):
        model = ConcreteModel()
        model.params = build_package("peng-robinson", ["nitrogen", "argon", "oxygen"], ["Liq","Vap"])

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 3
        for i in model.params.component_list:
            assert i in ["nitrogen", "argon", "oxygen"]
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 6
        for i in model.params._phase_component_set:
            assert i in [
                ("Liq", "nitrogen"),
                ("Liq", "argon"),
                ("Liq", "oxygen"),
                ("Vap", "nitrogen"),
                ("Vap", "argon"),
                ("Vap", "oxygen"),
            ]

        assert model.params.config.state_definition == FTPx

        assertStructuredAlmostEqual(
            model.params.config.state_bounds,
            {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (54.361, 150, 500, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            },
            item_callback=_as_quantity,
        )

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): SmoothVLE
        }

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 3
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2", "PE3"]

        assert model.params.phase_equilibrium_list == {
            "PE1": ["nitrogen", ("Vap", "Liq")],
            "PE2": ["argon", ("Vap", "Liq")],
            "PE3": ["oxygen", ("Vap", "Liq")],
        }

        assert_approx(model.params.pressure_ref.value, 101325, 0.2)
        assert_approx(model.params.temperature_ref.value, 298.15, 0.2)

        assert_approx(model.params.nitrogen.mw.value, 28.0135e-3, 0.2)
        assert_approx(model.params.nitrogen.pressure_crit.value, 34e5, 0.2)
        assert_approx(model.params.nitrogen.temperature_crit.value, 126.2, 0.2)

        assert_approx(model.params.argon.mw.value, 39.948e-3, 0.2)
        assert_approx(model.params.argon.pressure_crit.value, 48.98e5, 0.2)
        assert_approx(model.params.argon.temperature_crit.value, 150.86, 0.2)

        assert_approx(model.params.oxygen.mw.value, 31.999e-3, 0.2)
        assert_approx(model.params.oxygen.pressure_crit.value, 50.43e5, 0.2)
        assert_approx(model.params.oxygen.temperature_crit.value, 154.58, 0.2)

        assert_units_consistent(model)


@pytest.mark.skipif(not cubic_roots_available(), reason="Cubic functions not available")
class TestStateBlock(object):

    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = build_package("peng-robinson", ["nitrogen", "argon", "oxygen"], ["Liq", "Vap"])
        model.props = model.params.build_state_block([1], defined_state=True)
        return model

    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.props[1].flow_mol, Var)
        assert value(model.props[1].flow_mol) == 100
        assert model.props[1].flow_mol.ub == 1000
        assert model.props[1].flow_mol.lb == 0

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 1e5
        assert model.props[1].pressure.ub == 1e6
        assert model.props[1].pressure.lb == 5e4

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 150
        assert model.props[1].temperature.ub == 500
        assert model.props[1].temperature.lb == 54.361

        assert isinstance(model.props[1].mole_frac_comp, Var)
        assert len(model.props[1].mole_frac_comp) == 3
        for i in model.props[1].mole_frac_comp:
            assert value(model.props[1].mole_frac_comp[i]) == 1 / 3

        assert_units_consistent(model)

    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "Total Molar Flowrate",
                "Total Mole Fraction",
                "Temperature",
                "Pressure",
            ]

    def test_initialize(self, model):
        # Fix state
        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(85.00)
        model.props[1].pressure.fix(101325)
        model.props[1].mole_frac_comp["nitrogen"].fix(1 / 3)
        model.props[1].mole_frac_comp["argon"].fix(1 / 3)
        model.props[1].mole_frac_comp["oxygen"].fix(1 / 3)

        assert degrees_of_freedom(model.props[1]) == 0

        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(optarg={"tol": 1e-6})

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    def test_solution(self, model):
        # Check phase equilibrium results
        assert_approx(model.props[1].mole_frac_phase_comp[
            "Liq", "nitrogen"
        ].value, 0.1739, 2)
        assert_approx(model.props[1].mole_frac_phase_comp[
            "Vap", "nitrogen"
        ].value, 0.4221, 2)
        assert_approx(model.props[1].phase_frac["Vap"].value, 0.6422, 2)

    def test_report(self, model):
        model.props[1].report()