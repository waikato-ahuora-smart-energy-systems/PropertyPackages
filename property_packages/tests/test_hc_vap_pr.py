"""#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

import pytest
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Set,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.unittest import assertStructuredAlmostEqual

from idaes.core import (
    MaterialBalanceType,
    EnergyBalanceType,
    MaterialFlowBasis,
    Component,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)
from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from idaes.models.properties.modular_properties.state_definitions import FTPx

from idaes.models.properties.modular_properties.examples.HC_PR_vap import (
    configuration_vap,
)
from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from property_packages.build_package import build_package


# -----------------------------------------------------------------------------
# Get default solver for testing
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


# Test for configuration dictionaries with parameters from Properties of Gases
# and liquids 4th edition
class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()
        model.params = build_package("peng-robinson", ["methane", "hydrogen", "ethane", "propane", 
                                                  "n-butane", "isobutane", "ethylene", "propylene", 
                                                  "1-butene", "1-pentene", "1-hexene", "1-heptene", 
                                                  "1-octene"], ["Vap"])

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i in ["Vap"]
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 13
        for i in model.params.component_list:
            assert i in [
                "hydrogen",
                "methane",
                "ethane",
                "propane",
                "n-butane",
                "isobutane",
                "ethylene",
                "propylene",
                "1-butene",
                "1-pentene",
                "1-hexene",
                "1-heptene",
                "1-octene",
            ]

            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 13
        for i in model.params._phase_component_set:
            assert i in [
                ("Vap", "hydrogen"),
                ("Vap", "methane"),
                ("Vap", "ethane"),
                ("Vap", "propane"),
                ("Vap", "n-butane"),
                ("Vap", "isobutane"),
                ("Vap", "ethylene"),
                ("Vap", "propylene"),
                ("Vap", "1-butene"),
                ("Vap", "1-pentene"),
                ("Vap", "1-hexene"),
                ("Vap", "1-heptene"),
                ("Vap", "1-octene"),
            ]

        assert model.params.config.state_definition == FTPx

        assertStructuredAlmostEqual(
            model.params.config.state_bounds,
            {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (273.15, 300, 1500, pyunits.K),
                "pressure": (5e4, 1e5, 1e7, pyunits.Pa),
            },
            item_callback=_as_quantity,
        )

        assert model.params.pressure_ref.value == 101325
        assert model.params.temperature_ref.value == 298.15

        assert model.params.hydrogen.mw.value == 2.016e-3
        assert model.params.hydrogen.pressure_crit.value == 12.9e5
        assert model.params.hydrogen.temperature_crit.value == 33.2

        assert model.params.methane.mw.value == 16.043e-3
        assert model.params.methane.pressure_crit.value == 46e5
        assert model.params.methane.temperature_crit.value == 190.4

        assert model.params.ethane.mw.value == 30.070e-3
        assert model.params.ethane.pressure_crit.value == 48.8e5
        assert model.params.ethane.temperature_crit.value == 305.4

        assert model.params.propane.mw.value == 44.094e-3
        assert model.params.propane.pressure_crit.value == 42.5e5
        assert model.params.propane.temperature_crit.value == 369.8

        assert model.params.__dict__["n-butane"].mw.value == 58.124e-3
        assert model.params.__dict__["n-butane"].pressure_crit.value == 38.0e5
        assert model.params.__dict__["n-butane"].temperature_crit.value == 425.2

        assert model.params.isobutane.mw.value == 58.124e-3
        assert model.params.isobutane.pressure_crit.value == 36.5e5
        assert model.params.isobutane.temperature_crit.value == 408.2

        assert model.params.ethylene.mw.value == 28.054e-3
        assert model.params.ethylene.pressure_crit.value == 50.5e5
        assert model.params.ethylene.temperature_crit.value == 282.4

        assert model.params.propylene.mw.value == 42.081e-3
        assert model.params.propylene.pressure_crit.value == 46.2e5
        assert model.params.propylene.temperature_crit.value == 365.0

        assert model.params.__dict__["n-butene"].mw.value == 56.104e-3
        assert model.params.__dict__["n-butene"].pressure_crit.value == 40.2e5
        assert model.params.__dict__["n-butene"].temperature_crit.value == 419.3

        assert model.params.pentene.mw.value == 70.135e-3
        assert model.params.pentene.pressure_crit.value == 40.5e5
        assert model.params.pentene.temperature_crit.value == 464.7

        assert model.params.hexene.mw.value == 84.162e-3
        assert model.params.hexene.pressure_crit.value == 31.7e5
        assert model.params.hexene.temperature_crit.value == 504.0

        assert model.params.heptene.mw.value == 98.189e-3
        assert model.params.heptene.pressure_crit.value == 25.4e5
        assert model.params.heptene.temperature_crit.value == 537.2

        assert model.params.octene.mw.value == 112.216e-3
        assert model.params.octene.pressure_crit.value == 26.2e5
        assert model.params.octene.temperature_crit.value == 566.6
        assert_units_consistent(model)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**configuration_vap)

        model.props = model.params.build_state_block([1], defined_state=True)

        # Fix state
        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(295.00)
        model.props[1].pressure.fix(1e5)
        model.props[1].mole_frac_comp["hydrogen"].fix(0.077)
        model.props[1].mole_frac_comp["methane"].fix(0.077)
        model.props[1].mole_frac_comp["ethane"].fix(0.077)
        model.props[1].mole_frac_comp["propane"].fix(0.077)
        model.props[1].mole_frac_comp["n-butane"].fix(0.077)
        model.props[1].mole_frac_comp["isobutane"].fix(0.077)
        model.props[1].mole_frac_comp["ethylene"].fix(0.077)
        model.props[1].mole_frac_comp["propylene"].fix(0.077)
        model.props[1].mole_frac_comp["1-butene"].fix(0.077)
        model.props[1].mole_frac_comp["1-pentene"].fix(0.077)
        model.props[1].mole_frac_comp["1-hexene"].fix(0.077)
        model.props[1].mole_frac_comp["1-heptene"].fix(0.077)
        model.props[1].mole_frac_comp["1-octene"].fix(0.076)

        assert degrees_of_freedom(model.props[1]) == 0

        return model

    @pytest.mark.unit
    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.props[1].flow_mol, Var)
        assert value(model.props[1].flow_mol) == 1
        assert model.props[1].flow_mol.ub == 1000
        assert model.props[1].flow_mol.lb == 0

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 1e5
        assert model.props[1].pressure.ub == 1e7
        assert model.props[1].pressure.lb == 5e4

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 295
        assert model.props[1].temperature.ub == 1500
        assert model.props[1].temperature.lb == 273.15

        assert isinstance(model.props[1].mole_frac_comp, Var)
        assert len(model.props[1].mole_frac_comp) == 13
        for i in model.props[1].mole_frac_comp:
            assert value(model.props[1].mole_frac_comp[i]) == pytest.approx(
                0.077, abs=1e-2
            )

        # Check supporting variables
        assert isinstance(model.props[1].flow_mol_phase, Var)
        assert len(model.props[1].flow_mol_phase) == 1

        assert isinstance(model.props[1].mole_frac_phase_comp, Var)
        assert len(model.props[1].mole_frac_phase_comp) == 13

        assert isinstance(model.props[1].phase_frac, Var)
        assert len(model.props[1].phase_frac) == 1

        assert isinstance(model.props[1].total_flow_balance, Constraint)
        assert len(model.props[1].total_flow_balance) == 1

        assert isinstance(model.props[1].component_flow_balances, Constraint)
        assert len(model.props[1].component_flow_balances) == 13

        assert isinstance(model.props[1].phase_fraction_constraint, Constraint)
        assert len(model.props[1].phase_fraction_constraint) == 1

        assert_units_consistent(model)

    @pytest.mark.skipif(
        not cubic_roots_available(), reason="Cubic functions not available"
    )
    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert str(model.props[1].get_material_flow_terms(p, j)) == str(
                    model.props[1].flow_mol_phase_comp[p, j]
                )

    @pytest.mark.skipif(
        not cubic_roots_available(), reason="Cubic functions not available"
    )
    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert str(model.props[1].get_enthalpy_flow_terms(p)) == str(
                model.props[1]._enthalpy_flow_term[p]
            )

    @pytest.mark.skipif(
        not cubic_roots_available(), reason="Cubic functions not available"
    )
    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                if j == "hydrogen":
                    assert str(model.props[1].get_material_density_terms(p, j)) == str(
                        model.props[1]._material_density_term[p, j]
                    )

    @pytest.mark.skipif(
        not cubic_roots_available(), reason="Cubic functions not available"
    )
    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert str(model.props[1].get_energy_density_terms(p)) == str(
                model.props[1]._energy_density_term[p]
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.props[1].default_material_balance_type()
            == MaterialBalanceType.componentTotal
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.props[1].default_energy_balance_type()
            == EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_get_material_flow_basis(self, model):
        assert model.props[1].get_material_flow_basis() == MaterialFlowBasis.molar

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

    @pytest.mark.unit
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

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
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

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        # Check phase equilibrium results
        assert model.props[1].mole_frac_phase_comp[
            "Vap", "hydrogen"
        ].value == pytest.approx(0.077, abs=1e-4)
        assert model.props[1].mole_frac_phase_comp[
            "Vap", "propylene"
        ].value == pytest.approx(0.077, abs=1e-4)
        assert model.props[1].phase_frac["Vap"].value == pytest.approx(1.000, abs=1e-4)

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, model):
        model.props[1].report()"""