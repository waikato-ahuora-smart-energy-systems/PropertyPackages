#################################################################################
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

from math import exp
from pprint import pprint
import pytest
from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    Set,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
import pyomo.common.unittest as unittest

from idaes.core import Component
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)

from pyomo.environ import units as u

from idaes.core.solvers import get_solver
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available

from ..build_package import build_package
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver("ipopt")


"""
---------------------------------------------
Goal for absolute error margins
---------------------------------------------
- temp, pressure, specific volume : 0.5%
- compositional things : 2% (mole frac, flow, etc)
- thermodynamic values: 5% (entr, enth)
---------------------------------------------
"""

"""
Helper function to calculate the absolute error margin of a
value and asserts whether or not value lies within range
- Accepts: percent_error (as whole percent)
- Returns: None
"""

def assert_approx(value, expected_value, error_margin):
    percent_error = error_margin / 100
    tolerance = abs(percent_error * expected_value)
    assert pytest.approx(value, abs=tolerance) == expected_value

def _as_quantity(x):
    unit = pyunits.get_units(x)
    if unit is None:
        unit = pyunits.dimensionless
    return value(x) * unit._get_pint_unit()

def build_model():
  m = ConcreteModel()
  m.fs = FlowsheetBlock(dynamic=False)
  m.fs.props = build_package("peng-robinson", ["methane", "hydrogen", "ethane", "propane", "n-butane", 
                                               "isobutane", "ethylene", "propylene", "1-butene", "1-pentene", 
                                               "1-hexene", "1-heptene", "1-octene"])

  m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

  # iscale.calculate_scaling_factors(m.fs.props)
  # iscale.calculate_scaling_factors(m.fs.state[1])

  m.fs.state[1].flow_mol.fix(1)
  m.fs.state[1].temperature.fix(295)
  m.fs.state[1].pressure.fix(1e5)

  m.fs.state[1].mole_frac_comp["hydrogen"].fix(0.077)
  m.fs.state[1].mole_frac_comp["methane"].fix(0.077)
  m.fs.state[1].mole_frac_comp["ethane"].fix(0.077)
  m.fs.state[1].mole_frac_comp["propane"].fix(0.077)
  m.fs.state[1].mole_frac_comp["n-butane"].fix(0.077)
  m.fs.state[1].mole_frac_comp["isobutane"].fix(0.077)
  m.fs.state[1].mole_frac_comp["ethylene"].fix(0.077)
  m.fs.state[1].mole_frac_comp["propylene"].fix(0.077)
  m.fs.state[1].mole_frac_comp["1-butene"].fix(0.077)
  m.fs.state[1].mole_frac_comp["1-pentene"].fix(0.077)
  m.fs.state[1].mole_frac_comp["1-hexene"].fix(0.077)
  m.fs.state[1].mole_frac_comp["1-heptene"].fix(0.077)
  m.fs.state[1].mole_frac_comp["1-octene"].fix(0.076)

  return m

def initialize_model(model):
    model.fs.state.initialize(optarg={"tol": 1e-6})

def test_params():
  model = ConcreteModel()
  model.params = build_package("peng-robinson", ["methane", "hydrogen", "ethane", "propane", 
                                                  "n-butane", "isobutane", "ethylene", "propylene", 
                                                  "1-butene", "1-pentene", "1-hexene", "1-heptene", 
                                                  "1-octene"])

  assert isinstance(model.params.phase_list, Set)
  assert len(model.params.phase_list) == 2
  for i in model.params.phase_list:
      assert i in ["Liq", "Vap"]
  assert model.params.Liq.is_liquid_phase()
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
  assert len(model.params._phase_component_set) == 24
  for i in model.params._phase_component_set:
      assert i in [
          ("Liq", "ethane"),
          ("Vap", "hydrogen"),
          ("Vap", "methane"),
          ("Vap", "ethane"),
          ("Liq", "propane"),
          ("Liq", "n-butane"),
          ("Liq", "isobutane"),
          ("Vap", "propane"),
          ("Vap", "n-butane"),
          ("Vap", "isobutane"),
          ("Liq", "ethylene"),
          ("Liq", "propylene"),
          ("Liq", "1-butene"),
          ("Vap", "ethylene"),
          ("Vap", "propylene"),
          ("Vap", "1-butene"),
          ("Liq", "1-pentene"),
          ("Liq", "1-hexene"),
          ("Liq", "1-heptene"),
          ("Vap", "1-pentene"),
          ("Vap", "1-hexene"),
          ("Vap", "1-heptene"),
          ("Liq", "1-octene"),
          ("Vap", "1-octene"),
      ]

  assert model.params.config.state_definition == FTPx

  unittest.assertStructuredAlmostEqual(
      model.params.config.state_bounds,
      {
          "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
          "temperature": (10, 300, 500, pyunits.K),
          "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
      },
      item_callback=_as_quantity,
  )

  assert model.params.config.phase_equilibrium_state == {
      ("Vap", "Liq"): SmoothVLE
  }

  assert isinstance(model.params.phase_equilibrium_idx, Set)
  assert len(model.params.phase_equilibrium_idx) == 11
  for i in model.params.phase_equilibrium_idx:
      assert i in [
          "PE1",
          "PE2",
          "PE3",
          "PE4",
          "PE5",
          "PE6",
          "PE7",
          "PE8",
          "PE9",
          "PE10",
          "PE11",
      ]

  assert model.params.phase_equilibrium_list == {
      "PE1": {"ethane": ("Vap", "Liq")},
      "PE2": {"propane": ("Vap", "Liq")},
      "PE3": {"n-butane": ("Vap", "Liq")},
      "PE4": {"isobutane": ("Vap", "Liq")},
      "PE5": {"ethylene": ("Vap", "Liq")},
      "PE6": {"propylene": ("Vap", "Liq")},
      "PE7": {"1-butene": ("Vap", "Liq")},
      "PE8": {"1-pentene": ("Vap", "Liq")},
      "PE9": {"1-hexene": ("Vap", "Liq")},
      "PE10": {"1-heptene": ("Vap", "Liq")},
      "PE11": {"1-octene": ("Vap", "Liq")}
  }

  assert model.params.pressure_ref.value == 101325
  assert model.params.temperature_ref.value == 298.15

  assert_approx(model.params.hydrogen.mw.value, 2.016e-3, 0.5)
  assert_approx(model.params.hydrogen.pressure_crit.value, 13.13e5, 0.5)
  assert_approx(model.params.hydrogen.temperature_crit.value, 33.2, 0.5)

  assert_approx(model.params.methane.mw.value, 16.043e-3, 0.5)
  assert_approx(model.params.methane.pressure_crit.value, 46e5, 0.5)
  assert_approx(model.params.methane.temperature_crit.value, 190.4, 0.5)

  assert_approx(model.params.ethane.mw.value, 30.070e-3, 0.5)
  assert_approx(model.params.ethane.pressure_crit.value, 48.8e5, 0.5)
  assert_approx(model.params.ethane.temperature_crit.value, 305.4, 0.5)

  assert_approx(model.params.propane.mw.value, 44.094e-3, 0.5)
  assert_approx(model.params.propane.pressure_crit.value, 42.5e5, 0.5)
  assert_approx(model.params.propane.temperature_crit.value, 369.8, 0.5)

  assert_approx(model.params.__dict__["n-butane"].mw.value, 58.124e-3, 0.5)
  assert_approx(model.params.__dict__["n-butane"].pressure_crit.value, 38.0e5, 0.5)
  assert_approx(model.params.__dict__["n-butane"].temperature_crit.value, 425.2, 0.5)

  assert_approx(model.params.isobutane.mw.value, 58.124e-3, 0.5)
  assert_approx(model.params.isobutane.pressure_crit.value, 36.5e5, 0.5)
  assert_approx(model.params.isobutane.temperature_crit.value, 408.2, 0.5)

  assert_approx(model.params.ethylene.mw.value, 28.054e-3, 0.5)
  assert_approx(model.params.ethylene.pressure_crit.value, 50.5e5, 0.5)
  assert_approx(model.params.ethylene.temperature_crit.value, 282.4, 0.5)

  assert_approx(model.params.propylene.mw.value, 42.081e-3, 0.5)
  assert_approx(model.params.propylene.pressure_crit.value, 46.2e5, 0.5)
  assert_approx(model.params.propylene.temperature_crit.value, 365.0, 0.5)

  assert_approx(model.params.__dict__["1-butene"].mw.value, 56.104e-3, 0.5)
  assert_approx(model.params.__dict__["1-butene"].pressure_crit.value, 40.2e5, 0.5)
  assert_approx(model.params.__dict__["1-butene"].temperature_crit.value, 419.3, 0.5)

  assert_approx(model.params.__dict__["1-pentene"].mw.value, 70.135e-3, 0.5)
  assert_approx(model.params.__dict__["1-pentene"].pressure_crit.value, 35.6e5, 0.5)
  assert_approx(model.params.__dict__["1-pentene"].temperature_crit.value, 464.7, 0.5)

  assert_approx(model.params.__dict__["1-hexene"].mw.value, 84.162e-3, 0.5)
  assert_approx(model.params.__dict__["1-hexene"].pressure_crit.value, 31.43e5, 0.5)
  assert_approx(model.params.__dict__["1-hexene"].temperature_crit.value, 504.0, 0.5)

  assert_approx(model.params.__dict__["1-heptene"].mw.value, 98.189e-3, 0.5)
  assert_approx(model.params.__dict__["1-heptene"].pressure_crit.value, 29.2e5, 0.5)
  assert_approx(model.params.__dict__["1-heptene"].temperature_crit.value, 537.2, 0.5)

  assert_approx(model.params.__dict__["1-octene"].mw.value, 112.216e-3, 0.5)
  assert_approx(model.params.__dict__["1-octene"].pressure_crit.value, 26.8e5, 0.5)
  assert_approx(model.params.__dict__["1-octene"].temperature_crit.value, 566.6, 0.5)

  assert_units_consistent(model)

def test_build():

  model = build_model()

  # Check state variable values and bounds
  assert isinstance(model.fs.state[1].flow_mol, Var)
  assert value(model.fs.state[1].flow_mol) == 1
  assert model.fs.state[1].flow_mol.ub == 1000
  assert model.fs.state[1].flow_mol.lb == 0

  assert isinstance(model.fs.state[1].pressure, Var)
  assert value(model.fs.state[1].pressure) == 1e5
  assert model.fs.state[1].pressure.ub == 1e6 # changed from 1e7 (different bounds) TODO: check this
  assert model.fs.state[1].pressure.lb == 5e4

  assert isinstance(model.fs.state[1].temperature, Var)
  assert value(model.fs.state[1].temperature) == 295
  assert model.fs.state[1].temperature.ub == 500 # changed from 1500 (different bounds) TODO: check this
  assert model.fs.state[1].temperature.lb == 10

  assert isinstance(model.fs.state[1].mole_frac_comp, Var)
  assert len(model.fs.state[1].mole_frac_comp) == 13
  for i in model.fs.state[1].mole_frac_comp:
      assert value(model.fs.state[1].mole_frac_comp[i]) == pytest.approx(
          0.077, abs=1e-2
      )

def test_unit_consistency():
  model=build_model()
  assert_units_consistent(model)

def test_define_state_vars():
  model=build_model()
  sv = model.fs.state[1].define_state_vars()

  assert len(sv) == 4
  for i in sv:
      assert i in ["flow_mol", "mole_frac_comp", "temperature", "pressure"]

def test_define_display_vars():
  model=build_model()
  sv = model.fs.state[1].define_display_vars()

  assert len(sv) == 4
  for i in sv:
      assert i in [
          "Total Molar Flowrate",
          "Total Mole Fraction",
          "Temperature",
          "Pressure",
      ]

def test_initialize():
  model = build_model()
  
  assert degrees_of_freedom(model.fs.state[1]) == 0

  orig_fixed_vars = fixed_variables_set(model)
  orig_act_consts = activated_constraints_set(model)
  
  initialize_model(model)

  assert degrees_of_freedom(model) == 0

  fin_fixed_vars = fixed_variables_set(model)
  fin_act_consts = activated_constraints_set(model)

  assert len(fin_act_consts) == len(orig_act_consts)
  assert len(fin_fixed_vars) == len(orig_fixed_vars)

  for c in fin_act_consts:
      assert c in orig_act_consts
  
  for v in fin_fixed_vars:
      assert v in orig_fixed_vars

def test_solve_and_solution():
  # Solving the model
  model = build_model()
  results = solver.solve(model)
  assert_optimal_termination(results)
  # Check phase equilibrium results
  assert_approx(model.fs.state[1].mole_frac_phase_comp[
      "Vap", "hydrogen"
  ].value, 0.09996, 2)
  assert_approx(model.fs.state[1].mole_frac_phase_comp[
      "Liq", "propylene"
  ].value, 0.01056, 2)
  assert_approx(model.fs.state[1].mole_frac_phase_comp[
      "Vap", "propylene"
  ].value, 0.09681, 2)
  assert_approx(model.fs.state[1].phase_frac["Vap"].value, 0.77026, 2)

def test_report():
  model = build_model()
  model.fs.state[1].report()