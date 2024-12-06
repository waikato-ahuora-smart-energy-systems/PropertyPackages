# Build and solve a state block.
from ..build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, value

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver

from pyomo.util.check_units import assert_units_consistent

from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from pyomo.environ import  check_optimal_termination, ConcreteModel, Objective
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom
solver = get_solver(solver="ipopt")

import idaes.logger as idaeslog
SOUT = idaeslog.INFO

"""
Test Suite Sourced From IDAES

URL: https://github.com/IDAES/idaes-pse/blob/41bb3c9728ea227fa8bb0fa1bf35b8deec467783/idaes/models/
properties/modular_properties/examples/tests/test_BT_PR.py
"""

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
    assert approx(value, abs=tolerance) == expected_value

def get_m():
  m = ConcreteModel()
  m.fs = FlowsheetBlock(dynamic=False)
  m.fs.props = build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])
  m.fs.state = m.fs.props.build_state_block([1], defined_state=True)
  iscale.calculate_scaling_factors(m.fs.props)
  iscale.calculate_scaling_factors(sb)
  return m

def test_T_sweep():
  m = get_m()
  sb = m.fs.state[1]

  assert_units_consistent(m)

  m.fs.obj = Objective(expr=(sb.temperature - 510) ** 2)
  sb.temperature.setub(600)

  # Tests a variety of pressures, and makes sure that the benzene-toluene mixture
  # can move from liquid to vapor phase at each pressure
  for logP in [9.5, 10, 10.5, 11, 11.5, 12]:
      m.fs.obj.deactivate()

      sb.flow_mol.fix(100)
      sb.mole_frac_comp["benzene"].fix(0.5)
      sb.mole_frac_comp["toluene"].fix(0.5)
      sb.pressure.fix(10 ** (0.5 * logP))
      m.fs.state.initialize()
      sb.temperature.fix(300)

      results = solver.solve(m)
      assert sb.flow_mol_phase["Vap"].value <= 1e-2

      sb.temperature.unfix()
      m.fs.obj.activate()

      results = solver.solve(m)

      assert check_optimal_termination(results)
      assert sb.flow_mol_phase["Liq"].value <= 1e-2


def test_P_sweep():
  m = get_m()
  sb = m.fs.state[1]

  # Tests initiialization at a variety of temperatures.
  for T in range(370, 500, 25):
      sb.flow_mol.fix(100)
      sb.mole_frac_comp["benzene"].fix(0.5)
      sb.mole_frac_comp["toluene"].fix(0.5)
      sb.pressure.fix(1e5)
      print(T)
      m.fs.state.initialize()
      # Fix temperature (after initialisation, because otherwise the solver will complain about too few degrees of freedom)
      sb.temperature.fix(T)

      results = solver.solve(m)

      assert check_optimal_termination(results)

      while sb.pressure.value <= 1e6:

          results = solver.solve(m)
          assert check_optimal_termination(results)

          sb.pressure.value = sb.pressure.value + 1e5
      # unfix temperature for next iteration
      sb.temperature.unfix()

def test_T350_P1_x5():
  m = get_m()
  sb = m.fs.state[1]
  
  sb.flow_mol.fix(100)
  sb.mole_frac_comp["benzene"].fix(0.5)
  sb.mole_frac_comp["toluene"].fix(0.5)
  sb.temperature.fix(350)
  sb.pressure.fix(1e5)

  # Trigger build of enthalpy and entropy
  sb.enth_mol_phase
  sb.entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert_approx(value(sb._teq[("Vap", "Liq")]), 365, 0.5)
  assert_approx(value(sb.compress_fact_phase["Liq"]), 0.0035346, 0.5)
  assert_approx(value(sb.compress_fact_phase["Vap"]), 0.966749, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "benzene"]), 0.894676, 2) # Investigate
  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "toluene"]), 0.347566, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "benzene"]), 0.971072, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "toluene"]), 0.959791, 0.5)

  assert_approx(value(sb.mole_frac_phase_comp["Liq", "benzene"]), 0.5, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Liq", "toluene"]), 0.5, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "benzene"]), 0.70584, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "toluene"]), 0.29416, 2)

  assert_approx(value(sb.enth_mol_phase["Liq"]), 38942.8, 5)
  assert_approx(value(sb.enth_mol_phase["Vap"]), 78048.7, 5)
  assert_approx(value(sb.entr_mol_phase["Liq"]), -361.794, 5)
  assert_approx(value(sb.entr_mol_phase["Vap"]), -264.0181, 5)

def test_T350_P5_x5():
  m = get_m()
  sb = m.fs.state[1]
  
  sb.flow_mol.fix(100)
  sb.mole_frac_comp["benzene"].fix(0.5)
  sb.mole_frac_comp["toluene"].fix(0.5)
  sb.temperature.fix(350)
  sb.pressure.fix(5e5)

  # Trigger build of enthalpy and entropy
  sb.enth_mol_phase
  sb.entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert_approx(value(sb._teq[("Vap", "Liq")]), 431.47, 0.5)
  assert_approx(value(sb.compress_fact_phase["Liq"]), 0.01766, 0.5)
  assert_approx(value(sb.compress_fact_phase["Vap"]), 0.80245, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "benzene"]), 0.181229, 1.5) # Investigate
  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "toluene"]), 0.070601, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "benzene"]), 0.856523, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "toluene"]), 0.799237, 0.5)

  assert_approx(value(sb.mole_frac_phase_comp["Liq", "benzene"]), 0.5, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Liq", "toluene"]), 0.5, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "benzene"]), 0.65415, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "toluene"]), 0.34585, 2)

  assert_approx(value(sb.enth_mol_phase["Liq"]), 38966.9, 5)
  assert_approx(value(sb.enth_mol_phase["Vap"]), 75150.7, 5)
  assert_approx(value(sb.entr_mol_phase["Liq"]), -361.8433, 5)
  assert_approx(value(sb.entr_mol_phase["Vap"]), -281.9703, 5)

def test_T450_P1_x5():
  m = get_m()
  sb = m.fs.state[1]
  
  sb.flow_mol.fix(100)
  sb.mole_frac_comp["benzene"].fix(0.5)
  sb.mole_frac_comp["toluene"].fix(0.5)
  sb.temperature.fix(450)
  sb.pressure.fix(1e5)

  # Trigger build of enthalpy and entropy
  sb.enth_mol_phase
  sb.entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert_approx(value(sb._teq[("Vap", "Liq")]), 371.4, 0.5)
  assert_approx(value(sb.compress_fact_phase["Liq"]), 0.0033583, 0.5)
  assert_approx(value(sb.compress_fact_phase["Vap"]), 0.9821368, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "benzene"]), 8.069323, 1) # Investigate
  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "toluene"]), 4.304955, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "benzene"]), 0.985365, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "toluene"]), 0.979457, 0.5)

  assert_approx(value(sb.mole_frac_phase_comp["Liq", "benzene"]), 0.29861, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Liq", "toluene"]), 0.70139, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "benzene"]), 0.5, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "toluene"]), 0.5, 2)

  assert_approx(value(sb.enth_mol_phase["Liq"]), 49441.2, 5)
  assert_approx(value(sb.enth_mol_phase["Vap"]), 84175.1, 5)
  assert_approx(value(sb.entr_mol_phase["Liq"]), -328.766, 5)
  assert_approx(value(sb.entr_mol_phase["Vap"]), -241.622, 5)

def test_T450_P5_x5():
  m = get_m()
  sb = m.fs.state[1]
  
  sb.flow_mol.fix(100)
  sb.mole_frac_comp["benzene"].fix(0.5)
  sb.mole_frac_comp["toluene"].fix(0.5)
  sb.temperature.fix(450)
  sb.pressure.fix(5e5)

  # Trigger build of enthalpy and entropy
  sb.enth_mol_phase
  sb.entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert_approx(value(sb._teq[("Vap", "Liq")]), 436.93, 0.5)
  assert_approx(value(sb.compress_fact_phase["Liq"]), 0.0166181, 0.5)
  assert_approx(value(sb.compress_fact_phase["Vap"]), 0.9053766, 0.5)

  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "benzene"]), 1.63308, 1) # investigate
  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "toluene"]), 0.873213, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "benzene"]), 0.927534, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "toluene"]), 0.898324, 0.5)

  assert_approx(value(sb.mole_frac_phase_comp["Liq", "benzene"]), 0.3488737, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Liq", "toluene"]), 0.6511263, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "benzene"]), 0.5, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "toluene"]), 0.5, 2)

  assert_approx(value(sb.enth_mol_phase["Liq"]), 51095.2, 5)
  assert_approx(value(sb.enth_mol_phase["Vap"]), 83362.3, 5)
  assert_approx(value(sb.entr_mol_phase["Liq"]), -326.299, 5)
  assert_approx(value(sb.entr_mol_phase["Vap"]), -256.198, 5)

def test_T368_P1_x5():
  m = get_m()
  sb = m.fs.state[1]
  
  sb.flow_mol.fix(100)
  sb.mole_frac_comp["benzene"].fix(0.5)
  sb.mole_frac_comp["toluene"].fix(0.5)
  sb.temperature.fix(368)
  sb.pressure.fix(1e5)

  # Trigger build of enthalpy and entropy
  sb.enth_mol_phase
  sb.entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert_approx(value(sb._teq[("Vap", "Liq")]), 368, 0.5)
  assert_approx(value(sb.compress_fact_phase["Liq"]), 0.003504, 0.5)
  assert_approx(value(sb.compress_fact_phase["Vap"]), 0.97, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "benzene"]), 1.492049, 1.5) # Investigate
  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "toluene"]), 0.621563, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "benzene"]), 0.97469, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "toluene"]), 0.964642, 0.5)

  assert_approx(value(sb.mole_frac_phase_comp["Liq", "benzene"]), 0.4012128, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Liq", "toluene"]), 0.5987872, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "benzene"]), 0.6141738, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "toluene"]), 0.3858262, 2)

  assert_approx(value(sb.enth_mol_phase["Liq"]), 38235.1, 5)
  assert_approx(value(sb.enth_mol_phase["Vap"]), 77155.4, 5)
  assert_approx(value(sb.entr_mol_phase["Liq"]), -359.256, 5)
  assert_approx(value(sb.entr_mol_phase["Vap"]), -262.348, 5)

def test_T376_P1_x2():
  m = get_m()
  sb = m.fs.state[1]
  
  sb.flow_mol.fix(100)
  sb.mole_frac_comp["benzene"].fix(0.2)
  sb.mole_frac_comp["toluene"].fix(0.8)
  sb.temperature.fix(376)
  sb.pressure.fix(1e5)

  # Trigger build of enthalpy and entropy
  sb.enth_mol_phase
  sb.entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert_approx(value(sb._teq[("Vap", "Liq")]), 376, 0.5)
  assert_approx(value(sb.compress_fact_phase["Liq"]), 0.00361333, 0.5)
  assert_approx(value(sb.compress_fact_phase["Vap"]), 0.968749, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "benzene"]), 1.8394188, 1.5) # Investigate
  assert_approx(value(sb.fug_coeff_phase_comp["Liq", "toluene"]), 0.7871415, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "benzene"]), 0.9763608, 0.5)
  assert_approx(value(sb.fug_coeff_phase_comp["Vap", "toluene"]), 0.9663611, 0.5)

  assert_approx(value(sb.mole_frac_phase_comp["Liq", "benzene"]), 0.17342, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Liq", "toluene"]), 0.82658, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "benzene"]), 0.3267155, 2)
  assert_approx(value(sb.mole_frac_phase_comp["Vap", "toluene"]), 0.6732845, 2)

  assert_approx(value(sb.enth_mol_phase["Liq"]), 31535.8, 5)
  assert_approx(value(sb.enth_mol_phase["Vap"]), 69175.3, 5)
  assert_approx(value(sb.entr_mol_phase["Liq"]), -369.033, 5)
  assert_approx(value(sb.entr_mol_phase["Vap"]), -273.513, 5)

def test_basic_scaling():
  m = get_m()
  sb = m.fs.state[1]
  
  assert len(sb.scaling_factor) == 23
  assert sb.scaling_factor[sb.flow_mol] == 1e-2
  assert sb.scaling_factor[sb.flow_mol_phase["Liq"]] == 1e-2
  assert sb.scaling_factor[sb.flow_mol_phase["Vap"]] == 1e-2
  assert (
      sb.scaling_factor[
          sb.flow_mol_phase_comp["Liq", "benzene"]
      ]
      == 1e-2
  )
  assert (
      sb.scaling_factor[
          sb.flow_mol_phase_comp["Liq", "toluene"]
      ]
      == 1e-2
  )
  assert (
      sb.scaling_factor[
          sb.flow_mol_phase_comp["Vap", "benzene"]
      ]
      == 1e-2
  )
  assert (
      sb.scaling_factor[
          sb.flow_mol_phase_comp["Vap", "toluene"]
      ]
      == 1e-2
  )
  assert (
      sb.scaling_factor[sb.mole_frac_comp["benzene"]]
      == 1000
  )
  assert (
      sb.scaling_factor[sb.mole_frac_comp["toluene"]]
      == 1000
  )
  assert (
      sb.scaling_factor[
          sb.mole_frac_phase_comp["Liq", "benzene"]
      ]
      == 1000
  )
  assert (
      sb.scaling_factor[
          sb.mole_frac_phase_comp["Liq", "toluene"]
      ]
      == 1000
  )
  assert (
      sb.scaling_factor[
          sb.mole_frac_phase_comp["Vap", "benzene"]
      ]
      == 1000
  )
  assert (
      sb.scaling_factor[
          sb.mole_frac_phase_comp["Vap", "toluene"]
      ]
      == 1000
  )
  assert sb.scaling_factor[sb.pressure] == 1e-5
  assert sb.scaling_factor[sb.temperature] == 1e-2
  assert sb.scaling_factor[sb._teq["Vap", "Liq"]] == 1e-2
  assert sb.scaling_factor[sb._t1_Vap_Liq] == 1e-2

  assert (
      sb.scaling_factor[
          sb._mole_frac_tbub["Vap", "Liq", "benzene"]
      ]
      == 1000
  )
  assert (
      sb.scaling_factor[
          sb._mole_frac_tbub["Vap", "Liq", "toluene"]
      ]
      == 1000
  )
  assert (
      sb.scaling_factor[
          sb._mole_frac_tdew["Vap", "Liq", "benzene"]
      ]
      == 1000
  )
  assert (
      sb.scaling_factor[
          sb._mole_frac_tdew["Vap", "Liq", "toluene"]
      ]
      == 1000
  )
  assert (
      sb.scaling_factor[sb.temperature_bubble["Vap", "Liq"]]
      == 1e-2
  )
  assert (
      sb.scaling_factor[sb.temperature_dew["Vap", "Liq"]]
      == 1e-2
  )