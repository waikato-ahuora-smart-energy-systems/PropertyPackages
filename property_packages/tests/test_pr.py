# Build and solve a state block.
from pprint import pprint
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

def get_m():
  m = ConcreteModel()
  m.fs = FlowsheetBlock(dynamic=False)
  m.fs.props = build_package("peng-robinson", ["benzene", "toluene"])
  m.fs.state = m.fs.props.build_state_block([1], defined_state=True)
  iscale.calculate_scaling_factors(m.fs.props)
  iscale.calculate_scaling_factors(m.fs.state[1])
  return m

def test_T_sweep():
  m = get_m()

  assert_units_consistent(m)

  m.fs.obj = Objective(expr=(m.fs.state[1].temperature - 510) ** 2)
  m.fs.state[1].temperature.setub(600)

  for logP in [9.5, 10, 10.5, 11, 11.5, 12]:
      m.fs.obj.deactivate()

      m.fs.state[1].flow_mol.fix(100)
      m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
      m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
      m.fs.state[1].temperature.fix(300)
      m.fs.state[1].pressure.fix(10 ** (0.5 * logP))

      m.fs.state.initialize()

      m.fs.state[1].temperature.unfix()
      m.fs.obj.activate()

      results = solver.solve(m)

      assert check_optimal_termination(results)
      assert m.fs.state[1].flow_mol_phase["Liq"].value <= 1e-2

def test_P_sweep():
  m = get_m()

  for T in range(370, 500, 25):
      m.fs.state[1].flow_mol.fix(100)
      m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
      m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
      m.fs.state[1].temperature.fix(T)
      m.fs.state[1].pressure.fix(1e5)

      m.fs.state.initialize()

      results = solver.solve(m)

      assert check_optimal_termination(results)

      while m.fs.state[1].pressure.value <= 1e6:

          results = solver.solve(m)
          assert check_optimal_termination(results)

          m.fs.state[1].pressure.value = m.fs.state[1].pressure.value + 1e5

def test_T350_P1_x5():
  m = get_m()
  
  m.fs.state[1].flow_mol.fix(100)
  m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
  m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
  m.fs.state[1].temperature.fix(350)
  m.fs.state[1].pressure.fix(1e5)

  # Trigger build of enthalpy and entropy
  m.fs.state[1].enth_mol_phase
  m.fs.state[1].entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), abs=1) == 365
  assert 0.0035346 == approx(
      value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-2
  )
  assert 0.966749 == approx(
      value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-2
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-1
      )
      == 0.894676
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-2
      )
      == 0.347566
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.971072
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-2
      )
      == 0.959791
  )

  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-2
      )
      == 0.5
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-2
      )
      == 0.5
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.70584
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-1
      )
      == 0.29416
  )
  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-2) == 38942.8
  )
  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-2) == 78048.7
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-2) == -361.794
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-2) == -264.0181
  )

def test_T350_P5_x5():
  m = get_m()
  
  m.fs.state[1].flow_mol.fix(100)
  m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
  m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
  m.fs.state[1].temperature.fix(350)
  m.fs.state[1].pressure.fix(5e5)

  # Trigger build of enthalpy and entropy
  m.fs.state[1].enth_mol_phase
  m.fs.state[1].entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-2) == 431.47
  assert (
      approx(value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-2)
      == 0.01766
  )
  assert (
      approx(value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-2)
      == 0.80245
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-1
      )
      == 0.181229
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-2
      )
      == 0.070601
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.856523
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-2
      )
      == 0.799237
  )

  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-2
      )
      == 0.5
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-2
      )
      == 0.5
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.65415
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-2
      )
      == 0.34585
  )

  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-2) == 38966.9
  )
  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-2) == 75150.7
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-2) == -361.8433
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-2) == -281.9703
  )

def test_T450_P1_x5():
  m = get_m()
  
  m.fs.state[1].flow_mol.fix(100)
  m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
  m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
  m.fs.state[1].temperature.fix(450)
  m.fs.state[1].pressure.fix(1e5)

  # Trigger build of enthalpy and entropy
  m.fs.state[1].enth_mol_phase
  m.fs.state[1].entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-2) == 371.4
  assert 0.0033583 == approx(
      value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-2
  )
  assert 0.9821368 == approx(
      value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-2
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-2
      )
      == 8.069323
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-2
      )
      == 4.304955
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.985365
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-2
      )
      == 0.979457
  )

  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-2
      )
      == 0.29861
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-2
      )
      == 0.70139
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.5
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-2
      )
      == 0.5
  )

  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-2) == 49441.2
  )
  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-2) == 84175.1
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-2) == -328.766
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-2) == -241.622
  )

def test_T450_P5_x5():
  m = get_m()
  
  m.fs.state[1].flow_mol.fix(100)
  m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
  m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
  m.fs.state[1].temperature.fix(450)
  m.fs.state[1].pressure.fix(5e5)

  # Trigger build of enthalpy and entropy
  m.fs.state[1].enth_mol_phase
  m.fs.state[1].entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-2) == 436.93
  assert 0.0166181 == approx(
      value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-2
  )
  assert 0.9053766 == approx(
      value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-2
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-2
      )
      == 1.63308
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-2
      )
      == 0.873213
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.927534
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-2
      )
      == 0.898324
  )

  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-2
      )
      == 0.3488737
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-2
      )
      == 0.6511263
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.5
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-2
      )
      == 0.5
  )

  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-2) == 51095.2
  )
  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-2) == 83362.3
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-2) == -326.299
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-2) == -256.198
  )

def test_T368_P1_x5():
  m = get_m()
  
  m.fs.state[1].flow_mol.fix(100)
  m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
  m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
  m.fs.state[1].temperature.fix(368)
  m.fs.state[1].pressure.fix(1e5)

  # Trigger build of enthalpy and entropy
  m.fs.state[1].enth_mol_phase
  m.fs.state[1].entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-2) == 368
  assert 0.003504 == approx(
      value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-2
  )
  assert 0.97 == approx(
      value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-2
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-1
      )
      == 1.492049
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-2
      )
      == 0.621563
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.97469
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-2
      )
      == 0.964642
  )

  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-1
      )
      == 0.4012128
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-1
      )
      == 0.5987872
  )

  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.6141738
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-1
      )
      == 0.3858262
  )

  m.fs.state[1].mole_frac_phase_comp.display()
  m.fs.state[1].enth_mol_phase_comp.display()

  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-2) == 38235.1
  )
  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-2) == 77155.4
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-2) == -359.256
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-2) == -262.348
  )

def test_T376_P1_x2():
  m = get_m()
  
  m.fs.state[1].flow_mol.fix(100)
  m.fs.state[1].mole_frac_comp["benzene"].fix(0.2)
  m.fs.state[1].mole_frac_comp["toluene"].fix(0.8)
  m.fs.state[1].temperature.fix(376)
  m.fs.state[1].pressure.fix(1e5)

  # Trigger build of enthalpy and entropy
  m.fs.state[1].enth_mol_phase
  m.fs.state[1].entr_mol_phase

  m.fs.state.initialize(outlvl=SOUT)

  results = solver.solve(m)

  # Check for optimal solution
  assert check_optimal_termination(results)

  assert approx(value(m.fs.state[1]._teq[("Vap", "Liq")]), 1e-2) == 376
  assert 0.00361333 == approx(
      value(m.fs.state[1].compress_fact_phase["Liq"]), 1e-2
  )
  assert 0.968749 == approx(
      value(m.fs.state[1].compress_fact_phase["Vap"]), 1e-2
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "benzene"]), 1e-1
      )
      == 1.8394188
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Liq", "toluene"]), 1e-2
      )
      == 0.7871415
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.9763608
  )
  assert (
      approx(
          value(m.fs.state[1].fug_coeff_phase_comp["Vap", "toluene"]), 1e-2
      )
      == 0.9663611
  )

  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]), 1e-1
      )
      == 0.17342
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]), 1e-2
      )
      == 0.82658
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]), 1e-2
      )
      == 0.3267155
  )
  assert (
      approx(
          value(m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]), 1e-2
      )
      == 0.6732845
  )

  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Liq"]), 1e-2) == 31535.8
  )
  assert (
      approx(value(m.fs.state[1].enth_mol_phase["Vap"]), 1e-2) == 69175.3
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Liq"]), 1e-2) == -369.033
  )
  assert (
      approx(value(m.fs.state[1].entr_mol_phase["Vap"]), 1e-2) == -273.513
  )

def test_basic_scaling():
  m = get_m()
  
  assert len(m.fs.state[1].scaling_factor) == 23
  assert m.fs.state[1].scaling_factor[m.fs.state[1].flow_mol] == 1e-2
  assert m.fs.state[1].scaling_factor[m.fs.state[1].flow_mol_phase["Liq"]] == 1e-2
  assert m.fs.state[1].scaling_factor[m.fs.state[1].flow_mol_phase["Vap"]] == 1e-2
  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1].flow_mol_phase_comp["Liq", "benzene"]
      ]
      == 1e-2
  )
  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1].flow_mol_phase_comp["Liq", "toluene"]
      ]
      == 1e-2
  )
  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1].flow_mol_phase_comp["Vap", "benzene"]
      ]
      == 1e-2
  )
  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1].flow_mol_phase_comp["Vap", "toluene"]
      ]
      == 1e-2
  )
  assert (
      m.fs.state[1].scaling_factor[m.fs.state[1].mole_frac_comp["benzene"]]
      == 1000
  )
  assert (
      m.fs.state[1].scaling_factor[m.fs.state[1].mole_frac_comp["toluene"]]
      == 1000
  )
  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1].mole_frac_phase_comp["Liq", "benzene"]
      ]
      == 1000
  )
  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1].mole_frac_phase_comp["Liq", "toluene"]
      ]
      == 1000
  )
  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1].mole_frac_phase_comp["Vap", "benzene"]
      ]
      == 1000
  )
  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1].mole_frac_phase_comp["Vap", "toluene"]
      ]
      == 1000
  )
  assert m.fs.state[1].scaling_factor[m.fs.state[1].pressure] == 1e-5
  assert m.fs.state[1].scaling_factor[m.fs.state[1].temperature] == 1e-2
  assert m.fs.state[1].scaling_factor[m.fs.state[1]._teq["Vap", "Liq"]] == 1e-2
  assert m.fs.state[1].scaling_factor[m.fs.state[1]._t1_Vap_Liq] == 1e-2

  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1]._mole_frac_tbub["Vap", "Liq", "benzene"]
      ]
      == 1000
  )
  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1]._mole_frac_tbub["Vap", "Liq", "toluene"]
      ]
      == 1000
  )
  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1]._mole_frac_tdew["Vap", "Liq", "benzene"]
      ]
      == 1000
  )
  assert (
      m.fs.state[1].scaling_factor[
          m.fs.state[1]._mole_frac_tdew["Vap", "Liq", "toluene"]
      ]
      == 1000
  )
  assert (
      m.fs.state[1].scaling_factor[m.fs.state[1].temperature_bubble["Vap", "Liq"]]
      == 1e-2
  )
  assert (
      m.fs.state[1].scaling_factor[m.fs.state[1].temperature_dew["Vap", "Liq"]]
      == 1e-2
  )