# Build and solve a state block.
from pprint import pprint
from ..build_package import build_package
from pytest import approx


# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver

from pyomo.util.check_units import assert_units_consistent

from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from pyomo.environ import  check_optimal_termination, ConcreteModel, Objective
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
import idaes.core.util.scaling as iscale


from numpy import logspace

def test_gen():

  solver = get_solver(solver="ipopt")

  assert cubic_roots_available() == True

  m = ConcreteModel()
  m.fs = FlowsheetBlock(dynamic=False)
  m.fs.props = build_package("peng-robinson", ["benzene", "toluene"])
  m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

  iscale.calculate_scaling_factors(m.fs.props)
  iscale.calculate_scaling_factors(m.fs.state[1])

  assert_units_consistent(m)

  m.fs.obj = Objective(expr=(m.fs.state[1].temperature - 510) ** 2)
  m.fs.state[1].temperature.setub(600)

  for logP in (10.5, 11.5, 12.0):
      m.fs.obj.deactivate()

      m.fs.state[1].flow_mol.fix(100)
      m.fs.state[1].mole_frac_comp["Benzene"].fix(0.5)
      m.fs.state[1].mole_frac_comp["Toluene"].fix(0.5)
      m.fs.state[1].temperature.fix(300)
      m.fs.state[1].pressure.fix(10 ** (0.5 * logP))

      m.fs.state.initialize()
      
      m.fs.state[1].temperature.unfix()
      m.fs.obj.activate()

      results = solver.solve(m)

      assert check_optimal_termination(results)
      assert m.fs.state[1].flow_mol_phase["Liq"].value <= 1e-5