# Build and solve a state block.
from pprint import pprint
from ..build_package import build_package
from pytest import approx


# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.statejunction import StateJunction
from idaes.models.unit_models.heater import Heater
from idaes.core.solvers import get_solver

from pyomo.util.check_units import assert_units_consistent

from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from pyomo.environ import  check_optimal_termination, ConcreteModel, Objective, value
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

from ..modular.template_builder import build_config
from compounds.CompoundDB import get_compound
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom
from ..config import config


from numpy import logspace

def test_gen():

  solver = get_solver(solver="ipopt")

  assert cubic_roots_available() == True

  m = ConcreteModel()
  m.fs = FlowsheetBlock(dynamic=False)
  m.fs.props = GenericParameterBlock(**config)
  m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

  iscale.calculate_scaling_factors(m.fs.props)
  iscale.calculate_scaling_factors(m.fs.state[1])

  assert_units_consistent(m)

  m.fs.obj = Objective(expr=(m.fs.state[1].temperature - 510) ** 2)
  m.fs.state[1].temperature.setub(600)

  for logP in (10.5, 11.5, 12.0):
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
      assert m.fs.state[1].flow_mol_phase["Liq"].value <= 1e-5