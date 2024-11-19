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
from .config_pure import config
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

solver = get_solver(solver="ipopt")
solver.options['tol'] = 1e-6  # Adjust tolerance
solver.options['max_iter'] = 1000  # Allow more iteration

import idaes.logger as idaeslog
SOUT = idaeslog.INFO

def assert_approx(value, expected_value, error_margin):
    percent_error = error_margin / 100
    tolerance = abs(percent_error * expected_value)
    assert approx(value, abs=tolerance) == expected_value

def get_m():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = GenericParameterBlock(**config)
    m.fs.state = m.fs.props.build_state_block([1], defined_state=True)
    iscale.calculate_scaling_factors(m.fs.props)
    iscale.calculate_scaling_factors(m.fs.state[1])
    return m

def test_T_sweep():
    m = get_m()

    assert_units_consistent(m)

    m.fs.state[1].flow_mol.fix(100)
    m.fs.state[1].temperature.fix(300)
    m.fs.state[1].pressure.fix(100000)
    m.fs.state.initialize()


    # results = solver.solve(m)
    
    # assert check_optimal_termination(results)
    # assert m.fs.state[1].flow_mol_phase["Liq"].value <= 1e-2