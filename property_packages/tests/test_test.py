# Build and solve a heater block.
from ..build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.heater import Heater
from idaes.core.util.tables import _get_state_from_port

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, Constraint

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the mixer unit model
from idaes.models.unit_models import Mixer, MomentumMixingType

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the BTX_ideal property package to create a properties block for the flowsheet
from idaes.models.properties.activity_coeff_models import BTX_activity_coeff_VLE

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import _get_state_from_port

def assert_approx(value, expected_value, error_margin):
    percent_error = error_margin / 100
    tolerance = abs(percent_error * expected_value)
    assert approx(value, abs=tolerance) == expected_value

def test():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = build_package("peng-robinson", ["carbon dioxide", "benzene"], ["Liq", "Vap"])
    m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

    m.fs.state[1].flow_mol.fix(1)
    m.fs.state[1].mole_frac_comp["carbon dioxide"].fix(0.5)
    m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
    m.fs.state[1].pressure.fix(101325)
    m.fs.state[1].temperature.fix(300)

    m.fs.state.initialize()

    #m.fs.state[1].constraint = Constraint(expr=m.fs.state[1].enth_mol == 54807.613150693374)

    assert degrees_of_freedom(m) == 0

    solver = SolverFactory('ipopt')
    solver.options['tol'] = 1e-3
    solver.solve(m, tee=True)

    
    