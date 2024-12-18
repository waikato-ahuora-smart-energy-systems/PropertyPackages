from ..build_package import build_package
from pytest import approx

from compounds.CompoundDB import get_compound_names, get_compound

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, value

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver

from pyomo.util.check_units import assert_units_consistent

from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from pyomo.environ import  check_optimal_termination, ConcreteModel, Objective
import idaes.core.util.scaling as iscale
from pyomo.environ import assert_optimal_termination

from pprint import pprint

solver = get_solver(solver="ipopt")

from numpy import logspace
import idaes.logger as idaeslog
SOUT = idaeslog.INFO

def test():

    #pprint(get_compound("benzene"))

    #assert False
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])
    m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

    m.fs.state[1].flow_mol.fix(1)
    m.fs.state[1].mole_frac_comp["benzene"].fix(0.5)
    m.fs.state[1].mole_frac_comp["toluene"].fix(0.5)
    m.fs.state[1].temperature.fix(300)
    m.fs.state[1].pressure.fix(101325)

    assert degrees_of_freedom(m) == 0

    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)
