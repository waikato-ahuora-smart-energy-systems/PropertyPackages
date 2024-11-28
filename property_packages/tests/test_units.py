from compounds.CompoundDB import get_compound
from property_packages.types import States
from typing import List
from pyomo.environ import units as u
from ..build_package import build_package

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

solver = get_solver(solver="ipopt")

import idaes.logger as idaeslog
SOUT = idaeslog.INFO


def test_units():

    # Build list of compound objects
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.props = build_package("peng-robinson", ["benzene"], ["Liq", "Vap"])
    m.fs.state = m.fs.props.build_state_block([1], defined_state=True)
    iscale.calculate_scaling_factors(m.fs.props)
    iscale.calculate_scaling_factors(m.fs.state[1])

    assert_units_consistent(m)

    m.fs.state[1].enth_mol_phase
    m.fs.state[1].entr_mol_phase

    assert value(m.fs.state[1].enth_mol_phase["Vap"]) == 4

    

