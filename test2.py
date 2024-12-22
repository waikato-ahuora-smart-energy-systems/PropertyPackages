
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
from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.state_definitions import FTPx, FPhx
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE

from idaes.models.properties.modular_properties.eos.ceos import cubic_roots_available
from property_packages.build_package import build_package

solver = get_solver("ipopt")

model = ConcreteModel()
model.params = build_package("peng-robinson", ["methane", "hydrogen", "ethane", "propane", 
                                                "n-butane", "isobutane", "ethylene", "propylene", 
                                                "1-butene", "1-pentene", "1-hexene", "1-heptene", 
                                                "1-octene"], ["Liq", "Vap"])

model.props = model.params.build_state_block([1], defined_state=True)

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

model.props.initialize(outlvl=1)

solver.solve(model, tee=True)

assert_optimal_termination(model)