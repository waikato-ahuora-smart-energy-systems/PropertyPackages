# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, value

# Import the solver
from idaes.core.solvers import get_solver

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock

# Import the option to set the type of material balance 
from idaes.core import MaterialBalanceType

# Import the separator unit model
from idaes.models.unit_models import Separator

# Import the option to set the basis for splitting
from idaes.models.unit_models.separator import SplittingType

# Import idaes logger to set output levels
import idaes.logger as idaeslog

# Import the modular property package to create a property block for the flowsheet
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock

# Import the BT_Ideal property package to create a configuration file for the GenericParameterBlock
from idaes.models.properties.modular_properties.examples.CO2_H2O_Ideal_VLE import configuration

# Import the degrees_of_freedom function from the idaes.core.util.model_statistics package
# DOF = Number of Model Variables - Number of Model Constraints
from idaes.core.util.model_statistics import degrees_of_freedom

from pyomo.environ import ConcreteModel, SolverFactory, value, units

# Import the build function to create a property package
from property_packages.build_package import build_package

from pytest import approx

def test_separator():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_package("peng-robinson", ["water", "carbon dioxide"], ["Liq", "Vap"])

    m.fs.sep_1 = Separator(
        property_package=m.fs.properties,
        split_basis=SplittingType.totalFlow,
        outlet_list=["a1", "b1", "c1"],
        ideal_separation=False,
        has_phase_equilibrium=False,
    )

    assert degrees_of_freedom(m) == 7

    m.fs.sep_1.inlet.flow_mol.fix(10) # mol/s
    m.fs.sep_1.inlet.mole_frac_comp[0, "water"].fix(0.9)
    m.fs.sep_1.inlet.mole_frac_comp[0, "carbon dioxide"].fix(0.1)
    m.fs.sep_1.inlet.pressure.fix(101325) # Pa
    m.fs.sep_1.inlet.temperature.fix(353) # K

    m.fs.sep_1.split_fraction[0, "a1"].fix(0.2)
    m.fs.sep_1.split_fraction[0, "b1"].fix(0.5)

    assert degrees_of_freedom(m) == 0

    m.fs.sep_1.initialize()

    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)

    assert value(m.fs.sep_1.a1_state[0].flow_mol) == approx(2)
    assert value(m.fs.sep_1.b1_state[0].flow_mol) == approx(5)
    assert value(m.fs.sep_1.c1_state[0].flow_mol) == approx(3)

    assert value(m.fs.sep_1.a1_state[0].mole_frac_comp["water"]) == approx(0.9)
    assert value(m.fs.sep_1.a1_state[0].mole_frac_comp["carbon dioxide"]) == approx(0.1)

    assert value(m.fs.sep_1.b1_state[0].mole_frac_comp["water"]) == approx(0.9)
    assert value(m.fs.sep_1.b1_state[0].mole_frac_comp["carbon dioxide"]) == approx(0.1)

    assert value(m.fs.sep_1.c1_state[0].mole_frac_comp["water"]) == approx(0.9)
    assert value(m.fs.sep_1.c1_state[0].mole_frac_comp["carbon dioxide"]) == approx(0.1)

    assert value(m.fs.sep_1.a1_state[0].pressure) == approx(101325)
    assert value(m.fs.sep_1.b1_state[0].pressure) == approx(101325)
    assert value(m.fs.sep_1.c1_state[0].pressure) == approx(101325)

    assert value(m.fs.sep_1.a1_state[0].temperature) == approx(353)
    assert value(m.fs.sep_1.b1_state[0].temperature) == approx(353)
    assert value(m.fs.sep_1.c1_state[0].temperature) == approx(353)
