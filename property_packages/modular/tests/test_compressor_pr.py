# Build and solve a heater block.
from property_packages.build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units, Var

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent
from idaes.models.unit_models.heater import Heater
from idaes.models.unit_models.pressure_changer import Pump
from idaes.models.unit_models.pressure_changer import PressureChanger, ThermodynamicAssumption
from property_packages.utils import solver_graph
import idaes.core.util.scaling as iscale
from property_packages.utils.debug_block import debug_block
from pyomo.environ import (
    Objective,
    check_optimal_termination,
    Constraint,
    Block,
)
import json

import warnings
from pyomo.environ import *

# Suppress only Pyomo warnings
warnings.filterwarnings('ignore', category=UserWarning, module='pyomo')

def assert_approx(value, expected_value, error_margin):
    percent_error = error_margin / 100
    tolerance = abs(percent_error * expected_value)
    assert approx(value, abs=tolerance) == expected_value

def test_compressor_asu():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("peng-robinson", ["argon", "oxygen", "nitrogen"], ["Liq", "Vap"])

    m.fs.compressor = PressureChanger(
        property_package=m.fs.properties, 
        thermodynamic_assumption=ThermodynamicAssumption.isentropic,
        compressor=True
    )

    m.fs.compressor.inlet.flow_mol[0].fix(1 * units.mol/units.s)
    m.fs.compressor.inlet.pressure.fix(200000 * units.Pa)
    m.fs.compressor.inlet.temperature.fix((273.15+25) * units.K)
    m.fs.compressor.inlet.mole_frac_comp[0, "argon"].fix(1/3)
    m.fs.compressor.inlet.mole_frac_comp[0, "oxygen"].fix(1/3)
    m.fs.compressor.inlet.mole_frac_comp[0, "nitrogen"].fix(1/3)
    m.fs.compressor.outlet.pressure.fix(800000*units.Pa)
    m.fs.compressor.efficiency_isentropic.fix(0.75)
    
    m.fs.compressor.control_volume.properties_in[0].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.compressor.control_volume.properties_in[0].eps_z_Vap_Liq.set_value(1e-4)
    m.fs.compressor.control_volume.properties_out[0].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.compressor.control_volume.properties_out[0].eps_z_Vap_Liq.set_value(1e-4)

    iscale.calculate_scaling_factors(m.fs)
    iscale.calculate_scaling_factors(m.fs.compressor.control_volume)

    assert degrees_of_freedom(m) == 0

    m.fs.compressor.initialize()

    assert degrees_of_freedom(m) == 0

    solver = SolverFactory('ipopt')
    solver.options['max_iter'] = 500
    solver.solve(m, tee=True)

    assert_approx(value(m.fs.compressor.outlet.temperature[0]), 512.38, 0.1)

def test_expander_asu():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("peng-robinson", ["argon", "oxygen", "nitrogen"], ["Liq", "Vap"])
    
    m.fs.compressor = PressureChanger(
        property_package=m.fs.properties, 
        thermodynamic_assumption=ThermodynamicAssumption.adiabatic,
        compressor=False)
    
    m.fs.compressor.inlet.flow_mol[0].fix(1)
    m.fs.compressor.inlet.pressure.fix(1000000 * units.Pa)
    m.fs.compressor.inlet.temperature.fix((273.15+25) * units.K)
    m.fs.compressor.inlet.mole_frac_comp[0, "argon"].fix(1/3)
    m.fs.compressor.inlet.mole_frac_comp[0, "oxygen"].fix(1/3)
    m.fs.compressor.inlet.mole_frac_comp[0, "nitrogen"].fix(1/3)

    m.fs.compressor.outlet.pressure.fix(100000*units.Pa)

    m.fs.compressor.control_volume.properties_in[0].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.compressor.control_volume.properties_in[0].eps_z_Vap_Liq.set_value(1e-4)
    m.fs.compressor.control_volume.properties_out[0].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.compressor.control_volume.properties_out[0].eps_z_Vap_Liq.set_value(1e-4)

    iscale.calculate_scaling_factors(m.fs)
    iscale.calculate_scaling_factors(m.fs.compressor.control_volume)

    assert degrees_of_freedom(m) == 0

    m.fs.compressor.initialize()

    assert degrees_of_freedom(m) == 0

    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)

    assert_approx(int(value(m.fs.compressor.deltaP[0])), -900000, 0)

def test_solved_blocks():
    
    res = []

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.props = build_package("peng-robinson", ["argon", "oxygen", "nitrogen"], ["Liq", "Vap"])
    m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

    # Required setup of CCVLE
    m.fs.state[1].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.state[1].eps_z_Vap_Liq.set_value(1e-4)
    # iscale.calculate_scaling_factors(m.fs.props)
    # iscale.calculate_scaling_factors(m.fs.state[1])

    m.fs.state[1].flow_mol.fix(1)
    m.fs.state[1].temperature.fix(300)
    m.fs.state[1].mole_frac_comp["argon"].fix(1/3)
    m.fs.state[1].mole_frac_comp["oxygen"].fix(1/3)
    m.fs.state[1].mole_frac_comp["nitrogen"].fix(1/3)
    
    for p in range (50000, 800000, 50000): # fails on 750000

        # Fixes initialisation
        m.fs.state[1].pressure.fix(p)
        assert degrees_of_freedom(m) == 0
        m.fs.state.initialize()
        res.append(p)
        m.fs.state[1].pressure.unfix()
        
    json.dump(res, open("solved_blocks.json", 'w'))

# Initialising multiple chages the state of a block in a way that breaks the test?
# 


# Potential points of failure

# 1. Cannot be bad Teq guess (same).
# 2. Cannot be wrong T (fixed).
# 3. Cannot be s_vap (0)

# Issues 
# s_liq is 0 rather than ~T
# this is because F-liq is not zero?
# why is F-liq close but not zero?

# Investigate iscale???

def test_2():

    # for p in range (50000, 1000000, 50000):

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.props = build_package("peng-robinson", ["argon", "nitrogen", "oxygen"], ["Liq", "Vap"])
    m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

    # Required setup of CCVLE
    m.fs.state[1].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.state[1].eps_z_Vap_Liq.set_value(1e-4)
    
    # Fixes initialisation
    m.fs.state[1].flow_mol.fix(1)
    m.fs.state[1].pressure.fix(200000)
    m.fs.state[1].temperature.fix(300)
    m.fs.state[1].mole_frac_comp["argon"].fix(1/3)
    m.fs.state[1].mole_frac_comp["nitrogen"].fix(1/3)
    m.fs.state[1].mole_frac_comp["oxygen"].fix(1/3)

    # Look into scaling causes failures & passes
    iscale.calculate_scaling_factors(m.fs.props)
    iscale.calculate_scaling_factors(m.fs.state[1])

    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 0
    m.fs.state.initialize()

def test_3():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.props = build_package("peng-robinson", ["benzene", "toluene"], ["Liq", "Vap"])
    m.fs.state = m.fs.props.build_state_block([1], defined_state=True)

    # Required setup of CCVLE
    m.fs.state[1].eps_t_Vap_Liq.set_value(1e-4)
    m.fs.state[1].eps_z_Vap_Liq.set_value(1e-4)
    iscale.calculate_scaling_factors(m.fs.props)
    iscale.calculate_scaling_factors(m.fs.state[1])
    
    # Fixes initialisation
    m.fs.state[1].flow_mol.fix(100)
    m.fs.state[1].pressure.fix(350000)
    m.fs.state[1].temperature.fix(300)
    m.fs.state[1].mole_frac_comp["benzene"].fix(1/2)
    m.fs.state[1].mole_frac_comp["toluene"].fix(1/2)

    # m.fs.state[1].display()

    assert degrees_of_freedom(m) == 0
    m.fs.state.initialize()

    assert 1 == 0