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

import logging

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.phase_equil import (
    CubicComplementarityVLE,
)
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    LogBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.pure import RPP4

import warnings
from pyomo.environ import *

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from idaes.models.properties.modular_properties.examples.BT_PR import configuration

# Suppress only Pyomo warnings
warnings.filterwarnings('ignore', category=UserWarning, module='pyomo')

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


    assert degrees_of_freedom(m) == 0

    m.fs.compressor.initialize()

    solver = SolverFactory('ipopt')
    solver.options['max_iter'] = 500
    res = solver.solve(m, tee=True)

    assert value(m.fs.compressor.outlet.temperature[0]) == approx(512, rel=1e-2)
    assert_optimal_termination(res)

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
    
    assert degrees_of_freedom(m) == 0

    m.fs.compressor.initialize()

    solver = SolverFactory('ipopt')
    solver.options['max_iter'] = 500
    res = solver.solve(m, tee=True)

    assert int(value(m.fs.compressor.deltaP[0])) == approx(-900000, rel=1e-2)
    assert_optimal_termination(res)