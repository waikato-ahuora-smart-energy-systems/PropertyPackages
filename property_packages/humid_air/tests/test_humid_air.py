from pyomo.environ import ConcreteModel, value,units,Constraint
import idaes.logger as idaeslog
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Heater
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import MaterialBalanceType
from idaes.core.solvers import get_solver
from CoolProp.CoolProp import HAPropsSI
from property_packages.build_package import build_package

def testInit(): ##not expected to do anything at this stage
    vars = ["flow_mass","flow_vol","vol_mass","enth_mass","entr_mass","total_energy_flow","vapor_frac"]
    indexedVars = ["flow_mol_comp","flow_mass_comp","enth_mol_comp","enth_mass_comp","mass_frac_comp","entr_mol_comp"]
    # haha = HAPropsSI("Hha", "T", 416, "P", 5325647, "psi_w", 0.03)
    # print(haha)
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_package("humid_air", ["water", "air"])
    m.fs.sb = m.fs.properties.build_state_block()
    m.fs.sb.mole_frac_comp["water"].fix(0.00116041013882064)
    m.fs.sb.mole_frac_comp["air"].fix(0.998839589861179)
    m.fs.sb.flow_mol.fix(1)
    m.fs.sb.pressure.fix(108862.4193)
    m.fs.sb.temperature.fix(273.1966)
    m.fs.sb.initialize(outlvl=idaeslog.INFO_HIGH)

    solver = get_solver("ipopt")
    solver.solve(m)

def testConstraints(): ##not expected to do anything at this stage
    vars = ["flow_mass","flow_vol","vol_mass","enth_mass","entr_mass","total_energy_flow","vapor_frac"]
    indexedVars = ["flow_mol_comp","flow_mass_comp","enth_mol_comp","enth_mass_comp","mass_frac_comp","entr_mol_comp"]
    # haha = HAPropsSI("Hha", "T", 416, "P", 5325647, "psi_w", 0.03)
    # print(haha)
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_package("humid_air", ["water", "air"])
    m.fs.sb = m.fs.properties.build_state_block()
    m.fs.sb.mole_frac_comp["water"].fix(0.00116041013882064)
    m.fs.sb.mole_frac_comp["air"].fix(0.998839589861179)
    # m.fs.sb.flow_mol.fix(1)
    m.fs.sb.constrain_component(m.fs.sb.flow_mass, 1)
    m.fs.sb.pressure.fix(108862.4193)
    m.fs.sb.temperature.fix(273.1966)
    m.fs.sb.initialize(outlvl=idaeslog.INFO_HIGH)

    solver = get_solver("ipopt")
    solver.solve(m)