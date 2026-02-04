from pyomo.environ import ConcreteModel, value
import idaes.logger as idaeslog
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from ahuora_property_packages.build_package import build_package

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
    m.fs.sb.temperature_dry_bulb.fix(273.1966)
    m.fs.sb.initialize(outlvl=idaeslog.INFO_HIGH)

    solver = get_solver("ipopt")
    solver.solve(m)
    m.fs.sb.display()
    for var in vars:
        print(var +": ", value(m.fs.sb.__getattr__(var)))
    for indexedVar in indexedVars:
        print(indexedVar)
        for index in m.fs.sb.__getattr__(indexedVar):
            print(index+": ")#,value(m.fs.sb.__getattr__(indexedVar)[index]))
    # print(value(m.fs.sb.vol_mass))