"""
Try to test the property package at the bounds
of the state vars. For example, going to the upper
bound of the temperature or pressure.
"""

from pyomo.environ import ConcreteModel, SolverFactory, check_optimal_termination
from idaes.core import FlowsheetBlock
from idaes.models.unit_models.pressure_changer import PressureChanger
from idaes.models.unit_models.heater import Heater
from property_packages.build_package import build_package


def test_compressor():
    """
    Test a compressor with pressure change from lb to ub
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_package(
        "peng-robinson",
        ["methane", "ethane", "propane", "carbon dioxide", "nitrogen"],
        ["Liq", "Vap"],
    )
    m.fs.compressor = PressureChanger(
        property_package=m.fs.properties,
        thermodynamic_assumption="isentropic",
        compressor=True,
    )
    m.fs.compressor.efficiency_isentropic.fix(0.75)
    m.fs.compressor.inlet.flow_mol.fix(1)
    m.fs.compressor.inlet.pressure.fix(m.fs.compressor.inlet.pressure[0].lb)
    m.fs.compressor.inlet.temperature.fix(300)
    m.fs.compressor.inlet.mole_frac_comp[0, "methane"].fix(0.87)
    m.fs.compressor.inlet.mole_frac_comp[0, "ethane"].fix(0.06)
    m.fs.compressor.inlet.mole_frac_comp[0, "propane"].fix(0.03)
    m.fs.compressor.inlet.mole_frac_comp[0, "carbon dioxide"].fix(0.02)
    m.fs.compressor.inlet.mole_frac_comp[0, "nitrogen"].fix(0.02)
    m.fs.compressor.outlet.pressure.fix(m.fs.compressor.outlet.pressure[0].ub)

    m.fs.compressor.initialize(outlvl=1)
    opt = SolverFactory("ipopt")
    res = opt.solve(m, tee=True)
    assert check_optimal_termination(res)


def test_heater():
    """
    Test a heater with temperature change from lb to ub
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = build_package(
        "peng-robinson",
        ["nitrogen", "carbon dioxide", "water"],
        ["Liq", "Vap"],
    )
    m.fs.heater = Heater(
        property_package=m.fs.properties,
        has_pressure_change=False,
    )
    m.fs.heater.inlet.flow_mol.fix(1)
    m.fs.heater.inlet.pressure.fix(100000)
    # temperature currently doesn't have a good method for determining lb
    # we will do 150 K for now
    m.fs.heater.inlet.temperature.fix(150)
    m.fs.heater.inlet.mole_frac_comp[0, "nitrogen"].fix(0.79)
    m.fs.heater.inlet.mole_frac_comp[0, "carbon dioxide"].fix(0.2)
    m.fs.heater.inlet.mole_frac_comp[0, "water"].fix(0.01)
    m.fs.heater.outlet.temperature.fix(m.fs.heater.outlet.temperature[0].ub)

    m.fs.heater.initialize(outlvl=1)
    opt = SolverFactory("ipopt")
    res = opt.solve(m, tee=True)

    assert check_optimal_termination(res)
