# Build and solve a heater block.
from property_packages.build_package import build_package
from pytest import approx

# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.heater import Heater


def test_helmholtz():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("helmholtz", ["water"], ["Liq", "Vap"])
    
    m.fs.heater = Heater(property_package=m.fs.properties)
    m.fs.heater.heat_duty.fix(0)
    m.fs.heater.inlet.flow_mol.fix(1)
    #m.fs.heater.inlet.vapor_frac.fix(0)
    #m.fs.heater.inlet.temperature.fix(298) # room temperature in K
    m.fs.heater.inlet.enth_mol.fix(1878.71)
    m.fs.heater.inlet.pressure.fix(101325)
    assert degrees_of_freedom(m) == 0
    m.fs.heater.initialize()
    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)
    assert value(m.fs.heater.control_volume.properties_out[0].temperature) == approx(298)
    assert value(m.fs.heater.outlet.pressure[0]) == approx(101325)
    assert value(m.fs.heater.outlet.flow_mol[0]) == approx(1)

def test_vapor_fraction_state_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("helmholtz", ["water"], ["Liq", "Vap"])

    m.fs.sb = m.fs.properties.build_state_block([0])
    m.fs.sb[0].flow_mol.fix(100)
    m.fs.sb[0].constrain_component(m.fs.sb[0].temperature,373.15) # 100C
    m.fs.sb[0].constrain_component(m.fs.sb[0].vapor_frac,0.8)
    m.fs.sb[0].pressure.value = 101325 # 1 atm as defualt guess.
    m.fs.sb.initialize()
    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)
    assert value(m.fs.sb[0].temperature) == approx(373.15)
    assert value(m.fs.sb[0].flow_mol) == approx(100)
    assert value(m.fs.sb[0].vapor_frac) == approx(0.8)
    assert value(m.fs.sb[0].pressure) == approx(101273)

def test_ammonia():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("helmholtz", ["ammonia"], ["Liq", "Vap"])
    
    m.fs.heater = Heater(property_package=m.fs.properties)

    m.fs.heater.heat_duty.fix(0)
    m.fs.heater.inlet.flow_mol.fix(1)
    m.fs.heater.inlet.enth_mol.fix(1878.71)
    m.fs.heater.inlet.pressure.fix(101325)

    assert degrees_of_freedom(m) == 0

    m.fs.heater.initialize()

    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)

    assert value(m.fs.heater.control_volume.properties_out[0].temperature) == approx(220.850, rel=1e-3)
    assert value(m.fs.heater.outlet.pressure[0]) == approx(101325)
    assert value(m.fs.heater.outlet.flow_mol[0]) == approx(1)



def test_vapor_fraction_ammonia():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("helmholtz", ["ammonia"], ["Liq", "Vap"])
    m.fs.sb = m.fs.properties.build_state_block([0])
    m.fs.sb[0].flow_mol.fix(100)
    m.fs.sb[0].constrain_component(m.fs.sb[0].temperature,283.15) # 100C
    m.fs.sb[0].constrain_component(m.fs.sb[0].vapor_frac,0.8)
    m.fs.sb[0].pressure.value = 700000 # 1 atm as defualt guess.
    m.fs.sb.initialize()
    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)
    assert value(m.fs.sb[0].temperature) == approx(283.15)
    assert value(m.fs.sb[0].flow_mol) == approx(100)
    assert value(m.fs.sb[0].vapor_frac) == approx(0.8)
    assert value(m.fs.sb[0].pressure) == approx(614294.349)
    # check that temperature_bubble exists
    assert value(m.fs.sb[0].temperature_sat_liq) == approx(283.1266) 
    assert value(m.fs.sb[0].temperature_sat_vap) == approx(283.1266)

