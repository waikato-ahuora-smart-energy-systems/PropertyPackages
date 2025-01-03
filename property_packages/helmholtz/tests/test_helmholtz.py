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