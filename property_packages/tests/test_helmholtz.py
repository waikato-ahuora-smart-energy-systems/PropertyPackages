# Build and solve a state block.
from ..build_package import build_package
from pytest import approx


# Import objects from pyomo package 
from pyomo.environ import ConcreteModel, SolverFactory, value, units

# Import the main FlowsheetBlock from IDAES. The flowsheet block will contain the unit model
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.statejunction import StateJunction
from idaes.models.unit_models.heater import Heater

def test_helmholtz():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 
    m.fs.properties = build_package("helmholtz", ["h2o"])
    
    m.fs.state = Heater(property_package=m.fs.properties)
    m.fs.state.heat_duty.fix(0)
    m.fs.state.inlet.flow_mol.fix(1)
    m.fs.state.inlet.vapor_frac.fix(0)
    m.fs.state.inlet.temperature.fix(298)# room temperature in K
    m.fs.state.inlet.pressure.fix(101325)
    assert degrees_of_freedom(m) == 0
    solver = SolverFactory('ipopt')
    solver.solve(m, tee=True)
    assert value(m.fs.state.outlet.temperature[0]) == approx(298)
    assert value(m.fs.state.outlet.pressure[0]) == approx(101325)
    assert value(m.fs.state.outlet.flow_mol[0]) == approx(1)
    
    