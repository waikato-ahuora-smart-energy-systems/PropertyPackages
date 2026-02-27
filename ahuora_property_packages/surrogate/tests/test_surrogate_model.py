from ..surrogate_builder import build_surrogate_package
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock


m = ConcreteModel()
m.fs = FlowsheetBlock() 

