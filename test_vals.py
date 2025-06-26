from property_packages.helmholtz.helmholtz_builder import build_helmholtz_package
from pytest import approx
from pyomo.environ import ConcreteModel, SolverFactory, value, units, assert_optimal_termination
from pyomo.core.base.constraint import Constraint, ScalarConstraint
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Compressor
from idaes.core.util.model_statistics import degrees_of_freedom

def solve(m):
    solver = SolverFactory('ipopt')
    result = solver.solve(m, tee=True)
    assert_optimal_termination(result)


m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.properties = build_helmholtz_package(["i-butane"])
m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)
state = m.fs.state[0]

state.constrain("flow_mol",1)
m.fs.state.initialize()

    
sat_thermo_data = {
        1: {  # not that close to the triple point, but as low as I could get it to go
            "T": 244,
            "p": 48.482,
            "rhol": 732.52996,
            "hl": -699.97,
            "sl": -3.079,
            "rhov": 0.0000096,
            "hv": -225.91,
            "sv": 0.806,
        },
        2: {  # between critical and triple point
            "T": 266,
            "p": 120.89,
            "rhol": 588.69,
            "hl": -417.88,
            "sl": -15.79,
            "rhov": 3.3332,
            "hv": -56.92,
            "sv": -0.222,
        },
        3: {  # near critical point
            "T": 407,
            "p": 3580.1,
            "rhol": 276.83,
            "hl": 7.30,
            "sl": -0.354,
            "rhov": 173.46,
            "hv": 59.64,
            "sv": -0.225,
        },
        
    }


import json

for data in sat_thermo_data.values():
    state.constrain("temperature",data["T"])
    state.pressure.fix(data["p"]*1000)

    solve(m)

    data["rhol"] = value(state.dens_mass_phase["Liq"])
    data["rhov"] = value(state.dens_mass_phase["Vap"])
    data["hl"] = value(state.enth_mass_phase["Liq"]/1000)
    data["hv"] = value(state.enth_mass_phase["Vap"]/1000)
    data["sl"] = value(state.entr_mass_phase["Liq"]/1000)
    data["sv"] = value(state.entr_mass_phase["Vap"]/1000)

    


    state.constraints.del_component(state.constraints.temperature)

print(json.dumps(sat_thermo_data))
