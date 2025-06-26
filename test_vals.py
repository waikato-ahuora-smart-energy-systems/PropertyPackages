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
m.fs.properties = build_helmholtz_package(["n-butane"])
m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)
state = m.fs.state[0]

state.constrain("flow_mol",1)
m.fs.state.initialize()

    
sat_thermo_data = {
        1: {  # near triple point
            "T": 136,
            "p": 0.00082,
            "rhol": 733.9269,
            "hl": -719.46,
            "sl": -3.02,
            "rhov": 0.000042,
            "hv": -225.72,
            "sv": -0.640,
        },
        2: {  # between critical and triple point
            "T": 222,
            "p": 8.829,
            "rhol": 652.828,
            "hl": -544.84,
            "sl": -2.028,
            "rhov": 0.27998,
            "hv": -117.82,
            "sv": -0.104,
        },
        3: {  # near critical point
            "T": 425,
            "p": 3788.1,
            "rhol": 250.17,
            "hl": 50.78,
            "sl": -0.235,
            "rhov": 205.54,
            "hv": 73.93,
            "sv": -0.180,
        },
    }

MOLAR_WEIGHT = 0.0581222 # kg/mol

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
