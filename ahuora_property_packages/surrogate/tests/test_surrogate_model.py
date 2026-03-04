from ahuora_property_packages.surrogate.surrogate_builder import build_surrogate_package
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import SolverFactory


def test_surrogate_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock() 
    m.fs.pp = build_surrogate_package(["water"])
    m.fs.sb = m.fs.pp.build_state_block(m.fs.time)
    m.fs.sb[0].temperature.fix(300)
    m.fs.sb[0].pressure.fix(101325)
    m.fs.sb[0].flow_mol.fix(1)

    assert degrees_of_freedom(m) == 0

    solver = SolverFactory("ipopt")
    results = solver.solve(m, tee=True)
    assert results.solver.termination_condition == "optimal"

    print("Post-solve surrogate state block values:")
    print("Temperature:", m.fs.sb[0].temperature.value)
    print("Pressure:", m.fs.sb[0].pressure.value)
    print("Flow Mol:", m.fs.sb[0].flow_mol.value)
    print("Enthalpy:", m.fs.sb[0].enth_mol.value)
    print("Entropy:", m.fs.sb[0].entr_mol.value)
    print("Dynamic Viscosity:", m.fs.sb[0].dynamic_viscosity.value)
    print("Kinematic Viscosity:", m.fs.sb[0].kinematic_viscosity.value)
    print("Molar Mass:", m.fs.sb[0].flow_mass.value)
    print("Specific Volume:", m.fs.sb[0].specific_volume.value)


