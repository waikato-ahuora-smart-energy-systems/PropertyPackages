from pyomo.environ import ConcreteModel, SolverFactory, value, units
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.pressure_changer import PressureChanger, ThermodynamicAssumption
from property_packages.build_package import build_package

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False) 
m.fs.properties = build_package("peng-robinson", ["argon", "oxygen", "nitrogen"], ["Liq", "Vap"])

m.fs.compressor = PressureChanger(
    property_package=m.fs.properties, 
    thermodynamic_assumption=ThermodynamicAssumption.isentropic,
    compressor=True
)

m.fs.compressor.inlet.flow_mol[0].fix(2)
m.fs.compressor.inlet.pressure.fix(100000 * units.Pa)
m.fs.compressor.inlet.temperature.fix((273.15+25) * units.K)
m.fs.compressor.inlet.mole_frac_comp[0, "argon"].fix(1/3)
m.fs.compressor.inlet.mole_frac_comp[0, "oxygen"].fix(1/3)
m.fs.compressor.inlet.mole_frac_comp[0, "nitrogen"].fix(1/3)
m.fs.compressor.efficiency_isentropic.fix(0.75)

m.fs.compressor.outlet.pressure.fix(120000*units.Pa)

assert degrees_of_freedom(m) == 0

m.fs.compressor.initialize(outlvl=1)

assert degrees_of_freedom(m) == 0

solver = SolverFactory('ipopt')

solver.solve(m, tee=True)

# assert_approx(int(value(m.fs.compressor.deltaP[0])), -900000, 0)