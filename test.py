### Imports
from pyomo.environ import ConcreteModel, SolverFactory, SolverStatus, TerminationCondition, Block, TransformationFactory, Constraint, value
from pyomo.network import SequentialDecomposition, Port, Arc
from pyomo.core.base.units_container import _PyomoUnit, units as pyomo_units
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import report_statistics, degrees_of_freedom, fixed_variables_set, activated_block_component_generator
import idaes.logger as idaeslog
from idaes.models.properties.general_helmholtz import HelmholtzParameterBlock, PhaseType, StateVars, AmountBasis
from idaes.models.properties.iapws95 import Iapws95ParameterBlock
from idaes.models.unit_models.heater import Heater
 
 
### Utility Methods
def units(item: str) -> _PyomoUnit:
    ureg = pyomo_units._pint_registry
    pint_unit = getattr(ureg, item)
    return _PyomoUnit(pint_unit, ureg)
 
 
### Build Model
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
 
# Set up property packages
m.fs.properties = HelmholtzParameterBlock(
    pure_component="h2o",
    phase_presentation=PhaseType.MIX,
    state_vars=StateVars.PH,
    amount_basis=AmountBasis.MOLE,
)
# m.fs.properties = Iapws95ParameterBlock()
# Create unit models
 
# HT-772
m.fs.HT_772 = Heater(
    property_package=m.fs.properties,
    has_pressure_change=True,
)
inlet_state = m.fs.HT_772.control_volume.properties_in[0]
outlet_state = m.fs.HT_772.control_volume.properties_out[0]
inlet_state.flow_mol.fix(1.0 * units("mol/s"))
inlet_state.pressure.fix(100.0 * units("kPa"))
# inlet_state.temperature.fix(283.15 * units("K"))
# inlet_state.vapor_frac.fix(0.0)
# inlet_state.enth_mol.fix(80000 * units("J/mol"))
# m.fs.HT_772.deltaP[0].fix(-0.00001 * (inlet_state.enth_mol - outlet_state.enth_mol))
# outlet_state.pressure.fix(99.999 * units("kPa"))
# outlet_state.enth_mol.fix(80000 * units("J/mol"))
# outlet_state.temperature.fix(293.15 * units("K"))
# outlet_state.enth_mol.fix(1513.4 * units("J/mol"))
m.constraint0 = Constraint(expr=inlet_state.temperature == 383.15)  # K
m.contraint1 = Constraint(expr=outlet_state.temperature == 393.15)  # K
m.fs.HT_772.constraint2 = Constraint(expr=m.fs.HT_772.deltaP[0] == (0.0001 * (outlet_state.enth_mol - inlet_state.enth_mol)))
print(hasattr(m.fs.HT_772, "deltaT"))
 
 
### Check Model Status
report_statistics(m)
print("Degrees of freedom:", degrees_of_freedom(m))
 
 
### Initialize Model
m.fs.HT_772.initialize(outlvl=idaeslog.INFO)
 
outlet_state.enth_mol.unfix()
 
m.fs.HT_772.report()
 
fixed = fixed_variables_set(m)
for v in fixed:
    print(v, "fixed to", value(v))
constraints = activated_block_component_generator(m, Constraint)
for c in constraints:
    print(c, "expression:", c.expr)
 
 
### Solve
solver = SolverFactory("ipopt")
result = solver.solve(m, tee=True)
# check for optimal solution
if result.solver.status != SolverStatus.ok or result.solver.termination_condition != TerminationCondition.optimal:
    raise Exception("Solver did not converge to an optimal solution")
 
 
### Report
m.fs.HT_772.report()
 
#print(value(m.fs.HT_772.deltaP[0]))