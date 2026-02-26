from ahuora_property_packages.seawater.seawater_builder import build_seawater_package
import pyomo.environ as pe # Pyomo environment
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Heater
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
import matplotlib.pyplot as plt

# Create an empty flowsheet and steam property parameter block.
model = pe.ConcreteModel()
model.fs = FlowsheetBlock(dynamic=False)
model.fs.properties = build_seawater_package(['water', 'tds'])

# Add a Heater model to the flowsheet.
model.fs.heater = Heater(property_package=model.fs.properties)

# Setup the heater model by fixing the inputs and heat duty
model.fs.heater.inlet.temperature.fix(295)
model.fs.heater.inlet.pressure.fix(101325)

sb = model.fs.heater.control_volume.properties_in[0]

# model.fs.heater.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(9)
# model.fs.heater.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "TDS"].fix(0.1)

# model.fs.heater.control_volume.properties_in[0].constrain_component(
# model.fs.heater.control_volume.properties_in[0].mole_frac_comp['H2O']
# , 0.9936)

# model.fs.heater.control_volume.properties_in[0].constrain_component(
# model.fs.heater.control_volume.properties_in[0].mole_frac_comp['TDS']
# , 0.01)



sb.flow_mol_phase_expression = pe.Expression(expr=sb.mole_frac_comp['H2O'])
sb.constrain_component(
sb.flow_mol_phase_expression
, 0.9936)

sb.constrain_component(
sb.flow_mol
, 510)


print("Flow Mass Phase Comp H2O",pe.value(sb.flow_mass_phase_comp["Liq","H2O"]))
print("Flow Mass Phase Comp TDS",pe.value(sb.flow_mass_phase_comp["Liq","TDS"]))



model.fs.heater.heat_duty.fix(20e5)
# model.fs.heater.heat_duty.unfix()

# model.fs.heater.outlet.temperature.fix(450)


dt = DiagnosticsToolbox(model)
dt.display_underconstrained_set()

print("DOF:" , degrees_of_freedom(model))

# Initialize the model.
model.fs.heater.initialize()

print("DOF After:" , degrees_of_freedom(model))

# Touch the properties to make sure they are calculated and added to the graph 
# (since IDAES uses lazy evaluation for properties, they won't be added to the graph until they are accessed)
sb.flow_mol_phase_comp
sb.flow_mol
sb.mole_frac_comp






#results = []
#outlet_temperature = 295
# while outlet_temperature < 1000:

#     model.fs.heater.outlet.temperature.fix(outlet_temperature)

#     solver = pe.SolverFactory("ipopt")
#     solver.solve(model,tee=True)

#     results.append((outlet_temperature, pe.value(model.fs.heater.heat_duty[0])))

#     outlet_temperature += 5

# plt.plot([r[0] for r in results], [r[1] for r in results])
# plt.xlabel("Outlet Temperature (K)")
# plt.ylabel("Heat Duty (W)")
# plt.title("Heat Duty vs Outlet Temperature for Seawater Heater")
# plt.grid()
# plt.show()

solver = pe.SolverFactory("ipopt")
solver.solve(model,tee=True)

# duty = 12e5
# while duty < 10e6:
#     model.fs.heater.heat_duty.fix(duty)
#     solver.solve(model)
#     duty = duty + 5e6
# solver.solve(model,tee=True)


print("Flow Mol Phase Comp H2O",pe.value(sb.flow_mol_phase_comp["Liq","H2O"]))
print("Flow Mol Phase Comp TDS",pe.value(sb.flow_mol_phase_comp["Liq","TDS"]))
print("Flow Mol Total",pe.value(sb.flow_mol))
print("Mole Frac H2O",pe.value(sb.mole_frac_comp["H2O"]))
print("Mole Frac TDS",pe.value(sb.mole_frac_comp["TDS"]))

print("Flow Mass Phase Comp H2O",pe.value(model.fs.heater.control_volume.properties_out[0].flow_mass_phase_comp["Liq","H2O"]))



sb_out = model.fs.heater.control_volume.properties_out[0]
print("T out",pe.value(sb_out.temperature))
print("work",pe.value(model.fs.heater.heat_duty[0]))

