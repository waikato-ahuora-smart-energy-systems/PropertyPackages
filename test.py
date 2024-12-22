### Imports
from pyomo.environ import ConcreteModel, SolverFactory, SolverStatus, TerminationCondition, Block, TransformationFactory, Constraint, value, Var, log
from pyomo.network import SequentialDecomposition, Port, Arc
from pyomo.core.base.units_container import _PyomoUnit, units as pyomo_units
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import report_statistics, degrees_of_freedom, fixed_variables_set, activated_block_component_generator
import idaes.logger as idaeslog
from idaes.models.properties.general_helmholtz import HelmholtzParameterBlock, PhaseType, StateVars, AmountBasis
from idaes.models.properties.iapws95 import Iapws95ParameterBlock
from idaes.models.unit_models.heater import Heater
from property_packages.build_package import build_package
from idaes.core.util.scaling import get_scaling_factor, set_scaling_factor

print('imported everything')

def report_vars(blk, file):
    with open(file, "w") as f:
        for c in blk.component_data_objects(Var, descend_into=True, active=True):
            f.write(f"{c}: {c.value} {"fixed" if c.fixed else "unfixed"}\n")
    # with open("file8.txt", "r") as f:
    #     with open(f"{file}_diff", "w") as f2:
    #         for c in blk.component_data_objects(Var, descend_into=True):
    #             line = f.readline()
    #             key, value, _ = line.split(" ")
    #             if type(c.value) not in [float, int]:
    #                 diff = "None - " + type(c.value).__name__
    #             else:
    #                 diff = abs(float(value) - c.value)
    #             f2.write(f"{c}: {diff} {"fixed" if c.fixed else "unfixed"}\n")

def report_constraints(blk, file):
    with open(file, "w") as f:
        for c in blk.component_data_objects(Constraint, descend_into=True, active=True):
            f.write(f"{c}\n")

def add_report_vars_to_blk(blk):
    blk.report_vars = lambda x: report_vars(blk, x)
def set_bubble(blk, val):
    blk.temperature_bubble["Vap", "Liq"].value = val
def set_dew(blk, val):
    blk.temperature_dew["Vap", "Liq"].value = val
def add_set_bubble_to_blk(blk, val):
    blk.set_bubble = lambda: set_bubble(blk, val)
        # set_bubble(blk, val),
        # set_dew(blk, 92.46536639111942)
    # )

### Build Model
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

m.fs.properties = build_package("peng-robinson", ["oxygen", "nitrogen"])
# m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)
m.fs.heater = Heater(property_package=m.fs.properties)
sb1 = m.fs.heater.control_volume.properties_in[0]
sb2 = m.fs.heater.control_volume.properties_out[0]

# sb1.flow_mol.fix(1)
sb1.constrain("flow_mass", 1)
sb1.mole_frac_comp["oxygen"].fix(1/2)
sb1.mole_frac_comp["nitrogen"].fix(1/2)
sb1.pressure.fix(100000)
sb1.constrain("enth_mol", -10)
sb2.constrain("enth_mol", 2400)
# sb2.temperature.fix(400)

m.fs.heater.initialize()
solver = SolverFactory("ipopt")
solver.options["max_iter"] = 100
solver.solve(m, tee=True)

for c in sb1.component_data_objects(Block):
    print(c)

m.fs.heater.report()
exit()

sb = m.fs.state[0]
# add_report_vars_to_blk(sb)

# # report_vars(sb, "file21.txt")
# # report_constraints(sb, "file22.txt")
# # exit()
# sb.flow_mol.fix(1)
# sb.mole_frac_comp["oxygen"].fix(0.5)
# # sb.mole_frac_comp["argon"].fix(0.33)
# sb.mole_frac_comp["nitrogen"].fix(0.5)
# # sb.enth_mol.fix(30000)
# # sb.enth_mol.fix(-9)
# # sb.enth_mol.value = 100000
# sb.constrain("enth_mol", -10)
# # sb.enth_mol.fix(-10)
# # sb.temperature.fix(410)
# # sb.temperature.fix(404)
# # sb._t1_Vap_Liq.value = 280
# # sb._teq["Vap","Liq"].value = 280
# sb.pressure.fix(100000)

# m.fs.state.initialize()
# exit()
# sb.enth_mol.fix(60000)
# except:
    # sb.report_vars("file17.txt")
solver = SolverFactory("ipopt")
solver.options["max_iter"] = 100
solver.solve(m, tee=True)
print(value(sb.enth_mol))
sb.temperature.unfix()
sb.enth_constraint = Constraint(expr = sb.enth_mol == 4000)
solver.solve(m, tee=True)

sb.report_vars("file23.txt")
exit()

# for i, v in sb.mole_frac_comp.items():
#     sb.log_mole_frac_comp[i].value = log(v)

# for i, v in sb.log_mole_frac_comp.items():
#     print(i, v, v.value)

# sb.flow_mol_phase["Liq"].value = 0.5
# sb.flow_mol_phase["Vap"].value = 0.5
# sb.temperature_dew["Vap", "Liq"].value = 85.64626567281583
# sb.temperature_bubble["Vap", "Liq"].value = 7.251561660768924

sb.calculate_scaling_factors()
# print("flow_mol_scale", get_scaling_factor(sb.flow_mol))

sb.flow_mol.value = 1
sb.mole_frac_comp["oxygen"].value = 0.33
sb.mole_frac_comp["argon"].value = 0.33
sb.mole_frac_comp["nitrogen"].value = 0.33
sb.pressure.value = 200000
sb.enth_mol.value = 2800
sb.flow_mol_phase["Liq"].value = 0
sb.flow_mol_phase["Vap"].value = 1
sb.mole_frac_phase_comp["Liq","oxygen"].value = 0.47079352385149564
sb.mole_frac_phase_comp["Liq","argon"].value = 0.368908063480675
sb.mole_frac_phase_comp["Liq","nitrogen"].value = 0.1502984126678293
sb.mole_frac_phase_comp["Vap","oxygen"].value = 0.3218184161151369
sb.mole_frac_phase_comp["Vap","argon"].value = 0.3277390324749561
sb.mole_frac_phase_comp["Vap","nitrogen"].value = 0.340442551409907
# sb.temperature.value = 404.04401232619784
sb.phase_frac["Liq"].value = 0.05491913386860608
sb.phase_frac["Vap"].value = 0.945080866131394
sb._teq["Vap","Liq"].value = 92.4653654398964
# sb._t1_Vap_Liq.value = 404.04401241707967
sb.temperature_bubble["Vap","Liq"].value = 89.81904298694349
sb._mole_frac_tbub["Vap","Liq","oxygen"].value = 0.17962142058397806
sb._mole_frac_tbub["Vap","Liq","argon"].value = 0.2312975489647376
sb._mole_frac_tbub["Vap","Liq","nitrogen"].value = 0.5890810304512843
sb.log_mole_frac_comp["oxygen"].value = -1.1086626245216111
sb.log_mole_frac_comp["argon"].value = -1.1086626245216111
sb.log_mole_frac_comp["nitrogen"].value = -1.1086626245216111
sb.log_mole_frac_tbub["Vap","Liq","oxygen"].value = -1.7169038619409782
sb.log_mole_frac_tbub["Vap","Liq","argon"].value = -1.4640503065810992
sb.log_mole_frac_tbub["Vap","Liq","nitrogen"].value = -0.5291915318704526
sb.temperature_dew["Vap","Liq"].value = 92.46536544012616
sb._mole_frac_tdew["Vap","Liq","oxygen"].value = 0.4839403861534417
sb._mole_frac_tdew["Vap","Liq","argon"].value = 0.3714379049091942
sb._mole_frac_tdew["Vap","Liq","nitrogen"].value = 0.14462170893736398
sb.log_mole_frac_tdew["Vap","Liq","oxygen"].value = -0.7257935489559016
sb.log_mole_frac_tdew["Vap","Liq","argon"].value = -0.9903735757419181
sb.log_mole_frac_tdew["Vap","Liq","nitrogen"].value = -1.9336338495622085
sb.log_mole_frac_phase_comp["Liq","oxygen"].value = -0.7533356593161652
sb.log_mole_frac_phase_comp["Liq","argon"].value = -0.9972078164436171
sb.log_mole_frac_phase_comp["Liq","nitrogen"].value = -1.8951325433712636
sb.log_mole_frac_phase_comp["Vap","oxygen"].value = -1.1337678176115666
sb.log_mole_frac_phase_comp["Vap","argon"].value = -1.1155376199799802
sb.log_mole_frac_phase_comp["Vap","nitrogen"].value = -1.0775088859533726

assert degrees_of_freedom(sb) == 0
import time

# report_statistics(m)
# add_set_bubble_to_blk(sb, 5.25)
add_report_vars_to_blk(sb)
# try:
m.fs.state.initialize()
# except:
    # sb.report_vars("file17.txt")
solver = SolverFactory("ipopt")
solver.options["max_iter"] = 100
solver.solve(m, tee=True)
sb.report_vars("file18.txt")

print(value(sb.enth_mol))

# from idaes.core.util import DiagnosticsToolbox
# dt = DiagnosticsToolbox(sb)
# dt.report_structural_issues()
# exit()

import numpy as np
import matplotlib.pyplot as plt
import json

# Example setup
mean = 403.1255
std_dev = 0.002
temp_pts = list(np.random.normal(mean, std_dev, 30))
result_pts = {}

# Main logic
for pt in temp_pts:
    sb.temperature.set_value(pt)  # Set the temperature in the solver
    try:
        solver.solve(m, tee=True)  # Attempt to solve
        result_pts[pt] = {
            "status": "success",
            "enth_mol": value(sb.enth_mol)  # Fetch the result
        }
    except:
        result_pts[pt] = {
            "status": "failed",
            "enth_mol": None  # No valid result for failures
        }
    print("done", pt)

# Output the result dictionary for reference
print(json.dumps(result_pts, indent=4))

# Separate points based on status
success_points = [(data["enth_mol"], point) for point, data in result_pts.items() if data["status"] == "success"]
failed_points = [(data["enth_mol"], point) for point, data in result_pts.items() if data["status"] == "failed"]

# Sort success points by x (enthalpy)
success_points.sort()  # Sort by the first element (x)

# Unpack points and temperatures for plotting
success_x, success_y = zip(*success_points) if success_points else ([], [])
failed_x, failed_y = zip(*failed_points) if failed_points else ([], [])

# Plot the success points with a line
plt.plot(success_x, success_y, color="green", label="Success", marker="o")  # Line with markers
plt.scatter(failed_x, failed_y, color="red", label="Failed")  # Keep failed points as dots

# Add labels, legend, and grid
plt.xlabel("Enthalpy (Enth_mol)")
plt.ylabel("Temperature")
plt.title("Point Status with Temperature")
plt.legend()
plt.grid(True)

# Show the plot
plt.show()


# sb.report_vars()
exit()
# report_vars(sb, "file2.txt")

opt = SolverFactory("ipopt")
opt.options["tol"] = 1e-10
opt.solve(sb, tee=True)

sb.report_vars()
# print(value(sb.enth_mol))

# report_vars(sb, "file3.txt")

m.fs.state.report()