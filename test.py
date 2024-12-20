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

def report_vars(blk, file):
    with open(file, "w") as f:
        for c in blk.component_data_objects(Var, descend_into=True, active=True):
            f.write(f"{c}: {c.value} {"fixed" if c.fixed else "unfixed"}\n")
    with open("file5.txt", "r") as f:
        with open(f"{file}_diff", "w") as f2:
            for c in blk.component_data_objects(Var, descend_into=True):
                line = f.readline()
                key, value, _ = line.split(" ")
                if type(c.value) != float:
                    diff = "None"
                else:
                    diff = abs(float(value) - c.value)
                f2.write(f"{c}: {diff} {"fixed" if c.fixed else "unfixed"}\n")

def add_report_vars_to_blk(blk, file):
    blk.report_vars = lambda x = file: report_vars(blk, x)

### Build Model
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

m.fs.properties = build_package("peng-robinson", ["oxygen", "argon", "nitrogen"])
m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)

sb = m.fs.state[0]
sb.flow_mol.fix(1)
sb.mole_frac_comp["oxygen"].fix(0.33)
sb.mole_frac_comp["argon"].fix(0.33)
sb.mole_frac_comp["nitrogen"].fix(0.33)
# sb.enth_mol.fix(20000)
sb.enth_mol.fix(-9)
# sb.enth_mol.value = 100000
# sb.temperature.fix(298.15)
sb.pressure.fix(110000)

# for i, v in sb.mole_frac_comp.items():
#     sb.log_mole_frac_comp[i].value = log(v)

# for i, v in sb.log_mole_frac_comp.items():
#     print(i, v, v.value)

# sb.flow_mol_phase["Liq"].value = 0.5
# sb.flow_mol_phase["Vap"].value = 0.5
# sb.temperature_dew["Vap", "Liq"].value = 85.64626567281583
# sb.temperature_bubble["Vap", "Liq"].value = 7.251561660768924

# sb.calculate_scaling_factors()
# print("flow_mol_scale", get_scaling_factor(sb.flow_mol))

# sb.flow_mol.value = 1
# sb.mole_frac_comp["oxygen"].value = 0.33
# sb.mole_frac_comp["argon"].value = 0.33
# sb.mole_frac_comp["nitrogen"].value = 0.33
# sb.pressure.value = 100000
# sb.enth_mol.value = -8.6
# sb.flow_mol_phase["Liq"].value = 0.04492645320343111 
# sb.flow_mol_phase["Vap"].value = 0.955073546796569 
# sb.mole_frac_phase_comp["Liq","oxygen"].value = 0.4920469916779807 
# sb.mole_frac_phase_comp["Liq","argon"].value = 0.36636366805258 
# sb.mole_frac_phase_comp["Liq","nitrogen"].value = 0.13158934026943928 
# sb.mole_frac_phase_comp["Vap","oxygen"].value = 0.3223773445377108 
# sb.mole_frac_phase_comp["Vap","argon"].value = 0.3282894609147667 
# sb.mole_frac_phase_comp["Vap","nitrogen"].value = 0.33933319454752253 
# sb.temperature.value = 298.14976092335405 
# sb.phase_frac["Liq"].value = 0.04492645320343111 
# sb.phase_frac["Vap"].value = 0.9550735467965689 
# sb._teq["Vap","Liq"].value = 85.6462656725217 
# sb._t1_Vap_Liq.value = 298.1497610092947 
# sb.temperature_bubble["Vap","Liq"].value = 7.251481101358245 
# sb._mole_frac_tbub["Vap","Liq","oxygen"].value = 0.7837418157206659 
# sb._mole_frac_tbub["Vap","Liq","argon"].value = 0.10315729845781076 
# sb._mole_frac_tbub["Vap","Liq","nitrogen"].value = 0.11310135409682268 
# sb.log_mole_frac_comp["oxygen"].value = -1.1086652122703258 
# sb.log_mole_frac_comp["argon"].value = -1.1086634525257204 
# sb.log_mole_frac_comp["nitrogen"].value = -1.108663556211434 
# sb.log_mole_frac_tbub["Vap","Liq","oxygen"].value = -0.24367515217305877 
# sb.log_mole_frac_tbub["Vap","Liq","argon"].value = -2.2714918514296145 
# sb.log_mole_frac_tbub["Vap","Liq","nitrogen"].value = -2.1794632110268286 
# sb.temperature_dew["Vap","Liq"].value = 86.51987459210113 
# sb._mole_frac_tdew["Vap","Liq","oxygen"].value = 0.5019690696002069 
# sb._mole_frac_tdew["Vap","Liq","argon"].value = 0.368685898040895 
# sb._mole_frac_tdew["Vap","Liq","nitrogen"].value = 0.12934503373478523 
# sb.log_mole_frac_tdew["Vap","Liq","oxygen"].value = -0.6892167691974668 
# sb.log_mole_frac_tdew["Vap","Liq","argon"].value = -0.9978102109224087 
# sb.log_mole_frac_tdew["Vap","Liq","nitrogen"].value = -2.0452716854508997 
# sb.log_mole_frac_phase_comp["Liq","oxygen"].value = -0.7091810555101541 
# sb.log_mole_frac_phase_comp["Liq","argon"].value = -1.0041288103439043 
# sb.log_mole_frac_phase_comp["Liq","nitrogen"].value = -2.028069264363777 
# sb.log_mole_frac_phase_comp["Vap","oxygen"].value = -1.1320325421038222 
# sb.log_mole_frac_phase_comp["Vap","argon"].value = -1.1138595569857048 
# sb.log_mole_frac_phase_comp["Vap","nitrogen"].value = -1.0807727795353723 


# sb.temperature_bubble["Vap","Liq"].value = 7.618121925180626 
# sb._mole_frac_tbub["Vap","Liq","oxygen"].value = 0.4527322206422321 
# sb._mole_frac_tbub["Vap","Liq","argon"].value = 0.2685391710602624 
# sb._mole_frac_tbub["Vap","Liq","nitrogen"].value = 0.27872859468496625 
# sb.log_mole_frac_comp["oxygen"].value = -1.0816657621185595 
# sb.log_mole_frac_comp["argon"].value = -1.1086626017918755 
# sb.log_mole_frac_comp["nitrogen"].value = -1.108662599498811 
# sb.log_mole_frac_tbub["Vap","Liq","oxygen"].value = -0.7924544649616694 
# sb.log_mole_frac_tbub["Vap","Liq","argon"].value = -1.3147586692624962 
# sb.log_mole_frac_tbub["Vap","Liq","nitrogen"].value = -1.2775168739875298 
# sb.temperature_dew["Vap","Liq"].value = 91.78871730120795 
# sb._mole_frac_tdew["Vap","Liq","oxygen"].value = 0.4930555368058721 
# sb._mole_frac_tdew["Vap","Liq","argon"].value = 0.36607675681388535 
# sb._mole_frac_tdew["Vap","Liq","nitrogen"].value = 0.14086770638024265 
# sb.log_mole_frac_tdew["Vap","Liq","oxygen"].value = -0.7071334605622126 
# sb.log_mole_frac_tdew["Vap","Liq","argon"].value = -1.0049122494978062 
# sb.log_mole_frac_tdew["Vap","Liq","nitrogen"].value = -1.9599340816608668 




assert degrees_of_freedom(sb) == 0

# report_statistics(m)

add_report_vars_to_blk(sb, "file5.txt")
# try:
m.fs.state.initialize()
# except:
sb.report_vars()
exit()
# report_vars(sb, "file2.txt")

opt = SolverFactory("ipopt")
opt.options["tol"] = 1e-10
opt.solve(sb, tee=True)

sb.report_vars()
# print(value(sb.enth_mol))

# report_vars(sb, "file3.txt")

m.fs.state.report()