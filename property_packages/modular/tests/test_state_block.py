"""
This file tests initialising solving single state blocks
"""

from pyomo.environ import ConcreteModel, SolverFactory, units, assert_optimal_termination, Var, value, Constraint
from idaes.core import FlowsheetBlock
from property_packages.build_package import build_package
from idaes.core.util.model_statistics import degrees_of_freedom

import json


def report_vars(b, prefix="sb", filename=None):
    f = open(filename, 'w') if filename is not None else None

    def write_output(line):
        if f is not None:
            f.write(line + "\n")
        else:
            # use print if no file is provided
            print(line)

    for v in b.component_objects(Var):
        if v.is_indexed():
            for idx in v:
                if isinstance(idx, tuple):
                    idx2 = str(idx)
                else:
                    idx2 = f"'{idx}'"
                write_output(prefix + "." + v.local_name + "[" + str(idx2) + "].value = " + str(value(v[idx])))
        else:
            write_output(prefix + "." + v.local_name + ".value = " + str(value(v)))


m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.properties = build_package("peng-robinson", ["argon", "oxygen", "nitrogen"], ["Liq", "Vap"])
m.fs.state = m.fs.properties.build_state_block([0], defined_state=True)
sb = m.fs.state[0]

sb.flow_mol.fix(1)
sb.pressure.fix(800000 * units.Pa)
sb.temperature.fix((273.15+25) * units.K)
sb.mole_frac_comp["argon"].fix(1/3)
sb.mole_frac_comp["oxygen"].fix(1/3)
sb.mole_frac_comp["nitrogen"].fix(1/3)

assert degrees_of_freedom(m) == 0

sb.report_vars = lambda filename="vars.txt": report_vars(sb, filename=filename)

# m.fs.state.initialize()

# sb.flow_mol.value = 1
# sb.mole_frac_comp['argon'].value = 0.3333333333333333
# sb.mole_frac_comp['oxygen'].value = 0.3333333333333333
# sb.mole_frac_comp['nitrogen'].value = 0.3333333333333333
sb.pressure.value = 800000
# sb.temperature.value = 298.15
sb.flow_mol_phase['Liq'].value = 9.99999999995449e-06
sb.flow_mol_phase['Vap'].value = 0.99999
# sb.flow_mol_phase['Liq'].value = 0.00132944991115618
# sb.flow_mol_phase['Vap'].value = 0.9986705500888439

# sb.mole_frac_phase_comp[('Liq', 'argon')].value = 0.3720125371751377
# sb.mole_frac_phase_comp[('Liq', 'oxygen')].value = 0.43801074807815793
# sb.mole_frac_phase_comp[('Liq', 'nitrogen')].value = 0.18997671474670427
# sb.mole_frac_phase_comp[('Vap', 'argon')].value = 0.3332818428151573
# sb.mole_frac_phase_comp[('Vap', 'oxygen')].value = 0.33319398469656786
# sb.mole_frac_phase_comp[('Vap', 'nitrogen')].value = 0.33352417248827465

# sb.phase_frac['Liq'].value = 0.00132944991115618
# sb.phase_frac['Vap'].value = 0.9986705500888439
# sb._teq[('Vap', 'Liq')].value = 110.352631340186
# sb.s_Vap_Liq['Vap'].value = 0.2503328049252674
# sb.s_Vap_Liq['Liq'].value = 188.04770146473925
# sb.gp_Vap_Liq['Vap'].value = 4.251743987917596
# sb.gp_Vap_Liq['Liq'].value = 188.0477014587557
# sb.gn_Vap_Liq['Vap'].value = 0.2503328049252674
# sb.gn_Vap_Liq['Liq'].value = 184.0504938509725
# sb.log_mole_frac_phase_comp[('Liq', 'argon')].value = -0.98882772319317
# sb.log_mole_frac_phase_comp[('Liq', 'oxygen')].value = -0.8255118299154842
# sb.log_mole_frac_phase_comp[('Liq', 'nitrogen')].value = -1.6608537682964142
# sb.log_mole_frac_phase_comp[('Vap', 'argon')].value = -1.098766772154597
# sb.log_mole_frac_phase_comp[('Vap', 'oxygen')].value = -1.0990304219839582
# sb.log_mole_frac_phase_comp[('Vap', 'nitrogen')].value = -1.0980399350288836

sb.flow_mol.value = 1
sb.mole_frac_comp['argon'].value = 0.3333333333333333
sb.mole_frac_comp['oxygen'].value = 0.3333333333333333
sb.mole_frac_comp['nitrogen'].value = 0.3333333333333333
sb.pressure.value = 800000.0
sb.temperature.value = 298.15
sb.flow_mol_phase['Liq'].value = 9.99999999995449e-06
sb.flow_mol_phase['Vap'].value = 0.99999
sb.mole_frac_phase_comp[('Liq', 'argon')].value = 0.004135564045118999
sb.mole_frac_phase_comp[('Liq', 'oxygen')].value = 0.004001382229450751
sb.mole_frac_phase_comp[('Liq', 'nitrogen')].value = 0.0034676699341971654
sb.mole_frac_phase_comp[('Vap', 'argon')].value = 0.3333332919781064
sb.mole_frac_phase_comp[('Vap', 'oxygen')].value = 0.33333329331991113
sb.mole_frac_phase_comp[('Vap', 'nitrogen')].value = 0.3333332986569807
sb.phase_frac['Liq'].value = 9.99999999995449e-06
sb.phase_frac['Vap'].value = 0.99999
sb._teq[('Vap', 'Liq')].value = 298.15
sb.s_Vap_Liq['Vap'].value = 0.0001
sb.s_Vap_Liq['Liq'].value = 0.0001
sb.gp_Vap_Liq['Vap'].value = 3.9780552428143774
sb.gp_Vap_Liq['Liq'].value = 0.0001
sb.gn_Vap_Liq['Vap'].value = 0.0001
sb.gn_Vap_Liq['Liq'].value = 0.0001
sb.log_mole_frac_phase_comp[('Liq', 'argon')].value = -5.488131450865674
sb.log_mole_frac_phase_comp[('Liq', 'oxygen')].value = -5.521115271805621
sb.log_mole_frac_phase_comp[('Liq', 'nitrogen')].value = -5.664271697364069
sb.log_mole_frac_phase_comp[('Vap', 'argon')].value = -1.0986124127205792
sb.log_mole_frac_phase_comp[('Vap', 'oxygen')].value = -1.0986124086951645
sb.log_mole_frac_phase_comp[('Vap', 'nitrogen')].value = -1.0986123926839542

# sb.mole_frac_phase_comp[('Liq', 'argon')].value = 0.3720125371751377
# sb.mole_frac_phase_comp[('Liq', 'oxygen')].value = 0.43801074807815793
# sb.mole_frac_phase_comp[('Liq', 'nitrogen')].value = 0.18997671474670427
# sb.mole_frac_phase_comp[('Vap', 'argon')].value = 0.3332818428151573
# sb.mole_frac_phase_comp[('Vap', 'oxygen')].value = 0.33319398469656786
# sb.mole_frac_phase_comp[('Vap', 'nitrogen')].value = 0.33352417248827465

sum_molar_fracs=sum([sb.mole_frac_phase_comp[('Liq', i)].value for i in sb.component_list])
print(sum_molar_fracs)
# scale the mole fractions to sum to 1
for i in sb.component_list:
    sb.mole_frac_phase_comp[('Liq', i)].value /= sum_molar_fracs
    print(sb.mole_frac_phase_comp[('Liq', i)], sb.mole_frac_phase_comp[('Liq', i)].value)
# sb.phase_frac['Liq'].value = 0.00132944991115618
# sb.phase_frac['Vap'].value = 0.9986705500888439
# sb._teq[('Vap', 'Liq')].value = 110.352631340186
# sb.s_Vap_Liq['Vap'].value = 0.2503328049252674
# sb.s_Vap_Liq['Liq'].value = 188.04770146473925
# sb.gp_Vap_Liq['Vap'].value = 4.251743987917596
# sb.gp_Vap_Liq['Liq'].value = 188.0477014587557
# sb.gn_Vap_Liq['Vap'].value = 0.2503328049252674
# sb.gn_Vap_Liq['Liq'].value = 184.0504938509725
# sb.log_mole_frac_phase_comp[('Liq', 'argon')].value = -0.98882772319317
# sb.log_mole_frac_phase_comp[('Liq', 'oxygen')].value = -0.8255118299154842
# sb.log_mole_frac_phase_comp[('Liq', 'nitrogen')].value = -1.6608537682964142
# sb.log_mole_frac_phase_comp[('Vap', 'argon')].value = -1.098766772154597
# sb.log_mole_frac_phase_comp[('Vap', 'oxygen')].value = -1.0990304219839582
# sb.log_mole_frac_phase_comp[('Vap', 'nitrogen')].value = -1.0980399350288836

### Values that we are trying to reach

# opt = SolverFactory('ipopt')
# res = opt.solve(m, tee=True)

# assert_optimal_termination(res)

# for p in range(100000, 1000000, 50000):
#     sb.pressure.value = p
#     opt.solve(m, tee=False)
#     assert_optimal_termination(res)

# sb.pressure.value = 800000
# opt.solve(m, tee=False)

# report_vars(sb)

 

import matplotlib.pyplot as plt
def solver_graph(blk, step=1):
    # solve the block and see how the variables are changing per iteration
    initial_values = {}
    for v in blk.component_objects(Var):
        if v.is_indexed():
            for idx in v:
                if v[idx].is_fixed():
                    # skip fixed variables
                    continue
                initial_values[str(v[idx])] = value(v[idx])
        else:
            if v.is_fixed():
                # skip fixed variables
                continue
            initial_values[str(v)] = value(v)
    # add to this dictionary the values of the variables per iteration
    var_values = {}
    max_iter = 0
    keep_solving = True
    opt = SolverFactory("ipopt")
    res = None
    while keep_solving:
        # reset the values of the variables
        for v in blk.component_objects(Var):
            if v.is_indexed():
                for idx in v:
                    if v[idx].is_fixed():
                        continue
                    v[idx].value = initial_values[str(v[idx])]
            else:
                if v.is_fixed():
                    continue
                v.value = initial_values[str(v)]
        opt.options["max_iter"] = max_iter
        res = opt.solve(blk, tee=False)
        # add the values of the variables to the dictionary
        for v in blk.component_objects(Var):
            if v.is_indexed():
                for idx in v:
                    # skip fixed variables
                    if v[idx].is_fixed():
                        continue
                    var_values.setdefault(str(v[idx]), []).append(value(v[idx]))
            else:
                # skip fixed variables
                if v.is_fixed():
                    continue
                var_values.setdefault(str(v), []).append(value(v))
        if res.solver.termination_condition == "maxIterations":
            max_iter += step
        elif res.solver.termination_condition == "optimal":
            # we're done
            keep_solving = False
        else:
            # infeasible or some other condition
            keep_solving = False
    if res is not None:
        print(res.solver.termination_condition)
    print("Number of iterations: ", max_iter)
    # write the values to a file
    with open("solver_graph.json", 'w') as f:
        json.dump(var_values, f, indent=4)
    # plot the values of the variables per iteration
    for v in var_values:
        plt.plot(var_values[v], label=v)
    plt.legend(fontsize="xx-small")
    # label with the solver status
    plt.title(blk.name + ": " + res.solver.termination_condition)
    # label with number of iterations
    plt.xlabel("Iterations")
    # label with the value of the variables
    plt.ylabel("Variable values")
    plt.show()

solver_graph(sb)
 