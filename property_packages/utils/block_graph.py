import matplotlib.pyplot as plt
from pyomo.environ import Var, value
import json
from idaes.core.solvers import get_solver
import matplotlib

#matplotlib.use('Agg')

def solver_graph(blk, step=1):
    # solve the block and see how the variables are changing per iteration
    initial_values = {}
    for v in blk.component_objects(Var):
        if v.is_indexed():
            for idx in v:
                if v[idx].is_fixed():
                    # skip fixed variables
                    continue
                initial_values[str(v[idx])] = value(v[idx], exception=False)
        else:
            if v.is_fixed():
                # skip fixed variables
                continue
            initial_values[str(v)] = value(v)
    # add to this dictionary the values of the variables per iteration
    var_values = {}
    max_iter = 0
    keep_solving = True
    opt = get_solver("ipopt")
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
 
def normalize_variables(vv):
    var_values = {}
    for v in vv:
        if(max(vv[v]) < 0): # all values are negative
            var_values[v] = [(-1*x / min(vv[v])) for x in vv[v]]
        else:
            var_values[v] = [(x / max(vv[v])) for x in vv[v]]
        if(max(var_values[v]) > 1):
            print(f"Max in {v} is {max(vv[v])} normalised max is {max(var_values[v])}")
    with open("normalized.json", 'w') as f:
        json.dump(var_values, f, indent=4)
    return var_values

def plot_graph():
    var_values = json.load(open("solver_graph.json", 'r'))
    # var_values = normalize_variables(var_values)
    for v in var_values:
        plt.plot(var_values[v], label=v)
    plt.legend(fontsize="xx-small")
    plt.xlabel("Iterations")
    plt.ylabel("Variable values")
    plt.show()
    # plt.savefig('output_plot.png')


