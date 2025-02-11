from pyomo.environ import (
    Block,
    check_optimal_termination,
    Constraint,
    value,
)
from pyomo.network import Arc
from pyomo.dae import ContinuousSet
from pyomo.core.expr.visitor import identify_variables

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import (
    get_activity_dict,
    deactivate_model_at,
    deactivate_constraints_unindexed_by,
    fix_vars_unindexed_by,
    get_derivatives_at,
    copy_values_at_time,
    get_implicit_index_of_set,
)
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver

import matplotlib.pyplot as plt
from pyomo.environ import Var, value
import json
from idaes.core.solvers import get_solver
import matplotlib
import os
from datetime import datetime

step=1

def debug_block(blk):

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
    
    print("Initial values: ", initial_values)
    

    # add to this dictionary the values of the variables per iteration
    var_values = {}
    max_iter = 1
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
        
        if max_iter == 25:
            keep_solving = False
        
        max_iter += step
    if res is not None:
        print(res.solver.termination_condition)

    print("Number of iterations: ", max_iter)

    # write the values to a file
    with open(f"{datetime.now()}.json", 'w') as f:
        json.dump(var_values, f, indent=4)
