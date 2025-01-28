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

def debug_block(blks):
    """
    This method allows for solving of Indexed Block components as if they were
    a single Block. A temporary Block object is created which is populated with
    the contents of the objects in the blocks argument and then solved.

    Args:
        solver : a Pyomo solver object to use when solving the Indexed Block
        blocks : an object which inherits from Block, or a list of Blocks
        kwds : a dict of arguments to be passed to the solver

    Returns:
        A Pyomo solver results object
    """
    # We need to play with Pyomo internals for this
    # Check blocks argument, and convert to a list of Blocks
    if isinstance(blks, Block):
        blks = [blks]

    # Create a temporary Block
    tmp = Block(concrete=True)

    nBlocks = len(blks)

    # Iterate over indexed objects
    for i, blk in enumerate(blks):
        # Check that object is a Block
        if not isinstance(blk, Block):
            raise TypeError(
                "Trying to apply solve_indexed_blocks to "
                "object containing non-Block objects"
            )
        # Append components of BlockData to temporary Block
        try:
            tmp._decl["block_%s" % i] = i  # pylint: disable=protected-access
            tmp._decl_order.append(  # pylint: disable=protected-access
                (blk, i + 1 if i < nBlocks - 1 else None)
            )
        except Exception:
            # PYLINT-TODO
            # pylint: disable-next=broad-exception-raised
            raise Exception(
                "solve_indexed_blocks method failed adding "
                "components to temporary block."
            )

        # Set ctypes on temporary Block
        tmp._ctypes[Block] = [  # pylint: disable=protected-access
            0,
            nBlocks - 1,
            nBlocks,
        ]

        initial_values = {}

        for v in tmp.component_objects(Var):
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
            for v in tmp.component_objects(Var):
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
            res = opt.solve(tmp, tee=False)
            # add the values of the variables to the dictionary
            for v in tmp.component_objects(Var):
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
        with open(f"{datetime.now()}.json", 'w') as f:
            json.dump(var_values, f, indent=4)
