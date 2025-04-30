from idaes.core.util.exceptions import ConfigurationError

"""
This method is analogous to the fix_state_vars method in
idaes.core.util.initialization, but it allows us to handle
the case where we have constraints that define the state.
"""

def fix_state_vars(blk, state_args=None):
    """
    Method for fixing state variables within StateBlocks. Method takes an
    optional argument of values to use when fixing variables.

    Args:
        blk : An IDAES StateBlock object in which to fix the state variables
        state_args : a dict containing values to use when fixing state
                variables. Keys must match with names used in the
                define_state_vars method, and indices of any variables must
                agree.

    Returns:
        A dict keyed by block index, state variable name (as defined by
        define_state_variables) and variable index indicating the fixed status
        of each variable before the fix_state_vars method was applied.
    """
    if state_args is None:
        state_args = {}

    flags = {}
    for k, b in blk.items():
        flow_exprs = {
            "flow_mol",
            "flow_mass",
            "flow_vol"
        }
        other_exprs = {
            "pressure",
            "enth_mol",
            "enth_mass",
            "entr_mol",
            "entr_mass",
            "temperature",
            "total_energy_flow",
            "custom_vapor_frac",
            "vapor_frac",
        }
        for n, v in b.define_state_vars().items():
            fix_var = True
            if n in flow_exprs:
                for prop_name in flow_exprs:
                    if hasattr(b.constraints, prop_name):
                        fix_var = False
                        break
            elif n in other_exprs:
                if not v.is_fixed():
                    # check if any of the constraints exist
                    for prop_name in other_exprs:
                        if hasattr(b.constraints, prop_name):
                            # don't fix this variable - it is defined by a constraint
                            other_exprs.remove(prop_name)
                            fix_var = False
                            break
            for i in v:
                flags[k, n, i] = v[i].is_fixed()
                if fix_var and not v[i].is_fixed():
                    # If not fixed, fix at either guess provided or current value
                    if n in state_args:
                        # Try to get initial guess from state_args
                        try:
                            if i is None:
                                val = state_args[n]
                            else:
                                val = state_args[n][i]
                        except KeyError:
                            raise ConfigurationError(
                                "Indexes in state_args did not agree with "
                                "those of state variable {}. Please ensure "
                                "that indexes for initial guesses are correct.".format(
                                    n
                                )
                            )
                        v[i].fix(val)
                    else:
                        # No guess, try to use current value
                        if v[i].value is not None:
                            v[i].fix()
                        else:
                            # No initial value - raise Exception before this
                            # gets to a solver.
                            raise ConfigurationError(
                                "State variable {} does not have a value "
                                "assigned. This usually occurs when a Var "
                                "is not assigned an initial value when it is "
                                "created. Please ensure all variables have "
                                "valid values before fixing them.".format(v.name)
                            )
    return flags