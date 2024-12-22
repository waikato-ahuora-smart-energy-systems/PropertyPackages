from pyomo.environ import Block, Constraint
from pyomo.core.base.expression import ScalarExpression
from pyomo.core.base.var import ScalarVar, _GeneralVarData, VarData
from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import ConfigurationError
from idaes.models.properties.modular_properties.base.generic_property import (
    _GenericStateBlock,
    GenericParameterData,
    GenericStateBlockData,
)
import idaes.models.properties.modular_properties.base.utility as utility
from property_packages.utils.add_extra_expressions import add_extra_expressions

# increase max iterations for estimating values for bubble, dew, and
# critical point initialization (ie. temperature_bubble, temperature_dew)
utility.MAX_ITER = 1000


def fix_state_vars(blk, state_args=None):
    """
    Method for fixing state variables within StateBlocks. Method takes an
    optional argument of values to use when fixing variables.

    This method is analogous to the fix_state_vars method in
    idaes.core.util.initialization, but it allows us to handle
    the case where we have constraints that define the state.

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
        available_constraints = [
            "enth_mol",
            "enth_mass",
            "entr_mol",
            "entr_mass",
            "smooth_temperature",
        ]
        for n, v in b.define_state_vars().items():
            fix_var = True
            if n == "flow_mol":
                for prop_name in ["flow_mass", "flow_vol"]:
                    if hasattr(b.constraints, prop_name):
                        fix_var = False
                        break
            elif n in ["temperature", "pressure"]:
                if not v.is_fixed():
                    # check if any of the constraints exist
                    for prop_name in available_constraints:
                        if hasattr(b.constraints, prop_name):
                            # don't fix this variable - it is defined by a constraint
                            available_constraints.remove(prop_name)
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


class _ExtendedGenericStateBlock(_GenericStateBlock):

    def initialize(blk, *args, **kwargs):
        flag_dict = fix_state_vars(blk, kwargs.get("state_args", None))

        # Set state_vars_fixed to True to avoid fixing state variables
        # during the initialize method, since this would overdefine
        # the block if we are using constraints
        kwargs["state_vars_fixed"] = True

        # Call the base class initialize method
        super().initialize(*args, **kwargs)

        if kwargs.get("hold_state") is True:
            return flag_dict
        else:
            blk.release_state(flag_dict)


@declare_process_block_class(
    "GenericExtendedStateBlock", block_class=_ExtendedGenericStateBlock
)
class GenericExtendedStateBlockData(GenericStateBlockData):

    def build(self, *args):
        super().build(*args)
        # Add expressions for smooth_temperature, enthalpy in terms of mass, etc.
        add_extra_expressions(self)
        # Add a block for constraints, so we can disable or enable them in bulk
        self.constraints = Block()

    def constrain(self, name: str, value: float):
        # Value must be a float. TODO: Handle unit conversion.
        var = getattr(self, name)
        if type(var) == ScalarExpression:
            c = Constraint(expr=var == value)
            c.abcdef = True
            self.constraints.add_component(name, c)
        elif type(var) in (ScalarVar, _GeneralVarData, VarData):
            var.fix(value)
        else:
            raise Exception(
                f"Variable {self} {name} is not a Var or Expression: {type(var)}"
            )


@declare_process_block_class("GenericExtendedParameterBlock")
class GenericExtendedParameterData(GenericParameterData):
    def build(self):
        super().build()
        self._state_block_class = GenericExtendedStateBlock  # noqa: F821
