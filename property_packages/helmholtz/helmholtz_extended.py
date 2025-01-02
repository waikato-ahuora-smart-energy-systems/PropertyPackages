from pyomo.environ import SolverFactory

from idaes.models.properties.general_helmholtz.helmholtz_state import HelmholtzStateBlockData, _StateBlock
from idaes.models.properties.general_helmholtz.helmholtz_functions import HelmholtzParameterBlockData
from idaes.core import declare_process_block_class
from idaes.core.util.initialization import solve_indexed_blocks, revert_state_vars
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import InitializationError
import idaes.logger as idaeslog

from property_packages.base.state_block_constraints import StateBlockConstraints
from property_packages.utils.fix_state_vars import fix_state_vars


class _ExtendedStateBlock(_StateBlock):
    """
    This class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs) 

    def initialize(blk, *args, **kwargs):
        outlvl = kwargs.get("outlvl", idaeslog.NOTSET)
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="properties")

        flag_dict = fix_state_vars(blk, kwargs.get("state_args", None))

        dof = degrees_of_freedom(blk)
        if dof != 0:
            raise InitializationError(
                f"{blk.name} Unexpected degrees of freedom during "
                f"initialization at property initialization step: {dof}."
            )
        opt = SolverFactory('ipopt')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            try:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
            except ValueError as e:
                if str(e).startswith("No variables appear"):
                    # https://github.com/Pyomo/pyomo/pull/3445
                    pass
                else:
                    raise e

        if kwargs.get("hold_state") is True:
            return flag_dict
        else:
            blk.release_state(flag_dict)
    
    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        revert_state_vars(blk, flags)


@declare_process_block_class("HelmholtzExtendedStateBlock", block_class=_ExtendedStateBlock)
class HelmholtzExtendedStateBlockData(HelmholtzStateBlockData, StateBlockConstraints):

    def build(self, *args):
        HelmholtzStateBlockData.build(self, *args)
        StateBlockConstraints.build(self, *args)


@declare_process_block_class("HelmholtzExtendedParameterBlock")
class HelmholtzExtendedParameterBlockData(HelmholtzParameterBlockData):
    def build(self):
        super().build()
        self._state_block_class = HelmholtzExtendedStateBlock # noqa: F821
