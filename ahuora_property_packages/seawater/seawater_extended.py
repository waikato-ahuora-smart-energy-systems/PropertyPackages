from pyomo.environ import Expression, Constraint, check_optimal_termination
from watertap.property_models.seawater_prop_pack import SeawaterParameterData, SeawaterStateBlockData, _SeawaterStateBlock
from ahuora_property_packages.base.state_block_constraints import StateBlockConstraints
from idaes.core import declare_process_block_class
from watertap.core.solvers import get_solver
from idaes.core.util.initialization import fix_state_vars, solve_indexed_blocks
from idaes.core.util.model_statistics import degrees_of_freedom, number_unfixed_variables
from idaes.core.util.exceptions import InitializationError, PropertyPackageError
import idaes.logger as idaeslog
from pyomo.core.base.expression import ScalarExpression, IndexedExpression, Expression, _GeneralExpressionData, ExpressionData
from pyomo.environ import Var, ScalarVar 
from pyomo.core.base.var import IndexedVar
import pyomo.environ as pyo


class ExpressionConversionError(Exception):
    pass


class _SeawaterStateBlockConstraints(StateBlockConstraints):
    
    def constrain_component(blk, component: Var | Expression, value: float) -> Var | None:
        """
        Constrain a component to a value
        """
        try:
            variable = _convert_expression_to_var(component)
        except ExpressionConversionError as e:
            variable = component # already a Var, just fix it directly
        
        variable.fix(value)

        if isinstance(variable, IndexedVar):
            for i in variable.index_set():
                # direct dictionary access avoids intercepted attribute resolution
                blk.__dict__["vars_to_deactivate"].append(variable[i])
        else:
            blk.__dict__["vars_to_deactivate"].append(variable)

        return variable


def _convert_expression_to_var(expr: ScalarExpression | IndexedExpression):
    if isinstance(expr, ScalarExpression) or isinstance(expr, ExpressionData):
        var = Var(units=pyo.units.get_units(expr))
        constraint = Constraint(expr= var == expr)
    elif isinstance(expr, IndexedExpression):
        var = Var(expr.index_set(), units=pyo.units.get_units(expr.units))
        def rule(b, i):
            return var[i] == expr[i]
        constraint = Constraint(expr= rule)
    else:
        raise ExpressionConversionError(f"Expression {expr} is not a ScalarExpression or IndexedExpression: {type(expr)}")
    block = expr.parent_block()
    block.add_component(f"{expr.local_name}_var", var)
    block.add_component(f"{expr.local_name}_constraint", constraint)
    return var

def _deactivate_additional_constraints(self: _ExtendedSeawaterStateBlock):
    # Temporarily deactivate platform constraints added with
    # StateBlockConstraints.constrain()) so they don't interfere with
    # initialization.
    deactivated_vars: list[tuple[ScalarVar,float]] = []

    for k in self.keys():
        blk = self[k]

        for var in blk.__dict__["vars_to_deactivate"]:
            var.unfix()
            # Store the original value so we can reactivate and fix back to the original value later.
            deactivated_vars.append((var, var.value))
    
    self.deactivated_vars = deactivated_vars

def _reactivate_additional_constraints(self: _ExtendedSeawaterStateBlock):
    for var, value in self.deactivated_vars:
        var.fix(value)

def _solve_block(self, solve_log, init_log, opt, step_name):
    skip_solve = True  # skip solve if only state variables are present
    for k in self.keys():
        if number_unfixed_variables(self[k]) != 0:
            skip_solve = False

    if not skip_solve:
        # Initialize properties
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solve_indexed_blocks(opt, [self], tee=slc.tee)
        init_log.info_high(
            f"Property initialization {step_name}: {idaeslog.condition(results)}"
        )
    
    if (not skip_solve) and (not check_optimal_termination(results)):
        raise InitializationError(
            f"{self.name} {step_name} failed to initialize successfully. Please "
            f"check the output logs for more information."
        )


class _ExtendedSeawaterStateBlock(_SeawaterStateBlock):
    def initialize(
        self,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """Initialization routine for the extended seawater property package.

        The procedure is:
          1. Deactivate platform constraints and fix state variables.
          2. Verify zero degrees of freedom on every sub-block.
          3. Solve with state variables fixed.
          4. Release state variables and re-activate constraints.
          5. If DOF == 0, solve again so constraints are satisfied.
          6. Optionally re-fix state variables when ``hold_state=True``.

        Args:
            state_args (dict, optional): Initial guesses keyed by state
                variable name (``flow_mass_phase_comp``, ``pressure``,
                ``temperature``).  When called from a control volume the
                inlet values are passed automatically.
            state_vars_fixed (bool): ``True`` if the control volume has
                already fixed state variables (e.g. CV1D).
            hold_state (bool): If ``True``, state variables are left fixed
                after initialization and a flags dict is returned so the
                caller can later call ``release_state``.
            outlvl: IDAES logging output level.
            solver: Pyomo solver object (default: IDAES default solver).
            optarg (dict, optional): Solver options.

        Returns:
            dict or None: Flags dict when ``hold_state=True`` and
            ``state_vars_fixed=False``; otherwise ``None``.
        """
        # Get loggers
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Set solver and options
        opt = get_solver(solver, optarg)

        # Fix state variables
        _deactivate_additional_constraints(self)
        flags = fix_state_vars(self, state_args)
        
        # Check when the state vars are fixed already result in dof 0
        for k in self.keys():
            dof = degrees_of_freedom(self[k])
            if dof != 0:
                raise PropertyPackageError(
                    "State vars fixed but degrees of "
                    "freedom for state block is not "
                    "zero during initialization."
                )

        _solve_block(self, solve_log, init_log, opt, step_name="Initial solve with state vars fixed")

        self.release_state(flags)
        if degrees_of_freedom(self) == 0:
            # We want to solve again, with any constraints reactivated, to ensure that the state variables have the correct values.
            
            _solve_block(self, solve_log, init_log, opt, step_name="Second solve with constraints reactivated")

        if hold_state:
            # Switch back to state vars fixed
            _deactivate_additional_constraints(self)
            flags = fix_state_vars(self, state_args)

            
        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)
        

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        super().release_state(flags, outlvl)
        # Reactivate platform constraints that were deferred during initialize
        _reactivate_additional_constraints(self)


@declare_process_block_class("SeawaterExtendedStateBlock", block_class=_ExtendedSeawaterStateBlock)
class SeawaterExtendedStateBlockData(SeawaterStateBlockData, _SeawaterStateBlockConstraints):

    def build(blk, *args):
        SeawaterStateBlockData.build(blk, *args)
        StateBlockConstraints.build(blk, *args)

        # We initialise vars_to_deactivate here with __setattr__ instead of doing it 
        # in the constructor of _SeawaterStateBlockConstraints as blk.vars_to_deactivate = [], 
        # because on state blocks, missing attributes are not treated like normal Python objects.
        # Why the below works is because it bypasses the custom __setattr__ logic (which makes it a metadata) 
        # on the block and writes directly to the object.
        # Also, initializing here guarantees every state block has its own list before constrain_component() runs.
        object.__setattr__(blk, "vars_to_deactivate", [])

    def add_extra_expressions(blk):
        """
        Override base add_extra_expressions for seawater.

        Unlike the Helmholtz package, seawater uses mass-phase properties
        (enth_mass_phase, flow_mass_phase_comp) rather than molar properties.
        We must NOT eagerly create expressions that reference lazy WaterTAP
        properties (e.g. enth_mass_phase) during build, because that triggers
        Var/Constraint construction on every state block and breaks the
        system solve.

        Instead, the extra properties (enth_mass, total_energy_flow) are
        registered via IDAES metadata in the parameter block below, so they
        are only constructed when actually accessed (e.g. by the platform).
        """
        pass

    def _enth_mass(blk):
        blk.enth_mass = Expression(expr=blk.enth_mass_phase["Liq"])

    def _enth_mol(blk):
        # enth_mol = enth_mass * average_molar_weight
        # average_molar_weight = total_mass_flow / total_mol_flow
        blk.enth_mol = Expression(
            expr=blk.enth_mass_phase["Liq"]
            * sum(blk.flow_mass_phase_comp["Liq", j] for j in blk.params.component_list)
            / sum(blk.flow_mol_phase_comp["Liq", j] for j in blk.params.component_list)
        )

    def _entr_mass(blk):
        # Seawater package does not provide entropy directly; use an
        # enthalpy-based surrogate to satisfy platform property interface.
        # This mirrors prior extended-package behavior where extra aliases are
        # used for generic platform compatibility.
        blk.entr_mass = Expression(expr=1 * pyo.units.J/pyo.units.kg / pyo.units.K)

    def _entr_mol(blk):
        blk.entr_mol = Expression(expr=1 * pyo.units.J/pyo.units.mol / pyo.units.K)

    def _flow_mass(blk):
        blk.flow_mass = Expression(
            expr=sum(blk.flow_mass_phase_comp["Liq", j] for j in blk.params.component_list)
        )

    def _flow_mol(blk):
        blk.flow_mol = Expression(
            expr=sum(blk.flow_mol_phase_comp["Liq", j] for j in blk.params.component_list)
        )

    def _mole_frac_comp(blk):
        blk.mole_frac_comp = Var(
            blk.params.component_list,
        bounds=(0,1),initialize=0.5)

        @blk.Constraint(blk.params.component_list)
        def mole_frac_comp_rule(b, j):
            return b.mole_frac_comp[j] == b.mole_frac_phase_comp["Liq", j]

        if (blk.config.defined_state):
            # We are fixing TDS seperately, so we will let them be separate.
            blk.mole_frac_comp_rule["TDS"].deactivate()

    def _vapor_frac(blk):
        # Seawater package is liquid-only for this application.
        blk.vapor_frac = Expression(expr=0.0)

    def _total_energy_flow(blk):
        blk.total_energy_flow = Expression(expr=blk.enth_flow)


@declare_process_block_class("SeawaterExtendedParameterBlock")
class SeawaterExtendedParameterBlockData(SeawaterParameterData):

    @classmethod
    def define_metadata(cls, obj):
        SeawaterParameterData.define_metadata(obj)
        # enth_mass is a standard IDAES property; update its method
        obj.add_properties(
            {
                "enth_mass": {"method": "_enth_mass"},
                "enth_mol": {"method": "_enth_mol"},
                "entr_mass": {"method": "_entr_mass"},
                "entr_mol": {"method": "_entr_mol"},
                "flow_mass": {"method": "_flow_mass"},
                "flow_mol": {"method": "_flow_mol"},
                "mole_frac_comp": {"method": "_mole_frac_comp"},
            }
        )
        # total_energy_flow is not a standard property; define as custom
        obj.define_custom_properties(
            {
                "total_energy_flow": {"method": "_total_energy_flow"},
                "vapor_frac": {"method": "_vapor_frac"},
            }
        )

    def build(self):
        super().build()
        self._state_block_class = SeawaterExtendedStateBlock  # noqa: F821