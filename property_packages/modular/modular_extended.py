from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import ConfigurationError
from idaes.models.properties.modular_properties.base.generic_property import (
    _GenericStateBlock,
    GenericParameterData,
    GenericStateBlockData,
)
import idaes.models.properties.modular_properties.base.utility as utility

from property_packages.base.state_block_constraints import StateBlockConstraints
from property_packages.utils.add_extra_expressions import add_extra_expressions
from property_packages.utils.fix_state_vars import fix_state_vars


# TODO: remove these imports once https://github.com/IDAES/idaes-pse/pull/1554 is resolved
from pyomo.environ import (
    Block,
    check_optimal_termination,
    Constraint,
    exp,
    Expression,
    log,
    Set,
    Param,
    value,
    Var,
    units as pyunits,
    Reference,
)
from pyomo.common.config import ConfigBlock, ConfigDict, ConfigValue, In, Bool
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialFlowBasis,
    ElectrolytePropertySet,
)
from idaes.core.base.components import Component, __all_components__
from idaes.core.base.phases import (
    Phase,
    AqueousPhase,
    LiquidPhase,
    VaporPhase,
    __all_phases__,
)
from idaes.core.util.initialization import (
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_activated_constraints,
)
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyPackageError,
    PropertyNotSupportedError,
    InitializationError,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.initialization.initializer_base import InitializerBase

from idaes.models.properties.modular_properties.base.generic_reaction import (
    equil_rxn_config,
)
from idaes.models.properties.modular_properties.base.utility import (
    get_method,
    get_phase_method,
    GenericPropertyPackageError,
    StateIndex,
    identify_VL_component_list,
    estimate_Tbub,
    estimate_Tdew,
    estimate_Pbub,
    estimate_Pdew,
)
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    LogBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.henry import HenryType
from idaes.models.properties.modular_properties.base.generic_property import (
    _initialize_critical_props,
    _init_Tbub,
    _init_Tdew,
    _init_Pbub,
    _init_Pdew,
)


# increase max iterations for estimating values for bubble, dew, and
# critical point initialization (ie. temperature_bubble, temperature_dew)
utility.MAX_ITER = 1000


class _ExtendedGenericStateBlock(_GenericStateBlock):

    def initialize(blk, *args, **kwargs):
        flag_dict = fix_state_vars(blk, kwargs.get("state_args", None))

        # Set state_vars_fixed to True to avoid fixing state variables
        # during the initialize method, since this would overdefine
        # the block if we are using constraints
        kwargs["state_vars_fixed"] = True

        # TODO: replace this with
        # super().initialize(*args, **kwargs)
        # once https://github.com/IDAES/idaes-pse/pull/1554 has
        # been resolved and released (allows using constraints to
        # define state variables during initialization)
        blk._custom_super_initialize(*args, **kwargs)

        if kwargs.get("hold_state") is True:
            return flag_dict
        else:
            blk.release_state(flag_dict)
    

    def _custom_super_initialize(
        blk,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine for property package.
        Keyword Arguments:
            state_args : a dict of initial values for the state variables
                    defined by the property package.
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed by the
                                       control volume 1D. Control volume 0D
                                       does not fix the state vars, so will
                                       be False if this state block is used
                                       with 0D blocks.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         - True - states variables are not unfixed, and
                                 a dict of returned containing flags for
                                 which states were fixed during
                                 initialization.
                        - False - state variables are unfixed after
                                 initialization by calling the
                                 release_state method
        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="properties")

        init_log.info("Starting initialization")

        res = None

        for k in blk.values():
            # Deactivate the constraints specific for outlet block i.e.
            # when defined state is False
            if k.config.defined_state is False:
                try:
                    k.sum_mole_frac_out.deactivate()
                except AttributeError:
                    pass

                if hasattr(k, "inherent_equilibrium_constraint") and (
                    not k.params._electrolyte
                    or k.params.config.state_components == StateIndex.true
                ):
                    k.inherent_equilibrium_constraint.deactivate()

        # Fix state variables if not already fixed
        if state_vars_fixed is False:
            flag_dict = fix_state_vars(blk, state_args)
            # Confirm DoF for sanity
            for k in blk.values():
                if k.always_flash:
                    # If not always flash, DoF is probably less than zero
                    # We will handle this elsewhere
                    dof = degrees_of_freedom(k)
                    if dof != 0:
                        raise BurntToast(
                            "Degrees of freedom were not zero [{}] "
                            "after trying to fix state variables. "
                            "Something broke in the generic property "
                            "package code - please inform the IDAES "
                            "developers.".format(dof)
                        )
        else:
            # When state vars are fixed, check that DoF is 0
            for k in blk.values():
                if degrees_of_freedom(k) != 0:
                    # PYLINT-TODO
                    # pylint: disable-next=broad-exception-raised
                    raise Exception(
                        "State vars fixed but degrees of "
                        "freedom for state block is not zero "
                        "during initialization."
                    )

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # If present, initialize bubble, dew , and critical point calculations
        for k in blk.values():
            T_units = k.params.get_metadata().default_units.TEMPERATURE

            # List of bubble and dew point constraints
            cons_list = [
                "eq_pressure_dew",
                "eq_pressure_bubble",
                "eq_temperature_dew",
                "eq_temperature_bubble",
                "eq_mole_frac_tbub",
                "eq_mole_frac_tdew",
                "eq_mole_frac_pbub",
                "eq_mole_frac_pdew",
                "log_mole_frac_tbub_eqn",
                "log_mole_frac_tdew_eqn",
                "log_mole_frac_pbub_eqn",
                "log_mole_frac_pdew_eqn",
                "mole_frac_comp_eq",
                "log_mole_frac_comp_eqn",
            ]

            # Critical point
            with k.lock_attribute_creation_context():
                # Only need to look for one, as it is all-or-nothing
                if hasattr(k, "pressure_crit"):
                    # Initialize critical point properties
                    _initialize_critical_props(k)
                    # Add critical point constraints to cons_list
                    cons_list += k.list_critical_property_constraint_names()

            # Bubble temperature initialization
            if hasattr(k, "_mole_frac_tbub"):
                _init_Tbub(k, T_units)

            # Dew temperature initialization
            if hasattr(k, "_mole_frac_tdew"):
                _init_Tdew(k, T_units)

            # Bubble pressure initialization
            if hasattr(k, "_mole_frac_pbub"):
                _init_Pbub(k)

            # Dew pressure initialization
            if hasattr(k, "_mole_frac_pdew"):
                _init_Pdew(k)

            # Solve bubble, dew, and critical point constraints
            for c in k.component_objects(Constraint):
                # Deactivate all constraints not associated with bubble, dew,
                # or critical points
                if c.local_name not in cons_list:
                    c.deactivate()

        # If StateBlock has active constraints (i.e. has bubble, dew, or critical
        # point calculations), solve the block to converge these
        n_cons = 0
        dof = 0
        for k in blk.values():
            n_cons += number_activated_constraints(k)
            dof += degrees_of_freedom(k)
        if n_cons > 0:
            if dof > 0:
                raise InitializationError(
                    f"{blk.name} Unexpected degrees of freedom during "
                    f"initialization at bubble, dew, and critical point step: {dof}."
                )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
            init_log.info(
                "Bubble, dew, and critical point initialization: {}.".format(
                    idaeslog.condition(res)
                )
            )
        # ---------------------------------------------------------------------
        # Calculate _teq if required
        # Using iterator k outside of for loop - this should be OK as we just need
        # a valid StateBlockData an assume they are all the same.
        if k.params.config.phases_in_equilibrium is not None and (
            not k.config.defined_state or k.always_flash
        ):
            for k in blk.values():
                for pp in k.params._pe_pairs:
                    k.params.config.phase_equilibrium_state[pp].calculate_teq(k, pp)

            init_log.info("Equilibrium temperature initialization completed.")

        # ---------------------------------------------------------------------
        # Initialize flow rates and compositions
        for k in blk.values():

            k.params.config.state_definition.state_initialization(k)

            if k.params._electrolyte:
                if k.params.config.state_components == StateIndex.true:
                    # First calculate initial values for apparent species flows
                    for p, j in k.params.apparent_phase_component_set:
                        calculate_variable_from_constraint(
                            k.flow_mol_phase_comp_apparent[p, j],
                            k.true_to_appr_species[p, j],
                        )
                    # Need to calculate all flows before doing mole fractions
                    for p, j in k.params.apparent_phase_component_set:
                        sum_flow = sum(
                            k.flow_mol_phase_comp_apparent[p, jj]
                            for jj in k.params.apparent_species_set
                            if (p, jj) in k.params.apparent_phase_component_set
                        )
                        if value(sum_flow) == 0:
                            x = 1
                        else:
                            x = value(k.flow_mol_phase_comp_apparent[p, j] / sum_flow)
                        lb = k.mole_frac_phase_comp_apparent[p, j].lb
                        if lb is not None and x <= lb:
                            k.mole_frac_phase_comp_apparent[p, j].set_value(lb)
                        else:
                            k.mole_frac_phase_comp_apparent[p, j].set_value(x)
                elif k.params.config.state_components == StateIndex.apparent:
                    # First calculate initial values for true species flows
                    for p, j in k.params.true_phase_component_set:
                        calculate_variable_from_constraint(
                            k.flow_mol_phase_comp_true[p, j],
                            k.appr_to_true_species[p, j],
                        )
                    # Need to calculate all flows before doing mole fractions
                    for p, j in k.params.true_phase_component_set:
                        sum_flow = sum(
                            k.flow_mol_phase_comp_true[p, jj]
                            for jj in k.params.true_species_set
                            if (p, jj) in k.params.true_phase_component_set
                        )
                        if value(sum_flow) == 0:
                            x = 1
                        else:
                            x = value(k.flow_mol_phase_comp_true[p, j] / sum_flow)
                        lb = k.mole_frac_phase_comp_true[p, j].lb
                        if lb is not None and x <= lb:
                            k.mole_frac_phase_comp_true[p, j].set_value(lb)
                        else:
                            k.mole_frac_phase_comp_true[p, j].set_value(x)

            # If state block has phase equilibrium, use the average of all
            # _teq's as an initial guess for T
            if (
                k.params.config.phases_in_equilibrium is not None
                and isinstance(k.temperature, Var)
                and not k.temperature.fixed
            ):
                k.temperature.value = value(
                    sum(k._teq[i] for i in k.params._pe_pairs) / len(k.params._pe_pairs)
                )

        if outlvl > 0:  # TODO: Update to use logger Enum
            init_log.info("State variable initialization completed.")

        # ---------------------------------------------------------------------
        n_cons = 0
        dof = 0
        skip = False
        Tfix = {}  # In enth based state defs, need to also fix T until later
        for k, b in blk.items():
            if b.params.config.phase_equilibrium_state is not None and (
                not b.config.defined_state or b.always_flash
            ):
                if not b.temperature.fixed:
                    b.temperature.fix()
                    Tfix[k] = True
                for c in b.component_objects(Constraint):
                    # Activate common constraints
                    if c.local_name in (
                        "total_flow_balance",
                        "component_flow_balances",
                        "sum_mole_frac",
                        "phase_fraction_constraint",
                        "mole_frac_phase_comp_eq",
                        "mole_frac_comp_eq",
                    ):
                        c.activate()
                    if c.local_name == "log_mole_frac_phase_comp_eqn":
                        c.activate()
                        for p, j in b.params._phase_component_set:
                            calculate_variable_from_constraint(
                                b.log_mole_frac_phase_comp[p, j],
                                b.log_mole_frac_phase_comp_eqn[p, j],
                            )
                    elif c.local_name == "equilibrium_constraint":
                        # For systems where the state variables fully define the
                        # phase equilibrium, we cannot activate the equilibrium
                        # constraint at this stage.
                        if "flow_mol_phase_comp" not in b.define_state_vars():
                            c.activate()
                    elif getattr(c, "defining_state_var", False):
                        # Allow using a constraint to define a state var
                        # (rather than fixing it directly) by setting 
                        # c.defining_state_var = True. This can be used in
                        # conjunction with setting state_vars_fixed = True
                        c.activate()

                for pp in b.params._pe_pairs:
                    # Activate formulation specific constraints
                    b.params.config.phase_equilibrium_state[
                        pp
                    ].phase_equil_initialization(b, pp)

            n_cons += number_activated_constraints(b)
            dof += degrees_of_freedom(b)
            if degrees_of_freedom(b) < 0:
                # Skip solve if DoF < 0 - this is probably due to a
                # phase-component flow state with flash
                skip = True

        if n_cons > 0 and not skip:
            if dof > 0:
                raise InitializationError(
                    f"{blk.name} Unexpected degrees of freedom during "
                    f"initialization at phase equilibrium step: {dof}."
                )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
            init_log.info(
                "Phase equilibrium initialization: {}.".format(idaeslog.condition(res))
            )

        # ---------------------------------------------------------------------
        # Initialize other properties
        for k, b in blk.items():
            for c in b.component_objects(Constraint):
                # Activate all constraints except flagged do_not_initialize
                if c.local_name not in (
                    b.params.config.state_definition.do_not_initialize
                ):
                    c.activate()
            if k in Tfix:
                b.temperature.unfix()

            # Initialize log-form variables
            log_form_vars = [
                "act_phase_comp",
                "act_phase_comp_apparent",
                "act_phase_comp_true",
                "conc_mol_phase_comp",
                "conc_mol_phase_comp_apparent",
                "conc_mol_phase_comp_true",
                "mass_frac_phase_comp",
                "mass_frac_phase_comp_apparent",
                "mass_frac_phase_comp_true",
                "molality_phase_comp",
                "molality_phase_comp_apparent",
                "molality_phase_comp_true",
                "mole_frac_comp",  # Might have already been initialized
                "mole_frac_phase_comp",  # Might have already been initialized
                "mole_frac_phase_comp_apparent",
                "mole_frac_phase_comp_true",
                "pressure_phase_comp",
                "pressure_phase_comp_apparent",
                "pressure_phase_comp_true",
            ]

            for prop in log_form_vars:
                if b.is_property_constructed("log_" + prop):
                    comp = getattr(b, prop)
                    lcomp = getattr(b, "log_" + prop)
                    for k2, v in lcomp.items():
                        c = value(comp[k2])
                        if c <= 0:
                            c = 1e-8
                        lc = log(c)
                        v.set_value(value(lc))

        n_cons = 0
        dof = 0
        skip = False
        for k in blk.values():
            if degrees_of_freedom(k) < 0:
                # Skip solve if DoF < 0 - this is probably due to a
                # phase-component flow state with flash
                skip = True
            n_cons += number_activated_constraints(k)
            dof += degrees_of_freedom(k)
        if n_cons > 0 and not skip:
            if dof > 0:
                raise InitializationError(
                    f"{blk.name} Unexpected degrees of freedom during "
                    f"initialization at property initialization step: {dof}."
                )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
            init_log.info(
                "Property initialization: {}.".format(idaeslog.condition(res))
            )

        # ---------------------------------------------------------------------
        # Return constraints to initial state
        for k in blk.values():
            for c in k.component_objects(Constraint):
                if c.local_name in (k.params.config.state_definition.do_not_initialize):
                    c.activate()

        if res is not None and not check_optimal_termination(res):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        if state_vars_fixed is False:
            if hold_state is True:
                return flag_dict
            else:
                blk.release_state(flag_dict)

        init_log.info(
            "Property package initialization: {}.".format(idaeslog.condition(res))
        )


@declare_process_block_class(
    "GenericExtendedStateBlock", block_class=_ExtendedGenericStateBlock
)
class GenericExtendedStateBlockData(GenericStateBlockData, StateBlockConstraints):

    def build(self, *args):
        GenericStateBlockData.build(self, *args)
        StateBlockConstraints.build(self, *args)


@declare_process_block_class("GenericExtendedParameterBlock")
class GenericExtendedParameterData(GenericParameterData):
    def build(self):
        super().build()
        self._state_block_class = GenericExtendedStateBlock  # noqa: F821
