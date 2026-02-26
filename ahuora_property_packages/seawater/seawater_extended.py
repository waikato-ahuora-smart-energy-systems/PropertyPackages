from pyomo.environ import Expression, Constraint, check_optimal_termination
from watertap.property_models.seawater_prop_pack import SeawaterParameterData, SeawaterStateBlockData, _SeawaterStateBlock
from ahuora_property_packages.base.state_block_constraints import StateBlockConstraints
from idaes.core import declare_process_block_class
from watertap.core.solvers import get_solver
from idaes.core.util.initialization import fix_state_vars, solve_indexed_blocks
from idaes.core.util.model_statistics import degrees_of_freedom, number_unfixed_variables
from idaes.core.util.exceptions import InitializationError, PropertyPackageError
import idaes.logger as idaeslog


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
        """
        Initialization routine for property package.
        Keyword Arguments:
            state_args : Dictionary with initial guesses for the state vars
                         chosen. Note that if this method is triggered
                         through the control volume, and if initial guesses
                         were not provided at the unit model level, the
                         control volume passes the inlet values as initial
                         guess.The keys for the state_args dictionary are:
                         flow_mass_phase_comp : value at which to initialize
                                               phase component flows
                         pressure : value at which to initialize pressure
                         temperature : value at which to initialize temperature
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default={})
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed by the
                                       control volume 1D. Control volume 0D
                                       does not fix the state vars, so will
                                       be False if this state block is used
                                       with 0D blocks.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            solver : Solver object to use during initialization if None is provided
                     it will use the default solver for IDAES (default = None)
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
        # Get loggers
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Set solver and options
        opt = get_solver(solver, optarg)

        
        # Temporarily deactivate platform constraints added with
        # StateBlockConstraints.constrain()) so they don't interfere with
        # initialization.
        deactivated_constraints = {}
        for k in self.keys():
            blk = self[k]
            if hasattr(blk, "constraints"):
                deactivated = []
                for c in blk.constraints.component_objects(
                    Constraint, active=True, descend_into=False
                ):
                    c.deactivate()
                    deactivated.append(c)
                if deactivated:
                    deactivated_constraints[k] = deactivated

        # Fix state variables
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

        # ---------------------------------------------------------------------
        skip_solve = True  # skip solve if only state variables are present
        for k in self.keys():
            if number_unfixed_variables(self[k]) != 0:
                skip_solve = False

        if not skip_solve:
            # Initialize properties
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = solve_indexed_blocks(opt, [self], tee=slc.tee)
            init_log.info_high(
                "Property initialization: {}.".format(idaeslog.condition(results))
            )

        # Reactivate the platform constraints that were deactivated above.
        # If hold_state is True, defer reactivation until release_state(),
        # independent of state_vars_fixed. This avoids over-constraining
        # unit-model initialization solves when state vars are held fixed.
        if hold_state is True:
            self._deactivated_platform_constraints = deactivated_constraints
        else:
            for k, cons_list in deactivated_constraints.items():
                for c in cons_list:
                    c.activate()

        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)

        if (not skip_solve) and (not check_optimal_termination(results)):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please "
                f"check the output logs for more information."
            )

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        # Reactivate platform constraints that were deferred during initialize
        if hasattr(self, "_deactivated_platform_constraints"):
            for k, cons_list in self._deactivated_platform_constraints.items():
                for c in cons_list:
                    c.activate()
            del self._deactivated_platform_constraints
        super().release_state(flags, outlvl)


@declare_process_block_class("SeawaterExtendedStateBlock", block_class=_ExtendedSeawaterStateBlock)
class SeawaterExtendedStateBlockData(SeawaterStateBlockData, StateBlockConstraints):

    def build(blk, *args):
        SeawaterStateBlockData.build(blk, *args)
        StateBlockConstraints.build(blk, *args)

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
        blk.entr_mass = Expression(expr=blk.enth_mass)

    def _entr_mol(blk):
        blk.entr_mol = Expression(expr=blk.enth_mol)

    def _flow_mass(blk):
        blk.flow_mass = Expression(
            expr=sum(blk.flow_mass_phase_comp["Liq", j] for j in blk.params.component_list)
        )

    def _flow_mol(blk):
        blk.flow_mol = Expression(
            expr=sum(blk.flow_mol_phase_comp["Liq", j] for j in blk.params.component_list)
        )

    def _mole_frac_comp(blk):
        blk.mole_frac_comp = Expression(
            blk.params.component_list,
            rule=lambda b, j: b.mole_frac_phase_comp["Liq", j],
        )

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