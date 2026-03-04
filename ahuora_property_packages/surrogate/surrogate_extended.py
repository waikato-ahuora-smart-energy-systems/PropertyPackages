from pyomo.environ import (
    Expression,
    SolverFactory,
    check_optimal_termination,
    units,
    Constraint,
    Var,
)
from idaes.core import (
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
)
from ahuora_property_packages.utils.fix_state_vars import fix_state_vars
from ahuora_property_packages.base.state_block_constraints import StateBlockConstraints
from idaes.core import declare_process_block_class
from idaes.core.util.initialization import get_solver, fix_state_vars, solve_indexed_blocks
from idaes.core.util.model_statistics import degrees_of_freedom, number_unfixed_variables
from idaes.core.util.exceptions import InitializationError, PropertyPackageError
import idaes.logger as idaeslog
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate, PysmoPolyTrainer
from idaes.core.base.components import Solvent, Solute
from idaes.core.base.phases import LiquidPhase
from pyomo.environ import units as pyunits


def _deactivate_additional_constraints(self):
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
    return deactivated_constraints

def _reactivate_additional_constraints(self, deactivated_constraints):
    for k, cons_list in deactivated_constraints.items():
        for c in cons_list:
            c.activate()

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


class _ExtendedSurrogateStateBlock(StateBlock):
    def initialize(
        self,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """Initialization routine for the extended surrogate property package.

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
        deactivated_constraints = _deactivate_additional_constraints(self)
        self._deactivated_platform_constraints = deactivated_constraints
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
            deactivated_constraints = _deactivate_additional_constraints(self)
            self._deactivated_platform_constraints = deactivated_constraints
            flags = fix_state_vars(self, state_args)

            
        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)
        

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        # Reactivate platform constraints that were deferred during initialize
        if hasattr(self, "_deactivated_platform_constraints"):
            _reactivate_additional_constraints(self, self._deactivated_platform_constraints)
            del self._deactivated_platform_constraints
        super().release_state(flags, outlvl)



surrogate = PysmoSurrogate.load_from_file("coolprop_surrogate.json")

@declare_process_block_class("SurrogateExtendedStateBlock", block_class=_ExtendedSurrogateStateBlock)
class SurrogateExtendedStateBlockData(StateBlockData, StateBlockConstraints):

    def build(blk, *args):
        StateBlockData.build(blk, *args)
        blk.temperature = Var( units = pyunits.K)
        blk.pressure = Var( units = pyunits.Pa)
        blk.flow_mol = Var( units = pyunits.mol/pyunits.s)
        blk.flow_mass = Var( units = pyunits.kg/pyunits.s)
        blk.enth_mass = Var( units = pyunits.J/pyunits.kg)
        blk.enth_mol = Var( units = pyunits.J/pyunits.mol)
        blk.entr_mass = Var( units = pyunits.J/pyunits.kg/pyunits.K)
        blk.entr_mol = Var( units = pyunits.J/pyunits.mol/pyunits.K)
        blk.dynamic_viscosity = Var( units= pyunits.Pa*pyunits.s)
        blk.kinematic_viscosity = Var( units = pyunits.m**2/pyunits.s)
        blk.specific_volume = Var( units = pyunits.m**3/pyunits.kg)
        blk.mw = Var( units = pyunits.kg/pyunits.mol)
        blk.mw.fix(0.01801528) # No need to calculate as this is a constant, and makes the problem poorly defined.
        blk.unused_mw = Var( units = pyunits.kg/pyunits.mol)


        @blk.Constraint()
        def _rule_flow_mass(b, t):
            return b.flow_mass == b.flow_mol * b.mw

        @blk.Constraint()
        def _rule_enth_mass(b, t):
            return b.enth_mass == b.enth_mol / b.mw
        
        @blk.Constraint()
        def _rule_entr_mass(b, t):
            return b.entr_mass == b.entr_mol / b.mw

        # Add the surrogate model block
        blk.surrogate = SurrogateBlock()
        blk.surrogate.build_model(surrogate, input_vars=[blk.temperature, blk.pressure], 
                                  output_vars=[
                                      blk.entr_mass,
                                      blk.enth_mass,
                                      blk.dynamic_viscosity,
                                      blk.kinematic_viscosity,
                                      blk.unused_mw, # Molar mass
                                      blk.specific_volume
        ])

        StateBlockConstraints.build(blk, *args)
    
    def define_state_vars(self):
        """Define state vars."""
        return {
            "flow_mol": self.flow_mol,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }


@declare_process_block_class("SurrogateExtendedParameterBlock")
class SurrogateExtendedParameterBlockData(PhysicalParameterBlock):
    def build(self):
        super().build()

        # components
        self.H2O = Solvent()

        # phases
        self.Liq = LiquidPhase()
        self._state_block_class = SurrogateExtendedStateBlock # noqa: F821
    

    @classmethod
    def define_metadata(cls, pcm):
        """Define properties supported and units."""
        pcm.add_properties(
            {
                "flow_mass_phase_comp": {"method": None},
                "flow_mass": {"method": None},
                "flow_mol": {"method": None},
                "enth_mol": {"method": None},
                "temperature": {"method": None},
                "pressure": {"method": None},
                "mass_frac_phase_comp": {"method": "_mass_frac_phase_comp"},
                "dens_mass_phase": {"method": "_dens_mass_phase"},
                "flow_vol_phase": {"method": "_flow_vol_phase"},
                "flow_vol": {"method": "_flow_vol"},
                "conc_mass_phase_comp": {"method": "_conc_mass_phase_comp"},
                "flow_mol_phase_comp": {"method": "_flow_mol_phase_comp"},
                "mole_frac_phase_comp": {"method": "_mole_frac_phase_comp"},
                "molality_phase_comp": {"method": "_molality_phase_comp"},
                "visc_d_phase": {"method": "_visc_d_phase"},
                "pressure_osm_phase": {"method": "_pressure_osm_phase"},
                "energy_density_phase": {"method": "_energy_density_phase"},
                "enth_mass_phase": {"method": "_enth_mass_phase"},
                "pressure_sat": {"method": "_pressure_sat"},
                "cp_mass_phase": {"method": "_cp_mass_phase"},
                "therm_cond_phase": {"method": "_therm_cond_phase"},
                "diffus_phase_comp": {"method": "_diffus_phase_comp"},
            }
        )

        # pcm.define_custom_properties(
        #     {
        #         "dens_mass_solvent": {"method": "_dens_mass_solvent"},
                
        #     }
        # )

        pcm.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )
