from idaes.models.properties.activity_coeff_models.activity_coeff_prop_pack import (
    ActivityCoeffParameterData,
    ActivityCoeffStateBlockData,
    _ActivityCoeffStateBlock,
)

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
)

from pyomo.environ import (Expression, Param, Var, Constraint, units as pyunits)


@declare_process_block_class("ExtendedActivityCoeffParameterBlock")
class ExtendedActivityCoeffParameterData(ActivityCoeffParameterData):
    def build(self):
      super().build()
      self._state_block_class = ExtendedActivityCoeffStateBlock # noqa: F821

    @classmethod
    def define_metadata(cls, obj):
            ActivityCoeffParameterData.define_metadata(obj)
            """Define properties supported and units."""
            obj.add_properties(
                    {
                        "enth_mol_comp": {"method": "_enth_mol_comp", "units": "J/mol"},
                        "enth_mol": {"method": "_enth_mol", "units": "J/mol"},
                        "enth_mass_comp": {"method": "_enth_mass_comp", "units": "J/kg"},
                        "enth_mass": {"method": "_enth_mass", "units": "J/kg"},
                        "flow_mol_comp": {"method": "_flow_mol_comp", "units": "mol/s"},
                        "flow_mass_comp": {"method": "_flow_mass_comp", "units": "kg/s"},
                        "flow_mass": {"method": "_flow_mass", "units": "kg/s"},
                        "entr_mol_comp": {"method": "_entr_mol_comp", "units": "J/mol"},
                        "entr_mol": {"method": "_entr_mol", "units": "J/mol"},
                        "entr_mass_comp": {"method": "_entr_mass_comp", "units": "J/kg"},
                        "entr_mass": {"method": "_entr_mass", "units": "J/kg"},
                        "vapor_frac": {"method": "_vapor_frac", "units": None},
                        "flow_vol": {"method": "_flow_vol", "units": "m^3/s"},
                        "total_energy_flow": {"method": "_total_energy_flow", "units": "J/s"},
                        "phase_frac_mol": {"method": "_phase_frac_mol", "units": None},
                        "flow_mass_phase": {"method": "_flow_mass_phase_comp", "units": "kg/s"},
                        "phase_frac_mass": {"method": "_phase_frac_mass", "units": None},
                        "mass_frac_phase_comp": {"method": "_mass_frac_phase_comp", "units": None},
                    }
            )

@declare_process_block_class("ExtendedActivityCoeffStateBlock", block_class=_ActivityCoeffStateBlock)
class ExtendedActivityCoeffStateBlockData(ActivityCoeffStateBlockData):

    def build(self):
        super().build()
        # we want to always create these properties (We said the magic words)
        self.enth_mol_phase
        self.entr_mol_phase

        #Ahuora layout mol_frac_comp, flow_mol, temperature, pressure, flow_mass, vapor_frac, enth_mol , enth_mass, entr_mol, entr_mass, total_energy_flow, flow_vol
    def _flow_mol_comp(self):
      def _rule_flow_mol_comp(model, i):
          return model.mole_frac_comp[i] * model.flow_mol
      self.flow_mol_comp = Expression(
          self.params.component_list,
          rule=_rule_flow_mol_comp,
      )

    def _flow_mass_comp(self):
      def _rule_flow_mass_comp(model, i):
          return model.flow_mol_comp[i] * model.params.mw_comp[i]
      self.flow_mass_comp = Expression(
          self.params.component_list,
          rule=_rule_flow_mass_comp,
      )
    
    def _flow_mass(self):
      def _rule_flow_mass(model):
          return sum(model.flow_mass_comp[i] for i in model.params.component_list)
      self.flow_mass = Expression(
          rule=_rule_flow_mass,
      )

    def _vapor_frac(self):
        def _rule_vapor_frac(model):
            return model.flow_mol_phase['Vap'] / model.flow_mol
        self.vapor_frac = Expression(
            rule=_rule_vapor_frac,
        )
    def _phase_frac_mol(self):
        def _rule_phase_frac_mol(model, p):
            return model.flow_mol_phase[p] / model.flow_mol
        self.phase_frac_mol = Expression(
            self.params.phase_list,
            rule=_rule_phase_frac_mol,
        )

    def _enth_mol(self):
      def _rule_enth_mol(model):
          return sum(model.phase_frac_mol[p]*model.enth_mol_phase[p] for p in model.params.phase_list)
      self.enth_mol = Expression(
          rule=_rule_enth_mol
      )

    def _weighted_Molecular_weight(self):
        def _rule_weighted_Molecular_weight(model):
            return sum(model.mole_frac_comp[c] * model.params.mw_comp[c] for c in model.params.component_list)
        self.weighted_Molecular_weight = Expression(
            rule=_rule_weighted_Molecular_weight,
        )

    def _flow_mass_phase(self):
        def _rule_flow_mass_phase(model, p):
            return model.flow_mol_phase[p] * self.weighted_Molecular_weight
        self.flow_mass_phase = Expression(
            self.params.phase_list,
            rule=_rule_flow_mass_phase
        )
    def _phase_frac_mass(self):
        def _rule_phase_frac_mass(model, p):
            return model.flow_mass_phase[p] / model.flow_mass
        self.phase_frac_mass = Expression(
            self.params.phase_list,
            rule=_rule_phase_frac_mass
        )
    def _mass_frac_phase_comp(self):
        def _rule_mass_frac_phase_comp(model, p, i):
            return model.flow_mass_phase_comp[p, i] / model.flow_mass
        self.mass_frac_phase_comp = Expression(
            self.params.phase_list,
            self.params.component_list,
            rule=_rule_mass_frac_phase_comp,
        )

    def enth_mass_phase(self):
        def _rule_enth_mass_phase(model,p):
            return model.enth_mol_phase[p] / model.weighted_Molecular_weight
        self.enth_mass_phase = Expression(
            self.params.phase_list,
            rule=_rule_enth_mass_phase
        )

    def _enth_mass(self):
      def _rule_enth_mass(model):
          return sum(model.phase_frac_mass[p]*model.enth_mass_phase[p] for p in model.params.phase_list)
                
      self.enth_mass = Expression(
          rule=_rule_enth_mass
      )

    def _entr_mol(self):
      def _rule_entr_mol(model):
          return sum(model.phase_frac_mol[p]*model.entr_mol_phase[p] for p in model.params.phase_list)
      self.entr_mol = Expression(
          rule=_rule_entr_mol
      )

    def _entr_mass_phase(self):
        def _rule_entr_mass_phase(model,p):
            return model.entr_mol_phase[p] / model.weighted_Molecular_weight
        self.entr_mass_phase = Expression(
            self.params.phase_list,
            rule=_rule_entr_mass_phase
        )  

    def _entr_mass(self):
      def _rule_entr_mass(model):
          return sum(model.phase_frac_mass[p]*model.entr_mass_phase[p] for p in model.params.phase_list)
      self.entr_mass = Expression(
          rule=_rule_entr_mass
      )

    def _total_energy_flow(self):
        def _rule_total_energy_flow(model):
            return model.flow_mol * model.enth_mol
        self.total_energy_flow = Expression(
            rule=_rule_total_energy_flow,
        )
    
    def _flow_vol(self):
      def _rule_flow_vol(model):
          return 0 * pyunits.m**3 / pyunits.s # Placeholder until we have density data
      self.flow_vol = Expression(
          rule=_rule_flow_vol,
      )
