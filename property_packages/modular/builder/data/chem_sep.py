"""
Methods for calculating pure component properties from:

http://www.chemsep.org/downloads/docs/book.htm

IDAES naming conventions followed for compatibility with modular property packages
"""

from pyomo.environ import log, exp, units as pyunits
from idaes.core.util.misc import set_param_from_config
from compounds.CompoundDB import get_compound
from pyomo.environ import units as u
from pyomo.environ import Var, value
from scipy.integrate import quad

name_map = {
    "enth_mol_form_vap_comp_ref": "HeatOfFormation",
    "entr_mol_form_vap_comp_ref": "AbsEntropy",
    "enth_mol_form_liq_comp_ref": 0,
    "entr_mol_form_liq_comp_ref": 0,
    "mw": "MolecularWeight",
    "omega": "AcentricityFactor",
    "pressure_crit": "CriticalPressure",
    "temperature_crit": "CriticalTemperature",
}

class ChemSep(object):

    _cached_components = {}

    @staticmethod
    def get_parameter_value(name, param):
        map_prop = name_map[param]
        if map_prop is None:
            raise Exception("Unrecognized parameter in ChemSep wrapper")
        else:
            return ChemSep._get_property(name, map_prop)

    class cp_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            ChemSep._create_params(cobj, "cp_mol_ig_comp_coeff")
            ChemSep._set_params(cobj, "cp_mol_ig_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            T = pyunits.convert(T, to_units=pyunits.K)
            cp = (
                cobj.cp_mol_ig_comp_coeff_A +
                cobj.cp_mol_ig_comp_coeff_B*T +
                cobj.cp_mol_ig_comp_coeff_C*T**2 +
                cobj.cp_mol_ig_comp_coeff_D*T**3 +
                cobj.cp_mol_ig_comp_coeff_E*T**4
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(cp, units.HEAT_CAPACITY_MOLE)

    class enth_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            # Required for calculating entropy
            if not hasattr(cobj, "cp_mol_ig_comp_coeff_A"):
                ChemSep.cp_mol_ig_comp.build_parameters(cobj)

            if cobj.parent_block().config.include_enthalpy_of_formation:

                units = cobj.parent_block().get_metadata().derived_units

                cobj.enth_mol_form_vap_comp_ref = Var(
                    doc="Vapor phase molar heat of formation",
                    units=units.ENERGY_MOLE,
                )
                set_param_from_config(cobj, param="enth_mol_form_vap_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific enthalpy
            T = pyunits.convert(T, to_units=pyunits.K)
            Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

            units = b.params.get_metadata().derived_units

            h_form = (
                cobj.enth_mol_form_vap_comp_ref
                if b.params.config.include_enthalpy_of_formation
                else 0 * units.ENERGY_MOLE
            )

            h = (
                pyunits.convert(
                    (cobj.cp_mol_ig_comp_coeff_A * (T - Tr)
                     + (cobj.cp_mol_ig_comp_coeff_B / 2) * (T**2 - Tr**2)
                     + (cobj.cp_mol_ig_comp_coeff_C / 3) * (T**3 - Tr**3)
                     + (cobj.cp_mol_ig_comp_coeff_D / 4) * (T**4 - Tr**4)
                     + (cobj.cp_mol_ig_comp_coeff_E / 5) * (T**5 - Tr**5)),
                    units.ENERGY_MOLE,
                ) + h_form
            )

            return h

    class entr_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            # Required for calculating enthalpy
            if not hasattr(cobj, "cp_mol_ig_comp_coeff_A"):
                ChemSep.cp_mol_ig_comp.build_parameters(cobj)

            units = cobj.parent_block().get_metadata().derived_units

            cobj.entr_mol_form_vap_comp_ref = Var(
                doc="Vapor phase molar entropy of formation",
                units=units.ENTROPY_MOLE,
            )

            set_param_from_config(cobj, param="entr_mol_form_vap_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy
            T = pyunits.convert(T, to_units=pyunits.K)
            Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

            units = b.params.get_metadata().derived_units

            s = (
                pyunits.convert(
                    (cobj.cp_mol_ig_comp_coeff_A * log(T / Tr)
                     + cobj.cp_mol_ig_comp_coeff_B * (T - Tr)
                     + (cobj.cp_mol_ig_comp_coeff_C / 2) * (T ** 2 - Tr ** 2)
                     + (cobj.cp_mol_ig_comp_coeff_D / 3) * (T ** 3 - Tr ** 3)
                     + (cobj.cp_mol_ig_comp_coeff_E / 4) * (T ** 4 - Tr ** 4)),
                    units.ENTROPY_MOLE,
                ) + cobj.entr_mol_form_vap_comp_ref
            )

            return s

    class pressure_sat_comp:
        @staticmethod
        def build_parameters(cobj):
            ChemSep._create_params(cobj, "pressure_sat_comp_coeff")
            ChemSep._set_params(cobj, "pressure_sat_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T, dT=False):
            psat = (
                    exp(
                        cobj.pressure_sat_comp_coeff_A - cobj.pressure_sat_comp_coeff_B /
                        (pyunits.convert(T, to_units=pyunits.K) + cobj.pressure_sat_comp_coeff_C))
                    ) * pyunits.Pa

            units = b.params.get_metadata().derived_units
            return pyunits.convert(psat, to_units=units.PRESSURE)
    
    class cp_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            ChemSep._create_params(cobj, "cp_mol_liq_comp_coeff")
            ChemSep._set_params(cobj, "cp_mol_liq_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity
            T = pyunits.convert(T, to_units=pyunits.K)

            cp = (
                cobj.cp_mol_liq_comp_coeff_A +
                exp(
                    cobj.cp_mol_liq_comp_coeff_B / T
                    + cobj.cp_mol_liq_comp_coeff_C
                    + cobj.cp_mol_liq_comp_coeff_D * T
                    + cobj.cp_mol_liq_comp_coeff_E * T**2
                )
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(cp, units.HEAT_CAPACITY_MOLE)

    class enth_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_liq_comp_coeff_A"):
                ChemSep.cp_mol_liq_comp.build_parameters(cobj)

            if cobj.parent_block().config.include_enthalpy_of_formation:
                units = cobj.parent_block().get_metadata().derived_units

                cobj.enth_mol_form_liq_comp_ref = Var(
                    doc="Liquid phase molar heat of formation @ Tref",
                    units=units.ENERGY_MOLE,
                )
                set_param_from_config(cobj, param="enth_mol_form_liq_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific enthalpy
            T = pyunits.convert(T, to_units=pyunits.K)
            Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

            units = b.params.get_metadata().derived_units

            integrand = lambda T: ChemSep.cp_mol_liq_comp.return_expression(b, cobj, T)

            h = (
                pyunits.convert(
                    quad(integrand, Tr, T) + Tr,
                    units.ENERGY_MOLE,
                )
            )

            return h

    class entr_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_liq_comp_coeff_A"):
                ChemSep.cp_mol_liq_comp.build_parameters(cobj)

            units = cobj.parent_block().get_metadata().derived_units

            cobj.entr_mol_form_liq_comp_ref = Var(
                doc="Liquid phase molar entropy of formation @ Tref",
                units=units.ENTROPY_MOLE,
            )
            set_param_from_config(cobj, param="entr_mol_form_liq_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy
            T = pyunits.convert(T, to_units=pyunits.K)
            Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

            units = b.params.get_metadata().derived_units

            integrand = lambda T: ChemSep.cp_mol_liq_comp.return_expression(b, cobj, T)

            s = (
                pyunits.convert(
                    quad(integrand, Tr, T)/T + Tr,
                    units.ENERGY_MOLE,
                )
            )

            return s

    class dens_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            ChemSep._create_params(cobj, "dens_mol_liq_comp_coeff")
            ChemSep._set_params(cobj, "dens_mol_liq_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Molar density
            Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

            rho = (
                cobj.dens_mol_liq_comp_coeff_A * (1 - Tr) ^
                (
                    cobj.dens_mol_liq_comp_coeff_B
                    + cobj.dens_mol_liq_comp_coeff_C * Tr
                    + cobj.dens_mol_liq_comp_coeff_D * Tr**2
                    + cobj.dens_mol_liq_comp_coeff_E * Tr**3
                )
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(rho, units.MOLAR_DENSITY)
    
    #----------------------------
    # Internal Methods

    @staticmethod
    def _create_params(cobj, param):
        for c in "ABCDE":
            cobj.add_component(
                param + "_" + c,
                Var(
                    doc=f"{param} parameter {c}",
                    units=None,
                ),
            )
    
    @staticmethod
    def _set_params(cobj, param):
        for c in "ABCDE":
            set_param_from_config(cobj, param=param, index=c)

    @staticmethod
    def _get_property(comp_name, prop_name, index=None):
        print(comp_name)
        if comp_name not in ChemSep._cached_components:
            ChemSep._cached_components[comp_name] = get_compound(comp_name)
        if index is not None:
            return ChemSep._cached_components[comp_name][prop_name][index]
        else:
            return ChemSep._cached_components[comp_name][prop_name]
