"""
Methods for calculating pure component properties from:

http://www.chemsep.org/downloads/docs/book.htm

IDAES naming conventions followed for compatibility with modular property packages
"""

from idaes.core.util.misc import set_param_from_config
from pyomo.environ import log, exp, units as pyunits
from scipy.integrate import quad
from pyomo.environ import Var

class ChemSep(object):

    class cp_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.cp_mol_ig_comp_coeff_A = Var(
                doc="Parameter A for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.kilomol / pyunits.K,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="A")

            cobj.cp_mol_ig_comp_coeff_B = Var(
                doc="Parameter B for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.kilomol / pyunits.K ** 2,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="B")

            cobj.cp_mol_ig_comp_coeff_C = Var(
                doc="Parameter C for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.kilomol / pyunits.K ** 3,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="C")

            cobj.cp_mol_ig_comp_coeff_D = Var(
                doc="Parameter D for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.kilomol / pyunits.K ** 4,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="D")

            cobj.cp_mol_ig_comp_coeff_E = Var(
                doc="Parameter E for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.kilomol / pyunits.K ** 5,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="E")

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
    
    class enth_entr_mol_ig:
        @staticmethod
        def build_parameters(cobj):
            cobj.enth_entr_mol_ig_coeff_A = Var(
                doc="Parameter A for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.kilomol / pyunits.K,
            )
            set_param_from_config(cobj, param="enth_entr_mol_ig_coeff", index="A")

            cobj.enth_entr_mol_ig_coeff_B = Var(
                doc="Parameter B for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.kilomol / pyunits.K ** 2,
            )
            set_param_from_config(cobj, param="enth_entr_mol_ig_coeff", index="B")

            cobj.enth_entr_mol_ig_coeff_C = Var(
                doc="Parameter C for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.kilomol / pyunits.K ** 3,
            )
            set_param_from_config(cobj, param="enth_entr_mol_ig_coeff", index="C")

            cobj.enth_entr_mol_ig_coeff_D = Var(
                doc="Parameter D for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.kilomol / pyunits.K ** 4,
            )
            set_param_from_config(cobj, param="enth_entr_mol_ig_coeff", index="D")

            cobj.enth_entr_mol_ig_coeff_E = Var(
                doc="Parameter E for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.kilomol / pyunits.K ** 5,
            )
            set_param_from_config(cobj, param="enth_entr_mol_ig_coeff", index="E")
    
    class enth_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            # Required for calculating entropy
            if not hasattr(cobj, "enth_entr_mol_ig_coeff_A"):
                ChemSep.enth_entr_mol_ig.build_parameters(cobj)

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
                    (cobj.enth_entr_mol_ig_coeff_A * (T - Tr)
                     + (cobj.enth_entr_mol_ig_coeff_B / 2) * (T**2 - Tr**2)
                     + (cobj.enth_entr_mol_ig_coeff_C / 3) * (T**3 - Tr**3)
                     + (cobj.enth_entr_mol_ig_coeff_D / 4) * (T**4 - Tr**4)
                     + (cobj.enth_entr_mol_ig_coeff_E / 5) * (T**5 - Tr**5)),
                    units.ENERGY_MOLE,
                ) + h_form
            )

            return h

    class entr_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            # Required for calculating enthalpy
            if not hasattr(cobj, "enth_entr_mol_ig_coeff_A"):
                ChemSep.enth_entr_mol_ig.build_parameters(cobj)

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
                    (cobj.enth_entr_mol_ig_coeff_A * log(T / Tr)
                     + cobj.enth_entr_mol_ig_coeff_B * (T - Tr)
                     + (cobj.enth_entr_mol_ig_coeff_C / 2) * (T ** 2 - Tr ** 2)
                     + (cobj.enth_entr_mol_ig_coeff_D / 3) * (T ** 3 - Tr ** 3)
                     + (cobj.enth_entr_mol_ig_coeff_E / 4) * (T ** 4 - Tr ** 4)),
                    units.ENTROPY_MOLE,
                ) + cobj.entr_mol_form_vap_comp_ref
            )

            return s

    class pressure_sat_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.pressure_sat_comp_coeff_A = Var(
                doc="Antoine A coefficient for calculating P-sat",
                units=pyunits.dimensionless
            )
            set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="A")

            cobj.pressure_sat_comp_coeff_B = Var(
                doc="Antoine B coefficient for calculating P-sat",
                units=pyunits.K
            )
            set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="B")

            cobj.pressure_sat_comp_coeff_C = Var(
                doc="Antoine C coefficient for calculating P-sat",
                units=pyunits.K
            )
            set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="C")

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
            cobj.cp_mol_liq_comp_coeff_A = Var(
                doc="Parameter A for liquid phase molar heat capacity",
                units=pyunits.J * pyunits.kmol**-1 * pyunits.K**-1,
            )
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="A")

            cobj.cp_mol_liq_comp_coeff_B = Var(
                doc="Parameter B for liquid phase molar heat capacity",
                units=pyunits.J * pyunits.kmol**-1 * pyunits.K**-2,
            )
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="B")

            cobj.cp_mol_liq_comp_coeff_C = Var(
                doc="Parameter C for liquid phase molar heat capacity",
                units=pyunits.J * pyunits.kmol**-1 * pyunits.K**-3,
            )
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="C")

            cobj.cp_mol_liq_comp_coeff_D = Var(
                doc="Parameter D for liquid phase molar heat capacity",
                units=pyunits.J * pyunits.kmol**-1 * pyunits.K**-4,
            )
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="D")

            cobj.cp_mol_liq_comp_coeff_E = Var(
                doc="Parameter E for liquid phase molar heat capacity",
                units=pyunits.J * pyunits.kmol**-1 * pyunits.K**-5,
            )
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="E")

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
            cobj.dens_mol_liq_comp_coeff_A = Var(
                doc="Parameter A for liquid phase molar density",
                units=pyunits.kmol * pyunits.m**-3,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="A")

            cobj.dens_mol_liq_comp_coeff_B = Var(
                doc="Parameter B for liquid phase molar density",
                units=pyunits.dimensionless,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="B")

            cobj.dens_mol_liq_comp_coeff_C = Var(
                doc="Parameter C for liquid phase molar density",
                units=pyunits.dimensionless,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="C")

            cobj.dens_mol_liq_comp_coeff_D = Var(
                doc="Parameter D for liquid phase molar density",
                units=pyunits.dimensionless,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="D")

            cobj.dens_mol_liq_comp_coeff_E = Var(
                doc="Parameter E for liquid phase molar density",
                units=pyunits.dimensionless,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="E")

        @staticmethod
        def return_expression(b, cobj, T):

            T = pyunits.convert(T, to_units=pyunits.K)
            A = cobj.dens_mol_liq_comp_coeff_A
            B = cobj.dens_mol_liq_comp_coeff_B
            C = cobj.dens_mol_liq_comp_coeff_C
            D = cobj.dens_mol_liq_comp_coeff_D
            E = cobj.dens_mol_liq_comp_coeff_E

            rho = A * (1 - T)^(B + C * T + D * T^2 + E * T^3)

            units = b.params.get_metadata().derived_units

            return pyunits.convert(rho, units.DENSITY_MOLE)