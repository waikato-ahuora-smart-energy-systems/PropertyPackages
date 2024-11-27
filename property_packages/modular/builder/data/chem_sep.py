"""
Methods for calculating pure component properties from:

http://www.chemsep.org/downloads/docs/book.htm

IDAES naming conventions followed for compatibility with modular property packages
"""

from idaes.core.util.misc import set_param_from_config
from pyomo.environ import log, exp, units as pyunits
from scipy.integrate import quad
from pyomo.environ import Var, value, units
from abc import abstractmethod
from typing import List, Literal

class ChemSep(object):

    class cp_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            ChemSepEqn.build_parameters(
                prefix="cp_mol_ig_comp", 
                doc="ideal gas molar heat capacity", 
                cobj=cobj)

        @staticmethod
        def return_expression(b, cobj, T):
            return ChemSepEqn.return_expression(
                prefix="cp_mol_ig_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="regular", 
                end_units="HEAT_CAPACITY_MOLE"
            )

    class enth_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_ig_comp_coeff_A"):
                ChemSep.cp_mol_ig_comp.build_parameters(cobj)
            
            units = cobj.parent_block().get_metadata().derived_units

            cobj.enth_mol_form_vap_comp_ref = Var(
                doc="Vapor phase molar enthalpy of formation",
                units=units.ENERGY_MOLE,
            )

            set_param_from_config(cobj, param="enth_mol_form_vap_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            return ChemSepEqn.return_expression(
                prefix="cp_mol_ig_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="enth", 
                end_units="ENERGY_MOLE"
            )

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
            return ChemSepEqn.return_expression(
                prefix="cp_mol_ig_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="entr", 
                end_units="ENTROPY_MOLE"
            )
            
    class pressure_sat_comp:
        @staticmethod
        def build_parameters(cobj):
            ChemSepEqn.build_parameters(
                prefix="pressure_sat_comp", 
                doc="calculating P-sat", 
                cobj=cobj)

        @staticmethod
        def return_expression(b, cobj, T, dT=False):
            return ChemSepEqn.return_expression(
                prefix="pressure_sat_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="regular", 
                end_units="PRESSURE"
            )
    
    class cp_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            ChemSepEqn.build_parameters(
                prefix="cp_mol_liq_comp", 
                doc="liquid phase molar heat capacity", 
                cobj=cobj)

        @staticmethod
        def return_expression(b, cobj, T):
            return ChemSepEqn.return_expression(
                prefix="cp_mol_liq_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="regular", 
                end_units="HEAT_CAPACITY_MOLE"
            )

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
            return ChemSepEqn.return_expression(
                prefix="cp_mol_liq_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="enth", 
                end_units="ENERGY_MOLE"
            )

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
            return ChemSepEqn.return_expression(
                prefix="cp_mol_liq_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="entr", 
                end_units="ENTROPY_MOLE"
            )

    class dens_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            ChemSepEqn.build_parameters(
                prefix="dens_mol_liq_comp", 
                doc="liquid phase molar density", 
                cobj=cobj)

        @staticmethod
        def return_expression(b, cobj, T):
            return ChemSepEqn.return_expression(
                prefix="dens_mol_liq_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="regular", 
                end_units="DENSITY_MOLE"
            )

EqnType = Literal["entr"] | Literal["enth"] | Literal["regular"]

class ChemSepEqn:

    @staticmethod
    def return_expression(prefix, b, cobj, T, eqn_type, end_units):
        """
        Generic return expression method which calls the appropriate equation method

        Args:
            prefix : str
                Prefix of the associated parameters
            b : Block
                Block object
            cobj : Component
                Component object
            T : float
                Temperature equation is evaluated at
            eqn_str : str
                Equation string
            units : Any
                Units
        
        Returns:
            float : Result of the equation
        """
        eqn_number = str(ChemSepEqn.get_params(cobj, prefix)[5].value)
        if eqn_number is None:
            raise Exception(f"Equation number not found for {prefix}")
        units = b.params.get_metadata().derived_units.__getitem__(end_units)
        eqn_obj = None

        if eqn_number in equation_map:
            eqn_obj = equation_map[eqn_number]
        else:
            raise Exception(f"Equation {eqn_number} not found")
            
        res = None

        if eqn_type == "enth":
            res = eqn_obj.enth(prefix, b, cobj, T)
        elif eqn_type == "entr":
            res = eqn_obj.entr(prefix, b, cobj, T)
        else:
            res = eqn_obj.return_expression(prefix, b, cobj, T)
        
        print(res)
        return res * units

    @staticmethod
    def get_params(cobj, prefix):
        coeff_names = ['A', 'B', 'C', 'D', 'E', "eqno", "units"]
        params = {name: getattr(cobj, f"{prefix}_coeff_{name}", None) for name in coeff_names}
        return tuple(params[name] for name in coeff_names)
    
    @staticmethod
    def build_parameters(prefix, doc, cobj):
        coeff_names = ['A', 'B', 'C', 'D', 'E', "eqno", "units"]
        for name in coeff_names:
            setattr(cobj, f"{prefix}_coeff_{name}", Var(
                doc=f"Parameter {name} for {doc}",
                units=None,
            ))
            set_param_from_config(cobj, param=f"{prefix}_coeff", index=name)

class eqn_100:

    @staticmethod
    def return_expression(prefix, b, cobj, T):
        # Ensuring temperature is in Kelvin
        T = pyunits.convert(T, to_units=pyunits.K)
        # Retrieving A-E coefficients based on prefix
        A, B, C, D, E, eqno, units = ChemSepEqn.get_params(cobj, prefix)
        # Equation 100 taken from Chem Sep Book
        rho = (A + B * T + C * T^2 + D * T^3 + E * T^4) * units
        # Getting the final units
        return rho
    
    def enth(prefix, b, cobj, T):
        # Specific enthalpy (eq_100)
        T = pyunits.convert(T, to_units=pyunits.K)
        Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

        A, B, C, D, E, eqno, units = ChemSepEqn.get_params(cobj, prefix)

        h_form = (
            cobj.enth_mol_form_vap_comp_ref
            if b.params.config.include_enthalpy_of_formation
            else 0
        )

        h = (A * (T - Tr)
            + (B / 2) * (T**2 - Tr**2)
            + (C / 3) * (T**3 - Tr**3)
            + (D / 4) * (T**4 - Tr**4)
            + (E / 5) * (T**5 - Tr**5)
        ) + h_form

        return h

    def entr(prefix, b, cobj, T):
        # Specific entropy (eq_100)
        T = pyunits.convert(T, to_units=pyunits.K)
        Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

        A, B, C, D, E, eqno, units = ChemSepEqn.get_params(cobj, prefix)

        s = (A * log(T / Tr)
            + B * (T - Tr)
            + (C / 2) * (T ** 2 - Tr ** 2)
            + (D / 3) * (T ** 3 - Tr ** 3)
            + (E / 4) * (T ** 4 - Tr ** 4)
        ) + cobj.entr_mol_form_vap_comp_ref

        return s

class eqn_106:

    @staticmethod
    def return_expression(prefix, b, cobj, T):
        # Ensuring temperature is in Kelvin
        T = pyunits.convert(T, to_units=pyunits.K)
        # Retrieving A-E coefficients based on prefix
        A, B, C, D, E, eqno, units = ChemSepEqn.get_params(cobj, prefix)
        # Equation 106 taken from Chem Sep Book
        rho = params.A + params.B * T + params.C * T^2 + params.D * T^3 + params.E * T^4
        units = b.params.get_metadata().derived_units
        return pyunits.convert(rho, units.DENSITY_MOLE)

class eqn_4:

    @staticmethod
    def build_parameters(prefix, doc, cobj):
        pass

class eqn_16:

    @staticmethod
    def return_expression(prefix, b, cobj, T):
        # Ensuring temperature is in Kelvin
        T = pyunits.convert(T, to_units=pyunits.K)
        # Retrieving A-E coefficients based on prefix
        A, B, C, D, E, eqno, units = ChemSepEqn.get_params(cobj, prefix)
        # Equation 16 taken from Chem Sep Book
        eqn = (A + exp( B / T + C + D * T + E * T**2))
        return eqn
    
    def enth(prefix, b, cobj, T):
        # Specific enthalpy (eq_16)
        T = pyunits.convert(T, to_units=pyunits.K)
        Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)
        def integrand(T):
            return eqn_16.return_expression(prefix, b, cobj, T)
        return pyunits.convert(quad(integrand, Tr, T) + Tr, units.ENERGY_MOLE)

    def entr(prefix, b, cobj, T):
        # Specific entropy (eq_16)
        T = pyunits.convert(T, to_units=pyunits.K)
        Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)
        def integrand(T):
            return eqn_16.return_expression(prefix, b, cobj, T)
        return pyunits.convert(quad(integrand, Tr, T) / T + Tr, units.ENTROPY_MOLE)

class eqn_10:

    @staticmethod
    def return_expression(prefix, b, cobj, T):
        # Ensuring temperature is in Kelvin
        T = value(pyunits.convert(T, to_units=pyunits.K))
        # Retrieving A-E coefficients based on prefix
        A, B, C, D, E, eqno, units = ChemSepEqn.get_params(cobj, prefix)
        # Calculating final result
        eqn = (exp(A - B / (T + C)))
        return eqn

class eqn_105:

    @staticmethod
    def build_parameters(prefix, doc, cobj):
        pass

equation_map = {
    "100": eqn_100,
    "106": eqn_106,
    "4": eqn_4,
    "16": eqn_16,
    "10": eqn_10,
    "105": eqn_105,
}