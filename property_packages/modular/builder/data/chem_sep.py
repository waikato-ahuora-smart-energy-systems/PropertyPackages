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
            raise Exception("Does not need")
            return ChemSepEqn.return_expression(
                prefix="cp_mol_ig_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="heatcp", 
                units="HEAT_CAPACITY_MOLE"
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
                units="ENERGY_MOLE"
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
                units="ENTROPY_MOLE"
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
                eqn_type="heatcp", 
                units="PRESSURE"
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
            raise Exception("Does not need")
            return ChemSepEqn.return_expression(
                prefix="cp_mol_liq_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="heatcp", 
                units="HEAT_CAPACITY_MOLE"
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
            raise Exception("Does not need")
            return ChemSepEqn.return_expression(
                prefix="cp_mol_liq_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="enth", 
                units="ENERGY_MOLE"
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
            raise Exception("Does not need")
            return ChemSepEqn.return_expression(
                prefix="cp_mol_liq_comp",
                b=b, cobj=cobj, T=T, 
                eqn_type="entr", 
                units="ENTROPY_MOLE"
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
                eqn_type="heatcp", 
                units="DENSITY_MOLE"
            )

EqnType = Literal["entr"] | Literal["enth"] | Literal["heatcp"]
UnitsType = Literal["HEAT_CAPCITY_MOLE"] | Literal["ENERGY_MOLE"] | Literal["ENTROPY_MOLE"] | Literal["DENSITY_MOLE"]
PrefixType = Literal["cp_mol_ig_comp"] | ["enth_mol_ig_comp"] | ["entr_mol_ig_comp"] | ["pressure_sat_comp"] | ["cp_mol_liq_comp"] | ["enth_mol_liq_comp"] | ["entr_mol_liq_comp"] | ["dens_mol_liq_comp"]

class ChemSepEqn:

    @staticmethod
    def get_unit_from_string(unit_str):

        unit_dict = {
            "K": pyunits.K,
            "Pa": pyunits.Pa,
            "m3/kmol": pyunits.m**3 / pyunits.kmol,
            "kmol/m3": pyunits.kmol / pyunits.m**3,
            "kg/kmol": pyunits.kg / pyunits.kmol,
            "J/kmol": pyunits.J / pyunits.kmol,
            "J/kmol/K": pyunits.J / (pyunits.kmol * pyunits.K),
            "m": pyunits.m,
            "_": pyunits.dimensionless,
            "J0.5/m1.5": pyunits.J**0.5 / pyunits.m**1.5,
            "Coulomb.m": pyunits.C * pyunits.m, # TODO: update
            "m2/kmol": pyunits.m**2 / pyunits.kmol,
            "W/m/K": pyunits.W / (pyunits.m * pyunits.K),
            "N/m": pyunits.N / pyunits.m,
            "kg0.25.m3/s0.5/kmol": pyunits.kg**0.25 * pyunits.m**3 / (pyunits.s**0.5 * pyunits.kmol),
            "m3/kmol": pyunits.m**3 / pyunits.kmol,
            "Pa.s": pyunits.Pa * pyunits.s,
        }

        try:
            unit_str = value(unit_str)
        except:
            unit_str = unit_str

        if unit_str in unit_dict:
            return unit_dict[unit_str]
        else:
            raise Exception(f"Unit {unit_str} not found")
            return None

    @staticmethod
    def get_units_from_unittype(b, units: UnitsType):
        return b.params.get_metadata().derived_units.__getitem__(units)

    @staticmethod
    def return_expression(b, cobj, T, prefix: PrefixType, eqn_type: EqnType, units: UnitsType):
        """
        Generic return expression method which calls the appropriate equation method

        Args:
            prefix : PrefixType
                Prefix of the associated coefficients

            eqn_str : EqnType
                Identifier for the equation type
            
            units : UnitsType
                Final units for the expression
    
        Returns:
            Pyomo expression
        """

        A, B, C, D, E, eqno, _ = ChemSepEqn.get_params(cobj, prefix)
        
        if value(eqno) is None:
            raise Exception(f"Equation number not found for {prefix}")
        else:
            if str(value(eqno)) in equation_map:
                eqn_obj = equation_map[str(value(eqno))]
                if eqn_type == "enth":
                    return eqn_obj.enth(b, cobj, T, prefix, units)
                elif eqn_type == "entr":
                    return eqn_obj.entr(b, cobj, T, prefix, units)
                elif eqn_type == "heatcp":
                    return eqn_obj.return_expression(b, cobj, T, prefix, units)
                else:
                    raise Exception(f"Invalid equation type: {eqn_type}")
            else:
                raise Exception(f"Equation {eqno} not found")

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

"""

Equations taken from ChemSep book

http://www.chemsep.org/book/docs/book2.pdf

All implementations except should evaluate the equations as according to the book
using dimensionless coefficients A-E and temperature T. This value should then be
converted to the appropriate final units as defined within the database. Before, 
finally being converted to the desired units. Derived from the block using the
identifier passed in.

"""


class eqn_100:

    @staticmethod
    def return_expression(b, cobj, T, prefix, units):
        # Ensuring temperature is in Kelvin
        T = value(pyunits.convert(T, to_units=pyunits.K))

        # Retrieving A-E coefficients based on prefix
        A, B, C, D, E, eqno, eq_units = ChemSepEqn.get_params(cobj, prefix)

        # Equation 100 taken from Chem Sep Book
        eq_units = ChemSepEqn.get_unit_from_string(eq_units)
        converted_units = ChemSepEqn.get_units_from_unittype(b, units)
        rho = (A + B * T + C * T^2 + D * T^3 + E * T^4) * eq_units

        # Getting the final units
        return pyunits.convert(rho, converted_units)
    
    def enth(b, cobj, T, prefix, units):
        # Specific enthalpy (eq_100)
        T = value(pyunits.convert(T, to_units=pyunits.K))
        Tr = value(pyunits.convert(b.params.temperature_ref, to_units=pyunits.K))

        A, B, C, D, E, eqno, eq_units = ChemSepEqn.get_params(cobj, prefix)

        eq_units = ChemSepEqn.get_unit_from_string(eq_units)
        converted_units = ChemSepEqn.get_units_from_unittype(b, units)

        h_form = (
            cobj.enth_mol_form_vap_comp_ref
            if b.params.config.include_enthalpy_of_formation
            else 0 * converted_units
        )

        h = (A * (T - Tr)
            + (B / 2) * (T**2 - Tr**2)
            + (C / 3) * (T**3 - Tr**3)
            + (D / 4) * (T**4 - Tr**4)
            + (E / 5) * (T**5 - Tr**5)) * converted_units + h_form

        return h

    def entr(b, cobj, T, prefix, units):
        # Specific entropy (eq_100)
        T = value(pyunits.convert(T, to_units=pyunits.K))
        Tr = value(pyunits.convert(b.params.temperature_ref, to_units=pyunits.K))

        A, B, C, D, E, eqno, eq_units = ChemSepEqn.get_params(cobj, prefix)

        eq_units = ChemSepEqn.get_unit_from_string(eq_units)
        converted_units = ChemSepEqn.get_units_from_unittype(b, units)

        s = pyunits.convert((A * log(T / Tr)
            + B * (T - Tr)
            + (C / 2) * (T ** 2 - Tr ** 2)
            + (D / 3) * (T ** 3 - Tr ** 3)
            + (E / 4) * (T ** 4 - Tr ** 4)
        ) * eq_units, converted_units) + cobj.entr_mol_form_vap_comp_ref

        return s

class eqn_106:

    @staticmethod
    def return_expression(b, cobj, T, prefix, units):
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
    def return_expression(b, cobj, T, prefix, units):
        # Ensuring temperature is in Kelvin
        T = value(pyunits.convert(T, to_units=pyunits.K))

        # Retrieving A-E coefficients based on prefix
        A, B, C, D, E, eqno, eq_units = ChemSepEqn.get_params(cobj, prefix)
        eq_units = ChemSepEqn.get_unit_from_string(eq_units)
        converted_units = ChemSepEqn.get_units_from_unittype(b, units)

        # Equation 4 taken from Chem Sep Book
        eqn = (A + B * T + C * T**2 + D * T**3) * eq_units
        return pyunits.convert(eqn, converted_units)
    
    @staticmethod
    def entr(b, cobj, T, prefix, units):
        # Specific entropy (eq_4)
        T = value(pyunits.convert(T, to_units=pyunits.K))
        Tr = value(pyunits.convert(b.params.temperature_ref, to_units=pyunits.K))

        A, B, C, D, E, eqno, eq_units = ChemSepEqn.get_params(cobj, prefix)

        eq_units = ChemSepEqn.get_unit_from_string(eq_units)
        converted_units = ChemSepEqn.get_units_from_unittype(b, units)

        s = pyunits.convert(((D / 3) * (T**3 - Tr**3)
            + (C / 2) * (T**2 - Tr**2)
            + B * (T - Tr)
            + A * log(T / Tr)) * eq_units, converted_units) + cobj.entr_mol_form_vap_comp_ref

        return s
    
    @staticmethod
    def enth(b, cobj, T, prefix, units):
        # Specific enthalpy (eq_4)
        T = value(pyunits.convert(T, to_units=pyunits.K))
        Tr = value(pyunits.convert(b.params.temperature_ref, to_units=pyunits.K))

        A, B, C, D, E, eqno, eq_units = ChemSepEqn.get_params(cobj, prefix)

        eq_units = ChemSepEqn.get_unit_from_string(eq_units)
        converted_units = ChemSepEqn.get_units_from_unittype(b, units)

        aunits = b.params.get_metadata().derived_units

        raise Exception(aunits.ENERGY_MOLE)

        #raise Exception(eq_units)

        # h_form = (
        #     cobj.enth_mol_form_vap_comp_ref
        #     if b.params.config.include_enthalpy_of_formation
        #     else 0 * converted_units
        # )

        h = pyunits.convert((
                (D / 4) * (T**4 - Tr**4)
                + (C / 3) * (T**3 - Tr**3)
                + (B / 2) * (T**2 - Tr**2)
                + A * (T - Tr)) * eq_units, converted_units) + cobj.enth_mol_form_vap_comp_ref

        return h

class eqn_16:

    @staticmethod
    def return_expression(b, cobj, T, prefix, units):
        # Ensuring temperature is in Kelvin
        if pyunits.get_units(T) != pyunits.dimensionless: # This check is required for the enthalpy calculation
            T = value(pyunits.convert(T, to_units=pyunits.K))
        
        # Retrieving A-E coefficients based on prefix
        A, B, C, D, E, eqno, eq_units = ChemSepEqn.get_params(cobj, prefix)

        # Converting string to pyomo units
        eq_units = ChemSepEqn.get_unit_from_string(eq_units)
        converted_units = ChemSepEqn.get_units_from_unittype(b, units)

        # Equation 16 taken from Chem Sep Book
        eqn = (A + exp( B / T + C + D * T + E * T**2)) * eq_units

        if(prefix == "cp_mol_liq_comp"):
            # Equation 16 taken from Chem Sep Book
            raise Exception(eqn)

        return pyunits.convert(eqn, converted_units)

    @staticmethod
    def enth(b, cobj, T, prefix, units):
        # Specific enthalpy (eq_16)
        T = value(pyunits.convert(T, to_units=pyunits.K))
        Tr = value(pyunits.convert(b.params.temperature_ref, to_units=pyunits.K))
        integrand = lambda T: value(eqn_16.return_expression(b, cobj, T, prefix, "HEAT_CAPACITY_MOLE"))
        h = quad(integrand, Tr, T)[0] + Tr
        converted_units = ChemSepEqn.get_units_from_unittype(b, units)
        return h * converted_units

    @staticmethod
    def entr(b, cobj, T, prefix, units):
        # Specific entropy (eq_16)
        T = value(pyunits.convert(T, to_units=pyunits.K))
        Tr = value(pyunits.convert(b.params.temperature_ref, to_units=pyunits.K))
        def integrand(T):
            return value(eqn_16.return_expression(b, cobj, T, prefix, "HEAT_CAPACITY_MOLE"))
        s = quad(integrand, Tr, T)[0] / T + Tr
        converted_units = ChemSepEqn.get_units_from_unittype(b, units)
        return s * converted_units

class eqn_10:

    @staticmethod
    def return_expression(b, cobj, T, prefix, units):
        # Ensuring temperature is in Kelvin
        T = value(pyunits.convert(T, to_units=pyunits.K))
        # Retrieving A-E coefficients based on prefix
        A, B, C, D, E, eqno, units = ChemSepEqn.get_params(cobj, prefix)
        # Calculating final result
        eqn = (exp(A - B / (T + C)))
        return eqn

equation_map = {
    "100": eqn_100,
    "10": eqn_10,
    "16": eqn_16,
    "106": eqn_106,
    "4": eqn_4
}