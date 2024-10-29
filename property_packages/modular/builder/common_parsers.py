from typing import Any, Dict, List
from .base_parser import BuildBase
from compounds.Compound import Compound
from pyomo.environ import units as pyunits
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (LogBubbleDew)
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.models.properties.modular_properties.phase_equil import (SmoothVLE)
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.pure import RPP4, RPP3, Perrys
from property_packages.modular.builder.data.chem_sep import ChemSep

class base_units_parser(BuildBase):

    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return {
            'time': pyunits.s,
            'length': pyunits.m,
            'mass': pyunits.kg,
            'amount': pyunits.mol,
            'temperature': pyunits.K,
        }

class bubble_dew_method_parser(BuildBase):
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return LogBubbleDew

class components_parser(BuildBase):

    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
    
        def serialise_component(compound: Compound) -> Dict[str, Any]:

            # configuration default to all components 
            config = {
                "type": Component,
                "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                "parameter_data": {
                    "mw": (compound["MolecularWeight"].value, pyunits.kg/pyunits.kilomol),
                    "pressure_crit": (compound["CriticalPressure"].value, pyunits.Pa),
                    "temperature_crit": (compound["CriticalTemperature"].value, pyunits.K),
                    "omega": compound["AcentricityFactor"].value,
                }
            }

            # Energies of Formation
            if compound["HeatOfFormation"] is not None: # this does not work when passed
                config["parameter_data"].update({
                    "enth_mol_form_vap_comp_ref": (compound["HeatOfFormation"].value, pyunits.J/pyunits.kilomol)
                })
            if compound["AbsEntropy"] is not None:
                config["parameter_data"].update({
                    "entr_mol_form_vap_comp_ref": (-1 * compound["AbsEntropy"].value, pyunits.J/pyunits.kilomol/pyunits.K)
                })

            # Ideal Gas Molar Calculations
            # All three properties intrinsically linked together
            if compound["RPPHeatCapacityCp"] is not None:
                if compound["RPPHeatCapacityCp"]["eqno"] == 4:
                    config["enth_mol_ig_comp"] = RPP4
                    config["entr_mol_ig_comp"] = RPP4
                    config["parameter_data"].update({"cp_mol_ig_comp_coeff": {
                        "A": (float(compound["RPPHeatCapacityCp"]["A"]), pyunits.J / pyunits.kilomol / pyunits.K),
                        "B": (float(compound["RPPHeatCapacityCp"]["B"]), pyunits.J / pyunits.kilomol / pyunits.K**2),
                        "C": (float(compound["RPPHeatCapacityCp"]["C"]), pyunits.J / pyunits.kilomol / pyunits.K**3),
                        "D": (float(compound["RPPHeatCapacityCp"]["D"]), pyunits.J / pyunits.kilomol / pyunits.K**4),
                    }})
                elif compound["RPPHeatCapacityCp"]["eqno"] == 100:
                    config["enth_mol_ig_comp"] = ChemSep
                    config["entr_mol_ig_comp"] = ChemSep
                    config["parameter_data"].update({"cp_mol_ig_comp_coeff": {
                        "A": (compound["RPPHeatCapacityCp"]["A"], pyunits.J / pyunits.kilomol / pyunits.K),
                        "B": (compound["RPPHeatCapacityCp"]["B"], pyunits.J / pyunits.kilomol / pyunits.K**2),
                        "C": (compound["RPPHeatCapacityCp"]["C"], pyunits.J / pyunits.kilomol / pyunits.K**3),
                        "D": (compound["RPPHeatCapacityCp"]["D"], pyunits.J / pyunits.kilomol / pyunits.K**4),
                        "E": (compound["RPPHeatCapacityCp"]["E"], pyunits.J / pyunits.kilomol / pyunits.K**5),
                    }})

            # Saturation Pressure (Vapor)
            if compound["VaporPressure"] is not None:
                config["pressure_sat_comp"] = ChemSep
                config["parameter_data"].update({"pressure_sat_comp_coeff": {
                    "A": (compound["VaporPressure"]["A"], None),
                    "B": (compound["VaporPressure"]["B"], pyunits.K),
                    "C": (compound["VaporPressure"]["C"], pyunits.K),
                }})

            # Liquid Density
            if compound["LiquidDensity"] is not None:
                if compound["LiquidDensity"]["eqno"] == 105:
                    config["dens_mol_liq_comp"] = Perrys
                    config["parameter_data"].update({"dens_mol_liq_comp_coeff": {
                        "eqn_type": 1,
                        "1": (compound["LiquidDensity"]["A"], pyunits.kmol / pyunits.m**3),
                        "2": (compound["LiquidDensity"]["B"], None),
                        "3": (compound["LiquidDensity"]["C"], pyunits.K),
                        "4": (compound["LiquidDensity"]["D"], None),
                    }})

            # Liquid Heat Capacity & Entropy / Enthalpy
            if compound["LiquidHeatCapacityCp"] is not None:
                if compound["LiquidHeatCapacityCp"]["eqno"] == 100:
                    # Uses correct equations to calculate
                    config["enth_mol_liq_comp"] = Perrys
                    config["entr_mol_liq_comp"] = Perrys
                    config["parameter_data"].update({"cp_mol_liq_comp_coeff": {
                        "1": (compound["LiquidHeatCapacityCp"]["A"], pyunits.J / pyunits.kilomol / pyunits.K**1),
                        "2": (compound["LiquidHeatCapacityCp"]["B"], pyunits.J / pyunits.kilomol / pyunits.K**2),
                        "3": (compound["LiquidHeatCapacityCp"]["C"], pyunits.J / pyunits.kilomol / pyunits.K**3),
                        "4": (compound["LiquidHeatCapacityCp"]["D"], pyunits.J / pyunits.kilomol / pyunits.K**4),
                        "5": (compound["LiquidHeatCapacityCp"]["E"], pyunits.J / pyunits.kilomol / pyunits.K**5),
                    }})
                    
                    # ASSUMPTION: Molar heat of formation, liq is zero - given the semi-okay by Ben
                    config["parameter_data"].update({
                        "enth_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.kilomol)
                    })

                    config["parameter_data"].update({
                        "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.kilomol / pyunits.K)
                    })  

            return config
        
        components_output = {}
        for compound in compounds:
            components_output[compound["CompoundID"].value] = serialise_component(compound)
        return components_output

class phase_equilibrium_state_parser(BuildBase):
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return {("Vap", "Liq"): SmoothVLE} # TODO: Update to Smooth VLE_v2 https://github.com/IDAES/idaes-pse/blob/main/idaes/models/properties/modular_properties/phase_equil/smooth_VLE_2.py

class phases_parser(BuildBase):
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return {
            "Liq": {
                "type": LiquidPhase,
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR},
            },
            "Vap": {
                "type": VaporPhase,
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR},
            },
        }

class phases_in_equilibrium_parser(BuildBase):
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return [("Vap", "Liq")]

class pressure_ref_parser(BuildBase):
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return (101325, pyunits.Pa)

class state_bounds_parser(BuildBase):
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return {
            "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
            "temperature": (273.15, 300, 500, pyunits.K),
            "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
        }

class state_definition_parser(BuildBase):
    def serialise(compounds: List[Compound]) -> str:
        return FTPx

class temperature_ref_parser(BuildBase):
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        return (298.15, pyunits.K)


class pr_kappa_parser(BuildBase):
    def serialise(compounds: List[Compound]) -> Dict[str, Any]:
        kappa_parameters = {}
        for i, compound1 in enumerate(compounds):
            for j, compound2 in enumerate(compounds):
                kappa_parameters[(compound1["CompoundID"].value, compound2["CompoundID"].value)] = 0.000
                # Setting all interactions initially to zero
                # TODO: Pass over the correct values
        return {"PR_kappa": kappa_parameters}
