from math import floor
from typing import Any, Dict, List
from .base_parser import BuildBase
from compounds.Compound import Compound
from pyomo.environ import units as pyunits
from idaes.models.properties.modular_properties.state_definitions import FTPx, FPhx
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (LogBubbleDew)
from idaes.core import LiquidPhase, VaporPhase, Component, PhaseType as PT
from idaes.models.properties.modular_properties.phase_equil import (SmoothVLE)
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.pure import RPP4, RPP3, Perrys
from property_packages.modular.builder.data.chem_sep import ChemSep
from pyomo.common.fileutils import this_file_dir
from property_packages.types import States
from idaes.models.properties.modular_properties.eos.ideal import Ideal
import csv

class components_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
    
        def serialise_component(compound: Compound) -> Dict[str, Any]:

            # configuration default to all components 
            config = {
                "type": Component,
                "valid_phase_types": VaporPhase,
                "parameter_data": {
                    "mw": (compound["MolecularWeight"].value, pyunits.kg/pyunits.kilomol),
                    "pressure_crit": (compound["CriticalPressure"].value, pyunits.Pa),
                    "temperature_crit": (compound["CriticalTemperature"].value, pyunits.K),
                }
            }

            # Energies of Formation
            if compound["NISTVaporReference"] is not None: # this does not work when passed
                config["parameter_data"].update({
                    "enth_mol_form_vap_comp_ref": (compound["NISTVaporReference"].value, pyunits.J/pyunits.kilomol)
                })
            else:
                raise ValueError("No Heat of Formation Data")

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
                elif compound["RPPHeatCapacityCp"]["eqno"] == 100 or compound["RPPHeatCapacityCp"]["eqno"] == 5:
                    config["enth_mol_ig_comp"] = ChemSep
                    config["entr_mol_ig_comp"] = ChemSep
                    config["parameter_data"].update({"cp_mol_ig_comp_coeff": {
                        "A": (compound["RPPHeatCapacityCp"]["A"], pyunits.J / pyunits.kilomol / pyunits.K),
                        "B": (compound["RPPHeatCapacityCp"]["B"], pyunits.J / pyunits.kilomol / pyunits.K**2),
                        "C": (compound["RPPHeatCapacityCp"]["C"], pyunits.J / pyunits.kilomol / pyunits.K**3),
                        "D": (compound["RPPHeatCapacityCp"]["D"], pyunits.J / pyunits.kilomol / pyunits.K**4),
                        "E": (compound["RPPHeatCapacityCp"]["E"], pyunits.J / pyunits.kilomol / pyunits.K**5),
                    }})
                else:
                    raise ValueError(f"Invalid equation number for heat capacity {compound['RPPHeatCapacityCp']['eqno']}")
            else:
                raise ValueError("No Heat Capacity Data")

            # Saturation Pressure (Vapor)
            if compound["NISTVaporPressure"] is not None:
                config["pressure_sat_comp"] = ChemSep
                config["parameter_data"].update({"pressure_sat_comp_coeff": {
                    "A": (compound["NISTVaporPressure"]["A"], None),
                    "B": (compound["NISTVaporPressure"]["B"], pyunits.K),
                    "C": (compound["NISTVaporPressure"]["C"], pyunits.K),
                }})
            else:
                raise ValueError("No NIST Vapor Pressure Data")
            
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
                elif compound["LiquidHeatCapacityCp"]["eqno"] == 16:
                    # Uses correct equations to calculate
                    config["enth_mol_liq_comp"] = ChemSep
                    config["entr_mol_liq_comp"] = ChemSep
                    config["parameter_data"].update({"cp_mol_liq_comp_coeff": {
                        "A": (compound["LiquidHeatCapacityCp"]["A"], pyunits.J / pyunits.kilomol / pyunits.K**1),
                        "B": (compound["LiquidHeatCapacityCp"]["B"], pyunits.J / pyunits.kilomol / pyunits.K**2),
                        "C": (compound["LiquidHeatCapacityCp"]["C"], pyunits.J / pyunits.kilomol / pyunits.K**3),
                        "D": (compound["LiquidHeatCapacityCp"]["D"], pyunits.J / pyunits.kilomol / pyunits.K**4),
                        "E": (compound["LiquidHeatCapacityCp"]["E"], pyunits.J / pyunits.kilomol / pyunits.K**5),
                    }})

                    # ASSUMPTION: Molar heat of formation, liq is zero - given the semi-okay by Ben
                    config["parameter_data"].update({
                        "enth_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.kilomol)
                    })

                    config["parameter_data"].update({
                        "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.kilomol / pyunits.K)
                    }) 
                else:
                    raise ValueError(f"No Liquid Heat Capacity Equation Data {compound['LiquidHeatCapacityCp']['eqno']}")
            else: 
                # Compound only exists in vapor phase
                config["parameter_data"].update({"valid_phase_types": PT.vaporphase})

            return config
        
        def valid_phases(compound: Compound) -> PT:
            return PT.vaporPhase
        
        components_output = {}
        for compound in compounds:
            components_output[compound["CompoundID"].value] = serialise_component(compound)
        return components_output

class phases_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
      return {"Vap": {"type": VaporPhase, "equation_of_state": Ideal}},

class phases_in_equilibrium_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        return [("Vap", "Liq")]

class state_bounds_parser(BuildBase):
    """
    State bounds are used to find optimal solution
    """

    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        return {
            "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
            "temperature": (273.15, 300, 1600, pyunits.K),
            "pressure": (5e3, 1e5, 1e6, pyunits.Pa),
        }

class temperature_ref_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        return (300, pyunits.K)
    
class pressure_ref_parser(BuildBase):
    @staticmethod
    def serialise(compounds: List[Compound], valid_states: List[States]) -> Dict[str, Any]:
        return (1e5, pyunits.Pa)
