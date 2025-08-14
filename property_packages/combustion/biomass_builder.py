from .biomass_combustion_rp import BMCombReactionParameterBlock
from .biomass_comb_pp import configuration
from property_packages.modular.modular_extended import GenericExtendedParameterBlock
from idaes.models.properties.modular_properties import GenericParameterBlock
from typing import List

def build_biomass_and_flue_package(compound_list: List[str]):

    # Validate compound list

    for compound in compound_list:
        if compound not in ["biomass", "water", "carbon dioxide", "oxygen", "carbon monoxide", "nitrogen", "ash"]:
            raise ValueError(f"Compound {compound} is not valid for biomass and flue package.")

    return GenericParameterBlock(**configuration)


def build_biomass_combustion_reaction_package(compound_list: List[str], pp=None):

    # Validate compound list

    for compound in compound_list:
        if compound not in ["biomass", "water", "carbon dioxide", "oxygen", "carbon monoxide", "nitrogen"]:
            raise ValueError(f"Compound {compound} is not valid for biomass combustion reaction package.")

    return BMCombReactionParameterBlock(
        property_package=pp
    )