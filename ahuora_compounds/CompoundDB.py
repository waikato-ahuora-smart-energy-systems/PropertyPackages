from ahuora_compounds import deprecated
from ahuora_compounds.CompoundRegistry import CompoundRegistry
from ahuora_compounds.Compound import Compound
from ahuora_compounds.RegistrySearch import RegistrySearch

__author__ = "Mahaki Leach"

# Main registry instance
_registry = CompoundRegistry()

# Frontend registry view
db = RegistrySearch(_registry)

_registry._discover_loaders()

@deprecated("get_compound")
def get_compound(name: str) -> Compound | None:
    """
    Get a compound by its name.
    
    Args:
        name (str): Name of the compound to retrieve.
    
    Returns:
        compound (dict): the compound chem-sep source data if found, otherwise None.
    """
    return db.get_compound(name).get_source("chemsep")


@deprecated("get_compound_names")
def get_compound_names() -> list:
    """
    Returns a list of all compound names.

    Returns:
        list: List of compound names.
    """
    return db.get_compound_names()


@deprecated("search_compounds")
def search_compounds(query: str) -> list:
    """
    Searches for compounds based on the query.

    Args:
        query (str): The search query.

    Returns:
        list: List of compound names that match the query.
    """
    return db.search_compounds(query)

# Restricting what can be imported from this module
__all__ = ["db", "get_compound", "get_compound_names", "search_compounds"]