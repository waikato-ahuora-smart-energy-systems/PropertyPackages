from compounds import deprecated
from compounds.CompoundRegistry import CompoundRegistry
from compounds.Compound import Compound

__author__ = "Mahaki Leach"


#  Simplified wrapper for back-end interactions
class RegistryView:
    def __init__(self, compound_registry):
        self._registry = compound_registry

    def register_package(self, pkg):
        self._registry.queue_package(pkg)

    def register_compound(self, name, source, data):
        self._registry.queue_compound(name, source, data)

    def bind(self, compound_name, package_name):
        self._registry.queue_bind(compound_name, package_name)


# Simplified wrapper for front-end integration
class RegistrySearch:

    def __init__(self, compound_registry):
        self._registry = compound_registry

    def search_compounds(self, query, package_filters=None, filter_strict=False):
        return self._registry.search_compounds(query, package_filters, filter_strict)

    def get_compound_names(self):
        return self._registry._get_compound_names()

    def get_compound(self, name):
        return self._registry._get_compound(name)

    def get_supported_packages(self, compounds, strict=True):
        return self._registry.get_supported_packages(compounds, strict)

    def get_supported_compounds(self, packages, strict=True):
        return self._registry.get_supported_compounds(packages, strict)


# Main registry instance
_registry = CompoundRegistry()

# Backend registry view
registry_view = RegistryView(_registry)

# Frontend registry view
db = RegistrySearch(_registry)

_registry._discover_loaders()

def init_registry():
    ### this ensures that all loaders are imported and registered before the registry is built
    from compounds.loaders import loaders_list

    # Import all loaders
    for loader in loaders_list:
        loader(registry_view)

    # Building the registry
    _registry._build()


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

# To force build module on import
def __getattribute__(self, attr):
        if attr == "db":
            init_registry()
            return _registry
        else:
            return super().__getattribute__(attr)

# Restricting what can be imported from this module
__all__ = ["db", "get_compound", "get_compound_names", "search_compounds"]