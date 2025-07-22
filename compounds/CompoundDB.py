from compounds import deprecated
from CompoundRegistry import CompoundRegistry
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


registry = CompoundRegistry()
registry_view = RegistryView(registry)
registry_search = RegistrySearch(registry)

registry._discover_loaders()

def main():
    ### this ensures that all loaders are imported and registered before the registry is built
    from compounds.loaders import loaders_list

    for loader in loaders_list:
        print("test")
        loader(registry_view)

    registry._build()


@deprecated("get_compound")
def get_compound(name: str) -> Compound | None:
    """
    Get a compound by its name.
    
    Args:
        name (str): Name of the compound to retrieve.
    
    Returns:
        compound (dict): the compound chem-sep source data if found, otherwise None.
    """
    return registry_search.get_compound(name).get_source("chemsep")


@deprecated("get_compound_names")
def get_compound_names() -> list:
    """
    Returns a list of all compound names.

    Returns:
        list: List of compound names.
    """
    return registry_search.get_compound_names()


@deprecated("search_compounds")
def search_compounds(query: str) -> list:
    """
    Searches for compounds based on the query.

    Args:
        query (str): The search query.

    Returns:
        list: List of compound names that match the query.
    """
    return registry_search.search_compounds(query)

if __name__ == "__main__":
    main()