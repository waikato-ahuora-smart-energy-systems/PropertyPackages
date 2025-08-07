# Simplified wrapper for front-end integration
class RegistrySearch:

    def __init__(self, compound_registry):
        self._registry = compound_registry

    def search_compounds(self, query, package_filters=None, filter_strict=False):
        return self._registry._search_compounds(query, package_filters, filter_strict)

    def get_compound_names(self):
        return self._registry._get_compound_names()

    def get_compound(self, name):
        return self._registry._get_compound(name)

    def get_supported_packages(self, compounds, strict=True):
        return self._registry._get_supported_packages(compounds, strict)

    def get_supported_compounds(self, packages, strict=True):
        return self._registry._get_supported_compounds(packages, strict)