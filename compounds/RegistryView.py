#  Simplified wrapper for back-end interactions
class RegistryView:
    def __init__(self, compound_registry):
        self._registry = compound_registry

    def register_package(self, pkg):
        self._registry.queue_package(pkg)

    def register_compound(self, name, source, data):
        self._registry.queue_compound(name, source, data)

    def bind(self, compound_name, package_name):
        self._registry.queue_binding(compound_name, package_name)