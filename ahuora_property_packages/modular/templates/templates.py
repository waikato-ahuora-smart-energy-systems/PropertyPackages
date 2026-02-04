from enum import Enum
from .pengrobinson import template as pr_template

class PropertyPackage(Enum):
    PengRobinson = "peng-robinson"
    # ...
    # Add other packages

    def get_template(self):
        return {
            PropertyPackage.PengRobinson: pr_template
            # ...
            # Add other mappings
        }.get(self, None)

    @classmethod
    def from_string(cls, package_name: str):
        # Find the enum member based on the string value
        for member in cls:
            if member.value == package_name.lower():
                return member
        return None