
from pprint import pprint
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from compounds.CompoundDB import get_compound
from .templates.templates import PropertyPackage
from typing import List

"""
Creates an idaes property package, to follow the template found at:
https://idaes-pse.readthedocs.io/en/stable/explanations/components/property_package/general/index.html

"""

def build_config(property_package_name, compound_names: List[str]) -> dict[str,any]:

  # Build list of compound objects

  compounds = []

  for compound in compound_names:
    # Verifying compound name
    try:
        compounds.append(get_compound(compound))
    except:
       raise ValueError(f"Invalid compound name {compound}")

  # Verifying property package name

  package = PropertyPackage.from_string(property_package_name)

  if package is None:
      raise ValueError("Invalid property package name")

  # Retrieve property package template

  template = package.get_template()
  
  # Building template

  new_template = {}

  for key, obj in template.items():
    # Call the parse method on each object and update the template
    new_template[key] = obj.serialise(compounds)

  pprint(new_template)

  # Building property package and returning

  return GenericParameterBlock(**new_template)
