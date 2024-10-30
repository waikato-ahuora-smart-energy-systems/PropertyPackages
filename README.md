To install this:

```sh
pip install git+https://github.com/waikato-ahuora-smart-energy-systems/PropertyPackages.git
```

To use


```python
import property_packages
property_packages.build_package("helmholtz",["h2o"])

from compounds import CompoundDB as db

benzene = db.get_compound("benzene")
benzene["AbsEntropy"].value # "269300"
benzene["AbsEntropy"].unit # "J/kmol/K" 
```

# Next steps
Implement something to translate between apis, e.g for if the user specifies in Temperature and mass flow but the property package uses flow_mol and enthalpy

```py
property_package.build_state(
    "helmholtz",
    "h2o",
    {
        "temperature",
        "pressure",
        "mass_flow",
    }
) -> {
    "enthalpy":
    "pressure":
    "quality":
    "flow_mol":
}
```

Support different phases:

```py
property_packages.build_package("helmholtz",["h2o"],["vap","liq"])
```