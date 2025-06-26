from loaders import loaders
from typing import Dict

# Loads all the compounds, and enables some basic queries for/on the compounds.
# Handy to avoid extra reloading of the compounds, but does consume more memory.
# Right now, the compounds are stored as xml files, but we could eventually
# move to a database, hdf5, or use a redis cache or something else as the backend of this.
# https://chatgpt.com/share/67196214-92a0-8005-8d2d-68713f50c9d6


# compounds["source"]["name"]["parameter"]
all_compounds = {}

# iterate through all files in the data_files folder
# and create a Compound object for each file
for file in os.listdir(os.path.dirname(__file__) + "/data_files/chemsep"):
    if file.endswith(".xml"):
        compound_name = file[:-4]
        compound = load_compound(compound_name)
        all_compounds[compound_name] = compound

# Calling all loaders
for name, loader in loaders.items():
    try:
        all_compounds["name"] = loader()
    except Exception as e:
        print(f"Error loading compounds with loader '{loader}': {e}")

def get_compound(compound_name: str, ) -> Dict:
    """
    Returns a Compound object for the given compound name.

    Args:
        compound_name (str): Name of the compound.

    Returns:
        Compound: The Compound object.
    """
    return all_compounds[compound_name.lower()]

def get_compound_names() -> list:
    """
    Returns a list of all compound names.

    Returns:
        list: List of compound names.
    """
    return list(all_compounds.keys())

def search_compounds(query: str) -> list:
    """
    Searches for compounds based on the query.

    Args:
        query (str): The search query.

    Returns:
        list: List of compound names that match the query.
    """
    return [name for name in all_compounds.keys() if query.lower() in name.lower()]