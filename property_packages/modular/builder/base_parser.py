from property_packages.types import States
from typing import List, Dict

class BuildBase:

  def __init__(self) -> None:
    pass
  
  @staticmethod
  def serialise(compounds: List[Dict], valid_states: List[States]) -> Dict:
    """
    accepts
    - List of compounds

    returns
    - parsable segment of idaes dictionary
    """
    pass



