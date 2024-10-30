from typing import List, Dict
from compounds.Compound import load_compound

class BuildBase:

  def __init__(self) -> None:
    pass
  
  @staticmethod
  def serialise( compounds: List[Dict]) -> Dict:
    """
    accepts
    - List of compounds

    returns
    - parsable segment of idaes dictionary
    """
    pass



