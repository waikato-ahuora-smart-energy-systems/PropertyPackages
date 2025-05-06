from pyomo.common.fileutils import this_file_dir
import idaes
from idaes.models.properties.general_helmholtz.components.parameters import auto_register
from idaes.models.properties.general_helmholtz import helmholtz_functions 
from idaes.models.properties.general_helmholtz import helmholtz_state
import os

def register_compounds():
    idaes.cfg.properties.helmholtz.parameter_file_path = os.path.join(this_file_dir(), "")
    helmholtz_functions._data_dir = idaes.cfg.properties.helmholtz.parameter_file_path
    helmholtz_functions.helmholtz_data_dir = idaes.cfg.properties.helmholtz.parameter_file_path
    helmholtz_state._data_dir = idaes.cfg.properties.helmholtz.parameter_file_path 
    auto_register()