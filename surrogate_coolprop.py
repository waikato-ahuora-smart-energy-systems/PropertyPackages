from idaes.core.surrogate.pysmo.sampling import HammersleySampling
from pyfluids import Fluid, FluidsList, Input
from idaes.core.surrogate.pysmo.polynomial_regression import PolynomialRegression
from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer, PysmoSurrogate, PysmoPolyTrainer
import pandas as pd
from idaes.core.surrogate.sampling.data_utils import split_training_validation
from idaes.core.surrogate.plotting.sm_plotter import *
import time

bounds_info = [[273.15, 1000.0], [373.15, 1000000.0]]

init_data = HammersleySampling(data_input=bounds_info, number_of_samples=1000)

time_start = time.time()

result = {}

for [temp, press] in init_data.sample_points():
    # Units and descriptions for properties are https://github.com/portyanikhin/PyFluids?tab=readme-ov-file#properties-of-fluid-and-mixture-instances
    water = Fluid(FluidsList.Water).with_state(
        Input.pressure(press), Input.temperature(temp)
    )
    result[temp, press] = [
        temp,
        press,
        water.entropy, 
        water.enthalpy, 
        water.dynamic_viscosity, 
        water.kinematic_viscosity, 
        water.molar_mass,
        water.specific_volume,
    ]

df = pd.DataFrame(
    result.values(),
    columns=[
        "temperature",
        "pressure",
        "entropy", 
        "enthalpy", 
        "dynamic_viscosity", 
        "kinematic_viscosity", 
        "molar_mass",
        "specific_volume",
    ],
)

data_training, data_validation = split_training_validation(df, 0.8, seed=42) 


input_labels=["temperature", "pressure"]
output_labels=[
    "entropy", 
    "enthalpy", 
    "dynamic_viscosity", 
    "kinematic_viscosity", 
    "molar_mass",
    "specific_volume",
]

trainer = PysmoRBFTrainer(
    input_labels=input_labels,
    output_labels=output_labels,
    training_dataframe=data_training,
)

trained_surr = trainer.train_surrogate()

time_end = time.time()
print(f"Time taken to train surrogate: {time_end - time_start} seconds")

surr = PysmoSurrogate(trained_surr, input_labels, output_labels)
surr.save_to_file("coolprop_surrogate.json", overwrite=True)

surrogate_scatter2D(surr, data_validation, filename='scatter2D.pdf',show=False)
surrogate_parity(surr, data_validation, filename='parity.pdf',show=False)
surrogate_residual(surr, data_validation, filename='residual.pdf',show=False)

