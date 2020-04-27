""" Global parameters file for evaporation simulation """
import numpy as np

""" --------------------------------------- Constants """
eps0 = 8.8541878128E-12
kB = 1.380649E-23
u = 1.66053906660E-27
D2CM = 3.33564E-30

global_parameters = {
    # --------------------------------------- Simulation parameters """
    "tmax": 10,  # Time in ms
    "tau": 0.002,  # Time step in ms
    "equilibrationtime": 0.0,
    "writeevery": 0.1,  # Timestep to write the data in ms

    "n": 20,

    "bound": 150,
    "collisioncutoff": 0.05,

    "nonequilibrium": 0.0,

    # --------------------------------------- Particle parameters
    "m": 127 * u,  # Mass in kg
    "t": 300E-9,  # Initial temperature
    "inelastic": True,  # Set false to turn off inelastic collisions

    "elasticcoeff": [-7.525E-6, 3.095E-4, - 4.731E-3, 2.913E-2, - 1.469E-3],
    "reactivecoeff": [5.763E-5, - 7.05277E-5],

    # --------------------------------------- Trap parameters
    "depth": 2.5 * 1E-6 * kB,  # Trap depth in J
    "freq": 2 * np.pi * 35,

    # --------------------------------------- Evaporation parameters
    "a": 0.0,
    "b": 0.0,
    "evaporationramp": 20.0
}

meta_data = {"name": '',
             "comment": '',
             }

cross_sections = {"elastic": [],
                  "reactive": []
                  }

derived_parameters = {
    'sigtrapinv': 0,
    'fmax': 0,
    'tleninv': 0,
    'emax': 0,
    'sigmaVelocity': 0,
    'sigmaPosition': 0,
    'collisionProbabilityFactor': 0,
    'time': 0,
    'nT': 0,
    'writeEveryInv': 0
}
