""" Global parameters file for evaporation simulation """
import numpy as np

filename = "testo"
filepath = "./scratch/"

""" --------------------------------------- Constants """

eps0 = 8.8541878128E-12
kB = 1.380649E-23
u = 1.66053906660E-27
D2CM = 3.33564E-30

""" --------------------------------------- Simulation parameters """

tmax = 0.5  # Time in ms
tau = 0.002  # Time step in ms
writeEvery = 0.1  # Timestep to write the data in ms
N = 8  # Number of particles

bound = 150.
dipoleCutoff = 0.6
collisionCutoff = 0.05

equilibrationTime = 20.0
evaporationRamp = 20.0

""" --------------------------------------- Particle parameters """

m = 127 * u  # Mass in kg
d = 0.2 * D2CM  # Dipole moment in C m
T0 = 600E-9  # Initial temperature

reactiveCrossSection = 4.28E-4  # 2D Cross-section @ 0.2D in [um] PHYSICAL REVIEW A 81, 060701(R) (2010)
elasticCrossSection = reactiveCrossSection * 127.0


""" --------------------------------------- Trap parameters """

Ud = 2.5 * 1E-6 * kB  # Trap depth in J
omega = 2 * np.pi * 35


""" --------------------------------------- Evaporation parameters """

a = 0.1 * 0
b = 0.1 * 0

""" --------------------------------------- Derived parameters """

sigtrapinv = m * omega * omega / (2. * Ud)
fmax = m * omega * omega

tleninv = np.sqrt(fmax / Ud)
emax = tleninv * Ud

sigmaVelocity = np.sqrt(kB * T0 / m)
sigmaPosition = sigmaVelocity / omega

inelasticProbability = reactiveCrossSection / (reactiveCrossSection + elasticCrossSection)
totalCrossSection = reactiveCrossSection + elasticCrossSection

collisionProbabilityFactor = tau/(np.pi*collisionCutoff*collisionCutoff)

time = np.arange(0, tmax + tau, tau)
nT = len(time)
writeEveryInv = 1.0 / writeEvery