2D Evaporation Simulation
Luigi De Marco
############################################################################################

Simulate the evaporation of particles in a two dimensional Gaussian trap.
Both elastic and inelastic collisions are treated using Monte Carlo.
Trap dynamics are treated through Verlet integration of the equations of motion.




To run:

    python main.py --input file.in --output output.out

The output file is a text file with columns:

    time [ms]       mean kinetic energy [J]     mean potential energy [J]       number

A parameters file is also generated that summarizes the input parameters and reports the runtime of the simulation.

**NOTE: Running 2000 particles with no inelastic loss for 1s takes ~18 hours on a decent computer
        Runtime scales roughly quadratically with N and linearly with number of time steps




The input file should be a text file with lines of the form

    Parameter1: Value1
    Parameter2: Value2
    Parameter3: Value3
    etc...

The parameters are not case sensitive.

Possible parameters are:

-----------------------------------------------------------------------------------
Parameter               Value                                               Default
---------               -----                                               -------

Name                    A Name for your simulation
Comment                 A Comment

TMax                    Length of simulation [ms]                           10
Tau                     Time step [ms]                                      0.002
WriteEvery              Time step for writing data [ms]                     0.1
EquilibrationTime       Time before evaporation starts [ms]                 0.0

N                       Number of particles                                 20
M                       Mass of particle [amu]                              127
T                       Starting temperature [nK]                           300
Inelastic               Whether inelastic collisions occur                  True

Bound                   Boundary at which particles are removed [um]        150
CollisionCutOff         Cutoff at which a collision can occur [um]          0.05

Depth                   Trap depth [uK]                                     2.5
Freq                    Trap frequency [Hz]                                 35
A                       Evaporation gradient parameter                      0.0
B                       Evaporation curvature parameter                     0.0
EvaporationRamp         Timescale on which Grad./Curv. is turned on [ms]    20.0
---------------------------------------------------------------------------------

The elastic and inelastic collision cross sections are defined in global_parameters.py.
The functions take in a relative velocity [um/ms] and return a 2D cross section [um].

