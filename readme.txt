2D Evaporation Simulation

Written by:
Luigi De Marco
Kyle Matsuda

############################################################################################

Simulate the evaporation of particles in a two dimensional Gaussian trap.
Both elastic and inelastic collisions are treated using Monte Carlo.
Trap dynamics are treated through Verlet integration of the equations of motion.

Simulations are greatly sped up via a C++ extension to deal with collisions:
    C++ extension to speed up collision_check() in function_library.py.
    To build and install extension: "python setup.py"
    Need setuptools and distutils installed.
    For Windows, need Microsoft Visual C++ Compiler for Python 2.7.
    For Mac, need Xcode command line tools installed. May also need to run setup.py as root.
    After installing, should be able to "import krbcollision" in Python.

############################################################################################

To run:

    python main.py --input file.in --output output.out

The output file is a text file with columns:

    time [ms]       mean kinetic energy [J]     mean potential energy [J]       number       Elastic Collisions

However, if the "Nonequilibrium" parameter is specified, the kinetic energy for both the x and y axes are both contained
in the output file.

A parameters file is also generated that summarizes the input parameters and reports the runtime of the simulation.

############################################################################################

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
Nonequilibrium          Initial Temperature Scaling for X-Axis              0.0

Inelastic               Whether inelastic collisions occur                  True
ElasticCoeff            Elastic cross section polynomial coefficients       KRb @ 0.2 D
ReactiveCoeff           Reactive cross section polynomial coefficients      KRb @ 0.2 D
Collision               Specify "swave" or "differential" collision type    differential
DiffCSParams            Differential cross-section parameters
                        corresponding to Eq. 1 in ArXiV 1311.0429
                        Comma separated list: a, a', alpha, alpha'          KRb @ 0.2 D and 100 nK

Bound                   Boundary at which particles are removed [um]        150
CollisionCutOff         Cutoff at which a collision can occur [um]          0.05

Trap                    Specify: "Harmonic" or "Gaussian"                   Gaussian
Depth                   Trap depth [uK]                                     2.5
Freq                    Trap frequency [Hz]                                 35
A                       Evaporation gradient parameter                      0.0
B                       Evaporation curvature parameter                     0.0
EvaporationRamp         Timescale on which Grad./Curv. is turned on [ms]    20.0
---------------------------------------------------------------------------------

The elastic and inelastic collision cross sections are defined via a polynomial relation.
The input takes a series of comma separated coefficients, e.g:

    ElasticCoeff: -7.525E-6, 3.095E-4, -4.731E-3, 2.913E-2, -1.469E-3

and returns an Nth order polynomial function. For this example, the cross section would be

    -7.525E-6 * V ** 4 + 3.095E-4 * V ** 3 -4.731E-3 * V ** 2 + 2.913E-2 * V - 1.469E-3

The coefficients should be in units such that the relative velocity is specified in [um/ms]
and 2D cross section is returned in [um].

