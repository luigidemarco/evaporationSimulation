from global_parameters import *
import numpy as np
import sys
import getopt
from os import path


def parse_inputs(sysargv):
    try:
        opts, args = getopt.getopt(sysargv, "hi:o:", ["input=", "output="])
    except getopt.GetoptError:
        print "Unknown option specified. Specify -h for usage"
        sys.exit()

    input_file = ''
    output_file = ''

    for opt, arg in opts:
        if opt in "-h":
            print("main.py -i <input file> -o <output file>")
            sys.exit()
        elif opt in ["-i", "--input"]:
            input_file = arg
        elif opt in ["-o", "--output"]:
            output_file = arg
        else:
            pass

    if not input_file:
        print('Input file not specified. Specify an input file to run simulation.')
        sys.exit()
    elif not output_file:
        output_file = ".".join(input_file.split('.')[0:-1]) + ".out"
        params_file = ".".join(input_file.split('.')[0:-1]) + ".params"
        print('Output file not specified. Output will be written to {}'.format(output_file))
    else:
        params_file = ".".join(output_file.split('.')[0:-1]) + ".params"

    """ ------- Make sure not to overwrite file ------- """
    if path.exists(output_file):
        overwrite = raw_input("Are you sure you want to overwrite {}? Data will be lost! (y/N): ".format(output_file))

        if overwrite.lower() == 'y':
            print("Data in {} will be overwritten.".format(output_file))
        else:
            raise SystemExit("Simulation aborted. Change output file path.")
        pass
    return input_file, output_file, params_file


def set_global_parameters(input_file):
    f = open(input_file)
    global_parameter_keys = global_parameters.keys()
    meta_data_keys = meta_data.keys()

    for line in f:
        if line.strip():

            stripped_nl = line.replace('\n', "")
            line_split = stripped_nl.split(':')

            parameter = line_split[0].strip().lower()
            value = line_split[1].strip()

            try:
                value = float(value)
            except ValueError:
                pass

            if parameter not in global_parameter_keys + meta_data_keys:
                raise AttributeError('Found undefined parameter \"{}\" in {}.'.format(parameter, input_file))

            elif parameter in global_parameter_keys:
                if parameter == 'freq':
                    global_parameters[parameter] = 2.0*np.pi*value

                elif parameter == 'depth':
                    global_parameters[parameter] = value * kB * 1E-6

                elif parameter == 't':
                    global_parameters[parameter] = value * 1E-9

                elif parameter == 'm':
                    global_parameters[parameter] = value * u

                elif parameter == 'n':
                    global_parameters[parameter] = int(value)

                elif parameter == 'inelastic':

                    if isinstance(value, float):
                        global_parameters[parameter] = bool(value)
                    elif value.lower() in ('false', 'no', 'f'):
                        global_parameters[parameter] = False
                    elif value.lower() in ('true', 'yes', 't'):
                        global_parameters[parameter] = True
                    else:
                        raise ValueError('Expected a boolean input for \"inelastic\" parameter,'
                                         ' got {} instead.'.format(value))

                elif parameter == 'elasticcoeff' or parameter == 'reactivecoeff':

                    if type(value) == float:
                        global_parameters[parameter] = [value]
                    else:
                        value = value.replace(" ", "")
                        value = value.split(",")
                        value = map(float, value)
                        global_parameters[parameter] = value

                else:
                    global_parameters[parameter] = value

            elif parameter in meta_data_keys:
                meta_data[parameter] = value
    f.close()
    return 0


def calculate_derived_parameters():

    derived_parameters['sigtrapinv'] = global_parameters["m"] * global_parameters['freq'] * global_parameters['freq'] \
                 / (2. * global_parameters['depth'])

    derived_parameters['fmax'] = global_parameters['m'] * global_parameters['freq'] * global_parameters['freq']

    derived_parameters['tleninv'] = np.sqrt(derived_parameters['fmax'] / global_parameters['depth'])

    derived_parameters['emax'] = derived_parameters['tleninv'] * global_parameters['depth']

    derived_parameters['sigmaVelocity'] = np.sqrt(kB * global_parameters['t'] / global_parameters['m'])
    derived_parameters['sigmaPosition'] = derived_parameters['sigmaVelocity'] / global_parameters['freq']

    derived_parameters['collisionProbabilityFactor'] = global_parameters['tau'] / (np.pi * global_parameters['collisioncutoff']
                                                             * global_parameters['collisioncutoff'])

    derived_parameters['time'] = np.arange(0, global_parameters['tmax'] + global_parameters['tau'], global_parameters['tau'])
    derived_parameters['nT'] = len(derived_parameters['time'])
    derived_parameters['writeEveryInv'] = 1.0 / global_parameters['writeevery']

def kinetic_energy(v):
    x = v * v * 1E-6
    return 0.5 * global_parameters['m'] * (x[:, 0] + x[:, 1])


def potential_energy(r):
    # Positions should be in Microns
    # return -Ud * 0.5 * (
    #        np.exp(-sigtrapinv * r[:, 0] * r[:, 0] * 1E-12) + np.exp(-sigtrapinv * r[:, 1] * r[:, 1] * 1E-12)) + Ud

    return global_parameters['depth'] - global_parameters['depth'] \
           * np.exp(-derived_parameters['sigtrapinv'] * (r[:, 0] * r[:, 0] + r[:, 1] * r[:, 1]) * 1E-12)


def trap_force(r):
    # Accepts an N x 2 matrix of positions, and returns an N x 2 matrix of forces
    # Positions should be in Microns
    # The force is returned in units of [kg*um/ms^2]
    #     X = -fmax*np.array([r[:, 0]*1E-6*np.exp(-sigtrapinv*r[:, 0]*r[:, 0]*1E-12),
    #     r[:,1]*1E-6*np.exp(-sigtrapinv*r[:, 1]*r[:, 1]*1E-12)])

    X = -derived_parameters['fmax'] * np.multiply(np.exp(-derived_parameters['sigtrapinv'] * (r[:, 0] * r[:, 0] + r[:, 1] * r[:, 1]) * 1E-12), r.T * 1E-6)
    return X.T


def evap_force(r, t):
    if t <= global_parameters['equilibrationtime']:
        return r * 0
    else:
        f = -derived_parameters['emax'] * (global_parameters['a'] -
                                           2.0 * global_parameters['b'] * derived_parameters['tleninv'] * r * 1E-6)
        f[:, 1] = 0

    if t < global_parameters['equilibrationtime'] + global_parameters['evaporationramp']:
        return f * (t - global_parameters['equilibrationtime']) / global_parameters['evaporationramp']
    else:
        return f


def initialize_velocities(N):
    # Returns N initial velocities in [um/ms]
    return np.random.normal(0, derived_parameters['sigmaVelocity'], (N, 2)) * 1E3


def initialize_positions(N):
    # Returns N initial positions in [um]
    return np.random.normal(0, derived_parameters['sigmaPosition'], (N, 2)) * 1E6


def make_cross_section(coefficients):
    # Given a list of coefficients [a1, a2, a3, ..., aN] returns a function
    # that is the polynomial a1*x**N + a2*x**(N-1) + ... + aN
    def polynomial(v):
        result = 0
        for c in coefficients:
            result = result * v + c
        return result
    return polynomial


def collision_check(r, colcut=0.1):
    # Takes a N x 2 matrix of positions, and returns elastic and inelastic loss candidates [List of pairs]

    N = r.shape[0]

    ####Calculate all interparticle distances
    x1, x2 = np.meshgrid(r[:, 0], r[:, 0])
    dx = x2 - x1

    tx = np.argwhere((np.abs(dx) < colcut) * np.triu(np.ones((N, N)), 1))

    squareColIndex = []
    collIndex = []
    coords = []

    for k in tx:
        dy = r[k[0], 1] - r[k[1], 1]
        if abs(dy) < colcut:
            squareColIndex.append(k)
            coords.append([dx[k[0], k[1]], dy])

    nSqCol = len(squareColIndex)
    if nSqCol == 0:
        return []
    else:
        pass

    coords = np.array(coords)
    rij2 = coords[:, 0] * coords[:, 0] + coords[:, 1] * coords[:, 1]

    #### Calculate the dipole-dipole forces for particles separated by less than dipcut
    for k in range(nSqCol):
        if rij2[k] < colcut * colcut:
            collIndex.append(squareColIndex[k])

    return np.array(collIndex)


def collision_montecarlo(colList, V):
    n = len(colList)

    pEvI = np.random.random(n)  # Probabilities for reactive collision
    P = np.random.random(n)  # Probabilities for the montecarlo
    reactiveSuccess = []
    elasticSuccess = []

    j = 0
    for k in colList:
        dVx = V[k[0], 0] - V[k[1], 0]
        dVy = V[k[0], 1] - V[k[1], 1]

        VRel = np.sqrt(dVx * dVx + dVy * dVy)

        elasticCrossSection = cross_sections['elastic'](VRel)
        reactiveCrossSection = cross_sections['reactive'](VRel)

        totalCrossSection = reactiveCrossSection + elasticCrossSection
        inelasticProbability = reactiveCrossSection / (reactiveCrossSection + elasticCrossSection)

        # print(VRel * collisionProbabilityFactor * totalCrossSection)
        if VRel * derived_parameters['collisionProbabilityFactor'] * totalCrossSection > P[j]:

            if pEvI[j] > inelasticProbability:
                # Elastic collisions
                elasticSuccess.append(k)

                # Generate an angle between pi and pi/2 at random for the rotation of the velocities
                costheta = -np.random.random(1)[0]
                sintheta = np.sqrt(1 - costheta * costheta)

                VCOM = 0.5 * np.array([V[k[0], 0] + V[k[1], 0], V[k[0], 1] + V[k[1], 1]])
                VDIFROT = np.array([(dVx * costheta - dVy * sintheta), (dVx * sintheta + dVy * costheta)])

                V[k[0], :] = VCOM + 0.5 * VDIFROT
                V[k[1], :] = VCOM - 0.5 * VDIFROT

            else:
                # Reactive collisions
                if global_parameters['inelastic']:
                    reactiveSuccess.append(k)
        else:
            pass

        j += 1

    return V, np.array(elasticSuccess).flatten(), np.array(reactiveSuccess).flatten()


def particles_out_of_bounds(r, bound):
    # Returns the indicies of particles that are beyond bound
    x = np.sum(abs(r), 1) > 2.0 * bound
    return np.array(np.nonzero(x)[0])


def write_params_file(params_file):

    from datetime import datetime
    now = datetime.now()
    date_string = now.strftime("%Y/%m/%d")
    time_string = now.strftime("%H:%M:%S")

    f = open(params_file, 'w')
    if meta_data['name']:
        simulation_name = meta_data['name']
    else:
        simulation_name = params_file.split('/')[-1].split('.')[0]

    f.write('Simulation name: {}\n'.format(simulation_name))
    if meta_data['comment']:
        f.write('Comments: {}\n'.format(meta_data['comment']))

    f.write('Simulation started: {} at {}\n\n'.format(date_string, time_string))
    print('Simulation started: {} at {}\n\n'.format(date_string, time_string))

    f.write('########### Simulation Parameters ###########\n\n')
    f.write(
        'T_MAX: {} ms\nTAU: {} ms\nWRITE_TIME: {} ms\n'
        'EQUILIBRATION: {} ms\nN: {}\n\n'.format(global_parameters['tmax'],
                                                                   global_parameters['tau'],
                                                                   global_parameters['writeevery'],
                                                                   global_parameters['equilibrationtime'],
                                                                   global_parameters['n']))

    f.write('BOUNDARY: {} um\nCOLL_CUTOFF: {} nm\n\n'.format(global_parameters['bound'],
                                                             global_parameters['collisioncutoff'] * 1E3))

    f.write('########### Particle Parameters ###########\n\n')
    f.write(
        'MASS: {} AMU\nTEMP: {} nK\nINELASTIC_COLL: {}\n\n'.format(global_parameters['m'] / u,
                                                                      global_parameters['t'] * 1E9,
                                                                      global_parameters['inelastic']))
    f.write('ELASTICCOEFF: {}\nREACTIVECOEFF: {}\n\n'.format(global_parameters['elasticcoeff'],
                                                             global_parameters['reactivecoeff']))

    f.write('########### Trap Parameters ###########\n\n')
    f.write('DEPTH: {} uK\nFREQ: {} Hz\nA: {}\nB: {}\nEVAP_RAMP: {} ms\n\n'.format(
                                                                 global_parameters['depth'] / kB * 1E6,
                                                                 global_parameters['freq'] / (2.0 * np.pi),
                                                                 global_parameters['a'],
                                                                 global_parameters['b'],
                                                                 global_parameters['evaporationramp'],))

    f.close()


def trap_potential_for_plot(r):
    # Positions should be in Microns
    #     return -Ud*0.5*(np.exp(-sigtrapinv*r[0]*r[0]*1E-12)+np.exp(-sigtrapinv*r[1]*r[1]*1E-12))
    return -global_parameters['depth']  * np.exp(-derived_parameters['sigtrapinv'] * (r[0] * r[0] + r[1] * r[1]) * 1E-12) \
           + global_parameters['depth']  * (
            global_parameters['a'] * derived_parameters['tleninv'] * r[0] - global_parameters['b']
            * derived_parameters['tleninv'] * derived_parameters['tleninv'] * r[0] * r[0] * 1E-6) * 1E-6
