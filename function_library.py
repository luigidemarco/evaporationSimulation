from global_parameters import *
import numpy as np


def kinetic_energy(v):
    x = v * v * 1E-6
    return 0.5 * m * (x[:, 0] + x[:, 1])


def potential_energy(r):
    # Positions should be in Microns
    # return -Ud * 0.5 * (
    #        np.exp(-sigtrapinv * r[:, 0] * r[:, 0] * 1E-12) + np.exp(-sigtrapinv * r[:, 1] * r[:, 1] * 1E-12)) + Ud

    return Ud - Ud * np.exp(-sigtrapinv * (r[:, 0] * r[:, 0] + r[:, 1] * r[:, 1]) * 1E-12)


def trap_force(r):
    # Accepts an N x 2 matrix of positions, and returns an N x 2 matrix of forces
    # Positions should be in Microns
    # The force is returned in units of [kg*um/ms^2]
    #     X = -fmax*np.array([r[:, 0]*1E-6*np.exp(-sigtrapinv*r[:, 0]*r[:, 0]*1E-12),r[:,1]*1E-6*np.exp(-sigtrapinv*r[:, 1]*r[:, 1]*1E-12)])
    X = -fmax * np.multiply(np.exp(-sigtrapinv * (r[:, 0] * r[:, 0] + r[:, 1] * r[:, 1]) * 1E-12), r.T * 1E-6)
    return X.T


def evap_force(r, t):
    if t <= equilibrationTime:
        return r * 0
    else:
        f = -emax * (a - 2.0 * b * tleninv * r * 1E-6)
        f[:, 1] = 0

    if t < equilibrationTime + evaporationRamp:
        return f * (t - equilibrationTime) / evaporationRamp
    else:
        return f


def initialize_velocities(N):
    # Returns N initial velocities in [um/ms]
    return np.random.normal(0, sigmaVelocity, (N, 2)) * 1E3


def initialize_positions(N):
    # Returns N initial positions in [um]
    return np.random.normal(0, sigmaPosition, (N, 2)) * 1E6


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

        elasticCrossSection = elastic_cs(VRel)
        reactiveCrossSection = reactive_cs(VRel)

        totalCrossSection = reactiveCrossSection + elasticCrossSection
        inelasticProbability = reactiveCrossSection / (reactiveCrossSection + elasticCrossSection)

        # print(VRel * collisionProbabilityFactor * totalCrossSection)
        if VRel * collisionProbabilityFactor * totalCrossSection > P[j]:

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
                reactiveSuccess.append(k)
        else:
            pass

        j += 1

    return V, np.array(elasticSuccess).flatten(), np.array(reactiveSuccess).flatten()


def particles_out_of_bounds(r, bound):
    # Returns the indicies of particles that are beyond bound
    x = np.sum(abs(r), 1) > 2.0 * bound
    return np.array(np.nonzero(x)[0])


def write_params_file():
    savePath = filepath + filename + '.params'

    from datetime import datetime
    now = datetime.now()
    date_string = now.strftime("%Y/%m/%d")
    time_string = now.strftime("%H:%M:%S")

    f = open(savePath, 'w')
    f.write('Simulation name: {}\n'.format(filename))
    f.write('Simulation started: {} at {}\n\n'.format(date_string, time_string))
    print('Simulation started: {} at {}\n\n'.format(date_string, time_string))

    f.write('########### Simulation Parameters ###########\n\n')
    f.write(
        'T_MAX: {} ms\nTAU: {} ms\nWRITE_TIME: {} ms\nEQUILIBRATION: {} ms\nEVAP_RAMP: {} ms\nN: {}\n\n'.format(tmax,
                                                                                                                tau,
                                                                                                                writeEvery,
                                                                                                                equilibrationTime,
                                                                                                                evaporationRamp,
                                                                                                                N))
    f.write('BOUNDARY: {} um\nDIPOLE_CUTOFF: {} um\nCOLL_CUTOFF: {} nm\n\n'.format(bound, dipoleCutoff,
                                                                                   collisionCutoff * 1E3))

    f.write('########### Particle Parameters ###########\n\n')
    f.write(
        'MASS: {} AMU\nDIPOLE: {} D\nTEMP: {} nK\nINELASTIC_CS: {} um\nELASTIC_CS: {} um\n\n'.format(m / u, d / D2CM,
                                                                                                     T0 * 1E9,
                                                                                                     "na",
                                                                                                     "na"))

    f.write('########### Trap Parameters ###########\n\n')
    f.write('DEPTH: {} uK\nFREQ: {} Hz\nA: {}\nB: {}\n\n'.format(Ud / kB * 1E6, omega / (2.0 * np.pi), a, b))

    f.close()


def trap_potential_for_plot(r):
    # Positions should be in Microns
    #     return -Ud*0.5*(np.exp(-sigtrapinv*r[0]*r[0]*1E-12)+np.exp(-sigtrapinv*r[1]*r[1]*1E-12))
    return -Ud * np.exp(-sigtrapinv * (r[0] * r[0] + r[1] * r[1]) * 1E-12) + Ud * (
            a * tleninv * r[0] - b * tleninv * tleninv * r[0] * r[0] * 1E-6) * 1E-6
