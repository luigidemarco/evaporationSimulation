from function_library import *
import time as clock


write_params_file()
resFile = open(filepath + filename + '.out', 'w')

# Initialize the positions and velocities
R0 = initialize_positions(N)
V0 = initialize_velocities(N)

# R0 = np.array([[10,0.00],[-10,-0.00]])
# V0 = np.array([[-15,0],[15,0]])

# Remove particles initialized out of bounds or within 0.02 um of each other
particlesOOB = particles_out_of_bounds(R0, bound)
particlesCollision = collision_check(R0, 0.02)
particlesRemoval = np.unique(np.concatenate((particlesCollision, particlesOOB)))

R0 = np.delete(R0, particlesRemoval, 0)
V0 = np.delete(V0, particlesRemoval, 0)


""" ------- Calculated the zeroth time step ------- """
start = clock.time()
F0 = trap_force(R0) + evap_force(R0, 0.0)

T = np.mean(kinetic_energy(V0))
U = np.mean(potential_energy(R0))
pN = R0.shape[0]

resFile.write('{}\t{:.4e}\t{:.4e}\t{:0.0f}'.format(time[0], T, U, pN))

# trajectory = np.zeros((nT, N, 2))
# trajectory[0, :, :] = R0


""" ------- Verlet integration for the dynamics ------- """
for k in range(1, nT):
    # Calculate the Verlet Equations

    R = R0 + tau * V0 + tau * tau * F0 * 0.5 / m

    particlesCollision = collision_check(R, collisionCutoff)
    particlesOOB = particles_out_of_bounds(R, bound)

    V0, elasticPairs, inelasticPairs = collision_montecarlo(particlesCollision, V0)

    # Removal of particles

    particlesRemoval = np.unique(np.concatenate((inelasticPairs, particlesOOB)))
    #     particlesRemoval = []

    R = np.delete(R, particlesRemoval, 0)
    V0 = np.delete(V0, particlesRemoval, 0)
    F0 = np.delete(F0, particlesRemoval, 0)

    F = trap_force(R) + evap_force(R, time[k])
    V = V0 + tau * (F0 + F) * 0.5 / m

    # Verlet integration complete
    # Store important time step info at this point

    # trajectory[k, :, :] = R
    #     Note that the trajectories shouldn't be saved if N is larger than 10 or so

    if time[k] * writeEveryInv % 1 == 0:
        T = np.mean(kinetic_energy(V))
        U = np.mean(potential_energy(R))
        pN = R.shape[0]
        resFile.write('\n{}\t{:.4e}\t{:.4e}\t{:0.0f}'.format(time[k], T, U, pN))

    # Reinitialize the variables
    R0 = R
    V0 = V
    F0 = F

end = clock.time()
resFile.close()

with open(filepath + filename + '.params', 'a') as f:
    f.write('#############################################\n\n')
    from datetime import datetime

    now = datetime.now()
    date_string = now.strftime("%Y/%m/%d")
    time_string = now.strftime("%H:%M:%S")
    f.write('Simulation completed: {} at {}\nRuntime: {} s'.format(date_string, time_string, end - start))

print("Simulation complete.\nRuntime: {} s.".format(end - start))
