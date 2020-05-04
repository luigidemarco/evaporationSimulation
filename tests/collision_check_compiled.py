import numpy as np
import time

import sys
sys.path.append('..')
from function_library import collision_check

import krbcollision

def collision_check_compiled(r, colcut=0.1):
    # Takes a N x 2 matrix of positions, and returns elastic and inelastic loss candidates [List of pairs]

    # Sort r based on the x-coordinate
    # see https://stackoverflow.com/a/30623882
    s = np.lexsort(np.fliplr(r).T)
    xx = r[s]

    [squareColIndex, coords] = krbcollision.find_pairs(xx.tolist(), colcut)
    squareColIndex = [[min(s[sci[0]], s[sci[1]]), max(s[sci[0]], s[sci[1]])] for sci in squareColIndex]

    if not squareColIndex:
        return []

    squareColIndex = np.array(squareColIndex)
    coords = np.array(coords)

    collIndex = []
    rij2 = np.sum(np.square(coords), axis=-1)

    for r, sci in zip(rij2, squareColIndex):
        if r < colcut * colcut:
            collIndex.append(sci)

    return np.array(collIndex)

N = 2000
sigma = 20E-6
nreps = 1000
ncomps = 500

t1 = time.time()
for i in range(nreps):
    r = np.random.normal(0, sigma, (N, 2)) * 1E6
    collision_check(r)
t2 = time.time()
print "Original function, {} trials: {:.2f}s".format(nreps, t2 - t1)

t1 = time.time()
for i in range(nreps):
    r = np.random.normal(0, sigma, (N, 2)) * 1E6
    collision_check_compiled(r)
t2 = time.time()
print "C extension, {} trials: {:.2f}s".format(nreps, t2 - t1)

failed = False
for i in range(ncomps):
    r = np.random.normal(0, sigma, (N, 2)) * 1E6

    check = collision_check_compiled(r) == collision_check(r)
    # print check

    if not np.all(check):
        print "{} random comparisons failed! At test {}.\n".format(ncomps, i)
        failed = True
        break

if not failed:
    print "{} random comparisons succeeded!\n".format(ncomps)