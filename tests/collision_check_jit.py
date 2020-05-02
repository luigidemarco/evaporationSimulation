import numpy as np
from numba import jit
import time

import sys
sys.path.append('..')
from function_library import collision_check

# Comment @jit to disable compilation
@jit(nopython=True)
def find_pairs(xx, colcut=0.1):
    # Collect particles closer than colcut in x or in y
    squareColIndex = []
    coords = []
    for i, a in enumerate(xx):
        for j, b in enumerate(xx[(i+1):-1]):
            dx = np.abs(b[0] - a[0])
            if dx < colcut:
                dy = np.abs(b[1] - a[1])
                if dy < colcut:
                    squareColIndex.append([i, i+1+j])
                    coords.append([dx, dy])
            else:
                break
    return (squareColIndex, coords)

def collision_check_JIT(r, colcut=0.1):
    # Takes a N x 2 matrix of positions, and returns elastic and inelastic loss candidates [List of pairs]

    # Sort r based on the x-coordinate
    # see https://stackoverflow.com/a/30623882
    s = np.lexsort(np.fliplr(r).T)
    xx = r[s]

    (squareColIndex, coords) = find_pairs(xx, colcut)
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
print("Original function (Python 3), {} trials: {:.2f}s".format(nreps, t2 - t1))

t1 = time.time()
for i in range(nreps):
    r = np.random.normal(0, sigma, (N, 2)) * 1E6
    collision_check_JIT(r)
t2 = time.time()
print("Numba jit'ed function, {} trials: {:.2f}s".format(nreps, t2 - t1))

failed = False
for i in range(ncomps):
    r = np.random.normal(0, sigma, (N, 2)) * 1E6

    check = collision_check_JIT(r) == collision_check(r)
    # print(check)

    if not np.all(check):
        print("{} random comparisons failed! At test {}.\n".format(ncomps, i))
        failed = True
        break

if not failed:
    print("{} random comparisons succeeded!\n".format(ncomps))