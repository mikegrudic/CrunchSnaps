import numpy as np
from numpy import abs, copysign, sqrt
from numba import njit, vectorize

cubemap_directions = "forward", "left", "right", "up", "down", "backward"

def cubemapify(params):
    new_params = []
    for dir in cubemap_directions:
        d = params.copy()
        d.update({"cubemap_dir": dir})
        new_params.append(d)
    return new_params

@njit
def NormalizeVector(v):
    norm = 0
    for k in range(v.shape[0]):
        norm += v[k]*v[k]
    norm = 1./sqrt(norm)
    for k in range(v.shape[0]):
        v[k] *= norm

@vectorize
def NearestImage(x,boxsize):
    if abs(x) > boxsize/2: return -copysign(boxsize-abs(x),x)
    else: return x
