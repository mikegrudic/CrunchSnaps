import numpy as np
from numpy import abs, copysign, sqrt
from numba import njit, vectorize
from glob import glob
from natsort import natsorted
import h5py
from os.path import abspath, exists

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

def get_snapshot_time_dict(snaps):
    snaps = natsorted(snaps)
    all_snaps = natsorted(glob(snaps[0].split("snapshot")[0] + "snapshot*.hdf5")) # look for other snapshots in same directory
    if len(snaps) == 1: # if we only have one snapshot, keep it simple and just open the file. otherwise we will do some fancy stuff to avoid opening multiple files for the timeline
        return {snapnum_from_path(snaps[0]) : h5py.File(snaps[0],'r')["Header"].attrs["Time"]}
        
    # don't yet know what the snapshot times are - get the snapshot times in a prepass
    snaptimes = []
    snapnums = []

    snaptimes_path = abspath(snaps[0]).split("/snapshot")[0] + "/.snapshot_times"
    
    if exists(snaptimes_path):
        snapnums, snaptimes = np.atleast_2d(np.loadtxt(snaptimes_path).T)
    snaptime_dict = dict(zip(snapnums, snaptimes))

    do_snapshot_pass = False
    if np.any([not (snapnum_from_path(s) in snaptime_dict.keys()) for s in snaps]): # check if we have a snapshot missing from the dictionary, if so we must do a pass
        print("Sinkvis2 getting snapshot times...")
        for s in all_snaps:
#             snapnums.append(snapnum_from_path(s))
             with h5py.File(s, 'r') as F:
                 snaptime_dict[snapnum_from_path(s)] = F["Header"].attrs["Time"]
        np.savetxt(snaptimes_path, np.c_[[k for k in snaptime_dict.keys()], [v for v in snaptime_dict.values()]])

    return snaptime_dict

def snapnum_from_path(path):
    return int(path.split("snapshot_")[1].split(".")[0])
