import numpy as np
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


def _NormalizeVector_numpy(v):
    v /= np.sqrt(np.dot(v, v))


def _NearestImage_numpy(dx, boxsize):
    half = boxsize / 2
    return np.where(np.abs(dx) > half, -np.copysign(boxsize - np.abs(dx), dx), dx)


def _init_numba_funcs():
    """Lazy-initialize numba-accelerated versions of NormalizeVector and NearestImage."""
    global NormalizeVector, NearestImage
    try:
        from numba import njit, vectorize
        from numpy import abs, copysign, sqrt

        @njit(cache=True)
        def _NormalizeVector_numba(v):
            norm = 0
            for k in range(v.shape[0]):
                norm += v[k] * v[k]
            norm = 1.0 / sqrt(norm)
            for k in range(v.shape[0]):
                v[k] *= norm

        @vectorize(cache=True)
        def _NearestImage_numba(dx, boxsize):
            if abs(dx) > boxsize / 2:
                return -copysign(boxsize - abs(dx), dx)
            else:
                return dx

        NormalizeVector = _NormalizeVector_numba
        NearestImage = _NearestImage_numba
    except ImportError:
        NormalizeVector = _NormalizeVector_numpy
        NearestImage = _NearestImage_numpy


def NormalizeVector(v):
    _init_numba_funcs()
    NormalizeVector(v)


def NearestImage(dx, boxsize):
    _init_numba_funcs()
    return NearestImage(dx, boxsize)


def get_snapshot_time_dict(snaps, save_to_file=True):
    snaps = natsorted(snaps)
    all_snaps = natsorted(
        glob(snaps[0].split("snapshot")[0] + "snapshot*.hdf5")
    )  # look for other snapshots in same directory
    if (
        len(snaps) == 1
    ):  # if we only have one snapshot, keep it simple and just open the file. otherwise we will do some fancy stuff to avoid opening multiple files for the timeline
        return {snapnum_from_path(snaps[0]): h5py.File(snaps[0], "r")["Header"].attrs["Time"]}

    # don't yet know what the snapshot times are - get the snapshot times in a prepass
    snaptimes = []
    snapnums = []

    snaptimes_path = abspath(snaps[0]).split("/snapshot")[0] + "/.snapshot_times"

    if exists(snaptimes_path):
        snapnums, snaptimes = np.atleast_2d(np.loadtxt(snaptimes_path).T)
    snaptime_dict = dict(zip(snapnums, snaptimes))

    do_snapshot_pass = False
    if np.any(
        [not (snapnum_from_path(s) in snaptime_dict.keys()) for s in snaps]
    ):  # check if we have a snapshot missing from the dictionary, if so we must do a pass
        print("Sinkvis2 getting snapshot times...")
        for s in all_snaps:
            print(s)
            try:
                with h5py.File(s, "r") as F:
                    snaptime_dict[snapnum_from_path(s)] = F["Header"].attrs["Time"]
            except:
                pass
        if save_to_file:
            np.savetxt(snaptimes_path, np.c_[[k for k in snaptime_dict.keys()], [v for v in snaptime_dict.values()]])

    return snaptime_dict


def snapnum_from_path(path):
    try:
        return int(path.split("snapshot_")[1].split("_")[0].split(".")[0])
    except Exception as e:
        print(f"Exception {f} when reading {path}")
        return 0
