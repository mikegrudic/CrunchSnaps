import numpy as np
from glob import glob
from math import erf
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


def coords_in_box_frame(x, boxsize, tol=0.01):
    """Whether coordinates are expressed in the periodic box frame [0, BoxSize).

    Periodic GIZMO runs write coordinates wrapped into [0, BoxSize), so
    box-frame operations (wrapping, periodic trees) are safe.  Non-periodic
    runs generally center the domain on the origin while still writing a
    nonzero BoxSize in the header; wrapping those into [0, BoxSize) tears the
    domain into pieces at the box faces.  The tolerance allows for particles
    that have drifted slightly outside a periodic box.
    """
    if boxsize is None or not np.isfinite(boxsize) or boxsize <= 0:
        return False
    return bool(x.min() > -tol * boxsize and x.max() < (1 + tol) * boxsize)


def sigma_quantile_limits(map_arrays, log_scale=False, nsigma=3.0, max_samples=50_000, max_total=4_000_000):
    """(lo, hi) at the +/-nsigma Gaussian-equivalent quantiles of the pooled pixel values.

    For nsigma=3 these are the 0.135th and 99.865th percentiles, i.e. the range a
    normal distribution keeps within +/-3 sigma.  Pixels are pooled over every frame
    so that a movie gets one stable color scale.  With ``log_scale=True`` the
    quantiles are taken in log10 space (matching how the values are displayed) and
    returned in linear units.  ``log_scale="auto"`` uses log10 space unless some
    frame has a negative pixel, in which case the limits come back symmetric about
    zero for a diverging colormap.

    Args:
        map_arrays: iterable of 2D arrays (one per frame, e.g. surface density maps).
            May be a generator, so that long sequences never hold every map at once.
        log_scale: True, False, or "auto" (decide from the sign of the data).  In
            log space, non-positive pixels are ignored.
        nsigma: half-width of the kept range, in Gaussian sigmas.
        max_samples: per-frame subsample cap, to bound the cost on large maps.
        max_total: cap on pooled samples; once reached, the pool is thinned and
            later frames are sampled more sparsely, so memory stays bounded no
            matter how many frames there are.
    """
    auto = str(log_scale).lower() == "auto"
    rng = np.random.default_rng(0)
    samples = []
    n_pooled = 0
    for arr in map_arrays:
        vals = np.asarray(arr).ravel()
        vals = vals[np.isfinite(vals)]
        if log_scale is True:
            vals = np.log10(vals[vals > 0])
        # in auto mode the pool stays linear: the sign of the whole sequence is
        # only known once every frame has been seen
        if not len(vals):
            continue
        if len(vals) > max_samples:
            vals = vals[rng.choice(len(vals), size=max_samples, replace=False)]
        samples.append(vals)
        n_pooled += len(vals)
        if n_pooled > max_total:  # thin the pool and take less from each later frame
            samples = [s[::2] for s in samples]
            n_pooled = sum(len(s) for s in samples)
            max_samples = max(1000, max_samples // 2)

    if not samples:
        return (0.0, 1.0)

    pooled = np.concatenate(samples)
    q_hi = 0.5 * (1.0 + erf(nsigma / np.sqrt(2)))  # 0.99865 for 3 sigma
    if auto:
        if pooled.min() < 0:  # signed field: symmetric limits, linear scale
            lo, hi = np.quantile(pooled, [1.0 - q_hi, q_hi])
            vext = max(abs(lo), abs(hi))
            return -float(vext), float(vext)
        # non-negative: log limits from the pixels with signal, so that a mostly
        # empty frame does not drag the low end down to zero
        pooled = np.log10(pooled[pooled > 0])
        if not pooled.size:
            return (0.0, 1.0)
    lo, hi = np.quantile(pooled, [1.0 - q_hi, q_hi])
    if log_scale is True or auto:
        lo, hi = 10.0**lo, 10.0**hi
    return float(lo), float(hi)


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
        snapnums = snapnums.astype(int)  # loadtxt returns float64; keep snapnums as ints so index stays int
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
        print(f"Exception {e} when reading {path}")
        return 0
