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


def max_entropy_limits(vals, weights, n_bins=128, n_search=32, log_scale=False):
    """Find (lo, hi) that maximize mass-weighted Shannon entropy of the
    color histogram.

    The search runs on the *displayed* axis: when ``log_scale=True``, the
    objective is the entropy of the log10(value) histogram, and the
    returned limits are still in linear units.

    Ported from vizmo.field_ops.max_entropy_limits.
    """
    vals = np.asarray(vals).ravel()
    weights = np.asarray(weights).ravel()

    if log_scale:
        pos = vals > 0
        if not pos.any():
            return (float(vals.min()) if len(vals) else 0.0,
                    float(vals.max()) if len(vals) else 1.0)
        vals = np.log10(vals[pos])
        weights = weights[pos]

    if len(vals) > 100_000:
        rng = np.random.default_rng(0)
        idx = rng.choice(len(vals), size=100_000, replace=False)
        vals = vals[idx]
        weights = weights[idx]

    order = np.argsort(vals)
    sv = vals[order]
    sw = weights[order]
    cw = np.cumsum(sw)
    total = cw[-1]
    if total <= 0:
        lo_lin, hi_lin = float(sv[0]), float(sv[-1])
        if log_scale:
            return float(10 ** lo_lin), float(10 ** hi_lin)
        return lo_lin, hi_lin
    cw /= total

    lo_fracs = np.linspace(0.0, 0.45, n_search)
    hi_fracs = np.linspace(0.55, 1.0, n_search)
    lo_vals = np.interp(lo_fracs, cw, sv)
    hi_vals = np.interp(hi_fracs, cw, sv)

    t = np.arange(n_bins + 1, dtype=np.float64) / n_bins

    span = hi_vals[None, :, None] - lo_vals[:, None, None]
    edges = lo_vals[:, None, None] + t * span

    cdf_at_edges = np.interp(edges.ravel(), sv, cw).reshape(edges.shape)

    p = np.diff(cdf_at_edges, axis=2)
    p[:, :, 0] += cdf_at_edges[:, :, 0]
    p[:, :, -1] += 1.0 - cdf_at_edges[:, :, -1]

    with np.errstate(divide="ignore", invalid="ignore"):
        log_p = np.where(p > 0, np.log2(p), 0.0)
    H = -np.sum(p * log_p, axis=2)

    H[hi_vals[None, :] <= lo_vals[:, None]] = -np.inf

    i, j = np.unravel_index(H.argmax(), H.shape)
    lo_opt, hi_opt = float(lo_vals[i]), float(hi_vals[j])
    if log_scale:
        return float(10 ** lo_opt), float(10 ** hi_opt)
    return lo_opt, hi_opt


def max_entropy_limits_multi(map_arrays, n_bins=128, n_search=32, log_scale=False):
    """Find (lo, hi) that maximize the sum of per-frame Shannon entropies.

    Each frame's entropy is mass-weighted (pixel values used as weights),
    so limits are chosen to maximize information about where the mass is
    across the entire sequence.

    Args:
        map_arrays: list of 2D arrays (one per frame, e.g. surface density maps).
        n_bins: histogram bin count.
        n_search: grid resolution per axis (lo, hi).
        log_scale: when True, optimize in log10 space and return linear limits.
    """
    # Preprocess each frame: sort values, build mass-weighted CDFs
    frame_cdfs = []
    pooled_samples = []
    for arr in map_arrays:
        vals = np.asarray(arr).ravel()
        weights = vals.copy()
        if log_scale:
            pos = vals > 0
            if not pos.any():
                continue
            vals, weights = np.log10(vals[pos]), weights[pos]

        if len(vals) > 50_000:
            rng = np.random.default_rng(0)
            idx = rng.choice(len(vals), size=50_000, replace=False)
            vals, weights = vals[idx], weights[idx]

        order = np.argsort(vals)
        sv = vals[order]
        sw = weights[order]
        cw = np.cumsum(sw)
        total = cw[-1]
        if total <= 0:
            continue
        cw /= total
        frame_cdfs.append((sv, cw))
        pooled_samples.append(sv[:: max(1, len(sv) // 1000)])

    if not frame_cdfs:
        return (0.0, 1.0)

    # Build candidate lo/hi from pooled percentiles
    pooled = np.sort(np.concatenate(pooled_samples))
    pooled_cdf = np.linspace(0, 1, len(pooled))
    lo_fracs = np.linspace(0.0, 0.45, n_search)
    hi_fracs = np.linspace(0.55, 1.0, n_search)
    lo_vals = np.interp(lo_fracs, pooled_cdf, pooled)
    hi_vals = np.interp(hi_fracs, pooled_cdf, pooled)

    # Bin edges for all candidate pairs: shape (n_search, n_search, n_bins+1)
    t = np.arange(n_bins + 1, dtype=np.float64) / n_bins
    span = hi_vals[None, :, None] - lo_vals[:, None, None]
    edges = lo_vals[:, None, None] + t * span
    edges_flat = edges.ravel()

    # Sum per-frame entropies for each (lo, hi) pair
    H_total = np.zeros((n_search, n_search))
    for sv, cw in frame_cdfs:
        cdf_at_edges = np.interp(edges_flat, sv, cw).reshape(edges.shape)
        p = np.diff(cdf_at_edges, axis=2)
        p[:, :, 0] += cdf_at_edges[:, :, 0]
        p[:, :, -1] += 1.0 - cdf_at_edges[:, :, -1]
        with np.errstate(divide="ignore", invalid="ignore"):
            log_p = np.where(p > 0, np.log2(p), 0.0)
        H_total -= np.sum(p * log_p, axis=2)

    H_total[hi_vals[None, :] <= lo_vals[:, None]] = -np.inf

    i, j = np.unravel_index(H_total.argmax(), H_total.shape)
    lo_opt, hi_opt = float(lo_vals[i]), float(hi_vals[j])
    if log_scale:
        return float(10 ** lo_opt), float(10 ** hi_opt)
    return lo_opt, hi_opt


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
