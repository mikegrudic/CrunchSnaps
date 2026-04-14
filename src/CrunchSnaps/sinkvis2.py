"""
Usage:
SinkVis2 <args> ... [options]

Options:
    -h --help                    Show this screen.
    --tasks=<task1,task2...>     Comma-separated list of tasks (alternative to positional args)
                                 Built-in: SigmaGas, HubbleSHO, CoolMap
                                 Custom: SurfaceDensity(expr), Projection(expr),
                                         ProjectedAverage(expr), Slice(expr)
    --rmax=<pc>                  Half-width of the field of view; defaults to mass-weighted spatial extent of gas. FOV/2 in radians if camera_dist is <inf
    --no_stars                   Hide sink particles
    --no_overwrite               Skip rendering if the output file already exists
    --backend=<b>                matplotlib vs PIL [default: PIL]
    --id_mask=<file>            .npy file containing the gas particle IDs to be plotted
    --dir=<x,y,z>                Coordinate direction to orient the image along - x, y, or z. It also accepts vector values [default: z]
    --target_time=<f>            Render a single image at this simulation time, interpolating between bounding snapshots
    --limits=<min,max>           Dynamic range of surface density colormap
    --v_limits=<min,max>         Dynamic range of kinematic map in km/s
    --SHO_RGB_norm=<f>           Normalization constant for narrow band plot, set automatically by default. If a vector is provided, then each channel is normalized by the correponding component [default: 0.0]
    --camera_distance=<D>        Camera distance if perspective rendering is required [default: inf]
    --camera_up=<x,y,z>          Camera up vector [default: None]
    --pan=<deg>                  Pan angle in degrees (rotation about Y axis) [default: 0]
    --tilt=<deg>                 Tilt angle in degrees (rotation about X axis) [default: 0]
    --freeze_rotation=<num1,num2,...> Snapshot numbers at which to add a freeze-frame rotation [default: None]
    --cmap=<name>                Name of colormap to use [default: viridis]
    --cool_cmap=<name>           Name of colormap to use for plot_cool_map, defaults to same as cmap [default: magma]
    --interp_fac=<N>             Number of interpolating frames per snapshot [default: 1]
    --np=<N>                     Number of processors to run on [default: 1]
    --np_render=<N>              Number of threads for rendering (-1 uses all available cores divided by --np) [default: -1]
    --res=<N>                    Image resolution [default: 1024]
    --center=<s>                 Center of the image. Argument can include 3D comma-separated coordinates, a particle
                                 ID, 'densest' to center on the densest gas, 'median' to center on the median gas cell
                                 position, or 'massive' to center on the most-massive star [default: None]
    --outputfolder=<name>        Specifies the folder to save the images and movies to [default: .]
    --no_timestamp               Flag, if set no timestamp will be put on the images
    --no_size_scale              Flag, if set no size scale will be put on the images
    --no_colorbar                Flag, if set no colorbar will be put on surface density images
    --rescale_hsml=<f>           Factor by which the smoothing lengths of the particles are rescaled [default: 1]
    --sparse_snaps               Flag, if enabled then corrections are applied to the interpolation algorithm to make the movies from sensitive maps (e.g. SHO narrowband) less flickery
    --equal_frame_times          Ensure frames in render sequence are equally spaced, even if snapshot times are not
    --outflow_only               Only show gas moving away from the nearest star
    --realstars                  Use realistic PSFs for stars
    --realstars_lum_exp=<f>      Exponent p such that luminosities are rescaled by raising to that power [default: 1.0]
    --realstars_max_lum=<f>      Maximum stellar luminosity in realistic PSF rendering [default: 1.0e3]
    --realstars_opacity=<f>      Opacity scaling factor for realstars [default: 1.0]
    --supersample=<N>            Anti-aliasing supersampling factor for Slice renders [default: 2]
    --unit_length=<f>            Override UnitLength_In_CGS from snapshot header
    --unit_mass=<f>              Override UnitMass_In_CGS from snapshot header
    --unit_velocity=<f>          Override UnitVelocity_In_CGS from snapshot header
    --unit_B=<f>                 Magnetic field unit in Gauss (default: 1, i.e. B stored in Gauss)
    --make_movie                 Stitch rendered frames into a 24fps mp4 movie using ffmpeg
    --fps=<N>                    Frames per second for movie [default: 24]
"""

import os
import sys

if not os.getenv("PYTHONHASHSEED"):
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

from docopt import docopt
from natsort import natsorted
import numpy as np
import CrunchSnaps
from CrunchSnaps.misc_functions import *

taskdict = {
    "SigmaGas": CrunchSnaps.SinkVisSigmaGas,
    "HubbleSHO": CrunchSnaps.SinkVisNarrowbandComposite,
    "CoolMap": CrunchSnaps.SinkVisCoolMap,
}

# map task names to their output filename prefixes
task_prefix = {
    "SigmaGas": "SurfaceDensity",
    "HubbleSHO": "NarrowbandComposite",
    "CoolMap": "CoolMap",
}

import re


def _resolve_task(spec):
    """Resolve a task spec string to (task_class, extra_params, movie_prefix).

    Handles both built-in names ('SigmaGas') and custom field expressions
    like 'SurfaceDensity(Masses*InternalEnergy)'.
    """
    if spec in taskdict:
        return taskdict[spec], {}, task_prefix[spec]

    parsed = CrunchSnaps.parse_custom_task(spec)
    if parsed is not None:
        render_mode, field_expr = parsed
        safe_expr = re.sub(r"[^\w]", "_", field_expr)
        prefix = f"{render_mode}_{safe_expr}"
        extra = {"_render_mode": render_mode, "_field_expr": field_expr}
        return CrunchSnaps.SinkVisCustomField, extra, prefix

    raise ValueError(
        f"Unknown task '{spec}'. Use a built-in name ({', '.join(taskdict)}) "
        f"or a custom expression like SurfaceDensity(Masses*Temperature)"
    )


def parse_inputs_to_jobparams(input):  # parse input parameters to generate a list of job parameters
    arguments = input
    filenames = natsorted(arguments["<files>"])
    np_render = int(arguments["--np_render"])
    n_interp = int(arguments["--interp_fac"])
    tasks = arguments["--tasks"].split(",")
    N_tasks = len(tasks)
    direction = arguments["--dir"]
    if arguments["--no_overwrite"]:
        overwrite = False
    else:
        overwrite = True
    equal_frame_times = arguments["--equal_frame_times"]

    if "," in arguments["--SHO_RGB_norm"]:  # normalization by channel
        SHO_RGB_norm = np.array([float(c) for c in arguments["--SHO_RGB_norm"].split(",")])
    else:  # same normalization constant for each channel
        SHO_RGB_norm = float(arguments["--SHO_RGB_norm"])

    if arguments["--limits"]:
        limits = np.array([float(c) for c in arguments["--limits"].split(",")])

    # parameters that every single task will have in common
    common_params = {i.replace("--", ""): input[i] for i in input}
    for k in list(common_params):
        if k.startswith("<"):
            del common_params[k]
    for c, i in common_params.items():
        if not isinstance(i, str):
            continue
        if "," in i and c != "tasks":
            common_params[c] = np.array([float(k) for k in i.split(",")])
        elif i.replace(".", "").replace("e", "").isnumeric() or i == "inf":
            common_params[c] = float(i)

    common_params.update(
        {
            "res": int(input["--res"]),
            "limits": (limits if arguments["--limits"] else None),
            "no_timestamp": input["--no_timestamp"],
            "threads": np_render,
            "outputfolder": input["--outputfolder"],
            "SHO_RGB_norm": SHO_RGB_norm,
            "cool_cmap": input["--cool_cmap"],
            "center": input["--center"],
            "sparse_snaps": input["--sparse_snaps"],
            "backend": input["--backend"],
            "overwrite": overwrite,
            "id_mask": input["--id_mask"],
            "no_stars": input["--no_stars"],
            "no_colorbar": input["--no_colorbar"],
        }
    )

    # Unit overrides — injected into snapshot headers at load time
    unit_overrides = {}
    if input["--unit_length"]:
        unit_overrides["UnitLength_In_CGS"] = float(input["--unit_length"])
    if input["--unit_mass"]:
        unit_overrides["UnitMass_In_CGS"] = float(input["--unit_mass"])
    if input["--unit_velocity"]:
        unit_overrides["UnitVelocity_In_CGS"] = float(input["--unit_velocity"])
    if input["--unit_B"]:
        unit_overrides["UnitB_In_Gauss"] = float(input["--unit_B"])
    common_params["_unit_overrides"] = unit_overrides

    common_params["pan"] = float(input["--pan"])
    common_params["tilt"] = float(input["--tilt"])

    if direction == "x":
        common_params["camera_dir"] = np.array([1.0, 0, 0])
        common_params["camera_up"] = np.array([0, 0, 1.0])
    elif direction == "y":
        common_params["camera_dir"] = np.array([0, 1.0, 0])
        common_params["camera_up"] = np.array([0, 0, 1.0])
    elif isinstance(direction, np.ndarray):
        common_params["camera_dir"] = direction

    if input["--camera_up"] is not None and input["--camera_up"] != "None":
        common_params["camera_up"] = np.array([float(x) for x in input["--camera_up"].split(",")])

    if input["--rmax"] is None:
        common_params["rmax"] = None
    else:
        common_params["rmax"] = float(input["--rmax"])

    print(input["<files>"])
    snaptime_dict = get_snapshot_time_dict(input["<files>"])  # get times of snapshots
    snaptime_dict_inv = {v: k for v, k in zip(snaptime_dict.values(), snaptime_dict.keys())}
    snaptimes_orig = np.array(natsorted([snaptime_dict[snapnum_from_path(s)] for s in input["<files>"]]))

    target_time = input["--target_time"]
    if target_time is not None:
        # Single image at a specific simulation time
        target_time = float(target_time)
        # Find the bounding snapshot index for labeling
        idx = np.searchsorted(snaptimes_orig, target_time, side="right")
        idx = np.clip(idx, 0, len(snaptimes_orig) - 1)
        snapnum = snaptime_dict_inv[snaptimes_orig[idx]]
        d = common_params.copy()
        d["Time"] = target_time
        d["index"] = round(target_time / max(snaptimes_orig[1] - snaptimes_orig[0], 1e-30))
        p = [d]
    else:
        N_params = len(filenames) + (n_interp - 1) * (len(filenames) - 1)
        if n_interp > 1:  # get times of frames if we're doing an interpolated movie
            frametimes = np.interp(
                np.arange(N_params) / (N_params - 1), np.linspace(0, 1, len(snaptimes_orig)), snaptimes_orig
            )
        else:
            frametimes = snaptimes_orig

        if equal_frame_times:
            frametimes = np.linspace(frametimes.min(), frametimes.max(), N_params)

        p = []
        for i in range(N_params):
            d = common_params.copy()
            d["Time"] = frametimes[i]
            snapnum = snaptime_dict_inv[snaptimes_orig[i // n_interp]]
            d["index"] = snapnum * 10 + i % n_interp
            p.append(d.copy())

            if i % n_interp == 0 and (input["--freeze_rotation"] is not None):
                num_rotation_frames = 720
                if snapnum in [int(f) for f in input["--freeze_rotation"].split(",")]:  # add a rotation freeze
                    for k in range(1, num_rotation_frames):  # do a pan, skip k=0 (duplicate of original)
                        d = p[-1].copy()
                        d["index"] = snapnum * 10 + i % n_interp
                        d["pan"] = k * 360 / num_rotation_frames
                        d["tilt"] = 10 * np.sin(2 * np.pi * k / num_rotation_frames)  # add a bit of tilt for 3D look
                        p.append(d)

    params = [[d.copy() for d in p] for _ in range(N_tasks)]  # independent copy per task
    return params


def make_movie(outputfolder, prefix, fps):
    """Stitch all frames matching prefix in outputfolder into an mp4 movie."""
    import subprocess
    from glob import glob

    frames = natsorted(glob(os.path.join(outputfolder, prefix + "_*.png")))
    # exclude incomplete files
    frames = [f for f in frames if "_incomplete" not in f]
    if len(frames) < 2:
        print(f"Skipping movie for {prefix}: only {len(frames)} frame(s) found")
        return

    # write a temporary file list for ffmpeg concat demuxer
    listfile = os.path.join(outputfolder, f".{prefix}_framelist.txt")
    with open(listfile, "w") as f:
        for frame in frames:
            f.write(f"file {os.path.abspath(frame)}\n")
            f.write(f"duration {1/fps}\n")

    movie_path = os.path.join(outputfolder, f"{prefix}.mp4")
    cmd = [
        "ffmpeg", "-y", "-f", "concat", "-safe", "0", "-i", listfile,
        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-c:v", "libx264",  "-pix_fmt", "yuv420p",
        "-r", str(fps),
        movie_path,
    ]
    print(f"Encoding {movie_path} from {len(frames)} frames at {fps} fps...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    os.remove(listfile)
    if result.returncode != 0:
        print(f"ffmpeg error:\n{result.stderr}")
    else:
        print(f"Wrote {movie_path}")


def _compute_global_limits(outputfolder):
    """Load all cached sigma_gas maps and compute global max-entropy limits."""
    from glob import glob

    map_files = natsorted(glob(os.path.join(outputfolder, ".maps", "sigma_gas_*.npz")))
    if not map_files:
        return None
    maps = [np.load(f)["sigma_gas"] for f in map_files]
    print(f"Computing global max-entropy limits from {len(maps)} cached maps...")
    limits = np.array(max_entropy_limits_multi(maps, log_scale=True))
    print(f"Global limits: [{limits[0]:.3g}, {limits[1]:.3g}]")
    return limits


def main(input):
    task_specs = input["--tasks"].split(",")
    resolved = [_resolve_task(spec) for spec in task_specs]
    tasks = [r[0] for r in resolved]
    extra_params = [r[1] for r in resolved]
    movie_prefixes = [r[2] for r in resolved]

    nproc = int(input["--np"])
    np_render = int(input["--np_render"])
    params = parse_inputs_to_jobparams(input)
    snaps = natsorted(input["<files>"])
    user_set_limits = input["--limits"] is not None
    outputfolder = input["--outputfolder"]

    # Inject per-task extra params (e.g. _render_mode, _field_expr for custom tasks)
    for i, ep in enumerate(extra_params):
        if ep:
            for p in params[i]:
                p.update(ep)

    # First pass: render maps (and images with per-frame limits)
    CrunchSnaps.DoTasksForSimulation(
        snaps, task_types=tasks, task_params=params, nproc=nproc, nthreads=np_render, id_mask=input["--id_mask"]
    )

    if input["--make_movie"]:
        # If user didn't specify limits, re-render with global max-entropy limits
        if not user_set_limits:
            global_limits = _compute_global_limits(outputfolder)
            if global_limits is not None:
                for task_params_list in params:
                    for p in task_params_list:
                        p["limits"] = global_limits
                        p["overwrite"] = True
                # Second pass: maps are cached, only images are re-rendered
                CrunchSnaps.DoTasksForSimulation(
                    snaps, task_types=tasks, task_params=params,
                    nproc=nproc, nthreads=np_render, id_mask=input["--id_mask"],
                )

        fps = int(input["--fps"])
        for prefix in movie_prefixes:
            make_movie(outputfolder, prefix, fps)


def main_cli():
    arguments = docopt(__doc__)
    # Split positional args into snapshot files and task specs
    positional = arguments["<args>"]
    files = [a for a in positional if a.endswith(".hdf5")]
    pos_tasks = [a for a in positional if not a.endswith(".hdf5")]
    arguments["<files>"] = files
    if pos_tasks:
        arguments["--tasks"] = ",".join(pos_tasks)
    elif not arguments["--tasks"]:
        arguments["--tasks"] = "SigmaGas"
    main(arguments)


if __name__ == "__main__":
    main_cli()
