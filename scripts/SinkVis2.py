#!/usr/bin/env python
"""
Usage:
SinkVis2.py <files> ... [options]

Options:
    -h --help                    Show this screen.
    --tasks=<task1,task2...>     List of types of the plots you want to make for each frame [default: SigmaGas]
    --rmax=<pc>                  Maximum radius of plot window; defaults to box size/10. Note that this is FOV/2 in radians if camera_dist is <inf
    --no_stars                   Hide sink particles
    --no_overwrite               Overwrite existing files if they already exist
    --backend=<b>                matplotlib vs PIL [default: PIL]
    --id_mask=<file>            .npy file containing the gas particle IDs to be plotted
    --dir=<x,y,z>                Coordinate direction to orient the image along - x, y, or z. It also accepts vector values [default: z] 
    --target_time=<f>            If set to nonzero, SinkVis will try to make a single image by interpolating from the available files [default: 0.0] 
    --limits=<min,max>           Dynamic range of surface density colormap
    --v_limits=<min,max>         Dynamic range of kinematic map in km/s
    --SHO_RGB_norm=<f>           Normalization constant for narrow band plot, set automatically by default. If a vector is provided, then each channel is normalized by the correponding component [default: 0.0]
    --camera_distance=<D>        Camera distance if perspective rendering is required [default: inf]
    --freeze_rotation=<num1,num2,...> Snapshot numbers at which to add a freeze-frame rotation [default: None]
    --cmap=<name>                Name of colormap to use [default: viridis]
    --cool_cmap=<name>           Name of colormap to use for plot_cool_map, defaults to same as cmap [default: magma]
    --interp_fac=<N>             Number of interpolating frames per snapshot [default: 1]
    --np=<N>                     Number of processors to run on [default: 1]
    --np_render=<N>              Number of openMP threads to run rendering on (-1 uses all available cores divided by --np) [default: 1]
    --res=<N>                    Image resolution [default: 1024]
    --center_on_star             Center image on the Nth most-massive sink particle
    --center_on_ID=<ID>          Center image on sink particle with specific ID, does not center if zero [default: 0]
    --plot_cool_map              Plots surface density+kinematics map that looks cool
    --outputfolder=<name>        Specifies the folder to save the images and movies to [default: .]
    --no_timestamp               Flag, if set no timestamp will be put on the images
    --no_size_scale              Flag, if set no size scale will be put on the images
    --rescale_hsml=<f>           Factor by which the smoothing lengths of the particles are rescaled [default: 1]
    --highlight_wind=<f>         Factor by which to increase wind particle masses if you want to highlight them [default: 1]
    --sparse_snaps               Flag, if enabled then corrections are applied to the interpolation algorithm to make the movies from sensitive maps (e.g. SHO narrowband) less flickery
    --equal_frame_times          Ensure frames in render sequence are equally spaced, even if snapshot times are not
    --outflow_only               Only show gas moving away from the nearest star
    --realstars                  Use realistic PSFs for stars
    --realstars_lum_exp=<f>      Exponent p such that luminosities are rescaled by raising to that power [default: 1.0]
    --realstars_max_lum=<f>      Maximum stellar luminosity in realistic PSF rendering [default: 1.0e3]

"""

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
    print(common_params["realstars_max_lum"], common_params["realstars_lum_exp"])
    del common_params["<files>"]
    for c, i in common_params.items():
        if not isinstance(i, str):
            continue
        if "," in i and c != "tasks":
            common_params[c] = np.array([float(k) for k in i.split(",")])
        elif i.replace(".", "").replace("e", "").isnumeric() or i == "inf":
            common_params[c] = float(i)

    print(common_params["realstars_max_lum"], common_params["realstars_lum_exp"])
    common_params.update(
        {
            "res": int(input["--res"]),
            "limits": (limits if arguments["--limits"] else None),
            "no_timestamp": input["--no_timestamp"],
            "threads": np_render,
            "outputfolder": input["--outputfolder"],
            "SHO_RGB_norm": SHO_RGB_norm,
            "cool_cmap": input["--cool_cmap"],
            "center_on_star": int(input["--center_on_star"]),
            "center_on_ID": int(input["--center_on_ID"]),
            "sparse_snaps": input["--sparse_snaps"],
            "backend": input["--backend"],
            "overwrite": overwrite,
            "id_mask": input["--id_mask"],
            "no_stars": input["--no_stars"],
        }
    )

    if direction == "x":
        common_params["camera_dir"] = np.array([1.0, 0, 0])
        common_params["camera_up"] = np.array([0, 0, 1.0])
    elif direction == "y":
        common_params["camera_dir"] = np.array([0, 1.0, 0])
        common_params["camera_up"] = np.array([0, 0, 1.0])

    if input["--rmax"] is None:
        common_params["rmax"] = None
    else:
        common_params["rmax"] = float(input["--rmax"])

    N_params = len(filenames) + (n_interp - 1) * (len(filenames) - 1)
    print(input["<files>"])
    snaptime_dict = get_snapshot_time_dict(input["<files>"])  # get times of snapshots
    snaptime_dict_inv = {v: k for v, k in zip(snaptime_dict.values(), snaptime_dict.keys())}
    snaptimes_orig = np.array(natsorted([snaptime_dict[snapnum_from_path(s)] for s in input["<files>"]]))
    #    print([snapnum_from_path(s) for s in input["<files>"]])
    if n_interp > 1:  # get times of frames if we're doing an interpolated movie
        frametimes = np.interp(
            np.arange(N_params) / (N_params - 1), np.linspace(0, 1, len(snaptimes_orig)), snaptimes_orig
        )
    else:
        frametimes = snaptimes_orig

    #    print(frametimes)
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
                for k in range(num_rotation_frames):  # do a pan
                    d = p[-1].copy()
                    d["index"] = snapnum * 10 + i % n_interp
                    d["pan"] = k * 360 / num_rotation_frames
                    d["tilt"] = 10 * np.sin(2 * np.pi * k / num_rotation_frames)  # add a bit of tilt for 3D look
                    p.append(d)

    #        print(i,d["Time"],d["index"])
    params = N_tasks * [p]  # one list of params per task
    return params


def main(input):
    tasks = input["--tasks"].split(",")
    tasks = [taskdict[t] for t in tasks]
    nproc = int(input["--np"])
    np_render = int(input["--np_render"])
    params = parse_inputs_to_jobparams(input)
    snaps = natsorted(input["<files>"])
    CrunchSnaps.DoTasksForSimulation(
        snaps, task_types=tasks, task_params=params, nproc=nproc, nthreads=np_render, id_mask=input["--id_mask"]
    )


if __name__ == "__main__":
    arguments = docopt(__doc__)
    main(arguments)
