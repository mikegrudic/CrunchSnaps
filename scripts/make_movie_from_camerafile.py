#!/usr/bin/env python
"""
Usage:
make_movie_from_camerafile.py <camerafile> <simdir> ... [options]

Options:
    -h --help              Show this screen.
    --map_type=<type>      Map type: SigmaGas, CoolMap, SHOMap [default: SigmaGas]
    --realstars            Render stars with realistic PSFs
    --limits=<min,max>     Surface density limits
    --no_overwrite         Skip rendering if the output file already exists
    --res=<N>              Resolution [default: 256]
    --np=<N>               Number of renders to do in parallel [default: 1]
    --np_render=<N>        Number of cores per process for rendering [default: 1]
    --cubemap              Render 6 faces of a cubemap surrounding the camera
    --no_timestamp         Don't add timestamp
    --no_size_scale        Don't draw size scale
    --no_colorbar          Don't draw colorbar
    --cmap=<name>          Colormap [default: viridis]
    --SHO_RGB_norm=<f>     Normalization for narrow band plot [default: 0.0]
    --sparse_snaps         Corrections for flickery interpolation on sparse snapshots
    --id_mask=<f>          Path of the .npy file containing IDs of the particles to plot
    --backend=<b>          Backend [default: PIL]
    --outputfolder=<name>  Output directory [default: .]
    --make_movie           Stitch rendered frames into an mp4 movie
    --fps=<N>              Movie frame rate [default: 24]
"""

from docopt import docopt
import numpy as np
from natsort import natsorted
from CrunchSnaps import *
from CrunchSnaps.sinkvis2 import make_movie
from glob import glob

options = docopt(__doc__)
cubemap = options["--cubemap"]
res = int(options["--res"])
nproc = int(options["--np"])
np_render = int(options["--np_render"])
if "," in options["--SHO_RGB_norm"]:
    SHO_RGB_norm = np.array([float(c) for c in options["--SHO_RGB_norm"].split(",")])
else:
    SHO_RGB_norm = float(options["--SHO_RGB_norm"])
overwrite = not options["--no_overwrite"]
if options["--limits"]:
    limits = np.array([float(c) for c in options["--limits"].split(",")])
else:
    limits = None
outputfolder = options["--outputfolder"]

common_params = {
    "backend": options["--backend"],
    "realstars": options["--realstars"],
    "res": res,
    "limits": limits,
    "no_timestamp": options["--no_timestamp"],
    "no_size_scale": options["--no_size_scale"],
    "no_colorbar": options["--no_colorbar"],
    "cmap": options["--cmap"],
    "threads": np_render,
    "SHO_RGB_norm": SHO_RGB_norm,
    "overwrite": overwrite,
    "sparse_snaps": options["--sparse_snaps"],
    "outputfolder": outputfolder,
}

map_type_to_task = {
    "SigmaGas": SinkVisSigmaGas,
    "CoolMap": SinkVisCoolMap,
    "SHOMap": SinkVisNarrowbandComposite,
}

map_type_to_prefix = {
    "SigmaGas": "SurfaceDensity",
    "CoolMap": "CoolMap",
    "SHOMap": "NarrowbandComposite",
}

map_type = options["--map_type"]
if map_type not in map_type_to_task:
    raise ValueError(f"Unknown map type '{map_type}'. Choose from: {', '.join(map_type_to_task)}")
tasks = [map_type_to_task[map_type]]

camera_data = np.atleast_2d(np.loadtxt(options["<camerafile>"]))

params = []

match camera_data.shape[1]:
    case 1:  # just times
        time = camera_data[:, 0]
        for i in range(len(time)):
            params.append({"Time": time[i], "index": i})
    case 4:  # time, distance, pan, tilt
        time, camera_dist, pan, tilt = camera_data.T
        for i in range(len(time)):
            params.append(
                {"Time": time[i], "camera_distance": camera_dist[i], "pan": pan[i], "tilt": tilt[i], "index": i}
            )
    case 5:  # time, position(3), rmax
        time = camera_data[:, 0]
        camera_pos = camera_data[:, 1:4]
        rmax = camera_data[:, 4]
        for i in range(len(time)):
            params.append({"Time": time[i], "center": camera_pos[i], "rmax": rmax[i], "index": i})
    case 7:  # time, distance, position(3), pan, tilt
        time = camera_data[:, 0]
        camera_dist = camera_data[:, 1]
        camera_pos = camera_data[:, 2:5]
        pan, tilt = camera_data[:, -2:].T
        for i in range(len(time)):
            params.append(
                {
                    "Time": time[i],
                    "center": camera_pos[i],
                    "camera_distance": camera_dist[i],
                    "pan": pan[i],
                    "tilt": tilt[i],
                    "index": i,
                }
            )
    case 8:  # time, position(3), direction(3), distance
        time = camera_data[:, 0]
        camera_pos = camera_data[:, 1:4]
        camera_dir = camera_data[:, 4:7]
        camera_dist = camera_data[:, 7]
        for i in range(len(time)):
            params.append(
                {
                    "Time": time[i],
                    "center": camera_pos[i],
                    "camera_dir": camera_dir[i],
                    "camera_distance": camera_dist[i],
                    "index": i,
                }
            )
    case 11:  # time, position(3), direction(3), up(3), distance
        time = camera_data[:, 0]
        camera_pos = camera_data[:, 1:4]
        camera_dir = camera_data[:, 4:7]
        camera_up = camera_data[:, 7:10]
        camera_dist = camera_data[:, 10]
        for i in range(len(time)):
            params.append(
                {
                    "Time": time[i],
                    "center": camera_pos[i],
                    "camera_dir": camera_dir[i],
                    "camera_up": camera_up[i],
                    "camera_distance": camera_dist[i],
                    "index": i,
                }
            )
    case n:
        raise ValueError(f"Camera file has {n} columns; expected 1, 4, 5, 7, 8, or 11")


for p in params:
    p.update(common_params)

if cubemap:
    params = [cp for p in params for cp in cubemapify(p)]

sim_dir = options["<simdir>"][0]
snaps = natsorted(glob(sim_dir + "/snapshot*.hdf5"))

DoTasksForSimulation(
    snaps,
    task_types=tasks,
    task_params=[params],
    nproc=nproc,
    nthreads=np_render,
    id_mask=options["--id_mask"],
)

if options["--make_movie"]:
    fps = int(options["--fps"])
    prefix = map_type_to_prefix[map_type]
    make_movie(outputfolder, prefix, fps)
