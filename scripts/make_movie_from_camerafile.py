#!/usr/bin/env python
"""
Usage:
make_movie_from_camerafile.py <camerafile> <simdir> ... [options]

Options:
    -h --help              Show this screen.
    --map_type=<type>      Set the type of map to do, available options are SigmaGas,CoolMap and SHOMap [default: SigmaGas]
    --fresco_stars         Render stars with Fresco
    --extinct_stars        Calculate the extinction of stars to observers and attenuate their light, used only if --fresco_stars is set. Note: enabling this can make the calculation significantly slower
    --limits=<min,max>     Surface density limits
    --no_overwrite         Flag, if enabled existing figures are not overwritten
    --res=<N>              Resolution [default: 256]
    --np=<N>               Number of renders to do in parallel [default: 1]
    --np_render=<N>        Number of cores per process to run rundering calls on [default: 1]
    --cubemap              Render 6 faces of a cubemap surrounding the camera
    --no_timestamp         Don't add timestamp
    --no_size_scale        Don't draw size scale
    --SHO_RGB_norm=<f>     Normalization constant for narrow band plot, set automatically by default. If a vector is provided, then each channel is normalized by the correponding component [default: 0.0]
    --sparse_snaps         Flag, if enabled then corrections are applied to the interpolation algorithm to make the movies from sensitive maps (e.g. SHO narrowband) less flickery
    --id_mask=<f>          Path of the .npy file containing IDs of the particles to plot
    --backend=<b>          Backend [default: PIL]
    --sink_scale=<m>       Minimum sink mass to visualize in msun [default: 1e-3]
"""

from docopt import docopt
import numpy as np
from sys import argv
from natsort import natsorted
from CrunchSnaps import *
from glob import glob

options = docopt(__doc__)
cubemap = options["--cubemap"]
res = int(options["--res"])
nproc = int(options["--np"])
np_render = int(options["--np_render"])
if "," in options["--SHO_RGB_norm"]:  # normalization by channel
    SHO_RGB_norm = np.array([float(c) for c in options["--SHO_RGB_norm"].split(",")])
else:  # same normalization constant for each channel
    SHO_RGB_norm = float(options["--SHO_RGB_norm"])
if options["--no_overwrite"]:
    overwrite = False
else:
    overwrite = True
if options["--limits"]:
    limits = np.array([float(c) for c in options["--limits"].split(",")])
else:
    limits = None
id_mask = options["--id_mask"]
options["--sink_scale"] = float(options["--sink_scale"])

common_params = {
    "sink_scale": options["--sink_scale"],
    "backend": options["--backend"],
    "fresco_stars": options["--fresco_stars"],
    "res": res,
    "limits": limits,
    "no_timestamp": options["--no_timestamp"],
    "no_size_scale": options["--no_size_scale"],
    "threads": np_render,
    "SHO_RGB_norm": SHO_RGB_norm,
    "extinct_stars": options["--extinct_stars"],
    "overwrite": overwrite,
    "sparse_snaps": options["--sparse_snaps"],
}

camera_data = np.atleast_2d(np.loadtxt(options["<camerafile>"]))
sim_dir = options["<simdir>"][0]


# N = 1000
# time = np.linspace(0,5e-3,N)
# camera_dist = 3 + 30*np.exp(-time / 3e-4)
# angular_vel = 360 / (1e-3) * (3 / camera_dist)
# pan = cumtrapz(angular_vel,time,initial=0.0)
# tilt = 30*np.sin(2*np.pi * (time/1.3e-3))
# np.savetxt("camerafile.txt", np.c_[time, camera_dist, pan, tilt])

params = []


if options["--map_type"] == "SigmaGas":
    tasks = [SinkVisSigmaGas]
elif options["--map_type"] == "CoolMap":
    tasks = [SinkVisCoolMap]
elif options["--map_type"] == "SHOMap":
    tasks = [SinkVisNarrowbandComposite]
else:
    print("Map type %s not recognized, exiting..." % (options["--map_type"]))
    exit()


match camera_data.shape[1]:
    case 1:  # just times
        time = camera_data[:, 0]
        for i in range(len(time)):
            params.append({"Time": time[i]})
    case 4:  # columns will be time, distance, pan, tilt
        time, camera_dist, pan, tilt = camera_data.T
        for i in range(len(time)):
            params.append(
                {"Time": time[i], "camera_distance": camera_dist[i], "pan": pan[i], "tilt": tilt[i], "index": i}
            )
    case 5:  # time, position, rmax
        time = camera_data[:, 0]
        camera_pos = camera_data[:, 1:4]
        rmax = camera_data[:, 4]
        for i in range(len(time)):
            params.append({"Time": time[i], "center": camera_pos[i], "rmax": rmax[i], "index": i})
    case 7:  # columns will be time, distance, camera_pos, pan, tilt
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
                    "index": i,
                    "pan": pan[i],
                    "tilt": tilt[i],
                }
            )
    case 8:  # columns will be time, camera position, camera direction, camera distance
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
    case 11:  # full camera data: time, camera position, camera forward vector, camera up vector, camera distance
        time = camera_data[:, 0]
        camera_pos = camera_data[:, 1:4]
        camera_dir = camera_data[:, 4:7]
        camera_up = camera_data[:, 7:10]
        camera_dist = camera_data[:, 10]
    case _:
        raise ("camera file format not implemented :( do you have the right number of columns?")


for p in params:
    for k in common_params.keys():
        p[k] = common_params[k]

if cubemap:
    params_cubemap = []
    for p in params:
        params_cubemap += cubemapify(p)
    params = params_cubemap

snaps = natsorted(glob(sim_dir + "/snapshot*.hdf5"))

DoTasksForSimulation(
    snaps,
    task_types=tasks,
    task_params=[
        params,
    ],
    nproc=nproc,
    nthreads=np_render,
    id_mask=options["--id_mask"],
)
