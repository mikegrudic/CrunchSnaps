#!/usr/bin/env python
"""
Usage:
SinkVis.py <camerafile> <simdir> ... [options]

Options:
    -h --help              Show this screen.
    --fresco_stars         Render stars with Fresco
    --limits=<min,max>     Surface density limits  [default: 1,3e3]
    --res=<N>              Resolution [default: 256]
    --np=<N>               Number of renders to do in parallel [default: 1]
    --np_render=<N>        Number of cores per process to run rundering calls on [default: 1]
    --cubemap              Render 6 faces of a cubemap surrounding the camera
    --no_timestamp
"""

from docopt import docopt
import numpy as np
from sys import argv
from natsort import natsorted
from CrunchSnaps import *
from glob import glob
from scipy.integrate import cumtrapz

options = docopt(__doc__)
cubemap = options["--cubemap"]
res = int(options["--res"])
nproc = int(options["--np"])
np_render = int(options["--np_render"])

limits = np.array([float(c) for c in options["--limits"].split(',')])

common_params = {"fresco_stars": options["--fresco_stars"], "res": res, "limits": limits, "no_timestamp": options["--no_timestamp"], "threads": np_render}

camera_data = np.loadtxt(options["<camerafile>"])
sim_dir = options["<simdir>"][0] 


#N = 1000
#time = np.linspace(0,5e-3,N)
#camera_dist = 3 + 30*np.exp(-time / 3e-4)
#angular_vel = 360 / (1e-3) * (3 / camera_dist)
#pan = cumtrapz(angular_vel,time,initial=0.0)
#tilt = 30*np.sin(2*np.pi * (time/1.3e-3))
#np.savetxt("camerafile.txt", np.c_[time, camera_dist, pan, tilt])


params = []
tasks = [SinkVisCoolMap,]

if camera_data.shape[1] == 4: # columns will be time, distance, pan, tilt
    time, camera_dist, pan, tilt = camera_data.T
    for i in range(len(time)):
        params.append({"Time": time[i], "camera_distance": camera_dist[i], "pan": pan[i], "tilt": tilt[i], "index": i})
elif camera_data.shape[1] == 7: # columns will be time, camera_pos, distance, pan, tilt
    time = camera_data[:,0]
    camera_dist = camera_data[:,1]
    camera_pos = camera_data[:,1:4]
    pan, tilt = camera_data[:,-2:].T
elif camera_data.shape[1] == 8: # columns will be time, camera position, camera direction, camera distance
    time = camera_data[:,0]
    camera_pos = camera_data[:,1:4]
    camera_dir = camera_data[:,4:7]
    camera_dist = camera_data[:,7]
    for i in range(len(time)):
        params.append({"Time": time[i], "center": camera_pos[i], "camera_dir": camera_dir[i], "camera_distance": camera_dist[i], "index": i})
elif camera_data.shape[1] == 11:# full camera data: time, camera position, camera forward vector, camera up vector, camera distance
    time = camera_data[:,0]
    camera_pos = camera_data[:,1:4]
    camera_dir = camera_data[:,4:7]
    camera_up = camera_data[:,7:10]
    camera_dist = camera_data[:,10]
else: 
    raise("camera file format not implemented :( do you have the right number of columns?")


for p in params:
    for k in common_params.keys():
        p[k] = common_params[k]

if cubemap:
    params_cubemap = []
    for p in params:
        params_cubemap += cubemapify(p)
    params = params_cubemap

snaps = natsorted(glob(sim_dir + "/snapshot*.hdf5"))

DoTasksForSimulation(snaps, tasks=tasks, task_params=[params],nproc=nproc,nthreads=np_render)
