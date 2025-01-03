from glob import glob
from CrunchSnaps import *
from natsort import natsorted
import numpy as np
from sys import argv

tasks = [
    SinkVisCoolMap,
]  # select the tasks that you want to run
snaps = natsorted(
    glob(argv[1] + "/snapshot*.hdf5")
)  # get the list of snapshots to run on (should sort chronologically for best access pattern)

# Now here we generate a list of parameters for each frame we want to make, this will be done "by hand" and all unintialized parameters will adopt default values
Nangle = 1080
params = []
tilt = 0
pan = 0
pan_seconds_per_rotation = 30
fps = 24
pan_deg_per_frame = 360 / (fps * pan_seconds_per_rotation)
res = 512  # 1920 #3840
frame = 0
Tmax = 7.65e-3
center = np.float64([50, 50, 50])
tot_frames = Nangle + 10 * fps + 360 + pan_seconds_per_rotation * fps

centers = center + np.zeros(
    (tot_frames, 3)
)  # np.c_[8*np.cos(3*np.pi*np.arange(tot_frames)/tot_frames), np.zeros(tot_frames),8*np.sin(6*np.pi*np.arange(tot_frames)/tot_frames)]

for i in list(range(Nangle)):
    tilt = 25 * np.sin(2 * np.pi * i / Nangle)
    params.append(
        {
            "Time": 5e-3 * i / Nangle,
            "pan": pan,
            "tilt": tilt,
            "camera_distance": 8 * (1 + 4 * np.exp(-5 * i / Nangle)),
            "filename": "cool__%s.png" % str(frame).zfill(4),
            "res": res,
            "fresco_stars": True,
            #                   "center": centers[frame],
            "no_timestamp": True,
        }
    )
    frame += 1
    if i > fps * 5:
        pan += pan_deg_per_frame
for i in range(10 * fps):
    params.append(
        {
            "Time": 5e-3,
            "pan": pan,
            "tilt": tilt,
            "camera_distance": 8 * np.exp(0.5 * (-1 + np.cos(2 * np.pi * i / (10 * fps)))),
            "filename": "cool__%s.png" % str(frame).zfill(4),
            "res": res,
            "fresco_stars": True,
            "no_timestamp": True,
        }
    )
    frame += 1
    pan += pan_deg_per_frame

i0 = i
for i in range(i0 + 1, i0 + 360):
    params.append(
        {
            "Time": 5e-3 + (Tmax - 5e-3) * (i - i0) / 360,
            "pan": pan,
            "camera_distance": 8 + (i - i0) * 10 / 360,
            "filename": "cool__%s.png" % str(frame).zfill(4),
            "res": res,
            "fresco_stars": True,
            "no_timestamp": True,
        }
    )
    frame += 1

for i in range(pan_seconds_per_rotation * fps):
    pan += pan_deg_per_frame
    tilt = np.sin(3 * np.pi * i / (pan_seconds_per_rotation * fps))
    params.append(
        {
            "Time": Tmax,
            "pan": pan,
            "camera_distance": 18 * np.exp(np.sin(1.5 * np.pi * i / (pan_seconds_per_rotation * fps))),
            "filename": "cool__%s.png" % str(frame).zfill(4),
            "res": res,
            "tilt": tilt,
            "fresco_stars": True,
            "no_timestamp": True,
        }
    )
    frame += 1
print(len(params))
params = [
    params,
]  # format as a list of lists with one entry per task
DoTasksForSimulation(snaps=snaps, tasks=tasks, task_params=params, nproc=16, nthreads=1)
