from glob import glob
from CrunchSnaps import *
from tasks import *
from natsort import natsorted
import numpy as np
from sys import argv

tasks = [SinkVisCoolMap,] # select the tasks that you want to run
snaps = natsorted(glob(argv[1]+"/snapshot*.hdf5")) # get the list of snapshots to run on (should sort chronologically for best access pattern)

# Now here we generate a list of parameters for each frame we want to make, this will be done "by hand" and all unintialized parameters will adopt default values
N = 256
dists = np.logspace(np.log10(30),np.log10(1),N)
params = []
for i in list(range(N)):
    params.append({"Time": 3e-3,
                   "focal_distance": dists[i],
                   "FOV": 90,
                   "filename": "zoom_%s.png"%str(i).zfill(4),
                   "center_on_star": 1,
                   "res": 512}
    )
for th in range(360):
    i += 1
    params.append({"Time": 3e-3,
                   "focal_distance": dists[-1],
                   "FOV": 90,
                   "pan": th,
                   "tilt": 20*np.sin(1.5*th/180*np.pi),
                   "filename": "zoom_%s.png"%str(i).zfill(4),
                   "center_on_star": 1,
                   "res": 512}
    )

params = [params] # format as a list of lists with one entry per task
Nchunks = 16 
DoTasksForSimulation(snaps=snaps, tasks=tasks, task_params=params,nproc=16)
