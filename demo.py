from glob import glob
from CrunchSnaps import *
from tasks import *
from natsort import natsorted
import numpy as np
from sys import argv

tasks = [SinkVisCoolMap,] # select the tasks that you want to run
snaps = natsorted(glob(argv[1]+"/snapshot*.hdf5")) # get the list of snapshots to run on (should sort chronologically for best access pattern)

# Now here we generate a list of parameters for each frame we want to make, this will be done "by hand" and all unintialized parameters will adopt default values
Nangle = 1080
params = []
for i in list(range(Nangle)):
    params.append({"Time": 5e-3 * i/Nangle,
                   "pan": float(i)/Nangle * 360,
                   "tilt": 30*np.sin(2*np.pi*i/Nangle),                   
                   "rmax": 10*(1+np.exp(-5*i/Nangle)),
                   "filename": "cool__%s.png"%str(i).zfill(4),
                   "res": 1024}
    )
for i in range(Nangle-1,Nangle+360):
    params.append({"Time": 5e-3, "pan": float(i),  "rmax": 10, "filename": "cool__%s.png"%str(i).zfill(4),"res": 1024})
i0 = i
for i in range(i0+1,i0+360):
    params.append({"Time": 5e-3+(7e-3-5e-3)*(i-i0)/360, "pan": float(i0) + (i-i0)*0.3, "rmax": 10 + (i-i0)*10/360, "filename": "cool__%s.png"%str(i).zfill(4), "res": 1024})    
i0 = i
for i in range(i0+1,i0+1080):
    params.append({"Time": 7e-3, "pan": float(i0) + (i-i0)*0.3, "rmax": 10 + (i-i0)*10/1080, "filename": "cool__%s.png"%str(i).zfill(4), "res": 1024})    

params = [params] # format as a list of lists with one entry per task
Nchunks = 16 
DoTasksForSimulation(snaps=snaps, tasks=tasks, task_params=params,nproc=16)
