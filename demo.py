from glob import glob
from CrunchSnaps import *
from tasks import *
from multiprocessing import Pool, cpu_count
from tasks import *
from natsort import natsorted
import h5py
import itertools
import numpy as np
from sys import argv

tasks = [SinkVisCoolMap,]
snaps = glob(argv[1]+"/snapshot*.hdf5")
#Nparams = len(snaps)
Tmax = 7e-3
params = []
Nframes = 3000
times = np.linspace(0,Tmax,Nframes)

Nangle = 1080
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
    params.append({"Time": 5e-3+(7e-3-5e-3)*(i-i0)/360, "pan": float(i), "rmax": 10, "filename": "cool__%s.png"%str(i).zfill(4), "res": 1024})    

params = params

params = [params,]
Nchunks = 20 #cpu_count()
#params_chunks = [params[i*len(params)//Nchunks:(i+1)*len(params)//Nchunks] for i in range(Nchunks)]


#def DoStuff(params_chunk):
DoTasksForSimulation(snaps=snaps, tasks=tasks, task_params=params,nproc=Nchunks)

#[DoStuff(chunk) for chunk in params_chunks]
#Pool(Nchunks).map(DoStuff, params_chunks)
