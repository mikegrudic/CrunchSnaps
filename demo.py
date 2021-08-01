from glob import glob
from CrunchSnaps import *
from tasks import *
from multiprocessing import Pool, cpu_count
from tasks import *
from natsort import natsorted
import h5py
import itertools
import numpy as np

tasks = [SinkVisSigmaGas,]
snaps = glob("M2e3*/output/snap*.hdf5")

params = []
Nangle = 1080
for i in list(range(Nangle)):
    params.append({"Time": 5e-3 * i/Nangle,
                   "pan": float(i)/Nangle * 360,
                   "tilt": 30*np.sin(2*np.pi*i/Nangle),                   
                   "rmax": 3*np.exp(np.cos((5e-3*i/Nangle)/1e-3)**2),                   
                   "filename": "sigma_gas_%s.png"%str(i).zfill(4)}
    )
for i in range(Nangle-1,Nangle+360):
    params.append({"Time": 5e-3, "pan": float(i),  "rmax": 3*np.exp(np.cos((5e-3)/1e-3)**2), "filename": "sigma_gas_%s.png"%str(i).zfill(4)})
i0 = i
for i in range(i0+1,i0+360):
    params.append({"Time": 5e-3+(0.00983896001630692-5e-3)*(i-i0)/360, "pan": float(i0), "rmax": 3*np.exp(np.cos((5e-3)/1e-3)**2), "filename": "sigma_gas_%s.png"%str(i).zfill(4)})    

#params = [params,]
Nchunks = cpu_count()
params_chunks = [params[i*len(params)//Nchunks:(i+1)*len(params)//Nchunks] for i in range(Nchunks)]


def DoStuff(params_chunk):
    DoTasksForSimulation(snaps, tasks, [params_chunk,params_chunk])

#[DoStuff(chunk) for chunk in params_chunks]
Pool(Nchunks).map(DoStuff, params_chunks)
