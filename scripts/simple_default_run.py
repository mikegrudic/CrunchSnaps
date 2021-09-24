from glob import glob
from CrunchSnaps import *
from multiprocessing import Pool, cpu_count
from natsort import natsorted
import h5py
import itertools
import numpy as np
from sys import argv
from time import time

tasks = [SinkVisSigmaGas, SinkVisCoolMap]
snaps = natsorted(glob(argv[1]+"/snap*.hdf5"))
Nchunks = cpu_count()

params = {"res": 512}

t = time()
DoTasksForSimulation(snaps, task_params=params, interp_fac=1,tasks=tasks,nproc=4, nthreads=4)
print(time() - t)
