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
snaps = natsorted(glob(argv[1] + "/snap*.hdf5"))
np = cpu_count()

params = {"res": 512, "fresco_stars": True}

t = time()
DoTasksForSimulation(snaps, task_params=params, interp_fac=4, tasks=tasks, nproc=np, nthreads=1)
print(time() - t)
