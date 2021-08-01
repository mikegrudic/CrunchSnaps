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
from time import time

tasks = [SinkVisSigmaGas,]
snaps = natsorted(glob(argv[1]+"/snap*.hdf5"))[:250]
Nchunks = cpu_count()

t = time()
DoTasksForSimulation(snaps,interp_fac=5,tasks=tasks,nproc=Nchunks)
print time() -t 
