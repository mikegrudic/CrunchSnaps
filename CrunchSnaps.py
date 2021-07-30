from multiprocessing import Pool
from tasks import *
from natsort import natsorted
import h5py
import itertools
import numpy as np

def DoPassOnSimulation(snaps=[], tasks=[], task_params=[], N_threads=1):
    """
    Main CrunchSnaps routine, performs a list of tasks using (possibly interpolated) time series simulation data
    snaps - list of paths of simulation snapshots
    times - list of times at which to perform tasks 
    tasks - list of tasks to perform at each simulation time
    task_params - shape [N_tasks,N_params] array of dictionaries containing the parameters for each task and each time - if only (N_tasks,) list is provided this will be broadcast over the time dimension
    """

    ####################################################################################################
    # do a bunch of default behaviors, broadcasting of input parameters, figure out the timeline, etc.
    ####################################################################################################
    if len(snaps) == 0 or len(tasks) == 0: return # no work to do so just quit

    snaps = natsorted(snaps)
        
    N_tasks = len(tasks)

    if not task_params:
        N_params = len(snaps)
    else:
        N_params = len(task_params)    
    
    # don't yet know what the snapshot times are - get the snapshot times in a prepass
    snaptimes = []
    for s in snaps:
        with h5py.File(s, 'r') as F:
            snaptimes.append(F["Header"].attrs["Time"])
    snaptimes = np.array(snaptimes)
    print(snaptimes)
    snapdict = dict(zip(snaptimes, snaps))
    
    # note that params must be sorted by time!
    if not task_params:
        task_params = [[{"Time": s} for s in snaptimes] for t in tasks]

    snapdata_buffer = {} # should store data from at most 2 snapshots
    
    for i in range(N_params):
        ######################### initialization  ############################################
        # initialize the task objects at figure out what data the tasks need
        
        task_instances = [tasks[n](task_params[n][i]) for n in range(N_tasks)] # instantiate a task object for each task to be done
        time = task_params[0][i]["Time"]
        required_snapdata = set()
        for t in task_instances:
            fields = t.GetRequiredSnapdata() # get the set of required fields, e.g. PartType0/Coordinates
            for f in fields: required_snapdata.add(f)
            for ptype in "PartType0","PartType5": required_snapdata.add(ptype + "/ParticleIDs")

        ############################ IO ######################################################
        # OK now we do file I/O if we don't find what we need in the buffer
        t1, t2 = snaptimes[snaptimes<=time][-1], snaptimes[snaptimes>=time][0]
        # do a pass to delete anything that will no longer be needed
        print(list(snapdata_buffer.keys()))
        for k in list(snapdata_buffer.keys()):
            if k < t1: del snapdata_buffer[k]
        for t in t1, t2:
            if not t in snapdata_buffer.keys(): snapdata_buffer[t] = GetSnapData(snapdict[t], required_snapdata)

        ##########################  interpolation #############################################
        # now get the particle data we need for this time, interpolating if needed
        
        if time in snapdata_buffer.keys(): # if we have the snapshot for this exact time in the buffer, no 
            snapdata_for_thistime = snapdata_buffer[t]
        else:
            snapdata_for_thistime = SnapInterpolate(time, t1, t2, snapdata_buffer)

        
        ################# task execution  ####################################################
        data = [t.DoTask(snapdata_for_thistime) for t in task_instances] # actually do the task - each method can optionally return data to be compiled in the pass through the snapshots


def SnapInterpolate(t,t1,t2,snapdata_buffer):    
    stuff_to_interpolate_linearly = "PartType0/Coordinates", "PartType0/Velocities", "PartType0/MagneticField", "PartType5/Coordinates", "PartType5/Velocities"
    stuff_to_interp_logarithmically = "PartType0/SmoothingLength", "PartType0/InternalEnergy", "PartType0/Pressure", "PartType0/SoundSpeed", "PartType0/Density"

def GetSnapData(snappath, required_snapdata):
    print("loading " + snappath)
    snapdata = {}
    with h5py.File(snappath,'r') as F:
        snapdata["Header"] = dict(F["Header"].attrs)
        for field in required_snapdata:
            if field in F.keys(): snapdata[field] = np.array(F[field])
    return snapdata
        
from glob import glob

tasks = [SurfaceDensityMap,]
snaps = glob("M2e3*/output/snap*.hdf5")


DoPassOnSimulation(snaps, tasks)
