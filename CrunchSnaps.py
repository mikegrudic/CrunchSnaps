from multiprocessing import Pool, cpu_count
from tasks import *
from natsort import natsorted
import h5py
import itertools
import numpy as np
from functools import partial


def DoTasksForSimulation(snaps=[], tasks=[], task_params=[], interp_fac=1, nproc=1):
    """
    Main CrunchSnaps routine, performs a list of tasks using (possibly interpolated) time series simulation data
    snaps - list of paths of simulation snapshots
    tasks - list of tasks to perform at each simulation time
    task_params - shape [N_tasks,N_params] array of dictionaries containing the parameters for each task - if only (N_tasks,) list is provided this will be broadcast
    """

    ####################################################################################################
    # do a bunch of default behaviors, broadcasting of input parameters, figure out the timeline, etc.
    ####################################################################################################
    if len(snaps) == 0 or len(tasks) == 0: return # no work to do so just quit

    snaps = natsorted(snaps)
        
    N_tasks = len(tasks)

    if not task_params:
        N_params = len(snaps)*interp_fac
    else:
        N_params = len(task_params[0])    
    
    # don't yet know what the snapshot times are - get the snapshot times in a prepass
    snaptimes = []
    for s in snaps:
        with h5py.File(s, 'r') as F:
            snaptimes.append(F["Header"].attrs["Time"])
    snaptimes = np.array(snaptimes)    
    snapdict = dict(zip(snaptimes, snaps))
    
    if not task_params: # set up a default params with the desired interpolated times as parameters
        if interp_fac > 1: params_times = np.interp(np.arange(N_params)/interp_fac,np.arange(N_params//interp_fac),snaptimes)
        else: params_times = snaptimes
        task_params = [[{"Time": params_times[i]} for i in range(N_params)] for t in tasks]
    # note that params must be sorted by time!        

    index_chunks = np.array_split(np.arange(N_params), nproc)
    chunks=[(index_chunks[i], tasks, snaps, task_params, snapdict, snaptimes) for i in range(nproc)]            
    Pool(nproc).map(DoParamsPass, chunks) # this is where we fork into parallel tasks 

def DoParamsPass(chunk):
    task_chunk_indices, tasks, snaps, task_params, snapdict, snaptimes = chunk
    N_tasks = len(tasks)
    snapdata_buffer = {} # should store data from at most 2 snapshots        
    for i in task_chunk_indices:        
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
        for k in list(snapdata_buffer.keys()):
            if k < t1: del snapdata_buffer[k] # delete if not needed for interpolation (including old interpolants)
        for t in t1, t2:
            if not t in snapdata_buffer.keys(): snapdata_buffer[t] = GetSnapData(snapdict[t], required_snapdata)

        ##########################  interpolation #############################################
        # now get the particle data we need for this time, interpolating if needed

        if time in snapdata_buffer.keys(): # if we have the snapshot for this exact time in the buffer, no interpolation needed
            snapdata_for_thistime = snapdata_buffer[time]
        else:
            snapdata_for_thistime = SnapInterpolate(time, t1, t2, snapdata_buffer)
            snapdata_buffer[time] = snapdata_for_thistime

        ################# task execution  ####################################################
        # actually do the task - each method can optionally return data to be compiled in the pass through the snapshots

        data = [t.DoTask(snapdata_for_thistime) for t in task_instances]    


def SnapInterpolate(t,t1,t2,snapdata_buffer):
    stuff_to_interp_lin = "PartType0/Coordinates", "PartType0/Velocities", "PartType0/MagneticField", "PartType5/Coordinates", "PartType5/Velocities",
    stuff_to_interp_log = "PartType0/SmoothingLength", "PartType0/InternalEnergy", "PartType0/Pressure", "PartType0/SoundSpeed", "PartType0/Density", "PartType5/Masses", "PartType5/BH_Mass"
    interpolated_data = snapdata_buffer[t1].copy()

    idx1, idx2 = {}, {}
    for ptype in "PartType0", "PartType5":
        if ptype+"/ParticleIDs" in snapdata_buffer[t1].keys(): id1 = snapdata_buffer[t1][ptype+"/ParticleIDs"]
        else: id1 = np.array([])
        if ptype+"/ParticleIDs" in snapdata_buffer[t2].keys(): id2 = snapdata_buffer[t2][ptype+"/ParticleIDs"]
        else: id2 = np.array([])                    
        common_ids = np.intersect1d(id1,id2)
        idx1[ptype] = np.in1d(np.sort(id1),common_ids)
        idx2[ptype] = np.in1d(np.sort(id2),common_ids)    
    
    wt1, wt2 = (t2 - t)/(t2 - t1), (t - t1)/(t2 - t1)
    for field in stuff_to_interp_lin:
        if field in interpolated_data.keys():
            ptype = field.split("/")[0]
            interpolated_data[field] = snapdata_buffer[t1][field][idx1[ptype]] * wt1 + snapdata_buffer[t2][field][idx2[ptype]] * wt2
    for field in stuff_to_interp_log:        
        if field in interpolated_data.keys():
            ptype = field.split("/")[0]
            interpolated_data[field] = np.exp(np.log(snapdata_buffer[t1][field][idx1[ptype]]) * wt1 + np.log(snapdata_buffer[t2][field][idx2[ptype]]) * wt2)
    return interpolated_data
        

def GetSnapData(snappath, required_snapdata):
    wind_ids = np.array([1913298393, 1913298394])
    snapdata = {}
    with h5py.File(snappath,'r') as F:
        snapdata["Header"] = dict(F["Header"].attrs)
        for field in required_snapdata:
            if field in F.keys():
                snapdata[field] = np.array(F[field])
                if "ID" in field: snapdata[field] = np.int_(snapdata[field]) # cast to int for things that should be integers


    id_order = {}  # have to pre-sort everything by ID and fix the IDs of the wind particles
    for ptype in "PartType0", "PartType5": 
        if not ptype+"/ParticleIDs" in snapdata.keys(): continue
        ids = snapdata[ptype + "/ParticleIDs"]

        if ptype == "PartType0": # if we have to worry about wind IDs
            wind_idx1 = np.in1d(ids, wind_ids)
            if np.any(wind_idx1):
                progenitor_ids = snapdata["PartType0/ParticleIDGenerationNumber"][wind_idx1]
                child_ids = snapdata["PartType0/ParticleChildIDsNumber"][wind_idx1]
                wind_particle_ids = -((progenitor_ids << 16) + child_ids) # bit-shift the progenitor ID outside the plausible range for particle count, then add child ids to get a unique new id                                         
                ids[wind_idx1] = wind_particle_ids

        id_order[ptype] = ids.argsort()

    for field in snapdata.keys():
        if field == "Header": continue
        ptype = field.split("/")[0]
        snapdata[field] = np.take(snapdata[field], id_order[ptype],axis=0)
#        unique, counts = np.unique(id2, return_counts=True)
#        doubles = unique[counts>1]
#        id2[np.in1d(id2,doubles)]=-1                                                                                                                             
    return snapdata
        


