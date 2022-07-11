from multiprocessing import Pool, cpu_count
from .snapshot_tasks import *
from .misc_functions import *
from natsort import natsorted
import h5py
import itertools
import numpy as np
from functools import partial
from os.path import isfile
from numba import vectorize, njit
from math import copysign

def DoTasksForSimulation(snaps=[], tasks=[], task_params=[], interp_fac=1, nproc=1, nthreads=1, snaptime_dict=None):
    """
    Main CrunchSnaps routine, performs a list of tasks using (possibly interpolated) time series simulation data
    snaps - list of paths of simulation snapshots
    tasks - list of tasks to perform at each simulation time
    task_params - shape [N_tasks,N_params] list of dictionaries containing the parameters for each task - if only (N_tasks,) list is provided this will be broadcast
    """

    ####################################################################################################
    # do a bunch of default behaviors, broadcasting of input parameters, figure out the timeline, etc.
    ####################################################################################################
    if len(snaps) == 0 or len(tasks) == 0: return # no work to do so just quit

    snaps = natsorted(snaps)
        
    N_tasks = len(tasks)
        
    # broadcast task params across the different tasks if necessary
    if len(tasks) > 1 and task_params: # if we are doing more than one type of task and have given parameters
        for i in range(len(task_params)): # loop over parameters list
            if type(tasks[i]) == dict: # if we have a list of dicts
                tasks[i] = [tasks[i], tasks[i]]
            elif type(tasks[i]) == list: # if we have a list of lists
                if len(tasks[i]) < N_tasks: 
                    tasks[i] = [tasks[0], tasks[0]]            
        
    if snaptime_dict is None:
        snaptime_dict = get_snapshot_time_dict(snaps)

    snapnums = np.array([snapnum_from_path(s) for s in snaps])
    snaptimes = np.array([snaptime_dict[snapnum_from_path(s)] for s in snaps])
    snapdict = dict(zip(snaptimes, snaps))
#    print(snapnums, snaptimes, snapdict)

    if (not task_params) or (type(task_params) == dict): # if task_params is empty or a single dict that we must broadcast
        N_params = len(snaps)*interp_fac # default to just doing a pass on snapshots with optional interpolation and default parameters
        if interp_fac > 1: params_times = np.interp(np.arange(N_params)/interp_fac,np.arange(N_params//interp_fac),snaptimes)
        else: params_times = snaptimes
        if not task_params: task_params = [[{"Time": params_times[i], "index": i, "threads": nthreads} for i in range(N_params)] for t in tasks] # need to generate params from scratch
        else: # already have the dict of defaults; need to broadcast
            task_params_orig = task_params.copy()
            task_params = [[task_params_orig.copy() for i in range(N_params)] for t in tasks]
            [[task_params[j][i].update({"Time": params_times[i], "index": i, "threads": nthreads}) for i in range(N_params)] for j in range(N_tasks)] # add the other defaults
    else:
        N_params = len(task_params[0])    
    # note that params must be sorted by time!
    
    index_chunks = np.array_split(np.arange(N_params), nproc)
    chunks=[(index_chunks[i], tasks, snaps, task_params, snapdict, snaptimes, snapnums) for i in range(nproc)]
    if nproc > 1:
        Pool(nproc).map(DoParamsPass, chunks,chunksize=1) # this is where we fork into parallel tasks
    else:
        [DoParamsPass(c) for c in chunks]

def DoParamsPass(chunk):
    task_chunk_indices, tasks, snaps, task_params, snapdict, snaptimes, snapnums = chunk # unpack chunk data
    N_tasks = len(tasks)
    snapdata_buffer = {} # should store data from at most 2 snapshots
    for i in task_chunk_indices:        
        ######################### initialization  ############################################
        # initialize the task objects and figure out what data the tasks need
        task_instances = [tasks[n](task_params[n][i]) for n in range(N_tasks)] # instantiate a task object for each task to be done
        if np.all([t.TaskDone for t in task_instances]): continue
        time = task_params[0][i]["Time"]
        required_snapdata = set()
        for t in task_instances:
            fields = t.GetRequiredSnapdata() # get the set of required fields, e.g. PartType0/Coordinates
            for f in fields: required_snapdata.add(f)            
            for ptype in "PartType0","PartType5":
                for f in fields:
                    if ptype in f: required_snapdata.add(ptype + "/ParticleIDs")

        ############################ IO ######################################################
        # OK now we do file I/O if we don't find what we need in the buffer
        if len(snaptimes) >= 2:
            if time < snaptimes.max(): t1, t2 = snaptimes[snaptimes<=time][-1], snaptimes[snaptimes>=time][0] # if the time is within our timeline
            else: t1, t2 = snaptimes[-2:] # else if we have >1 snapshot, let those be the two times we interpolate/extrapolate from
        else: t1=t2=snaptimes[0] # otherwise we have exactly 1 snapshot, so set t1=t2 and ignore all time dependence
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
            if t1 != t2:
                snapdata_for_thistime = SnapInterpolate(time, t1, t2, snapdata_buffer)
            else:
                snapdata_for_thistime = snapdata_buffer[t1]
            snapdata_buffer[time] = snapdata_for_thistime

        ################# task execution  ####################################################
        # actually do the task - each method can optionally return data to be compiled in the pass through the snapshots

        data = [t.DoTask(snapdata_for_thistime) for t in task_instances]    


def SnapInterpolate(t,t1,t2,snapdata_buffer, min_val=1e-40):
    stuff_to_skip = 'Header', 'PartType0/ParticleIDs', 'PartType5/ParticleIDs', 'PartType0/ParticleIDGenerationNumber', 'PartType0/ParticleChildIDsNumber'
    #stuff_to_interp_lin = "PartType0/Coordinates", "PartType0/Velocities", "PartType0/MagneticField", "PartType5/Coordinates", "PartType5/Velocities",
    stuff_to_interp_log = "PartType0/SmoothingLength", "PartType0/InternalEnergy", "PartType0/Pressure", "PartType0/SoundSpeed", "PartType0/Density", "PartType5/Masses", "PartType5/BH_Mass", "PartType0/Masses", "PartType0/ElectronAbundance", "PartType0/HII", "PartType0/Temperature",
    interpolated_data = snapdata_buffer[t1].copy()

    idx1, idx2 = {}, {}
    for ptype in "PartType0","PartType5":
        if ptype+"/ParticleIDs" in snapdata_buffer[t1].keys(): id1 = snapdata_buffer[t1][ptype+"/ParticleIDs"]
        else: id1 = np.array([])
        if ptype+"/ParticleIDs" in snapdata_buffer[t2].keys(): id2 = snapdata_buffer[t2][ptype+"/ParticleIDs"]
        else: id2 = np.array([])
        common_ids = np.intersect1d(id1,id2)
        idx1[ptype] = np.in1d(np.sort(id1),common_ids)
        idx2[ptype] = np.in1d(np.sort(id2),common_ids)    
        
    
    wt1, wt2 = (t2 - t)/(t2 - t1), (t - t1)/(t2 - t1)
    for field in snapdata_buffer[t1].keys():
        if not (field in stuff_to_skip):
            ptype = field.split("/")[0]
            f1, f2 = snapdata_buffer[t1][field][idx1[ptype]], snapdata_buffer[t2][field][idx2[ptype]]
            if field in stuff_to_interp_log:        
                f1 += min_val; f2 += min_val; #nonzero min_val added to avoid issues
                interpolated_data[field] = np.exp(np.log(f1) * wt1 + np.log(f2) * wt2)
            else: #interpolate everything else linearily
                if "Coordinates" in field: # special behaviour to handle periodic BCs
                    dx = f2 - f1
                    dx = NearestImage(dx,snapdata_buffer[t1]["Header"]["BoxSize"])
                    interpolated_data[field] = f1 + wt2 * dx
                else:                
                    interpolated_data[field] = f1 * wt1 + f2 * wt2
                                       
    return interpolated_data
        

def GetSnapData(snappath, required_snapdata):
    ptypes_toread = set()
    for s in required_snapdata:
        for i in range(6):
            if "PartType%d"%i in s: ptypes_toread.add("PartType%d"%i)
        
    wind_ids = np.array([1913298393, 1913298394])
    snapdata = {}
    print("opening ", snappath)
    with h5py.File(snappath,'r') as F:
        snapdata["Header"] = dict(F["Header"].attrs)
        header = snapdata["Header"]
        time = header["Time"]
        z = header["Redshift"]
        hubble = header["HubbleParam"]
        # attempt to guess if this is a cosmological simulation from the agreement or lack thereof between time and redshift. note at t=1,z=0, even if non-cosmological, this won't do any harm
        if(np.abs(time*(1.+z)-1.) < 1.e-6): 
            cosmological=True; ascale=time;
        else:
            cosmological = False; ascale=1;

        for field in required_snapdata:
            if field in F.keys():
                snapdata[field] = F[field][:]
                if "ID" in field: snapdata[field] = np.int_(snapdata[field]) # cast to int for things that should be signed integers

                if cosmological:
                    if "Coordinates" in field or "SmoothingLength" in field: 
                        snapdata[field] *= time/hubble
                    if "Mass" in field: snapdata[field] /= hubble
                    if "Velocities" in field: snapdata[field] *= time**0.5
                    if "Density" in field or "Pressure" in field: snapdata[field] *= 1/hubble / (time/hubble)**3

    id_order = {}  # have to pre-sort everything by ID and fix the IDs of the wind particles
    for ptype in ptypes_toread:
        if not ptype+"/ParticleIDs" in snapdata.keys(): continue
        ids = snapdata[ptype + "/ParticleIDs"]

        if ptype == "PartType0": # if we have to worry about wind IDs and splitting
            child_ids = snapdata["PartType0/ParticleChildIDsNumber"]
            wind_idx1 = np.in1d(ids, wind_ids)
            ids[np.invert(wind_idx1)] = ((child_ids << 32) + ids)[np.invert(wind_idx1)]
            if np.any(wind_idx1):
                progenitor_ids = snapdata["PartType0/ParticleIDGenerationNumber"][wind_idx1]
                child_ids = child_ids[wind_idx1]
                wind_particle_ids = -((progenitor_ids << 16) + child_ids) # bit-shift the progenitor ID outside the plausible range for particle count, then add child ids to get a unique new id                                         
                ids[wind_idx1] = wind_particle_ids
            

        unique, counts = np.unique(ids, return_counts=True)
        doubles = unique[counts>1]
        ids[np.in1d(ids,doubles)]=-1

        id_order[ptype] = ids.argsort()

    for field in snapdata.keys():
        if field == "Header": continue
        ptype = field.split("/")[0]
        snapdata[field] = np.take(snapdata[field], id_order[ptype],axis=0)
    return snapdata
        
