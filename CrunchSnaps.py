from multiprocessing import Pool, cpu_count
from tasks import *
from natsort import natsorted
import h5py
import itertools
import numpy as np

def DoTasksForSimulation(snaps=[], tasks=[], task_params=[]):
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
        N_params = len(snaps)
    else:
        N_params = len(task_params[0])    
    
    # don't yet know what the snapshot times are - get the snapshot times in a prepass
    snaptimes = []
    for s in snaps:
        with h5py.File(s, 'r') as F:
            snaptimes.append(F["Header"].attrs["Time"])
    snaptimes = np.array(snaptimes)
#    print(snaptimes)
    snapdict = dict(zip(snaptimes, snaps))
    
    # note that params must be sorted by time!
    if not task_params:
        task_params = [[{"Time": s} for s in snaptimes] for t in tasks]

    snapdata_buffer = {} # should store data from at most 2 snapshots
    
    for i in range(N_params): # this loop cannot be parallelized without sacrificing IO performance because it has a persistent IO buffer modified on the fly
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
#        print(list(snapdata_buffer.keys()))
        for k in list(snapdata_buffer.keys()):
            if k < t1: del snapdata_buffer[k] # delete if not needed for interpolation (including old interpolants)
        for t in t1, t2:
            if not t in snapdata_buffer.keys(): snapdata_buffer[t] = GetSnapData(snapdict[t], required_snapdata)

        ##########################  interpolation #############################################
        # now get the particle data we need for this time, interpolating if needed

#        print(time, [k for k in snapdata_buffer.keys()])
        if time in snapdata_buffer.keys(): # if we have the snapshot for this exact time in the buffer, no interpolation needed
            snapdata_for_thistime = snapdata_buffer[time]
        else:
            print("interpolating to time ", time)
            snapdata_for_thistime = SnapInterpolate(time, t1, t2, snapdata_buffer)
            snapdata_buffer[time] = snapdata_for_thistime

        ################# task execution  ####################################################
        # actually do the task - each method can optionally return data to be compiled in the pass through the snapshots
        
        data = [t.DoTask(snapdata_for_thistime) for t in task_instances] 


def SnapInterpolate(t,t1,t2,snapdata_buffer):
    wind_ids = np.array([1913298393, 1913298394])
    
    stuff_to_interp_lin = "PartType0/Coordinates", "PartType0/Velocities", "PartType0/MagneticField", "PartType5/Coordinates", "PartType5/Velocities"
    stuff_to_interp_log = "PartType0/SmoothingLength", "PartType0/InternalEnergy", "PartType0/Pressure", "PartType0/SoundSpeed", "PartType0/Density", "PartType5/Masses"
    interpolated_data = snapdata_buffer[t1].copy()
    print(interpolated_data["PartType0/Coordinates"][0])
    # idx1, idx2, id1_order, id2_order = {}, {}, {}, {}
    # for ptype in "PartType0", "PartType5":
    #     id1 = snapdata_buffer[t1][ptype + "/ParticleIDs"]
    #     id2 = snapdata_buffer[t2][ptype + "/ParticleIDs"]

    #     if ptype == "PartType0": # if we have to worry about wind IDs
    #         wind_idx1 = np.in1d(id1, wind_ids)
    #         if np.any(wind_idx1):
    #             progenitor_ids = np.int_(snapdata_buffer[t1]["PartType0/ParticleIDGenerationNumber"])[wind_idx1]
    #             child_ids = np.int_(snapdata_buffer[t1]["PartType0/ParticleChildIDsNumber"])[wind_idx1]
    #             wind_particle_ids = -((progenitor_ids << 16) + child_ids) # bit-shift the progenitor ID outside the plausible range for particle count, then add child ids to get a unique new id                                         
    #             id1[wind_idx1] = wind_particle_ids

    #         wind_idx2 = np.in1d(id2, wind_ids)
    #         if np.any(wind_idx2):
    #             progenitor_ids = np.int_(snapdata_buffer[t2]["PartType0/ParticleIDGenerationNumber"])[wind_idx2]
    #             child_ids = np.int_(snapdata_buffer[t2]["PartType0/ParticleChildIDsNumber"])[wind_idx2]
    #             wind_particle_ids = -((progenitor_ids << 16) + child_ids) # bit-shift the progenitor ID outside the plausible range for particle count, then add child ids to get a unique new id                                         
    #             id2[wind_idx2] = wind_particle_ids        
        
    #     common_ids = np.intersect1d(id1,id2)
    #     id1_order[ptype] = id1.argsort()
    #     id2_order[ptype] = id2.argsort()
    #     idx1[ptype] = np.in1d(np.sort(id1),common_ids)
    #     idx2[ptype] = np.in1d(np.sort(id2),common_ids)
        
    # for field in stuff_to_interp_lin:
    #     if field in interpolated_data.keys():
    #         if field == "Header": continue
    #         ptype = field.split("/")[0]
    #         interpolated_data[field] = snapdata_buffer[t1][field][id1_order[ptype]][idx1[ptype]] * (t2 - t)/(t2 - t1) + snapdata_buffer[t2][field][id2_order[ptype]][idx2[ptype]] * (t - t1)/(t2 - t1)
    # for field in stuff_to_interp_log:        
    #     if field in interpolated_data.keys():
    #         if field == "Header": continue
    #         ptype = field.split("/")[0]
    #         interpolated_data[field] = np.exp(np.log(snapdata_buffer[t1][field][id1_order[ptype]][idx1[ptype]]) * (t2 - t)/(t2 - t1) + np.log(snapdata_buffer[t2][field][id2_order[ptype]][idx2[ptype]]) * (t - t1)/(t2 - t1))
    return interpolated_data
        

def GetSnapData(snappath, required_snapdata):
    print("loading " + snappath)
    snapdata = {}
    with h5py.File(snappath,'r') as F:
        snapdata["Header"] = dict(F["Header"].attrs)
        for field in required_snapdata:
            if field in F.keys(): snapdata[field] = np.array(F[field])
    return snapdata
        
from glob import glob

tasks = [SimVis,SimVis2]
snaps = glob("M2e3*/output/snap*.hdf5")

params = []
Nangle = 360
for i in range(Nangle):
    params.append({"Time": 0.001, "pan": float(i)/Nangle * 360})
#params = [params,]
Nchunks = cpu_count()
params_chunks = [params[i*len(params)//Nchunks:(i+1)*len(params)//Nchunks] for i in range(Nchunks)]


def DoStuff(params_chunk):
    DoTasksForSimulation(snaps, tasks, [params_chunk,params_chunk])


Pool(1).map(DoStuff, params_chunks)


