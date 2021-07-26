from multiprocessing import Pool
from .tasks import *
from natsort import natsorted

def DoPassOnSimulation(snaps=[], times=[], tasks=[], task_params=None, N_threads=1):
    """
    Main CrunchSnaps routine, performs a list of tasks using (possibly interpolated) time series simulation data
    snaps - list of paths of simulation snapshots
    times - list of times at which to perform tasks 
    tasks - list of tasks to perform at each simulation time
    task_params - shape [N_times, N_tasks] array of dictionaries containing the parameters for each task and each time - if only (N_tasks,) list is provided this will be broadcast over the time dimension
    """

    ####################################################################################################
    # do a bunch of default behaviors, broadcasting of input parameters, figure out the timeline, etc.
    ####################################################################################################
    if len(snaps) == 0 or len(tasks) == 0: return # no work to do so just quit

    snaps = natsorted(snaps)

    if len(times) == 0:
        N_times = len(snaps)
    else:
        N_times = len(times)
        
    N_tasks = len(tasks)
    
    if task_params == None: # defaults for everything, make empty dicts to pass to tasks
        task_params = np.array([[{} for i in range(len(tasks))] for j in range(N_times)])
    elif len(task_params) == N_tasks: # broadcast over the times
        task_params = np.repeat(task_params, N_times)
    # set up assumed broadcasting over tasks as well?        
    

    if len(times) > 0: # we'be been given a list of times and don't yet know what the snapshot times are - get the snapshot times in a prepass
        snaptimes = []
        for s in snaps:
            with h5py.File(f, 'r') as F:
                snaptimes.append(F["Header"].attrs["Time"])
    else:
        snaptimes = []

    alldata = []
    for i in range(N_times):
        task_instances = [t(task_params[j]) for t in tasks] # instantiate a task object for each task to be done
        required_snapdata = set(sum([t.GetRequiredSnapdata() for t in tasks])) # get the set of required fields, e.g. PartType0/Coordinates
        snapdata = GetSnapData(i, times, snaps, snaptimes, required_snapdata)
        data = [t.DoTask(fields) for t in tasks] # actually do the task - each method can optionally return data to be compiled in the pass through the snapshots
        alldata.append(data)

        

def GetSnapData(i, times, snaps, snaptimes, required_snapdata):
    snapdata = {}
    # case 1: the desired time is precisely equal to the desired time - no interpolation required
    if len(times) > 0:
        if times[i] in snaptimes:
            num_snaps_required = 1
    elif len(times) == 0:
        num_snaps_required = 1
    else:
        num_snaps_required = 2

    snap1 = h5py.File(snaps[i],'r')
    if num_snaps_required == 2:
        snap2 = h5py.File(snaps[i+1])

        
    for key in required_snapdata:
        
        

        
                
# MAIN LOOP:
# 1. initialize task objects with parameters and get required fields and constraints from them
# 2. do IO to get the required fields
# 3. do interpolation to current time if needed
# 4. do tasks, provided field data
# 5. retrieve data
