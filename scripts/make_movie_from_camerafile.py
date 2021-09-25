import numpy as np
from sys import argv
from natsort import natsorted
from CrunchSnaps import *
from glob import glob

#N = 100
#time = np.linspace(0,5e-3,N)
#camera_pos = np.c_[np.zeros(N), np.zeros(N), np.logspace(2,0,N)] + np.array([15,15,15])
#camera_pos = np.concatenate([camera_pos, 
#camera_dir = np.gradient(camera_pos,axis=0) 
#camera_dist = np.zeros(N)
#np.savetxt("camerafile.txt",np.c_[time, camera_pos, camera_dir, camera_dist])

camera_data = np.loadtxt(argv[1])
sim_dir = argv[2]

time = camera_data[:,0]
camera_pos = camera_data[:,1:4]
camera_dir = camera_data[:,4:7]
camera_dist = camera_data[:,7]

params = []
tasks = [SinkVisCoolMap,]
for i in range(len(time)):
    params.append({"Time": time[i], "center": camera_pos[i], "camera_dir": camera_dir[i], "camera_distance": camera_dist[i], "index": i})

j = i
for i in range(100):
    j += 1
    params.append({"Time": time[-1], "center": camera_pos[-1], "camera_dir": camera_dir[-1], "camera_distance": float(i), "index": j})  

snaps = natsorted(glob(sim_dir + "/snapshot*.hdf5"))

DoTasksForSimulation(snaps, tasks=tasks, task_params=[params],nproc=16)
