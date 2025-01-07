import h5py
from CrunchSnaps import *
from sys import argv

task = SinkVisCoolMap

for f in argv[1:]:
    snapnum = f.split("snapshot_")[1].split(".hdf5")[0]
    params = {
        "camera_distance": 0,
        "Time": 0,
        "fresco_stars": True,
        "index": int(snapnum),
        "pan": 0,
        "no_timestamp": True,
        "threads": -1,
        "no_stars": False,
        "res": 512,
        "cmap": "inferno",
        "limits": [10, 1e4],
        "tilt": -90,
        "threads": 1,
    }
    params = cubemapify(params)
    snapdata = None
    for p in params:
        print(p)
        t = task(p)
        if not snapdata:
            fields = t.GetRequiredSnapdata()
            snapdata = GetSnapData(f, fields)
        t.DoTask(snapdata)
