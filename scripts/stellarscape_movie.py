import numpy as np

Tmax = 0.00917
fps = 60
total_length_seconds = 10 * 60

total_frames = fps * total_length_seconds

pan_seconds_per_rotation = 120
tilt_period_seconds = 150
# res = 1024

freeze_times = [6.0 / 979, Tmax]
freeze_intervals = [30, 30]  # [30, 30, 30]
tot_timelapse_time = total_length_seconds - sum(freeze_intervals)
tot_timelapse_frames = fps * tot_timelapse_time
frames_per_sim_time = tot_timelapse_time / Tmax * fps
frame_dt = 1 / frames_per_sim_time

sim_times = []
t0 = 0
for i in range(len(freeze_times)):
    dt = freeze_times[i] - t0
    Nframes_timelapse = round(dt * frames_per_sim_time)

    sim_times.append(np.arange(t0, freeze_times[i], frame_dt))
    sim_times.append(np.repeat(freeze_times[i], round(freeze_intervals[i] * fps)))
    t0 = freeze_times[i] + frame_dt

sim_times = np.concatenate(sim_times)
real_times = np.linspace(0, total_length_seconds, total_frames)
if len(sim_times) < len(real_times):
    sim_times = np.append(sim_times, sim_times[-1])


pan = (360 / pan_seconds_per_rotation * real_times) % 360
tilt = 25 * np.sin(2 * np.pi * real_times / tilt_period_seconds)

camera_dist = (
    7
    * (1 + 4 * np.exp(-real_times / 10) + 2 * np.exp((sim_times * 979 - 7.5) / 1))
    * np.exp((np.cos(4 * np.pi * real_times / total_length_seconds) - 1) / 2)
)


data = np.c_[sim_times, camera_dist, pan, tilt]
np.savetxt("camerafile.txt", data)
Nchunks = 8
chunks = np.array_split(data, Nchunks)

for i in range(Nchunks):
    np.savetxt("%d/camerafile_%d.txt" % (i, i), chunks[i])
