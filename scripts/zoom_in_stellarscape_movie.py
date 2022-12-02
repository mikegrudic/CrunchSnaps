import numpy as np
from numpy.polynomial import Polynomial
import os
from scipy import interpolate
from scipy import optimize


#Parameters
Tstart=0.0030 #in code units, units are roughly 1 Gyr
Tend = 0.0070186 #ends at a snapshot, to make the zoom in easier (278 in this case)
fps = 10 #60
total_length_seconds = 1.5*60 #10/1.5*60
Nchunks = 1 #16 #number of folders to distribute the files, i.e. number of separate slurm jobs to render
pan_seconds_per_rotation = 80 #80
tilt_period_seconds = 100 #100

#Zooming (extra segment added at the end)
zoom_in_time_seconds = 10 # set 0 to turn it off
zoom_dt = 0.0002 #in code units, how much time elapses in the simulation during the zoom
zoom_in_coord = np.array([47.34047, 48.086853, 48.908130]) #coordinates to zoom in on at the end
boxsize =100

center_coord = np.ones(3)/2*boxsize #center of the box
total_length_seconds_nozoom = total_length_seconds - zoom_in_time_seconds #non zoom part
zoom_in_coord = np.array(zoom_in_coord)

#Segments before zoom, this is et up for a 2 segment video + zoom in at the end
T1 = 0.0045 #first segment end time
segment1_realtime_end = 10 #length of initial video segment in seconds

segment_sim_times = np.array([Tstart, T1, Tend]) #Bookmarks for video segments (e.g., frozen or slowed down parts of the video). The times in code units, use [Tstart, Tend] if there are no segments.
segment_real_times = np.array([0, segment1_realtime_end, total_length_seconds_nozoom]) #real time of video segment bookmarks
segment_frames_per_sim_time = fps * np.diff(segment_real_times)/np.diff(segment_sim_times)
Nsegments = len(segment_sim_times) - 1

def angle_positive_map(x):
    return (x+2*np.pi)%(2*np.pi)

def Spherical_to_Euclidean(r,theta,phi):
    z = r*np.cos(theta)
    y = r*np.sin(theta)*np.sin(phi)
    x = r*np.sin(theta)*np.cos(phi)
    return x,y,z
    
def Euclidean_to_Spherical(x,y,z):
    r = np.sqrt( x**2 + y**2 + z**2 )
    theta = angle_positive_map(np.arccos(z/r))
    phi = angle_positive_map(np.arctan2(y,x))
    return r,theta,phi

def poly_exp_decay_func(x,params):
    #f(x) = polynom(x)*exp(x/scale) + conv_val
    poly_params = params[:-2]; scale = params[-2]; conv_val = params[-1] #unpack
    return ( Polynomial(poly_params)(x) * np.exp(-x/scale) + conv_val )
    
def poly_exp_decay_func_deriv(x,params):
    #Derivative of f(x) = polynom(x)*exp(x/scale) + conv_val
    poly_params = params[:-2]; scale = params[-2]; conv_val = params[-1] #unpack
    return ( np.exp(-x/scale) * ( Polynomial(poly_params).deriv(1)(x) - Polynomial(poly_params)(x)/scale ) )
    
def poly_exp_decay_func_deriv2(x,params):
    #Second derivative of f(x) = polynom(x)*exp(x/scale) + conv_val
    poly_params = params[:-2]; scale = params[-2]; conv_val = params[-1] #unpack
    return ( np.exp(-x/scale) * ( Polynomial(poly_params).deriv(2)(x) - Polynomial(poly_params).deriv(1)(x)/scale ) -  poly_exp_decay_func_deriv(x,params)/scale)


def fitfunc(params,args,regularization_param=1e-4):
    y0, v0, a0, y1, v1, dx = args #unpack target values
    tot_err = 0 #total error
    weights = [100,10,1,100,1] #relative weigts of constraints
    for target, func, dist, wt in zip( [y0, v0, a0, y1, v1, dx], [poly_exp_decay_func,poly_exp_decay_func_deriv,poly_exp_decay_func_deriv2, poly_exp_decay_func, poly_exp_decay_func_deriv], [0,0,0,dx,dx], weights ):
        # if target:
            # tot_err += wt*( 1.0 - func(dist,params)/target )**2 #match BC
        # else: #need to handle 0
            # tot_err += wt*func(dist,params)**2 
        tot_err += wt*(func(dist,params)-target)**2 
    #regularization error
    tot_err += regularization_param * np.sum(params**2)
    return tot_err

def zoom_between(y,y1,n, func='exp_poly2', verbose=False, errtol=1e-2): #movement of camera during zoom, default is linear
    if hasattr(y, "__iter__"): #spline interpolation
        x = np.append(np.arange(len(y))-len(y)+1,[n])
        x_interp = np.linspace(x[-2],x[-1],num=n+1)
        if ( (len(y)<3) or func=='spline' ):
            if len(y)==1: 
                kind = 'slinear'; 
            if len(y)<4: 
                kind = 'quadratic'; 
            else:
                kind = 'cubic'
            y_full = np.append(y,[y1]);
            return interpolate.interp1d(x,y_full,kind=kind)(x_interp)
        elif func=='exp_poly2':
            y0 = y[-1]; v0 = np.diff(y)[-1]; a0 = np.diff(np.diff(y))[-1];
            v1=0 #let's have zero final velocity
            init_guess = [a0/2 + 4*(y0-y1)/(n**2) + 4*(v0+2*(y0-y1)/n)/n,v0+(y0-y1)/(n/2),y0-y1,n/2,y1] #initial guess for second order case, using partial solutions for first 3 constraints
            BC_params = [y0, v0, a0, y1, v1, n]
            if verbose:
                print("Boundary conditions parameters params ", BC_params)
                print("Initial guess for parameters ", init_guess)
            opt = optimize.minimize(fitfunc, init_guess, args=BC_params, method = 'BFGS', options={"disp":verbose, 'eps': 1e-8, 'maxiter': 1000}) #need tgo reduce eps if it gives inaccurate fits
            if verbose: 
                print("Fitted params ", opt.x)
                print("Fit error %g"%(opt.fun))
                print("At position %g:"%(x[-2]))
                print("\t y0: %g fit %g"%(y0,poly_exp_decay_func(x[-2],opt.x)))
                print("\t v0: %g fit %g"%(v0,poly_exp_decay_func_deriv(x[-2],opt.x)))
                print("\t a0: %g fit %g"%(a0,poly_exp_decay_func_deriv2(x[-2],opt.x)))
                print("At position %g:"%(x[-1]))
                print("\t y1: %g fit %g"%(y1,poly_exp_decay_func(x[-1],opt.x)))
                print("\t v1: %g fit %g"%(v1,poly_exp_decay_func_deriv(x[-1],opt.x)))
            if np.sqrt( (y0-poly_exp_decay_func(x[-2],opt.x))**2 + (y1-poly_exp_decay_func(x[-1],opt.x))**2 ) > max(errtol,np.abs(y0)*errtol,np.abs(y1)*errtol) :
                raise Exception("Fit does not match endpoints within error tolerance!")
            x_interp = np.linspace(x[-2],x[-1],num=n+1)[1:-1]#skip last frame which would be on top of target
            return poly_exp_decay_func(x_interp,opt.x)
        else:
            print("No known function called "+func)
    else: #simple linear
        return np.linspace(y,y1,num=n+1)[1:-1] #skip last frame which would be on top of target
        


if np.any(np.diff(segment_sim_times)<=0): 
    print("Error! segment_sim_times not monotonically increasing"); exit()
if np.any(np.diff(segment_real_times)<=0): 
    print("Error! segment_real_times not monotonically increasing"); exit()

sim_times = np.array([])
for i in range(Nsegments):
    sim_times = np.append(sim_times, np.arange(segment_sim_times[i], segment_sim_times[i+1], 1/segment_frames_per_sim_time[i]) )
sim_times = np.array(sim_times)
real_times = np.linspace(0, total_length_seconds_nozoom, len(sim_times))

# tot_timelapse_time = total_length_seconds_nozoom - sum(freeze_intervals)
# tot_timelapse_frames = fps * tot_timelapse_time
# frames_per_sim_time = tot_timelapse_time / (Tmax-Tstart) * fps
# frame_dt = 1/frames_per_sim_time

# sim_times = []
# t0 = Tstart
# for i in range(len(freeze_times)):
    # dt = freeze_times[i] - t0
    # sim_times.append(np.arange(t0, freeze_times[i], frame_dt))
    # if len(freeze_intervals):
        # sim_times.append(np.repeat(freeze_times[i], round(freeze_intervals[i] * fps)))
    # t0 = freeze_times[i]+frame_dt
# sim_times = np.concatenate(sim_times)
# real_times = np.linspace(0,total_length_seconds_nozoom,int(total_frames))
# if len(sim_times) < len(real_times):
    # sim_times = np.append(sim_times,sim_times[-1])

#Camera track before zooming in
pan = (360 / pan_seconds_per_rotation * real_times) % 360
tilt = 26 + 25 * np.sin(2*np.pi*real_times / tilt_period_seconds)
#r = 7*(1+4*np.exp(-20*real_times / total_length_seconds_nozoom) + 2*np.exp((sim_times*979  - 7.5) / 1)) * np.exp((np.cos(4*np.pi * real_times / total_length_seconds_nozoom) - 1)/2) #more complicated path, better for longer movies
r = 7*(1+4*np.exp(-20*real_times / total_length_seconds_nozoom) ) * np.exp((np.cos(2*np.pi * real_times / total_length_seconds_nozoom) - 1)/2)

#Track coordinates
theta = angle_positive_map(tilt/180*np.pi)
phi = angle_positive_map(pan/180*np.pi)
dx,dy,dz = Spherical_to_Euclidean(r,theta,phi)
directions = -np.c_[dx/r,dy/r,dz/r] #camera directions os far

#Zoom in
if zoom_in_time_seconds:
    #dir_0 = directions[-1,:]; dir_0 /= np.linalg.norm(dir_0) #last viewing direction at the end
    dir_final = zoom_in_coord-center_coord-np.array([dx[-1],dy[-1],dz[-1]]); dir_final /= np.linalg.norm(dir_final)
    zoom_in_frame_number = int(zoom_in_time_seconds*fps)
    #Append zoom in to data so far
    dx=np.append(dx,zoom_between(dx[:-1],zoom_in_coord[0]-center_coord[0],zoom_in_frame_number))
    dy=np.append(dy,zoom_between(dy[:-1],zoom_in_coord[1]-center_coord[1],zoom_in_frame_number))
    dz=np.append(dz,zoom_between(dz[:-1],zoom_in_coord[2]-center_coord[2],zoom_in_frame_number))
    r,theta,phi = Euclidean_to_Spherical(dx,dy,dz)
    new_dir = np.c_[zoom_between(directions[:-1,0],dir_final[0],zoom_in_frame_number),zoom_between(directions[:-1,1],dir_final[1],zoom_in_frame_number),zoom_between(directions[:-1,2],dir_final[2],zoom_in_frame_number)]
    directions = np.vstack( (directions,new_dir) )
    #normalize just to be sure
    directions = directions/np.linalg.norm(directions,axis=1)[:,None]
    real_times = np.append(real_times, real_times[-1]+(np.arange(zoom_in_frame_number-1)+1)/fps)
    sim_times = np.append(sim_times, sim_times[-1]+(np.arange(zoom_in_frame_number-1)+1)/zoom_in_frame_number*zoom_dt)
camera_pos = np.c_[dx,dy,dz] + center_coord
camera_dist = np.zeros_like(sim_times) #since we are using explicit camera coordinates instead of pan/tilt 

#Save camerafile
data = np.c_[sim_times, camera_pos, directions, camera_dist]
np.savetxt("camerafile.txt",data)
chunks = np.array_split(data,Nchunks)
for i in range(Nchunks):
    if not os.path.exists("%d"%(i)):
        os.mkdir("%d"%(i))
    np.savetxt("%d/camerafile_%d.txt"%(i,i), chunks[i])
                        


################
#Plotting
norm_time = real_times/np.max(real_times)
xlim=ylim=[-5,5]
zlim=[-1,15]
N_frames_tot = len(real_times)
v = np.c_[np.diff(dx),np.diff(dy),np.diff(dz)]/np.diff(real_times)[:,None]
vel = np.linalg.norm(v,axis=1)
a = np.c_[np.diff(dx,2),np.diff(dy,2),np.diff(dz,2)]/(np.diff(real_times)[:-1,None]**2)
acc = np.linalg.norm(a,axis=1)
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

fig = plt.figure()
for i in range(N_frames_tot-1):
    plt.plot(real_times[i:i+2], r[i:i+2], c = plt.get_cmap('viridis')(norm_time[i]))
plt.xlabel('T [s]'); plt.ylabel('R [pc]')
fig.savefig('Path_Dist_t.png', dpi=200)

fig = plt.figure()
for i in range(len(vel)-1):
    plt.plot(real_times[i:i+2], vel[i:i+2], c = plt.get_cmap('viridis')(norm_time[i]))
plt.xlabel('T [s]'); plt.ylabel('v [pc/s]')
fig.savefig('Path_vel_t.png', dpi=200)

fig = plt.figure()
for i in range(len(acc)-1):
    plt.plot(real_times[i:i+2], acc[i:i+2], c = plt.get_cmap('viridis')(norm_time[i]))
plt.xlabel('T [s]'); plt.ylabel('Acceleration [$\mathrm{pc/s^2}$]')
plt.ylim([0,np.percentile(acc,99)])
fig.savefig('Path_acc_t.png', dpi=200)

fig = plt.figure()
for i in range(N_frames_tot-1):
    plt.plot(real_times[i:i+2], dz[i:i+2], c = plt.get_cmap('viridis')(norm_time[i]))
plt.xlabel('T [s]'); plt.ylabel('Z [pc]')
fig.savefig('Path_z_t.png', dpi=200)

fig = plt.figure()
for i in range(N_frames_tot-1):
    plt.plot(real_times[i:i+2], dx[i:i+2], c = plt.get_cmap('viridis')(norm_time[i]))
plt.xlabel('T [s]'); plt.ylabel('X [pc]')
fig.savefig('Path_x_t.png', dpi=200)

# fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# for i in range(N_frames_tot-1):
    # ax.plot((phi[i:i+2]+2*np.pi)%(2*np.pi), r[i:i+2], c = plt.get_cmap('viridis')(norm_time[i]))
# ax.set_rmax(15)
# fig.savefig('Path_polar_2D.png', dpi=200)

fig = plt.figure()
for i in range(N_frames_tot-1):
    plt.plot(dx[i:i+2], dy[i:i+2], c = plt.get_cmap('viridis')(norm_time[i]))
plt.scatter(zoom_in_coord[0]-center_coord[0],zoom_in_coord[1]-center_coord[1],marker='x',c='r')
plt.scatter(0,0,marker='o',c='b')
plt.xlim(xlim); plt.ylim(xlim); plt.xlabel('X [pc]'); plt.ylabel('Y [pc]')
fig.savefig('Path_xy.png', dpi=200)

fig = plt.figure()
for i in range(N_frames_tot-1):
    plt.plot(dx[i:i+2], dz[i:i+2], c = plt.get_cmap('viridis')(norm_time[i]))
plt.scatter(zoom_in_coord[0]-center_coord[0],zoom_in_coord[2]-center_coord[2],marker='x',c='r')
plt.scatter(0,0,marker='o',c='b')
plt.xlim(xlim); plt.ylim(zlim); plt.xlabel('X [pc]'); plt.ylabel('Z [pc]')
fig.savefig('Path_xz.png', dpi=200)

fig = plt.figure()
for i in range(N_frames_tot-1):
    plt.plot(dy[i:i+2], dz[i:i+2], c = plt.get_cmap('viridis')(norm_time[i]))
plt.scatter(zoom_in_coord[0]-center_coord[0],zoom_in_coord[2]-center_coord[2],marker='x',c='r')
plt.scatter(0,0,marker='o',c='b')
plt.xlim(xlim); plt.ylim(zlim); plt.xlabel('Y [pc]'); plt.ylabel('Z [pc]')
fig.savefig('Path_yz.png', dpi=200)

fig = plt.figure()
for i in range(N_frames_tot-1):
    plt.plot(dy[i:i+2], dz[i:i+2], c = plt.get_cmap('viridis')(norm_time[i]))
plt.scatter(zoom_in_coord[1]-center_coord[1],zoom_in_coord[2]-center_coord[2],marker='x',c='r')
plt.scatter(0,0,marker='o',c='b')
plt.xlim(ylim); plt.ylim(zlim); plt.xlabel('Y [pc]'); plt.ylabel('Z [pc]')
fig.savefig('Path_yz.png', dpi=200)

fig = plt.figure()
ax = plt.axes(projection='3d')
for i in range(N_frames_tot-1):
    ax.plot3D(dx[i:i+2], dy[i:i+2],dz[i:i+2], c = plt.get_cmap('viridis')(norm_time[i]))
ax.scatter(zoom_in_coord[0]-center_coord[0], zoom_in_coord[1]-center_coord[1], zoom_in_coord[2]-center_coord[2], marker='x',c='r')
ax.scatter(0, 0, 0, marker='o',c='b')
ax.set_xlim(xlim[0], xlim[1])
ax.set_ylim(ylim[0], ylim[1])
ax.set_zlim(zlim[0], zlim[1])
fig.savefig('Path_3D.png', dpi=200)