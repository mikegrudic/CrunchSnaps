#!/usr/bin/env python
#Used by SinkVis to create observation-like images using amuse-fresco (https://pypi.org/project/amuse-fresco/). See install.txt for instructions!

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from amuse.datamodel import Particles
from amuse.units import units, nbody_system
#from amuse.community.sse.interface import SSE
#from amuse.ext.masc import make_a_star_cluster
#from amuse.ext import masc
from amuse.plot.fresco import make_fresco_image
import h5py
from scipy.spatial.distance import pdist,squareform
from scipy.spatial import cKDTree
#from numba import njit




#Calculate the luminosity of a main sequence star in code units: ORION version: fitting formulas of Tout et al (1996) 
def lum_MS(m):
    m = np.array(m,dtype=np.float64)
    lum = ( 0.39704170*(m**5.5) + 8.52762600*(m**11) ) / (0.00025546+(m**3)+5.43288900*(m**5) + 5.56357900*(m**7) + 0.78866060*(m**8)+0.00586685*(m**9.5))
    lum[m<0.1] = 0 #BD are not on MS
    return lum 

#Calculate the radius of a main sequence star in solar units: ORION version: fitting formulas of Tout et al (1996) 
def rad_MS(m):
    m = np.array(m,dtype=np.float64)
    rad =  (1.71535900*(m**2.5)+6.59778800*(m**6.5)+10.08855000*(m**11)+1.01249500*(m**19)+0.07490166*(m**19.5)) / (0.01077422+3.08223400*(m**2)+17.84778000*(m**8.5)+(m**18.5)+0.00022582*(m**19.5));
    return rad 

def inspectvar(var):
    import inspect
    attributes = inspect.getmembers(var, lambda a:not(inspect.isroutine(a)))
    att_list = [a for a in attributes if not(a[0].startswith('__') and a[0].endswith('__'))]
    print('Printing attributes')
    for a in att_list:
        print(a)
    return

def Cartesian_to_Spherical(x_cart, map_to_positive=False):
    if not (type(x_cart)==np.ndarray): x_cart=np.array(x_cart)
    if len(x_cart.shape)==1: x_cart.shape=(1,3)
    x_sph = np.zeros(x_cart.shape)
    r = x_cart[:,0]**2 + x_cart[:,1]**2
    x_sph[:,0] = np.sqrt(r + x_cart[:,2]**2)
    x_sph[:,1] = np.arctan2(np.sqrt(r), x_cart[:,2]) # for elevation angle defined from Z-axis down
    x_sph[:,2] = np.arctan2(x_cart[:,1], x_cart[:,0])
    if map_to_positive:
        neg_phi = x_sph[:,2]<0
        x_sph[neg_phi,2] += 2*np.pi #instead of -Pi,Pi we map to 0,2*Pi
    return x_sph


def optical_depth_for_targets(target_coords,x_obs,x,m,h,kappa,nthreads=1):  
    #Calculate the optical depth between a set of targets and an observer
    optical_depths = np.zeros(len(target_coords),dtype=np.float32) #init
    #Shorthand, roughly the optical length through a cell
    d_tau = kappa*m/(h*h)
    #Center everything on observer
    dx = x - x_obs; target_coords = target_coords - x_obs
    ##################
    #Convert to Spherical, centered on observer and roughly looking towards where most stars are
    x_target_mean = np.mean(target_coords,axis=0)
    #Get normalized vectors
    e1= x_target_mean/np.linalg.norm(x_target_mean)
    if np.abs(e1[2])<1:
        e2 = np.cross(e1, [0,0,1])
    else:
        e2 = np.cross(e1, [0,1,0]) #exception if the observer is looking directly at the target
    e2 /= np.linalg.norm(e2)
    e3 = np.cross(e1,e2); e3 /= np.linalg.norm(e3)
    #Transform into new Cartesian coordinates
    T = np.vstack( (e1,e2,e3) ).T
    dx = dx@T
    target_coords = target_coords@T
    #Transform to spherical
    dx_sph = Cartesian_to_Spherical(dx,map_to_positive=True)
    target_coords_sph = Cartesian_to_Spherical(target_coords,map_to_positive=True)
    ##############
    #Build a tree in a fake 3D space that we will query for gas cells along lines of sights, based on https://github.com/mikegrudic/ComputeColumns/blob/master/ComputeColumns.py
    R = dx_sph[:,0]
    r_eff = h/R #The angular size of the softening
    K = r_eff.max() * (1+1e-6) # largest possible gas cell radius
    w = np.sqrt(K**2 - r_eff**2) # this is the fake z coordinate in the 3D space we're embedding the tree in
    x_fake = np.c_[dx_sph[:,1],dx_sph[:,2],w] # construct the fake 3D coordinates
    #print("\t Coordinate transform for %d sources and %d gas done, building KDTree..."%(target_coords_sph.shape[0],len(R)))
    tree3d = cKDTree(x_fake, leafsize=64) # construct the fake 3D tree
    #print("\t KDTree built, start neighbor queries...")
    #Query tree for each target       
    target_coords_fake = np.vstack((target_coords_sph[:,1],target_coords_sph[:,2],np.zeros(target_coords_sph.shape[0]))).T
    ind_all = tree3d.query_ball_point(target_coords_fake,K,workers=nthreads) # particles within distance K in the fake 3D space will be within their actual radius in the 2D space
    #print("\t Neighbor queries done, integrating along field lines...")
    for i,ind in enumerate(ind_all):
        #Choose only those between observer and target
        between_ind = R[ind]<=target_coords_sph[i,0]
        #Find distance to sightline normalized by cell size (<1 by definition)
        delta = dx_sph[ind][between_ind] - target_coords_sph[i]
        q = np.sqrt(delta[:,1]**2+delta[:,2]**2)/r_eff[ind][between_ind]
        #Calculate contributions for each cell
        kernel = np.zeros_like(q) #numerical factor for how much the particle overlap with our path
        ind1 = (q<=0.5); ind2 = (q>0.5) & (q<1.0)
        kernel[ind1] = 1 - 6*q[ind1]*q[ind1] * (1-q[ind1])
        kernel[ind2] = 2.0 * (1-q[ind2])**3
        optical_depths[i] = np.sum(1.8189136353359467 * kernel * d_tau[ind][between_ind]) #This means we take either 100% or zero of a gas cell's d_tau, with a hard transition at the softening length 
    return  optical_depths

def make_amuse_fresco_stars_only(x,mstar,age_yr,L,lum=None,stage=None,res=512,p=5e-4,mass_limits=[0,0],mass_rescale=1.,filename=None,vmax=None,optical_depth=None, MS_only=False):
    number_of_stars = len(mstar)

    if optical_depth is None:
        optical_depth = np.zeros_like(mstar)
    if (mass_limits[0]!=0.0) or (mass_limits[1]!=0.0):
        if (mass_limits[0]==0.0): mass_limits[0]=np.min(mstar)
        if (mass_limits[1]==0.0): mass_limits[1]=np.max(mstar)
        mstar_new = np.clip(mstar,mass_limits[0],mass_limits[1])
        #small correction to make them stand apart (e.g. instea of clipping both 100 and 50 msun to 50, we get 50 and 60  so they don't look identical)
        mstar_new = mstar_new + (mstar-mstar_new)*0.1
        ##limits masses of star
        #logm = np.log10(mstar); logm0 = np.max(logm) + np.min(logm);
        #mstar_new = 10**( (logm - logm0)/mass_limits + logm0 ) 
    mstar_new = mstar**mass_rescale
    new_stars = Particles(number_of_stars)
    new_stars.age = age_yr | units.yr
    new_stars.mass = mstar_new  | units.MSun
    new_stars.position = x | units.pc
    stars = new_stars
    gas = Particles()
    
    #use luminsoity if provided
    lum_final = lum_MS(mstar_new)
    if (lum is not None): #take the lower of the rescaled MS and the current luminosity
        ind = (lum_final > lum)
        lum_final[ind] = lum[ind]
    if MS_only and (stage is not None):
        lum_final[stage<5] = 0
    
    #inspectvar(units)
    stars.luminosity = ( np.exp(-optical_depth)*lum_final ) | units.LSun
    #stars.main_sequence_lifetime = [450.0, 420.0] | units.Myr
    stars.radius = rad_MS(mstar_new) | units.RSun
    #stars.spin = [4700, 4700] | units.yr**-1
    # stars.stellar_type = [1, 1] | units.stellar_type
    
    # se = SSE()
    # se.particles.add_particles(stars)
    # from_se = se.particles.new_channel_to(stars)
    # from_se.copy()
    # inspectvar(stars)
    image, _ = make_fresco_image( stars, gas, return_vmax=True,\
        image_width=[L | units.pc,L | units.pc], image_size=[res,res],percentile=1-p,vmax=vmax)
    #image: (2048,2048,3) RGB
    if not(filename is None):
        #Save image to file
        plt.imshow(image[::-1],extent=(-L/2.0,L/2.0,-L/2.0,L/2.0))
        plt.xlim(-L/2.0,L/2.0)
        plt.ylim(-L/2.0,L/2.0)
        plt.imsave(filename,image[::-1])
    
    return image[::-1]
