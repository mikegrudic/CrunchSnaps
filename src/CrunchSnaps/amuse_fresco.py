#!/usr/bin/env python
#Used by SinkVis to create observation-like images using amuse-fresco (https://pypi.org/project/amuse-fresco/). See install.txt for instructions!

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from amuse.datamodel import Particles
from amuse.units import units, nbody_system
from amuse.community.sse.interface import SSE
#from amuse.ext.masc import make_a_star_cluster
from amuse.ext import masc
from amuse.ext.fresco import make_fresco_image
import h5py



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

def make_amuse_fresco_stars_only(x,mstar,age_yr,L,res=512,p=5e-4,mass_limits=[0,0],mass_rescale=1.,filename=None,vmax=None):
#    print(x,mstar)
    number_of_stars = len(mstar)

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
    
    
    #inspectvar(units)
    stars.luminosity = lum_MS(mstar_new) | units.LSun
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
