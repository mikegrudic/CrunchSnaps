import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LightSource
from meshoid import GridSurfaceDensity, GridRadTransfer
import aggdraw
from skimage.color import rgb2hsv, hsv2rgb
from PIL import Image, ImageDraw, ImageFont, ImageChops
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from .amuse_fresco import *
from numba import get_num_threads, set_num_threads
from .misc_functions import *
from os.path import isfile
import json
import os
import sys
hashseed = os.getenv('PYTHONHASHSEED')
if not hashseed:
    os.environ['PYTHONHASHSEED'] = '0'
    os.execv(sys.executable, [sys.executable] + sys.argv)

class Task:
    """Class containing generic routines common to all tasks, and assigns default (null/empty) attributes that any task should have"""
    def __init__(self,params):
        self.RequiredSnapdata = []
        self.params = params

    def GetRequiredSnapdata(self):
        return self.RequiredSnapdata

    def AssignDefaultParams(self):
        for k in self.default_params.keys():
            if not k in self.params.keys(): self.params[k] = self.default_params[k]

    def AssignDefaultParamsFromSnapdata(self,snapdata):
        return

    

class SinkVis(Task):
    def __init__(self, params):
        """Class containing methods for coordinate transformations, rendering, etc. for a generic SinkVis-type map plot"""
        super().__init__(params)        

        self.default_params = {"Time": 0,
                               "res": 512,
                               "rmax": None,
                               "limits": [10,3e3],
                               "center": None,
                               "pan": 0,
                               "tilt": 0,
                               "no_timestamp": False,
                               "no_size_scale": False,
                               "filename": None,
                               "sink_scale": 1,
                               "cmap": "viridis",
                               "backend": "PIL",
                               "rescale_hsml": False,
                               "FOV": 90,
                               "camera_distance": np.inf,
                               "center_on_star": False,
                               "fresco_stars": False,
                               "fresco_param": 0.001,
                               "fresco_mass_limits": [0,0],
                               "fresco_mass_rescale": 0.3,
                               "threads": 1,
                               "cubemap_dir": "forward",
                               "camera_dir": None,
                               "camera_right": None,
                               "camera_up": None,
                               "index": None,
                               "no_stars": False,
                               "overwrite": False,
                               "unit_scalefac": 1,
                               "outputfolder": ".",
                               "SHO_RGB_norm": 0
        }


        self.AssignDefaultParams()
        basename, ext = os.path.splitext(self.params["filename"])
        self.params["filename_incomplete"] = basename+"_incomplete"+ext
        self.params_that_affect_maps = ["Time", "res", "rmax", "center", "pan", "tilt", "FOV", "camera_distance", "center_on_star", "cubemap_dir", "camera_dir", "camera_right", "camera_up", "rescale_hsml", "res"]
        dump = {}
        for k in self.params_that_affect_maps:
            if type(self.params[k]) == np.ndarray:
                dump[k] = self.params[k].tolist()
            else:
                dump[k] = self.params[k]
        self.params_hash = str(hash(json.dumps(dump,sort_keys=True)))

        mapdir = self.params["outputfolder"]+"/.maps"
        if not os.path.isdir(mapdir): os.mkdir(mapdir)
        self.map_files = dict([(m, mapdir+"/" + m + "_" + self.params_hash) for m in self.required_maps]) # filename for the saved maps will by MAPNAME_(hash # of input params)
        self.maps = {}

        self.DetermineRequiredSnapdata()

        if self.params["threads"] != 1:
            self.parallel = True
            if self.params["threads"] > 0: # if negative, just use all available threads, otherwise set to desired value
                set_num_threads(self.params["threads"])
        else: self.parallel = False

        if isfile(self.params["filename"]) and not self.params["overwrite"]:
            self.RequiredSnapdata = []
            self.TaskDone = True
        else:
            self.TaskDone = False

    def DetermineRequiredSnapdata(self):
        self.RequiredSnapdata = ["PartType5/Coordinates","PartType5/Masses","PartType5/ParticleIDs", "PartType5/BH_Mass"]

    def AssignDefaultParams(self):
        super().AssignDefaultParams()
        if self.params["index"] is None:
            self.params["index"] = round(self.params["Time"]/1e-6)            
        self.params["filename_suffix"] = "%s_%s_%s.png"%(str(self.params["index"]).zfill(5), str(round(self.params["pan"]*10)).zfill(4), self.params["cubemap_dir"]) 


    def CoordinateTransform(self,x,m=None,h=None, contravariant=False):
        # center on the designated center coordinate
        if not contravariant: x[:] -= self.params["center"]

        if self.params["camera_dir"] is None: # without a specified camera direction, we just use a simple tilt/pan scheme
            tilt, pan = self.params["tilt"], self.params["pan"]
            if contravariant: tilt, pan = -tilt, -pan
            # first pan
            cosphi, sinphi = np.cos(np.pi*pan/180), np.sin(np.pi*pan/180)
            x[:] = np.c_[cosphi*x[:,0] + sinphi*x[:,2],x[:,1], -sinphi*x[:,0] + cosphi*x[:,2]]
            # then tilt
            costheta, sintheta = np.cos(np.pi*tilt/180), np.sin(np.pi*tilt/180)
            x[:] = np.c_[x[:,0], costheta*x[:,1] + sintheta*x[:,2], -sintheta*x[:,1] + costheta*x[:,2]]
        else: # we have a camera position and coordinate basis
            if contravariant: x[:] = (self.camera_matrix_vectors @ x.T).T # note that @ performs matrix multiplication
            else: x[:] = (self.camera_matrix @ x.T).T

        if self.params["camera_distance"] != np.inf and not contravariant:
            # transform so camera is at z=0:
            x[:,2] += self.params["camera_distance"]
        
        # shuffle the axes to get the desired cubemap direction
        cubedir = self.params["cubemap_dir"]
        if cubedir != "forward":
            if cubedir == "right": x[:] = np.c_[-x[:,2],x[:,1],x[:,0]]
            elif cubedir == "left": x[:] = np.c_[x[:,2],x[:,1],-x[:,0]]
            elif cubedir == "up": x[:] = np.c_[x[:,0],-x[:,2],x[:,1]]
            elif cubedir == "down": x[:] = np.c_[x[:,0],x[:,2],-x[:,1]]
            elif cubedir == "backward": x[:] = np.c_[-x[:,0],x[:,1],-x[:,2]]                
        
        # then do projection if desired
        if self.params["camera_distance"] != np.inf:
            if not contravariant:
                # now transform from 3D to angular system
                r = np.sum(x*x,axis=1)**0.5 # distance from camera
                x[:,:2] = x[:,:2] / x[:,2][:,None] # homogeneous coordinates
                r = np.abs(x[:,2])
                self.r = r
                if h is not None:
                    h[:] = h / r # kernel lengths are now angular (divide by distance)
                    h[x[:,2]<0] = 0             # assign 0 weight/size to anything behind the camera
                if m is not None:
                    m[:] /= r**2 # rescale mass weights so that integrated surface density remains the same
                    m[x[:,2]<0] = 0

            else: # dealing with a contravariant vector such as velocity - want the [:,2] component to correspond to line-of-sight value
                global_coords = np.copy(self.pos) # this would have been converted to angular by now - let's convery back to real space
                global_coords[:,:2] *= global_coords[:,2][:,None] # multiply by z, now we're in the rotated real space frame
                x[:,2] = np.sum(x*global_coords,axis=1)/np.sum(global_coords**2,axis=1)**0.5 # get the radial component
                
            
            

    def SetupCoordsAndWeights(self, snapdata):
        res = self.params["res"]
        if "PartType0/Coordinates" in snapdata.keys():
            self.pos, self.mass, self.hsml = np.copy(snapdata["PartType0/Coordinates"]), np.copy(snapdata["PartType0/Masses"]), np.copy(snapdata["PartType0/SmoothingLength"]) # copy these because we don't want to modify them
            if self.params["rescale_hsml"]: self.hsml *= self.params["rescale_hsml"]

        # Setting up coordinate basis
        if self.params["camera_dir"] is not None:
            self.camera_dir = self.params["camera_dir"]
            NormalizeVector(self.params["camera_dir"])
            if not self.params["camera_up"]: self.camera_up = np.array([0,1.,0]) # default "up" direction is +y, we will project it out if the camera is tilted
            else: self.camera_up = self.params["camera_up"]
            
            # if we've specified an up direction, project out the component parallel to the forward direction and normalize
            self.camera_up -= sum(self.camera_dir * self.camera_up).sum() * self.camera_dir
            NormalizeVector(self.camera_up)
            # now get the "right" vector as the cross product of forward x up. this will be normalized to machine precision
            self.camera_right = np.cross(self.camera_up, self.camera_dir)

            self.camera_matrix = np.c_[self.camera_right, self.camera_up, self.camera_dir].T # matrix of coordinate vectors - operate this on coordinates to apply transformation - operates on COORDINATES not vectors
            self.camera_matrix_vectors = self.camera_matrix.T # since vector fields are contravariant, this is the operator for transforming v and B (note that this is an orthogonal matrix so the transpose is the inverse)
            
        if "PartType0/Coordinates" in snapdata.keys():             
            self.CoordinateTransform(self.pos,self.mass,self.hsml)
            self.hsml = np.clip(self.hsml,2*self.params["rmax"]/res, 1e100)


    def GenerateMaps(self,snapdata):
        return

    def SaveImage(self):
        print("saving ", self.params["filename"])        
        if self.params["backend"] == "matplotlib":
            rmax = self.params["rmax"]
            self.ax.set(xlim=[-rmax,rmax],ylim=[-rmax,rmax])
            plt.savefig(self.params["filename_incomplete"],bbox_inches='tight',dpi=200)            
            plt.close()
        os.rename(self.params["filename_incomplete"], self.params["filename"])            

    def MakeImages(self,snapdata):
        if not self.params["no_stars"]: self.AddStarsToImage(snapdata)
        self.AddSizeScaleToImage()
        self.AddTimestampToImage()
        self.SaveImage()
        

    def AddTimestampToImage(self):
        if self.params["no_timestamp"]: return
        fname = self.params["filename_incomplete"]
        time = self.params["Time"]
        if (time*979>=1e-2):
            time_text="%3.2gMyr"%(time*979)
        elif(time*979>=1e-4):
            time_text="%3.2gkyr"%(time*979*1e3)
        else:
            time_text="%3.2gyr"%(time*979*1e6)
            
        if self.params["backend"]=="PIL":
            F = Image.open(fname)
            gridres = F.size[0]
            draw = ImageDraw.Draw(F)
            font = ImageFont.truetype("LiberationSans-Regular.ttf", gridres//12)
            draw.text((gridres/16, gridres/24), time_text, font=font)
            F.save(fname)
            F.close()
        elif self.params["backend"]=="matplotlib":
            self.ax.text(-self.params["rmax"]*0.85, self.params["rmax"]*0.85,time_text,color="#FFFFFF")


    def AddSizeScaleToImage(self):
#        if self.params["camera_distance"] < np.inf: return
        if self.params["backend"]=="matplotlib": return # matplotlib will have axis ticks for scale
        pc_to_AU = 206265.0
        if self.params["no_size_scale"]: return
        fname = self.params["filename_incomplete"]
        F = Image.open(fname)
        draw = ImageDraw.Draw(F)
        gridres = self.params["res"]
        font = ImageFont.truetype("LiberationSans-Regular.ttf", gridres//12)
        r = self.params["rmax"] * self.params["unit_scalefac"]
        if self.params["camera_distance"] < np.inf:
            r = self.params["rmax"] * self.params["camera_distance"] * self.params["unit_scalefac"]
#        print("r=%g\n"%r)
        if (r*2>1000):
            scale_kpc=10**np.round(np.log10(r*0.5/1000))
            size_scale_text="%3.3gkpc"%(scale_kpc)
            size_scale_ending=gridres/16+gridres*(scale_kpc*1000)/(2*r)
        elif (r>1e-2):
            scale_pc=10**np.round(np.log10(r*0.5))
            size_scale_text="%3.3gpc"%(scale_pc)
            size_scale_ending=gridres/16+gridres*(scale_pc)/(2*r)
        else:
            scale_AU=10**np.round(np.log10(r*0.5*pc_to_AU))
            size_scale_text="%3.4gAU"%(scale_AU)
            size_scale_ending=gridres/16+gridres*(scale_AU)/(2*r*pc_to_AU)
        draw.line(((gridres/16, 7*gridres/8), (size_scale_ending, 7*gridres/8)), fill="#FFFFFF", width=6)
        draw.text((gridres/16, 7*gridres/8 + 5), size_scale_text, font=font)
        F.save(fname)
        F.close()                       
            
    def AddStarsToImage(self,snapdata):
#        print([k for k in snapdata.keys()])
        if not "PartType5/Coordinates" in snapdata.keys(): return
        X_star = np.copy(snapdata["PartType5/Coordinates"])
        m_star = snapdata["PartType5/BH_Mass"]

        self.CoordinateTransform(X_star, np.ones(len(X_star)), np.ones(len(X_star)))
        
        if self.params["backend"]=="PIL":
            fname = self.params["filename_incomplete"]
            if self.params["fresco_stars"]: # use fresco for stellar images
                if self.params["camera_distance"] < np.inf:
                    X_star, m_star = X_star[X_star[:,2]>0], m_star[X_star[:,2]>0]
#                    m_star /= X_star[:,2]**2
                if len(X_star) == 0: return
                data_stars_fresco = make_amuse_fresco_stars_only(X_star,m_star, np.zeros_like(m_star),2*self.params["rmax"],res=self.params["res"],vmax=self.params["fresco_param"],mass_rescale=self.params["fresco_mass_rescale"],mass_limits=self.params["fresco_mass_limits"])
                img = plt.imread(fname)
                plt.imsave(fname,np.clip(img[:,:,:3]+data_stars_fresco,0,1))
            else: # use derpy PIL circles
                F = Image.open(fname)
                gridres = F.size[0]
                draw = ImageDraw.Draw(F)
                d = aggdraw.Draw(F)
                pen = aggdraw.Pen(self.Star_Edge_Color(),1) #gridres/800
                sink_relscale = 0.0025
                X_star ,m_star = X_star[m_star.argsort()[::-1]], np.sort(m_star)[::-1]

                for j in np.arange(len(X_star))[m_star>1e-2]:                    
                    X = X_star[j]
                    ms = m_star[j]
                    star_size = max(1,gridres * sink_relscale * (np.log10(ms/self.params["sink_scale"]) + 1))
#                    if self.params["camera_distance"] < np.inf:
                        # make 100msun ~ 0.03pc, scale down from there
#                        if X[2] < 0: continue
#                        star_size = gridres * 0.03 / dist_to_camera / self.params["rmax"] * (ms/100)**(1./3)                        
                    star_size = max(1,star_size)
                    p = aggdraw.Brush(self.GetStarColor(ms))
                    norm_coords = (X[:2]+self.params["rmax"])/(2*self.params["rmax"])*gridres
                    #Pillow puts the origin in th top left corner, so we need to flip the y axis
                    norm_coords[1] = gridres - norm_coords[1]
                    coords = np.concatenate([norm_coords-star_size, norm_coords+star_size])
                    d.ellipse(coords, pen, p)#, fill=(155, 176, 255))
                    d.flush()
                F.save(fname)
                F.close()
        elif self.params["backend"]=="matplotlib":
            star_size = np.log10(m_star/self.params["sink_scale"])+2
            colors = np.array([self.GetStarColor(m) for m in m_star])/255
 
            self.ax.scatter(X_star[:,0], X_star[:,1],s=star_size*5,edgecolor=self.Star_Edge_Color(),lw=0.1,facecolor=colors,marker='*')
            
    def Star_Edge_Color(self):
        if self.params["cmap"] in ('afmhot', 'inferno', "Blues"):
            return 'black'
        else:
            return 'white'
        
    def GetStarColor(self, mass_in_msun):
        if self.params["cmap"] in ('afmhot', 'inferno', "Blues"):
            star_colors = np.array([[255, 100, 60],[120, 200, 150],[75, 80, 255]]) #alternate colors, red-green-blue, easier to see on a bright color map
        else:
            star_colors = np.array([[255, 203, 132],[255, 243, 233],[155, 176, 255]]) #default colors, reddish for small ones, yellow-white for mid sized and blue for large
        colors = np.int_([np.interp(np.log10(mass_in_msun),[-1,0,1],star_colors[:,i]) for i in range(3)])
        return (colors[0],colors[1],colors[2])# if len(colors)==1 else colors)

    def AssignDefaultParamsFromSnapdata(self,snapdata):
        if self.params["center"] is None:
            if self.params["center_on_star"]:
                if "PartType5/Coordinates" in snapdata.keys():
                    self.params["center"] = snapdata["PartType5/Coordinates"][snapdata["PartType5/BH_Mass"].argsort()[::-1]][self.params["center_on_star"]-1] # center on the n'th most massive star
                else: # otherwise center on the densest gas cell
                    self.params["center"] = snapdata["PartType0/Coordinates"][snapdata["PartType0/Density"].argmax()]
            else:
                self.params["center"] = np.repeat(snapdata["Header"]["BoxSize"]*0.5,3)
            center = self.params["center"]
        else: center = self.params["center"]
        if self.params["rmax"] is None:
            if self.params["camera_distance"] < np.inf:
                self.params["rmax"] = self.params["FOV"]/90 # angular width
            else:
                self.params["rmax"] = snapdata["Header"]["BoxSize"]/10
#            if self.params["camera_distance"] < np.inf and self.params["FOV"] is None:
#                self.params["rmax"] /= self.params["camera_distance"] # convert to angular assuming rmax is real-space half-width at the focal distance

    
    def DoTask(self, snapdata):
        if self.TaskDone: return 
        self.AssignDefaultParamsFromSnapdata(snapdata)
        if set(self.required_maps) != set(self.maps.keys()): # if we don't already have the maps we need
            self.SetupCoordsAndWeights(snapdata)        
        self.GenerateMaps(snapdata)
        self.MakeImages(snapdata)
        return self.maps


class SinkVisSigmaGas(SinkVis):
    def __init__(self,params):
        self.required_maps = ["sigma_gas"]                                 
        super().__init__(params)
        if self.TaskDone: return
        self.AssignDefaultParams()

    def AssignDefaultParams(self):
        super().AssignDefaultParams()
#         self.params["filename"] = "SurfaceDensity_%s_%s.png"%(str(self.params["index"]).zfill(4), str(round(self.params["pan"])).zfill(4))
#        else:
        if self.params["filename"] is None:
            self.params["filename"] = self.params["outputfolder"]+"/"+"SurfaceDensity_" + self.params["filename_suffix"]

    def DetermineRequiredSnapdata(self):
        super().DetermineRequiredSnapdata()
        # check if we have sigma_gas map already saved
        if isfile(self.map_files["sigma_gas"] + ".npz"):
            print("Loading sigma_gas map from %s"%(self.map_files["sigma_gas"]))
            self.maps["sigma_gas"] = np.load(self.map_files["sigma_gas"] + ".npz")['sigma_gas']            
        else:
            self.RequiredSnapdata += ["PartType0/Coordinates","PartType0/Masses","PartType0/ParticleIDs", "PartType0/BH_Mass", "PartType0/SmoothingLength","PartType0/ParticleChildIDsNumber","PartType0/ParticleIDGenerationNumber"]
        
    def GenerateMaps(self,snapdata):
        if not "sigma_gas" in self.maps.keys():
            self.maps["sigma_gas"] = GridSurfaceDensity(self.mass, self.pos, self.hsml, np.zeros(3), 2*self.params["rmax"], res=self.params["res"],parallel=self.parallel).T        
            np.savez_compressed(self.map_files["sigma_gas"], sigma_gas=self.maps["sigma_gas"])        

    def MakeImages(self,snapdata):
        vmin, vmax = self.params["limits"]
        f = (np.log10(self.maps["sigma_gas"])-np.log10(vmin))/(np.log10(vmax)-np.log10(vmin))

        if self.params["backend"]=="PIL":
            plt.imsave(self.params["filename_incomplete"], plt.get_cmap(self.params["cmap"])(np.flipud(f))) # NOTE - we invert this to get the coordinate system right
        elif self.params["backend"]=="matplotlib":
            self.fig, self.ax = plt.subplots(figsize=(4,4))
            X = Y = np.linspace(-self.params["rmax"], self.params["rmax"], self.params["res"])
            X, Y = np.meshgrid(X, Y)
            p = self.ax.pcolormesh(X, Y, self.maps["sigma_gas"], norm=matplotlib.colors.LogNorm(vmin=self.params["limits"][0],vmax=self.params["limits"][1]),cmap=self.params["cmap"])
            self.ax.set_aspect('equal')
 
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.0)
            self.fig.colorbar(p,label=r"$\Sigma_{\rm gas}$ $(\rm M_\odot\,pc^{-2})$",cax=cax)
            if self.params["camera_distance"] == np.inf:
                self.ax.set_xlabel("X (pc)")
                self.ax.set_ylabel("Y (pc)")
            else:
                self.ax.set_xlabel("X (rad)")
                self.ax.set_ylabel("Y (rad)")

        super().MakeImages(snapdata)




class SinkVisCoolMap(SinkVis):
    def __init__(self,params):
        self.required_maps = ["sigma_gas", "sigma_1D"] #physical rendered quantities that can get saved and reused
        super().__init__(params)
        if self.TaskDone: return        
        self.default_params["cool_cmap"] = 'magma'
        self.AssignDefaultParams()

    def DetermineRequiredSnapdata(self):
        super().DetermineRequiredSnapdata()
        # check if we have sigma_gas map already saved
        if isfile(self.map_files["sigma_gas"] + ".npz"):
            print("Loading sigma_gas map from %s"%(self.map_files["sigma_gas"]))
            self.maps["sigma_gas"] = np.load(self.map_files["sigma_gas"]+".npz")["sigma_gas"]            
        else:
            self.RequiredSnapdata += ["PartType0/Coordinates","PartType0/Masses","PartType0/ParticleIDs", "PartType0/BH_Mass", "PartType0/SmoothingLength","PartType0/ParticleChildIDsNumber","PartType0/ParticleIDGenerationNumber"]

        if isfile(self.map_files["sigma_1D"] + ".npz"):
            print("Loading sigma_1D map from %s"%(self.map_files["sigma_1D"]))
            self.maps["sigma_1D"] = np.load(self.map_files["sigma_1D"] + ".npz")["sigma_1D"]
        else:
            self.RequiredSnapdata += ["PartType0/Velocities"]
            self.RequiredSnapdata += ["PartType0/Coordinates","PartType0/Masses","PartType0/ParticleIDs", "PartType0/BH_Mass", "PartType0/SmoothingLength","PartType0/ParticleChildIDsNumber","PartType0/ParticleIDGenerationNumber"]
#        print(self.RequiredSnapdata)
        

    def AssignDefaultParams(self):
        super().AssignDefaultParams()
#        if self.params["filename"] is None: self.params["filename"] = "CoolMap_%s_%s.png"%(str(self.params["Index"]).zfill(4), str(round(self.params["pan"])).zfill(4))
#        else:
        if self.params["filename"] is None:
            self.params["filename"] = self.params["outputfolder"]+"/"+"CoolMap_" + self.params["filename_suffix"]
#            self.params["filename_incomplete"] = self.params["filename"].replace(".png",".incomplete.png")
        
        
    def GenerateMaps(self,snapdata):
        super().GenerateMaps(snapdata)
                
        if not "sigma_gas" in self.maps.keys():
            self.maps["sigma_gas"] = GridSurfaceDensity(self.mass, self.pos, self.hsml, np.zeros(3), 2*self.params["rmax"], res=self.params["res"],parallel=self.parallel).T        
            np.savez_compressed(self.map_files["sigma_gas"], sigma_gas=self.maps["sigma_gas"])
        if not "sigma_1D" in self.maps.keys():
            # need to apply coordinate transforms to z-velocity
            v = np.copy(snapdata["PartType0/Velocities"])
            self.CoordinateTransform(v,contravariant=True)
            sigma_1D = GridSurfaceDensity(self.mass * v[:,2]**2, self.pos, self.hsml, np.zeros(3), 2*self.params["rmax"], res=self.params["res"],parallel=self.parallel).T/self.maps["sigma_gas"]
            v_avg = GridSurfaceDensity(self.mass * v[:,2], self.pos, self.hsml, np.zeros(3), 2*self.params["rmax"], res=self.params["res"],parallel=self.parallel).T/self.maps["sigma_gas"]
            self.maps["sigma_1D"] = np.sqrt(sigma_1D - v_avg**2)/1e3
            np.savez_compressed(self.map_files["sigma_1D"], sigma_1D=self.maps["sigma_1D"])

    def MakeImages(self,snapdata):
        fgas = (np.log10(self.maps["sigma_gas"])-np.log10(self.params["limits"][0]))/np.log10(self.params["limits"][1]/self.params["limits"][0])
        fgas = np.clip(fgas,0,1)
        ls = LightSource(azdeg=315, altdeg=45)
        #lightness = ls.hillshade(z, vert_exag=4)
        mapcolor = plt.get_cmap(self.params["cool_cmap"])(np.log10(self.maps["sigma_1D"]/0.1)/2)
        cool_data = ls.blend_hsv(mapcolor[:,:,:3], fgas[:,:,None])
        self.maps["coolmap"] = cool_data
    
        plt.imsave(self.params["filename_incomplete"], np.flipud(self.maps["coolmap"])) # NOTE - we invert this to get the coordinate system right
        super().MakeImages(snapdata)


class SinkVisNarrowbandComposite(SinkVis):
    def __init__(self,params):
        self.required_maps = ["SHO_RGB"] #RGB map of SII, Halpha, and OIII
        super().__init__(params)
        if self.TaskDone: return
        self.AssignDefaultParams()

    def DetermineRequiredSnapdata(self):
        super().DetermineRequiredSnapdata()
        # check if we have maps already saved
        render_maps = False
        for mapname in self.required_maps:
            if not isfile(self.map_files[mapname]+".npz"): render_maps = True
        if render_maps:
            self.RequiredSnapdata += ["PartType0/Coordinates","PartType0/Masses","PartType0/ParticleIDs", "PartType0/Temperature", "PartType0/ElectronAbundance", "PartType0/SmoothingLength","PartType0/ParticleChildIDsNumber","PartType0/ParticleIDGenerationNumber", "PartType0/Density", "PartType0/HII"]
        else: # load pre-existing
            for mapname in self.required_maps:
                print("Loading %s map from %s"%(mapname, self.map_files[mapname]))
                self.maps[mapname] = np.load(self.map_files[mapname]+".npz")[mapname]            
                

    def AssignDefaultParams(self):
        super().AssignDefaultParams()
        if self.params["filename"] is None:
            self.params["filename"] = self.params["outputfolder"]+"/"+"NarrowbandComposite_" + self.params["filename_suffix"]
#            self.params["filename_incomplete"] = self.params["filename"].replace(".png",".incomplete.png")
        
        
    def GenerateMaps(self,snapdata):
        super().GenerateMaps(snapdata)
             
        if not "SHO_RGB" in self.maps.keys():
            rho = snapdata["PartType0/Density"]
            T = snapdata["PartType0/Temperature"]
            fe = snapdata["PartType0/ElectronAbundance"]
            hii = snapdata["PartType0/HII"]
            nH = rho * 30
            ne = nH * fe
            
            ne = np.clip(ne,None,np.percentile(ne,100*(1.0-100/len(ne))) ) #clip by 100th largest value in case we have few rogue cells with extremely large values

            T4 = T/1e4
            j_B_Ha = 1.24e-25 * (T4)**(-0.942-0.031 * np.log(T4)) * 2.86 * nH*hii * ne

            wavelength = 5007
            ncrit = 1e3
            j_OIII = nH*hii*ne * np.exp(-1.5/T4) / (1 + (ne/ncrit) * T4**-0.5) * np.exp(-14388 / wavelength / T4) * T4**-0.5 * np.exp(-14388 / wavelength / T4)

            wavelength = 6584
            ncrit = 1e3
            j_NII = nH*hii*ne / (1 + (ne/ncrit) * T4**-0.5) * np.exp(-14388 / wavelength/T4) * T4**-0.5
            
            wavelength = 6716
            ncrit = 1e3
            j_SII = nH*hii*ne  * np.exp(-1/T4)/ (1 + (ne/ncrit) * T4**-0.5) * np.exp(-14388 / wavelength/T4) * T4**-0.5

            pc_to_cm = 3.08e18
            msun_to_g = 2e33
            

            lum = np.c_[j_B_Ha,j_OIII,j_SII] * pc_to_cm**3 *  (snapdata["PartType0/Masses"]/rho)[:,None]
            lum_sum = np.sum(lum,axis=0); full_sum = np.sum(lum_sum)
            if full_sum:
                if lum_sum[0]: lum[:,0] *= full_sum/lum_sum[0]
                if lum_sum[1]: lum[:,1] *= full_sum/lum_sum[1]
                if lum_sum[2]: lum[:,2] *= full_sum/lum_sum[2]
            

            def get_color_matrix(rot):
                a = np.eye(3)
                a = rgb2hsv(a)
                a[:,0] += rot
                a[:,0] = a[:,0]%1
                a = hsv2rgb(a)
                return a

            color_matrix = get_color_matrix(-0.3)
            lum = lum @ color_matrix
            kappa = np.array([200,300,400]) * np.ones_like(self.mass)[:,None] / (pc_to_cm**2 / msun_to_g)
            if self.params["camera_distance"] < np.inf:
                lum[:] /= self.r[:,None]**2 # have to convert here because smoothing lengths are now in angular units
                lum[self.pos[:,2]<0] = 0  # ignore stuff behind the camera 
                kappa[self.pos[:,2]<0] = 0
                
            self.maps["SHO_RGB"] = GridRadTransfer(np.copy(lum), np.copy(self.mass), np.copy(kappa),  np.copy(self.pos), np.copy(self.hsml), self.params["res"], 2*self.params["rmax"]).swapaxes(0,1)
            np.savez_compressed(self.map_files["SHO_RGB"], SHO_RGB=self.maps["SHO_RGB"])

    def MakeImages(self,snapdata):
        #sigmoid = lambda x: x/np.sqrt(1+x*x) # tapering function to soften the saturation
        sigmoid = lambda x: x/(1+x) # tapering function to soften the saturation
        ha_map = np.copy(self.maps["SHO_RGB"])
        
        if self.params["SHO_RGB_norm"] == 0: 
            norm = [np.percentile(ha_map[:,:,2],99), np.percentile(ha_map[:,:,0],99), np.percentile(ha_map[:,:,1],99) ]
        else:
            norm = [self.params["SHO_RGB_norm"], self.params["SHO_RGB_norm"], self.params["SHO_RGB_norm"]]
        print("Using SHO_RGB normalizations %g %g %g"%(norm[0],norm[1],norm[2]))
        self.maps["SHO_RGB"][:,:,0] = sigmoid(ha_map[:,:,2]/norm[0])
        self.maps["SHO_RGB"][:,:,1] = sigmoid(ha_map[:,:,0]/norm[1])
        self.maps["SHO_RGB"][:,:,2] = sigmoid(ha_map[:,:,1]/norm[2])
    
        plt.imsave(self.params["filename_incomplete"], self.maps["SHO_RGB"][::-1]) # NOTE - we invert this to get the coordinate system right
        super().MakeImages(snapdata)
        #self.AddStarsToImage(snapdata)        
#        self.AddSizeScaleToImage()
#        self.AddTimestampToImage()
        
