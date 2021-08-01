import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LightSource
from Meshoid import GridSurfaceDensityMultigrid as GridSurfaceDensity
import aggdraw
from PIL import Image, ImageDraw, ImageFont, ImageChops
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

class Task:
    """Class containing generic routines common to all tasks, and assigns default (null/empty) attributes that any task should have"""
    def __init__(self,params):
        self.RequiredSnapdata = []
        self.params = params

    def GetRequiredSnapdata(self):
        return self.RequiredSnapdata

    def AssignDefaultParams(self):
        for p in self.default_params.keys():
            if not p in self.params.keys(): self.params[p] = self.default_params[p]

    def AssignDefaultParamsFromSnapdata(self,snapdata):
        return

    

class SinkVis(Task):
    def __init__(self, params):
        """Class containing methods for coordinate transformations, rendering, etc. for a generic SinkVis-type map plot"""
        super().__init__(params)
        
        self.RequiredSnapdata = ["PartType0/Coordinates", "PartType0/Masses", "PartType0/SmoothingLength","PartType0/ParticleIDGenerationNumber","PartType0/ParticleChildIDsNumber", "PartType0/ParticleIDs", "PartType5/Coordinates","PartType5/Masses","PartType5/ParticleIDs", "PartType5/BH_Mass"]

        self.default_params = {"res": 512,
                          "rmax": None,
                          "limits": [1,2e3],
                          "center": None,
                          "pan": 0,
                          "tilt": 0,
                          "no_timestamp": False,
                          "no_size_scale": False,
                          "filename": None,
                          "sink_scale": 1,
                          "cmap": "inferno",
                           "backend": "PIL",
                          "rescale_hsml": False
                          }

        super().AssignDefaultParams()

    def CoordinateTransform(self,x,m=None,h=None):
        tilt, pan = self.params["tilt"], self.params["pan"]
        # first pan
        cosphi, sinphi = np.cos(np.pi*pan/180), np.sin(np.pi*pan/180)
        x[:] = np.c_[cosphi*x[:,0] + sinphi*x[:,2],x[:,1], -sinphi*x[:,0] + cosphi*x[:,2]]
        # then tilt
        costheta, sintheta = np.cos(np.pi*tilt/180), np.sin(np.pi*tilt/180)
        x[:] = np.c_[x[:,0], costheta*x[:,1] + sintheta*x[:,2], -sintheta*x[:,1] + costheta*x[:,2]]
        # then do projection if desired
        
        

    def SetupCoordsAndWeights(self, snapdata):
        res = self.params["res"]
        self.pos, self.mass, self.hsml = np.copy(snapdata["PartType0/Coordinates"]), np.copy(snapdata["PartType0/Masses"]), np.copy(snapdata["PartType0/SmoothingLength"]) # copy these because we don't want to modify them
        self.pos -= self.params["center"]
        self.CoordinateTransform(self.pos,self.mass,self.hsml)
        if self.params["rescale_hsml"]: self.hsml *= self.params["rescale_hsml"]
        self.hsml = np.clip(self.hsml,2*self.params["rmax"]/res, 1e100)

    def GenerateMaps(self,snapdata):
        self.maps = {}

    def SaveImage(self):
        if self.params["backend"] == "matplotlib":
            rmax = self.params["rmax"]
            self.ax.set(xlim=[-rmax,rmax],ylim=[-rmax,rmax])
            plt.savefig(self.params["filename"],bbox_inches='tight',dpi=200)            
            plt.close()

    def MakeImages(self,snapdata):
        self.AddStarsToImage(snapdata)
        self.AddSizeScaleToImage()
        self.AddTimestampToImage()
        self.SaveImage()

    def AddTimestampToImage(self):
        if self.params["no_timestamp"]: return
        fname = self.params["filename"]
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
        if self.params["backend"]=="matplotlib": return # matplotlib will have axis ticks for scale
        pc_to_AU = 206265.0
        if self.params["no_size_scale"]: return
        fname = self.params["filename"]
        F = Image.open(fname)
        draw = ImageDraw.Draw(F)
        gridres = self.params["res"]
        font = ImageFont.truetype("LiberationSans-Regular.ttf", gridres//12)
        r = self.params["rmax"]
        if (r>1000):
            scale_kpc=10**np.round(np.log10(r*0.5/1000))
            size_scale_text="%3.3gkpc"%(scale_kpc)
            size_scale_ending=gridres/16+gridres*(scale_kpc*1000)/(2*r)
        if (r>1e-2):
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
        if not "PartType5/Coordinates" in snapdata.keys(): return
        X_star = snapdata["PartType5/Coordinates"] - self.params["center"]
        m_star = snapdata["PartType5/BH_Mass"]
        self.CoordinateTransform(X_star, np.ones(len(X_star)), np.ones(len(X_star)))        
        if self.params["backend"]=="PIL":
            fname = self.params["filename"]
            F = Image.open(fname)
            gridres = F.size[0]
            draw = ImageDraw.Draw(F)
            d = aggdraw.Draw(F)
            pen = aggdraw.Pen('white',1) #gridres/800
            sink_relscale = 0.0025
            X_star ,m_star = X_star[m_star.argsort()[::-1]], np.sort(m_star)[::-1]
            for j in np.arange(len(X_star))[m_star>1e-2]:
                X = X_star[j]
                ms = m_star[j]
                star_size = gridres* sink_relscale * (np.log10(ms/self.params["sink_scale"]) + 1)
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
#            print(star_size)
            colors = np.array([self.GetStarColor(m) for m in m_star])/255
 
            self.ax.scatter(X_star[:,0], X_star[:,1],s=star_size*5,edgecolor='white',lw=0.1,facecolor=colors,marker='*')
            

    def GetStarColor(self, mass_in_msun):
        if self.params["cmap"] in ('afmhot', 'inferno', "Blues"):
            star_colors = np.array([[255, 100, 60],[120, 200, 150],[75, 80, 255]]) #alternate colors, red-green-blue, easier to see on a bright color map
        else:
            star_colors = np.array([[255, 203, 132],[255, 243, 233],[155, 176, 255]]) #default colors, reddish for small ones, yellow-white for mid sized and blue for large
        colors = np.int_([np.interp(np.log10(mass_in_msun),[-1,0,1],star_colors[:,i]) for i in range(3)])
        return (colors[0],colors[1],colors[2])# if len(colors)==1 else colors)

    def AssignDefaultParamsFromSnapdata(self,snapdata):
        if self.params["center"] is None:
            self.params["center"] = np.repeat(snapdata["Header"]["BoxSize"]*0.5,3)
            center = self.params["center"]
        else: center = self.params["center"]
        if self.params["rmax"] is None:
            self.params["rmax"] = snapdata["Header"]["BoxSize"]/10
        if self.params["filename"] is None:
            self.params["filename"] = "sigma_gas_%s_%s.png"%(str(round(self.params["Time"]/1e-6)).zfill(4), str(round(self.params["pan"])).zfill(4))

    
    def DoTask(self, snapdata):
        self.AssignDefaultParamsFromSnapdata(snapdata)
        self.SetupCoordsAndWeights(snapdata)        
        self.GenerateMaps(snapdata)
        self.MakeImages(snapdata)
        return self.maps


class SinkVisSigmaGas(SinkVis):
    def GenerateMaps(self,snapdata):
        super().GenerateMaps(snapdata)
        self.maps["sigma_gas"] = GridSurfaceDensity(self.mass, self.pos, self.hsml, np.zeros(3), 2*self.params["rmax"], res=self.params["res"]).T

    def MakeImages(self,snapdata):
        vmin, vmax = self.params["limits"]
                              
        f = (np.log10(self.maps["sigma_gas"])-np.log10(vmin))/(np.log10(vmax)-np.log10(vmin))

        if self.params["backend"]=="PIL":
            plt.imsave(self.params["filename"], plt.get_cmap(self.params["cmap"])(np.flipud(f))) # NOTE - we invert this to get the coordinate system right
        elif self.params["backend"]=="matplotlib":
            self.fig, self.ax = plt.subplots(figsize=(4,4))
            X = Y = np.linspace(-self.params["rmax"], self.params["rmax"], self.params["res"])
            X, Y = np.meshgrid(X, Y)
            p = self.ax.pcolormesh(X, Y, self.maps["sigma_gas"], norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e3),cmap=self.params["cmap"])
            self.ax.set_aspect('equal')
 
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.0)
            self.fig.colorbar(p,label=r"$\Sigma_{\rm gas}$ $(\rm M_\odot\,pc^{-2})$",cax=cax)
            self.ax.set_xlabel("X (pc)")
            self.ax.set_ylabel("Y (pc)")

        super().MakeImages(snapdata)


class SinkVisCoolMap(SinkVis):
    def __init__(self,params):
        super().__init__(params)
        self.RequiredSnapdata.append("PartType0/Velocities")
        self.default_params["cool_cmap"] = 'magma'
        super().AssignDefaultParams()
        
    def GenerateMaps(self,snapdata):
        super().GenerateMaps(snapdata)
        # need to apply coordinate transforms to z-velocity
        v = np.copy(snapdata["PartType0/Velocities"])
        self.CoordinateTransform(v)
        sigma_gas = GridSurfaceDensity(self.mass, self.pos, self.hsml, np.zeros(3), 2*self.params["rmax"], res=self.params["res"]).T
        sigma_1D = GridSurfaceDensity(self.mass * v[:,2]**2, self.pos, self.hsml, np.zeros(3), 2*self.params["rmax"], res=self.params["res"]).T/sigma_gas
        v_avg = GridSurfaceDensity(self.mass * v[:,2], self.pos, self.hsml, np.zeros(3), 2*self.params["rmax"], res=self.params["res"]).T/sigma_gas
        sigma_1D = np.sqrt(sigma_1D - v_avg**2)/1e3
        fgas = (np.log10(sigma_gas)-np.log10(self.params["limits"][0]))/np.log10(self.params["limits"][1]/self.params["limits"][0])
        fgas = np.clip(fgas,0,1)
        ls = LightSource(azdeg=315, altdeg=45)
        #lightness = ls.hillshade(z, vert_exag=4)
        mapcolor = plt.get_cmap(self.params["cool_cmap"])(np.log10(sigma_1D/0.1)/2)
        cool_data = ls.blend_hsv(mapcolor[:,:,:3], fgas[:,:,None])
        self.maps["coolmap"] = cool_data

    def MakeImages(self,snapdata):
        vmin, vmax = self.params["limits"]
        if self.params["backend"]=="PIL":
            plt.imsave(self.params["filename"], self.maps["coolmap"]) # NOTE - we invert this to get the coordinate system right            
        
