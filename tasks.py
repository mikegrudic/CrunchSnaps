import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Meshoid import GridSurfaceDensity
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

class SimVis:
    def __init__(self, params):
        self.RequiredSnapdata = ["PartType0/Coordinates", "PartType0/Masses", "PartType0/SmoothingLength","PartType0/ParticleIDGenerationNumber","PartType0/ParticleChildIDsNumber"]
        self.params = params

        default_params = {"res": 1024,
                          "rmax": None,
                          "limits": None,
                          "center": None,
                          "pan": None,
                          "tilt": None
                          }

        for p in default_params.keys():
            if not p in params.keys(): self.params[p] = default_params[p]

    def GetRequiredSnapdata(self):
        return self.RequiredSnapdata            

    def CoordinateTransform(self, x, m, h):
        tilt, pan = self.params["tilt"], self.params["pan"]
#        print(tilt, pan)
        # first pan
        cosphi, sinphi = np.cos(np.pi*pan/180), np.sin(np.pi*pan/180)
        x[:] = np.c_[cosphi*x[:,0] + sinphi*x[:,2],x[:,1], -sinphi*x[:,0] + cosphi*x[:,2]]
        # then tilt
        costheta, sintheta = np.cos(np.pi*tilt/180), np.sin(np.pi*tilt/180)
        x[:] = np.c_[x[:,0], costheta*x[:,1] + sintheta*x[:,2], -sintheta*x[:,1] + costheta*x[:,2]]
        # then do projection
            
    
    def DoTask(self, snapdata):        
        if self.params["center"] is None:
            center = np.repeat(snapdata["Header"]["BoxSize"]*0.5,3)
        else: center = self.params["center"]
        if self.params["rmax"] is None:
            rmax = snapdata["Header"]["BoxSize"]/10
        else: rmax = self.params["rmax"]
        if self.params["pan"] is None:
            self.params["pan"] = 0
        if self.params["tilt"] is None:
            self.params["tilt"] = 0
            
        res = self.params["res"]
        pos, mass, hsml = np.copy(snapdata["PartType0/Coordinates"]), np.copy(snapdata["PartType0/Masses"]), np.copy(snapdata["PartType0/SmoothingLength"]) # copy these because we don't want to modify them
        pos -= center
        print(pos[0])
        self.CoordinateTransform(pos,mass,hsml)
        
        hsml = np.clip(hsml,2*rmax/res, 1e100)
        sigma_gas = GridSurfaceDensity(mass, pos, hsml, np.zeros(3), 2*rmax, res=res).T
        f = (np.log10(sigma_gas)-np.log10(1))/(np.log10(1e3)-np.log10(1))
        plt.imsave("sigma_gas_%s_%s.png"%(str(round(self.params["pan"])).zfill(4), str(round(self.params["tilt"])).zfill(4)), plt.get_cmap('viridis')(f))
        # fig, ax = plt.subplots(figsize=(6,6))
        # X = Y = np.linspace(-rmax, rmax, res)
        # X, Y = np.meshgrid(X, Y)
        # p = ax.pcolormesh(X, Y, sigma_gas, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e3))
        # ax.set_aspect('equal')

        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes("right", size="5%", pad=0.0)
        # fig.colorbar(p,label=r"$\Sigma_{\rm gas}$ $(\rm M_\odot\,pc^{-2})$",cax=cax)
        # ax.set_xlabel("X (pc)")
        # ax.set_ylabel("Y (pc)")
        # plt.savefig("sigma_gas_%s_%s.png"%(str(round(self.params["pan"])).zfill(4), str(round(self.params["tilt"])).zfill(4)),bbox_inches='tight',dpi=300)
        # plt.close()
#        print("done")
        return sigma_gas    

class SimVis2:
    def __init__(self, params):
        self.RequiredSnapdata = ["PartType0/Coordinates", "PartType0/Masses", "PartType0/SmoothingLength","PartType0/ParticleIDGenerationNumber","PartType0/ParticleChildIDsNumber"]
        self.params = params

        default_params = {"res": 256,
                          "rmax": None,
                          "limits": None,
                          "center": None,
                          "pan": None,
                          "tilt": None
                          }

        for p in default_params.keys():
            if not p in params.keys(): self.params[p] = default_params[p]

    def GetRequiredSnapdata(self):
        return self.RequiredSnapdata            

    def CoordinateTransform(self, x, m, h):
        tilt, pan = self.params["tilt"], self.params["pan"]
#        print(tilt, pan)
        # first pan
        cosphi, sinphi = np.cos(np.pi*pan/180), np.sin(np.pi*pan/180)
        x[:] = np.c_[cosphi*x[:,0] + sinphi*x[:,2],x[:,1], -sinphi*x[:,0] + cosphi*x[:,2]]
        # then tilt
        costheta, sintheta = np.cos(np.pi*tilt/180), np.sin(np.pi*tilt/180)
        x[:] = np.c_[x[:,0], costheta*x[:,1] + sintheta*x[:,2], -sintheta*x[:,1] + costheta*x[:,2]]
        # then do projection
            
    
    def DoTask(self, snapdata):        
        if self.params["center"] is None:
            center = np.repeat(snapdata["Header"]["BoxSize"]*0.5,3)
        else: center = self.params["center"]
        if self.params["rmax"] is None:
            rmax = snapdata["Header"]["BoxSize"]/10
        else: rmax = self.params["rmax"]
        if self.params["pan"] is None:
            self.params["pan"] = 0
        if self.params["tilt"] is None:
            self.params["tilt"] = 0
            
        res = self.params["res"]
        pos, mass, hsml = np.copy(snapdata["PartType0/Coordinates"]), np.copy(snapdata["PartType0/Masses"]), np.copy(snapdata["PartType0/SmoothingLength"]) # copy these because we don't want to modify them
        pos -= center
        print(pos[0])
        self.CoordinateTransform(pos,mass,hsml)
        
        hsml = np.clip(hsml,2*rmax/res, 1e100)
        sigma_gas = GridSurfaceDensity(mass, pos, hsml, np.zeros(3), 2*rmax, res=res).T
        f = (np.log10(sigma_gas)-np.log10(1))/(np.log10(1e3)-np.log10(1))
        plt.imsave("sigma_gas2_%s_%s.png"%(str(round(self.params["pan"])).zfill(4), str(round(self.params["tilt"])).zfill(4)), plt.get_cmap('inferno')(f))
        # fig, ax = plt.subplots(figsize=(6,6))
        # X = Y = np.linspace(-rmax, rmax, res)
        # X, Y = np.meshgrid(X, Y)
        # p = ax.pcolormesh(X, Y, sigma_gas, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e3),cmap='inferno')
        # ax.set_aspect('equal')

        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes("right", size="5%", pad=0.0)
        # fig.colorbar(p,label=r"$\Sigma_{\rm gas}$ $(\rm M_\odot\,pc^{-2})$",cax=cax)
        # ax.set_xlabel("X (pc)")
        # ax.set_ylabel("Y (pc)")
        # plt.savefig("sigma_gas2_%s_%s.png"%(str(round(self.params["pan"])).zfill(4), str(round(self.params["tilt"])).zfill(4)),bbox_inches='tight',dpi=300)
        # plt.close()
#        print("done")
        return sigma_gas    
