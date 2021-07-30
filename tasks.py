import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Meshoid import GridSurfaceDensity
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

class SurfaceDensityMap:
    def __init__(self, params):
        
        
        self.RequiredSnapdata = ["PartType0/Coordinates", "PartType0/Masses", "PartType0/SmoothingLength"]
        self.params = params
        default_params = {"res": 128,
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

    def DoTask(self, snapdata):        
        if self.params["center"] is None:
            center = np.repeat(snapdata["Header"]["BoxSize"]*0.5,3)
        else: center = self.params["center"]
        if self.params["rmax"] is None:
            rmax = snapdata["Header"]["BoxSize"]/10
        else: rmax = self.params["rmax"]
        if self.params["pan"] is None:
            pan = 0
        else: pan = self.params["pan"]
        if self.params["tilt"] is None:
            tilt = 0
        else:
            tilt = self.params["tilt"]
            
        res = self.params["res"]
        pos, mass, hsml = np.copy(snapdata["PartType0/Coordinates"]), np.copy(snapdata["PartType0/Masses"]), np.copy(snapdata["PartType0/SmoothingLength"]) # copy these because we don't want to modify them
        hsml = np.clip(hsml,2*rmax/res, 1e100)
        pos -= center
        sigma_gas = GridSurfaceDensity(mass, pos, hsml, np.zeros(3), 2*rmax, res=res).T        
        fig, ax = plt.subplots(figsize=(6,6))
        X = Y = np.linspace(-rmax, rmax, res)
        X, Y = np.meshgrid(X, Y)
        p = ax.pcolormesh(X, Y, sigma_gas, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e3))
        ax.set_aspect('equal')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.0)
        fig.colorbar(p,label=r"$\Sigma_{\rm gas}$ $(\rm M_\odot\,pc^{-2})$",cax=cax)
        ax.set_xlabel("X (pc)")
        ax.set_ylabel("Y (pc)")
        plt.savefig("sigma_gas_%g.png"%snapdata["Header"]["Time"],bbox_inches='tight',dpi=300)
        plt.close()
        print("done")
        return sigma_gas    

