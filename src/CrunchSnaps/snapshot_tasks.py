import numpy as np
import re
from astropy import constants as _ac, units as _au
from .misc_functions import *
from os.path import isfile
import json
import os



def _get_font(size):
    """Get a PIL ImageFont at the requested pixel size, with robust fallbacks."""
    import matplotlib.font_manager as fm
    from PIL import ImageFont

    # Try matplotlib's font manager to find a sans-serif system font
    try:
        font_path = fm.findfont(fm.FontProperties(family="sans-serif"))
        if font_path and os.path.isfile(font_path):
            return ImageFont.truetype(font_path, size)
    except Exception:
        pass
    # Try Pillow >= 10.1 load_default with size
    try:
        return ImageFont.load_default(size=size)
    except TypeError:
        pass
    # Last resort: unsized default bitmap font
    return ImageFont.load_default()


class Task:
    """Class containing generic routines common to all tasks, and assigns default (null/empty) attributes that any task should have"""

    def __init__(self, params):
        self.RequiredSnapdata = []
        self.params = params.copy()

    def GetRequiredSnapdata(self):
        return list(np.unique(self.RequiredSnapdata))

    def AssignDefaultParams(self):
        for k in self.default_params.keys():
            if not k in self.params.keys():
                self.params[k] = self.default_params[k]

    def AssignDefaultParamsFromSnapdata(self, snapdata):
        return

    def DoTask(self, snapdata):
        return


class SinkVis(Task):
    def __init__(self, params):
        """Class containing methods for coordinate transformations, rendering, etc. for a generic SinkVis-type map plot"""
        super().__init__(params)

        self.default_params = {
            "Time": 0,
            "res": 512,
            "rmax": None,
            "limits": None,  # [10,3e3],
            "center": None,
            "pan": 0,
            "tilt": 0,
            "no_timestamp": False,
            "no_size_scale": False,
            "filename": None,
            "sink_scale": 1,
            "cmap": "magma",
            "backend": "PIL",
            "rescale_hsml": False,
            "FOV": 90,
            "camera_distance": np.inf,
            "center_on_star": False,
            "center_on_ID": False,
            "center_on_densest": False,
            "realstars": True,
            "realstars_opacity": 1.0,
            "realstars_max_lum": 1e7,
            "realstars_lum_exp": 1.0,
            "realstars_background": 0,
            "threads": -1,
            "cubemap_dir": "forward",
            "camera_dir": None,
            "camera_right": None,
            "camera_up": None,
            "index": None,
            "no_stars": False,
            "overwrite": True,
            "unit_scalefac": 1,
            "outputfolder": ".",
            "SHO_RGB_norm": 0,
            "outflow_only": False,
            "no_colorbar": False,
        }

        self.AssignDefaultParams()
        basename, ext = os.path.splitext(self.params["filename"])
        self.params["filename_incomplete"] = basename + "_incomplete" + ext
        self.params_that_affect_maps = [
            "Time",
            "res",
            "rmax",
            "center",
            "pan",
            "tilt",
            "FOV",
            "camera_distance",
            "center",
            "center_on_star",
            "center_on_densest",
            "cubemap_dir",
            "camera_dir",
            "camera_right",
            "camera_up",
            "rescale_hsml",
            "res",
            "realstars_max_lum",
            "realstars_lum_exp",
            "realstars_opacity",
        ]
        dump = {}
        for k in self.params_that_affect_maps:
            if type(self.params[k]) == np.ndarray:
                dump[k] = self.params[k].tolist()
            else:
                dump[k] = self.params[k]
        self.params_hash = str(hash(json.dumps(dump, sort_keys=True)))

        mapdir = self.params["outputfolder"] + "/.maps"
        while not os.path.isdir(mapdir):
            try:
                os.mkdir(mapdir)
            except Exception as e:
                print(f"Exception when loading map: {e}")
                continue

        if self.params["realstars"]:
            self.required_maps.add("realstars")

        # filename for the saved maps will by MAPNAME_(hash # of input params)
        self.map_files = dict([(m, mapdir + "/" + m + "_" + self.params_hash) for m in self.required_maps])
        self.maps = {}
        if self.params["threads"] != 1:
            self.parallel = True
            # if negative, just use all available threads, otherwise set to desired value
            if self.params["threads"] > 0:
                try:
                    from numba import set_num_threads as _numba_set_num_threads
                    _numba_set_num_threads(self.params["threads"])
                except (ImportError, RuntimeError):
                    pass  # Numba pool already launched with fewer threads
        else:
            self.parallel = False

        if isfile(self.params["filename"]) and not self.params["overwrite"]:
            self.RequiredSnapdata = []
            print("%s already exists, skipping figure..." % (self.params["filename"]))
            self.TaskDone = True
        else:
            self.TaskDone = False
            self.DetermineRequiredSnapdata()

    def DetermineRequiredSnapdata(self):
        self.RequiredSnapdata = [
            "PartType5/Coordinates",
            "PartType5/Masses",
            "PartType5/ParticleIDs",
            "PartType5/Sink_Mass",
        ]

        if self.params["realstars"]:
            self.RequiredSnapdata += [
                "PartType5/ProtoStellarRadius_inSolar",
                "PartType5/StarLuminosity_Solar",
                "PartType0/Masses",
                "PartType0/Coordinates",
                "PartType0/KernelMaxRadius",
            ]
        if any((self.params[k] for k in ("center", "center_on_star", "center_on_ID", "center_on_densest"))):
            self.RequiredSnapdata += ["PartType0/Density", "PartType0/Coordinates"]
        if self.params["outflow_only"]:
            self.RequiredSnapdata += ["PartType0/Velocities", "PartType5/Velocities"]
        # check if we have maps already saved
        self.render_maps = False
        for mapname in self.required_maps:
            if not isfile(self.map_files[mapname] + ".npz"):
                self.render_maps = True
            else:
                print("Loading %s map from %s" % (mapname, self.map_files[mapname]))
                try:
                    self.maps[mapname] = np.load(self.map_files[mapname] + ".npz")[mapname]
                except Exception as e:
                    print(
                        "Error when loading file %s! Removing potentially corrupted file..." % self.map_files[mapname]
                    )
                    os.remove(self.map_files[mapname] + ".npz")
                    self.render_maps = True

    def AssignDefaultParams(self):
        super().AssignDefaultParams()
        if self.params["index"] is None:
            self.params["index"] = round(self.params["Time"] / 1e-6)
        self.params["filename_suffix"] = "%s_%s_%s.png" % (
            str(self.params["index"]).zfill(5),
            str(round(self.params["pan"] * 10)).zfill(4),
            self.params["cubemap_dir"],
        )

    def DoCoordinateTransform(self, x, m=None, h=None, contravariant=False, update_r=True):
        # center on the designated center coordinate
        if not contravariant:
            x[:] -= self.params["center"]

        # without a specified camera direction, we just use a simple tilt/pan scheme
        if self.params["camera_dir"] is None:
            tilt, pan = self.params["tilt"], self.params["pan"]
            if contravariant:
                tilt, pan = -tilt, -pan
            # first pan
            cosphi, sinphi = np.cos(np.pi * pan / 180), np.sin(np.pi * pan / 180)
            x[:] = np.stack([cosphi * x[:, 0] + sinphi * x[:, 2], x[:, 1], -sinphi * x[:, 0] + cosphi * x[:, 2]], 1)
            # then tilt
            costheta, sintheta = np.cos(np.pi * tilt / 180), np.sin(np.pi * tilt / 180)
            x[:] = np.stack(
                [x[:, 0], costheta * x[:, 1] + sintheta * x[:, 2], -sintheta * x[:, 1] + costheta * x[:, 2]], 1
            )
        else:  # we have a camera position and coordinate basis
            if contravariant:
                x[:] = (self.camera_matrix_vectors @ x.T).T  # note that @ performs matrix multiplication
            else:
                x[:] = (self.camera_matrix @ x.T).T

        if self.params["camera_distance"] != np.inf and not contravariant:
            # transform so camera is at z=0:
            x[:, 2] = x[:, 2] - self.params["camera_distance"]

        # shuffle the axes to get the desired cubemap direction
        cubedir = self.params["cubemap_dir"]
        if cubedir != "forward":
            if cubedir == "right":
                x[:] = np.c_[-x[:, 2], x[:, 1], x[:, 0]]
            elif cubedir == "left":
                x[:] = np.c_[x[:, 2], x[:, 1], -x[:, 0]]
            elif cubedir == "up":
                x[:] = np.c_[x[:, 0], -x[:, 2], x[:, 1]]
            elif cubedir == "down":
                x[:] = np.c_[x[:, 0], x[:, 2], -x[:, 1]]
            elif cubedir == "backward":
                x[:] = np.c_[-x[:, 0], x[:, 1], -x[:, 2]]

        # then do projection if desired
        if self.params["camera_distance"] != np.inf:
            if not contravariant:
                # now transform from 3D to angular system
                with np.errstate(divide="ignore", invalid="ignore"):
                    r = np.abs(x[:, 2])
                    x[:, :2] = x[:, :2] / (-x[:, 2][:, None])  # homogeneous coordinates
                    if update_r:
                        self.r = r
                    behind = x[:, 2] >= 0  # at or behind the camera
                    if h is not None:
                        h[:] = h / r  # kernel lengths are now angular (divide by distance)
                        h[behind] = 0
                    if m is not None:
                        m[:] /= r**2  # rescale mass weights so that integrated surface density remains the same
                        m[behind] = 0

            else:  # dealing with a contravariant vector such as velocity - want the [:,2] component to correspond to line-of-sight value
                # this would have been converted to angular by now - let's convert back to real space
                global_coords = np.copy(self.pos)
                # multiply by z, now we're in the rotated real space frame
                global_coords[:, :2] *= -global_coords[:, 2][:, None]
                # get the radial component
                x[:, 2] = np.sum(x * global_coords, axis=1) / np.sum(global_coords**2, axis=1) ** 0.5

    def SetupCoordsAndWeights(self, snapdata):
        from meshoid import Meshoid

        res = self.params["res"]
        if "PartType0/Coordinates" in snapdata.keys():
            if "PartType0/KernelMaxRadius" not in snapdata:
                snapdata["PartType0/KernelMaxRadius"] = Meshoid(
                    snapdata["PartType0/Coordinates"], boxsize=snapdata["Header"]["BoxSize"]
                ).SmoothingLength()
            self.pos, self.mass, self.hsml = (
                np.copy(snapdata["PartType0/Coordinates"]),
                np.copy(snapdata["PartType0/Masses"]),
                np.copy(snapdata["PartType0/KernelMaxRadius"]),
            )  # copy these because we don't want to modify them
            if self.params["rescale_hsml"]:
                self.hsml *= self.params["rescale_hsml"]

        if self.params["outflow_only"]:
            if "PartType5/Coordinates" in snapdata.keys():
                from scipy.spatial import KDTree
                # find nearest star to each gas cell
                _, ngb = KDTree(snapdata["PartType5/Coordinates"]).query(self.pos)
                dx = self.pos - snapdata["PartType5/Coordinates"][ngb]
                dv = snapdata["PartType0/Velocities"] - snapdata["PartType5/Velocities"][ngb]
                self.mass *= (dx * dv).sum(1) > 0

        # Setting up coordinate basis
        if self.params["camera_dir"] is not None:
            self.camera_dir = self.params["camera_dir"]
            NormalizeVector(self.params["camera_dir"])
            if self.params["camera_up"] is None:
                self.camera_up = np.array(
                    [0, 1.0, 0]
                )  # default "up" direction is +y, we will project it out if the camera is tilted
            else:
                self.camera_up = self.params["camera_up"]

            # if we've specified an up direction, project out the component parallel to the forward direction and normalize
            self.camera_up -= sum(self.camera_dir * self.camera_up).sum() * self.camera_dir
            NormalizeVector(self.camera_up)
            # now get the "right" vector as the cross product of forward x up. this will be normalized to machine precision
            self.camera_right = np.cross(self.camera_up, self.camera_dir)

            # matrix of coordinate vectors - operate this on coordinates to apply transformation - operates on COORDINATES not vectors
            self.camera_matrix = np.c_[self.camera_right, self.camera_up, self.camera_dir].T
            # since vector fields are contravariant, this is the operator for transforming v and B (note that this is an orthogonal matrix so the transpose is the inverse)
            self.camera_matrix_vectors = self.camera_matrix.T

        if "PartType0/Coordinates" in snapdata.keys():
            self.DoCoordinateTransform(self.pos, self.mass, self.hsml)
            # cull particles that can't contribute to the grid (behind camera, outside FOV)
            rmax = self.params["rmax"]
            keep = (self.hsml > 0) & (self.mass != 0)
            keep &= (np.abs(self.pos[:, 0]) - self.hsml < rmax) & (np.abs(self.pos[:, 1]) - self.hsml < rmax)
            if keep.sum() < len(keep):
                self.pos = self.pos[keep]
                self.mass = self.mass[keep]
                self.hsml = self.hsml[keep]
            self._keep_mask = keep  # save for subclasses that need to cull field arrays
            self.hsml = np.clip(self.hsml, 2 * rmax / res, np.inf)

    def GenerateMaps(self, snapdata):
        return

    def SaveImage(self):
        print("Saving ", self.params["filename"])
        if self.params["backend"] == "matplotlib":
            from matplotlib import pyplot as plt
            rmax = self.params["rmax"]
            self.ax.set(xlim=[-rmax, rmax], ylim=[-rmax, rmax])
            plt.savefig(self.params["filename_incomplete"], bbox_inches="tight", dpi=400)
            plt.close()
        elif self.params["backend"] == "PIL":
            self._pil_image.convert("RGB").save(self.params["filename_incomplete"])
            self._pil_image = None
        os.rename(self.params["filename_incomplete"], self.params["filename"])

    def MakeImages(self, snapdata):
        if not self.params["no_stars"]:
            self.AddStarsToImage(snapdata)
        self.AddSizeScaleToImage(snapdata["Header"])
        self.AddTimestampToImage(snapdata["Header"])
        self.SaveImage()

    def AddTimestampToImage(self, header):
        if self.params["no_timestamp"]:
            return
        if "UnitLength_In_CGS" in header.keys():
            unit_time_in_s = header["UnitLength_In_CGS"] / header["UnitVelocity_In_CGS"]
        else:
            unit_time_in_s = _au.kpc.to(_au.cm) / (_au.km.to(_au.cm))  # kpc/(km/s) in seconds
        time_Myr = self.params["Time"] * unit_time_in_s / _au.Myr.to(_au.s)
        if time_Myr >= 1e-2:
            time_text = "%3.2gMyr" % (time_Myr)
        elif time_Myr >= 1e-4:
            time_text = "%3.2gkyr" % (time_Myr * 1e3)
        else:
            time_text = "%3.2gyr" % (time_Myr * 1e6)

        if self.params["backend"] == "PIL":
            from PIL import ImageDraw
            gridres = self._pil_image.size[0]
            draw = ImageDraw.Draw(self._pil_image)
            font = _get_font(gridres // 12)
            draw.text((gridres / 16, gridres / 24), time_text, font=font)
        elif self.params["backend"] == "matplotlib":
            self.ax.text(-self.params["rmax"] * 0.85, self.params["rmax"] * 0.85, time_text, color="#FFFFFF")

    def AddSizeScaleToImage(self, header):
        #        if self.params["camera_distance"] < np.inf: return
        if self.params["backend"] == "matplotlib":
            return  # matplotlib will have axis ticks for scale
        pc_to_AU = _au.pc.to(_au.AU)
        if self.params["no_size_scale"]:
            return
        from PIL import ImageDraw
        draw = ImageDraw.Draw(self._pil_image)
        gridres = self.params["res"]
        font = _get_font(gridres // 12)
        if "UnitLength_In_CGS" in header.keys():
            r = self.params["rmax"] * header["UnitLength_In_CGS"] / _au.pc.to(_au.cm)  # in pc
        else:
            r = self.params["rmax"]
        if self.params["camera_distance"] < np.inf:
            r = self.params["rmax"] * self.params["camera_distance"]  # * self.params["unit_scalefac"]
        if r * 2 > 1000:
            scale_kpc = 10 ** np.round(np.log10(r * 0.5 / 1000))
            size_scale_text = "%3.3gkpc" % (scale_kpc)
            size_scale_ending = gridres / 16 + gridres * (scale_kpc * 1000) / (2 * r)
        elif r > 1e-2:
            scale_pc = 10 ** np.round(np.log10(r * 0.5))
            size_scale_text = "%3.3gpc" % (scale_pc)
            size_scale_ending = gridres / 16 + gridres * (scale_pc) / (2 * r)
        else:
            scale_AU = 10 ** np.round(np.log10(r * 0.5 * pc_to_AU))
            size_scale_text = "%3.4gAU" % (scale_AU)
            size_scale_ending = gridres / 16 + gridres * (scale_AU) / (2 * r * pc_to_AU)
        draw.line(((gridres / 16, 7 * gridres / 8), (size_scale_ending, 7 * gridres / 8)), fill="#FFFFFF", width=6)
        draw.text((gridres / 16, 7 * gridres / 8 + 5), size_scale_text, font=font)

    def AddStarsToImage(self, snapdata):
        if "PartType5/Coordinates" not in snapdata.keys():
            return
        if not len(snapdata["PartType5/Coordinates"]):
            return
        X_star = np.copy(snapdata["PartType5/Coordinates"])
        m_star = snapdata["PartType5/Sink_Mass"]

        if self.params["backend"] == "PIL":
            from PIL import Image
            if self.params["realstars"]:  # use realstars for stellar images
                if "realstars" not in self.maps:
                    from .realistic_stars import make_stars_image
                    self.maps["realstars"] = make_stars_image(
                        self,
                        snapdata,
                        lum_max_solar=self.params["realstars_max_lum"],
                        lum_renorm_exponent=self.params["realstars_lum_exp"],
                        IMG_RES=self.params["res"],
                        IMG_SIZE=2 * self.params["rmax"],
                        opacity_scalefac=self.params["realstars_opacity"],
                        I_background=self.params["realstars_background"],
                        threads=self.params["threads"],
                    )
                np.savez_compressed(self.map_files["realstars"], realstars=self.maps["realstars"])
                # blend realstars into the in-memory image
                img_arr = np.array(self._pil_image.convert("RGB")).astype(np.float32) / 255
                img_arr = np.clip(img_arr + self.maps["realstars"], 0, 1)
                self._pil_image = Image.fromarray((img_arr * 255).astype(np.uint8), "RGB").convert("RGBA")

            else:  # use derpy PIL circles
                import aggdraw
                self.DoCoordinateTransform(X_star, np.ones(len(X_star)), np.ones(len(X_star)))
                gridres = self._pil_image.size[0]
                d = aggdraw.Draw(self._pil_image)
                pen = aggdraw.Pen(self.Star_Edge_Color(), 1)  # gridres/800
                sink_relscale = 0.0025
                X_star, m_star = X_star[m_star.argsort()[::-1]], np.sort(m_star)[::-1]

                for j in np.arange(len(X_star)):
                    X = X_star[j]
                    ms = m_star[j]
                    if ms == 0:
                        continue
                    star_size = max(1, gridres * sink_relscale * (np.log10(ms / self.params["sink_scale"]) + 1))
                    if self.params["camera_distance"] < np.inf:
                        if X[2] > 0:
                            continue
                    star_size = max(1, star_size)
                    p = aggdraw.Brush(self.GetStarColor(ms))
                    norm_coords = (X[:2] + self.params["rmax"]) / (2 * self.params["rmax"]) * gridres
                    # Pillow puts the origin in the top left corner, so we need to flip the y axis
                    norm_coords[1] = gridres - norm_coords[1]
                    coords = np.concatenate([norm_coords - star_size, norm_coords + star_size])
                    d.ellipse(coords, pen, p)
                d.flush()
        elif self.params["backend"] == "matplotlib":
            star_size = np.log10(m_star / self.params["sink_scale"]) + 2
            star_size[m_star < self.params["sink_scale"]] = 0
            colors = np.array([self.GetStarColor(m) for m in m_star]) / 255

            self.ax.scatter(
                X_star[:, 0],
                X_star[:, 1],
                s=star_size * 5,
                edgecolor=self.Star_Edge_Color(),
                lw=0.1,
                facecolor=colors,
                marker="*",
            )

    def Star_Edge_Color(self):
        if self.params["cmap"] in ("afmhot", "inferno", "Blues"):
            return "black"
        else:
            return "white"

    def GetStarColor(self, mass_in_msun):
        if self.params["cmap"] in ("afmhot", "inferno", "Blues"):
            # alternate colors, red-green-blue, easier to see on a bright color map
            star_colors = np.array([[255, 100, 60], [120, 200, 150], [75, 80, 255]])
        else:
            # default colors, reddish for small ones, yellow-white for mid sized and blue for large
            star_colors = np.array([[255, 203, 132], [255, 243, 233], [155, 176, 255]])
        if mass_in_msun > 1e4:  # assume a black hole
            colors = [np.zeros_like(mass_in_msun) for i in range(3)]
        else:
            colors = np.int_([np.interp(np.log10(mass_in_msun), [-1, 0, 1], star_colors[:, i]) for i in range(3)])
        return (colors[0], colors[1], colors[2])  # if len(colors)==1 else colors)

    def assign_center(self, snapdata):
        """Assign the center of the image"""

        if isinstance(self.params["center"], np.ndarray):
            if len(self.params["center"]) == 3:
                return

        match self.params["center"]:
            case "densest":
                rho = snapdata["PartType0/Masses"] / snapdata["PartType0/KernelMaxRadius"] ** 3
                self.params["center"] = snapdata["PartType0/Coordinates"][rho.argmax()]
            case "median":
                self.params["center"] = np.median(snapdata["PartType0/Coordinates"], axis=0)
            case x if "massive" in x:
                if "PartType5/Coordinates" in snapdata.keys():
                    if "=" in self.params["center"]:
                        num = int(self.params["center"].split("=")[1])
                    else:
                        num = 1
                    self.params["center"] = snapdata["PartType5/Coordinates"][
                        snapdata["PartType5/Sink_Mass"].argsort()[::-1]
                    ][num - 1]
                else:
                    rho = snapdata["PartType0/Masses"] / snapdata["PartType0/KernelMaxRadius"] ** 3
                    self.params["center"] = snapdata["PartType0/Coordinates"][rho.argmax()]
            case x if "ID" in x:
                if "=" in self.params["center"]:
                    id = int(self.params["center"].split("=")[1])
                else:
                    raise ValueError("Must specify particle ID for center as 'ID=value'")
                for k, data in snapdata.items():
                    if "IDs" in k:
                        ids = data
                        coords = snapdata[k.replace("ParticleIDs", "Coordinates")]
                        if id in ids:
                            self.params["center"] = coords[ids == id]

            case x if "," in x:
                self.params["center"] = np.array([float(s) for s in self.params["center"].split(",")])
        #            case _:

        if not isinstance(self.params["center"], np.ndarray):
            self.params["center"] = np.repeat(snapdata["Header"]["BoxSize"] * 0.5, 3)  # default

    def _compute_default_rmax(self, snapdata):
        """Compute default rmax from mass-weighted 2D variance in the viewing plane.

        Scaled so that a uniform-density cube returns rmax = BoxSize/2.
        """
        if "PartType0/Coordinates" not in snapdata or "PartType0/Masses" not in snapdata:
            return snapdata["Header"]["BoxSize"] / 2

        pos = snapdata["PartType0/Coordinates"]
        mass = snapdata["PartType0/Masses"]
        total_mass = mass.sum()
        if total_mass == 0:
            return snapdata["Header"]["BoxSize"] / 2

        dx = pos - self.params["center"]

        # line-of-sight direction in the original frame
        if self.params["camera_dir"] is not None:
            los = np.array(self.params["camera_dir"], dtype=float)
            los /= np.sqrt(np.dot(los, los))
        else:
            pan_rad = np.pi * self.params["pan"] / 180
            tilt_rad = np.pi * self.params["tilt"] / 180
            los = np.array([
                np.sin(pan_rad) * np.cos(tilt_rad),
                np.sin(tilt_rad),
                np.cos(pan_rad) * np.cos(tilt_rad),
            ])

        # Var_2D = Tr(Cov) - los^T @ Cov @ los, without forming the full 3x3
        trace_cov = np.sum(mass * np.sum(dx ** 2, axis=1)) / total_mass
        los_proj = dx @ los
        var_los = np.sum(mass * los_proj ** 2) / total_mass
        var_2d = max(trace_cov - var_los, 0)

        # For uniform cube of side L: Var_2D = L^2/6, want rmax = L/2
        # => rmax = sqrt(6)/2 * sqrt(Var_2D)
        return np.sqrt(6) / 2 * np.sqrt(var_2d)

    def AssignDefaultParamsFromSnapdata(self, snapdata):
        self.assign_center(snapdata)
        if self.params["rmax"] is None:
            if self.params["camera_distance"] < np.inf:
                self.params["rmax"] = self.params["FOV"] / 90  # angular width
            else:
                self.params["rmax"] = self._compute_default_rmax(snapdata)

    #            if self.params["camera_distance"] < np.inf and self.params["FOV"] is None:
    #                self.params["rmax"] /= self.params["camera_distance"] # convert to angular assuming rmax is real-space half-width at the focal distance
    def DoTask(self, snapdata):
        if self.TaskDone:
            return
        self.AssignDefaultParamsFromSnapdata(snapdata)
        if not self.has_required_maps():
            self.SetupCoordsAndWeights(snapdata)
            self.GenerateMaps(snapdata)
        self.MakeImages(snapdata)
        return self.maps

    @staticmethod
    def _render_latex_label(text, fontsize, color, rotation=0, dpi=200):
        """Render a LaTeX string via matplotlib and return as a RGBA PIL Image."""
        import io
        import matplotlib
        from matplotlib import pyplot as plt
        from PIL import Image

        matplotlib.use("Agg")
        tmp_fig = plt.figure(figsize=(0.01, 0.01), dpi=dpi)
        tmp_fig.text(
            0, 0, text,
            fontsize=fontsize, color=color,
            ha="left", va="bottom",
            rotation=rotation, rotation_mode="anchor",
        )
        buf = io.BytesIO()
        tmp_fig.savefig(buf, format="png", bbox_inches="tight", transparent=True, pad_inches=0.05, dpi=dpi)
        plt.close(tmp_fig)
        buf.seek(0)
        return Image.open(buf).convert("RGBA")

    def _add_colorbar_to_image(self, vmin, vmax, label=None):
        """Overlay a horizontal colorbar in the bottom-right of the in-memory PIL image."""
        from PIL import Image, ImageDraw
        from matplotlib import pyplot as plt

        img = self._pil_image
        W, H = img.size
        cmap = plt.get_cmap(self.params["cmap"])
        font_size = max(6, W // 80)
        tick_gap = 3

        # Compute tick values — integer powers of 10 between vmin and vmax,
        # dropping any that would overlap with the endpoint labels
        log_vmin, log_vmax = np.log10(vmin), np.log10(vmax)
        log_range = log_vmax - log_vmin
        tick_values = [vmin, vmax]
        tick_exp_min = int(np.ceil(log_vmin))
        tick_exp_max = int(np.floor(log_vmax))
        # Require at least 15% of the bar between any two ticks
        min_gap = 0.15
        for e in range(tick_exp_min, tick_exp_max + 1):
            tv = 10.0**e
            frac = (e - log_vmin) / log_range
            if frac > min_gap and frac < (1 - min_gap):
                tick_values.append(tv)
        tick_values.sort()

        # Format ticks
        def _format_tick(tv):
            exp = int(np.floor(np.log10(np.abs(tv))))
            coeff = tv / 10**exp
            if abs(coeff - 1.0) < 0.01:
                return r"$10^{%d}$" % exp
            return r"$%.2g\times10^{%d}$" % (coeff, exp)

        tick_labels_text = [_format_tick(tv) for tv in tick_values]

        # Measure tick labels to size the bar
        sample = self._render_latex_label(r"$10^{3}$", font_size, "#FFFFFF")
        tick_h = sample.size[1]
        max_tick_w = max(
            self._render_latex_label(t, font_size, "#FFFFFF").size[0]
            for t in tick_labels_text
        )

        # Horizontal bar geometry — bottom right
        bar_h = max(int(H * 0.02), 6)
        bar_w = int(W * 0.35)
        margin_r = max(int(W * 0.04), 8)
        margin_b = max(int(H * 0.04), 8)
        bar_x2 = W - margin_r
        bar_x1 = bar_x2 - bar_w
        bar_y2 = H - margin_b
        bar_y1 = bar_y2 - bar_h

        # Sample image brightness to pick contrasting color
        crop = (max(0, bar_x1), max(0, bar_y1 - tick_h - 10),
                min(W, bar_x2), min(H, bar_y2 + tick_h + 10))
        region = np.array(img.crop(crop))
        luminance = 0.299 * region[:, :, 0] + 0.587 * region[:, :, 1] + 0.114 * region[:, :, 2]
        text_color = "#FFFFFF" if luminance.mean() < 140 else "#000000"

        # Draw horizontal gradient bar (left = vmin, right = vmax)
        gradient = cmap(np.linspace(0, 1, bar_w))[:, :3]
        gradient_row = (gradient * 255).astype(np.uint8)[None, :, :]
        bar_arr = np.tile(gradient_row, (bar_h, 1, 1))
        bar_img = Image.fromarray(bar_arr, "RGB").convert("RGBA")
        img.paste(bar_img, (bar_x1, bar_y1))

        # Border
        draw = ImageDraw.Draw(img)
        draw.rectangle([bar_x1, bar_y1, bar_x2, bar_y2], outline=text_color, width=1)

        # Tick labels below the bar
        for tv, label_text in zip(tick_values, tick_labels_text):
            frac = (np.log10(tv) - log_vmin) / (log_vmax - log_vmin)
            frac = np.clip(frac, 0, 1)
            tick_x = bar_x1 + int(frac * bar_w)
            draw.line([(tick_x, bar_y2), (tick_x, bar_y2 + tick_gap)],
                      fill=text_color, width=1)
            tick_img = self._render_latex_label(label_text, font_size, text_color)
            tw, th = tick_img.size
            # Scale down if too wide
            max_tw = bar_w // max(len(tick_values) - 1, 1)
            if tw > max_tw:
                th = int(th * max_tw / tw)
                tw = max_tw
                tick_img = tick_img.resize((tw, th), Image.LANCZOS)
            paste_x = np.clip(tick_x - tw // 2, 0, W - tw)
            paste_y = bar_y2 + tick_gap + 1
            img.paste(tick_img, (paste_x, paste_y), tick_img)

        # Title label above the bar
        if label:
            title_fontsize = max(6, W // 90)
            title_latex = "$" + label + "$"
            label_img = self._render_latex_label(title_latex, title_fontsize, text_color)
            lw, lh = label_img.size
            if lw > bar_w:
                lh = int(lh * bar_w / lw)
                lw = bar_w
                label_img = label_img.resize((lw, lh), Image.LANCZOS)
            label_x = bar_x1 + bar_w // 2 - lw // 2
            label_y = bar_y1 - lh - 2
            img.paste(label_img, (label_x, label_y), label_img)

    def has_required_maps(self):
        return np.all([i in self.maps for i in self.required_maps])


class SinkVisSigmaGas(SinkVis):
    def __init__(self, params):
        self.required_maps = set(["sigma_gas"])
        super().__init__(params)
        if self.TaskDone:
            return
        self.AssignDefaultParams()

    def AssignDefaultParams(self):
        super().AssignDefaultParams()
        if self.params["filename"] is None:
            self.params["filename"] = (
                self.params["outputfolder"] + "/" + "SurfaceDensity_" + self.params["filename_suffix"]
            )

    def DetermineRequiredSnapdata(self):
        super().DetermineRequiredSnapdata()
        if self.render_maps:
            # gas data for surface density
            self.RequiredSnapdata += [
                "PartType0/Coordinates",
                "PartType0/Masses",
                "PartType0/ParticleIDs",
                "PartType0/KernelMaxRadius",
                "PartType0/ParticleChildIDsNumber",
                "PartType0/ParticleIDGenerationNumber",
            ]

    def GenerateMaps(self, snapdata):
        if "sigma_gas" not in self.maps.keys():
            from meshoid import GridSurfaceDensity
            self.maps["sigma_gas"] = GridSurfaceDensity(
                self.mass,
                self.pos,
                self.hsml,
                np.zeros(3),
                2 * self.params["rmax"],
                res=self.params["res"],
                parallel=self.parallel,
            ).T.clip(1e-100)
            # self.maps["sigma_gas"] = self.maps["sigma_gas"]
            np.savez_compressed(self.map_files["sigma_gas"], sigma_gas=self.maps["sigma_gas"])

    def MakeImages(self, snapdata):
        if self.params["limits"] is None:
            # Mass-weighted 1st/99th percentiles for surface density
            sigma_flat = np.sort(self.maps["sigma_gas"].ravel())
            cw = sigma_flat.cumsum() / sigma_flat.sum()
            self.params["limits"] = np.interp([0.01, 0.99], cw, sigma_flat)

        vmin, vmax = self.params["limits"]
        if vmax > vmin:
            f = (np.log10(self.maps["sigma_gas"]) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))
        else:
            f = np.zeros_like(self.maps["sigma_gas"])

        import matplotlib
        from matplotlib import pyplot as plt

        if self.params["backend"] == "PIL":
            from PIL import Image
            rgba = plt.get_cmap(self.params["cmap"])(np.flipud(f))
            self._pil_image = Image.fromarray((rgba * 255).astype(np.uint8), "RGBA")
            if not self.params["no_colorbar"]:
                self._add_colorbar_to_image(
                    vmin, vmax,
                    label=r"\Sigma_{\rm gas}\,\left(M_\odot\,\rm pc^{-2}\right)",
                )
        elif self.params["backend"] == "matplotlib":
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            matplotlib.use("Agg")
            self.fig, self.ax = plt.subplots(figsize=(4, 4))
            X = Y = np.linspace(-self.params["rmax"], self.params["rmax"], self.params["res"])
            X, Y = np.meshgrid(X, Y)
            p = self.ax.pcolormesh(
                X,
                Y,
                self.maps["sigma_gas"],
                norm=matplotlib.colors.LogNorm(vmin=self.params["limits"][0], vmax=self.params["limits"][1]),
                cmap=self.params["cmap"],
            )
            self.ax.set_aspect("equal")

            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.0)
            self.fig.colorbar(p, label=r"$\Sigma_{\rm gas}$ $(\rm M_\odot\,pc^{-2})$", cax=cax)
            if self.params["camera_distance"] == np.inf:
                self.ax.set_xlabel("X (pc)")
                self.ax.set_ylabel("Y (pc)")
            else:
                self.ax.set_xlabel("X (rad)")
                self.ax.set_ylabel("Y (rad)")

        super().MakeImages(snapdata)


class SinkVisCoolMap(SinkVis):
    def __init__(self, params):
        self.required_maps = set(
            ["sigma_gas", "sigma_1D"]
        )  # physical rendered quantities that can get saved and reused
        super().__init__(params)
        if self.TaskDone:
            return
        self.default_params["cool_cmap"] = "magma"
        self.default_params["v_limits"] = None
        self.AssignDefaultParams()

    def DetermineRequiredSnapdata(self):
        super().DetermineRequiredSnapdata()
        if self.render_maps:
            # gas data for surface density
            self.RequiredSnapdata += [
                "PartType0/Coordinates",
                "PartType0/Masses",
                "PartType0/ParticleIDs",
                "PartType0/KernelMaxRadius",
                "PartType0/ParticleChildIDsNumber",
                "PartType0/ParticleIDGenerationNumber",
            ]
            # extra velocity data for kinemtic map
            self.RequiredSnapdata += ["PartType0/Velocities"]

    def AssignDefaultParams(self):
        super().AssignDefaultParams()
        if self.params["filename"] is None:
            self.params["filename"] = self.params["outputfolder"] + "/" + "CoolMap_" + self.params["filename_suffix"]

    def GenerateMaps(self, snapdata):
        super().GenerateMaps(snapdata)
        from meshoid import GridSurfaceDensity

        if not "sigma_gas" in self.maps.keys():
            self.maps["sigma_gas"] = GridSurfaceDensity(
                self.mass,
                self.pos,
                self.hsml,
                np.zeros(3),
                2 * self.params["rmax"],
                res=self.params["res"],
                parallel=self.parallel,
            ).T
            np.savez_compressed(self.map_files["sigma_gas"], sigma_gas=self.maps["sigma_gas"])
        if not "sigma_1D" in self.maps.keys():
            # need to apply coordinate transforms to z-velocity
            v = np.copy(snapdata["PartType0/Velocities"])
            if hasattr(self, "_keep_mask"):
                v = v[self._keep_mask]
            self.DoCoordinateTransform(v, contravariant=True)
            sigma_1D = (
                GridSurfaceDensity(
                    self.mass * v[:, 2] ** 2,
                    self.pos,
                    self.hsml,
                    np.zeros(3),
                    2 * self.params["rmax"],
                    res=self.params["res"],
                    parallel=self.parallel,
                ).T
                / self.maps["sigma_gas"]
            )
            v_avg = (
                GridSurfaceDensity(
                    self.mass * v[:, 2],
                    self.pos,
                    self.hsml,
                    np.zeros(3),
                    2 * self.params["rmax"],
                    res=self.params["res"],
                    parallel=self.parallel,
                ).T
                / self.maps["sigma_gas"]
            )
            self.maps["sigma_1D"] = np.sqrt(sigma_1D - v_avg**2) / 1e3
            np.savez_compressed(self.map_files["sigma_1D"], sigma_1D=self.maps["sigma_1D"])

    def MakeImages(self, snapdata):
        if self.params["limits"] is None:
            # Mass-weighted 1st/99th percentiles for surface density
            sigma_flat = np.sort(self.maps["sigma_gas"].ravel())
            cw = sigma_flat.cumsum() / sigma_flat.sum()
            self.params["limits"] = np.interp([0.01, 0.99], cw, sigma_flat)
        if self.params["v_limits"] is None:
            # Energy-weighted percentiles: pixels with more kinetic energy count more
            order = self.maps["sigma_1D"].ravel().argsort()
            sigma_1D_sorted = self.maps["sigma_1D"].ravel()[order]
            Ekin = (self.maps["sigma_gas"] * self.maps["sigma_1D"] ** 2).ravel()[order]
            cw = Ekin.cumsum() / Ekin.sum()
            self.params["v_limits"] = np.interp([0.0, 0.99], cw, sigma_1D_sorted)
        fgas = (np.log10(self.maps["sigma_gas"]) - np.log10(self.params["limits"][0])) / np.log10(
            self.params["limits"][1] / self.params["limits"][0]
        )
        fgas = np.clip(fgas, 0, 1)
        from matplotlib import pyplot as plt
        from matplotlib.colors import rgb_to_hsv, hsv_to_rgb
        mapcolor = plt.get_cmap(self.params["cool_cmap"])(
            np.log10(self.maps["sigma_1D"] / self.params["v_limits"][0])
            / np.log10(self.params["v_limits"][1] / self.params["v_limits"][0])
        )
        # blend HSV: use fgas as intensity to modulate saturation and value
        hsv = rgb_to_hsv(mapcolor[:, :, :3])
        intensity = 2 * fgas - 1  # remap [0,1] -> [-1,1]
        hue, sat, val = np.moveaxis(hsv, -1, 0)
        bright = intensity > 0
        dark = intensity < 0
        nontrivial_sat = np.abs(sat) > 1e-10
        np.putmask(sat, nontrivial_sat & bright, (1 - intensity) * sat)
        np.putmask(sat, nontrivial_sat & dark, (1 + intensity) * sat - intensity)
        np.putmask(val, bright, (1 - intensity) * val + intensity)
        np.putmask(val, dark, (1 + intensity) * val)
        np.clip(hsv[:, :, 1:], 0, 1, out=hsv[:, :, 1:])
        cool_data = hsv_to_rgb(hsv)
        self.maps["coolmap"] = cool_data

        if self.params["backend"] == "PIL":
            from PIL import Image
            rgb = np.flipud(self.maps["coolmap"])
            self._pil_image = Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8), "RGB").convert("RGBA")
        else:
            plt.imsave(
                self.params["filename_incomplete"], np.flipud(self.maps["coolmap"])
            )
        super().MakeImages(snapdata)


class SinkVisNarrowbandComposite(SinkVis):
    def __init__(self, params):
        self.required_maps = set(["SHO_RGB"])  # RGB map of SII, Halpha, and OIII
        super().__init__(params)
        if self.TaskDone:
            return
        self.AssignDefaultParams()

    def DetermineRequiredSnapdata(self):
        super().DetermineRequiredSnapdata()
        if self.render_maps:
            self.RequiredSnapdata += [
                "PartType0/Coordinates",
                "PartType0/Masses",
                "PartType0/ParticleIDs",
                "PartType0/Temperature",
                "PartType0/ElectronAbundance",
                "PartType0/KernelMaxRadius",
                "PartType0/ParticleChildIDsNumber",
                "PartType0/ParticleIDGenerationNumber",
                "PartType0/Density",
                "PartType0/HII",
            ]

    def AssignDefaultParams(self):
        super().AssignDefaultParams()
        if self.params["filename"] is None:
            self.params["filename"] = (
                self.params["outputfolder"] + "/" + "NarrowbandComposite_" + self.params["filename_suffix"]
            )

    #            self.params["filename_incomplete"] = self.params["filename"].replace(".png",".incomplete.png")

    def GenerateMaps(self, snapdata):
        super().GenerateMaps(snapdata)

        if not "SHO_RGB" in self.maps.keys():
            from skimage.color import rgb2hsv, hsv2rgb
            from meshoid.radiation import radtransfer
            # print("Generating SHO map...")
            _k = self._keep_mask if hasattr(self, "_keep_mask") else slice(None)
            rho = snapdata["PartType0/Density"][_k]
            T = snapdata["PartType0/Temperature"][_k]
            fe = snapdata["PartType0/ElectronAbundance"][_k]
            hii = snapdata["PartType0/HII"][_k]
            UnitDensity = snapdata["Header"].get("UnitMass_In_CGS", 1.989e33) / snapdata["Header"].get("UnitLength_In_CGS", 3.086e18)**3
            nH = rho * UnitDensity / _ac.m_p.cgs.value
            ne = nH * fe

            # ne = np.clip(ne,None,np.percentile(ne,100*(1.0-100/len(ne))) ) #clip by 100th largest value in case we have few rogue cells with extremely large values
            #            ne = np.clip(ne,None,np.percentile(ne,99)) #clip by 99th percentile, this removes some too bright pixels in Ha, usually at interfaces with dense regions, which are hard to interpolate anyway

            wavelength = 6562
            T4 = T / 1e4
            j_B_Ha = 1.24e-25 * (T4) ** (-0.942 - 0.031 * np.log(T4)) * 2.86 * nH * hii * ne

            wavelength = 5007
            ncrit = 1e3
            j_OIII = (
                nH
                * hii
                * ne
                * np.exp(-1.5 / T4)
                / (1 + (ne / ncrit) * T4**-0.5)
                * np.exp(-14388 / wavelength / T4)
                * T4**-0.5
                * np.exp(-14388 / wavelength / T4)
            )

            wavelength = 6584
            ncrit = 1e3
            j_NII = nH * hii * ne / (1 + (ne / ncrit) * T4**-0.5) * np.exp(-14388 / wavelength / T4) * T4**-0.5

            wavelength = 6716
            ncrit = 1e3
            j_SII = (
                nH
                * hii
                * ne
                * np.exp(-1 / T4)
                / (1 + (ne / ncrit) * T4**-0.5)
                * np.exp(-14388 / wavelength / T4)
                * T4**-0.5
            )

            pc_to_cm = _au.pc.to(_au.cm)
            msun_to_g = _ac.M_sun.cgs.value

            lum = np.c_[j_B_Ha, j_OIII, j_SII] * pc_to_cm**3 * (self.mass / rho)[:, None]
            #            lum = np.c_[j_NII,j_OIII,j_SII] * pc_to_cm**3 *  (snapdata["PartType0/Masses"]/rho)[:,None] #NII behaves much better for interpolation than Ha, mostly because it does not diverge in the limit of ne->inf, similar to OIII and SII

            # Here we normalize the emission, but it is done relative to the current emissions, so each snapshot has a different absolute normalization. If an absolute normalization is desired across snapshots, this part should be commented out and a vector should be used for SHO_RGB_norm to normalize individual channels
            lum_sum = np.sum(lum, axis=0)
            full_sum = np.sum(lum_sum)
            if full_sum:
                if lum_sum[0]:
                    lum[:, 0] *= full_sum / lum_sum[0]
                if lum_sum[1]:
                    lum[:, 1] *= full_sum / lum_sum[1]
                if lum_sum[2]:
                    lum[:, 2] *= full_sum / lum_sum[2]

            def get_color_matrix(rot):
                a = np.eye(3)
                a = rgb2hsv(a)
                a[:, 0] += rot
                a[:, 0] = a[:, 0] % 1
                a = hsv2rgb(a)
                return a

            color_matrix = get_color_matrix(-0.3)
            lum = lum @ color_matrix
            kappa = np.array([159, 233, 283]) * np.ones_like(self.mass)[:, None] / (pc_to_cm**2 / msun_to_g)
            if self.params["camera_distance"] < np.inf:
                lum[:] /= (
                    self.r[:, None] ** 2
                )  # have to convert here because smoothing lengths are now in angular units
                lum[self.pos[:, 2] < 0] = 0  # ignore stuff behind the camera
                kappa[self.pos[:, 2] < 0] = 0

            self.maps["SHO_RGB"] = radtransfer(
                np.copy(lum),
                np.copy(self.mass),
                np.copy(kappa),
                np.copy(self.pos),
                np.copy(self.hsml),
                self.params["res"],
                2 * self.params["rmax"],
            ).swapaxes(0, 1)

            np.savez_compressed(self.map_files["SHO_RGB"], SHO_RGB=self.maps["SHO_RGB"])

    def MakeImages(self, snapdata):
        # sigmoid = lambda x: x/np.sqrt(1+x*x) # tapering function to soften the saturation
        sigmoid = lambda x: x / (1 + x)  # tapering function to soften the saturation
        ha_map = np.copy(self.maps["SHO_RGB"])

        if hasattr(self.params["SHO_RGB_norm"], "__iter__"):  # normalization constant per channel provided
            norm = self.params["SHO_RGB_norm"]
        elif self.params["SHO_RGB_norm"] == 0:  # guess normalization for each channel
            norm = [
                np.percentile(ha_map[:, :, 2], 99),
                np.percentile(ha_map[:, :, 0], 99),
                np.percentile(ha_map[:, :, 1], 99),
            ]
            print("Using SHO_RGB normalizations %g %g %g" % (norm[0], norm[1], norm[2]))
        else:  # use the same normlization constant for each channel
            norm = [self.params["SHO_RGB_norm"], self.params["SHO_RGB_norm"], self.params["SHO_RGB_norm"]]
        self.maps["SHO_RGB"][:, :, 0] = sigmoid(ha_map[:, :, 2] / norm[0])
        self.maps["SHO_RGB"][:, :, 1] = sigmoid(ha_map[:, :, 0] / norm[1])
        self.maps["SHO_RGB"][:, :, 2] = sigmoid(ha_map[:, :, 1] / norm[2])

        if self.params["backend"] == "PIL":
            from PIL import Image
            rgb = self.maps["SHO_RGB"][::-1]
            self._pil_image = Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8), "RGB").convert("RGBA")
        else:
            from matplotlib import pyplot as plt
            plt.imsave(
                self.params["filename_incomplete"], self.maps["SHO_RGB"][::-1]
            )
        super().MakeImages(snapdata)


# ---------- Custom field task infrastructure ----------

# Physical constants and unit conversions (CGS, from astropy)
_CONSTANTS = {
    "pi": np.pi,
    "k_B": _ac.k_B.cgs.value,
    "m_p": _ac.m_p.cgs.value,
    "m_e": _ac.m_e.cgs.value,
    "c_light": _ac.c.cgs.value,
    "G": _ac.G.cgs.value,
    "Msun": _ac.M_sun.cgs.value,
    "Lsun": _ac.L_sun.cgs.value,
    "pc": _au.pc.to(_au.cm),
    "AU": _au.AU.to(_au.cm),
    "yr": _au.yr.to(_au.s),
    "Myr": _au.Myr.to(_au.s),
    "eV": _au.eV.to(_au.erg),
}

# Registry of derived fields: name -> expression string.
# Expressions can reference snapshot fields, constants, other derived fields,
# and the builtins (abs, sqrt, norm, log, log10, exp, ...).
DERIVED_FIELDS = {}

# Fallback expressions for snapshot fields that may not exist.
# Used only when PartType0/<name> is absent from the snapshot.
FIELD_FALLBACKS = {}


def register_derived_field(name, expr):
    """Register a named derived field.

    >>> register_derived_field("MagneticPressure", "norm(MagneticField)**2 / (8*pi)")
    """
    DERIVED_FIELDS[name] = expr


def register_field_fallback(name, expr):
    """Register a fallback expression for a snapshot field.

    The fallback is used only when PartType0/<name> is not present in the
    snapshot.  If the field exists on disk, the snapshot value is used.

    >>> register_field_fallback("Pressure", "(5./3 - 1) * Density * InternalEnergy")
    """
    FIELD_FALLBACKS[name] = expr


# Built-in derived fields (Gaussian CGS conventions matching GIZMO)
register_derived_field("MagneticPressure", "norm(MagneticField)**2 / (8*pi)")
register_derived_field("PlasmaBeta", "Pressure / MagneticPressure")
register_derived_field("AlfvenSpeed", "norm(MagneticField) / sqrt(4*pi*Density)")
register_derived_field("MachNumber", "norm(Velocities) / SoundSpeed")
register_derived_field("JeansLength", "SoundSpeed / sqrt(G * Density)")
register_derived_field("ThermalEnergy", "Masses * InternalEnergy")
register_derived_field("KineticEnergy", "0.5 * Masses * norm(Velocities)**2")
register_derived_field("MagneticEnergy", "norm(MagneticField)**2 / (8*pi) * Masses / Density")
register_derived_field("NumberDensity", "Density / m_p")
register_derived_field("Entropy", "Pressure / Density**(5./3)")

# Photon energy density fields (eV/cm^3) — col(PhotonEnergy, band) extracts
# band 0=EUV, 1=FUV, 2=NUV, 3=ONIR, 4=FIR from the (N,5) PhotonEnergy array.
# PhotonEnergy is specific (per unit mass) in code units, so:
#   energy_density [code] = PhotonEnergy * Density / Masses
#   energy_density [CGS]  = energy_density [code] * UnitEnergyDensity_In_CGS
#   energy_density [eV/cm^3] = energy_density [CGS] / eV
register_derived_field("PhotonEnergyDensity_EUV", "col(PhotonEnergy, 0) * Density / Masses * UnitEnergyDensity_In_CGS / eV")
register_derived_field("PhotonEnergyDensity_FUV", "col(PhotonEnergy, 1) * Density / Masses * UnitEnergyDensity_In_CGS / eV")
register_derived_field("PhotonEnergyDensity_NUV", "col(PhotonEnergy, 2) * Density / Masses * UnitEnergyDensity_In_CGS / eV")
register_derived_field("PhotonEnergyDensity_ONIR", "col(PhotonEnergy, 3) * Density / Masses * UnitEnergyDensity_In_CGS / eV")
register_derived_field("PhotonEnergyDensity_FIR", "col(PhotonEnergy, 4) * Density / Masses * UnitEnergyDensity_In_CGS / eV")

# G0: FUV photon energy density in Habing units (1 Habing = 5.29e-14 erg/cm^3)
register_derived_field("G0", "col(PhotonEnergy, 1) * Density / Masses * UnitEnergyDensity_In_CGS / 5.29e-14")

# LaTeX symbols for colorbar labels.  Keys can be snapshot field names,
# derived field names, or full expression strings.
FIELD_SYMBOLS = {}


def register_field_symbol(name, latex):
    r"""Register a LaTeX symbol for a field or expression.

    >>> register_field_symbol("Temperature", r"T\;\mathrm{(K)}")
    """
    FIELD_SYMBOLS[name] = latex


# Built-in symbols
register_field_symbol("Density", r"\rho")
register_field_symbol("Temperature", r"T\;\mathrm{(K)}")
register_field_symbol("Pressure", r"P")
register_field_symbol("InternalEnergy", r"u")
register_field_symbol("Masses", r"M")
register_field_symbol("SoundSpeed", r"c_s")
register_field_symbol("MagneticPressure", r"P_B")
register_field_symbol("PlasmaBeta", r"\beta")
register_field_symbol("AlfvenSpeed", r"v_A")
register_field_symbol("MachNumber", r"\mathcal{M}")
register_field_symbol("JeansLength", r"\lambda_J")
register_field_symbol("NumberDensity", r"n\;\mathrm{(cm^{-3})}")
register_field_symbol("ThermalEnergy", r"E_\mathrm{th}")
register_field_symbol("KineticEnergy", r"E_\mathrm{kin}")
register_field_symbol("MagneticEnergy", r"E_B")
register_field_symbol("Entropy", r"P/\rho^{\gamma}")
register_field_symbol("PhotonEnergyDensity_EUV", r"u_\mathrm{EUV}\;\mathrm{(eV\,cm^{-3})}")
register_field_symbol("PhotonEnergyDensity_FUV", r"u_\mathrm{FUV}\;\mathrm{(eV\,cm^{-3})}")
register_field_symbol("PhotonEnergyDensity_NUV", r"u_\mathrm{NUV}\;\mathrm{(eV\,cm^{-3})}")
register_field_symbol("PhotonEnergyDensity_ONIR", r"u_\mathrm{ONIR}\;\mathrm{(eV\,cm^{-3})}")
register_field_symbol("PhotonEnergyDensity_FIR", r"u_\mathrm{FIR}\;\mathrm{(eV\,cm^{-3})}")
register_field_symbol("G0", r"G_0")


# Built-in fallbacks for fields that GIZMO may or may not write
register_field_fallback("Pressure", "(5./3 - 1) * Density * InternalEnergy")
register_field_fallback("SoundSpeed", "sqrt(5./3 * (5./3 - 1) * InternalEnergy)")
register_field_fallback("Temperature", "(5./3 - 1) * InternalEnergy * m_p / k_B")


# Regex to extract identifier tokens, skipping the 'e'/'E' in scientific notation
_TOKEN_RE = re.compile(r"(?<!\d)[A-Za-z_]\w*")

# Tokens that are builtins, NOT field names
_UNIT_NAMES = {
    "UnitLength_In_CGS", "UnitMass_In_CGS", "UnitVelocity_In_CGS",
    "UnitEnergyDensity_In_CGS", "UnitDensity_In_CGS",
    "UnitTime_In_CGS", "UnitEnergy_In_CGS",
}

_EXPR_BUILTINS = {
    "abs", "sqrt", "norm", "col", "log", "log2", "log10", "exp",
    "sin", "cos", "tan", "minimum", "maximum", "clip", "where",
} | set(_CONSTANTS.keys()) | _UNIT_NAMES


def _extract_field_names(expr):
    """Extract field name tokens from an expression, resolving derived fields
    and fallbacks recursively to their base PartType0 snapshot fields.

    Returns all snapshot fields that *might* be needed — both the primary
    field and any fields its fallback expression requires.
    """
    tokens = _TOKEN_RE.findall( expr)
    raw = set(t for t in tokens if t not in _EXPR_BUILTINS)

    base_fields = set()
    seen = set()

    def _resolve(names):
        for name in names:
            if name in seen:
                continue
            seen.add(name)
            if name in DERIVED_FIELDS:
                sub = _TOKEN_RE.findall( DERIVED_FIELDS[name])
                _resolve(t for t in sub if t not in _EXPR_BUILTINS)
            else:
                # This is a snapshot field (or has a fallback).  Request
                # both the field itself and any fields the fallback needs,
                # since we won't know until runtime which is available.
                base_fields.add(name)
                if name in FIELD_FALLBACKS:
                    sub = _TOKEN_RE.findall( FIELD_FALLBACKS[name])
                    _resolve(t for t in sub if t not in _EXPR_BUILTINS)

    _resolve(raw)
    return sorted(base_fields)


def _eval_field_expr(expr, snapdata, _cache=None):
    """Evaluate a field expression against loaded PartType0 snapshot data.

    Field names (e.g. 'Masses', 'Temperature') are resolved to their
    PartType0 arrays.  Derived fields are evaluated recursively and cached
    within the call tree.  Numpy ufuncs and physical constants are available.
    """
    if _cache is None:
        _cache = {}

    def _norm(x):
        return np.sqrt(np.sum(np.asarray(x) ** 2, axis=-1))

    def _col(arr, i):
        return np.asarray(arr)[:, int(i)]

    ns = {
        "np": np, "abs": np.abs, "sqrt": np.sqrt, "norm": _norm, "col": _col,
        "log": np.log, "log2": np.log2, "log10": np.log10,
        "exp": np.exp, "sin": np.sin, "cos": np.cos, "tan": np.tan,
        "minimum": np.minimum, "maximum": np.maximum,
        "clip": np.clip, "where": np.where,
    }
    ns.update(_CONSTANTS)

    # Add code-unit conversion factors from snapshot header
    header = snapdata.get("Header", {})
    UL = header.get("UnitLength_In_CGS", 1.0)
    UM = header.get("UnitMass_In_CGS", 1.0)
    UV = header.get("UnitVelocity_In_CGS", 1.0)
    ns["UnitLength_In_CGS"] = UL
    ns["UnitMass_In_CGS"] = UM
    ns["UnitVelocity_In_CGS"] = UV
    ns["UnitEnergyDensity_In_CGS"] = UM * UV**2 / UL**3
    ns["UnitDensity_In_CGS"] = UM / UL**3
    ns["UnitTime_In_CGS"] = UL / UV
    ns["UnitEnergy_In_CGS"] = UM * UV**2

    tokens = _TOKEN_RE.findall(expr)
    for name in set(tokens) - _EXPR_BUILTINS:
        if name in ns:
            continue
        if name in _cache:
            ns[name] = _cache[name]
            continue

        # 1) Explicit derived field — always computed from expression
        if name in DERIVED_FIELDS:
            _cache[name] = _eval_field_expr(DERIVED_FIELDS[name], snapdata, _cache)
            ns[name] = _cache[name]
            continue

        # 2) Snapshot field — use if present
        key = "PartType0/" + name
        if key in snapdata and snapdata[key] is not None:
            ns[name] = snapdata[key]
            continue

        # 3) Fallback expression — used when snapshot field is missing
        if name in FIELD_FALLBACKS:
            _cache[name] = _eval_field_expr(FIELD_FALLBACKS[name], snapdata, _cache)
            ns[name] = _cache[name]
            continue

        raise KeyError(
            f"'{name}' is not a snapshot field (PartType0/{name}), "
            f"a derived field, a fallback, or a builtin"
        )
    return eval(expr, {"__builtins__": {}}, ns)


def parse_custom_task(spec):
    """Parse a task spec like 'SurfaceDensity(Masses*Temperature)'.

    Returns (render_mode, field_expr) or None if not a custom task.
    """
    m = re.match(r"^(SurfaceDensity|Projection|ProjectedAverage|Slice)\((.+)\)$", spec)
    if m:
        return m.group(1), m.group(2)
    return None


class SinkVisCustomField(SinkVis):
    """Generic task that renders an arbitrary field expression.

    Supports four render modes:
      - SurfaceDensity(expr): integral of expr along the line of sight  (Σ)
      - Projection(expr):     same as SurfaceDensity (alias)
      - ProjectedAverage(expr): mass-weighted projected average of expr
      - Slice(expr):          midplane slice of expr
    """

    def __init__(self, params):
        self._render_mode = params["_render_mode"]
        self._field_expr = params["_field_expr"]
        self._map_key = f"{self._render_mode}_{self._field_expr}"
        self.required_maps = set([self._map_key])
        super().__init__(params)
        if self.TaskDone:
            return
        self.AssignDefaultParams()

    def AssignDefaultParams(self):
        super().AssignDefaultParams()
        if self.params["filename"] is None:
            safe_expr = re.sub(r"[^\w]", "_", self._field_expr)
            self.params["filename"] = (
                self.params["outputfolder"]
                + "/"
                + f"{self._render_mode}_{safe_expr}_"
                + self.params["filename_suffix"]
            )

    def DetermineRequiredSnapdata(self):
        super().DetermineRequiredSnapdata()
        # Always need coordinates, masses, smoothing lengths for rendering
        self.RequiredSnapdata += [
            "PartType0/Coordinates",
            "PartType0/Masses",
            "PartType0/KernelMaxRadius",
            "PartType0/ParticleIDs",
            "PartType0/ParticleChildIDsNumber",
            "PartType0/ParticleIDGenerationNumber",
        ]
        # Add any fields referenced in the expression
        for name in _extract_field_names(self._field_expr):
            self.RequiredSnapdata.append("PartType0/" + name)

    def _colorbar_label(self):
        """Return a LaTeX string for the colorbar title."""
        expr = self._field_expr
        # Check for exact match on the expression or field name
        if expr in FIELD_SYMBOLS:
            return FIELD_SYMBOLS[expr]
        # For simple single-field expressions, look up the field
        tokens = [t for t in _TOKEN_RE.findall( expr) if t not in _EXPR_BUILTINS]
        if len(tokens) == 1 and tokens[0] == expr and expr in FIELD_SYMBOLS:
            return FIELD_SYMBOLS[expr]
        # Fall back to rendering the raw expression in mathtt
        return r"\mathtt{" + expr.replace("_", r"\_") + "}"

    def GenerateMaps(self, snapdata):
        if self._map_key in self.maps:
            return
        from meshoid import Meshoid

        f = _eval_field_expr(self._field_expr, snapdata)

        # Handle vector fields: if f has multiple columns, take magnitude
        if f.ndim > 1:
            f = np.sqrt(np.sum(f ** 2, axis=1))

        # Apply the same cull mask that SetupCoordsAndWeights applied to pos/mass/hsml
        if hasattr(self, "_keep_mask"):
            f = f[self._keep_mask]

        res = self.params["res"]
        rmax = self.params["rmax"]
        # Build a Meshoid from the already-transformed coordinates
        n_jobs = self.params["threads"] if self.params["threads"] > 0 else -1
        M = Meshoid(self.pos, self.mass, self.hsml, n_jobs=n_jobs)

        if self._render_mode == "SurfaceDensity":
            # Surface density: f is an extensive/conserved quantity (e.g. Masses)
            result = M.SurfaceDensity(
                f, center=np.zeros(3), size=2 * rmax, res=res,
            ).T
        elif self._render_mode == "Projection":
            # Projection: f is a volume density / intensive quantity (e.g. Density)
            # computes the line integral ∫ f dz
            result = M.Projection(
                f, center=np.zeros(3), size=2 * rmax, res=res,
            ).T
        elif self._render_mode == "ProjectedAverage":
            result = M.ProjectedAverage(
                f, center=np.zeros(3), size=2 * rmax, res=res,
            ).T
        elif self._render_mode == "Slice":
            # Supersample and downsample to anti-alias Voronoi edges
            ss = int(self.params.get("supersample", 2))
            # For positive quantities, slice in log space to guarantee positivity
            positive = np.all(f > 0)
            slice_f = np.log(f) if positive else f
            hi_res = M.Slice(
                slice_f, center=np.zeros(3), size=2 * rmax, res=res * ss,
                order=1, slope_limiter=True,
            ).T
            if positive:
                hi_res = np.exp(hi_res)
            if ss > 1:
                if positive:
                    result = np.exp(
                        np.log(hi_res).reshape(res, ss, res, ss).mean(axis=(1, 3))
                    )
                else:
                    result = hi_res.reshape(res, ss, res, ss).mean(axis=(1, 3))
                self._slice_hires = hi_res
            else:
                result = hi_res
        else:
            raise ValueError(f"Unknown render mode: {self._render_mode}")

        self.maps[self._map_key] = result
        np.savez_compressed(self.map_files[self._map_key], **{self._map_key: result})

    def MakeImages(self, snapdata):
        data = self.maps[self._map_key]
        positive = data > 0

        if self.params["limits"] is None:
            # Use hi-res slice data for limits if available (before AA smoothing)
            limit_data = getattr(self, "_slice_hires", data)
            if self._render_mode in ("SurfaceDensity", "Projection"):
                # Mass-weighted 1st/99th percentiles for integral quantities
                flat = np.sort(limit_data.ravel())
                cw = flat.cumsum() / flat.sum()
                self.params["limits"] = np.interp([0.01, 0.99], cw, flat)
            else:
                # Raw 1st/99th percentiles for slice/projected average
                self.params["limits"] = np.percentile(limit_data, [1, 99])

        vmin, vmax = self.params["limits"]
        if vmax <= vmin:
            f = np.zeros_like(data)
        elif positive.all():
            f = (np.log10(data) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))
        else:
            f = (data - vmin) / (vmax - vmin)
        f = np.clip(f, 0, 1)

        from matplotlib import pyplot as plt

        if self.params["backend"] == "PIL":
            from PIL import Image
            rgba = plt.get_cmap(self.params["cmap"])(np.flipud(f))
            self._pil_image = Image.fromarray((rgba * 255).astype(np.uint8), "RGBA")
            if not self.params["no_colorbar"]:
                self._add_colorbar_to_image(vmin, vmax, label=self._colorbar_label())
        elif self.params["backend"] == "matplotlib":
            import matplotlib
            matplotlib.use("Agg")
            self.fig, self.ax = plt.subplots(figsize=(4, 4))
            X = Y = np.linspace(-self.params["rmax"], self.params["rmax"], self.params["res"])
            X, Y = np.meshgrid(X, Y)
            if positive.all():
                norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
            else:
                norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            p = self.ax.pcolormesh(X, Y, data, norm=norm, cmap=self.params["cmap"])
            self.ax.set_aspect("equal")
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.0)
            cb_label = "$" + self._colorbar_label() + "$"
            self.fig.colorbar(p, label=cb_label, cax=cax)
            if self.params["camera_distance"] == np.inf:
                self.ax.set_xlabel("X (pc)")
                self.ax.set_ylabel("Y (pc)")
            else:
                self.ax.set_xlabel("X (rad)")
                self.ax.set_ylabel("Y (rad)")

        super().MakeImages(snapdata)
