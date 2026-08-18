import numpy as np
import re
import inspect as _inspect
from astropy import constants as _ac, units as _au
from .misc_functions import *
from os.path import isfile
import json
import hashlib
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


# smallest sink marker size, in the units of the log-mass scale used by AddStarsToImage
_MIN_STAR_SIZE = 1.0


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
    # maps whose values set a color-limit param, as {param name: map name}.  Used to
    # compute one color scale for a whole sequence of frames; tasks that don't take
    # color limits (e.g. RGB composites) leave this empty.
    color_limit_maps = {}

    # code length -> matplotlib axis unit; set by _set_axis_length_unit
    _axis_unit_factor = 1.0

    @property
    def splat_kwargs(self):
        """meshoid deposition backend/precision for this run; see splat_kwargs()."""
        return splat_kwargs(self.params.get("gpu", "off"))

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
            "backend": "matplotlib",
            "rescale_hsml": False,
            "recompute_hsml": False,
            "FOV": 90,
            "camera_distance": np.inf,
            "center_on_star": False,
            "center_on_ID": False,
            "center_on_densest": False,
            "realstars": True,
            "realstars_opacity": 1.0,
            "time_offset": 0,
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
            "star_legend": False,
            "draw_axes": False,
            "no_map_cache": False,
            "gpu": "off",
            "tight_bbox": True,
            "overwrite": True,
            "unit_scalefac": 1,
            "outputfolder": ".",
            "SHO_RGB_norm": 0,
            "outflow_only": False,
            "no_colorbar": False,
            "order": 0,
            # extra label inserted into the output filename, so that several
            # renders of the same snapshot (e.g. one per sink) can coexist in
            # one output folder instead of overwriting each other
            "filename_tag": "",
        }

        self.AssignDefaultParams()
        if self.params["camera_distance"] < np.inf and self.params["backend"] == "matplotlib":
            self.params["backend"] = "PIL"
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
            "recompute_hsml",
            "res",
            "realstars_max_lum",
            "realstars_lum_exp",
            "realstars_opacity",
            "order",
        ]
        dump = {}
        for k in self.params_that_affect_maps:
            if type(self.params[k]) == np.ndarray:
                dump[k] = self.params[k].tolist()
            else:
                dump[k] = self.params[k]
        # Use md5 rather than the built-in hash(): Python's hash() of strings is
        # PYTHONHASHSEED-salted per interpreter, so cache filenames would differ
        # between runs (callers used to work around this with an os.execv re-exec
        # setting PYTHONHASHSEED=0). md5 is deterministic across runs and Python
        # versions, removing that requirement.
        self.params_hash = hashlib.md5(
            json.dumps(dump, sort_keys=True).encode()
        ).hexdigest()

        mapdir = self.params["outputfolder"] + "/.maps"
        # makedirs creates the output folder itself if needed, and exist_ok
        # tolerates another worker process winning the race to create it.
        # (os.mkdir in a retry loop spun forever when outputfolder didn't exist.)
        os.makedirs(mapdir, exist_ok=True)

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
            self.RequiredSnapdata += [
                "PartType0/Density", "PartType0/Coordinates",
                "PartType0/Masses", "PartType0/KernelMaxRadius",
            ]
        if self.params["outflow_only"]:
            self.RequiredSnapdata += ["PartType0/Velocities", "PartType5/Velocities"]
        # check if we have maps already saved
        self.render_maps = False
        if self.params["no_map_cache"]:
            self.render_maps = True
            return
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
        # the tag goes in front of the frame index so that everything belonging
        # to one tag sorts together and shares a movie prefix
        tag = self.params.get("filename_tag") or ""
        self.params["filename_suffix"] = "%s%s_%s_%s.png" % (
            tag + "_" if tag else "",
            str(self.params["index"]).zfill(5),
            str(round(self.params["pan"] * 10)).zfill(4),
            self.params["cubemap_dir"],
        )

    def DoCoordinateTransform(self, x, m=None, h=None, contravariant=False, update_r=True, pos=None):
        # center on the designated center coordinate
        if not contravariant:
            x[:] -= self.params["center"]

        # orient along the camera basis if one is specified (--dir / camera_dir);
        # same matrix for coordinates and vectors: components along the camera
        # axes are (v.right, v.up, v.dir) for any vector
        if self.params["camera_dir"] is not None:
            x[:] = (self.camera_matrix @ x.T).T
        # pan/tilt: about the world axes without a camera basis, else composed
        # in the camera frame so a pan orbits the chosen view (a freeze rotation
        # with --dir used to silently render the same angle every frame).
        # Vectors rotate with the same matrix as coordinate displacements:
        # the transformed v is d/dt of the transformed x.
        tilt, pan = self.params["tilt"], self.params["pan"]
        # first pan
        cosphi, sinphi = np.cos(np.pi * pan / 180), np.sin(np.pi * pan / 180)
        x[:] = np.stack([cosphi * x[:, 0] + sinphi * x[:, 2], x[:, 1], -sinphi * x[:, 0] + cosphi * x[:, 2]], 1)
        # then tilt
        costheta, sintheta = np.cos(np.pi * tilt / 180), np.sin(np.pi * tilt / 180)
        x[:] = np.stack(
            [x[:, 0], costheta * x[:, 1] + sintheta * x[:, 2], -sintheta * x[:, 1] + costheta * x[:, 2]], 1
        )

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
                # pos overrides self.pos when the vector field is evaluated on a
                # different particle set (e.g. the unculled arrays Div/Curl need)
                global_coords = np.copy(self.pos if pos is None else pos)
                # multiply by z, now we're in the rotated real space frame
                global_coords[:, :2] *= -global_coords[:, 2][:, None]
                # get the radial component
                x[:, 2] = np.sum(x * global_coords, axis=1) / np.sum(global_coords**2, axis=1) ** 0.5

    def SetupCameraBasis(self):
        """Build the camera coordinate basis from camera_dir/camera_up.

        Depends only on params, not on any snapshot data, and is idempotent.  It has
        to be available even when the maps come from cache and SetupCoordsAndWeights
        never runs, because anything drawn on top of the map (sinks, vectors) still
        needs to be transformed into the camera frame.
        """
        if self.params["camera_dir"] is None:
            return  # no camera basis: DoCoordinateTransform uses the pan/tilt path instead
        self.camera_dir = self.params["camera_dir"]
        NormalizeVector(self.params["camera_dir"])
        if self.params["camera_up"] is None:
            # default "up" direction is +y, we will project it out if the camera is tilted
            self.camera_up = np.array([0, 1.0, 0])
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

    def SetupCoordsAndWeights(self, snapdata):
        from meshoid import Meshoid

        res = self.params["res"]
        if "PartType0/Coordinates" in snapdata.keys():
            have_radii = snapdata.get("PartType0/KernelMaxRadius") is not None
            recompute = self.params["recompute_hsml"] and not snapdata.get("_hsml_recomputed")
            if not have_radii or recompute:
                # a periodic tree is only valid if the coordinates are in the box frame
                coords = snapdata["PartType0/Coordinates"]
                boxsize = snapdata["Header"]["BoxSize"]
                snapdata["PartType0/KernelMaxRadius"] = Meshoid(
                    coords, boxsize=boxsize if coords_in_box_frame(coords, boxsize) else None
                ).SmoothingLength()
                # sticky for this loaded snapshot so later tasks reuse the tree's result
                snapdata["_hsml_recomputed"] = True
            self.pos, self.mass, self.hsml = (
                np.copy(snapdata["PartType0/Coordinates"]),
                np.copy(snapdata["PartType0/Masses"]),
                np.copy(snapdata["PartType0/KernelMaxRadius"]),
            )  # copy these because we don't want to modify them
            if self.params["rescale_hsml"]:
                self.hsml *= self.params["rescale_hsml"]
        else:
            # sink-only snapshot (e.g. pure N-body): keep the gas arrays defined but empty so
            # the splatting routines return blank maps instead of the task dying on a missing
            # attribute.  The stars still render.
            self.pos = np.zeros((0, 3))
            self.mass = np.zeros(0)
            self.hsml = np.zeros(0)
            self._keep_mask = np.zeros(0, dtype=bool)

        if self.params["outflow_only"]:
            if "PartType5/Coordinates" in snapdata.keys() and len(self.mass):
                from scipy.spatial import KDTree
                # find nearest star to each gas cell
                _, ngb = KDTree(snapdata["PartType5/Coordinates"]).query(self.pos)
                dx = self.pos - snapdata["PartType5/Coordinates"][ngb]
                dv = snapdata["PartType0/Velocities"] - snapdata["PartType5/Velocities"][ngb]
                self.mass *= (dx * dv).sum(1) > 0

        self.SetupCameraBasis()

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

    def GenerateBlankMaps(self):
        """Stand-in gas maps for a snapshot with no gas.

        Every renderer reads gas fields straight out of the snapshot, so rather than teach each
        one to guard, hand them empty maps and let them draw the stars over a blank field.  Not
        cached -- there is nothing here worth reusing.
        """
        for name in self.required_maps:
            if name != "realstars" and name not in self.maps:  # realstars renders from the sinks
                self.maps[name] = np.zeros(self.blank_map_shape(name))

    def blank_map_shape(self, name):
        return self.params["res"], self.params["res"]

    @staticmethod
    def _mass_weighted_limits(sigma, quantiles=(0.01, 0.99)):
        """Mass-weighted percentile color limits for a surface density map.

        A frame with no gas in it makes the weighted CDF 0/0, and any limits are arbitrary then,
        so fall back to a positive non-degenerate pair the log color scales can accept.  Nothing
        is read off it: a uniformly zero map clips to the bottom of the scale everywhere.
        """
        flat = np.sort(np.asarray(sigma).ravel())
        total = flat.sum()
        if not np.isfinite(total) or total <= 0:
            return np.array([1e-100, 1e-99])
        return np.interp(quantiles, flat.cumsum() / total, flat)

    def SaveImage(self):
        if self.params.get("_return_figure"):
            rmax = self.params["rmax"] * self._axis_unit_factor
            self.ax.set(xlim=[-rmax, rmax], ylim=[-rmax, rmax])
            self._returned_fig = self.fig
            return
        print("Saving ", self.params["filename"])
        if self.params["backend"] == "matplotlib":
            from matplotlib import pyplot as plt
            rmax = self.params["rmax"] * self._axis_unit_factor
            self.ax.set(xlim=[-rmax, rmax], ylim=[-rmax, rmax])
            # When rendering a movie the per-frame PNG dimensions and the
            # axes position within them must stay constant; bbox_inches="tight"
            # crops to the rendered content, so tick-label width changes
            # (e.g. "0.005" -> "5e-05") produce different image sizes. Opt
            # out via tight_bbox=False to preserve fixed dimensions.
            bbox = "tight" if self.params["tight_bbox"] else None
            plt.savefig(self.params["filename_incomplete"], bbox_inches=bbox, dpi=400)
            plt.close()
        elif self.params["backend"] == "PIL":
            self._pil_image.convert("RGB").save(self.params["filename_incomplete"])
            self._pil_image = None
        os.rename(self.params["filename_incomplete"], self.params["filename"])

    def MakeImages(self, snapdata):
        if not self.params["no_stars"]:
            self.AddStarsToImage(snapdata)
        if self.params["draw_axes"]:
            self.AddAxesToImage()
        self.AddSizeScaleToImage(snapdata["Header"])
        self.AddTimestampToImage(snapdata["Header"])
        self.SaveImage()

    def _save_map(self, name):
        """Write a rendered map to the .maps cache, unless caching is disabled."""
        if not self.params["no_map_cache"]:
            np.savez_compressed(self.map_files[name], **{name: self.maps[name]})

    def _orientation_triad(self):
        """World x/y/z unit vectors in the displayed camera frame (rows).

        Rotation + cubemap shuffle only: the triad shows on-screen directions,
        so the perspective LOS projection is bypassed by transforming as if
        orthographic."""
        camera_distance = self.params["camera_distance"]
        self.params["camera_distance"] = np.inf
        try:
            axes = np.eye(3)
            self.DoCoordinateTransform(axes, contravariant=True)
        finally:
            self.params["camera_distance"] = camera_distance
        return axes

    # x/y/z triad colors
    _AXIS_COLORS = ("#e05659", "#63c74d", "#4d9be6")

    def AddAxesToImage(self):
        """Draw an x/y/z orientation triad in the top-right corner."""
        axes = self._orientation_triad()
        order = np.argsort(axes[:, 2])  # back-to-front: -z is away from viewer
        labels = "xyz"

        if self.params["backend"] == "PIL":
            from PIL import ImageDraw

            W = self._pil_image.size[0]
            length = 0.06 * W
            cx, cy = W - 2.2 * length, 2.2 * length
            draw = ImageDraw.Draw(self._pil_image)
            font = _get_font(max(10, W // 24))
            for i in order:
                dx, dy = axes[i, 0], axes[i, 1]  # PIL y runs downward
                draw.line(
                    [(cx, cy), (cx + dx * length, cy - dy * length)],
                    fill=self._AXIS_COLORS[i], width=max(1, W // 256),
                )
                label_pos = (cx + 1.35 * dx * length, cy - 1.35 * dy * length)
                try:
                    draw.text(label_pos, labels[i], fill=self._AXIS_COLORS[i], font=font, anchor="mm")
                except ValueError:  # bitmap fallback fonts don't support anchor
                    draw.text(label_pos, labels[i], fill=self._AXIS_COLORS[i], font=font)
        elif self.params["backend"] == "matplotlib" and hasattr(self, "ax"):
            s = 0.07  # arrow length in axes fraction
            cx, cy = 0.87, 0.85
            for i in order:
                dx, dy = axes[i, 0] * s, axes[i, 1] * s
                self.ax.annotate(
                    "", xy=(cx + dx, cy + dy), xytext=(cx, cy),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(
                        arrowstyle="-|>", color=self._AXIS_COLORS[i], lw=1.5,
                        shrinkA=0, shrinkB=0,
                    ),
                )
                self.ax.text(
                    cx + 1.45 * dx, cy + 1.45 * dy, labels[i],
                    transform=self.ax.transAxes, color=self._AXIS_COLORS[i],
                    fontsize=9, fontweight="bold", ha="center", va="center",
                )

    def AddTimestampToImage(self, header):
        if self.params["no_timestamp"]:
            return
        # Detect cosmological simulations: Time is scale factor, show redshift
        z = header.get("Redshift", 0)
        time = self.params["Time"]
        if "Redshift" in header and abs(time * (1 + z) - 1.0) < 1e-6:
            if z < 10:
                time_text = "z=%.2g" % z
            else:
                time_text = "z=%.3g" % z
        else:
            if "UnitLength_In_CGS" in header.keys():
                unit_time_in_s = header["UnitLength_In_CGS"] / header["UnitVelocity_In_CGS"]
            else:
                unit_time_in_s = _au.kpc.to(_au.cm) / (_au.km.to(_au.cm))
            time_Myr = time * unit_time_in_s / _au.Myr.to(_au.s)
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
            import matplotlib.patheffects as path_effects
            txt = self.ax.text(
                0.05, 0.95, time_text,
                transform=self.ax.transAxes,
                color="white", fontsize=10, fontweight="bold",
                ha="left", va="top",
            )
            txt.set_path_effects([
                path_effects.Stroke(linewidth=1.5, foreground="black"),
                path_effects.Normal(),
            ])

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
            UL = header["UnitLength_In_CGS"]
        else:
            UL = _au.kpc.to(_au.cm)  # GADGET default: kpc
        r = self.params["rmax"] * UL / _au.pc.to(_au.cm)  # in pc
        if self.params["camera_distance"] < np.inf:
            r = self.params["rmax"] * self.params["camera_distance"]  # * self.params["unit_scalefac"]
        if r * 2 > 1000:
            scale_kpc = 10 ** np.round(np.log10(r * 0.5 / 1000))
            size_scale_text = "%3.3gkpc" % (scale_kpc)
            size_scale_ending = gridres / 16 + gridres * (scale_kpc * 1000) / (2 * r)
        elif r * 2 > 0.1:  # below 0.1pc across, AU reads better than fractions of a pc
            scale_pc = 10 ** np.round(np.log10(r * 0.5))
            size_scale_text = "%3.3gpc" % (scale_pc)
            size_scale_ending = gridres / 16 + gridres * (scale_pc) / (2 * r)
        else:
            scale_AU = 10 ** np.round(np.log10(r * 0.5 * pc_to_AU))
            size_scale_text = "%gAU" % (scale_AU)  # plain digits, not 1e+04
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
                    self._save_map("realstars")  # only when freshly computed, not on cache hits
                # blend realstars into the in-memory image; nan_to_num so a bad map
                # (possibly loaded from an old cache) can't NaN out the whole frame
                img_arr = np.array(self._pil_image.convert("RGB")).astype(np.float32) / 255
                img_arr = np.clip(img_arr + np.nan_to_num(self.maps["realstars"]), 0, 1)
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
            self.DoCoordinateTransform(X_star, np.ones(len(X_star)), np.ones(len(X_star)))
            # Floor the marker size instead of zeroing it, so that sinks below
            # sink_scale (e.g. protostars in a sub-solar collapse run) are still
            # visible.  The floor continues the formula down to sink_scale/10,
            # below which the size would shrink to nothing.
            star_size = np.maximum(np.log10(m_star / self.params["sink_scale"]) + 2, _MIN_STAR_SIZE)
            colors = np.array([self.GetStarColor(m) for m in m_star]) / 255

            self.ax.scatter(
                X_star[:, 0] * self._axis_unit_factor,
                X_star[:, 1] * self._axis_unit_factor,
                s=star_size * 5,
                edgecolor=self.Star_Edge_Color(),
                lw=0.1,
                facecolor=colors,
                marker="*",
            )

            if self.params["star_legend"]:
                # Stellar-mass key in the bottom-left, styled after
                # starforge_tools.plots.star_markers.plot_star_legend but using
                # CrunchSnaps' own size/color scheme so the legend markers
                # match the rendered stars.
                for m_dummy in (1, 10, 100, 1000):
                    s_dummy = max(np.log10(m_dummy / self.params["sink_scale"]) + 2, _MIN_STAR_SIZE) * 5
                    self.ax.scatter(
                        [np.inf], [np.inf],
                        s=s_dummy,
                        facecolor=np.array(self.GetStarColor(m_dummy)) / 255,
                        edgecolor=self.Star_Edge_Color(),
                        lw=0.1,
                        marker="*",
                        label=r"$%g\,M_\odot$" % m_dummy,
                    )
                ledge = self.ax.legend(
                    loc=3, frameon=True, facecolor="black",
                    labelspacing=0.1, fontsize=6, edgecolor="white",
                )
                ledge.get_frame().set_linewidth(0.5)
                ledge.get_frame().set_alpha(0.5)
                for txt in ledge.get_texts():
                    txt.set_color("white")

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

    @staticmethod
    def _default_center(snapdata):
        """Center of the simulation domain: BoxSize/2 in the box frame, else the origin."""
        boxsize = snapdata["Header"]["BoxSize"]
        coords = snapdata.get("PartType0/Coordinates")
        if coords is None or not len(coords):
            coords = snapdata.get("PartType5/Coordinates")  # sink-only snapshot
        if coords is not None and len(coords) and not coords_in_box_frame(coords, boxsize):
            return np.zeros(3)  # non-periodic runs are generally centered on the origin
        return np.repeat(boxsize * 0.5, 3)

    @staticmethod
    def _has_gas(snapdata):
        """Whether the snapshot has any gas cells to render."""
        mass = snapdata.get("PartType0/Masses")
        return mass is not None and len(mass) > 0

    def assign_center(self, snapdata):
        """Assign the center of the image"""

        center = self.params["center"]
        if isinstance(center, (np.ndarray, list, tuple)) and len(center) == 3:
            self.params["center"] = np.asarray(center, dtype=float)
            return

        if center is None:
            self.params["center"] = self._default_center(snapdata)
            return

        # the gas-based modes fall through to the default center when there is no gas
        have_gas = self._has_gas(snapdata)
        match self.params["center"]:
            case "densest" if have_gas:
                if "PartType0/Density" in snapdata and snapdata["PartType0/Density"] is not None:
                    rho = snapdata["PartType0/Density"]
                else:
                    rho = snapdata["PartType0/Masses"] / snapdata["PartType0/KernelMaxRadius"] ** 3
                self.params["center"] = snapdata["PartType0/Coordinates"][rho.argmax()]
            case "median" if have_gas:
                self.params["center"] = np.median(snapdata["PartType0/Coordinates"], axis=0)
            case x if "massive" in x and (have_gas or "PartType5/Coordinates" in snapdata.keys()):
                if "PartType5/Coordinates" in snapdata.keys():
                    if "=" in self.params["center"]:
                        num = int(self.params["center"].split("=")[1])
                    else:
                        num = 1
                    self.params["center"] = snapdata["PartType5/Coordinates"][
                        snapdata["PartType5/Sink_Mass"].argsort()[::-1]
                    ][num - 1]
                else:
                    if "PartType0/Density" in snapdata and snapdata["PartType0/Density"] is not None:
                        rho = snapdata["PartType0/Density"]
                    else:
                        rho = snapdata["PartType0/Masses"] / snapdata["PartType0/KernelMaxRadius"] ** 3
                    self.params["center"] = snapdata["PartType0/Coordinates"][rho.argmax()]
            case x if x.startswith("sink_ID="):
                # sink-only lookup: sinks inherit the ID of the gas cell they
                # formed from, so a plain "ID=" search can match a recycled gas
                # ID first and silently center somewhere else entirely
                sink_id = int(x.split("=")[1])
                ids = snapdata.get("PartType5/ParticleIDs")
                idx = np.flatnonzero(ids == sink_id) if ids is not None else []
                if not len(idx):
                    raise ValueError(f"No sink particle with ID {sink_id} in this snapshot")
                self.params["center"] = snapdata["PartType5/Coordinates"][idx[0]]
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
            self.params["center"] = self._default_center(snapdata)  # default

    def _line_of_sight(self):
        """Unit viewing direction in the original (untransformed) frame."""
        if self.params["camera_dir"] is not None:
            los = np.array(self.params["camera_dir"], dtype=float)
            return los / np.sqrt(np.dot(los, los))
        pan_rad = np.pi * self.params["pan"] / 180
        tilt_rad = np.pi * self.params["tilt"] / 180
        return np.array([
            np.sin(pan_rad) * np.cos(tilt_rad),
            np.sin(tilt_rad),
            np.cos(pan_rad) * np.cos(tilt_rad),
        ])

    def _compute_default_rmax(self, snapdata):
        """Compute default rmax from mass-weighted 2D variance in the viewing plane.

        Scaled so that a uniform-density cube returns rmax = BoxSize/2.
        """
        pos = snapdata.get("PartType0/Coordinates")
        mass = snapdata.get("PartType0/Masses")
        if pos is None or mass is None or not len(mass) or mass.sum() == 0:
            return self._sink_extent_rmax(snapdata)

        total_mass = mass.sum()
        dx = pos - self.params["center"]
        los = self._line_of_sight()

        # Var_2D = Tr(Cov) - los^T @ Cov @ los, without forming the full 3x3
        trace_cov = np.sum(mass * np.sum(dx ** 2, axis=1)) / total_mass
        los_proj = dx @ los
        var_los = np.sum(mass * los_proj ** 2) / total_mass
        var_2d = max(trace_cov - var_los, 0)

        # For uniform cube of side L: Var_2D = L^2/6, want rmax = L/2
        # => rmax = sqrt(6)/2 * sqrt(Var_2D)
        return np.sqrt(6) / 2 * np.sqrt(var_2d)

    def _sink_extent_rmax(self, snapdata):
        """Field of view for a snapshot with no gas: frame the sinks instead.

        The variance estimator above is calibrated for a continuous mass distribution and badly
        under-covers a handful of point masses, so use their actual projected radius.  Only a set
        of sinks projecting onto the center itself carries no scale at all; that falls back to
        the box.
        """
        boxsize = snapdata["Header"]["BoxSize"]
        pos = snapdata.get("PartType5/Coordinates")
        if pos is None or not len(pos):
            return boxsize / 2
        dx = np.asarray(pos) - self.params["center"]
        los = self._line_of_sight()
        r2d = np.sqrt(np.maximum(np.sum(dx ** 2, axis=1) - (dx @ los) ** 2, 0))
        rmax = 1.2 * r2d.max()  # pad so the outermost star is not sitting on the frame edge
        return rmax if rmax > 0 else boxsize / 2

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
        # always needed by MakeImages (e.g. to place sinks), even on the cached-map path
        self.SetupCameraBasis()
        if not self.has_required_maps():
            self.SetupCoordsAndWeights(snapdata)
            if self._has_gas(snapdata):
                self.GenerateMaps(snapdata)
            else:
                self.GenerateBlankMaps()
        self.MakeImages(snapdata)
        if self.params.get("_return_figure"):
            return getattr(self, "_returned_fig", None)
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

    def _add_colorbar_to_image(self, vmin, vmax, label=None, log_scale=True):
        """Overlay a horizontal colorbar in the bottom-right of the in-memory PIL image."""
        from PIL import Image, ImageDraw
        from matplotlib import pyplot as plt

        if not (np.isfinite(vmin) and np.isfinite(vmax) and vmax > vmin):
            return
        if log_scale and vmin <= 0:
            log_scale = False

        img = self._pil_image
        W, H = img.size
        cmap = plt.get_cmap(self.params["cmap"])
        font_size = max(6, W // 80)
        tick_gap = 3

        # Tick values: integer powers of 10 on a log scale, evenly spaced
        # values on a linear one (signed data, e.g. LOS components)
        if log_scale:
            log_vmin, log_vmax = np.log10(vmin), np.log10(vmax)
            tick_values = [vmin, vmax]
            for e in range(int(np.ceil(log_vmin)), int(np.floor(log_vmax)) + 1):
                tv = 10.0**e
                if tv > vmin * 1.1 and tv < vmax * 0.9:
                    tick_values.append(tv)
            tick_values.sort()
            tick_frac = lambda tv: (np.log10(tv) - log_vmin) / (log_vmax - log_vmin)
        else:
            tick_values = list(np.linspace(vmin, vmax, 5))
            tick_frac = lambda tv: (tv - vmin) / (vmax - vmin)

        def _format_tick(tv):
            if tv == 0:
                return r"$0$"
            exp = int(np.floor(np.log10(np.abs(tv))))
            if not log_scale and -2 <= exp < 4:
                return r"$%g$" % np.round(tv, 6)
            coeff = tv / 10**exp
            if abs(coeff - 1.0) < 0.01:
                return r"$10^{%d}$" % exp
            return r"$%.2g\times10^{%d}$" % (coeff, exp)

        # Measure a sample label to size the contrast-sampling region
        sample = self._render_latex_label(r"$10^{3}$", font_size, "#FFFFFF")
        tick_h = sample.size[1]

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
        font_size = max(8, W // 55)
        for tv in tick_values:
            frac = tick_frac(tv)
            if not np.isfinite(frac):
                continue
            label_text = _format_tick(tv)
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

    @staticmethod
    def _length_unit_label(header):
        """Return a human-readable label for the code length unit."""
        if "UnitLength_In_CGS" in header:
            UL = header["UnitLength_In_CGS"]
        else:
            UL = _au.kpc.to(_au.cm)  # GADGET default
        ratio = UL / _au.pc.to(_au.cm)
        if abs(ratio - 1.0) < 0.01:
            return "pc"
        elif abs(ratio - 1e3) < 10:
            return "kpc"
        elif abs(ratio - 1e6) < 1e4:
            return "Mpc"
        elif abs(ratio * _au.pc.to(_au.AU) - 1.0) < 0.01:
            return "AU"
        return "code length"

    def _set_axis_length_unit(self, header):
        """Pick the matplotlib axis unit: code units, or AU once the frame is
        under 0.1pc across, where fractions of a pc stop being readable.  Matches
        the switch in AddSizeScaleToImage, which draws the PIL backend's scale bar.

        Sets self._axis_unit_factor, by which every data coordinate handed to the
        axes must be multiplied, and returns the label for it.
        """
        self._axis_unit_factor = 1.0
        lu = self._length_unit_label(header)
        if self.params["camera_distance"] < np.inf:
            return lu  # axes are in radians; nothing to convert
        UL = header["UnitLength_In_CGS"] if "UnitLength_In_CGS" in header else _au.kpc.to(_au.cm)
        if lu != "AU" and 2 * self.params["rmax"] * UL / _au.pc.to(_au.cm) < 0.1:
            self._axis_unit_factor = UL / _au.AU.to(_au.cm)
            return "AU"
        return lu

    @staticmethod
    def _velocity_unit_label(header):
        """Return a human-readable label for the code velocity unit."""
        if "UnitVelocity_In_CGS" in header:
            UV = header["UnitVelocity_In_CGS"]
        else:
            UV = _au.km.to(_au.cm)  # GADGET default
        ratio = UV / _au.km.to(_au.cm)
        if abs(ratio - 1.0) < 0.01:
            return r"\mathrm{km\,s^{-1}}"
        ratio_ms = UV / 100.0
        if abs(ratio_ms - 1.0) < 0.01:
            return r"\mathrm{m\,s^{-1}}"
        return r"\mathrm{code\;vel}"

    @staticmethod
    def _surface_density_unit_label(header):
        """Return a label for code mass / code length^2."""
        if "UnitMass_In_CGS" in header:
            UM = header["UnitMass_In_CGS"]
        else:
            UM = 1e10 * _ac.M_sun.cgs.value  # GADGET default
        ratio = UM / _ac.M_sun.cgs.value
        if abs(ratio - 1.0) < 0.01:
            mass_label = r"M_\odot"
        elif abs(ratio - 1e10) < 1e8:
            mass_label = r"10^{10}\,M_\odot"
        else:
            mass_label = r"\mathrm{code\;mass}"
        length_label = SinkVis._length_unit_label(header)
        return r"\Sigma_{\rm gas}\;\left(" + mass_label + r"\," + length_label + r"^{-2}\right)"

    def has_required_maps(self):
        return np.all([i in self.maps for i in self.required_maps])


class SinkVisSigmaGas(SinkVis):
    color_limit_maps = {"limits": "sigma_gas"}

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
                **self.splat_kwargs,
            ).T.clip(1e-100)
            # self.maps["sigma_gas"] = self.maps["sigma_gas"]
            self._save_map("sigma_gas")

    def MakeImages(self, snapdata):
        # a snapshot with no gas gets a blank field for the stars to sit on: the map is all
        # zeros, so its color scale would be pure noise and its colorbar meaningless
        no_gas = not self._has_gas(snapdata)
        if self.params["limits"] is None:
            # Mass-weighted 1st/99th percentiles for surface density
            self.params["limits"] = self._mass_weighted_limits(self.maps["sigma_gas"])

        vmin, vmax = self.params["limits"]
        # empty pixels are exactly zero on a blank map -> -inf -> clipped to the bottom of the scale
        with np.errstate(divide="ignore"):
            if vmax > vmin:
                f = (np.log10(self.maps["sigma_gas"]) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))
            else:
                f = np.zeros_like(self.maps["sigma_gas"])

        import matplotlib
        from matplotlib import pyplot as plt

        if self.params["backend"] == "PIL":
            from PIL import Image
            if no_gas:  # a blank field reads better than the bottom of the colormap
                self._pil_image = Image.new("RGBA", (f.shape[1], f.shape[0]), (0, 0, 0, 255))
            else:
                # NaN pixels (no gas) pass through the colormap as the transparent
                # bad color; the cast warnings that triggers are expected
                with np.errstate(invalid="ignore", over="ignore"):
                    rgba = plt.get_cmap(self.params["cmap"])(np.flipud(f))
                    self._pil_image = Image.fromarray((rgba * 255).astype(np.uint8), "RGBA")
            if not self.params["no_colorbar"] and not no_gas:
                self._add_colorbar_to_image(
                    vmin, vmax,
                    label=self._surface_density_unit_label(snapdata["Header"]),
                )
        elif self.params["backend"] == "matplotlib":
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            self.fig, self.ax = plt.subplots(figsize=(4, 4))
            # Explicit margins so the Y label, colorbar, and the colorbar's own
            # "Σ_gas" label all fit without bbox_inches="tight" needing to crop
            # to content (which would make per-frame PNG dimensions vary).
            self.fig.subplots_adjust(left=0.16, right=0.82, top=0.95, bottom=0.12)
            lu = self._set_axis_length_unit(snapdata["Header"])
            rmax = self.params["rmax"] * self._axis_unit_factor
            # imshow rather than pcolormesh: the axes is almost never an exact
            # integer number of device pixels per map cell, and Agg snaps every
            # quad of a QuadMesh to pixel boundaries, so neighbouring cells get
            # merged in an irregular pattern that shows up as blocky square
            # artifacts.  imshow resamples the array properly.  extent is also
            # the cell *edges*, so the map no longer sits half a cell off the
            # axis coordinates the way center-shaded pcolormesh does.
            if no_gas:
                self.ax.set_xlim(-rmax, rmax)
                self.ax.set_ylim(-rmax, rmax)
                self.ax.set_facecolor("black")
            else:
                p = self.ax.imshow(
                    self.maps["sigma_gas"],
                    extent=[-rmax, rmax, -rmax, rmax],
                    origin="lower",
                    interpolation="antialiased",
                    norm=matplotlib.colors.LogNorm(vmin=self.params["limits"][0], vmax=self.params["limits"][1]),
                    cmap=self.params["cmap"],
                )
                divider = make_axes_locatable(self.ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                self.fig.colorbar(p, label="$" + self._surface_density_unit_label(snapdata["Header"]) + "$", cax=cax)
            self.ax.set_aspect("equal")
            if self.params["camera_distance"] == np.inf:
                self.ax.set_xlabel(f"X ({lu})")
                self.ax.set_ylabel(f"Y ({lu})")
            else:
                self.ax.set_xlabel("X (rad)")
                self.ax.set_ylabel("Y (rad)")

        super().MakeImages(snapdata)


class SinkVisCoolMap(SinkVis):
    color_limit_maps = {"limits": "sigma_gas", "v_limits": "sigma_1D"}

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
        self.params["backend"] = "PIL"  # the HSV blend + colorbar path is PIL-only; matplotlib's is a thin fallback

    def GenerateMaps(self, snapdata):
        super().GenerateMaps(snapdata)
        from meshoid import GridSurfaceDensity

        need_sigma = "sigma_gas" not in self.maps
        need_disp = "sigma_1D" not in self.maps
        if need_disp:
            # need to apply coordinate transforms to z-velocity
            v = np.copy(snapdata["PartType0/Velocities"])
            if hasattr(self, "_keep_mask"):
                v = v[self._keep_mask]
            self.DoCoordinateTransform(v, contravariant=True)

        # when all three gas maps must be rendered, deposit them in one
        # particle traversal (kernel evaluated once) if meshoid supports it
        GridSurfaceDensityMulti = None
        if need_sigma and need_disp:
            try:
                from meshoid import GridSurfaceDensityMulti
            except ImportError:
                pass  # older meshoid: fall back to one splat per map

        splat_args = dict(res=self.params["res"], parallel=self.parallel, **self.splat_kwargs)
        if GridSurfaceDensityMulti is not None:
            stacked = GridSurfaceDensityMulti(
                np.column_stack([self.mass, self.mass * v[:, 2], self.mass * v[:, 2] ** 2]),
                self.pos,
                self.hsml,
                np.zeros(3),
                2 * self.params["rmax"],
                **splat_args,
            )
            # slice the channel BEFORE transposing: .T on the 3D array would
            # reverse all three axes and scramble channels into the y-axis
            self.maps["sigma_gas"] = stacked[:, :, 0].T
            mom1, mom2 = stacked[:, :, 1].T, stacked[:, :, 2].T
            self._save_map("sigma_gas")
        else:
            if need_sigma:
                self.maps["sigma_gas"] = GridSurfaceDensity(
                    self.mass,
                    self.pos,
                    self.hsml,
                    np.zeros(3),
                    2 * self.params["rmax"],
                    **splat_args,
                ).T
                self._save_map("sigma_gas")
            if need_disp:
                mom2 = GridSurfaceDensity(
                    self.mass * v[:, 2] ** 2,
                    self.pos,
                    self.hsml,
                    np.zeros(3),
                    2 * self.params["rmax"],
                    **splat_args,
                ).T
                mom1 = GridSurfaceDensity(
                    self.mass * v[:, 2],
                    self.pos,
                    self.hsml,
                    np.zeros(3),
                    2 * self.params["rmax"],
                    **splat_args,
                ).T
        if need_disp:
            # pixels with no gas divide 0/0; NaN marks them as empty.  The
            # maximum() guards against <v^2> - <v>^2 rounding negative.
            with np.errstate(invalid="ignore", divide="ignore"):
                sigma_1D = mom2 / self.maps["sigma_gas"]
                v_avg = mom1 / self.maps["sigma_gas"]
                dispersion_sq = np.maximum(sigma_1D - v_avg**2, 0)
            # convert code velocity to km/s; snapshots without unit metadata
            # are assumed to follow the STARFORGE km/s convention
            _unit_ov = self.params.get("_unit_overrides") or {}
            uv = _unit_ov.get("UnitVelocity_In_CGS", snapdata["Header"].get("UnitVelocity_In_CGS", 1e5))
            v_to_kms = uv / 1e5
            self.maps["sigma_1D"] = np.sqrt(dispersion_sq) * v_to_kms
            self._save_map("sigma_1D")

    def MakeImages(self, snapdata):
        if self.params["limits"] is None:
            # Mass-weighted 1st/99th percentiles for surface density
            self.params["limits"] = self._mass_weighted_limits(self.maps["sigma_gas"])
        if self.params["v_limits"] is None:
            # Energy-weighted percentiles: pixels with more kinetic energy count more.
            # Empty pixels are NaN in sigma_1D and would poison the cumsum.
            s1_flat = self.maps["sigma_1D"].ravel()
            Ekin = (self.maps["sigma_gas"].ravel() * s1_flat**2)
            good = np.isfinite(s1_flat) & np.isfinite(Ekin)
            s1_flat, Ekin = s1_flat[good], Ekin[good]
            if s1_flat.size and Ekin.sum() > 0:
                order = s1_flat.argsort()
                cw = Ekin[order].cumsum() / Ekin.sum()
                vmin, vmax = np.interp([0.0, 0.99], cw, s1_flat[order])
            else:
                vmax = 0.0  # no kinetic energy anywhere in frame
            if vmax <= 0:
                vmin, vmax = 1e-3, 1.0  # arbitrary but finite; e.g. zero-dispersion scenes
            elif vmin <= 0:
                vmin = 1e-3 * vmax  # log color scale needs a positive floor
            self.params["v_limits"] = np.array([vmin, vmax])
        # log10(0) -> -inf is fine here: empty pixels clip to fgas=0 (black) and
        # zero-dispersion pixels map to the bottom of the colormap
        with np.errstate(invalid="ignore", divide="ignore"):
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

    def blank_map_shape(self, name):
        if name == "SHO_RGB":
            return self.params["res"], self.params["res"], 3
        return super().blank_map_shape(name)

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
        # like CoolMap, this renders a finished RGB array rather than a colormapped scalar, so
        # there is no matplotlib axes for the star markers to be drawn on
        self.params["backend"] = "PIL"

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

            self._save_map("SHO_RGB")

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
            # a channel with no signal at all would otherwise normalize 0/0 to NaN
            norm = [n if n > 0 else 1.0 for n in norm]
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

# Helium mass fraction assumed when the snapshot does not track He (GIZMO's
# HYDROGEN_MASSFRAC = 0.76)
_PRIMORDIAL_HELIUM_MASSFRAC = 0.24

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
# For callable derived fields, stores the list of snapshot fields they may need.
DERIVED_FIELD_DEPS = {}

# Fallback expressions for snapshot fields that may not exist.
# Used only when PartType0/<name> is absent from the snapshot.
FIELD_FALLBACKS = {}


def register_derived_field(name, expr, deps=None):
    """Register a named derived field.

    expr may be a string expression or a callable(snapdata, _cache, _extra_ns).
    When expr is callable, deps lists the snapshot fields it may access so that
    _extract_field_names can request them for pre-loading.

    >>> register_derived_field("MagneticPressure", "norm(MagneticField)**2 / (8*pi)")
    """
    DERIVED_FIELDS[name] = expr
    if callable(expr):
        DERIVED_FIELD_DEPS[name] = list(deps) if deps else []


def register_field_fallback(name, expr):
    """Register a fallback expression for a snapshot field.

    The fallback is used only when PartType0/<name> is not present in the
    snapshot.  If the field exists on disk, the snapshot value is used.

    >>> register_field_fallback("Pressure", "(5./3 - 1) * Density * InternalEnergy")
    """
    FIELD_FALLBACKS[name] = expr


def _number_density_field(snapdata, _cache, _extra_ns, unit_overrides=None):
    density = _eval_field_expr("Density", snapdata, _cache, unit_overrides, _extra_ns)
    header = snapdata.get("Header", {})
    ov = unit_overrides or {}
    unit_mass = ov.get("UnitMass_In_CGS", header.get("UnitMass_In_CGS", 1.989e33))      # default: 1 M_sun
    unit_length = ov.get("UnitLength_In_CGS", header.get("UnitLength_In_CGS", 3.086e18))  # default: 1 pc
    unit_density = unit_mass / unit_length**3
    key = "PartType0/Metallicity"
    if key in snapdata and snapdata[key] is not None:
        met = np.asarray(snapdata[key])
        if met.ndim < 2:  # NUM_METAL_SPECIES=1: GIZMO writes a rank-1 dataset
            met = met[:, None]
        if met.shape[1] > 1:
            # Metallicity[:,0] = total metal fraction (Z), [:,1] = He fraction (Y)
            xH = 1.0 - met[:, 0] - met[:, 1]
        else:
            # only total Z is tracked, so assume primordial He
            xH = 1.0 - met[:, 0] - _PRIMORDIAL_HELIUM_MASSFRAC
    else:
        xH = 1.0
    return xH * density * unit_density / _CONSTANTS["m_p"]


# GIZMO writes two different quantities to the dataset PartType0/CoolingRate
# depending on which flag is set, and records neither flag in the snapshot:
#   OUTPUT_COOLRATE_DETAIL: Lambda, the cooling rate per n_H^2 in cgs
#                           (erg cm^3 s^-1), positive when the gas is cooling
#   OUTPUT_COOLRATE:        du/dt in code units (specific, signed)
# "auto" tells them apart by the companion rate fields, which only
# OUTPUT_COOLRATE_DETAIL writes.  Set to "lambda" or "dudt" to force one.
COOLRATE_CONVENTION = "auto"

# Written under the same #if as CoolingRate by OUTPUT_COOLRATE_DETAIL
_COOLRATE_DETAIL_COMPANIONS = ("HeatingRate", "NetHeatingRateQ", "HydroHeatingRate")

# t_cool = e/|de/dt| in seconds.  Lambda*n_H^2 is a volumetric rate, so the
# energy must be volumetric too: e = rho_cgs * u_cgs.
_COOLING_TIME_LAMBDA_EXPR = (
    "InternalEnergy * UnitVelocity_In_CGS**2 * Density * UnitDensity_In_CGS"
    " / abs(CoolingRate * nH**2)"
)
# du/dt is already specific, so only the code-to-cgs time conversion is needed
_COOLING_TIME_DUDT_EXPR = "InternalEnergy * UnitTime_In_CGS / abs(CoolingRate)"


def _cooling_time_field(snapdata, _cache, _extra_ns, unit_overrides=None):
    """Cooling time in seconds, t_cool = e_thermal / |de/dt|.

    Requires PartType0/CoolingRate; the interpretation of that dataset
    depends on :data:`COOLRATE_CONVENTION` (see above).  The rate is used in
    absolute value, so this is the thermal timescale whether the gas is
    cooling or being heated, and cells with zero rate come back as inf.
    """
    if snapdata.get("PartType0/CoolingRate") is None:
        raise KeyError(
            "'CoolingTime' requires PartType0/CoolingRate, which this snapshot "
            "does not have - rerun GIZMO with OUTPUT_COOLRATE_DETAIL (preferred) "
            "or OUTPUT_COOLRATE"
        )

    convention = COOLRATE_CONVENTION
    if convention == "auto":
        has_detail = any(
            snapdata.get("PartType0/" + f) is not None for f in _COOLRATE_DETAIL_COMPANIONS
        )
        convention = "lambda" if has_detail else "dudt"
    if convention == "lambda":
        expr = _COOLING_TIME_LAMBDA_EXPR
    elif convention == "dudt":
        expr = _COOLING_TIME_DUDT_EXPR
    else:
        raise ValueError(
            f"COOLRATE_CONVENTION must be 'auto', 'lambda', or 'dudt', not {COOLRATE_CONVENTION!r}"
        )

    with np.errstate(divide="ignore", invalid="ignore"):
        return _eval_field_expr(expr, snapdata, _cache, unit_overrides, _extra_ns)


# Built-in derived fields (Gaussian CGS conventions matching GIZMO)
register_derived_field("MagneticPressure", "norm(MagneticField)**2 / (8*pi)")
register_derived_field("PlasmaBeta", "Pressure / MagneticPressure")
register_derived_field("AlfvenSpeed", "norm(MagneticField) / sqrt(4*pi*Density)")
register_derived_field("MachNumber", "norm(Velocities) / SoundSpeed")
register_derived_field("JeansLength", "SoundSpeed / sqrt(G * Density)")
register_derived_field("ThermalEnergy", "Masses * InternalEnergy")
register_derived_field("KineticEnergy", "0.5 * Masses * norm(Velocities)**2")
register_derived_field("MagneticEnergy", "norm(MagneticField)**2 / (8*pi) * Masses / Density")
register_derived_field("NumberDensity", _number_density_field, deps=["Density", "Metallicity"])
register_derived_field("nH", "NumberDensity")
register_derived_field("n_H", "NumberDensity")
register_derived_field(
    "CoolingTime",
    _cooling_time_field,
    deps=["CoolingRate", "InternalEnergy", "Density", "nH"] + list(_COOLRATE_DETAIL_COMPANIONS),
)
register_derived_field("t_cool", "CoolingTime")
register_derived_field("dx", "cbrt(Masses/Density)")
register_derived_field("DivergenceError", "abs(DivergenceOfMagneticField)*cbrt(Masses/Density)/norm(MagneticField)")

# LaTeX symbols for colorbar labels.  Keys can be snapshot field names,
# derived field names, or full expression strings.
FIELD_SYMBOLS = {}


def register_field_symbol(name, latex, unit_func=None):
    r"""Register a LaTeX symbol for a field or expression.

    Parameters
    ----------
    name : str
        Field name or expression.
    latex : str
        LaTeX symbol (without $ delimiters).
    unit_func : callable, optional
        A function ``f(header) -> str`` returning a LaTeX unit string
        (e.g. ``r"\mathrm{g\,cm^{-3}}"``).  If provided, the unit is
        appended in parentheses.

    >>> register_field_symbol("Temperature", r"T", lambda h: r"\mathrm{K}")
    """
    FIELD_SYMBOLS[name] = latex
    if unit_func is not None:
        _FIELD_UNIT_FUNCS[name] = unit_func


_FIELD_UNIT_FUNCS = {}


def _density_unit(header):
    """Return density unit label from header."""
    if "UnitMass_In_CGS" in header:
        UM = header["UnitMass_In_CGS"]
    else:
        UM = 1e10 * _ac.M_sun.cgs.value
    if "UnitLength_In_CGS" in header:
        UL = header["UnitLength_In_CGS"]
    else:
        UL = _au.kpc.to(_au.cm)
    ratio_m = UM / _ac.M_sun.cgs.value
    ratio_l = UL / _au.pc.to(_au.cm)
    if abs(ratio_m - 1) < 0.01 and abs(ratio_l - 1) < 0.01:
        return r"M_\odot\,\mathrm{pc}^{-3}"
    if abs(ratio_m - 1e10) < 1e8 and abs(ratio_l / 1e3 - 1) < 0.01:
        return r"10^{10}\,M_\odot\,\mathrm{kpc}^{-3}"
    return r"\mathrm{code\;density}"


def _velocity_unit(header):
    if "UnitVelocity_In_CGS" in header:
        UV = header["UnitVelocity_In_CGS"]
    else:
        UV = _au.km.to(_au.cm)
    if abs(UV / _au.km.to(_au.cm) - 1) < 0.01:
        return r"\mathrm{km\,s^{-1}}"
    if abs(UV / 100.0 - 1) < 0.01:
        return r"\mathrm{m\,s^{-1}}"
    return r"\mathrm{code\;vel}"


# Built-in symbols
register_field_symbol("Density", r"\rho", _density_unit)
register_field_symbol("Temperature", r"T", lambda h: r"\mathrm{K}")
register_field_symbol("Pressure", r"P")
register_field_symbol("InternalEnergy", r"u")
register_field_symbol("Masses", r"M")
register_field_symbol("SoundSpeed", r"c_s", _velocity_unit)
register_field_symbol("MagneticPressure", r"P_B")
register_field_symbol("PlasmaBeta", r"\beta")
register_field_symbol("AlfvenSpeed", r"v_A")
register_field_symbol("MachNumber", r"\mathcal{M}")
register_field_symbol("JeansLength", r"\lambda_J")
register_field_symbol("NumberDensity", r"n_\mathrm{H}\;\mathrm{(cm^{-3})}")
register_field_symbol("nH", r"n_\mathrm{H}\;\mathrm{(cm^{-3})}")
register_field_symbol("n_H", r"n_\mathrm{H}\;\mathrm{(cm^{-3})}")
register_field_symbol("ThermalEnergy", r"E_\mathrm{th}")
register_field_symbol("KineticEnergy", r"E_\mathrm{kin}")
register_field_symbol("MagneticEnergy", r"E_B")
register_field_symbol("CoolingTime", r"t_\mathrm{cool}\;\mathrm{(s)}")
register_field_symbol("t_cool", r"t_\mathrm{cool}\;\mathrm{(s)}")
register_field_symbol("CoolingTime/Myr", r"t_\mathrm{cool}\;\mathrm{(Myr)}")
register_field_symbol("CoolingTime/yr", r"t_\mathrm{cool}\;\mathrm{(yr)}")
register_field_symbol("Entropy", r"P/\rho^{\gamma}")
register_field_symbol("Sigma1D", r"\sigma_\mathrm{1D}")
register_field_symbol("PhotonEnergyDensity_EUV", r"u_\mathrm{EUV}\;\mathrm{(eV\,cm^{-3})}")
register_field_symbol("PhotonEnergyDensity_FUV", r"u_\mathrm{FUV}\;\mathrm{(eV\,cm^{-3})}")
register_field_symbol("PhotonEnergyDensity_NUV", r"u_\mathrm{NUV}\;\mathrm{(eV\,cm^{-3})}")
register_field_symbol("PhotonEnergyDensity_ONIR", r"u_\mathrm{ONIR}\;\mathrm{(eV\,cm^{-3})}")
register_field_symbol("PhotonEnergyDensity_FIR", r"u_\mathrm{FIR}\;\mathrm{(eV\,cm^{-3})}")
register_field_symbol("G0", r"G_0")


# Built-in fallbacks for fields that GIZMO may or may not write
# Pressure fallback stays in code units (Density * velocity^2)
register_field_fallback("Pressure", "(5./3 - 1) * Density * InternalEnergy")
# SoundSpeed: InternalEnergy is in code velocity^2, cs is in code velocity
register_field_fallback("SoundSpeed", "sqrt(5./3 * (5./3 - 1) * InternalEnergy)")
# Temperature must be in Kelvin: convert InternalEnergy from code to CGS first
register_field_fallback("Temperature", "(5./3 - 1) * InternalEnergy * UnitVelocity_In_CGS**2 * m_p / k_B")
register_field_fallback("DivergenceOfMagneticField", "Div(MagneticField)")


# Regex to extract identifier tokens, skipping the 'e'/'E' in scientific notation
_TOKEN_RE = re.compile(r"(?<!\d)[A-Za-z_]\w*")

# Tokens that are builtins, NOT field names
_UNIT_NAMES = {
    "UnitLength_In_CGS", "UnitMass_In_CGS", "UnitVelocity_In_CGS",
    "UnitEnergyDensity_In_CGS", "UnitDensity_In_CGS",
    "UnitTime_In_CGS", "UnitEnergy_In_CGS", "UnitB_In_Gauss",
}

_EXPR_BUILTINS = {
    "np",
    "abs", "sqrt", "cbrt", "norm", "col", "log", "log2", "log10", "exp",
    "sin", "cos", "tan", "minimum", "maximum", "clip", "where",
    "Div", "Curl", "LOS",
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
                df = DERIVED_FIELDS[name]
                sub = DERIVED_FIELD_DEPS.get(name, []) if callable(df) else re.findall(r"[A-Za-z_]\w*", df)
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


def _op_needs_meshoid(name):
    def _op(*args, **kwargs):
        raise RuntimeError(
            f"'{name}()' requires particle positions; use it inside a "
            f"Slice/SurfaceDensity/Projection/ProjectedAverage render task."
        )
    return _op


def _op_needs_camera(name):
    def _op(*args, **kwargs):
        raise RuntimeError(
            f"'{name}()' requires a camera frame; use it inside a "
            f"Slice/SurfaceDensity/Projection/ProjectedAverage render task."
        )
    return _op


def _call_derived_callable(df, snapdata, _cache, unit_overrides, _extra_ns):
    """Invoke a callable derived field.

    The protocol is ``df(snapdata, _cache, _extra_ns)``; callables that also
    declare a ``unit_overrides`` parameter additionally receive the CLI unit
    overrides, which they need if their result depends on code-unit scalings.
    """
    if "unit_overrides" in _inspect.signature(df).parameters:
        return df(snapdata, _cache, _extra_ns, unit_overrides=unit_overrides)
    return df(snapdata, _cache, _extra_ns)


def _mask_snapdata(snapdata, mask):
    """Return a shallow copy of snapdata with PartType0 particle arrays filtered by mask."""
    n_all = mask.size
    result = {}
    for k, v in snapdata.items():
        if (k.startswith("PartType0/") and v is not None
                and isinstance(v, np.ndarray)
                and v.ndim >= 1 and len(v) == n_all):
            result[k] = v[mask]
        else:
            result[k] = v
    return result


def _eval_field_expr(expr, snapdata, _cache=None, unit_overrides=None, _extra_ns=None):
    """Evaluate a field expression against loaded PartType0 snapshot data.

    Field names (e.g. 'Masses', 'Temperature') are resolved to their
    PartType0 arrays.  Derived fields are evaluated recursively and cached
    within the call tree.  Numpy ufuncs and physical constants are available.

    Pass ``_extra_ns`` to inject additional names (e.g. ``Div`` and ``Curl``
    bound to a :class:`meshoid.Meshoid` instance).
    """
    if _cache is None:
        _cache = {}

    def _norm(x):
        return np.sqrt(np.sum(np.asarray(x) ** 2, axis=-1))

    def _col(arr, i):
        return np.asarray(arr)[:, int(i)]

    ns = {
        "np": np, "abs": np.abs, "sqrt": np.sqrt, "cbrt": np.cbrt,
        "norm": _norm, "col": _col,
        "log": np.log, "log2": np.log2, "log10": np.log10,
        "exp": np.exp, "sin": np.sin, "cos": np.cos, "tan": np.tan,
        "minimum": np.minimum, "maximum": np.maximum,
        "clip": np.clip, "where": np.where,
        "Div": _op_needs_meshoid("Div"),
        "Curl": _op_needs_meshoid("Curl"),
        "LOS": _op_needs_camera("LOS"),
    }
    ns.update(_CONSTANTS)
    if _extra_ns:
        ns.update(_extra_ns)

    # Add code-unit conversion factors from snapshot header, with CLI overrides
    header = snapdata.get("Header", {})
    ov = unit_overrides or {}
    # GADGET defaults: kpc, 10^10 Msun, km/s
    _default_UL = _au.kpc.to(_au.cm)
    _default_UM = 1e10 * _ac.M_sun.cgs.value
    _default_UV = _au.km.to(_au.cm)
    UL = ov.get("UnitLength_In_CGS", header.get("UnitLength_In_CGS", _default_UL))
    UM = ov.get("UnitMass_In_CGS", header.get("UnitMass_In_CGS", _default_UM))
    UV = ov.get("UnitVelocity_In_CGS", header.get("UnitVelocity_In_CGS", _default_UV))
    ns["UnitLength_In_CGS"] = UL
    ns["UnitMass_In_CGS"] = UM
    ns["UnitVelocity_In_CGS"] = UV
    ns["UnitEnergyDensity_In_CGS"] = UM * UV**2 / UL**3
    ns["UnitDensity_In_CGS"] = UM / UL**3
    ns["UnitTime_In_CGS"] = UL / UV
    ns["UnitEnergy_In_CGS"] = UM * UV**2
    ns["UnitB_In_Gauss"] = ov.get("UnitB_In_Gauss", header.get("UnitB_In_Gauss", 1.0))

    tokens = _TOKEN_RE.findall(expr)
    for name in set(tokens) - _EXPR_BUILTINS:
        if name in ns:
            continue
        if name in _cache:
            ns[name] = _cache[name]
            continue

        # 1) Explicit derived field — always computed from expression or callable
        if name in DERIVED_FIELDS:
            df = DERIVED_FIELDS[name]
            if callable(df):
                _cache[name] = _call_derived_callable(df, snapdata, _cache, unit_overrides, _extra_ns)
            else:
                _cache[name] = _eval_field_expr(df, snapdata, _cache, unit_overrides, _extra_ns)
            ns[name] = _cache[name]
            continue

        # 2) Snapshot field — use if present
        key = "PartType0/" + name
        if key in snapdata and snapdata[key] is not None:
            ns[name] = snapdata[key]
            continue

        # 3) Fallback expression — used when snapshot field is missing
        if name in FIELD_FALLBACKS:
            _cache[name] = _eval_field_expr(
                FIELD_FALLBACKS[name], snapdata, _cache, unit_overrides, _extra_ns
            )
            ns[name] = _cache[name]
            continue

        raise KeyError(
            f"'{name}' is not a snapshot field (PartType0/{name}), "
            f"a derived field, a fallback, or a builtin"
        )
    return eval(expr, {"__builtins__": {}}, ns)


def _expr_needs_meshoid(expr, snapdata=None, _seen=None):
    """Return True if expr (after expanding derived fields and fallbacks) contains Div or Curl.

    When snapdata is provided, fallback expressions are only followed for fields
    that are actually absent from the snapshot — if the snapshot already has the
    field, no Meshoid is needed to compute it.
    """
    if re.search(r'\b(Div|Curl)\b', expr):
        return True
    if _seen is None:
        _seen = set()
    for name in re.findall(r"[A-Za-z_]\w*", expr):
        if name in _seen:
            continue
        _seen.add(name)
        if name in DERIVED_FIELDS:
            df = DERIVED_FIELDS[name]
            if not callable(df) and _expr_needs_meshoid(df, snapdata, _seen):
                return True
        elif name in FIELD_FALLBACKS:
            # Only follow the fallback if the snapshot doesn't supply the field directly
            if snapdata is None or f"PartType0/{name}" not in snapdata:
                if _expr_needs_meshoid(FIELD_FALLBACKS[name], snapdata, _seen):
                    return True
    return False


def parse_custom_task(spec):
    """Parse a task spec like 'SurfaceDensity(Masses*Temperature)'.

    Returns (render_mode, field_expr) or None if not a custom task.
    """
    m = re.match(r"^(SurfaceDensity|Projection|ProjectedAverage|WeightedVariance|Sigma1D|Slice)\((.+)\)$", spec)
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

    Expressions may wrap any (N,3) vector field in LOS() to take its signed
    line-of-sight component in the camera frame (positive = receding), e.g.
    ProjectedAverage(LOS(Velocities)) or SurfaceDensity(Masses*LOS(MagneticField)).
    """

    def __init__(self, params):
        self._render_mode = params["_render_mode"]
        self._field_expr = params["_field_expr"]
        _safe = re.sub(r"[^\w]", "_", self._field_expr)
        self._map_key = f"{self._render_mode}_{_safe}"
        self.required_maps = set([self._map_key])
        self.color_limit_maps = {"limits": self._map_key}
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

    def _los_operator(self, snapdata, full_arrays):
        """LOS(v): signed line-of-sight component of an (N,3) vector field,
        positive = receding from the camera in both camera modes.

        full_arrays says whether expressions are evaluated on the complete
        particle set (the Div/Curl path) rather than the culled one matching
        self.pos; the perspective projection then needs the transformed
        positions of that full set, computed lazily on first use.
        """
        state = {}

        def LOS(v):
            v = np.array(v, dtype=np.float64)
            if v.ndim != 2 or v.shape[1] != 3:
                raise ValueError("LOS() requires an (N,3) vector field")
            pos = None
            if full_arrays and self.params["camera_distance"] != np.inf:
                if "pos" not in state:
                    p = np.copy(snapdata["PartType0/Coordinates"]).astype(np.float64)
                    self.DoCoordinateTransform(p, update_r=False)
                    state["pos"] = p
                pos = state["pos"]
            self.DoCoordinateTransform(v, contravariant=True, pos=pos)
            if self.params["camera_distance"] == np.inf:
                # orthographic v_z is positive toward the viewer; flip so the
                # sign convention matches the perspective branch (+ = receding)
                return -v[:, 2]
            return v[:, 2]

        return LOS

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
        if self._render_mode == "Sigma1D":
            # dispersion of the LOS velocity, whatever the expression says
            self.RequiredSnapdata.append("PartType0/Velocities")

    def _colorbar_label(self, header=None):
        """Return a LaTeX string for the colorbar title."""
        # Sigma1D has its own symbol with velocity units
        if self._render_mode == "Sigma1D":
            sym = FIELD_SYMBOLS.get("Sigma1D", r"\sigma_\mathrm{1D}")
            if header:
                vel_unit = self._velocity_unit_label(header)
                sym += r"\;(" + vel_unit + ")"
            return sym

        expr = self._field_expr
        # Check for exact match on the expression or field name
        if expr in FIELD_SYMBOLS:
            inner = FIELD_SYMBOLS[expr]
        elif expr in [t for t in _TOKEN_RE.findall(expr) if t not in _EXPR_BUILTINS]:
            inner = FIELD_SYMBOLS.get(expr, r"\mathtt{" + expr.replace("_", r"\_") + "}")
        else:
            inner = r"\mathtt{" + expr.replace("_", r"\_") + "}"

        # Wrap in σ(...) for WeightedVariance mode
        if self._render_mode == "WeightedVariance":
            inner = r"\sigma(" + inner + ")"

        # Append unit label if registered and header available
        if header and expr in _FIELD_UNIT_FUNCS:
            unit = _FIELD_UNIT_FUNCS[expr](header)
            inner += r"\;(" + unit + ")"

        return inner

    def GenerateMaps(self, snapdata):
        if self._map_key in self.maps:
            return
        from meshoid import Meshoid

        # Div/Curl need the full particle neighbourhood — culling to the view
        # box removes neighbours and zeroes out kernel sums at the boundary.
        # When the expression contains them, build a Meshoid from ALL particles
        # in the transformed frame and evaluate on the full (unmasked) arrays,
        # then mask the scalar result down to the render region afterward.
        _unit_ov = self.params.get("_unit_overrides")
        if _expr_needs_meshoid(self._field_expr, snapdata):
            # Div/Curl are coordinate-invariant: compute in the original frame
            # so that field vectors (e.g. MagneticField) and position offsets
            # are in the same basis.  Applying DoCoordinateTransform to positions
            # but not to the vector field would mix frames and zero the kernel sums.
            # Use a larger des_ngb than the default so the gradient regression
            # always has enough well-distributed neighbours (avoids singular matrices).
            M_diff = Meshoid(
                snapdata["PartType0/Coordinates"],
                snapdata["PartType0/Masses"],
                des_ngb=64,
            )
            f = _eval_field_expr(
                self._field_expr, snapdata,
                unit_overrides=_unit_ov,
                _extra_ns={"Div": M_diff.Div, "Curl": M_diff.Curl,
                           "LOS": self._los_operator(snapdata, full_arrays=True)},
            )
            if f.ndim > 1:
                f = np.sqrt(np.sum(f ** 2, axis=1))
            if hasattr(self, "_keep_mask"):
                f = f[self._keep_mask]
        else:
            # No differential operators: evaluate on the already-masked
            # particle set for efficiency.
            mask = getattr(self, "_keep_mask", None)
            eval_snapdata = _mask_snapdata(snapdata, mask) if mask is not None else snapdata
            f = _eval_field_expr(
                self._field_expr, eval_snapdata, unit_overrides=_unit_ov,
                _extra_ns={"LOS": self._los_operator(eval_snapdata, full_arrays=False)},
            )
            if f.ndim > 1:
                f = np.sqrt(np.sum(f ** 2, axis=1))

        res = self.params["res"]
        rmax = self.params["rmax"]
        # Render Meshoid always uses the culled, transformed particle set
        M = Meshoid(self.pos, self.mass, self.hsml)

        from meshoid import GridSurfaceDensity

        # every Meshoid.* grid method floors the kernel radius at two pixels;
        # do the same here so bypassing them changes nothing but the backend
        hsml = np.clip(self.hsml, 4 * rmax / res, None)

        def splat(field):
            """Σ field_i W_ij over the grid of sightlines."""
            return GridSurfaceDensity(
                field, self.pos, hsml, np.zeros(3), 2 * rmax,
                res=res, parallel=self.parallel, **self.splat_kwargs,
            ).T

        weights = []  # Σ W_ij, the denominator shared by every projected average

        def projected_average(field):
            """<field> along sightlines.  Meshoid.ProjectedAverage does this in one
            traversal but has no GPU path and no CPU threading; the ratio of two
            depositions is identical to roundoff and faster on both backends."""
            if not weights:
                weights.append(splat(np.ones_like(self.mass)))
            with np.errstate(invalid="ignore", divide="ignore"):
                return splat(field) / weights[0]

        if self._render_mode == "SurfaceDensity":
            # Surface density: f is an extensive/conserved quantity (e.g. Masses)
            result = splat(f)
        elif self._render_mode == "Projection":
            # Projection: f is a volume density / intensive quantity (e.g. Density)
            # computes the line integral ∫ f dz
            result = splat(f * M.vol)
        elif self._render_mode == "ProjectedAverage":
            result = projected_average(f)
        elif self._render_mode == "WeightedVariance":
            # σ(f) = sqrt(<f²> - <f>²)
            mean_f = projected_average(f)
            mean_f2 = projected_average(f ** 2)
            result = np.sqrt(np.maximum(mean_f2 - mean_f ** 2, 0))
        elif self._render_mode == "Sigma1D":
            # Line-of-sight velocity dispersion: sqrt(<v_z²> - <v_z>²)
            # Velocities have already been coordinate-transformed by
            # SetupCoordsAndWeights, so z is the LOS direction
            v_los = np.copy(snapdata["PartType0/Velocities"])
            if hasattr(self, "_keep_mask"):
                v_los = v_los[self._keep_mask]
            self.DoCoordinateTransform(v_los, contravariant=True)
            vz = v_los[:, 2]
            mean_vz = projected_average(vz)
            mean_vz2 = projected_average(vz ** 2)
            result = np.sqrt(np.maximum(mean_vz2 - mean_vz ** 2, 0))
        elif self._render_mode == "Slice":
            # Supersample and downsample to anti-alias Voronoi edges
            ss = int(self.params.get("supersample", 2))
            # For positive quantities, slice in log space to guarantee positivity
            positive = np.all(f > 0)
            slice_f = np.log(f) if positive else f
            recon_order = int(self.params.get("order", 0))
            hi_res = M.Slice(
                slice_f, center=np.zeros(3), size=2 * rmax, res=res * ss,
                order=recon_order, slope_limiter=True,
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
        self._save_map(self._map_key)

    def MakeImages(self, snapdata):
        data = self.maps[self._map_key]
        # Fields like CoolingTime or PlasmaBeta are legitimately infinite where
        # their denominator vanishes; those pixels would otherwise poison the
        # percentiles (and the cumsum) and blank out the image.
        finite_all = data[np.isfinite(data)]
        # Signedness is a property of the field, not of one frame.  A non-negative
        # map keeps its log scale even when part of it is exactly zero, so a
        # sequence cannot flip scales the moment its empty pixels fill in.
        signed = bool(finite_all.size and finite_all.min() < 0)

        if self.params["limits"] is None:
            # Use hi-res slice data for limits if available (before AA smoothing)
            limit_data = getattr(self, "_slice_hires", data)
            finite = limit_data[np.isfinite(limit_data)]
            if not signed:
                # zeros carry no weight and would pin vmin at 0, which has no log
                # scale; take the limits from the pixels that have signal
                finite = finite[finite > 0]
            if finite.size == 0:
                self.params["limits"] = np.zeros(2)
            elif signed:
                # symmetric limits so a diverging colormap puts zero at the center
                vext = np.abs(np.percentile(finite, [0.1, 99.9])).max()
                self.params["limits"] = np.array([-vext, vext])
            elif self._render_mode in ("SurfaceDensity", "Projection"):
                # Mass-weighted 1st/99th percentiles for integral quantities
                flat = np.sort(finite.ravel())
                cw = flat.cumsum() / flat.sum()
                self.params["limits"] = np.interp([0.001, 0.999], cw, flat)
            else:
                # Raw 1st/99th percentiles for slice/projected average
                self.params["limits"] = np.percentile(finite, [0.1, 99.9])

        vmin, vmax = self.params["limits"]
        log_scale = vmin > 0 and vmax > 0
        if vmax <= vmin:
            f = np.zeros_like(data)
        elif log_scale:
            # non-positive pixels -> -inf -> clipped to the bottom of the scale;
            # NaN (no data) passes through to the colormap's transparent bad color
            with np.errstate(divide="ignore"):
                f = (np.log10(np.maximum(data, 0)) - np.log10(vmin)) / (np.log10(vmax) - np.log10(vmin))
        else:
            f = (data - vmin) / (vmax - vmin)
        f = np.clip(f, 0, 1)

        from matplotlib import pyplot as plt

        if self.params["backend"] == "PIL":
            from PIL import Image
            # NaN pixels (no gas) pass through the colormap as the transparent
            # bad color; the cast warnings that triggers are expected
            with np.errstate(invalid="ignore", over="ignore"):
                rgba = plt.get_cmap(self.params["cmap"])(np.flipud(f))
                self._pil_image = Image.fromarray((rgba * 255).astype(np.uint8), "RGBA")
            if not self.params["no_colorbar"]:
                self._add_colorbar_to_image(
                    vmin, vmax, label=self._colorbar_label(snapdata.get("Header")),
                    log_scale=log_scale,
                )
        elif self.params["backend"] == "matplotlib":
            import matplotlib
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            self.fig, self.ax = plt.subplots(figsize=(4, 4))
            self.fig.subplots_adjust(left=0.16, right=0.82, top=0.95, bottom=0.12)
            lu = self._set_axis_length_unit(snapdata["Header"])
            rmax = self.params["rmax"] * self._axis_unit_factor
            if log_scale:
                # clip so empty (zero) pixels land on the bottom color instead of
                # being masked out, matching the PIL path
                norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax, clip=True)
            else:
                norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            # imshow rather than pcolormesh - see the comment in SinkVisSigmaGas
            p = self.ax.imshow(
                data,
                extent=[-rmax, rmax, -rmax, rmax],
                origin="lower",
                interpolation="antialiased",
                norm=norm,
                cmap=self.params["cmap"],
            )
            self.ax.set_aspect("equal")
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cb_label = "$" + self._colorbar_label(snapdata.get("Header")) + "$"
            self.fig.colorbar(p, label=cb_label, cax=cax)
            if self.params["camera_distance"] == np.inf:
                self.ax.set_xlabel(f"X ({lu})")
                self.ax.set_ylabel(f"Y ({lu})")
            else:
                self.ax.set_xlabel("X (rad)")
                self.ax.set_ylabel("Y (rad)")

        super().MakeImages(snapdata)
