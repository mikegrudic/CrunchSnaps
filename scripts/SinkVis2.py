#!/usr/bin/env python
"""
Usage:
SinkVis2.py <files> ... [options]

Options:
    -h --help                  Show this screen.
    --tasks=<task1,task2...>   List of names of the plots you want to make for each snapshot [default: SigmaGas]
    --rmax=<pc>                Maximum radius of plot window; defaults to box size/10. Note that this is FOV/2 in radians if FOV_plot is enabled
    --overwrite                Overwrite existing files if they already exist
    --FOV_plot                 Flag, if enables an image is created for an observer at the coordinates defined by c, looking in direction dir, with FOV of 2*rmax, projection options 'spherical', 'frustum', the default is 'frustum'
    --backend=<b>              matplotlib vs PIL [default: PIL]
    --dir=<x,y,z>              Coordinate direction to orient the image along - x, y, or z. It also accepts vector values [default: z] 
    --full_box                 Sets the plot to the entire box, overrides rmax
    --target_time=<f>          If set to nonzero, SinkVis will try to make a single image by interpolating from the available files [default: 0.0] 
    --c=<cx,cy,cz>             Coordinates of plot window center relative to box center [default: 0.0,0.0,0.0]
    --limits=<min,max>         Dynamic range of surface density colormap
    --Tlimits=<min,max>        Dynamic range of temperature colormap in K [default: 0,0]
    --SHO_RGB_norm=<f>         Normalization constant for narrow band plot, set automatically by default. If a vector is provided, then each channel is normalized by the correponding component [default: 0.0]
    --energy_limits=<min,max>  Dynamic range of kinetic energy colormap in code units [default: 0,0]
    --ecmap=<name>             Name of colormap to use for kinetic energy [default: viridis]
    --Tcmap=<name>             Name of colormap to use for temperature [default: inferno]
    --cmap=<name>              Name of colormap to use [default: viridis]
    --cmap_fresco=<name>       Name of colormap to use for plot_fresco_stars, defaults to same as cmap [default: same]
    --cool_cmap=<name>         Name of colormap to use for plot_cool_map, defaults to same as cmap [default: magma]
    --abundance_map=<N>        Will plot the surface density of metal species N (so P[:].Metallicity[N] in GIZMO),off by default [default: -1]
    --interp_fac=<N>           Number of interpolating frames per snapshot [default: 1]
    --np=<N>                   Number of processors to run on [default: 1]
    --np_render=<N>            Number of openMP threads to run rendering on (-1 uses all available cores divided by --np) [default: 1]
    --res=<N>                  Image resolution [default: 512]
    --v_res=<N>                Resolution for overplotted velocity field if plot_v_map is on [default: 32]
    --vector_quiver_map        If enabled the velocity map will be quivers, if not then field line maps
    --velocity_scale=<f>       Scale for the quivers when using plot_v_map with vector_quiver_map, in m/s [default: 1000]
    --arrow_color=<name>       Color of the velocity arrows if plot_v_map is enabled with vector_quiver_map, [default: white]
    --slice_height=<pc>        Calculation is only done on particles within a box of 2*slice_height size around the center (mostly for zoom-ins), no slicing if set to zero [default: 0]
    --keep_only_movie          Only the movie is saved, the images are removed at the end
    --no_movie                 Does not create a movie, only makes images (legacy, default behavior now is not to make a movie)
    --make_movie               Also makes movie
    --make_movie_only          Flag, if set only the movie is made from premade images
    --fps=<fps>                Frame per second for movie [default: 24]
    --movie_name=<name>        Filename of the output movie file without format [default: sink_movie]
    --sink_type=<N>            Particle type of sinks [default: 5]
    --sink_scale=<msun>        Sink particle mass such that apparent sink size is 1 pixel for that and all asses below [default: 0.1]
    --sink_relscale=<f>        Relative size scale of a sink particles at 10xsink_scale to the entire picture, e.g. 0.01 means these stars will be 1% of the entire plotting area, [default: 0.0025]
    --center_on_star           Center image on the N_high most massive sink particles
    --center_on_densest        Center image on the N_high sinks with the densest gas nearby
    --N_high=<N>               Number of sinks to center on using the center_on_star or center_on_densest flags [default: 1]
    --center_on_ID=<ID>        Center image on sink particle with specific ID, does not center if zero [default: 0]
    --rotating_images          If set SinkVis will create a set of images for a single snashot by rotating the system around a pre-specified rotation_axis
    --rotation_init=<f>        Rotation angle for the first image around the rotation axis [default: 0]
    --rotation_max=<f>         Rotation angle for the final image around the rotation axis [default: 6.2831853]
    --rotation_steps=<N>       Number rotational steps (i.e. number of images to be made), spanning from the initial to the maximum rotation [default: 4]
    --rotation_axis=<x,y,z>    Vector defining the axis around whic the rotated images are made [default: 0.0,1.0,0.0]
    --galunits                 Use default GADGET units
    --plot_T_map               Plots both surface density and average temperature maps
    --plot_B_map               Overplots magnetic field map on plots
    --plot_v_map               Overplots velocity map on plots
    --plot_energy_map          Plots kinetic energy map
    --plot_cool_map            Plots cool map that looks cool
    --sharpen_LIC_map          If enabled SinkVis will try to sharpen the field lines in line-integral convolution maps (i.e. when using plot_B_map). May produce artifacts in the image.
    --LIC_map_max_alpha=<f>        Sets the maximum opacity of the field line map as it is blended with the background (i.e. surface density map). [default: 0.5]
    --calculate_all_maps       Calculates all data for the pickle files, even if they won't get plotted
    --plot_fresco_stars        Plots surface density map with Hubble-like PSFs for the stars
    --plot_cool_map_fresco     Plots cool map that uses Hubble-like PSFs for the stars
    --fresco_param=<f>         Parameter that sets the vmax parameter of amuse-fresco, the larger the value the more extended stellar PSFs are [default: 0.002]
    --fresco_mass_limits=<min,max>  Parameter that determines how masses are rescaled for fresco. Stellar masses are roughly clipped between min and max values, useful to define a max as massive stars are extremely luminous and dominate the image [default: 0,0]
    --fresco_mass_rescale=<f>  Rescale masses plugged into Fresco mass-luminosity relation by raising masses to this power [default: 0.3]
    --energy_v_scale=<v0>      Scale in the weighting of kinetic energy (w=m*(1+(v/v0)^2)), [default: 1000.0]
    --outputfolder=<name>      Specifies the folder to save the images and movies to [default: .]
    --name_addition=<name>     Extra string to be put after the name of the ouput files, defaults to empty string
    --no_pickle                Flag, if set no pickle file is created to make replots faster
    --no_timestamp             Flag, if set no timestamp will be put on the images
    --no_size_scale            Flag, if set no size scale will be put on the images
    --draw_axes                Flag, if set the coordinate axes are added to the figure
    --remake_only              Flag, if set SinkVis will only used already calculated pickle files, used to remake plots
    --rescale_hsml=<f>         Factor by which the smoothing lengths of the particles are rescaled [default: 1]
    --highlight_wind=<f>       Factor by which to increase wind particle masses if you want to highlight them [default: 1]
    --smooth_center=<Ns>       If not 0 and SinkVis is supposed to center on a particle (e.g. with center_on_ID) then the center coordinates are smoothed across Ns snapshots, [default: 0]
    --disable_multigrid        Disables GridSurfaceDensityMultigrid froms meshoid, uses slower GridSurfaceDensity instead
    --extinct_stars            Flag on whether to extinct stars
    --sparse_snaps             Flag, if enabled then corrections are applied to the interpolation algorithm to make the movies from sensitive maps (e.g. SHO narrowband) less flickery
"""

#Example
# python SinkVis.py /panfs/ds08/hopkins/guszejnov/GMC_sim/Tests/200msun/MHD_isoT_2e6/output/snapshot*.hdf5 --np=24 --keep_only_movie --movie_name=200msun_MHD_isoT_2e6

from docopt import docopt
from natsort import natsorted
import numpy as np
import CrunchSnaps
from CrunchSnaps.misc_functions import *
import h5py
from os.path import exists, abspath
from glob import glob

taskdict = {"SigmaGas": CrunchSnaps.SinkVisSigmaGas, "HubbleSHO": CrunchSnaps.SinkVisNarrowbandComposite, "CoolMap": CrunchSnaps.SinkVisCoolMap}


def parse_inputs_to_jobparams(input): # parse input parameters to generate a list of job parameters
    arguments=input
    filenames = natsorted(arguments["<files>"])
    nproc = int(arguments["--np"])
    np_render = int(arguments["--np_render"])
    n_interp = int(arguments["--interp_fac"])
    tasks = arguments["--tasks"].split(",")
    N_tasks = len(tasks)    
    res = int(arguments["--res"])
    direction = arguments["--dir"]
    if ',' in arguments["--SHO_RGB_norm"]: #normalization by channel
        SHO_RGB_norm = np.array([float(c) for c in arguments["--SHO_RGB_norm"].split(',')])
    else: #same normalization constant for each channel
        SHO_RGB_norm = float(arguments["--SHO_RGB_norm"])
        
    if arguments["--limits"]:
        limits = np.array([float(c) for c in arguments["--limits"].split(',')])

    # parameters that every single task will have in common
    common_params = {i.replace("--",""): input[i] for i in input}
    del common_params["<files>"]
    for c, i in common_params.items():
        if type(i) != str: continue
        if "," in i:
            common_params[c] = np.array([float(k) for k in i.split(",")])
        elif i.replace(".","").isnumeric():
            common_params[c] = float(i)
    common_params.update({"fresco_stars": input["--plot_fresco_stars"], "res": int(input["--res"]), "limits": (limits if arguments["--limits"] else None), "no_timestamp": input["--no_timestamp"], "threads": np_render, "outputfolder": input["--outputfolder"], "SHO_RGB_norm": SHO_RGB_norm, "cool_cmap": input["--cool_cmap"], "center_on_star": int(input["--center_on_star"]), "extinct_stars": int(input["--extinct_stars"]), "sparse_snaps": input["--sparse_snaps"], "backend": input["--backend"]})


    if direction=='x':
        common_params["camera_dir"] = np.array([1.,0,0])
        common_params["camera_up"] = np.array([0,0,1.])
    elif direction=='y':
        common_params["camera_dir"] = np.array([0,1.,0])
        common_params["camera_up"] = np.array([0,0,1.])

    if input["--rmax"] is None:
        common_params["rmax"] = None
    else:
        common_params["rmax"] = float(input["--rmax"])

    N_params = len(filenames)*n_interp
    
    snaptime_dict = get_snapshot_time_dict(input["<files>"]) # get times of snapshots
    snaptime_dict_inv = {v:k for v, k in zip(snaptime_dict.values(),snaptime_dict.keys())}
    snaptimes_orig = np.array([snaptime_dict[snapnum_from_path(s)] for s in input["<files>"]])

    if n_interp>1: # get times of frames if we're doing an interpolated movie
        snaptimes = np.interp(np.arange(n_interp*len(filenames))/n_interp, np.arange(len(filenames)), snaptimes_orig)
    else:
        snaptimes = snaptimes_orig
    
    params = []
    for j in range(N_tasks):
        p = []
        for i in range(N_params):
            d = common_params.copy()
            d["Time"] = snaptimes[i]
            d["index"] = snaptime_dict_inv[snaptimes_orig[i]] * 10 + i%n_interp
            p.append(d.copy())
        params.append(p)
    
    return params
    
def main(input):
    tasks = input["--tasks"].split(",")
    tasks = [taskdict[t] for t in tasks]
    nproc = int(input["--np"])
    np_render = int(input["--np_render"])
    params = parse_inputs_to_jobparams(input)
    snaps = natsorted(input["<files>"])
    CrunchSnaps.DoTasksForSimulation(snaps, tasks=tasks, task_params=params,nproc=nproc,nthreads=np_render)

if __name__ == "__main__":
    arguments = docopt(__doc__)
    main(arguments)
    
