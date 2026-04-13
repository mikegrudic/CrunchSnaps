Camera-File Movies
==================

``make_movie_from_camerafile.py`` renders frames from a camera trajectory
file, giving full control over the camera position, direction, and distance
at each point in time.  This is the tool to use for fly-through movies,
zoom sequences, and orbit animations.

Basic Usage
-----------

::

    make_movie_from_camerafile.py trajectory.txt /path/to/output/ --make_movie

where ``trajectory.txt`` is a whitespace-delimited ASCII file defining the
camera state at each frame, and ``/path/to/output/`` is the simulation
directory containing ``snapshot_*.hdf5`` files.

Camera File Formats
-------------------

The script auto-detects the format from the number of columns:

1 column: times only
^^^^^^^^^^^^^^^^^^^^

::

    # time
    0.000
    0.001
    0.002

Uses default centering and field of view at each time.

4 columns: time, distance, pan, tilt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # time  camera_dist  pan(deg)  tilt(deg)
    0.000   100.0        0.0       0.0
    0.001    80.0       10.0       5.0
    0.002    60.0       20.0      10.0

Perspective projection from ``camera_dist`` with rotation angles.

5 columns: time, position, rmax
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # time  cx     cy     cz     rmax
    0.000   0.5    0.5    0.5    10.0
    0.001   0.5    0.5    0.5     8.0

Orthographic projection centered at (cx, cy, cz) with half-width ``rmax``.
Useful for zoom sequences.

7 columns: time, distance, position, pan, tilt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # time  dist   cx     cy     cz    pan   tilt
    0.000   50.0   0.5    0.5    0.5   0.0   0.0
    0.001   40.0   0.5    0.5    0.5  15.0   5.0

Perspective projection with explicit center and rotation.

8 columns: time, position, direction, distance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # time  cx   cy   cz   dx   dy   dz   dist
    0.000  0.5  0.5  0.5  0.0  0.0  1.0  100.0

Full camera specification with arbitrary viewing direction vector.

11 columns: time, position, direction, up, distance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # time  cx   cy   cz   dx   dy   dz   ux   uy   uz   dist
    0.000  0.5  0.5  0.5  0.0  0.0  1.0  0.0  1.0  0.0  100.0

Full camera specification with explicit up vector, for complete control
over camera roll.

Generating Camera Files
-----------------------

Camera files are plain text, easily generated with numpy::

    import numpy as np

    N = 500
    time = np.linspace(0, 5e-3, N)

    # Zoom in while orbiting
    camera_dist = 3 + 30 * np.exp(-time / 3e-4)
    pan = np.cumsum(360 / 1e-3 * (3 / camera_dist) * np.gradient(time))
    tilt = 30 * np.sin(2 * np.pi * time / 1.3e-3)

    np.savetxt("trajectory.txt", np.c_[time, camera_dist, pan, tilt])

The companion script ``stellarscape_movie.py`` generates camera trajectory
files with smooth interpolated paths for complex fly-throughs.

Map Types
---------

Select the visualization type with ``--map_type``:

``SigmaGas`` (default)
    Gas surface density.

``CoolMap``
    Surface density modulated by velocity dispersion.

``SHOMap``
    Narrowband SII/H-alpha/OIII composite.

Cubemap Rendering
-----------------

Use ``--cubemap`` to render all six faces of a cubemap at each frame,
producing files with suffixes ``forward``, ``left``, ``right``, ``up``,
``down``, ``backward``::

    make_movie_from_camerafile.py traj.txt sim/ --cubemap --res=512

This is useful for VR or 360-degree video.

Examples
--------

Simple fly-through with movie output::

    make_movie_from_camerafile.py traj.txt sim/ --make_movie --res=1024 --cmap=inferno

High-resolution with realistic stars::

    make_movie_from_camerafile.py traj.txt sim/ --realstars --res=2048 --np=4 --np_render=4

Narrowband composite::

    make_movie_from_camerafile.py traj.txt sim/ --map_type=SHOMap --make_movie

Full Option Reference
---------------------

.. code-block:: text

    Usage:
    make_movie_from_camerafile.py <camerafile> <simdir> ... [options]

    Options:
        -h --help              Show this screen.
        --map_type=<type>      Map type: SigmaGas, CoolMap, SHOMap [default: SigmaGas]
        --realstars            Render stars with realistic PSFs
        --limits=<min,max>     Surface density limits
        --no_overwrite         Skip existing files
        --res=<N>              Resolution [default: 256]
        --np=<N>               Parallel workers [default: 1]
        --np_render=<N>        Render threads per worker [default: 1]
        --cubemap              Render 6 cubemap faces
        --no_timestamp         Hide timestamp
        --no_size_scale        Hide size scale bar
        --no_colorbar          Hide colorbar
        --cmap=<name>          Colormap [default: viridis]
        --SHO_RGB_norm=<f>     Narrowband normalization [default: 0.0]
        --sparse_snaps         Reduce flicker for sparse snapshots
        --id_mask=<f>          Particle ID filter (.npy file)
        --backend=<b>          PIL or matplotlib [default: PIL]
        --outputfolder=<name>  Output directory [default: .]
        --make_movie           Encode frames into mp4
        --fps=<N>              Movie frame rate [default: 24]
