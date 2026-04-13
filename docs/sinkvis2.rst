SinkVis2 Command-Line Tool
==========================

``SinkVis2`` renders images and movies from GIZMO simulation snapshots.

Basic Usage
-----------

Tasks (what to plot) are given as positional arguments after the snapshot
files.  Anything ending in ``.hdf5`` is treated as a snapshot; everything
else is a task.  If no tasks are given, ``SigmaGas`` is the default.

::

    # Default surface density
    SinkVis2 output/snapshot_*.hdf5

    # Surface density + temperature slice
    SinkVis2 output/snapshot_*.hdf5 SigmaGas Slice(Temperature)

    # Just a temperature slice, higher resolution, inferno colormap
    SinkVis2 output/snapshot_*.hdf5 Slice(Temperature) --res=2048 --cmap=inferno

    # Make a movie
    SinkVis2 output/snapshot_*.hdf5 SigmaGas --make_movie --fps=30

.. note::

   Shells like zsh interpret parentheses as glob qualifiers.  Quote task
   names that contain parentheses::

       SinkVis2 output/snapshot_*.hdf5 'Slice(Temperature)'

   Or use the ``--tasks`` flag, which avoids the issue::

       SinkVis2 output/snapshot_*.hdf5 --tasks=Slice(Temperature)

Built-in Tasks
--------------

The following built-in task names are available:

``SigmaGas`` (default)
    Gas surface density map with colorbar.

``CoolMap``
    Surface density modulated by velocity dispersion in HSV space.

``HubbleSHO``
    Narrowband SII/H-alpha/OIII composite resembling Hubble SHO palette
    images.

Example::

    SinkVis2 snapshot_*.hdf5 SigmaGas CoolMap

Custom Field Tasks
------------------

Beyond the built-in tasks, you can render arbitrary field expressions using
four render modes:

``SurfaceDensity(expr)``
    Surface density of *expr* --- the line-of-sight integral of the volume
    density of the given quantity.  E.g. ``SurfaceDensity(Masses)`` gives
    mass surface density :math:`\Sigma = \int \rho\, dz`.

``Projection(expr)``
    Alias for ``SurfaceDensity``.

``ProjectedAverage(expr)``
    Mass-weighted average of *expr* along the line of sight.

``Slice(expr)``
    Midplane slice of *expr*, with order-1 linear reconstruction in log
    space (for positive quantities) and anti-aliasing via supersampling
    (default 2x, controllable with ``--supersample``).

Expressions can reference any ``PartType0`` field by name (``Masses``,
``Temperature``, ``Density``, ``Velocities``, etc.), use arithmetic
operators, and call numpy functions.

Examples::

    # Temperature slice
    SinkVis2 snapshot_*.hdf5 Slice(Temperature)

    # Column density of thermal energy
    SinkVis2 snapshot_*.hdf5 SurfaceDensity(Masses*InternalEnergy)

    # Projected average of log density
    SinkVis2 snapshot_*.hdf5 ProjectedAverage(log10(Density))

    # Multiple tasks at once
    SinkVis2 snapshot_*.hdf5 SigmaGas Slice(Temperature) Slice(PlasmaBeta)

Available Functions in Expressions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following numpy functions are available in expressions:

- ``abs``, ``sqrt``, ``norm`` (vector magnitude)
- ``log``, ``log2``, ``log10``, ``exp``
- ``sin``, ``cos``, ``tan``
- ``minimum``, ``maximum``, ``clip``, ``where``

Physical constants (CGS, from astropy) are also available:

.. list-table::
   :header-rows: 1
   :widths: 15 50

   * - Name
     - Value
   * - ``pi``
     - 3.14159...
   * - ``k_B``
     - Boltzmann constant (erg/K)
   * - ``m_p``
     - Proton mass (g)
   * - ``m_e``
     - Electron mass (g)
   * - ``G``
     - Gravitational constant (cm\ :sup:`3` g\ :sup:`-1` s\ :sup:`-2`)
   * - ``c_light``
     - Speed of light (cm/s)
   * - ``Msun``
     - Solar mass (g)
   * - ``Lsun``
     - Solar luminosity (erg/s)
   * - ``pc``
     - Parsec (cm)
   * - ``AU``
     - Astronomical unit (cm)
   * - ``yr``, ``Myr``
     - Year / megayear (s)
   * - ``eV``
     - Electron-volt (erg)

Derived Fields
^^^^^^^^^^^^^^

Named derived fields are pre-defined expressions that can be used by name in
any task.  They resolve recursively, so a derived field can depend on other
derived fields.

.. list-table::
   :header-rows: 1
   :widths: 25 45

   * - Name
     - Expression
   * - ``MagneticPressure``
     - ``norm(MagneticField)**2 / (8*pi)``
   * - ``PlasmaBeta``
     - ``Pressure / MagneticPressure``
   * - ``AlfvenSpeed``
     - ``norm(MagneticField) / sqrt(4*pi*Density)``
   * - ``MachNumber``
     - ``norm(Velocities) / SoundSpeed``
   * - ``JeansLength``
     - ``SoundSpeed / sqrt(G * Density)``
   * - ``ThermalEnergy``
     - ``Masses * InternalEnergy``
   * - ``KineticEnergy``
     - ``0.5 * Masses * norm(Velocities)**2``
   * - ``MagneticEnergy``
     - ``norm(MagneticField)**2 / (8*pi) * Masses / Density``
   * - ``NumberDensity``
     - ``Density / m_p``

Example::

    SinkVis2 snapshot_*.hdf5 Slice(PlasmaBeta) Slice(MachNumber)

Field Fallbacks
^^^^^^^^^^^^^^^

Some fields may not be present in every snapshot.  Fallback expressions are
used automatically when a field is missing:

.. list-table::
   :header-rows: 1
   :widths: 20 40

   * - Field
     - Fallback expression
   * - ``Pressure``
     - ``(5/3 - 1) * Density * InternalEnergy``
   * - ``SoundSpeed``
     - ``sqrt(5/3 * (5/3 - 1) * InternalEnergy)``
   * - ``Temperature``
     - ``(5/3 - 1) * InternalEnergy * m_p / k_B``

If the field exists in the snapshot, its stored value is used instead.

Registering Custom Fields
^^^^^^^^^^^^^^^^^^^^^^^^^

You can register your own derived fields, fallbacks, and colorbar symbols
from Python before calling ``main``::

    from CrunchSnaps.snapshot_tasks import (
        register_derived_field,
        register_field_fallback,
        register_field_symbol,
    )

    register_derived_field("CoolingTime", "InternalEnergy / CoolingRate")
    register_field_fallback("ElectronAbundance", "1.0")
    register_field_symbol("CoolingTime", r"t_{\rm cool}\;\mathrm{(s)}")

Camera and Projection
---------------------

Viewing direction
^^^^^^^^^^^^^^^^^

By default the image is projected along the z-axis.  Use ``--dir`` to
change::

    SinkVis2 snapshot_*.hdf5 --dir=x       # view along x
    SinkVis2 snapshot_*.hdf5 --dir=y       # view along y

Arbitrary directions can be specified as a vector::

    SinkVis2 snapshot_*.hdf5 --dir=1,1,0   # view along (1,1,0)

Camera up vector
^^^^^^^^^^^^^^^^

By default the up direction is inferred from the viewing direction.  Use
``--camera_up`` to set it explicitly (it will be orthogonalized against
the viewing direction automatically)::

    # Edge-on disk with rotation axis along z
    SinkVis2 snapshot_*.hdf5 --dir=0,1,0 --camera_up=0,0,1

Pan and tilt
^^^^^^^^^^^^

Apply rotation about the Y and X axes::

    SinkVis2 snapshot_*.hdf5 --pan=45 --tilt=10

Perspective projection
^^^^^^^^^^^^^^^^^^^^^^

Set a finite camera distance for perspective (non-orthographic) rendering::

    SinkVis2 snapshot_*.hdf5 --camera_distance=100

Centering
^^^^^^^^^

Control the center of the image::

    SinkVis2 snapshot_*.hdf5 --center=densest       # densest gas cell
    SinkVis2 snapshot_*.hdf5 --center=massive        # most massive star
    SinkVis2 snapshot_*.hdf5 --center=massive=3      # 3rd most massive star
    SinkVis2 snapshot_*.hdf5 --center=median          # median gas position
    SinkVis2 snapshot_*.hdf5 --center=0.5,0.5,0.5    # explicit coordinates
    SinkVis2 snapshot_*.hdf5 --center=ID=12345        # specific particle ID

When not specified, the default center is the box center, and the default
field of view (``--rmax``) is determined from the mass-weighted spatial
extent of the gas.

Movies
------

``--make_movie`` renders all frames, then stitches them into an mp4 using
ffmpeg::

    SinkVis2 snapshot_*.hdf5 --make_movie --fps=24

When ``--limits`` is not explicitly set, the movie pipeline automatically
computes **global colormap limits** that maximize the sum of per-frame
Shannon entropies across the sequence, ensuring consistent and
information-rich coloring throughout the movie.

Interpolated movies
^^^^^^^^^^^^^^^^^^^

Use ``--interp_fac`` to generate interpolated frames between snapshots for
smoother movies::

    SinkVis2 snapshot_*.hdf5 --make_movie --interp_fac=4

Freeze-frame rotation
^^^^^^^^^^^^^^^^^^^^^

Add a 360-degree rotation at specific snapshot numbers::

    SinkVis2 snapshot_*.hdf5 --make_movie --freeze_rotation=100,200

Parallelism
-----------

``--np``
    Number of worker processes for processing different snapshots in parallel.

``--np_render``
    Number of threads per rendering call.  Defaults to ``-1``, which uses
    all available cores divided by ``--np``.

Example: 4 snapshot workers, each using 4 render threads::

    SinkVis2 snapshot_*.hdf5 --np=4 --np_render=4

Stars
-----

Sink particles are rendered on top of the image by default.  Use
``--no_stars`` to hide them.

``--realstars``
    Render stars with realistic PSFs based on their luminosity and
    temperature, using the AMUSE fresco algorithm.

::

    SinkVis2 snapshot_*.hdf5 --realstars --realstars_max_lum=1e4

Full Option Reference
---------------------

.. code-block:: text

    Usage:
    SinkVis2 <args> ... [options]

    Positional arguments are split automatically: .hdf5 files are snapshots,
    everything else is a task (e.g. SigmaGas, Slice(Temperature)).
    If no tasks are given, defaults to SigmaGas.

    Options:
        -h --help                         Show this screen.
        --tasks=<task1,task2...>           Alternative to positional task args
        --rmax=<pc>                        Half-width of the field of view
        --res=<N>                          Image resolution in pixels [default: 1024]
        --cmap=<name>                      Colormap name [default: viridis]
        --limits=<min,max>                 Colormap dynamic range
        --center=<s>                       Image center [default: None]
        --dir=<x,y,z>                      Viewing direction [default: z]
        --camera_up=<x,y,z>                Camera up vector [default: None]
        --pan=<deg>                        Pan angle in degrees [default: 0]
        --tilt=<deg>                       Tilt angle in degrees [default: 0]
        --camera_distance=<D>              Camera distance for perspective [default: inf]
        --target_time=<f>                  Render single image at this simulation time
        --np=<N>                           Number of worker processes [default: 1]
        --np_render=<N>                    Render threads (-1 = auto) [default: -1]
        --interp_fac=<N>                   Interpolated frames per snapshot [default: 1]
        --supersample=<N>                  Anti-aliasing factor for Slice [default: 2]
        --make_movie                       Encode frames into an mp4 movie
        --fps=<N>                          Movie frame rate [default: 24]
        --outputfolder=<name>              Output directory [default: .]
        --backend=<b>                      PIL or matplotlib [default: PIL]
        --no_stars                         Hide sink particles
        --no_timestamp                     Hide timestamp overlay
        --no_size_scale                    Hide size scale bar
        --no_colorbar                      Hide colorbar
        --no_overwrite                     Skip existing files
        --realstars                        Realistic stellar PSFs
        --realstars_lum_exp=<f>            Luminosity exponent [default: 1.0]
        --realstars_max_lum=<f>            Max stellar luminosity [default: 1.0e3]
        --realstars_opacity=<f>            Star opacity scaling [default: 1.0]
        --rescale_hsml=<f>                 Smoothing length scale factor [default: 1]
        --freeze_rotation=<num1,num2,...>  Snapshots for freeze-frame rotation
        --equal_frame_times               Equally-spaced frame times
        --sparse_snaps                    Reduce flicker for sparse snapshots
        --outflow_only                    Show only outflowing gas
        --cool_cmap=<name>                CoolMap colormap [default: magma]
        --v_limits=<min,max>              CoolMap velocity dispersion range
        --SHO_RGB_norm=<f>                Narrowband normalization [default: 0.0]
        --id_mask=<file>                  Particle ID filter (.npy file)
