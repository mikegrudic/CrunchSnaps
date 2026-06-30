Custom Tasks In Depth
=====================

The custom task system lets you visualize any quantity computable from
snapshot fields without writing a new task class.

Expression Evaluation
---------------------

Expressions are evaluated per-particle against the loaded ``PartType0`` data.
The result is then projected onto a 2D grid using the selected render mode.

Resolution order for names in an expression:

1. **Derived fields** -- looked up in the ``DERIVED_FIELDS`` registry and
   evaluated recursively.
2. **Snapshot fields** -- loaded from ``PartType0/<name>`` in the HDF5 file.
3. **Field fallbacks** -- if the snapshot field is missing, the
   ``FIELD_FALLBACKS`` expression is evaluated instead.
4. **Error** -- if none of the above match.

Results are cached within each frame, so referencing the same derived field
in multiple places does not cause redundant computation.

Render Modes
------------

``SurfaceDensity(expr)``
    Computes :math:`\int \mathtt{expr}\, dz` along each sightline using
    SPH kernel integration (``meshoid.GridSurfaceDensity``).  Colormap
    limits are chosen with mass-weighted max-entropy optimization.

``Projection(expr)``
    Alias for ``SurfaceDensity``.

``ProjectedAverage(expr)``
    Computes the mass-weighted average
    :math:`\int \rho\, \mathtt{expr}\, dz \;/\; \int \rho\, dz`
    via ``Meshoid.ProjectedAverage``.  Colormap limits use uniform
    (per-pixel) weighting.

``Slice(expr)``
    Evaluates *expr* on the midplane (z = 0 in the camera frame) using
    ``Meshoid.Slice``.  The slice is rendered at 4x supersampling and
    downsampled to anti-alias Voronoi cell boundaries.  Colormap limits
    use uniform weighting.

Colormap Behavior
-----------------

When ``--limits`` is not specified, limits are chosen automatically using
maximum-entropy optimization:

- **Surface density / projection** tasks use mass-weighted entropy
  (pixels with more accumulated mass count more).
- **Slice / projected average** tasks use uniform-weighted entropy
  (all pixels count equally).
- **Log scale** is used automatically when all values are positive;
  otherwise linear scale is used.

Colorbar Labels
---------------

Each field can have a LaTeX symbol registered for colorbar display.  If no
symbol is registered, the raw expression is shown in monospace.

Built-in symbols include:

.. list-table::
   :header-rows: 1
   :widths: 25 25

   * - Field
     - Symbol
   * - ``Temperature``
     - :math:`T\;\mathrm{(K)}`
   * - ``Density``
     - :math:`\rho`
   * - ``Pressure``
     - :math:`P`
   * - ``PlasmaBeta``
     - :math:`\beta`
   * - ``MachNumber``
     - :math:`\mathcal{M}`
   * - ``MagneticPressure``
     - :math:`P_B`
   * - ``NumberDensity``
     - :math:`n\;\mathrm{(cm^{-3})}`

Register custom symbols::

    from CrunchSnaps.snapshot_tasks import register_field_symbol
    register_field_symbol("Masses*InternalEnergy", r"E_\mathrm{th}")

Built-in Functions
------------------

.. list-table::
   :header-rows: 1
   :widths: 20 40

   * - Name
     - Description
   * - ``abs(x)``
     - Absolute value
   * - ``sqrt(x)``
     - Square root
   * - ``cbrt(x)``
     - Cube root (:math:`x^{1/3}`)
   * - ``norm(v)``
     - Euclidean magnitude of a vector field (shape ``(N, 3)``)
   * - ``log(x)``, ``log2(x)``, ``log10(x)``
     - Natural, base-2, and base-10 logarithm
   * - ``exp(x)``
     - Exponential
   * - ``sin(x)``, ``cos(x)``, ``tan(x)``
     - Trigonometric functions
   * - ``minimum(a,b)``, ``maximum(a,b)``
     - Element-wise min/max
   * - ``clip(x,a,b)``
     - Clamp values to ``[a, b]``
   * - ``where(cond,a,b)``
     - Element-wise conditional selection

Physical constants (``pi``, ``G``, ``k_B``, ``m_p``, ``c_light``,
``Msun``, ``pc``, ``yr``, ``Myr``, ``eV``, …) are also available.
Use ``np`` for any other NumPy function (e.g. ``np.arctan2(y, x)``).

Built-in Aliases
----------------

``dx``
    Inter-particle spacing estimate, ``cbrt(Masses/Density)``.  Useful
    as a local resolution scale::

        --tasks='Slice(Div(MagneticField) * dx)'

Differential Operators
----------------------

``Div(v)`` and ``Curl(v)`` compute the SPH-kernel divergence and curl of
a vector field ``v`` (shape ``(N, 3)``) using
:meth:`meshoid.Meshoid.Div` and :meth:`meshoid.Meshoid.Curl`.  They are
evaluated on the particles in the render volume and return a scalar
(shape ``(N,)``) or vector (shape ``(N, 3)``) array respectively.

These operators require the particle positions and are therefore only
available inside a rendering task (``Slice``, ``SurfaceDensity``,
``Projection``, ``ProjectedAverage``).  Using them in a bare expression
outside a render context raises ``RuntimeError``.

Example — divergence of :math:`\mathbf{B}` scaled by the local cell size::

    --tasks='Slice(Div(MagneticField) * cbrt(Masses/Density)/norm(MagneticField))'

The same expression using the ``dx`` alias::

    --tasks='Slice(Div(MagneticField) * dx/norm(MagneticField))'

Vector Fields
-------------

If an expression evaluates to a vector (shape ``(N, 3)``), the magnitude
is taken automatically.  Use ``norm()`` explicitly if combining with other
terms::

    # These are equivalent for a vector field:
    --tasks=Slice(Velocities)
    --tasks=Slice(norm(Velocities))

    # But for expressions you need norm():
    --tasks=Slice(norm(Velocities)/SoundSpeed)

Defining New Derived Fields
---------------------------

From Python::

    from CrunchSnaps.snapshot_tasks import register_derived_field

    register_derived_field("SpecificAngularMomentum",
                           "norm(cross(Coordinates, Velocities))")

Derived fields can depend on other derived fields::

    register_derived_field("MagneticPressure",
                           "norm(MagneticField)**2 / (8*pi)")
    register_derived_field("PlasmaBeta",
                           "Pressure / MagneticPressure")

Field Fallbacks
---------------

Fallbacks provide computed alternatives when a snapshot field is missing.
The snapshot value is always preferred when available::

    from CrunchSnaps.snapshot_tasks import register_field_fallback

    # Use ideal-gas EOS when Pressure isn't stored
    register_field_fallback("Pressure",
                            "(5./3 - 1) * Density * InternalEnergy")
