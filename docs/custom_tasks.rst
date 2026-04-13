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
    Computes the surface density (line-of-sight integral of the volume
    density) of the quantity given by *expr*.  For example,
    ``SurfaceDensity(Masses)`` gives :math:`\Sigma = \int \rho\, dz`
    and ``SurfaceDensity(Masses*InternalEnergy)`` gives the thermal
    energy surface density :math:`\int \rho\, u\, dz`.  Uses SPH kernel
    integration via ``meshoid.GridSurfaceDensity``.  Colormap limits are
    chosen with mass-weighted max-entropy optimization.

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
