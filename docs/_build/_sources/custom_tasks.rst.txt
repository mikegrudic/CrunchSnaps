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
    Computes the surface density of an extensive (conserved) quantity
    :math:`\int (f/V)\, dz` via ``Meshoid.SurfaceDensity``.
    *expr* should be a per-particle quantity (e.g. ``Masses``,
    ``Masses*InternalEnergy``).  For example, ``SurfaceDensity(Masses)``
    gives :math:`\Sigma = \int \rho\, dz`.  Colormap limits are chosen
    with mass-weighted percentiles.

``Projection(expr)``
    Computes the line-of-sight integral of a volume density / intensive
    quantity :math:`\int f\, dz` via ``Meshoid.Projection``.
    *expr* should be a volumetric quantity (e.g. ``Density``,
    ``NumberDensity``).  For example, ``Projection(Density)`` also gives
    :math:`\Sigma`.  Colormap limits use mass-weighted percentiles.

``ProjectedAverage(expr)``
    Computes the mass-weighted average
    :math:`\int \rho\, \mathtt{expr}\, dz \;/\; \int \rho\, dz`
    via ``Meshoid.ProjectedAverage``.  Colormap limits use uniform
    (per-pixel) weighting.

``WeightedVariance(expr)``
    Computes the mass-weighted projected variance
    :math:`\sigma(f) = \sqrt{\langle f^2 \rangle - \langle f \rangle^2}`.
    Uses two ``Meshoid.ProjectedAverage`` calls internally.  Colormap
    limits use uniform weighting.

``Sigma1D(expr)``
    Line-of-sight velocity dispersion.  Computes
    :math:`\sigma_\mathrm{1D} = \sqrt{\langle v_z^2 \rangle - \langle v_z \rangle^2}`
    where :math:`v_z` is the line-of-sight velocity after camera
    coordinate transformation.  The *expr* argument is unused.

``Slice(expr)``
    Evaluates *expr* on the midplane (z = 0 in the camera frame) using
    ``Meshoid.Slice`` with order-1 linear reconstruction.  Positive
    quantities are reconstructed in log space to guarantee positivity.
    Anti-aliased via supersampling (default 2x, set with
    ``--supersample``).  Colormap limits
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
