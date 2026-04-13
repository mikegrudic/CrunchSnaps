API Reference
=============

Field Registration Functions
----------------------------

.. autofunction:: CrunchSnaps.snapshot_tasks.register_derived_field

.. autofunction:: CrunchSnaps.snapshot_tasks.register_field_fallback

.. autofunction:: CrunchSnaps.snapshot_tasks.register_field_symbol

.. autofunction:: CrunchSnaps.snapshot_tasks.parse_custom_task

Expression Evaluation
---------------------

.. autofunction:: CrunchSnaps.snapshot_tasks.max_entropy_limits

.. autofunction:: CrunchSnaps.snapshot_tasks.max_entropy_limits_multi

Task Classes
------------

.. autoclass:: CrunchSnaps.snapshot_tasks.SinkVis
   :members: DoTask, SetupCoordsAndWeights, GenerateMaps, MakeImages

.. autoclass:: CrunchSnaps.snapshot_tasks.SinkVisSigmaGas
   :members:
   :show-inheritance:

.. autoclass:: CrunchSnaps.snapshot_tasks.SinkVisCoolMap
   :members:
   :show-inheritance:

.. autoclass:: CrunchSnaps.snapshot_tasks.SinkVisNarrowbandComposite
   :members:
   :show-inheritance:

.. autoclass:: CrunchSnaps.snapshot_tasks.SinkVisCustomField
   :members:
   :show-inheritance:

Main Entry Point
----------------

.. autofunction:: CrunchSnaps.CrunchSnaps.DoTasksForSimulation
