# CrunchSnaps

CrunchSnaps is a Python package for performing analysis tasks on GIZMO simulations as efficiently as possible, with an efficient IO access pattern when performing passes over all snapshots and options for parallelization. It is implemented in a very generic way: the base Task class can be *any* thing that you want to compute for some number of snapshots. But included are routines in the SinkVis class including the coordinate transforms, projection routines, and plotting routines for making the STARFORGE Project's simuilation movies.

# Installation

Run the setup.py to install the package:
``python setup.py install``

# Examples

Various scripts that call CrunchSnaps routines can be found in `scripts`, including the very flexible script `SinkVis2.py`, which has a command-line interface for generating visualizations from snapshots.
