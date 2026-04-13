Installation
============

Install from source::

    git clone https://github.com/mikegrudic/CrunchSnaps
    cd CrunchSnaps
    pip install -e .

This installs the ``CrunchSnaps`` Python package and the ``SinkVis2``
command-line tool.  You may need to ensure that your pip scripts directory
is on your ``PATH``::

    export PATH="$PATH:$(python -m site --user-base)/bin"

Requirements
------------

- Python >= 3.6
- meshoid_
- numpy, scipy, matplotlib, numba, h5py, astropy
- Pillow, aggdraw, scikit-image
- docopt, natsort
- ffmpeg (optional, for ``--make_movie``)

.. _meshoid: https://github.com/mikegrudic/meshoid
