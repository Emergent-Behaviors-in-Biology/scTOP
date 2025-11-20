============
Installation
============

Requirements
============

Python 3.8 or higher is required.

Install from PyPI
=================

The easiest way to install scTOP is using pip::

    pip install sctop

This will install scTOP and all required dependencies.

Install from Source
===================

To install the latest development version from GitHub::

    git clone https://github.com/Emergent-Behaviors-in-Biology/scTOP.git
    cd scTOP
    pip install -e .

The ``-e`` flag installs the package in editable mode, which is useful for development.

Verify Installation
===================

You can verify your installation by importing the package::

    import sctop as top
    print(top.__version__)

Dependencies
============

scTOP automatically installs the following required packages:

Core Dependencies
-----------------

* **numpy**: Numerical computing
* **pandas**: Data manipulation
* **scipy**: Scientific computing (linear algebra, sparse matrices)
* **anndata**: Annotated data structures for scRNA-seq
* **scikit-learn**: Machine learning utilities (metrics, feature selection)
* **tables** (PyTables): HDF5 file support
* **numba**: JIT compilation for performance

Visualization
-------------

* **matplotlib**: Plotting library
* **seaborn**: Statistical visualizations

Utilities
---------

* **tqdm**: Progress bars
* **requests**: HTTP library for downloading remote bases

Optional Components
===================

For Jupyter Notebooks
---------------------

If you plan to use scTOP in Jupyter notebooks::

    pip install jupyter

For Development
---------------

Additional packages for development and testing::

    pip install pytest
    pip install sphinx
    pip install sphinx-autoapi

Troubleshooting
===============

PyTables Installation Issues
-----------------------------

If you encounter issues installing ``tables`` (PyTables), you may need to install HDF5 first:

On macOS with Homebrew::

    brew install hdf5

On Ubuntu/Debian::

    sudo apt-get install libhdf5-dev

On Windows, consider using Anaconda which includes HDF5.

Numba Installation Issues
--------------------------

Numba requires LLVM. If installation fails, you can install scTOP without numba, but performance will be reduced. The package will fall back to non-JIT implementations.

Memory Issues
-------------

For large datasets, ensure you have sufficient RAM. scTOP is optimized for memory efficiency with chunked processing, but large atlases (>100k cells, >20k genes) may require 16GB+ RAM.

Getting Help
============

If you encounter installation issues:

1. Check the `GitHub Issues <https://github.com/Emergent-Behaviors-in-Biology/scTOP/issues>`_
2. Make sure you're using Python 3.8+
3. Try installing in a fresh virtual environment
4. Open a new issue with your error message and system details
