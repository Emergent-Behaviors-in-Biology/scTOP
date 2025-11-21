=========================================================================================
Single-cell Type Order Parameters: identifying cell phenotype from gene expression data
=========================================================================================

scTOP is a Python package for identifying cell phenotypes from single-cell RNA-sequencing data by projecting samples onto reference cell type bases. The package uses Hopfield-inspired order parameters to define the coordinates of cell fate space, where differentiation trajectories can be studied.

The paper explaining the theory and mathematical details of this model is `scTOP: physics-inspired order parameters for cellular identification and visualization <https://doi.org/10.1242/dev.201873>`_ by Maria Yampolskaya, Michael J. Herriges, Laertis Ikonomou, Darrell Kotton, and Pankaj Mehta.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   api/index
   tutorials
   theory

Installation
============

Install via pip::

    pip install sctop

Or from source::

    git clone https://github.com/Emergent-Behaviors-in-Biology/scTOP.git
    cd scTOP
    pip install -e .

Quick Start
===========

Basic Usage
-----------

Load a pre-made reference basis::

    import sctop as top
    
    # List available bases
    print(top.list_available_bases())
    
    # Load a basis
    basis, metadata = top.load_basis(basis_key="MCKO legacy")

Process and score your data::

    # Process raw count data
    processed_sample = top.process(raw_counts)
    
    # Score against basis
    projections = top.score(basis, processed_sample)

Create a custom basis from your own data::

    import anndata as ad
    
    # Load your annotated data
    adata = ad.read_h5ad("your_data.h5ad")
    
    # Create basis with validation
    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=100,
        do_anova=False
    )
    
    # Access the basis and metrics
    basis = results['basis']
    print(f"Accuracy: {results['metrics']['accuracy']:.3f}")

Sources for Reference Databases
================================

* `Mouse Cell Atlas <http://bis.zju.edu.cn/MCA/>`_
* `Atlas of Mouse Lung Development <https://journals.biologists.com/dev/article-abstract/148/24/dev199512/273783/A-single-cell-atlas-of-mouse-lung-development?redirectedFrom=fulltext>`_

Dependencies
============

Core Dependencies
-----------------

* NumPy
* Pandas
* SciPy
* AnnData
* scikit-learn
* tables (PyTables)

Visualization
-------------

* Matplotlib
* Seaborn

Performance
-----------

* Numba (for fast rank computation)
* tqdm (progress bars)

Optional
--------

* requests (for downloading remote bases)

Citation
========

If you use scTOP in your research, please cite:

[Citation information to be added]

License
=======

scTOP is released under the MIT License.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
