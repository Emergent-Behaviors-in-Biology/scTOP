===========
Quick Start
===========

This guide will help you get started with scTOP quickly.

Loading Pre-made Reference Bases
=================================

scTOP provides pre-computed reference bases from published atlases.

List Available Bases
--------------------

::

    import sctop as top
    
    # See what's available
    available = top.list_available_bases()
    print(available)
    # ['MCKO legacy']

Load a Basis
------------

::

    # Load the Mouse Cell Atlas Kotton lab basis
    basis, metadata = top.load_basis(
        basis_key="MCKO legacy",
        cache_dir="./bases_cache"  # Optional: cache downloaded files
    )
    
    print(f"Loaded basis with {basis.shape[0]} genes and {basis.shape[1]} cell types")
    print(f"Cell types: {list(basis.columns)}")

The basis is a pandas DataFrame with:

* **Index**: Gene names
* **Columns**: Cell type names
* **Values**: Processed expression values representing each cell type

Scoring Samples Against a Basis
================================

Once you have a basis, you can score your samples to identify cell types.

Process Your Data
-----------------

Start with raw count data (genes × samples)::

    import pandas as pd
    
    # Your raw count matrix
    # Rows = genes, Columns = samples/cells
    raw_counts = pd.DataFrame(...)
    
    # Process: normalize → log → rank → z-score
    processed_sample = top.process(raw_counts)

Score Against Basis
-------------------

::

    # Project onto cell type basis
    projections = top.score(basis, processed_sample)
    
    # projections is a DataFrame: cell_types × samples
    # Each column shows how well that sample matches each cell type

Interpret Results
-----------------

::

    # For a single sample
    sample_scores = projections['sample_1'].sort_values(ascending=False)
    
    print("Top 5 cell type matches:")
    print(sample_scores.head(5))
    
    # Visualize top matches
    top.plot_highest(sample_scores, n=10)

Creating a Custom Basis
=======================

You can create your own reference basis from annotated scRNA-seq data.

Prepare Your Data
-----------------

Your data should be in AnnData format with cell type annotations::

    import anndata as ad
    
    # Load your annotated data
    adata = ad.read_h5ad("your_atlas.h5ad")
    
    # Check that you have cell type annotations
    print(adata.obs.columns)  # Should include your cell type column
    print(adata.obs['cell_type'].value_counts())

Create Basis with Validation
-----------------------------

::

    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',  # Column name in adata.obs
        threshold=100,                 # Minimum cells per type
        test_size=0.2,                 # 20% held out for testing
        random_state=42,
        do_anova=False,                # Optional: feature selection
        cv_folds=None,                  # Or use cross-validation
        plot_results=True              # Plot performance summary
    )
    
    # Access results
    basis = results['basis']

Review Performance
------------------

::

    # Per cell type accuracy
    per_type = results['per_cell_type']
    print("\nWorst performing cell types:")
    print(per_type.nsmallest(10, 'accuracy')[['accuracy', 'total']])
    
    # Confusion matrix
    cm = results['confusion_matrix']
    labels = results['confusion_matrix_labels']

Save Your Basis
---------------

::

    # Save as HDF5 for reuse
    import anndata as ad
    
    # Create AnnData object from basis
    adata_basis = ad.AnnData(
        X=basis.T.values,
        obs=pd.DataFrame(index=basis.columns),
        var=pd.DataFrame(index=basis.index)
    )
    
    adata_basis.write_h5ad("my_custom_basis.h5ad")

Advanced: Cross-Validation
===========================

For more robust validation, use k-fold cross-validation::

    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=100,
        cv_folds=5,  # 5-fold cross-validation
        n_jobs=-1,   # Use all CPU cores
        random_state=42
    )
    
    # Cross-validation results
    cv_results = results['cv_results']
    avg_metrics = results['cv_avg_metrics']
    
    print(f"Average accuracy: {avg_metrics['accuracy_mean']:.3f} ± {avg_metrics['accuracy_std']:.3f}")

Advanced: Feature Selection with ANOVA
=======================================

For large datasets, you can select informative genes::

    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=100,
        do_anova=True,
        n_features=5000,  # Select top 5000 genes
        # OR use percentile:
        # anova_percentile=25  # Keep top 25% of genes
    )
    
    selected_genes = results['selected_genes']
    print(f"Selected {len(selected_genes)} informative genes")

Gene Contribution Analysis
===========================

Understand which genes drive cell type assignments::

    # Analyze contributions
    contributions = top.compute_gene_contributions(
        data=raw_counts,
        basis=basis,
        cell_types=['T cell', 'B cell', 'Macrophage']
    )
    
    # Find top genes for each cell type
    for cell_type, contrib in contributions.items():
        top_genes = top.find_top_contributing_genes(contrib, n_genes=20)
        print(f"\n{cell_type} - Top 20 genes:")
        print(top_genes)

Visualization
=============

Basic Plots
-----------

::

    import matplotlib.pyplot as plt
    
    # Plot top cell type matches
    top.plot_highest(projections['sample_1'], n=15)
    plt.tight_layout()
    plt.show()
    
    # Plot expression distribution of top genes
    top.plot_expression_distribution(processed_sample, n=10)
    plt.show()

2D Projections
--------------

::

    # Compare two cell types
    top.plot_two(
        projections, 
        celltype1='T cell',
        celltype2='B cell',
        alpha=0.5
    )
    plt.xlabel('T cell score')
    plt.ylabel('B cell score')
    plt.show()

Next Steps
==========

* Read the :doc:`api/index` for detailed function documentation
* See :doc:`tutorials` for more complex workflows
* Understand the :doc:`theory` behind the method
