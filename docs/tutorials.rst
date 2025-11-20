=========
Tutorials
=========

This page contains detailed tutorials for common scTOP workflows.

Tutorial 1: Basic Cell Type Identification
===========================================

**Goal**: Score samples against a reference basis to identify cell types.

Step 1: Load Reference Basis
-----------------------------

::

    import sctop as top
    import pandas as pd
    
    # List available references
    print(top.list_available_bases())
    
    # Load Mouse Cell Atlas basis
    basis, metadata = top.load_basis("MCKO legacy")
    
    print(f"Basis: {basis.shape[0]} genes × {basis.shape[1]} cell types")
    print(f"Cell types: {list(basis.columns)[:10]}...")  # First 10

Step 2: Load Your Data
-----------------------

::

    import anndata as ad
    
    # Load your scRNA-seq data
    # Can be from: Scanpy, Cell Ranger, CSV, etc.
    adata = ad.read_h5ad("your_sample.h5ad")
    
    # OR load from CSV
    # counts = pd.read_csv("counts.csv", index_col=0)  # genes × cells
    
    print(f"Data: {adata.shape[0]} cells × {adata.shape[1]} genes")

Step 3: Process and Score
--------------------------

::

    # Convert to DataFrame if using AnnData
    if 'adata' in locals():
        counts = pd.DataFrame(
            adata.X.T if hasattr(adata.X, 'toarray') else adata.X.T,
            index=adata.var_names,
            columns=adata.obs_names
        )
    
    # Process: normalize → log → rank → z-score
    processed = top.process(counts)
    
    # Score against basis
    scores = top.score(basis, processed)
    
    print(f"Scores: {scores.shape}")  # cell_types × cells

Step 4: Interpret Results
--------------------------

::

    import matplotlib.pyplot as plt
    
    # Single cell example
    cell_id = scores.columns[0]
    cell_scores = scores[cell_id].sort_values(ascending=False)
    
    print(f"\n{cell_id} - Top 10 matches:")
    print(cell_scores.head(10))
    
    # Visualize
    top.plot_highest(cell_scores, n=15)
    plt.title(f'Cell Type Scores: {cell_id}')
    plt.tight_layout()
    plt.show()

Step 5: Visualize trajectories 
--------------------------------

::
    import matplotlib.pyplot as plt

    top.plot_two(
        scores,
        cell_x='Alveolar Type 1',
        cell_y='Alveolar Type 2',
        gene='Nkx2-1'
    )


Tutorial 2: Creating a Custom Reference
========================================

**Goal**: Build your own reference basis from annotated data.

Step 1: Prepare Annotated Data
-------------------------------

::

    import anndata as ad
    import sctop as top
    
    # Load your atlas with cell type labels
    adata = ad.read_h5ad("annotated_atlas.h5ad")
    
    # Verify annotations
    print("Cell type column:", 'cell_type' in adata.obs.columns)
    print("\nCell type counts:")
    print(adata.obs['cell_type'].value_counts().head(20))

Step 2: Quality Check
----------------------

::

    import matplotlib.pyplot as plt
    
    # Check cell counts per type
    counts = adata.obs['cell_type'].value_counts()
    
    plt.figure(figsize=(12, 6))
    counts.head(30).plot(kind='barh')
    plt.xlabel('Number of Cells')
    plt.title('Cells per Type (Top 30)')
    plt.tight_layout()
    plt.show()
    
    # Decide on threshold (e.g., 100 cells minimum)
    threshold = 100
    types_above = counts[counts >= threshold]
    print(f"\n{len(types_above)} types with ≥{threshold} cells")

Step 3: Create Basis
---------------------

::

    # Create basis with validation
    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=100,
        test_size=0.2,
        random_state=42,
        n_jobs=-1,
        plot_results=True
    )
    
    # Extract results
    basis = results['basis']

Step 4: Review Performance
---------------------------

::

    # Per-type accuracy
    per_type = results['per_cell_type']
    
    print("\nBest performing (Top 10):")
    print(per_type.head(10)[['accuracy', 'total']])
    
    print("\nWorst performing (Bottom 10):")
    print(per_type.tail(10)[['accuracy', 'total']])
    
    # Check confusion
    cm = results['confusion_matrix']
    labels = results['confusion_matrix_labels']
    
    # Find most confused pairs
    import numpy as np
    np.fill_diagonal(cm, 0)
    max_confusion_idx = np.unravel_index(cm.argmax(), cm.shape)
    print(f"\nMost confused pair:")
    print(f"  {labels[max_confusion_idx[0]]} → {labels[max_confusion_idx[1]]}")
    print(f"  {cm[max_confusion_idx]} misclassifications")

Step 5: Improve Basis (Optional)
---------------------------------

::

    # If performance is poor, try:
    
    # Option 1: Merge similar types (especially recommended for immune and stromal cells, or cell types that only differ over/under-expression of 1 gene)
    # (Do this to your adata before create_basis)
    mapping = {
        'T cell CD4+': 'T cell',
        'T cell CD8+': 'T cell',
        'T cell regulatory': 'T cell'
    }
    adata.obs['cell_type'] = adata.obs['cell_type'].replace(mapping)

    # Option 2: Drop questionably-annotated types
    # (Do this to your adata before create_basis)
    to_drop = ['Unknown', 'Doublet', 'Mesenchymal Stem Cell']
    adata = adata[~adata.obs['cell_type'].isin(to_drop), :]

    # Option 3: Increase threshold for number of cells required per cell type
    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=200,  # More stringent
        test_size=0.2,
        random_state=42
    )

    # Option 4: Use ANOVA feature selection (last resort)
    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=100,
        do_anova=True,
        n_features=5000,
        test_size=0.2,
    )

Step 6: Save Basis
-------------------

::

    # Save as HDF5
    basis_adata = ad.AnnData(
        X=basis.T.values,
        obs=pd.DataFrame(index=basis.columns),
        var=pd.DataFrame(index=basis.index)
    )
    basis_adata.write_h5ad("my_custom_basis.h5ad")
    
    # Save as CSV (for sharing)
    basis.to_csv("my_basis.csv")
    
    print(f"Saved basis: {basis.shape}")

Tutorial 3: Gene Contribution Analysis
=======================================

**Goal**: Understand which genes drive cell type assignments.

Step 1: Compute Contributions
------------------------------

::

    import sctop as top
    
    # From Tutorial 1, we have:
    # - basis: cell type reference
    # - counts: sample expression data
    
    # Compute contributions for specific types
    contributions = top.compute_gene_contributions(
        data=counts,
        basis=basis,
        cell_types=['T cell', 'B cell', 'Macrophage'],
        process_data=True
    )
    
    print("Computed contributions for:")
    for ct in contributions:
        print(f"  {ct}: {contributions[ct].shape}")

Step 2: Find Top Genes
-----------------------

::

    # For each cell type
    for cell_type in ['Alveolar Type 1 Cell', 'Alveolar Type 2 Cell', 'Basal Cell']:
        contrib = contributions[cell_type]
        
        top_genes = top.find_top_contributing_genes(
            contrib, 
            n_genes=20,
            aggregate='mean'
        )
        
        print(f"\n{cell_type} - Top 20 marker genes:")
        for gene, score in top_genes.items():
            print(f"  {gene}: {score:.4f}")

Step 3: Visualize Contributions
--------------------------------

::

    # Single cell type heatmap
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    contrib = contributions['T cell']
    top_genes = top.find_top_contributing_genes(contrib, n_genes=30).index
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        contrib.loc[top_genes],
        cmap='YlOrRd',
        cbar_kws={'label': 'Contribution'}
    )
    plt.title('Top 30 T Cell Gene Contributions', fontsize=16)
    plt.ylabel('Genes', fontsize=14)
    plt.xlabel('Cells', fontsize=14)
    plt.tight_layout()
    plt.show()


Tutorial 4: Cross-Validation Workflow
======================================

**Goal**: Robustly evaluate basis quality using cross-validation.

Step 1: 5-Fold Cross-Validation
--------------------------------

::

    import sctop as top
    
    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=100,
        cv_folds=5,
        random_state=42,
        n_jobs=-1,
        plot_results=False  # We'll plot manually
    )
    
    # Cross-validation summary
    cv_metrics = results['cv_avg_metrics']
    
    print("Cross-Validation Results:")
    print(f"  Accuracy: {cv_metrics['accuracy_mean']:.3f} ± {cv_metrics['accuracy_std']:.3f}")
    print(f"  F1 (macro): {cv_metrics['f1_macro_mean']:.3f} ± {cv_metrics['f1_macro_std']:.3f}")

Step 2: Plot Fold-by-Fold Performance
--------------------------------------

::

    import matplotlib.pyplot as plt
    import numpy as np
    
    cv_results = results['cv_results']
    
    # Extract metrics per fold
    folds = [r['fold'] for r in cv_results]
    accuracies = [r['metrics']['accuracy'] for r in cv_results]
    f1_scores = [r['metrics']['f1_weighted'] for r in cv_results]
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    ax1.bar(folds, accuracies, color='steelblue', alpha=0.7)
    ax1.axhline(np.mean(accuracies), color='red', linestyle='--', label='Mean')
    ax1.set_xlabel('Fold')
    ax1.set_ylabel('Accuracy')
    ax1.set_title('Accuracy per Fold')
    ax1.legend()
    ax1.set_ylim(0, 1)
    
    ax2.bar(folds, f1_scores, color='coral', alpha=0.7)
    ax2.axhline(np.mean(f1_scores), color='red', linestyle='--', label='Mean')
    ax2.set_xlabel('Fold')
    ax2.set_ylabel('F1 Score (Weighted)')
    ax2.set_title('F1 Score per Fold')
    ax2.legend()
    ax2.set_ylim(0, 1)
    
    plt.tight_layout()
    plt.show()

Step 3: Overall Confusion Matrix
---------------------------------

::

    from sklearn.metrics import ConfusionMatrixDisplay
    import matplotlib.pyplot as plt
    
    true_labels = results['all_true_labels']
    pred_labels = results['all_predicted_labels']
    
    # Plot confusion matrix (subset for readability)
    unique_labels = sorted(set(true_labels))[:20]  # Top 20 types
    
    fig, ax = plt.subplots(figsize=(14, 12))
    ConfusionMatrixDisplay.from_predictions(
        true_labels,
        pred_labels,
        labels=unique_labels,
        normalize='true',
        xticks_rotation='vertical',
        ax=ax,
        values_format='.2f'
    )
    plt.title('Confusion Matrix (Top 20 Cell Types)', fontsize=16)
    plt.tight_layout()
    plt.show()

Best Practices
==============

Data Quality
------------
* **High-Quality Annotations**: Ensure cell type labels are accurate and consistent.
* **Sufficient Cell Counts**: Aim for at least 100 cells per type; more is better.
* **Use raw counts**: scTOP expects raw counts, not normalized

Basis Quality
-------------

Check your basis quality by considering Top-1 accuracy and F1 scores. Seeing which cell types are confused for others is very informative, so make good use of the confusion matrix.

If the accuracy is bad:

* Increase cell count threshold
* Merge similar cell types
* Remove rare or ambiguous types
* Consider ANOVA feature selection as a last resort

Memory Management
-----------------

For large datasets (>100k cells, >20k genes)::

    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=100,
        n_jobs=4,  # Don't use all cores
        n_scoring_jobs=2,
        inner_chunk_size=500,  # Smaller chunks
        outer_chunks=100,  # More chunks
        do_anova=True,  # Reduce gene space
        n_features=5000
    )

Processing speed vs memory tradeoff::

    # Fast but memory-intensive
    n_jobs=-1, n_scoring_jobs=8, inner_chunk_size=2000
    
    # Slow but memory-efficient
    n_jobs=2, n_scoring_jobs=1, inner_chunk_size=500
