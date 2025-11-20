==================
Analysis Functions
==================

Utilities for computing metrics, gene contributions, and understanding basis performance.

compute_predictivity
====================

.. autofunction:: sctop.compute_predictivity

**Purpose**: Compute the predictivity matrix showing how each gene contributes to each cell type score.

**Mathematical Background**:

The predictivity matrix :math:`\\eta` is defined as:

.. math::

   \\eta = (\\xi^T \\xi)^{-1} \\xi^T

where :math:`\\xi` is the basis matrix (genes × cell types).

Each entry :math:`\\eta_{ct,g}` represents how much gene :math:`g` contributes to the score for cell type :math:`ct`.

**Parameters**:

* **basis** (*pd.DataFrame*): Cell type basis (genes × cell types)

**Returns**:

* **pd.DataFrame**: Predictivity matrix (cell types × genes)

**Example**::

    # Compute predictivity from basis
    predictivity = top.compute_predictivity(basis)
    
    # See which genes contribute most to T cell score
    t_cell_pred = predictivity.loc['T cell'].sort_values(ascending=False)
    print("Top 10 T cell predictor genes:")
    print(t_cell_pred.head(10))
    
    # Positive values: gene expression increases cell type score
    # Negative values: gene expression decreases cell type score

**Notes**:

* Predictivity is computed once and can be reused
* Used internally by contribution analysis functions
* Reveals gene importance for each cell type

compute_gene_contributions
===========================

.. autofunction:: sctop.compute_gene_contributions

**Purpose**: Compute gene-level contributions to cell type scores for specific samples.

**How It Works**:

For each cell type, contribution is computed as:

.. math::

   \\text{contribution}_{g,s} = \\text{expression}_{g,s} \\times \\text{predictivity}_{ct,g}

where:

* :math:`g` = gene
* :math:`s` = sample
* :math:`ct` = cell type

**Parameters**:

* **data** (*pd.DataFrame or array*): Expression data (genes × samples)
* **basis** (*pd.DataFrame*): Cell type basis
* **predictivity** (*pd.DataFrame*, optional): Precomputed predictivity (computed if None)
* **cell_types** (*list*, optional): Cell types to analyze (default: all)
* **process_data** (*bool*, default=True): Whether to process raw counts first

**Returns**:

* **dict**: Maps cell_type → contribution_matrix (genes × samples)

**Example**::

    # Compute contributions for all cell types
    contributions = top.compute_gene_contributions(
        data=sample_data,
        basis=basis,
        process_data=True
    )
    
    # Analyze T cell contributions
    t_cell_contrib = contributions['T cell']
    
    # Find genes driving T cell score in sample 1
    sample1_contrib = t_cell_contrib['sample_1'].sort_values(ascending=False)
    print("Top 20 genes driving T cell assignment:")
    print(sample1_contrib.head(20))

**Use Cases**:

1. **Identify marker genes**: Which genes drive cell type assignments?
2. **Validate assignments**: Do expected markers contribute highly?
3. **Compare samples**: How do contributions differ across conditions?
4. **Quality control**: Are unexpected genes contributing?

find_top_contributing_genes
============================

.. autofunction:: sctop.find_top_contributing_genes

**Purpose**: Identify the top contributing genes from a contribution matrix.

**Parameters**:

* **contributions** (*pd.DataFrame*): Gene contributions (genes × samples)
* **n_genes** (*int*, default=20): Number of top genes to return
* **aggregate** (*str*, default='mean'): How to aggregate across samples

  * ``'mean'``: Average contribution
  * ``'median'``: Median contribution
  * ``'max'``: Maximum contribution

**Returns**:

* **pd.Series**: Top genes with their aggregated contribution scores

**Example**::

    # Get contributions for a cell type
    contrib = contributions['Macrophage']
    
    # Find top 30 genes by mean contribution
    top_genes = top.find_top_contributing_genes(
        contributions=contrib,
        n_genes=30,
        aggregate='mean'
    )
    
    print("Top macrophage marker genes:")
    for gene, score in top_genes.items():
        print(f"  {gene}: {score:.4f}")

**Aggregation Methods**:

* **mean**: Good for consistent markers across samples
* **median**: Robust to outliers
* **max**: Finds genes with strongest contribution in any sample

perform_anova_selection
=======================

.. autofunction:: sctop.perform_anova_selection

**Purpose**: Select informative genes using one-way ANOVA F-test.

**Description**:

IMPORTANT: The best practice is to carefully curate your basis by seriously considering what is a biologically meaningful cell type and combining similar cell types (e.g. although epithelial cells are very specialized, many stromal and immune cell types are functionally identical). Merging similar cell types, dropping cell types with very few cells, and dropping questionably-annotated cell types should ALWAYS be done first before considering feature selection. Feature selection should be a last resort to improve performance after careful curation of cell types.

ANOVA feature selection identifies genes that differ significantly across cell types. This selects genes that discriminate between cell types the most. These may or may not be biologically meaningful.

**Parameters**:

* **basis** (*pd.DataFrame*): The basis matrix
* **adata** (*ad.AnnData*): Full annotated dataset
* **training_IDs** (*np.ndarray*): Training cell IDs
* **cell_type_column** (*str*): Column with cell type labels
* **n_features** (*int*, default=2000): Number of genes to keep
* **percentile** (*float*, optional): Keep top percentile (overrides n_features)
* **standardize** (*bool*, default=True): Whether to standardize basis after selection

**Returns**:

* **basis_selected** (*pd.DataFrame*): Basis with only selected genes
* **selected_genes** (*np.ndarray*): Array of selected gene names

**Example**::

    from sctop.utils import perform_anova_selection
    
    # Select top 5000 genes
    basis_filtered, genes = perform_anova_selection(
        basis=basis,
        adata=adata,
        training_IDs=train_ids,
        cell_type_column='cell_type',
        n_features=5000
    )
    
    print(f"Reduced from {basis.shape[0]} to {basis_filtered.shape[0]} genes")

**Notes**:

* Higher F-scores indicate better discrimination
* May remove cell-type-specific markers with moderate expression
* Most useful for large gene sets (>20k genes)
* Standardization ensures basis vectors have unit norm

calculate_metrics
=================

.. autofunction:: sctop.calculate_metrics

**Purpose**: Calculate comprehensive classification metrics from predictions.

**Parameters**:

* **true_labels** (*list*): Ground truth cell type labels
* **predicted_labels** (*list*): Predicted cell type labels
* **total_cells** (*int*): Total number of cells
* **accuracies** (*dict*): Dictionary with counts for top1, top3, unspecified

**Returns**:

* **dict**: Metrics including:

  * ``accuracy``: Top-1 accuracy (correct predictions / total)
  * ``top3_accuracy``: Top-3 accuracy
  * ``unspecified_rate``: Fraction below confidence threshold
  * ``f1_macro``: Macro-averaged F1 score
  * ``f1_weighted``: Weighted F1 score
  * ``precision_macro``, ``precision_weighted``
  * ``recall_macro``, ``recall_weighted``
  * ``total_cells``

**Example**::

    from sctop.utils import calculate_metrics
    
    metrics = calculate_metrics(
        true_labels=true_types,
        predicted_labels=pred_types,
        total_cells=len(test_ids),
        accuracies={'top1': 850, 'top3': 920, 'unspecified': 30}
    )
    
    print(f"Accuracy: {metrics['accuracy']:.3f}")
    print(f"F1 (weighted): {metrics['f1_weighted']:.3f}")

calculate_per_cell_type_accuracy
=================================

.. autofunction:: sctop.calculate_per_cell_type_accuracy

**Purpose**: Compute accuracy metrics for each individual cell type.

**Parameters**:

* **cell_accuracies** (*dict*): Per-cell accuracy information

**Returns**:

* **pd.DataFrame**: Per-cell-type metrics with columns:

  * ``correct``: Number of correctly classified cells
  * ``total``: Total cells of this type
  * ``accuracy``: Fraction correct
  * ``top3_correct``: Correct within top 3 predictions
  * ``top3_accuracy``: Top-3 accuracy
  * ``unspecified_count``: Cells below confidence threshold
  * ``unspecified_rate``: Fraction unspecified

**Example**::

    from sctop.utils import calculate_per_cell_type_accuracy
    
    per_type = calculate_per_cell_type_accuracy(cell_accs)
    
    # Find worst performers
    worst = per_type.nsmallest(10, 'accuracy')
    print("\nWorst performing cell types:")
    print(worst[['accuracy', 'total']])
    
    # Find types with many unspecified
    high_unspec = per_type.nsmallest(10, 'unspecified_rate')

**Notes**:

* Sorted by accuracy (best to worst)
* Useful for identifying problematic cell types
* Consider merging types with low accuracy and high confusion

run_scoring_parallel
====================

.. autofunction:: sctop.run_scoring_parallel

**Purpose**: Score test cells against basis in parallel (internal function used by ``create_basis``).

**Key Features**:

* Thread-based parallelism for shared-memory efficiency
* Automatic chunking of test set
* Progress bar via tqdm
* Returns detailed per-cell metrics

**Parameters**:

* **adata** (*ad.AnnData*): Full dataset
* **basis** (*pd.DataFrame*): Cell type basis
* **test_IDs** (*np.ndarray*): Test cell IDs to score
* **cell_type_column** (*str*): Cell type column name
* **spec_value** (*float*): Threshold for unspecified predictions
* **outer_chunks** (*int*): Number of chunks
* **inner_chunk_size** (*int*): Chunk size for internal processing
* **n_jobs** (*int*, default=4): Number of parallel workers

**Returns**:

Tuple of:

* **cell_accuracies** (*dict*): Per-cell results
* **true_labels** (*list*): Ground truth labels
* **predicted_labels** (*list*): Predicted labels
* **accuracies** (*dict*): Aggregate counts

**Example**::

    from sctop.utils import run_scoring_parallel
    
    cell_accs, true_labs, pred_labs, accs = run_scoring_parallel(
        adata=adata,
        basis=basis,
        test_IDs=test_ids,
        cell_type_column='cell_type',
        spec_value=0.1,
        outer_chunks=20,
        inner_chunk_size=500,
        n_jobs=4
    )
    
    print(f"Top-1 accuracy: {accs['top1'] / len(test_ids):.3f}")

**Performance Tuning**:

For fast scoring::

    n_jobs=8, outer_chunks=50, inner_chunk_size=1000

For memory-constrained::

    n_jobs=2, outer_chunks=100, inner_chunk_size=500

plot_performance_summary
========================

.. autofunction:: sctop.plot_performance_summary

**Purpose**: Generate comprehensive visualization of classification performance.

**Creates**:

1. **Confusion Matrix**: Normalized by true labels (recall)
2. **F1 Score Bar Plot**: Per-cell-type F1 scores

**Parameters**:

* **true_labels** (*list*): True cell type labels
* **predicted_labels** (*list*): Predicted labels
* **figsize_base** (*int*, default=10): Base figure size (scales with #types)
* **f1_df** (*pd.DataFrame*, optional): Precomputed F1 scores

**Example**::

    from sctop.utils import plot_performance_summary
    
    plot_performance_summary(
        true_labels=results['true_labels'],
        predicted_labels=results['predicted_labels'],
        f1_df=results['f1_scores']
    )

**Notes**:

* Automatically called by ``create_basis`` if ``plot_results=True``
* Figure size scales with number of cell types
* Confusion matrix is normalized (shows recall per type)
* Useful for identifying confused cell type pairs

print_metrics
=============

.. autofunction:: sctop.print_metrics

**Purpose**: Pretty-print metrics dictionary.

**Parameters**:

* **metrics** (*dict*): Metrics dictionary from ``calculate_metrics``

**Example**::

    from sctop.utils import print_metrics
    
    print_metrics(results['metrics'])

Output::

    Accuracy (Top-1): 0.8723
    Top-3 Accuracy: 0.9541
    Unspecified Rate: 0.0234
    F1 Score (Macro): 0.8456
    F1 Score (Weighted): 0.8701
    Precision (Macro): 0.8532
    Precision (Weighted): 0.8745
    Recall (Macro): 0.8512
    Recall (Weighted): 0.8723
