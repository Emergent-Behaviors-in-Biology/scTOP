====================
Refrence Basis Management
====================

Functions for loading pre-made reference bases and creating custom reference bases from annotated data.

list_available_bases
====================

.. autofunction:: sctop.list_available_bases

**Purpose**: List all pre-computed reference bases available for download.

**Returns**:

* **list**: Names of available basis keys

**Example**::

    import sctop as top
    
    available = top.list_available_bases()
    print("Available bases:", available)
    # ['MCKO legacy']

**Notes**:

* Bases are hosted on Figshare
* New bases may be added over time
* Use the returned keys with ``load_basis()``

load_basis
==========

.. autofunction:: sctop.load_basis

**Purpose**: Download and load a pre-made reference basis from online storage.

**Key Features**:

* Automatic download with progress bar
* Local caching to avoid re-downloading
* Returns both basis and metadata

**Parameters**:

* **basis_key** (*str*): Name of the basis (from ``list_available_bases()``)
* **cache_dir** (*str or Path*, optional): Directory to cache downloaded files (default: system temp)
* **force_download** (*bool*, default=False): Re-download even if cached

**Returns**:

* **basis** (*pd.DataFrame*): Cell type basis (genes × cell types)
* **metadata** (*pd.DataFrame*): Cell type metadata

**Example**::

    # Load with default cache location
    basis, metadata = top.load_basis("MCKO legacy")
    
    # Use custom cache directory
    basis, metadata = top.load_basis(
        basis_key="MCKO legacy",
        cache_dir="./my_bases"
    )
    
    # Force re-download
    basis, metadata = top.load_basis(
        basis_key="MCKO legacy",
        force_download=True
    )
    
    print(f"Loaded basis: {basis.shape[0]} genes × {basis.shape[1]} cell types")
    print("Cell types:", list(basis.columns))

**Notes**:

* First download may take time depending on connection
* Subsequent loads are instant (cached)
* Basis is already processed (ready for ``score()``)
* Check metadata for cell type counts and other info

**Available Bases**:

MCKO legacy
-----------

* **Source**: Mouse Cell Atlas (Kotton lab variant)
* **Organism**: Mouse (*Mus musculus*)
* **Cell types**: Over 100 major cell types
* **Format**: Pre-processed, ready to use

create_basis
============

.. autofunction:: sctop.create_basis

**Purpose**: Create a custom cell type reference basis from annotated scRNA-seq data with automatic validation.

**Key Features**:

* Train-test split or k-fold cross-validation
* Optional ANOVA feature selection
* Parallel processing for speed
* Comprehensive performance metrics
* Confusion matrix and per-cell-type accuracy
* Memory-efficient chunking
* Automatic visualization of results

**Parameters**:

* **adata** (*ad.AnnData*): Annotated data object with cell type labels
* **cell_type_column** (*str*): Column name in ``adata.obs`` containing cell types
* **threshold** (*int*): Minimum number of cells required per cell type
* **test_size** (*float*, default=0.2): Fraction of data for testing (if cv_folds=None)
* **random_state** (*int*, default=42): Random seed for reproducibility
* **n_jobs** (*int*, default=-1): Number of parallel jobs (-1 = all cores)
* **do_anova** (*bool*, default=False): Whether to perform feature selection
* **n_features** (*int*, default=20000): Number of features to select (if do_anova=True)
* **anova_percentile** (*float*, optional): Keep top percentile of genes (overrides n_features)
* **spec_value** (*float*, default=0.1): Threshold for "unspecified" predictions
* **outer_chunks** (*int*, default=10): Number of chunks for parallel scoring
* **inner_chunk_size** (*int*, default=1000): Chunk size for internal processing
* **n_scoring_jobs** (*int*, default=4): Number of parallel jobs for scoring
* **cv_folds** (*int*, optional): If specified, use k-fold cross-validation instead of train-test
* **plot_results** (*bool*, default=True): Whether to generate performance plots

**Returns**:

Dictionary containing:

* **basis** (*pd.DataFrame*): Final cell type basis (genes × cell types)
* **selected_genes** (*np.ndarray* or None): Selected genes if ANOVA was used
* **training_IDs** (*np.ndarray*): Cell IDs used for training
* **test_IDs** (*np.ndarray*): Cell IDs used for testing
* **metrics** (*dict*): Performance metrics including:
  
  * ``accuracy``: Top-1 accuracy
  * ``top3_accuracy``: Top-3 accuracy
  * ``unspecified_rate``: Fraction of low-confidence predictions
  * ``f1_macro``, ``f1_weighted``: F1 scores
  * ``precision_macro``, ``precision_weighted``
  * ``recall_macro``, ``recall_weighted``

* **confusion_matrix** (*np.ndarray*): Confusion matrix
* **confusion_matrix_labels** (*list*): Cell type labels for confusion matrix
* **per_cell_type** (*pd.DataFrame*): Per-cell-type accuracy metrics
* **f1_scores** (*pd.DataFrame*): F1 scores for each cell type
* **true_labels** (*list*): True cell type labels (test set)
* **predicted_labels** (*list*): Predicted cell type labels (test set)
* **cv_results** (*list*, optional): Cross-validation fold results (if cv_folds specified)
* **cv_avg_metrics** (*dict*, optional): Average CV metrics (if cv_folds specified)

**Example 1: Basic Usage**::

    import anndata as ad
    import sctop as top
    
    # Load annotated data
    adata = ad.read_h5ad("mouse_atlas.h5ad")
    
    # Create basis with validation
    # plot_results = True will display performance metrics like accuracy and F1 scores
    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=100,  # At least 100 cells per type
        test_size=0.2,
        plot_results=True
    )
    
    # Get the basis
    basis = results['basis']

**Example 2: Memory-Efficient for Large Datasets**::

    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=100,
        n_jobs=4,  # Limit parallel jobs
        n_scoring_jobs=2,
        inner_chunk_size=500,  # Smaller chunks
        outer_chunks=50  # More chunks
    )


**Example 3: With Feature Selection**::
    IMPORTANT: The best practice is to carefully curate your basis by seriously considering what is a biologically meaningful cell type and combining similar cell types (e.g. although epithelial cells are very specialized, many stromal and immune cell types are functionally identical). Merging similar cell types, dropping cell types with very few cells, and dropping questionably-annotated cell types should ALWAYS be done first before considering feature selection. Feature selection should be a last resort to improve performance after careful curation of cell types.

    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=100,
        do_anova=True,
        n_features=5000,  # Select top 5000 genes
        random_state=42
    )
    
    selected_genes = results['selected_genes']
    print(f"Selected {len(selected_genes)} informative genes")

**Example 4: Cross-Validation**::

    results = top.create_basis(
        adata=adata,
        cell_type_column='cell_type',
        threshold=100,
        cv_folds=5,  # 5-fold cross-validation
        n_jobs=-1,
        random_state=42
    )
    
    # Check CV results
    cv_metrics = results['cv_avg_metrics']
    print(f"CV Accuracy: {cv_metrics['accuracy_mean']:.3f} ± {cv_metrics['accuracy_std']:.3f}")

**Notes**:

* **Threshold**: Cell types with fewer than ``threshold`` cells are excluded
* **Unspecified**: Cells with max score < ``spec_value`` are marked "unspecified"
* **Parallelization**: Uses thread pools for shared-memory parallel processing
* **Memory**: Adjust ``inner_chunk_size`` and ``n_scoring_jobs`` for large datasets
* **ANOVA**: Reduces features but may remove important genes
* **Cross-validation**: More robust but slower than single train-test split

**Performance Tips**:

For small datasets (<10k cells)::

    results = create_basis(
        adata, 'cell_type', 100,
        n_jobs=-1,  # Use all cores
        n_scoring_jobs=8
    )

For large datasets (>100k cells)::

    results = create_basis(
        adata, 'cell_type', 100,
        n_jobs=4,  # Conservative
        n_scoring_jobs=2,
        inner_chunk_size=500,
        outer_chunks=100,
    )

**Troubleshooting**:

* Check per-cell-type metrics to find problematic types
* Examine confusion matrix for similar types
* Consider merging similar cell types
* Increase ``threshold`` to require more cells per type
* Try ANOVA feature selection as a last resort

analyze_sample_contributions
=============================

.. autofunction:: sctop.analyze_sample_contributions

**Purpose**: Analyze which genes contribute most to cell type scores across multiple samples.

**Key Features**:

* Computes gene-level contributions for each cell type
* Identifies top contributing genes per sample
* Handles multiple samples/clusters simultaneously
* Optional data processing

**Parameters**:

* **sample_data_dict** (*dict*): Maps sample_name → expression_data (DataFrame or array)
* **basis** (*pd.DataFrame*): Cell type basis
* **cell_types** (*list*, optional): Cell types to analyze (default: all in basis)
* **n_top_genes** (*int*, default=20): Number of top genes to identify
* **process_data** (*bool*, default=True): Whether to process raw counts

**Returns**:

Dictionary with structure::

    {
        cell_type: {
            'contributions': {sample_name: contribution_matrix},
            'top_genes': {sample_name: [gene1, gene2, ...]},
            'expressions': {sample_name: expression_matrix}
        }
    }

**Example**::

    # Prepare sample dictionary
    sample_dict = {
        'cluster_1': cluster1_counts,
        'cluster_2': cluster2_counts,
        'cluster_3': cluster3_counts
    }
    
    # Analyze contributions
    results = top.analyze_sample_contributions(
        sample_data_dict=sample_dict,
        basis=basis,
        cell_types=['T cell', 'B cell', 'Macrophage'],
        n_top_genes=20
    )
    
    # Get top genes for T cells in cluster 1
    t_cell_genes = results['T cell']['top_genes']['cluster_1']
    print("Top T cell genes in cluster 1:", t_cell_genes[:10])
    
    # Get contribution matrix
    t_cell_contrib = results['T cell']['contributions']['cluster_1']

**Notes**:

* Contribution = gene expression × predictivity
* Higher contribution = gene more responsible for that cell type score
* Useful for identifying markers and understanding assignments
