==============
Core Functions
==============

These are the fundamental data processing and scoring functions.

process
=======

.. autofunction:: sctop.process

**Purpose**: Transform raw count data into normalized, rank-based z-scores suitable for projection.

**Key Features**:

* Column-wise normalization (library size)
* Log transformation
* Rank-based conversion to normal scores
* Handles both individual samples and averaged profiles
* Optional chunking for memory efficiency

**Pipeline**:

1. Normalize each sample: ``data / sum(data)``
2. Log transform: ``log2(normalized + 1)``
3. Rank genes within each sample
4. Convert ranks to z-scores via inverse normal CDF

**Parameters**:

* **df_in** (*pd.DataFrame, np.ndarray, or sparse matrix*): Input data with genes as rows, samples as columns
* **average** (*bool*, default=False): If True, averages across samples (for creating cell type profiles)
* **chunk_size** (*int*, optional): Process data in chunks to reduce memory usage

**Returns**:

* **pd.DataFrame**: Processed data with same shape as input (or single column if average=True)

**Example**::

    import pandas as pd
    import sctop as top
    
    # Raw count matrix: genes × cells
    raw_counts = pd.DataFrame(...)
    
    # Process individual cells
    processed = top.process(raw_counts)
    
    # Or average across samples (for basis creation)
    averaged = top.process(raw_counts, average=True)
    
    # For large datasets, use chunking
    processed_chunked = top.process(raw_counts, chunk_size=500)

**Notes**:

* Input should be raw counts (not normalized or log-transformed)
* The rank-based transformation is robust to technical variation and outliers
* Sparse matrices are converted to dense arrays internally

score
=====

.. autofunction:: sctop.score

**Purpose**: Project processed samples onto a cell type basis to compute similarity scores.

**Key Features**:

* Finds the projections onto cell types
* Optional chunking for large sample sets
* Returns scores for all cell types simultaneously
* Can return full decomposition (projection matrix, correlation matrix, etc.)

**Parameters**:

* **basis** (*pd.DataFrame*): Cell type basis (genes × cell types)
* **sample** (*pd.DataFrame*): Processed sample data (genes × samples)
* **full_output** (*bool*, default=False): If True, returns [projections, correlation_matrix, eta]
* **chunk_size** (*int*, optional): Process samples in chunks

**Returns**:

* **pd.DataFrame**: Cell type scores (cell types × samples), or
* **list**: [projections, A, eta] if full_output=True, where:
  
  * **projections**: Cell type scores
  * **A**: Correlation matrix (:math:`B^T B / G`)
  * **eta**: Projection matrix

**Example**::

    # Score samples against basis
    scores = top.score(basis, processed_sample)
    
    # scores.loc['Basophil', 'sample_1'] gives Basophil score for sample_1
    
    # For large sample sets, use chunking
    scores = top.score(basis, processed_sample, chunk_size=1000)
    
    # Get full output for analysis
    projections, A, eta = top.score(basis, processed_sample, full_output=True)

**Notes**:

* Only genes present in both basis and sample are used
* Basis and sample must be **processed** (not raw counts)
* Typical workflow: ``process()`` → ``score()``

rank_zscore_fast
================

.. autofunction:: sctop.rank_zscore_fast

**Purpose**: Convert values to rank-based z-scores (internal function).

**Description**:

This function implements the rank-based normal scores transformation:

1. Rank values within each column
2. Handle ties by assigning average rank
3. Convert ranks to percentiles
4. Apply inverse normal CDF (probit transform)

**Parameters**:

* **arr** (*np.ndarray*): Input array (genes × samples)

**Returns**:

* **np.ndarray**: Z-scored array with same shape

**Example**::

    import numpy as np
    from sctop.processing import rank_zscore_fast
    
    # Raw expression values
    data = np.random.lognormal(size=(1000, 100))
    
    # Convert to z-scores
    z_scores = rank_zscore_fast(data)
    
    # Z-scores are approximately N(0,1) distributed
    print(f"Mean: {z_scores.mean():.3f}, Std: {z_scores.std():.3f}")

**Notes**:

* Uses Numba JIT compilation for speed
* Handles ties correctly (average rank method)
* Parallel processing across samples
* Used internally by ``process()``

Technical Details
=================

Data Processing Pipeline
------------------------

The complete processing pipeline in ``process()`` is:

.. code-block:: python

    # 1. Normalize by library size
    normalized = data / data.sum(axis=0)
    
    # 2. Log transform
    log_data = np.log2(normalized + 1)
    
    # 3. Rank-based z-score
    z_scored = rank_zscore_fast(log_data)

This transformation:

* Removes library size effects
* Stabilizes variance
* Produces normal-like distributions

Projection Mathematics
----------------------

The scoring step solves:

.. math::

   \\min_{\\mathbf{a}} \\|B \\mathbf{a} - \\mathbf{s}\\|^2

where:

* :math:`B`: basis matrix (genes × cell types)
* :math:`\\mathbf{s}`: processed sample vector
* :math:`\\mathbf{a}`: cell type scores (output)

The solution is:

.. math::

   \\mathbf{a} = (B^T B)^{-1} B^T \\mathbf{s} / G

where :math:`G` is the number of genes.

Memory Optimization
-------------------

Both ``process()`` and ``score()`` support chunking:

**process()** chunks:

* Splits samples into chunks
* Processes each chunk independently
* Concatenates results
* Reduces peak memory by factor of ``n_samples / chunk_size``

**score()** chunks:

* Computes projection matrix once: :math:`\\eta = (B^T B)^{-1} B^T / G`
* Applies to sample chunks: :math:`\\mathbf{a}_i = \\eta \\cdot \\mathbf{s}_i`
* Minimal memory overhead

Example: memory-efficient workflow::

    # For 100k cells × 20k genes
    processed = top.process(
        large_dataset, 
        chunk_size=500  # Process 500 cells at a time
    )
    
    scores = top.score(
        basis, 
        processed, 
        chunk_size=1000  # Score 1000 cells at a time
    )
