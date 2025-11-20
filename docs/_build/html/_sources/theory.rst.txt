======
Theory
======

This page explains the theoretical foundation and mathematical details of scTOP.

Overview
========

scTOP (Single-cell Type Order Parameters) identifies cell phenotypes by projecting single-cell expression profiles onto a reference basis of known cell types. The method uses order parameters for Hopfield networks as defined by Kanter and Sompolinsky.

Theoretical Background
======================

The approach is based on the work by Lang et al. (2014) in `Epigenetic Landscapes Explain Partially Reprogrammed Cells and Identify Key Reprogramming Genes <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003734>`_.

Key Concepts
------------

**Order Parameters**: In physics, order parameters are a macroscopic variable describing the phase of a system. For cell types, Hopfield-inspired order parameters measure cellular identity.

**Cell Type as Attractor State**: This approach is inspired by the Waddington landscape and treating cell types as attractors in a Hopfield network. Each cell type corresponds to a stable attractor in gene expression space.

**Projection onto Basis**: By creating a basis from reference cell types, we can decompose any sample's expression profile into contributions from each cell type.

Mathematical Framework
======================

Notation
--------

* :math:`\\xi`: Basis matrix (genes × cell types)
* :math:`\\mathbf{s}`: Processed sample vector (genes × 1)
* :math:`\\mathbf{a}`: Cell type scores (cell types × 1)
* :math:`G`: Number of genes
* :math:`C`: Number of cell types

Data Processing
===============

Step 1: Normalization
----------------------

Each sample is independently normalized by library size:

.. math::

   x'_{g,i} = \\frac{x_{g,i}}{\\sum_{g'} x_{g',i}}

where :math:`x_{g,i}` is the raw count for gene :math:`g` in sample :math:`i`.

**Purpose**: Remove technical variation due to sequencing depth.

Step 2: Log Transformation
---------------------------

.. math::

   y_{g,i} = \\log_2(x'_{g,i} + 1)

**Purpose**: Stabilize variance and reduce the influence of highly expressed genes.

Step 3: Rank-Based Z-Score Transform
-------------------------------------

For each sample :math:`i`:

1. **Rank genes**: Assign ranks :math:`r_{g,i} \\in [1, G]` to genes based on :math:`y_{g,i}`

2. **Handle ties**: Use average rank for tied values

3. **Convert to percentiles**:

   .. math::

      p_{g,i} = \\frac{r_{g,i}}{G + 1}

4. **Apply inverse normal CDF** (probit transform):

   .. math::

      z_{g,i} = \\Phi^{-1}(p_{g,i}) = \\sqrt{2} \\cdot \\text{erf}^{-1}(2p_{g,i} - 1)

**Purpose**: Ranking is valuable because there is a lot of batch-to-bath variability in scRNA-seq data. Ranking converts all gene expression to relative ordering instead of absolute expression, helping ameliorate batch effects.Rank-based methods are more stable across experiments


Basis Construction
==================

Creating the Reference Basis
-----------------------------

Given annotated reference data, the basis is constructed by:

1. **Group cells by type**: Collect all cells :math:`\\{\\mathbf{s}_1, \\mathbf{s}_2, \\ldots, \\mathbf{s}_N\\}` for each cell type

2. **Process each cell**: Apply the full processing pipeline to cells individually, avoiding cross-sample operations

3. **Average within type**:

   .. math::

      \\mathbf{b}_{ct} = \\frac{1}{N_{ct}} \\sum_{i \\in ct} \\mathbf{s}_i

4. **Assemble basis**: :math:`\\xi = [\\mathbf{b}_1, \\mathbf{b}_2, \\ldots, \\mathbf{b}_C]`

**Result**: Each column of :math:`\\xi` represents the average processed expression profile for one cell type.

Basis Properties
----------------

**Non-orthogonality**: The gene expression profiles of cell types are correlated, so :math:`\\xi^T \\xi \\neq I`. The basis captures these relationships.

**Similarity Matrix**: :math:`A = \\xi^T \\xi / G` encodes cell type similarities:

* :math:`A_{ii} \\approx 1`: Self-similarity
* :math:`A_{ij}`: Similarity between types :math:`i` and :math:`j`

Scoring (Projection)
====================

Non-Orthogonal Projection
-------------------------

To score a sample :math:`\\mathbf{s}` against basis :math:`\\xi`, we solve:

.. math::

   \\mathbf{a}^* = (\\xi^T \\xi)^{-1} \\xi^T \\mathbf{s}
   
This gives the coefficients :math:`\\mathbf{a}^*` that best reconstruct :math:`\\mathbf{s}` as a linear combination of basis vectors.

Interpretation of Scores
-------------------------

* **Higher score** → better match to that cell type
* **Scores can be negative**: Sample anti-correlates with that type
* **Scores can be very low**: Because of dropout, scores for single cells may be low across all types. Try using pseudo-bulk (i.e. average across populations) for more robust scores.
* **Use endogenous samples as controls**: To interpret scores, compare to known endogenous samples processed the same way whenever possible. This gives a baseline for expected score ranges given dropout, technical variability, and data quality.

Predictivity Matrix
====================

The matrix :math:`\\eta` is called the **predictivity matrix**:

.. math::

   \\eta = (\\xi^T \\xi)^{-1} \\xi^T / G

Each entry :math:`\\eta_{ct,g}` represents:

**How much does expression of gene** :math:`g` **contribute to the score for cell type** :math:`ct`?

Gene Contributions
------------------

For a sample :math:`\\mathbf{s}` and cell type :math:`ct`:

.. math::

   a_{ct} = \\sum_g \\eta_{ct,g} \\cdot s_g

The **contribution of gene** :math:`g` to type :math:`ct` is:

.. math::

   c_{g,ct} = \\eta_{ct,g} \\cdot s_g

**Interpretation**:

* :math:`c_{g,ct} > 0`: High expression of gene :math:`g` increases score for type :math:`ct`. Low expression of gene :math:`g` decreases score.
* :math:`c_{g,ct} < 0`: High expression of gene :math:`g` decreases score for type :math:`ct`. Low expression of gene :math:`g` increases score.
* :math:`|c_{g,ct}|` large: Gene :math:`g` strongly influences assignment

**Use Case**: Identify which genes drive cell type assignments, validate with known markers.

Performance Metrics
===================

Top-1 Accuracy
--------------

.. math::

   \\text{Acc}_1 = \\frac{\\#\\{\\text{argmax}(\\mathbf{a}) = \\text{true type}\\}}{N}

The fraction of cells assigned to the correct type.

Top-3 Accuracy
--------------

.. math::

   \\text{Acc}_3 = \\frac{\\#\\{\\text{true type} \\in \\text{top 3 of } \\mathbf{a}\\}}{N}

The fraction where true type is in top 3 predictions.

F1 Score
--------

Harmonic mean of precision and recall:

.. math::

   F_1 = 2 \\cdot \\frac{\\text{Precision} \\cdot \\text{Recall}}{\\text{Precision} + \\text{Recall}}

Computed per cell type and aggregated (macro or weighted average).

Unspecified Rate
----------------

Fraction of predictions with low confidence:

.. math::

   \\text{Unspec} = \\frac{\\#\\{\\max(\\mathbf{a}) < \\tau\\}}{N}

where :math:`\\tau` is the specification threshold (default: 0.1).

Feature Selection
=================

ANOVA-Based Selection
---------------------

scTOP can select informative genes using one-way ANOVA:

1. **Compute F-statistic** for each gene:

   .. math::

      F_g = \\frac{\\text{MS}_{\\text{between}}}{\\text{MS}_{\\text{within}}}

   where MS = mean square.

2. **Select top genes**: Keep top :math:`k` genes or top percentile by F-score.

3. **Standardize basis**: Optionally standardize selected features to unit norm.

**Disclaimer**:
May discard important and relevant genes. May be misleading in defining what distinct cell types are.

Cross-Validation
================

To robustly evaluate basis quality, scTOP supports k-fold cross-validation:

1. **Split data** into :math:`k` folds
2. **For each fold**:

   * Train basis on :math:`k-1` folds
   * Test on held-out fold
   * Compute metrics

3. **Average metrics** across folds:

   .. math::

      \\mu_{\\text{metric}} = \\frac{1}{k} \\sum_{i=1}^k \\text{metric}_i

4. **Create final basis** on all data

**Advantage**: More reliable performance estimate than single train-test split.

Best Practices
==============

For Optimal Results
-------------------

1. **Use high-quality reference data**: The entire method relies on the idea that each cell type is an attractor basin. Although it can distinguish highly-correlated cell types very well, it is very important to carefully curate a reference basis. The quality of scTOP's output depends heavily on the quality of the reference basis.
* Use sufficient cells per type (at least 100, preferably much more)
* Ensure accurate and consistent annotations (make sure there are no duplicates or highly-similar types, e.g. "T cell CD4+" and "T cell")

2. **Common gene space**: Ensure query and reference share many genes

3. **Consistent processing**: Do not pre-normalize or log-transform; scTOP handles this

4. **Validate basis**: Check accuracy, confusion matrix, and F1 scores before use. It's good to iteratively create a basis by merging similar cell types and dropping suspicious ones.

For Interpretation
------------------

1. **Examine contributions**: Use gene contributions to validate and understand assignments

2. **Compare to markers**: Check if known markers have high contributions

3. **Visualize cell trajectories**: Plot scores in 2D or 3D to explore differentiation paths. scTOP is not just for annotating cells; it provides a set of cell fate coordinates to study the trajectories of differentiation. This is the most powerful use of scTOP. 

4. **Check confusion matrix**: Understand which types are commonly confused


References
==========

1. Kanter, I., & Sompolinsky, H. (1987). Associative recall of memory without errors. *Physical Review A*, 35(1), 380.

2. Yampolskaya, Maria, et al. "scTOP: physics-inspired order parameters for cellular identification and visualization." Development 150.21 (2023): dev201873.

3. Yampolskaya, Maria, and Pankaj Mehta. "Hopfield Networks as Models of Emergent Function in Biology." arXiv preprint arXiv:2506.13076 (2025).

4. Lang, A. H., Li, H., Collins, J. J., & Mehta, P. (2014). `Epigenetic Landscapes Explain Partially Reprogrammed Cells and Identify Key Reprogramming Genes <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003734>`_. *PLoS Computational Biology*, 10(8), e1003734.
