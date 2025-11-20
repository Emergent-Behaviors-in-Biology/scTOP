==============
Visualization
==============

Functions for visualizing cell type scores, gene expression, and analysis results.

plot_highest
============

.. autofunction:: sctop.plot_highest

**Purpose**: Plot horizontal bar chart of top N cell type scores for a sample.

**Parameters**:

* **projections** (*pd.Series or pd.DataFrame*): Cell type scores (if DataFrame, uses first column)
* **n** (*int*, default=10): Number of top cell types to show
* **ax** (*matplotlib.axes.Axes*, optional): Axes to plot on (creates new if None)
* **color** (*str*, default='olive'): Bar color
* **fontsize** (*int*, default=40): Font size for labels
* ****kwargs**: Additional arguments passed to ``plot.barh()``

**Returns**:

* **matplotlib.axes.Axes**: The axes object

**Example**::

    import matplotlib.pyplot as plt
    import sctop as top
    
    # Score a sample
    scores = top.score(basis, processed_sample)
    
    # Plot top 15 matches
    fig, ax = plt.subplots(figsize=(10, 8))
    top.plot_highest(
        scores['sample_1'], 
        n=15, 
        ax=ax,
        color='steelblue',
        fontsize=20
    )
    plt.title('Top Cell Type Matches', fontsize=24)
    plt.tight_layout()
    plt.show()

**Notes**:

* X-axis fixed to [0, 1] for consistency across plots
* Scores > 1 indicate an issue with normalization
* Useful for quick sample interpretation

plot_expression_distribution
=============================

.. autofunction:: sctop.plot_expression_distribution

**Purpose**: Plot boxplots of gene expression across samples for top expressed genes.

**Parameters**:

* **scores** (*pd.DataFrame*): Expression data (genes × samples)
* **n** (*int*, default=10): Number of top genes to show
* **ax** (*matplotlib.axes.Axes*, optional): Axes to plot on
* **box_color** (*str*, default='skyblue'): Box color
* **fontsize** (*int*, default=30): Font size
* ****kwargs**: Additional arguments passed to ``sns.boxplot()``

**Returns**:

* **matplotlib.axes.Axes**: The axes object

**Example**::

    import matplotlib.pyplot as plt
    import sctop as top
    
    # Show expression of top 20 genes
    fig, ax = plt.subplots(figsize=(14, 6))
    top.plot_expression_distribution(
        processed_sample,
        n=20,
        ax=ax,
        box_color='coral',
        fontsize=16
    )
    plt.title('Top Expressed Genes', fontsize=20)
    plt.tight_layout()
    plt.show()

**Notes**:

* Genes ordered by median expression
* Shows distribution across all samples
* Useful for QC and identifying highly variable genes

plot_two
========

.. autofunction:: sctop.plot_two

**Purpose**: Create 2D scatter plot comparing scores for two cell types.

**Parameters**:

* **projections** (*pd.DataFrame*): Cell type scores (cell types × samples)
* **celltype1** (*str*): First cell type (x-axis)
* **celltype2** (*str*): Second cell type (y-axis)
* **gene** (*str*, optional): Gene name to color points by
* **gene_expressions** (*pd.DataFrame*, optional): Gene expression data (genes × samples)
* **ax** (*matplotlib.axes.Axes*, optional): Axes to plot on
* ****kwargs**: Additional arguments passed to ``sns.scatterplot()``

**Returns**:

* None (modifies axes in place)

**Example 1: Basic 2D plot**::

    import matplotlib.pyplot as plt
    import sctop as top
    
    # Compare T cells vs B cells
    fig, ax = plt.subplots(figsize=(8, 8))
    top.plot_two(
        projections,
        celltype1='T cell',
        celltype2='B cell',
        ax=ax,
        alpha=0.5,
        s=20
    )
    plt.xlabel('T cell score', fontsize=16)
    plt.ylabel('B cell score', fontsize=16)
    plt.title('T cell vs B cell Scores', fontsize=18)
    plt.show()

**Example 2: Color by gene expression**::

    # Color points by CD3D expression
    fig, ax = plt.subplots(figsize=(10, 8))
    top.plot_two(
        projections,
        celltype1='T cell',
        celltype2='B cell',
        gene='CD3D',
        gene_expressions=processed_sample,
        ax=ax,
        s=30
    )
    plt.xlabel('T cell score', fontsize=16)
    plt.ylabel('B cell score', fontsize=16)
    plt.title('T vs B cells (colored by CD3D)', fontsize=18)
    plt.show()

**Use Cases**:

* **Compare related cell types**: See separation between similar types
* **Validate markers**: Check if marker genes correlate with expected type
* **Identify doublets**: Cells with high scores for multiple incompatible types
* **Explore transitions**: Identify intermediate states

plot_all_contributions
=======================

.. autofunction:: sctop.plot_all_contributions

**Purpose**: Generate and save gene contribution scatter plots for all cell types and samples.

**Key Features**:

* Automatically creates organized directory structure
* Plots expression vs contribution for each gene
* Highlights top contributing genes
* Optional highlighting of specific genes (e.g., known markers)
* Saves high-resolution images

**Parameters**:

* **results** (*dict*): Results from ``analyze_sample_contributions``
* **sample_names** (*list*): Sample names to plot
* **output_dir** (*str*, optional): Base directory for saving (default: current)
* **highlight_genes** (*dict*, optional): Dict mapping cell_type → [genes_to_highlight]
* **dpi** (*int*, default=150): Image resolution
* ****plot_kwargs**: Additional styling arguments:

  * ``figsize``: Figure size (default: (15, 8))
  * ``fontsize_labels``: Axis label font size (default: 20)
  * ``fontsize_title``: Title font size (default: 30)
  * ``fontsize_legend``: Legend font size (default: 14)
  * ``fontsize_annotations``: Gene label font size (default: 18)

**Returns**:

* None (saves plots to disk)

**Example**::

    import sctop as top
    
    # Analyze contributions
    contrib_results = top.analyze_sample_contributions(
        sample_data_dict={'cluster1': data1, 'cluster2': data2},
        basis=basis,
        cell_types=['T cell', 'B cell', 'Macrophage'],
        n_top_genes=20
    )
    
    # Highlight known markers
    markers = {
        'T cell': ['CD3D', 'CD3E', 'CD8A'],
        'B cell': ['CD19', 'MS4A1', 'CD79A'],
        'Macrophage': ['CD68', 'CD14', 'FCGR3A']
    }
    
    # Generate all plots
    top.plot_all_contributions(
        results=contrib_results,
        sample_names=['cluster1', 'cluster2'],
        output_dir='./gene_contribution_plots',
        highlight_genes=markers,
        dpi=200,
        figsize=(16, 10),
        fontsize_title=24
    )

**Output Structure**::

    output_dir/
    ├── gene_contributions_T cell/
    │   ├── cluster1.png
    │   └── cluster2.png
    ├── gene_contributions_B cell/
    │   ├── cluster1.png
    │   └── cluster2.png
    └── gene_contributions_Macrophage/
        ├── cluster1.png
        └── cluster2.png

**Interpreting Plots**:

* **X-axis**: Mean gene expression across samples
* **Y-axis**: Mean contribution to cell type score
* **Gray points**: All genes
* **Blue points**: Top N contributing genes (with labels)
* **Red points**: Highlighted genes (if specified)

High contribution + high expression → strong positive marker

High contribution + low expression → potential issue (check basis)

High expression + low contribution → not discriminative for this type

create_colorbar
===============

.. autofunction:: sctop.create_colorbar

**Purpose**: Create a colorbar for continuous data visualization (helper function).

**Parameters**:

* **data** (*array-like*): Data values to create colorbar for
* **label** (*str*): Colorbar label
* **colormap** (*str*, default='rocket_r'): Matplotlib colormap name
* **ax** (*matplotlib.axes.Axes*, optional): Axes to add colorbar to

**Returns**:

* **matplotlib.colors.Colormap**: The colormap object

**Example**::

    import matplotlib.pyplot as plt
    import seaborn as sns
    from sctop.visualization import create_colorbar
    
    fig, ax = plt.subplots()
    
    # Create heatmap with colorbar
    cmap = create_colorbar(
        data=values,
        label='Expression Level',
        colormap='viridis',
        ax=ax
    )
    
    # Use cmap for plotting
    scatter = ax.scatter(x, y, c=values, cmap=cmap)

**Notes**:

* Used internally by other plotting functions
* Supports all matplotlib colormaps
* Automatically normalizes data range

Visualization Tips
==================

Figure Layout
-------------

For publications::

    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Set style
    sns.set_style('whitegrid')
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 12
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot on each axis
    top.plot_highest(scores1, ax=axes[0,0], fontsize=14)
    top.plot_highest(scores2, ax=axes[0,1], fontsize=14)
    # ...
    
    plt.tight_layout()
    plt.savefig('figure.pdf', dpi=300, bbox_inches='tight')

