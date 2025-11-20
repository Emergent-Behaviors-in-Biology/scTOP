import matplotlib.pyplot as plt
import seaborn as sns
import os
from typing import Optional, List, Dict
from .utils import *
from .processing import *

def create_colorbar(data, label, colormap='rocket_r', ax = None):
    ax = ax or plt.gca()
    
    cmap = plt.get_cmap(colormap)
    scalarmap = plt.cm.ScalarMappable(norm=plt.Normalize(min(data), max(data)),
                               cmap=cmap)
    scalarmap.set_array([])
    plt.colorbar(scalarmap, label=label, ax = ax)
    
    return cmap
    
def plot_highest(projections, n=10, ax=None, color="olive", fontsize=40, **kwargs):
    """
    Plots a horizontal bar chart of the top N projections with a fixed x-axis scale.
    """
    ax = ax or plt.gca()
    
    projections_sorted = projections.sort_values(by=projections.columns[0])
    projections_top_n = projections_sorted.iloc[-n:]
    projections_top_n.plot.barh(ax=ax, color=color, legend=False, **kwargs)

    # --- Adjustments for Presentation ---
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    
    xlabel = projections.columns[0]
    ylabel = projections.index.name or 'Items'
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    
    ax.set_xlim(0.0, 1.0)

    ax.grid(axis='x', linestyle='--', alpha=0.6)
    return ax
    
def plot_expression_distribution(scores, n=10, ax=None, box_color="skyblue", fontsize=30, **kwargs):
    """
    Plots boxplots of expression for top genes with a fixed y-axis scale.
    """
    ax = ax or plt.gca()

    gene_meds = scores.median(axis=1).sort_values(ascending=False)
    top_n_genes = gene_meds.head(n).index

    data_to_plot_genes = scores.loc[top_n_genes].T
    melted_data = data_to_plot_genes.melt(var_name="", value_name="Expression")
    order_genes = top_n_genes.tolist()

    sns.boxplot(
        x="",
        y="Expression",
        data=melted_data,
        ax=ax,
        color=box_color,
        order=order_genes,
        showfliers=True,
        **kwargs
    )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right', fontsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)

    ax.set_ylabel("Expression", fontsize=fontsize)
    
    
    return ax
    
def plot_two(projections, celltype1, celltype2, 
             gene=None, gene_expressions=None, ax=None, **kwargs):
    
    ax = ax or plt.gca()

    if gene:
        palette = create_colorbar(gene_expressions.loc[gene],
                                  '{} expression'.format(gene), ax = ax)
        
        plot = sns.scatterplot(x = projections.loc[celltype1],
                        y = projections.loc[celltype2],
                        hue = gene_expressions.loc[gene],
                        palette = palette,
                        alpha = 0.5,
                        ax = ax,
                        **kwargs
                       )
        plot.legend_.remove()
        
    else:
        sns.scatterplot(x = projections.loc[celltype1],
                        y = projections.loc[celltype2],
                        alpha=0.5,
                        ax=ax,
                        **kwargs
                       )
        
def plot_all_contributions(
    results: Dict[str, Dict],
    sample_names: List[str],
    output_dir: Optional[str] = None,
    highlight_genes: Optional[Dict[str, List[str]]] = None,
    dpi: int = 150,
    **plot_kwargs
) -> None:
    """
    Create and save contribution plots for all cell types and samples.
    
    Parameters
    ----------
    results : dict
        Results from analyze_sample_contributions
    sample_names : list
        List of sample names to plot
    output_dir : str, optional
        Base directory for saving plots. If None, uses current directory
    highlight_genes : dict, optional
        Dictionary mapping cell_type -> [genes_to_highlight]
    dpi : int
        DPI for saved images
    **plot_kwargs
        Additional kwargs passed to plot_gene_contribution_scatter
    """
    if output_dir is None:
        output_dir = "."
    
    # Compute predictivity from first cell type (should be same for all)
    basis_needed = False  # We'll get predictivity from results
    
    for cell_type, cell_data in results.items():
        # Create directory for this cell type
        cell_type_dir = os.path.join(output_dir, f"gene_contributions_{cell_type}")
        os.makedirs(cell_type_dir, exist_ok=True)
        
        print(f"Creating plots for {cell_type}...")
        
        for sample_name in sample_names:
            if sample_name not in cell_data['expressions']:
                print(f"  Warning: {sample_name} not found in results, skipping")
                continue
            
            # Get data
            expression = cell_data['expressions'][sample_name]
            contributions = cell_data['contributions'][sample_name]
            top_genes = cell_data['top_genes'][sample_name]
            
            # Compute predictivity from contributions and expression
            # predictivity = contributions / expression (approximately)
            # But better to pass it separately
            common_genes = contributions.index
            mean_expression = expression.loc[common_genes].mean(axis=1)
            mean_contribution = contributions.mean(axis=1)
            
            # Approximate predictivity for this cell type
            predictivity_approx = mean_contribution / (mean_expression + 1e-10)
            
            # Get highlight genes for this cell type
            highlight = None
            if highlight_genes and cell_type in highlight_genes:
                highlight = highlight_genes[cell_type]
            
            # Create plot
            fig, ax = plt.subplots(figsize=plot_kwargs.get('figsize', (15, 8)))
            
            # Use mean expression and mean contribution for plotting
            ax.scatter(mean_expression, mean_contribution, 
                      color='gray', alpha=0.6, label='All Genes', s=3)
            
            # Highlight top genes
            ax.scatter(mean_expression[top_genes], mean_contribution[top_genes],
                      color='blue', label=f'Top {len(top_genes)} Genes', s=10)
            
            # Annotate
            texts = []
            for gene in top_genes:
                x = mean_expression[gene]
                y = mean_contribution[gene]
                texts.append(ax.text(x, y, gene, 
                                   fontsize=plot_kwargs.get('fontsize_annotations', 18)))
            
            
            # Highlight special genes
            if highlight:
                for gene in highlight:
                    if gene in mean_expression.index:
                        x = mean_expression[gene]
                        y = mean_contribution[gene]
                        ax.scatter(x, y, color='red', s=20, zorder=10)
                        ax.text(x, y, gene, fontsize=18, color='red', weight='bold')
            
            # Formatting
            fontsize_labels = plot_kwargs.get('fontsize_labels', 20)
            fontsize_title = plot_kwargs.get('fontsize_title', 30)
            fontsize_legend = plot_kwargs.get('fontsize_legend', 14)
            
            ax.set_xlabel('Mean Gene Expression', fontsize=fontsize_labels)
            ax.set_ylabel(f'Mean Contribution to {cell_type} Score', fontsize=fontsize_labels)
            ax.set_title(f'{sample_name}', fontsize=fontsize_title)
            ax.legend(fontsize=fontsize_legend)
            ax.grid(True, alpha=0.3)
            
            # Save
            
            filename = os.path.join(cell_type_dir, f"{sample_name}.png")
            plt.savefig(filename, dpi=dpi, bbox_inches='tight')
            plt.show()
            plt.close(fig)
        
        print(f"  Saved plots to {cell_type_dir}/")