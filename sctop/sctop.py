import numpy as np
import pandas as pd
import anndata as ad
import tempfile
import requests
from pathlib import Path

from typing import Optional, Union, List, Tuple, Dict
from .utils import *
from .processing import *

# Define URLs for premade bases
BASIS_URLS = {
    "MCKO legacy": "https://figshare.com/ndownloader/files/59736809",
}

def list_available_bases() -> List[str]:
    """
    List available premade bases that can be loaded.
    
    Returns:
    --------
    basis_keys : list
        List of available basis keys
    """
    return list(BASIS_URLS.keys())

def load_basis(
    basis_key: str,
    cache_dir: Optional[Union[str, Path]] = None,
    force_download: bool = False
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load a basis from an h5ad file hosted online.
    
    Parameters:
    -----------
    basis_key : str
        Name/key of the basis to load (e.g., "MCKO legacy")
    cache_dir : str, optional
        Directory to cache downloaded files. If None, uses system temp directory.
    force_download : bool
        If True, re-downloads even if cached file exists
    
    Returns:
    --------
    basis : pd.DataFrame
        Basis matrix (genes x cell types)
    metadata : pd.DataFrame
        Metadata for the basis (cell types x attributes)
    
    Example:
    --------
    >>> basis, metadata = load_basis(
    ...     basis_key="MCKO legacy"
    ... )
    """
    # Setup cache directory
    if cache_dir is None:
        cache_dir = Path(tempfile.gettempdir()) / "sctop_cache"
    else:
        cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    # Create filename from URL and basis key
    url = BASIS_URLS.get(basis_key)
    if url is None:
        raise ValueError(f"Unknown basis_key '{basis_key}'. Available keys: {list(BASIS_URLS.keys())}")

    cache_file = cache_dir / f"basis_{basis_key}.h5ad"
    
    # Download if needed
    if force_download or not cache_file.exists():
        print(f"Downloading basis from {url}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        with open(cache_file, 'wb') as f:
            if total_size > 0:
                downloaded = 0
                chunk_size = 8192
                for chunk in response.iter_content(chunk_size=chunk_size):
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        progress = (downloaded / total_size) * 100
                        print(f"\rProgress: {progress:.1f}%", end='', flush=True)
                print()  # New line after progress
            else:
                f.write(response.content)
        print(f"Cached to {cache_file}")
    else:
        print(f"Loading from cache: {cache_file}")
    
    # Load the h5ad file
    adata_basis = ad.read_h5ad(cache_file)
    
    # Extract basis matrix and metadata
    # Assuming the basis is stored as a layer or in X with cell types as obs
    if 'cell_type' in adata_basis.obs.columns:
        # Basis stored with cell types as observations
        basis = pd.DataFrame(
            adata_basis.X.T if hasattr(adata_basis.X, 'toarray') else adata_basis.X.T,
            index=adata_basis.var_names,
            columns=adata_basis.obs['cell_type']
        )
        metadata = adata_basis.obs.copy()
    elif 'Cell Type' in adata_basis.obs.columns:
        # Alternative column name
        basis = pd.DataFrame(
            adata_basis.X.T if hasattr(adata_basis.X, 'toarray') else adata_basis.X.T,
            index=adata_basis.var_names,
            columns=adata_basis.obs['Cell Type']
        )
        metadata = adata_basis.obs.copy()
    
    print(f"\nLoaded basis '{basis_key}':")
    print(f"  Genes: {basis.shape[0]:,}")
    print(f"  Cell types: {basis.shape[1]:,}")
    
    return basis, metadata

def create_basis(
    adata: ad.AnnData,
    cell_type_column: str,
    threshold: int,
    test_size: float = 0.2,
    random_state: int = 42,
    n_jobs: int = -1,
    do_anova: bool = False,
    n_features: int = 20000,
    anova_percentile: Optional[float] = None,
    spec_value: float = 0.1,
    outer_chunks: int = 10,
    inner_chunk_size: int = 1000,
    n_scoring_jobs: int = 4,
    cv_folds: Optional[int] = None,
    plot_results: bool = True
) -> Dict:
    """
    Create basis and evaluate with optional ANOVA selection and cross-validation.
    
    Parameters:
    -----------
    adata : ad.AnnData
        Annotated data object
    cell_type_column : str
        Column name for cell types in adata.obs
    threshold : int
        Minimum number of cells per cell type
    test_size : float
        Fraction of data to use for testing (if cv_folds is None)
    random_state : int
        Random seed
    n_jobs : int
        Number of parallel jobs for basis creation
    do_anova : bool
        Whether to perform ANOVA feature selection
    n_features : int
        Number of features to select with ANOVA
    anova_percentile : float, optional
        Percentile of features to keep (overrides n_features)
    spec_value : float
        Threshold for unspecified predictions
    outer_chunks : int
        Number of chunks for parallel scoring
    inner_chunk_size : int
        Chunk size for internal processing
    n_scoring_jobs : int
        Number of parallel jobs for scoring
    cv_folds : int, optional
        Number of cross-validation folds. If None, uses single train-test split
    
    Returns:
    --------
    results : dict
        Dictionary containing:
        - 'basis': final basis
        - 'selected_genes': selected genes (if ANOVA)
        - 'metrics': performance metrics
        - 'cv_results': cross-validation results (if cv_folds is not None)
        - 'confusion_matrix': confusion matrix
        - 'per_cell_type': per cell type accuracy
    """
    
    if cv_folds is not None:
        return _create_basis_with_cv(
            adata=adata,
            cell_type_column=cell_type_column,
            threshold=threshold,
            random_state=random_state,
            n_jobs=n_jobs,
            do_anova=do_anova,
            n_features=n_features,
            anova_percentile=anova_percentile,
            spec_value=spec_value,
            outer_chunks=outer_chunks,
            inner_chunk_size=inner_chunk_size,
            n_scoring_jobs=n_scoring_jobs,
            cv_folds=cv_folds,
            plot_results=plot_results
        )
    
    # Single train-test split
    print("\n" + "="*60)
    print("Creating basis with single train-test split")
    print("="*60)
    
    # Create basis with training data
    basis, training_IDs, test_IDs = create_basis_optimized(
        adata=adata,
        cell_type_column=cell_type_column,
        threshold=threshold,
        test_size=test_size,
        random_state=random_state,
        n_jobs=n_jobs
    )
    
    selected_genes = None
    
    # Optional ANOVA selection
    if do_anova:
        basis, selected_genes = perform_anova_selection(
            basis=basis,
            adata=adata,
            training_IDs=training_IDs,
            cell_type_column=cell_type_column,
            n_features=n_features,
            percentile=anova_percentile
        )
    
    # Score on test set
    cell_accuracies, true_labels, predicted_labels, accuracies = run_scoring_parallel(
        adata=adata,
        basis=basis,
        test_IDs=test_IDs,
        cell_type_column=cell_type_column,
        spec_value=spec_value,
        outer_chunks=outer_chunks,
        inner_chunk_size=inner_chunk_size,
        n_jobs=n_scoring_jobs
    )
    
    # Calculate metrics
    metrics = calculate_metrics(true_labels, predicted_labels, len(test_IDs), accuracies)
    
    # Per cell type accuracy
    per_cell_type = calculate_per_cell_type_accuracy(cell_accuracies)
    
    # Confusion matrix
    unique_labels = sorted(set(true_labels + predicted_labels))
    from sklearn.metrics import confusion_matrix, f1_score
    cm = confusion_matrix(true_labels, predicted_labels, labels=unique_labels)

    # F1 scores per cell type
    f1_scores = f1_score(true_labels, predicted_labels, average=None, labels=unique_labels)
    
    # Create DataFrame for easier sorting/plotting
    f1_df = pd.DataFrame({
        'Cell Type': unique_labels,
        'F1 Score': f1_scores
    }).sort_values('F1 Score', ascending=True)
    
    print("\n" + "="*60)
    print("FINAL RESULTS")
    print("="*60)
    print_metrics(metrics)
    
    if plot_results:
        print("\nGenerating performance plots...")
        plot_performance_summary(true_labels, predicted_labels, f1_df)
        
    return {
        'basis': basis,
        'selected_genes': selected_genes,
        'training_IDs': training_IDs,
        'test_IDs': test_IDs,
        'metrics': metrics,
        'confusion_matrix': cm,
        'confusion_matrix_labels': unique_labels,
        'per_cell_type': per_cell_type,
        'f1_scores': f1_df,
        'cell_accuracies': cell_accuracies,
        'true_labels': true_labels,
        'predicted_labels': predicted_labels
    }


def analyze_sample_contributions(
    sample_data_dict: Dict[str, Union[pd.DataFrame, np.ndarray]],
    basis: pd.DataFrame,
    cell_types: Optional[List[str]] = None,
    n_top_genes: int = 20,
    process_data: bool = True,
) -> Dict[str, Dict]:
    """
    Analyze gene contributions for multiple samples/clusters.
    
    Parameters
    ----------
    sample_data_dict : dict
        Dictionary mapping sample_name -> expression_data
    basis : pd.DataFrame
        Basis matrix
    cell_types : list, optional
        Cell types to analyze. If None, uses all
    n_top_genes : int
        Number of top genes to identify per sample
    process_data : bool
        Whether to process the data
        
    Returns
    -------
    results : dict
        Nested dictionary with structure:
        {cell_type: {
            'contributions': {sample_name: contribution_matrix},
            'top_genes': {sample_name: [gene1, gene2, ...]},
            'expressions': {sample_name: expression_matrix}
        }}
    """
    # Compute predictivity once
    predictivity = compute_predictivity(basis)
    
    if cell_types is None:
        cell_types = basis.columns.tolist()
    
    results = {}
    
    for cell_type in cell_types:
        contributions = {}
        top_genes = {}
        expressions = {}
        
        for sample_name, data in sample_data_dict.items():
            # Process data
            if process_data:
                if not isinstance(data, pd.DataFrame):
                    data = pd.DataFrame(data)
                expression = process(data)
            else:
                expression = data if isinstance(data, pd.DataFrame) else pd.DataFrame(data)
            
            # Compute contributions for this cell type
            common_genes = np.intersect1d(expression.index, predictivity.columns)
            contrib = expression.loc[common_genes].multiply(
                predictivity.loc[cell_type, common_genes],
                axis=0
            )
            
            # Find top genes
            top_gene_list = find_top_contributing_genes(
                contrib, 
                n_genes=n_top_genes
            ).index.tolist()
            
            # Store results
            contributions[sample_name] = contrib
            top_genes[sample_name] = top_gene_list
            expressions[sample_name] = expression
        
        results[cell_type] = {
            'contributions': contributions,
            'top_genes': top_genes,
            'expressions': expressions
        }
    
    return results
