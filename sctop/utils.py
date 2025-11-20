import numpy as np
import pandas as pd
from typing import Tuple, Optional, Dict, List
from sklearn.model_selection import train_test_split, KFold
from sklearn.feature_selection import f_classif
from sklearn.metrics import f1_score, precision_score, recall_score, ConfusionMatrixDisplay
from sklearn.preprocessing import StandardScaler
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import multiprocessing as mp
import anndata as ad
from functools import partial
import matplotlib.pyplot as plt
import seaborn as sns
from .processing import *

from typing import Optional, Union, List, Tuple

    
def perform_anova_selection(
    basis: pd.DataFrame,
    adata: ad.AnnData,
    training_IDs: np.ndarray,
    cell_type_column: str,
    n_features: int = 2000,
    percentile: Optional[float] = None,
    standardize: bool = True
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Perform ANOVA feature selection on the basis and optionally standardize.
    
    Parameters:
    -----------
    basis : pd.DataFrame
        The basis matrix (genes x cell types)
    adata : ad.AnnData
        The AnnData object
    training_IDs : np.ndarray
        Training cell IDs
    cell_type_column : str
        Column name for cell types
    n_features : int
        Number of top features to select (if percentile is None)
    percentile : float, optional
        Percentile of features to keep (overrides n_features)
    standardize : bool
        Whether to standardize the basis after selection (default: True)
    
    Returns:
    --------
    basis_selected : pd.DataFrame
        Basis with selected features only (and standardized if requested)
    selected_genes : np.ndarray
        Array of selected gene names
    """
    print(f"\nPerforming ANOVA feature selection...")
    
    # Get training data
    adata_train = adata[training_IDs, :]
    
    # Convert to dense if sparse
    if hasattr(adata_train.X, 'toarray'):
        X_train = adata_train.X.toarray()
    else:
        X_train = adata_train.X
    
    # Get labels
    y_train = adata_train.obs[cell_type_column].values
    
    # Compute ANOVA F-statistic
    print("Computing F-statistics...")
    f_scores, p_values = f_classif(X_train, y_train)
    
    # Handle NaN values (can occur with zero-variance features)
    f_scores = np.nan_to_num(f_scores, nan=0.0)
    
    # Select features
    if percentile is not None:
        threshold = np.percentile(f_scores, 100 - percentile)
        selected_mask = f_scores >= threshold
        n_selected = selected_mask.sum()
    else:
        n_selected = min(n_features, len(f_scores))
        top_indices = np.argsort(f_scores)[-n_selected:]
        selected_mask = np.zeros(len(f_scores), dtype=bool)
        selected_mask[top_indices] = True
    
    selected_genes = adata.var.index[selected_mask].values
    
    print(f"Selected {n_selected} genes out of {len(f_scores)}")
    print(f"F-score range: [{f_scores[selected_mask].min():.2f}, {f_scores[selected_mask].max():.2f}]")
    
    # Filter basis
    basis_selected = basis.loc[selected_genes]
    
    # Standardize if requested
    if standardize:
        print("Standardizing basis (scaling to unit norm per column)...")
        
        # Apply StandardScaler per column (cell type)
        # This ensures each cell type basis vector has mean=0 and std=1
        scaler = StandardScaler()
        basis_scaled = scaler.fit_transform(basis_selected.values)
        
        # Convert back to DataFrame
        basis_selected = pd.DataFrame(
            basis_scaled,
            index=basis_selected.index,
            columns=basis_selected.columns
        )
        
        # Re-normalize each column to unit length
        # This ensures basis.T @ basis has 1s on diagonal
        norms = np.linalg.norm(basis_selected.values, axis=0, keepdims=True)
        basis_selected = pd.DataFrame(
            basis_selected.values / norms,
            index=basis_selected.index,
            columns=basis_selected.columns
        )
        
        # Verify normalization
        gram_matrix = basis_selected.T @ basis_selected
        diagonal_values = np.diag(gram_matrix.values)
        print(f"Diagonal of basis.T @ basis: min={diagonal_values.min():.6f}, max={diagonal_values.max():.6f}")
        print(f"All diagonal values ≈ 1.0: {np.allclose(diagonal_values, 1.0, atol=1e-5)}")
    
    return basis_selected, selected_genes

def _create_basis_with_cv(
    adata: ad.AnnData,
    cell_type_column: str,
    threshold: int,
    random_state: int,
    n_jobs: int,
    do_anova: bool,
    n_features: int,
    anova_percentile: Optional[float],
    spec_value: float,
    outer_chunks: int,
    inner_chunk_size: int,
    n_scoring_jobs: int,
    cv_folds: int,
    plot_results: bool = True
) -> Dict:
    """
    Create basis with cross-validation.
    """
    print("\n" + "="*60)
    print(f"Creating basis with {cv_folds}-fold cross-validation")
    print("="*60)
    
    atlas_metadata = adata.obs
    type_counts = atlas_metadata[cell_type_column].value_counts().sort_index()
    types_above_threshold = type_counts[type_counts > threshold].index.tolist()
    
    # Get all cells for cell types above threshold
    valid_cells = atlas_metadata[atlas_metadata[cell_type_column].isin(types_above_threshold)].index.values
    
    print(f"\nTotal cells to use: {len(valid_cells)}")
    print(f"Cell types: {len(types_above_threshold)}")
    
    kf = KFold(n_splits=cv_folds, shuffle=True, random_state=random_state)
    
    cv_results = []
    all_true_labels = []
    all_predicted_labels = []
    
    for fold_idx, (train_idx, test_idx) in enumerate(kf.split(valid_cells)):
        print(f"\n{'='*60}")
        print(f"Fold {fold_idx + 1}/{cv_folds}")
        print(f"{'='*60}")
        
        train_IDs = valid_cells[train_idx]
        test_IDs = valid_cells[test_idx]
        
        print(f"Training cells: {len(train_IDs)}")
        print(f"Test cells: {len(test_IDs)}")
        
        # Create basis for this fold
        fold_basis = _create_fold_basis(
            adata=adata,
            train_IDs=train_IDs,
            cell_type_column=cell_type_column,
            types_above_threshold=types_above_threshold,
            n_jobs=n_jobs
        )
        
        # Optional ANOVA
        if do_anova:
            fold_basis, _ = perform_anova_selection(
                basis=fold_basis,
                adata=adata,
                training_IDs=train_IDs,
                cell_type_column=cell_type_column,
                n_features=n_features,
                percentile=anova_percentile
            )
        
        # Score
        cell_accuracies, true_labels, predicted_labels, accuracies = run_scoring_parallel(
            adata=adata,
            basis=fold_basis,
            test_IDs=test_IDs,
            cell_type_column=cell_type_column,
            spec_value=spec_value,
            outer_chunks=outer_chunks,
            inner_chunk_size=inner_chunk_size,
            n_jobs=n_scoring_jobs
        )
        
        # Calculate metrics for this fold
        fold_metrics = calculate_metrics(true_labels, predicted_labels, len(test_IDs), accuracies)
        
        cv_results.append({
            'fold': fold_idx + 1,
            'metrics': fold_metrics,
            'train_size': len(train_IDs),
            'test_size': len(test_IDs)
        })
        
        all_true_labels.extend(true_labels)
        all_predicted_labels.extend(predicted_labels)
        
        print(f"\nFold {fold_idx + 1} Results:")
        print_metrics(fold_metrics)
    
    # Aggregate CV results
    print("\n" + "="*60)
    print("CROSS-VALIDATION SUMMARY")
    print("="*60)
    
    avg_metrics = {}
    metric_keys = cv_results[0]['metrics'].keys()
    
    for key in metric_keys:
        values = [fold['metrics'][key] for fold in cv_results]
        avg_metrics[f'{key}_mean'] = np.mean(values)
        avg_metrics[f'{key}_std'] = np.std(values)
    
    print("\nAverage Metrics Across Folds:")
    for key in ['accuracy', 'top3_accuracy', 'f1_macro', 'f1_weighted']:
        if f'{key}_mean' in avg_metrics:
            print(f"{key}: {avg_metrics[f'{key}_mean']:.4f} ± {avg_metrics[f'{key}_std']:.4f}")
        
    # Create final basis on all data
    print("\n" + "="*60)
    print("Creating final basis on all data")
    print("="*60)
    
    final_basis = _create_fold_basis(
        adata=adata,
        train_IDs=valid_cells,
        cell_type_column=cell_type_column,
        types_above_threshold=types_above_threshold,
        n_jobs=n_jobs
    )
    
    selected_genes = None
    if do_anova:
        final_basis, selected_genes = perform_anova_selection(
            basis=final_basis,
            adata=adata,
            training_IDs=valid_cells,
            cell_type_column=cell_type_column,
            n_features=n_features,
            percentile=anova_percentile
        )
    
    # Overall confusion matrix
    unique_labels = sorted(set(all_true_labels + all_predicted_labels))
    from sklearn.metrics import confusion_matrix
    cm = confusion_matrix(all_true_labels, all_predicted_labels, labels=unique_labels)

    # F1 scores per cell type
    f1_scores = f1_score(all_true_labels, all_predicted_labels, average=None, labels=unique_labels)
    
    # Create DataFrame for easier sorting/plotting
    f1_df = pd.DataFrame({
        'Cell Type': unique_labels,
        'F1 Score': f1_scores
    }).sort_values('F1 Score', ascending=True)

    if plot_results:
        print("\nGenerating performance plots...")
        plot_performance_summary(true_labels, predicted_labels, f1_df)
    
    return {
        'basis': final_basis,
        'selected_genes': selected_genes,
        'cv_results': cv_results,
        'cv_avg_metrics': avg_metrics,
        'confusion_matrix': cm,
        'confusion_matrix_labels': unique_labels,
        'all_true_labels': all_true_labels,
        'all_predicted_labels': all_predicted_labels,
        'f1_scores': f1_df
    }


def _create_fold_basis(
    adata: ad.AnnData,
    train_IDs: np.ndarray,
    cell_type_column: str,
    types_above_threshold: List[str],
    n_jobs: int
) -> pd.DataFrame:
    """
    Helper function to create basis for a single fold.
    """
    from functools import partial
    
    if n_jobs == -1:
        n_jobs = mp.cpu_count()
    
    def process_cell_type(cell_type):
        # Get only training cells of this type
        cell_IDs = adata.obs.loc[train_IDs][adata.obs.loc[train_IDs, cell_type_column] == cell_type].index
        
        if len(cell_IDs) == 0:
            return None
        
        adata_train = adata[cell_IDs, :]
        cell_data_sparse = adata_train.X.T
        
        cell_data_df = pd.DataFrame(
            cell_data_sparse.toarray(),
            index=adata.var.index,
            columns=adata_train.obs.index
        )
        
        # Import process function (assuming it's in scope)
        processed = process(cell_data_df, average=True)
        
        return {
            'cell_type': cell_type,
            'processed': processed
        }
    
    with ThreadPoolExecutor(max_workers=min(n_jobs, len(types_above_threshold))) as executor:
        results = list(tqdm(
            executor.map(process_cell_type, types_above_threshold),
            total=len(types_above_threshold),
            desc="Creating basis"
        ))
    
    # Filter out None results and combine
    results = [r for r in results if r is not None]
    basis_list = [r['processed'] for r in results]
    cell_types = [r['cell_type'] for r in results]
    
    basis = pd.concat(basis_list, axis=1)
    basis.columns = cell_types
    
    return basis


def calculate_metrics(
    true_labels: List,
    predicted_labels: List,
    total_cells: int,
    accuracies: Dict
) -> Dict:
    """
    Calculate comprehensive metrics.
    """
    metrics = {
        'accuracy': accuracies['top1'] / total_cells,
        'top3_accuracy': accuracies['top3'] / total_cells,
        'unspecified_rate': accuracies['unspecified'] / total_cells,
        'f1_macro': f1_score(true_labels, predicted_labels, average='macro', zero_division=0),
        'f1_weighted': f1_score(true_labels, predicted_labels, average='weighted', zero_division=0),
        'precision_macro': precision_score(true_labels, predicted_labels, average='macro', zero_division=0),
        'precision_weighted': precision_score(true_labels, predicted_labels, average='weighted', zero_division=0),
        'recall_macro': recall_score(true_labels, predicted_labels, average='macro', zero_division=0),
        'recall_weighted': recall_score(true_labels, predicted_labels, average='weighted', zero_division=0),
        'total_cells': total_cells
    }
    
    return metrics


def calculate_per_cell_type_accuracy(cell_accuracies: Dict) -> pd.DataFrame:
    """
    Calculate per cell type accuracy.
    """
    cell_acc_data = []
    for cell_id, data in cell_accuracies.items():
        cell_acc_data.append({
            'cell_id': cell_id,
            'true': data['true'],
            'predicted': data['predicted'],
            'top1': data['top1'],
            'top3': data['top3'],
            'unspecified': data['unspecified']
        })
    
    df = pd.DataFrame(cell_acc_data)
    
    per_type = df.groupby('true').agg({
        'top1': ['sum', 'count', 'mean'],
        'top3': ['sum', 'mean'],
        'unspecified': ['sum', 'mean']
    })
    
    per_type.columns = ['_'.join(col).strip() for col in per_type.columns.values]
    per_type = per_type.rename(columns={
        'top1_sum': 'correct',
        'top1_count': 'total',
        'top1_mean': 'accuracy',
        'top3_sum': 'top3_correct',
        'top3_mean': 'top3_accuracy',
        'unspecified_sum': 'unspecified_count',
        'unspecified_mean': 'unspecified_rate'
    })
    
    per_type = per_type.sort_values('accuracy', ascending=False)
    
    return per_type


def print_metrics(metrics: Dict):
    """
    Pretty print metrics.
    """
    print(f"\nAccuracy (Top-1): {metrics['accuracy']:.4f}")
    print(f"Top-3 Accuracy: {metrics['top3_accuracy']:.4f}")
    print(f"Unspecified Rate: {metrics['unspecified_rate']:.4f}")
    print(f"F1 Score (Macro): {metrics['f1_macro']:.4f}")
    print(f"F1 Score (Weighted): {metrics['f1_weighted']:.4f}")
    print(f"Precision (Macro): {metrics['precision_macro']:.4f}")
    print(f"Precision (Weighted): {metrics['precision_weighted']:.4f}")
    print(f"Recall (Macro): {metrics['recall_macro']:.4f}")
    print(f"Recall (Weighted): {metrics['recall_weighted']:.4f}")


def create_basis_optimized(
    adata: ad.AnnData,
    cell_type_column: str,
    threshold: int,
    test_size: float,
    random_state: int,
    n_jobs: int = -1
) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """
    Original function - kept for backwards compatibility.
    """
    print("\nBuilding cell type basis (parallelized)...")
    atlas_metadata = adata.obs
    type_counts = atlas_metadata[cell_type_column].value_counts().sort_index()
    types_above_threshold = type_counts[type_counts > threshold].index.tolist()

    if n_jobs == -1:
        n_jobs = mp.cpu_count()
    
    def process_cell_type(cell_type):
        cell_IDs = atlas_metadata[atlas_metadata[cell_type_column] == cell_type].index
        X_train, X_test = train_test_split(cell_IDs, test_size=test_size, random_state=random_state)
        
        adata_train = adata[X_train, :]
        cell_data_sparse = adata_train.X.T
        
        cell_data_df = pd.DataFrame(
            cell_data_sparse.toarray(),
            index=adata.var.index,
            columns=adata_train.obs.index
        )
        
        processed = process(cell_data_df, average=True)
        
        return {
            'cell_type': cell_type,
            'processed': processed,
            'X_train': X_train,
            'X_test': X_test
        }
    
    with ThreadPoolExecutor(max_workers=min(n_jobs, len(types_above_threshold))) as executor:
        results = list(tqdm(
            executor.map(process_cell_type, types_above_threshold),
            total=len(types_above_threshold),
            desc="Creating basis (parallel)"
        ))
    
    basis_list = [r['processed'] for r in results]
    training_IDs = np.concatenate([r['X_train'] for r in results])
    test_IDs = np.concatenate([r['X_test'] for r in results])
    
    basis = pd.concat(basis_list, axis=1)
    basis.columns = types_above_threshold
    
    print(f"Basis created: {basis.shape}")
    print(f"Total training cells: {len(training_IDs)}")
    print(f"Total test cells: {len(test_IDs)}")
    
    return basis, training_IDs, test_IDs

def run_scoring_parallel(
    adata: ad.AnnData,
    basis: pd.DataFrame,
    test_IDs: np.ndarray,
    cell_type_column: str,
    spec_value: float,
    outer_chunks: int,
    inner_chunk_size: int,
    n_jobs: int = 4
) -> Tuple[dict, list, list, dict]:
    """
    OPTIMIZED: Parallel scoring of test cells.
    Uses ThreadPoolExecutor for shared-memory parallel processing.
    """
    print(f"\nScoring {len(test_IDs)} test cells (parallel, {n_jobs} workers)...")
    
    split_IDs = np.array_split(test_IDs, outer_chunks)
    
    score_func = partial(
        score_chunk_optimized,
        adata,  
        basis,  
        cell_type_column=cell_type_column,
        spec_value=spec_value,
        inner_chunk_size=inner_chunk_size
    )
    

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        results = list(tqdm(
            executor.map(score_func, split_IDs),
            total=len(split_IDs),
            desc="Scoring chunks (parallel)"
        ))
    
    # Aggregate results
    all_cell_accuracies = {}
    all_true_labels = []
    all_predicted_labels = []
    total_accuracies = {'top1': 0, 'top3': 0, 'unspecified': 0}
    
    for cell_acc, true_lab, pred_lab, acc in results:
        all_cell_accuracies.update(cell_acc)
        all_true_labels.extend(true_lab)
        all_predicted_labels.extend(pred_lab)
        for key in total_accuracies:
            total_accuracies[key] += acc[key]
    
    return all_cell_accuracies, all_true_labels, all_predicted_labels, total_accuracies

def score_chunk_optimized(
    adata: ad.AnnData,
    basis: pd.DataFrame,
    sample_IDs: np.ndarray,
    cell_type_column: str,
    spec_value: float,
    inner_chunk_size: int
) -> Tuple[dict, list, list, dict]:
    """
    OPTIMIZED: Score a single chunk of cells.
    Extracted for parallel processing.
    """
    if len(sample_IDs) == 0:
        return {}, [], [], {'top1': 0, 'top3': 0, 'unspecified': 0}
    
    # Slice sparse data
    adata_test_chunk = adata[sample_IDs, :]
    test_data_sparse = adata_test_chunk.X.T
    
    # Convert to dense DataFrame
    test_data_df = pd.DataFrame(
        test_data_sparse.toarray(),
        index=adata.var.index,
        columns=adata_test_chunk.obs.index
    )
    
    # Process and score with internal chunking
    test_processed = process(test_data_df, chunk_size=inner_chunk_size)
    test_projections = score(basis, test_processed, chunk_size=inner_chunk_size)
    
    # Calculate metrics
    cell_accuracies = {}
    predicted_labels = []
    true_labels = []
    accuracies = {'top1': 0, 'top3': 0, 'unspecified': 0}
    
    atlas_metadata = adata.obs
    
    for sample_id, sample_projections in test_projections.items():
        types_sorted = sample_projections.sort_values(ascending=False).index
        true_type = atlas_metadata.loc[sample_id, cell_type_column]
        
        true_labels.append(true_type)
        top_type = types_sorted[0]
        predicted_labels.append(top_type)
        
        is_unspecified = sample_projections.max() < spec_value
        is_top1 = top_type == true_type
        is_top3 = true_type in types_sorted[:3]
        
        if is_unspecified:
            accuracies['unspecified'] += 1
        if is_top1:
            accuracies['top1'] += 1
        if is_top3:
            accuracies['top3'] += 1
        
        cell_accuracies[sample_id] = {
            'true': true_type,
            'predicted': top_type,
            'top1': int(is_top1),
            'top3': int(is_top3),
            'unspecified': int(is_unspecified)
        }
    
    return cell_accuracies, true_labels, predicted_labels, accuracies

def plot_performance_summary(
    true_labels: List, 
    predicted_labels: List, 
    figsize_base: int = 10,
    f1_df: Optional[pd.DataFrame] = None
):
    """
    Generates and displays a Confusion Matrix and a Per-Cell-Type F1 Score plot.
    """
    unique_labels = sorted(list(set(true_labels) | set(predicted_labels)))
    n_classes = len(unique_labels)
    
    # Dynamic figure size based on number of cell types
    # Scale factor ensures labels don't overlap for large datasets
    scale_factor = max(1, n_classes / 20) 
    fig, (ax1, ax2) = plt.subplots(
        1, 2, 
        figsize=(figsize_base * 2 * scale_factor, figsize_base * scale_factor),
        gridspec_kw={'width_ratios': [1, 1]}
    )

    # --- 1. Confusion Matrix ---
    ConfusionMatrixDisplay.from_predictions(
        true_labels, 
        predicted_labels, 
        labels=unique_labels,
        xticks_rotation='vertical',
        normalize='true', # Normalize over true labels (Recall)
        values_format=".2f",
        cmap='Blues',
        ax=ax1,
        colorbar=False 
    )
    ax1.set_title(f"Normalized Confusion Matrix ({n_classes} Cell Types)", fontsize=16)
    ax1.grid(False)

    # --- 2. F1 Scores per Cell Type ---

    # Plot
    sns.barplot(data=f1_df, x='F1 Score', y='Cell Type', ax=ax2)
    
    ax2.set_title("F1 Score per Cell Type", fontsize=16)
    ax2.set_xlabel("F1 Score", fontsize=14)
    ax2.set_xlim(0, 1.05)
    ax2.grid(axis='x', linestyle='--', alpha=0.7)
    
    # Add value labels to the end of bars
    for i, v in enumerate(f1_df['F1 Score']):
        ax2.text(v + 0.01, i, f"{v:.2f}", va='center', fontsize=10)

    plt.tight_layout()
    plt.show()

def compute_predictivity(basis: pd.DataFrame) -> pd.DataFrame:
    """
    Compute predictivity matrix from basis.
    
    The predictivity shows how each gene contributes to each cell type's score.
    Formula: predictivity = inv(B^T @ B) @ B^T
    
    Parameters
    ----------
    basis : pd.DataFrame
        Basis matrix (genes x cell_types)
        
    Returns
    -------
    predictivity : pd.DataFrame
        Predictivity matrix (cell_types x genes)
        Shows how each gene contributes to each cell type score
    """
    
    correlation_matrix = basis.T @ basis
    
    eta = np.linalg.inv(correlation_matrix.values) @ basis.T.values
    
    predictivity = pd.DataFrame(
        eta,
        index=basis.columns,  # Cell types
        columns=basis.index   # Genes
    )
    
    return predictivity

def compute_gene_contributions(
    data: Union[pd.DataFrame, np.ndarray],
    basis: pd.DataFrame,
    predictivity: Optional[pd.DataFrame] = None,
    cell_types: Optional[List[str]] = None,
    process_data: bool = True,
) -> Dict[str, pd.DataFrame]:
    """
    Compute gene-level contributions to cell type scores.
    
    For each cell type, computes: contribution = expression * predictivity
    
    Parameters
    ----------
    data : DataFrame or array
        Expression data (genes x samples)
    basis : pd.DataFrame
        Basis matrix
    predictivity : pd.DataFrame, optional
        Precomputed predictivity matrix. If None, computed from basis
    cell_types : list, optional
        Cell types to compute contributions for. If None, uses all
    process_data : bool
        Whether to process the data first (default: True)
        
    Returns
    -------
    contributions : dict
        Dictionary mapping cell_type -> contribution_matrix (genes x samples)
    """
    # Compute predictivity if not provided
    if predictivity is None:
        predictivity = compute_predictivity(basis)
    
    # Process data if requested
    if process_data:
        if not isinstance(data, pd.DataFrame):
            data = pd.DataFrame(data)
        expression = process(data)
    else:
        expression = data if isinstance(data, pd.DataFrame) else pd.DataFrame(data)
    
    # Determine cell types
    if cell_types is None:
        cell_types = basis.columns.tolist()
    
    # Find common genes
    common_genes = np.intersect1d(expression.index, predictivity.columns)
    
    if len(common_genes) == 0:
        raise ValueError("No common genes between data and predictivity matrix")
    
    # Compute contributions for each cell type
    contributions = {}
    for cell_type in cell_types:
        # contribution = expression * predictivity
        contrib = expression.loc[common_genes].multiply(
            predictivity.loc[cell_type, common_genes],
            axis=0
        )
        contributions[cell_type] = contrib
    
    return contributions

def find_top_contributing_genes(
    contributions: pd.DataFrame,
    n_genes: int = 20,
    aggregate: str = 'mean'
) -> pd.Series:
    """
    Find top contributing genes from contribution matrix.
    
    Parameters
    ----------
    contributions : pd.DataFrame
        Gene contributions (genes x samples)
    n_genes : int
        Number of top genes to return
    aggregate : str
        How to aggregate across samples: 'mean', 'median', 'max'
        
    Returns
    -------
    top_genes : pd.Series
        Top contributing genes with their aggregated scores
    """
    if aggregate == 'mean':
        gene_scores = contributions.mean(axis=1)
    elif aggregate == 'median':
        gene_scores = contributions.median(axis=1)
    elif aggregate == 'max':
        gene_scores = contributions.max(axis=1)
    else:
        raise ValueError(f"Unknown aggregate method: {aggregate}")
    
    return gene_scores.sort_values(ascending=False).head(n_genes)
