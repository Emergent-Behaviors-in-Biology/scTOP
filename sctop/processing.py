import numpy as np
import pandas as pd
from typing import Optional, Union, List, Tuple
from scipy import special
import scipy.linalg as la
from numba import njit, prange
from scipy import sparse

def process(
    df_in: Union[pd.DataFrame, np.ndarray, sparse.spmatrix],
    average: bool = False,
    chunk_size: Optional[int] = None,
) -> pd.DataFrame:
    """Process scRNA-seq data with optional chunking."""
    if sparse.issparse(df_in):
        arr_in = df_in.toarray()
    elif isinstance(df_in, pd.DataFrame):
        arr_in = df_in.values
        genes = df_in.index
        cols = df_in.columns
    else:
        arr_in = np.asarray(df_in)

    if not isinstance(df_in, pd.DataFrame):
        genes = [f"gene_{i}" for i in range(arr_in.shape[0])]
        cols = [f"sample_{i}" for i in range(arr_in.shape[1])]

    if average:
        arr = arr_in.astype(np.float32, copy=False)
        z = _process_array(arr, average=True)
        return pd.DataFrame(z, index=genes, columns=["avg"])

    n_samples = arr_in.shape[1]
    if chunk_size is None or chunk_size >= n_samples:
        arr = arr_in.astype(np.float32, copy=False)
        z = _process_array(arr, average=False)
        return pd.DataFrame(z, index=genes, columns=cols)
    
    pieces = []
    starts = list(range(0, n_samples, chunk_size))
    
    for start in starts:
        end = min(start + chunk_size, n_samples)
        arr_chunk = arr_in[:, start:end].astype(np.float32, copy=False)
        z_chunk = _process_array(arr_chunk, average=False)
        pieces.append(z_chunk)

    z_all = np.hstack(pieces)
    out_df = pd.DataFrame(z_all, index=genes, columns=cols)
    return out_df

def score(
    basis: pd.DataFrame,
    sample: pd.DataFrame,
    full_output: bool = False,
    chunk_size: Optional[int] = None,
) -> Union[pd.DataFrame, List]:
    """Project sample onto basis with optional chunking."""
    common = np.intersect1d(basis.index.values, sample.index.values)
    if common.size == 0:
        raise ValueError("Basis and sample have no genes in common.")

    basis_sub = basis.loc[common].to_numpy(dtype=np.float32, copy=False)
    if basis_sub.ndim != 2:
        basis_sub = basis_sub.reshape(-1, 1)

    G = float(basis_sub.shape[0])
    A = (basis_sub.T @ basis_sub) / G
    BT = basis_sub.T
    X = la.solve(A, BT, assume_a="sym")
    eta = (X / G).astype(np.float32)

    n_samples = sample.shape[1]
    
    if chunk_size is None or chunk_size >= n_samples:
        sample_sub = sample.loc[common].to_numpy(dtype=np.float32, copy=False)
        if sample_sub.ndim == 1:
            sample_sub = sample_sub.reshape(-1, 1)
        a = eta @ sample_sub
    else:
        starts = list(range(0, n_samples, chunk_size))
        proj_pieces = []
        
        for start in starts:
            end = min(start + chunk_size, n_samples)
            sample_chunk_df = sample.iloc[:, start:end]
            sample_chunk_sub = sample_chunk_df.loc[common].to_numpy(dtype=np.float32, copy=False)
            
            if sample_chunk_sub.ndim == 1:
                sample_chunk_sub = sample_chunk_sub.reshape(-1, 1)
            
            a_chunk = eta @ sample_chunk_sub
            proj_pieces.append(a_chunk)
            
        a = np.hstack(proj_pieces)

    proj_df = pd.DataFrame(a, index=basis.columns, columns=sample.columns)
    if full_output:
        return [proj_df, A, eta]
    return proj_df
    
@njit(parallel=True, fastmath=True)
def _compute_ranks_parallel(arr: np.ndarray) -> np.ndarray:
    """Compute average ranks per column with tie-handling."""
    G, S = arr.shape
    ranks = np.empty((G, S), dtype=np.float32)

    for j in prange(S):
        col = arr[:, j]
        order = np.argsort(col)
        i = 0
        while i < G:
            start = i
            idx0 = order[i]
            val = col[idx0]
            i += 1
            while i < G:
                idxi = order[i]
                if col[idxi] == val:
                    i += 1
                else:
                    break
            end = i - 1
            avg_rank = 0.5 * (start + 1 + end + 1)
            for k in range(start, end + 1):
                ranks[order[k], j] = avg_rank
    return ranks

def rank_zscore_fast(arr: np.ndarray) -> np.ndarray:
    """Rank-based normal scores transformation."""
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    G, S = arr.shape
    if G <= 1:
        return np.zeros_like(arr, dtype=np.float32)
        
    ranks = _compute_ranks_parallel(arr)
    p = ranks.astype(np.float64) / (G + 1.0)
    z = np.sqrt(2.0) * special.erfinv(2.0 * p - 1.0)
    return z.astype(np.float32)

def _process_array(arr: np.ndarray, average: bool = False) -> np.ndarray:
    """Core processing on numpy array."""
    arr = arr.astype(np.float32, copy=False)
    
    col_sums = arr.sum(axis=0)
    zero_mask = col_sums == 0.0
    if zero_mask.any():
        col_sums = col_sums.copy()
        col_sums[zero_mask] = 1.0

    arr = arr / col_sums[np.newaxis, :]

    if average:
        arr = arr.mean(axis=1, keepdims=True)

    arr = np.log2(arr + 1.0).astype(np.float32, copy=False)
    z = rank_zscore_fast(arr)
    z_final = z/np.linalg.norm(z, axis=0, keepdims=True)
    return z_final