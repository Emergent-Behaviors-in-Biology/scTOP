"""
scTOP: Single-Cell Type Order Parameters
-----

    A Python module for calculating cell phenotype order parameters from single-cell RNA sequencing data.

    @author: Maria Yampolskaya (with modified code from R. A. Marsland)
    
    Copyright (C) 2022  Maria Yampolskaya, Robert A. Marsland III

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

import numpy as np
import pandas as pd
import scipy.stats as sps

def load_basis(basis_name, minimum_cells=100, source_dir=__file__[:-9] + '/databases/'):
    """
    Return the desired basis and corresponding metadata.
    
    
    Parameters
    ----------
    basis_name : {'MC20-K022', 'LD21'}
        'MC20' : Mouse Cell Atlas
        'LD21' : Atlas of Lung Development
    minimum_cells : int, optional
        Only include cell types in the basis with a number of cells greater than or
        equal to ``minimum_cells``.
    source_dir : string, optional
        Path to where data_``basis_name``.h5 and metadata_``basis_name``.csv are located.
        
    Returns
    -------
    out : tuple of DataFrames
        Tuple of the form ``(basis, metadata)``.
        basis : (G, C) DataFrame
            Table with ``G`` genes and ``C`` cell types containing processed scRNA-seq data.
        metadata : (C, T) DataFrame
            Table with ``C`` cell types and ``T`` tags containing metadata information.
            
    """
    
    acceptable_bases = ['MC-KO', 'MC', 'LD']
    if basis_name not in acceptable_bases:
        raise ValueError('Provided basis_name must refer to one of the existing bases: {}'
                         .format(acceptable_bases))
    
    if basis_name == 'MC':
        basis_name = 'MC20_apr22'
    elif basis_name == 'MC-KO':
        basis_name = 'MCKO_apr22'
    elif basis_name == 'LD':
        basis_name = 'LD_apr22'
    
    data = pd.read_hdf(source_dir + 'data_' + basis_name + '.h5')
    metadata = pd.read_csv(source_dir + 'metadata_' + basis_name + '.csv', index_col=0)
    
    # if Mouse Cell Atlas, drop cultured mesenchymal and trophoblast stem cells
    # ALso drop embroynic cells, because we only want an adult basis
    if 'MC' in basis_name:
        exceptions = [exception for exception in data.columns 
                      if ('Trophoblast' in exception) or ('Cultured' in exception) or ('E14.5' in exception) or ('Embryonic' in exception)]
        data = data.drop(columns = exceptions)
        metadata = metadata.drop(index = exceptions)
        
    types_above_minimum = metadata['Cell Count'] > minimum_cells
    
    return (data.loc[:, types_above_minimum].dropna(axis=0), metadata.loc[types_above_minimum])

def score(basis, sample, full_output=False):
    """
    Return the cell type order parameters of the given sample projected onto the given basis.
    Optionally return correlation matrix and predictivity.
    
    
    Parameters
    ----------
    basis : (G_B, C) DataFrame
        Table with ``G_B`` genes and ``C`` cell types containing processed scRNA-seq data.
        ``basis.index``` should be a list of gene names.
    sample : (G_D, S) DataFrame
        Table with ``G_D`` genes and ``S`` samples of processed scRNA-seq data.
        ``sample.index`` should be a list of gene names.
    full_output : bool, optional
        if False: only return projections, ``a``
        if True: return array containing projections, correlation matrix, and predictivity, ``[a, A, eta]``
    
    Returns
    -------
    out : DataFrame (or list of DataFrame and arrays, if full_output)
        Projections of the sample on the given basis.
        If ``full_output``, return a list containing the projections, correlation matrix, and predictivity.
        
    """

    common_genes = np.intersect1d(basis.index, sample.index)
    
    if len(common_genes) == 0:
        raise ValueError('Basis and sample have no genes in common.')
        
    basis_values = basis.loc[common_genes].values
    sample_values = sample.loc[common_genes].values
    
    if len(sample_values.shape) == 1:
        sample_values = sample_values.reshape((sample_values.size,1))

    A = np.dot(basis_values.T, basis_values) / basis_values.shape[0]
    eta = np.linalg.solve(A, basis_values.T) / basis_values.shape[0]
    a = pd.DataFrame(np.dot(eta, sample_values), index=basis.columns, columns=sample.columns)

    if full_output:
        return [a, A, eta]
    else:
        return a
    
def rank_zscore(data):
    """ 
    Return the normally-distributed z-scores of the given data. Each column is individually 
    ranked against itself, then the ranked data is z-scored assuming a normal distribution.
    
    
    Parameters
    ----------
    data : (G, S) array
        Array with ``G`` genes and ``S`` independent samples containing scRNA-seq data.
    
    Returns
    -------
    out : (G, S) array
        Array with ``G`` genes and ``S`` independent samples containing ranked then z-scored scRNA-seq data.
    
    """

    if len(data.shape) == 1:
        data = data.reshape((data.size,1))
        
    G, S = data.shape
    
    # Calculate rank for each sample
    rank_data = np.zeros(data.shape) 
    
    for i in range(S):
        # by default, ties method is average
        rank_data[:,i] = sps.rankdata(data[:,i]) 

    # Convert rank_data to percentiles according to a normal distribution
    output_data = sps.norm.ppf(rank_data / (G+1))

    return output_data

def process(data, average=False):
    """
    Return the processed scRNA-seq data. Each sample is independently normalized, then the fitted z-score
    is calculated. Optionally, an average is taken across samples before z-scoring.
    
    Parameters
    ----------
    data : (G, S) DataFrame
        Table with ``G`` genes and ``S`` independent samples containing raw counts from scRNA-seq.
    average : bool, optional
        if True: take the average across samples after normalizing and before z-scoring
        
    Returns
    -------
    out : (G, S) DataFrame (or (G, 1) if ``average``)
        Table with ``G`` genes and ``S`` independent samples containing processed scRNA-seq data.
    
    """
    
    # Normalize each sample independently
    data_normalized = data/np.sum(data,axis=0)
    
    # Average across samples, if requested
    if average:
        data_normalized = data_normalized.mean(axis=1)
    
    # Find the normal-distribution z-scores of log(expression + 1)
    data_zscored = rank_zscore(np.log2(data_normalized.values + 1))
    
    if average:
        return pd.DataFrame(data_zscored, index=data.index)
    else:
        return pd.DataFrame(data_zscored, index=data.index, columns=data.columns)