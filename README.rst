=========================================================================================
scTOP: Single-cell Type Order Parameters
=========================================================================================

A Python package for finding projections onto known cell phenotyes, given matrices of raw RNA count data. 
A manuscript containing the technical details of applying this method is available `here <https://github.com/Emergent-Behaviors-in-Biology/scTOP-manuscript/tree/main>`_. It has been published in `Development <https://doi.org/10.1242/dev.201873>`_.

Installation
=============

To install using pip: ``pip install sctop``
To install from github: clone or download this repo and import the module from ``sctop/sctop.py``. See ``tutorial/tutorial.ipynb`` for an example.

`(In progress!) Documentation available via Read the Docs <https://sctop.readthedocs.io/>`_ 

Please report bugs or request new features in the Github "Issues" tab. You can also contact me directly at ``mariay@bu.edu``.

Dependencies
-------------
* NumPy
* Pandas
* SciPy

Sources for reference databases
=================================
* `Mouse Cell Atlas <http://bis.zju.edu.cn/MCA/>`_
* `Herriges et al. <https://doi.org/10.1101/2022.07.26.501591>`
* `Atlas of Lung Development <https://doi.org/10.1101/2021.01.21.427641>`

Tutorial on creating your own reference basis using an existing atlas

I suggest making your own reference basis to analyze data that's relevant to your biological context. You can find a few examples of reference basis creation in the [scTOP manuscript repository "notebooks" folder](https://github.com/Emergent-Behaviors-in-Biology/scTOP-manuscript/tree/main/notebooks) -- particularly, you can looks at [create MCA reference basis.ipynb](https://github.com/Emergent-Behaviors-in-Biology/scTOP-manuscript/blob/main/notebooks/create%20MCA%20reference%20basis.ipynb) or [cellbench, human lung atlas.ipynb](https://github.com/Emergent-Behaviors-in-Biology/scTOP-manuscript/blob/main/notebooks/cellbench%2C%20human%20lung%20atlas.ipynb). 

Here is a general version of code I've used to make reference bases. This code assumes the original data source is an h5ad file (for this example, I'm using the LungMAP human lung atlas). ::

    # Load the h5ad file corresponding to the atlas
    atlas = sc.read_h5ad('LungMAP_HumanLung_CellRef.v1.1.h5ad')
    
    # Create a dataframe corresponding to the raw RNA counts
    atlas_df = pd.DataFrame(atlas.X.toarray(), index = atlas.obs.index , columns = atlas.var['_index']).T
    atlas_metadata = atlas.obs
    
    # Choose the metadata column to use as the reference types
    cell_type_column = 'celltype_level3_fullname'
    
    # Count the number of cells per type
    type_counts = atlas_metadata[cell_type_column].value_counts()
    
    # Using fewer than 150-200 cells leads to nonsensical results, due to noise. More cells -> less sampling error 
    threshold = 200 # only use cell types with at least this many cells (but use all cells for training)
    types_above_threshold = type_counts[type_counts > threshold].index
    
    basis_list = []
    training_IDs = []
    
    rng = np.random.default_rng()
    
    for cell_type in tqdm(types_above_threshold):
        cell_IDs = atlas_metadata[atlas_metadata[cell_type_column] == cell_type].index
    #     current_IDs = rng.choice(cell_IDs, size=threshold, replace=False) # This line is for using only the threshold number of cells for the reference basis. This can be useful for testing the accuracy of the basis, but it performs notably worse in accuracy metrics compared to using all possible cells.
        current_IDs = cell_IDs
        
        cell_data = atlas_df[current_IDs]
        training_IDs += [current_IDs] # Keep track of training_IDs so that you can exclude them if you want to test the accuracy
        
        # Average across the cells and process them using the scTOP processing method
        processed = top.process(cell_data, average=True)
        basis_list += [processed]
        
    training_IDs = np.concatenate(training_IDs)
    basis = pd.concat(basis_list, axis=1)
    basis.columns = types_above_threshold

It's also important to validate the reference basis by finding the accuracy. As you can see in the paper, top1 accuracies ranging from 80 to 100% are generally acceptable. Certain less-specified cell types, especially immune cells, may lower the accuracy without affecting the efficacy of identifying other cell types.

Here's a sample code for finding the top1 and top3 accuracies: ::

    test_IDs = np.setdiff1d(atlas_df.columns, training_IDs)
    split_IDs = np.array_split(test_IDs, 10) # I split this test dataset because it's very large and took up a lot of memory -- you don't need to do this if you have enough memory to test the entire dataset at once
    
    # cells with maximum projection under this value are considered "unspecified"
    specification_value = 0.1
    
    accuracies = {'top1': 0,
                  'top3': 0,
                  'unspecified': 0
                 }
    
    predicted_labels = []
    predicted_labels_specified = []
    true_labels = []
    
    for sample_IDs in tqdm(split_IDs):
        test_data = atlas_df[sample_IDs]
        test_processed = top.process(test_data)
        test_projections = top.score(basis, test_processed)
    
        for sample_id, sample_projections in test_projections.iteritems():
            types_sorted_by_projections = sample_projections.sort_values(ascending=False).index
            true_type = atlas_metadata.loc[sample_id, cell_type_column]
    
            true_labels += [true_type]
            top_type = types_sorted_by_projections[0]
            predicted_labels += [top_type]
    
            if sample_projections.max() < specification_value:
                predicted_labels_specified += ['Unspecified']
                accuracies['unspecified'] += 1
            else:
                predicted_labels_specified += [top_type]
    
            if top_type == true_type:
                accuracies['top1'] += 1
            if true_type in types_sorted_by_projections[:3]:
                accuracies['top3'] += 1
                
        del test_data
        del test_processed
        del test_projections

Then you simply check the accuracy like so: ::

    for key, value in accuracies.items():
        print("{}: {}".format(key, value/len(test_IDs)))

For support in applying scTOP, feel free to open an issue ticket on this repository or email me at mariay@bu.edu.
