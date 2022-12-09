=========================================================================================
scTOP: Single-cell Type Order Parameters
=========================================================================================

A Python package for finding projections onto known cell phenotyes, given matrices of raw RNA count data. 
The theoretical background for this project can be found in `Epigenetic Landscapes Explain Partially Reprogrammed Cells and Identify Key Reprogramming Genes <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003734>`_ by Alex H. Lang, Hu Li, James J. Collins, and Pankaj Mehta. 

A manuscript containing the technical details of applying our method to single-cell RNA-seq will be available in January 2023.

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