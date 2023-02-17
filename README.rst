=========================================================================================
scTOP: Single-cell Type Order Parameters
=========================================================================================

A Python package for finding projections onto known cell phenotyes, given matrices of raw RNA count data. 
A manuscript containing the technical details of applying this method is available `here <https://github.com/Emergent-Behaviors-in-Biology/scTOP-manuscript/tree/main>`_.

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
