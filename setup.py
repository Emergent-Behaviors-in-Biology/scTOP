import setuptools

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scTOP",
    version="2.0.2",
    author="Maria Yampolskaya, Huan Souza",
    author_email="mariay@bu.edu, hsouza@bu.edu",
    description="A package for analyzing cell type using single-cell RNA-sequencing data.",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/Emergent-Behaviors-in-Biology/scTOP",
    project_urls={
        "Bug Tracker": "https://github.com/Emergent-Behaviors-in-Biology/scTOP/issues",
        "Docs": "https://sctop.readthedocs.io/"
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    packages=["sctop"],
    license="GPL-3.0-or-later", 
    python_requires=">=3.8",
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "tables",
        "anndata",
        "scikit-learn",
        "matplotlib",
        "seaborn",
        "tqdm",
        "numba",
        "requests"
    ]
)
