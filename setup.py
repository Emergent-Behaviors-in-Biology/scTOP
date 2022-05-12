import setuptools

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scTOP",
    version="0.0.4",
    author="Maria Yampolskaya",
    author_email="mariay@bu.edu",
    description="A package for identifying cell phenotype from single-cell RNA-sequencing data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Emergent-Behaviors-in-Biology/scTOP",
    project_urls={
        "Bug Tracker": "https://github.com/Emergent-Behaviors-in-Biology/scTOP/issues",
        "Docs": "https://sctop.readthedocs.io/"
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "sctop"},
    packages=setuptools.find_packages(where="sctop"),
    python_requires=">=3.6",
    include_package_data=True,
)
