# Building the Documentation

This directory contains the Sphinx documentation for scTOP.

## Prerequisites

Install the documentation dependencies:

```bash
pip install -r requirements.txt
```

## Building HTML Documentation

From the `docs/` directory, run:

```bash
make html
```

Or on Windows:

```bash
make.bat html
```

The built documentation will be in `_build/html/`. Open `_build/html/index.html` in your browser.

## Building for ReadTheDocs

The documentation is configured to build automatically on ReadTheDocs when pushed to GitHub. The configuration is in:

- `docs/conf.py` - Sphinx configuration
- `docs/requirements.txt` - Python dependencies for building

## Documentation Structure

```
docs/
├── index.rst              # Main landing page
├── installation.rst       # Installation instructions
├── quickstart.rst         # Quick start guide
├── tutorials.rst          # Detailed tutorials
├── theory.rst            # Mathematical theory
├── api/                  # API Reference
│   ├── index.rst         # API overview
│   ├── core.rst          # Core functions (process, score)
│   ├── basis.rst         # Basis management
│   ├── analysis.rst      # Analysis utilities
│   └── visualization.rst # Plotting functions
└── requirements.txt      # Dependencies for building
```

## Clean Build

To start fresh:

```bash
make clean
make html
```

## Viewing Locally

After building, you can serve the documentation locally:

```bash
cd _build/html
python -m http.server 8000
```

Then open http://localhost:8000 in your browser.

## Notes

- The API documentation is auto-generated from docstrings using `sphinx-autoapi`
- Mathematical formulas use MathJax for rendering
- The theme is ReadTheDocs (sphinx_rtd_theme)
- Documentation follows NumPy/Google docstring conventions
