# Building the documentation

Install the documentation dependencies and the package from the repository root:

```bash
python -m pip install -e "./python[docs]"
```

Build the JupyterBook site:

```bash
jupyter-book build docs -vnW --keep-going
```

The generated site is in `docs/_build/html/`. The GitHub Actions workflow builds
the book for pull requests and publishes it to GitHub Pages when changes are
pushed to `main`.
