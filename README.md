# SampleDisc documentation site

This directory now contains the MkDocs source for the SampleDisc documentation website.

## Local preview

```bash
cd /users/hjiang/GenoDistance/website/SampleDisco_tutorial
python -m venv .venv
source .venv/bin/activate
pip install -r requirements-docs.txt
mkdocs serve
```

Then open `http://127.0.0.1:8000`.

## Build static files

```bash
mkdocs build
```

The generated site is written to `site/`.

## Deploy on GitHub Pages

A GitHub Actions workflow is included at `.github/workflows/deploy-docs.yml`.
Push the repository to GitHub, enable Pages with **GitHub Actions** as the source,
and each push to `main` or `master` will rebuild the docs site.
