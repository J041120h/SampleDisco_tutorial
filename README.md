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
Push the repository to GitHub, then use one of these two deployment modes:

### Preferred: GitHub Actions

1. Open the repository on GitHub.
2. Go to `Settings` -> `Pages`.
3. Under **Build and deployment**, set **Source** to `GitHub Actions`.
4. Push to `main` or `master`.

The included workflow will build MkDocs and publish the generated HTML, so the
documentation appears directly at your GitHub Pages URL.

### Fallback: Branch deployment

If you prefer **Deploy from a branch**:

1. Run `mkdocs build`.
2. Commit the generated `site/` folder.
3. Keep the repository root `index.html` file.

That root page now automatically redirects to `./site/` when the built docs are
present, so GitHub Pages can still display the generated documentation.
