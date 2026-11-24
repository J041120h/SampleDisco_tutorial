# SampleDisco documentation site

MkDocs source for the SampleDisco documentation: <https://j041120h.github.io/SampleDisco_tutorial/>.
The package itself lives at <https://github.com/J041120h/SampleDisco>.

## Local preview

```bash
git clone https://github.com/J041120h/SampleDisco_tutorial.git
cd SampleDisco_tutorial
python -m venv .venv && source .venv/bin/activate
pip install -r requirements-docs.txt
mkdocs serve          # http://127.0.0.1:8000
```

`mkdocs build --strict` checks every link and anchor; the build goes to `site/` (not committed).

## Deployment

Every push to `main` builds the site with `mkdocs build` and publishes it with
GitHub Actions (`.github/workflows/deploy-docs.yml`). The repository's Pages
source must be set to **GitHub Actions** (`Settings` → `Pages`).
