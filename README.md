## Sample Disco Tutorial Website

This repository contains the source for the **Sample Disco** tutorial website, used as the project page for the GenoDistance / Sample Disco method on GitHub.

The site is meant to:
- **Introduce** the Sample Disco method at a high level
- **Show** small, concrete examples (toy datasets, basic commands)
- **Link** to the main paper, code, and datasets

### Live demo (example)

If you deploy this with GitHub Pages, your site URL will look like:

`https://<your-github-username>.github.io/SampleDisco_tutorial/`

Update this line once you have the real URL:

- **Live site**: _coming soon_  

### Repository structure

- `index.html`: **Landing page** with a short method overview and links
- `assets/` or `static/`: **CSS, JS, and images** for the website (create as needed)
- `examples/`: **Small tutorial examples** (e.g., example input/output, small figures)

You can adapt this structure to match your actual files.

### How to view locally

From the `SampleDisco_tutorial` folder:

```bash
# Option 1: Python simple server
python -m http.server 8000

# Then open in a browser:
# http://localhost:8000
```

Or open `index.html` directly in a browser (right-click → "Open With" → your browser).

### How to publish on GitHub Pages

1. **Push** this folder as a GitHub repository (or as a subdirectory if you prefer).
2. In your GitHub repo, go to **Settings → Pages**.
3. Under **Source**, select the branch (e.g., `main`) and the root folder (or `/docs` if you move the site there).
4. Click **Save**. After a few minutes, GitHub will give you a public URL.

### Customization ideas

- Add a **figure** that illustrates the Sample Disco pipeline.
- Include a **"Quick Start"** section with 3–5 minimal commands.
- Add a **"Cite us"** section with BibTeX for the related paper.
- Embed small **interactive plots** or demo scripts if relevant.

Feel free to edit any of this text to better match your actual project.