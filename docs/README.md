# Documentation

The primary documentation site is now built with Sphinx from `docs/source/`.

## Build Locally

```bash
python -m pip install -r docs/requirements.txt
python docs/scripts/render_schema_docs.py
sphinx-build -b html docs/source docs/build/html
```

Open `docs/build/html/index.html` after the build completes.

## Legacy Guides

These Markdown guides are still useful while the Sphinx site is filled out:

- [QUICK_START.md](QUICK_START.md)
- [INDEX.md](INDEX.md)
- [MODULES.md](MODULES.md)
- [PATTERNS.md](PATTERNS.md)
- [workflows.md](workflows.md)
- [CONFIGURATION.md](CONFIGURATION.md)
