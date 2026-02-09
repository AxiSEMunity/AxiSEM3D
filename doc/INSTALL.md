For build and environment information see:
https://www.sphinx-doc.org/en/master/tutorial/getting-started.html


To install:

    python -m venv .venv
    source .venv/bin/activate
    python -m pip install sphinx

You may need to also install:

    pip install sphinxcontrib-bibtex
    pip install sphinx-rtd-theme
    pip install sphinx-book-theme
    pip install myst-parser

verify environment:

     sphinx-build --version

Build:

    sphinx-build -M html docs/source/ docs/build/
    