"""
Configuration file for the Sphinx documentation builder.
"""

import datetime
import os
from pathlib import Path

from sunpy_sphinx_theme import PNG_ICON

from irispy import __version__

# -- Project information -----------------------------------------------------
project = "irispy-lmsal"
author = "IRIS Instrument Team"
copyright = f"{datetime.datetime.now(datetime.timezone.utc).year}, {author}"  # NOQA: A001
release = __version__
is_development = ".dev" in __version__

# -- General configuration ---------------------------------------------------
extensions = [
    "hoverxref.extension",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_automodapi.automodapi",
    "sphinx_automodapi.smart_resolver",
    "sphinx_changelog",
    "sphinx_gallery.gen_gallery",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
]
automodapi_toctreedirnm = "generated/api"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = ".rst"
master_doc = "index"
default_role = "obj"

# -- Options for hoverxref -----------------------------------------------------
if os.environ.get("READTHEDOCS"):
    hoverxref_api_host = "https://readthedocs.org"
    if os.environ.get("PROXIED_API_ENDPOINT"):
        # Use the proxied API endpoint
        # - A RTD thing to avoid a CSRF block when docs are using a
        #   custom domain
        hoverxref_api_host = "/_"

hoverxref_tooltip_maxwidth = 600  # RTD main window is 696px
hoverxref_auto_ref = True
hoverxref_mathjax = True
# hoverxref has to be applied to these
hoverxref_domains = ["py"]
hoverxref_role_types = {
    # roles with py domain
    "attr": "tooltip",
    "class": "tooltip",
    "const": "tooltip",
    "data": "tooltip",
    "exc": "tooltip",
    "func": "tooltip",
    "meth": "tooltip",
    "mod": "tooltip",
    "obj": "tooltip",
    # roles with std domain
    "confval": "tooltip",
    "hoverxref": "tooltip",
    "ref": "tooltip",
    "term": "tooltip",
}

# -- Options for sphinx-copybutton ---------------------------------------------
# Python Repl + continuation, Bash, ipython and qtconsole + continuation, jupyter-console + continuation
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True

# -- Options for intersphinx extension ---------------------------------------
intersphinx_mapping = {
    "python": (
        "https://docs.python.org/3/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/python3.inv"),
    ),
    "numpy": (
        "https://numpy.org/doc/stable/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/numpy.inv"),
    ),
    "scipy": (
        "https://docs.scipy.org/doc/scipy/reference/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/scipy.inv"),
    ),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "aiapy": ("https://aiapy.readthedocs.io/en/stable/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "astroquery": ("https://astroquery.readthedocs.io/en/latest/", None),
    "drms": ("https://docs.sunpy.org/projects/drms/en/stable/", None),
    "mpl_animators": ("https://docs.sunpy.org/projects/mpl-animators/en/stable/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "parfive": ("https://parfive.readthedocs.io/en/stable/", None),
    "reproject": ("https://reproject.readthedocs.io/en/stable/", None),
    "skimage": ("https://scikit-image.org/docs/stable/", None),
    "sqlalchemy": ("https://docs.sqlalchemy.org/en/latest/", None),
    "sunkit_image": ("https://docs.sunpy.org/projects/sunkit-image/en/stable/", None),
    "sunkit_instruments": (
        "https://docs.sunpy.org/projects/sunkit-instruments/en/stable/",
        None,
    ),
    "asdf": ("https://asdf.readthedocs.io/en/stable/", None),
    "sunpy": ("https://docs.sunpy.org/en/stable/", None),
    "sunraster": ("https://docs.sunpy.org/projects/sunraster/en/stable/", None),
    "zeep": ("https://docs.python-zeep.org/en/stable/", None),
}

# -- Options for HTML output -------------------------------------------------
html_theme = "sunpy"
graphviz_output_format = "svg"
sphinx_gallery_conf = {
    "backreferences_dir": Path("generated") / "modules",
    "filename_pattern": "^((?!skip_).)*$",
    "examples_dirs": Path("..") / "examples",
    "within_subsection_order": "ExampleTitleSortKey",
    "gallery_dirs": Path("generated") / "gallery",
    "default_thumb_file": PNG_ICON,
    "abort_on_example_error": False,
    "plot_gallery": "True",
    "remove_config_comments": True,
    "doc_module": ("sunpy"),
    "only_warn_on_example_error": True,
    "matplotlib_animations": True,
}
