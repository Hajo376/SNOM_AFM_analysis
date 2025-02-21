# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# import os
# import sys
# sys.path.insert(0, os.path.abspath('../src'))
import sys
from pathlib import Path

sys.path.insert(0, str(Path('..', 'src').resolve()))

project = 'SNOM Analysis'
copyright = '2025, Hans-Joachim Schill'
author = 'Hans-Joachim Schill'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc', 
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    "sphinx_design",
    # "sphinx_togglebutton",
    "numpydoc",
    'sphinx.ext.napoleon', 
    'sphinx.ext.viewcode', 
    # 'sphinx.ext.graphviz',
    # 'sphinx.ext.inheritance_diagram',
    # 'sphinx.ext.ifconfig',
    # 'sphinx.ext.autosectionlabel',
    ]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autosectionlabel_prefix_document = True
autosectionlabel_maxdepth = 4

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']


# -- Options for python code -------------------------------------------------

pygments_style = 'sphinx'

# -- Options for autodoc -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#confval-autodoc_member_order
autodoc_member_order = 'bysource'