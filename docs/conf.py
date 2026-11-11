# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import re

project = 'discoal'
copyright = '2024-2026, Andrew D. Kern'
author = 'Andrew D. Kern'

# The version comes from src/core/version.h, the one place that sets it.
with open(os.path.join(os.path.dirname(__file__), '..', 'src', 'core',
                       'version.h')) as _f:
    release = re.search(r'#define DISCOAL_VERSION_STRING "([^"]+)"',
                        _f.read()).group(1)
version = '.'.join(release.split('.')[:2])

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
    'myst_parser',
]

# The YAML and demes pages are Markdown, so that they also read well on
# GitHub; myst_parser builds them with the rest.
myst_heading_anchors = 3

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Options for LaTeX output ------------------------------------------------
latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
}