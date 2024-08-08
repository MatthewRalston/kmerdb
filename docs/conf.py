#!/bin/env python3

import sys
import os

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'kmerdb'

release = '0.8.6'
author = 'Matt Ralston <mralston.development@gmail.com>, and the kmerdb project contributors'
copyright = '2020, Matt Ralston <mralston.development@gmail.com>'
language = 'en'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']


#html_theme = 'sphinx_rtd_theme'
#html_context = {}





extensions = [
    #'sphinx.ext.duration',
    #'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
#    'sphinx_rtd_theme',
]

templates_path = ['_templates']
source_suffix = '.rst'
exclude_patterns = []
locale_dirs = ['locale/']
gettext_compact = False

master_doc = 'index'
suppress_warnings = ['image.nonlocal_uri']
pygments_style = 'default'

if sys.version_info < (3, 0):
    tags.add("python2")
else:
    tags.add("python3")

intersphinx_mapping = {
    'numpy': ('https://numpy.readthedocs.io/en/latest', None),
    'rtd': ('https://docs.readthedocs.io/en/stable/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}


