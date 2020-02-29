# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'pymatmc2'
copyright = '2020, Eugene J. Ragasa'
author = 'Eugene J. Ragasa'

# The full version, including alpha/beta/rc tags
release = '0.1.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    'sphinx.ext.doctest',
    'sphinx.ext.napoleon',
    'sphinx.ext.imgmath',
    'sphinxcontrib.bibtexwww'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Setting up image math ---------------------------------------------------
default_role = 'math'
pngmath_divpng_args = ['-gamma 1.5','-D 110']
#pngmath_divpng_args = ['-gamma', '1.5', '-D', '110', '-bg', 'Transparent']
imgmath_latex_preamble =  '\\usepackage{amsmath}\n'+\
                          '\\usepackage{mathtools}\n'+\
                          '\\usepackage{amsfonts}\n'+\
                          '\\usepackage{amssymb}\n'+\
                          '\\usepackage{dsfont}\n'+\
                          '\\def\\Z{\\mathbb{Z}}\n'+\
                          '\\def\\R{\\mathbb{R}}\n'+\
                          '\\def\\bX{\\mathbf{X}}\n'+\
                          '\\def\\X{\\mathbf{X}}\n'+\
                          '\\def\\By{\\mathbf{y}}\n'+\
                          '\\def\\Bbeta{\\boldsymbol{\\beta}}\n'+\
                          '\\def\\U{\\mathbf{U}}\n'+\
                          '\\def\\V{\\mathbf{V}}\n'+\
                          '\\def\\V1{\\mathds{1}}\n'+\
                          '\\def\\hU{\\mathbf{\hat{U}}}\n'+\
                          '\\def\\hS{\\mathbf{\hat{\Sigma}}}\n'+\
                          '\\def\\hV{\\mathbf{\hat{V}}}\n'+\
                          '\\def\\E{\\mathbf{E}}\n'+\
                          '\\def\\F{\\mathbf{F}}\n'+\
                          '\\def\\x{\\mathbf{x}}\n'+\
                          '\\def\\h{\\mathbf{h}}\n'+\
                          '\\def\\v{\\mathbf{v}}\n'+\
                          '\\def\\nv{\\mathbf{v^{{\bf -}}}}\n'+\
                          '\\def\\nh{\\mathbf{h^{{\bf -}}}}\n'+\
                          '\\def\\s{\\mathbf{s}}\n'+\
                          '\\def\\b{\\mathbf{b}}\n'+\
                          '\\def\\c{\\mathbf{c}}\n'+\
                          '\\def\\W{\\mathbf{W}}\n'+\
                          '\\def\\C{\\mathbf{C}}\n'+\
                          '\\def\\P{\\mathbf{P}}\n'+\
                          '\\def\\T{{\\bf \\mathcal T}}\n'+\
                          '\\def\\B{{\\bf \\mathcal B}}\n'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']