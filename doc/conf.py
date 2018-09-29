import sphinx_bootstrap_theme

project = 'Piecewise Polynomial Curves'
author = 'Matthias Geier'
copyright = '2018, ' + author

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'nbsphinx',
]

highlight_language = 'none'
html_sourcelink_suffix = ''

intersphinx_mapping = {'https://docs.python.org/': None}

master_doc = 'index'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']
pygments_style = 'sphinx'

nbsphinx_prolog = r"""
{% set docname = env.doc2path(env.docname, base='doc') %}

.. only:: html

    .. role:: raw-html(raw)
        :format: html

    .. nbinfo::

        This page was generated from `{{ docname }}`__.
        Interactive online version:
        :raw-html:`<a href="https://mybinder.org/v2/gh/AudioSceneDescriptionFormat/splines/{{ env.config.release }}?filepath={{ docname }}"><img alt="Binder badge" src="https://mybinder.org/badge.svg" style="vertical-align:text-bottom"></a>`

    __ https://github.com/AudioSceneDescriptionFormat/splines/blob/
        {{ env.config.release }}/{{ docname }}

.. raw:: latex

    \vfil\penalty-1\vfilneg
    \vspace{\baselineskip}
    \textcolor{gray}{The following section was generated from
    \texttt{\strut{}{{ docname }}}\\[-0.5\baselineskip]
    \noindent\rule{\textwidth}{0.4pt}}
    \vspace{-2\baselineskip}
"""

# This is processed by Jinja2 and inserted after each notebook
nbsphinx_epilog = r"""
.. raw:: latex

    \textcolor{gray}{\noindent\rule{\textwidth}{0.4pt}\\
    \hbox{}\hfill End of
    \texttt{\strut{}{{ env.doc2path(env.docname, base='doc') }}}}
    \vfil\penalty-1\vfilneg
"""
# -- Get version information and date from Git ----------------------------

try:
    from subprocess import check_output
    release = check_output(['git', 'describe', '--tags', '--always'])
    release = release.decode().strip()
    today = check_output(['git', 'show', '-s', '--format=%ad', '--date=short'])
    today = today.decode().strip()
except Exception:
    release = '<unknown>'
    today = '<unknown date>'

# -- Options for HTML output -------------------------------------------------

html_title = project + ', version ' + release
html_theme = 'bootstrap'
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()
html_theme_options = {
    'navbar_title': 'splines',
    'navbar_site_name': 'Pages',
    'navbar_pagenav_name': 'This Page',
    'navbar_fixed_top': True,
    'source_link_position': 'footer',
    #'bootswatch_theme': 'cosmo',
    'bootswatch_theme': 'yeti',
}
html_domain_indices = False
html_show_sourcelink = True

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    'papersize': 'a4paper',
    'printindex': '',
    'preamble': r"""
\usepackage[sc,osf]{mathpazo}
\linespread{1.05}  % see http://www.tug.dk/FontCatalogue/urwpalladio/
\renewcommand{\sfdefault}{pplj}  % Palatino instead of sans serif
\IfFileExists{zlmtt.sty}{
    \usepackage[light,scaled=1.05]{zlmtt}  % light typewriter font from lmodern
}{
    \renewcommand{\ttdefault}{lmtt}  % typewriter font from lmodern
}
""",
}

latex_show_urls = 'footnote'
latex_domain_indices = False

latex_documents = [
    (master_doc, 'PiecewisePolynomialCurves.tex', project, author, 'howto'),
]

# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright
epub_exclude_files = ['search.html']
