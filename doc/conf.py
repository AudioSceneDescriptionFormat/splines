project = 'Splines in Euclidean Space and Beyond'
author = 'Matthias Geier'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'nbsphinx',
    'sphinxcontrib.bibtex',
    'sphinx_last_updated_by_git',
]

bibtex_bibfiles = ['references.bib']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/', None),
}

autosummary_generate = ['python-module']
autoclass_content = 'init'
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}

default_role = 'any'

linkcheck_ignore = [
    # Anchors with line numbers don't seem to work with linkcheck builder
    'https://github.com/scipy/scipy/blob/',
]

nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'svg', 'pdf'}",
]


def if_docname(text):
    return r"""
{% if not env.docname.endswith('index') and env.docname not in [
    'euclidean/hermite',
    'euclidean/natural',
    'euclidean/bezier',
    'euclidean/catmull-rom',
    'euclidean/kochanek-bartels',
    'euclidean/end-conditions',
] %}
{% set docname = 'doc/' + env.doc2path(env.docname, base=None) %}
{% set latex_href = ''.join([
    '\href{https://github.com/AudioSceneDescriptionFormat/splines/blob/',
    env.config.release,
    '/',
    docname | escape_latex,
    '}{\sphinxcode{\sphinxupquote{',
    docname | escape_latex,
    '}}}',
]) %}
""" + text + r"""
{% endif %}
"""


nbsphinx_prolog = if_docname(r"""
.. raw:: html

    <div class="admonition note">
      This page was generated from
      <a class="reference external" href="https://github.com/AudioSceneDescriptionFormat/splines/blob/{{ env.config.release|e }}/{{ docname|e }}">{{ docname|e }}</a>.
      Interactive online version:
      <span style="white-space: nowrap;"><a href="https://mybinder.org/v2/gh/AudioSceneDescriptionFormat/splines/{{ env.config.release|e }}?filepath={{ docname|e }}"><img alt="Binder badge" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a>.</span>
    </div>

.. raw:: latex

    \nbsphinxstartnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{The following section was generated from {{ latex_href }}
    \dotfill}}
""")

nbsphinx_epilog = if_docname(r"""
.. raw:: latex

    \nbsphinxstopnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{\dotfill\ {{ latex_href }} ends here.}}
""")

# -- Work-around to get LaTeX References at the same place as HTML --------

# See https://github.com/mcmtroffaes/sphinxcontrib-bibtex/issues/156

import docutils
import sphinx.builders.latex

class DummyTransform(docutils.transforms.Transform):

    default_priority = 0

    def apply(self):
        pass

sphinx.builders.latex.BibliographyTransform = DummyTransform

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

html_title = 'splines, version ' + release
html_theme = 'insipid'
html_domain_indices = False
html_favicon = 'favicon.svg'
html_copy_source = False
html_add_permalinks = '\N{SECTION SIGN}'
html_show_copyright = False

mathjax_config = {
    'TeX': {'equationNumbers': {'autoNumber': 'AMS', 'useLabelIds': True}},
}

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    'papersize': 'a4paper',
    'printindex': '',
    'sphinxsetup': r"""
        VerbatimColor={HTML}{F5F5F5},
        VerbatimBorderColor={HTML}{E0E0E0},
        noteBorderColor={HTML}{E0E0E0},
        noteborder=1.5pt,
        warningBorderColor={HTML}{E0E0E0},
        warningborder=1.5pt,
        warningBgColor={HTML}{FBFBFB},
    """,
    'preamble': r"""
\usepackage[sc,osf]{mathpazo}
\linespread{1.05}  % see http://www.tug.dk/FontCatalogue/urwpalladio/
\renewcommand{\sfdefault}{pplj}  % Palatino instead of sans serif
\IfFileExists{zlmtt.sty}{
    \usepackage[light,scaled=1.05]{zlmtt}  % light typewriter font from lmodern
}{
    \renewcommand{\ttdefault}{lmtt}  % typewriter font from lmodern
}
\usepackage{mathrsfs}  % for \mathscr{}
""",
}

latex_show_urls = 'footnote'
latex_show_pagerefs = True
latex_domain_indices = False

latex_documents = [
    ('index', 'splines.tex', project, author, 'howto'),
]

# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = ''
epub_exclude_files = ['search.html']
