project = "GauXC"
version = "1.1"
release = version


extensions = [
    "breathe",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.extlinks",
    "sphinxcontrib.bibtex",
]

breathe_default_project = project
breathe_projects = {
    project: "_doxygen/xml",
}
breathe_domain_by_extension = {
    "h": "c",
}
breathe_show_include = True

bibtex_bibfiles = ["_static/bib/references.bib"]

cpp_id_attributes = [
    "HOST_DEVICE_ACCESSIBLE",
    "GAUXC_MPI_CODE",
]

templates_path = ["_templates"]
source_suffix = [".rst"]
master_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
pygments_style = "default"


html_theme = "sphinx_book_theme"
html_title = project

html_theme_options = {
    "repository_url": "https://github.com/wavefunction91/gauxc",
    "repository_branch": "master",
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_download_button": False,
    "path_to_docs": "docs",
}

html_static_path = ["_static"]


def run_doxygen(folder):
    """Run the doxygen make command in the designated folder"""

    import subprocess, sys

    try:
        retcode = subprocess.call("set -ex; cd %s; doxygen Doxyfile" % folder, shell=True)
        if retcode < 0:
            sys.stderr.write("doxygen terminated by signal %s" % (-retcode))
    except OSError as e:
        sys.stderr.write("doxygen execution failed: %s" % e)


def generate_doxygen_xml(app):
    """Run the doxygen make commands if we're on the ReadTheDocs server"""

    import os, os.path

    read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

    if read_the_docs_build:

        _dir = os.path.dirname(os.path.abspath(__file__))
        run_doxygen(_dir)


def setup(app):

    # Add hook for building doxygen xml when needed
    app.connect("builder-inited", generate_doxygen_xml)
