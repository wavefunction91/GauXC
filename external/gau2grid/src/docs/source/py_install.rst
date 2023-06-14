Python installation
===================

You can install gau2grid with ``conda`` or by installing from source.

Conda
-----

You can update gau2grid using `conda <https://www.anaconda.com/download/>`_::

    conda install pygau2grid -c psi4

This installs gau2grid and the NumPy dependancy.


Install from Source
-------------------

To install gau2grid from source, clone the repository from `github
<https://github.com/dgasmith/gau2grid>`_::

    git clone https://github.com/dgasmith/gau2grid.git
    cd gau2grid
    python setup.py install


Test
----

Test gau2grid with ``py.test``::

    cd gau2grid
    py.test
