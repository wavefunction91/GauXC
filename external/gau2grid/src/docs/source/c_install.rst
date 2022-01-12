C installation
==============

You can install gau2grid with ``conda`` or by installing from source.

Conda
-----

You can update gau2grid using `conda <https://www.anaconda.com/download/>`_::

    conda install gau2grid -c psi4

This installs the gau2grid library.


Install from Source
-------------------

Gau2grid uses the CMake build system to compile and configure options. To begin, clone the repository:

.. code-block:: bash

    git clone https://github.com/dgasmith/gau2grid.git
    cd gau2grid

A basic CMake build can then be executed with:

.. code-block:: bash

    cmake -H. -Bobjdir
    cd objdir
    make
    make install

CMake Options
-------------
Gau2grid can be compiled with the following CMake options:

 - ``CMAKE_INSTALL_PREFIX`` - The path to install the library to (default, ``/usr/local``)
 - ``CMAKE_INSTALL_LIBDIR`` - Directory to which libraries installed
 - ``MAX_AM`` - The maximum gaussian angular momentum to compile (default, ``8``)
 - ``CMAKE_BUILD_TYPE`` - Build type (Release or Debug) (default, ``Release``)
 - ``ENABLE_XHOST`` - Enables processor-specific optimization (default, ``ON``)
 - ``BUILD_FPIC`` - Libraries will be compiled with position independent code (default, ``ON``)
 - ``BUILD_SHARED_LIBS`` - Build final library as shared, not static (default, ``ON``)
 - ``ENABLE_GENERIC`` - Enables mostly static linking of system libraries for shared library (default, ``OFF``)

CMake options should be prefixed with ``-D``, for example:

.. code-block:: bash

    cmake -H. -Bobjdir -DCMAKE_INSTALL_PREFIX=~/installs
