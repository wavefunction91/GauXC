Collocation Example
===================

Single Basis Functions
----------------------

A collocation grid between a single basis and a Cartesian grid can be computed
with the :c:func:`~gg_collocation` function. For example, we will use a grid
starting at the origin along the ``z`` axis and a ``S`` shell at the origin:

.. code-block:: C

  #include <stdio.h>
  #include "gau2grid.h"

  int main() {
      // Generate grid
      long int npoints = 5;
      double xyz[15] = {0, 0, 0, 0, 0, // x components
                        0, 0, 0, 0, 0}; // y components
                        0, 1, 2, 3, 4}; // z components
      long int xyz_stride = 1; // This is a contiguous format

      // Gaussian data
      int nprim = 1;
      double coef[1] = {1};
      double exp[1] = {1};
      double center[3] = {0, 0, 0};
      int order = GG_CARTESIAN_CCA; // Use cartesian components

      double s_output[5] = {0};
      gg_collocation(0,                                // The angular momentum
                     npoints, xyz, xyz_stride,          // Grid data
                     nprim, coef, exp, center, order,  // Gaussian data
                     s_output);                        // Output

      // Print output to stdout
      for (int i = 0; i < npoints; i += 1) {
          printf("%lf  ", s_output[i]);
      }
      printf("\n");
  }

The resulting output should be:

.. code-block:: bash

  1.000000  0.367879  0.018316  0.000123  0.000000

For higher angular momentum functions that output size should ``ncomponents x
npoints`` in size. Where each component is on a unique row or the ``X``
component starts at position ``0``, the ``Y`` component starts at position
``5``, and the ``Z`` component starts at position ``10`` as out grid is of
length ``5``. See :ref:`Gaussian Component Orders <gpo_order>` for more details or order output.

The xyz input shape can either be organized contiguously in each dimension like
the above or packed in a xyz, xyz, ... fashion. If the ``xyz_stride`` is not 1,
the shape refers to the strides per row. For example, if the data is packed as
xyzw, xyzw, ... (where w could be a DFT grid weight) the ``xyz_stride`` should
be 4.

.. code-block:: C

      long int xyz_stride = 3;
      double xyz[15] = {0, 0, 0,
                        0, 0, 1,
                        0, 0, 2,
                        0, 0, 3,
                        0, 0, 4}; // xyz, xyz, ... format


      gg_collocation(0,                                // The angular momentum
                     npoints, xyz, xyz_stride,          // Grid data
                     nprim, coef, exp, center, order,  // Gaussian data
                     s_output);                        // Output

Multiple Basis Functions
------------------------

Often collocation matrices are computed for multiple basis functions at once.
The below is an example of usage:

.. code-block:: C

  #include <stdio.h>
  #include "gau2grid.h"

  int main() {
      // Generate grid
      long int npoints = 5;
      double xyz[15] = {0, 0, 0, 0, 0, // x components
                        0, 0, 0, 0, 0}; // y components
                        0, 1, 2, 3, 4}; // z components
      long int xyz_stride = 1;

      // Gaussian data
      int nprim = 1;
      double coef[1] = {1};
      double exp[1] = {1};
      double center[3] = {0, 0, 0};
      int order = GG_SPHERICAL_CCA; // Use cartesian components

      // Size ncomponents * npoints, (1 + 3 + 5) * 5
      double output[45] = {0};
      int row = 0;
      for (int L = 0; L < 3; L++) {
          gg_collocation(L,                                 // The angular momentum
                         npoints, xyz, xyz_stride            // Grid data
                         nprim, coef, exp, center, order,   // Gaussian data
                         output + (row * npoints));         // Output, shift pointer

          row += gg_ncomponents(L, spherical); // Increment rows skipped
      }

      // Print out by row
      for (int i = 0; i < row; i += 1) {
          for (int j = 0; j < npoints; j += 1) {
              printf("%lf  ", output[i * npoints + j]);
          }
          printf("\n");
      }
  }

The resulting output should be:

.. code-block:: bash

  1.000000  0.367879  0.018316  0.000123  0.000000 // S
  0.000000  0.367879  0.036631  0.000370  0.000000 // P_0
  0.000000  0.000000  0.000000  0.000000  0.000000 // P^+_0
  0.000000  0.000000  0.000000  0.000000  0.000000 // P^-_0
  0.000000  0.367879  0.073263  0.001111  0.000002 // D_0
  0.000000  0.000000  0.000000  0.000000  0.000000 // D^+_1
  0.000000  0.000000  0.000000  0.000000  0.000000 // D^-_1
  0.000000  0.000000  0.000000  0.000000  0.000000 // D^+_2
  0.000000  0.000000  0.000000  0.000000  0.000000 // D^-_2
