#pragma once
/*
  Copyright (c) 2015-2017, Norbert Juffa
  All rights reserved.

  Redistribution and use in source and binary forms, with or without 
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright 
     notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* Compute exponential function. maximum ulp error observed = 0.89028 */
__device__ __noinline__ double my_exp (double a)
{
    const double ln2_hi = 6.9314718055829871e-01;
    const double ln2_lo = 1.6465949582897082e-12;
    const double l2e = 1.4426950408889634; // log2(e)
    const double cvt = 6755399441055744.0; // 3 * 2**51
    double f, j, p, r;
    int i;

    // exp(a) = 2**i * exp(f); i = rint (a / log(2))
    j = fma (l2e, a, cvt);
    i = __double2loint (j);
    j = j - cvt;
    f = fma (j, -ln2_hi, a);
    f = fma (j, -ln2_lo, f);

    // approximate p = exp(f) on interval [-log(2)/2, +log(2)/2]
    p =            2.5022018235176802e-8;
    p = fma (p, f, 2.7630903481118922e-7);
    p = fma (p, f, 2.7557514543922205e-6);
    p = fma (p, f, 2.4801491039429033e-5);
    p = fma (p, f, 1.9841269589083001e-4);
    p = fma (p, f, 1.3888888945916664e-3);
    p = fma (p, f, 8.3333333334557492e-3);
    p = fma (p, f, 4.1666666666519782e-2);
    p = fma (p, f, 1.6666666666666477e-1);
    p = fma (p, f, 5.0000000000000122e-1);
    p = fma (p, f, 1.0000000000000000e+0);
    p = fma (p, f, 1.0000000000000000e+0);

    // exp(a) = 2**i * exp(f);
    int rlo = __double2loint (p);
    int rhi = (i << 20) + __double2hiint (p);
    r = __hiloint2double (rhi, rlo);

    // handle special cases
    int ia = __double2hiint (a);
    int ib = __double2hiint (708.0); // |a| >= 708 requires double scaling
    int ic = __double2hiint (746.0); // |a| >= 746 severe overflow / underflow
    float fa = __int_as_float (ia);
    float fb = __int_as_float (ib);
    float fc = __int_as_float (ic);

    if (! (fabsf (fa) < fb)) { // !(|a| < 708)
        i = (i > 0) ?  0 : 0x80300000;
        r = __hiloint2double (0x7fe00000 + i, 0);
        r = r * __hiloint2double (rhi - i - 0x3ff00000, rlo);
        if (! (fabsf (fa) < fc)) { // !(|a| < 746)
            r = __hiloint2double ((ia > 0) ? 0x7ff00000 : 0, 0); // +INF, +0
            if (isnan (a)) r = a + a;
        }
    }
    return r;
}
