/*
 * File: sum.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 14:00:11
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Plot_U.h"
#include "sum.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T x[64]
 * Return Type  : creal_T
 */
creal_T sum(const creal_T x[64])
{
  creal_T y;
  int k;
  y = x[0];
  for (k = 0; k < 63; k++) {
    y.re += x[k + 1].re;
    y.im += x[k + 1].im;
  }

  return y;
}

/*
 * File trailer for sum.c
 *
 * [EOF]
 */
