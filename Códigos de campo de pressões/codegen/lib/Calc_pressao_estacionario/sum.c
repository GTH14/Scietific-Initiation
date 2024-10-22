/*
 * File: sum.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao_estacionario.h"
#include "sum.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T x[128]
 * Return Type  : creal_T
 */
creal_T sum(const creal_T x[128])
{
  creal_T y;
  int k;
  y = x[0];
  for (k = 0; k < 127; k++) {
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
