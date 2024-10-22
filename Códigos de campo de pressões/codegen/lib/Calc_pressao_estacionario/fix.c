/*
 * File: fix.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao_estacionario.h"
#include "fix.h"

/* Function Definitions */

/*
 * Arguments    : double *x
 * Return Type  : void
 */
void b_fix(double *x)
{
  if (*x < 0.0) {
    *x = ceil(*x);
  } else {
    *x = floor(*x);
  }
}

/*
 * File trailer for fix.c
 *
 * [EOF]
 */
