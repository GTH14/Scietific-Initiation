/*
 * File: cosh.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao_estacionario.h"
#include "cosh.h"

/* Function Definitions */

/*
 * Arguments    : creal_T *x
 * Return Type  : void
 */
void b_cosh(creal_T *x)
{
  double x_re;
  double x_im;
  if (x->im == 0.0) {
    x->re = cosh(x->re);
    x->im = 0.0;
  } else {
    x_re = x->re;
    x_im = x->im;
    x->re = cosh(x->re) * cos(x->im);
    x->im = sinh(x_re) * sin(x_im);
  }
}

/*
 * File trailer for cosh.c
 *
 * [EOF]
 */
