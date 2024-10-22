/*
 * File: exp.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao_estacionario.h"
#include "exp.h"

/* Function Definitions */

/*
 * Arguments    : creal_T *x
 * Return Type  : void
 */
void b_exp(creal_T *x)
{
  double r;
  double x_im;
  if (x->im == 0.0) {
    x->re = exp(x->re);
    x->im = 0.0;
  } else if (rtIsInf(x->im) && rtIsInf(x->re) && (x->re < 0.0)) {
    x->re = 0.0;
    x->im = 0.0;
  } else {
    r = exp(x->re / 2.0);
    x_im = x->im;
    x->re = r * (r * cos(x->im));
    x->im = r * (r * sin(x_im));
  }
}

/*
 * File trailer for exp.c
 *
 * [EOF]
 */
