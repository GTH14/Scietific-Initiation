/*
 * File: log1.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 08-Aug-2022 21:44:19
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao.h"
#include "log1.h"
#include "cmlri.h"
#include "Calc_pressao_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : creal_T *x
 * Return Type  : void
 */
void b_log(creal_T *x)
{
  double x_im;
  double x_re;
  if ((x->im == 0.0) && rtIsNaN(x->re)) {
  } else if ((fabs(x->re) > 8.9884656743115785E+307) || (fabs(x->im) >
              8.9884656743115785E+307)) {
    x_im = x->im;
    x_re = x->re;
    x->re = log(rt_hypotd_snf(x->re / 2.0, x->im / 2.0)) + 0.69314718055994529;
    x->im = rt_atan2d_snf(x_im, x_re);
  } else {
    x_im = x->im;
    x_re = x->re;
    x->re = log(rt_hypotd_snf(x->re, x->im));
    x->im = rt_atan2d_snf(x_im, x_re);
  }
}

/*
 * File trailer for log1.c
 *
 * [EOF]
 */
