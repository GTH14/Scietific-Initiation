/*
 * File: calc_v.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 12:47:38
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Plot_Uz.h"
#include "calc_v.h"

/* Function Definitions */

/*
 * Arguments    : double k
 *                double a
 *                const creal_T theta
 * Return Type  : creal_T
 */
creal_T calc_v(double k, double a, const creal_T theta)
{
  creal_T v;
  double x_re;
  double x_im;
  if (theta.im == 0.0) {
    x_re = sin(theta.re);
    x_im = 0.0;
  } else if (theta.re == 0.0) {
    x_re = 0.0;
    x_im = sinh(theta.im);
  } else {
    x_re = sin(theta.re) * cosh(theta.im);
    x_im = cos(theta.re) * sinh(theta.im);
  }

  v.re = k * a * x_re + 1.0E-12;
  v.im = k * a * x_im;
  return v;
}

/*
 * File trailer for calc_v.c
 *
 * [EOF]
 */
