/*
 * File: calc_p.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "f_Simulacao_dinamica_v3.h"
#include "calc_p.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T Pax
 *                const creal_T H
 *                double omega
 *                double k
 *                const creal_T r
 *                double phase
 * Return Type  : creal_T
 */
creal_T calc_p(const creal_T Pax, const creal_T H, double omega, double k, const
               creal_T r, double phase)
{
  creal_T p;
  double b_r;
  double omega_im;
  double x_re;
  double x_im;
  b_r = (omega * 0.0 - k * r.re) + phase;
  omega_im = 0.0 - k * r.im;
  x_re = 0.0 * b_r - omega_im;
  x_im = 0.0 * omega_im + b_r;
  if (x_im == 0.0) {
    x_re = exp(x_re);
    x_im = 0.0;
  } else if (rtIsInf(x_im) && rtIsInf(x_re) && (x_re < 0.0)) {
    x_re = 0.0;
    x_im = 0.0;
  } else {
    b_r = exp(x_re / 2.0);
    x_re = b_r * (b_r * cos(x_im));
    x_im = b_r * (b_r * sin(x_im));
  }

  omega_im = Pax.re * H.re - Pax.im * H.im;
  b_r = Pax.re * H.im + Pax.im * H.re;
  p.re = omega_im * x_re - b_r * x_im;
  p.im = omega_im * x_im + b_r * x_re;
  return p;
}

/*
 * File trailer for calc_p.c
 *
 * [EOF]
 */
