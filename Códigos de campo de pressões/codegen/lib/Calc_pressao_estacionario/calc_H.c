/*
 * File: calc_H.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao_estacionario.h"
#include "calc_H.h"
#include "cbesj.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T v
 * Return Type  : creal_T
 */
creal_T calc_H(const creal_T v)
{
  creal_T H;
  int ierr;
  creal_T w;
  int unusedU0;
  double ar;
  double ai;
  double brm;
  double bim;
  double d;
  ierr = 0;
  if (rtIsNaN(v.re) || rtIsNaN(v.im)) {
    w.re = rtNaN;
    w.im = 0.0;
  } else {
    cbesj(v, 1.0, 1, &w, &unusedU0, &ierr);
  }

  if (ierr == 5) {
    w.re = rtNaN;
    w.im = 0.0;
  } else {
    if (ierr == 2) {
      w.re = rtInf;
      w.im = 0.0;
    }
  }

  if (v.im == 0.0) {
    w.im = 0.0;
  }

  ar = 2.0 * w.re;
  ai = 2.0 * w.im;
  if (v.im == 0.0) {
    if (ai == 0.0) {
      H.re = ar / v.re;
      H.im = 0.0;
    } else if (ar == 0.0) {
      H.re = 0.0;
      H.im = ai / v.re;
    } else {
      H.re = ar / v.re;
      H.im = ai / v.re;
    }
  } else if (v.re == 0.0) {
    if (ar == 0.0) {
      H.re = ai / v.im;
      H.im = 0.0;
    } else if (ai == 0.0) {
      H.re = 0.0;
      H.im = -(ar / v.im);
    } else {
      H.re = ai / v.im;
      H.im = -(ar / v.im);
    }
  } else {
    brm = fabs(v.re);
    bim = fabs(v.im);
    if (brm > bim) {
      bim = v.im / v.re;
      d = v.re + bim * v.im;
      H.re = (ar + bim * ai) / d;
      H.im = (ai - bim * ar) / d;
    } else if (bim == brm) {
      if (v.re > 0.0) {
        bim = 0.5;
      } else {
        bim = -0.5;
      }

      if (v.im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }

      H.re = (ar * bim + ai * d) / brm;
      H.im = (ai * bim - ar * d) / brm;
    } else {
      bim = v.re / v.im;
      d = v.im + bim * v.re;
      H.re = (bim * ar + ai) / d;
      H.im = (bim * ai - ar) / d;
    }
  }

  return H;
}

/*
 * File trailer for calc_H.c
 *
 * [EOF]
 */
