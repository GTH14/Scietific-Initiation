/*
 * File: casyi.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 08-Aug-2022 21:44:19
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao.h"
#include "casyi.h"
#include "cmlri.h"
#include "sqrt1.h"
#include "Calc_pressao_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int nin
 *                creal_T *y
 *                double rl
 *                double tol
 *                double elim
 * Return Type  : int
 */
int casyi(const creal_T z, double fnu, int kode, int nin, creal_T *y, double rl,
          double tol, double elim)
{
  int nz;
  int n;
  double x;
  creal_T ak1;
  double brm;
  double acz;
  double s;
  double cz_re;
  double cz_im;
  double dnu2;
  double fdn;
  double ez_re;
  double ez_im;
  double aez;
  double b_dnu2;
  int jl;
  double p1_re;
  int inu;
  double bk;
  boolean_T exitg1;
  double sqk;
  double b_atol;
  double sgn;
  double cs1_re;
  double cs1_im;
  double cs2_re;
  double cs2_im;
  double ak;
  double aa;
  double bb;
  double dk_re;
  double dk_im;
  boolean_T errflag;
  boolean_T exitg2;
  double ni;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  nz = 0;
  x = z.re;
  if (z.im == 0.0) {
    ak1.re = 0.15915494309189535 / z.re;
    ak1.im = 0.0;
  } else if (z.re == 0.0) {
    ak1.re = 0.0;
    ak1.im = -(0.15915494309189535 / z.im);
  } else {
    brm = fabs(z.re);
    acz = fabs(z.im);
    if (brm > acz) {
      s = z.im / z.re;
      acz = z.re + s * z.im;
      ak1.re = (0.15915494309189535 + s * 0.0) / acz;
      ak1.im = (0.0 - s * 0.15915494309189535) / acz;
    } else if (acz == brm) {
      if (z.re > 0.0) {
        acz = 0.5;
      } else {
        acz = -0.5;
      }

      if (z.im > 0.0) {
        dnu2 = 0.5;
      } else {
        dnu2 = -0.5;
      }

      ak1.re = (0.15915494309189535 * acz + 0.0 * dnu2) / brm;
      ak1.im = (0.0 * acz - 0.15915494309189535 * dnu2) / brm;
    } else {
      s = z.re / z.im;
      acz = z.im + s * z.re;
      ak1.re = s * 0.15915494309189535 / acz;
      ak1.im = (s * 0.0 - 0.15915494309189535) / acz;
    }
  }

  b_sqrt(&ak1);
  if (kode == 2) {
    cz_re = 0.0;
    cz_im = z.im;
    acz = 0.0;
  } else {
    cz_re = z.re;
    cz_im = z.im;
    acz = z.re;
  }

  if (fabs(acz) > elim) {
    nz = -1;
    y->re = rtNaN;
    y->im = 0.0;
  } else {
    dnu2 = fnu + fnu;
    if (cz_im == 0.0) {
      cz_re = exp(cz_re);
      cz_im = 0.0;
    } else if (rtIsInf(cz_im) && rtIsInf(cz_re) && (cz_re < 0.0)) {
      cz_re = 0.0;
      cz_im = 0.0;
    } else {
      acz = exp(cz_re / 2.0);
      cz_re = acz * (acz * cos(cz_im));
      cz_im = acz * (acz * sin(cz_im));
    }

    acz = ak1.re;
    ak1.re = ak1.re * cz_re - ak1.im * cz_im;
    ak1.im = acz * cz_im + ak1.im * cz_re;
    fdn = 0.0;
    if (dnu2 > 4.7170688552396617E-153) {
      fdn = dnu2 * dnu2;
    }

    ez_re = 8.0 * z.re;
    ez_im = 8.0 * z.im;
    aez = 8.0 * rt_hypotd_snf(z.re, z.im);
    s = tol / aez;
    dnu2 = rl + rl;
    if (dnu2 < 0.0) {
      b_dnu2 = ceil(dnu2);
    } else {
      b_dnu2 = floor(dnu2);
    }

    jl = (int)b_dnu2 + 2;
    if (z.im != 0.0) {
      inu = (int)fnu;
      acz = (fnu - (double)inu) * 3.1415926535897931;
      dnu2 = sin(acz);
      bk = cos(acz);
      if (z.im < 0.0) {
        bk = -bk;
      }

      if ((inu & 1) != 0) {
        p1_re = -(-dnu2);
        bk = -bk;
      } else {
        p1_re = -dnu2;
      }
    } else {
      p1_re = 0.0;
      bk = 0.0;
    }

    inu = 1;
    exitg1 = false;
    while ((!exitg1) && (inu <= n)) {
      sqk = fdn - 1.0;
      b_atol = s * fabs(fdn - 1.0);
      sgn = 1.0;
      cs1_re = 1.0;
      cs1_im = 0.0;
      cs2_re = 1.0;
      cs2_im = 0.0;
      cz_re = 1.0;
      cz_im = 0.0;
      ak = 0.0;
      aa = 1.0;
      bb = aez;
      dk_re = ez_re;
      dk_im = ez_im;
      errflag = true;
      inu = 1;
      exitg2 = false;
      while ((!exitg2) && (inu <= jl)) {
        cz_re *= sqk;
        cz_im *= sqk;
        acz = cz_re;
        if (dk_im == 0.0) {
          if (cz_im == 0.0) {
            cz_re /= dk_re;
            cz_im = 0.0;
          } else if (cz_re == 0.0) {
            cz_re = 0.0;
            cz_im /= dk_re;
          } else {
            cz_re /= dk_re;
            cz_im /= dk_re;
          }
        } else if (dk_re == 0.0) {
          if (cz_re == 0.0) {
            cz_re = cz_im / dk_im;
            cz_im = 0.0;
          } else if (cz_im == 0.0) {
            cz_re = 0.0;
            cz_im = -(acz / dk_im);
          } else {
            cz_re = cz_im / dk_im;
            cz_im = -(acz / dk_im);
          }
        } else {
          brm = fabs(dk_re);
          acz = fabs(dk_im);
          if (brm > acz) {
            dnu2 = dk_im / dk_re;
            acz = dk_re + dnu2 * dk_im;
            ni = cz_im - dnu2 * cz_re;
            cz_re = (cz_re + dnu2 * cz_im) / acz;
            cz_im = ni / acz;
          } else if (acz == brm) {
            if (dk_re > 0.0) {
              acz = 0.5;
            } else {
              acz = -0.5;
            }

            if (dk_im > 0.0) {
              dnu2 = 0.5;
            } else {
              dnu2 = -0.5;
            }

            ni = cz_im * acz - cz_re * dnu2;
            cz_re = (cz_re * acz + cz_im * dnu2) / brm;
            cz_im = ni / brm;
          } else {
            dnu2 = dk_re / dk_im;
            acz = dk_im + dnu2 * dk_re;
            ni = dnu2 * cz_im - cz_re;
            cz_re = (dnu2 * cz_re + cz_im) / acz;
            cz_im = ni / acz;
          }
        }

        cs2_re += cz_re;
        cs2_im += cz_im;
        sgn = -sgn;
        cs1_re += cz_re * sgn;
        cs1_im += cz_im * sgn;
        dk_re += ez_re;
        dk_im += ez_im;
        aa = aa * fabs(sqk) / bb;
        bb += aez;
        ak += 8.0;
        sqk -= ak;
        if (aa <= b_atol) {
          errflag = false;
          exitg2 = true;
        } else {
          inu++;
        }
      }

      if (errflag) {
        nz = -2;
        exitg1 = true;
      } else {
        if (x + x < elim) {
          cz_re = -2.0 * z.re;
          cz_im = -2.0 * z.im;
          if (cz_im == 0.0) {
            cz_re = exp(cz_re);
            cz_im = 0.0;
          } else if (rtIsInf(cz_im) && rtIsInf(cz_re) && (cz_re < 0.0)) {
            cz_re = 0.0;
            cz_im = 0.0;
          } else {
            acz = exp(cz_re / 2.0);
            cz_re = acz * (acz * cos(cz_im));
            cz_im = acz * (acz * sin(cz_im));
          }

          acz = cz_re;
          cz_re = cz_re * cs2_re - cz_im * cs2_im;
          cz_im = acz * cs2_im + cz_im * cs2_re;
          acz = cz_re;
          cz_re = cz_re * p1_re - cz_im * bk;
          cz_im = acz * bk + cz_im * p1_re;
          cs1_re += cz_re;
          cs1_im += cz_im;
        }

        fdn = (fdn + 8.0 * fnu) + 4.0;
        p1_re = -p1_re;
        bk = -bk;
        y->re = cs1_re * ak1.re - cs1_im * ak1.im;
        y->im = cs1_re * ak1.im + cs1_im * ak1.re;
        inu = 2;
      }
    }
  }

  return nz;
}

/*
 * File trailer for casyi.c
 *
 * [EOF]
 */
