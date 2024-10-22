/*
 * File: casyi.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "f_Simulacao_dinamica_v3.h"
#include "casyi.h"
#include "Calc_Gorkov.h"
#include "sqrt1.h"
#include "f_Simulacao_dinamica_v3_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                creal_T *y
 *                double elim
 * Return Type  : int
 */
int casyi(const creal_T z, double fnu, int kode, creal_T *y, double elim)
{
  int nz;
  double x;
  double cz_re;
  double cz_im;
  double acz;
  double r;
  creal_T ck;
  double brm;
  double s;
  double ak1_re;
  double ak1_im;
  double fdn;
  double sgnbi;
  double ez_re;
  double ez_im;
  double aez;
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
  boolean_T errflag;
  boolean_T exitg2;
  double ck_re;
  nz = 0;
  x = z.re;
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
    acz = fnu + fnu;
    if (cz_im == 0.0) {
      cz_re = exp(cz_re);
      cz_im = 0.0;
    } else if (rtIsInf(cz_im) && rtIsInf(cz_re) && (cz_re < 0.0)) {
      cz_re = 0.0;
      cz_im = 0.0;
    } else {
      r = exp(cz_re / 2.0);
      cz_re = r * (r * cos(cz_im));
      cz_im = r * (r * sin(cz_im));
    }

    if (z.im == 0.0) {
      ck.re = 0.15915494309189535 / z.re;
      ck.im = 0.0;
    } else if (z.re == 0.0) {
      ck.re = 0.0;
      ck.im = -(0.15915494309189535 / z.im);
    } else {
      brm = fabs(z.re);
      r = fabs(z.im);
      if (brm > r) {
        s = z.im / z.re;
        r = z.re + s * z.im;
        ck.re = (0.15915494309189535 + s * 0.0) / r;
        ck.im = (0.0 - s * 0.15915494309189535) / r;
      } else if (r == brm) {
        if (z.re > 0.0) {
          r = 0.5;
        } else {
          r = -0.5;
        }

        if (z.im > 0.0) {
          sgnbi = 0.5;
        } else {
          sgnbi = -0.5;
        }

        ck.re = (0.15915494309189535 * r + 0.0 * sgnbi) / brm;
        ck.im = (0.0 * r - 0.15915494309189535 * sgnbi) / brm;
      } else {
        s = z.re / z.im;
        r = z.im + s * z.re;
        ck.re = s * 0.15915494309189535 / r;
        ck.im = (s * 0.0 - 0.15915494309189535) / r;
      }
    }

    b_sqrt(&ck);
    ak1_re = ck.re * cz_re - ck.im * cz_im;
    ak1_im = ck.re * cz_im + ck.im * cz_re;
    fdn = 0.0;
    if (acz > 4.7170688552396617E-153) {
      fdn = acz * acz;
    }

    ez_re = 8.0 * z.re;
    ez_im = 8.0 * z.im;
    aez = 8.0 * rt_hypotd_snf(z.re, z.im);
    s = 2.2204460492503131E-16 / aez;
    if (z.im != 0.0) {
      inu = (int)fnu;
      acz = (fnu - (double)inu) * 3.1415926535897931;
      r = sin(acz);
      bk = cos(acz);
      if (z.im < 0.0) {
        bk = -bk;
      }

      if ((inu & 1) != 0) {
        p1_re = -(-r);
        bk = -bk;
      } else {
        p1_re = -r;
      }
    } else {
      p1_re = 0.0;
      bk = 0.0;
    }

    inu = 1;
    exitg1 = false;
    while ((!exitg1) && (inu <= 1)) {
      sqk = fdn - 1.0;
      b_atol = s * fabs(fdn - 1.0);
      sgn = 1.0;
      cs1_re = 1.0;
      cs1_im = 0.0;
      cs2_re = 1.0;
      cs2_im = 0.0;
      ck.re = 1.0;
      ck.im = 0.0;
      ak = 0.0;
      aa = 1.0;
      bb = aez;
      cz_re = ez_re;
      cz_im = ez_im;
      errflag = true;
      inu = 1;
      exitg2 = false;
      while ((!exitg2) && (inu <= 45)) {
        ck.re *= sqk;
        ck.im *= sqk;
        ck_re = ck.re;
        if (cz_im == 0.0) {
          if (ck.im == 0.0) {
            ck.re /= cz_re;
            ck.im = 0.0;
          } else if (ck.re == 0.0) {
            ck.re = 0.0;
            ck.im /= cz_re;
          } else {
            ck.re /= cz_re;
            ck.im /= cz_re;
          }
        } else if (cz_re == 0.0) {
          if (ck.re == 0.0) {
            ck.re = ck.im / cz_im;
            ck.im = 0.0;
          } else if (ck.im == 0.0) {
            ck.re = 0.0;
            ck.im = -(ck_re / cz_im);
          } else {
            ck.re = ck.im / cz_im;
            ck.im = -(ck_re / cz_im);
          }
        } else {
          brm = fabs(cz_re);
          r = fabs(cz_im);
          if (brm > r) {
            acz = cz_im / cz_re;
            r = cz_re + acz * cz_im;
            ck.re = (ck.re + acz * ck.im) / r;
            ck.im = (ck.im - acz * ck_re) / r;
          } else if (r == brm) {
            if (cz_re > 0.0) {
              r = 0.5;
            } else {
              r = -0.5;
            }

            if (cz_im > 0.0) {
              sgnbi = 0.5;
            } else {
              sgnbi = -0.5;
            }

            ck.re = (ck.re * r + ck.im * sgnbi) / brm;
            ck.im = (ck.im * r - ck_re * sgnbi) / brm;
          } else {
            acz = cz_re / cz_im;
            r = cz_im + acz * cz_re;
            ck.re = (acz * ck.re + ck.im) / r;
            ck.im = (acz * ck.im - ck_re) / r;
          }
        }

        cs2_re += ck.re;
        cs2_im += ck.im;
        sgn = -sgn;
        cs1_re += ck.re * sgn;
        cs1_im += ck.im * sgn;
        cz_re += ez_re;
        cz_im += ez_im;
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
            r = exp(cz_re / 2.0);
            cz_re = r * (r * cos(cz_im));
            cz_im = r * (r * sin(cz_im));
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
        y->re = cs1_re * ak1_re - cs1_im * ak1_im;
        y->im = cs1_re * ak1_im + cs1_im * ak1_re;
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
