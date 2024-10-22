/*
 * File: cseri.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 08-Aug-2022 21:44:19
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao.h"
#include "cseri.h"
#include "cuchk.h"
#include "gammaln.h"
#include "cmlri.h"
#include "Calc_pressao_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int nin
 *                creal_T *y
 *                double tol
 *                double elim
 *                double alim
 * Return Type  : int
 */
int cseri(const creal_T z, double fnu, int kode, int nin, creal_T *y, double tol,
          double elim, double alim)
{
  int nz;
  int n;
  double az;
  double crsc_re;
  boolean_T iflag;
  double hz_re;
  double hz_im;
  double cz_re;
  double cz_im;
  double acz;
  creal_T ck;
  double ak;
  double ascle;
  double aa;
  double coef_re;
  double coef_im;
  double b_atol;
  int i;
  boolean_T exitg1;
  double s1_re;
  double s1_im;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  nz = 0;
  az = rt_hypotd_snf(z.re, z.im);
  if (az == 0.0) {
    if (fnu == 0.0) {
      y->re = 1.0;
      y->im = 0.0;
    } else {
      y->re = 0.0;
      y->im = 0.0;
    }
  } else {
    crsc_re = 1.0;
    iflag = false;
    if (az < 2.2250738585072014E-305) {
      nz = n;
      if (fnu == 0.0) {
        nz = n - 1;
        y->re = 1.0;
        y->im = 0.0;
      } else {
        y->re = 0.0;
        y->im = 0.0;
      }
    } else {
      hz_re = 0.5 * z.re;
      hz_im = 0.5 * z.im;
      if (az > 4.7170688552396617E-153) {
        cz_re = hz_re * hz_re - hz_im * hz_im;
        cz_im = hz_re * hz_im + hz_im * hz_re;
        acz = rt_hypotd_snf(cz_re, cz_im);
      } else {
        cz_re = 0.0;
        cz_im = 0.0;
        acz = 0.0;
      }

      ck.re = hz_re;
      ck.im = hz_im;
      if ((hz_im == 0.0) && rtIsNaN(hz_re)) {
      } else if ((fabs(hz_re) > 8.9884656743115785E+307) || (fabs(hz_im) >
                  8.9884656743115785E+307)) {
        ck.re = log(rt_hypotd_snf(hz_re / 2.0, hz_im / 2.0)) +
          0.69314718055994529;
        ck.im = rt_atan2d_snf(hz_im, hz_re);
      } else {
        ck.re = log(rt_hypotd_snf(hz_re, hz_im));
        ck.im = rt_atan2d_snf(hz_im, hz_re);
      }

      az = (fnu + (double)n) - 1.0;
      ak = az + 1.0;
      gammaln(&ak);
      ck.re = ck.re * az - ak;
      ck.im *= az;
      if (kode == 2) {
        ck.re -= z.re;
      }

      if (ck.re > -elim) {
        ascle = 0.0;
        if (ck.re <= -alim) {
          iflag = true;
          crsc_re = tol;
          ascle = 2.2250738585072014E-305 / tol;
        }

        aa = exp(ck.re);
        if (iflag) {
          aa /= tol;
        }

        coef_re = aa * cos(ck.im);
        coef_im = aa * sin(ck.im);
        b_atol = tol * acz / (az + 1.0);
        i = 1;
        exitg1 = false;
        while ((!exitg1) && (i <= n)) {
          s1_re = 1.0;
          s1_im = 0.0;
          if (!(acz < tol * (fnu + 1.0))) {
            ck.re = 1.0;
            ck.im = 0.0;
            ak = (fnu + 1.0) + 2.0;
            az = fnu + 1.0;
            aa = 2.0;
            do {
              hz_re = 1.0 / az;
              hz_im = ck.re;
              ck.re = hz_re * (ck.re * cz_re - ck.im * cz_im);
              ck.im = hz_re * (hz_im * cz_im + ck.im * cz_re);
              s1_re += ck.re;
              s1_im += ck.im;
              az += ak;
              ak += 2.0;
              aa = aa * acz * hz_re;
            } while (!!(aa > b_atol));
          }

          ck.re = s1_re * coef_re - s1_im * coef_im;
          ck.im = s1_re * coef_im + s1_im * coef_re;
          if (iflag && (cuchk(ck, ascle, tol) != 0)) {
            exitg1 = true;
          } else {
            y->re = ck.re * crsc_re - ck.im * 0.0;
            y->im = ck.re * 0.0 + ck.im * crsc_re;
            i = 2;
          }
        }
      } else {
        nz = 1;
        y->re = 0.0;
        y->im = 0.0;
        if (acz > az) {
          nz = -1;
        }
      }
    }
  }

  return nz;
}

/*
 * File trailer for cseri.c
 *
 * [EOF]
 */
