/*
 * File: cmlri.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao_estacionario.h"
#include "cmlri.h"
#include "gammaln.h"
#include "Calc_pressao_estacionario_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int nin
 *                creal_T *y
 *                double tol
 * Return Type  : int
 */
int cmlri(const creal_T z, double fnu, int kode, int nin, creal_T *y, double tol)
{
  int nz;
  double az;
  double flooraz;
  int iaz;
  double fixfnu;
  int ifnu;
  int icounter;
  int inu;
  double ck_re;
  double tst;
  double ck_im;
  double flam;
  double rz_re;
  double rz_im;
  double rho2;
  double ack;
  double p1_re;
  double p1_im;
  double p2_re;
  double p2_im;
  boolean_T earlyExit;
  int i;
  boolean_T exitg1;
  int kcounter;
  double pt_re;
  boolean_T guard1 = false;
  double pt_im;
  int itime;
  nz = 0;
  az = rt_hypotd_snf(z.re, z.im);
  flooraz = floor(az);
  iaz = (int)flooraz;
  if (fnu < 0.0) {
    fixfnu = ceil(fnu);
  } else {
    fixfnu = floor(fnu);
  }

  ifnu = (int)fixfnu;
  if (nin < 1) {
    icounter = nin;
  } else {
    icounter = 1;
  }

  inu = (ifnu + icounter) - 1;
  if (z.im == 0.0) {
    ck_re = (flooraz + 1.0) / z.re;
    ck_im = 0.0;
    rz_re = 2.0 / z.re;
    rz_im = 0.0;
  } else {
    if (z.re == 0.0) {
      ck_re = 0.0;
      ck_im = -((flooraz + 1.0) / z.im);
    } else {
      tst = fabs(z.re);
      flam = fabs(z.im);
      if (tst > flam) {
        rho2 = z.im / z.re;
        ack = z.re + rho2 * z.im;
        ck_re = ((flooraz + 1.0) + rho2 * 0.0) / ack;
        ck_im = (0.0 - rho2 * (flooraz + 1.0)) / ack;
      } else if (flam == tst) {
        if (z.re > 0.0) {
          rho2 = 0.5;
        } else {
          rho2 = -0.5;
        }

        if (z.im > 0.0) {
          ack = 0.5;
        } else {
          ack = -0.5;
        }

        ck_re = ((flooraz + 1.0) * rho2 + 0.0 * ack) / tst;
        ck_im = (0.0 * rho2 - (flooraz + 1.0) * ack) / tst;
      } else {
        rho2 = z.re / z.im;
        ack = z.im + rho2 * z.re;
        ck_re = rho2 * (flooraz + 1.0) / ack;
        ck_im = (rho2 * 0.0 - (flooraz + 1.0)) / ack;
      }
    }

    if (z.re == 0.0) {
      rz_re = 0.0;
      rz_im = -(2.0 / z.im);
    } else {
      tst = fabs(z.re);
      flam = fabs(z.im);
      if (tst > flam) {
        rho2 = z.im / z.re;
        ack = z.re + rho2 * z.im;
        rz_re = (2.0 + rho2 * 0.0) / ack;
        rz_im = (0.0 - rho2 * 2.0) / ack;
      } else if (flam == tst) {
        if (z.re > 0.0) {
          rho2 = 0.5;
        } else {
          rho2 = -0.5;
        }

        if (z.im > 0.0) {
          ack = 0.5;
        } else {
          ack = -0.5;
        }

        rz_re = (2.0 * rho2 + 0.0 * ack) / tst;
        rz_im = (0.0 * rho2 - 2.0 * ack) / tst;
      } else {
        rho2 = z.re / z.im;
        ack = z.im + rho2 * z.re;
        rz_re = rho2 * 2.0 / ack;
        rz_im = (rho2 * 0.0 - 2.0) / ack;
      }
    }
  }

  p1_re = 0.0;
  p1_im = 0.0;
  p2_re = 1.0;
  p2_im = 0.0;
  ack = ((flooraz + 1.0) + 1.0) / az;
  ack += sqrt(ack * ack - 1.0);
  rho2 = ack * ack;
  tst = (rho2 + rho2) / ((rho2 - 1.0) * (ack - 1.0)) / tol;
  flooraz++;
  earlyExit = true;
  icounter = 1;
  i = 1;
  exitg1 = false;
  while ((!exitg1) && (i < 81)) {
    icounter++;
    pt_re = p2_re;
    pt_im = p2_im;
    rho2 = ck_re * p2_im + ck_im * p2_re;
    p2_re = p1_re - (ck_re * p2_re - ck_im * p2_im);
    p2_im = p1_im - rho2;
    p1_re = pt_re;
    p1_im = pt_im;
    ck_re += rz_re;
    ck_im += rz_im;
    if (rt_hypotd_snf(p2_re, p2_im) > tst * flooraz * flooraz) {
      earlyExit = false;
      exitg1 = true;
    } else {
      flooraz++;
      i++;
    }
  }

  if (earlyExit) {
    nz = -2;
  } else {
    kcounter = 1;
    guard1 = false;
    if (inu >= iaz) {
      p1_re = 0.0;
      p1_im = 0.0;
      p2_re = 1.0;
      p2_im = 0.0;
      if (z.im == 0.0) {
        ck_re = ((double)inu + 1.0) / z.re;
        ck_im = 0.0;
      } else if (z.re == 0.0) {
        if ((double)inu + 1.0 == 0.0) {
          ck_re = 0.0 / z.im;
          ck_im = 0.0;
        } else {
          ck_re = 0.0;
          ck_im = -(((double)inu + 1.0) / z.im);
        }
      } else {
        tst = fabs(z.re);
        flam = fabs(z.im);
        if (tst > flam) {
          rho2 = z.im / z.re;
          ack = z.re + rho2 * z.im;
          ck_re = (((double)inu + 1.0) + rho2 * 0.0) / ack;
          ck_im = (0.0 - rho2 * ((double)inu + 1.0)) / ack;
        } else if (flam == tst) {
          if (z.re > 0.0) {
            rho2 = 0.5;
          } else {
            rho2 = -0.5;
          }

          if (z.im > 0.0) {
            ack = 0.5;
          } else {
            ack = -0.5;
          }

          ck_re = (((double)inu + 1.0) * rho2 + 0.0 * ack) / tst;
          ck_im = (0.0 * rho2 - ((double)inu + 1.0) * ack) / tst;
        } else {
          rho2 = z.re / z.im;
          ack = z.im + rho2 * z.re;
          ck_re = rho2 * ((double)inu + 1.0) / ack;
          ck_im = (rho2 * 0.0 - ((double)inu + 1.0)) / ack;
        }
      }

      tst = sqrt(((double)inu + 1.0) / az / tol);
      itime = 1;
      earlyExit = true;
      i = 1;
      exitg1 = false;
      while ((!exitg1) && (i < 81)) {
        kcounter++;
        pt_re = p2_re;
        pt_im = p2_im;
        rho2 = ck_re * p2_im + ck_im * p2_re;
        p2_re = p1_re - (ck_re * p2_re - ck_im * p2_im);
        p2_im = p1_im - rho2;
        p1_re = pt_re;
        p1_im = pt_im;
        ck_re += rz_re;
        ck_im += rz_im;
        rho2 = rt_hypotd_snf(p2_re, p2_im);
        if (rho2 >= tst * flooraz * flooraz) {
          if (itime == 2) {
            earlyExit = false;
            exitg1 = true;
          } else {
            ack = rt_hypotd_snf(ck_re, ck_im);
            flam = ack + sqrt(ack * ack - 1.0);
            ack = rho2 / rt_hypotd_snf(pt_re, pt_im);
            if ((flam < ack) || rtIsNaN(ack)) {
              ack = flam;
            }

            tst *= sqrt(ack / (ack * ack - 1.0));
            itime = 2;
            i++;
          }
        } else {
          i++;
        }
      }

      if (earlyExit) {
        nz = -2;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      itime = icounter + iaz;
      icounter = kcounter + inu;
      if (itime > icounter) {
        icounter = itime;
      }

      az = icounter;
      p1_re = 0.0;
      p1_im = 0.0;
      p2_re = 2.2250738585072014E-305 / tol;
      p2_im = 0.0;
      fixfnu = fnu - fixfnu;
      flooraz = fixfnu + fixfnu;
      rho2 = ((double)icounter + flooraz) + 1.0;
      gammaln(&rho2);
      ack = (double)icounter + 1.0;
      gammaln(&ack);
      flam = flooraz + 1.0;
      gammaln(&flam);
      flam = exp((rho2 - ack) - flam);
      ck_re = 0.0;
      ck_im = 0.0;
      icounter -= inu;
      for (i = 1; i <= icounter; i++) {
        pt_re = p2_re;
        pt_im = p2_im;
        rho2 = (az + fixfnu) * rz_re;
        ack = (az + fixfnu) * rz_im;
        tst = rho2 * p2_im + ack * p2_re;
        p2_re = p1_re + (rho2 * p2_re - ack * p2_im);
        p2_im = p1_im + tst;
        p1_re = pt_re;
        p1_im = pt_im;
        ack = flam * (1.0 - flooraz / (az + flooraz));
        ck_re += (ack + flam) * pt_re;
        ck_im += (ack + flam) * pt_im;
        flam = ack;
        az--;
      }

      y->re = p2_re;
      y->im = p2_im;
      if (ifnu > 0) {
        for (i = 1; i <= ifnu; i++) {
          pt_re = p2_re;
          pt_im = p2_im;
          rho2 = (az + fixfnu) * rz_re;
          ack = (az + fixfnu) * rz_im;
          tst = rho2 * p2_im + ack * p2_re;
          p2_re = p1_re + (rho2 * p2_re - ack * p2_im);
          p2_im = p1_im + tst;
          p1_re = pt_re;
          p1_im = pt_im;
          ack = flam * (1.0 - flooraz / (az + flooraz));
          ck_re += (ack + flam) * pt_re;
          ck_im += (ack + flam) * pt_im;
          flam = ack;
          az--;
        }
      }

      pt_re = z.re;
      pt_im = z.im;
      if (kode == 2) {
        pt_re = z.re - z.re;
        pt_im = z.im;
      }

      if ((rz_im == 0.0) && rtIsNaN(rz_re)) {
      } else if ((fabs(rz_re) > 8.9884656743115785E+307) || (fabs(rz_im) >
                  8.9884656743115785E+307)) {
        rho2 = rz_re;
        rz_re = log(rt_hypotd_snf(rz_re / 2.0, rz_im / 2.0)) +
          0.69314718055994529;
        rz_im = rt_atan2d_snf(rz_im, rho2);
      } else {
        rho2 = rz_re;
        rz_re = log(rt_hypotd_snf(rz_re, rz_im));
        rz_im = rt_atan2d_snf(rz_im, rho2);
      }

      rho2 = rz_re;
      rz_re = -fixfnu * rz_re - -0.0 * rz_im;
      rz_im = -fixfnu * rz_im + -0.0 * rho2;
      rho2 = 1.0 + fixfnu;
      gammaln(&rho2);
      pt_re = (rz_re + pt_re) - rho2;
      pt_im += rz_im;
      p2_re += ck_re;
      p2_im += ck_im;
      p1_re = 1.0 / rt_hypotd_snf(p2_re, p2_im);
      if (pt_im == 0.0) {
        rz_re = exp(pt_re);
        rz_im = 0.0;
      } else if (rtIsInf(pt_im) && rtIsInf(pt_re) && (pt_re < 0.0)) {
        rz_re = 0.0;
        rz_im = 0.0;
      } else {
        rho2 = exp(pt_re / 2.0);
        rz_re = rho2 * (rho2 * cos(pt_im));
        rz_im = rho2 * (rho2 * sin(pt_im));
      }

      rho2 = rz_re * p1_re - rz_im * 0.0;
      rz_im = rz_re * 0.0 + rz_im * p1_re;
      ack = p2_re * p1_re + p2_im * 0.0;
      p2_im = p2_re * 0.0 - p2_im * p1_re;
      rz_re = rho2 * ack - rz_im * p2_im;
      rz_im = rho2 * p2_im + rz_im * ack;
      rho2 = y->re;
      ack = y->im;
      y->re = rho2 * rz_re - ack * rz_im;
      y->im = rho2 * rz_im + ack * rz_re;
    }
  }

  return nz;
}

/*
 * File trailer for cmlri.c
 *
 * [EOF]
 */
