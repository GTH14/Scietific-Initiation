/*
 * File: cmlri.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "f_Simulacao_dinamica_v3.h"
#include "cmlri.h"
#include "Calc_Gorkov.h"
#include "gammaln.h"
#include "f_Simulacao_dinamica_v3_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int nin
 *                creal_T *y
 * Return Type  : int
 */
int cmlri(const creal_T z, double fnu, int kode, int nin, creal_T *y)
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
  double ap;
  double rz_re;
  double rz_im;
  double ack;
  double rho;
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
  double fnf;
  double s_im;
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
      ap = fabs(z.im);
      if (tst > ap) {
        ack = z.im / z.re;
        rho = z.re + ack * z.im;
        ck_re = ((flooraz + 1.0) + ack * 0.0) / rho;
        ck_im = (0.0 - ack * (flooraz + 1.0)) / rho;
      } else if (ap == tst) {
        if (z.re > 0.0) {
          ack = 0.5;
        } else {
          ack = -0.5;
        }

        if (z.im > 0.0) {
          rho = 0.5;
        } else {
          rho = -0.5;
        }

        ck_re = ((flooraz + 1.0) * ack + 0.0 * rho) / tst;
        ck_im = (0.0 * ack - (flooraz + 1.0) * rho) / tst;
      } else {
        ack = z.re / z.im;
        rho = z.im + ack * z.re;
        ck_re = ack * (flooraz + 1.0) / rho;
        ck_im = (ack * 0.0 - (flooraz + 1.0)) / rho;
      }
    }

    if (z.re == 0.0) {
      rz_re = 0.0;
      rz_im = -(2.0 / z.im);
    } else {
      tst = fabs(z.re);
      ap = fabs(z.im);
      if (tst > ap) {
        ack = z.im / z.re;
        rho = z.re + ack * z.im;
        rz_re = (2.0 + ack * 0.0) / rho;
        rz_im = (0.0 - ack * 2.0) / rho;
      } else if (ap == tst) {
        if (z.re > 0.0) {
          ack = 0.5;
        } else {
          ack = -0.5;
        }

        if (z.im > 0.0) {
          rho = 0.5;
        } else {
          rho = -0.5;
        }

        rz_re = (2.0 * ack + 0.0 * rho) / tst;
        rz_im = (0.0 * ack - 2.0 * rho) / tst;
      } else {
        ack = z.re / z.im;
        rho = z.im + ack * z.re;
        rz_re = ack * 2.0 / rho;
        rz_im = (ack * 0.0 - 2.0) / rho;
      }
    }
  }

  p1_re = 0.0;
  p1_im = 0.0;
  p2_re = 1.0;
  p2_im = 0.0;
  ack = ((flooraz + 1.0) + 1.0) / az;
  rho = ack + sqrt(ack * ack - 1.0);
  ack = rho * rho;
  tst = (ack + ack) / ((ack - 1.0) * (rho - 1.0)) / 2.2204460492503131E-16;
  flooraz++;
  earlyExit = true;
  icounter = 1;
  i = 1;
  exitg1 = false;
  while ((!exitg1) && (i < 81)) {
    icounter++;
    pt_re = p2_re;
    pt_im = p2_im;
    ack = ck_re * p2_im + ck_im * p2_re;
    p2_re = p1_re - (ck_re * p2_re - ck_im * p2_im);
    p2_im = p1_im - ack;
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
        ap = fabs(z.im);
        if (tst > ap) {
          ack = z.im / z.re;
          rho = z.re + ack * z.im;
          ck_re = (((double)inu + 1.0) + ack * 0.0) / rho;
          ck_im = (0.0 - ack * ((double)inu + 1.0)) / rho;
        } else if (ap == tst) {
          if (z.re > 0.0) {
            ack = 0.5;
          } else {
            ack = -0.5;
          }

          if (z.im > 0.0) {
            rho = 0.5;
          } else {
            rho = -0.5;
          }

          ck_re = (((double)inu + 1.0) * ack + 0.0 * rho) / tst;
          ck_im = (0.0 * ack - ((double)inu + 1.0) * rho) / tst;
        } else {
          ack = z.re / z.im;
          rho = z.im + ack * z.re;
          ck_re = ack * ((double)inu + 1.0) / rho;
          ck_im = (ack * 0.0 - ((double)inu + 1.0)) / rho;
        }
      }

      tst = sqrt(((double)inu + 1.0) / az / 2.2204460492503131E-16);
      itime = 1;
      earlyExit = true;
      i = 1;
      exitg1 = false;
      while ((!exitg1) && (i < 81)) {
        kcounter++;
        pt_re = p2_re;
        pt_im = p2_im;
        ack = ck_re * p2_im + ck_im * p2_re;
        p2_re = p1_re - (ck_re * p2_re - ck_im * p2_im);
        p2_im = p1_im - ack;
        p1_re = pt_re;
        p1_im = pt_im;
        ck_re += rz_re;
        ck_im += rz_im;
        ap = rt_hypotd_snf(p2_re, p2_im);
        if (ap >= tst * flooraz * flooraz) {
          if (itime == 2) {
            earlyExit = false;
            exitg1 = true;
          } else {
            ack = rt_hypotd_snf(ck_re, ck_im);
            ack += sqrt(ack * ack - 1.0);
            rho = ap / rt_hypotd_snf(pt_re, pt_im);
            if ((ack < rho) || rtIsNaN(rho)) {
              rho = ack;
            }

            tst *= sqrt(rho / (rho * rho - 1.0));
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
      p2_re = 1.0020841800044864E-289;
      p2_im = 0.0;
      fnf = fnu - fixfnu;
      flooraz = fnf + fnf;
      ack = ((double)icounter + flooraz) + 1.0;
      gammaln(&ack);
      rho = (double)icounter + 1.0;
      gammaln(&rho);
      ap = flooraz + 1.0;
      gammaln(&ap);
      tst = exp((ack - rho) - ap);
      fixfnu = 0.0;
      s_im = 0.0;
      icounter -= inu;
      for (i = 1; i <= icounter; i++) {
        pt_re = p2_re;
        pt_im = p2_im;
        ack = (az + fnf) * rz_re;
        ap = (az + fnf) * rz_im;
        rho = ack * p2_im + ap * p2_re;
        p2_re = p1_re + (ack * p2_re - ap * p2_im);
        p2_im = p1_im + rho;
        p1_re = pt_re;
        p1_im = pt_im;
        ack = tst * (1.0 - flooraz / (az + flooraz));
        fixfnu += (ack + tst) * pt_re;
        s_im += (ack + tst) * pt_im;
        tst = ack;
        az--;
      }

      y->re = p2_re;
      y->im = p2_im;
      if (ifnu > 0) {
        for (i = 1; i <= ifnu; i++) {
          pt_re = p2_re;
          pt_im = p2_im;
          ack = (az + fnf) * rz_re;
          ap = (az + fnf) * rz_im;
          rho = ack * p2_im + ap * p2_re;
          p2_re = p1_re + (ack * p2_re - ap * p2_im);
          p2_im = p1_im + rho;
          p1_re = pt_re;
          p1_im = pt_im;
          ack = tst * (1.0 - flooraz / (az + flooraz));
          fixfnu += (ack + tst) * pt_re;
          s_im += (ack + tst) * pt_im;
          tst = ack;
          az--;
        }
      }

      pt_re = z.re;
      pt_im = z.im;
      if (kode == 2) {
        pt_re = z.re - z.re;
        pt_im = z.im;
      }

      if (rz_im == 0.0) {
        if (rz_re < 0.0) {
          ck_re = log(fabs(rz_re));
          ck_im = 3.1415926535897931;
        } else {
          ck_re = log(rz_re);
          ck_im = 0.0;
        }
      } else if ((fabs(rz_re) > 8.9884656743115785E+307) || (fabs(rz_im) >
                  8.9884656743115785E+307)) {
        ck_re = log(rt_hypotd_snf(rz_re / 2.0, rz_im / 2.0)) +
          0.69314718055994529;
        ck_im = rt_atan2d_snf(rz_im, rz_re);
      } else {
        ck_re = log(rt_hypotd_snf(rz_re, rz_im));
        ck_im = rt_atan2d_snf(rz_im, rz_re);
      }

      ack = ck_re;
      ck_re = -fnf * ck_re - -0.0 * ck_im;
      ck_im = -fnf * ck_im + -0.0 * ack;
      ap = 1.0 + fnf;
      gammaln(&ap);
      pt_re = (ck_re + pt_re) - ap;
      pt_im += ck_im;
      p2_re += fixfnu;
      p2_im += s_im;
      p1_re = 1.0 / rt_hypotd_snf(p2_re, p2_im);
      if (pt_im == 0.0) {
        ck_re = exp(pt_re);
        ck_im = 0.0;
      } else if (rtIsInf(pt_im) && rtIsInf(pt_re) && (pt_re < 0.0)) {
        ck_re = 0.0;
        ck_im = 0.0;
      } else {
        ack = exp(pt_re / 2.0);
        ck_re = ack * (ack * cos(pt_im));
        ck_im = ack * (ack * sin(pt_im));
      }

      ack = ck_re * p1_re - ck_im * 0.0;
      ck_im = ck_re * 0.0 + ck_im * p1_re;
      rho = p2_re * p1_re + p2_im * 0.0;
      p2_im = p2_re * 0.0 - p2_im * p1_re;
      ck_re = ack * rho - ck_im * p2_im;
      ck_im = ack * p2_im + ck_im * rho;
      ack = y->re;
      rho = y->im;
      y->re = ack * ck_re - rho * ck_im;
      y->im = ack * ck_im + rho * ck_re;
    }
  }

  return nz;
}

/*
 * File trailer for cmlri.c
 *
 * [EOF]
 */
