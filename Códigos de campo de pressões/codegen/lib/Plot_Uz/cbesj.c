/*
 * File: cbesj.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 12:47:38
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Plot_Uz.h"
#include "cbesj.h"
#include "cmlri.h"
#include "casyi.h"
#include "cuchk.h"
#include "gammaln.h"
#include "Plot_Uz_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T z
 *                creal_T *cy
 *                int *nz
 *                int *ierr
 * Return Type  : void
 */
void cbesj(const creal_T z, creal_T *cy, int *nz, int *ierr)
{
  double az;
  int csgn_im;
  creal_T zn;
  boolean_T guard1 = false;
  int nw;
  double b_az;
  double crsc_re;
  boolean_T iflag;
  creal_T hz;
  double cz_re;
  double c_az;
  double acz;
  double b_atol;
  double hz_re;
  double ak;
  double ascle;
  double aa;
  double coef_re;
  double coef_im;
  int i;
  boolean_T exitg1;
  double s1_re;
  double s1_im;
  double s;
  double rs;
  *ierr = 0;
  az = rt_hypotd_snf(z.re, z.im);
  if (az > 1.0737418235E+9) {
    *ierr = 4;
  } else {
    if (az > 32767.999992370605) {
      *ierr = 3;
    }
  }

  csgn_im = 1;
  zn.re = -z.re * 0.0 - (-z.im);
  zn.im = -z.re + -z.im * 0.0;
  if (z.im < 0.0) {
    zn.re = -zn.re;
    zn.im = -zn.im;
    csgn_im = -1;
  }

  cy->re = 0.0;
  cy->im = 0.0;
  az = rt_hypotd_snf(zn.re, zn.im);
  guard1 = false;
  if ((az <= 2.0) || (az * az * 0.25 <= 2.0)) {
    nw = 0;
    b_az = rt_hypotd_snf(zn.re, zn.im);
    if (!(b_az == 0.0)) {
      crsc_re = 1.0;
      iflag = false;
      if (b_az < 2.2250738585072014E-305) {
        nw = 1;
      } else {
        hz.re = 0.5 * zn.re;
        hz.im = 0.5 * zn.im;
        if (b_az > 4.7170688552396617E-153) {
          cz_re = hz.re * hz.re - hz.im * hz.im;
          b_az = hz.re * hz.im + hz.im * hz.re;
          acz = rt_hypotd_snf(cz_re, b_az);
        } else {
          cz_re = 0.0;
          b_az = 0.0;
          acz = 0.0;
        }

        if (hz.im == 0.0) {
          if (hz.re < 0.0) {
            hz.re = log(fabs(hz.re));
            hz.im = 3.1415926535897931;
          } else {
            hz.re = log(hz.re);
            hz.im = 0.0;
          }
        } else if ((fabs(hz.re) > 8.9884656743115785E+307) || (fabs(hz.im) >
                    8.9884656743115785E+307)) {
          hz_re = hz.re;
          hz.re = log(rt_hypotd_snf(hz.re / 2.0, hz.im / 2.0)) +
            0.69314718055994529;
          hz.im = rt_atan2d_snf(hz.im, hz_re);
        } else {
          hz_re = hz.re;
          hz.re = log(rt_hypotd_snf(hz.re, hz.im));
          hz.im = rt_atan2d_snf(hz.im, hz_re);
        }

        ak = 2.0;
        gammaln(&ak);
        hz.re -= ak;
        if (hz.re > -701.61506577445994) {
          ascle = 0.0;
          if (hz.re <= -665.56491761372422) {
            iflag = true;
            crsc_re = 2.2204460492503131E-16;
            ascle = 1.0020841800044864E-289;
          }

          aa = exp(hz.re);
          if (iflag) {
            aa /= 2.2204460492503131E-16;
          }

          coef_re = aa * cos(hz.im);
          coef_im = aa * sin(hz.im);
          b_atol = 2.2204460492503131E-16 * acz / 2.0;
          i = 1;
          exitg1 = false;
          while ((!exitg1) && (i <= 1)) {
            s1_re = 1.0;
            s1_im = 0.0;
            if (!(acz < 4.4408920985006262E-16)) {
              hz.re = 1.0;
              hz.im = 0.0;
              ak = 4.0;
              s = 2.0;
              aa = 2.0;
              do {
                rs = 1.0 / s;
                hz_re = hz.re;
                hz.re = rs * (hz.re * cz_re - hz.im * b_az);
                hz.im = rs * (hz_re * b_az + hz.im * cz_re);
                s1_re += hz.re;
                s1_im += hz.im;
                s += ak;
                ak += 2.0;
                aa = aa * acz * rs;
              } while (!!(aa > b_atol));
            }

            hz.re = s1_re * coef_re - s1_im * coef_im;
            hz.im = s1_re * coef_im + s1_im * coef_re;
            if (iflag && (cuchk(hz, ascle) != 0)) {
              exitg1 = true;
            } else {
              cy->re = hz.re * crsc_re - hz.im * 0.0;
              cy->im = hz.re * 0.0 + hz.im * crsc_re;
              i = 2;
            }
          }
        } else {
          nw = 1;
          if (acz > 1.0) {
            nw = -1;
          }
        }
      }
    }

    if (nw < 0) {
      *nz = 1;
    } else {
      *nz = nw;
    }

    if ((1 - *nz == 0) || (nw >= 0)) {
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    if (az < 21.784271729432426) {
      nw = cmlri(zn, 1.0, 1, 1, cy);
      if (nw < 0) {
        if (nw == -2) {
          *nz = -2;
        } else {
          *nz = -1;
        }
      } else {
        *nz = 0;
      }
    } else {
      nw = casyi(zn, 1.0, 1, cy, 701.61506577445994);
      if (nw < 0) {
        if (nw == -2) {
          *nz = -2;
        } else {
          *nz = -1;
        }
      } else {
        *nz = 0;
      }
    }
  }

  if (*nz < 0) {
    if (*nz == -2) {
      *nz = 0;
      *ierr = 5;
      cy->re = rtNaN;
      cy->im = rtNaN;
    } else {
      *nz = 0;
      *ierr = 2;
      cy->re = rtInf;
      cy->im = 0.0;
    }
  } else {
    if (*nz != 1) {
      zn = *cy;
      b_az = fabs(cy->re);
      cz_re = fabs(cy->im);
      if ((b_az > cz_re) || rtIsNaN(cz_re)) {
        c_az = b_az;
      } else {
        c_az = cz_re;
      }

      if (c_az <= 1.0020841800044864E-289) {
        zn.re = 4.503599627370496E+15 * cy->re;
        zn.im = 4.503599627370496E+15 * cy->im;
        b_atol = 2.2204460492503131E-16;
      } else {
        b_atol = 1.0;
      }

      cy->re = b_atol * (zn.re * 6.123233995736766E-17 - zn.im * (double)csgn_im);
      cy->im = b_atol * (zn.re * (double)csgn_im + zn.im * 6.123233995736766E-17);
    }
  }
}

/*
 * File trailer for cbesj.c
 *
 * [EOF]
 */
