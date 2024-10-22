/*
 * File: cairy.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 08-Aug-2022 21:44:19
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao.h"
#include "cairy.h"
#include "cbknu.h"
#include "cmlri.h"
#include "casyi.h"
#include "cuchk.h"
#include "gammaln.h"
#include "sqrt1.h"
#include "Calc_pressao_rtwutil.h"

/* Function Declarations */
static int cacai(const creal_T z, double fnu, int kode, int mr, creal_T *y,
                 double rl, double tol, double elim, double alim);

/* Function Definitions */

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int mr
 *                creal_T *y
 *                double rl
 *                double tol
 *                double elim
 *                double alim
 * Return Type  : int
 */
static int cacai(const creal_T z, double fnu, int kode, int mr, creal_T *y,
                 double rl, double tol, double elim, double alim)
{
  int nz;
  creal_T zn;
  double az;
  boolean_T guard1 = false;
  int nw;
  double crsc_re;
  boolean_T iflag;
  double hz_re;
  double hz_im;
  double cz_re;
  double cz_im;
  double acz;
  creal_T cy[2];
  double ck_re;
  double aa;
  double ak1_re;
  double ak1_im;
  creal_T s2;
  double ascle;
  double coef_im;
  double b_atol;
  int exitg1;
  double s;
  double rs;
  double b_ak1_re;
  nz = 0;
  zn.re = -z.re;
  zn.im = -z.im;
  az = rt_hypotd_snf(z.re, z.im);
  guard1 = false;
  if ((!(az <= 2.0)) && (az * az * 0.25 > ((fnu + 1.0) - 1.0) + 1.0)) {
    if (az < rl) {
      nw = cmlri(zn, fnu, kode, 1, y, tol);
    } else {
      nw = casyi(zn, fnu, kode, 1, y, rl, tol, elim);
    }

    if (nw < 0) {
      if (nw == -2) {
        nz = -2;
      } else {
        nz = -1;
      }
    } else {
      guard1 = true;
    }
  } else {
    az = rt_hypotd_snf(-z.re, -z.im);
    if (az == 0.0) {
      if (fnu == 0.0) {
        y->re = 1.0;
        y->im = 0.0;
      } else {
        y->re = 0.0;
        y->im = 0.0;
      }

      guard1 = true;
    } else {
      crsc_re = 1.0;
      iflag = false;
      if (az < 2.2250738585072014E-305) {
        if (fnu == 0.0) {
          y->re = 1.0;
          y->im = 0.0;
        } else {
          y->re = 0.0;
          y->im = 0.0;
        }

        guard1 = true;
      } else {
        hz_re = 0.5 * -z.re;
        hz_im = 0.5 * -z.im;
        if (az > 4.7170688552396617E-153) {
          cz_re = hz_re * hz_re - hz_im * hz_im;
          cz_im = hz_re * hz_im + hz_im * hz_re;
          acz = rt_hypotd_snf(cz_re, cz_im);
        } else {
          cz_re = 0.0;
          cz_im = 0.0;
          acz = 0.0;
        }

        ck_re = hz_re;
        if ((hz_im == 0.0) && rtIsNaN(hz_re)) {
        } else if ((fabs(hz_re) > 8.9884656743115785E+307) || (fabs(hz_im) >
                    8.9884656743115785E+307)) {
          ck_re = log(rt_hypotd_snf(hz_re / 2.0, hz_im / 2.0)) +
            0.69314718055994529;
          hz_im = rt_atan2d_snf(hz_im, hz_re);
        } else {
          ck_re = log(rt_hypotd_snf(hz_re, hz_im));
          hz_im = rt_atan2d_snf(hz_im, hz_re);
        }

        az = ((fnu + 1.0) - 1.0) + 1.0;
        gammaln(&az);
        ak1_re = ck_re * ((fnu + 1.0) - 1.0) - az;
        ak1_im = hz_im * ((fnu + 1.0) - 1.0);
        if (kode == 2) {
          ak1_re -= -z.re;
        }

        if (ak1_re > -elim) {
          ascle = 0.0;
          if (ak1_re <= -alim) {
            iflag = true;
            crsc_re = tol;
            ascle = 2.2250738585072014E-305 / tol;
          }

          aa = exp(ak1_re);
          if (iflag) {
            aa /= tol;
          }

          hz_re = aa * cos(ak1_im);
          coef_im = aa * sin(ak1_im);
          b_atol = tol * acz / (((fnu + 1.0) - 1.0) + 1.0);
          nw = 1;
          do {
            exitg1 = 0;
            if (nw <= 1) {
              ck_re = 1.0;
              hz_im = 0.0;
              if (!(acz < tol * (fnu + 1.0))) {
                ak1_re = 1.0;
                ak1_im = 0.0;
                az = (fnu + 1.0) + 2.0;
                s = fnu + 1.0;
                aa = 2.0;
                do {
                  rs = 1.0 / s;
                  b_ak1_re = ak1_re;
                  ak1_re = ak1_re * cz_re - ak1_im * cz_im;
                  ak1_im = b_ak1_re * cz_im + ak1_im * cz_re;
                  ak1_re *= rs;
                  ak1_im *= rs;
                  ck_re += ak1_re;
                  hz_im += ak1_im;
                  s += az;
                  az += 2.0;
                  aa = aa * acz * rs;
                } while (!!(aa > b_atol));
              }

              s2.re = ck_re * hz_re - hz_im * coef_im;
              s2.im = ck_re * coef_im + hz_im * hz_re;
              if (iflag && (cuchk(s2, ascle, tol) != 0)) {
                guard1 = true;
                exitg1 = 1;
              } else {
                y->re = s2.re * crsc_re - s2.im * 0.0;
                y->im = s2.re * 0.0 + s2.im * crsc_re;
                nw = 2;
              }
            } else {
              guard1 = true;
              exitg1 = 1;
            }
          } while (exitg1 == 0);
        } else {
          y->re = 0.0;
          y->im = 0.0;
          guard1 = true;
        }
      }
    }
  }

  if (guard1) {
    for (nw = 0; nw < 2; nw++) {
      cy[nw].re = 0.0;
      cy[nw].im = 0.0;
    }

    nw = b_cbknu(zn, fnu, kode, 1, cy, tol, elim, alim);
    if (nw != 0) {
      if (nw == -2) {
        nz = -2;
      } else {
        nz = -1;
      }
    } else {
      if (mr < 0) {
        aa = 3.1415926535897931;
      } else {
        aa = -3.1415926535897931;
      }

      ck_re = 0.0;
      hz_im = aa;
      if (kode != 1) {
        az = cos(-(-z.im));
        hz_re = sin(-(-z.im));
        ck_re = 0.0 * az - aa * hz_re;
        hz_im = 0.0 * hz_re + aa * az;
      }

      nw = (int)fnu;
      az = (fnu - (double)nw) * aa;
      s2.re = cos(az);
      s2.im = sin(az);
      if ((nw & 1) != 0) {
        s2.re = -s2.re;
        s2.im = -s2.im;
      }

      ak1_re = cy[0].re;
      ak1_im = cy[0].im;
      if (kode != 1) {
        ascle = 2.2250738585072014E-305 / tol;
        az = rt_hypotd_snf(cy[0].re, cy[0].im);
        if (az > 0.0) {
          if ((-(-z.re) - (-z.re)) + log(az) < -alim) {
            ak1_re = 0.0;
            ak1_im = 0.0;
            az = 0.0;
          } else {
            if ((cy[0].im == 0.0) && rtIsNaN(cy[0].re)) {
            } else if ((fabs(cy[0].re) > 8.9884656743115785E+307) || (fabs(cy[0]
                         .im) > 8.9884656743115785E+307)) {
              ak1_re = log(rt_hypotd_snf(cy[0].re / 2.0, cy[0].im / 2.0)) +
                0.69314718055994529;
              ak1_im = rt_atan2d_snf(cy[0].im, cy[0].re);
            } else {
              ak1_re = log(rt_hypotd_snf(cy[0].re, cy[0].im));
              ak1_im = rt_atan2d_snf(cy[0].im, cy[0].re);
            }

            ak1_re -= -z.re;
            ak1_im -= -z.im;
            ak1_re -= -z.re;
            ak1_im -= -z.im;
            if (ak1_im == 0.0) {
              ak1_re = exp(ak1_re);
              ak1_im = 0.0;
            } else if (rtIsInf(ak1_im) && rtIsInf(ak1_re) && (ak1_re < 0.0)) {
              ak1_re = 0.0;
              ak1_im = 0.0;
            } else {
              az = exp(ak1_re / 2.0);
              ak1_re = az * (az * cos(ak1_im));
              ak1_im = az * (az * sin(ak1_im));
            }

            az = rt_hypotd_snf(ak1_re, ak1_im);
          }
        }

        if ((az > ascle) || (rt_hypotd_snf(y->re, y->im) > ascle)) {
          nz = 0;
        } else {
          ak1_re = 0.0;
          ak1_im = 0.0;
          y->re = 0.0;
          y->im = 0.0;
          nz = 1;
        }
      }

      az = y->re;
      hz_re = y->im;
      y->re = ck_re * az - hz_im * hz_re;
      y->im = ck_re * hz_re + hz_im * az;
      y->re += s2.re * ak1_re - s2.im * ak1_im;
      y->im += s2.re * ak1_im + s2.im * ak1_re;
    }
  }

  return nz;
}

/*
 * Arguments    : const creal_T z
 *                int id
 *                int kode
 * Return Type  : creal_T
 */
creal_T cairy(const creal_T z, int id, int kode)
{
  creal_T ai;
  double az;
  double s1_re;
  double az3;
  double s1_im;
  double r;
  creal_T s2;
  creal_T trm2;
  double aa;
  int iflag;
  creal_T trm1;
  double ak;
  double atrm;
  boolean_T guard1 = false;
  double z3_re;
  boolean_T guard2 = false;
  double z3_im;
  boolean_T guard3 = false;
  double d2;
  double ad;
  int i;
  double bk;
  boolean_T exitg1;
  int b_z;
  double b_z3_re;
  double b_z3_im;
  ai.re = 0.0;
  ai.im = 0.0;
  az = rt_hypotd_snf(z.re, z.im);
  if (az > 1.0) {
    az3 = (1.0 + (double)id) / 3.0;
    r = log(az);
    trm2 = z;
    b_sqrt(&trm2);
    s2.re = 0.66666666666666663 * (z.re * trm2.re - z.im * trm2.im);
    s2.im = 0.66666666666666663 * (z.re * trm2.im + z.im * trm2.re);
    iflag = 0;
    aa = 1.0;
    ak = s2.im;
    if (z.re < 0.0) {
      s2.re = -fabs(s2.re);
    }

    if ((z.im != 0.0) || (z.re > 0.0)) {
    } else {
      s2.re = 0.0;
      s2.im = ak;
    }

    guard1 = false;
    guard2 = false;
    guard3 = false;
    if ((s2.re >= 0.0) && (z.re > 0.0)) {
      if ((kode != 2) && (s2.re >= 664.87164553371019)) {
        iflag = 2;
        aa = 4.503599627370496E+15;
        if (-s2.re - 0.25 * r < -700.92179369444591) {
        } else {
          guard3 = true;
        }
      } else {
        guard3 = true;
      }
    } else if ((kode != 2) && (s2.re <= -664.87164553371019)) {
      iflag = 1;
      aa = 2.2204460492503131E-16;
      if (-s2.re + 0.25 * r > 700.92179369444591) {
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard3) {
      cbknu(s2, az3, kode, 664.87164553371019, &trm1, &i);
      guard1 = true;
    }

    if (guard2) {
      trm1.re = 0.0;
      trm1.im = 0.0;
      if (z.im < 0.0) {
        b_z = -1;
      } else {
        b_z = 1;
      }

      i = cacai(s2, az3, kode, b_z, &trm1, 21.784271729432426,
                2.2204460492503131E-16, 700.92179369444591, 664.87164553371019);
      if (i < 0) {
      } else {
        guard1 = true;
      }
    }

    if (guard1) {
      s1_re = 0.18377629847393068 * trm1.re;
      s1_im = 0.18377629847393068 * trm1.im;
      if (iflag != 0) {
        s1_re *= aa;
        s1_im *= aa;
        if (id == 1) {
          s1_re = -s1_re;
          s1_im = -s1_im;
          r = s1_re;
          s1_re = s1_re * z.re - s1_im * z.im;
          s1_im = r * z.im + s1_im * z.re;
        } else {
          r = s1_re;
          s1_re = s1_re * trm2.re - s1_im * trm2.im;
          s1_im = r * trm2.im + s1_im * trm2.re;
        }

        ai.re = 1.0 / aa * s1_re;
        ai.im = 1.0 / aa * s1_im;
      } else if (id == 1) {
        ai.re = -z.re * s1_re - -z.im * s1_im;
        ai.im = -z.re * s1_im + -z.im * s1_re;
      } else {
        ai.re = trm2.re * s1_re - trm2.im * s1_im;
        ai.im = trm2.re * s1_im + trm2.im * s1_re;
      }
    }
  } else {
    s1_re = 1.0;
    s1_im = 0.0;
    s2.re = 1.0;
    s2.im = 0.0;
    if (az < 2.2204460492503131E-16) {
      s1_re = 0.0;
      s1_im = 0.0;
      if (id == 1) {
        if (az > 4.7170688552396617E-153) {
          s1_re = 0.5 * (z.re * z.re - z.im * z.im);
          s1_im = 0.5 * (z.re * z.im + z.im * z.re);
        }

        s1_re *= 0.35502805388781722;
        s1_im *= 0.35502805388781722;
        ai.re = -0.25881940379280682 + s1_re;
        ai.im = s1_im;
      } else {
        if (az > 2.2250738585072014E-305) {
          s1_re = 0.25881940379280682 * z.re;
          s1_im = 0.25881940379280682 * z.im;
        }

        ai.re = 0.35502805388781722 - s1_re;
        ai.im = 0.0 - s1_im;
      }
    } else {
      aa = az * az;
      if (aa >= 2.2204460492503131E-16 / az) {
        trm1.re = 1.0;
        trm1.im = 0.0;
        trm2.re = 1.0;
        trm2.im = 0.0;
        atrm = 1.0;
        r = z.re * z.re - z.im * z.im;
        az3 = z.re * z.im + z.im * z.re;
        z3_re = r * z.re - az3 * z.im;
        z3_im = r * z.im + az3 * z.re;
        az3 = az * aa;
        aa = (2.0 + (double)id) * ((3.0 + (double)id) + (double)id);
        d2 = ((3.0 - (double)id) - (double)id) * (4.0 - (double)id);
        if (aa < d2) {
          ad = aa;
        } else {
          ad = d2;
        }

        ak = 24.0 + 9.0 * (double)id;
        bk = 30.0 - 9.0 * (double)id;
        i = 1;
        exitg1 = false;
        while ((!exitg1) && (i < 26)) {
          if (z3_im == 0.0) {
            b_z3_re = z3_re / aa;
            b_z3_im = 0.0;
          } else if (z3_re == 0.0) {
            b_z3_re = 0.0;
            b_z3_im = z3_im / aa;
          } else {
            b_z3_re = z3_re / aa;
            b_z3_im = z3_im / aa;
          }

          r = trm1.re;
          trm1.re = trm1.re * b_z3_re - trm1.im * b_z3_im;
          trm1.im = r * b_z3_im + trm1.im * b_z3_re;
          s1_re += trm1.re;
          s1_im += trm1.im;
          if (z3_im == 0.0) {
            b_z3_re = z3_re / d2;
            b_z3_im = 0.0;
          } else if (z3_re == 0.0) {
            b_z3_re = 0.0;
            b_z3_im = z3_im / d2;
          } else {
            b_z3_re = z3_re / d2;
            b_z3_im = z3_im / d2;
          }

          r = trm2.re;
          trm2.re = trm2.re * b_z3_re - trm2.im * b_z3_im;
          trm2.im = r * b_z3_im + trm2.im * b_z3_re;
          s2.re += trm2.re;
          s2.im += trm2.im;
          atrm = atrm * az3 / ad;
          aa += ak;
          d2 += bk;
          if ((aa < d2) || rtIsNaN(d2)) {
            ad = aa;
          } else {
            ad = d2;
          }

          if (atrm < 2.2204460492503131E-16 * ad) {
            exitg1 = true;
          } else {
            ak += 18.0;
            bk += 18.0;
            i++;
          }
        }
      }

      if (id == 1) {
        ai.re = -0.25881940379280682 * s2.re;
        ai.im = -0.25881940379280682 * s2.im;
        if (az > 2.2204460492503131E-16) {
          r = z.re * z.re - z.im * z.im;
          az3 = z.re * z.im + z.im * z.re;
          ai.re += (r * s1_re - az3 * s1_im) * 0.17751402694390861;
          ai.im += (r * s1_im + az3 * s1_re) * 0.17751402694390861;
        }

        if (kode == 1) {
        } else {
          s2 = z;
          b_sqrt(&s2);
          r = s2.re * z.im + s2.im * z.re;
          s2.re = 0.66666666666666663 * (s2.re * z.re - s2.im * z.im);
          s2.im = 0.66666666666666663 * r;
          if (s2.im == 0.0) {
            s2.re = exp(s2.re);
            s2.im = 0.0;
          } else if (rtIsInf(s2.im) && rtIsInf(s2.re) && (s2.re < 0.0)) {
            s2.re = 0.0;
            s2.im = 0.0;
          } else {
            r = exp(s2.re / 2.0);
            s2.re = r * (r * cos(s2.im));
            s2.im = r * (r * sin(s2.im));
          }

          r = ai.re;
          ai.re = ai.re * s2.re - ai.im * s2.im;
          ai.im = r * s2.im + ai.im * s2.re;
        }
      } else {
        ai.re = s1_re * 0.35502805388781722 - (z.re * s2.re - z.im * s2.im) *
          0.25881940379280682;
        ai.im = s1_im * 0.35502805388781722 - (z.re * s2.im + z.im * s2.re) *
          0.25881940379280682;
        if (kode == 1) {
        } else {
          s2 = z;
          b_sqrt(&s2);
          r = s2.re * z.im + s2.im * z.re;
          s2.re = 0.66666666666666663 * (s2.re * z.re - s2.im * z.im);
          s2.im = 0.66666666666666663 * r;
          if (s2.im == 0.0) {
            s2.re = exp(s2.re);
            s2.im = 0.0;
          } else if (rtIsInf(s2.im) && rtIsInf(s2.re) && (s2.re < 0.0)) {
            s2.re = 0.0;
            s2.im = 0.0;
          } else {
            r = exp(s2.re / 2.0);
            s2.re = r * (r * cos(s2.im));
            s2.im = r * (r * sin(s2.im));
          }

          r = ai.re;
          ai.re = ai.re * s2.re - ai.im * s2.im;
          ai.im = r * s2.im + ai.im * s2.re;
        }
      }
    }
  }

  return ai;
}

/*
 * File trailer for cairy.c
 *
 * [EOF]
 */
