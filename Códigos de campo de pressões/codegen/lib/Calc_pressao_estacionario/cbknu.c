/*
 * File: cbknu.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao_estacionario.h"
#include "cbknu.h"
#include "exp.h"
#include "cuchk.h"
#include "cmlri.h"
#include "log1.h"
#include "cosh.h"
#include "gammaln.h"
#include "sinh.h"
#include "fix.h"
#include "sqrt1.h"
#include "Calc_pressao_estacionario_rtwutil.h"

/* Function Declarations */
static int b_ckscl(const creal_T zr, int nin, creal_T y[2], double ascle, double
                   tol, double elim);
static creal_T calccs(const creal_T zr, const creal_T s1, double tol);
static int ckscl(const creal_T zr, creal_T *y, double elim);

/* Function Definitions */

/*
 * Arguments    : const creal_T zr
 *                int nin
 *                creal_T y[2]
 *                double ascle
 *                double tol
 *                double elim
 * Return Type  : int
 */
static int b_ckscl(const creal_T zr, int nin, creal_T y[2], double ascle, double
                   tol, double elim)
{
  int nz;
  int n;
  int ic;
  creal_T s1;
  if (nin < 2) {
    n = nin;
  } else {
    n = 2;
  }

  ic = 0;
  nz = 0;
  if (n >= 1) {
    s1 = y[0];
    nz = 1;
    y[0].re = 0.0;
    y[0].im = 0.0;
    if (-zr.re + log(rt_hypotd_snf(s1.re, s1.im)) >= -elim) {
      s1 = calccs(zr, s1, tol);
      if (cuchk(s1, ascle, tol) == 0) {
        y[0] = s1;
        nz = 0;
        ic = 1;
      }
    }
  }

  if (n == 2) {
    if (nz == 1) {
      nz = 2;
    } else {
      nz = 1;
    }

    s1 = y[1];
    y[1].re = 0.0;
    y[1].im = 0.0;
    if (-zr.re + log(rt_hypotd_snf(s1.re, s1.im)) >= -elim) {
      s1 = calccs(zr, s1, tol);
      if (cuchk(s1, ascle, tol) == 0) {
        y[1] = s1;
        nz = 0;
        ic = 2;
      }
    }

    if (ic < 2) {
      y[0].re = 0.0;
      y[0].im = 0.0;
      nz = 2;
    }
  }

  return nz;
}

/*
 * Arguments    : const creal_T zr
 *                const creal_T s1
 *                double tol
 * Return Type  : creal_T
 */
static creal_T calccs(const creal_T zr, const creal_T s1, double tol)
{
  creal_T cs;
  double tmp_re;
  double tmp_im;
  tmp_re = s1.re;
  tmp_im = s1.im;
  if ((s1.im == 0.0) && rtIsNaN(s1.re)) {
  } else if ((fabs(s1.re) > 8.9884656743115785E+307) || (fabs(s1.im) >
              8.9884656743115785E+307)) {
    tmp_re = log(rt_hypotd_snf(s1.re / 2.0, s1.im / 2.0)) + 0.69314718055994529;
    tmp_im = rt_atan2d_snf(s1.im, s1.re);
  } else {
    tmp_re = log(rt_hypotd_snf(s1.re, s1.im));
    tmp_im = rt_atan2d_snf(s1.im, s1.re);
  }

  cs.re = -zr.re + tmp_re;
  cs.im = -zr.im + tmp_im;
  cs.re = exp(cs.re) / tol * cos(cs.im);
  cs.im = exp(-zr.re + tmp_re) / tol * sin(cs.im);
  return cs;
}

/*
 * Arguments    : const creal_T zr
 *                creal_T *y
 *                double elim
 * Return Type  : int
 */
static int ckscl(const creal_T zr, creal_T *y, double elim)
{
  int nz;
  creal_T s1;
  s1 = *y;
  nz = 1;
  y->re = 0.0;
  y->im = 0.0;
  if (-zr.re + log(rt_hypotd_snf(s1.re, s1.im)) >= -elim) {
    s1 = calccs(zr, s1, 2.2204460492503131E-16);
    if (cuchk(s1, 1.0020841800044864E-289, 2.2204460492503131E-16) == 0) {
      *y = s1;
      nz = 0;
    }
  }

  return nz;
}

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int nin
 *                creal_T y[2]
 *                double tol
 *                double elim
 *                double alim
 * Return Type  : int
 */
int b_cbknu(const creal_T z, double fnu, int kode, int nin, creal_T y[2], double
            tol, double elim, double alim)
{
  int nz;
  int n;
  double yy;
  double caz;
  double cssr[3];
  double csrr[3];
  double bry[3];
  int iflag;
  double rz_re;
  double fks;
  double rz_im;
  double etest;
  int inu;
  double a1;
  double dnu;
  double rk;
  boolean_T goto_mw110;
  double dnu2;
  double ak;
  double fhs;
  double s1_re;
  double s1_im;
  double s2_re;
  double s2_im;
  creal_T zd;
  double ck_re;
  double ck_im;
  int inub;
  boolean_T goto_mw225;
  boolean_T goto_mw240;
  boolean_T goto_mw270;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T goto_mw210;
  creal_T cs;
  creal_T p2;
  double coef_re;
  creal_T fmu;
  double coef_im;
  int kflag;
  double a2;
  double t1;
  int i;
  boolean_T exitg1;
  creal_T f;
  static const double dv4[8] = { 0.57721566490153287, -0.042002635034095237,
    -0.042197734555544333, 0.0072189432466631, -0.00021524167411495098,
    -2.0134854780788239E-5, 1.1330272319816959E-6, 6.1160951044814161E-9 };

  creal_T p1;
  int k;
  int j;
  creal_T cy[2];
  double b_etest;
  boolean_T b_guard1 = false;
  if (nin < 2) {
    n = nin;
  } else {
    n = 2;
  }

  yy = z.im;
  caz = rt_hypotd_snf(z.re, z.im);
  cssr[0] = 1.0 / tol;
  cssr[1] = 1.0;
  cssr[2] = tol;
  csrr[0] = tol;
  csrr[1] = 1.0;
  csrr[2] = 1.0 / tol;
  bry[0] = 2.2250738585072014E-305 / tol;
  bry[1] = tol / 2.2250738585072014E-305;
  bry[2] = 1.7976931348623157E+308;
  nz = 0;
  iflag = 0;
  if (z.im == 0.0) {
    rz_re = 2.0 / z.re;
    rz_im = 0.0;
  } else if (z.re == 0.0) {
    rz_re = 0.0;
    rz_im = -(2.0 / z.im);
  } else {
    fks = fabs(z.re);
    etest = fabs(z.im);
    if (fks > etest) {
      a1 = z.im / z.re;
      rk = z.re + a1 * z.im;
      rz_re = (2.0 + a1 * 0.0) / rk;
      rz_im = (0.0 - a1 * 2.0) / rk;
    } else if (etest == fks) {
      if (z.re > 0.0) {
        rk = 0.5;
      } else {
        rk = -0.5;
      }

      if (z.im > 0.0) {
        a1 = 0.5;
      } else {
        a1 = -0.5;
      }

      rz_re = (2.0 * rk + 0.0 * a1) / fks;
      rz_im = (0.0 * rk - 2.0 * a1) / fks;
    } else {
      a1 = z.re / z.im;
      rk = z.im + a1 * z.re;
      rz_re = a1 * 2.0 / rk;
      rz_im = (a1 * 0.0 - 2.0) / rk;
    }
  }

  inu = (int)(fnu + 0.5);
  dnu = fnu - (double)inu;
  goto_mw110 = (fabs(dnu) == 0.5);
  if ((!goto_mw110) && (fabs(dnu) > tol)) {
    dnu2 = dnu * dnu;
  } else {
    dnu2 = 0.0;
  }

  ak = 1.0;
  fhs = 0.0;
  s1_re = 0.0;
  s1_im = 0.0;
  s2_re = 0.0;
  s2_im = 0.0;
  zd.re = 0.0;
  zd.im = 0.0;
  ck_re = 1.0;
  ck_im = 0.0;
  inub = 1;
  goto_mw225 = false;
  goto_mw240 = false;
  goto_mw270 = false;
  if (goto_mw110 || (caz > 2.0)) {
    goto_mw110 = true;
  } else {
    goto_mw110 = false;
  }

  guard1 = false;
  guard2 = false;
  if (!goto_mw110) {
    rk = 1.0;
    p2.re = rz_re;
    p2.im = rz_im;
    b_log(&p2);
    fmu.re = dnu * p2.re;
    fmu.im = dnu * p2.im;
    if (dnu != 0.0) {
      rk = dnu * 3.1415926535897931;
      rk /= sin(rk);
      cs = fmu;
      b_sinh(&cs);
      p2.re = 1.0 / dnu * cs.re;
      p2.im = 1.0 / dnu * cs.im;
    }

    etest = 1.0 + dnu;
    gammaln(&etest);
    a2 = exp(-etest);
    t1 = 1.0 / (a2 * rk);
    if (fabs(dnu) > 0.1) {
      etest = (t1 - a2) / (dnu + dnu);
    } else {
      a1 = 0.57721566490153287;
      i = 2;
      exitg1 = false;
      while ((!exitg1) && (i < 9)) {
        ak *= dnu2;
        s1_re = dv4[i - 1] * ak;
        a1 += s1_re;
        if (fabs(s1_re) < tol) {
          exitg1 = true;
        } else {
          i++;
        }
      }

      etest = -a1;
    }

    etest *= rk;
    cs = fmu;
    b_cosh(&cs);
    f.re = p2.re * (0.5 * (t1 + a2) * rk) + etest * cs.re;
    f.im = p2.im * (0.5 * (t1 + a2) * rk) + etest * cs.im;
    b_exp(&fmu);
    p1.re = 0.5 / a2 * fmu.re;
    p1.im = 0.5 / a2 * fmu.im;
    ak = 0.5 / t1;
    if (fmu.im == 0.0) {
      cs.re = ak / fmu.re;
      cs.im = 0.0;
    } else if (fmu.re == 0.0) {
      if (ak == 0.0) {
        cs.re = 0.0 / fmu.im;
        cs.im = 0.0;
      } else {
        cs.re = 0.0;
        cs.im = -(ak / fmu.im);
      }
    } else {
      fks = fabs(fmu.re);
      etest = fabs(fmu.im);
      if (fks > etest) {
        a1 = fmu.im / fmu.re;
        rk = fmu.re + a1 * fmu.im;
        cs.re = (ak + a1 * 0.0) / rk;
        cs.im = (0.0 - a1 * ak) / rk;
      } else if (etest == fks) {
        if (fmu.re > 0.0) {
          rk = 0.5;
        } else {
          rk = -0.5;
        }

        if (fmu.im > 0.0) {
          a1 = 0.5;
        } else {
          a1 = -0.5;
        }

        cs.re = (ak * rk + 0.0 * a1) / fks;
        cs.im = (0.0 * rk - ak * a1) / fks;
      } else {
        a1 = fmu.re / fmu.im;
        rk = fmu.im + a1 * fmu.re;
        cs.re = a1 * ak / rk;
        cs.im = (a1 * 0.0 - ak) / rk;
      }
    }

    s1_re = f.re;
    s1_im = f.im;
    s2_re = p1.re;
    s2_im = p1.im;
    ak = 1.0;
    a1 = 1.0;
    a2 = 1.0 - dnu2;
    if ((inu > 0) || (n > 1)) {
      goto_mw110 = true;
    } else {
      goto_mw110 = false;
    }

    if (!goto_mw110) {
      if (caz >= tol) {
        coef_re = 0.25 * (z.re * z.re - z.im * z.im);
        coef_im = 0.25 * (z.re * z.im + z.im * z.re);
        t1 = 0.25 * caz * caz;
        do {
          f.re *= ak;
          f.im *= ak;
          f.re += p1.re;
          f.im += p1.im;
          f.re += cs.re;
          f.im += cs.im;
          f.re *= 1.0 / a2;
          f.im *= 1.0 / a2;
          p1.re *= 1.0 / (ak - dnu);
          p1.im *= 1.0 / (ak - dnu);
          cs.re *= 1.0 / (ak + dnu);
          cs.im *= 1.0 / (ak + dnu);
          etest = ck_re;
          ck_re = ck_re * coef_re - ck_im * coef_im;
          ck_im = etest * coef_im + ck_im * coef_re;
          ck_re *= 1.0 / ak;
          ck_im *= 1.0 / ak;
          s1_re += ck_re * f.re - ck_im * f.im;
          s1_im += ck_re * f.im + ck_im * f.re;
          a1 = a1 * t1 / ak;
          a2 = ((a2 + ak) + ak) + 1.0;
          ak++;
        } while (!!(a1 > tol));
      }

      y[0].re = s1_re;
      y[0].im = s1_im;
      if (kode != 1) {
        y[0] = z;
        b_exp(&y[0]);
        etest = y[0].re;
        rk = y[0].im;
        y[0].re = etest * s1_re - rk * s1_im;
        y[0].im = etest * s1_im + rk * s1_re;
      }
    } else {
      if (caz >= tol) {
        coef_re = 0.25 * (z.re * z.re - z.im * z.im);
        coef_im = 0.25 * (z.re * z.im + z.im * z.re);
        t1 = 0.25 * caz * caz;
        do {
          f.re *= ak;
          f.im *= ak;
          f.re += p1.re;
          f.im += p1.im;
          f.re += cs.re;
          f.im += cs.im;
          f.re *= 1.0 / a2;
          f.im *= 1.0 / a2;
          p1.re *= 1.0 / (ak - dnu);
          p1.im *= 1.0 / (ak - dnu);
          cs.re *= 1.0 / (ak + dnu);
          cs.im *= 1.0 / (ak + dnu);
          etest = ck_re;
          ck_re = ck_re * coef_re - ck_im * coef_im;
          ck_im = etest * coef_im + ck_im * coef_re;
          ck_re *= 1.0 / ak;
          ck_im *= 1.0 / ak;
          s1_re += ck_re * f.re - ck_im * f.im;
          s1_im += ck_re * f.im + ck_im * f.re;
          etest = p1.re - f.re * ak;
          rk = p1.im - f.im * ak;
          s2_re += etest * ck_re - rk * ck_im;
          s2_im += etest * ck_im + rk * ck_re;
          a1 = a1 * t1 / ak;
          a2 = ((a2 + ak) + ak) + 1.0;
          ak++;
        } while (!!(a1 > tol));
      }

      kflag = 1;
      if ((fnu + 1.0) * fabs(p2.re) > alim) {
        kflag = 2;
      }

      etest = cssr[kflag] * s2_re;
      s2_im *= cssr[kflag];
      s2_re = etest * rz_re - s2_im * rz_im;
      s2_im = etest * rz_im + s2_im * rz_re;
      s1_re *= cssr[kflag];
      s1_im *= cssr[kflag];
      if (kode != 1) {
        f = z;
        b_exp(&f);
        etest = s1_re;
        s1_re = s1_re * f.re - s1_im * f.im;
        s1_im = etest * f.im + s1_im * f.re;
        etest = s2_re;
        s2_re = s2_re * f.re - s2_im * f.im;
        s2_im = etest * f.im + s2_im * f.re;
      }

      goto_mw210 = true;
      guard1 = true;
    }
  } else {
    goto_mw210 = false;
    cs = z;
    b_sqrt(&cs);
    if (cs.im == 0.0) {
      coef_re = 1.2533141373155001 / cs.re;
      coef_im = 0.0;
    } else if (cs.re == 0.0) {
      coef_re = 0.0;
      coef_im = -(1.2533141373155001 / cs.im);
    } else {
      fks = fabs(cs.re);
      etest = fabs(cs.im);
      if (fks > etest) {
        a1 = cs.im / cs.re;
        rk = cs.re + a1 * cs.im;
        coef_re = (1.2533141373155001 + a1 * 0.0) / rk;
        coef_im = (0.0 - a1 * 1.2533141373155001) / rk;
      } else if (etest == fks) {
        if (cs.re > 0.0) {
          rk = 0.5;
        } else {
          rk = -0.5;
        }

        if (cs.im > 0.0) {
          a1 = 0.5;
        } else {
          a1 = -0.5;
        }

        coef_re = (1.2533141373155001 * rk + 0.0 * a1) / fks;
        coef_im = (0.0 * rk - 1.2533141373155001 * a1) / fks;
      } else {
        a1 = cs.re / cs.im;
        rk = cs.im + a1 * cs.re;
        coef_re = a1 * 1.2533141373155001 / rk;
        coef_im = (a1 * 0.0 - 1.2533141373155001) / rk;
      }
    }

    kflag = 1;
    if (!(kode == 2)) {
      if (z.re > alim) {
        iflag = 1;
      } else {
        a1 = exp(-z.re);
        fmu.re = cos(z.im) * a1;
        fmu.im = sin(z.im) * -a1;
        a2 = coef_re;
        coef_re = coef_re * fmu.re - coef_im * fmu.im;
        coef_im = a2 * fmu.im + coef_im * fmu.re;
      }
    }

    if (fabs(dnu) == 0.5) {
      s1_re = coef_re;
      s1_im = coef_im;
      s2_re = coef_re;
      s2_im = coef_im;
      goto_mw210 = true;
    } else {
      ak = fabs(cos(3.1415926535897931 * dnu));
      if (ak == 0.0) {
        s1_re = coef_re;
        s1_im = coef_im;
        s2_re = coef_re;
        s2_im = coef_im;
        goto_mw210 = true;
      } else {
        fhs = fabs(0.25 - dnu2);
        if (fhs == 0.0) {
          s1_re = coef_re;
          s1_im = coef_im;
          s2_re = coef_re;
          s2_im = coef_im;
          goto_mw210 = true;
        }
      }
    }

    if (!goto_mw210) {
      if (z.re != 0.0) {
        t1 = fabs(atan(z.im / z.re));
      } else {
        t1 = 1.5707963267948966;
      }

      if (28.66666665740641 > caz) {
        ak = 1.8976999933151775 * ak / (tol * sqrt(sqrt(caz)));
        ak = (log(ak) + caz * cos(3.0 * t1 / (1.0 + caz)) / (1.0 + 0.008 * caz))
          / cos(14.7 * t1 / (28.0 + caz));
        ak = 0.12125 * ak * ak / caz + 1.5;
        guard2 = true;
      } else {
        etest = ak / (3.1415926535897931 * caz * tol);
        ak = 1.0;
        if (etest >= 1.0) {
          fks = 2.0;
          rk = (caz + caz) + 2.0;
          a1 = 0.0;
          a2 = 1.0;
          goto_mw110 = true;
          i = 1;
          exitg1 = false;
          while ((!exitg1) && (i < 31)) {
            s1_re = a2;
            a2 = rk / (ak + 1.0) * a2 - fhs / fks * a1;
            a1 = s1_re;
            rk += 2.0;
            fks = ((fks + ak) + ak) + 2.0;
            fhs = (fhs + ak) + ak;
            ak++;
            if (etest < fabs(a2) * ak) {
              goto_mw110 = false;
              exitg1 = true;
            } else {
              i++;
            }
          }

          if (goto_mw110) {
            nz = -2;
          } else {
            ak += 1.909859317102744 * t1 * sqrt(28.66666665740641 / caz);
            fhs = fabs(0.25 - dnu2);
            guard2 = true;
          }
        } else {
          guard2 = true;
        }
      }
    } else {
      guard1 = true;
    }
  }

  if (guard2) {
    b_fix(&ak);
    k = (int)ak;
    ak = k;
    fks = (double)k * (double)k;
    p1.re = 0.0;
    p1.im = 0.0;
    p2.re = tol;
    p2.im = 0.0;
    cs = p2;
    for (i = 1; i <= k; i++) {
      a1 = fks - ak;
      rk = 2.0 / (ak + 1.0);
      fmu = p2;
      a2 = (ak + z.re) * rk;
      etest = yy * rk;
      rk = p2.re;
      p2.re = p2.re * a2 - p2.im * etest;
      p2.im = rk * etest + p2.im * a2;
      p2.re -= p1.re;
      p2.im -= p1.im;
      p2.re *= (fks + ak) / (a1 + fhs);
      p2.im *= (fks + ak) / (a1 + fhs);
      p1 = fmu;
      cs.re += p2.re;
      cs.im += p2.im;
      fks = (a1 - ak) + 1.0;
      ak--;
    }

    fmu.re = 1.0 / rt_hypotd_snf(cs.re, cs.im);
    cs.im = -cs.im;
    etest = cs.re;
    cs.re = cs.re * fmu.re - cs.im * 0.0;
    cs.im = etest * 0.0 + cs.im * fmu.re;
    rk = fmu.re * p2.re - 0.0 * p2.im;
    etest = fmu.re * p2.im + 0.0 * p2.re;
    a2 = coef_re * rk - coef_im * etest;
    coef_im = coef_re * etest + coef_im * rk;
    s1_re = a2 * cs.re - coef_im * cs.im;
    s1_im = a2 * cs.im + coef_im * cs.re;
    if ((inu > 0) || (n > 1)) {
      goto_mw110 = true;
    } else {
      goto_mw110 = false;
    }

    if (!goto_mw110) {
      zd = z;
      if (iflag == 1) {
        goto_mw270 = true;
      } else {
        goto_mw240 = true;
      }
    } else {
      fmu.re = 1.0 / rt_hypotd_snf(p2.re, p2.im);
      etest = p1.re;
      p1.re = p1.re * fmu.re - p1.im * 0.0;
      p1.im = etest * 0.0 + p1.im * fmu.re;
      p2.im = -p2.im;
      rk = p2.re;
      p2.re = p2.re * fmu.re - p2.im * 0.0;
      p2.im = rk * 0.0 + p2.im * fmu.re;
      rk = p1.re * p2.im + p1.im * p2.re;
      ak = (dnu + 0.5) - (p1.re * p2.re - p1.im * p2.im);
      a2 = 0.0 - (p1.re * p2.im + p1.im * p2.re);
      if (z.im == 0.0) {
        if (0.0 - rk == 0.0) {
          etest = ak / z.re;
          rk = 0.0;
        } else if (ak == 0.0) {
          etest = 0.0;
          rk = (0.0 - rk) / z.re;
        } else {
          etest = ak / z.re;
          rk = (0.0 - rk) / z.re;
        }
      } else if (z.re == 0.0) {
        if (ak == 0.0) {
          etest = (0.0 - rk) / z.im;
          rk = 0.0;
        } else if (0.0 - rk == 0.0) {
          etest = 0.0;
          rk = -(ak / z.im);
        } else {
          etest = (0.0 - rk) / z.im;
          rk = -(ak / z.im);
        }
      } else {
        fks = fabs(z.re);
        etest = fabs(z.im);
        if (fks > etest) {
          a1 = z.im / z.re;
          rk = z.re + a1 * z.im;
          etest = (ak + a1 * a2) / rk;
          rk = (a2 - a1 * ak) / rk;
        } else if (etest == fks) {
          if (z.re > 0.0) {
            rk = 0.5;
          } else {
            rk = -0.5;
          }

          if (z.im > 0.0) {
            a1 = 0.5;
          } else {
            a1 = -0.5;
          }

          etest = (ak * rk + a2 * a1) / fks;
          rk = (a2 * rk - ak * a1) / fks;
        } else {
          a1 = z.re / z.im;
          rk = z.im + a1 * z.re;
          etest = (a1 * ak + a2) / rk;
          rk = (a1 * a2 - ak) / rk;
        }
      }

      etest++;
      s2_re = s1_re * etest - s1_im * rk;
      s2_im = s1_re * rk + s1_im * etest;
      goto_mw210 = true;
    }

    guard1 = true;
  }

  if (guard1) {
    if (goto_mw240 || goto_mw270) {
    } else if (goto_mw210) {
      ck_re = (dnu + 1.0) * rz_re;
      ck_im = (dnu + 1.0) * rz_im;
      if (n == 1) {
        inu--;
      }

      if (inu > 0) {
        if (iflag == 1) {
          a2 = 0.5 * elim;
          coef_re = exp(-elim);
          zd = z;
          rk = z.re;
          k = inu;
          j = 1;
          for (i = 0; i < 2; i++) {
            cy[i].re = 0.0;
            cy[i].im = 0.0;
          }

          i = 1;
          exitg1 = false;
          while ((!exitg1) && (i <= inu)) {
            f.re = s2_re;
            f.im = s2_im;
            etest = s2_re;
            s2_re = s2_re * ck_re - s2_im * ck_im;
            s2_im = etest * ck_im + s2_im * ck_re;
            s2_re += s1_re;
            s2_im += s1_im;
            s1_re = f.re;
            s1_im = f.im;
            ck_re += rz_re;
            ck_im += rz_im;
            etest = log(rt_hypotd_snf(s2_re, s2_im));
            b_guard1 = false;
            if (-rk + etest >= -elim) {
              cs.re = s2_re;
              cs.im = s2_im;
              b_log(&cs);
              p2.re = cs.re + -zd.re;
              p2.im = cs.im + -zd.im;
              p1.re = exp(p2.re) / tol * cos(p2.im);
              p1.im = exp(p2.re) / tol * sin(p2.im);
              if (cuchk(p1, bry[0], tol) == 0) {
                j = 1 - j;
                cy[j] = p1;
                if (k == i - 1) {
                  kflag = 0;
                  inub = i + 1;
                  s2_re = cy[j].re;
                  s2_im = cy[j].im;
                  j = 1 - j;
                  s1_re = cy[j].re;
                  s1_im = cy[j].im;
                  if (i + 1 <= inu) {
                    goto_mw225 = true;
                  } else {
                    if (n == 1) {
                      s1_re = s2_re;
                      s1_im = s2_im;
                    }

                    goto_mw240 = true;
                  }

                  exitg1 = true;
                } else {
                  k = i;
                  i++;
                }
              } else {
                b_guard1 = true;
              }
            } else {
              b_guard1 = true;
            }

            if (b_guard1) {
              if (etest >= a2) {
                rk -= elim;
                s1_re = f.re * coef_re - f.im * 0.0;
                s1_im = f.re * 0.0 + f.im * coef_re;
                etest = s2_re;
                s2_re = s2_re * coef_re - s2_im * 0.0;
                s2_im = etest * 0.0 + s2_im * coef_re;
                zd.re = rk;
                zd.im = yy;
              }

              i++;
            }
          }

          if (goto_mw225 || goto_mw240) {
            goto_mw110 = true;
          } else {
            goto_mw110 = false;
          }

          if (!goto_mw110) {
            if (n == 1) {
              s1_re = s2_re;
              s1_im = s2_im;
            }

            goto_mw270 = true;
          }
        } else {
          goto_mw225 = true;
        }
      }

      if (goto_mw225 || goto_mw240 || goto_mw270) {
        goto_mw110 = true;
      } else {
        goto_mw110 = false;
      }

      if (!goto_mw110) {
        if (n == 1) {
          s1_re = s2_re;
          s1_im = s2_im;
        }

        zd = z;
        if (iflag != 1) {
          goto_mw240 = true;
        }
      }
    } else {
      goto_mw225 = true;
    }

    if (goto_mw225 || goto_mw240) {
      if (goto_mw225) {
        p1.re = csrr[kflag];
        a2 = bry[kflag];
        while (inub <= inu) {
          f.re = s2_re;
          f.im = s2_im;
          etest = s2_re;
          s2_re = s2_re * ck_re - s2_im * ck_im;
          s2_im = etest * ck_im + s2_im * ck_re;
          s2_re += s1_re;
          s2_im += s1_im;
          s1_re = f.re;
          s1_im = f.im;
          ck_re += rz_re;
          ck_im += rz_im;
          if (kflag + 1 < 3) {
            p2.re = s2_re * p1.re - s2_im * 0.0;
            p2.im = s2_re * 0.0 + s2_im * p1.re;
            etest = fabs(p2.re);
            rk = fabs(p2.im);
            if ((etest > rk) || rtIsNaN(rk)) {
              b_etest = etest;
            } else {
              b_etest = rk;
            }

            if (b_etest > a2) {
              kflag++;
              a2 = bry[kflag];
              s1_re = cssr[kflag] * (f.re * p1.re - f.im * 0.0);
              s1_im = cssr[kflag] * (f.re * 0.0 + f.im * p1.re);
              s2_re = cssr[kflag] * p2.re;
              s2_im = cssr[kflag] * p2.im;
              p1.re = csrr[kflag];
            }
          }

          inub++;
        }

        if (n == 1) {
          s1_re = s2_re;
          s1_im = s2_im;
        }
      }

      y[0].re = csrr[kflag] * s1_re;
      y[0].im = csrr[kflag] * s1_im;
      if (n == 2) {
        y[1].re = csrr[kflag] * s2_re;
        y[1].im = csrr[kflag] * s2_im;
      }
    } else {
      y[0].re = s1_re;
      y[0].im = s1_im;
      if (n != 1) {
        y[1].re = s2_re;
        y[1].im = s2_im;
      }

      nz = b_ckscl(zd, n, y, bry[0], tol, elim);
      if (n <= nz) {
      } else {
        y[nz].re *= tol;
        y[nz].im *= tol;
        if ((n >= 0) && (nz < n - MAX_int32_T)) {
          k = MAX_int32_T;
        } else if ((n < 0) && (nz > n - MIN_int32_T)) {
          k = MIN_int32_T;
        } else {
          k = n - nz;
        }

        if (k >= 2) {
          y[nz + 1].re *= tol;
          y[nz + 1].im *= tol;
        }
      }
    }
  }

  return nz;
}

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                double alim
 *                creal_T *y
 *                int *nz
 * Return Type  : void
 */
void cbknu(const creal_T z, double fnu, int kode, double alim, creal_T *y, int
           *nz)
{
  double yy;
  double caz;
  int iflag;
  double rz_re;
  double fks;
  double rz_im;
  double rk;
  int inu;
  double dnu;
  double etest;
  boolean_T goto_mw110;
  double dnu2;
  double tm;
  double ak;
  double fhs;
  double s1_re;
  double s1_im;
  double s2_re;
  double s2_im;
  creal_T zd;
  double ck_re;
  double ck_im;
  int inub;
  boolean_T goto_mw225;
  boolean_T goto_mw240;
  boolean_T goto_mw270;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  double a2;
  creal_T coef;
  creal_T smu;
  creal_T fmu;
  int kflag;
  double t1;
  int i;
  boolean_T exitg1;
  creal_T p1;
  static const double dv0[8] = { 0.57721566490153287, -0.042002635034095237,
    -0.042197734555544333, 0.0072189432466631, -0.00021524167411495098,
    -2.0134854780788239E-5, 1.1330272319816959E-6, 6.1160951044814161E-9 };

  double p2_re;
  double p2_im;
  int k;
  double cs_re;
  double a1;
  static const double dv1[3] = { 2.2204460492503131E-16, 1.0,
    4.503599627370496E+15 };

  static const double dv2[3] = { 1.0020841800044864E-289, 9.9792015476736E+288,
    1.7976931348623157E+308 };

  boolean_T earlyExit;
  int j;
  creal_T cy[2];
  static const double dv3[3] = { 4.503599627370496E+15, 1.0,
    2.2204460492503131E-16 };

  double b_etest;
  boolean_T b_guard1 = false;
  y->re = 0.0;
  y->im = 0.0;
  yy = z.im;
  caz = rt_hypotd_snf(z.re, z.im);
  *nz = 0;
  iflag = 0;
  if (z.im == 0.0) {
    rz_re = 2.0 / z.re;
    rz_im = 0.0;
  } else if (z.re == 0.0) {
    rz_re = 0.0;
    rz_im = -(2.0 / z.im);
  } else {
    fks = fabs(z.re);
    rk = fabs(z.im);
    if (fks > rk) {
      fks = z.im / z.re;
      etest = z.re + fks * z.im;
      rz_re = (2.0 + fks * 0.0) / etest;
      rz_im = (0.0 - fks * 2.0) / etest;
    } else if (rk == fks) {
      if (z.re > 0.0) {
        etest = 0.5;
      } else {
        etest = -0.5;
      }

      if (z.im > 0.0) {
        tm = 0.5;
      } else {
        tm = -0.5;
      }

      rz_re = (2.0 * etest + 0.0 * tm) / fks;
      rz_im = (0.0 * etest - 2.0 * tm) / fks;
    } else {
      fks = z.re / z.im;
      etest = z.im + fks * z.re;
      rz_re = fks * 2.0 / etest;
      rz_im = (fks * 0.0 - 2.0) / etest;
    }
  }

  inu = (int)(fnu + 0.5);
  dnu = fnu - (double)inu;
  goto_mw110 = (fabs(dnu) == 0.5);
  if ((!goto_mw110) && (fabs(dnu) > 2.2204460492503131E-16)) {
    dnu2 = dnu * dnu;
  } else {
    dnu2 = 0.0;
  }

  ak = 1.0;
  fhs = 0.0;
  s1_re = 0.0;
  s1_im = 0.0;
  s2_re = 0.0;
  s2_im = 0.0;
  zd.re = 0.0;
  zd.im = 0.0;
  ck_re = 1.0;
  ck_im = 0.0;
  inub = 1;
  goto_mw225 = false;
  goto_mw240 = false;
  goto_mw270 = false;
  if (goto_mw110 || (caz > 2.0)) {
    goto_mw110 = true;
  } else {
    goto_mw110 = false;
  }

  guard1 = false;
  guard2 = false;
  if (!goto_mw110) {
    a2 = 1.0;
    smu.re = rz_re;
    smu.im = rz_im;
    b_log(&smu);
    fmu.re = dnu * smu.re;
    fmu.im = dnu * smu.im;
    if (dnu != 0.0) {
      a2 = dnu * 3.1415926535897931;
      a2 /= sin(a2);
      coef = fmu;
      b_sinh(&coef);
      smu.re = 1.0 / dnu * coef.re;
      smu.im = 1.0 / dnu * coef.im;
    }

    etest = 1.0 + dnu;
    gammaln(&etest);
    rk = exp(-etest);
    t1 = 1.0 / (rk * a2);
    if (fabs(dnu) > 0.1) {
      etest = (t1 - rk) / (dnu + dnu);
    } else {
      fks = 0.57721566490153287;
      i = 2;
      exitg1 = false;
      while ((!exitg1) && (i < 9)) {
        ak *= dnu2;
        tm = dv0[i - 1] * ak;
        fks += tm;
        if (fabs(tm) < 2.2204460492503131E-16) {
          exitg1 = true;
        } else {
          i++;
        }
      }

      etest = -fks;
    }

    etest *= a2;
    coef = fmu;
    b_cosh(&coef);
    p1.re = smu.re * (0.5 * (t1 + rk) * a2) + etest * coef.re;
    p1.im = smu.im * (0.5 * (t1 + rk) * a2) + etest * coef.im;
    b_exp(&fmu);
    p2_re = 0.5 / rk * fmu.re;
    p2_im = 0.5 / rk * fmu.im;
    ak = 0.5 / t1;
    if (fmu.im == 0.0) {
      cs_re = ak / fmu.re;
      tm = 0.0;
    } else if (fmu.re == 0.0) {
      if (ak == 0.0) {
        cs_re = 0.0 / fmu.im;
        tm = 0.0;
      } else {
        cs_re = 0.0;
        tm = -(ak / fmu.im);
      }
    } else {
      fks = fabs(fmu.re);
      rk = fabs(fmu.im);
      if (fks > rk) {
        fks = fmu.im / fmu.re;
        etest = fmu.re + fks * fmu.im;
        cs_re = (ak + fks * 0.0) / etest;
        tm = (0.0 - fks * ak) / etest;
      } else if (rk == fks) {
        if (fmu.re > 0.0) {
          etest = 0.5;
        } else {
          etest = -0.5;
        }

        if (fmu.im > 0.0) {
          tm = 0.5;
        } else {
          tm = -0.5;
        }

        cs_re = (ak * etest + 0.0 * tm) / fks;
        tm = (0.0 * etest - ak * tm) / fks;
      } else {
        fks = fmu.re / fmu.im;
        etest = fmu.im + fks * fmu.re;
        cs_re = fks * ak / etest;
        tm = (fks * 0.0 - ak) / etest;
      }
    }

    s1_re = p1.re;
    s1_im = p1.im;
    s2_re = p2_re;
    s2_im = p2_im;
    ak = 1.0;
    a1 = 1.0;
    etest = 1.0 - dnu2;
    if (!(inu > 0)) {
      if (caz >= 2.2204460492503131E-16) {
        coef.re = 0.25 * (z.re * z.re - z.im * z.im);
        coef.im = 0.25 * (z.re * z.im + z.im * z.re);
        t1 = 0.25 * caz * caz;
        do {
          p1.re *= ak;
          p1.im *= ak;
          p1.re += p2_re;
          p1.im += p2_im;
          p1.re += cs_re;
          p1.im += tm;
          p1.re *= 1.0 / etest;
          p1.im *= 1.0 / etest;
          p2_re *= 1.0 / (ak - dnu);
          p2_im *= 1.0 / (ak - dnu);
          cs_re *= 1.0 / (ak + dnu);
          tm *= 1.0 / (ak + dnu);
          rk = ck_re;
          ck_re = ck_re * coef.re - ck_im * coef.im;
          ck_im = rk * coef.im + ck_im * coef.re;
          ck_re *= 1.0 / ak;
          ck_im *= 1.0 / ak;
          s1_re += ck_re * p1.re - ck_im * p1.im;
          s1_im += ck_re * p1.im + ck_im * p1.re;
          a1 = a1 * t1 / ak;
          etest = ((etest + ak) + ak) + 1.0;
          ak++;
        } while (!!(a1 > 2.2204460492503131E-16));
      }

      y->re = s1_re;
      y->im = s1_im;
      if (kode != 1) {
        coef = z;
        b_exp(&coef);
        y->re = coef.re * s1_re - coef.im * s1_im;
        y->im = coef.re * s1_im + coef.im * s1_re;
      }
    } else {
      if (caz >= 2.2204460492503131E-16) {
        coef.re = 0.25 * (z.re * z.re - z.im * z.im);
        coef.im = 0.25 * (z.re * z.im + z.im * z.re);
        t1 = 0.25 * caz * caz;
        do {
          p1.re *= ak;
          p1.im *= ak;
          p1.re += p2_re;
          p1.im += p2_im;
          p1.re += cs_re;
          p1.im += tm;
          p1.re *= 1.0 / etest;
          p1.im *= 1.0 / etest;
          p2_re *= 1.0 / (ak - dnu);
          p2_im *= 1.0 / (ak - dnu);
          cs_re *= 1.0 / (ak + dnu);
          tm *= 1.0 / (ak + dnu);
          rk = ck_re;
          ck_re = ck_re * coef.re - ck_im * coef.im;
          ck_im = rk * coef.im + ck_im * coef.re;
          ck_re *= 1.0 / ak;
          ck_im *= 1.0 / ak;
          s1_re += ck_re * p1.re - ck_im * p1.im;
          s1_im += ck_re * p1.im + ck_im * p1.re;
          a2 = p2_re - p1.re * ak;
          rk = p2_im - p1.im * ak;
          s2_re += a2 * ck_re - rk * ck_im;
          s2_im += a2 * ck_im + rk * ck_re;
          a1 = a1 * t1 / ak;
          etest = ((etest + ak) + ak) + 1.0;
          ak++;
        } while (!!(a1 > 2.2204460492503131E-16));
      }

      kflag = 1;
      if ((fnu + 1.0) * fabs(smu.re) > alim) {
        kflag = 2;
      }

      etest = dv3[kflag] * s2_re;
      s2_im *= dv3[kflag];
      s2_re = etest * rz_re - s2_im * rz_im;
      s2_im = etest * rz_im + s2_im * rz_re;
      s1_re *= dv3[kflag];
      s1_im *= dv3[kflag];
      if (kode != 1) {
        p1 = z;
        b_exp(&p1);
        etest = s1_re;
        s1_re = s1_re * p1.re - s1_im * p1.im;
        s1_im = etest * p1.im + s1_im * p1.re;
        etest = s2_re;
        s2_re = s2_re * p1.re - s2_im * p1.im;
        s2_im = etest * p1.im + s2_im * p1.re;
      }

      goto_mw110 = true;
      guard1 = true;
    }
  } else {
    goto_mw110 = false;
    coef = z;
    b_sqrt(&coef);
    if (coef.im == 0.0) {
      coef.re = 1.2533141373155001 / coef.re;
      coef.im = 0.0;
    } else if (coef.re == 0.0) {
      coef.re = 0.0;
      coef.im = -(1.2533141373155001 / coef.im);
    } else {
      fks = fabs(coef.re);
      rk = fabs(coef.im);
      if (fks > rk) {
        fks = coef.im / coef.re;
        etest = coef.re + fks * coef.im;
        coef.re = (1.2533141373155001 + fks * 0.0) / etest;
        coef.im = (0.0 - fks * 1.2533141373155001) / etest;
      } else if (rk == fks) {
        if (coef.re > 0.0) {
          etest = 0.5;
        } else {
          etest = -0.5;
        }

        if (coef.im > 0.0) {
          tm = 0.5;
        } else {
          tm = -0.5;
        }

        coef.re = (1.2533141373155001 * etest + 0.0 * tm) / fks;
        coef.im = (0.0 * etest - 1.2533141373155001 * tm) / fks;
      } else {
        fks = coef.re / coef.im;
        etest = coef.im + fks * coef.re;
        coef.re = fks * 1.2533141373155001 / etest;
        coef.im = (fks * 0.0 - 1.2533141373155001) / etest;
      }
    }

    kflag = 1;
    if (!(kode == 2)) {
      if (z.re > alim) {
        iflag = 1;
      } else {
        etest = exp(-z.re);
        fmu.re = cos(z.im) * etest;
        fmu.im = sin(z.im) * -etest;
        rk = coef.re;
        coef.re = coef.re * fmu.re - coef.im * fmu.im;
        coef.im = rk * fmu.im + coef.im * fmu.re;
      }
    }

    if (fabs(dnu) == 0.5) {
      s1_re = coef.re;
      s1_im = coef.im;
      s2_re = coef.re;
      s2_im = coef.im;
      goto_mw110 = true;
    } else {
      ak = fabs(cos(3.1415926535897931 * dnu));
      if (ak == 0.0) {
        s1_re = coef.re;
        s1_im = coef.im;
        s2_re = coef.re;
        s2_im = coef.im;
        goto_mw110 = true;
      } else {
        fhs = fabs(0.25 - dnu2);
        if (fhs == 0.0) {
          s1_re = coef.re;
          s1_im = coef.im;
          s2_re = coef.re;
          s2_im = coef.im;
          goto_mw110 = true;
        }
      }
    }

    if (!goto_mw110) {
      if (z.re != 0.0) {
        t1 = fabs(atan(z.im / z.re));
      } else {
        t1 = 1.5707963267948966;
      }

      if (28.66666665740641 > caz) {
        ak = 1.8976999933151775 * ak / (2.2204460492503131E-16 * sqrt(sqrt(caz)));
        ak = (log(ak) + caz * cos(3.0 * t1 / (1.0 + caz)) / (1.0 + 0.008 * caz))
          / cos(14.7 * t1 / (28.0 + caz));
        ak = 0.12125 * ak * ak / caz + 1.5;
        guard2 = true;
      } else {
        etest = ak / (3.1415926535897931 * caz * 2.2204460492503131E-16);
        ak = 1.0;
        if (etest >= 1.0) {
          fks = 2.0;
          rk = (caz + caz) + 2.0;
          a1 = 0.0;
          a2 = 1.0;
          earlyExit = true;
          i = 1;
          exitg1 = false;
          while ((!exitg1) && (i < 31)) {
            tm = a2;
            a2 = rk / (ak + 1.0) * a2 - fhs / fks * a1;
            a1 = tm;
            rk += 2.0;
            fks = ((fks + ak) + ak) + 2.0;
            fhs = (fhs + ak) + ak;
            ak++;
            if (etest < fabs(a2) * ak) {
              earlyExit = false;
              exitg1 = true;
            } else {
              i++;
            }
          }

          if (earlyExit) {
            *nz = -2;
          } else {
            ak += 1.909859317102744 * t1 * sqrt(28.66666665740641 / caz);
            fhs = fabs(0.25 - dnu2);
            guard2 = true;
          }
        } else {
          guard2 = true;
        }
      }
    } else {
      guard1 = true;
    }
  }

  if (guard2) {
    b_fix(&ak);
    k = (int)ak;
    ak = k;
    fks = (double)k * (double)k;
    p1.re = 0.0;
    p1.im = 0.0;
    p2_re = 2.2204460492503131E-16;
    p2_im = 0.0;
    cs_re = 2.2204460492503131E-16;
    tm = 0.0;
    for (i = 1; i <= k; i++) {
      a1 = fks - ak;
      rk = 2.0 / (ak + 1.0);
      fmu.re = p2_re;
      fmu.im = p2_im;
      etest = (ak + z.re) * rk;
      rk *= yy;
      a2 = p2_re;
      p2_re = p2_re * etest - p2_im * rk;
      p2_im = a2 * rk + p2_im * etest;
      p2_re -= p1.re;
      p2_im -= p1.im;
      p2_re *= (fks + ak) / (a1 + fhs);
      p2_im *= (fks + ak) / (a1 + fhs);
      p1 = fmu;
      cs_re += p2_re;
      tm += p2_im;
      fks = (a1 - ak) + 1.0;
      ak--;
    }

    fmu.re = 1.0 / rt_hypotd_snf(cs_re, tm);
    tm = -tm;
    etest = cs_re;
    cs_re = cs_re * fmu.re - tm * 0.0;
    tm = etest * 0.0 + tm * fmu.re;
    etest = fmu.re * p2_re - 0.0 * p2_im;
    a2 = fmu.re * p2_im + 0.0 * p2_re;
    rk = coef.re * etest - coef.im * a2;
    etest = coef.re * a2 + coef.im * etest;
    s1_re = rk * cs_re - etest * tm;
    s1_im = rk * tm + etest * cs_re;
    if (!(inu > 0)) {
      zd = z;
      if (iflag == 1) {
        goto_mw270 = true;
      } else {
        goto_mw240 = true;
      }
    } else {
      fmu.re = 1.0 / rt_hypotd_snf(p2_re, p2_im);
      etest = p1.re;
      p1.re = p1.re * fmu.re - p1.im * 0.0;
      p1.im = etest * 0.0 + p1.im * fmu.re;
      p2_im = -p2_im;
      a2 = p2_re;
      p2_re = p2_re * fmu.re - p2_im * 0.0;
      p2_im = a2 * 0.0 + p2_im * fmu.re;
      etest = p1.re * p2_im + p1.im * p2_re;
      ak = (dnu + 0.5) - (p1.re * p2_re - p1.im * p2_im);
      a2 = 0.0 - (p1.re * p2_im + p1.im * p2_re);
      if (z.im == 0.0) {
        if (0.0 - etest == 0.0) {
          rk = ak / z.re;
          etest = 0.0;
        } else if (ak == 0.0) {
          rk = 0.0;
          etest = (0.0 - etest) / z.re;
        } else {
          rk = ak / z.re;
          etest = (0.0 - etest) / z.re;
        }
      } else if (z.re == 0.0) {
        if (ak == 0.0) {
          rk = (0.0 - etest) / z.im;
          etest = 0.0;
        } else if (0.0 - etest == 0.0) {
          rk = 0.0;
          etest = -(ak / z.im);
        } else {
          rk = (0.0 - etest) / z.im;
          etest = -(ak / z.im);
        }
      } else {
        fks = fabs(z.re);
        rk = fabs(z.im);
        if (fks > rk) {
          fks = z.im / z.re;
          etest = z.re + fks * z.im;
          rk = (ak + fks * a2) / etest;
          etest = (a2 - fks * ak) / etest;
        } else if (rk == fks) {
          if (z.re > 0.0) {
            etest = 0.5;
          } else {
            etest = -0.5;
          }

          if (z.im > 0.0) {
            tm = 0.5;
          } else {
            tm = -0.5;
          }

          rk = (ak * etest + a2 * tm) / fks;
          etest = (a2 * etest - ak * tm) / fks;
        } else {
          fks = z.re / z.im;
          etest = z.im + fks * z.re;
          rk = (fks * ak + a2) / etest;
          etest = (fks * a2 - ak) / etest;
        }
      }

      rk++;
      s2_re = s1_re * rk - s1_im * etest;
      s2_im = s1_re * etest + s1_im * rk;
      goto_mw110 = true;
    }

    guard1 = true;
  }

  if (guard1) {
    if (goto_mw240 || goto_mw270) {
    } else if (goto_mw110) {
      ck_re = (dnu + 1.0) * rz_re;
      ck_im = (dnu + 1.0) * rz_im;
      inu--;
      if (inu > 0) {
        if (iflag == 1) {
          zd = z;
          a2 = z.re;
          k = inu;
          j = 1;
          for (i = 0; i < 2; i++) {
            cy[i].re = 0.0;
            cy[i].im = 0.0;
          }

          i = 1;
          exitg1 = false;
          while ((!exitg1) && (i <= inu)) {
            cs_re = s2_re;
            tm = s2_im;
            etest = s2_re;
            s2_re = s2_re * ck_re - s2_im * ck_im;
            s2_im = etest * ck_im + s2_im * ck_re;
            s2_re += s1_re;
            s2_im += s1_im;
            s1_re = cs_re;
            s1_im = tm;
            ck_re += rz_re;
            ck_im += rz_im;
            etest = log(rt_hypotd_snf(s2_re, s2_im));
            b_guard1 = false;
            if (-a2 + etest >= -700.92179369444591) {
              coef.re = s2_re;
              coef.im = s2_im;
              b_log(&coef);
              p2_re = coef.re + -zd.re;
              p2_im = coef.im + -zd.im;
              p1.re = exp(p2_re) / 2.2204460492503131E-16 * cos(p2_im);
              p1.im = exp(p2_re) / 2.2204460492503131E-16 * sin(p2_im);
              if (cuchk(p1, 1.0020841800044864E-289, 2.2204460492503131E-16) ==
                  0) {
                j = 1 - j;
                cy[j] = p1;
                if (k == i - 1) {
                  kflag = 0;
                  inub = i + 1;
                  s2_re = cy[j].re;
                  s2_im = cy[j].im;
                  j = 1 - j;
                  s1_re = cy[j].re;
                  s1_im = cy[j].im;
                  if (i + 1 <= inu) {
                    goto_mw225 = true;
                  } else {
                    s1_re = s2_re;
                    s1_im = s2_im;
                    goto_mw240 = true;
                  }

                  exitg1 = true;
                } else {
                  k = i;
                  i++;
                }
              } else {
                b_guard1 = true;
              }
            } else {
              b_guard1 = true;
            }

            if (b_guard1) {
              if (etest >= 350.46089684722295) {
                a2 -= 700.92179369444591;
                s1_re = cs_re * 3.9222272510438042E-305 - tm * 0.0;
                s1_im = cs_re * 0.0 + tm * 3.9222272510438042E-305;
                etest = s2_re;
                s2_re = s2_re * 3.9222272510438042E-305 - s2_im * 0.0;
                s2_im = etest * 0.0 + s2_im * 3.9222272510438042E-305;
                zd.re = a2;
                zd.im = yy;
              }

              i++;
            }
          }

          if (goto_mw225 || goto_mw240) {
            goto_mw110 = true;
          } else {
            goto_mw110 = false;
          }

          if (!goto_mw110) {
            s1_re = s2_re;
            s1_im = s2_im;
            goto_mw270 = true;
          }
        } else {
          goto_mw225 = true;
        }
      }

      if (goto_mw225 || goto_mw240 || goto_mw270) {
        goto_mw110 = true;
      } else {
        goto_mw110 = false;
      }

      if (!goto_mw110) {
        s1_re = s2_re;
        s1_im = s2_im;
        zd = z;
        if (iflag != 1) {
          goto_mw240 = true;
        }
      }
    } else {
      goto_mw225 = true;
    }

    if (goto_mw225 || goto_mw240) {
      if (goto_mw225) {
        p1.re = dv1[kflag];
        rk = dv2[kflag];
        while (inub <= inu) {
          cs_re = s2_re;
          tm = s2_im;
          etest = s2_re;
          s2_re = s2_re * ck_re - s2_im * ck_im;
          s2_im = etest * ck_im + s2_im * ck_re;
          s2_re += s1_re;
          s2_im += s1_im;
          s1_re = cs_re;
          s1_im = tm;
          ck_re += rz_re;
          ck_im += rz_im;
          if (kflag + 1 < 3) {
            p2_re = s2_re * p1.re - s2_im * 0.0;
            p2_im = s2_re * 0.0 + s2_im * p1.re;
            etest = fabs(p2_re);
            a2 = fabs(p2_im);
            if ((etest > a2) || rtIsNaN(a2)) {
              b_etest = etest;
            } else {
              b_etest = a2;
            }

            if (b_etest > rk) {
              kflag++;
              rk = dv2[kflag];
              s1_re = dv3[kflag] * (cs_re * p1.re - tm * 0.0);
              s1_im = dv3[kflag] * (cs_re * 0.0 + tm * p1.re);
              s2_re = dv3[kflag] * p2_re;
              s2_im = dv3[kflag] * p2_im;
              p1.re = dv1[kflag];
            }
          }

          inub++;
        }

        s1_re = s2_re;
        s1_im = s2_im;
      }

      y->re = dv1[kflag] * s1_re;
      y->im = dv1[kflag] * s1_im;
    } else {
      y->re = s1_re;
      y->im = s1_im;
      *nz = ckscl(zd, y, 700.92179369444591);
      if (1 <= *nz) {
      } else {
        y->re *= 2.2204460492503131E-16;
        y->im *= 2.2204460492503131E-16;
        if (*nz < -2147483646) {
          k = MAX_int32_T;
        } else {
          k = 1 - *nz;
        }

        if (k >= 2) {
          y->re *= 2.2204460492503131E-16;
          y->im *= 2.2204460492503131E-16;
        }
      }
    }
  }
}

/*
 * File trailer for cbknu.c
 *
 * [EOF]
 */
