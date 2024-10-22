/*
 * File: calc_th.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 12:47:38
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Plot_Uz.h"
#include "calc_th.h"
#include "cmlri.h"
#include "complexTimes.h"
#include "sqrt1.h"
#include "Plot_Uz_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : double x0
 *                double b_y0
 *                double z0
 *                double x
 *                double y
 *                double z
 *                const double normal[3]
 *                const creal_T r
 * Return Type  : creal_T
 */
creal_T calc_th(double x0, double b_y0, double z0, double x, double y, double z,
                const double normal[3], const creal_T r)
{
  creal_T theta;
  double ar;
  double brm;
  double t3;
  creal_T v;
  double t4;
  creal_T u;
  double ci;
  boolean_T xneg;
  double sbi;
  double scaleB;
  ar = (normal[0] * (x - x0) + normal[1] * (y - b_y0)) + normal[2] * (z - z0);
  if (r.im == 0.0) {
    theta.re = ar / r.re;
    theta.im = 0.0;
  } else if (r.re == 0.0) {
    if (ar == 0.0) {
      theta.re = 0.0 / r.im;
      theta.im = 0.0;
    } else {
      theta.re = 0.0;
      theta.im = -(ar / r.im);
    }
  } else {
    brm = fabs(r.re);
    t3 = fabs(r.im);
    if (brm > t3) {
      t3 = r.im / r.re;
      t4 = r.re + t3 * r.im;
      theta.re = (ar + t3 * 0.0) / t4;
      theta.im = (0.0 - t3 * ar) / t4;
    } else if (t3 == brm) {
      if (r.re > 0.0) {
        t3 = 0.5;
      } else {
        t3 = -0.5;
      }

      if (r.im > 0.0) {
        t4 = 0.5;
      } else {
        t4 = -0.5;
      }

      theta.re = (ar * t3 + 0.0 * t4) / brm;
      theta.im = (0.0 * t3 - ar * t4) / brm;
    } else {
      t3 = r.re / r.im;
      t4 = r.im + t3 * r.re;
      theta.re = t3 * ar / t4;
      theta.im = (t3 * 0.0 - ar) / t4;
    }
  }

  if ((theta.im == 0.0) && (!(fabs(theta.re) > 1.0))) {
    theta.re = acos(theta.re);
    theta.im = 0.0;
  } else {
    v.re = 1.0 + theta.re;
    v.im = theta.im;
    b_sqrt(&v);
    u.re = 1.0 - theta.re;
    u.im = 0.0 - theta.im;
    b_sqrt(&u);
    if ((-v.im == 0.0) && (u.im == 0.0)) {
      ci = 0.0;
    } else {
      ci = v.re * u.im + -v.im * u.re;
      if ((!((!rtIsInf(ci)) && (!rtIsNaN(ci)))) && (!rtIsNaN(v.re)) && (!rtIsNaN
           (-v.im)) && (!rtIsNaN(u.re)) && (!rtIsNaN(u.im))) {
        t3 = v.re;
        t4 = -v.im;
        brm = rescale(&t3, &t4);
        ar = u.re;
        sbi = u.im;
        scaleB = rescale(&ar, &sbi);
        if ((!rtIsInf(brm)) && (!rtIsNaN(brm)) && ((!rtIsInf(scaleB)) &&
             (!rtIsNaN(scaleB)))) {
          xneg = true;
        } else {
          xneg = false;
        }

        if (rtIsNaN(ci) || (rtIsInf(ci) && xneg)) {
          ci = t3 * sbi + t4 * ar;
          if (ci != 0.0) {
            ci = ci * brm * scaleB;
          } else {
            if ((rtIsInf(brm) && ((u.re == 0.0) || (u.im == 0.0))) || (rtIsInf
                 (scaleB) && ((v.re == 0.0) || (-v.im == 0.0)))) {
              t3 = v.re * u.im;
              t4 = -v.im * u.re;
              if (rtIsNaN(t3)) {
                t3 = 0.0;
              }

              if (rtIsNaN(t4)) {
                t4 = 0.0;
              }

              ci = t3 + t4;
            }
          }
        }
      }
    }

    xneg = (ci < 0.0);
    if (xneg) {
      ci = -ci;
    }

    if (ci >= 2.68435456E+8) {
      ci = log(ci) + 0.69314718055994529;
    } else if (ci > 2.0) {
      ci = log(2.0 * ci + 1.0 / (sqrt(ci * ci + 1.0) + ci));
    } else {
      t3 = ci * ci;
      ci += t3 / (1.0 + sqrt(1.0 + t3));
      t3 = fabs(ci);
      if ((t3 > 4.503599627370496E+15) || (!((!rtIsInf(ci)) && (!rtIsNaN(ci)))))
      {
        ci = log(1.0 + ci);
      } else {
        if (!(t3 < 2.2204460492503131E-16)) {
          ci = log(1.0 + ci) * (ci / ((1.0 + ci) - 1.0));
        }
      }
    }

    if (xneg) {
      ci = -ci;
    }

    theta.re = 2.0 * rt_atan2d_snf(u.re, v.re);
    theta.im = ci;
  }

  return theta;
}

/*
 * File trailer for calc_th.c
 *
 * [EOF]
 */
