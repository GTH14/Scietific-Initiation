/*
 * File: cuoik.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 08-Aug-2022 21:44:19
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao.h"
#include "cuoik.h"
#include "cuchk.h"
#include "cmlri.h"
#include "cunhj.h"
#include "cunik.h"
#include "Calc_pressao_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int ikflg
 *                int nin
 *                creal_T y[2]
 *                double tol
 *                double elim
 *                double alim
 * Return Type  : int
 */
int b_cuoik(const creal_T z, double fnu, int kode, int ikflg, int nin, creal_T
            y[2], double tol, double elim, double alim)
{
  int nuf;
  int n;
  creal_T zr;
  int iform;
  double gnn;
  double gnu;
  double aarg;
  creal_T arg;
  creal_T an;
  creal_T phi;
  creal_T cz;
  creal_T zeta2;
  boolean_T guard1 = false;
  if (nin < 2) {
    n = nin;
  } else {
    n = 2;
  }

  nuf = 0;
  if (z.re < 0.0) {
    zr.re = -z.re;
    zr.im = -z.im;
  } else {
    zr = z;
  }

  if (fabs(zr.im) > fabs(z.re) * 1.7321) {
    iform = 2;
  } else {
    iform = 1;
  }

  if (ikflg == 1) {
    if (fnu < 1.0) {
      gnu = 1.0;
    } else {
      gnu = fnu;
    }
  } else {
    gnn = (fnu + (double)n) - 1.0;
    gnu = n;
    if (gnn > gnu) {
      gnu = gnn;
    }
  }

  aarg = 0.0;
  if (iform == 2) {
    if (zr.im <= 0.0) {
      an.re = -zr.im;
      an.im = -zr.re;
    } else {
      an.re = zr.im;
      an.im = -zr.re;
    }

    cunhj(an, gnu, tol, &phi, &arg, &cz, &zeta2);
    cz.re = zeta2.re - cz.re;
    cz.im = zeta2.im - cz.im;
    aarg = rt_hypotd_snf(arg.re, arg.im);
  } else {
    arg.re = 0.0;
    arg.im = 0.0;
    cunik(zr, gnu, ikflg, 1, tol, 0, &phi, &cz, &zeta2);
    cz.re = zeta2.re - cz.re;
    cz.im = zeta2.im - cz.im;
  }

  if (kode == 2) {
    cz.re -= zr.re;
    cz.im -= zr.im;
  }

  if (ikflg == 2) {
    cz.re = -cz.re;
    cz.im = -cz.im;
  }

  gnn = rt_hypotd_snf(phi.re, phi.im);
  if (cz.re >= alim) {
    gnn = cz.re + log(gnn);
    if (iform == 2) {
      gnn = (gnn - 0.25 * log(aarg)) - 1.2655121234846454;
    }

    if (gnn > elim) {
      nuf = -1;
    }
  } else {
    guard1 = false;
    if (cz.re >= -elim) {
      if (cz.re > -alim) {
      } else {
        gnn = cz.re + log(gnn);
        if (iform == 2) {
          gnn = (gnn - 0.25 * log(aarg)) - 1.2655121234846454;
        }

        if (gnn > -elim) {
          if ((phi.im == 0.0) && rtIsNaN(phi.re)) {
          } else if ((fabs(phi.re) > 8.9884656743115785E+307) || (fabs(phi.im) >
                      8.9884656743115785E+307)) {
            phi.im = rt_atan2d_snf(phi.im, phi.re);
          } else {
            phi.im = rt_atan2d_snf(phi.im, phi.re);
          }

          cz.im += phi.im;
          if (iform == 2) {
            phi = arg;
            if ((arg.im == 0.0) && rtIsNaN(arg.re)) {
            } else if ((fabs(arg.re) > 8.9884656743115785E+307) || (fabs(arg.im)
                        > 8.9884656743115785E+307)) {
              phi.im = rt_atan2d_snf(arg.im, arg.re);
            } else {
              phi.im = rt_atan2d_snf(arg.im, arg.re);
            }

            phi.im *= 0.25;
            phi.im = cz.im - phi.im;
            cz.im = phi.im;
          }

          gnn = exp(gnn) / tol;
          zr.re = gnn * cos(cz.im);
          zr.im = gnn * sin(cz.im);
          if (cuchk(zr, 2.2250738585072014E-305 / tol, tol) != 1) {
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      for (iform = 1; iform <= n; iform++) {
        y[iform - 1].re = 0.0;
        y[iform - 1].im = 0.0;
      }

      nuf = n;
    }
  }

  return nuf;
}

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int ikflg
 *                int nin
 *                creal_T *y
 *                double tol
 *                double elim
 *                double alim
 * Return Type  : int
 */
int cuoik(const creal_T z, double fnu, int kode, int ikflg, int nin, creal_T *y,
          double tol, double elim, double alim)
{
  int nuf;
  int n;
  creal_T zr;
  int iform;
  double gnn;
  double gnu;
  double aarg;
  creal_T arg;
  creal_T an;
  creal_T phi;
  creal_T cz;
  creal_T zeta2;
  boolean_T guard1 = false;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  nuf = 0;
  if (z.re < 0.0) {
    zr.re = -z.re;
    zr.im = -z.im;
  } else {
    zr = z;
  }

  if (fabs(zr.im) > fabs(z.re) * 1.7321) {
    iform = 2;
  } else {
    iform = 1;
  }

  if (ikflg == 1) {
    if (fnu < 1.0) {
      gnu = 1.0;
    } else {
      gnu = fnu;
    }
  } else {
    gnn = (fnu + (double)n) - 1.0;
    gnu = n;
    if (gnn > gnu) {
      gnu = gnn;
    }
  }

  aarg = 0.0;
  if (iform == 2) {
    if (zr.im <= 0.0) {
      an.re = -zr.im;
      an.im = -zr.re;
    } else {
      an.re = zr.im;
      an.im = -zr.re;
    }

    cunhj(an, gnu, tol, &phi, &arg, &cz, &zeta2);
    cz.re = zeta2.re - cz.re;
    cz.im = zeta2.im - cz.im;
    aarg = rt_hypotd_snf(arg.re, arg.im);
  } else {
    arg.re = 0.0;
    arg.im = 0.0;
    cunik(zr, gnu, ikflg, 1, tol, 0, &phi, &cz, &zeta2);
    cz.re = zeta2.re - cz.re;
    cz.im = zeta2.im - cz.im;
  }

  if (kode == 2) {
    cz.re -= zr.re;
    cz.im -= zr.im;
  }

  if (ikflg == 2) {
    cz.re = -cz.re;
    cz.im = -cz.im;
  }

  gnn = rt_hypotd_snf(phi.re, phi.im);
  if (cz.re >= alim) {
    gnn = cz.re + log(gnn);
    if (iform == 2) {
      gnn = (gnn - 0.25 * log(aarg)) - 1.2655121234846454;
    }

    if (gnn > elim) {
      nuf = -1;
    }
  } else {
    guard1 = false;
    if (cz.re >= -elim) {
      if (cz.re > -alim) {
      } else {
        gnn = cz.re + log(gnn);
        if (iform == 2) {
          gnn = (gnn - 0.25 * log(aarg)) - 1.2655121234846454;
        }

        if (gnn > -elim) {
          if ((phi.im == 0.0) && rtIsNaN(phi.re)) {
          } else if ((fabs(phi.re) > 8.9884656743115785E+307) || (fabs(phi.im) >
                      8.9884656743115785E+307)) {
            phi.im = rt_atan2d_snf(phi.im, phi.re);
          } else {
            phi.im = rt_atan2d_snf(phi.im, phi.re);
          }

          cz.im += phi.im;
          if (iform == 2) {
            phi = arg;
            if ((arg.im == 0.0) && rtIsNaN(arg.re)) {
            } else if ((fabs(arg.re) > 8.9884656743115785E+307) || (fabs(arg.im)
                        > 8.9884656743115785E+307)) {
              phi.im = rt_atan2d_snf(arg.im, arg.re);
            } else {
              phi.im = rt_atan2d_snf(arg.im, arg.re);
            }

            phi.im *= 0.25;
            phi.im = cz.im - phi.im;
            cz.im = phi.im;
          }

          gnn = exp(gnn) / tol;
          zr.re = gnn * cos(cz.im);
          zr.im = gnn * sin(cz.im);
          if (cuchk(zr, 2.2250738585072014E-305 / tol, tol) != 1) {
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      iform = 1;
      while (iform <= n) {
        y->re = 0.0;
        y->im = 0.0;
        iform = 2;
      }

      nuf = n;
    }
  }

  return nuf;
}

/*
 * File trailer for cuoik.c
 *
 * [EOF]
 */
