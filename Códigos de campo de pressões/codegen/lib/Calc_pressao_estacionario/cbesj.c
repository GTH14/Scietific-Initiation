/*
 * File: cbesj.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao_estacionario.h"
#include "cbesj.h"
#include "cmlri.h"
#include "cbinu.h"
#include "cuoik.h"
#include "casyi.h"
#include "cseri.h"
#include "Calc_pressao_estacionario_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                creal_T *cy
 *                int *nz
 *                int *ierr
 * Return Type  : void
 */
void cbesj(const creal_T z, double fnu, int kode, creal_T *cy, int *nz, int
           *ierr)
{
  double az;
  int inuh;
  double csgn_re;
  double csgn_im;
  creal_T zn;
  int nn;
  double dfnu;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  int nw;
  creal_T cw[2];
  int b_dfnu;
  double b_az;
  *ierr = 0;
  az = rt_hypotd_snf(z.re, z.im);
  if ((az > 1.0737418235E+9) || (fnu > 1.0737418235E+9)) {
    *ierr = 4;
  } else {
    if ((az > 32767.999992370605) || (fnu > 32767.999992370605)) {
      *ierr = 3;
    }
  }

  inuh = (int)fnu >> 1;
  az = (fnu - (double)(inuh << 1)) * 1.5707963267948966;
  csgn_re = cos(az);
  csgn_im = sin(az);
  if ((inuh & 1) == 1) {
    csgn_re = -csgn_re;
    csgn_im = -csgn_im;
  }

  zn.re = -z.re * 0.0 - (-z.im);
  zn.im = -z.re + -z.im * 0.0;
  if (z.im < 0.0) {
    zn.re = -zn.re;
    zn.im = -zn.im;
    csgn_im = -csgn_im;
  }

  cy->re = 0.0;
  cy->im = 0.0;
  *nz = 0;
  az = rt_hypotd_snf(zn.re, zn.im);
  nn = 1;
  dfnu = fnu;
  guard1 = false;
  guard2 = false;
  guard3 = false;
  if ((az <= 2.0) || (az * az * 0.25 <= fnu + 1.0)) {
    nw = cseri(zn, fnu, kode, 1, cy, 2.2204460492503131E-16, 701.61506577445994,
               665.56491761372422);
    if (nw <= MIN_int32_T) {
      inuh = MAX_int32_T;
    } else {
      inuh = -nw;
    }

    if (nw < 0) {
      *nz = inuh;
    } else {
      *nz = nw;
    }

    nn = 1 - *nz;
    if ((1 - *nz == 0) || (nw >= 0)) {
    } else {
      dfnu = (fnu + (1.0 - (double)*nz)) - 1.0;
      guard3 = true;
    }
  } else {
    guard3 = true;
  }

  if (guard3) {
    if (az < 21.784271729432426) {
      if (dfnu <= 1.0) {
        nw = cmlri(zn, fnu, kode, nn, cy, 2.2204460492503131E-16);
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
        guard2 = true;
      }
    } else if ((dfnu <= 1.0) || (!(2.0 * az < dfnu * dfnu))) {
      nw = casyi(zn, fnu, kode, nn, cy, 21.784271729432426,
                 2.2204460492503131E-16, 701.61506577445994);
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
      guard2 = true;
    }
  }

  if (guard2) {
    nw = cuoik(zn, fnu, kode, 1, nn, cy, 2.2204460492503131E-16,
               701.61506577445994, 665.56491761372422);
    if (nw < 0) {
      if (nw == -2) {
        *nz = -2;
      } else {
        *nz = -1;
      }
    } else {
      *nz += nw;
      nn -= nw;
      if (nn == 0) {
      } else {
        dfnu = (fnu + (double)nn) - 1.0;
        if ((dfnu > 85.921358647162123) || (az > 85.921358647162123)) {
          if (dfnu > 86.921358647162123) {
            b_dfnu = 0;
          } else {
            b_dfnu = (int)((85.921358647162123 - dfnu) + 1.0);
          }

          cbuni(zn, fnu, kode, nn, cy, b_dfnu, 85.921358647162123,
                2.2204460492503131E-16, 701.61506577445994, 665.56491761372422,
                &inuh, &nw);
          if (nw < 0) {
            if (nw == -2) {
              *nz = -2;
            } else {
              *nz = -1;
            }
          } else {
            *nz += nw;
            if (inuh == 0) {
            } else {
              nn = inuh;
              guard1 = true;
            }
          }
        } else {
          guard1 = true;
        }
      }
    }
  }

  if (guard1) {
    if (az > 21.784271729432426) {
      for (inuh = 0; inuh < 2; inuh++) {
        cw[inuh].re = 0.0;
        cw[inuh].im = 0.0;
      }

      nw = b_cuoik(zn, fnu, kode, 2, 2, cw, 2.2204460492503131E-16,
                   701.61506577445994, 665.56491761372422);
      if (nw < 0) {
        *nz = nn;
        inuh = 1;
        while (inuh <= nn) {
          cy->re = 0.0;
          cy->im = 0.0;
          inuh = 2;
        }
      } else if (nw > 0) {
        *nz = -1;
      } else {
        nw = cwrsk(zn, fnu, kode, nn, cy, cw, 2.2204460492503131E-16,
                   701.61506577445994, 665.56491761372422);
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
    } else {
      nw = cmlri(zn, fnu, kode, nn, cy, 2.2204460492503131E-16);
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
      az = fabs(cy->re);
      dfnu = fabs(cy->im);
      if ((az > dfnu) || rtIsNaN(dfnu)) {
        b_az = az;
      } else {
        b_az = dfnu;
      }

      if (b_az <= 1.0020841800044864E-289) {
        zn.re = 4.503599627370496E+15 * cy->re;
        zn.im = 4.503599627370496E+15 * cy->im;
        az = 2.2204460492503131E-16;
      } else {
        az = 1.0;
      }

      cy->re = az * (zn.re * csgn_re - zn.im * csgn_im);
      cy->im = az * (zn.re * csgn_im + zn.im * csgn_re);
    }
  }
}

/*
 * File trailer for cbesj.c
 *
 * [EOF]
 */
