/*
 * File: cbinu.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 08-Aug-2022 21:44:19
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao.h"
#include "cbinu.h"
#include "cmlri.h"
#include "cbknu.h"
#include "cuchk.h"
#include "cairy.h"
#include "cunhj.h"
#include "cunik.h"
#include "cuoik.h"
#include "Calc_pressao_rtwutil.h"

/* Function Declarations */
static void b_cuni2(const creal_T z, double fnu, int kode, int nin, creal_T y[2],
                    double fnul, double tol, double elim, double alim, int
                    *nlast, int *nz);
static void cuni1(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                  double tol, double elim, double alim, int *nlast, int *nz);
static void cuni2(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                  double tol, double elim, double alim, int *nlast, int *nz);

/* Function Definitions */

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int nin
 *                creal_T y[2]
 *                double fnul
 *                double tol
 *                double elim
 *                double alim
 *                int *nlast
 *                int *nz
 * Return Type  : void
 */
static void b_cuni2(const creal_T z, double fnu, int kode, int nin, creal_T y[2],
                    double fnul, double tol, double elim, double alim, int
                    *nlast, int *nz)
{
  int n;
  double cssr[3];
  double csrr[3];
  double bry1;
  double yy;
  creal_T zn;
  double zb_re;
  double zb_im;
  signed char cid_im;
  double ffnu;
  int inu;
  double c2_re;
  double c2_im;
  double zar_re;
  double zar_im;
  int in;
  static const cint8_T icv1[4] = { { 1,/* re */
      0                                /* im */
    }, { 0,                            /* re */
      1                                /* im */
    }, { -1,                           /* re */
      0                                /* im */
    }, { 0,                            /* re */
      -1                               /* im */
    } };

  creal_T dai;
  creal_T s1;
  creal_T zeta1;
  creal_T zeta2;
  creal_T ai;
  creal_T unusedU7;
  double br;
  double bi;
  double rs1;
  double brm;
  int exitg1;
  boolean_T goto_mw120;
  int i;
  double b_br;
  boolean_T exitg2;
  double b_bi;
  double fn;
  creal_T phi;
  creal_T asum;
  creal_T bsum;
  boolean_T guard1 = false;
  if (nin < 2) {
    n = nin;
  } else {
    n = 2;
  }

  *nz = 0;
  *nlast = 0;
  cssr[0] = 1.0 / tol;
  cssr[1] = 1.0;
  cssr[2] = tol;
  csrr[0] = tol;
  csrr[1] = 1.0;
  csrr[2] = 1.0 / tol;
  bry1 = 2.2250738585072014E-305 / tol;
  yy = z.im;
  zn.re = z.im;
  zn.im = -z.re;
  zb_re = z.re;
  zb_im = z.im;
  cid_im = -1;
  if (fnu < 0.0) {
    ffnu = ceil(fnu);
  } else {
    ffnu = floor(fnu);
  }

  inu = (int)ffnu - 1;
  ffnu = 1.5707963267948966 * (fnu - ffnu);
  c2_re = cos(ffnu);
  c2_im = sin(ffnu);
  zar_re = c2_re;
  zar_im = c2_im;
  in = inu + n;
  in -= (in >> 2) << 2;
  ffnu = c2_re;
  c2_re = c2_re * (double)icv1[in].re - c2_im * (double)icv1[in].im;
  c2_im = ffnu * (double)icv1[in].im + c2_im * (double)icv1[in].re;
  if (z.im <= 0.0) {
    zn.re = -z.im;
    zn.im = -z.re;
    zb_re = z.re;
    zb_im = -z.im;
    cid_im = 1;
    c2_im = -c2_im;
  }

  if (fnu > 1.0) {
    ffnu = fnu;
  } else {
    ffnu = 1.0;
  }

  b_cunhj(zn, ffnu, 1, tol, &dai, &s1, &zeta1, &zeta2, &ai, &unusedU7);
  if (kode == 1) {
    s1.re = zeta2.re - zeta1.re;
  } else {
    br = zb_re + zeta2.re;
    bi = zb_im + zeta2.im;
    if (bi == 0.0) {
      ffnu = fnu / br;
    } else if (br == 0.0) {
      if (fnu == 0.0) {
        ffnu = 0.0 / bi;
      } else {
        ffnu = 0.0;
      }
    } else {
      brm = fabs(br);
      ffnu = fabs(bi);
      if (brm > ffnu) {
        brm = bi / br;
        ffnu = (fnu + brm * 0.0) / (br + brm * bi);
      } else if (ffnu == brm) {
        if (br > 0.0) {
          b_br = 0.5;
        } else {
          b_br = -0.5;
        }

        if (bi > 0.0) {
          b_bi = 0.5;
        } else {
          b_bi = -0.5;
        }

        ffnu = (fnu * b_br + 0.0 * b_bi) / brm;
      } else {
        brm = br / bi;
        ffnu = brm * fnu / (bi + brm * br);
      }
    }

    s1.re = fnu * ffnu - zeta1.re;
  }

  rs1 = s1.re;
  if (fabs(s1.re) > elim) {
    if (s1.re > 0.0) {
      *nz = -1;
    } else {
      *nz = n;
      for (i = 1; i <= n; i++) {
        y[i - 1].re = 0.0;
        y[i - 1].im = 0.0;
      }
    }
  } else {
    do {
      exitg1 = 0;
      goto_mw120 = false;
      in = -1;
      i = 1;
      exitg2 = false;
      while ((!exitg2) && (i <= n)) {
        fn = fnu + (double)(n - i);
        b_cunhj(zn, fn, 0, tol, &phi, &unusedU7, &zeta1, &zeta2, &asum, &bsum);
        if (kode == 1) {
          s1.re = zeta2.re - zeta1.re;
          s1.im = zeta2.im - zeta1.im;
        } else {
          br = zb_re + zeta2.re;
          bi = zb_im + zeta2.im;
          if (bi == 0.0) {
            rs1 = fn / br;
            ffnu = 0.0;
          } else if (br == 0.0) {
            if (fn == 0.0) {
              rs1 = 0.0 / bi;
              ffnu = 0.0;
            } else {
              rs1 = 0.0;
              ffnu = -(fn / bi);
            }
          } else {
            brm = fabs(br);
            ffnu = fabs(bi);
            if (brm > ffnu) {
              brm = bi / br;
              ffnu = br + brm * bi;
              rs1 = (fn + brm * 0.0) / ffnu;
              ffnu = (0.0 - brm * fn) / ffnu;
            } else if (ffnu == brm) {
              if (br > 0.0) {
                br = 0.5;
              } else {
                br = -0.5;
              }

              if (bi > 0.0) {
                ffnu = 0.5;
              } else {
                ffnu = -0.5;
              }

              rs1 = (fn * br + 0.0 * ffnu) / brm;
              ffnu = (0.0 * br - fn * ffnu) / brm;
            } else {
              brm = br / bi;
              ffnu = bi + brm * br;
              rs1 = brm * fn / ffnu;
              ffnu = (brm * 0.0 - fn) / ffnu;
            }
          }

          s1.re = (fn * rs1 - zeta1.re) + fabs(yy) * 0.0;
          s1.im = (fn * ffnu - zeta1.im) + fabs(yy);
        }

        rs1 = s1.re;
        if (fabs(s1.re) > elim) {
          goto_mw120 = true;
          exitg2 = true;
        } else {
          if (i == 1) {
            in = 1;
          }

          guard1 = false;
          if (fabs(s1.re) >= alim) {
            rs1 = ((s1.re + log(rt_hypotd_snf(phi.re, phi.im))) - 0.25 * log
                   (rt_hypotd_snf(unusedU7.re, unusedU7.im))) -
              1.2655121234846454;
            if (fabs(rs1) > elim) {
              goto_mw120 = true;
              exitg2 = true;
            } else {
              if (i == 1) {
                in = 0;
              }

              if ((rs1 >= 0.0) && (i == 1)) {
                in = 2;
              }

              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            ai = cairy(unusedU7, 0, 2);
            dai = cairy(unusedU7, 1, 2);
            ffnu = (ai.re * asum.re - ai.im * asum.im) + (dai.re * bsum.re -
              dai.im * bsum.im);
            brm = (ai.re * asum.im + ai.im * asum.re) + (dai.re * bsum.im +
              dai.im * bsum.re);
            bi = ffnu * phi.re - brm * phi.im;
            brm = ffnu * phi.im + brm * phi.re;
            ffnu = exp(s1.re) * cssr[in] * cos(s1.im);
            br = exp(s1.re) * cssr[in] * sin(s1.im);
            dai.re = bi * ffnu - brm * br;
            dai.im = bi * br + brm * ffnu;
            if ((in + 1 == 1) && (cuchk(dai, bry1, tol) != 0)) {
              goto_mw120 = true;
              exitg2 = true;
            } else {
              if (yy <= 0.0) {
                dai.im = -dai.im;
              }

              y[n - i].re = csrr[in] * (dai.re * c2_re - dai.im * c2_im);
              y[n - i].im = csrr[in] * (dai.re * c2_im + dai.im * c2_re);
              ffnu = c2_re;
              c2_re = c2_re * 0.0 - c2_im * (double)cid_im;
              c2_im = ffnu * (double)cid_im + c2_im * 0.0;
              i++;
            }
          }
        }
      }

      if (!goto_mw120) {
        exitg1 = 1;
      } else if (rs1 > 0.0) {
        *nz = -1;
        exitg1 = 1;
      } else {
        y[n - 1].re = 0.0;
        y[n - 1].im = 0.0;
        (*nz)++;
        if (n - 1 == 0) {
          exitg1 = 1;
        } else {
          in = b_cuoik(z, fnu, kode, 1, 1, y, tol, elim, alim);
          if (in < 0) {
            *nz = -1;
            exitg1 = 1;
          } else {
            n = 1 - in;
            *nz += in;
            if (1 - in == 0) {
              exitg1 = 1;
            } else if ((fnu + (1.0 - (double)in)) - 1.0 < fnul) {
              *nlast = 1 - in;
              exitg1 = 1;
            } else {
              in = (inu - in) + 1;
              i = icv1[in - ((in >> 2) << 2)].re;
              in = icv1[in - ((in >> 2) << 2)].im;
              c2_re = zar_re * (double)i - zar_im * (double)in;
              c2_im = zar_re * (double)in + zar_im * (double)i;
              if (yy <= 0.0) {
                c2_im = -c2_im;
              }
            }
          }
        }
      }
    } while (exitg1 == 0);
  }
}

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int nin
 *                creal_T *y
 *                double tol
 *                double elim
 *                double alim
 *                int *nlast
 *                int *nz
 * Return Type  : void
 */
static void cuni1(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                  double tol, double elim, double alim, int *nlast, int *nz)
{
  int n;
  double cssr[3];
  double csrr[3];
  double bry1;
  double fn;
  int iflag;
  creal_T cwrk[16];
  creal_T s1;
  creal_T zeta1;
  creal_T zeta2;
  creal_T unusedU1;
  double fn_re;
  double bi;
  double rs1;
  double brm;
  boolean_T goto_mw110;
  boolean_T exitg1;
  double b_fn_re;
  creal_T summ;
  double b_bi;
  boolean_T guard1 = false;
  double sgnbr;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  *nz = 0;
  *nlast = 0;
  cssr[0] = 1.0 / tol;
  cssr[1] = 1.0;
  cssr[2] = tol;
  csrr[0] = tol;
  csrr[1] = 1.0;
  csrr[2] = 1.0 / tol;
  bry1 = 2.2250738585072014E-305 / tol;
  if (fnu > 1.0) {
    fn = fnu;
  } else {
    fn = 1.0;
  }

  iflag = 0;
  memset(&cwrk[0], 0, sizeof(creal_T) << 4);
  b_cunik(z, fn, 1, 1, tol, &iflag, cwrk, &s1, &zeta1, &zeta2, &unusedU1);
  if (kode == 1) {
    s1.re = zeta2.re - zeta1.re;
  } else {
    fn_re = z.re + zeta2.re;
    bi = z.im + zeta2.im;
    if (bi == 0.0) {
      fn_re = fn / fn_re;
    } else if (fn_re == 0.0) {
      fn_re = 0.0;
    } else {
      brm = fabs(fn_re);
      rs1 = fabs(bi);
      if (brm > rs1) {
        brm = bi / fn_re;
        fn_re = (fn + brm * 0.0) / (fn_re + brm * bi);
      } else if (rs1 == brm) {
        if (fn_re > 0.0) {
          b_fn_re = 0.5;
        } else {
          b_fn_re = -0.5;
        }

        if (bi > 0.0) {
          b_bi = 0.5;
        } else {
          b_bi = -0.5;
        }

        fn_re = (fn * b_fn_re + 0.0 * b_bi) / brm;
      } else {
        brm = fn_re / bi;
        fn_re = brm * fn / (bi + brm * fn_re);
      }
    }

    s1.re = fn * fn_re - zeta1.re;
  }

  rs1 = s1.re;
  if (fabs(s1.re) > elim) {
    if (s1.re > 0.0) {
      *nz = -1;
    } else {
      *nz = n;
      iflag = 1;
      while (iflag <= n) {
        y->re = 0.0;
        y->im = 0.0;
        iflag = 2;
      }
    }
  } else {
    goto_mw110 = false;
    iflag = 1;
    exitg1 = false;
    while ((!exitg1) && (iflag <= n)) {
      iflag = 0;
      b_cunik(z, fnu, 1, 0, tol, &iflag, cwrk, &unusedU1, &zeta1, &zeta2, &summ);
      if (kode == 1) {
        s1.re = zeta2.re - zeta1.re;
        s1.im = zeta2.im - zeta1.im;
      } else {
        fn_re = z.re + zeta2.re;
        bi = z.im + zeta2.im;
        if (bi == 0.0) {
          fn_re = fnu / fn_re;
          rs1 = 0.0;
        } else if (fn_re == 0.0) {
          if (fnu == 0.0) {
            fn_re = 0.0 / bi;
            rs1 = 0.0;
          } else {
            fn_re = 0.0;
            rs1 = -(fnu / bi);
          }
        } else {
          brm = fabs(fn_re);
          rs1 = fabs(bi);
          if (brm > rs1) {
            brm = bi / fn_re;
            rs1 = fn_re + brm * bi;
            fn_re = (fnu + brm * 0.0) / rs1;
            rs1 = (0.0 - brm * fnu) / rs1;
          } else if (rs1 == brm) {
            if (fn_re > 0.0) {
              sgnbr = 0.5;
            } else {
              sgnbr = -0.5;
            }

            if (bi > 0.0) {
              rs1 = 0.5;
            } else {
              rs1 = -0.5;
            }

            fn_re = (fnu * sgnbr + 0.0 * rs1) / brm;
            rs1 = (0.0 * sgnbr - fnu * rs1) / brm;
          } else {
            brm = fn_re / bi;
            rs1 = bi + brm * fn_re;
            fn_re = brm * fnu / rs1;
            rs1 = (brm * 0.0 - fnu) / rs1;
          }
        }

        s1.re = fnu * fn_re - zeta1.re;
        s1.im = (fnu * rs1 - zeta1.im) + z.im;
      }

      rs1 = s1.re;
      if (fabs(s1.re) > elim) {
        goto_mw110 = true;
        exitg1 = true;
      } else {
        iflag = 1;
        guard1 = false;
        if (fabs(s1.re) >= alim) {
          rs1 = s1.re + log(rt_hypotd_snf(unusedU1.re, unusedU1.im));
          if (fabs(rs1) > elim) {
            goto_mw110 = true;
            exitg1 = true;
          } else {
            iflag = 0;
            if (rs1 >= 0.0) {
              iflag = 2;
            }

            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          fn_re = unusedU1.re * summ.re - unusedU1.im * summ.im;
          brm = unusedU1.re * summ.im + unusedU1.im * summ.re;
          sgnbr = exp(s1.re) * cssr[iflag] * cos(s1.im);
          bi = exp(s1.re) * cssr[iflag] * sin(s1.im);
          s1.re = fn_re * sgnbr - brm * bi;
          s1.im = fn_re * bi + brm * sgnbr;
          if ((iflag + 1 == 1) && (cuchk(s1, bry1, tol) != 0)) {
            goto_mw110 = true;
            exitg1 = true;
          } else {
            y->re = csrr[iflag] * s1.re;
            y->im = csrr[iflag] * s1.im;
            iflag = 2;
          }
        }
      }
    }

    if (goto_mw110) {
      if (rs1 > 0.0) {
        *nz = -1;
      } else {
        y->re = 0.0;
        y->im = 0.0;
        *nz = 1;
      }
    }
  }
}

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int nin
 *                creal_T *y
 *                double tol
 *                double elim
 *                double alim
 *                int *nlast
 *                int *nz
 * Return Type  : void
 */
static void cuni2(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                  double tol, double elim, double alim, int *nlast, int *nz)
{
  int n;
  double cssr[3];
  double csrr[3];
  double bry1;
  double yy;
  creal_T zn;
  double zb_re;
  double zb_im;
  signed char cid_im;
  double ffnu;
  double ang;
  int in;
  double re;
  double c2_re;
  static const cint8_T icv0[4] = { { 1,/* re */
      0                                /* im */
    }, { 0,                            /* re */
      1                                /* im */
    }, { -1,                           /* re */
      0                                /* im */
    }, { 0,                            /* re */
      -1                               /* im */
    } };

  double c2_im;
  creal_T dai;
  creal_T s1;
  creal_T zeta1;
  creal_T zeta2;
  creal_T ai;
  creal_T unusedU7;
  double bi;
  double rs1;
  boolean_T goto_mw120;
  boolean_T exitg1;
  double b_ang;
  creal_T phi;
  creal_T asum;
  creal_T bsum;
  double b_bi;
  boolean_T guard1 = false;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  *nz = 0;
  *nlast = 0;
  cssr[0] = 1.0 / tol;
  cssr[1] = 1.0;
  cssr[2] = tol;
  csrr[0] = tol;
  csrr[1] = 1.0;
  csrr[2] = 1.0 / tol;
  bry1 = 2.2250738585072014E-305 / tol;
  yy = z.im;
  zn.re = z.im;
  zn.im = -z.re;
  zb_re = z.re;
  zb_im = z.im;
  cid_im = -1;
  if (fnu < 0.0) {
    ffnu = ceil(fnu);
  } else {
    ffnu = floor(fnu);
  }

  ang = 1.5707963267948966 * (fnu - ffnu);
  in = ((int)ffnu + n) - 1;
  in -= (in >> 2) << 2;
  re = cos(ang);
  ffnu = sin(ang);
  c2_re = re * (double)icv0[in].re - ffnu * (double)icv0[in].im;
  c2_im = re * (double)icv0[in].im + ffnu * (double)icv0[in].re;
  if (z.im <= 0.0) {
    zn.re = -z.im;
    zn.im = -z.re;
    zb_re = z.re;
    zb_im = -z.im;
    cid_im = 1;
    c2_im = -c2_im;
  }

  if (fnu > 1.0) {
    ffnu = fnu;
  } else {
    ffnu = 1.0;
  }

  b_cunhj(zn, ffnu, 1, tol, &dai, &s1, &zeta1, &zeta2, &ai, &unusedU7);
  if (kode == 1) {
    s1.re = zeta2.re - zeta1.re;
  } else {
    ang = zb_re + zeta2.re;
    bi = zb_im + zeta2.im;
    if (bi == 0.0) {
      ffnu = fnu / ang;
    } else if (ang == 0.0) {
      if (fnu == 0.0) {
        ffnu = 0.0 / bi;
      } else {
        ffnu = 0.0;
      }
    } else {
      re = fabs(ang);
      ffnu = fabs(bi);
      if (re > ffnu) {
        re = bi / ang;
        ffnu = (fnu + re * 0.0) / (ang + re * bi);
      } else if (ffnu == re) {
        if (ang > 0.0) {
          b_ang = 0.5;
        } else {
          b_ang = -0.5;
        }

        if (bi > 0.0) {
          b_bi = 0.5;
        } else {
          b_bi = -0.5;
        }

        ffnu = (fnu * b_ang + 0.0 * b_bi) / re;
      } else {
        re = ang / bi;
        ffnu = re * fnu / (bi + re * ang);
      }
    }

    s1.re = fnu * ffnu - zeta1.re;
  }

  rs1 = s1.re;
  if (fabs(s1.re) > elim) {
    if (s1.re > 0.0) {
      *nz = -1;
    } else {
      *nz = n;
      in = 1;
      while (in <= n) {
        y->re = 0.0;
        y->im = 0.0;
        in = 2;
      }
    }
  } else {
    goto_mw120 = false;
    in = 1;
    exitg1 = false;
    while ((!exitg1) && (in <= n)) {
      b_cunhj(zn, fnu, 0, tol, &phi, &unusedU7, &zeta1, &zeta2, &asum, &bsum);
      if (kode == 1) {
        s1.re = zeta2.re - zeta1.re;
        s1.im = zeta2.im - zeta1.im;
      } else {
        ang = zb_re + zeta2.re;
        bi = zb_im + zeta2.im;
        if (bi == 0.0) {
          rs1 = fnu / ang;
          ffnu = 0.0;
        } else if (ang == 0.0) {
          if (fnu == 0.0) {
            rs1 = 0.0 / bi;
            ffnu = 0.0;
          } else {
            rs1 = 0.0;
            ffnu = -(fnu / bi);
          }
        } else {
          re = fabs(ang);
          ffnu = fabs(bi);
          if (re > ffnu) {
            re = bi / ang;
            ffnu = ang + re * bi;
            rs1 = (fnu + re * 0.0) / ffnu;
            ffnu = (0.0 - re * fnu) / ffnu;
          } else if (ffnu == re) {
            if (ang > 0.0) {
              ang = 0.5;
            } else {
              ang = -0.5;
            }

            if (bi > 0.0) {
              ffnu = 0.5;
            } else {
              ffnu = -0.5;
            }

            rs1 = (fnu * ang + 0.0 * ffnu) / re;
            ffnu = (0.0 * ang - fnu * ffnu) / re;
          } else {
            re = ang / bi;
            ffnu = bi + re * ang;
            rs1 = re * fnu / ffnu;
            ffnu = (re * 0.0 - fnu) / ffnu;
          }
        }

        s1.re = (fnu * rs1 - zeta1.re) + fabs(yy) * 0.0;
        s1.im = (fnu * ffnu - zeta1.im) + fabs(yy);
      }

      rs1 = s1.re;
      if (fabs(s1.re) > elim) {
        goto_mw120 = true;
        exitg1 = true;
      } else {
        in = 1;
        guard1 = false;
        if (fabs(s1.re) >= alim) {
          rs1 = ((s1.re + log(rt_hypotd_snf(phi.re, phi.im))) - 0.25 * log
                 (rt_hypotd_snf(unusedU7.re, unusedU7.im))) - 1.2655121234846454;
          if (fabs(rs1) > elim) {
            goto_mw120 = true;
            exitg1 = true;
          } else {
            in = 0;
            if (rs1 >= 0.0) {
              in = 2;
            }

            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          ai = cairy(unusedU7, 0, 2);
          dai = cairy(unusedU7, 1, 2);
          ffnu = (ai.re * asum.re - ai.im * asum.im) + (dai.re * bsum.re -
            dai.im * bsum.im);
          ang = (ai.re * asum.im + ai.im * asum.re) + (dai.re * bsum.im + dai.im
            * bsum.re);
          bi = ffnu * phi.re - ang * phi.im;
          ang = ffnu * phi.im + ang * phi.re;
          re = exp(s1.re) * cssr[in] * cos(s1.im);
          ffnu = exp(s1.re) * cssr[in] * sin(s1.im);
          dai.re = bi * re - ang * ffnu;
          dai.im = bi * ffnu + ang * re;
          if ((in + 1 == 1) && (cuchk(dai, bry1, tol) != 0)) {
            goto_mw120 = true;
            exitg1 = true;
          } else {
            if (yy <= 0.0) {
              dai.im = -dai.im;
            }

            y->re = csrr[in] * (dai.re * c2_re - dai.im * c2_im);
            y->im = csrr[in] * (dai.re * c2_im + dai.im * c2_re);
            ffnu = c2_re;
            c2_re = c2_re * 0.0 - c2_im * (double)cid_im;
            c2_im = ffnu * (double)cid_im + c2_im * 0.0;
            in = 2;
          }
        }
      }
    }

    if (goto_mw120) {
      if (rs1 > 0.0) {
        *nz = -1;
      } else {
        y->re = 0.0;
        y->im = 0.0;
        *nz = 1;
      }
    }
  }
}

/*
 * Arguments    : const creal_T z
 *                double fnu
 *                int kode
 *                int nin
 *                creal_T *y
 *                int nui
 *                double fnul
 *                double tol
 *                double elim
 *                double alim
 *                int *nlast
 *                int *nz
 * Return Type  : void
 */
void cbuni(const creal_T z, double fnu, int kode, int nin, creal_T *y, int nui,
           double fnul, double tol, double elim, double alim, int *nlast, int
           *nz)
{
  int n;
  int iform;
  double fnui;
  double dfnu;
  int nw;
  double gnu;
  creal_T cy[2];
  int nd;
  double cssr[3];
  double rs1;
  double bry0;
  double csrr[3];
  double bry1;
  int iflag;
  double fn;
  double dscl;
  creal_T cwrk[16];
  creal_T s2;
  creal_T s1;
  creal_T zeta2;
  creal_T rz;
  double dscr;
  double brm;
  int exitg1;
  boolean_T goto_mw110;
  double b_bry0;
  boolean_T exitg2;
  double b_rs1;
  int unusedU3;
  boolean_T guard1 = false;
  double c_rs1;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  *nz = 0;
  iform = 1;
  if (fabs(z.im) > fabs(z.re) * 1.7321) {
    iform = 2;
  }

  if (nui == 0) {
    if (iform == 2) {
      cuni2(z, fnu, kode, n, y, tol, elim, alim, nlast, &nw);
    } else {
      cuni1(z, fnu, kode, n, y, tol, elim, alim, nlast, &nw);
    }

    if (nw < 0) {
      *nz = -1;
      if (nw == -2) {
        *nz = -2;
      }
    } else {
      *nz = nw;
    }
  } else {
    fnui = nui;
    dfnu = (fnu + (double)n) - 1.0;
    gnu = dfnu + (double)nui;
    if (iform == 2) {
      for (iform = 0; iform < 2; iform++) {
        cy[iform].re = 0.0;
        cy[iform].im = 0.0;
      }

      b_cuni2(z, gnu, kode, 2, cy, fnul, tol, elim, alim, nlast, &nw);
    } else {
      for (iform = 0; iform < 2; iform++) {
        cy[iform].re = 0.0;
        cy[iform].im = 0.0;
      }

      nw = 0;
      nd = 2;
      *nlast = 0;
      cssr[0] = 1.0 / tol;
      cssr[1] = 1.0;
      cssr[2] = tol;
      csrr[0] = tol;
      csrr[1] = 1.0;
      csrr[2] = 1.0 / tol;
      bry1 = 2.2250738585072014E-305 / tol;
      if (gnu > 1.0) {
        fn = gnu;
      } else {
        fn = 1.0;
      }

      iform = 0;
      memset(&cwrk[0], 0, sizeof(creal_T) << 4);
      b_cunik(z, fn, 1, 1, tol, &iform, cwrk, &s2, &s1, &zeta2, &rz);
      if (kode == 1) {
        s1.re = zeta2.re - s1.re;
      } else {
        bry0 = z.re + zeta2.re;
        rs1 = z.im + zeta2.im;
        if (rs1 == 0.0) {
          dscr = fn / bry0;
        } else if (bry0 == 0.0) {
          dscr = 0.0;
        } else {
          brm = fabs(bry0);
          dscl = fabs(rs1);
          if (brm > dscl) {
            dscl = rs1 / bry0;
            dscr = (fn + dscl * 0.0) / (bry0 + dscl * rs1);
          } else if (dscl == brm) {
            if (bry0 > 0.0) {
              b_bry0 = 0.5;
            } else {
              b_bry0 = -0.5;
            }

            if (rs1 > 0.0) {
              b_rs1 = 0.5;
            } else {
              b_rs1 = -0.5;
            }

            dscr = (fn * b_bry0 + 0.0 * b_rs1) / brm;
          } else {
            dscl = bry0 / rs1;
            dscr = dscl * fn / (rs1 + dscl * bry0);
          }
        }

        s1.re = fn * dscr - s1.re;
      }

      rs1 = s1.re;
      if (fabs(s1.re) > elim) {
        if (s1.re > 0.0) {
          nw = -1;
        } else {
          nw = 2;
          for (iform = 0; iform < 2; iform++) {
            cy[iform].re = 0.0;
            cy[iform].im = 0.0;
          }
        }
      } else {
        do {
          exitg1 = 0;
          iflag = -1;
          goto_mw110 = false;
          iform = 1;
          exitg2 = false;
          while ((!exitg2) && (iform <= nd)) {
            fn = gnu + (double)(nd - iform);
            unusedU3 = 0;
            b_cunik(z, fn, 1, 0, tol, &unusedU3, cwrk, &s2, &s1, &zeta2, &rz);
            if (kode == 1) {
              s1.re = zeta2.re - s1.re;
              s1.im = zeta2.im - s1.im;
            } else {
              bry0 = z.re + zeta2.re;
              rs1 = z.im + zeta2.im;
              if (rs1 == 0.0) {
                dscr = fn / bry0;
                rs1 = 0.0;
              } else if (bry0 == 0.0) {
                if (fn == 0.0) {
                  dscr = 0.0 / rs1;
                  rs1 = 0.0;
                } else {
                  dscr = 0.0;
                  rs1 = -(fn / rs1);
                }
              } else {
                brm = fabs(bry0);
                dscl = fabs(rs1);
                if (brm > dscl) {
                  dscl = rs1 / bry0;
                  bry0 += dscl * rs1;
                  dscr = (fn + dscl * 0.0) / bry0;
                  rs1 = (0.0 - dscl * fn) / bry0;
                } else if (dscl == brm) {
                  if (bry0 > 0.0) {
                    dscl = 0.5;
                  } else {
                    dscl = -0.5;
                  }

                  if (rs1 > 0.0) {
                    bry0 = 0.5;
                  } else {
                    bry0 = -0.5;
                  }

                  dscr = (fn * dscl + 0.0 * bry0) / brm;
                  rs1 = (0.0 * dscl - fn * bry0) / brm;
                } else {
                  dscl = bry0 / rs1;
                  bry0 = rs1 + dscl * bry0;
                  dscr = dscl * fn / bry0;
                  rs1 = (dscl * 0.0 - fn) / bry0;
                }
              }

              s1.re = fn * dscr - s1.re;
              s1.im = fn * rs1 - s1.im;
              s1.im += z.im;
            }

            rs1 = s1.re;
            if (fabs(s1.re) > elim) {
              goto_mw110 = true;
              exitg2 = true;
            } else {
              if (iform == 1) {
                iflag = 1;
              }

              guard1 = false;
              if (fabs(s1.re) >= alim) {
                rs1 = s1.re + log(rt_hypotd_snf(s2.re, s2.im));
                if (fabs(rs1) > elim) {
                  goto_mw110 = true;
                  exitg2 = true;
                } else {
                  if (iform == 1) {
                    iflag = 0;
                  }

                  if ((rs1 >= 0.0) && (iform == 1)) {
                    iflag = 2;
                  }

                  guard1 = true;
                }
              } else {
                guard1 = true;
              }

              if (guard1) {
                brm = s2.re * rz.re - s2.im * rz.im;
                bry0 = s2.re * rz.im + s2.im * rz.re;
                dscl = exp(s1.re) * cssr[iflag] * cos(s1.im);
                dscr = exp(s1.re) * cssr[iflag] * sin(s1.im);
                s2.re = brm * dscl - bry0 * dscr;
                s2.im = brm * dscr + bry0 * dscl;
                if ((iflag + 1 == 1) && (cuchk(s2, bry1, tol) != 0)) {
                  goto_mw110 = true;
                  exitg2 = true;
                } else {
                  cy[nd - iform].re = csrr[iflag] * s2.re;
                  cy[nd - iform].im = csrr[iflag] * s2.im;
                  iform++;
                }
              }
            }
          }

          if (!goto_mw110) {
            exitg1 = 1;
          } else if (rs1 > 0.0) {
            nw = -1;
            exitg1 = 1;
          } else {
            cy[nd - 1].re = 0.0;
            cy[nd - 1].im = 0.0;
            nw++;
            if (nd - 1 == 0) {
              exitg1 = 1;
            } else {
              iform = b_cuoik(z, gnu, kode, 1, 1, cy, tol, elim, alim);
              if (iform < 0) {
                nw = -1;
                exitg1 = 1;
              } else {
                nd = 1 - iform;
                nw += iform;
                if (1 - iform == 0) {
                  exitg1 = 1;
                } else {
                  if (!((gnu + (1.0 - (double)iform)) - 1.0 >= fnul)) {
                    *nlast = 1 - iform;
                    exitg1 = 1;
                  }
                }
              }
            }
          }
        } while (exitg1 == 0);
      }
    }

    if (nw < 0) {
      *nz = -1;
      if (nw == -2) {
        *nz = -2;
      }
    } else if (nw != 0) {
      *nlast = n;
    } else {
      rs1 = rt_hypotd_snf(cy[0].re, cy[0].im);
      bry0 = 2.2250738585072014E-305 / tol;
      bry1 = 1.0 / bry0;
      cssr[0] = bry0;
      cssr[1] = bry1;
      cssr[2] = bry1;
      iflag = 2;
      fn = 1.0;
      dscl = 1.0;
      if (rs1 > bry0) {
        if (rs1 >= bry1) {
          iflag = 3;
          fn = tol;
          dscl = tol;
        }
      } else {
        iflag = 1;
        bry1 = bry0;
        fn = 1.0 / tol;
        dscl = fn;
      }

      dscr = 1.0 / fn;
      s1.re = dscl * cy[1].re;
      s1.im = dscl * cy[1].im;
      s2.re = dscl * cy[0].re;
      s2.im = dscl * cy[0].im;
      if (z.im == 0.0) {
        rz.re = 2.0 / z.re;
        rz.im = 0.0;
      } else if (z.re == 0.0) {
        rz.re = 0.0;
        rz.im = -(2.0 / z.im);
      } else {
        brm = fabs(z.re);
        dscl = fabs(z.im);
        if (brm > dscl) {
          dscl = z.im / z.re;
          bry0 = z.re + dscl * z.im;
          rz.re = (2.0 + dscl * 0.0) / bry0;
          rz.im = (0.0 - dscl * 2.0) / bry0;
        } else if (dscl == brm) {
          if (z.re > 0.0) {
            dscl = 0.5;
          } else {
            dscl = -0.5;
          }

          if (z.im > 0.0) {
            bry0 = 0.5;
          } else {
            bry0 = -0.5;
          }

          rz.re = (2.0 * dscl + 0.0 * bry0) / brm;
          rz.im = (0.0 * dscl - 2.0 * bry0) / brm;
        } else {
          dscl = z.re / z.im;
          bry0 = z.im + dscl * z.re;
          rz.re = dscl * 2.0 / bry0;
          rz.im = (dscl * 0.0 - 2.0) / bry0;
        }
      }

      for (iform = 1; iform <= nui; iform++) {
        *y = s2;
        brm = s2.re;
        s2.re = s2.re * rz.re - s2.im * rz.im;
        s2.im = brm * rz.im + s2.im * rz.re;
        s2.re *= dfnu + fnui;
        s2.im *= dfnu + fnui;
        s2.re += s1.re;
        s2.im += s1.im;
        s1 = *y;
        fnui--;
        if (!(iflag >= 3)) {
          y->re = dscr * s2.re;
          y->im = dscr * s2.im;
          rs1 = fabs(y->re);
          bry0 = fabs(y->im);
          if ((rs1 > bry0) || rtIsNaN(bry0)) {
            c_rs1 = rs1;
          } else {
            c_rs1 = bry0;
          }

          if (!(c_rs1 <= bry1)) {
            iflag++;
            bry1 = cssr[iflag - 1];
            s1.re *= dscr;
            s1.im *= dscr;
            fn *= tol;
            dscr = 1.0 / fn;
            s1.re *= fn;
            s1.im *= fn;
            s2.re = fn * y->re;
            s2.im = fn * y->im;
          }
        }
      }

      y->re = dscr * s2.re;
      y->im = dscr * s2.im;
    }
  }
}

/*
 * Arguments    : const creal_T zr
 *                double fnu
 *                int kode
 *                int nin
 *                creal_T *y
 *                const creal_T cw[2]
 *                double tol
 *                double elim
 *                double alim
 * Return Type  : int
 */
int cwrsk(const creal_T zr, double fnu, int kode, int nin, creal_T *y, const
          creal_T cw[2], double tol, double elim, double alim)
{
  int nz;
  int n;
  int i;
  int itime;
  creal_T b_cw[2];
  double p1_re;
  double p1_im;
  int id;
  double flam;
  double test;
  double pt_re;
  int k;
  double rz_re;
  double brm;
  double rz_im;
  double p2_re;
  double p2_im;
  double t1_re;
  double t1_im;
  double ap2;
  double test1;
  double ap1;
  double pt_im;
  if (nin < 1) {
    n = nin;
  } else {
    n = 1;
  }

  nz = 0;
  for (i = 0; i < 2; i++) {
    b_cw[i] = cw[i];
  }

  itime = b_cbknu(zr, fnu, kode, 2, b_cw, tol, elim, alim);
  if (itime != 0) {
    nz = -1;
    if (itime == -2) {
      nz = -2;
    }
  } else {
    if (!(n < 1)) {
      n = (int)fnu;
      i = (int)rt_hypotd_snf(zr.re, zr.im);
      if (i > 2147483646) {
        id = MAX_int32_T;
      } else {
        id = i + 1;
      }

      if (n > id) {
        id = 0;
      } else {
        id -= n;
      }

      k = 1;
      if (zr.im == 0.0) {
        rz_re = 2.0 / zr.re;
        rz_im = 0.0;
      } else if (zr.re == 0.0) {
        rz_re = 0.0;
        rz_im = -(2.0 / zr.im);
      } else {
        brm = fabs(zr.re);
        flam = fabs(zr.im);
        if (brm > flam) {
          test = zr.im / zr.re;
          flam = zr.re + test * zr.im;
          rz_re = (2.0 + test * 0.0) / flam;
          rz_im = (0.0 - test * 2.0) / flam;
        } else if (flam == brm) {
          if (zr.re > 0.0) {
            test = 0.5;
          } else {
            test = -0.5;
          }

          if (zr.im > 0.0) {
            flam = 0.5;
          } else {
            flam = -0.5;
          }

          rz_re = (2.0 * test + 0.0 * flam) / brm;
          rz_im = (0.0 * test - 2.0 * flam) / brm;
        } else {
          test = zr.re / zr.im;
          flam = zr.im + test * zr.re;
          rz_re = test * 2.0 / flam;
          rz_im = (test * 0.0 - 2.0) / flam;
        }
      }

      flam = (double)i + 1.0;
      test = n;
      if (flam > test) {
        test = flam;
      }

      t1_re = test * rz_re;
      flam = (double)i + 1.0;
      test = n;
      if (flam > test) {
        test = flam;
      }

      t1_im = test * rz_im;
      p2_re = -t1_re;
      p2_im = -t1_im;
      t1_re += rz_re;
      t1_im += rz_im;
      ap2 = rt_hypotd_snf(p2_re, p2_im);
      test1 = sqrt((ap2 + ap2) / tol);
      test = test1;
      p1_re = 1.0;
      p1_im = 0.0;
      itime = 1;
      while (itime <= 2) {
        k++;
        ap1 = ap2;
        pt_re = p2_re;
        pt_im = p2_im;
        brm = p2_re;
        p2_re = t1_re * p2_re - t1_im * p2_im;
        p2_im = t1_re * p2_im + t1_im * brm;
        p2_re = p1_re - p2_re;
        p2_im = p1_im - p2_im;
        p1_re = pt_re;
        p1_im = pt_im;
        t1_re += rz_re;
        t1_im += rz_im;
        ap2 = rt_hypotd_snf(p2_re, p2_im);
        if (!(ap1 <= test)) {
          if (itime == 1) {
            test = rt_hypotd_snf(t1_re, t1_im) * 0.5;
            flam = test + sqrt(test * test - 1.0);
            brm = ap2 / ap1;
            if ((brm < flam) || rtIsNaN(flam)) {
              flam = brm;
            }

            test = test1 * sqrt(flam / (flam * flam - 1.0));
          }

          itime++;
        }
      }

      itime = (k + id) + 1;
      t1_re = itime;
      p1_re = 1.0 / ap2;
      p1_im = 0.0;
      p2_re = 0.0;
      p2_im = 0.0;
      for (i = 1; i <= itime; i++) {
        pt_re = p1_re;
        pt_im = p1_im;
        flam = ((fnu + 1.0) - 1.0) + t1_re;
        ap2 = rz_re * flam - rz_im * 0.0;
        test = rz_re * 0.0 + rz_im * flam;
        brm = ap2 * p1_im + test * p1_re;
        p1_re = (ap2 * p1_re - test * p1_im) + p2_re;
        p1_im = brm + p2_im;
        p2_re = pt_re;
        p2_im = pt_im;
        t1_re--;
      }

      if ((p1_re == 0.0) && (p1_im == 0.0)) {
        p1_re = tol;
        p1_im = tol;
      }

      if (p1_im == 0.0) {
        if (p2_im == 0.0) {
          y->re = p2_re / p1_re;
          y->im = 0.0;
        } else if (p2_re == 0.0) {
          y->re = 0.0;
          y->im = p2_im / p1_re;
        } else {
          y->re = p2_re / p1_re;
          y->im = p2_im / p1_re;
        }
      } else if (p1_re == 0.0) {
        if (p2_re == 0.0) {
          y->re = p2_im / p1_im;
          y->im = 0.0;
        } else if (p2_im == 0.0) {
          y->re = 0.0;
          y->im = -(p2_re / p1_im);
        } else {
          y->re = p2_im / p1_im;
          y->im = -(p2_re / p1_im);
        }
      } else {
        brm = fabs(p1_re);
        flam = fabs(p1_im);
        if (brm > flam) {
          test = p1_im / p1_re;
          flam = p1_re + test * p1_im;
          y->re = (p2_re + test * p2_im) / flam;
          y->im = (p2_im - test * p2_re) / flam;
        } else if (flam == brm) {
          if (p1_re > 0.0) {
            test = 0.5;
          } else {
            test = -0.5;
          }

          if (p1_im > 0.0) {
            flam = 0.5;
          } else {
            flam = -0.5;
          }

          y->re = (p2_re * test + p2_im * flam) / brm;
          y->im = (p2_im * test - p2_re * flam) / brm;
        } else {
          test = p1_re / p1_im;
          flam = p1_im + test * p1_re;
          y->re = (test * p2_re + p2_im) / flam;
          y->im = (test * p2_im - p2_re) / flam;
        }
      }
    }

    if (kode == 1) {
      p1_re = 1.0;
      p1_im = 0.0;
    } else {
      p1_re = cos(zr.im);
      p1_im = sin(zr.im);
    }

    flam = rt_hypotd_snf(b_cw[1].re, b_cw[1].im);
    test = 2.2250738585072014E-305 / tol;
    if (flam > test) {
      test = 1.0 / test;
      if (flam >= test) {
        pt_re = tol;
      } else {
        pt_re = 1.0;
      }
    } else {
      pt_re = 1.0 / tol;
    }

    flam = b_cw[0].re * pt_re - b_cw[0].im * 0.0;
    test = b_cw[0].re * 0.0 + b_cw[0].im * pt_re;
    brm = (y->re * flam - y->im * test) + (b_cw[1].re * pt_re - b_cw[1].im * 0.0);
    flam = (y->re * test + y->im * flam) + (b_cw[1].re * 0.0 + b_cw[1].im *
      pt_re);
    p2_re = brm * zr.re - flam * zr.im;
    p2_im = brm * zr.im + flam * zr.re;
    flam = 1.0 / rt_hypotd_snf(p2_re, p2_im);
    p2_im = -p2_im;
    p2_re *= flam;
    p2_im *= flam;
    p1_re *= flam;
    p1_im *= flam;
    flam = p1_re;
    p1_re = p1_re * p2_re - p1_im * p2_im;
    p1_im = flam * p2_im + p1_im * p2_re;
    y->re = p1_re * pt_re - p1_im * 0.0;
    y->im = p1_re * 0.0 + p1_im * pt_re;
  }

  return nz;
}

/*
 * File trailer for cbinu.c
 *
 * [EOF]
 */
