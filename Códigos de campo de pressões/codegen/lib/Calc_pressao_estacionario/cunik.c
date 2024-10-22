/*
 * File: cunik.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao_estacionario.h"
#include "cunik.h"
#include "sqrt1.h"
#include "cmlri.h"
#include "Calc_pressao_estacionario_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T zr
 *                double fnu
 *                int ikflg
 *                int ipmtr
 *                double tol
 *                int *init
 *                creal_T cwrk[16]
 *                creal_T *phi
 *                creal_T *zeta1
 *                creal_T *zeta2
 *                creal_T *summ
 * Return Type  : void
 */
void b_cunik(const creal_T zr, double fnu, int ikflg, int ipmtr, double tol, int
             *init, creal_T cwrk[16], creal_T *phi, creal_T *zeta1, creal_T
             *zeta2, creal_T *summ)
{
  double rfn;
  double ac;
  boolean_T guard1 = false;
  double t_re;
  double s_re;
  double s_im;
  creal_T sr;
  int i;
  double cfn_im;
  double cfn_re;
  double br;
  double crfn_im;
  double crfn_re;
  int l;
  int exitg1;
  int j;
  static const double C[120] = { 1.0, -0.20833333333333334, 0.125,
    0.3342013888888889, -0.40104166666666669, 0.0703125, -1.0258125964506173,
    1.8464626736111112, -0.8912109375, 0.0732421875, 4.6695844234262474,
    -11.207002616222994, 8.78912353515625, -2.3640869140625, 0.112152099609375,
    -28.212072558200244, 84.636217674600729, -91.818241543240021,
    42.534998745388457, -7.3687943594796321, 0.22710800170898438,
    212.57013003921713, -765.25246814118168, 1059.9904525279999,
    -699.57962737613252, 218.19051174421159, -26.491430486951554,
    0.57250142097473145, -1919.4576623184071, 8061.7221817373093,
    -13586.550006434138, 11655.393336864534, -5305.646978613403,
    1200.9029132163525, -108.09091978839466, 1.7277275025844574,
    20204.291330966149, -96980.598388637518, 192547.00123253153,
    -203400.17728041555, 122200.46498301746, -41192.65496889755,
    7109.5143024893641, -493.915304773088, 6.074042001273483,
    -242919.18790055133, 1.3117636146629772E+6, -2.9980159185381066E+6,
    3.7632712976564039E+6, -2.8135632265865342E+6, 1.2683652733216248E+6,
    -331645.17248456361, 45218.768981362729, -2499.8304818112097,
    24.380529699556064, 3.2844698530720379E+6, -1.9706819118432228E+7,
    5.0952602492664643E+7, -7.4105148211532652E+7, 6.6344512274729028E+7,
    -3.7567176660763353E+7, 1.3288767166421818E+7, -2.7856181280864547E+6,
    308186.40461266239, -13886.08975371704, 110.01714026924674,
    -4.932925366450996E+7, 3.2557307418576574E+8, -9.394623596815784E+8,
    1.55359689957058E+9, -1.6210805521083372E+9, 1.1068428168230145E+9,
    -4.9588978427503031E+8, 1.4206290779753309E+8, -2.447406272573873E+7,
    2.2437681779224495E+6, -84005.433603024081, 551.33589612202059,
    8.1478909611831212E+8, -5.8664814920518475E+9, 1.8688207509295826E+10,
    -3.4632043388158775E+10, 4.1280185579753975E+10, -3.3026599749800724E+10,
    1.79542137311556E+10, -6.5632937926192846E+9, 1.5592798648792574E+9,
    -2.2510566188941526E+8, 1.7395107553978164E+7, -549842.32757228869,
    3038.0905109223841, -1.4679261247695616E+10, 1.144982377320258E+11,
    -3.9909617522446649E+11, 8.1921866954857727E+11, -1.0983751560812233E+12,
    1.0081581068653821E+12, -6.4536486924537646E+11, 2.8790064990615057E+11,
    -8.786707217802327E+10, 1.7634730606834969E+10, -2.1671649832237949E+9,
    1.4315787671888897E+8, -3.8718334425726128E+6, 18257.755474293175,
    2.86464035717679E+11, -2.4062979000285039E+12, 9.1093411852398984E+12,
    -2.0516899410934437E+13, 3.056512551993532E+13, -3.166708858478516E+13,
    2.334836404458184E+13, -1.2320491305598287E+13, 4.6127257808491318E+12,
    -1.1965528801961816E+12, 2.0591450323241E+11, -2.1822927757529224E+10,
    1.2470092935127103E+9, -2.9188388122220814E+7, 118838.42625678325 };

  rfn = 1.0 / fnu;
  ac = fnu * 2.2250738585072014E-305;
  guard1 = false;
  if (*init == 0) {
    summ->re = 0.0;
    summ->im = 0.0;
    if ((fabs(zr.re) > ac) || (fabs(zr.im) > ac)) {
      t_re = zr.re * rfn - zr.im * 0.0;
      ac = zr.re * 0.0 + zr.im * rfn;
      s_re = (t_re * t_re - ac * ac) + 1.0;
      s_im = t_re * ac + ac * t_re;
      sr.re = s_re;
      sr.im = s_im;
      b_sqrt(&sr);
      cfn_im = 1.0 + sr.re;
      cfn_re = sr.im;
      br = t_re;
      if (ac == 0.0) {
        if (cfn_re == 0.0) {
          t_re = cfn_im / t_re;
          ac = 0.0;
        } else if (cfn_im == 0.0) {
          t_re = 0.0;
          ac = cfn_re / br;
        } else {
          t_re = cfn_im / t_re;
          ac = cfn_re / br;
        }
      } else if (t_re == 0.0) {
        if (cfn_im == 0.0) {
          t_re = cfn_re / ac;
          ac = 0.0;
        } else if (cfn_re == 0.0) {
          t_re = 0.0;
          ac = -(cfn_im / ac);
        } else {
          t_re = cfn_re / ac;
          ac = -(cfn_im / ac);
        }
      } else {
        crfn_im = fabs(t_re);
        br = fabs(ac);
        if (crfn_im > br) {
          crfn_re = ac / t_re;
          br = t_re + crfn_re * ac;
          t_re = (cfn_im + crfn_re * cfn_re) / br;
          ac = (cfn_re - crfn_re * cfn_im) / br;
        } else if (br == crfn_im) {
          if (t_re > 0.0) {
            crfn_re = 0.5;
          } else {
            crfn_re = -0.5;
          }

          if (ac > 0.0) {
            br = 0.5;
          } else {
            br = -0.5;
          }

          t_re = (cfn_im * crfn_re + cfn_re * br) / crfn_im;
          ac = (cfn_re * crfn_re - cfn_im * br) / crfn_im;
        } else {
          crfn_re = t_re / ac;
          br = ac + crfn_re * t_re;
          t_re = (crfn_re * cfn_im + cfn_re) / br;
          ac = (crfn_re * cfn_re - cfn_im) / br;
        }
      }

      if ((ac == 0.0) && rtIsNaN(t_re)) {
      } else if ((fabs(t_re) > 8.9884656743115785E+307) || (fabs(ac) >
                  8.9884656743115785E+307)) {
        br = t_re;
        t_re = log(rt_hypotd_snf(t_re / 2.0, ac / 2.0)) + 0.69314718055994529;
        ac = rt_atan2d_snf(ac, br);
      } else {
        br = t_re;
        t_re = log(rt_hypotd_snf(t_re, ac));
        ac = rt_atan2d_snf(ac, br);
      }

      zeta1->re = fnu * t_re - 0.0 * ac;
      zeta1->im = fnu * ac + 0.0 * t_re;
      zeta2->re = fnu * sr.re - 0.0 * sr.im;
      zeta2->im = fnu * sr.im + 0.0 * sr.re;
      if (sr.im == 0.0) {
        sr.re = rfn / sr.re;
        sr.im = 0.0;
      } else if (sr.re == 0.0) {
        if (rfn == 0.0) {
          sr.re = 0.0 / sr.im;
          sr.im = 0.0;
        } else {
          sr.re = 0.0;
          sr.im = -(rfn / sr.im);
        }
      } else {
        crfn_im = fabs(sr.re);
        br = fabs(sr.im);
        if (crfn_im > br) {
          crfn_re = sr.im / sr.re;
          br = sr.re + crfn_re * sr.im;
          sr.re = (rfn + crfn_re * 0.0) / br;
          sr.im = (0.0 - crfn_re * rfn) / br;
        } else if (br == crfn_im) {
          if (sr.re > 0.0) {
            crfn_re = 0.5;
          } else {
            crfn_re = -0.5;
          }

          if (sr.im > 0.0) {
            br = 0.5;
          } else {
            br = -0.5;
          }

          sr.re = (rfn * crfn_re + 0.0 * br) / crfn_im;
          sr.im = (0.0 * crfn_re - rfn * br) / crfn_im;
        } else {
          crfn_re = sr.re / sr.im;
          br = sr.im + crfn_re * sr.re;
          sr.re = crfn_re * rfn / br;
          sr.im = (crfn_re * 0.0 - rfn) / br;
        }
      }

      cwrk[15] = sr;
      b_sqrt(&cwrk[15]);
      phi->re = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) *
        cwrk[15].re;
      phi->im = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) *
        cwrk[15].im;
      if (ipmtr != 0) {
      } else {
        if (s_im == 0.0) {
          cfn_re = 1.0 / s_re;
          cfn_im = 0.0;
        } else if (s_re == 0.0) {
          cfn_re = 0.0;
          cfn_im = -(1.0 / s_im);
        } else {
          crfn_im = fabs(s_re);
          br = fabs(s_im);
          if (crfn_im > br) {
            crfn_re = s_im / s_re;
            br = s_re + crfn_re * s_im;
            cfn_re = (1.0 + crfn_re * 0.0) / br;
            cfn_im = (0.0 - crfn_re) / br;
          } else if (br == crfn_im) {
            if (s_re > 0.0) {
              crfn_re = 0.5;
            } else {
              crfn_re = -0.5;
            }

            if (s_im > 0.0) {
              br = 0.5;
            } else {
              br = -0.5;
            }

            cfn_re = (crfn_re + 0.0 * br) / crfn_im;
            cfn_im = (0.0 * crfn_re - br) / crfn_im;
          } else {
            crfn_re = s_re / s_im;
            br = s_im + crfn_re * s_re;
            cfn_re = crfn_re / br;
            cfn_im = (crfn_re * 0.0 - 1.0) / br;
          }
        }

        cwrk[0].re = 1.0;
        cwrk[0].im = 0.0;
        crfn_re = 1.0;
        crfn_im = 0.0;
        ac = 1.0;
        l = 0;
        *init = 15;
        i = 1;
        do {
          exitg1 = 0;
          if (i < 15) {
            s_re = 0.0;
            s_im = 0.0;
            for (j = 0; j <= i; j++) {
              l++;
              br = s_re;
              s_re = s_re * cfn_re - s_im * cfn_im;
              s_im = br * cfn_im + s_im * cfn_re;
              s_re += C[l];
            }

            br = crfn_re;
            crfn_re = crfn_re * sr.re - crfn_im * sr.im;
            crfn_im = br * sr.im + crfn_im * sr.re;
            cwrk[i].re = crfn_re * s_re - crfn_im * s_im;
            cwrk[i].im = crfn_re * s_im + crfn_im * s_re;
            ac *= rfn;
            if ((ac < tol) && (fabs(cwrk[i].re) + fabs(cwrk[i].im) < tol)) {
              *init = i + 1;
              guard1 = true;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            guard1 = true;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      }
    } else {
      zeta1->re = 1402.9773265065639 + fnu;
      zeta1->im = 0.0;
      zeta2->re = fnu;
      zeta2->im = 0.0;
      phi->re = 1.0;
      phi->im = 0.0;
    }
  } else {
    zeta1->re = 0.0;
    zeta1->im = 0.0;
    zeta2->re = 0.0;
    zeta2->im = 0.0;
    guard1 = true;
  }

  if (guard1) {
    if (ikflg == 2) {
      s_re = 0.0;
      s_im = 0.0;
      t_re = 1.0;
      ac = 0.0;
      for (i = 1; i <= *init; i++) {
        s_re += t_re * cwrk[i - 1].re - ac * cwrk[i - 1].im;
        s_im += t_re * cwrk[i - 1].im + ac * cwrk[i - 1].re;
        t_re = -t_re;
        ac = -ac;
      }

      summ->re = s_re;
      summ->im = s_im;
      phi->re = 1.2533141373155003 * cwrk[15].re;
      phi->im = 1.2533141373155003 * cwrk[15].im;
    } else {
      s_re = 0.0;
      s_im = 0.0;
      for (i = 1; i <= *init; i++) {
        s_re += cwrk[i - 1].re;
        s_im += cwrk[i - 1].im;
      }

      summ->re = s_re;
      summ->im = s_im;
      phi->re = 0.3989422804014327 * cwrk[15].re;
      phi->im = 0.3989422804014327 * cwrk[15].im;
    }
  }
}

/*
 * Arguments    : const creal_T zr
 *                double fnu
 *                int ikflg
 *                int ipmtr
 *                double tol
 *                int init
 *                creal_T *phi
 *                creal_T *zeta1
 *                creal_T *zeta2
 * Return Type  : void
 */
void cunik(const creal_T zr, double fnu, int ikflg, int ipmtr, double tol, int
           init, creal_T *phi, creal_T *zeta1, creal_T *zeta2)
{
  int i;
  double rfn;
  creal_T cwrk[16];
  double ac;
  boolean_T guard1 = false;
  double t_re;
  double s_re;
  double s_im;
  creal_T sr;
  double cfn_im;
  double cfn_re;
  double br;
  double crfn_im;
  double crfn_re;
  int l;
  int exitg1;
  int j;
  static const double C[120] = { 1.0, -0.20833333333333334, 0.125,
    0.3342013888888889, -0.40104166666666669, 0.0703125, -1.0258125964506173,
    1.8464626736111112, -0.8912109375, 0.0732421875, 4.6695844234262474,
    -11.207002616222994, 8.78912353515625, -2.3640869140625, 0.112152099609375,
    -28.212072558200244, 84.636217674600729, -91.818241543240021,
    42.534998745388457, -7.3687943594796321, 0.22710800170898438,
    212.57013003921713, -765.25246814118168, 1059.9904525279999,
    -699.57962737613252, 218.19051174421159, -26.491430486951554,
    0.57250142097473145, -1919.4576623184071, 8061.7221817373093,
    -13586.550006434138, 11655.393336864534, -5305.646978613403,
    1200.9029132163525, -108.09091978839466, 1.7277275025844574,
    20204.291330966149, -96980.598388637518, 192547.00123253153,
    -203400.17728041555, 122200.46498301746, -41192.65496889755,
    7109.5143024893641, -493.915304773088, 6.074042001273483,
    -242919.18790055133, 1.3117636146629772E+6, -2.9980159185381066E+6,
    3.7632712976564039E+6, -2.8135632265865342E+6, 1.2683652733216248E+6,
    -331645.17248456361, 45218.768981362729, -2499.8304818112097,
    24.380529699556064, 3.2844698530720379E+6, -1.9706819118432228E+7,
    5.0952602492664643E+7, -7.4105148211532652E+7, 6.6344512274729028E+7,
    -3.7567176660763353E+7, 1.3288767166421818E+7, -2.7856181280864547E+6,
    308186.40461266239, -13886.08975371704, 110.01714026924674,
    -4.932925366450996E+7, 3.2557307418576574E+8, -9.394623596815784E+8,
    1.55359689957058E+9, -1.6210805521083372E+9, 1.1068428168230145E+9,
    -4.9588978427503031E+8, 1.4206290779753309E+8, -2.447406272573873E+7,
    2.2437681779224495E+6, -84005.433603024081, 551.33589612202059,
    8.1478909611831212E+8, -5.8664814920518475E+9, 1.8688207509295826E+10,
    -3.4632043388158775E+10, 4.1280185579753975E+10, -3.3026599749800724E+10,
    1.79542137311556E+10, -6.5632937926192846E+9, 1.5592798648792574E+9,
    -2.2510566188941526E+8, 1.7395107553978164E+7, -549842.32757228869,
    3038.0905109223841, -1.4679261247695616E+10, 1.144982377320258E+11,
    -3.9909617522446649E+11, 8.1921866954857727E+11, -1.0983751560812233E+12,
    1.0081581068653821E+12, -6.4536486924537646E+11, 2.8790064990615057E+11,
    -8.786707217802327E+10, 1.7634730606834969E+10, -2.1671649832237949E+9,
    1.4315787671888897E+8, -3.8718334425726128E+6, 18257.755474293175,
    2.86464035717679E+11, -2.4062979000285039E+12, 9.1093411852398984E+12,
    -2.0516899410934437E+13, 3.056512551993532E+13, -3.166708858478516E+13,
    2.334836404458184E+13, -1.2320491305598287E+13, 4.6127257808491318E+12,
    -1.1965528801961816E+12, 2.0591450323241E+11, -2.1822927757529224E+10,
    1.2470092935127103E+9, -2.9188388122220814E+7, 118838.42625678325 };

  for (i = 0; i < 16; i++) {
    cwrk[i].re = 0.0;
    cwrk[i].im = 0.0;
  }

  rfn = 1.0 / fnu;
  ac = fnu * 2.2250738585072014E-305;
  guard1 = false;
  if (init == 0) {
    if ((fabs(zr.re) > ac) || (fabs(zr.im) > ac)) {
      t_re = zr.re * rfn - zr.im * 0.0;
      ac = zr.re * 0.0 + zr.im * rfn;
      s_re = (t_re * t_re - ac * ac) + 1.0;
      s_im = t_re * ac + ac * t_re;
      sr.re = s_re;
      sr.im = s_im;
      b_sqrt(&sr);
      cfn_im = 1.0 + sr.re;
      cfn_re = sr.im;
      br = t_re;
      if (ac == 0.0) {
        if (cfn_re == 0.0) {
          t_re = cfn_im / t_re;
          ac = 0.0;
        } else if (cfn_im == 0.0) {
          t_re = 0.0;
          ac = cfn_re / br;
        } else {
          t_re = cfn_im / t_re;
          ac = cfn_re / br;
        }
      } else if (t_re == 0.0) {
        if (cfn_im == 0.0) {
          t_re = cfn_re / ac;
          ac = 0.0;
        } else if (cfn_re == 0.0) {
          t_re = 0.0;
          ac = -(cfn_im / ac);
        } else {
          t_re = cfn_re / ac;
          ac = -(cfn_im / ac);
        }
      } else {
        crfn_im = fabs(t_re);
        br = fabs(ac);
        if (crfn_im > br) {
          crfn_re = ac / t_re;
          br = t_re + crfn_re * ac;
          t_re = (cfn_im + crfn_re * cfn_re) / br;
          ac = (cfn_re - crfn_re * cfn_im) / br;
        } else if (br == crfn_im) {
          if (t_re > 0.0) {
            crfn_re = 0.5;
          } else {
            crfn_re = -0.5;
          }

          if (ac > 0.0) {
            br = 0.5;
          } else {
            br = -0.5;
          }

          t_re = (cfn_im * crfn_re + cfn_re * br) / crfn_im;
          ac = (cfn_re * crfn_re - cfn_im * br) / crfn_im;
        } else {
          crfn_re = t_re / ac;
          br = ac + crfn_re * t_re;
          t_re = (crfn_re * cfn_im + cfn_re) / br;
          ac = (crfn_re * cfn_re - cfn_im) / br;
        }
      }

      if ((ac == 0.0) && rtIsNaN(t_re)) {
      } else if ((fabs(t_re) > 8.9884656743115785E+307) || (fabs(ac) >
                  8.9884656743115785E+307)) {
        br = t_re;
        t_re = log(rt_hypotd_snf(t_re / 2.0, ac / 2.0)) + 0.69314718055994529;
        ac = rt_atan2d_snf(ac, br);
      } else {
        br = t_re;
        t_re = log(rt_hypotd_snf(t_re, ac));
        ac = rt_atan2d_snf(ac, br);
      }

      zeta1->re = fnu * t_re - 0.0 * ac;
      zeta1->im = fnu * ac + 0.0 * t_re;
      zeta2->re = fnu * sr.re - 0.0 * sr.im;
      zeta2->im = fnu * sr.im + 0.0 * sr.re;
      if (sr.im == 0.0) {
        sr.re = rfn / sr.re;
        sr.im = 0.0;
      } else if (sr.re == 0.0) {
        if (rfn == 0.0) {
          sr.re = 0.0 / sr.im;
          sr.im = 0.0;
        } else {
          sr.re = 0.0;
          sr.im = -(rfn / sr.im);
        }
      } else {
        crfn_im = fabs(sr.re);
        br = fabs(sr.im);
        if (crfn_im > br) {
          crfn_re = sr.im / sr.re;
          br = sr.re + crfn_re * sr.im;
          sr.re = (rfn + crfn_re * 0.0) / br;
          sr.im = (0.0 - crfn_re * rfn) / br;
        } else if (br == crfn_im) {
          if (sr.re > 0.0) {
            crfn_re = 0.5;
          } else {
            crfn_re = -0.5;
          }

          if (sr.im > 0.0) {
            br = 0.5;
          } else {
            br = -0.5;
          }

          sr.re = (rfn * crfn_re + 0.0 * br) / crfn_im;
          sr.im = (0.0 * crfn_re - rfn * br) / crfn_im;
        } else {
          crfn_re = sr.re / sr.im;
          br = sr.im + crfn_re * sr.re;
          sr.re = crfn_re * rfn / br;
          sr.im = (crfn_re * 0.0 - rfn) / br;
        }
      }

      cwrk[15] = sr;
      b_sqrt(&cwrk[15]);
      phi->re = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) *
        cwrk[15].re;
      phi->im = (0.3989422804014327 + 0.8543718569140677 * (double)(ikflg - 1)) *
        cwrk[15].im;
      if (ipmtr != 0) {
      } else {
        if (s_im == 0.0) {
          cfn_re = 1.0 / s_re;
          cfn_im = 0.0;
        } else if (s_re == 0.0) {
          cfn_re = 0.0;
          cfn_im = -(1.0 / s_im);
        } else {
          crfn_im = fabs(s_re);
          br = fabs(s_im);
          if (crfn_im > br) {
            crfn_re = s_im / s_re;
            br = s_re + crfn_re * s_im;
            cfn_re = (1.0 + crfn_re * 0.0) / br;
            cfn_im = (0.0 - crfn_re) / br;
          } else if (br == crfn_im) {
            if (s_re > 0.0) {
              crfn_re = 0.5;
            } else {
              crfn_re = -0.5;
            }

            if (s_im > 0.0) {
              br = 0.5;
            } else {
              br = -0.5;
            }

            cfn_re = (crfn_re + 0.0 * br) / crfn_im;
            cfn_im = (0.0 * crfn_re - br) / crfn_im;
          } else {
            crfn_re = s_re / s_im;
            br = s_im + crfn_re * s_re;
            cfn_re = crfn_re / br;
            cfn_im = (crfn_re * 0.0 - 1.0) / br;
          }
        }

        cwrk[0].re = 1.0;
        cwrk[0].im = 0.0;
        crfn_re = 1.0;
        crfn_im = 0.0;
        ac = 1.0;
        l = 0;
        i = 1;
        do {
          exitg1 = 0;
          if (i < 15) {
            s_re = 0.0;
            s_im = 0.0;
            for (j = 0; j <= i; j++) {
              l++;
              br = s_re;
              s_re = s_re * cfn_re - s_im * cfn_im;
              s_im = br * cfn_im + s_im * cfn_re;
              s_re += C[l];
            }

            br = crfn_re;
            crfn_re = crfn_re * sr.re - crfn_im * sr.im;
            crfn_im = br * sr.im + crfn_im * sr.re;
            cwrk[i].re = crfn_re * s_re - crfn_im * s_im;
            cwrk[i].im = crfn_re * s_im + crfn_im * s_re;
            ac *= rfn;
            if ((ac < tol) && (fabs(cwrk[i].re) + fabs(cwrk[i].im) < tol)) {
              guard1 = true;
              exitg1 = 1;
            } else {
              i++;
            }
          } else {
            guard1 = true;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      }
    } else {
      zeta1->re = 1402.9773265065639 + fnu;
      zeta1->im = 0.0;
      zeta2->re = fnu;
      zeta2->im = 0.0;
      phi->re = 1.0;
      phi->im = 0.0;
    }
  } else {
    zeta1->re = 0.0;
    zeta1->im = 0.0;
    zeta2->re = 0.0;
    zeta2->im = 0.0;
    guard1 = true;
  }

  if (guard1) {
    if (ikflg == 2) {
      phi->re = 1.2533141373155003 * cwrk[15].re;
      phi->im = 1.2533141373155003 * cwrk[15].im;
    } else {
      phi->re = 0.3989422804014327 * cwrk[15].re;
      phi->im = 0.3989422804014327 * cwrk[15].im;
    }
  }
}

/*
 * File trailer for cunik.c
 *
 * [EOF]
 */
