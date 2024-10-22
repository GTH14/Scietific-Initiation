/*
 * File: Calc_pressao.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 12:47:38
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Plot_Uz.h"
#include "Calc_pressao.h"
#include "sum.h"
#include "calc_p.h"
#include "calc_H.h"
#include "calc_v.h"
#include "calc_th.h"
#include "calc_r.h"

/* Function Definitions */

/*
 * Arguments    : double x
 *                double y
 *                double z
 *                double transdutor_A
 *                double transdutor_omega
 *                double transdutor_k
 *                double transdutor_a
 *                const double transdutor_posx[64]
 *                const double transdutor_posy[64]
 *                const double transdutor_posz[64]
 *                const double transdutor_normal[3]
 *                const double transdutor_phase[64]
 * Return Type  : creal_T
 */
creal_T Calc_pressao(double x, double y, double z, double transdutor_A, double
                     transdutor_omega, double transdutor_k, double transdutor_a,
                     const double transdutor_posx[64], const double
                     transdutor_posy[64], const double transdutor_posz[64],
                     const double transdutor_normal[3], const double
                     transdutor_phase[64])
{
  int i;
  creal_T p[64];
  creal_T r[64];
  creal_T theta[64];
  creal_T Pax[64];
  double brm;
  double bim;
  creal_T v[64];
  creal_T H[64];
  double d;

  /*  FUNÇÕES */
  for (i = 0; i < 64; i++) {
    r[i].re = calc_r(transdutor_posx[i], transdutor_posy[i], transdutor_posz[i],
                     x, y, z);
    r[i].im = 0.0;
    theta[i] = calc_th(transdutor_posx[i], transdutor_posy[i], transdutor_posz[i],
                       x, y, z, transdutor_normal, r[i]);
    if (r[i].im == 0.0) {
      Pax[i].re = transdutor_A / r[i].re;
      Pax[i].im = 0.0;
    } else if (r[i].re == 0.0) {
      if (transdutor_A == 0.0) {
        Pax[i].re = 0.0 / r[i].im;
        Pax[i].im = 0.0;
      } else {
        Pax[i].re = 0.0;
        Pax[i].im = -(transdutor_A / r[i].im);
      }
    } else {
      brm = fabs(r[i].re);
      bim = fabs(r[i].im);
      if (brm > bim) {
        bim = r[i].im / r[i].re;
        d = r[i].re + bim * r[i].im;
        Pax[i].re = (transdutor_A + bim * 0.0) / d;
        Pax[i].im = (0.0 - bim * transdutor_A) / d;
      } else if (bim == brm) {
        if (r[i].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (r[i].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        Pax[i].re = (transdutor_A * bim + 0.0 * d) / brm;
        Pax[i].im = (0.0 * bim - transdutor_A * d) / brm;
      } else {
        bim = r[i].re / r[i].im;
        d = r[i].im + bim * r[i].re;
        Pax[i].re = bim * transdutor_A / d;
        Pax[i].im = (bim * 0.0 - transdutor_A) / d;
      }
    }

    v[i] = calc_v(transdutor_k, transdutor_a, theta[i]);
    H[i] = calc_H(v[i]);
    p[i] = calc_p(Pax[i], H[i], transdutor_omega, transdutor_k, r[i],
                  transdutor_phase[i]);
  }

  return sum(p);
}

/*
 * File trailer for Calc_pressao.c
 *
 * [EOF]
 */
