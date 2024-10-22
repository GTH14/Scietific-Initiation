/*
 * File: Calc_Gorkov.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "f_Simulacao_dinamica_v3.h"
#include "Calc_Gorkov.h"
#include "Calc_pressao.h"
#include "f_Simulacao_dinamica_v3_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T x_data[]
 *                const creal_T y_data[]
 *                const creal_T z_data[]
 *                double transdutor_c
 *                double transdutor_A
 *                double transdutor_rho0
 *                double transdutor_omega
 *                double transdutor_k
 *                double transdutor_a
 *                double transdutor_disc
 *                const double transdutor_posx[64]
 *                const double transdutor_posy[64]
 *                const double transdutor_posz[64]
 *                const double transdutor_normal[3]
 *                const double transdutor_phase[64]
 *                double transdutor_R
 * Return Type  : creal_T
 */
creal_T Calc_Gorkov(const creal_T x_data[], const creal_T y_data[], const
                    creal_T z_data[], double transdutor_c, double transdutor_A,
                    double transdutor_rho0, double transdutor_omega, double
                    transdutor_k, double transdutor_a, double transdutor_disc,
                    const double transdutor_posx[64], const double
                    transdutor_posy[64], const double transdutor_posz[64], const
                    double transdutor_normal[3], const double transdutor_phase
                    [64], double transdutor_R)
{
  creal_T U;
  double re;
  creal_T b_x_data[1];
  creal_T c_x_data[1];
  double A;
  double br;
  double bi;
  double vx_re;
  double vx_im;
  double vy_re;
  double vy_im;
  double vz_re;
  re = Calc_pressao(x_data, y_data, z_data, transdutor_A, transdutor_omega,
                    transdutor_k, transdutor_a, transdutor_posx, transdutor_posy,
                    transdutor_posz, transdutor_normal, transdutor_phase);

  /*      x_v1 = [(x-disc/2) (x+disc/2)]; */
  /*      y_v1 = [y y]; */
  /*      z_v1 = [z z]; */
  /*      x_v2 = [x x]; */
  /*      y_v2 = [(y-disc/2) (y+disc/2)]; */
  /*      z_v2 = [z z]; */
  /*      x_v3 = [x x]; */
  /*      y_v3 = [y y]; */
  /*      z_v3 = [(z-disc/2) (z+disc/2)]; */
  /*      p1 = Calc_pressao(x_v1, y_v1, z_v1, transdutor); */
  /*      p2 = Calc_pressao(x_v2, y_v2, z_v2, transdutor); */
  /*      p3 = Calc_pressao(x_v3, y_v3, z_v3, transdutor); */
  /*      potential_v1 = p1./(j*omega*rho0); */
  /*      potential_v2 = p2./(j*omega*rho0); */
  /*      potential_v3 = p3./(j*omega*rho0); */
  b_x_data[0].re = x_data[0].re + transdutor_disc / 2.0;
  b_x_data[0].im = x_data[0].im;
  c_x_data[0].re = x_data[0].re - transdutor_disc / 2.0;
  c_x_data[0].im = x_data[0].im;
  A = Calc_pressao(b_x_data, y_data, z_data, transdutor_A, transdutor_omega,
                   transdutor_k, transdutor_a, transdutor_posx, transdutor_posy,
                   transdutor_posz, transdutor_normal, transdutor_phase) -
    Calc_pressao(c_x_data, y_data, z_data, transdutor_A, transdutor_omega,
                 transdutor_k, transdutor_a, transdutor_posx, transdutor_posy,
                 transdutor_posz, transdutor_normal, transdutor_phase);
  br = transdutor_rho0 * (transdutor_omega * 0.0) * transdutor_disc;
  bi = transdutor_rho0 * transdutor_omega * transdutor_disc;
  if (bi == 0.0) {
    vx_re = A / br;
    vx_im = 0.0;
  } else if (br == 0.0) {
    if (A == 0.0) {
      vx_re = 0.0 / bi;
      vx_im = 0.0;
    } else {
      vx_re = 0.0;
      vx_im = -(A / bi);
    }
  } else {
    vx_re = rtNaN;
    vx_im = rtNaN;
  }

  b_x_data[0].re = y_data[0].re + transdutor_disc / 2.0;
  b_x_data[0].im = y_data[0].im;
  c_x_data[0].re = y_data[0].re - transdutor_disc / 2.0;
  c_x_data[0].im = y_data[0].im;
  A = Calc_pressao(x_data, b_x_data, z_data, transdutor_A, transdutor_omega,
                   transdutor_k, transdutor_a, transdutor_posx, transdutor_posy,
                   transdutor_posz, transdutor_normal, transdutor_phase) -
    Calc_pressao(x_data, c_x_data, z_data, transdutor_A, transdutor_omega,
                 transdutor_k, transdutor_a, transdutor_posx, transdutor_posy,
                 transdutor_posz, transdutor_normal, transdutor_phase);
  br = transdutor_rho0 * (transdutor_omega * 0.0) * transdutor_disc;
  bi = transdutor_rho0 * transdutor_omega * transdutor_disc;
  if (bi == 0.0) {
    vy_re = A / br;
    vy_im = 0.0;
  } else if (br == 0.0) {
    if (A == 0.0) {
      vy_re = 0.0 / bi;
      vy_im = 0.0;
    } else {
      vy_re = 0.0;
      vy_im = -(A / bi);
    }
  } else {
    vy_re = rtNaN;
    vy_im = rtNaN;
  }

  b_x_data[0].re = z_data[0].re + transdutor_disc / 2.0;
  b_x_data[0].im = z_data[0].im;
  c_x_data[0].re = z_data[0].re - transdutor_disc / 2.0;
  c_x_data[0].im = z_data[0].im;
  A = Calc_pressao(x_data, y_data, b_x_data, transdutor_A, transdutor_omega,
                   transdutor_k, transdutor_a, transdutor_posx, transdutor_posy,
                   transdutor_posz, transdutor_normal, transdutor_phase) -
    Calc_pressao(x_data, y_data, c_x_data, transdutor_A, transdutor_omega,
                 transdutor_k, transdutor_a, transdutor_posx, transdutor_posy,
                 transdutor_posz, transdutor_normal, transdutor_phase);
  br = transdutor_rho0 * (transdutor_omega * 0.0) * transdutor_disc;
  bi = transdutor_rho0 * transdutor_omega * transdutor_disc;
  if (bi == 0.0) {
    vz_re = A / br;
    A = 0.0;
  } else if (br == 0.0) {
    if (A == 0.0) {
      vz_re = 0.0 / bi;
      A = 0.0;
    } else {
      vz_re = 0.0;
      A = -(A / bi);
    }
  } else {
    vz_re = rtNaN;
    A = rtNaN;
  }

  br = rt_hypotd_snf(re, 0.0);
  U.re = 6.2831853071795862 * rt_powd_snf(transdutor_R, 3.0) * (0.5 * (br * br) /
    (3.0 * transdutor_rho0 * (transdutor_c * transdutor_c)) - transdutor_rho0 /
    2.0 * (0.5 * (((vx_re * vx_re - vx_im * -vx_im) + (vy_re * vy_re - vy_im *
    -vy_im)) + (vz_re * vz_re - A * -A))));
  U.im = 6.2831853071795862 * rt_powd_snf(transdutor_R, 3.0) * (0.0 -
    transdutor_rho0 / 2.0 * (0.5 * (((vx_re * -vx_im + vx_im * vx_re) + (vy_re *
    -vy_im + vy_im * vy_re)) + (vz_re * -A + A * vz_re))));
  return U;
}

/*
 * File trailer for Calc_Gorkov.c
 *
 * [EOF]
 */
