/*
 * File: f_Simulacao_dinamica_v3.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "f_Simulacao_dinamica_v3.h"
#include "Calc_Gorkov.h"
#include "f_Simulacao_dinamica_v3_emxutil.h"
#include "f_Simulacao_dinamica_v3_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : const struct0_T *transdutor
 *                const struct1_T *bolinha
 *                emxArray_creal_T *x
 * Return Type  : void
 */
void f_Simulacao_dinamica_v3(const struct0_T *transdutor, const struct1_T
  *bolinha, emxArray_creal_T *x)
{
  double h;
  int i0;
  double m;
  double y;
  double b_y;
  double c_y;
  double d_y;
  double e_y;
  double f_y;
  int i;
  creal_T y_data[1];
  creal_T z_data[1];
  creal_T bolinha_data[1];
  creal_T b_bolinha_data[1];
  creal_T dc0;
  creal_T dc1;
  double A_re;
  double A_im;
  h = bolinha->h;
  i0 = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = 5001;
  emxEnsureCapacity_creal_T(x, i0);
  for (i0 = 0; i0 < 5001; i0++) {
    x->data[i0].re = 0.0;
    x->data[i0].im = 0.0;
  }

  x->data[0].re = bolinha->pos_i[0];
  x->data[0].im = 0.0;
  m = bolinha->rho * 4.0 / 3.0 * 3.1415926535897931 * rt_powd_snf(bolinha->R,
    3.0);

  /*  Teste inicial com integração de Euler Explicito */
  y = transdutor->disc / 2.0;
  b_y = transdutor->disc / 2.0;
  c_y = transdutor->disc / 2.0;
  d_y = transdutor->disc / 2.0;
  e_y = transdutor->disc / 2.0;
  f_y = transdutor->disc / 2.0;
  for (i = 0; i < 5000; i++) {
    i0 = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = 1;
    emxEnsureCapacity_creal_T(x, i0);
    x->data[0].re = bolinha->x[i];
    x->data[0].im = 0.0;
    y_data[0].re = bolinha->y[i];
    y_data[0].im = 0.0;
    z_data[0].re = bolinha->z[i];
    z_data[0].im = 0.0;
    bolinha_data[0].re = bolinha->x[i] + y;
    bolinha_data[0].im = 0.0;
    b_bolinha_data[0].re = bolinha->x[i] - b_y;
    b_bolinha_data[0].im = 0.0;
    dc0 = Calc_Gorkov(bolinha_data, y_data, z_data, transdutor->c, transdutor->A,
                      transdutor->rho0, transdutor->omega, transdutor->k,
                      transdutor->a, transdutor->disc, transdutor->posx,
                      transdutor->posy, transdutor->posz, transdutor->normal,
                      transdutor->phase, transdutor->R);
    dc1 = Calc_Gorkov(b_bolinha_data, y_data, z_data, transdutor->c,
                      transdutor->A, transdutor->rho0, transdutor->omega,
                      transdutor->k, transdutor->a, transdutor->disc,
                      transdutor->posx, transdutor->posy, transdutor->posz,
                      transdutor->normal, transdutor->phase, transdutor->R);
    A_re = -(dc0.re - dc1.re);
    A_im = -(dc0.im - dc1.im);
    bolinha_data[0].re = bolinha->y[i] + c_y;
    bolinha_data[0].im = 0.0;
    Calc_Gorkov(x->data, bolinha_data, z_data, transdutor->c, transdutor->A,
                transdutor->rho0, transdutor->omega, transdutor->k,
                transdutor->a, transdutor->disc, transdutor->posx,
                transdutor->posy, transdutor->posz, transdutor->normal,
                transdutor->phase, transdutor->R);
    bolinha_data[0].re = bolinha->y[i] - d_y;
    bolinha_data[0].im = 0.0;
    Calc_Gorkov(x->data, bolinha_data, z_data, transdutor->c, transdutor->A,
                transdutor->rho0, transdutor->omega, transdutor->k,
                transdutor->a, transdutor->disc, transdutor->posx,
                transdutor->posy, transdutor->posz, transdutor->normal,
                transdutor->phase, transdutor->R);
    bolinha_data[0].re = bolinha->z[i] + e_y;
    bolinha_data[0].im = 0.0;
    Calc_Gorkov(x->data, y_data, bolinha_data, transdutor->c, transdutor->A,
                transdutor->rho0, transdutor->omega, transdutor->k,
                transdutor->a, transdutor->disc, transdutor->posx,
                transdutor->posy, transdutor->posz, transdutor->normal,
                transdutor->phase, transdutor->R);
    bolinha_data[0].re = bolinha->z[i] - f_y;
    bolinha_data[0].im = 0.0;
    Calc_Gorkov(x->data, y_data, bolinha_data, transdutor->c, transdutor->A,
                transdutor->rho0, transdutor->omega, transdutor->k,
                transdutor->a, transdutor->disc, transdutor->posx,
                transdutor->posy, transdutor->posz, transdutor->normal,
                transdutor->phase, transdutor->R);
    if (A_im == 0.0) {
      A_re /= transdutor->disc;
      A_im = 0.0;
    } else if (A_re == 0.0) {
      A_re = 0.0;
      A_im /= transdutor->disc;
    } else {
      A_re /= transdutor->disc;
      A_im /= transdutor->disc;
    }

    if (A_im == 0.0) {
      A_re /= m;
      A_im = 0.0;
    } else if (A_re == 0.0) {
      A_re = 0.0;
      A_im /= m;
    } else {
      A_re /= m;
      A_im /= m;
    }

    x->data[i + 1].re = bolinha->x[i] + (bolinha->vx[i] + h * A_re) * h;
    x->data[i + 1].im = h * A_im * h;
  }
}

/*
 * File trailer for f_Simulacao_dinamica_v3.c
 *
 * [EOF]
 */
