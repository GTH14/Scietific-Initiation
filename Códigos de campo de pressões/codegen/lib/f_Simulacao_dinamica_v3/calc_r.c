/*
 * File: calc_r.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "f_Simulacao_dinamica_v3.h"
#include "calc_r.h"
#include "sqrt1.h"
#include "power.h"

/* Function Definitions */

/*
 * Arguments    : double x0
 *                double b_y0
 *                double z0
 *                const creal_T x_data[]
 *                const creal_T y_data[]
 *                const creal_T z_data[]
 *                creal_T r_data[]
 *                int r_size[2]
 * Return Type  : void
 */
void calc_r(double x0, double b_y0, double z0, const creal_T x_data[], const
            creal_T y_data[], const creal_T z_data[], creal_T r_data[], int
            r_size[2])
{
  creal_T b_x_data[1];
  creal_T tmp_data[1];
  int tmp_size[2];
  creal_T b_tmp_data[1];
  b_x_data[0].re = x_data[0].re - x0;
  b_x_data[0].im = x_data[0].im;
  power(b_x_data, r_data, r_size);
  b_x_data[0].re = y_data[0].re - b_y0;
  b_x_data[0].im = y_data[0].im;
  power(b_x_data, tmp_data, tmp_size);
  b_x_data[0].re = z_data[0].re - z0;
  b_x_data[0].im = z_data[0].im;
  power(b_x_data, b_tmp_data, tmp_size);
  r_size[0] = 1;
  r_size[1] = 1;
  r_data[0].re = (r_data[0].re + tmp_data[0].re) + b_tmp_data[0].re;
  r_data[0].im = (r_data[0].im + tmp_data[0].im) + b_tmp_data[0].im;
  b_sqrt(&r_data[0]);
}

/*
 * File trailer for calc_r.c
 *
 * [EOF]
 */
