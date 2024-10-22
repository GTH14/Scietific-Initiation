/*
 * File: power.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "f_Simulacao_dinamica_v3.h"
#include "power.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T a_data[]
 *                creal_T y_data[]
 *                int y_size[2]
 * Return Type  : void
 */
void power(const creal_T a_data[], creal_T y_data[], int y_size[2])
{
  y_size[0] = 1;
  y_size[1] = 1;
  y_data[0].re = a_data[0].re * a_data[0].re - a_data[0].im * a_data[0].im;
  y_data[0].im = a_data[0].re * a_data[0].im + a_data[0].im * a_data[0].re;
}

/*
 * File trailer for power.c
 *
 * [EOF]
 */
