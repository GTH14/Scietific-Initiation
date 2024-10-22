/*
 * File: calc_th.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

#ifndef CALC_TH_H
#define CALC_TH_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "f_Simulacao_dinamica_v3_types.h"

/* Function Declarations */
extern void calc_th(double x0, double b_y0, double z0, const creal_T x_data[],
                    const creal_T y_data[], const creal_T z_data[], const double
                    normal[3], const creal_T r, creal_T theta_data[], int
                    theta_size[2]);

#endif

/*
 * File trailer for calc_th.h
 *
 * [EOF]
 */
