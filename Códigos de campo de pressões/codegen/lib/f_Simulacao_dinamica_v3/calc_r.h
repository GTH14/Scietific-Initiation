/*
 * File: calc_r.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

#ifndef CALC_R_H
#define CALC_R_H

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
extern void calc_r(double x0, double b_y0, double z0, const creal_T x_data[],
                   const creal_T y_data[], const creal_T z_data[], creal_T
                   r_data[], int r_size[2]);

#endif

/*
 * File trailer for calc_r.h
 *
 * [EOF]
 */
