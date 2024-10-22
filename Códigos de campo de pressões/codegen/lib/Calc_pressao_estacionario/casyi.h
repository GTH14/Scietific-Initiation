/*
 * File: casyi.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

#ifndef CASYI_H
#define CASYI_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "Calc_pressao_estacionario_types.h"

/* Function Declarations */
extern int casyi(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                 double rl, double tol, double elim);

#endif

/*
 * File trailer for casyi.h
 *
 * [EOF]
 */
