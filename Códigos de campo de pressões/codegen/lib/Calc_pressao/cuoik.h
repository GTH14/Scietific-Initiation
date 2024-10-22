/*
 * File: cuoik.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 08-Aug-2022 21:44:19
 */

#ifndef CUOIK_H
#define CUOIK_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "Calc_pressao_types.h"

/* Function Declarations */
extern int b_cuoik(const creal_T z, double fnu, int kode, int ikflg, int nin,
                   creal_T y[2], double tol, double elim, double alim);
extern int cuoik(const creal_T z, double fnu, int kode, int ikflg, int nin,
                 creal_T *y, double tol, double elim, double alim);

#endif

/*
 * File trailer for cuoik.h
 *
 * [EOF]
 */
