/*
 * File: cunhj.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

#ifndef CUNHJ_H
#define CUNHJ_H

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
extern void b_cunhj(const creal_T z, double fnu, int ipmtr, double tol, creal_T *
                    phi, creal_T *arg, creal_T *zeta1, creal_T *zeta2, creal_T
                    *asum, creal_T *bsum);
extern void cunhj(const creal_T z, double fnu, double tol, creal_T *phi, creal_T
                  *arg, creal_T *zeta1, creal_T *zeta2);

#endif

/*
 * File trailer for cunhj.h
 *
 * [EOF]
 */
