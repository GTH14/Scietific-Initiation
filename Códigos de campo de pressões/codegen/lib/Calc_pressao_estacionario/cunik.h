/*
 * File: cunik.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

#ifndef CUNIK_H
#define CUNIK_H

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
extern void b_cunik(const creal_T zr, double fnu, int ikflg, int ipmtr, double
                    tol, int *init, creal_T cwrk[16], creal_T *phi, creal_T
                    *zeta1, creal_T *zeta2, creal_T *summ);
extern void cunik(const creal_T zr, double fnu, int ikflg, int ipmtr, double tol,
                  int init, creal_T *phi, creal_T *zeta1, creal_T *zeta2);

#endif

/*
 * File trailer for cunik.h
 *
 * [EOF]
 */
