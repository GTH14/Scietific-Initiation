/*
 * File: cbinu.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

#ifndef CBINU_H
#define CBINU_H

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
extern void cbuni(const creal_T z, double fnu, int kode, int nin, creal_T *y,
                  int nui, double fnul, double tol, double elim, double alim,
                  int *nlast, int *nz);
extern int cwrsk(const creal_T zr, double fnu, int kode, int nin, creal_T *y,
                 const creal_T cw[2], double tol, double elim, double alim);

#endif

/*
 * File trailer for cbinu.h
 *
 * [EOF]
 */
