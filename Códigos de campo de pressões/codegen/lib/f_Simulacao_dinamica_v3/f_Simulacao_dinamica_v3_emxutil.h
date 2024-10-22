/*
 * File: f_Simulacao_dinamica_v3_emxutil.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

#ifndef F_SIMULACAO_DINAMICA_V3_EMXUTIL_H
#define F_SIMULACAO_DINAMICA_V3_EMXUTIL_H

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
extern void emxEnsureCapacity_creal_T(emxArray_creal_T *emxArray, int oldNumel);
extern void emxFree_creal_T(emxArray_creal_T **pEmxArray);
extern void emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions);

#endif

/*
 * File trailer for f_Simulacao_dinamica_v3_emxutil.h
 *
 * [EOF]
 */
