/*
 * File: f_Simulacao_dinamica_v3_emxAPI.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

#ifndef F_SIMULACAO_DINAMICA_V3_EMXAPI_H
#define F_SIMULACAO_DINAMICA_V3_EMXAPI_H

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
extern emxArray_creal_T *emxCreateND_creal_T(int numDimensions, int *size);
extern emxArray_creal_T *emxCreateWrapperND_creal_T(creal_T *data, int
  numDimensions, int *size);
extern emxArray_creal_T *emxCreateWrapper_creal_T(creal_T *data, int rows, int
  cols);
extern emxArray_creal_T *emxCreate_creal_T(int rows, int cols);
extern void emxDestroyArray_creal_T(emxArray_creal_T *emxArray);
extern void emxInitArray_creal_T(emxArray_creal_T **pEmxArray, int numDimensions);

#endif

/*
 * File trailer for f_Simulacao_dinamica_v3_emxAPI.h
 *
 * [EOF]
 */
