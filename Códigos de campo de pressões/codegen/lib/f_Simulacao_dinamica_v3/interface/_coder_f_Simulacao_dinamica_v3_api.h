/*
 * File: _coder_f_Simulacao_dinamica_v3_api.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

#ifndef _CODER_F_SIMULACAO_DINAMICA_V3_API_H
#define _CODER_F_SIMULACAO_DINAMICA_V3_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_f_Simulacao_dinamica_v3_api.h"

/* Type Definitions */
#ifndef struct_emxArray_creal_T
#define struct_emxArray_creal_T

struct emxArray_creal_T
{
  creal_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_creal_T*/

#ifndef typedef_emxArray_creal_T
#define typedef_emxArray_creal_T

typedef struct emxArray_creal_T emxArray_creal_T;

#endif                                 /*typedef_emxArray_creal_T*/

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  real_T freq;
  real_T c;
  real_T A;
  real_T rho0;
  real_T lambda;
  real_T omega;
  real_T k;
  real_T a;
  real_T disc;
  real_T posx[64];
  real_T posy[64];
  real_T posz[64];
  real_T normal[3];
  real_T foco[3];
  real_T phase[64];
  real_T R;
} struct0_T;

#endif                                 /*typedef_struct0_T*/

#ifndef typedef_struct1_T
#define typedef_struct1_T

typedef struct {
  real_T rho;
  real_T R;
  real_T pos_i[3];
  real_T v_i[3];
  real_T h;
  real_T t_v[5001];
  real_T x[50001];
  real_T y[50001];
  real_T z[50001];
  real_T vx[50001];
  real_T vy[50001];
  real_T vz[50001];
  real_T ti;
  real_T tf;
} struct1_T;

#endif                                 /*typedef_struct1_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void f_Simulacao_dinamica_v3(struct0_T *transdutor, struct1_T *bolinha,
  emxArray_creal_T *x);
extern void f_Simulacao_dinamica_v3_api(const mxArray * const prhs[2], const
  mxArray *plhs[1]);
extern void f_Simulacao_dinamica_v3_atexit(void);
extern void f_Simulacao_dinamica_v3_initialize(void);
extern void f_Simulacao_dinamica_v3_terminate(void);
extern void f_Simulacao_dinamica_v3_xil_terminate(void);

#endif

/*
 * File trailer for _coder_f_Simulacao_dinamica_v3_api.h
 *
 * [EOF]
 */
