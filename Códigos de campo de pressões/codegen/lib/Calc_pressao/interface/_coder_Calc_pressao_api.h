/*
 * File: _coder_Calc_pressao_api.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 08-Aug-2022 21:44:19
 */

#ifndef _CODER_CALC_PRESSAO_API_H
#define _CODER_CALC_PRESSAO_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_Calc_pressao_api.h"

/* Type Definitions */
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
  real_T nu;
} struct0_T;

#endif                                 /*typedef_struct0_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern creal_T Calc_pressao(real_T x, real_T y, real_T z, struct0_T *transdutor);
extern void Calc_pressao_api(const mxArray * const prhs[4], const mxArray *plhs
  [1]);
extern void Calc_pressao_atexit(void);
extern void Calc_pressao_initialize(void);
extern void Calc_pressao_terminate(void);
extern void Calc_pressao_xil_terminate(void);

#endif

/*
 * File trailer for _coder_Calc_pressao_api.h
 *
 * [EOF]
 */
