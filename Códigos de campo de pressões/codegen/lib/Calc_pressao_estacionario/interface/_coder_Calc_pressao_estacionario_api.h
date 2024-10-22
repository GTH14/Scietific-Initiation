/*
 * File: _coder_Calc_pressao_estacionario_api.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

#ifndef _CODER_CALC_PRESSAO_ESTACIONARIO_API_H
#define _CODER_CALC_PRESSAO_ESTACIONARIO_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_Calc_pressao_estacionario_api.h"

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
  real_T posx[128];
  real_T posy[128];
  real_T posz[128];
  real_T normal[3];
  real_T foco[3];
  real_T phase[128];
  real_T R;
  real_T nu;
} struct0_T;

#endif                                 /*typedef_struct0_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern creal_T Calc_pressao_estacionario(real_T x, real_T y, real_T z, struct0_T
  *transdutor);
extern void Calc_pressao_estacionario_api(const mxArray * const prhs[4], const
  mxArray *plhs[1]);
extern void Calc_pressao_estacionario_atexit(void);
extern void Calc_pressao_estacionario_initialize(void);
extern void Calc_pressao_estacionario_terminate(void);
extern void Calc_pressao_estacionario_xil_terminate(void);

#endif

/*
 * File trailer for _coder_Calc_pressao_estacionario_api.h
 *
 * [EOF]
 */
