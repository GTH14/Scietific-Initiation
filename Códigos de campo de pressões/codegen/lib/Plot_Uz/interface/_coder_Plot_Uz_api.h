/*
 * File: _coder_Plot_Uz_api.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 12:47:38
 */

#ifndef _CODER_PLOT_UZ_API_H
#define _CODER_PLOT_UZ_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_Plot_Uz_api.h"

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
  real_T R;
  real_T nu;
  real_T foco[3];
  real_T phase[64];
} struct0_T;

#endif                                 /*typedef_struct0_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void Plot_Uz(struct0_T *transdutor);
extern void Plot_Uz_api(const mxArray * const prhs[1]);
extern void Plot_Uz_atexit(void);
extern void Plot_Uz_initialize(void);
extern void Plot_Uz_terminate(void);
extern void Plot_Uz_xil_terminate(void);

#endif

/*
 * File trailer for _coder_Plot_Uz_api.h
 *
 * [EOF]
 */
