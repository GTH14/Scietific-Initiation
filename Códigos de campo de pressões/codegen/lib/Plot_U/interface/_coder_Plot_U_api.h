/*
 * File: _coder_Plot_U_api.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 14:00:11
 */

#ifndef _CODER_PLOT_U_API_H
#define _CODER_PLOT_U_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_Plot_U_api.h"

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
extern real_T Plot_U(struct0_T *transdutor);
extern void Plot_U_api(const mxArray * const prhs[1], const mxArray *plhs[1]);
extern void Plot_U_atexit(void);
extern void Plot_U_initialize(void);
extern void Plot_U_terminate(void);
extern void Plot_U_xil_terminate(void);

#endif

/*
 * File trailer for _coder_Plot_U_api.h
 *
 * [EOF]
 */
