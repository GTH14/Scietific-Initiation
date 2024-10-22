/*
 * File: _coder_Plot_Uz_mex.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 12:47:38
 */

/* Include Files */
#include "_coder_Plot_Uz_api.h"
#include "_coder_Plot_Uz_mex.h"

/* Function Declarations */
static void Plot_Uz_mexFunction(int32_T nlhs, int32_T nrhs, const mxArray *prhs
  [1]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                int32_T nrhs
 *                const mxArray *prhs[1]
 * Return Type  : void
 */
static void Plot_Uz_mexFunction(int32_T nlhs, int32_T nrhs, const mxArray *prhs
  [1])
{
  const mxArray *inputs[1];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 1, 4, 7,
                        "Plot_Uz");
  }

  if (nlhs > 0) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 7,
                        "Plot_Uz");
  }

  /* Temporary copy for mex inputs. */
  if (0 <= nrhs - 1) {
    memcpy((void *)&inputs[0], (void *)&prhs[0], (uint32_T)(nrhs * (int32_T)
            sizeof(const mxArray *)));
  }

  /* Call the function. */
  Plot_Uz_api(inputs);

  /* Module termination. */
  Plot_Uz_terminate();
}

/*
 * Arguments    : int32_T nlhs
 *                const mxArray * const plhs[]
 *                int32_T nrhs
 *                const mxArray * const prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  (void)plhs;
  mexAtExit(Plot_Uz_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  Plot_Uz_initialize();

  /* Dispatch the entry-point. */
  Plot_Uz_mexFunction(nlhs, nrhs, prhs);
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_Plot_Uz_mex.c
 *
 * [EOF]
 */
