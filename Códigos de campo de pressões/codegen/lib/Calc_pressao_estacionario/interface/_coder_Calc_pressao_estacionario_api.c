/*
 * File: _coder_Calc_pressao_estacionario_api.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_Calc_pressao_estacionario_api.h"
#include "_coder_Calc_pressao_estacionario_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131450U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "Calc_pressao_estacionario",         /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId);
static void c_emlrt_marshallIn(const mxArray *transdutor, const char_T
  *identifier, struct0_T *y);
static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, struct0_T *y);
static void e_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[128]);
static real_T emlrt_marshallIn(const mxArray *x, const char_T *identifier);
static const mxArray *emlrt_marshallOut(const creal_T u);
static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[3]);
static void g_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[3]);
static real_T h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId);
static void i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[128]);
static void j_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[3]);
static void k_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[3]);

/* Function Definitions */

/*
 * Arguments    : const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId)
{
  real_T y;
  y = h_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const mxArray *transdutor
 *                const char_T *identifier
 *                struct0_T *y
 * Return Type  : void
 */
static void c_emlrt_marshallIn(const mxArray *transdutor, const char_T
  *identifier, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(emlrtAlias(transdutor), &thisId, y);
  emlrtDestroyArray(&transdutor);
}

/*
 * Arguments    : const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                struct0_T *y
 * Return Type  : void
 */
static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[17] = { "freq", "c", "A", "rho0", "lambda",
    "omega", "k", "a", "disc", "posx", "posy", "posz", "normal", "foco", "phase",
    "R", "nu" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(emlrtRootTLSGlobal, parentId, u, 17, fieldNames, 0U,
    &dims);
  thisId.fIdentifier = "freq";
  y->freq = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "freq")), &thisId);
  thisId.fIdentifier = "c";
  y->c = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u,
    0, "c")), &thisId);
  thisId.fIdentifier = "A";
  y->A = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u,
    0, "A")), &thisId);
  thisId.fIdentifier = "rho0";
  y->rho0 = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "rho0")), &thisId);
  thisId.fIdentifier = "lambda";
  y->lambda = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a
    (emlrtRootTLSGlobal, u, 0, "lambda")), &thisId);
  thisId.fIdentifier = "omega";
  y->omega = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a
    (emlrtRootTLSGlobal, u, 0, "omega")), &thisId);
  thisId.fIdentifier = "k";
  y->k = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u,
    0, "k")), &thisId);
  thisId.fIdentifier = "a";
  y->a = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u,
    0, "a")), &thisId);
  thisId.fIdentifier = "disc";
  y->disc = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "disc")), &thisId);
  thisId.fIdentifier = "posx";
  e_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "posx")), &thisId, y->posx);
  thisId.fIdentifier = "posy";
  e_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "posy")), &thisId, y->posy);
  thisId.fIdentifier = "posz";
  e_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "posz")), &thisId, y->posz);
  thisId.fIdentifier = "normal";
  f_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "normal")), &thisId, y->normal);
  thisId.fIdentifier = "foco";
  g_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "foco")), &thisId, y->foco);
  thisId.fIdentifier = "phase";
  e_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u, 0,
    "phase")), &thisId, y->phase);
  thisId.fIdentifier = "R";
  y->R = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal, u,
    0, "R")), &thisId);
  thisId.fIdentifier = "nu";
  y->nu = b_emlrt_marshallIn(emlrtAlias(emlrtGetFieldR2013a(emlrtRootTLSGlobal,
    u, 0, "nu")), &thisId);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real_T y[128]
 * Return Type  : void
 */
static void e_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[128])
{
  i_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const mxArray *x
 *                const char_T *identifier
 * Return Type  : real_T
 */
static real_T emlrt_marshallIn(const mxArray *x, const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(emlrtAlias(x), &thisId);
  emlrtDestroyArray(&x);
  return y;
}

/*
 * Arguments    : const creal_T u
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const creal_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxCOMPLEX);
  *mxGetPr(m0) = u.re;
  *mxGetPi(m0) = u.im;
  emlrtFreeImagIfZero(m0);
  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real_T y[3]
 * Return Type  : void
 */
static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[3])
{
  j_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real_T y[3]
 * Return Type  : void
 */
static void g_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T y[3])
{
  k_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 0U,
    &dims);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real_T ret[128]
 * Return Type  : void
 */
static void i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[128])
{
  static const int32_T dims[1] = { 128 };

  int32_T i0;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 1U,
    dims);
  for (i0 = 0; i0 < 128; i0++) {
    ret[i0] = (*(real_T (*)[128])mxGetData(src))[i0];
  }

  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real_T ret[3]
 * Return Type  : void
 */
static void j_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[3])
{
  static const int32_T dims[2] = { 1, 3 };

  int32_T i1;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 2U,
    dims);
  for (i1 = 0; i1 < 3; i1++) {
    ret[i1] = (*(real_T (*)[3])mxGetData(src))[i1];
  }

  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real_T ret[3]
 * Return Type  : void
 */
static void k_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T ret[3])
{
  static const int32_T dims[1] = { 3 };

  int32_T i2;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", false, 1U,
    dims);
  for (i2 = 0; i2 < 3; i2++) {
    ret[i2] = (*(real_T (*)[3])mxGetData(src))[i2];
  }

  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const mxArray * const prhs[4]
 *                const mxArray *plhs[1]
 * Return Type  : void
 */
void Calc_pressao_estacionario_api(const mxArray * const prhs[4], const mxArray *
  plhs[1])
{
  real_T x;
  real_T y;
  real_T z;
  struct0_T transdutor;
  creal_T p_total;

  /* Marshall function inputs */
  x = emlrt_marshallIn(emlrtAliasP(prhs[0]), "x");
  y = emlrt_marshallIn(emlrtAliasP(prhs[1]), "y");
  z = emlrt_marshallIn(emlrtAliasP(prhs[2]), "z");
  c_emlrt_marshallIn(emlrtAliasP(prhs[3]), "transdutor", &transdutor);

  /* Invoke the target function */
  p_total = Calc_pressao_estacionario(x, y, z, &transdutor);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(p_total);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void Calc_pressao_estacionario_atexit(void)
{
  mexFunctionCreateRootTLS();
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  Calc_pressao_estacionario_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void Calc_pressao_estacionario_initialize(void)
{
  mexFunctionCreateRootTLS();
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, false, 0U, 0);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void Calc_pressao_estacionario_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_Calc_pressao_estacionario_api.c
 *
 * [EOF]
 */
