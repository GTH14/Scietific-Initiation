/*
 * File: f_Simulacao_dinamica_v3_types.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

#ifndef F_SIMULACAO_DINAMICA_V3_TYPES_H
#define F_SIMULACAO_DINAMICA_V3_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef struct_emxArray_creal_T
#define struct_emxArray_creal_T

struct emxArray_creal_T
{
  creal_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
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
  double freq;
  double c;
  double A;
  double rho0;
  double lambda;
  double omega;
  double k;
  double a;
  double disc;
  double posx[64];
  double posy[64];
  double posz[64];
  double normal[3];
  double foco[3];
  double phase[64];
  double R;
} struct0_T;

#endif                                 /*typedef_struct0_T*/

#ifndef typedef_struct1_T
#define typedef_struct1_T

typedef struct {
  double rho;
  double R;
  double pos_i[3];
  double v_i[3];
  double h;
  double t_v[5001];
  double x[50001];
  double y[50001];
  double z[50001];
  double vx[50001];
  double vy[50001];
  double vz[50001];
  double ti;
  double tf;
} struct1_T;

#endif                                 /*typedef_struct1_T*/
#endif

/*
 * File trailer for f_Simulacao_dinamica_v3_types.h
 *
 * [EOF]
 */
