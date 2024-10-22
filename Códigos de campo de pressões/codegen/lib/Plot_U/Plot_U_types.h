/*
 * File: Plot_U_types.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 14:00:11
 */

#ifndef PLOT_U_TYPES_H
#define PLOT_U_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

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
  double R;
  double nu;
  double foco[3];
  double phase[64];
} struct0_T;

#endif                                 /*typedef_struct0_T*/
#endif

/*
 * File trailer for Plot_U_types.h
 *
 * [EOF]
 */
