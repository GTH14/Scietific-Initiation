/*
 * File: Calc_pressao_estacionario_types.h
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

#ifndef CALC_PRESSAO_ESTACIONARIO_TYPES_H
#define CALC_PRESSAO_ESTACIONARIO_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
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
  double posx[128];
  double posy[128];
  double posz[128];
  double normal[3];
  double foco[3];
  double phase[128];
  double R;
  double nu;
} struct0_T;

#endif                                 /*typedef_struct0_T*/
#endif

/*
 * File trailer for Calc_pressao_estacionario_types.h
 *
 * [EOF]
 */
