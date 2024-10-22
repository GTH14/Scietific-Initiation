/*
 * File: Calc_pressao.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 12:47:38
 */

#ifndef CALC_PRESSAO_H
#define CALC_PRESSAO_H

/* Include Files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "Plot_Uz_types.h"

/* Function Declarations */
extern creal_T Calc_pressao(double x, double y, double z, double transdutor_A,
  double transdutor_omega, double transdutor_k, double transdutor_a, const
  double transdutor_posx[64], const double transdutor_posy[64], const double
  transdutor_posz[64], const double transdutor_normal[3], const double
  transdutor_phase[64]);

#endif

/*
 * File trailer for Calc_pressao.h
 *
 * [EOF]
 */
