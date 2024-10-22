/*
 * File: Calc_pressao.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
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
#include "f_Simulacao_dinamica_v3_types.h"

/* Function Declarations */
extern double Calc_pressao(const creal_T x_data[], const creal_T y_data[], const
  creal_T z_data[], double transdutor_A, double transdutor_omega, double
  transdutor_k, double transdutor_a, const double transdutor_posx[64], const
  double transdutor_posy[64], const double transdutor_posz[64], const double
  transdutor_normal[3], const double transdutor_phase[64]);

#endif

/*
 * File trailer for Calc_pressao.h
 *
 * [EOF]
 */
