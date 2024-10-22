/*
 * File: Calc_Gorkov.h
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

#ifndef CALC_GORKOV_H
#define CALC_GORKOV_H

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
extern creal_T Calc_Gorkov(const creal_T x_data[], const creal_T y_data[], const
  creal_T z_data[], double transdutor_c, double transdutor_A, double
  transdutor_rho0, double transdutor_omega, double transdutor_k, double
  transdutor_a, double transdutor_disc, const double transdutor_posx[64], const
  double transdutor_posy[64], const double transdutor_posz[64], const double
  transdutor_normal[3], const double transdutor_phase[64], double transdutor_R);

#endif

/*
 * File trailer for Calc_Gorkov.h
 *
 * [EOF]
 */
