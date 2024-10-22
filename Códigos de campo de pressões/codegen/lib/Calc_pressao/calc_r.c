/*
 * File: calc_r.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 08-Aug-2022 21:44:19
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao.h"
#include "calc_r.h"

/* Function Definitions */

/*
 * Arguments    : double x0
 *                double b_y0
 *                double z0
 *                double x
 *                double y
 *                double z
 * Return Type  : double
 */
double calc_r(double x0, double b_y0, double z0, double x, double y, double z)
{
  double a;
  double b_a;
  double c_a;
  a = x - x0;
  b_a = y - b_y0;
  c_a = z - z0;
  return sqrt((a * a + b_a * b_a) + c_a * c_a);
}

/*
 * File trailer for calc_r.c
 *
 * [EOF]
 */
