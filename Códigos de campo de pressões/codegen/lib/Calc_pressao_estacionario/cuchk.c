/*
 * File: cuchk.c
 *
 * MATLAB Coder version            : 3.3
 * C/C++ source code generated on  : 30-Mar-2022 16:35:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Calc_pressao_estacionario.h"
#include "cuchk.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T y
 *                double ascle
 *                double tol
 * Return Type  : int
 */
int cuchk(const creal_T y, double ascle, double tol)
{
  int nz;
  double yr;
  double yi;
  double smallpart;
  yr = fabs(y.re);
  yi = fabs(y.im);
  if (yr > yi) {
    smallpart = yi;
    yi = yr;
  } else {
    smallpart = yr;
  }

  if ((smallpart <= ascle) && (yi < smallpart / tol)) {
    nz = 1;
  } else {
    nz = 0;
  }

  return nz;
}

/*
 * File trailer for cuchk.c
 *
 * [EOF]
 */
