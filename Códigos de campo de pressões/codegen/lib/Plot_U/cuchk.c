/*
 * File: cuchk.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 14:00:11
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Plot_U.h"
#include "cuchk.h"

/* Function Definitions */

/*
 * Arguments    : const creal_T y
 *                double ascle
 * Return Type  : int
 */
int cuchk(const creal_T y, double ascle)
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

  if ((smallpart <= ascle) && (yi < smallpart / 2.2204460492503131E-16)) {
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
