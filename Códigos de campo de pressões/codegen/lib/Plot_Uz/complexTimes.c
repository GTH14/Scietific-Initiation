/*
 * File: complexTimes.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 12:47:38
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Plot_Uz.h"
#include "complexTimes.h"

/* Function Definitions */

/*
 * Arguments    : double *re
 *                double *im
 * Return Type  : double
 */
double rescale(double *re, double *im)
{
  double scale;
  double absim;
  scale = fabs(*re);
  absim = fabs(*im);
  if (scale > absim) {
    if (*re < 0.0) {
      *re = -1.0;
    } else {
      *re = 1.0;
    }

    *im /= scale;
  } else if (absim > scale) {
    *re /= absim;
    if (*im < 0.0) {
      *im = -1.0;
    } else {
      *im = 1.0;
    }

    scale = absim;
  } else {
    if (*re < 0.0) {
      *re = -1.0;
    } else {
      *re = 1.0;
    }

    if (*im < 0.0) {
      *im = -1.0;
    } else {
      *im = 1.0;
    }
  }

  return scale;
}

/*
 * File trailer for complexTimes.c
 *
 * [EOF]
 */
