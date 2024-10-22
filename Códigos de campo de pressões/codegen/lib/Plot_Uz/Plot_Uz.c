/*
 * File: Plot_Uz.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 12:47:38
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Plot_Uz.h"
#include "Calc_pressao.h"
#include "Plot_Uz_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const struct0_T *transdutor
 * Return Type  : void
 */
void Plot_Uz(const struct0_T *transdutor)
{
  emxArray_real_T *y;
  int ibtile;
  double ndbl;
  double apnd;
  double cdiff;
  emxArray_real_T *b_y;
  int n;
  emxArray_real_T *xs;
  int nm1d2;
  emxArray_real_T *zs;
  emxArray_real_T *a;
  int k;
  int s[2];

  /*  Parâmetros de entrada */
  emxInit_real_T(&y, 2);
  if (rtIsNaN(transdutor->disc)) {
    ibtile = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 1;
    emxEnsureCapacity_real_T(y, ibtile);
    y->data[0] = rtNaN;
  } else if ((transdutor->disc == 0.0) || (transdutor->disc < 0.0)) {
    ibtile = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 0;
    emxEnsureCapacity_real_T(y, ibtile);
  } else if (rtIsInf(transdutor->disc)) {
    ibtile = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 1;
    emxEnsureCapacity_real_T(y, ibtile);
    y->data[0] = -0.05;
  } else {
    ndbl = floor(0.1 / transdutor->disc + 0.5);
    apnd = -0.05 + ndbl * transdutor->disc;
    if (transdutor->disc > 0.0) {
      cdiff = apnd - 0.05;
    } else {
      cdiff = 0.05 - apnd;
    }

    if (fabs(cdiff) < 2.2204460492503132E-17) {
      ndbl++;
      apnd = 0.05;
    } else if (cdiff > 0.0) {
      apnd = -0.05 + (ndbl - 1.0) * transdutor->disc;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }

    ibtile = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n;
    emxEnsureCapacity_real_T(y, ibtile);
    if (n > 0) {
      y->data[0] = -0.05;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nm1d2 = (n - 1) / 2;
        for (k = 1; k < nm1d2; k++) {
          ndbl = (double)k * transdutor->disc;
          y->data[k] = -0.05 + ndbl;
          y->data[(n - k) - 1] = apnd - ndbl;
        }

        if (nm1d2 << 1 == n - 1) {
          y->data[nm1d2] = (-0.05 + apnd) / 2.0;
        } else {
          ndbl = (double)nm1d2 * transdutor->disc;
          y->data[nm1d2] = -0.05 + ndbl;
          y->data[nm1d2 + 1] = apnd - ndbl;
        }
      }
    }
  }

  emxInit_real_T(&b_y, 2);
  if (rtIsNaN(transdutor->disc)) {
    ibtile = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = 1;
    emxEnsureCapacity_real_T(b_y, ibtile);
    b_y->data[0] = rtNaN;
  } else if ((transdutor->disc == 0.0) || (transdutor->disc < 0.0)) {
    ibtile = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = 0;
    emxEnsureCapacity_real_T(b_y, ibtile);
  } else if (rtIsInf(transdutor->disc)) {
    ibtile = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = 1;
    emxEnsureCapacity_real_T(b_y, ibtile);
    b_y->data[0] = 0.0;
  } else if (floor(transdutor->disc) == transdutor->disc) {
    ibtile = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = (int)floor(0.1 / transdutor->disc) + 1;
    emxEnsureCapacity_real_T(b_y, ibtile);
    nm1d2 = (int)floor(0.1 / transdutor->disc);
    for (ibtile = 0; ibtile <= nm1d2; ibtile++) {
      b_y->data[b_y->size[0] * ibtile] = transdutor->disc * (double)ibtile;
    }
  } else {
    ndbl = floor(0.1 / transdutor->disc + 0.5);
    apnd = ndbl * transdutor->disc;
    if (transdutor->disc > 0.0) {
      cdiff = apnd - 0.1;
    } else {
      cdiff = 0.1 - apnd;
    }

    if (fabs(cdiff) < 4.4408920985006264E-17) {
      ndbl++;
      apnd = 0.1;
    } else if (cdiff > 0.0) {
      apnd = (ndbl - 1.0) * transdutor->disc;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }

    ibtile = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = n;
    emxEnsureCapacity_real_T(b_y, ibtile);
    if (n > 0) {
      b_y->data[0] = 0.0;
      if (n > 1) {
        b_y->data[n - 1] = apnd;
        nm1d2 = (n - 1) / 2;
        for (k = 1; k < nm1d2; k++) {
          ndbl = (double)k * transdutor->disc;
          b_y->data[k] = ndbl;
          b_y->data[(n - k) - 1] = apnd - ndbl;
        }

        if (nm1d2 << 1 == n - 1) {
          b_y->data[nm1d2] = apnd / 2.0;
        } else {
          ndbl = (double)nm1d2 * transdutor->disc;
          b_y->data[nm1d2] = ndbl;
          b_y->data[nm1d2 + 1] = apnd - ndbl;
        }
      }
    }
  }

  emxInit_real_T(&xs, 2);
  emxInit_real_T(&zs, 2);
  emxInit_real_T(&a, 2);
  if ((y->size[1] == 0) || (b_y->size[1] == 0)) {
    ibtile = xs->size[0] * xs->size[1];
    xs->size[0] = 0;
    xs->size[1] = 0;
    emxEnsureCapacity_real_T(xs, ibtile);
    ibtile = zs->size[0] * zs->size[1];
    zs->size[0] = 0;
    zs->size[1] = 0;
    emxEnsureCapacity_real_T(zs, ibtile);
  } else {
    nm1d2 = y->size[1];
    ibtile = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = nm1d2;
    emxEnsureCapacity_real_T(a, ibtile);
    for (ibtile = 0; ibtile < nm1d2; ibtile++) {
      a->data[a->size[0] * ibtile] = y->data[ibtile];
    }

    s[0] = b_y->size[1];
    nm1d2 = a->size[1];
    ibtile = xs->size[0] * xs->size[1];
    xs->size[0] = s[0];
    xs->size[1] = nm1d2;
    emxEnsureCapacity_real_T(xs, ibtile);
    nm1d2 = s[0];
    for (n = 0; n + 1 <= a->size[1]; n++) {
      ibtile = n * nm1d2;
      for (k = 1; k <= nm1d2; k++) {
        xs->data[(ibtile + k) - 1] = a->data[n];
      }
    }

    s[1] = y->size[1];
    nm1d2 = b_y->size[1];
    ibtile = zs->size[0] * zs->size[1];
    zs->size[0] = nm1d2;
    zs->size[1] = s[1];
    emxEnsureCapacity_real_T(zs, ibtile);
    nm1d2 = b_y->size[1];
    for (n = 1; n <= s[1]; n++) {
      ibtile = (n - 1) * nm1d2;
      for (k = 1; k <= nm1d2; k++) {
        zs->data[(ibtile + k) - 1] = b_y->data[k - 1];
      }
    }
  }

  emxFree_real_T(&a);
  emxFree_real_T(&b_y);
  emxFree_real_T(&y);

  /*  zs = 0:1e-4:.3; */
  /*  n = length(zs); */
  for (ibtile = 0; ibtile < 2; ibtile++) {
    s[ibtile] = xs->size[ibtile];
  }

  for (n = 0; n < s[0]; n++) {
    for (nm1d2 = 0; nm1d2 < s[1]; nm1d2++) {
      Calc_pressao(xs->data[n + xs->size[0] * nm1d2], 0.0, zs->data[n + zs->
                   size[0] * nm1d2], transdutor->A, transdutor->omega,
                   transdutor->k, transdutor->a, transdutor->posx,
                   transdutor->posy, transdutor->posz, transdutor->normal,
                   transdutor->phase);

      /*      x_v1 = [(x-disc/2) (x+disc/2)]; */
      /*      y_v1 = [y y]; */
      /*      z_v1 = [z z]; */
      /*      x_v2 = [x x]; */
      /*      y_v2 = [(y-disc/2) (y+disc/2)]; */
      /*      z_v2 = [z z]; */
      /*      x_v3 = [x x]; */
      /*      y_v3 = [y y]; */
      /*      z_v3 = [(z-disc/2) (z+disc/2)]; */
      /*      p1 = Calc_pressao(x_v1, y_v1, z_v1, transdutor); */
      /*      p2 = Calc_pressao(x_v2, y_v2, z_v2, transdutor); */
      /*      p3 = Calc_pressao(x_v3, y_v3, z_v3, transdutor); */
      /*      potential_v1 = p1./(j*omega*rho0); */
      /*      potential_v2 = p2./(j*omega*rho0); */
      /*      potential_v3 = p3./(j*omega*rho0); */
      Calc_pressao(xs->data[n + xs->size[0] * nm1d2] + transdutor->disc / 2.0,
                   0.0, zs->data[n + zs->size[0] * nm1d2], transdutor->A,
                   transdutor->omega, transdutor->k, transdutor->a,
                   transdutor->posx, transdutor->posy, transdutor->posz,
                   transdutor->normal, transdutor->phase);
      Calc_pressao(xs->data[n + xs->size[0] * nm1d2] - transdutor->disc / 2.0,
                   0.0, zs->data[n + zs->size[0] * nm1d2], transdutor->A,
                   transdutor->omega, transdutor->k, transdutor->a,
                   transdutor->posx, transdutor->posy, transdutor->posz,
                   transdutor->normal, transdutor->phase);
      Calc_pressao(xs->data[n + xs->size[0] * nm1d2], transdutor->disc / 2.0,
                   zs->data[n + zs->size[0] * nm1d2], transdutor->A,
                   transdutor->omega, transdutor->k, transdutor->a,
                   transdutor->posx, transdutor->posy, transdutor->posz,
                   transdutor->normal, transdutor->phase);
      Calc_pressao(xs->data[n + xs->size[0] * nm1d2], 0.0 - transdutor->disc /
                   2.0, zs->data[n + zs->size[0] * nm1d2], transdutor->A,
                   transdutor->omega, transdutor->k, transdutor->a,
                   transdutor->posx, transdutor->posy, transdutor->posz,
                   transdutor->normal, transdutor->phase);
      Calc_pressao(xs->data[n + xs->size[0] * nm1d2], 0.0, zs->data[n + zs->
                   size[0] * nm1d2] + transdutor->disc / 2.0, transdutor->A,
                   transdutor->omega, transdutor->k, transdutor->a,
                   transdutor->posx, transdutor->posy, transdutor->posz,
                   transdutor->normal, transdutor->phase);
      Calc_pressao(xs->data[n + xs->size[0] * nm1d2], 0.0, zs->data[n + zs->
                   size[0] * nm1d2] - transdutor->disc / 2.0, transdutor->A,
                   transdutor->omega, transdutor->k, transdutor->a,
                   transdutor->posx, transdutor->posy, transdutor->posz,
                   transdutor->normal, transdutor->phase);
      Calc_pressao(xs->data[n + xs->size[0] * nm1d2], 0.0, zs->data[n + zs->
                   size[0] * nm1d2], transdutor->A, transdutor->omega,
                   transdutor->k, transdutor->a, transdutor->posx,
                   transdutor->posy, transdutor->posz, transdutor->normal,
                   transdutor->phase);
    }
  }

  emxFree_real_T(&zs);
  emxFree_real_T(&xs);
}

/*
 * File trailer for Plot_Uz.c
 *
 * [EOF]
 */
