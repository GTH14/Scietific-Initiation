/*
 * File: main.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 03-Jun-2021 14:42:05
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "rt_nonfinite.h"
#include "f_Simulacao_dinamica_v3.h"
#include "main.h"
#include "f_Simulacao_dinamica_v3_terminate.h"
#include "f_Simulacao_dinamica_v3_emxAPI.h"
#include "f_Simulacao_dinamica_v3_initialize.h"

/* Function Declarations */
static void argInit_1x3_real_T(double result[3]);
static void argInit_1x50001_real_T(double result[50001]);
static void argInit_1x5001_real_T(double result[5001]);
static void argInit_3x1_real_T(double result[3]);
static void argInit_64x1_real_T(double result[64]);
static double argInit_real_T(void);
static void argInit_struct0_T(struct0_T *result);
static void argInit_struct1_T(struct1_T *result);
static void main_f_Simulacao_dinamica_v3(void);

/* Function Definitions */

/*
 * Arguments    : double result[3]
 * Return Type  : void
 */
static void argInit_1x3_real_T(double result[3])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 3; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[50001]
 * Return Type  : void
 */
static void argInit_1x50001_real_T(double result[50001])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 50001; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[5001]
 * Return Type  : void
 */
static void argInit_1x5001_real_T(double result[5001])
{
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx1 = 0; idx1 < 5001; idx1++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx1] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[3]
 * Return Type  : void
 */
static void argInit_3x1_real_T(double result[3])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[64]
 * Return Type  : void
 */
static void argInit_64x1_real_T(double result[64])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 64; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : struct0_T *result
 * Return Type  : void
 */
static void argInit_struct0_T(struct0_T *result)
{
  /* Set the value of each structure field.
     Change this value to the value that the application requires. */
  result->freq = argInit_real_T();
  result->c = argInit_real_T();
  result->A = argInit_real_T();
  result->rho0 = argInit_real_T();
  result->lambda = argInit_real_T();
  result->omega = argInit_real_T();
  result->k = argInit_real_T();
  result->a = argInit_real_T();
  result->disc = argInit_real_T();
  argInit_64x1_real_T(result->posx);
  argInit_64x1_real_T(result->posy);
  argInit_64x1_real_T(result->posz);
  argInit_1x3_real_T(result->normal);
  argInit_3x1_real_T(result->foco);
  argInit_64x1_real_T(result->phase);
  result->R = argInit_real_T();
}

/*
 * Arguments    : struct1_T *result
 * Return Type  : void
 */
static void argInit_struct1_T(struct1_T *result)
{
  /* Set the value of each structure field.
     Change this value to the value that the application requires. */
  result->rho = argInit_real_T();
  result->R = argInit_real_T();
  argInit_3x1_real_T(result->pos_i);
  argInit_3x1_real_T(result->v_i);
  result->h = argInit_real_T();
  argInit_1x5001_real_T(result->t_v);
  argInit_1x50001_real_T(result->x);
  argInit_1x50001_real_T(result->y);
  argInit_1x50001_real_T(result->z);
  argInit_1x50001_real_T(result->vx);
  argInit_1x50001_real_T(result->vy);
  argInit_1x50001_real_T(result->vz);
  result->ti = argInit_real_T();
  result->tf = argInit_real_T();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_f_Simulacao_dinamica_v3(void)
{
  emxArray_creal_T *x;
  struct0_T r0;
  static struct1_T r1;
  emxInitArray_creal_T(&x, 2);

  /* Initialize function 'f_Simulacao_dinamica_v3' input arguments. */
  /* Initialize function input argument 'transdutor'. */
  /* Initialize function input argument 'bolinha'. */
  /* Call the entry-point 'f_Simulacao_dinamica_v3'. */
  argInit_struct0_T(&r0);
  argInit_struct1_T(&r1);
  f_Simulacao_dinamica_v3(&r0, &r1, x);
  emxDestroyArray_creal_T(x);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  f_Simulacao_dinamica_v3_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_f_Simulacao_dinamica_v3();

  /* Terminate the application.
     You do not need to do this more than one time. */
  f_Simulacao_dinamica_v3_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
