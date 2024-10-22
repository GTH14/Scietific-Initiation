/*
 * File: main.c
 *
 * MATLAB Coder version            : 3.4
 * C/C++ source code generated on  : 17-Jun-2021 12:47:38
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
#include "Plot_Uz.h"
#include "main.h"
#include "Plot_Uz_terminate.h"
#include "Plot_Uz_initialize.h"

/* Function Declarations */
static void argInit_1x3_real_T(double result[3]);
static void argInit_3x1_real_T(double result[3]);
static void argInit_64x1_real_T(double result[64]);
static double argInit_real_T(void);
static void argInit_struct0_T(struct0_T *result);
static void main_Plot_Uz(void);

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
  result->R = argInit_real_T();
  result->nu = argInit_real_T();
  argInit_3x1_real_T(result->foco);
  argInit_64x1_real_T(result->phase);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_Plot_Uz(void)
{
  struct0_T r0;

  /* Initialize function 'Plot_Uz' input arguments. */
  /* Initialize function input argument 'transdutor'. */
  /* Call the entry-point 'Plot_Uz'. */
  argInit_struct0_T(&r0);
  Plot_Uz(&r0);
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
  Plot_Uz_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_Plot_Uz();

  /* Terminate the application.
     You do not need to do this more than one time. */
  Plot_Uz_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
