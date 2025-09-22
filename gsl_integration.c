#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "glafic.h"

double gsl_qgaus_func(double x, void *p);
double gsl_qgaus2_func(double x, void *p);
static double (*func_qgaus_sav)(double);
static double (*func_qgaus2_sav)(double);

double gsl_romberg1_func(double x, void *p);
double gsl_romberg2_func(double x, void *p);
double gsl_romberg3_func(double x, void *p);
static double (*func_romberg1_sav)(double);
static double (*func_romberg2_sav)(double);
static double (*func_romberg3_sav)(double);

/*--------------------------------------------------------------
  Gauss integral, fixed points (21 points)
*/

double gsl_qgaus_func(double x, void *p)
{
  return (*func_qgaus_sav)(x);
}

double gsl_qgaus(double (*func)(double), double a, double b)
{
  double integral, err;
  size_t neval;
  gsl_function F;

  func_qgaus_sav = func;
  F.function = &gsl_qgaus_func;
  
  gsl_integration_qng(&F, a, b, TOL_QGAUS, TOL_QGAUS, &integral, &err, &neval);

  return integral;
}

double gsl_qgaus2_func(double x, void *p)
{
  return (*func_qgaus2_sav)(x);
}

double gsl_qgaus2(double (*func)(double), double a, double b)
{
  double integral, err;
  size_t neval;
  gsl_function F;

  func_qgaus2_sav = func;
  F.function = &gsl_qgaus2_func;
  
  gsl_integration_qng(&F, a, b, TOL_QGAUS, TOL_QGAUS, &integral, &err, &neval);

  return integral;
}

/*--------------------------------------------------------------
  Romberg integral
*/

double gsl_romberg1_func(double x, void *p)
{
  return (*func_romberg1_sav)(x);
}

double gsl_romberg2_func(double x, void *p)
{
  return (*func_romberg2_sav)(x);
}

double gsl_romberg3_func(double x, void *p)
{
  return (*func_romberg3_sav)(x);
}

double gsl_romberg1(double (*func)(double), double a, double b, double eps)
{
  int gslstatus;
  double integral;
  size_t neval;
  gsl_integration_romberg_workspace *workspace = NULL; 
  gsl_function F;

  func_romberg1_sav = func;
  F.function = &gsl_romberg1_func;

  workspace = gsl_integration_romberg_alloc( GSL_ROMBERG_N ); 
  
  if(workspace == NULL){
    fprintf(stderr, "failed to create workspace in gsl_romberg1\n");
    exit(EXIT_FAILURE);
  } else {
    gslstatus = gsl_integration_romberg(&F, a, b, 0.0, eps, &integral, &neval, workspace);

    /* if (gslstatus != GSL_SUCCESS){ 
       fprintf(stderr, "integration failed in gsl_romberg1\n");
       exit(EXIT_FAILURE);
       } */
  }

  gsl_integration_romberg_free(workspace); 

  return integral;
}

double gsl_romberg2(double (*func)(double), double a, double b, double eps)
{
  int gslstatus;
  double integral;
  size_t neval;
  gsl_integration_romberg_workspace *workspace = NULL; 
  gsl_function F;

  func_romberg2_sav = func;
  F.function = &gsl_romberg2_func;

  workspace = gsl_integration_romberg_alloc( GSL_ROMBERG_N ); 
  
  if(workspace == NULL){
    fprintf(stderr, "failed to create workspace in gsl_romberg2\n");
    exit(EXIT_FAILURE);
  } else {
    gslstatus = gsl_integration_romberg(&F, a, b, 0.0, eps, &integral, &neval, workspace); 
    /* if (gslstatus != GSL_SUCCESS){
       fprintf(stderr, "integration failed in gsl_romberg2\n");
       exit(EXIT_FAILURE);
       } */
  }

  gsl_integration_romberg_free(workspace); 

  return integral;
}

double gsl_romberg3(double (*func)(double), double a, double b, double eps)
{
  int gslstatus;
  double integral;
  size_t neval;
  gsl_integration_romberg_workspace *workspace = NULL; 
  gsl_function F;

  func_romberg3_sav = func;
  F.function = &gsl_romberg3_func;

  workspace = gsl_integration_romberg_alloc( GSL_ROMBERG_N ); 
  
  if(workspace == NULL){
    fprintf(stderr, "failed to create workspace in gsl_romberg3\n");
    exit(EXIT_FAILURE);
  } else {
    gslstatus = gsl_integration_romberg(&F, a, b, 0.0, eps, &integral, &neval, workspace); 
    /* if (gslstatus != GSL_SUCCESS){
       fprintf(stderr, "integration failed in gsl_romberg3\n");
       exit(EXIT_FAILURE);
      } */
  }

  gsl_integration_romberg_free(workspace); 

  return integral;
}

