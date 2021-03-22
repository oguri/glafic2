#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "glafic.h"

double zbrent_func(double x, void *p);
static double (*func_sav)(double);

double gsl_zbrent(double (*func)(double), double x_lo, double x_hi, double tol)
{
  int status;
  int iter = 0;
  const gsl_root_fsolver_type *solver_type;
  gsl_root_fsolver *solver;
  gsl_function F;
  double r;

  func_sav = func;
  F.function = &zbrent_func;

  solver_type = gsl_root_fsolver_brent;
  solver = gsl_root_fsolver_alloc(solver_type);
  gsl_root_fsolver_set(solver, &F, x_lo, x_hi);
 
  do{
    iter++;
    status = gsl_root_fsolver_iterate(solver);
    r = gsl_root_fsolver_root(solver);
    x_lo = gsl_root_fsolver_x_lower(solver);
    x_hi = gsl_root_fsolver_x_upper(solver);
    status = gsl_root_test_interval(x_lo, x_hi, 0, tol);
  }while((status == GSL_CONTINUE) && (iter < GSL_ZBRENT_MAXITER));
  
  gsl_root_fsolver_free(solver);
  
  return r;
}

double zbrent_func(double x, void *p)
{
  return (*func_sav)(x);
}
