#define GLOBAL_SET
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../glafic.h"
#include <gsl/gsl_randist.h>

void do_mock(void);

double gene_e(void);
double gene_gam(void);
double gene_ang(void);

gsl_rng *ran;

/*--------------------------------------------------------------
   main function
*/

int main(int argc, char **argv)
{
  ran = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(ran, -342116);

  do_mock();

  return 0;
}

/*--------------------------------------------------------------
  generate mock
*/

void do_mock(void)
{
  int i, j, k, n_mock, mlev, f, ni;
  double zl, zs, sig, m, ein, rt_range, pix;
  double e, e_pa, g, g_pa, xs, ys, s, sep;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  FILE *fptr_w1;
  FILE *fptr_w2;

  /* number of trials */
  n_mock = 100000;
  
  /* fixed lens and source parameters so that theta_Ein = 1 arcsec */
  zl = 0.5;
  zs = 2.0;
  sig = 233.605333;

  /* dummy values */
  f = 0;
  m = 99.0;
  
  /* parameters for glafic */
  rt_range = 4.0;
  pix = 0.2;
  mlev = 5;

  /* write results to these files */
  fptr_w1 = fopen("out_siex_result.dat", "w");
  fptr_w2 = fopen("out_siex_log.dat", "w");
  
  /* run glafic */
  glafic_init(0.3, 0.7, -1.0, 0.7, "out", (-1.0) * 2.0 * rt_range, (-1.0) * 2.0 * rt_range, 2.0 * rt_range, 2.0 * rt_range, pix, pix, mlev, 0, 0);

  glafic_startup_setnum(2, 0, 1);

  glafic_set_lens(1, "sie",  zl, sig, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  glafic_set_lens(2, "pert", zl,  zs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  glafic_set_point(1, zs, 0.0, 0.0);

  glafic_model_init(0);
  ein = glafic_calcein_i(1, zs);

  for(i=0;i<n_mock;i++){
    e = gene_e();
    e_pa = gene_ang();
    g = gene_gam();
    g_pa = gene_ang();
    
    xs = (gsl_rng_uniform(ran) - 0.5) * rt_range * 2.0;
    ys = (gsl_rng_uniform(ran) - 0.5) * rt_range * 2.0;
    
    glafic_set_lens(1, "sie",  zl, sig, 0.0, 0.0, e, e_pa, 0.0, 0.0);
    glafic_set_lens(2, "pert", zl,  zs, 0.0, 0.0, g, g_pa, 0.0, 0.0);
    
    glafic_model_init(0);
    glafic_point_solve(zs, xs, ys, &ni, rr, 0);

    if(ni > 1){
      sep = 0.0;
      for(j=0;j<ni-1;j++){
	for(k=ni-1;k>j;k--){
	  s = (rr[j][0] - rr[k][0]) * (rr[j][0] - rr[k][0]) + (rr[j][1] - rr[k][1]) * (rr[j][1] - rr[k][1]);

	  if(s > sep){ sep = s; }
	}	
      }
      sep = sqrt(sep);
      
      fprintf(fptr_w1,"%d %e %e %e %e %e %e %13e %13e %13e %13e %13e %13e %d %d\n", ni, zl, sig, zs, m, m, sep, e, e_pa, g, g_pa, xs, ys, f, i);
      fflush(stdout);

      fprintf(fptr_w2, "%d %13e %13e %d %d\n", ni, xs, ys, f, i);
      for(j=0;j<ni;j++){
	fprintf(fptr_w2," %13e %13e %13e %13e\n", rr[j][0], rr[j][1], rr[j][2], rr[j][3]);
      }
    }
  }
  
  fclose(fptr_w1);
  fclose(fptr_w2);
  
  glafic_quit();
    
  return;
}

/*--------------------------------------------------------------
  generate variables
*/

double gene_e(void)
{
  double e, em, se;

  em = 0.3;
  se = 0.16;

  do{
    e = gsl_ran_gaussian(ran, 1.0) * se + em;
  }while((e < 0.0) || (e > 0.9));

  return e;
}

double gene_gam(void)
{
  double lg, lgm, slg, g;

  lgm = log10(0.05);
  slg = 0.2;

  do{
    lg = gsl_ran_gaussian(ran, 1.0) * slg + lgm;
    g = pow(10.0, lg);
  }while(g > 0.9);

  return g;
}

double gene_ang(void)
{
  return (gsl_rng_uniform(ran) - 0.5) * 360.0;
}

