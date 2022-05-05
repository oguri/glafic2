#define GLOBAL_SET
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../glafic.h"

void test_example(void);
void test_quadhost(void);
void test_point(void);
void test_readexample(void);

/*--------------------------------------------------------------
   main function
*/

int main(int argc, char **argv)
{
  test_example();
  test_quadhost();
  test_point();
  test_readexample();
}

void test_example(void)
{
  int verb = 0;
  int i, j, ni;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  double pout[NPAR_LMODEL];

  printf("### test_example\n");
  
  glafic_init(0.3, 0.7, -1.0, 0.7, "out", -60.0, -60.0, 60.0, 60.0, 0.2, 3.0, 5, 0, verb);

  glafic_set_secondary("chi2_splane    1", verb);

  glafic_startup_setnum(2, 2, 1);

  glafic_set_lens(1, "anfw", 0.3, 7.2e14, 0.0, 0.0, 0.5, -45.0, 6.0, 0.0);
  glafic_set_lens(2, "sie", 0.5, 300.0,  2.0, 2.0, 0.2, -20.0, 0.02, 0.0);
  glafic_set_extend(1, "sersic", 1.5, 150.0, -1.0, -1.5, 0.3, 90.0, 0.8, 1.0);
  glafic_set_extend(2, "gauss", 2.0, 150.0,  1.2,  1.0, 0.2, 10.0, 0.6, 0.0);
  glafic_set_point(1, 2.5, 1.0, 0.5);

  glafic_model_init(verb);

  glafic_point_solve(2.5, 1.0, 0.5, &ni, rr, verb);
  
  for(i=0;i<ni;i++){
    printf("image %d: %13e %13e %13e %13e\n", i + 1, rr[i][0], rr[i][1], rr[i][2], rr[i][3]);
  }

  glafic_calcimage(2.5, 1.0, -1.5, pout, 0, verb);
  printf("kappa = %e\n", pout[3]);
  printf("ein(id:1) = %e \n", glafic_calcein_i(1, 2.5));
  printf("ein(id:2) = %e \n", glafic_calcein_i(2, 2.5));
  printf("ein2 = %e\n",glafic_calcein2(0, 2.5, 0.0, 0.0));
  printf("kappa_ave = %e\n",glafic_kappa_ave(0, 2.5, 16.0, 0.0, 0.0));
  printf("kappa_cum = %e\n",glafic_kappa_cum(0, 2.5, 16.0, 0.0, 0.0));


  // glafic_findimg();
  // glafic_writelens(2.5);
  // glafic_writecrit(2.5);
  // glafic_writemesh(2.5);
  // glafic_lenscenter(2.5);

  glafic_set_array_extend(0, 0.0, 0.0, 0);
  i = 300;
  j = 250;
  printf("extend(i=%d,j=%d) = %e\n", i, j, glafic_extend_array_ij(i, j));
  glafic_unset_array_extend();

  glafic_quit();
}

void test_quadhost(void)
{
  int verb = 0;

  printf("### test_quadhost\n");

  glafic_init(0.3, 0.7, -1.0, 0.7, "out", -5.0, -5.0, 5.0, 5.0, 0.1, 0.5, 6, 0, verb);

  glafic_set_secondary("chi2_restart   -1", verb);
  glafic_set_secondary("obs_gain       3.0", verb);
  glafic_set_secondary("obs_ncomb      10", verb);
  glafic_set_secondary("psfconv_size   3.0", verb);
  glafic_set_secondary("flag_extnorm   1", verb);

  glafic_startup_setnum(2, 3, 0);

  glafic_set_lens(1, "sie", 0.3, 315.0, 0.0, 0.0, 0.3,  37.0, 0.0, 0.0);
  glafic_set_lens(2, "pert", 0.3,  2.0, 0.0, 0.0, 0.042, 60.0, 0.0, 0.0);
  glafic_set_extend(1,"sersic", 0.3, 4.9e5,  0.0,  0.0, 0.3,   37.0, 1.0, 4.0);
  glafic_set_extend(2, "sersic", 2.0, 1.2e4,  0.21, -0.09, 0.3, -10.0, 0.5, 2.0);
  glafic_set_extend(3, "point",  2.0, 2.1e4,  0.21, -0.09, 0.0,   0.0, 0.0, 0.0);
  glafic_set_psf(0.22, 0.07, 15.0, 3.0, 0.8, 0.04, -20.0, 2.5, 0.35);

  // glafic_setopt_lens(1, 0, 1, 0, 0, 0, 0, 0, 0);
  // glafic_setopt_lens(2, 0, 0, 0, 0, 1, 1, 0, 0); 
  // glafic_setopt_extend(1, 0, 1, 0, 0, 0, 0, 1, 0); 
  // glafic_setopt_extend(2, 0, 1, 1, 1, 0, 0, 1, 0); 
  // glafic_setopt_extend(3, 0, 1, 1, 1, 0, 0, 0, 0); 
  // glafic_setopt_psf(1, 0, 0, 0, 1, 0, 0, 0, 1); 

  glafic_model_init(verb);

  glafic_readobs_extend("../samples/obs_quadhost.fits", "N", verb);
  // glafic_optimize(verb);
  printf("chi2 = %e\n", glafic_c2calc());
  
  glafic_set_lens(1, "sie", 0.3, 320.0, 0.0, 0.0, 0.3,  37.0, 0.0, 0.0);
  glafic_model_init(verb);

  printf("chi2 = %e\n", glafic_c2calc());

  glafic_quit();
}

void test_point(void)
{
  int verb = 0;
  
  printf("### test_point\n");

  glafic_init(0.3, 0.7, -1.0, 0.7, "out", -5.0, -5.0, 5.0, 5.0, 0.02, 0.5, 5, 0, verb);

  glafic_set_secondary("chi2_splane    0", verb);
  glafic_set_secondary("chi2_checknimg 1", verb);
  glafic_set_secondary("chi2_restart   -1", verb);
  glafic_set_secondary("chi2_usemag    0", verb);
  glafic_set_secondary("hvary          1", verb);
  
  glafic_startup_setnum(2, 0, 1);

  glafic_set_lens(1, "sie", 0.5, 300.0,  0.0, 0.0, 0.35,   0.0, 0.0, 0.0);
  glafic_set_lens(2, "pert", 0.5,   2.0,  0.0, 0.0, 0.05,  60.0, 0.0, 0.0);
  glafic_set_point(1, 2.0, -0.15, 0.05);

  glafic_setopt_lens(1, 0, 1, 1, 1, 1, 1, 0, 0);
  glafic_setopt_lens(2, 0, 0, 0, 0, 1, 1, 0, 0);
  glafic_setopt_point(1, 0, 1, 1);

  glafic_model_init(verb);

  glafic_readobs_point("../samples/obs_point.dat", verb);
  glafic_parprior("../samples/prior_point.dat", verb);
  // glafic_optimize(verb);
  printf("chi2 = %e\n", glafic_c2calc());

  glafic_quit();
}

void test_readexample(void)
{
  int verb = 0;
  int i, j, ni;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  double pout[NPAR_LMODEL];

  printf("### test_readexample\n");

  glafic_init_file("../samples/example.input", 0);

  glafic_point_solve(2.5, 1.0, 0.5, &ni, rr, verb);
  
  for(i=0;i<ni;i++){
    printf("image %d: %13e %13e %13e %13e\n", i + 1, rr[i][0], rr[i][1], rr[i][2], rr[i][3]);
  }

  glafic_calcimage(2.5, 1.0, -1.5, pout, 0, verb);
  printf("kappa = %e\n", pout[3]);
  printf("ein(id:1) = %e \n", glafic_calcein_i(1, 2.5));
  printf("ein(id:2) = %e \n", glafic_calcein_i(2, 2.5));
  printf("ein2 = %e\n",glafic_calcein2(0, 2.5, 0.0, 0.0));
  printf("kappa_ave = %e\n",glafic_kappa_ave(0, 2.5, 16.0, 0.0, 0.0));
  printf("kappa_cum = %e\n",glafic_kappa_cum(0, 2.5, 16.0, 0.0, 0.0));


  // glafic_findimg();
  // glafic_writelens(2.5);
  // glafic_writecrit(2.5);
  // glafic_writemesh(2.5);
  // glafic_lenscenter(2.5);

  glafic_set_array_extend(0, 0.0, 0.0, 0);
  i = 300;
  j = 250;
  printf("extend(i=%d,j=%d) = %e\n", i, j, glafic_extend_array_ij(i, j));
  glafic_unset_array_extend();

  glafic_quit();
}

