#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "glafic.h"

static int mcmc_calc_ndim;
static double mcmc_calc_sigma[NDIMMAX + 1];

/*--------------------------------------------------------------
  Markov Chain Monte Carlo
*/

void mcmc_calc(int n)
{
  int i, j, k;
  double c0, c1, c10, p10, xx;
  double par[NDIMMAX + 1], par2[NDIMMAX + 1];
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  fprintf(stderr, "######## Markov Chain Monte Carlo\n");
  fprintf(stderr, " number of evaluations = %d\n", n);

  sprintf(fname, "%s_mcmc.dat", file_prefix);
  fprintf(stderr, " output file name = %s\n\n", fname);
    
  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (mcmc)");

  opt_lens_static(0);
  mcmc_calc_ndim = opt_lens_calcndim();

  if(mcmc_calc_ndim > NDIMMAX) terminator("too many lens parameters");

  paratopar(par);

  c0 = chi2calc(par);
  /* p0 = exp((-0.5) * c0); */
  paratopar(par);

  if(flag_mcmcall != 0){ fprintf(fptr, "accept "); }
  fprintf(fptr, "%e ", c0);
  for(j=1;j<=mcmc_calc_ndim;j++){
    fprintf(fptr, "%13e ", par[j]);
  }
  fprintf(fptr, "\n");

  k = 0;
  for(i=0;i<n;i++){
    mcmc_calc_step(par, par2);
    c1 = chi2calc(par2);
    /* p1 = exp((-0.5) * c1); */
    paratopar(par2);
    c10 = c1 - c0;
    if(c10 > CUT_CHI2_EXP){ c10 = CUT_CHI2_EXP; }
    if(c10 < ((-1.0) * CUT_CHI2_EXP) ){ c10 = (-1.0) * CUT_CHI2_EXP; }
    p10 = exp((-0.5) * c10);

    /* fprintf(stderr, "%8d: %e %e ", i + 1, p0, p1); */
    fprintf(stderr, "%8d: %e %e %e ", i + 1, c0, c1, p10);
    xx = gsl_rng_uniform(ran_gsl);

    if(xx < p10){
      if(flag_mcmcall != 0){ fprintf(fptr, "accept "); }
      fprintf(fptr, "%e ", c1);
      for(j=1;j<=mcmc_calc_ndim;j++){
	par[j] = par2[j];
	fprintf(fptr, "%13e ", par[j]);
      }
      c0 = c1;
      /* p0 = p1; */
      fprintf(fptr, "\n");
      fflush(fptr);
      k++;
    } else {
      if(flag_mcmcall != 0){ 
	fprintf(fptr, "reject ");
	fprintf(fptr, "%e ", c1);
	for(j=1;j<=mcmc_calc_ndim;j++){
	  fprintf(fptr, "%13e ", par2[j]);
	}
	fprintf(fptr, "\n");
	fflush(fptr);
      }
      fprintf(stderr, "[rejected] ");
    }
    fprintf(stderr, "\n");
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Acceptance rate: %5.3f\n\n", ((double)k) / ((double)n));
  fclose(fptr);
  
  return;
}

/* {lens} -> {cosmo} -> {ext} -> sky -> {point} -> {psf} */
void paratopar(double par[NDIMMAX + 1])
{
  int i, j, ndim, nn;
  
  ndim = 0;
  for(i=0;i<num_len;i++){
    for(j=0;j<NPAR_LEN;j++){
      if(flag_para_lens[i][j] == 1){
	ndim++;
	par[ndim] = para_lens[i][j];
      }
    }
  }
  if(ovary != 0){ ndim++; par[ndim] = omega;  }
  if(lvary != 0){ ndim++; par[ndim] = lambda; }
  if(wvary != 0){ ndim++; par[ndim] = weos;   }
  if(hvary != 0){ ndim++; par[ndim] = hubble; }

  nn = 0;
  for(i=0;i<num_ext;i++){
    for(j=0;j<NPAR_EXT;j++){
      if(flag_para_ext[i][j] == 1){
	nn++;
	ndim++;
	par[ndim] = para_ext[i][j];
      }
    }
  }
  if((skyfix == 0) && (nn > 0)){
    ndim++;
    par[ndim] = skymed;
  }
  
  for(i=0;i<num_poi;i++){
    if(flag_para_poi[i][0] == 1){
      ndim++;
      par[ndim] = para_poi[i][0];
    }
  }

  for(j=0;j<NPAR_PSF;j++){
    if((nn > 0) && (flag_seeing == 1) && (flag_para_psf[j] == 1)){
      ndim++;
      par[ndim] = para_psf[j];
    }
  }
  if(ndim != mcmc_calc_ndim) terminator("invalid number of MCMC parameters");

  return;
}

void mcmc_calc_step(double par[NDIMMAX + 1], double par2[NDIMMAX + 1])
{
  int i;

  for(i=1;i<=mcmc_calc_ndim;i++){
    par2[i] = mcmc_calc_par(i, par[i]);
  }

  return;
}

double mcmc_calc_par(int i, double x0)
{
  double x;

  if(mcmc_calc_sigma[i] >= 0.0){
    x = x0 + gsl_ran_gaussian(ran_gsl, 1.0) * mcmc_calc_sigma[i];
  } else {
    x = x0 * pow(10.0, gsl_ran_gaussian(ran_gsl, 1.0) * fabs(mcmc_calc_sigma[i]));
  }

  return x;
}

/*--------------------------------------------------------------
  set sigma for proposal function
*/

void mcmc_calc_init(char *infile)
{
  int i, j, f, n, nn;
  char buffer[INPUT_MAXCHAR];
  FILE* fptr;
  
  fprintf(stderr, "######## reading sigma file for Markov Chain Monte Carlo\n");
  fprintf(stderr, " input file name = %s \n\n", infile);

  fptr = fopen(infile, "r");
  
  if(fptr == NULL) terminator("failed at fopen (mcmc_sigma)");

  n = opt_lens_calcndim();
  
  f = 0;
  j = 0;
  while((fgets(buffer, INPUT_MAXCHAR, fptr)) && (f == 0)){
    if(buffer[0] != '#'){
      nn = sscanf(buffer, "%d", &mcmc_calc_ndim);
      if(mcmc_calc_ndim != n) terminator("invalid number of parameters (mcmc_sigma)");
      if(mcmc_calc_ndim > NDIMMAX) terminator("too many MCMC parameters");
      if (nn != EOF){
	f = 1;
	for(i=1;i<=mcmc_calc_ndim;i++){
	  if(fgets(buffer, INPUT_MAXCHAR, fptr)){
	    nn = sscanf(buffer, "%lf", &mcmc_calc_sigma[i]);
	    if(nn == 1) j++;
	  }
	}
	if(mcmc_calc_ndim != j) terminator("invalid format (mcmc_sigma)");
      }
    }
  }

  fprintf(stderr, "read %d sigmas\n\n", j);

  fclose(fptr);

  return;
}

/*--------------------------------------------------------------
  read chains, convert it to some quantities
*/

void mcmc_out_kappa_rad(char *infile, double nd, double zs, double x0, double y0, double r1, double r2, int n, int lensid)
{
  int i, f, nn, ndim, flag;
  double c2;
  double kapbin[NMAX_KAPBIN + 1], rbin[NMAX_KAPBIN + 1];
  double par[NDIMMAX + 1];
  char fname[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  FILE* fptr_r;
  FILE* fptr;

  fprintf(stderr, "######## kappa radial profiles from mcmc chain\n");
  fprintf(stderr, " zs = %e,  lens id = %d,  step = %d\n", zs, lensid, n);
  fprintf(stderr, " center = (%e, %e), range = [%e, %e]\n\n", x0, y0, r1, r2);
  sprintf(fname, "%s_mcmc_kaprad.dat", file_prefix);
  fprintf(stderr, " input file name  = %s\n", infile);
  fprintf(stderr, " outfut file name = %s\n\n", fname);

  if(n < 0){ nn = n * (-1); } else { nn = n; }

  ndim = opt_lens_calcndim();

  if(nd != ndim) terminator("invalid number of MCMC parameters");

  fptr_r = fopen(infile, "r");
  fptr = fopen(fname, "w");
  if((fptr_r == NULL) || (fptr == NULL)) terminator("failed at fopen (mcmc_kaprad)");

  flag = 0;
  /* check file type */
  fscanf(fptr_r, "%s ", keyword);
  if(strcmp(keyword, "accept") == 0){ flag = 1;
  } else if(strcmp(keyword, "reject") == 0){ flag = -1;
  } else { rewind(fptr_r); }

  f = 0;
  while(fscanf(fptr_r, "%lf ", &c2) != EOF){
    for(i=1;i<=nd;i++){
      fscanf(fptr_r, "%lf ", &par[i]);
    }
    if(flag >= 0){
      partopara(par);

      kappa_rad(zs, x0, y0, r1, r2, n, lensid, kapbin, rbin, 0);
      if(f == 0){
	f = 1;
	fprintf(fptr, "# zs = %e,  x0 = %e,  y0 = %e,  id = %d \n# r = ", zs, x0, y0, lensid);
	for(i=0;i<=nn;i++) fprintf(fptr, "%e ", rbin[i]);
	fprintf(fptr, "\n");
      }
      
      fprintf(fptr, "%e ", c2);
      for(i=0;i<=nn;i++) fprintf(fptr, "%e ", kapbin[i]);
      fprintf(fptr, "\n");
      fflush(fptr);
    }
    if(flag != 0){
      fscanf(fptr_r, "%s ", keyword);
      if(strcmp(keyword, "accept") == 0){ flag = 1;
      } else if(strcmp(keyword, "reject") == 0){ flag = -1; } 
    }
  }

  fclose(fptr_r);
  fclose(fptr);

}

void mcmc_out_ein(char *infile, double nd, double zs, int lensid)
{
  int i, ndim, flag;
  double c2, ein;
  double par[NDIMMAX + 1];
  char fname[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  FILE* fptr_r;
  FILE* fptr;

  fprintf(stderr, "######## Einstein radii from mcmc chain\n");
  fprintf(stderr, " zs = %e,  lens id = %d\n\n", zs, lensid);
  sprintf(fname, "%s_mcmc_ein.dat", file_prefix);
  fprintf(stderr, " input file name  = %s\n", infile);
  fprintf(stderr, " outfut file name = %s\n\n", fname);

  set_distance_lpl_zs(zs);
  
  ndim = opt_lens_calcndim();

  if(nd != ndim) terminator("invalid number of MCMC parameters");

  fptr_r = fopen(infile, "r");
  fptr = fopen(fname, "w");
  if((fptr_r == NULL) || (fptr == NULL)) terminator("failed at fopen (mcmc_ein)");

  flag = 0;
  /* check file type */
  fscanf(fptr_r, "%s ", keyword);
  if(strcmp(keyword, "accept") == 0){ flag = 1;
  } else if(strcmp(keyword, "reject") == 0){ flag = -1;
  } else { rewind(fptr_r); }

  fprintf(fptr, "# zs = %e,  id = %d \n", zs, lensid);
  while(fscanf(fptr_r, "%lf ", &c2) != EOF){
    for(i=1;i<=nd;i++){
      fscanf(fptr_r, "%lf ", &par[i]);
    }
    if(flag >= 0){
      partopara(par);

      fprintf(fptr, "%e ", c2);
      
      for(i=0;i<num_len;i++){
	if((lensid == (i + 1)) || (lensid == 0)){
	  ein = calcein_i(i);
	  if(ein < 0.0){ ein = 0.0; }
	  fprintf(fptr, "%e ", ein);
	}
      }
      fprintf(fptr, "\n");
      fflush(fptr);
    }
    if(flag != 0){
      fscanf(fptr_r, "%s ", keyword);
      if(strcmp(keyword, "accept") == 0){ flag = 1;
      } else if(strcmp(keyword, "reject") == 0){ flag = -1; } 
    }
  }

  fclose(fptr_r);
  fclose(fptr);

}

void mcmc_out_ein2(char *infile, double nd, double zs, double x0, double y0, int lensid)
{
  int i, ndim, flag;
  double c2, ein;
  double par[NDIMMAX + 1];
  char fname[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  FILE* fptr_r;
  FILE* fptr;

  fprintf(stderr, "######## Einstein radii from mcmc chain (from average kappa)\n");
  fprintf(stderr, " zs = %e,  lens id = %d\n", zs, lensid);
  fprintf(stderr, " center = (%e, %e)\n\n", x0, y0);
  sprintf(fname, "%s_mcmc_ein.dat", file_prefix);
  fprintf(stderr, " input file name  = %s\n", infile);
  fprintf(stderr, " outfut file name = %s\n\n", fname);

  set_distance_lpl_zs(zs);
  
  ndim = opt_lens_calcndim();

  if(nd != ndim) terminator("invalid number of MCMC parameters");

  fptr_r = fopen(infile, "r");
  fptr = fopen(fname, "w");
  if((fptr_r == NULL) || (fptr == NULL)) terminator("failed at fopen (mcmc_ein2)");

  flag = 0;
  /* check file type */
  fscanf(fptr_r, "%s ", keyword);
  if(strcmp(keyword, "accept") == 0){ flag = 1;
  } else if(strcmp(keyword, "reject") == 0){ flag = -1;
  } else { rewind(fptr_r); }

  fprintf(fptr, "# zs = %e, x0 = %e, y0 = %e, id = %d \n", zs, x0, y0, lensid);
  while(fscanf(fptr_r, "%lf ", &c2) != EOF){
    for(i=1;i<=nd;i++){
      fscanf(fptr_r, "%lf ", &par[i]);
    }
    if(flag >= 0){
      partopara(par);
      
      ein = calcein2_calc(zs, x0, y0, lensid);
      
      fprintf(fptr, "%e %e\n", c2, ein);
      fflush(fptr);
    }
    if(flag != 0){
      fscanf(fptr_r, "%s ", keyword);
      if(strcmp(keyword, "accept") == 0){ flag = 1;
      } else if(strcmp(keyword, "reject") == 0){ flag = -1; } 
    }
  }

  fclose(fptr_r);
  fclose(fptr);

}

void mcmc_out_calcim(char *infile, double nd, double zs, double x0, double y0)
{
  int i, ndim, flag;
  double c2;
  double par[NDIMMAX + 1];
  char fname[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  double pout[NPAR_LMODEL];
  FILE* fptr_r;
  FILE* fptr;

  fprintf(stderr, "######## image properties from mcmc chain\n");
  fprintf(stderr, " zs = %e,  x = %e,  y = %e\n", zs, x0, y0);
  sprintf(fname, "%s_mcmc_calcim.dat", file_prefix);
  fprintf(stderr, " input file name  = %s\n", infile);
  fprintf(stderr, " outfut file name = %s\n\n", fname);

  set_distance_lpl_zs(zs);
  
  ndim = opt_lens_calcndim();

  if(nd != ndim) terminator("invalid number of MCMC parameters");

  fptr_r = fopen(infile, "r");
  fptr = fopen(fname, "w");
  if((fptr_r == NULL) || (fptr == NULL)) terminator("failed at fopen (mcmc_calcim)");

  flag = 0;
  /* check file type */
  fscanf(fptr_r, "%s ", keyword);
  if(strcmp(keyword, "accept") == 0){ flag = 1;
  } else if(strcmp(keyword, "reject") == 0){ flag = -1;
  } else { rewind(fptr_r); }

  fprintf(fptr, "# zs = %e, x = %e, y = %e\n", zs, x0, y0);
  fprintf(fptr, "# chi^2 mag kappa gamma_1 gamma_2 gamma alpha_x alpha_y rotation\n");
  while(fscanf(fptr_r, "%lf ", &c2) != EOF){
    for(i=1;i<=nd;i++){
      fscanf(fptr_r, "%lf ", &par[i]);
    }
    if(flag >= 0){
      partopara(par);

      lensmodel(x0, y0, pout, 0, 0);

      fprintf(fptr, "%e %e %e %e %e %e %e %e %e\n", c2, 1.0 / (pout[6] + imag_ceil), pout[3], pout[4], pout[5], sqrt(pout[4] * pout[4] + pout[5] * pout[5]), pout[0], pout[1], pout[7]);
      fflush(fptr);
    }
     if(flag != 0){
      fscanf(fptr_r, "%s ", keyword);
      if(strcmp(keyword, "accept") == 0){ flag = 1;
      } else if(strcmp(keyword, "reject") == 0){ flag = -1; } 
    }
  }   
   
  fclose(fptr_r);
  fclose(fptr);
 
}
