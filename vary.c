#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "glafic.h"

/*--------------------------------------------------------------
  change parameters
*/

void varyone(int i, int j, double pmin, double pmax, int n, int flag)
{
  int ff, k, nn;
  double p, hh, c2;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  fprintf(stderr, "######## compute chi^2 with changing parameters\n");
  fprintf(stderr, " lens id = %d,  par no. = %d,  step = %d\n", i + 1, j + 1, n);
  
  sprintf(fname, "%s_vary.dat", file_prefix);
  fprintf(stderr, " output file name = %s\n\n", fname);
    
  if((i >= num_len) || (j >= NPAR_LEN) || (i < 0) || (j < 0))
    terminator("lens id irrelevant (varyone)"); 

  if((n < 0) && (pmin <= 0.0))
    terminator("parameter must be positive for log spacing (varyone)"); 
  
  if(n < 0){
    nn = n * (-1);
  } else if(n > 0){
    nn = n;
  } else {
    nn = n;
    terminator("invalid step number (varyone)"); 
  }

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (varyone)");

  ff = flag_para_lens[i][j];
  flag_para_lens[i][j] = 0;

  if(n > 0){
    hh = (pmax - pmin) / ((double)nn);
  } else {
    hh = log10(pmax / pmin) / ((double)nn);
  }

  chi2_dumplimit = C2DUMPLIMIT_SET;

  for(k=0;k<=nn;k++){
    if(n>0){
      p = pmin + ((double)k) * hh;
    } else {
      p = pmin * pow(10.0, ((double)k) * hh);
    }
    if(k == 0){ p = pmin; }
    if(k == nn){ p = pmax; }
    para_lens[i][j] = p;
    
    c2 = opt_lens(flag, -1);
    fprintf(stderr, "%e %e\n", p, c2);
    fprintf(fptr, "%e %e\n", p, c2);
    fflush(fptr);
  }
  
  fprintf(stderr, "\n");
  flag_para_lens[i][j] = ff;

  fclose(fptr);

  return;
}

void varytwo(int i1, int j1, double pmin1, double pmax1, int n1, int i2, int j2, double pmin2, double pmax2, int n2, int flag)
{
  int ff1, ff2, k, l, nn1, nn2;
  double p1, p2, hh1, hh2, c2;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  fprintf(stderr, "######## compute chi^2 with changing parameters\n");
  fprintf(stderr, "lens id = %d,  par no. = %d,  step = %d\n", i1 + 1, j1 + 1, n1);
  fprintf(stderr, "lens id = %d,  par no. = %d,  step = %d\n", i2 + 1, j2 + 1, n2);
  
  sprintf(fname, "%s_vary.dat", file_prefix);
  fprintf(stderr, " output file name = %s\n\n", fname);
    
  if((i1 >= num_len) || (j1 >= NPAR_SRC) || (i1 < 0) || (j1 < 0) || (i2 >= num_len) || (j2 >= NPAR_SRC) || (i2 < 0) || (j2 < 0))
    terminator("lens id irrelevant (varytwo)"); 

  if(((n1 < 0) && (pmin1 <= 0.0)) || ((n2 < 0) && (pmin2 <= 0.0)))
    terminator("parameter must be positive for log spacing (varytwo)"); 
  
  if(n1 < 0){
    nn1 = n1 * (-1);
  } else if(n1 > 0){
    nn1 = n1;
  } else {
    nn1 = n1;
    terminator("invalid step number (varytwo)"); 
  }

  if(n2 < 0){
    nn2 = n2 * (-1);
  } else if(n2 > 0){
    nn2 = n2;
  } else {
    nn2 = n2;
    terminator("invalid step number (varytwo)"); 
  }

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (varyone)");

  ff1 = flag_para_lens[i1][j1];
  ff2 = flag_para_lens[i2][j2];
  flag_para_lens[i1][j1] = 0;
  flag_para_lens[i2][j2] = 0;

  if(n1 > 0){
    hh1 = (pmax1 - pmin1) / ((double)nn1);
  } else {
    hh1 = log10(pmax1 / pmin1) / ((double)nn1);
  }

  if(n2 > 0){
    hh2 = (pmax2 - pmin2) / ((double)nn2);
  } else {
    hh2 = log10(pmax2 / pmin2) / ((double)nn2);
  }

  chi2_dumplimit = C2DUMPLIMIT_SET;

  for(k=0;k<=nn1;k++){
    if(n1 > 0){
      p1 = pmin1 + ((double)k) * hh1;
    } else {
      p1 = pmin1 * pow(10.0, ((double)k) * hh1);
    }
    if(k == 0){ p1 = pmin1; }
    if(k == nn1){ p1 = pmax1; }
    para_lens[i1][j1] = p1;
    
    for(l=0;l<=nn2;l++){
      if(n2 > 0){
	p2 = pmin2 + ((double)l) * hh2;
      } else {
	p2 = pmin2 * pow(10.0, ((double)l) * hh2);
      }
      if(l == 0){ p2 = pmin2; }
      if(l == nn2){ p2 = pmax2; }
      para_lens[i2][j2] = p2;
    
      c2 = opt_lens(flag, -1);
      fprintf(stderr, "%e %e %e\n", p1, p2, c2);
      fprintf(fptr, "%e %e %e\n", p1, p2, c2);
      fflush(fptr);
    }
  }
  
  fprintf(stderr, "\n");
  flag_para_lens[i1][j1] = ff1;
  flag_para_lens[i2][j2] = ff2;

  fclose(fptr);

  return;
}

void varyzs_extend(int i, double pmin, double pmax, int n, int flag)
{
  int ff, k, nn;
  double p, hh, c2;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  fprintf(stderr, "######## compute chi^2 with changing extend zs\n");
  fprintf(stderr, " extend id = %d,  step = %d\n", i + 1, n);
  
  sprintf(fname, "%s_vary.dat", file_prefix);
  fprintf(stderr, " output file name = %s\n\n", fname);
    
  if((i >= num_ext) || (i < 0))
    terminator("extend id irrelevant (varyzs_extend)"); 

  if((n < 0) && (pmin <= 0.0))
    terminator("parameter must be positive for log spacing (varyzs_extend)"); 
  
  if(n < 0){
    nn = n * (-1);
  } else if(n > 0){
    nn = n;
  } else {
    nn = n;
    terminator("invalid step number (varyzs_extend)"); 
  }

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (varyzs_extend)");

  ff = flag_para_ext[i][0];
  flag_para_ext[i][0] = 0;

  if(n > 0){
    hh = (pmax - pmin) / ((double)nn);
  } else {
    hh = log10(pmax / pmin) / ((double)nn);
  }

  chi2_dumplimit = C2DUMPLIMIT_SET;

  for(k=0;k<=nn;k++){
    if(n > 0){
      p = pmin + ((double)k) * hh;
    } else {
      p = pmin * pow(10.0, ((double)k) * hh);
    }
    if(k == 0){ p = pmin; }
    if(k == nn){ p = pmax; }
    para_ext[i][0] = p;
    
    c2 = opt_lens(flag, -1);
    fprintf(stderr, "%e %e\n", p, c2);
    fprintf(fptr, "%e %e\n", p, c2);
    fflush(fptr);
  }
  
  fprintf(stderr, "\n");
  flag_para_ext[i][0] = ff;

  fclose(fptr);

  return;
}

void varyzs_point(int i, double pmin, double pmax, int n, int flag)
{
  int ff, k, nn;
  double p, hh, c2;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  fprintf(stderr, "######## compute chi^2 with changing point zs\n");
  fprintf(stderr, " point id = %d,  step = %d\n", i + 1, n);
  
  sprintf(fname, "%s_vary.dat", file_prefix);
  fprintf(stderr, " output file name = %s\n\n", fname);
    
  if((i >= num_poi) || (i < 0))
    terminator("point id irrelevant (varyzs_point)"); 

  if((n < 0) && (pmin <= 0.0))
    terminator("parameter must be positive for log spacing (varyzs_point)"); 
  
  if(n < 0){
    nn = n * (-1);
  } else if(n > 0){
    nn = n;
  } else {
    nn = n;
    terminator("invalid step number (varyzs_point)"); 
  }

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (varyzs_point)");

  ff = flag_para_poi[i][0];
  flag_para_poi[i][0] = 0;

  if(n > 0){
    hh = (pmax - pmin) / ((double)nn);
  } else {
    hh = log10(pmax / pmin) / ((double)nn);
  }

  chi2_dumplimit = C2DUMPLIMIT_SET;

  for(k=0;k<=nn;k++){
    if(n > 0){
      p = pmin + ((double)k) * hh;
    } else {
      p = pmin * pow(10.0, ((double)k) * hh);
    }
    if(k == 0){ p = pmin; }
    if(k == nn){ p = pmax; }
    para_poi[i][0] = p;
    
    c2 = opt_lens(flag, -1);
    fprintf(stderr, "%e %e\n", p, c2);
    fprintf(fptr, "%e %e\n", p, c2);
    fflush(fptr);
  }
  
  fprintf(stderr, "\n");
  flag_para_poi[i][0] = ff;

  fclose(fptr);

  return;
}

void varycosmo(char *keyword, double pmin, double pmax, int n, int flag)
{
  int ff, k, nn;
  double p, hh, c2;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  fprintf(stderr, "######## compute chi^2 with changing cosmological parameters\n");
  fprintf(stderr, " name = %s,  step = %d\n", keyword, n);
  
  if(!((strcmp(keyword, "omega") == 0) || (strcmp(keyword, "lambda") == 0) || (strcmp(keyword, "weos") == 0) || (strcmp(keyword, "hubble") == 0))){
    terminator("parameter name irrelevant (varycosmo)");
  }
 
  sprintf(fname, "%s_vary.dat", file_prefix);
  fprintf(stderr, " output file name = %s\n\n", fname);
    
  if((n < 0) && (pmin <= 0.0))
    terminator("parameter must be positive for log spacing (varycosmo)"); 
  
  if(n < 0){
    nn = n * (-1);
  } else if(n > 0){
    nn = n;
  } else {
    nn = n;
    terminator("invalid step number (varycosmo)"); 
  }

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (varycosmo)");

  ff = 0;
  if(strcmp(keyword, "omega") == 0){
    ff = ovary;
    ovary = 0;
  } else if(strcmp(keyword, "lambda") == 0){
    ff = lvary;
    lvary = 0;
  } else if(strcmp(keyword, "weos") == 0){
    ff = wvary;
    wvary = 0;
  } else if(strcmp(keyword, "hubble") == 0){
    ff = hvary;
    hvary = 0;
  }

  if(n > 0){
    hh = (pmax - pmin) / ((double)nn);
  } else {
    hh = log10(pmax / pmin) / ((double)nn);
  }

  chi2_dumplimit = C2DUMPLIMIT_SET;

  for(k=0;k<=nn;k++){
    if(n > 0){
      p = pmin + ((double)k) * hh;
    } else {
      p = pmin * pow(10.0, ((double)k) * hh);
    }
    if(k == 0){ p = pmin; }
    if(k == nn){ p = pmax; }

    if(strcmp(keyword, "omega") == 0){
      omega = p;
      if(flatfix != 0 ){ lambda = 1.0 - omega; }
    } else if(strcmp(keyword, "lambda") == 0){
      lambda = p;
      if(flatfix != 0 ){ omega = 1.0 - lambda; }
    } else if(strcmp(keyword, "weos") == 0){
      weos = p;
    } else if(strcmp(keyword, "hubble") == 0){
      hubble = p;
    }
    
    c2 = opt_lens(flag, -1);
    fprintf(stderr, "%e %e\n", p, c2);
    fprintf(fptr, "%e %e\n", p, c2);
    fflush(fptr);
  }
  
  fprintf(stderr, "\n");
  
  if(strcmp(keyword, "omega") == 0){
    ovary = ff;
  } else if(strcmp(keyword, "lambda") == 0){
    lvary = ff;
  } else if(strcmp(keyword, "weos") == 0){
    wvary = ff;
  } else if(strcmp(keyword, "hubble") == 0){
    hvary = ff;
  }

  fclose(fptr);

  return;
}

/*--------------------------------------------------------------
  repeat randomize+optimize to find global minimum
*/

void opt_explore(int nt, double c2lim, int flag)
{
  int i;
  double c2, c2min = C2DUMPLIMIT_SET;
  double bak_cosmo[NPAR_COSMO];
  double bak_para_lens[NMAX_LEN][NPAR_LEN];
  double bak_para_ext[NMAX_EXT][NPAR_EXT];
  double bak_para_poi[NMAX_POI][NPAR_POITAB];
  double bak_para_psf[NPAR_PSF];
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  fprintf(stderr, "######## explore chi^2\n");
  fprintf(stderr, " n_tot = %d,  chi2lim. = %e\n", nt, c2lim);
  
  sprintf(fname, "%s_explore.dat", file_prefix);
  fprintf(stderr, " output file name = %s\n\n", fname);
    
  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (opt_explore)");

  explore_parin(bak_cosmo, bak_para_lens, bak_para_ext, bak_para_poi, bak_para_psf);

  chi2_dumplimit = c2lim;
  for(i=0;i<nt;i++){
    c2 = opt_lens(flag, -1);
    fprintf(stderr, "run %5d: chi^2 = %e", i + 1, c2);
    fprintf(fptr, "run %5d: chi^2 = %e", i + 1, c2);
    if(c2 <= chi2_dumplimit){
      fprintf(stderr, "\n");
      fprintf(fptr, "\n");
      if(c2 < c2min){
	c2min = c2;
	explore_parin(bak_cosmo, bak_para_lens, bak_para_ext, bak_para_poi, bak_para_psf);
      }
    } else {
      fprintf(stderr, " [exceeds limit]\n");
      fprintf(fptr, " [exceeds limit]\n");
    }
    fflush(fptr);
    randomize(0);
  }

  fprintf(stderr, "\n");
  fprintf(fptr, "\n");
  fprintf(stderr, "best-fit chi^2 = %e\n", c2min);
  fprintf(fptr, "best-fit chi^2 = %e\n", c2min);

  explore_parout(bak_cosmo, bak_para_lens, bak_para_ext, bak_para_poi, bak_para_psf);
  dump_model(stderr);
  dump_model(fptr);

  fprintf(stderr, "\n");
  
  fclose(fptr);

  return;
}

void explore_parin(double bak_cosmo[NPAR_COSMO], double bak_para_lens[NMAX_LEN][NPAR_LEN], double bak_para_ext[NMAX_EXT][NPAR_EXT], double bak_para_poi[NMAX_POI][NPAR_POITAB], double bak_para_psf[NPAR_PSF])
{
  int j, k;

  bak_cosmo[0] = omega;
  bak_cosmo[1] = lambda;
  bak_cosmo[2] = weos;
  bak_cosmo[3] = hubble;
  for(j=0;j<num_len;j++){
    for(k=0;k<NPAR_LEN;k++) bak_para_lens[j][k] = para_lens[j][k];
  }
  for(j=0;j<num_ext;j++){
    for(k=0;k<NPAR_EXT;k++) bak_para_ext[j][k] = para_ext[j][k];
  }
  for(j=0;j<num_poi;j++){
    for(k=0;k<NPAR_POITAB;k++) bak_para_poi[j][k] = para_poi[j][k];
  }
  for(k=0;k<NPAR_PSF;k++) bak_para_psf[k] = para_psf[k];

  return;
}

void explore_parout(double bak_cosmo[NPAR_COSMO], double bak_para_lens[NMAX_LEN][NPAR_LEN], double bak_para_ext[NMAX_EXT][NPAR_EXT], double bak_para_poi[NMAX_POI][NPAR_POITAB], double bak_para_psf[NPAR_PSF])
{
  int j, k;

  omega = bak_cosmo[0];
  lambda = bak_cosmo[1];
  weos = bak_cosmo[2];
  hubble = bak_cosmo[3];
  for(j=0;j<num_len;j++){
    for(k=0;k<NPAR_LEN;k++) para_lens[j][k] = bak_para_lens[j][k];
  }
  for(j=0;j<num_ext;j++){
    for(k=0;k<NPAR_EXT;k++) para_ext[j][k] = bak_para_ext[j][k];
  }
  for(j=0;j<num_poi;j++){
    for(k=0;k<NPAR_POITAB;k++) para_poi[j][k] = bak_para_poi[j][k];
  }
  for(k=0;k<NPAR_PSF;k++) para_psf[k] = bak_para_psf[k];

  return;
}

/*--------------------------------------------------------------
  randomize lens parameters 
*/

void randomize(int verb)
{
  int i, j, nn;
  
  if(verb > 0) fprintf(stderr, "######## randomize parameters\n\n");

  for(i=0;i<num_len;i++){
    for(j=0;j<NPAR_LEN;j++){
      if(flag_para_lens[i][j] == 1){
	para_lens[i][j] = para_lens_min[i][j] + gsl_rng_uniform(ran_gsl) * (para_lens_max[i][j] - para_lens_min[i][j]);
      }
    }
  }

  if(ovary != 0){
    omega = omega_min + gsl_rng_uniform(ran_gsl) * (omega_max - omega_min);
    if(flatfix != 0 ){ lambda = 1.0 - omega; }
  }
  if(lvary != 0){
    lambda = lambda_min + gsl_rng_uniform(ran_gsl) * (lambda_max - lambda_min);
    if(flatfix != 0 ){ omega = 1.0 - lambda; }
  }
  if(wvary != 0){
    weos = weos_min + gsl_rng_uniform(ran_gsl) * (weos_max - weos_min);
  }
  if(hvary != 0){
    hubble = hubble_min + gsl_rng_uniform(ran_gsl) * (hubble_max - hubble_min);
  }

  nn = 0;
  for(i=0;i<num_ext;i++){
    for(j=0;j<NPAR_EXT;j++){
      if(flag_para_ext[i][j] == 1){
	nn++;
	para_ext[i][j] = para_ext_min[i][j] + gsl_rng_uniform(ran_gsl) * (para_ext_max[i][j] - para_ext_min[i][j]);
      }
    }
  }

  for(i=0;i<num_poi;i++){
    if(flag_para_poi[i][0] == 1){
      para_poi[i][0] = para_poi_min[i][0] + gsl_rng_uniform(ran_gsl) * (para_poi_max[i][0] - para_poi_min[i][0]);
    }
  }

  for(j=0;j<NPAR_PSF;j++){
    if((nn > 0) && (flag_seeing == 1) && (flag_para_psf[j] == 1)){
      para_psf[j] = para_psf_min[j] + gsl_rng_uniform(ran_gsl) * (para_psf_max[j] - para_psf_min[j]);
    }
  }

  return;
}
