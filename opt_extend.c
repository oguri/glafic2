#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

static int ndim;

/*--------------------------------------------------------------
  set def, optimize source 
*/

double chi2calc_opt_extend_out(void)
{
  double chi2min[NPAR_CHI2MIN];
  
  return chi2calc_opt_extend(chi2min, 1, 0);
}

double chi2calc_opt_extend(double chi2min[NPAR_CHI2MIN], int verb, int flag_reset)
{
  int i, j, k, nd, nfunc;

  double r, xx, yy;
  static double p[NDIMMAX + 2][NDIMMAX + 1];
  double y[NDIMMAX + 2];
  double par[NDIMMAX + 1], parori[NDIMMAX + 1], dp[NDIMMAX + 1];
  char fname[INPUT_MAXCHAR];

  if(flag_arrayobs == 0) terminator("you need to set obs file first");
  
  set_distance_lpl_init();

  ndim = 0;
  for(i=0;i<num_ext;i++){
    for(j=1;j<NPAR_EXT;j++){
      if(flag_para_ext[i][j] == 1) ndim++;
    }
  }
  if(skyfix == 0) ndim++;

  for(j=0;j<NPAR_PSF;j++){
      if((flag_seeing == 1) && (flag_para_psf[j] == 1)) ndim++;
  }

  if(ndim > NDIMMAX) terminator("too many source parameters");
  
  if(verb > 0){
    fprintf(stderr, "######## optimizing extended sources\n");
    fprintf(stderr, " number of parameters = %d\n", ndim);
  }

  flag_computeall = 0;
  i_ext_fid = -1;
  
  if(ndim == 0){ 
    r = chi2calc_ext_func(par);
    if(verb > 0){ 
      fprintf(stderr, "chi^2 = %e \n\n", r);
      dump_model(stderr);
      return r; 
    }
  }

  ndim = 0;
  for(i=0;i<num_ext;i++){

    ext_set_table_all(i);
    
    for(j=1;j<NPAR_EXT;j++){
      if(flag_para_ext[i][j] == 1){
	ndim++;
	par[ndim] = para_ext[i][j];

	switch(j){
	case 1:  /* norm */
	  if(obs_ext_prior[i] >= 0) par[ndim] = array_obs[obs_ext_prior[i]];
	  dp[ndim] = par[ndim] * amoeba_dp_mass;
	  break;
	  
	case 2:  /* x position */
	  if(obs_ext_prior[i] >= 0){ 
	    ktoxy_ext(obs_ext_prior[i], &xx, &yy);
	    par[ndim] = xx - array_ext_def[obs_ext_prior[i]] * dis_fac_ext[i];
	  }
	  dp[ndim] = pix_ext * amoeba_dp_xy;
	  break;

	case 3:  /* y position */
	  if(obs_ext_prior[i] >= 0){ 
	    ktoxy_ext(obs_ext_prior[i], &xx, &yy);
	    par[ndim] = yy - array_ext_def[obs_ext_prior[i] + nx_ext * ny_ext] * dis_fac_ext[i];
	  }
	  dp[ndim] = pix_ext * amoeba_dp_xy;
	  break;

	case 4:  /* e */
	  dp[ndim] = (1.0 - par[ndim]) * amoeba_dp_e;
	  break;

	case 5:  /* theta_e */
	  dp[ndim] = amoeba_dp_ang;
	  break;

	case 6:  /* size */
	  dp[ndim] = pix_ext * amoeba_dp_r; 
	  break;
	  
	case 7:  /* n */
	  dp[ndim] = amoeba_dp_n;
	  break;
	
	}

      }
    }
  }
  if(skyfix == 0){
    ndim++;
    par[ndim] = skymed;
    dp[ndim] = skysigma;
  }

  for(j=0;j<NPAR_PSF;j++){
    if((flag_seeing == 1) && (flag_para_psf[j] == 1)){
      ndim++;
      par[ndim] = para_psf[j];
	
      switch(j){
      case 0:  /* FWHM1 */
	dp[ndim] = amoeba_dp_psfw;
	break;
	  
      case 1:  /* e1 */
	dp[ndim] = amoeba_dp_psfe;
	break;
	
      case 2:  /* PA1 */
	dp[ndim] = amoeba_dp_psfpa;
	break;
	
      case 3:  /* beta1 */
	dp[ndim] = amoeba_dp_psfb;
	break;
	
      case 4:  /* FWHM2 */
	dp[ndim] = amoeba_dp_psfw;
	break;
	  
      case 5:  /* e2 */
	dp[ndim] = amoeba_dp_psfe;
	break;
	  
      case 6:  /* PA2 */
	dp[ndim] = amoeba_dp_psfpa;
	break;
	
      case 7:  /* beta2 */
	dp[ndim] = amoeba_dp_psfb;
	break;
	
      case 8:  /* frac */
	dp[ndim] = amoeba_dp_psff; 
	break;
	  
      }
    }
  }
  
  for(j=1;j<=ndim;j++){ parori[j] = par[j]; }
  
  /* input parameter set */
  for(i=1;i<=(ndim+1);i++){
    for(j=1;j<=ndim;j++){
      p[i][j] = par[j];
      if(i == j){ p[i][j] = p[i][j] + dp[i]; }
    }
  }

  for(i=1;i<=(ndim+1);i++){
    for(j=1;j<=ndim;j++){
      par[j] = p[i][j];
    }
    y[i] = chi2calc_ext_func(par);
  } 
  
  simplex(p, y, ndim, tol_amoeba, chi2calc_ext_func, &nfunc, nmax_amoeba, verb); 

  ndim = 0;
  for(i=0;i<num_ext;i++){
    for(j=1;j<NPAR_EXT;j++){
      if(flag_para_ext[i][j] == 1){
	ndim++;
	para_ext[i][j] = p[1][ndim]; 
      }
    }
  }
  if(skyfix == 0){
    ndim++;
     skymed = p[1][ndim]; 
  }

  for(j=0;j<NPAR_PSF;j++){
    if((flag_seeing == 1) && (flag_para_psf[j] == 1)){
      ndim++;
      para_psf[j] = p[1][ndim]; 
    }
  }

  r = chi2calc_extend(chi2min);

  if(verb > 0){ 
    nd = 0;
    for(i=0;i<(nx_ext*ny_ext);i++){
      if(array_ext_mask[i] == 0) nd++;
    }

    fprintf(stderr, "\n");
    dump_optextend(NULL, r, ndim, nd, nfunc, chi2min);
    sprintf(fname, "%s_optresult_extend.dat", file_prefix);
    fprintf(stderr, "\nwriting result to %s\n\n", fname);
    dump_optextend(fname, r, ndim, nd, nfunc, chi2min);
  } 

  if(flag_reset == 1) {
    k = 0;
    for(i=0;i<num_ext;i++){
      for(j=1;j<NPAR_EXT;j++){
	if(flag_para_ext[i][j] == 1){
	  k++;
	  para_ext[i][j] = parori[k];
	}
      }
    }
    if(skyfix == 0){
      k++;
      skymed = parori[k];
    }
    for(j=0;j<NPAR_PSF;j++){
      if((flag_seeing == 1) && (flag_para_psf[j] == 1)){
	k++;
	para_psf[j] = parori[k];
      }
    }
  }

  flag_computeall = 1;
  
  return r;
}

void dump_optextend(char *infile, double c2, int ndim, int nd, int nfunc, double chi2min[NPAR_CHI2MIN])
{
  static int flag;
  FILE* fptr;
  
  if(infile != NULL){
    if(flag != 1){
      fptr = fopen(infile, "w");
      flag = 1;
    } else {
      fptr = fopen(infile, "a");
    }
    if(fptr == NULL) terminator("failed at fopen (opt)");
  } else {
    fptr = stderr;
  }
  
  fprintf(fptr, "------------------------------------------\n");
  fprintf(fptr, "optextend ndim=%d\n", ndim);
  fprintf(fptr, "%d models calculated\n", nfunc);
  fprintf(fptr, "chi^2 = %e  [N_data: %d]\n", c2, nd);
  fprintf(fptr, " [image = %e; prior = %e]\n\n", chi2min[1], chi2min[2]);
  fprintf(fptr, "sky = %e  sigma(sky) = %e\n", skymed, skysigma);
  fprintf(fptr, "omega = %6.4f  lambda = %6.4f  weos = %7.4f  hubble = %6.4f\n", omega, lambda, weos, hubble);

  dump_model(fptr);

  fprintf(fptr, "------------------------------------------\n");

  if(infile != NULL) fclose(fptr);

  return;
}

/*--------------------------------------------------------------
  compute chi^2
*/

double chi2calc_ext_func(double par[NDIMMAX + 1])
{
  int i, j, k;
  double chi2min[NPAR_CHI2MIN];

  k = 0;
  for(i=0;i<num_ext;i++){
    for(j=1;j<NPAR_EXT;j++){
      if(flag_para_ext[i][j] == 1){
	k++;
	para_ext[i][j] = par[k];
      }
    }
  }
  if(skyfix == 0){
    k++;
    skymed = par[k];
  }
  for(j=0;j<NPAR_PSF;j++){
    if((flag_seeing == 1) && (flag_para_psf[j] == 1)){
      k++;
      para_psf[j] = par[k];
    }
  }

  return chi2calc_extend(chi2min);
}

double chi2calc_extend(double chi2min[NPAR_CHI2MIN])
{
  int i, k, nn;
  double fo, f, s;

  parmatch_ext();
  
  if(check_para_ext_all() > 0){ 
    chi2min[0] = chi2min[2] = chi2pen_range;
    chi2min[1] = 0.0;
    return chi2pen_range; 
  } 

  ext_set_image(0, 0, 0);
  nn = nx_ext * ny_ext;
    
  if(skyfix == 1) skymed = skyfix_value;

  chi2min[1] = 0.0;
  for(k=0;k<nn;k++){
    if(array_ext_mask[k] == 0){
      fo = array_obs[k];
      if(flag_obssig == 0){
	s = array_obs_noise[k];
      } else {
	s = array_obsnoise_file[k];
      }
      f = 0.0;
      for(i=0;i<num_ext;i++){
	f = f + array_ext_img[k + i * nn];
      }
      f = f + skymed;
      chi2min[1] = chi2min[1] + (f - fo) * (f - fo) / (s * s);
    }
  }
  

  chi2min[2] = chi2prior_ext();
  chi2min[0] = chi2min[1] + chi2min[2];

  return chi2min[0];

}

/*--------------------------------------------------------------
  check priors
*/

int check_para_ext(int i, int j)
{
  if((para_ext[i][j] < para_ext_min[i][j]) || (para_ext[i][j] > para_ext_max[i][j])){
    return 1;
  } else {
    return 0;
  }
    
}

int check_para_psf(int j)
{
  if((para_psf[j] < para_psf_min[j]) || (para_psf[j] > para_psf_max[j])){
    return 1;
  } else {
    return 0;
  }
    
}

int check_para_ext_all(void)
{
  int i, j;
  int r = 0;

  for(i=0;i<num_ext;i++){
    for(j=1;j<NPAR_EXT;j++){
      r = r + check_para_ext(i, j);
    }
  }
  
  if(flag_seeing == 1){
    for(j=0;j<NPAR_PSF;j++){
      r = r + check_para_psf(j);
    }
  }

  return r;
}

double chi2prior_ext(void)
{
  int i, j;
  double c2, lp, lm;
  
  c2 = 0.0;
  for(i=0;i<num_ext;i++){
    for(j=1;j<NPAR_EXT;j++){
      if(para_ext_sig[i][j] > 0.0){
	c2 = c2 + (para_ext[i][j] - para_ext_med[i][j]) * (para_ext[i][j] - para_ext_med[i][j]) / (para_ext_sig[i][j] * para_ext_sig[i][j]);
      } else if(para_ext_sig[i][j] < 0.0){
	if((para_ext[i][j] < 0.0) || (para_ext_med[i][j] < 0.0)) return chi2pen_range; 
	lp = log10(para_ext[i][j]);
	lm = log10(para_ext_med[i][j]);
	c2 = c2 + (lp - lm) * (lp - lm) / (para_ext_sig[i][j] * para_ext_sig[i][j]);
      }
      if((para_ext_rai[i][j] != i) || (para_ext_raj[i][j] != j)){
	if(para_ext_ras[i][j] > 0.0){
	  c2 = c2 + (para_ext[i][j] - para_ext_rat[i][j] * para_ext[para_ext_rai[i][j]][para_ext_raj[i][j]]) * (para_ext[i][j] - para_ext_rat[i][j] * para_ext[para_ext_rai[i][j]][para_ext_raj[i][j]]) / (para_ext_ras[i][j] * para_ext_ras[i][j]);
	} else if(para_ext_ras[i][j] < 0.0){
	  if((para_ext[i][j] < 0.0) || (para_ext_rat[i][j] < 0.0) || (para_ext[para_ext_rai[i][j]][para_ext_raj[i][j]] < 0.0)) return chi2pen_range; 
	  lp = log10(para_ext[i][j]);
	  lm = log10(para_ext_rat[i][j] * para_ext[para_ext_rai[i][j]][para_ext_raj[i][j]]);
	  c2 = c2 + (lp - lm) * (lp - lm) / (para_ext_ras[i][j] * para_ext_ras[i][j]);
	}
      }
      if((para_extlen_rai[i][j] >= 0) || (para_extlen_raj[i][j] >= 0)){
	if(para_extlen_ras[i][j] > 0.0){
	  c2 = c2 + (para_ext[i][j] - para_extlen_rat[i][j] * para_lens[para_extlen_rai[i][j]][para_extlen_raj[i][j]]) * (para_ext[i][j] - para_extlen_rat[i][j] * para_lens[para_extlen_rai[i][j]][para_extlen_raj[i][j]]) / (para_extlen_ras[i][j] * para_extlen_ras[i][j]);
	} else if(para_extlen_ras[i][j] < 0.0){
	  if((para_ext[i][j] < 0.0) || (para_extlen_rat[i][j] < 0.0) || (para_lens[para_extlen_rai[i][j]][para_extlen_raj[i][j]] < 0.0)) return chi2pen_range; 
	  lp = log10(para_ext[i][j]);
	  lm = log10(para_extlen_rat[i][j] * para_lens[para_extlen_rai[i][j]][para_extlen_raj[i][j]]);
	  c2 = c2 + (lp - lm) * (lp - lm) / (para_extlen_ras[i][j] * para_extlen_ras[i][j]);
	}
      }
    }
  }
   
  if(flag_seeing == 1){
    for(j=0;j<NPAR_PSF;j++){
      if(para_psf_sig[j] > 0.0){
	c2 = c2 + (para_psf[j] - para_psf_med[j]) * (para_psf[j] - para_psf_med[j]) / (para_psf_sig[j] * para_psf_sig[j]);
      } else if(para_psf_sig[j] < 0.0){
	if((para_psf[j] < 0.0) || (para_psf_med[j] < 0.0)) return chi2pen_range; 
	lp = log10(para_psf[j]);
	lm = log10(para_psf_med[j]);
	c2 = c2 + (lp - lm) * (lp - lm) / (para_psf_sig[j] * para_psf_sig[j]);
      }
      if(para_psf_raj[j] != j){
	if(para_psf_ras[j] > 0.0){
	  c2 = c2 + (para_psf[j] - para_psf_rat[j] * para_psf[para_psf_raj[j]]) * (para_psf[j] - para_psf_rat[j] * para_psf[para_psf_raj[j]]) / (para_psf_ras[j] * para_psf_ras[j]);
	} else if(para_psf_ras[j] < 0.0){
	  if((para_psf[j] < 0.0) || (para_psf_rat[j] < 0.0) || (para_psf[para_psf_raj[j]] < 0.0)) return chi2pen_range; 
	  lp = log10(para_psf[j]);
	  lm = log10(para_psf_rat[j] * para_psf[para_psf_raj[j]]);
	  c2 = c2 + (lp - lm) * (lp - lm) / (para_psf_ras[j] * para_psf_ras[j]);
	}
      }
    }
  }
 
  return c2;
}

void parmatch_ext(void)
{
  int i, j;
  
  for(i=0;i<num_ext;i++){
    for(j=1;j<NPAR_EXT;j++){
      if((para_ext_rai[i][j] != i) || (para_ext_raj[i][j] != j)){
	if(para_ext_ras[i][j] == 0.0){
	  para_ext[i][j] = para_ext_rat[i][j] * para_ext[para_ext_rai[i][j]][para_ext_raj[i][j]];
	}
      }
      if((para_extlen_rai[i][j] >= 0) && (para_extlen_raj[i][j] >= 0)){
	if(para_extlen_ras[i][j] == 0.0){
	  para_ext[i][j] = para_extlen_rat[i][j] * para_lens[para_extlen_rai[i][j]][para_extlen_raj[i][j]];
	}
      }
    }
  }
  
  for(j=0;j<NPAR_PSF;j++){
    if(para_psf_raj[j] != j){
      if(para_psf_ras[j] == 0.0){
	para_psf[j] = para_psf_rat[j] * para_psf[para_psf_raj[j]];
      }
    }
  }
  
  return ;
}
