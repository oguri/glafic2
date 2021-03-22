#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

static int ndim, nl, np, ne;
static int flag_sav;

/*--------------------------------------------------------------
  optimize lens new
*/

double opt_lens(int flag, int verb)
{
  /* flag - 0: point+extend,  1: extend,  2: point */
  int i, j, nd, res, nfunc;
  static double p[NDIMMAX + 2][NDIMMAX + 1]; 
  double c2, c2old, del, dmax, xx, yy;
  double y[NDIMMAX + 2]; 
  double par[NDIMMAX + 1], dp[NDIMMAX + 1];
  double chi2min_point[NMAX_POI][NPAR_CHI2];
  double chi2min_extend[NPAR_CHI2MIN];
  char fname[INPUT_MAXCHAR];

  opt_lens_static(flag);

  for(i=0;i<NMAX_POI;i++){
    for(j=0;j<NPAR_CHI2;j++) chi2min_point[i][j] = 0.0;
  }
  for(j=0;j<NPAR_CHI2MIN;j++) chi2min_extend[j] = 0.0;

  if(verb > 0){
    fprintf(stderr, "######## optimizing lens\n");
    fprintf(stderr, " number of parameters = %d (lens) ", nl);
    if(flag != 2) fprintf(stderr, "+ %d (extend) ", ne);
    if(flag != 1) fprintf(stderr, "+ %d (point) ", np);
    fprintf(stderr, "\n\n");
  }
  sprintf(fname, "%s_optresult.dat", file_prefix);

  res = -1;
  c2 = 0.0;
  do{
    res++;

    if(np > 0){ 
      dmax = 0.0;
      for(i=0;i<num_poi;i++){
	del = calcdelta_i(i);
	if(del > dmax) dmax = del;
      }
      del = dmax; 
    } else { 
      del = amoeba_delmax; 
    }
    
    if(verb > 0) printf("amoeba_delta = %e\n\n", del);
    
    if(ne > 0){
      flag_computeall = 0;
      i_ext_fid = -1;
    }

    ndim = 0;
    for(i=0;i<num_len;i++){
      for(j=0;j<NPAR_LEN;j++){
	if(flag_para_lens[i][j] == 1){
	  ndim++;
	  par[ndim] = para_lens[i][j];
	  
	  switch(j){
	  case 0:  /* zl */
	    dp[ndim] = par[ndim] * amoeba_dp_z * del;
	    break;
	    
	  case 1:  /* mass */
	    dp[ndim] = par[ndim] * amoeba_dp_mass * del;
	    break;
	    
	  case 2:  /* x position */
	    dp[ndim] = pix_ext * amoeba_dp_xy * del;
	    break;
	    
	  case 3:  /* y position */
	    dp[ndim] = pix_ext * amoeba_dp_xy * del;
	    break;
	    
	  case 4:  /* e */
	    dp[ndim] = (1.0 - par[ndim]) * amoeba_dp_e * del;
	    break;
	    
	  case 5:  /* theta_e */
	    dp[ndim] = amoeba_dp_ang * del;
	    break;
	    
	  case 6:  /* r0 */
	    dp[ndim] = pix_ext * amoeba_dp_r * del; 
	    break;
	    
	  case 7:  /* n */
	    dp[ndim] = amoeba_dp_n * del;
	    break;
	    
	  }
	  
	}
      }
    }
    if(ovary != 0){ ndim++; par[ndim] = omega;  dp[ndim] = amoeba_dp_cosmo * del; }
    if(lvary != 0){ ndim++; par[ndim] = lambda; dp[ndim] = amoeba_dp_cosmo * del; }
    if(wvary != 0){ ndim++; par[ndim] = weos;   dp[ndim] = amoeba_dp_cosmo * del; }
    if(hvary != 0){ ndim++; par[ndim] = hubble; dp[ndim] = amoeba_dp_cosmo * del; }
    
    for(i=0;i<num_ext;i++){
      if((res == 0) && (obs_ext_prior[i] >= 0)) ext_set_table_all(i);
      for(j=0;j<NPAR_EXT;j++){
	if(flag_para_ext[i][j] == 1){
	  ndim++;
	  par[ndim] = para_ext[i][j];

	  switch(j){
	  case 0:  /* zs */
	    dp[ndim] = amoeba_dp_z;
	    break;
	    
	  case 1:  /* norm */
	    if((res == 0) && (obs_ext_prior[i] >= 0)) par[ndim] = array_obs[obs_ext_prior[i]];
	    dp[ndim] = par[ndim] * amoeba_dp_mass;
	    break;
	    
	  case 2:  /* x position */
	    if((res == 0) && (obs_ext_prior[i] >= 0)){ 
	      ktoxy_ext(obs_ext_prior[i], &xx, &yy);
	      par[ndim] = xx - array_ext_def[obs_ext_prior[i]] * dis_fac_ext[i];
	    }
	    dp[ndim] = pix_ext * amoeba_dp_xy;
	    break;
	    
	  case 3:  /* y position */
	  if((res == 0) && (obs_ext_prior[i] >= 0)){ 
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
    if((skyfix == 0) && (ne > 0)){
      ndim++;
      par[ndim] = skymed;
      dp[ndim] = skysigma;
    }
   
    for(i=0;i<num_poi;i++){
      if(flag_para_poi[i][0] == 1){
	ndim++;
	par[ndim] = para_poi[i][0];
	dp[ndim] = amoeba_dp_z * del;  /* zs */
      }
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
      y[i] = chi2calc(par);
    } 
    if(ndim > 0) simplex(p, y, ndim, tol_amoeba_lens, chi2calc, &nfunc, nmax_amoeba, verb); 
 
    for(j=1;j<=ndim;j++){ par[j] = p[1][j]; } 

    partopara(par);
    parmatch_lens();
    parmatch_ext();
    set_distance_lpl_init();

    c2old = c2;
    c2 = chi2tot(chi2min_point, chi2min_extend);
    
    poi_unset_table();

    if(verb > 0){
      nd = 0;
      if(ne > 0){
	for(i=0;i<(nx_ext*ny_ext);i++){
	  if(array_ext_mask[i] == 0) nd++;
	}
      }
      
      fprintf(stderr, "\n");
      dump_opt(NULL, c2, nl + np + ne, res + 1, nfunc, nd, chi2min_point, chi2min_extend);
      fprintf(stderr, "\nwriting result to %s\n\n", fname);
      dump_opt(fname, c2, nl + np + ne, res + 1, nfunc, nd, chi2min_point, chi2min_extend);
    }

  }while((nfunc > 0) && ((res < chi2_restart) || ((chi2_restart < 0) && (fabs((c2 - c2old) / c2) > (0.5 * tol_amoeba_lens)))) && (res < chi2_restart_max));

  if((verb < 0) && (c2 <= chi2_dumplimit)){
    nd = 0;
    if(ne > 0){
      for(i=0;i<(nx_ext*ny_ext);i++){
	if(array_ext_mask[i] == 0) nd++;
      }
    }
    dump_opt(fname, c2, nl + np + ne, res + 1, nfunc, nd, chi2min_point, chi2min_extend);
  }

  return c2;
}

/*--------------------------------------------------------------
  compute chi^2
*/

double chi2calc(double par[NDIMMAX + 1])
{
  double chi2min_point[NMAX_POI][NPAR_CHI2];
  double chi2min_extend[NPAR_CHI2MIN];

  partopara(par);
  
  parmatch_lens();
  parmatch_ext();

  set_distance_lpl_init();
  
  return chi2tot(chi2min_point, chi2min_extend);
}

double chi2tot(double chi2min_point[NMAX_POI][NPAR_CHI2],  double chi2min_extend[NPAR_CHI2MIN])
{
  double r;

  if(check_para_lens_all() > 0) return chi2pen_range; 
  
  r = 0.0;
  
  if((flag_sav != 2) && (ne > 0)){
    if(check_para_ext_all() > 0) return chi2pen_range; 
    flag_computeall = 0;
    i_ext_fid = -1;
    r = r + chi2calc_extend(chi2min_extend);
  }

  if((flag_sav != 1) && (np > 0)) r = r + chi2calc_opt_point(chi2min_point, 0);
  
  r = r + chi2prior_lens() + chi2prior_map();
  
  return r;
}

/* {lens} -> {cosmo} -> {ext} -> sky -> {point} -> {psf} */

void partopara(double par[NDIMMAX + 1])
{
  int i, j, k, nn;
  
  k = 0;
  for(i=0;i<num_len;i++){
    for(j=0;j<NPAR_LEN;j++){
      if(flag_para_lens[i][j] == 1){
	k++;
	para_lens[i][j] = par[k];
      }
    }
  }
  if(ovary != 0){ k++; omega = par[k];  }
  if(lvary != 0){ k++; lambda = par[k]; }
  if(wvary != 0){ k++; weos = par[k];   }
  if(hvary != 0){ k++; hubble = par[k]; }
  
  nn = 0;
  for(i=0;i<num_ext;i++){
    for(j=0;j<NPAR_EXT;j++){
      if(flag_para_ext[i][j] == 1){
	k++;
	nn++;
	para_ext[i][j] = par[k];
      }
    }
  }
  if((skyfix == 0) && (nn > 0)){
    k++;
    skymed = par[k];
  }
  
  for(i=0;i<num_poi;i++){
    if(flag_para_poi[i][0] == 1){
      k++;
      para_poi[i][0] = par[k];
    }
  }

  for(j=0;j<NPAR_PSF;j++){
    if((nn > 0) && (flag_seeing == 1) && (flag_para_psf[j] == 1)){
      k++;
      para_psf[j] = par[k];
    }
  }

  return;
}

int opt_lens_calcndim(void)
{
  int i, j, n, nn;

  n = 0;
  nn = 0;
  nl = 0;
  ne = 0;
  np = 0;
  for(i=0;i<num_len;i++){
    for(j=0;j<NPAR_LEN;j++){
      if(flag_para_lens[i][j] == 1){ n++; nl++; }
    }
  }
  if(ovary != 0){ n++; nl++; }
  if(lvary != 0){ n++; nl++; }
  if(wvary != 0){ n++; nl++; }
  if(hvary != 0){ n++; nl++; }

  nn = 0;
  for(i=0;i<num_ext;i++){
    for(j=0;j<NPAR_EXT;j++){
      if(flag_para_ext[i][j] == 1){ n++; nn++; ne++; }
    }
  }
  if((skyfix == 0) && (nn > 0)){
    n++; ne++;
  }
  
  for(i=0;i<num_poi;i++){
    if(flag_para_poi[i][0] == 1){ n++; np++; }
    for(j=1;j<NPAR_POI;j++){
      if(flag_para_poi[i][j] == 1){ np++; }
    }
  }
  
  if((nn > 0) && (flag_seeing == 1)){
    for(j=0;j<NPAR_PSF;j++){
      if(flag_para_psf[j] == 1){ n++; ne++; }
    }
  }

  if((ne > 0) && (flag_arrayobs == 0)) terminator("you need to set obs file first");

  return n;
}

void opt_lens_static(int flag)
{
  flag_sav = flag;
  ndim = opt_lens_calcndim();

  if(ndim > NDIMMAX) terminator("too many lens parameters");
  
  if((ne > NDIMMAX) && (flag != 2)) terminator("too many extend parameters");
  
  if((np > NDIMMAX) && (flag != 1)) terminator("too many point parameters");
  
  return;
}

void dump_opt(char *infile, double c2, int ndim, int res, int nfunc, int nd, double chi2min_point[NMAX_POI][NPAR_CHI2], double chi2min_extend[NPAR_CHI2MIN])
{
  static int flag;
  int i, j;
  FILE* fptr;
  
  if(infile != NULL){
    if(flag != 1){
      fptr = fopen(infile, "w");
      flag = 1;
    } else {
      fptr = fopen(infile, "a");
    }
    if(fptr == NULL) terminator("failed at fopen (optimize)");
  } else {
    fptr = stderr;
  }
  
  fprintf(fptr, "------------------------------------------\n");
  fprintf(fptr, "optimize ndim=%d\n", ndim);
  fprintf(fptr, "run %d: %d lens models calculated\n", res, nfunc);
  fprintf(fptr, "chi^2 = %e  [N_data(extend): %d]\n", c2, nd);
  fprintf(fptr, " extend     : %e %e %e\n", chi2min_extend[0], chi2min_extend[1], chi2min_extend[2]);
  for(i=0;i<num_poi;i++){
    fprintf(fptr, " point no %-2d:", i + 1);
    for(j=0;j<NPAR_CHI2;j++) fprintf(fptr, " %e", chi2min_point[i][j]);
    fprintf(fptr, "\n");
  }
  fprintf(fptr, " lens prior : %e\n", chi2prior_lens());
  fprintf(fptr, " map prior  : %e\n\n", chi2prior_map());
  if(ne > 0){
    fprintf(fptr, "sky = %e  sigma(sky) = %e\n", skymed, skysigma);
  }
  fprintf(fptr, "omega = %6.4f  lambda = %6.4f  weos = %7.4f  hubble = %6.4f\n", omega, lambda, weos, hubble);
  
  dump_model(fptr);

  fprintf(fptr, "------------------------------------------\n");

  if(infile != NULL) fclose(fptr);

  return;
}

/*--------------------------------------------------------------
  check parameter range, 0 = ok, 1 = par out of range
*/

int check_para_cosmo(void)
{
  int r;

  r = 0;

  if((omega < omega_min) || (omega > omega_max)) r++;
  if((lambda < lambda_min) || (lambda > lambda_max)) r++;
  if((weos < weos_min) || (weos > weos_max)) r++;
  if((hubble < hubble_min) || (hubble > hubble_max)) r++;

  return r;
}

int check_para_lens(int i, int j)
{
  if((para_lens[i][j] < para_lens_min[i][j]) || (para_lens[i][j] > para_lens_max[i][j])){
    return 1;
  } else {
    return 0;
  }
    
}

int check_para_lens_all(void)
{
  int i, j;
  int r = 0;

  for(i=0;i<num_len;i++){
    for(j=0;j<NPAR_LEN;j++){
      r = r + check_para_lens(i, j);
    }
  }
  
  for(i=0;i<num_ext;i++) r = r + check_para_ext(i, 0);  
  for(i=0;i<num_poi;i++) r = r + check_para_poi_zs(i);

  r = r + check_para_cosmo();

  return r;
}

double chi2prior_lens(void)
{
  int i, j;
  double c2, lp, lm;
  
  if(check_para_lens_all() > 0) return chi2pen_range; 

  c2 = 0.0;
  for(i=0;i<num_len;i++){
    for(j=0;j<NPAR_LEN;j++){
      if(para_lens_sig[i][j] > 0.0){
	c2 = c2 + (para_lens[i][j] - para_lens_med[i][j]) * (para_lens[i][j] - para_lens_med[i][j]) / (para_lens_sig[i][j] * para_lens_sig[i][j]);
      } else if(para_lens_sig[i][j] < 0.0){
	if((para_lens[i][j] < 0.0) || (para_lens_med[i][j] < 0.0)) return chi2pen_range; 
	lp = log10(para_lens[i][j]);
	lm = log10(para_lens_med[i][j]);
	c2 = c2 + (lp - lm) * (lp - lm) / (para_lens_sig[i][j] * para_lens_sig[i][j]);
      }
      if((para_lens_rai[i][j] != i) || (para_lens_raj[i][j] != j)){
	if(para_lens_ras[i][j] > 0.0){
	  c2 = c2 + (para_lens[i][j] - para_lens_rat[i][j] * para_lens[para_lens_rai[i][j]][para_lens_raj[i][j]]) * (para_lens[i][j] - para_lens_rat[i][j] * para_lens[para_lens_rai[i][j]][para_lens_raj[i][j]]) / (para_lens_ras[i][j] * para_lens_ras[i][j]);
	} else if(para_lens_ras[i][j] < 0.0){
	  if((para_lens[i][j] < 0.0) || (para_lens_rat[i][j] < 0.0) || (para_lens[para_lens_rai[i][j]][para_lens_raj[i][j]] < 0.0)) return chi2pen_range; 
	  lp = log10(para_lens[i][j]);
	  lm = log10(para_lens_rat[i][j] * para_lens[para_lens_rai[i][j]][para_lens_raj[i][j]]);
	  c2 = c2 + (lp - lm) * (lp - lm) / (para_lens_ras[i][j] * para_lens_ras[i][j]);
	}
      }
    }
  }
  
  if(oerror > 0.0){ c2 = c2 + (omega - omedian) * (omega - omedian) / (oerror * oerror); 
  } else if(oerror < 0.0){
    if((omega < 0.0) || (omedian < 0.0)) return chi2pen_range; 
    lp = log10(omega); lm = log10(omedian);
    c2 = c2 + (lp - lm) * (lp - lm) / (oerror * oerror); 
  }
  
  if(lerror > 0.0){ c2 = c2 + (lambda - lmedian) * (lambda - lmedian) / (lerror * lerror); 
  } else if(lerror < 0.0){
    if((lambda < 0.0) || (lmedian < 0.0)) return chi2pen_range; 
    lp = log10(lambda); lm = log10(lmedian);
    c2 = c2 + (lp - lm) * (lp - lm) / (lerror * lerror); 
  }
  
  if(werror > 0.0){ c2 = c2 + (weos - wmedian) * (weos - wmedian) / (werror * werror); 
  } else if(werror < 0.0){
    if((weos < 0.0) || (wmedian < 0.0)) return chi2pen_range; 
    lp = log10(weos); lm = log10(wmedian);
    c2 = c2 + (lp - lm) * (lp - lm) / (werror * werror); 
  }
  
  if(herror > 0.0){ c2 = c2 + (hubble - hmedian) * (hubble - hmedian) / (herror * herror); 
  } else if(herror < 0.0){
    if((hubble < 0.0) || (hmedian < 0.0)) return chi2pen_range; 
    lp = log10(hubble); lm = log10(hmedian);
    c2 = c2 + (lp - lm) * (lp - lm) / (herror * herror); 
  }
  
  for(i=0;i < num_ext;i++){
    if(para_ext_sig[i][0] > 0.0){
      c2 = c2 + (para_ext[i][0] - para_ext_med[i][0]) * (para_ext[i][0] - para_ext_med[i][0]) / (para_ext_sig[i][0] * para_ext_sig[i][0]);
    } else if(para_ext_sig[i][0] < 0.0){
      if((para_ext[i][0] < 0.0) || (para_ext_med[i][0] < 0.0)) return chi2pen_range; 
      lp = log10(para_ext[i][0]);
      lm = log10(para_ext_med[i][0]);
      c2 = c2 + (lp - lm) * (lp - lm) / (para_ext_sig[i][0] * para_ext_sig[i][0]);
    }
    if(para_ext_rai[i][0] != i){
      if(para_ext_ras[i][0] > 0.0){
	c2 = c2 + (para_ext[i][0] - para_ext_rat[i][0] * para_ext[para_ext_rai[i][0]][0]) * (para_ext[i][0] - para_ext_rat[i][0] * para_ext[para_ext_rai[i][0]][0]) / (para_ext_ras[i][0] * para_ext_ras[i][0]);
      } else if(para_ext_ras[i][0] < 0.0){
	if((para_ext[i][0] < 0.0) || (para_ext_rat[i][0] < 0.0) || (para_ext[para_ext_rai[i][0]][0] < 0.0)) return chi2pen_range; 
	lp = log10(para_ext[i][0]);
	lm = log10(para_ext_rat[i][0] * para_ext[para_ext_rai[i][0]][0]);
	c2 = c2 + (lp - lm) * (lp - lm) / (para_ext_ras[i][0] * para_ext_ras[i][0]);
      }
    }
  }

  for(i=0;i<num_poi;i++){
    if(para_poi_sig[i][0] > 0.0){
      c2 = c2 + (para_poi[i][0] - para_poi_med[i][0]) * (para_poi[i][0] - para_poi_med[i][0]) / (para_poi_sig[i][0] * para_poi_sig[i][0]);
    } else if(para_poi_sig[i][0] < 0.0){
      if((para_poi[i][0] < 0.0) || (para_poi_med[i][0] < 0.0)) return chi2pen_range; 
      lp = log10(para_poi[i][0]);
      lm = log10(para_poi_med[i][0]);
      c2 = c2 + (lp - lm) * (lp - lm) / (para_poi_sig[i][0] * para_poi_sig[i][0]);
    }
    if(para_poi_rai[i][0] != i){
      if(para_poi_ras[i][0] > 0.0){
	c2 = c2 + (para_poi[i][0] - para_poi_rat[i][0] * para_poi[para_poi_rai[i][0]][0]) * (para_poi[i][0] - para_poi_rat[i][0] * para_poi[para_poi_rai[i][0]][0]) / (para_poi_ras[i][0] * para_poi_ras[i][0]);
      } else if(para_poi_ras[i][0] < 0.0){
	if((para_poi[i][0] < 0.0) || (para_poi_rat[i][0] < 0.0) || (para_poi[para_poi_rai[i][0]][0] < 0.0)) return chi2pen_range; 
	lp = log10(para_poi[i][0]);
	lm = log10(para_poi_rat[i][0] * para_poi[para_poi_rai[i][0]][0]);
	c2 = c2 + (lp - lm) * (lp - lm) / (para_poi_ras[i][0] * para_poi_ras[i][0]);
      }
    }
  }

  return c2;
}

void parmatch_lens(void)
{
  int i, j;
  
  for(i=0;i<num_len;i++){
    for(j=0;j<NPAR_LEN;j++){
      if((para_lens_rai[i][j] != i) || (para_lens_raj[i][j] != j)){
	if(para_lens_ras[i][j] == 0.0){
	  para_lens[i][j] = para_lens_rat[i][j] * para_lens[para_lens_rai[i][j]][para_lens_raj[i][j]];
	}
      }
    }
  }

  for(i=0;i<num_poi;i++){
    if(para_poi_rai[i][0] != i){
      if(para_poi_ras[i][0] == 0.0){
	para_poi[i][0] = para_poi_rat[i][0] * para_poi[para_poi_rai[i][0]][0];
      }
    }
    if((para_poilen_rai[i][0] >= 0) && (para_poilen_raj[i][0] >= 0)){
      if(para_poilen_ras[i][0] == 0.0){
	para_poi[i][0] = para_poilen_rat[i][0] * para_lens[para_poilen_rai[i][0]][para_poilen_raj[i][0]];
      }
    }
  }
  
  for(i=0;i<num_ext;i++){
    if(para_ext_rai[i][0] != i){
      if(para_ext_ras[i][0] == 0.0){
	para_ext[i][0] = para_ext_rat[i][0] * para_ext[para_ext_rai[i][0]][0];
      }
    }
    if((para_extlen_rai[i][0] >= 0) && (para_extlen_raj[i][0] >= 0)){
      if(para_extlen_ras[i][0] == 0.0){
	para_ext[i][0] = para_extlen_rat[i][0] * para_lens[para_extlen_rai[i][0]][para_extlen_raj[i][0]];
      }
    }
  }
  
  return ;
}

double chi2prior_map(void)
{
  int i;
  double pout[NPAR_LMODEL];
  double c2, f, par, sig, lp, lm;
  static double zszs;

  zszs = 0.0;
  c2 = 0.0;
  for(i=0;i<num_mapprior;i++){
    if(zszs != para_mapprior[i][0]){
      set_distance_lpl_zs(para_mapprior[i][0]);
      zszs = para_mapprior[i][0];
    }
    /* alponly == 0 -> no potential */
    lensmodel(para_mapprior[i][1], para_mapprior[i][2], pout, 0, 0);
    par = para_mapprior[i][3];
    sig = para_mapprior[i][4];
    
    f = 0.0;
    switch(flag_para_mapprior[i]){
    case 1:
      f = pout[3];                    /* kappa */
      break;
      
    case 2:
      f = 1.0 / fabs(pout[6] + imag_ceil);  /* mag */
      break;
      
    case 3:
      f = pout[4];                    /* gamma1 */
      break;

    case 4:
      f = pout[5];                    /* gamma2 */
      break;
      
    case 5:
      f = pout[4] / (1.0 - pout[3]);         /* g1 */
      break;

    case 6:
      f = pout[5] / (1.0 - pout[3]);         /* g2 */
      break;
    }

    if(sig > 0.0){
      c2 = c2 + (f - par) * (f - par) / (sig * sig);
    } else if(sig < 0.0){
      if((f < 0.0) || (par < 0.0)) return chi2pen_range; 
      lp = log10(f);
      lm = log10(par);
      c2 = c2 + (lp - lm) * (lp - lm) / (sig * sig);
    }
  }

  return c2;
}

