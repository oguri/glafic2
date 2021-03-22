#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

static int ntot, ntot2, ndim[NMAX_POI], nt[NMAX_POI];
static int i_sav;
static int flag_scalc;
static double umod[2];

/*--------------------------------------------------------------
  optimizing point sources
*/

double chi2calc_opt_point_out(void)
{
  double c2min[NMAX_POI][NPAR_CHI2];

  return chi2calc_opt_point(c2min, 1);
}

double chi2calc_opt_point(double c2min[NMAX_POI][NPAR_CHI2], int verb)
{
  int i;
  double r;
  char fname[INPUT_MAXCHAR];

  if(flag_pointobs == 0) terminator("you need to set obs file first");

  set_distance_lpl_init();
  poi_unset_table();

  if(chi2_point_splane == 0){
    r = chi2calc_opt_iplane(c2min, verb);
  } else {
    r = chi2calc_opt_splane(c2min, verb);
  }

  if(verb > 0){ 
    ntot = 0;
    ntot2 = 0;
    for(i=0;i<num_poi;i++){
      ntot = ntot + ndim[i];
      ntot2 = ntot2 + nt[i];
    }

    fprintf(stderr, "\n");
    dump_optpoint(NULL, r, ntot, ntot2, c2min);
    sprintf(fname, "%s_optresult_point.dat", file_prefix);
    fprintf(stderr, "\nwriting result to %s\n\n", fname);
    dump_optpoint(fname, r, ntot, ntot2, c2min);
  }

  return r;
}

void dump_optpoint(char *infile, double c2, int ntot, int nfunc, double c2min[NMAX_POI][NPAR_CHI2])
{
  static int flag;
  int i;
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
  fprintf(fptr, "optpoint ndim=%d\n", ntot);
  fprintf(fptr, "%d lens models calculated\n", nfunc);
  fprintf(fptr, "chi^2 = %e \n\n", c2);
  for(i=0;i<num_poi;i++){
    fprintf(fptr, " point id = %d, number of parameters = %d\n", i + 1, ndim[i]);
    fprintf(fptr, " %d models calculated\n", nt[i]);
    fprintf(fptr, " flux = %e, t0 = %e\n", para_poi[i][3], para_poi[i][4]);
    fprintf(fptr, " chi^2 = %e\n", c2min[i][0]);
    fprintf(fptr, "  [pos = %e; flux = %e; td = %e; prior = %e]\n\n", c2min[i][1], c2min[i][2], c2min[i][3], c2min[i][4]);
  }
  fprintf(fptr, "omega = %6.4f  lambda = %6.4f  weos = %7.4f  hubble = %6.4f\n", omega, lambda, weos, hubble);
    
  dump_model(fptr);

  fprintf(fptr, "------------------------------------------\n");

  if(infile != NULL) fclose(fptr);

  return;
}

/*--------------------------------------------------------------
  lens plane chi2
*/

double chi2calc_opt_iplane(double c2min[NMAX_POI][NPAR_CHI2], int verb)
{
  int i, k, l, f, nfunc;
  double t1, t2;
  static double p[NDIMMAX + 2][NDIMMAX + 1]; 
  double xs[2 * NMAX_POIMG], ys[2 * NMAX_POIMG], c2m, c2tot;
  double c2[NPAR_CHI2];
  double pout[NPAR_LMODEL];
  double y[NDIMMAX + 2];
  double parf[NDIMMAX + 1];

  if(verb > 0){ 
    fprintf(stderr, "######## optimizing point sources (image plane)\n\n");
  }

  c2tot = 0.0;
  t1 = 0.0;
  t2 = 0.0;

  for(i=0;i<num_poi;i++){
    ndim[i] = flag_para_poi[i][1] + flag_para_poi[i][2];
    i_sav = i;
    nt[i] = 0;
    set_distance_lpl_zs(para_poi[i][0]);

    f = parmatch_poi(i);
    
    if((ndim[i] > 0) && (obs_numimg[i] >= 2) && (f == 0)){
      c2m = CHI2_MIN_SET;
      for(k=0;k<obs_numimg[i];k++){
	lensmodel(tab_obs[i][k][0], tab_obs[i][k][1], pout, 1, 0);
	xs[k] = tab_obs[i][k][0] - pout[0];
	ys[k] = tab_obs[i][k][1] - pout[1];
	if(k > 0){
	  xs[obs_numimg[i] + k - 1] = 0.5 * (xs[k - 1] + xs[k]);
	  ys[obs_numimg[i] + k - 1] = 0.5 * (ys[k - 1] + ys[k]);
	} 
      }
      
      for(k=0;k<=((obs_numimg[i]-1)/(ndim[i]+1));k++){ 
	if(ndim[i] == 2){
	  p[1][1] = xs[k * 3];
	  p[1][2] = ys[k * 3];
	  p[2][1] = xs[k * 3 + 1];
	  p[2][2] = ys[k * 3 + 1];
	  p[3][1] = xs[k * 3 + 2];
	  p[3][2] = ys[k * 3 + 2];
	  for(l=1;l<=3;l++) y[l] = chi2calc_point_iplane(i, p[l][1], p[l][2], c2); 
	}
	
	if((ndim[i] == 1) && (flag_para_poi[i][1] == 1)){
	  p[1][1] = xs[k * 2];
	  p[2][1] = xs[k * 2 + 1];
	  for(l=1;l<=2;l++) y[l] = chi2calc_point_iplane(i, p[l][1], para_poi[i][2], c2);  
	}
	
	if((ndim[i] == 1) && (flag_para_poi[i][2] == 1)){
	  p[1][1] = ys[k * 2];
	  p[2][1] = ys[k * 2 + 1];
	  for(l=1;l<=2;l++) y[l] = chi2calc_point_iplane(i, para_poi[i][1], p[l][1], c2);  
	}
	
	simplex(p, y, ndim[i], tol_amoeba, chi2calc_opt_func, &nfunc, nmax_amoeba_point, verb);
	
	nt[i] = nt[i] + nfunc;

       if(y[1]<c2m){ 
          c2m=y[1]; 
          for(l=1;l<=ndim[i];l++) parf[l]=p[1][l];
        }
      }
      
      if(flag_para_poi[i][1] == 1){
	para_poi[i][1] = parf[1];
	if(flag_para_poi[i][2] == 1){
	  para_poi[i][2] = parf[2];
	} 
      } else {
	if(flag_para_poi[i][2] == 1){
	  para_poi[i][2] = parf[1];
	} 
      }
      
    } else if((obs_numimg[i] == 1) && (f == 0)){
      lensmodel(tab_obs[i][0][0], tab_obs[i][0][1], pout, 1, 0);
      if(flag_para_poi[i][1] == 1){ para_poi[i][1] = tab_obs[i][0][0] - pout[0]; }
      if(flag_para_poi[i][2] == 1){ para_poi[i][2] = tab_obs[i][0][1] - pout[1]; }
    } 

    chi2calc_point_iplane(i, para_poi[i][1], para_poi[i][2], c2);
    nt[i]++;

    for(l=0;l<NPAR_CHI2;l++) c2min[i][l] = c2[l];
    c2tot = c2tot + c2[0];
  }
  
  return c2tot;
}

double chi2calc_point_iplane(int i, double xs, double ys, double c2[NPAR_CHI2])
{
  int j, k, ni, fp;
  int flag[NMAX_POIMG], kk[NMAX_POIMG];
  double dm, dis2;
  double f1, f2;
  double t1, t2;
  double rr[NMAX_POIMG][NPAR_IMAGE];

  for(j=0;j<NPAR_CHI2;j++) c2[j] = 0.0;

  if(obs_numimg[i] == 0){ return 0.0; }
  if(check_para_poi(i, xs, ys) > 0){ 
    c2[0] = c2[4] = chi2pen_range;
    c2[1] = c2[2] = c2[3] = 0.0;
    return chi2pen_range; 
  }

  for(j=0;j<NMAX_POIMG;j++) flag[j] = 0;

  findimg(xs, ys, para_poi[i][0], &ni, rr, 0);

  if(((ni != obs_numimg[i]) && (chi2_checknimg == 1)) || (ni < obs_numimg[i])){ 
    c2[0] = c2[1] = chi2pen_nimg;
    c2[2] = c2[3] = c2[4] = 0.0;
    return chi2pen_nimg; 
  }

  f1 = f2 = 0.0;
  t1 = t2 = 0.0;
  para_poi[i][3] = 1.0;
  para_poi[i][4] = 0.0;
  for(j=0;j<obs_numimg[i];j++){
    dm = DIS2_MIN_SET;
    fp = 0;
    for(k=0;k<ni;k++){
      if(flag[k] == 0){
	if((obs_parity[i][j] == 0) || (((double)obs_parity[i][j]) * rr[k][2]) > 0.0){
	  dis2 = (rr[k][0] - tab_obs[i][j][0]) * (rr[k][0] - tab_obs[i][j][0]) + (rr[k][1] - tab_obs[i][j][1]) * (rr[k][1] - tab_obs[i][j][1]);
	  if(dis2 < dm){ dm = dis2; kk[j] = k; fp = 1; }
	}
      }
    }
    if(fp == 0){
      c2[0] = c2[1] = chi2pen_nimg;
      c2[2] = c2[3] = c2[4] = 0.0;
      return chi2pen_nimg; 
    }
    flag[kk[j]] = 1;
    if(tab_obs[i][j][3] > 0.0) c2[1] = c2[1] + dm / (tab_obs[i][j][3] * tab_obs[i][j][3]);
    
    if(tab_obs[i][j][4] > 0.0){
      if(chi2_usemag == 0){
	f1 = f1 + fabs(tab_obs[i][j][2] * rr[kk[j]][2]) / (tab_obs[i][j][4] * tab_obs[i][j][4]);
	f2 = f2 + (rr[kk[j]][2] * rr[kk[j]][2]) / (tab_obs[i][j][4] * tab_obs[i][j][4]);
      } else {
	f1 = f1 + (tab_obs[i][j][2] + 2.5 * log10(fabs(rr[kk[j]][2]))) / (tab_obs[i][j][4] * tab_obs[i][j][4]);
	f2 = f2 + 1.0 / (tab_obs[i][j][4] * tab_obs[i][j][4]);
      }
    }
    
    if(tab_obs[i][j][6] > 0.0){
      t1 = t1 + (tab_obs[i][j][5] - rr[kk[j]][3]) / (tab_obs[i][j][6] * tab_obs[i][j][6]);
      t2 = t2 + 1.0 / (tab_obs[i][j][6] * tab_obs[i][j][6]);
    }
  }
  
  if(f2 > 0.0) para_poi[i][3] = f1 / f2;
  if(t2 > 0.0) para_poi[i][4] = t1 / t2;
  
  fp = 0;
  for(j=0;j<obs_numimg[i];j++){
    if((obs_parity[i][j] != 0) && ((((double)obs_parity[i][j]) * rr[kk[j]][2]) < 0.0)) fp = 1;
    if(tab_obs[i][j][4] > 0.0){
      if(chi2_usemag == 0){
	c2[2] = c2[2] + (fabs(tab_obs[i][j][2]) - fabs(rr[kk[j]][2]) * para_poi[i][3]) * (fabs(tab_obs[i][j][2]) - fabs(rr[kk[j]][2]) * para_poi[i][3]) / (tab_obs[i][j][4] * tab_obs[i][j][4]);
      } else {
	c2[2] = c2[2] + (tab_obs[i][j][2] + 2.5 * log10(fabs(rr[kk[j]][2])) - para_poi[i][3]) * (tab_obs[i][j][2] + 2.5 * log10(fabs(rr[kk[j]][2])) - para_poi[i][3]) / (tab_obs[i][j][4] * tab_obs[i][j][4]);
      }
    }
    if(tab_obs[i][j][6] > 0.0){ 
      c2[3] = c2[3] + (tab_obs[i][j][5] - rr[kk[j]][3] - para_poi[i][4]) * (tab_obs[i][j][5] - rr[kk[j]][3] - para_poi[i][4]) / (tab_obs[i][j][6] * tab_obs[i][j][6]);
    }
  }

  if(fp == 1) c2[2] = chi2pen_parity;
  c2[4] = chi2prior_point(i, xs, ys);

  c2[0] = c2[1] + c2[2] + c2[3] + c2[4];

  return c2[0];
}

double chi2calc_opt_func(double par[NDIMMAX + 1])
{
  double xs, ys;
  double c2[NPAR_CHI2];

  if(flag_para_poi[i_sav][1] == 1){
    xs = par[1];
    if(flag_para_poi[i_sav][2] == 1){
      ys = par[2];
    } else {
      ys = para_poi[i_sav][2];
    }
  } else {
    xs = para_poi[i_sav][1];
    if(flag_para_poi[i_sav][2] == 1){
      ys = par[1];
    } else {
      ys = para_poi[i_sav][2];
    }
  }
  
  if(chi2_point_splane == 0){
    return chi2calc_point_iplane(i_sav, xs, ys, c2);
  } else {
    return chi2calc_point_splane(i_sav, xs, ys, c2);
  }
}

/*--------------------------------------------------------------
  source plane chi2
*/

double chi2calc_opt_splane(double c2min[NMAX_POI][NPAR_CHI2], int verb)
{
  int i, l, f, nfunc;
  static double p[NDIMMAX + 2][NDIMMAX + 1];
  double c2tot, c2[NPAR_CHI2];
  double hh;
  double pout[NPAR_LMODEL];
  double y[NDIMMAX + 2]; 

  if(verb > 0){ 
    fprintf(stderr, "######## optimizing point sources (source plane)\n\n");
  }
  
  c2tot = 0.0;

  hh = dp_lev(maxlev - 1);

  for(i=0;i<num_poi;i++){
    ndim[i]=flag_para_poi[i][1] + flag_para_poi[i][2];
    i_sav = i;
    set_distance_lpl_zs(para_poi[i][0]);
    flag_scalc = 0;

    f = parmatch_poi(i);
 
    if((ndim[i] > 0) && (obs_numimg[i] >= 2) && (f == 0)){

      chi2calc_point_splane(i, 0.0, 0.0, c2);
    
      if(ndim[i] == 2){
	p[1][1] = umod[0];
	p[1][2] = umod[1];
	p[2][1] = umod[0] + hh;
	p[2][2] = umod[1];
	p[3][1] = umod[0];
	p[3][2] = umod[1] + hh;
	for(l=1;l<=3;l++) y[l] = chi2calc_point_splane(i, p[l][1], p[l][2], c2); 
      }
      
      if((ndim[i] == 1) && (flag_para_poi[i][1] == 1)){
	p[1][1] = umod[0];
	p[2][1] = umod[0] + hh;
	for(l=1;l<=2;l++) y[l] = chi2calc_point_splane(i, p[l][1], para_poi[i][2], c2); 
      }
      
      if((ndim[i] == 1) && (flag_para_poi[i][2] == 1)){
	p[1][1] = umod[1];
	p[2][1] = umod[1] + hh;
	for(l=1;l<=2;l++) y[l] = chi2calc_point_splane(i, para_poi[i][1], p[l][1], c2); 
      }
      
      simplex(p, y, ndim[i], tol_amoeba, chi2calc_opt_func, &nfunc, nmax_amoeba_point, verb); 
      
      nt[i] = nfunc;
      
      if(flag_para_poi[i][1] == 1){
	para_poi[i][1] = p[1][1]; 
	if(flag_para_poi[i][2] == 1){
	  para_poi[i][2] = p[1][2]; 
	} 
      } else {
	if(flag_para_poi[i][2] == 1){
	  para_poi[i][2] = p[1][1]; 
	} 
      }
    } else if((obs_numimg[i] == 1) && (f == 0)){
      lensmodel(tab_obs[i][0][0], tab_obs[i][0][1], pout, 1, 0);
      if(flag_para_poi[i][1] == 1){ para_poi[i][1] = tab_obs[i][0][0] - pout[0]; }
      if(flag_para_poi[i][2] == 1){ para_poi[i][2] = tab_obs[i][0][1] - pout[1]; }
    }

    chi2calc_point_splane(i, para_poi[i][1], para_poi[i][2], c2);
    nt[i]++;
    
    for(l=0;l<NPAR_CHI2;l++) c2min[i][l] = c2[l];
    c2tot = c2tot + c2[0];
  }

  return c2tot;
}

double chi2calc_point_splane(int i, double xs, double ys, double c2[NPAR_CHI2])
{
  int j, k, l, fp;
  double cc, f1, f2, t1, t2, dx, dy;
  double pout[NPAR_LMODEL];
  double m1, m2, mumod[NMAX_POIMG], tdmod[NMAX_POIMG];
  double aa[2][2], bb[2];
  static double a_mat[NMAX_POIMG][2][2];
  static double mu_mat[NMAX_POIMG][2][2];
  static double uobs[NMAX_POIMG][2];
  static double dmudx[NMAX_POIMG], dmudy[NMAX_POIMG];
  static double rr[NMAX_POIMG][4];
  static double def[NMAX_POIMG][2];
  
  double hh;

  hh = dp_lev(maxlev - 1) * 0.1;

  for(j=0;j<NPAR_CHI2;j++) c2[j] = 0.0;
    
  ndim[i] = flag_para_poi[i][1] + flag_para_poi[i][2];
    
  if(check_para_poi(i, xs, ys) > 0){ 
    c2[0] = c2[4] = chi2pen_range;
    c2[1] = c2[2] = c2[3] = 0.0;
    return chi2pen_range; 
  }
  
  if(flag_scalc == 0){
    flag_scalc = 1;
    for(k=0;k<obs_numimg[i];k++){
      lensmodel(tab_obs[i][k][0], tab_obs[i][k][1], pout, -1, 0);
      set_matrix(pout, a_mat[k], mu_mat[k]);
      uobs[k][0] = tab_obs[i][k][0] - pout[0];
      uobs[k][1] = tab_obs[i][k][1] - pout[1];
      rr[k][2] = 1.0 / (pout[6] + imag_ceil);
      rr[k][3] = pout[2];
      /* deflection angle for time delay */
      def[k][0] = def_lpl[nlp - 1][0];
      def[k][1] = def_lpl[nlp - 1][1];
      /*  computing derivative of mu */
      lensmodel(tab_obs[i][k][0] + 0.5 * hh, tab_obs[i][k][1], pout, 0, 0);
      m1 = pout[6];
      lensmodel(tab_obs[i][k][0] - 0.5 * hh, tab_obs[i][k][1], pout, 0, 0);
      m2 = pout[6];
      dmudx[k] = (-1.0) * rr[k][2] * rr[k][2] * (m1 - m2) / hh; 
      lensmodel(tab_obs[i][k][0], tab_obs[i][k][1] + 0.5 * hh, pout, 0, 0);
      m1 = pout[6];
      lensmodel(tab_obs[i][k][0], tab_obs[i][k][1] - 0.5 * hh, pout, 0, 0);
      m2 = pout[6];
      dmudy[k] = (-1.0) * rr[k][2] * rr[k][2] * (m1 - m2) / hh; 
    }
 
    for(k=0;k<2;k++) for(l=0;l<2;l++) aa[k][l] = 0.0;
    for(k=0;k<2;k++) bb[k] = 0.0;

    for(k=0;k<obs_numimg[i];k++){
      if(tab_obs[i][k][3] > 0.0){
	aa[0][0] = aa[0][0] + a_mat[k][0][0] / (tab_obs[i][k][3] * tab_obs[i][k][3]);
	aa[0][1] = aa[0][1] + a_mat[k][0][1] / (tab_obs[i][k][3] * tab_obs[i][k][3]);
	aa[1][0] = aa[1][0] + a_mat[k][1][0] / (tab_obs[i][k][3] * tab_obs[i][k][3]);
	aa[1][1] = aa[1][1] + a_mat[k][1][1] / (tab_obs[i][k][3] * tab_obs[i][k][3]);
	bb[0] = bb[0] + (a_mat[k][0][0] * uobs[k][0] + a_mat[k][0][1] * uobs[k][1]) / (tab_obs[i][k][3] * tab_obs[i][k][3]);
	bb[1] = bb[1] + (a_mat[k][1][0] * uobs[k][0] + a_mat[k][1][1] * uobs[k][1]) / (tab_obs[i][k][3] * tab_obs[i][k][3]);
      }
    }
    
    if(para_poi_sig[i][1] > 0.0){
      aa[0][0] = aa[0][0] + 1.0 / (para_poi_sig[i][1] * para_poi_sig[i][1]);
      bb[0] = bb[0] + para_poi_med[i][1] / (para_poi_sig[i][1] * para_poi_sig[i][1]);
    }

    if(para_poi_sig[i][2] > 0.0){
      aa[1][1] = aa[1][1] + 1.0 / (para_poi_sig[i][2] * para_poi_sig[i][2]);
      bb[1] = bb[1] + para_poi_med[i][2] / (para_poi_sig[i][2] * para_poi_sig[i][2]);
    }
    
    umod[0] = (aa[1][1] * bb[0] - aa[0][1] * bb[1]) / (aa[0][0] * aa[1][1] - aa[0][1] * aa[1][0]);
    umod[1] = (aa[0][0] * bb[1] - aa[1][0] * bb[0]) / (aa[0][0] * aa[1][1] - aa[0][1] * aa[1][0]);
  }
    
  f1 = f2 = 0.0;
  t1 = t2 = 0.0;
  for(k=0;k<obs_numimg[i];k++){
    if(tab_obs[i][k][3] > 0.0){
      cc = (a_mat[k][0][0] * (uobs[k][0] - xs) * (uobs[k][0] - xs) + (a_mat[k][0][1] + a_mat[k][1][0]) * (uobs[k][0] - xs) * (uobs[k][1] - ys) + a_mat[k][1][1] * (uobs[k][1] - ys) * (uobs[k][1] - ys)) / (tab_obs[i][k][3] * tab_obs[i][k][3]);
      if(cc < 0.0) cc = chi2pen_nimg;
      c2[1] = c2[1] + cc;
    }
    
    dx = mu_mat[k][0][0] * (uobs[k][0] - xs) + mu_mat[k][0][1] * (uobs[k][1] - ys);
    dy = mu_mat[k][1][0] * (uobs[k][0] - xs) + mu_mat[k][1][1] * (uobs[k][1] - ys);
    
    mumod[k] = rr[k][2] - (dmudx[k] * dx + dmudy[k] * dy);
    /* single lens plane */
    /* tdmod[k] = rr[k][3] - tdelay_fac(zl_lpl[0], dis_os, dis_ol_lpl[0], dis_ls_lpl[0]) * ((uobs[k][0] - tab_obs[i][k][0]) * (uobs[k][0] - xs) + (uobs[k][1] - tab_obs[i][k][1]) * (uobs[k][1] - ys)); */
    /* multiple lens plane */
    tdmod[k] = rr[k][3] - dis_tdelay[nlp - 1][nlp] * ((uobs[k][0] - (tab_obs[i][k][0] - def[k][0])) * (uobs[k][0] - xs) + (uobs[k][1] - (tab_obs[i][k][1] - def[k][1])) * (uobs[k][1] - ys));
      
    if(tab_obs[i][k][4] > 0.0){
      if(chi2_usemag == 0){
	f1 = f1 + fabs(tab_obs[i][k][2] * mumod[k]) / (tab_obs[i][k][4] * tab_obs[i][k][4]);
	f2 = f2 + (mumod[k] * mumod[k]) / (tab_obs[i][k][4] * tab_obs[i][k][4]);
      } else {
	f1 = f1 + (tab_obs[i][k][2] + 2.5 * log10(fabs(mumod[k]))) / (tab_obs[i][k][4] * tab_obs[i][k][4]);
	f2 = f2 + 1.0 / (tab_obs[i][k][4] * tab_obs[i][k][4]);
      }
    }
    
    if(tab_obs[i][k][6] > 0.0){
      t1 = t1 + (tab_obs[i][k][5] - tdmod[k]) / (tab_obs[i][k][6] * tab_obs[i][k][6]);
      t2 = t2 + 1.0 / (tab_obs[i][k][6] * tab_obs[i][k][6]);
    }
  }
    
  if(f2 > 0.0) para_poi[i][3] = f1 / f2;
  if(t2 > 0.0) para_poi[i][4] = t1 / t2;
  
  fp = 0;
  for(k=0;k<obs_numimg[i];k++){
    if((obs_parity[i][k] != 0) && ((((double)obs_parity[i][k]) * mumod[k]) < 0.0)) fp = 1;
    if(tab_obs[i][k][4] > 0.0){
      if(chi2_usemag == 0){
	c2[2] = c2[2] + (fabs(tab_obs[i][k][2]) - fabs(mumod[k]) * para_poi[i][3]) * (fabs(tab_obs[i][k][2]) - fabs(mumod[k]) * para_poi[i][3]) / (tab_obs[i][k][4] * tab_obs[i][k][4]);
      } else {
	c2[2] = c2[2] + (tab_obs[i][k][2] + 2.5 * log10(fabs(mumod[k])) - para_poi[i][3]) * (tab_obs[i][k][2] + 2.5 * log10(fabs(mumod[k])) - para_poi[i][3]) / (tab_obs[i][k][4] * tab_obs[i][k][4]);
      }
    } 
    if(tab_obs[i][k][6] > 0.0) c2[3] = c2[3] + (tab_obs[i][k][5] - tdmod[k] - para_poi[i][4]) * (tab_obs[i][k][5] - tdmod[k] - para_poi[i][4]) / (tab_obs[i][k][6] * tab_obs[i][k][6]);
  }
  
  if(fp == 1) c2[2] = chi2pen_parity;
  c2[4] = chi2prior_point(i, xs, ys);

  c2[0] = c2[1] + c2[2] + c2[3] + c2[4];
  
  return c2[0];

}

void set_matrix(double pout[NPAR_LMODEL], double mat_a[2][2], double mat_mu[2][2])
{
  double norm;

  norm = 1.0 / (pout[6] + imag_ceil);
  
  mat_a[0][0] = norm * norm * ((1.0 - pout[3] + pout[4]) * (1.0 - pout[3] + pout[4]) + (pout[5] - pout[7]) * (pout[5] - pout[7]));
  mat_a[0][1] = 2.0 * norm * norm * (pout[5] * (1.0 - pout[3]) + pout[4] * pout[7]);
  mat_a[1][0] = 2.0 * norm * norm * (pout[5] * (1.0 - pout[3]) + pout[4] * pout[7]);
  mat_a[1][1] = norm * norm * ((1.0 - pout[3] - pout[4]) * (1.0 - pout[3] - pout[4]) + (pout[5] + pout[7]) * (pout[5] + pout[7]));
  mat_mu[0][0] = norm * (1.0 - pout[3] + pout[4]);
  mat_mu[0][1] = norm * (pout[5] + pout[7]);
  mat_mu[1][0] = norm * (pout[5] - pout[7]);
  mat_mu[1][1] = norm * (1.0 - pout[3] - pout[4]);

  return;
}

/*--------------------------------------------------------------
  calculate scale length for amoeba
*/

double calcdelta_i(int i)
{
  int k, l;
  double c2, cc, sig, del, odel, dd;
  double pout[NPAR_LMODEL];
  double aa[2][2], bb[2];
  double a_mat[NMAX_POIMG][2][2];
  double mu_mat[NMAX_POIMG][2][2];
  double uobs[NMAX_POIMG][2];

  set_distance_lpl_init();
  set_distance_lpl_zs(para_poi[i][0]);

  if((flag_para_poi[i][1] + flag_para_poi[i][2]) > 0){

    for(k=0;k<obs_numimg[i];k++){
      lensmodel(tab_obs[i][k][0], tab_obs[i][k][1], pout, 0, 0);
      set_matrix(pout, a_mat[k], mu_mat[k]);
      uobs[k][0] = tab_obs[i][k][0] - pout[0];
      uobs[k][1] = tab_obs[i][k][1] - pout[1];
    }
 
    for(k=0;k<2;k++) for(l=0;l<2;l++) aa[k][l] = 0.0;
    for(k=0;k<2;k++) bb[k] = 0.0;

    for(k=0;k<obs_numimg[i];k++){
      if(tab_obs[i][k][3] > 0.0){
	aa[0][0] = aa[0][0] + a_mat[k][0][0] / (tab_obs[i][k][3] * tab_obs[i][k][3]);
	aa[0][1] = aa[0][1] + a_mat[k][0][1] / (tab_obs[i][k][3] * tab_obs[i][k][3]);
	aa[1][0] = aa[1][0] + a_mat[k][1][0] / (tab_obs[i][k][3] * tab_obs[i][k][3]);
	aa[1][1] = aa[1][1] + a_mat[k][1][1] / (tab_obs[i][k][3] * tab_obs[i][k][3]);
	bb[0] = bb[0] + (a_mat[k][0][0] * uobs[k][0] + a_mat[k][0][1] * uobs[k][1]) / (tab_obs[i][k][3] * tab_obs[i][k][3]);
	bb[1] = bb[1] + (a_mat[k][1][0] * uobs[k][0] + a_mat[k][1][1] * uobs[k][1]) / (tab_obs[i][k][3] * tab_obs[i][k][3]);
      }
    }
    
    umod[0] = (aa[1][1] * bb[0] - aa[0][1] * bb[1]) / (aa[0][0] * aa[1][1] - aa[0][1] * aa[1][0]);
    umod[1] = (aa[0][0] * bb[1] - aa[1][0] * bb[0]) / (aa[0][0] * aa[1][1] - aa[0][1] * aa[1][0]);

    c2 = 0.0;
    sig = 0.0;
    for(k=0;k<obs_numimg[i];k++){
      if(tab_obs[i][k][3] > 0.0){
	cc = (a_mat[k][0][0] * (uobs[k][0] - umod[0]) * (uobs[k][0] - umod[0]) + (a_mat[k][0][1] + a_mat[k][1][0]) * (uobs[k][0] - umod[0]) * (uobs[k][1] - umod[1]) + a_mat[k][1][1] * (uobs[k][1] - umod[1]) * (uobs[k][1] - umod[1]));
	if(cc < 0.0) cc = chi2pen_nimg;
	c2 = c2 + cc;
	sig = sig + 1.0;
      }
    }
    
    del = sqrt(c2 / (sig + OFFSET_DEL));
    odel = calcdelta_obs(i);


    dd = del / odel;

    if(dd < amoeba_delmin) dd = amoeba_delmin;
    if(dd > amoeba_delmax) dd = amoeba_delmax;
  } else {
    dd = amoeba_delmin;
  }
    
  return dd;
}

double calcdelta_obs(int i)
{
  int k, l;
  double d2, d2max;

  d2max = 0.0;
  for(k=0;k<obs_numimg[i];k++){
      for(l=k+1;l<obs_numimg[i];l++){
	d2 = (tab_obs[i][k][0] - tab_obs[i][l][0]) * (tab_obs[i][k][0] - tab_obs[i][l][0]) + (tab_obs[i][k][1] - tab_obs[i][l][1]) * (tab_obs[i][k][1] - tab_obs[i][l][1]);
	if(d2 > d2max) d2max = d2;
      }
    }
  
  return sqrt(d2max) / 2.0;
    
}

/*--------------------------------------------------------------
  prior
*/

int check_para_poi(int i, double xs, double ys)
{
  int r = 0;

  if((xs < para_poi_min[i][1]) || (xs > para_poi_max[i][1])) r++;
  if((ys < para_poi_min[i][2]) || (ys > para_poi_max[i][2])) r++;

  return r;
}

int check_para_poi_zs(int i)
{
  if((para_poi[i][0] < para_poi_min[i][0]) || (para_poi[i][0] > para_poi_max[i][0])){
    return 1;
  } else {
    return 0;
  }
    
}

int check_para_poi_all(void)
{
  int i;
  int r = 0;

  for(i=0;i<num_poi;i++){
    r = r + check_para_poi(i, para_poi[i][1], para_poi[i][2]);
  }
  
  return r;
}

double chi2prior_point(int i, double xs, double ys)
{
  double c2;

  c2 = 0.0;
  if(para_poi_sig[i][1] > 0.0){
    c2 = c2 + (para_poi_med[i][1] - xs) * (para_poi_med[i][1] - xs) / (para_poi_sig[i][1] * para_poi_sig[i][1]);
  }
  if(para_poi_sig[i][2] > 0.0){
    c2 = c2 + (para_poi_med[i][2] - ys) * (para_poi_med[i][2] - ys) / (para_poi_sig[i][2] * para_poi_sig[i][2]);
  }
  
  return c2;
}

int parmatch_poi(int i)
{
  int j, f;

  f = 0;
  for(j=1;j<NPAR_POI;j++){
    /* 
      if((para_poi_rai[i][j] != i) || (para_poi_raj[i][j] != j)){
      if(para_poi_ras[i][j] == 0.0){
      para_poi[i][j] = para_poi_rat[i][j] * para_poi[para_poi_rai[i][j]][para_poi_raj[i][j]];
      }
      }
    */
    if((para_poilen_rai[i][j] >= 0) && (para_poilen_raj[i][j] >= 0)){
      if(para_poilen_ras[i][j] == 0.0){
	para_poi[i][j] = para_poilen_rat[i][j] * para_lens[para_poilen_rai[i][j]][para_poilen_raj[i][j]];
	f++;
      }
    }
  }

  if(f == 1) terminator("match poilens should be used for both x and y");
  
  return f;
}
