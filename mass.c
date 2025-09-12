#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>

#include "glafic.h"

#define NM 25
static char lmodelname[NM][10]={
"gals",     /* 1 */
"nfwpot",   /* 2 */
"sie",      /* 3 */
"jaffe",    /* 4 */
"point",    /* 5 */
"pert",     /* 6 */
"clus3",    /* 7 */
"mpole",    /* 8 */
"hernpot",  /* 9 */
"nfw",      /* 10 */
"hern",     /* 11 */
"powpot",   /* 12 */
"pow",      /* 13 */
"gnfwpot",  /* 14 */
"gnfw",     /* 15 */
"serspot",  /* 16 */
"sers",     /* 17 */
"tnfwpot",  /* 18 */
"tnfw",     /* 19 */
"einpot",   /* 20 */
"ein",      /* 21 */
"anfw",     /* 22 */
"ahern",    /* 23 */
"crline",   /* 24 */
"gaupot"    /* 25 */
};

static double alpha_gnfw_sav, x_gnfw_sav;
static double alpha_ein_sav, x_ein_sav;
static double gam_pow_sav;
static double n_sers_sav;
static double tau_tnfw_sav;

static int n_integ_sav;
static double x_integ_sav, y_integ_sav, q_integ_sav;
static int ell_u_lnint;
static double ell_u_lnmin;
static double (*func_sav)(double);

/*--------------------------------------------------------------
  lensmodel

  pout[0]: ax
  pout[1]: ay
  pout[2]: tdelay
  pout[3]: kap
  pout[4]: gam1
  pout[5]: gam2
  pout[6]: muinv
  pout[7]: rot
*/

void lensmodel(double tx, double ty, double pout[NPAR_LMODEL], int alponly, int lensid)
{
  double par_multi[NMAX_LPL + 1][NPAR_MULTI];
  /* npar_multi: 0 - ty  1 - tx  2 - alp1  3 - alp2  
                 4 - U_11  5 - U_12  6 - U_21  7 - U_22  
                 8 - A_11  9 - A_12  10 - A_21  11 - A_22  12 - phi */
  double ax, ay, phi, kap, gam1, gam2;
  
  if((lensid > num_len) || (lensid < 0)) terminator("wrong lens id"); 
  
  if(lensid != 0){
    lensmodel_get_i(lensid - 1, tx, ty, &kap, &gam1, &gam2, &phi, &ax, &ay, alponly);
    pout[0] = ax;
    pout[1] = ay;
    pout[2] = tdelay_fac(zl_ext, dis_os, dis_ol, dis_ls) * (0.5 * (ax * ax + ay * ay) - phi);
    pout[3] = kap;
    pout[4] = gam1;
    pout[5] = gam2;
    pout[6] = (1.0 - kap) * (1.0 - kap) - (gam1 * gam1 + gam2 * gam2);
    pout[7] = 0.0;
  } else {
    multi_calc_all(tx, ty, par_multi, alponly);
    multi_set_pout(par_multi, pout);
  }
  
  return;
}

void multi_calc_all(double tx, double ty, double par_multi[NMAX_LPL + 1][NPAR_MULTI], int alponly)
{
  int i, j;
  double tax, tay, tk, tg1, tg2, tph;
  double ax, ay, phi, kap, gam1, gam2;
  
  par_multi[0][0] = tx;
  par_multi[0][1] = ty;
  par_multi[0][8]  = 1.0;
  par_multi[0][9]  = 0.0;
  par_multi[0][10] = 0.0;
  par_multi[0][11] = 1.0;
  def_lpl[0][0] = 0.0;
  def_lpl[0][1] = 0.0;
  def_lpl[0][2] = 1.0;
  def_lpl[0][3] = 0.0;
  def_lpl[0][4] = 0.0;
  def_lpl[0][5] = 1.0;
  
  for(j=0;j<nlp;j++){
    lensmodel_sum_init(&ax, &ay, &phi, &kap, &gam1, &gam2);
    
    for(i=0;i<num_len;i++){
      if(lens_lpl_id[i] == j){
	lensmodel_get_i(i, par_multi[j][0], par_multi[j][1], &tk, &tg1, &tg2, &tph, &tax, &tay, alponly);
	lensmodel_sum(&ax, &ay, &phi, &kap, &gam1, &gam2, tax, tay, tph, tk, tg1, tg2, alponly);
      }
    }
    
    update_par_multi(ax, ay, phi, kap, gam1, gam2, j, par_multi);
  }

  return;
}

void multi_set_pout(double par_multi[NMAX_LPL + 1][NPAR_MULTI], double pout[NPAR_LMODEL])
{
  int i, k;
  double dd;
  double mat[4], par_out[4];
  
  for(i=0;i<NPAR_LMODEL;i++)  pout[i] = 0.0; 
  for(k=0;k<4;k++) mat[k] = 0.0;

  for(i=0;i<nlp;i++){
    pout[0] = pout[0] + par_multi[i][2];
    pout[1] = pout[1] + par_multi[i][3];
  }
  
  par_multi[nlp][0] = par_multi[0][0] - pout[0];
  par_multi[nlp][1] = par_multi[0][1] - pout[1];
  
  for(i=0;i<nlp;i++){
    dd = (par_multi[i][0] - par_multi[i+1][0]) * (par_multi[i][0] - par_multi[i+1][0]) + (par_multi[i][1] - par_multi[i+1][1]) * (par_multi[i][1] - par_multi[i+1][1]);
    pout[2] = pout[2] + dis_tdelay[i][i+1] * (0.5 * dd - dis_beta[i][i+1] * par_multi[i][12]);
  }
  
  for(i=0;i<nlp;i++){
    matrix_ua_calc(par_multi, i, par_out);
    for(k=0;k<4;k++) mat[k] = mat[k] + par_out[k];
  }

  pout[3] = 0.5 * (mat[0] + mat[3]);
  pout[4] = 0.5 * (mat[0] - mat[3]);
  pout[5] = 0.5 * (mat[1] + mat[2]);
  pout[6] = (1.0 - mat[0]) * (1.0 - mat[3]) - mat[1] * mat[2];
  /* rotation */
  pout[7] = 0.5 * (mat[1] - mat[2]); 
     
  return;
}

void update_par_multi(double ax, double ay, double phi, double kap, double gam1, double gam2, int j, double par_multi[NMAX_LPL + 1][NPAR_MULTI])
{
  int i, k;
  double alp1, alp2;
  double mat[4], par_out[4];
  
  par_multi[j][2] = ax;
  par_multi[j][3] = ay;
  par_multi[j][4] = kap + gam1;
  par_multi[j][5] = gam2;
  par_multi[j][6] = gam2;
  par_multi[j][7] = kap - gam1;
  par_multi[j][12] = phi;

  if(j < (nlp - 1)){
    alp1 = 0.0;
    alp2 = 0.0;
    for(k=0;k<4;k++) mat[k] = 0.0;

    for(i=0;i<=j;i++){
      alp1 = alp1 + dis_beta[i][j+1] * par_multi[i][2];
      alp2 = alp2 + dis_beta[i][j+1] * par_multi[i][3];
      
      matrix_ua_calc(par_multi, i, par_out);
      for(k=0;k<4;k++) mat[k] = mat[k] + dis_beta[i][j+1] * par_out[k];
    }

    par_multi[j + 1][0] = par_multi[0][0] - alp1;
    par_multi[j + 1][1] = par_multi[0][1] - alp2 ;
    par_multi[j + 1][8]  = 1.0 - mat[0];
    par_multi[j + 1][9]  = (-1.0) *  mat[1];
    par_multi[j + 1][10] = (-1.0) *  mat[2];
    par_multi[j + 1][11] = 1.0 - mat[3];

    def_lpl[j + 1][0] = alp1;
    def_lpl[j + 1][1] = alp2;
    def_lpl[j + 1][2] = par_multi[j + 1][8];
    def_lpl[j + 1][3] = par_multi[j + 1][9];
    def_lpl[j + 1][4] = par_multi[j + 1][10];
    def_lpl[j + 1][5] = par_multi[j + 1][11];
  }
  
  return;
}

void matrix_ua_calc(double par_multi[NMAX_LPL + 1][NPAR_MULTI], int i, double par_out[4])
{
  par_out[0] = par_multi[i][4] * par_multi[i][8] +  par_multi[i][5] * par_multi[i][10];
  par_out[1] = par_multi[i][4] * par_multi[i][9] +  par_multi[i][5] * par_multi[i][11];
  par_out[2] = par_multi[i][6] * par_multi[i][8] +  par_multi[i][7] * par_multi[i][10];
  par_out[3] = par_multi[i][6] * par_multi[i][9] +  par_multi[i][7] * par_multi[i][11];

  return;
}

void lensmodel_get_i(int i, double tx, double ty, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  double cx, cy;
  
  *kap  = 0.0;
  *gam1 = 0.0;
  *gam2 = 0.0;
  *phi  = 0.0;
  *ax   = 0.0;
  *ay   = 0.0;
  
  set_distance_lpl_i(lens_lpl_id[i]);

  cx = para_lens[i][2];
  cy = para_lens[i][3];

  switch(model_lens[i]){
    
  case 1:
    kapgam_gals(tx, ty, para_lens[i][1], para_lens[i][6], para_lens[i][7], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
    
  case 2:
    kapgam_nfwpot(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
    
  case 3:
    kapgam_sie(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
    
  case 4:
    kapgam_jaffe(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
    
  case 5:
    kapgam_point(tx, ty, cx, cy, para_lens[i][1], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
    
  case 6:
    kapgam_pert(tx, ty, cx, cy, para_lens[i][1], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
	
  case 7:
    kapgam_clus3(tx, ty, cx, cy, para_lens[i][1], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
    
  case 8:
    kapgam_mpole(tx, ty, cx, cy, para_lens[i][1], para_lens[i][4], para_lens[i][5], para_lens[i][6], para_lens[i][7], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
	
  case 9:
    kapgam_hernpot(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
	
  case 10:
    kapgam_nfw(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;

  case 11:
    kapgam_hern(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
	
  case 12:
    kapgam_powpot(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
	
  case 13:
    kapgam_pow(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
	
  case 14:
    kapgam_gnfwpot(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
	
  case 15:
    kapgam_gnfw(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;

  case 16:
    kapgam_serspot(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;

  case 17:
    kapgam_sers(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;

  case 18:
    kapgam_tnfwpot(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;

  case 19:
    kapgam_tnfw(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;

  case 20:
    kapgam_einpot(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
	
  case 21:
    kapgam_ein(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;

  case 22:
    kapgam_anfw(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;

  case 23:
    kapgam_ahern(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;

  case 24:
    kapgam_crline(tx, ty, cx, cy, para_lens[i][1], para_lens[i][7], para_lens[i][5], para_lens[i][6], kap, gam1, gam2, phi, ax, ay, alponly);
    break;

  case 25:
    kapgam_gaupot(tx, ty, cx, cy, para_lens[i][1], para_lens[i][6], para_lens[i][7], para_lens[i][4], para_lens[i][5], kap, gam1, gam2, phi, ax, ay, alponly);
    break;
  }

  return;
}

void lensmodel_sum(double *ax, double *ay, double *phi, double *kap, double *gam1, double *gam2, double tax, double tay, double tph, double tk, double tg1, double tg2, int alponly)
{
  (*ax) = (*ax) + tax;
  (*ay) = (*ay) + tay;

  if(alponly != 1){
    (*kap) = (*kap) + tk;
    (*gam1) = (*gam1) + tg1;
    (*gam2) = (*gam2) + tg2;
    if(alponly < 0) (*phi) = (*phi) + tph;
  }

  return;
}

void lensmodel_sum_init(double *ax, double *ay, double *phi, double *kap, double *gam1, double *gam2)
{
  (*ax)   = 0.0;
  (*ay)   = 0.0;
  (*kap)  = 0.0;
  (*gam1) = 0.0;
  (*gam2) = 0.0;
  (*phi)  = 0.0;

  return;
}

int lmodeltoint(char *model)
{
  int i;

  for(i=0;i<NM;i++){
    if(strcmp(model, lmodelname[i]) == 0){
    return i + 1;
    } 
  }
   
  return 0;
}

char* inttolmodel(int i)
{
  return lmodelname[i - 1];
}

/*--------------------------------------------------------------
  external perturbation
*/

void kapgam_pert(double tx, double ty, double tx0, double ty0, double zs_fid, double k, double g, double tg, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double co, si, tgtg;
  double fac;

  fac = fac_pert(zs_fid);

  if((ff != 1) || (tgtg != tg)){
    ff = 1;
    tgtg = tg;
    co = cos(2.0 * (tg - 90.0) * M_PI / 180.0);
    si = sin(2.0 * (tg - 90.0) * M_PI / 180.0);
  }

  *ax = fac * ((tx - tx0) * k - (tx - tx0) * g * co - (ty - ty0) * g * si);
  *ay = fac * ((ty - ty0) * k + (ty - ty0) * g * co - (tx - tx0) * g * si);

  if(alponly != 1){
    *kap = k * fac;
    *gam1 = (-1.0) * g * fac * co;
    *gam2 = (-1.0) * g * fac * si;
    if(alponly < 0) *phi = 0.5 * ((tx - tx0) * (tx - tx0) + (ty - ty0) * (ty - ty0)) * (*kap) + 0.5 * ((tx - tx0) * (tx - tx0) - (ty - ty0) * (ty - ty0)) * (*gam1) + (tx - tx0) * (ty - ty0) * (*gam2);
  }

  return;
}

void kapgam_clus3(double tx, double ty, double tx0, double ty0, double zs_fid, double g, double tg, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double cog, sig, cog3, sig3, tgtg;
  double fac, r, cox, six, cox3, six3, co, si, co3, si3, pxx, pyy, pxy;
  
  fac = fac_pert(zs_fid);

  r = sqrt((tx - tx0) * (tx - tx0) + (ty - ty0) * (ty - ty0))
     + smallcore; /* avoid error at the center */
  /* t = atan2(ty - ty0, tx - tx0 + smallcore) * 180.0 / M_PI; */
  cox = (tx - tx0) / r;
  six = (ty - ty0) / r;
  cox3 = cox * (cox * cox - 3.0 * six * six);
  six3 = (3.0 * cox * cox - six * six) * six;

  if((ff != 1) || (tgtg != tg)){
    ff = 1;
    tgtg = tg;
    cog = cos(tg * M_PI / 180.0);
    sig = sin(tg * M_PI / 180.0);
    /* gravlens uses the following angle (?) */
    /* cog = cos((tg + 90.0) * M_PI / 180.0);
       sig = sin((tg + 90.0) * M_PI / 180.0); */
    cog3 = cog * (cog * cog - 3.0 * sig * sig);
    sig3 = (3.0 * cog * cog - sig * sig) * sig;
  }

  co = cox * cog + six * sig;
  si = six * cog - cox * sig;
  co3 = cox3 * cog3 + six3 * sig3;
  si3 = six3 * cog3 - cox3 * sig3;

  *ax = fac * (g / 4.0) * (3.0 * r * (tx - tx0) * (si + si3) - r * (ty - ty0) * (co + 3.0 * co3));
  *ay = fac * (g / 4.0) * (3.0 * r * (ty - ty0) * (si + si3) + r * (tx - tx0) * (co + 3.0 * co3));

  if(alponly != 1){
    pxx = fac * (g / 4.0) * (3.0 * (r + (tx - tx0) * (tx - tx0) / r) * (si + si3) - 4.0 * ((tx - tx0) * (ty - ty0) / r) * (co + 3.0 * co3) - ((ty - ty0) * (ty - ty0) / r) * (si + 9.0 * si3));
    pyy = fac * (g / 4.0) * (3.0 * (r + (ty - ty0) * (ty - ty0) / r) * (si + si3) + 4.0 * ((tx - tx0) * (ty - ty0) / r) * (co + 3.0 * co3) - ((tx - tx0) * (tx - tx0) / r) * (si + 9.0 * si3));
    pxy = fac * (g / 4.0) * (3.0 * ((tx - tx0) * (ty - ty0) / r) * (si + si3) + 2.0 * (((tx - tx0) * (tx - tx0) - (ty - ty0) * (ty - ty0)) / r) * (co + 3.0 * co3) + ((tx - tx0) * (ty - ty0) / r) * (si + 9.0 * si3));

    *kap = 0.5 * (pxx + pyy);
    *gam1 = 0.5 * (pxx - pyy);
    *gam2 = pxy;
    if(alponly < 0) *phi = fac * (g / 4.0) * r * r * r * (si + si3);
  }

  return;
}

void kapgam_mpole(double tx, double ty, double tx0, double ty0, double zs_fid, double g, double tg, double m, double n, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  double fac, f2, r, t, co, si, pxx, pyy, pxy;

  checkmodelpar_min(m, 0.0);

  fac = fac_pert(zs_fid);

  r = sqrt((tx - tx0) * (tx - tx0) + (ty - ty0) * (ty - ty0))
     + smallcore; /* avoid error at the center */
  t = atan2(ty - ty0, tx - tx0 + smallcore) * 180.0 / M_PI;

  co = cos(m * (t - tg - 90.0) * M_PI / 180.0);
  si = sin(m * (t - tg - 90.0) * M_PI / 180.0);

  f2 = (-1.0) * g * pow(r, n - 2.0) / m;

  *ax = fac * f2 * (n * (tx - tx0) * co + m * (ty - ty0) * si);
  *ay = fac * f2 * (n * (ty - ty0) * co - m * (tx - tx0) * si);

  if(alponly != 1){
    pxx = fac * f2 * ((n - 2.0) * (tx - tx0) * (n * (tx - tx0) * co + m * (ty - ty0) * si) + n * r * r * co + n * m * (tx - tx0) * (ty - ty0) * si - m * m * (ty - ty0) * (ty - ty0) * co) / (r * r);
    pyy = fac * f2 * ((n - 2.0) * (ty - ty0) * (n * (ty - ty0) * co - m * (tx - tx0) * si) + n * r * r * co - n * m * (tx - tx0) * (ty - ty0) * si - m * m * (tx - tx0) * (tx - tx0) * co) / (r * r);
    pxy = fac * f2 * ((n * (n - 2.0) + m * m) * (tx - tx0) * (ty - ty0) * co + m * (n - 1.0) * ((ty - ty0) * (ty - ty0) - (tx - tx0) * (tx - tx0)) * si) / (r * r);

    *kap = 0.5 * (pxx + pyy);
    *gam1 = 0.5 * (pxx - pyy);
    *gam2 = pxy;
    if(alponly < 0) *phi = fac * f2 * r * r * co;
  }

  return;
}

double fac_pert(double zs_fid)
{
  static double fac, zsf, ll, dd;

  if((zsf != zs_fid) || (dd != dis_os) || (ll != dis_ls)){
    if(zl_ext >= zs_fid) terminator("wrong zs_fid\n");
    dd = dis_os;
    ll = dis_ls;
    zsf = zs_fid;
    fac = (dis_ls / dis_os) / (dis_angulard(zl_ext, zs_fid) / dis_angulard(0.0, zs_fid));
  }

  return fac;
}

/*--------------------------------------------------------------
  straight critical line
*/

void kapgam_crline(double tx, double ty, double tx0, double ty0, double zs_fid, double k, double pa, double epsilon, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double g, dx, dy;
  double fac;

  fac = fac_pert(zs_fid);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  g = 1.0 - k;
  
  dx = (tx - tx0) * co - (ty - ty0) * si;
  dy = (tx - tx0) * si + (ty - ty0) * co;
  
  *ax = fac * (dx * k + dx * g - 0.5 * epsilon * dx * dx);
  *ay = fac * (dy * k - dy * g);

  if(alponly != 1){
    *kap = fac * (k - 0.5 * epsilon * dx);
    *gam1 = fac * (g - 0.5 * epsilon * dx);
    *gam2 = 0.0;
    if(alponly < 0) *phi = fac * (0.5 * (dx * dx + dy * dy) * k + 0.5 * (dx * dx - dy * dy) * g - (1.0 / 6.0) * epsilon * dx * dx * dx);
  }

  return;
}

/*--------------------------------------------------------------
  galaxies (sum of jaffe)
*/

void kapgam_gals(double tx, double ty, double sig, double a, double alpha, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  int i;
  double dx, dy, ax1, ax2, ay1, ay2, ak1, ak2, ag1, ag2, ar1, ar2, ap1, ap2;
  static double bb[NMAX_GALS], aa[NMAX_GALS], si[NMAX_GALS], co[NMAX_GALS], qq[NMAX_GALS];
  static double sss, aaa, dd, ll;
  
  if((sss != sig) || (aaa != a) || (dd != dis_os) || (ll != dis_ls)){
    dd = dis_os;
    ll = dis_ls;
    sss = sig;
    aaa = a;
    for(i=0;i<num_gal;i++){
      qq[i] = 1.0 - para_gals[i][3];
      bb[i] = b_sie(pow(para_gals[i][2], 0.25) * sig, qq[i]);
      aa[i] = pow(para_gals[i][2], alpha) * a;
      si[i] = sin((-1.0) * (para_gals[i][4] - 90.0) * M_PI / 180.0);
      co[i] = cos((-1.0) * (para_gals[i][4] - 90.0) * M_PI / 180.0);
    }
  }

  *ax = 0.0;
  *ay = 0.0;
  
  for(i=0;i<num_gal;i++){
    /* dx = tx - (para_gals[i][0] - gals_offset_x[i]); 
       dy = ty - (para_gals[i][1] - gals_offset_y[i]); */
    dx = tx - para_gals[i][0];
    dy = ty - para_gals[i][1];
    
    alpha_sie_bq(dx, dy, bb[i], smallcore, qq[i], si[i], co[i], &ax1, &ay1);
    alpha_sie_bq(dx, dy, bb[i], aa[i], qq[i], si[i], co[i], &ax2, &ay2);
    
    (*ax) = (*ax) + (ax1 - ax2);
    (*ay) = (*ay) + (ay1 - ay2);
  }

  if(alponly != 1){
    *kap = 0.0;
    *gam1 = 0.0;
    *gam2 = 0.0;
    *phi = 0.0;
    
    for(i=0;i<num_gal;i++){
      /* dx = tx - (para_gals[i][0] - gals_offset_x[i]);
	 dy = ty - (para_gals[i][1] - gals_offset_y[i]); */
      dx = tx - para_gals[i][0];
      dy = ty - para_gals[i][1];
      
      kapgam_sie_bq(dx, dy, bb[i], smallcore, qq[i], si[i], co[i], &ak1, &ag1, &ar1, &ap1, alponly);
      kapgam_sie_bq(dx, dy, bb[i], aa[i], qq[i], si[i], co[i], &ak2, &ag2, &ar2, &ap2, alponly);
      
      (*kap) = (*kap) + (ak1 - ak2);
      (*gam1) = (*gam1) + (ag1 - ag2);
      (*gam2) = (*gam2) + (ar1 - ar2);
      if(alponly < 0) (*phi) = (*phi) + (ap1 - ap2);
    }
  }

  return;
}


/*--------------------------------------------------------------
  point mass
*/

void kapgam_point(double tx, double ty, double tx0, double ty0, double m, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  double r2, rr, dx, dy, re2;

  checkmodelpar_mineq(m, 0.0);

  re2 = re2_point(m);
  dx = tx - tx0;
  dy = ty - ty0;
  r2 = dx * dx + dy * dy;

  rr = re2 / (r2 + smallcore * smallcore);
  *ax = rr * dx;
  *ay = rr * dy;

  if(alponly != 1){
    if(r2 < (smallcore * smallcore)){
      r2 = smallcore * smallcore;
      *kap = 0.0;
      *gam1 = re2 / (2.0 * smallcore * smallcore);
      *gam2 = re2 / (2.0 * smallcore * smallcore);
    } else {
      *kap = 0.0;
      *gam1 = (re2 / (r2 * r2)) * (dy * dy - dx * dx);
      *gam2 = (re2 / (r2 * r2)) * ((-2.0) * dx * dy);
    }
    if(alponly < 0) *phi = re2 * 0.5 * log(r2);
  }

  return;
}

double re2_point(double m)
{
  double d;

  d = dis_ls / ( COVERH_MPCH * dis_ol * dis_os );

  return ( 2.0 * ( R_SCHWARZ * m / MPC2METER ) * d ) / ( ARCSEC2RADIAN * ARCSEC2RADIAN );
}

/*--------------------------------------------------------------
  pseudo-Jaffe profile
*/

void kapgam_jaffe(double tx, double ty, double tx0, double ty0, double sig, double a, double rco, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double si, co, papa;
  double bb, q, dx, dy, ax1, ay1, ax2, ay2, k1, k2, g1, g2, r1, r2, p1, p2;

  checkmodelpar_mineq(sig, 0.0);
  checkmodelpar_mineq(rco, 0.0);
  checkmodelpar_min(a, 0.0);

  if(rco < smallcore) rco = smallcore;

  /* checkmodelpar_min(a, rco); */
  if(a > rco){
    checkmodelpar_mineq(e, 0.0);
    checkmodelpar_max(e, 1.0);
    
    q = 1.0 - e;
    bb = b_sie(sig, q);
    
    dx = tx - tx0;
    dy = ty - ty0;
    
    if((ff != 1) || (papa != pa)){
      ff = 1;
      papa = pa;
      si = sin((-1.0) * (pa - 90.0) * M_PI / 180.0);
      co = cos((-1.0) * (pa - 90.0) * M_PI / 180.0);
    }
    
    alpha_sie_bq(dx, dy, bb, rco, q, si, co, &ax1, &ay1);
    alpha_sie_bq(dx, dy, bb, a, q, si, co, &ax2, &ay2);
    
    *ax = ax1 - ax2;
    *ay = ay1 - ay2;

    if(alponly != 1){
      kapgam_sie_bq(dx, dy, bb, rco, q, si, co, &k1, &g1, &r1, &p1, alponly);
      kapgam_sie_bq(dx, dy, bb, a, q, si, co, &k2, &g2, &r2, &p2, alponly);
      
      *kap = k1 - k2;
      *gam1 = g1 - g2;
      *gam2 = r1 - r2;
      if(alponly < 0) *phi = p1 - p2;
    }
  } else {
    *kap  = 0.0;
    *gam1 = 0.0;
    *gam2 = 0.0;
    *phi  = 0.0;
    *ax   = 0.0;
    *ay   = 0.0;
  }
  return;

}

/*--------------------------------------------------------------
  SIE profile
*/

void kapgam_sie(double tx, double ty, double tx0, double ty0, double sig, double s, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double si, co, papa;
  double q, bb;
  
  checkmodelpar_mineq(sig, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_mineq(s, 0.0);
  
  q = 1.0 - e;
  bb = b_sie(sig, q);
  if(s < smallcore) s = smallcore;  /* avoid error at the center */

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * (pa - 90.0) * M_PI / 180.0);
    co = cos((-1.0) * (pa - 90.0) * M_PI / 180.0);
  }

  alpha_sie_bq(tx - tx0, ty - ty0, bb, s, q, si, co, ax, ay);
  
  if(alponly != 1){
    kapgam_sie_bq(tx - tx0, ty - ty0, bb, s, q, si, co, kap, gam1, gam2, phi, alponly);
  }

  return;
}

void kapgam_sie_bq(double dx, double dy, double bb, double s, double q, double si, double co, double *kap, double *gam1, double *gam2, double *phi, int alponly)
{
  double ss;
  double ddx, ddy;
  double pxx, pyy, pxy, rpxx, rpyy, rpxy, p;

  ss = s * facq_sie(q);

  ddx = co * dx - si * dy;
  ddy = si * dx + co * dy;

  if(alponly < 0){ 
    p = phi_sie_dl(ddx, ddy, ss, q);
    *phi = bb * p;
  }
  ddphi_sie_dl(ddx, ddy, ss, q, &pxx, &pxy, &pyy);

  rpxx = co * co * pxx + 2.0 * co * si * pxy + si * si * pyy;
  rpyy = si * si * pxx - 2.0 * co * si * pxy + co * co * pyy;
  rpxy = si * co * (pyy - pxx) + (co * co - si * si) * pxy;

  *kap = 0.5 * bb * (rpxx + rpyy);
  *gam1 = 0.5 * bb * (rpxx - rpyy);
  *gam2 = bb * rpxy;

  return;
}

void alpha_sie_bq(double dx, double dy, double bb, double s, double q, double si, double co, double *ax, double *ay)
{
  double ss;
  double ddx, ddy, aax, aay;

  ss = s * facq_sie(q);

  ddx = co * dx - si * dy;
  ddy = si * dx + co * dy;

  alpha_sie_dl(ddx, ddy, ss, q, &aax, &aay);

  *ax = bb * (aax * co + aay * si);
  *ay = bb * (aax * (-1.0) * si + aay * co);
  
  return;
}

double phi_sie_dl(double x, double y, double s, double q)
{
  double psi, px, py, aa;

  alpha_sie_dl(x, y, s, q, &px, &py);
  
  psi = sqrt(q * q * (s * s + x * x) + y * y);
  aa = log((1.0 + q) * s) - log(sqrt((psi + s) * (psi + s) + (1.0 - q * q) * x * x));

  return x * px + y * py + q * s * aa;
  
}

void ddphi_sie_dl(double x, double y, double s, double q, double *pxx, double *pxy, double *pyy)
{
  double psi, f;

  psi = sqrt(q * q * (s * s + x * x) + y * y);
  f = (1.0 + q * q) * s * s + 2.0 * psi * s + x * x + y * y;

  *pxx = (q / psi) * (q * q * s * s + y * y + s * psi) / f;
  *pyy = (q / psi) * (s * s + x * x + s * psi) / f;
  *pxy = (q / psi) * ((-1.0) * x * y) / f;

  return;
}

void alpha_sie_dl(double x, double y, double s, double q, double *ax, double *ay)
{
  double sq, psi;

  if((1.0 - q) > 1.0e-5){
    sq = sqrt(1.0 - q * q);
    psi = sqrt(q * q * (s * s + x * x) + y * y);
     
    *ax = (q / sq) * atan(sq * x / (psi + s));
    *ay = (q / sq) * atanh(sq * y / (psi + q * q * s));
  } else {
    psi = sqrt(s * s + x * x + y * y);
    *ax = x / (psi + s);
    *ay = y / (psi + s);
  }

  return;
}

double b_sie(double sig, double q)
{
  double ss;

  ss = sig / C_LIGHT_KMS;

  return facq_sie(q) * (4.0 * M_PI * ss * ss * dis_ls / dis_os) / ARCSEC2RADIAN;
}

double facq_sie(double q)
{
  static int ff;
  static double qq, fac;

  if((ff != 1) || (qq != q)){
    ff = 1;
    qq = q;
    
    /* my normalization */
    fac = 1.0 / sqrt(q);
    /* normalization in gravlens */
    /* fac = sqrt((1.0 + q * q) / (2.0 * q)) / sqrt(q); */
  }

  return fac;
}

/*--------------------------------------------------------------
  NFW potential
*/

void kapgam_nfwpot(double tx, double ty, double tx0, double ty0, double m, double c, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double bb, tt, a, b, pxx, pxy, pyy;
  double u[6];

  checkmodelpar_min(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(c, 0.0);

  calc_bbtt_nfw(m, c, &bb, &tt);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  u_calc((tx - tx0) / tt, (ty - ty0) / tt, e, si, co, u);
  
  a = bb * dphi_nfw_dl(u[0]);
  *ax = a * u[1] * tt;
  *ay = a * u[2] * tt;

  if(alponly != 1){
    b = bb * ddphi_nfw_dl(u[0]);
    
    pxx = b * u[1] * u[1] + a * u[3];
    pxy = b * u[1] * u[2] + a * u[4];
    pyy = b * u[2] * u[2] + a * u[5];
    
    *kap = 0.5 * (pxx + pyy);
    *gam1 = 0.5 * (pxx - pyy);
    *gam2 = pxy;
    if(alponly < 0) *phi = bb * phi_nfw_dl(u[0]) * tt * tt;
  }

  return;
} 

void calc_bbtt_nfw(double m, double c, double *bb, double *tt)
{
  double cc;

  if(nfw_users == 0){
    *bb = b_func(m, c);
    *tt = rtotheta(rs(m, c));
  } else {
    *tt = c;
    cc = rs(m, 1.0) / thetator(c);
    *bb = b_func(m, cc);
  }

  return;
}

double ddphi_nfw_dl(double x)
{
  return ddphi_kappadphi(kappa_nfw_dl(x), dphi_nfw_dl(x), x);
}

double ddphi_kappadphi(double kappa, double dphi, double x)
{
  return 2.0 * kappa - (dphi / x);
}

double kappa_nfw_dl(double x)
{
  if(x > (1.0 + 1.0e-6)){
    return 0.5 * (1.0 - 2.0 * atan(sqrt((x - 1.0) / (x + 1.0))) / sqrt(x * x - 1.0)) / (x * x - 1.0);
  } else if(x < (1.0 - 1.0e-6)){
    return 0.5 * (2.0 * atanh(sqrt((1.0 - x) / (x + 1.0))) / sqrt(1.0 - x * x) - 1.0) / (1.0 - x * x);
  } else {
    return 0.5 / 3.0;
  }
}

double dkappa_nfw_dl(double x)
{
   if((x > (1.0 + 1.0e-5)) || (x < (1.0 - 1.0e-5))){
     return 0.5 * (3.0 * x * x * func_hern_dl(x) - 2.0 * x * x - 1.0) / (x * (x * x - 1.0) * (x * x - 1.0));
   } else {
     return -0.2;
   }
}

double dphi_nfw_dl(double x)
{
  if(x > (1.0 + 1.0e-9)){
    return 2.0 * atan(sqrt((x - 1.0) / (x + 1.0))) / (x * sqrt(x * x - 1.0)) + log(0.5 * x) / x;
  } else if(x < (1.0 - 1.0e-9)){
    if(x > 1.0e-5){
      return 2.0 * atanh(sqrt((1.0 - x) / (x + 1.0))) / (x * sqrt(1.0 - x * x)) + log(0.5 * x) / x;
    } else {
      return 0.5 * x * log(2.0 / x);
    } 
  } else {
    return 1.0 + log(0.5);
  }
}

double phi_nfw_dl(double x)
{
  double a, b, bb;

  if(x > 1.0e-6){
    a = log(0.5 * x);
  
    if(x < (1.0 - 1.0e-9)){
      b = acosh(1.0 / x);
      bb = (-1.0) * b * b;
    } else if(x > (1.0 + 1.0e-9)){
      b = acos(1.0 / x);
      bb = b * b;
    } else {
      bb = 0.0;
    }

    return 0.5 * (a * a + bb); 
  } else {
    return 0.25 * x * x * log(2.0 / x);
  }
}

/* in units of Mpc/h */

double rs(double m, double c)
{
  static int ff;
  static double rr, mm, cc, dd;
  
  if((ff != 1) || (mm != m) || (cc != c) || (dd != delome)){
    ff = 1;
    mm = m;
    cc = c;
    dd = delome;
    rr = NFW_RS_NORM * pow(m / delome, 1.0 / 3.0) / c;
  }

  return rr;
}

double b_func(double m, double c)
{
  static int ff;
  static double bb, mm, cc, dd, ss, ll;

  if((ff != 1) || (mm != m) || (cc != c) || (dd != delome) || (ss != dis_os) || (ll != dis_ol)){
    ff = 1;
    mm = m;
    cc = c;
    dd = delome;
    ss = dis_os;
    ll = dis_ol;
    
    bb = (NFW_B_NORM) * dis_ol * dis_ls * pow(delome * delome * m, 1.0 / 3.0) * (c * c / hnfw(c)) / dis_os;
  }

  return bb;
}

double hnfw(double c)
{
  if(c < 1.0e-6){
    return 0.5 * c * c;
  } else {
    return log(1.0 + c) - (c / (1.0 + c));
  }
}

/*--------------------------------------------------------------
  NFW density
*/

void kapgam_nfw(double tx, double ty, double tx0, double ty0, double m, double c, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double q, bb, tt, uu, j0, j1;
  double bx, by, px, py, bpx, bpy, pxx, pxy, pyy, bpxx, bpxy, bpyy;

  checkmodelpar_min(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(c, 0.0);

  q = 1.0 - e;

  calc_bbtt_nfw(m, c, &bb, &tt);
  tt = tt / sqrt(q);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  bx = (co * (tx - tx0) - si * (ty - ty0)) / tt;
  by = (si * (tx - tx0) + co * (ty - ty0)) / tt;

  x_integ_sav = bx;
  y_integ_sav = by;
  q_integ_sav = q;

  uu = 1.0 / ell_xi2(1.0, 1.0);
  
  if(uu > 0.1){
    ell_u_lnint = 0;
  } else {
    ell_u_lnint = 1;
    ell_u_lnmin = (1.0e-4) * uu;
  }
  
  j1 = ell_integ_j(kappa_nfw_dl, 1);
  bpx = q * bx * j1;
  j0 = ell_integ_j(kappa_nfw_dl, 0);
  bpy = q * by * j0;

  ell_pxpy(bpx, bpy, si, co, &px, &py);

  *ax = bb * tt * px;
  *ay = bb * tt * py;

  if(alponly != 1){
   bpxx = 2.0 * q * bx * bx * ell_integ_k(dkappa_nfw_dl, 2) + q * j1;
   bpyy = 2.0 * q * by * by * ell_integ_k(dkappa_nfw_dl, 0) + q * j0;
   bpxy = 2.0 * q * bx * by * ell_integ_k(dkappa_nfw_dl, 1);

   ell_pxxpyy(bpxx, bpyy, bpxy, si, co, &pxx, &pyy, &pxy);

   *kap = 0.5 * bb * (pxx + pyy);
   *gam1 = 0.5 * bb * (pxx - pyy);
   *gam2 = bb * pxy;
   if(alponly < 0) *phi = 0.5 * q * bb * ell_integ_i(dphi_nfw_dl) * tt * tt;
  }

  return;
} 

/*--------------------------------------------------------------
  NFW density, approximation using CSE
*/

void kapgam_anfw(double tx, double ty, double tx0, double ty0, double m, double c, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double q, bb, tt, dx, dy, ddx, ddy, aax, aay;
  double pxx, pyy, pxy, rpxx, rpyy, rpxy;

  checkmodelpar_min(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(c, 0.0);

  q = 1.0 - e;

  calc_bbtt_nfw(m, c, &bb, &tt);
  tt = tt / sqrt(q);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * (pa - 90.0) * M_PI / 180.0);
    co = cos((-1.0) * (pa - 90.0) * M_PI / 180.0);
  }

  dx = (tx - tx0) / tt;
  dy = (ty - ty0) / tt;

  ddx = co * dx - si * dy;
  ddy = si * dx + co * dy;

  alpha_anfw_dl(ddx, ddy, q, &aax, &aay);

  *ax = bb * (aax * co + aay * si) * tt;
  *ay = bb * (aax * (-1.0) * si + aay * co) * tt;
  
  if(alponly != 1){
    ddphi_anfw_dl(ddx, ddy, q, &pxx, &pxy, &pyy);

    rpxx = co * co * pxx + 2.0 * co * si * pxy + si * si * pyy;
    rpyy = si * si * pxx - 2.0 * co * si * pxy + co * co * pyy;
    rpxy = si * co * (pyy - pxx) + (co * co - si * si) * pxy;

    *kap = 0.5 * bb * (rpxx + rpyy);
    *gam1 = 0.5 * bb * (rpxx - rpyy);
    *gam2 = bb * rpxy;

    if(alponly < 0){
      *phi = bb * phi_anfw_dl(ddx, ddy, q) * tt * tt;
    }
  }

  return;
} 

/*--------------------------------------------------------------
  generalized NFW potential
*/

void kapgam_gnfwpot(double tx, double ty, double tx0, double ty0, double m, double c, double alpha, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double bb, tt, a, b, pxx, pxy, pyy;
  double u[6];

  checkmodelpar_min(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(c, 0.0);
  checkmodelpar_mineq(alpha, 0.0);
  checkmodelpar_maxeq(alpha, 2.0);

  if(gnfw_usetab > 0) gnfw_maketable();

  alpha_gnfw_sav = alpha;
  calc_bbtt_gnfw(m, c, alpha, &bb, &tt);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  u_calc((tx - tx0) / tt, (ty - ty0) / tt, e, si, co, u);
  
  a = bb * dphi_gnfw_dl(u[0]);
  *ax = a * u[1] * tt;
  *ay = a * u[2] * tt;

  if(alponly != 1){
    b = bb * ddphi_gnfw_dl(u[0]);
    
    pxx = b * u[1] * u[1] + a * u[3];
    pxy = b * u[1] * u[2] + a * u[4];
    pyy = b * u[2] * u[2] + a * u[5];
    
    *kap = 0.5 * (pxx + pyy);
    *gam1 = 0.5 * (pxx - pyy);
    *gam2 = pxy;
    if(alponly < 0) *phi = bb * phi_gnfw_dl(u[0]) * tt * tt;
  }

  return;
} 

void calc_bbtt_gnfw(double m, double c, double alpha, double *bb, double *tt)
{
  double cc;

  if(nfw_users == 0){
    cc = (2.0 - alpha) * c;  /* c_-2 -> c_vir */
    *bb = b_func_gnfw(m, cc, alpha);
    *tt = rtotheta(rs(m, cc));
  } else {
    *tt = c;
    cc = rs(m, 1.0) / thetator(c);
    *bb = b_func_gnfw(m, cc, alpha);
  }

  return;
}

double ddphi_gnfw_dl(double x)
{
  return ddphi_kappadphi(kappa_gnfw_dl(x), dphi_gnfw_dl(x), x);
}

double kappa_gnfw_dl(double x)
{
  if(gnfw_usetab > 0){
    return kappa_gnfw_dl_tab(x, alpha_gnfw_sav);
  } else {
    x_gnfw_sav = x;
    return 0.5 * gsl_romberg1(kappa_gnfw_dl_func, -25.0, 25.0, TOL_ROMBERG_GNFW);
  }
}

double kappa_gnfw_dl_func(double x)
{
  double a, aa, xx, pp;
  
  xx = exp(x);
  a = sqrt(x_gnfw_sav * x_gnfw_sav + xx * xx);
  aa = 1.0 + a;
  pp = pow(a / aa, alpha_gnfw_sav);
  
  return xx / (aa * aa * aa * pp);
}
 
double dkappa_gnfw_dl(double x)
{
  if(gnfw_usetab > 0){
    return dkappa_gnfw_dl_tab(x, alpha_gnfw_sav);
  } else {
    x_gnfw_sav = x;
    return 0.5 * gsl_romberg1(dkappa_gnfw_dl_func, -25.0, 25.0, TOL_ROMBERG_GNFW);
  }
}

double dkappa_gnfw_dl_func(double x)
{
  double a, aa, xx, pp;
  
  xx = exp(x);
  a = sqrt(x_gnfw_sav * x_gnfw_sav + xx * xx);
  aa = 1.0 + a;
  pp = pow(a / aa, alpha_gnfw_sav);
  
  return x_gnfw_sav * xx * ((-3.0) * a - alpha_gnfw_sav) / (a * a * pp * aa * aa * aa * aa);
}
 
double dphi_gnfw_dl(double x)
{
  if(gnfw_usetab > 0){
    return dphi_gnfw_dl_tab(x, alpha_gnfw_sav);
  } else {
    double f;
    
    x_gnfw_sav = x;
    
    if(x > 1.0e-3){
      f = gsl_romberg1(dphi_gnfw_dl_func, smallcore, 1.0, TOL_ROMBERG_GNFW);
    } else {
      f = gsl_romberg1(dphi_gnfw_dl_funcln, -30.0, 0.0, TOL_ROMBERG_GNFW);
    }
    return (hgnfw(x) / x) + pow(x, 2.0 - alpha_gnfw_sav) * f;
  }
}

double dphi_gnfw_dl_func(double x)
{
  return pow(x + x_gnfw_sav, alpha_gnfw_sav - 3.0) * (1.0 - sqrt(1.0 - x * x)) / x;
}

double dphi_gnfw_dl_funcln(double lx)
{
  double x;
  
  x = exp(lx);

  return x * dphi_gnfw_dl_func(x);
}

double phi_gnfw_dl(double x)
{
  return gsl_qgaus(phi_gnfw_dl_func, -25.0, log(x));
}

double phi_gnfw_dl_func(double lx)
{
  double x;

  x = exp(lx);

  return x * dphi_gnfw_dl(x);
}

double hgnfw(double x)
{
  if(x < 1.0e3){
    return gsl_romberg1(hgnfw_func, 0.0, x, TOL_ROMBERG_GNFW);
  } else {
    return gsl_romberg1(hgnfw_funcln, -10.0, log(x), TOL_ROMBERG_GNFW);
  }
}

double hgnfw_func(double x)
{
  return pow(x / (1.0 + x), 2.0 - alpha_gnfw_sav) / (1.0 + x);
}

double hgnfw_funcln(double lx)
{
  double x;
  
  x = exp(lx);

  return x * hgnfw_func(x);
}

double b_func_gnfw(double m, double c, double alpha)
{
  static int ff;
  static double bb, mm, cc, dd, ss, ll, aa;

   alpha_gnfw_sav = alpha;

  if((ff != 1) || (mm != m) || (cc != c) || (dd != delome) || (ss != dis_os) || (ll != dis_ol) || (aa != alpha)){
    ff = 1;
    mm = m;
    cc = c;
    aa = alpha;
    dd = delome;
    ss = dis_os;
    ll = dis_ol;
    
    bb = (NFW_B_NORM) * dis_ol * dis_ls * pow(delome * delome * m, 1.0 / 3.0) * (c * c / hgnfw(c)) / dis_os;
  }

  return bb;
}

/*--------------------------------------------------------------
  generalized NFW density
*/

void kapgam_gnfw(double tx, double ty, double tx0, double ty0, double m, double c, double alpha, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double bb, tt, q, uu, j0, j1;
  double bx, by, px, py, bpx, bpy, pxx, pxy, pyy, bpxx, bpxy, bpyy;

  checkmodelpar_min(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(c, 0.0);
  checkmodelpar_mineq(alpha, 0.0);
  checkmodelpar_maxeq(alpha, 2.0);

  if(gnfw_usetab > 0) gnfw_maketable();

  alpha_gnfw_sav = alpha;
  q = 1.0 - e;
  calc_bbtt_gnfw(m, c, alpha, &bb, &tt);
  tt = tt / sqrt(q);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  bx = (co * (tx - tx0) - si * (ty - ty0)) / tt;
  by = (si * (tx - tx0) + co * (ty - ty0)) / tt;

  x_integ_sav = bx;
  y_integ_sav = by;
  q_integ_sav = q;

  uu = 1.0 / ell_xi2(1.0, 1.0);
  
  if(alpha <= 1.0){
    if(uu > 0.1){
      ell_u_lnint = 0;
    } else {
      ell_u_lnint = 1;
      ell_u_lnmin = (1.0e-4) * uu;
    }
  } else {
    ell_u_lnint = 1;
    ell_u_lnmin = smallcore * smallcore / ell_xi2(1.0, 1.0);
  }
  
  j1 = ell_integ_j(kappa_gnfw_dl, 1);
  bpx = q * bx * j1;
  j0 = ell_integ_j(kappa_gnfw_dl, 0);
  bpy = q * by * j0;

  ell_pxpy(bpx, bpy, si, co, &px, &py);

  *ax = bb * tt * px;
  *ay = bb * tt * py;

  if(alponly != 1){
   bpxx = 2.0 * q * bx * bx * ell_integ_k(dkappa_gnfw_dl, 2) + q * j1;
   bpyy = 2.0 * q * by * by * ell_integ_k(dkappa_gnfw_dl, 0) + q * j0;
   bpxy = 2.0 * q * bx * by * ell_integ_k(dkappa_gnfw_dl, 1);

   ell_pxxpyy(bpxx, bpyy, bpxy, si, co, &pxx, &pyy, &pxy);

   *kap = 0.5 * bb * (pxx + pyy);
   *gam1 = 0.5 * bb * (pxx - pyy);
   *gam2 = bb * pxy;
   if(alponly < 0) *phi = 0.5 * q * bb * ell_integ_i(dphi_gnfw_dl) * tt * tt;
  }

  return;
} 

/*--------------------------------------------------------------
  Hernquist potential
*/

void kapgam_hernpot(double tx, double ty, double tx0, double ty0, double m, double rb, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double bb, tt, a, b, pxx, pxy, pyy;
  double u[6];

  checkmodelpar_mineq(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(rb, 0.0);

  bb = b_func_hern(m, rb);
  tt = rb;

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  u_calc((tx - tx0) / tt, (ty - ty0) / tt, e, si, co, u);
  
  a = bb * dphi_hern_dl(u[0]);
  *ax = a * u[1] * tt;
  *ay = a * u[2] * tt;

  if(alponly != 1){
    b = bb * ddphi_hern_dl(u[0]);
    
    pxx = b * u[1] * u[1] + a * u[3];
    pxy = b * u[1] * u[2] + a * u[4];
    pyy = b * u[2] * u[2] + a * u[5];
    
    *kap = 0.5 * (pxx + pyy);
    *gam1 = 0.5 * (pxx - pyy);
    *gam2 = pxy;
    if(alponly < 0) *phi = bb * phi_hern_dl(u[0]) * tt * tt;
  }

  return;
} 

double b_func_hern(double m, double rb)
{
  double rr;

  rr = thetator(rb);
  return (m * inv_sigma_crit() / (rr * rr)) / (2.0 * M_PI);
}

double ddphi_hern_dl(double x)
{
  return ddphi_kappadphi(kappa_hern_dl(x), dphi_hern_dl(x), x);
}

double kappa_hern_dl(double x)
{
  if((x > (1.0 + 1.0e-5)) || (x < (1.0 - 1.0e-5))){
    return ((2.0 + x * x) * func_hern_dl(x) - 3.0) / ((x * x - 1.0) * (x * x - 1.0));
  } else {
    return 4.0 / 15.0;
  }
}

double dkappa_hern_dl(double x)
{
  if((x > (1.0 + 1.0e-4)) || (x < (1.0 - 1.0e-4))){
    return (2.0 + 13.0 * x * x - 3.0 * x * x * (x * x + 4.0) * func_hern_dl(x)) / (x * (x * x - 1.0) * (x * x - 1.0) * (x * x - 1.0)); 
  } else {
    return -0.4571429;
  }
}

double dphi_hern_dl(double x)
{
  if((x > (1.0 + 1.0e-7)) || (x < (1.0 - 1.0e-7))){
    return 2.0 * x * (1.0 - func_hern_dl(x)) / (x * x - 1.0); 
  } else {
    return 2.0 / 3.0;
  }
}

double phi_hern_dl(double x)
{
  if(x > 1.0e-3){
    return log(x * x / 4.0) + 2.0 * func_hern_dl(x); 
  } else {
    return x * x * (log(2.0 / x) - 0.5);
  }
}

double func_hern_dl(double x)
{
  double xx;
  if(x > (1.0 + 1.0e-9)){
    xx = sqrt(x * x - 1.0);
    return atan(xx) / xx;
  } else if(x < (1.0 - 1.0e-9)){
    if(x > 1.0e-5){
      xx = sqrt(1.0 - x * x);
      return atanh(xx) / xx;
    } else {
      return log(2.0 / x);
    }
  } else {
    return 1.0;
  }
}

/*--------------------------------------------------------------
  Hernquist density
*/

void kapgam_hern(double tx, double ty, double tx0, double ty0, double m, double rb, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double q, bb, tt, uu, j0, j1;
  double bx, by, px, py, bpx, bpy, pxx, pxy, pyy, bpxx, bpxy, bpyy;

  checkmodelpar_mineq(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(rb, 0.0);
  
  q = 1.0 - e;

  bb = b_func_hern(m, rb);
  tt = rb / sqrt(q);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  bx = (co * (tx - tx0) - si * (ty - ty0)) / tt;
  by = (si * (tx - tx0) + co * (ty - ty0)) / tt;

  x_integ_sav = bx;
  y_integ_sav = by;
  q_integ_sav = q;

  uu = 1.0 / ell_xi2(1.0, 1.0);
  
  if(uu > 0.1){
    ell_u_lnint = 0;
  } else {
    ell_u_lnint = 1;
    ell_u_lnmin = (1.0e-4) * uu;
  }
  
  j1 = ell_integ_j(kappa_hern_dl, 1);
  bpx = q * bx * j1;
  j0 = ell_integ_j(kappa_hern_dl, 0);
  bpy = q * by * j0;

  ell_pxpy(bpx, bpy, si, co, &px, &py);

  *ax = bb * tt * px;
  *ay = bb * tt * py;

  if(alponly != 1){
   bpxx = 2.0 * q * bx * bx * ell_integ_k(dkappa_hern_dl, 2) + q * j1;
   bpyy = 2.0 * q * by * by * ell_integ_k(dkappa_hern_dl, 0) + q * j0;
   bpxy = 2.0 * q * bx * by * ell_integ_k(dkappa_hern_dl, 1);

   ell_pxxpyy(bpxx, bpyy, bpxy, si, co, &pxx, &pyy, &pxy);

   *kap = 0.5 * bb * (pxx + pyy);
   *gam1 = 0.5 * bb * (pxx - pyy);
   *gam2 = bb * pxy;
   if(alponly < 0) *phi = 0.5 * q * bb * ell_integ_i(dphi_hern_dl) * tt * tt;
  }

  return;
} 

/*--------------------------------------------------------------
  Hernquist density, approximation using CSE
*/

void kapgam_ahern(double tx, double ty, double tx0, double ty0, double m, double rb, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double q, bb, tt, dx, dy, ddx, ddy, aax, aay;
  double pxx, pyy, pxy, rpxx, rpyy, rpxy;

  checkmodelpar_min(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(rb, 0.0);

  q = 1.0 - e;

  bb = b_func_hern(m, rb);
  tt = rb / sqrt(q);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * (pa - 90.0) * M_PI / 180.0);
    co = cos((-1.0) * (pa - 90.0) * M_PI / 180.0);
  }

  dx = (tx - tx0) / tt;
  dy = (ty - ty0) / tt;

  ddx = co * dx - si * dy;
  ddy = si * dx + co * dy;

  alpha_ahern_dl(ddx, ddy, q, &aax, &aay);

  *ax = bb * (aax * co + aay * si) * tt;
  *ay = bb * (aax * (-1.0) * si + aay * co) * tt;
  
  if(alponly != 1){
    ddphi_ahern_dl(ddx, ddy, q, &pxx, &pxy, &pyy);

    rpxx = co * co * pxx + 2.0 * co * si * pxy + si * si * pyy;
    rpyy = si * si * pxx - 2.0 * co * si * pxy + co * co * pyy;
    rpxy = si * co * (pyy - pxx) + (co * co - si * si) * pxy;

    *kap = 0.5 * bb * (rpxx + rpyy);
    *gam1 = 0.5 * bb * (rpxx - rpyy);
    *gam2 = bb * rpxy;

    if(alponly < 0){
      *phi = bb * phi_ahern_dl(ddx, ddy, q) * tt * tt;
    }
  }

  return;
} 

/*--------------------------------------------------------------
  power-law potential
*/

void kapgam_powpot(double tx, double ty, double tx0, double ty0, double zs_fid, double re, double gam, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double fac, a, b, pxx, pxy, pyy;
  double u[6];

  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(re, 0.0);
  checkmodelpar_min(gam, 1.0);
  checkmodelpar_max(gam, 3.0);

  fac = fac_pert(zs_fid);
  
  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }
  
  gam_pow_sav = gam;

  u_calc((tx - tx0) / re, (ty - ty0) / re, e, si, co, u);
  
  a = fac * dphi_pow_dl(u[0]);
  *ax = a * u[1] * re;
  *ay = a * u[2] * re;
  
  if(alponly != 1){
    b = fac * ddphi_pow_dl(u[0]);
    
    pxx = b * u[1] * u[1] + a * u[3];
    pxy = b * u[1] * u[2] + a * u[4];
    pyy = b * u[2] * u[2] + a * u[5];
    
    *kap = 0.5 * (pxx + pyy);
    *gam1 = 0.5 * (pxx - pyy);
    *gam2 = pxy;
    if(alponly < 0) *phi = fac * phi_pow_dl(u[0]) * re * re;
  }

  return;
} 

double ddphi_pow_dl(double x)
{
  return (2.0 - gam_pow_sav) * pow(x, 1.0 - gam_pow_sav);
}

double kappa_pow_dl(double x)
{
  return 0.5 * (3.0 - gam_pow_sav) * pow(x, 1.0 - gam_pow_sav);
}

double dkappa_pow_dl(double x)
{
  return 0.5 * (3.0 - gam_pow_sav) * (1.0 - gam_pow_sav) * pow(x, (-1.0) * gam_pow_sav);
}

double dphi_pow_dl(double x)
{
  return pow(x, 2.0 - gam_pow_sav);
}

double phi_pow_dl(double x)
{
  return pow(x, 3.0 - gam_pow_sav) / (3.0 - gam_pow_sav);
}

/*--------------------------------------------------------------
  power-law density
*/

void kapgam_pow(double tx, double ty, double tx0, double ty0, double zs_fid, double re, double gam, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  if(flag_pow_tm15 == 0){
    kapgam_pow_direct(tx, ty, tx0, ty0, zs_fid, re, gam, e, pa, kap, gam1, gam2, phi, ax, ay, alponly);
  } else {
    kapgam_pow_tm15(tx, ty, tx0, ty0, zs_fid, re, gam, e, pa, kap, gam1, gam2, phi, ax, ay, alponly);
  }
  return;
}

/*--------------------------------------------------------------
  power-law density, direct integration
*/

void kapgam_pow_direct(double tx, double ty, double tx0, double ty0, double zs_fid, double re, double gam, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double fac, q, tt, j0, j1;
  double bx, by, px, py, bpx, bpy, pxx, pxy, pyy, bpxx, bpxy, bpyy;

  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(re, 0.0);
  checkmodelpar_min(gam, 1.0);
  checkmodelpar_max(gam, 3.0);

  q = 1.0 - e;
  fac = fac_pert(zs_fid);
  
  gam_pow_sav = gam;
  tt = re / sqrt(q);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  bx = (co * (tx - tx0) - si * (ty - ty0)) / tt;
  by = (si * (tx - tx0) + co * (ty - ty0)) / tt;

  x_integ_sav = bx;
  y_integ_sav = by;
  q_integ_sav = q;

  if(gam < 1.7){
    ell_u_lnint = 0;
  } else {
    ell_u_lnint = 1;
    ell_u_lnmin = smallcore * smallcore / ell_xi2(1.0, 1.0);
  }
  
  j1 = ell_integ_j(kappa_pow_dl, 1);
  bpx = q * bx * j1;
  j0 = ell_integ_j(kappa_pow_dl, 0);
  bpy = q * by * j0;

  ell_pxpy(bpx, bpy, si, co, &px, &py);

  *ax = fac * tt * px;
  *ay = fac * tt * py;

  if(alponly != 1){
    /* if((gam <= 2.0) || ((bx * bx + by * by) > smallcore * smallcore)){ */
      bpxx = 2.0 * q * bx * bx * ell_integ_k(dkappa_pow_dl, 2) + q * j1;
      bpyy = 2.0 * q * by * by * ell_integ_k(dkappa_pow_dl, 0) + q * j0;
      bpxy = 2.0 * q * bx * by * ell_integ_k(dkappa_pow_dl, 1);
      
      ell_pxxpyy(bpxx, bpyy, bpxy, si, co, &pxx, &pyy, &pxy);
      
      *kap = 0.5 * fac * (pxx + pyy);
      *gam1 = 0.5 * fac * (pxx - pyy);
      *gam2 = fac * pxy;
      /* } else {
      *kap = fac * kappa_pow_dl(smallcore);
      *gam1 = fac * kappa_pow_dl(smallcore);
      *gam2 = fac * kappa_pow_dl(smallcore);
      } */
    if(alponly < 0) *phi = 0.5 * q * fac * ell_integ_i(dphi_pow_dl) * tt * tt;
  }

  return;
} 

/*--------------------------------------------------------------
  power-law density, using Tessore & Metcalf 2015
*/

void kapgam_pow_tm15(double tx, double ty, double tx0, double ty0, double zs_fid, double re, double gam, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double fac, q, tt;
  double ddx, ddy, r, psi, aa, aax, aay, dx, dy;
  double g1, g2, c1, s1, c2, s2;
  double ome[2];
  
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(re, 0.0);
  checkmodelpar_min(gam, 1.0);
  checkmodelpar_max(gam, 3.0);

  q = 1.0 - e;
  fac = fac_pert(zs_fid);
  tt = re * sqrt(q);
  dx = tx - tx0;
  dy = ty - ty0;
  
  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * (pa - 90.0) * M_PI / 180.0);
    co = cos((-1.0) * (pa - 90.0) * M_PI / 180.0);
  }

  ddx = co * dx - si * dy;
  ddy = si * dx + co * dy;

  r  = sqrt(q * q * ddx * ddx + ddy * ddy) + smallcore;
  psi = atan2(ddy, q * ddx);
  aa = 2.0 * tt * pow(tt / r, gam - 2.0) / (1.0 + q);
  pow_tm15_omega(psi, q, gam, ome);
 
  aax = aa * ome[0];
  aay = aa * ome[1];
  
  *ax = fac * (aax * co + aay * si);
  *ay = fac * (aax * (-1.0) * si + aay * co);
 
  if(alponly != 1){
    *kap = fac * 0.5 * (3.0 - gam) * pow(r / tt, 1.0 - gam);
    r = sqrt(ddx * ddx + ddy * ddy) + smallcore;
    c1 = ddx / r;
    s1 = ddy / r;
    c2 = c1 * c1 - s1 * s1;
    s2 = 2.0 * c1 * s1;
    g1 = (-1.0) * c2 * (*kap) + fac * (2.0 - gam) * (c1 * aax - s1 * aay) / r;
    g2 = (-1.0) * s2 * (*kap) + fac * (2.0 - gam) * (c1 * aay + s1 * aax) / r;
    c2 = co * co - si * si;
    s2 = 2.0 * co * si;
    *gam1 = c2 * g1 + s2 * g2;
    *gam2 = (-1.0) * s2 * g1 + c2 * g2;
    if(alponly < 0) *phi = fac * (ddx * aax + ddy * aay) / (3.0 - gam);
  }
  
  return;
}

void pow_tm15_omega(double psi, double q, double gam, double ome[])
{
  int k;
  double t, f, kk, fac, c1, s1, c2, s2;
  double a[2], b[2];
  
  t = 3.0 - gam;
  f = (1.0 - q) / (1.0 + q);
    
  c1 = cos(psi);
  s1 = sin(psi);
  c2 = c1 * c1 - s1 * s1;
  s2 = 2.0 * c1 * s1;

  a[0] = c1;
  a[1] = s1;
  
  ome[0] = a[0];
  ome[1] = a[1];

  k=0;
  do{
    k++;
    kk = (double)k;
    fac = (-1.0) * f * (2.0 * kk - t) / (2.0 * kk + t);

    b[0] = fac * (c2 * a[0] - s2 * a[1]);
    b[1] = fac * (s2 * a[0] + c2 * a[1]);
      
    a[0] = b[0];
    a[1] = b[1];

    ome[0] += a[0];
    ome[1] += a[1];

  } while((fabs(a[0]) > tol_pow_tm15) || (fabs(a[1]) > tol_pow_tm15));

  return;
}

/*--------------------------------------------------------------
  Sersic potential
*/

void kapgam_serspot(double tx, double ty, double tx0, double ty0, double m, double re, double n, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double bb, tt, a, b, pxx, pxy, pyy;
  double u[6];

  checkmodelpar_mineq(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(re, 0.0);
  checkmodelpar_mineq(n, SERSIC_N_MIN);
  checkmodelpar_maxeq(n, SERSIC_N_MAX);

  n_sers_sav = n;

  tt = re * bnn_sers(n);
  bb = b_func_sers(m, tt, n);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  u_calc((tx - tx0) / tt, (ty - ty0) / tt, e, si, co, u);
  
  a = bb * dphi_sers_dl(u[0]);
  *ax = a * u[1] * tt;
  *ay = a * u[2] * tt;

  if(alponly != 1){
    b = bb * ddphi_sers_dl(u[0]);
    
    pxx = b * u[1] * u[1] + a * u[3];
    pxy = b * u[1] * u[2] + a * u[4];
    pyy = b * u[2] * u[2] + a * u[5];
    
    *kap = 0.5 * (pxx + pyy);
    *gam1 = 0.5 * (pxx - pyy);
    *gam2 = pxy;
    if(alponly < 0) *phi = bb * phi_sers_dl(u[0]) * tt * tt;
  }

  return;
} 

double b_func_sers(double m, double rs, double n)
{
  double rr, gam;

  n_sers_sav = n;

  rr = thetator(rs);
  gam = gam2n1_sers(n);

  return (m * inv_sigma_crit() / (M_PI * rr * rr)) / gam;
}

double bn_sers(double n)
{
  double n2, n3, n4, bn;
  
  n2 = n * n;
  n3 = n2 * n;
  n4 = n3 * n;
  if(n > 0.36){
    /*  see Ciotti & Bertin astro-ph/9911078 */
    bn = 2.0 * n - (1.0 / 3.0) + (4.0 / 405.0) / n + (46.0 / 25515.0) / n2 + (131.0 / 1148175.0) / n3 - (2194697.0 / 30690717750.0) / n4;
  } else {
    /* see MacArthur et al. astro-ph/0208404 */
    bn = 0.01945 - 0.8902 * n + 10.95 * n2 - 19.67 * n3 + 13.43 * n4;
  }
  return bn;
}

double bnn_sers(double n)
{
  static int ff;
  static double nn, bn, bnn;

  n_sers_sav = n;

  if((ff != 1) || (nn != n)){
    ff = 1;
    nn = n;
    bn = bn_sers(n);
    bnn = pow(bn, (-1.0) * n);
  }

  return bnn;
}

double gam2n1_sers(double n)
{
  static int ff;
  static double nn, gam;

  if((ff != 1) || (nn != n)){
    ff = 1;
    nn = n;
    gam = exp(gsl_sf_lngamma(2.0 * n + 1.0));
  }

  return gam;
}

double ddphi_sers_dl(double x)
{
  return ddphi_kappadphi(kappa_sers_dl(x), dphi_sers_dl(x), x);
}

double kappa_sers_dl(double x)
{
  double xx;

  xx = pow(x, 1.0 / n_sers_sav);

  return exp((-1.0) * xx);
}

double dkappa_sers_dl(double x)
{
  double xx;

  xx = pow(x, 1.0 / n_sers_sav);

  return (-1.0) * xx * exp((-1.0) * xx) / (x * n_sers_sav);
}

double dphi_sers_dl(double x)
{
  double xx;

  xx = pow(x, 1.0 / n_sers_sav);

  return gam2n1_sers(n_sers_sav) * gsl_sf_gamma_inc_P(2.0 * n_sers_sav, xx) / x;
}

double phi_sers_dl(double x)
{
  return gsl_qgaus(phi_sers_dl_func, 1.0e-12, x);
}

double phi_sers_dl_func(double x)
{
  return dphi_sers_dl(x);
}

/*--------------------------------------------------------------
  Sersic density
*/

void kapgam_sers(double tx, double ty, double tx0, double ty0, double m, double re, double n, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double q, bb, tt, uu, j0, j1;
  double bx, by, px, py, bpx, bpy, pxx, pxy, pyy, bpxx, bpxy, bpyy;

  checkmodelpar_mineq(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(re, 0.0);
  checkmodelpar_mineq(n, SERSIC_N_MIN);
  checkmodelpar_maxeq(n, SERSIC_N_MAX);

  q = 1.0 - e;
  n_sers_sav = n;

  tt = re * bnn_sers(n);
  bb = b_func_sers(m, tt, n);
  tt = tt / sqrt(q);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  bx = (co * (tx - tx0) - si * (ty - ty0)) / tt;
  by = (si * (tx - tx0) + co * (ty - ty0)) / tt;

  x_integ_sav = bx;
  y_integ_sav = by;
  q_integ_sav = q;

  uu = 1.0 / ell_xi2(1.0, 1.0);
  
  if(uu > 0.1){
    ell_u_lnint = 0;
  } else {
    ell_u_lnint = 1;
    ell_u_lnmin = (1.0e-4) * uu;
  }
  
  j1 = ell_integ_j(kappa_sers_dl, 1);
  bpx = q * bx * j1;
  j0 = ell_integ_j(kappa_sers_dl, 0);
  bpy = q * by * j0;

  ell_pxpy(bpx, bpy, si, co, &px, &py);

  *ax = bb * tt * px;
  *ay = bb * tt * py;

  if(alponly != 1){
   bpxx = 2.0 * q * bx * bx * ell_integ_k(dkappa_sers_dl, 2) + q * j1;
   bpyy = 2.0 * q * by * by * ell_integ_k(dkappa_sers_dl, 0) + q * j0;
   bpxy = 2.0 * q * bx * by * ell_integ_k(dkappa_sers_dl, 1);

   ell_pxxpyy(bpxx, bpyy, bpxy, si, co, &pxx, &pyy, &pxy);

   *kap = 0.5 * bb * (pxx + pyy);
   *gam1 = 0.5 * bb * (pxx - pyy);
   *gam2 = bb * pxy;
   if(alponly < 0) *phi = 0.5 * q * bb * ell_integ_i(dphi_sers_dl) * tt * tt;
  }

  return;
} 

/*--------------------------------------------------------------
  truncated NFW potential
*/

void kapgam_tnfwpot(double tx, double ty, double tx0, double ty0, double m, double c, double t, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double bb, tt, cc, a, b, pxx, pxy, pyy;
  double u[6];

  checkmodelpar_min(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(c, 0.0);
  checkmodelpar_min(t, 0.0);

  calc_bbtt_tnfw(m, c, &bb, &tt, &cc);
  tnfw_set_tau(t * cc);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  u_calc((tx - tx0) / tt, (ty - ty0) / tt, e, si, co, u);
  
  a = bb * dphi_tnfw_dl(u[0]);
  *ax = a * u[1] * tt;
  *ay = a * u[2] * tt;

  if(alponly != 1){
    b = bb * ddphi_tnfw_dl(u[0]);
    
    pxx = b * u[1] * u[1] + a * u[3];
    pxy = b * u[1] * u[2] + a * u[4];
    pyy = b * u[2] * u[2] + a * u[5];
    
    *kap = 0.5 * (pxx + pyy);
    *gam1 = 0.5 * (pxx - pyy);
    *gam2 = pxy;
    if(alponly < 0) *phi = bb * phi_tnfw_dl(u[0]) * tt * tt;
  }

  return;
} 

void calc_bbtt_tnfw(double m, double c, double *bb, double *tt, double *cc)
{
  if(nfw_users == 0){
    *cc = c;
    *bb = b_func(m, c);
    *tt = rtotheta(rs(m, c));
  } else {
    *tt = c;
    *cc = rs(m, 1.0) / thetator(c);
    *bb = b_func(m, *cc);
  }

  return;
}

double ddphi_tnfw_dl(double x)
{
  return ddphi_kappadphi(kappa_tnfw_dl(x), dphi_tnfw_dl(x), x);
}

double kappa_tnfw_dl(double x)
{
  double f, t2, t4, u;
  double xa, ua, xb, ub, fa, fb;

  t2 = tau_tnfw_sav * tau_tnfw_sav;
  t4 = t2 * t2;
  u = x * x;

  if(t2 > 1.0e-4){
    if(u < ((1.0e5) * t2)){ 
      f = kappa_tnfw_dl_func1(x, u, t2, t4);
    } else {
      f = kappa_tnfw_dl_func2(x, u, t2, t4);
    }
  } else {
    if(u > (1.0e4)){
       f = kappa_tnfw_dl_func2(x, u, t2, t4);
    } else if(u < ((1.0e6) * t2)){
      f = kappa_tnfw_dl_func1(x, u, t2, t4);
    } else {
      /* derive values from interpolation ...  */
      xa = 1.0e2;
      ua = xa * xa;
      fa = kappa_tnfw_dl_func2(xa, ua, t2, t4);
      xb = (1.0e3) * tau_tnfw_sav;
      ub = xb * xb;
      fb = kappa_tnfw_dl_func1(xb, ub, t2, t4);
      f = exp((log(fa) - log(fb)) / (log(xa) - log(xb)) * (log(x) - log(xb)) + log(fb));
    }
  }

  return f;
}

double kappa_tnfw_dl_func1(double x, double u, double t2, double t4)
{
  double f, fx, st2u, st1;
  
  if((x > (1.0 - 1.0e-7)) && (x < (1.0 + 1.0e-7))){ 
    st1 = sqrt(1.0 + t2);
    f = (tau_tnfw_sav / (24.0 * (1.0 + t2) * (1.0 + t2) * (1.0 + t2) * (1.0 + t2))) * (tau_tnfw_sav * (-6.0 + 62.0 * t4 + 4.0 * t4 * t2 + t2 * (52.0 - 30.0 * M_PI * st1)) + 3.0 * st1 * (1.0 + 7.0 * t2 - 4.0 * t4) * log(1.0 + 2.0 * tau_tnfw_sav * (tau_tnfw_sav + st1)));
  } else {
    st2u = sqrt(t2 + u);
    fx = func_hern_dl(x);
    f = (t4 / (4.0 * (t2 + 1.0) * (t2 + 1.0) * (t2 + 1.0))) * (2.0 * (t2 + 1.0) * (1.0 - fx) / (u - 1.0) + 8.0 * fx + (t4 - 1.0) / (t2 * (t2 + u)) - M_PI * (4.0 * (t2 + u) + t2 + 1.0) / (st2u * (t2 + u)) + func_tnfw_dl(x) * (t2 * (t4 - 1.0) + (t2 + u) * (3.0 * t4 - 6.0 * t2 - 1.0)) / (tau_tnfw_sav * t2 * st2u * (t2 + u)));
  }

  return f;
}

double kappa_tnfw_dl_func2(double x, double u, double t2, double t4)
{
  double f, u3;
  
  u3 = u * u * u;
  f = 4.0 * t4 / (15.0 * u3) - 5.0 * M_PI * t4 / (32.0 * u3 * x) - 8.0 * (t4 * (-3.0 + 2.0 * t2)) / (35.0 * u3 * u) + 35.0 * t4 * M_PI * (t2 - 1.0) / (128.0 * u3 * u * x);

  return f;
}

double dkappa_tnfw_dl(double x)
{
  double f, t2, t4, u;
  double xa, ua, xb, ub, fa, fb;
    
  t2 = tau_tnfw_sav * tau_tnfw_sav;
  t4 = t2 * t2;
  u = x * x;
  
  if(t2 > 1.0e-4){
    if(u < ((1.0e5) * t2)){ 
      f = dkappa_tnfw_dl_func1(x, u, t2, t4);
    } else {
      f = dkappa_tnfw_dl_func2(x, u, t2, t4);
    }
  } else {
    if(u > (1.0e4)){
       f = dkappa_tnfw_dl_func2(x, u, t2, t4);
    } else if(u < ((1.0e6) * t2)){
      f = dkappa_tnfw_dl_func1(x, u, t2, t4);
    } else {
      /* derive values from interpolation ... */
      xa = 1.0e2;
      ua = xa * xa;
      fa = dkappa_tnfw_dl_func2(xa, ua, t2, t4);
      xb = (1.0e3) * tau_tnfw_sav;
      ub = xb * xb;
      fb = dkappa_tnfw_dl_func1(xb, ub, t2, t4);
      f = exp((log(fa) - log(fb)) / (log(xa) - log(xb)) * (log(x) - log(xb)) + log(fb));
    }
  }

  return f;
}

double dkappa_tnfw_dl_func1(double x, double u, double t2, double t4)
{
  double f, fx, st2u, t2u2, t2u52, st1;

  if((x > (1.0 - 1.0e-6)) && (x < (1.0 + 1.0e-6))){ 
    st1 = sqrt(1.0 + t2);
    f = ((-1.0) * tau_tnfw_sav / (120.0 * (1.0 + t2) * (1.0 + t2) * (1.0 + t2) * (1.0 + t2) * (1.0 + t2))) * (2.0 * tau_tnfw_sav * (-15.0 + 271.0 * t4 + 56.0 * t4 * t2 + 12.0 * t4 * t4 + t2 * (212.0 - 105.0 * M_PI * st1)) - 15.0 * st1 * (-1.0 - 9.0 * t2 + 6.0 * t4) * log(1.0 + 2.0 * tau_tnfw_sav * (tau_tnfw_sav + st1)));
  } else {
    st2u = sqrt(t2 + u);
    t2u2 = (t2 + u) * (t2 + u);
    t2u52 = t2u2 * st2u;
    fx = func_hern_dl(x);

    f = (t4 / (4.0 * (t2 + 1.0) * (t2 + 1.0) * (t2 + 1.0))) * (2.0 * (t2 + 1.0) * (3.0 * u * fx - 2.0 * u - 1.0) / (x * (u - 1.0) * (u - 1.0)) + 8.0 * (1.0 - u * fx) / (x * (u - 1.0)) - 2.0 * x * (t4 - 1.0) / (t2 * t2u2) + M_PI * x * (4.0 * (t2 + u) + 3.0 * (t2 + 1.0)) / t2u52 + (t2 * (t4 - 1.0) + (t2 + u) * (3.0 * t4 - 6.0 * t2 - 1.0)) / (x * t2 * t2u2) - func_tnfw_dl(x) * x * (3.0 * t2 * (t4 - 1.0) + (t2 + u) * (3.0 * t4 - 6.0 * t2 - 1.0)) / (tau_tnfw_sav * t2 * t2u52));
  }

  return f;
}

double dkappa_tnfw_dl_func2(double x, double u, double t2, double t4)
{
  double f, u3;
  
  u3 = u * u * u;
  f = (-8.0 * t4) / (5.0 * u3 * x) + 35.0 * M_PI * t4 / (32.0 * u3 * u) + 64.0 * (t4 * (-3.0 + 2.0 * t2)) / (35.0 * u3 * u * x) - 315.0 * t4 * M_PI * (t2 - 1.0) / (128.0 * u3 * u * u);

  return f;
}

double dphi_tnfw_dl(double x)
{
  double f, t2, t4, u, st2u, fx, lx;
  double lnt, ln2;

  t2 = tau_tnfw_sav * tau_tnfw_sav;
  t4 = t2 * t2;
  u = x * x;
  lnt = log(tau_tnfw_sav);
  
  if(((x > 3.0e-3) && (tau_tnfw_sav > 3.0e-3)) || (u > (t2 * 1.0e-3))){
    st2u = sqrt(t2 + u);
    fx = func_hern_dl(x);
    lx = func_tnfw_dl(x);

    f = (t4 / (2.0 * (t2 + 1.0) * (t2 + 1.0) * (t2 + 1.0) * x)) * (2.0 * (t2 + 4.0 * u - 3.0) * fx + (M_PI * (3.0 * t2 - 1.0) + 2.0 * tau_tnfw_sav * (t2 - 3.0) * lnt) / tau_tnfw_sav + (-t2 * tau_tnfw_sav * M_PI * (4.0 * u + 3.0 * t2 - 1.0) + (2.0 * t4 * t2 - 6.0 * t4 + u * (3.0 * t4 - 6.0 * t2 - 1.0)) * lx) / (t2 * tau_tnfw_sav * st2u));
  } else {
    ln2 = log(2.0);

    if(tau_tnfw_sav > 3.0e-3){
      f = (1.0 / (4.0 * ((t2 + 1.0) * (t2 + 1.0) * (t2 + 1.0)))) * (3.0 * t2 + (-1.0) * M_PI * tau_tnfw_sav - 5.0 * M_PI * tau_tnfw_sav * t2 + t4 * t2 * (-1.0 + 2.0 * ln2) + t4 * (2.0 + 10.0 * ln2) + (2.0 + 6.0 * t2 - 4.0 * t4) * (ln2 + lnt) - 2.0 * (t2 + 1.0) * (t2 + 1.0) * (t2 + 1.0) * log(x)) * x;
    } else {
      f = (x * (2.0 * ln2 + 2.0 * lnt - 2.0 * log(x) - (M_PI * tau_tnfw_sav)) / 4.0) + 3.0 * x * t2 / 4.0;
    }
  }

  return f;
}

double phi_tnfw_dl(double x)
{
  double f, t2, t4, u, st2u, fx, lx;
  double ln2, lnt;

  t2 = tau_tnfw_sav * tau_tnfw_sav;
  t4 = t2 * t2;
  u = x * x;
  ln2 = log(2.0);
  lnt = log(tau_tnfw_sav);

  if(((x > 3.0e-3) && (tau_tnfw_sav > 3.0e-3)) || (u > (t2 * 1.0e-3))){
    st2u = sqrt(t2 + u);
    fx = func_hern_dl(x);
    lx = func_tnfw_dl(x);
    
    f = (1.0 / (4.0 * (t2 + 1.0) * (t2 + 1.0) * (t2 + 1.0))) * (2.0 * t2 * tau_tnfw_sav * M_PI * (-4.0 * tau_tnfw_sav * st2u + (3.0 * t2 - 1.0) * log(st2u + tau_tnfw_sav)) + 2.0 * tau_tnfw_sav * (3.0 * t4 - 6.0 * t2 - 1.0) * st2u * lx + 2.0 * t4 * (t2 - 3.0) * lx * lx + 16.0 * t4 * (u - 1.0) * fx + 2.0 * t4 * (t2 - 3.0) * (u - 1.0) * fx * fx + t2 * (2.0 * t2 * (t2 - 3.0) * lnt - 3.0 * t4 - 2.0 * t2 + 1.0) * log(u) + 2.0 * t2 * (t2 * (4.0 * tau_tnfw_sav * M_PI + (t2 - 3.0) * ln2 * ln2 + 8.0 * ln2) - (ln2 + lnt) * (1.0 + 6.0 * t2 - 3.0 * t4 + t2 * (t2 - 3.0) * (ln2 + lnt)) - tau_tnfw_sav * M_PI * (3.0 * t2 - 1.0) * (ln2 + lnt)));
  } else {
    if(tau_tnfw_sav > 3.0e-3){
      f = (1.0 / (8.0 * ((t2 + 1.0) * (t2 + 1.0) * (t2 + 1.0)))) * (1.0 + tau_tnfw_sav * ((-1.0) * M_PI + tau_tnfw_sav * (6.0 + tau_tnfw_sav * ((-5.0) * M_PI + tau_tnfw_sav * (5.0 + (5.0 + t2) * 2.0 * ln2)))) + (2.0 + 6.0 * t2 - 4.0 * t4) * (ln2 + lnt) - 2.0 * (t2 + 1.0) * (t2 + 1.0) * (t2 + 1.0) * log(x)) * u;
    } else {
      f = u * (1.0 + 2.0 * ln2 + 2.0 * lnt - 2.0 * log(x) - (M_PI * tau_tnfw_sav)) / 8.0;
    }
  }

  return f;
}

double tnfw_mtot(double c, double t)
{
  double f, t2, g;

  t2 = t * t + 1.0;

  f = t * t / (2.0 * hnfw(c) * t2 * t2 * t2);

  g = 2.0 * t * t * (t * t - 3.0) * log(t) - (3.0 * t * t - 1.0) * (t * t + 1.0 - M_PI * t);
  
  return f * g;
}

double func_tnfw_dl(double x)
{
  double xx;

  xx = x / (sqrt(tau_tnfw_sav * tau_tnfw_sav + x * x) + tau_tnfw_sav);

  return log(xx);
}

void tnfw_set_tau(double tau)
{
  tau_tnfw_sav = tau;

  return;
}

/*--------------------------------------------------------------
  truncated NFW density
*/

void kapgam_tnfw(double tx, double ty, double tx0, double ty0, double m, double c, double t, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double q, bb, tt, cc, uu, j0, j1;
  double bx, by, px, py, bpx, bpy, pxx, pxy, pyy, bpxx, bpxy, bpyy;

  checkmodelpar_min(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(c, 0.0);
  checkmodelpar_min(t, 0.0);

  q = 1.0 - e;

  calc_bbtt_tnfw(m, c, &bb, &tt, &cc);
  tt = tt / sqrt(q);
  tnfw_set_tau(t * cc);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  bx = (co * (tx - tx0) - si * (ty - ty0)) / tt;
  by = (si * (tx - tx0) + co * (ty - ty0)) / tt;

  x_integ_sav = bx;
  y_integ_sav = by;
  q_integ_sav = q;

  uu = 1.0 / ell_xi2(1.0, 1.0);
  
  if(uu > 0.1){
    ell_u_lnint = 0;
  } else {
    ell_u_lnint = 1;
    ell_u_lnmin = (1.0e-4) * uu;
  }
  
  j1 = ell_integ_j(kappa_tnfw_dl, 1);
  bpx = q * bx * j1;
  j0 = ell_integ_j(kappa_tnfw_dl, 0);
  bpy = q * by * j0;

  ell_pxpy(bpx, bpy, si, co, &px, &py);

  *ax = bb * tt * px;
  *ay = bb * tt * py;

  if(alponly != 1){
   bpxx = 2.0 * q * bx * bx * ell_integ_k(dkappa_tnfw_dl, 2) + q * j1;
   bpyy = 2.0 * q * by * by * ell_integ_k(dkappa_tnfw_dl, 0) + q * j0;
   bpxy = 2.0 * q * bx * by * ell_integ_k(dkappa_tnfw_dl, 1);

   ell_pxxpyy(bpxx, bpyy, bpxy, si, co, &pxx, &pyy, &pxy);

   *kap = 0.5 * bb * (pxx + pyy);
   *gam1 = 0.5 * bb * (pxx - pyy);
   *gam2 = bb * pxy;
   if(alponly < 0) *phi = 0.5 * q * bb * ell_integ_i(dphi_tnfw_dl) * tt * tt;
  }

  return;
} 

/*--------------------------------------------------------------
  Einasto potential
*/

void kapgam_einpot(double tx, double ty, double tx0, double ty0, double m, double c, double alpha, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double bb, tt, a, b, pxx, pxy, pyy;
  double u[6];

  checkmodelpar_min(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(c, 0.0);
  checkmodelpar_mineq(alpha, 0.02);
  checkmodelpar_maxeq(alpha, 1.0);

  if(ein_usetab > 0) ein_maketable();

  alpha_ein_sav = alpha;

  calc_bbtt_ein(m, c, alpha, &bb, &tt);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  u_calc((tx - tx0) / tt, (ty - ty0) / tt, e, si, co, u);
  
  a = bb * dphi_ein_dl(u[0]);
  *ax = a * u[1] * tt;
  *ay = a * u[2] * tt;

  if(alponly != 1){
    b = bb * ddphi_ein_dl(u[0]);
    
    pxx = b * u[1] * u[1] + a * u[3];
    pxy = b * u[1] * u[2] + a * u[4];
    pyy = b * u[2] * u[2] + a * u[5];
    
    *kap = 0.5 * (pxx + pyy);
    *gam1 = 0.5 * (pxx - pyy);
    *gam2 = pxy;
    if(alponly < 0) *phi = bb * phi_ein_dl(u[0]) * tt * tt;
  }

  return;
} 

void calc_bbtt_ein(double m, double c, double alpha, double *bb, double *tt)
{
  double cc;

  if(nfw_users == 0){
    *bb = b_func_ein(m, c, alpha);
    *tt = rtotheta(rs(m, c));
  } else {
    *tt = c;
    cc = rs(m, 1.0) / thetator(c);
    *bb = b_func_ein(m, cc, alpha);
  }

  return;
}

double ddphi_ein_dl(double x)
{
  return ddphi_kappadphi(kappa_ein_dl(x), dphi_ein_dl(x), x);
}

double kappa_ein_dl(double x)
{
  if(ein_usetab > 0){
    return kappa_ein_dl_tab(x, alpha_ein_sav);
  } else { 
    x_ein_sav = x;
    return 0.5 * gsl_romberg1(kappa_ein_dl_func, -25.0, 25.0, TOL_ROMBERG_EIN);
  }
}

double kappa_ein_dl_func(double x)
{
  double xx, yy, ff;
  
  xx = exp(x);
  yy = (x_ein_sav * x_ein_sav + xx * xx);
  ff = (-2.0) * pow(yy, 0.5 * alpha_ein_sav) / alpha_ein_sav;

  if(ff > (-300.0)){
    return exp(ff) * xx;
  } else {
    return 0.0;
  }
}
 
double dkappa_ein_dl(double x)
{
   if(ein_usetab > 0){
     return dkappa_ein_dl_tab(x, alpha_ein_sav);
   } else { 
     x_ein_sav = x;
     return 0.5 * gsl_romberg1(dkappa_ein_dl_func, -25.0, 25.0, TOL_ROMBERG_EIN);
   }
}

double dkappa_ein_dl_func(double x)
{
  double xx, yy, pp, ff;
  
  xx = exp(x);
  yy = (x_ein_sav * x_ein_sav + xx * xx);
  pp = pow(yy, 0.5 * alpha_ein_sav);
  ff = (-2.0) * pp / alpha_ein_sav;

  if(ff > (-300.0)){
    return exp(ff) * xx * (-2.0) * (pp / yy) * x_ein_sav;
  } else {
    return 0.0;
  }
}
 
double dphi_ein_dl(double x)
{
  if(ein_usetab > 0){
    return dphi_ein_dl_tab(x, alpha_ein_sav);
  } else { 
    double f;
    
    x_ein_sav = x;
    
    /*  f = gsl_romberg1(dphi_ein_dl_func, smallcore, 1.0, TOL_ROMBERG_EIN); */
     f = gsl_romberg1(dphi_ein_dl_funcln, -30.0, 0.0, TOL_ROMBERG_EIN);
    
    return f / x;
  }
}

double dphi_ein_dl_func(double x)
{
  double xx;

  xx = sqrt(1.0 + x * x);

  return (hein(x_ein_sav * xx) + x * hein(x_ein_sav * xx / x)) / ((1.0 + x * x) * xx);
}

double dphi_ein_dl_funcln(double lx)
{
  double x;
  
  x = exp(lx);

  return x * dphi_ein_dl_func(x);
}

double phi_ein_dl(double x)
{
  return gsl_qgaus(phi_ein_dl_func, -25.0, log(x));
}

double phi_ein_dl_func(double lx)
{
  double x;

  x = exp(lx);

  return x * dphi_ein_dl(x);
}

double hein(double x)
{
  double gam;

  gam = gsl_sf_gamma_inc_P(3.0 / alpha_ein_sav, 2.0 * pow(x, alpha_ein_sav) / alpha_ein_sav) * exp(gsl_sf_lngamma(3.0 / alpha_ein_sav));
  
  return pow(0.5 * alpha_ein_sav, 3.0 / alpha_ein_sav) * gam / alpha_ein_sav;

}

double b_func_ein(double m, double c, double alpha)
{
  static int ff;
  static double bb, mm, cc, dd, ss, ll, aa;

   alpha_ein_sav = alpha;

  if((ff != 1) || (mm != m) || (cc != c) || (dd != delome) || (ss != dis_os) || (ll != dis_ol) || (aa != alpha)){
    ff = 1;
    mm = m;
    cc = c;
    aa = alpha;
    dd = delome;
    ss = dis_os;
    ll = dis_ol;
    
    bb = (NFW_B_NORM) * dis_ol * dis_ls * pow(delome * delome * m, 1.0 / 3.0) * (c * c / hein(c)) / dis_os;
  }

  return bb;
}

/*--------------------------------------------------------------
  Einasto density
*/

void kapgam_ein(double tx, double ty, double tx0, double ty0, double m, double c, double alpha, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double bb, tt, q, uu, j0, j1;
  double bx, by, px, py, bpx, bpy, pxx, pxy, pyy, bpxx, bpxy, bpyy;

  checkmodelpar_min(m, 0.0);
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(c, 0.0);
  checkmodelpar_mineq(alpha, 0.02);
  checkmodelpar_maxeq(alpha, 1.0);

  if(ein_usetab > 0) ein_maketable();

  alpha_ein_sav = alpha;
  q = 1.0 - e;

  calc_bbtt_ein(m, c, alpha, &bb, &tt);
  tt = tt / sqrt(q);

  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  bx = (co * (tx - tx0) - si * (ty - ty0)) / tt;
  by = (si * (tx - tx0) + co * (ty - ty0)) / tt;

  x_integ_sav = bx;
  y_integ_sav = by;
  q_integ_sav = q;

  uu = 1.0 / ell_xi2(1.0, 1.0);
  
  if(alpha <= 1.0){
    if(uu > 0.1){
      ell_u_lnint = 0;
    } else {
      ell_u_lnint = 1;
      ell_u_lnmin = (1.0e-4) * uu;
    }
  } else {
    ell_u_lnint = 1;
    ell_u_lnmin = smallcore * smallcore / ell_xi2(1.0, 1.0);
  }
  
  j1 = ell_integ_j(kappa_ein_dl, 1);
  bpx = q * bx * j1;
  j0 = ell_integ_j(kappa_ein_dl, 0);
  bpy = q * by * j0;

  ell_pxpy(bpx, bpy, si, co, &px, &py);

  *ax = bb * tt * px;
  *ay = bb * tt * py;

  if(alponly != 1){
   bpxx = 2.0 * q * bx * bx * ell_integ_k(dkappa_ein_dl, 2) + q * j1;
   bpyy = 2.0 * q * by * by * ell_integ_k(dkappa_ein_dl, 0) + q * j0;
   bpxy = 2.0 * q * bx * by * ell_integ_k(dkappa_ein_dl, 1);

   ell_pxxpyy(bpxx, bpyy, bpxy, si, co, &pxx, &pyy, &pxy);

   *kap = 0.5 * bb * (pxx + pyy);
   *gam1 = 0.5 * bb * (pxx - pyy);
   *gam2 = bb * pxy;
   if(alponly < 0) *phi = 0.5 * q * bb * ell_integ_i(dphi_ein_dl) * tt * tt;
  }

  return;
} 

/*--------------------------------------------------------------
  Gaussian potential
*/

void kapgam_gaupot(double tx, double ty, double tx0, double ty0, double zs_fid, double sig, double kap0, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly)
{
  static int ff;
  static double papa, si, co;
  double fac, a, b, pxx, pxy, pyy;
  double u[6];

  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(sig, 0.0);

  fac = fac_pert(zs_fid) * kap0;
  
  if((ff != 1) || (papa != pa)){
    ff = 1;
    papa = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }
  
  u_calc((tx - tx0) / sig, (ty - ty0) / sig, e, si, co, u);
  
  a = fac * dphi_gau_dl(u[0]);
  *ax = a * u[1] * sig;
  *ay = a * u[2] * sig;
  
  if(alponly != 1){
    b = fac * ddphi_gau_dl(u[0]);
    
    pxx = b * u[1] * u[1] + a * u[3];
    pxy = b * u[1] * u[2] + a * u[4];
    pyy = b * u[2] * u[2] + a * u[5];
    
    *kap = 0.5 * (pxx + pyy);
    *gam1 = 0.5 * (pxx - pyy);
    *gam2 = pxy;
    if(alponly < 0) *phi = fac * phi_gau_dl(u[0]) * sig * sig;
  }

  return;
} 

double ddphi_gau_dl(double x)
{
  if(x < 1.0e-4){
    return 1.0 - 3.0 * x * x / 4.0;
  } else {
    return 2.0 * exp((-0.5) * x * x) - (dphi_gau_dl(x) / x);
  }
}

double dphi_gau_dl(double x)
{
  if(x < 1.0e-4){
    return x - (x * x * x) / 4.0;
  } else {
    return 2.0 * (1.0 - exp((-0.5) * x * x)) / x;
  }
}

double phi_gau_dl(double x)
{
  if(x < 9.0e-5){
    return log(2.0) - 0.57721566490153286 + 0.5 * x * x;
  } else if(x < 1.0e1){
    return 2.0 * log(x) - gsl_sf_expint_Ei((-0.5) * x * x);
  } else {
    return 2.0 * log(x);
  }
}

/*--------------------------------------------------------------
  for elliptical potential
*/

void u_calc(double dx, double dy, double e, double si, double co, double u[6])
{
  double ep, si2, co2, ddx, ddy;

  /* ep = e / 2.27; */
  ep = e;

  si2 = 2.0 * si * co;
  co2 = co * co - si * si;

  ddx = co * dx - si * dy;
  ddy = si * dx + co * dy;

  /* u */
  u[0] = sqrt((1.0 + ep) * ddx * ddx + (1.0 - ep) * ddy * ddy)
     + smallcore; /* avoid error at the center */
  /* u_x */
  u[1] = (dx + ep * (dx * co2 - dy * si2)) / u[0];
  /* u_y */
  u[2] = (dy - ep * (dy * co2 + dx * si2)) / u[0];
  /* u_xx */
  u[3] = (1.0 + ep * co2 - u[1] * u[1]) / u[0];
  /* u_xy */
  u[4] = (( - 1.0) * ep * si2 - u[1] * u[2]) / u[0];
  /* u_yy */
  u[5] = (1.0 - ep * co2 - u[2] * u[2]) / u[0];
  
  return;
}

/*--------------------------------------------------------------
  for elliptical dentity
*/

double ell_integ_k(double (*func)(double), int n)
{
  func_sav = func;
  n_integ_sav = n;

  if(ell_u_lnint == 0){
    return gsl_romberg2(ell_integ_k_func, ULIM_JHK, 1.0, TOL_ROMBERG_JHK);
  } else {
    return gsl_romberg2(ell_integ_k_funcln, log(ell_u_lnmin), 0.0, TOL_ROMBERG_JHK);
  }
}

double ell_integ_k_func(double u)
{
  double se, equ, f;

  equ = ell_qu(q_integ_sav, u);
  se = sqrt(ell_xi2(u, equ));
  f = 2.0 * se * ell_nhalf(equ, n_integ_sav);

  return u * func_sav(se) / f;
}

double ell_integ_k_funcln(double lu)
{
  double u;

  u = exp(lu);

  return u * ell_integ_k_func(u);
}

double ell_integ_j(double (*func)(double), int n)
{
  func_sav = func;
  n_integ_sav = n;

  if(ell_u_lnint == 0){
    return gsl_romberg2(ell_integ_j_func, ULIM_JHK, 1.0, TOL_ROMBERG_JHK);
  } else {
    return gsl_romberg2(ell_integ_j_funcln, log(ell_u_lnmin), 0.0, TOL_ROMBERG_JHK);
  }
}

double ell_integ_j_func(double u)
{
  double se, equ, f;

  equ = ell_qu(q_integ_sav, u);
  se = sqrt(ell_xi2(u, equ));
  f = ell_nhalf(equ, n_integ_sav);

  return func_sav(se) / f;
}

double ell_integ_j_funcln(double lu)
{
  double u;

  u = exp(lu);

  return u * ell_integ_j_func(u);
}

double ell_integ_i(double (*func)(double))
{
  func_sav = func;

  if(ell_u_lnint == 0){
    return gsl_romberg2(ell_integ_i_func, ULIM_JHK, 1.0, TOL_ROMBERG_JHK);
  } else {
    return gsl_romberg2(ell_integ_i_funcln, log(ell_u_lnmin), 0.0, TOL_ROMBERG_JHK);
  }
}

double ell_integ_i_func(double u)
{
  double se, equ, f;

  equ = ell_qu(q_integ_sav, u);
  se = sqrt(ell_xi2(u, equ));
  f = u * ell_nhalf(equ, 0);

  return se * func_sav(se) / f;
}

double ell_integ_i_funcln(double lu)
{
  double u;

  u = exp(lu);

  return u * ell_integ_i_func(u);
}

double ell_xi2(double u, double equ)
{
  return u * (y_integ_sav * y_integ_sav + x_integ_sav * x_integ_sav / equ + smallcore * smallcore);
}

double ell_qu(double q, double u)
{
  return 1.0 - (1.0 - q * q) * u;
}

double ell_nhalf(double x, int n)
{
  int i;
  double f;

  f = sqrt(x);
  for(i=0;i<n;i++) f = f * x;

  return f;
}

void ell_pxpy(double bpx, double bpy, double si, double co, double *px, double *py)
{
  *px =          bpx * co + bpy * si;
  *py = (-1.0) * bpx * si + bpy * co;

  return;
}

void ell_pxxpyy(double bpxx, double bpyy, double bpxy, double si, double co, double *pxx, double *pyy, double *pxy)
{
  *pxx = bpxx * co * co          + bpxy * 2.0 * co * si       + bpyy * si * si;
  *pxy = bpxx * (-1.0) * si * co + bpxy * (co * co - si * si) + bpyy * co * si;
  *pyy = bpxx * si * si          + bpxy * (-2.0) * co * si    + bpyy * co * co;
}

#undef NM
