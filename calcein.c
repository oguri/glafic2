#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

static double a_sav, b_sav, s_sav;
static double x0_calkap_sav, y0_calkap_sav, r_calkap_sav;
static int  lensid_calkap_sav;

/*--------------------------------------------------------------
  calculate Einstein radius of individual lens ( e -> 0 limit ) 
*/

void calcein(double zs)
{
  int i;
  double ein, mein;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;
  
  sprintf(fname, "%s_ein.dat", file_prefix);

  fprintf(stderr, "######## calculating Einstein radii (e=0 limit)\n");
  fprintf(stderr, " zs = %e\n", zs);
  fprintf(stderr, " output file name = %s \n\n", fname);

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (calcein)");

  set_distance_lpl_zs(zs);
  
  for(i=0;i<num_len;i++){
    ein = calcein_i_calc(i, zs);
    
    fprintf(stderr, "lens id = %-2d [%-7s] : ", i + 1, inttolmodel(model_lens[i]));
    if(ein > 0.0){ 
      mein = sigma_crit() * M_PI * thetator(ein) * thetator(ein);
      if(outformat_exp == 0){
	fprintf(stderr, "%8.4f [arcsec]  %e [M_Sun/h]\n", ein, mein); 
      } else {
	fprintf(stderr, "%e [arcsec]  %e [M_Sun/h]\n", ein, mein); 
      }

      fprintf(fptr, "%e %3d %e %e\n", zs, i + 1, ein, mein); 
    } else { fprintf(stderr, "N/A\n"); }
  }


  fprintf(stderr, "\n");

  fclose(fptr);
  
  return;
}

double calcein_i_calc(int i, double zs)
{
  double ein;
  
  if(para_lens[i][0] >= zs){
    ein = CALCEIN_NAN;
  } else {
    set_distance_lpl_i(lens_lpl_id[i]);
    ein = calcein_i(i);
  }

  return ein;
}

double calcein_i(int i)
{
  double ein, x;

  ein = CALCEIN_NAN;

  switch(model_lens[i]){
  case 1:
    /* Einstein radius for L/L_* = 1, a/a_* = 1 member */
    ein = calcein_jaffe(i, 0.0);
    break;
    
  case 2:
    ein = calcein_nfw(i);
      break;
      
  case 3:
    x = b_sie(para_lens[i][1], 1.0) * (b_sie(para_lens[i][1], 1.0) - 2.0 * para_lens[i][6]);
    if(x  > 0.0) ein = sqrt(x);
    break;
    
  case 4:
    ein = calcein_jaffe(i, para_lens[i][7]);
    break;
    
  case 5:
    ein = sqrt(re2_point(para_lens[i][1]));
    break;
    
  case 9:
    ein = calcein_hern(i);
    break;
    
  case 10:
    ein = calcein_nfw(i);
    break;
    
  case 11:
    ein = calcein_hern(i);
    break;
    
  case 12:
    ein = calcein_pow(i);
    break;
    
  case 13:
    ein = calcein_pow(i);
    break;
    
  case 14:
    ein = calcein_gnfw(i);
    break;
    
  case 15:
    ein = calcein_gnfw(i);
    break;
    
  case 16:
    ein = calcein_sers(i);
    break;
    
  case 17:
    ein = calcein_sers(i);
    break;
    
  case 18:
    ein = calcein_tnfw(i);
    break;
    
  case 19:
    ein = calcein_tnfw(i);
    break;
    
  case 20:
    ein = calcein_ein(i);
    break;
    
  case 21:
    ein = calcein_ein(i);
    break;

  case 22:
    ein = calcein_nfw(i);
    break;

  case 23:
    ein = calcein_hern(i);
    break;
  }

  return ein;
}

double calcein_jaffe(int i, double rco)
{
  double r;

  a_sav = para_lens[i][6];
  s_sav = rco;
  b_sav = b_sie(para_lens[i][1], 1.0);
  
  if((calcein_jaffe_func(smallcore) > 0.0) && (calcein_jaffe_func(XMAX_CALCEIN) < 0.0)){
    r = gsl_zbrent(calcein_jaffe_func, smallcore, XMAX_CALCEIN, TOL_ZBRENT_CALCEIN);
  } else {
    r = CALCEIN_NAN;
  }

  return r;
}

double calcein_nfw(int i)
{
  double x, r, c;

  if(nfw_users == 0){
    b_sav = b_func(para_lens[i][1], para_lens[i][6]);
    x = rtotheta(rs(para_lens[i][1], para_lens[i][6]));
  } else {
    x = para_lens[i][6];
    c = rs(para_lens[i][1], 1.0) / thetator(para_lens[i][6]);
    b_sav = b_func(para_lens[i][1], c);
  }
  if((calcein_nfw_func(smallcore) > 0.0) && (calcein_nfw_func(XMAX_CALCEIN) < 0.0)){
    r = gsl_zbrent(calcein_nfw_func, smallcore, XMAX_CALCEIN, TOL_ZBRENT_CALCEIN) * x;
  } else {
    r = CALCEIN_NAN;
  }

  return r;
}

double calcein_gnfw(int i)
{
  double x, r, c;

  if(gnfw_usetab == 1) gnfw_maketable();

  if(nfw_users == 0){
    c = (2.0 - para_lens[i][7]) * para_lens[i][6];
    b_sav = b_func_gnfw(para_lens[i][1], c, para_lens[i][7]);
    x = rtotheta(rs(para_lens[i][1], c));
  } else {
    x = para_lens[i][6];
    c = rs(para_lens[i][1], 1.0) / thetator(para_lens[i][6]);
    b_sav = b_func_gnfw(para_lens[i][1], c, para_lens[i][7]);
  }

  if((calcein_gnfw_func(smallcore) > 0.0) && (calcein_gnfw_func(XMAX_CALCEIN) < 0.0)){
    r = gsl_zbrent(calcein_gnfw_func, smallcore, XMAX_CALCEIN, TOL_ZBRENT_CALCEIN) * x;
  } else {
    r = CALCEIN_NAN;
  }

  return r;
}

double calcein_tnfw(int i)
{
  double x, r, c;

  if(nfw_users == 0){
    b_sav = b_func(para_lens[i][1], para_lens[i][6]);
    x = rtotheta(rs(para_lens[i][1], para_lens[i][6]));
  } else {
    x = para_lens[i][6];
    c = rs(para_lens[i][1], 1.0) / thetator(para_lens[i][6]);
    b_sav = b_func(para_lens[i][1], c);
  }
  tnfw_set_tau(para_lens[i][7] * para_lens[i][6]);
  if((calcein_tnfw_func(smallcore) > 0.0) && (calcein_tnfw_func(XMAX_CALCEIN) < 0.0)){
    r = gsl_zbrent(calcein_tnfw_func, smallcore, XMAX_CALCEIN, TOL_ZBRENT_CALCEIN) * x;
  } else {
    r = CALCEIN_NAN;
  }

  return r;
}

double calcein_hern(int i)
{
  double x, r;

  b_sav = b_func_hern(para_lens[i][1], para_lens[i][6]);
  x = para_lens[i][6];
  if((calcein_hern_func(smallcore) > 0.0) && (calcein_hern_func(XMAX_CALCEIN) < 0.0)){
    r = gsl_zbrent(calcein_hern_func, smallcore, XMAX_CALCEIN, TOL_ZBRENT_CALCEIN) * x;
  } else {
    r = CALCEIN_NAN;
  }

  return r;
}

double calcein_sers(int i)
{
  double x, r;

  x = para_lens[i][6] * bnn_sers(para_lens[i][7]);
  b_sav = b_func_sers(para_lens[i][1], x, para_lens[i][7]);
  if((calcein_sers_func(smallcore) > 0.0) && (calcein_sers_func(XMAX_CALCEIN) < 0.0)){
    r = gsl_zbrent(calcein_sers_func, smallcore, XMAX_CALCEIN, TOL_ZBRENT_CALCEIN) * x;
  } else {
    r = CALCEIN_NAN;
  }

  return r;
}

double calcein_pow(int i)
{
  double r;

  if((para_lens[i][1] > para_lens[i][0]) && (para_lens[i][7] > 1.0)){
    r = pow(1.0 / fac_pert(para_lens[i][1]), 1.0 / (1.0 - para_lens[i][7])) * para_lens[i][6];
  } else {
    r = CALCEIN_NAN;
  }

  return r;
}

double calcein_ein(int i)
{
  double x, r, c;

  if(ein_usetab == 1) ein_maketable();

  if(nfw_users == 0){
    b_sav = b_func_ein(para_lens[i][1], para_lens[i][6], para_lens[i][7]);
    x = rtotheta(rs(para_lens[i][1], para_lens[i][6]));
  } else {
    x = para_lens[i][6];
    c = rs(para_lens[i][1], 1.0) / thetator(para_lens[i][6]);
    b_sav = b_func_ein(para_lens[i][1], c, para_lens[i][7]);
  }
  if((calcein_ein_func(smallcore) > 0.0) && (calcein_ein_func(XMAX_CALCEIN) < 0.0)){
    r = gsl_zbrent(calcein_ein_func, smallcore, XMAX_CALCEIN, TOL_ZBRENT_CALCEIN) * x;
  } else {
    r = CALCEIN_NAN;
  }

  return r;
}

/*--------------------------------------------------------------
  for solving equations
*/

double calcein_jaffe_func(double x)
{
  double r, f1, f2;

  f1 = sqrt(a_sav * a_sav + x * x) + a_sav;
  f2 = sqrt(s_sav * s_sav + x * x) + s_sav;
  r = b_sav * (f1 - f2) / (f1 * f2) - 1.0;
  return r;
}

double calcein_nfw_func(double x)
{
  return b_sav * dphi_nfw_dl(x) / x - 1.0;
}

double calcein_gnfw_func(double x)
{
  return b_sav * dphi_gnfw_dl(x) / x - 1.0;
}

double calcein_tnfw_func(double x)
{
  return b_sav * dphi_tnfw_dl(x) / x - 1.0;
}

double calcein_hern_func(double x)
{
  return b_sav * dphi_hern_dl(x) / x - 1.0;
}

double calcein_sers_func(double x)
{
  return b_sav * dphi_sers_dl(x) / x - 1.0;
}

double calcein_ein_func(double x)
{
  return b_sav * dphi_ein_dl(x) / x - 1.0;
}

/*--------------------------------------------------------------
  calculate circular-average of kappa
*/

void kappa_rad_out(double zs, double x0, double y0, double r1, double r2, int n, int lensid)
{
  double kapbin[NMAX_KAPBIN + 1], rbin[NMAX_KAPBIN + 1];

  return kappa_rad(zs, x0, y0, r1, r2, n, lensid, kapbin, rbin, 1);
}

void kappa_cum_out(double zs, double x0, double y0, double r1, double r2, int n, int lensid)
{
  double kapbin[NMAX_KAPBIN + 1], rbin[NMAX_KAPBIN + 1];

  return kappa_cum(zs, x0, y0, r1, r2, n, lensid, kapbin, rbin, 1);
}

void kappa_rad(double zs, double x0, double y0, double r1, double r2, int n, int lensid, double kapbin[NMAX_KAPBIN + 1], double rbin[NMAX_KAPBIN + 1], int verb)
{
  int k, nn;
  double r, hh;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  if(verb == 1){
    fprintf(stderr, "######## computing radial profile of kappa\n");
    fprintf(stderr, " zs = %e,  lens id = %d,  step = %d\n", zs, lensid, n);
    fprintf(stderr, " center = (%e, %e),  range = [%e, %e]\n", x0, y0, r1, r2);
    
    sprintf(fname, "%s_kaprad.dat", file_prefix);
    fprintf(stderr, " output file name = %s\n\n", fname);
    
    fptr = fopen(fname, "w");
    if(fptr == NULL) terminator("failed at fopen (kappa_rad)");
  }

  set_distance_lpl_zs(zs);

  if(n < 0){
    nn = n * (-1);
  } else if(n > 0){
    nn = n;
  } else {
    nn = n;
  }

  if((r1 <= 0.0) || (r1 >= r2) || (nn > NMAX_KAPBIN))
    terminator("invalid parameter (kappa_rad)"); 

  if(n > 0){
    hh = (r2 - r1) / ((double)nn);
  } else if(n < 0) {
    hh = log10(r2 / r1) / ((double)nn);
  } else {
    hh = 0.0;
  }
  
  for(k=0;k<=nn;k++){
    if(n > 0){
      r = r1 + ((double)k) * hh;
    } else if(n < 0){
      r = r1 * pow(10.0, ((double)k) * hh);
    } else {
      r = r1;
    }
    if(k == nn){ r = r2; }
    if(k == 0){ r = r1; }

    rbin[k] = r;
    kapbin[k] = calc_kappa_ave(r, x0, y0, lensid);

    if(verb == 1){
      fprintf(stderr, "%e %e\n", r, kapbin[k]);
      fprintf(fptr, "%e %e\n", r, kapbin[k]);
    }
  }

  if(verb == 1){
    fprintf(stderr, "\n");
    fclose(fptr);
  }

  return;
}

void kappa_cum(double zs, double x0, double y0, double r1, double r2, int n, int lensid, double kapbin[NMAX_KAPBIN + 1], double rbin[NMAX_KAPBIN + 1], int verb)
{
  int k, nn;
  double r, hh;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  if(verb == 1){
    fprintf(stderr, "######## computing cumulative profile of kappa\n");
    fprintf(stderr, " zs = %e,  lens id = %d,  step = %d\n", zs, lensid, n);
    fprintf(stderr, " center = (%e, %e), range = [%e, %e]\n", x0, y0, r1, r2);
    
    sprintf(fname, "%s_kapcum.dat", file_prefix);
    fprintf(stderr, " output file name = %s\n\n", fname);
    
    fptr = fopen(fname, "w");
    if(fptr == NULL) terminator("failed at fopen (kappa_cum)");
  }

  set_distance_lpl_zs(zs);

  if(n < 0){
    nn = n * (-1);
  } else if(n > 0){
    nn = n;
  } else {
    nn = n;
  }

  if((r1 <= 0.0) || (r1 >= r2) || (nn > NMAX_KAPBIN))
    terminator("invalid parameter (kappa_rad)"); 

  if(n > 0){
    hh = (r2 - r1) / ((double)nn);
  } else if(n < 0) {
    hh = log10(r2 / r1) / ((double)nn);
  } else {
    hh = 0.0;
  }
  
  for(k=0;k<=nn;k++){
    if(n > 0){
      r = r1 + ((double)k) * hh;
    } else if(n < 0){
      r = r1 * pow(10.0, ((double)k) * hh);
    } else {
      r = r1;
    }
    if(k == nn){ r = r2; }
    if(k == 0){ r = r1; }

    rbin[k] = r;
    kapbin[k] = calc_kappa_cum(r, x0, y0, lensid);

    if(verb == 1){
      fprintf(stderr, "%e %e\n", r, kapbin[k]);
      fprintf(fptr, "%e %e\n", r, kapbin[k]);
    }
  }

  if(verb == 1){
    fprintf(stderr, "\n");
    fclose(fptr);
  }

  return;
}

double calc_kappa_cum(double r, double x0, double y0, int lensid)
{
  
  lensid_calkap_sav = lensid;
  x0_calkap_sav = x0;
  y0_calkap_sav = y0;

  return 2.0 * gsl_qgaus(calc_kappa_cum_func, 0.0, r) / (r * r);
}

double calc_kappa_cum_func(double r)
{
  return r * calc_kappa_ave(r, x0_calkap_sav, y0_calkap_sav, lensid_calkap_sav);
}

double calc_kappa_ave(double r, double x0, double y0, int lensid)
{
  
  lensid_calkap_sav = lensid;
  x0_calkap_sav = x0;
  y0_calkap_sav = y0;
  r_calkap_sav = r;

  return gsl_romberg3(calc_kappa_ave_func, 0.0, 2.0 * M_PI, TOL_ROMBERG_AVE) / (2.0 * M_PI);
}

double calc_kappa_ave_func(double t)
{
  double x, y;
  double pout[NPAR_LMODEL];
	   
  x = x0_calkap_sav + r_calkap_sav * cos(t);
  y = y0_calkap_sav + r_calkap_sav * sin(t);

  lensmodel(x, y, pout, 0, lensid_calkap_sav);
  
  return pout[3];
}

/*--------------------------------------------------------------
  calculate Einstein radius in another way
*/

void calcein2(double zs, double x0, double y0, int lensid)
{
  double ein, mein;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  sprintf(fname, "%s_ein2.dat", file_prefix);

  fprintf(stderr, "######## calculating Einstein radius (from average kappa)\n");
  fprintf(stderr, " zs = %e,  lens id = %d\n", zs, lensid);
  fprintf(stderr, " center = (%e, %e)\n", x0, y0);
  fprintf(stderr, " output file name = %s \n\n", fname);

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (calcein2)");

  set_distance_lpl_zs(zs);
  
  ein = calcein2_calc(zs, x0, y0, lensid);

  if(ein > 0.0){ 
    mein = sigma_crit() * M_PI * thetator(ein) * thetator(ein);
    if(outformat_exp == 0){
      fprintf(stderr, "%8.4f [arcsec]  %e [M_Sun/h]\n", ein, mein); 
    } else {
      fprintf(stderr, "%e [arcsec]  %e [M_Sun/h]\n", ein, mein); 
    }
    fprintf(fptr, "%e %3d %e %e %e %e\n", zs, lensid, x0, y0, ein, mein); 
  } else { fprintf(stderr, "ein = N/A\n"); }

  fprintf(stderr, "\n");

  fclose(fptr);
  
  return;
}

double calcein2_calc(double zs, double x0, double y0, int lensid)
{
  double ein, rmax, rmin;

  lensid_calkap_sav = lensid;
  x0_calkap_sav = x0;
  y0_calkap_sav = y0;

  /* infer maximum from the box size */
  rmax = (xmax - xmin) + (ymax - ymin);
  rmin = rmax * TOL_ZBRENT_CALCEIN;
  
  if((calcein2_func(rmax) < 0.0) && (calcein2_func(rmin) > 0.0)){
    ein = gsl_zbrent(calcein2_func, rmin, rmax, TOL_ZBRENT_CALCEIN);
  } else {
    ein = CALCEIN_NAN;
  } 

  return ein;
}

double calcein2_func(double r)
{
  return calc_kappa_cum(r, x0_calkap_sav, y0_calkap_sav, lensid_calkap_sav) - 1.0;
}

/*--------------------------------------------------------------
  calculate masses, virial radii, etc.
*/

void calcmr(void)
{
  int i;
  double mt, md, rd;

  fprintf(stderr, "######## calculating masses and radii\n");
  fprintf(stderr, " Delta = ");
  if(flag_hodensity == 1){
    fprintf(stderr, "%9.3f times mean density\n", hodensity);
  } else if(flag_hodensity == 2){
    fprintf(stderr, "%9.3f times critical density\n", hodensity);
  } else {
    fprintf(stderr, "Delta_vir\n");
  }

  fprintf(stderr, "\n");

  for(i=0;i<num_len;i++){

    set_distance_lpl_i(lens_lpl_id[i]);
    calcmr_i(i, &mt, &md, &rd);
    
    fprintf(stderr, "lens id = %-2d [%-7s]\n", i + 1, inttolmodel(model_lens[i]));
    if(mt > 0.0){
      fprintf(stderr, " Mtot = %12e [M_Sun/h]\n", mt);
    } else {
      fprintf(stderr, " Mtot = N/A\n");
    }
    if(md > 0.0){
      fprintf(stderr, " Mdel = %12e [M_Sun/h]\n", md);
    } else {
      fprintf(stderr, " Mdel = N/A\n");
    }
    if(rd > 0.0){
      if(outformat_exp == 0){
	fprintf(stderr, " rdel = %9.5f [Mpc/h] = %9.4f [arcsec]\n", rd, rtotheta(rd));
      } else {
	fprintf(stderr, " rdel = %e [Mpc/h] = %e [arcsec]\n", rd, rtotheta(rd));
      }
    } else {
      fprintf(stderr, " rdel = N/A\n");
    }
    
    fprintf(stderr, "\n");

  }

  return;
}

void calcmr_i(int i, double *mtot, double *mdel, double *rdel)
{
  double c;

  *mtot = CALCEIN_NAN;
  *mdel = CALCEIN_NAN;
  *rdel = CALCEIN_NAN;

  switch(model_lens[i]){
  case 1:
    /* total mass for L/L_* = 1, a/a_* = 1 member */
    *mtot = ( FAC_SIE_CALCEIN ) * M_PI * (para_lens[i][1] / C_LIGHT_KMS) * (para_lens[i][1] / C_LIGHT_KMS) * thetator(para_lens[i][6] - para_lens[i][7]);
    *mdel = CALCEIN_NAN;
    *rdel = CALCEIN_NAN;
    break;
    
  case 2:
    *mtot = CALCEIN_NAN;
    *mdel = para_lens[i][1];
    *rdel = rs(para_lens[i][1], 1.0);
    break;
      
  case 4:
    *mtot = ( FAC_SIE_CALCEIN ) * M_PI * (para_lens[i][1] / C_LIGHT_KMS) * (para_lens[i][1] / C_LIGHT_KMS) * thetator(para_lens[i][6] - para_lens[i][7]);
    *mdel = CALCEIN_NAN;
    *rdel = CALCEIN_NAN;
    break;
    
  case 5:
    *mtot = para_lens[i][1];
    *mdel = CALCEIN_NAN;
    *rdel = CALCEIN_NAN;
    break;
    
  case 9:
    *mtot = para_lens[i][1];
    *mdel = CALCEIN_NAN;
    *rdel = CALCEIN_NAN;
    break;
    
  case 10:
    *mtot = CALCEIN_NAN;
    *mdel = para_lens[i][1];
    *rdel = rs(para_lens[i][1], 1.0);
    break;
    
  case 11:
    *mtot = para_lens[i][1];
    *mdel = CALCEIN_NAN;
    *rdel = CALCEIN_NAN;
    break;
    
  case 14:
    *mtot = CALCEIN_NAN;
    *mdel = para_lens[i][1];
    *rdel = rs(para_lens[i][1], 1.0);
    break;
    
  case 15:
    *mtot = CALCEIN_NAN;
    *mdel = para_lens[i][1];
    *rdel = rs(para_lens[i][1], 1.0);
    break;
    
  case 16:
    *mtot = para_lens[i][1];
    *mdel = CALCEIN_NAN;
    *rdel = CALCEIN_NAN;
    break;
    
  case 17:
    *mtot = para_lens[i][1];
    *mdel = CALCEIN_NAN;
    *rdel = CALCEIN_NAN;
    break;
    
  case 18:
    if(nfw_users == 0){
      c = para_lens[i][6];
    } else {
      c = rs(para_lens[i][1], 1.0) / thetator(para_lens[i][6]);
      printf("%e\n", c);
    }
    *mtot = para_lens[i][1] * tnfw_mtot(c, c * para_lens[i][7]);
    *mdel = para_lens[i][1];
    *rdel = rs(para_lens[i][1], 1.0);
    break;
    
  case 19:
    if(nfw_users == 0){
      c = para_lens[i][6];
    } else {
      c = rs(para_lens[i][1], 1.0) / thetator(para_lens[i][6]);
    }
    *mtot = para_lens[i][1] * tnfw_mtot(c, c * para_lens[i][7]);
    *mdel = para_lens[i][1];
    *rdel = rs(para_lens[i][1], 1.0);
    break;
    
  case 20:
    *mtot = para_lens[i][1];
    *mdel = CALCEIN_NAN;
    *rdel = CALCEIN_NAN;
    break;
    
  case 21:
    *mtot = para_lens[i][1];
    *mdel = CALCEIN_NAN;
    *rdel = CALCEIN_NAN;
    break;
    
  case 22:
    *mtot = CALCEIN_NAN;
    *mdel = para_lens[i][1];
    *rdel = rs(para_lens[i][1], 1.0);
    break;

  case 23:
    *mtot = para_lens[i][1];
    *mdel = CALCEIN_NAN;
    *rdel = CALCEIN_NAN;
    break;
  }

  return;
}

