#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

#define NM 7
static char smodelname[NM][10] = {"gauss", "sersic", "tophat", "moffat", "jaffe", "srcs", "point"};

static double xsmin, xsmax, ysmin, ysmax;
static int nxtab, nytab;
static int nmax_srcsbin;
static int *nsrcstab;
static int *isrcstab;
static int tflag;

/*--------------------------------------------------------------
  source model
*/

double sourcemodel(double x, double y, int i, double pxx, double pxy, double pyx, double pyy, double pix)
{
  double f = 0.0;

  switch(model_ext[i]){

  case 1:
    f = para_ext[i][1] * source_all(1, x, y, para_ext[i][2], para_ext[i][3], para_ext[i][4], para_ext[i][5], para_ext[i][6], para_ext[i][7], pxx, pxy, pyx, pyy, pix);
    break;

  case 2:
    f = para_ext[i][1] * source_all(2, x, y, para_ext[i][2], para_ext[i][3], para_ext[i][4], para_ext[i][5], para_ext[i][6], para_ext[i][7], pxx, pxy, pyx, pyy, pix);
    break;

  case 3:
    f = para_ext[i][1] * source_all(3, x, y, para_ext[i][2], para_ext[i][3], para_ext[i][4], para_ext[i][5], para_ext[i][6], para_ext[i][7], pxx, pxy, pyx, pyy, pix);
    break;

  case 4:
    f = para_ext[i][1] * source_all(4, x, y, para_ext[i][2], para_ext[i][3], para_ext[i][4], para_ext[i][5], para_ext[i][6], para_ext[i][7], pxx, pxy, pyx, pyy, pix);
    break;

  case 5:
    f = para_ext[i][1] * source_all(5, x, y, para_ext[i][2], para_ext[i][3], para_ext[i][4], para_ext[i][5], para_ext[i][6], para_ext[i][7], pxx, pxy, pyx, pyy, pix);
    break;

  case 6:
    f = calc_srcs(x, y, i, 0, pxx, pxy, pyx, pyy, pix);
    break;

  case 7:
    f = 0.0;
    break;
  }

  return f;
}

int emodeltoint(char *model)
{
  int i;

  for(i=0;i<NM;i++){
    if(strcmp(model, smodelname[i]) == 0){
    return i + 1;
    } 
  }
   
  return 0;
}

char* inttoemodel(int i)
{
  return smodelname[i - 1];
}

/*--------------------------------------------------------------
  sources
*/

double calc_srcs(double x, double y, int i, int flag_psf, double pxx, double pxy, double pyx, double pyy, double pix)
{
  int j, jj, ix, iy, k;
  double f = 0.0;

  if(flag_srcsbin == 0){
    if(flag_psf == 0){
      for(j=0;j<num_src;j++){
	f = f + calc_srcs_obj(j, x, y, pxx, pxy, pyx, pyy, pix);
      }
    } else {
      for(j=0;j<num_src;j++){
	f = f + calc_srcs_psf(j, x, y);
      }
    }
  } else {
    if(tflag != TFLAG_VALUE) tab_calc_src();
    
    ix = (int)((x - xsmin) / srcsbinsize);
    iy = (int)((y - ysmin) / srcsbinsize);
    k = ix + iy * nxtab;

    if((k >= 0) && (k < (nxtab * nytab))){
      if(flag_psf == 0){
	for(jj=0;jj<nsrcstab[k];jj++){
	  j = isrcstab[jj + k * nmax_srcsbin];
	  f = f + calc_srcs_obj(j, x, y, pxx, pxy, pyx, pyy, pix);
	}
      } else {
	for(jj=0;jj<nsrcstab[k];jj++){
	  j = isrcstab[jj + k * nmax_srcsbin];
	  f = f + calc_srcs_psf(j, x, y);
	}
      }
    }
  }
  
  f = f * para_ext[i][1];
    
  return f;
}

double calc_srcs_obj(int j, double x, double y, double pxx, double pxy, double pyx, double pyy, double pix)
{
  double f;

  if((para_srcs[5 + j * NPAR_SRC] > 0.0) && (para_srcs[6 + j * NPAR_SRC] > 0.0)){
    f = para_srcs[j * NPAR_SRC] * source_all(2, x, y, para_srcs[1 + j * NPAR_SRC], para_srcs[2 + j * NPAR_SRC], para_srcs[3 + j * NPAR_SRC], para_srcs[4 + j * NPAR_SRC], para_srcs[5 + j * NPAR_SRC], para_srcs[6 + j * NPAR_SRC], pxx, pxy, pyx, pyy, pix);
  } else {
    f = 0.0;
  }

  return f;
}

double calc_srcs_psf(int j, double x, double y)
{
  double f;
  
  if((para_srcs[5 + j * NPAR_SRC] <= 0.0) || (para_srcs[6 + j * NPAR_SRC] <= 0.0)){
    f = para_srcs[j * NPAR_SRC] * source_psf_pix(x, y, para_srcs[1 + j * NPAR_SRC], para_srcs[2 + j * NPAR_SRC], pix_ext);
  } else {
    f = 0.0;
  }

  return f;
}

void tab_calc_src(void)
{
  int j, ix, iy, k, ix2, iy2, ixs, iys, nn;

  xsmin = XYSMIN_SET;
  xsmax = (-1.0) * XYSMIN_SET;
  ysmin = XYSMIN_SET;
  ysmax = (-1.0) * XYSMIN_SET;

  for(j=0;j<num_src;j++){
    if(xsmin > para_srcs[1 + j * NPAR_SRC]) xsmin = para_srcs[1 + j * NPAR_SRC];
    if(xsmax < para_srcs[1 + j * NPAR_SRC]) xsmax = para_srcs[1 + j * NPAR_SRC];
    if(ysmin > para_srcs[2 + j * NPAR_SRC]) ysmin = para_srcs[2 + j * NPAR_SRC];
    if(ysmax < para_srcs[2 + j * NPAR_SRC]) ysmax = para_srcs[2 + j * NPAR_SRC];
  }

  nxtab = (int)((xsmax - xsmin) / srcsbinsize) + 1;
  nytab = (int)((ysmax - ysmin) / srcsbinsize) + 1;

  nn = nxtab * nytab;

  nmax_srcsbin = (100 * num_src) / nn;
  if(nmax_srcsbin < 10) nmax_srcsbin = 10;
  
  if((nn * nmax_srcsbin) > NMAX_PIXEL_POINT)
    terminator("srcs binning failed"); 

  nsrcstab = (int*)malloc(sizeof(int) * (nn));
  if(nsrcstab == NULL) terminator("memory allocation failed");

  isrcstab = (int*)malloc(sizeof(int) * (nn * nmax_srcsbin));
  if(isrcstab == NULL) terminator("memory allocation failed");
  
  for(j=0;j<(nxtab*nytab);j++){
    nsrcstab[j] = 0;
  }

  for(j=0;j<num_src;j++){
    /* big sources in all bins */
    if(para_srcs[5 + j * NPAR_SRC] > (srcsbinsize / 10.0)){
      for(k=0;k<nn;k++){
	if(nsrcstab[k] >= nmax_srcsbin) terminator("srcs binning failed");
	isrcstab[nsrcstab[k] + k * nmax_srcsbin] = j;
	nsrcstab[k]++;
      }
    } else {
      ix = (int)((para_srcs[1 + j * NPAR_SRC] - xsmin) / srcsbinsize);
      iy = (int)((para_srcs[2 + j * NPAR_SRC] - ysmin) / srcsbinsize);
      
      for(ix2=(-1);ix2<=1;ix2++){
	ixs = ix + ix2;
	if((ixs >= 0) && (ixs < nxtab)){
	  for(iy2=(-1);iy2<=1;iy2++){
	    iys = iy + iy2;
	    if((iys >= 0) && (iys < nytab)){
	      k = ixs + iys * nxtab;
	      if(nsrcstab[k] >= nmax_srcsbin) terminator("srcs binning failed");
	      isrcstab[nsrcstab[k] + k * nmax_srcsbin] = j;
	      nsrcstab[k]++;
	    }
	  }
	}
      }
    }
  }

  tflag = TFLAG_VALUE;

  return;
}

void unset_tab_calc_src(void)
{
  if(tflag == TFLAG_VALUE){
    free(nsrcstab);
    free(isrcstab);
    tflag = 0;
  }

  return;
}

/*--------------------------------------------------------------
  all sources (1 = Gaussian, 2 = Sersic, 3 = tophat, 4 = moffat, 5 = jaffe)
*/

double source_all(int id, double x, double y, double x0, double y0, double e, double pa, double r0, double n, double pxx, double pxy, double pyx, double pyy, double pix)
{
  int i, j, np;
  double f, dx, dy, dsx, dsy, h, hh, muinv, dr2, dr2s;

  f = 0.0;

  muinv = fabs((1.0 - pxx) * (1.0 - pyy) - pxy * pyx + imag_ceil);
  dr2s = (x - x0) * (x - x0) + (y - y0) * (y - y0);
  dr2 = dr2s / muinv;
  
  if(((dr2s < (source_refr0 * source_refr0 * r0 * r0)) || (dr2 < (4.0 * pix_ext * pix_ext))) && (id != 3) && (flag_extref != 0) && (pix > 0.0)){
  
    if(dr2<(pix_ext * pix_ext)){
      np = num_pixint * 4;
    } else {
      np = num_pixint;
    }

    hh = 1.0 / ((double)np);
    h = pix / ((double)np);
    
    for(i=0;i<np;i++){
      for(j=0;j<np;j++){
	dx = (-0.5) * pix + (((double)i) + 0.5) * h;
	dy = (-0.5) * pix + (((double)j) + 0.5) * h;
	dsx = (1.0 - pxx) * dx - pxy * dy;
	dsy = (1.0 - pyy) * dy - pyx * dx;
	if(id == 1) f = f + source_gauss(x + dsx, y + dsy, x0, y0, e, pa, r0) * hh * hh;
	if(id == 2) f = f + source_sersic(x + dsx, y + dsy, x0, y0, e, pa, r0, n) * hh * hh;
	if(id == 4) f = f + source_moffat(x + dsx, y + dsy, x0, y0, e, pa, r0, n) * hh * hh;
	if(id == 5) f = f + source_jaffe(x + dsx, y + dsy, x0, y0, e, pa, r0, n) * hh * hh;
      }
    }
    
  } else {
    if(id == 1) f = source_gauss(x, y, x0, y0, e, pa, r0);
    if(id == 2) f = source_sersic(x, y, x0, y0, e, pa, r0, n);
    if(id == 3) f = source_tophat(x, y, x0, y0, e, pa, r0);
    if(id == 4) f = source_moffat(x, y, x0, y0, e, pa, r0, n);
    if(id == 5) f = source_jaffe(x, y, x0, y0, e, pa, r0, n);
  }

  if(flag_extnorm != 0) f = f / source_all_norm(id, r0, n, pix_ext);

  return f;
}

double source_all_norm(int id, double r0, double n, double pix)
{
  double ftot = 1.0;

  switch(id){

  case 1:
    ftot = 2.0 * M_PI * r0 * r0 / (pix * pix);
    break;

  case 2:
    ftot = (M_PI * r0 * r0 / (pix * pix)) * calc_norm_sersic(n);
    break;

  case 3:
    ftot = M_PI * r0 * r0 / (pix * pix);
    break;

  case 4:
    ftot = M_PI * r0 * r0 / ((n - 1.0) * (pix * pix));
    break;

  case 5:
    if(n < smallcore) n = smallcore;
    ftot = 2.0 * M_PI * r0 * n / (pix * pix);
    break;

  }

  return ftot;
}

double calc_norm_sersic(double n)
{
  static int ff;
  static double norm, nn, bn, ga;

  if((n != nn) || (ff != 1)){
    ff = 1;
    nn = n;
    bn = bnn_sers(n);
    ga = gam2n1_sers(n);
    norm = bn * bn * ga;
  }

  return norm;
}

double source_psf_pix(double x, double y, double x0, double y0, double pix)
{
  int i, j, np;
  double r02, f, dx, dy, h, hh, dr2;

  if(flag_seeing == 1){
    
    f = 0.0;
    
    dr2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
    r02 = moffat_fwhmtoa(para_psf[0], para_psf[3]) * moffat_fwhmtoa(para_psf[4], para_psf[7]);
    
    if(((dr2 < (source_refr0 * source_refr0 * r02)) || (dr2 < (4.0 * pix_ext * pix_ext))) && (pix > 0.0)){
       
      if(dr2 < (pix_ext * pix_ext)){
	np = num_pixint * 4;
      } else {
	np = num_pixint;
      }

      hh = 1.0 / ((double)np);
      h = pix / ((double)np);
      
      for(i=0;i<np;i++){
	for(j=0;j<np;j++){
	  dx = (-0.5) * pix + (((double)i) + 0.5) * h;
	  dy = (-0.5) * pix + (((double)j) + 0.5) * h;
	  f = f + source_psf(x + dx, y + dy, x0, y0, pix) * hh * hh;
	}
      }
      
    } else {
      f = source_psf(x, y, x0, y0, pix);
    }

  } else {
    f = source_psf(x, y, x0, y0, pix);
  }

  return f;
}

double source_psf(double x, double y, double x0, double y0, double pix)
{
  int i, j, nnk;
  double fw1, e1, pa1, b1, fw2, e2, pa2, b2, f, x1, x2, pix2;
  double xx, yy, p, q;

  if(flag_seeing == 1){
    fw1 = para_psf[0];
    e1 = para_psf[1];
    pa1 = para_psf[2];
    b1 = para_psf[3];
    fw2 = para_psf[4];
    e2 = para_psf[5];
    pa2 = para_psf[6];
    b2 = para_psf[7];
    f = para_psf[8];
    
    x1 = source_moffat_psf(x, y, x0, y0, e1, pa1, fw1, b1); 
    /* x1 = source_gauss(x, y, x0, y0, e1, pa1, fw1); */
    x2 = source_moffat_psf(x, y, x0, y0, e2, pa2, fw2, b2); 
    
    return (pix * pix) * (f * x1 + (1.0 - f) * x2);
  } else if(flag_seeing == (-1)){
    pix2 = pix / ((double)seeing_sub);
    nnk = 2 * fpsf_nk + 1;
    xx = ((x - x0) / pix2) + fpsf_nk;
    yy = ((y - y0) / pix2) + fpsf_nk;
    i = (int)xx;
    j = (int)yy;
    p = xx - (double)i;
    q = yy - (double)j;

    if((i < 0) || (i > (nnk - 2)) || (j < 0) || (j > (nnk - 2))){ return 0.0; }
    if((p < 0.0) || (p > 1.0) || (q < 0.0) || (q > 1.0)){ return 0.0; }
   
    return ((1.0 - p) * (1.0 - q) * array_fpsf[i + j * nnk] + p * (1.0 - q) * array_fpsf[i + 1 + j * nnk]
	     + (1.0 - p) * q * array_fpsf[i + (j + 1) * nnk] + p * q * array_fpsf[i + 1 + (j + 1) * nnk]) * ((double)(seeing_sub * seeing_sub));
  } else {
    terminator("psf must be set for point source images");
  }

  return 0.0;
}

/*--------------------------------------------------------------
  Gaussian
*/

double source_gauss(double x, double y, double x0, double y0, double e, double pa, double r0)
{
  double u;

  if(source_checkdis(x, y, x0, y0, r0) > 0) return 0.0;

  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(r0, 0.0);

  u = ucalc(x - x0, y - y0, e, pa);

  return exp((-0.5) * u / (r0 * r0));
}

/*--------------------------------------------------------------
  sersic
*/

double source_sersic(double x, double y, double x0, double y0, double e, double pa, double r0, double n)
{
  static int ff;
  double u;
  static double nn, bn;

  if(source_checkdis(x, y, x0, y0, r0) > 0) return 0.0;

  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(r0, 0.0);
  checkmodelpar_mineq(n, SERSIC_N_MIN);
  checkmodelpar_maxeq(n, SERSIC_N_MAX);

  if((n != nn) || (ff != 1)){
    ff = 1;
    nn = n;
    bn = bn_sers(n);
  }

  u = ucalc(x - x0, y - y0, e, pa);

  return exp(bn * (-1.0) * pow(u / (r0 * r0), 0.5 / n));
}

/*--------------------------------------------------------------
  tophat
*/

double source_tophat(double x, double y, double x0, double y0, double e, double pa, double r0)
{
  double u, r2;
  
  if(source_checkdis(x, y, x0, y0, r0) > 0) return 0.0;

  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(r0, 0.0);

  u = ucalc(x - x0, y - y0, e, pa);
  r2 = r0 * r0;

  if(u <= r2){
    return 1.0;
  } else {
    return 0.0;
  }
}

/*--------------------------------------------------------------
  moffat
*/

double source_moffat_psf(double x, double y, double x0, double y0, double e, double pa, double fwhm, double b)
{
  double a, u, r2, ftotinv;
  
  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(b, 1.0);
  a = moffat_fwhmtoa(fwhm, b);
  checkmodelpar_min(a, 0.0);

  u = ucalc(x - x0, y - y0, e, pa);
  r2 = a * a;
  ftotinv = (b - 1.0) / (M_PI * a * a);

  return ftotinv * pow(1.0 + (u / r2), (-1.0) * b);
}

double source_moffat(double x, double y, double x0, double y0, double e, double pa, double a, double b)
{
  double u, r2;
  
  if(source_checkdis(x, y, x0, y0, a) > 0) return 0.0;

  checkmodelpar_mineq(e, 0.0);
  checkmodelpar_max(e, 1.0);
  checkmodelpar_min(a, 0.0);
  checkmodelpar_min(b, 1.0);

  u = ucalc(x - x0, y - y0, e, pa);
  r2 = a * a;

  return pow(1.0 + (u / r2), (-1.0) * b);
}

double moffat_fwhmtoa(double fwhm, double b)
{
  static double bb, aa;

  if(bb != b){
    bb = b;
    aa = 0.5 / sqrt(pow(2.0, 1.0 / b) - 1.0);
  }

  return fwhm * aa;
}

/*--------------------------------------------------------------
  sersic
*/

double source_jaffe(double x, double y, double x0, double y0, double e, double pa, double a, double rco)
{
  double u, f1, f2;
  
  if(source_checkdis(x, y, x0, y0, a) > 0) return 0.0;

  if(a > rco){
    checkmodelpar_mineq(e, 0.0);
    checkmodelpar_max(e, 1.0);
    checkmodelpar_min(a, 0.0);
    checkmodelpar_mineq(rco, 0.0);
    
    if(rco < smallcore) rco = smallcore;
    
    u = ucalc(x - x0, y - y0, e, pa);
    
    f1 = 1.0 / sqrt(rco * rco + u);
    f2 = 1.0 / sqrt(a * a + u);
    
    return (f1 - f2) / ((1.0 / rco) - (1.0 / a));
  } else {
    return 0.0;
  } 
}

/*--------------------------------------------------------------
  misc
*/

int source_checkdis(double x, double y, double x0, double y0, double r0)
{
  int i;
  double dx, dy, rr;

  i = 0;

  dx = fabs(x - x0);
  dy = fabs(y - y0);
  rr = source_calcr0 * r0;
  
  if(dx > rr) i++;
  if(dy > rr) i++;

  return i;
}

double ucalc(double dx, double dy, double e, double pa)
{
  static int ff;
  double ddx, ddy;
  static double si, co, pp;

  if((pa != pp) || (ff != 1)){
    ff = 1;
    pp = pa;
    si = sin((-1.0) * pa * M_PI / 180.0);
    co = cos((-1.0) * pa * M_PI / 180.0);
  }

  ddx = co * dx - si * dy;
  ddy = si * dx + co * dy;

  return ddx * ddx / (1.0 - e) + (1.0 - e) * ddy * ddy;
}


#undef NM
