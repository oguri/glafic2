#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#include "glafic.h"

/*--------------------------------------------------------------
  set def angle and/or allocate matrix for images
*/

void ext_set_table_all(int i)
{
  if(num_lpl == 1){
    if(i_ext_fid < 0){   /* set table only once for a single lens plane */
      ext_set_table_single();
    }
  } else {
    ext_set_table_lpl(i);
  }

  return;
}

void ext_set_table_single(void)
{
  set_distance_facext();
  ext_set_table(i_ext_fid);
  
  return;
}

void ext_set_table_lpl(int i)
{
  set_distance_facext();
  ext_set_table(i);

  return;
}

void ext_set_table(int i)
{
  int k, nn, ntot, ntot2, ntot3, flag;
  double tx, ty, pout[NPAR_LMODEL];
  int ix, iy, ixp, ixm, iyp, iym;
  float axxp, axxm, ayyp, ayym, axyp, axym, ayxp, ayxm, ddx, ddy;

  if(num_ext < 1) terminator("setting deflections failed");

  set_distance_lpl_zs(para_ext[i][0]);

  nn = nx_ext * ny_ext;
  ntot = nn * 2;
  ntot2 = nn * (num_ext + 1);
  ntot3 = nn * 4;

  if(flag_set_array == 0){
    flag_set_array = 1;
    
    array_ext_def = (float*)malloc(sizeof(float) * ntot);
    array_ext_img = (float*)malloc(sizeof(float) * ntot2);
    array_ext_img_ori = (float*)malloc(sizeof(float) * ntot2);
    array_ext_mask = (int*)malloc(sizeof(int) * nn);
    array_ext_mag = (float*)malloc(sizeof(float) * ntot3);
    
    if((array_ext_def == NULL) || (array_ext_mag == NULL) || (array_ext_img == NULL) || (array_ext_mask == NULL))
      terminator("memory allocation failed");
  }

  for(k=0;k<nn;k++){
    flag = 1;
    if((flag_computeall == 0) && (flag_obsmask == 1)){
      if(array_obs_mask[k] > 0) flag = 0;
    }

    if(flag == 1){ 
      ktoxy_ext(k, &tx, &ty);
      lensmodel(tx, ty, pout, 1, 0);
      array_ext_def[k] = pout[0];
      array_ext_def[k + nn] = pout[1];
      array_ext_mask[k] = 0;
    } else {
      array_ext_def[k] = 0.0;
      array_ext_def[k + nn] = 0.0;
      array_ext_mask[k] = 1;
    }
  }

  for(k=0;k<nn;k++){
    ix = k % nx_ext;
    iy = (k - ix) / nx_ext;

    ixm = ix - 1;
    ixp = ix + 1;
    iym = iy - 1;
    iyp = iy + 1;

    ddx = 0.0;
    ddy = 0.0;

    if(ixm >= 0){ 
      axxm = array_ext_def[ixm + iy * nx_ext]; 
      ayxm = array_ext_def[ixm + iy * nx_ext + nn]; 
      ddx = ddx + pix_ext;
    } else { 
      axxm = array_ext_def[ix + iy * nx_ext]; 
      ayxm = array_ext_def[ix + iy * nx_ext + nn]; 
    }
    if(ixp < nx_ext){ 
      axxp = array_ext_def[ixp + iy * nx_ext]; 
      ayxp = array_ext_def[ixp + iy * nx_ext + nn]; 
      ddx = ddx + pix_ext;
    } else { 
      axxp = array_ext_def[ix + iy * nx_ext]; 
      ayxp = array_ext_def[ix + iy * nx_ext + nn]; 
    }
    if(iym >= 0){ 
      axym = array_ext_def[ix + iym * nx_ext]; 
      ayym = array_ext_def[ix + iym * nx_ext + nn]; 
      ddy = ddy + pix_ext;
    } else { 
      axym = array_ext_def[ix + iy * nx_ext]; 
      ayym = array_ext_def[ix + iy * nx_ext + nn]; 
    }
    if(iyp < ny_ext){ 
      axyp = array_ext_def[ix + iyp * nx_ext]; 
      ayyp = array_ext_def[ix + iyp * nx_ext + nn]; 
      ddy = ddy + pix_ext;
    } else { 
      axyp = array_ext_def[ix + iy * nx_ext]; 
      ayyp = array_ext_def[ix + iy * nx_ext + nn]; 
    }

    /* phi_xx */
    array_ext_mag[k] = (axxp - axxm) / ddx;
    /* phi_xy */
    array_ext_mag[k + nn] = (ayxp - ayxm) / ddy;
    /* phi_yx */
    array_ext_mag[k + 2 * nn] = (axyp - axym) / ddx;
    /* phi_yy */
    array_ext_mag[k + 3 * nn] = (ayyp - ayym) / ddy;
    
  }

  return;
}

void ext_unset_table(void)
{
  if(flag_set_array == 1){
    free(array_ext_def);
    free(array_ext_mag);
    free(array_ext_img);
    free(array_ext_img_ori);
    free(array_ext_mask);
    flag_set_array = 0;
    i_ext_fid = -1;
  }

  return;
}

/*--------------------------------------------------------------
  compute images
  isrc is id for source, set 0 if you want to compute images for all sources
*/

void ext_set_image(int isrc, int flag_source, int verb)
{
  int ii, k, nn, i, j, j0, j1, i2, j2, ni, l;
  int nx, ny, dny1, dny2, ddny, nk, nnk, npadx, npady, nnp;
  double sx, sy, x, y, pxx, pyy, pxy, pyx, pix_psf;
  double ff, re, im, ff2;
  double *gau;
  double *img;
  double rr[NMAX_POIMG][NPAR_IMAGE];
 
  fftw_complex *fft_img, *fft_gau;
  fftw_plan fft_p;

  nn = nx_ext * ny_ext;
  
  pix_psf = calc_pixpsf();

  for(ii=0;ii<num_ext;ii++){
    if(((ii + 1) == isrc) || (isrc == 0)){

      ext_set_table_all(ii);
      
      for(k=0;k<nn;k++){
	if(array_ext_mask[k] == 0){
	  ktoxy_ext(k, &x, &y);
	  
	  if(flag_source == 0){
	    sx = x - array_ext_def[k] * dis_fac_ext[ii];
	    sy = y - array_ext_def[k + nn] * dis_fac_ext[ii];
	    pxx = array_ext_mag[k] * dis_fac_ext[ii];
	    pxy = array_ext_mag[k + nn] * dis_fac_ext[ii];
	    pyx = array_ext_mag[k + 2 * nn] * dis_fac_ext[ii];
	    pyy = array_ext_mag[k + 3 * nn] * dis_fac_ext[ii];
	  } else {
	    sx = x;
	    sy = y;
	    pxx = 0.0;
	    pxy = 0.0;
	    pyx = 0.0;
	    pyy = 0.0;
	  }
	  array_ext_img[k + ii * nn] = (float)sourcemodel(sx, sy, ii, pxx, pxy, pyx, pyy, pix_psf);
	  if(flag_seeing!=0){
	    array_ext_img_ori[k + ii * nn] = array_ext_img[k + ii * nn];
	    array_ext_img[k + ii * nn] = 0.0;
	  }
	} else {
	  if(flag_seeing!= 0){
	    array_ext_img_ori[k + ii * nn] = 0.0;
	  }
	}
      }
    }
  }

  if(flag_seeing!=0){

    if(seeing_sub <= 0) terminator("invalid parameter");

    ff2 = 1.0 / ((double)(seeing_sub * seeing_sub));
    
    nk = calc_psfnk();
    nnk = nk * 2 + 1;
    
    nx = seeing_sub * nx_ext;
    npadx = nx + nk;
    
    ny = (nmax_fft / npadx) - nk;

    if(ny > (seeing_sub * ny_ext)){
      ny = seeing_sub * ny_ext;
      dny1 = ny_ext;
      dny2 = ny_ext;
      ddny = 0;
    } else {
      ny = ny - ny % seeing_sub;
      dny1 = ny / seeing_sub;
      dny2 = (ny - 2 * nk) / seeing_sub;
      ddny = nk / seeing_sub;
    }

    npady = ny + nk;

    if(ny <= (3 * nk)) terminator("PSF convolution failed");

    nnp = npadx * npady;

    ff = 1.0 / ((double)nnp);

    /*------------------------------------
      making PSF images (real and fourier)
    */
    gau = (double*)malloc(sizeof(double) * nnk * nnk);
    if(gau == NULL) terminator("memory allocation failed");

    set_psfimage(gau, nk, pix_psf);

    fft_gau = (fftw_complex *)fftw_malloc(nnp * sizeof(fftw_complex));

    for(k=0;k<nnp;k++){
      fft_gau[k][0] = 0.0;
      fft_gau[k][1] = 0.0;
    }
    
    for(j=0;j<=nk;j++){
      for(i=0;i<=nk;i++) fft_gau[i + npadx * j][0] = gau[(i + nk) + (j + nk) * nnk];
      for(i=(npadx - nk);i<npadx;i++) fft_gau[i + npadx * j][0] = gau[(i + nk - npadx) + (j + nk) * nnk];
    }
    
    for(j=(npady-nk);j<npady;j++){
      for(i=0;i<=nk;i++) fft_gau[i + npadx * j][0] = gau[(i + nk) + (j + nk - npady) * nnk];
      for(i=(npadx - nk);i<npadx;i++) fft_gau[i + npadx * j][0] = gau[(i + nk - npadx) + (j + nk - npady) * nnk];
    }
    
    fft_p = fftw_plan_dft_2d(npady, npadx, fft_gau, fft_gau, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(fft_p);

    fftw_destroy_plan(fft_p);

    /* normalize */
    for(k=0;k<nnp;k++){
      fft_gau[k][0] = fft_gau[k][0] * ff;
      fft_gau[k][1] = fft_gau[k][1] * ff;
    }

    /*------------------------------------
      end of PSF making
    */
    for(ii=0;ii<num_ext;ii++){
      if(((ii + 1) == isrc) || (isrc == 0)){
	
	if(verb == 1) fprintf(stderr, "   convolving PSF ( id = %d) ...\n", ii + 1);

	for(j0=0;j0<ny_ext;j0=j0+dny2){
	  
	  if(verb == 1){
	    if(ny_ext < (j0 + dny2)){ j1 = ny_ext; } else { j1 = j0 + dny2; }
	    fprintf(stderr, "   y-range : %5d -- %5d \n", j0 + 1, j1);
	  }
    
	  img = (double*)malloc(sizeof(double) * nx * ny);
	  if(img == NULL) terminator("memory allocation failed");
	  
	  set_psfconv_img(img, 0, nx_ext - 1, j0 - ddny, j0 - ddny + dny1 - 1, ii * nn);
	  
	  fft_img = (fftw_complex *)fftw_malloc(nnp * sizeof(fftw_complex));
	  
	  for(j=0;j<npady;j++){
	    for(i=0;i<npadx;i++){
	      if((i < nx) && (j < ny)){
		fft_img[i + npadx * j][0] = img[i + nx * j];
		fft_img[i + npadx * j][1] = 0.0;
	      } else {
		fft_img[i + npadx * j][0] = 0.0;
		fft_img[i + npadx * j][1] = 0.0;
	      }
	    }
	  }
	  
	  fft_p = fftw_plan_dft_2d(npady, npadx, fft_img, fft_img, FFTW_FORWARD, FFTW_ESTIMATE);
	  
	  fftw_execute(fft_p);
	
	  fftw_destroy_plan(fft_p);
	  
	  /* convolve with PSF in fourier space */
	  for(k=0;k<nnp;k++){
	    re = fft_img[k][0] * fft_gau[k][0] - fft_img[k][1] * fft_gau[k][1];
	    im = fft_img[k][0] * fft_gau[k][1] + fft_img[k][1] * fft_gau[k][0];
	    fft_img[k][0] = re;
	    fft_img[k][1] = im;
	  }
	  
	  fft_p = fftw_plan_dft_2d(npady, npadx, fft_img, fft_img, FFTW_BACKWARD, FFTW_ESTIMATE);
	  
	  fftw_execute(fft_p);
	  
	  fftw_destroy_plan(fft_p);
	  
	  for(j=0;j<(dny2 * seeing_sub);j++){
	    j2 = j0 + j / seeing_sub;
	    if(j2<ny_ext){
	      for(i=0;i<nx;i++){
		i2 = i / seeing_sub;
		array_ext_img[i2 + j2 * nx_ext + ii * nn] = array_ext_img[i2 + j2 * nx_ext + ii * nn] + fft_img[i + npadx * (j + ddny * seeing_sub)][0] * ff2;
	      }
	    }
	  }
	  
	  free(img);
	  fftw_free(fft_img);
	}
	
	if(verb == 1) fprintf(stderr, "\n");

      }
      
    }
    
    free(gau);
    fftw_free(fft_gau);
    
    /* adding PSFs */
    for(ii=0;ii<num_ext;ii++){
      if(((ii + 1) == isrc) || (isrc == 0)){
	
	/* for ``srcs'' */
	if(model_ext[ii] == emodeltoint("srcs")){
	  
	  for(k=0;k<nn;k++){
	    if(array_ext_mask[k] == 0){
	      ktoxy_ext(k, &x, &y);

	      array_ext_img[k + ii * nn] = array_ext_img[k + ii * nn] + (float)calc_srcs(x, y, ii, 1, 0.0, 0.0, 0.0, 0.0, pix_ext);
	    }
	  }
	}
	
	/* for ``point'' */
	if(model_ext[ii] == emodeltoint("point")){
	  
	  if(flag_source == 0){
	    poi_unset_table();
	    findimg(para_ext[ii][2], para_ext[ii][3], para_ext[ii][0], &ni, rr, 0);
	  } else {
	    ni = 1; 
	    rr[0][0] = para_ext[ii][2]; 
	    rr[0][1] = para_ext[ii][3]; 
	    rr[0][2] = 1.0;
	  }
	  
	  for(l=0;l<ni;l++){
	    for(k=0;k<nn;k++){
	      if(array_ext_mask[k] == 0){
		ktoxy_ext(k, &x, &y);
		
		array_ext_img[k + ii * nn] = array_ext_img[k + ii * nn] + para_ext[ii][1] * fabs(rr[l][2]) * source_psf_pix(x, y, rr[l][0], rr[l][1], pix_ext);
	      }
	    }
	  }
	}
	
      }
    }
  }
  
  return;
}

double calc_pixpsf(void)
{
  double pix_psf;

  if((seeing_sub > 1) && (flag_seeing != 0)){
    pix_psf = pix_ext / ((double)seeing_sub);
  } else {
    pix_psf = pix_ext;
  }

  return pix_psf;
}

int calc_psfnk(void)
{
  double sig;
  
  if(flag_seeing == 1){
    
    if(para_psf[0] > para_psf[4]){ 
      sig = INVSIG2FWHM * para_psf[0]; 
    } else { 
      sig = INVSIG2FWHM * para_psf[4]; 
    }
    /* convolve within seeing_calcsig */
    /* nk = ((int)((seeing_calcsig * sig / pix_ext) + 1)) * seeing_sub; */
    return ((int)((psfconv_size / pix_ext) + 1)) * seeing_sub;
  } else if(flag_seeing == (-1)){
    return fpsf_nk;
  } else {
    terminator("PSF not set");
  }

  return 0;
}

void set_psfimage(double *psf, int nk, double pix_psf)
{
  int i, j, k, nnk;
  double f1, f2, x2, y2;

  nnk = nk * 2 + 1;
  
  f2 = 0.0;
  if(flag_seeing == 1){
    for(j=0;j<nnk;j++){
      for(i=0;i<nnk;i++){
	x2 = ((double)(i - nk)) * pix_psf;
	y2 = ((double)(j - nk)) * pix_psf;
	f1 = source_psf_pix(x2, y2, 0.0, 0.0, pix_psf);
	psf[i + j * nnk] = f1;
	f2 = f2 + f1;
      }
    }
  } else if(flag_seeing == (-1)){
    for(j=0;j<nnk;j++){
      for(i=0;i<nnk;i++){
	f1 = array_fpsf[i + j * nnk];
	psf[i + j * nnk] = f1;
	f2 = f2 + f1;
      }
    }
  } else {
    terminator("PSF not set");
  }
  
  for(k=0;k<(nnk * nnk);k++) psf[k] = psf[k] / f2;

  return;
}

void set_psfconv_img(double *img, int ix0, int ix1, int iy0, int iy1, int in)
{
  int i, j, i2, j2, i3, j3, fb;
  int nx, ny;
  double x2, y2, sx, sy, ff;
  double xx[4], yy[4];

  nx = (ix1 - ix0 + 1) * seeing_sub;
  ny = (iy1 - iy0 + 1) * seeing_sub;
  
  fb = 0;

  for(j2=0;j2<ny;j2++){
    y2 = NPIX_SMALL_OFFSET + iy0 - 1.0 + (2.0 * j2 + seeing_sub + 1.0) / (2.0 * seeing_sub);
    j = (int)y2;
    sy = y2 - j;
    
    if((seeing_sub > 1) && (j >= 1) && (j <= (ny_ext - 3))){
      fb = 1;
      yy[0] = bicub_func2(sy + 1.0);
      yy[1] = bicub_func1(sy);
      yy[2] = bicub_func1(1.0 - sy);
      yy[3] = bicub_func2(2.0 - sy);
    }

    for(i2=0;i2<nx;i2++){
      x2 = NPIX_SMALL_OFFSET + ix0 - 1.0 + (2.0 * i2 + seeing_sub + 1.0) / (2.0 * seeing_sub);
      i = (int)x2;
      sx = x2 - i;

      if((fb == 1) && (i >= 1) && (i <= (nx_ext - 3))){
       
	/* bilinear */
	/* img[i2 + j2 * n] = (1.0 - sx) * (1.0 - sy) * array_ext_img_ori[i + j * nx_ext + in] + sx * (1.0 - sy) * array_ext_img_ori[i + 1 + j * nx_ext + in] + (1.0 - sx) * sy * array_ext_img_ori[i + (j + 1) * nx_ext + in] + sx * sy * array_ext_img_ori[i + 1 + (j + 1) * nx_ext + in]; */
	
	/* bicubic */
	xx[0] = bicub_func2(sx + 1.0);
	xx[1] = bicub_func1(sx);
	xx[2] = bicub_func1(1.0 - sx);
	xx[3] = bicub_func2(2.0 - sx);
	
	ff = 0.0;

	for(j3=0;j3<4;j3++){
	  for(i3=0;i3<4;i3++){
	    ff = ff + xx[i3] * yy[j3] * array_ext_img_ori[i + i3 - 1 + (j + j3 - 1) * nx_ext + in];
	  }
	}

	img[i2 + j2 * nx] = ff;
	
      } else if((i >=0 ) && (i <= (nx_ext - 1)) && (j >= 0) && (j <= (ny_ext - 1))){

	img[i2 + j2 * nx] = array_ext_img_ori[i + j * nx_ext + in];

      } else {

	img[i2 + j2 * nx] = 0.0;

      }
    }
    fb = 0;
  } 	 
  
  return;
}

double bicub_func1(double x)
{
  return ((bicub_a + 2.0) * x - (bicub_a + 3.0)) * x * x + 1.0;
}

double bicub_func2(double x)
{
  return (bicub_a) * (((x - 5.0) * x + 8.0) * x - 4.0);
}

/*--------------------------------------------------------------
  calculate image properties
  isrc is id for source, set 0 if you want to compute images for all sources
*/

void ext_est_image(int isrc, int flag_source, double sbth, double lwlim)
{
  int narc;
  double flux, peak, px, py, area, flux2, lwmax;
  char fname[INPUT_MAXCHAR], fname2[INPUT_MAXCHAR];
  FILE* fptr;
  FILE* fptr_a;

  if(flag_source == 0){
    fprintf(stderr, "######## computing extend image properties (lensed)\n");
  } else {
    fprintf(stderr, "######## computing extend image properties (original)\n");
  }
  sprintf(fname, "%s_extend.dat", file_prefix);
  fprintf(stderr, " output file name = %s \n", fname);
  sprintf(fname2, "%s_extend_arc.dat", file_prefix);
  fprintf(stderr, " output file name = %s \n\n", fname2);
  
  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (calcext)");

  fptr_a = fopen(fname2, "w");
  if(fptr_a == NULL) terminator("failed at fopen (calcext)");

  i_ext_fid = -1;
  ext_set_image(isrc, flag_source, 1);

  if(isrc == 0){
    fprintf(stderr, "extend id = all ( 1 -- %2d )\n", num_ext);
  } else {
    fprintf(stderr, "extend id = %-2d [%-7s] : \n", isrc, inttoemodel(model_ext[isrc - 1]));
  }
  
  ext_est_image_i(isrc, lwlim, &sbth, &flux, &peak, &px, &py, &area, &flux2, &lwmax, &narc, fptr_a);
  
  fprintf(stderr, "  flux = %e\n  peak = %e at (%e, %e)\n", flux, peak, px, py); 
  fprintf(stderr, "  area = %e [arcsec^2] and  flux = %e for threshold = %e\n", area, flux2, sbth); 
  fprintf(stderr, "  l/w(max) = %e, %d arcs with l/w > %e\n\n", lwmax, narc, lwlim); 
  
  fprintf(fptr, "%-2d %13e %13e %13e %13e %13e %13e %13e %13e %d\n", isrc, flux, peak, px, py, sbth, area, flux2, lwmax, narc);
  
  fclose(fptr); 
  fclose(fptr_a); 

  return;
}

void ext_est_image_i(int isrc, double lwlim, double *sbth, double *flux, double *peak, double *px, double *py, double *area, double *flux2, double *lwmax, int *narc, FILE* fptr)
{
  int i, nn, k, pk, kmin, kmax, l, nl;
  double aa, th, apix, lw;
  double xx, yy, rr, r1, r2, len, are, wid, x0, y0;
  int np[NMAX_ARC];
  int akpeak[NMAX_ARC], akedge1[NMAX_ARC], akedge2[NMAX_ARC];
  double afpeak[NMAX_ARC], adismax[NMAX_ARC], afarea[NMAX_ARC], afflux[NMAX_ARC];
  float *aimg;
  int *simg;
  int *aptr;

  nn = nx_ext * ny_ext;
  apix = pix_ext * pix_ext;

  *flux = 0.0;
  *peak = EXTEND_INIT_PEAK;
  pk = 0;

  aimg = (float*)malloc(sizeof(float) * nn);
  if(aimg == NULL) terminator("memory allocation failed");

  simg = (int*)malloc(sizeof(int) * nn);
  if(simg == NULL) terminator("memory allocation failed");

  aptr = (int*)malloc(sizeof(int) * nn);
  if(aptr == NULL) terminator("memory allocation failed");

  if(isrc == 0){
    for(k=0;k<nn;k++) aimg[k] = 0.0;
    for(i=0;i<num_ext;i++){
      for(k=0;k<nn;k++){
	aimg[k] = aimg[k] + array_ext_img[k + i * nn];
      }
    }
  } else {
    for(k=0;k<nn;k++){
      aimg[k] = array_ext_img[k + (isrc - 1) * nn];
    }
  }

  /* total flux, peak */
  for(k=0;k<nn;k++){
    aa = aimg[k];
    *flux = (*flux) + aa;
    if((*peak) < aa){
      *peak = aa;
      pk = k;
    }
  }
  
  /* total flux, peak */
  if((*sbth) > 0.0){
    th = (*sbth); 
  } else if(((*sbth) < 0.0) && ((*sbth) > (-1.0))){
    th = (-1.0) * (*sbth) * (*peak);
  } else {
    th = 0.5 * (*peak);
  }
  
  *sbth = th;

  ktoxy_ext(pk, px, py);

  /* flux above threshold */
  *flux2 = 0.0;
  *area = 0.0;
  kmin = nn;
  kmax = -1;

  for(k=0;k<nn;k++){
    aa = aimg[k];
    if(aa >= th){
      *flux2 = (*flux2) + aa;
      *area = (*area) + apix;
      simg[k] = 1;
      if(k > kmax){ kmax = k; }
      if(k < kmin){ kmin = k; }
    } else {
      simg[k] = 0;
    }
    aptr[k] = 0;
  }

  /* arc id, l/w */
  l = 1;
  for(k=kmin;k<=kmax;k++){
    np[l] = arc_id(k, l, simg, aptr);
    if(np[l] > 0){ 
      l++; 
    }
    if(l >= NMAX_ARC) break;
  }

  nl = l - 1;

  for(l=1;l<=nl;l++) afpeak[l] = 0.0;

  for(k=kmin;k<=kmax;k++){
    l = aptr[k];
    if(aimg[k] > afpeak[l]){
      afpeak[l] = aimg[k];
      akpeak[l] = k;
    }
  }

  for(l=1;l<=nl;l++) adismax[l] = 0.0;

  for(k=kmin;k<=kmax;k++){
    l = aptr[k];
    rr = r2_extpix(akpeak[l], k);
    if(rr > adismax[l]){
      adismax[l] = rr;
      akedge1[l] = k;
    }
  }
  
  for(l=1;l<=nl;l++){
    adismax[l] = 0.0;
    afarea[l] = 0.0;
    afflux[l] = 0.0;
  }

  for(k=kmin;k<=kmax;k++){
    l = aptr[k];
    rr = r2_extpix(akedge1[l], k);
    afarea[l] = afarea[l] + apix;
    afflux[l] = afflux[l] + aimg[k];
    if(rr > adismax[l]){
      adismax[l] = rr;
      akedge2[l] = k;
    }
  }

  *lwmax = 0.0;
  *narc = 0;
  for(l=1;l<=nl;l++){
    if(np[l] >= NMIN_ARCANA){
      r1 = sqrt(r2_extpix(akpeak[l], akedge1[l]));
      r2 = sqrt(r2_extpix(akpeak[l], akedge2[l]));
      len = r1 + r2;
      are = ((double)(np[l])) * apix;
      wid = 4.0 * are / (M_PI * len);
      lw = len / wid;
      ktoxy_ext(akedge2[l], &xx, &yy);
      if(lw > (*lwmax)) *lwmax = lw;
      if(lw > lwlim){
	(*narc)++;
	ktoxy_ext(akpeak[l], &x0, &y0);
	fprintf(fptr, "%13e %13e %13e %13e %13e\n", x0, y0, lw, afarea[l], afflux[l]);
      }
    }
  }

  free(aimg);
  free(simg);
  free(aptr);

  return;
}

int arc_id(int k, int aid, int *im, int *ia)
{
  int i;
  
  if((k < 0) || (k >= (nx_ext * ny_ext))){
    return 0;
  }

  if(((*(im + k)) != 1) || ((*(ia + k)) != 0)){
    return 0;
  }
  
  i = 1;
  *(ia + k) = aid;

  if((k % nx_ext) != 0) i = i + arc_id(k - 1, aid, im, ia);
  if((k % nx_ext) != (nx_ext - 1)) i = i + arc_id(k + 1, aid, im, ia);
  i = i + arc_id(k - nx_ext, aid, im, ia);
  i = i + arc_id(k + nx_ext, aid, im, ia);

  return i;
}

double r2_extpix(int k1, int k2)
{
  double x1, y1, x2, y2;

  ktoxy_ext(k1, &x1, &y1);
  ktoxy_ext(k2, &x2, &y2);

  return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

/*--------------------------------------------------------------
  convert int to/from x and y
*/

void ktoxy_ext(int k, double *x, double *y)
{
  int i, j;
  
  i = k % nx_ext;
  j = (k - i) / nx_ext;
  *x = xmin + pix_ext * (i + 0.5);
  *y = ymin + pix_ext * (j + 0.5);
 
  return;
}

int xytok_ext(double x, double y)
{
  int i, j;
  
  i = (int)(((x - xmin) / pix_ext) - 0.5);
  j = (int)(((y - ymin) / pix_ext) - 0.5);
 
  return i + j * nx_ext;
}

