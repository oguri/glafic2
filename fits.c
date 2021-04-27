#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <gsl/gsl_randist.h>

#include "glafic.h"

void fits_header(fitsfile *fptr);

/*--------------------------------------------------------------
  write kappa, gamma, mag map
*/

void writelens(double zs)
{
  long num[3];
  int status = 0;
  long fpixel[] = {1, 1, 1};
  float *array;
  int k, nn, ntot;
  double tx, ty, pout[NPAR_LMODEL];
  char fname[INPUT_MAXCHAR];

  sprintf(fname, "!%s_lens.fits", file_prefix);

  fprintf(stderr, "######## writing lens properties\n");
  fprintf(stderr, " zs = %e \n", zs);
  fprintf(stderr, " output file name = %s_lens.fits \n\n", file_prefix);
  
  set_distance_lpl_zs(zs);

  num[0] = nx_ext;
  num[1] = ny_ext;
  nn = num[0] * num[1];
  num[2] = 8;
  ntot = num[0] * num[1] * num[2];

  array = (float*)malloc(sizeof(float) * ntot);
  if(array == NULL) terminator("memory allocation failed");

  for(k=0;k<nn;k++){
    ktoxy_ext(k, &tx, &ty);
    
    lensmodel(tx, ty, pout, -1, 0);
    array[k] = pout[0];
    array[k + nn] = pout[1];
    array[k + nn * 2] = pout[2];
    array[k + nn * 3] = pout[3];
    array[k + nn * 4] = pout[4];
    array[k + nn * 5] = pout[5];
    array[k + nn * 6] = pout[6];
    array[k + nn * 7] = pout[7];
  }

  fitsfile *fptr;
  fits_create_file(&fptr, fname, &status);
  fits_create_img(fptr, -32, 3, num, &status);

  fits_write_pix(fptr, TFLOAT, fpixel, ntot,  array, &status);

  fits_header(fptr);
  
  fits_write_key(fptr, TDOUBLE, "ZS", (void*)&zs, "source redshift", &status);
  fits_write_comment(fptr, "cube 1: alpha_x", &status);
  fits_write_comment(fptr, "cube 2: alpha_y", &status);
  fits_write_comment(fptr, "cube 3: time delay", &status);
  fits_write_comment(fptr, "cube 4: kappa", &status);
  fits_write_comment(fptr, "cube 5: gamma1", &status);
  fits_write_comment(fptr, "cube 6: gamma2", &status);
  fits_write_comment(fptr, "cube 7: mu^-1", &status);
  fits_write_comment(fptr, "cube 8: rotation", &status);

  fits_close_file(fptr, &status);

  free(array);

  return;
} 

void fits_header(fitsfile *fptr)
{
  int i, j, ii, jj;
  double tx, ty, x, y, x0, y0;
  int status = 0;

  if(flag_addwcs == 0){
    ktoxy_ext(0, &tx, &ty);
    i = j = 1;
    fits_write_key(fptr, TSTRING, "CTYPE1", (void*)"LINEAR", "", &status);
    fits_write_key(fptr, TSTRING, "CTYPE2", (void*)"LINEAR", "", &status);
    fits_write_key(fptr, TINT, "CRPIX1", (void*)&i, "", &status);
    fits_write_key(fptr, TINT, "CRPIX2", (void*)&j, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL1", (void*)&tx, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", (void*)&ty, "", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT1", (void*)&pix_ext, "", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT2", (void*)&pix_ext, "", &status);
  } else {
    tx = (-1.0) * pix_ext / 3600.0;
    ty = pix_ext / 3600.0;
    i = (int)(((0.0 - xmin) / pix_ext) - 0.5) + 1;
    j = (int)(((0.0 - ymin) / pix_ext) - 0.5) + 1;
    /* ktoxy_ext(xytok_ext(0.0, 0.0), &x, &y); */
    ii = (int)(((0.0 - xmin) / pix_ext) - 0.5);
    jj = (int)(((0.0 - ymin) / pix_ext) - 0.5);
    x = xmin + pix_ext * (ii + 0.5);
    y = ymin + pix_ext * (jj + 0.5);
    x0 = wcs_ra0 - x / (3600.0 * cos((M_PI / 180.0) * wcs_dec0));
    y0 = wcs_dec0 + y / 3600.0;
    fits_write_key(fptr, TSTRING, "CTYPE1", (void*)"RA---CAR", "", &status);
    fits_write_key(fptr, TSTRING, "CTYPE2", (void*)"DEC--CAR", "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL1", (void*)&x0, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", (void*)&y0, "", &status);
    fits_write_key(fptr, TINT, "CRPIX1", (void*)&i, "", &status);
    fits_write_key(fptr, TINT, "CRPIX2", (void*)&j, "", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT1", (void*)&tx, "", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT2", (void*)&ty, "", &status); 
  }

  fits_write_key(fptr, TSTRING, "IN_FILE", (void*)fname_input, "input file name", &status);

  return;
}

/*--------------------------------------------------------------
  write source plane map of magnifications and the number of images 
*/

void writelens_splane(double zs, double sxmin, double sxmax, double symin, double symax, double spix)
{
  long num[3];
  int status = 0;
  long fpixel[] = {1, 1, 1};
  float *array;
  int i, j, k, l, ni, nx, ny, nn, ntot, ii, jj;
  double xs, ys, magtot, magmax, tx, ty, x, y, x0, y0;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  char fname[INPUT_MAXCHAR];

  sprintf(fname, "!%s_lens_splane.fits", file_prefix);

  fprintf(stderr, "######## writing lens properties in the source plane\n");
  fprintf(stderr, " zs = %e  pixel = %e\n", zs, spix);
  fprintf(stderr, " sxmin = %e  sxmax = %e  symin = %e  symax = %e \n", sxmin, sxmax, symin, symax);
  
  if((sxmax <= sxmin) || (symax <= symin) || (spix <= 0))
    terminator("invalid parameter (writelens_splane)");
  
  nx = (int)((sxmax - sxmin + NPIX_SMALL_OFFSET) / spix);
  ny = (int)((symax - symin + NPIX_SMALL_OFFSET) / spix);
  
  if((nx > NMAX_PIXEL) || (ny > NMAX_PIXEL))
    terminator("pixel number exceeds the limit (writelens_splane)");
  
  fprintf(stderr, " nx = %d  ny = %d\n", nx, ny);
  fprintf(stderr, " output file name = %s_lens_splane.fits \n\n", file_prefix);

  num[0] = nx;
  num[1] = ny;
  nn = num[0] * num[1];
  num[2] = 3;
  ntot = nn * num[2];

  array = (float*)malloc(sizeof(float) * ntot);
  if(array == NULL) terminator("memory allocation failed");

  for(k=0;k<nn;k++){
    i = k % nx;
    j = (k - i) / nx;
    xs = sxmin + spix * (i + 0.5);
    ys = symin + spix * (j + 0.5);

    findimg(xs, ys, zs, &ni, rr, 0);

    magtot = 0.0;
    magmax = 0.0;
    for(l=0;l<ni;l++){
      magtot = magtot + fabs(rr[l][2]);
      if(fabs(rr[l][2]) > magmax) magmax = fabs(rr[l][2]);
    }
    
    array[k]          = magtot;
    array[k + nn]     = magmax;
    array[k + 2 * nn] = (double)ni;
  }
     
     fitsfile *fptr;
     fits_create_file(&fptr, fname, &status);

  fits_create_img(fptr, -32, 3, num, &status);
  fits_write_pix(fptr, TFLOAT, fpixel, ntot, array, &status);

  if(flag_addwcs == 0){
    tx = sxmin + spix * 0.5;
    ty = symin + spix * 0.5;
    i = j = 1;
    fits_write_key(fptr, TSTRING, "CTYPE1", (void*)"LINEAR", "", &status);
    fits_write_key(fptr, TSTRING, "CTYPE2", (void*)"LINEAR", "", &status);
    fits_write_key(fptr, TINT, "CRPIX1", (void*)&i, "", &status);
    fits_write_key(fptr, TINT, "CRPIX2", (void*)&j, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL1", (void*)&tx, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", (void*)&ty, "", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT1", (void*)&spix, "", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT2", (void*)&spix, "", &status);
  } else {
    tx = (-1.0) * spix / 3600.0;
    ty = spix / 3600.0;
    i = (int)(((0.0 - sxmin) / spix) - 0.5) + 1;
    j = (int)(((0.0 - symin) / spix) - 0.5) + 1;
    /* ktoxy_ext(xytok_ext(0.0, 0.0), &x, &y); */
    ii = (int)(((0.0 - sxmin) / spix) - 0.5);
    jj = (int)(((0.0 - symin) / spix) - 0.5);
    x = sxmin + spix * (ii + 0.5);
    y = symin + spix * (jj + 0.5);
    x0 = wcs_ra0 - x / (3600.0 * cos((M_PI / 180.0) * wcs_dec0));
    y0 = wcs_dec0 + y / 3600.0;
    fits_write_key(fptr, TSTRING, "CTYPE1", (void*)"RA---CAR", "", &status);
    fits_write_key(fptr, TSTRING, "CTYPE2", (void*)"DEC--CAR", "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL1", (void*)&x0, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", (void*)&y0, "", &status);
    fits_write_key(fptr, TINT, "CRPIX1", (void*)&i, "", &status);
    fits_write_key(fptr, TINT, "CRPIX2", (void*)&j, "", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT1", (void*)&tx, "", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT2", (void*)&ty, "", &status); 
  }

  fits_write_key(fptr, TSTRING, "IN_FILE", (void*)fname_input, "input file name", &status);
     
  fits_write_key(fptr, TDOUBLE, "ZS", (void*)&zs, "source redshift", &status);
  fits_write_comment(fptr, "cube 1: total manification", &status);
  fits_write_comment(fptr, "cube 2: maximum magnification", &status);
  fits_write_comment(fptr, "cube 3: number of images", &status);

  fits_close_file(fptr, &status);

  free(array);

  return;
} 

/*--------------------------------------------------------------
  write lensed images
*/

void writeimage(double sky, double sigma, int flag_source)
{
  long num[3];
  int status = 0;
  long fpixel[] = {1, 1, 1};
  int i, j, k, nn, ntot;
  double x, p, tx, ty, f;
  char fname[INPUT_MAXCHAR];

  long num2[2];
  long fpixel2[] = {1, 1};

  if(flag_source == 0){
    sprintf(fname, "!%s_image.fits", file_prefix);
    
    fprintf(stderr, "######## writing extended images (lensed)\n");
    fprintf(stderr, " sky = %e,  noise = %e\n", sky, sigma);
    fprintf(stderr, " output file name = %s_image.fits \n\n", file_prefix);
  } else {
    sprintf(fname, "!%s_source.fits", file_prefix);
    
    fprintf(stderr, "######## writing extended images (original)\n");
    fprintf(stderr, " sky = %e,  noise = %e\n", sky, sigma);
    fprintf(stderr, " output file name = %s_source.fits \n\n", file_prefix);
    
  }

  flag_computeall = 1;
  i_ext_fid = -1;
  ext_set_image(0, flag_source, 1);
  
  if(num_ext > 1){ 
    num[0] = nx_ext;
    num[1] = ny_ext;
    nn = num[0] * num[1];
    num[2] = num_ext + 1; 
    ntot = nn * num[2];
  } else {
    num2[0] = nx_ext; 
    num2[1] = ny_ext;
    nn = num2[0] * num2[1];
    ntot = nn;
  }
  
  for(k=0;k<nn;k++){
    ktoxy_ext(k, &tx, &ty);
    
    f = 0.0;
    for(i=0;i<num_ext;i++){
      x = array_ext_img[k + i * nn];
      f = f + x;
      p = calc_pix_noise(x, sky, sigma);
      array_ext_img[k + i * nn] = (float)p;
    }

    if(num_ext > 1){
      p = calc_pix_noise(f, sky, sigma);
      array_ext_img[k + num_ext * nn] = (float)p;
    }
  }
  
  ktoxy_ext(0, &tx, &ty);
  i = j = 1;

  fitsfile *fptr;
  fits_create_file(&fptr, fname, &status);
  if(num_ext > 1){ 
    fits_create_img(fptr, -32, 3, num, &status);
    fits_write_pix(fptr, TFLOAT, fpixel, ntot, array_ext_img, &status);
  } else {
    fits_create_img(fptr, -32, 2, num2, &status);
    fits_write_pix(fptr, TFLOAT, fpixel2, ntot, array_ext_img, &status);
  }
 
  fits_header(fptr);

  if(num_ext > 1){
    fits_write_comment(fptr, "cube 1--N: individual extend images", &status);
    fits_write_comment(fptr, "cube N+1 : all extend images", &status);
  }
  
  fits_close_file(fptr, &status);

  return;
}

double calc_pix_noise(double flux, double sky, double sigma)
{
  double s, p;

  if(sigma > 0.0){
    /* use Gaussian distribution */
    s = sqrt(fabs(flux) / (obs_gain * ((double)obs_ncomb)) + sigma * sigma);
    p = flux + sky + gsl_ran_gaussian(ran_gsl, 1.0) * s;
  } else if(sigma < 0.0){
    /* use Poisson distribution */
    s = ((double)obs_ncomb) * ((flux + sky) * obs_gain + obs_readnoise * obs_readnoise);
    if(s < FLOOR_POISSON) s = FLOOR_POISSON;
    p = flux + sky + ((double)(gsl_ran_poisson(ran_gsl, s)) - s) / (obs_gain * ((double)obs_ncomb));
  } else {
    s = 0.0;
    p = flux + sky;
  }

  return p;
}

/*--------------------------------------------------------------
  write lensed images, both point and extend
*/

void writeimageall(double sky, double sigma, int flag_source)
{
  long num[3];
  int status = 0;
  long fpixel[] = {1, 1, 1};
  int i, j, k, ns, nn, ntot;
  double x, p, tx, ty, f;
  char fname[INPUT_MAXCHAR];

  long num2[2];
  long fpixel2[] = {1, 1};

  if(flag_source == 0){
    sprintf(fname, "!%s_image.fits", file_prefix);
    
    fprintf(stderr, "######## writing point and extended images (lensed)\n");
    fprintf(stderr, " sky = %e,  noise = %e\n", sky, sigma);
    fprintf(stderr, " output file name = %s_image.fits \n\n", file_prefix);
  } else {
    sprintf(fname, "!%s_source.fits", file_prefix);
    
    fprintf(stderr, "######## writing point and extended images (original)\n");
    fprintf(stderr, " sky = %e,  noise = %e\n", sky, sigma);
    fprintf(stderr, " output file name = %s_source.fits \n\n", file_prefix);
    
  }

  ns = num_ext + num_poi;
  if(ns < 1) terminator("no source included");

  if(ns == 1){
    num2[0] = nx_ext; 
    num2[1] = ny_ext;
    nn = num2[0] * num2[1];
    ntot = nn;
  } else {
    num[0] = nx_ext;
    num[1] = ny_ext;
    nn = num[0] * num[1];
    num[2] = ns + 1; 
    ntot = nn * num[2];
  }

  array_imageall = (float*)malloc(sizeof(float) * ntot);
  if(array_imageall == NULL) terminator("memory allocation failed");

  if(num_poi > 0){
    poimg_set_table(flag_source);
  }

  if(num_ext > 0){
    flag_computeall = 1;
    i_ext_fid = -1;
    ext_set_image(0, flag_source, 1);
  }
  
  for(k=0;k<nn;k++){

    f = 0.0;
    for(i=0;i<num_poi;i++){
      x = array_poimg[k + i * nn];
      f = f + x;
      p = calc_pix_noise(x, sky, sigma);
      array_imageall[k + i * nn] = (float)p;
    }
    for(i=0;i<num_ext;i++){
      j = num_poi + i;
      x = array_ext_img[k + i * nn];
      f = f + x;
      p = calc_pix_noise(x, sky, sigma);
      array_imageall[k + j * nn] = (float)p;
    }
    
    if(ns > 1){
      p = calc_pix_noise(f, sky, sigma);
      array_imageall[k + ns * nn] = (float)p;
    }
  }
  
  ktoxy_ext(0, &tx, &ty);
  i = j = 1;

  fitsfile *fptr;
  fits_create_file(&fptr, fname, &status);
  if(ns > 1){ 
    fits_create_img(fptr, -32, 3, num, &status);
    fits_write_pix(fptr, TFLOAT, fpixel, ntot, array_imageall, &status);
  } else {
    fits_create_img(fptr, -32, 2, num2, &status);
    fits_write_pix(fptr, TFLOAT, fpixel2, ntot, array_imageall, &status);
  }
 
  fits_header(fptr);

  if(ns > 1){
    fits_write_comment(fptr, "cube 1--N: individual point and extend images", &status);
    fits_write_comment(fptr, "cube N+1 : all point and extend images", &status);
  }

  fits_close_file(fptr, &status);
  free(array_imageall);

  return;
}

/*--------------------------------------------------------------
  write time delay surface
*/

void writetd_ext(void)
{
  long num[3];
  int status = 0;
  long fpixel[] = {1, 1, 1};
  float *array;
  int ii, k, nn, ntot;
  double xs, ys, tx, ty, pout[NPAR_LMODEL];
  char fname[INPUT_MAXCHAR];

  long num2[2];
  long fpixel2[] = {1, 1};

  sprintf(fname, "!%s_td_extend.fits", file_prefix);

  fprintf(stderr, "######## writing time delay surfaces for extended sources\n");
  fprintf(stderr, " output file name = %s_td_extend.fits \n\n", file_prefix);
  
  if(num_ext > 1){ 
    num[0] = nx_ext;
    num[1] = ny_ext;
    nn = num[0] * num[1];
    num[2] = num_ext; 
    ntot = nn * num[2];
  } else {
    num2[0] = nx_ext; 
    num2[1] = ny_ext;
    nn = num2[0] * num2[1];
    ntot = nn;
  }

  array = (float*)malloc(sizeof(float) * ntot);
  if(array == NULL) terminator("memory allocation failed");

  for(ii=0;ii<num_ext;ii++){
    set_distance_lpl_zs(para_ext[ii][0]);
    if(nlp > 1) terminator("writetd available only for a single lens plane");
   
    xs = para_ext[ii][2];
    ys = para_ext[ii][3];

    for(k=0;k<nn;k++){
      ktoxy_ext(k, &tx, &ty);
      
      lensmodel(tx, ty, pout, -1, 0);
      
      array[k + ii * nn] = pout[2] - tdelay_fac(zl_ext, dis_os, dis_ol, dis_ls) * 0.5 * (pout[0] * pout[0] + pout[1] * pout[1] - (tx - xs) * (tx - xs) - (ty - ys) * (ty - ys));
    }
  }
  
  fitsfile *fptr;
  fits_create_file(&fptr, fname, &status);

  if(num_ext > 1){ 
    fits_create_img(fptr, -32, 3, num, &status);
    fits_write_pix(fptr, TFLOAT, fpixel, ntot, array, &status);
  } else {
    fits_create_img(fptr, -32, 2, num2, &status);
    fits_write_pix(fptr, TFLOAT, fpixel2, ntot, array, &status);
  }

  fits_header(fptr);

  fits_close_file(fptr, &status);

  free(array);

  return;
} 

void writetd_poi(void)
{
  long num[3];
  int status = 0;
  long fpixel[] = {1, 1, 1};
  float *array;
  int ii, k, nn, ntot;
  double xs, ys, tx, ty, pout[NPAR_LMODEL];
  char fname[INPUT_MAXCHAR];

  long num2[2];
  long fpixel2[] = {1, 1};

  sprintf(fname, "!%s_td_point.fits", file_prefix);

  fprintf(stderr, "######## writing time delay surfaces for point sources\n");
  fprintf(stderr, " output file name = %s_td_point.fits \n\n", file_prefix);
  
  if(num_poi > 1){ 
    num[0] = nx_ext;
    num[1] = ny_ext;
    nn = num[0] * num[1];
    num[2] = num_poi; 
    ntot = nn * num[2];
  } else {
    num2[0] = nx_ext; 
    num2[1] = ny_ext;
    nn = num2[0] * num2[1];
    ntot = nn;
  }

  array = (float*)malloc(sizeof(float) * ntot);
  if(array == NULL) terminator("memory allocation failed");

  for(ii=0;ii<num_poi;ii++){
    set_distance_lpl_zs(para_poi[ii][0]);
    if(nlp > 1) terminator("writetd available only for a single lens plane");
			   
    xs = para_poi[ii][1];
    ys = para_poi[ii][2];

    for(k=0;k<nn;k++){
      ktoxy_ext(k, &tx, &ty);
      
      lensmodel(tx, ty, pout, -1, 0);
      
      array[k + ii * nn] = pout[2] - tdelay_fac(zl_ext, dis_os, dis_ol, dis_ls) * 0.5 * (pout[0] * pout[0] + pout[1] * pout[1] - (tx - xs) * (tx - xs) - (ty - ys) * (ty - ys));
    }
  }
  
  fitsfile *fptr;
  fits_create_file(&fptr, fname, &status);

  if(num_poi > 1){ 
    fits_create_img(fptr, -32, 3, num, &status);
    fits_write_pix(fptr, TFLOAT, fpixel, ntot, array, &status);
  } else {
    fits_create_img(fptr, -32, 2, num2, &status);
    fits_write_pix(fptr, TFLOAT, fpixel2, ntot, array, &status);
  }

  fits_header(fptr);

  fits_close_file(fptr, &status);

  free(array);

  return;
} 

/*--------------------------------------------------------------
  read obsimage (extend)
*/

void readobs_extend(char *infile, int verb)
{
  long num[2], num2[2];
  int status = 0;
  int nkeys, it, naxis, anynul;
  long fpixel[] = {1, 1};
  int nn;

  if(verb == 1){
    fprintf(stderr, "######## reading obs fits file for extend\n");
    fprintf(stderr, " input file name = %s \n\n", infile);
  }
  
  if(flag_arrayobs == 0){
    flag_arrayobs = 1;

    num[0] = nx_ext;
    num[1] = ny_ext;
    nn = num[0] * num[1];
    
    array_obs = (float*)malloc(sizeof(float) * nn);
    if(array_obs == NULL){
      flag_arrayobs = 0;
      terminator("memory allocation failed");
    }

    fitsfile *fptr;
    fits_open_file(&fptr, infile, READONLY, &status);
    if(status > 100) terminator("failed at open fits file");

    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    fits_get_img_type(fptr, &it, &status);
    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, 2, num2, &status);
    
    if(it != -32) terminator("input obs fits must be float");

    if((naxis != 2) || (num[0] != num2[0]) || (num[1] != num2[1]))
      terminator("invalid obs fits size");

    fits_read_pix(fptr, TFLOAT, fpixel, nn, NULL, array_obs, &anynul, &status);

    fits_close_file(fptr, &status);
  } else {
    if(verb == 1){
      fprintf(stderr, " [skipped because the obs file already read]\n\n");
    }
  }

  return;
}

/*--------------------------------------------------------------
  read noise image
*/

void readnoise_extend(char *infile, int verb)
{
  long num[2], num2[2];
  int status = 0;
  int nkeys, it, naxis, anynul;
  long fpixel[] = {1, 1};
  int nn;
  
  if(verb == 1){
     fprintf(stderr, "######## reading noise fits file for extend\n");
     fprintf(stderr, " input file name = %s \n\n", infile);
   }
   
  if(flag_obssig == 0){
    flag_obssig = 1;
   
    num[0] = nx_ext;
    num[1] = ny_ext;
    nn = num[0] * num[1];
    
    array_obsnoise_file = (float*)malloc(sizeof(float) * nn);
    if(array_obsnoise_file == NULL){
      flag_obssig = 0;
      terminator("memory allocation failed");
    }
    
    fitsfile *fptr;
    fits_open_file(&fptr, infile, READONLY, &status);
    if(status > 100) terminator("failed at open fits file");
    
    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    fits_get_img_type(fptr, &it, &status);
    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, 2, num2, &status);
    
    if(it != -32) terminator("input noise fits must be float");
    
    if((naxis != 2) || (num[0] != num2[0]) || (num[1] != num2[1]))
      terminator("invalid noise fits size");
   
    fits_read_pix(fptr, TFLOAT, fpixel, nn, NULL, array_obsnoise_file, &anynul, &status);
    
    fits_close_file(fptr, &status);
 
  } else {
    if(verb == 1){
      fprintf(stderr, " [skipped because the noise file already read]\n\n");
    }
  }
 
  return;
}

void readmask(char *infile, int verb)
{
  long num[2], num2[2];
  int status = 0;
  int nkeys, it, naxis, anynul;
  long fpixel[] = {1, 1};
  int nn;
  
  if(verb == 1){
    fprintf(stderr, "######## reading mask file for extend\n");
    fprintf(stderr, " input file name = %s \n\n", infile);
  }
  
  if(flag_obsmask == 0){
    flag_obsmask = 1;

    num[0] = nx_ext;
    num[1] = ny_ext;
    nn = num[0] * num[1];

    array_obs_mask = (int*)malloc(sizeof(int) * nn);

    if(array_obs_mask == NULL){
      flag_obsmask = 0;
      terminator("memory allocation failed");
    }

    fitsfile *fptr;
    fits_open_file(&fptr, infile, READONLY, &status);
    if(status > 100) terminator("failed at open fits file");

    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    fits_get_img_type(fptr, &it, &status);
    fits_get_img_dim(fptr, &naxis, &status);
    fits_get_img_size(fptr, 2, num2, &status);
 
    if(it != 32) terminator("input mask fits must be int");

    if((naxis != 2) || (num[0] != num2[0]) || (num[1] != num2[1]))
      terminator("invalid mask fits size");
    
    fits_read_pix(fptr, TINT, fpixel, nn, NULL, array_obs_mask, &anynul, &status);
    
    fits_close_file(fptr, &status);
  } else {
    if(verb == 1){
      fprintf(stderr, " [skipped because the mask file already read]\n\n");
    }
  }
  
  return;
}

void calc_obsnoise(void)
{
  int i, j, k, ko, n, nn;
  float *atmp;
  double med, sig, med2, sig2, x, hh;
  
  /* hh-sigma clip */
  hh = noise_clip;

  /* first compute sky noise from the data */

  nn = nx_ext * ny_ext;
  atmp = (float*)malloc(sizeof(float) * nn);
  if(atmp == NULL) terminator("memory allocation failed");

  n = 0;
  med = 0.0;
  sig = 0.0;
  for(i=0;i<nn;i++){
    j = 1;
    if(flag_obsmask == 1){
      if(array_obs_mask[i] > 0) j = 0;
    }

    if(j == 1){
      atmp[n] = array_obs[i];
      med = med + atmp[n];
      sig = sig + atmp[n] * atmp[n];
      n++;
    }
  }
  
  med = med / ((double)n);
  sig = sqrt(sig / ((double)n) - med * med);

  k = n;
  do{
    ko = k;
    k = 0;
    med2 = 0.0;
    sig2 = 0.0;
    for(j=0;j<n;j++){
      if((atmp[j] > (med - hh * sig)) && (atmp[j] < (med + hh * sig))){
	k++;
	med2 = med2 + atmp[j];
	sig2 = sig2 + atmp[j] * atmp[j];
      }
    }
    med = med2 / ((double)k);
    sig = sqrt(sig2 / ((double)k) - med * med);
  }while(k < ko);

  free(atmp);
  
  if(skyfix_value > (0.1 * DEF_SKYFIX_VALUE)) skyfix_value = med;
  skymed = skyfix_value;
  skysigma = sig;

  /* start making noise array */

  array_obs_noise = (float*)malloc(sizeof(float) * nn);
  if(array_obs_noise == NULL) terminator("memory allocation failed");
  
  for(i=0;i<nn;i++){
    j = 1;
    if(flag_obsmask == 1){
      if(array_obs_mask[i] > 0) j = 0;
    }
    
    if(j == 1){
      x = fabs(array_obs[i] - med) / (sig * sig * (obs_gain * ((float)obs_ncomb)));
      array_obs_noise[i] = sig * sqrt(1.0 + x);
    } else {
      array_obs_noise[i] = 0.0;
    }
  }

  return;
}

void obs_unset_table(void)
{
  if(flag_arrayobs == 1){
    free(array_obs);
    free(array_obs_noise);
    flag_arrayobs = 0;
  }

  if(flag_obsmask == 1){
    free(array_obs_mask);
    flag_obsmask = 0;
  }

  if(flag_obssig == 1){
    free(array_obsnoise_file);
    flag_obssig = 0;
  }

  return;
}

/*--------------------------------------------------------------
  write obs noise file
*/

void writenoise(void)
{
  long num[2];
  int status = 0;
  long fpixel[] = {1, 1};
  int nn;
  char fname[INPUT_MAXCHAR];

  fprintf(stderr, "######## writing obs noise fits file\n");

  if(flag_arrayobs == 1){

    sprintf(fname, "!%s_obsnoise.fits", file_prefix);

    
    fprintf(stderr, " output file name = %s_obsnoise.fits \n\n", file_prefix);
    
    num[0] = nx_ext;
    num[1] = ny_ext;
    nn = num[0] * num[1];
    
    fitsfile *fptr;
    fits_create_file(&fptr, fname, &status);
    fits_create_img(fptr, -32, 2, num, &status);
    
    fits_write_pix(fptr, TFLOAT, fpixel, nn, array_obs_noise, &status);
    
    fits_header(fptr);
    
    fits_close_file(fptr, &status);
  } else {
    fprintf(stderr, "\n");
  }

  return;
}

/*--------------------------------------------------------------
  add noise to obs file
*/

void addnoise(double sky, double sigma)
{
  long num[2];
  int status = 0;
  long fpixel[] = {1, 1};
  int i, j, nn;
  double tx, ty, x, p;
  char fname[INPUT_MAXCHAR];
  float *atmp;
 
  fprintf(stderr, "######## add noise to obs file\n");
  fprintf(stderr, " sky = %e,   noise = %e\n", sky, sigma);

  if(flag_arrayobs == 1){

    sprintf(fname, "!%s_addnoise.fits", file_prefix);

    fprintf(stderr, " output file name = %s_addnoise.fits \n\n", file_prefix);
    
    num[0] = nx_ext;
    num[1] = ny_ext;
    nn = num[0] * num[1];

    atmp = (float*)malloc(sizeof(float) * nn);
    if(atmp == NULL) terminator("memory allocation failed");

    for(i=0;i<nn;i++){
      x = array_obs[i];
      p = calc_pix_noise(x, sky, sigma);
      atmp[i] = (float)p;
    }

    ktoxy_ext(0, &tx, &ty);
    i = j = 1;
   
    fitsfile *fptr;
    fits_create_file(&fptr, fname, &status);
    fits_create_img(fptr, -32, 2, num, &status);
    
    fits_write_pix(fptr, TFLOAT, fpixel, nn, atmp, &status);
    
    fits_header(fptr);
     
    fits_close_file(fptr, &status);
  } else {
    fprintf(stderr, "\n");
  }

  return;
}

/*--------------------------------------------------------------
  dump PSF fits file
*/

void writepsf(void)
{
  long num[2];
  long fpixel[] = {1, 1};
  int i, j, nk, nnk, ntot, status = 0;
  double pix_psf;
  char fname[INPUT_MAXCHAR];
  double *gau;
  static float *bout;
  
  sprintf(fname, "!%s_psf.fits", file_prefix);

  fprintf(stderr, "######## writing PSF image\n"); 
  fprintf(stderr, " output file name = %s\n\n", fname); 

  if(flag_seeing == 0) terminator("no PSF set");

  nk = calc_psfnk();
  nnk = nk * 2 + 1;

  bout = (float*)malloc(sizeof(float) * nnk * nnk);
  if(bout == NULL) terminator("memory allocation failed");
  
  gau = (double*)malloc(sizeof(double) * nnk * nnk);
  if(gau == NULL) terminator("memory allocation failed");
  
  pix_psf = calc_pixpsf();
  set_psfimage(gau, nk, pix_psf);

  for(i=0;i<nnk;i++){
    for(j=0;j<nnk;j++){
      bout[i + j * nnk] = gau[i + j * nnk];
    }
  }

  num[0] = nnk;
  num[1] = nnk;
  ntot = num[0] * num[1];

  fitsfile *fptr;
  fits_create_file(&fptr, fname, &status);
  fits_create_img(fptr, -32, 2, num, &status);

  fits_write_pix(fptr, TFLOAT, fpixel, ntot, bout, &status);

  fits_write_key(fptr, TSTRING, "CTYPE1", (void*)"LINEAR", "", &status);
  fits_write_key(fptr, TSTRING, "CTYPE2", (void*)"LINEAR", "", &status);

  fits_close_file(fptr, &status);
  free(bout);

  return;
}

/*--------------------------------------------------------------
  read PSF fits file
*/

void read_psffits(char *infile, int verb)
{
  long num[2];
  int status = 0;
  int nkeys, it, naxis, anynul, i, j, fpsf_xsize, fpsf_ysize, nnfpsf;
  long fpixel[] = {1, 1};
  double f;
  fitsfile *fptr;

  if(verb == 1){
    fprintf(stderr, "######## reading PSF file for extend\n");
    fprintf(stderr, " input file name = %s \n\n", infile);
  }
  
  fits_open_file(&fptr, infile, READONLY, &status);
  if(status > 100) terminator("failed at open PSF fits file");

  fits_get_hdrspace(fptr, &nkeys, NULL, &status);
  fits_get_img_type(fptr, &it, &status);
  if(it != -32) terminator("input PSF fits must be float");

  fits_get_img_dim(fptr, &naxis, &status);
  if(naxis != 2) terminator("input PSF fits file has wrong dimension");
  
  fits_get_img_size(fptr, 2, num, &status);

  fpsf_xsize = num[0];
  fpsf_ysize = num[1];
  nnfpsf = fpsf_xsize * fpsf_ysize;

  if(fpsf_xsize != fpsf_ysize) terminator("input PSF fits must be square");
  if((fpsf_xsize%2) == 0) terminator("input PSF fits pixel number must be odd");

  array_fpsf = (float*)malloc(sizeof(float) * nnfpsf);
  if(array_fpsf == NULL) terminator("memory allocation failed");

  fits_read_pix(fptr, TFLOAT, fpixel, nnfpsf, NULL, array_fpsf, &anynul, &status);
  fits_close_file(fptr, &status);

  fpsf_nk = (fpsf_xsize - 1) / 2;
  flag_seeing = -1;

  f = 0.0;
  for(j=0;j<fpsf_ysize;j++){
    for(i=0;i<fpsf_xsize;i++){
      f = f + array_fpsf[i + j * fpsf_xsize];
    }
  }
  
  for(i=0;i<nnfpsf;i++) array_fpsf[i] = array_fpsf[i] / f;
  
  return;
}
