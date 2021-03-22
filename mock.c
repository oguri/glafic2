#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "glafic.h"

/*--------------------------------------------------------------
  mock image samples
*/

void mock1(int n, double zs, double x1, double x2, double y1, double y2)
{
  int i, j, k, ni;
  int num[NMAX_POIMG_MOCK+1];
  double x, y, area;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  sprintf(fname, "%s_mock.dat", file_prefix);
  
  fprintf(stderr, "######## generating mock sample\n");
  fprintf(stderr, " n_src = %d  zs = %e \n", n, zs);
  fprintf(stderr, " region = rectangle [ %e < x < %e  %e < y < %e ]\n", x1, x2, y1, y2); 
  fprintf(stderr, " output file name = %s \n\n", fname);

  if((x2 < x1) || (y2 < y1) || (n <= 0))
    terminator("invalid range (mock1)");

  area = (x2 - x1) * (y2 - y1);
  for(j=0;j<=NMAX_POIMG_MOCK;j++) num[j] = 0;

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (mock1)");
  
  fprintf(fptr, "# %d %e %e\n", n, zs, area);
  
  k = 1;
  for(i=0;i<n;i++){
    if((i + 1) >= ((k * n) / NUM_MIDREPORT_MOCK)){ 
      fprintf(stderr, "%2d/%2d done\n", k, NUM_MIDREPORT_MOCK);
      k++;
    }
    x = x1 + gsl_rng_uniform(ran_gsl) * (x2 - x1);
    y = y1 + gsl_rng_uniform(ran_gsl) * (y2 - y1);
    findimg(x, y, zs, &ni, rr, 0);
    if(ni <= NMAX_POIMG_MOCK) num[ni]++;
    if(outformat_exp == 0){
      fprintf(fptr, "%d %8.4f %9.4f %9.4f\n", ni, zs, x, y);
    } else{
      fprintf(fptr, "%d %13e %13e %13e\n", ni, zs, x, y);
    }
    for(j=0;j<ni;j++){
      if(outformat_exp == 0){
	fprintf(fptr, "%9.4f %9.4f %9.4f %9.3f\n", rr[j][0], rr[j][1], rr[j][2], rr[j][3]);
      } else {
	fprintf(fptr, "%13e %13e %13e %13e\n", rr[j][0], rr[j][1], rr[j][2], rr[j][3]);
      }
    }
  }
  fprintf(stderr, "\n");

  printf("# Ni   Nevent crosssection\n");
  for(j=0;j<=NMAX_POIMG_MOCK;j++){
    printf("  %2d %8d %e\n", j, num[j], area * ((double)num[j]) / ((double)n));
  }
  fprintf(stderr, "\n");

  fclose(fptr); 

  return;
}

void mock2(int n, double zs, double rmax, double x0, double y0)
{
  int i, j, k, ni;
  int num[NMAX_POIMG_MOCK + 1];
  double r, sr, tt, x, y, area;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  sprintf(fname, "%s_mock.dat", file_prefix);

  fprintf(stderr, "######## generating mock sample\n");
  fprintf(stderr, " n_src = %d  zs = %e \n", n, zs);
  fprintf(stderr, " region = circle [ center = (%e, %e)  r < %e  ]\n", x0, y0, rmax); 
  fprintf(stderr, " output file name = %s \n\n", fname);

  if((rmax <= 0.0) || (n <= 0)) terminator("invalid range (mock2)");

  area = M_PI * rmax * rmax;
  for(j=0;j<=NMAX_POIMG_MOCK;j++) num[j] = 0;

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (mock1)");
  
  fprintf(fptr, "# %d %e %e\n", n, zs, area);
  
  k = 1;
  for(i=0;i<n;i++){
    if((i + 1) >= ((k * n) / NUM_MIDREPORT_MOCK)){ 
      fprintf(stderr, "%2d/%2d done\n", k, NUM_MIDREPORT_MOCK);
      k++;
    }
    r = gsl_rng_uniform(ran_gsl);
    sr = rmax * sqrt(r);
    tt = 2.0 * M_PI * gsl_rng_uniform(ran_gsl);
    x = sr * cos(tt) + x0;
    y = sr * sin(tt) + y0;
    findimg(x, y, zs, &ni, rr, 0);
    if(ni <= NMAX_POIMG_MOCK) num[ni]++;
    if(outformat_exp == 0){
      fprintf(fptr, "%d %8.4f %9.4f %9.4f\n", ni, zs, x, y);
    } else{
      fprintf(fptr, "%d %13e %13e %13e\n", ni, zs, x, y);
    }
    for(j=0;j<ni;j++){
      if(outformat_exp == 0){
	fprintf(fptr, "%9.4f %9.4f %9.4f %9.3f\n", rr[j][0], rr[j][1], rr[j][2], rr[j][3]);
      } else {
	fprintf(fptr, "%13e %13e %13e %13e\n", rr[j][0], rr[j][1], rr[j][2], rr[j][3]);
      }
    }
  }
  fprintf(stderr, "\n");

  printf("# Ni   Nevent crosssection\n");
  for(j=0;j<=NMAX_POIMG_MOCK;j++){
    printf("  %2d %8d %e\n", j, num[j], area * ((double)num[j]) / ((double)n));
  }
  fprintf(stderr, "\n");

  fclose(fptr); 
    
  return;
}

void mock3(int n, double zs, double fac)
{
  int i, j, k, ni;
  int num[NMAX_POIMG_MOCK + 1];
  double x1, x2, y1, y2, x, y, area;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  double causize[NPAR_CAUSIZE];
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  sprintf(fname, "%s_mock.dat", file_prefix);
  
  fprintf(stderr, "######## generating mock sample\n");
  fprintf(stderr, " n_src = %d  zs = %e  size factor = %e\n", n, zs, fac);

  writecrit(zs, causize, 1, 0);

  x1 = 0.5 * (causize[1] + causize[0]) - 0.5 * fac * (causize[1] - causize[0]);
  x2 = 0.5 * (causize[1] + causize[0]) + 0.5 * fac * (causize[1] - causize[0]);
  y1 = 0.5 * (causize[3] + causize[2]) - 0.5 * fac * (causize[3] - causize[2]);
  y2 = 0.5 * (causize[3] + causize[2]) + 0.5 * fac * (causize[3] - causize[2]);

  fprintf(stderr, " region = rectangle [ %e < x < %e  %e < y < %e ]\n", x1, x2, y1, y2); 
  fprintf(stderr, " output file name = %s \n\n", fname);

  if((x2 < x1) || (y2 < y1) || (n <= 0))
    terminator("invalid range (mock3)");

  area = (x2 - x1) * (y2 - y1);
  for(j=0;j<=NMAX_POIMG_MOCK;j++) num[j] = 0;

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (mock3)");
  
  fprintf(fptr, "# %d %e %e\n", n, zs, area);
  
  k = 1;
  for(i=0;i<n;i++){
    if((i + 1) >= ((k * n) / NUM_MIDREPORT_MOCK)){ 
      fprintf(stderr, "%2d/%2d done\n", k, NUM_MIDREPORT_MOCK);
      k++;
    }
    x = x1 + gsl_rng_uniform(ran_gsl) * (x2 - x1);
    y = y1 + gsl_rng_uniform(ran_gsl) * (y2 - y1);
    findimg(x, y, zs, &ni, rr, 0);
    if(ni <= NMAX_POIMG_MOCK) num[ni]++;
    if(outformat_exp == 0){
      fprintf(fptr, "%d %8.4f %9.4f %9.4f\n", ni, zs, x, y);
    } else{
      fprintf(fptr, "%d %13e %13e %13e\n", ni, zs, x, y);
    }
    for(j=0;j<ni;j++){
      if(outformat_exp == 0){
	fprintf(fptr, "%9.4f %9.4f %9.4f %9.3f\n", rr[j][0], rr[j][1], rr[j][2], rr[j][3]);
      } else {
	fprintf(fptr, "%13e %13e %13e %13e\n", rr[j][0], rr[j][1], rr[j][2], rr[j][3]);
      }
    }
  }
  fprintf(stderr, "\n");

  printf("# Ni   Nevent crosssection\n");
  for(j=0;j<=NMAX_POIMG_MOCK;j++){
    printf("  %2d %8d %e\n", j, num[j], area * ((double)num[j]) / ((double)n));
  }
  fprintf(stderr, "\n");

  fclose(fptr); 

  return;
}

void mockline(int n, double zs, double x1, double x2, double y1, double y2, int flag_full)
{
  int i, j, ni;
  double x, y, dx, dy, mag;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  sprintf(fname, "%s_mockline.dat", file_prefix);
  
  fprintf(stderr, "######## generating mock sample\n");
  fprintf(stderr, " n_src = %d  zs = %e \n", n, zs);
  fprintf(stderr, " region = line connecting (%e, %e) and (%e, %e)\n", x1, y1, x2, y2); 
  fprintf(stderr, " output file name = %s \n\n", fname);

  if(n <= 0) terminator("invalid range (mockline)");

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (mockline)");
  
  fprintf(fptr, "# %d %e %e %e %e %e\n", n, zs, x1, x2, y1, y2);

  dx = (x2 - x1) / ((double)n);
  dy = (y2 - y1) / ((double)n);
  
  for(i=0;i<=n;i++){
    x = x1 + ((double)i) * dx;
    y = y1 + ((double)i) * dy;
    findimg(x, y, zs, &ni, rr, 0);
    if(flag_full == 0){
      mag = 0.0;
      for(j=0;j<ni;j++){
	mag = mag + fabs(rr[j][2]);
      }
      if(outformat_exp == 0){
	fprintf(fptr, "%5d %9.4f %9.4f %2d %9.4f\n", i + 1, x, y, ni, mag);
      } else{
	fprintf(fptr, "%5d %13e %13e %2d %13e\n", i + 1, x, y, ni, mag);
      }
    } else {
      if(outformat_exp == 0){
	fprintf(fptr, "%d %8.4f %9.4f %9.4f\n", ni, zs, x, y);
      } else{
	fprintf(fptr, "%d %13e %13e %13e\n", ni, zs, x, y);
      }
      for(j=0;j<ni;j++){
	if(outformat_exp == 0){
	  fprintf(fptr, "%9.4f %9.4f %9.4f %9.3f\n", rr[j][0], rr[j][1], rr[j][2], rr[j][3]);
	} else {
	  fprintf(fptr, "%13e %13e %13e %13e\n", rr[j][0], rr[j][1], rr[j][2], rr[j][3]);
	}
      }
    }
  }

  fclose(fptr); 

  return;
}

/*--------------------------------------------------------------
  mock samples for extended sources
*/

void mockext1(int n, int id, double x1, double x2, double y1, double y2, double e1, double e2, double sbth, double lwlim)
{
  int i, k, narc;
  double x, y, e, t, xx0, yy0, ee0, tt0, area, lwmax;
  double fl, pe, px, py, ar, ff, sb;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;
  char fname2[INPUT_MAXCHAR];
  FILE* fptr_a;

  sprintf(fname, "%s_mockext.dat", file_prefix);
  sprintf(fname2, "%s_mockext_arc.dat", file_prefix);
  
  fprintf(stderr, "######## generating mock sample (extend source)\n");
  fprintf(stderr, " n_src = %d  id = %d \n", n, id);
  fprintf(stderr, " region = rectangle [ %e < x < %e  %e < y < %e ]\n", x1, x2, y1, y2); 
  fprintf(stderr, " ellipticity between %e and  %e\n", e1, e2); 
  fprintf(stderr, " output file name = %s \n", fname);
  fprintf(stderr, " output file name = %s \n\n", fname2);

  if((x2 < x1) || (y2 < y1) || (e2 < e1) || (n <= 0) || (id <= 0) || (id > num_ext))
    terminator("invalid range (mockext1)");

  area = (x2 - x1) * (y2 - y1);

  xx0 = para_ext[id - 1][2];
  yy0 = para_ext[id - 1][3];
  ee0 = para_ext[id - 1][4];
  tt0 = para_ext[id - 1][5];

  i_ext_fid = -1;
  ext_set_table_all(id - 1);
  
  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (mockext1)");
  fptr_a = fopen(fname2, "w");
  if(fptr_a == NULL) terminator("failed at fopen (mockext1)");
  
  fprintf(fptr, "# %d %e %-7s %e %e %e %e %e\n", n, area, inttoemodel(model_ext[id - 1]), para_ext[id - 1][0], para_ext[id - 1][1], para_ext[id - 1][6], para_ext[id - 1][7], lwlim);
  
  k = 1;
  for(i=0;i<n;i++){
    if((i + 1) >= ((k * n) / NUM_MIDREPORT_MOCK)){ 
      fprintf(stderr, "%2d/%2d done\n", k, NUM_MIDREPORT_MOCK);
      k++;
    }
    x = x1 + gsl_rng_uniform(ran_gsl) * (x2 - x1);
    y = y1 + gsl_rng_uniform(ran_gsl) * (y2 - y1);
    e = e1 + gsl_rng_uniform(ran_gsl) * (e2 - e1);
    t = (-180.0) + gsl_rng_uniform(ran_gsl) * 360.0;
    
    para_ext[id - 1][2] = x;
    para_ext[id - 1][3] = y;
    para_ext[id - 1][4] = e;
    para_ext[id - 1][5] = t;

    sb = sbth;

    /* image */
    ext_set_image(id, 0, 0);
    ext_est_image_i(id, lwlim, &sb, &fl, &pe, &px, &py, &ar, &ff, &lwmax, &narc, fptr_a);
    
    fprintf(fptr, "%13e %13e %13e %13e %13e %13e %13e %13e %13e %13e %13e %13e %d\n", x, y, e, t, fl, pe, px, py, sb, ar, ff, lwmax, narc);
  }
  fprintf(stderr, "\n");


  fclose(fptr); 
  fclose(fptr_a); 

  para_ext[id - 1][2] = xx0;
  para_ext[id - 1][3] = yy0;
  para_ext[id - 1][4] = ee0;
  para_ext[id - 1][5] = tt0;

  return;
}

void mockext2(int n, int id,  double rmax, double x0, double y0, double e1, double e2, double sbth, double lwlim)
{
  int i, k, narc;
  double x, y, e, t, r, sr, tt, xx0, yy0, ee0, tt0, area, lwmax;
  double fl, pe, px, py, ar, ff, sb;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;
  char fname2[INPUT_MAXCHAR];
  FILE* fptr_a;

  sprintf(fname, "%s_mockext.dat", file_prefix);
  sprintf(fname2, "%s_mockext_arc.dat", file_prefix);
  
  fprintf(stderr, "######## generating mock sample (extend source)\n");
  fprintf(stderr, " n_src = %d  id = %d \n", n, id);
  fprintf(stderr, " region = circle [ center = (%e, %e)  r < %e  ]\n", x0, y0, rmax); 
  fprintf(stderr, " ellipticity between %e and  %e\n", e1, e2); 
  fprintf(stderr, " output file name = %s \n", fname);
  fprintf(stderr, " output file name = %s \n\n", fname2);

  if((rmax <= 0.0) || (e2 < e1) || (n <= 0) || (id <= 0) || (id > num_ext))
    terminator("invalid range (mockext2)");

  area = M_PI * rmax * rmax;

  xx0 = para_ext[id - 1][2];
  yy0 = para_ext[id - 1][3];
  ee0 = para_ext[id - 1][4];
  tt0 = para_ext[id - 1][5];

  i_ext_fid = -1;
  ext_set_table_all(id - 1);

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (mockext2)");
  fptr_a = fopen(fname2, "w");
  if(fptr_a == NULL) terminator("failed at fopen (mockext2)");
  
  fprintf(fptr, "# %d %e %-7s %e %e %e %e %e\n", n, area, inttoemodel(model_ext[id - 1]), para_ext[id - 1][0], para_ext[id - 1][1], para_ext[id - 1][6], para_ext[id - 1][7], lwlim);
  
  k = 1;
  for(i=0;i<n;i++){
    if((i + 1) >= ((k * n) / NUM_MIDREPORT_MOCK)){ 
      fprintf(stderr, "%2d/%2d done\n", k, NUM_MIDREPORT_MOCK);
      k++;
    }
    r = gsl_rng_uniform(ran_gsl);
    sr = rmax * sqrt(r);
    tt = 2.0 * M_PI * gsl_rng_uniform(ran_gsl);
    x = sr * cos(tt) + x0;
    y = sr * sin(tt) + y0;
    e = e1 + gsl_rng_uniform(ran_gsl) * (e2 - e1);
    t = (-180.0) + gsl_rng_uniform(ran_gsl) * 360.0;
    
    para_ext[id - 1][2] = x;
    para_ext[id - 1][3] = y;
    para_ext[id - 1][4] = e;
    para_ext[id - 1][5] = t;

    sb = sbth;

    /* image */
    ext_set_image(id, 0, 0);
    ext_est_image_i(id, lwlim, &sb, &fl, &pe, &px, &py, &ar, &ff, &lwmax, &narc, fptr_a);
    
    fprintf(fptr, "%13e %13e %13e %13e %13e %13e %13e %13e %13e %13e %13e %13e %d\n", x, y, e, t, fl, pe, px, py, sb, ar, ff, lwmax, narc);
  }
  fprintf(stderr, "\n");


  fclose(fptr); 
  fclose(fptr_a); 

  para_ext[id - 1][2] = xx0;
  para_ext[id - 1][3] = yy0;
  para_ext[id - 1][4] = ee0;
  para_ext[id - 1][5] = tt0;

  return;
}

void mockext3(int n, int id, double fac, double e1, double e2, double sbth, double lwlim)
{
  int i, k, narc;
  double x1, x2, y1, y2, x, y, e, t, xx0, yy0, ee0, tt0, area, lwmax;
  double fl, pe, px, py, ar, ff, sb;
  double causize[NPAR_CAUSIZE];
  char fname[INPUT_MAXCHAR];
  FILE* fptr;
  char fname2[INPUT_MAXCHAR];
  FILE* fptr_a;

  sprintf(fname, "%s_mockext.dat", file_prefix);
  sprintf(fname2, "%s_mockext_arc.dat", file_prefix);
  
  fprintf(stderr, "######## generating mock sample (extend source)\n");
  fprintf(stderr, " n_src = %d  id = %d  size factor = %e\n", n, id, fac);

  writecrit(para_ext[id - 1][0], causize, 1, 0);

  x1 = 0.5 * (causize[1] + causize[0]) - 0.5 * fac * (causize[1] - causize[0]);
  x2 = 0.5 * (causize[1] + causize[0]) + 0.5 * fac * (causize[1] - causize[0]);
  y1 = 0.5 * (causize[3] + causize[2]) - 0.5 * fac * (causize[3] - causize[2]);
  y2 = 0.5 * (causize[3] + causize[2]) + 0.5 * fac * (causize[3] - causize[2]);

  fprintf(stderr, " region = rectangle [ %e < x < %e  %e < y < %e ]\n", x1, x2, y1, y2); 
  fprintf(stderr, " ellipticity between %e and  %e\n", e1, e2); 
  fprintf(stderr, " output file name = %s \n", fname);
  fprintf(stderr, " output file name = %s \n\n", fname2);

  if((x2 < x1) || (y2 < y1) || (e2 < e1) || (n <= 0) || (id <= 0) || (id > num_ext))
    terminator("invalid range (mockext3)");

  area = (x2 - x1) * (y2 - y1);

  xx0 = para_ext[id - 1][2];
  yy0 = para_ext[id - 1][3];
  ee0 = para_ext[id - 1][4];
  tt0 = para_ext[id - 1][5];

  i_ext_fid = -1;
  ext_set_table_all(id - 1);

  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (mockext3)");
  fptr_a = fopen(fname2, "w");
  if(fptr_a == NULL) terminator("failed at fopen (mockext3)");
  
  fprintf(fptr, "# %d %e %-7s %e %e %e %e %e\n", n, area, inttoemodel(model_ext[id - 1]), para_ext[id - 1][0], para_ext[id - 1][1], para_ext[id - 1][6], para_ext[id - 1][7], lwlim);
  
  k = 1;
  for(i=0;i<n;i++){
    if((i + 1) >= ((k * n) / NUM_MIDREPORT_MOCK)){ 
      fprintf(stderr, "%2d/%2d done\n", k, NUM_MIDREPORT_MOCK);
      k++;
    }
    x = x1 + gsl_rng_uniform(ran_gsl) * (x2 - x1);
    y = y1 + gsl_rng_uniform(ran_gsl) * (y2 - y1);
    e = e1 + gsl_rng_uniform(ran_gsl) * (e2 - e1);
    t = (-180.0) + gsl_rng_uniform(ran_gsl) * 360.0;
    
    para_ext[id - 1][2] = x;
    para_ext[id - 1][3] = y;
    para_ext[id - 1][4] = e;
    para_ext[id - 1][5] = t;

    sb = sbth;

    /* image */
    ext_set_image(id, 0, 0);
    ext_est_image_i(id, lwlim, &sb, &fl, &pe, &px, &py, &ar, &ff, &lwmax, &narc, fptr_a);
    
    fprintf(fptr, "%13e %13e %13e %13e %13e %13e %13e %13e %13e %13e %13e %13e %d\n", x, y, e, t, fl, pe, px, py, sb, ar, ff, lwmax, narc);
  }
  fprintf(stderr, "\n");


  fclose(fptr); 
  fclose(fptr_a); 

  para_ext[id - 1][2] = xx0;
  para_ext[id - 1][3] = yy0;
  para_ext[id - 1][4] = ee0;
  para_ext[id - 1][5] = tt0;

  return;
}

