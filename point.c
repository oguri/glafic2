#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

static double zs_point;

/*--------------------------------------------------------------
  compute lensing properties
*/

void calcimg(double x, double y, double zs)
{
  int i;
  double pout[NPAR_LMODEL];

  set_distance_lpl_zs(zs);

  lensmodel(x, y, pout, -1, 0);

  fprintf(stderr, "######## calculating lens property\n");
  fprintf(stderr, " zs = %e  x = %e  y = %e\n\n", zs, x, y);
  for(i=0;i<nlp;i++){
    printf("sigma_crit (zl = %f) = %e [h M_Sun Mpc^-2 (physical density)]\n", zl_lpl[i], sigma_crit_dis(dis_angulard(0.0, zs), dis_angulard(0.0, zl_lpl[i]), dis_angulard(zl_lpl[i], zs)));
  }
  printf("\n");
  if(outformat_exp == 0){
    printf("alpha_x  = %f\n", pout[0]);
    printf("alpha_y  = %f\n", pout[1]);
    printf("kappa    = %f\n", pout[3]);
    printf("gamma_1  = %f\n", pout[4]);
    printf("gamma_2  = %f\n", pout[5]);
    printf("gamma    = %f\n", sqrt(pout[4] * pout[4] + pout[5] * pout[5]));
    printf("mag      = %f\n", 1.0 / (pout[6] + imag_ceil));
    printf("td[day]  = %f\n", pout[2]);
    printf("rotation = %f\n", pout[7]);
    printf("x_src    = %f\n", x - pout[0]);
    printf("y_src    = %f\n", y - pout[1]);
    fprintf(stderr, "\n");
  } else {
    printf("alpha_x  = %13e\n", pout[0]);
    printf("alpha_y  = %13e\n", pout[1]);
    printf("kappa    = %13e\n", pout[3]);
    printf("gamma_1  = %13e\n", pout[4]);
    printf("gamma_2  = %13e\n", pout[5]);
    printf("gamma    = %13e\n", sqrt(pout[4] * pout[4] + pout[5] * pout[5]));
    printf("mag      = %13e\n", 1.0 / (pout[6] + imag_ceil));
    printf("td[day]  = %13e\n", pout[2]);
    printf("rotation = %13e\n", pout[7]);
    printf("x_src    = %13e\n", x - pout[0]);
    printf("y_src    = %13e\n", y - pout[1]);
    fprintf(stderr, "\n");
  }

  return;
}

/*--------------------------------------------------------------
  find source position and other images
*/

void findsrcimg(int i, double x, double y)
{
  double pout[NPAR_LMODEL], xs, ys;

  fprintf(stderr, "######## find source position, and the other images \n\n");
  fprintf(stderr, "point id = %d,  x = %e,  y = %e\n\n", i, x, y);

  if((i > num_poi) || (i <= 0)) terminator("wrong id (findsrcimg)");
  
  set_distance_lpl_zs(para_poi[i - 1][0]);
  lensmodel(x, y, pout, 1, 0);

  xs = x - pout[0];
  ys = y - pout[1];

  para_poi[i - 1][1] = xs;
  para_poi[i - 1][2] = ys;

  fprintf(stderr, "source moved to: (%e, %e)\n\n", xs, ys);

  if(check_para_poi_all() > 0)
    terminator("invalid input parameter (findsrcimg)");
  
  findimg_i(i);

  return;
}

/*--------------------------------------------------------------
  set def angle and/or allocate matrix for images
*/

void poi_set_table(double zs, int flag_setlcenter, int verb)
{
  int lev, k, kk, i, j, l, nn, flag;
  int stot;
  double dp, tx, ty, ttx, tty, dd, pp, im, pout[NPAR_LMODEL];
  float unit[9][3];
  float *array_tmp;

  if(flag_set_point == 0){
    set_distance_lpl_zs(zs);
    flag_set_point = 1;

    /* assignments of lens centers for 1st multiple lens plane */
    if(flag_setlcenter != 0) set_lens_center_npl0();

    array_poi_nbox[0] = (nx_poi - 1) * (ny_poi - 1);
    nn = nx_poi * ny_poi;
    
    array_tmp = (float*)malloc(sizeof(float) * (3 * nn));
    array_poi_defx[0] = (float*)malloc(sizeof(float) * (4 * array_poi_nbox[0]));
    array_poi_defy[0] = (float*)malloc(sizeof(float) * (4 * array_poi_nbox[0]));
    array_poi_smag[0] = (float*)malloc(sizeof(float) * (4 * array_poi_nbox[0]));
    array_poi_flag[0] = (int*)malloc(sizeof(int) * array_poi_nbox[0]);
    
    if((array_tmp == NULL) || (array_poi_defx[0] == NULL) || (array_poi_defy[0] == NULL) || (array_poi_smag[0] == NULL) || (array_poi_flag[0] == NULL)) 
      terminator("memory allocation failed");
    
    for(k=0;k<nn;k++){
      ktoxy_poi_init(k, &tx, &ty);
      lensmodel(tx, ty, pout, 0, 0);
      
      array_tmp[k] = pout[0];
      array_tmp[k + nn] = pout[1];
      array_tmp[k + 2 * nn] = pout[6];

      /* approximate assignments of lens centers on the fly in case of multiple lens planes */
      if(flag_setlcenter != 0){
	for(j=1;j<nlp;j++){
	  ttx = tx - def_lpl[j][0];
	  tty = ty - def_lpl[j][1];
	  for(i=0;i<num_len;i++){
	    if(lens_lpl_id[i] == j){
	      if(model_lens[i] != 1){
		dd = mag_matrix_dr2(def_lpl[j][2], def_lpl[j][3], def_lpl[j][4], def_lpl[j][5], ttx - para_lens[i][2], tty - para_lens[i][3]);
		pp = dr_lens_center * dp_lev(0);
		if(dd < (pp * pp)){
		  if(num_lcent >= NMAX_LENSCENTER) terminator("too many lens centers");
		  lens_center[num_lcent][0] = tx;
		  lens_center[num_lcent][1] = ty;
		  lens_center_id[num_lcent] = i;
		  num_lcent++;
		}
	      } else {
		for(l=0;l<num_gal;l++){
		  dd = mag_matrix_dr2(def_lpl[j][2], def_lpl[j][3], def_lpl[j][4], def_lpl[j][5], ttx - para_gals[l][0], tty - para_lens[l][1]);
		  
		  pp = dr_lens_center * dp_lev(0);
		  if(dd < (pp * pp)){
		    if(num_lcent >= NMAX_LENSCENTER) terminator("too many lens centers");
		    lens_center[num_lcent][0] = tx;
		    lens_center[num_lcent][1] = ty;
		    lens_center_id[num_lcent] = i;
		    num_lcent++;
		  }
		}
	      }
	    }
	  }
	}
      }
    }   
    
    for(k=0;k<array_poi_nbox[0];k++){
      i = k % (nx_poi - 1);
      j = (k - i) / (nx_poi - 1);
      
      kk = i + j * nx_poi;
      
      for(i=0;i<4;i++){
	j = (i % 2) + nx_poi * (i - (i % 2)) / 2;
	array_poi_defx[0][k + i * array_poi_nbox[0]] = array_tmp[kk + j];
	array_poi_defy[0][k + i * array_poi_nbox[0]] = array_tmp[kk + nn + j];
	array_poi_smag[0][k + i * array_poi_nbox[0]] = array_tmp[kk + 2 * nn + j];
      }
    }
    
    free(array_tmp);
    
    if(verb == 1){
      fprintf(stderr, "   initializing the adaptive mesh...\n");
      if(outformat_exp == 0){
	fprintf(stderr, "   N_box ( lev =  1  pix_size = %8.5f ) = %d\n", pix_poi, array_poi_nbox[0]);
      } else {
	fprintf(stderr, "   N_box ( lev =  1  pix_size = %e ) = %d\n", pix_poi, array_poi_nbox[0]);
      }
    }

    for(k=0;k<array_poi_nbox[0];k++) array_poi_flag[0][k] = 0;
    
    
    for(lev=1;lev<maxlev;lev++){
      dp = dp_lev(lev);
      nn = 0;
      for(k=0;k<array_poi_nbox[lev - 1];k++){
	flag = 0;
	stot = magsigntot(k, lev - 1);
	if((stot != 0) && (stot != 15)){ 
	  flag = 1;
	}  else {
	  for(i=0;i<4;i++){
	    im = fabs(array_poi_smag[lev - 1][k + i * array_poi_nbox[lev - 1]]);
	    if((im > poi_imag_max) || (im < poi_imag_min)) flag = 1;
	  }
	  if(check_mesh_center(k, lev - 1, -1.0, 2.0) == 1) flag = 1;
	} 
	if(flag == 1){ 
	  array_poi_flag[lev - 1][k] = 1;
	  nn++;
	}
      } 
      
      array_poi_nbox[lev] = nn * 4;
      if(verb == 1){
	if(outformat_exp == 0){
	  fprintf(stderr, "   N_box ( lev = %2d  pix_size = %8.5f ) = %d\n", lev + 1, dp, array_poi_nbox[lev]);
	} else {
	  fprintf(stderr, "   N_box ( lev = %2d  pix_size = %e ) = %d\n", lev + 1, dp, array_poi_nbox[lev]);
	}
      }

      if(array_poi_nbox[lev] > NMAX_PIXEL_POINT)
	terminator("pixel number exceeds the limit");
      
      array_poi_defx[lev] = (float*)malloc(sizeof(float) * (4 * array_poi_nbox[lev]));
      array_poi_defy[lev] = (float*)malloc(sizeof(float) * (4 * array_poi_nbox[lev]));
      array_poi_smag[lev] = (float*)malloc(sizeof(float) * (4 * array_poi_nbox[lev]));
      array_poi_flag[lev] = (int*)malloc(sizeof(int) * array_poi_nbox[lev]);
      array_poi_kref[lev] = (int*)malloc(sizeof(int) * array_poi_nbox[lev]);
      
      if((array_poi_defx[lev] == NULL) || (array_poi_defy[lev] == NULL) || (array_poi_smag[lev] == NULL || (array_poi_flag[lev] == NULL) || (array_poi_kref[lev] == NULL)))
	terminator("memory allocation failed");
     
      for(k=0;k<array_poi_nbox[lev];k++) array_poi_flag[lev][k] = 0;
      
      nn = 0;
      for(k=0;k<array_poi_nbox[lev - 1];k++){
	if(array_poi_flag[lev - 1][k] == 1){
	  for(i=0;i<4;i++) array_poi_kref[lev][nn + i] = k;
	  nn = nn + 4;
	}
      }
      
      for(k=0;k<array_poi_nbox[lev];k=k+4){
	
	kk = array_poi_kref[lev][k];
	ktoxy_poi(k, lev, &tx, &ty);
	
	int num[] = {0, 2, 6, 8};
	for(l=0;l<4;l++){
	  unit[num[l]][0] = array_poi_defx[lev - 1][kk + l * array_poi_nbox[lev - 1]];
	  unit[num[l]][1] = array_poi_defy[lev - 1][kk + l * array_poi_nbox[lev - 1]];
	  unit[num[l]][2] = array_poi_smag[lev - 1][kk + l * array_poi_nbox[lev - 1]];
	}
	
	lensmodel(tx + dp, ty, pout, 0, 0);
	unit[1][0] = pout[0];
	unit[1][1] = pout[1];
	unit[1][2] = pout[6];
	lensmodel(tx, ty + dp, pout, 0, 0);
	unit[3][0] = pout[0];
	unit[3][1] = pout[1];
	unit[3][2] = pout[6];
	lensmodel(tx + dp, ty + dp, pout, 0, 0);
	unit[4][0] = pout[0];
	unit[4][1] = pout[1];
	unit[4][2] = pout[6];
	lensmodel(tx + 2.0 * dp, ty + dp, pout, 0, 0);
	unit[5][0] = pout[0];
	unit[5][1] = pout[1];
	unit[5][2] = pout[6];
	lensmodel(tx + dp, ty + 2.0 * dp, pout, 0, 0);
	unit[7][0] = pout[0];
	unit[7][1] = pout[1];
	unit[7][2] = pout[6];
     
	int num2[4][4] = {{0, 1, 3, 4}, {1, 2, 4, 5}, {3, 4, 6, 7}, {4, 5, 7, 8}};
	for(j=0;j<4;j++){
	  for(l=0;l<4;l++){
	    array_poi_defx[lev][k + j + l * array_poi_nbox[lev]] = unit[num2[j][l]][0];
	    array_poi_defy[lev][k + j + l * array_poi_nbox[lev]] = unit[num2[j][l]][1];
	    array_poi_smag[lev][k + j + l * array_poi_nbox[lev]] = unit[num2[j][l]][2];
	  }
	}
      }
    }
    if(verb == 1) fprintf(stderr, "\n");

  }
  
  return;
}

double mag_matrix_dr2(double m11, double m12, double m21, double m22, double dx, double dy)
{
  double det, ddx, ddy;

  det = m11 * m22 - m12 * m21 + imag_ceil;
  
  ddx = (m22 * dx - m12 * dy) / det;
  ddy = ((-1.0) * m12 * dx + m11 * dy) / det;

  return ddx * ddx + ddy * ddy;
}

void poi_unset_table(void)
{
  int lev;
  
  if(flag_set_point == 1){
    flag_set_point = 0;
    for(lev=0;lev<maxlev;lev++){
      free(array_poi_defx[lev]);
      free(array_poi_defy[lev]);
      free(array_poi_smag[lev]);
      free(array_poi_flag[lev]);
      free(array_poi_kref[lev]);
    }
  }
  
  return;
}

/*--------------------------------------------------------------
  check critical curves
*/

int magsigntot(int k, int lev)
{
  int i, l, s, stot;
  float mm;
  
  stot = 0;
  l = 1;
  for(i=0;i<4;i++){ 
    mm = array_poi_smag[lev][k + i * array_poi_nbox[lev]];
    if(mm > 0.0){ s = 0; } else { s = 1; }
    stot = stot + l * s;
    l = l * 2;
  }

  return stot;
}

/*--------------------------------------------------------------
  compute positions of images
*/

void findimg_i(int i)
{
  int ii, ni;
  double zs, xs, ys;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  char fname[INPUT_MAXCHAR];
  FILE* fptr;
  
  if(i > num_poi) terminator("wrong id (findimg)");

  sprintf(fname, "%s_point.dat", file_prefix);
  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (findimg)");

  fclose(fptr); 
  
  for(ii=0;ii<num_poi;ii++){
    zs = para_poi[ii][0];
    xs = para_poi[ii][1];
    ys = para_poi[ii][2];
    
    if((i == 0) || (ii == (i - 1))){
      fprintf(stderr, "######## finding images for point %d\n", ii + 1); 
      fprintf(stderr, " zs = %e  xs = %e  ys = %e\n", zs, xs, ys); 
      findimg(xs, ys, zs, &ni, rr, 1);
    }
  }

  return;
}

void findimg(double xs, double ys, double zs, int *ni, double rr[][NPAR_IMAGE], int verb)
{
  double pout[NPAR_LMODEL];
  int ii, i, k, kk[4], ite, lev;
  int fff[NMAX_POIMG];
  int nimg, nf;
  double xi[NMAX_POIMG], yi[NMAX_POIMG], dpi[NMAX_POIMG];
  double mag[NMAX_POIMG], td[NMAX_POIMG];
  double d1[2], d2[2], d3[2], xx, yy, dp;
  double d12, d23, d31;
  double pxx, pyy, pxy, pyx, ff, gg, mm, dx, dy, xo, yo, dis2, tdmin;
  /* for gap check, 0: bottom, 1: right, 2: top, 3: left */
  double gxs1, gxs2, gxs3, gxs4, gys1, gys2, gys3, gys4;
  double gdx[4] = { 0.4, 2.4, 0.4, -0.4};
  double gdy[4] = {-0.4, 0.4, 2.4,  0.4};
  double gdxs1[4] = {0.0, 2.0, 0.0, 0.0};
  double gdys1[4] = {0.0, 0.0, 2.0, 0.0};
  double gdxs2[4] = {2.0, 2.0, 2.0, 0.0};
  double gdys2[4] = {0.0, 2.0, 2.0, 2.0};
  double gdxs3[4] = {1.0, 2.0, 1.0, 0.0};
  double gdys3[4] = {0.0, 1.0, 2.0, 1.0};
  int gj1[4] = {0, 1, 2, 0};
  int gl1[4] = {0, 1, 2, 0};
  int gj2[4] = {1, 3, 3, 2};
  int gl2[4] = {1, 3, 3, 2};
  int gj3[4] = {0, 1, 2, 0};
  int gl3[4] = {1, 3, 3, 2};

  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  if(verb == 1){
    sprintf(fname, "%s_point.dat", file_prefix);
    fprintf(stderr, " output file name = %s \n\n", fname);
    fptr = fopen(fname, "a");
    
    if(fptr == NULL) terminator("failed at fopen (findimg)");
  }

  if((zs != zs_point) || (flag_set_point == 0)){
    zs_point = zs;
    poi_unset_table();
    poi_set_table(zs, 1, verb);
  }

  nimg = 0;
  
  for(lev=0;lev<maxlev;lev++){
    dp = dp_lev(lev);
    for(k=0;k<array_poi_nbox[lev];k++){
      /* examining gap */
      if((lev > 0) && ((k % 4) == 0)){
	ktoxy_poi(k, lev, &xx, &yy);
	for(i=0;i<4;i++){
	  if(lev > xytolev(xx + gdx[i] * dp, yy + gdy[i] * dp)){
	    gxs1 = xx + gdxs1[i] * dp - array_poi_defx[lev][k + gj1[i] + gl1[i] * array_poi_nbox[lev]];
	    gys1 = yy + gdys1[i] * dp - array_poi_defy[lev][k + gj1[i] + gl1[i] * array_poi_nbox[lev]];
	    gxs2 = xx + gdxs2[i] * dp - array_poi_defx[lev][k + gj2[i] + gl2[i] * array_poi_nbox[lev]];
	    gys2 = yy + gdys2[i] * dp - array_poi_defy[lev][k + gj2[i] + gl2[i] * array_poi_nbox[lev]];
	    gxs3 = xx + gdxs3[i] * dp - array_poi_defx[lev][k + gj3[i] + gl3[i] * array_poi_nbox[lev]];
	    gys3 = yy + gdys3[i] * dp - array_poi_defy[lev][k + gj3[i] + gl3[i] * array_poi_nbox[lev]];
	    gxs4 = xx + dp - array_poi_defx[lev][k + 3 * array_poi_nbox[lev]];
	    gys4 = yy + dp - array_poi_defy[lev][k + 3 * array_poi_nbox[lev]];
	    d1[0] = gxs2 - gxs1;
	    d1[1] = gys2 - gys1;
	    d2[0] = gxs3 - gxs1;
	    d2[1] = gys3 - gys1;
	    d3[0] = gxs4 - gxs1;
	    d3[1] = gys4 - gys1;
	    if((vec_product(d1, d2) * vec_product(d1, d3)) > 0.0){
	      d1[0] = xs - gxs1;
	      d1[1] = ys - gys1;
	      d2[0] = xs - gxs2;
	      d2[1] = ys - gys2;
	      d3[0] = xs - gxs3;
	      d3[1] = ys - gys3;
	      d12 = vec_product(d1, d2);
	      d23 = vec_product(d2, d3);
	      d31 = vec_product(d3, d1);
	      
	      if(((d12 >= 0.0) && (d23 >= 0.0) && (d31 >= 0.0)) || ((d12 <= 0.0) && (d23 <= 0.0) && (d31 <= 0.0))){
		if(nimg < NMAX_POIMG){
		  nimg++;
		  xi[nimg - 1] = xx + gdxs3[i] * dp;
		  yi[nimg - 1] = yy + gdys3[i] * dp;
		  dpi[nimg - 1] = dp;
		}
	      }
	    }
	  }
	}
      }
      /* end examining gap */

      /* standard image checking from here */
      if(array_poi_flag[lev][k] == 0){
	ktoxy_poi(k, lev, &xx, &yy);

	for(i=0;i<4;i++) kk[i] = k + i * array_poi_nbox[lev];

	d1[0] = xs - (xx - array_poi_defx[lev][kk[0]]);
	d1[1] = ys - (yy - array_poi_defy[lev][kk[0]]);
	d2[0] = xs - (xx + dp - array_poi_defx[lev][kk[3]]);
	d2[1] = ys - (yy + dp - array_poi_defy[lev][kk[3]]);
	d3[0] = xs - (xx + dp - array_poi_defx[lev][kk[1]]);
	d3[1] = ys - (yy - array_poi_defy[lev][kk[1]]);
	
	d12 = vec_product(d1, d2);
	d23 = vec_product(d2, d3);
	d31 = vec_product(d3, d1);

	if(((d12 >= 0.0) && (d23 >= 0.0) && (d31 >= 0.0)) || ((d12 <= 0.0) && (d23 <= 0.0) && (d31 <= 0.0))){
	  if(nimg < NMAX_POIMG){
	    nimg++;
	    xi[nimg - 1] = xx + 0.667 * dp;
	    yi[nimg - 1] = yy + 0.333 * dp;
	    dpi[nimg - 1] = dp;
	  }
	}

	d3[0] = xs - (xx - array_poi_defx[lev][kk[2]]);
	d3[1] = ys - (yy + dp - array_poi_defy[lev][kk[2]]);
	
	d23 = vec_product(d2, d3);
	d31 = vec_product(d3, d1);

	if(((d12 > 0.0) && (d23 > 0.0) && (d31 > 0.0)) || ((d12<0.0) && (d23<0.0) && (d31<0.0))){
	  if(nimg < NMAX_POIMG){
	    nimg++;
	    xi[nimg - 1] = xx + 0.333 * dp;
	    yi[nimg - 1] = yy + 0.667 * dp;
	    dpi[nimg - 1] = dp;
	  }
	}
      }
    }
  }
  
  /* improve solutions */
  for(ii=0;ii<nimg;ii++){
    ite = 0;
    xo = xi[ii];
    yo = yi[ii];
    do{
      ite++;
      lensmodel(xi[ii], yi[ii], pout, 0, 0);
      pxx = pout[3] + pout[4];
      pyy = pout[3] - pout[4];
      pxy = pout[5] + pout[7];
      pyx = pout[5] - pout[7];
      ff = xs - xi[ii] + pout[0];
      gg = ys - yi[ii] + pout[1];
      mm = (1.0 - pxx) * (1.0 - pyy) - pxy * pyx;
      dx = ((1.0 - pyy) * ff + pxy * gg) / mm;
      dy = ((1.0 - pxx) * gg + pyx * ff) / mm;
      xi[ii] = xi[ii] + dx;
      yi[ii] = yi[ii] + dy;
    }while(((fabs(ff) > max_poi_tol) || (fabs(gg) > max_poi_tol)) && (ite <= nmax_poi_ite));
    if(((xo - xi[ii]) * (xo - xi[ii]) + (yo - yi[ii]) * (yo - yi[ii])) > (2.0 * dpi[ii] * dpi[ii])){
      xi[ii] = xo;
      yi[ii] = yo;
      fff[ii] = 1;
    } else {
      fff[ii] = 0;
    }
    lensmodel(xi[ii], yi[ii], pout, -1, 0);
    mag[ii] = 1.0 / (pout[6] + imag_ceil);
    td[ii] = pout[2];
  } 
  
  /* remove same images */
  nf = nimg;
  for(ii=0;ii<nimg;ii++){
    for(i=ii+1;i<nimg;i++){
      mm = fabs(mag[ii] * mag[i]);
      dis2 = ((xi[ii] - xi[i]) * (xi[ii] - xi[i]) + (yi[ii] - yi[i]) * (yi[ii] - yi[i])) / mm;
      if((dis2 <= (10.0 * max_poi_tol * max_poi_tol)) && (fff[ii] >= 0)){
	fff[ii] = -1;
	nf--;
      }
    }
  }
  
  /* remove `iteration failed' images */
  for(ii=0;ii<nimg;ii++){
    if(fff[ii] == 1){
	fff[ii] = -1;
	nf--;
    }
  }

  /* calculate time delay zero-point */
  tdmin = TDMIN_SET;
  for(ii=0;ii<nimg;ii++){
    if((fff[ii] >= 0) && (td[ii] < tdmin)) tdmin = td[ii];
  }

  for(ii=0;ii<nimg;ii++){
    td[ii] = td[ii] - tdmin;
  } 
 
  if(verb == 1){
    printf("n_img = %d\n", nf);
    if(outformat_exp == 0){
      fprintf(fptr, "%d %8.4f %9.4f %9.4f\n", nf, zs, xs, ys);
    } else{
      fprintf(fptr, "%d %13e %13e %13e\n", nf, zs, xs, ys);
    }
    for(ii=0;ii<nimg;ii++){
      if(fff[ii] >= 0){
	if(outformat_exp == 0){
	  printf("x = %9.4f   y = %9.4f   mag = %9.4f [%8.3f]   td[day] = %9.3f", xi[ii], yi[ii], mag[ii], (-2.5) * log10(fabs(mag[ii]) + OFFSET_LOG), td[ii]);
	  fprintf(fptr, "%9.4f %9.4f %9.4f %9.3f\n", xi[ii], yi[ii], mag[ii], td[ii]);
	} else {
	  printf("x = %13e   y = %13e   mag = %13e [%13e]   td[day] = %13e", xi[ii], yi[ii], mag[ii], (-2.5) * log10(fabs(mag[ii]) + OFFSET_LOG), td[ii]);
	  fprintf(fptr, "%13e %13e %13e %13e\n", xi[ii], yi[ii], mag[ii], td[ii]);
	}
	if(fff[ii] == 1){
	  printf("  [iteration failed]\n");
	} else {
	  printf("\n");
	}
      }
    }
    fprintf(stderr, "\n");
    fclose(fptr);
  }

  (*ni) = nf;
  i = 0;
  for(ii=0;ii<nimg;ii++){
    if(fff[ii] >= 0){
      rr[i][0] = xi[ii];
      rr[i][1] = yi[ii];
      rr[i][2] = mag[ii];
      rr[i][3] = td[ii];
      i++;
    }
  }

  return;
}

double vec_product(double d1[2], double d2[2])
{
  return d1[0] * d2[1] - d1[1] * d2[0];
}

/*--------------------------------------------------------------
  get full list of lens centers by solving lens equations
*/

void full_lens_center(double zs)
{
  int i, j, k, l, ni;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  int tmp_nlp, tmp_num_lcent;
  double tmp_zl_lpl[NMAX_LPL];
  double tmp_lens_lpl_id[NMAX_LEN];
  double tmp_lens_center[NMAX_LENSCENTER][2];
  int tmp_lens_center_id[NMAX_LENSCENTER];

  set_distance_lpl_zs(zs);
  set_lens_center_npl0();
  
  tmp_num_lcent = num_lcent;
  for(i=0;i<tmp_num_lcent;i++){
    tmp_lens_center[i][0] = lens_center[i][0];
    tmp_lens_center[i][1] = lens_center[i][1];
    tmp_lens_center_id[i] = lens_center_id[i];
  }
  tmp_nlp = nlp;
  for(j=0;j<nlp;j++) tmp_zl_lpl[j] = zl_lpl[j];
  for(i=0;i<num_len;i++) tmp_lens_lpl_id[i] = lens_lpl_id[i];
  
  poi_unset_table();

  for(j=1;j<tmp_nlp;j++){
    for(i=0;i<num_len;i++){
      if(tmp_lens_lpl_id[i] == j){
	if(model_lens[i] != 1){
	  findimg(para_lens[i][2], para_lens[i][3], tmp_zl_lpl[j], &ni, rr, 0);
	  for(k=0;k<ni;k++){
	    if(tmp_num_lcent >= NMAX_LENSCENTER) terminator("too many lens centers");
	    tmp_lens_center[tmp_num_lcent][0] = rr[k][0];
	    tmp_lens_center[tmp_num_lcent][1] = rr[k][1];
	    tmp_lens_center_id[tmp_num_lcent] = i;
	    tmp_num_lcent++;
	  }
	} else {
	  for(l=0;l<num_gal;l++){
	    findimg(para_gals[l][0], para_gals[l][1], tmp_zl_lpl[j], &ni, rr, 0);
	    for(k=0;k<ni;k++){
	      if(tmp_num_lcent >= NMAX_LENSCENTER) terminator("too many lens centers");
	      tmp_lens_center[tmp_num_lcent][0] = rr[k][0];
	      tmp_lens_center[tmp_num_lcent][1] = rr[k][1];
	      tmp_lens_center_id[tmp_num_lcent] = i;
	      tmp_num_lcent++;
	    }
	  }
	}
      }
    }
  }

  num_lcent = tmp_num_lcent;
  for(i=0;i<num_lcent;i++){
    lens_center[i][0] = tmp_lens_center[i][0];
    lens_center[i][1] = tmp_lens_center[i][1];
    lens_center_id[i] = tmp_lens_center_id[i];
  }
  
  return;
}

void lenscenter(double zs)
{
  int i;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  fprintf(stderr, "######## computing centers of all lenses\n");
  fprintf(stderr, " zs = %e\n", zs);
  
  sprintf(fname, "%s_lenscenter.dat", file_prefix);
  fprintf(stderr, " output file name = %s\n\n", fname);
  
  fptr = fopen(fname, "w");
  if(fptr == NULL) terminator("failed at fopen (lenscenter)");

  full_lens_center(zs);

  for(i=0;i<num_lcent;i++){
    fprintf(stderr, "%d %f %13e %13e\n", lens_center_id[i] + 1, zl_lpl[lens_lpl_id[lens_center_id[i]]], lens_center[i][0], lens_center[i][1]);
    fprintf(fptr, "%d %f %13e %13e\n", lens_center_id[i] + 1, zl_lpl[lens_lpl_id[lens_center_id[i]]], lens_center[i][0], lens_center[i][1]);
  }

  fprintf(stderr, "\n");
  fclose(fptr);

  return;
} 

/*--------------------------------------------------------------
  write critcurves and caustics
*/

void writecrit(double zs, double causize[NPAR_CAUSIZE ], int flag_newcent, int verb)
{
  int lev, i, ii, k, stot, kk[4];
  int l[4] = {0, 0, 0, 0};
  double dp;
  /* 0: -  1: |  2: \ (low) 3: / (low) 4: / (up) 5: \ (up) */
  int flag[6];
  double pout[NPAR_LMODEL], xc[4], yc[4], xx, yy;
  double x1, y1, x2, y2, sx1, sx2, sy1, sy2;
  double t, dt, cx, cy, s1, s2, s3;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;

  if((flag_newcent != 0) && (num_lpl > 1)){
    full_lens_center(zs);
    zs_point = zs;
    poi_unset_table();
    poi_set_table(zs, 0, 1);
  } else if((zs != zs_point) || (flag_set_point == 0)){
    zs_point = zs;
    poi_unset_table();
    poi_set_table(zs, 1, 1);
  }

  causize[0] = causize[2] = CAUSIZE_SET;
  causize[1] = causize[3] = (-1.0) * CAUSIZE_SET;
  
  if(verb == 1){
    sprintf(fname, "%s_crit.dat", file_prefix);
  
    fprintf(stderr, "######## writing critical curve\n");
    fprintf(stderr, " zs = %e \n", zs);
    fprintf(stderr, " output file name = %s \n\n", fname);
  
    fptr = fopen(fname, "w");
    if(fptr == NULL) terminator("failed at fopen (writecrit)");
  }

  for(lev=0;lev<maxlev;lev++){
    dp = dp_lev(lev);
    for(k=0;k<array_poi_nbox[lev];k++){
      if((array_poi_flag[lev][k] == 0) && (check_mesh_center(k, lev, -0.001, 1.001) == 0)){
	
	for(ii=0;ii<6;ii++){ flag[ii] = 0;  } 
	ktoxy_poi(k, lev, &xx, &yy);
	
	for(i=0;i<4;i++) kk[i] = k + i * array_poi_nbox[lev];

	stot = magsigntot(k, lev);
	
	switch(stot){
	case 0: 
	  break;
	  
	case 1:
	  flag[2]++;
	  break;
	  
	case 2:
	  flag[3]++;
	  break;
	  
	case 3:
	  flag[0]++;
	  break;
	  
	case 4:
	  flag[4]++;
	  break;
	  
	case 5:
	  flag[1]++;
	  break;
	  
	case 6:
	  lensmodel(xx + 0.5 * dp, yy + 0.5 * dp, pout, 0, 0);
	  if(pout[6] > 0.0){ 
	    flag[3]++;
	    flag[4]++; 
	  } else { 
	    flag[2]++;
	    flag[5]++;
	  }
	  break;
	  
	case 7:
	  flag[5]++;
	  break;
	  
	case 8:
	  flag[5]++;
	  break;
	  
	case 9:
	  lensmodel(xx + 0.5 * dp, yy + 0.5 * dp, pout, 0, 0);
	  if(pout[6] > 0.0){ 
	    flag[2]++;
	    flag[5]++;
	  } else { 
	    flag[3]++;
	    flag[4]++; 
	  }
	  break;

	case 10:
	  flag[1]++;
	  break;
	  
	case 11:
	  flag[4]++;
	  break;
	  
	case 12:
	  flag[0]++;
	  break;
	  
	case 13:
	  flag[3]++;
	  break;
	  
	case 14:
	  flag[2]++;
	  break;

	case 15:
	  break;
	}
	
	for(ii=0;ii<6;ii++){
	  if(flag[ii] == 1){
	    xc[0] = xx;    yc[0] = yy;
	    xc[1] = xx + dp; yc[1] = yy;
	    xc[2] = xx;    yc[2] = yy + dp;
	    xc[3] = xx + dp; yc[3] = yy + dp;

	    switch(ii){
	    case 0:   
	      l[0] = 0; l[1] = 2; l[2] = 1; l[3] = 3;
	      break;
	      
	    case 1:
	      l[0] = 0; l[1] = 1; l[2] = 2; l[3] = 3;
	      break;
	      
	    case 2:
	      l[0] = 0; l[1] = 1; l[2] = 0; l[3] = 2;
	      break;
	      
	    case 3:
	      l[0] = 0; l[1] = 1; l[2] = 1; l[3] = 3;
	      break;
	      
	    case 4:
	      l[0] = 0; l[1] = 2; l[2] = 2; l[3] = 3;
	      break;
	      
	    case 5:
	      l[0] = 1; l[1] = 3; l[2] = 2; l[3] = 3;
	      break;
	    }
	    
	    x1 = 0.5 * (xc[l[0]] + xc[l[1]]);
	    y1 = 0.5 * (yc[l[0]] + yc[l[1]]);
	    x2 = 0.5 * (xc[l[2]] + xc[l[3]]);
	    y2 = 0.5 * (yc[l[2]] + yc[l[3]]);
	    
	    sx1 = x1 - 0.5 * (array_poi_defx[lev][kk[l[0]]] + array_poi_defx[lev][kk[l[1]]]);
	    sy1 = y1 - 0.5 * (array_poi_defy[lev][kk[l[0]]] + array_poi_defy[lev][kk[l[1]]]);
	    sx2 = x2 - 0.5 * (array_poi_defx[lev][kk[l[2]]] + array_poi_defx[lev][kk[l[3]]]);
	    sy2 = y2 - 0.5 * (array_poi_defy[lev][kk[l[2]]] + array_poi_defy[lev][kk[l[3]]]);
	    
	    if(sx1 < causize[0]){ causize[0] = sx1; }
	    if(sx1 > causize[1]){ causize[1] = sx1; }
	    if(sy1 < causize[2]){ causize[2] = sy1; }
	    if(sy1 > causize[3]){ causize[3] = sy1; }

	    if(verb == 1){
	      fprintf(fptr, "%e %e %e %e %e %e %e %e\n", x1, y1, sx1, sy1, x2, y2, sx2, sy2);
	    }
	  }
	}
      } 
    }
  }

  /* check centers of galaxies
     for the case inner crit curves are degenerate at the center
  */
  dp = dp_lev(maxlev - 1);
  i = 0;
  dt = 2.0 / ((double)center_ang_step);
  do{
    cx = lens_center[i][0];
    cy = lens_center[i][1];
    i++;
    for(t=0.0;t<(2.0-(0.5*dt));t=t+dt){
      x1 = dp * cos(t * M_PI) + cx;
      y1 = dp * sin(t * M_PI) + cy;
      x2 = dp * cos((t + dt) * M_PI) + cx;
      y2 = dp * sin((t + dt) * M_PI) + cy;
      lensmodel(x1, y1, pout, 0, 0);
      s1 = pout[6];
      lensmodel(x2, y2, pout, 0, 0);
      s2 = pout[6];
      lensmodel(cx, cy, pout, 0, 0);
      s3 = pout[6];
      if(((s1 > 0.0) && (s2 > 0.0) && (s3 <= 0.0)) || ((s1 < 0.0) && (s2 < 0.0) && (s3 >= 0.0))){
	x1 = 0.5 * dp * cos(t * M_PI) + cx;
	y1 = 0.5 * dp * sin(t * M_PI) + cy;
	x2 = 0.5 * dp * cos((t + dt) * M_PI) + cx;
	y2 = 0.5 * dp * sin((t + dt) * M_PI) + cy;
	lensmodel(x1, y1, pout, 1, 0);
	sx1 = x1 - pout[0];
	sy1 = y1 - pout[1];
	lensmodel(x2, y2, pout, 1, 0);
	sx2 = x2 - pout[0];
	sy2 = y2 - pout[1];

	if(sx1 < causize[0]){ causize[0] = sx1; }
	if(sx1 > causize[1]){ causize[1] = sx1; }
	if(sy1 < causize[2]){ causize[2] = sy1; }
	if(sy1 > causize[3]){ causize[3] = sy1; }

	if(verb == 1){
	  fprintf(fptr, "%e %e %e %e %e %e %e %e\n", cx, cy, sx1, sy1, cx, cy, sx2, sy2);
	}
      }
    }
  }while(i < num_lcent);

  if(verb == 1){
    fclose(fptr);
  }

  return;
}

/*--------------------------------------------------------------
  write the strcture of the adaptove mesh
*/

void writemesh(double zs)
{
  int lev, k, i, kk[4];
  double x, y, xx[4], yy[4], sx[4], sy[4], dp;
  char fname[INPUT_MAXCHAR];
  FILE* fptr;
  int ii[4] = {0, 1, 3, 2};
  int jj[4] = {1, 3, 2, 0};

  if((zs != zs_point) || (flag_set_point == 0)){
    zs_point = zs;
    poi_unset_table();
    poi_set_table(zs, 1, 1);
  }
  
  sprintf(fname, "%s_mesh.dat", file_prefix);
  
  fprintf(stderr, "######## writing mesh for point sources\n");
  fprintf(stderr, " zs = %e \n", zs);
  fprintf(stderr, " output file name = %s \n\n", fname);
    
  fptr = fopen(fname, "w");

  if(fptr == NULL) terminator("failed at fopen (writemesh)");

  for(lev=0;lev<maxlev;lev++){
    dp = dp_lev(lev);
    for(k=0;k<array_poi_nbox[lev];k++){
      if(array_poi_flag[lev][k] == 0){
	ktoxy_poi(k, lev, &x, &y);
	xx[0] = x;    yy[0] = y;    kk[0] = k;
	xx[1] = x + dp; yy[1] = y;    kk[1] = k + array_poi_nbox[lev];
	xx[2] = x;    yy[2] = y + dp; kk[2] = k + 2 * array_poi_nbox[lev];
	xx[3] = x + dp; yy[3] = y + dp; kk[3] = k + 3 * array_poi_nbox[lev];
	for(i = 0;i<4;i++){
	  sx[i] = xx[i] - array_poi_defx[lev][kk[i]];
	  sy[i] = yy[i] - array_poi_defy[lev][kk[i]];
	}
	for(i=0;i<4;i++){
	  fprintf(fptr, "%e %e %e %e %e %e %e %e\n", xx[ii[i]], yy[ii[i]], sx[ii[i]], sy[ii[i]], xx[jj[i]], yy[jj[i]], sx[jj[i]], sy[jj[i]]);
	}
      }
    }
  }

  fclose(fptr);
  
  return;
}

/*--------------------------------------------------------------
  convert int to x and y
*/

void ktoxy_poi(int k, int lev, double *x, double *y)
{
  int i, j, l, km;
  double x0, y0, dx, dy, dp;
  
  if(lev == 0){
    i = k % (nx_poi - 1);
    j = (k - i) / (nx_poi - 1);
    
    *x = xmin + pix_poi * (i + 0.5);
    *y = ymin + pix_poi * (j + 0.5);
  } else {
    dp = dp_lev(lev);

    l = k % 4;
    
    dx = 0.0;
    dy = 0.0;
    if((l == 1) || (l == 3)){ dx = dp; }
    if((l == 2) || (l == 3)){ dy = dp; }

    km = array_poi_kref[lev][k];
    ktoxy_poi(km, lev - 1, &x0, &y0);

    *x = x0 + dx;
    *y = y0 + dy;
  }

  return;
}

void ktoxy_poi_init(int k, double *x, double *y)
{
  int i, j;
  
  i = k % nx_poi;
  j = (k - i) / nx_poi;
  
  *x = xmin + pix_poi * (i + 0.5);
  *y = ymin + pix_poi * (j + 0.5);
  
  return;
}

/*--------------------------------------------------------------
  check subgrid level at each point (x, y)
*/

int xytolev(double x, double y)
{
  int r, i, j, k0, k1;
  double p0, dx, dy;

  i = (int)(((x - xmin) / pix_poi) - 0.5);
  j = (int)(((y - ymin) / pix_poi) - 0.5);
  k0 = i + j * (nx_poi - 1);
  k1 = 0;

  if((k0 < 0) || (k0 >= array_poi_nbox[0])){
    r = -1;
  } else if(array_poi_flag[0][k0] == 0){
    r = 0;
  } else {
    dx = x - (xmin + pix_poi * (i + 0.5));
    dy = y - (ymin + pix_poi * (j + 0.5));
    p0 = pix_poi;
    r = 0;
    do{
      r++;
      p0 = p0 * 0.5;
      i = 0;
      if(dx >= p0){ i = i + 1; dx = dx - p0; }
      if(dy >= p0){ i = i + 2; dy = dy - p0; }
      for(j=0;j<array_poi_nbox[r];j=j+4){
	if(array_poi_kref[r][j] == k0) k1 = j + i;
      }
      k0 = k1;
    }while(array_poi_flag[r][k1] != 0);
  }

  return r;
}

/*--------------------------------------------------------------
  check if a grid is in the centers of lens objects or not
*/

int check_mesh_center(int k, int lev, double r1, double r2)
{
  int i, f;
  double tx, ty, x, y, dp;

  ktoxy_poi(k, lev, &tx, &ty);
  dp = dp_lev(lev);

  f = 0;

  i = 0;
  do{
    x = lens_center[i][0];
    y = lens_center[i][1];
    i++;
    if((x >= (tx + r1 * dp)) && (x <= (tx + r2 * dp)) && (y >= (ty + r1 * dp)) && (y <= (ty + r2 * dp))) f = 1;
  }while(i < num_lcent);

  return f;
}

void set_lens_center_npl0(void)
{
  int i, l;

  num_lcent = 0;
  for(i=0;i<num_len;i++){
    if(lens_lpl_id[i] == 0){
      if(model_lens[i] != 1){
	if(num_lcent >= NMAX_LENSCENTER) terminator("too many lens centers");
	lens_center[num_lcent][0] = para_lens[i][2];
	lens_center[num_lcent][1] = para_lens[i][3];
	lens_center_id[num_lcent] = i;
	num_lcent++;
      } else {
	for(l=0;l<num_gal;l++){
	  if(num_lcent >= NMAX_LENSCENTER) terminator("too many lens centers");
	  lens_center[num_lcent][0] = para_gals[l][0];
	  lens_center[num_lcent][1] = para_gals[l][1];
	  lens_center_id[num_lcent] = i;
	  num_lcent++;
	}
      }
    }
  }

  return;  
}

/*--------------------------------------------------------------
  mesh size at each level
*/

double dp_lev(int lev)
{
  int i;
  double dp;
  
  dp = pix_poi;
  for(i=0;i<lev;i++) dp = dp * 0.5;

  return dp;
}

/*--------------------------------------------------------------
  set image table for point sources
*/

void poimg_set_table(int flag_source)
{
  int i, ni, nn, ntot, j, k;
  double zs, xs, ys, x, y, f;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  
  if(num_poi < 1) terminator("setting table failed");
  if(flag_seeing == 0) terminator("no PSF set");

  nn = nx_ext * ny_ext;
  ntot = nn * num_poi;

  if(flag_set_poimg == 0){
    flag_set_poimg = 1;
    array_poimg = (float*)malloc(sizeof(float) * ntot);
    if(array_poimg == NULL) terminator("memory allocation failed");
 }

  for(i=0;i<num_poi;i++){
    zs = para_poi[i][0];
    xs = para_poi[i][1];
    ys = para_poi[i][2];

    if(flag_source == 0){
      findimg(xs, ys, zs, &ni, rr, 0);
    } else {
      ni = 1;
      rr[0][0] = xs;
      rr[0][1] = ys;
      rr[0][2] = 1.0;
    }

    for(k=0;k<nn;k++){
      ktoxy_ext(k, &x, &y);
      f = 0.0;
      for(j=0;j<ni;j++){
	f = f + flux_poi[i] * fabs(rr[j][2]) * source_psf_pix(x, y, rr[j][0], rr[j][1], pix_ext);
      } 
      array_poimg[k + i * nn] = f;
    }
  }

  return;
}

void poimg_unset_table(void)
{
  if(flag_set_poimg == 1){
    free(array_poimg);
    flag_set_poimg = 0;
  }

  return;
}

void read_poimg_flux(char *infile)
{
  int i, j, f, nn, np;
  char buffer[INPUT_MAXCHAR];
  FILE* fptr;
  
  fprintf(stderr, "######## reading flux file for writeimageall\n");
  fprintf(stderr, " input file name = %s \n\n", infile);

  fptr = fopen(infile, "r");
  
  if(fptr == NULL) terminator("failed at fopen (point_flux)");

  f = 0;
  j = 0;
  while((fgets(buffer, INPUT_MAXCHAR, fptr)) && (f == 0)){
    if(buffer[0] != '#'){
      nn = sscanf(buffer, "%d", &np);
      if(np != num_poi) terminator("invalid format (point_flux)");
      if (nn != EOF){
	f = 1;
	for(i=1;i<=np;i++){
	  if(fgets(buffer, INPUT_MAXCHAR, fptr)){
	    nn = sscanf(buffer, "%lf", &flux_poi[i - 1]);
	    if(nn == 1) j++;
	  }
	}
	if(np != j) terminator("invalid format (point_flux)");
      }
    }
  }

  fprintf(stderr, "read %d fluxes\n\n", j);

  fclose(fptr);

  return;
}

