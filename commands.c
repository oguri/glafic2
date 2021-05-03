#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

/*--------------------------------------------------------------
  execute command
*/

int do_command(char *buffer)
{
  int i, j, k, l, n, r, f;
  double zs, zl, x, y, x1, y1, x2, y2, r1, r2, rr, e1, e2;
  double csize[NPAR_CAUSIZE];
  char keyword[INPUT_MAXCHAR];
  char cdummy[INPUT_MAXCHAR];
  char fname[INPUT_MAXCHAR];
  char fname2[INPUT_MAXCHAR];

  r = 0;

  if(sscanf(buffer, "%s", keyword) != EOF){
    if(keyword[0] != '#'){
      
      /*------------------------------------
        debug
      */
      if(strcmp(keyword, "debug") == 0){
	deb();
      }
      
      /*------------------------------------
        properties of lens
      */
      if(strcmp(keyword, "calcimage") == 0){
	zs = DEF_COMMAND_ZS;
	x = DEF_COMMAND_XY; y = DEF_COMMAND_XY;
	sscanf(buffer, "%s %lf %lf %lf", cdummy, &zs, &x, &y);
	calcimg(x, y, zs);
      }
      
      if(strcmp(keyword, "calcein") == 0){
	zs = DEF_COMMAND_ZS;
	sscanf(buffer, "%s %lf", cdummy, &zs);
	calcein(zs);
      }
      
      if(strcmp(keyword, "calcein2") == 0){
	zs = DEF_COMMAND_ZS;
	x = DEF_COMMAND_XY; y = DEF_COMMAND_XY; i = DEF_COMMAND_I;
	sscanf(buffer, "%s %lf %lf %lf %d", cdummy, &zs, &x, &y, &i);
	calcein2(zs, x, y, i);
      }

      if(strcmp(keyword, "calcmr") == 0){
	sscanf(buffer, "%s", cdummy);
	calcmr();
      }
      
      if(strcmp(keyword, "writelens") == 0){
	zs = DEF_COMMAND_ZS;
	sscanf(buffer, "%s %lf", cdummy, &zs);
	writelens(zs);
      }

      if(strcmp(keyword, "kapparad") == 0){
	zs = DEF_COMMAND_ZS;
	x = DEF_COMMAND_XY; y = DEF_COMMAND_XY;
	r1 = DEF_COMMAND_R1; r2 = DEF_COMMAND_R2; n = DEF_COMMAND_N;
	i = DEF_COMMAND_I;
	sscanf(buffer, "%s %lf %lf %lf %lf %lf %d %d", cdummy, &zs, &x, &y, &r1, &r2, &n, &i);
	kappa_rad_out(zs, x, y, r1, r2, n, i);
      }
      
      if(strcmp(keyword, "kappacum") == 0){
	zs = DEF_COMMAND_ZS;
	x = DEF_COMMAND_XY; y = DEF_COMMAND_XY;
	r1 = DEF_COMMAND_R1; r2 = DEF_COMMAND_R2; n = DEF_COMMAND_N;
	i = DEF_COMMAND_I;
	sscanf(buffer, "%s %lf %lf %lf %lf %lf %d %d", cdummy, &zs, &x, &y, &r1, &r2, &n, &i);
	kappa_cum_out(zs, x, y, r1, r2, n, i);
      }
       
      if(strcmp(keyword, "lenscenter") == 0){
	zs = DEF_COMMAND_ZS;
	sscanf(buffer, "%s %lf", cdummy, &zs);
	lenscenter(zs);
      }
       
      /*------------------------------------
	extended sources
      */
      if(strcmp(keyword, "writeimage") == 0){
	x = DEF_COMMAND_SKY; y = DEF_COMMAND_NOISE;
	sscanf(buffer, "%s %lf %lf", cdummy, &x, &y);
	writeimage(x, y, 0);
      }
      
      if(strcmp(keyword, "writeimage_ori") == 0){
	x = DEF_COMMAND_SKY; y = DEF_COMMAND_NOISE;
	sscanf(buffer, "%s %lf %lf", cdummy, &x, &y);
	writeimage(x, y, 1);
      }
       
      if(strcmp(keyword, "calcextend") == 0){
	i = DEF_COMMAND_I;
	x = DEF_COMMAND_FTH;
	y = DEF_COMMAND_RLIM;
	sscanf(buffer, "%s %d %lf %lf", cdummy, &i, &x, &y);
	ext_est_image(i, 0, x, y);
      }
       
      if(strcmp(keyword, "calcextend_ori") == 0){
	i = DEF_COMMAND_I;
	x = DEF_COMMAND_FTH;
	y = DEF_COMMAND_RLIM;
	sscanf(buffer, "%s %d %lf %lf", cdummy, &i, &x, &y);
	ext_est_image(i, 1, x, y);
      }
       
      if(strcmp(keyword, "writetd_extend") == 0){
	writetd_ext();
      }
      
      if(strcmp(keyword, "writepsf") == 0){
	sscanf(buffer, "%s", cdummy);
	writepsf();
      }
      
      if(strcmp(keyword, "readpsf") == 0){
	sscanf(buffer, "%s %s", cdummy, fname);
	read_psffits(fname, 1);
       }
      
      /*------------------------------------
	point sources
      */
      if(strcmp(keyword, "findimg") == 0){
	i = DEF_COMMAND_I;
	sscanf(buffer, "%s %d", cdummy, &i);
	findimg_i(i);
      }
      
      if(strcmp(keyword, "findsrcimg") == 0){
	i = DEF_COMMAND_I; x1 = DEF_COMMAND_XY; y1 = DEF_COMMAND_XY;
	sscanf(buffer, "%s %d %lf %lf", cdummy, &i, &x1, &y1);
	findsrcimg(i, x1, y1);
      }
      
      if(strcmp(keyword, "writecrit") == 0){
	zs = DEF_COMMAND_ZS;
	sscanf(buffer, "%s %lf", cdummy, &zs);
	writecrit(zs, csize, 1, 1);
      }
      
      if(strcmp(keyword, "writemesh") == 0){
	zs = DEF_COMMAND_ZS;
	sscanf(buffer, "%s %lf", cdummy, &zs);
	writemesh(zs);
      }
      
      if(strcmp(keyword, "writetd_point") == 0){
	writetd_poi();
      }
      
      if(strcmp(keyword, "writelens_splane") == 0){
	zs = DEF_COMMAND_ZS;
	x1 = (-1.0) * DEF_COMMAND_R2; x2 = DEF_COMMAND_R2;
	y1 = (-1.0) * DEF_COMMAND_R2; y2 = DEF_COMMAND_R2;
	rr = DEF_COMMAND_R1;
	sscanf(buffer, "%s %lf %lf %lf %lf %lf %lf", cdummy, &zs, &x1, &x2, &y1, &y2, &rr);
	writelens_splane(zs, x1, x2, y1, y2, rr);
      }

      /*------------------------------------
	write both point and extended sources
      */
      if(strcmp(keyword, "point_flux") == 0){
	sscanf(buffer, "%s %s", cdummy, fname);
	read_poimg_flux(fname);
      }
      
      if(strcmp(keyword, "writeimageall") == 0){
	x = DEF_COMMAND_SKY; y = DEF_COMMAND_NOISE;
	sscanf(buffer, "%s %lf %lf", cdummy, &x, &y);
	writeimageall(x, y, 0);
      }
       
      if(strcmp(keyword, "writeimageall_ori") == 0){
	x = DEF_COMMAND_SKY; y = DEF_COMMAND_NOISE;
	sscanf(buffer, "%s %lf %lf", cdummy, &x, &y);
	writeimageall(x, y, 1);
      }
       
      /*------------------------------------
	mock
      */
      if(strcmp(keyword, "mock1") == 0){
	zs = DEF_COMMAND_ZS;
	i = DEF_COMMAND_MOCK;
	x1 = (-1.0) * DEF_COMMAND_MX; x2 = DEF_COMMAND_MX;
	y1 = (-1.0) * DEF_COMMAND_MY; y2 = DEF_COMMAND_MY;
	sscanf(buffer, "%s %d %lf %lf %lf %lf %lf", cdummy, &i, &zs, &x1, &x2, &y1, &y2);
	mock1(i, zs, x1, x2, y1, y2);
      }
      
      if(strcmp(keyword, "mock2") == 0){
	zs = DEF_COMMAND_ZS;
	i = DEF_COMMAND_MOCK;
	rr = DEF_COMMAND_MR; x1 = DEF_COMMAND_XY; y1 = DEF_COMMAND_XY;
	sscanf(buffer, "%s %d %lf %lf %lf %lf", cdummy, &i, &zs, &rr, &x1, &y1);
	mock2(i, zs, rr, x1, y1);
      }
      
      if(strcmp(keyword, "mock3") == 0){
	zs = DEF_COMMAND_ZS;
	i = DEF_COMMAND_MOCK;
	rr = DEF_COMMAND_FREC;
	sscanf(buffer, "%s %d %lf %lf", cdummy, &i, &zs, &rr);
	mock3(i, zs, rr);
      }

      if(strcmp(keyword, "mockline") == 0){
	zs = DEF_COMMAND_ZS;
	i = DEF_COMMAND_MOCK;
	j = DEF_COMMAND_FLAG_OUT;
	x1 = (-1.0) * DEF_COMMAND_MX; x2 = DEF_COMMAND_MX;
	y1 = (-1.0) * DEF_COMMAND_MY; y2 = DEF_COMMAND_MY;
	sscanf(buffer, "%s %d %lf %lf %lf %lf %lf %d", cdummy, &i, &zs, &x1, &x2, &y1, &y2, &j);
	mockline(i, zs, x1, x2, y1, y2, j);
      }

      if(strcmp(keyword, "mockext1") == 0){
	i = DEF_COMMAND_MOCK; j = DEF_COMMAND_EI;
	x1 = (-1.0) * DEF_COMMAND_MX; x2 = DEF_COMMAND_MX;
	y1 = (-1.0) * DEF_COMMAND_MY; y2 = DEF_COMMAND_MY;
	e1 = DEF_COMMAND_ELLIP; e2 = DEF_COMMAND_ELLIP;
	x = DEF_COMMAND_FTH; y = DEF_COMMAND_RLIM; 
	sscanf(buffer, "%s %d %d %lf %lf %lf %lf %lf %lf %lf %lf", cdummy, &i, &j, &x1, &x2, &y1, &y2, &e1, &e2, &x, &y);
	mockext1(i, j, x1, x2, y1, y2, e1, e2, x, y);
      }

      if(strcmp(keyword, "mockext2") == 0){
	i = DEF_COMMAND_MOCK; j = DEF_COMMAND_EI;
	rr = DEF_COMMAND_MR; x1 = DEF_COMMAND_XY; y1 = DEF_COMMAND_XY;
	e1 = DEF_COMMAND_ELLIP; e2 = DEF_COMMAND_ELLIP;
	x = DEF_COMMAND_FTH; y = DEF_COMMAND_RLIM; 
	sscanf(buffer, "%s %d %d %lf %lf %lf %lf %lf %lf %lf", cdummy, &i, &j, &rr, &x1, &y1, &e1, &e2, &x, &y);
	mockext2(i, j, rr, x1, y1, e1, e2, x, y);
      }
       
      if(strcmp(keyword, "mockext3") == 0){
	i = DEF_COMMAND_MOCK; j = DEF_COMMAND_EI;
	rr = DEF_COMMAND_FREC; 
	e1 = DEF_COMMAND_ELLIP; e2 = DEF_COMMAND_ELLIP;
	x = DEF_COMMAND_FTH; y = DEF_COMMAND_RLIM;  
	sscanf(buffer, "%s %d %d %lf %lf %lf %lf %lf", cdummy, &i, &j, &rr, &e1, &e2, &x, &y);
	mockext3(i, j, rr, e1, e2, x, y);
      }

      /*------------------------------------
	reading (writing) opt-related data
      */
      if(strcmp(keyword, "readobs_extend") == 0){
	i = sscanf(buffer, "%s %s %s", cdummy, fname, fname2);
	readobs_extend(fname, 1);
	if(i == 3) readmask(fname2, 1);
	calc_obsnoise();
      }
      
      if(strcmp(keyword, "readnoise_extend") == 0){
	i = sscanf(buffer, "%s %s", cdummy, fname);
	readnoise_extend(fname, 1);
      }
      
      if(strcmp(keyword, "writenoise") == 0){
	sscanf(buffer, "%s", cdummy);
	writenoise();
      }
      
      if(strcmp(keyword, "addnoise") == 0){
	sscanf(buffer, "%s %lf %lf", cdummy, &x1, &x2);
	addnoise(x1, x2);
      }
      
      if(strcmp(keyword, "readobs_point") == 0){
	sscanf(buffer, "%s %s", cdummy, fname);
	readobs_point(fname, 1);
      }
      
      if(strcmp(keyword, "parprior") == 0){
	sscanf(buffer, "%s %s", cdummy, fname);
	parprior(fname, 1);
      }
      
      if(strcmp(keyword, "mapprior") == 0){
	sscanf(buffer, "%s %s", cdummy, fname);
	mapprior(fname, 1);
      }

      /*------------------------------------
	optimizations
      */
      if(strcmp(keyword, "optpoint") == 0){
	sscanf(buffer, "%s", cdummy);
	chi2calc_opt_point_out();
      }
      
      if(strcmp(keyword, "optextend") == 0){
	sscanf(buffer, "%s", cdummy);
	chi2calc_opt_extend_out();
      }
      
      if(strcmp(keyword, "optimize") == 0){
	sscanf(buffer, "%s %s", cdummy, keyword);
	f = command_srcflag(keyword);
	opt_lens(f, 1);
      }
      
      if(strcmp(keyword, "varyone") == 0){
	sscanf(buffer, "%s %d %d %lf %lf %d %s", cdummy, &i, &j, &x1, &x2, &n, keyword);
	 f = command_srcflag(keyword);
	 varyone(i - 1, j - 1, x1, x2, n, f);
      }

      if(strcmp(keyword, "varytwo") == 0){
	sscanf(buffer, "%s %d %d %lf %lf %d %d %d %lf %lf %d %s", cdummy, &i, &j, &x1, &x2, &n, &k, &l, &y1, &y2, &r, keyword);
	f = command_srcflag(keyword);
	varytwo(i - 1, j - 1, x1, x2, n, k - 1, l - 1, y1, y2, r, f);
      }

      if(strcmp(keyword, "varyzs_extend") == 0){
	sscanf(buffer, "%s %d %lf %lf %d %s", cdummy, &i, &x1, &x2, &n, keyword);
	f = command_srcflag(keyword);
	varyzs_extend(i - 1, x1, x2, n, f);
       }
      
      if(strcmp(keyword, "varyzs_point") == 0){
	sscanf(buffer, "%s %d %lf %lf %d %s", cdummy, &i, &x1, &x2, &n, keyword);
	f = command_srcflag(keyword);
	varyzs_point(i - 1, x1, x2, n, f);
      }

      if(strcmp(keyword, "varycosmo") == 0){
	sscanf(buffer, "%s %s %lf %lf %d %s", cdummy, fname, &x1, &x2, &n, keyword);
	f = command_srcflag(keyword);
	varycosmo(fname, x1, x2, n, f);
      }

      if(strcmp(keyword, "randomize") == 0){
	sscanf(buffer, "%s", cdummy);
	randomize(1);
      }

      if(strcmp(keyword, "opt_explore") == 0){
	i = DEF_COMMAND_MOCK;
	x1 = DEF_COMMAND_C2LIMIT;
	sscanf(buffer, "%s %d %lf %s", cdummy, &i, &x1, keyword);
	f = command_srcflag(keyword);
	opt_explore(i, x1, f);
      }

      /*------------------------------------
	Markov Chain Monte Carlo
      */
      if(strcmp(keyword, "mcmc_sigma") == 0){
	sscanf(buffer, "%s %s", cdummy, fname);
	mcmc_calc_init(fname);
      }
      
      if(strcmp(keyword, "mcmc") == 0){
	 i = DEF_COMMAND_MCMC;
	 sscanf(buffer, "%s %d", cdummy, &i);
	 mcmc_calc(i);
      }
      
      if(strcmp(keyword, "mcmc_kaprad") == 0){
	zs = DEF_COMMAND_ZS;
	x = DEF_COMMAND_XY; y = DEF_COMMAND_XY;
	r1 = DEF_COMMAND_R1; r2 = DEF_COMMAND_R2; n = DEF_COMMAND_N;
	i = DEF_COMMAND_I;
	sscanf(buffer, "%s %s %d %lf %lf %lf %lf %lf %d %d", cdummy, fname, &j, &zs, &x, &y, &r1, &r2, &n, &i);
	mcmc_out_kappa_rad(fname, j, zs, x, y, r1, r2, n, i);
      }

      if(strcmp(keyword, "mcmc_ein") == 0){
	zs = DEF_COMMAND_ZS;
	i = DEF_COMMAND_I;
	sscanf(buffer, "%s %s %d %lf %d", cdummy, fname, &j, &zs, &i);
	mcmc_out_ein(fname, j, zs, i);
      }
      
      if(strcmp(keyword, "mcmc_ein2") == 0){
	zs = DEF_COMMAND_ZS;
	i = DEF_COMMAND_I;
	x = DEF_COMMAND_XY; y = DEF_COMMAND_XY; 
	sscanf(buffer, "%s %s %d %lf %lf %lf %d", cdummy, fname, &j, &zs, &x, &y, &i);
	 mcmc_out_ein2(fname, j, zs, x, y, i);
      }
      
      if(strcmp(keyword, "mcmc_calcim") == 0){
	zs = DEF_COMMAND_ZS;
	i = DEF_COMMAND_I;
	x = DEF_COMMAND_XY; y = DEF_COMMAND_XY; 
	sscanf(buffer, "%s %s %d %lf %lf %lf", cdummy, fname, &j, &zs, &x, &y);
	mcmc_out_calcim(fname, j, zs, x, y);
      }
      
      /*------------------------------------
	reset parameters
      */
      if(strcmp(keyword, "reset_par") == 0){
	strncpy(cdummy, buffer + 9, strlen(buffer) - 9);
	strcpy(buffer, cdummy);
	sscanf(buffer, "%s", keyword);
	fprintf(stderr, "######## re-setting parameter \n\n");
	
	init_para_body(keyword, buffer, 1);
	init_para2_body(keyword, buffer, 1);
	fprintf(stderr, "\n");
	
	if(strcmp(keyword, "ran_seed") == 0){
	  gsl_rng_set(ran_gsl, ran_seed);
	}
	
	if(strcmp(keyword, "galfile") == 0) readgals();
	if(strcmp(keyword, "srcfile") == 0) readsrcs();
	
	if(strcmp(keyword, "srcsbinsize") == 0) unset_tab_calc_src();
	
	if(check_para_cosmo() > 0) 
	  terminator("invalid input parameter (reset_par)"); 
	
	if((strcmp(keyword, "omega") == 0) || (strcmp(keyword, "lambda") == 0)
	   || (strcmp(keyword, "weos") == 0) || (strcmp(keyword, "hubble") == 0)
	   || (strcmp(keyword, "xmin") == 0) || (strcmp(keyword, "xmax") == 0)
	   || (strcmp(keyword, "ymin") == 0) || (strcmp(keyword, "ymax") == 0)
	   || (strcmp(keyword, "pix_ext") == 0) || (strcmp(keyword, "pix_poi") == 0)
	   || (strcmp(keyword, "maxlev") == 0)){
	  nx_ext = (int)((xmax - xmin + NPIX_SMALL_OFFSET) / pix_ext);
	  ny_ext = (int)((ymax - ymin + NPIX_SMALL_OFFSET) / pix_ext);
	  nx_poi = (int)((xmax - xmin + NPIX_SMALL_OFFSET) / pix_poi);
	  ny_poi = (int)((ymax - ymin + NPIX_SMALL_OFFSET) / pix_poi);
	  set_distance_lpl_init();
	  ext_unset_table();
	  poi_unset_table();
	  obs_unset_table();
	  poimg_unset_table();
	}
      }
      
      if(strcmp(keyword, "reset_lens") == 0){
	i = DEF_COMMAND_RESET_I; j = DEF_COMMAND_RESET_J; x = DEF_COMMAND_RESET_P;
	sscanf(buffer, "%s %d %d %lf", cdummy, &i, &j, &x);
	fprintf(stderr, "######## re-setting lens parameter \n\n");
	fprintf(stderr, "lens id = %d,  par no. = %d,  value = %e\n\n", i, j, x);
	
	if((i > num_len) || (j > NPAR_LEN) || (i <= 0) || (j <= 0))
	  terminator("id irrelevant (reset_lens)"); 
	if(check_para_lens_all() > 0)
	  terminator("invalid input parameter (reset_lens)");
	 
	para_lens[i - 1][j - 1] = x;
	if(j == 1) set_distance_lpl_init();
	poi_unset_table();
	ext_unset_table();
	poimg_unset_table();
      }
       
      if(strcmp(keyword, "reset_extend") == 0){
	i = DEF_COMMAND_RESET_I; j = DEF_COMMAND_RESET_J; x = DEF_COMMAND_RESET_P;
	sscanf(buffer, "%s %d %d %lf", cdummy, &i, &j, &x);
	fprintf(stderr, "######## re-setting extend parameter \n\n");
	fprintf(stderr, "extend id = %d,  par no. = %d,  value = %e\n\n", i, j, x);
	 
	if((i > num_ext) || (j > NPAR_EXT) || (i <= 0) || (j <= 0))
	  terminator("id irrelevant (reset_extend)"); 
	if((check_para_lens_all() > 0) || (check_para_ext_all() > 0))
	  terminator("invalid input parameter (reset_extend)");
	
	para_ext[i - 1][j - 1] = x;
	if(j == 1) ext_unset_table();
      }
       
      if(strcmp(keyword, "reset_point") == 0){
	i = DEF_COMMAND_RESET_I; j = DEF_COMMAND_RESET_J; x = DEF_COMMAND_RESET_P;
	sscanf(buffer, "%s %d %d %lf", cdummy, &i, &j, &x);
	fprintf(stderr, "######## re-setting point parameter \n\n");
	fprintf(stderr, "point id = %d,  par no. = %d,  value = %e\n\n", i, j, x);
	 
	if((i > num_poi) || (j > NPAR_POI) || (i <= 0) || (j <= 0))
	  terminator("id irrelevant (reset_point)"); 
	if((check_para_lens_all() > 0) || (check_para_poi_all() > 0))
	  terminator("invalid input parameter (reset_point)");
	
	para_poi[i - 1][j - 1] = x;
	if(j == 1){
	  poi_unset_table();
	  poimg_unset_table();
	}
      }
       
      if(strcmp(keyword, "reset_psf") == 0){
	j = DEF_COMMAND_RESET_J; x = DEF_COMMAND_RESET_P;
	sscanf(buffer, "%s %d %lf", cdummy, &j, &x);
	fprintf(stderr, "######## re-setting psf parameter \n\n");
	fprintf(stderr, "psf par no. = %d,  value = %e\n\n", j, x);

	if((flag_seeing == 0) || (j > NPAR_PSF) || (j <= 0))
	  terminator("id irrelevant (reset_psf)"); 
	if(check_para_ext_all() > 0)
	  terminator("invalid input parameter (reset_psf)");
	
	para_psf[j - 1] = x;
      }
       
      if(strcmp(keyword, "resetopt_lens") == 0){
	i = DEF_COMMAND_RESET_I; j = DEF_COMMAND_RESET_J; k = DEF_COMMAND_RESET_F;
	sscanf(buffer, "%s %d %d %d", cdummy, &i, &j, &k);
	fprintf(stderr, "######## re-setting lens optimization flag\n\n");
	fprintf(stderr, "lens id = %d,  par no. = %d,  value = %d\n\n", i, j, k);
	
	if((i > num_len) || (j > NPAR_LEN) || (i <= 0) || (j <= 0))
	  terminator("id irrelevant (resetopt_lens)"); 
	 
	flag_para_lens[i - 1][j - 1] = k;
      }

      if(strcmp(keyword, "resetopt_extend") == 0){
	i = DEF_COMMAND_RESET_I; j = DEF_COMMAND_RESET_J; k = DEF_COMMAND_RESET_F;
	sscanf(buffer, "%s %d %d %d", cdummy, &i, &j, &k);
	fprintf(stderr, "######## re-setting extend optimization flag\n\n");
	fprintf(stderr, "extend id = %d,  par no. = %d,  value = %d\n\n", i, j, k);
	
	if((i > num_ext) || (j > NPAR_EXT) || (i <= 0) || (j <= 0))
	  terminator("id irrelevant (resetopt_extend)"); 
	
	flag_para_ext[i - 1][j - 1] = k;
      }

      if(strcmp(keyword, "resetopt_point") == 0){
	i = DEF_COMMAND_RESET_I; j = DEF_COMMAND_RESET_J; k = DEF_COMMAND_RESET_F;
	sscanf(buffer, "%s %d %d %d", cdummy, &i, &j, &k);
	fprintf(stderr, "######## re-setting point optimization flag\n\n");
	fprintf(stderr, "point id = %d,  par no. = %d,  value = %d\n\n", i, j, k);
	 
	if((i > num_poi) || (j > NPAR_POI) || (i <= 0) || (j <= 0))
	  terminator("id irrelevant (resetopt_point)"); 
	
	flag_para_poi[i - 1][j - 1] = k;
      }

      if(strcmp(keyword, "resetopt_psf") == 0){
	j = DEF_COMMAND_RESET_J; k = DEF_COMMAND_RESET_F;
	sscanf(buffer, "%s %d %d", cdummy, &j, &k);
	fprintf(stderr, "######## re-setting psf optimization flag\n\n");
	fprintf(stderr, "psf par no. = %d,  value = %d\n\n", j, k);
	
	if((flag_seeing == 0) || (j > NPAR_PSF) || (j <= 0))
	  terminator("id irrelevant (resetopt_psf)"); 
	
	flag_para_psf[j - 1] = k;
      }

      if(strcmp(keyword, "mv_lens") == 0){
	i = DEF_COMMAND_RESET_I; x = DEF_COMMAND_RESET_P; y = DEF_COMMAND_RESET_P;
	sscanf(buffer, "%s %d %lf %lf", cdummy, &i, &x, &y);
	fprintf(stderr, "######## re-setting lens position \n\n");
	fprintf(stderr, "lens id = %d,  x = %e,  y = %e\n\n", i, x, y);
	
	if((i > num_len) || (i <= 0))
	  terminator("id irrelevant (mv_lens)"); 
	
	para_lens[i - 1][2] = x;
	para_lens[i - 1][3] = y;
	poi_unset_table();
	ext_unset_table();
	poimg_unset_table();

	if(check_para_lens_all() > 0)
	  terminator("invalid input parameter (mv_lens)");
      }

      if(strcmp(keyword, "mv_extend") == 0){
	i = DEF_COMMAND_RESET_I; x = DEF_COMMAND_RESET_P; y = DEF_COMMAND_RESET_P;
	sscanf(buffer, "%s %d %lf %lf", cdummy, &i, &x, &y);
	fprintf(stderr, "######## re-setting extend position \n\n");
	fprintf(stderr, "extend id = %d,  x = %e,  y = %e\n\n", i, x, y);
	
	if((i > num_ext) || (i <= 0))
	  terminator("id irrelevant (mv_extend)"); 
	 
	para_ext[i - 1][2] = x;
	para_ext[i - 1][3] = y;

	if(check_para_ext_all() > 0)
	  terminator("invalid input parameter (mv_extend)");
      }

      if(strcmp(keyword, "mv_point") == 0){
	i = DEF_COMMAND_RESET_I; x = DEF_COMMAND_RESET_P; y = DEF_COMMAND_RESET_P;
	sscanf(buffer, "%s %d %lf %lf", cdummy, &i, &x, &y);
	fprintf(stderr, "######## re-setting point position \n\n");
	fprintf(stderr, "point id = %d,  x = %e,  y = %e\n\n", i, x, y);
	
	if((i > num_poi) || (i <= 0))
	  terminator("id irrelevant (mv_point)"); 
	
	para_poi[i - 1][1] = x;
	para_poi[i - 1][2] = y;
	 
	if(check_para_poi_all() > 0)
	  terminator("invalid input parameter (mv_point)");
      }

      /*------------------------------------
        some utils
      */
      if(strcmp(keyword, "printmodel") == 0){
	fprintf(stderr, "######## printing current model\n");
	dump_model(stderr);
	fprintf(stderr, "\n");
      }
      
      if(strcmp(keyword, "printopt") == 0){
	fprintf(stderr, "######## printing flags for optimization\n");
	dump_opt_flag(stderr);
	fprintf(stderr, "\n");
      }

      if(strcmp(keyword, "printlensplane") == 0){
	fprintf(stderr, "######## printing lens plane\n");
	gen_lensplane(1);
	fprintf(stderr, "\n");
      }

      if(strcmp(keyword, "arcsec2mpc") == 0){
	zs = DEF_COMMAND_ZL;
	x = DEF_COMMAND_R1;
	sscanf(buffer, "%s %lf %lf", cdummy, &zs, &x);
	fprintf(stderr, "######## arcsec to Mpc/h for z = %e\n\n", zs);
	fprintf(stderr, "%e [arcsec] = %e [Mpc/h]\n", x, thetator_dis(x, dis_angulard(0.0, zs)));
	fprintf(stderr, "\n");
      }

      if(strcmp(keyword, "mpc2arcsec") == 0){
	zs = DEF_COMMAND_ZL;
	x = DEF_COMMAND_R1;
	sscanf(buffer, "%s %lf %lf", cdummy, &zs, &x);
	fprintf(stderr, "######## Mpc/h to arcsec for z = %e\n\n", zs);
	fprintf(stderr, "%e [Mpc/h] = %e [arcsec]\n", x, rtotheta_dis(x, dis_angulard(0.0, zs)));
	fprintf(stderr, "\n");
      }

      if(strcmp(keyword, "critdens") == 0){
	zs = DEF_COMMAND_ZS;
	zl = DEF_COMMAND_ZL;
	sscanf(buffer, "%s %lf %lf", cdummy, &zl, &zs);
	fprintf(stderr, "######## critical density for  zl = %e and zs = %e \n\n", zl, zs);
	
	if((zl <= 0.0) || (zs <= 0.0) || (zl >= zs))
	  terminator("redshift irrelevant (critdens)"); 

	fprintf(stderr, "Sigma_crit = %e [h M_Sun Mpc^-2]\n", sigma_crit_dis(dis_angulard(0.0, zs), dis_angulard(0.0, zl), dis_angulard(zl, zs)));
	fprintf(stderr, "\n");
      }
      
      if(strcmp(keyword, "dismod") == 0){
	zs = DEF_COMMAND_ZL;
	sscanf(buffer, "%s %lf", cdummy, &zs);
	fprintf(stderr, "######## distance modulus for zs = %e\n\n", zs);
	fprintf(stderr, "DM = %e\n", dis_mod(zs) - 5.0 * log10(hubble));
	fprintf(stderr, "DM = %e [hubble=1]\n", dis_mod(zs));
	fprintf(stderr, "\n");
      }

      /*------------------------------------
	quit
      */
      if(strcmp(keyword, "quit") == 0){
	r = 1;
      }
      
      /*------------------------------------
       */
    }
  }
  
  return r;
}

int command_srcflag(char *keyword)
{
  int f = 0;

  if(strcmp(keyword, "extend") == 0) f = 1;
  if(strcmp(keyword, "point") == 0) f = 2;

  return f;
}

/*--------------------------------------------------------------
  interactive mode
*/

void interactive(void)
{
  int r;
  char buffer[INPUT_MAXCHAR];

  do{
    fprintf(stderr, "glafic> ");
    fgets(buffer, INPUT_MAXCHAR, stdin);
    fprintf(stderr, "\n");
    r = do_command(buffer);
  }while(r == 0);

  return;
}

/*--------------------------------------------------------------
  for debug
*/

void deb(void)
{
  return;
}

