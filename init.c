#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

#define PAREADXE(s1, s2, f)						\
  if(strcmp(keyword, #s1) == 0){ sscanf(buffer, "%s %lf", cdummy, &s2); \
     if(f == 1){ fprintf(stderr, "%-15s =  %e\n", cdummy, s2); } }
#define PAREADXD(s1, s2, f)					       \
  if(strcmp(keyword, #s1) == 0){ sscanf(buffer, "%s %d", cdummy, &s2); \
     if(f == 1){ fprintf(stderr, "%-15s =  %d\n", cdummy, s2); } }
#define PAREADXS(s1, s2, f)					      \
  if(strcmp(keyword, #s1) == 0){ sscanf(buffer, "%s %s", cdummy, s2); \
     if(f == 1){ fprintf(stderr, "%-15s =  %s\n", cdummy, s2); } }

static int order_opt[NMAX_LEN + NMAX_EXT + NMAX_POI + 2];

/*--------------------------------------------------------------
  read input file
*/

void init_para(char *infile, int verb)
{
  char buffer[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  FILE* fptr;
  
  def_parameters();
  
  fptr = fopen(infile, "r");
  
  if(fptr == NULL) terminator("failed at fopen (input file)");
  
  while(fgets(buffer, INPUT_MAXCHAR, fptr)){
    if(sscanf(buffer, "%s", keyword) != EOF){
      if(keyword[0] != '#'){

	init_para_body(keyword, buffer, 0);

	/* error if version 1 input file is fed */
	if(strcmp(keyword, "zl") == 0)
	  terminator("wrong input file (primary parameter zl not used in ver. 2)");
	  
	if(strcmp(keyword, "startup") == 0) break;

      }
    }
  }
  
  fclose(fptr);
  
  set_npix();
  if(verb == 1) out_para();
  
  return;
}

void set_npix(void)
{
  nx_ext = (int)((xmax - xmin + NPIX_SMALL_OFFSET) / pix_ext);
  ny_ext = (int)((ymax - ymin + NPIX_SMALL_OFFSET) / pix_ext);
  /* xmax = xmin + ((double)nx_ext) * pix_ext;
     ymax = ymin + ((double)ny_ext) * pix_ext; */
  
  nx_poi = (int)((xmax - xmin + NPIX_SMALL_OFFSET) / pix_poi);
  ny_poi = (int)((ymax - ymin + NPIX_SMALL_OFFSET) / pix_poi);
  /* xmax = xmin + ((double)nx_poi) * pix_poi;
     ymax = ymin + ((double)ny_poi) * pix_poi; */
  
  if(maxlev < 1) maxlev = 1;
  if(maxlev > NMAX_MAXLEV) maxlev = NMAX_MAXLEV;

  if((nx_ext > NMAX_PIXEL) || (ny_ext > NMAX_PIXEL) || ((nx_poi * ny_poi) > NMAX_PIXEL_POINT)){
    terminator("pixel number exceeds the limit");
  }

  if((nx_ext <= 0) || (ny_ext <= 0) || (nx_poi <= 0) || (ny_poi <= 0)){
    terminator("pixel number invalid");
  }

  return;
}

void init_para_body(char *keyword, char *buffer, int verb)
{
  char cdummy[INPUT_MAXCHAR];

  PAREADXE(omega, omega, verb)
  PAREADXE(lambda, lambda, verb)
  PAREADXE(weos, weos, verb)
  PAREADXE(hubble, hubble, verb)
  PAREADXS(prefix, file_prefix, verb)
  PAREADXE(xmin, xmin, verb)
  PAREADXE(ymin, ymin, verb)
  PAREADXE(xmax, xmax, verb)
  PAREADXE(ymax, ymax, verb)
  PAREADXE(pix_ext, pix_ext, verb)
  PAREADXE(pix_poi, pix_poi, verb)
  PAREADXD(maxlev, maxlev, verb)

  return;
}

void init_para2(char *infile, int verb)
{
  char buffer[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  FILE* fptr;

  fptr = fopen(infile, "r");

  if(fptr == NULL) terminator("failed at fopen (input file)");
  
  if(verb == 1) fprintf(stderr, "\n");

  while(fgets(buffer, INPUT_MAXCHAR, fptr)){
    if(sscanf(buffer, "%s", keyword) != EOF){
      if(keyword[0] != '#'){
	
	/* secondary parameters */
	init_para2_body(keyword, buffer, verb);

	/* end of parameter setting */
	if(strcmp(keyword, "startup") == 0) break;

      }
    }
  }
  
  fclose(fptr);

  if(skyfix_value < (0.1 * DEF_SKYFIX_VALUE)) skymed = skyfix_value;

  return;
}

void init_para2_body(char *keyword, char *buffer, int verb)
{
  char cdummy[INPUT_MAXCHAR];

  PAREADXD(ovary, ovary, verb)
  PAREADXD(lvary, lvary, verb)
  PAREADXD(wvary, wvary, verb)
  PAREADXD(hvary, hvary, verb)

  PAREADXS(galfile, file_gal, verb)
  PAREADXS(srcfile, file_src, verb)
  PAREADXD(ran_seed, ran_seed, verb)
  PAREADXD(flag_hodensity, flag_hodensity, verb)
  PAREADXE(hodensity, hodensity, verb)
  PAREADXD(gnfw_usetab, gnfw_usetab, verb)
  PAREADXD(ein_usetab, ein_usetab, verb)
  PAREADXD(nfw_users, nfw_users, verb)
  PAREADXD(nmax_poi_ite, nmax_poi_ite, verb)
  PAREADXE(max_poi_tol, max_poi_tol, verb)
  PAREADXE(poi_imag_max, poi_imag_max, verb)
  PAREADXE(poi_imag_min, poi_imag_min, verb)
  PAREADXD(ang_step, center_ang_step, verb)
  PAREADXE(imag_ceil, imag_ceil, verb)
  PAREADXE(smallcore, smallcore, verb)
  PAREADXD(outformat_exp, outformat_exp, verb)
  PAREADXD(flatfix, flatfix, verb)
  PAREADXD(flag_outpot, flag_outpot, verb)
  
  PAREADXE(amoeba_tol1, tol_amoeba_lens, verb)
  PAREADXE(amoeba_tol2, tol_amoeba, verb)
  PAREADXD(nmax_amoeba, nmax_amoeba, verb)
  PAREADXD(nmax_amoeba_point, nmax_amoeba_point, verb)
  PAREADXE(amoeba_dp_mass, amoeba_dp_mass, verb)
  PAREADXE(amoeba_dp_xy, amoeba_dp_xy, verb)
  PAREADXE(amoeba_dp_e, amoeba_dp_e, verb)
  PAREADXE(amoeba_dp_ang, amoeba_dp_ang, verb)
  PAREADXE(amoeba_dp_r, amoeba_dp_r, verb)
  PAREADXE(amoeba_dp_n, amoeba_dp_n, verb)
  PAREADXE(amoeba_dp_z, amoeba_dp_z, verb)
  PAREADXE(amoeba_dp_cos, amoeba_dp_cosmo, verb)
  PAREADXE(amoeba_delmin, amoeba_delmin, verb)
  PAREADXE(amoeba_delmax, amoeba_delmax, verb)

  PAREADXE(amoeba_dp_psfw, amoeba_dp_psfw, verb)
  PAREADXE(amoeba_dp_psfe, amoeba_dp_psfe, verb)
  PAREADXE(amoeba_dp_psfpa, amoeba_dp_psfpa, verb)
  PAREADXE(amoeba_dp_psfb, amoeba_dp_psfb, verb)
  PAREADXE(amoeba_dp_psff, amoeba_dp_psff, verb)

  PAREADXD(chi2_splane, chi2_point_splane, verb)
  PAREADXD(chi2_checknimg, chi2_checknimg, verb)
  PAREADXD(chi2_usemag, chi2_usemag, verb)
  PAREADXD(chi2_restart, chi2_restart, verb)
  PAREADXE(chi2pen_range, chi2pen_range, verb)
  PAREADXE(chi2pen_nimg, chi2pen_nimg, verb)
  PAREADXE(chi2pen_parity, chi2pen_parity, verb)
  PAREADXD(chi2_restart_max, chi2_restart_max, verb)
  PAREADXE(obs_gain, obs_gain, verb)
  PAREADXD(obs_ncomb, obs_ncomb, verb)
  PAREADXE(obs_readnoise, obs_readnoise, verb)
  PAREADXE(noise_clip, noise_clip, verb)
  PAREADXD(skyfix, skyfix, verb)
  PAREADXE(skyfix_value, skyfix_value, verb)
  PAREADXE(psfconv_size, psfconv_size, verb)
  PAREADXD(seeing_sub, seeing_sub, verb)
  PAREADXE(bicub_a, bicub_a, verb)
  PAREADXE(source_calcr0, source_calcr0, verb)
  PAREADXD(flag_extref, flag_extref, verb)
  PAREADXE(source_refr0, source_refr0, verb)
  PAREADXD(num_pixint, num_pixint, verb)
  PAREADXD(flag_srcsbin, flag_srcsbin, verb)
  PAREADXE(srcsbinsize, srcsbinsize, verb)
  PAREADXD(flag_extnorm, flag_extnorm, verb)
  PAREADXD(nmax_srcs, nmax_srcs, verb)
  PAREADXD(nmax_fft, nmax_fft, verb)
  PAREADXD(addwcs, flag_addwcs, verb)
  PAREADXD(flag_mcmcall, flag_mcmcall, verb)
  PAREADXE(wcs_ra0, wcs_ra0, verb)
  PAREADXE(wcs_dec0, wcs_dec0, verb)
  PAREADXE(dr_lens_center, dr_lens_center, verb)
  PAREADXD(flag_pow_tm15, flag_pow_tm15, verb)
  PAREADXE(tol_pow_tm15, tol_pow_tm15, verb)
    
  return;
}

/*--------------------------------------------------------------
  startup
*/

void startup(char *infile, int verb)
{
  int n0, n1, n2, n3, n4;
  char buffer[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  char cdummy[INPUT_MAXCHAR];
  FILE* fptr;

  if(verb == 1) fprintf(stderr, "\n######## startup");

  fptr = fopen(infile, "r");

  if(fptr == NULL) terminator("failed at fopen (input file)");

  do{
    if(fgets(buffer, INPUT_MAXCHAR, fptr) == NULL)
      terminator("no startup found");
    sscanf(buffer, "%s", keyword);
  }while(strcmp(keyword, "startup") != 0);

  sscanf(buffer, "%s %d %d %d", cdummy, &num_len, &num_ext, &num_poi);
  
  if((num_len > NMAX_LEN) || (num_ext > NMAX_EXT) || (num_poi > NMAX_POI))
    terminator("number exceeds maximum");
  if(num_len == 0)
    terminator("need at least one lens");
  
  n0 = 0;
  n1 = 0;
  n2 = 0;
  n3 = 0;
  n4 = 0;

  do{
    fgets(buffer, INPUT_MAXCHAR, fptr);
    if(sscanf(buffer, "%s", keyword) != EOF){
      if(keyword[0] != '#'){
	
	if(strcmp(keyword, "lens") == 0){
	  order_opt[n0] = 1;
	  n0++;
	  startup_lens(buffer, n1);
	  n1++;
	}
	
	if(strcmp(keyword, "extend") == 0){
	  order_opt[n0] = 2;
	  n0++;
	  startup_extend(buffer, n2);
	  n2++;
	}
	
	if(strcmp(keyword, "point") == 0){
	  order_opt[n0] = 3;
	  n0++;
	  startup_point(buffer, n3);
	  n3++;
	}

	if(strcmp(keyword, "psf") == 0){
	  order_opt[n0] = 4;
	  n0++;
	  startup_psf(buffer);
	  n4++;
	}
      }
    }
  }while(strcmp(keyword, "end_startup") != 0);

  order_opt[n0] = 0;
  
  if((num_len != n1) || (num_ext != n2) || (num_poi != n3) || (n4 > 1)){
    terminator("startup failed (number mismatch)");
  }

  if((check_para_lens_all() > 0) || (check_para_ext_all() > 0) || (check_para_poi_all() > 0)){
    terminator("startup failed (invalid input parameter)");
  }

  fclose(fptr);

  if(verb == 1) dump_model(stderr);

  return;
}

void startup_lens(char *buffer, int n1)
{
  int nn;
  char cdummy[INPUT_MAXCHAR];
  char cdummy2[INPUT_MAXCHAR];

  nn = sscanf(buffer, "%s %s %lf %lf %lf %lf %lf %lf %lf %lf", cdummy, cdummy2, &para_lens[n1][0], &para_lens[n1][1], &para_lens[n1][2], &para_lens[n1][3], &para_lens[n1][4], &para_lens[n1][5], &para_lens[n1][6], &para_lens[n1][7]);
  if(nn != (NPAR_LEN + 2)) terminator("startup failed (invalid 'lens' format)"); 
  model_lens[n1] = lmodeltoint(cdummy2);
  if(model_lens[n1] == 0) terminator("invalid lens model name"); 
  if(strcmp(cdummy2, "gals") == 0) readgals(); 

  return;
}

void startup_extend(char *buffer, int n2)
{
  int nn;
  char cdummy[INPUT_MAXCHAR];
  char cdummy2[INPUT_MAXCHAR];

  nn = sscanf(buffer, "%s %s %lf %lf %lf %lf %lf %lf %lf %lf", cdummy, cdummy2, &para_ext[n2][0], &para_ext[n2][1], &para_ext[n2][2], &para_ext[n2][3], &para_ext[n2][4], &para_ext[n2][5], &para_ext[n2][6], &para_ext[n2][7]);
  if(nn != (NPAR_EXT + 2)) terminator("startup failed (invalid 'extend' format)"); 
  model_ext[n2] = emodeltoint(cdummy2);
  if(model_ext[n2] == 0) terminator("invalid extend model name"); 
  if(strcmp(cdummy2, "srcs") == 0) readsrcs(); 

  return;
}

void startup_point(char *buffer, int n3)
{
  int nn;
  char cdummy[INPUT_MAXCHAR];

  nn = sscanf(buffer, "%s %lf %lf %lf", cdummy, &para_poi[n3][0], &para_poi[n3][1], &para_poi[n3][2]);
  if(nn != (NPAR_POI + 1)) terminator("startup failed (invalid 'point' format)"); 

  return;
}

void startup_psf(char *buffer)
{
  int nn;
  char cdummy[INPUT_MAXCHAR];

  nn = sscanf(buffer, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf", cdummy, &para_psf[0], &para_psf[1], &para_psf[2], &para_psf[3], &para_psf[4], &para_psf[5], &para_psf[6], &para_psf[7], &para_psf[8]);
  if(nn != (NPAR_PSF + 1)) terminator("startup failed (invalid 'psf' format)"); 
  
  flag_seeing = 1;

  return;
}

/*--------------------------------------------------------------
  gen_lensplane
*/

void gen_lensplane(int verb)
{
  int i, j, k, flag;
  double zl_tmp;
  
  num_lpl = 0;
  zl_lpl[num_lpl] = para_lens[0][0];
  lens_lpl_id[0] = num_lpl;
  num_lpl++;
  for(i=1;i<num_len;i++){
    flag = 0;
    for(j=0;j<i;j++){
      if(fabs(para_lens[i][0] - para_lens[j][0]) < TOL_ZLPL){
	lens_lpl_id[i] = lens_lpl_id[j];
	flag = 1;
	break;
      }
    }
    if(flag == 0){
      if(num_lpl >= NMAX_LPL) terminator("too many lens planes");
      zl_lpl[num_lpl] = para_lens[i][0];
      lens_lpl_id[i] = num_lpl;
      num_lpl++;
    }
  }

  for(i=0;i<(num_lpl-1);i++){
    for(j=i+1;j<num_lpl;j++){
      if(zl_lpl[i] > zl_lpl[j]){
	zl_tmp = zl_lpl[i];
	zl_lpl[i] = zl_lpl[j];
	zl_lpl[j] = zl_tmp;
	for(k=0;k<num_len;k++){
	  if(lens_lpl_id[k] == i){
	    lens_lpl_id[k] = j;
	  } else if(lens_lpl_id[k] == j){
	    lens_lpl_id[k] = i;
	  }
	}
      }
    }
  }

  if(verb == 1) dump_lensplane(stderr);

  return;
}

/*--------------------------------------------------------------
  setopt
*/

void setopt(char *infile, int verb)
{
  int i, j, n0, n1, n2, n3, n4;
  char buffer[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  FILE* fptr;

  fptr = fopen(infile, "r");

  if(fptr == NULL) terminator("failed at fopen (input file)");

  for(j=0;j<NPAR_LEN;j++){
    for(i=0;i<NMAX_LEN;i++) flag_para_lens[i][j] = 0; 
  }
  for(j=0;j<NPAR_EXT;j++){
    for(i=0;i<NMAX_EXT;i++) flag_para_ext[i][j] = 0; 
  }
  for(j=0;j<NPAR_POITAB;j++){
    for(i=0;i<NMAX_POI;i++) flag_para_poi[i][j] = 0; 
  }
  for(j=0;j<NPAR_PSF;j++){
    flag_para_psf[j] = 0;
  }

  while(fgets(buffer, INPUT_MAXCHAR, fptr)){
    sscanf(buffer, "%s", keyword);
    if(strcmp(keyword, "start_setopt") == 0){ 
      
      n0 = 0;
      n1 = 0;
      n2 = 0;
      n3 = 0;
      n4 = 0;
      
      do{
	fgets(buffer, INPUT_MAXCHAR, fptr);
	if(sscanf(buffer, "%s", keyword) != EOF){
	  if(keyword[0] != '#'){

	    switch(order_opt[n0]){
	      
	    case 1:
	      setopt_lens(buffer, n1);
	      n1++;
	      n0++;
	      break;

	    case 2:
	      setopt_extend(buffer, n2);
	      n2++;
	      n0++;
	      break;
	      
	    case 3:
	      setopt_point(buffer, n3);
	      n3++;
	      n0++;
	      break;
	      
	    case 4:
	      setopt_psf(buffer);
	      n4++;
	      n0++;
	      break;
	    }
	  }
	}
      }while(strcmp(keyword, "end_setopt") != 0);
      
      if((num_len != n1) || (num_ext != n2) || (num_poi != n3) || (n4 != flag_seeing)){
	terminator("setopt failed (number mismatch)");
      }

    }
  }

  fclose(fptr);

  if(verb == 1) dump_opt_flag(stderr);

  return;
}

void setopt_lens(char *buffer, int n1)
{
  int nn;
  
  nn = sscanf(buffer, "%d %d %d %d %d %d %d %d", &flag_para_lens[n1][0], &flag_para_lens[n1][1], &flag_para_lens[n1][2], &flag_para_lens[n1][3], &flag_para_lens[n1][4], &flag_para_lens[n1][5], &flag_para_lens[n1][6], &flag_para_lens[n1][7]);
  if(nn != NPAR_LEN) terminator("setopt failed (invalid format)"); 

  return;
}

void setopt_extend(char *buffer, int n2)
{
  int nn;

  nn = sscanf(buffer, "%d %d %d %d %d %d %d %d", &flag_para_ext[n2][0], &flag_para_ext[n2][1], &flag_para_ext[n2][2], &flag_para_ext[n2][3], &flag_para_ext[n2][4], &flag_para_ext[n2][5], &flag_para_ext[n2][6], &flag_para_ext[n2][7]);
  if(nn != NPAR_EXT) terminator("setopt failed (invalid format)"); 

  return;
}

void setopt_point(char *buffer, int n3)
{
  int nn;

  nn = sscanf(buffer, "%d %d %d", &flag_para_poi[n3][0], &flag_para_poi[n3][1], &flag_para_poi[n3][2]);
  if(nn != NPAR_POI) terminator("setopt failed (invalid format)"); 
  
  return;
}

void setopt_psf(char *buffer)
{
  int nn;

  nn = sscanf(buffer, "%d %d %d %d %d %d %d %d %d", &flag_para_psf[0], &flag_para_psf[1], &flag_para_psf[2], &flag_para_psf[3], &flag_para_psf[4], &flag_para_psf[5], &flag_para_psf[6], &flag_para_psf[7], &flag_para_psf[8]);
  if(nn != NPAR_PSF) terminator("setopt failed (invalid format)"); 
  
  return;
}

/*--------------------------------------------------------------
  default parameters
*/

void def_parameters(void)
{
  int i, j, k;

  omega = DEF_OMEGA;
  lambda = DEF_LAMBDA;
  weos = DEF_WEOS;
  hubble = DEF_HUBBLE;
  sprintf(file_prefix, DEF_PREFIX);
  xmin = DEF_XMIN;
  ymin = DEF_YMIN;
  xmax = DEF_XMAX;
  ymax = DEF_YMAX;
  pix_ext = DEF_PIX_EXT;
  pix_poi = DEF_PIX_POI;
  maxlev = DEF_MAXLEV;

  ovary = lvary = wvary = hvary = 0;
  o_error = l_error = w_error = h_error = 0.0;
  omedian = omega;
  lmedian = lambda;
  wmedian = weos;
  hmedian = hubble;

  sprintf(file_gal, DEF_GALFILE);
  sprintf(file_src, DEF_SRCFILE);
  ran_seed = DEF_RAN_SEED;
  flag_hodensity = DEF_FLAG_HODENSITY;
  hodensity = DEF_HODENSITY;
  gnfw_usetab = DEF_GNFW_USETAB;
  nfw_users = DEF_NFW_USERS;
  ein_usetab = DEF_EIN_USETAB;
  nmax_poi_ite = DEF_NMAX_POI_ITE;
  max_poi_tol = DEF_MAX_POI_TOL;
  poi_imag_max = DEF_POI_IMAG_MAX;
  poi_imag_min = DEF_POI_IMAG_MIN;
  center_ang_step = DEF_CENTER_ANG_STEP;
  imag_ceil = DEF_IMAG_CEIL;
  smallcore = DEF_SMALLCORE;
  outformat_exp = DEF_OUTFORMAT_EXP;
  flatfix = DEF_FLATFIX;
  flag_outpot = DEF_FLAG_OUTPOT;
  
  tol_amoeba_lens = DEF_TOL_AMOEBA_LENS;
  tol_amoeba = DEF_TOL_AMOEBA;
  nmax_amoeba = DEF_NMAX_AMOEBA;
  nmax_amoeba_point = DEF_NMAX_AMOEBA_POINT;
  amoeba_dp_mass = DEF_AMOEBA_DP_MASS;
  amoeba_dp_xy = DEF_AMOEBA_DP_XY;
  amoeba_dp_e = DEF_AMOEBA_DP_E;
  amoeba_dp_ang = DEF_AMOEBA_DP_ANG;
  amoeba_dp_r = DEF_AMOEBA_DP_R;
  amoeba_dp_n = DEF_AMOEBA_DP_N;
  amoeba_dp_z = DEF_AMOEBA_DP_Z;
  amoeba_dp_cosmo = DEF_AMOEBA_DP_COSMO;
  amoeba_delmin = DEF_AMOEBA_DELMIN;
  amoeba_delmax = DEF_AMOEBA_DELMAX;

  amoeba_dp_psfw = DEF_AMOEBA_DP_PSFW;
  amoeba_dp_psfe = DEF_AMOEBA_DP_PSFE;
  amoeba_dp_psfpa = DEF_AMOEBA_DP_PSFPA;
  amoeba_dp_psfb = DEF_AMOEBA_DP_PSFB;
  amoeba_dp_psff = DEF_AMOEBA_DP_PSFF;

  chi2_point_splane = DEF_CHI2_POINT_SPLANE;
  chi2_checknimg = DEF_CHI2_CHECKNIMG;
  chi2_usemag = DEF_CHI2_USEMAG;
  chi2_restart = DEF_CHI2_RESTART;
  chi2pen_range = DEF_CHI2PEN_RANGE;
  chi2pen_nimg = DEF_CHI2PEN_NIMG;
  chi2pen_parity = DEF_CHI2PEN_PARITY;
  chi2_restart_max = DEF_CHI2_RESTART_MAX;
  obs_gain = DEF_OBS_GAIN;
  obs_ncomb = DEF_OBS_NCOMB;
  obs_readnoise = DEF_OBS_READNOISE;
  noise_clip = DEF_NOISE_CLIP;
  skyfix = DEF_SKYFIX;
  skyfix_value = DEF_SKYFIX_VALUE;
  psfconv_size = DEF_PSFCONV_SIZE;
  seeing_sub = DEF_SEEING_SUB;
  bicub_a = DEF_BICUB_A;
  source_calcr0 = DEF_SOURCE_CALCR0;
  flag_extref = DEF_FLAG_EXTREF;
  source_refr0 = DEF_SOURCE_REFR0;
  num_pixint = DEF_NUM_PIXINT;
  flag_srcsbin = DEF_FLAG_SRCSBIN;
  srcsbinsize = DEF_SRCSBINSIZE;
  flag_extnorm = DEF_FLAG_EXTNORM;
  nmax_srcs = DEF_NMAX_SRCS;
  nmax_fft = DEF_NMAX_FFT;
  flag_mcmcall = DEF_FLAG_MCMCALL;
  flag_addwcs = DEF_FLAG_ADDWCS;
  wcs_ra0 = DEF_WCS_RA0;
  wcs_dec0 = DEF_WCS_DEC0;
  dr_lens_center = DEF_DR_LENS_CENTER;
  flag_pow_tm15 = DEF_FLAG_POW_TM15;
  tol_pow_tm15 = DEF_TOL_POW_TM15;
    
  for(i=0;i<NMAX_LEN;i++){
    para_lens_min[i][0] = INIT_ZMIN;   /* zl */
    para_lens_max[i][0] = INIT_ZMAX; 
    para_lens_min[i][1] = INIT_MMIN;   /* mass */
    para_lens_max[i][1] = INIT_MMAX; 
    para_lens_min[i][2] = INIT_XYMIN;  /* x */
    para_lens_max[i][2] = INIT_XYMAX; 
    para_lens_min[i][3] = INIT_XYMIN;  /* y */
    para_lens_max[i][3] = INIT_XYMAX; 
    para_lens_min[i][4] = INIT_EMIN;   /* e */
    para_lens_max[i][4] = INIT_EMAX; 
    para_lens_min[i][5] = INIT_ANGMIN; /* theta_e */
    para_lens_max[i][5] = INIT_ANGMAX;   
    para_lens_min[i][6] = INIT_R0MIN;  /* r0 */
    para_lens_max[i][6] = INIT_R0MAX;
    para_lens_min[i][7] = INIT_NMIN;   /* n */
    para_lens_max[i][7] = INIT_NMAX;
    for(j=0;j<NPAR_LEN;j++){
      para_lens_med[i][j] = 0.0;
      para_lens_sig[i][j] = 0.0;
      para_lens_rai[i][j] = i;
      para_lens_raj[i][j] = j;
      para_lens_rat[i][j] = 1.0;
      para_lens_ras[i][j] = 0.0;
    }
  }

  for(i=0;i<NMAX_EXT;i++){
    para_ext_min[i][0] = INIT_ZMIN;   /* zs */
    para_ext_max[i][0] = INIT_ZMAX; 
    para_ext_min[i][1] = INIT_MMIN;   /* norm */
    para_ext_max[i][1] = INIT_MMAX;
    para_ext_min[i][2] = INIT_XYMIN;  /* x */
    para_ext_max[i][2] = INIT_XYMAX;
    para_ext_min[i][3] = INIT_XYMIN;  /* y */
    para_ext_max[i][3] = INIT_XYMAX;
    para_ext_min[i][4] = INIT_EMIN;   /* e */
    para_ext_max[i][4] = INIT_EMAX;
    para_ext_min[i][5] = INIT_ANGMIN; /* theta_e */
    para_ext_max[i][5] = INIT_ANGMAX;
    para_ext_min[i][6] = INIT_R0MIN;  /* r0 */
    para_ext_max[i][6] = INIT_R0MAX;
    para_ext_min[i][7] = INIT_NMIN;   /* n */
    para_ext_max[i][7] = INIT_NMAX;
    for(j=0;j<NPAR_EXT;j++){
      para_ext_med[i][j] = 0.0;
      para_ext_sig[i][j] = 0.0;
      para_ext_rai[i][j] = i;
      para_ext_raj[i][j] = j;
      para_ext_rat[i][j] = 1.0;
      para_ext_ras[i][j] = 0.0;
      para_extlen_rai[i][j] = -1;
      para_extlen_raj[i][j] = -1;
      para_extlen_rat[i][j] = 1.0;
      para_extlen_ras[i][j] = 0.0;
    }
    obs_ext_prior[i] = -1;
  }

  for(i=0;i<NMAX_POI;i++){
    para_poi_min[i][0] = INIT_ZMIN;  /* zs */
    para_poi_max[i][0] = INIT_ZMAX;  
    para_poi_min[i][1] = INIT_XYMIN; /* x */
    para_poi_max[i][1] = INIT_XYMAX; 
    para_poi_min[i][2] = INIT_XYMIN; /* y */
    para_poi_max[i][2] = INIT_XYMAX; 
    obs_numimg[i] = 0;
    for(j=0;j<NPAR_POITAB;j++){
      para_poi_med[i][j] = 0.0;
      para_poi_sig[i][j] = 0.0;
      para_poi_rai[i][j] = i;
      para_poi_raj[i][j] = j;
      para_poi_rat[i][j] = 1.0;
      para_poi_ras[i][j] = 0.0;
      para_poilen_rai[i][j] = -1;
      para_poilen_raj[i][j] = -1;
      para_poilen_rat[i][j] = 1.0;
      para_poilen_ras[i][j] = 0.0;
    }
    flux_poi[i] = 100.0;
  }

  para_psf_min[0] = INIT_R0MIN;   /* FWHM1 */
  para_psf_max[0] = INIT_R0MAX; 
  para_psf_min[1] = INIT_EMIN;    /* e1 */
  para_psf_max[1] = INIT_EMAX; 
  para_psf_min[2] = INIT_ANGMIN;  /* PA1 */
  para_psf_max[2] = INIT_ANGMAX; 
  para_psf_min[3] = INIT_BETAMIN; /* beta1 */
  para_psf_max[3] = INIT_BETAMAX; 
  para_psf_min[4] = INIT_R0MIN;  /* FWHM2 */
  para_psf_max[4] = INIT_R0MAX; 
  para_psf_min[5] = INIT_EMIN;   /* e2 */
  para_psf_max[5] = INIT_EMAX; 
  para_psf_min[6] = INIT_ANGMIN; /* PA2 */
  para_psf_max[6] = INIT_ANGMAX; 
  para_psf_min[7] = INIT_BETAMIN; /* beta2 */
  para_psf_max[7] = INIT_BETAMAX; 
  para_psf_min[8] = INIT_FMIN;    /* frac */
  para_psf_max[8] = INIT_FMAX; 

  para_psf[0] = DEF_SEEING;
  para_psf[1] = DEF_SEEING_E;
  para_psf[2] = DEF_SEEING_PA;
  para_psf[3] = DEF_SEEING_BETA;
  para_psf[4] = DEF_SEEING;
  para_psf[5] = DEF_SEEING_E;
  para_psf[6] = DEF_SEEING_PA;
  para_psf[7] = DEF_SEEING_BETA;
  para_psf[8] = 1.0;

  for(j=0;j<NPAR_PSF;j++){ 
    para_psf_med[j] = 0.0;
    para_psf_sig[j] = 0.0;
    para_psf_raj[j] = j;
    para_psf_rat[j] = 1.0;
    para_psf_ras[j] = 0.0;
  }

  omega_min = INIT_OMMIN;
  omega_max = INIT_OMMAX;
  lambda_min = INIT_OMMIN;
  lambda_max = INIT_OMMAX;
  weos_min = INIT_WMIN;
  weos_max = INIT_WMAX;
  hubble_min = INIT_HMIN;
  hubble_max = INIT_HMAX;

  for(i=0;i<NMAX_POI;i++){
    obs_numimg[i] = 0;
    for(j=0;j<NMAX_POIMG;j++){
      obs_parity[i][j] = 0;
      for(k=0;k<NPAR_READOBS;k++){
	tab_obs[i][j][k] = 0.0;
      }
    }
  } 
    
  return;
}

void init_flags(void)
{
  flag_set_array = 0;
  flag_set_point = 0;
  flag_arrayobs = 0;
  flag_obsmask = 0;
  flag_pointobs = 0;
  flag_computeall = 1;
  flag_set_srcs = 0;
  flag_obssig = 0;
  flag_seeing = 0;

  num_gal = 0;
  num_src = 0;
  num_mapprior = 0;
  num_lcent = 0;
  i_ext_fid = -1;

  return;
}

/*--------------------------------------------------------------
  dump current model
*/

void dump_model(FILE* fptr)
{
  int i, j;
  
  fprintf(fptr, "\n");
  for(i=0;i<num_len;i++){
    fprintf(fptr, "lens   %-7s %6.4f ", inttolmodel(model_lens[i]), para_lens[i][0]);
    for(j=1;j<NPAR_LEN;j++) fprintf(fptr, "%13e ", para_lens[i][j]);
    fprintf(fptr, "\n");
  }

  for(i=0;i<num_ext;i++){
    fprintf(fptr, "extend %-7s %6.4f ", inttoemodel(model_ext[i]), para_ext[i][0]);
    for(j=1;j<NPAR_EXT;j++) fprintf(fptr, "%13e ", para_ext[i][j]);
    fprintf(fptr, "\n");
  }

  for(i=0;i<num_poi;i++){
    fprintf(fptr, "point  %6.4f ", para_poi[i][0]);
    for(j=1;j<NPAR_POI;j++) fprintf(fptr, "%13e ", para_poi[i][j]);
    fprintf(fptr, "\n");
  }

  for(i=0;i<flag_seeing;i++){
    fprintf(fptr, "psf %13e ", para_psf[0]);
    for(j=1;j<NPAR_PSF;j++) fprintf(fptr, "%13e ", para_psf[j]);
    fprintf(fptr, "\n");
  }

  return;
}

void dump_opt_flag(FILE* fptr)
{
  int i, j;
  
  fprintf(fptr, "\n");
  for(i=0;i<num_len;i++){
    for(j=0;j<NPAR_LEN;j++) fprintf(fptr, "%d ", flag_para_lens[i][j]);
    fprintf(fptr, "\n");
  }

  for(i=0;i<num_ext;i++){
    for(j=0;j<NPAR_EXT;j++) fprintf(fptr, "%d ", flag_para_ext[i][j]);
    fprintf(fptr, "\n");
  }

  for(i=0;i<num_poi;i++){
    for(j=0;j<NPAR_POI;j++) fprintf(fptr, "%d ", flag_para_poi[i][j]);
    fprintf(fptr, "\n");
  }
  
  for(i=0;i<flag_seeing;i++){
    for(j=0;j<NPAR_PSF;j++) fprintf(fptr, "%d ", flag_para_psf[j]);
    fprintf(fptr, "\n");
  }

  return;
}

void dump_lensplane(FILE* fptr)
{
  int i;

  fprintf(fptr, "\n");
  for(i=0;i<num_lpl;i++){
    fprintf(fptr, "lens plane %2d : z = %f\n", i + 1, zl_lpl[i]);
  }

  return;
}

/*--------------------------------------------------------------
  output parameters
*/

void out_para(void)
{
  fprintf(stderr, "######## parameter definition\n");
  fprintf(stderr, "omega    =  %e\n", omega);
  fprintf(stderr, "lambda   =  %e\n", lambda);
  fprintf(stderr, "weos     =  %e\n", weos);
  fprintf(stderr, "hubble   =  %e\n", hubble);
  fprintf(stderr, "prefix   =  %s\n", file_prefix);
  fprintf(stderr, "xmin     =  %e\n", xmin);
  fprintf(stderr, "ymin     =  %e\n", ymin);
  fprintf(stderr, "xmax     =  %e\n", xmax);
  fprintf(stderr, "ymax     =  %e\n", ymax);
  fprintf(stderr, "pix_ext  =  %e\n", pix_ext);
  fprintf(stderr, "nx_ext   =  %d\n", nx_ext);
  fprintf(stderr, "ny_ext   =  %d\n", ny_ext);
  fprintf(stderr, "pix_poi  =  %e\n", pix_poi);
  fprintf(stderr, "nx_poi   =  %d\n", nx_poi);
  fprintf(stderr, "ny_poi   =  %d\n", ny_poi);
  fprintf(stderr, "maxlev   =  %d\n", maxlev);

  return;
}

/*--------------------------------------------------------------
  set parameter range
*/

void parprior(char *infile, int verb)
{
  int i, j, k, ii, jj, n, nn;
  double xx, yy, min, max, med, rat, sig;
  char buffer[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  char ptype[INPUT_MAXCHAR];
  FILE* fptr;

  fptr = fopen(infile, "r");

  if(verb == 1){
    fprintf(stderr, "######## reading parameter prior file\n");
    fprintf(stderr, " input file name = %s \n\n", infile);
  }
  
  if(fptr == NULL) terminator("failed at fopen (parprior)");
  
  n = 0;

  while(fgets(buffer, INPUT_MAXCHAR, fptr)){
    nn = sscanf(buffer, "%s %s", ptype, keyword);
    if(nn != EOF){
      if(ptype[0] != '#'){
	if(strcmp(ptype, "range") == 0){
	  if((strcmp(keyword, "lens") == 0) || (strcmp(keyword, "extend") == 0) || (strcmp(keyword, "point") == 0)){
	    nn = sscanf(buffer, "%s %s %d %d %lf %lf", ptype, keyword, &i, &j, &min, &max);
	    if(nn != 6) terminator("input file format irrelevant (parprior)"); 
	      
	    if(strcmp(keyword, "lens") == 0){
	      if((i > num_len) || (j > NPAR_LEN) || (i < 1) || (j < 1)){ terminator("lens id irrelevant (parprior)"); }
	      para_lens_min[i - 1][j - 1] = min;
	      para_lens_max[i - 1][j - 1] = max;
	      n++;
	    } else if(strcmp(keyword, "extend") == 0){
	      if((i > num_ext) || (j > NPAR_EXT) || (i < 1) || (j < 1)){ terminator("extend id irrelevant (parprior)"); }
	      para_ext_min[i - 1][j - 1] = min;
	      para_ext_max[i - 1][j - 1] = max;
	      n++;
	    } else if(strcmp(keyword, "point") == 0){
	      if((i > num_poi) || (j > NPAR_POI) || (i < 1) || (j < 1)){ terminator("point id irrelevant (parprior)"); }
	      para_poi_min[i - 1][j - 1] = min;
	      para_poi_max[i - 1][j - 1] = max;
	      n++;
	    } 
	  } else if((strcmp(keyword, "psf") == 0)){
	    nn = sscanf(buffer, "%s %s %d %lf %lf", ptype, keyword, &j, &min, &max);
	    if(nn != 5) terminator("input file format irrelevant (parprior)"); 
	    
	    if((j > NPAR_PSF) || (j < 1)){ terminator("psf id irrelevant (parprior)"); }
	    para_psf_min[j - 1] = min;
	    para_psf_max[j - 1] = max;
	    n++;
	  } else if((strcmp(keyword, "omega") == 0) || (strcmp(keyword, "lambda") == 0) || (strcmp(keyword, "weos") == 0) || (strcmp(keyword, "hubble") == 0)){
	    nn = sscanf(buffer, "%s %s %lf %lf", ptype, keyword, &min, &max);
	    if(nn != 4) terminator("input file format irrelevant (parprior)"); 
	    
	    if(strcmp(keyword, "omega") == 0){
	      omega_min = min;
	      omega_max = max;
	      n++;
	    } else if(strcmp(keyword, "lambda") == 0){
	      lambda_min = min;
	      lambda_max = max;
	      n++;
	    } else if(strcmp(keyword, "weos") == 0){
	      weos_min = min;
	      weos_max = max;
	      n++;
	    } else if(strcmp(keyword, "hubble") == 0){
	      hubble_min = min;
	      hubble_max = max;
	      n++;
	    } 
	  }
	}
	
	if(strcmp(ptype, "gauss") == 0){
	  if((strcmp(keyword, "lens") == 0) || (strcmp(keyword, "extend") == 0) || (strcmp(keyword, "point") == 0)){
	    
	    nn = sscanf(buffer, "%s %s %d %d %lf %lf", ptype, keyword, &i, &j, &med, &sig);
	    if(nn != 6) terminator("input file format irrelevant (parprior)"); 
	    
	    if(strcmp(keyword, "lens") == 0){
	      if((i > num_len) || (j > NPAR_LEN) || (i < 1) || (j < 1)){ terminator("lens id irrelevant (parprior)"); }
	      para_lens_med[i - 1][j - 1] = med;
	      para_lens_sig[i - 1][j - 1] = sig;
	      n++;
	    } else if(strcmp(keyword, "extend") == 0){
	      if((i > num_ext) || (j > NPAR_EXT) || (i < 1) || (j < 1)){ terminator("extend id irrelevant (parprior)"); }
	      para_ext_med[i - 1][j - 1] = med;
	      para_ext_sig[i - 1][j - 1] = sig;
	      n++;
	    } else if(strcmp(keyword, "point") == 0){
	      if((i > num_poi) || (j > NPAR_POI) || (i < 1) || (j < 1)){ terminator("point id irrelevant (parprior)"); }
	      para_poi_med[i - 1][j - 1] = med;
	      para_poi_sig[i - 1][j - 1] = sig;
	      n++;
	    } 	  
	  } else if((strcmp(keyword, "psf") == 0)){
	    nn = sscanf(buffer, "%s %s %d %lf %lf", ptype, keyword, &j, &med, &sig);
	    if(nn != 5) terminator("input file format irrelevant (parprior)"); 
	    
	    if((j > NPAR_PSF) || (j < 1)){ terminator("psf id irrelevant (parprior)"); }
	    para_psf_med[j - 1] = med;
	    para_psf_sig[j - 1] = sig;
	    n++;
	  } else if((strcmp(keyword, "omega") == 0) || (strcmp(keyword, "lambda") == 0) || (strcmp(keyword, "weos") == 0) || (strcmp(keyword, "hubble") == 0)){
	    nn = sscanf(buffer, "%s %s %lf %lf", ptype, keyword, &med, &sig);
	    if(nn != 4) terminator("input file format irrelevant (parprior)"); 
	    
	    if(strcmp(keyword, "omega") == 0){
	      omedian = med;
	      o_error = sig;
	      n++;
	    } else if(strcmp(keyword, "lambda") == 0){
	      lmedian = med;
	      l_error = sig;
	      n++;
	    } else if(strcmp(keyword, "weos") == 0){
	      wmedian = med;
	      w_error = sig;
	      n++;
	    } else if(strcmp(keyword, "hubble") == 0){
	      hmedian = med;
	      h_error = sig;
	      n++;
	    }
	  }
	}

	if(strcmp(ptype, "match") == 0){
	  nn = sscanf(buffer, "%s %s %d %d %d %d %lf %lf", ptype, keyword, &i, &j, &ii, &jj, &rat, &sig);
	  if(nn != 8) terminator("input file format irrelevant (parprior)"); 
	  
	  if(strcmp(keyword, "lens") == 0){
	    if((i > num_len) || (j > NPAR_LEN) || (i < 1) || (j < 1)){ terminator("lens id irrelevant (parprior)"); }
	    if((ii > num_len) || (jj > NPAR_LEN) || (ii < 1) || (jj < 1)){ terminator("lens id irrelevant (parprior)"); }
	    para_lens_rai[i - 1][j - 1] = ii - 1;
	    para_lens_raj[i - 1][j - 1] = jj - 1;
	    para_lens_rat[i - 1][j - 1] = rat;
	    para_lens_ras[i - 1][j - 1] = sig;
	    n++;
	  } else if(strcmp(keyword, "extend") == 0){
	    if((i > num_ext) || (j > NPAR_EXT) || (i < 1) || (j < 1)){ terminator("extend id irrelevant (parprior)"); }
	    if((ii > num_ext) || (jj > NPAR_EXT) || (ii < 1) || (jj < 1)){ terminator("extend id irrelevant (parprior)"); }
	    para_ext_rai[i - 1][j - 1] = ii - 1;
	    para_ext_raj[i - 1][j - 1] = jj - 1;
	    para_ext_rat[i - 1][j - 1] = rat;
	    para_ext_ras[i - 1][j - 1] = sig;
	    n++;
	  } else if(strcmp(keyword, "point") == 0){
	    if((i > num_poi) || (i < 1) || (j != 1)){ terminator("point id irrelevant (parprior)"); }
	    if((ii > num_poi) || (ii < 1) || (jj != 1)){ terminator("point id irrelevant (parprior)"); }
	    para_poi_rai[i - 1][j - 1] = ii - 1;
	    para_poi_raj[i - 1][j - 1] = jj - 1;
	    para_poi_rat[i - 1][j - 1] = rat;
	    para_poi_ras[i - 1][j - 1] = sig;
	    n++;
	  } else if(strcmp(keyword, "psf") == 0){
	    nn = sscanf(buffer, "%s %s %d %d %lf %lf", ptype, keyword, &j, &jj, &rat, &sig);
	    if(nn != 6) terminator("input file format irrelevant (parprior)"); 
	    if((j > NPAR_PSF) || (j < 1) || (jj > NPAR_PSF) || (jj < 1)){ terminator("psf id irrelevant (parprior)"); }
	    para_psf_raj[j - 1] = jj - 1;
	    para_psf_rat[j - 1] = rat;
	    para_psf_ras[j - 1] = sig;
	    n++;
	  } else if(strcmp(keyword, "extlens") == 0){
	    if((i > num_ext) || (j > NPAR_EXT) || (i < 1) || (j < 1)){ terminator("extend id irrelevant (parprior)"); }
	    if((ii > num_len) || (jj > NPAR_LEN) || (ii < 1) || (jj < 1)){ terminator("lens id irrelevant (parprior)"); }
	    para_extlen_rai[i - 1][j - 1] = ii - 1;
	    para_extlen_raj[i - 1][j - 1] = jj - 1;
	    para_extlen_rat[i - 1][j - 1] = rat;
	    para_extlen_ras[i - 1][j - 1] = sig;
	    n++;
	  } else if(strcmp(keyword, "poilens") == 0){
	    if((i > num_poi) || (j > NPAR_POI) || (i < 1) || (j < 1)){ terminator("point id irrelevant (parprior)"); }
	    if((ii > num_len) || (jj > NPAR_LEN) || (ii < 1) || (jj < 1)){ terminator("lens id irrelevant (parprior)"); }
	    para_poilen_rai[i - 1][j - 1] = ii - 1;
	    para_poilen_raj[i - 1][j - 1] = jj - 1;
	    para_poilen_rat[i - 1][j - 1] = rat;
	    para_poilen_ras[i - 1][j - 1] = sig;
	    n++;
	  } 
	} 
	
	if(strcmp(ptype, "obsext") == 0){
	  nn = sscanf(buffer, "%s %d %lf %lf", ptype, &i, &xx, &yy);
	  if(nn != 4) terminator("input file format irrelevant (parprior)"); 
	  if((i > NMAX_EXT) || (i < 1)){
	    terminator("extobs id out of range (parprior)"); 
	  }
	  k = xytok_ext(xx, yy);
	  if((k < 0) || (k >= (nx_ext * ny_ext))){
	    terminator("obsext position out of field (parprior)");
	  }
	  obs_ext_prior[i - 1] = k;
	  n++;
	}
	
      }
    }
  }

  fclose(fptr);

  if(verb == 1){
    fprintf(stderr, "read %d priors\n\n", n);
  }
  
  return;
}

/*--------------------------------------------------------------
  read galaxy file
*/

void mapprior(char *infile, int verb)
{
  int n, nn;
  double x, y, z, p, e;
  char buffer[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  FILE* fptr;

  fptr = fopen(infile, "r");

  if(verb == 1){
    fprintf(stderr, "######## reading map prior file\n");
    fprintf(stderr, " input file name = %s \n\n", infile);
  }
  
  if(fptr == NULL) terminator("failed at fopen (mapprior)");
  
  n = 0;
  while(fgets(buffer, INPUT_MAXCHAR, fptr)){
    if(sscanf(buffer, "%s", keyword) != EOF){
      if(keyword[0] != '#'){
	nn = sscanf(buffer, "%s %lf %lf %lf %lf %lf", keyword, &z, &x, &y, &p, &e);
	if((nn != 6) || (fabs(e) < TOL_ERROR_MAPPRIOR)){
	  terminator("input file format irrelevant (mapprior)"); 
	}
	if(n >= NMAX_MAPPRIOR) terminator("mapprior number exceeds maximum"); 
	flag_para_mapprior[n] = 0;
	if(strcmp(keyword, "kappa") == 0){ flag_para_mapprior[n] = 1; }
	else if(strcmp(keyword, "mag") == 0){ flag_para_mapprior[n] = 2; }
	else if(strcmp(keyword, "gamma1") == 0){ flag_para_mapprior[n] = 3; }
	else if(strcmp(keyword, "gamma2") == 0){ flag_para_mapprior[n] = 4; }
	else if(strcmp(keyword, "g1") == 0){ flag_para_mapprior[n] = 5; }
	else if(strcmp(keyword, "g2") == 0){ flag_para_mapprior[n] = 6; }
	if(flag_para_mapprior[n] == 0){ 
	  terminator("input file format irrelevant (mapprior)"); 
	}
	para_mapprior[n][0] = z;
	para_mapprior[n][1] = x;
	para_mapprior[n][2] = y;
	para_mapprior[n][3] = p;
	para_mapprior[n][4] = e;
	n++;
      }
    }
  }
  num_mapprior = n;

  fclose(fptr);

  if(verb == 1){
    fprintf(stderr, "read %d priors\n\n", n);
  }
  
  return;
}

/*--------------------------------------------------------------
  read galaxy file
*/

void readgals(void)
{
  int i, nn;
  double x, y, l, e, t;
  char buffer[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  FILE* fptr;

  fptr = fopen(file_gal, "r");

  if(fptr == NULL) terminator("failed at fopen (readgals)"); 

  i = 0;
  while(fgets(buffer, INPUT_MAXCHAR, fptr)){
    e = 0.0;
    t = 0.0;
    if(sscanf(buffer, "%s", keyword) != EOF){
      if(keyword[0] != '#'){
	nn = sscanf(buffer, "%lf %lf %lf %lf %lf", &x, &y, &l, &e, &t);
	if((nn != 5) && (nn != 3)){
	  terminator("input file format irrelevant (readgals)"); 
	}
	if(i >= NMAX_GALS) terminator("galaxy number exceeds maximum"); 
	if((l < 0.0) || (e < INIT_EMIN) || (e >= INIT_EMAX)){ 
	  terminator("input file format irrelevant (readgals)"); 
	}
	para_gals[i][0] = x;
	para_gals[i][1] = y;
	para_gals[i][2] = l;
	para_gals[i][3] = e;
	para_gals[i][4] = t;
	i++;
      }
    }
  }
  num_gal = i;

  fclose(fptr);

  return;
}

void readsrcs(void)
{
  int i, nn;
  double f, x, y, e, t, r, n;
  char buffer[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  FILE* fptr;

  unset_srcs();
  unset_tab_calc_src();

  /* check num_src */

  fptr = fopen(file_src, "r");
  if(fptr == NULL) terminator("failed at fopen (readsrcs)"); 

  i = 0;
  while(fgets(buffer, INPUT_MAXCHAR, fptr)){
    if(sscanf(buffer, "%s", keyword) != EOF){
      if(keyword[0] != '#'){
	if(i >= nmax_srcs) terminator("source number exceeds maximum"); 
	i++;
      }
    }
  }
  num_src = i;

  fclose(fptr);

  /* memory allocation */

  para_srcs = (float*)malloc(sizeof(float) * (num_src * NPAR_SRC));
  if(para_srcs == NULL) terminator("memory allocation failed");

  /* read values */

  fptr = fopen(file_src, "r");
  if(fptr == NULL) terminator("failed at fopen (readsrcs)"); 

  i = 0;
  while(fgets(buffer, INPUT_MAXCHAR, fptr)){
    e = 0.0;
    t = 0.0;
    if(sscanf(buffer, "%s", keyword) != EOF){
      if(keyword[0] != '#'){
	nn = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf", &f, &x, &y, &e, &t, &r, &n);
	if((nn != NPAR_SRC)){
	  terminator("input file format irrelevant (readsrcs)"); 
	}
	if((f <= 0.0) || (e < INIT_EMIN) || (e >= INIT_EMAX)){ 
	  terminator("input file format irrelevant (readsrcs)"); 
	}
	para_srcs[i * NPAR_SRC] = f;
	para_srcs[1 + i * NPAR_SRC] = x;
	para_srcs[2 + i * NPAR_SRC] = y;
	para_srcs[3 + i * NPAR_SRC] = e;
	para_srcs[4 + i * NPAR_SRC] = t;
	para_srcs[5 + i * NPAR_SRC] = r;
	para_srcs[6 + i * NPAR_SRC] = n;
	i++;
      }
    }
  }

  fclose(fptr);

  flag_set_srcs = 1;

  return;
}

void unset_srcs(void)
{
  if(flag_set_srcs == 1){
    free(para_srcs);
    flag_set_srcs = 0;
  }

  return;
}

/*--------------------------------------------------------------
  read obs file (point)
*/

void readobs_point(char *infile, int verb)
{
  int i, j, k, ii, jj, jjj, nn;
  double zs, zserr;
  char buffer[INPUT_MAXCHAR];
  FILE* fptr;

  if(verb == 1){
    fprintf(stderr, "######## reading obs file for point\n");
    fprintf(stderr, " input file name = %s \n\n", infile);
  }
  
  fptr = fopen(infile, "r");

  if(fptr == NULL) terminator("failed at fopen (readobs_point)");
  
  flag_pointobs = 1;

  ii = 0;
  jjj = 0;
  while(fgets(buffer, INPUT_MAXCHAR, fptr)){
    i = 0;
    j = 0;
    zs = -1.0;
    zserr = 0.0;
    if(buffer[0] != '#'){
      nn = sscanf(buffer, "%d %d %lf %lf", &i, &j, &zs, &zserr);
      if(nn != EOF){
	if((i < 1) || (i > num_poi)){
	  terminator("id out of range (readobs_point)"); 
	}
	if(j > NMAX_POIMG) terminator("too many images (readobs_point)"); 
	
	obs_numimg[i - 1] = j;
	
	if(zs >= 0.0){
	  if((fabs(zs - para_poi[i - 1][0]) > TOL_ZS) && (zserr == 0.0)){
	    terminator("wrong zs (readobs_point)"); 
	  }
	} else {
	  zserr = 0.0;
	}

	para_poi_med[i - 1][0] = zs;
	para_poi_sig[i - 1][0] = zserr;

	ii++;
	jj = 0;
	for(k=0;k<j;k++){
	  if(fgets(buffer, INPUT_MAXCHAR, fptr)){
	      nn = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf %d", &tab_obs[i - 1][k][0], &tab_obs[i - 1][k][1], &tab_obs[i - 1][k][2], &tab_obs[i - 1][k][3], &tab_obs[i - 1][k][4], &tab_obs[i - 1][k][5], &tab_obs[i - 1][k][6], &obs_parity[i - 1][k]);
	      if(nn != 8) terminator("wrong format (readobs_point)"); 
	      jj++;
	  }
	}
	if(j != jj) terminator("wrong format (readobs_point)"); 
	jjj = jjj + jj;
      }
    }
  }

  if(verb == 1){
    fprintf(stderr, "read %d sources, %d images in total\n\n", ii, jjj);
  }
 
  fclose(fptr);

  return;
}

void reset_obs_point(int i, int j, int k, double p)
{
  if((i > num_poi) || (i <= 0)  || (j > obs_numimg[i - 1]) || (j <= 0) || (k > (NPAR_READOBS + 1)) || (k <= 0))
    terminator("id irrelevant (reset_obs_point)"); 

  if(k == (NPAR_READOBS + 1)){
    obs_parity[i - 1][j - 1] = (int)p;
  } else {
    tab_obs[i - 1][j - 1][k - 1] = p;
  }
  
  return;
}

#undef PAREADXE
#undef PAREADXD
#undef PAREADXS
