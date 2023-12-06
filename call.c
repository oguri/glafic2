#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

static float *array_call;
static int flag_array_call;

/*--------------------------------------------------------------
  initialization, quit
*/

void glafic_init(double in_omega, double in_lambda, double in_weos, double in_hubble, char *in_file_prefix, double in_xmin,  double in_ymin, double in_xmax,double in_ymax, double in_pix_ext, double in_pix_poi, int in_maxlev, int in_ran_seed, int verb)
{
  init_flags();
  def_parameters();
 
  glafic_set_primary(in_omega, in_lambda, in_weos, in_hubble, in_file_prefix, in_xmin, in_ymin, in_xmax, in_ymax, in_pix_ext, in_pix_poi, in_maxlev, verb);

  /* use default seed when in_ran_seed == 0 */
  if(in_ran_seed != 0){
    ran_seed = in_ran_seed;
  }

  if(verb == 1) fprintf(stderr, "ran_seed =  %d\n", ran_seed);
  
  ran_gsl = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(ran_gsl, ran_seed);

  return;
}

void glafic_init_file(char *infile, int verb)
{
  init_flags();

  /* read parameters */
  init_para(infile, verb);
  init_para2(infile, verb);
  /* startup */
  startup(infile, verb);
  /* set para for opt */
  setopt(infile, verb);
  /* prepare lens planes */
  gen_lensplane(verb);
  
  ran_gsl = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(ran_gsl, ran_seed);

  set_distance_lpl_init();

  return;
}

void glafic_set_primary(double in_omega, double in_lambda, double in_weos, double in_hubble, char *in_file_prefix, double in_xmin, double in_ymin, double in_xmax, double in_ymax, double in_pix_ext, double in_pix_poi, int in_maxlev, int verb)
{
  sprintf(file_prefix, "%s", in_file_prefix);
  xmin = in_xmin;
  ymin = in_ymin;
  xmax = in_xmax;
  ymax = in_ymax;
  pix_ext = in_pix_ext;
  pix_poi = in_pix_poi;
  maxlev = in_maxlev;
  
  glafic_set_cosmo(in_omega, in_lambda, in_weos, in_hubble);

  set_npix();
  glafic_unset_array_extend();
  if(verb == 1) out_para();
  
  return;
}

void glafic_set_cosmo(double in_omega, double in_lambda, double in_weos, double in_hubble)
{
  omega = in_omega;
  lambda = in_lambda;
  weos = in_weos;
  hubble = in_hubble;

  ext_unset_table();
  poi_unset_table();
  poimg_unset_table();
  
  return;
}

void glafic_quit(void)
{
  ext_unset_table();
  obs_unset_table();
  poi_unset_table();
  poimg_unset_table();
  unset_srcs();
  unset_tab_calc_src();
  glafic_unset_array_extend();
  gsl_rng_free(ran_gsl);

  return;
}

void glafic_set_secondary(char *buffer, int verb)
{
  char keyword[INPUT_MAXCHAR];

  if(sscanf(buffer, "%s", keyword) != EOF){
    if(keyword[0] != '#'){
      init_para2_body(keyword, buffer, verb);
    }
  }

  return;
}

/*--------------------------------------------------------------
  setting lens model
*/

void glafic_startup_setnum(int in_num_len, int in_num_ext, int in_num_poi)
{
  num_len = in_num_len;
  num_ext = in_num_ext;
  num_poi = in_num_poi;

  return;
}

void glafic_set_lens(int id, char *model, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8)
{
  if((id <= 0) || (id > num_len))
    terminator("id out of range (glafic_set_lens)");

  model_lens[id - 1] = lmodeltoint(model);
  if((strcmp(model, "gals") == 0) && (num_gal == 0)) readgals(); 
  
  para_lens[id - 1][0] = p1;
  para_lens[id - 1][1] = p2;
  para_lens[id - 1][2] = p3;
  para_lens[id - 1][3] = p4;
  para_lens[id - 1][4] = p5;
  para_lens[id - 1][5] = p6;
  para_lens[id - 1][6] = p7;
  para_lens[id - 1][7] = p8;

  return;
}

void glafic_set_extend(int id, char *model, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8)
{
  if((id <= 0) || (id > num_ext))
    terminator("id out of range (glafic_set_extend)");

  model_ext[id - 1] = emodeltoint(model);
  if((strcmp(model, "srcs") == 0) && (num_src == 0)) readsrcs(); 
  
  para_ext[id - 1][0] = p1;
  para_ext[id - 1][1] = p2;
  para_ext[id - 1][2] = p3;
  para_ext[id - 1][3] = p4;
  para_ext[id - 1][4] = p5;
  para_ext[id - 1][5] = p6;
  para_ext[id - 1][6] = p7;
  para_ext[id - 1][7] = p8;
  
  return;
}

void glafic_set_point(int id, double p1, double p2, double p3)
{
  if((id <= 0) || (id > num_poi))
    terminator("id out of range (glafic_set_point)");

  para_poi[id - 1][0] = p1;
  para_poi[id - 1][1] = p2;
  para_poi[id - 1][2] = p3;
  
  return;
}

void glafic_set_psf(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9)
{
  para_psf[0] = p1;
  para_psf[1] = p2;
  para_psf[2] = p3;
  para_psf[3] = p4;
  para_psf[4] = p5;
  para_psf[5] = p6;
  para_psf[6] = p7;
  para_psf[7] = p8;
  para_psf[8] = p9;

  flag_seeing = 1;
  
  return;
}

void glafic_setopt_lens(int id, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8)
{
  if((id <= 0) || (id > num_len))
    terminator("id out of range (glafic_setopt_lens)");

  flag_para_lens[id - 1][0] = p1;
  flag_para_lens[id - 1][1] = p2;
  flag_para_lens[id - 1][2] = p3;
  flag_para_lens[id - 1][3] = p4;
  flag_para_lens[id - 1][4] = p5;
  flag_para_lens[id - 1][5] = p6;
  flag_para_lens[id - 1][6] = p7;
  flag_para_lens[id - 1][7] = p8;
  
  return;
}

void glafic_setopt_extend(int id, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8)
{
  if((id <= 0) || (id > num_ext))
    terminator("id out of range (glafic_setopt_extend)");

  flag_para_ext[id - 1][0] = p1;
  flag_para_ext[id - 1][1] = p2;
  flag_para_ext[id - 1][2] = p3;
  flag_para_ext[id - 1][3] = p4;
  flag_para_ext[id - 1][4] = p5;
  flag_para_ext[id - 1][5] = p6;
  flag_para_ext[id - 1][6] = p7;
  flag_para_ext[id - 1][7] = p8;
  
  return;
}

void glafic_setopt_point(int id, int p1, int p2, int p3)
{
  if((id <= 0) || (id > num_poi))
    terminator("id out of range (glafic_setopt_point)");

  flag_para_poi[id - 1][0] = p1;
  flag_para_poi[id - 1][1] = p2;
  flag_para_poi[id - 1][2] = p3;
  
  return;
}

void glafic_setopt_psf(int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8, int p9)
{
  flag_para_psf[0] = p1;
  flag_para_psf[1] = p2;
  flag_para_psf[2] = p3;
  flag_para_psf[3] = p4;
  flag_para_psf[4] = p5;
  flag_para_psf[5] = p6;
  flag_para_psf[6] = p7;
  flag_para_psf[7] = p8;
  flag_para_psf[8] = p9;
  
  return;
}

void glafic_model_init(int verb)
{
  poi_unset_table();
  ext_unset_table();
  poimg_unset_table();

  if(verb == 1){
    dump_model(stderr);
    dump_opt_flag(stderr);
  }
  gen_lensplane(verb);
 
  set_distance_lpl_init();

  if(verb == 1) fprintf(stderr, "\n");

  return;
}

/*--------------------------------------------------------------
  lens properties
*/

void glafic_calcimage(double zs, double x, double y, double pout[NPAR_LMODEL], int alponly, int verb)
{
  set_distance_lpl_zs(zs);
  lensmodel(x, y, pout, alponly, 0);

  if(verb == 1){
    printf("alpha_x  = %13e\n", pout[0]);
    printf("alpha_y  = %13e\n", pout[1]);
    printf("td[day]  = %13e\n", pout[2]);
    printf("kappa    = %13e\n", pout[3]);
    printf("gamma_1  = %13e\n", pout[4]);
    printf("gamma_2  = %13e\n", pout[5]);
    printf("mag^-1   = %13e\n", pout[6]);
    printf("rotation = %13e\n", pout[7]);
    fprintf(stderr, "\n");
  }
  
  return;
}

double glafic_calcein_i(int id, double zs)
{
  double ein;
  
  set_distance_lpl_zs(zs);

  ein = calcein_i_calc(id - 1, zs);

  return ein;
}

double glafic_calcein2(int id, double zs, double x0, double y0)
{
  double ein;
  
  set_distance_lpl_zs(zs);

  ein = calcein2_calc(zs, x0, y0, id);

  return ein;
}

double glafic_kappa_ave(int id, double zs, double r, double x0, double y0)
{
  double kap;
  
  set_distance_lpl_zs(zs);

  kap = calc_kappa_ave(r, x0, y0, id);

  return kap;
}

double glafic_kappa_cum(int id, double zs, double r, double x0, double y0)
{
  double kap;
  
  set_distance_lpl_zs(zs);

  kap = calc_kappa_cum(r, x0, y0, id);

  return kap;
}

/*--------------------------------------------------------------
  solve lens equation for point source
*/

/* rr[i][0]: x position of i-th image
   rr[i][1]: y position of i-th image
   rr[i][2]: (signed) magnification of i-th image
   rr[i][3]: time delay of i-th image
 */
void glafic_point_solve(double zs, double x, double y, int *ni, double rr[NMAX_POIMG][NPAR_IMAGE], int verb)
{
  findimg(x, y, zs, ni, rr, verb);

  return;
}

void glafic_findimg_i(int id, int *ni, double rr[NMAX_POIMG][NPAR_IMAGE], int verb)
{
  if((id <= 0) || (id > num_poi))
    terminator("id out of range (glafic_findimg_i)");

  glafic_point_solve(para_poi[id - 1][0], para_poi[id - 1][1], para_poi[id - 1][2], ni, rr, verb);
  
  return;
}

void glafic_findimg(void)
{
  findimg_i(0);

  return;
}

void glafic_writelens(double zs)
{
  writelens(zs);

  return;
}

void glafic_writecrit(double zs)
{
  double csize[NPAR_CAUSIZE];

  writecrit(zs, csize, 1, 1);

  return;
}

void glafic_writemesh(double zs)
{
  writemesh(zs);

  return;
}

void glafic_lenscenter(double zs)
{
  lenscenter(zs);

  return;
}

/*--------------------------------------------------------------
  solve lens equation for extended source
*/

void glafic_set_array_extend(int id, double sky, double noise, int flag_source)
{
  int l, k, nn;
  double f, p;

  if(flag_array_call != TFLAG_VALUE){
    flag_array_call = TFLAG_VALUE;
    array_call = (float*)malloc(sizeof(float) * nx_ext * ny_ext);
    if(array_call == NULL) terminator("memory allocation failed");
  }
  
  flag_computeall = 1;
  i_ext_fid = -1;
  ext_set_image(id, flag_source, 1);

  nn = nx_ext * ny_ext;

  for(k=0;k<nn;k++){
    if(id == 0){
      f = 0.0;
      for(l=0;l<num_ext;l++){
	f = f + array_ext_img[k + l * nn];
      }
    } else {
      f = array_ext_img[k + (id - 1) * nn];
    }
    p = calc_pix_noise(f, sky, noise);
    array_call[k] = (float)p;
  }
  
  return;
}

void glafic_unset_array_extend(void)
{
  if(flag_array_call == TFLAG_VALUE){
    flag_array_call = 0;
    free(array_call);
  }

  return;
}

void glafic_readpsf(char *fname, int verb)
{
  read_psffits(fname, verb);
  
  return;
}

double glafic_extend_array_ij(int i, int j)
{
  return array_call[i + j * nx_ext];
}

double glafic_extend_array_k(int k)
{
  return array_call[k];
}

void glafic_extend_array_ktoxy(int k, double *x, double *y)
{
  ktoxy_ext(k, x, y);

  return;
}

void glafic_writepsf(void)
{
  writepsf();

  return;
}

/*------------------------------------
  reading (writing) opt-related data
*/

void glafic_readobs_extend(char *fname, char *fname_mask, int verb)
{
  readobs_extend(fname, verb);
  
  if(strcmp(fname_mask, "N") != 0){
    readmask(fname_mask, verb);
  }
  
  calc_obsnoise();

  return;
}

void glafic_readnoise_extend(char *fname, int verb)
{
  readnoise_extend(fname, verb);

  return;
}
  
void glafic_readobs_point(char *fname, int verb)
{
  readobs_point(fname, verb);

  return;
}
  
void glafic_parprior(char *fname, int verb)
{
  parprior(fname, verb);

  return;
}
  
void glafic_mapprior(char *fname, int verb)
{
  mapprior(fname, verb);

  return;
}
  
/*------------------------------------
  optimizations
*/

void glafic_optimize(int verb)
{
  opt_lens(0, verb);

  return;
}

void glafic_optpoint(int verb)
{
  double c2min[NMAX_POI][NPAR_CHI2];

  chi2calc_opt_point(c2min, verb);

  return;
}

void glafic_optextend(int verb)
{
  double chi2min[NPAR_CHI2MIN];
  
  chi2calc_opt_extend(chi2min, verb, 0);

  return;
}

double glafic_c2calc(void)
{
  opt_lens_static(-1);
  
  return chi2calc_nopar();
}

void glafic_reset_obs_point(int i, int j, int k, double p)
{
  reset_obs_point(i, j, k, p);

  return;
}

/*------------------------------------
  get parameters
*/

double glafic_getpar_lens(int id, int ip)
{
  if((id <= 0) || (id > num_len) || (ip <= 0) || (ip > NPAR_LEN)) 
    terminator("id out of range (glafic_getpar_lens)");

  return para_lens[id - 1][ip - 1];
}

double glafic_getpar_extend(int id, int ip)
{
  if((id <= 0) || (id > num_ext) || (ip <= 0) || (ip > NPAR_EXT))
    terminator("id out of range (glafic_getpar_extend)");

  return para_ext[id - 1][ip - 1];
}

double glafic_getpar_point(int id, int ip)
{
  if((id <= 0) || (id > num_poi) || (ip <= 0) || (ip > NPAR_POI))
    terminator("id out of range (glafic_getpar_point)");

  return para_poi[id - 1][ip - 1];
}

double glafic_getpar_psf(int ip)
{

  if((flag_seeing == 0) || (ip <= 0) || (ip > NPAR_PSF))
    terminator("id out of range (glafic_getpar_psf)");

  return para_psf[ip - 1];
}

double glafic_getpar_omega(void)
{
  return omega;
}

double glafic_getpar_lambda(void)
{
  return lambda;
}

double glafic_getpar_hubble(void)
{
  return hubble;
}

double glafic_getpar_weos(void)
{
  return weos;
}

double glafic_getpar_sky(void)
{
  return skymed;
}
