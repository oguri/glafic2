#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "glafic.h"

/*--------------------------------------------------------------
  functions using cosmo_para are defined here
*/

typedef struct {
  double o;
  double l;
  double w;
} cosmo_para;

double dis_chi(double a, double b, cosmo_para *cpar);
double dis_comoving(double a, double b, cosmo_para *cpar);
double dis_angulard_cpar(double a, double b, cosmo_para *cpar);
double dis_luminosity_cpar(double a, double b, cosmo_para *cpar);
double derivative_dis_proper_cpar(double z, cosmo_para *cpar);
double inv_hubble_aint(double a, void *p);
double hubble_ez(double z, cosmo_para *cpar);
double hubble_ez2(double z, cosmo_para *cpar);
void set_cosmo_para(cosmo_para *cpar);

/*--------------------------------------------------------------
  radial distance
*/

double dis_chi(double a, double b, cosmo_para *cpar)
{
  int gslstatus;
  double integral;
  gsl_integration_cquad_workspace *workspace = NULL;
  /* size_t neval;
     cgsl_integration_romberg_workspace *workspace = NULL; */
  gsl_function F;

  F.function = &inv_hubble_aint;
  F.params = cpar;

  workspace = gsl_integration_cquad_workspace_alloc( GSL_CQUAD_ITERATION );
  /* workspace = gsl_integration_romberg_alloc( GSL_ROMBERG_N_DISTANCE ); */
  
  if(workspace == NULL){
    fprintf(stderr, "failed to create workspace in dis_chi\n");
    exit(EXIT_FAILURE);
  } else {
    gslstatus = gsl_integration_cquad(&F, 1.0 / (1.0 + b), 1.0 / (1.0 + a), 0.0, GSL_EPSREL_DISTANCE, workspace, &integral, NULL, NULL);
    /* gslstatus = gsl_integration_romberg(&F, 1.0 / (1.0 + b), 1.0 / (1.0 + a), 0.0, GSL_EPSREL_DISTANCE, &integral, &neval, workspace); */
    if (gslstatus != GSL_SUCCESS){
      fprintf(stderr, "integration failed in dis_chi\n");
      exit(EXIT_FAILURE);
    }
  }

  gsl_integration_cquad_workspace_free(workspace);
  /* gsl_integration_romberg_free(workspace); */

  return integral;
}

/*--------------------------------------------------------------
  comoving distance
*/

double dis_comoving(double a, double b, cosmo_para *cpar)
{
  double k, chi;

  if(a >= b){ return 0.0; }
  
  chi = dis_chi(a, b, cpar);
  
  k = cpar->o + cpar->l - 1.0;

  if(fabs(k) < TOL_CURVATURE){
    return chi;
  } else if(k > 0.0){
    return sin(chi * sqrt(k)) / sqrt(k);
  } else {
    return sinh(chi * sqrt(-k)) / sqrt(-k);
  }
}

/*--------------------------------------------------------------
  angular diameter distance
*/

double dis_angulard(double a, double b)
{
  cosmo_para cpar;

  set_cosmo_para(&cpar);

  return dis_angulard_cpar(a, b, &cpar);
}

double dis_angulard_cpar(double a, double b, cosmo_para *cpar)
{
  return dis_comoving(a, b, cpar) / (1.0 + b);
}

/*--------------------------------------------------------------
  luminosity distance
*/

double dis_luminosity(double a, double b)
{
  cosmo_para cpar;

  set_cosmo_para(&cpar);

  return dis_luminosity_cpar(a, b, &cpar);
}

double dis_luminosity_cpar(double a, double b, cosmo_para *cpar)
{
  return dis_comoving(a, b, cpar) * (1.0 + b) / (1.0 + a);
}

/*--------------------------------------------------------------
  z derivative of proper distance
*/

double derivative_dis_proper(double z)
{
  cosmo_para cpar;

  set_cosmo_para(&cpar);

  return derivative_dis_proper_cpar(z, &cpar);
}

double derivative_dis_proper_cpar(double z, cosmo_para *cpar)
{
  return (1.0 / hubble_ez(z, cpar)) / (1.0 + z);
}

/*--------------------------------------------------------------
  distance moduls w/o hubble
*/

double dis_mod(double z)
{
  double dl;

  dl = dis_luminosity(0.0, z);

  return 25.0 + 5.0 * log10(COVERH_MPCH * dl);
}

/*--------------------------------------------------------------
  Hubble parameter
*/

double inv_hubble_aint(double a, void *p)
{
  cosmo_para *cpar = (cosmo_para *)p;
  
  return  1.0 / (a * a * hubble_ez(1.0 / a - 1.0, cpar));
}

double hubble_ez(double z, cosmo_para *cpar)
{
  double h2;
  
  h2 = hubble_ez2(z, cpar);
  
  return sqrt(h2);
}

double hubble_ez2(double z, cosmo_para *cpar)
{
  double o, l, w, h2;

  o = cpar->o;
  l = cpar->l;
  w = cpar->w;
  
  h2 = (1.0 + o * z - l) * (1.0 + z) * (1.0 + z) + l * pow(1.0 + z, 3.0 * (1.0 + w));
  
  return h2;
}

/*--------------------------------------------------------------
  set distance
*/

void set_cosmo_para(cosmo_para *cpar)
{
  cpar->o = omega;
  cpar->l = lambda; 
  cpar->w = weos;   
}

void set_distance_lpl_zs(double zs)
{
  int i, j;

  nlp = 0;
  for(i=0;i<num_lpl;i++){
    if(zl_lpl[i] < zs) nlp = i + 1;
  }

  dis_os = dis_angulard(0.0, zs);

  for(i=0;i<nlp;i++){
    dis_ls_lpl[i] = dis_angulard(zl_lpl[i], zs);
  }
    
  for(j=1;j<nlp;j++){
    for(i=0;i<j;i++){
      dis_beta[i][j] = dis_beta_fac[i][j] * dis_os / dis_ls_lpl[i];
      dis_tdelay[i][j] = dis_tdelay_fac[i][j];
    }
  }
  for(i=0;i<nlp;i++){
    dis_beta[i][nlp] = 1.0;
    dis_tdelay[i][nlp] = tdelay_fac(zl_lpl[nlp - 1], dis_os, dis_ol_lpl[nlp - 1], dis_ls_lpl[nlp - 1]);
  }
		
  return;
}

void set_distance_lpl_i(int i)
{
  delome = delome_lpl[i];
  dis_ol = dis_ol_lpl[i];
  dis_ls = dis_ls_lpl[i];
  zl_ext = zl_lpl[i];
}

void set_distance_lpl_init(void)
{
  int i, j;
  double dij, dj;

  gen_lensplane(0);

  for(i=0;i<num_lpl;i++){
    dis_ol_lpl[i] = dis_angulard(0.0, zl_lpl[i]);
    delome_lpl[i] = deltaomega(zl_lpl[i]);
    for(j=0;j<num_lpl;j++) dis_beta_fac[i][j] = 0.0;
  }

  for(j=1;j<num_lpl;j++){
    for(i=0;i<j;i++){
      dij = dis_angulard(zl_lpl[i], zl_lpl[j]);
      dj  = dis_angulard(0.0, zl_lpl[j]);
      dis_beta_fac[i][j] = dij / dj;
      dis_tdelay_fac[i][j] = tdelay_fac(zl_lpl[i], dj, dis_ol_lpl[i], dij);
    }
  }

  return;
}

void set_distance_raw_zlzs(double zl, double zs)
{
  set_distance_raw_zl(zl);
  
  dis_os = dis_angulard(0.0, zs);
  dis_ls = dis_angulard(zl, zs);
  
  return;
}

void set_distance_raw_zl(double zl)
{
  if(zl <= 0.0) terminator("irrelevant lens redshift");

  delome = deltaomega(zl);
  dis_ol = dis_angulard(0.0, zl);
  zl_ext = zl;
  
  return;
}

void set_distance_facext(void)
{
  int i;
  double zsmax = INIT_ZMIN;
  
  if(num_lpl == 1){
    i_ext_fid = 0;
    for(i=0;i<num_ext;i++){
      if(para_ext[i][0] > zsmax){
	zsmax = para_ext[i][0];
	i_ext_fid = i;
      }
    }
    for(i=0;i<num_ext;i++){
      dis_fac_ext[i] = disratio(zl_lpl[0], para_ext[i_ext_fid][0], para_ext[i][0]);
    }
  } else {
    i_ext_fid = -1;
    for(i=0;i<num_ext;i++){
      dis_fac_ext[i] = 1.0;
    }
  }

  return;
}

double disratio(double zl, double zs_fid, double zs)
{
  if(zs_fid <= zl){
    /* terminator("irrelevant source redshift"); */
    return 0.0;
  } else {
    double d1, d2;
    
    d1 = dis_angulard(zl, zs_fid) / dis_angulard(0.0, zs_fid);
    d2 = dis_angulard(zl, zs) / dis_angulard(0.0, zs);
  
    return d2 / d1;
  }
} 

/*--------------------------------------------------------------
  for time delay [day]
*/

double tdelay_fac(double zl, double dos, double dol, double dls)
{
  double ddd;

  ddd = (1.0 + zl) * dos * dol / (dls + OFFSET_TDELAY_FAC);
 
  return FAC_TDELAY_DAY * ARCSEC2RADIAN * ARCSEC2RADIAN * ddd / hubble;
}

/*--------------------------------------------------------------
  arcsec <-> Mpc/h conversion
*/

double thetator(double theta)
{
  return thetator_dis(theta, dis_ol);
}

double thetator_dis(double theta, double dis)
{
  return COVERH_MPCH * dis * ARCSEC2RADIAN * theta;
}

double rtotheta(double r)
{
  return rtotheta_dis(r, dis_ol);
}

double rtotheta_dis(double r, double dis)
{
  return r / (COVERH_MPCH * dis * ARCSEC2RADIAN);
}

/*--------------------------------------------------------------
  critical Sigma in units of h * M_sun Mpc^-2 (physical)
*/

double inv_sigma_crit(void)
{
  return inv_sigma_crit_dis(dis_os, dis_ol, dis_ls);
}

double inv_sigma_crit_dis(double dos, double dol, double dls)
{
  return dol * dls / (FAC_CRITDENS * dos);
}

double sigma_crit(void)
{
  return sigma_crit_dis(dis_os, dis_ol, dis_ls);
}

double sigma_crit_dis(double dos, double dol, double dls)
{
  return FAC_CRITDENS * dos / (dol * dls);
}

/*--------------------------------------------------------------
  cosmology
*/

double delta_vir(double zl)
{
  double omega_vir, d_vir, omega_vir_2, omega_vir_3, eta_vir, f;
  
  omega_vir = omegaz(zl);
  omega_vir_2 = 1.0 / omega_vir - 1.0;
  
  if((lambda > 0.0) && (lambda < 1.0) && (omega > 0.0) && (omega < 1.0)){
    /* omega+lambda=1 */
    d_vir = 177.6528 * (1.0 + 0.40929 * pow(omega_vir_2, 0.90524)); 
  } else if((omega < 1.0) && (lambda < 1.0e-10) && (lambda > (-1.0e-10))){
    omega_vir_3 = 1.0 + (2.0 * omega_vir_2);
    eta_vir = log(omega_vir_3 + sqrt((omega_vir_3 * omega_vir_3) - 1.0));
    f = ((exp(eta_vir) - 1.0 / exp(eta_vir)) / 2.0) - eta_vir;
    /* omega<1,  lambda=0 */
    d_vir = 315.827 * omega_vir_2 * omega_vir_2 * omega_vir_2 / (f * f); 
  } else {
    /* omega=1,  lambda=0 or something else */
    d_vir = 177.6528; 
  }
  
  return d_vir;
}

double deltaomega(double zl)
{
  if(flag_hodensity == 1){
    /* Delta time mean matter density at z */
    return hodensity * omega * (1.0 + zl) * (1.0 + zl) * (1.0 + zl);
  } else if(flag_hodensity == 2){
    /* Delta time critical density at z */
    return hodensity * critdensz(zl);
  } else {
    /* overdensity from spherical collapse */
    return delta_vir(zl) * omega * (1.0 + zl) * (1.0 + zl) * (1.0 + zl);
  }
}

double omegaz(double zl)
{
  return omega * (1.0 + zl) * (1.0 + zl) * (1.0 + zl) / critdensz(zl);
}

double critdensz(double zl)
{
  cosmo_para cpar;

  set_cosmo_para(&cpar);

  return hubble_ez2(zl, &cpar);
}

