/*--------------------------------------------------------------
  for global parameters
*/
#ifdef	GLOBAL_SET
#define	GLOBAL
#else
#define	GLOBAL	extern
#endif

/*--------------------------------------------------------------
  gsl
*/

#include <gsl/gsl_rng.h>

/*--------------------------------------------------------------
  glafic version
*/

#define VERSION "2.0b3"
#define RELEASE_DATE "2021.4.28"

/*--------------------------------------------------------------
  primary default parameters
  can be set in an input file
*/

#define DEF_OMEGA 0.3
#define DEF_LAMBDA 0.7
#define DEF_WEOS -1.0
#define DEF_HUBBLE 0.7
#define DEF_PREFIX "out"
#define DEF_XMIN -60.0
#define DEF_YMIN -60.0
#define DEF_XMAX 60.0
#define DEF_YMAX 60.0
#define DEF_PIX_EXT 0.2
#define DEF_PIX_POI 3.0
#define DEF_MAXLEV 5

/*--------------------------------------------------------------
  secondary default parameters
*/

/* for 'gals' */
#define DEF_GALFILE "galfile.dat"
/* for 'srcs' */
#define DEF_SRCFILE "srcfile.dat"
/* seed for random number */
#define DEF_RAN_SEED -1234
/* halo overdenity 0 -> vir, 1 -> Del x mean, 2 -> Del x critical  */
#define DEF_FLAG_HODENSITY 0
/* Fixed halo overdensity Delta  */
#define DEF_HODENSITY 200.0
/* use table for gnfw profile? */
#define DEF_GNFW_USETAB 1
/* use table for Einasto profile?  */
#define DEF_EIN_USETAB 1
/* use rs instead of concentration parameter? */
#define DEF_NFW_USERS 0
/* maximum number of iteration for point source image finding */
#define DEF_NMAX_POI_ITE 10
/* accuracy of point source image finding */
#define DEF_MAX_POI_TOL 1.0e-10
/* maximum absolute value of mu^-1 beyond which subgrids are made */
#define DEF_POI_IMAG_MAX 5.0
/* minimum absolute value of mu^-1 below which subgrids are made */
#define DEF_POI_IMAG_MIN 0.001
/* for examining galaxy centers */
#define DEF_CENTER_ANG_STEP 90
/* limit of mu^-1 to avoid too large magnifications */
#define DEF_IMAG_CEIL 1.0e-10
/* small core at the center of each lens, to avoid divergence */
#define DEF_SMALLCORE 1.0e-10

/* parameters for amoeba */

/* minimum finding accuracy for amoeba (lens) */
#define DEF_TOL_AMOEBA_LENS 1.0e-4
/* minimum finding accuracy for amoeba (source) */
#define DEF_TOL_AMOEBA 1.0e-5
/* mamimum number of amoeba */
#define DEF_NMAX_AMOEBA 10000
/* mamimum number of amoeba for point source */
#define DEF_NMAX_AMOEBA_POINT 100
/* dp for mass/norm (fractional) */
#define DEF_AMOEBA_DP_MASS -0.5
/* dp for x and y in units of pixel size (extend) */
#define DEF_AMOEBA_DP_XY 10.0
/* dp for ellipticity, (1-e)*DP */
#define DEF_AMOEBA_DP_E 0.5
/* dp for angle */
#define DEF_AMOEBA_DP_ANG 45.0
/* dp for scale radius in units of pixel size (extend) */
#define DEF_AMOEBA_DP_R 10.0
/* dp for `n' parameter */
#define DEF_AMOEBA_DP_N 1.0
/* dp for redshift */
#define DEF_AMOEBA_DP_Z 0.2
/* dp for cosmological parameters */
#define DEF_AMOEBA_DP_COSMO 0.2
/* minimum length scale for amoeba */
#define DEF_AMOEBA_DELMIN 0.01
/* maximum length scale for amoeba */
#define DEF_AMOEBA_DELMAX 1.0
/* dp for PSF FWHM  */
#define DEF_AMOEBA_DP_PSFW 0.1
/* dp for PSF ellipticity */
#define DEF_AMOEBA_DP_PSFE 0.02
/* dp for PSF PA */
#define DEF_AMOEBA_DP_PSFPA 10.0
/* dp for PSF beta */
#define DEF_AMOEBA_DP_PSFB 0.5
/* dp for PSF fraction of two Moffat */
#define DEF_AMOEBA_DP_PSFF 0.05

/* parameters for optimize 
   0 -> use image plane chi^2,  1 -> use source plane chi^2 
   (for point sources) */
#define DEF_CHI2_POINT_SPLANE 0
/* check number of images or not? */
#define DEF_CHI2_CHECKNIMG 1
/* use mag instead of flux for flux constraints? */
#define DEF_CHI2_USEMAG 0
/* number of optimization */
#define DEF_CHI2_RESTART 0
/* chi2 penalty */
#define DEF_CHI2PEN_RANGE 1.0e30
/* chi2 penalty */
#define DEF_CHI2PEN_NIMG 1.0e30
/* chi2 penalty */
#define DEF_CHI2PEN_PARITY 1.0e30
/* maximum number of restarts (for chi2_restart<0) */
#define DEF_CHI2_RESTART_MAX 1000
/* gain for obs fits file [e-/ADU] */
#define DEF_OBS_GAIN 3.0
/* number of frame combine for obs fits file */
#define DEF_OBS_NCOMB 1
/* readout noise in e- */
#define DEF_OBS_READNOISE 10.0
/* sigma rejection for noise calculation */
#define DEF_NOISE_CLIP 3.0
/* fix sky value in fitting? */
#define DEF_SKYFIX 0
/* fixed sky value 
   taken from the obs data if it is not set */
#define DEF_SKYFIX_VALUE 1.0e10
/* printf format: exponential form or not? */
#define DEF_OUTFORMAT_EXP 0
/* seeing FWHM value in units of arcsec */
#define DEF_SEEING 1.0
/* range of PSF convolution (half the box size, arcsec) */
#define DEF_PSFCONV_SIZE 4.0
/* PSF ellipticity */
#define DEF_SEEING_E 0.0
/* PSF position angle */
#define DEF_SEEING_PA 0.0
/* PSF beta for Moffat, use Gaussian if <1 */
#define DEF_SEEING_BETA 3.0
/* PSF convolution subpixel */
#define DEF_SEEING_SUB 1
/* bicubic a for subpixel interpolation */
#define DEF_BICUB_A -0.5
/* range of calculating extended sources (in units of r0) */
#define DEF_SOURCE_CALCR0 20.0
/* refine extened source calculations near the center? */
#define DEF_FLAG_EXTREF 1
/* range of refining extended sources (in units of r0) */
#define DEF_SOURCE_REFR0 3.0
/* pixel integral numer in refining */
#define DEF_NUM_PIXINT 5
/* binning srcs model for speed-up */
#define DEF_FLAG_SRCSBIN 1
/* size of srcs bin, in units of arcsec */
#define DEF_SRCSBINSIZE 20.0
/* normalization of extended sources, 1 -> use total counts */
#define DEF_FLAG_EXTNORM 0
/* maximum number of sources for model `srcs' */
#define DEF_NMAX_SRCS 5000000
/* maximum number of pixels for PSF convolution */
#define DEF_NMAX_FFT 20000000
/* flag for MCMC output (include rejected points or not) */
#define DEF_FLAG_MCMCALL 0
/* flag for adding WCS in fits output */
#define DEF_FLAG_ADDWCS 0
/* RA value at x=0 */
#define DEF_WCS_RA0 150.0
/* Dec value at y=0 */
#define DEF_WCS_DEC0 30.0
/* radius to search for lens center for multiple lens planes */
#define DEF_DR_LENS_CENTER 0.9

/*--------------------------------------------------------------
  numerical parameters fixed throughout the code
*/

/* number of cosmological parameters */
#define NPAR_COSMO 4
/* number of parameters for lens, point, extend, psf, gals, srcs */
#define NPAR_LEN 8
#define NPAR_POI 3
#define NPAR_EXT 8
#define NPAR_PSF 9
#define NPAR_GAL 5
#define NPAR_SRC 7
/* number of parameters including flux and time delay zero point */
#define NPAR_POITAB 5
/* number of parameters for chi2 */
#define NPAR_CHI2 5
/* number of parameters for chi2min */
#define NPAR_CHI2MIN 5
/* number of parameters for mapprior */
#define NPAR_MPRIOR 5
/* number of parameters for lensmodel table */
#define NPAR_LMODEL 8
/* number of parameters for image table */
#define NPAR_IMAGE 4
/* number of parameters for readobs table */
#define NPAR_READOBS 7
/* number of parameters for causize in writecrit */
#define NPAR_CAUSIZE 4
/* number of parameters for def_lpl */
#define NPAR_DEFLPL 6
/* number of parameters for multi-plane lensing calculations */
#define NPAR_MULTI 13
/* maximum number of lenses, extended sources and point sources */
#define NMAX_LEN 2000
#define NMAX_EXT 1000
#define NMAX_POI 1000
/* maximum number of lens planes */
#define NMAX_LPL 20
/* tolerance of redshift difference when defining lens planes */
#define TOL_ZLPL 1.0e-3
/* maximum characters in each column of the input file */
#define INPUT_MAXCHAR 200
/* max prefix for output files */
#define PREFIX_MAXCHAR 40
/* maximum number of pixels (nx_ext or ny_ext) */
#define NMAX_PIXEL 20000
/* maximum number of total pixels for each sublevel 
   of the adaptive mesh for point sources */
#define NMAX_PIXEL_POINT 1000000000
/* maximum number of images for point source */
#define NMAX_POIMG 50
/* maximum number of images examined in mock.c */
#define NMAX_POIMG_MOCK 10
/* maxim number of maxlev */
#define NMAX_MAXLEV 20
/* maximum number of dimension for amoeba */
#define NDIMMAX 1000
/* maximum number of galaxies for model `gals' */
#define NMAX_GALS 10000
/* maximum number of bins for kappa radial profile */
#define NMAX_KAPBIN 1000
/* conversion factor between Gaussian FWHM and sigma */
#define INVSIG2FWHM 0.42466
/* maximum number of arcs for analyze */
#define NMAX_ARC 5000
/* minimum number of pixels for arc analyze */
#define NMIN_ARCANA 5
/* maximum number of mapprior */
#define NMAX_MAPPRIOR 10000
/* maximum number of lens centers */
#define NMAX_LENSCENTER 30000
/* offset in definting the number of pixels */
#define NPIX_SMALL_OFFSET 1.0e-10
/* offset to avoid diverge in del */
#define OFFSET_DEL 1.0e-30
/* offset to avoid diverge of log */
#define OFFSET_LOG 1.0e-300
/* offset to avoid diverge of tdelay_fac */
#define OFFSET_TDELAY_FAC 1.0e-300
/* intial value in peak search in extend */
#define EXTEND_INIT_PEAK -1.0e30
/* minimum error allowed in mapprior (for judging errors) */
#define TOL_ERROR_MAPPRIOR 1.0e-30
/* for zs comparison */
#define TOL_ZS 1.0e-6
/* cutoff of exp from chi2 to probability */
#define CUT_CHI2_EXP 1000.0
/* number of midway reports in mock.c */
#define NUM_MIDREPORT_MOCK 10
/* initial value of minimum chi2 in point.c */
#define CHI2_MIN_SET 1.0e60
/* initial value of minimum chi2 in point.c */
#define DIS2_MIN_SET 1.0e60
/* initial value of time delay zero point in point.c */
#define TDMIN_SET 1.0e30
/* initial value of the size of caustic in point.c */
#define CAUSIZE_SET 1.0e30
/* value of tflag in source.c */
#define TFLAG_VALUE 1234
/* initial value of xsmin etc in source.c */
#define XYSMIN_SET 1.0e30
/* initial value of chi2_dumplimit in vary.c */
#define C2DUMPLIMIT_SET 1.0e60
/* floor for computing Poisson distribution */
#define FLOOR_POISSON 1.0e-6
/* allowed range of Sersic n */
#define SERSIC_N_MIN 0.06
#define SERSIC_N_MAX 20.0
/* number when calcein fails */
#define CALCEIN_NAN -1.0

/* initial range of mass */
#define INIT_MMIN 0.0
#define INIT_MMAX 1.0e30
/* initial range of x and y */
#define INIT_XYMIN -1.0e30
#define INIT_XYMAX 1.0e30
/* initial range of ellipticity */
#define INIT_EMIN 0.0
#define INIT_EMAX 1.0
/* initial range of angles [deg] */
#define INIT_ANGMIN -360.0
#define INIT_ANGMAX  360.0
/* initial range of r0 */
#define INIT_R0MIN 0.0
#define INIT_R0MAX 1.0e30
/* initial range of n */
#define INIT_NMIN -1.0e30
#define INIT_NMAX 1.0e30
/* initial range of z */
#define INIT_ZMIN 1.0e-6
#define INIT_ZMAX 1.0e4
/* initial range of beta */
#define INIT_BETAMIN 1.0
#define INIT_BETAMAX 1.0e30
/* initial range of frac */
#define INIT_FMIN 0.0
#define INIT_FMAX 1.0
/* initial range of omega, lambda */
#define INIT_OMMIN -1.0
#define INIT_OMMAX 3.0
/* initial range of weos */
#define INIT_WMIN -5.0
#define INIT_WMAX 1.0
/* initial range of hubble */
#define INIT_HMIN 0.0
#define INIT_HMAX 3.0

/*--------------------------------------------------------------
  parameters for GSL
*/

/* parameters for calcein.c */
#define TOL_ROMBERG_AVE 1.0e-3
#define TOL_ZBRENT_CALCEIN 1.0e-7
#define XMAX_CALCEIN 1.0e4
/* parameters for distance.c */
#define GSL_CQUAD_ITERATION 100
#define GSL_ROMBERG_N_DISTANCE 20
#define TOL_CURVATURE 1.0e-6
#define GSL_EPSREL_DISTANCE 1.0e-6
/* parameters for gsl_integration.c */
#define GSL_ROMBERG_N 16
#define TOL_QGAUS 1.0e100
/* parameter for gsl_zbrent.c */
#define GSL_ZBRENT_MAXITER 100
/* parameter for mass.c */
#define TOL_ROMBERG_JHK 5.0e-4
#define TOL_ROMBERG_GNFW 3.0e-4
#define TOL_ROMBERG_EIN  1.0e-3
#define ULIM_JHK 1.0e-8

/*--------------------------------------------------------------
  default parameters for commands
*/

#define DEF_COMMAND_ZS 3.0
#define DEF_COMMAND_ZL 0.5
#define DEF_COMMAND_XY 0.0
#define DEF_COMMAND_I 0
#define DEF_COMMAND_R1 1.0
#define DEF_COMMAND_R2 10.0
#define DEF_COMMAND_N -10
#define DEF_COMMAND_SKY 0.0
#define DEF_COMMAND_NOISE 0.0
#define DEF_COMMAND_FTH 0.0
#define DEF_COMMAND_RLIM 7.5
#define DEF_COMMAND_MOCK 10
#define DEF_COMMAND_MX 1.0
#define DEF_COMMAND_MY 1.0
#define DEF_COMMAND_MR 1.0
#define DEF_COMMAND_FREC 1.0
#define DEF_COMMAND_FLAG_OUT 0
#define DEF_COMMAND_EI 1
#define DEF_COMMAND_ELLIP 0.0
#define DEF_COMMAND_C2LIMIT 1.0e30
#define DEF_COMMAND_MCMC 10000
#define DEF_COMMAND_RESET_I 1
#define DEF_COMMAND_RESET_J 1
#define DEF_COMMAND_RESET_P 1.0
#define DEF_COMMAND_RESET_F 0

/*--------------------------------------------------------------
  useful numbers
*/

#define ARCSEC2RADIAN 0.00000484813681109536
#define COVERH_MPCH 2997.92458
#define FAC_CRITDENS 5.5467248e+14
#define FAC_TDELAY_DAY 3.571386089689082e12
#define MPC2METER 3.085677581e22
#define R_SCHWARZ 2953.339382
#define C_LIGHT_KMS 2.99792458e5
#define NFW_RS_NORM 0.00009510361
#define NFW_B_NORM 6.34482175e-8
#define FAC_SIE_CALCEIN 2.08955538e19

/*--------------------------------------------------------------
  glafic.c
*/

/*--------------------------------------------------------------
  amoeba_opt.c
*/

double simplex(double v[][NDIMMAX + 1], double f[], int n, double ftol, double (*func)(double []), int *nfunc, int nmax, int verb);

/*--------------------------------------------------------------
  calcein.c
*/

void calcein(double zs);
double calcein_i_calc(int i, double zs);
double calcein_i(int i);
double calcein_jaffe(int i, double rco);
double calcein_nfw(int i);
double calcein_gnfw(int i);
double calcein_tnfw(int i);
double calcein_hern(int i);
double calcein_sers(int i);
double calcein_pow(int i);
double calcein_ein(int i);
double calcein_jaffe_func(double x);
double calcein_nfw_func(double x);
double calcein_gnfw_func(double x);
double calcein_tnfw_func(double x);
double calcein_hern_func(double x);
double calcein_sers_func(double x);
double calcein_ein_func(double x);

void kappa_rad_out(double zs, double x0, double y0, double r1, double r2, int n, int lensid);
void kappa_cum_out(double zs, double x0, double y0, double r1, double r2, int n, int lensid);
void kappa_rad(double zs, double x0, double y0, double r1, double r2, int n, int lensid, double kapbin[], double rbin[], int verb);
void kappa_cum(double zs, double x0, double y0, double r1, double r2, int n, int lensid, double kapbin[], double rbin[], int verb);
double calc_kappa_cum(double r, double x0, double y0, int lensid);
double calc_kappa_cum_func(double r);
double calc_kappa_ave(double r, double x0, double y0, int lensid);
double calc_kappa_ave_func(double t);

void calcein2(double zs, double x0, double y0, int lensid);
double calcein2_calc(double zs, double x0, double y0, int lensid);
double calcein2_func(double r);

void calcmr(void);
void calcmr_i(int i, double *mtot, double *mdel, double *rdel);

/*--------------------------------------------------------------
  call.c 
*/

void glafic_init(double in_omega, double in_lambda, double in_weos, double in_hubble, char *in_file_prefix, double in_xmin, double in_ymin, double in_xmax, double in_ymax, double in_pix_ext, double in_pix_poi, int in_maxlev, int in_ran_seed, int verb);
void glafic_set_primary(double in_omega, double in_lambda, double in_weos, double in_hubble, char *in_file_prefix, double in_xmin, double in_ymin, double in_xmax, double in_ymax, double in_pix_ext, double in_pix_poi, int in_maxlev, int verb);
void glafic_set_cosmo(double in_omega, double in_lambda, double in_weos, double in_hubble);
void glafic_quit(void);
void glafic_set_secondary(char *buffer, int verb);

void glafic_startup_setnum(int in_num_len, int in_num_ext, int in_num_poi);
void glafic_set_lens(int id, char *model, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
void glafic_set_extend(int id, char *model, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
void glafic_set_point(int id, double p1, double p2, double p3);
void glafic_set_psf(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8, double p9);
void glafic_setopt_lens(int id, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8);
void glafic_setopt_extend(int id, int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8);
void glafic_setopt_point(int id, int p1, int p2, int p3);
void glafic_setopt_psf(int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8, int p9);
void glafic_model_init(int verb);

void glafic_calcimage(double zs, double x, double y, double pout[NPAR_LMODEL], int alponly, int verb);
double glafic_calcein_i(int id, double zs);
double glafic_calcein2(int id, double zs, double x0, double y0);
double glafic_kappa_ave(int id, double zs, double r, double x0, double y0);
double glafic_kappa_cum(int id, double zs, double r, double x0, double y0);

void glafic_point_solve(double zs, double x, double y, int *ni, double rr[NMAX_POIMG][NPAR_IMAGE], int verb);
void glafic_findimg_i(int id, int *ni, double rr[NMAX_POIMG][NPAR_IMAGE], int verb);
void glafic_findimg(void);
void glafic_writelens(double zs);
void glafic_writecrit(double zs);
void glafic_writemesh(double zs);
void glafic_lenscenter(double zs);

void glafic_set_array_extend(int id, double sky, double noise, int flag_source);
void glafic_unset_array_extend(void);
void glafic_readpsf(char *fname, int verb);
double glafic_extend_array_ij(int i, int j);
double glafic_extend_array_k(int k);
void glafic_extend_array_ktoxy(int k, double *x, double *y);
void glafic_writepsf(void);

void glafic_readobs_extend(char *fname, char *fname_mask, int verb);
void glafic_readnoise_extend(char *fname, int verb);
void glafic_readobs_point(char *fname, int verb);
void glafic_parprior(char *fname, int verb);
void glafic_mapprior(char *fname, int verb);

void glafic_optimize(int verb);
double glafic_c2calc(void);

/*--------------------------------------------------------------
  commands.c
*/

int do_command(char *buffer);
int command_srcflag(char *keyword);
void interactive(void);
void deb(void);

/*--------------------------------------------------------------
  distance.c
*/

double dis_angulard(double a, double b);
double dis_luminosity(double a, double b);
double derivative_dis_proper(double z);
double dis_mod(double z);
void set_distance_lpl_zs(double zs);
void set_distance_lpl_i(int i);
void set_distance_lpl_init(void);
void set_distance_raw_zlzs(double zl, double zs);
void set_distance_raw_zl(double zl);
void set_distance_facext(void);
double disratio(double zl, double zs_fid, double zs);
double tdelay_fac(double zl, double dos, double dol, double dls);
double thetator(double theta);
double thetator_dis(double theta, double dis);
double rtotheta(double r);
double rtotheta_dis(double r, double dis);
double inv_sigma_crit(void);
double inv_sigma_crit_dis(double dos, double dol, double dls);
double sigma_crit(void);
double sigma_crit_dis(double dos, double dol, double dls);
double delta_vir(double zl);
double deltaomega(double zl);
double omegaz(double zl);
double critdensz(double zl);

/*--------------------------------------------------------------
  ein_tab.c
*/

void ein_maketable(void);
double kappa_ein_dl_tab(double x, double alpha);
double dkappa_ein_dl_tab(double x, double alpha);
double dphi_ein_dl_tab(double x, double alpha);

/*--------------------------------------------------------------
  example.c
*/

void example_infile(void);

/*--------------------------------------------------------------
  extend.c
*/

void ext_set_table_all(int i);
void ext_set_table_single(void);
void ext_set_table_lpl(int i);
void ext_set_table(int i);
void ext_unset_table(void);
void ext_set_image(int isrc, int flag_source, int verb);
double calc_pixpsf(void);
int calc_psfnk(void);
void set_psfimage(double *psf, int nk, double pix_psf);
void ext_est_image(int isrc, int flag_source, double sbth, double lwlim);
void ext_est_image_i(int isrc, double lwlim, double *sbth, double *flux, double *peak, double *px, double *py, double *area, double *flux2, double *lwmax, int *narc, FILE* fptr);
int arc_id(int k, int aid, int *im, int *ia);
double r2_extpix(int k1, int k2);
void set_psfconv_img(double *img, int ix0, int ix1, int iy0, int iy1, int in);
double bicub_func1(double x);
double bicub_func2(double x);
void ktoxy_ext(int k, double *x, double *y);
int xytok_ext(double x, double y);

/*--------------------------------------------------------------
  fits.c
*/

void writelens(double zs);
void writelens_splane(double zs, double sxmin, double sxmax, double symin, double symax, double spix);
void writeimage(double sky, double sigma, int flag_source);
double calc_pix_noise(double flux, double sky, double sigma);
void writeimageall(double sky, double sigma, int flag_source);
void writetd_ext(void);
void writetd_poi(void);
void readobs_extend(char *infile, int verb);
void readnoise_extend(char *infile, int verb);
void readmask(char *infile, int verb);
void calc_obsnoise(void);
void obs_unset_table(void);
void writenoise(void);
void addnoise(double sky, double sigma);
void writepsf(void);
void read_psffits(char *infile, int verb);

/*--------------------------------------------------------------
  gnfw_tab.c
*/

void gnfw_maketable(void);
double kappa_gnfw_dl_tab(double x, double alpha);
double dkappa_gnfw_dl_tab(double x, double alpha);
double dphi_gnfw_dl_tab(double x, double alpha);

/*--------------------------------------------------------------
  gsl_integration.c
*/

double gsl_qgaus(double (*func)(double), double a, double b);
double gsl_romberg1(double (*func)(double), double a, double b, double eps);
double gsl_romberg2(double (*func)(double), double a, double b, double eps);
double gsl_romberg3(double (*func)(double), double a, double b, double eps);

/*--------------------------------------------------------------
  gsl_zbrent.c
*/

double gsl_zbrent(double (*func)(double), double x_lo, double x_hi, double tol);

/*--------------------------------------------------------------
  init.c
*/

void init_para(char *infile);
void set_npix(void);
void init_para_body(char *keyword, char *buffer, int verb);
void init_para2(char *infile);
void init_para2_body(char *keyword, char *buffer, int verb);
void startup(char *infile);
void startup_lens(char *buffer, int n1);
void startup_extend(char *buffer, int n2);
void startup_point(char *buffer, int n3);
void startup_psf(char *buffer);
void gen_lensplane(int verb);
void setopt(char *infile);
void setopt_lens(char *buffer, int n1);
void setopt_extend(char *buffer, int n2);
void setopt_point(char *buffer, int n3);
void def_parameters(void);
void init_flags(void);
void setopt_psf(char *buffer);
void dump_model(FILE* fptr);
void dump_opt_flag(FILE* fptr);
void dump_lensplane(FILE* fptr);
void out_para(void);
void parprior(char *infile, int verb);
void mapprior(char *infile, int verb);
void readgals(void);
void readsrcs(void);
void unset_srcs(void);
void readobs_point(char *infile, int verb);

/*--------------------------------------------------------------
  mass.c
*/

void lensmodel(double tx, double ty, double pout[], int alponly, int lensid);
void multi_calc_all(double tx, double ty, double par_multi[][NPAR_MULTI], int alponly);
void multi_set_pout(double par_multi[][NPAR_MULTI], double pout[]);
void update_par_multi(double ax, double ay, double phi, double kap, double gam1, double gam2, int j, double par_multi[][NPAR_MULTI]);
void matrix_ua_calc(double par_multi[][NPAR_MULTI], int i, double par_out[]);
void lensmodel_get_i(int i, double tx, double ty, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
void lensmodel_sum(double *ax, double *ay, double *phi, double *kap, double *gam1, double *gam2, double tax, double tay, double tph, double tk, double tg1, double tg2, int alponly);
void lensmodel_sum_init(double *ax, double *ay, double *phi, double *kap, double *gam1, double *gam2);
int lmodeltoint(char *model);
char* inttolmodel(int i);

void kapgam_pert(double tx, double ty, double tx0, double ty0, double zs_fid, double k, double g, double tg, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
void kapgam_clus3(double tx, double ty, double tx0, double ty0, double zs_fid, double g, double tg, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
void kapgam_mpole(double tx, double ty, double tx0, double ty0, double zs_fid, double g, double tg, double m, double n, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
double fac_pert(double zs_fid);

void kapgam_gals(double tx, double ty, double sig, double a, double alpha, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);

void kapgam_point(double tx, double ty, double tx0, double ty0, double m, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
double re2_point(double m);

void kapgam_jaffe(double tx, double ty, double tx0, double ty0, double sig, double a, double rco, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);

void kapgam_sie(double tx, double ty, double tx0, double ty0, double sig, double s, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
void kapgam_sie_bq(double dx, double dy, double bb, double s, double q, double si, double co, double *kap, double *gam1, double *gam2, double *phi, int alponly);
void alpha_sie_bq(double dx, double dy, double bb, double s, double q, double si, double co, double *ax, double *ay);
double phi_sie_dl(double x, double y, double s, double q);
void ddphi_sie_dl(double x, double y, double s, double q, double *pxx, double *pxy, double *pyy);
void alpha_sie_dl(double x, double y, double s, double q, double *ax, double *ay);
double b_sie(double sig, double q);
double facq_sie(double q);

void kapgam_nfwpot(double tx, double ty, double tx0, double ty0, double m, double c, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
void calc_bbtt_nfw(double m, double c, double *bb, double *tt);
double ddphi_nfw_dl(double x);
double ddphi_kappadphi(double kappa, double dphi, double x);
double kappa_nfw_dl(double x);
double dkappa_nfw_dl(double x);
double dphi_nfw_dl(double x);
double phi_nfw_dl(double x);
double rs(double m, double c);
double b_func(double m, double c);
double hnfw(double c);

void kapgam_nfw(double tx, double ty, double tx0, double ty0, double m, double c, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);

void kapgam_hernpot(double tx, double ty, double tx0, double ty0, double m, double rb, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
double b_func_hern(double m, double rb);
double ddphi_hern_dl(double x);
double kappa_hern_dl(double x);
double dkappa_hern_dl(double x);
double dphi_hern_dl(double x);
double phi_hern_dl(double x);
double func_hern_dl(double x);

void kapgam_hern(double tx, double ty, double tx0, double ty0, double m, double rb, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);

void kapgam_gnfwpot(double tx, double ty, double tx0, double ty0, double m, double c, double alpha, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
void calc_bbtt_gnfw(double m, double c, double alpha, double *bb, double *tt);
double ddphi_gnfw_dl(double x);
double kappa_gnfw_dl(double x);
double kappa_gnfw_dl_func(double t);
double dkappa_gnfw_dl(double x);
double dkappa_gnfw_dl_func(double x);
double dphi_gnfw_dl(double x);
double dphi_gnfw_dl_func(double x);
double dphi_gnfw_dl_funcln(double lx);
double phi_gnfw_dl(double x);
double phi_gnfw_dl_func(double lx);
double hgnfw(double x);
double hgnfw_func(double x);
double hgnfw_funcln(double lx);
double b_func_gnfw(double m, double c, double alpha);

void kapgam_gnfw(double tx, double ty, double tx0, double ty0, double m, double c, double alpha, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);

void kapgam_powpot(double tx, double ty, double tx0, double ty0, double zs_fid, double re, double gam, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
double ddphi_pow_dl(double x);
double kappa_pow_dl(double x);
double dkappa_pow_dl(double x);
double dphi_pow_dl(double x);
double phi_pow_dl(double x);

void kapgam_pow(double tx, double ty, double tx0, double ty0, double zs_fid, double re, double gam, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);

void kapgam_serspot(double tx, double ty, double tx0, double ty0, double m, double re, double n, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
double b_func_sers(double m, double rs, double n);
double bn_sers(double n);
double bnn_sers(double n);
double gam2n1_sers(double n);
double ddphi_sers_dl(double x);
double kappa_sers_dl(double x);
double dkappa_sers_dl(double x);
double dphi_sers_dl(double x);
double phi_sers_dl(double x);
double phi_sers_dl_func(double x);

void kapgam_sers(double tx, double ty, double tx0, double ty0, double m, double re, double n, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);

void kapgam_tnfwpot(double tx, double ty, double tx0, double ty0, double m, double c, double t, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
void calc_bbtt_tnfw(double m, double c, double *bb, double *tt, double *cc);
double ddphi_tnfw_dl(double x);
double kappa_tnfw_dl(double x);
double kappa_tnfw_dl_func1(double x, double u, double t2, double t4);
double kappa_tnfw_dl_func2(double x, double u, double t2, double t4);
double dkappa_tnfw_dl(double x);
double dkappa_tnfw_dl_func1(double x, double u, double t2, double t4);
double dkappa_tnfw_dl_func2(double x, double u, double t2, double t4);
double dphi_tnfw_dl(double x);
double phi_tnfw_dl(double x);
double tnfw_mtot(double c, double t);
double func_tnfw_dl(double x);
void tnfw_set_tau(double tau);

void kapgam_tnfw(double tx, double ty, double tx0, double ty0, double m, double c, double t, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);

void kapgam_einpot(double tx, double ty, double tx0, double ty0, double m, double c, double alpha, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);
void calc_bbtt_ein(double m, double c, double alpha, double *bb, double *tt);
double ddphi_ein_dl(double x);
double kappa_ein_dl(double x);
double kappa_ein_dl_func(double t);
double dkappa_ein_dl(double x);
double dkappa_ein_dl_func(double x);
double dphi_ein_dl(double x);
double dphi_ein_dl_func(double x);
double dphi_ein_dl_funcln(double lx);
double phi_ein_dl(double x);
double phi_ein_dl_func(double lx);
double hein(double x);
double b_func_ein(double m, double c, double alpha);

void kapgam_ein(double tx, double ty, double tx0, double ty0, double m, double c, double alpha, double e, double pa, double *kap, double *gam1, double *gam2, double *phi, double *ax, double *ay, int alponly);

void u_calc(double dx, double dy, double e, double si, double co, double u[]);

double ell_integ_k(double (*func)(double), int n);
double ell_integ_k_func(double u);
double ell_integ_k_funcln(double lu);
double ell_integ_j(double (*func)(double), int n);
double ell_integ_j_func(double u);
double ell_integ_j_funcln(double lu);
double ell_integ_i(double (*func)(double));
double ell_integ_i_func(double u);
double ell_integ_i_funcln(double lu);
double ell_xi2(double u, double equ);
double ell_qu(double q, double u);
double ell_nhalf(double x, int n);
void ell_pxpy(double bpx, double bpy, double si, double co, double *px, double *by);
void ell_pxxpyy(double bpxx, double bpyy, double bpxy, double si, double co, double *pxx, double *pyy, double *pxy);

/*--------------------------------------------------------------
  mcmc.c
*/

void mcmc_calc(int n);
void paratopar(double par[]);
void mcmc_calc_step(double par[], double par2[]);
double mcmc_calc_par(int i, double x0);
void mcmc_calc_init(char *infile);

void mcmc_out_kappa_rad(char *infile, double nd, double zs, double x0, double y0, double r1, double r2, int n, int lensid);
void mcmc_out_ein(char *infile, double nd, double zs, int lensid);
void mcmc_out_ein2(char *infile, double nd, double zs, double x0, double y0, int lensid);
void mcmc_out_calcim(char *infile, double nd, double zs, double x0, double y0);

/*--------------------------------------------------------------
  mock.c
*/

void mock1(int n, double zs, double x1, double x2, double y1, double y2);
void mock2(int n, double zs, double rmax, double x0, double y0);
void mock3(int n, double zs, double fac);
void mockline(int n, double zs, double x1, double x2, double y1, double y2, int flag_full);
void mockext1(int n, int id, double x1, double x2, double y1, double y2, double e1, double e2, double sbth, double lwlim);
void mockext2(int n, int id,  double rmax, double x0, double y0, double e1, double e2, double sbth, double lwlim);
void mockext3(int n, int id, double fac, double e1, double e2, double sbth, double lwlim);

/*--------------------------------------------------------------
  opt_lens.c
*/

double opt_lens(int flag, int verb);
double chi2calc(double par[]);
double chi2calc_nopar(void);
double chi2tot(double chi2min_point[][NPAR_CHI2], double chi2min_extend[]);
void partopara(double par[]);
int opt_lens_calcndim(void);
void opt_lens_static(int flag);
void dump_opt(char *infile, double c2, int ndim, int res, int nfunc, int nd, double chi2min_point[][NPAR_CHI2], double chi2min_extend[]);
int check_para_cosmo(void);
int check_para_lens(int i, int j);
int check_para_lens_all(void);
double chi2prior_lens(void);
void parmatch_lens(void);
double chi2prior_map(void);


/*--------------------------------------------------------------
  opt_point.c
*/

double chi2calc_opt_point_out(void);
double chi2calc_opt_point(double c2min[][NPAR_CHI2], int verb);
void dump_optpoint(char *infile, double c2, int ntot, int nfunc, double c2min[][NPAR_CHI2]);
double chi2calc_opt_iplane(double c2min[][NPAR_CHI2], int verb);
double chi2calc_opt_func(double par[]);
double chi2calc_point_iplane(int i, double xs, double ys, double c2[]);
double chi2calc_opt_splane(double c2min[][NPAR_CHI2], int verb);
double chi2calc_point_splane(int i, double xs, double ys, double c2[]);
void set_matrix(double pout[NPAR_LMODEL], double mat_a[2][2], double mat_mu[2][2]);
double calcdelta_i(int i);
double calcdelta_obs(int i);
int check_para_poi(int i, double xs, double ys);
int check_para_poi_zs(int i);
int check_para_poi_all(void);
double chi2prior_point(int i, double xs, double ys);
int parmatch_poi(int i);

/*--------------------------------------------------------------
  opt_extend.c
*/

double chi2calc_opt_extend_out(void);
double chi2calc_opt_extend(double chi2min[], int verb, int flag_reset);
void dump_optextend(char *infile, double c2, int ndim, int nd, int nfunc, double chi2min[]);
double chi2calc_ext_func(double par[]);
double chi2calc_extend(double chi2min[]);
int check_para_ext(int i, int j);
int check_para_psf(int j);
int check_para_ext_all(void);
double chi2prior_ext(void);
void parmatch_ext(void);

/*--------------------------------------------------------------
  point.c
*/

void calcimg(double x, double y, double zs);
void findsrcimg(int i, double x, double y);
void poi_set_table(double zs, int flag_setlcenter, int verb);
double mag_matrix_dr2(double m11, double m12, double m21, double m22, double dx, double dy);
void poi_unset_table(void);
int magsigntot(int k, int lev);
void findimg_i(int i);
void findimg(double xs, double ys, double zs, int *ni, double rr[][NPAR_IMAGE], int verb);
double vec_product(double d1[], double d2[]);
void ktoxy_poi(int k, int lev, double *x, double *y);
void ktoxy_poi_init(int k, double *x, double *y);
void full_lens_center(double zs);
void lenscenter(double zs);
void writecrit(double zs, double causize[], int flag_newcent, int verb);
void writemesh(double zs);
int xytolev(double x, double y);
int check_mesh_center(int k, int lev, double r1, double r2);
void set_lens_center_npl0(void);
double dp_lev(int lev);
void poimg_set_table(int flag_source);
void poimg_unset_table(void);
void read_poimg_flux(char *infile);

/*--------------------------------------------------------------
  source.c
*/

double sourcemodel(double x, double y, int i, double pxx, double pxy, double pyx, double pyy, double pix);
int emodeltoint(char *model);
char* inttoemodel(int i);

double calc_srcs(double x, double y, int i, int flag_psf, double pxx, double pxy, double pyx, double pyy, double pix);
double calc_srcs_obj(int j, double x, double y, double pxx, double pxy, double pyx, double pyy, double pix);
double calc_srcs_psf(int j, double x, double y);
void tab_calc_src(void);
void unset_tab_calc_src(void);

double source_all(int id, double x, double y, double x0, double y0, double e, double pa, double r0, double n, double pxx, double pxy, double pyx, double pyy, double pix);
double source_all_norm(int id, double r0, double n, double pix);
double calc_norm_sersic(double n);
double source_psf_pix(double x, double y, double x0, double y0, double pix);
double source_psf(double x, double y, double x0, double y0, double pix);

double source_gauss(double x, double y, double x0, double y0, double e, double pa, double r0);
double source_sersic(double x, double y, double x0, double y0, double e, double pa, double r0, double n);
double source_tophat(double x, double y, double x0, double y0, double e, double pa, double r0);
double source_moffat_psf(double x, double y, double x0, double y0, double e, double pa, double fwhm, double b);
double source_moffat(double x, double y, double x0, double y0, double e, double pa, double a, double b);
double moffat_fwhmtoa(double fwhm, double b);
int source_checkdis(double x, double y, double x0, double y0, double r0);
double ucalc(double dx, double dy, double e, double pa);

/*--------------------------------------------------------------
  util.c
*/

void checkmodelpar_min(double p, double pmin);
void checkmodelpar_max(double p, double pmax);
void checkmodelpar_mineq(double p, double pmin);
void checkmodelpar_maxeq(double p, double pmax);

void terminator(char *message);

/*--------------------------------------------------------------
  vary.c
*/

void varyone(int i, int j, double pmin, double pmax, int n, int flag);
void varytwo(int i1, int j1, double pmin1, double pmax1, int n1, int i2, int j2, double pmin2, double pmax2, int n2, int flag);
void varyzs_extend(int i, double pmin, double pmax, int n, int flag);
void varyzs_point(int i, double pmin, double pmax, int n, int flag);
void varycosmo(char *keyword, double pmin, double pmax, int n, int flag);
void opt_explore(int nt, double c2lim, int flag);
void explore_parin(double bak_cosmo[], double bak_para_lens[][NPAR_LEN], double bak_para_ext[][NPAR_EXT], double bak_para_poi[][NPAR_POITAB], double bak_para_psf[NPAR_PSF]);
void explore_parout(double bak_cosmo[], double bak_para_lens[][NPAR_LEN], double bak_para_ext[][NPAR_EXT], double bak_para_poi[][NPAR_POITAB], double bak_para_psf[NPAR_PSF]);
void randomize(int verb);

/*--------------------------------------------------------------
  global parameters
*/

GLOBAL char *fname_input;

GLOBAL double delome, dis_ol, dis_os, dis_ls, zl_ext;
GLOBAL double omega, lambda, weos, hubble;
GLOBAL int num_len, num_ext, num_poi, num_gal, num_src;

GLOBAL int num_lpl, nlp; 
GLOBAL double zl_lpl[NMAX_LPL];
GLOBAL double delome_lpl[NMAX_LPL];
GLOBAL int lens_lpl_id[NMAX_LEN];
GLOBAL double dis_ol_lpl[NMAX_LPL];
GLOBAL double dis_ls_lpl[NMAX_LPL];
GLOBAL double dis_beta[NMAX_LPL][NMAX_LPL + 1];
GLOBAL double dis_tdelay[NMAX_LPL][NMAX_LPL + 1];
GLOBAL double dis_beta_fac[NMAX_LPL][NMAX_LPL];
GLOBAL double dis_tdelay_fac[NMAX_LPL][NMAX_LPL];
GLOBAL double def_lpl[NMAX_LPL + 1][NPAR_DEFLPL];
GLOBAL int i_ext_fid;

GLOBAL int ovary, lvary, wvary, hvary;
GLOBAL double omedian, lmedian, wmedian, hmedian;
GLOBAL double oerror, lerror, werror, herror;
GLOBAL double omega_min, lambda_min, weos_min, hubble_min;
GLOBAL double omega_max, lambda_max, weos_max, hubble_max;

GLOBAL double xmin, xmax, ymin, ymax;
GLOBAL double pix_ext;
GLOBAL double pix_poi;
GLOBAL int nx_ext, ny_ext;
GLOBAL int nx_poi, ny_poi;
GLOBAL int maxlev;

GLOBAL int ran_seed;
GLOBAL int flag_hodensity;
GLOBAL double hodensity;
GLOBAL int gnfw_usetab;
GLOBAL int nfw_users;
GLOBAL int ein_usetab;
GLOBAL int nmax_poi_ite;
GLOBAL double max_poi_tol;
GLOBAL double poi_imag_max;
GLOBAL double poi_imag_min;
GLOBAL int center_ang_step;
GLOBAL double imag_ceil;
GLOBAL double smallcore;
GLOBAL int outformat_exp;

GLOBAL double tol_amoeba_lens;
GLOBAL double tol_amoeba;
GLOBAL int nmax_amoeba;
GLOBAL int nmax_amoeba_point;
GLOBAL double amoeba_dp_mass;
GLOBAL double amoeba_dp_xy;
GLOBAL double amoeba_dp_e;
GLOBAL double amoeba_dp_ang;
GLOBAL double amoeba_dp_r;
GLOBAL double amoeba_dp_n;
GLOBAL double amoeba_dp_z;
GLOBAL double amoeba_dp_cosmo;
GLOBAL double amoeba_delmin;
GLOBAL double amoeba_delmax;

GLOBAL double amoeba_dp_psfw;
GLOBAL double amoeba_dp_psfe;
GLOBAL double amoeba_dp_psfpa;
GLOBAL double amoeba_dp_psfb;
GLOBAL double amoeba_dp_psff;

GLOBAL int chi2_point_splane;
GLOBAL int chi2_checknimg;
GLOBAL int chi2_usemag;
GLOBAL int chi2_restart;
GLOBAL double chi2pen_range;
GLOBAL double chi2pen_nimg;
GLOBAL double chi2pen_parity;
GLOBAL int chi2_restart_max;
GLOBAL double obs_gain;
GLOBAL int obs_ncomb;
GLOBAL double obs_readnoise;
GLOBAL double noise_clip;
GLOBAL int skyfix;
GLOBAL double skyfix_value;
GLOBAL double skymed;
GLOBAL double skysigma;
GLOBAL int flag_seeing;
GLOBAL double psfconv_size;
GLOBAL int seeing_sub;
GLOBAL double bicub_a;
GLOBAL double source_calcr0;
GLOBAL int flag_extref;
GLOBAL double source_refr0;
GLOBAL int num_pixint;
GLOBAL int flag_srcsbin;
GLOBAL double srcsbinsize;
GLOBAL int flag_extnorm;
GLOBAL int nmax_srcs;
GLOBAL int nmax_fft;
GLOBAL int flag_mcmcall;
GLOBAL double dr_lens_center;

GLOBAL double para_lens[NMAX_LEN][NPAR_LEN];
GLOBAL int model_lens[NMAX_LEN];
GLOBAL int flag_para_lens[NMAX_LEN][NPAR_LEN];
GLOBAL double para_lens_min[NMAX_LEN][NPAR_LEN];
GLOBAL double para_lens_max[NMAX_LEN][NPAR_LEN];
GLOBAL double para_lens_med[NMAX_LEN][NPAR_LEN];
GLOBAL double para_lens_sig[NMAX_LEN][NPAR_LEN];
GLOBAL double para_lens_rat[NMAX_LEN][NPAR_LEN];
GLOBAL double para_lens_ras[NMAX_LEN][NPAR_LEN];
GLOBAL int para_lens_rai[NMAX_LEN][NPAR_LEN];
GLOBAL int para_lens_raj[NMAX_LEN][NPAR_LEN];

GLOBAL double para_ext[NMAX_EXT][NPAR_EXT];
GLOBAL int model_ext[NMAX_EXT];
GLOBAL int flag_para_ext[NMAX_EXT][NPAR_EXT];
GLOBAL double dis_fac_ext[NMAX_EXT];
GLOBAL double para_ext_min[NMAX_EXT][NPAR_EXT];
GLOBAL double para_ext_max[NMAX_EXT][NPAR_EXT];
GLOBAL double para_ext_med[NMAX_EXT][NPAR_EXT];
GLOBAL double para_ext_sig[NMAX_EXT][NPAR_EXT];
GLOBAL double para_ext_rat[NMAX_EXT][NPAR_EXT];
GLOBAL double para_ext_ras[NMAX_EXT][NPAR_EXT];
GLOBAL int para_ext_rai[NMAX_EXT][NPAR_EXT];
GLOBAL int para_ext_raj[NMAX_EXT][NPAR_EXT];
GLOBAL double para_extlen_rat[NMAX_EXT][NPAR_EXT];
GLOBAL double para_extlen_ras[NMAX_EXT][NPAR_EXT];
GLOBAL int para_extlen_rai[NMAX_EXT][NPAR_EXT];
GLOBAL int para_extlen_raj[NMAX_EXT][NPAR_EXT];

GLOBAL double para_psf[NPAR_PSF];
GLOBAL int flag_para_psf[NPAR_PSF];
GLOBAL double para_psf_min[NPAR_PSF];
GLOBAL double para_psf_max[NPAR_PSF];
GLOBAL double para_psf_med[NPAR_PSF];
GLOBAL double para_psf_sig[NPAR_PSF];
GLOBAL double para_psf_rat[NPAR_PSF];
GLOBAL double para_psf_ras[NPAR_PSF];
GLOBAL int para_psf_raj[NPAR_PSF];

GLOBAL double para_poi[NMAX_POI][NPAR_POITAB];
GLOBAL int flag_para_poi[NMAX_POI][NPAR_POITAB];
GLOBAL double para_poi_min[NMAX_POI][NPAR_POITAB];
GLOBAL double para_poi_max[NMAX_POI][NPAR_POITAB];
GLOBAL double para_poi_med[NMAX_POI][NPAR_POITAB];
GLOBAL double para_poi_sig[NMAX_POI][NPAR_POITAB];
GLOBAL double para_poi_rat[NMAX_POI][NPAR_POITAB];
GLOBAL double para_poi_ras[NMAX_POI][NPAR_POITAB];
GLOBAL int para_poi_rai[NMAX_POI][NPAR_POITAB];
GLOBAL int para_poi_raj[NMAX_POI][NPAR_POITAB];
GLOBAL double para_poilen_rat[NMAX_POI][NPAR_POITAB];
GLOBAL double para_poilen_ras[NMAX_POI][NPAR_POITAB];
GLOBAL int para_poilen_rai[NMAX_POI][NPAR_POITAB];
GLOBAL int para_poilen_raj[NMAX_POI][NPAR_POITAB];

GLOBAL float para_gals[NMAX_GALS][NPAR_GAL];
GLOBAL float *para_srcs;
GLOBAL int flag_set_srcs;

GLOBAL char file_prefix[PREFIX_MAXCHAR];
GLOBAL double obs_sigma;
GLOBAL char file_gal[PREFIX_MAXCHAR];
GLOBAL char file_src[PREFIX_MAXCHAR];
GLOBAL int flag_obssig;

GLOBAL float *array_ext_def;
GLOBAL float *array_ext_mag;
GLOBAL float *array_ext_img;
GLOBAL float *array_ext_img_ori;
GLOBAL int *array_ext_mask;
GLOBAL int flag_set_array;
GLOBAL double flux_poi[NMAX_POI];

GLOBAL float *array_poi_defx[NMAX_MAXLEV];
GLOBAL float *array_poi_defy[NMAX_MAXLEV];
GLOBAL float *array_poi_smag[NMAX_MAXLEV];
GLOBAL int *array_poi_flag[NMAX_MAXLEV];
GLOBAL int *array_poi_kref[NMAX_MAXLEV];
GLOBAL int array_poi_nbox[NMAX_MAXLEV];
GLOBAL int flag_set_point;
GLOBAL float *array_poimg;
GLOBAL int flag_set_poimg;
GLOBAL float *array_imageall;

GLOBAL float *array_obs;
GLOBAL float *array_obs_noise;
GLOBAL float *array_obsnoise_file;
GLOBAL int flag_arrayobs;
GLOBAL int *array_obs_mask;
GLOBAL int flag_obsmask;
GLOBAL int obs_ext_prior[NMAX_EXT];
GLOBAL int num_mapprior;
GLOBAL double para_mapprior[NMAX_MAPPRIOR][NPAR_MPRIOR];
GLOBAL int flag_para_mapprior[NMAX_MAPPRIOR];

GLOBAL int num_lcent;
GLOBAL double lens_center[NMAX_LENSCENTER][2];
GLOBAL int lens_center_id[NMAX_LENSCENTER];

GLOBAL int obs_numimg[NMAX_POI];
GLOBAL double tab_obs[NMAX_POI][NMAX_POIMG][NPAR_READOBS];
GLOBAL int obs_parity[NMAX_POI][NMAX_POIMG];
GLOBAL int flag_pointobs;
GLOBAL double chi2_dumplimit;

GLOBAL int flag_computeall;
GLOBAL gsl_rng *ran_gsl;

GLOBAL int fpsf_nk;
GLOBAL float *array_fpsf;

GLOBAL int flag_addwcs;
GLOBAL double wcs_ra0;
GLOBAL double wcs_dec0;
