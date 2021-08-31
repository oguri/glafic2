#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

void example_infile(void)
{
  printf("## example input file (for ver. %s)\n", VERSION);
  printf("## generated by glafic -d\n");
  printf("\n");
  printf("## setting primary parameters\n");
  printf("omega     %f\n", DEF_OMEGA);
  printf("lambda	  %f\n", DEF_LAMBDA); 
  printf("weos	  %f\n", DEF_WEOS); 
  printf("hubble	  %f\n", DEF_HUBBLE);
  printf("prefix	  %s\n", DEF_PREFIX);
  printf("xmin	  %f\n", DEF_XMIN);
  printf("ymin	  %f\n", DEF_YMIN);
  printf("xmax	  %f\n", DEF_XMAX);
  printf("ymax	  %f\n", DEF_YMAX);
  printf("pix_ext   %f\n", DEF_PIX_EXT);
  printf("pix_poi   %f\n", DEF_PIX_POI);
  printf("maxlev	  %d\n", DEF_MAXLEV);
  printf("\n");
  printf("## examples of secondary parameters\n");
  printf("#galfile        %s\n", DEF_GALFILE);
  printf("#srcfile        %s\n", DEF_SRCFILE);
  printf("#ran_seed       %d\n", DEF_RAN_SEED);
  printf("#outformat_exp  %d\n", DEF_OUTFORMAT_EXP);
  printf("#flag_hodensity %d\n", DEF_FLAG_HODENSITY);
  printf("#hodensity      %f\n", DEF_HODENSITY);
  printf("#gnfw_usetab    %d\n", DEF_GNFW_USETAB);
  printf("#ein_usetab     %d\n", DEF_EIN_USETAB);
  printf("#nfw_users      %d\n", DEF_NFW_USERS);
  printf("#flag_extnorm   %d\n", DEF_FLAG_EXTNORM);
  printf("#flatfix        %d\n", DEF_FLATFIX);
  printf("#chi2_checknimg %d\n", DEF_CHI2_CHECKNIMG);
  printf("#chi2_splane    %d\n", DEF_CHI2_POINT_SPLANE);
  printf("#chi2_usemag    %d\n", DEF_CHI2_USEMAG);
  printf("#chi2_restart   %d\n", DEF_CHI2_RESTART);
  printf("#obs_gain       %f\n", DEF_OBS_GAIN);
  printf("#obs_ncomb      %d\n", DEF_OBS_NCOMB);
  printf("#obs_readnoise  %f\n", DEF_OBS_READNOISE);
  printf("#skyfix         %d\n", DEF_SKYFIX);
  printf("#skyfix_value   %e\n", DEF_SKYFIX_VALUE);
  printf("#psfconv_size   %e\n", DEF_PSFCONV_SIZE);
  printf("#seeing_sub     %d\n", DEF_SEEING_SUB);
  printf("#flag_srcsbin   %d\n", DEF_FLAG_SRCSBIN);
  printf("#srcsbinsize    %f\n", DEF_SRCSBINSIZE);
  printf("#flag_mcmcall   %d\n", DEF_FLAG_MCMCALL);
  printf("#addwcs         %d\n", DEF_FLAG_ADDWCS);
  printf("#wcs_ra0        %f\n", DEF_WCS_RA0);
  printf("#wcs_dec0       %f\n", DEF_WCS_DEC0);
  printf("#ovary          0\n");
  printf("#lvary          0\n");
  printf("#wvary          0\n");
  printf("#hvary          0\n");
  printf("\n");
  printf("## define lenses and sources\n");
  printf("startup 2 2 1\n");
  printf("lens anfw   0.3 7.2e14 0.0 0.0 0.5 -45.0 6.0  0.0\n");
  printf("lens sie    0.5 300.0  2.0 2.0 0.2 -20.0 0.02 0.0\n");
  printf("extend sersic 1.5 150.0 -1.0 -1.5 0.3 90.0 0.8 1.0\n");
  printf("extend gauss  2.0 150.0  1.2  1.0 0.2 10.0 0.6 0.0\n");
  printf("point 2.5 1.0 0.5\n");
  printf("#psf 0.8 0.02 60.0 5.0 1.2 0.02 -30.0 3.0 0.7\n");
  printf("end_startup\n");
  printf("\n");
  printf("## for optimizations\n");
  printf("## can be ignored unless you do opts\n");
  printf("start_setopt\n");
  printf("0 1 0 0 1 1 0 0 \n");
  printf("0 1 0 0 1 1 0 0 \n");
  printf("0 1 1 1 1 1 0 0 \n");
  printf("0 1 1 1 1 1 0 0 \n");
  printf("0 1 1\n");
  printf("#0 0 0 0 0 0 0 0 0 \n");
  printf("end_setopt\n");
  printf("\n");
  printf("## execute commands\n");
  printf("start_command\n");
  printf("\n");
  printf("calcimage 2.5 1.0 -1.5\n");
  printf("calcein2 2.5 0.0 0.0\n");
  printf("\n");
  printf("lenscenter 2.5\n");
  printf("\n");
  printf("writelens 2.5\n");
  printf("writeimage 0.0 10.0\n");
  printf("#writepsf\n");
  printf("\n");
  printf("findimg\n");
  printf("calcextend\n");
  printf("writecrit 2.5\n");
  printf("writemesh 2.5\n");
  printf("\n");
  printf("quit\n");
  printf("\n");

  return;
}

