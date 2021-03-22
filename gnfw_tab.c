#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

#define NUM_LNX 801
#define LNX_MIN -23.0
#define DLNX 0.0575 
#define NUM_ALP 101
#define ALP_MIN 0.0
#define DALP 0.02
#define NUM_FLAG 13

double intpol_gnfw(double tab[NUM_ALP][NUM_LNX], double lnx, double alp);
double intpol_gnfw_lin(double tab[NUM_ALP][NUM_LNX], double lnx, double alp);

static double tab_alp[NUM_ALP];
static double tab_lnx[NUM_LNX];
static double tab_lnkappa_gnfw[NUM_ALP][NUM_LNX];
static double tab_lndkappa_gnfw[NUM_ALP][NUM_LNX];
static double tab_lndphi_gnfw[NUM_ALP][NUM_LNX];

/*--------------------------------------------------------------
  make tables
*/

void gnfw_maketable(void)
{
  int i, j, f;
  double x, lx, alp;
  static int flag;

 if(flag != NUM_FLAG){
    flag = NUM_FLAG;

    fprintf(stderr, " making tables for gnfw... [%d x %d = %d points]\n\n", NUM_ALP, NUM_LNX, NUM_ALP * NUM_LNX);

    f = gnfw_usetab;
    gnfw_usetab = 0;

    for(i=0;i<NUM_ALP;i++){
      tab_alp[i] = (ALP_MIN) + ((double)i) * (DALP);
    }
    
    for(i=0;i<NUM_LNX;i++){
      tab_lnx[i] = (LNX_MIN) + ((double)i) * (DLNX);
    }
    
    for(i=0;i<NUM_ALP;i++){
      alp = tab_alp[i];
      b_func_gnfw(1.0, 1.0, alp); /*  set alpha value  */
      for(j=0;j<NUM_LNX;j++){
	lx = tab_lnx[j];
	x = exp(lx);
	tab_lnkappa_gnfw[i][j] = log(kappa_gnfw_dl(x) + OFFSET_LOG);
	tab_lndkappa_gnfw[i][j] = log((-1.0) * dkappa_gnfw_dl(x) + OFFSET_LOG);
	tab_lndphi_gnfw[i][j] = log(dphi_gnfw_dl(x) + OFFSET_LOG);
      }
    }

    gnfw_usetab = f;
  }

  return;
}

/*--------------------------------------------------------------
  estimate values from the tables
*/

double kappa_gnfw_dl_tab(double x, double alpha)
{
  double lx, ll;

  lx = log(x);
  ll = intpol_gnfw(tab_lnkappa_gnfw, lx, alpha);

  return exp(ll);
}

double dkappa_gnfw_dl_tab(double x, double alpha)
{
  double lx, ll;

  lx = log(x);
  ll = intpol_gnfw(tab_lndkappa_gnfw, lx, alpha);

  return (-1.0) * exp(ll);
}

double dphi_gnfw_dl_tab(double x, double alpha)
{
  double lx, ll;

  lx = log(x);
  ll = intpol_gnfw(tab_lndphi_gnfw, lx, alpha);

  return exp(ll);
}

double intpol_gnfw(double tab[NUM_ALP][NUM_LNX], double lnx, double alp)
{
  return intpol_gnfw_lin(tab, lnx, alp);
  /* return intpol_gnfw_bcu(tab, lnx, alp); */
}

double intpol_gnfw_lin(double tab[NUM_ALP][NUM_LNX], double lnx, double alp)
{
  int i, j;
  double t, u;

  i = (int)((alp - (ALP_MIN)) / DALP);
  j = (int)((lnx - (LNX_MIN)) / DLNX);

  if(i < 0){ i = 0; alp = tab_alp[i]; }
  if(j < 0){ j = 0; lnx = tab_lnx[j]; }
  if(i >= NUM_ALP){ i = NUM_ALP - 1; alp = tab_alp[i]; }
  if(j >= NUM_LNX){ j = NUM_LNX - 1; lnx = tab_lnx[j]; }
  
  t = (alp - tab_alp[i]) / (tab_alp[i + 1] - tab_alp[i]);
  u = (lnx - tab_lnx[j]) / (tab_lnx[j + 1] - tab_lnx[j]);
  
  return (1.0 - t) * (1.0 - u) * tab[i][j] + t * (1.0 - u) * tab[i + 1][j] + t * u * tab[i + 1][j + 1] + (1.0 - t) * u * tab[i][j + 1];
}

#undef NUM_LNX 
#undef LNX_MIN 
#undef DLNX 
#undef NUM_ALP 
#undef ALP_MIN 
#undef DALP 
#undef NUM_FLAG 

