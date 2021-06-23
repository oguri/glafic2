#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

/*--------------------------------------------------------------
  parameters for approx. NFW and Hernquist profiles derived in M. Oguri, arXiv:2106.11464
*/

static int n_cse_nfw = 44;
static double par_cse_nfw[]={
1.082411e-06,
1.648988e-18,
8.786566e-06,
6.274458e-16,
3.292868e-06,
3.646620e-17,
1.860019e-05,
3.459206e-15,
3.274231e-05,
2.457389e-14,
6.232485e-05,
1.059319e-13,
9.256333e-05,
4.211597e-13,
1.546762e-04,
1.142832e-12,
2.097321e-04,
4.391215e-12,
3.391140e-04,
1.556500e-11,
5.178790e-04,
6.951271e-11,
8.636736e-04,
3.147466e-10,
1.405152e-03,
1.379109e-09,
2.193855e-03,
3.829778e-09,
3.179572e-03,
1.384858e-08,
4.970987e-03,
5.370951e-08,
7.631970e-03,
1.804384e-07,
1.119413e-02,
5.788608e-07,
1.827267e-02,
3.205256e-06,
2.945251e-02,
1.102422e-05,
4.562723e-02,
4.093971e-05,
6.782509e-02,
1.282206e-04,
1.596987e-01,
4.575541e-04,
1.127751e-01,
7.995270e-04,
2.169469e-01,
5.013701e-03,
3.423835e-01,
1.403508e-02,
5.194527e-01,
5.230727e-02,
8.623185e-01,
1.898907e-01,
1.382737e+00,
3.643448e-01,
2.034929e+00,
7.203734e-01,
3.402979e+00,
1.717667e+00,
5.594276e+00,
2.217566e+00,
8.052345e+00,
3.187447e+00,
1.349045e+01,
8.194898e+00,
2.603825e+01,
1.765210e+01,
4.736823e+01,
1.974319e+01,
6.559320e+01,
2.783688e+01,
1.087932e+02,
4.482311e+01,
1.477673e+02,
5.598897e+01,
2.495341e+02,
1.426485e+02,
4.305999e+02,
2.279833e+02,
7.760206e+02,
5.401335e+02,
2.143057e+03,
9.743682e+02,
1.935749e+03,
1.775124e+03
};

static int n_cse_hern = 41;
static double par_cse_hern[]={
1.199110e-06,
9.200445e-18,
3.751762e-06,
2.184724e-16,
9.927207e-06,
3.548079e-15,
2.206076e-05,
2.823716e-14,
3.781528e-05,
1.091876e-13,
6.659808e-05,
6.998697e-13,
1.154366e-04,
3.142264e-12,
1.924150e-04,
1.457280e-11,
3.040440e-04,
4.472783e-11,
4.683051e-04,
2.042079e-10,
7.745084e-04,
8.708137e-10,
1.175953e-03,
2.423649e-09,
1.675459e-03,
7.353440e-09,
2.801948e-03,
5.470738e-08,
9.712807e-03,
2.445878e-07,
5.469589e-03,
4.541672e-07,
1.104654e-02,
3.227611e-06,
1.893893e-02,
1.110690e-05,
2.792864e-02,
3.725101e-05,
4.152834e-02,
1.056271e-04,
6.640398e-02,
6.531501e-04,
1.107083e-01,
2.121330e-03,
1.648028e-01,
8.285518e-03,
2.839601e-01,
4.084190e-02,
4.129439e-01,
5.760942e-02,
8.239115e-01,
1.788945e-01,
6.031726e-01,
2.092774e-01,
1.145604e+00,
3.697750e-01,
1.401895e+00,
3.440555e-01,
2.512223e+00,
5.792737e-01,
2.038025e+00,
2.325935e-01,
4.644014e+00,
5.227961e-01,
9.301590e+00,
3.079968e-01,
2.039273e+01,
1.633456e-01,
4.896534e+01,
7.410900e-02,
1.252311e+02,
3.123329e-02,
3.576766e+02,
1.292488e-02,
2.579464e+04,
2.156527e+00,
2.944679e+04,
1.652553e-02,
2.834717e+03,
2.314934e-02,
5.931328e+04,
3.992313e-01
};

/*--------------------------------------------------------------
  approximated nfw
*/

double phi_anfw_dl(double x, double y, double q)
{
  int i;
  double f;

  f = 0.0;
  for(i=0;i<n_cse_nfw;i++){
    f = f + phi_cse_dl(x, y, par_cse_nfw[2 * i], q) * par_cse_nfw[2 * i + 1];
  }

  return f;
}

void ddphi_anfw_dl(double x, double y, double q, double *pxx, double *pxy, double *pyy)
{
  int i;
  double p1, p2, p3, pp1, pp2, pp3;

  pp1 = 0.0;
  pp2 = 0.0;
  pp3 = 0.0;

  for(i=0;i<n_cse_nfw;i++){
    ddphi_cse_dl(x, y, par_cse_nfw[2 * i], q, &p1, &p2, &p3);
    pp1 = pp1 + par_cse_nfw[2 * i + 1] * p1;
    pp2 = pp2 + par_cse_nfw[2 * i + 1] * p2;
    pp3 = pp3 + par_cse_nfw[2 * i + 1] * p3;
  }

  *pxx = pp1;
  *pxy = pp2;
  *pyy = pp3;

  return;
}

void alpha_anfw_dl(double x, double y, double q, double *ax, double *ay)
{
  int i;
  double p1, p2, pp1, pp2;

  pp1 = 0.0;
  pp2 = 0.0;

  for(i=0;i<n_cse_nfw;i++){
    alpha_cse_dl(x, y, par_cse_nfw[2 * i], q, &p1, &p2);
    pp1 = pp1 + par_cse_nfw[2 * i + 1] * p1;
    pp2 = pp2 + par_cse_nfw[2 * i + 1] * p2;
  }

  *ax = pp1;
  *ay = pp2;

  return;
}

/*--------------------------------------------------------------
  approximated hern
*/

double phi_ahern_dl(double x, double y, double q)
{
  int i;
  double f;

  f = 0.0;
  for(i=0;i<n_cse_hern;i++){
    f = f + phi_cse_dl(x, y, par_cse_hern[2 * i], q) * par_cse_hern[2 * i + 1];
  }

  return f;
}

void ddphi_ahern_dl(double x, double y, double q, double *pxx, double *pxy, double *pyy)
{
  int i;
  double p1, p2, p3, pp1, pp2, pp3;

  pp1 = 0.0;
  pp2 = 0.0;
  pp3 = 0.0;

  for(i=0;i<n_cse_hern;i++){
    ddphi_cse_dl(x, y, par_cse_hern[2 * i], q, &p1, &p2, &p3);
    pp1 = pp1 + par_cse_hern[2 * i + 1] * p1;
    pp2 = pp2 + par_cse_hern[2 * i + 1] * p2;
    pp3 = pp3 + par_cse_hern[2 * i + 1] * p3;
  }

  *pxx = pp1;
  *pxy = pp2;
  *pyy = pp3;

  return;
}

void alpha_ahern_dl(double x, double y, double q, double *ax, double *ay)
{
  int i;
  double p1, p2, pp1, pp2;

  pp1 = 0.0;
  pp2 = 0.0;

  for(i=0;i<n_cse_hern;i++){
    alpha_cse_dl(x, y, par_cse_hern[2 * i], q, &p1, &p2);
    pp1 = pp1 + par_cse_hern[2 * i + 1] * p1;
    pp2 = pp2 + par_cse_hern[2 * i + 1] * p2;
  }

  *ax = pp1;
  *ay = pp2;

  return;
}

/*--------------------------------------------------------------
  cored steep ellipsoid model
*/

double phi_cse_dl(double x, double y, double s, double q)
{
  double psi, aa;

  psi = sqrt(q * q * (s * s + x * x) + y * y);

  aa = 0.5 * log((psi + s) * (psi + s) + (1.0 - q * q) * x * x) - log((1.0 + q) * s);

  return (q / s) * aa;
  
}

void ddphi_cse_dl(double x, double y, double s, double q, double *pxx, double *pxy, double *pyy)
{
  double psi, f, pi, fi, qs;

  psi = sqrt(q * q * (s * s + x * x) + y * y);
  f = (psi + s) * (psi + s) + (1.0 - q * q) * x * x;

  pi = 1.0 / psi;
  fi = 1.0 / f;
  qs = q / s;
  
  *pxx = qs * fi * (1.0 + pi * pi * pi * q * q * s * (q * q * s * s + y * y) - 2.0 * pi * pi * fi * x * x * (psi + q * q * s) * (psi + q * q * s));
  *pyy = qs * fi * (1.0 + pi * pi * pi * q * q * s * (s * s + x * x) - 2.0 * pi * pi * fi * y * y * (psi + s) * (psi + s));
  *pxy = (-1.0) * qs * x * y * fi * (pi * pi * pi * q * q * s + 2.0 * pi * pi * fi * (psi + q * q * s) * (psi + s));

  return;
}

void alpha_cse_dl(double x, double y, double s, double q, double *ax, double *ay)
{
  double psi, f, fac;

  psi = sqrt(q * q * (s * s + x * x) + y * y);
  f = (psi + s) * (psi + s) + (1.0 - q * q) * x * x;

  fac = q / (s * psi * f);
  
  *ax = fac * x * (psi + q * q * s);
  *ay = fac * y * (psi + s);

  return;
}
