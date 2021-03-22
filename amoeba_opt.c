#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

/* 
 * Program: nmsimplex.c
 * Author : Michael F. Hutt
 * http://www.mikehutt.com
 * Nov. 3, 1997
 *
 * An implementation of the Nelder-Mead simplex method.
 *
 * Copyright (c) 1997 Michael F. Hutt
 * Released under the MIT License
 *
 * Modified by M. Oguri, March 23, 2021
 *
 */

#define ALPHA 1.0  /* reflection coefficient */
#define BETA  0.5  /* contraction coefficient */
#define GAMMA 2.0  /* expansion coefficient */

int vg_index(double f[], int vg, int n);
int vs_index(double f[], int vs, int n);
int vh_index(double f[], int vh, int vg, int n);
void centroid(double vm[], double v[][NDIMMAX + 1], int n, int vg);

double simplex(double v[][NDIMMAX + 1], double f[], int n, double ftol, double (*func)(double []), int *nfunc, int nmax, int verb)
{
  int vs;    /* vertex with smallest value */
  int vh;    /* vertex with next smallest value */
  int vg;    /* vertex with largest value */
	
  int j, row;
  int k;     /* track the number of function evaluations */
  int itr;   /* track the number of iterations */
	
  double fr;            /* value of function at reflection point */
  double fe;            /* value of function at expansion point */
  double fc;            /* value of function at contraction point */
  double vr[NDIMMAX + 1]; /* reflection - coordinates */
  double ve[NDIMMAX + 1]; /* expansion - coordinates */
  double vc[NDIMMAX + 1]; /* contraction - coordinates */
  double vm[NDIMMAX + 1]; /* centroid - coordinates */
  double min, max, rtol;
	
  double fsum, favg, s;
	
  k = 0;

  /* begin the main loop of the minimization */
  for(itr=1;itr<=nmax;itr++) {     
    /* find the index of the largest value */
    vg = vg_index(f, 1, n);

    /* find the index of the smallest value */
    vs = vs_index(f, 1, n);
		
    /* find the index of the second largest value */
    vh = vh_index(f, vs, vg, n);
		
    /* calculate the centroid */
    centroid(vm, v, n, vg);

    /* reflect vg to new vertex vr */
    for(j=1;j<=n;j++){
      /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
      vr[j] = vm[j] + ALPHA * (vm[j] - v[vg][j]);
    }
    
    fr = func(vr);
    k++;
    
    if((fr < f[vh]) && (fr >= f[vs])){
      for(j=1;j<=n;j++){
	v[vg][j] = vr[j];
      }
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    if(fr <  f[vs]){
      for(j=1;j<=n;j++){
	/*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
	ve[j] = vm[j] + GAMMA * (vr[j] - vm[j]);
      }
      fe = func(ve);
      k++;
			
      /* 
	 by making fe < fr as opposed to fe < f[vs], 			   
	 Rosenbrocks function takes 63 iterations as opposed 
	 to 64 when using double variables. 
      */
			
      if(fe < fr){
	for(j=1;j<=n;j++){
	  v[vg][j] = ve[j];
	}
	f[vg] = fe;
      } else{
	for(j=1;j<=n;j++){
	  v[vg][j] = vr[j];
	}
	f[vg] = fr;
      }
    }
		
    /* check to see if a contraction is necessary */
    if(fr >= f[vh]){
      if((fr < f[vg]) && (fr >= f[vh])){
	/* perform outside contraction */
	for (j=1;j<=n;j++) {
	  /*vc[j] = BETA * v[vg][j] + (1 - BETA) * vm[j];*/
	  vc[j] = vm[j] + BETA * (vr[j] - vm[j]);
	}
	fc = func(vc);
	k++;
      } else{
	/* perform inside contraction */
	for(j=1;j<=n;j++){
	  /*vc[j] = BETA * v[vg][j] + (1 - BETA) * vm[j];*/
	  vc[j] = vm[j] - BETA * (vm[j] - v[vg][j]);
	}
	fc = func(vc);
	k++;
      }
			
			
      if(fc < f[vg]){
	for(j=1;j<=n;j++){
	  v[vg][j] = vc[j];
	}
	f[vg] = fc;
      } else {
	/* 
	   at this point the contraction is not successful,
	   we must halve the distance from vs to all the 
	   vertices of the simplex and then continue.
	   1997-10-31 - modified to account for ALL vertices. 
	*/
	
	for(row=1;row<=(n+1);row++){
	  if(row != vs){
	    for(j=1;j<=n;j++){
	      v[row][j] = v[vs][j] + (v[row][j] - v[vs][j]) / 2.0;
	    }
	  }
	}

	/* re-evaluate all the vertices */
	for(j=1;j<=(n+1);j++){
	  f[j] = func(v[j]);
	}
	
	/* find the index of the largest value */
	vg = vg_index(f, 1, n);
	
	/* find the index of the smallest value */
	vs = vs_index(f, 1, n);
	
	/* find the index of the second largest value */
	vh = vh_index(f, vs, vg, n);

	f[vg] = func(v[vg]);
	k++;

	f[vh] = func(v[vh]);
	k++;
      }
    }
		
    /* test for convergence */
    fsum = 0.0;
    for(j=1;j<=(n+1);j++){
      fsum += f[j];
    }
    favg = fsum / (n + 1);
    s = 0.0;
    for(j=1;j<=(n+1);j++){
      s += pow((f[j] - favg), 2.0) / (n);
    }
    s = sqrt(s);

    max = f[vg_index(f, 1, n)];
    min = f[vs_index(f, 1, n)];
    rtol = 2.0 * fabs(max - min) / (fabs(max) + fabs(min));

    if(verb > 0) fprintf(stderr, "amoeba:  n = %5d  y = %e  tol = %e\n", itr, min, rtol);

    /* if(s < ftol) break; */
    if((rtol < ftol) || (fabs(max) < (0.1 * ftol))) break;
  }
  /* end main loop of the minimization */
	
  /* find the index of the smallest value */
  vs = vs_index(f, 1, n);
	
  /*printf("The minimum was found at\n"); */
  /* set minimum to i = 1 */
  for(j=1;j<(n+1);j++){
    /*printf("%e\n", v[vs][j]);*/
    v[1][j] = v[vs][j];
  }
  min = func(v[vs]);
  f[1] = min;
  k++;

  *nfunc = itr;
  
  return min;
}

/* find the index of the largest value */
int vg_index(double f[], int vg, int n)
{
  int j;
  
  for(j=1;j<=(n+1);j++){
    if(f[j] > f[vg]){
      vg = j;
    }
  }
  
  return vg;
}

/* find the index of the smallest value */
int vs_index(double f[], int vs, int n)
{
  int j;
    
  for(j=1;j<=(n+1);j++){
    if(f[j] < f[vs]){
      vs = j;
    }
  }
  
  return vs;
}

/* find the index of the second largest value */
int vh_index(double f[], int vh, int vg, int n)
{
  int j;
  
  for(j=1;j<=(n+1);j++){
    if((f[j] > f[vh]) && (f[j] < f[vg])){
      vh = j;
    }
  }
  
  return vh;
}

/* calculate the centroid */
void centroid(double vm[], double v[][NDIMMAX + 1], int n, int vg)
{
  int j, m;
  double cent;
  
  for(j=1;j<=n;j++){
    cent = 0.0;
    for(m=1;m<=(n+1);m++){
      if(m != vg){
	cent += v[m][j];
      }
    }
    vm[j] = cent / n;
  }

  return;
}

#undef ALPHA 
#undef BETA  
#undef GAMMA 
