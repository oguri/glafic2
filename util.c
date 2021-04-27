#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

/*--------------------------------------------------------------
  check parameter range
*/

void checkmodelpar_min(double p, double pmin)
{
  if(p > pmin){
    return;
  } else {
    terminator("model parameter out of range");
  }
}

void checkmodelpar_max(double p, double pmax)
{
  if(p < pmax){
    return;
  } else {
    terminator("model parameter out of range");
  }
}

void checkmodelpar_mineq(double p, double pmin)
{
  if(p >= pmin){
    return;
  } else {
    terminator("model parameter out of range");
  }
}

void checkmodelpar_maxeq(double p, double pmax)
{
  if(p <= pmax){
    return;
  } else {
    terminator("model parameter out of range");
  }
}

/*--------------------------------------------------------------
  damn it! 
*/

void terminator(char *message)
{
  fprintf(stderr, "Error: %s\n", message);
  
  ext_unset_table();
  obs_unset_table();
  poi_unset_table();
  poimg_unset_table();
  gsl_rng_free(ran_gsl);
  
  exit(EXIT_FAILURE);

  return;
}

