#define GLOBAL_SET
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

/*--------------------------------------------------------------
   main function
*/

int main(int argc, char **argv)
{
  int i;
  char buffer[INPUT_MAXCHAR];
  char keyword[INPUT_MAXCHAR];
  FILE* fptr;
  
  init_flags();

  fprintf(stderr, "\n");
  fprintf(stderr, "glafic  ver. %s (%s) by Masamune Oguri\n", VERSION, RELEASE_DATE);
  fprintf(stderr, "\n");
  
  if(argc != 2) terminator("please specify an input file \n > glafic file.input \nexample of file.input can be printed by -d option \n > glafic -d\n");
  
  if(strcmp(argv[1], "-d") == 0){
    example_infile();
    return 0;
  }

  fname_input = argv[1];
  
  /* read parameters */
  init_para(argv[1], 1);
  init_para2(argv[1], 1);
  /* startup */
  startup(argv[1], 1);
  /* set para for opt */
  setopt(argv[1], 1);
  /* prepare lens planes */
  gen_lensplane(1);
  
  ran_gsl = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(ran_gsl, ran_seed);

  set_distance_lpl_init();
  
  fprintf(stderr, "\n######## run commands\n\n");

  fptr = fopen(argv[1], "r");
  
  i = 0;

  while(fgets(buffer, INPUT_MAXCHAR, fptr)){
    sscanf(buffer, "%s", keyword);
    if(strcmp(keyword, "start_command") == 0) break;
  }
  
  while(fgets(buffer, INPUT_MAXCHAR, fptr)){
    if(sscanf(buffer, "%s",keyword) != EOF){
      if(keyword[0] != '#'){
	fprintf(stderr, "glafic> %s\n", buffer);
	i = do_command(buffer);
      }
    }
    if(i == 1) break;
  }
  
  if(i == 0) interactive();

  ext_unset_table();
  obs_unset_table();
  poi_unset_table();
  poimg_unset_table();
  unset_srcs();
  unset_tab_calc_src();
  gsl_rng_free(ran_gsl);

  fclose(fptr);

  return 0;
}

#undef GLOBAL

