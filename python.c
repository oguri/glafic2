#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define GLOBAL_SET
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glafic.h"

/*--------------------------------------------------------------
  initialization, quit
*/

PyObject* python_init(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"omgea", "lambda", "weos", "hubble", "prefix", "xmin", "ymin", "xmax", "ymax", "pix_ext", "pix_poi", "maxlev", "ran_seed", "verb", NULL};
  double in_omega, in_lambda, in_weos, in_hubble;
  char *in_file_prefix;
  double in_xmin, in_ymin, in_xmax, in_ymax, in_pix_ext, in_pix_poi;
  int in_maxlev;
  int in_ran_seed = 0;
  int verb = 1;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ddddsddddddi|ii", argnames, &in_omega, &in_lambda, &in_weos, &in_hubble, &in_file_prefix, &in_xmin, &in_ymin, &in_xmax, &in_ymax, &in_pix_ext, &in_pix_poi, &in_maxlev, &in_ran_seed, &verb))
    return NULL;

  glafic_init(in_omega, in_lambda, in_weos, in_hubble, in_file_prefix, in_xmin, in_ymin, in_xmax, in_ymax, in_pix_ext, in_pix_poi, in_maxlev, in_ran_seed, verb);

  return Py_BuildValue("");
}

PyObject* python_set_primary(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"omgea", "lambda", "weos", "hubble", "prefix", "xmin", "ymin", "xmax", "ymax", "pix_ext", "pix_poi", "maxlev", "verb", NULL};
  double in_omega, in_lambda, in_weos, in_hubble;
  char *in_file_prefix;
  double in_xmin, in_ymin, in_xmax, in_ymax, in_pix_ext, in_pix_poi;
  int in_maxlev;
  int verb = 1;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ddddsddddddi|i", argnames, &in_omega, &in_lambda, &in_weos, &in_hubble, &in_file_prefix, &in_xmin, &in_ymin, &in_xmax, &in_ymax, &in_pix_ext, &in_pix_poi, &in_maxlev, &verb))
    return NULL;

  glafic_set_primary(in_omega, in_lambda, in_weos, in_hubble, in_file_prefix, in_xmin, in_ymin, in_xmax, in_ymax, in_pix_ext, in_pix_poi, in_maxlev, verb);

  return Py_BuildValue("");
}

PyObject* python_set_cosmo(PyObject* self, PyObject* args)
{
  double in_omega, in_lambda, in_weos, in_hubble;
  
  if(!PyArg_ParseTuple(args, "dddd", &in_omega, &in_lambda, &in_weos, &in_hubble))
    return NULL;

  glafic_set_cosmo(in_omega, in_lambda, in_weos, in_hubble);

  return Py_BuildValue("");
}

PyObject* python_quit(PyObject* self, PyObject* args)
{
  glafic_quit();

  return Py_BuildValue("");
}

PyObject* python_set_secondary(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"buffer", "verb", NULL};
  char *buffer;
  int verb = 1;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|i", argnames, &buffer, &verb))
    return NULL;

  glafic_set_secondary(buffer, verb);

  return Py_BuildValue("");
}

PyObject* python_get_nxy_ext(PyObject* self, PyObject* args)
{
  return Py_BuildValue("ii", nx_ext, ny_ext);
}

/*--------------------------------------------------------------
  setting lens model
*/

PyObject* python_startup_setnum(PyObject* self, PyObject* args)
{
  int in_num_len, in_num_ext, in_num_poi;
  
  if(!PyArg_ParseTuple(args, "iii", &in_num_len, &in_num_ext, &in_num_poi))
    return NULL;

  glafic_startup_setnum(in_num_len, in_num_ext, in_num_poi);

  return Py_BuildValue("");
}

PyObject* python_set_lens(PyObject* self, PyObject* args)
{
  int i;
  char *model;
  double p1, p2, p3, p4, p5, p6, p7, p8;
  
  if(!PyArg_ParseTuple(args, "isdddddddd", &i, &model, &p1, &p2, &p3, &p4, &p5, &p6, &p7, &p8))
    return NULL;

  glafic_set_lens(i, model, p1, p2, p3, p4, p5, p6, p7, p8);

  return Py_BuildValue("");
}

PyObject* python_set_extend(PyObject* self, PyObject* args)
{
  int i;
  char *model;
  double p1, p2, p3, p4, p5, p6, p7, p8;
  
  if(!PyArg_ParseTuple(args, "isdddddddd", &i, &model, &p1, &p2, &p3, &p4, &p5, &p6, &p7, &p8))
    return NULL;

  glafic_set_extend(i, model, p1, p2, p3, p4, p5, p6, p7, p8);

  return Py_BuildValue("");
}

PyObject* python_set_point(PyObject* self, PyObject* args)
{
  int i;
  double p1, p2, p3;
  
  if(!PyArg_ParseTuple(args, "iddd", &i, &p1, &p2, &p3))
    return NULL;

  glafic_set_point(i, p1, p2, p3);

  return Py_BuildValue("");
}

PyObject* python_set_psf(PyObject* self, PyObject* args)
{
  double p1, p2, p3, p4, p5, p6, p7, p8, p9;
  
  if(!PyArg_ParseTuple(args, "ddddddddd", &p1, &p2, &p3, &p4, &p5, &p6, &p7, &p8, &p9))
    return NULL;

  glafic_set_psf(p1, p2, p3, p4, p5, p6, p7, p8, p9);

  return Py_BuildValue("");
}

PyObject* python_setopt_lens(PyObject* self, PyObject* args)
{
  int i, p1, p2, p3, p4, p5, p6, p7, p8;
  
  if(!PyArg_ParseTuple(args, "iiiiiiiii", &i, &p1, &p2, &p3, &p4, &p5, &p6, &p7, &p8))
    return NULL;

  glafic_setopt_lens(i, p1, p2, p3, p4, p5, p6, p7, p8);

  return Py_BuildValue("");
}

PyObject* python_setopt_extend(PyObject* self, PyObject* args)
{
  int i, p1, p2, p3, p4, p5, p6, p7, p8;
  
  if(!PyArg_ParseTuple(args, "iiiiiiiii", &i, &p1, &p2, &p3, &p4, &p5, &p6, &p7, &p8))
    return NULL;

  glafic_setopt_extend(i, p1, p2, p3, p4, p5, p6, p7, p8);

  return Py_BuildValue("");
}

PyObject* python_setopt_point(PyObject* self, PyObject* args)
{
  int i, p1, p2, p3;
  
  if(!PyArg_ParseTuple(args, "iiii", &i, &p1, &p2, &p3))
    return NULL;

  glafic_setopt_point(i, p1, p2, p3);

  return Py_BuildValue("");
}

PyObject* python_setopt_psf(PyObject* self, PyObject* args)
{
  int p1, p2, p3, p4, p5, p6, p7, p8, p9;
  
  if(!PyArg_ParseTuple(args, "iiiiiiiii", &p1, &p2, &p3, &p4, &p5, &p6, &p7, &p8, &p9))
    return NULL;

  glafic_setopt_psf(p1, p2, p3, p4, p5, p6, p7, p8, p9);

  return Py_BuildValue("");
}

PyObject* python_model_init(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"verb", NULL};
  int verb = 1;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "|i", argnames, &verb))
    return NULL;

  glafic_model_init(verb);

  return Py_BuildValue("");
}

/*--------------------------------------------------------------
  lens properties
*/

PyObject* python_calcimage(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"zs", "x", "y", "alponly", "verb", NULL};
  double zs, x, y;
  int alponly = -1;
  int verb = 0;
  double pout[NPAR_LMODEL];
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ddd|ii", argnames, &zs, &x, &y, &alponly, &verb))
    return NULL;

  glafic_calcimage(zs, x, y, pout, alponly, verb);

  return Py_BuildValue("dddddddd", pout[0], pout[1], pout[2], pout[3], pout[4], pout[5], pout[6], pout[7]);
}

PyObject* python_calcein_i(PyObject* self, PyObject* args, PyObject* kwds)
{
 static char* argnames[] = {"zs", "id", NULL};
 double zs, ein;
  int id = 1;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "d|i", argnames, &zs, &id))
    return NULL;

  ein = glafic_calcein_i(id, zs);

  return Py_BuildValue("d", ein);
}

PyObject* python_calcein2(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"zs", "x0", "y0", "id", NULL};
  double zs, x0, y0, ein;
  int id = 0;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ddd|i", argnames, &zs, &x0, &y0, &id))
    return NULL;

  ein = glafic_calcein2(id, zs, x0, y0);

  return Py_BuildValue("d", ein);
}

PyObject* python_kappa_ave(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"zs", "x0", "y0", "r", "id", NULL};
  double zs, x0, y0, r, kap;
  int id = 0;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "dddd|i", argnames, &zs, &x0, &y0, &r, &id))
    return NULL;

  kap = glafic_kappa_ave(id, zs, r, x0, y0);

  return Py_BuildValue("d", kap);
}

PyObject* python_kappa_cum(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"zs", "x0", "y0", "r", "id", NULL};
  double zs, x0, y0, r, kap;
  int id = 0;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "dddd|i", argnames, &zs, &x0, &y0, &r, &id))
    return NULL;

  kap = glafic_kappa_cum(id, zs, r, x0, y0);

  return Py_BuildValue("d", kap);
}

/*--------------------------------------------------------------
  solve lens equation for point source
*/

PyObject* python_image_tuple(int ni, double rr[NMAX_POIMG][NPAR_IMAGE])
{
  int i, j;
  PyObject *tuple, *list_image;

  list_image = PyTuple_New(ni);
  if(list_image == NULL) terminator("point_solve failed");
  
  for(i=0;i<ni;i++){
    tuple = PyTuple_New(NPAR_IMAGE);
    if(tuple == NULL) terminator("point_solve failed");
    for(j=0;j<NPAR_IMAGE;j++){
      PyTuple_SET_ITEM(tuple, j, PyFloat_FromDouble(rr[i][j]));
    }
    PyTuple_SET_ITEM(list_image, i, tuple);
  } 
  
  return list_image;
}

PyObject* python_point_solve(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"zs", "x", "y", "verb", NULL};
  double zs, x, y;
  int verb = 0;
  int ni;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "ddd|i", argnames, &zs, &x, &y, &verb))
    return NULL;

  findimg(x, y, zs, &ni, rr, verb);

  return python_image_tuple(ni, rr);
}

PyObject* python_findimg_i(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"id", "verb", NULL};
  int verb = 0;
  int id, ni;
  double rr[NMAX_POIMG][NPAR_IMAGE];
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "i|i", argnames, &id, &verb))
    return NULL;

  if((id <= 0) || (id > num_poi))
    terminator("findimg_i failed");

  findimg(para_poi[id - 1][1], para_poi[id - 1][2], para_poi[id - 1][0], &ni, rr, verb);

  return python_image_tuple(ni, rr);
}

PyObject* python_findimg(PyObject* self, PyObject* args)
{
  glafic_findimg();

  return Py_BuildValue("");
}

PyObject* python_writelens(PyObject* self, PyObject* args)
{
  double zs;
  
  if(!PyArg_ParseTuple(args, "d", &zs))
    return NULL;

  glafic_writelens(zs);

  return Py_BuildValue("");
}

PyObject* python_writecrit(PyObject* self, PyObject* args)
{
  double zs;
  
  if(!PyArg_ParseTuple(args, "d", &zs))
    return NULL;

  glafic_writecrit(zs);

  return Py_BuildValue("");
}

PyObject* python_writemesh(PyObject* self, PyObject* args)
{
  double zs;
  
  if(!PyArg_ParseTuple(args, "d", &zs))
    return NULL;

  glafic_writemesh(zs);

  return Py_BuildValue("");
}

PyObject* python_lenscenter(PyObject* self, PyObject* args)
{
  double zs;
  
  if(!PyArg_ParseTuple(args, "d", &zs))
    return NULL;

  glafic_lenscenter(zs);

  return Py_BuildValue("");
}

/*--------------------------------------------------------------
  solve lens equation for extended source
*/

PyObject* python_writeimage(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"id", "sky", "noise", "ori", NULL};
  int id = 0;
  double sky = 0.0;
  double noise = 0.0;
  int flag_source = 0;
  int i, j, l, k, nn;
  double f, p;
  PyObject *tuple, *list_image;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "|iddi", argnames, &id, &sky, &noise, &flag_source))
    return NULL;

  if((id < 0) || (id > num_ext))
    terminator("writeimage failed");

  flag_computeall = 1;
  i_ext_fid = -1;
  ext_set_image(id, flag_source, 1);

  nn = nx_ext * ny_ext;
  
  list_image = PyTuple_New(ny_ext);
  if(list_image == NULL) terminator("writeimage failed");
  
  for(i=0;i<ny_ext;i++){
    tuple = PyTuple_New(nx_ext);
    if(tuple == NULL) terminator("writeimage failed");
    for(j=0;j<nx_ext;j++){
      k = j + i * nx_ext;
      if(id == 0){
	f = 0.0;
	for(l=0;l<num_ext;l++){
	  f = f + array_ext_img[k + l * nn];
	}
      } else {
	f = array_ext_img[k + (id - 1) * nn];
      }
      p = calc_pix_noise(f, sky, noise);
      
      PyTuple_SET_ITEM(tuple, j, PyFloat_FromDouble(p));
    }
    
    PyTuple_SET_ITEM(list_image, i, tuple);
  }
  
  return list_image;
}

PyObject* python_readpsf(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"fname", "verb", NULL};
  int verb = 0;
  char *fname;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|i", argnames, &fname, &verb))
    return NULL;

  glafic_readpsf(fname, verb);

  return Py_BuildValue("");
}

PyObject* python_writepsf(PyObject* self, PyObject* args)
{
  glafic_writepsf();

  return Py_BuildValue("");
}

/*------------------------------------
  reading (writing) opt-related data
*/

PyObject* python_readobs_extend(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"fname", "mask", "verb", NULL};
  int verb = 0;
  char *fname;
  char fname_mask[INPUT_MAXCHAR];

  sprintf(fname_mask, "N");
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|si", argnames, &fname, &fname_mask, &verb))
    return NULL;

  glafic_readobs_extend(fname, fname_mask, verb);

  return Py_BuildValue("");
}

PyObject* python_readnoise_extend(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"fname", "verb", NULL};
  int verb = 0;
  char *fname;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|i", argnames, &fname, &verb))
    return NULL;

  glafic_readnoise_extend(fname, verb);

  return Py_BuildValue("");
}

PyObject* python_readobs_point(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"fname", "verb", NULL};
  int verb = 0;
  char *fname;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|i", argnames, &fname, &verb))
    return NULL;

  glafic_readobs_point(fname, verb);

  return Py_BuildValue("");
}

PyObject* python_parprior(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"fname", "verb", NULL};
  int verb = 0;
  char *fname;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|i", argnames, &fname, &verb))
    return NULL;

  glafic_parprior(fname, verb);

  return Py_BuildValue("");
}

PyObject* python_mapprior(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"fname", "verb", NULL};
  int verb = 0;
  char *fname;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "s|i", argnames, &fname, &verb))
    return NULL;

  glafic_mapprior(fname, verb);

  return Py_BuildValue("");
}

/*------------------------------------
  optimizations
*/

PyObject* python_optimize(PyObject* self, PyObject* args, PyObject* kwds)
{
  static char* argnames[] = {"verb", NULL};
  int verb = 1;
  
  if(!PyArg_ParseTupleAndKeywords(args, kwds, "|i", argnames, &verb))
    return NULL;

  glafic_optimize(verb);

  return Py_BuildValue("");
}

PyObject* python_c2calc(PyObject* self, PyObject* args)
{
  double c2;

  c2 = glafic_c2calc();

  return Py_BuildValue("d", c2);
}

/*--------------------------------------------------------------
  define methods
*/

static PyMethodDef methods[] = {
  {"init", (PyCFunction)python_init, METH_VARARGS|METH_KEYWORDS},
  {"set_primary", (PyCFunction)python_set_primary, METH_VARARGS|METH_KEYWORDS},
  {"set_cosmo", python_set_cosmo, METH_VARARGS},
  {"quit", python_quit, METH_VARARGS},
  {"set_secondary", (PyCFunction)python_set_secondary, METH_VARARGS|METH_KEYWORDS},
  {"get_nxy_ext", python_get_nxy_ext, METH_VARARGS},
  {"startup_setnum", python_startup_setnum, METH_VARARGS},
  {"set_lens", python_set_lens, METH_VARARGS},
  {"set_extend", python_set_extend, METH_VARARGS},
  {"set_point", python_set_point, METH_VARARGS},
  {"set_psf", python_set_psf, METH_VARARGS},
  {"setopt_lens", python_setopt_lens, METH_VARARGS},
  {"setopt_extend", python_setopt_extend, METH_VARARGS},
  {"setopt_point", python_setopt_point, METH_VARARGS},
  {"setopt_psf", python_setopt_psf, METH_VARARGS},
  {"model_init", (PyCFunction)python_model_init, METH_VARARGS|METH_KEYWORDS},
  {"calcimage", (PyCFunction)python_calcimage, METH_VARARGS|METH_KEYWORDS},
  {"calcein_i", (PyCFunction)python_calcein_i, METH_VARARGS|METH_KEYWORDS},
  {"calcein2", (PyCFunction)python_calcein2, METH_VARARGS|METH_KEYWORDS},
  {"kappa_ave", (PyCFunction)python_kappa_ave, METH_VARARGS|METH_KEYWORDS},
  {"kappa_cum", (PyCFunction)python_kappa_cum, METH_VARARGS|METH_KEYWORDS},
  {"point_solve", (PyCFunction)python_point_solve, METH_VARARGS|METH_KEYWORDS},
  {"findimg_i", (PyCFunction)python_findimg_i, METH_VARARGS|METH_KEYWORDS},
  {"findimg", python_findimg, METH_VARARGS},
  {"writelens", python_writelens, METH_VARARGS},
  {"writecrit", python_writecrit, METH_VARARGS},
  {"writemesh", python_writemesh, METH_VARARGS},
  {"lenscenter", python_lenscenter, METH_VARARGS},
  {"writeimage", (PyCFunction)python_writeimage, METH_VARARGS|METH_KEYWORDS},
  {"readpsf", (PyCFunction)python_readpsf, METH_VARARGS|METH_KEYWORDS},
  {"writepsf", (PyCFunction)python_writepsf, METH_VARARGS|METH_KEYWORDS},
  {"readobs_extend", (PyCFunction)python_readobs_extend, METH_VARARGS|METH_KEYWORDS},
  {"readnoise_extend", (PyCFunction)python_readnoise_extend, METH_VARARGS|METH_KEYWORDS},
  {"readobs_point", (PyCFunction)python_readobs_point, METH_VARARGS|METH_KEYWORDS},
  {"parprior", (PyCFunction)python_parprior, METH_VARARGS|METH_KEYWORDS},
  {"mapprior", (PyCFunction)python_mapprior, METH_VARARGS|METH_KEYWORDS},
  {"optimize", (PyCFunction)python_optimize, METH_VARARGS|METH_KEYWORDS},
  {"c2calc", python_c2calc, METH_VARARGS},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef glafic =
{
  PyModuleDef_HEAD_INIT,
  "glafic",
  "",
  -1,
  methods
};

PyMODINIT_FUNC PyInit_glafic(void)
{
    return PyModule_Create(&glafic);
}

#undef GLOBAL
