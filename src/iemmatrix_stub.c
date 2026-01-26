/*
 *  iemmatrix_stub
 *
 *  loaders for 3rd party dependencies
 *
 * Copyright (c) IOhannes m zmölnig, forum::für::umläute
 * IEM, Graz, Austria
 *
 * For information on usage and redistribution, and for a DISCLAIMER OF ALL
 * WARRANTIES, see the file, "LICENSE.txt," in this distribution.
 *
 */

#ifdef _WIN32
#include <windows.h>
typedef HINSTANCE module_t;
#define snprintf _snprintf
#else
#define _GNU_SOURCE
#include <dlfcn.h>
typedef void *module_t;
#endif

#include "iemmatrix_stub.h"
#include "m_pd.h"
#include <math.h>
#include <string.h>

#if defined(__GNUC__)
#warning iemmatrix_stub
#endif

/* dlopen() well-known libraries and get symbols */
const char *extension =
#ifdef SHARED_LIBRARY_EXTENSION
    SHARED_LIBRARY_EXTENSION
#elif defined _WIN32
    "dll"
#else
    "so"
#endif
    ;

static module_t getmodule(const char *modulename, const char *path)
{
  module_t mod = 0;
#ifdef _WIN32
  if (path)
    SetDllDirectory(path);
  mod = LoadLibrary(modulename);
#else
  // search recursively, starting from the main program
  mod = dlopen(modulename, RTLD_NOW);
  if (!mod && path)
  {
    char fullmodulename[MAXPDSTRING];
    if (snprintf(fullmodulename, MAXPDSTRING, "%s/%s", path, modulename) < MAXPDSTRING)
      mod = dlopen(fullmodulename, RTLD_NOW);
  }
#endif
  return mod;
}

static module_t getstubmodule(const char *basename, const char *path)
{
  char modulename[MAXPDSTRING];
  snprintf(modulename, MAXPDSTRING - 1, "libiemmatrixStub_%s.%s", basename, extension);
  modulename[MAXPDSTRING - 1] = 0;
  return getmodule(modulename, path);
}

typedef void *(*getfun_t)(void);
static getfun_t getfun(module_t module, const char *name)
{
  static module_t module0 = 0;
  getfun_t fun = 0;
  (void)module0;
#ifdef _WIN32
  /* fall back to pd binary */
  if (!module0)
    module0 = GetModuleHandle(0);
  if (!module)
    module = module0;

  fun = (void *)GetProcAddress(module, name);
#else
  /* fall back to pd binary */
  if (!module)
  {
#ifdef _GNU_SOURCE
    module = RTLD_DEFAULT;
#else
    if (!module0)
      module0 = dlopen(0, RTLD_NOW);
    module = module0;
#endif
  }
  fun = dlsym(module, name);
#endif
  // post("%s(%p, %s) -> %p", __FUNCTION__, module, name, fun);
  return fun;
}

/* macros for retrieving the actual functions */

#define DECLARE_STUB(fun) \
  static void *iemmatrix_stub_##fun = 0
#define MAKE_STUB(fun, module, function)                   \
  do                                                       \
  {                                                        \
    if (!iemmatrix_stub_##fun)                             \
    {                                                      \
      getfun_t f = getfun(module, "iemmatrix_" #function); \
      if (!f)                                              \
        f = getfun(module, #function);                     \
      if (f)                                               \
        iemmatrix_stub_##fun = f();                        \
    }                                                      \
  } while (0)
#define GET_STUB(module, fun, name, path) \
  do                                      \
  {                                       \
    if (!strcmp(name, #fun))              \
    {                                     \
      if (!iemmatrix_stub_##fun)          \
        setup_##module(path);             \
      return iemmatrix_stub_##fun;        \
    }                                     \
  } while (0)

/* GSL */
DECLARE_STUB(jn);
DECLARE_STUB(yn);
DECLARE_STUB(gsl_sf_bessel_Jn);
DECLARE_STUB(gsl_sf_bessel_Yn);
DECLARE_STUB(gsl_eigen_nonsymm);
DECLARE_STUB(gsl_eigen_nonsymm_alloc);
DECLARE_STUB(gsl_eigen_nonsymm_free);
DECLARE_STUB(gsl_eigen_nonsymmv);
DECLARE_STUB(gsl_eigen_nonsymmv_alloc);
DECLARE_STUB(gsl_eigen_nonsymmv_free);
DECLARE_STUB(gsl_linalg_QR_decomp);
DECLARE_STUB(gsl_linalg_QR_unpack);
DECLARE_STUB(gsl_linalg_SV_decomp);
DECLARE_STUB(gsl_matrix_alloc);
DECLARE_STUB(gsl_matrix_calloc);
DECLARE_STUB(gsl_matrix_get);
DECLARE_STUB(gsl_matrix_set);
DECLARE_STUB(gsl_matrix_complex_alloc);
DECLARE_STUB(gsl_matrix_complex_free);
DECLARE_STUB(gsl_matrix_free);
DECLARE_STUB(gsl_vector_alloc);
DECLARE_STUB(gsl_vector_complex_alloc);
DECLARE_STUB(gsl_vector_complex_free);
DECLARE_STUB(gsl_vector_free);

static void setup_gsl(const char *path)
{
  static module_t mod = 0;
  if (!mod)
    mod = getstubmodule("gsl", path);

  MAKE_STUB(jn, mod, gsl_sf_bessel_Jn);
  MAKE_STUB(yn, mod, gsl_sf_bessel_Yn);

  MAKE_STUB(gsl_sf_bessel_Jn, mod, gsl_sf_bessel_Jn);
  MAKE_STUB(gsl_sf_bessel_Yn, mod, gsl_sf_bessel_Yn);
  MAKE_STUB(gsl_eigen_nonsymm, mod, gsl_eigen_nonsymm);
  MAKE_STUB(gsl_eigen_nonsymm_alloc, mod, gsl_eigen_nonsymm_alloc);
  MAKE_STUB(gsl_eigen_nonsymm_free, mod, gsl_eigen_nonsymm_free);
  MAKE_STUB(gsl_eigen_nonsymmv, mod, gsl_eigen_nonsymmv);
  MAKE_STUB(gsl_eigen_nonsymmv_alloc, mod, gsl_eigen_nonsymmv_alloc);
  MAKE_STUB(gsl_eigen_nonsymmv_free, mod, gsl_eigen_nonsymmv_free);
  MAKE_STUB(gsl_linalg_QR_decomp, mod, gsl_linalg_QR_decomp);
  MAKE_STUB(gsl_linalg_QR_unpack, mod, gsl_linalg_QR_unpack);
  MAKE_STUB(gsl_linalg_SV_decomp, mod, gsl_linalg_SV_decomp);
  MAKE_STUB(gsl_matrix_alloc, mod, gsl_matrix_alloc);
  MAKE_STUB(gsl_matrix_calloc, mod, gsl_matrix_calloc);
  MAKE_STUB(gsl_matrix_get, mod, gsl_matrix_get);
  MAKE_STUB(gsl_matrix_set, mod, gsl_matrix_set);
  MAKE_STUB(gsl_matrix_complex_alloc, mod, gsl_matrix_complex_alloc);
  MAKE_STUB(gsl_matrix_complex_free, mod, gsl_matrix_complex_free);
  MAKE_STUB(gsl_matrix_free, mod, gsl_matrix_free);
  MAKE_STUB(gsl_vector_alloc, mod, gsl_vector_alloc);
  MAKE_STUB(gsl_vector_complex_alloc, mod, gsl_vector_complex_alloc);
  MAKE_STUB(gsl_vector_complex_free, mod, gsl_vector_complex_free);
  MAKE_STUB(gsl_vector_free, mod, gsl_vector_free);

  /* default to math.h */
  if (!iemmatrix_stub_gsl_sf_bessel_Jn)
    iemmatrix_stub_gsl_sf_bessel_Jn = jn;
  if (!iemmatrix_stub_gsl_sf_bessel_Yn)
    iemmatrix_stub_gsl_sf_bessel_Yn = yn;
  iemmatrix_stub_jn = iemmatrix_stub_gsl_sf_bessel_Jn;
  iemmatrix_stub_yn = iemmatrix_stub_gsl_sf_bessel_Yn;
}

/* sndfile */
DECLARE_STUB(sf_open);
DECLARE_STUB(sf_close);
DECLARE_STUB(sf_readf_float);

static void setup_sndfile(const char *path)
{
  static module_t mod = 0;
  if (!mod)
    mod = getstubmodule("sndfile", path);

  MAKE_STUB(sf_open, mod, sf_open);
  MAKE_STUB(sf_close, mod, sf_close);
  MAKE_STUB(sf_readf_float, mod, sf_readf_float);
}

/* fftw */
DECLARE_STUB(fftw_malloc);
DECLARE_STUB(fftw_free);
DECLARE_STUB(fftw_plan_dft_1d);
DECLARE_STUB(fftw_plan_dft_c2r_1d);
DECLARE_STUB(fftw_plan_dft_r2c_1d);
DECLARE_STUB(fftw_execute);
DECLARE_STUB(fftw_destroy_plan);

static void setup_fftw(const char *path)
{
  static module_t mod = 0;
  if (!mod)
    mod = getstubmodule("fftw", path);

  MAKE_STUB(fftw_malloc, mod, fftw_malloc);
  MAKE_STUB(fftw_free, mod, fftw_free);
  MAKE_STUB(fftw_plan_dft_1d, mod, fftw_plan_dft_1d);
  MAKE_STUB(fftw_plan_dft_c2r_1d, mod, fftw_plan_dft_c2r_1d);
  MAKE_STUB(fftw_plan_dft_r2c_1d, mod, fftw_plan_dft_r2c_1d);
  MAKE_STUB(fftw_execute, mod, fftw_execute);
  MAKE_STUB(fftw_destroy_plan, mod, fftw_destroy_plan);
}

/* fftwf */
DECLARE_STUB(fftwf_malloc);
DECLARE_STUB(fftwf_free);
DECLARE_STUB(fftwf_plan_dft_1d);
DECLARE_STUB(fftwf_plan_dft_c2r_1d);
DECLARE_STUB(fftwf_plan_dft_r2c_1d);
DECLARE_STUB(fftwf_execute);
DECLARE_STUB(fftwf_destroy_plan);

static void setup_fftwf(const char *path)
{
  static module_t mod = 0;
  if (!mod)
    mod = getstubmodule("fftwf", path);

  MAKE_STUB(fftwf_malloc, mod, fftwf_malloc);
  MAKE_STUB(fftwf_free, mod, fftwf_free);
  MAKE_STUB(fftwf_plan_dft_1d, mod, fftwf_plan_dft_1d);
  MAKE_STUB(fftwf_plan_dft_c2r_1d, mod, fftwf_plan_dft_c2r_1d);
  MAKE_STUB(fftwf_plan_dft_r2c_1d, mod, fftwf_plan_dft_r2c_1d);
  MAKE_STUB(fftwf_execute, mod, fftwf_execute);
  MAKE_STUB(fftwf_destroy_plan, mod, fftwf_destroy_plan);
}

void *iemmatrix_get_stub(const char *name, struct _class *cls)
{
  /* gsl */
  const char *path = cls ? class_gethelpdir(cls) : 0;
  GET_STUB(gsl, jn, name, path);
  GET_STUB(gsl, yn, name, path);
  GET_STUB(gsl, gsl_sf_bessel_Jn, name, path);
  GET_STUB(gsl, gsl_sf_bessel_Yn, name, path);
  GET_STUB(gsl, gsl_eigen_nonsymm, name, path);
  GET_STUB(gsl, gsl_eigen_nonsymm_alloc, name, path);
  GET_STUB(gsl, gsl_eigen_nonsymm_free, name, path);
  GET_STUB(gsl, gsl_eigen_nonsymmv, name, path);
  GET_STUB(gsl, gsl_eigen_nonsymmv_alloc, name, path);
  GET_STUB(gsl, gsl_eigen_nonsymmv_free, name, path);
  GET_STUB(gsl, gsl_linalg_QR_decomp, name, path);
  GET_STUB(gsl, gsl_linalg_QR_unpack, name, path);
  GET_STUB(gsl, gsl_linalg_SV_decomp, name, path);
  GET_STUB(gsl, gsl_matrix_alloc, name, path);
  GET_STUB(gsl, gsl_matrix_calloc, name, path);
  GET_STUB(gsl, gsl_matrix_get, name, path);
  GET_STUB(gsl, gsl_matrix_set, name, path);
  GET_STUB(gsl, gsl_matrix_complex_alloc, name, path);
  GET_STUB(gsl, gsl_matrix_complex_free, name, path);
  GET_STUB(gsl, gsl_matrix_free, name, path);
  GET_STUB(gsl, gsl_vector_alloc, name, path);
  GET_STUB(gsl, gsl_vector_complex_alloc, name, path);
  GET_STUB(gsl, gsl_vector_complex_free, name, path);
  GET_STUB(gsl, gsl_vector_free, name, path);

  /* sndfile */
  GET_STUB(sndfile, sf_open, name, path);
  GET_STUB(sndfile, sf_close, name, path);
  GET_STUB(sndfile, sf_readf_float, name, path);

  /* fftw */
  GET_STUB(fftw, fftw_malloc, name, path);
  GET_STUB(fftw, fftw_free, name, path);
  GET_STUB(fftw, fftw_plan_dft_1d, name, path);
  GET_STUB(fftw, fftw_plan_dft_c2r_1d, name, path);
  GET_STUB(fftw, fftw_plan_dft_r2c_1d, name, path);
  GET_STUB(fftw, fftw_execute, name, path);
  GET_STUB(fftw, fftw_destroy_plan, name, path);

  /* fftwf */
  GET_STUB(fftwf, fftwf_malloc, name, path);
  GET_STUB(fftwf, fftwf_free, name, path);
  GET_STUB(fftwf, fftwf_plan_dft_1d, name, path);
  GET_STUB(fftwf, fftwf_plan_dft_c2r_1d, name, path);
  GET_STUB(fftwf, fftwf_plan_dft_r2c_1d, name, path);
  GET_STUB(fftwf, fftwf_execute, name, path);
  GET_STUB(fftwf, fftwf_destroy_plan, name, path);
  return 0;
}

void iemmatrix_stub_setup()
{
  setup_gsl(0);
  setup_sndfile(0);
  setup_fftw(0);
  setup_fftwf(0);
}
