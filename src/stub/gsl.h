#ifndef _iemmatrix_stub_gsl_h_
#define _iemmatrix_stub_gsl_h_
#ifndef HAVE_LIBGSL


typedef struct _gsl_complext_private_tag {
  double dat[2];
} gsl_complex;
typedef struct _gsl_matrix_private_tag {
  size_t size1;
  size_t size2;
  size_t tda;
  double *data;
  void *block;
  int owner;
} gsl_matrix;
typedef struct _gsl_vector_private_tag {
  size_t size;
  size_t stride;
  double *data;
  void *block;
  int owner;
} gsl_vector;
typedef struct _gsl_matrix_complex_private_tag {
  size_t size1;
  size_t size2;
  size_t tda;
  double *data;
  void *block;
  int owner;
} gsl_matrix_complex;
typedef struct _gsl_vector_complex_private_tag {
  size_t size;
  size_t stride;
  double *data;
  void *block;
  int owner;
} gsl_vector_complex;
#ifndef GSL_VECTOR_REAL
# define GSL_VECTOR_REAL(z, i)  ((z)->data[2*(i)*(z)->stride + 0])
#endif
#ifndef GSL_VECTOR_IMAG
# define GSL_VECTOR_IMAG(z, i)  ((z)->data[2*(i)*(z)->stride + 1])
#endif

typedef struct _gsl_eigen_nonsymm_workspace_private_tag gsl_eigen_nonsymm_workspace;
typedef struct _gsl_eigen_nonsymmv_workspace_private_tag gsl_eigen_nonsymmv_workspace;


#endif /* GSL */
#endif /* _iemmatrix_stub_gsl_h_ */
