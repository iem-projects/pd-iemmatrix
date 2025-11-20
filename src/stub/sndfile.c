#ifdef HAVE_SNDFILE_H
# include <sndfile.h>

# undef STUB
# define STUB(x) \
  void* iemmatrix_##x() { return x; }

STUB(sf_open);
STUB(sf_close);
STUB(sf_readf_float);

#endif
