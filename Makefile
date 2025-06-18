# Makefile to build class 'iemmatrix' for Pure Data.
# Needs Makefile.pdlibbuilder as helper makefile for platform-dependent build
# settings and rules.

# library name
lib.name = iemmatrix
# use 'lib.version' to set the library version from the cmdline
#lib.version =

# mtx_*~ and friends make problems...
make-lib-executable=yes

cflags =
ldlibs =
ifneq ($(strip $(lib.version)),)
cflags += -DVERSION='"$(lib.version)"'
endif

cflags += -I. -Isrc

#####################################################################
## external dependencies (used by STUB LIBRARIES at the end)

# set any of these to 'no' to disable the library
# typically this is not needed, as
# - dependencies are only used if found
# - we do "weak" linking (so the external can still be used
#   if the dependencies are not present on the host systems
with-sndfile := yes
with-gsl := yes
with-fftw := yes

# make pkg-config overridable (for cross-building)
PKG_CONFIG ?= pkg-config

# libsndfile
ifneq ($(with-sndfile),no)
have-sndfile := $(shell $(PKG_CONFIG) --exists sndfile && echo yes || echo no)
endif
ifeq ($(have-sndfile),yes)
$(info ++++ info: found libsndfile)
SNDFILE_CFLAGS = $(shell $(PKG_CONFIG) --cflags sndfile)
SNDFILE_LIBS = $(shell $(PKG_CONFIG) --libs sndfile)
else
$(info ++++ info: missing libsndfile)
endif
ifneq ($(SNDFILE_LIBS),)
SNDFILE_CFLAGS += -DUSE_SNDFILE=1 -DHAVE_SNDFILE_H=1
else
have-sndfile=no
endif
cflags += $(SNDFILE_CFLAGS)

# GNU scientific library
ifneq ($(with-gsl),no)
have-gsl := $(shell $(PKG_CONFIG) --exists gsl && echo yes || echo no)
endif
ifeq ($(have-gsl),yes)
$(info ++++ info: found GSL)
GSL_CFLAGS = $(shell $(PKG_CONFIG) --cflags gsl)
GSL_LIBS = $(shell $(PKG_CONFIG) --libs gsl)
else
$(info ++++ info: missing GSL)
endif
ifneq ($(GSL_LIBS),)
have-gsl=yes
GSL_CFLAGS += -DHAVE_LIBGSL=1 -DHAVE_GSL_EIGEN_NONSYMM=1 -DHAVE_GSL_BESSEL=1
else
have-gsl=no
endif
cflags += $(GSL_CFLAGS)
# afaik, all libm implementations have a bessel function
cflags += -DHAVE_MATH_BESSEL=1

# FFTW
ifneq ($(with-fftw),no)
have-fftw3 := $(shell $(PKG_CONFIG) --exists fftw3 && echo yes || echo no)
endif
ifeq ($(have-fftw3),yes)
$(info ++++ info: found FFTW3)
FFTW_CFLAGS = $(shell $(PKG_CONFIG) --cflags fftw3)
FFTW_LIBS = $(shell $(PKG_CONFIG) --libs fftw3)
FFTWF_CFLAGS = $(shell $(PKG_CONFIG) --cflags fftw3f)
FFTWF_LIBS = $(shell $(PKG_CONFIG) --libs fftw3f)
else
$(info ++++ info: missing FFTW3)
endif
ifneq ($(FFTW_LIBS),)
have-fftw=yes
FFTW_CFLAGS += -DHAVE_FFTW=1
else
have-fftw=no
endif
ifneq ($(FFTWF_LIBS),)
have-fftwf=yes
FFTWF_CFLAGS += -DHAVE_FFTWF=1
else
have-fftwf=no
endif
cflags += $(FFTW_CFLAGS) $(FFTWF_CFLAGS)

##
#####################################################################

lib.setup.sources = \
	src/iemmatrix.c

common.sources = \
	src/iemmatrix_binops.c \
	src/iemmatrix_utility.c \
	src/iemmatrix_stub.c \
	$(empty)

# input source file (class name == source file basename)
class.sources = \
	src/matrix.c \
	src/mtx_abs.c \
	src/mtx_add.c \
	src/mtx_and.c \
	src/mtx_atan.c \
	src/mtx_atan2.c \
	src/mtx_bessel.c \
	src/mtx_bitand.c \
	src/mtx_bitleft.c \
	src/mtx_bitor.c \
	src/mtx_bitright.c \
	src/mtx_bspline.c \
	src/mtx_check.c \
	src/mtx_cholesky.c \
	src/mtx_col.c \
	src/mtx_colon.c \
	src/mtx_concat.c \
	src/mtx_conv.c \
	src/mtx_cos.c \
	src/mtx_cumprod.c \
	src/mtx_cumsum.c \
	src/mtx_dbtopow.c \
	src/mtx_dbtorms.c \
	src/mtx_decay.c \
	src/mtx_diag.c \
	src/mtx_diegg.c \
	src/mtx_diff.c \
	src/mtx_dispersive_dline.c \
	src/mtx_distance2.c \
	src/mtx_egg.c \
	src/mtx_eig.c \
	src/mtx_ei~.c \
	src/mtx_element.c \
	src/mtx_eq.c \
	src/mtx_exp.c \
	src/mtx_eye.c \
	src/mtx_fft.c \
	src/mtx_fill.c \
	src/mtx_find.c \
	src/mtx_gauss.c \
	src/mtx_ge.c \
	src/mtx_gt.c \
	src/mtx_ifft.c \
	src/mtx_index.c \
	src/mtx_int.c \
	src/mtx_inverse.c \
	src/mtx_isequal.c \
	src/mtx_le.c \
	src/mtx_log.c \
	src/mtx_lt.c \
	src/mtx_max2.c \
	src/mtx_mean.c \
	src/mtx_min2.c \
	src/mtx_minmax.c \
	src/mtx_mul.c \
	src/mtx_mul~.c \
	src/mtx_neq.c \
	src/mtx_not.c \
	src/mtx_ones.c \
	src/mtx_or.c \
	src/mtx_pack~.c \
	src/mtx_pivot.c \
	src/mtx_pow.c \
	src/mtx_powtodb.c \
	src/mtx_print.c \
	src/mtx_prod.c \
	src/mtx_qr.c \
	src/mtx_rand.c \
	src/mtx_repmat.c \
	src/mtx_resize.c \
	src/mtx_reverse.c \
	src/mtx_rfft.c \
	src/mtx_rifft.c \
	src/mtx_rmstodb.c \
	src/mtx_roll.c \
	src/mtx_row.c \
	src/mtx_scroll.c \
	src/mtx_sin.c \
	src/mtx_size.c \
	src/mtx_slice.c \
	src/mtx_sndfileread.c \
	src/mtx_sort.c \
	src/mtx_sub.c \
	src/mtx_sum.c \
	src/mtx_svd.c \
	src/mtx_trace.c \
	src/mtx_transpose.c \
	src/mtx_unpack~.c \
	src/mtx_zeros.c \
	$(empty)

class.sources += \
	src/mtx_qhull.c
mtx_qhull.class.sources = \
	src/mtx_qhull/list.c \
	src/mtx_qhull/vectors.c \
	src/mtx_qhull/zhull.c \
	$(empty)

class.sources += \
	src/mtx_spherical_radial.c
mtx_spherical_radial.class.sources = \
	src/mtx_spherical_harmonics/sph_radial.c \
	$(empty)

class.sources += \
	src/mtx_spherical_harmonics.c
mtx_spherical_harmonics.class.sources = \
	src/mtx_spherical_harmonics/chebyshev12.c \
	src/mtx_spherical_harmonics/legendre_a.c \
	src/mtx_spherical_harmonics/sharmonics.c \
	src/mtx_spherical_harmonics/sharmonics_normalization.c \
	$(empty)

class.sources += \
	src/mtx_convolver~.c
mtx_convolver~.class.sources = \
	src/mtx_convolver/array.c \
	src/mtx_convolver/convolver.c \
	$(empty)


#ifneq ($(make-lib-executable),yes)
#class.sources += \
#	alias/matrix_mul_line~.c \
#	alias/matrix_mul~.c \
#	alias/matrix~.c \
#	alias/mtx.c \
#	alias/mtx_0x21.c \
#	alias/mtx_0x210x3d.c \
#	alias/mtx_0x26.c \
#	alias/mtx_0x260x26.c \
#	alias/mtx_0x2a.c \
#	alias/mtx_0x2a0x7e.c \
#	alias/mtx_0x2b.c \
#	alias/mtx_0x2d.c \
#	alias/mtx_0x2e0x2a.c \
#	alias/mtx_0x2e0x2f.c \
#	alias/mtx_0x2e0x5e.c \
#	alias/mtx_0x2f.c \
#	alias/mtx_0x3a.c \
#	alias/mtx_0x3c.c \
#	alias/mtx_0x3c0x3c.c \
#	alias/mtx_0x3c0x3d.c \
#	alias/mtx_0x3d0x3d.c \
#	alias/mtx_0x3e.c \
#	alias/mtx_0x3e0x3d.c \
#	alias/mtx_0x3e0x3e.c \
#	alias/mtx_0x7c.c \
#	alias/mtx_0x7c0x7c.c \
#	alias/mtx_div.c \
#	$(empty)
#endif

# all extra files to be included in binary distribution of the library
datafiles = \
	AUTHORS.txt \
	BUGS.txt \
	GnuGPL.txt \
	LICENSE.txt \
	README.txt \
	$(empty)
datafiles = \
	VERSION.txt \
	iemmatrix-meta.pd \
	$(empty)

datafiles += \
	$(wildcard abs/*.pd) \
	$(wildcard reference/*.pd) \
	$(empty)


# include Makefile.pdlibbuilder from submodule directory 'pd-lib-builder'
PDLIBBUILDER_DIR=pd-lib-builder/
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder


src/iemmatrix.c: iemmatrix_sources.h VERSION.txt iemmatrix-meta.pd

iemmatrix_sources.h:
	src/makesource.sh $(sort $(class.sources)) >$@

%: %.in
	sed -e 's|@PACKAGE_NAME@|$(lib.name)|g' -e 's|@PACKAGE_VERSION@|$(lib.version)|g' $< > $@


.PHONY: check
check:
	-make -C tests


# build stub libraries
-include Make.stublibs
