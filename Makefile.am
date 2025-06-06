AUTOMAKE_OPTIONS = foreign subdir-objects
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = abs doc

iemmatrixdir=$(pkglibdir)
dist_iemmatrix_DATA = \
	AUTHORS.txt  BUGS.txt  GnuGPL.txt  INSTALL.txt  LICENSE.txt  VERSION.txt \
	iemmatrix-meta.pd

pkglib_LTLIBRARIES=

AM_CPPFLAGS=-I$(top_srcdir)/src @IEMMATRIX_CPPFLAGS@
AM_CFLAGS =
AM_LDFLAGS= -module -avoid-version -shared -shrext .@EXTENSION@

BUILT_SOURCES=
CLEANFILES   =

if WINDOWS
AM_LDFLAGS+= -no-undefined
endif

AM_LDFLAGS+= @EXTRA_LTFLAGS@

# RTE
AM_CFLAGS +=
AM_LDFLAGS+=


noinst_LTLIBRARIES = \
	libiemmatrixCore.la \
	libiemmatrixSphericalHarmonics.la libiemmatrixSphericalRadial.la \
	libiemmatrixQhull.la \
	libiemmatrixConvolver.la \
	$(empty)

libiemmatrixCore_la_SOURCES=\
	src/iemmatrix.h \
	src/mtx_binop_generic.h \
	src/iemmatrix_binops.c \
	src/iemmatrix_utility.c \
	$(empty)

libiemmatrixSphericalHarmonics_la_SOURCES =  \
	src/mtx_spherical_harmonics/chebyshev12.c src/mtx_spherical_harmonics/chebyshev12.h \
	src/mtx_spherical_harmonics/legendre_a.c src/mtx_spherical_harmonics/legendre_a.h \
	src/mtx_spherical_harmonics/sharmonics.c src/mtx_spherical_harmonics/sharmonics.h \
	src/mtx_spherical_harmonics/sharmonics_normalization.c src/mtx_spherical_harmonics/sharmonics_normalization.h \
	$(empty)

libiemmatrixSphericalRadial_la_SOURCES =  \
	src/mtx_spherical_harmonics/sph_radial.c src/mtx_spherical_harmonics/sph_radial.h \
	$(empty)


## qhull-stuff
libiemmatrixQhull_la_SOURCES = \
	src/mtx_qhull/entry.h \
	src/mtx_qhull/list.h \
        src/mtx_qhull/list.c \
	src/mtx_qhull/vectors.h \
	src/mtx_qhull/vectors.c \
	src/mtx_qhull/zhull.h \
	src/mtx_qhull/zhull.c \
	$(empty)

## convolver~-stuff
libiemmatrixConvolver_la_SOURCES = \
	src/mtx_convolver/convolver.c \
	src/mtx_convolver/convolver.h \
	src/mtx_convolver/array.h \
        src/mtx_convolver/array.c \
	$(empty)

libiemmatrixConvolver_la_CFLAGS = $(FFTW3F_CFLAGS)
libiemmatrixConvolver_la_LIBADD = $(FFTW3F_LIBS)

iemmatrix_la_CPPFLAGS=-DIEMMATRIX_LIBRARY $(AM_CPPFLAGS)
iemmatrix_la_CFLAGS =$(AM_CFLAGS)
iemmatrix_la_CFLAGS+=$(SNDFILE_CFLAGS) $(FFTW3_CFLAGS) $(GSL_CFLAGS)
iemmatrix_la_LIBADD = libiemmatrixCore.la
iemmatrix_la_LIBADD+= libiemmatrixSphericalHarmonics.la libiemmatrixSphericalRadial.la
iemmatrix_la_LIBADD+= libiemmatrixQhull.la
iemmatrix_la_LIBADD+= libiemmatrixConvolver.la



if LIBRARY
pkglib_LTLIBRARIES+=iemmatrix.la
endif LIBRARY
if SINGLEOBJECTS
pkglib_LTLIBRARIES += \
	matrix.la mtx_abs.la mtx_add.la \
	mtx_and.la mtx_atan.la mtx_atan2.la \
	mtx_bitand.la mtx_bitleft.la \
	mtx_bitor.la mtx_bitright.la mtx_bspline.la \
	mtx_check.la mtx_cholesky.la mtx_col.la \
	mtx_colon.la mtx_concat.la mtx_conv.la \
	mtx_convoler~.la \
	mtx_cos.la mtx_cumprod.la mtx_cumsum.la \
	mtx_dbtopow.la mtx_dbtorms.la mtx_decay.la \
	mtx_diag.la mtx_diegg.la mtx_diff.la \
	mtx_dispersive_dline.la mtx_distance2.la mtx_egg.la \
	mtx_ei~.la mtx_element.la \
	mtx_eq.la mtx_exp.la mtx_eye.la \
	mtx_fft.la mtx_fill.la mtx_find.la \
	mtx_gauss.la mtx_ge.la mtx_gt.la \
	mtx_ifft.la mtx_index.la mtx_int.la \
	mtx_inverse.la mtx_isequal.la mtx_le.la \
	mtx_log.la mtx_lt.la mtx_max2.la \
	mtx_mean.la mtx_min2.la mtx_minmax.la \
	mtx_mul.la mtx_mul~.la mtx_neq.la \
	mtx_not.la mtx_ones.la mtx_or.la \
	mtx_pack~.la mtx_pivot.la mtx_pow.la \
	mtx_powtodb.la mtx_print.la mtx_prod.la \
	mtx_qhull.la mtx_rand.la \
	mtx_repmat.la mtx_resize.la mtx_reverse.la \
	mtx_rmstodb.la \
	mtx_roll.la mtx_row.la mtx_scroll.la \
	mtx_sin.la mtx_size.la mtx_slice.la \
	mtx_sort.la mtx_spherical_harmonics.la \
	mtx_spherical_radial.la mtx_sub.la mtx_sum.la \
	mtx_trace.la mtx_transpose.la \
	mtx_unpack~.la mtx_zeros.la matrix_mul_line~.la \
	matrix_mul~.la matrix~.la mtx.la \
	mtx_0x21.la mtx_0x210x3d.la mtx_0x26.la \
	mtx_0x260x26.la mtx_0x2a.la mtx_0x2a0x7e.la \
	mtx_0x2b.la mtx_0x2d.la mtx_0x2e0x2a.la \
	mtx_0x2e0x2f.la mtx_0x2e0x5e.la mtx_0x2f.la \
	mtx_0x3a.la mtx_0x3c.la mtx_0x3c0x3c.la \
	mtx_0x3c0x3d.la mtx_0x3d0x3d.la mtx_0x3e.la \
	mtx_0x3e0x3d.la mtx_0x3e0x3e.la mtx_0x7c.la \
	mtx_0x7c0x7c.la mtx_div.la \
	$(empty)
endif SINGLEOBJECTS

# these objects are always built separately
# hmm, these can either use fftw3 or Pd's mayerFFT
# LATER provide both binaries, preferring fftw3 if available at runtime
pkglib_LTLIBRARIES += mtx_rfft.la mtx_rifft.la mtx_convolver~.la

# depends on libsndfile
pkglib_LTLIBRARIES += mtx_sndfileread.la
# depends on libgsl
pkglib_LTLIBRARIES += \
	mtx_bessel.la \
	mtx_eig.la \
	mtx_qr.la \
	mtx_svd.la \
	$(empty)

iemmatrix_la_SOURCES = \
	src/iemmatrix.c \
	iemmatrix_sources.h \
	$(empty)
iemmatrix_la_SOURCES += \
	src/matrix.c src/mtx_abs.c src/mtx_add.c \
	src/mtx_and.c src/mtx_atan.c src/mtx_atan2.c \
	src/mtx_bitand.c src/mtx_bitleft.c \
	src/mtx_bitor.c src/mtx_bitright.c src/mtx_bspline.c \
	src/mtx_check.c src/mtx_cholesky.c src/mtx_col.c \
	src/mtx_colon.c src/mtx_concat.c src/mtx_conv.c \
	src/mtx_cos.c src/mtx_cumprod.c src/mtx_cumsum.c \
	src/mtx_dbtopow.c src/mtx_dbtorms.c src/mtx_decay.c \
	src/mtx_diag.c src/mtx_diegg.c src/mtx_diff.c \
	src/mtx_dispersive_dline.c src/mtx_distance2.c src/mtx_egg.c \
	src/mtx_ei~.c src/mtx_element.c \
	src/mtx_eq.c src/mtx_exp.c src/mtx_eye.c \
	src/mtx_fft.c src/mtx_fill.c src/mtx_find.c \
	src/mtx_gauss.c src/mtx_ge.c src/mtx_gt.c \
	src/mtx_ifft.c src/mtx_index.c src/mtx_int.c \
	src/mtx_inverse.c src/mtx_isequal.c src/mtx_le.c \
	src/mtx_log.c src/mtx_lt.c src/mtx_max2.c \
	src/mtx_mean.c src/mtx_min2.c src/mtx_minmax.c \
	src/mtx_mul.c src/mtx_mul~.c src/mtx_neq.c \
	src/mtx_not.c src/mtx_ones.c src/mtx_or.c \
	src/mtx_pack~.c src/mtx_pivot.c src/mtx_pow.c \
	src/mtx_powtodb.c src/mtx_print.c src/mtx_prod.c \
	src/mtx_qhull.c src/mtx_rand.c \
	src/mtx_repmat.c src/mtx_resize.c src/mtx_reverse.c \
	src/mtx_rmstodb.c \
	src/mtx_roll.c src/mtx_row.c src/mtx_scroll.c \
	src/mtx_sin.c src/mtx_size.c src/mtx_slice.c \
	src/mtx_sort.c src/mtx_spherical_harmonics.c \
	src/mtx_spherical_radial.c src/mtx_sub.c src/mtx_sum.c \
	src/mtx_trace.c src/mtx_transpose.c \
	src/mtx_unpack~.c src/mtx_zeros.c \
	$(empty)

BUILT_SOURCES+=iemmatrix_sources.h
CLEANFILES+=iemmatrix_sources.h
iemmatrix_sources.h: $(iemmatrix_la_SOURCES)
	@$(srcdir)/src/makesource.sh $^ > $@

matrix_la_SOURCES	= src/matrix.c
matrix_la_LIBADD	= libiemmatrixCore.la
mtx_abs_la_SOURCES	= src/mtx_abs.c
mtx_abs_la_LIBADD	= libiemmatrixCore.la
mtx_add_la_SOURCES	= src/mtx_add.c
mtx_add_la_LIBADD	= libiemmatrixCore.la
mtx_and_la_SOURCES	= src/mtx_and.c
mtx_and_la_LIBADD	= libiemmatrixCore.la
mtx_atan_la_SOURCES	= src/mtx_atan.c
mtx_atan_la_LIBADD	= libiemmatrixCore.la
mtx_atan2_la_SOURCES	= src/mtx_atan2.c
mtx_atan2_la_LIBADD	= libiemmatrixCore.la
mtx_bessel_la_SOURCES	= src/mtx_bessel.c
mtx_bessel_la_CFLAGS	= $(AM_CFLAGS) $(GSL_CFLAGS)
mtx_bessel_la_LIBADD	= libiemmatrixCore.la $(GSL_LIBS)
mtx_bitand_la_SOURCES	= src/mtx_bitand.c
mtx_bitand_la_LIBADD	= libiemmatrixCore.la
mtx_bitleft_la_SOURCES	= src/mtx_bitleft.c
mtx_bitleft_la_LIBADD	= libiemmatrixCore.la
mtx_bitor_la_SOURCES	= src/mtx_bitor.c
mtx_bitor_la_LIBADD	= libiemmatrixCore.la
mtx_bitright_la_SOURCES	= src/mtx_bitright.c
mtx_bitright_la_LIBADD	= libiemmatrixCore.la
mtx_bspline_la_SOURCES	= src/mtx_bspline.c
mtx_bspline_la_LIBADD	= libiemmatrixCore.la
mtx_check_la_SOURCES	= src/mtx_check.c
mtx_check_la_LIBADD	= libiemmatrixCore.la
mtx_cholesky_la_SOURCES	= src/mtx_cholesky.c
mtx_cholesky_la_LIBADD	= libiemmatrixCore.la
mtx_col_la_SOURCES	= src/mtx_col.c
mtx_col_la_LIBADD	= libiemmatrixCore.la
mtx_colon_la_SOURCES	= src/mtx_colon.c
mtx_colon_la_LIBADD	= libiemmatrixCore.la
mtx_concat_la_SOURCES	= src/mtx_concat.c
mtx_concat_la_LIBADD	= libiemmatrixCore.la
mtx_conv_la_SOURCES	= src/mtx_conv.c
mtx_conv_la_LIBADD	= libiemmatrixCore.la
mtx_convolver__la_CFLAGS=  $(AM_CFLAGS) $(FFTW3F_CFLAGS)
mtx_convolver__la_SOURCES= src/mtx_convolver~.c
mtx_convolver__la_LIBADD = libiemmatrixCore.la libiemmatrixConvolver.la $(FFTW3F_LIBS)
mtx_cos_la_SOURCES	= src/mtx_cos.c
mtx_cos_la_LIBADD	= libiemmatrixCore.la
mtx_cumprod_la_SOURCES	= src/mtx_cumprod.c
mtx_cumprod_la_LIBADD	= libiemmatrixCore.la
mtx_cumsum_la_SOURCES	= src/mtx_cumsum.c
mtx_cumsum_la_LIBADD	= libiemmatrixCore.la
mtx_dbtopow_la_SOURCES	= src/mtx_dbtopow.c
mtx_dbtopow_la_LIBADD	= libiemmatrixCore.la
mtx_dbtorms_la_SOURCES	= src/mtx_dbtorms.c
mtx_dbtorms_la_LIBADD	= libiemmatrixCore.la
mtx_decay_la_SOURCES	= src/mtx_decay.c
mtx_decay_la_LIBADD	= libiemmatrixCore.la
mtx_diag_la_SOURCES	= src/mtx_diag.c
mtx_diag_la_LIBADD	= libiemmatrixCore.la
mtx_diegg_la_SOURCES	= src/mtx_diegg.c
mtx_diegg_la_LIBADD	= libiemmatrixCore.la
mtx_diff_la_SOURCES	= src/mtx_diff.c
mtx_diff_la_LIBADD	= libiemmatrixCore.la
mtx_dispersive_dline_la_SOURCES	= src/mtx_dispersive_dline.c
mtx_dispersive_dline_la_LIBADD	= libiemmatrixCore.la
mtx_distance2_la_SOURCES	= src/mtx_distance2.c
mtx_distance2_la_LIBADD	= libiemmatrixCore.la
mtx_egg_la_SOURCES	= src/mtx_egg.c
mtx_egg_la_LIBADD	= libiemmatrixCore.la
mtx_eig_la_SOURCES	= src/mtx_eig.c
mtx_eig_la_CFLAGS	= $(AM_CFLAGS) $(GSL_CFLAGS)
mtx_eig_la_LIBADD	= libiemmatrixCore.la $(GSL_LIBS)
mtx_ei__la_SOURCES	= src/mtx_ei~.c
mtx_ei__la_LIBADD	= libiemmatrixCore.la
mtx_element_la_SOURCES	= src/mtx_element.c
mtx_element_la_LIBADD	= libiemmatrixCore.la
mtx_eq_la_SOURCES	= src/mtx_eq.c
mtx_eq_la_LIBADD	= libiemmatrixCore.la
mtx_exp_la_SOURCES	= src/mtx_exp.c
mtx_exp_la_LIBADD	= libiemmatrixCore.la
mtx_eye_la_SOURCES	= src/mtx_eye.c
mtx_eye_la_LIBADD	= libiemmatrixCore.la
mtx_fft_la_SOURCES	= src/mtx_fft.c
mtx_fft_la_LIBADD	= libiemmatrixCore.la
mtx_fill_la_SOURCES	= src/mtx_fill.c
mtx_fill_la_LIBADD	= libiemmatrixCore.la
mtx_find_la_SOURCES	= src/mtx_find.c
mtx_find_la_LIBADD	= libiemmatrixCore.la
mtx_gauss_la_SOURCES	= src/mtx_gauss.c
mtx_gauss_la_LIBADD	= libiemmatrixCore.la
mtx_ge_la_SOURCES	= src/mtx_ge.c
mtx_ge_la_LIBADD	= libiemmatrixCore.la
mtx_gt_la_SOURCES	= src/mtx_gt.c
mtx_gt_la_LIBADD	= libiemmatrixCore.la
mtx_ifft_la_SOURCES	= src/mtx_ifft.c
mtx_ifft_la_LIBADD	= libiemmatrixCore.la
mtx_index_la_SOURCES	= src/mtx_index.c
mtx_index_la_LIBADD	= libiemmatrixCore.la
mtx_int_la_SOURCES	= src/mtx_int.c
mtx_int_la_LIBADD	= libiemmatrixCore.la
mtx_inverse_la_SOURCES	= src/mtx_inverse.c
mtx_inverse_la_LIBADD	= libiemmatrixCore.la
mtx_isequal_la_SOURCES	= src/mtx_isequal.c
mtx_isequal_la_LIBADD	= libiemmatrixCore.la
mtx_le_la_SOURCES	= src/mtx_le.c
mtx_le_la_LIBADD	= libiemmatrixCore.la
mtx_log_la_SOURCES	= src/mtx_log.c
mtx_log_la_LIBADD	= libiemmatrixCore.la
mtx_lt_la_SOURCES	= src/mtx_lt.c
mtx_lt_la_LIBADD	= libiemmatrixCore.la
mtx_max2_la_SOURCES	= src/mtx_max2.c
mtx_max2_la_LIBADD	= libiemmatrixCore.la
mtx_mean_la_SOURCES	= src/mtx_mean.c
mtx_mean_la_LIBADD	= libiemmatrixCore.la
mtx_min2_la_SOURCES	= src/mtx_min2.c
mtx_min2_la_LIBADD	= libiemmatrixCore.la
mtx_minmax_la_SOURCES	= src/mtx_minmax.c
mtx_minmax_la_LIBADD	= libiemmatrixCore.la
mtx_mul_la_SOURCES	= src/mtx_mul.c
mtx_mul_la_LIBADD	= libiemmatrixCore.la
mtx_mul__la_SOURCES	= src/mtx_mul~.c
mtx_mul__la_LIBADD	= libiemmatrixCore.la
mtx_neq_la_SOURCES	= src/mtx_neq.c
mtx_neq_la_LIBADD	= libiemmatrixCore.la
mtx_not_la_SOURCES	= src/mtx_not.c
mtx_not_la_LIBADD	= libiemmatrixCore.la
mtx_ones_la_SOURCES	= src/mtx_ones.c
mtx_ones_la_LIBADD	= libiemmatrixCore.la
mtx_or_la_SOURCES	= src/mtx_or.c
mtx_or_la_LIBADD	= libiemmatrixCore.la
mtx_pack__la_SOURCES	= src/mtx_pack~.c
mtx_pack__la_LIBADD	= libiemmatrixCore.la
mtx_pivot_la_SOURCES	= src/mtx_pivot.c
mtx_pivot_la_LIBADD	= libiemmatrixCore.la
mtx_pow_la_SOURCES	= src/mtx_pow.c
mtx_pow_la_LIBADD	= libiemmatrixCore.la
mtx_powtodb_la_SOURCES	= src/mtx_powtodb.c
mtx_powtodb_la_LIBADD	= libiemmatrixCore.la
mtx_print_la_SOURCES	= src/mtx_print.c
mtx_print_la_LIBADD	= libiemmatrixCore.la
mtx_prod_la_SOURCES	= src/mtx_prod.c
mtx_prod_la_LIBADD	= libiemmatrixCore.la
mtx_qhull_la_SOURCES	= src/mtx_qhull.c
mtx_qhull_la_LIBADD	= libiemmatrixCore.la libiemmatrixQhull.la
mtx_qr_la_SOURCES	= src/mtx_qr.c
mtx_qr_la_CFLAGS	= $(AM_CFLAGS) $(GSL_CFLAGS)
mtx_qr_la_LIBADD	= libiemmatrixCore.la $(GSL_LIBS)
mtx_rand_la_SOURCES	= src/mtx_rand.c
mtx_rand_la_LIBADD	= libiemmatrixCore.la
mtx_repmat_la_SOURCES	= src/mtx_repmat.c
mtx_repmat_la_LIBADD	= libiemmatrixCore.la
mtx_resize_la_SOURCES	= src/mtx_resize.c
mtx_resize_la_LIBADD	= libiemmatrixCore.la
mtx_reverse_la_SOURCES	= src/mtx_reverse.c
mtx_reverse_la_LIBADD	= libiemmatrixCore.la
mtx_rfft_la_SOURCES	= src/mtx_rfft.c
mtx_rfft_la_CFLAGS      = $(AM_CFLAGS) $(FFTW3_CFLAGS)
mtx_rfft_la_LIBADD	= libiemmatrixCore.la $(FFTW3_LIBS)
mtx_rifft_la_SOURCES	= src/mtx_rifft.c
mtx_rifft_la_CFLAGS     = $(AM_CFLAGS) $(FFTW3_CFLAGS)
mtx_rifft_la_LIBADD	= libiemmatrixCore.la $(FFTW3_LIBS)
mtx_rmstodb_la_SOURCES	= src/mtx_rmstodb.c
mtx_rmstodb_la_LIBADD	= libiemmatrixCore.la
mtx_roll_la_SOURCES	= src/mtx_roll.c
mtx_roll_la_LIBADD	= libiemmatrixCore.la
mtx_row_la_SOURCES	= src/mtx_row.c
mtx_row_la_LIBADD	= libiemmatrixCore.la
mtx_scroll_la_SOURCES	= src/mtx_scroll.c
mtx_scroll_la_LIBADD	= libiemmatrixCore.la
mtx_sin_la_SOURCES	= src/mtx_sin.c
mtx_sin_la_LIBADD	= libiemmatrixCore.la
mtx_size_la_SOURCES	= src/mtx_size.c
mtx_size_la_LIBADD	= libiemmatrixCore.la
mtx_slice_la_SOURCES	= src/mtx_slice.c
mtx_slice_la_LIBADD	= libiemmatrixCore.la
mtx_sndfileread_la_SOURCES	= src/mtx_sndfileread.c
mtx_sndfileread_la_CFLAGS	= $(AM_CFLAGS) $(SNDFILE_CFLAGS)
mtx_sndfileread_la_LIBADD	= libiemmatrixCore.la $(SNDFILE_LIBS)
mtx_sort_la_SOURCES	= src/mtx_sort.c
mtx_sort_la_LIBADD	= libiemmatrixCore.la
mtx_spherical_harmonics_la_SOURCES	= src/mtx_spherical_harmonics.c
mtx_spherical_harmonics_la_LIBADD	= libiemmatrixCore.la libiemmatrixSphericalHarmonics.la
mtx_spherical_radial_la_SOURCES	= src/mtx_spherical_radial.c
mtx_spherical_radial_la_LIBADD	= libiemmatrixCore.la libiemmatrixSphericalRadial.la
mtx_sub_la_SOURCES	= src/mtx_sub.c
mtx_sub_la_LIBADD	= libiemmatrixCore.la
mtx_sum_la_SOURCES	= src/mtx_sum.c
mtx_sum_la_LIBADD	= libiemmatrixCore.la
mtx_svd_la_SOURCES	= src/mtx_svd.c
mtx_svd_la_CFLAGS	= $(AM_CFLAGS) $(GSL_CFLAGS)
mtx_svd_la_LIBADD	= libiemmatrixCore.la $(GSL_LIBS)
mtx_trace_la_SOURCES	= src/mtx_trace.c
mtx_trace_la_LIBADD	= libiemmatrixCore.la
mtx_transpose_la_SOURCES= src/mtx_transpose.c
mtx_transpose_la_LIBADD	= libiemmatrixCore.la
mtx_unpack__la_SOURCES	= src/mtx_unpack~.c
mtx_unpack__la_LIBADD	= libiemmatrixCore.la
mtx_zeros_la_SOURCES	= src/mtx_zeros.c
mtx_zeros_la_LIBADD	= libiemmatrixCore.la
matrix_mul_line__la_SOURCES	= alias/matrix_mul_line~.c
matrix_mul_line__la_LIBADD	= libiemmatrixCore.la
matrix_mul__la_SOURCES	= alias/matrix_mul~.c
matrix_mul__la_LIBADD	= libiemmatrixCore.la
matrix__la_SOURCES	= alias/matrix~.c
matrix__la_LIBADD	= libiemmatrixCore.la
mtx_la_SOURCES		= alias/mtx.c
mtx_la_LIBADD		= libiemmatrixCore.la
mtx_0x21_la_SOURCES	= alias/mtx_0x21.c
mtx_0x21_la_LIBADD	= libiemmatrixCore.la
mtx_0x210x3d_la_SOURCES	= alias/mtx_0x210x3d.c
mtx_0x210x3d_la_LIBADD	= libiemmatrixCore.la
mtx_0x26_la_SOURCES	= alias/mtx_0x26.c
mtx_0x26_la_LIBADD	= libiemmatrixCore.la
mtx_0x260x26_la_SOURCES	= alias/mtx_0x260x26.c
mtx_0x260x26_la_LIBADD	= libiemmatrixCore.la
mtx_0x2a_la_SOURCES	= alias/mtx_0x2a.c
mtx_0x2a_la_LIBADD	= libiemmatrixCore.la
mtx_0x2a0x7e_la_SOURCES	= alias/mtx_0x2a0x7e.c
mtx_0x2a0x7e_la_LIBADD	= libiemmatrixCore.la
mtx_0x2b_la_SOURCES	= alias/mtx_0x2b.c
mtx_0x2b_la_LIBADD	= libiemmatrixCore.la
mtx_0x2d_la_SOURCES	= alias/mtx_0x2d.c
mtx_0x2d_la_LIBADD	= libiemmatrixCore.la
mtx_0x2e0x2a_la_SOURCES	= alias/mtx_0x2e0x2a.c
mtx_0x2e0x2a_la_LIBADD	= libiemmatrixCore.la
mtx_0x2e0x2f_la_SOURCES	= alias/mtx_0x2e0x2f.c
mtx_0x2e0x2f_la_LIBADD	= libiemmatrixCore.la
mtx_0x2e0x5e_la_SOURCES	= alias/mtx_0x2e0x5e.c
mtx_0x2e0x5e_la_LIBADD	= libiemmatrixCore.la
mtx_0x2f_la_SOURCES	= alias/mtx_0x2f.c
mtx_0x2f_la_LIBADD	= libiemmatrixCore.la
mtx_0x3a_la_SOURCES	= alias/mtx_0x3a.c
mtx_0x3a_la_LIBADD	= libiemmatrixCore.la
mtx_0x3c_la_SOURCES	= alias/mtx_0x3c.c
mtx_0x3c_la_LIBADD	= libiemmatrixCore.la
mtx_0x3c0x3c_la_SOURCES	= alias/mtx_0x3c0x3c.c
mtx_0x3c0x3c_la_LIBADD	= libiemmatrixCore.la
mtx_0x3c0x3d_la_SOURCES	= alias/mtx_0x3c0x3d.c
mtx_0x3c0x3d_la_LIBADD	= libiemmatrixCore.la
mtx_0x3d0x3d_la_SOURCES	= alias/mtx_0x3d0x3d.c
mtx_0x3d0x3d_la_LIBADD	= libiemmatrixCore.la
mtx_0x3e_la_SOURCES	= alias/mtx_0x3e.c
mtx_0x3e_la_LIBADD	= libiemmatrixCore.la
mtx_0x3e0x3d_la_SOURCES	= alias/mtx_0x3e0x3d.c
mtx_0x3e0x3d_la_LIBADD	= libiemmatrixCore.la
mtx_0x3e0x3e_la_SOURCES	= alias/mtx_0x3e0x3e.c
mtx_0x3e0x3e_la_LIBADD	= libiemmatrixCore.la
mtx_0x7c_la_SOURCES	= alias/mtx_0x7c.c
mtx_0x7c_la_LIBADD	= libiemmatrixCore.la
mtx_0x7c0x7c_la_SOURCES	= alias/mtx_0x7c0x7c.c
mtx_0x7c0x7c_la_LIBADD	= libiemmatrixCore.la
mtx_div_la_SOURCES	= alias/mtx_div.c
mtx_div_la_LIBADD	= libiemmatrixCore.la
