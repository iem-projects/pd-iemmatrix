dnl Process this file with autoconf to produce a configure script.
AC_PREREQ([2.60])

AC_INIT([iemmatrix], [0.2], [zmoelnig@iem.at], [iemmatrix], [http://puredata.info/downloads/iemmatrix])
AM_INIT_AUTOMAKE([1.10 foreign])
m4_ifdef([LT_INIT], 
	[LT_INIT([disable-static win32-dll])],
	[AC_ENABLE_SHARED
	 AC_DISABLE_STATIC
	 AC_LIBTOOL_WIN32_DLL])

AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([doc/Makefile abs/Makefile])
AC_CONFIG_FILES([iemmatrix-meta.pd])
AC_CONFIG_FILES([VERSION.txt])

AC_CONFIG_SRCDIR([src/iemmatrix.c])

IEM_OPERATING_SYSTEM
ARCH=$(uname -m)
KERN=$(uname -s)

dnl m4_ifdef([AC_LIB_RPATH], [AC_LIB_RPATH])

dnl Checks for programs.
AC_PROG_CC
AC_PROG_LN_S
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL

AC_ARG_WITH([pd],       [  --with-pd=</path/to/pd> where to look for pd-headers and and -libs])
AC_ARG_WITH([extension],[  --with-extension=<ext>  enforce a certain extension for the dynamic library (e.g. dll)])
AC_ARG_ENABLE([library],[  --enable-library        build iemmatrix as multi-object library (default, "yes"), as single-object externals ("no") or as both ("both")])
AM_CONDITIONAL([LIBRARY], [test "x${enable_library}" != "xno"])
AM_CONDITIONAL([SINGLEOBJECTS], [test "x${enable_library}" = "xno" || test "x${enable_library}" = "xboth"])


if test "x${libdir}" = "x\${exec_prefix}/lib"; then
 libdir='${exec_prefix}/lib/pd/extra'
fi


if test "x$with_pd" != "x"; then
 if test -d "${with_pd}"; then
   PDPATH=${with_pd}
 fi
 if test -d "${PDPATH}/src"; then
  AC_MSG_RESULT([adding ${PDPATH}/src to CPPFLAGS])
  CPPFLAGS="-I${PDPATH}/src ${CPPFLAGS}"

  AC_MSG_RESULT([adding ${PDPATH}/src to LDFLAGS])
  LDFLAGS="-L${PDPATH}/src ${LDFLAGS}"
 else
  if test -d "${PDPATH}"; then
   AC_MSG_RESULT([adding ${PDPATH} to CPPFLAGS])
   CPPFLAGS="-I${PDPATH} ${CPPFLAGS}"
  fi
 fi
 if test -d "${PDPATH}/bin"; then
  AC_MSG_RESULT([adding ${PDPATH}/bin to LDFLAGS])
  LDFLAGS="-L${PDPATH}/bin ${LDFLAGS}"
 else
  if test -d "${PDPATH}"; then
   AC_MSG_RESULT([adding ${PDPATH} to LDFLAGS])
   LDFLAGS="-L${PDPATH} ${LDFLAGS}"
  fi
 fi
fi

dnl Checks for libraries.
dnl Replace `main' with a function in -lc:
AC_CHECK_LIB([c], [main])
AC_CHECK_LIB([crtdll], [fclose])

AC_CHECK_LIB([coldname], [lseek])
AC_CHECK_LIB([kernel32], [main])
AC_CHECK_LIB([user32], [main])
AC_CHECK_LIB([gdi32], [main])

dnl Replace `main' with a function in -lm:
AC_CHECK_LIB([m], [main])
dnl Replace `main' with a function in -lpthread:
dnl AC_CHECK_LIB([pthread], [main])
dnl Replace `main' with a function in -lstk:
dnl AC_CHECK_LIB([stk], [main], [STK=yes])

tmp_pddll_LIBS=$LIBS
AC_CHECK_LIB([:pd.dll], [nullfn], have_pddll="yes", have_pddll="no")
LIBS=$tmp_pddll_LIBS

if test "x$have_pddll" = "xyes"
then
  EXTRA_LTFLAGS="$EXTRA_LTFLAGS -Xlinker -l:pd.dll"
else
  AC_CHECK_LIB([pd], [nullfn])
fi

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([io.h math.h memory.h stdio.h stdlib.h string.h unistd.h])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_TIME

dnl Checks for library functions.
AC_FUNC_MMAP

IEM_CHECK_FAT(ARCH_FLAG)

IEM_CHECK_CFLAGS([-mms-bitfields])
IEM_CHECK_SIMD

AC_MSG_CHECKING([extension])
if test "x$with_extension" != "x"
then
 EXTENSION=$with_extension
else
 case "$KERN" in
 Darwin)
   EXTENSION=pd_darwin
   ;;
 Linux)
   EXTENSION=pd_linux
   ;;
 IRIX64)
   EXTENSION=pd_irix6
   ;;
 IRIX32)
   EXTENSION=pd_irix5
   ;;
 IRIX64)
   EXTENSION=pd_irix6
   ;;
 *CYGWIN*)
   EXTENSION=dll
   ;;
 *MINGW*)
   EXTENSION=dll
   ;;
 *)
   EXTENSION=pd_linux
   ;;
 esac
fi
AC_MSG_RESULT([$EXTENSION])


have_pd="no"
AC_CHECK_HEADERS([pd/m_pd.h m_pd.h], [have_pd="yes"])
if test "x${have_pd}" = "xno"
then
AC_ERROR([m_pd.h is desperately needed!
        install pd and/or use
        "--with-pd=</path/to/pd/>"])
fi

AC_ARG_WITH([fftw], AC_HELP_STRING([--without-fftw], [compile without FFTW]))
if test "xno" != "x${with_fftw}" ; then
have_fftw="no"
have_fftw3="no"
AC_CHECK_LIB([fftw3],
	     [fftw_execute],
	     [AC_CHECK_HEADER([fftw3.h],
		[have_fftw3="yes"],
		)])
fi
if test "xyes" = "x${have_fftw3}"; then
   have_fftw="yes"
   AC_DEFINE(HAVE_LIBFFTW3, 1, [Define to 1 if you have the `fftw3' library (-lfftw3).])
fi
if test "xyes" = "x${have_fftw}"; then
dnl   AC_DEFINE(HAVE_LIBFFTW,1,[Define to 1 if you have the `fftw' library (e.g. -lfftw3).])
   AC_DEFINE(USE_FFTW, 1, [can we use fftw?])
fi

AM_CONDITIONAL(HAVE_FFTW, [test "x${have_fftw}" = "xyes"])

dnl for soundfile reading (and writing, if we do that...)
AC_ARG_WITH([sndfile], AC_HELP_STRING([--without-sndfile], [compile without libsndfile support]))
have_sndfile="no"
if test "xno" != "x${with_sndfile}" ; then
   AC_CHECK_LIB([sndfile], [sf_close], [have_sndfile="yes"])
fi
if test "xyes" = "x${have_sndfile}"; then
   AC_DEFINE(HAVE_LIBSNDFILE,1,[Define to 1 if you have the `sndfile' library (-lsndfile).])
fi
AM_CONDITIONAL(HAVE_SNDFILE, [test "x${have_sndfile}" = "xyes"])



dnl for gnu scientific library -lgsl:
AC_ARG_WITH([gsl], AC_HELP_STRING([--without-gsl], [compile without GNU Scientific Library]))
have_gsl="no"
have_gsl_cblas="no"
have_gsl_eigen="no"
have_gsl_bessel="no"
if test "xno" != "x${with_gsl}" ; then
   AC_CHECK_LIB([gslcblas], [cblas_dgemm],[have_gsl_cblas="yes"])
   AC_CHECK_LIB([gsl], [gsl_blas_dgemm],[have_gsl="yes"])
   AC_CHECK_LIB([gsl], [gsl_eigen_nonsymm],[have_gsl_eigen="yes"])
   AC_CHECK_LIB([gsl], [gsl_sf_bessel_Jn],[have_gsl_eigen="yes"])
fi
AM_CONDITIONAL(HAVE_GSL, [test "x${have_gsl}" = "xyes"])
AM_CONDITIONAL(HAVE_GSL_CBLAS, [test "x${have_gsl_cblas}" = "xyes"])
if test "xyes" = "x${have_gsl}"; then
   AC_DEFINE(HAVE_LIBGSL,1,[Define to 1 if you have the `gsl' library (-lgsl).])
fi
if test "xyes" = "x${have_gsl_cblas}"; then
   AC_DEFINE(HAVE_LIBGSLCBLAS,1,[Define to 1 if you have the `gslcblas' library (-lgslcblas).])
fi
if test "xyes" = "x${have_gsl_eigen}"; then
   AC_DEFINE(HAVE_GSL_EIGEN_NONSYMM,1,[can gsl compute complex eigenvalues?])
fi
AM_CONDITIONAL(HAVE_GSL_EIGEN, [test "x${have_gsl_eigen}" = "xyes"])
if test "xyes" = "x${have_gsl_bessel}"; then
   AC_DEFINE(HAVE_GSL_BESSEL,1,[can gsl compute Bessel functions?])
fi
AM_CONDITIONAL(HAVE_GSL_BESSEL, [test "x${have_gsl_bessel}" = "xyes"])


dnl for math.h Bessel/Neumann functions
have_math_bessel="no"
AC_CHECK_LIB([m], [jn],have_math_bessel="yes")
if test "xyes" = "x${have_math_bessel}"; then
   AC_DEFINE(HAVE_MATH_BESSEL,1,[can math compute Bessel functions?])
fi
AM_CONDITIONAL(HAVE_MATH_BESSEL, [test "x${have_math_bessel}" = "xyes"])



AC_SUBST([CPPFLAGS])
AC_SUBST([CFLAGS])
AC_SUBST([DEFS])
AC_SUBST([ARCH_FLAG])
AC_SUBST([LIBTOOL_DEPS])
AC_SUBST([EXTRA_LTFLAGS])

AC_SUBST(EXTENSION)

AC_OUTPUT()

