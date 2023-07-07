dnl Process this file with autoconf to produce a configure script.
AC_PREREQ([2.60])

AC_INIT([iemmatrix], [0.3.2], [zmoelnig@iem.at], [iemmatrix], [https://git.iem.at/pd/iemmatrix])
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
KERN=$(uname -s)

dnl m4_ifdef([AC_LIB_RPATH], [AC_LIB_RPATH])

dnl Checks for programs.
AC_PROG_CC
AC_PROG_LN_S
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL


dnl Checks for libraries.
dnl Replace `main' with a function in -lc:
AC_CHECK_LIB([c], [main])

dnl Replace `main' with a function in -lm:
AC_CHECK_LIB([m], [cos])
dnl for math.h Bessel/Neumann functions
AC_CHECK_LIB([m], [jn])
AS_IF([test "xyes" = "x${ac_cv_lib_m_jn}"],[
   AC_DEFINE(HAVE_MATH_BESSEL,1,[can math compute Bessel functions?])
])

dnl FFT
AC_ARG_WITH([fftw], AC_HELP_STRING([--without-fftw], [compile without FFTW]))
AS_IF([test "xno" != "x${with_fftw}" ],[
   PKG_CHECK_MODULES([FFTW3], [fftw3], [have_fftw3=yes], [have_fftw3=no])
])
AS_IF([test "xyes" = "x${have_fftw3}" ],[
   AC_DEFINE(USE_FFTW, 1, [can we use fftw?])
])

dnl for soundfile reading (and writing, if we do that...)
AC_ARG_WITH([sndfile], AC_HELP_STRING([--without-sndfile], [compile without libsndfile support]))
have_sndfile="no"
AS_IF([test "xno" != "x${with_sndfile}"],[
   PKG_CHECK_MODULES([SNDFILE], [sndfile], [have_sndfile=yes], [have_sndfile=no])
])
AS_IF([test "xyes" = "x${have_sndfile}"],[
   AC_DEFINE([HAVE_LIBSNDFILE],1,[Define to 1 if you have the 'sndfile' library (-lsndfile).])
   AC_DEFINE([HAVE_SNDFILE_H],1,[Define to 1 if you have the <sndfile.h> header file.])
])


dnl for gnu scientific library -lgsl:
AC_ARG_WITH([gsl], AC_HELP_STRING([--without-gsl], [compile without GNU Scientific Library]))
AS_IF([test "xno" != "x${with_gsl}"],[
   AC_CHECK_LIB([gslcblas], [cblas_dgemm])
   AC_CHECK_LIB([gsl], [gsl_blas_dgemm])
   AS_IF([test "x${ac_cv_lib_gsl_gsl_blas_dgemm}" = "xyes"],[
        AC_CHECK_LIB([gsl], [gsl_eigen_nonsymm], [:])
        AC_CHECK_LIB([gsl], [gsl_sf_bessel_Jn], [:])
   ])
])
AS_IF([test "x${ac_cv_lib_gsl_gsl_blas_dgemm}" = "xyes"],[
   AC_DEFINE(HAVE_LIBGSL,1,[Define to 1 if you have the 'gsl' library (-lgsl).])
])
AS_IF([test "x${ac_cv_lib_gslcblas_cblas_dgemm}" = "xyes"],[
   AC_DEFINE(HAVE_LIBGSLCBLAS,1,[Define to 1 if you have the 'gslcblas' library (-lgslcblas).])
])
AS_IF([test "x${ac_cv_lib_gsl_gsl_eigen_nonsymm}" = "xyes"],[
   AC_DEFINE(HAVE_GSL_EIGEN_NONSYMM,1,[can gsl compute complex eigenvalues?])
])

AS_IF([test "x${ac_cv_lib_gsl_gsl_sf_bessel_Jn}" = "xyes"],[
   AC_DEFINE(HAVE_GSL_BESSEL,1,[can gsl compute Bessel functions?])
])


###############################################

AS_IF([test "x${libdir}" = "x\${exec_prefix}/lib"],
 [libdir='${exec_prefix}/lib/pd/extra']
)


# handle Pd path
AC_ARG_WITH([pd],
  AC_HELP_STRING([--with-pd=</path/to/pd>],
                 [where to look for Pd-headers and and -libs]))
AS_IF([test -d "${with_pd}"], [PDPATH=${with_pd}])
AS_IF([test -d "${PDPATH}"], [
 AS_IF([test -d "${PDPATH}/src"], [
   AC_MSG_RESULT([adding ${PDPATH}/src to CPPFLAGS])
   CPPFLAGS="-I${PDPATH}/src ${CPPFLAGS}"

   AC_MSG_RESULT([adding ${PDPATH}/src to LDFLAGS])
   LDFLAGS="-L${PDPATH}/src ${LDFLAGS}"
 ],[
   AC_MSG_RESULT([adding ${PDPATH} to CPPFLAGS])
   CPPFLAGS="-I${PDPATH} ${CPPFLAGS}"
 ])
 AS_IF([test -d "${PDPATH}/bin"], [
   AC_MSG_RESULT([adding ${PDPATH}/bin to LDFLAGS])
   LDFLAGS="-L${PDPATH}/bin ${LDFLAGS}"
 ],[
   AC_MSG_RESULT([adding ${PDPATH} to LDFLAGS])
   LDFLAGS="-L${PDPATH} ${LDFLAGS}"
 ])
])

AC_ARG_ENABLE([library],
  AC_HELP_STRING([--enable-library],
                 [build iemmatrix as multi-object library (default, "yes"), as single-object externals ("no") or as both ("both")]))
AM_CONDITIONAL([LIBRARY], [test "x${enable_library}" != "xno"])
AM_CONDITIONAL([SINGLEOBJECTS], [test "x${enable_library}" = "xno" || test "x${enable_library}" = "xboth"])

AC_ARG_WITH([extension],
  AC_HELP_STRING([--with-extension=<ext>],
                 [enforce a certain extension for the dynamic library (e.g. dll)]))
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


tmp_pddll_LIBS=$LIBS
AC_CHECK_LIB([:pd.dll], [nullfn], [have_pddll="yes"], [have_pddll="no"])
LIBS=$tmp_pddll_LIBS

AS_IF([test "x$have_pddll" = "xyes"],[
  EXTRA_LTFLAGS="$EXTRA_LTFLAGS -Xlinker -l:pd.dll"
  ],[
  AC_CHECK_LIB([pd], [nullfn])
])

IEM_CHECK_FAT(ARCH_FLAG)

IEM_CHECK_CFLAGS([-mms-bitfields])
IEM_CHECK_SIMD

have_pd="no"
AC_CHECK_HEADERS([pd/m_pd.h m_pd.h], [have_pd="yes"])
AS_IF([test "x${have_pd}" = "xno"],[
AC_ERROR([m_pd.h is desperately needed!
        install pd and/or use
        "--with-pd=</path/to/pd/>"])
])

dnl Reproducible builds
AC_ARG_WITH([build-date],
  AC_HELP_STRING([--with-build-date=<str>],
                 [alternative date string]))
AC_MSG_CHECKING([build-date])
AS_IF([test "x$with_build_date" = "xyes" ],[with_build_date=""])
AS_IF([test "x$with_build_date" = "xno"  ],[with_build_date=""])
AS_IF([test "x$with_build_date" != "x"  ],[ AC_DEFINE_UNQUOTED([BUILD_DATE], ["$with_build_date"], [user-defined build-date string])])
## urgh, brackets need to be escaped with @<:@ and @:>@...
AS_IF([test "x$with_build_date" != "x"  ],[ with_build_date="@<:@default@:>@" ])
AC_MSG_RESULT([$with_build_date])

AC_SUBST([CPPFLAGS])
AC_SUBST([CFLAGS])
AC_SUBST([LIBS])
AC_SUBST([DEFS])
AC_SUBST([ARCH_FLAG])
AC_SUBST([LIBTOOL_DEPS])
AC_SUBST([EXTRA_LTFLAGS])
AC_SUBST([EXTENSION])

AC_OUTPUT()

