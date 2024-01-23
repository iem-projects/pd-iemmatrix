iemmatrix - matrix objects for Pure Data
========================================

Homepage: https://git.iem.at/pd/iemmatrix



Installation/Compilation instructions for "iemmatrix"
=====================================================

iemmatrix uses the *autoconf* build system.
For generic build instructions, please see INSTALL.txt

Dependencies
============

in order to compile iemmatrix, you will need some libraries/applications installed.
Some of these dependencies are optional.
If they are present during compilation, iemmatrix will have certain functionality
enabled, which might be missing otherwise.

 - puredata
     http://puredata.info/downloads
     make sure you also have the development files (headers)

 - GNU scientific library (aka gsl) [optional]
     http://www.gnu.org/software/gsl
     needed for higher maths, including eigenvalues and singular-value decomposition
 - FFTW (Fast Fourier Transform) [optional]
     http://www.fftw.org
     fast, high-precision FFTs; if not present, iemmatrix will use Pd's internal FFT
 - sndfile (reading/writing audio files) [optional]
     http://www.mega-nerd.com/libsndfile/
     for reading soundfiles into matrices (many soundfiles are supported)


 - autotools, libtool, gettext
     that's the basic build-system

 - a compiler, linker,...
     tested with
      - GCC (the GNU compiler collection)
      - Clang

Building
========

Make sure that each step succeeds, before proceeding to the next one!

 $ ./autogen.sh
 $ ./configure
 $ make


./configure will need to find the Pd-headers (and on some platforms the Pd-library).
If these files are in a non-standard location (e.g. on W32 and OSX),
you have to manually tell it where to look for them, using the `--with-pd` flag
and pointing it to the directory that contains the 'bin/' and 'src/' (or 'include/')
folders of your Pd-distribution.

E.g.

 $ ./configure --with-pd=/Applications/Pd-0.45-4.app/Contents/Resources

Once `make` has succeeded, you will have a `iemmatrix` binary in the (hidden)
`.libs` folder.


Installing
==========

The ordinary way to install, is by running the following with the proper privileges
(e.g. as root):

 # make install

This will install the entire iemmatrix into /usr/local/lib/pd/extra/iemmatrix.

On systems that have no standard filesystem layout for Pd-externals (e.g. W32 and OSX),
this is not exactly what you want.
Instead, you can use the following to collect all installation data into a single directory:

 $ make install DESTDIR=$(pwd)/ pkglibdir=iemmatrix

This will create a new directory `iemmatrix` (in your current directory),
containing all binaries and abstractions needed.

You can then take this directory, and put it into a place, where Pd will look for it:

E.g.

 - Linux: `~/pd-externals` (e.g. `/home/frodo/pd-externals`)
 - OSX  : `~/Library/Pd`   (e.g. `/Users/frodo/Pd`)
 - W32  : `%AppData%\Pd`   (e.g. `C:\Documents and Settings\frodo\Application Data\Pd`)

A full list of default search paths for externals, can be found at
  https://puredata.info/docs/faq/how-do-i-install-externals-and-help-files/


System-specific instructions
============================

Linux (Debian-based)
--------------------

`iemmatrix` is available as Debian package

  # aptitude install pd-iemmatrix

To install all dependencies for compiling iemmatrix yourself, use:

 # apt-get build-dep pd-iemmatrix

or

 # aptitude install automake1.11 autoconf puredata libfftw3-dev libsndfile1-dev libgsl0-dev


mac OS-X
--------

Most dependencies are available via `brew`

   http://brew.sh

Install them using

   $ brew install gsl --universal
   $ brew install fftw --universal
   $ brew install libsndfile --universal

When running iemmatrix' configure, make sure to specify the path to Pd.
If you want to build universal binaries (e.g. both 32bit (i386) and 64bit (x86_64)),
you can specify the `--enable-fat-binary` flag.
E.g.

 $ ./configure  --with-pd=/Applications/Pd-0.45-4.app/Contents/Resources --enable-fat-binary=i386,x86_64

Note: when building fat binaries, all dependencies must be fat too.

Microsoft W32
-------------

The only sane way to build iemmatrix for W32 using autotools,
is currently by using `MinGW-w64` as a cross-compilation environment on Linux.

Debian (and derivates) provides packages for this:

 # aptitude install mingw-w64 mingw-w64-i686-dev mingw-w64-tools binutils-mingw-w64

Run "./configure" and specify the path to your W32 installation of Pd
(where you unzipped the W32 package of Pd) using the `--with-pd` flag.
Make sure that the Pd-sources are in PDPATH/src and the compiled pd-binaries in PDPATH/bin.
If they are scattered across you filesystem you can alternatively give explicitly the paths
to your "m_pd.h" (with `--includedir`) and to your "pd.lib" (with `--libdir`).
Don't forget to override the default extension ("pd_linux" on Linux-systems)
with the more appropriate "dll".
On bash this looks like:

 $ ./configure --with-extension=dll --host=i686-w64-mingw32 --with-pd=/home/frodo/W32/pd-0.46-6

Now run "make" and enjoy.
	

W32 using VisualStudio
----------------------
There are outdated VisualStudio project files in src/
Use them at your own risk.
