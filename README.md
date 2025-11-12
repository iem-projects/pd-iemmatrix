iemmatrix - matrix objects for Pure Data
========================================

Homepage: https://git.iem.at/pd/iemmatrix

`iemmatrix` is a collection of objects for Pure Data (Pd)
that allow manipulation of simple matrices.


## Found a bug? Miss a feature?

Please use our issue-tracker at
https://git.iem.at/pd/iemmatrix/-/issues

## Installation instructions

### Simple installation

`iemmatrix` is available via deken (Pd's built in package manager):

For Pd>=0.56
- Open the `Tools` menu
- Select `Find externals...`
- Type **`iemmatrix`** and hit `Search`
- Select the first entry and click on `Install`

For older Pd versions, `Find externals...` can be found in the `Help` menu.


### Compiling iemmatrix

If deken does not offer any downloads for your system
(e.g. because you are stuck on an old OSX system with a PowerPC processor,
or you want to use a not-yet released version),
you can locally compile `iemmatrix` instead.

#### Dependencies

in order to compile iemmatrix, you will need some libraries/applications installed.
Some of these dependencies are optional.
If they are present during compilation, iemmatrix will have certain functionality
enabled, which might be missing otherwise.

 - Pure Data
     https://puredata.info/downloads
     make sure you also have the development files (headers)

 - GNU scientific library (aka gsl) [optional]
     https://www.gnu.org/software/gsl/
     needed for higher maths, including eigenvalues and singular-value decomposition
 - FFTW (Fast Fourier Transform) [optional]
     https://fftw.org
     fast, high-precision FFTs; if not present, iemmatrix will use Pd's internal FFT
 - sndfile (reading/writing audio files) [optional]
     https://libsndfile.github.io/libsndfile/
     for reading soundfiles into matrices (many soundfiles are supported)

  - the `pkg-config` tool (to detect) the optional dependencies

 - a compiler, linker,...
     tested with
      - GCC (the GNU compiler collection)
      - Clang

#### Building

iemmatrix uses the [pd-lib-builder](https://github.com/pure-data/pd-lib-builder) build system.

just run:

```sh
make
```

To get some basic help with the build-system run

```sh
make help
```

or read the [online documentation for pd-lib-builder](https://github.com/pure-data/pd-lib-builder).

#### Installing

The ordinary way to install, is by running the following with the proper privileges
(e.g. as root):

```sh
make install
```

This will install the entire iemmatrix into `/usr/local/lib/pd-externals/iemmatrix/`

On systems that have no standard filesystem layout for Pd-externals (e.g. W32 and macOS),
this is not exactly what you want.
Instead, you can use the following to collect all installation data into a single directory:

```sh
make install DESTDIR=$(pwd) pkglibdir=
```

This will create a new directory `iemmatrix` (in your current directory),
containing all binaries and abstractions needed.

You can then take this directory, and put it into a place, where Pd will look for it:

E.g.

| OS      | path                   | example
|---------|------------------------|------------------------------------------------------
| Linux   |`~/.local/lib/pd/extra` | `/home/frodo/.local/lib/pd/extra/iemmatrix`
| macOS   | `~/Documents/Pd/extra` | `/Users/frodo/Documents/Pd/extra/iemmatrix`
| Windows | `%AppData%\Pd`         | `C:\Users\frodo\AppData\Roaming\Pd\iemmatrix`

A full list of default search paths for externals, can be found at
  https://puredata.info/docs/faq/how-do-i-install-externals-and-help-files/


## System-specific instructions

### Linux (Debian-based)

`iemmatrix` is available as a Debian package

```sh
apt-get install pd-iemmatrix
```

To get all the build-dependencies, use:

```sh
apt-get build-dep pd-iemmatrix
```

or manually with:
```sh
apt-get install build-essential puredata libfftw3-dev libsndfile1-dev libgsl0-dev
```

### macOS


The additional dependencies are available via [`brew`](https://brew.sh)

```
brew install pkgconf gsl fftw libsndfile
```

Then proceed with the [building instructions](#building)

### Windows


You will need the [MSYS2](https://www.msys2.org/) environment with [MinGW](https://www.mingw-w64.org/) installed.

Open an MinGW shell (`MSYS2 MinGW 64-bit`; for older 32bit Windows systems use `MSYS Mingw 32-bit`)
and install the dependencies with:

```sh
pacman -Suy \
  ${MINGW_PACKAGE_PREFIX}-gcc \
  ${MINGW_PACKAGE_PREFIX}-pkgconf \
  ${MINGW_PACKAGE_PREFIX}-ntldd \
  ${MINGW_PACKAGE_PREFIX}-libsndfile \
  ${MINGW_PACKAGE_PREFIX}-fftw \
  ${MINGW_PACKAGE_PREFIX}-gsl
```

Then proceed with the [building instructions](#building)


## License
This program is free software; you can redistribute it and/or
modify it under the terms of the **GNU General Public License**
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Authors

`iemmatrix` is being developped at the
[Institute of Electronic Music and Acoustics](https://iem.at),
which is part of the [University of Music and Performing Arts, Graz/Austria](https://www.kug.ac.at).

For a list of developers, see the `AUTHORS.txt` file.
