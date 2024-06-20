# Makefile to build class 'mtx_convolver~' for Pure Data.
# Needs Makefile.pdlibbuilder as helper makefile for platform-dependent build
# settings and rules.

# library name
lib.name = mtx_convolver

# input source file (class name == source file basename)
mtx_convolver~.class.sources = mtx_convolver~.c ../src/array.c ../src/convolver.c
cflags=-I../include
ldlibs = -L. -L../bin/ -lfftw3f-3

# all extra files to be included in binary distribution of the library
datafiles =

CPPFLAGS=-I.

# include Makefile.pdlibbuilder from submodule directory 'pd-lib-builder'
PDLIBBUILDER_DIR=../pd-lib-builder/
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
