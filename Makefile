# makefile for mcmax (with comments!)
# Tested on MacOSX 10.6 with ifort 11.1.080
# Tested on Fedora Core 8 with ifort 10.1.015

# compiler= FC, flags = FFlags
# linker= LINKER, flags= LDFLAGS, libraries=LIBS
FC	      = ifort
LINKER	      = ifort

# Platform specific compilation options
FLAG_ALL      = -O3 -extend-source -traceback -zero -prec-div
FLAG_LINUX    = -msse3 -prefetch
FLAG_MAC      = -mssse3 -opt-prefetch

ifeq ($(shell uname),Linux)
  FFLAGS   = $(FLAG_ALL) $(FLAG_LINUX) -diag-disable vec
  LDFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) -fpp Version.f
  LIBS     = -lm -lfftw3 -lcfitsio -I/sw/include -L/home/sw-astro/cfitsio64/lib
else
  FFLAGS  = $(FLAG_ALL) $(FLAG_MAC)
  LDFLAGS = $(FLAG_ALL) $(FLAG_MAC) -fpp Version.f
  LIBS	  = -lm -lfftw3 -lcfitsio -I/sw/include
endif

# use a suffix in file name (i.e. static, test etc.)
# cl> make suffix=32 ...
ifneq ($(name),)
  SUFFIX = -$(name)
endif

# files to make
OBJS	      = Modules.o \
		Main.o \
		DiskStructure.o \
		RegridR.o \
		Init.o \
		MCRadiation.o \
		Subroutines.o \
		Radiation.o \
		TraceMono.o \
		Telescope.o \
		Visibility.o \
		Diffusion.o \
		RandomWalk.o \
		Dtable.o \
		Lapack.o \
		Stochastic_MCMax.o \
		Settling.o \
		MPSettling.o \
		fit_module.o \
		GrainsizeDistribution.o \
		Topac.o \
		MRN.o \
		GasOpacity.o

# program name and install location
PROGRAM       = MCMax$(SUFFIX)-$(shell date +%d-%m-%Y)
DEST	      = ${HOME}/bin

# make actions 
all:		$(PROGRAM)
clean:;		rm -f $(OBJS) $(PROGRAM)
install:	$(PROGRAM)
		mv $(PROGRAM) $(DEST)
echo:;		@echo $(SUFFIX)

# special rule to make fortran 90 files
fit_module.o:	fit_module.f90
		${FC} $(FFLAGS) -c fit_module.f90 -o fit_module.o

# how to compile program 
$(PROGRAM):     $(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

# recompile everything if Modules.f has changed 
$(OBJS):	Modules.f
