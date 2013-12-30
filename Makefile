# makefile for mcmax (with comments!)
# Tested on Fedora Core 8 with ifort 10.1.015 (20/12/2012)
# Tested on MacOSX 10.6 with ifort 11.1.080 (20/12/2012)
# Tested on MacOSX 10.8 with ifort 14.0.1 (30/12/2013)
# Tested on MacOSX 10.9 with ifort 12.0.0 (20/12/2013)
	
#GITVERSION = $(echo "#define gitversion = \"$(shell git rev-parse HEAD)\"" > gitversion.h)
GITVERSION = $(ls -l)

# compiler= FC, flags = FFlags
# linker= LINKER, flags= LDFLAGS, libraries=LIBS
FC	      = ifort
LINKER	      = ifort

# enforce single core compilation with:
# cl> make multi=false
ifneq ($(multi),false)
  MULTICORE = -openmp -fp-model strict
endif

# array boundary check
ifeq ($(debug),true)
  DEBUGGING = -debug -check bounds -ftrapuv #-fpe0 -O0
  #DEBUGGING = -debug -ftrapuv -g -check all -fp-stack-check
  #DEBUGGING = -heap-arrays
  #DEBUGGING = -gen-interfaces -warn interfaces
endif

# Platform specific compilation options
FLAG_ALL      = -O3 -extend-source -g -traceback -zero -prec-div $(MULTICORE) $(DEBUGGING)
FLAG_LINUX    = -msse3 -prefetch
FLAG_MAC      = -mssse3 -opt-prefetch -static-intel


ifeq ($(shell uname),Linux)
  FFLAGS   = $(FLAG_ALL) $(FLAG_LINUX) -diag-disable vec
  LDFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) -fpp Version.f
  LIBS     = -lm -lfftw3 -lcfitsio -I/sw/include -L/home/sw-astro/cfitsio/lib
else
  FFLAGS  = $(FLAG_ALL) $(FLAG_MAC) -diag-disable 8290,8291
  LDFLAGS = $(FLAG_ALL) $(FLAG_MAC) -fpp Version.f -Wl,-macosx_version_min,10.6
  #LIBS	  =  -L/sw/lib -lm -lfftw3 -lcfitsio -I/sw/include
  LIBS	  =  -lm -lfftw3 -L/usr/local/lib -lcfitsio -L/opt/local/lib
endif

# use a suffix in file name (i.e. static, test etc.)
# cl> make name=test
ifneq ($(name),)
  SUFFIX = -$(name)
endif

# files to make
OBJS	      = Modules.o \
		Main.o \
		DiskStructure.o \
		RegridR.o \
		Init.o \
		ImportProdimo.o \
		ExportProdimo.o \
		Keywordlist.o \
		CreateLRF.o \
		MCRadiation.o \
		Subroutines.o \
		Radiation.o \
		TraceMono.o \
		TraceLine.o \
		Telescope.o \
		Visibility.o \
		Diffusion.o \
		RandomWalk.o \
		Dtable.o \
		Lapack.o \
		Stochastic_MCMax.o \
		PAH_MCMax.o \
		Ionization.o \
		Settling.o \
		MPSettling.o \
		fit_module.o \
		GrainsizeDistribution.o \
		Topac.o \
		MRN.o \
		GasOpacity.o \
		ComputeLTE.o \
		Reddening.o \
		Smooth.o \
		ReadParticleFits.o \
		InputOutput.o \
		ComputePart.o \
		ComputePAH.o \
		DeadZone.o \
		RadialStruct.o

# program name and install location
PROGRAM       = MCMax$(SUFFIX)-$(shell date +%d-%m-%Y)
DEST	      = ${HOME}/bin

# make actions 
all:		version $(PROGRAM)
version:;	echo "#define gitversion \"$(shell git rev-parse HEAD)\"" > gitversion.h
clean:;		rm -f $(OBJS) $(PROGRAM)
install:	$(PROGRAM)
		mv $(PROGRAM) $(DEST)
echo:;		@echo $(SUFFIX)

# special rule to make fortran 90 files
fit_module.o:	fit_module.f90
		${FC} $(FFLAGS) -c fit_module.f90 -o fit_module.o

# special rule to make Regrid.f files for fast compilation
#RegridR.o:	RegridR.f
#		${FC} $(FFLAGS) -O0 -c RegridR.f -o RegridR.o
#Main.o:	Main.f
#		${FC} $(FFLAGS) -O0 -c Main.f -o Main.o
#Init.o:	Init.f
#		${FC} $(FFLAGS) -O0 -c Init.f -o Init.o

# how to compile program 
$(PROGRAM):     $(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

# recompile everything if Modules.f has changed 
$(OBJS):	Modules.f
