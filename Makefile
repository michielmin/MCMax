# makefile for mcmax (with comments!)
# Tested on Fedora Core 8 with ifort 10.1.015 (20/12/2012)
# Tested on MacOSX 10.6 with ifort 11.1.080 (20/12/2012)
# Tested on MacOSX 10.8 with ifort 14.0.1 (30/12/2013)
# Tested on MacOSX 10.9 with ifort 12.0.0 (20/12/2013)
# Tested on Scientific Linux 7.7 with with GNU Fortran (GCC) 6.3.1 (14/02/2020)
	
#GITVERSION = $(echo "#define gitversion = \"$(shell git rev-parse HEAD)\"" > gitversion.h)
GITVERSION = $(ls -l)


# compiler= FC, flags = FFlags
# linker= LINKER, flags= LDFLAGS, libraries=LIBS
#################################################################
# Default compiler path (ifort)
ifneq ($(compiler),gfortran)
	# compiler= FC, flags = FFlags
	# linker= LINKER, flags= LDFLAGS, libraries=LIBS
	FC	      = ifort
	LINKER	      = ifort

	# enforce single core compilation with:
	# cl> make multi=false
	ifneq ($(multi),false)
	  ifeq ($(shell uname),Linux)
	    MULTICORE = -qopenmp -fp-model strict
	  else
	    MULTICORE = -openmp -fp-model strict  
	  endif
	endif

	# array boundary check
	ifeq ($(debug),true)
	  DEBUGGING = -debug -check bounds -ftrapuv -fpe0 -O0 -check all -fp-stack-check -CB -C
	  #DEBUGGING = -debug -ftrapuv -g -check all -fp-stack-check
	  #DEBUGGING = -heap-arrays
	  #DEBUGGING = -gen-interfaces -warn interfaces
	endif

	# Platform specific compilation options
	FLAG_ALL      = -O3 -extend-source -zero -prec-div $(MULTICORE) $(DEBUGGING)
	FLAG_LINUX    = -msse3 #-prefetch
	FLAG_MAC      = -xHOST -static-intel #-opt-prefetch 


	ifeq ($(shell uname),Linux)
	  FFLAGS   = $(FLAG_ALL) $(FLAG_LINUX) -diag-disable vec
	  LDFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) -fpp Version.f
	  LIBS     = -lm -lfftw3 -lcfitsio -I/sw/include -L/home/sw-astro/cfitsio/lib -L$(HOME)/lib
	else
	  FFLAGS  = $(FLAG_ALL) $(FLAG_MAC) -diag-disable 8290,8291
	  LDFLAGS = $(FLAG_ALL) $(FLAG_MAC) -fpp Version.f
	  #LIBS	  =  -L/sw/lib -lm -lfftw3 -lcfitsio -I/sw/include
	  LIBS	  =  -L/sw/lib -lm -lfftw3 -L/usr/local/lib -lcfitsio -L/opt/local/lib -L$(HOME)/lib
	endif

#################################################################
# GFOTRAN path 
else
	FC	      = gfortran
	LINKER	  = gfortran
	
  	# cl> make multi=false
  # maybe remove this as it does  not compile without openmp anyway
  ifneq ($(multi),false)
  	MULTICORE = -fopenmp
  endif
  
  # array boundary check
  ifeq ($(debug),true)
    DEBUGGING =  -fcheck=all
  endif
  
  # -fno-whole-file deactivates some checks (e.g. explicit interface required for optional arguments)
  #FLAG_ALL = -O1 -g -fbacktrace -fcheck=all -finit-local-zero -ffixed-line-length-none $(MULTICORE) $(DEBUGGING)
  FLAG_ALL = -O3 -g -fbacktrace -finit-local-zero -ffixed-line-length-none $(MULTICORE) $(DEBUGGING)
    
  FFLAGS   = $(FLAG_ALL)
  LDFLAGS  = $(FLAG_ALL) -cpp Version.f

# add here library paths if required  
	ifeq ($(shell uname),Linux)
		LIBS   = -lm -lfftw3 -lcfitsio
	else
		LIBS   = -lm -lfftw3 -lcfitsio -L/opt/local/lib
	endif
endif

# use a suffix in file name (i.e. static, test etc.)
# cl> make name=test
ifneq ($(name),)
  SUFFIX = -$(name)
endif

# files to make
OBJS	      = Modules.o \
		InputOutput.o \
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
		ComputePart.o \
		ComputePAH.o \
		DeadZone.o \
		RadialStruct.o \
		KuruczData.o \
		RefIndData.o \
		TraceOpticalDepth.o \
		ISRFData.o

# program name and install location
PROGRAM       = MCMax$(SUFFIX)-$(shell date +%d-%m-%Y)
DEST	      = ${HOME}/bin

# make actions 
all:		version $(PROGRAM)
version:;	echo "#define gitversion \"$(shell git rev-parse HEAD)\"" > gitversion.h
clean:;		rm -f $(OBJS) $(PROGRAM) *.mod *.i
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

# recompile everything if InputOutput.f has changed 
$(OBJS):	InputOutput.f



