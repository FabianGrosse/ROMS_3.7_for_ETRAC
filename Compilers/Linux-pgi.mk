# svn $Id: Linux-pgi.mk 834 2017-01-25 18:49:17Z arango $
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2017 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for PGI Fortran compiler on Linux
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF libary directory
# LD             Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := mpif90
           FFLAGS :=
              CPP := /usr/bin/cpp
         CPPFLAGS := -P -traditional
               LD := $(FC)
          LDFLAGS :=
               AR := ar
          ARFLAGS := r
            MKDIR := mkdir -p
               RM := rm -f
           RANLIB := ranlib
             PERL := perl
             TEST := test

        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)

#
# Perform floating-point operations in strict conformance with the
# IEEE standard. This may slow down computations because some
# optimizations are disabled.  However, we noticed a speed-up.
# The user may want to uncomment this option to allow similar,
# if not identical solutions between different of the PGI compiler.

#           FFLAGS += -Kieee

#
# Library locations, can be overridden by environment variables.
#

ifdef USE_NETCDF4
    NETCDF_INCDIR ?= /opt/include
    NETCDF_LIBDIR ?= /opt/lib
      HDF5_LIBDIR := /opt/lib
else
    NETCDF_INCDIR ?= /opt/include
    NETCDF_LIBDIR ?= /opt/lib
endif
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf -lmpi -lmpi_f77 -lmpi_f90
ifdef USE_NETCDF4
             LIBS += -L$(HDF5_LIBDIR) -lhdf5_hl -lhdf5 -lz
endif

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /opt/pgisoft/PARPACK
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /opt/pgisoft/PARPACK
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
 ifdef USE_MPIF90
               FC := mpif90
 else
             LIBS += -L$(NETCDF_LIBDIR) -lnetcdf  -lmpi -lmpi_f77 -lmpi_f90
 endif
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
endif

# According to the PGI manual, the -u -Bstatic flags initializes
# the symbol table with -Bstatic, which is undefined for the linker.
# An undefined symbol triggers loading of the first member of an
# archive library. The -u flag fails with version 7.x of the compiler
# because it expects an argument.

ifdef USE_DEBUG
#          FFLAGS += -g -C -Mchkstk -Mchkfpstk
           FFLAGS += -g -C
#          FFLAGS += -gopt -C
#          FFLAGS += -g
else
#          FFLAGS += -Bstatic -fastsse -Mipa=fast
           FFLAGS += -O3 -fastsse -tp k8-64e -Mipa=fast
endif

# Save compiler flags without the MCT or ESMF libraries additions
# to keep the string (MY_FFLAGS) in "mod_strings.o" short. Otherwise,
# it will exceed the maximum number of characters allowed for
# free-format compilation.

ifdef USE_MCT
       MCT_INCDIR ?= /opt/pgisoft/mct/include
       MCT_LIBDIR ?= /opt/pgisoft/mct/lib
           FFLAGS += -I$(MCT_INCDIR)
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
endif

ifdef USE_ESMF
      ESMF_MK_DIR ?= /opt/pgisoft/esmf/lib/libO/Linux.pgi.64.mpich.default 
                     include $(ESMF_MK_DIR)/esmf.mk
           FFLAGS += $(ESMF_F90COMPILEPATHS)
             LIBS += $(ESMF_F90LINKPATHS) -lesmf -lC
endif

       clean_list += ifc* work.pc*

#
# Use full path of compiler.
#
               FC := $(shell which ${FC})
               LD := $(FC)

#
# Set free form format in source files to allow long string for
# local directory and compilation flags inside the code.
#

$(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/mod_strings.o: FFLAGS := $(MY_FFLAGS) -Mfree
$(SCRATCH_DIR)/analytical.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/biology.o: FFLAGS += -Mfree
ifdef USE_ADJOINT
$(SCRATCH_DIR)/ad_biology.o: FFLAGS += -Mfree
endif
ifdef USE_REPRESENTER
$(SCRATCH_DIR)/rp_biology.o: FFLAGS += -Mfree
endif
ifdef USE_TANGENT
$(SCRATCH_DIR)/tl_biology.o: FFLAGS += -Mfree
endif

#
# Supress free format in SWAN source files since there are comments
# beyond column 72.
#

ifdef USE_SWAN

$(SCRATCH_DIR)/ocpcre.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/ocpids.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/ocpmix.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom1.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom2.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom3.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom4.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swancom5.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanmain.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanout1.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanout2.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanparll.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swanser.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swmod1.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/swmod2.o: FFLAGS += -Mnofree
$(SCRATCH_DIR)/m_constants.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/m_fileio.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/mod_xnl4v5.o: FFLAGS += -Mfree
$(SCRATCH_DIR)/serv_xnl4v5.o: FFLAGS += -Mfree

endif
