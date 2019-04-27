#!/bin/csh -f
#
# svn $Id: build.sh 835 2017-01-25 18:59:53Z arango $
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::: John Wilkin :::
# Copyright (c) 2002-2017 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS/TOMS Compiling Script                                            :::
#                                                                       :::
# Script to compile an user application where the application-specific  :::
# files are kept separate from the ROMS source code.                    :::
#                                                                       :::
# Q: How/why does this script work?                                     :::
#                                                                       :::
# A: The ROMS makefile configures user-defined options with a set of    :::
#    flags such as ROMS_APPLICATION. Browse the makefile to see these.  :::
#    If an option in the makefile uses the syntax ?= in setting the     :::
#    default, this means that make will check whether an environment    :::
#    variable by that name is set in the shell that calls make. If so   :::
#    the environment variable value overrides the default (and the      :::
#    user need not maintain separate makefiles, or frequently edit      :::
#    the makefile, to run separate applications).                       :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build.sh [options]                                               :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]      Compile in parallel using N CPUs                       :::
#                  omit argument for all available CPUs                 :::
#    -noclean    Do not clean already compiled objects                  :::
#                                                                       :::
# Notice that sometimes the parallel compilation fail to find MPI       :::
# include file "mpif.h".                                                :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

set parallel = 0
set clean = 1

while ( ($#argv) > 0 )
  switch ($1)
    case "-noclean"
      shift
      set clean = 0
    breaksw

    case "-j"
      shift
      set parallel = 1
      if (`echo $1 | grep '^[0-9]\+$'` != "" ) then
        set NCPUS = "-j $1"
        shift
      else
        set NCPUS = "-j"
      endif
    breaksw

    case "-*":
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]      Compile in parallel using N CPUs"
      echo "              omit argument for all avaliable CPUs"
      echo "-noclean    Do not clean already compiled objects"
      echo ""
      exit 1
    breaksw

  endsw
end

# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application
# CPP definitions.

setenv ROMS_APPLICATION      MCH_BIO_TBNT

# Set a local environmental variable to define the path to the directories
# where all this project's files are kept.

setenv MY_ROOT_DIR           ${HOME}/ROMS/roms854_TBNT
setenv MY_PROJECT_DIR        ${PWD}

# The path to the user's local current ROMS source code.
#
# If using svn locally, this would be the user's Working Copy Path (WCPATH).
# Note that one advantage of maintaining your source code locally with svn
# is that when working simultaneously on multiple machines (e.g. a local
# workstation, a local cluster and a remote supercomputer) you can checkout
# the latest release and always get an up-to-date customized source on each
# machine. This script is designed to more easily allow for differing paths
# to the code and inputs on differing machines.

#setenv MY_ROMS_SRC          ${MY_ROOT_DIR}/branches/arango
 setenv MY_ROMS_SRC          ${MY_ROOT_DIR}

# Set path of the directory containing makefile configuration (*.mk) files.
# The user has the option to specify a customized version of these files
# in a different directory than the one distributed with the source code,
# ${MY_ROMS_SRC}/Compilers. If this is the case, the you need to keep
# these configurations files up-to-date.

 setenv COMPILERS            ${MY_ROMS_SRC}/Compilers
#setenv COMPILERS            ${HOME}/Compilers

# Set tunable CPP options.
#
# Sometimes it is desirable to activate one or more CPP options to run
# different variants of the same application without modifying its header
# file. If this is the case, specify each options here using the -D syntax.
# Notice also that you need to use shell's quoting syntax to enclose the
# definition.  Both single or double quotes work. For example,
#
#    setenv MY_CPP_FLAGS "-DAVERAGES"
#    setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DDEBUGGING"
#
# can be use to write time-averaged fields. Notice that you can have as
# many definitions as you want by appending values.

setenv MY_CPP_FLAGS ""
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DDIAGNOSTICS_TS"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DDIAGNOSTICS_UV"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DSTATIONS"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DFLOATS"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DFLOAT_OYSTER"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DFLOAT_VWALK"

#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DDEBUGGING"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DOUT_DOUBLE"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DPOSITIVE_ZERO"

#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DUV_VIS2"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DUV_VIS4"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DMIX_S_UV"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DMIX_GEO_UV"

#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DTS_DIF2"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DMIX_S_TS"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DMIX_GEO_TS"

#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DUV_LOGDRAG"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DUV_LDRAG"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DUV_QDRAG"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DUV_DRAG_GRID"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DANA_DRAG"

#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DTS_U3HADVECTION"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DTS_C4VADVECTION"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DTS_MPDATA"

#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DBIO_FENNEL"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DECOSIM"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNEMURO"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNPZD_FRANKS"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNPZD_IRON"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNPZD_POWELL"

#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DGLS_MIXING"
#setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DMY25_MIXING"

# Set deprecated lateral boundary conditions CPP flags for backward
# compatibility with older versions of the code.

# setenv BACK_COMPATIBILITY  on          # needed for ROMS 3.4 or older

if ($?BACK_COMPATIBILITY) then
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DNORTHERN_WALL"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DSOUTHERN_WALL"
 setenv MY_CPP_FLAGS "${MY_CPP_FLAGS} -DEW_PERIODIC"
endif

# Other user defined environmental variables. See the ROMS makefile for
# details on other options the user might want to set here. Be sure to
# leave the switches meant to be off set to an empty string or commented
# out. Any string value (including off) will evaluate to TRUE in
# conditional if-statements.

 setenv USE_MPI             on          # distributed-memory parallelism
 setenv USE_MPIF90          on          # compile with mpif90 script
#setenv which_MPI           mpich       # compile with MPICH library
#setenv which_MPI           mpich2      # compile with MPICH2 library
 setenv which_MPI           openmpi     # compile with OpenMPI library

#setenv USE_OpenMP          on          # shared-memory parallelism

 setenv FORT                ifort
#setenv FORT                gfortran
#setenv FORT                pgi

#setenv USE_DEBUG           on          # use Fortran debugging flags
 setenv USE_LARGE           on          # activate 64-bit compilation
 setenv USE_NETCDF4         on          # compile with NetCDF-4 library
#setenv USE_PARALLEL_IO     on          # Parallel I/O with NetCDF-4/HDF5

#setenv USE_MY_LIBS         on          # use my library paths below

# There are several MPI libraries available. Here, we set the desired
# "mpif90" script to use during compilation. This only works if the make
# configuration file (say, Linux-pgi.mk) in the "Compilers" directory
# has the following definition for FC (Fortran Compiler) in the USE_MPI
# section:
#
#              FC := mpif90
#
# that is, "mpif90" defined without any path. Notice that the path
# where the MPI library is installed is computer dependent. Recall
# that you still need to use the appropriate "mpirun" to execute.

if ($?USE_MPIF90) then
  switch ($FORT)

    case "ifort"
      if ($which_MPI == "mpich" ) then
        setenv PATH /opt/intelsoft/mpich/bin:$PATH
      else if ($which_MPI == "mpich2" ) then
        setenv PATH /opt/intelsoft/mpich2/bin:$PATH
      else if ($which_MPI == "openmpi" ) then
#        setenv PATH /opt/intelsoft/openmpi/bin:$PATH
      endif
    breaksw

    case "pgi"
      if ($which_MPI == "mpich" ) then
        setenv PATH /opt/pgisoft/mpich/bin:$PATH
      else if ($which_MPI == "mpich2" ) then
        setenv PATH /opt/pgisoft/mpich2/bin:$PATH
      else if ($which_MPI == "openmpi" ) then
        setenv PATH /opt/pgisoft/openmpi/bin:$PATH
      endif
    breaksw

    case "gfortran"
      if ($which_MPI == "mpich2" ) then
        setenv PATH /opt/gfortransoft/mpich2/bin:$PATH
      else if ($which_MPI == "openmpi" ) then
        setenv PATH /opt/gfortransoft/openmpi/bin:$PATH
      endif
    breaksw

  endsw
endif

# If the USE_MY_LIBS is activated above, the path of the libraries
# required by ROMS can be set here using environmental variables
# which take precedence to the values specified in the make macro
# definitions file (Compilers/*.mk). For most applications, only
# the location of the NetCDF library is needed during compilation.
#
# Notice that when the USE_NETCDF4 macro is activated, we need the
# serial or parallel version of the NetCDF-4/HDF5 library. The
# configuration script NF_CONFIG (available since NetCDF 4.0.1)
# is used to set up all the required libraries according to the
# installed options (openDAP, netCDF4/HDF5 file format). The
# parallel library uses the MPI-I/O layer (usually available
# in MPICH2 and OpenMPI) requiring compiling with the selected
# MPI library.
#
# In ROMS distributed-memory applications, you may use either the
# serial or parallel version of the NetCDF-4/HDF5 library. The
# parallel version is required when parallel I/O is activated
# (ROMS cpp option PARALLEL_IO and HDF5).
#
# However, in serial or shared-memory ROMS applications, we need
# to use the serial version of the NetCDF-4/HDF5 to avoid conflicts
# with the compiler. We cannot activate MPI constructs in serial
# or shared-memory ROMS code. Hybrid parallelism is not possible.
#
# Recall also that the MPI library comes in several flavors:
# MPICH, MPICH2, OpenMPI, etc.

if ($?USE_MY_LIBS) then
  switch ($FORT)

    case "ifort"
      setenv ESMF_OS              Linux
      setenv ESMF_COMPILER        ifort
      setenv ESMF_BOPT            O
      setenv ESMF_ABI             64
      setenv ESMF_COMM            mpich
      setenv ESMF_SITE            default

      setenv ARPACK_LIBDIR        /opt/intelsoft/serial/ARPACK
      if ($?USE_MPI) then
        if ($which_MPI == "mpich" ) then
          setenv ESMF_DIR         /opt/intelsoft/mpich/esmf
          setenv MCT_INCDIR       /opt/intelsoft/mpich/mct/include
          setenv MCT_LIBDIR       /opt/intelsoft/mpich/mct/lib
          setenv PARPACK_LIBDIR   /opt/intelsoft/mpich/PARPACK
        else if ($which_MPI == "mpich2" ) then
          setenv ESMF_DIR         /opt/intelsoft/mpich2/esmf
          setenv MCT_INCDIR       /opt/intelsoft/mpich2/mct/include
          setenv MCT_LIBDIR       /opt/intelsoft/mpich2/mct/lib
          setenv PARPACK_LIBDIR   /opt/intelsoft/mpich2/PARPACK
        else if ($which_MPI == "openmpi" ) then
          setenv ESMF_DIR         /opt/intelsoft/openmpi/esmf
          setenv MCT_INCDIR       /opt/intelsoft/openmpi/mct/include
          setenv MCT_LIBDIR       /opt/intelsoft/openmpi/mct/lib
          setenv PARPACK_LIBDIR   /opt/intelsoft/openmpi/PARPACK
        endif
      endif

      if ($?USE_NETCDF4) then
        if ($?USE_PARALLEL_IO && $?USE_MPI) then
          if ($which_MPI == "mpich" ) then
            setenv NF_CONFIG      /opt/intelsoft/mpich/netcdf4/bin/nf-config
            setenv NETCDF_INCDIR  /opt/intelsoft/mpich/netcdf4/include
          else if ($which_MPI == "mpich2" ) then
            setenv NF_CONFIG      /opt/intelsoft/mpich2/netcdf4/bin/nf-config
            setenv NETCDF_INCDIR  /opt/intelsoft/mpich2/netcdf4/include
          else if ($which_MPI == "openmpi" ) then
            setenv NF_CONFIG      /opt/intelsoft/openmpi/netcdf4/bin/nf-config
            setenv NETCDF_INCDIR  /opt/intelsoft/openmpi/netcdf4/include
          endif
        else
          setenv NF_CONFIG        /opt/intelsoft/serial/netcdf4/bin/nf-config
          setenv NETCDF_INCDIR    /opt/intelsoft/serial/netcdf4/include
        endif
      else
        setenv NETCDF_INCDIR      /opt/intelsoft/serial/netcdf3/include
        setenv NETCDF_LIBDIR      /opt/intelsoft/serial/netcdf3/lib
      endif
    breaksw

    case "pgi"
      setenv ESMF_OS              Linux
      setenv ESMF_COMPILER        pgi
      setenv ESMF_BOPT            O
      setenv ESMF_ABI             64
      setenv ESMF_COMM            mpich
      setenv ESMF_SITE            default

      setenv ARPACK_LIBDIR        /opt/pgisoft/serial/ARPACK
      if ($?USE_MPI) then
        if ($which_MPI == "mpich" ) then
          setenv ESMF_DIR         /opt/pgisoft/mpich/esmf
          setenv MCT_INCDIR       /opt/pgisoft/mpich/mct/include
          setenv MCT_LIBDIR       /opt/pgisoft/mpich/mct/lib
          setenv PARPACK_LIBDIR   /opt/pgisoft/mpich/PARPACK
        else if ($which_MPI == "mpich2" ) then
          setenv ESMF_DIR         /opt/pgisoft/mpich2/esmf
          setenv MCT_INCDIR       /opt/pgisoft/mpich2/mct/include
          setenv MCT_LIBDIR       /opt/pgisoft/mpich2/mct/lib
          setenv PARPACK_LIBDIR   /opt/pgisoft/mpich2/PARPACK
        else if ($which_MPI == "openmpi" ) then
          setenv ESMF_DIR         /opt/pgisoft/openmpi/esmf
          setenv MCT_INCDIR       /opt/pgisoft/openmpi/mct/include
          setenv MCT_LIBDIR       /opt/pgisoft/openmpi/mct/lib
          setenv PARPACK_LIBDIR   /opt/pgisoft/openmpi/PARPACK
        endif
      endif

      if ($?USE_NETCDF4) then
        if ($?USE_PARALLEL_IO && $?USE_MPI) then
          if ($which_MPI == "mpich" ) then
            setenv NF_CONFIG      /opt/pgisoft/mpich/netcdf4/bin/nf-config
            setenv NETCDF_INCDIR  /opt/pgisoft/mpich/netcdf4/include
          else if ($which_MPI == "mpich2" ) then
            setenv NF_CONFIG      /opt/pgisoft/mpich2/netcdf4/bin/nf-config
            setenv NETCDF_INCDIR  /opt/pgisoft/mpich2/netcdf4/include
          else if ($which_MPI == "openmpi" ) then
            setenv NF_CONFIG      /opt/pgisoft/openmpi/netcdf4/bin/nf-config
            setenv NETCDF_INCDIR  /opt/pgisoft/openmpi/netcdf4/include
          endif
        else
          setenv NF_CONFIG        /opt/pgisoft/serial/netcdf4/bin/nf-config
          setenv NETCDF_INCDIR    /opt/pgisoft/serial/netcdf4/include
        endif
      else
        setenv NETCDF_INCDIR      /opt/pgisoft/serial/netcdf3/include
        setenv NETCDF_LIBDIR      /opt/pgisoft/serial/netcdf3/lib
      endif
    breaksw

    case "gfortran"
      setenv ESMF_OS              Linux
      setenv ESMF_COMPILER        gfortran
      setenv ESMF_BOPT            O
      setenv ESMF_ABI             64
      setenv ESMF_COMM            mpich
      setenv ESMF_SITE            default

      setenv ARPACK_LIBDIR        /opt/gfortransoft/serial/ARPACK
      if ($?USE_MPI) then
        if ($which_MPI == "mpich2" ) then
          setenv ESMF_DIR         /opt/gfortransoft/mpich2/esmf
          setenv MCT_INCDIR       /opt/gfortransoft/mpich2/mct/include
          setenv MCT_LIBDIR       /opt/gfortransoft/mpich2/mct/lib
          setenv PARPACK_LIBDIR   /opt/gfortransoft/mpich2/PARPACK
        else if ($which_MPI == "openmpi" ) then
          setenv ESMF_DIR         /opt/gfortransoft/openmpi/esmf
          setenv MCT_INCDIR       /opt/gfortransoft/openmpi/mct/include
          setenv MCT_LIBDIR       /opt/gfortransoft/openmpi/mct/lib
          setenv PARPACK_LIBDIR   /opt/gfortransoft/openmpi/PARPACK
        endif
      endif

      if ($?USE_NETCDF4) then
        if ($?USE_PARALLEL_IO && $?USE_MPI) then
          if ($which_MPI == "mpich2" ) then
            setenv NF_CONFIG      /opt/gfortransoft/mpich2/netcdf4/bin/nf-config
            setenv NETCDF_INCDIR  /opt/gfortransoft/mpich2/netcdf4/include
          else if ($which_MPI == "openmpi" ) then
            setenv NF_CONFIG      /opt/gfortransoft/openmpi/netcdf4/bin/nf-config
            setenv NETCDF_INCDIR  /opt/gfortransoft/openmpi/netcdf4/include
          endif
        else
          setenv NF_CONFIG        /opt/gfortransoft/serial/netcdf4/bin/nf-config
          setenv NETCDF_INCDIR    /opt/gfortransoft/serial/netcdf4/include
        endif
      else
        setenv NETCDF_INCDIR      /opt/gfortransoft/serial/netcdf3/include
        setenv NETCDF_LIBDIR      /opt/gfortransoft/serial/netcdf3/lib
      endif
    breaksw

  endsw
endif

# The rest of this script sets the path to the users header file and
# analytical source files, if any. See the templates in User/Functionals.
#
# If applicable, use the MY_ANALYTICAL_DIR directory to place your
# customized biology model header file (like fennel.h, nemuro.h, ecosim.h,
# etc).

 setenv MY_HEADER_DIR       ${MY_PROJECT_DIR}

 setenv MY_ANALYTICAL_DIR   ${MY_ROOT_DIR}/Functionals

# Put the binary to execute in the following directory.

 setenv BINDIR              ${MY_PROJECT_DIR}

# Put the f90 files in a project specific Build directory to avoid conflict
# with other projects.

 setenv SCRATCH_DIR         ${MY_PROJECT_DIR}/Build

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

 cd ${MY_ROMS_SRC}

# Stop if activating both MPI and OpenMP at the same time.

if ( ${?USE_MPI} & ${?USE_OpenMP} ) then
  echo "You cannot activate USE_MPI and USE_OpenMP at the same time!"
  exit 1
endif

# Remove build directory.

if ( $clean == 1 ) then
  make clean
endif

# Compile (the binary will go to BINDIR set above).

if ( $parallel == 1 ) then
  make $NCPUS
else
  make
endif
