/*
** svn $Id: ias.h 8 2007-02-06 19:00:29Z arango $
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for South Atlantic Bight and Gulf of Mexico Application.
**
*/

# undef SABGOM
# define DEFLATE
# define HDF5

# define UV_ADV
# define DJ_GRADPS
# define UV_COR
# define UV_QDRAG
# define UV_VIS2
# define MIX_S_UV
# define TS_MPDATA
# undef TS_U3HADVECTION
# define TS_DIF2
# define MIX_GEO_TS

# define SOLVE3D
# define SALINITY
# define NONLIN_EOS
# define CURVGRID
# define MASKING
# define SPLINES_VDIFF
# define SPLINES_VVISC
# define RI_SPLINES

# define AVERAGES
# undef AVERAGES_DETIDE
# define PERFECT_RESTART

/* surface forcing */
# define BULK_FLUXES
# ifdef BULK_FLUXES
# define DIURNAL_SRF
# define SOLAR_SOURCE
# define EMINUSP
# else
#  define ANA_SSFLUX
#  define ANA_SMFLUX
#  define ANA_STFLUX
# endif

/* Select turbulence mixing scheme */

# define MY25_MIXING
# ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
# endif


/* Select  Biological model option */

# undef RIVER_BIOLOGY

# undef  BIO_FENNEL
#ifdef BIO_FENNEL
# define CARBON
!# define PCO2_RZ
# define TALK_NONCONSERV
# define DENITRIFICATION
# define BIO_SEDIMENT
# define DIAGNOSTICS_BIO
# define OXYGEN
#endif

# define  BIO_HYDROCARBON
#ifdef BIO_HYDROCARBON
# define OXYGEN
# define DIAGNOSTICS_BIO
#endif
# define OIL_SPILL

/* Tide */
# define SSH_TIDES
# define ADD_FSOBC
# define UV_TIDES
# define ADD_M2OBC
# define RAMP_TIDES

# define RADIATION_2D

# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX

# undef ANA_NUDGCOEF

