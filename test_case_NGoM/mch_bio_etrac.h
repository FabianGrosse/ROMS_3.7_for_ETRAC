/*
**-----------------------------------------------------------------------------
**  MAIN CONFIGURATION OPTIONS
**-----------------------------------------------------------------------------
*/
#define BIO_FENNEL
#undef BIO_OXYSINK
#undef TXLA_MODEL
#define MCH_MODEL

#define MCH_BIO

#define PERFECT_RESTART

/*
**-----------------------------------------------------------------------------
**  OPTIONS FOR TBNT OUTPUT
**-----------------------------------------------------------------------------
*/
#define TBNT_OUT
#define TBNT_NFDOUBLE
#define PSRC_UV_BUGFIX

/*
**-----------------------------------------------------------------------------
**  OPTIONS FOR netCDF HANDLING
**-----------------------------------------------------------------------------
*/
#define DEFLATE
#define HDF5

/*
**-----------------------------------------------------------------------------
**  OPTIONS FOR BIOLOGICAL MODULE (FENNEL ET AL, 2006)
**-----------------------------------------------------------------------------
*/
#ifdef BIO_FENNEL
# define OXYGEN
# define RW14_OXYGEN_SC
# undef CARBON
# ifdef CARBON
#  define RW14_CO2_SC
#  define TALK_NONCONSERV
#  define PCO2AIR_DATA
#  undef PCO2AIR_SECULAR
# endif
# define PO4
# define BIO_SEDIMENT
# undef BIO_SOD
# undef SED_TEMP
# ifdef SED_TEMP
#  undef SED_TEMP_ADJUSTED
#  define SED_TEMP_KDT
#  undef SED_TEMP_KDT_FUTURE
# endif
# define DENITRIFICATION
# define RIVER_PON
# define DIAGNOSTICS_BIO
#endif

/*
**-----------------------------------------------------------------------------
**  OPTIONS FOR SIMPLE OXYGEN MODULE (YU ET AL, 2015)
**-----------------------------------------------------------------------------
*/
#ifdef BIO_OXYSINK
# undef BIO_WR
# define BIO_WR_MM
# undef BIO_SOD_TFIX
# undef BIO_SOD
# define SED_TEMP
# define SED_TEMP_KDT
# define DIAGNOSTICS_BIO
#endif

/*
**-----------------------------------------------------------------------------
**  OPTIONS FOR CIRCULATION MODEL
**-----------------------------------------------------------------------------
*/
#define DIAGNOSTICS_TS
#undef DIAGNOSTICS_UV
#undef AVERAGES
#define SOLAR_SOURCE
#define SWFRAC_KD
#define DIURNAL_SRFLUX
#define ANA_CLOUD

#define UV_ADV
#define UV_COR
#define UV_LOGDRAG
#define UV_VIS2
#define MIX_S_UV
#define VISC_GRID

#define TS_MPDATA
#undef TS_U3HADVECTION
#ifdef TS_U3HADVECTION
# define UPWIND_LIMITER
# define UPL_MASK
#endif
#define NONLIN_EOS
#define SALINITY
#define TS_DIF2
#define MIX_GEO_TS
#define DIFF_GRID

#define MASKING
#define CURVGRID
#define SOLVE3D
#define DJ_GRADPS

#define QCORRECTION
#define MY25_MIXING
#define N2S2_HORAVG

#define SPLINES_VVISC
#define SPLINES_VDIFF
#define RI_SPLINES

#define BULK_FLUXES
#define LONGWAVE

#define ANA_BSFLUX
#define ANA_BTFLUX

#define ANA_BPFLUX
#define ANA_SPFLUX

/*
**-----------------------------------------------------------------------------
**  OPTIONS FOR MCH MODEL
**-----------------------------------------------------------------------------
*/
#ifdef MCH_MODEL

# define  ANA_M2OBC
# define  ANA_FSOBC
# undef  GLS_MIXING
/*
**-----------------------------------------------------------------------------
**  OPTIONS FOR TXLA MODEL
**-----------------------------------------------------------------------------
*/
#elif defined TXLA_MODEL

# undef  ANA_M2OBC
# undef  ANA_FSOBC

# define EMINUSP

# undef ANA_SSFLUX

# undef TIDES
# ifdef TIDES
#  define SSH_TIDES
#  define UV_TIDES
#  define RAMP_TIDES
#  define ADD_FSOBC
#  define ADD_M2OBC
# endif
#endif
