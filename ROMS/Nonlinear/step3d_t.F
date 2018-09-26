#include "cppdefs.h"
      MODULE step3d_t_mod
#if !defined TS_FIXED && (defined NONLINEAR && defined SOLVE3D)
!
!svn $Id: step3d_t.F 854 2017-07-18 23:28:45Z arango $
!=======================================================================
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This routine time-steps tracer equations.  Notice that advective    !
!  and diffusive terms are time-stepped differently. It applies the    !
!  corrector time-step for horizontal/vertical advection,  vertical    !
!  diffusion, nudging if necessary, and lateral boundary conditions.   !
!                                                                      !
!  Notice that at input the tracer arrays have:                        !
!                                                                      !
!    t(:,:,:,nnew,:)   m Tunits  n+1     horizontal/vertical diffusion !
!                                        terms plus source/sink terms  !
!                                        (biology, sediment), if any   !
!                                                                      !
!    t(:,:,:,3   ,:)   Tunits    n+1/2   advective terms and vertical  !
!                                        diffusion predictor step      !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: step3d_t
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE step3d_t (ng, tile)
!***********************************************************************
!
      USE mod_param
# ifdef DIAGNOSTICS_TS
      USE mod_diags
# endif
# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION
      USE mod_tbnt
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 35, __LINE__, __FILE__)
# endif
      CALL step3d_t_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nrhs(ng), nstp(ng), nnew(ng),                 &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
     &                    GRID(ng) % umask,                             &
     &                    GRID(ng) % vmask,                             &
# endif
# ifdef UPL_MASK
     &                    GRID(ng) % uplmask,                           &
# endif
# if (defined TS_MPDATA || defined TS_HSIMT) && defined WET_DRY
     &                    GRID(ng) % rmask_wet,                         &
     &                    GRID(ng) % umask_wet,                         &
     &                    GRID(ng) % vmask_wet,                         &
# endif
# if defined TS_MPDATA || \
     (defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION)
     &                    GRID(ng) % omn,                               &
# endif
# ifdef TS_MPDATA
     &                    GRID(ng) % om_u,                              &
     &                    GRID(ng) % om_v,                              &
     &                    GRID(ng) % on_u,                              &
     &                    GRID(ng) % on_v,                              &
# endif
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
     &                    GRID(ng) % Huon,                              &
     &                    GRID(ng) % Hvom,                              &
     &                    GRID(ng) % z_r,                               &
     &                    MIXING(ng) % Akt,                             &
     &                    OCEAN(ng) % W,                                &
# if defined FLOATS && defined FLOAT_VWALK
     &                    MIXING(ng) % dAktdz,                          &
# endif
# ifdef DIAGNOSTICS_TS
     &                    DIAGS(ng) % DiaTwrk,                          &
# endif
# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION
     &                    TBNT(ng) % TBNToce2tbnt,                      &
     &                    TBNT(ng) % TBNTphysWRK,                       &
     &                    TBNT(ng) % TBNTpsrcWRK,                       &
# endif
# if defined UPWIND_LIMITER_DIAG
     &                    MIXING(ng) % LimP,                            &
     &                    MIXING(ng) % LimN,                            &
     &                    MIXING(ng) % LmD1,                            &
     &                    MIXING(ng) % LmD2,                            &
     &                    MIXING(ng) % LmD3,                            &
     &                    MIXING(ng) % LmD4,                            &
# endif
     &                    OCEAN(ng) % t)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 35, __LINE__, __FILE__)
# endif
      RETURN
      END SUBROUTINE step3d_t
!
!***********************************************************************
      SUBROUTINE step3d_t_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nrhs, nstp, nnew,                       &
# ifdef MASKING
     &                          rmask, umask, vmask,                    &
# endif
# ifdef UPL_MASK
     &                          uplmask,                                & 
# endif
# if (defined TS_MPDATA || defined TS_HSIMT) && defined WET_DRY
     &                          rmask_wet, umask_wet, vmask_wet,        &
# endif
# if defined TS_MPDATA || \
     (defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION)
     &                          omn,                                    &
# endif
# ifdef TS_MPDATA
     &                          om_u, om_v, on_u, on_v,                 &
# endif
     &                          pm, pn,                                 &
     &                          Hz, Huon, Hvom,                         &
     &                          z_r,                                    &
     &                          Akt,                                    &
     &                          W,                                      &
# if defined FLOATS && defined FLOAT_VWALK
     &                          dAktdz,                                 &
# endif
# ifdef DIAGNOSTICS_TS
     &                          DiaTwrk,                                &
# endif
# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION
     &                          TBNToce2tbnt,                           &
     &                          TBNTphysWRK,                            &
     &                          TBNTpsrcWRK,                            &
# endif
# if defined UPWIND_LIMITER_DIAG
                                LimP,LimN,LmD1,LmD2,LmD3,LmD4,          &
# endif
     &                          t)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_ncparam
# if defined NESTING && !defined ONE_WAY
      USE mod_nesting
# endif
      USE mod_scalars
      USE mod_sources
!
      USE exchange_3d_mod, ONLY : exchange_r3d_tile
# ifdef DISTRIBUTE
#  if defined FLOATS && defined FLOAT_VWALK
      USE mp_exchange_mod, ONLY : mp_exchange3d
#  endif
      USE mp_exchange_mod, ONLY : mp_exchange4d
# endif
# ifdef TS_MPDATA
      USE mpdata_adiff_mod
# endif
# ifdef TS_HSIMT
      USE hsimt_tvd_mod
# endif
# ifdef NESTING
      USE nesting_mod, ONLY : bry_fluxes
# endif
      USE t3dbc_mod, ONLY : t3dbc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs, nstp, nnew
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef UPL_MASK
      real(r8), intent(in) :: uplmask(LBi:,LBj:)
#  endif
#  if (defined TS_MPDATA || defined TS_HSIMT) && defined WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:,LBj:)
      real(r8), intent(in) :: umask_wet(LBi:,LBj:)
      real(r8), intent(in) :: vmask_wet(LBi:,LBj:)
#  endif
#  if defined TS_MPDATA || \
      (defined TBNT_OUT && defined BIO_FENNEL && \
       defined BIO_SEDIMENT && defined DENITRIFICATION)
      real(r8), intent(in) :: omn(LBi:,LBj:)
#  endif
#  ifdef TS_MPDATA
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
#  ifdef SUN
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
#  else
      real(r8), intent(in) :: Akt(LBi:,LBj:,0:,:)
#  endif
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:,LBj:,:,:,:)
#  endif
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
      integer , intent(in)    :: TBNToce2tbnt(NT(ng))
      real(r8), intent(inout) :: TBNTphysWRK(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: TBNTpsrcWRK(LBi:,LBj:,:,:,:)
#  endif
#  if defined UPWIND_LIMITER_DIAG
      real(r8), intent(inout) :: LimP(LBi:,LBj:,:)
      real(r8), intent(inout) :: LimN(LBi:,LBj:,:)
      real(r8), intent(inout) :: LmD1(LBi:,LBj:,:)
      real(r8), intent(inout) :: LmD2(LBi:,LBj:,:)
      real(r8), intent(inout) :: LmD3(LBi:,LBj:,:)
      real(r8), intent(inout) :: LmD4(LBi:,LBj:,:)
#  endif
#  ifdef SUN
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#  else
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#  endif
#  if defined FLOATS && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:,LBj:,:)
#  endif

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef UPL_MASK
      real(r8), intent(in) :: uplmask(LBi:UBi,LBj:UBj)
#  endif
#  if (defined TS_MPDATA || defined TS_HSIMT) && defined WET_DRY
      real(r8), intent(in) :: rmask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask_wet(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask_wet(LBi:UBi,LBj:UBj)
#  endif
#  if defined TS_MPDATA || \
      (defined TBNT_OUT && defined BIO_FENNEL && \
       defined BIO_SEDIMENT && defined DENITRIFICATION)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
#  endif
#  ifdef TS_MPDATA
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Huon(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Hvom(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
      real(r8), intent(in) :: W(LBi:UBi,LBj:UBj,0:N(ng))

#  ifdef DIAGNOSTICS_TS
      real(r8), intent(inout) :: DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),  &
     &                                   NDT)
#  endif
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
      integer , intent(in)    :: TBNToce2tbnt(NT(ng))
      real(r8), intent(inout) :: TBNTphysWRK(LBi:UBi,LBj:UBj,N(ng),     &
     &                                       NTBNTTRC(ng),nTBNTphys)
      real(r8), intent(inout) :: TBNTpsrcWRK(LBi:UBi,LBj:UBj,N(ng),     &
     &                                       NTBNTTRC(ng),nTBNTpsrc)
#  endif
#  if defined UPWIND_LIMITER_DIAG
      real(r8), intent(inout) :: LimP(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: LimN(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: LmD1(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: LmD2(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: LmD3(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: LmD4(LBi:UBi,LBj:UBj,N(ng))
#  endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

#  if defined FLOATS && defined FLOAT_VWALK
      real(r8), intent(out) :: dAktdz(LBi:UBi,LBj:UBj,N(ng))
#  endif
# endif
!
!  Local variable declarations.
!
# ifdef NESTING
      integer :: ILB, IUB, JLB, JUB
      integer :: dg, cr, rg
# endif
      integer :: Isrc, Jsrc
      integer :: i, ic, is, itrc, j, k, ltrc
# if defined AGE_MEAN && defined T_PASSIVE
      integer :: iage
# endif
# ifdef DIAGNOSTICS_TS
      integer :: idiag
# endif
      real(r8), parameter :: eps = 1.0E-16_r8

# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION
      real(r8), parameter :: TBNTeps = 1.0E-10_r8
# endif

      real(r8) :: cff, cff1, cff2, cff3

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: oHz
# ifdef TS_MPDATA
#  ifdef DIAGNOSTICS_TS
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,3) :: Dhadv
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng),NT(ng)) :: Dvadv
#  endif
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng),NT(ng)) :: Ta
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: Ua
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: Va
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: Wa
# endif

# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION

      logical, dimension(IminS:ImaxS,JminS:JmaxS,N(ng),                 &
     &                               NTBNTTRC(ng)) :: psrcMaskX
      logical, dimension(IminS:ImaxS,JminS:JmaxS,N(ng),                 &
     &                               NTBNTTRC(ng)) :: psrcMaskE

      integer  :: itbnt
      real(r8) :: psrcIn
      real(r8), dimension(IminS:ImaxS) :: vdif_bot
      real(r8), dimension(IminS:ImaxS,0:N(ng),2) :: TBNT_FC
#  ifdef TS_MPDATA
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: TBNTadh
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng),                &
     &                                 NTBNTTRC(ng)) :: TBNTadv
#  endif
# endif
# ifdef TS_HSIMT
      real(r8):: rd, rkad, a1, b1, betad, sw, ru, rkau
      real(r8):: betau, epson, cc1, cc2, cc3
      real(r8), dimension(0:N(ng)) :: kaz, kaz_inverse, grad_k
# endif
# ifdef UPWIND_LIMITER
      real(r8), parameter :: alpha_limiter=5.0_r8
      real(r8), parameter :: uplim = 1

      real(r8) :: alpha
      real(r8) :: FX_tmp, FX_med_pos, FX_pos, FX_med_neg, FX_neg
      real(r8) :: FE_tmp, FE_med_pos, FE_pos, FE_med_neg, FE_neg

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX_EE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX_E
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX_W
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX_WW

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE_NN
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE_N
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE_S
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE_SS
# endif

# include "set_bounds.h"

# ifdef NESTING
!
!  Notice that the trace flux boundary arrays are dimensioned with the
!  global dimensions of grid to facilitate processing.
!
      ILB=BOUNDS(ng)%LBi(-1)
      IUB=BOUNDS(ng)%UBi(-1)
      JLB=BOUNDS(ng)%LBj(-1)
      JUB=BOUNDS(ng)%UBj(-1)
# endif
!
!-----------------------------------------------------------------------
!  Time-step horizontal advection term.
!-----------------------------------------------------------------------
!
!  Compute inverse thickness.
!
# if defined TS_MPDATA || defined TS_HSIMT
#  define I_RANGE Istrm2,Iendp2
#  define J_RANGE Jstrm2,Jendp2
# else
#  define I_RANGE Istr,Iend
#  define J_RANGE Jstr,Jend
# endif

      DO k=1,N(ng)
        DO j=J_RANGE
          DO i=I_RANGE
            oHz(i,j,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
      END DO
# undef I_RANGE
# undef J_RANGE
# if defined TS_MPDATA || defined TS_HSIMT
!
!  The MPDATA algorithm requires a three-point footprint, so exchange
!  boundary data on t(:,:,:,nnew,:) so other processes computed earlier
!  (horizontal diffusion, biology, or sediment) are accounted.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        DO itrc=1,NT(ng)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END DO
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,nnew,:))
#  endif
# endif
!
!  Compute horizontal tracer advection fluxes.
!
      T_LOOP : DO itrc=1,NT(ng)
# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION
        itbnt = TBNToce2tbnt(itrc)
# endif
        K_LOOP : DO k=1,N(ng)

# if defined TS_C2HADVECTION
!
!  Second-order, centered differences horizontal advective fluxes.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              FX(i,j)=Huon(i,j,k)*                                      &
     &                0.5_r8*(t(i-1,j,k,3,itrc)+                        &
     &                        t(i  ,j,k,3,itrc))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              FE(i,j)=Hvom(i,j,k)*                                      &
     &                0.5_r8*(t(i,j-1,k,3,itrc)+                        &
     &                        t(i,j  ,k,3,itrc))
            END DO
          END DO
# elif defined TS_HSIMT
!
!  Third-order HSIMT-TVD horizontal advective fluxes.
!
          CALL hsimt_tvd_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
#  ifdef MASKING
     &                         rmask, umask, vmask,                     &
#  endif
#  ifdef WET_DRY
     &                         rmask_wet, umask_wet, vmask_wet,         &
#  endif
     &                         pm, pn,                                  &
     &                         Huon(:,:,k), Hvom(:,:,k),                &
     &                         oHz(:,:,k), t(:,:,k,3,itrc),             &
     &                         FX, FE)

# elif defined TS_MPDATA
!
!  First-order, upstream differences horizontal advective fluxes.
!
          DO j=JstrVm2,Jendp2i
            DO i=IstrUm2,Iendp3
              cff1=MAX(Huon(i,j,k),0.0_r8)
              cff2=MIN(Huon(i,j,k),0.0_r8)
              FX(i,j)=cff1*t(i-1,j,k,3,itrc)+                           &
     &                cff2*t(i  ,j,k,3,itrc)
            END DO
          END DO
          DO j=JstrVm2,Jendp3
            DO i=IstrUm2,Iendp2i
              cff1=MAX(Hvom(i,j,k),0.0_r8)
              cff2=MIN(Hvom(i,j,k),0.0_r8)
              FE(i,j)=cff1*t(i,j-1,k,3,itrc)+                           &
     &                cff2*t(i,j  ,k,3,itrc)
            END DO
          END DO

# else
!
#  if defined TS_U3HADVECTION
!  Third-order, uptream-biased horizontal advective fluxes.
#   if defined UPWIND_LIMITER
!  Prepare the limiting values East/West
!  Note that valid values need to span i-2:i+1 cells for each cell face
!            write(*,*)'Jstr, Jend, Istr, Iend',Jstr,Jend,Istr,Iend
!            write(*,*)'IminS,ImaxS,JminS,JmaxS',IminS,ImaxS,JminS,JmaxS
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
               FX_E(i,j) = t(i  ,j,k,3,itrc)*umask(i,j)
               FX_EE(i,j)= t(i+1,j,k,3,itrc)*umask(i,j)
               FX_W(i,j) = t(i-1,j,k,3,itrc)*umask(i,j)
               FX_WW(i,j)= t(i-2,j,k,3,itrc)*umask(i,j)
            END DO
          END DO
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
               FE_N(i,j) = t(i,j  ,k,3,itrc)*vmask(i,j)
               FE_NN(i,j)= t(i,j+1,k,3,itrc)*vmask(i,j)
               FE_S(i,j) = t(i,j-1,k,3,itrc)*vmask(i,j)
               FE_SS(i,j)= t(i,j-2,k,3,itrc)*vmask(i,j)
            END DO
          END DO
#   endif
#  elif defined TS_A4HADVECTION
!  Fourth-order, Akima horizontal advective fluxes.
#  else
!  Fourth-order, centered differences horizontal advective fluxes.
#  endif
!
          DO j=Jstr,Jend
            DO i=Istrm1,Iendp2
              FX(i,j)=t(i  ,j,k,3,itrc)-                                &
     &                t(i-1,j,k,3,itrc)
#  ifdef MASKING
              FX(i,j)=FX(i,j)*umask(i,j)
#  endif
            END DO
          END DO
          IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Western_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Istr-1,j)=FX(Istr,j)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
            IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
              DO j=Jstr,Jend
                FX(Iend+2,j)=FX(Iend+1,j)
              END DO
            END IF
          END IF
!
          DO j=Jstr,Jend
            DO i=Istr-1,Iend+1
#  if defined TS_U3HADVECTION
              curv(i,j)=FX(i+1,j)-FX(i,j)
#  elif defined TS_A4HADVECTION
              cff=2.0_r8*FX(i+1,j)*FX(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#  else
              grad(i,j)=0.5_r8*(FX(i+1,j)+FX(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
#  ifdef TS_U3HADVECTION
#   ifdef UPWIND_LIMITER
#    ifdef UPL_MASK
              IF (uplmask(i,j).eq.1) THEN
#    endif
              FX(i,j)=0.5_r8*                                           &
     &                (t(i-1,j,k,3,itrc)+                               &
     &                 t(i  ,j,k,3,itrc))-                              &
     &         cff1*(curv(i-1,j)*MAX(SIGN(1.0_r8,Huon(i,j,k)),0.0_r8)-  &
     &                 curv(i,j)*MIN(SIGN(1.0_r8,Huon(i,j,k)),0.0_r8))
              FX_tmp=FX_WW(i,j)+alpha_limiter*(FX_W(i,j)-FX_WW(i,j))
              FX_med_pos=MAX(MIN(FX_W(i,j),FX_E(i,j)),                  &
     &                       MIN(FX_E(i,j),FX_tmp),                     &
     &                       MIN(FX_tmp,FX_W(i,j)))
              FX_pos=MAX(MIN(FX_W(i,j),FX(i,j)),                        &
     &                   MIN(FX(i,j),FX_med_pos),                       &
     &                   MIN(FX_med_pos,FX_W(i,j)))
!FX Negative limit
              FX_tmp=FX_EE(i,j)+alpha_limiter*(FX_E(i,j)-FX_EE(i,j))
              FX_med_neg=MAX(MIN(FX_E(i,j),FX_W(i,j)),                  &
     &                       MIN(FX_W(i,j),FX_tmp),                     &
     &                       MIN(FX_tmp,FX_E(i,j)))
              FX_neg=MAX(MIN(FX_E(i,j),FX(i,j)),                        &
     &                   MIN(FX(i,j),FX_med_neg),                       &
     &                   MIN(FX_med_neg,FX_E(i,j)))

              FX(i,j)= MAX(Huon(i,j,k),0.0_r8)*FX_pos +                 &
     &                 MIN(Huon(i,j,k),0.0_r8)*FX_neg
#    ifdef UPL_MASK
              ELSE IF (uplmask(i,j).eq.0) THEN
#    endif 
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                   (t(i-1,j,k,3,itrc)+                            &
     &                    t(i,j  ,k,3,itrc))-                           &
     &                cff1*(curv(i-1,j)*MAX(Huon(i,j,k),0.0_r8)+        &
     &                      curv(i,j  )*MIN(Huon(i,j,k),0.0_r8))
#    ifdef UPL_MASK
              END IF
#    endif
#   else !END UPWIND LIMITER
!
!  ORIGINAL UPWINDING SCHEME
!
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,3,itrc)+                               &
     &                 t(i  ,j,k,3,itrc))-                              &
     &                cff1*(curv(i-1,j)*MAX(Huon(i,j,k),0.0_r8)+        &
     &                      curv(i  ,j)*MIN(Huon(i,j,k),0.0_r8))
#   endif
#  else
              FX(i,j)=Huon(i,j,k)*0.5_r8*                               &
     &                (t(i-1,j,k,3,itrc)+                               &
     &                 t(i  ,j,k,3,itrc)-                               &
     &                 cff2*(grad(i  ,j)-                               &
     &                       grad(i-1,j)))
#  endif
            END DO
          END DO
!
          DO j=Jstrm1,Jendp2
            DO i=Istr,Iend
              FE(i,j)=t(i,j  ,k,3,itrc)-                                &
     &                t(i,j-1,k,3,itrc)
#  ifdef MASKING
              FE(i,j)=FE(i,j)*vmask(i,j)
#  endif
            END DO
          END DO
          IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jstr-1)=FE(i,Jstr)
              END DO
            END IF
          END IF
          IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
            IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
              DO i=Istr,Iend
                FE(i,Jend+2)=FE(i,Jend+1)
              END DO
            END IF
          END IF
!
          DO j=Jstr-1,Jend+1
            DO i=Istr,Iend
#  if defined TS_U3HADVECTION
              curv(i,j)=FE(i,j+1)-FE(i,j)
#  elif defined TS_A4HADVECTION
              cff=2.0_r8*FE(i,j+1)*FE(i,j)
              IF (cff.gt.eps) THEN
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              ELSE
                grad(i,j)=0.0_r8
              END IF
#  else
              grad(i,j)=0.5_r8*(FE(i,j+1)+FE(i,j))
#  endif
            END DO
          END DO
!
          cff1=1.0_r8/6.0_r8
          cff2=1.0_r8/3.0_r8
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
#  ifdef TS_U3HADVECTION
#   ifdef UPWIND_LIMITER
#    ifdef UPL_MASK
              IF (uplmask(i,j).eq.1) THEN
#    endif 
              FE(i,j)=0.5_r8*                                           &
     &                (t(i,j-1,k,3,itrc)+                               &
     &                 t(i,j  ,k,3,itrc))-                              &
     &         cff1*(curv(i,j-1)*MAX(SIGN(1.0_r8,Hvom(i,j,k)),0.0_r8)-  &
     &               curv(i,j  )*MIN(SIGN(1.0_r8,Hvom(i,j,k)),0.0_r8))
!
! apply limiters
! Compute median(a,b,c)=MAX(MIN(a,b),MIN(b,c), MIN(c,a))
! Postive limit
! FE_med_pos=median(FE_S,FE_N,FE_SS+alpha*(FE_S-FE_SS))
              FE_tmp=FE_SS(i,j)+alpha_limiter*(FE_S(i,j)-FE_SS(i,j))

              FE_med_pos=MAX(MIN(FE_S(i,j),FE_N(i,j)),                  &
     &                       MIN(FE_N(i,j),FE_tmp),                     &
     &                       MIN(FE_tmp,FE_S(i,j)))
! FE_pos=median(FE_S,FE,FE_med_pos)
              FE_pos=MAX(MIN(FE_S(i,j),FE(i,j)),                        &
                         MIN(FE(i,j),FE_med_pos),                       &
     &                   MIN(FE_med_pos,FE_S(i,j)))
! Negative limit
! FE_med_neg=median(FE_N,FE_S,FE_NN+alpha*(FE_N-FE_NN))
              FE_tmp=FE_NN(i,j)+alpha_limiter*( FE_N(i,j)-FE_NN(i,j))
              FE_med_neg=MAX(MIN(FE_N(i,j),FE_S(i,j)),                  &
     &                       MIN(FE_S(i,j),FE_tmp),                     &
     &                       MIN(FE_tmp,FE_N(i,j)))
! FX_neg=median(FE_N,FE,FE_med_neg)
              FE_neg=MAX(MIN(FE_N(i,j),FE(i,j)),                        &
     &                   MIN(FE(i,j),FE_med_neg),                       &
     &                   MIN(FE_med_neg,FE_N(i,j)))

              FE(i,j)= MAX(Hvom(i,j,k),0.0_r8)*FE_pos +                 &
     &                 MIN(Hvom(i,j,k),0.0_r8)*FE_neg
#    ifdef UPL_MASK
              ELSE IF (uplmask(i,j).eq.0) THEN
#    endif
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
&                        (t(i,j-1,k,3,itrc)+                            &
&                         t(i,j  ,k,3,itrc))-                           &
&                     cff1*(curv(i,j-1)*MAX(Hvom(i,j,k),0.0_r8)+        &
&                           curv(i,j  )*MIN(Hvom(i,j,k),0.0_r8))
#    ifdef UPL_MASK
              END IF
#    endif
#   else 
!
!  ORIGINAL UPWINDING SCHEME
!
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,3,itrc)+                               &
     &                 t(i,j  ,k,3,itrc))-                              &
     &                cff1*(curv(i,j-1)*MAX(Hvom(i,j,k),0.0_r8)+        &
     &                      curv(i,j  )*MIN(Hvom(i,j,k),0.0_r8))
#   endif
#  else
              FE(i,j)=Hvom(i,j,k)*0.5_r8*                               &
     &                (t(i,j-1,k,3,itrc)+                               &
     &                 t(i,j  ,k,3,itrc)-                               &
     &                 cff2*(grad(i,j  )-                               &
     &                       grad(i,j-1)))
#  endif
            END DO
          END DO
# endif
!
!  Apply tracers point sources to the horizontal advection terms,
!  if any.
!
# if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
#  if defined TS_MPDATA || defined TS_HSIMT
#    define I_RANGE IstrUm2,Iendp2i
#    define J_RANGE JstrVm2,Jendp2i
#  else
#    define I_RANGE Istr,Iend
#    define J_RANGE Jstr,Jend
#  endif
!
!  Initialize masks for separation of point source fluxes from advection
!
          IF (itbnt>0) THEN
            DO j=J_RANGE
              DO i=I_RANGE
                psrcMaskX(i,j,k,itbnt) = .false.
                psrcMaskE(i,j,k,itbnt) = .false.
                TBNTpsrcWRK(i,j,k,itbnt,iTBNT_psrcx) = 0.0_r8
                TBNTpsrcWRK(i,j,k,itbnt,iTBNT_psrce) = 0.0_r8
                TBNTpsrcWRK(i,j,k,itbnt,iTBNT_psrcs) = 0.0_r8
              END DO
            END DO
          END IF
# endif
          IF (LuvSrc(ng).and.ANY(LtracerSrc(:,ng))) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
# if defined TS_MPDATA || defined TS_HSIMT
                IF (((IstrUm2.le.Isrc).and.(Isrc.le.Iendp3)).and.       &
     &              ((JstrVm2.le.Jsrc).and.(Jsrc.le.Jendp2i))) THEN
# else
                IF (((Istr.le.Isrc).and.(Isrc.le.Iend+1)).and.          &
     &              ((Jstr.le.Jsrc).and.(Jsrc.le.Jend))) THEN
# endif
                  IF (LtracerSrc(itrc,ng)) THEN
                    FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                    &
     &                            SOURCES(ng)%Tsrc(is,k,itrc)
# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION && defined MASKING
                    IF (itbnt>0.and.ABS(FX(Isrc,Jsrc))>TBNTeps) THEN
                      psrcMaskX(Isrc,Jsrc,k,itbnt) = .true.
                      IF (FX(Isrc,Jsrc)<0.0_r8) THEN
                        TBNTpsrcWRK(Isrc-1,Jsrc,k,itbnt,iTBNT_psrcx) =  &
     &                   TBNTpsrcWRK(Isrc-1,Jsrc,k,itbnt,iTBNT_psrcx) - &
     &                   FX(Isrc,Jsrc) * dt(ng)
                      ELSE
                        TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrcx) =    &
     &                   TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrcx) +   &
     &                   FX(Isrc,Jsrc) * dt(ng)
                      END IF
!                      IF ((rmask(Isrc  ,Jsrc).eq.0.0_r8).and.           &
!     &                    (rmask(Isrc-1,Jsrc).eq.1.0_r8)) THEN
!                        TBNTpsrcWRK(Isrc-1,Jsrc,k,itbnt,iTBNT_psrcx) =  &
!     &                   TBNTpsrcWRK(Isrc-1,Jsrc,k,itbnt,iTBNT_psrcx) + &
!     &                   ABS(FX(Isrc,Jsrc)) * dt(ng)
!                      ELSE IF ((rmask(Isrc-1,Jsrc).eq.0.0_r8).and.      &
!     &                         (rmask(Isrc  ,Jsrc).eq.1.0_r8)) THEN
!                        TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrcx) =    &
!     &                   TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrcx) +   &
!     &                   ABS(FX(Isrc,Jsrc)) * dt(ng)
!                      END IF
                    END IF
# endif
# ifdef MASKING
                  ELSE
                    IF ((rmask(Isrc  ,Jsrc).eq.0.0_r8).and.             &
     &                  (rmask(Isrc-1,Jsrc).eq.1.0_r8)) THEN
                      FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                  &
     &                              t(Isrc-1,Jsrc,k,3,itrc)
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
                      IF (itbnt>0.and.ABS(FX(Isrc,Jsrc))>TBNTeps) THEN
                        psrcMaskX(Isrc,Jsrc,k,itbnt) = .true.
                        TBNTpsrcWRK(Isrc-1,Jsrc,k,itbnt,iTBNT_psrcx) =  &
     &                   TBNTpsrcWRK(Isrc-1,Jsrc,k,itbnt,iTBNT_psrcx) - &
     &                   FX(Isrc,Jsrc) * dt(ng)
                      END IF
#  endif
                    ELSE IF ((rmask(Isrc  ,Jsrc).eq.1.0_r8).and.        &
     &                       (rmask(Isrc-1,Jsrc).eq.0.0_r8)) THEN
                      FX(Isrc,Jsrc)=Huon(Isrc,Jsrc,k)*                  &
     &                              t(Isrc  ,Jsrc,k,3,itrc)
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
                      IF (itbnt>0.and.ABS(FX(Isrc,Jsrc))>TBNTeps) THEN
                        psrcMaskX(Isrc,Jsrc,k,itbnt) = .true.
                        TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrcx) =    &
     &                   TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrcx) +   &
     &                   FX(Isrc,Jsrc) * dt(ng)
                      END IF
#  endif
                    END IF
# endif
                  END IF
                END IF
              ELSE IF (INT(SOURCES(ng)%Dsrc(is)).eq.1) THEN
# if defined TS_MPDATA || defined TS_HSIMT
                IF (((IstrUm2.le.Isrc).and.(Isrc.le.Iendp2i)).and.      &
     &              ((JstrVm2.le.Jsrc).and.(Jsrc.le.Jendp3))) THEN
# else
                IF (((Istr.le.Isrc).and.(Isrc.le.Iend)).and.            &
     &              ((Jstr.le.Jsrc).and.(Jsrc.le.Jend+1))) THEN
# endif
                  IF (LtracerSrc(itrc,ng)) THEN
                    FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                    &
     &                            SOURCES(ng)%Tsrc(is,k,itrc)
# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION && defined MASKING
                    IF (itbnt>0.and.ABS(FE(Isrc,Jsrc))>TBNTeps) THEN
                      psrcMaskE(Isrc,Jsrc,k,itbnt) = .true.
                      IF (FE(Isrc,Jsrc)<0.0_r8) THEN
                        TBNTpsrcWRK(Isrc,Jsrc-1,k,itbnt,iTBNT_psrce) =  &
     &                   TBNTpsrcWRK(Isrc,Jsrc-1,k,itbnt,iTBNT_psrce) - &
     &                   FE(Isrc,Jsrc) * dt(ng)
                      ELSE
                        TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrce) =    &
     &                   TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrce) +   &
     &                   FE(Isrc,Jsrc) * dt(ng)
                      END IF

!                      IF ((rmask(Isrc,Jsrc  ).eq.0.0_r8).and.           &
!     &                    (rmask(Isrc,Jsrc-1).eq.1.0_r8)) THEN
!                        TBNTpsrcWRK(Isrc,Jsrc-1,k,itbnt,iTBNT_psrce) =  &
!     &                   TBNTpsrcWRK(Isrc,Jsrc-1,k,itbnt,iTBNT_psrce)   &
!     &                   +ABS(FE(Isrc,Jsrc)) * dt(ng)
!                      ELSE IF ((rmask(Isrc,Jsrc-1).eq.0.0_r8).and.      &
!     &                         (rmask(Isrc,Jsrc  ).eq.1.0_r8)) THEN
!                        TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrce) =    &
!     &                   TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrce)     &
!     &                   +ABS(FE(Isrc,Jsrc)) * dt(ng)
!                      END IF
                    END IF
# endif
# ifdef MASKING
                  ELSE
                    IF ((rmask(Isrc,Jsrc  ).eq.0.0_r8).and.             &
     &                  (rmask(Isrc,Jsrc-1).eq.1.0_r8)) THEN
                      FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                  &
     &                              t(Isrc,Jsrc-1,k,3,itrc)
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
                      IF (itbnt>0.and.ABS(FE(Isrc,Jsrc))>TBNTeps) THEN
                        psrcMaskE(Isrc,Jsrc,k,itbnt) = .true.
                        TBNTpsrcWRK(Isrc,Jsrc-1,k,itbnt,iTBNT_psrce) =  &
     &                   TBNTpsrcWRK(Isrc,Jsrc-1,k,itbnt,iTBNT_psrce)   &
     &                   -FE(Isrc,Jsrc) * dt(ng)
                      END IF
#  endif
                    ELSE IF ((rmask(Isrc,Jsrc  ).eq.1.0_r8).and.        &
     &                       (rmask(Isrc,Jsrc-1).eq.0.0_r8)) THEN
                      FE(Isrc,Jsrc)=Hvom(Isrc,Jsrc,k)*                  &
     &                              t(Isrc,Jsrc  ,k,3,itrc)
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
                      IF (itbnt>0.and.ABS(FE(Isrc,Jsrc))>TBNTeps) THEN
                        psrcMaskE(Isrc,Jsrc,k,itbnt) = .true.
                        TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrce) =    &
     &                   TBNTpsrcWRK(Isrc,Jsrc,k,itbnt,iTBNT_psrce)     &
     &                   +FE(Isrc,Jsrc) * dt(ng)
                      END IF
#  endif
                    END IF
# endif
                  END IF
                END IF
              END IF
            END DO
          END IF

# if defined NESTING && !defined ONE_WAY
!
!  If refinement grids, extract tracer horizontal advection fluxes
!  (Hz*u*T/n, Hz*v*T/m) at the grid contact boundary (physical
!  domain period) to be used in two-way nesting.
!
         IF (RefinedGrid(ng)) THEN
           DO cr=1,Ncontact
             dg=Rcontact(cr)%donor_grid
             rg=Rcontact(cr)%receiver_grid
             IF (ng.eq.rg) THEN
               CALL bry_fluxes (dg, rg, cr, iNLM, tile,                 &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          ILB, IUB, JLB, JUB,                     &
     &                          dt(ng), FX, FE,                         &
     &                          BRY_CONTACT(iwest, cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(ieast, cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(isouth,cr)%Tflux(:,k,itrc), &
     &                          BRY_CONTACT(inorth,cr)%Tflux(:,k,itrc))
             END IF
           END DO
         END IF
# endif
!
# ifdef TS_MPDATA
!  Time-step horizontal advection for intermediate diffusive tracer, Ta.
!  Advective fluxes have units of Tunits m3/s.  The new tracer has
!  units of m Tunits.
# else
!  Time-step horizontal advection term.  Advective fluxes have units
!  of Tunits m3/s.  The new tracer has units of m Tunits.
# endif
!
# ifdef TS_MPDATA
#  define I_RANGE IstrUm2,Iendp2i
#  define J_RANGE JstrVm2,Jendp2i
# else
#  define I_RANGE Istr,Iend
#  define J_RANGE Jstr,Jend
# endif
          DO j=J_RANGE
            DO i=I_RANGE
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j)-FX(i,j))
              cff2=cff*(FE(i,j+1)-FE(i,j))
              cff3=cff1+cff2
# ifdef TS_MPDATA
              Ta(i,j,k,itrc)=t(i,j,k,nnew,itrc)-cff3
# else
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff3
# endif
!              IF (i==82.and.j==20.and.k==1.and.itrc==3) THEN
!                write(222,'(4(a7,E20.10))')                             &
!     &           '   FX_W',FX(i,j)*dt(ng),'   FX_E',FX(i+1,j)*dt(ng),   &
!     &           ' Huon_W',Huon(i,j,k)   ,' Huon_E',Huon(i+1,j,k)
!              END IF
# ifdef DIAGNOSTICS_TS
#  ifdef TS_MPDATA
              Dhadv(i,j,iTxadv)=-cff1
              Dhadv(i,j,iTyadv)=-cff2
              Dhadv(i,j,iThadv)=-cff3
#  else
              DiaTwrk(i,j,k,itrc,iTxadv)=-cff1
              DiaTwrk(i,j,k,itrc,iTyadv)=-cff2
              DiaTwrk(i,j,k,itrc,iThadv)=-cff3
#  endif
# endif
# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION
              IF (itbnt>0) THEN
                IF (.not.psrcMaskX(i,j,k,itbnt)) THEN
#  ifdef TS_MPDATA
                  TBNTadh(i,j,1) = FX(i,j) * dt(ng)
#  else
                  TBNTphysWRK(i,j,k,itbnt,iTBNT_hAdx) = FX(i,j) * dt(ng)
#  endif
                END IF
                IF (.not.psrcMaskE(i,j,k,itbnt)) THEN
#  ifdef TS_MPDATA
                  TBNTadh(i,j,2) = FE(i,j) * dt(ng)
#  else
                  TBNTphysWRK(i,j,k,itbnt,iTBNT_hAde) = FE(i,j) * dt(ng)
#  endif
                END IF
              END IF
# endif


            END DO
          END DO
# if defined TS_MPDATA && (defined DIAGNOSTICS_TS || \
     (defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION))
          DO j=Jstr,Jend
            DO i=Istr,Iend
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTxadv)=Dhadv(i,j,iTxadv)
              DiaTwrk(i,j,k,itrc,iTyadv)=Dhadv(i,j,iTyadv)
              DiaTwrk(i,j,k,itrc,iThadv)=Dhadv(i,j,iThadv)
#  endif
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
              IF (itbnt>0) THEN
                IF (.not.psrcMaskX(i,j,k,itbnt)) THEN
                  TBNTphysWRK(i,j,k,itbnt,iTBNT_hAdx) = TBNTadh(i,j,1)
                END IF
                IF (.not.psrcMaskE(i,j,k,itbnt)) THEN
                  TBNTphysWRK(i,j,k,itbnt,iTBNT_hAde) = TBNTadh(i,j,2)
                END IF
              END IF
#  endif
            END DO
          END DO
# endif
        END DO K_LOOP
      END DO T_LOOP
!
!-----------------------------------------------------------------------
!  Time-step vertical advection term.
!-----------------------------------------------------------------------
!
# ifdef TS_HSIMT
      cc1=0.25_r8
      cc2=0.5_r8
      cc3=1.0_r8/12.0_r8
      epson=1.0E-12_r8
# endif
      DO j=J_RANGE
        DO itrc=1,NT(ng)
# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION
          itbnt = TBNToce2tbnt(itrc)
# endif

# if defined TS_SVADVECTION
!
!  Build conservative parabolic splines for the vertical derivatives
!  "FC" of the tracer.  Then, the interfacial "FC" values are
!  converted to vertical advective flux.
!
          DO i=Istr,Iend
#  ifdef NEUMANN
            FC(i,0)=1.5_r8*t(i,j,1,3,itrc)
            CF(i,1)=0.5_r8
#  else
            FC(i,0)=2.0_r8*t(i,j,1,3,itrc)
            CF(i,1)=1.0_r8
#  endif
          END DO
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(2.0_r8*Hz(i,j,k)+                             &
     &                    Hz(i,j,k+1)*(2.0_r8-CF(i,k)))
              CF(i,k+1)=cff*Hz(i,j,k)
              FC(i,k)=cff*(3.0_r8*(Hz(i,j,k  )*t(i,j,k+1,3,itrc)+       &
     &                             Hz(i,j,k+1)*t(i,j,k  ,3,itrc))-      &
     &                     Hz(i,j,k+1)*FC(i,k-1))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef NEUMANN
            FC(i,N(ng))=(3.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/     &
     &                  (2.0_r8-CF(i,N(ng)))
#  else
            FC(i,N(ng))=(2.0_r8*t(i,j,N(ng),3,itrc)-FC(i,N(ng)-1))/     &
     &                  (1.0_r8-CF(i,N(ng)))
#  endif
          END DO
          DO k=N(ng)-1,0,-1
            DO i=Istr,Iend
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)
              FC(i,k+1)=W(i,j,k+1)*FC(i,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,0)=0.0_r8
          END DO

# elif defined TS_A4VADVECTION
!
!  Fourth-order, Akima vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=t(i,j,k+1,3,itrc)-                                &
     &                t(i,j,k  ,3,itrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=FC(i,1)
            FC(i,N(ng))=FC(i,N(ng)-1)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=2.0_r8*FC(i,k)*FC(i,k-1)
              IF (cff.gt.eps) THEN
                CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
              ELSE
                CF(i,k)=0.0_r8
              END IF
            END DO
          END DO
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,3,itrc)+                        &
     &                        t(i,j,k+1,3,itrc)-                        &
     &                        cff1*(CF(i,k+1)-CF(i,k)))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
          END DO
# elif defined TS_C2VADVECTION
!
!  Second-order, central differences vertical advective flux.
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                0.5_r8*(t(i,j,k  ,3,itrc)+                        &
     &                        t(i,j,k+1,3,itrc))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
          END DO
# elif defined TS_MPDATA
!
!  First_order, upstream differences vertical advective flux.
!
          DO i=I_RANGE
            DO k=1,N(ng)-1
              cff1=MAX(W(i,j,k),0.0_r8)
              cff2=MIN(W(i,j,k),0.0_r8)
              FC(i,k)=cff1*t(i,j,k  ,3,itrc)+                           &
     &                cff2*t(i,j,k+1,3,itrc)
            END DO
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,N(ng))=0.0_r8
          END DO
# elif defined TS_HSIMT
          kaz(0)=0.0_r8
          kaz(N(ng))=0.0_r8
          kaz_inverse(0)=0.0_r8
          kaz_inverse(N(ng))=0.0_r8
          grad_k(0)=0.0_r8
          grad_k(N(ng))=0.0_r8
          DO i=I_RANGE
            DO k=1,N(ng)-1
              kaz(k)=1.0_r8-abs(w(i,j,k)*pm(i,j)*pn(i,j)*dt(ng)/        &
     &                         (z_r(i,j,k+1)-z_r(i,j,k)))
              kaz_inverse(k)=1.0_r8/(kaz(k))
              grad_k(k)=t(i,j,k+1,3,itrc)-t(i,j,k,3,itrc)
            END DO
            DO k=1,N(ng)-1
              IF (k.eq.1.and.W(i,j,k).ge.0.0_r8) THEN
                FC(i,k)=W(i,j,k)*t(i,j,k,3,itrc)
              ELSE IF (k.eq.N(ng)-1.and.W(i,j,k).lt.0.0_r8) THEN
                FC(i,k)=W(i,j,k)*t(i,j,k+1,3,itrc)
              ELSE
                IF (W(i,j,k).ge.0) THEN
                  IF (abs(grad_k(k)).le.epson) THEN
                    rd=0.0_r8
                    rkad=0.0_r8
                  ELSE
                    rd=grad_k(k-1)/grad_k(k)
                    rkad=kaz(k-1)*kaz_inverse(k)
!                   rkad=kaz(k-1)*kaz_inverse(k)*W(i,j,k-1)/W(i,j,k)
                  END IF
                  a1= cc1*kaz(k)+cc2-cc3*kaz_inverse(k)
                  b1=-cc1*kaz(k)+cc2+cc3*kaz_inverse(k)
                  betad=a1+b1*rd
                  sw=t(i,j,k,3,itrc)+                                   &
     &               0.5_r8*MAX(0.0_r8,                                 &
     &                          MIN(2.0_r8,2.0_r8*rd*rkad,betad))*      &
     &                grad_k(k)*kaz(k)
                ELSE
                  IF (abs(grad_k(k)).le.epson) THEN
                    ru=0.0_r8
                    rkau=0.0_r8
                  ELSE
                    ru=grad_k(k+1)/grad_k(k)
                    rkau=kaz(k+1)*kaz_inverse(k)
!                   rkau=kaz(k+1)*kaz_inverse(k)*W(i,j,k+1)/W(i,j,k)
                  END IF
                  a1= cc1*kaz(k)+cc2-cc3*kaz_inverse(k)
                  b1=-cc1*kaz(k)+cc2+cc3*kaz_inverse(k)
                  betau=a1+b1*ru
                  sw=t(i,j,k+1,3,itrc)-                                 &
     &               0.5_r8*MAX(0.0_r8,                                 &
     &                          MIN(2.0_r8,2.0_r8*ru*rkau,betau))*      &
     &               grad_k(k)*kaz(k)
                END IF 
                FC(i,k)=W(i,j,k)*sw
              END IF
            END DO
#   ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*t(i,j,1,3,itrc)
#   else
            FC(i,0)=0.0_r8
#   endif
            FC(i,N(ng))=0.0_r8
          END DO
# else
!
!  Fourth-order, central differences vertical advective flux.
!
          cff1=0.5_r8
          cff2=7.0_r8/12.0_r8
          cff3=1.0_r8/12.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=W(i,j,k)*                                         &
     &                (cff2*(t(i,j,k  ,3,itrc)+                         &
     &                       t(i,j,k+1,3,itrc))-                        &
     &                 cff3*(t(i,j,k-1,3,itrc)+                         &
     &                       t(i,j,k+2,3,itrc)))
            END DO
          END DO
          DO i=Istr,Iend
#  ifdef SED_MORPH
            FC(i,0)=W(i,j,0)*2.0_r8*                                    &
     &              (cff2*t(i,j,1,3,itrc)-                              &
     &               cff3*t(i,j,2,3,itrc))
#  else
            FC(i,0)=0.0_r8
#  endif
            FC(i,1)=W(i,j,1)*                                           &
     &              (cff1*t(i,j,1,3,itrc)+                              &
     &               cff2*t(i,j,2,3,itrc)-                              &
     &               cff3*t(i,j,3,3,itrc))
            FC(i,N(ng)-1)=W(i,j,N(ng)-1)*                               &
     &                    (cff1*t(i,j,N(ng)  ,3,itrc)+                  &
     &                     cff2*t(i,j,N(ng)-1,3,itrc)-                  &
     &                     cff3*t(i,j,N(ng)-2,3,itrc))
            FC(i,N(ng))=0.0_r8
          END DO
# endif
!
# ifdef TS_MPDATA
!  Time-step vertical advection for intermediate diffusive tracer, Ta
!  (Tunits).
# else
#  ifdef SPLINES_VDIFF
!  Time-step vertical advection term (Tunits).
#  else
!  Time-step vertical advection term (m Tunits).
#  endif
# endif
# ifdef DIAGNOSTICS_TS
!  Convert units of tracer diagnostic terms to Tunits.
# endif
!
          DO i=I_RANGE
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION
            DO k=0,N(ng)
              TBNT_FC(i,k,1) = FC(i,k)
              TBNT_FC(i,k,2) = 0.0_r8
            END DO
# endif
          END DO
!
!  Apply mass point sources (volume vertical influx), if any.
!
          IF (LwSrc(ng).and.ANY(LtracerSrc(:,ng))) THEN
            DO is=1,Nsrc(ng)
              Isrc=SOURCES(ng)%Isrc(is)
              Jsrc=SOURCES(ng)%Jsrc(is)
              IF (LtracerSrc(itrc,ng).and.                              &
# if defined TS_MPDATA || defined TS_HSIMT
     &            ((IstrUm2.le.Isrc).and.(Isrc.le.Iendp2i)).and.        &
# else
     &            ((IstrR.le.Isrc).and.(Isrc.le.IendR)).and.            &
# endif
     &            (j.eq.Jsrc)) THEN
                DO k=1,N(ng)-1
# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION
                  IF (itbnt>0) THEN
                    psrcIn = 0.5_r8*                                    &
     &                       (SOURCES(ng)%Qsrc(is,k  )*                 &
     &                        SOURCES(ng)%Tsrc(is,k  ,itrc)+            &
     &                        SOURCES(ng)%Qsrc(is,k+1)*                 &
     &                        SOURCES(ng)%Tsrc(is,k+1,itrc))
                    IF (ABS(psrcIn)>TBNTeps) THEN
                      TBNT_FC(Isrc,k,2) = psrcIn
                    ELSE
                      TBNT_FC(Isrc,k,1) = TBNT_FC(Isrc,k,1) + psrcIn
                    END IF
                  END IF
# endif
                  FC(Isrc,k)=FC(Isrc,k)+0.5_r8*                         &
     &                       (SOURCES(ng)%Qsrc(is,k  )*                 &
     &                        SOURCES(ng)%Tsrc(is,k  ,itrc)+            &
     &                        SOURCES(ng)%Qsrc(is,k+1)*                 &
     &                        SOURCES(ng)%Tsrc(is,k+1,itrc))
                END DO
              END IF
            END DO
          END IF
!
          DO k=1,N(ng)
            DO i=I_RANGE
              cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
# ifdef TS_MPDATA
              Ta(i,j,k,itrc)=(Ta(i,j,k,itrc)-cff1)*oHz(i,j,k)
#  ifdef DIAGNOSTICS_TS
              Dvadv(i,j,k,itrc)=-cff1
#  endif
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
              IF (itbnt>0) THEN
                TBNTadv(i,j,k,itbnt) = TBNT_FC(i,k-1,1) * dt(ng)
                TBNTpsrcWRK(i,j,k,itbnt,iTBNT_psrcs) =                  &
     &                                 TBNT_FC(i,k-1,2) * dt(ng)
              END IF
#  endif
# else
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
#  ifdef SPLINES_VDIFF
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*oHz(i,j,k)
#  endif
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvadv)=-cff1
              DO idiag=1,NDT
                DiaTwrk(i,j,k,itrc,idiag)=DiaTwrk(i,j,k,itrc,idiag)*    &
     &                                    oHz(i,j,k)
              END DO
#  endif
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
              IF (itbnt>0) THEN
                TBNTphysWRK(i,j,k,itbnt,iTBNT_vAds) =                   &
     &                                   TBNT_FC(i,k-1,1) * dt(ng)
                TBNTpsrcWRK(i,j,k,itbnt,iTBNT_psrcs) =                  &
     &                                   TBNT_FC(i,k-1,2) * dt(ng)
              END IF
#  endif
# endif
            END DO
          END DO
        END DO
# undef I_RANGE
# undef J_RANGE
# ifdef TS_MPDATA
      END DO
!
!-----------------------------------------------------------------------
!  Compute anti-diffusive velocities to corrected advected tracers
!  using MPDATA recursive method.  Notice that pipelined J-loop ended.
!-----------------------------------------------------------------------
!
      DO itrc=1,NT(ng)     
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
        itbnt = TBNToce2tbnt(itrc)
#  endif

        CALL mpdata_adiff_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
#  ifdef MASKING
     &                          rmask, umask, vmask,                    &
#  endif
#  ifdef WET_DRY
     &                          rmask_wet, umask_wet, vmask_wet,        &
#  endif
     &                          pm, pn, omn, om_u, on_v,                &
     &                          z_r, oHz,                               &
     &                          Huon, Hvom, W,                          &
     &                          t(:,:,:,3,itrc),                        &
     &                          Ta(:,:,:,itrc),  Ua, Va, Wa)
!
!  Compute anti-diffusive corrected advection fluxes.
!
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff1=MAX(Ua(i,j,k),0.0_r8)
              cff2=MIN(Ua(i,j,k),0.0_r8)
              FX(i,j)=(cff1*Ta(i-1,j,k,itrc)+                           &
     &                 cff2*Ta(i  ,j,k,itrc))*                          &
     &                0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff1=MAX(Va(i,j,k),0.0_r8)
              cff2=MIN(Va(i,j,k),0.0_r8)
              FE(i,j)=(cff1*Ta(i,j-1,k,itrc)+                           &
     &                 cff2*Ta(i,j  ,k,itrc))*                          &
     &                0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
            END DO
          END DO
!
!  Time-step corrected horizontal advection (Tunits m).
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              cff=dt(ng)*pm(i,j)*pn(i,j)
              cff1=cff*(FX(i+1,j)-FX(i,j))
              cff2=cff*(FE(i,j+1)-FE(i,j))
              cff3=cff1+cff2
              t(i,j,k,nnew,itrc)=Ta(i,j,k,itrc)*Hz(i,j,k)-cff3
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTxadv)=DiaTwrk(i,j,k,itrc,iTxadv)-    &
     &                                   cff1
              DiaTwrk(i,j,k,itrc,iTyadv)=DiaTwrk(i,j,k,itrc,iTyadv)-    &
     &                                   cff2
              DiaTwrk(i,j,k,itrc,iThadv)=DiaTwrk(i,j,k,itrc,iThadv)-    &
     &                                   cff3
#  endif
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
              IF (itbnt>0) THEN
                IF (.not.psrcMaskX(i,j,k,itbnt)) THEN
                  TBNTphysWRK(i,j,k,itbnt,iTBNT_hAdx) =                 &
     &                           TBNTphysWRK(i,j,k,itbnt,iTBNT_hAdx) +  &
     &                           FX(i,j) * dt(ng)
                ELSE
                  IF (FX(i,j)>0.0_r8) THEN
                    TBNTpsrcWRK(i,j,k,itbnt,iTBNT_psrcx) =              &
     &                           TBNTpsrcWRK(i,j,k,itbnt,iTBNT_psrcx) + &
     &                           FX(i,j) * dt(ng)
                  ELSE
                    TBNTpsrcWRK(i-1,j,k,itbnt,iTBNT_psrcx) =            &
     &                         TBNTpsrcWRK(i-1,j,k,itbnt,iTBNT_psrcx) - &
     &                         FX(i,j) * dt(ng)
                  END IF
                END IF
                IF (.not.psrcMaskE(i,j,k,itbnt)) THEN
                  TBNTphysWRK(i,j,k,itbnt,iTBNT_hAde) =                 &
     &                           TBNTphysWRK(i,j,k,itbnt,iTBNT_hAde) +  &
                                 FE(i,j) * dt(ng)
                ELSE
                  IF (FE(i,j)>0.0_r8) THEN
                    TBNTpsrcWRK(i,j,k,itbnt,iTBNT_psrce) =              &
     &                           TBNTpsrcWRK(i,j,k,itbnt,iTBNT_psrce) + &
     &                           FE(i,j) * dt(ng)
                  ELSE
                    TBNTpsrcWRK(i,j-1,k,itbnt,iTBNT_psrce) =            &
     &                         TBNTpsrcWRK(i,j-1,k,itbnt,iTBNT_psrce) - &
     &                         FE(i,j) * dt(ng)
                  END IF
                END IF
              END IF
#  endif
            END DO
          END DO
        END DO
!
!  Compute anti-diffusive corrected vertical advection flux.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=MAX(Wa(i,j,k),0.0_r8)
              cff2=MIN(Wa(i,j,k),0.0_r8)
              FC(i,k)=cff1*Ta(i,j,k  ,itrc)+                            &
     &                cff2*Ta(i,j,k+1,itrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Time-step corrected vertical advection (Tunits).
#  ifdef DIAGNOSTICS_TS
!  Convert units of tracer diagnostic terms to Tunits.
#  endif
!
          DO i=Istr,Iend
            CF(i,0)=dt(ng)*pm(i,j)*pn(i,j)
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=CF(i,0)*(FC(i,k)-FC(i,k-1))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-cff1
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvadv)=Dvadv(i,j,k,itrc)-             &
     &                                   cff1
              DO idiag=1,NDT
                DiaTwrk(i,j,k,itrc,idiag)=DiaTwrk(i,j,k,itrc,idiag)*    &
     &                                    oHz(i,j,k)
              END DO
#  endif
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
              IF (itbnt>0) THEN
                TBNTphysWRK(i,j,k,itbnt,iTBNT_vAds) =                   &
     &                                           TBNTadv(i,j,k,itbnt) + &
     &                                           FC(i,k-1) * dt(ng)
              END IF
#  endif
            END DO
          END DO
        END DO
      END DO
!
!  Start pipelined J-loop.
!
      DO j=Jstr,Jend
# endif /* TS_MPDATA */
!
!-----------------------------------------------------------------------
!  Time-step vertical diffusion term.
!-----------------------------------------------------------------------
!
        DO itrc=1,NT(ng)
# if defined TBNT_OUT && defined BIO_FENNEL && \
     defined BIO_SEDIMENT && defined DENITRIFICATION
          itbnt = TBNToce2tbnt(itrc)
# endif
          ltrc=MIN(NAT,itrc)

# if defined SPLINES_VDIFF && !defined TS_MPDATA
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hz(i,j,k  )-                                 &
     &                dt(ng)*Akt(i,j,k-1,ltrc)*oHz(i,j,k  )
              CF(i,k)=cff1*Hz(i,j,k+1)-                                 &
     &                dt(ng)*Akt(i,j,k+1,ltrc)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
            DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hz(i,j,k)+Hz(i,j,k+1))+                     &
     &                dt(ng)*Akt(i,j,k,ltrc)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
              DC(i,k)=cff*(t(i,j,k+1,nnew,itrc)-t(i,j,k,nnew,itrc)-     &
     &                     FC(i,k)*DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
            DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=DC(i,k)*Akt(i,j,k,ltrc)
              cff1=dt(ng)*oHz(i,j,k)*(DC(i,k)-DC(i,k-1))
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+cff1
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvdif)=DiaTwrk(i,j,k,itrc,iTvdif)+    &
     &                                   cff1
#  endif
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
              IF (itbnt>0.and.k>1) THEN
               ! NOTE: "k>1" forces bottom mixing to be ZERO
               !       => DC(i,0) = 0, i.e., no error
                TBNTphysWRK(i,j,k,itbnt,iTBNT_vDis) =                   &
     &                            TBNTphysWRK(i,j,k,itbnt,iTBNT_vDis) - &
     &                            DC(i,k-1) * omn(i,j) * dt(ng)
              END IF
#  endif
            END DO
          END DO
# else
!
!  Compute off-diagonal coefficients FC [lambda*dt*Akt/Hz] for the
!  implicit vertical diffusion terms at future time step, located
!  at horizontal RHO-points and vertical W-points.
!  Also set FC at the top and bottom levels.
!
          cff=-dt(ng)*lambda
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              cff1=1.0_r8/(z_r(i,j,k+1)-z_r(i,j,k))
              FC(i,k)=cff*cff1*Akt(i,j,k,ltrc)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,0)=0.0_r8
            FC(i,N(ng))=0.0_r8
          END DO
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tracer equation into DC.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=Hz(i,j,k)-FC(i,k)-FC(i,k-1)
              DC(i,k)=t(i,j,k,nnew,itrc)
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
#  if defined DIAGNOSTICS_TS || \
      (defined TBNT_OUT && defined BIO_FENNEL && \
       defined BIO_SEDIMENT && defined DENITRIFICATION)
             cff1=t(i,j,N(ng),nnew,itrc)*oHz(i,j,N(ng))
#  endif
             DC(i,N(ng))=(DC(i,N(ng))-FC(i,N(ng)-1)*DC(i,N(ng)-1))/     &
     &                   (BC(i,N(ng))-FC(i,N(ng)-1)*CF(i,N(ng)-1))
             t(i,j,N(ng),nnew,itrc)=DC(i,N(ng))
#  ifdef DIAGNOSTICS_TS
             DiaTwrk(i,j,N(ng),itrc,iTvdif)=                            &
     &                              DiaTwrk(i,j,N(ng),itrc,iTvdif)+     &
     &                              t(i,j,N(ng),nnew,itrc)-cff1
#  endif
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
             IF (itbnt>0) THEN
               ! vertical diffusion across bottom of surface grid cell
               vdif_bot(i) = (t(i,j,N(ng),nnew,itrc)-cff1) *            &
     &                       Hz(i,j,N(ng)) * omn(i,j)
               TBNTphysWRK(i,j,N(ng),itbnt,iTBNT_vDis) =                &
     &           TBNTphysWRK(i,j,N(ng),itbnt,iTBNT_vDis) + vdif_bot(i)
             END IF
#  endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
#  if defined DIAGNOSTICS_TS || \
      (defined TBNT_OUT && defined BIO_FENNEL && \
       defined BIO_SEDIMENT && defined DENITRIFICATION)
              cff1=t(i,j,k,nnew,itrc)*oHz(i,j,k)
#  endif
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              t(i,j,k,nnew,itrc)=DC(i,k)
#  ifdef DIAGNOSTICS_TS
              DiaTwrk(i,j,k,itrc,iTvdif)=DiaTwrk(i,j,k,itrc,iTvdif)+    &
     &                                   t(i,j,k,nnew,itrc)-cff1
#  endif
#  if defined TBNT_OUT && defined BIO_FENNEL && \
      defined BIO_SEDIMENT && defined DENITRIFICATION
              IF (itbnt>0.and.k>1) THEN
                ! vertical diffusion across bottom of grid cell calculated as
                ! change in concentration in grid cell PLUS flux at bottom of above grid cell
                ! NOTE: "k>1" forces bottom mixing to be ZERO
                !       => relative error in order of 1.0E-11 => acceptable
                vdif_bot(i) = (t(i,j,k,nnew,itrc)-cff1) *               &
                              Hz(i,j,k)*omn(i,j) + vdif_bot(i)
                TBNTphysWRK(i,j,k,itbnt,iTBNT_vDis) =                   &
     &            TBNTphysWRK(i,j,k,itbnt,iTBNT_vDis) + vdif_bot(i)
              END IF
#  endif
            END DO
          END DO
# endif
        END DO
      END DO

# if defined AGE_MEAN && defined T_PASSIVE
!
!-----------------------------------------------------------------------
!  If inert passive tracer and Mean Age, compute age concentration (even
!  inert index) forced by the right-hand-side term that is concentration
!  of an associated conservative passive tracer (odd inert index). Mean
!  Age is age concentration divided by conservative passive tracer
!  concentration. Code implements NPT/2 mean age tracer pairs.
!
!  Implemented and tested by W.G. Zhang and J. Wilkin. See following
!  reference for details.
!
!   Zhang et al. (2010): Simulation of water age and residence time in
!      the New York Bight, JPO, 40,965-982, doi:10.1175/2009JPO4249.1
!-----------------------------------------------------------------------
!
      DO itrc=1,NPT,2
        iage=inert(itrc+1)                     ! even inert tracer index
        DO k=1,N(ng)
          DO j=Jstr,Jend
            DO i=Istr,Iend
              t(i,j,k,nnew,iage)=t(i,j,k,nnew,iage)+                    &
     &                           dt(ng)*                                &
#  ifdef TS_MPDATA
     &                           t(i,j,k,nnew,inert(itrc))
#  else
     &                           t(i,j,k,3,inert(itrc))
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions and, if appropriate, nudge
!  to tracer data and apply Land/Sea mask.
!-----------------------------------------------------------------------
!
!  Initialize tracer counter index. The "tclm" array is only allocated
!  to the NTCLM fields that need to be processed. This is done to
!  reduce memory.
!
      ic=0
!
      DO itrc=1,NT(ng)
!
!  Set compact reduced memory tracer index for nudging coefficients and
!  climatology arrays.
!
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
        END IF
!
!  Set lateral boundary conditions.
!
        CALL t3dbc_tile (ng, tile, itrc, ic,                            &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp, nnew,                                    &
     &                   t)
!
!  Nudge towards tracer climatology.
!
        IF (LtracerCLM(itrc,ng).and.LnudgeTCLM(itrc,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+                  &
     &                             dt(ng)*                              &
     &                             CLIMA(ng)%Tnudgcof(i,j,k,ic)*        &
     &                             (CLIMA(ng)%tclm(i,j,k,ic)-           &
     &                              t(i,j,k,nnew,itrc))
              END DO
            END DO
          END DO
        END IF

# ifdef MASKING
!
!  Apply Land/Sea mask.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*rmask(i,j)
            END DO
          END DO
        END DO
# endif
# ifdef DIAGNOSTICS_TS
!
!  Compute time-rate-of-change diagnostic term.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              DiaTwrk(i,j,k,itrc,iTrate)=t(i,j,k,nnew,itrc)-            &
     &                                   DiaTwrk(i,j,k,itrc,iTrate)
!!            DiaTwrk(i,j,k,itrc,iTrate)=t(i,j,k,nnew,itrc)-            &
!!   &                                   t(i,j,k,nstp,itrc)
            END DO
          END DO
        END DO
# endif
!
!  Apply periodic boundary conditions.
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            t(:,:,:,nnew,itrc))
        END IF
      END DO
# ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    t(:,:,:,nnew,:))
# endif
# if defined FLOATS && defined FLOAT_VWALK
!
!-----------------------------------------------------------------------
!  Compute vertical gradient in vertical T-diffusion coefficient for
!  floats random walk.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          DO k=1,N(ng)
            dAktdz(i,j,k)=(Akt(i,j,k,1)-Akt(i,j,k-1,1))/Hz(i,j,k)
          END DO
        END DO
      END DO
!
!  Apply periodic boundary conditions.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          dAktdz)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange3d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    dAktdz)
#  endif
# endif

      RETURN
      END SUBROUTINE step3d_t_tile
#endif
      END MODULE step3d_t_mod
