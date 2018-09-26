#include "cppdefs.h"
#if defined TBNT_OUT && defined BIO_FENNEL && \
    defined BIO_SEDIMENT && defined DENITRIFICATION && \
    defined SOLVE3D
      SUBROUTINE set_tbnt (ng, tile)
!
!svn $Id: set_tbnt.F 26 2007-04-11 03:51:31Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-integrated     !
!  fluxes and snapshopt state variable fields for TBNT application.    !
!  Due to synchronization, the TBNT fields are computed in delayed     !
!  mode. All values are accumulated at the beginning of the next       !
!  time-step.                                                          !
!                                                                      !
!  NOTE: written by Fabian Grosse based on set_diags (original ROMS3)  !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
!
      implicit none
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
      CALL wclock_on (ng, iNLM, 5)
# endif
      CALL set_tbnt_tile (ng, tile,                                    &
     &                    LBi, UBi, LBj, UBj,                          &
     &                    IminS, ImaxS, JminS, JmaxS,                  &
     &                    kstp(ng), nrhs(ng))
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_tbnt
!
!***********************************************************************
      SUBROUTINE set_tbnt_tile (ng, tile,                              &
     &                          LBi, UBi, LBj, UBj,                    &
     &                          IminS, ImaxS, JminS, JmaxS,            &
     &                          kout, nrhs)
!***********************************************************************
!
      USE mod_param
      USE mod_tbnt
      USE mod_grid
      USE mod_scalars
!
      USE bc_2d_mod
      USE bc_3d_mod
# ifdef WET_DRY
      USE exchange_2d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d,        &
     &                            mp_exchange4d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kout, nrhs
!
!  Local variable declarations.
!
      integer :: i, it, j, k
      integer :: itbnt

      real(r8) :: fac

      real(r8) :: rfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: ufac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: vfac(IminS:ImaxS,JminS:JmaxS)

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-integrating window is zero.
!-----------------------------------------------------------------------
!
      IF (nTBNT(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
! Initialize time-integrated TBNT arrays when appropriate
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsTBNT(ng)).and.                                &
     &     (MOD(iic(ng)-1,nTBNT(ng)).eq.1)).or.                         &
     &    ((iic(ng).ge.ntsTBNT(ng)).and.(nTBNT(ng).eq.1)).or.           &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
        DO it=1,NTBNTTRC(ng)
          DO itbnt=1,nTBNTphys
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  TBNT(ng)%TBNTphys(i,j,k,it,itbnt)=                    &
#   ifdef WET_DRY
     &                     GRID(ng)%rmask_full(i,j)*                    &
#   endif
     &                     TBNT(ng)%TBNTphysWRK(i,j,k,it,itbnt)
                END DO
              END DO
            END DO
          END DO
          DO itbnt=1,nTBNTpsrc
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  TBNT(ng)%TBNTpsrc(i,j,k,it,itbnt)=                    &
#   ifdef WET_DRY
     &                     GRID(ng)%rmask_full(i,j)*                    &
#   endif
     &                     TBNT(ng)%TBNTpsrcWRK(i,j,k,it,itbnt)
                END DO
              END DO
            END DO
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Accumulate time-integrated fields.
!  (only physical and point source fluxes; biological fluxes are updated
!   in biology module)
!-----------------------------------------------------------------------
!
      ELSEIF ((iic(ng).gt.ntsTBNT(ng)).and.                             &
              (iic(ng).gt.ntstart(ng))) THEN
        DO it=1,NTBNTTRC(ng)
          DO itbnt=1,nTBNTphys
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  TBNT(ng)%TBNTphys(i,j,k,it,itbnt)=                    &
     &                     TBNT(ng)%TBNTphys(i,j,k,it,itbnt)+           &
#   ifdef WET_DRY
     &                     GRID(ng)%rmask_full(i,j)*                    &
#   endif
     &                     TBNT(ng)%TBNTphysWRK(i,j,k,it,itbnt)
                END DO
              END DO
            END DO
          END DO
          DO itbnt=1,nTBNTpsrc
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  TBNT(ng)%TBNTpsrc(i,j,k,it,itbnt)=                    &
     &                     TBNT(ng)%TBNTpsrc(i,j,k,it,itbnt)+           &
#   ifdef WET_DRY
     &                     GRID(ng)%rmask_full(i,j)*                    &
#   endif
     &                     TBNT(ng)%TBNTpsrcWRK(i,j,k,it,itbnt)
                END DO
              END DO
            END DO
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions and land-mask
!  for output purposes.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsTBNT(ng)).and.                                &
     &     (MOD(iic(ng)-1,nTBNT(ng)).eq.0).and.                         &
     &     ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))).or.          &
          ((iic(ng).ge.ntsTBNT(ng)).and.(nTBNT(ng).eq.1))) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (nTBNT(ng).eq.1) THEN
            TBNTtime(ng)=time(ng)
          ELSE
            TBNTtime(ng)=TBNTtime(ng)+REAL(nTBNT(ng),r8)*dt(ng)
          END IF
        END IF
!
!  Physical fluxes
!
        DO itbnt=1,nTBNTphys
          DO it=1,NTBNTTRC(ng)
            CALL bc_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        TBNT(ng)%TBNTphys(:,:,:,it,itbnt))
          END DO
# ifdef DISTRIBUTE
          CALL mp_exchange4d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        1, NTBNTTRC(ng), NghostPoints,            &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        TBNT(ng)%TBNTphys(:,:,:,:,itbnt))
# endif
        END DO
!
!  Point source fluxes
!
        DO itbnt=1,nTBNTpsrc
          DO it=1,NTBNTTRC(ng)
            CALL bc_r3d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        TBNT(ng)%TBNTpsrc(:,:,:,it,itbnt))
          END DO
#  ifdef DISTRIBUTE
          CALL mp_exchange4d (ng, tile, iNLM, 1,                        &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        1, NTBNTTRC(ng), NghostPoints,            &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        TBNT(ng)%TBNTpsrc(:,:,:,:,itbnt))
#  endif
        END DO
!
!  2D biological fluxes
!        
        DO itbnt=1,nTBNTflx2
          CALL bc_r2d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      TBNT(ng)%TBNTflx2(:,:,itbnt))
        END DO
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, iNLM, 1,                          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      1, nTBNTflx2, NghostPoints,                 &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      TBNT(ng)%TBNTflx2)
#  endif
!
!  3D biological fluxes
!        
        DO itbnt=1,nTBNTflx3
          CALL bc_r3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      TBNT(ng)%TBNTflx3(:,:,:,itbnt))
        END DO
#  ifdef DISTRIBUTE
        CALL mp_exchange4d (ng, tile, iNLM, 1,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      1, nTBNTflx3, NghostPoints,                 &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      TBNT(ng)%TBNTflx3)
#  endif
!
!  Sinking fluxes
!        
        DO itbnt=1,nTBNTsink
          CALL bc_r3d_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      TBNT(ng)%TBNTsink(:,:,:,itbnt))
        END DO
#  ifdef DISTRIBUTE
        CALL mp_exchange4d (ng, tile, iNLM, 1,                          &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      1, nTBNTsink, NghostPoints,                 &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      TBNT(ng)%TBNTsink)
#  endif
      END IF
      
      RETURN
      END SUBROUTINE set_tbnt_tile
#else
      SUBROUTINE set_tbnt
      RETURN
      END SUBROUTINE set_tbnt
#endif
