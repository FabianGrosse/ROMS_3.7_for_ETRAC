#include "cppdefs.h"
      MODULE mod_tbnt
#if defined TBNT_OUT && defined BIO_FENNEL &&          \
    defined BIO_SEDIMENT && defined DENITRIFICATION && \
    defined SOLVE3D
!
!svn $Id: mod_tbnt.F 8 2007-02-06 19:00:29Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  TBNT fields for output.                                             !
!                                                                      !
!  TBNToce2tbnt TBNT tracer pointer from OCEAN array to TBNT array     !
!  TBNTtbnt2oce TBNT tracer pointer from TBNT array to OCEAN array     !
!  TBNTflx2     TBNT 2D fluxes (sediment/air-sea)                      !
!  TBNTflx3     TBNT 3D biochemical fluxes                             !
!  TBNTsink     TBNT 3D sinking fluxes (for tracers with wSink>0)      !
!  TBNTphys     TBNT 3D physical exchange fluxes (advection/diffusion) !
!  TBNTpsrc     TBNT 3D point source fluxes                             !
!                                                                      !
!  TBNT fields work arrays.                                            !
!                                                                      !
!  TBNTphysWRK  TBNT work array for 3D physical exchange fluxes        !
!  TBNTpsrcWRK  TBNT work array for 3D point source fluxes                    !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_TBNT
          
          integer, allocatable :: TBNToce2tbnt(:)
          integer, allocatable :: TBNTtbnt2oce(:)
          real(r8), pointer :: TBNTgvol(:,:,:)        
          real(r8), pointer :: TBNTflx2(:,:,:)
          real(r8), pointer :: TBNTflx3(:,:,:,:)
          real(r8), pointer :: TBNTsink(:,:,:,:)
          real(r8), pointer :: TBNTphys(:,:,:,:,:)
          real(r8), pointer :: TBNTphysWRK(:,:,:,:,:)
          real(r8), pointer :: TBNTpsrc(:,:,:,:,:)
          real(r8), pointer :: TBNTpsrcWRK(:,:,:,:,:)

        END TYPE T_TBNT

        TYPE (T_TBNT), allocatable :: TBNT(:)

      CONTAINS

      SUBROUTINE allocate_tbnt (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1 ) allocate ( TBNT(Ngrids) )
!
      allocate ( TBNT(ng) % TBNToce2tbnt(NT(ng)) )
      allocate ( TBNT(ng) % TBNTtbnt2oce(NTBNTTRC(ng)) )
      allocate ( TBNT(ng) % TBNTgvol(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( TBNT(ng) % TBNTflx2(LBi:UBi,LBj:UBj,nTBNTflx2) )
      allocate ( TBNT(ng) % TBNTflx3(LBi:UBi,LBj:UBj,N(ng),nTBNTflx3) )
      allocate ( TBNT(ng) % TBNTsink(LBi:UBi,LBj:UBj,N(ng),nTBNTsink) )
      allocate ( TBNT(ng) % TBNTphys(LBi:UBi,LBj:UBj,N(ng),             &
                                     NTBNTTRC(ng),nTBNTphys) )
      allocate ( TBNT(ng) % TBNTphysWRK(LBi:UBi,LBj:UBj,N(ng),          &
     &                                  NTBNTTRC(ng),nTBNTphys) )
      allocate ( TBNT(ng) % TBNTpsrc(LBi:UBi,LBj:UBj,N(ng),             &
     &                               NTBNTTRC(ng),nTBNTpsrc) )
      allocate ( TBNT(ng) % TBNTpsrcWRK(LBi:UBi,LBj:UBj,N(ng),          &
     &                                  NTBNTTRC(ng),nTBNTpsrc) )

      RETURN
      END SUBROUTINE allocate_tbnt

      SUBROUTINE initialize_tbnt (ng, tile)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_biology
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: itbnt, itrc, i, j, k

      real(r8), parameter :: IniVal = 0.0_r8
!
# include "set_bounds.h"
!
!  Set array initialization range.
!
# ifdef DISTRIBUTE
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
# else
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      ! tracer pointer from OCEAN array to TBNT array and vice versa
      TBNT(ng) % TBNToce2tbnt(:) = 0
      TBNT(ng) % TBNTtbnt2oce(:) = 0
      DO itbnt = 1,NTBNTTRC(ng)
        SELECT CASE ( itbnt )
          CASE ( iTBNT_NO3 )
             TBNT(ng) % TBNToce2tbnt(iNO3_) = itbnt
             TBNT(ng) % TBNTtbnt2oce(itbnt) = iNO3_
          CASE ( iTBNT_NH4 )
             TBNT(ng) % TBNToce2tbnt(iNH4_) = itbnt
             TBNT(ng) % TBNTtbnt2oce(itbnt) = iNH4_
          CASE ( iTBNT_PhyN )
             TBNT(ng) % TBNToce2tbnt(iPhyt) = itbnt
             TBNT(ng) % TBNTtbnt2oce(itbnt) = iPhyt
          CASE ( iTBNT_ZooN )
             TBNT(ng) % TBNToce2tbnt(iZoop) = itbnt
             TBNT(ng) % TBNTtbnt2oce(itbnt) = iZoop
          CASE ( iTBNT_SDeN )
             TBNT(ng) % TBNToce2tbnt(iSDeN) = itbnt
             TBNT(ng) % TBNTtbnt2oce(itbnt) = iSDeN
          CASE ( iTBNT_LDeN )
             TBNT(ng) % TBNToce2tbnt(iLDeN) = itbnt
             TBNT(ng) % TBNTtbnt2oce(itbnt) = iLDeN
# ifdef RIVER_PON
          CASE ( iTBNT_RDeN )
             TBNT(ng) % TBNToce2tbnt(iRDeN) = itbnt
             TBNT(ng) % TBNTtbnt2oce(itbnt) = iRDeN
# endif
        END SELECT
      END DO

      ! initalize TBNT output fields
      DO j=Jmin,Jmax
        ! grid cell volume
        DO k=1,N(ng)
          DO i=Imin,Imax
            TBNT(ng) % TBNTgvol(i,j,k) = IniVal
          END DO
        END DO
        ! 2D biochemical fluxes
        DO itbnt=1,nTBNTflx2
          DO i=Imin,Imax
            TBNT(ng) % TBNTflx2(i,j,itbnt) = IniVal
          END DO
        END DO
        ! 3D biochemical fluxes
        DO itbnt=1,nTBNTflx3
          DO k=1,N(ng)
            DO i=Imin,Imax
              TBNT(ng) % TBNTflx3(i,j,k,itbnt) = IniVal
            END DO
          END DO
        END DO
        ! sinking fluxes
        DO itbnt=1,nTBNTsink
          DO k=1,N(ng)
            DO i=Imin,Imax
              TBNT(ng) % TBNTsink(i,j,k,itbnt) = IniVal
            END DO
          END DO
        END DO
        ! physical fluxes
        DO itrc=1,NTBNTTRC(ng)
          ! 3D physical fluxes
          DO itbnt=1,nTBNTphys
            DO k=1,N(ng)
              DO i=Imin,Imax
                TBNT(ng) % TBNTphys(i,j,k,itrc,itbnt) = IniVal
                TBNT(ng) % TBNTphysWRK(i,j,k,itrc,itbnt) = IniVal
              END DO
            END DO
          END DO
          ! point source fluxes
          DO itbnt=1,nTBNTpsrc
            DO k=1,N(ng)
              DO i=Imin,Imax
                TBNT(ng) % TBNTpsrc(i,j,k,itrc,itbnt) = IniVal
                TBNT(ng) % TBNTpsrcWRK(i,j,k,itrc,itbnt) = IniVal
              END DO
            END DO
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE initialize_tbnt
#endif
      END MODULE mod_tbnt
