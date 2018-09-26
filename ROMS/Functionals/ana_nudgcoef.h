      SUBROUTINE ana_nudgcoef (ng, tile, model)
!
!! svn $Id: ana_nudgcoef.h 830 2017-01-24 21:21:11Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine sets spatially varying nudging coefficients time-      !
!  scales (1/s). They are used for nuding to climatology in the        !
!  governing equations.                                                !
!                                                                      !
!  It is HIGHLY recommended to write these nudging coefficients into   !
!  input NetCDF NUDNAME instead of using analytical expressions        !
!  below.  It is very easy to introduce parallel bugs.  Also, Users    !
!  can plot their spatial distribution and fine tune their values      !
!  during the pre-proccessing stage for a particular application.      !
!                                                                      !
!  REMARK:  Nudging of free-surface in the vertically integrated       !
!  ======   continuity equation is NOT allowed because it VIOLATES     !
!  mass/volume conservation. If such nudging effects are required,     !
!  it needs to be specified on the momentum equations for (u,v)        !
!  and/or (ubar,vbar). If done on (u,v) only, its effects enter        !
!  the 2D momentum equations via the residual vertically integrated    !
!  forcing term.                                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
#include "tile.h"
!
      CALL ana_nudgcoef_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(16)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_nudgcoef
!
!***********************************************************************
      SUBROUTINE ana_nudgcoef_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_clima
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_collect
      USE mp_exchange_mod, ONLY : mp_exchange2d
# ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE mp_exchange_mod, ONLY : mp_exchange4d
# endif
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Iwrk, i, itrc, j, k
#if defined TXLA_MODEL || defined MCH_MODEL
      integer :: ncells
#endif

      real(r8) :: cff1, cff2, cff3

      real(r8), parameter :: IniVal = 0.0_r8

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set up nudging towards data time-scale coefficients (1/s).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          wrk(i,j)=0.0_r8
#if defined TXLA_MODEL || defined MCH_MODEL
          wrk_south(i,j)=0.0_r8
          wrk_east(i,j)=0.0_r8
          wrk_west(i,j)=0.0_r8
#endif
        END DO
      END DO

#if defined TXLA_MODEL || defined MCH_MODEL
!
!  Set tracer nudging coefficients in the southern and western edges following
!  the same weighting used by Barth et al. 2008, except relaxing over cff3 grid 
!  points away from boundary instead of 10.  Equation looks like as follows:
!
!    dT/dt = ...+ c (T - T_OCGM)  where        T = nested ROMS tracer
!                                         T_OCGM = tracer from outer model
!
!    Note that this formulation is designed for a relaxation zone of 46 grid cells
!
!    c_j = ( (1 / dt) * ( exp ( - (j+10) / cff3 ) ) ) - cff1 for j <= ncells
!    c_j =                0                 otherwise
!
!    Here, j = grid index perpendicular to nesting boundary (j=1 @ boundary)
!    tau = time step, and ncells = the maximum number of relaxation grid points. 
!    The term cff1 is subtracted from the function to allow for the function to
!    equal zero at 46 grid cells in from the boundary
!
      cff1=1.0_r8/(0.3333_r8*86400.0_r8)
      ncells=6

      ! southern boundary
      DO j=JstrT,MIN(JendT,ncells)
        DO i=IstrT,IendT
           wrk_south(i,j)=cff1*(2.0_r8+2.0_r8*                          &
     &      cos(pi*REAL(j,r8)/REAL(ncells*3+1,r8)+2.0_r8/3.0_r8*pi))
        END DO
      END DO

      ! eastern boundary
      DO j=JstrT,JendT
        DO i=MAX(Lm(ng)+1-ncells,IstrT),IendT
           wrk_east(i,j)=cff1*(2.0_r8+2.0_r8*                           &
     &      cos(pi*REAL(Lm(ng)-i+1,r8)/REAL(ncells*3+1,r8)+             &
     &      2.0_r8/3.0_r8*pi))
        END DO
      END DO

      ! western boundary
      DO j=JstrT,JendT
        DO i=IstrT,MIN(IendT,ncells)
           wrk_west(i,j)=cff1*(2.0_r8+2.0_r8*                           &
     &      cos(pi*REAL(i,r8)/REAL(ncells*3+1,r8)+2.0_r8/3.0_r8*pi))
        END DO
      END DO

      DO j=JstrT,JendT
        DO i=IstrT,IendT
          wrk(i,j) = MAX(wrk_south(i,j),wrk_west(i,j),wrk_east(i,j))
        END DO
      END DO

      IF (ANY(LnudgeTCLM(:,ng))) THEN
        DO itrc=1,NTCLM(ng)
          DO k=1,N(ng)
            DO j=JstrT,JendT
              DO i=IstrT,IendT
                CLIMA(ng)%Tnudgcof(i,j,k,itrc)=wrk(i,j)
              END DO
            END DO
          END DO
        END DO
      END IF

      IF (LnudgeM3CLM(ng)) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              CLIMA(ng)%M3nudgcof(i,j,k)=wrk(i,j)
            END DO
          END DO
        END DO
      END IF

      IF (LnudgeM2CLM(ng)) THEN
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            CLIMA(ng)%M2nudgcof(i,j)=wrk(i,j)
          END DO
        END DO
      END IF

#elif defined ADRIA02
!
!  Set tracer nudging coefficients in the southern edges from a three
!  days time scale at the boundary point to decrease linearly to 30 days
!  six grids points away from the boundary.
!
      cff1=1.0_r8/(3.0_r8*86400.0_r8)
      cff2=1.0_r8/(30.0_r8*86400.0_r8)
      DO j=JstrT,MIN(6,JendT)
        DO i=IstrT,IendT
          wrk(i,j)=cff2+REAL(6-j,r8)*(cff1-cff2)/6.0_r8
        END DO
      END DO

      IF (ANY(LnudgeTCLM(:,ng))) THEN
        DO itrc=1,NTCLM(ng)
          DO k=1,N(ng)
            DO j=JstrT,JendT
              DO i=IstrT,IendT
                CLIMA(ng)%Tnudgcof(i,j,k,itrc)=wrk(i,j)
              END DO
            END DO
          END DO
        END DO
      END IF
 
#elif defined SABGOM
!
!  Set tracer nudging coefficients in the southern edges from a 1
!  days time scale at the boundary point to decrease linearly to 20 days
!  30 grids points away from the boundary.
!
      cff1=1.0_r8/(1.0_r8*86400.0_r8)
      cff2=1.0_r8/(20.0_r8*86400.0_r8)
      cff3=1.0_r8/(60.0_r8*86400.0_r8)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
           wrk(i,j)=cff2
        END DO
      END DO

!south boundary

      DO i=IstrT,IendT
        DO j=JstrT,MIN(30,JendT)
          wrk(i,j)=cff2+REAL(30-j,r8)*(cff1-cff2)/30.0_r8
        END DO
      END DO

!east boundary
       DO i=MAX(IstrT,Lm(ng)-30),IendR
        DO j=JstrT,JendT
          wrk(i,j)=cff1+REAL(Lm(ng)-i,r8)*(cff2-cff1)/30.0_r8
        END DO
      END DO

       DO k=1,N(ng)
         DO j=JstrT,JendT
           DO i=IstrT,IendT
             CLIMA(ng)%Tnudgcof(i,j,k,itemp)=wrk(i,j)
             CLIMA(ng)%Tnudgcof(i,j,k,isalt)=wrk(i,j)
           END DO
          END DO
        END DO
      
!      DO k=1,N(ng) 
!        DO j=JstrT,JendT
!          DO i=IstrT,IendT
!            CLIMA(ng)%M3nudgcof(i,j,k)=wrk(i,j)
!          END DO
!        END DO
!      END DO       

#elif defined DAMEE_4
!
!  Set tracer nudging coefficients in the southern and northern edges
!  from a five days time scale at the boundary point to decrease
!  linearly to 60 days seven grids points away from the boundary.
!
      cff1=1.0_r8/(5.0_r8*86400.0_r8)
      cff2=1.0_r8/(60.0_r8*86400.0_r8)
      cff3=(7.0_r8*cff1-cff2)/6.0_r8
      DO j=JstrT,MIN(8,JendT)
        DO i=IstrT,IendT
          wrk(i,j)=cff2+REAL(8-j,r8)*(cff1-cff2)/7.0_r8
        END DO
      END DO

      DO j=MAX(JstrT,Mm(ng)-7),JendT
        DO i=IstrT,IendT
          wrk(i,j)=cff1+REAL(Mm(ng)-j,r8)*(cff2-cff1)/7.0_r8
        END DO
      END DO

      DO j=MAX(JstrT,74),MIN(80,JendT)
        DO i=MAX(IstrT,102),MIN(108,IendT)
          cff1=SQRT(REAL((i-109)*(i-109)+(j-77)*(j-77),r8))
          wrk(i,j)=MAX(0.0_r8,(cff3+cff1*(cff2-cff3)/6.0_r8))
        END DO
      END DO

      IF (ANY(LnudgeTCLM(:,ng))) THEN
        DO itrc=1,NTCLM(ng)
          DO k=1,N(ng)
            DO j=JstrT,JendT
              DO i=IstrT,IendT
                CLIMA(ng)%Tnudgcof(i,j,k,itrc)=wrk(i,j)
              END DO
            END DO
          END DO
        END DO
      END IF

#else
!
!  Default nudging coefficients.  Set nudging coefficients uniformly to
!  the values specified in the standard input file.
!
      IF (LnudgeM2CLM(ng)) THEN
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            CLIMA(ng)%M2nudgcof(i,j)=M2nudg(ng)
          END DO
        END DO
      END IF

# ifdef SOLVE3D
!
      IF (LnudgeM3CLM(ng)) THEN
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              CLIMA(ng)%M3nudgcof(i,j,k)=M3nudg(ng)
            END DO
          END DO
        END DO
      END IF
!
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        DO itrc=1,NTCLM(ng)
          DO k=1,N(ng)
            DO j=JstrT,JendT
              DO i=IstrT,IendT
                CLIMA(ng)%Tnudgcof(i,j,k,itrc)=Tnudg(itrc,ng)
              END DO
            END DO
          END DO
        END DO
      END IF
# endif
#endif
#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Exchage nudging coefficients information.
!-----------------------------------------------------------------------
!
      IF (LnudgeM2CLM(ng)) THEN
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      CLIMA(ng)%M2nudgcof)
      END IF

# ifdef SOLVE3D
!
      IF (LnudgeM3CLM(ng)) THEN
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      CLIMA(ng)%M3nudgcof)
      END IF
!
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        CALL mp_exchange4d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng), 1, NTCLM(ng), &
     &                      NghostPoints, .FALSE., .FALSE.,             &
     &                      CLIMA(ng)%Tnudgcof)
      END IF
# endif
#endif

      RETURN
      END SUBROUTINE ana_nudgcoef_tile