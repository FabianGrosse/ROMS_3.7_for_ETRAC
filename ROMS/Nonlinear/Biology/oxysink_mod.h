!
!svn $Id: oxysink_mod.h 645 2013-01-22 23:21:54Z arango $
! Modified from fennel_mod.h
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for the simple oxygen sink model:                        !
!                                                                      !
!   BioIter    Maximum number of iterations to achieve convergence     !
!              of the nonlinear solution.                              !
!   OxyCCR   Constant respiration rate [mmol_O2/(m3 day)].             !
!   OxyCCR_hmin Minimum depth for decrease in respiration rate [m].    !
!   OxyCCR_hmax Depth at which respiration rate goes to zero [m].      !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:)  ! Biological tracers
      integer :: iOxyg                  ! Dissolved oxygen concentration

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Biological 2D diagnostic variable IDs.
!
      integer, allocatable :: iDbio2(:)       ! 2D biological terms

      integer  :: iO2fx                       ! air-sea O2 flux
!
!  Biological 3D diagnostic variable IDs.
!
      integer, allocatable :: iDbio3(:)       ! 3D biological terms

      integer  :: iPPro = 1                   ! primary productivity
#endif
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)

      real(r8), allocatable :: OxyCCR(:)             ! mmol/(m-2 day)
      real(r8), allocatable :: OxyCCR_hmin(:)        ! m
      real(r8), allocatable :: OxyCCR_hmax(:)        ! m

      CONTAINS

      SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic
!
!-----------------------------------------------------------------------
!  Determine number of biological tracers.
!-----------------------------------------------------------------------
!
      NBT=1

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
!  Set sources and sinks biology diagnostic parameters.
!-----------------------------------------------------------------------
!
!  Set number of diagnostics terms.
!
      NDbio3d=0
      NDbio2d=1
!
!  Initialize biology diagnostic indices.
!
      ic=0
      iO2fx=ic+1
#endif
!
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(OxyCCR)) THEN
        allocate ( OxyCCR(Ngrids) )
      END IF
      IF (.not.allocated(OxyCCR_hmin)) THEN
        allocate ( OxyCCR_hmin(Ngrids) )
      END IF
      IF (.not.allocated(OxyCCR_hmax)) THEN
        allocate ( OxyCCR_hmax(Ngrids) )
      END IF
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Allocate biological diagnostics vectors
!
      IF (.not.allocated(iDbio2)) THEN
        allocate ( iDbio2(NDbio2d) )
      END IF
      IF (.not.allocated(iDbio3)) THEN
        allocate ( iDbio3(NDbio3d) )
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iOxyg=ic+1
      ic=ic+1

      RETURN
      END SUBROUTINE initialize_biology
