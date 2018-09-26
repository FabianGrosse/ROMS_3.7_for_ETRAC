!
!svn $Id: hydrocarbon_mod.h 645 2013-01-22 23:21:54Z arango $
! Modified from fennel_mod.h
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for the hydrocarbon model:                               !
!                                                                      !
!   BioIter    Maximum number of iterations to achieve convergence     !
!              of the nonlinear solution.                              !
!   HydroCDR   Hydrocarbon decay rate [day-1].                         !
!   Oxy2HC     ratio of oxygen consumption to hydroC decay             !
!               [mmol_O2/(mmol_hydroC)].                               ! 
!   nday_startoil  starting day of oil spill                           !
!   nday_endoil    ending day of oil spill                             !
!   iloc_oil       longitude (i) index of the model grid where oil     !
!                  spill is located                                    !
!   jloc_oil       latitude (j) index of the model grid where oil      !
!                  spill is located                                    !
!   kloc_oil       vertical (k) index of the model grid where oil      !
!                  spill is located                                    !
!   density_oil    density of the oil added to the grid cell at        !
!                  each time step [kg/m3]                              !
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
      integer :: iHydroC                ! Hydrocarbon concentration      

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

      real(r8), allocatable :: HydroCDR(:)             ! day-1
      real(r8), allocatable :: Oxy2HC(:)               ! mmol_O2/mmol_hydroC
      real(r8), allocatable :: nday_startoil(:)         
      real(r8), allocatable :: nday_endoil(:)
      real(r8), allocatable :: iloc_oil(:)
      real(r8), allocatable :: jloc_oil(:)
      real(r8), allocatable :: kloc_oil(:)
      real(r8), allocatable :: density_oil(:)
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
      NBT=2

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
      IF (.not.allocated(HydroCDR)) THEN
        allocate ( HydroCDR(Ngrids) )
      END IF
      IF (.not.allocated(Oxy2HC)) THEN
        allocate ( Oxy2HC(Ngrids) )
      END IF
      IF (.not.allocated(nday_startoil)) THEN
        allocate ( nday_startoil(Ngrids) )
      END IF
      IF (.not.allocated(nday_endoil)) THEN
        allocate ( nday_endoil(Ngrids) )
      END IF
      IF (.not.allocated(iloc_oil)) THEN
        allocate ( iloc_oil(Ngrids) )
      END IF
      IF (.not.allocated(jloc_oil)) THEN
        allocate ( jloc_oil(Ngrids) )
      END IF 
      IF (.not.allocated(kloc_oil)) THEN
        allocate ( kloc_oil(Ngrids) )
      END IF 
      IF (.not.allocated(density_oil)) THEN
        allocate ( density_oil(Ngrids) )
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
      iHydroC=ic+2
      ic=ic+2

      RETURN
      END SUBROUTINE initialize_biology
