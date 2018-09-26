#include "cppdefs.h"
#if defined TBNT_OUT && defined BIO_FENNEL && \
    defined BIO_SEDIMENT && defined DENITRIFICATION && \
    defined SOLVE3D
      SUBROUTINE wrt_tbnt (ng)
!
!svn $Id: wrt_tbnt.F 75 2007-06-07 01:01:18Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes TBNT variables to TBNT NetCDF file           !
!   - fluxes: time-integrated over output interval                     !
!   - variables: snapshots at output time                              !
!                                                                      !
!  NOTE: written by Fabian Grosse, based on wrt_his                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_tbnt
      USE mod_ocean
      USE mod_biology
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
      USE omega_mod, ONLY : scale_omega
      USE uv_rotate_mod, ONLY : uv_rotate2d
      USE uv_rotate_mod, ONLY : uv_rotate3d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Fcount, gfactor, gtype, status
      integer :: ifield, itrc, ivar
      
      real(r8) :: scale, dtBIOfac
      
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
#ifdef SOLVE3D
      real(r8), allocatable :: Ur3d(:,:,:)
      real(r8), allocatable :: Vr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
#endif
!
      SourceFile='wrt_tbnt.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out time-integrated TBNT fields when appropriate.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
        gfactor=-1
# else
        gfactor=1
# endif
!
!  Set time and time-record index.
!
      TBNT_NC(ng)%Rindex=TBNT_NC(ng)%Rindex+1
      Fcount=TBNT_NC(ng)%Fcount
      TBNT_NC(ng)%Nrec(Fcount)=TBNT_NC(ng)%Nrec(Fcount)+1
!
!  Write out time.
!
      CALL netcdf_put_fvar (ng, iNLM, TBNT_nc(ng)%name,                 &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/TBNT_NC(ng)%Rindex/), (/1/),              &
     &                      ncid = TBNT_NC(ng)%ncid,                    &
     &                      varid = TBNT_NC(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out grid cell volume (m3)
!
      IF (TBNTout(idTBNTgvol,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, TBNT_NC(ng)%ncid,                  &
     &                     TBNT_NC(ng)%Vid(idTBNTgvol),                 &
     &                     TBNT_NC(ng)%Rindex, gtype,                   &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#endif
	                   TBNT(ng) % TBNTgvol,                         &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTBNTgvol)),                &
     &                        TBNT_NC(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF      
!
!  Write out 3D biological tracer fields.
!  (tracer fields from OCEAN structure are used)
!
!
      DO ivar=1,NTBNTTRC(ng)
        ifield=idTBNTtrc3(ivar)
        itrc=TBNT(ng)%TBNTtbnt2oce(ivar)
        IF (TBNTout(ifield,ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, TBNT_NC(ng)%ncid,                &
     &                       TBNT_NC(ng)%Vid(ifield),                   &
     &                       TBNT_NC(ng)%Rindex, gtype,                 &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       OCEAN(ng) % t(:,:,:,NOUT,itrc),            &
     &                       SetFillVal = .FALSE.)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)),                  &
     &                          TBNT_NC(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out 2D biological flux fields.
!
!      dtBIOfac = REAL(BioIter(ng),r8)/(dt(ng)*sec2day*nTBNT(ng))

      DO ivar=1,nTBNTflx2
        ifield=idTBNTflx2(ivar)
        IF (TBNTout(ifield,ng)) THEN
          scale=1.0_r8!*dtBIOfac
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, TBNT_NC(ng)%ncid,                &
     &                       TBNT_NC(ng)%Vid(ifield),                   &
     &                       TBNT_NC(ng)%Rindex, gtype,                 &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       TBNT(ng) % TBNTflx2(:,:,ivar),             &
     &                       SetFillVal = .FALSE.)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)),                  &
     &                          TBNT_NC(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out 3D biological flux fields.
!
      DO ivar=1,nTBNTflx3
        ifield=idTBNTflx3(ivar)
        IF (TBNTout(ifield,ng)) THEN
          scale=1.0_r8!*dtBIOfac
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, TBNT_NC(ng)%ncid,                &
     &                       TBNT_NC(ng)%Vid(ifield),                   &
     &                       TBNT_NC(ng)%Rindex, gtype,                 &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       TBNT(ng) % TBNTflx3(:,:,:,ivar),           &
     &                       SetFillVal = .FALSE.)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)),                  &
     &                          TBNT_NC(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out 3D sinking flux fields.
!
      DO ivar=1,nTBNTsink
        ifield=idTBNTsink(ivar)
        IF (TBNTout(ifield,ng)) THEN
          scale=1.0_r8!*dtBIOfac
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, TBNT_NC(ng)%ncid,                &
     &                       TBNT_NC(ng)%Vid(ifield),                   &
     &                       TBNT_NC(ng)%Rindex, gtype,                 &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       TBNT(ng) % TBNTsink(:,:,:,ivar),           &
     &                       SetFillVal = .FALSE.)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)),                  &
     &                          TBNT_NC(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out physical exchange flux fields
!
      DO ivar=1,nTBNTphys
        DO itrc=1,NTBNTTRC(ng)
          ifield=idTBNTphys(itrc,ivar)
          IF (TBNTout(ifield,ng)) THEN
            scale=1.0_r8
            gtype=gfactor*r3dvar
            status=nf_fwrite3d(ng, iNLM, TBNT_NC(ng)%ncid,              &
     &                         TBNT_NC(ng)%Vid(ifield),                 &
     &                         TBNT_NC(ng)%Rindex, gtype,               &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
# endif
     &                         TBNT(ng) % TBNTphys(:,:,:,itrc,ivar),    &
     &                         SetFillVal = .FALSE.)
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,ifield)),                &
     &                            TBNT_NC(ng)%Rindex
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END DO
      END DO
!
!  Write out 3D point source flux fields.
!
      DO ivar=1,nTBNTpsrc
        DO itrc=1,NTBNTTRC(ng)
          ifield=idTBNTpsrc(itrc,ivar)
          IF (TBNTout(ifield,ng)) THEN
            scale=1.0_r8
            gtype=gfactor*r3dvar
            status=nf_fwrite3d(ng, iNLM, TBNT_NC(ng)%ncid,              &
     &                         TBNT_NC(ng)%Vid(ifield),                 &
     &                         TBNT_NC(ng)%Rindex, gtype,               &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
# endif
     &                         TBNT(ng) % TBNTpsrc(:,:,:,itrc,ivar),    &
     &                         SetFillVal = .FALSE.)
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,ifield)),                &
     &                            TBNT_NC(ng)%Rindex
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END DO
      END DO
!
!  Synchronize TBNT NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!
      CALL netcdf_sync(ng, iNLM, TBNT_NC(ng)%name, TBNT_NC(ng)%ncid)
# ifdef NESTING
      IF (Master) WRITE (stdout,20) KOUT, NOUT, TBNT_NC(ng)%Rindex, ng
# else
      IF (Master) WRITE (stdout,20) KOUT, NOUT, TBNT_NC(ng)%Rindex
# endif
!
  10  FORMAT (/,' WRT_TBNT    - error while writing variable: ',a,/,11x,  &
     &        'into TBNT NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_TBNT    - wrote TBNT', t39,                         &
     &        'fields (Index=',i1,',',i1,') in record = ',i7.7,t92,i2.2)
#else
      SUBROUTINE wrt_tbnt
#endif
      RETURN
      END SUBROUTINE wrt_tbnt
      
