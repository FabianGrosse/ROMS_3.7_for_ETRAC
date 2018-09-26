#include "cppdefs.h"
#if defined TBNT_OUT && defined BIO_FENNEL && \
    defined BIO_SEDIMENT && defined DENITRIFICATION && \
    defined SOLVE3D
      SUBROUTINE def_tbnt (ng,ldef)
!
!svn $Id: def_tbnt.F 75 2007-06-07 01:01:18Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates TBNT NetCDF file, it defines its dimensions,   !
!  attributes, and variables.                                          !
!                                                                      !
!  NOTE: written by Fabian Grosse, based on def_diags                  !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE def_var_mod, ONLY : def_var
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)
      
      integer, parameter :: Natt = 25
# ifdef TBNT_NFDOUBLE
      integer, parameter :: NF_RTYPE = NF_TYPE
# else
      integer, parameter :: NF_RTYPE = NF_FOUT
# endif

      integer :: i, ifield, itrc, ivar, j, nvd3, nvd4
      integer :: recdim, status
# if defined WRITE_WATER && defined MASKING
      integer :: xy_pdim, xyz_pdim
# endif
      integer :: DimIDs(32), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)

      integer :: def_dim

      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)

      real(r8) :: Aval(6)

      character (len=200) :: Vinfo(Natt)
      character (len=256) :: fname, ncname
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=TBNT_NC(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new TBNT NetCDF file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), TBNT_NC(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define the grid dimensions - NOTE: all variables set to rho points
!-----------------------------------------------------------------------
!
        DimIDs=0
        
        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xi_rho',     &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN
      
        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xi_u',       &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN
      
        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xi_v',       &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xi_psi',     &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN
      
        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'eta_rho',    &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN
      
        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'eta_u',      &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN
      
        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'eta_v',      &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN
      
        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'eta_psi',    &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xy_psi',     &
     &                 IOBOUNDS(ng)%xy_psi, xy_pdim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xy_rho',     &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xy_u',       &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xy_v',       &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (exit_flag.ne.NoError) RETURN
# endif

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xyz_psi',    &
     &                 IOBOUNDS(ng)%xy_psi*N(ng), xyz_pdim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xyz_rho',    &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xyz_u',      &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xyz_v',      &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'xyz_w',      &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 's_rho',      &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 's_w',        &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'tracer',     &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'boundary',   &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname, 'Nstate',     &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, TBNT_NC(ng)%ncid, ncname,               &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                   &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
#if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
#else
        nvd3=3
        nvd4=4
#endif
!
!  Define dimension vectors for staggered tracer type variables.
!
# if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
# else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
# endif
!
!  Define dimension vectors for staggered u-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
# else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
# endif
!
!  Define dimension vectors for staggered v-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
# else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
# endif
!
!  Define dimension vector for staggered w-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
# else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
# endif
!
!  Initialize unlimited time record dimension.
!
        TBNT_NC(ng)%Rindex=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(i))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, iNLM, TBNT_NC(ng)%ncid, ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,a)') 'averaged ', TRIM(Vname(2,idtime))
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, iNLM, TBNT_NC(ng)%ncid,                       &
                       TBNT_NC(ng)%Vid(idtime), NF_TYPE,                 &
     &                 1, (/recdim/) ,Aval, Vinfo, ncname,               &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define grid cell volume.
!
        Vinfo( 1)=Vname(1,idTBNTgvol)
        Vinfo( 2)=Vname(2,idTBNTgvol)
        Vinfo( 3)=Vname(3,idTBNTgvol)
        Vinfo(14)=Vname(4,idTBNTgvol)
        Vinfo(16)=Vname(1,idTBNTgvol)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(r3dvar,r8)
        status=def_var(ng, iNLM, TBNT_NC(ng)%ncid,                       &
     &                 TBNT_NC(ng)%Vid(idTBNTgvol), NF_RTYPE,            &
     &                 nvd4, t3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 3D biochemical tracer fields.
!
        DO ivar=1,nTBNTtrc3
          ifield=idTBNTtrc3(ivar)
          IF (TBNTout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, iNLM, TBNT_NC(ng)%ncid,                   &
     &                     TBNT_NC(ng)%Vid(ifield), NF_RTYPE,            &
     &                     nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define 2D biochemical flux fields
!
        DO ivar=1,nTBNTflx2
          ifield=idTBNTflx2(ivar)
          IF (TBNTout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng, iNLM, TBNT_NC(ng)%ncid,                   &
     &                     TBNT_NC(ng)%Vid(ifield), NF_RTYPE,            &
     &                     nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define 3D biochemical flux fields.
!
        DO ivar=1,nTBNTflx3
          ifield=idTBNTflx3(ivar)
          IF (TBNTout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, iNLM, TBNT_NC(ng)%ncid,                   &
     &                     TBNT_NC(ng)%Vid(ifield), NF_RTYPE,            &
     &                     nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define 3D sinking fields.
!
        DO ivar=1,nTBNTsink
          ifield=idTBNTsink(ivar)
          IF (TBNTout(ifield,ng)) THEN
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, iNLM, TBNT_NC(ng)%ncid,                   &
     &                     TBNT_NC(ng)%Vid(ifield), NF_RTYPE,            &
     &                     nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define physical flux fields (advection/diffusion)
!
        DO itrc=1,NTBNTTRC(ng)
          DO ivar=1,nTBNTphys
            ifield=idTBNTphys(itrc,ivar)
            IF (TBNTout(ifield,ng)) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
              Vinfo(20)='mask_rho'
# endif
              Vinfo(22)='coordinates'
              Aval(5)=REAL(r3dvar,r8)
              status=def_var(ng, iNLM, TBNT_NC(ng)%ncid,                 &
     &                       TBNT_NC(ng)%Vid(ifield), NF_RTYPE,          &
     &                       nvd4, t3dgrd, Aval, Vinfo, ncname)
              IF (exit_flag.ne.NoError) RETURN
            END IF
          END DO
!
!  Define point source flux fields.
!
          DO ivar=1,nTBNTpsrc
            ifield=idTBNTpsrc(itrc,ivar)
            IF (TBNTout(ifield,ng)) THEN
              Vinfo( 1)=Vname(1,ifield)
              Vinfo( 2)=Vname(2,ifield)
              Vinfo( 3)=Vname(3,ifield)
              Vinfo(14)=Vname(4,ifield)
              Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
              Vinfo(20)='mask_rho'
#  endif
              Vinfo(22)='coordinates'
              Aval(5)=REAL(r3dvar,r8)
              status=def_var(ng, iNLM, TBNT_NC(ng)%ncid,                 &
     &                       TBNT_NC(ng)%Vid(ifield), NF_RTYPE,          &
     &                       nvd4, t3dgrd, Aval, Vinfo, ncname)
              IF (exit_flag.ne.NoError) RETURN
            END IF
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef(ng, iNLM, ncname, TBNT_NC(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info(ng, iNLM, TBNT_NC(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing TBNT file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=TBNT_NC(ng)%name
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open TBNT file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, TBNT_NC(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout, 40) TRIM(ncname)
          RETURN
        END IF
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  TBNT variables. Get variable IDs.
!
        DO i=1,n_var
          ! time
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            TBNT_NC(ng)%Vid(idtime)=var_id(i)
          END IF
          ! volume
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTBNTgvol))) THEN
            got_var(idTBNTgvol)=.TRUE.
            TBNT_NC(ng)%Vid(idTBNTgvol)=var_id(i)
          END IF
          ! 3D tracer variables
          DO ivar=1,nTBNTtrc3
            ifield=idTBNTtrc3(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              TBNT_NC(ng)%Vid(ifield)=var_id(i)
            END IF
          END DO
          ! 2D flux variables
          DO ivar=1,nTBNTflx2
            ifield=idTBNTflx2(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              TBNT_NC(ng)%Vid(ifield)=var_id(i)
            END IF
          END DO
          ! 3D flux variables
          DO ivar=1,nTBNTflx3
            ifield=idTBNTflx3(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              TBNT_NC(ng)%Vid(ifield)=var_id(i)
            END IF
          END DO
          ! sinking flux variables
          DO ivar=1,nTBNTsink
            ifield=idTBNTsink(ivar)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
              got_var(ifield)=.TRUE.
              TBNT_NC(ng)%Vid(ifield)=var_id(i)
            END IF
          END DO
          ! 3D physical fluxes variables
          DO itrc=1,NTBNTTRC(ng)
            DO ivar=1,nTBNTphys
              ifield=idTBNTphys(itrc,ivar)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                TBNT_NC(ng)%Vid(ifield)=var_id(i)
              END IF
            END DO
            ! point source input variables
            DO ivar=1,nTBNTpsrc
              ifield=idTBNTpsrc(itrc,ivar)
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,ifield))) THEN
                got_var(ifield)=.TRUE.
                TBNT_NC(ng)%Vid(ifield)=var_id(i)
              END IF
            END DO
          END DO
        END DO
!
!  Check if TBNT variables are available in input NetCDF file.
!
        ! time
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idtime)),         &
                                        TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        ! volume
        IF (.not.got_var(idTBNTgvol).and.TBNTout(idTBNTgvol,ng)) THEN
          IF (Master) WRITE (stdout,50) TRIM(Vname(1,idTBNTgvol)),     &
                                        TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        ! 3D tracer variables
        DO ivar=1,nTBNTtrc3
          ifield=idTBNTtrc3(ivar)
          IF (.not.got_var(ifield).and.TBNTout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,50) TRIM(Vname(1,ifield)),       &
                                          TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        ! 2D flux variables
        DO ivar=1,nTBNTflx2
          ifield=idTBNTflx2(ivar)
          IF (.not.got_var(ifield).and.TBNTout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,50) TRIM(Vname(1,ifield)),       &
                                          TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        ! 3D flux variables
        DO ivar=1,nTBNTflx3
          ifield=idTBNTflx3(ivar)
          IF (.not.got_var(ifield).and.TBNTout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,50) TRIM(Vname(1,ifield)),       &
                                          TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        ! sinking flux variables
        DO ivar=1,nTBNTsink
          ifield=idTBNTsink(ivar)
          IF (.not.got_var(ifield).and.TBNTout(ifield,ng)) THEN
            IF (Master) WRITE (stdout,50) TRIM(Vname(1,ifield)),       &
                                          TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        ! 3D physical flux variables
        DO itrc=1,NTBNTTRC(ng)
          DO ivar=1,nTBNTphys
            ifield=idTBNTphys(itrc,ivar)
            IF (.not.got_var(ifield).and.TBNTout(ifield,ng)) THEN
              IF (Master) WRITE (stdout,50) TRIM(Vname(1,ifield)),     &
                                            TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END DO
          ! point source flux variables
          DO ivar=1,nTBNTpsrc
            ifield=idTBNTpsrc(itrc,ivar)
            IF (.not.got_var(ifield).and.TBNTout(ifield,ng)) THEN
              IF (Master) WRITE (stdout,50) TRIM(Vname(1,ifield)),     &
                                            TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END DO
        END DO
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (ndefTBNT(ng).gt.0) THEN
          IF ((ntstart(ng)-1).lt.ndefTBNT(ng)) THEN
            TBNT_NC(ng)%Rindex=1+((ntstart(ng)-1)- ndefTBNT(ng)*       &
     &                            ((ntstart(ng)-1)/ndefTBNT(ng)))/     &
     &                           nTBNT(ng)
          ELSE
            TBNT_NC(ng)%Rindex=((ntstart(ng)-1)- ndefTBNT(ng)*         &
     &                          ((ntstart(ng)-1)/ndefTBNT(ng)))/       &
     &                         nTBNT(ng)
          END IF
        ELSE
          TBNT_NC(ng)%Rindex=(ntstart(ng)-1)/nTBNT(ng)
        END IF
        TBNT_NC(ng)%Rindex=MIN(TBNT_NC(ng)%Rindex,rec_size)
      END IF QUERY
!
  10  FORMAT (6x,'DEF_TBNT    - creating TBNT file', t43,              &
     &        ' file, Grid ',i2.2,': ', a)
  20  FORMAT (6x,'DEF_TBNT    - inquiring TBNT file', t43,             &
     &        ' file, Grid ',i2.2,': ', a)
  30  FORMAT (/,' DEF_TBNT    - unable to create TBNT NetCDF file: ',a)
  40  FORMAT (/,' DEF_TBNT    - unable to open TBNT NetCDF file: ',a)
  50  FORMAT (/,' DEF_TBNT    - unable to find variable: ',a,2x,       &
     &        ' in TBNT NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_tbnt
#else
      SUBROUTINE def_tbnt
      RETURN
      END SUBROUTINE def_tbnt
#endif
