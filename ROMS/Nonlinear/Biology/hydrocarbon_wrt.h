/*
** svn $Id: hydrocarbon_wrt.h 645 2013-01-22 23:21:54Z arango
** Modified from fennel_wrt.h
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2013 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Writes hydrocarbon model input parameters into                    **
**  output NetCDF files. It is included in routine "wrt_info.F".      **
**                                                                    **
************************************************************************
*/

!
!  Write out simple oxygen model parameters.
!
      CALL netcdf_put_ivar (ng, model, ncname, 'BioIter',               &
     &                      BioIter(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'HydroCDR',              &
     &                      HydroCDR(ng), (/0/), (/0/),                 &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'Oxy2HC',                &
     &                      Oxy2HC(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'nday_startoil',                &
     &                      nday_startoil(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      
      CALL netcdf_put_fvar (ng, model, ncname, 'nday_endoil',                &
     &                      nday_endoil(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      
      CALL netcdf_put_fvar (ng, model, ncname, 'iloc_oil',                &
     &                      iloc_oil(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN
      
      CALL netcdf_put_fvar (ng, model, ncname, 'jloc_oil',                &
     &                      jloc_oil(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, model, ncname, 'kloc_oil',                &
     &                      kloc_oil(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN                              

      CALL netcdf_put_fvar (ng, model, ncname, 'density_oil',                &
     &                      density_oil(ng), (/0/), (/0/),                   &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN    