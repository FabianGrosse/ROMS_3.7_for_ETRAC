/*
** svn $Id: hydrocarbon_def.h 645 2013-01-22 23:21:54Z arango $
** Modified from fennel_def.h
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2013 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines the hydrocarbon model input parameters in                 **
**  output NetCDF files. It is included in routine "def_info.F".      **
**                                                                    **
************************************************************************
*/

!
!  Define hydrocarbon model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HydroCDR'
      Vinfo( 2)='hydrocarbon decay rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Oxy2HC'
      Vinfo( 2)='ratio of oxygen consumption to hydrocarbon decay'
      Vinfo( 3)='millimole_oxygen millimole_hydrocarbon-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='nday_startoil'
      Vinfo( 2)='starting day of oil spill'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      
      Vinfo( 1)='nday_endoil'
      Vinfo( 2)='ending day of oil spill'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
      
      Vinfo( 1)='iloc_oil'
      Vinfo( 2)='longitude (i) index of the model grid where oil spill is located'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN            

      Vinfo( 1)='jloc_oil'
      Vinfo( 2)='latitude (j) index of the model grid where oil spill is located'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN      
      
      Vinfo( 1)='kloc_oil'
      Vinfo( 2)='vertical (k) index of the model grid where oil spill is located'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN 
      
      Vinfo( 1)='density_oil'
      Vinfo( 2)='density of the oil added to the grid cell at each time step'
      Vinfo( 3)='kg m-3'      
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN       
                  