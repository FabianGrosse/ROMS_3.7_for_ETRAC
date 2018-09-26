/*
** svn $Id: oxysink_var.h 645 2013-01-22 23:21:54Z arango $
** Modified from fennel_var.h
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2013 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the simple oxygen sink               **
**  model variables that are used in input and output NetCDF files.   **
**  The metadata information is read from "varinfo.dat".              **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/
              CASE ('idTvar(iOxyg)')
                idTvar(iOxyg)=varid

/*
**  Adjoint sensitivity state biological tracers.
*/

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
              CASE ('idTads(iOxyg)')
                idTads(iOxyg)=varid
#endif

/*
**  Biological tracers open boundary conditions.
*/
              CASE ('idTbry(iwest,iOxyg)')
                idTbry(iwest,iOxyg)=varid
              CASE ('idTbry(ieast,iOxyg)')
                idTbry(ieast,iOxyg)=varid
              CASE ('idTbry(isouth,iOxyg)')
                idTbry(isouth,iOxyg)=varid
              CASE ('idTbry(inorth,iOxyg)')
                idTbry(inorth,iOxyg)=varid

#ifdef TS_PSOURCE

/*
**  Biological tracers point Source/Sinks (river runoff).
*/
              CASE ('idRtrc(iOxyg)')
                idRtrc(iOxyg)=varid
#endif

#ifdef DIAGNOSTICS_BIO

/*
**  Biological tracers term diagnostics.
*/
              CASE ('iDbio2(iO2fx)')
                iDbio2(iO2fx)=varid
#endif
