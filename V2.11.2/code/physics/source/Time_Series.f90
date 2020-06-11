!
! Copyright 2013 RBINS-MUMM
!
! Licensed under the EUPL, Version 1.1 or - as soon they will be approved by
! the European Commission - subsequent versions of the EUPL (the "Licence");
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!
! http://ec.europa.eu/idabc/eupl
!
! Unless required by the applicable law or agreed to in writing, software
! distributed under the Licence is distributed on an "AS IS" basis,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.

SUBROUTINE time_series
!************************************************************************
!
! *time_series* Time series output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Time_Series.f90  V2.11
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
!
! External calls - read_cif_params, usrdef_tsr0d_vals, usrdef_tsr2d_vals,
!                  usrdef_tsr3d_vals
!
! Internal calls - time_series_grid, time_series_init 
!
! Module calls - close_filepars, combine_write_stats_loc, combine_write_submod,
!                define_out0d_vals, define_out2d_vals, define_out3d_vals,
!                error_alloc, loop_index, mask_array, write_vars
!
!************************************************************************
!
USE datatypes
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE switches
USE syspars
USE timepars
USE check_model, ONLY: check_out_filepars, check_out_gpars, check_statlocs, &
                     & check_variables
USE datatypes_init, ONLY: filepars_init, outgpars_init, statlocs_init, &
                        & varatts_init
USE default_model, ONLY: default_out_files, default_out_gpars
USE error_routines, ONLY: error_abort, error_alloc, error_alloc_struc, &
                        & error_diff_vals_arrlist, error_limits_arr, &
                        & warning_reset_arr_struc
USE grid_routines, ONLY: global_mask, local_proc, mask_array, Zcoord_var
USE inout_paral, ONLY: combine_write_stats_loc, combine_write_submod
USE inout_routines, ONLY: close_filepars, inout_atts_out, write_time, &
                        & write_vars
USE model_output, ONLY: define_out0d_vals, define_out2d_vals, define_out3d_vals
USE paral_comms, ONLY: combine_mod
USE paral_utilities, ONLY: sum_vars
USE reset_model, ONLY: reset_out_files, reset_out_gpars, reset_out_stats, &
                     & reset_out_vars
USE time_routines, ONLY: date_to_year, diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index

IMPLICIT none

!
!* Local variables
!
LOGICAL, SAVE :: last_call
LOGICAL :: flag, fmask, gridded, packing, time_grid
CHARACTER (LEN=lentime) :: refdate
INTEGER, SAVE :: maxstats, novars0d, novars2d, novars3d
INTEGER :: i, ii, iivar, imax, iproc, iset, istat, ivar, i0, i1, i2, i3, j, jj, &
         & jmax, j1, j2, j3, k, kk, l, lloc, maxstatsloc, maxvars, ncout, &
         & nocoords0d, nocoords2d, nocoords3d, nodat, nodim, nostats, &
         & noutvars0d, noutvars2d, noutvars3d, nowetlocout, nowetout, npcc, &
         & nrank, nrout, nzout, time_format, vcoord, zvarid
REAL :: fill_value
REAL (KIND=kndrlong) :: rtime
TYPE (FileParams) :: filepars

INTEGER, DIMENSION(3) :: tlims, xlims, xlimsglb, ylims, ylimsglb, zlims
INTEGER, DIMENSION(novarstsr) :: vecids
INTEGER, DIMENSION(2,2,nprocs) :: lprocs
REAL, DIMENSION(novarstsr) :: outdat

LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: maskvals
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: lmask, nclocout, nostatsloc, &
                                          & nrlocout, ntout 
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: istatpos, ivarstsr0d, &
                                            & ivarstsr2d, ivarstsr3d, &
                                            & jstatpos, nostatsprocs
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: limloc, lstatsprocs
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: limprocs
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb
REAL, ALLOCATABLE, DIMENSION(:) :: gsigout, out1dsub 
REAL, ALLOCATABLE, DIMENSION(:,:) :: depout, out2dsub, xout, yout, zetout
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: gsout, out3dsub, zout
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: out4dsub
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: outvars

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*last_call*   LOGICAL .TRUE. if the routine has been called for the last time
!*istatpos*    INTEGER X-index of local output station
!*ivarstsr0d*  INTEGER Array of indices (per file) of 0-D output variables
!*ivarstsr2d*  INTEGER Array of indices (per file) of 2-D output variables
!*ivarstsr3d*  INTEGER Array of indices (per file) of 3-D output variables
!*jstatpos*    INTEGER Y-index of local output station
!*limloc*      INTEGER Start/End/Increment values of local (regular) output
!                      grid in X- and Y-directions
!*limprocs*    INTEGER Start/End indices of local array section within global
!                      array (regular grid)
!*lstatprocs*  INTEGER Output labels of local stations for each process domain
!*nclocout*    INTEGER X-dimension of local (regular) output grid
!*nostatsloc*  INTEGER Number of local output stations (per file)
!*nostatprocs* INTEGER Number of local output stations for each process domain 
!*novars0d*    INTEGER Total number of 0-D output variables
!*novars2d*    INTEGER Total number of 2-D output variables
!*novars3d*    INTEGER Total number of 3-D output variables
!*nroutloc*    INTEGER Y-dimension of local (regular) output grid
!*ntout*       INTEGER Output time counter
!*depout*      REAL    Output depth array                                   [m]
!*xout*        REAL    Output X-coordinates                  [degree_east or m]
!*yout*        REAL    Output Y-coordinates                 [degree_north or m]
!*zetout*      REAL    Output surface elevations                            [m]
!
!------------------------------------------------------------------------------
!

IF (nt.EQ.0) last_call=.FALSE.
IF (last_call) RETURN

procname(pglev+1) = 'time_series'
CALL log_timer_in(npcc)

!
!1. Initialisation on first call
!-------------------------------
!    

IF (nt.EQ.0) CALL time_series_init

!
!2. Write output grid
!--------------------
!

CALL time_series_grid

!
!3. Write time series output
!---------------------------
!

iset_300:DO iset=1,nosetstsr

   tlims = tsrgpars(iset)%tlims
   IF (loop_index(tlims,nt)) THEN

!
!3.1 Initialise
!--------------
!

      gridded = tsrgpars(iset)%gridded
      nodim = tsrgpars(iset)%nodim
      ncout = tsrgpars(iset)%ncout
      nrout = tsrgpars(iset)%nrout
      nzout = tsrgpars(iset)%nzout
      nowetout = tsrgpars(iset)%nowetout
      nostats = tsrgpars(iset)%nostats
      xlims = limloc(:,1,iset)
      ylims = limloc(:,2,iset)
      zlims = tsrgpars(iset)%zlims
      ntout(iset) = ntout(iset) + 1
      lprocs = limprocs(:,:,:,iset)
      fill_value = double_fill

      noutvars0d = tsr0d(iset)%novars
      noutvars2d = tsr2d(iset)%novars
      noutvars3d = tsr3d(iset)%novars

      nocoords0d = tsr0d(iset)%nocoords
      nocoords2d = tsr2d(iset)%nocoords
      nocoords3d = tsr3d(iset)%nocoords

!
!3.2 Create mask array
!---------------------
!

      IF (gridded) THEN
         ALLOCATE (maskvals(nclocout(iset),nrlocout(iset)),STAT=errstat)
         CALL error_alloc('maskvals',2,(/nclocout(iset),nrlocout(iset)/),&
                        & kndlog)
         IF (SIZE(maskvals).GT.0) THEN
            fmask = .TRUE.
            maskvals = mask_array(nclocout(iset),nrlocout(iset),&
                                & limloc(:,:,iset),'C')
         ELSE
            fmask = .FALSE.
         ENDIF
      ENDIF

!
!3.3 Allocate arrays for temporary data storage
!----------------------------------------------
!

!     ---0-D
      ALLOCATE (out1dsub(noutvars0d),STAT=errstat)
      CALL error_alloc('out1dsub',1,(/noutvars0d/),kndrtype)
      IF (noutvars0d.GT.0) out1dsub = 0.0

!     ---regular output
      IF (gridded) THEN
         ALLOCATE (out3dsub(nclocout(iset),nrlocout(iset),noutvars2d),&
                 & STAT=errstat)
         CALL error_alloc('out3dsub',3,(/nclocout(iset),nrlocout(iset),&
                                       & noutvars2d/),kndrtype)
         IF (noutvars2d.GT.0) out3dsub = 0.0
         ALLOCATE (out4dsub(nclocout(iset),nrlocout(iset),nzout,&
                          & noutvars3d),STAT=errstat)
         CALL error_alloc('out4dsub',4,(/nclocout(iset),nrlocout(iset),&
                                       & nzout,noutvars3d/),kndrtype)
         IF (noutvars3d.GT.0) out4dsub = 0.0

!     ---station output
      ELSE
         ALLOCATE (out2dsub(nostatsloc(iset),noutvars2d),STAT=errstat)
         CALL error_alloc('out2dsub',2,(/nostatsloc(iset),noutvars2d/),&
                         & kndrtype)
         IF (noutvars2d.GT.0) out2dsub = 0.0
         ALLOCATE (out3dsub(nostatsloc(iset),nzout,noutvars3d),STAT=errstat)
         CALL error_alloc('out3dsub',3,&
                       & (/nostatsloc(iset),nzout,noutvars3d/),kndrtype)
         IF (noutvars3d.GT.0) out3dsub = 0.0
      ENDIF

!
!3.4 Store data into subarrays
!-----------------------------
!
!3.4.1 0-D data
!--------------
!

      IF (tsr0d(iset)%defined) THEN
         outvars(1:noutvars0d) = tsrvars(ivarstsr0d(iset,1:noutvars0d)) 
         CALL define_out0d_vals(out1dsub,noutvars0d,&
                              & outvars=outvars(1:noutvars0d))
         IF (ANY(outvars(1:noutvars0d)%ivarid.EQ.0)) THEN
            CALL usrdef_tsr0d_vals(outdat(1:novars0d),novars0d)
            WHERE (outvars(1:noutvars0d)%ivarid.EQ.0)
               out1dsub = outdat(ivarstsr0d(iset,1:noutvars0d))
            END WHERE
         ENDIF
      ENDIF

!
!3.4.2 2-D data
!--------------
!

      IF (nodim.GT.0.AND.tsr2d(iset)%defined) THEN
         outvars(1:noutvars2d) = &
                              & tsrvars(novars0d+ivarstsr2d(iset,1:noutvars2d))
         flag = ANY(outvars(1:noutvars2d)%ivarid.EQ.0)

!        ---regular data
         IF (gridded) THEN
            i_3421:DO i=xlims(1),xlims(2),xlims(3)
            j_3421:DO j=ylims(1),ylims(2),ylims(3)
               ii = (i-xlims(1))/xlims(3) + 1
               jj = (j-ylims(1))/ylims(3) + 1
               CALL define_out2d_vals(out3dsub(ii,jj,:),i,j,noutvars2d,&
                                    & outvars=outvars(1:noutvars2d))
               IF (flag) THEN
                  CALL usrdef_tsr2d_vals(outdat(1:novars2d),i,j,novars2d)
                  WHERE (outvars(1:noutvars2d)%ivarid.EQ.0)
                     out3dsub(ii,jj,:)  = outdat(ivarstsr2d(iset,1:noutvars2d))
                  END WHERE
               ENDIF
            ENDDO j_3421
            ENDDO i_3421

!        ---irregular data
         ELSE
            l_3422:DO l=1,nostatsloc(iset)
               i = istatpos(l,iset); j = jstatpos(l,iset)
               CALL define_out2d_vals(out2dsub(l,:),i,j,noutvars2d,&
                                    & outvars=outvars(1:noutvars2d))
               IF (flag) THEN
                  CALL usrdef_tsr2d_vals(outdat(1:novars2d),i,j,novars2d)
                  WHERE (outvars(1:noutvars2d)%ivarid.EQ.0)
                     out2dsub(l,:) = outdat(ivarstsr2d(iset,1:noutvars2d))
                  END WHERE
               ENDIF
            ENDDO l_3422
         ENDIF

      ENDIF

!
!3.4.3 3-D data
!--------------
!

      IF (nodim.EQ.3.AND.tsr3d(iset)%defined) THEN
         outvars(1:noutvars3d) = &
                     & tsrvars(novars0d+novars2d+ivarstsr3d(iset,1:noutvars3d))
         flag = ANY(outvars(1:noutvars3d)%ivarid.EQ.0)

!        ---regular data
         IF (gridded) THEN
            i_3431:DO i=xlims(1),xlims(2),xlims(3)
            j_3431:DO j=ylims(1),ylims(2),ylims(3)
            k_3431:DO k=zlims(1),zlims(2),zlims(3)
               ii = (i-xlims(1))/xlims(3) + 1
               jj = (j-ylims(1))/ylims(3) + 1
               kk = (k-zlims(1))/zlims(3) + 1
               CALL define_out3d_vals(out4dsub(ii,jj,kk,:),i,j,k,noutvars3d,&
                                    & outvars=outvars(1:noutvars3d))
               IF (flag) THEN
                  CALL usrdef_tsr3d_vals(outdat(1:novars3d),i,j,k,novars3d)
                  WHERE (outvars(1:noutvars3d)%ivarid.EQ.0)
                     out4dsub(ii,jj,kk,:) = &
                                         & outdat(ivarstsr3d(iset,1:noutvars3d))
                  END WHERE
               ENDIF
            ENDDO k_3431
            ENDDO j_3431
            ENDDO i_3431

!        ---irregular data
         ELSE
            l_3432:DO l=1,nostatsloc(iset)
               i = istatpos(l,iset); j = jstatpos(l,iset)
               k_34231:DO k=zlims(1),zlims(2),zlims(3)
                  kk = (k-zlims(1))/zlims(3) + 1
                  CALL define_out3d_vals(out3dsub(l,kk,:),i,j,k,noutvars3d,&
                                       & outvars=outvars(1:noutvars3d))
                  IF (flag) THEN
                     CALL usrdef_tsr3d_vals(outdat(1:novars3d),i,j,k,novars3d)
                     WHERE (outvars(1:noutvars3d)%ivarid.EQ.0)
                        out3dsub(l,kk,:) = outdat(ivarstsr3d(iset,1:noutvars3d))
                     END WHERE
                  ENDIF
               ENDDO k_34231
            ENDDO l_3432
         ENDIF

      ENDIF

!         
!3.4.4 Insert fill values
!------------------------
!

      IF (gridded.AND.fmask) THEN
         IF (tsr2d(iset)%defined.AND.tsr2d(iset)%fill) THEN
            ivar_3441: DO ivar=1,noutvars2d
               WHERE (.NOT.maskvals)
                  out3dsub(:,:,ivar) = fill_value
               END WHERE
            ENDDO ivar_3441
         ENDIF
         IF (tsr3d(iset)%defined.AND.tsr3d(iset)%fill) THEN
            ivar_3442: DO ivar=1,noutvars3d
            k_3442: DO k=1,nzout
               WHERE (.NOT.maskvals)
                  out4dsub(:,:,k,ivar) = fill_value
               END WHERE
            ENDDO k_3442
            ENDDO ivar_3442
         ENDIF
      ENDIF

!
!3.5 Write data
!--------------
!
!3.5.1 0-D output
!----------------
!

      IF (master.AND.tsr0d(iset)%defined) THEN
         outvars(1:noutvars0d) = tsrvars(ivarstsr0d(iset,1:noutvars0d)) 
         vecids(1:noutvars0d) = (/(nocoords0d+ivar,ivar=1,noutvars0d)/)
         CALL write_vars(out1dsub,tsr0d(iset),0,varatts=outvars(1:noutvars0d),&
                       & vecids=vecids(1:noutvars0d))
      ENDIF

!
!3.5.2 2-D output
!----------------
!

      IF (tsr2d(iset)%defined) THEN
         outvars(1:noutvars2d) = &
                              & tsrvars(novars0d+ivarstsr2d(iset,1:noutvars2d))
         vecids(1:noutvars2d) = (/(nocoords2d+ivar,ivar=1,noutvars2d)/)
         IF (gridded) THEN
            CALL combine_write_submod(out3dsub,tsr2d(iset),0,&
                                   & (/ncout,nrout,noutvars2d/),lprocs,&
                                   & nowet=nowetout,&
                                   & varatts=outvars(1:noutvars2d),&
                                   & vecids=vecids(1:noutvars2d))
         ELSE
            CALL combine_write_stats_loc(out2dsub,tsr2d(iset),0,maxstats,&
                                   & nostats,nostatsprocs(:,iset),&
                                   & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                   & varatts=outvars(1:noutvars2d),&
                                   & vecids=vecids(1:noutvars2d))
         ENDIF
      ENDIF

!
!3.5.3 3-D output
!----------------
!

      IF (tsr3d(iset)%defined) THEN
         outvars(1:noutvars3d) = &
                     & tsrvars(novars0d+novars2d+ivarstsr3d(iset,1:noutvars3d))
         vecids(1:noutvars3d) = (/(nocoords3d+ivar,ivar=1,noutvars3d)/)
         IF (gridded) THEN
            CALL combine_write_submod(out4dsub,tsr3d(iset),0,&
                                  & (/ncout,nrout,nzout,noutvars3d/),lprocs,&
                                  & nowet=nowetout,&
                                  & varatts=outvars(1:noutvars3d),&
                                  & vecids=vecids(1:noutvars3d))
         ELSE
            CALL combine_write_stats_loc(out3dsub,tsr3d(iset),0,maxstats,&
                                   & nostats,nostatsprocs(:,iset),&
                                   & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                   & varatts=outvars(1:noutvars3d),&
                                   & vecids=vecids(1:noutvars3d))
         ENDIF
      ENDIF

!
!3.6 Deallocate arrays
!---------------------
!

      IF (ALLOCATED(maskvals)) DEALLOCATE (maskvals)
      DEALLOCATE(out1dsub)
      IF (gridded) THEN
         DEALLOCATE (out3dsub,out4dsub)
      ELSE
         DEALLOCATE (out2dsub,out3dsub)
      ENDIF
   ENDIF

!
!3.7 Close files after last output
!---------------------------------
!

   IF (master) THEN
      IF (nt.EQ.tlims(2).OR.cold_start) THEN
         IF (tsr0d(iset)%defined) CALL close_filepars(tsr0d(iset))
         IF (tsr2d(iset)%defined) CALL close_filepars(tsr2d(iset))
         IF (tsr3d(iset)%defined) CALL close_filepars(tsr3d(iset))
      ENDIF
   ENDIF
 
ENDDO iset_300

!
!4. Deallocate arrays at last output time
!----------------------------------------
!
!---check if last output has been written
last_call = .TRUE.
IF (.NOT.cold_start) THEN
   iset_410: DO iset=1,nosetstsr
      IF (.NOT.(nt.GE.tsrgpars(iset)%tlims(2))) last_call = .FALSE.
   ENDDO iset_410
ENDIF

!---deallocate
IF (last_call) THEN
   DEALLOCATE (ivarstsr0d,ivarstsr2d,ivarstsr3d,outvars)
   DEALLOCATE (nclocout,nrlocout,limprocs,limloc)
   DEALLOCATE (nostatsprocs,nostatsloc,lstatsprocs,istatpos,jstatpos)
   DEALLOCATE (ntout)
ENDIF

CALL log_timer_out(npcc,itm_user)


RETURN

CONTAINS

!========================================================================

SUBROUTINE time_series_init
!************************************************************************
!
! *time_series_init* Initialise parameters for time series output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Time_Series.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - time_series
!
! External calls - read_cif_params, usrdef_tsr_params
!
! Module calls - check_out_filepars, check_out_gpars, check_statlocs,
!   check_variables, default_out_files, default_out_gpars, error_abort,
!   error_alloc, error_alloc_struc, error_diff_vals_arrlist, error_limits_arr,
!   filepars_init, inout_atts_out, lim_dims, local_proc, outgpars_init,
!   reset_out_files, reset_out_gpars, reset_out_stats, reset_out_vars,
!   statlocs_init, sum_vars, varatts_init, warning_reset_arr_struc
!
!************************************************************************
!


procname(pglev+1) = 'time_series_init'
CALL log_timer_in()

!
!1. Parameters for time series output
!------------------------------------
!
!1.1 Allocate and initialise
!---------------------------
!
!---variables
ALLOCATE (tsrvars(novarstsr),STAT=errstat)
CALL error_alloc_struc('tsrvars',1,(/novarstsr/),'VariableAtts')
CALL varatts_init(tsrvars)

!---file parameters
ALLOCATE (tsr0d(nosetstsr),STAT=errstat)
CALL error_alloc_struc('tsr0d',1,(/nosetstsr/),'FileParams')
CALL filepars_init(tsr0d)
ALLOCATE (tsr2d(nosetstsr),STAT=errstat)
CALL error_alloc_struc('tsr2d',1,(/nosetstsr/),'FileParams')
CALL filepars_init(tsr2d)
ALLOCATE (tsr3d(nosetstsr),STAT=errstat)
CALL error_alloc_struc('tsr3d',1,(/nosetstsr/),'FileParams')
CALL filepars_init(tsr3d)

!---output specifiers
ALLOCATE (tsrgpars(nosetstsr),STAT=errstat)
CALL error_alloc_struc('tsrgpars',1,(/nosetstsr/),'OutFileParams')
CALL outgpars_init(tsrgpars)

!---output stations
ALLOCATE (lstatstsr(nosetstsr,nostatstsr),STAT=errstat)
CALL error_alloc('lstatstsr',2,(/nosetstsr,nostatstsr/),kndint)
IF (SIZE(lstatstsr).GT.0) lstatstsr = 0
ALLOCATE (tsrstatlocs(nostatstsr),STAT=errstat)
CALL error_alloc_struc('tsrstatsloc',1,(/nostatstsr/),'StationLocs')
CALL statlocs_init(tsrstatlocs)
 
!---variable indices
ALLOCATE (ivarstsr(nosetstsr,novarstsr),STAT=errstat)
CALL error_alloc('ivarstsr',2,(/nosetstsr,novarstsr/),kndint)
ivarstsr = 0

!
!1.2 Defaults
!------------
!
!---file parameters
CALL default_out_files(tsr0d,'TS')
CALL default_out_files(tsr2d,'TS')
CALL default_out_files(tsr3d,'TS')

!---output parameters
CALL default_out_gpars(tsrgpars)

!
!1.3 Define
!----------
!

IF (ciffile%status.EQ.'R') THEN
   CALL read_cif_params(icif_tsout)
ELSE
   CALL usrdef_tsr_params
ENDIF

!
!1.4 Reset
!---------
!
!---file parameters
CALL reset_out_files(tsr0d,'tsr0d',tsrgpars,.TRUE.)
CALL reset_out_files(tsr2d,'tsr2d',tsrgpars,.TRUE.)
CALL reset_out_files(tsr3d,'tsr3d',tsrgpars,.TRUE.)

!---variable atrributes
CALL reset_out_vars(tsrvars)

!---output grid dimension
iset_141: DO iset=1,nosetstsr
   nodim = 0
   ivar_1411: DO ivar=1,novarstsr
      iivar = ivarstsr(iset,ivar)
      IF (iivar.GT.0) THEN
         IF (tsrvars(iivar)%nrank.EQ.3) THEN
            nodim = 3; EXIT ivar_1411
         ENDIF
      ENDIF
   ENDDO ivar_1411
   IF (nodim.EQ.0) THEN
      ivar_1412: DO ivar=1,novarstsr
         iivar = ivarstsr(iset,ivar)
         IF (iivar.GT.0) THEN
            IF (tsrvars(iivar)%nrank.EQ.2) THEN
               nodim = 2; EXIT ivar_1412
            ENDIF
         ENDIF
      ENDDO ivar_1412
   ENDIF
   tsrgpars(iset)%nodim = nodim
   IF (iopt_grid_nodim.EQ.2.AND.tsrgpars(iset)%nodim.EQ.3) THEN
      CALL warning_reset_arr_struc(tsrgpars(iset)%nodim,'tsrgpars','%nodim',2,&
                                 & 1,(/iset/))
   ENDIF
ENDDO iset_141

!---output grid parameters
CALL reset_out_gpars(tsrgpars,'tsrgpars',.TRUE.,.TRUE.)

!---station attributes
CALL reset_out_stats(tsrstatlocs)

!---file parameters
tsr0d%nodim = 0; tsr2d%nodim = 2; tsr3d%nodim = 3
iset_142: DO iset=1,nosetstsr
   IF (tsrgpars(iset)%nodim.LT.3) THEN
      CALL warning_reset_arr_struc(tsr3d(iset)%defined,'tsr3d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
   IF (tsrgpars(iset)%nodim.EQ.0) THEN
      CALL warning_reset_arr_struc(tsr2d(iset)%defined,'tsr2d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
ENDDO iset_142

!
!1.5 Check
!---------
!
!---output variables
CALL check_variables(tsrvars,'tsrvars')

!---file parameters
CALL check_out_filepars(tsr0d,'tsr0d')
CALL check_out_filepars(tsr2d,'tsr2d')
CALL check_out_filepars(tsr3d,'tsr3d')

!---output parameters
CALL check_out_gpars(tsrgpars,'tsrgpars',nostatstsr)

!---output stations
iset_151: DO iset=1,nosetstsr
istat_151: DO istat=1,tsrgpars(iset)%nostats
   CALL error_limits_arr(lstatstsr(iset,istat),'lstatstsr',1,&
                       & nostatstsr,2,(/iset,istat/))
ENDDO istat_151
ENDDO iset_151
CALL check_statlocs(tsrstatlocs,'tsrstatlocs')

!---variable indices
iset_152: DO iset=1,nosetstsr
   ivar_1521: DO ivar=1,novarstsr
      CALL error_limits_arr(ivarstsr(iset,ivar),'ivarstsr',0,novarstsr,2,&
                         & (/iset,ivar/))
   ENDDO ivar_1521
   CALL error_diff_vals_arrlist(ivarstsr(iset,:),'ivarstsr',2,2,(/iset,1/),&
                             & .TRUE.)
ENDDO iset_152

!
!2. Variable indices
!-------------------
!
!2.1 Number of variables per dimension
!-------------------------------------
!

novars0d = COUNT(tsrvars%nrank.EQ.0)
novars2d = COUNT(tsrvars%nrank.EQ.2)
novars3d = COUNT(tsrvars%nrank.EQ.3)

!
!2.2 Allocate and initialise
!---------------------------
!
 
ALLOCATE (ivarstsr0d(nosetstsr,novars0d),STAT=errstat)
CALL error_alloc('ivarstsr0d',2,(/nosetstsr,novars0d/),kndint)
IF (novars0d.GT.0) ivarstsr0d = 0 
ALLOCATE (ivarstsr2d(nosetstsr,novars2d),STAT=errstat)
CALL error_alloc('ivarstsr2d',2,(/nosetstsr,novars2d/),kndint)
IF (novars2d.GT.0) ivarstsr2d = 0 
ALLOCATE (ivarstsr3d(nosetstsr,novars3d),STAT=errstat)
CALL error_alloc('ivarstsr3d',2,(/nosetstsr,novars3d/),kndint)
IF (novars3d.GT.0) ivarstsr3d = 0 
maxvars = MAX(novars0d,novars2d,novars3d)
ALLOCATE (outvars(maxvars),STAT=errstat)
CALL error_alloc_struc('outvars',1,(/maxvars/),'VariableAtts')

!
!2.3 Define
!----------
!

iset_230: DO iset=1,nosetstsr
   i0 = 0; i2 = 0; i3 = 0
   ivar_231: DO ivar=1,novarstsr
      l = ivarstsr(iset,ivar)
      IF (l.GT.0) THEN
         nrank = tsrvars(l)%nrank
         SELECT CASE (nrank)
         CASE (0)
            i0 = i0 + 1
            ivarstsr0d(iset,i0) = l
         CASE (2)
            i2 = i2 + 1
            ivarstsr2d(iset,i2) = l
         CASE (3)
            i3 = i3 + 1
            ivarstsr3d(iset,i3) = l
         END SELECT
      ENDIF
   ENDDO ivar_231
   tsr0d(iset)%novars = i0; tsr2d(iset)%novars = i2; tsr3d(iset)%novars = i3
   ivarstsr2d(iset,:) = ivarstsr2d(iset,:) - novars0d
   ivarstsr3d(iset,:) = ivarstsr3d(iset,:) - novars0d - novars2d
ENDDO iset_230

!
!2.4 Write output only if number of variables is greater than zero
!-----------------------------------------------------------------
!

iset_240: DO iset=1,nosetstsr
   IF (tsr0d(iset)%novars.EQ.0) THEN
      CALL warning_reset_arr_struc(tsr0d(iset)%defined,'tsr0d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
   IF (tsr2d(iset)%novars.EQ.0) THEN
      CALL warning_reset_arr_struc(tsr2d(iset)%defined,'tsr2d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
   IF (tsr3d(iset)%novars.EQ.0) THEN
      CALL warning_reset_arr_struc(tsr3d(iset)%defined,'tsr3d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
ENDDO iset_240

!
!3. Space output specifiers
!--------------------------
!
!3.1 Allocate
!------------
!
!---regular output
ALLOCATE (nclocout(nosetstsr),STAT=errstat)
CALL error_alloc('nclocout',1,(/nosetstsr/),kndint)
nclocout = 0
ALLOCATE (nrlocout(nosetstsr),STAT=errstat)
CALL error_alloc('nrlocout',1,(/nosetstsr/),kndint)
nrlocout = 0
ALLOCATE (limprocs(2,2,nprocs,nosetstsr),STAT=errstat)
CALL error_alloc('limprocs',4,(/2,2,nprocs,nosetstsr/),kndint)
limprocs = 0
ALLOCATE (limloc(3,2,nosetstsr),STAT=errstat)
CALL error_alloc('limloc',3,(/3,2,nosetstsr/),kndint)
limloc = 0

!---station output
ALLOCATE (nostatsprocs(nprocs,nosetstsr),STAT=errstat)
CALL error_alloc('nostatsprocs',2,(/nprocs,nosetstsr/),kndint)
nostatsprocs = 0
ALLOCATE (nostatsloc(nosetstsr),STAT=errstat)
CALL error_alloc('nostatsloc',1,(/nosetstsr/),kndint)
nostatsloc = 0

!
!3.2 Regular output
!------------------
!

iset_320: DO iset=1,nosetstsr

   IF (tsrgpars(iset)%nodim.GT.0.AND.tsrgpars(iset)%gridded) THEN

      tsrgpars(iset)%ncout = lim_dims(tsrgpars(iset)%xlims)
      tsrgpars(iset)%nrout = lim_dims(tsrgpars(iset)%ylims)
      tsrgpars(iset)%nzout = lim_dims(tsrgpars(iset)%zlims)

      iproc_321: DO iproc=1,nprocs

!        ---X-direction
         xlims = tsrgpars(iset)%xlims
         nodat = 0; limprocs(1,2,iproc,iset) = 0; ii = 0
         IF (idloc.EQ.idprocs(iproc)) nclocout(iset) = 0
         i_3211: DO i=xlims(1),xlims(2),xlims(3)
            ii = ii + 1
            IF (i.GE.nc1procs(iproc).AND.i.LE.nc2procs(iproc)) THEN
               limprocs(1,2,iproc,iset) = ii
               nodat = nodat + 1
               IF (idloc.EQ.idprocs(iproc)) THEN
                  nclocout(iset) = nclocout(iset) + 1
                  i2 = i
               ENDIF
            ENDIF
         ENDDO i_3211
         limprocs(1,1,iproc,iset) = limprocs(1,2,iproc,iset) - nodat + 1
         IF (idloc.EQ.idprocs(iproc)) THEN
            IF (nclocout(iset).GT.0) THEN
               i2 = i2 - nc1loc + 1
               i1 = i2 + (1-nclocout(iset))*xlims(3)
               i3 = xlims(3)
            ELSE
               i1 = 1; i2 = 0; i3 = 1
            ENDIF
            limloc(:,1,iset) = (/i1,i2,i3/)
         ENDIF
            
!        ---Y-direction
         ylims = tsrgpars(iset)%ylims
         nodat = 0; limprocs(2,2,iproc,iset) = 0; jj = 0
         IF (idloc.EQ.idprocs(iproc)) nrlocout(iset) = 0
         j_3212: DO j=ylims(1),ylims(2),ylims(3)
            jj = jj + 1
            IF (j.GE.nr1procs(iproc).AND.j.LE.nr2procs(iproc)) THEN
               limprocs(2,2,iproc,iset) = jj
               nodat = nodat + 1
               IF (idloc.EQ.idprocs(iproc)) THEN
                  nrlocout(iset) = nrlocout(iset) + 1
                  j2 = j
               ENDIF
            ENDIF
         ENDDO j_3212
         limprocs(2,1,iproc,iset) = limprocs(2,2,iproc,iset) - nodat + 1
         IF (idloc.EQ.idprocs(iproc)) THEN
            IF (nrlocout(iset).GT.0) THEN
               j2 = j2 - nr1loc + 1
               j1 = j2 + (1-nrlocout(iset))*ylims(3)
               j3 = ylims(3)
            ELSE
               j1 = 1; j2 = 0; j3 = 1
            ENDIF
            limloc(:,2,iset) = (/j1,j2,j3/)
         ENDIF

      ENDDO iproc_321

   ENDIF

ENDDO iset_320
            
!
!3.3 Station output
!------------------
!
!3.3.1 Number of local stations
!----------------------------
!

iset_331: DO iset=1,nosetstsr
   IF (tsrgpars(iset)%nodim.GT.0.AND.(.NOT.tsrgpars(iset)%gridded)) THEN
      istat_3311: DO istat=1,tsrgpars(iset)%nostats
         l = lstatstsr(iset,istat)
         i = tsrstatlocs(l)%ipos; j = tsrstatlocs(l)%jpos
         iproc_33111: DO iproc=1,nprocs
            IF (local_proc(i,j,iproc=iproc)) THEN
               nostatsprocs(iproc,iset) = nostatsprocs(iproc,iset) + 1
               IF (idloc.EQ.idprocs(iproc)) THEN
                  nostatsloc(iset) = nostatsloc(iset) + 1
               ENDIF
            ENDIF
         ENDDO iproc_33111
      ENDDO istat_3311
   ENDIF
ENDDO iset_331

!
!3.3.2 Allocate
!--------------
!

maxstats = MAXVAL(nostatsprocs)
ALLOCATE (lstatsprocs(maxstats,nprocs,nosetstsr),STAT=errstat)
CALL error_alloc('lstatsprocs',3,(/maxstats,nprocs,nosetstsr/),kndint)
IF (maxstats.GT.0) lstatsprocs = 0

maxstatsloc = MAXVAL(nostatsloc)
ALLOCATE (istatpos(maxstatsloc,nosetstsr),STAT=errstat)
CALL error_alloc('istatpos',2,(/maxstatsloc,nosetstsr/),kndint)
IF (maxstatsloc.GT.0) istatpos = 0
ALLOCATE (jstatpos(maxstatsloc,nosetstsr),STAT=errstat)
CALL error_alloc('jstatpos',2,(/maxstatsloc,nosetstsr/),kndint)
IF (maxstatsloc.GT.0) jstatpos = 0

!
!3.3.3 Define
!------------
!

iproc_333: DO iproc=1,nprocs
iset_333: DO iset=1,nosetstsr
   IF (tsrgpars(iset)%nodim.GT.0.AND.(.NOT.tsrgpars(iset)%gridded)) THEN
      lloc = 0
      istat_3331: DO istat=1,tsrgpars(iset)%nostats
         l = lstatstsr(iset,istat)
         i = tsrstatlocs(l)%ipos; j = tsrstatlocs(l)%jpos
         IF (local_proc(i,j,iproc=iproc)) THEN
            lloc = lloc + 1
            lstatsprocs(lloc,iproc,iset) = l
            IF (idloc.EQ.idprocs(iproc)) THEN
               istatpos(lloc,iset) = i - nc1loc + 1
               jstatpos(lloc,iset) = j - nr1loc + 1
            ENDIF
         ENDIF
      ENDDO istat_3331
   ENDIF
ENDDO iset_333
ENDDO iproc_333

!
!3.4 Number of wet points
!------------------------
!

iset_340: DO iset=1,nosetstsr
   IF (tsrgpars(iset)%packing) THEN
      imax = lim_dims(limloc(:,1,iset))
      jmax = lim_dims(limloc(:,2,iset))
      ALLOCATE (maskvals(imax,jmax),STAT=errstat)
      CALL error_alloc('maskvals',2,(/imax,jmax/),kndlog)
      maskvals = mask_array(imax,jmax,limloc(:,:,iset),'C')
      nowetlocout = COUNT(maskvals)
      CALL sum_vars(nowetlocout,nowetout,0,commall=.TRUE.)
      tsrgpars(iset)%nowetout = nowetout   
      DEALLOCATE (maskvals)
   ELSE
      tsrgpars(iset)%nowetout = 0
   ENDIF
ENDDO iset_340

!
!4. Time specifiers
!-------------------
!

ALLOCATE (ntout(nosetstsr),STAT=errstat)
CALL error_alloc('ntout',1,(/nosetstsr/),kndint)
ntout = 0

!
!5. Abort if necessary
!----------------------
!

CALL error_abort('time_series_init',ierrno_inival)

!
!6. Write file headers and initial arrays
!----------------------------------------
!

CALL inout_atts_out('TS')

CALL log_timer_out()


RETURN

END SUBROUTINE time_series_init

!========================================================================

SUBROUTINE time_series_grid
!************************************************************************
!
! *time_series_grid* Write output grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Time_Series.f90  V2.10.1
!
! Description -
!
! Reference -
!
! Calling program - time_series
!
! External calls -
!
! Module calls - combine_mod, combine_write_stats_loc, combine_write_submod,
!                date_to_year, diff_dates, error_alloc, global_mask,
!                loop_index, write_vars, write_time, Zcoord_var
!
!************************************************************************
!
!1. Log info
!-----------
!

flag = .FALSE.
iset_110: DO iset=1,nosetstsr
   IF (loop_index(tsrgpars(iset)%tlims,nt)) THEN
      flag = .TRUE.
   ENDIF
ENDDO iset_110
IF (flag) THEN
   procname(pglev+1) = 'time_series_grid'
   CALL log_timer_in()
ELSE
   RETURN
ENDIF

!
!2. First output time
!--------------------
!

iset_200: DO iset=1,nosetstsr

   nodim = tsrgpars(iset)%nodim
   IF (nodim.GT.0.AND.nt.EQ.tsrgpars(iset)%tlims(1)) THEN

      gridded = tsrgpars(iset)%gridded
      packing = tsrgpars(iset)%packing
      time_grid = tsrgpars(iset)%time_grid.AND.nodim.EQ.3
      ncout = tsrgpars(iset)%ncout
      nrout = tsrgpars(iset)%nrout
      nzout = tsrgpars(iset)%nzout
      nostats = tsrgpars(iset)%nostats
      nowetout = tsrgpars(iset)%nowetout
      vcoord = tsrgpars(iset)%vcoord
      fill_value = MERGE(double_fill,0.0_kndrlong,packing)
      xlimsglb = tsrgpars(iset)%xlims
      ylimsglb = tsrgpars(iset)%ylims
      xlims = limloc(:,1,iset)
      ylims = limloc(:,2,iset)
      zlims = tsrgpars(iset)%zlims
      lprocs = limprocs(:,:,:,iset)
      fill_value = double_fill

!
!2.1 Allocate arrays
!-------------------
!
!2.1.1 Land mask arrays
!----------------------
!      

      IF (packing) THEN
         ALLOCATE(maskglb(nc,nr),STAT=errstat)
         CALL error_alloc('maskglb',2,(/nc,nr/),kndlog)
         ALLOCATE(lmask(nowetout),STAT=errstat)
         CALL error_alloc('lmask',1,(/nowetout/),kndint)
      ENDIF
 
!
!2.1.2 Gridded
!-------------
!

      IF (gridded) THEN
         
!        ---horizontal grid
         ALLOCATE (xout(ncout,nrout),STAT=errstat)
         CALL error_alloc('xout',2,(/ncout,nrout/),kndrtype)
         ALLOCATE (yout(ncout,nrout),STAT=errstat)
         CALL error_alloc('yout',2,(/ncout,nrout/),kndrtype)
         
!        ---bathymetry
         ALLOCATE (depout(nclocout(iset),nrlocout(iset)),STAT=errstat)
         CALL error_alloc('depout',2,(/nclocout(iset),nrlocout(iset)/),kndrtype)
         IF (SIZE(depout).GT.0) depout = fill_value

!        ---vertical grid
         IF (nodim.EQ.3) THEN
            IF (vcoord.EQ.1.AND.(.NOT.time_grid)) THEN
               ALLOCATE (zout(nclocout(iset),nrlocout(iset),nzout),STAT=errstat)
               CALL error_alloc('zout',3,&
                             & (/nclocout(iset),nrlocout(iset),nzout/),kndrtype)
               IF (SIZE(zout).GT.0) zout = fill_value
            ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.LE.2) THEN
               ALLOCATE (gsigout(nzout),STAT=errstat)
               CALL error_alloc('gsigout',1,(/nzout/),kndrtype)
            ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
               ALLOCATE (gsout(nclocout(iset),nrlocout(iset),nzout),&
                       & STAT=errstat)
               CALL error_alloc('gsout',3,&
                             & (/nclocout(iset),nrlocout(iset),nzout/),kndrtype)
               IF (SIZE(gsout).GT.0) gsout = fill_value
            ENDIF
         ENDIF
 
!
!2.1.3 Non-gridded
!-----------------      
!
         
      ELSE

!        ---horizontal grid         
         ALLOCATE (xout(nostats,1),STAT=errstat)
         CALL error_alloc('xout',2,(/nostats,1/),kndrtype)
         ALLOCATE (yout(nostats,1),STAT=errstat)
         CALL error_alloc('yout',2,(/nostats,1/),kndrtype)

!        ---bathymetry         
         ALLOCATE (depout(nostatsloc(iset),1),STAT=errstat)
         CALL error_alloc('depout',2,(/nostatsloc(iset),1/),kndrtype)

!        ---vertical grid
         IF (nodim.EQ.3) THEN
            IF (vcoord.EQ.1.AND.(.NOT.time_grid)) THEN
               ALLOCATE (zout(nostatsloc(iset),1,nzout),STAT=errstat)
               CALL error_alloc('zout',3,(/nostatsloc(iset),1,nzout/),kndrtype)
            ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.LE.2) THEN
               ALLOCATE (gsigout(nzout),STAT=errstat)
               CALL error_alloc('gsigout',1,(/nzout/),kndrtype)
            ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
               ALLOCATE (gsout(nostatsloc(iset),1,nzout),STAT=errstat)
               CALL error_alloc('gsout',3,(/nostatsloc(iset),1,nzout/),kndrtype)
            ENDIF
         ENDIF

      ENDIF

!
!2.2 Output grid
!---------------
!
!2.2.1 Index vecor of sea points
!-------------------------------         
!

      IF (packing) THEN
         CALL global_mask(maskglb,'C  ',1)            
         l = 0
         j_221: DO j=ylimsglb(1),ylimsglb(2),ylimsglb(3)
         i_221: DO i=xlimsglb(1),xlimsglb(2),xlimsglb(3)
            IF (maskglb(i,j)) THEN 
               l = l + 1
               lmask(l) = (j-1)*ncout + i
            ENDIF
         ENDDO i_221
         ENDDO j_221
      ENDIF

!         
!2.2.2 Gridded
!-------------
!      

      IF (gridded) THEN

!        ---horizontal grid         
         j_2221: DO j=ylimsglb(1),ylimsglb(2),ylimsglb(3)
         i_2221: DO i=xlimsglb(1),xlimsglb(2),xlimsglb(3)
            ii = (i-xlimsglb(1))/xlimsglb(3) + 1
            jj = (j-ylimsglb(1))/ylimsglb(3) + 1
            xout(ii,jj) = gxcoordglbatc(i,j)
            yout(ii,jj) = gycoordglbatc(i,j)
         ENDDO i_2221
         ENDDO j_2221
            
!        ---bathymetry and vertical grid         
         j_2222: DO j=ylims(1),ylims(2),ylims(3)
         i_2222: DO i=xlims(1),xlims(2),xlims(3)
            ii = (i-xlims(1))/xlims(3) + 1
            jj = (j-ylims(1))/ylims(3) + 1
            depout(ii,jj) = MERGE(depmeanatc(i,j),fill_value,maskatc_int(i,j))
            IF (nodim.EQ.3) THEN
               IF (maskatc_int(i,j)) THEN
                  k_22221: DO k=zlims(1),zlims(2),zlims(3)
                     kk = (k-zlims(1))/zlims(3) + 1
                     IF (vcoord.EQ.1.AND.(.NOT.time_grid)) THEN
                        zout(ii,jj,kk) = Zcoord_var(i,j,k,'C  ',.TRUE.)
                     ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
                        gsout(ii,jj,kk) = gscoordatc(i,j,k) - 1.0
                     ENDIF
                  ENDDO k_22221
               ENDIF
            ENDIF
         ENDDO i_2222
         ENDDO j_2222

!
!2.2.3 Non-gridded
!-----------------
!         

      ELSE

!        ---horizontal grid
         l_2231: DO l=1,nostats
            i = tsrstatlocs(l)%ipos; j = tsrstatlocs(l)%jpos
            xout(l,1) = gxcoordglbatc(i,j); yout(l,1) = gycoordglbatc(i,j)
         ENDDO l_2231
         
!        ---bathymetry and vertical grid 
         l_2232: DO l=1,nostatsloc(iset)
            i = istatpos(l,iset); j = jstatpos(l,iset)
            depout(l,1) = depmeanatc(i,j)
            IF (nodim.EQ.3) THEN
               k_22321: DO k=zlims(1),zlims(2),zlims(3)
                  kk = (k-zlims(1))/zlims(3) + 1
                  IF (vcoord.EQ.1.AND.(.NOT.time_grid)) THEN
                     zout(l,1,kk) = Zcoord_var(i,j,k,'C  ',.TRUE.)
                  ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
                     gsout(l,1,kk) = gscoordatc(i,j,k) - 1.0
                  ENDIF
               ENDDO k_22321
            ENDIF
         ENDDO l_2232

      ENDIF

!
!2.2.4 Sigma coordinate
!----------------------
!      

      IF (nodim.EQ.3.AND.vcoord.EQ.2.AND.iopt_grid_vtype.LE.2) THEN
         k_224: DO k=zlims(1),zlims(2),zlims(3)
            kk = (k-zlims(1))/zlims(3) + 1
            gsigout(kk)= gsigcoordatc(k) - 1.0
         ENDDO k_224
      ENDIF
      
!
!2.3 Write
!---------
!

      l_230: DO l=2,3
         SELECT CASE (l)
            CASE (2)
               filepars = tsr2d(iset)
            CASE (3)
               filepars = tsr3d(iset)
         END SELECT
         IF (.NOT.filepars%defined) CYCLE l_230
         ivar = 1
            
!
!2.3.1 Gridded
!--------------            
!

         IF (gridded) THEN

!              ---horizontal grid
            IF (master) THEN
               IF (iopt_grid_htype.LE.2) THEN
                  CALL write_vars(xout(:,1),filepars,ivar+1)
                  CALL write_vars(yout(1,:),filepars,ivar+2)
               ELSE
                  CALL write_vars(xout,filepars,ivar+1)
                  CALL write_vars(yout,filepars,ivar+2)
               ENDIF
            ENDIF
            ivar = ivar + 2

!              ---index vector of sea points               
            IF (packing) THEN
               IF (master) CALL write_vars(lmask,filepars,ivar+1)
               ivar = ivar + 1
            ENDIF

!           ---bathymetry
            CALL combine_write_submod(depout,filepars,ivar+1,&
                                    & (/ncout,nrout/),lprocs,nowet=nowetout)
            ivar = ivar + 1

!           ---vertical grid
            IF (filepars%nodim.EQ.3) THEN
               IF (vcoord.EQ.1.AND.(.NOT.time_grid)) THEN
                  CALL combine_write_submod(zout,filepars,ivar+1,&
                                  & (/ncout,nrout,nzout/),lprocs,nowet=nowetout)
                  ivar = ivar + 1
               ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
                  CALL combine_write_submod(gsout,filepars,ivar+1,&
                                  & (/ncout,nrout,nzout/),lprocs,nowet=nowetout)
                  ivar = ivar + 1
               ENDIF
            ENDIF
            
!
!2.3.2 Non-gridded
!-----------------            
!

               
         ELSE
!           ---horizontal grid
            IF (master) THEN
               CALL write_vars(xout(:,1),filepars,ivar+1)
               CALL write_vars(yout(:,1),filepars,ivar+2)
            ENDIF
            ivar = ivar + 2
!           ---bathymetry
            CALL combine_write_stats_loc(depout(:,1),filepars,&
                                & ivar+1,maxstats,nostats,nostatsprocs(:,iset),&
                                & lstatsprocs(:,:,iset))
            ivar = ivar + 1
!           ---vertical grid
            IF (filepars%nodim.EQ.3) THEN               
               IF (vcoord.EQ.1.AND.(.NOT.time_grid)) THEN
                  CALL combine_write_stats_loc(zout(:,1,:),filepars,&
                                & ivar+1,maxstats,nostats,nostatsprocs(:,iset),&
                                & lstatsprocs(:,:,iset),reduced=.FALSE.)
                  ivar = ivar + 1
               ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
                  CALL combine_write_stats_loc(gsout(:,1,:),filepars,&
                                & ivar+1,maxstats,nostats,nostatsprocs(:,iset),&
                                & lstatsprocs(:,:,iset),reduced=.FALSE.)
                  ivar = ivar + 1
               ENDIF
            ENDIF

         ENDIF

!
!2.3.3 Sigma coordinate
!----------------------
!
            
         IF (filepars%nodim.EQ.3.AND.vcoord.EQ.2.AND.&
           & iopt_grid_vtype.LE.2) THEN
            IF (master) CALL write_vars(gsigout,filepars,ivar+1)
            ivar = ivar + 1
         ENDIF

      ENDDO l_230

!
!2.4 Deallocate
!--------------
!

      DEALLOCATE (xout,yout,depout)
      IF (packing) DEALLOCATE (maskglb,lmask)
      IF (nodim.EQ.3) THEN
         IF (vcoord.EQ.1.AND.(.NOT.time_grid)) THEN
            DEALLOCATE (zout)
         ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.LE.2) THEN
            DEALLOCATE (gsigout)
         ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
            DEALLOCATE (gsout)
         ENDIF
      ENDIF

   ENDIF
      
ENDDO iset_200

!
!3. Time grid
!------------
!

iset_300: DO iset=1,nosetstsr

   IF (loop_index(tsrgpars(iset)%tlims,nt)) THEN
      time_format = tsrgpars(iset)%time_format
      ncout = tsrgpars(iset)%ncout
      nrout = tsrgpars(iset)%nrout
      nzout = tsrgpars(iset)%nzout
      vcoord = tsrgpars(iset)%vcoord
      fill_value = MERGE(double_fill,0.0_kndrlong,packing)
      time_grid = tsrgpars(iset)%time_grid.AND.tsrgpars(iset)%nodim.EQ.3
      refdate = tsrgpars(iset)%refdate

!
!3.1 Time coordinate
!-------------------
!

      IF (time_format.EQ.8) THEN
         CALL date_to_year(CDateTime,rtime)
      ELSEIF (time_format.GT.0) THEN
         CALL diff_dates(refdate,CDateTime,time_format,rtime=rtime)
      ENDIF

!
!3.2 Elevations
!--------------
!

      IF (time_grid) THEN
         IF (tsrgpars(iset)%gridded) THEN
            IF (vcoord.EQ.1) THEN
               ALLOCATE (zout(nclocout(iset),nrlocout(iset),nzout),STAT=errstat)
               CALL error_alloc('zout',3,&
                    & (/nclocout(iset),nrlocout(iset),nzout/),kndrtype)
               IF (SIZE(zout).GT.0) zout = fill_value
            ELSE
               ALLOCATE (zetout(nclocout(iset),nrlocout(iset)),STAT=errstat)
               CALL error_alloc('zetout',2,(/nclocout(iset),nrlocout(iset)/),&
                    & kndrtype)
               IF (SIZE(zetout).GT.0) zetout = fill_value
            ENDIF
            xlims = limloc(:,1,iset); ylims = limloc(:,2,iset)
            zlims = tsrgpars(iset)%zlims
            i_321: DO i=xlims(1),xlims(2),xlims(3)
            j_321: DO j=ylims(1),ylims(2),ylims(3)
               ii = (i-xlims(1))/xlims(3) + 1
               jj = (j-ylims(1))/ylims(3) + 1
               IF (vcoord.EQ.1) THEN
                  IF (maskatc_int(i,j)) THEN
                     k_3211: DO k=zlims(1),zlims(2),zlims(3)
                        kk = (k-zlims(1))/zlims(3) + 1
                        zout(ii,jj,kk) = Zcoord_var(i,j,k,'C  ',.FALSE.)
                     ENDDO k_3211
                  ELSE
                     zout(ii,jj,:) = 0.0
                  ENDIF
               ELSE
                  zetout(ii,jj) = MERGE(zeta(i,j),0.0,maskatc_int(i,j))        
               ENDIF
            ENDDO j_321
            ENDDO i_321
         ELSE
            IF (vcoord.EQ.1) THEN
               ALLOCATE (zout(nostatsloc(iset),1,nzout),STAT=errstat)
               CALL error_alloc('zout',3,(/nostatsloc,1,nzout/),kndrtype)
            ELSE
               ALLOCATE (zetout(nostatsloc(iset),1),STAT=errstat)
               CALL error_alloc('zetout',2,(/nostatsloc,1/),kndrtype)
            ENDIF
            l_322: DO l=1,nostatsloc(iset)
               i = istatpos(l,iset); j = jstatpos(l,iset)
               IF (vcoord.EQ.1) THEN
                  k_3221: DO k=zlims(1),zlims(2),zlims(3)
                     kk = (k-zlims(1))/zlims(3) + 1
                     zout(l,1,k) = Zcoord_var(i,j,k,'C  ',.FALSE.)
                  ENDDO k_3221
               ELSE
                  zetout(l,1) = zeta(i,j)
               ENDIF
            ENDDO l_322
         ENDIF
      ENDIF

!
!3.3 Write
!---------
!

      l_330: DO l=0,3
         SELECT CASE (l)
            CASE (0)
               filepars = tsr0d(iset)
            CASE (2)
               filepars = tsr2d(iset)
            CASE (3)
               filepars = tsr3d(iset)
         END SELECT
         IF (.NOT.filepars%defined.OR.l.EQ.1) CYCLE l_330
!        ---time coordinate
         IF (time_format.EQ.0) THEN
            CALL write_time(CDateTime,filepars)
         ELSE
            CALL write_time(rtime,filepars)
         ENDIF
!        ---elevations
         IF (time_grid.AND.filepars%nodim.EQ.3) THEN
            zvarid = filepars%zvarid
            nostats = tsrgpars(iset)%nostats
            nowetout = tsrgpars(iset)%nowetout
            lprocs = limprocs(:,:,:,iset)
            IF (tsrgpars(iset)%gridded) THEN
               IF (vcoord.EQ.1) THEN
                  CALL combine_write_submod(zout,filepars,zvarid,&
                                  & (/ncout,nrout,nzout/),lprocs,nowet=nowetout)
               ELSE
                  CALL combine_write_submod(zetout,filepars,zvarid,&
                                  & (/ncout,nrout/),lprocs,nowet=nowetout)
               ENDIF
            ELSE
               IF (vcoord.EQ.1) THEN
                  CALL combine_write_stats_loc(zout(:,1,:),filepars,zvarid,&
                                   & maxstats,nostats,nostatsprocs(:,iset),&
                                   & lstatsprocs(:,:,iset))
               ELSE
                  CALL combine_write_stats_loc(zetout(:,1),filepars,zvarid,&
                                   & maxstats,nostats,nostatsprocs(:,iset),&
                                   & lstatsprocs(:,:,iset))
               ENDIF
            ENDIF
            IF (vcoord.EQ.1) THEN
               DEALLOCATE (zout)
            ELSE
               DEALLOCATE (zetout)
            ENDIF
         ENDIF
!        ---store changes in file structure
         SELECT CASE (l)
            CASE (0)
               tsr0d(iset) = filepars
            CASE (2)
               tsr2d(iset) = filepars
            CASE (3)
               tsr3d(iset) = filepars
         END SELECT
      ENDDO l_330
      
   ENDIF

ENDDO iset_300

CALL log_timer_out()


RETURN

END SUBROUTINE time_series_grid


END SUBROUTINE time_series
