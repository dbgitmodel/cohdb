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

SUBROUTINE time_averages
!************************************************************************
!
! *time_averages* Time averaged output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Time_Averages.f90  V2.11
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
! External calls - read_cif_params
!
! Internal calls - time_averages_grid, time_averages_init, time_averages_reset,
!                  time_averages_update, time_averages_write
!
! Module calls  - close_filepars, error_alloc
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
USE time_routines, ONLY: add_secs_to_date, date_to_year, diff_dates, &
                       & log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index

IMPLICIT none

!
!* Local variables
!
LOGICAL, SAVE :: last_call
LOGICAL :: first, flag, fmask, gridded, newvals, oldvals, packing, second
CHARACTER (LEN=lentime) :: ctime, refdate
INTEGER, SAVE :: maxstats, novals0dnew, novals0dold, novals2dnew, &
               & novals2dold, novals3dnew, novals3dold, novars0d, novars2d, &
               & novars3d
INTEGER :: i, ii, iivar, imax, iproc, iset, istat, ivar, i0, i1, i2, i3, j, jj, &
         & jmax, j1, j2, j3, k, kk, l, lbuf0d, lbuf2d, lbuf3d, lloc, lnew, &
         & lold, l0dat, l2dat, l3dat, maxstatsloc, maxvars, ncout, nocoords0d, &
         & nocoords2d, nocoords3d, nodat, nodim, nostats, noutvars0d, &
         & noutvars2d, noutvars3d, nowetlocout, nowetout, npcc, nrank, nrout, &
         & nzout, time_format, vcoord
REAL :: fill_value
REAL (KIND=kndrlong) :: rtime
TYPE (FileParams) :: filepars

INTEGER, DIMENSION(3) :: tlims, xlims, xlimsglb, ylims, ylimsglb, zlims
INTEGER, DIMENSION(novarsavr) :: vecids
INTEGER, DIMENSION(2,2,nprocs) :: lprocs
REAL, DIMENSION(novarsavr) :: outdat, outsub

LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb, maskvals
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: lmask, navrtot, nclocout, ndat0d, &
                                          & ndat2d, ndat3d, nostatsloc, &
                                          & nrlocout, ntout
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: istatpos, ivarsavr0d, &
                                            & ivarsavr2d, ivarsavr3d, &
                                            & jstatpos, nostatsprocs
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: limloc, lstatsprocs
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: limprocs
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: buffold, buff0d, buff2d, buff3d, gsigout
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: depout, out2dsub, xout, yout
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: gsout, out3dsub, zout
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: out4dsub
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: outvars

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*last_call*   LOGICAL .TRUE. if the routine has been called for the last time
!*istatpos*    INTEGER X-index of local output station
!*ivarsavr0d*  INTEGER Array of indices (per file) of 0-D output variables
!*ivarsavr2d*  INTEGER Array of indices (per file) of 2-D output variables
!*ivarsavr3d*  INTEGER Array of indices (per file) of 3-D output variables
!*jstatpos*    INTEGER Y-index of local output station
!*limloc*      INTEGER Start/End/Increment values of local (regular) output
!                      grid in X- and Y-directions
!*limprocs*    INTEGER Start/End indices of local array section within global
!                      array (regular grid)
!*lstatprocs*  INTEGER Output labels of local stations for each process domain
!*navrtot*     INTEGER Number of time steps within one averaging period (+1)
!*nclocout*    INTEGER X-dimension of local (regular) output grid
!*ndat0d*      INTEGER Total number (per file) of 0-D data
!*ndat2d*      INTEGER Total number (per file) of 2-D data
!*ndat3d*      INTEGER Total number (per file) of 3-D data
!*nostatsloc*  INTEGER Number of local output stations (per file)
!*nostatprocs* INTEGER Number of local output stations for each process domain 
!*novals0dnew* INTEGER Total number (over all files) of 0-D data at new time
!*novals0dold* INTEGER Total number (over all files) of 0-D data at old time
!*novals2dnew* INTEGER Total number (over all files) of 2-D data at new time
!*novals2dold* INTEGER Total number (over all files) of 2-D data at old time
!*novals3dnew* INTEGER Total number (over all files) of 3-D data at new time
!*novals3dold* INTEGER Total number (over all files) of 3-D data at old time
!*novars0d*    INTEGER Total number of 0-D output variables
!*novars2d*    INTEGER Total number of 2-D output variables
!*novars3d*    INTEGER Total number of 3-D output variables
!*nroutloc*    INTEGER Y-dimension of local (regular) output grid
!*ntout*       INTEGER Output time counter
!*buff0d*      REAL    Buffer for temporary storage of 0-D data
!*buff2d*      REAL    Buffer for temporary storage of 2-D data
!*buff3d*      REAL    Buffer for temporary storage of 3-D data
!*depout*      REAL    Output depth array                                   [m]
!*xout*        REAL    Output X-coordinates                  [degree_east or m]
!*yout*        REAL    Output Y-coordinates                 [degree_north or m]
!*zout*        REAL    Output Z-coordinates                                 [m]
!
!------------------------------------------------------------------------------
!

IF (nt.EQ.0) last_call=.FALSE.
IF (last_call) RETURN

procname(pglev+1) = 'time_averages'
CALL log_timer_in(npcc)

!
!1. Initialisation on first call
!-------------------------------
!    

IF (nt.EQ.0) CALL time_averages_init

!
!2. Write output grid
!--------------------
!

CALL time_averages_grid

!
!3. Create buffers for data storage
!----------------------------------
!
!3.1 Buffer sizes
!----------------
!

iset_310: DO iset=1,nosetsavr
   tlims = avrgpars(iset)%tlims
   IF ((nt.GE.tlims(1)).AND.(nt.LE.tlims(2))) THEN
      IF (avr0d(iset)%defined) ndat0d(iset) = avr0d(iset)%novars
      IF (avrgpars(iset)%gridded) THEN
         IF (avr2d(iset)%defined) THEN
            ndat2d(iset) = avr2d(iset)%novars*nclocout(iset)*nrlocout(iset)
         ENDIF
         IF (avr3d(iset)%defined) THEN
            ndat3d(iset) = avr3d(iset)%novars*nclocout(iset)*nrlocout(iset)*&
                         & avrgpars(iset)%nzout
         ENDIF
      ELSE
         IF (avr2d(iset)%defined) THEN
            ndat2d(iset) = avr2d(iset)%novars*nostatsloc(iset)
         ENDIF
         IF (avr3d(iset)%defined) THEN
            ndat3d(iset) = avr3d(iset)%novars*avrgpars(iset)%nzout*&
                         & nostatsloc(iset)
         ENDIF
      ENDIF
   ENDIF
ENDDO iset_310

!
!3.2 Refresh buffers (if necessary)
!----------------------------------
!
!3.2.1 0-D buffer
!----------------
!

novals0dnew = SUM(ndat0d)
IF ((novals0dold.EQ.0).AND.(novals0dnew.GT.0)) THEN
   ALLOCATE (buff0d(novals0dnew),STAT=errstat)
   CALL error_alloc('buff0d',1,(/novals0dnew/),kndrtype)
ELSEIF ((novals0dold.GT.0).AND.(novals0dnew.EQ.0)) THEN
   DEALLOCATE (buff0d)
ELSEIF (novals0dold.NE.novals0dnew) THEN
   ALLOCATE (buffold(novals0dold),STAT=errstat)
   CALL error_alloc('buffold',1,(/novals0dold/),kndrtype)
   buffold = buff0d
   DEALLOCATE (buff0d)
   ALLOCATE (buff0d(novals0dnew),STAT=errstat)
   CALL error_alloc('buff0d',1,(/novals0dnew/),kndrtype)

   lold = 0; lnew = 0
   iset_321: DO iset=1,nosetsavr
      tlims = avrgpars(iset)%tlims
      IF (avr0d(iset)%defined) THEN
         nodat = ndat0d(iset)
         oldvals = ((nt-ic3d).GE.tlims(1)).AND.((nt-ic3d).LE.tlims(2))
         newvals = (nt.GE.tlims(1)).AND.(nt.LE.tlims(2))
         IF (oldvals.AND.newvals) THEN
            buff0d(lnew+1:lnew+nodat) = buffold(lold+1:lold+nodat)
            lold = lold + nodat; lnew = lnew + nodat 
         ELSEIF (oldvals.AND.(.NOT.newvals)) THEN
            lold = lold + nodat
         ELSEIF ((.NOT.oldvals).AND.newvals) THEN
            lnew = lnew + nodat
         ENDIF
      ENDIF
   ENDDO iset_321

   DEALLOCATE (buffold)

ENDIF

novals0dold = novals0dnew

!
!3.2.2 2-D buffer
!----------------
!

novals2dnew = SUM(ndat2d)
IF ((novals2dold.EQ.0).AND.(novals2dnew.GT.0)) THEN
   ALLOCATE (buff2d(novals2dnew),STAT=errstat)
   CALL error_alloc('buff2d',1,(/novals2dnew/),kndrtype)
ELSEIF ((novals2dold.GT.0).AND.(novals2dnew.EQ.0)) THEN
   DEALLOCATE (buff2d)
ELSEIF (novals2dold.NE.novals2dnew) THEN
   ALLOCATE (buffold(novals2dold),STAT=errstat)
   CALL error_alloc('buffold',1,(/novals2dold/),kndrtype)
   buffold = buff2d
   DEALLOCATE (buff2d)
   ALLOCATE (buff2d(novals2dnew),STAT=errstat)
   CALL error_alloc('buff2d',1,(/novals2dnew/),kndrtype)

   lold = 0; lnew = 0
   iset_322: DO iset=1,nosetsavr
      tlims = avrgpars(iset)%tlims
      IF (avr2d(iset)%defined) THEN
         nodat = ndat2d(iset)
         oldvals = ((nt-ic3d).GE.tlims(1)).AND.((nt-ic3d).LE.tlims(2))
         newvals = (nt.GE.tlims(1)).AND.(nt.LE.tlims(2))
         IF (oldvals.AND.newvals) THEN
            buff2d(lnew+1:lnew+nodat) = buffold(lold+1:lold+nodat)
            lold = lold + nodat; lnew = lnew + nodat 
         ELSEIF (oldvals.AND.(.NOT.newvals)) THEN
            lold = lold + nodat
         ELSEIF ((.NOT.oldvals).AND.newvals) THEN
            lnew = lnew + nodat
         ENDIF
      ENDIF
   ENDDO iset_322

   DEALLOCATE (buffold)

ENDIF

novals2dold = novals2dnew

!
!3.2.3 3-D buffer
!----------------
!

novals3dnew = SUM(ndat3d)
IF ((novals3dold.EQ.0).AND.(novals3dnew.GT.0)) THEN
   ALLOCATE (buff3d(novals3dnew),STAT=errstat)
   CALL error_alloc('buff3d',1,(/novals3dnew/),kndrtype)
ELSEIF ((novals3dold.GT.0).AND.(novals3dnew.EQ.0)) THEN
   DEALLOCATE (buff3d)
ELSEIF (novals3dold.NE.novals3dnew) THEN
   ALLOCATE (buffold(novals3dold),STAT=errstat)
   CALL error_alloc('buffold',1,(/novals3dold/),kndrtype)
   buffold = buff3d
   DEALLOCATE (buff3d)
   ALLOCATE (buff3d(novals3dnew),STAT=errstat)
   CALL error_alloc('buff3d',1,(/novals3dnew/),kndrtype)

   lold = 0; lnew = 0
   iset_323: DO iset=1,nosetsavr
      tlims = avrgpars(iset)%tlims
      IF (avr3d(iset)%defined) THEN
         nodat = ndat3d(iset)
         oldvals = ((nt-ic3d).GE.tlims(1)).AND.((nt-ic3d).LE.tlims(2))
         newvals = (nt.GE.tlims(1)).AND.(nt.LE.tlims(2))
         IF (oldvals.AND.newvals) THEN
            buff3d(lnew+1:lnew+nodat) = buffold(lold+1:lold+nodat)
            lold = lold + nodat; lnew = lnew + nodat 
         ELSEIF (oldvals.AND.(.NOT.newvals)) THEN
            lold = lold + nodat
         ELSEIF ((.NOT.oldvals).AND.newvals) THEN
            lnew = lnew + nodat
         ENDIF
      ENDIF
   ENDDO iset_323

   DEALLOCATE (buffold)

ENDIF

novals3dold = novals3dnew

!
!4. Time-averaged output
!-----------------------
!
!4.1 Initialise buffers
!----------------------
!

first = .TRUE.; second = .FALSE.
CALL time_averages_reset

!
!4.2 Update buffers
!------------------
!

CALL time_averages_update

!
!4.3 Write time-averaged output
!------------------------------
!

CALL time_averages_write

!
!4.4 Reinitialise/update buffers
!-------------------------------
!

first = .FALSE.; second = .TRUE.
CALL time_averages_reset
CALL time_averages_update

!
!5. Finalise after last output
!-----------------------------
!
!5.1 Close data files
!--------------------
!

last_call = .TRUE.
iset_510: DO iset=1,nosetsavr
   IF (.NOT.cold_start.AND.nt.LT.avrgpars(iset)%tlims(2)) last_call = .FALSE.
   IF (master.AND.(cold_start.OR.nt.EQ.avrgpars(iset)%tlims(2))) THEN
      IF (avr0d(iset)%defined) CALL close_filepars(avr0d(iset))
      IF (avr2d(iset)%defined) CALL close_filepars(avr2d(iset))
      IF (avr3d(iset)%defined) CALL close_filepars(avr3d(iset))
   ENDIF
ENDDO iset_510

!
!5.2 Deallocate
!--------------
!

IF (last_call) THEN

   DEALLOCATE (ivarsavr0d,ivarsavr2d,ivarsavr3d,outvars)
   DEALLOCATE (nclocout,nrlocout,limprocs,limloc)
   DEALLOCATE (nostatsprocs,nostatsloc,lstatsprocs,istatpos,jstatpos)
   DEALLOCATE (ndat0d,ndat2d,ndat3d)
   DEALLOCATE (navrtot,ntout)

   IF (ALLOCATED(buff0d)) DEALLOCATE (buff0d)
   IF (ALLOCATED(buff2d)) DEALLOCATE (buff2d)
   IF (ALLOCATED(buff3d)) DEALLOCATE (buff3d)

ENDIF

CALL log_timer_out(npcc,itm_user)


RETURN


CONTAINS

!========================================================================

SUBROUTINE time_averages_init
!************************************************************************
!
! *time_averages_init* Initialise parameters for time averaging
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Time_Averages.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - time_averages
!
! External calls - read_cif_params, usrdef_avr_params
!
! Module calls - check_out_filepars, check_out_gpars, check_statslocs,
!   check_variables, combine_mod, default_out_files, default_out_gpars,
!   error_abort, error_alloc, error_alloc_struc, error_diff_vals_arrlist,
!   error_limits_arr, filepars_init, inout_atts_out, lim_dims, local_proc,
!   mask_array, outgpars_init, reset_out_files, reset_out_gpars,
!   reset_out_stats, reset_out_vars, statlocs_init, sum_vars, varatts_init,
!   warning_reset_arr_struc
!
!************************************************************************
!


procname(pglev+1) = 'time_averages_init'
CALL log_timer_in()

!
!1. Parameters for time averaging
!--------------------------------
!
!1.1 Allocate and initialise
!---------------------------
!
!---variables and vectors
ALLOCATE (avrvars(novarsavr),STAT=errstat)
CALL error_alloc_struc('avrvars',1,(/novarsavr/),'VariableAtts')
CALL varatts_init(avrvars)

!---file parameters
ALLOCATE (avr0d(nosetsavr),STAT=errstat)
CALL error_alloc_struc('avr0d',1,(/nosetsavr/),'FileParams')
CALL filepars_init(avr0d)
ALLOCATE (avr2d(nosetsavr),STAT=errstat)
CALL error_alloc_struc('avr2d',1,(/nosetsavr/),'FileParams')
CALL filepars_init(avr2d)
ALLOCATE (avr3d(nosetsavr),STAT=errstat)
CALL error_alloc_struc('avr3d',1,(/nosetsavr/),'FileParams')
CALL filepars_init(avr3d)

!---output specifiers
ALLOCATE (avrgpars(nosetsavr),STAT=errstat)
CALL error_alloc_struc('avrgpars',1,(/nosetsavr/),'OutFileParams')
CALL outgpars_init(avrgpars)

!---output stations
ALLOCATE (lstatsavr(nosetsavr,nostatsavr),STAT=errstat)
CALL error_alloc('lstatsavr',2,(/nosetsavr,nostatsavr/),kndint)
lstatsavr = 0
ALLOCATE (avrstatlocs(nostatsavr),STAT=errstat)
CALL error_alloc_struc('avrstatsloc',1,(/nostatsavr/),'StationLocs')
CALL statlocs_init(avrstatlocs)

!---variable indices
ALLOCATE (ivarsavr(nosetsavr,novarsavr),STAT=errstat)
CALL error_alloc('ivarsavr',2,(/nosetsavr,novarsavr/),kndint)
ivarsavr = 0

!
!1.2 Defaults
!------------
!
!---file parameters
CALL default_out_files(avr0d,'TA')
CALL default_out_files(avr2d,'TA')
CALL default_out_files(avr3d,'TA')

!---output specifiers
CALL default_out_gpars(avrgpars)

!
!1.3 Define
!---------
!

IF (ciffile%status.EQ.'R') THEN
   CALL read_cif_params(icif_avrgd)
ELSE
   CALL usrdef_avr_params
ENDIF

!
!1.4 Reset
!---------
!
!---file parameters
CALL reset_out_files(avr0d,'avr0d',avrgpars,.FALSE.)
CALL reset_out_files(avr2d,'avr2d',avrgpars,.FALSE.)
CALL reset_out_files(avr3d,'avr3d',avrgpars,.FALSE.)

!---variable atrributes
CALL reset_out_vars(avrvars)

!---output grid dimension
iset_141: DO iset=1,nosetsavr
   nodim = 0
   ivar_1411: DO ivar=1,novarsavr
      iivar = ivarsavr(iset,ivar)
      IF (iivar.GT.0) THEN
         IF (avrvars(iivar)%nrank.EQ.3) THEN
            nodim = 3; EXIT ivar_1411
         ENDIF
      ENDIF
   ENDDO ivar_1411
   IF (nodim.EQ.0) THEN
      ivar_1412: DO ivar=1,novarsavr
         iivar = ivarsavr(iset,ivar)
         IF (iivar.GT.0) THEN
            IF (avrvars(iivar)%nrank.EQ.2) THEN
               nodim = 2; EXIT ivar_1412
            ENDIF
         ENDIF
      ENDDO ivar_1412
   ENDIF
   avrgpars(iset)%nodim = nodim
   IF (iopt_grid_nodim.EQ.2.AND.avrgpars(iset)%nodim.EQ.3) THEN
      CALL warning_reset_arr_struc(avrgpars(iset)%nodim,'avrgpars','%nodim',2,&
                                 & 1,(/iset/))
   ENDIF
ENDDO iset_141

!---output grid parameters
CALL reset_out_gpars(avrgpars,'avrgpars',.FALSE.,.FALSE.)

!---station attributes
CALL reset_out_stats(avrstatlocs)

!---file parameters
avr0d%nodim = 0; avr2d%nodim = 2; avr3d%nodim = 3
iset_142: DO iset=1,nosetsavr
   IF (avrgpars(iset)%nodim.LT.3) THEN
      CALL warning_reset_arr_struc(avr3d(iset)%defined,'avr3d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
   IF (avrgpars(iset)%nodim.EQ.0) THEN
      CALL warning_reset_arr_struc(avr2d(iset)%defined,'avr2d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
ENDDO iset_142

!
!1.5 Check
!---------
!
!---output variables
CALL check_variables(avrvars,'avrvars')

!---file parameters
CALL check_out_filepars(avr0d,'avr0d')
CALL check_out_filepars(avr2d,'avr2d')
CALL check_out_filepars(avr3d,'avr3d')

!---output parameters
CALL check_out_gpars(avrgpars,'avrgpars',nostatsavr)

!---output stations
iset_151: DO iset=1,nosetsavr
   istat_1151: DO istat=1,avrgpars(iset)%nostats
      CALL error_limits_arr(lstatsavr(iset,istat),'lstatsavr',1,nostatsavr,2,&
                          & (/iset,istat/))
   ENDDO istat_1151
ENDDO iset_151
CALL check_statlocs(avrstatlocs,'avrstatlocs')

!---variable indices
iset_152: DO iset=1,nosetsavr
   ivar_1521: DO ivar=1,novarsavr
      CALL error_limits_arr(ivarsavr(iset,ivar),'ivarsavr',0,novarsavr,2,&
                          & (/iset,ivar/))
   ENDDO ivar_1521
   CALL error_diff_vals_arrlist(ivarsavr(iset,:),'ivarsavr',2,2,(/iset,1/),&
                             & .TRUE.)
ENDDO iset_152

!
!2. Variable indices
!-------------------
!
!2.1 Number of variables per dimension
!-------------------------------------
!

novars0d = COUNT(avrvars%nrank.EQ.0)
novars2d = COUNT(avrvars%nrank.EQ.2)
novars3d = COUNT(avrvars%nrank.EQ.3)

!
!2.2 Allocate and initialise
!---------------------------
!

ALLOCATE (ivarsavr0d(nosetsavr,novars0d),STAT=errstat)
CALL error_alloc('ivarsavr0d',2,(/nosetsavr,novars0d/),kndint)
IF (novars0d.GT.0) ivarsavr0d = 0
ALLOCATE (ivarsavr2d(nosetsavr,novars2d),STAT=errstat)
CALL error_alloc('ivarsavr2d',2,(/nosetsavr,novars2d/),kndint)
IF (novars2d.GT.0) ivarsavr2d = 0
ALLOCATE (ivarsavr3d(nosetsavr,novars3d),STAT=errstat)
CALL error_alloc('ivarsavr3d',2,(/nosetsavr,novars3d/),kndint)
IF (novars3d.GT.0) ivarsavr3d = 0
maxvars = MAX(novars0d,novars2d,novars3d)
ALLOCATE (outvars(maxvars),STAT=errstat)
CALL error_alloc_struc('outvars',1,(/maxvars/),'VariableAtts')

!
!2.3 Define
!----------
!

iset_230: DO iset=1,nosetsavr
   i0 = 0; i2 = 0; i3 = 0
   ivar_231: DO ivar=1,novarsavr
      l = ivarsavr(iset,ivar)
      IF (l.GT.0) THEN
         nrank = avrvars(l)%nrank
         SELECT CASE (nrank)
         CASE (0)
            i0 = i0 + 1
            ivarsavr0d(iset,i0) = l
         CASE (2)
            i2 = i2 + 1
            ivarsavr2d(iset,i2) = l
         CASE (3)
            i3 = i3 + 1
            ivarsavr3d(iset,i3) = l
         END SELECT
      ENDIF
   ENDDO ivar_231
   avr0d(iset)%novars = i0; avr2d(iset)%novars = i2; avr3d(iset)%novars = i3
   ivarsavr2d(iset,:) = ivarsavr2d(iset,:) - novars0d
   ivarsavr3d(iset,:) = ivarsavr3d(iset,:) - novars0d - novars2d
ENDDO iset_230

!
!2.4 Write output only if number of variables is greater than zero
!-----------------------------------------------------------------
!

iset_240: DO iset=1,nosetsavr
   IF (avr0d(iset)%novars.EQ.0) THEN
      CALL warning_reset_arr_struc(avr0d(iset)%defined,'avr0d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
   IF (avr2d(iset)%novars.EQ.0) THEN
      CALL warning_reset_arr_struc(avr2d(iset)%defined,'avr2d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
   IF (avr3d(iset)%novars.EQ.0) THEN
      CALL warning_reset_arr_struc(avr3d(iset)%defined,'avr3d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
ENDDO iset_240

!
!3. Space output specifiers
!--------------------------
!
!3.1 Allocate and initialise
!---------------------------
!
!---regular output
ALLOCATE (nclocout(nosetsavr),STAT=errstat)
CALL error_alloc('nclocout',1,(/nosetsavr/),kndint)
nclocout = 0
ALLOCATE (nrlocout(nosetsavr),STAT=errstat)
CALL error_alloc('nrlocout',1,(/nosetsavr/),kndint)
nrlocout = 0
ALLOCATE (limprocs(2,2,nprocs,nosetsavr),STAT=errstat)
CALL error_alloc('limprocs',4,(/2,2,nprocs,nosetsavr/),kndint)
limprocs = 0
ALLOCATE (limloc(3,2,nosetsavr),STAT=errstat)
CALL error_alloc('limloc',3,(/3,2,nosetsavr/),kndint)
limloc = 0

!---station output
ALLOCATE (nostatsprocs(nprocs,nosetsavr),STAT=errstat)
CALL error_alloc('nostatsprocs',2,(/nprocs,nosetsavr/),kndint)
nostatsprocs= 0
ALLOCATE (nostatsloc(nosetsavr),STAT=errstat)
CALL error_alloc('nostatsloc',1,(/nosetsavr/),kndint)
nostatsloc = 0

!
!3.2 Regular output
!------------------
!

iset_320: DO iset=1,nosetsavr

   IF (avrgpars(iset)%nodim.GT.0.AND.avrgpars(iset)%gridded) THEN

      avrgpars(iset)%ncout = lim_dims(avrgpars(iset)%xlims)
      avrgpars(iset)%nrout = lim_dims(avrgpars(iset)%ylims)
      avrgpars(iset)%nzout = lim_dims(avrgpars(iset)%zlims)

      iproc_321: DO iproc=1,nprocs

!        ---X-direction
         xlims = avrgpars(iset)%xlims
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
         ylims = avrgpars(iset)%ylims
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

iset_331: DO iset=1,nosetsavr
   IF (avrgpars(iset)%nodim.GT.0.AND.(.NOT.avrgpars(iset)%gridded)) THEN
      istat_3311: DO istat=1,avrgpars(iset)%nostats
         l = lstatsavr(iset,istat)
         i = avrstatlocs(l)%ipos; j = avrstatlocs(l)%jpos
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
ALLOCATE (lstatsprocs(nprocs,maxstats,nosetsavr),STAT=errstat)
CALL error_alloc('lstatsprocs',3,(/nprocs,maxstats,nosetsavr/),kndint)
IF (maxstats.GT.0) lstatsprocs = 0

maxstatsloc = MAXVAL(nostatsloc)
ALLOCATE (istatpos(maxstatsloc,nosetsavr),STAT=errstat)
CALL error_alloc('istatpos',2,(/maxstatsloc,nosetsavr/),kndint)
IF (maxstatsloc.GT.0) istatpos = 0
ALLOCATE (jstatpos(maxstatsloc,nosetsavr),STAT=errstat)
CALL error_alloc('jstatpos',2,(/maxstatsloc,nosetsavr/),kndint)
IF (maxstatsloc.GT.0) jstatpos = 0

!
!3.3.3 Define
!------------
!

iproc_333: DO iproc=1,nprocs
iset_333: DO iset=1,nosetsavr
   IF (avrgpars(iset)%nodim.GT.0.AND.(.NOT.avrgpars(iset)%gridded)) THEN
      lloc = 0
      istat_3331: DO istat=1,avrgpars(iset)%nostats
         l = lstatsavr(iset,istat)
         i = avrstatlocs(l)%ipos; j = avrstatlocs(l)%jpos
         IF (local_proc(i,j,iproc=iproc)) THEN
            lloc = lloc + 1
            lstatsprocs(iproc,lloc,iset) = l
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

iset_340: DO iset=1,nosetsavr
   IF (avrgpars(iset)%packing) THEN
      imax = lim_dims(limloc(:,1,iset))
      jmax = lim_dims(limloc(:,2,iset))
      ALLOCATE (maskvals(imax,jmax),STAT=errstat)
      CALL error_alloc('maskvals',2,(/imax,jmax/),kndlog)
      maskvals = mask_array(imax,jmax,limloc(:,:,iset),'C')
      nowetlocout = COUNT(maskvals)
      CALL sum_vars(nowetlocout,nowetout,0,commall=.TRUE.)
      avrgpars(iset)%nowetout = nowetout   
      DEALLOCATE (maskvals)
   ELSE
      avrgpars(iset)%nowetout = 0
   ENDIF
ENDDO iset_340

!
!3.5 Number of data
!------------------

ALLOCATE (ndat0d(nosetsavr),STAT=errstat)
CALL error_alloc('ndat0d',1,(/nosetsavr/),kndint)
ndat0d = 0
ALLOCATE (ndat2d(nosetsavr),STAT=errstat)
CALL error_alloc('ndat2d',1,(/nosetsavr/),kndint)
ndat2d = 0
ALLOCATE (ndat3d(nosetsavr),STAT=errstat)
CALL error_alloc('ndat3d',1,(/nosetsavr/),kndint)
ndat3d = 0

novals0dold = 0; novals2dold = 0; novals3dold = 0      

!
!4. Time specifiers
!------------------
!
!4.1 Allocate and initialise
!---------------------------
!

ALLOCATE (navrtot(nosetsavr),STAT=errstat)
CALL error_alloc('navrtot',1,(/nosetsavr/),kndint)
ALLOCATE (ntout(nosetsavr),STAT=errstat)
CALL error_alloc('ntout',1,(/nosetsavr/),kndint)
ntout = 0

!
!4.2 Initialise
!--------------
!

iset_420: DO iset=1,nosetsavr
   navrtot(iset) = avrgpars(iset)%tlims(3)/ic3d + 1
ENDDO iset_420

!
!5. Abort if necessary
!----------------------
!

CALL error_abort('time_averages_init',ierrno_inival)

!
!6. Write file headers and initial arrays
!----------------------------------------
!

CALL inout_atts_out('TA')

CALL log_timer_out()


RETURN

END SUBROUTINE time_averages_init

!========================================================================

SUBROUTINE time_averages_grid
!************************************************************************
!
! *time_averages_grid* Write output grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Time_Averages.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - time_averages
!
! External calls -
!
! Module calls - add_secs_to_date, date_to_year, diff_dates, error_alloc,
!                global_mask, loop_index, write_time, Zcoord_var
!
!************************************************************************
!
!1. Log info
!-----------
!

flag = .FALSE.
iset_110: DO iset=1,nosetsavr
   IF (loop_index(avrgpars(iset)%tlims,nt)) THEN
      flag = .TRUE.
   ENDIF
ENDDO iset_110
IF (flag) THEN
   procname(pglev+1) = 'time_averages_grid'
   CALL log_timer_in()
ELSE
   RETURN
ENDIF

!
!2. First output time
!--------------------
!

iset_200: DO iset=1,nosetsavr

   nodim = avrgpars(iset)%nodim
   IF (nodim.GT.0.AND.nt.EQ.avrgpars(iset)%tlims(1)) THEN

      gridded = avrgpars(iset)%gridded
      packing = avrgpars(iset)%packing
      ncout = avrgpars(iset)%ncout
      nrout = avrgpars(iset)%nrout
      nzout = avrgpars(iset)%nzout
      nostats = avrgpars(iset)%nostats
      nowetout = avrgpars(iset)%nowetout
      vcoord = avrgpars(iset)%vcoord
      fill_value = MERGE(double_fill,0.0_kndrlong,packing)
      xlimsglb = avrgpars(iset)%xlims
      ylimsglb = avrgpars(iset)%ylims
      xlims = limloc(:,1,iset)
      ylims = limloc(:,2,iset)
      zlims = avrgpars(iset)%zlims
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
            IF (vcoord.EQ.1) THEN
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
            IF (vcoord.EQ.1) THEN
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
                     IF (vcoord.EQ.1) THEN
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
            i = avrstatlocs(l)%ipos; j = avrstatlocs(l)%jpos
            xout(l,1) = gxcoordglbatc(i,j); yout(l,1) = gycoordglbatc(i,j)
         ENDDO l_2231
         
!        ---bathymetry and vertical grid 
         l_2232: DO l=1,nostatsloc(iset)
            i = istatpos(l,iset); j = jstatpos(l,iset)
            depout(l,1) = depmeanatc(i,j)
            IF (nodim.EQ.3) THEN
               k_22321: DO k=zlims(1),zlims(2),zlims(3)
                  kk = (k-zlims(1))/zlims(3) + 1
                  IF (vcoord.EQ.1) THEN
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
               filepars = avr2d(iset)
            CASE (3)
               filepars = avr3d(iset)
         END SELECT
         IF (.NOT.filepars%defined) CYCLE l_230
         ivar = 1
            
!
!2.3.1 Gridded
!--------------            
!

         IF (gridded) THEN
!           ---horizontal grid
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
!           ---index vector of sea points               
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
               IF( iopt_grid_vtype.EQ.3) THEN
                  CALL combine_write_submod(gsout,filepars,ivar+1,&
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
               IF (vcoord.EQ.1) THEN
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
         IF (vcoord.EQ.1) THEN
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
!3. Time data
!------------
!

iset_300: DO iset=1,nosetsavr

   IF (loop_index(avrgpars(iset)%tlims,nt).AND.&
     & nt.GT.avrgpars(iset)%tlims(1)) THEN
      time_format = avrgpars(iset)%time_format
      refdate = avrgpars(iset)%refdate
      CALL add_secs_to_date(CDateTime,ctime,-avrgpars(iset)%tlims(3),&
                          & 0.5*delt2d)

!
!3.1 Time coordinate
!-------------------
!

      IF (time_format.EQ.8) THEN
         CALL date_to_year(ctime,rtime)
      ELSEIF (time_format.GT.0) THEN
         CALL diff_dates(refdate,ctime,time_format,rtime=rtime)
      ENDIF

!
!3.2 Write
!---------
!

      l_320: DO l=0,3
         SELECT CASE (l)
            CASE (0)
               filepars = avr0d(iset)
            CASE (2)
               filepars = avr2d(iset)
            CASE (3)
               filepars = avr3d(iset)
         END SELECT
         IF (.NOT.filepars%defined.OR.l.EQ.1) CYCLE l_320
         IF (time_format.EQ.0) THEN
            CALL write_time(ctime,filepars)
         ELSE
            CALL write_time(rtime,filepars)
         ENDIF
!        ---store changes in file structure
         SELECT CASE (l)
            CASE (0)
               avr0d(iset) = filepars
            CASE (2)
               avr2d(iset) = filepars
            CASE (3)
               avr3d(iset) = filepars
         END SELECT
      ENDDO l_320

   ENDIF

ENDDO iset_300

CALL log_timer_out()


RETURN

END SUBROUTINE time_averages_grid

!========================================================================

SUBROUTINE time_averages_reset
!************************************************************************
!
! *time_averages_reset* Initialise buffers with time-averaged data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Time_Averages.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - time_averages
!
! Module calls - loop_index 
!
!************************************************************************
!


procname(pglev+1) = 'time_averages_reset'
CALL log_timer_in()

lbuf0d = 0; lbuf2d = 0; lbuf3d = 0

iset_1000: DO iset=1,nosetsavr
   tlims = avrgpars(iset)%tlims
   l0dat = ndat0d(iset); l2dat = ndat2d(iset); l3dat = ndat3d(iset)
   IF ((first.AND.(nt.EQ.tlims(1))).OR.&
     & (second.AND.loop_index(tlims,nt).AND.&
      & nt.GT.tlims(1).AND.nt.LT.tlims(2))) THEN
      IF (l0dat.GT.0) buff0d(lbuf0d+1:lbuf0d+l0dat) = 0.0
      IF (l2dat.GT.0) buff2d(lbuf2d+1:lbuf2d+l2dat) = 0.0
      IF (l3dat.GT.0) buff3d(lbuf3d+1:lbuf3d+l3dat) = 0.0
   ENDIF
   lbuf0d = lbuf0d + l0dat
   lbuf2d = lbuf2d + l2dat
   lbuf3d = lbuf3d + l3dat
ENDDO iset_1000

CALL log_timer_out()
  

RETURN

END SUBROUTINE time_averages_reset

!========================================================================

SUBROUTINE time_averages_update
!************************************************************************
!
! *time_averages_update* Update buffers with time-averaged data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Time_Averages.f90  V2.1.2
!
! Description -
!
! Reference -
!
! Calling program - time_averages
!
! External calls - usrdef_avr0d_vals, usrdef_avr2d_vals, usrdef_avr3d_vals
!
! Module calls - define_out0d_vals, define_out2d_vals, define_out3d_vals,
!                loop_index
!
!************************************************************************
!


procname(pglev+1) = 'time_averages_update'
CALL log_timer_in()

lbuf0d = 0; lbuf2d = 0; lbuf3d = 0

iset_1000: DO iset=1,nosetsavr

   tlims = avrgpars(iset)%tlims
   nodim = avrgpars(iset)%nodim

   IF ((first.AND.nt.GE.tlims(1).AND.nt.LE.tlims(2)).OR.&
     & (second.AND.loop_index(tlims,nt).AND.&
      & nt.GT.tlims(1).AND.nt.LT.tlims(2))) THEN

!
!1. Initialise
!-------------
!
!     ---resolution
      xlims = limloc(:,1,iset)
      ylims = limloc(:,2,iset)
      zlims = avrgpars(iset)%zlims

!     ---number of variables
      noutvars0d = avr0d(iset)%novars
      noutvars2d = avr2d(iset)%novars
      noutvars3d = avr3d(iset)%novars

!
!2. 0-D data
!-----------
!

      IF (avr0d(iset)%defined) THEN
         outvars(1:noutvars0d) = avrvars(ivarsavr0d(iset,1:noutvars0d)) 
         CALL define_out0d_vals(outdat(1:noutvars0d),noutvars0d,&
                              & outvars=outvars(1:noutvars0d))
         WHERE (outvars(1:noutvars0d)%ivarid.GT.0)
            buff0d(lbuf0d+1:lbuf0d+noutvars0d) = outdat/navrtot(iset)
         END WHERE
         IF (ANY(outvars(1:noutvars0d)%ivarid.EQ.0)) THEN
            CALL usrdef_avr0d_vals(outdat(1:novars0d),novars0d)
            WHERE (outvars(1:noutvars0d)%ivarid.EQ.0)
               buff0d(lbuf0d+1:lbuf0d+noutvars0d) = &
                          & outdat(ivarsavr0d(iset,1:noutvars0d))/navrtot(iset)
            END WHERE
         ENDIF
      ENDIF

!
!3. 2-D data
!-----------
!

      IF (nodim.GT.0.AND.avr2d(iset)%defined) THEN
         outvars(1:noutvars2d) = &
                             & avrvars(novars0d+ivarsavr2d(iset,1:noutvars2d))
         flag = ANY(outvars(1:noutvars2d)%ivarid.EQ.0)

!        ---regular data
         IF (avrgpars(iset)%gridded) THEN
            i_302: DO i=xlims(1),xlims(2),xlims(3)
            j_302: DO j=ylims(1),ylims(2),ylims(3)
               ii = (i-xlims(1))/xlims(3) + 1
               jj = (j-ylims(1))/ylims(3) + 1
               CALL define_out2d_vals(outsub(1:noutvars2d),i,j,noutvars2d,&
                                    & outvars=outvars(1:noutvars2d))
               IF (flag) THEN
                  CALL usrdef_tsr2d_vals(outdat(1:novars2d),i,j,novars2d)
                  WHERE (outvars(1:noutvars2d)%ivarid.EQ.0)
                     outsub(1:noutvars2d) = &
                                        & outdat(ivarsavr2d(iset,1:noutvars2d))
                  END WHERE
               ENDIF
               buff2d(lbuf2d+1:lbuf2d+noutvars2d) = &
                                        & buff2d(lbuf2d+1:lbuf2d+noutvars2d) + &
                                        & outsub(1:noutvars2d)/navrtot(iset)
               lbuf2d = lbuf2d + noutvars2d
            ENDDO j_302
            ENDDO i_302

!        ---irregular data
         ELSE
            l_303: DO l=1,nostatsloc(iset)
               i = istatpos(l,iset); j = jstatpos(l,iset)
               CALL define_out2d_vals(outsub(1:noutvars2d),i,j,noutvars2d,&
                                    & outvars=outvars(1:noutvars2d))
               IF (flag) THEN
                  CALL usrdef_tsr2d_vals(outdat(1:novars2d),i,j,novars2d)
                  WHERE (outvars(1:noutvars2d)%ivarid.EQ.0)
                     outsub(1:noutvars2d) = &
                                         & outdat(ivarsavr2d(iset,1:noutvars2d))
                  END WHERE
               ENDIF
               buff2d(lbuf2d+1:lbuf2d+noutvars2d) = &
                                        & buff2d(lbuf2d+1:lbuf2d+noutvars2d) + &
                                        & outsub(1:noutvars2d)/navrtot(iset)
               lbuf2d = lbuf2d + noutvars2d
            ENDDO l_303
         ENDIF

      ENDIF

!
!4. 3-D data
!-----------
!

      IF (nodim.EQ.3.AND.avr3d(iset)%defined) THEN
         outvars(1:noutvars3d) = &
                      & avrvars(novars0d+novars2d+ivarsavr3d(iset,1:noutvars3d))
         flag = ANY(outvars(1:noutvars3d)%ivarid.EQ.0)

!        ---regular data
         IF (avrgpars(iset)%gridded) THEN
            i_401: DO i=xlims(1),xlims(2),xlims(3)
            j_401: DO j=ylims(1),ylims(2),ylims(3)
               ii = (i-xlims(1))/xlims(3) + 1
               jj = (j-ylims(1))/ylims(3) + 1
               k_4011: DO k=zlims(1),zlims(2),zlims(3)
                  kk = (k-zlims(1))/zlims(3) + 1
                  CALL define_out3d_vals(outsub(1:noutvars3d),i,j,k,noutvars3d,&
                                       & outvars=outvars(1:noutvars3d))
                  IF (flag) THEN
                     CALL usrdef_tsr3d_vals(outdat(1:novars3d),i,j,k,novars3d)
                     WHERE (outvars(1:noutvars3d)%ivarid.EQ.0)
                        outsub(1:noutvars3d) = &
                                         & outdat(ivarsavr3d(iset,1:noutvars3d))
                     END WHERE
                  ENDIF
                  buff3d(lbuf3d+1:lbuf3d+noutvars3d) = &
                                        & buff3d(lbuf3d+1:lbuf3d+noutvars3d) + &
                                        & outsub(1:noutvars3d)/navrtot(iset)
                  lbuf3d = lbuf3d + noutvars3d
               ENDDO k_4011
            ENDDO j_401
            ENDDO i_401

!        ---irregular data
         ELSE
            l_402: DO l=1,nostatsloc(iset)
               i = istatpos(l,iset); j = jstatpos(l,iset)
               k_4021: DO k=zlims(1),zlims(2),zlims(3)
                  kk = (k-zlims(1))/zlims(3) + 1
                  CALL define_out3d_vals(outsub(1:noutvars3d),i,j,k,noutvars3d,&
                                       & outvars=outvars(1:noutvars3d))
                  IF (flag) THEN
                     CALL usrdef_tsr3d_vals(outdat(1:novars3d),i,j,k,novars3d)
                     WHERE (outvars(1:noutvars3d)%ivarid.EQ.0)
                        outsub(1:noutvars3d) = &
                                         & outdat(ivarsavr3d(iset,1:noutvars3d))
                     END WHERE
                  ENDIF
                  buff3d(lbuf3d+1:lbuf3d+noutvars3d) = &
                                        & buff3d(lbuf3d+1:lbuf3d+noutvars3d) + &
                                        & outsub(1:noutvars3d)/navrtot(iset)
                  lbuf3d = lbuf3d + noutvars3d
               ENDDO k_4021
            ENDDO l_402
         ENDIF

      ENDIF

   ELSE

      lbuf0d = lbuf0d + ndat0d(iset)
      lbuf2d = lbuf2d + ndat2d(iset)
      lbuf3d = lbuf3d + ndat3d(iset)

   ENDIF

ENDDO iset_1000

CALL log_timer_out()


RETURN

END SUBROUTINE time_averages_update

!========================================================================

SUBROUTINE time_averages_write
!************************************************************************
!
! *time_averages_write* Write time-averaged output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Time_Averages.f90  V2.10.1
!
! Description -
!
! Reference -
!
! Calling program - time_averages
!
! Module calls - combine_write_stats_loc, combine_write_submod, error_alloc,
!                loop_index, mask_array, write_vars
!
!************************************************************************
!
!1. Log info
!-----------
!

flag = .FALSE.
iset_110: DO iset=1,nosetsavr
   tlims = avrgpars(iset)%tlims
   IF (loop_index(tlims,nt).AND.(nt.GT.tlims(1))) THEN
      flag = .TRUE.
   ENDIF
ENDDO iset_110
IF (flag) THEN
   procname(pglev+1) = 'time_averages_write'
   CALL log_timer_in()
ELSE
   RETURN
ENDIF

!
!2. Time-averaged output
!-----------------------
!

lbuf0d = 0; lbuf2d = 0; lbuf3d = 0

iset_2000: DO iset=1,nosetsavr
   tlims = avrgpars(iset)%tlims
   IF (loop_index(tlims,nt).AND.(nt.GT.tlims(1))) THEN

!
!2.1 Initialise
!--------------
!

      gridded = avrgpars(iset)%gridded
      ncout = avrgpars(iset)%ncout
      nrout = avrgpars(iset)%nrout
      nzout = avrgpars(iset)%nzout
      nowetout = avrgpars(iset)%nowetout
      nostats = avrgpars(iset)%nostats
      lprocs = limprocs(:,:,:,iset)
      ntout(iset) = ntout(iset) + 1
      fill_value = double_fill

      noutvars0d = avr0d(iset)%novars
      noutvars2d = avr2d(iset)%novars
      noutvars3d = avr3d(iset)%novars

      nocoords0d = avr0d(iset)%nocoords
      nocoords2d = avr2d(iset)%nocoords
      nocoords3d = avr3d(iset)%nocoords

!
!2.2 Create mask array
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
!2.3 0-D output
!--------------
!

      IF (master.AND.avr0d(iset)%defined) THEN
         outvars(1:noutvars0d) = avrvars(ivarsavr0d(iset,1:noutvars0d)) 
         vecids(1:noutvars0d) = (/(nocoords0d+ivar,ivar=1,noutvars0d)/)
         CALL write_vars(buff0d(lbuf0d+1:lbuf0d+noutvars0d),avr0d(iset),0,&
                       & varatts=outvars(1:noutvars0d),&
                       & vecids=vecids(1:noutvars0d))
      ENDIF

!
!2.4 2-D output
!--------------
!

      IF (avr2d(iset)%defined) THEN

!
!2.4.1 Create output array
!-------------------------
!

         IF (gridded) THEN
            ALLOCATE (out3dsub(nclocout(iset),nrlocout(iset),noutvars2d),&
                    & STAT=errstat)
            CALL error_alloc('out3dsub',3,(/nclocout(iset),nrlocout(iset),&
                                          & noutvars2d/),kndrtype)
            nodat = ndat2d(iset)
            out3dsub = RESHAPE(buff2d(lbuf2d+1:lbuf2d+nodat),&
                         & SHAPE=(/nclocout(iset),nrlocout(iset),noutvars2d/),&
                         & ORDER=(/3,2,1/))
         ELSE
            ALLOCATE (out2dsub(nostatsloc(iset),noutvars2d),STAT=errstat)
            CALL error_alloc('out2dsub',2,(/nostatsloc(iset),noutvars2d/),&
                            & kndrtype)
            nodat = ndat2d(iset)
            out2dsub = RESHAPE(buff2d(lbuf2d+1:lbuf2d+nodat),&
                         & SHAPE=(/nostatsloc(iset),noutvars2d/),ORDER=(/2,1/))
         ENDIF

!
!2.4.2 Apply land mask
!---------------------
!

         IF (gridded.AND.fmask.AND.avr2d(iset)%fill) THEN
            ivar_242: DO ivar=1,noutvars2d
               WHERE (.NOT.maskvals)
                  out3dsub(:,:,ivar) = fill_value
               END WHERE
            ENDDO ivar_242
         ENDIF

!
!2.4.3 Write data
!----------------
!

         outvars(1:noutvars2d) = &
                             & avrvars(novars0d+ivarsavr2d(iset,1:noutvars2d))
         vecids(1:noutvars2d) = (/(nocoords2d+ivar,ivar=1,noutvars2d)/)
         IF (gridded) THEN
            CALL combine_write_submod(out3dsub,avr2d(iset),0,&
                                   & (/ncout,nrout,noutvars2d/),lprocs,&
                                   & nowet=nowetout,&
                                   & varatts=outvars(1:noutvars2d),&
                                   & vecids=vecids(1:noutvars2d))
         ELSE
            CALL combine_write_stats_loc(out2dsub,avr2d(iset),0,maxstats,&
                                   & nostats,nostatsprocs(:,iset),&
                                   & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                   & varatts=outvars(1:noutvars2d),&
                                   & vecids=vecids(1:noutvars2d))
         ENDIF

!
!2.4.4 Deallocate data array
!---------------------------
!

         IF (gridded) THEN
            DEALLOCATE (out3dsub)
         ELSE
            DEALLOCATE (out2dsub)
         ENDIF

!
!2.4.5 Update position in buffer
!-------------------------------
!

         lbuf2d = lbuf2d + nodat

      ENDIF

!
!2.5 3-D output
!--------------
!

      IF (avr3d(iset)%defined) THEN

!
!2.5.1 Create output array
!-------------------------
!

         IF (gridded) THEN
            ALLOCATE (out4dsub(nclocout(iset),nrlocout(iset),&
                             & nzout,noutvars3d),STAT=errstat)
            CALL error_alloc('out4dsub',4,(/nclocout(iset),nrlocout(iset),&
                                          & nzout,noutvars3d/),kndrtype)
            nodat = ndat3d(iset)
            out4dsub = RESHAPE(buff3d(lbuf3d+1:lbuf3d+nodat),&
                & SHAPE=(/nclocout(iset),nrlocout(iset),nzout,noutvars3d/),&
                & ORDER=(/4,3,2,1/))
         ELSE
            ALLOCATE (out3dsub(nostatsloc(iset),nzout,noutvars3d),STAT=errstat)
            CALL error_alloc('out3dsub',3,(/nostatsloc(iset),nzout,&
                                          & noutvars3d/),kndrtype)
            nodat = ndat3d(iset)
            out3dsub = RESHAPE(buff3d(lbuf3d+1:lbuf3d+nodat),&
                & SHAPE=(/nostatsloc(iset),nzout,noutvars3d/),ORDER=(/3,2,1/))
         ENDIF

!
!2.5.2 Apply land mask
!---------------------
!

         IF (gridded.AND.fmask.AND.avr3d(iset)%fill) THEN
            ivar_252: DO ivar=1,noutvars3d
            k_252: DO k=1,nzout
               WHERE (.NOT.maskvals)
                  out4dsub(:,:,k,ivar) = fill_value
               END WHERE
            ENDDO k_252
            ENDDO ivar_252
         ENDIF

!
!2.5.3 Write data
!----------------
!

         outvars(1:noutvars3d) = &
                      & avrvars(novars0d+novars2d+ivarsavr3d(iset,1:noutvars3d))
         vecids(1:noutvars3d) = (/(nocoords3d+ivar,ivar=1,noutvars3d)/)
         IF (gridded) THEN
            CALL combine_write_submod(out4dsub,avr3d(iset),0,&
                                   & (/ncout,nrout,nzout,noutvars3d/),lprocs,&
                                   & nowet=nowetout,&
                                   & varatts=outvars(1:noutvars3d),&
                                   & vecids=vecids(1:noutvars3d))
         ELSE
            CALL combine_write_stats_loc(out3dsub,avr3d(iset),0,maxstats,&
                                   & nostats,nostatsprocs(:,iset),&
                                   & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                   & varatts=outvars(1:noutvars3d),&
                                   & vecids=vecids(1:noutvars3d))
         ENDIF

!
!2.5.4 Deallocate data array
!---------------------------
!

         IF (gridded) THEN
            DEALLOCATE (out4dsub)
         ELSE
            DEALLOCATE (out3dsub)
         ENDIF

!
!2.5.5 Update position in buffer
!-------------------------------
!

         lbuf3d = lbuf3d + nodat

      ENDIF

!
!2.6 Deallocate mask arrays
!--------------------------
!

      IF (ALLOCATED(maskvals)) DEALLOCATE(maskvals)

   ELSE

!     ---update position in buffer in absence of output
      lbuf0d = lbuf0d + ndat0d(iset)
      lbuf2d = lbuf2d + ndat2d(iset)
      lbuf3d = lbuf3d + ndat3d(iset)

   ENDIF

ENDDO iset_2000

CALL log_timer_out()


RETURN

END SUBROUTINE time_averages_write

END SUBROUTINE time_averages
