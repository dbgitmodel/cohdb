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

!************************************************************************
!
! *Model_Finalisation* Series of routines for model finalisation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Finalisation.f90  V2.11
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - 
!
! Reference -
!
! Routines - coherens_end, simulation_end, timer_report, write_phsics
!
!************************************************************************
!

!========================================================================

SUBROUTINE coherens_end
!************************************************************************
!
! *coherens_end*  Finalize program
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Finalisation.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - comms_barrier, comms_finalize, error_abort
!
!************************************************************************
!
USE iopars
USE paralpars
USE comms_MPI, ONLY: comms_barrier, comms_finalize
USE error_routines, ONLY: error_abort

IMPLICIT NONE


IF (parallel_set) THEN
   CALL comms_barrier(comm_world_MPI)
   CALL comms_finalize
ENDIF
CALL error_abort('coherens_end',ierrno_MPI)

DEALLOCATE(idprocs,idprocsglb)
  

RETURN

END SUBROUTINE coherens_end

!========================================================================

SUBROUTINE simulation_end
!************************************************************************
!
! *simulation_end* Finalise current simulation before starting next one
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Finalisation.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - coh2wav_end, comms_barrier, deallocate_bio_arrays,
!                  deallocate_mg_arrays, deallocate_mod_arrays,
!                  deallocate_morph_arrays, deallocate_part_arrays,
!                  deallocate_sed_arrays, inquire_unit, timer_report
!
! Module calls - close_file, get_unit, rng_finalize
!
!************************************************************************
!
USE grid
USE iopars
USE paralpars
USE switches
USE syspars
USE timepars
USE wavevars
USE comms_MPI, ONLY: comms_barrier
USE inout_routines, ONLY: close_file, close_filepars, get_unit, inquire_unit
USE rng_library, ONLY: rng_finalize
USE time_routines, ONLY: log_timer_in

IMPLICIT NONE
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: form
INTEGER :: idesc, ifil, iostat, iotype, iunit


procname(pglev+1) = 'simulation_end'
CALL log_timer_in()

!
!1. Close opened files
!---------------------
!

idesc_110: DO idesc=1,MaxIOTypes
iotype_110: DO iotype=1,2
ifil_110: DO ifil=1,maxdatafiles(idesc,iotype)
   iunit = modfiles(idesc,ifil,iotype)%iunit
   iostat = modfiles(idesc,ifil,iotype)%iostat
   IF (iostat.GT.0.AND.iunit.NE.int_fill) THEN
      form = modfiles(idesc,ifil,iotype)%form
      flag = inquire_unit(iunit,form)
      IF (flag) CALL close_filepars(modfiles(idesc,ifil,iotype))
   ENDIF
ENDDO ifil_110
ENDDO iotype_110
ENDDO idesc_110

!
!2. Close random generators
!--------------------------
!

CALL rng_finalize

!
!3. Write timer report
!---------------------
!

IF (timer) CALL timer_report

!
!4. Deallocate arrays
!--------------------
!

IF (modelid.EQ.modelidcoh) THEN
   IF (iopt_part_model.LT.3) THEN
      IF (iopt_hydro_impl.EQ.1) CALL deallocate_mg_arrays
      CALL deallocate_mod_arrays
      IF (iopt_sed.GT.0) CALL deallocate_sed_arrays
      IF (iopt_biolgy.GT.0) CALL deallocate_bio_arrays
      IF (iopt_morph.GT.0) CALL deallocate_morph_arrays
      IF (iopt_dar.GT.0) CALL deallocate_dar_arrays
   ENDIF
   IF (iopt_waves_couple.GT.0.AND.(.NOT.cold_start)) CALL coh2wav_end
ELSEIF (.NOT.cold_start.AND.modelid.EQ.modelidwav) THEN
   DEALLOCATE (gxcoordglb,gycoordglb,seapointglb)
   DEALLOCATE (gxcoordglbwav,gycoordglbwav,maskglbwav)
   DEALLOCATE (comprocs)
ENDIF
IF (iopt_part_model.GT.0) CALL deallocate_part_arrays

!
!5. Write exit code
!------------------
!

IF (loglev2.GT.0) THEN
   IF (.NOT.cold_start.AND.modelid.EQ.modelidpart) THEN
      WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Particle model run terminated'
   ENDIF
      IF (.NOT.cold_start.AND.modelid.EQ.modelidwav) THEN
      WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Wave model run terminated'
   ENDIF
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, &
                          & 'End of simulation: '//TRIM(runtitle)
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, &
                          & 'End of simulation: '//TRIM(runtitle)
   ENDIF
ENDIF

!
!6. Close and reset log and error files, reset parameters
!--------------------------------------------------------
!
!---warning
IF (warnflag) CALL close_file(iowarn,'A',filename=warlog_file)
warning = .TRUE.; warnflag = .FALSE.; iowarn = 0

!---monitoring
IF (monflag) CALL close_file(iomon,'A',filename=monlog_file)
monlog = 0; monflag = .FALSE.; iomon = 0 
IF (sedflag) CALL close_file(iosed,'A',filename=sedlog_file)
sedlog = .TRUE.; sedflag = .FALSE.; iosed = 0  

!---error
INQUIRE (FILE=errlog_file,OPENED=flag)
IF (flag) CALL close_file(ioerr,'A',filename=errlog_file,fildel=.TRUE.)
errchk = master
ioerr = 0; maxerrors = MaxErrMesgs
DEALLOCATE (levprocs_err)

!---inilog
INQUIRE (FILE=inilog_file,OPENED=flag)
IF (flag) CALL close_file(iolog,'A',filename=inilog_file)
DEALLOCATE (levprocs_ini)

!---runlog
INQUIRE (FILE=runlog_file,OPENED=flag)
IF (flag) CALL close_file(iolog,'A',filename=runlog_file)
exitlog = .TRUE.; iolog = 0; loglev1 = 0; loglev2 = 0; nerrs = 0
DEALLOCATE (levprocs_run)

!---time
nt = 0

!
!7. Synchronise processors
!-------------------------
!

IF (parallel_set) CALL comms_barrier(comm_world_MPI)

pglev = pglev - 1


RETURN

END SUBROUTINE simulation_end

!========================================================================

SUBROUTINE timer_report
!************************************************************************
!
! *timer_report* Write timing report file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Finalisation.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - simulation_end
!
! External calls -
!
! Module calls - clock_date_time, close_file, collect_vars, comms_recv_int,
!                comms_send_int, diff_dates, error_alloc, open_file
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE syspars
USE timepars
USE comms_MPI, ONLY: comms_recv_real, comms_send_real
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_file, open_file
USE paral_comms, ONLY: collect_vars
USE time_routines, ONLY: clock_date_time, diff_dates, log_timer_in, &
                       & log_timer_out

!
!* Local variables
!
CHARACTER (LEN=lentime) :: clocktime1, clocktime2
CHARACTER (LEN=36) :: exectime
INTEGER :: iotime, iproc, it, n, nodat, np, nomsecs
INTEGER (KIND=kndilong) :: nohours, nomins, nosecs
REAL (KIND=kndrlong) :: exsecs
REAL :: ptmaster, ptmax, ptmean, ptmin
CHARACTER (LEN=22), DIMENSION(5) :: ctime
INTEGER (KIND=kndilong), DIMENSION(5) :: itime
REAL, ALLOCATABLE, DIMENSION(:) :: ptpcc
REAL, ALLOCATABLE, DIMENSION(:,:) :: ptpccprocs


procname(pglev+1) = 'timer_report'
CALL log_timer_in()


!
!1. Open file
!------------
!

timer = .FALSE.
IF (master) CALL open_file(iotime,timing_file,'OUT','A')

!
!2. Evaluate and write execution time
!------------------------------------
!
!2.1 Execution time from real-time clock
!---------------------------------------
!

clocktime1 = ClockTime
CALL clock_date_time(clocktime2)
CALL diff_dates(clocktime1,clocktime2,0,nosecs=nosecs,millisecs=nomsecs)

!
!2.2 Write execution time in appropriate format
!----------------------------------------------
!
!---milliseconds
itime(1) = nomsecs
WRITE (ctime(1),'(I3.3)') itime(1); ctime(1) = ADJUSTL(ctime(1))
!---seconds
If (timer_format.EQ.1) THEN
   itime(2) = nosecs
   WRITE (ctime(2),'(I22)') itime(2); ctime(2) = ADJUSTL(ctime(2))
ELSE
   itime(2) = MOD(nosecs,60_kndilong)
   WRITE (ctime(2),'(I2.2)') itime(2); ctime(2) = ADJUSTL(ctime(2))
   nomins = nosecs/60
ENDIF
!---minutes
IF (timer_format.EQ.2) THEN
   itime(3) = nomins
   WRITE (ctime(3),'(I22)') itime(3); ctime(3) = ADJUSTL(ctime(3))
ELSEIF (timer_format.GT.2) THEN
   itime(3) = MOD(nomins,60_kndilong)
   WRITE (ctime(3),'(I2.2)') itime(3); ctime(3) = ADJUSTL(ctime(3))
   nohours = nomins/60
ENDIF
!---hours and days
IF (timer_format.EQ.3) THEN
   itime(4) = nohours
   WRITE (ctime(4),'(I22)') itime(4); ctime(4) = ADJUSTL(ctime(4))
ELSEIF (timer_format.EQ.4) THEN
   itime(4) = MOD(nohours,24_kndilong)
   WRITE (ctime(4),'(I2.2)') itime(4); ctime(4) = ADJUSTL(ctime(4))
   itime(5) = nohours/24
   WRITE (ctime(5),'(I22)') itime(5); ctime(5) = ADJUSTL(ctime(5))
ENDIF

!---format
SELECT CASE(timer_format)
CASE(1); exectime = TRIM(ctime(2))//'s.'//TRIM(ctime(1))
CASE(2); exectime = TRIM(ctime(3))//'m'//TRIM(ctime(2))//'s.'//TRIM(ctime(1))
CASE(3); exectime = TRIM(ctime(4))//'h'//TRIM(ctime(3))//'m'//&
                  & TRIM(ctime(2))//'s.'//TRIM(ctime(1))
CASE(4); exectime = TRIM(ctime(5))//'d'//TRIM(ctime(4))//'h'//&
                  & TRIM(ctime(3))//'m'//TRIM(ctime(2))//'s.'//TRIM(ctime(1))
END SELECT

!---write to time report
IF (master) THEN
   WRITE (iotime,'(A)') TRIM(runtitle)//': '//TRIM(exectime)
ENDIF
IF (levtimer.EQ.1) GOTO 1000

!
!3. Time report for model components
!-----------------------------------
!
!3.1 Reorder terms
!-----------------
!

SELECT CASE (iopt_grid_nodim)
CASE (1)
   nopcc(itm_1dmode) = nopcc(itm_1dmode) + nopcc(itm_3dmode)
CASE (2)
   nopcc(itm_2dmode) = nopcc(itm_2dmode) + nopcc(itm_astro)
CASE (3)
   nopcc(itm_3dmode) = nopcc(itm_3dmode) + nopcc(itm_astro)
END SELECT

IF (iopt_dens.GT.0) THEN
   nopcc(itm_dens) = nopcc(itm_dens) + nopcc(itm_sal) + nopcc(itm_temp) + &
                   & nopcc(itm_phgrad)
ENDIF

IF (iopt_MPI.EQ.0) THEN
   nopcc(itm_com_coll) = 0; nopcc(itm_com_comb) = 0
   nopcc(itm_com_copy) = 0; nopcc(itm_com_dist) = 0
   nopcc(itm_com_exch) = 0; nopcc(itm_com_util) = 0
   nopcc(itm_coms) = 0
ELSEIF (iopt_MPI.EQ.1) THEN
   nopcc(itm_coms) = nopcc(itm_com_coll) + nopcc(itm_com_comb) + &
                   & nopcc(itm_com_copy) + nopcc(itm_com_dist) + &
                   & nopcc(itm_com_exch) + nopcc(itm_com_util)
ENDIF

nopcc(itm_inout) = nopcc(itm_input) + nopcc(itm_output)

IF (modelid.EQ.modelidwav) THEN
   it_310: DO it=1,MaxTimers
      IF (it.NE.itm_MCT.AND.it.NE.itm_wavemod) nopcc(it) = 0
   ENDDO it_310
ENDIF

!
!3.2 Collect
!-----------
!
!---allocate
ALLOCATE (ptpcc(MaxTimers),STAT=errstat)
CALL error_alloc('ptpcc',1,(/MaxTimers/),kndrtype)
ptpcc = 0.0
ALLOCATE (ptpccprocs(MaxTimers,npworld),STAT=errstat)
CALL error_alloc('ptpccprocs',2,(/MaxTimers,npworld/),kndrtype)
ptpccprocs = 0.0

!--- % of execution time
exsecs = nosecs + 0.001*nomsecs
ptpcc = ((nopcc/npcc_rate+MOD(nopcc,INT(npcc_rate,KIND=kndilong))/&
             & REAL(npcc_rate))/exsecs)*100.0
ptpcc = NINT(1000.0*ptpcc)/1000.0

!---collect from particle model
np = nprocscoh + 1
IF (modelid.EQ.modelidpart) THEN
   ptpccprocs(:,np) = ptpcc
   CALL collect_vars(ptpcc,ptpccprocs(:,np:np),1,0,&
                   & comm=comm_model,commtype=5)
   IF (iopt_part_model.EQ.2.AND.mastermod) THEN
      nodat = MaxTimers
      CALL comms_send_real(ptpccprocs(:,np),nodat,idmaster,1001,comm_world_MPI,&
                         & 1,0)
   ENDIF
ELSEIF (modelid.EQ.modelidcoh) THEN
   CALL collect_vars(ptpcc,ptpccprocs(:,1:nprocs),nprocs,0)
   IF (iopt_part_model.EQ.2.AND.master) THEN
      nodat = MaxTimers
      CALL comms_recv_real(ptpccprocs(:,np),nodat,idmasterpart,1001,&
                         & comm_world_MPI,0)
   ENDIF 
ENDIF

!---collect from wave model
np = nprocscoh+nprocspart+1
IF (modelid.EQ.modelidwav) THEN 
   CALL collect_vars(ptpcc,ptpccprocs(:,np:npworld),nprocswav,0,&
                   & comm=comm_model,commtype=5)
   IF (iopt_waves_model.GT.0.AND.mastermod) THEN
      nodat = MaxTimers*nprocswav
      CALL comms_send_real(ptpccprocs(:,np:npworld),nodat,idmaster,1002,&
                         & comm_world_MPI,1,0)
   ENDIF
ELSEIF (modelid.EQ.modelidcoh) THEN
   CALL collect_vars(ptpcc,ptpccprocs(:,1:nprocs),nprocs,0)
   IF (iopt_waves_model.GT.0.AND.master) THEN
      nodat = MaxTimers*nprocswav
      CALL comms_recv_real(ptpccprocs(:,np:npworld),nodat,idmasterwav,1002,&
                         & comm_world_MPI,0)
   ENDIF 
ENDIF

!
!3.4 Write to time report
!------------------------
!

IF (master) THEN
   it_340: DO it=1,MaxTimers
      IF (ANY(ptpccprocs(it,:).GT.0.0)) THEN
         IF (iopt_MPI.EQ.0) THEN
            WRITE (iotime,9001) desctimer(it), ptpccprocs(it,1)
         ELSE
            n = 0; ptmax = 0.0; ptmin = 100.0; ptmean = 0.0; ptmaster = 0.0
            iproc_341: DO iproc=1,npworld
               IF (idprocsglb(iproc).EQ.idmaster) ptmaster = ptpccprocs(it,iproc)
               IF (ptpccprocs(it,iproc).GT.0.0) THEN
                  n = n + 1
                  ptmean = ptmean + ptpccprocs(it,iproc)
                  ptmin = MIN(ptmin,ptpccprocs(it,iproc))
                  ptmax = MAX(ptmax,ptpccprocs(it,iproc))
               ENDIF
            ENDDO iproc_341
            ptmean = ptmean/REAL(n)
            WRITE (iotime,9001) desctimer(it), ptmax, ptmin, ptmean, ptmaster
            IF (levtimer.EQ.3) WRITE (iotime,9002) ptpccprocs(it,:)
         ENDIF
      ENDIF
   ENDDO it_340
ENDIF

!
!3.5 Deallocate
!--------------
!

DEALLOCATE (ptpcc,ptpccprocs)

!
!4 Close file
!------------
!

IF (master) CALL close_file(iotime,'A')

1000 CALL log_timer_out()


RETURN

9001 FORMAT (A,': ',4F7.3)
9002 FORMAT (2X,10F7.3)

END SUBROUTINE timer_report

!========================================================================

SUBROUTINE write_phsics
!************************************************************************
!
! *write_phsics* Write physical initial conditions in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Model_Finalisation.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - coherens_main
!
! External calls -
!
! Module calls - close_filepars, combine_write_mod, combine_write_stats_glb,
!                error_alloc_struc, open_filepars, set_modfiles_atts,
!                set_modvars_atts, write_atts_mod, write_time, write_vars
!
!************************************************************************
!
USE currents
USE datatypes
USE density
USE depths
USE fluxes
USE grid
USE gridpars
USE iopars
USE obconds
USE paralpars
USE structures
USE switches
USE tide
USE timepars
USE turbulence
USE error_routines, ONLY: error_alloc_struc
USE inout_paral, ONLY: combine_write_mod, combine_write_stats_glb
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_time, write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
INTEGER :: numvars, varid
INTEGER, SAVE :: itrestart
TYPE(FileParams) :: filepars
INTEGER, DIMENSION(3) :: lbounds, ubounds
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


IF (norestarts.EQ.0.OR.modfiles(io_fincon,ics_phys,2)%status.EQ.'0') RETURN
IF (nt.EQ.0) itrestart = 0
IF (itrestart.GT.norestarts) RETURN

procname(pglev+1) = 'write_phsics'
flag = .FALSE.
filepars = modfiles(io_fincon,ics_phys,2)

!
!1. Write file header on first call
!----------------------------------
!

IF (nt.EQ.ntrestart(1)) THEN

   CALL log_timer_in()
   flag = .TRUE.

!  ---file attributes
   CALL set_modfiles_atts(io_fincon,ics_phys,2)
   filepars = modfiles(io_fincon,ics_phys,2)
   numvars = filepars%novars + 1

!  ---open file
   IF (master) CALL open_filepars(filepars)

!  ---variable attributes
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
   CALL set_modvars_atts(io_fincon,ics_phys,2,varatts,numvars)

!  ---write attributes
   IF (master) CALL write_atts_mod(filepars,varatts,numvars)

   itrestart = 1

ENDIF

!
!2. Write physical data
!---------------------
!

IF (itrestart.GT.0.AND.itrestart.LE.norestarts) THEN

   IF (nt.EQ.ntrestart(itrestart)) THEN
      IF (itrestart.GT.1) THEN
         CALL log_timer_in()
         flag = .TRUE.
      ENDIF

!
!2.1 Date/time
!-------------
!

      CALL write_time(CDateTime,filepars)
      varid = 1

!
!2.2 2-D mode
!------------
!

      IF (iopt_mode_2D.GT.0) THEN

!        ---X-current
         varid = varid + 1
         lbounds(1:2) = (/1-nhalo,1-nhalo/)
         ubounds(1:2) = (/nc+nhalo,nr+nhalo/)
         CALL combine_write_mod(udvel,filepars,varid,lbounds(1:2),ubounds(1:2),&
                              & 'U  ',varatts=varatts)

!        ---Y-current
         varid = varid + 1
         CALL combine_write_mod(vdvel,filepars,varid,lbounds(1:2),ubounds(1:2),&
                              & 'V  ',varatts=varatts)

!        ---surface elevation
         varid = varid + 1
         lbounds(1:2) = 0; ubounds(1:2) = (/nc+1,nr+1/)
         CALL combine_write_mod(zeta,filepars,varid,lbounds(1:2),ubounds(1:2),&
                              & 'C  ',varatts=varatts)

      ENDIF

!
!2.3 3-D currents
!----------------
!

      IF (iopt_mode_3D.GT.0) THEN

!        ---X-current
         varid = varid + 1
         lbounds = (/1-nhalo,1-nhalo,1/); ubounds = (/nc+nhalo,nr+nhalo,nz/)
         CALL combine_write_mod(uvel,filepars,varid,lbounds,ubounds,'U  ',&
                              & varatts=varatts)

!        ---Y-current
         varid = varid + 1
         CALL combine_write_mod(vvel,filepars,varid,lbounds,ubounds,'V  ',&
                              & varatts=varatts)

      ENDIF

!     ---vertical current
      IF (iopt_mode_3D.GT.0.AND.iopt_grid_nodim.EQ.3) THEN
         varid = varid + 1
         lbounds = (/0,0,1/); ubounds = (/nc,nr,nz+1/)
         CALL combine_write_mod(wvel,filepars,varid,lbounds,ubounds,'W  ',&
                              & varatts=varatts)
      ENDIF

!
!2.4  Density arrays
!-------------------
!

      lbounds = (/1-nhalo,1-nhalo,1/); ubounds = (/nc+nhalo,nr+nhalo,nz/)
!     ---temperature
      IF (iopt_temp.GT.0) THEN
         varid = varid + 1
         CALL combine_write_mod(temp,filepars,varid,lbounds,ubounds,'C  ',&
                              & varatts=varatts)
      ENDIF

!     ---salinity
      IF (iopt_sal.GT.0) THEN
         varid = varid + 1
         CALL combine_write_mod(sal,filepars,varid,lbounds,ubounds,'C  ',&
                              & varatts=varatts)
      ENDIF

!
!2.5 Turbulence arrays
!---------------------
!

      IF (iopt_vdif_coef.EQ.3) THEN
         ubounds = (/nc+nhalo,nr+nhalo,nz+1/)
         
!        ---turbulence energy
         varid = varid + 1
         CALL combine_write_mod(tke,filepars,varid,lbounds,ubounds,'W  ',&
                              & varatts=varatts)

!        ---mixing length
         IF (iopt_turb_ntrans.EQ.2.AND.iopt_turb_param.EQ.1) THEN
            varid = varid + 1
            CALL combine_write_mod(zlmix,filepars,varid,lbounds,ubounds,'W  ',&
                                 & varatts=varatts)

!        ---dissipation rate
         ELSEIF (iopt_turb_ntrans.EQ.2.AND.iopt_turb_param.EQ.2) THEN
            varid = varid + 1
            CALL combine_write_mod(dissip,filepars,varid,lbounds,ubounds,'W  ',&
                                 & varatts=varatts)
         ENDIF

      ENDIF

!
!2.6 Bottom stress arrays
!------------------------
!

      IF (iopt_bstres_form.EQ.2) THEN
         lbounds(1:2) = 1; ubounds(1:2) = (/nc,nr/)
         IF (iopt_bstres_drag.EQ.2) THEN

!           ---bottom drag coefficient
            varid = varid + 1
            CALL combine_write_mod(bdragcoefatc(1:ncloc,1:nrloc),filepars,&
                                 & varid,lbounds(1:2),ubounds(1:2),'C  ',&
                                 & varatts=varatts)
            
         ELSEIF (iopt_bstres_drag.EQ.4) THEN

!           ---roughness length
            varid = varid + 1
            CALL combine_write_mod(zroughatc(1:ncloc,1:nrloc),filepars,varid,&
                                 & lbounds(1:2),ubounds(1:2),'C  ',&
                                 & varatts=varatts)
         ENDIF
      ENDIF

!
!2.7 Tidal arrays
!----------------
!

      IF (nconobc.GT.0) THEN 
         varid = varid + 1
         CALL write_vars(fnode_obc,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(phase_obc,filepars,varid,varatts=varatts)
      ENDIF
      IF (iopt_astro_tide.EQ.1) THEN
         varid = varid + 1
         CALL write_vars(fnode_astro,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(phase_astro,filepars,varid,varatts=varatts)
      ENDIF

!
!2.8 Energy loss terms at weirs
!------------------------------
!

      IF (iopt_weibar.EQ.1) THEN
         IF (numwbaru.GT.0) THEN
            varid = varid + 1
            CALL combine_write_stats_glb(wbarelossu,filepars,varid,numwbaru,&
                                       & nowbaruprocs,indexwbaruprocs,&
                                       & varatts=varatts)
         ENDIF
         IF (numwbarv.GT.0) THEN
            varid = varid + 1
            CALL combine_write_stats_glb(wbarelossv,filepars,varid,numwbarv,&
                                       & nowbarvprocs,indexwbarvprocs,&
                                       & varatts=varatts)
         ENDIF
      ENDIF

!
!2.9 Open boundary arrays
!------------------------
!

      IF (iopt_obc_sal.EQ.1) THEN
         IF (nobu.GT.0) THEN
            varid = varid + 1
            CALL combine_write_stats_glb(obcsalatu,filepars,varid,nobu,&
                                       & nobuprocs,indexobuprocs,&
                                       & varatts=varatts)
         ENDIF
         IF (nobv.GT.0) THEN
            varid = varid + 1
            CALL combine_write_stats_glb(obcsalatv,filepars,varid,nobv,&
                                       & nobvprocs,indexobvprocs,&
                                       & varatts=varatts)
         ENDIF
      ENDIF
      IF (iopt_obc_temp.EQ.1) THEN
         IF (nobu.GT.0) THEN
            varid = varid + 1
            CALL combine_write_stats_glb(obctmpatu,filepars,varid,nobu,&
                                       & nobuprocs,indexobuprocs,&
                                       & varatts=varatts)
         ENDIF
         IF (nobv.GT.0) THEN
            varid = varid + 1
            CALL combine_write_stats_glb(obctmpatv,filepars,varid,nobv,&
                                       & nobvprocs,indexobvprocs,&
                                       & varatts=varatts)
         ENDIF
      ENDIF
      IF (iopt_obc_2D.EQ.1) THEN
         IF (nobu.GT.0) THEN
            varid = varid + 1
            CALL combine_write_stats_glb(obc2uvatu,filepars,varid,nobu,&
                                       & nobuprocs,indexobuprocs,&
                                       & varatts=varatts)
         ENDIF
         IF (nobv.GT.0) THEN
            varid = varid + 1
            CALL combine_write_stats_glb(obc2uvatv,filepars,varid,nobv,&
                                       & nobvprocs,indexobvprocs,&
                                       & varatts=varatts)
         ENDIF
      ENDIF
      IF (iopt_obc_3D.EQ.1) THEN
         IF (nobu.GT.0) THEN
            varid = varid + 1
            CALL combine_write_stats_glb(obc3uvatu,filepars,varid,nobu,&
                                       & nobuprocs,indexobuprocs,&
                                       & varatts=varatts)
         ENDIF
         IF (nobv.GT.0) THEN
            varid = varid + 1
            CALL combine_write_stats_glb(obc3uvatv,filepars,varid,nobv,&
                                       & nobvprocs,indexobvprocs,&
                                       & varatts=varatts)
         ENDIF
      ENDIF

!
!2.10 Next output index
!----------------------
!

      itrestart = itrestart + 1

   ENDIF

ENDIF

!
!3. Finalise
!-----------
!

IF ((cold_start.AND.ntrestart(1).EQ.0).OR.&
  & (itrestart.GT.norestarts)) THEN
   IF (master) CALL close_filepars(filepars)
   DEALLOCATE (varatts)
ENDIF

modfiles(io_fincon,ics_phys,2) = filepars

IF (flag) CALL log_timer_out() 


RETURN

END SUBROUTINE write_phsics
