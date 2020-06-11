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

SUBROUTINE usrdef_output
!************************************************************************
!
! *usrdef_output* User-formatted output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.11.1
!
! $Date: 2018-01-25 15:49:18 +0100 (Thu, 25 Jan 2018) $
!
! $Revision: 1079 $
!
! Description - output parameters for test case wcinlet
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_stats_loc, define_out2d_vals, error_alloc,
!                local_proc, mult_index, open_file, Uvar_at_C, Vvar_at_C
!
!************************************************************************
!
USE gridpars
USE iopars
USE modids
USE paralpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: local_proc
USE inout_routines, ONLY: close_file, open_file
USE model_output, ONLY: define_out2d_vals
USE paral_comms, ONLY: combine_stats_loc
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=1) :: cl 
INTEGER, SAVE :: icout, iunit, nostatsloc, novars
INTEGER :: i, iproc, istat, ivar,j, l, lloc, mtime
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: varids
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: istatpos, jstatpos, nostatsprocs
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: lstatsprocs
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: outdat, outdatglb


!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,1)%status = 'R'
   modfiles(io_modgrd,1,1)%form = modfiles(io_modgrd,1,2)%form
   modfiles(io_modgrd,1,1)%filename = modfiles(io_modgrd,1,2)%filename
   modfiles(io_modgrd,1,2)%status = '0'
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1,1)%status = 'R'
   modfiles(io_2uvobc,1,1)%form = modfiles(io_2uvobc,1,2)%form
   modfiles(io_2uvobc,1,1)%filename = modfiles(io_2uvobc,1,2)%filename
   modfiles(io_2uvobc,1,2)%status = '0'
ENDIF

IF (cold_start) RETURN

!
!1. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

!
!1.1 Open output file
!--------------------
!

   IF (master) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case wcinlet: '&
                      & //'simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF

!
!1.2 Number of local stations
!----------------------------
!

   ALLOCATE (nostatsprocs(nprocs),STAT=errstat)
   CALL error_alloc('nostatsprocs',2,(/nprocs/),int_type)
   nostatsprocs = 0

   nostatsloc = 0
   istat_120: DO istat=1,nostatstsr
      l = lstatstsr(2,istat)
      i = tsrstatlocs(l)%ipos; j = tsrstatlocs(l)%jpos
      iproc_121: DO iproc=1,nprocs
         IF (local_proc(i,j,iproc=iproc)) THEN
            nostatsprocs(iproc) = nostatsprocs(iproc) + 1
            IF (idloc.EQ.idprocs(iproc)) THEN
               nostatsloc = nostatsloc + 1
            ENDIF
         ENDIF
      ENDDO iproc_121
   ENDDO istat_120

!
!1.3 Allocate arrays
!-------------------
!

   novars = MERGE(4,15,iopt_waves.EQ.0)
   ALLOCATE (varids(novars),STAT=errstat)
   CALL error_alloc('varids',1,(/novars/),kndint)
   
   ALLOCATE (outdat(nostatsloc,novars),STAT=errstat)
   CALL error_alloc('outdat',2,(/nostatsloc,novars/),kndrtype)
   outdat = 0.0
   ALLOCATE (outdatglb(nostatstsr,novars),STAT=errstat)
   CALL error_alloc('outdatglb',2,(/nostatstsr,novars/),kndrtype)
   outdatglb = 0.0
   ALLOCATE (lstatsprocs(nostatstsr,nprocs),STAT=errstat)
   CALL error_alloc('lstatsprocs',2,(/nostatstsr,nprocs/),kndint)
   lstatsprocs = 0
   ALLOCATE (istatpos(nostatsloc),STAT=errstat)
   CALL error_alloc('istatpos',1,(/nostatsloc/),kndint)
   istatpos = 0
   ALLOCATE (jstatpos(nostatsloc),STAT=errstat)
   CALL error_alloc('jstatpos',1,(/nostatsloc/),kndint)
   jstatpos = 0

!
!1.4 Mapping array and local station locations
!---------------------------------------------
!

   iproc_140: DO iproc=1,nprocs
      lloc = 0
      istat_141: DO istat=1,nostatstsr
         l = lstatstsr(2,istat)
         i = tsrstatlocs(l)%ipos; j = tsrstatlocs(l)%jpos
         IF (local_proc(i,j,iproc=iproc)) THEN
            lloc = lloc + 1
            lstatsprocs(lloc,iproc) = l
            IF (idloc.EQ.idprocs(iproc)) THEN
               istatpos(lloc) = i - nc1loc + 1
               jstatpos(lloc) = j - nr1loc + 1
            ENDIF
         ENDIF
      ENDDO istat_141
   ENDDO iproc_140

!
!1.5 Output variables
!--------------------
!   

   IF (iopt_waves.EQ.0) THEN
      varids = (/iarr_zeta,iarr_umvel,iarr_vmvel,iarr_bstresatc/)
   ELSE
      varids = (/iarr_zeta,iarr_umvel,iarr_vmvel,iarr_umstokesatc,&
        & iarr_vmstokesatc,iarr_umveltot, iarr_vmveltot,iarr_waveheight,&
        & iarr_waveperiod,iarr_wavedir,iarr_bstresatc,iarr_bstresatc_max,&
        & iarr_bstresatc_wav,iarr_wavethickatc,iarr_zaroughatc/)
   ENDIF

!
!1.5 Output frequency
!--------------------
!

   icout = 240

   GOTO 1000

ENDIF

!
!3. Store parameters
!-------------------
!

IF (mult_index(nt,icout)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

   mtime = nosecsrun/3600.0

!
!3.1 Store data
!--------------
!
!  ---local
   outdat = 0.0
   l_310: DO l=1,nostatsloc
      i = istatpos(l); j = jstatpos(l)
      CALL define_out2d_vals(outdat(l,:),i,j,novars,ivarid=varids)
   ENDDO l_310

!  ---combine globally
   CALL combine_stats_loc(outdatglb,outdat,nostatstsr,nostatsprocs,&
                        & lstatsprocs,0)

!
!3.2 Write data
!--------------
!

   IF (master) THEN
      WRITE (iunit,9001) mtime
      l_320: DO l=1,nostatstsr
         WRITE (iunit,*)
         cl = tsrstatlocs(l)%name(1:1)
         WRITE (iunit,9002) 'zeta'//TRIM(cl), 100.0*outdatglb(l,1)
         WRITE (iunit,9002) 'umvel'//TRIM(cl), 100.0*outdatglb(l,2)
         WRITE (iunit,9002) 'vmvel'//TRIM(cl), 100.0*outdatglb(l,3)
         IF (iopt_waves.GT.0) THEN
            WRITE (iunit,9002) 'umstokes'//TRIM(cl), 100.0*outdatglb(l,4)
            WRITE (iunit,9002) 'vmstokes'//TRIM(cl), 100.0*outdatglb(l,5)
            WRITE (iunit,9002) 'umveltot'//TRIM(cl), 100.0*outdatglb(l,6)
            WRITE (iunit,9002) 'vmveltot'//TRIM(cl), 100.0*outdatglb(l,7)
            WRITE (iunit,9002) 'waveheight'//TRIM(cl), outdatglb(l,8)
            WRITE (iunit,9002) 'waveperiod'//TRIM(cl), outdatglb(l,9)
            WRITE (iunit,9002) 'wavedir'//TRIM(cl), outdatglb(l,10)
         ENDIF
         ivar = MERGE(4,11,iopt_waves.EQ.0)
         WRITE (iunit,9002) 'bstresatc'//TRIM(cl), outdatglb(l,ivar)
         IF (iopt_waves.GT.0) THEN
            WRITE (iunit,9002) 'bstresatc_max'//TRIM(cl), outdatglb(l,12)
            WRITE (iunit,9002) 'bstresatc_wav'//TRIM(cl), outdatglb(l,13)
            WRITE (iunit,9002) 'wavethick'//TRIM(cl), 100.0*outdatglb(l,14)
            WRITE (iunit,9002) 'zaroughatc'//TRIM(cl), 100.0*outdatglb(l,15)
         ENDIF
      ENDDO l_320
      IF (nt.LT.nstep) WRITE (iunit,*)
   ENDIF

ENDIF

!
!4. Finalise at last time step
!-----------------------------
!

IF (nt.EQ.nstep) THEN

!  ---close output file
   IF (master) CALL close_file(iunit,'A')

!  ---deallocate arrays
   DEALLOCATE(istatpos,jstatpos,lstatsprocs,nostatsprocs,outdat,outdatglb,&
            & varids)

ENDIF

1000 IF (mult_index(nt,icout)) CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',I2,' hours')
9002 FORMAT (T3,A,T18,': ',G12.4E3)

END SUBROUTINE usrdef_output
