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
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.11.1
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - output parameters for test case flocvprof
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, mult_indesx, open_file
!
!************************************************************************
!
USE currents
USE fluxes
USE gridpars
USE iopars
USE physpars
USE sedarrays
USE syspars
USE timepars
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, SAVE :: icount, iunit
REAL :: rhour



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
!  ---surface forcing
   modfiles(io_1uvsur,1:2,1)%status = 'R'
   modfiles(io_1uvsur,1:2,1)%form = &
 & modfiles(io_1uvsur,1:2,2)%form
   modfiles(io_1uvsur,1:2,1)%filename = modfiles(io_1uvsur,1:2,2)%filename
   modfiles(io_1uvsur,1:2,2)%status = '0'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,1)%status = 'R'
   modfiles(io_sedspc,1,1)%form =   modfiles(io_sedspc,1,2)%form
   modfiles(io_sedspc,1,1)%filename = modfiles(io_sedspc,1,2)%filename
   modfiles(io_sedspc,1,2)%status = '0'
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (cold_start.OR.runtitle(10:10).EQ.'0') RETURN

IF (nt.EQ.0) THEN
!  ---open output file
   CALL open_file(iunit,runtitle(1:9)//runtitle(11:11)//'.'//suffix,'OUT','A')
   WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
   WRITE (iunit,'(A)') 'Output parameters for test case flocvprof: '&
                        & //TRIM(runtitle)
   WRITE (iunit,*)
   
!  ---time counter
   icount = 75
   
ENDIF

!
!3. Write output parameters
!--------------------------
!

IF (nt.GT.0.AND.mult_index(nt,icount)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

   rhour = nosecsrun/3600.0
   WRITE (iunit,9001) rhour
   WRITE (iunit,9002) 'umean', umvel(1,1)
   WRITE (iunit,9002) 'ubot', uvel(1,1,1)
   WRITE (iunit,9002) 'usur', uvel(1,1,nz)
   WRITE (iunit,9002) 'bstres', density_ref*bstresatc(1,1)
   WRITE (iunit,9002) 'floc_P_mean', rhos(1)*SUM(cvol(1,1,:,1))/REAL(nz)
   WRITE (iunit,9002) 'floc_F_mean', &
                     & SUM(floc_dens(1,1,1)*cvol(1,1,:,2))/REAL(nz)
   WRITE (iunit,9002) 'floc_dia_mean', 1.0E+06*SUM(floc_dia(1,1,:))/REAL(nz)
   WRITE (iunit,9002) 'floc_dens_mean', SUM(floc_dens(1,1,:))/REAL(nz)
   WRITE (iunit,9002) 'floc_ws_mean', SUM(wfall(1,1,:,2))/REAL(nz)
   WRITE (iunit,9002) 'floc_P_bot', rhos(1)*cvol(1,1,1,1)
   WRITE (iunit,9002) 'floc_F_bot', floc_dens(1,1,1)*cvol(1,1,1,2)
   WRITE (iunit,9002) 'floc_dia_bot', 1.0E+06*floc_dia(1,1,1)
   WRITE (iunit,9002) 'floc_dens_bot', floc_dens(1,1,1)
   WRITE (iunit,9002) 'floc_ws_bot', wfall(1,1,1,2)
   WRITE (iunit,9002) 'floc_P_sur', rhos(1)*cvol(1,1,nz,1)
   WRITE (iunit,9002) 'floc_F_sur', floc_dens(1,1,1)*cvol(1,1,nz,2)
   WRITE (iunit,9002) 'floc_dia_sur', 1.0E+06*floc_dia(1,1,nz)
   WRITE (iunit,9002) 'floc_dens_sur', floc_dens(1,1,nz)
   WRITE (iunit,9002) 'floc_ws_sur', wfall(1,1,nz,2)

   CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (nt.EQ.nstep) CALL close_file(iunit,'A')


RETURN

9001 FORMAT ('Time: ',F5.2,' hours')
9002 FORMAT (T3,A,T19,': ',G12.4E3)

END SUBROUTINE usrdef_output





