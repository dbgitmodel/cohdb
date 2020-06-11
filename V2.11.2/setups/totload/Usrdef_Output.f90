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
! Version - @(COHERENS)Usrdef_Output.f90  V2.10.3
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - output parameters for test case totload
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, mult_index, open_file
!
!************************************************************************
!
USE currents
USE fluxes
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
INTEGER, PARAMETER :: notests = 7
INTEGER :: iunit2, nn
INTEGER, SAVE :: icount, itest, iunit, n
REAL, SAVE, DIMENSION(10) :: bstresatc_cst, vvel_cst
REAL, SAVE, DIMENSION(10,notests) :: qtotvals

!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
   modfiles(io_sedspc,1,1)%status = 'R'
   modfiles(io_sedspc,1,1)%form = modfiles(io_sedspc,1,2)%form
   modfiles(io_sedspc,1,1)%filename = modfiles(io_sedspc,1,2)%filename
   modfiles(io_sedspc,1,2)%status = '0'
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (cold_start) RETURN

IF (nt.EQ.0) THEN
   icount = modfiles(io_inicon,ics_phys,1)%tlims(3)
!  ---open output file
   CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
   WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
   WRITE (iunit,'(A)') 'Output parameters for test case totload: '&
                      & //'simulation '//TRIM(runtitle)
   WRITE (iunit,*)
!  ---counters
   n = 0
   SELECT CASE (runtitle(8:8))
      CASE ('A'); itest = 1
      CASE ('B'); itest = 2
      CASE ('C'); itest = 3
      CASE ('D'); itest = 4
      CASE ('E'); itest = 5
      CASE ('F'); itest = 6
      CASE ('G'); itest = 7
   END SELECT

ENDIF

!
!3. Write output parameters
!--------------------------
!

IF (mult_index(nt+1,icount)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

!  ---write output
   WRITE (iunit,9001) nosecsrun
   WRITE (iunit,9002) 'vmean', vdvel(2,1)/depmean_cst
   WRITE (iunit,9002) 'bstres', density_ref*bstresatc(1,1)
   WRITE (iunit,9002) 'qtot', rhos(1)*qtotatv(2,1,1)

!  ---parameters for table format
   n = n + 1
   IF (itest.EQ.1) THEN
      vvel_cst(n) = vdvel(2,1)/depmean_cst
      bstresatc_cst(n) = density_ref*bstresatc(1,1)
   ENDIF
   qtotvals(n,itest) = rhos(1)*qtotatv(2,1,1)

   CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (nt.EQ.nstep) THEN 
   CALL close_file(iunit,'A')
   IF (isimul.EQ.nosimul) THEN
      CALL open_file(iunit2,runtitle(1:7)//'.dat','OUT','A')
      WRITE (iunit2,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit2,'(A)') 'Output data for test case totload'
      WRITE (iunit2,*)
      WRITE (iunit2,'(A)') 'vmean bstres      EH1         EH2         AW'//&
                         & '          WU          VR03        VR07        GM'
      WRITE (iunit2,*)
      nn_410: DO nn=1,10
         WRITE (iunit2,9003) vvel_cst(nn), bstresatc_cst(nn), qtotvals(nn,:)
      ENDDO nn_410
      CALL close_file(iunit2,'A')
   ENDIF
ENDIF

RETURN

9001 FORMAT ('Time: ',I3,' seconds') 
9002 FORMAT (T3,A,T12,': ',G12.4E3)
9003 format (F4.1,8(2X,G10.4))

END SUBROUTINE usrdef_output
