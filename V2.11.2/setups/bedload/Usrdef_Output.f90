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
! Version - @(COHERENS)Usrdef_Output.f90  V2.5
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - output parameters for test case bedload
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
INTEGER, PARAMETER :: notests = 6
INTEGER :: iunit2, nn
INTEGER, SAVE :: icount, itest, iunit, n
REAL, SAVE, DIMENSION(10) :: bstresatc_cst, uvel_cst
REAL, SAVE, DIMENSION(10,notests) :: qbedvals


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
   WRITE (iunit,'(A)') 'Output parameters for test case bedload: simulation '&
                     & //TRIM(runtitle)
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
   WRITE (iunit,9003) nosecsrun
   WRITE (iunit,9001) 'umean', udvel(2,1)/depmean_cst
   WRITE (iunit,9001) 'bstres', density_ref*bstresatc(1,1)
   WRITE (iunit,9001) 'qbed', rhos(1)*qbedatu(2,1,1)

!  ---parameters for table format
   n = n + 1
   IF (itest.EQ.1) THEN
      uvel_cst(n) = udvel(2,1)/depmean_cst
      bstresatc_cst(n) = density_ref*bstresatc(1,1)
   ENDIF
   qbedvals(n,itest) = rhos(1)*qbedatu(2,1,1)

   CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (nt.EQ.nstep) THEN 
   CALL close_file(iunit,'A')
   IF (itest.EQ.notests) THEN
      CALL open_file(iunit2,runtitle(1:7)//'.dat','OUT','A')
      WRITE (iunit2,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit2,'(A)') 'Output data for test case bedload'
      WRITE (iunit2,*)
      WRITE (iunit2,'(A)') 'umean bstres      MPM         EF          VR84'//&
                         & '        WU          SO          VR03'
      WRITE (iunit2,*)
      nn_410: DO nn=1,10
         WRITE (iunit2,9002) uvel_cst(nn), bstresatc_cst(nn), qbedvals(nn,:)
      ENDDO nn_410
      CALL close_file(iunit2,'A')
   ENDIF
ENDIF

RETURN

9001 FORMAT (T3,A,T12,': ',G12.4E3)
9002 FORMAT (F4.1,8(2X,G10.4))
9003 FORMAT ('Time: ',I3,' seconds') 

END SUBROUTINE usrdef_output
