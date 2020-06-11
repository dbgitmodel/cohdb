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
! Description - output parameters for test case wavload
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, mult_index, open_file, vector_mag_var_atc
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
USE wavevars
USE inout_routines, ONLY: close_file, open_file
USE math_library, ONLY: complex_polar, vector_mag_var_atc
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, PARAMETER :: notests = 4
INTEGER :: iunit2, iunit3, nn
INTEGER, SAVE :: icount, itest, iunit, n
REAL :: qloaddir, qloadmag
REAL, SAVE, DIMENSION(900) :: bstresatc_cst, uvel_cst, wdir, wheight, wperiod
REAL, SAVE, DIMENSION(900,4) :: qdirvals, qmagvals


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
   icount = modfiles(io_wavsur,1,1)%tlims(3)
!  ---open output file
   CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
   WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
   WRITE (iunit,'(A)') 'Output parameters for test case wavload: '&
                   & //'simulation '//TRIM(runtitle)
!  ---counters
   n = 0
   SELECT CASE (runtitle(8:8))
      CASE ('A'); itest = 1
      CASE ('B'); itest = 2
      CASE ('C'); itest = 3
      CASE ('D'); itest = 4
   END SELECT

ENDIF

!
!3. Write output parameters
!--------------------------
!

IF (mult_index(nt+1,3)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

!  ---write output
   WRITE (iunit,*)
   WRITE (iunit,9001) 'umean', umvel(2,1)
   WRITE (iunit,9001) 'bstres', density_ref*bstresatc(1,1)
   WRITE (iunit,9001) 'wheight', waveheight(1,1)
   WRITE (iunit,9001) 'wperiod', waveperiod(1,1)
   WRITE (iunit,9001) 'wdir', radtodeg*wavedir(1,1)

   SELECT CASE (runtitle(8:8))
      CASE ('A')
         CALL complex_polar(qbedatu(2,1,1),qbedatv(1,2,1),xamp=qloadmag,&
                          & xpha=qloaddir)
         WRITE (iunit,9001) 'qloadmag', qloadmag
         WRITE (iunit,9001) 'qloaddir', radtodeg*qloaddir
         WRITE (iunit,9001) 'qbedu', rhos(1)*qbedatu(2,1,1)
         WRITE (iunit,9001) 'qbedv', rhos(1)*qbedatv(1,2,1)
      CASE ('B','C','D')
         CALL complex_polar(qtotatu(2,1,1),qtotatv(1,2,1),xamp=qloadmag,&
                          & xpha=qloaddir)
         WRITE (iunit,9001) 'qloadmag', qloadmag
         WRITE (iunit,9001) 'qloaddir', radtodeg*qloaddir
         WRITE (iunit,9001) 'qtotu', rhos(1)*qtotatu(2,1,1)
         WRITE (iunit,9001) 'qtotv', rhos(1)*qtotatv(1,2,1)
   END SELECT

!  ---parameters for table format
   n = n + 1
   IF (itest.EQ.1) THEN
      uvel_cst(n) = umvel(2,1)
      bstresatc_cst(n) = density_ref*bstresatc(1,1)
      wheight(n) = waveheight(1,1)
      wperiod(n) = waveperiod(1,1)
      wdir(n) = radtodeg*wavedir(1,1)
   ENDIF
   qmagvals(n,itest) = rhos(1)*qloadmag
   qdirvals(n,itest) = radtodeg*qloaddir

   CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (nt.EQ.nstep) then
   CALL close_file(iunit,'A')
   IF (isimul.EQ.nosimul) THEN
      CALL open_file(iunit2,runtitle(1:7)//'_mag.dat','OUT','A')
      WRITE (iunit2,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit2,'(A)') 'Load magnitudes for test case wavload'
      WRITE (iunit2,*)
      WRITE (iunit2,'(A)') 'umean bstres      wheight     wperiod     wdir'//&
                         & '        SO          MG          VR03        VR07'
      WRITE (iunit2,*)
      CALL open_file(iunit3,runtitle(1:7)//'_dir.dat','OUT','A')
      WRITE (iunit3,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit3,'(A)') 'Load directions for test case wavload'
      WRITE (iunit3,*)
      WRITE (iunit3,'(A)') 'umean bstres      wheight     wperiod     wdir'//&
                         & '        SO          MG          VR03        VR07'
      WRITE (iunit3,*)
      nn_410: DO nn=1,900
         WRITE (iunit2,9002) uvel_cst(nn), bstresatc_cst(nn), wheight(nn), &
                           & wperiod(nn), wdir(nn), qmagvals(nn,:)
         WRITE (iunit3,9002) uvel_cst(nn), bstresatc_cst(nn), wheight(nn), &
                           & wperiod(nn), wdir(nn), qdirvals(nn,:)
      ENDDO nn_410
      CALL close_file(iunit2,'A')
      CALL close_file(iunit3,'A')
   ENDIF
ENDIF


RETURN

9001 FORMAT (T3,A,T12,': ',G12.4E3)
9002 format (F4.1,8(2X,G10.4))

END SUBROUTINE usrdef_output
