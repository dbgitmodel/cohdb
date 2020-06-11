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
! Version - @(COHERENS)Usrdef_Output.f90  V2.11.2
!
! $Date: 2017-06-22 11:10:44 +0200 (Thu, 22 Jun 2017) $
!
! $Revision: 1037 $
!
! Description - output parameters for test case papa
!
! Reference -
!
! Calling program - initialse_model, coherens_main
!
! Module calls - cf90_get_var, cf90_inquire, cf90_inquire_variable,
!                close_file, error_alloc, open_file,
!
!************************************************************************
!
USE gridpars  
USE iopars
USE paralpars
USE syspars
USE timepars
USE cf90_routines, ONLY: cf90_get_var, cf90_inquire, cf90_inquire_variable
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=lenname) :: varname
INTEGER, PARAMETER :: norecs = 48, noutvars2d = 7
INTEGER, SAVE :: iunit
INTEGER :: iout2, iout3, irec, ivar, novars
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: out2ddat, out3ddat


IF (.NOT.master.OR.cold_start) RETURN

!
!1. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN

!  ---open output file
   CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
   WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
   WRITE (iunit,'(A)') 'Output parameters for test case papa: simulation '&
                        & //TRIM(runtitle)
   WRITE (iunit,*)
ENDIF

IF (nt.LT.nstep) RETURN

procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

!
!2. Allocate
!-----------
!

ALLOCATE (out2ddat(1,1,noutvars2d),STAT=errstat)
CALL error_alloc('out2ddat',3,(/1,1,noutvars2d/),real_type)
ALLOCATE (out3ddat(1,1,nz),STAT=errstat)
CALL error_alloc('out3ddat',3,(/1,1,nz/),real_type)

!
!3. Open data files
!------------------
!
!---2-D file
CALL open_file(iout2,avr2d(2)%filename,'IN','N')
!---3-D file
CALL open_file(iout3,avr3d(2)%filename,'IN','N')

!
!4. Read/store data
!------------------
!

irec_410: DO irec=1,norecs

!   
!4.1 Read data
!-------------
!   
!4.1.1 2-D data
!--------------
!   

   CALL cf90_inquire(iout2,nvariables=novars)
   ivar_411: DO ivar=1,novars
      CALL cf90_inquire_variable(iout2,ivar,varname)
      SELECT CASE (TRIM(varname))
         CASE ('tmld')
            CALL cf90_get_var(iout2,ivar,out2ddat(:,:,1),irec)
         CASE ('qrad')
            CALL cf90_get_var(iout2,ivar,out2ddat(:,:,2),irec)
         CASE ('qlwave')
            CALL cf90_get_var(iout2,ivar,out2ddat(:,:,3),irec)
         CASE ('qlatent')
            CALL cf90_get_var(iout2,ivar,out2ddat(:,:,4),irec)
         CASE ('qsensible')
            CALL cf90_get_var(iout2,ivar,out2ddat(:,:,5),irec)
         CASE ('qnonsol')
            CALL cf90_get_var(iout2,ivar,out2ddat(:,:,6),irec)
         CASE ('qtot')
            CALL cf90_get_var(iout2,ivar,out2ddat(:,:,7),irec)
      END SELECT
   ENDDO ivar_411

!
!4.1.2 3-D data
!--------------
!   

   CALL cf90_inquire(iout3,nvariables=novars)
   ivar_412: DO ivar=1,novars
      CALL cf90_inquire_variable(iout3,ivar,varname)
      SELECT CASE (TRIM(varname))
         CASE ('temp')
            CALL cf90_get_var(iout3,ivar,out3ddat,irec)
      END SELECT
   ENDDO ivar_412

!
!4.2 Write test case parameters
!------------------------------  
!  

   WRITE (iunit,9001) irec
   WRITE (iunit,9002) 'tsur', out3ddat(1,1,nz)
   WRITE (iunit,9002) 'temp40', out3ddat(1,1,260)
   WRITE (iunit,9002) 'temp120', out3ddat(1,1,180)
   WRITE (iunit,9002) 'tbot', out3ddat(1,1,1)
   WRITE (iunit,9002) 'depmix', out2ddat(1,1,1)
   WRITE (iunit,9002) 'qrad', out2ddat(1,1,2)
   WRITE (iunit,9002) 'qlwave', out2ddat(1,1,3)
   WRITE (iunit,9002) 'qlatent', out2ddat(1,1,4)
   WRITE (iunit,9002) 'qsensible', out2ddat(1,1,5)
   WRITE (iunit,9002) 'qnonsol', out2ddat(1,1,6)
   WRITE (iunit,9002) 'qtot', out2ddat(1,1,7)

ENDDO irec_410

!
!5. Close files
!---------------------
!

CALL close_file(iunit,'A')
CALL close_file(iout2,'N')
CALL close_file(iout3,'N')

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (out2ddat,out3ddat)

CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',I2,' months') 
9002 FORMAT (T3,A,T13,': ',G12.4E3)

END SUBROUTINE usrdef_output
