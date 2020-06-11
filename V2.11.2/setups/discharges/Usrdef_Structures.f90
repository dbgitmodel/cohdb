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
!

!************************************************************************
!
! *Usrdef_Structures* User-defined model setup for structures
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! $Date: 2017-02-13 14:59:12 +0100 (Mon, 13 Feb 2017) $
!
! $Revision: 1004 $
!
! Description - test case discharges 
!
! Reference -
!
! Subroutines - usrdef_dry_cells, usrdef_thin_dams, usrdef_weirs,
!               usrdef_dischr_spec, usrdef_dischr_data
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_dry_cells
!************************************************************************
!
! *usrdef_dry_cells* Define locations of dry cells
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! Description - 
!
! Reference -
!
! Calling program - define_global_grid
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_dry_cells

!========================================================================

SUBROUTINE usrdef_thin_dams
!************************************************************************
!
! *usrdef_thin_dams* Define locations of thin dams
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_thin_dams

!========================================================================

SUBROUTINE usrdef_weirs
!************************************************************************
!
! *usrdef_weirs* Define locations and arrays for weirs and barriers
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_weirs

!========================================================================

SUBROUTINE usrdef_dischr_spec
!************************************************************************
!
! *usrdef_discharges* Define specifier arrays for discharge module
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! Description - test case discharges
!
! Reference -
!
! Calling program - define_dischr_spec
!
!************************************************************************
!
USE iopars
USE structures
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'usrdef_dischr_spec'
CALL log_timer_in()

SELECT CASE (runtitle(11:11))
   CASE ('A','B')
      kdistype = (/0,1,1,2/)
   CASE ('C','D')
      kdistype = 1
END SELECT

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_dischr_spec

!========================================================================

SUBROUTINE usrdef_dischr_data(iddesc,ifil,ciodatetime,disdata,nodat,novars)
!************************************************************************
!
! *usrdef_discharges* Define discharge data
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Structures.f90  V2.6
!
! Description - test case discharges
!
! Reference -
!
! Calling program - define_dischr_data
!
! Module calls - add_secs_to_data, lim_dims
! 
!************************************************************************
!
USE iopars
USE physpars
USE syspars
USE timepars
USE time_routines, ONLY: add_secs_to_date, log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(OUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
REAL, INTENT(INOUT), DIMENSION(nodat,novars) :: disdata

!
! Name         Type Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*disdata*     REAL    Input data
!*nodat*       INTEGER Number of discharge locations
!*novars*      INTEGER Number of input variables
!
!------------------------------------------------------------------------------
!
!* Local variables
!
LOGICAL :: time_series
INTEGER, SAVE :: icloc, icsal, icvol
INTEGER :: icount
REAL, SAVE :: deltdat, speed
REAL :: disvol1, disvol2, period, time


procname(pglev+1) = 'usrdef_dischr_data'
CALL log_timer_in()

time_series =  lim_dims(ABS(modfiles(iddesc,ifil,1)%tlims)).GT.1

!
!1. Initialise parameters on first call
!--------------------------------------
!

IF (modfiles(iddesc,ifil,1)%iostat.EQ.0) THEN
   modfiles(iddesc,ifil,1)%iostat = 1
   SELECT CASE (iddesc)
      CASE (io_disloc); icloc = 0
      CASE (io_disvol); icvol = 0
      CASE (io_dissal); icsal = 0
   END SELECT
   GOTO 1000
ENDIF

!
!2. Date/time
!------------
!

IF (time_series) THEN
   SELECT CASE (iddesc)
      CASE (io_disloc)
         icount = icloc
         deltdat = 60.0
      CASE (io_disvol)
         icount = icvol
         deltdat = 60.0
      CASE (io_dissal)
         icount = icsal
         deltdat = delt2d*(nstep/2)
   END SELECT
   CALL add_secs_to_date(CStartDateTime,ciodatetime,icount,deltdat)
ELSE
   ciodatetime = CStartDateTime
ENDIF

!
!3. Data values
!--------------
!

SELECT CASE (iddesc)
!  ---discharge locations
   CASE (io_disloc)
      disdata(:,1) = (/7500.0,12500.0,7500.0,12500.0/)
      disdata(:,2) = (/7500.0,7500.0,12500.0,12500.0/)
      SELECT CASE (runtitle(11:11))
         CASE ('A','B')
            disdata(:,3) = (/0.0,1.0,20.0,5.0/)
         CASE ('C','D')
            disdata(:,3) = 20
      END SELECT
      IF (time_series) THEN
         speed = 50.0/3.0
         time = icloc*deltdat
         disdata(1,1) = disdata(1,1) + time*speed
         disdata(2,2) = disdata(2,2) + time*speed
         disdata(3,2) = disdata(3,2) - time*speed
         disdata(4,1) = disdata(4,1) - time*speed
         icloc = icloc + 1
      ENDIF
!  ---discharge rates
   CASE (io_disvol)
      IF (time_series) THEN
         time = icloc*deltdat
         period = nstep*delt2d
         disvol1 = 50.0
         disvol2 = 150.0
         disdata(:,1) = ((period-time)*disvol1+time*disvol2)/period
         icvol = icvol + 1
      ELSE
         disdata(:,1) = 100.0
      ENDIF
!  ---discharge area and direction
   CASE (io_discur)
      disdata(:,1) = 10.0
      disdata(:,2) = (/0.0,halfpi,3.0*halfpi,pi/)
!  ---salinity discharge
   CASE (io_dissal)
      IF (time_series) THEN
         time = icsal*deltdat
         disdata(:,1) = MERGE(sal_ref,0.0,icsal.EQ.0)
         icsal = icsal + 1
      ELSE
         disdata(:,1) = sal_ref
      ENDIF
END SELECT

1000 CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_dischr_data
