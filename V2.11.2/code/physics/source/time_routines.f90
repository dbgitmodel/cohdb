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

MODULE time_routines
!************************************************************************
!
! *time_routines* Date/time routines
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.11
!
! Generic routines - add_secs_to_date, check_date, convert_date, date_to_year,
!                    day_number, diff_dates, num_time_steps, year_to_date
!
! Defined operators - .earlier., .later., .noearlier., .nolater.
!
! Routines - add_secs_to_phase, clock_date_time, diff_clock,
!            error_lbound_var_date, error_ubound_var_date, initialise_time,
!            leap_year, log_timer_in, log_timer_out, read_clock, suspend_proc,
!            update_time
!
!************************************************************************
!
USE iopars
USE syspars
USE timepars

IMPLICIT NONE

!
! Interfaces
!
INTERFACE add_secs_to_date
   MODULE PROCEDURE add_secs_to_date_char, add_secs_to_date_char_d, &
                  & add_secs_to_date_int, add_secs_to_date_int_d
END INTERFACE

INTERFACE check_date
   MODULE PROCEDURE check_date_char, check_date_int
END INTERFACE

INTERFACE convert_date
   MODULE PROCEDURE convert_date_to_char, convert_date_to_int
END INTERFACE

INTERFACE date_to_year
   MODULE PROCEDURE date_to_year_char, date_to_year_int
END INTERFACE

INTERFACE day_number
   MODULE PROCEDURE day_number_char, day_number_int
END INTERFACE

INTERFACE diff_dates
   MODULE PROCEDURE diff_dates_char, diff_dates_int
END INTERFACE

INTERFACE num_time_steps
   MODULE PROCEDURE num_time_steps_char, num_time_steps_int
END INTERFACE

INTERFACE year_to_date
   MODULE PROCEDURE year_to_date_char, year_to_date_int
END INTERFACE

INTERFACE OPERATOR (.earlier.)
   MODULE PROCEDURE earlier_char, earlier_int
END INTERFACE

INTERFACE OPERATOR (.later.)
   MODULE PROCEDURE later_char, later_int
END INTERFACE

INTERFACE OPERATOR (.noearlier.)
   MODULE PROCEDURE noearlier_char, noearlier_int
END INTERFACE

INTERFACE OPERATOR (.nolater.)
   MODULE PROCEDURE nolater_char, nolater_int
END INTERFACE

CONTAINS

!========================================================================

SUBROUTINE add_secs_to_date_char(chardate1,chardate2,numsteps,dsecs)
!************************************************************************
!
! *add_secs_to_date_char* Add 'numsteps' time steps of 'dsecs' secs to
!                         input date/time
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.11
!
! Description - returns new date/time in character format
!             - dsecs given in single precision  
!
! Module calls - add_secs_to_date_int, convert_date
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: chardate1
CHARACTER (LEN=lentime), INTENT(OUT) :: chardate2
INTEGER, INTENT(IN) :: numsteps
REAL (KIND=kndreal), INTENT(IN) :: dsecs

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*chardate1* CHAR    Old date/time
!*chardate2* CHAR    New date/time
!*numsteps*  INTEGER Number of time steps
!*dsecs*     REAL    Number of secs in one time step                        [s]
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER, DIMENSION(7) :: intdate1, intdate2


procname(pglev+1) = 'add_secs_to_date_char'
CALL log_timer_in()

!
!1. Convert to integer format
!----------------------------
!

intdate1 = convert_date(chardate1)

!
!2. New date/time
!----------------
!

CALL add_secs_to_date_int(intdate1,intdate2,numsteps,dsecs)

!
!3. Reconvert to character format
!--------------------------------
!

chardate2 = convert_date(intdate2)

CALL log_timer_out()


RETURN

END SUBROUTINE add_secs_to_date_char

!========================================================================

SUBROUTINE add_secs_to_date_char_d(chardate1,chardate2,numsteps,dsecs)
!************************************************************************
!
! *add_secs_to_date_char_d* Add 'numsteps' time steps of 'dsecs' secs to
!                           input date/time
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.11
!
! Description - returns new date/time in character format
!             - dsecs given in double precision  
!
! Module calls - add_secs_to_date_int, convert_date
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: chardate1
CHARACTER (LEN=lentime), INTENT(OUT) :: chardate2
INTEGER, INTENT(IN) :: numsteps
REAL (KIND=kndrlong), INTENT(IN) :: dsecs

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*chardate1* CHAR    Old date/time
!*chardate2* CHAR    New date/time
!*numsteps*  INTEGER Number of time steps
!*dsecs*     DOUBLE  Number of secs in one time step                        [s]
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER, DIMENSION(7) :: intdate1, intdate2


procname(pglev+1) = 'add_secs_to_date_char_d'
CALL log_timer_in()

!
!1. Convert to integer format
!----------------------------
!

intdate1 = convert_date(chardate1)

!
!2. New date/time
!----------------
!

CALL add_secs_to_date_int_d(intdate1,intdate2,numsteps,dsecs)

!
!3. Reconvert to character format
!--------------------------------
!

chardate2 = convert_date(intdate2)

CALL log_timer_out()


RETURN

END SUBROUTINE add_secs_to_date_char_d

!========================================================================

SUBROUTINE add_secs_to_date_int(intdate1,intdate2,numsteps,dsecs)
!************************************************************************
!
! *add_secs_to_date_int* Add 'numsteps' time steps of 'dsecs' secs to
!                        input date/time
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.11
!
! Description - returns new date/time in integer format
!
! Module calls - leap_year
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(IN) :: numsteps
REAL (KIND=kndreal), INTENT(IN) :: dsecs
INTEGER, INTENT(IN), DIMENSION(7) :: intdate1
INTEGER, INTENT(OUT), DIMENSION(7) :: intdate2

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*intdate1*  INTEGER  Old date/time
!*intdate2*  INTEGER  New date/time
!*numsteps*  INTEGER  Number of time steps
!*dsecs*     REAL     Number of secs in one time step                       [s]
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: ileap, l, m, nday2, nhours2, nmins2, nmont2, nodays, &
         & nohours, nomins, nosecs, nsecs2, nyear2
INTEGER (KIND=kndilong) :: msecs2 


procname(pglev+1) = 'add_secs_to_date_int'
CALL log_timer_in()

nosecs = numsteps*INT(dsecs)

!
!1. Zero addition
!----------------
!

IF (numsteps.EQ.0.OR.dsecs.EQ.0) THEN
   intdate2 = intdate1

!
!2. Add seconds
!--------------
!

ELSEIF (numsteps.GT.0) THEN

!  ---milliseconds
   msecs2 = intdate1(7) + numsteps*INT(1000*(dsecs-INT(dsecs)),KIND=kndilong)
   IF (msecs2.GE.1000) THEN
      m = MOD(msecs2,1000_kndilong)
      nosecs = nosecs + (msecs2-m)/1000
      msecs2 = m
   ENDIF

!  ---seconds
   nsecs2 = intdate1(6) + nosecs
   nomins = 0
   IF (nsecs2.GE.60) THEN
      m = MOD(nsecs2,60)
      nomins = (nsecs2-m)/60
      nsecs2 = m
   ENDIF

!  ---minutes
   nmins2 = intdate1(5) + nomins
   nohours = 0
   IF (nmins2.GE.60) THEN
      m = MOD(nmins2,60)
      nohours = (nmins2-m)/60
      nmins2 = m
   ENDIF

!  ---hours
   nhours2 = intdate1(4) + nohours
   nodays = 0
   IF (nhours2.GE.24) THEN
      m = MOD(nhours2,24)
      nodays = (nhours2-m)/24
      nhours2 = m
   ENDIF

!  ---year
   nyear2 = intdate1(1)
   ileap = leap_year(nyear2)
   nodays = nodays + monthdays(intdate1(2)) + intdate1(3) - 1
   IF (ileap.EQ.1.AND.intdate1(2).GT.2) nodays = nodays + 1
   DO WHILE (nodays.GE.(365+ileap))
      nodays = nodays - 365 - ileap
      nyear2 = nyear2 + 1
      ileap = leap_year(nyear2)
   ENDDO

!  ---month
   l_210: DO l=2,13
      IF (nodays.LT.(monthdays(l)+ileap)) EXIT l_210
   ENDDO l_210
   nmont2 = l - 1
   IF (ileap.EQ.1.AND.nodays.EQ.31) nmont2 = 2

!  ---day
   IF (ileap.EQ.1.AND.nmont2.GT.2) THEN
      nday2 = nodays + 1 - monthdays(nmont2) - ileap
   ELSE
      nday2 = nodays + 1 - monthdays(nmont2)
   ENDIF

!  ---store in date/time array
   intdate2(1:6) = (/nyear2,nmont2,nday2,nhours2,nmins2,nsecs2/)
   intdate2(7) = msecs2

!
!3. Substract seconds
!--------------------
!

ELSEIF (numsteps.LT.0) THEN

!  ---milliseconds
   msecs2 = intdate1(7) + numsteps*INT(1000*(dsecs-INT(dsecs)),KIND=kndilong)
   IF (msecs2.LT.0) THEN
      m = MOD(msecs2,1000_kndilong)
      IF (m.EQ.0) THEN
         nosecs = nosecs + msecs2/1000
         msecs2 = 0
      ELSE
         nosecs = nosecs + (msecs2-m)/1000 - 1
         msecs2 = 1000 + m
      ENDIF
   ENDIF

!  ---seconds
   nsecs2 = intdate1(6) + nosecs
   nomins = 0
   IF (nsecs2.LT.0) THEN
      m = MOD(nsecs2,60)
      IF (m.EQ.0) THEN
         nomins = nsecs2/60
         nsecs2 = 0
      ELSE
         nomins = (nsecs2-m)/60 - 1
         nsecs2 = 60 + m
      ENDIF
   ENDIF

!  ---minutes
   nmins2 = intdate1(5) + nomins
   nohours = 0
   IF (nmins2.LT.0) THEN
      m = MOD(nmins2,60)
      IF (m.EQ.0) THEN
         nohours = nmins2/60
         nmins2 = 0
      ELSE
         nohours = (nmins2-m)/60 - 1
         nmins2 = 60 + m
      ENDIF
   ENDIF

!  ---hours
   nhours2 = intdate1(4) + nohours
   nodays = 0
   IF (nhours2.LT.0) THEN
      m = MOD(nhours2,24)
      IF (m.EQ.0) THEN
         nodays = nhours2/24
         nhours2 = 0
      ELSE
         nodays = (nhours2-m)/24 - 1
         nhours2 = 24 + m
      ENDIF
   ENDIF

!  ---year
   nyear2 = intdate1(1)
   ileap = leap_year(nyear2)
   nodays = nodays + monthdays(intdate1(2)) + intdate1(3) - 1
   IF (ileap.EQ.1.AND.intdate1(2).GT.2) nodays = nodays + 1
   DO WHILE (nodays.LT.0)
      ileap = leap_year(nyear2-1)
      nodays = nodays + 365 + ileap
      nyear2 = nyear2 - 1
   ENDDO

!  ---month
   l_310: DO l=2,13
      IF (nodays.LT.(monthdays(l)+ileap)) EXIT l_310
   ENDDO l_310
   nmont2 = l - 1
   IF ((ileap.EQ.1).AND.(nodays.EQ.31)) nmont2 = 2

!  ---day
   IF ((ileap.EQ.1).AND.(nmont2.GT.2)) THEN
      nday2 = nodays + 1 - monthdays(nmont2) - ileap
   ELSE
      nday2 = nodays + 1 - monthdays(nmont2)
   ENDIF

!  ---store in date/time array
   intdate2(1:6) = (/nyear2,nmont2,nday2,nhours2,nmins2,nsecs2/)
   intdate2(7) = msecs2

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE add_secs_to_date_int

!========================================================================

SUBROUTINE add_secs_to_date_int_d(intdate1,intdate2,numsteps,dsecs)
!************************************************************************
!
! *add_secs_to_date_int_d* Add 'numsteps' time steps of 'dsecs' secs to
!                        input date/time
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.11
!
! Description - returns new date/time in integer format
!             - dsecs given in double precision  
!
! Module calls - leap_year
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(IN) :: numsteps
REAL (KIND=kndrlong), INTENT(IN) :: dsecs
INTEGER, INTENT(IN), DIMENSION(7) :: intdate1
INTEGER, INTENT(OUT), DIMENSION(7) :: intdate2

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*intdate1*  INTEGER  Old date/time
!*intdate2*  INTEGER  New date/time
!*numsteps*  INTEGER  Number of time steps
!*dsecs*     REAL     Number of secs in one time step                       [s]
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: ileap, l, m, nday2, nhours2, nmins2, nmont2, nodays, &
         & nohours, nomins, nosecs, nsecs2, nyear2
INTEGER (KIND=kndilong) :: msecs2 


procname(pglev+1) = 'add_secs_to_date_int'
CALL log_timer_in()

nosecs = numsteps*INT(dsecs)

!
!1. Zero addition
!----------------
!

IF (numsteps.EQ.0.OR.dsecs.EQ.0) THEN
   intdate2 = intdate1

!
!2. Add seconds
!--------------
!

ELSEIF (numsteps.GT.0) THEN

!  ---milliseconds
   msecs2 = intdate1(7) + numsteps*INT(1000*(dsecs-INT(dsecs)),KIND=kndilong)
   IF (msecs2.GE.1000) THEN
      m = MOD(msecs2,1000_kndilong)
      nosecs = nosecs + (msecs2-m)/1000
      msecs2 = m
   ENDIF

!  ---seconds
   nsecs2 = intdate1(6) + nosecs
   nomins = 0
   IF (nsecs2.GE.60) THEN
      m = MOD(nsecs2,60)
      nomins = (nsecs2-m)/60
      nsecs2 = m
   ENDIF

!  ---minutes
   nmins2 = intdate1(5) + nomins
   nohours = 0
   IF (nmins2.GE.60) THEN
      m = MOD(nmins2,60)
      nohours = (nmins2-m)/60
      nmins2 = m
   ENDIF

!  ---hours
   nhours2 = intdate1(4) + nohours
   nodays = 0
   IF (nhours2.GE.24) THEN
      m = MOD(nhours2,24)
      nodays = (nhours2-m)/24
      nhours2 = m
   ENDIF

!  ---year
   nyear2 = intdate1(1)
   ileap = leap_year(nyear2)
   nodays = nodays + monthdays(intdate1(2)) + intdate1(3) - 1
   IF (ileap.EQ.1.AND.intdate1(2).GT.2) nodays = nodays + 1
   DO WHILE (nodays.GE.(365+ileap))
      nodays = nodays - 365 - ileap
      nyear2 = nyear2 + 1
      ileap = leap_year(nyear2)
   ENDDO

!  ---month
   l_210: DO l=2,13
      IF (nodays.LT.(monthdays(l)+ileap)) EXIT l_210
   ENDDO l_210
   nmont2 = l - 1
   IF (ileap.EQ.1.AND.nodays.EQ.31) nmont2 = 2

!  ---day
   IF (ileap.EQ.1.AND.nmont2.GT.2) THEN
      nday2 = nodays + 1 - monthdays(nmont2) - ileap
   ELSE
      nday2 = nodays + 1 - monthdays(nmont2)
   ENDIF

!  ---store in date/time array
   intdate2(1:6) = (/nyear2,nmont2,nday2,nhours2,nmins2,nsecs2/)
   intdate2(7) = msecs2

!
!3. Substract seconds
!--------------------
!

ELSEIF (numsteps.LT.0) THEN

!  ---milliseconds
   msecs2 = intdate1(7) + numsteps*INT(1000*(dsecs-INT(dsecs)),KIND=kndilong)
   IF (msecs2.LT.0) THEN
      m = MOD(msecs2,1000_kndilong)
      IF (m.EQ.0) THEN
         nosecs = nosecs + msecs2/1000
         msecs2 = 0
      ELSE
         nosecs = nosecs + (msecs2-m)/1000 - 1
         msecs2 = 1000 + m
      ENDIF
   ENDIF

!  ---seconds
   nsecs2 = intdate1(6) + nosecs
   nomins = 0
   IF (nsecs2.LT.0) THEN
      m = MOD(nsecs2,60)
      IF (m.EQ.0) THEN
         nomins = nsecs2/60
         nsecs2 = 0
      ELSE
         nomins = (nsecs2-m)/60 - 1
         nsecs2 = 60 + m
      ENDIF
   ENDIF

!  ---minutes
   nmins2 = intdate1(5) + nomins
   nohours = 0
   IF (nmins2.LT.0) THEN
      m = MOD(nmins2,60)
      IF (m.EQ.0) THEN
         nohours = nmins2/60
         nmins2 = 0
      ELSE
         nohours = (nmins2-m)/60 - 1
         nmins2 = 60 + m
      ENDIF
   ENDIF

!  ---hours
   nhours2 = intdate1(4) + nohours
   nodays = 0
   IF (nhours2.LT.0) THEN
      m = MOD(nhours2,24)
      IF (m.EQ.0) THEN
         nodays = nhours2/24
         nhours2 = 0
      ELSE
         nodays = (nhours2-m)/24 - 1
         nhours2 = 24 + m
      ENDIF
   ENDIF

!  ---year
   nyear2 = intdate1(1)
   ileap = leap_year(nyear2)
   nodays = nodays + monthdays(intdate1(2)) + intdate1(3) - 1
   IF (ileap.EQ.1.AND.intdate1(2).GT.2) nodays = nodays + 1
   DO WHILE (nodays.LT.0)
      ileap = leap_year(nyear2-1)
      nodays = nodays + 365 + ileap
      nyear2 = nyear2 - 1
   ENDDO

!  ---month
   l_310: DO l=2,13
      IF (nodays.LT.(monthdays(l)+ileap)) EXIT l_310
   ENDDO l_310
   nmont2 = l - 1
   IF ((ileap.EQ.1).AND.(nodays.EQ.31)) nmont2 = 2

!  ---day
   IF ((ileap.EQ.1).AND.(nmont2.GT.2)) THEN
      nday2 = nodays + 1 - monthdays(nmont2) - ileap
   ELSE
      nday2 = nodays + 1 - monthdays(nmont2)
   ENDIF

!  ---store in date/time array
   intdate2(1:6) = (/nyear2,nmont2,nday2,nhours2,nmins2,nsecs2/)
   intdate2(7) = msecs2

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE add_secs_to_date_int_d

!========================================================================

SUBROUTINE add_secs_to_phase(phasein,phaseout,numsteps,dsecs,freq)
!************************************************************************
!
! *add_secs_to_phase* Calculate phasein+MOD(numsteps*dsecs*freq,twopi)
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - algorithm avoids rounding errors
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(IN) :: numsteps
REAL, INTENT(IN) :: dsecs, freq, phasein
REAL, INTENT(OUT) :: phaseout

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*phasein*   REAL     Old phase
!*phaseout*  REAL     New phase
!*numsteps*  INTEGER  Number of time steps
!*dsecs*     REAL     Number of secs in one time step                       [s]
!*freq*      REAL     Frequency                                         [rad/s]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iperiod, isign, msecs, n, noperiods, nosecs
REAL :: dmsecs


msecs = INT(numsteps)*NINT(1000*dsecs)
nosecs = msecs/1000
dmsecs = 0.001*(msecs-1000*nosecs)
isign = SIGN(1,INT(nosecs))
nosecs = ABS(nosecs)
iperiod = twopi/freq
noperiods = nosecs/iperiod
phaseout = 0.0
n_110: DO n=1,noperiods
   phaseout = phaseout + freq*iperiod
   phaseout = MOD(phaseout,twopi)
ENDDO n_110
nosecs = MOD(nosecs,iperiod)
phaseout = phaseout + freq*(nosecs+dmsecs)
phaseout = MOD(phaseout,twopi)
phaseout = MOD(phasein+isign*phaseout,twopi)
IF (phaseout.LT.0.0) phaseout = phaseout + twopi


RETURN

END SUBROUTINE add_secs_to_phase
  
!========================================================================

SUBROUTINE check_date_char(chardate,varname)
!************************************************************************
!
! *check_date_char* Check date/time variable
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.1.1
!
! Description - input date/time is in character format
!
! Module calls - convert_date, error_lbound_arr, error_limits_arr,
!                leap_year
!
!************************************************************************
!
USE error_routines, ONLY: error_lbound_arr, error_limits_arr

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
CHARACTER (LEN=lentime), INTENT(IN) :: chardate

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*chardate*  CHAR     Date/time
!*varname*   CHAR     Name of date/time variable
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ileap, ndays
INTEGER, DIMENSION(7) :: intdate


procname(pglev+1) = 'check_date_char'
CALL log_timer_in()

intdate = convert_date(chardate)
CALL error_lbound_arr(intdate(1),varname,0,.FALSE.,1,indx=(/1/))
CALL error_limits_arr(intdate(2),varname,1,12,1,indx=(/2/))
ileap = leap_year(intdate(1))
ndays = days_in_month(intdate(2))
IF (ileap.EQ.1.AND.intdate(2).GT.1) ndays = ndays + 1
CALL error_limits_arr(intdate(3),varname,1,ndays,1,indx=(/3/))
CALL error_limits_arr(intdate(4),varname,0,23,1,indx=(/4/))
CALL error_limits_arr(intdate(5),varname,0,59,1,indx=(/5/))
CALL error_limits_arr(intdate(6),varname,0,59,1,indx=(/6/))
CALL error_limits_arr(intdate(7),varname,0,999,1,indx=(/7/))

CALL log_timer_out()


RETURN

END SUBROUTINE check_date_char

!========================================================================

SUBROUTINE check_date_int(intdate,varname)
!************************************************************************
!
! *check_date_int* Check date/time variable
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.1.1
!
! Description - input date/time is in integer format
!
! Module calls - error_lbound_arr, error_limits_arr, leap_year
!
!************************************************************************
!
USE error_routines, ONLY: error_lbound_arr, error_limits_arr

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN), DIMENSION(7) :: intdate

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*intdate*   INTEGER  Date/time
!*varname*   CHAR     Name of date/time variable
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ileap, ndays


procname(pglev+1) = 'check_date_int'
CALL log_timer_in()

CALL error_lbound_arr(intdate(1),varname,0,.FALSE.,1,indx=(/1/))
CALL error_limits_arr(intdate(2),varname,1,12,1,indx=(/2/))
ileap = leap_year(intdate(1))
ndays = days_in_month(intdate(2))
IF (ileap.EQ.1.AND.intdate(2).GT.1) ndays = ndays + 1
CALL error_limits_arr(intdate(3),varname,1,ndays,1,indx=(/3/))
CALL error_limits_arr(intdate(4),varname,0,23,1,indx=(/4/))
CALL error_limits_arr(intdate(5),varname,0,59,1,indx=(/5/))
CALL error_limits_arr(intdate(6),varname,0,59,1,indx=(/6/))
CALL error_limits_arr(intdate(7),varname,0,999,1,indx=(/7/))

CALL log_timer_out()
  

RETURN

END SUBROUTINE check_date_int

!========================================================================

SUBROUTINE clock_date_time(charclock)
!************************************************************************
!
! *clock_date_time* Returns date/time from systems's real-time clock
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description -
!
! Module calls -
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(OUT) :: charclock

!
! Name       Type   Purpose
!------------------------------------------------------------------------------
!*charclock* CHAR   Date/time from system clock   
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=8) :: cdate
CHARACTER (LEN=10) :: ctime


CALL DATE_AND_TIME(DATE=cdate,TIME=ctime)
charclock = cdate(1:4)//'/'//cdate(5:6)//'/'//cdate(7:8)//';'//ctime(1:2)//&
          & ':'//ctime(3:4)//':'//ctime(5:6)//','//ctime(8:10)
  

RETURN

END SUBROUTINE clock_date_time

!========================================================================

FUNCTION convert_date_to_char(intdate,seps) RESULT(chardate)
!************************************************************************
!
! *convert_date_to_char* Convert date/time from integer to character format
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.1.0
!
! Description -
!
! Module calls - error_abort
!
!************************************************************************
!
USE error_routines, ONLY: error_abort

!
!* Arguments
!
INTEGER, INTENT(IN), DIMENSION(7) :: intdate
CHARACTER (LEN=4), INTENT(IN), OPTIONAL :: seps
CHARACTER (LEN=lentime) :: chardate

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intdate*   INTEGER Date/time in integer format
!*seps*      CHAR    Separators used in the date/time string
!*chardate*  CHAR    Date/time in character format
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=4) :: sepsx
INTEGER :: istat


procname(pglev+1) = 'convert_date_to_char'
CALL log_timer_in()

!
!1. Separators in date/time string
!---------------------------------
!

IF (PRESENT(seps)) THEN
   sepsx = seps
ELSE
   sepsx = '/;::'
ENDIF

!
!2. Convert
!----------
!

istat = 0
WRITE (chardate,9001,IOSTAT=istat) intdate
IF (istat.NE.0) THEN
   nerrs = 1
   IF (errchk) WRITE (ioerr,9002) 'Invalid Date/Time: ', intdate
   CALL error_abort('convert_date_to_char',ierrno_write)
ENDIF
chardate(5:5) = sepsx(1:1); chardate(8:8) = sepsx(1:1)
chardate(11:11) = sepsx(2:2); chardate(14:14) = sepsx(3:3)
chardate(17:17) = sepsx(3:3); chardate(20:20) = sepsx(4:4)

CALL log_timer_out()


RETURN

9001 FORMAT(I4.4,1X,I2.2,1X,I2.2,1X,I2.2,1X,I2.2,1X,I2.2,1X,I3.3)
9002 FORMAT(A,7(I6,1X))

END FUNCTION convert_date_to_char

!========================================================================

FUNCTION convert_date_to_int(chardate) RESULT(intdate)
!************************************************************************
!
! *convert_date_to_int* Convert date/time from character to integer format
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description -
!
! Module calls - error_abort
!
!************************************************************************
!
USE error_routines, ONLY: error_abort

!
!* Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: chardate
INTEGER, DIMENSION(7) :: intdate

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*chardate*  CHAR    Date/time in character format
!*intdate*   INTEGER Date/time in integer format
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: istat


procname(pglev+1) = 'convert_date_to_int'
CALL log_timer_in()

READ (chardate,'(I4,5(TR1,I2),TR1,I3)',IOSTAT=istat) intdate
IF (istat.NE.0) THEN
   nerrs = 1
   IF (errchk) WRITE (ioerr,'(A)') 'Invalid Date/Time: '//chardate
   CALL error_abort('convert_date_to_char',ierrno_read)
ENDIF

CALL log_timer_out()


RETURN

END FUNCTION convert_date_to_int

!========================================================================

SUBROUTINE date_to_year_char(chardate,yearnum)
!************************************************************************
!
! *date_to_year_char* Convert date/time in year date
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.7.1
!
! Description - input date is in character format
!
! Module calls - convert_date, leap_year
!
!************************************************************************
!
!* Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: chardate
REAL (KIND=kndrlong), INTENT(OUT) :: yearnum

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*chardate*  CHAR    Date/time
!*yearnum*   REAL    Time in years
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ileap
REAL :: daynum
INTEGER, DIMENSION(7) :: intdate


procname(pglev+1) = 'date_to_year_char'
CALL log_timer_in()

intdate = convert_date(chardate)
ileap = leap_year(intdate(1))
daynum = monthdays(intdate(2))+intdate(3) - 1 +&
                & (intdate(4)+(intdate(5)+(intdate(6)+&
                &  intdate(7)/1000.0)/60.0)/60.0)/24.0 + ileap
yearnum = intdate(1) + daynum/(365+ileap)

CALL log_timer_out()


RETURN

END SUBROUTINE date_to_year_char

!========================================================================

SUBROUTINE date_to_year_int(intdate,yearnum)
!************************************************************************
!
! *date_to_year_int* Convert date/time in year date
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.7.1
!
! Description - input date is in integer format
!
! Module calls - leap_year
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(IN), DIMENSION(7) :: intdate
REAL (KIND=kndrlong), INTENT(OUT) :: yearnum

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intdate*   INTEGER Date/time
!*yearnum*   REAL    Time in years
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ileap
REAL :: daynum


procname(pglev+1) = 'date_to_year_int'
CALL log_timer_in()

ileap = leap_year(intdate(1))
daynum = monthdays(intdate(2))+intdate(3) - 1 +&
                & (intdate(4)+(intdate(5)+(intdate(6)+&
                &  intdate(7)/1000.0)/60.0)/60.0)/24.0 + ileap
yearnum = intdate(1) + daynum/(365+ileap)

CALL log_timer_out()


RETURN

END SUBROUTINE date_to_year_int

!========================================================================

SUBROUTINE day_number_char(chardate,daynum)
!************************************************************************
!
! *day_number_char* Return day number of the year
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - result is between 0 and 364 or 365 (leap year)
!             - input date/time in character format
!
! Module calls - convert_date, day_number_int
!
!************************************************************************
!
!* Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: chardate
REAL, INTENT(OUT) :: daynum

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*chardate*  CHAR    Date/time
!*daynum*    REAL    Day number
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER, DIMENSION(7) :: intdate


procname(pglev+1) = 'day_number_char'
CALL log_timer_in()

!
!1. Convert date to integer format
!---------------------------------
!

intdate = convert_date(chardate)

!
!2. Day number
!-------------
!

CALL day_number_int(intdate,daynum)

CALL log_timer_out()


RETURN

END SUBROUTINE day_number_char

!========================================================================

SUBROUTINE day_number_int(intdate,daynum)
!************************************************************************
!
! *day_number_int* Return day number of the year
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - result is between 0 and 364 or 365 (leap year)
!             - input date/time in integer format
!
! Module calls - leap_year
!
!************************************************************************
!
!* Arguments
!
REAL, INTENT(OUT) :: daynum
INTEGER, INTENT(IN), DIMENSION(7) :: intdate

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intdate*   INTEGER Date/time
!*daynum*    REAL    Day number
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'day_number_int'
CALL log_timer_in()

daynum = monthdays(intdate(2))+intdate(3) - 1 +&
                & (intdate(4)+(intdate(5)+(intdate(6)+&
                &  intdate(7)/1000.0)/60.0)/60.0)/24.0
IF ((leap_year(intdate(1)).EQ.1).AND.(intdate(2).GT.2)) daynum = daynum+1

CALL log_timer_out()


RETURN

END SUBROUTINE day_number_int

!========================================================================

FUNCTION diff_clock(npccold)
!************************************************************************
!
! *diff_clock* Return number of clock counts since a previous call to
!              process clock
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description -
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(IN) :: npccold
INTEGER :: diff_clock

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*npccold*    INTEGER Old clock count
!*diff_clock* INTEGER New minus old clock count 
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: npccnew


CALL SYSTEM_CLOCK(COUNT=npccnew)
IF (npccnew.GE.npccold) THEN
   diff_clock = npccnew - npccold
ELSE
   diff_clock = npcc_max + npccnew - npccold
ENDIF


RETURN

END FUNCTION diff_clock

!========================================================================

SUBROUTINE diff_dates_char(chardate1,chardate2,tunit,nosecs,millisecs,rtime)
!************************************************************************
!
! *diff_dates_char* Return time difference between two character dates
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.7.1
!
! Description - if tunit = 0, result is integer and returned in seconds
!                (nosecs) and (optionally) milliseconds
!             - if tunit > 0, result is real and returned in rtime
!
! Module calls - convert_date, diff_dates_int
!
!************************************************************************
!
!* Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: chardate1, chardate2
INTEGER, INTENT(IN) :: tunit
INTEGER, INTENT(OUT), OPTIONAL :: millisecs
INTEGER (KIND=kndilong), INTENT(OUT), OPTIONAL :: nosecs
REAL (KIND=kndrlong), INTENT(OUT), OPTIONAL :: rtime

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*chardate1* CHAR   Start date
!*chardate2* CHAR   End date
!*tunit*     INTEGER Type of result
!              = 0   => seconds and milliseconds
!              = 1   => seconds
!              = 2   => minutes
!              = 3   => hours
!              = 4   => days
!              = 5   => months
!              = 6   => years
!*nosecs*    LONGINT Result in seconds
!*millisecs* INTEGER Residual milleseconds
!*rtime*     DOUBLE  Result in appropriate units
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER, DIMENSION(7) :: intdate1, intdate2


procname(pglev+1) = 'diff_dates_char'
CALL log_timer_in()

!
!1. Convert dates to integer format
!----------------------------------
!

intdate1 = convert_date(chardate1)
intdate2 = convert_date(chardate2)

!
!2. Number of seconds
!--------------------
!

IF (tunit.EQ.0) THEN
   CALL diff_dates_int(intdate1,intdate2,0,nosecs,millisecs)
ELSE
   CALL diff_dates_int(intdate1,intdate2,tunit,rtime=rtime)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE diff_dates_char

!========================================================================

SUBROUTINE diff_dates_int(intdate1,intdate2,tunit,nosecs,millisecs,rtime)
!************************************************************************
!
! *diff_dates_int* Return time difference between two integer dates
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.11
!
! Description - if tunit = 0, result is integer and returned in seconds
!               (nosecs) and (optionally) milliseconds
!             - if tunit > 0, result is real and returned in rtime
!
! Module calls - day_number, .earlier., .later., leap_year
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(IN) :: tunit
INTEGER, INTENT(OUT), OPTIONAL :: millisecs
INTEGER (KIND=kndilong), INTENT(OUT), OPTIONAL :: nosecs
REAL (KIND=kndrlong), INTENT(OUT), OPTIONAL :: rtime
INTEGER, INTENT(IN), DIMENSION(7) :: intdate1, intdate2

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intdate1*  INTEGER Start date
!*intdate2*  INTEGER End date
!*tunit*     INTEGER Type of result
!              = 0   => seconds and milliseconds
!              = 1   => seconds
!              = 2   => minutes
!              = 3   => hours
!              = 4   => days
!              = 5   => weeks
!              = 6   => months
!              = 7   => years
!*nosecs*    LONGINT Result in seconds
!*millisecs* INTEGER Residual milleseconds
!*rtime*     DOUBLE  Result in appropriate units
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: isign, iyear, iyear1, iyear2, ndays1, ndays2, tconv
INTEGER (KIND=kndilong) :: nodays, numsecs
REAL :: daynum1, daynum2
INTEGER, DIMENSION(7) :: idate1, idate2


procname(pglev+1) = 'diff_dates_int'
CALL log_timer_in()

!
!1. Determine which date comes first
!-----------------------------------
!

IF (intdate1.earlier.intdate2) THEN
   idate1 = intdate1; idate2 = intdate2; isign = 1
ELSEIF (intdate1.later.intdate2) THEN
   idate1 = intdate2; idate2 = intdate1; isign = -1
ELSE
   isign = 0
ENDIF

IF (isign.EQ.0) THEN
   IF (PRESENT(nosecs)) nosecs = 0
   IF (PRESENT(millisecs)) millisecs = 0
   IF (PRESENT(rtime)) rtime = 0.0
   GOTO 1000
ENDIF

!
!2. Number of days between start/end dates
!-----------------------------------------
!

iyear1 = idate1(1)
iyear2 = idate2(1)
CALL day_number(idate1,daynum1)
CALL day_number(idate2,daynum2)
ndays1 = daynum1; ndays2 = daynum2
IF (iyear1.EQ.iyear2) THEN
   nodays = ndays2 - ndays1
ELSE
   nodays = 365 + leap_year(iyear1) - ndays1
   iyear_110: DO iyear=iyear1+1,iyear2-1
      nodays = nodays + 365 + leap_year(iyear)
   ENDDO iyear_110
   nodays = nodays + ndays2
ENDIF

!
!3. Number of seconds
!--------------------
!

numsecs = nodays*86400_kndilong-((idate1(4)*60+idate1(5))*60+&
        & idate1(6))+(idate2(4)*60+idate2(5))*60+idate2(6)

!
!4. Return result in appropriate format
!--------------------------------------
!

IF (tunit.EQ.0) THEN
   nosecs = isign*numsecs
   IF (PRESENT(millisecs)) THEN
      millisecs = idate2(7) - idate1(7)
      IF (millisecs.LT.0) THEN
         millisecs = millisecs + 1000
         nosecs = nosecs + isign
      ENDIF
      millisecs = isign*millisecs
   ENDIF
ELSE
   tconv = REAL(time_convert(tunit))
   rtime =  numsecs/tconv + MOD(numsecs,INT(tconv,KIND=kndilong))/REAL(tconv) +&
          & 0.001*(idate2(7)-idate1(7))/REAL(tconv,KIND=kndrlong)
   rtime = isign*rtime
ENDIF

1000 CALL log_timer_out()


RETURN

END SUBROUTINE diff_dates_int

!========================================================================

FUNCTION earlier_char(chardate1,chardate2)
!************************************************************************
!
! *earlier_char* Checks whether date 'chardate1' is earlier than 'chardate2'
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - input date/times in character format
!
! Module calls - convert_date
!
!************************************************************************
!
!* Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: chardate1, chardate2
LOGICAL :: earlier_char

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*chardate1*    CHAR    First date/time
!*chardate2*    CHAR    Second date/time
!*earlier_char* LOGICAL .TRUE. if chardate1 is earlier than chardate2
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: l
INTEGER, DIMENSION(7) :: intdate1, intdate2


intdate1 = convert_date(chardate1)
intdate2 = convert_date(chardate2)

earlier_char = .FALSE.
l_110: DO l=1,7
   IF (intdate1(l).NE.intdate2(l)) THEN
      earlier_char = MERGE(.TRUE.,.FALSE.,intdate1(l).LT.intdate2(l))
      EXIT l_110
   ENDIF
ENDDO l_110


RETURN

END FUNCTION earlier_char

!========================================================================

FUNCTION earlier_int(intdate1,intdate2)
!************************************************************************
!
! *earlier_int* checks whether date 'intdate1' is earlier than 'intdate2'
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - input date/times in integer format
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(IN), DIMENSION(7) :: intdate1, intdate2
LOGICAL :: earlier_int

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*intdate1*    INTEGER First date/time
!*intdate2*    INTEGER Second date/time
!*earlier_int* LOGICAL .TRUE. if intdate1 is earlier than intdate2
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: l


earlier_int = .FALSE.
l_110: DO l=1,7
   IF (intdate1(l).NE.intdate2(l)) THEN
      earlier_int = MERGE(.TRUE.,.FALSE.,intdate1(l).LT.intdate2(l))
      EXIT l_110
   ENDIF
ENDDO l_110


RETURN

END FUNCTION earlier_int

!========================================================================

SUBROUTINE error_lbound_var_date(cdate,varname,cdatemin,matchmin)
!************************************************************************
!
! *error_lbound_var_date* Checks whether date/time 'cdate' is not earlier than
!                         'cdatemin' if 'matchmin' is .TRUE. or is later than
!                         'cdatemin' if 'matchmin' is .FALSE.
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0  
!
! Description - displays error message if needed
!
! Module calls - .earlier., .nolater.
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmin
CHARACTER (LEN=*), INTENT(IN) :: varname
CHARACTER (LEN=lentime), INTENT(IN) :: cdate, cdatemin

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*cdate*    CHAR    Input date/time
!*varname*  CHAR    Name of date/time variable
!*cdatemin* CHAR    Date/time used for comparison
!*matchmin* LOGICAL Allows matching with cdatemin if .TRUE.
!
!------------------------------------------------------------------------------
!

IF (matchmin) THEN
   IF (cdate.earlier.cdatemin) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (ioerr,'(A)') 'Invalid value for date/time parameter '&
                           & //TRIM(varname)//': '//cdate
         WRITE (ioerr,'(A)') 'Must be not earlier than: '//cdatemin
      ENDIF
   ENDIF
ELSEIF (cdate.nolater.cdatemin) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (ioerr,'(A)') 'Invalid value for date/time parameter '&
                        & //TRIM(varname)//': '//cdate
      WRITE (ioerr,'(A)') 'Must be later than: '//cdatemin
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_lbound_var_date

!========================================================================

SUBROUTINE error_ubound_var_date(cdate,varname,cdatemax,matchmax)
!************************************************************************
!
! *error_ubound_var_date* Checks whether date/time cdate is not later than
!                         'cdatemax' if 'matchmax' is .TRUE. or is earlier
!                         than 'cdatemax' if 'matchmax' is .FALSE.
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - displays error message if needed
!
! Module calls - .later., .noearlier.
!
!************************************************************************
!
!*Arguments
!
LOGICAL, INTENT(IN) :: matchmax
CHARACTER (LEN=*), INTENT(IN) :: varname
CHARACTER (LEN=lentime), INTENT(IN) :: cdate, cdatemax

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*cdate*    CHAR    Input date/time
!*varname*  CHAR    Name of date/time variable
!*cdatemax* CHAR    Date/time used for comparison
!*matchmax* LOGICAL Allows matching with cdatemin if .TRUE.
!
!------------------------------------------------------------------------------
!

IF (matchmax) THEN
   IF (cdate.later.cdatemax) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (ioerr,'(A)') 'Invalid value for date/time parameter '&
                           & //TRIM(varname)//': '//cdate
         WRITE (ioerr,'(A)') 'Must be not later than: '//cdatemax
      ENDIF
   ENDIF
ELSEIF (cdate.noearlier.cdatemax) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (ioerr,'(A)') 'Invalid value for date/time parameter '&
                        & //TRIM(varname)//': '//cdate
      WRITE (ioerr,'(A)') 'Must be earlier than: '//cdatemax
   ENDIF
ENDIF


RETURN

END SUBROUTINE error_ubound_var_date

!========================================================================

SUBROUTINE initialise_time
!************************************************************************
!
! *initialise_time* Initialise date/time parameters
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - convert_date, day_number
!
!************************************************************************
!
USE switches

!
!*Local variables
!
REAL :: daynum


procname(pglev+1) = 'initialise_time'
CALL log_timer_in()

!---time counter
nt = 0

!---date/time
CDateTime = CStartDateTime
IStartDateTime = convert_date(CStartDateTime)
IDateTime = IStartDateTime
IEndDateTime = convert_date(CEndDateTime)
CALL day_number(IStartDateTime,daynum)
julianday = daynum + 1
nosecsrun = 0

!---time steps
metstep = iopt_meteo.EQ.1
wavestep = iopt_waves.EQ.1
physinit = modfiles(io_inicon,ics_phys,1)%status.NE.'0'

CALL log_timer_out()


RETURN

END SUBROUTINE initialise_time

!========================================================================

FUNCTION later_char(chardate1,chardate2)
!************************************************************************
!
! *later_char* Checks whether date 'chardate1' is later than 'chardate2'
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - input date/times in character format
!
! Module calls - convert_date
!
!************************************************************************
!
!* Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: chardate1, chardate2
LOGICAL :: later_char

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*chardate1*  CHAR    First date/time
!*chardate2*  CHAR    Second date/time
!*later_char* LOGICAL .TRUE. if chardate1 is later than chardate2
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: l
INTEGER, DIMENSION(7) :: intdate1, intdate2


intdate1 = convert_date(chardate1)
intdate2 = convert_date(chardate2)

later_char = .FALSE.
l_110: DO l=1,7
   IF (intdate1(l).NE.intdate2(l)) THEN
      later_char = MERGE(.TRUE.,.FALSE.,intdate1(l).GT.intdate2(l))
      EXIT l_110
   ENDIF
ENDDO l_110


RETURN

END FUNCTION later_char

!========================================================================

FUNCTION later_int(intdate1,intdate2)
!************************************************************************
!
! *later_int* Checks whether date 'intdate1' is later than 'intdate2'
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - input date/times in integer format
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(IN), DIMENSION(7) :: intdate1, intdate2
LOGICAL :: later_int

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*intdate1*   INTEGER First date/time
!*intdate2*   INTEGER Second date/time
!*later_int*  LOGICAL .TRUE. if intdate1 is later than intdate2
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: l


later_int = .FALSE.
l_110: DO l=1,7
   IF (intdate1(l).NE.intdate2(l)) THEN
      later_int = MERGE(.TRUE.,.FALSE.,intdate1(l).GT.intdate2(l))
      EXIT l_110
   ENDIF
ENDDO l_110


RETURN

END FUNCTION later_int

!========================================================================

FUNCTION leap_year(iyear)
!************************************************************************
!
! *leap_year* Return 1 or 0 if 'iyear' is a leap year or not 
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description -
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(IN) :: iyear
INTEGER :: leap_year

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*iyear*      INTEGER Year
!*leap_year*  INTEGER Number of leap days
!
!------------------------------------------------------------------------------
!


IF (MOD(iyear,4).EQ.0.AND.(MOD(iyear,100).NE.0.OR.MOD(iyear,400).EQ.0)) THEN
   leap_year = 1
ELSE
   leap_year = 0
ENDIF


RETURN

END FUNCTION leap_year

!========================================================================

SUBROUTINE log_timer_in(npcc,ivarid,logname,numvar,info)
!************************************************************************
!
! *log_timer_in* Write log info and initialise timer on entry of a procedure
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.10.1
!
! Description -
!
! Module calls - inquire_var, read_clock
!
!************************************************************************
!
USE modvars_routines, ONLY: inquire_var

!
!* Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: info
CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: logname
INTEGER, INTENT(IN), OPTIONAL :: ivarid, numvar
INTEGER, INTENT(OUT), OPTIONAL :: npcc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*npcc*       INTEGER Clock count
!*ivarid*     INTEGER Variable id
!*logname*    CHAR    String used in log output if present, otherwise name of
!                     calling program is taken
!*numvar*     INTEGER Variable number in case of a multi-variable key id
!*info*       LOGICAL No log info is written if present and .FALSE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: infox
CHARACTER (LEN=lenname) :: f90_name
CHARACTER (LEN=12) :: cnum
CHARACTER (LEN=200) :: pname
INTEGER :: ivaridx, l


pglev = pglev + 1

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(info)) THEN
   infox = info
ELSE
   infox = .TRUE.
ENDIF

IF (PRESENT(ivarid)) THEN
   ivaridx = ivarid
ELSE
   ivaridx = 0
ENDIF

IF (PRESENT(logname)) THEN
   l = LEN(logname)
   pname = logname(1:l)
ELSE
   pname = procname(pglev)
ENDIF

IF (PRESENT(numvar)) THEN
   WRITE (cnum,'(I12)') numvar
   cnum = ADJUSTL(cnum)
ELSE
   cnum = ''
ENDIF

!
!2. Write info
!-------------
!

IF (infox.AND.loglev1.GE.pglev) THEN
   IF (ivaridx.GT.0) THEN
      CALL inquire_var(ivaridx,f90_name=f90_name)
      IF (pglev.LT.10) THEN
         WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, &
                      & TRIM(pname)//': '//TRIM(f90_name)//TRIM(' '//TRIM(cnum))
      ELSE
         WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, &
                      & TRIM(pname)//': '//TRIM(f90_name)//TRIM(' '//TRIM(cnum))
      ENDIF
   ELSE
      IF (pglev.LT.10) THEN
         WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, TRIM(pname)
      ELSE
         WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, TRIM(pname)
      ENDIF
   ENDIF
ENDIF

!
!3. Timer
!--------
!

IF (timer) THEN
   IF (PRESENT(npcc)) npcc = read_clock()
ENDIF


RETURN

END SUBROUTINE log_timer_in

!========================================================================

SUBROUTINE log_timer_out(npcc,itm_type,info)
!************************************************************************
!
! *log_timer_out* Write log info and evaluate timer on exit of a procedure
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.10.1
!
! Description -
!
! Module calls - diff_clock
!
!************************************************************************
!
!* Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: info
INTEGER, INTENT(IN), OPTIONAL :: itm_type, npcc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*npcc*       INTEGER Old clock count
!*itm_type*   INTEGER Timer type
!*info*       LOGICAL No log info is written if present and .FALSE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: infox


IF (timer) THEN
   IF (PRESENT(npcc).AND.PRESENT(itm_type)) THEN
      nopcc(itm_type) = nopcc(itm_type) + diff_clock(npcc)
   ENDIF
ENDIF

IF (PRESENT(info)) THEN
   infox = info
ELSE
   infox = .TRUE.
ENDIF
   
IF (infox.AND.loglev2.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, logexit
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, logexit
   ENDIF
ENDIF

pglev = pglev - 1


RETURN

END SUBROUTINE log_timer_out

!========================================================================

FUNCTION noearlier_char(chardate1,chardate2)
!************************************************************************
!
! *noearlier_char* Checks whether date 'chardate1' is not earlier than
!                  'chardate2'
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - input date/times in character format
!
! Module calls - convert_date
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: chardate1, chardate2
LOGICAL :: noearlier_char

!
! Name            Type    Purpose
!------------------------------------------------------------------------------
!*chardate1*      CHAR    First date/time
!*chardate2*      CHAR    Second date/time
!*noearlier_char* LOGICAL .TRUE. if chardate1 is not earlier than chardate2
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: l
INTEGER, DIMENSION(7) :: intdate1, intdate2


intdate1 = convert_date(chardate1)
intdate2 = convert_date(chardate2)

noearlier_char = .TRUE.
l_110: DO l=1,7
   IF (intdate1(l).NE.intdate2(l)) THEN
      noearlier_char = MERGE(.TRUE.,.FALSE.,intdate1(l).GT.intdate2(l))
      EXIT l_110
   ENDIF
ENDDO l_110


RETURN

END FUNCTION noearlier_char

!========================================================================

FUNCTION noearlier_int(intdate1,intdate2)
!************************************************************************
!
! *noearlier_int* Checks whether date 'intdate1' is not earlier than 'intdate2'
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - input date/times in integer format
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(IN), DIMENSION(7) :: intdate1, intdate2
LOGICAL :: noearlier_int

!
! Name            Type    Purpose
!------------------------------------------------------------------------------
!*intdate1*       INTEGER First date/time
!*intdate2*       INTEGER Second date/time
!*noearlier_int*  LOGICAL .TRUE. if intdate1 is not earlier than intdate2
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: l


noearlier_int = .TRUE.
l_110: DO l=1,7
   IF (intdate1(l).NE.intdate2(l)) THEN
      noearlier_int = MERGE(.TRUE.,.FALSE.,intdate1(l).GT.intdate2(l))
      EXIT l_110
   ENDIF
ENDDO l_110


RETURN

END FUNCTION noearlier_int

!========================================================================

FUNCTION nolater_char(chardate1,chardate2)
!************************************************************************
!
! *nolater_char* Checks whether date 'chardate1' is not later than 'chardate2'
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - input date/times in character format
!
! Module calls - convert_date
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: chardate1, chardate2
LOGICAL :: nolater_char

!
! Name            Type    Purpose
!------------------------------------------------------------------------------
!*chardate1*      CHAR    First date/time
!*chardate2*      CHAR    Second date/time
!*nolater_char*   LOGICAL .TRUE. if chardate1 is not later than chardate2
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: l
INTEGER, DIMENSION(7) :: intdate1, intdate2


intdate1 = convert_date(chardate1)
intdate2 = convert_date(chardate2)

nolater_char = .TRUE.
l_110: DO l=1,7
   IF (intdate1(l).NE.intdate2(l)) THEN
      nolater_char = MERGE(.TRUE.,.FALSE.,intdate1(l).LT.intdate2(l))
      EXIT l_110
   ENDIF
ENDDO l_110


RETURN

END FUNCTION nolater_char

!========================================================================

FUNCTION nolater_int(intdate1,intdate2)
!************************************************************************
!
! *nolater_int* Checks whether date 'intdate1' is not later than 'intdate2'
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - input date/times in integer format
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(IN), DIMENSION(7) :: intdate1, intdate2
LOGICAL :: nolater_int

!
! Name            Type    Purpose
!------------------------------------------------------------------------------
!*intdate1*       INTEGER First date/time
!*intdate2*       INTEGER Second date/time
!*nolater_int*    LOGICAL .TRUE. if intdate1 is not later than intdate2
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: l


nolater_int = .TRUE.
l_110: DO l=1,7
   IF (intdate1(l).NE.intdate2(l)) THEN
      nolater_int = MERGE(.TRUE.,.FALSE.,intdate1(l).LT.intdate2(l))
      EXIT l_110
   ENDIF
ENDDO l_110


RETURN

END FUNCTION nolater_int

!========================================================================

SUBROUTINE num_time_steps_char(chardate1,chardate2,dsecs,numsteps)
!************************************************************************
!
! *num_time_steps_char* Return the number of time steps between two dates
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - input date/times in character format
!
! Module calls - convert_date, num_time_steps_int
!
!************************************************************************
!
!* Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: chardate1, chardate2
INTEGER, INTENT(OUT) :: numsteps
REAL, INTENT(IN) :: dsecs

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*chardate1* CHAR    First date/time
!*chardate2* CHAR    Second date/time
!*dsecs*     REAL    Number of secs in one time step                        [s]
!*numsteps*  INTEGER Number of time steps
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER, DIMENSION(7) :: intdate1, intdate2


procname(pglev+1) = 'num_time_steps_char'
CALL log_timer_in()

!
!1. Convert dates to integer format
!----------------------------------
!

intdate1 = convert_date(chardate1)
intdate2 = convert_date(chardate2)

!
!2. Number of time steps
!-----------------------
!

CALL num_time_steps_int(intdate1,intdate2,dsecs,numsteps)

CALL log_timer_out()


RETURN

END SUBROUTINE num_time_steps_char

!========================================================================

SUBROUTINE num_time_steps_int(intdate1,intdate2,dsecs,numsteps)
!************************************************************************
!
! *num_time_steps_int* Return the number of time steps between two dates
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description - input date/times in integer format
!
! Module calls - day_number, diff_dates, .earlier., .later., leap_year
!
!************************************************************************
!
!* Arguments
!
INTEGER, INTENT(OUT) :: numsteps
REAL, INTENT(IN) :: dsecs
INTEGER, INTENT(IN), DIMENSION(7) :: intdate1, intdate2

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*intdate1*  INTEGER First date/time
!*intdate2*  INTEGER Second date/time
!*dsecs*     REAL    Number of secs in one time step                        [s]
!*numsteps*  INTEGER Number of time steps
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: isign, iyear, iyear1, iyear2, msecs, ndays1, ndays2,&
         & nodays, nomsecs, nosecs
INTEGER (KIND=kndilong) :: numsecs
REAL :: daynum1, daynum2
INTEGER, DIMENSION(7) :: idate1, idate2


procname(pglev+1) = 'num_time_steps_int'
CALL log_timer_in()

!
!1. Determine which date comes first
!-----------------------------------
!

IF (intdate1.earlier.intdate2) THEN
   idate1 = intdate1; idate2 = intdate2; isign = 1
ELSEIF (intdate1.later.intdate2) THEN
   idate1 = intdate2; idate2 = intdate1; isign = -1
ELSE
   numsteps = 0
   GOTO 1000
ENDIF

!
!2. Number of time steps
!-----------------------
!

nosecs = dsecs
msecs = MERGE(INT((dsecs-nosecs)*1000),0,dsecs.LT.1000.0)

!---without milliseconds
IF (msecs.EQ.0.AND.idate1(7).EQ.0.AND.idate2(7).EQ.0) THEN
   CALL diff_dates(idate1,idate2,0,numsecs)
   numsteps = isign*numsecs/nosecs

!---with milliseconds
ELSE
!  --number of days between start and end dates
   iyear1 = idate1(1)
   iyear2 = idate2(1)
   CALL day_number(idate1,daynum1)
   CALL day_number(idate2,daynum2)
   ndays1 = daynum1; ndays2 = daynum2
   IF (iyear1.EQ.iyear2) THEN
      nodays = ndays2 - ndays1
   ELSE
      nodays = 365 + leap_year(iyear1) - ndays1
      iyear_110:DO iyear=iyear1+1,iyear2-1
         nodays = nodays + 365 + leap_year(iyear)
      ENDDO iyear_110
      nodays = nodays + ndays2
   ENDIF
!  --number of time steps
   nomsecs = 1000*nosecs+msecs
   numsteps = isign*(1000_kndilong*(86400_kndilong*nodays+&
                   & 3600*(idate2(4)-idate1(4))+&
                   & 60*(idate2(5)-idate1(5))+idate2(6)-idate1(6))+&
                   & idate2(7)-idate1(7))/nomsecs
ENDIF

1000 CALL log_timer_out()


RETURN

END SUBROUTINE num_time_steps_int

!========================================================================

FUNCTION read_clock()
!************************************************************************
!
! *read_clock* Read clock count from process clock
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description -
!
!************************************************************************
!
!* Arguments
!
INTEGER :: read_clock

!
!* Local variables
!
INTEGER :: n


CALL SYSTEM_CLOCK(COUNT=n)
read_clock = n


RETURN

END FUNCTION read_clock

!========================================================================

SUBROUTINE suspend_proc(nosecswait)
!************************************************************************
!
! *suspend_proc* Suspend program execution on calling process by a number of
!                seconds 
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - diff_clock, read_clock
!
!************************************************************************
!
!*Arguments
!
INTEGER, INTENT(IN) :: nosecswait

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*nosecswait* INTEGER Number of seconds to suspend program execution
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: maxcounts, npcc, numcounts


IF (nosecswait.EQ.0) RETURN
npcc = read_clock()
numcounts = 0; maxcounts = nosecswait*npcc_rate
DO WHILE (numcounts.LE.maxcounts)
   numcounts = diff_clock(npcc)
ENDDO

nopcc(itm_wait) = nopcc(itm_wait) + maxcounts


RETURN

END SUBROUTINE suspend_proc

!========================================================================

SUBROUTINE update_time
!************************************************************************
!
! *update_time* Update date/time
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - add_secs_to_date, convert_date, day_number, loop_index
!
!************************************************************************
!
USE paralpars
USE switches
USE utility_routines, ONLY: loop_index

!
!*Local variables
!
LOGICAL :: flag
INTEGER :: icmeteo, icwaves, isign
REAL :: daynum
INTEGER, DIMENSION(7) :: IDateTimeold


procname(pglev+1) = 'update_time'
CALL log_timer_in()

!---current date/time
IDateTimeold = IDateTime
isign = MERGE(1,-1,iopt_part_model.LT.4)
CALL add_secs_to_date(IDateTimeold,IDateTime,isign,delt2d)
CDateTime = convert_date(IDateTime)

!---number of since start
nosecsrun = INT(nt,KIND=kndilong)*INT(delt2d,KIND=kndilong) + &
          & INT(nt,KIND=kndilong)*(delt2d-INT(delt2d))

!---Julian day
CALL day_number(IDateTime,daynum)
julianday = daynum + 1

!---type of time step
IF (iopt_hydro_impl.EQ.0) THEN
   predstep = ((nt-1)/ic3d)*ic3d.EQ.(nt-1)
   corrstep = (nt/ic3d)*ic3d.EQ.nt
ELSE
   predstep = .TRUE.
   corrstep = .TRUE.
ENDIF

!---meteo input time
metstep = .FALSE.
IF (iopt_meteo.EQ.1) THEN
   icmeteo = modfiles(io_metsur,1,1)%tlims(3)
   IF (icmeteo.NE.0) THEN
      IF ((nt/icmeteo)*icmeteo.EQ.nt) metstep = .TRUE.
   ENDIF
ENDIF

!---wave input time
wavestep = .FALSE.
IF (iopt_waves.EQ.1) THEN
   icwaves = modfiles(io_wavsur,1,1)%tlims(3)
   IF (icwaves.NE.0) THEN
      IF ((nt/icwaves)*icwaves.EQ.nt) wavestep = .TRUE.
   ENDIF
ENDIF

!---redefine initial conditions
physinit = modfiles(io_inicon,ics_phys,1)%status.NE.'0'.AND.&
         & loop_index(modfiles(io_inicon,ics_phys,1)%tlims,nt)

!---write current date/time
flag = loglev1.GT.0.AND.(iopt_part_model.NE.2.OR.&
     & (modelid.EQ.modelidcoh.OR.(modelid.EQ.modelidpart.AND.corrstep)))
IF (flag) WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//CDateTime

!---date/time flag
sedlog_date = master

CALL log_timer_out()


RETURN

END SUBROUTINE update_time

!========================================================================

SUBROUTINE write_time_log(iunit,flag)
!************************************************************************
!
! *write_time_log* write date/time to log file if flag is set
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.4
!
! Description -
!
! Reference -
!
!************************************************************************
!

USE timepars

!
!*Arguments
!
INTEGER, INTENT(IN) :: iunit
LOGICAL, INTENT(INOUT) :: flag

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iunit*     INTEGER Unitb of the log file
!*flag*      LOGICAL Write date/time if TRUE. Ser to FALSE on exit
!------------------------------------------------------------------------------
!


IF (flag) THEN
   WRITE (iunit,'(A)') CDateTime
   flag = .FALSE.
ENDIF


RETURN

END SUBROUTINE write_time_log

!========================================================================

SUBROUTINE year_to_date_char(yearnum,chardate)
!************************************************************************
!
! *year_to_date_char* Convert date in absolute years to date/time in character
!                     format
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description -
!
! Module calls - add_secs_to_date, convert_date, leap_year
!
!************************************************************************
!
!* Arguments
!
CHARACTER (LEN=lentime), INTENT(OUT) :: chardate
REAL (KIND=kndrlong), INTENT(IN) :: yearnum

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*yearnum*   REAL    Time in years
!*chardate*  CHAR    Date/time
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ileap, nosecs
INTEGER, DIMENSION(7) :: intdate1, intdate2 


intdate1(1) = INT(yearnum); intdate1(2:7) = 0
ileap = leap_year(intdate1(1))
nosecs = (365+ileap)*86400*(yearnum-INT(yearnum))
CALL add_secs_to_date(intdate1,intdate2,nosecs,1.0)
chardate = convert_date(intdate2)


RETURN

END SUBROUTINE year_to_date_char

!========================================================================

SUBROUTINE year_to_date_int(yearnum,intdate)
!************************************************************************
!
! *year_to_date_int* Convert date in absolute years to date/time in integer
!                    format
!
! Author - Patrick Luyten
!
! Last update - @(COHERENS)time_routines.f90  V2.0
!
! Description -
!
! Module calls - add_secs_to_date, leap_year
!
!************************************************************************
!
!* Arguments
!
REAL (KIND=kndrlong), INTENT(IN) :: yearnum
INTEGER, INTENT(OUT), DIMENSION(7) :: intdate

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*yearnum*   REAL    Time in years
!*intdate*   INTEGER Date/time
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ileap, nosecs
INTEGER, DIMENSION(7) :: intdate1


intdate1(1) = INT(yearnum); intdate1(2:7) = 0
ileap = leap_year(intdate1(1))
nosecs = (365+ileap)*86400*(yearnum-INT(yearnum))
CALL add_secs_to_date(intdate1,intdate,nosecs,1.0)


RETURN

END SUBROUTINE year_to_date_int


END MODULE time_routines
