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

MODULE timepars
!************************************************************************
!
! *timepars* Time parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)timepars.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

LOGICAL :: corrstep, metstep, predstep, physinit, wavestep
CHARACTER (LEN=lentime) :: CDateTime, CEndDateTime = cdatetime_undef,&
                         & CStartDateTime = cdatetime_undef, ClockTime
INTEGER :: julianday, norestarts = 0, nstep = 0, nt = 0, ntobcrlx = 0
INTEGER :: iccvt = 0, ic3d = 1, icnodal = 0
INTEGER (KIND=kndilong) :: nosecsrun = 0
REAL :: delt2d, delt3d, hydro_end = 0, time_zone = 0.0
INTEGER, DIMENSION(7) :: time_convert = (/1,60,3600,86400,604800,2629800,&
                                        & 31557600/)
INTEGER, DIMENSION(7) :: IDateTime, IEndDateTime, IStartDateTime
INTEGER, DIMENSION(12) :: days_in_month = &
                        & (/31,28,31,30,31,30,31,31,30,31,30,31/) 
INTEGER, DIMENSION(13) :: monthdays = &
                        & (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
INTEGER, DIMENSION(MaxRestarts) :: ntrestart = 0

SAVE

!
! Name           Type     Purpose
!------------------------------------------------------------------------------
!*CDateTime*     CHAR     Current date/time
!*CEndDateTime*  CHAR     End date/time
!*ClockTime*     CHAR     Real date/time at start of simulation
!*CStartDateTime*CHAR     Start date/time
!*corrstep*      LOGICAL  .TRUE. at corrector time steps
!*days_in_month* INTEGER  Number of days in each month
!*delt2d*        REAL     Time step for 2-D mode                           [s]
!*delt3d*        REAL     Time step for 3-D mode and scalars               [s]
!*hydro_end*     INTEGER  Time step at the end of the last hydrodynamical cycle
!                         (used in case of tidal acceleration)
!*iccvt*         INTER    Counter for convective adjustment
!*icnodal*       INTEGER  Counter for update of nodal factors and phases
!*ic3d*          INTEGER  Counter for 3-D mode calculations
!*IDateTime*     INTEGER  Current date/time
!*IEndDateTime*  INTEGER  End date/time
!*IStartDateTime*INTEGER  Start date/time
!*julianday*     INTEGER  Julian day (year day between 1 and 365/366)
!*metstep*       LOGICAL  .TRUE. for meteo input time steps
!*monthdays*     INTEGER  Day number of first day in each month
!*norestarts*    INTEGER  Number of restarts
!*nosecsrun*     LONGINT  Number of seconds since start of simulation
!*nstep*         INTEGER  Number of 2-D time steps for simulation
!*nt*            INTEGER  (2-D) time step counter
!*ntobcrlx*      INTEGER  Period for relaxation of 2-D open boundary condition
!                         from zero initial state, measured in 2-D time steps
!                         (mode splitting scheme only)
!*ntrestart*     INTEGER   Time indices for writing to restart file
!                         zone                                         [hours]
!*physinit*      LOGICAL  .TRUE. for (re)initialisation of physical conditions
!*predstep*      LOGICAL  .TRUE. at predictor time steps
!*time_convert*  INTEGER  Time conversion units                            [s]
!*time_zone*     REAL     Time difference between local and Greenwich time 
!*wavestep*      LOGICAL  .TRUE. for wave input time steps
!
!************************************************************************
!
!integer Date/Time format : (year,month,day,hour,min,sec,millisec)
!character Date/Time format :yyyy/mm/dd;hh:mm:ss,xxx


END MODULE timepars
