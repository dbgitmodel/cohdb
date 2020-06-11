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
! Version - @(COHERENS)Usrdef_Output.f90  V2.1.1
!
! $Date: 2018-07-23 10:21:02 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1167 $
!
! Description - output parameters for test case csnsp
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - 
!
! Module calls - close_file, open_file
!
!************************************************************************
!
USE density
USE depths
USE fluxes
USE grid
USE gridpars
USE iopars
USE optics
USE timepars
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: day_number, log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, SAVE :: iunit
INTEGER :: k, kc, klo, kup
REAL :: daynum, dtdz, tbot, tdep, tgrad, tlo, tmean, tsur, tup, twidth, zlo, zup
REAL :: rflag = -999.9

!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---initial conditions
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = modfiles(io_fincon,ics_phys,2)%form
   modfiles(io_inicon,ics_phys,1)%filename = &
 & modfiles(io_fincon,ics_phys,2)%filename
   modfiles(io_fincon,ics_phys,2)%status = '0'
!  ---surface forcing
   modfiles(io_1uvsur,1,1)%status = 'R'
   modfiles(io_1uvsur,1,1)%form = modfiles(io_1uvsur,1,2)%form
   modfiles(io_1uvsur,1,1)%filename = modfiles(io_1uvsur,1,2)%filename
   modfiles(io_1uvsur,1,2)%status = '0'
ENDIF

IF (cold_start) RETURN

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN
!  ---open output file
   CALL open_file(iunit,runtitle(1:6)//'.'//suffix,'OUT','A')
   WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
   WRITE (iunit,'(A)') 'Output parameters for test case csnsp: simulation '//&
                      & TRIM(runtitle)
   WRITE (iunit,*)
ENDIF

!
!3. Store parameters
!-------------------
!
!3.1 Evaluate parameters
!-----------------------
!

IF (nt.GT.0.AND.(mult_index(nt+36,504))) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

!  ---current day   
   CALL day_number(IdateTime,daynum)
   
!  ---surface temperature
   tsur = temp(1,1,nz)

!  ---bottom temperature
   tbot = temp(1,1,1)

!  ---mean temperature
   tmean = SUM(temp(1,1,:))/REAL(nz)

!  ---thermocline depth
   kc = 0
   k_311: DO k=2,nz
      IF ((kc.EQ.0).AND.((temp(1,1,k)-temp(1,1,1)).GT.1.0)) THEN
         kc = k
      ENDIF
   ENDDO k_311
   IF (kc.EQ.0) THEN
      tdep = rflag
   ELSE
      tdep = depmeanatc(1,1)*(1.0-gscoordatc(1,1,kc-1))-delzatw(1,1,kc)*&
           & (temp(1,1,1)+1.0-temp(1,1,kc-1))/(temp(1,1,kc)-temp(1,1,kc-1))
   ENDIF

!  ---maximum temperature gradient
   tgrad = ABS(temp(1,1,nz)-temp(1,1,nz-1))/delzatw(1,1,nz)
   k_312: DO k=nz-1,2,-1
      dtdz = ABS(temp(1,1,k)-temp(1,1,k-1))/delzatw(1,1,k)
      IF (dtdz.GT.tgrad) tgrad = dtdz
   ENDDO k_312

!  ---thermocline width
   klo = 0; kup = 0
   tlo = 8.5; tup = 12.0
   k_313: DO k=1,nz
      IF ((klo.EQ.0).AND.(temp(1,1,k).GT.tlo)) THEN
         klo = k
      ENDIF
      IF ((kup.EQ.0).AND.(temp(1,1,k).GT.tup)) then
         kup = k
      ENDIF
   ENDDO k_313
      
   IF ((klo*kup).GT.0) THEN
      IF (klo.EQ.1) THEN
         zlo = depmeanatc(1,1)*gscoordatc(1,1,klo)
      ELSE
         zlo = depmeanatc(1,1)*gscoordatc(1,1,klo-1)+&
             & (tlo-temp(1,1,klo-1))*delzatw(1,1,klo)&
             & /(temp(1,1,klo)-temp(1,1,klo-1))
      ENDIF
      IF (kup.EQ.1) then
         zup = depmeanatc(1,1)*gscoordatc(1,1,kup)
      ELSE
         zup = depmeanatc(1,1)*gscoordatc(1,1,kup-1)+&
             & (tup-temp(1,1,kup-1))*delzatw(1,1,kup)&
             & /(temp(1,1,kup)-temp(1,1,kup-1))
      ENDIF
      twidth = zup - zlo
   ELSE
      twidth = rflag
   ENDIF

!
!3.2 Write output
!----------------
!

   WRITE (iunit,9002) INT(daynum)
   WRITE (iunit,9001) 'tsur', tsur
   WRITE (iunit,9001) 'tbot', tbot
   WRITE (iunit,9001) 'tmean', tmean
   WRITE (iunit,9001) 'tdep', tdep
   WRITE (iunit,9001) 'tgrad', tgrad
   WRITE (iunit,9001) 'twidth', twidth
   WRITE (iunit,9001) 'qrad', qrad(1,1)
   WRITE (iunit,9001) 'qlwave', qlwave(1,1)

   CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (nt.EQ.nstep) CALL close_file(iunit,'A')


RETURN

9001 FORMAT (T3,A,T12,': ',G12.4E3)
9002 FORMAT ('Time: ',I3,' days')

END SUBROUTINE usrdef_output
