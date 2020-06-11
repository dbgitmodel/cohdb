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
! Description - output parameters for test case convtest
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, max_vars, min_vars, mult_index,
!                open_file
!
!************************************************************************
!
USE currents  
USE density
USE depths
USE grid  
USE gridpars
USE iopars
USE modids
USE paralpars
use physpars
USE syspars
USE timepars
USE inout_routines, ONLY: close_file, open_file
USE paral_comms, ONLY: combine_mod
USE paral_utilities, ONLY: max_vars, min_vars 
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=12) :: cl
INTEGER, PARAMETER :: nostats = 9
INTEGER, SAVE :: icout, iunit
INTEGER :: i, k, l
INTEGER, SAVE, DIMENSION(nostats) :: ipos
REAL :: betaw, uvelmax, uvelmin, wphysmax, wphysmin
REAL, DIMENSION(nostats) :: fbuo2mean, fbuo2max, fbuo2min
REAL, DIMENSION(2:nz,nostats) :: bgrad
REAL, DIMENSION(0:nc+1,0:nr+1) :: gaccatcglb
REAL, DIMENSION(0:nc+1,0:nr+1,nz) :: densglb
REAL, DIMENSION(0:nc+1,0:nr+1,nz+1) :: delzatwglb


!
!1. Reset forcing attributes for forcing files
!---------------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---initial conditions
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = modfiles(io_fincon,ics_phys,2)%form
   modfiles(io_inicon,ics_phys,1)%filename = &
 & modfiles(io_fincon,ics_phys,2)%filename
   modfiles(io_inicon,ics_phys,2)%status = '0'
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN

!  ---open output file
   IF (master) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case convtest: '//&
                        & 'simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF

!  ---output frequency
   icout = 2880

!   ---station locations
   ipos = (/5,10,15,20,25,30,35,40,45/)

ENDIF

!
!3. Store parameters
!-------------------
!

IF (mult_index(nt,icout,matchzero=.TRUE.)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

!
!3.1 Combine
!-----------
!

   CALL combine_mod(densglb,dens,(/0,0,1/),iarr_dens,0.0,commall=.TRUE.)
   CALL combine_mod(gaccatcglb,gaccatc,(/0,0/),iarr_gaccatc,0.0,commall=.TRUE.)
   CALL combine_mod(delzatwglb,delzatw,(/0,0,1/),iarr_delzatw,0.0,commall=.TRUE.)

!   
!3.2 Current data 
!----------------
!
   
   CALL max_vars(uvel(1:ncloc,1:nrloc,:),uvelmax,iarr_uvel,commall=.TRUE.,&
               & mask=node2du(1:ncloc,1:nrloc).GT.1)
   CALL min_vars(uvel(1:ncloc,1:nrloc,:),uvelmin,iarr_uvel,commall=.TRUE.,&
               & mask=node2du(1:ncloc,1:nrloc).GT.1)
   CALL max_vars(wphys,wphysmax,iarr_wphys,commall=.TRUE.,mask=maskatc_int)
   CALL min_vars(wphys,wphysmin,iarr_wphys,commall=.TRUE.,mask=maskatc_int)
   
!
!3.3 Station data
!----------------
!   

   l_330: DO l=1,nostats
      i = ipos(l)
      k_331: DO k=2,nz
         bgrad(k,l)= -gaccatcglb(i,1)*(densglb(i,1,k)-densglb(i,1,k-1))/&
                    & delzatwglb(i,1,k)/density_ref
      ENDDO k_331
      fbuo2mean(l) = SUM(bgrad(:,l))/REAL(nz-1)
      fbuo2max(l) = MAXVAL(bgrad(:,l))
      fbuo2min(l) = MINVAL(bgrad(:,l))
   ENDDO l_330

!
!3.4 Write output
!----------------
!

   IF (master) THEN
      
      WRITE (iunit,9001) nosecsrun/86400.0
      WRITE (iunit,*)
      WRITE (iunit,'(A)') 'Global parameters'
      WRITE (iunit,9003) 'uvelmax', uvelmax
      WRITE (iunit,9003) 'uvelmin', uvelmin
      WRITE (iunit,9003) 'wphysmax', wphysmax
      WRITE (iunit,9003) 'wphysmin', wphysmin
      WRITE (iunit,*)
   
      l_340: DO l=1,nostats
         WRITE (iunit,9002) 'Station', l
         WRITE (cl,'(I12)') l; cl = ADJUSTL(cl)
         WRITE (iunit,9003) 'fbuo2mean'//TRIM(cl), fbuo2mean(l)
         WRITE (iunit,9003) 'fbuo2max'//TRIM(cl), fbuo2max(l)
         WRITE (iunit,9003) 'fbuo2min'//TRIM(cl), fbuo2min(l)
      ENDDO l_340
      
      WRITE (iunit,*)
      
   ENDIF

   CALL log_timer_out()
   
ENDIF

!
!4. Close output file
!--------------------
!

IF (master.AND.(cold_start.OR.nt.EQ.nstep)) CALL close_file(iunit,'A')


RETURN

9001 FORMAT ('Time: ',F5.1,' days')
9002 FORMAT (A,': ',I2)
9003 FORMAT (T3,A,T14,': ',G12.4E3)

END SUBROUTINE usrdef_output
