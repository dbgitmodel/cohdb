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
! Version - @(COHERENS)Usrdef_Output.f90  V2.4
!
! $Date: 2018-07-23 10:21:02 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1167 $
!
! Description - output parameters for test case fredy
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - 
!
! Module calls - close_file, energy_0d, enstrophy_0d, max_vars, min_vars,
!                mult_index, open_file, sum_vars, sum2_vars
!
!************************************************************************
!
USE density
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE syspars
USE timepars
USE diagnostic_routines, ONLY: energy_0d, enstrophy_0d
USE inout_routines, ONLY: close_file, open_file
USE paral_utilities, ONLY: max_vars, min_vars, sum_vars, sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, SAVE :: iunit, knt, mtot, ntday
INTEGER :: iday, i1pt, i1ptloc, k
REAL, SAVE :: stotin, sumcos, sumres
REAL :: a1pt, ekin, enstr, epot, salint, salmax, salmin, scrit, sdev, stot, &
      & sum1, sum2, theta, vdens
INTEGER, DIMENSION(4) :: nhdims = 0
REAL, DIMENSION(5) :: xvec
REAL, DIMENSION(ncloc,nrloc,nz) :: array3dc


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
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN
!  ---open output file
   IF (.NOT.cold_start.AND.master) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case fredy: simulation '&
                        & //TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF
!  ---initial value of volume averaged salinity
   CALL sum2_vars(delzatc(1:ncloc,1:nrloc,:),vdens,nhdims,'C  ',0)
   array3dc = sal(1:ncloc,1:nrloc,:)*delzatc(1:ncloc,1:nrloc,:)
   CALL sum2_vars(array3dc,salint,nhdims,'C  ',0)
   IF (master) stotin = salint/vdens
!  ---initialise parameters
   sumres = 0.0
   sumcos = 0.0
   mtot = 288
   knt = -mtot
   ntday = 2880
ENDIF

IF (.NOT.(nt.EQ.0.OR.corrstep)) RETURN

CALL log_timer_in()
procname(pglev+1) = 'usrdef_output'
 
!
!3. Evaluate parameters
!----------------------
!
!---day number
iday = nt/ntday

!---energy components
CALL energy_0d(3,xvec)
IF (master) THEN
   ekin = 0.001*xvec(1) 
   epot = 0.001*(xvec(2)+xvec(3)+xvec(4))
ENDIF

!---enstrophy 
enstr =  enstrophy_0d(3)
IF (master) enstr = 1.0E+06*enstr

!---area bounded by surface 34.839 psu contour
scrit = 34.839
i1ptloc = COUNT(sal(1:ncloc,1:nrloc,nz).LT.scrit.AND.maskatc_int)
CALL sum_vars(i1ptloc,i1pt,0)
IF (master) a1pt = 0.01*i1pt

!---minimum and maximum salinity
CALL min_vars(sal(1:ncloc,1:nrloc,:),salmin,iarr_sal,mask=maskatc_int)
CALL max_vars(sal(1:ncloc,1:nrloc,:),salmax,iarr_sal,mask=maskatc_int)

!---salinity deviation
CALL sum2_vars(delzatc(1:ncloc,1:nrloc,:),vdens,nhdims,'C  ',0)
array3dc = sal(1:ncloc,1:nrloc,:)*delzatc(1:ncloc,1:nrloc,:)
CALL sum2_vars(array3dc,salint,nhdims,'C  ',0)
IF (master) THEN
   stot = salint/vdens
   sdev = 1.0E+06*(stot-stotin)/stotin
ENDIF

!
!4. Harmonic analysis
!--------------------
!

IF (master.AND.iday.GE.4) THEN
   sumres = sumres + ekin/epot
   sumcos = sumcos + COS(coriolatu(1,1)*knt*delt3d)*ekin/epot
   knt = knt + 1
ENDIF

!
!5. Write output
!---------------
!

IF (master.AND.(nt.GT.0).AND.(mult_index(nt,ntday))) THEN
   WRITE (iunit,9001) iday
   WRITE (iunit,9002) 'ekin', ekin
   WRITE (iunit,9002) 'epot', epot
   WRITE (iunit,9002) 'enstr', enstr
   WRITE (iunit,9002) 'a1pt', a1pt
   WRITE (iunit,9002) 'salmin', salmin
   WRITE (iunit,9002) 'salmax', salmax
   WRITE (iunit,9002) 'sdev', sdev
ENDIF

!
!6. Write harmonic output
!------------------------
!

IF (master.AND.nt.EQ.nstep) THEN

   sum1 = SIN((mtot+0.5)*coriolatu(1,1)*delt3d)/SIN(0.5*coriolatu(1,1)*delt3d)
   sum2 = 0.5*SIN((2*mtot+1)*coriolatu(1,1)*delt3d)/&
        & SIN(coriolatu(1,1)*delt3d) + mtot + 0.5
   theta = (sumres*sum2-sum1*sumcos)/((2*mtot+1)*sum2-sum1*sum1)       

   WRITE (iunit,*)
   WRITE (iunit,'(A)') 'Harmonic parameters'
   WRITE (iunit,9002) 'theta', theta

ENDIF

!
!7. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')

CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',I1,' days')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
