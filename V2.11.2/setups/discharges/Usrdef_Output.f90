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

SUBROUTINE usrdef_output
!************************************************************************
!
! *usrdef_output* User-formatted output
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.6
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - output parameters for test case discharges
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, complex_polar, energy_0d, error_alloc,
!                mult_index, open_file, sum2_vars
!
!************************************************************************
!
USE currents
USE depths
USE gridpars
USE iopars
USE modids
USE paralpars
USE physpars
USE structures
USE syspars
USE timepars
USE diagnostic_routines, ONLY: energy_0d
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: complex_polar
USE paral_comms, ONLY: combine_mod
USE paral_utilities, ONLY: sum2_vars
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: outflag
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=12) :: cn
INTEGER, SAVE :: icout, iunit
INTEGER :: i, j, n
INTEGER, SAVE, DIMENSION(4) :: nhdims
REAL, SAVE :: balance, delxdat, delydat, distot, dvolume, volume_ini
REAL :: ekin, epot
REAL, DIMENSION(5) :: xvec
REAL, DIMENSION(numdis) :: veldird, velmagd, zetad
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: umvelglb, vmvelglb, zetaglb


procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN

!  ---discharge specifiers
   modfiles(io_disspc,1,1)%status = 'R'
   modfiles(io_disspc,1,1)%form = modfiles(io_disspc,1,2)%form
   modfiles(io_disspc,1,1)%filename = modfiles(io_disspc,1,2)%filename
   modfiles(io_disspc,1,2)%status = '0'
!  ---discharge locations
   IF (LLT(runtitle(11:11),'C')) THEN
      modfiles(io_disloc,1,1)%status = 'R'
      modfiles(io_disloc,1,1)%form = modfiles(io_disloc,1,2)%form
      modfiles(io_disloc,1,1)%filename = modfiles(io_disloc,1,2)%filename
      modfiles(io_disloc,1,2)%status = '0'
   ENDIF
!  ---discharge rates
   IF (LLT(runtitle(11:11),'D')) THEN
      modfiles(io_disvol,1,1)%status = 'R'
      modfiles(io_disvol,1,1)%form = modfiles(io_disvol,1,2)%form
      modfiles(io_disvol,1,1)%filename = modfiles(io_disvol,1,2)%filename
      modfiles(io_disvol,1,2)%status = '0'
   ENDIF
!  ---momentum discharge
   IF (LGT(runtitle(11:11),'A')) THEN
      modfiles(io_discur,1,1)%status = 'R'
      modfiles(io_discur,1,1)%form = modfiles(io_discur,1,2)%form
      modfiles(io_discur,1,1)%filename = modfiles(io_discur,1,2)%filename
      modfiles(io_discur,1,2)%status = '0'
   ENDIF
!  ---salinity discharge
   IF (LLT(runtitle(11:11),'D')) THEN
      modfiles(io_dissal,1,1)%status = 'R'
      modfiles(io_dissal,1,1)%form = modfiles(io_dissal,1,2)%form
      modfiles(io_dissal,1,1)%filename = modfiles(io_dissal,1,2)%filename
      modfiles(io_dissal,1,2)%status = '0'
   ENDIF

   GOTO 1001

ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN
!  ---open output file
   IF (.NOT.cold_start.AND.master) THEN
      CALL open_file(iunit,TRIM(outtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case dischg:'//&
                          ' simulation '//TRIM(runtitle)
      WRITE (iunit,*)

   ENDIF
!  ---initialise parameters
   icout = 12
   nhdims = 1
   delxdat = surfacegrids(igrd_model,1)%delxdat
   delydat = surfacegrids(igrd_model,1)%delydat
   volume_ini = (nc-1)*(nr-1)*depmean_cst
   distot = 0.0
ENDIF

!
!3. Evaluate parameters
!----------------------
!
!3.1 Update total volume increase and amount of discharge
!--------------------------------------------------------
!

CALL sum2_vars(zeta,dvolume,nhdims,'C  ',iarr_zeta,commall=.TRUE.)
IF (nt.GT.0) THEN
   distot = distot + delt2d*SUM(disvol)/(delxdat*delydat)
ENDIF

!
!3.2 Evaluate output parameters
!------------------------------
!
 
outflag = nt.GT.0.AND.mult_index(nt,icout)

IF (outflag) THEN

!  ---volume increase and balance
   dvolume = 100.0*dvolume/volume_ini
   balance = 100.0*(dvolume-distot)/volume_ini

!  ---energy components
   CALL energy_0d(3,xvec)
   ekin = 0.001*xvec(1) 
   epot = 0.001*(xvec(2)+xvec(3)+xvec(4))

!  ---allocate
   ALLOCATE(zetaglb(0:nc+1,0:nr+1),STAT=errstat)
   CALL error_alloc('zetaglb',2,(/nc+2,nr+2/),real_type)
   ALLOCATE(umvelglb(1-nhalo:nc+nhalo,0:nr+1),STAT=errstat)
   CALL error_alloc('umvelglb',2,(/nc+2*nhalo,nr+2/),real_type)
   ALLOCATE(vmvelglb(0:nc+1,1-nhalo:nr+nhalo),STAT=errstat)
   CALL error_alloc('vmvelglb',2,(/nc+2,nr+2*nhalo/),real_type)

!  ---combine local arrays
   CALL combine_mod(zetaglb,zeta,(/0,0/),iarr_zeta,0.0,commall=.TRUE.)
   CALL combine_mod(umvelglb,umvel,(/1-nhalo,0/),iarr_umvel,0.0,commall=.TRUE.)
   CALL combine_mod(vmvelglb,vmvel,(/0,1-nhalo/),iarr_vmvel,0.0,commall=.TRUE.)

!  ---parameters at discharge locations
   n_320: DO n=1,numdis
      i = idis(n); j = jdis(n)
      zetad(n) = zetaglb(i,j)
      CALL complex_polar(umvelglb(i,j),vmvelglb(i,j),xamp=velmagd(n),&
                       & xpha=veldird(n))
   ENDDO n_320

!  ---deallocate
   DEALLOCATE (zetaglb,umvelglb,vmvelglb)

ENDIF

!
!4. Write output parameters
!--------------------------
!

IF (master.AND.outflag) THEN

   WRITE (iunit,9001) nosecsrun
   WRITE (iunit,9002) 'dvolume', dvolume
   WRITE (iunit,9002) 'balance', balance
   WRITE (iunit,9002) 'ekin', ekin
   WRITE (iunit,9002) 'epot', epot
   n_410: DO n=1,numdis
      WRITE (cn,'(I12)') n; cn = ADJUSTL(cn)
      WRITE (iunit,9002) 'zetad'//TRIM(cn), 100.0*zetad(n)
      WRITE (iunit,9002) 'velmagd'//TRIM(cn), 100.0*velmagd(n)
      WRITE (iunit,9002) 'veldird'//TRIM(cn), radtodeg*veldird(n)
   ENDDO n_410

ENDIF

!
!5. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')

1001 CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',I3,' seconds')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
