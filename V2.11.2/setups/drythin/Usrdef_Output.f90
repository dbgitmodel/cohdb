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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.8
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - output parameters for test case drythin
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, define_out0d_vals, max_vars, min_vars, mult_index,
!                open_file
!
!************************************************************************
!
USE currents
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE switches
USE syspars
USE timepars
USE paral_utilities, ONLY: max_vars, min_vars
USE inout_routines, ONLY: close_file, open_file
USE model_output, ONLY: define_out0d_vals
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: outflag
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, SAVE :: icout, iunit
REAL :: uvelmax, uvelmin, vvelmax, vvelmin, wphysmax, wphysmin
REAL, DIMENSION(4) :: out0ddat


procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,1)%status = 'R'
   modfiles(io_modgrd,1,1)%form = modfiles(io_modgrd,1,2)%form
   modfiles(io_modgrd,1,1)%filename = modfiles(io_modgrd,1,2)%filename
   modfiles(io_modgrd,1,2)%status = '0'
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1,1)%status = 'R'
   modfiles(io_2uvobc,1,1)%form = modfiles(io_2uvobc,1,2)%form
   modfiles(io_2uvobc,1,1)%filename = modfiles(io_2uvobc,1,2)%filename
   modfiles(io_2uvobc,1,2)%status = '0'
!  ---dry cells
   IF (iopt_drycel.EQ.1) THEN
      modfiles(io_drycel,1,1)%status = 'R'
      modfiles(io_drycel,1,1)%form = modfiles(io_drycel,1,2)%form
      modfiles(io_drycel,1,1)%filename = modfiles(io_drycel,1,2)%filename
      modfiles(io_drycel,1,2)%status = '0'
   ENDIF
!  ---thin dams
   IF (iopt_thndam.EQ.1) THEN
      modfiles(io_thndam,1,1)%status = 'R'
      modfiles(io_thndam,1,1)%form = modfiles(io_thndam,1,2)%form
      modfiles(io_thndam,1,1)%filename = modfiles(io_thndam,1,2)%filename
      modfiles(io_thndam,1,2)%status = '0'
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
      WRITE (iunit,'(A)') 'Output parameters for test case drythin:'//&
                          ' simulation '//TRIM(runtitle)
      WRITE (iunit,*)

   ENDIF
!  ---initialise parameters
   icout = 1440
ENDIF

!
!3. Evaluate parameters
!----------------------
!

outflag = nt.GT.0.AND.mult_index(nt,icout)

IF (outflag) THEN

!
!3.1 Integral parameters
!-----------------------
!

   CALL define_out0d_vals(out0ddat,4,ivarid=(/iarr_ekin0d,iarr_epot0d,&
                        & iarr_edissip0d,iarr_enstr0d/))

!
!3.2 Velocity extrema
!--------------------
!

   CALL max_vars(uvel(1:ncloc,1:nrloc,:),uvelmax,iarr_uvel,commall=.TRUE.,&
                  & mask=node2du(1:ncloc,1:nrloc).GT.1)
   CALL min_vars(uvel(1:ncloc,1:nrloc,:),uvelmin,iarr_uvel,commall=.TRUE.,&
                  & mask=node2du(1:ncloc,1:nrloc).GT.1)
   CALL max_vars(vvel(1:ncloc,1:nrloc,:),vvelmax,iarr_vvel,&
                     & commall=.TRUE.,mask=node2dv(1:ncloc,1:nrloc).GT.1)
   CALL min_vars(vvel(1:ncloc,1:nrloc,:),vvelmin,iarr_vvel,&
                     & commall=.TRUE.,mask=node2dv(1:ncloc,1:nrloc).GT.1)
   CALL max_vars(wphys,wphysmax,iarr_wphys,commall=.TRUE.,mask=maskatc_int)
   CALL min_vars(wphys,wphysmin,iarr_wphys,commall=.TRUE.,mask=maskatc_int)

ENDIF

!
!4. Write output
!---------------
!

IF (master.AND.outflag) THEN
   WRITE (iunit,9001) nosecsrun/3600.0
   WRITE (iunit,9002) 'ekin', 1.0E-09*out0ddat(1)
   WRITE (iunit,9002) 'epot', 1.0E-09*out0ddat(2)
   WRITE (iunit,9002) 'edissip', 0.001*out0ddat(3)
   WRITE (iunit,9002) 'enstr', out0ddat(4)
   WRITE (iunit,9002) 'uvelmax', 100.0*uvelmax
   WRITE (iunit,9002) 'uvelmin', 100.0*uvelmin
   WRITE (iunit,9002) 'vvelmax', 100.0*vvelmax
   WRITE (iunit,9002) 'vvelmin', 100.0*vvelmin
   WRITE (iunit,9002) 'wphysmax', 1000.0*wphysmax
   WRITE (iunit,9002) 'wphysmin', 1000.0*wphysmin
ENDIF

!
!5. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')

1001 CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',F4.1,' hours')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
