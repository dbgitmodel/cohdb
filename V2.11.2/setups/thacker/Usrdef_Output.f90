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
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.5
!
! $Date: 2018-07-23 10:21:02 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1167 $
!
! Description - test case thacker
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - define_out0d_vals
! 
!************************************************************************
!
USE gridpars
USE iopars
USE modids
USE paralpars
USE sedids
USE sedarrays
USE syspars
USE timepars
USE inout_routines, ONLY: close_file, open_file
USE model_output, ONLY: define_out0d_vals
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER :: itime
INTEGER, SAVE :: iunit
REAL :: bflxtot, bstresmean, ceqmean, dryarea, sedmax, sedmin, sedtot
REAL, DIMENSION(5) :: outdat

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
!  ---initial conditions
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = modfiles(io_fincon,ics_phys,2)%form
   modfiles(io_inicon,ics_phys,1)%filename = modfiles(io_fincon,1,2)%filename
   modfiles(io_fincon,ics_phys,2)%status = '0'
!  ---sediment particle attributes
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
!  ---open output file
   IF (master) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case thacker: '&
                      & //'simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF
ENDIF

!
!3. Store parameters
!-------------------
!

IF (mult_index(nt,300)) THEN

   procname(pglev+1) = 'usrdef_output'
   CALL log_timer_in()

   itime = nosecsrun

!  ---max, min, total concentrations, bottom flux
   CALL define_out0d_vals(outdat,5,&
             & ivarid=(/iarr_bstresatc,iarr_ctot,iarr_ctot,iarr_ctot,&
             & iarr_bottom_sed_flux/),&
             & numvar=(/0,0,0,0,1/),&
             & oopt=(/oopt_mean,oopt_max,oopt_min,oopt_int,oopt_int/))
   bstresmean = outdat(1)
   sedmax = outdat(2); sedmin = outdat(3)
   sedtot = rhos(1)*outdat(4); bflxtot = rhos(1)*outdat(5)

!  ---dry area
   dryarea = 1.0-nowetatc/REAL(noseaatc)

!  ---write output
   IF (master) THEN
      WRITE (iunit,9001) itime
      WRITE (iunit,9002) 'bstresmean', bstresmean
      WRITE (iunit,9002) 'sedmax', sedmax
      WRITE (iunit,9002) 'sedmin', sedmin
      WRITE (iunit,9002) 'sedtot', sedtot
      WRITE (iunit,9002) 'bflxtot', bflxtot
      WRITE (iunit,9002) 'dryarea', dryarea
   ENDIF

   CALL log_timer_out()

ENDIF

!
!4. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')


RETURN

9001 FORMAT ('Time: ',I4,' seconds')
9002 FORMAT (T3,A,T15,': ',G12.4E3)

END SUBROUTINE usrdef_output
