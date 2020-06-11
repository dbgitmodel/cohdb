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
! Version - @(COHERENS)Usrdef_Output.f90  V2.5
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - test case flicest
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, define_out0d_vals, error_alloc,
!                open_file, mult_index
! 
!************************************************************************
!
USE gridpars
USE iopars
USE paralpars
USE sedids
USE sedarrays
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_file, open_file
USE model_output, ONLY: define_out0d_vals
USE paral_comms, ONLY: combine_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=12) :: cl
INTEGER, PARAMETER :: nostats = 9
INTEGER, SAVE :: icount, iunit
INTEGER :: i, j, k, l
INTEGER, SAVE, DIMENSION(nostats) :: istat, kstat
INTEGER, DIMENSION(4) :: lbounds
REAL, DIMENSION(3) :: outdat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: floc_densglb, floc_diaglb, &
                                           & floc_ncglb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: cvolglb, wfallglb 

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
   modfiles(io_2uvobc,1:2,1)%status = 'R'
   modfiles(io_2uvobc,1:2,1)%form = modfiles(io_2uvobc,1:2,2)%form
   modfiles(io_2uvobc,1,1)%filename = modfiles(io_2uvobc,1,2)%filename
   modfiles(io_2uvobc,2,1)%filename = modfiles(io_2uvobc,2,2)%filename
   modfiles(io_2uvobc,1:2,2)%status = '0'
!  ---open boundary conditions (3-D)
   modfiles(io_3uvobc,1:2,1)%status = 'R'
   modfiles(io_3uvobc,1:2,1)%form = modfiles(io_3uvobc,1:2,2)%form
   modfiles(io_3uvobc,1,1)%filename = modfiles(io_3uvobc,1,2)%filename
   modfiles(io_3uvobc,2,1)%filename = modfiles(io_3uvobc,2,2)%filename
   modfiles(io_3uvobc,1:2,2)%status = '0'
 !  ---open boundary conditions (salinity)
   modfiles(io_salobc,1:2,1)%status = 'R'
   modfiles(io_salobc,1:2,1)%form = modfiles(io_salobc,1:2,2)%form
   modfiles(io_salobc,1,1)%filename = modfiles(io_salobc,1,2)%filename
   modfiles(io_salobc,2,1)%filename = modfiles(io_salobc,2,2)%filename
   modfiles(io_salobc,1:2,2)%status = '0'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,1)%status = 'R'
   modfiles(io_sedspc,1,1)%form = modfiles(io_sedspc,1,2)%form
   modfiles(io_sedspc,1,1)%filename = modfiles(io_sedspc,1,2)%filename
   modfiles(io_sedspc,1,2)%status = '0'
!  ---sediment initial conditions
   IF (runtitle(8:8).EQ.'B') THEN    
      modfiles(io_inicon,ics_sed,1)%status = 'R'
      modfiles(io_inicon,ics_sed,1)%form = modfiles(io_fincon,ics_sed,2)%form
      modfiles(io_inicon,ics_sed,1)%filename = &
    & modfiles(io_fincon,ics_sed,2)%filename
      modfiles(io_fincon,ics_sed,2)%status = '0'
   ENDIF

ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (cold_start.OR.runtitle(8:8).EQ.'0') RETURN

IF (nt.EQ.0) THEN

!  ---open output file
   IF (master) THEN
      CALL open_file(iunit,TRIM(runtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case flocest: '&
                      & //'simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF

!  ---station indices
   istat = (/25,25,25,50,50,50,75,75,75/)
   kstat = (/1,10,20,1,10,20,1,10,20/)

!  ---allocate arrays
   ALLOCATE (cvolglb(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz,3),STAT=errstat)
   CALL error_alloc('cvolglb',4,(/nc+2*nhalo,nr+2*nhalo,nz,3/),kndrtype)
   ALLOCATE (floc_densglb(0:nc+1,0:nr+1,nz),STAT=errstat)
   CALL error_alloc('floc_densglb',3,(/nc+2,nr+2,nz/),kndrtype)
   ALLOCATE (floc_diaglb(nc,nr,nz),STAT=errstat)
   CALL error_alloc('floc_diaglb',3,(/nc,nr,nz/),kndrtype)
   ALLOCATE (floc_ncglb(0:nc+1,0:nr+1,nz),STAT=errstat)
   CALL error_alloc('floc_ncglb',3,(/nc+2,nr+2,nz/),kndrtype)
   ALLOCATE (wfallglb(0:nc,0:nr,nz+1,3),STAT=errstat)
   CALL error_alloc('wfallglb',4,(/nc+1,nr+1,nz+1,3/),kndrtype)
   
!  ---output frequency   
   icount = 2400
   
ENDIF

IF (.NOT.mult_index(nt,icount)) RETURN

!
!3. Store parameters
!-------------------
!

procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

!
!3.1 Domain integrated masses
!----------------------------
!

CALL define_out0d_vals(outdat(1:2),2,ivarid=(/iarr_floc_P,iarr_floc_F/),&
                     & oopt=(/oopt_int,oopt_int/))
outdat(3) = outdat(1) + outdat(2)
outdat = 0.001*outdat

!
!3.2 Combine into global arrays
!------------------------------
!   

lbounds = (/1-nhalo,1-nhalo,1,1/)
CALL combine_mod(cvolglb,cvol,lbounds,iarr_cvol,0.0,commall=.TRUE.)
lbounds(1:3) = (/0,0,1/)
CALL combine_mod(floc_densglb,floc_dens,lbounds(1:3),iarr_floc_dia,0.0,&
               & commall=.TRUE.)
lbounds(1:3) = 1
CALL combine_mod(floc_diaglb,floc_dia,lbounds(1:3),iarr_floc_dia,0.0,&
               & commall=.TRUE.)
lbounds(1:3) = (/0,0,1/)
CALL combine_mod(floc_ncglb,floc_nc,lbounds(1:3),iarr_floc_dia,0.0,&
               & commall=.TRUE.)
lbounds = (/0,0,1,1/)
CALL combine_mod(wfallglb,wfall,lbounds,iarr_wfall,0.0,commall=.TRUE.)

!
!4. Write parameters
!-------------------
!

IF (master) THEN

!  ---domain integrated parameters    
   WRITE (iunit,9001) nosecsrun/3600.0
   WRITE (iunit,9002) 'Pintmass', outdat(1)
   WRITE (iunit,9002) 'Fintmass', outdat(2)
   WRITE (iunit,9002) 'Sedintmass', outdat(3)

!  ---station data
   j = 5
   l_410: DO l=1,nostats
      i = istat(l); k = kstat(l)
      WRITE (iunit,9003) 'Station', l
      WRITE (cl,'(I12)') l; cl = ADJUSTL(cl)
      WRITE (iunit,9002) 'Pdens_'//TRIM(cl), rhos(1)*cvolglb(i,j,k,1)
      WRITE (iunit,9002) 'Fdens_'//TRIM(cl), floc_densglb(i,j,k)*&
                        & cvolglb(i,j,k,2)
      WRITE (iunit,9002) 'floc_dens_'//TRIM(cl), floc_densglb(i,j,k)
      WRITE (iunit,9002) 'floc_dia'//TRIM(cl), 1.0E+06*floc_diaglb(i,j,k)
      WRITE (iunit,9002) 'floc_nc_'//TRIM(cl), floc_ncglb(i,j,k)
      WRITE (iunit,9002) 'Fwfall_'//TRIM(cl), 100.0*wfallglb(i,j,k,2)
   ENDDO l_410
      
ENDIF

CALL log_timer_out()

!
!5. Finalise
!-----------
!

IF (nt.EQ.nstep) THEN

!  ---close output file   
   IF (master) CALL close_file(iunit,'A')

!  ---deallocate
   DEALLOCATE (cvolglb,floc_densglb,floc_diaglb,floc_ncglb,wfallglb)

ENDIF


RETURN

9001 FORMAT ('Time: ',F5.2,' hours')
9002 FORMAT (T3,A,T15,': ',G12.4E3)
9003 FORMAT (A,': ',I3)

END SUBROUTINE usrdef_output
