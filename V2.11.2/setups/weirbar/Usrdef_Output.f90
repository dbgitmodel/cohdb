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
! Version - @(COHERENS)Usrdef_Output.f90  V2.6
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - output parameters for test case weirbar
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_mod, combine_stats_glb, error_alloc,
!                max_vars, min_vars, mult_index, open_file
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE structures
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: combine_mod, combine_stats_glb
USE paral_utilities, ONLY: max_vars, min_vars
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: outflag
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
CHARACTER (LEN=12) :: cl
INTEGER :: i, ii, l
INTEGER, SAVE :: icout, iunit, numweirs
REAL :: uvelmax, uvelmin, vvelmax, vvelmin, wphysmax, wphysmin
REAL, DIMENSION(7) :: eloss, wdeptot, wmean
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: deptotatuglb, umvelglb


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
   modfiles(io_2uvobc,1,1)%form =modfiles(io_2uvobc,1,2)%form
   modfiles(io_2uvobc,1,1)%filename = modfiles(io_2uvobc,1,2)%filename
   modfiles(io_2uvobc,1,2)%status = '0'
!  ---structures
   modfiles(io_weibar,1,1)%status = 'R'
   modfiles(io_weibar,1,1)%form = modfiles(io_weibar,1,2)%form
   modfiles(io_weibar,1,1)%filename = modfiles(io_weibar,1,2)%filename
   modfiles(io_weibar,1,2)%status = '0'
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
      WRITE (iunit,'(A)') 'Output parameters for test case weirbar:'//&
                          ' simulation '//TRIM(runtitle)
      WRITE (iunit,*)
   ENDIF
!  ---initialise parameters
   numweirs = MERGE(1,7,numwbaru.EQ.1)
   icout = 1440
ENDIF

!
!3. Evaluate parameters
!----------------------
!

outflag = nt.GT.0.AND.mult_index(nt,icout)

IF (outflag) THEN

!
!3.1 Allocate
!------------
!

   ALLOCATE(deptotatuglb(1-nhalo:nc+nhalo,0:nr+1),STAT=errstat)
   CALL error_alloc('deptotatuglb',2,(/nc+2*nhalo,nr+2/),real_type)
   ALLOCATE(umvelglb(1-nhalo:nc+nhalo,0:nr+1),STAT=errstat)
   CALL error_alloc('umvelglb',2,(/nc+2*nhalo,nr+2/),real_type)

!
!3.2 Combine local arrays
!------------------------
!

   CALL combine_mod(deptotatuglb,deptotatu,(/1-nhalo,0/),iarr_deptotatu,0.0)
   CALL combine_mod(umvelglb,umvel,(/1-nhalo,0/),iarr_umvel,0.0)
   IF (iopt_MPI.EQ.1) THEN
      CALL combine_stats_glb(wbarelossu,numwbaru,nowbaruprocs,indexwbaruprocs,&
                           & iarr_wbarelossu,commall=.TRUE.)
   ENDIF

!
!3.3 Parameters at weirs
!-----------------------
!

   IF (numweirs.EQ.1) THEN
      i = iwbaru(1)
      wdeptot(1) = deptotatuglb(i,1)
      wmean(1) = umvelglb(i,1)
      eloss(1) = wbarelossu(1)
   ELSE
      l_310: DO l=1,7
         ii = 5*l-4; i = iwbaru(ii)
         wdeptot(l) = SUM(deptotatuglb(i,1:5))/5.0
         wmean(l) = SUM(umvelglb(i,1:5))/5.0
         eloss(l) = SUM(wbarelossu(ii:ii+4))/5.0
      ENDDO l_310
   ENDIF

!
!3.4 Velocity extrema
!--------------------
!

   IF (iopt_grid_nodim.EQ.3) THEN
      CALL max_vars(uvel(1:ncloc,1:nrloc,:),uvelmax,iarr_uvel,commall=.TRUE.,&
                  & mask=node2du(1:ncloc,1:nrloc).GT.1)
      CALL min_vars(uvel(1:ncloc,1:nrloc,:),uvelmin,iarr_uvel,commall=.TRUE.,&
                  & mask=node2du(1:ncloc,1:nrloc).GT.1)
      IF (nr.GT.2) THEN
         CALL max_vars(vvel(1:ncloc,1:nrloc,:),vvelmax,iarr_vvel,&
                     & commall=.TRUE.,mask=node2dv(1:ncloc,1:nrloc).GT.1)
         CALL min_vars(vvel(1:ncloc,1:nrloc,:),vvelmin,iarr_vvel,&
                     & commall=.TRUE.,mask=node2dv(1:ncloc,1:nrloc).GT.1)
      ENDIF
      CALL max_vars(wphys,wphysmax,iarr_wphys,commall=.TRUE.,mask=maskatc_int)
      CALL min_vars(wphys,wphysmin,iarr_wphys,commall=.TRUE.,mask=maskatc_int)
   ELSE
      CALL max_vars(umvel(1:ncloc,1:nrloc),uvelmax,iarr_umvel,commall=.TRUE.,&
                  & mask=node2du(1:ncloc,1:nrloc).GT.1)
      CALL min_vars(umvel(1:ncloc,1:nrloc),uvelmin,iarr_umvel,commall=.TRUE.,&
                  & mask=node2du(1:ncloc,1:nrloc).GT.1)
      IF (nr.GT.2) THEN
         CALL max_vars(vmvel(1:ncloc,1:nrloc),vvelmax,iarr_vmvel,&
                     & commall=.TRUE.,mask=node2dv(1:ncloc,1:nrloc).GT.1)
         CALL min_vars(vmvel(1:ncloc,1:nrloc),vvelmin,iarr_vmvel,&
                     & commall=.TRUE.,mask=node2dv(1:ncloc,1:nrloc).GT.1)
      ENDIF
   ENDIF

!
!3.5 Deallocate
!--------------
!

   DEALLOCATE (deptotatuglb,umvelglb)

ENDIF

!
!4. Write output
!---------------
!

IF (master.AND.outflag) THEN
   WRITE (iunit,9001) nosecsrun/3600.0
   l_410: DO l=1,numweirs
      WRITE (cl,'(I12)') l; cl = ADJUSTL(cl)
      WRITE (iunit,9002) 'wdeptot'//TRIM(cl), wdeptot(l)
      WRITE (iunit,9002) 'wmean'//TRIM(cl), wmean(l)
      WRITE (iunit,9002) 'eloss'//TRIM(cl), eloss(l)
   ENDDO l_410
   WRITE (iunit,9002) 'uvelmax', 100.0*uvelmax
   WRITE (iunit,9002) 'uvelmin', 100.0*uvelmin
   IF (nr.GT.2) THEN
      WRITE (iunit,9002) 'vvelmax', 100.0*vvelmax
      WRITE (iunit,9002) 'vvelmin', 100.0*vvelmin
   ENDIF
   IF (iopt_grid_nodim.EQ.3) THEN
      WRITE (iunit,9002) 'wphysmax', 1000.0*wphysmax
      WRITE (iunit,9002) 'wphysmin', 1000.0*wphysmin
   ENDIF
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
