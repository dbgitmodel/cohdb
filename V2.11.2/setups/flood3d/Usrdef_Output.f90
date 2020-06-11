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
! Author - ANTEA/MUMM
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.4
!
! $Date: 2018-07-23 10:21:02 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1167 $
!
! Description - output parameters for test case flood3d
!             - flooding/drying experiment over a 3-D obstacle
!
! Reference -
!
! Calling program - coherens_main
!
! Module calls - close_file, combine_stats_glb, define_out0d_vals, error_alloc,
!                max_vars, min_vars, mult_index, open_file, sum2_vars
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
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE model_output, ONLY: define_out0d_vals
USE paral_comms, ONLY: combine_stats_glb
USE paral_utilities, ONLY: max_vars, min_vars, sum2_vars
USE inout_routines, ONLY: close_file, open_file
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: outflag
CHARACTER (LEN=*), PARAMETER :: suffix = 'tst'
INTEGER, PARAMETER :: novars = 14
INTEGER, SAVE :: icout, iunit
INTEGER :: i, ii, iiloc, j
INTEGER, SAVE, DIMENSION(4) :: nhdims
REAL, SAVE :: dhun, volerr
REAL :: ddepint, depmaxerr, depmeanerr, depthsum, deptotint, dryarea, hour, &
      & transobu, uvelmax, uvelmin, vvelmax, vvelmin, wphysmax, wphysmin
REAL, DIMENSION(4) :: ecomps
REAL, DIMENSION(nobu) :: disobu
REAL, DIMENSION(novars) :: out0ddat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: dep_old
REAL, DIMENSION(ncloc,nrloc) :: ddeptot


procname(pglev+1) = 'usrdef_output'
CALL log_timer_in()

!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,1)%status = 'R'
   modfiles(io_modgrd,1,1)%form = 'U'
   modfiles(io_modgrd,1,1)%filename = modfiles(io_modgrd,1,2)%filename
   modfiles(io_modgrd,1,2)%status = '0'
!  ---initial conditions
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = modfiles(io_fincon,ics_phys,2)%form
   modfiles(io_inicon,ics_phys,1)%filename = &
 & modfiles(io_fincon,ics_phys,2)%filename
   modfiles(io_fincon,ics_phys,2)%status = '0'
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1,1)%status = 'R'
   modfiles(io_2uvobc,1,1)%form = modfiles(io_2uvobc,1,2)%form 
   modfiles(io_2uvobc,1,1)%filename = modfiles(io_2uvobc,1,2)%filename
   modfiles(io_2uvobc,1,2)%status = '0'
   GOTO 1001
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (nt.EQ.0) THEN

!
!2.1 Open files, write headers
!-----------------------------
!

   IF (.NOT.cold_start.AND.master) THEN

!     ---.tst file
      CALL open_file(iunit,TRIM(outtitle)//'.'//suffix,'OUT','A')
      WRITE (iunit,'(A)') 'Version: '//TRIM(model_version)
      WRITE (iunit,'(A)') 'Output parameters for test case flood3d:'//&
                          ' simulation '//TRIM(runtitle)
      WRITE (iunit,*)

   ENDIF

!
!2.2 Allocate
!------------
!

   ALLOCATE (dep_old(ncloc,nrloc),STAT=errstat)
   CALL error_alloc ('dep_old',2,(/ncloc,nrloc/),real_type)

!
!2.3 Initialise parameters and arrays
!------------------------------------
!

   icout = 600
   nhdims = 0
   dhun = surfacegrids(igrd_model,1)%delydat
   volerr = 0.0
   dep_old = deptotatc(1:ncloc,1:nrloc)

ENDIF

!
!3. Update parameters
!--------------------
!

outflag = mult_index(nt,icout).AND.nt.GT.0

!
!3.1 Volume error
!----------------
!
!---volume
CALL sum2_vars(deptotatc(1:ncloc,1:nrloc),deptotint,nhdims,'C  ',0,&
             & commall=.TRUE.,mask=seapoint(1:ncloc,1:nrloc))

!---open boundary transport term
iiloc_310: DO iiloc=1,nobuloc
   i = iobuloc(iiloc); j = jobuloc(iiloc)
   ii = indexobu(iiloc)
   disobu(ii) = 0.5*(udvel_old(i,j)+udvel(i,j))
ENDDO iiloc_310
IF (iopt_MPI.EQ.1) THEN
   CALL combine_stats_glb(disobu,nobu,nobuprocs,indexobuprocs,0,commall=.TRUE.)
ENDIF
transobu = SUM(disobu)

!---volume error
IF (nt.GT.0) THEN
   CALL sum2_vars(deptotatc(1:ncloc,1:nrloc),depthsum,nhdims,'C  ',0,&
                & commall=.TRUE.,mask=seapoint(1:ncloc,1:nrloc))
   ddeptot = deptotatc(1:ncloc,1:nrloc) - dep_old
   CALL sum2_vars(ddeptot,ddepint,nhdims,'C  ',0,commall=.TRUE.,&
                & mask=seapoint(1:ncloc,1:nrloc))
   volerr = volerr + ddepint - delt2d*transobu/dhun
   dep_old = deptotatc(1:ncloc,1:nrloc)
ENDIF

!
!3.2 Depth error
!---------------
!
IF (outflag) THEN
   CALL sum2_vars(deptotatc_err,depmeanerr,nhdims,'C  ',0,commall=.TRUE.,&
                & mask=seapoint(1:ncloc,1:nrloc))
   depmeanerr = depmeanerr/REAL(noseaatc)
   CALL max_vars(deptotatc_err,depmaxerr,0,commall=.TRUE.,&
               & mask=seapoint(1:ncloc,1:nrloc)) 
ENDIF   

!
!3.3 Dry area
!------------
!

IF (outflag) dryarea = 1.0-nowetatc/REAL(noseaatc)

!
!3.4 Energy
!----------
!

CALL define_out0d_vals(ecomps,4,ivarid=(/iarr_ekin0d,iarr_epot0d,iarr_etot0d,&
                     & iarr_edissip0d/))

!
!3.5 Current extrema
!-------------------
!

IF (outflag) THEN
   IF (iopt_grid_nodim.EQ.3) THEN
      CALL max_vars(uvel(1:ncloc,1:nrloc,:),uvelmax,iarr_uvel,commall=.TRUE.,&
                  & mask=node2du(1:ncloc,1:nrloc).GT.1)
      CALL min_vars(uvel(1:ncloc,1:nrloc,:),uvelmin,iarr_uvel,commall=.TRUE.,&
                  & mask=node2du(1:ncloc,1:nrloc).GT.1)
      CALL max_vars(vvel(1:ncloc,1:nrloc,:),vvelmax,iarr_vvel,commall=.TRUE.,&
                  & mask=node2dv(1:ncloc,1:nrloc).GT.1)
      CALL min_vars(vvel(1:ncloc,1:nrloc,:),vvelmin,iarr_vvel,commall=.TRUE.,&
                  & mask=node2dv(1:ncloc,1:nrloc).GT.1)
      CALL max_vars(wphys,wphysmax,iarr_wphys,commall=.TRUE.,mask=maskatc_int)
      CALL min_vars(wphys,wphysmin,iarr_wphys,commall=.TRUE.,mask=maskatc_int)
   ELSE
      CALL max_vars(umvel(1:ncloc,1:nrloc),uvelmax,iarr_umvel,commall=.TRUE.,&
                  & mask=node2du(1:ncloc,1:nrloc).GT.1)
      CALL min_vars(umvel(1:ncloc,1:nrloc),uvelmin,iarr_umvel,commall=.TRUE.,&
                  & mask=node2du(1:ncloc,1:nrloc).GT.1)
      CALL max_vars(vmvel(1:ncloc,1:nrloc),vvelmax,iarr_vmvel,commall=.TRUE.,&
                  & mask=node2dv(1:ncloc,1:nrloc).GT.1)
      CALL min_vars(vmvel(1:ncloc,1:nrloc),vvelmin,iarr_vmvel,commall=.TRUE.,&
                  & mask=node2dv(1:ncloc,1:nrloc).GT.1)
   ENDIF
ENDIF

!
!3.6 Store for output
!--------------------
!

IF (outflag) THEN
   out0ddat(1)  = 100.0*volerr/depthsum
   out0ddat(2)  = 100.0*depmeanerr
   out0ddat(3)  = 100.0*depmaxerr
   out0ddat(4)  = 100.0*dryarea
   out0ddat(5)  = 1.0E-09*ecomps(1)
   out0ddat(6)  = 1.0E-09*ecomps(2)
   out0ddat(7)  = 1.0E-09*ecomps(3)
   out0ddat(8)  = 0.0001*ecomps(4)
   out0ddat(9)  = 100.0*uvelmax
   out0ddat(10) = 100.0*uvelmin
   out0ddat(11) = 100.0*vvelmax
   out0ddat(12) = 100.0*vvelmin
   IF (iopt_grid_nodim.EQ.3) THEN
      out0ddat(13) = 1000.0*wphysmax
      out0ddat(14) = 1000.0*wphysmin
   ENDIF
ENDIF

!
!4. Write output
!---------------
!

IF (master.AND.outflag) THEN
   hour = nosecsrun/3600.0
   WRITE (iunit,9001) hour
   WRITE (iunit,9002) 'volerr',    out0ddat(1)
!   WRITE (iunit,9002) 'depmeanerr',out0ddat(2)
!   WRITE (iunit,9002) 'depmaxerr', out0ddat(3)
   WRITE (iunit,9002) 'dryarea',   out0ddat(4)
   WRITE (iunit,9002) 'ekin',      out0ddat(5)
   WRITE (iunit,9002) 'epot',      out0ddat(6)
   WRITE (iunit,9002) 'etot',      out0ddat(7)
   WRITE (iunit,9002) 'edissip',   out0ddat(8)
   WRITE (iunit,9002) 'uvelmax',   out0ddat(9)
   WRITE (iunit,9002) 'uvelmin',   out0ddat(10)
   WRITE (iunit,9002) 'vvelmax',   out0ddat(11)
   WRITE (iunit,9002) 'vvelmin',   out0ddat(12)
   IF (iopt_grid_nodim.EQ.3) THEN
      WRITE (iunit,9002) 'wphysmax',  out0ddat(13)
      WRITE (iunit,9002) 'wphysmin',  out0ddat(14)
   ENDIF
ENDIF

!
!5. Close output file
!--------------------
!

IF (master.AND.nt.EQ.nstep) CALL close_file(iunit,'A')

!
!6. Deallocate
!-------------
!

IF (cold_start.OR.nt.EQ.nstep) DEALLOCATE (dep_old)

1001 CALL log_timer_out()


RETURN

9001 FORMAT ('Time: ',F4.1,' hours')
9002 FORMAT (T3,A,T12,': ',G12.4E3)

END SUBROUTINE usrdef_output
