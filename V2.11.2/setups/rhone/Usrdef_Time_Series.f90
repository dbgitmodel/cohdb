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

!************************************************************************
!
! *Usrdef_Time_Series* Define time series output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.7.1
!
! $Date: 2016-10-10 16:36:15 +0200 (Mon, 10 Oct 2016) $
!
! $Revision: 973 $
!
! Description - test case rhone
!
! Reference -
!
! Routines - usrdef_tsr_params, usrdef_tsr0d_vals, usrdef_tsr2d_vals,
!            usrdef_tsr3d_vals
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_tsr_params
!************************************************************************
!
! *usrdef_out_params* Specifiers for time series output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.7.1
!
! Description - test case rhone
!
! Reference -
!
! Calling program - time_series_init
!
!************************************************************************
!
USE gridpars
USE iopars
USE modids
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: icount

procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---fortran name
tsrvars(1)%f90_name = 'area340d'
tsrvars(2)%f90_name = 'area220d'
tsrvars(3)%f90_name = 'ekin0d'
tsrvars(4)%f90_name = 'epot0d'
tsrvars(5)%f90_name = 'etot0d'
tsrvars(6)%f90_name = 'bdissip0d'
tsrvars(7)%f90_name = 'enstr0d'
tsrvars(8)%f90_name = 'curmax0d'
tsrvars(9)%f90_name = 'wmax0d'
tsrvars(10)%f90_name = 'wmin0d'
tsrvars(11)%f90_name = 'd340d'

!---standard name
tsrvars(1)%standard_name = 'area_bounded_by_surface_34_psu_contour'
tsrvars(2)%standard_name = 'area_bounded_by_surface_22_psu_contour'
tsrvars(3)%standard_name = 'kinetic_energy_in_sea_water'
tsrvars(4)%standard_name = 'available_potential_energy_in_sea_water'
tsrvars(5)%standard_name = 'total_energy_in_sea_water'
tsrvars(6)%standard_name = 'energy_dissipation_in_sea_water'
tsrvars(7)%standard_name = 'enstrophy_in_sea_water'
tsrvars(8)%standard_name = 'maximum_magnitude_of_horizontal_sea_water_velocity'
tsrvars(9)%standard_name = 'maximum_physical_upward_sea_water_velocity'
tsrvars(10)%standard_name = 'minimum_physical_upward_sea_water_velocity'
tsrvars(11)%standard_name = 'distance_of_34_psu_contour_to_coast_line'

!---comment
tsrvars(1:11)%comment = &
            & 'standard name constructed following the netCDF CF-1.6 guidelines'

!---long name
tsrvars(1)%long_name = 'Area between 34 contour'
tsrvars(2)%long_name = 'Area between 22 contour'
tsrvars(3)%long_name = 'Kinetic energy'
tsrvars(4)%long_name = 'Available potential energy'
tsrvars(5)%long_name = 'Total energy'
tsrvars(6)%long_name = 'Tidal dissipation'
tsrvars(7)%long_name = 'Enstrophy'
tsrvars(8)%long_name = 'Maximum horizontal current'
tsrvars(9)%long_name = 'Maximum vertical current'
tsrvars(10)%long_name = 'Minimum vertical current'
tsrvars(11)%long_name = 'Distance of 34 PSU line'

!---units
tsrvars(1)%units = 'km2'
tsrvars(2)%units = 'km2'
tsrvars(3)%units = 'GJ'
tsrvars(4)%units = 'GJ'
tsrvars(5)%units = 'GJ'
tsrvars(6)%units = 'MW'
tsrvars(7)%units = 'm3 s-2'
tsrvars(8)%units = 'cm s-1'
tsrvars(9)%units = 'mm s-1'
tsrvars(10)%units = 'mm s-1'
tsrvars(11)%units = 'km'

!---varids and ranks
tsrvars%ivarid = (/0,0,0,0,0,0,0,0,0,0,0,iarr_zeta,iarr_umvel,iarr_vmvel,&
                 & iarr_uvel,iarr_vvel,iarr_wphys,iarr_sal/)
tsrvars%nrank = (/0,0,0,0,0,0,0,0,0,0,0,2,2,2,3,3,3,3/)

!
!2. Variable indices
!-------------------
!

IF (iopt_verif.EQ.0) THEN
   ivarstsr(1,1:7) = (/12,13,14,15,16,17,18/)
   ivarstsr(2,1:11) = (/1,2,3,4,5,6,7,8,9,10,11/)
ELSE
   ivarstsr(1,:) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)
ENDIF

!
!3. File parameters
!------------------
!

IF (iopt_verif.EQ.0) THEN
    tsr2d(1)%defined = .TRUE.
    tsr3d(1)%defined = .TRUE.
    tsr0d(2)%defined = .TRUE.
ELSE
    tsr0d(1)%defined = .TRUE.
    tsr2d(1)%defined = .TRUE.
    tsr3d(1)%defined = .TRUE.
ENDIF

!
!4. Output grid
!--------------
!

IF (iopt_verif.EQ.0) THEN
!  ---first file
   icount = MERGE(4320,240,iopt_hydro_impl.EQ.0)
   tsrgpars(1)%tlims = (/icount,nstep,icount/)
!  ---second file
   icount = MERGE(180,10,iopt_hydro_impl.EQ.0)
   tsrgpars(2)%tlims = (/0,nstep,icount/)
ELSE
   icount = MERGE(2160,120,iopt_hydro_impl.EQ.0)
   tsrgpars(1)%tlims = (/0,nstep,icount/) 
ENDIF

tsrgpars(1:nosetstsr)%time_format = 3

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_tsr_params

!========================================================================

SUBROUTINE usrdef_tsr0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_tsr0d_vals* 0-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description - test case rhone
!
! Reference -
!
! Calling program - time_series
!
! Module calls - combine_mod, define_out0d_vals, max_vars, sum2_vars
!
!************************************************************************
!
USE currents
USE density
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE physpars
USE syspars
USE model_output, ONLY: define_out0d_vals
USE paral_comms, ONLY: combine_mod
USE paral_utilities, ONLY: max_vars, sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n0vars
REAL, INTENT(OUT), DIMENSION(n0vars) :: out0ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out0ddat* REAL    output data
!*n0vars*   INTEGER number of 0-D variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL, DIMENSION(nc,nr) :: maskglb
INTEGER :: i, j, j1, j2
INTEGER, DIMENSION(4) :: nhdims = 0
REAL :: fac
REAL, DIMENSION(nc,nr) :: gyglb
REAL, DIMENSION(nc,nr) :: salglb
REAL, DIMENSION(nc-1) :: dist
REAL, DIMENSION(ncloc,nrloc) :: area, array2dc1, array2dc2
REAL, DIMENSION(ncloc,nrloc,nz) :: array3dc


procname(pglev+1) = 'usrdef_tsr0d_vals'
CALL log_timer_in()

!---initialise
out0ddat = 0
WHERE (maskatc_int)
   area = delxatc(1:ncloc,1:nrloc)*delyatc(1:ncloc,1:nrloc)
ELSEWHERE
   array2dc1 = 0.0
   array2dc2 = 0.0
END WHERE

!---plume areas
i_210: DO i=1,ncloc
j_210: DO j=1,nrloc
   IF (maskatc_int(i,j)) THEN
      array2dc1(i,j) = MERGE (area(i,j),0.0,sal(i,j,nz).LE.34.0)
      array2dc2(i,j) = MERGE (area(i,j),0.0,sal(i,j,nz).LE.22.0)
   ENDIF
ENDDO j_210
ENDDO i_210
CALL sum2_vars(array2dc1,out0ddat(1),nhdims,'C  ',0)
CALL sum2_vars(array2dc2,out0ddat(2),nhdims,'C  ',0)
IF (master) out0ddat(1:2) = 1.0E-06*out0ddat(1:2)

!---energy components
CALL define_out0d_vals(out0ddat(3:7),5,ivarid=(/iarr_ekin0d,iarr_epot0d,&
                                    & iarr_etot0d,iarr_edissip0d,iarr_enstr0d/))
IF (master) THEN
   out0ddat(3:5) = 0.001*out0ddat(3:5)
   out0ddat(7) = 1.0E+06*out0ddat(7)
ENDIF

!---maximum horizontal current magnitude
array3dc = uvel(1:ncloc,1:nrloc,:)**2 + vvel(1:ncloc,1:nrloc,:)**2
CALL max_vars(array3dc,out0ddat(8),0,mask=maskatc_int)
out0ddat(8) = 100.0*SQRT(out0ddat(8))

!---maximum/minimum vertical current
CALL define_out0d_vals(out0ddat(9:10),2,&
                     & ivarid=(/iarr_wvel,iarr_wvel/),&
                     & oopt=(/oopt_max,oopt_min/))
IF (master) out0ddat(9:10) = 1000.0*out0ddat(9:10)

!---distance of 34 PSU line
CALL combine_mod(gyglb,gycoord(1:ncloc,1:nrloc),(/1,1/),iarr_gycoord,0.0)
CALL combine_mod(maskglb,maskatc_int,(/1,1/),0,.FALSE.)
CALL combine_mod(salglb,sal(1:ncloc,1:nrloc,nz),(/1,1/),iarr_sal,0.0)
IF (master) THEN
   i_220: DO i=1,nc-1
      dist(i) = 0.0
      j1 = nr
      DO WHILE (.NOT.maskglb(i,j1))
         j1 = j1 - 1
      ENDDO
      j1 = j1 + 1
      j2 = 1
      DO WHILE (salglb(i,j2).GT.34.0)
         j2 = j2 + 1
      ENDDO
      IF (j2.LT.j1) THEN
         fac = 0.001*Rearth*degtorad
         dist(i) = fac*(gyglb(i,j1)-gyglb(i,j2))
      ENDIF
   ENDDO i_220
   out0ddat(11) = MAXVAL(dist)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_tsr0d_vals

!========================================================================

SUBROUTINE usrdef_tsr2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_tsr2d_vals* 2-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description - test case rhone
!
! Reference -
!
! Calling program - time_series
!
!************************************************************************
!

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, n2vars
REAL, INTENT(OUT), DIMENSION(n2vars) :: out2ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out2ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*n2vars*   INTEGER number of 2-D variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_tsr2d_vals

!========================================================================

SUBROUTINE usrdef_tsr3d_vals(out3ddat,i,j,k,n3vars)
!************************************************************************
!
! *usrdef_tsr3d_vals* 3-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description - test case rhone
!
! Reference -
!
! Calling program - time_series
!
!************************************************************************
!

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, k, n3vars
REAL, INTENT(OUT), DIMENSION(n3vars) :: out3ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out3ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*k*        INTEGER Vertical index of output location
!*n3vars*   INTEGER number of 2-D variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_tsr3d_vals
