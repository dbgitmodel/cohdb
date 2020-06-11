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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.7.1
!
! $Date: 2014-12-05 17:20:26 +0100 (Fri, 05 Dec 2014) $
!
! $Revision: 780 $
!
! Description - test case obcest 
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
! *usrdef_tsr_params* Specifiers for time series output
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.7.1
!
! Description - test case obcest
!
! Reference -
!
! Calling program - time_series_init
!
!************************************************************************
!
USE iopars
USE modids
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: l

procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

!
!1. Output variables
!-------------------
!
!---fortran name
tsrvars(1)%f90_name = 'd30max0d'
tsrvars(2)%f90_name = 'dist300d'
tsrvars(3)%f90_name = 'd18max0d'
tsrvars(4)%f90_name = 'dist180d'
tsrvars(5)%f90_name = 'pdep10d'
tsrvars(6)%f90_name = 'pdep20d'
tsrvars(7)%f90_name = 'pdep30d'
tsrvars(8)%f90_name = 'area300d'
tsrvars(9)%f90_name = 'area150d'
tsrvars(10)%f90_name = 'umax0d'
tsrvars(11)%f90_name = 'wmax0d'
tsrvars(12)%f90_name = 'wmin0d'
tsrvars(13)%f90_name = 'ekin0d'
tsrvars(14)%f90_name = 'epot0d'
tsrvars(15)%f90_name = 'enstr0d'

!---long name
tsrvars(1)%standard_name = 'distance_of_30_psu_contour_to_western_boundary'
tsrvars(2)%standard_name = 'distance_of_30_psu_contour_to_southern_boundary'
tsrvars(3)%standard_name = 'distance_of_15_psu_contour_to_eastern_boundary'
tsrvars(4)%standard_name = 'distance_of_15_psu_contour_to_southern_boundary'
tsrvars(5)%standard_name = 'Halocline_depth at 2 km from western_boundary'
tsrvars(6)%standard_name = 'Halocline depth at 5 km from western_boundary'
tsrvars(7)%standard_name = 'Halocline depth at 8 km from western_boundary'
tsrvars(8)%standard_name = 'surface_area_fraction_above_30 PSU'
tsrvars(9)%standard_name = 'surface_area_fraction_below_15_PSU'
tsrvars(10)%standard_name = 'maximum_magnitude_of_horizontal_sea_water_velocity'
tsrvars(11)%standard_name = 'maximum_physical_upward_sea_water_velocity'
tsrvars(12)%standard_name = 'minimum_physical_upward_sea_water_velocity'
tsrvars(13)%standard_name = 'kinetic_energy_in_sea_water'
tsrvars(14)%standard_name = 'potential_energy_in_sea_water'
tsrvars(15)%standard_name = 'enstrophy_in_sea_water'

!---comment
tsrvars(1:15)%comment = &
            & 'standard name constructed following the netCDF CF-1.6 guidelines'

!---long name
tsrvars(1)%long_name = 'Distance of 30 PSU contour to western boundary'
tsrvars(2)%long_name = 'Distance of 30 PSU contour to southern boundary'
tsrvars(3)%long_name = 'Distance of 15 PSU contour to eastern boundary'
tsrvars(4)%long_name = 'Distance of 15 PSU contour to southern boundary'
tsrvars(5)%long_name = 'Halocline depth at 2 km from western boundary'
tsrvars(6)%long_name = 'Halocline depth at 5 km from western boundary'
tsrvars(7)%long_name = 'Halocline depth at 8 km from western boundary'
tsrvars(8)%long_name = 'Fraction of surface area above 30 PSU'
tsrvars(9)%long_name = 'Fraction of surface area below 15 PSU'
tsrvars(10)%long_name = 'Maximum horizontal current'
tsrvars(11)%long_name = 'Maximum vertical current'
tsrvars(12)%long_name = 'Minimum vertical current'
tsrvars(13)%long_name = 'Kinetic energy'
tsrvars(14)%long_name = 'Potential energy'
tsrvars(15)%long_name = 'Enstrophy'

!---units
tsrvars(1)%units = 'km'
tsrvars(2)%units = 'km'
tsrvars(3)%units = 'km'
tsrvars(4)%units = 'km'
tsrvars(5)%units = 'm'
tsrvars(6)%units = 'm'
tsrvars(7)%units = 'm'
tsrvars(8)%units = '1'
tsrvars(9)%units = '1'
tsrvars(10)%units = 'cm s-1'
tsrvars(11)%units = 'mm s-1'
tsrvars(12)%units = 'mm s-1'
tsrvars(13)%units = 'GJ'
tsrvars(14)%units = 'GJ'
tsrvars(15)%units = '1e6 m3 s-2'

!---varids and ranks
tsrvars(1:15)%ivarid = 0
tsrvars(16:24)%ivarid = (/iarr_zeta,iarr_umvel,iarr_vmvel,iarr_hmvelmag,&
                        & iarr_uvel,iarr_vvel,iarr_wphys,iarr_sal,iarr_hvelmag/)
tsrvars(1:15)%nrank = 0
tsrvars(16:24)%nrank = (/2,2,2,2,3,3,3,3,3/)

!
!2. Variable indices
!-------------------
!

ivarstsr(1,1:9) = (/(l,l=16,24)/)
ivarstsr(2,1:15) = (/(l,l=1,15)/)

!
!3. File parameters
!------------------
!

tsr0d(2)%defined = .TRUE.
tsr2d(1)%defined = .TRUE.
tsr3d(1)%defined = .TRUE.
 
!
!4. Output grid
!--------------
!
!---first file
tsrgpars(1)%time_format = 3
tsrgpars(1)%time_grid = .TRUE.
tsrgpars(1)%tlims = (/0,nstep,3600/)     

!---second file
tsrgpars(2)%time_format = 3
tsrgpars(2)%tlims = (/0,nstep,100/)

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
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.7
!
! Description - test case obcest
!
! Reference -
!
! Calling program - time_series
!
! Module calls - combine_mod, define_out0d_vals
!
!************************************************************************
!
USE density
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE model_output, ONLY: define_out0d_vals
USE paral_comms, ONLY: combine_mod
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
INTEGER :: i, ic, ivar, j, jc, k, kc
INTEGER, DIMENSION(1) :: jmax
INTEGER, DIMENSION(3) :: ipd
INTEGER, DIMENSION(nr-1) :: imax
REAL :: delh, hest, hriv, htot, scrit, x
LOGICAL, DIMENSION(nc,nr) :: maskglb
REAL, DIMENSION(3) :: pdep
REAL, DIMENSION(nr-1) :: dist
REAL, DIMENSION(nc,nr) :: array2dc1, array2dc2, deptotglb
REAL, DIMENSION(nc,nr,nz) :: salglb


procname(pglev+1) = 'usrdef_tsr0d_vals'
CALL log_timer_in()

!
!1. Initialise parameters 
!------------------------
!
!---domain sizes
delh = surfacegrids(igrd_model,1)%delxdat
htot = (nc-1)*delh
hest = 50*delh
hriv = htot - hest
ipd = (/20,50,80/)

!---global arrays
CALL combine_mod(maskglb,maskatc_int,(/1,1/),0,.FALSE.,commall=.TRUE.)
CALL combine_mod(deptotglb,deptotatc(1:ncloc,1:nrloc),(/1,1/),iarr_deptotatc,&
               & 0.0,commall=.TRUE.)
CALL combine_mod(salglb,sal(1:ncloc,1:nrloc,:),(/1,1,1/),iarr_sal,0.0,&
               & commall=.TRUE.)

!
!2. Maximum distance of plume contours
!-------------------------------------
!
!---30 PSU contour
scrit = 30.0; dist = 0.0
j_210: DO j=1,nr-1
   ic = 0
   i_211: DO i=1,nc
      IF (ic.EQ.0.AND.(.NOT.maskglb(i,j).OR.(maskglb(i,j).AND.&
        & salglb(i,j,nz).LE.scrit))) THEN
         ic = i
      ENDIF
   ENDDO i_211
   SELECT CASE(ic)
      CASE (1)
         dist(j) = 0.0; imax(j) = 0
      CASE (51)
         dist(j) = hest; imax(j) = 0
      CASE (0,101)
         dist(j) = htot; imax(j) = 0
      CASE DEFAULT
         x = (scrit-salglb(ic-1,j,nz))/(salglb(ic,j,nz)-salglb(ic-1,j,nz))
         dist(j) = delh*(ic-1.5+x)
         imax(j) = MERGE(ic-1,ic,x.LE.0.5)
   END SELECT
ENDDO j_210

out0ddat(1) = 0.001*MAXVAL(dist)
IF (ANY(imax.GT.0)) THEN
   jmax = MAXLOC(dist,dist.GT.0.0)
   out0ddat(2) = 0.001*delh*(jmax(1)-0.5)
ELSE
   out0ddat(2) = 0.0
ENDIF

!---18 PSU contour
scrit = 18.0; dist = 0.0
j_220: DO j=1,nr-1
   IF (ALL(salglb(1:nc-1,j,nz).GE.scrit.OR.(.NOT.maskglb(1:nc-1,j)))) THEN
      dist(j) = 0.0; imax(j) = 0
   ELSEIF (ALL(salglb(1:nc-1,j,nz).LT.scrit.OR.(.NOT.maskglb(1:nc-1,j)))) THEN
      dist(j) = htot; imax(j) = 0
   ELSE
      ic = 0
      i_221: DO i=nc-1,1,-1
         IF (ic.EQ.0.AND.maskglb(i,j).AND.salglb(i,j,nz).GE.scrit)  THEN
            ic = i + 1
         ENDIF
      ENDDO i_221
      x = (scrit-salglb(ic-1,j,nz))/(salglb(ic,j,nz)-salglb(ic-1,j,nz))
      dist(j) = delh*(nc+0.5-ic-x)
      imax(j) = MERGE(ic-1,ic,x.LE.0.5)
   ENDIF
ENDDO j_220

IF (ANY(imax.GT.0)) THEN
   out0ddat(3) = 0.001*MAXVAL(dist,dist.GT.0.0)
   jmax = MAXLOC(dist,dist.GT.0.0)
   out0ddat(4) = 0.001*delh*(jmax(1)-0.5)
ELSE
   out0ddat(3) = 0.0
   out0ddat(4) = 0.0
ENDIF

!
!3. Plume areas
!--------------
!

j_310: DO j=1,nr
i_310: DO i=1,nc
   IF (maskglb(i,j)) THEN
      array2dc1(i,j) = MERGE (1.0,0.0,salglb(i,j,nz).GE.30.0)
      array2dc2(i,j) = MERGE (1.0,0.0,salglb(i,j,nz).LE.18.0)
   ELSE
      array2dc1(i,j) = 0.0
   ENDIF
ENDDO i_310
ENDDO j_310
out0ddat(5) = 100.0*SUM(array2dc1,MASK=maskglb)/REAL(noseaatc)
out0ddat(6) = 100.0*SUM(array2dc2,MASK=maskglb)/REAL(noseaatc)

!
!4. Halocline depths
!-------------------
!

scrit = 30.0; jc = 5
ivar_410: DO ivar=1,3
   ic = ipd(ivar); kc = 0
   k_411: DO k=nz,1,-1
      IF ((kc.EQ.0).AND.(salglb(ic,jc,k).GT.scrit)) kc = k
   ENDDO k_411
   IF (kc.EQ.0) THEN
      pdep(ivar) = deptotglb(ic,jc)
   ELSEIF  (kc.EQ.nz) THEN
      pdep(ivar) = 0.0
   ELSE
      pdep(ivar) = (deptotglb(ic,jc)/REAL(nz))*&
                 & (nz-kc-0.5+(scrit-salglb(ic,jc,kc+1))/&
                 & (salglb(ic,jc,kc)-salglb(ic,jc,kc+1)))
   ENDIF
ENDDO ivar_410
out0ddat(7) = pdep(1); out0ddat(8) = pdep(2); out0ddat(9) = pdep(3)

!
!5. Current maxima/minima
!------------------------
!

CALL define_out0d_vals(out0ddat(10:12),3,&
                     & ivarid=(/iarr_hvelmag,iarr_wvel,iarr_wphys/),&
                     & oopt=(/oopt_max,oopt_max,oopt_min/))
out0ddat(10) = 100.0*out0ddat(10)
out0ddat(11:12) = 1000.0*out0ddat(11:12)

!
!6. Energy components
!--------------------
!

CALL define_out0d_vals(out0ddat(13:15),5,&
                     & ivarid=(/iarr_ekin0d,iarr_epot0d,iarr_enstr0d/))
out0ddat(13:14) = 0.001*out0ddat(13:14)
out0ddat(15) = 1.0E+06*out0ddat(15)

CALL log_timer_out()

      
RETURN

END SUBROUTINE usrdef_tsr0d_vals

!========================================================================

SUBROUTINE usrdef_tsr2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_tsr2d_vals* 2-D output data
!
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description -
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
! Author - ANTEA
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description -
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

