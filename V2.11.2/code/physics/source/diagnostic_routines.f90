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

MODULE diagnostic_routines
!************************************************************************
!
! *diagnostic_routines* Diagnostic tools for user output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.11.2
!
! $Date: 2018-01-19 09:43:56 +0100 (Fri, 19 Jan 2018) $
!
! $Revision: 1074 $
!
! Description -
!
! Routines - energy_0d, energy_2d, energy_3d, energy_flux_2d, energy_flux_3d,
!            energy_force_0d, energy_force_2d, energy_force_3d, enstrophy_0d,
!            vorticity_2d, vorticity_3d
!
!************************************************************************
!

IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE energy_0d(nodim,ecomps)
!************************************************************************
!
! *energy_0d* Domain integrated energy components
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.11.2
!
! Description -
!
! Module calls - error_abort, error_alloc, error_limits_var, sum2_vars,
!                Uarr_at_C, Varr_at_C, Zcoord_arr
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
USE physpars
USE switches
USE array_interp, ONLY: Uarr_at_C, Varr_at_C
USE error_routines, ONLY: error_abort, error_alloc, error_limits_var
USE grid_routines, ONLY: Zcoord_arr
USE paral_utilities, ONLY: sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
INTEGER, INTENT(IN) :: nodim
REAL, INTENT(OUT), DIMENSION(5) :: ecomps 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*nodim*   INTEGER Dimension of current field
!                  = 2 => depth-mean currents
!                  = 3 => 3-D currents
!*ecomps*  REAL    Energy components                                       [MJ]
!                  1: kinetic energy
!                  2: potential energy (barotropic)
!                  3: potential energy (baroclinic)
!                  4: baroclinic energy
!                  5: total energy
!
!------------------------------------------------------------------------------
!*Local variables
!
INTEGER :: k
REAL :: energy0d
INTEGER, DIMENSION(4) :: nhdims
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: area2d, area3d, array2dc, &
                                         & energy2d, rhodev2d, uvel2datc, &
                                         & vvel2datc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: energy3d, rhodev3d, &
                                           & uvel3datc, vvel3datc, zcoord


procname(pglev+1) = 'energy_0d'
CALL log_timer_in()

!
!1. Check argument
!-----------------
!

CALL error_limits_var(nodim,'nodim',2,3)
CALL error_abort('energy_0d',ierrno_arg)

!
!2. Allocate arrays
!------------------
!

ALLOCATE (area2d(ncloc,nrloc),STAT=errstat)
CALL error_alloc('area2d',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (area3d(ncloc,nrloc),STAT=errstat)
CALL error_alloc('area3d',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (energy2d(ncloc,nrloc),STAT=errstat)
CALL error_alloc('energy2d',2,(/ncloc,nrloc/),kndrtype)
IF (nodim.EQ.2) THEN
   ALLOCATE (uvel2datc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('uvel2datc',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (vvel2datc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('vvel2datc',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (nodim.EQ.3) THEN
   ALLOCATE (energy3d(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('energy3d',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (uvel3datc(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('uvel3datc',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (vvel3datc(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('vvel3datc',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF
IF (nodim.EQ.3.AND.iopt_dens.GT.0) THEN
   ALLOCATE (rhodev2d(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('rhodev2d',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (rhodev3d(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('rhodev3d',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (zcoord(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('zcoord',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

!
!3. Initialise
!-------------
!
!---parameters
ecomps = 0.0
nhdims = 0

!---domain area
WHERE (maskatc_int)
   area2d = 1.0E-06*garea
   area3d = area2d*deptotatc(1:ncloc,1:nrloc)
END WHERE

!---work space arrays
energy2d = 0.0
IF (nodim.EQ.3) energy3d = 0.0

!
!4. Kinetic energy
!-----------------
!
!4.1 2-D case
!------------
!

IF (nodim.EQ.2) THEN

   CALL Uarr_at_C(umvel(1:ncloc+1,1:nrloc),uvel2datc,1,1,(/1,1,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_umvel,.TRUE.)
   CALL Varr_at_C(vmvel(1:ncloc,1:nrloc+1),vvel2datc,1,1,(/1,1,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vmvel,.TRUE.)
   WHERE (maskatc_int)
      energy2d = area3d*(uvel2datc**2+vvel2datc**2)
   END WHERE
   CALL sum2_vars(energy2d,energy0d,nhdims,'C  ',0,commall=.TRUE.)
   ecomps(1) = 0.5*density_ref*energy0d

!
!4.2 3-D case
!------------
!

ELSEIF (nodim.EQ.3) THEN

   CALL Uarr_at_C(uvel(1:ncloc+1,1:nrloc,:),uvel3datc,1,1,(/1,1,1/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_uvel,.TRUE.)
   CALL Varr_at_C(vvel(1:ncloc,1:nrloc+1,:),vvel3datc,1,1,(/1,1,1/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vvel,.TRUE.)
   k_420: DO k=1,nz
      WHERE (maskatc_int)
         energy3d(:,:,k) = area2d*delzatc(1:ncloc,1:nrloc,k)*&
                         & (uvel3datc(:,:,k)**2+vvel3datc(:,:,k)**2)
      END WHERE
   ENDDO k_420
   CALL sum2_vars(energy3d,energy0d,nhdims,'C  ',0,commall=.TRUE.)
   ecomps(1) = 0.5*density_ref*energy0d

ENDIF

!
!5. Potential energy
!-------------------
!
!---barotropic component
WHERE (maskatc_int)
   array2dc = area2d*gaccatc(1:ncloc,1:nrloc)
   energy2d = array2dc*zeta(1:ncloc,1:nrloc)**2
END WHERE
CALL sum2_vars(energy2d,energy0d,nhdims,'C  ',0,commall=.TRUE.)
ecomps(2) = 0.5*density_ref*energy0d


!---baroclinic component
IF (nodim.EQ.3.AND.iopt_dens.GT.0) THEN
   WHERE (maskatc_int)
      rhodev2d = SUM(delzatc(1:ncloc,1:nrloc,:)*&
                 & (dens(1:ncloc,1:nrloc,:)/density_ref),DIM=3)&
                 & /deptotatc(1:ncloc,1:nrloc) - 1.0
      energy2d = array2dc*rhodev2d*&
              & (zeta(1:ncloc,1:nrloc)**2-depmeanatc(1:ncloc,1:nrloc)**2)
   END WHERE
   CALL sum2_vars(energy2d,energy0d,nhdims,'C  ',0,commall=.TRUE.)
   ecomps(3) = 0.5*density_ref*energy0d
ENDIF

!
!6. Baroclinic energy
!--------------------
!

IF (nodim.EQ.3.AND.iopt_dens.GT.0) THEN

!  ---vertical density deviation
   IF (iopt_sal.GT.0) THEN
      WHERE (maskatc_int)
         array2dc = SUM(delzatc(1:ncloc,1:nrloc,:)*sal(1:ncloc,1:nrloc,:),&
                      & DIM=3)/deptotatc(1:ncloc,1:nrloc)
      END WHERE
      k_610: DO k=1,nz
         WHERE (maskatc_int)
            rhodev3d(:,:,k) = beta_sal(1:ncloc,1:nrloc,k)*&
                            & (sal(1:ncloc,1:nrloc,k)-array2dc)
         END WHERE
      ENDDO k_610
   ELSE
      rhodev3d = 0.0
   ENDIF
   IF (iopt_temp.GT.0) THEN
      WHERE (maskatc_int)
         array2dc = SUM(delzatc(1:ncloc,1:nrloc,:)*temp(1:ncloc,1:nrloc,:),&
                      & DIM=3)/deptotatc(1:ncloc,1:nrloc)
      END WHERE
      k_620: DO k=1,nz
         WHERE (maskatc_int)
            rhodev3d(:,:,k) = rhodev3d(:,:,k) - beta_temp(1:ncloc,1:nrloc,k)*&
                            & (temp(1:ncloc,1:nrloc,k)-array2dc)
         END WHERE
      ENDDO k_620
   ENDIF

!  ---Z-coordinates
   CALL Zcoord_arr(zcoord,(/1,1,1/),(/ncloc,nrloc,nz/),'C  ',.FALSE.)

!  ---baroclinic energy
   WHERE (maskatc_int)
      array2dc = area2d*gaccatc(1:ncloc,1:nrloc)
   END WHERE
   k_630: DO k=1,nz
      WHERE (maskatc_int)
         energy3d(:,:,k) = array2dc*delzatc(1:ncloc,1:nrloc,k)*&
                         & rhodev3d(:,:,k)*zcoord(:,:,k)
      END WHERE
   ENDDO k_630
   CALL sum2_vars(energy3d,energy0d,nhdims,'C  ',0,commall=.TRUE.)
   ecomps(4) = density_ref*energy0d

ENDIF

!
!7. Total energy
!---------------
!

ecomps(5) = SUM(ecomps(1:4))

!
!8. Deallocate arrays
!--------------------
!

DEALLOCATE (area2d,area3d,array2dc,energy2d)
IF (nodim.EQ.2) DEALLOCATE (uvel2datc,vvel2datc)
IF (nodim.EQ.3) DEALLOCATE (energy3d,uvel3datc,vvel3datc)
IF (nodim.EQ.3.AND.iopt_dens.GT.0) THEN
   DEALLOCATE (rhodev2d,rhodev3d,zcoord)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE energy_0d

!========================================================================

SUBROUTINE energy_2d(i,j,nodim,ecomps)
!************************************************************************
!
! *energy_2d* Vertically integrated energy components at grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.0
!
! Description -
!
! Module calls - Uarr_at_C, Uvar_at_C, Varr_at_C, Vvar_at_C, Zcoord_arr
!
!************************************************************************
!
USE currents
USE density
USE depths
USE grid
USE gridpars
USE modids
USE physpars
USE switches
USE array_interp, ONLY: Uarr_at_C, Uvar_at_C, Varr_at_C, Vvar_at_C
USE grid_routines, ONLY: Zcoord_arr

!
!*Arguments
!
INTEGER, INTENT(IN) :: i,j, nodim
REAL, INTENT(OUT), DIMENSION(4) :: ecomps 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*i*       INTEGER X-index of C-node grid point
!*j*       INTEGER Y-index of C-node grid point
!*nodim*   INTEGER Dimension of current field
!                  = 2 => depth-mean currents
!                  = 3 => 3-D currents
!*ecomps*  REAL    Energy components                                    [J/m^2]
!                  1: kinetic energy
!                  2: potential energy
!                  3: baroclinic energy
!                  4: total energy
!
!------------------------------------------------------------------------------
!*Local variables
!
REAL :: halfdens, uvel2datc, vvel2datc
REAL, DIMENSION(nz) :: curmag, energy3d, saldef, tempdef, uvel3datc, &
                     & vvel3datc, zcoord


!
!1. Initialise
!-------------
!
!---parameters
halfdens = 0.5*density_ref
ecomps = 0.0
IF (nodeatc(i,j).EQ.0) RETURN

!
!2. Kinetic energy
!-----------------
!
!2.1 2-D case
!------------
!

IF (nodim.EQ.2) THEN

   uvel2datc =  Uvar_at_C(umvel(i:i+1,j),i,j,nz,1,1)
   vvel2datc =  Vvar_at_C(vmvel(i,j:j+1),i,j,nz,1,1)
   ecomps(1) = halfdens*deptotatc(i,j)*(uvel2datc**2+vvel2datc**2)

!
!2.2 3-D case
!------------
!

ELSEIF (nodim.EQ.3) THEN

   CALL Uarr_at_C(uvel(i:i+1,j,:),uvel3datc,1,1,(/i,j,1/),(/i+1,j,nz/),1,&
                & iarr_uvel,.FALSE.)
   CALL Varr_at_C(vvel(i,j:j+1,:),vvel3datc,1,1,(/i,j,1/),(/i,j+1,nz/),1,&
                & iarr_vvel,.FALSE.)
   curmag = delzatc(i,j,:)*(uvel3datc**2+vvel3datc**2)
   ecomps(1) = halfdens*SUM(curmag)

ENDIF

!
!3. Potential energy
!-------------------
!

ecomps(2) = halfdens*gaccatc(i,j)*zeta(i,j)**2

!
!4. Baroclinic energy
!--------------------
!

IF (nodim.EQ.3.AND.iopt_dens.GT.0) THEN
   IF (iopt_sal.GT.0) THEN
      saldef = sal(i,j,:) - SUM(delzatc(i,j,:)*sal(i,j,:))/deptotatc(i,j)
   ELSE
      saldef = 0.0
   ENDIF
   IF (iopt_temp.GT.0) THEN
      tempdef = temp(i,j,:) - SUM(delzatc(i,j,:)*temp(i,j,:))/deptotatc(i,j)
   ELSE
      tempdef = 0.0
   ENDIF
   CALL Zcoord_arr(zcoord,(/i,j,1/),(/i,j,nz/),'C  ',.FALSE.)
   energy3d = delzatc(i,j,:)*&
            & (beta_sal(i,j,:)*saldef-beta_temp(i,j,:)*tempdef)*zcoord
   ecomps(3) = density_ref*gaccatc(i,j)*SUM(energy3d)

ENDIF

!
!5. Total energy
!---------------
!

ecomps(4) = SUM(ecomps(1:3))


RETURN

END SUBROUTINE energy_2d

!========================================================================

SUBROUTINE energy_3d(i,j,k,ecomps)
!************************************************************************
!
! *energy_3d* Energy components at grid point (i,j,k)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.0
!
! Description -
!
! Module calls - Uvar_at_C, Vvar_at_C, Zcoord_var
!
!************************************************************************
!
USE currents
USE density
USE depths
USE grid
USE physpars
USE switches
USE array_interp, ONLY: Uvar_at_C, Vvar_at_C
USE grid_routines, ONLY: Zcoord_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: i,j, k
REAL, INTENT(OUT), DIMENSION(3) :: ecomps 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*i*       INTEGER X-index of C-node grid point
!*j*       INTEGER Y-index of C-node grid point
!*k*       INTEGER Z-index of C-node grid point
!*ecomps*  REAL    Energy components                                    [J/m^3]
!                  1: kinetic energy
!                  2: baroclinic energy
!                  3: total energy
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: saldef, tempdef, uvelatc, vvelatc, zcoord


!
!1. Initialise
!-------------
!
!---parameters
ecomps = 0.0
IF (nodeatc(i,j).EQ.0) RETURN

!
!2. Kinetic energy
!-----------------
!

uvelatc =  Uvar_at_C(uvel(i:i+1,j,k),i,j,k,1,1)
vvelatc =  Vvar_at_C(vvel(i,j:j+1,k),i,j,k,1,1)
ecomps(1) = 0.5*density_ref*(uvelatc**2+vvelatc**2)

!
!3. Baroclinic energy
!--------------------
!

IF (iopt_dens.GT.0) THEN
   IF (iopt_sal.GT.0) THEN
      saldef = sal(i,j,k) - SUM(delzatc(i,j,:)*sal(i,j,:))/deptotatc(i,j)
   ELSE
      saldef = 0.0
   ENDIF
   IF (iopt_temp.GT.0) THEN
      tempdef = temp(i,j,k) - SUM(delzatc(i,j,:)*temp(i,j,:))/deptotatc(i,j)
   ELSE
      tempdef = 0.0
   ENDIF
   zcoord = Zcoord_var(i,j,k,'C  ',.FALSE.)
   ecomps(2) = density_ref*gaccatc(i,j)*(beta_sal(i,j,k)*saldef-&
                                       & beta_temp(i,j,k)*tempdef)*zcoord
ENDIF

!
!4. Total energy
!---------------
!

ecomps(3) = ecomps(1) + ecomps(2)


RETURN

END SUBROUTINE energy_3d

!========================================================================

SUBROUTINE energy_flux_2d(i,j,nodim,fcomps)
!************************************************************************
!
! *energy_flux_2d* Vertically integrated energy flux components at grid point
!                  (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.0
!
! Description -
!
! Module calls - Uarr_at_C, Uvar_at_C, Varr_at_C, Vvar_at_C, Zcoord_arr
!
!************************************************************************
!
USE currents
USE density
USE depths
USE grid
USE gridpars
USE modids
USE physpars
USE switches
USE array_interp, ONLY: Uarr_at_C, Uvar_at_C, Varr_at_C, Vvar_at_C
USE grid_routines, ONLY: Zcoord_arr

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, nodim
REAL, INTENT(OUT), DIMENSION(4,2) :: fcomps

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*i*       INTEGER X-index of C-node grid point
!*j*       INTEGER Y-index of C-node grid point
!*nodim*   INTEGER Dimension of current field
!                  = 2 => depth-mean currents
!                  = 3 => 3-D currents
!*fcomps*  REAL    Energy flux components                                 [W/m]
!                  (1,*): kinetic energy
!                  (2,*): potential energy
!                  (3,*): baroclinic energy
!                  (4,*): total energy
!                  (*,1:2) : X-, Y-components
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k
REAL :: fmag, halfdens, uvel2datc, vvel2datc 
REAL, DIMENSION(nz) :: array1d, rhodef, saldef, tempdef, uvel3datc, &
                     & vvel3datc , zcoord


!
!1. Initialise
!-------------
!
!---parameters
halfdens = 0.5*density_ref
fcomps = 0.0
IF (nodeatc(i,j).EQ.0) RETURN

!---currents
IF (nodim.EQ.2) THEN
   uvel2datc =  Uvar_at_C(umvel(i:i+1,j),i,j,nz,1,1)
   vvel2datc =  Vvar_at_C(vmvel(i,j:j+1),i,j,nz,1,1)
ELSEIF (nodim.EQ.3) THEN
   CALL Uarr_at_C(uvel(i:i+1,j,:),uvel3datc,1,1,(/i,j,1/),(/i+1,j,nz/),1,&
                & iarr_uvel,.FALSE.)
   CALL Varr_at_C(vvel(i,j:j+1,:),vvel3datc,1,1,(/i,j,1/),(/i,j+1,nz/),1,&
                & iarr_vvel,.FALSE.)
ENDIF

!
!2. Kinetic energy flux
!----------------------
!
!2.1 2-D case
!------------
!

IF (nodim.EQ.2) THEN

   fmag = halfdens*deptotatc(i,j)*(uvel2datc**2+vvel2datc**2)
   fcomps(1,1) = fmag*uvel2datc
   fcomps(1,2) = fmag*vvel2datc

!
!2.2 3-D case
!------------
!

ELSEIF (nodim.EQ.3) THEN

   array1d = delzatc(i,j,:)*(uvel3datc**2+vvel3datc**2)
   fcomps(1,1) = halfdens*SUM(array1d*uvel3datc)
   fcomps(1,2) = halfdens*SUM(array1d*vvel3datc)

ENDIF

!
!3. Potential energy flux
!------------------------
!
!3.1 2-D case
!------------
!

IF (nodim.EQ.2) THEN
   
   fmag = density_ref*deptotatc(i,j)*gaccatc(i,j)*zeta(i,j)
   fcomps(2,1) = fmag*uvel2datc
   fcomps(2,2) = fmag*vvel2datc

!
!3.2 3-D case
!------------
!

ELSEIF (nodim.EQ.3) THEN

   fmag = density_ref*gaccatc(i,j)*zeta(i,j)
   fcomps(2,1) = fmag*SUM(delzatc(i,j,:)*uvel3datc)
   fcomps(2,2) = fmag*SUM(delzatc(i,j,:)*vvel3datc)

ENDIF

!
!4. Baroclinic energy flux
!-------------------------
!

IF (nodim.EQ.3.AND.iopt_dens.GT.0) THEN
   
   IF (iopt_sal.GT.0) THEN
      saldef = sal(i,j,:) - SUM(delzatc(i,j,:)*sal(i,j,:))/deptotatc(i,j)
   ELSE
      saldef = 0.0
   ENDIF
   IF (iopt_temp.GT.0) THEN
      tempdef = temp(i,j,:) - SUM(delzatc(i,j,:)*temp(i,j,:))/deptotatc(i,j)
   ELSE
      tempdef = 0.0
   ENDIF
   rhodef = beta_sal(i,j,:)*(sal(i,j,:)-saldef)-beta_temp(i,j,:)*&
                          & (temp(i,j,:)-tempdef)
   CALL Zcoord_arr(zcoord,(/i,j,1/),(/i,j,nz/),'C  ',.FALSE.)
   array1d(nz) = 0.0
   k_410: DO k=nz-1,1,-1
      array1d(k) = array1d(k+1) + rhodef(k+1)*delzatc(i,j,k+1)
   ENDDO k_410
   array1d = delzatc(i,j,:)*(array1d+rhodef*(zcoord+0.5*delzatc(i,j,:)))
   fcomps(3,1) = density_ref*gaccatc(i,j)*SUM(array1d*uvel3datc)
   fcomps(3,2) = density_ref*gaccatc(i,j)*SUM(array1d*vvel3datc)

ENDIF

!
!5. Total energy flux
!--------------------
!

fcomps(4,:) = SUM(fcomps(1:3,:),DIM=1)


RETURN

END SUBROUTINE energy_flux_2d

!========================================================================

SUBROUTINE energy_flux_3d(i,j,k,fcomps)
!************************************************************************
!
! *energy_flux_3d* Energy flux components at grid point (i,j,k)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.0
!
! Description -
!
! Module calls - Uvar_at_C, Vvar_at_C, Wvar_at_C, Zcoord_var
!
!************************************************************************
!
USE currents
USE density
USE depths
USE grid
USE gridpars
USE physpars
USE switches
USE array_interp, ONLY: Uvar_at_C, Vvar_at_C, Wvar_at_C
USE grid_routines, ONLY: Zcoord_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, k
REAL, INTENT(OUT), DIMENSION(4,3) :: fcomps

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*i*       INTEGER X-index of C-node grid point
!*j*       INTEGER Y-index of C-node grid point
!*k*       INTEGER Z-index of C-node grid point
!*fcomps*  REAL    Energy flux components                               [W/m^2]
!                  (1,*): kinetic energy
!                  (2,*): potential energy
!                  (3,*): baroclinic energy
!                  (4,*): total energy
!                  (*,1:3) : X-, Y-, Z-components
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: l
REAL :: fmag, uvelatc, vvelatc, wvelatc, zcoord
REAL, DIMENSION(nz) :: array1d, rhodef, saldef, tempdef


!
!1. Initialise
!-------------
!
!---parameters
fcomps = 0.0
IF (nodeatc(i,j).EQ.0) RETURN

!---currents
uvelatc =  Uvar_at_C(uvel(i:i+1,j,k),i,j,k,1,1)
vvelatc =  Vvar_at_C(vvel(i,j:j+1,k),i,j,k,1,1)

!
!2. Kinetic energy flux
!----------------------
!

wvelatc = Wvar_at_C(wvel(i,j,k:k+1),i,j)
fmag = 0.5*density_ref*(uvelatc**2+vvelatc**2)
fcomps(1,1) = fmag*uvelatc
fcomps(1,2) = fmag*vvelatc
fcomps(1,3) = fmag*wvelatc

!
!3. Potential energy flux
!------------------------
!

fmag = density_ref*gaccatc(i,j)*zeta(i,j)
fcomps(2,1) = fmag*uvelatc
fcomps(2,2) = fmag*vvelatc
fcomps(2,3) = fmag*wvelatc

!
!4. Baroclinic energy flux
!-------------------------
!

IF (iopt_dens.GT.0) THEN
   
   IF (iopt_sal.GT.0) THEN
      saldef(k:nz) = sal(i,j,k:nz) - SUM(delzatc(i,j,:)*sal(i,j,:))&
                                      & /deptotatc(i,j)
   ELSE
      saldef(k:nz) = 0.0
   ENDIF
   IF (iopt_temp.GT.0) THEN
      tempdef(k:nz) = temp(i,j,k:nz) - SUM(delzatc(i,j,:)*temp(i,j,:))&
                                        & /deptotatc(i,j)
   ELSE
      tempdef(k:nz) = 0.0
   ENDIF
   rhodef(k:nz) = beta_sal(i,j,k:nz)*(sal(i,j,k:nz)-saldef(k:nz))-&
                & beta_temp(i,j,k:nz)*(temp(i,j,k:nz)-tempdef(k:nz))
   zcoord =  Zcoord_var(i,j,k,'C  ',.FALSE.)
   array1d(nz) = 0.0
   l_510: DO l=nz-1,k,-1
      array1d(l) = array1d(l+1) + rhodef(l+1)*delzatc(i,j,l+1)
   ENDDO l_510
   fmag = density_ref*gaccatc(i,j)*deptotatc(i,j)*(array1d(k)+rhodef(k)*&
                                 & (zcoord+0.5*delzatc(i,j,k)))
   fcomps(3,1) = fmag*uvelatc
   fcomps(3,2) = fmag*vvelatc
   fcomps(3,3) = fmag*wvelatc

ENDIF

!
!6. Total energy flux
!--------------------
!

fcomps(4,:) = SUM(fcomps(1:3,:),DIM=1)


RETURN

END SUBROUTINE energy_flux_3d

!========================================================================

SUBROUTINE energy_force_0d(nodim,wcomps)
!************************************************************************
!
! *energy_force_0d* Domain integrated forcing terms in energy balance
!                   equation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.11.2
!
! Description -
!
! Module calls - error_abort, error_alloc, error_limits_var, sum2_vars,
!                Uarr_at_C, Varr_at_C, Zcoord_arr
!
!************************************************************************
!
USE currents
USE density
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE optics
USE physpars
USE switches
USE tide
USE array_interp, ONLY: Uarr_at_C, Varr_at_C
USE error_routines, ONLY: error_abort, error_alloc, error_limits_var
USE grid_routines, ONLY: Zcoord_arr
USE paral_utilities, ONLY: sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
INTEGER, INTENT(IN) :: nodim
REAL, INTENT(OUT), DIMENSION(6) :: wcomps 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*nodim*   INTEGER Dimension of current field
!                  = 2 => depth-mean currents
!                  = 3 => 3-D currents
!*wcomps*  REAL    Forcing terms in the energy balance equation            [MW]
!                  1: work done by tidal force
!                  2: work done by surface stress
!                  3: work done by bottom stress
!                  4: work done by internal shear stresses
!                  5: work done by buoyancy fluxes
!                  6: total work
!
!------------------------------------------------------------------------------
!*Local variables
!
INTEGER :: k
REAL :: xfac, work0d
INTEGER, DIMENSION(4) :: nhdims
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: area, uvel2datc, vvel2datc, work2d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: flux3d, uvel3datc, vvel3datc, &
                                           & zcoord


procname(pglev+1) = 'energy_force_0d'
CALL log_timer_in()

!
!1. Check argument
!-----------------
!

CALL error_limits_var(nodim,'nodim',2,3)
CALL error_abort('energy_force_0d',ierrno_arg)

!
!2. Allocate arrays
!------------------
!

ALLOCATE (area(ncloc,nrloc),STAT=errstat)
CALL error_alloc('area',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (work2d(ncloc,nrloc),STAT=errstat)
CALL error_alloc('work2d',2,(/ncloc,nrloc/),kndrtype)
IF (nodim.EQ.2.OR.iopt_astro_tide.EQ.1) THEN
   ALLOCATE (uvel2datc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('uvel2datc',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (vvel2datc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('vvel2datc',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (nodim.EQ.3) THEN
   ALLOCATE (flux3d(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('flux3d',3,(/ncloc,nrloc,nz+1/),kndrtype)
   ALLOCATE (uvel3datc(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('uvel3datc',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (vvel3datc(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('vvel3datc',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF
IF (nodim.EQ.3.AND.iopt_dens.GT.0) THEN
   ALLOCATE (zcoord(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('zcoord',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

!
!3. Initialise
!-------------
!
!---parameters
wcomps = 0
nhdims = 0

!---domain area
WHERE (maskatc_int)
   area = 1.0E-06*garea
END WHERE

!---currents
IF (nodim.EQ.2.OR.iopt_astro_tide.EQ.1) THEN
   CALL Uarr_at_C(umvel(1:ncloc+1,1:nrloc),uvel2datc,1,1,(/1,1,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_umvel,.TRUE.)
   CALL Varr_at_C(vmvel(1:ncloc,1:nrloc+1),vvel2datc,1,1,(/1,1,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vmvel,.TRUE.)
ENDIF
IF (nodim.EQ.3) THEN
   CALL Uarr_at_C(uvel(1:ncloc+1,1:nrloc,:),uvel3datc,1,1,(/1,1,1/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_uvel,.TRUE.)
   CALL Varr_at_C(vvel(1:ncloc,1:nrloc+1,:),vvel3datc,1,1,(/1,1,1/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vvel,.TRUE.)
ENDIF

!---Z-coordinates
IF (nodim.EQ.3.AND.iopt_dens.GT.0) THEN
   xfac = density_ref*specheat
   CALL Zcoord_arr(zcoord,(/1,1,1/),(/ncloc,nrloc,nz/),'C  ',.FALSE.)
ENDIF

!---work space
work2d = 0.0

!
!4. Work done by tidal force
!---------------------------
!

IF (iopt_astro_tide.EQ.1) THEN

   WHERE (maskatc_int)
      work2d = area*deptotatc(1:ncloc,1:nrloc)*&
             & (uvel2datc*fxastro(1:ncloc,:)+vvel2datc*fyastro(:,1:nrloc))
   END WHERE
   CALL sum2_vars(work2d,work0d,nhdims,'C  ',0,commall=.TRUE.)
   wcomps(1) = density_ref*work0d
ENDIF

!
!5. Work done by surface stress
!------------------------------
!

IF (iopt_meteo_stres.GT.0) THEN

!
!5.1 2-D case
!------------
!

   IF (nodim.EQ.2) THEN
   
      WHERE (maskatc_int)
         work2d = area*(uvel2datc*usstresatc(1:ncloc,1:nrloc) + &
                & vvel2datc*vsstresatc(1:ncloc,1:nrloc))
      END WHERE
      CALL sum2_vars(work2d,work0d,nhdims,'C  ',0,commall=.TRUE.)
      wcomps(2) = density_ref*work0d

!
!5.2 3-D case
!------------
!

   ELSEIF (nodim.EQ.3) THEN

      WHERE (maskatc_int)
         work2d = area*(uvel3datc(:,:,nz)*usstresatc(1:ncloc,1:nrloc) + &
                      & vvel3datc(:,:,nz)*vsstresatc(1:ncloc,1:nrloc))
      END WHERE
      CALL sum2_vars(work2d,work0d,nhdims,'C  ',0,commall=.TRUE.)
      wcomps(2) = density_ref*work0d

   ENDIF

ENDIF

!
!6. Work done by bottom stress
!-----------------------------
!

IF (iopt_bstres_form.EQ.2) THEN

!
!6.1 2-D case
!------------
!

   IF (nodim.EQ.2) THEN

      WHERE (maskatc_int)
         work2d = bdragcoefatc(1:ncloc,1:nrloc)*area*&
                & (uvel2datc**2+vvel2datc**2)**1.5
      END WHERE
      CALL sum2_vars(work2d,work0d,nhdims,'C  ',0,commall=.TRUE.)
      wcomps(3) = -density_ref*work0d

!
!6.2 3-D case
!------------
!

   ELSEIF (nodim.EQ.3) THEN

      WHERE (maskatc_int)
         work2d = bdragcoefatc(1:ncloc,1:nrloc)*area*&
                & (uvel3datc(:,:,1)**2+vvel3datc(:,:,1)**2)**1.5
      END WHERE
      CALL sum2_vars(work2d,work0d,nhdims,'C  ',0,commall=.TRUE.)
      wcomps(3) = -density_ref*work0d

   ENDIF

ENDIF

!
!7. Work done by internal shear stresses
!---------------------------------------
!

IF (nodim.EQ.3) THEN

   WHERE (maskatc_int)
      flux3d(:,:,1) = 0.0
      flux3d(:,:,nz+1) = 0.0
      work2d = 0.0
   END WHERE

!
!7.1 X-component
!---------------
!

   k_711: DO k=2,nz
      WHERE (maskatc_int)
         flux3d(:,:,k) = vdifcoefmom(1:ncloc,1:nrloc,k)*&
                      & (uvel3datc(:,:,k)-uvel3datc(:,:,k-1))&
                      & /delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_711
   k_712: DO k=1,nz
      WHERE (maskatc_int)
         work2d = work2d + uvel3datc(:,:,k)*(flux3d(:,:,k+1)-flux3d(:,:,k))
      END WHERE
   ENDDO k_712

!
!7.2 Y-component
!---------------
!

   k_721: DO k=2,nz
      WHERE (maskatc_int)
         flux3d(:,:,k) = vdifcoefmom(1:ncloc,1:nrloc,k)*&
                       & (vvel3datc(:,:,k)-vvel3datc(:,:,k-1))&
                       & /delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_721
   k_722: DO k=1,nz
      WHERE (maskatc_int)
         work2d = work2d + vvel3datc(:,:,k)*(flux3d(:,:,k+1)-flux3d(:,:,k))
      END WHERE
   ENDDO k_722

!
!7.3 Domain average
!------------------
!

   WHERE (maskatc_int)
      work2d = work2d*area
   END WHERE
   CALL sum2_vars(work2d,work0d,nhdims,'C  ',0,commall=.TRUE.)
   wcomps(4) = density_ref*work0d

ENDIF

!
!8. Work done by buoyancy fluxes
!-------------------------------
!

IF (nodim.EQ.3.AND.iopt_dens.GT.0) THEN
   

   WHERE (maskatc_int)
      work2d = 0.0
   END WHERE

!
!8.1 Salinity fluxes
!-------------------
!

   IF (iopt_sal.GT.0) THEN
      WHERE (maskatc_int)
         flux3d(:,:,nz+1) = beta_sal(1:ncloc,1:nrloc,nz)*ssalflux
      END WHERE
      k_811: DO k=2,nz
         WHERE (maskatc_int)
            flux3d(:,:,k) = vdifcoefscal(1:ncloc,1:nrloc,k)*&
                          & beta_sal(1:ncloc,1:nrloc,k)*&
                          & (sal(1:ncloc,1:nrloc,k)-sal(1:ncloc,1:nrloc,k-1))&
                          & /delzatw(1:ncloc,1:nrloc,k)
         END WHERE
      ENDDO k_811
      k_812: DO k=1,nz
         WHERE (maskatc_int)
            work2d = work2d - zcoord(:,:,k)*(flux3d(:,:,k+1)-flux3d(:,:,k))
         END WHERE
      ENDDO k_812
   ENDIF

!
!8.2 Temperature fluxes
!----------------------
!

   IF (iopt_temp.GT.0) THEN

      WHERE (maskatc_int)
         flux3d(:,:,nz+1) = beta_temp(1:ncloc,1:nrloc,nz)*&
                          & (qrad-qnonsol)/xfac
      END WHERE
      k_821: DO k=2,nz
         WHERE (maskatc_int)
            flux3d(:,:,k) = vdifcoefscal(1:ncloc,1:nrloc,k)*&
                          & beta_temp(1:ncloc,1:nrloc,k)*&
                         & (temp(1:ncloc,1:nrloc,k)-temp(1:ncloc,1:nrloc,k-1))&
                         & /delzatw(1:ncloc,1:nrloc,k)
         END WHERE
      ENDDO k_821
      k_822: DO k=1,nz
         WHERE (maskatc_int)
            work2d = work2d + zcoord(:,:,k)*(flux3d(:,:,k+1)-flux3d(:,:,k))
         END WHERE
      ENDDO k_822
   ENDIF

!
!8.3 Domain average
!------------------
!

   WHERE (maskatc_int)
      work2d = work2d*area*gaccatc(1:ncloc,1:nrloc)
   END WHERE
   CALL sum2_vars(work2d,work0d,nhdims,'C  ',0,commall=.TRUE.)
   wcomps(5) = density_ref*work0d

ENDIF

!
!9. Total work
!-------------
!

wcomps(6) = SUM(wcomps(1:5))

!
!10. Deallocate arrays
!---------------------
!

DEALLOCATE (area,work2d)
IF (nodim.EQ.2.OR.iopt_astro_tide.EQ.1) DEALLOCATE (uvel2datc,vvel2datc)
IF (nodim.EQ.3) DEALLOCATE (flux3d,uvel3datc,vvel3datc)
IF (nodim.EQ.3.AND.iopt_dens.GT.0) DEALLOCATE (zcoord)

CALL log_timer_out()


RETURN

END SUBROUTINE energy_force_0d

!========================================================================

SUBROUTINE energy_force_2d(i,j,nodim,wcomps)
!************************************************************************
!
! *energy_force_2d* Vertically integrated forcing terms in energy balance
!                   equation at grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.0
!
! Description -
!
! Module calls - Uarr_at_C, Uvar_at_C, Varr_at_C, Vvar_at_C, Zcoord_arr
!
!************************************************************************
!
USE currents
USE density
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE modids
USE obconds
USE optics
USE physpars
USE switches
USE tide
USE array_interp, ONLY: Uarr_at_C, Uvar_at_C, Varr_at_C, Vvar_at_C
USE grid_routines, ONLY: Zcoord_arr

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, nodim
REAL, INTENT(OUT), DIMENSION(6) :: wcomps 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*i*       INTEGER X-index of C-node grid point
!*j*       INTEGER Y-index of C-node grid point
!*nodim*   INTEGER Dimension of current field
!                  = 2 => depth-mean currents
!                  = 3 => 3-D currents
!*wcomps*  REAL    Forcing terms in the energy balance equation         [W/m^2]
!                  1: work done by tidal force
!                  2: work done by surface stress
!                  3: work done by bottom stress
!                  4: work done by internal shear stresses
!                  5: work done by density fluxes
!                  6: total work
!
!------------------------------------------------------------------------------
!*Local variables
!
INTEGER :: k
REAL :: uvel2datc, vvel2datc, work
REAL, DIMENSION(nz) :: uvel3datc, vvel3datc, zcoord
REAL, DIMENSION(nz+1) :: flux3d


!
!1. Initialise
!-------------
!
!---parameters
wcomps = 0
IF (nodeatc(i,j).EQ.0) RETURN

!---currents
IF (nodim.EQ.2.OR.iopt_astro_tide.EQ.1) THEN
   uvel2datc =  Uvar_at_C(umvel(i:i+1,j),i,j,nz,1,1)
   vvel2datc =  Vvar_at_C(vmvel(i,j:j+1),i,j,nz,1,1)
ENDIF
IF (nodim.EQ.3) THEN
   CALL Uarr_at_C(uvel(i:i+1,j,:),uvel3datc,1,1,(/i,j,1/),(/i+1,j,nz/),1,&
                & iarr_uvel,.FALSE.)
   CALL Varr_at_C(vvel(i,j:j+1,:),vvel3datc,1,1,(/i,j,1/),(/i,j+1,nz/),1,&
                & iarr_vvel,.FALSE.)
ENDIF

!---Z-coordinates
IF (nodim.EQ.3.AND.iopt_dens.GT.0) THEN
   CALL Zcoord_arr(zcoord,(/i,j,1/),(/i,j,nz/),'C  ',.FALSE.)
ENDIF

!
!2. Work done by tidal force
!---------------------------
!

IF (iopt_astro_tide.EQ.1) THEN
   wcomps(1) = density_ref*deptotatc(i,j)*&
             & (uvel2datc*fxastro(i,j)+vvel2datc*fyastro(i,j))
ENDIF

!
!3. Work done by surface stress
!------------------------------
!

IF (iopt_meteo_stres.GT.0) THEN

!
!3.1 2-D case
!------------
!

   IF (nodim.EQ.2) THEN
      wcomps(2) = density_ref*(uvel2datc*usstresatc(i,j) + &
                             & vvel2datc*vsstresatc(i,j))

!
!3.2 3-D case
!------------
!

   ELSEIF (nodim.EQ.3) THEN
      wcomps(2) = density_ref*(uvel3datc(nz)*usstresatc(i,j) + &
                             & vvel3datc(nz)*vsstresatc(i,j))
   ENDIF

ENDIF

!
!4. Work done by bottom stress
!-----------------------------
!

IF (iopt_bstres_form.EQ.2) THEN

!
!4.1 2-D case
!------------
!

   IF (nodim.EQ.2) THEN
      wcomps(3) = -density_ref*bdragcoefatc(i,j)*&
                & (uvel2datc**2+vvel2datc**2)**1.5

!
!4.2 3-D case
!------------
!

   ELSEIF (nodim.EQ.3) THEN
      wcomps(3) = -density_ref*bdragcoefatc(i,j)*&
                & (uvel3datc(1)**2+vvel3datc(1)**2)**1.5

   ENDIF

ENDIF

!
!5. Work done by internal shear stresses
!---------------------------------------
!

IF (nodim.EQ.3) THEN

   flux3d(1) = 0.0
   flux3d(nz+1) = 0.0
   work = 0.0

!
!5.1 X-component
!---------------
!

   k_511: DO k=2,nz
      flux3d(k) = vdifcoefmom(i,j,k)*(uvel3datc(k)-uvel3datc(k-1))&
                & /delzatw(i,j,k)
   ENDDO k_511
   k_512: DO k=1,nz
      work = work + uvel3datc(k)*(flux3d(k+1)-flux3d(k))
   ENDDO k_512

!
!5.2 Y-component
!---------------
!

   k_521: DO k=2,nz
      flux3d(k) = vdifcoefmom(i,j,k)*(vvel3datc(k)-vvel3datc(k-1))&
                & /delzatw(i,j,k)
   ENDDO k_521
   k_522: DO k=1,nz
      work = work + vvel3datc(k)*(flux3d(k+1)-flux3d(k))
   ENDDO k_522

   wcomps(4) = density_ref*work

ENDIF

!
!6. Work done by buoyancy fluxes
!-------------------------------
!

IF (nodim.EQ.3.AND.iopt_dens.GT.0) THEN

   work = 0.0

!
!6.1 Salinity fluxes
!-------------------
!

   IF (iopt_sal.GT.0) THEN
      flux3d(nz+1) = beta_sal(i,j,nz)*ssalflux(i,j)
      k_611: DO k=2,nz
         flux3d(k) = vdifcoefscal(i,j,k)*beta_sal(i,j,k)*&
                   & (sal(i,j,k)-sal(i,j,k-1))/delzatw(i,j,k)
      ENDDO k_611
      k_612: DO k=1,nz
         work = work - zcoord(k)*(flux3d(k+1)-flux3d(k))
      ENDDO k_612
   ENDIF

!
!6.2 Temperature fluxes
!----------------------
!

   IF (iopt_temp.GT.0) THEN
      flux3d(nz+1) = beta_temp(i,j,nz)*(qrad(i,j)-qnonsol(i,j))&
                   & /(density_ref*specheat)
      k_621: DO k=2,nz
         flux3d(k) = vdifcoefscal(i,j,k)*beta_temp(i,j,k)*&
                   & (temp(i,j,k)-temp(i,j,k-1))/delzatw(i,j,k)
      ENDDO k_621
      k_622: DO k=1,nz
         work = work + zcoord(k)*(flux3d(k+1)-flux3d(k))
      ENDDO k_622
   ENDIF

   wcomps(5) = density_ref*work

ENDIF

!
!7. Total work
!-------------
!

wcomps(6) = SUM(wcomps(1:5))


RETURN

END SUBROUTINE energy_force_2d

!========================================================================

SUBROUTINE energy_force_3d(i,j,k,wcomps)
!************************************************************************
!
! *energy_force_3d* Forcing terms in energy balance equation at grid point
!                   (i,j,k)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.2
!
! Description -
!
! Module calls - Uarr_at_C, Varr_at_C, Zcoord_var
!
!************************************************************************
!
USE currents
USE density
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE modids
USE obconds
USE optics
USE physpars
USE switches
USE tide
USE array_interp, ONLY: Uarr_at_C, Varr_at_C
USE grid_routines, ONLY: Zcoord_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, k
REAL, INTENT(OUT), DIMENSION(4) :: wcomps 

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*i*       INTEGER X-index of C-node grid point
!*j*       INTEGER Y-index of C-node grid point
!*k*       INTEGER Z-index of C-node grid point
!*wcomps*  REAL    Forcing terms in the energy balance equation         [W/m^2]
!                  1: work done by tidal force
!                  2: work done by shear stress
!                  3: work done by density flux
!                  4: total work
!
!------------------------------------------------------------------------------
!*Local variables
!
INTEGER :: klo, kup
REAL :: fac, zcoord
REAL, DIMENSION(2) :: bflux, uflux, vflux
REAL, DIMENSION(nz) :: uvelatc, vvelatc


!
!1. Initialise
!-------------
!
!---parameters
wcomps = 0
IF (nodeatc(i,j).EQ.0) RETURN

!---currents
klo = MAX(1,k-1); kup = MIN(nz,k+1)
CALL Uarr_at_C(uvel(i:i+1,j,klo:kup),uvelatc(klo:kup),1,1,(/i,j,klo/),&
             & (/i+1,j,kup/),1,iarr_uvel,.FALSE.)
CALL Varr_at_C(vvel(i,j:j+1,klo:kup),vvelatc(klo:kup),1,1,(/i,j,klo/),&
            & (/i,j+1,kup/),1,iarr_vvel,.FALSE.)

!---Z-coordinates
IF (iopt_dens.GT.0) zcoord = Zcoord_var(i,j,k,'C  ',.FALSE.)

!
!2. Work done by tidal force
!---------------------------
!

IF (iopt_astro_tide.EQ.1) THEN
   wcomps(1) = density_ref*(uvelatc(k)*fxastro(i,j)+vvelatc(k)*fyastro(i,j))
ENDIF

!
!3. Work done by shear stress
!-----------------------------
!
!---lower fluxes
IF (k.EQ.1) THEN
   fac = bdragcoefatc(i,j)*SQRT(uvelatc(1)**2+vvelatc(1)**2)
   uflux(1) = fac*uvelatc(1)
   vflux(1) = fac*vvelatc(1)
ELSE
   uflux(1) = vdifcoefmom(i,j,k)*&
            & (uvelatc(k)-uvelatc(k-1))/delzatw(i,j,k)
   vflux(1) = vdifcoefmom(i,j,k)*&
            & (vvelatc(k)-vvelatc(k-1))/delzatw(i,j,k)
ENDIF

!---upper fluxes
IF (k.EQ.nz) THEN
   uflux(2) = usstresatc(i,j)
   vflux(2) = vsstresatc(i,j)
ELSE
   uflux(2) = vdifcoefmom(i,j,k+1)*&
            & (uvelatc(k+1)-uvelatc(k))/delzatw(i,j,k+1)
   vflux(2) = vdifcoefmom(i,j,k+1)*&
            & (vvelatc(k+1)-vvelatc(k))/delzatw(i,j,k+1)
ENDIF

!---work by shear stress
wcomps(2) = density_ref*(uvelatc(k)*(uflux(2)-uflux(1))+&
                       & vvelatc(k)*(vflux(2)-vflux(1)))/delzatc(i,j,k)

!
!4. Work done by buoyancy flux
!-----------------------------
!

IF (iopt_dens.GT.0) THEN

   bflux = 0.0

!
!4.1 Salinity fluxes
!-------------------
!

   IF (iopt_sal.GT.0) THEN
!     ---lower flux
      IF (k.GT.1) THEN
         bflux(1) = beta_sal(i,j,nz)*vdifcoefscal(i,j,k)*&
                  & (sal(i,j,k)-sal(i,j,k-1))/delzatw(i,j,k)
      ENDIF
!     ---upper flux
      IF (k.EQ.nz) THEN
         bflux(2) = beta_sal(i,j,nz)*ssalflux(i,j)
      ELSE
         bflux(2) = beta_sal(i,j,k)*vdifcoefscal(i,j,k+1)*&
                  & (sal(i,j,k+1)-sal(i,j,k))/delzatw(i,j,k+1)
      ENDIF
   ENDIF

!
!4.2 Temperature fluxes
!----------------------
!

   IF (iopt_temp.GT.0) THEN
!     ---lower flux
      IF (k.GT.1) THEN
         bflux(1) = bflux(1) - beta_temp(i,j,k)*vdifcoefscal(i,j,k)*&
                  & (temp(i,j,k)-temp(i,j,k-1))/delzatw(i,j,k)
      ENDIF
!     ---upper flux
      IF (k.EQ.nz) THEN
         bflux(2) = bflux(2) - beta_temp(i,j,nz)*(qrad(i,j)-qnonsol(i,j))&
                  & /(density_ref*specheat)
      ELSE
         bflux(2) = bflux(2) - beta_temp(i,j,k)*vdifcoefscal(i,j,k+1)*&
                  & (temp(i,j,k+1)-temp(i,j,k))/delzatw(i,j,k+1)
      ENDIF
   ENDIF

!
!4.3 Work by buoyancy flux
!-------------------------
!

   wcomps(3) = density_ref*zcoord*(bflux(2)-bflux(1))/delzatc(i,j,k)

ENDIF

!
!5. Total work
!-------------
!

wcomps(4) = SUM(wcomps(1:3))


RETURN

END SUBROUTINE energy_force_3d

!========================================================================

FUNCTION enstrophy_0d(nodim)
!************************************************************************
!
! *enstrophy_0d* Domain integrated enstrophy
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.11.2
!
! Description -
!
! Module calls - error_abort, error_alloc, error_limits_var, sum2_vars
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE error_routines, ONLY: error_abort, error_alloc, error_limits_var
USE paral_utilities, ONLY: sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
INTEGER, INTENT(IN) :: nodim
REAL :: enstrophy_0d

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*nodim*        INTEGER Dimension of current field
!                       = 2 => depth-mean currents
!                       = 3 => 3-D currents
!*enstrophy_0d* REAL    Enstrophy                                     [m^3/s^2]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k
REAL :: enstr0d
INTEGER, DIMENSION(4) :: nhdims
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: area, array2dc, enstr2d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: enstr3d


procname(pglev+1) = 'enstrophy_0d'
CALL log_timer_in()

!
!1. Check argument
!-----------------
!

CALL error_limits_var(nodim,'nodim',2,3)
CALL error_abort('enstrophy_0d',ierrno_arg)

!
!2. Allocate arrays
!------------------
!

ALLOCATE (area(ncloc,nrloc),STAT=errstat)
CALL error_alloc('area',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
IF (nodim.EQ.2) THEN
   ALLOCATE (enstr2d(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('enstr2d',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (nodim.EQ.3) THEN
   ALLOCATE (enstr3d(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('enstr3d',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

!
!3. Initialise
!-------------
!
!---parameters
enstrophy_0d = 0.0
nhdims = 0

!---arrays
IF (nodim.EQ.2) THEN
   enstr2d = 0.0
ELSEIF (nodim.EQ.3) THEN
   enstr3d = 0.0
ENDIF
WHERE (node2duv(1:ncloc,1:nrloc).EQ.1)
   array2dc = 0.0
   area = 1.0E-06*delxatuv(1:ncloc,1:nrloc)*delyatuv(1:ncloc,1:nrloc)
END WHERE

!
!4. 2-D case
!-----------
!

IF (nodim.EQ.2) THEN
   WHERE (node2duv(1:ncloc,1:nrloc).EQ.4)
      array2dc = (vmvel(1:ncloc,1:nrloc)-vmvel(0:ncloc-1,1:nrloc))&
               & /delxatuv(1:ncloc,1:nrloc) - &
               & (umvel(1:ncloc,1:nrloc)-umvel(1:ncloc,0:nrloc-1))&
               & /delyatuv(1:ncloc,1:nrloc)
      enstr2d = area*deptotatuv(1:ncloc,1:nrloc)*array2dc*array2dc
   END WHERE
   CALL sum2_vars(enstr2d,enstr0d,nhdims,'UV ',0,commall=.TRUE.)
   enstrophy_0d = enstr0d

!
!5. 3-D case
!-----------
!

ELSEIF (nodim.EQ.3) THEN
   k_510: DO k=1,nz
      WHERE (nodeatuv(1:ncloc,1:nrloc,k).EQ.1)
         array2dc = (vvel(1:ncloc,1:nrloc,k)-vvel(0:ncloc-1,1:nrloc,k))&
                  & /delxatuv(1:ncloc,1:nrloc) - &
                  & (uvel(1:ncloc,1:nrloc,k)-uvel(1:ncloc,0:nrloc-1,k))&
                  & /delyatuv(1:ncloc,1:nrloc)
         enstr3d(:,:,k) = area*delzatuv(1:ncloc,1:nrloc,k)*array2dc*array2dc
      END WHERE
   ENDDO k_510
   CALL sum2_vars(enstr3d,enstr0d,nhdims,'UV ',0,commall=.TRUE.)
   enstrophy_0d = enstr0d
ENDIF

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (area,array2dc)
IF (nodim.EQ.2) DEALLOCATE (enstr2d)
IF (nodim.EQ.3) DEALLOCATE (enstr3d)

CALL log_timer_out()


RETURN

END FUNCTION enstrophy_0d

!========================================================================

FUNCTION vorticity_2d(i,j,nodim)
!************************************************************************
!
! *vorticity_2d* Vertically-integrated vorticity at grid point (i,j)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.0
!
! Description -
!
! Module calls -
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, nodim
REAL :: vorticity_2d

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*i*            INTEGER X-index of E-node grid point
!*j*            INTEGER Y-index of E-node grid point
!*nodim*        INTEGER Dimension of current field
!                       = 2 => depth-mean currents
!                       = 3 => 3-D currents
!*vorticity_2d* REAL    Vorticity                                         [m/s]
!
!------------------------------------------------------------------------------
!
!1. Initialise
!-------------
!
!---parameters
vorticity_2d = 0.0
IF (node2duv(i,j).EQ.0) RETURN

!
!2. 2-D case
!-----------
!

IF (nodim.EQ.2) THEN
   vorticity_2d = deptotatuv(i,j)*((vmvel(i,j)-vmvel(i-1,j))/delxatuv(i,j)- &
                                 & (umvel(i,j)-umvel(i,j-1))/delyatuv(i,j))

!
!3. 3-D case
!-----------
!

ELSEIF (nodim.EQ.3) THEN
   vorticity_2d = SUM(delzatuv(i,j,:)*&
                   & (vvel(i,j,:)-vvel(i-1,j,:))/delxatuv(i,j) - &
                   & (uvel(i,j,:)-uvel(i,j-1,:))/delyatuv(i,j))
ENDIF


RETURN

END FUNCTION vorticity_2d

!========================================================================

FUNCTION vorticity_3d(i,j,k)
!************************************************************************
!
! *vorticity_3d* Vorticity at grid point (i,j,k)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)diagnostic_routines.f90  V2.0
!
! Description -
!
! Module calls -
!
!************************************************************************
!
USE currents
USE grid

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, k
REAL :: vorticity_3d

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*i*            INTEGER X-index of E-node grid point
!*j*            INTEGER Y-index of E-node grid point
!*k*            INTEGER Z-index of E-node grid point
!*vorticity_3d* REAL    Vorticity                                         [1/s]
!
!------------------------------------------------------------------------------
!
!1. Initialise
!-------------
!
!---parameters
vorticity_3d = 0.0

!
!2. Vorticity
!------------
!

IF (nodeatuv(i,j,k).EQ.1) THEN
   vorticity_3d = (vvel(i,j,k)-vvel(i-1,j,k))/delxatuv(i,j) - &
                & (uvel(i,j,k)-uvel(i,j-1,k))/delyatuv(i,j)
ENDIF


RETURN

END FUNCTION vorticity_3d


END MODULE diagnostic_routines
