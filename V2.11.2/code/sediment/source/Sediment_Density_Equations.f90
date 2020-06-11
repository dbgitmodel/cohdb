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
! *Sediment_Density_Equations* Sediment contributions to the main coherens
!                              code dealing with density
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Density_Equations.f90  V2.11.2
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description -
!
! Routines - baroclinic_gradient_sed_cubic, baroclinic_gradient_sed_sigma,
!            baroclinic_gradient_sed_z, buoyancy_frequency_sed,
!            equation_of_state_sed
!
!************************************************************************
!

!========================================================================

SUBROUTINE baroclinic_gradient_sed_cubic(zcoord,dzx,dzy,dzz,cdir)
!************************************************************************
!
! *baroclinic_gradient_sed_cubic* Baroclinic density gradient due to sediments
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Density_Equations.f90  V2.11.2
!
! Description - uses cubic method
!
! Calling program - baroclinic_gradient
!
! External calls - hcube_deriv, hcube_fluxes
!
! Module calls - Carr_at_U, Carr_at_UW, Carr_at_V, Carr_at_VW, error_alloc
!
!************************************************************************
!
USE currents
USE grid
USE gridpars
USE iopars
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE syspars
USE array_interp, ONLY: Carr_at_U, Carr_at_UW, Carr_at_V, Carr_at_VW
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE
!
!*Arguments
!
CHARACTER (LEN=1) :: cdir
REAL, INTENT(IN), DIMENSION(0:ncloc+1,0:nrloc+1,nz) ::  dzx, dzy, dzz
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,&
                          & 1-nhalo:nrloc+nhalo,nz) :: zcoord

!
! Name      Type    Purpose
!-------------------------------------------------------------------------------
!*zcoord*   REAL    Z cordinate at W points                                  [m]
!*dzx*      REAL    Harmonic derivative of x coordinate
!*dzy*      REAL    Harmonic derivative of y coordinate
!*dzz*      REAL    Harmonic derivative of z coordinate
!*cdir*     CHAR    Direction of the derivative ('X','Y')
!
!*******************************************************************************
!
!*Local variables
!
INTEGER :: f, k
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: mask3d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: beta, dxycvol, dzcvol, &
                                           & fcvolatvel, fcvolatw


IF (iopt_sed_dens_grad.EQ.0) RETURN

procname(pglev+1) = 'baroclinic_gradient_sed_cubic '
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (mask3d(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('mask3d',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (array2d1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2d2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (beta(ncloc,nrloc,2:nz+1),STAT=errstat)
CALL error_alloc('beta',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (dxycvol(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('dxycvol',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
ALLOCATE (dzcvol(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('dzcvol',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
ALLOCATE (fcvolatvel(0:ncloc,0:nrloc,nz),STAT=errstat)
CALL error_alloc('fcvolatvel',3,(/ncloc+1,nrloc+1,nz/),kndrtype)
ALLOCATE (fcvolatw(0:ncloc,0:nrloc,nz),STAT=errstat)
CALL error_alloc('fcvolatw',3,(/ncloc+1,nrloc+1,nz/),kndrtype)

!
!2. Mask arrays
!--------------
!

IF (cdir.EQ.'X') THEN
   mask3d = nodeatu(1:ncloc,1:nrloc,:).EQ.2
ELSEIF (cdir.EQ.'Y') THEN
   mask3d = nodeatv(1:ncloc,1:nrloc,:).EQ.2
ENDIF

!
!3. Contributions to baroclinic gradient
!---------------------------------------
!

f_300: DO f=1,nf

!  ---Z-derivative and fluxes
   CALL hcube_deriv(cvol(:,:,:,f),dzcvol,'Z')
   CALL hcube_fluxes(cvol(:,:,:,f),zcoord,dzcvol,dzz,fcvolatw,'Z')

!
!3.1 X-direction
!----------------
!

   IF (cdir.EQ.'X') THEN

!     ---interpolate arrays
      CALL Carr_at_U(beta_state_sed(0:ncloc,1:nrloc,nz,f),beta(:,:,nz+1),1,2,&
                  & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)
      CALL Carr_at_UW(beta_state_sed(0:ncloc,1:nrloc,:,f),beta(:,:,2:nz),2,&
                   & (/0,1,1/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)

!     ---X-derivative and fluxes
      CALL hcube_deriv(cvol(:,:,:,f),dxycvol,'X')
      CALL hcube_fluxes(cvol(:,:,:,f),zcoord,dxycvol,dzx,fcvolatvel,'X')

!     ---surface
      WHERE (mask3d(:,:,nz))
         array2d1 = gaccatu(1:ncloc,:)/delxatu(1:ncloc,1:nrloc)
         array2d2 = 0.5*array2d1*delzatu(1:ncloc,:,nz)*beta(:,:,nz+1)*&
                  & (cvol(1:ncloc,1:nrloc,nz,f)-cvol(0:ncloc-1,1:nrloc,nz,f))
         p3dbcgradatu(:,:,nz) = p3dbcgradatu(:,:,nz) - array2d2
      END WHERE

!     ---interior points
      k_310: DO k=nz,2,-1
         WHERE (mask3d(:,:,k))
            array2d2 = array2d2 + array2d1*beta(:,:,k)*&
               & (fcvolatw(1:ncloc,1:nrloc,k)-fcvolatw(0:ncloc-1,1:nrloc,k)&
               & +fcvolatvel(1:ncloc,1:nrloc,k-1)-fcvolatvel(1:ncloc,1:nrloc,k))
            p3dbcgradatu(:,:,k-1) = p3dbcgradatu(:,:,k-1) - array2d2
         END WHERE
      ENDDO k_310

!
!3.2 Y-direction
!---------------
!

   ELSEIF (cdir.EQ.'Y') THEN

!     ---interpolate arrays
      CALL Carr_at_V(beta_state_sed(1:ncloc,0:nrloc,nz,f),beta(:,:,nz+1),1,2,&
                  & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)
      CALL Carr_at_VW(beta_state_sed(1:ncloc,0:nrloc,:,f),beta(:,:,2:nz),2,&
                   & (/1,0,1/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)

!     ---X-derivative and fluxes
      CALL hcube_deriv(cvol(:,:,:,f),dxycvol,'Y')
      CALL hcube_fluxes(cvol(:,:,:,f),zcoord,dxycvol,dzy,fcvolatvel,'Y')

!     ---surface
      WHERE (mask3d(:,:,nz))
         array2d1 = gaccatv(:,1:nrloc)/delyatv(1:ncloc,1:nrloc)
         array2d2 = 0.5*array2d1*delzatv(:,1:nrloc,nz)*beta(:,:,nz+1)*&
                  & (cvol(1:ncloc,1:nrloc,nz,f)-cvol(1:ncloc,0:nrloc-1,nz,f))
         p3dbcgradatv(:,:,nz) = p3dbcgradatv(:,:,nz) - array2d2
      END WHERE

!     ---interior points
      k_320: DO k=nz,2,-1
         WHERE (mask3d(:,:,k))
            array2d2 = array2d2 + array2d1*beta(:,:,k)*&
               & (fcvolatw(1:ncloc,1:nrloc,k)-fcvolatw(1:ncloc,0:nrloc-1,k)&
               & +fcvolatvel(1:ncloc,1:nrloc,k-1)-fcvolatvel(1:ncloc,1:nrloc,k))
            p3dbcgradatv(:,:,k-1) = p3dbcgradatv(:,:,k-1) - array2d2
         END WHERE
      ENDDO k_320

   ENDIF

ENDDO f_300

!
!4. Deallocate
!-------------
!

DEALLOCATE (mask3d)
DEALLOCATE (array2d1,array2d2,beta,dxycvol,dzcvol,fcvolatvel,fcvolatw)

CALL log_timer_out()


RETURN

END SUBROUTINE baroclinic_gradient_sed_cubic

!========================================================================

SUBROUTINE baroclinic_gradient_sed_sigma(zcoord,cdir)
!************************************************************************
!
! *baroclinic_grad_sed_sigma* Baroclinic density gradient due to sediments
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Density_Equations.f90  V2.11.2
!
! Description - uses second order sigma method
!
! Calling program - baroclinic_gradient
!
! Module calls - Carr_at_U, Carr_at_UW, Carr_at_V, Carr_at_VW, Carr_at_W,
!                error_alloc
!
!************************************************************************
!
USE currents
USE grid
USE gridpars
USE iopars
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE syspars
USE array_interp, ONLY: Carr_at_U, Carr_at_UW, Carr_at_V, Carr_at_VW, Carr_at_W
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=1), INTENT(IN) :: cdir
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,&
                          & 1-nhalo:nrloc+nhalo,nz) :: zcoord

!
! Name      Type    Purpose
!-------------------------------------------------------------------------------
!*zcoord*   REAL    Z-coordinates at W-nodes                                 [m]
!*cdir*     CHAR    Direction of the derivative ('X','Y')
!
!*******************************************************************************
!
!*Local variables
!
INTEGER :: f, k
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: mask3d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: beta, cvolatvel, cvolatw


IF (iopt_sed_dens_grad.EQ.0) RETURN

procname(pglev+1) = 'baroclinic_gradient_sed_sigma'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (mask3d(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('mask3d',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (array2d1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2d2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (beta(ncloc,nrloc,2:nz+1),STAT=errstat)
CALL error_alloc('beta',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (cvolatvel(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('cvolatvel',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (cvolatw(0:ncloc,0:nrloc,2:nz),STAT=errstat)
CALL error_alloc('cvolatw',3,(/ncloc+1,nrloc+1,nz-1/),kndrtype)

!
!2. Mask arrays
!--------------
!

IF (cdir.EQ.'X') THEN
   mask3d = nodeatu(1:ncloc,1:nrloc,:).EQ.2
ELSEIF (cdir.EQ.'Y') THEN
   mask3d = nodeatv(1:ncloc,1:nrloc,:).EQ.2
ENDIF

!
!3. Contributions to baroclinic gradient
!---------------------------------------
!

f_300: DO f=1,nf

!
!3.1 Interpolate total concentration at W-nodes
!---------------------------------------------
!

   CALL Carr_at_W(cvol(0:ncloc,0:nrloc,:,f),cvolatw,(/0,0,1/),&
               & (/ncloc,nrloc,nz/),1,iarr_cvol,.TRUE.)

!
!3.2 X-direction
!----------------
!

   IF (cdir.EQ.'X') THEN

!     ---interpolate arrays
      CALL Carr_at_U(beta_state_sed(0:ncloc,1:nrloc,nz,f),beta(:,:,nz+1),1,2,&
                  & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)
      CALL Carr_at_UW(beta_state_sed(0:ncloc,1:nrloc,:,f),beta(:,:,2:nz),2,&
                   & (/0,1,1/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)
      CALL Carr_at_U(cvol(0:ncloc,1:nrloc,:,f),cvolatvel,1,2,(/0,1,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_cvol,.TRUE.)

!     ---surface
      WHERE (mask3d(:,:,nz))
         array2d1 = gaccatu(1:ncloc,:)/delxatu(1:ncloc,1:nrloc)
         array2d2 = 0.5*array2d1*delzatu(1:ncloc,:,nz)*beta(:,:,nz+1)*&
                  & (cvol(1:ncloc,1:nrloc,nz,f)-cvol(0:ncloc-1,1:nrloc,nz,f))
         p3dbcgradatu(:,:,nz) = p3dbcgradatu(:,:,nz) - array2d2
      END WHERE

!     ---interior points
      k_320: DO k=nz,2,-1
         WHERE (mask3d(:,:,k))
            array2d2 = array2d2 + array2d1*beta(:,:,k)*&
                & ((cvolatw(1:ncloc,1:nrloc,k)-cvolatw(0:ncloc-1,1:nrloc,k))*&
                & delzatuw(1:ncloc,:,k)-&
                & (cvolatvel(:,:,k)-cvolatvel(:,:,k-1))*&
                & (zcoord(1:ncloc,1:nrloc,k)-zcoord(0:ncloc-1,1:nrloc,k)))
            p3dbcgradatu(:,:,k-1) = p3dbcgradatu(:,:,k-1) - array2d2
         END WHERE
      ENDDO k_320

!
!3.3 Y-direction
!---------------
!

   ELSEIF (cdir.EQ.'Y') THEN

!     ---interpolate arrays
      CALL Carr_at_V(beta_state_sed(1:ncloc,0:nrloc,nz,f),beta(:,:,nz+1),1,2,&
                  & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)
      CALL Carr_at_VW(beta_state_sed(1:ncloc,0:nrloc,:,f),beta(:,:,2:nz),2,&
                   & (/1,0,1/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)
      CALL Carr_at_V(cvol(1:ncloc,0:nrloc,:,f),cvolatvel,1,2,(/1,0,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_cvol,.TRUE.)

!     ---surface
      WHERE (mask3d(:,:,nz))
         array2d1 = gaccatv(:,1:nrloc)/delyatv(1:ncloc,1:nrloc)
         array2d2 = 0.5*array2d1*delzatv(:,1:nrloc,nz)*beta(:,:,nz+1)*&
                  & (cvol(1:ncloc,1:nrloc,nz,f)-cvol(1:ncloc,0:nrloc-1,nz,f))
         p3dbcgradatv(:,:,nz) = p3dbcgradatv(:,:,nz) - array2d2
      END WHERE

!     ---interior points
      k_330: DO k=nz,2,-1
         WHERE (mask3d(:,:,k))
            array2d2 = array2d2 + array2d1*beta(:,:,k)*&
                & ((cvolatw(1:ncloc,1:nrloc,k)-cvolatw(1:ncloc,0:nrloc-1,k))*&
                & delzatvw(:,1:nrloc,k)-&
                & (cvolatvel(:,:,k)-cvolatvel(:,:,k-1))*&
                & (zcoord(1:ncloc,1:nrloc,k)-zcoord(1:ncloc,0:nrloc-1,k)))
            p3dbcgradatv(:,:,k-1) = p3dbcgradatv(:,:,k-1) - array2d2
         END WHERE
      ENDDO k_330

   ENDIF

ENDDO f_300

!
!4. Deallocate
!-------------
!

DEALLOCATE (mask3d) 
DEALLOCATE (array2d1,array2d2,beta,cvolatvel,cvolatw)

CALL log_timer_out()


RETURN

END SUBROUTINE baroclinic_gradient_sed_sigma

!========================================================================

SUBROUTINE baroclinic_gradient_sed_z(sigint,kint,cdir)
!************************************************************************
!
! *baroclinic_grad_sed_z* Baroclinic density gradient due to sediments
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Density_Equations.f90  V2.11.2
!
! Description - uses Z-level method
!
! Calling program - baroclinic_gradient
!
! Module calls - Carr_at_U, Carr_at_UW, Carr_at_V, Carr_at_VW, error_alloc
!
!************************************************************************
!
USe currents
USE grid
USE gridpars
USE iopars
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE syspars
USE array_interp, ONLY: Carr_at_U, Carr_at_UW, Carr_at_V, Carr_at_VW
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE
!
!*Arguments
!
CHARACTER (LEN=1), INTENT(IN) :: cdir
INTEGER, INTENT(IN), DIMENSION(ncloc,nrloc,2:nz+1,2) :: kint
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,2:nz+1,2) :: sigint

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*sigint*   REAL    Interpolated sigma coordinates
!*kint*     INTEGER Counter for the interpolated coordinates
!*cdir*     CHAR    Direction of the derivative ('X','Y')
!
!******************************************************************************
!
!*Local variables
!
INTEGER :: f, i, ii, kc, j, jj, k, l
REAL :: sigk
REAL, DIMENSION(nz+1,2) :: cvolint
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: mask3d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: beta, dhcvol


IF (iopt_sed_dens_grad.EQ.0) RETURN

procname(pglev+1) = 'baroclinic_gradient_sed_z'
CALL log_timer_in()

ALLOCATE (mask3d(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('mask3d',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (array2d1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2d2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (beta(ncloc,nrloc,2:nz+1),STAT=errstat)
CALL error_alloc('beta',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (dhcvol(ncloc,nrloc,2:nz+1),STAT=errstat)
CALL error_alloc('dhcvol',3,(/ncloc,nrloc,nz/),kndrtype)

!
!2. Mask arrays
!--------------
!

IF (cdir.EQ.'X') THEN
   mask3d = nodeatu(1:ncloc,1:nrloc,:).EQ.2
ELSEIF (cdir.EQ.'Y') THEN
   mask3d = nodeatv(1:ncloc,1:nrloc,:).EQ.2
ENDIF

!
!3. Contributions to baroclinic gradient
!---------------------------------------
!

f_300: DO f=1,nf

!
!3.1 X-direction
!---------------
!

   IF (cdir.EQ.'X') THEN

!     ---horizontal concentration difference
      j_311: DO j=1,nrloc
      i_311: DO i=1,ncloc
      k_311: DO k=2,nz+1
         IF ((k.EQ.(nz+1).AND.mask3d(i,j,nz)).OR.nodeatuw(i,j,k).EQ.2) THEN
            ii_3111: DO ii=i-1,i
               l = ii-i+2
               kc = kint(i,j,k,l)
               IF (kc.EQ.0) THEN
                  cvolint(k,l) = cvol(ii,j,1,f)
               ELSEIF (kc.EQ.nz) THEN
                  cvolint(k,l) = cvol(ii,j,nz,f)
               ELSE
                  sigk = sigint(i,j,k,l)
                  cvolint(k,l) = (1.0-sigk)*cvol(ii,j,kc,f)+&
                                & sigk*cvol(ii,j,kc+1,f)
               ENDIF
            ENDDO ii_3111
            dhcvol(i,j,k) = cvolint(k,2) - cvolint(k,1)
         ENDIF
      ENDDO k_311
      ENDDO i_311
      ENDDO j_311

!     ---interpolate arrays
      CALL Carr_at_U(beta_state_sed(0:ncloc,1:nrloc,nz,f),beta(:,:,nz+1),1,2,&
                  & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)
      CALL Carr_at_UW(beta_state_sed(0:ncloc,1:nrloc,:,f),beta(:,:,2:nz),2,&
                   & (/0,1,1/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)

!     ---surface
      WHERE (mask3d(:,:,nz))
         array2d1 = gaccatu(1:ncloc,:)/delxatu(1:ncloc,1:nrloc)
         array2d2 = 0.5*array2d1*delzatu(1:ncloc,:,nz)*beta(:,:,nz+1)*&
                  & dhcvol(:,:,nz+1)
         p3dbcgradatu(:,:,nz) = p3dbcgradatu(:,:,nz) - array2d2
      END WHERE
      
!     ---interior points
      k_312: DO k=nz,2,-1
         WHERE (mask3d(:,:,k))
            array2d2 = array2d2 + array2d1*beta(:,:,k)*dhcvol(:,:,k)*&
                     & delzatuw(1:ncloc,:,k)
            p3dbcgradatu(:,:,k-1) = p3dbcgradatu(:,:,k-1) - array2d2
         END WHERE
      ENDDO k_312

!
!3.2 Y-direction
!---------------
!

   ELSEIF (cdir.EQ.'Y') THEN

!     ---horizontal concentration difference
      j_321: DO j=1,nrloc
      i_321: DO i=1,ncloc
      k_321: DO k=2,nz+1
         IF ((k.EQ.(nz+1).AND.mask3d(i,j,nz)).OR.nodeatvw(i,j,k).EQ.2) THEN
            jj_3211: DO jj=j-1,j
               l = jj-j+2
               kc = kint(i,j,k,l)
               IF (kc.EQ.0) THEN
                  cvolint(k,l) = cvol(i,jj,1,f)
               ELSEIF (kc.EQ.nz) THEN
                  cvolint(k,l) = cvol(i,jj,nz,f)
               ELSE
                  sigk = sigint(i,j,k,l)
                  cvolint(k,l) = (1.0-sigk)*cvol(i,jj,kc,f)+&
                               & sigk*cvol(i,jj,kc+1,f)
               ENDIF
            ENDDO jj_3211
            dhcvol(i,j,k) = cvolint(k,2) - cvolint(k,1)
         ENDIF
      ENDDO k_321
      ENDDO i_321
      ENDDO j_321

!     ---interpolate arrays
      CALL Carr_at_V(beta_state_sed(1:ncloc,0:nrloc,nz,f),beta(:,:,nz+1),1,2,&
                  & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)
      CALL Carr_at_VW(beta_state_sed(1:ncloc,0:nrloc,:,f),beta(:,:,2:nz),2,&
                   & (/1,0,1/),(/ncloc,nrloc,nz/),1,iarr_beta_state_sed,.TRUE.)

!     ---surface
      WHERE (mask3d(:,:,nz))
         array2d1 = gaccatv(:,1:nrloc)/delyatv(1:ncloc,1:nrloc)
         array2d2 = 0.5*array2d1*delzatv(:,1:nrloc,nz)*beta(:,:,nz+1)*&
                  & dhcvol(:,:,nz+1)
         p3dbcgradatv(:,:,nz) = p3dbcgradatv(:,:,nz) - array2d2
      END WHERE
      
!     ---interior points
      k_322: DO k=nz,2,-1
         WHERE (mask3d(:,:,k))
            array2d2 = array2d2 + array2d1*beta(:,:,k)*dhcvol(:,:,k)*&
                     & delzatvw(:,1:nrloc,k)
            p3dbcgradatv(:,:,k-1) = p3dbcgradatv(:,:,k-1) - array2d2
         END WHERE
      ENDDO k_322

   ENDIF

ENDDO f_300

!
!4. Deallocate
!-------------
!

DEALLOCATE (mask3d)
DEALLOCATE (array2d1,array2d2,beta,dhcvol)

CALL log_timer_out()


RETURN

END SUBROUTINE baroclinic_gradient_sed_z

!========================================================================

SUBROUTINE buoyancy_frequency_sed(bgrad)
!************************************************************************
!
! *buoyancy_frequency_sed* Squared buoyancy frequency due to sediments
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Density_Equations.f90  V2.11.2
!
! Description - defined as vertical gradient of buoyancy
!             - uses interpolation of neighbouring cells
!
! Calling program - buoyancy_frequency
!
! Module calls - Carr_at_W, error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE syspars
USE array_interp, ONLY: Carr_at_W
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
REAL, INTENT(INOUT), DIMENSION(0:ncloc+1,0:nrloc+1,2:nz) :: bgrad

!
! Name      Type    Purpose
!-------------------------------------------------------------------------------
!*bgrad*       REAL    Non-averaged buoyancy gradient                    [1/s^2]
!
!*******************************************************************************
!
!*Local variables
!
INTEGER :: f, k
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: beta


IF (iopt_sed_dens_grad.EQ.0) RETURN

procname(pglev+1) = 'buoyancy_frequency_sediment'
CALL log_timer_in()

!
!1. Allocate
!-----------
!

ALLOCATE (beta(0:ncloc+1,0:nrloc+1,2:nz),STAT=errstat)
CALL error_alloc('beta',3,(/ncloc+2,nrloc+2,nz-1/),kndrtype)

!
!2. Add sediment contribution
!----------------------------
!

f_210: DO f=1,nf
   CALL Carr_at_W(beta_state_sed(:,:,:,f),beta,(/0,0,1/),&
               & (/ncloc+1,nrloc+1,nz/),1,iarr_beta_state_sed,.TRUE.)
   k_211: DO k=2,nz
      WHERE (nodeatc(0:ncloc+1,0:nrloc+1).GT.0)
         bgrad(:,:,k) = bgrad(:,:,k) - gaccatc*beta(:,:,k)*&
                     & (cvol(0:ncloc+1,0:nrloc+1,k,f)-&
                      & cvol(0:ncloc+1,0:nrloc+1,k-1,f))/delzatw(:,:,k)
      END WHERE
   ENDDO k_211
ENDDO f_210

!
!2.Deallocate
!------------
!

DEALLOCATE (beta)

CALL log_timer_out()


RETURN

END SUBROUTINE buoyancy_frequency_sed

!======================================================================

SUBROUTINE equation_of_state_sed
!************************************************************************
!
! *equation_of_state_sed* Density and expansion coefficients of a
!                         water-sediment mixture
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Density_Equations.f90  V2.11.1
!
! Description - 
!
! Calling program - equation_of_state
!
! External calls - equation_of_state_floc
!
! Module calls -
!
!************************************************************************
!
USE density
USE grid
USE gridpars
USE iopars
USE physpars
USE sedarrays
USE sedids
USE sedpars
USE switches
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: k, f
INTEGER, DIMENSION(4) :: nhexch
LOGICAL, DIMENSION(0:ncloc+1,0:nrloc+1) :: maskatc
REAL, DIMENSION(0:ncloc+1,0:nrloc+1) :: array2dc

   
procname(pglev+1) = 'equation_of_state_sed'
CALL log_timer_in()

!
!1. Initialise
!-------------
!
!---mask array
maskatc = nodeatc(0:ncloc+1,0:nrloc+1).GT.0

!---"pure" water density
densw = dens

!
!2. Add sediment contributions
!-----------------------------
!
!2.1 Without flocculation
!------------------------
!

IF (iopt_sed.EQ.1) THEN

   k_210: DO k=1,nz

!
!2.1.1 Density
!-------------
!

      array2dc = 0.0
      f_211: DO f=1,nf
         WHERE (maskatc)
            array2dc = array2dc + cvol(0:ncloc+1,0:nrloc+1,k,f)*rhos(f)
         END WHERE
      ENDDO f_211
      IF (iopt_dens.EQ.0) THEN
         WHERE (maskatc)
            dens(:,:,k) = (1.0-ctot(:,:,k))*density_ref + array2dc
         END WHERE
      ELSE
         WHERE (maskatc)
            dens(:,:,k) = (1.0-ctot(:,:,k))*densw(:,:,k) + array2dc
         END WHERE
      ENDIF

!
!2.1.2 Expansion coefficients
!----------------------------
!

      f_212: DO f=1,nf
         IF (iopt_dens.EQ.0) THEN
            WHERE (maskatc)
               beta_state_sed(:,:,k,f) = (rhos(f)-density_ref)/dens(:,:,k)
            END WHERE
         ELSE
            WHERE (maskatc)
               beta_state_sed(:,:,k,f) = (rhos(f)-densw(:,:,k))/dens(:,:,k)
            END WHERE
         ENDIF
      ENDDO f_212
      
   ENDDO k_210
   
!
!2.2 With flocculation
!---------------------
!

ELSEIF (iopt_sed.EQ.2) THEN

   k_220: DO k=1,nz

!
!2.2.1 Floc density
!------------------
!
 
   WHERE (maskatc_int)
      array2dc(1:ncloc,1:nrloc) = floc_nc(1:ncloc,1:nrloc,k)**(1.0-3.0/nfrdim)
      floc_dens(1:ncloc,1:nrloc,k) = densw(1:ncloc,1:nrloc,k)*&
                                   & (1.0-array2dc(1:ncloc,1:nrloc)) + &
                                   & rhos(1)*array2dc(1:ncloc,1:nrloc)
   END WHERE

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      nhexch = 1
      CALL exchange_mod(floc_dens,(/0,0,1/),nhexch,iarr_floc_dens)
   ENDIF
   
!
!2.2.2 Density      
!-------------
!      
     
      array2dc = 0.0
      WHERE (maskatc)
         array2dc = rhos(1)*(cvol(0:ncloc+1,0:nrloc+1,k,1) + &
                           & cvol(0:ncloc+1,0:nrloc+1,k,2)*floc_nc(:,:,k))
      END WHERE
      IF (iopt_dens.EQ.0) THEN
         WHERE (maskatc)
            dens(:,:,k) = (1.0-ctot(:,:,k))*density_ref + array2dc
         END WHERE
      ELSE
         WHERE (maskatc)
            dens(:,:,k) = (1.0-ctot(:,:,k))*densw(:,:,k) + array2dc
         END WHERE
      ENDIF
      
!
!2.2.3 Expansion coefficients
!----------------------------
!

      IF (iopt_dens.EQ.0) THEN
         WHERE (maskatc)
            beta_state_sed(:,:,k,1) = (rhos(1)-density_ref)/dens(:,:,k)
            beta_state_sed(:,:,k,2) = (floc_dens(:,:,k)-density_ref)/dens(:,:,k)
         END WHERE
      ELSE
         WHERE (maskatc)
            beta_state_sed(:,:,k,1) = (rhos(1)-densw(:,:,k))/dens(:,:,k)
            beta_state_sed(:,:,k,2) = (floc_dens(:,:,k)-densw(:,:,k))&
                                    & /dens(:,:,k)
         END WHERE
      ENDIF
      
   ENDDO k_220
   
ENDIF
   
CALL log_timer_out()


RETURN

END SUBROUTINE equation_of_state_sed
