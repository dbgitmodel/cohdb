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
! *Diffusion_Coefficients* Diffusion coeffcients and slope factors for
!                          geophysical and isoneutral diffusion
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
! Reference -
!
! Routines - horizontal_diff_coeefs, geoslopes, isoslopes,
!            kinematic_viscosity, vertical_diff_coefs, vertical_diff_rot
!
!************************************************************************
!

!========================================================================

SUBROUTINE geoslopes
!************************************************************************
!
! *geoslopes* Geographical slope factors 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Coefficients.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
! Reference -
!
! Calling program -
!
! External calls -
!
! Module calls - error_alloc, Carr_at_W, Zcoord_arr
!
!************************************************************************
!
USE diffusion
USE grid
USE gridpars
USE iopars
USE physpars
USE switches
USE syspars
USE timepars
USE array_interp, ONLY: Carr_at_W
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: Zcoord_arr
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, k
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: delta
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: zcoord


procname(pglev+1) = 'geoslopes'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatu(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc+1,nrloc/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc+1/),kndlog)
ALLOCATE (zcoord(0:ncloc+1,0:nrloc+1,nz+1),STAT=errstat)
CALL error_alloc('zcoord',3,(/ncloc+2,nrloc+2,nz+1/),kndrtype)
IF (iopt_hdif_scal.EQ.2) THEN
   ALLOCATE (delta(ncloc+1,nrloc+1),STAT=errstat)
   CALL error_alloc('delta',2,(/ncloc+1,nrloc+1/),kndrtype)
ENDIF

!
!2. Initialise
!-------------
!   
!---mask arrays
maskatu = node2du(1:ncloc+1,1:nrloc).EQ.2
maskatv = node2dv(1:ncloc,1:nrloc+1).EQ.2
   
!
!2. Slopes in the X-direction
!----------------------------
!
!---Z-coordinates
CALL Zcoord_arr(zcoord(0:ncloc+1,1:nrloc,1:nz),(/0,1,1/),&
             & (/ncloc+1,nrloc,nz/),'C  ',.FALSE.)

!---U-nodes
k_210: DO k=1,nz
   WHERE (maskatu)
      xslopeatu_geo(1:ncloc+1,:,k) = (zcoord(1:ncloc+1,1:nrloc,k)-&
                                    & zcoord(0:ncloc,1:nrloc,k))&
                                    & /delxatu(1:ncloc+1,1:nrloc)
   END WHERE
ENDDO k_210

!---W-nodes
IF (iopt_hdif_scal.EQ.2) THEN
   CALL Zcoord_arr(zcoord(1:ncloc+1,1:nrloc,2:nz),(/1,1,2/),&
                & (/ncloc+1,nrloc,nz/),'UW ',.FALSE.)
   k_220: DO k=2,nz
   i_220: DO i=1,ncloc
   j_220: DO j=1,nrloc   
      IF (maskatu(i,j).AND.maskatu(i+1,j)) THEN
         xslopeatw_geo(i,j,k) = (zcoord(i+1,j,k)-zcoord(i,j,k))/delxatc(i,j)
      ELSEIF (maskatu(i,j)) THEN
         xslopeatw_geo(i,j,k) = 0.5*(xslopeatu_geo(i,j,k)+&
                                   & xslopeatu_geo(i,j,k-1))
      ELSEIF (maskatu(i+1,j)) THEN
         xslopeatw_geo(i,j,k) = 0.5*(xslopeatu_geo(i+1,j,k)+&
                                   & xslopeatu_geo(i+1,j,k-1))
      ELSE
         xslopeatw_geo(i,j,k) = 0.0
      ENDIF
   ENDDO j_220
   ENDDO i_220
   ENDDO k_220
ENDIF

!
!3. Slopes in the Y-direction
!------------------------
!
!---Z-coordinates
CALL Zcoord_arr(zcoord(1:ncloc,0:nrloc+1,1:nz),(/0,1,1/),&
             & (/ncloc,nrloc+1,nz/),'C  ',.FALSE.)

!---V-nodes
k_310: DO k=1,nz
   WHERE (maskatv)
      yslopeatv_geo(:,1:nrloc+1,k) = (zcoord(1:ncloc,1:nrloc+1,k)-&
                                    & zcoord(1:ncloc,0:nrloc,k))&
                                    & /delyatv(1:ncloc,1:nrloc+1)
   END WHERE
ENDDO k_310

!---W-nodes
IF (iopt_hdif_scal.EQ.2) THEN
   CALL Zcoord_arr(zcoord(1:ncloc,1:nrloc+1,2:nz),(/1,1,2/),&
                & (/ncloc,nrloc+1,nz/),'VW ',.FALSE.)
   k_320: DO k=2,nz
   i_320: DO i=1,ncloc
   j_320: DO j=1,nrloc   
      IF (maskatv(i,j).AND.maskatv(i,j+1)) THEN
         yslopeatw_geo(i,j,k) = (zcoord(i,j+1,k)-zcoord(i,j,k))/delyatc(i,j)
      ELSEIF (maskatv(i,j)) THEN
         yslopeatw_geo(i,j,k) = 0.5*(yslopeatv_geo(i,j,k)+&
                                   & yslopeatv_geo(i,j,k-1))
      ELSEIF (maskatv(i,j+1)) THEN
         yslopeatw_geo(i,j,k) = 0.5*(yslopeatv_geo(i,j+1,k)+&
                                   & yslopeatv_geo(i,j+1,k-1))
      ELSE
         yslopeatw_geo(i,j,k) = 0.0
      ENDIF
   ENDDO j_320
   ENDDO i_320
   ENDDO k_320
ENDIF
   
IF (iopt_hdif_scal.EQ.3) GOTO 1000

!
!4. Contibution to vertical diffusion coefficient
!------------------------------------------------
!

k_410: DO k=2,nz
   WHERE (maskatc_int)
      vdifcoefscal_rot(:,:,k) = hdifscal_fac*hdifcoef3datw(:,:,k)*&
                             & (xslopeatw_geo(:,:,k)**2+yslopeatw_geo(:,:,k)**2)
   END WHERE
ENDDO k_410

!
!5. Limiting factors for stability
!---------------------------------
!
!5.1 U-nodes
!-----------
!

k_510: DO k=1,nz
   
!  ---"delta" factor"
   WHERE (maskatu)
      delta(:,1:nrloc) = 0.5*delxatu(1:ncloc+1,1:nrloc)*delzatu(1:ncloc+1,:,k)&
                       & /(hdifscal_fac*hdifcoef3datu(:,:,k)*delt3d)
      xslopeatu_geo(:,:,k) = MIN(delta(:,1:nrloc),xslopeatu_geo(:,:,k))
   END WHERE

ENDDO k_510
   
!
!5.2 V-nodes
!-----------
!

k_520: DO k=1,nz
   
   WHERE (maskatv)
      delta(1:ncloc,:) = delyatv(1:ncloc,1:nrloc+1)*delzatv(:,1:nrloc+1,k)&
                       & /(hdifscal_fac*hdifcoef3datv(:,:,k)*delt3d)
      yslopeatv_geo(:,:,k) = 0.5*MIN(delta(1:ncloc,:),yslopeatv_geo(:,:,k))
   END WHERE

ENDDO k_520

!
!5.3 W-nodes
!-----------
!

k_530: DO k=2,nz

!
!5.3.1 X-direction
!-----------------
!   

   WHERE (maskatc_int)
      delta(1:ncloc,1:nrloc) = delxatc(1:ncloc,1:nrloc)*&
                       & delzatc(1:ncloc,1:nrloc,k)&
                       & /(hdifscal_fac*hdifcoef3datc(1:ncloc,1:nrloc,k)*delt3d)
      xslopeatw_geo(:,:,k) = 0.5*MIN(delta(1:ncloc,1:nrloc),&
                           & xslopeatw_geo(:,:,k))
   END WHERE


!
!5.3.2 Y-direction
!-----------------
!   

   WHERE (maskatc_int)
      delta(1:ncloc,1:nrloc) = delyatc(1:ncloc,1:nrloc)*&
                       & delzatc(1:ncloc,1:nrloc,k)&
                       & /(hdifscal_fac*hdifcoef3datc(1:ncloc,1:nrloc,k)*delt3d)
      yslopeatw_geo(:,:,k) = 0.5*MIN(delta(1:ncloc,1:nrloc),&
                           & yslopeatw_geo(:,:,k))
   END WHERE
   
ENDDO k_530

!
!6. Deallocate
!-------------
!

1000 CONTINUE

DEALLOCATE (maskatu,maskatv,zcoord)
IF (iopt_hdif_scal.EQ.2) DEALLOCATE (delta)
   
CALL log_timer_out()


RETURN

END SUBROUTINE geoslopes
   
!========================================================================

SUBROUTINE horizontal_diff_coefs
!************************************************************************
!
! *horizontal_diff_coefs* Horizontal diffusion coefficients
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Coefficients.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
! Reference -
!
  ! Calling program - coherens_main, initialise_model
!
! External calls - geoslopes, isoslopes
!
! Module calls - error_alloc, exchange_mod, Carr_at_U, Carr_at_UV, Carr_at_V,
!                Carr_at_W, UVarr_at_C, UVarr_at_U, UVarr_at_V,
!                warning_reset_var
!
!************************************************************************
!
USE currents
USE depths
USE diffusion
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE physpars
USE switches
USE timepars
USE array_interp, ONLY: Carr_at_U, Carr_at_UV, Carr_at_V, Carr_at_W, &
                      & UVarr_at_C, UVarr_at_U, UVarr_at_V
USE error_routines, ONLY: error_alloc, warning_reset_var
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: info = .FALSE.
INTEGER :: i, ic, ii, iiloc, j, jc, jj, jjloc, k, npcc
REAL :: cflmax, xfac
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhdims
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu, maskatv, maskatuv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc1, array2dc2, array2du1, &
                              & array2du2, array2du3, array2dv1, array2dv2, &
                              & array2dv3, array2duv1, array2duv2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: DT_tension, DS_strain

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!
!*DT_tension*  REAL     Horizontal tension                              [1/s]
!*DS_strain*   REAL     Horizontal shearing strain                      [1/s]
!
!------------------------------------------------------------------------------
!

procname(pglev+1) = 'horizontal_diff_coefs'
CALL log_timer_in(npcc)


!
!1. Time-independent arrays on first call
!----------------------------------------
!

IF (nt.EQ.0) THEN

!
!1.1 Uniform (3-D) diffusion coefficients
!----------------------------------------
!   

   IF (iopt_hdif_coef.EQ.1) THEN

!     ---reset to maximum value by stability constraint
      cflmax = MIN(delxatc_min,delyatc_min)**2/delt3d
      IF (hdifmom_cst.GT.cflmax) THEN
         CALL warning_reset_var(hdifmom_cst,'hdifmom_cst',cflmax)
      ENDIF
      IF (hdifscal_cst.GT.cflmax) THEN
         CALL warning_reset_var(hdifscal_cst,'hdifscal_cst',cflmax)
      ENDIF
   
!     ---diffusion coeffient arrays
      WHERE (nodeatc(0:ncloc,0:nrloc).GT.0)
         hdifcoef2datc = deptotatc(0:ncloc,0:nrloc)
      END WHERE
      WHERE (node2duv(1:ncloc+1,1:nrloc+1).EQ.1)
         hdifcoef2datuv = deptotatuv
      END WHERE
      IF (iopt_grid_nodim.EQ.3) THEN
         k_111: DO k=1,nz
            WHERE (nodeatc(0:ncloc,0:nrloc).GT.0)
               hdifcoef3datc(:,:,k) = 1.0
            END WHERE
         ENDDO k_111
         IF (iopt_hdif_scal.GE.2) THEN
            k_112: DO k=2,nz
               WHERE (maskatc_int)
                  hdifcoef3datw(:,:,k) = 1.0
               END WHERE
            ENDDO k_112
         ENDIF
         WHERE (nodeatu(1:ncloc+1,1:nrloc,:).GT.1)
            hdifcoef3datu = 1.0
         END WHERE
         WHERE (nodeatv(1:ncloc,1:nrloc+1,:).GT.1) 
            hdifcoef3datv = 1.0
         END WHERE
         WHERE (nodeatuv(1:ncloc+1,1:nrloc+1,:).EQ.1)
            hdifcoef3datuv = 1.0
         END WHERE
      ENDIF

   ENDIF


!
!1.2 Slope factors z-level diffusion
!-----------------------------------
!

   IF (iopt_hdif_scal.GT.1) CALL geoslopes
   IF (iopt_hdif_scal.EQ.3) CALL isoslopes

   GOTO 1000
   
ENDIF

!
!2. Allocate arrays
!------------------
!

IF (iopt_hdif_coef.EQ.2) THEN

   ALLOCATE (maskatc(0:ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('maskatc',2,(/ncloc+1,nrloc+1/),kndlog)
   ALLOCATE (maskatu(0:ncloc+1,0:nrloc,nz),STAT=errstat)
   CALL error_alloc('maskatu',3,(/ncloc+2,nrloc+1,nz/),kndlog)
   ALLOCATE (maskatv(0:ncloc,0:nrloc+1,nz),STAT=errstat)
   CALL error_alloc('maskatv',3,(/ncloc+1,nrloc+2,nz/),kndlog)
   ALLOCATE (maskatuv(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('maskatuv',3,(/ncloc,nrloc,nz/),kndlog)
   ALLOCATE (DT_tension(0:ncloc,0:nrloc,nz),STAT=errstat)
   CALL error_alloc('DT_tension',3,(/ncloc+1,nrloc+1,nz/),kndrtype)
   ALLOCATE (DS_strain(ncloc+1,nrloc+1,nz),STAT=errstat)
   CALL error_alloc('DT_strain',3,(/ncloc+1,nrloc+1,nz/),kndrtype)
   ALLOCATE (array2dc1(0:ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('array2dc1',2,(/ncloc+1,nrloc+1/),kndrtype)
   ALLOCATE (array2dc2(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc2',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (array2du1(0:ncloc+1,0:nrloc),STAT=errstat)
   CALL error_alloc('array2du1',2,(/ncloc+2,nrloc+1/),kndrtype)
   ALLOCATE (array2du2(ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('array2du2',2,(/ncloc,nrloc+1/),kndrtype)
   ALLOCATE (array2dv1(0:ncloc,0:nrloc+1),STAT=errstat)
   CALL error_alloc('array2dv1',2,(/ncloc+1,nrloc+2/),kndrtype)
   ALLOCATE (array2dv2(0:ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dv2',2,(/ncloc+1,nrloc/),kndrtype)
   ALLOCATE (array2duv1(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2duv1',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (array2duv2(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2duv2',2,(/ncloc,nrloc/),kndrtype)
   IF (iopt_grid_nodim.EQ.3) THEN
      ALLOCATE (array2du3(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('array2du3',2,(/ncloc,nrloc/),kndrtype)
      ALLOCATE (array2dv3(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('array2dv3',2,(/ncloc,nrloc/),kndrtype)
   ENDIF

ENDIF
   
!
!3. Uniform coefficients
!-----------------------
!

IF (iopt_hdif_coef.EQ.1) THEN

   WHERE (nodeatc(0:ncloc,0:nrloc).GT.0)
      hdifcoef2datc = deptotatc(0:ncloc,0:nrloc)
   END WHERE
   WHERE (node2duv(1:ncloc+1,1:nrloc+1).EQ.1)
      hdifcoef2datuv = deptotatuv
   END WHERE

   IF (iopt_hdif_scal.EQ.3) CALL isoslopes
   
   GOTO 1000

ENDIF

!
!4. Initialise
!-------------
!
!---mask arrays
maskatc = nodeatc(0:ncloc,0:nrloc).GT.0
maskatu = nodeatu(0:ncloc+1,0:nrloc,:).GT.1
maskatv = nodeatv(0:ncloc,0:nrloc+1,:).GT.1
maskatuv = nodeatuv(1:ncloc,1:nrloc,:).GT.0

!---horizontal tension
k_411: DO k=1,nz
   WHERE (.NOT.maskatc)
      DT_tension(:,:,k) = 0.0
   END WHERE
ENDDO k_411

!---horizontal shearing strain
DS_strain = 0.0

!
!5. 2-D grid
!-----------
!

IF (iopt_grid_nodim.EQ.2) THEN

!
!5.1 Horizontal tension
!----------------------
!

   WHERE (maskatc)
      array2dc1 = delyatc(0:ncloc,0:nrloc)/delxatc(0:ncloc,0:nrloc)
   END WHERE
   WHERE (maskatu(:,:,nz)) 
      array2du1 = umvel(0:ncloc+1,0:nrloc)/delyatu(0:ncloc+1,0:nrloc)
   ELSEWHERE
      array2du1 = 0.0
   END WHERE
   WHERE (maskatv(:,:,nz))
      array2dv1 = vmvel(0:ncloc,0:nrloc+1)/delxatv(0:ncloc,0:nrloc+1)
   ELSEWHERE
      array2dv1 = 0.0
   END WHERE
   WHERE (maskatc)
      DT_tension(:,:,1) = array2dc1*&
                      & (array2du1(1:ncloc+1,:)-array2du1(0:ncloc,:)) - &
                      & (array2dv1(:,1:nrloc+1)-array2dv1(:,0:nrloc))/array2dc1
   END WHERE

!
!5.2 Horizontal shearing strain
!------------------------------
!

   WHERE (maskatuv(:,:,nz))
      array2duv1 = delxatuv(1:ncloc,1:nrloc)/delyatuv(1:ncloc,1:nrloc)
   END WHERE
   WHERE (maskatu(1:ncloc,:,nz))
      array2du2 = umvel(1:ncloc,0:nrloc)/delxatu(1:ncloc,0:nrloc)
   ELSEWHERE
      array2du2 = 0.0
   END WHERE
   WHERE (maskatv(:,1:nrloc,nz))
      array2dv2 = vmvel(0:ncloc,1:nrloc)/delyatv(0:ncloc,1:nrloc)
   ELSEWHERE
      array2dv2 = 0.0
   END WHERE
   WHERE (maskatuv(:,:,nz))
      DS_strain(1:ncloc,1:nrloc,1) = array2duv1*&
                     & (array2du2(:,1:nrloc)-array2du2(:,0:nrloc-1)) + &
                     & (array2dv2(1:ncloc,:)-array2dv2(0:ncloc-1,:))/array2duv1
   END WHERE

!  ---X-open boundaries
   iiloc_521: DO iiloc=1,nobxloc
      i = iobxloc(iiloc); j = jobxloc(iiloc)
      ii = indexobx(iiloc); ic = MERGE(i,i-1,westobx(ii))
      xfac = delxatv(ic,j)/delyatuv(i,j)
      DS_strain(i,j,1) = xfac*(umvel(i,j)/delxatc(ic,j)-&
                             & umvel(i,j-1)/delxatc(ic,j-1))
   ENDDO iiloc_521

!  ---Y-open boundaries
   jjloc_522: DO jjloc=1,nobyloc
      i = iobyloc(jjloc); j = jobyloc(jjloc)
      jj = indexoby(jjloc); jc = MERGE(j,j-1,soutoby(jj))
      xfac = delyatu(i,jc)/delxatuv(i,j)
      DS_strain(i,j,1) = xfac*(vmvel(i,j)/delyatc(i,jc)-&
                             & vmvel(i-1,j)/delyatc(i-1,jc))
   ENDDO jjloc_522

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = 1; nhdims = (/0,1,0,1/)
      CALL exchange_mod(DS_strain,lbounds,nhdims,0)
   ENDIF

!
!5.3 Diffusion coefficients at C-nodes
!-------------------------------------
!

   CALL UVarr_at_C(DS_strain(:,:,1),array2dc2,0,1,(/1,1,nz/),&
                 & (/ncloc+1,nrloc+1,nz/),1,0,.TRUE.)
   WHERE (maskatc_int)
      hdifcoef2datc(1:ncloc,1:nrloc) = garea*deptotatc(1:ncloc,1:nrloc)*&
                          & SQRT(DT_tension(1:ncloc,1:nrloc,1)**2+array2dc2**2)
   END WHERE

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds(1:2) = 0; nhdims = (/1,0,1,0/)
      CALL exchange_mod(hdifcoef2datc,lbounds(1:2),nhdims,iarr_hdifcoef2datc)
   ENDIF

!
!5.4 Diffusion coefficients at UV-nodes
!--------------------------------------
!
!  ---interior
   CALL Carr_at_UV(DT_tension(:,:,1),array2duv2,1,1,(/0,0,nz/),&
                & (/ncloc,nrloc,nz/),1,0,.TRUE.)
   WHERE (maskatuv(:,:,nz))
      hdifcoef2datuv(1:ncloc,1:nrloc) = delxatuv(1:ncloc,1:nrloc)*&
                        & delyatuv(1:ncloc,1:nrloc)*&
                        & deptotatuv(1:ncloc,1:nrloc)*&
                        & SQRT(array2duv2**2+DS_strain(1:ncloc,1:nrloc,1)**2)
   END WHERE

!  ---X-open boundaries
   iiloc_541: DO iiloc=1,nobxloc
      i = iobxloc(iiloc); j = jobxloc(iiloc)
      ii = indexobx(iiloc); ic = MERGE(i,i-1,westobx(ii))
      hdifcoef2datuv(i,j) = 0.5*(hdifcoef2datc(ic,j-1)+hdifcoef2datc(ic,j))
   ENDDO iiloc_541

!  ---Y-open boundaries
   jjloc_542: DO jjloc=1,nobyloc
      i = iobyloc(jjloc); j = jobyloc(jjloc)
      jj = indexoby(jjloc); jc = MERGE(j,j-1,soutoby(jj))
      hdifcoef2datuv(i,j) = 0.5*(hdifcoef2datc(i-1,jc)+hdifcoef2datc(i,jc))
   ENDDO jjloc_542

!
!6. 3-D grid
!-----------
!

ELSE

!
!6.1 Horizontal tension
!----------------------
!

   WHERE (maskatc)
      array2dc1 = delyatc(0:ncloc,0:nrloc)/delxatc(0:ncloc,0:nrloc)
   END WHERE
   WHERE (.NOT.maskatu(:,:,nz)) array2du1 = 0.0
   WHERE (.NOT.maskatv(:,:,nz)) array2dv1 = 0.0
   k_610: DO k=1,nz
      WHERE (maskatu(:,:,k)) 
         array2du1 = uvel(0:ncloc+1,0:nrloc,k)/delyatu(0:ncloc+1,0:nrloc)
      END WHERE
      WHERE (maskatv(:,:,k))
         array2dv1 = vvel(0:ncloc,0:nrloc+1,k)/delxatv(0:ncloc,0:nrloc+1)
      END WHERE
      WHERE (maskatc)
         DT_tension(:,:,k) = array2dc1*&
                      & (array2du1(1:ncloc+1,:)-array2du1(0:ncloc,:)) - &
                      & (array2dv1(:,1:nrloc+1)-array2dv1(:,0:nrloc))/array2dc1
      END WHERE
   ENDDO k_610

!
!6.2 Horizontal shearing strain
!-----------------------------
!

   WHERE (maskatuv(:,:,nz))
      array2duv1 = delxatuv(1:ncloc,1:nrloc)/delyatuv(1:ncloc,1:nrloc)
   END WHERE
   WHERE (.NOT.maskatu(1:ncloc,:,nz)) array2du2 = 0.0
   WHERE (.NOT.maskatv(:,1:nrloc,nz)) array2dv2 = 0.0
   k_621: DO k=1,nz
      WHERE (maskatu(1:ncloc,:,k))
         array2du2 = uvel(1:ncloc,0:nrloc,k)/delxatu(1:ncloc,0:nrloc)
      END WHERE
      WHERE (maskatv(:,1:nrloc,k))
         array2dv2 = vvel(0:ncloc,1:nrloc,k)/delyatv(0:ncloc,1:nrloc)
      END WHERE
      WHERE (maskatuv(:,:,k))
         DS_strain(1:ncloc,1:nrloc,k) = array2duv1*&
                     & (array2du2(:,1:nrloc)-array2du2(:,0:nrloc-1)) + &
                     & (array2dv2(1:ncloc,:)-array2dv2(0:ncloc-1,:))/array2duv1
      END WHERE
   ENDDO k_621

!  ---X-open boundaries
   iiloc_622: DO iiloc=1,nobxloc
      i = iobxloc(iiloc); j = jobxloc(iiloc)
      ii = indexobx(iiloc); ic = MERGE(i,i-1,westobx(ii))
      xfac = delxatv(ic,j)/delyatuv(i,j)
      DS_strain(i,j,:) = xfac*(uvel(i,j,:)/delxatc(ic,j)-&
                             & uvel(i,j-1,:)/delxatc(ic,j-1))
   ENDDO iiloc_622

!  ---Y-open boundaries
   jjloc_623: DO jjloc=1,nobyloc
      i = iobyloc(jjloc); j = jobyloc(jjloc)
      jj = indexoby(jjloc); jc = MERGE(j,j-1,soutoby(jj))
      xfac = delyatu(i,jc)/delxatuv(i,j)
      DS_strain(i,j,:) = xfac*(vvel(i,j,:)/delyatc(i,jc)-&
                             & vvel(i-1,j,:)/delyatc(i-1,jc))
   ENDDO jjloc_623

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = 1; nhdims = (/0,1,0,1/)
      CALL exchange_mod(DS_strain,lbounds,nhdims,0)
   ENDIF

!
!6.3 Diffusion coefficients at C-nodes
!-------------------------------------
!

   WHERE (maskatc_int) 
      hdifcoef2datc(1:ncloc,1:nrloc) = 0.0
   END WHERE
   k_630: DO k=1,nz
      CALL UVarr_at_C(DS_strain(:,:,k),array2dc2,0,1,(/1,1,k/),&
                   & (/ncloc+1,nrloc+1,k/),1,0,info)
      WHERE (maskatc_int)
         hdifcoef3datc(1:ncloc,1:nrloc,k) = garea*&
                          & SQRT(DT_tension(1:ncloc,1:nrloc,k)**2+array2dc2**2)
         hdifcoef2datc(1:ncloc,1:nrloc) = hdifcoef2datc(1:ncloc,1:nrloc) +&
                & hdifcoef3datc(1:ncloc,1:nrloc,k)*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_630

!   ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds  = (/0,0,1/); nhdims = (/1,0,1,0/)
      CALL exchange_mod(hdifcoef2datc,lbounds(1:2),nhdims,iarr_hdifcoef2datc)
      CALL exchange_mod(hdifcoef3datc,lbounds,nhdims,iarr_hdifcoef3datc)
   ENDIF

!
!6.4 Diffusion coefficients at U-nodes
!-------------------------------------
!

   WHERE (maskatu(1:ncloc,1:nrloc,nz))
      array2du1(1:ncloc,1:nrloc) = delxatu(1:ncloc,1:nrloc)*&
                                 & delyatu(1:ncloc,1:nrloc)
   END WHERE
   k_640: DO k=1,nz
      CALL Carr_at_U(DT_tension(:,1:nrloc,k),array2du2(:,1:nrloc),1,1,&
                  & (/0,1,k/),(/ncloc,nrloc,k/),1,0,info)
      CALL UVarr_at_U(DS_strain(1:ncloc,:,k),array2du3,0,3,(/1,1,k/),&
                   & (/ncloc,nrloc+1,k/),1,0,info)
      WHERE (maskatu(1:ncloc,1:nrloc,k))
         hdifcoef3datu(1:ncloc,:,k) = array2du1(1:ncloc,1:nrloc)*&
                                   & SQRT(array2du2(:,1:nrloc)**2+array2du3**2)
      END WHERE
   ENDDO k_640

!
!6.5 Diffusion coefficients at V-nodes
!-------------------------------------
!

   WHERE (maskatv(1:ncloc,1:nrloc,nz)) 
      array2dv1(1:ncloc,1:nrloc) = delxatv(1:ncloc,1:nrloc)*&
                                 & delyatv(1:ncloc,1:nrloc)
   END WHERE
   k_650: DO k=1,nz
      CALL Carr_at_V(DT_tension(1:ncloc,:,k),array2dv2(1:ncloc,:),1,1,&
                  & (/1,0,k/),(/ncloc,nrloc,k/),1,0,info)
      CALL UVarr_at_V(DS_strain(:,1:nrloc,k),array2dv3,0,3,(/1,1,k/),&
                   & (/ncloc+1,nrloc,k/),1,0,info)
      WHERE (maskatv(1:ncloc,1:nrloc,k))
         hdifcoef3datv(:,1:nrloc,k) = array2dv1(1:ncloc,1:nrloc)*&
                                  & SQRT(array2dv2(1:ncloc,:)**2+array2dv3**2)
      END WHERE
   ENDDO k_650

!
!6.6 Diffusion coefficients at UV-nodes
!--------------------------------------
!
!  ---interior
   WHERE (maskatuv(:,:,nz)) 
      array2duv1 = delxatuv(1:ncloc,1:nrloc)*delyatuv(1:ncloc,1:nrloc)
      hdifcoef2datuv(1:ncloc,1:nrloc) = 0.0
   END WHERE
   k_661: DO k=1,nz
      CALL Carr_at_UV(DT_tension(:,:,k),array2duv2,1,1,(/0,0,k/),&
                   & (/ncloc,nrloc,k/),1,0,info)
      WHERE (maskatuv(:,:,k))
         hdifcoef3datuv(1:ncloc,1:nrloc,k) = array2duv1*&
                          & SQRT(array2duv2**2+DS_strain(1:ncloc,1:nrloc,k)**2)
         hdifcoef2datuv(1:ncloc,1:nrloc) = hdifcoef2datuv(1:ncloc,1:nrloc) + &
              & hdifcoef3datuv(1:ncloc,1:nrloc,k)*delzatuv(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_661

!   ---X-open boundaries
   iiloc_662: DO iiloc=1,nobxloc
      i = iobxloc(iiloc); j = jobxloc(iiloc)
      ii = indexobx(iiloc); ic = MERGE(i,i-1,westobx(ii))
      hdifcoef2datuv(i,j) = 0.5*(hdifcoef2datc(ic,j-1)+hdifcoef2datc(ic,j))
      hdifcoef3datuv(i,j,:) = 0.5*(hdifcoef3datc(ic,j-1,:)+&
                                 & hdifcoef3datc(ic,j,:))
   ENDDO iiloc_662

!   ---Y-open boundaries
   jjloc_663: DO jjloc=1,nobyloc
      i = iobyloc(jjloc); j = jobyloc(jjloc)
      jj = indexoby(jjloc); jc = MERGE(j,j-1,soutoby(jj))
      hdifcoef2datuv(i,j) = 0.5*(hdifcoef2datc(i-1,jc)+hdifcoef2datc(i,jc))
      hdifcoef3datuv(i,j,:) = 0.5*(hdifcoef3datc(i-1,jc,:)+&
                                 & hdifcoef3datc(i,jc,:))
   ENDDO jjloc_663

!
!6.7 Diffusion coefficients at W-nodes
!-------------------------------------
!

   IF (iopt_hdif_scal.GE.2) THEN
      CALL Carr_at_W(hdifcoef3datc(1:ncloc,1:nrloc,:),hdifcoef3datw(:,:,2:nz),&
           & (/1,1,1/),(/ncloc,nrloc,nz/),1,iarr_hdifcoef3datw,.TRUE.)
   ENDIF
   
ENDIF
   
!
!7. Exchange halo sections
!---------------------------
!

IF (iopt_MPI.EQ.1) THEN

!
!7.1 3-D arrays
!--------------
!

   IF (iopt_grid_nodim.EQ.3) THEN
      lbounds = 1
      nhdims = (/0,1,0,0/)
      CALL exchange_mod(hdifcoef3datu,lbounds,nhdims,iarr_hdifcoef3datu)
      nhdims = (/0,0,0,1/)
      CALL exchange_mod(hdifcoef3datv,lbounds,nhdims,iarr_hdifcoef3datv)
      nhdims = (/0,1,0,1/)
      CALL exchange_mod(hdifcoef3datuv,lbounds,nhdims,iarr_hdifcoef3datuv)
   ENDIF

!
!7.2 2-D arrays
!--------------
!

   lbounds(1:2) = 1; nhdims = (/0,1,0,1/)
   CALL exchange_mod(hdifcoef2datuv,lbounds(1:2),nhdims,iarr_hdifcoef2datuv)

ENDIF

!
!8. Isoneutral slope factors in case of rotated diffusion
!--------------------------------------------------------
!

IF (iopt_hdif_scal.EQ.3) CALL isoslopes

!
!8. Deallocate arrays
!--------------------
!

IF (iopt_hdif_coef.EQ.2) THEN
   DEALLOCATE (maskatc,maskatu,maskatv,maskatuv)
   DEALLOCATE (DT_tension,DS_strain)
   DEALLOCATE (array2dc1,array2dc2,array2du1,array2du2,array2dv1,array2dv2, &
             & array2duv1,array2duv2)
   IF (iopt_grid_nodim.EQ.3) DEALLOCATE (array2du3,array2dv3)
ENDIF

1000 CALL log_timer_out(npcc,itm_hdif)
 

RETURN

END SUBROUTINE horizontal_diff_coefs

!========================================================================

SUBROUTINE isoslopes
!************************************************************************
!
! *isoslopes* Isoneutral slope factors 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Coefficients.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
! Reference -
!
! Calling program -
!
! External calls -
!
! Module calls - error_alloc, Zcoord_arr
!
!************************************************************************
!
USE density
USE diffusion  
USE grid
USE gridpars
USE iopars
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: Zcoord_arr
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: salflag, tempflag
INTEGER :: i, ip, ipp, j, k, kk, kkp, kp
REAL, PARAMETER :: slopemin_iso = 1.0E-25
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc, maskatu, maskatv
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: kmix
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, array2du1, array2du2, &
                                         & array2dv1, array2dv2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: deltau, deltav, zcoordatw
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: xdensslope, xslopeatu_mld, &
                                             & ydensslope, yslopeatv_mld, &
                                             & zdensslope


procname(pglev+1) = 'isoslopes'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatc(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc+2,nrloc+2/),kndlog)
ALLOCATE (maskatu(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc+2,nrloc/),kndlog)
ALLOCATE (maskatv(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc+2/),kndlog)
ALLOCATE (array2dc(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc+2,nrloc+2/),kndrtype)
ALLOCATE (array2du1(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2du1',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (array2du2(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2du2',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (array2dv1(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2dv1',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (array2dv2(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2dv2',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (deltau(ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('deltau',3,(/ncloc+1,nrloc,nz/),kndrtype)
ALLOCATE (deltav(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('deltav',3,(/ncloc,nrloc+1,nz/),kndrtype)
ALLOCATE (kmix(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('kmix',2,(/ncloc+2,nrloc+2/),kndint)
ALLOCATE (xslopeatu_mld(0:ncloc+1,nrloc,0:1,0:1),STAT=errstat)
CALL error_alloc('xslopeatu_mld',4,(/ncloc+2,nrloc,1,1/),kndrtype)
ALLOCATE (xdensslope(0:ncloc+1,nrloc,nz,0:1),STAT=errstat)
CALL error_alloc('xdensslope',3,(/ncloc+2,nrloc,nz/),kndrtype)
ALLOCATE (ydensslope(ncloc,0:nrloc+1,nz,0:1),STAT=errstat)
CALL error_alloc('ydensslope',3,(/ncloc,nrloc+2,nz/),kndrtype)
ALLOCATE (yslopeatv_mld(ncloc,0:nrloc+1,0:1,0:1),STAT=errstat)
CALL error_alloc('yslopeatvu_mld',4,(/ncloc,nrloc+2,1,1/),kndrtype)
ALLOCATE (zcoordatw(0:ncloc+1,0:nrloc+1,nz+1),STAT=errstat)
CALL error_alloc('zcoordatw',3,(/ncloc+2,nrloc,nz+1/),kndrtype)
ALLOCATE (zdensslope(0:ncloc+1,0:nrloc+1,nz,0:1),STAT=errstat)
CALL error_alloc('zdensslope',3,(/ncloc+2,nrloc+2,nz/),kndrtype)

!
!2. Initialise
!-------------
!   
!---mask arrays
maskatc = nodeatc(0:ncloc+1,0:nrloc+1).EQ.1
maskatu = node2du(1:ncloc+1,1:nrloc).EQ.2
maskatv = node2dv(1:ncloc,1:nrloc+1).EQ.2   

!---Z-coordinate at W-nodes
CALL Zcoord_arr(zcoordatw,(/0,0,1/),(/ncloc+1,nrloc+1,nz+1/),'W  ',.FALSE.)

!---density slopes
xdensslope = 0.0; ydensslope = 0.0; zdensslope = 0.0

!---("triad") slope factors
xslopeatu_siso = 0.0; yslopeatv_siso = 0.0
xslopeatu_ziso = 0.0; yslopeatv_ziso = 0.0
xslopeatu_mld = 0.0; yslopeatv_mld = 0.0

!---vertical diffusion coefficient
vdifcoefscal_rot = 0.0

!---tracer flags
tempflag = iopt_temp.EQ.2.OR.(iopt_temp.EQ.1.AND.nt.EQ.0)
salflag = iopt_sal.EQ.2.OR.(iopt_sal.EQ.1.AND.nt.EQ.0)

!
!3. Delta factor for limiting conditions
!---------------------------------------
!

IF (iopt_hdif_lim.GT.1) THEN

   k_310: DO k=1,nz
   
!     ---X-direction
      WHERE (maskatu)
         deltau(:,:,k) = 0.5*delxatu(1:ncloc+1,1:nrloc)*delzatu(1:ncloc+1,:,k)&
                       & /(hdifscal_fac*hdifcoef3datu(:,:,k)*delt3d)
      END WHERE

!     ---Y-direction
      WHERE (maskatv)
         deltav(:,:,k) = 0.5*delyatv(1:ncloc,1:nrloc+1)*delzatv(:,1:nrloc+1,k)&
                       & /(hdifscal_fac*hdifcoef3datv(:,:,k)*delt3d)
      END WHERE

   ENDDO k_310
      
ENDIF

!
!4. Density slopes
!-----------------
!
!4.1 X-direction
!---------------
!

ip_410: DO ip=0,1
k_410: DO k=1,nz
      
   array2du1 = 0.0; array2du2 = 0.0
      
!  ---temperature
   IF (tempflag) THEN
      WHERE (maskatu)
         array2du1 = temp(1:ncloc+1,1:nrloc,k) - temp(0:ncloc,1:nrloc,k)
      END WHERE
      WHERE (maskatc(ip:ip+ncloc,1:nrloc)) 
         array2du2 = -beta_temp(ip:ip+ncloc,1:nrloc,k)*array2du1
      END WHERE
   ENDIF
      
!  ---salinity
   IF (salflag) THEN
      WHERE (maskatu)
         array2du1 = sal(1:ncloc+1,1:nrloc,k) - sal(0:ncloc,1:nrloc,k)
      END WHERE
      WHERE (maskatc(ip:ip+ncloc,1:nrloc))
         array2du2 = array2du2 + beta_sal(ip:ip+ncloc,1:nrloc,k)*array2du1
      END WHERE
   ENDIF

!  ---density slope
   WHERE (maskatu)
      array2du2 = array2du2/delxatu(1:ncloc+1,1:nrloc)
      xdensslope(ip:ip+ncloc,:,k,1-ip) = SIGN(MAX(slopemin_iso,&
                                            & ABS(array2du2)),array2du2)
   END WHERE
   
ENDDO k_410
ENDDO ip_410

!
!4.2 Y-direction
!---------------
!

ip_420: DO ip=0,1
k_420: DO k=1,nz

   array2dv1 = 0.0; array2dv2 = 0.0      
!  ---temperature
   IF (tempflag) THEN
      WHERE (maskatv)
         array2dv1 = temp(1:ncloc,1:nrloc+1,k) - temp(1:ncloc,0:nrloc,k)
      END WHERE
      WHERE (maskatc(1:ncloc,ip:ip+nrloc))
         array2dv2 = -beta_temp(1:ncloc,ip:ip+nrloc,k)*array2dv1
      END WHERE
   ENDIF
      
!  ---salinity
   IF (salflag) THEN
      WHERE (maskatv)
         array2dv1 = sal(1:ncloc,1:nrloc+1,k) - sal(1:ncloc,0:nrloc,k)
      END WHERE

      WHERE (maskatc(1:ncloc,ip:ip+nrloc))
         array2dv2 = array2dv2 + beta_sal(1:ncloc,ip:ip+nrloc,k)*array2dv1
      END WHERE
   ENDIF

!  ---density slope
   WHERE (maskatv)
      array2dv2 = array2dv2/delyatv(1:ncloc,1:nrloc+1)
      ydensslope(:,ip:ip+nrloc,k,1-ip) = SIGN(MAX(slopemin_iso,&
                                            & ABS(array2dv2)),array2dv2)
   END WHERE
   
ENDDO k_420
ENDDO ip_420

!
!4.3 Z-direction
!--------------
!

kp_430: DO kp=0,1
k_430: DO k=1,nz
      
   array2dc = 0.0
   kkp = k+1-kp
   IF (kkp.GE.2.AND.kkp.LE.nz) THEN

!     ---temperature
      IF (tempflag) THEN
         WHERE (maskatc)
            array2dc = -beta_temp(0:ncloc+1,0:nrloc+1,k)*&
                        & (temp(0:ncloc+1,0:nrloc+1,kkp)-&
                        &  temp(0:ncloc+1,0:nrloc+1,kkp-1))
         END WHERE
      ENDIF

!     ---salinity
      IF (salflag) THEN
         WHERE (maskatc)
            array2dc = array2dc + beta_sal(0:ncloc+1,0:nrloc+1,k)*&
                        & (sal(0:ncloc+1,0:nrloc+1,kkp)-&
                        &  sal(0:ncloc+1,0:nrloc+1,kkp-1))
         END WHERE
      ENDIF

!     ---density slope
      WHERE (maskatc)
         array2dc = array2dc/delzatw(0:ncloc+1,0:nrloc+1,kkp)
         zdensslope(:,:,k,kp) = MIN(-slopemin_iso,array2dc)
      END WHERE
         
   ENDIF

ENDDO k_430
ENDDO kp_430

WHERE (maskatc)
   zdensslope(:,:,nz,0) = -slopemin_iso
END WHERE
WHERE (maskatc)
   zdensslope(:,:,1,1) = -slopemin_iso
END WHERE

!
!5. Density slope wrt z-coordinate surfaces
!------------------------------------------
!
!5.1 Without limits
!------------------
!

kp_510: DO kp=0,1
ip_510: DO ip=0,1
   array2du1 = 0.0; array2dv1 = 0.0

   k_511: DO k=1,nz
      
!     ---X-direction   
      WHERE (maskatu)
         xslopeatu_siso(ip:ip+ncloc,:,k,1-ip,kp) = &
                                      & xdensslope(ip:ip+ncloc,:,k,1-ip)&
                                     & /zdensslope(ip:ip+ncloc,1:nrloc,k,kp)
         array2du1 = xslopeatu_siso(ip:ip+ncloc,:,k,1-ip,kp) - &
                   & xslopeatu_geo(:,:,k)
         xslopeatu_ziso(ip:ip+ncloc,:,k,1-ip,kp) = SIGN(MIN(slopemax_iso,&
                                                 & ABS(array2du1)),array2du1)
      END WHERE
!     ---Y-direction   
      WHERE (maskatv)
         yslopeatv_siso(:,ip:ip+nrloc,k,1-ip,kp) = &
                                      & ydensslope(:,ip:ip+nrloc,k,1-ip)&
                                     & /zdensslope(1:ncloc,ip:ip+nrloc,k,kp)
         array2dv1 = yslopeatv_siso(:,ip:ip+nrloc,k,1-ip,kp) - &
                   & yslopeatv_geo(:,:,k)
         yslopeatv_ziso(:,ip:ip+nrloc,k,1-ip,kp) = SIGN(MIN(slopemax_iso,&
                                                 & ABS(array2dv1)),array2dv1)
     END WHERE
     
  ENDDO k_511
  
ENDDO ip_510
ENDDO kp_510

!
!5.2 Linear interpolation in surface mixed layer
!----------------------------------------------
!
!5.2.1 Vertical location of mixed layer
!--------------------------------------
!

kmix = 1
j_521: DO j=1,nrloc+1
i_521: DO i=1,ncloc+1
   IF (maskatc(i,j)) THEN
      k_5211: DO k=1,nz-1
         IF ((dens(i,j,k)-dens(i,j,nz)).GT.rho_crit_iso) kmix(i,j) = k
      ENDDO k_5211
   ENDIF
ENDDO i_521
ENDDO j_521

!
!5.2.2 Slope at mixed layer depth
!--------------------------------
!

kp_522: DO kp=0,1
ip_522: DO ip=0,1

!   ---X-direction      
   j_5221: DO j=1,nrloc
   i_5221: DO i=1,ncloc+1
      kk = kmix(i,j)
      IF (maskatu(i,j)) THEN
         xslopeatu_mld(i+ip-1,j,1-ip,kp) = &
       & xslopeatu_ziso(i+ip-1,j,kk+kp,1-ip,kp)
      ENDIF
   ENDDO i_5221
   ENDDO j_5221

!  ---Y-direction 
   j_5222: DO j=1,nrloc+1
   i_5222: DO i=1,ncloc   
      IF (maskatv(i,j)) THEN
         yslopeatv_mld(i,j+ip-1,1-ip,kp) = &
       & yslopeatv_ziso(i,j+ip-1,kk+kp,1-ip,kp)
      ENDIF
   ENDDO i_5222
   ENDDO j_5222

ENDDO ip_522
ENDDO kp_522

!
!5.2.3 Interpolate
!-----------------
!

kp_523: DO kp=0,1
ip_523: DO ip=0,1

!  ---X-direction
   j_5231: DO j=1,nrloc
   i_5231: DO i=1,ncloc+1
      kk = kmix(i,j)
      IF (maskatu(i,j)) THEN
         kk_52311: DO k=kk,nz
            xslopeatu_ziso(i+ip-1,j,k,1-ip,kp) = &
                                & xslopeatu_mld(i+ip-1,j,1-ip,kp)*&
                                & (zcoordatw(i+ip-1,j,k)/zcoordatw(i+ip-1,j,kk))
         ENDDO kk_52311
      ENDIF
   ENDDO i_5231
   ENDDO j_5231
   
!  ---Y-direction
   j_5232: DO j=1,nrloc+1
   i_5232: DO i=1,ncloc
      kk = kmix(i,j)
      IF (maskatv(i,j)) THEN
         kk_52321: DO k=kk,nz
            yslopeatv_ziso(i,j+ip-1,k,1-ip,kp) = &
                                & yslopeatv_mld(i,j+ip-1,1-ip,kp)*&
                                & (zcoordatw(i,j+ip-1,k)/zcoordatw(i,j+ip-1,kk))
         ENDDO kk_52321
      ENDIF
   ENDDO i_5232
   ENDDO j_5232
   
ENDDO ip_523
ENDDO kp_523

!
!6. Density slope wrt s-coordinate surfaces
!------------------------------------------
!

kp_610: DO kp=0,1
ip_610: DO ip=0,1
k_610: DO k=1,nz
   
!
!6.1 Without limiting condition
!-------------------------------
!

!  ---X-direction
   WHERE (maskatu)
      xslopeatu_siso(ip:ip+ncloc,:,k,1-ip,kp) = &
                & xslopeatu_ziso(ip:ip+ncloc,:,k,1-ip,kp) + xslopeatu_geo(:,:,k)
   END WHERE
!  ---Y-direction
   WHERE (maskatv)
      yslopeatv_siso(:,ip:ip+nrloc,k,1-ip,kp) = &
                & yslopeatv_ziso(:,ip:ip+nrloc,k,1-ip,kp) + yslopeatv_geo(:,:,k)
   END WHERE

!
!6.2 NEMO limiting condition
!---------------------------
!

   IF (iopt_hdif_lim.EQ.1) THEN

!     ---X-direction
      WHERE (maskatu)
         array2du1 = xslopeatu_siso(ip:ip+ncloc,:,k,1-ip,kp)
         array2du2 = array2du1**2/xslopeatu_siso(ip:ip+ncloc,:,k,1-ip,kp)
         xslopeatu_siso(ip:ip+ncloc,:,k,1-ip,kp) = &
                      & SIGN(MIN(ABS(array2du1),ABS(array2du2)),array2du1) 
      END WHERE
!     ---Y-direction
      WHERE (maskatv)
         array2dv1 = yslopeatv_siso(:,ip:ip+nrloc,k,1-ip,kp)
         array2dv2 = array2dv1**2/yslopeatv_siso(:,ip:ip+nrloc,k,1-ip,kp)
         yslopeatv_siso(:,ip:ip+nrloc,k,1-ip,kp) = &
                     & SIGN(MIN(ABS(array2dv1),ABS(array2dv2)),array2dv1) 
      END WHERE
      
!
!6.3 Gerdes et al. limiting condition
!------------------------------------
!
      
   ELSEIF (iopt_hdif_lim.EQ.2) THEN

!     ---X-direction
      WHERE (maskatu)
         array2du1 = xslopeatu_siso(ip:ip+ncloc,:,k,1-ip,kp)
      END WHERE
      WHERE (maskatu.AND.deltau(:,:,k).LT.ABS(array2du1))
         xslopeatu_siso(ip:ip+ncloc,:,k,1-ip,kp) = deltau(:,:,k)**2/array2du1
      END WHERE

!     ---Y-direction
      WHERE (maskatv)
         array2dv1 = yslopeatv_siso(:,ip:ip+nrloc,k,1-ip,kp)
      END WHERE
      WHERE (maskatv.AND.deltav(:,:,k).LT.ABS(array2dv1))
         yslopeatv_siso(:,ip:ip+nrloc,k,1-ip,kp) = deltav(:,:,k)**2/array2dv1
      END WHERE

   ENDIF
   
ENDDO k_610
ENDDO ip_610
ENDDO kp_610

!
!7. Isoneutral vertical diffusion coefficient
!--------------------------------------------
!
   
kp_710: DO kp=0,1
ip_710: DO ip=0,1

   ipp = 1-ip
   WHERE (maskatu(ipp+1:ipp+ncloc,:))
      array2du1(1:ncloc,:) = delxatu(ipp+1:ipp+ncloc,1:nrloc)*&
                           & delyatu(ipp+1:ipp+ncloc,1:nrloc)/garea
   END WHERE
   WHERE (maskatv(:,ipp+1:ipp+nrloc))
      array2dv1(:,1:nrloc) = delxatv(1:ncloc,ipp+1:ipp+nrloc)*&
                           & delyatv(1:ncloc,ipp+1:ipp+nrloc)/garea
   END WHERE
      
   k_711: DO k=2,nz
      kkp = k+kp-1

!     ---X-direction
      WHERE (maskatu(ipp+1:ipp+ncloc,:))
         array2du2(1:ncloc,:) = 0.25*(hdifscal_fac*&
                              & hdifcoef3datu(ipp+1:ipp+ncloc,:,k)*&
                              & array2du1(1:ncloc,:)*&
                              & xslopeatu_ziso(1:ncloc,:,kkp,1-ip,kp)**2)
      END WHERE
      WHERE (maskatc_int)
         vdifcoefscal_rot(:,:,k) = vdifcoefscal_rot(:,:,k) + &
                                 & array2du2(1:ncloc,:)
      END WHERE
!     ---Y-direction
      WHERE (maskatv(:,ipp+1:ipp+nrloc))
         array2dv2(:,1:nrloc) = 0.25*(hdifscal_fac*&
                              & hdifcoef3datv(:,ipp+1:ipp+nrloc,k)*&
                              & array2dv1(:,1:nrloc)*&
                              & yslopeatv_ziso(:,1:nrloc,kkp,1-ip,kp)**2)
      END WHERE
      WHERE (maskatc_int)
         vdifcoefscal_rot(:,:,k) = vdifcoefscal_rot(:,:,k) + &
                                 & array2dv2(:,1:nrloc)
      END WHERE
      
   ENDDO k_711

ENDDO ip_710
ENDDO kp_710

!
!7. Deallocate
!-------------
!

DEALLOCATE (maskatc,maskatu,maskatv,array2dc,array2du1,array2du2, array2dv1,&
          & array2dv2,deltau,deltav,kmix,xdensslope,xslopeatu_mld, ydensslope,&
          & yslopeatv_mld,zcoordatw,zdensslope)

CALL log_timer_out()


RETURN

END SUBROUTINE isoslopes

!========================================================================

SUBROUTINE kinematic_viscosity
!************************************************************************
!
! *kinematic viscosity* kinematic viscosity of sea water
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Diffusion_Coefficients.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - vertical_diff_coefs
!
! External calls - 
!
! Module calls - error_alloc, Carr_at_W
!
!************************************************************************
!
USE density
USE diffusion
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE array_interp, ONLY: Carr_at_W
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: k
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: tempatw


procname(pglev+1) = 'kinematic_viscosity'
CALL log_timer_in()

!
!1. Allocate
!-----------
!

ALLOCATE (array2dc(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc+1,nrloc+1/),kndrtype)
ALLOCATE (tempatw(0:ncloc,0:nrloc,nz+1),STAT=errstat)
CALL error_alloc('tempatw',3,(/ncloc+1,nrloc+1,nz+1/),kndrtype)

!
!2. Constant value
!-----------------
!

IF (iopt_kinvisc.EQ.0) THEN
   kinvisc = kinvisc_cst

!
!3. Temperature dependent value
!------------------------------
!

ELSEIF (iopt_kinvisc.EQ.1) THEN
   CALL Carr_at_W(temp(0:ncloc,0:nrloc,:),tempatw(:,:,1:nz),(/0,0,0/),&
               & (/ncloc,nrloc,nz/),1,iarr_temp,.TRUE.)
   WHERE (nodeatc(0:ncloc,0:nrloc).EQ.1)
      tempatw(:,:,1) = tempatw(:,:,2)
      tempatw(:,:,nz+1) = tempatw(:,:,nz)
   END WHERE
   k_310: DO k=1,nz+1
      WHERE (nodeatc(0:ncloc,0:nrloc).EQ.1)
         array2dc = tempatw(:,:,k)-1.0
         kinvisc(:,:,k) = 1.0E-06*(1.7688+(0.659E-03*array2dc-0.05076)*array2dc)
      END WHERE
   ENDDO k_310
ENDIF

!
!4. Deallocate
!-------------
!

DEALLOCATE (array2dc,tempatw)

CALL log_timer_out()


RETURN

END SUBROUTINE kinematic_viscosity

!========================================================================

SUBROUTINE vertical_diff_coefs
!************************************************************************
!
! *vertical_diff_coefs* Vertical diffusion coefficients
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Coefficients.f90  V2.5
!
! Description - 
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! External calls - buoyancy_frequency, dissipation_equation, dissip_lim_conds,
!                  dissip_to_zlmix, eddy_coefs_alg, eddy_coefs_tc,
!                  init_turbulence, kinematic_viscosity, kl_equation,
!                  mixing_length, shear_frequency, tke_equation,
!                  tke_equilibrium, zlmix_lim_conds, zlmix_to_dissip
!
! Module calls -  
!                 
!************************************************************************
!
USE iopars
USE physpars
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc

procname(pglev+1) = 'vertical_diff_coefs'
CALL log_timer_in(npcc)

!
!1. Initialise on first call
!---------------------------
!

IF (nt.EQ.0) THEN

!  ---initialise turbulence parameters
   IF (iopt_vdif_coef.EQ.3) CALL init_turbulence

!  ---kinematic viscoity
   CALL kinematic_viscosity

!  ---shear and buoyancy frequencies
   CALL shear_frequency
   IF (iopt_dens.GT.0.OR.iopt_sed.GT.0) CALL buoyancy_frequency

!  ---algebraic case 
   IF (iopt_vdif_coef.EQ.2) CALL eddy_coefs_alg

!  ---mixing length, dissipation rate
   IF (iopt_vdif_coef.EQ.3) THEN
      IF (iopt_turb_ntrans.LT.2) CALL mixing_length
      IF ((iopt_turb_ntrans.LT.2).OR.(iopt_turb_param.EQ.1)) THEN
         CALL zlmix_to_dissip
      ENDIF
   ENDIF

!  ---diffusion coefficients
   IF (iopt_vdif_coef.EQ.3) CALL eddy_coefs_tc

   GOTO 1000

ENDIF

!
!2. Forcing frequencies
!----------------------
!
!---shear
CALL shear_frequency
!--stratification
IF (iopt_dens.GT.0.OR.iopt_sed.GT.0) CALL buoyancy_frequency

!
!3. Kinematic viscosity
!----------------------
!

IF (iopt_kinvisc.EQ.1.AND.iopt_temp.EQ.2) CALL kinematic_viscosity

!
!4. Algebraic case
!-----------------
!
 
IF (iopt_vdif_coef.EQ.2) THEN
   CALL eddy_coefs_alg
   GOTO 1000
ENDIF

!
!5. Turbulence closure schemes
!-----------------------------
!
!5.1 Equilibrium models
!----------------------
!

IF (iopt_turb_ntrans.EQ.0) THEN

!  ---turbulence energy
   CALL tke_equilibrium

!  ---mixing length
   CALL mixing_length

!  ---dissipation rate
   CALL zlmix_to_dissip

!
!5.2 Non-equilibrium models
!--------------------------
!

ELSE

!  ---turbulence energy
   IF (.NOT.physinit) CALL tke_equation

!
!5.2.1 k-l case
!--------------
!

   IF (iopt_turb_param.EQ.1) THEN

!     ---mixing length
      IF (iopt_turb_ntrans.EQ.1) THEN
         CALL mixing_length
      ELSEIF (.NOT.physinit.AND.iopt_turb_ntrans.EQ.2) THEN
         CALL kl_equation
      ENDIF

!     ---limiting conditions
      IF (iopt_turb_iwlim.EQ.1) CALL zlmix_lim_conds

!     ---dissipation rate
      CALL zlmix_to_dissip

!
!5.2.2 k-epsilon case
!--------------------
!

   ELSEIF (iopt_turb_param.EQ.2) THEN

!     ---dissipation rate/mixing length
      IF (iopt_turb_ntrans.EQ.1) THEN
         CALL mixing_length
!        --limiting conditions
         IF (iopt_turb_iwlim.EQ.1) CALL zlmix_lim_conds
         CALL zlmix_to_dissip
      ELSEIF (.NOT.physinit.AND.iopt_turb_ntrans.EQ.2) THEN
         CALL dissipation_equation
!        --limiting conditions
         IF (iopt_turb_iwlim.EQ.1) CALL dissip_lim_conds
         CALL dissip_to_zlmix
      ENDIF

   ENDIF

ENDIF

!
!6. Eddy coefficients
!--------------------
!

CALL eddy_coefs_tc

1000 CONTINUE

CALL log_timer_out(npcc,itm_vdif)


RETURN

END SUBROUTINE vertical_diff_coefs

!========================================================================

SUBROUTINE vertical_diff_rot
!************************************************************************
!
! *vertical_diff_rot* Add contributions from rotated diffusion to
!                     scalar diffusion coefficient
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Coefficients.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! External calls - 
!
! Module calls -  
!                 
!************************************************************************
!
USE diffusion
USE grid
USE gridpars
USE iopars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: k, npcc


procname(pglev+1) = 'vertical_diff_rot'
CALL log_timer_in(npcc)

IF (iopt_vdif_coef.EQ.1) THEN
   k_110: DO k=2,nz
      WHERE (maskatc_int)
         vdifcoefscal(:,:,k) = vdifcoefscal_rot(:,:,k)
      END WHERE
   ENDDO k_110
ELSE
   k_120: DO k=2,nz
      WHERE (maskatc_int)
         vdifcoefscal(:,:,k) = vdifcoefscal(:,:,k) + vdifcoefscal_rot(:,:,k)
      END WHERE
   ENDDO k_120
ENDIF

CALL log_timer_out(npcc,itm_vdif)


RETURN

END SUBROUTINE vertical_diff_rot
