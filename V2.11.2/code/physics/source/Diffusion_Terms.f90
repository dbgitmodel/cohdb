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
! *Diffusion_Terms* Diffusion terms in transport equations
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
! Routines - Xdif_at_C, Ydif_at_C, Zdif_at_C, Xdif_at_U_3d, Xdif_at_U_2d,
!            Ydif_at_U_3d, Ydif_at_U_2d, Zdif_at_U, Xdif_at_V_3d,
!            Xdif_at_V_2d, Ydif_at_V_3d, Ydif_at_V_2d, Zdif_at_V,
!            Xdif_at_W, Ydif_at_W, Zdif_at_W  
!
!************************************************************************
!

!========================================================================

SUBROUTINE Xdif_at_C(psic,tridcfc,novars,klo,kup)
!************************************************************************
!
! *Xdif_at_C* Diffusion term in the Y-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2
!
! External calls - Xhdif_at_C_lap, Xhdif_at_C_geo, Xhdif_at_C_iso
! 
! Module calls -
!
!************************************************************************
!
USE gridpars
USE switches

IMPLICIT NONE
  
!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                     [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*     INTEGER Number of variables
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!

SELECT CASE (iopt_hdif_scal)
   CASE (1); CALL Xhdif_at_C_lap(psic,tridcfc,novars,klo,kup)
   CASE (2); CALL Xhdif_at_C_geo(psic,tridcfc,novars,klo,kup)
   CASE (3); CALL Xhdif_at_C_iso(psic,tridcfc,novars,klo,kup)
END SELECT


RETURN

END SUBROUTINE Xdif_at_C
  
!========================================================================

SUBROUTINE Xhdif_at_C_lap(psic,tridcfc,novars,klo,kup)
!************************************************************************
!
! *Xhdif_at_C_lap* Laplacian diffusion term in the X-direction for a quantity
!                  at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - laplacian diffusion
!
! Reference -
!
! Calling program - Xdif_at_C
!
! Module calls - error_alloc
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
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                     [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*     INTEGER Number of variables
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar, k, npcc
REAL :: fac
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2du
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2du
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3du, denom, difflux

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*difflux*    REAL    Diffusive flux (times factor) at U-nodes       
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xhdif_at_C_lap'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2du(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (array3du(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('array3du',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (difflux(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('difflux',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (maskatu(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc+1,nrloc,kup-klo+1/),kndlog)
ALLOCATE (mask2du(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('mask2du',2,(/ncloc+1,nrloc/),kndlog)

!
!2. Initialise
!-------------
!
!2.1 Mask array
!--------------
!

maskatu = nodeatu(1:ncloc+1,1:nrloc,klo:kup).EQ.2
mask2du = node2du(1:ncloc+1,1:nrloc).EQ.2

!
!2.2 Fluxes
!----------
!

fac = delt3d*hdifscal_fac
difflux = 0.0

!
!2.3 Work space arrays
!---------------------
!

IF (dYregX) THEN
   k_231: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delxatc(1:ncloc,1:nrloc)*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_231
ELSE
   k_232: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_232
ENDIF

IF (dYregX) THEN
   WHERE (mask2du)
      array2du = fac/delxatu(1:ncloc+1,1:nrloc)
   END WHERE
ELSE
   WHERE (mask2du)
      array2du = fac*delyatu(1:ncloc+1,1:nrloc)&
               & /delxatu(1:ncloc+1,1:nrloc)
   END WHERE
ENDIF

k_233: DO k=klo,kup
   WHERE (maskatu(:,:,k))
      array3du(:,:,k) = array2du*delzatu(1:ncloc+1,:,k)*hdifcoef3datu(:,:,k)
   END WHERE
ENDDO k_233

!
!3. Explicit terms
!-----------------
!

ivar_300: DO ivar=1,novars

!
!3.1 Fluxes
!----------
! 

   WHERE (maskatu)
      difflux = array3du*(psic(1:ncloc+1,1:nrloc,klo:kup,ivar)-&
                        & psic(0:ncloc,1:nrloc,klo:kup,ivar))
   END WHERE

!
!3.2 Add flux divergence to explicit term
!----------------------------------------
!

   k_320: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) + &
                              & (difflux(2:ncloc+1,:,k)-difflux(1:ncloc,:,k))&
                              & /denom(:,:,k)
      END WHERE
   ENDDO k_320

ENDDO ivar_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array2du,array3du,difflux,denom,maskatu,mask2du)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Xhdif_at_C_lap
  
!========================================================================

SUBROUTINE Xhdif_at_C_geo(psic,tridcfc,novars,klo,kup)
!************************************************************************
!
! *Xhdif_at_C_geo* Isolevel contribution to the diffusion term in the
!                  X-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - Xdif_at_C
!
! Module calls - Carr_at_U, error_alloc
!
!************************************************************************
!
USE diffusion
USE grid
USE gridpars
USE iopars
USE physpars
USE syspars
USE timepars
USE array_interp, ONLY: Carr_at_U
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                     [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*     INTEGER Number of variables
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar, k, npcc
REAL :: fac
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2du
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2du
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3du, denom, difflux, &
                                           & psicatu, zdiffatw

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*difflux*    REAL    Diffusive flux (times factor) at U-nodes       
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xhdif_at_C_geo'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2du(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (array3du(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('array3du',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (difflux(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('difflux',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (maskatu(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc+1,nrloc,kup-klo+1/),kndlog)
ALLOCATE (mask2du(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('mask2du',2,(/ncloc+1,nrloc/),kndlog)
ALLOCATE (psicatu(ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('psicatuw',3,(/ncloc+1,nrloc,nz/),kndrtype)
ALLOCATE (zdiffatw(ncloc+1,nrloc,nz+1),STAT=errstat)
CALL error_alloc('zdiffatw',3,(/ncloc+1,nrloc,nz+1/),kndrtype)

!
!2. Initialise
!-------------
!
!2.1 Mask array
!--------------
!

maskatu = nodeatu(1:ncloc+1,1:nrloc,klo:kup).EQ.2
mask2du = node2du(1:ncloc+1,1:nrloc).EQ.2

!
!2.2 Fluxes
!----------
!

fac = delt3d*hdifscal_fac
difflux = 0.0

!
!2.3 Work space arrays
!---------------------
!

IF (dYregX) THEN
   k_231: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delxatc(1:ncloc,1:nrloc)*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_231
ELSE
   k_232: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_232
ENDIF

IF (dYregX) THEN
   k_2331: DO k=klo,kup 
      WHERE (maskatu(:,:,k))
         array3du(:,:,k) = fac*delzatu(1:ncloc+1,:,k)*hdifcoef3datu(:,:,k)
      END WHERE
   ENDDO k_2331
ELSE
   k_2332: DO k=klo,kup
      WHERE (maskatu(:,:,k))
         array3du(:,:,k) = (fac*delyatu(1:ncloc+1,1:nrloc))*&
                          & delzatu(1:ncloc+1,:,k)*hdifcoef3datu(:,:,k)
      END WHERE
   ENDDO k_2332
ENDIF

!
!3. Explicit terms
!-----------------
!

ivar_300: DO ivar=1,novars
   
!
!3.1 Tracer gradients
!--------------------
!
!  ---tracer at U-nodes
   CALL Carr_at_U(psic(0:ncloc+1,1:nrloc,:,ivar),psicatu,1,2,(/0,1,1/),&
               & (/ncloc+1,nrloc,nz/),1,0,.TRUE.)

!  ---vertical tracer gradient at U-nodes
   k_310: DO k=2,nz
      WHERE (mask2du)
         zdiffatw(:,:,k) = (psicatu(:,:,k)-psicatu(:,:,k-1))&
                         & /delzatuw(1:ncloc+1,:,k)
      END WHERE
   ENDDO k_310
   WHERE (mask2du)
      zdiffatw(:,:,1) = zdiffatw(:,:,2)
      zdiffatw(:,:,nz+1) = zdiffatw(:,:,nz)
   END WHERE
   
!
!3.2 Diffusive flux at U-nodes
!-----------------------------
!

   k_320: DO k=klo,kup
      WHERE (maskatu(:,:,k))
         array2du = (psic(1:ncloc+1,1:nrloc,k,ivar)-&
                   & psic(0:ncloc,1:nrloc,k,ivar))/delxatu(1:ncloc+1,1:nrloc)
         difflux(:,:,k) = array3du(:,:,k)*(array2du-xslopeatu_geo(:,:,k)*&
                                      & 0.5*(zdiffatw(:,:,k)+zdiffatw(:,:,k+1)))
      END WHERE
   ENDDO k_320
   
!   
!3.3 Add flux divergence to explicit term
!----------------------------------------
!

   k_330: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) + &
                              & (difflux(2:ncloc+1,:,k)-difflux(1:ncloc,:,k))&
                              & /denom(:,:,k)
      END WHERE
   ENDDO k_330

ENDDO ivar_300


!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array2du,array3du,denom,difflux,maskatu,mask2du,psicatu,zdiffatw)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Xhdif_at_C_geo

!========================================================================

SUBROUTINE Xhdif_at_C_iso(psic,tridcfc,novars,klo,kup)
!************************************************************************
!
! *Xhdif_at_C_iso* Isoneutral contribution to the diffusion term in the
!                  X-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - Xdif_at_C
!
! Module calls - error_alloc
!
!************************************************************************
!
USE diffusion
USE grid
USE gridpars
USE iopars
USE physpars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                     [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*     INTEGER Number of variables
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ip, ivar, k, kkp, kp, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2du
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3du, denom, difflux, zdiff
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: slopefac


procname(pglev+1) = 'Xdif_at_C_iso'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array3du(ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('array3du',3,(/ncloc+1,nrloc,nz/),kndrtype)
ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (difflux(ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('difflux',3,(/ncloc+1,nrloc,nz/),kndrtype)
ALLOCATE (mask2du(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('mask2du',2,(/ncloc+1,nrloc/),kndlog)
ALLOCATE (slopefac(ncloc+1,nrloc,nz,0:1,0:1),STAT=errstat)
CALL error_alloc('slopefac',5,(/ncloc+1,nrloc,nz,2,2/),kndrtype)
ALLOCATE (zdiff(ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('zdiff',3,(/ncloc+1,nrloc,nz/),kndrtype)

!
!2. Initialise
!-------------
!
!---mask arrays
mask2du = node2du(1:ncloc+1,1:nrloc).EQ.2

!---work space arrays
IF (dYregX) THEN
   k_211: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delxatc(1:ncloc,1:nrloc)*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_211
ELSE
   k_212: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_212
ENDIF
k_213: DO k=1,nz
   IF (dYregX) THEN
      WHERE (mask2du)
         array3du(:,:,k) = delt3d*delzatu(1:ncloc+1,:,k)
      END WHERE
   ELSE
      WHERE (mask2du)
         array3du(:,:,k) = delt3d*delyatu(1:ncloc+1,1:nrloc)*&
                          & delzatu(1:ncloc+1,:,k)
      END WHERE
   ENDIF
END DO k_213

!
!3. Explicit terms
!-----------------
!

ivar_300: DO ivar=1,novars

!
!3.1 Vertical (off-diagonal) gradient term
!----------------------------------------
!   
!3.1.1 Initialise
!----------------
!   

   difflux = 0.0
   
   kp_311: DO kp=0,1
   ip_311: DO ip=0,1

      IF (ivar.EQ.1) THEN
         k_3111: DO k=1,nz
            kkp = k+1-kp
            IF (kkp.GT.1.AND.kkp.LE.nz) THEN
               WHERE (mask2du)
                  slopefac(:,:,k,ip,kp) = 0.25*delzatuw(:,:,k)*&
                       & (hdifscal_fac*hdifcoef3datu(:,:,k)*&
                       & xslopeatu_siso(ip:ip+ncloc,:,k,1-ip,kp)&
                       & -skewdiff_cst*xslopeatu_ziso(ip:ip+ncloc,:,k,1-ip,kp))&
                       & /delzatw(ip:ip+ncloc,1:nrloc,kkp)
               END WHERE
            ENDIF
         ENDDO k_3111
      ENDIF
               
!
!3.1.2 Add "triad" contributions
!------------------------------
!      

      k_3112: DO k=1,nz
         kkp = k+1-kp
         IF (kkp.GT.1.AND.kkp.LE.nz) THEN
            WHERE (mask2du)
               zdiff(:,:,k) = slopefac(:,:,k,ip,kp)*&
                           & (psic(ip:ip+ncloc,1:nrloc,kkp,ivar)-&
                           &  psic(ip:ip+ncloc,1:nrloc,kkp-1,ivar))
               difflux(:,:,k) = difflux(:,:,k) - zdiff(:,:,k)
            END WHERE
         ENDIF
      ENDDO k_3112
      IF (kp.EQ.0) THEN
         WHERE (mask2du)
            difflux(:,:,nz) = difflux(:,:,nz) - zdiff(:,:,nz-1)
         END WHERE
      ELSE
         WHERE (mask2du)
            difflux(:,:,1) = difflux(:,:,1) - zdiff(:,:,2)
         END WHERE
      ENDIF

   ENDDO ip_311
   ENDDO kp_311

!
!3.2 Horizontal (diagonal) gradient term
!---------------------------------------
!

   k_320: DO k=1,nz
      WHERE (mask2du)
         difflux(:,:,k) = difflux(:,:,k) + &
                & hdifscal_fac*array3du(:,:,k)*hdifcoef3datu(:,:,k)*&
                & (psic(1:ncloc+1,1:nrloc,k,ivar)-psic(0:ncloc,1:nrloc,k,ivar))&
                & /delxatu(1:ncloc+1,1:nrloc)
      END WHERE
   ENDDO k_320
      
!   
!3.3 Add flux divergence to explicit term
!----------------------------------------   
!

   k_330: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) + &
                               & (difflux(2:ncloc+1,:,k)-difflux(1:ncloc,:,k))&
                               & /denom(:,:,k)
      END WHERE
   ENDDO k_330

ENDDO ivar_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array3du,denom,difflux,mask2du,slopefac,zdiff)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Xhdif_at_C_iso

!========================================================================

SUBROUTINE Ydif_at_C(psic,tridcfc,novars,klo,kup)
!************************************************************************
!
! *Ydif_at_C* Diffusion term in the Y-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2
!
! External calls - Yhdif_at_C_lap, Yhdif_at_C_geo, Yhdif_at_C_iso
! 
! Module calls -
!
!************************************************************************
!
USE gridpars
USE switches  

IMPLICIT NONE
  
!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                     [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*     INTEGER Number of variables
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!

SELECT CASE (iopt_hdif_scal)
   CASE (1); CALL Yhdif_at_C_lap(psic,tridcfc,novars,klo,kup)
   CASE (2); CALL Yhdif_at_C_geo(psic,tridcfc,novars,klo,kup)
   CASE (3); CALL Yhdif_at_C_iso(psic,tridcfc,novars,klo,kup)
END SELECT


RETURN

END SUBROUTINE Ydif_at_C

!========================================================================

SUBROUTINE Yhdif_at_C_lap(psic,tridcfc,novars,klo,kup)
!************************************************************************
!
! *Yhdif_at_C_lap* Laplacian diffusion term in the Y-direction for a quantity
!                  at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - Ydif_at_C
!
! Module calls - error_alloc
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
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                     [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*     INTEGER Number of variables
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar, k, npcc
REAL :: fac
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2dv
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dv, denom, difflux


procname(pglev+1) = 'Ydif_at_C_lap'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2dv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (array3dv(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('array3dv',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)
ALLOCATE (difflux(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('difflux',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)
ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (maskatv(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc+1,kup-klo+1/),kndlog)
ALLOCATE (mask2dv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('mask2dv',2,(/ncloc,nrloc+1/),kndlog)

!
!2. Initialise
!-------------
!
!2.1 Mask array
!--------------
!

maskatv = nodeatv(1:ncloc,1:nrloc+1,klo:kup).EQ.2
mask2dv = node2dv(1:ncloc,1:nrloc+1).EQ.2

!
!2.2 Fluxes
!----------
!

fac = delt3d*hdifscal_fac
difflux = 0.0

!
!2.3 Work space arrays
!---------------------
!

IF (dXregY) THEN
   k_231: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delyatc(1:ncloc,1:nrloc)*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_231
ELSE
   k_232: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_232
ENDIF

IF (dXregY) THEN
   WHERE (mask2dv)
      array2dv = fac/delyatv(1:ncloc,1:nrloc+1)
   END WHERE
ELSE
   WHERE (mask2dv)
      array2dv = fac*delxatv(1:ncloc,1:nrloc+1)&
               & /delyatv(1:ncloc,1:nrloc+1)
   END WHERE
ENDIF

k_233: DO k=klo,kup
   WHERE (maskatv(:,:,k))
      array3dv(:,:,k) = array2dv*delzatv(:,1:nrloc+1,k)*hdifcoef3datv(:,:,k)
   END WHERE
END DO k_233

!
!3. Explicit terms
!-----------------
!

ivar_300: DO ivar=1,novars

!
!3.1 Fluxes
!----------
!
!

   WHERE (maskatv)
      difflux = array3dv*(psic(1:ncloc,1:nrloc+1,klo:kup,ivar)-&
                        & psic(1:ncloc,0:nrloc,klo:kup,ivar))
   END WHERE

!
!3.2 Add flux divergence to explicit term
!----------------------------------------
!

   k_320: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) + &
                               & (difflux(:,2:nrloc+1,k)-difflux(:,1:nrloc,k))&
                               & /denom(:,:,k)
      END WHERE
   ENDDO k_320

ENDDO ivar_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array2dv,array3dv,difflux,denom,maskatv,mask2dv)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Yhdif_at_C_lap
  
!========================================================================

SUBROUTINE Yhdif_at_C_geo(psic,tridcfc,novars,klo,kup)
!************************************************************************
!
! *Yhdif_at_C_geo* Isolevel contribution to the diffusion term in the
!                  Y-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - diffusion along horizontal surfaces
!
! Reference -
!
! Calling program - Ydif_at_C
!
! Module calls - Carr_at_V, error_alloc
!
!************************************************************************
!
USE diffusion
USE grid
USE gridpars
USE iopars
USE physpars
USE syspars
USE timepars
USE array_interp, ONLY: Carr_at_V
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                     [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*     INTEGER Number of variables
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar, k, npcc
REAL :: fac
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2dv
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dv, denom, difflux, &
                                           & psicatv, zdiffatw

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*difflux*    REAL    Diffusive flux (times factor) at U-nodes       
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Yhdif_at_C_geo'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2dv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (array3dv(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('array3dv',3,(/ncloc,nrloc+1,nz/),kndrtype)
ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (difflux(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('difflux',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)
ALLOCATE (maskatv(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc+1,kup-klo+1/),kndlog)
ALLOCATE (mask2dv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('mask2dv',2,(/ncloc,nrloc+1/),kndlog)
ALLOCATE (psicatv(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('psicatvw',3,(/ncloc,nrloc+1,nz/),kndrtype)
ALLOCATE (zdiffatw(ncloc,nrloc+1,nz+1),STAT=errstat)
CALL error_alloc('zdiffatw',3,(/ncloc,nrloc+1,nz+1/),kndrtype)

!
!2. Initialise
!-------------
!
!2.1 Mask array
!--------------
!

maskatv = nodeatv(1:ncloc,1:nrloc+1,klo:kup).EQ.2
mask2dv = node2dv(1:ncloc,1:nrloc+1).EQ.2

!
!2.2 Fluxes
!----------
!

fac = delt3d*hdifscal_fac
difflux = 0.0

!
!2.3 Work space arrays
!---------------------
!

IF (dXregY) THEN
   k_231: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delyatc(1:ncloc,1:nrloc)*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_231
ELSE
   k_232: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_232
ENDIF

IF (dXregY) THEN
   k_2331: DO k=klo,kup 
      WHERE (maskatv(:,:,k))
         array3dv(:,:,k) = fac*delzatv(:,1:nrloc+1,k)*hdifcoef3datv(:,:,k)
      END WHERE
   ENDDO k_2331
ELSE
   k_2332: DO k=klo,kup
      WHERE (maskatv(:,:,k))
         array3dv(:,:,k) = (fac*delxatv(1:ncloc,1:nrloc+1))*&
                          & delzatv(:,1:nrloc+1,k)*hdifcoef3datv(:,:,k)
      END WHERE
   ENDDO k_2332
ENDIF
   
!
!3. Explicit terms
!-----------------
!

ivar_300: DO ivar=1,novars
   
!
!3.1 Tracer gradients
!--------------------
!
!  ---tracer at V-nodes
   CALL Carr_at_V(psic(1:ncloc,0:nrloc+1,:,ivar),psicatv,1,2,(/1,0,1/),&
                & (/ncloc,nrloc+1,nz/),1,0,.TRUE.)

!  ---vertical tracer gradient at V-nodes
   k_310: DO k=2,nz
      WHERE (mask2dv)
         zdiffatw(:,:,k) = (psicatv(:,:,k)-psicatv(:,:,k-1))&
                         & /delzatvw(:,1:nrloc+1,k)
      END WHERE
   ENDDO k_310
   WHERE (mask2dv)
      zdiffatw(:,:,1) = zdiffatw(:,:,2)
      zdiffatw(:,:,nz+1) = zdiffatw(:,:,nz)
   END WHERE
   
!
!3.2 Diffusive flux at V-nodes
!-----------------------------
!

   k_320: DO k=klo,kup
      WHERE (maskatv(:,:,k))
         array2dv = (psic(1:ncloc,1:nrloc+1,k,ivar)-&
                   & psic(1:ncloc,0:nrloc,k,ivar))/delyatv(1:ncloc,1:nrloc+1)
         difflux(:,:,k) = array3dv(:,:,k)*(array2dv-yslopeatv_geo(:,:,k)*&
                                     & 0.5*(zdiffatw(:,:,k)+zdiffatw(:,:,k+1)))
      END WHERE
   ENDDO k_320
   
!   
!3.3 Add flux divergence to explicit term
!----------------------------------------
!

   k_330: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) + &
                              & (difflux(:,2:nrloc+1,k)-difflux(:,1:nrloc,k))&
                              & /denom(:,:,k)
      END WHERE
   ENDDO k_330

ENDDO ivar_300


!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array2dv,array3dv,denom,difflux,maskatv,mask2dv,psicatv,zdiffatw)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Yhdif_at_C_geo
  
!========================================================================

SUBROUTINE Yhdif_at_C_iso(psic,tridcfc,novars,klo,kup)
!************************************************************************
!
! *Yhdif_at_C_iso* Isoneutral contribution to the diffusion term in the
!                  Y-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - diffusion takes place along geopotential surfaces
!
! Reference -
!
! Calling program - Ydif_at_C
!
! Module calls - error_alloc
!
!************************************************************************
!
USE diffusion
USE grid
USE gridpars
USE iopars
USE physpars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                     [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*     INTEGER Number of variables
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ip, ivar, k, kkp, kp, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dv, denom, difflux, zdiff
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: slopefac


procname(pglev+1) = 'Ydif_at_C_iso'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array3dv(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('array3dv',3,(/ncloc,nrloc+1,nz/),kndrtype)
ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (difflux(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('difflux',3,(/ncloc,nrloc+1,nz/),kndrtype)
ALLOCATE (mask2dv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('mask2dv',2,(/ncloc,nrloc+1/),kndlog)
ALLOCATE (slopefac(ncloc,nrloc+1,nz,0:1,0:1),STAT=errstat)
CALL error_alloc('slopefac',5,(/ncloc,nrloc+1,nz,2,2/),kndrtype)
ALLOCATE (zdiff(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('zdiff',3,(/ncloc,nrloc+1,nz/),kndrtype)

!
!2. Initialise
!-------------
!
!---mask arrays
mask2dv = node2dv(1:ncloc,1:nrloc+1).EQ.2

!---work space arrays
IF (dXregY) THEN
   k_211: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delyatc(1:ncloc,1:nrloc)*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_211
ELSE
   k_212: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_212
ENDIF
k_213: DO k=1,nz
   IF (dXregY) THEN
      WHERE (mask2dv)
         array3dv(:,:,k) = delt3d*delzatv(:,1:nrloc+1,k)
      END WHERE
   ELSE
      WHERE (mask2dv)
         array3dv(:,:,k) = delt3d*delxatv(1:nrloc,1:nrloc+1)*&
                         & delzatv(:,1:nrloc+1,k)
      END WHERE
   ENDIF
END DO k_213

!
!3. Explicit terms
!-----------------
!

ivar_300: DO ivar=1,novars

!
!3.1 Vertical (off-diagonal) gradient term
!----------------------------------------
!   
!3.1.1 Initialise
!----------------
 !
   
   difflux = 0.0
   
   kp_311: DO kp=0,1
   ip_311: DO ip=0,1

      IF (ivar.EQ.1) THEN
         k_3111: DO k=1,nz
            kkp = k+1-kp
            IF (kkp.GT.1.AND.kkp.LE.nz) THEN
               WHERE (mask2dv)
                  slopefac(:,:,k,ip,kp) = 0.25*delzatvw(:,:,k)*&
                       & (hdifscal_fac*hdifcoef3datv(:,:,k)*&
                       & yslopeatv_siso(:,ip:ip+nrloc,k,1-ip,kp)&
                       &  -skewdiff_cst*yslopeatv_ziso(:,ip:ip+nrloc,k,1-ip,kp))
               END WHERE
            ENDIF
         ENDDO k_3111
      ENDIF
      
!
!3.1.2 Add "triad" contributions
!------------------------------
!      

      k_3112: DO k=1,nz
         kkp = k+1-kp
         IF (kkp.GT.1.AND.kkp.LE.nz) THEN
            WHERE (mask2dv)
               zdiff(:,:,k) = slopefac(:,:,k,ip,kp)*&
                           & (psic(1:ncloc,ip:ip+nrloc,kkp,ivar)-&
                            & psic(1:ncloc,ip:ip+nrloc,kkp-1,ivar))
               difflux(:,:,k) = difflux(:,:,k) - zdiff(:,:,k)
            END WHERE
         ENDIF
      ENDDO k_3112
      IF (kp.EQ.0) THEN
         WHERE (mask2dv)
            difflux(:,:,nz) = difflux(:,:,nz) - zdiff(:,:,nz-1)
         END WHERE
         
      ELSE
         WHERE (mask2dv)
            difflux(:,:,1) = difflux(:,:,1) - zdiff(:,:,2)
         END WHERE
      ENDIF

   ENDDO ip_311
   ENDDO kp_311

!
!3.2 Horizontal (diagonal) gradient term
!---------------------------------------
!

   k_320: DO k=1,nz
      WHERE (mask2dv)
         difflux(:,:,k) = difflux(:,:,k) + &
                & hdifscal_fac*array3dv(:,:,k)*hdifcoef3datv(:,:,k)*&
                & (psic(1:ncloc,1:nrloc+1,k,ivar)-psic(1:ncloc,0:nrloc,k,ivar))&
                & /delyatv(1:ncloc,1:nrloc+1)
      END WHERE
   ENDDO k_320
      
!   
!3.3 Add flux divergence to explicit term
!----------------------------------------   
!

   k_330: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) + &
                               & (difflux(:,2:nrloc+1,k)-difflux(:,1:nrloc,k))&
                               & /denom(:,:,k)
      END WHERE
   ENDDO k_330

ENDDO ivar_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array3dv,denom,difflux,mask2dv,slopefac,zdiff)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Yhdif_at_C_iso

!========================================================================

SUBROUTINE Zdif_at_C(psic,tridcfc,vdifcoefatw,novars,ibcsur,ibcbot,nbcs,nbcb,&
                   & bcsur,bcbot,klo,kup)
!************************************************************************
!
! *Zdif_at_C* Diffusion term in the  vertical direction for a quantity
!             at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2
!
! Module calls -
!
!************************************************************************
!
USE gridpars
USE switches

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ibcbot, ibcsur, klo, kup, nbcb, nbcs, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1,novars) :: vdifcoefatw
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcs,novars) :: bcsur
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcb,novars) :: bcbot

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                    [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit vertical solution
!*vdifcoefatw*REAL    Diffusion coefficient at W-nodes                 [m^2/s]
!*novars*     INTEGER Number of variables
!*ibcsur*     INTEGER Type of surface boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!                 = 3 => Dirichlet at the surface
!*ibcbot*     INTEGER Type of bottom boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!                 = 3 => Dirichlet at the bottom
!*nbcs*       INTEGER Last dimension of array bcsur
!*nbcb*       INTEGER Last dimension of array bcbot
!*bcsur*      REAL    Data for surface boundary condition
!              ibcsur = 1 => (:,:,1,:) = prescribed surface flux      [psic m/s]
!              ibcsur = 2 => (:,:,1,:) = prescribed surface value         [psic]
!                            (:,:,2,:) = transfer velocity                 [m/s]
!              ibcsur = 3 => (:,:,1,:) = explicit part of the surface flux
!                                                                     [psic m/s]
!                            (:,:,2,:) = implicit part of the surface flux
!                                                                     [psic m/s]
!*bcbot*      REAL    Data for bottom boundary condition
!              ibcbot = 1 => (:,:,1,:) = prescribed bottom flux       [psic m/s]
!              ibcbot = 2 => (:,:,1,:) = prescribed bottom value          [psic]
!                            (:,:,2,:) = transfer velocity                 [m/s]
!              ibcbot = 3 => (:,:,1,:) = explicit part of the bottom flux
!                                                                     [psic m/s]
!                            (:,:,2,:) = implicit part of the bottom flux
!                                                                     [psic m/s]
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: vdiff


SELECT CASE (iopt_hdif_scal)
   CASE(2); CALL Zhdif_at_C_geo(psic,tridcfc,novars,klo,kup)
   CASE(3); CALL Zhdif_at_C_iso(psic,tridcfc,novars,klo,kup)
END SELECT

vdiff = iopt_vdif_coef.GT.0.OR.&
      & (iopt_grid_nodim.EQ.2.AND.(ibcbot.GT.0.OR.ibcsur.GT.0))
vdiff = vdiff.AND.iopt_vdif_impl.GT.0

IF (vdiff) THEN  
   CALL Zvdif_at_C(psic,tridcfc,vdifcoefatw,novars,ibcsur,ibcbot,nbcs,nbcb,&
                 & bcsur,bcbot,klo,kup)
ENDIF


RETURN

END SUBROUTINE Zdif_at_C
  
!========================================================================

SUBROUTINE Zhdif_at_C_geo(psic,tridcfc,novars,klo,kup)
!************************************************************************
!
! *Zhdif_at_C_geo* Isolevel contribution to the diffusion term in the
!                  Z-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - Zdif_at_C
!
! Module calls - Carr_at_UW, UWarr_at_W, VWarr_at_W, error_alloc
!
!************************************************************************
!
USE diffusion
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE syspars
USE timepars
USE array_interp, ONLY: Carr_at_W, UWarr_at_W, VWarr_at_W
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                     [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*     INTEGER Number of variables
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar, k, npcc
REAL :: fac
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dx, array3dy, difflux, &
                                           & psiatw, xdiff, ydiff, wdiff


procname(pglev+1) = 'Zhdif_at_C_geo'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array3dx(ncloc,nrloc,nz-1),STAT=errstat)
CALL error_alloc('array3dx',3,(/ncloc,nrloc,nz-1/),kndrtype)
ALLOCATE (array3dy(ncloc,nrloc,nz-1),STAT=errstat)
CALL error_alloc('array3dy',3,(/ncloc,nrloc,nz-1/),kndrtype)
ALLOCATE (difflux(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('difflux',3,(/ncloc,nrloc,nz+1/),kndrtype)
ALLOCATE (psiatw(0:ncloc+1,0:nrloc+1,nz-1),STAT=errstat)
CALL error_alloc('psiatuw',3,(/ncloc+2,nrloc+2,nz-1/),kndrtype)
ALLOCATE (xdiff(ncloc+1,nrloc,nz+1),STAT=errstat)
CALL error_alloc('xdiff',3,(/ncloc+1,nrloc,nz+1/),kndrtype)
ALLOCATE (ydiff(ncloc,nrloc+1,nz+1),STAT=errstat)
CALL error_alloc('ydiff',3,(/ncloc,nrloc+1,nz+1/),kndrtype)
ALLOCATE (wdiff(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('wdiff',3,(/ncloc,nrloc,nz+1/),kndrtype)

!
!2. Initialise
!-------------
!
!---work space arrays
difflux = 0.0; psiatw = 0.0; xdiff = 0.0; ydiff = 0.0; wdiff = 0.0
fac = delt3d*hdifscal_fac

!---work space array
k_210: DO k=2,nz
   WHERE (maskatc_int)
      array3dx(:,:,k-1) = fac*hdifcoef3datw(:,:,k)*&
                         & xslopeatw_geo(:,:,k)
      array3dy(:,:,k-1) = fac*hdifcoef3datw(:,:,k)*&
                         & yslopeatw_geo(:,:,k)
   END WHERE
ENDDO k_210

!
!3. Explicit terms
!-----------------
!

ivar_300: DO ivar=1,novars
   
!
!3.1 X-direction
!---------------
!   
!  ---tracer at W-nodes
   CALL Carr_at_W(psic(0:ncloc+1,0:nrloc+1,:,ivar),psiatw,(/0,0,1/),&
               & (/ncloc+1,nrloc+1,nz/),1,0,.TRUE.)
!  ---X-gradient   
   k_311: DO k=2,nz
      WHERE (node2du(1:ncloc+1,1:nrloc).EQ.2)
         xdiff(:,:,k) = (psiatw(1:ncloc+1,1:nrloc,k-1)-&
                       & psiatw(0:ncloc,1:nrloc,k-1))&
                       & /delxatu(1:ncloc+1,1:nrloc)
      END WHERE
   ENDDO k_311
   CALL UWarr_at_W(xdiff,wdiff,2,(/1,1,1/),(/ncloc+1,nrloc,nz+1/),1,0,.TRUE.)
      
!  ---add to diffusive flux
   k_312: DO k=2,nz
      WHERE (maskatc_int)
         difflux(:,:,k) = -array3dx(:,:,k-1)*wdiff(:,:,k)
      END WHERE
   ENDDO k_312

!
!3.2 Y-direction
!---------------
!   
!  ---Y-gradient
   k_321: DO k=2,nz
      WHERE (node2dv(1:ncloc,1:nrloc+1).EQ.2)
         ydiff(:,:,k) = (psiatw(1:ncloc,1:nrloc+1,k-1)-&
                       & psiatw(1:ncloc,0:nrloc,k-1))&
                       & /delyatv(1:ncloc,1:nrloc+1)
      END WHERE
   ENDDO k_321
   CALL VWarr_at_W(ydiff,wdiff,2,(/1,1,1/),(/ncloc,nrloc+1,nz+1/),1,0,.TRUE.)
   
!  ---add to diffusive flux
   k_322: DO k=2,nz
      WHERE (maskatc_int)
         difflux(:,:,k) = difflux(:,:,k) - array3dy(:,:,k-1)*wdiff(:,:,k)
      END WHERE
   ENDDO k_322
   
!
!3.3 Add flux divergence to explicit term
!----------------------------------------   
!

   k_330: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) + &
                  & (difflux(:,:,k+1)-difflux(:,:,k))/delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_330

ENDDO ivar_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array3dx,array3dy,difflux,psiatw,xdiff,ydiff,wdiff)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Zhdif_at_C_geo
  
!========================================================================

SUBROUTINE Zhdif_at_C_iso(psic,tridcfc,novars,klo,kup)
!************************************************************************
!
! *Zhdif_at_C_iso* Isoneutral contribution to the diffusion term in the
!                  Z-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - Zdif_at_C
!
! Module calls - error_alloc
!
!************************************************************************
!
USE diffusion
USE grid
USE gridpars
USE iopars
USE physpars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                      [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*     INTEGER Number of variables
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!
!*Local variables
!
INTEGER :: ip, ipp, ivar, k, kkp, kp, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2du, mask2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2du, array2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: denom, difflux
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: slopefacx, slopefacy


procname(pglev+1) = 'Zhdif_at_C_iso'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2du(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (denom(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (difflux(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('difflux',3,(/ncloc,nrloc,nz+1/),kndrtype)
ALLOCATE (mask2du(ncloc,nrloc),STAT=errstat)
CALL error_alloc('mask2du',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (mask2dv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('mask2dv',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (slopefacx(ncloc,nrloc,nz-1,0:1,0:1),STAT=errstat)
CALL error_alloc('slopefacx',5,(/ncloc,nrloc,nz-1,2,2/),kndrtype)
ALLOCATE (slopefacy(ncloc,nrloc,nz-1,0:1,0:1),STAT=errstat)
CALL error_alloc('slopefacy',5,(/ncloc,nrloc,nz-1,2,2/),kndrtype)

!
!2. Initialise
!-------------
!
!---mask arrays
mask2du = node2du(1:ncloc,1:nrloc).EQ.2
mask2dv = node2dv(1:ncloc,1:nrloc).EQ.2

!---work space array
k_210: DO k=1,nz
   WHERE (maskatc_int)
      denom(:,:,k) = garea*delzatc(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_210

!
!3. Explicit terms
!-----------------
!

ivar_300: DO ivar=1,novars

!
!3.1 Initialise
!--------------
!

   IF (ivar.EQ.1) THEN
      
      kp_310: DO kp=0,1
      ip_310: DO ip=0,1
         ipp = 1-ip
         array2du = 0.25*delxatu(ipp+1:ipp+ncloc,1:nrloc)*&
                 &  delyatu(ipp+1:ipp+ncloc,1:nrloc)
         array2dv = 0.25*delxatv(1:ncloc,ipp+1:ipp+nrloc)*&
                 &  delyatv(1:ncloc,ipp+1:ipp+nrloc)
         k_311: DO k=2,nz
            kkp = k+kp-1
            WHERE (mask2du)
               slopefacx(:,:,k-1,ip,kp) = array2du*&
                          & (hdifscal_fac*hdifcoef3datu(ipp+1:ipp+ncloc,:,k)*&
                          & xslopeatu_siso(1:ncloc,:,kkp,1-ip,kp)&
                          & +skewdiff_cst*xslopeatu_ziso(1:ncloc,:,kkp,1-ip,kp))
            END WHERE
            WHERE (mask2dv)
               slopefacy(:,:,k-1,ip,kp) = array2dv*&
                          & (hdifscal_fac*hdifcoef3datv(:,ipp+1:ipp+nrloc,k)*&
                          & yslopeatv_siso(:,1:nrloc,kkp,1-ip,kp)&
                          & +skewdiff_cst*yslopeatv_ziso(:,1:nrloc,kkp,1-ip,kp))
            END WHERE
         ENDDO k_311

      ENDDO ip_310
      ENDDO kp_310
      
   ENDIF

!
!3.2 Diffusive flux
!------------------
!

   difflux = 0.0
   
   kp_320: DO kp=0,1
   ip_320: DO ip=0,1
      
      k_321: DO k=2,nz
         ipp = 1-ip; kkp = k+kp-1

!        ---X-direction
         WHERE (mask2du)
            array2du = (psic(ipp+1:ipp+ncloc,1:nrloc,kkp,ivar)-&
                      & psic(ipp:ipp+ncloc-1,1:nrloc,kkp,ivar))&
                      & /delxatu(ipp+1:ipp+ncloc,1:nrloc)
            difflux(:,:,k) = difflux(:,:,k) - slopefacx(:,:,k-1,ip,kp)*array2du
         END WHERE

!        ---y-direction
         WHERE (mask2dv)
            array2dv = (psic(1:ncloc,ipp+1:ipp+nrloc,kkp,ivar)-&
                      & psic(1:ncloc,ipp:ipp+nrloc-1,kkp,ivar))&
                      & /delyatv(1:ncloc,ipp+1:ipp+nrloc)
            difflux(:,:,k) = difflux(:,:,k) - slopefacy(:,:,k-1,ip,kp)*array2dv
         END WHERE

      ENDDO k_321

   ENDDO ip_320
   ENDDO kp_320

!
!3.3 Add flux divergence to explicit term
!----------------------------------------   
!

   k_330: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) + &
                          & delt3d*(difflux(:,:,k+1)-difflux(:,:,k))&
                          & /(garea*delzatc(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_330

ENDDO ivar_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array2du,array2dv,denom,difflux,mask2du,mask2dv,slopefacx,&
          & slopefacy)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Zhdif_at_C_iso
   
!========================================================================

SUBROUTINE Zvdif_at_C(psic,tridcfc,vdifcoefatw,novars,ibcsur,ibcbot,nbcs,nbcb,&
                    & bcsur,bcbot,klo,kup)
!************************************************************************
!
! *Zvdif_at_C* Vertical diffusion term for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2
!
! Module calls - error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ibcbot, ibcsur, klo, kup, nbcb, nbcs, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1,novars) :: vdifcoefatw
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcs,novars) :: bcsur
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcb,novars) :: bcbot

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psic*       REAL    C-node quantity to be diffused                    [psic]
!*tridcfc*    REAL    Tridiagonal matrix for implicit vertical solution
!*vdifcoefatw*REAL    Diffusion coefficient at W-nodes                 [m^2/s]
!*novars*     INTEGER Number of variables
!*ibcsur*     INTEGER Type of surface boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!                 = 3 => Dirichlet at the surface
!*ibcbot*     INTEGER Type of bottom boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!                 = 3 => Dirichlet at the bottom
!*nbcs*       INTEGER Last dimension of array bcsur
!*nbcb*       INTEGER Last dimension of array bcbot
!*bcsur*      REAL    Data for surface boundary condition
!              ibcsur = 1 => (:,:,1,:) = prescribed surface flux      [psic m/s]
!              ibcsur = 2 => (:,:,1,:) = prescribed surface value         [psic]
!                            (:,:,2,:) = transfer velocity                 [m/s]
!              ibcsur = 3 => (:,:,1,:) = explicit part of the surface flux
!                                                                     [psic m/s]
!                            (:,:,2,:) = implicit part of the surface flux
!                                                                     [psic m/s]
!*bcbot*      REAL    Data for bottom boundary condition
!              ibcbot = 1 => (:,:,1,:) = prescribed bottom flux       [psic m/s]
!              ibcbot = 2 => (:,:,1,:) = prescribed bottom value          [psic]
!                            (:,:,2,:) = transfer velocity                 [m/s]
!              ibcbot = 3 => (:,:,1,:) = explicit part of the bottom flux
!                                                                     [psic m/s]
!                            (:,:,2,:) = implicit part of the bottom flux
!                                                                     [psic m/s]
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ivar, k, kmax, kmin, npcc
REAL :: theta_vdif1, xexp, ximp
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dw, difflux

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*difflux*   REAL    Diffusive flux (times factor) at W-nodes       [m^2/psic]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Zvdif_at_C'
CALL log_timer_in(npcc)

!
!1. Set parameters
!------------------
!

xexp = delt3d*(1.0-theta_vdif)
ximp = delt3d*theta_vdif
theta_vdif1 = 1.0-theta_vdif

!
!2. Allocate arrays
!------------------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array3dw(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('array3dw',3,(/ncloc,nrloc,nz-1/),kndrtype)
IF (iopt_vdif_impl.NE.2) THEN
   ALLOCATE (difflux(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('difflux',3,(/ncloc,nrloc,nz+1/),kndrtype)
ENDIF

!
!3. Vertical diffusion term
!--------------------------
!

ivar_300: DO ivar=1,novars

!
!3.1 Work space array
!--------------------
!

   k_310: DO k=2,nz
      WHERE (maskatc_int)
         array3dw(:,:,k) = vdifcoefatw(:,:,k,ivar)/delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_310

!
!3.2 Explicit fluxes
!-------------------
!

   IF (iopt_vdif_impl.NE.2) THEN
      k_320: DO k=2,nz
         WHERE (maskatc_int)
            difflux(:,:,k) = xexp*array3dw(:,:,k)*&
                           & (psic(1:ncloc,1:nrloc,k,ivar)-&
                            & psic(1:ncloc,1:nrloc,k-1,ivar))
         END WHERE
      ENDDO k_320
   ENDIF

!
!3.3 Explicit terms
!------------------
!

   IF (iopt_vdif_impl.NE.2) THEN
      k_330: DO k=klo,kup
         WHERE (maskatc_int)
            tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) + &
                             & (difflux(:,:,k+1)-difflux(:,:,k))&
                             & /delzatc(1:ncloc,1:nrloc,k)
         END WHERE
      ENDDO k_330
   ENDIF

!
!3.4 Implicit terms
!------------------
!

   IF (iopt_vdif_impl.NE.0) THEN

!     ---lower flux
      kmax = MERGE(nz,nz-1,ibcsur.LE.2)
      k_341: DO k=2,kmax
         WHERE (maskatc_int)
            array2dc = ximp*array3dw(:,:,k)/delzatc(1:ncloc,1:nrloc,k)
            tridcfc(:,:,k,1,ivar) = tridcfc(:,:,k,1,ivar) - array2dc
            tridcfc(:,:,k,2,ivar) = tridcfc(:,:,k,2,ivar) + array2dc
         END WHERE
      ENDDO k_341

!     ---upper flux
      kmin = MERGE(1,2,ibcbot.LE.2)
      k_342: DO k=kmin,nz-1
         WHERE (maskatc_int)
            array2dc = ximp*array3dw(:,:,k+1)/delzatc(1:ncloc,1:nrloc,k)
            tridcfc(:,:,k,2,ivar) = tridcfc(:,:,k,2,ivar) + array2dc
            tridcfc(:,:,k,3,ivar) = tridcfc(:,:,k,3,ivar) - array2dc
         END WHERE
      ENDDO k_342

   ENDIF

!
!3.5 Boundary conditions
!----------------------
!
!3.5.1 Surface
!-------------
!
!  ---Neumann (prescribed flux)
   IF (ibcsur.EQ.1) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,nz,4,ivar) = tridcfc(:,:,nz,4,ivar) + delt3d*&
                         & bcsur(:,:,1,ivar)/delzatc(1:ncloc,1:nrloc,nz)
      END WHERE

!  ---Neumann (using transfer velocity)
   ELSEIF (ibcsur.EQ.2) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,nz,2,ivar) = tridcfc(:,:,nz,2,ivar) + &
                         & ximp*bcsur(:,:,2,ivar)/delzatc(1:ncloc,1:nrloc,nz)
         tridcfc(:,:,nz,4,ivar) = tridcfc(:,:,nz,4,ivar) - &
                         & bcsur(:,:,2,ivar)*&
                         & (xexp*psic(1:ncloc,1:nrloc,nz,ivar)-&
                         & delt3d*bcsur(:,:,1,ivar))/delzatc(1:ncloc,1:nrloc,nz)
      END WHERE
      
!  ---Dirichlet
   ELSEIF (ibcsur.EQ.3) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,nz,2,ivar) = tridcfc(:,:,nz,2,ivar) - &
                         & ximp*bcsur(:,:,2,ivar)/delzatc(1:ncloc,1:nrloc,nz)
         tridcfc(:,:,nz,4,ivar) = tridcfc(:,:,nz,4,ivar) + &
                         & xexp*bcsur(:,:,1,ivar)/delzatc(1:ncloc,1:nrloc,nz)
      END WHERE
   ENDIF

!
!3.5.2 Bottom
!------------
!
!  ---Neumann (prescribed flux)
   IF (ibcbot.EQ.1) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,1,4,ivar) = tridcfc(:,:,1,4,ivar) - &
                          & delt3d*bcbot(:,:,1,ivar)/delzatc(1:ncloc,1:nrloc,1)
      END WHERE

!  ---Neumann (using transfer velocity)
   ELSEIF (ibcbot.EQ.2) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,1,2,ivar) = tridcfc(:,:,1,2,ivar) + &
                          & ximp*bcbot(:,:,2,ivar)/delzatc(1:ncloc,1:nrloc,1)
         tridcfc(:,:,1,4,ivar) = tridcfc(:,:,1,4,ivar) - &
                          & bcbot(:,:,2,ivar)*&
                          & (xexp*psic(1:ncloc,1:nrloc,1,ivar)-&
                          & delt3d*bcbot(:,:,1,ivar))/delzatc(1:ncloc,1:nrloc,1)
      END WHERE

!  ---Dirichlet
   ELSEIF (ibcbot.EQ.3) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,1,2,ivar) = tridcfc(:,:,1,2,ivar) + &
                          & delt3d*bcbot(:,:,2,ivar)/delzatc(1:ncloc,1:nrloc,1)
         tridcfc(:,:,1,4,ivar) = tridcfc(:,:,1,4,ivar) - &
                          & delt3d*bcbot(:,:,1,ivar)/delzatc(1:ncloc,1:nrloc,1)
      END WHERE
   ENDIF

ENDDO ivar_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array2dc,array3dw)
IF (iopt_vdif_impl.NE.2) DEALLOCATE (difflux)

CALL log_timer_out(npcc,itm_vdif)


RETURN

END SUBROUTINE Zvdif_at_C

!========================================================================

SUBROUTINE Xdif_at_U_3d(psiu,xdifatu3d,vintflag)
!************************************************************************
!
! *Xdif_at_U_3d* Diffusion term in X-direction for the 3-D current at U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_U_3d
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE diffusion  
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: vintflag
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiu
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,nz) :: xdifatu3d

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psiu*        REAL    U-node quantity to be diffused                    [psiu]
!*xdifatu3d*   REAL    Diffusion term (times delz)                   [m*psiu/s]
!*vintflag*    LOGICAL Update vertically integrated diffusive term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc1, array2dc2, denom, &
                                         & difflux, unorm, vnorm

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*difflux*  REAL    Diffusive flux at C-nodes
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xdif_at_U_3d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatc(0:ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc+1,nrloc/),kndlog)
ALLOCATE (maskatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (array2dc1(0:ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc1',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (array2dc2(0:ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc2',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (difflux(0:ncloc,nrloc),STAT=errstat)
CALL error_alloc('difflux',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (unorm(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('unorm',2,(/ncloc+2,nrloc/),kndrtype)
ALLOCATE (vnorm(0:ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('vnorm',2,(/ncloc+1,nrloc+1/),kndrtype)

!
!2. Initialise
!-------------
!
!2.1 Mask arrays
!---------------
!

maskatu = nodeatu(1:ncloc,1:nrloc,:).EQ.2
maskatc = nodeatc(0:ncloc,1:nrloc).GT.0

!
!2.2 Fluxes
!----------
!

WHERE (.NOT.maskatc) difflux = 0.0

!
!2.3 Work space arrays
!---------------------
!

IF (dYregX) THEN
   denom = delxatu(1:ncloc,1:nrloc)
ELSE
   denom = delxatu(1:ncloc,1:nrloc)*delyatu(1:ncloc,1:nrloc)**2
ENDIF
WHERE (maskatc)
   array2dc1 = delyatc(0:ncloc,1:nrloc)/delxatc(0:ncloc,1:nrloc)
END WHERE

!
!3. Explicit terms
!-----------------
!

k_300: DO k=1,nz

!
!3.1 Fluxes
!----------
!
!  ---work space arrays
   WHERE (nodeatu(0:ncloc+1,1:nrloc,k).GT.1)
      unorm = psiu(0:ncloc+1,1:nrloc,k)/delyatu(0:ncloc+1,1:nrloc)
   ELSEWHERE
      unorm = 0.0
   END WHERE
   WHERE (nodeatv(0:ncloc,1:nrloc+1,k).GT.1)
      vnorm = vvel_old(0:ncloc,1:nrloc+1,k)/delxatv(0:ncloc,1:nrloc+1)
   ELSEWHERE
      vnorm = 0.0
   END WHERE
   IF (dYregX) THEN
      WHERE (maskatc)
         array2dc2 = hdifmom_fac*hdifcoef3datc(:,1:nrloc,k)*&
                   & delzatc(0:ncloc,1:nrloc,k)
      END WHERE
   ELSE
      WHERE (maskatc)
         array2dc2 = hdifmom_fac*hdifcoef3datc(:,1:nrloc,k)*&
                   & delzatc(0:ncloc,1:nrloc,k)*delyatc(0:ncloc,1:nrloc)**2
      END WHERE
   ENDIF
   WHERE (maskatc)
      difflux = array2dc2*(array2dc1*(unorm(1:ncloc+1,:)-unorm(0:ncloc,:))-&
                        & (vnorm(:,2:nrloc+1)-vnorm(:,1:nrloc))/array2dc1)
   END WHERE

!
!3.3 Add flux divergence to explicit term
!----------------------------------------
!

   WHERE (maskatu(:,:,k))
      xdifatu3d(:,:,k) = (difflux(1:ncloc,:)-difflux(0:ncloc-1,:))/denom
   END WHERE

!
!3.4 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

   IF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatu(:,:,k))
         xdifatu3d(:,:,k) = alphatu_fld*xdifatu3d(:,:,k)
      END WHERE
   ENDIF

!
!3.5 Update depth-integrated diffusive term
!------------------------------------------
!

   IF (vintflag) THEN
      WHERE (maskatu(:,:,k))
         udevint = udevint + xdifatu3d(:,:,k)
      END WHERE
   ENDIF

ENDDO k_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatc,maskatu)
DEALLOCATE (array2dc1,array2dc2,denom,difflux,unorm,vnorm)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Xdif_at_U_3d

!========================================================================

SUBROUTINE Xdif_at_U_2d(psiu,xdifatu2d,vintflag)
!************************************************************************
!
! *Xdif_at_U_2d* Diffusion term in X-direction for the 2-D current at U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_U_2d
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE depths
USE diffusion
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(INOUT) :: vintflag
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo) :: psiu
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: xdifatu2d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psiu*        REAL    U-node quantity to be diffused                    [psiu]
!*xdifatu2d*   REAL    Diffusion term                                  [psiu/s]
!*hdifcoefatc* REAL    Diffusion coefficient at C-nodes                 [m^2/s]
!*vintflag*    LOGICAL Update vertically integrated diffusive term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc, maskatu
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, difflux, unorm, vnorm

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*difflux*   REAL    Diffusive flux at C-nodes
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xdif_at_U_2d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatc(0:ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc+1,nrloc/),kndlog)
ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (array2dc(0:ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (difflux(0:ncloc,nrloc),STAT=errstat)
CALL error_alloc('difflux',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (unorm(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('unorm',2,(/ncloc+2,nrloc/),kndrtype)
ALLOCATE (vnorm(0:ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('vnorm',2,(/ncloc+1,nrloc+1/),kndrtype)

!
!2. Initialise
!-------------
!
!2.1 Mask arrays
!---------------
!

maskatu = node2du(1:ncloc,1:nrloc).EQ.2
maskatc = nodeatc(0:ncloc,1:nrloc).GT.0

!
!2.2 Fluxes
!----------
!

WHERE (.NOT.maskatc)
   difflux = 0.0
END WHERE

!
!2.3 Work space arrays
!---------------------
!

WHERE (node2du(0:ncloc+1,1:nrloc).GT.1)
   unorm = psiu(0:ncloc+1,1:nrloc)&
         & /(delyatu(0:ncloc+1,1:nrloc)*deptotatu(0:ncloc+1,1:nrloc))
ELSEWHERE
   unorm = 0.0
END WHERE
WHERE (node2dv(0:ncloc,1:nrloc+1).GT.1)
   vnorm = vmvel_old(0:ncloc,1:nrloc+1)/delxatv(0:ncloc,1:nrloc+1)
ELSEWHERE
   vnorm = 0.0
END WHERE
WHERE (maskatc) 
   array2dc = delyatc(0:ncloc,1:nrloc)/delxatc(0:ncloc,1:nrloc)
END WHERE

!
!3. Explicit terms
!-----------------
!
!3.1 Fluxes
!----------
!
!---regular grid
IF (dYregX) THEN
   WHERE (maskatc)
      difflux = hdifmom_fac*hdifcoef2datc(:,1:nrloc)*&
              & (array2dc*(unorm(1:ncloc+1,:)-unorm(0:ncloc,:))-&
              & (vnorm(:,2:nrloc+1)-vnorm(:,1:nrloc))/array2dc)
   END WHERE

!---irregular grid
ELSE
   WHERE (maskatc)
      difflux = hdifmom_fac*hdifcoef2datc(:,1:nrloc)*&
              & (delyatc(0:ncloc,1:nrloc)**2)*&
              & (array2dc*(unorm(1:ncloc+1,:)-unorm(0:ncloc,:))-&
              & (vnorm(:,2:nrloc+1)-vnorm(:,1:nrloc))/array2dc)
   END WHERE
ENDIF

!
!3.2 Add flux divergence to explicit term
!----------------------------------------
!

IF (dYregX) THEN
   WHERE (maskatu)
      xdifatu2d = (difflux(1:ncloc,:)-difflux(0:ncloc-1,:))&
                & /delxatu(1:ncloc,1:nrloc)
   END WHERE
ELSE
   WHERE (maskatu)
      xdifatu2d = (difflux(1:ncloc,:)-difflux(0:ncloc-1,:))&
                & /(delxatu(1:ncloc,1:nrloc)*delyatu(1:ncloc,1:nrloc)**2)
   END WHERE
ENDIF

!
!3.3 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN
   WHERE (maskatu)
      xdifatu2d = alphatu_fld*xdifatu2d
   END WHERE
ENDIF

!
!3.4 Update depth-integrated diffusive term
!------------------------------------------
!

IF (vintflag) THEN
   WHERE (maskatu)
      udevint = udevint - xdifatu2d
   END WHERE
   vintflag = .FALSE.
ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatc,maskatu)
DEALLOCATE (array2dc,difflux,unorm,vnorm)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Xdif_at_U_2d

!========================================================================

SUBROUTINE Ydif_at_U_3d(psiu,ydifatu3d,vintflag)
!************************************************************************
!
! *Ydif_at_U_3d* Diffusion term in Y-direction for the 3-D current at U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_U_3d
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE diffusion  
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: vintflag
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiu
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,nz) :: ydifatu3d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psiu*        REAL    U-node quantity to be diffused                    [psiu]
!*ydifatu3d*   REAL    Diffusion term (times delz)                   [m*psiu/s]
!*vintflag*    LOGICAL Update vertically integrated diffusive term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, jj, jjloc, ju, jv, k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2duv
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu, maskatuv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2duv1, array2duv2, array2duv3, &
                                         & denom, difflux, unorm, vnorm

!
! Name       Type   Purpose
!------------------------------------------------------------------------------
!*difflux*   REAL  Diffusive flux at UV-nodes
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Ydif_at_U_3d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (maskatuv(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('maskatuv',3,(/ncloc,nrloc+1,nz/),kndlog)
ALLOCATE (mask2duv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('mask2duv',2,(/ncloc,nrloc+1/),kndlog)
ALLOCATE (array2duv1(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2duv1',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (array2duv2(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2duv2',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (difflux(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('difflux',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (unorm(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('unorm',2,(/ncloc,nrloc+2/),kndrtype)
ALLOCATE (vnorm(0:ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('vnorm',2,(/ncloc+1,nrloc+1/),kndrtype)
IF (.NOT.dXregY) THEN
   ALLOCATE (array2duv3(ncloc,nrloc+1),STAT=errstat)
   CALL error_alloc('array2duv3',2,(/ncloc,nrloc+1/),kndrtype)
ENDIF

!
!2. Initialise
!-------------
!
!2.1 Mask arrays
!---------------
!

maskatu = nodeatu(1:ncloc,1:nrloc,:).EQ.2
maskatuv = nodeatuv(1:ncloc,1:nrloc+1,:).GT.0
mask2duv = node2duv(1:ncloc,1:nrloc+1).GT.0

!
!2.2 Work space arrays
!---------------------
!

IF (dXregY) THEN
   denom = delyatu(1:ncloc,1:nrloc)
ELSE
   denom = delyatu(1:ncloc,1:nrloc)*delxatu(1:ncloc,1:nrloc)**2
   WHERE (mask2duv)
      array2duv3 = delxatuv(1:ncloc,1:nrloc+1)**2
   END WHERE
ENDIF
WHERE (mask2duv)
   array2duv1 = delxatuv(1:ncloc,1:nrloc+1)/delyatuv(1:ncloc,1:nrloc+1)
END WHERE

!
!3. Diffusion term
!-----------------
!

k_300: DO k=1,nz

!
!3.1 Fluxes at interior nodes
!----------------------------
!

   WHERE (nodeatu(1:ncloc,0:nrloc+1,k).GT.1) 
      unorm = psiu(1:ncloc,0:nrloc+1,k)/delxatu(1:ncloc,0:nrloc+1)
   ELSEWHERE
      unorm = 0.0
   END WHERE
   WHERE (nodeatv(0:ncloc,1:nrloc+1,k).GT.1)
      vnorm = vvel_old(0:ncloc,1:nrloc+1,k)/delyatv(0:ncloc,1:nrloc+1)
   ELSEWHERE
      vnorm = 0.0
   END WHERE
   IF (dXregY) THEN
      WHERE (maskatuv(:,:,k))
         array2duv2 = hdifmom_fac*hdifcoef3datuv(1:ncloc,:,k)*&
                    & delzatuv(1:ncloc,:,k)
      END WHERE
   ELSE
      WHERE (maskatuv(:,:,k))
         array2duv2 = hdifmom_fac*hdifcoef3datuv(1:ncloc,:,k)*&
                    & delzatuv(1:ncloc,:,k)*array2duv3
      END WHERE
   ENDIF
   WHERE (maskatuv(:,:,k))
      difflux = array2duv2*(array2duv1*(unorm(:,1:nrloc+1)-unorm(:,0:nrloc))+&
                         & (vnorm(1:ncloc,:)-vnorm(0:ncloc-1,:))/array2duv1)
   ELSEWHERE
      difflux = 0.0
   END WHERE

!
!3.2 Open boundary fluxes
!------------------------
!

   jjloc_320: DO jjloc=1,nobyloc_ext
      i = iobyloc(jjloc); j = jobyloc(jjloc)
      jj = indexoby(jjloc); ju = MERGE(j-1,j,soutoby(jj))
      IF (nodeatu(i,ju,k).GT.2.OR.nodeatv(i-1,j,k).LE.1.OR.&
        & nodeatv(i,j,k).LE.1) THEN
         difflux(i,j) = array2duv2(i,j)*(array2duv1(i,j)*&
                      & (unorm(i,j)-unorm(i,j-1))+&
                      & (vnorm(i,j)-vnorm(i-1,j))/array2duv1(i,j))
      ELSE
         jv = MERGE(j+1,j-1,soutoby(jj))
         IF (jv.GE.1.AND.jv.LE.nrloc+1.AND.nodeatu(i,ju,k).EQ.2) THEN
            difflux(i,j) = difflux(i,jv)
         ENDIF
      ENDIF
   ENDDO jjloc_320

!
!3.3 Add flux divergence to explicit term
!----------------------------------------
!

   WHERE (maskatu(:,:,k))
      ydifatu3d(:,:,k) = (difflux(:,2:nrloc+1)-difflux(:,1:nrloc))/denom
   END WHERE

!
!3.4 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

   IF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatu(:,:,k))
         ydifatu3d(:,:,k) = alphatu_fld*ydifatu3d(:,:,k)
      END WHERE
   ENDIF

!
!3.5 Update depth-integrated diffusive term
!------------------------------------------
!

   IF (vintflag) THEN
      WHERE (maskatu(:,:,k))
         udevint = udevint + ydifatu3d(:,:,k)
      END WHERE
   ENDIF

ENDDO k_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatu,maskatuv,mask2duv)
DEALLOCATE (array2duv1,array2duv2,denom,difflux,unorm,vnorm)
IF (.NOT.dXregY) DEALLOCATE (array2duv3)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Ydif_at_U_3d

!========================================================================

SUBROUTINE Ydif_at_U_2d(psiu,ydifatu2d,vintflag)
!************************************************************************
!
! *Ydif_at_U_2d* Diffusion term in Y-direction for the 2-D current at U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_U_2d
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE depths
USE diffusion
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(INOUT) :: vintflag
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo) :: psiu
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: ydifatu2d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psiu*        REAL    U-node quantity to be diffused                    [psiu]
!*ydifatu2d*   REAL    Diffusion term                                  [psiu/s]
!*vintflag*    LOGICAL Update vertically integrated diffusive term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, jj, jjloc, ju, jv, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatuv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2duv, difflux, unorm, vnorm

!
! Name      Type   Purpose
!------------------------------------------------------------------------------
!*difflux*  REAL  Diffusive flux at Y-nodes
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Ydif_at_U_2d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (maskatuv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('maskatuv',2,(/ncloc,nrloc+1/),kndlog)
ALLOCATE (array2duv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2duv',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (difflux(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('difflux',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (unorm(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('unorm',2,(/ncloc,nrloc+2/),kndrtype)
ALLOCATE (vnorm(0:ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('vnorm',2,(/ncloc+1,nrloc+1/),kndrtype)

!
!2. Initialise
!-------------
!
!2.1 Mask arrays
!---------------
!

maskatu = node2du(1:ncloc,1:nrloc).EQ.2
maskatuv = node2duv(1:ncloc,1:nrloc+1).GT.0

!
!2.2 Fluxes
!----------
!

WHERE (.NOT.maskatuv)
   difflux = 0.0
END WHERE

!
!2.3 Work space arrays
!---------------------
!

WHERE (maskatuv)
   array2duv = delxatuv(1:ncloc,1:nrloc+1)/delyatuv(1:ncloc,1:nrloc+1)
END WHERE
WHERE (node2du(1:ncloc,0:nrloc+1).GT.1)
   unorm = psiu(1:ncloc,0:nrloc+1)&
         & /(delxatu(1:ncloc,0:nrloc+1)*deptotatu(1:ncloc,:))
ELSEWHERE
   unorm = 0.0
END WHERE
WHERE (node2dv(0:ncloc,1:nrloc+1).GT.1)
   vnorm = vmvel_old(0:ncloc,1:nrloc+1)/delyatv(0:ncloc,1:nrloc+1)
ELSEWHERE
   vnorm = 0.0
END WHERE

!
!3. Explicit terms
!-----------------
!
!3.1 Fluxes
!----------
!

IF (dXregY) THEN
   WHERE (maskatuv)
      difflux = hdifmom_fac*hdifcoef2datuv(1:ncloc,:)*(array2duv*&
             & (unorm(:,1:nrloc+1)-unorm(:,0:nrloc))+&
             & (vnorm(1:ncloc,:)-vnorm(0:ncloc-1,:))/array2duv)
   END WHERE
ELSE
   WHERE (maskatuv)
      difflux = hdifmom_fac*hdifcoef2datuv(1:ncloc,:)*&
              & (delxatuv(1:ncloc,1:nrloc+1)**2)*&
               & (array2duv*(unorm(:,1:nrloc+1)-unorm(:,0:nrloc))+&
                        & (vnorm(1:ncloc,:)-vnorm(0:ncloc-1,:))/array2duv)
   END WHERE
ENDIF

!
!3.2 Open boundary fluxes
!------------------------
!

jjloc_320: DO jjloc=1,nobyloc_ext
   i = iobyloc(jjloc); j = jobyloc(jjloc)
   jj = indexoby(jjloc); ju = MERGE(j-1,j,soutoby(jj))
   IF (node2du(i,ju).GT.2.OR.node2dv(i-1,j).LE.1.OR.node2dv(i,j).LE.1) THEN
      IF (dXregY) THEN
         difflux(i,j) = hdifmom_fac*hdifcoef2datuv(i,j)*&
                      & (array2duv(i,j)*(unorm(i,j)-unorm(i,j-1))+&
                      & (vnorm(i,j)-vnorm(i-1,j))/array2duv(i,j))
      ELSE
         difflux(i,j) = hdifmom_fac*hdifcoef2datuv(i,j)*(delxatuv(i,j)**2)*&
                        (array2duv(i,j)*(unorm(i,j)-unorm(i,j-1))+&
                      & (vnorm(i,j)-vnorm(i-1,j))/array2duv(i,j))
      ENDIF
   ELSE
      jv = MERGE(j+1,j-1,soutoby(jj))
      IF (jv.GE.1.AND.jv.LE.nrloc+1.AND.node2du(i,ju).EQ.2) THEN
         difflux(i,j) = difflux(i,jv)
      ENDIF
   ENDIF
ENDDO jjloc_320

!
!3.3 Add flux divergence to explicit term
!----------------------------------------
!

IF (dXregY) THEN
   WHERE (maskatu)
      ydifatu2d = (difflux(:,2:nrloc+1)-difflux(:,1:nrloc))&
                & /delyatu(1:ncloc,1:nrloc)
   END WHERE
ELSE
   WHERE (maskatu)
      ydifatu2d = (difflux(:,2:nrloc+1)-difflux(:,1:nrloc))&
                & /(delyatu(1:ncloc,1:nrloc)*delxatu(1:ncloc,1:nrloc)**2)
   END WHERE
ENDIF

!
!3.4 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN
   WHERE (maskatu)
      ydifatu2d = alphatu_fld*ydifatu2d
   END WHERE
ENDIF

!
!3.5 Update depth-integrated diffusive term
!------------------------------------------
!

IF (vintflag) THEN
   WHERE (maskatu)
      udevint = udevint - ydifatu2d
   END WHERE
   vintflag = .FALSE.
ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatu,maskatuv)
DEALLOCATE (array2duv,difflux,unorm,vnorm)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Ydif_at_U_2d

!========================================================================

SUBROUTINE Zdif_at_U(psiu,tridcfu,vdifcoefatuw,ibcsur,ibcbot,nbcs,nbcb,&
                   & bcsur,bcbot)
!************************************************************************
!
! *Zdif_at_U* Vertical diffusion term for the 3-D current at U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - transport_at_U
!
! Module calls - error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ibcsur, ibcbot, nbcs, nbcb
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiu
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4) :: tridcfu
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1) :: vdifcoefatuw
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcs) :: bcsur
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcb) :: bcbot

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psiu*        REAL    U-node quantity to be diffused                  [psiu]
!*tridcfu*     REAL    Tridiagonal matrix for implicit vertical solution
!*vdifcoefatuw*REAL    Diffusion coefficient at UW-nodes              [m^2/s]
!*ibcsur*      INTEGER Type of surface boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!*ibcbot*      INTEGER Type of bottom boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!*nbcs*        INTEGER Last dimension of array bcsur
!*nbcb*        INTEGER Last dimension of array bcbot
!*bcsur*       REAL    Data for surface boundary condition
!              (:,:,1) => prescribed surface flux or surface value
!              (:,:,2) => transfer velocity                              [m/s]
!*bcbot*       REAL    Data for bottom boundary condition
!              (:,:,1) => prescribed surface flux or surface value
!              (:,:,2) => transfer velocity                              [m/s]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
REAL :: theta_vdif1, xexp, ximp
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu, maskatuw 
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3du, array3duw, difflux

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*difflux*    REAL    Diffusive flux (times factor) at UW-nodes      [m^2/psiu]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Zdif_at_U'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array3du(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('array3du',3,(/ncloc,nrloc,nz+1/),kndrtype)
ALLOCATE (array3duw(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('array3uw',3,(/ncloc,nrloc,nz-1/),kndrtype)
ALLOCATE (maskatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (maskatuw(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('maskatuw',3,(/ncloc,nrloc,nz-1/),kndlog)
IF (iopt_vdif_impl.NE.2) THEN
   ALLOCATE (difflux(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('difflux',3,(/ncloc,nrloc,nz+1/),kndrtype)
ENDIF

!
!2. Initialise
!-------------
!
!2.1 Mask arrays
!---------------
!

maskatu = nodeatu(1:ncloc,1:nrloc,:).EQ.2
maskatuw = nodeatuw(1:ncloc,1:nrloc,2:nz).EQ.2

!
!2.2 Fluxes
!----------
!

IF (iopt_vdif_impl.NE.2) difflux = 0.0

!
!2.3 Work space arrays
!---------------------
!

k_230: DO k=2,nz
   WHERE (maskatuw(:,:,k))
      array3duw(:,:,k) = vdifcoefatuw(:,:,k)/delzatuw(1:ncloc,:,k)
   ELSEWHERE
      array3duw(:,:,k) = 0.0
   END WHERE
ENDDO k_230

!---time factors
xexp = delt3d*(1.0-theta_vdif)
ximp = delt3d*theta_vdif
theta_vdif1 = 1.0-theta_vdif

!
!3. Explicit fluxes
!------------------
!

IF (iopt_vdif_impl.NE.2) THEN
   WHERE (maskatuw)
      difflux(:,:,2:nz) = xexp*array3duw*&
                     & (psiu(1:ncloc,1:nrloc,2:nz)-psiu(1:ncloc,1:nrloc,1:nz-1))
   END WHERE
ENDIF

!
!4. Explicit terms
!-----------------
!

IF (iopt_vdif_impl.NE.2) THEN
   WHERE (maskatu)
      tridcfu(:,:,:,4) = tridcfu(:,:,:,4) + &
                       & (difflux(:,:,2:nz+1)-difflux(:,:,1:nz))&
                       & /delzatu(1:ncloc,:,:)
   END WHERE
ENDIF

!
!5. Implicit terms
!-----------------
!

IF (iopt_vdif_impl.NE.0) THEN

!  ---lower flux
   WHERE (maskatu(:,:,2:nz))
      array3du(:,:,2:nz) = ximp*array3duw/delzatu(1:ncloc,:,2:nz)
      tridcfu(:,:,2:nz,1) = tridcfu(:,:,2:nz,1) - array3du(:,:,2:nz)
      tridcfu(:,:,2:nz,2) = tridcfu(:,:,2:nz,2) + array3du(:,:,2:nz)
   END WHERE

!  ---upper flux
   WHERE (maskatu(:,:,1:nz-1))
      array3du(:,:,1:nz-1) = ximp*array3duw/delzatu(1:ncloc,:,1:nz-1)
      tridcfu(:,:,1:nz-1,2) = tridcfu(:,:,1:nz-1,2) + array3du(:,:,1:nz-1)
      tridcfu(:,:,1:nz-1,3) = tridcfu(:,:,1:nz-1,3) - array3du(:,:,1:nz-1)
   END WHERE

ENDIF

!
!6. Boundary conditions
!----------------------
!
!6.1 Surface
!-----------
!
!---Neumann (prescribed flux)
IF (ibcsur.EQ.1) THEN
   WHERE (maskatu(:,:,nz))
      tridcfu(:,:,nz,4) = tridcfu(:,:,nz,4) + delt3d*bcsur(:,:,1)&
                       & /delzatu(1:ncloc,:,nz)
   END WHERE

!---Neumann (using transfer velocity)
ELSEIF (ibcsur.EQ.2) THEN
   WHERE (maskatu(:,:,nz))
      tridcfu(:,:,nz,2) = tridcfu(:,:,nz,2) + ximp*bcsur(:,:,2)&
                        & /delzatu(1:ncloc,:,nz)
      tridcfu(:,:,nz,4) = tridcfu(:,:,nz,4) - delt3d*bcsur(:,:,2)*&
                        & (theta_vdif1*psiu(1:ncloc,1:nrloc,nz)-bcsur(:,:,1))&
                        & /delzatu(1:ncloc,:,nz)
   END WHERE
ENDIF

!
!6.2 Bottom
!----------
!
!---Neumann (prescribed flux)
IF (ibcbot.EQ.1) THEN
   WHERE (maskatu(:,:,1))
      tridcfu(:,:,1,4) = tridcfu(:,:,1,4) - delt3d*bcbot(:,:,1)&
                       & /delzatu(1:ncloc,:,1)
   END WHERE

!---Neumann (using transfer velocity)
ELSEIF (ibcbot.EQ.2) THEN
   WHERE (maskatu(:,:,1))
      tridcfu(:,:,1,2) = tridcfu(:,:,1,2) + ximp*bcbot(:,:,2)&
                       & /delzatu(1:ncloc,:,1)
      tridcfu(:,:,1,4) = tridcfu(:,:,1,4) - delt3d*bcbot(:,:,2)*&
                       & (theta_vdif1*psiu(1:ncloc,1:nrloc,1)-bcbot(:,:,1))&
                       & /delzatu(1:ncloc,:,1)
   END WHERE
ENDIF

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE (array3du,array3duw,maskatu,maskatuw)
IF (iopt_vdif_impl.NE.2) DEALLOCATE (difflux)

CALL log_timer_out(npcc,itm_vdif)


RETURN

END SUBROUTINE Zdif_at_U

!========================================================================

SUBROUTINE Xdif_at_V_3d(psiv,xdifatv3d,vintflag)
!************************************************************************
!
! *Xdif_at_V_3d* Diffusion term in X-direction for the 3-D current at V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_V_3d
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE diffusion  
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: vintflag
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiv
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,nz) :: xdifatv3d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psiv*        REAL    V-node quantity to be diffused                   [psiv]
!*xdifatv3d*   REAL    Diffusion term (times delz)                  [m*psiv/s]
!*vintflag*    LOGICAL Update vertically integrated diffusive term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ii, iiloc, iu, iv, j, k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2duv
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatv, maskatuv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2duv1, array2duv2, array2duv3, &
                                         & denom, difflux, unorm, vnorm

!
! Name      Type   Purpose
!------------------------------------------------------------------------------
!*difflux*  REAL  Diffusive flux at UV-nodes
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xdif_at_V_3d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (maskatuv(ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatuv',3,(/ncloc+1,nrloc,nz/),kndlog)
ALLOCATE (mask2duv(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('mask2duv',2,(/ncloc+1,nrloc/),kndlog)
ALLOCATE (array2duv1(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2duv1',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (array2duv2(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2duv2',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (difflux(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('difflux',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (unorm(ncloc+1,0:nrloc),STAT=errstat)
CALL error_alloc('unorm',2,(/ncloc+1,nrloc+1/),kndrtype)
ALLOCATE (vnorm(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('vnorm',2,(/ncloc+2,nrloc/),kndrtype)
IF (.NOT.dYregX) THEN
   ALLOCATE (array2duv3(ncloc+1,nrloc),STAT=errstat)
   CALL error_alloc('array2duv3',2,(/ncloc+1,nrloc/),kndrtype)
ENDIF

!
!2. Initialise
!-------------
!
!2.1 Mask arrays
!---------------
!

maskatv = nodeatv(1:ncloc,1:nrloc,:).EQ.2
maskatuv = nodeatuv(1:ncloc+1,1:nrloc,:).GT.0
mask2duv = node2duv(1:ncloc+1,1:nrloc).GT.0

!
!2.2 Work space arrays
!---------------------
!

IF (dYregX) THEN
   denom = delxatv(1:ncloc,1:nrloc)
ELSE
   denom = delxatv(1:ncloc,1:nrloc)*delyatv(1:ncloc,1:nrloc)**2
   WHERE (mask2duv)
      array2duv3 = delyatuv(1:ncloc+1,1:nrloc)**2
   END WHERE
ENDIF
WHERE (mask2duv)
   array2duv1 = delyatuv(1:ncloc+1,1:nrloc)/delxatuv(1:ncloc+1,1:nrloc)
END WHERE

!
!3. Diffusion term
!-----------------
!

k_300: DO k=1,nz

!
!3.1 Fluxes at interior nodes
!----------------------------
!

   WHERE (nodeatu(1:ncloc+1,0:nrloc,k).GT.1)
      unorm = uvel_old(1:ncloc+1,0:nrloc,k)/delxatu(1:ncloc+1,0:nrloc)
   ELSEWHERE
      unorm = 0.0
   END WHERE
   WHERE (nodeatv(0:ncloc+1,1:nrloc,k).GT.1) 
      vnorm = psiv(0:ncloc+1,1:nrloc,k)/delyatv(0:ncloc+1,1:nrloc)
   ELSEWHERE
      vnorm = 0.0
   END WHERE
   IF (dYregX) THEN
      WHERE (maskatuv(:,:,k))
         array2duv2 = hdifmom_fac*hdifcoef3datuv(:,1:nrloc,k)*&
                    & delzatuv(:,1:nrloc,k)
      END WHERE
   ELSE
      WHERE (maskatuv(:,:,k))
         array2duv2 = hdifmom_fac*hdifcoef3datuv(:,1:nrloc,k)*&
                    & delzatuv(:,1:nrloc,k)*array2duv3
      END WHERE
   ENDIF
   WHERE (maskatuv(:,:,k))
      difflux = array2duv2*(array2duv1*(vnorm(1:ncloc+1,:)-vnorm(0:ncloc,:))+&
                         & (unorm(:,1:nrloc)-unorm(:,0:nrloc-1))/array2duv1)
   ELSEWHERE
      difflux = 0.0
   END WHERE

!
!3.2 Open boundary fluxes
!------------------------
!

   iiloc_320: DO iiloc=1,nobxloc_ext
      i = iobxloc(iiloc); j = jobxloc(iiloc)
      ii = indexobx(iiloc); iv = MERGE(i-1,i,westobx(ii))
      IF (nodeatv(iv,j,k).GT.2.OR.nodeatu(i,j-1,k).LE.1.OR.&
        & nodeatu(i,j,k).LE.1) THEN
         difflux(i,j) = array2duv2(i,j)*(array2duv1(i,j)*&
                      & (vnorm(i,j)-vnorm(i-1,j))+&
                      & (unorm(i,j)-unorm(i,j-1))/array2duv1(i,j))
      ELSE
         iu = MERGE(i+1,i-1,westobx(ii))
         IF (iu.GE.1.AND.iu.LE.ncloc+1.AND.nodeatv(iv,j,k).EQ.2) THEN
            difflux(i,j) = difflux(iu,j)
         ENDIF
      ENDIF
   ENDDO iiloc_320

!
!3.3 Add flux divergence to explicit term
!----------------------------------------
!

   WHERE (maskatv(:,:,k))
      xdifatv3d(:,:,k) = (difflux(2:ncloc+1,:)-difflux(1:ncloc,:))/denom
   END WHERE

!
!3.4 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

   IF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatv(:,:,k))
         xdifatv3d(:,:,k) = alphatv_fld*xdifatv3d(:,:,k)
      END WHERE
   ENDIF

!
!3.5 Update depth-integrated diffusive term
!------------------------------------------
!

   IF (vintflag) THEN
      WHERE (maskatv(:,:,k))
         vdevint = vdevint + xdifatv3d(:,:,k)
      END WHERE
   ENDIF

ENDDO k_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatv,maskatuv,mask2duv)
DEALLOCATE (array2duv1,array2duv2,denom,difflux,unorm,vnorm)
IF (.NOT.dYregX) DEALLOCATE (array2duv3)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Xdif_at_V_3d

!========================================================================

SUBROUTINE Xdif_at_V_2d(psiv,xdifatv2d,vintflag)
!************************************************************************
!
! *Xdif_at_V_2d* Diffusion term in X-direction for the 2-D current at V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_V_2d
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE depths
USE diffusion
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(INOUT) :: vintflag
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo) :: psiv
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: xdifatv2d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psiv*        REAL    V-node quantity to be diffused                    [psiv]
!*xdifatv2d*   REAL    Diffusion term                                  [psiv/s]
!*vintflag*    LOGICAL Update vertically integrated diffusive term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ii, iiloc, iu, iv, j, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatv, maskatuv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2duv, difflux, unorm, vnorm

!
! Name      Type   Purpose
!------------------------------------------------------------------------------
!*difflux*  REAL  Diffusive flux at X-nodes
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xdif_at_V_2d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (maskatuv(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('maskatuv',2,(/ncloc+1,nrloc/),kndlog)
ALLOCATE (array2duv(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2duv',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (difflux(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('difflux',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (unorm(ncloc+1,0:nrloc),STAT=errstat)
CALL error_alloc('unorm',2,(/ncloc+1,nrloc+1/),kndrtype)
ALLOCATE (vnorm(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('vnorm',2,(/ncloc+2,nrloc/),kndrtype)

!
!2. Initialise
!-------------
!
!2.1 Mask arrays
!---------------
!

maskatv = node2dv(1:ncloc,1:nrloc).EQ.2
maskatuv = node2duv(1:ncloc+1,1:nrloc).GT.0

!
!2.2 Fluxes
!----------
!

WHERE (.NOT.maskatuv)
   difflux = 0.0
END WHERE

!
!2.3 Work space arrays
!---------------------
!

WHERE (maskatuv)
   array2duv = delyatuv(1:ncloc+1,1:nrloc)/delxatuv(1:ncloc+1,1:nrloc)
END WHERE
WHERE (node2du(1:ncloc+1,0:nrloc).GT.1)
   unorm = umvel_old(1:ncloc+1,0:nrloc)/delxatu(1:ncloc+1,0:nrloc)
ELSEWHERE
   unorm = 0.0
END WHERE
WHERE (node2dv(0:ncloc+1,1:nrloc).GT.1)
   vnorm = psiv(0:ncloc+1,1:nrloc)&
         & /(delyatv(0:ncloc+1,1:nrloc)*deptotatv(:,1:nrloc))
ELSEWHERE
   vnorm = 0.0
END WHERE

!
!3. Explicit terms
!-----------------
!
!3.1 Fluxes
!----------
!

IF (dYregX) THEN
   WHERE (maskatuv)
      difflux = hdifmom_fac*hdifcoef2datuv(:,1:nrloc)*(array2duv*&
             & (vnorm(1:ncloc+1,:)-vnorm(0:ncloc,:))+&
             & (unorm(:,1:nrloc)-unorm(:,0:nrloc-1))/array2duv)
   END WHERE
ELSE
   WHERE (maskatuv)
      difflux = hdifmom_fac*hdifcoef2datuv(:,1:nrloc)*&
              & (delyatuv(1:ncloc+1,1:nrloc)**2)*&
              & (array2duv*(vnorm(1:ncloc+1,:)-vnorm(0:ncloc,:))+&
                        & (unorm(:,1:nrloc)-unorm(:,0:nrloc-1))/array2duv)
   END WHERE
ENDIF

!
!3.2 Open boundary fluxes
!------------------------
!

iiloc_320: DO iiloc=1,nobxloc_ext
   i = iobxloc(iiloc); j = jobxloc(iiloc)
   ii = indexobx(iiloc); iv = MERGE(i-1,i,westobx(ii))
   IF (node2dv(iv,j).GT.2.OR.node2du(i,j-1).LE.1.OR.node2du(i,j).LE.1) THEN
      IF (dYregX) THEN
         difflux(i,j) = hdifmom_fac*hdifcoef2datuv(i,j)*&
                      & (array2duv(i,j)*(vnorm(i,j)-vnorm(i-1,j))+&
                      & (unorm(i,j)-unorm(i,j-1))/array2duv(i,j))
      ELSE
         difflux(i,j) = hdifmom_fac*hdifcoef2datuv(i,j)*(delyatuv(i,j)**2)*&
                        (array2duv(i,j)*(vnorm(i,j)-vnorm(i-1,j))+&
                      & (unorm(i,j)-unorm(i,j-1))/array2duv(i,j))
      ENDIF
   ELSE
      iu = MERGE(i+1,i-1,westobx(ii))
      IF (iu.GE.1.AND.iu.LE.ncloc+1.AND.node2dv(iv,j).EQ.2) THEN
         difflux(i,j) = difflux(iu,j)
      ENDIF
   ENDIF
ENDDO iiloc_320

!
!3.3 Add flux divergence to explicit term
!----------------------------------------
!

IF (dYregX) THEN
   WHERE (maskatv)
      xdifatv2d = (difflux(2:ncloc+1,:)-difflux(1:ncloc,:))&
                & /delxatv(1:ncloc,1:nrloc)
   END WHERE
ELSE
   WHERE (maskatv)
      xdifatv2d = (difflux(2:ncloc+1,:)-difflux(1:ncloc,:))&
                & /(delxatv(1:ncloc,1:nrloc)*delyatv(1:ncloc,1:nrloc)**2)
   END WHERE
ENDIF

!
!3.4 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN
   WHERE (maskatv)
      xdifatv2d = alphatv_fld*xdifatv2d
   END WHERE
ENDIF

!
!3.5 Update depth-integrated diffusive term
!------------------------------------------
!

IF (vintflag) THEN
   WHERE (maskatv)
      vdevint = vdevint - xdifatv2d
   END WHERE
   vintflag = .FALSE.
ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatv,maskatuv)
DEALLOCATE (array2duv,difflux,unorm,vnorm)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Xdif_at_V_2d

!========================================================================

SUBROUTINE Ydif_at_V_3d(psiv,ydifatv3d,vintflag)
!************************************************************************
!
! *Ydif_at_V_3d* Diffusion term in Y-direction for the 3-D current at V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_V_3d
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE diffusion  
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: vintflag
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiv
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,nz) :: ydifatv3d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psiv*        REAL    V-node quantity to be diffused                    [psiv]
!*ydifatv3d*   REAL    Diffusion term (times delz)                   [m*psiu/s]
!*vintflag*    LOGICAL Update vertically integrated diffusive term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc1, array2dc2, denom, &
                                         & difflux, unorm, vnorm

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*difflux*  REAL    Diffusive flux at C-nodes
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Ydif_at_V_3d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatc(ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc,nrloc+1/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (array2dc1(ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('array2dc1',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (array2dc2(ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('array2dc2',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (difflux(ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('difflux',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (unorm(ncloc+1,0:nrloc),STAT=errstat)
CALL error_alloc('unorm',2,(/ncloc+1,nrloc+1/),kndrtype)
ALLOCATE (vnorm(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('vnorm',2,(/ncloc,nrloc+2/),kndrtype)

!
!2. Initialise
!-------------
!
!2.1 Mask arrays
!---------------
!

maskatv = nodeatv(1:ncloc,1:nrloc,:).EQ.2
maskatc = nodeatc(1:ncloc,0:nrloc).GT.0

!
!2.2 Fluxes
!----------
!

WHERE (.NOT.maskatc) difflux = 0.0

!
!2.3 Work space arrays
!---------------------
!

IF (dXregY) THEN
   denom = delyatv(1:ncloc,1:nrloc)
ELSE
   denom = delyatv(1:ncloc,1:nrloc)*delxatv(1:ncloc,1:nrloc)**2
ENDIF
WHERE (maskatc)
   array2dc1 = delxatc(1:ncloc,0:nrloc)/delyatc(1:ncloc,0:nrloc)
END WHERE

!
!3. Explicit terms
!-----------------
!

k_300: DO k=1,nz

!
!3.1 Fluxes
!----------
!
!  ---work space arrays
   WHERE (nodeatu(1:ncloc+1,0:nrloc,k).GT.1)
      unorm = uvel_old(1:ncloc+1,0:nrloc,k)/delyatu(1:ncloc+1,0:nrloc)
      ELSEWHERE
         unorm = 0.0
   END WHERE
   WHERE (nodeatv(1:ncloc,0:nrloc+1,k).GT.1)
      vnorm = psiv(1:ncloc,0:nrloc+1,k)/delxatv(1:ncloc,0:nrloc+1)
   ELSEWHERE
      vnorm = 0.0
   END WHERE
   IF (dXregY) THEN
      WHERE (maskatc)
         array2dc2 = hdifmom_fac*hdifcoef3datc(1:ncloc,:,k)*&
                   & delzatc(1:ncloc,0:nrloc,k)
      END WHERE
   ELSE
      WHERE (maskatc)
         array2dc2 = hdifmom_fac*hdifcoef3datc(1:ncloc,:,k)*&
                   & delzatc(1:ncloc,0:nrloc,k)*delxatc(1:ncloc,0:nrloc)**2
      END WHERE
   ENDIF
   WHERE (maskatc)
      difflux = array2dc2*(array2dc1*(vnorm(:,1:nrloc+1)-vnorm(:,0:nrloc))-&
                        & (unorm(2:ncloc+1,:)-unorm(1:ncloc,:))/array2dc1)
   END WHERE

!
!3.2 Add flux divergence to explicit term
!----------------------------------------
!

   WHERE (maskatv(:,:,k))
      ydifatv3d(:,:,k) = (difflux(:,1:nrloc)-difflux(:,0:nrloc-1))/denom
   END WHERE

!
!3.3 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

   IF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatv(:,:,k))
         ydifatv3d(:,:,k) = alphatv_fld*ydifatv3d(:,:,k)
      END WHERE
   ENDIF

!
!3.4 Update depth-integrated diffusive term
!------------------------------------------
!

   IF (vintflag) THEN
      WHERE (maskatv(:,:,k))
         vdevint = vdevint + ydifatv3d(:,:,k)
      END WHERE
   ENDIF

ENDDO k_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatc,maskatv)
DEALLOCATE (array2dc1,array2dc2,denom,difflux,unorm,vnorm)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Ydif_at_V_3d

!========================================================================

SUBROUTINE Ydif_at_V_2d(psiv,ydifatv2d,vintflag)
!************************************************************************
!
! *Ydif_at_V_2d* Diffusion term in Y-direction for the 2-D current at V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_V_2d
!
! Module calls - error_alloc
!
!************************************************************************
!
USE currents
USE depths
USE diffusion
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(INOUT) :: vintflag
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo) :: psiv
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: ydifatv2d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psiv*        REAL    V-node quantity to be diffused                    [psiv]
!*ydifatv2d*   REAL    Diffusion term                                  [psiv/s]
!*vintflag*    LOGICAL Update vertically integrated diffusive term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, difflux, unorm, vnorm

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*difflux*   REAL    Diffusive flux at C-nodes
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Ydif_at_V_2d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatc(ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc,nrloc+1/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (array2dc(ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (difflux(ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('difflux',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (unorm(ncloc+1,0:nrloc),STAT=errstat)
CALL error_alloc('unorm',2,(/ncloc+1,nrloc+1/),kndrtype)
ALLOCATE (vnorm(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('vnorm',2,(/ncloc,nrloc+2/),kndrtype)

!
!2. Initialise
!-------------
!
!2.1 Mask arrays
!---------------
!

maskatv = node2dv(1:ncloc,1:nrloc).EQ.2
maskatc = nodeatc(1:ncloc,0:nrloc).GT.0

!
!2.2 Fluxes
!----------
!

WHERE (.NOT.maskatc)
   difflux = 0.0
END WHERE

!
!2.3 Work space arrays
!---------------------
!

WHERE (node2du(1:ncloc+1,0:nrloc).GT.1)
   unorm = umvel_old(1:ncloc+1,0:nrloc)/delyatu(1:ncloc+1,0:nrloc)
ELSEWHERE
   unorm = 0.0
END WHERE
WHERE (node2dv(1:ncloc,0:nrloc+1).GT.1)
   vnorm = psiv(1:ncloc,0:nrloc+1)&
         & /(delxatv(1:ncloc,0:nrloc+1)*deptotatv(1:ncloc,0:nrloc+1))
ELSEWHERE
   vnorm = 0.0
END WHERE
WHERE (maskatc) 
   array2dc = delxatc(1:ncloc,0:nrloc)/delyatc(1:ncloc,0:nrloc)
END WHERE

!
!3. Explicit terms
!-----------------
!
!3.1 Fluxes
!----------
!
!---regular grid
IF (dXregY) THEN
   WHERE (maskatc)
      difflux = hdifmom_fac*hdifcoef2datc(1:ncloc,:)*&
              & (array2dc*(vnorm(:,1:nrloc+1)-vnorm(:,0:nrloc))-&
              & (unorm(2:ncloc+1,:)-unorm(1:ncloc,:))/array2dc)
   END WHERE

!---irregular grid
ELSE
   WHERE (maskatc)
      difflux = hdifmom_fac*hdifcoef2datc(1:ncloc,:)*&
              & (delxatc(1:ncloc,0:nrloc)**2)*&
              & (array2dc*(vnorm(:,1:nrloc+1)-vnorm(:,0:nrloc))-&
              & (unorm(2:ncloc+1,:)-unorm(1:ncloc,:))/array2dc)
   END WHERE
ENDIF

!
!3.2 Add flux divergence to explicit term
!----------------------------------------
!

IF (dXregY) THEN
   WHERE (maskatv)
      ydifatv2d = (difflux(:,1:nrloc)-difflux(:,0:nrloc-1))&
                & /delyatv(1:ncloc,1:nrloc)
   END WHERE
ELSE
   WHERE (maskatv)
      ydifatv2d = (difflux(:,1:nrloc)-difflux(:,0:nrloc-1))&
                & /(delyatv(1:ncloc,1:nrloc)*delxatv(1:ncloc,1:nrloc)**2)
   END WHERE
ENDIF

!
!3.3 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN
   WHERE (maskatv)
      ydifatv2d = alphatv_fld*ydifatv2d
   END WHERE
ENDIF

!
!3.4 Update depth-integrated diffusive term
!------------------------------------------
!

IF (vintflag) THEN
   WHERE (maskatv)
      vdevint = vdevint - ydifatv2d
   END WHERE
   vintflag = .FALSE.
ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatc,maskatv)
DEALLOCATE (array2dc,difflux,unorm,vnorm)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Ydif_at_V_2d

!========================================================================

SUBROUTINE Zdif_at_V(psiv,tridcfv,vdifcoefatvw,ibcsur,ibcbot,nbcs,nbcb,&
                   & bcsur,bcbot)
!************************************************************************
!
! *Zdif_at_V* Vertical diffusion term for the 3-D current at V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_V
!
! Module calls - error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ibcsur, ibcbot, nbcs, nbcb
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiv
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4) :: tridcfv
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1) :: vdifcoefatvw
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcs) :: bcsur
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcb) :: bcbot

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psiv*        REAL    V-node quantity to be diffused                   [psiv]
!*tridcfv*     REAL    Tridiagonal matrix for implicit vertical solution
!*vdifcoefatvw*REAL    Diffusion coefficient at VW-nodes               [m^2/s]
!*ibcsur*      INTEGER Type of surface boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!*ibcbot*      INTEGER Type of bottom boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!*nbcs*        INTEGER Last dimension of array bcsur
!*nbcb*        INTEGER Last dimension of array bcbot
!*bcsur*       REAL    Data for surface boundary condition
!              (:,:,1) => prescribed surface flux or surface value
!              (:,:,2) => transfer velocity                              [m/s]
!*bcbot*       REAL    Data for bottom boundary condition
!              (:,:,1) => prescribed surface flux or surface value
!              (:,:,2) => transfer velocity                              [m/s]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
REAL :: theta_vdif1, xexp, ximp
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatv, maskatvw 
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dv, array3dvw, difflux

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*difflux*    REAL    Diffusive flux (times factor) at VW-nodes      [m^2/psiv]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Zdif_at_V'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array3dv(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('array3dvw',3,(/ncloc,nrloc,nz+1/),kndrtype)
ALLOCATE (array3dvw(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('array3dvw',3,(/ncloc,nrloc,nz-1/),kndrtype)
ALLOCATE (maskatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (maskatvw(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('maskatvw',3,(/ncloc,nrloc,nz-1/),kndlog)
IF (iopt_vdif_impl.NE.2) THEN
   ALLOCATE (difflux(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('difflux',3,(/ncloc,nrloc,nz+1/),kndrtype)
ENDIF

!
!2. Initialise
!-------------
!
!2.1 Mask arrays
!---------------
!

maskatv = nodeatv(1:ncloc,1:nrloc,:).EQ.2
maskatvw = nodeatvw(1:ncloc,1:nrloc,2:nz).EQ.2

!
!2.2 Fluxes
!----------
!

IF (iopt_vdif_impl.NE.2) THEN
   difflux(:,:,nz+1) = 0.0
   WHERE (.NOT.maskatvw) difflux = 0.0
ENDIF

!
!2.3 Work space arrays
!---------------------
!

k_230: DO k=2,nz
   WHERE (maskatvw(:,:,k))
      array3dvw(:,:,k) = vdifcoefatvw(:,:,k)/delzatvw(:,1:nrloc,k)
   ELSEWHERE
      array3dvw(:,:,k) = 0.0
   END WHERE
ENDDO k_230

!---time factors
xexp = delt3d*(1.0-theta_vdif)
ximp = delt3d*theta_vdif
theta_vdif1 = 1.0-theta_vdif

!
!3. Explicit fluxes
!------------------
!
IF (iopt_vdif_impl.NE.2) THEN
   WHERE (maskatvw)
      difflux(:,:,2:nz) = xexp*array3dvw*&
                     & (psiv(1:ncloc,1:nrloc,2:nz)-psiv(1:ncloc,1:nrloc,1:nz-1))
      END WHERE
ENDIF

!
!4. Explicit terms
!-----------------
!

IF (iopt_vdif_impl.NE.2) THEN
   WHERE (maskatv)
      tridcfv(:,:,:,4) = tridcfv(:,:,:,4) + &
                       & (difflux(:,:,2:nz+1)-difflux(:,:,1:nz))&
                       & /delzatv(:,1:nrloc,:)
   END WHERE
ENDIF

!
!5. Implicit terms
!-----------------
!

IF (iopt_vdif_impl.NE.0) THEN

!  ---lower flux
   WHERE (maskatv(:,:,2:nz))
      array3dv(:,:,2:nz) = ximp*array3dvw/delzatv(:,1:nrloc,2:nz)
      tridcfv(:,:,2:nz,1) = tridcfv(:,:,2:nz,1) - array3dv(:,:,2:nz)
      tridcfv(:,:,2:nz,2) = tridcfv(:,:,2:nz,2) + array3dv(:,:,2:nz)
   END WHERE

!  ---upper flux
   WHERE (maskatv(:,:,1:nz-1))
      array3dv(:,:,1:nz-1) = ximp*array3dvw/delzatv(:,1:nrloc,1:nz-1)
      tridcfv(:,:,1:nz-1,2) = tridcfv(:,:,1:nz-1,2) + array3dv(:,:,1:nz-1)
      tridcfv(:,:,1:nz-1,3) = tridcfv(:,:,1:nz-1,3) - array3dv(:,:,1:nz-1)
   END WHERE
ENDIF

!
!6. Boundary conditions
!----------------------
!
!6.1 Surface
!-----------
!
!---Neumann (prescribed flux)
IF (ibcsur.EQ.1) THEN
   WHERE (maskatv(:,:,nz))
      tridcfv(:,:,nz,4) = tridcfv(:,:,nz,4) + delt3d*bcsur(:,:,1)&
                       & /delzatv(:,1:nrloc,nz)
   END WHERE

!---Neumann (using transfer velocity)
ELSEIF (ibcsur.EQ.2) THEN
   WHERE (maskatv(:,:,nz))
      tridcfv(:,:,nz,2) = tridcfv(:,:,nz,2) + ximp*bcsur(:,:,2)&
                        & /delzatv(:,1:nrloc,nz)
      tridcfv(:,:,nz,4) = tridcfv(:,:,nz,4) - delt3d*bcsur(:,:,2)*&
                        & (theta_vdif1*psiv(1:ncloc,1:nrloc,nz)-bcsur(:,:,1))&
                        & /delzatv(:,1:nrloc,nz)
   END WHERE
ENDIF

!
!6.2 Bottom
!----------
!
!---Neumann (prescribed flux)
IF (ibcbot.EQ.1) THEN
   WHERE (maskatv(:,:,1))
      tridcfv(:,:,1,4) = tridcfv(:,:,1,4) - delt3d*bcbot(:,:,1)&
                       & /delzatv(:,1:nrloc,1)
   END WHERE

!---Neumann (using transfer velocity)
ELSEIF (ibcbot.EQ.2) THEN
   WHERE (maskatv(:,:,1))
      tridcfv(:,:,1,2) = tridcfv(:,:,1,2) + ximp*bcbot(:,:,2)&
                       & /delzatv(:,1:nrloc,1)
      tridcfv(:,:,1,4) = tridcfv(:,:,1,4) - delt3d*bcbot(:,:,2)*&
                       & (theta_vdif1*psiv(1:ncloc,1:nrloc,1)-bcbot(:,:,1))&
                       & /delzatv(:,1:nrloc,1)
   END WHERE
ENDIF

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE (array3dv,array3dvw,maskatv,maskatvw)
IF (iopt_vdif_impl.NE.2) DEALLOCATE (difflux)

CALL log_timer_out(npcc,itm_vdif)


RETURN

END SUBROUTINE Zdif_at_V

!========================================================================

SUBROUTINE Xdif_at_W(psiw,tridcfw,hdifcoefatuw,klo,kup)
!************************************************************************
!
! *Xdif_at_W* Diffusion term in X-direction for a quantity at W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_W
!
! Module calls - error_alloc
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1)::&
                          & psiw
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz+1,4) :: tridcfw
REAL, INTENT(IN), DIMENSION(ncloc+1,nrloc,2:nz) :: hdifcoefatuw

!
! Name         Type    Purpose
!-----------------------------------------------------------------------------
!*psiw*        REAL    W-node quantity to be diffused                   [psiw]
!*tridcfw*     REAL    Tridiagonal matrix for implicit (vertical) solution
!*hdifcoefatuw*REAL    Diffusion coefficient at UW-nodes               [m^2/s]
!*klo*         INTEGER Lower vertical array bound
!*kup*         INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2duw
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatuw
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, array2du
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3duw, denom, difflux

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*difflux*    REAL    Diffusive flux (times factor) at UW-nodes       
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xdif_at_W'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2du(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (array3duw(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('array3duw',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (difflux(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('difflux',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (maskatuw(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('maskatuw',3,(/ncloc+1,nrloc,kup-klo+1/),kndlog)
ALLOCATE (mask2duw(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('mask2duw',3,(/ncloc+1,nrloc/),kndlog)

!
!2. Initialise
!-------------
!
!2.1 Mask array
!--------------
!

maskatuw = nodeatuw(1:ncloc+1,1:nrloc,klo:kup).EQ.2
mask2duw = node2du(1:ncloc+1,1:nrloc).EQ.2

!
!2.2 Fluxes
!----------
!

WHERE (.NOT.maskatuw)
   difflux = 0.0
END WHERE

!
!2.3 Work space arrays
!---------------------
!

IF (dYregX) THEN
   k_231: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delxatc(1:ncloc,1:nrloc)*delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_231
ELSE
   k_232: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_232
ENDIF

IF (dYregX) THEN
   WHERE (mask2duw)
      array2du = delt3d/delxatu(1:ncloc+1,1:nrloc)
   END WHERE
ELSE
   WHERE (mask2duw)
      array2du = delt3d*delyatu(1:ncloc+1,1:nrloc)/delxatu(1:ncloc+1,1:nrloc)
   END WHERE
ENDIF
k_233: DO k=klo,kup
   WHERE (maskatuw(:,:,k))
      array3duw(:,:,k) = hdifscal_fac*array2du*delzatuw(:,:,k)*&
                       & hdifcoefatuw(:,:,k)
   END WHERE
ENDDO k_233

!
!3. Explicit terms
!-----------------
!
!3.1 Fluxes
!----------
! 

WHERE (maskatuw)
   difflux = array3duw*(psiw(1:ncloc+1,1:nrloc,klo:kup)-&
                      & psiw(0:ncloc,1:nrloc,klo:kup))
END WHERE

!
!3.2 Add flux divergence to explicit term
!----------------------------------------
!

k_320: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfw(:,:,k,4) = tridcfw(:,:,k,4) + &
                       & (difflux(2:ncloc+1,:,k)-difflux(1:ncloc,:,k))&
                       & /denom(:,:,k)
   END WHERE
ENDDO k_320

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array2du,array3duw,difflux,denom,maskatuw,mask2duw)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Xdif_at_W

!========================================================================

SUBROUTINE Ydif_at_W(psiw,tridcfw,hdifcoefatvw,klo,kup)
!************************************************************************
!
! *Ydif_at_W* Diffusion term in Y-direction for a quantity at W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_W
!
! Module calls - error_alloc
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1)::psiw
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz+1,4) :: tridcfw
REAL, INTENT(IN), DIMENSION(ncloc,nrloc+1,2:nz) :: hdifcoefatvw

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psiw*        REAL    W-node quantity to be diffused                  [psiw]
!*tridcfw*     REAL    Tridiagonal matrix for implicit (vertical) solution
!*hdifcoefatvw*REAL    Diffusion coefficient at VW-nodes              [m^2/s]
!*klo*         INTEGER Lower vertical array bound
!*kup*         INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2dvw
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatvw
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dvw, denom, difflux

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*difflux*    REAL    Diffusive flux (times factor) at VW-nodes       
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Ydif_at_W'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2dv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (array3dvw(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('array3dvw',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)
ALLOCATE (difflux(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('difflux',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)
ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
ALLOCATE (maskatvw(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('maskatvw',3,(/ncloc,nrloc+1,kup-klo+1/),kndlog)
ALLOCATE (mask2dvw(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('mask2dvw',3,(/ncloc,nrloc+1/),kndlog)

!
!2. Initialise
!-------------
!
!2.1 Mask array
!--------------
!

maskatvw = nodeatvw(1:ncloc,1:nrloc+1,klo:kup).EQ.2
mask2dvw = node2dv(1:ncloc,1:nrloc+1).EQ.2

!
!2.2 Fluxes
!----------
!

WHERE (.NOT.maskatvw)
   difflux = 0.0
END WHERE

!
!2.3 Work space arrays
!---------------------
!

IF (dXregY) THEN
   k_231: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delyatc(1:ncloc,1:nrloc)*delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_231
ELSE
   k_232: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_232
ENDIF

IF (dYregX) THEN
   WHERE (mask2dvw)
      array2dv = delt3d/delyatv(1:ncloc,1:nrloc+1)
   END WHERE
ELSE
   WHERE (mask2dvw)
      array2dv = delt3d*delxatv(1:ncloc,1:nrloc+1)/delyatv(1:ncloc,1:nrloc+1)
   END WHERE
ENDIF
k_233: DO k=klo,kup
   WHERE (maskatvw(:,:,k))
      array3dvw(:,:,k) = hdifscal_fac*array2dv*delzatvw(:,:,k)*&
                       & hdifcoefatvw(:,:,k)
   END WHERE
ENDDO k_233

!
!3. Explicit terms
!-----------------
!
!3.1 Fluxes
!----------
! 

WHERE (maskatvw)
   difflux = array3dvw*(psiw(1:ncloc,1:nrloc+1,klo:kup)-&
                      & psiw(1:ncloc,0:nrloc,klo:kup))
END WHERE

!
!3.2 Add flux divergence to explicit term
!----------------------------------------
!

k_320: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfw(:,:,k,4) = tridcfw(:,:,k,4) + &
                       & (difflux(:,2:nrloc+1,k)-difflux(:,1:nrloc,k))&
                       & /denom(:,:,k)
   END WHERE
ENDDO k_320

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array2dv,array3dvw,difflux,denom,maskatvw,mask2dvw)

CALL log_timer_out(npcc,itm_hdif)


RETURN

END SUBROUTINE Ydif_at_W

!========================================================================

SUBROUTINE Zdif_at_W(psiw,tridcfw,vdifcoefatc,ibcsur,ibcbot,bcsur,bcbot,&
                   & klo,kup)
!************************************************************************
!
! *Zdif_at_W* Vertical diffusion term for a quantity at W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Diffusion_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_W
!
! Module calls - error_alloc
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ibcsur, ibcbot, klo, kup
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1)::psiw
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz+1,4) :: tridcfw
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz) :: vdifcoefatc
REAL, INTENT(IN), DIMENSION(ncloc,nrloc) :: bcsur, bcbot

!
! Name             Type    Purpose
!------------------------------------------------------------------------------
!*psiw*       REAL    W-node quantity to be diffused                   [psiw]
!*tridcfw*    REAL    Tridiagonal matrix for implicit vertical solution
!*vdifcoefatc*REAL    Diffusion coefficient at C-nodes                [m^2/s]
!*ibcsur*     INTEGER Type of surface boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux at the surface)
!                 = 2 => Neumann (prescribed flux at first C-node below the
!                        surface)
!                 = 3 => Dirichlet at the surface
!                 = 4 => Dirichlet at the first W-node below the surface
!*ibcbot*     INTEGER Type of bottom boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux at the bottom)
!                 = 2 => Neumann (prescribed flux at first C-node above the
!                        bottom)
!                 = 3 => Dirichlet at the bottom
!                 = 4 => Dirichlet at the first W-node above the bottom
!*bcsur*      REAL    Imposed surface flux or surface value
!*bcbot*      REAL    Imposed bottom flux or bottom value
!*klo*        INTEGER Lower vertical array bound
!*kup*        INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, kmax, kmin, npcc
REAL :: xexp, ximp
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc1, array2dc2, array2dc3
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3d, difflux

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*difflux*    REAL    Diffusive flux (times factor) at C-nodes       [m^2/psiw]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Zdif_at_W'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2dc1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dc2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dc3(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc3',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array3d(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('array3d',3,(/ncloc,nrloc,nz/),kndrtype)
IF (iopt_vdif_impl.NE.2) THEN
   ALLOCATE (difflux(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('difflux',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

!
!2. Initialise
!-------------
!
!2.1 Upper and lower bounds
!--------------------------
!

kmax = MERGE(nz,nz-1,ibcsur.EQ.3)
kmin = MERGE(1,2,ibcbot.EQ.3)

!
!2.2 Fluxes
!----------
!

IF (iopt_vdif_impl.NE.2) THEN
   WHERE (maskatc_int)
      difflux(:,:,1) = 0.0
      difflux(:,:,nz) = 0.0
   END WHERE
ENDIF

!
!2.3 Work space array
!--------------------
!

k_230: DO k=kmin,kmax
   WHERE (maskatc_int)
      array3d(:,:,k) = vdifcoefatc(:,:,k)/delzatc(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_230

!---time factors
xexp = delt3d*(1.0-theta_vdif)
ximp = delt3d*theta_vdif

!
!3. Explicit fluxes
!------------------
!

IF (iopt_vdif_impl.NE.2) THEN
   k_301: DO k=klo,kup
      WHERE (maskatc_int)
         difflux(:,:,k) = xexp*array3d(:,:,k)*&
                        & (psiw(1:ncloc,1:nrloc,k+1)-psiw(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_301
ENDIF

!
!4. Explicit terms
!-----------------
!

IF (iopt_vdif_impl.NE.2) THEN
   k_410: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfw(:,:,k,4) = tridcfw(:,:,k,4) + &
                          & (difflux(:,:,k)-difflux(:,:,k-1))&
                          & /delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_410
ENDIF

!
!5. Implicit terms
!-----------------
!

IF (iopt_vdif_impl.NE.0) THEN
   WHERE (.NOT.maskatc_int) array2dc1 = 0.0
   k_510: DO k=kmin,kmax
      WHERE (maskatc_int) array2dc1 = ximp*array3d(:,:,k)
!     ---lower fluxes
      IF (k.LT.kup) THEN
         WHERE (maskatc_int)
            array2dc2 = array2dc1/delzatw(1:ncloc,1:nrloc,k+1)
            tridcfw(:,:,k+1,1) = tridcfw(:,:,k+1,1) - array2dc2
            tridcfw(:,:,k+1,2) = tridcfw(:,:,k+1,2) + array2dc2
         END WHERE
      ENDIF
!     ---upper fluxes
      IF (k.GE.klo) THEN
         WHERE (maskatc_int)
            array2dc2 = array2dc1/delzatw(1:ncloc,1:nrloc,k)
            tridcfw(:,:,k,2) = tridcfw(:,:,k,2) + array2dc2
            tridcfw(:,:,k,3) = tridcfw(:,:,k,3) - array2dc2
         END WHERE
      ENDIF
   ENDDO k_510
ENDIF

!
!6. Boundary conditions
!----------------------
!

xexp = 0.5*xexp
ximp = 0.5*ximp

!
!6.1 Surface
!-----------
!
!---Neumann at the surface
IF (ibcsur.EQ.1) THEN
   WHERE (maskatc_int)
      array2dc1 = delzatw(1:ncloc,1:nrloc,nz) + 0.5*delzatc(1:ncloc,1:nrloc,nz)
      array2dc2 = delzatc(1:ncloc,1:nrloc,nz)*vdifcoefatc(:,:,nz-1)&
                & /(delzatw(1:ncloc,1:nrloc,nz)*array2dc1*&
                  & delzatc(1:ncloc,1:nrloc,nz-1))
      array2dc3 = ximp*array2dc2
      tridcfw(:,:,nz,1) = tridcfw(:,:,nz,1) + array2dc3
      tridcfw(:,:,nz,2) = tridcfw(:,:,nz,2) - array2dc3
      tridcfw(:,:,nz,4) = tridcfw(:,:,nz,4) + delt3d*bcsur/array2dc1 +&
                  & xexp*array2dc2*&
                  & (psiw(1:ncloc,1:nrloc,nz)-psiw(1:ncloc,1:nrloc,nz-1))
   END WHERE

!---Neumann at first C-node below the surface
ELSEIF (ibcsur.EQ.2) THEN
   WHERE (maskatc_int)
      tridcfw(:,:,nz,4) = tridcfw(:,:,nz,4) + delt3d*bcsur&
                        & /delzatw(1:ncloc,1:nrloc,nz)
   END WHERE
ENDIF

!
!6.2 Bottom
!----------
!
!---Neumann at the bottom
IF (ibcbot.EQ.1) THEN
   WHERE (maskatc_int)
      array2dc1 = delzatw(1:ncloc,1:nrloc,2)+0.5*delzatc(1:ncloc,1:nrloc,1)
      array2dc2 = delzatc(1:ncloc,1:nrloc,1)*vdifcoefatc(:,:,2)&
                & /(delzatw(1:ncloc,1:nrloc,2)*array2dc1*&
                  & delzatc(1:ncloc,1:nrloc,2))
      array2dc3 = ximp*array2dc2
      tridcfw(:,:,2,2) = tridcfw(:,:,2,2) - array2dc3
      tridcfw(:,:,2,3) = tridcfw(:,:,2,3) + array2dc3
      tridcfw(:,:,2,4) = tridcfw(:,:,2,4) - delt3d*bcbot/array2dc1 -&
           & xexp*array2dc2*(psiw(1:ncloc,1:nrloc,3)-psiw(1:ncloc,1:nrloc,2))
   END WHERE

!---Neumann at first C-node above the bottom
ELSEIF (ibcbot.EQ.2) THEN
   WHERE (maskatc_int)
      tridcfw(:,:,2,4) = tridcfw(:,:,2,4) - delt3d*bcbot&
                       & /delzatw(1:ncloc,1:nrloc,2)
   END WHERE
ENDIF

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE (array2dc1,array2dc2,array2dc3,array3d)
IF (iopt_vdif_impl.NE.2) DEALLOCATE (difflux)

CALL log_timer_out(npcc,itm_vdif)


RETURN

END SUBROUTINE Zdif_at_W
