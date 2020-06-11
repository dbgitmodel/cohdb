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
! *Advection_Terms* Advective terms in transport equations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! $Date: 2017-06-22 11:10:44 +0200 (Thu, 22 Jun 2017) $
!
! $Revision: 1037 $
!
! Description - 
!
! Reference -
!
! Routines - Xadv_at_C, Yadv_at_C, Zadv_at_C, Xadv_at_U_3d, Xadv_at_U_2d,
!            Yadv_at_U_3d, Yadv_at_U_2d, Zadv_at_U, Xadv_at_V_3d,
!            Xadv_at_V_2d, Yadv_at_V_3d, Yadv_at_V_2d, Zadv_at_V,
!            Xadv_at_W, Yadv_at_W, Zadv_at_W  
!
!************************************************************************
!

!========================================================================


SUBROUTINE Xadv_at_C(psic,tridcfc,novars,iopt_adv,nh,klo,kup,psiobu,uadvvel)
!************************************************************************
!
! *Xadv_at_C* Advection term in X-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2
!
! Module calls - error_alloc, tvd_limiter
!
!************************************************************************
!
USE currents
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
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: iopt_adv, klo, kup, nh, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc
REAL, INTENT(IN), DIMENSION(nobu,nz,novars) :: psiobu
REAL, INTENT(IN), DIMENSION(2-nhalo:ncloc+nhalo,nrloc,nz,novars) :: uadvvel

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psic*      REAL    C-node quantity to be advected                     [psic]
!*uadvvel*   REAL    Horizontal velocity (X-dir) for advection         [m/s] 
!*tridcfc*   REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*    INTEGER Number of variables
!*iopt_adv*  INTEGER Switch to select advection scheme
!*nh*        INTEGER Halo size of local arrays
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!*psiobu*    REAL    External open boundary profiles                    [psic]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ic, ii, iiloc, isign, ivar, j, k, npcc
REAL :: xexp
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2du
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2du
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advflux, advflux_lw, advflux_up, &
                             & array3du1, array3du2, cflobu, denom, dflux, omega
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: cfl, usign

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at U-nodes          [psic]
!*advflux_lw* REAL    Lax-Wendroff flux at U-nodes                [psic]*[m/s]
!*advflux_up* REAL    Upwind flux at U-nodes                      [psic]*[m/s]
!*omega*      REAL    TVD limiting factor
!*usign*      REAL    Upwind factor
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xadv_at_C'
CALL log_timer_in(npcc)

!
!1. Set parameters
!-----------------
!

xexp = 0.5*delt3d

!
!2. Allocate arrays
!------------------
!

ALLOCATE (advflux(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('advflux',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)

ALLOCATE (array2du(2-nh:ncloc+nh,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc+2*nh-1,nrloc/),kndrtype)

ALLOCATE (array3du1(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('array3du1',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)

ALLOCATE (cflobu(nobuloc_ext,klo:kup,novars),STAT=errstat)
CALL error_alloc('cflobu',3,(/nobuloc_ext,kup-klo+1,novars/),kndrtype)

ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)

ALLOCATE (maskatu(2-nh:ncloc+nh,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc+2*nh-1,nrloc,kup-klo+1/),kndlog)

ALLOCATE (mask2du(2-nh:ncloc+nh,nrloc),STAT=errstat)
CALL error_alloc('mask2du',2,(/ncloc+2*nh-1,nrloc/),kndlog)

ALLOCATE (usign(2-nh:ncloc+nh,nrloc,klo:kup,novars),STAT=errstat)
CALL error_alloc('usign',4,(/ncloc+2*nh-1,nrloc,kup-klo+1,novars/),kndrtype)

IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) THEN
   ALLOCATE (advflux_up(2-nh:ncloc+nh,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('advflux_up',3,(/ncloc+2*nh-1,nrloc,kup-klo+1/),kndrtype)
ENDIF

IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN
   ALLOCATE (advflux_lw(2-nh:ncloc+nh,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('advflux_lw',3,(/ncloc+2*nh-1,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (cfl(2-nh:ncloc+nh,nrloc,klo:kup,novars),STAT=errstat)
   CALL error_alloc('cfl',4,(/ncloc+2*nh-1,nrloc,kup-klo+1,novars/),kndrtype)
ENDIF

IF (iopt_adv.EQ.3) THEN
   ALLOCATE (array3du2(ncloc+1,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('array3du2',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (dflux(2-nh:ncloc+nh,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('dflux',3,(/ncloc+2*nh-1,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (omega(ncloc+1,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('omega',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)
ENDIF

!
!3. Initialise arrays
!--------------------
!
!3.1 Mask array
!--------------
!

maskatu = nodeatu(2-nh:ncloc+nh,1:nrloc,klo:kup).EQ.2
mask2du = node2du(2-nh:ncloc+nh,1:nrloc).EQ.2

!
!3.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) THEN 
   advflux_up = 0.0
ENDIF
IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!3.3 Metric factor in advective flux
!-----------------------------------
!

IF (dYregX) THEN
   WHERE (maskatu(1:ncloc+1,:,:))
      array3du1 = xexp*delzatu(1:ncloc+1,:,klo:kup)
   END WHERE
ELSE
   WHERE (mask2du(1:ncloc+1,:))
      array2du(1:ncloc+1,:) = xexp*delyatu(1:ncloc+1,1:nrloc)
   END WHERE
   k_330: DO k=klo,kup
      WHERE (maskatu(1:ncloc+1,:,k))
         array3du1(:,:,k) = array2du(1:ncloc+1,:)*delzatu(1:ncloc+1,:,k)
      END WHERE
   ENDDO k_330
ENDIF

!
!3.4 Sign of velocity
!--------------------
!

ivar_340: DO ivar=1,novars
   WHERE (nodeatu(2-nh:ncloc+nh,1:nrloc,klo:kup).GT.1)
      usign(:,:,:,ivar) = SIGN(1.0,uadvvel(2-nh:ncloc+nh,:,klo:kup,ivar))
   END WHERE
ENDDO ivar_340

!
!3.5 CFL number
!--------------
!

IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN
   WHERE (mask2du)
      array2du = delt3d/delxatu(2-nh:ncloc+nh,1:nrloc)
   END WHERE
   ivar_320: DO ivar=1,novars
   k_320: DO k=klo,kup
      WHERE (maskatu(:,:,k))
         cfl(:,:,k,ivar) = array2du*uadvvel(2-nh:ncloc+nh,:,k,ivar)
      END WHERE
   ENDDO k_320
   ENDDO ivar_320
ENDIF

!
!3.6 Work space arrays
!---------------------
!

IF (dYregX) THEN
   k_361: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delxatc(1:ncloc,1:nrloc)*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_361
ELSE
   k_362: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_362
ENDIF

!
!3.7 Open boundary conditions
!----------------------------
!

ivar_370: DO ivar=1,novars
iiloc_370: DO iiloc=1,nobuloc_ext
   i = iobuloc(iiloc); j = jobuloc(iiloc)
   ii = indexobu(iiloc); ic = MERGE(i,i-1,westobu(ii))
   IF (ic.GE.1.AND.ic.LE.ncloc) THEN
      isign = MERGE(-1,1,westobu(ii))
      IF (dYregX) THEN
         cflobu(iiloc,:,ivar) = isign*xexp*delzatu(i,j,klo:kup)*&
                              & uadvvel(i,j,klo:kup,ivar)/denom(ic,j,:)
      ELSE
         cflobu(iiloc,:,ivar) = isign*xexp*delyatu(i,j)*delzatu(i,j,klo:kup)*&
                              & uadvvel(i,j,klo:kup,ivar)/denom(ic,j,:)
      ENDIF
   ENDIF
ENDDO iiloc_370
ENDDO ivar_370

!
!4. Explicit terms
!-----------------
!

ivar_400: DO ivar=1,novars

!
!4.1 Upwind fluxes
!-----------------
!

   IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) THEN

      k_410:DO k=klo,kup
         WHERE (maskatu(:,:,k))
            advflux_up(:,:,k) = uadvvel(2-nh:ncloc+nh,:,k,ivar)*&
                 & ((1.0+usign(:,:,k,ivar))*&
                 & psic(1-nh:ncloc+nh-1,1:nrloc,k,ivar)&
                 & +(1.0-usign(:,:,k,ivar))*psic(2-nh:ncloc+nh,1:nrloc,k,ivar))
         END WHERE
      ENDDO k_410

!     ---non-TVD case
      IF (iopt_adv.EQ.1) THEN
         WHERE (maskatu(1:ncloc+1,:,:))
            advflux = array3du1*advflux_up(1:ncloc+1,:,:)
         END WHERE
      ENDIF

   ENDIF

!
!4.2 Lax-Wendroff fluxes
!-----------------------
!

   IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN


      k_420: DO k=klo,kup
         WHERE (maskatu(:,:,k))
            advflux_lw(:,:,k) = uadvvel(2-nh:ncloc+nh,:,k,ivar)*&
               & ((1.0+cfl(:,:,k,ivar))*psic(1-nh:ncloc+nh-1,1:nrloc,k,ivar)&
               & +(1.0-cfl(:,:,k,ivar))*psic(2-nh:ncloc+nh,1:nrloc,k,ivar))
         END WHERE
      ENDDO k_420

!     ---non-TVD case
      IF (iopt_adv.EQ.2) THEN
         WHERE (maskatu(1:ncloc+1,:,:))
            advflux = array3du1*advflux_lw(1:ncloc+1,:,:)
         END WHERE
      ENDIF

   ENDIF

!
!4.3 TVD scheme
!--------------
!

   IF (iopt_adv.EQ.3) THEN
   
!     ---flux differences
      WHERE (maskatu)
         dflux = advflux_lw - advflux_up
      END WHERE

!     ---flux ratios
      WHERE (maskatu(1:ncloc+1,:,:)) array3du2 = 0.0
      k_430: DO k=klo,kup
         WHERE (ABS(dflux(1:ncloc+1,:,k)).GT.eps_adv)
            array3du2(:,:,k) = 0.5*(&
                             & (1.0+usign(1:ncloc+1,:,k,ivar))*&
                             & dflux(0:ncloc,:,k)+&
                             & (1.0-usign(1:ncloc+1,:,k,ivar))*&
                             & dflux(2:ncloc+2,:,k))/dflux(1:ncloc+1,:,k)
         END WHERE
      ENDDO k_430

!     ---TVD limiter
      omega = tvd_limiter(array3du2,maskatu(1:ncloc+1,:,:))

!     ---fluxes
      WHERE (maskatu(1:ncloc+1,:,:))
         advflux = array3du1*(advflux_up(1:ncloc+1,:,:) + &
                            & omega*dflux(1:ncloc+1,:,:))
      END WHERE

   ENDIF

!
!4.4 Advective term
!------------------
!

   k_440: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) - &
                    & (advflux(2:ncloc+1,:,k)-advflux(1:ncloc,:,k))/denom(:,:,k)
      END WHERE
   ENDDO k_440

ENDDO ivar_400

!
!5. Open boundary fluxes
!-----------------------
!

ivar_510: DO ivar=1,novars
iiloc_510: DO iiloc=1,nobuloc_ext
   i = iobuloc(iiloc); j = jobuloc(iiloc)
   ii = indexobu(iiloc); ic = MERGE(i,i-1,westobu(ii))
   IF (ic.GE.1.AND.ic.LE.ncloc) THEN
      isign = MERGE(-1,1,westobu(ii))
      WHERE (ABS(psiobu(ii,klo:kup,ivar)-float_fill).GT.float_fill_eps)
         tridcfc(ic,j,klo:kup,4,ivar) = tridcfc(ic,j,klo:kup,4,ivar) - &
                                      & cflobu(iiloc,:,ivar)*&
                                      & (1.0+isign*usign(i,j,:,ivar))*&
                                      & psic(ic,j,klo:kup,ivar)
         tridcfc(ic,j,klo:kup,4,ivar) = tridcfc(ic,j,klo:kup,4,ivar) - &
                                      & cflobu(iiloc,:,ivar)*&
                                      & (1.0-isign*usign(i,j,:,ivar))*&
                                      & psiobu(ii,klo:kup,ivar)
      ELSEWHERE
         tridcfc(ic,j,klo:kup,4,ivar) = tridcfc(ic,j,klo:kup,4,ivar) - &
                                      & 2.0*cflobu(iiloc,:,ivar)*&
                                      & psic(ic,j,klo:kup,ivar)
      END WHERE
   ENDIF
ENDDO iiloc_510
ENDDO ivar_510

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatu,mask2du)
DEALLOCATE (advflux,array2du,array3du1,cflobu,denom,usign)
IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) DEALLOCATE (advflux_up)
IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) DEALLOCATE (advflux_lw,cfl)
IF (iopt_adv.EQ.3) DEALLOCATE (array3du2,dflux,omega)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Xadv_at_C

!========================================================================

SUBROUTINE Yadv_at_C(psic,tridcfc,novars,iopt_adv,nh,klo,kup,psiobv,vadvvel)
!************************************************************************
!
! *Yadv_at_C* Advection term in Y-direction for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2
!
! Module calls - error_alloc, tvd_limiter
!
!************************************************************************
!
USE currents
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
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iopt_adv, klo, kup, nh, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc
REAL, INTENT(IN), DIMENSION(nobv,nz,novars) :: psiobv
REAL, INTENT(IN),DIMENSION(ncloc,2-nhalo:nrloc+nhalo,nz,novars) ::  vadvvel

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psic*      REAL    C-node quantity to be advected                     [psic]
!*vadvvel*   REAL    Horizontal velocity (Y-dir) for advection           [m/s] 
!*tridcfc*   REAL    Tridiagonal matrix for implicit (vertical) solution
!*novars*    INTEGER Number of variables
!*iopt_adv*  INTEGER Switch to select advection scheme
!*nh*        INTEGER Halo size of local arrays
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!*psiobv*    REAL    External open boundary profiles                    [psic]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ivar, j, jc, jj, jjloc, jsign, k, npcc
REAL :: xexp
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2dv
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advflux, advflux_lw, advflux_up, &
                             & array3dv1, array3dv2, cflobv, denom, dflux, omega
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: cfl, vsign

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at V-nodes          [psic]
!*advflux_lw* REAL    Lax_Wendroff flux at V-nodes                [psic]*[m/s]
!*advflux_up* REAL    Upwind flux at V-nodes                      [psic]*[m/s]
!*omega*      REAL    TVD limiting factor
!*vsign*      REAL    Upwind factor
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Yadv_at_C'
CALL log_timer_in(npcc)

!
!1. Set parameters
!-----------------
!

xexp = 0.5*delt3d

!
!2. Allocate arrays
!------------------
!

ALLOCATE (advflux(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('advflux',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)

ALLOCATE (array2dv(ncloc,2-nh:nrloc+nh),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc,nrloc+2*nh-1/),kndrtype)

ALLOCATE (array3dv1(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('array3dv1',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)

ALLOCATE (cflobv(nobvloc_ext,klo:kup,novars),STAT=errstat)
CALL error_alloc('cflobv',3,(/nobvloc_ext,kup-klo+1,novars/),kndrtype)

ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)

ALLOCATE (mask2dv(ncloc,2-nh:nrloc+nh),STAT=errstat)
CALL error_alloc('mask2dv',2,(/ncloc,nrloc+2*nh-1/),kndlog)

ALLOCATE (maskatv(ncloc,2-nh:nrloc+nh,klo:kup),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc+2*nh-1,kup-klo+1/),kndlog)

ALLOCATE (vsign(ncloc,2-nh:nrloc+nh,klo:kup,novars),STAT=errstat)
CALL error_alloc('vsign',4,(/ncloc,nrloc+2*nh-1,kup-klo+1,novars/),kndrtype)

IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) THEN
   ALLOCATE (advflux_up(ncloc,2-nh:nrloc+nh,klo:kup),STAT=errstat)
   CALL error_alloc('advflux_up',3,(/ncloc,nrloc+2*nh-1,kup-klo+1/),kndrtype)
ENDIF

IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN
   ALLOCATE (advflux_lw(ncloc,2-nh:nrloc+nh,klo:kup),STAT=errstat)
   CALL error_alloc('advflux_lw',3,(/ncloc,nrloc+2*nh-1,kup-klo+1/),kndrtype)
   ALLOCATE (cfl(ncloc,2-nh:nrloc+nh,klo:kup,novars),STAT=errstat)
   CALL error_alloc('cfl',4,(/ncloc,nrloc+2*nh-1,kup-klo+1,novars/),kndrtype)
ENDIF

IF (iopt_adv.EQ.3) THEN
   ALLOCATE (array3dv2(ncloc,nrloc+1,klo:kup),STAT=errstat)
   CALL error_alloc('array3dv2',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)
   ALLOCATE (dflux(ncloc,2-nh:nrloc+nh,klo:kup),STAT=errstat)
   CALL error_alloc('dflux',3,(/ncloc,nrloc+2*nh-1,kup-klo+1/),kndrtype)
   ALLOCATE (omega(ncloc,nrloc+1,klo:kup),STAT=errstat)
   CALL error_alloc('omega',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)
ENDIF

!
!3. Initialise arrays
!--------------------
!
!3.1 Mask array
!--------------
!

maskatv = nodeatv(1:ncloc,2-nh:nrloc+nh,klo:kup).EQ.2
mask2dv = node2dv(1:ncloc,2-nh:nrloc+nh).EQ.2

!
!3.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) THEN
   advflux_up = 0.0
ENDIF
IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!3.3 Metric factor in advective flux
!-----------------------------------
!


IF (dXregY) THEN
   WHERE (maskatv(:,1:nrloc+1,:))
      array3dv1 = xexp*delzatv(:,1:nrloc+1,klo:kup)
   END WHERE
ELSE
   WHERE (mask2dv(:,1:nrloc+1))
      array2dv(:,1:nrloc+1) = xexp*delxatv(1:ncloc,1:nrloc+1)
   END WHERE
   k_330: DO k=klo,kup
      WHERE (maskatv(:,1:nrloc+1,k))
         array3dv1(:,:,k) = array2dv(:,1:nrloc+1)*delzatv(:,1:nrloc+1,k)
      END WHERE
   ENDDO k_330
ENDIF

!
!3.4 Sign of velocity
!--------------------
!
ivar_340: DO ivar=1,novars
   WHERE (nodeatv(1:ncloc,2-nh:nrloc+nh,klo:kup).GT.1)
      vsign(:,:,:,ivar) = SIGN(1.0,vadvvel(:,2-nh:nrloc+nh,klo:kup,ivar))
   END WHERE
ENDDO ivar_340

!
!3.5 CFL number
!--------------
!

IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN
   WHERE (mask2dv)
      array2dv = delt3d/delyatv(1:ncloc,2-nh:nrloc+nh)
   END WHERE
   ivar_350: DO ivar=1,novars
   k_350: DO k=klo,kup
      WHERE (maskatv(:,:,k))
         cfl(:,:,k,ivar) = array2dv*vadvvel(:,2-nh:nrloc+nh,k,ivar)
      END WHERE
   ENDDO k_350
   ENDDO ivar_350
ENDIF

!
!3.6 Work space arrays
!---------------------
!

IF (dXregY) THEN
   k_361: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delyatc(1:ncloc,1:nrloc)*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_361
ELSE
   k_362: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_362
ENDIF

!
!3.7 Open boundary conditions
!----------------------------
!

ivar_370: DO ivar=1,novars
jjloc_370: DO jjloc=1,nobvloc_ext
   i = iobvloc(jjloc); j = jobvloc(jjloc)
   jj = indexobv(jjloc); jc = MERGE(j,j-1,soutobv(jj))
   IF (jc.GE.1.AND.jc.LE.nrloc) THEN
      jsign = MERGE(-1,1,soutobv(jj))
      IF (dXregY) THEN
         cflobv(jjloc,:,ivar) = jsign*xexp*delzatv(i,j,klo:kup)*&
                              & vadvvel(i,j,klo:kup,ivar)/denom(i,jc,:)
      ELSE
         cflobv(jjloc,:,ivar) = jsign*xexp*delxatv(i,j)*delzatv(i,j,klo:kup)*&
                              & vadvvel(i,j,klo:kup,ivar)/denom(i,jc,:)
      ENDIF
   ENDIF
ENDDO jjloc_370
ENDDO ivar_370

!
!4. Explicit terms
!-----------------
!

ivar_400: DO ivar=1,novars

!
!4.1 Upwind fluxes
!-----------------
!

   IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) THEN

      k_410:DO k=klo,kup
         WHERE (maskatv(:,:,k))
            advflux_up(:,:,k) = vadvvel(:,2-nh:nrloc+nh,k,ivar)*&
                 & ((1.0+vsign(:,:,k,ivar))*&
                 & psic(1:ncloc,1-nh:nrloc+nh-1,k,ivar)&
                 & +(1.0-vsign(:,:,k,ivar))*psic(1:ncloc,2-nh:nrloc+nh,k,ivar))
         END WHERE
      END DO k_410

!     ---non-TVD case
      IF (iopt_adv.EQ.1) THEN
         WHERE (maskatv(:,1:nrloc+1,:))
            advflux = array3dv1*advflux_up(:,1:nrloc+1,:)
         END WHERE
      ENDIF

   ENDIF

!
!4.2 Lax-Wendroff fluxes
!-----------------------
!

   IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN

      k_420: DO k=klo,kup
         WHERE (maskatv(:,:,k))
            advflux_lw(:,:,k) = vadvvel(:,2-nh:nrloc+nh,k,ivar)*&
               & ((1.0+cfl(:,:,k,ivar))*psic(1:ncloc,1-nh:nrloc+nh-1,k,ivar)&
               & +(1.0-cfl(:,:,k,ivar))*psic(1:ncloc,2-nh:nrloc+nh,k,ivar))
         END WHERE
      ENDDO k_420

!     ---non-TVD case
      IF (iopt_adv.EQ.2) THEN
         WHERE (maskatv(:,1:nrloc+1,:))
            advflux = array3dv1*advflux_lw(:,1:nrloc+1,:)
         END WHERE
      ENDIF

   ENDIF

!
!4.3 TVD scheme
!--------------
!

   IF (iopt_adv.EQ.3) THEN

!     ---flux differences
      WHERE (maskatv)
         dflux = advflux_lw - advflux_up
      END WHERE
   
!     ---flux ratios
      WHERE (maskatv(:,1:nrloc+1,:)) array3dv2 = 0.0
      k_430: DO k=klo,kup
         WHERE (ABS(dflux(:,1:nrloc+1,k)).GT.eps_adv)
            array3dv2(:,:,k) = 0.5*(&
                             & (1.0+vsign(:,1:nrloc+1,k,ivar))*&
                             & dflux(:,0:nrloc,k)+&
                             & (1.0-vsign(:,1:nrloc+1,k,ivar))*&
                             & dflux(:,2:nrloc+2,k))/dflux(:,1:nrloc+1,k)
         END WHERE
      ENDDO k_430

!     ---TVD limiter
      omega = tvd_limiter(array3dv2,maskatv(:,1:nrloc+1,:))

!     ---fluxes
      WHERE (maskatv(:,1:nrloc+1,:))
         advflux = array3dv1*(advflux_up(:,1:nrloc+1,:) + &
                            & omega*dflux(:,1:nrloc+1,:))
      END WHERE

   ENDIF

!
!4.4 Advective term
!-------------------
!

   k_440: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) - &
                    & (advflux(:,2:nrloc+1,k)-advflux(:,1:nrloc,k))/denom(:,:,k)
      END WHERE
   ENDDO k_440

ENDDO ivar_400

!
!5. Open boundary fluxes
!-----------------------
!

ivar_510: DO ivar=1,novars
jjloc_510: DO jjloc=1,nobvloc_ext
   i = iobvloc(jjloc); j = jobvloc(jjloc)
   jj = indexobv(jjloc); jc = MERGE(j,j-1,soutobv(jj))
   IF (jc.GE.1.AND.jc.LE.nrloc) THEN
      jsign = MERGE(-1,1,soutobv(jj))
      WHERE (ABS(psiobv(jj,klo:kup,ivar)-float_fill).GT.float_fill_eps)
         tridcfc(i,jc,klo:kup,4,ivar) = tridcfc(i,jc,klo:kup,4,ivar) - &
                                      & cflobv(jjloc,:,ivar)*&
                                      & (1.0+jsign*vsign(i,j,:,ivar))*&
                                      & psic(i,jc,klo:kup,ivar)
         tridcfc(i,jc,klo:kup,4,ivar) = tridcfc(i,jc,klo:kup,4,ivar) - &
                                      & cflobv(jjloc,:,ivar)*&
                                      & (1.0-jsign*vsign(i,j,:,ivar))*&
                                      & psiobv(jj,klo:kup,ivar)
      ELSEWHERE
         tridcfc(i,jc,klo:kup,4,ivar) = tridcfc(i,jc,klo:kup,4,ivar) - &
                                      & 2.0*cflobv(jjloc,:,ivar)*&
                                      & psic(i,jc,klo:kup,ivar)
      END WHERE
   ENDIF
ENDDO jjloc_510
ENDDO ivar_510

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatv,mask2dv)
DEALLOCATE (advflux,array2dv,array3dv1,cflobv,denom,vsign)
IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) DEALLOCATE (advflux_up)
IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) DEALLOCATE (advflux_lw,cfl)
IF (iopt_adv.EQ.3) DEALLOCATE (array3dv2,dflux,omega)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Yadv_at_C

!========================================================================

SUBROUTINE Zadv_at_C(psic,tridcfc,wadvatw,novars,iopt_adv,klo,kup)
!************************************************************************
!
! *Zadv_at_C* Vertical advection term for a quantity at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2
!
! Module calls - error_alloc, tvd_limiter
!
!************************************************************************
!
USE currents
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
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: iopt_adv, klo, kup, novars
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,&
                          & nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4,novars) :: tridcfc
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1,novars) :: wadvatw

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psic*      REAL    C-node quantity to be advected                     [psic]
!*tridcfc*   REAL    Tridiagonal matrix for implicit vertical solution
!*wadvatw*   REAL    Advective velocity at W-nodes                       [m/s]
!*novars*    INTEGER Number of variables
!*iopt_adv*  INTEGER Switch to select advection scheme
!*kbounds*   INTEGER Vertical bounds
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iopt_impl, ivar, k, npcc
REAL :: xexp, ximp
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc1, array2dc2, cfl, upwind
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advflux, advflux_ce, advflux_lw, &
                                           & advflux_up, dflux, omega, wsign

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at W-nodes          [psic]
!*advflux_ce* REAL    Central flux at W-nodes                     [psic]*[m/s]
!*advflux_lw* REAL    Lax_Wendroff flux at W-nodes                [psic]*[m/s]
!*advflux_up* REAL    Upwind flux at W-nodes                      [psic]*[m/s]
!*omega*      REAL    TVD limiting factor
!*wsign*      REAL    Upwind factor
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Zadv_at_C'
CALL log_timer_in(npcc)

!
!1. Set parameters
!-----------------
!

IF (iopt_adv.EQ.2) THEN
   xexp = 0.5*delt3d
   ximp = 0.0
   iopt_impl = 0
ELSE
   xexp = 0.5*delt3d*(1.0-theta_vadv)
   ximp = 0.5*delt3d*theta_vadv
   iopt_impl = iopt_vadv_impl
ENDIF

!
!2. Allocate arrays
!------------------
!

IF (iopt_impl.NE.2) THEN
   ALLOCATE (advflux(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('advflux',3,(/ncloc,nrloc,nz+1/),kndrtype)
ENDIF

IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) THEN
   ALLOCATE (wsign(ncloc,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('wsign',3,(/ncloc,nrloc,nz-1/),kndrtype)
   IF (iopt_impl.NE.2) THEN
      ALLOCATE (advflux_up(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('advflux_up',3,(/ncloc,nrloc,nz+1/),kndrtype)
   ENDIF
ENDIF

IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      ALLOCATE (advflux_ce(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('advflux_ce',3,(/ncloc,nrloc,nz+1/),kndrtype)
   ELSE
      ALLOCATE (advflux_lw(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('advflux_lw',2,(/ncloc,nrloc,nz+1/),kndrtype)
      ALLOCATE (cfl(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('cfl',2,(/ncloc,nrloc/),kndrtype)
   ENDIF
ENDIF

IF (iopt_adv.EQ.3) THEN
   ALLOCATE (dflux(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('dflux',2,(/ncloc,nrloc,nz+1/),kndrtype)
   ALLOCATE (omega(ncloc,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('omega',2,(/ncloc,nrloc,nz-1/),kndrtype)
ENDIF

IF (iopt_impl.GT.0) THEN
   ALLOCATE (upwind(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('upwind',2,(/ncloc,nrloc/),kndrtype)
ENDIF


IF (iopt_adv.EQ.3.OR.iopt_impl.GT.0) THEN
   ALLOCATE (array2dc1(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc1',2,(/ncloc,nrloc/),kndrtype)
ENDIF

IF (iopt_scal_depos.EQ.2) THEN
   ALLOCATE (array2dc2(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc2',2,(/ncloc,nrloc/),kndrtype)
ENDIF
   
!
!3. Initialise arrays
!--------------------
!
!3.1 Fluxes
!----------
!

advflux = 0.0
IF (iopt_impl.NE.2) THEN
   IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) THEN
      advflux_up = 0.0
   ENDIF
ENDIF

IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      advflux_ce = 0.0
   ELSE
      advflux_lw = 0.0
   ENDIF
ENDIF

IF (iopt_adv.EQ.3) dflux = 0.0

!
!3.2 Work space arrays
!---------------------
!

IF (iopt_scal_depos.EQ.2) THEN
   WHERE (maskatc_int)
      array2dc2 = 0.5*delzatc(1:ncloc,1:nrloc,1)/delzatw(1:ncloc,1:nrloc,2)
   END WHERE
ENDIF

!
!4. Advective terms
!------------------
!

ivar_400: DO ivar=1,novars

!
!4.1 Sign of velocity
!--------------------
!

   IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) THEN
      k_410: DO k=2,nz
         WHERE (maskatc_int)
            wsign(:,:,k) = SIGN(1.0,wadvatw(:,:,k,ivar))
         END WHERE
      ENDDO k_410
   ENDIF

!
!4.2 Explicit terms
!-----------------
!
!4.2.1 Upwind fluxes
!-------------------
!

   IF ((iopt_adv.EQ.1.AND.iopt_impl.NE.2).OR.iopt_adv.EQ.3) THEN

      k_4211: DO k=2,nz
         WHERE (maskatc_int)
            advflux_up(:,:,k) = wadvatw(:,:,k,ivar)*&
                 & ((1.0+wsign(:,:,k))*psic(1:ncloc,1:nrloc,k-1,ivar)&
                 & +(1.0-wsign(:,:,k))*psic(1:ncloc,1:nrloc,k,ivar))
         END WHERE
      ENDDO k_4211

!     ---non-TVD case
      IF (iopt_adv.NE.3) THEN
         k_4212: DO k=2,nz
            WHERE (maskatc_int)
               advflux(:,:,k) = xexp*advflux_up(:,:,k)
            END WHERE
         ENDDO k_4212
      ENDIF

   ENDIF

!
!4.2.2 Central fluxes
!--------------------
!

   IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN

      IF (iopt_impl.GT.0) THEN

         k_4221: DO k=2,nz
            WHERE (maskatc_int)
               advflux_ce(:,:,k) = wadvatw(:,:,k,ivar)*&
                                & (psic(1:ncloc,1:nrloc,k-1,ivar)&
                                & +psic(1:ncloc,1:nrloc,k,ivar))
            END WHERE
         ENDDO k_4221

!        ---non-TVD case
         IF (iopt_adv.EQ.2.AND.iopt_impl.NE.2) THEN
            k_4222: DO k=2,nz
               WHERE (maskatc_int) 
                  advflux(:,:,k) = xexp*advflux_ce(:,:,k)
               END WHERE
            ENDDO k_4222
         ENDIF

!
!4.2.3 Lax-Wendroff fluxes
!-------------------------
!

      ELSE

         k_4231: DO k=2,nz
            WHERE (maskatc_int)
               cfl = delt3d*wadvatw(:,:,k,ivar)/delzatw(1:ncloc,1:nrloc,k)
               advflux_lw(:,:,k) = wadvatw(:,:,k,ivar)*&
                          & ((1.0+cfl)*psic(1:ncloc,1:nrloc,k-1,ivar)&
                          & +(1.0-cfl)*psic(1:ncloc,1:nrloc,k,ivar))
            END WHERE
         ENDDO k_4231

!        ---non-TVD case
         IF (iopt_adv.EQ.2) THEN
            k_4232: DO k=2,nz
               WHERE (maskatc_int)
                  advflux(:,:,k) = xexp*advflux_lw(:,:,k)
               END WHERE
            ENDDO k_4232
         ENDIF
      ENDIF

   ENDIF

!
!4.2.4 TVD scheme
!----------------
!

   IF (iopt_adv.EQ.3) THEN

!     ---flux differences
      k_4241: DO k=2,nz
         IF (iopt_impl.NE.0) THEN
            WHERE (maskatc_int)
               dflux(:,:,k) = advflux_ce(:,:,k) - advflux_up(:,:,k)
            END WHERE
         ELSE
            WHERE (maskatc_int)
               dflux(:,:,k) = advflux_lw(:,:,k) - advflux_up(:,:,k)
            END WHERE
         ENDIF
      ENDDO k_4241

      k_4242: DO k=2,nz
!        --flux ratios
         WHERE (maskatc_int) array2dc1 = 0.0
         WHERE (ABS(dflux(:,:,k)).GT.eps_adv)
            array2dc1 = 0.5*((1.0+wsign(:,:,k))*dflux(:,:,k-1)+&
                           & (1.0-wsign(:,:,k))*dflux(:,:,k+1))&
                           & /dflux(:,:,k)
         END WHERE
!        --TVD limiter
         omega(:,:,k) = tvd_limiter(array2dc1,maskatc_int) 
!        --fluxes
         IF (iopt_impl.NE.2) THEN
            WHERE (maskatc_int)
               advflux(:,:,k) = xexp*(advflux_up(:,:,k)+&
                                    & omega(:,:,k)*dflux(:,:,k))
            END WHERE
         ENDIF
      ENDDO k_4242

   ENDIF

!
!4.2.5 Bottom flux
!-----------------
!

   IF (iopt_impl.NE.2.AND.klo.EQ.1) THEN

!     ---first order
      IF (iopt_scal_depos.EQ.1) THEN
        WHERE (wadvatw(:,:,1,ivar).LT.0.0)
           advflux(:,:,1) = (2.0*xexp)*wadvatw(:,:,1,ivar)*&
                          & psic(1:ncloc,1:nrloc,1,ivar)
        END WHERE

!     ---second order
      ELSEIF (iopt_scal_depos.EQ.2) THEN
        WHERE (wadvatw(:,:,1,ivar).LT.0.0)
           advflux(:,:,1) = (2.0*xexp)*wadvatw(:,:,1,ivar)*&
                          & (psic(1:ncloc,1:nrloc,1,ivar)*(array2dc1+1.0) - &
                          &  psic(1:ncloc,1:nrloc,2,ivar)*array2dc1)
        END WHERE
      ENDIF
   ENDIF

!
!4.2.6 Advective term
!------------------
!

   IF (iopt_impl.NE.2) THEN
      k_426: DO k=klo,kup
         WHERE (maskatc_int)
            tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) - &
                                 & (advflux(:,:,k+1)-advflux(:,:,k))&
                                 & /delzatc(1:ncloc,1:nrloc,k)
         END WHERE
      ENDDO k_426
   ENDIF

!
!4.3. Implicit terms
!-------------------
!

   IF (iopt_impl.GT.0) THEN

      k_430: DO k=2,nz

!        ---upwind factor
         IF (iopt_adv.EQ.1) THEN
            WHERE (maskatc_int) upwind = wsign(:,:,k)
         ELSEIF (iopt_adv.EQ.2) THEN
            WHERE (maskatc_int) upwind = 0.0
         ELSEIF (iopt_adv.EQ.3) THEN
            WHERE (maskatc_int) upwind = (1.0-omega(:,:,k))*wsign(:,:,k)
         ENDIF
         WHERE (maskatc_int)
            array2dc1 = ximp*wadvatw(:,:,k,ivar)
         END WHERE

!        ---lower fluxes
         IF (k.LE.kup) THEN
            WHERE (maskatc_int)
               tridcfc(:,:,k,1,ivar) = tridcfc(:,:,k,1,ivar) - array2dc1*&
                                     & (1.0+upwind)/delzatc(1:ncloc,1:nrloc,k)
               tridcfc(:,:,k,2,ivar) = tridcfc(:,:,k,2,ivar) - array2dc1*&
                                     & (1.0-upwind)/delzatc(1:ncloc,1:nrloc,k)
            END WHERE
         ENDIF

!        ---upper fluxes
         IF (k.GT.klo) THEN
            WHERE (maskatc_int)
               tridcfc(:,:,k-1,2,ivar) = tridcfc(:,:,k-1,2,ivar) + array2dc1*&
                                     & (1.0+upwind)/delzatc(1:ncloc,1:nrloc,k-1)
               tridcfc(:,:,k-1,3,ivar) = tridcfc(:,:,k-1,3,ivar) + array2dc1*&
                                     & (1.0-upwind)/delzatc(1:ncloc,1:nrloc,k-1)
            END WHERE
         ENDIF

      ENDDO k_430

!     ---bottom flux
      IF (klo.EQ.1) THEN
!        --first order
         IF (iopt_scal_depos.EQ.1) THEN
            WHERE (wadvatw(:,:,1,ivar).LT.0.0)
               tridcfc(:,:,1,2,ivar) = tridcfc(:,:,1,2,ivar) - &
                     & (2.0*ximp)*wadvatw(:,:,1,ivar)/delzatc(1:ncloc,1:nrloc,1)
            END WHERE

!        --second order
         ELSEIF (iopt_scal_depos.EQ.2) THEN
            WHERE (wadvatw(:,:,1,ivar).LT.0.0)
               array2dc1 = (2.0*ximp)*wadvatw(:,:,1,ivar)/&
                          & delzatc(1:ncloc,1:nrloc,1)
               tridcfc(:,:,1,2,ivar) = tridcfc(:,:,1,2,ivar) - &
                                     & array2dc1*(1.0+array2dc2)
               tridcfc(:,:,1,3,ivar) = tridcfc(:,:,1,3,ivar) + &
                                     & array2dc1*array2dc2
             END WHERE
         ENDIF
      ENDIF

   ENDIF

ENDDO ivar_400

!
!5. Deallocate arrays
!--------------------
!

IF (iopt_impl.NE.2) DEALLOCATE (advflux)
IF (iopt_adv.EQ.1.OR.iopt_adv.EQ.3) THEN
   DEALLOCATE (wsign)
   IF (iopt_impl.NE.2) DEALLOCATE (advflux_up)
ENDIF
IF (iopt_adv.EQ.2.OR.iopt_adv.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      DEALLOCATE (advflux_ce)
   ELSE
      DEALLOCATE (advflux_lw,cfl)
   ENDIF
ENDIF
IF (iopt_adv.EQ.3) DEALLOCATE (dflux,omega)
IF (iopt_impl.GT.0) DEALLOCATE (upwind)
IF (iopt_adv.EQ.3.OR.iopt_impl.GT.0) DEALLOCATE (array2dc1)
IF (iopt_scal_depos.EQ.2) DEALLOCATE (array2dc2)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Zadv_at_C

!========================================================================

SUBROUTINE Xadv_at_U_3d(psiu,xadvatu3d,nh,vintflag)
!************************************************************************
!
! *Xadv_at_U_3d* Advection term in X-direction for the 3-D current at U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_U_3d
!
! Module calls - error_alloc, tvd_limiter, Carr_at_U, Uarr_at_C, Varr_at_U
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Carr_at_U, Uarr_at_C, Varr_at_U
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: vintflag
INTEGER, INTENT(IN) :: nh
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiu
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,nz) :: xadvatu3d

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiu*      REAL    U-node quantity to be advected                    [psiu]
!*xadvatu3d* REAL    Advective term (times delz)                   [m*psiu/s]
!*nh*        INTEGER Halo size of local arrays
!*vintflag*  LOGICAL Update vertically integrated advective term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: info = .FALSE.
INTEGER :: k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: advflux, advflux_lw, advflux_up,&
        & array2dc1, array2dc2, cfl, denom, dflux, fcurvatu, omega,  usign, &
        & vadvatu, vstokesatu
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: uadvatc

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at C-nodes          [psiu]
!*advflux_lw* REAL    Lax_Wendroff flux at C-nodes                [psiu]*[m/s]
!*advflux_up* REAL    Upwind flux at C-nodes                      [psiu]*[m/s]
!*omega*      REAL    TVD limiting factor
!*usign*      REAL    Upwind factor
!*uadvatc*    REAL    Advective velocity at C-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xadv_at_U_3d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatc(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc+2*nh-1,nrloc/),kndlog)

ALLOCATE (maskatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc,nrloc,nz/),kndlog)

ALLOCATE (advflux(0:ncloc,nrloc),STAT=errstat)
CALL error_alloc('advflux',2,(/ncloc+1,nrloc/),kndrtype)

ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE (uadvatc(1-nh:ncloc+nh-1,nrloc,nz),STAT=errstat)
CALL error_alloc('uadvatc',2,(/ncloc+2*nh-1,nrloc/),kndrtype)

IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (advflux_up(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
   CALL error_alloc('advflux_up',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
   ALLOCATE (usign(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
   CALL error_alloc('usign',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
ENDIF

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (advflux_lw(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
   CALL error_alloc('advflux_lw',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
   ALLOCATE (array2dc1(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
   CALL error_alloc('array2dc1',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
   ALLOCATE (cfl(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
   CALL error_alloc('cfl',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
ENDIF

IF (iopt_adv_3D.EQ.3) THEN
   ALLOCATE (array2dc2(0:ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc1',2,(/ncloc+1,nrloc/),kndrtype)
   ALLOCATE (dflux(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
   CALL error_alloc('dflux',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
   ALLOCATE (omega(0:ncloc,nrloc),STAT=errstat)
   CALL error_alloc('omega',2,(/ncloc+1,nrloc/),kndrtype)
ENDIF

IF (.NOT.dYregX) THEN
   ALLOCATE (fcurvatu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fcurvatu',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (vadvatu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('vadvatu',2,(/ncloc,nrloc/),kndrtype)
   IF (iopt_waves_curr.EQ.1) THEN
      ALLOCATE (vstokesatu(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('vstokesatu',2,(/ncloc,nrloc/),kndrtype)
   ENDIF
ENDIF

!
!2. Initialise arrays
!--------------------
!
!2.1 Mask arrays
!---------------
!

maskatu = nodeatu(1:ncloc,1:nrloc,:).EQ.2
maskatc = nodeatc(1-nh:ncloc+nh-1,1:nrloc).GT.0

!
!2.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   advflux_up = 0.0
ENDIF
IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv_3D.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!2.3 Work space arrays
!---------------------
!

IF (dYregX) THEN
   denom = 2.0*delxatu(1:ncloc,1:nrloc)
ELSE
   denom = 2.0*delyatu(1:ncloc,1:nrloc)*delxatu(1:ncloc,1:nrloc)
   fcurvatu = (delyatc(1:ncloc,1:nrloc)-delyatc(0:ncloc-1,1:nrloc))/denom
ENDIF

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   WHERE (maskatc)
      array2dc1 = delt3d/delxatc(1-nh:ncloc+nh-1,1:nrloc)
   END WHERE
ENDIF

!
!3. Explicit terms
!-----------------
!
!3.1 Advective velocity
!----------------------
!

CALL Uarr_at_C(uvel_old(1-nh:ncloc+nh,1:nrloc,:),uadvatc,1,1,(/1-nh,1,1/),&
            & (/ncloc+nh,nrloc,nz/),1,iarr_uvel_old,info)
IF (iopt_waves_curr.EQ.1) THEN
   uadvatc = uadvatc + ustokesatc(1-nh:ncloc+nh-1,1:nrloc,:)
ENDIF

k_300: DO k=1,nz

!
!
!3.2 Upwind fluxes
!-----------------
!

   IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN

      WHERE (maskatc)
         usign = SIGN(1.0,uadvatc(:,:,k))
         advflux_up = uadvatc(:,:,k)*&
                            & ((1.0+usign)*psiu(1-nh:ncloc+nh-1,1:nrloc,k)&
                            & +(1.0-usign)*psiu(2-nh:ncloc+nh,1:nrloc,k))
      END WHERE

!     ---non-TVD case
      IF (iopt_adv_3D.EQ.1) THEN
         IF (dYregX) THEN
            WHERE (maskatc(0:ncloc,:)) 
               advflux = delzatc(0:ncloc,1:nrloc,k)*advflux_up(0:ncloc,:)
            END WHERE
         ELSE
            WHERE (maskatc(0:ncloc,:)) 
               advflux = delyatc(0:ncloc,1:nrloc)*delzatc(0:ncloc,1:nrloc,k)*&
                       & advflux_up(0:ncloc,:)
            END WHERE
         ENDIF
      ENDIF
      
   ENDIF

!
!3.3 Lax-Wendroff fluxes
!-----------------------
!

   IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN

      WHERE (maskatc)
         cfl = array2dc1*uadvatc(:,:,k)
         advflux_lw = uadvatc(:,:,k)*&
                            & ((1.0+cfl)*psiu(1-nh:ncloc+nh-1,1:nrloc,k)&
                            & +(1.0-cfl)*psiu(2-nh:ncloc+nh,1:nrloc,k))
      END WHERE

!     ---non-TVD case
      IF (iopt_adv_3D.EQ.2) THEN
         IF (dYregX) THEN
            WHERE (maskatc(0:ncloc,:))
               advflux = delzatc(0:ncloc,1:nrloc,k)*advflux_lw(0:ncloc,:)
            END WHERE
         ELSE
            WHERE (maskatc(0:ncloc,:))
               advflux = delyatc(0:ncloc,1:nrloc)*delzatc(0:ncloc,1:nrloc,k)*&
                       & advflux_lw(0:ncloc,:)
            END WHERE
         ENDIF
      ENDIF

   ENDIF

!
!3.4 TVD scheme
!--------------
!

   IF (iopt_adv_3D.EQ.3) THEN

!     ---flux differences
      WHERE (maskatc)
         dflux = advflux_lw - advflux_up
      END WHERE
   
!     ---flux ratios
      WHERE (maskatc(0:ncloc,:)) array2dc2 = 0.0
      WHERE (ABS(dflux(0:ncloc,:)).GT.eps_adv)
         array2dc2 = 0.5*((1.0+usign(0:ncloc,:))*dflux(-1:ncloc-1,:)+&
                        & (1.0-usign(0:ncloc,:))*dflux(1:ncloc+1,:))&
                        & /dflux(0:ncloc,:)
      END WHERE

!     ---TVD limiter
      omega = tvd_limiter(array2dc2,maskatc(0:ncloc,:))

!     ---fluxes
      IF (dYregX) THEN
         WHERE (maskatc(0:ncloc,:))
            advflux = delzatc(0:ncloc,1:nrloc,k)*&
                    & (advflux_up(0:ncloc,:)+omega*dflux(0:ncloc,:))
         END WHERE
      ELSE
         WHERE (maskatc(0:ncloc,:))
            advflux = delyatc(0:ncloc,1:nrloc)*delzatc(0:ncloc,1:nrloc,k)*&
                    & (advflux_up(0:ncloc,:)+omega*dflux(0:ncloc,:))
         END WHERE
      ENDIF

   ENDIF

!
!3.5 Advective term
!------------------
!

   WHERE (maskatu(:,:,k))
      xadvatu3d(:,:,k) = (advflux(1:ncloc,:)-advflux(0:ncloc-1,:))/denom
   END WHERE

!
!3.6 Curvature term
!------------------
!

   IF (.NOT.dYregX) THEN
      CALL Varr_at_U(vvel_old(0:ncloc,1:nrloc+1,k),vadvatu,1,2,(/0,1,k/),&
                           & (/ncloc,nrloc+1,k/),1,iarr_vvel_old,info)
      IF (iopt_waves_curr.EQ.0) THEN
         WHERE (maskatu(:,:,k))
            xadvatu3d(:,:,k) = xadvatu3d(:,:,k) - 2.0*delzatu(1:ncloc,:,k)*&
                                             & vadvatu*vadvatu*fcurvatu
         END WHERE
      ELSE
         CALL Carr_at_U(vstokesatc(0:ncloc,1:nrloc,k),vstokesatu,1,1,&
                     & (/0,1,k/),(/ncloc,nrloc,k/),1,iarr_vstokesatc,info)
         WHERE (maskatu(:,:,k))
            xadvatu3d(:,:,k) = xadvatu3d(:,:,k) - 2.0*delzatu(1:ncloc,:,k)*&
                                        & vadvatu*(vadvatu+vstokesatu)*fcurvatu
         END WHERE
      ENDIF
   ENDIF

!
!3.7 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

   IF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatu(:,:,k))
         xadvatu3d(:,:,k) = alphatu_fld*xadvatu3d(:,:,k)
      END WHERE
   ENDIF

!
!3.8 Apply relaxation scheme
!---------------------------
!

   IF (iopt_obc_advrlx.EQ.1) THEN
      WHERE (maskatu(:,:,k))
         xadvatu3d(:,:,k) = rlxobcatu*xadvatu3d(:,:,k)
      END WHERE
   ENDIF

!
!3.9 Update depth-integrated advective term
!------------------------------------------
!

   IF (vintflag) THEN
      WHERE (maskatu(:,:,k))
         udevint = udevint - xadvatu3d(:,:,k)
      END WHERE
   ENDIF

ENDDO k_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatc,maskatu)
DEALLOCATE (advflux,denom,uadvatc)
IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) DEALLOCATE (advflux_up,usign)
IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) DEALLOCATE (advflux_lw,array2dc1,cfl)
IF (iopt_adv_3D.EQ.3) DEALLOCATE (array2dc2,dflux,omega)
IF (.NOT.dYregX) DEALLOCATE (fcurvatu,vadvatu)
IF (.NOT.dYregX.AND.iopt_waves_curr.EQ.1) DEALLOCATE (vstokesatu)


CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Xadv_at_U_3d

!========================================================================

SUBROUTINE Xadv_at_U_2d(psiu,xadvatu2d,nh,vintflag)
!************************************************************************
!
! *Xadv_at_U_2d* Advection term in X-direction for the 2-D current at U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_U_2d
!
! Module calls - error_alloc, tvd_limiter, Carr_at_U, Uarr_at_C, Varr_at_U
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Carr_at_U, Uarr_at_C, Varr_at_U
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(INOUT) :: vintflag
INTEGER, INTENT(IN) :: nh
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo) :: psiu
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: xadvatu2d

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiu*      REAL    U-node quantity to be advected                    [psiu]
!*xadvatu2d* REAL    Advective term                                  [psiu/s]
!*nh*        INTEGER Halo size of local arrays
!*vintflag*  LOGICAL Update vertically integrated advective term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc, maskatu
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: advflux, advflux_lw, advflux_up,&
                  & array2dc, cfl, denom, dflux, omega, uadvatc, usign, &
                  & vadvatu, vmstokesatu

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at C-nodes          [psiu]
!*advflux_lw* REAL    Lax_Wendroff flux at C-nodes                [psiu]*[m/s]
!*advflux_up* REAL    Upwind flux at C-nodes                      [psiu]*[m/s]
!*omega*      REAL    TVD limiting factor
!*usign*      REAL    Upwind factor
!*uadvatc*    REAL    Advective velocity at C-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xadv_at_U_2d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatc(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc+2*nh-1,nrloc/),kndlog)

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)

ALLOCATE (advflux(0:ncloc,nrloc),STAT=errstat)
CALL error_alloc('advflux',2,(/ncloc+1,nrloc/),kndrtype)

ALLOCATE (uadvatc(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
CALL error_alloc('uadvatc',2,(/ncloc+2*nh-1,nrloc/),kndrtype)

IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN
   ALLOCATE (advflux_up(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
   CALL error_alloc('advflux_up',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
   ALLOCATE (usign(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
   CALL error_alloc('usign',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
ENDIF

IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN
   ALLOCATE (advflux_lw(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
   CALL error_alloc('advflux_lw',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
   ALLOCATE (cfl(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
   CALL error_alloc('cfl',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
ENDIF

IF (iopt_adv_2D.EQ.3) THEN
   ALLOCATE (array2dc(0:ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc',2,(/ncloc+1,nrloc/),kndrtype)
   ALLOCATE (dflux(1-nh:ncloc+nh-1,nrloc),STAT=errstat)
   CALL error_alloc('dflux',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
   ALLOCATE (omega(0:ncloc,nrloc),STAT=errstat)
   CALL error_alloc('omega',2,(/ncloc+1,nrloc/),kndrtype)
ENDIF

IF (.NOT.dYregX) THEN
   ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (vadvatu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('vadvatu',2,(/ncloc,nrloc/),kndrtype)
   IF (iopt_waves_curr.EQ.1) THEN
      ALLOCATE (vmstokesatu(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('vmstokesatu',2,(/ncloc,nrloc/),kndrtype)
   ENDIF
ENDIF

!
!2. Initialise arrays
!--------------------
!
!2.1 Mask arrays
!---------------
!

maskatu = node2du(1:ncloc,1:nrloc).EQ.2
maskatc = nodeatc(1-nh:ncloc+nh-1,1:nrloc).GT.0

!
!2.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN
   advflux_up = 0.0
ENDIF
IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv_2D.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!2.3 Work space array
!--------------------
!

IF (.NOT.dYregX) THEN
   WHERE (maskatu)
      denom = delxatu(1:ncloc,1:nrloc)*delyatu(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!2.4 Advecting velocity
!----------------------
!

CALL Uarr_at_C(umvel_old(1-nh:ncloc+nh,1:nrloc),uadvatc,1,1,(/1-nh,1,nz/),&
            & (/ncloc+nh,nrloc,nz/),1,iarr_umvel_old,.TRUE.)
IF (iopt_waves_curr.EQ.1) THEN
   uadvatc = uadvatc + umstokesatc(1-nh:ncloc+nh-1,1:nrloc)
ENDIF

!
!3. Explicit terms
!-----------------
!
!3.1 Upwind fluxes
!-----------------
!

IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN

   WHERE (maskatc)
      usign = SIGN(1.0,uadvatc)
      advflux_up = uadvatc*((1.0+usign)*psiu(1-nh:ncloc+nh-1,1:nrloc)&
                         & +(1.0-usign)*psiu(2-nh:ncloc+nh,1:nrloc))
   END WHERE

!  ---non-TVD case
   IF (iopt_adv_2D.EQ.1) THEN
      IF (dYregX) THEN
         WHERE (maskatc(0:ncloc,:))
            advflux = advflux_up(0:ncloc,:)
         END WHERE
      ELSE
         WHERE (maskatc(0:ncloc,:))
            advflux = delyatc(0:ncloc,1:nrloc)*advflux_up(0:ncloc,:)
         END WHERE
      ENDIF
   ENDIF

ENDIF

!
!3.2 Lax-Wendroff fluxes
!-----------------------
!

IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN

   WHERE (maskatc)
      cfl = delt2d*uadvatc/delxatc(1-nh:ncloc+nh-1,1:nrloc)
      advflux_lw = uadvatc*((1.0+cfl)*psiu(1-nh:ncloc+nh-1,1:nrloc)&
                         & +(1.0-cfl)*psiu(2-nh:ncloc+nh,1:nrloc))
   END WHERE

!  ---non-TVD case
   IF (iopt_adv_2D.EQ.2) THEN
      IF (dYregX) THEN
         WHERE (maskatc(0:ncloc,:))
            advflux = advflux_lw(0:ncloc,:)
         END WHERE
      ELSE
         WHERE (maskatc(0:ncloc,:))
            advflux = delyatc(0:ncloc,1:nrloc)*advflux_lw(0:ncloc,:)
         END WHERE
      ENDIF
   ENDIF

ENDIF

!
!3.3 TVD scheme
!--------------
!

IF (iopt_adv_2D.EQ.3) THEN

!  ---flux differences
   WHERE (maskatc)
      dflux = advflux_lw - advflux_up
   END WHERE
   
!  ---flux ratios
   WHERE (maskatc(0:ncloc,:)) array2dc = 0.0
   WHERE (ABS(dflux(0:ncloc,:)).GT.eps_adv)
      array2dc = 0.5*((1.0+usign(0:ncloc,:))*dflux(-1:ncloc-1,:)+&
                    & (1.0-usign(0:ncloc,:))*dflux(1:ncloc+1,:))&
                    & /dflux(0:ncloc,:)
   END WHERE

!  ---TVD limiter
   omega = tvd_limiter(array2dc,maskatc(0:ncloc,:))

!  ---fluxes
   IF (dYregX) THEN
      WHERE (maskatc(0:ncloc,:))
         advflux = advflux_up(0:ncloc,:)+ omega*dflux(0:ncloc,:)
      END WHERE
   ELSE
      WHERE (maskatc(0:ncloc,:))
         advflux = delyatc(0:ncloc,1:nrloc)*(advflux_up(0:ncloc,:)+&
                                           & omega*dflux(0:ncloc,:))
      END WHERE
   ENDIF

ENDIF

!
!3.4 Advective term
!------------------
!
!---regular grid
IF (dYregX) THEN
   WHERE (maskatu)
      xadvatu2d = 0.5*(advflux(1:ncloc,:)-advflux(0:ncloc-1,:))&
                & /delxatu(1:ncloc,1:nrloc)
   END WHERE

!---irregular grid
ELSE
   WHERE (maskatu)
      xadvatu2d = 0.5*(advflux(1:ncloc,:)-advflux(0:ncloc-1,:))/denom
   END WHERE
ENDIF

!
!3.5 Curvature term
!-------------------
!

IF (.NOT.dYregX) THEN
   CALL Varr_at_U(vmvel_old(0:ncloc,1:nrloc+1),vadvatu,1,2,(/0,1,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vmvel_old,.TRUE.)
   IF (iopt_waves_curr.EQ.0) THEN
      WHERE (maskatu)
         xadvatu2d = xadvatu2d - deptotatu(1:ncloc,1:nrloc)*vadvatu*vadvatu*&
                  & (delyatc(1:ncloc,1:nrloc)-delyatc(0:ncloc-1,1:nrloc))/denom
      END WHERE
   ELSE
      CALL Carr_at_U(vmstokesatc(0:ncloc,1:nrloc),vmstokesatu,1,1,&
                  & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_vmstokesatc,.TRUE.)
      WHERE (maskatu)
         xadvatu2d = xadvatu2d - deptotatu(1:ncloc,1:nrloc)*vadvatu*&
                   & (vadvatu+vmstokesatu)*&
                   & (delyatc(1:ncloc,1:nrloc)-delyatc(0:ncloc-1,1:nrloc))/denom
      END WHERE
   ENDIF
ENDIF

!
!3.6 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN
   WHERE (maskatu)
      xadvatu2d = alphatu_fld*xadvatu2d
   END WHERE
ENDIF

!
!3.7 Apply relaxation scheme
!---------------------------
!

IF (iopt_obc_advrlx.EQ.1) THEN
   WHERE (maskatu)
      xadvatu2d = rlxobcatu*xadvatu2d
   END WHERE
ENDIF

!
!3.8 Update depth-integrated advective term
!------------------------------------------
!

IF (vintflag) THEN
   WHERE (maskatu)
      udevint = udevint + xadvatu2d
   END WHERE
   vintflag = .FALSE.
ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatc,maskatu)
DEALLOCATE (advflux,uadvatc)
IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) DEALLOCATE (advflux_up,usign)
IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) DEALLOCATE (advflux_lw,cfl)
IF (iopt_adv_2D.EQ.3) DEALLOCATE (array2dc,dflux,omega)
IF (.NOT.dYregX) DEALLOCATE (denom,vadvatu)
IF (.NOT.dYregX.AND.iopt_waves_curr.EQ.1) DEALLOCATE (vmstokesatu)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Xadv_at_U_2d

!========================================================================

SUBROUTINE Yadv_at_U_3d(psiu,yadvatu3d,nh,vintflag)
!************************************************************************
!
! *Yadv_at_U_3d* Advection term in Y-direction for the 3-D current at U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_U_3d
!
! Module calls - error_alloc, tvd_limiter, Varr_at_U, Varr_at_UV
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Varr_at_U, Varr_at_UV
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: vintflag
INTEGER, INTENT(IN) :: nh
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiu
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,nz) :: yadvatu3d

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiu*      REAL    U-node quantity to be advected                     [psiu]
!*yadvatu3d* REAL    Advective term (times delz)                    [m*psiu/s]
!*nh*        INTEGER Halo size of local arrays
!*vintflag*  LOGICAL Update vertically integrated advective term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: info = .FALSE.
INTEGER :: i, j, jj, jjloc, k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2duv
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu, maskatuv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2duv, denom, fcurvatu
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advflux, advflux_lw, advflux_up, &
              & array3duv1, array3duv2, cfl, dflux, omega, vadvatu, vadvatuv, &
              & vsign, vstokesatuv

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at UV-nodes         [psiu]
!*advflux_lw* REAL    Lax_Wendroff flux at UV- nodes              [psiu]*[m/s]
!*advflux_up* REAL    Upwind flux at UV nodes                     [psiu]*[m/s]
!*omega*      REAL    TVD limiting factor
!*vsign*      REAL    Upwind factor
!*vadvatu*    REAL    Advective velocity at U-nodes                      [m/s]
!*vadvatuv*   REAL    Advective velocity at UV-nodes                     [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Yadv_at_U_3d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (advflux(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('advflux',3,(/ncloc,nrloc+1,nz/),kndrtype)

ALLOCATE (array2duv(ncloc,2-nh:nrloc+nh),STAT=errstat)
CALL error_alloc('array2duv',2,(/ncloc,nrloc+2*nh-1/),kndrtype)

ALLOCATE (array3duv1(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('array3duv1',3,(/ncloc,nrloc+1,nz/),kndrtype)

ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE (maskatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc,nrloc,nz/),kndlog)

ALLOCATE (maskatuv(ncloc,2-nh:nrloc+nh,nz),STAT=errstat)
CALL error_alloc('maskatuv',3,(/ncloc,nrloc+2*nh-1,nz/),kndlog)

ALLOCATE (mask2duv(ncloc,2-nh:nrloc+nh),STAT=errstat)
CALL error_alloc('mask2duv',2,(/ncloc,nrloc+2*nh-1/),kndlog)

ALLOCATE (vadvatuv(ncloc,2-nh:nrloc+nh,nz),STAT=errstat)
CALL error_alloc('vadvatuv',3,(/ncloc,nrloc+2*nh-1,nz/),kndrtype)

IF (iopt_waves_curr.EQ.1) THEN
   ALLOCATE (vstokesatuv(ncloc,2-nh:nrloc+nh,nz),STAT=errstat)
   CALL error_alloc('vstokesatuv',3,(/ncloc,nrloc+2*nh-1,nz/),kndrtype)
ENDIF

ALLOCATE (vsign(ncloc,2-nh:nrloc+nh,nz),STAT=errstat)
CALL error_alloc('vsign',3,(/ncloc,nrloc+2*nh-1,nz/),kndrtype)

IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (advflux_up(ncloc,2-nh:nrloc+nh,nz),STAT=errstat)
   CALL error_alloc('advflux_up',3,(/ncloc,nrloc+2*nh-1,nz/),kndrtype)
ENDIF

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (advflux_lw(ncloc,2-nh:nrloc+nh,nz),STAT=errstat)
   CALL error_alloc('advflux_lw',3,(/ncloc,nrloc+2*nh-1,nz/),kndrtype)
   ALLOCATE (cfl(ncloc,2-nh:nrloc+nh,nz),STAT=errstat)
   CALL error_alloc('cfl',3,(/ncloc,nrloc+2*nh-1,nz/),kndrtype)
ENDIF

IF (iopt_adv_3D.EQ.3) THEN
   ALLOCATE (array3duv2(ncloc,nrloc+1,nz),STAT=errstat)
   CALL error_alloc('array3duv2',3,(/ncloc,nrloc+1,nz/),kndrtype)
   ALLOCATE (dflux(ncloc,2-nh:nrloc+nh,nz),STAT=errstat)
   CALL error_alloc('dflux',3,(/ncloc,nrloc+2*nh-1,nz/),kndrtype)
   ALLOCATE (omega(ncloc,nrloc+1,nz),STAT=errstat)
   CALL error_alloc('omega',3,(/ncloc,nrloc+1,nz/),kndrtype)
ENDIF

IF (.NOT.dXregY) THEN
   ALLOCATE (fcurvatu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fcurvatu',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (vadvatu(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('vadvatu',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

!
!2. Initialise arrays
!--------------------
!
!2.1 Mask arrays
!---------------
!

maskatu = nodeatu(1:ncloc,1:nrloc,:).EQ.2
maskatuv = nodeatuv(1:ncloc,2-nh:nrloc+nh,:).GT.0
mask2duv = node2duv(1:ncloc,2-nh:nrloc+nh).GT.0

!
!2.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   advflux_up = 0.0
ENDIF
IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv_3D.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!2.3 Metric factor in advective flux
!-----------------------------------
!

IF (dXregY) THEN
   k_231: DO k=1,nz
      WHERE (maskatuv(:,1:nrloc+1,k))
         array3duv1(:,:,k) = delzatuv(1:ncloc,:,k)
      END WHERE
   ENDDO k_231
ELSE
   k_232: DO k=1,nz
      WHERE (maskatuv(:,1:nrloc+1,k))
         array3duv1(:,:,k) = delxatuv(1:ncloc,1:nrloc+1)*delzatuv(1:ncloc,:,k)
      END WHERE
   ENDDO k_232
ENDIF

!
!2.4 Advective velocity
!----------------------
!

CALL Varr_at_UV(vvel_old(0:ncloc,2-nh:nrloc+nh,:),vadvatuv,1,1,&
             & (/0,2-nh,1/),(/ncloc,nrloc+nh,nz/),1,iarr_vvel_old,info)
IF (iopt_waves_curr.EQ.1) THEN
   CALL Varr_at_UV(vstokesatv(0:ncloc,2-nh:nrloc+nh,:),vstokesatuv,1,1,&
        & (/0,2-nh,1/),(/ncloc,nrloc+nh,nz/),1,iarr_vstokesatv,info)
   vadvatuv = vadvatuv + vstokesatuv
ENDIF

!
!2.5 Sign of velocity
!--------------------
!

WHERE (maskatuv)
   vsign = SIGN(1.0,vadvatuv)
END WHERE

!
!2.6 CFL number
!--------------
!

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   WHERE (mask2duv)
      array2duv = delt3d/delyatuv(1:ncloc,2-nh:nrloc+nh)
   END WHERE
   k_260: DO k=1,nz
      WHERE (maskatuv(:,:,k))
         cfl(:,:,k) = array2duv*vadvatuv(:,:,k)
      END WHERE
   ENDDO k_260
ENDIF

!
!2.7 Work space arrays
!---------------------
!

IF (dXregY) THEN
   denom = 2.0*delyatu(1:ncloc,1:nrloc)
ELSE
   denom = 2.0*delxatu(1:ncloc,1:nrloc)*delyatu(1:ncloc,1:nrloc)
   fcurvatu = (delxatuv(1:ncloc,2:nrloc+1)-delxatuv(1:ncloc,1:nrloc))/denom
ENDIF

!
!3. Explicit terms
!-----------------
!

maskatuv = nodeatuv(1:ncloc,2-nh:nrloc+nh,:).EQ.1.OR.&
         & nodeatuv(1:ncloc,2-nh:nrloc+nh,:).EQ.2

!
!3.1 Upwind fluxes
!-----------------
!

IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN

   k_310: DO k=1,nz
      WHERE (maskatuv(:,:,k))
         advflux_up(:,:,k) = vadvatuv(:,:,k)*&
                          & ((1.0+vsign(:,:,k))*psiu(1:ncloc,1-nh:nrloc+nh-1,k)&
                          & +(1.0-vsign(:,:,k))*psiu(1:ncloc,2-nh:nrloc+nh,k))
      END WHERE
   ENDDO k_310

!  ---non-TVD case
   IF (iopt_adv_3D.EQ.1) THEN
      WHERE (maskatuv(:,1:nrloc+1,:))
         advflux = array3duv1*advflux_up(:,1:nrloc+1,:)
      END WHERE
   ENDIF

ENDIF

!
!3.2 Lax-Wendroff fluxes
!-----------------------
!

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN

   k_320: DO k=1,nz
      WHERE (maskatuv(:,:,k))
         advflux_lw(:,:,k) = vadvatuv(:,:,k)*&
                           & ((1.0+cfl(:,:,k))*psiu(1:ncloc,1-nh:nrloc+nh-1,k)&
                           & +(1.0-cfl(:,:,k))*psiu(1:ncloc,2-nh:nrloc+nh,k))
      END WHERE
   ENDDO k_320

!  ---non-TVD case
   IF (iopt_adv_3D.EQ.2) THEN
      WHERE (maskatuv(:,1:nrloc+1,:))
         advflux = array3duv1*advflux_lw(:,1:nrloc+1,:)
      END WHERE
   ENDIF

ENDIF

!
!3.3 TVD scheme
!--------------
!

IF (iopt_adv_3D.EQ.3) THEN

!  ---flux differences
   WHERE (maskatuv)
      dflux = advflux_lw - advflux_up
   END WHERE

!  ---flux ratios
   WHERE (maskatuv(:,1:nrloc+1,:)) array3duv2 = 0.0
   k_330: DO k=1,nz
      WHERE (ABS(dflux(:,1:nrloc+1,k)).GT.eps_adv)
         array3duv2(:,:,k) = 0.5*(&
            & (1.0+vsign(:,1:nrloc+1,k))*dflux(:,0:nrloc,k)+&
            & (1.0-vsign(:,1:nrloc+1,k))*dflux(:,2:nrloc+2,k))&
            & /dflux(:,1:nrloc+1,k)
      END WHERE
   ENDDO k_330

!  ---TVD limiter
   omega = tvd_limiter(array3duv2,maskatuv(:,1:nrloc+1,:)) 

!  ---fluxes
   WHERE (maskatuv(:,1:nrloc+1,:))
      advflux = array3duv1*(advflux_up(:,1:nrloc+1,:)+&
                          & omega*dflux(:,1:nrloc+1,:))
   END WHERE
   
ENDIF

!
!3.4 Open boundary fluxes
!------------------------
!

jjloc_340: DO jjloc=1,nobyloc_ext
   i = iobyloc(jjloc); j = jobyloc(jjloc)
   jj = indexoby(jjloc)

   SELECT CASE (itypoby(jj))

!
!3.4.1 Zero gradient condition
!-----------------------------
!

      CASE (0)
         IF (soutoby(jj)) THEN
            IF (j.LT.nrloc+1) advflux(i,j,:) = advflux(i,j+1,:)
         ELSE
            IF (j.GT.1) advflux(i,j,:) = advflux(i,j-1,:)
         ENDIF

!
!3.4.2 Upwind scheme
!-------------------
!

      CASE (1)
         IF (soutoby(jj)) THEN
            k_3421: DO k=1,nz
               IF (vsign(i,j,k).LT.0.0) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*vadvatuv(i,j,k)*&
                                 & psiu(i,j,k)
               ELSEIF (nodeatu(i,j-1,k).GT.2.OR.nodeatv(i-1,j,k).LE.1.OR.&
                     & nodeatv(i,j,k).LE.1) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*vadvatuv(i,j,k)*&
                                 & psiu(i,j-1,k)
               ELSE
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*vadvatuv(i,j,k)*&
                                 & psiu(i,j,k)
               ENDIF
            ENDDO k_3421
         ELSE
            k_3422: DO k=1,nz
               IF (vsign(i,j,k).GT.0.0) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*vadvatuv(i,j,k)*&
                                 & psiu(i,j-1,k)
               ELSEIF (nodeatu(i,j,k).GT.2.OR.nodeatv(i-1,j,k).LE.1.OR.&
                     & nodeatv(i,j,k).LE.1) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*vadvatuv(i,j,k)*&
                                 & psiu(i,j,k)
               ELSE
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*vadvatuv(i,j,k)*&
                                 & psiu(i,j-1,k)
               ENDIF
            ENDDO k_3422
         ENDIF

!
!3.4.3 Tangential condition with external values
!-----------------------------------------------
!

      CASE (2)
         IF (soutoby(jj)) THEN
            k_3431: DO k=1,nz
               IF (vsign(i,j,k).LT.0.0) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*vadvatuv(i,j,k)*&
                                 & psiu(i,j,k)
               ELSE
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*vadvatuv(i,j,k)*&
                                 & uveloby(jj,k)  
               ENDIF
            ENDDO k_3431
         ELSE
            k_3432: DO k=1,nz
               IF (vsign(i,j,k).GT.0.0) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*vadvatuv(i,j,k)*&
                                 & psiu(i,j-1,k)
               ELSE
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*vadvatuv(i,j,k)*&
                                 & uveloby(jj,k)  
               ENDIF
            ENDDO k_3432
         ENDIF

   END SELECT

ENDDO jjloc_340

!
!3.5 Advective term
!------------------
!

k_350: DO k=1,nz
   WHERE (maskatu(:,:,k))
      yadvatu3d(:,:,k) = (advflux(:,2:nrloc+1,k)-advflux(:,1:nrloc,k))/denom
   END WHERE
ENDDO k_350

!
!3.6 Curvature term
!------------------
!

IF (.NOT.dXregY) THEN
   CALL Varr_at_U(vvel_old(0:ncloc,1:nrloc+1,:),vadvatu,1,2,(/0,1,1/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vvel_old,info)
   IF (iopt_waves_curr.EQ.0) THEN
      k_361: DO k=1,nz
         WHERE (maskatu(:,:,k))
            yadvatu3d(:,:,k) = yadvatu3d(:,:,k) + 2.0*delzatu(1:ncloc,:,k)*&
                             & vadvatu(:,:,k)*psiu(1:ncloc,1:nrloc,k)*fcurvatu
         END WHERE
      ENDDO k_361
   ELSE
      k_362: DO k=1,nz
         WHERE (maskatu(:,:,k))
            yadvatu3d(:,:,k) = yadvatu3d(:,:,k) + 2.0*delzatu(1:ncloc,:,k)*&
                             & vadvatu(:,:,k)*&
                             & (psiu(1:ncloc,1:nrloc,k)+&
                             & ustokesatu(1:ncloc,1:nrloc,k))*fcurvatu
         END WHERE
      ENDDO k_362
   ENDIF
ENDIF

!
!3.7 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN
   k_370: DO k=1,nz
      WHERE (maskatu(:,:,k))
         yadvatu3d(:,:,k) = alphatu_fld*yadvatu3d(:,:,k)
      END WHERE
   ENDDO k_370
ENDIF

!
!3.8 Apply relaxation scheme
!---------------------------
!

IF (iopt_obc_advrlx.EQ.1) THEN
   k_380: DO k=1,nz
      WHERE (maskatu(:,:,k))
         yadvatu3d(:,:,k) = rlxobcatu*yadvatu3d(:,:,k)
      END WHERE
   ENDDO k_380
ENDIF

!
!3.9 Update depth-integrated advective term
!------------------------------------------
!

IF (vintflag) THEN
   k_390: DO k=1,nz
      WHERE (maskatu(:,:,k))
         udevint = udevint - yadvatu3d(:,:,k)
      END WHERE
   ENDDO k_390
ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatu,maskatuv,mask2duv)
DEALLOCATE (advflux,array2duv,array3duv1,denom,vadvatuv,vsign)
IF (iopt_waves_curr.EQ.1) DEALLOCATE (vstokesatuv)   
IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) DEALLOCATE (advflux_up)
IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) DEALLOCATE (advflux_lw,cfl)
IF (iopt_adv_3D.EQ.3) DEALLOCATE (array3duv2,dflux,omega)
IF (.NOT.dXregY) DEALLOCATE (fcurvatu,vadvatu)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Yadv_at_U_3d

!========================================================================

SUBROUTINE Yadv_at_U_2d(psiu,yadvatu2d,nh,vintflag)
!************************************************************************
!
! *Yadv_at_U_2d* Advection term in Y-direction for the 2-D current at U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_U_2d
!
! Module calls - error_alloc, tvd_limiter, Varr_at_U, Varr_at_UV
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Varr_at_U, Varr_at_UV
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(INOUT) :: vintflag
INTEGER, INTENT(IN) :: nh
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo) :: psiu
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: yadvatu2d

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiu*      REAL    U-node quantity to be advected                     [psiu]
!*yadvatu2d* REAL    Advective term                                   [psiu/s]
!*nh*        INTEGER Halo size of local arrays
!*vintflag*  LOGICAL Update vertically integrated advective term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, jj, jjloc, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatuv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: advflux, advflux_lw, advflux_up,&
       & array2duv, cfl, dflux, denom, omega, vadvatu, vadvatuv, &
       & vmstokesatuv, vsign

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at UV-nodes         [psiu]
!*advflux_lw* REAL    Lax_Wendroff flux at UV-nodes               [psiu]*[m/s]
!*advflux_up* REAL    Upwind flux at UV-nodes                     [psiu]*[m/s]
!*omega*      REAL    TVD limiting factor
!*vsign*      REAL    Upwind factor
!*vadvatu*    REAL    Advective velocity at U-nodes                      [m/s]
!*vadvatuv*   REAL   Advective velocity at UV-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Yadv_at_U_2d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (advflux(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('advflux',2,(/ncloc,nrloc+1/),kndrtype)

ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)

ALLOCATE (maskatuv(ncloc,2-nh:nrloc+nh),STAT=errstat)
CALL error_alloc('maskatuv',2,(/ncloc,nrloc+2*nh-1/),kndlog)

ALLOCATE (vadvatuv(ncloc,2-nh:nrloc+nh),STAT=errstat)
CALL error_alloc('vadvatuv',2,(/ncloc,nrloc+2*nh-1/),kndrtype)

IF (iopt_waves_curr.EQ.1) THEN
   ALLOCATE (vmstokesatuv(ncloc,2-nh:nrloc+nh),STAT=errstat)
   CALL error_alloc('vmstokesatuv',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
ENDIF

ALLOCATE (vsign(ncloc,2-nh:nrloc+nh),STAT=errstat)
CALL error_alloc('vsign',2,(/ncloc,nrloc+2*nh-1/),kndrtype)

IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN
   ALLOCATE (advflux_up(ncloc,2-nh:nrloc+nh),STAT=errstat)
   CALL error_alloc('advflux_up',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
ENDIF

IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN
   ALLOCATE (advflux_lw(ncloc,2-nh:nrloc+nh),STAT=errstat)
   CALL error_alloc('advflux_lw',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
   ALLOCATE (cfl(ncloc,2-nh:nrloc+nh),STAT=errstat)
   CALL error_alloc('cfl',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
ENDIF

IF (iopt_adv_2D.EQ.3) THEN
   ALLOCATE (array2duv(ncloc,nrloc+1),STAT=errstat)
   CALL error_alloc('array2duv',2,(/ncloc,nrloc+1/),kndrtype)
   ALLOCATE (dflux(ncloc,2-nh:nrloc+nh),STAT=errstat)
   CALL error_alloc('dflux',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
   ALLOCATE (omega(ncloc,nrloc+1),STAT=errstat)
   CALL error_alloc('omega',2,(/ncloc,nrloc+1/),kndrtype)
ENDIF

IF (.NOT.dXregY) THEN
   ALLOCATE (vadvatu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('vadvatu',2,(/ncloc,nrloc/),kndrtype)
ENDIF

!
!2. Initialise arrays
!--------------------
!
!2.1 Mask arrays
!---------------
!

maskatu = node2du(1:ncloc,1:nrloc).EQ.2
maskatuv = node2duv(1:ncloc,2-nh:nrloc+nh).GT.0

!
!2.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN
   advflux_up = 0.0
ENDIF
IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv_2D.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!2.3 Advective velocities
!------------------------
!

CALL Varr_at_UV(vmvel_old(0:ncloc,2-nh:nrloc+nh),vadvatuv,1,1,(/0,2-nh,nz/),&
             & (/ncloc,nrloc+nh,nz/),1,iarr_vmvel_old,.TRUE.)
IF (iopt_waves_curr.EQ.1) THEN
   CALL Varr_at_UV(vmstokesatv(0:ncloc,2-nh:nrloc+nh),vmstokesatuv,1,1,&
        & (/0,2-nh,nz/),(/ncloc,nrloc+nh,nz/),1,iarr_vmstokesatv,.TRUE.)
   vadvatuv = vadvatuv + vmstokesatuv
ENDIF

!
!2.4 Sign of velocity
!--------------------
!

WHERE (maskatuv)
   vsign = SIGN(1.0,vadvatuv)
END WHERE

!
!2.5 Work space arrays
!---------------------
!

IF (dXregY) THEN
   WHERE (maskatu)
      denom = delyatu(1:ncloc,1:nrloc)
   END WHERE
ELSE
   WHERE (maskatu)
      denom = delxatu(1:ncloc,1:nrloc)*delyatu(1:ncloc,1:nrloc)
   END WHERE

ENDIF

!
!3. Explicit terms
!-----------------
!

maskatuv = node2duv(1:ncloc,2-nh:nrloc+nh).EQ.1.OR.&
         & node2duv(1:ncloc,2-nh:nrloc+nh).EQ.2

!
!3.1 Upwind fluxes
!-----------------
!

IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN

   WHERE (maskatuv)
      advflux_up = vadvatuv*((1.0+vsign)*psiu(1:ncloc,1-nh:nrloc+nh-1)&
                          & +(1.0-vsign)*psiu(1:ncloc,2-nh:nrloc+nh))
   END WHERE

!  ---non-TVD case
   IF (iopt_adv_2D.EQ.1) THEN
      IF (dXregY) THEN
         WHERE (maskatuv(:,1:nrloc+1))
            advflux = advflux_up(:,1:nrloc+1)
         END WHERE
      ELSE
         WHERE (maskatuv(:,1:nrloc+1))
            advflux = delxatuv(1:ncloc,1:nrloc+1)*advflux_up(:,1:nrloc+1)
         END WHERE
      ENDIF
   ENDIF

ENDIF

!
!3.2 Lax-Wendroff fluxes
!-----------------------
!

IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN

   WHERE (maskatuv)
      cfl = delt2d*vadvatuv/delyatuv(1:ncloc,2-nh:nrloc+nh)
      advflux_lw = vadvatuv*((1.0+cfl)*psiu(1:ncloc,1-nh:nrloc+nh-1)&
                          & +(1.0-cfl)*psiu(1:ncloc,2-nh:nrloc+nh))
   END WHERE

!  ---non-TVD case
   IF (iopt_adv_2D.EQ.2) THEN
      IF (dXregY) THEN
         WHERE (maskatuv(:,1:nrloc+1))
            advflux = advflux_lw(:,1:nrloc+1)
         END WHERE
      ELSE
         WHERE (maskatuv(:,1:nrloc+1))
            advflux = delxatuv(1:ncloc,1:nrloc+1)*advflux_lw(:,1:nrloc+1)
         END WHERE
      ENDIF
   ENDIF

ENDIF

!
!3.3 TVD scheme
!--------------
!

IF (iopt_adv_2D.EQ.3) THEN

!  ---flux differences
   WHERE (maskatuv)
      dflux = advflux_lw - advflux_up
   END WHERE

!  ---flux ratios
   WHERE (maskatuv(:,1:nrloc+1)) array2duv = 0.0
   WHERE (ABS(dflux(:,1:nrloc+1)).GT.eps_adv)
      array2duv = 0.5*(&
                  & (1.0+vsign(:,1:nrloc+1))*dflux(:,0:nrloc)+&
                  & (1.0-vsign(:,1:nrloc+1))*dflux(:,2:nrloc+2))&
                  & /dflux(:,1:nrloc+1)
   END WHERE

!  ---TVD limiter
   omega = tvd_limiter(array2duv,maskatuv(:,1:nrloc+1))

!  ---fluxes
   IF (dXregY) THEN
      WHERE (maskatuv(:,1:nrloc+1))
         advflux = advflux_up(:,1:nrloc+1)+omega*dflux(:,1:nrloc+1)
      END WHERE
   ELSE
      WHERE (maskatuv(:,1:nrloc+1))
         advflux = delxatuv(1:ncloc,1:nrloc+1)*&
                 & (advflux_up(:,1:nrloc+1)+omega*dflux(:,1:nrloc+1))
      END WHERE
   ENDIF

ENDIF

!
!3.4 Open boundary fluxes
!------------------------
!

jjloc_341: DO jjloc=1,nobyloc_ext
   i = iobyloc(jjloc); j = jobyloc(jjloc)
   jj = indexoby(jjloc)

   SELECT CASE (ityp2doby(jj))

!
!3.4.1 Zero gradient condition
!-----------------------------
!

      CASE (0)
         IF (soutoby(jj)) THEN
            IF (j.LT.nrloc+1) advflux(i,j) = advflux(i,j+1)
         ELSE
            IF (j.GT.1) advflux(i,j) = advflux(i,j-1)
         ENDIF

!
!3.4.2 Upwind scheme
!-------------------
!

      CASE (1)
         IF (soutoby(jj)) THEN
            IF (vsign(i,j).LT.0.0) THEN
               advflux(i,j) = 2.0*vadvatuv(i,j)*psiu(i,j)
            ELSEIF (node2du(i,j-1).GT.2.OR.node2dv(i-1,j).LE.1.OR.&
                  & node2dv(i,j).LE.1) THEN
               advflux(i,j) = 2.0*vadvatuv(i,j)*psiu(i,j-1)
            ELSE
               advflux(i,j) = 2.0*vadvatuv(i,j)*psiu(i,j)
            ENDIF
         ELSE
            IF (vsign(i,j).GT.0.0) THEN
               advflux(i,j) = 2.0*vadvatuv(i,j)*psiu(i,j-1)
            ELSEIF (node2du(i,j).GT.2.OR.node2dv(i-1,j).LE.1.OR.&
                  & node2dv(i,j).LE.1) THEN
               advflux(i,j) = 2.0*vadvatuv(i,j)*psiu(i,j)
            ELSE
               advflux(i,j) = 2.0*vadvatuv(i,j)*psiu(i,j-1)
            ENDIF
         ENDIF
         IF (.NOT.dXregY) advflux(i,j) = delxatuv(i,j)*advflux(i,j)

!
!3.4.3 Tangential condition with external values
!-----------------------------------------------
!

      CASE (2)
         IF (soutoby(jj)) THEN
            IF (vsign(i,j).LT.0.0) THEN
               advflux(i,j) = 2.0*vadvatuv(i,j)*psiu(i,j)
            ELSE
               advflux(i,j) = 2.0*vadvatuv(i,j)*udatoby(jj,2)
            ENDIF
         ELSE
            IF (vsign(i,j).GT.0.0) THEN
               advflux(i,j) = 2.0*vadvatuv(i,j)*psiu(i,j-1)
            ELSE
               advflux(i,j) = 2.0*vadvatuv(i,j)*udatoby(jj,2)
            ENDIF
         ENDIF
         IF (.NOT.dXregY) advflux(i,j) = delxatuv(i,j)*advflux(i,j)

   END SELECT

ENDDO jjloc_341

!
!3.5 Advective term
!------------------
!

WHERE (maskatu)
   yadvatu2d = 0.5*(advflux(:,2:nrloc+1)-advflux(:,1:nrloc))/denom
END WHERE

!
!3.6 Curvature term
!-------------------
!

IF (.NOT.dXregY) THEN
   CALL Varr_at_U(vmvel_old(0:ncloc,1:nrloc+1),vadvatu,1,2,(/0,1,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vmvel_old,.TRUE.)
   IF (iopt_waves_curr.EQ.0) THEN
      WHERE (maskatu)
         yadvatu2d = yadvatu2d + vadvatu*psiu(1:ncloc,1:nrloc)*&
                 & (delxatuv(1:ncloc,2:nrloc+1)-delxatuv(1:ncloc,1:nrloc))/denom
      END WHERE
   ELSE
      WHERE (maskatu)
         yadvatu2d = yadvatu2d + vadvatu*&
                 & (psiu(1:ncloc,1:nrloc)+umstokesatu(1:ncloc,1:nrloc))*&
                 & (delxatuv(1:ncloc,2:nrloc+1)-delxatuv(1:ncloc,1:nrloc))/denom
      END WHERE
   ENDIF
ENDIF

!
!3.7 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN
   WHERE (maskatu)
      yadvatu2d = alphatu_fld*yadvatu2d
   END WHERE
ENDIF

!
!3.8 Apply relaxation scheme
!---------------------------
!

IF (iopt_obc_advrlx.EQ.1) THEN
   WHERE (maskatu)
      yadvatu2d = rlxobcatu*yadvatu2d
   END WHERE
ENDIF

!
!3.9 Update depth-integrated advective term
!------------------------------------------
!

IF (vintflag) THEN
   WHERE (maskatu)
      udevint = udevint + yadvatu2d
   END WHERE
   vintflag = .FALSE.
ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatu,maskatuv)
DEALLOCATE (advflux,denom,vadvatuv,vsign)
IF (iopt_waves_curr.EQ.1) DEALLOCATE (vmstokesatuv)
IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) DEALLOCATE (advflux_up)
IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) DEALLOCATE (advflux_lw,cfl)
IF (iopt_adv_2D.EQ.3) DEALLOCATE (array2duv,dflux,omega)
IF (.NOT.dXregY) DEALLOCATE (vadvatu)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Yadv_at_U_2d

!========================================================================

SUBROUTINE Zadv_at_U(psiu,tridcfu,wadvatuw)
!************************************************************************
!
! *Zadv_at_U* Vertical advection term for the 3-D current at U-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_U_3d
!
! Module calls - error_alloc, tvd_limiter
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiu
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4) :: tridcfu
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1) :: wadvatuw

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiu*      REAL    U-node quantity to be advected                     [psiu]
!*tridcfu*   REAL    Tridiagonal matrix for implicit vertical solution
!*wadvatuw*  REAL    Advective velocity at UW-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iopt_impl, k, npcc
REAL :: xexp, ximp
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu, maskatuw
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2du
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advflux, advflux_ce, advflux_lw, &
                     & advflux_up, array3duw, cfl, dflux, omega, upwind, wsign

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at UW-nodes          [psiu]
!*advflux_lw* REAL    Lax_Wendroff flux at UW-nodes                [psiu]*[m/s]
!*advflux_up* REAL    Upwind flux at UW-nodes                      [psiu]*[m/s]
!*omega*      REAL    TVD limiting factor
!*wsign*      REAL    Upwind factor
!*wadvatu*    REAL    Advective velocity at U-nodes                       [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Zadv_at_U'
CALL log_timer_in(npcc)

!
!1. Set parameters
!-----------------
!

IF (iopt_adv_3D.EQ.2) THEN
   xexp = 0.5*delt3d
   ximp = 0.0
   iopt_impl = 0
ELSE
   xexp = 0.5*delt3d*(1.0-theta_vadv)
   ximp = 0.5*delt3d*theta_vadv
   iopt_impl = iopt_vadv_impl
ENDIF

!
!2. Allocate arrays
!------------------
!

ALLOCATE (maskatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc,nrloc,nz/),kndlog)

ALLOCATE (maskatuw(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('maskatuw',3,(/ncloc,nrloc,nz-1/),kndlog)

IF (iopt_impl.NE.2) THEN
   ALLOCATE (advflux(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('advflux',3,(/ncloc,nrloc,nz+1/),kndrtype)
ENDIF

IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (wsign(ncloc,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('wsign',3,(/ncloc,nrloc,nz-1/),kndrtype)
   IF (iopt_impl.NE.2) THEN
      ALLOCATE (advflux_up(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('advflux_up',3,(/ncloc,nrloc,nz+1/),kndrtype)
   ENDIF
ENDIF

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      ALLOCATE (advflux_ce(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('advflux_ce',3,(/ncloc,nrloc,nz+1/),kndrtype)
   ELSE
      ALLOCATE (advflux_lw(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('advflux_lw',3,(/ncloc,nrloc,nz+1/),kndrtype)
      ALLOCATE (cfl(ncloc,nrloc,2:nz),STAT=errstat)
      CALL error_alloc('cfl',3,(/ncloc,nrloc,nz-1/),kndrtype)
   ENDIF
ENDIF

IF (iopt_adv_3D.EQ.3) THEN
   ALLOCATE (dflux(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('dflux',3,(/ncloc,nrloc,nz+1/),kndrtype)
   ALLOCATE (omega(ncloc,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('omega',3,(/ncloc,nrloc,nz-1/),kndrtype)
ENDIF

IF (iopt_impl.GT.0) THEN
   ALLOCATE (upwind(ncloc,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('upwind',3,(/ncloc,nrloc,nz-1/),kndrtype)
   IF (iopt_fld_alpha.GT.0) THEN
      ALLOCATE (array2du(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('array2du',2,(/ncloc,nrloc/),kndrtype)
   ENDIF
ENDIF

IF (iopt_impl.GT.0.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (array3duw(ncloc,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('array3duw',3,(/ncloc,nrloc,nz-1/),kndrtype)
ENDIF

!
!3. Initialise arrays
!--------------------
!
!3.1 Mask arrays
!---------------
!

maskatu = nodeatu(1:ncloc,1:nrloc,:).EQ.2
maskatuw = nodeatuw(1:ncloc,1:nrloc,2:nz).EQ.2

!
!3.2 Fluxes
!----------
!

IF (iopt_impl.NE.2) THEN
   advflux = 0.0
   IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
      advflux_up = 0.0
   ENDIF
ENDIF

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      advflux_ce = 0.0
   ELSE
      advflux_lw = 0.0
   ENDIF
ENDIF
IF (iopt_adv_3D.EQ.3) dflux = 0.0

!
!3.3 CFL number
!--------------
!

IF (iopt_impl.EQ.0.AND.(iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3)) THEN
   k_330: DO k=2,nz
      WHERE (maskatuw(:,:,k))
         cfl(:,:,k) = delt3d*wadvatuw(:,:,k)/delzatuw(1:ncloc,:,k)
      END WHERE
   ENDDO k_330
ENDIF

!
!3.4 Work space arrays
!---------------------
!

IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   WHERE (maskatuw)
      wsign = SIGN(1.0,wadvatuw(:,:,2:nz))
   END WHERE
ENDIF

!
!4. Explicit terms
!-----------------
!
!4.1 Upwind fluxes
!-----------------
!

IF ((iopt_adv_3D.EQ.1.AND.iopt_impl.NE.2).OR.iopt_adv_3D.EQ.3) THEN

   WHERE (maskatuw)
      advflux_up(:,:,2:nz) = wadvatuw(:,:,2:nz)*&
                           & ((1.0+wsign)*psiu(1:ncloc,1:nrloc,1:nz-1)&
                           & +(1.0-wsign)*psiu(1:ncloc,1:nrloc,2:nz))
   END WHERE

!  ---non-TVD case
   IF (iopt_adv_3D.NE.3) THEN
      WHERE (maskatuw)
         advflux(:,:,2:nz) = xexp*advflux_up(:,:,2:nz)
      END WHERE
   ENDIF

ENDIF

!
!4.2 Central fluxes
!------------------
!

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   
   IF (iopt_impl.GT.0) THEN
      WHERE (maskatuw)
            advflux_ce(:,:,2:nz) = wadvatuw(:,:,2:nz)*&
                                 & (psiu(1:ncloc,1:nrloc,1:nz-1)&
                                 & +psiu(1:ncloc,1:nrloc,2:nz))
         END WHERE

!     ---non-TVD case
      IF (iopt_adv_3D.EQ.2.AND.iopt_impl.NE.2) THEN
         WHERE (maskatuw)
            advflux(:,:,2:nz) = xexp*advflux_ce(:,:,2:nz)
         END WHERE
      ENDIF

!
!4.3 Lax-Wendroff fluxes
!-----------------------
!

   ELSE
      WHERE (maskatuw)
         advflux_lw(:,:,2:nz) = wadvatuw(:,:,2:nz)*&
                              & ((1.0+cfl)*psiu(1:ncloc,1:nrloc,1:nz-1)&
                              & +(1.0-cfl)*psiu(1:ncloc,1:nrloc,2:nz))
      END WHERE

!     ---non-TVD case
      IF (iopt_adv_3D.EQ.2) THEN
         WHERE (maskatuw)
            advflux(:,:,2:nz) = xexp*advflux_lw(:,:,2:nz)
         END WHERE
      ENDIF
   ENDIF

ENDIF

!
!4.4 TVD scheme
!--------------
!

IF (iopt_adv_3D.EQ.3) THEN

!  ---flux differences
   IF (iopt_impl.NE.0) THEN
      WHERE (maskatuw)
         dflux(:,:,2:nz) = advflux_ce(:,:,2:nz) - advflux_up(:,:,2:nz)
      END WHERE
   ELSE
      WHERE (maskatuw)
         dflux(:,:,2:nz) = advflux_lw(:,:,2:nz) - advflux_up(:,:,2:nz)
      END WHERE
   ENDIF

!  ---flux ratios
   WHERE (maskatuw) array3duw = 0.0
   WHERE (ABS(dflux(:,:,2:nz)).GT.eps_adv)
      array3duw = 0.5*((1.0+wsign)*dflux(:,:,1:nz-1)+&
                    & (1.0-wsign)*dflux(:,:,3:nz+1))/dflux(:,:,2:nz)
   END WHERE

!  ---TVD limiter
   omega = tvd_limiter(array3duw,maskatuw)

!  ---fluxes
   IF (iopt_impl.NE.2) THEN
      WHERE (maskatuw)
         advflux(:,:,2:nz) = xexp*(advflux_up(:,:,2:nz)+omega*dflux(:,:,2:nz))
      END WHERE
   ENDIF

ENDIF

!
!4.5 Advective term
!------------------
!

IF (iopt_impl.NE.2) THEN
   IF (iopt_fld_alpha.EQ.0) THEN
      WHERE (maskatu)
         tridcfu(:,:,:,4) = tridcfu(:,:,:,4) - &
                          & (advflux(:,:,2:nz+1)-advflux(:,:,1:nz))&
                          & /delzatu(1:ncloc,:,:)
      END WHERE
   ELSEIF (iopt_fld_alpha.GT.0) THEN
      k_452: DO k=1,nz
         WHERE (maskatu(:,:,k))
            tridcfu(:,:,k,4) = tridcfu(:,:,k,4) - alphatu_fld*&
                            & (advflux(:,:,k+1)-advflux(:,:,k))&
                            & /delzatu(1:ncloc,:,k)
         END WHERE
      ENDDO k_452
   ENDIF
ENDIF

!
!5. Implicit terms
!-----------------
!

IF (iopt_impl.GT.0) THEN

   IF (iopt_fld_alpha.EQ.0) THEN
      WHERE (maskatuw)
         array3duw = ximp*wadvatuw(:,:,2:nz)
      END WHERE
   ELSEIF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatuw(:,:,nz))
         array2du = ximp*alphatu_fld
      END WHERE
      k_512: DO k=2,nz
         WHERE (maskatuw(:,:,k))
            array3duw(:,:,k) = array2du*wadvatuw(:,:,k)
         END WHERE
      ENDDO k_512
   ENDIF

   IF (iopt_adv_3D.EQ.1) THEN
      WHERE (maskatuw) upwind = wsign
   ELSEIF (iopt_adv_3D.EQ.2) THEN
      WHERE (maskatuw) upwind = 0.0
   ELSEIF (iopt_adv_3D.EQ.3) THEN
      WHERE (maskatuw) upwind = (1.0-omega)*wsign
   ENDIF

!  ---lower fluxes
   WHERE (maskatu(:,:,2:nz))
      tridcfu(:,:,2:nz,1) = tridcfu(:,:,2:nz,1) - array3duw*&
                          & (1.0+upwind)/delzatu(1:ncloc,:,2:nz)
      tridcfu(:,:,2:nz,2) = tridcfu(:,:,2:nz,2) - array3duw*&
                          & (1.0-upwind)/delzatu(1:ncloc,:,2:nz)
   END WHERE

!  ---upper fluxes
   WHERE (maskatu(:,:,1:nz-1))
      tridcfu(:,:,1:nz-1,2) = tridcfu(:,:,1:nz-1,2) + array3duw*&
                            & (1.0+upwind)/delzatu(1:ncloc,:,1:nz-1)
      tridcfu(:,:,1:nz-1,3) = tridcfu(:,:,1:nz-1,3) + array3duw*&
                            & (1.0-upwind)/delzatu(1:ncloc,:,1:nz-1)
   END WHERE

ENDIF

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatu,maskatuw)
IF (iopt_impl.NE.2) DEALLOCATE (advflux)
IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   DEALLOCATE (wsign)
   IF (iopt_impl.NE.2) DEALLOCATE (advflux_up)
ENDIF
IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      DEALLOCATE (advflux_ce)
   ELSE
      DEALLOCATE (advflux_lw,cfl)
   ENDIF
ENDIF
IF (iopt_adv_3D.EQ.3) DEALLOCATE (dflux,omega)
IF (iopt_impl.GT.0) DEALLOCATE (upwind)
IF (iopt_impl.GT.0.AND.iopt_fld_alpha.GT.0) DEALLOCATE (array2du)
IF (iopt_impl.GT.0.OR.iopt_adv_3D.EQ.3) DEALLOCATE (array3duw)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Zadv_at_U

!========================================================================

SUBROUTINE Xadv_at_V_3d(psiv,xadvatv3d,nh,vintflag)
!************************************************************************
!
! *Xadv_at_V_3d* Advection term in X-direction for the 3-D current at V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_V_3d
!
! Module calls - error_alloc, tvd_limiter, Uarr_at_UV, Uarr_at_V
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Uarr_at_UV, Uarr_at_V
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: vintflag
INTEGER, INTENT(IN) :: nh
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiv
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,nz) :: xadvatv3d

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiv*      REAL    V-node quantity to be advected                     [psiv]
!*xadvatv3d* REAL    Advective term (times delz)                    [m*psiv/s]
!*nh*        INTEGER Halo size of local arrays
!*vintflag*  LOGICAL Update vertically integrated advective term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: info = .FALSE.
INTEGER :: i, ii, iiloc, j, k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2duv
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatv, maskatuv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2duv, denom, fcurvatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advflux, advflux_lw, advflux_up, &
              & array3duv1, array3duv2, cfl, dflux, omega, uadvatv, uadvatuv, &
              & usign, ustokesatuv

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at UV-nodes         [psiv]
!*advflux_lw* REAL    Lax_Wendroff flux at UV- nodes              [psiv]*[m/s]
!*advflux_up* REAL    Upwind flux at UV nodes                     [psiv]*[m/s]
!*omega*      REAL    TVD limiting factor
!*usign*      REAL    Upwind factor
!*uadvatv*    REAL    Advective velocity at V-nodes                      [m/s]
!*uadvatuv*   REAL    Advective velocity at UV-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xadv_at_V_3d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (advflux(ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('advflux',3,(/ncloc+1,nrloc,nz/),kndrtype)

ALLOCATE (array2duv(2-nh:ncloc+nh,nrloc),STAT=errstat)
CALL error_alloc('array2duv',2,(/ncloc+2*nh-1,nrloc/),kndrtype)

ALLOCATE (array3duv1(ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('array3duv1',3,(/ncloc+1,nrloc,nz/),kndrtype)

ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE (maskatuv(2-nh:ncloc+nh,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatuv',3,(/ncloc+2*nh-1,nrloc,nz/),kndlog)

ALLOCATE (mask2duv(2-nh:ncloc+nh,nrloc),STAT=errstat)
CALL error_alloc('mask2duv',2,(/ncloc+2*nh-1,nrloc/),kndlog)

ALLOCATE (maskatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc,nrloc,nz/),kndlog)

ALLOCATE (uadvatuv(2-nh:ncloc+nh,nrloc,nz),STAT=errstat)
CALL error_alloc('uadvatuv',3,(/ncloc+2*nh-1,nrloc,nz/),kndrtype)

IF (iopt_waves_curr.EQ.1) THEN
   ALLOCATE (ustokesatuv(2-nh:ncloc+nh,nrloc,nz),STAT=errstat)
   CALL error_alloc('ustokesatuv',3,(/ncloc+2*nh-1,nrloc,nz/),kndrtype)
ENDIF

ALLOCATE (usign(2-nh:ncloc+nh,nrloc,nz),STAT=errstat)
CALL error_alloc('usign',3,(/ncloc+2*nh-1,nrloc,nz/),kndrtype)

IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (advflux_up(2-nh:ncloc+nh,nrloc,nz),STAT=errstat)
   CALL error_alloc('advflux_up',3,(/ncloc+2*nh-1,nrloc,nz/),kndrtype)
ENDIF

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (advflux_lw(2-nh:ncloc+nh,nrloc,nz),STAT=errstat)
   CALL error_alloc('advflux_lw',3,(/ncloc+2*nh-1,nrloc,nz/),&
        & kndrtype)
   ALLOCATE (cfl(2-nh:ncloc+nh,nrloc,nz),STAT=errstat)
   CALL error_alloc('cfl',3,(/ncloc+2*nh-1,nrloc,nz/),kndrtype)
ENDIF

IF (iopt_adv_3D.EQ.3) THEN
   ALLOCATE (array3duv2(ncloc+1,nrloc,nz),STAT=errstat)
   CALL error_alloc('array3duv2',3,(/ncloc+1,nrloc,nz/),kndrtype)
   ALLOCATE (dflux(2-nh:ncloc+nh,nrloc,nz),STAT=errstat)
   CALL error_alloc('dflux',3,(/ncloc+2*nh-1,nrloc,nz/),kndrtype)
   ALLOCATE (omega(ncloc+1,nrloc,nz),STAT=errstat)
   CALL error_alloc('omega',3,(/ncloc+1,nrloc,nz/),kndrtype)
ENDIF

IF (.NOT.dYregX) THEN
   ALLOCATE (fcurvatv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fcurvatu',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (uadvatv(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('uadvatv',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

!
!2. Initialise arrays
!--------------------
!
!2.1 Mask arrays
!---------------
!

maskatv = nodeatv(1:ncloc,1:nrloc,:).EQ.2
maskatuv = nodeatuv(2-nh:ncloc+nh,1:nrloc,:).GT.0
mask2duv = node2duv(2-nh:ncloc+nh,1:nrloc).GT.0

!
!2.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   advflux_up = 0.0
ENDIF
IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv_3D.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!2.3 Metric factor in advective flux
!-----------------------------------
!

IF (dYregX) THEN
   k_231: DO k=1,nz
      WHERE (maskatuv(1:ncloc+1,:,k))
         array3duv1(:,:,k) = delzatuv(:,1:nrloc,k)
      END WHERE
   ENDDO k_231
ELSE
   k_232: DO k=1,nz
      WHERE (maskatuv(1:ncloc+1,:,k))
         array3duv1(:,:,k) = delyatuv(1:ncloc+1,1:nrloc)*delzatuv(:,1:nrloc,k)
      END WHERE
   ENDDO k_232
ENDIF

!
!2.4 Advective velocity
!----------------------
!

CALL Uarr_at_UV(uvel_old(2-nh:ncloc+nh,0:nrloc,:),uadvatuv,1,1,&
             & (/2-nh,0,1/),(/ncloc+nh,nrloc,nz/),1,iarr_uvel_old,info)
IF (iopt_waves_curr.EQ.1) THEN
   CALL Uarr_at_UV(ustokesatu(2-nh:ncloc+nh,0:nrloc,:),ustokesatuv,1,1,&
                & (/2-nh,0,1/),(/ncloc+nh,nrloc,nz/),1,iarr_ustokesatu,info)
   uadvatuv = uadvatuv + ustokesatuv
ENDIF

!
!2.5 Sign of velocity
!--------------------
!

WHERE (maskatuv)
   usign = SIGN(1.0,uadvatuv)
END WHERE

!
!2.6 CFL number
!--------------
!

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   WHERE (mask2duv)
      array2duv = delt3d/delxatuv(2-nh:ncloc+nh,1:nrloc)
   END WHERE
   k_260: DO k=1,nz
      WHERE (maskatuv(:,:,k))
         cfl(:,:,k) = array2duv*uadvatuv(:,:,k)
      END WHERE
   ENDDO k_260
ENDIF

!
!2.7 Work space arrays
!---------------------
!

IF (dYregX) THEN
   denom = 2.0*delxatv(1:ncloc,1:nrloc)
ELSE
   denom = 2.0*delxatv(1:ncloc,1:nrloc)*delyatv(1:ncloc,1:nrloc)
   fcurvatv = (delyatuv(2:ncloc+1,1:nrloc)-delyatuv(1:ncloc,1:nrloc))/denom
ENDIF

!
!3. Explicit terms
!-----------------
!

maskatuv = nodeatuv(2-nh:ncloc+nh,1:nrloc,:).EQ.1.OR.&
         & nodeatuv(2-nh:ncloc+nh,1:nrloc,:).EQ.3

!
!3.1 Upwind fluxes
!-----------------
!

IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN

   k_310: DO k=1,nz
      WHERE (maskatuv(:,:,k))
         advflux_up(:,:,k) = uadvatuv(:,:,k)*&
                          & ((1.0+usign(:,:,k))*psiv(1-nh:ncloc+nh-1,1:nrloc,k)&
                          & +(1.0-usign(:,:,k))*psiv(2-nh:ncloc+nh,1:nrloc,k))
      END WHERE
   ENDDO k_310

!  ---non-TVD case
   IF (iopt_adv_3D.EQ.1) THEN
      WHERE (maskatuv(1:ncloc+1,:,:))
         advflux = array3duv1*advflux_up(1:ncloc+1,:,:)
      END WHERE
   ENDIF

ENDIF

!
!3.2 Lax-Wendroff fluxes
!-----------------------
!

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN

   k_320: DO k=1,nz
      WHERE (maskatuv(:,:,k))
         advflux_lw(:,:,k) = uadvatuv(:,:,k)*&
                           & ((1.0+cfl(:,:,k))*psiv(1-nh:ncloc+nh-1,1:nrloc,k)&
                           & +(1.0-cfl(:,:,k))*psiv(2-nh:ncloc+nh,1:nrloc,k))
      END WHERE
   ENDDO k_320

!  ---non-TVD case
   IF (iopt_adv_3D.EQ.2) THEN
      WHERE (maskatuv(1:ncloc+1,:,:))
         advflux = array3duv1*advflux_lw(1:ncloc+1,:,:)
      END WHERE
   ENDIF

ENDIF

!
!3.3 TVD scheme
!--------------
!

IF (iopt_adv_3D.EQ.3) THEN

!  ---flux differences
   WHERE (maskatuv)
      dflux = advflux_lw - advflux_up
   END WHERE
   
!  ---flux ratios
   WHERE (maskatuv(1:ncloc+1,:,:)) array3duv2 = 0.0
   k_330: DO k=1,nz
      WHERE (ABS(dflux(1:ncloc+1,:,k)).GT.eps_adv)
         array3duv2(:,:,k) = 0.5*(&
            & (1.0+usign(1:ncloc+1,:,k))*dflux(0:ncloc,:,k)+&
            & (1.0-usign(1:ncloc+1,:,k))*dflux(2:ncloc+2,:,k))&
            & /dflux(1:ncloc+1,:,k)
      END WHERE
   ENDDO k_330

!  ---TVD limiter
   omega = tvd_limiter(array3duv2,maskatuv(1:ncloc+1,:,:))

!  ---fluxes
   WHERE (maskatuv(1:ncloc+1,:,:))
      advflux = array3duv1*(advflux_up(1:ncloc+1,:,:)+&
                          & omega*dflux(1:ncloc+1,:,:))
   END WHERE

ENDIF

!
!3.4 Open boundary fluxes
!------------------------
!

iiloc_340: DO iiloc=1,nobxloc_ext
   i = iobxloc(iiloc); j = jobxloc(iiloc)
   ii = indexobx(iiloc)

   SELECT CASE (itypobx(ii))

!
!3.4.1 Zero gradient condition
!-----------------------------
!

      CASE (0)
         IF (westobx(ii)) THEN
            IF (i.LT.ncloc+1) advflux(i,j,:) = advflux(i+1,j,:)
         ELSE
            IF (i.GT.1) advflux(i,j,:) = advflux(i-1,j,:)
         ENDIF

!
!3.4.2 Upwind scheme
!-------------------
!

      CASE (1)
         IF (westobx(ii)) THEN
            k_3421: DO k=1,nz
               IF (usign(i,j,k).LT.0.0) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*uadvatuv(i,j,k)*&
                                 & psiv(i,j,k)
               ELSEIF (nodeatv(i-1,j,k).GT.2.OR.nodeatu(i,j-1,k).LE.1.OR.&
                     & nodeatu(i,j,k).LE.1) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*uadvatuv(i,j,k)*&
                                 & psiv(i-1,j,k)
               ELSE
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*uadvatuv(i,j,k)*&
                                 & psiv(i,j,k)
               ENDIF
            ENDDO k_3421
         ELSE
            k_3422: DO k=1,nz
               IF (usign(i,j,k).GT.0.0) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*uadvatuv(i,j,k)*&
                                 & psiv(i-1,j,k)
               ELSEIF (nodeatv(i,j,k).GT.2.OR.nodeatu(i,j-1,k).LE.1.OR.&
                     & nodeatu(i,j,k).LE.1) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*uadvatuv(i,j,k)*&
                                 & psiv(i,j,k)
               ELSE
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*uadvatuv(i,j,k)*&
                                 & psiv(i-1,j,k)
               ENDIF
            ENDDO k_3422
         ENDIF

!
!3.4.3 Tangential condition with external values
!-----------------------------------------------
!

      CASE (2)
         IF (westobx(ii)) THEN
            k_3431: DO k=1,nz
               IF (usign(i,j,k).LT.0.0) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*uadvatuv(i,j,k)*&
                                 & psiv(i,j,k)
               ELSE
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*uadvatuv(i,j,k)*&
                                 & vvelobx(ii,k)  
               ENDIF
            ENDDO k_3431
         ELSE
            k_3432: DO k=1,nz
               IF (usign(i,j,k).GT.0.0) THEN
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*uadvatuv(i,j,k)*&
                                 & psiv(i-1,j,k)
               ELSE
                  advflux(i,j,k) = 2.0*array3duv1(i,j,k)*uadvatuv(i,j,k)*&
                                 & vvelobx(ii,k)  
               ENDIF
            ENDDO k_3432
         ENDIF

   END SELECT

ENDDO iiloc_340

!
!3.5 Advective term
!------------------
!

k_350: DO k=1,nz
   WHERE (maskatv(:,:,k))
      xadvatv3d(:,:,k) = (advflux(2:ncloc+1,:,k)-advflux(1:ncloc,:,k))/denom
   END WHERE
ENDDO k_350

!
!3.6 Curvature term
!------------------
!

IF (.NOT.dYregX) THEN
   CALL Uarr_at_V(uvel_old(1:ncloc+1,0:nrloc,:),uadvatv,1,2,(/1,0,1/),&
                  & (/ncloc+1,nrloc,nz/),1,iarr_uvel_old,info)
   IF (iopt_waves_curr.EQ.0) THEN
      k_361: DO k=1,nz
         WHERE (maskatv(:,:,k))
            xadvatv3d(:,:,k) = xadvatv3d(:,:,k) + 2.0*delzatv(:,1:nrloc,k)*&
                             & uadvatv(:,:,k)*psiv(1:ncloc,1:nrloc,k)*fcurvatv
         END WHERE
      ENDDO k_361
   ELSE
      k_362: DO k=1,nz
         WHERE (maskatv(:,:,k))
            xadvatv3d(:,:,k) = xadvatv3d(:,:,k) + 2.0*delzatv(:,1:nrloc,k)*&
                             & uadvatv(:,:,k)*&
                             & (psiv(1:ncloc,1:nrloc,k)+&
                             & vstokesatv(1:ncloc,1:nrloc,k))*fcurvatv
         END WHERE
      ENDDO k_362
   ENDIF
ENDIF

!
!3.7 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN
   k_370: DO k=1,nz
      WHERE (maskatv(:,:,k))
         xadvatv3d(:,:,k) = alphatv_fld*xadvatv3d(:,:,k)
      END WHERE
   ENDDO k_370
ENDIF

!
!3.8 Apply relaxation scheme
!---------------------------
!

IF (iopt_obc_advrlx.EQ.1) THEN
   k_380: DO k=1,nz
      WHERE (maskatv(:,:,k))
         xadvatv3d(:,:,k) = rlxobcatv*xadvatv3d(:,:,k)
      END WHERE
   ENDDO k_380
ENDIF

!
!3.9 Update depth-integrated advective term
!-----------------------------------------
!

IF (vintflag) THEN
   k_390: DO k=1,nz
      WHERE (maskatv(:,:,k))
         vdevint = vdevint - xadvatv3d(:,:,k)
      END WHERE
   ENDDO k_390
ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatv,maskatuv,mask2duv)
DEALLOCATE (advflux,array2duv,array3duv1,denom,uadvatuv,usign)
IF (iopt_waves_curr.EQ.1) DEALLOCATE (ustokesatuv)
IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) DEALLOCATE (advflux_up)
IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) DEALLOCATE (advflux_lw,cfl)
IF (iopt_adv_3D.EQ.3) DEALLOCATE (array3duv2,dflux,omega)
IF (.NOT.dYregX) DEALLOCATE (fcurvatv,uadvatv)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Xadv_at_V_3d

!========================================================================

SUBROUTINE Xadv_at_V_2d(psiv,xadvatv2d,nh,vintflag)
!************************************************************************
!
! *Xadv_at_V_2d* Advective term in X-direction for the 2-D current at V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_V_2d
!
! Module calls - error_alloc, tvd_limiter, Uarr_at_UV, Uarr_at_V
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Uarr_at_UV, Uarr_at_V
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(INOUT) :: vintflag
INTEGER, INTENT(IN) :: nh
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo) :: psiv
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: xadvatv2d

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiv*      REAL    V-node quantity to be advected                     [psiv]
!*xadvatv2d* REAL    Advective term                                   [psiv/s]
!*nh*        INTEGER Halo size of local arrays
!*vintflag*  LOGICAL Update vertically integrated advective term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ii, iiloc, j, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatv, maskatuv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: advflux, advflux_lw, advflux_up,&
                 & array2duv, cfl, dflux, denom, omega, uadvatv, uadvatuv, &
                 & umstokesatuv, usign

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at UV-nodes         [psiv]
!*advflux_lw* REAL    Lax_Wendroff flux at UV-nodes               [psiv]*[m/s]
!*advflux_up* REAL    Upwind flux at UV-nodes                     [psiv]*[m/s]
!*omega*      REAL    TVD limiting factor
!*usign*      REAL    Upwind factor
!*uadvatv*    REAL    Advective velocity at V-nodes                      [m/s]
!*uadvatuv*   REAL    Advective velocity at UV-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xadv_at_V_2d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (advflux(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('advflux',2,(/ncloc+1,nrloc/),kndrtype)

ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE (maskatuv(2-nh:ncloc+nh,nrloc),STAT=errstat)
CALL error_alloc('maskatuv',2,(/ncloc+2*nh-1,nrloc/),kndlog)

ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc,nrloc/),kndlog)

ALLOCATE (uadvatuv(2-nh:ncloc+nh,nrloc),STAT=errstat)
CALL error_alloc('uadvatuv',2,(/ncloc+2*nh-1,nrloc/),kndrtype)

IF (iopt_waves_curr.EQ.1) THEN
   ALLOCATE (umstokesatuv(2-nh:ncloc+nh,nrloc),STAT=errstat)
   CALL error_alloc('umstokesatuv',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
ENDIF

ALLOCATE (usign(2-nh:ncloc+nh,nrloc),STAT=errstat)
CALL error_alloc('usign',2,(/ncloc+2*nh-1,nrloc/),kndrtype)

IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN
   ALLOCATE (advflux_up(2-nh:ncloc+nh,nrloc),STAT=errstat)
   CALL error_alloc('advflux_up',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
ENDIF

IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN
   ALLOCATE (advflux_lw(2-nh:ncloc+nh,nrloc),STAT=errstat)
   CALL error_alloc('advflux_lw',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
   ALLOCATE (cfl(2-nh:ncloc+nh,nrloc),STAT=errstat)
   CALL error_alloc('cfl',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
ENDIF

IF (iopt_adv_2D.EQ.3) THEN
   ALLOCATE (array2duv(ncloc+1,nrloc),STAT=errstat)
   CALL error_alloc('array2duv',2,(/ncloc+1,nrloc/),kndrtype)
   ALLOCATE (dflux(2-nh:ncloc+nh,nrloc),STAT=errstat)
   CALL error_alloc('dflux',2,(/ncloc+2*nh-1,nrloc/),kndrtype)
   ALLOCATE (omega(ncloc+1,nrloc),STAT=errstat)
   CALL error_alloc('omega',2,(/ncloc+1,nrloc/),kndrtype)
ENDIF

IF (.NOT.dYregX) THEN
   ALLOCATE (uadvatv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('vadvatu',2,(/ncloc,nrloc/),kndrtype)
ENDIF

!
!2. Initialise arrays
!--------------------
!
!2.1 Mask arrays
!---------------
!

maskatv = node2dv(1:ncloc,1:nrloc).EQ.2
maskatuv = node2duv(2-nh:ncloc+nh,1:nrloc).GT.0

!
!2.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN
   advflux_up = 0.0
ENDIF
IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv_2D.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!2.3 Advective velocities
!------------------------
!

CALL Uarr_at_UV(umvel_old(2-nh:ncloc+nh,0:nrloc),uadvatuv,1,1,(/2-nh,0,nz/),&
             & (/ncloc+nh,nrloc,nz/),1,iarr_umvel_old,.TRUE.)
IF (iopt_waves_curr.EQ.1) THEN
   CALL Uarr_at_UV(umstokesatu(2-nh:ncloc+nh,0:nrloc),umstokesatuv,1,1,&
                & (/2-nh,0,nz/),(/ncloc+nh,nrloc,nz/),1,iarr_umstokesatu,.TRUE.)
   uadvatuv = uadvatuv + umstokesatuv
ENDIF

!
!2.4 Sign of velocity
!--------------------
!

WHERE (maskatuv)
   usign = SIGN(1.0,uadvatuv)
END WHERE

!
!2.5 Work space arrays
!---------------------
!

IF (dYregX) THEN
   WHERE (maskatv)
      denom = delxatv(1:ncloc,1:nrloc)
   END WHERE
ELSE
   WHERE (maskatv)
      denom = delxatv(1:ncloc,1:nrloc)*delyatv(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!3. Explicit terms
!-----------------
!

maskatuv = node2duv(2-nh:ncloc+nh,1:nrloc).EQ.1.OR.&
         & node2duv(2-nh:ncloc+nh,1:nrloc).EQ.3

!
!3.1 Upwind fluxes
!-----------------
!

IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN

   WHERE (maskatuv)
      advflux_up = uadvatuv*((1.0+usign)*psiv(1-nh:ncloc+nh-1,1:nrloc)&
                          & +(1.0-usign)*psiv(2-nh:ncloc+nh,1:nrloc))
   END WHERE
!  ---non-TVD case
   IF (iopt_adv_2D.EQ.1) THEN
      IF (dYregX) THEN
         WHERE (maskatuv(1:ncloc+1,:))
            advflux = advflux_up(1:ncloc+1,:)
         END WHERE
      ELSE
         WHERE (maskatuv(1:ncloc+1,:))
            advflux = delyatuv(1:ncloc+1,1:nrloc)*advflux_up(1:ncloc+1,:)
         END WHERE
      ENDIF
   ENDIF

ENDIF

!
!3.2 Lax-Wendroff fluxes
!-----------------------
!

IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN

   WHERE (maskatuv)
      cfl = delt2d*uadvatuv/delxatuv(2-nh:ncloc+nh,1:nrloc)
      advflux_lw = uadvatuv*((1.0+cfl)*psiv(1-nh:ncloc+nh-1,1:nrloc)&
                          & +(1.0-cfl)*psiv(2-nh:ncloc+nh,1:nrloc))
   END WHERE

!  ---non-TVD case
   IF (iopt_adv_2D.EQ.2) THEN
      IF (dYregX) THEN
         WHERE (maskatuv(1:ncloc+1,:))
            advflux = advflux_lw(1:ncloc+1,:)
         END WHERE
      ELSE
         WHERE (maskatuv(1:ncloc+1,:))
            advflux = delyatuv(1:ncloc+1,1:nrloc)*advflux_lw(1:ncloc+1,:)
         END WHERE
      ENDIF
   ENDIF

ENDIF

!
!3.3 TVD scheme
!--------------
!

IF (iopt_adv_2D.EQ.3) THEN

!  ---flux differences
   WHERE (maskatuv)
      dflux = advflux_lw - advflux_up
   END WHERE

!  ---flux ratios
   WHERE (maskatuv(1:ncloc+1,:)) array2duv = 0.0
   WHERE (ABS(dflux(1:ncloc+1,:)).GT.eps_adv)
      array2duv = 0.5*(&
                  & (1.0+usign(1:ncloc+1,:))*dflux(0:ncloc,:)+&
                  & (1.0-usign(1:ncloc+1,:))*dflux(2:ncloc+2,:))&
                  & /dflux(1:ncloc+1,:)
   END WHERE

!  ---TVD limiter
   omega = tvd_limiter(array2duv,maskatuv(1:ncloc+1,:))

!  ---fluxes
   IF (dYregX) THEN
      WHERE (maskatuv(1:ncloc+1,:))
         advflux = advflux_up(1:ncloc+1,:)+omega*dflux(1:ncloc+1,:)
      END WHERE
   ELSE
      WHERE (maskatuv(1:ncloc+1,:))
         advflux = delyatuv(1:ncloc+1,1:nrloc)*&
                 & (advflux_up(1:ncloc+1,:)+omega*dflux(1:ncloc+1,:))
      END WHERE
   ENDIF

ENDIF

!
!3.4 Open boundary fluxes
!------------------------
!

iiloc_340: DO iiloc=1,nobxloc_ext
   i = iobxloc(iiloc); j = jobxloc(iiloc)
   ii = indexobx(iiloc)

   SELECT CASE (ityp2dobx(ii))

!
!3.4.1 Zero gradient condition
!-----------------------------
!

      CASE (0)
         IF (westobx(ii)) THEN
            IF (i.LT.ncloc+1) advflux(i,j) = advflux(i+1,j)
         ELSE
            IF (i.GT.1) advflux(i,j) = advflux(i-1,j)
         ENDIF

!
!3.4.2 Upwind scheme
!-------------------
!

      CASE (1)
         IF (westobx(ii)) THEN
            IF (usign(i,j).LT.0.0) THEN
               advflux(i,j) = 2.0*uadvatuv(i,j)*psiv(i,j)
            ELSEIF (node2dv(i-1,j).GT.2.OR.node2du(i,j-1).LE.1.OR.&
                  & node2du(i,j).LE.1) THEN
               advflux(i,j) = 2.0*uadvatuv(i,j)*psiv(i-1,j)
            ELSE
               advflux(i,j) = 2.0*uadvatuv(i,j)*psiv(i,j)
            ENDIF
         ELSE
            IF (usign(i,j).GT.0.0) THEN
               advflux(i,j) = 2.0*uadvatuv(i,j)*psiv(i-1,j)
            ELSEIF (node2dv(i,j).GT.2.OR.node2du(i,j-1).LE.1.OR.&
                  & node2du(i,j).LE.1) THEN
               advflux(i,j) = 2.0*uadvatuv(i,j)*psiv(i,j)
            ELSE
               advflux(i,j) = 2.0*uadvatuv(i,j)*psiv(i-1,j)
            ENDIF
         ENDIF
         IF (.NOT.dYregX) advflux(i,j) = delyatuv(i,j)*advflux(i,j)

!
!3.4.3 Tangential condition with external values
!-----------------------------------------------
!

      CASE (2)
         IF (westobx(ii)) THEN
            IF (usign(i,j).LT.0.0) THEN
               advflux(i,j) = 2.0*uadvatuv(i,j)*psiv(i,j)
            ELSE
               advflux(i,j) =  uadvatuv(i,j)*vdatobx(ii,2)
            ENDIF
         ELSE
            IF (usign(i,j).GT.0.0) THEN
               advflux(i,j) = 2.0*uadvatuv(i,j)*psiv(i-1,j)
            ELSE
               advflux(i,j) = uadvatuv(i,j)*vdatobx(ii,2)
            ENDIF
         ENDIF
         IF (.NOT.dYregX) advflux(i,j) = delyatuv(i,j)*advflux(i,j)

   END SELECT

ENDDO iiloc_340

!
!3.5 Advective term
!------------------
!

WHERE (maskatv)
   xadvatv2d = 0.5*(advflux(2:ncloc+1,:)-advflux(1:ncloc,:))/denom
END WHERE

!
!3.6 Curvature term
!------------------
!

IF (.NOT.dYregX) THEN
   CALL Uarr_at_V(umvel_old(1:ncloc+1,0:nrloc),uadvatv,1,2,(/1,0,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_umvel_old,.TRUE.)
   IF (iopt_waves_curr.EQ.0) THEN
      WHERE (maskatv)
         xadvatv2d = xadvatv2d + uadvatv*psiv(1:ncloc,1:nrloc)*&
                 & (delyatuv(2:ncloc+1,1:nrloc)-delyatuv(1:ncloc,1:nrloc))/denom
      END WHERE
   ELSE
      WHERE (maskatv)
         xadvatv2d = xadvatv2d + uadvatv*&
                 & (psiv(1:ncloc,1:nrloc)+vmstokesatv(1:ncloc,1:nrloc))*&
                 & (delyatuv(2:ncloc+1,1:nrloc)-delyatuv(1:ncloc,1:nrloc))/denom
      END WHERE
   ENDIF
ENDIF

!
!3.7 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN
   WHERE (maskatv)
      xadvatv2d = alphatv_fld*xadvatv2d
   END WHERE
ENDIF

!
!3.8 Apply relaxation scheme
!---------------------------
!

IF (iopt_obc_advrlx.EQ.1) THEN
   WHERE (maskatv)
      xadvatv2d = rlxobcatv*xadvatv2d
   END WHERE
ENDIF

!
!3.9 Update depth-integrated advective term
!------------------------------------------
!

IF (vintflag) THEN
   WHERE (maskatv)
      vdevint = vdevint + xadvatv2d
   END WHERE
   vintflag = .FALSE.
ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatuv,maskatv)
DEALLOCATE (advflux,denom,uadvatuv,usign)
IF (iopt_waves_curr.EQ.1) DEALLOCATE (umstokesatuv)
IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) DEALLOCATE (advflux_up)
IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) DEALLOCATE (advflux_lw,cfl)
IF (iopt_adv_2D.EQ.3) DEALLOCATE (array2duv,dflux,omega)
IF (.NOT.dYregX) DEALLOCATE (uadvatv)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Xadv_at_V_2d

!========================================================================

SUBROUTINE Yadv_at_V_3d(psiv,yadvatv3d,nh,vintflag)
!************************************************************************
!
! *Yadv_at_V_3d* Advection term in Y-direction for the 3-D current at V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_V_3d
!
  ! Module calls - error_alloc, tvd_limiter, Carr_at_V, Uarr_at_V, Varr_at_C
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Carr_at_V, Uarr_at_V, Varr_at_C
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: vintflag
INTEGER, INTENT(IN) :: nh
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiv
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,nz) :: yadvatv3d

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiv*      REAL    V-node quantity to be advected                     [psiv]
!*yadvat3d*  REAL    Advective term (times delz)                    [m*psiv/s]
!*nh*        INTEGER Halo size of local arrays
!*vintflag*  LOGICAL Update vertically integrated advective term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: info = .FALSE.
INTEGER :: k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: advflux, advflux_lw, advflux_up,&
               & array2dc1, array2dc2, cfl, denom, dflux, fcurvatv, omega, &
               & uadvatv, ustokesatv, vsign
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: vadvatc

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at C-nodes          [psiv]
!*advflux_lw* REAL    Lax_Wendroff flux at C-nodes                [psiv]*[m/s]
!*advflux_up* REAL    Upwind flux at C-nodes                      [psiv]*[m/s]
!*omega*      REAL    TVD limiting factor
!*vsign*      REAL    Upwind factor
!*vadvatc*    REAL    Advective velocity at C-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Yadv_at_V_3d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatc(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc,nrloc+2*nh-1/),kndlog)

ALLOCATE (maskatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc,nz/),kndlog)

ALLOCATE (advflux(ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('advflux',2,(/ncloc,nrloc+1/),kndrtype)

ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE (vadvatc(ncloc,1-nh:nrloc+nh-1,nz),STAT=errstat)
CALL error_alloc('vadvatc',3,(/ncloc,nrloc+2*nh-1,nz/),kndrtype)

IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (advflux_up(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
   CALL error_alloc('advflux_up',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
   ALLOCATE (vsign(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
   CALL error_alloc('vsign',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
ENDIF

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (advflux_lw(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
   CALL error_alloc('advflux_lw',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
   ALLOCATE (array2dc1(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
   CALL error_alloc('array2dc1',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
   ALLOCATE (cfl(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
   CALL error_alloc('cfl',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
ENDIF

IF (iopt_adv_3D.EQ.3) THEN
   ALLOCATE (array2dc2(ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('array2dc1',2,(/ncloc,nrloc+1/),kndrtype)
   ALLOCATE (dflux(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
   CALL error_alloc('dflux',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
   ALLOCATE (omega(ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('omega',2,(/ncloc,nrloc+1/),kndrtype)
ENDIF

IF (.NOT.dXregY) THEN
   ALLOCATE (fcurvatv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fcurvatv',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (uadvatv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('uadvatv',2,(/ncloc,nrloc/),kndrtype)
   IF (iopt_waves_curr.EQ.1) THEN
      ALLOCATE (ustokesatv(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('ustokesatv',2,(/ncloc,nrloc/),kndrtype)
   ENDIF
ENDIF

!
!2. Initialise arrays
!--------------------
!
!2.1 Mask arrays
!---------------
!

maskatv = nodeatv(1:ncloc,1:nrloc,:).EQ.2
maskatc = nodeatc(1:ncloc,1-nh:nrloc+nh-1).GT.0

!
!2.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   advflux_up = 0.0
ENDIF
IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv_3D.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!2.3 Work space arrays
!---------------------
!

IF (dXregY) THEN
   denom = 2.0*delyatv(1:ncloc,1:nrloc)
ELSE
   denom = 2.0*delxatv(1:ncloc,1:nrloc)*delyatv(1:ncloc,1:nrloc)
   fcurvatv = (delxatc(1:ncloc,1:nrloc)-delxatc(1:ncloc,0:nrloc-1))/denom
ENDIF

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   WHERE (maskatc)
      array2dc1 = delt3d/delyatc(1:ncloc,1-nh:nrloc+nh-1)
   END WHERE
ENDIF

!
!3. Explicit terms
!-----------------
!

!
!3.1 Advective velocity
!----------------------
!

CALL Varr_at_C(vvel_old(1:ncloc,1-nh:nrloc+nh,:),vadvatc,1,1,(/1,1-nh,1/),&
     & (/ncloc,nrloc+nh,nz/),1,iarr_vvel_old,info)
IF (iopt_waves_curr.EQ.1) THEN
   vadvatc = vadvatc + vstokesatc(1:ncloc,1-nh:nrloc+nh-1,:)
ENDIF

k_300: DO k=1,nz


!
!3.2 Upwind fluxes
!-----------------
!

   IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN

      WHERE (maskatc)
         vsign = SIGN(1.0,vadvatc(:,:,k))
         advflux_up = vadvatc(:,:,k)*&
                            & ((1.0+vsign)*psiv(1:ncloc,1-nh:nrloc+nh-1,k)&
                            & +(1.0-vsign)*psiv(1:ncloc,2-nh:nrloc+nh,k))
      END WHERE

!     ---non-TVD case
      IF (iopt_adv_3D.EQ.1) THEN
         IF (dXregY) THEN
            WHERE (maskatc(:,0:nrloc))
               advflux = delzatc(1:ncloc,0:nrloc,k)*advflux_up(:,0:nrloc)
            END WHERE
         ELSE
            WHERE (maskatc(:,0:nrloc))
               advflux = delxatc(1:ncloc,0:nrloc)*delzatc(1:ncloc,0:nrloc,k)*&
                       & advflux_up(:,0:nrloc)
            END WHERE
         ENDIF
      ENDIF

   ENDIF

!
!3.3 Lax-Wendroff fluxes
!------------------------
!

   IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN

      WHERE (maskatc)
         cfl = array2dc1*vadvatc(:,:,k)
         advflux_lw = vadvatc(:,:,k)*((1.0+cfl)*psiv(1:ncloc,1-nh:nrloc+nh-1,k)&
                            & +(1.0-cfl)*psiv(1:ncloc,2-nh:nrloc+nh,k))
      END WHERE

!     ---non-TVD case
      IF (iopt_adv_3D.EQ.2) THEN
         IF (dXregY) THEN
            WHERE (maskatc(:,0:nrloc))
               advflux = delzatc(0:ncloc,1:nrloc,k)*advflux_lw(:,0:nrloc)
            END WHERE
         ELSE
            WHERE (maskatc(:,0:nrloc))
               advflux = delxatc(1:ncloc,0:nrloc)*delzatc(0:ncloc,1:nrloc,k)*&
                       & advflux_lw(:,0:nrloc)
            END WHERE
         ENDIF
      ENDIF

   ENDIF

!
!3.4 TVD scheme
!--------------
!

   IF (iopt_adv_3D.EQ.3) THEN

!     ---flux differences
      WHERE (maskatc)
         dflux = advflux_lw - advflux_up
      END WHERE
   
!     ---flux ratios
      WHERE (maskatc(:,0:nrloc)) array2dc2 = 0.0
      WHERE (ABS(dflux(:,0:nrloc)).GT.eps_adv)
         array2dc2 = 0.5*((1.0+vsign(:,0:nrloc))*dflux(:,-1:nrloc-1)+&
                        & (1.0-vsign(:,0:nrloc))*dflux(:,1:nrloc+1))&
                        & /dflux(:,0:nrloc)
      END WHERE

!     ---TVD limiter
      omega = tvd_limiter(array2dc2,maskatc(:,0:nrloc)) 

!     ---fluxes
      IF (dXregY) THEN
         WHERE (maskatc(:,0:nrloc))
            advflux = delzatc(1:ncloc,0:nrloc,k)*&
                    & (advflux_up(:,0:nrloc)+omega*dflux(:,0:nrloc))
         END WHERE
      ELSE
         WHERE (maskatc(:,0:nrloc))
            advflux = delxatc(1:ncloc,0:nrloc)*delzatc(1:ncloc,0:nrloc,k)*&
                    & (advflux_up(:,0:nrloc)+omega*dflux(:,0:nrloc))
         END WHERE
      ENDIF

   ENDIF

!
!3.5 Advective term
!------------------
!

   WHERE (maskatv(:,:,k))
      yadvatv3d(:,:,k) = (advflux(:,1:nrloc)-advflux(:,0:nrloc-1))/denom
   END WHERE

!
!3.6 Curvature term
!------------------
!

   IF (.NOT.dXregY) THEN
      CALL Uarr_at_V(uvel_old(1:ncloc+1,0:nrloc,k),uadvatv,1,2,(/1,0,k/),&
                        & (/ncloc+1,nrloc,k/),1,iarr_uvel_old,info)
      IF (iopt_waves_curr.EQ.0) THEN
         WHERE (maskatv(:,:,k))
            yadvatv3d(:,:,k) = yadvatv3d(:,:,k) - 2.0*delzatv(:,1:nrloc,k)*&
                                                & uadvatv*uadvatv*fcurvatv
         END WHERE
      ELSE
         CALL Carr_at_V(ustokesatc(1:ncloc,0:nrloc,k),ustokesatv,1,1,&
                     & (/1,0,k/),(/ncloc,nrloc,k/),1,iarr_ustokesatc,info)
         WHERE (maskatv(:,:,k))
            yadvatv3d(:,:,k) = yadvatv3d(:,:,k) - 2.0*delzatv(:,1:nrloc,k)*&
                                        & uadvatv*(uadvatv+ustokesatv)*fcurvatv
         END WHERE
      ENDIF
   ENDIF

!
!3.7 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

   IF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatv(:,:,k))
         yadvatv3d(:,:,k) = alphatv_fld*yadvatv3d(:,:,k)
      END WHERE
   ENDIF

!
!3.8 Apply relaxation scheme
!---------------------------
!

   IF (iopt_obc_advrlx.EQ.1) THEN
      WHERE (maskatv(:,:,k))
         yadvatv3d(:,:,k) = rlxobcatv*yadvatv3d(:,:,k)
      END WHERE
   ENDIF

!
!3.9 Update depth-integrated advective term
!------------------------------------------
!

   IF (vintflag) THEN
      WHERE (maskatv(:,:,k))
         vdevint = vdevint - yadvatv3d(:,:,k)
      END WHERE
   ENDIF

ENDDO k_300

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatc,maskatv)
DEALLOCATE (advflux,denom,vadvatc)
IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) DEALLOCATE (advflux_up,vsign)
IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   DEALLOCATE (advflux_lw,array2dc1,cfl)
ENDIF
IF (iopt_adv_3D.EQ.3) DEALLOCATE (array2dc2,dflux,omega)
IF (.NOT.dXregY) DEALLOCATE (fcurvatv,uadvatv)
IF (.NOT.dXregY.AND.iopt_waves_curr.EQ.1) DEALLOCATE (ustokesatv)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Yadv_at_V_3d

!========================================================================

SUBROUTINE Yadv_at_V_2d(psiv,yadvatv2d,nh,vintflag)
!************************************************************************
!
! *Yadv_at_V_2d* Advection term in Y-direction for the 2-D current at V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_V_2d
!
! Module calls - error_alloc, tvd_limiter, Carr_at_V, Uarr_at_V, Varr_at_C
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE obconds
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Carr_at_V, Uarr_at_V, Varr_at_C
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(INOUT) :: vintflag
INTEGER, INTENT(IN) :: nh
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo) :: psiv
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: yadvatv2d

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiv*      REAL    V-node quantity to be advected                    [psiv]
!*yadvatv2d* REAL    Advective term                                  [psiv/s]
!*nh*        INTEGER Halo size of local arrays
!*vintflag*  LOGICAL Update vertically integrated advective term if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: advflux, advflux_lw, advflux_up,&
               & array2dc, cfl, denom, dflux, omega, uadvatv, umstokesatv, &
               & vadvatc, vsign

!
! Name        Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at C-nodes          [psiv]
!*advflux_lw* REAL    Lax_Wendroff flux at C-nodes                [psiv]*[m/s]
!*advflux_up* REAL    Upwind flux at C-nodes                      [psiv]*[m/s]
!*omega*      REAL    TVD limiting factor
!*vsign*      REAL    Upwind factor
!*vadvatc*    REAL    Advective velocity at C-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Yadv_at_V_at_2d'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatc(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc,nrloc+2*nh-1/),kndlog)

ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)

ALLOCATE (advflux(ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('advflux',2,(/ncloc,nrloc+1/),kndrtype)

ALLOCATE (vadvatc(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
CALL error_alloc('vadvatc',2,(/ncloc,nrloc+2*nh-1/),kndrtype)

IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN
   ALLOCATE (advflux_up(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
   CALL error_alloc('advflux_up',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
   ALLOCATE (vsign(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
   CALL error_alloc('vsign',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
ENDIF

IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN
   ALLOCATE (advflux_lw(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
   CALL error_alloc('advflux_lw',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
   ALLOCATE (cfl(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
   CALL error_alloc('cfl',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
ENDIF

IF (iopt_adv_2D.EQ.3) THEN
   ALLOCATE (array2dc(ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('array2dc',2,(/ncloc,nrloc+1/),kndrtype)
   ALLOCATE (dflux(ncloc,1-nh:nrloc+nh-1),STAT=errstat)
   CALL error_alloc('dflux',2,(/ncloc,nrloc+2*nh-1/),kndrtype)
   ALLOCATE (omega(ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('omega',2,(/ncloc,nrloc+1/),kndrtype)
ENDIF

IF (.NOT.dXregY) THEN
   ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (uadvatv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('uadvatv',2,(/ncloc,nrloc/),kndrtype)
   IF (iopt_waves_curr.EQ.1) THEN
      ALLOCATE (umstokesatv(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('umstokesatv',2,(/ncloc,nrloc/),kndrtype)
   ENDIF
ENDIF

!
!2. Initialise arrays
!--------------------
!
!2.1 Mask arrays
!---------------
!

maskatv = node2dv(1:ncloc,1:nrloc).EQ.2
maskatc = nodeatc(1:ncloc,1-nh:nrloc+nh-1).GT.0

!
!2.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN
   advflux_up = 0.0
ENDIF
IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv_2D.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!2.3 Work space array
!--------------------
!

IF (.NOT.dXregY) THEN
   WHERE (maskatv)
      denom = delxatv(1:ncloc,1:nrloc)*delyatv(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!2.4 Advecting velocity
!----------------------
!

CALL Varr_at_C(vmvel_old(1:ncloc,1-nh:nrloc+nh),vadvatc,1,1,(/1,1-nh,nz/),&
            & (/ncloc,nrloc+nh,nz/),1,iarr_vmvel_old,.TRUE.)
IF (iopt_waves_curr.EQ.1) THEN
   vadvatc = vadvatc + vmstokesatc(1:ncloc,1-nh:nrloc+nh-1)
ENDIF

!
!3. Explicit terms
!-----------------
!
!3.1 Upwind fluxes
!-----------------
!

IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) THEN

   WHERE (maskatc)
      vsign = SIGN(1.0,vadvatc)
      advflux_up = vadvatc*((1.0+vsign)*psiv(1:ncloc,1-nh:nrloc+nh-1)&
                         & +(1.0-vsign)*psiv(1:ncloc,2-nh:nrloc+nh))
   END WHERE

!  ---non-TVD case
   IF (iopt_adv_2D.EQ.1) THEN
      IF (dXregY) THEN
         WHERE (maskatc(:,0:nrloc))
            advflux = advflux_up(:,0:nrloc)
         END WHERE
      ELSE
         WHERE (maskatc(:,0:nrloc))
            advflux = delxatc(1:ncloc,0:nrloc)*advflux_up(:,0:nrloc)
         END WHERE
      ENDIF
   ENDIF

ENDIF

!
!3.2 Lax-Wendroff fluxes
!-----------------------
!

IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) THEN

   WHERE (maskatc)
      cfl = delt2d*vadvatc/delyatc(1:ncloc,1-nh:nrloc+nh-1)
      advflux_lw = vadvatc*((1.0+cfl)*psiv(1:ncloc,1-nh:nrloc+nh-1)&
                         & +(1.0-cfl)*psiv(1:ncloc,2-nh:nrloc+nh))
   END WHERE

!  ---non-TVD case
   IF (iopt_adv_2D.EQ.2) THEN
      IF (dXregY) THEN
         WHERE (maskatc(:,0:nrloc))
            advflux = advflux_lw(:,0:nrloc)
         END WHERE
      ELSE
         WHERE (maskatc(:,0:nrloc))
            advflux = delxatc(1:ncloc,0:nrloc)*advflux_lw(:,0:nrloc)
         END WHERE
      ENDIF
   ENDIF

ENDIF

!
!3.3 TVD scheme
!--------------
!

IF (iopt_adv_2D.EQ.3) THEN

!  ---flux differences
   WHERE (maskatc)
      dflux = advflux_lw - advflux_up
   END WHERE
   
!  ---flux ratios
   WHERE (maskatc(:,0:nrloc)) array2dc = 0.0
   WHERE (ABS(dflux(:,0:nrloc)).GT.eps_adv)
      array2dc = 0.5*((1.0+vsign(:,0:nrloc))*dflux(:,-1:nrloc-1)+&
                    & (1.0-vsign(:,0:nrloc))*dflux(:,1:nrloc+1))&
                    & /dflux(:,0:nrloc)
   END WHERE

!  ---TVD limiter
   omega = tvd_limiter(array2dc,maskatc(:,0:nrloc))

!  ---fluxes
   IF (dXregY) THEN
      WHERE (maskatc(:,0:nrloc))
         advflux = advflux_up(:,0:nrloc)+omega*dflux(:,0:nrloc)
      END WHERE
   ELSE
      WHERE (maskatc(:,0:nrloc))
         advflux = delxatc(1:ncloc,0:nrloc)*&
                 & (advflux_up(:,0:nrloc)+omega*dflux(:,0:nrloc))
      END WHERE
   ENDIF

ENDIF

!
!3.4 Advective term
!------------------
!
!---regular grid
IF (dXregY) THEN
   WHERE (maskatv)
      yadvatv2d = 0.5*(advflux(:,1:nrloc)-advflux(:,0:nrloc-1))&
                & /delyatv(1:ncloc,1:nrloc)
   END WHERE

!---irregular grid
ELSE
   WHERE (maskatv)
      yadvatv2d = 0.5*(advflux(:,1:nrloc)-advflux(:,0:nrloc-1))/denom
   END WHERE
ENDIF

!
!3.5 Curvature term
!-------------------
!

IF (.NOT.dXregY) THEN
   CALL Uarr_at_V(umvel_old(1:ncloc+1,0:nrloc),uadvatv,1,2,(/1,0,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_umvel_old,.TRUE.)
   IF (iopt_waves_curr.EQ.0) THEN
      WHERE (maskatv)
         yadvatv2d = yadvatv2d - deptotatv(1:ncloc,1:nrloc)*uadvatv*uadvatv*&
            & (delxatc(1:ncloc,1:nrloc)-delxatc(1:ncloc,0:nrloc-1))/denom
      END WHERE
   ELSE
      CALL Carr_at_V(umstokesatc(1:ncloc,0:nrloc),umstokesatv,1,1,&
                  & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_umstokesatc,.TRUE.)
      WHERE (maskatv)
         yadvatv2d = yadvatv2d - deptotatv(1:ncloc,1:nrloc)*uadvatv*&
                   & (uadvatv+umstokesatv)*&
                   & (delxatc(1:ncloc,1:nrloc)-delxatc(1:ncloc,0:nrloc-1))/denom
      END WHERE
   ENDIF
ENDIF

!
!3.6 Apply reduction factor for flooding/drying scheme
!-----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN
   WHERE (maskatv)
      yadvatv2d = alphatv_fld*yadvatv2d
   END WHERE
ENDIF

!
!3.7 Apply relaxation scheme
!---------------------------
!

IF (iopt_obc_advrlx.EQ.1) THEN
   WHERE (maskatv)
      yadvatv2d = rlxobcatv*yadvatv2d
   END WHERE
ENDIF

!
!3.8 Update depth-integrated advective term
!------------------------------------------
!

IF (vintflag) THEN
   WHERE (maskatv)
      vdevint = vdevint + yadvatv2d
   END WHERE
   vintflag = .FALSE.
ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatc,maskatv)
DEALLOCATE (advflux,vadvatc)
IF (iopt_adv_2D.EQ.1.OR.iopt_adv_2D.EQ.3) DEALLOCATE (advflux_up,vsign)
IF (iopt_adv_2D.EQ.2.OR.iopt_adv_2D.EQ.3) DEALLOCATE (advflux_lw,cfl)
IF (iopt_adv_2D.EQ.3) DEALLOCATE (array2dc,dflux,omega)
IF (.NOT.dXregY) DEALLOCATE (denom,uadvatv)
IF (.NOT.dXregY.AND.iopt_waves_curr.EQ.1) DEALLOCATE (umstokesatv)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Yadv_at_V_2d

!========================================================================

SUBROUTINE Zadv_at_V(psiv,tridcfv,wadvatvw)
!************************************************************************
!
! *Zadv_at_V* Vertical advection term for the 3-D current at V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_V_3d
!
! Module calls - error_alloc, tvd_limiter
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psiv
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,4) :: tridcfv
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1) :: wadvatvw

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiv*      REAL    V-node quantity to be advected                     [psiv]
!*tridcfv*   REAL    Tridiagonal matrix for implicit vertical solution
!*wadvatvw*  REAL    Advective velocity at VW-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iopt_impl, k, npcc
REAL :: xexp, ximp
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatv, maskatvw
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advflux, advflux_ce, advflux_lw, &
                                  & advflux_up, array3dvw, cfl, dflux, omega, &
                                  & upwind, wsign

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at VW-nodes          [psiv]
!*advflux_ce* REAL    Central flux at VW-nodes                     [psiv]*[m/s]
!*advflux_lw* REAL    Lax_Wendroff flux at VW-nodes                [psiv]*[m/s]
!*advflux_up* REAL    Upwind flux at VW-nodes                      [psiv]*[m/s]
!*omega*      REAL    TVD limiting factor
!*wsign*      REAL    Upwind factor
!*wadvatv*    REAL    Advective velocity at V-nodes                       [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Zadv_at_V'
CALL log_timer_in(npcc)

!
!1. Set parameters
!-----------------
!

IF (iopt_adv_3D.EQ.2) THEN
   xexp = 0.5*delt3d
   ximp = 0.0
   iopt_impl = 0
ELSE
   xexp = 0.5*delt3d*(1.0-theta_vadv)
   ximp = 0.5*delt3d*theta_vadv
   iopt_impl = iopt_vadv_impl
ENDIF

!
!2. Allocate arrays
!------------------
!

ALLOCATE (array2dv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE (maskatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc,nz/),kndlog)

ALLOCATE (maskatvw(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('maskatvw',3,(/ncloc,nrloc,nz-1/),kndlog)

IF (iopt_impl.NE.2) THEN
   ALLOCATE (advflux(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('advflux',3,(/ncloc,nrloc,nz+1/),kndrtype)
ENDIF

IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (wsign(ncloc,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('wsign',3,(/ncloc,nrloc,nz-1/),kndrtype)
   IF (iopt_impl.NE.2) THEN
      ALLOCATE (advflux_up(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('advflux_up',3,(/ncloc,nrloc,nz+1/),kndrtype)
   ENDIF
ENDIF

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      ALLOCATE (advflux_ce(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('advflux_ce',3,(/ncloc,nrloc,nz+1/),kndrtype)
   ELSE
      ALLOCATE (advflux_lw(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('advflux_lw',3,(/ncloc,nrloc,nz+1/),kndrtype)
      ALLOCATE (cfl(ncloc,nrloc,2:nz),STAT=errstat)
      CALL error_alloc('cfl',3,(/ncloc,nrloc,nz-1/),kndrtype)
   ENDIF
ENDIF

IF (iopt_adv_3D.EQ.3) THEN
   ALLOCATE (dflux(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('dflux',3,(/ncloc,nrloc,nz+1/),kndrtype)
   ALLOCATE (omega(ncloc,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('omega',3,(/ncloc,nrloc,nz-1/),kndrtype)
ENDIF

IF (iopt_impl.GT.0) THEN
   ALLOCATE (upwind(ncloc,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('upwind',3,(/ncloc,nrloc,nz-1/),kndrtype)
ENDIF

IF (iopt_impl.GT.0.OR.iopt_adv_3D.EQ.3) THEN
   ALLOCATE (array3dvw(ncloc,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('array3dvw',3,(/ncloc,nrloc,nz-1/),kndrtype)
ENDIF

!
!3. Initialise arrays
!--------------------
!
!3.1 Mask arrays
!---------------
!

maskatv = nodeatv(1:ncloc,1:nrloc,:).EQ.2
maskatvw = nodeatvw(1:ncloc,1:nrloc,2:nz).EQ.2

!
!3.2 Fluxes
!----------
!

IF (iopt_impl.NE.2) THEN
   advflux = 0.0
   IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
      advflux_up = 0.0
   ENDIF
ENDIF

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      advflux_ce = 0.0
   ELSE
      advflux_lw = 0.0
   ENDIF
ENDIF

IF (iopt_adv_3D.EQ.3) dflux = 0.0

!
!3.3 CFL number
!--------------
!

IF (iopt_impl.EQ.0.AND.(iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3)) THEN
   k_330: DO k=2,nz
      WHERE (maskatvw(:,:,k))
         cfl(:,:,k) = delt3d*wadvatvw(:,:,k)/delzatvw(:,1:nrloc,k)
      END WHERE
   ENDDO k_330
ENDIF

!
!3.4 Work space arrays
!---------------------
!

IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   WHERE (maskatvw)
      wsign = SIGN(1.0,wadvatvw(:,:,2:nz))
   END WHERE
ENDIF

!
!4. Explicit terms
!-----------------
!
!4.1 Upwind fluxes
!-----------------
!

IF ((iopt_adv_3D.EQ.1.AND.iopt_impl.NE.2).OR.iopt_adv_3D.EQ.3) THEN

   WHERE (maskatvw)
      advflux_up(:,:,2:nz) = wadvatvw(:,:,2:nz)*&
                          & ((1.0+wsign)*psiv(1:ncloc,1:nrloc,1:nz-1)&
                          & +(1.0-wsign)*psiv(1:ncloc,1:nrloc,2:nz))
   END WHERE

!  ---non-TVD case
   IF (iopt_adv_3D.NE.3) THEN
      WHERE (maskatvw)
         advflux(:,:,2:nz) = xexp*advflux_up(:,:,2:nz)
      END WHERE
   ENDIF

ENDIF

!
!4.2 Central fluxes
!------------------
!

IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN

   IF (iopt_impl.GT.0) THEN
      WHERE (maskatvw)
         advflux_ce(:,:,2:nz) = wadvatvw(:,:,2:nz)*&
                              & (psiv(1:ncloc,1:nrloc,1:nz-1)&
                              & +psiv(1:ncloc,1:nrloc,2:nz))
      END WHERE

!     ---non-TVD case
      IF (iopt_adv_3D.EQ.2.AND.iopt_impl.NE.2) THEN
         WHERE (maskatvw)
            advflux(:,:,2:nz) = xexp*advflux_ce(:,:,2:nz)
         END WHERE
      ENDIF

!
!4.3 Lax-Wendroff fluxes
!-----------------------
!

   ELSE
      WHERE (maskatvw)
         advflux_lw(:,:,2:nz) = wadvatvw(:,:,2:nz)*&
                              & ((1.0+cfl)*psiv(1:ncloc,1:nrloc,1:nz-1)&
                              & +(1.0-cfl)*psiv(1:ncloc,1:nrloc,2:nz))
      END WHERE

!     ---non-TVD case
      IF (iopt_adv_3D.EQ.2) THEN
         WHERE (maskatvw)
            advflux(:,:,2:nz) = xexp*advflux_lw(:,:,2:nz)
         END WHERE
      ENDIF
   ENDIF

ENDIF

!
!4.4 TVD scheme
!--------------
!

IF (iopt_adv_3D.EQ.3) THEN

!  ---flux differences
   IF (iopt_impl.NE.0) THEN
      WHERE (maskatvw)
         dflux(:,:,2:nz) = advflux_ce(:,:,2:nz) - advflux_up(:,:,2:nz)
      END WHERE
   ELSE
      WHERE (maskatvw)
         dflux(:,:,2:nz) = advflux_lw(:,:,2:nz) - advflux_up(:,:,2:nz)
      END WHERE
   ENDIF

!  ---flux ratios
   WHERE (maskatvw) array3dvw = 0.0
   WHERE (ABS(dflux(:,:,2:nz)).GT.eps_adv)
      array3dvw = 0.5*((1.0+wsign)*dflux(:,:,1:nz-1)+&
                     & (1.0-wsign)*dflux(:,:,3:nz+1))/dflux(:,:,2:nz)
   END WHERE

!  ---TVD limiter
   omega = tvd_limiter(array3dvw,maskatvw)

!  ---fluxes
   IF (iopt_impl.NE.2) THEN
      WHERE (maskatvw)
         advflux(:,:,2:nz) = xexp*(advflux_up(:,:,2:nz)+omega*dflux(:,:,2:nz))
      END WHERE
   ENDIF

ENDIF

!
!4.5 Advective term
!------------------
!

IF (iopt_impl.NE.2) THEN
   IF (iopt_fld_alpha.EQ.0) THEN
      WHERE (maskatv)
         tridcfv(:,:,:,4) = tridcfv(:,:,:,4)- &
                         & (advflux(:,:,2:nz+1)-advflux(:,:,1:nz))&
                         & /delzatv(:,1:nrloc,1:nz)
      END WHERE
   ELSEIF (iopt_fld_alpha.GT.0) THEN
      k_452: DO k=1,nz
         WHERE (maskatv(:,:,k))
            tridcfv(:,:,k,4) = tridcfv(:,:,k,4)- alphatv_fld*&
                            & (advflux(:,:,k+1)-advflux(:,:,k))&
                            & /delzatv(:,1:nrloc,k)
         END WHERE
      ENDDO k_452
   ENDIF
ENDIF

!
!5. Implicit terms
!-----------------
!

IF (iopt_impl.GT.0) THEN

   IF (iopt_fld_alpha.EQ.0) THEN
      WHERE (maskatvw)
         array3dvw = ximp*wadvatvw(:,:,2:nz)
      END WHERE
   ELSEIF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatvw(:,:,nz))
         array2dv = ximp*alphatv_fld
      END WHERE
      k_512: DO k=2,nz
         WHERE (maskatvw(:,:,k))
            array3dvw(:,:,k) = array2dv*wadvatvw(:,:,k)
         END WHERE
      ENDDO k_512
   ENDIF

   IF (iopt_adv_3D.EQ.1) THEN
      WHERE (maskatvw) upwind = wsign
   ELSEIF (iopt_adv_3D.EQ.2) THEN
      WHERE (maskatvw) upwind = 0.0
   ELSEIF (iopt_adv_3D.EQ.3) THEN
      WHERE (maskatvw) upwind = (1.0-omega)*wsign
   ENDIF

!  ---lower fluxes
   WHERE (maskatv(:,:,2:nz))
      tridcfv(:,:,2:nz,1) = tridcfv(:,:,2:nz,1) - array3dvw*&
                          & (1.0+upwind)/delzatv(:,1:nrloc,2:nz)
      tridcfv(:,:,2:nz,2) = tridcfv(:,:,2:nz,2) - array3dvw*&
                          & (1.0-upwind)/delzatv(:,1:nrloc,2:nz)
   END WHERE

!  ---upper fluxes
   WHERE (maskatv(:,:,1:nz-1))
      tridcfv(:,:,1:nz-1,2) = tridcfv(:,:,1:nz-1,2) + array3dvw*&
                            & (1.0+upwind)/delzatv(:,1:nrloc,1:nz-1)
      tridcfv(:,:,1:nz-1,3) = tridcfv(:,:,1:nz-1,3) + array3dvw*&
                            & (1.0-upwind)/delzatv(:,1:nrloc,1:nz-1)
   END WHERE

ENDIF

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatv,maskatvw)
DEALLOCATE (array2dv)
IF (iopt_impl.NE.2) DEALLOCATE (advflux)
IF (iopt_adv_3D.EQ.1.OR.iopt_adv_3D.EQ.3) THEN
   DEALLOCATE (wsign)
   IF (iopt_impl.NE.2) DEALLOCATE (advflux_up)
ENDIF
IF (iopt_adv_3D.EQ.2.OR.iopt_adv_3D.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      DEALLOCATE (advflux_ce)
   ELSE
      DEALLOCATE (advflux_lw,cfl)
   ENDIF
ENDIF
IF (iopt_adv_3D.EQ.3) DEALLOCATE (dflux,omega)
IF (iopt_impl.GT.0) DEALLOCATE (upwind)
IF (iopt_impl.GT.0.OR.iopt_adv_3D.EQ.3) DEALLOCATE (array3dvw)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Zadv_at_V

!========================================================================

SUBROUTINE Xadv_at_W(psiw,tridcfw,uadvatuw,nh,klo,kup)
!************************************************************************
!
! *Xadv_at_W* Advection term in X-direction for a quantity at W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_W
!
! Module calls - error_alloc, tvd_limiter
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
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, nh
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1)::psiw
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz+1,4) :: tridcfw
REAL, INTENT(IN), DIMENSION(2-nhalo:ncloc+nhalo,nrloc,2:nz) :: uadvatuw

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiw*      REAL    W-node quantity to be advected                     [psiw]
!*tridcfw*   REAL    Tridiagonal matrix for implicit (vertical) solution
!*uadvatuw*  REAL    Advective velocity at UW-nodes                      [m/s]
!*nh*        INTEGER Halo size of local arrays
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ic, ii, iiloc, i1, j, k, npcc
REAL :: xexp
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2duw
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatuw
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2du
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advflux, advflux_lw, advflux_up, &
                      & array3duw1, array3duw2, cfl, denom, dflux, omega, usign

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at UW-nodes         [psiw]
!*advflux_lw* REAL    Lax_Wendroff flux at UW-nodes               [psiw]*[m/s]
!*advflux_up* REAL    Upwind flux at UW-nodes                     [psiw]*[m/s]
!*omega*      REAL    TVD limiting factor
!*usign*      REAL    Upwind factor
!*uadvatw*    REAL    Advective velocity at W-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Xadv_at_W'
CALL log_timer_in(npcc)

!
!1. Set parameters
!-----------------
!

xexp = 0.5*delt3d

!
!2. Allocate arrays
!------------------
!

ALLOCATE (advflux(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('advflux',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)

ALLOCATE (array2du(2-nh:ncloc+nh,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc+2*nh-1,nrloc/),kndrtype)

ALLOCATE (array3duw1(ncloc+1,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('array3duw1',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)

ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)

ALLOCATE (maskatuw(2-nh:ncloc+nh,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('maskatuw',3,(/ncloc+2*nh-1,nrloc,kup-klo+1/),kndlog)

ALLOCATE (mask2duw(2-nh:ncloc+nh,nrloc),STAT=errstat)
CALL error_alloc('maskatuw',2,(/ncloc+2*nh-1,nrloc/),kndlog)

ALLOCATE (usign(2-nh:ncloc+nh,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('usign',3,(/ncloc+2*nh-1,nrloc,kup-klo+1/),kndrtype)

IF (iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3) THEN
   ALLOCATE (advflux_up(2-nh:ncloc+nh,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('advflux_up',3,(/ncloc+2*nh-1,nrloc,kup-klo+1/),kndrtype)
ENDIF

IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN
   ALLOCATE (advflux_lw(2-nh:ncloc+nh,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('advflux_lw',3,(/ncloc+2*nh-1,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (cfl(2-nh:ncloc+nh,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('cfl',3,(/ncloc+2*nh-1,nrloc,kup-klo+1/),kndrtype)
ENDIF

IF (iopt_adv_turb.EQ.3) THEN
   ALLOCATE (array3duw2(ncloc+1,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('array3duw2',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (dflux(2-nh:ncloc+nh,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('dflux',3,(/ncloc+2*nh-1,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (omega(ncloc+1,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('omega',3,(/ncloc+1,nrloc,kup-klo+1/),kndrtype)
ENDIF

!
!3. Initialise arrays
!--------------------
!
!
!3.1 Mask array
!--------------
!

maskatuw = nodeatuw(2-nh:ncloc+nh,1:nrloc,klo:kup).EQ.2
mask2duw = node2du(2-nh:ncloc+nh,1:nrloc).EQ.2

!
!3.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3) THEN
   advflux_up = 0.0
ENDIF
IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv_turb.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!3.3 Metric factor in advective flux
!-----------------------------------
!

IF (dYregX) THEN
   WHERE (maskatuw(1:ncloc+1,:,:))
      array3duw1 = xexp*delzatuw(1:ncloc+1,:,klo:kup)
   END WHERE
ELSE
   WHERE (mask2duw(1:ncloc+1,:))
      array2du(1:ncloc+1,:) = xexp*delyatu(1:ncloc+1,1:nrloc)
   END WHERE
   k_330: DO k=klo,kup
      WHERE (maskatuw(1:ncloc+1,:,k))
         array3duw1(:,:,k) = array2du(1:ncloc+1,:)*delzatuw(1:ncloc+1,:,k)
      END WHERE
   ENDDO k_330
ENDIF

!
!3.4 Sign of velocity
!--------------------
!

WHERE (nodeatuw(2-nh:ncloc+nh,1:nrloc,klo:kup).GT.1)
   usign = SIGN(1.0,uadvatuw(2-nh:ncloc+nh,:,klo:kup))
END WHERE

!
!3.5 CFL number
!--------------
!

IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN
   WHERE (mask2duw)
      array2du = delt3d/delxatu(2-nh:ncloc+nh,1:nrloc)
   END WHERE
   k_320: DO k=klo,kup
      WHERE (maskatuw(:,:,k))
         cfl(:,:,k) = array2du*uadvatuw(2-nh:ncloc+nh,:,k)
      END WHERE
   ENDDO k_320
ENDIF

!
!3.6 Work space arrays
!---------------------
!

IF (dYregX) THEN
   k_361: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delxatc(1:ncloc,1:nrloc)*delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_361
ELSE
   k_362: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_362
ENDIF

!
!4. Explicit terms
!-----------------
!
!4.1 Upwind fluxes
!-----------------
!

IF (iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3) THEN

   k_410: DO k=klo,kup
      WHERE (maskatuw(:,:,k))
         advflux_up(:,:,k) = uadvatuw(2-nh:ncloc+nh,:,k)*&
                          & ((1.0+usign(:,:,k))*psiw(1-nh:ncloc+nh-1,1:nrloc,k)&
                          & +(1.0-usign(:,:,k))*psiw(2-nh:ncloc+nh,1:nrloc,k))
      END WHERE
   ENDDO k_410

!  ---non-TVD case
   IF (iopt_adv_turb.EQ.1) THEN
      WHERE (maskatuw(1:ncloc+1,:,:))
         advflux = array3duw1*advflux_up(1:ncloc+1,:,:)
      END WHERE
   ENDIF
   
ENDIF

!
!4.2 Lax-Wendroff fluxes
!-----------------------
!

IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN

   k_420: DO k=klo,kup
      WHERE (maskatuw(:,:,k))
         advflux_lw(:,:,k) = uadvatuw(2-nh:ncloc+nh,:,k)*&
                           & ((1.0+cfl(:,:,k))*psiw(1-nh:ncloc+nh-1,1:nrloc,k)&
                           & +(1.0-cfl(:,:,k))*psiw(2-nh:ncloc+nh,1:nrloc,k))
      END WHERE
   ENDDO k_420

!  ---non-TVD case
   IF (iopt_adv_turb.EQ.2) THEN
      WHERE (maskatuw(1:ncloc+1,:,:))
         advflux = array3duw1*advflux_lw(1:ncloc+1,:,:)
      END WHERE
      
   ENDIF

ENDIF

!
!4.3 TVD scheme
!--------------
!

IF (iopt_adv_turb.EQ.3) THEN
   
!  ---flux differences
   WHERE (maskatuw)
      dflux = advflux_lw - advflux_up
   END WHERE

!  ---flux ratios
   WHERE (maskatuw(1:ncloc+1,:,:)) array3duw2 = 0.0
   k_430: DO k=klo,kup
      WHERE (ABS(dflux(1:ncloc+1,:,k)).GT.eps_adv)
         array3duw2(:,:,k) = 0.5*(&
            & (1.0+usign(1:ncloc+1,:,k))*dflux(0:ncloc,:,k)+&
            & (1.0-usign(1:ncloc+1,:,k))*dflux(2:ncloc+2,:,k))&
            & /dflux(1:ncloc+1,:,k)
      END WHERE
   ENDDO k_430

!  ---TVD limiter
   omega = tvd_limiter(array3duw2,maskatuw(1:ncloc+1,:,:))

!  ---fluxes
   WHERE (maskatuw(1:ncloc+1,:,:))
      advflux = array3duw1*(advflux_up(1:ncloc+1,:,:)+&
                          & omega*dflux(1:ncloc+1,:,:))
   END WHERE

ENDIF

!
!4.4 Open boundary fluxes
!------------------------
!

iiloc_440: DO iiloc=1,nobuloc_ext
   i = iobuloc(iiloc); j = jobuloc(iiloc)
   ii = indexobu(iiloc); ic = MERGE(i,i-1,westobu(ii))
   IF (ic.GE.1.AND.ic.LE.ncloc) THEN
      i1 = MERGE(i+1,i-1,westobu(ii))
      advflux(i,j,:) = advflux(i1,j,:)
   ENDIF
ENDDO iiloc_440

!
!4.5 Advective term
!------------------
!

k_450: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfw(:,:,k,4) = tridcfw(:,:,k,4) - &
                       & (advflux(2:ncloc+1,:,k)-advflux(1:ncloc,:,k))&
                       & /denom(:,:,k)
   END WHERE
ENDDO k_450

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatuw,mask2duw)
DEALLOCATE (advflux,array2du,array3duw1,denom,usign)
IF (iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3) DEALLOCATE (advflux_up)
IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) DEALLOCATE (advflux_lw,cfl)
IF (iopt_adv_turb.EQ.3) DEALLOCATE (array3duw2,dflux,omega)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Xadv_at_W

!========================================================================

SUBROUTINE Yadv_at_W(psiw,tridcfw,vadvatvw,nh,klo,kup)
!************************************************************************
!
! *Yadv_at_W* Advection term in Y-direction for a quantity at W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_W
!
! Module calls - error_alloc, tvd_limiter
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
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup, nh
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1)::psiw
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz+1,4) :: tridcfw
REAL, INTENT(IN), DIMENSION(ncloc,2-nhalo:nrloc+nhalo,2:nz) :: vadvatvw

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*psiw*      REAL    W-node quantity to be advected                     [psiw]
!*tridcfw*   REAL    Tridiagonal matrix for implicit (vertical) solution
!*vadvatvw*  REAL    Advective velocity at VW-nodes                      [m/s]
!*nh*        INTEGER Halo size of local arrays
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, j, jc, jj, jjloc, j1, k, npcc
REAL :: xexp
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2dvw
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatvw
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advflux, advflux_lw, advflux_up, &
                     & array3dvw1, array3dvw2, cfl, denom, dflux, omega, vsign

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at V-nodes          [psiw]
!*advflux_lw* REAL    Lax_Wendroff flux at VW-nodes               [psiw]*[m/s]
!*advflux_up* REAL    Upwind flux at VW-nodes                     [psiw]*[m/s]
!*omega*      REAL    TVD limiting factor
!*vsign*      REAL    Upwind factor
!*vadvatw*    REAL    Advective velocity at W-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Yadv_at_W'
CALL log_timer_in(npcc)

!
!1. Set parameters
!-----------------
!

xexp = 0.5*delt3d

!
!2. Allocate arrays
!------------------
!

ALLOCATE (advflux(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('advflux',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)

ALLOCATE (array2dv(ncloc,2-nh:nrloc+nh),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc+2*nh-1/),kndrtype)

ALLOCATE (array3dvw1(ncloc,nrloc+1,klo:kup),STAT=errstat)
CALL error_alloc('array3dvw1',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)

ALLOCATE (denom(ncloc,nrloc,klo:kup),STAT=errstat)
CALL error_alloc('denom',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)

ALLOCATE (maskatvw(ncloc,2-nh:nrloc+nh,klo:kup),STAT=errstat)
CALL error_alloc('maskatvw',3,(/ncloc,nrloc+2*nh-1,kup-klo+1/),kndlog)

ALLOCATE (mask2dvw(ncloc,2-nh:nrloc+nh),STAT=errstat)
CALL error_alloc('mask2dvw',3,(/ncloc,nrloc+2*nh-1/),kndlog)

ALLOCATE (vsign(ncloc,2-nh:nrloc+nh,klo:kup),STAT=errstat)
CALL error_alloc('vsign',3,(/ncloc,nrloc+2*nh-1,kup-klo+1/),kndrtype)

IF (iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3) THEN
   ALLOCATE (advflux_up(ncloc,2-nh:nrloc+nh,klo:kup),STAT=errstat)
   CALL error_alloc('advflux_up',3,(/ncloc,nrloc+2*nh-1,kup-klo+1/),kndrtype)
ENDIF

IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN
   ALLOCATE (advflux_lw(ncloc,2-nh:nrloc+nh,klo:kup),STAT=errstat)
   CALL error_alloc('advflux_lw',3,(/ncloc,nrloc+2*nh-1,kup-klo+1/),kndrtype)
   ALLOCATE (cfl(ncloc,2-nh:nrloc+nh,klo:kup),STAT=errstat)
   CALL error_alloc('cfl',3,(/ncloc,nrloc+2*nh-1,kup-klo+1/),kndrtype)
ENDIF

IF (iopt_adv_turb.EQ.3) THEN
   ALLOCATE (array3dvw2(ncloc,nrloc+1,klo:kup),STAT=errstat)
   CALL error_alloc('array3dvw2',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)
   ALLOCATE (dflux(ncloc,2-nh:nrloc+nh,klo:kup),STAT=errstat)
   CALL error_alloc('dflux',3,(/ncloc,nrloc+2*nh-1,kup-klo+1/),kndrtype)
   ALLOCATE (omega(ncloc,nrloc+1,klo:kup),STAT=errstat)
   CALL error_alloc('omega',3,(/ncloc,nrloc+1,kup-klo+1/),kndrtype)
ENDIF

!
!3. Initialise arrays
!--------------------
!
!3.1 Mask array
!--------------
!

maskatvw = nodeatvw(1:ncloc,2-nh:nrloc+nh,klo:kup).EQ.2
mask2dvw = node2dv(1:ncloc,2-nh:nrloc+nh).EQ.2

!
!3.2 Fluxes
!----------
!

advflux = 0.0
IF (iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3) THEN
   advflux_up = 0.0
ENDIF
IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN
   advflux_lw = 0.0
ENDIF
IF (iopt_adv_turb.EQ.3) THEN
   dflux = 0.0
ENDIF

!
!3.3 Metric factor in advective flux
!-----------------------------------
!

IF (dXregY) THEN
   WHERE (maskatvw(:,1:nrloc+1,:))
      array3dvw1 = xexp*delzatvw(:,1:nrloc+1,klo:kup)
   END WHERE
ELSE
   WHERE (mask2dvw(:,1:nrloc+1))
      array2dv(:,1:nrloc+1) = xexp*delxatv(1:ncloc,1:nrloc+1)
   END WHERE
   k_330: DO k=klo,kup
      WHERE (maskatvw(:,1:nrloc+1,k))
         array3dvw1(:,:,k) = array2dv(:,1:nrloc+1)*delzatvw(:,1:nrloc+1,k)
      END WHERE
   ENDDO k_330
ENDIF

!
!3.4 Sign of velocity
!--------------------
!

WHERE (nodeatvw(1:ncloc,2-nh:nrloc+nh,klo:kup).GT.1)
   vsign = SIGN(1.0,vadvatvw(:,2-nh:nrloc+nh,klo:kup))
END WHERE

!
!3.5 CFL number
!--------------
!

IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN
   WHERE (mask2dvw)
      array2dv = delt3d/delyatc(1:ncloc,2-nh:nrloc+nh)
   END WHERE
   k_320: DO k=klo,kup
      WHERE (maskatvw(:,:,k))
         cfl(:,:,k) = array2dv*vadvatvw(:,2-nh:nrloc+nh,k)
      END WHERE
   ENDDO k_320
ENDIF

!
!3.6 Work space arrays
!---------------------
!

IF (dXregY) THEN
   k_361: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = delyatc(1:ncloc,1:nrloc)*delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_361
ELSE
   k_362: DO k=klo,kup
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_362
ENDIF

!
!4. Explicit terms
!-----------------
!
!4.1 Upwind fluxes
!-----------------
!

IF (iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3) THEN
   
   k_410: DO k=klo,kup
      WHERE (maskatvw(:,:,k))
         advflux_up(:,:,k) = vadvatvw(:,2-nh:nrloc+nh,k)*&
                          & ((1.0+vsign(:,:,k))*psiw(1:ncloc,1-nh:nrloc+nh-1,k)&
                          & +(1.0-vsign(:,:,k))*psiw(1:ncloc,2-nh:nrloc+nh,k))
      END WHERE
   ENDDO k_410

!  ---non-TVD case
   IF (iopt_adv_turb.EQ.1) THEN
      WHERE (maskatvw(:,1:nrloc+1,:))
         advflux = array3dvw1*advflux_up(:,1:nrloc+1,:)
      END WHERE
   ENDIF

ENDIF

!
!4.2 Lax-Wendroff fluxes
!-----------------------
!

IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN

   k_420: DO k=klo,kup
      WHERE (maskatvw(:,:,k))
         advflux_lw(:,:,k) = vadvatvw(:,2-nh:nrloc+nh,k)*&
                           & ((1.0+cfl(:,:,k))*psiw(1:ncloc,1-nh:nrloc+nh-1,k)&
                           & +(1.0-cfl(:,:,k))*psiw(1:ncloc,2-nh:nrloc+nh,k))
      END WHERE
   ENDDO k_420

!  ---non-TVD case
   IF (iopt_adv_turb.EQ.2) THEN
      WHERE (maskatvw(:,1:nrloc+1,:))
         advflux = array3dvw1*advflux_lw(:,1:nrloc+1,:)
      END WHERE
   ENDIF

ENDIF

!
!4.3 TVD scheme
!-------------
!

IF (iopt_adv_turb.EQ.3) THEN

!  ---flux differences
   WHERE (maskatvw)
      dflux = advflux_lw - advflux_up
   END WHERE

!  ---flux ratios
   WHERE (maskatvw(:,1:nrloc+1,:)) array3dvw2 = 0.0
   k_430: DO k=klo,kup
      WHERE (ABS(dflux(:,1:nrloc+1,k)).GT.eps_adv)
         array3dvw2(:,:,k) = 0.5*(&
            & (1.0+vsign(:,1:nrloc+1,k))*dflux(:,0:nrloc,k)+&
            & (1.0-vsign(:,1:nrloc+1,k))*dflux(:,2:nrloc+2,k))&
            & /dflux(:,1:nrloc+1,k)
      END WHERE
   ENDDO k_430

!  ---TVD limiter
   omega = tvd_limiter(array3dvw2,maskatvw(:,1:nrloc+1,:))

!  ---fluxes
   WHERE (maskatvw(:,1:nrloc+1,:))
      advflux = array3dvw1*(advflux_up(:,1:nrloc+1,:)+&
                          & omega*dflux(:,1:nrloc+1,:))
   END WHERE

ENDIF

!
!4.4 Open boundary fluxes
!------------------------
!

jjloc_440: DO jjloc=1,nobvloc_ext
   i = iobvloc(jjloc); j = jobvloc(jjloc)
   jj = indexobv(jjloc); jc = MERGE(j,j-1,soutobv(jj))
   IF (jc.GE.1.AND.jc.LE.nrloc) THEN
      j1 = MERGE(j+1,j-1,soutobv(jj))
      advflux(i,j,:) = advflux(i,j1,:)
   ENDIF
ENDDO jjloc_440

!
!4.5 Advective term
!-------------------
!

k_450: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfw(:,:,k,4) = tridcfw(:,:,k,4) - &
                       & (advflux(:,2:nrloc+1,k)-advflux(:,1:nrloc,k))&
                       & /denom(:,:,k)
   END WHERE
ENDDO k_450

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatvw,mask2dvw)
DEALLOCATE (advflux,array2dv,array3dvw1,denom,vsign)
IF (iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3) DEALLOCATE (advflux_up)
IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) DEALLOCATE (advflux_lw,cfl)
IF (iopt_adv_turb.EQ.3) DEALLOCATE (array3dvw2,dflux,omega)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Yadv_at_W

!========================================================================

SUBROUTINE Zadv_at_W(psiw,tridcfw,wadvatc,klo,kup)
!************************************************************************
!
! *Zadv_at_W* Vertical advection term for a quantity at W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Advection_Terms.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - transport_at_W
!
! Module calls - error_alloc, tvd_limiter
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
USE utility_routines, ONLY: tvd_limiter

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: klo, kup
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1)::psiw
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz+1,4) :: tridcfw
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz) :: wadvatc

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*psiw*      REAL    W-node quantity to be advected                     [psic]
!*tridcfw*   REAL    Tridiagonal matrix for implicit vertical solution
!*wadvatw*   REAL    Advective velocity at W-nodes                       [m/s]
!*klo*       INTEGER Lower vertical array bound
!*kup*       INTEGER Upper vertical array bound
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: iopt_impl, k, npcc
REAL :: xexp, ximp
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, cfl
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advflux, advflux_ce, advflux_lw, &
                            & advflux_up, array3dw1, array3dw2, dflux, omega, &
                            & upwind,wsign

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*advflux*    REAL    Advective flux (times factor) at C-nodes          [psiw]
!*advflux_ce* REAL    Central flux at C-nodes                     [psiw]*[m/s]
!*advflux_lw* REAL    Lax_Wendroff flux at C-nodes                [psiw]*[m/s]
!*advflux_up* REAL    Upwind flux at C-nodes                      [psiw]*[m/s]
!*omega*      REAL    TVD limiting factor
!*wsign*      REAL    Upwind factor
!*wadvatc*    REAL    Advective velocity at C-nodes                      [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'Zadv_at_w'
CALL log_timer_in(npcc)

!
!1. Set parameters
!-----------------
!

IF (iopt_adv_turb.EQ.2) THEN
   xexp = 0.5*delt3d
   ximp = 0.0
   iopt_impl = 0
ELSE
   xexp = 0.5*delt3d*(1.0-theta_vadv)
   ximp = 0.5*delt3d*theta_vadv
   iopt_impl = iopt_vadv_impl
ENDIF

!
!2. Allocate arrays
!------------------
!

IF (iopt_impl.NE.2) THEN
   ALLOCATE (advflux(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('advflux',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

IF (iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3) THEN
   ALLOCATE (wsign(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('wsign',3,(/ncloc,nrloc,nz/),kndrtype)
   IF (iopt_impl.NE.2) THEN
      ALLOCATE (advflux_up(ncloc,nrloc,0:nz+1),STAT=errstat)
      CALL error_alloc('advflux_up',3,(/ncloc,nrloc,nz+2/),kndrtype)
   ENDIF
ENDIF

IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      ALLOCATE (advflux_ce(ncloc,nrloc,0:nz+1),STAT=errstat)
      CALL error_alloc('advflux_ce',3,(/ncloc,nrloc,nz+2/),kndrtype)
   ELSE
      ALLOCATE (advflux_lw(ncloc,nrloc,0:nz+1),STAT=errstat)
      CALL error_alloc('advflux_lw',2,(/ncloc,nrloc,nz+2/),kndrtype)
      ALLOCATE (cfl(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('cfl',2,(/ncloc,nrloc/),kndrtype)
   ENDIF
ENDIF

IF (iopt_adv_turb.EQ.3) THEN
   ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc1',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (dflux(ncloc,nrloc,0:nz+1),STAT=errstat)
   CALL error_alloc('dflux',2,(/ncloc,nrloc,nz+2/),kndrtype)
   ALLOCATE (omega(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('omega',2,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

IF (iopt_impl.GT.0) THEN
   ALLOCATE (array3dw1(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('array3dw1',2,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (array3dw2(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('array2dw2',2,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (upwind(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('upwind',2,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

!
!3. Initialise arrays
!--------------------
!
!3.1 Fluxes
!----------
!

IF (iopt_impl.NE.2.AND.(iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3)) THEN
   advflux_up = 0.0
ENDIF

IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      advflux_ce = 0.0
   ELSE
      advflux_lw = 0.0
   ENDIF
ENDIF

IF (iopt_adv_turb.EQ.3) dflux = 0.0

!
!3.2 Work space arrays
!---------------------
!

IF (iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3) THEN
   k_321: DO k=klo-1,kup
      WHERE (maskatc_int)
         wsign(:,:,k) = SIGN(1.0,wadvatc(:,:,k))
      END WHERE
   ENDDO k_321
ENDIF

!
!4. Explicit terms
!-----------------
!
!4.1 Upwind fluxes
!-----------------
!

IF ((iopt_adv_turb.EQ.1.AND.iopt_impl.NE.2).OR.iopt_adv_turb.EQ.3) THEN

   k_411: DO k=klo-1,kup
      WHERE (maskatc_int)
         advflux_up(:,:,k) = wadvatc(:,:,k)*&
                           & ((1.0+wsign(:,:,k))*psiw(1:ncloc,1:nrloc,k)&
                           & +(1.0-wsign(:,:,k))*psiw(1:ncloc,1:nrloc,k+1))
      END WHERE
   ENDDO k_411

!  ---non-TVD case
   IF (iopt_adv_turb.EQ.1) THEN
      k_412: DO k=klo-1,kup
         WHERE (maskatc_int)
            advflux(:,:,k) = xexp*advflux_up(:,:,k)
         END WHERE
      ENDDO k_412
   ENDIF

ENDIF

!
!4.2 Central fluxes
!------------------
!

IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN

   IF (iopt_impl.GT.0) THEN

      k_421: DO k=klo-1,kup
         WHERE (maskatc_int)
            advflux_ce(:,:,k) = wadvatc(:,:,k)*(psiw(1:ncloc,1:nrloc,k)+&
                                              & psiw(1:ncloc,1:nrloc,k+1))
         END WHERE
      ENDDO k_421

!     ---non-TVD case
      IF (iopt_adv_turb.EQ.2.AND.iopt_impl.NE.2) THEN
         k_422: DO k=klo-1,kup
            WHERE (maskatc_int)
               advflux(:,:,k) = xexp*advflux_ce(:,:,k)
            END WHERE
         ENDDO k_422
      ENDIF

!
!4.3 Lax-Wendroff fluxes
!-----------------------
!

   ELSE

      k_431: DO k=klo-1,kup
         WHERE (maskatc_int)
            cfl = delt3d*wadvatc(:,:,k)/delzatc(1:ncloc,1:nrloc,k)
            advflux_lw(:,:,k) = wadvatc(:,:,k)*&
                              & ((1.0+cfl)*psiw(1:ncloc,1:nrloc,k)&
                              & +(1.0-cfl)*psiw(1:ncloc,1:nrloc,k+1))
         END WHERE
      ENDDO k_431
      
!     ---non-TVD case
      IF (iopt_adv_turb.EQ.2) THEN
         k_332: DO k=klo-1,kup
            WHERE (maskatc_int)
               advflux(:,:,k) = xexp*advflux_lw(:,:,k)
            END WHERE
         ENDDO k_332
      ENDIF
   ENDIF

ENDIF

!
!4.4 TVD scheme
!--------------
!

IF (iopt_adv_turb.EQ.3) THEN

!  ---flux differences
   k_441: DO k=klo-1,kup
      IF (iopt_impl.NE.0) THEN
         WHERE (maskatc_int)
            dflux(:,:,k) = advflux_ce(:,:,k) - advflux_up(:,:,k)
         END WHERE
      ELSE
         WHERE (maskatc_int)
            dflux(:,:,k) = advflux_lw(:,:,k) - advflux_up(:,:,k)
         END WHERE
      ENDIF
   ENDDO k_441

   k_442: DO k=klo-1,kup

!     --flux ratios
      WHERE (maskatc_int) array2dc = 0.0
      WHERE (ABS(dflux(:,:,k)).GT.eps_adv)
         array2dc = 0.5*((1.0+wsign(:,:,k))*dflux(:,:,k-1)+&
                       & (1.0-wsign(:,:,k))*dflux(:,:,k+1))/dflux(:,:,k)
      END WHERE

!     --TVD limiter
      omega(:,:,k) = tvd_limiter(array2dc,maskatc_int)

!     --fluxes
      IF (iopt_impl.NE.2) THEN
         WHERE (maskatc_int)
            advflux(:,:,k) = xexp*(advflux_up(:,:,k)+omega(:,:,k)*dflux(:,:,k))
         END WHERE
      ENDIF

   ENDDO k_442

ENDIF

!
!4.5 Advective term
!------------------
!

IF (iopt_impl.NE.2) THEN
   k_451: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfw(:,:,k,4) = tridcfw(:,:,k,4) + &
                          & (advflux(:,:,k)-advflux(:,:,k-1))&
                          & /delzatw(1:ncloc,1:nrloc,k)
         END WHERE
      ENDDO k_451
ENDIF

!
!5. Implicit terms
!-----------------
!

IF (iopt_impl.GT.0) THEN

!
!5.1 Upwind factor
!-----------------
!

   k_511: DO k=klo-1,kup
      IF (iopt_adv_turb.EQ.1) THEN
         WHERE (maskatc_int) upwind(:,:,k) = wsign(:,:,k)
      ELSEIF (iopt_adv_turb.EQ.2) THEN
         WHERE (maskatc_int) upwind(:,:,k) = 0.0
      ELSEIF (iopt_adv_turb.EQ.3) THEN
         WHERE (maskatc_int) upwind(:,:,k) = (1.0-omega(:,:,k))*wsign(:,:,k)
      ENDIF
   ENDDO k_511
   
!
!5.2 Downward terms
!------------------
!

   WHERE (maskatc_int)
      array3dw2(:,:,klo-1) = wadvatc(:,:,klo-1)*(1.0+upwind(:,:,klo-1))
   END WHERE
   k_521: DO k=klo,kup
      WHERE (maskatc_int)
         array3dw1(:,:,k) = ximp/delzatw(1:ncloc,1:nrloc,k)
         array3dw2(:,:,k) = wadvatc(:,:,k)*(1.0+upwind(:,:,k))
      END WHERE
   ENDDO k_521
   k_522: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfw(:,:,k,1) = tridcfw(:,:,k,1) - &
                          & array3dw1(:,:,k)*array3dw2(:,:,k-1)
         tridcfw(:,:,k,2) = tridcfw(:,:,k,2) + &
                          & array3dw1(:,:,k)*array3dw2(:,:,k)
      END WHERE
   ENDDO k_522

!
!5.3 Upward terms
!----------------
!

   k_531: DO k=klo-1,kup
      WHERE (maskatc_int)
         array3dw2(:,:,k) = wadvatc(:,:,k)*(1.0-upwind(:,:,k))
      END WHERE
   ENDDO k_531
   k_532: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfw(:,:,k,2) = tridcfw(:,:,k,2) - &
                          & array3dw1(:,:,k)*array3dw2(:,:,k-1)
         tridcfw(:,:,k,3) = tridcfw(:,:,k,3) + &
                          & array3dw1(:,:,k)*array3dw2(:,:,k)
      END WHERE
   ENDDO k_532

ENDIF

!
!6. Deallocate arrays
!--------------------
!

IF (iopt_impl.NE.2) DEALLOCATE (advflux)
IF (iopt_adv_turb.EQ.1.OR.iopt_adv_turb.EQ.3) THEN
   DEALLOCATE (wsign)
   IF (iopt_impl.NE.2) DEALLOCATE (advflux_up)
ENDIF
IF (iopt_adv_turb.EQ.2.OR.iopt_adv_turb.EQ.3) THEN
   IF (iopt_impl.GT.0) THEN
      DEALLOCATE (advflux_ce)
   ELSE
      DEALLOCATE (advflux_lw,cfl)
   ENDIF
ENDIF
IF (iopt_adv_turb.EQ.3) DEALLOCATE (array2dc,dflux,omega)
IF (iopt_impl.GT.0) DEALLOCATE (array3dw1,array3dw2,upwind)

CALL log_timer_out(npcc,itm_adv)


RETURN

END SUBROUTINE Zadv_at_W
