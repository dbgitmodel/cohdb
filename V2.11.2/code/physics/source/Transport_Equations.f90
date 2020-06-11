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
! *Transport_Equations* Solve transport equations at various nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Transport_Equations.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
! Reference -
!
! Routines - transport_at_C_3d, transport_at_C_4d1, transport_at_C_4d2,
!            transport_at_U_3d, transport_at_U_2d, transport_at_V_3d,
!            transport_at_V_2d, transport_at_W
!
!************************************************************************
!

!========================================================================

SUBROUTINE transport_at_C_3d(psic,source,vdifcoefatw,iopt_adv,iopt_hdif,psiobu,&
                           & psiobv,iprofrlx,ibcsur,ibcbot,nbcs,nbcb,bcsur,&
                           & bcbot,ivarid)
!************************************************************************
!
! *transport_at_C_3d* Solve advection-diffusion equation for a 3-D quantity at
!                     C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Transport_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - salinity_equation, temperature_equation
!
! External calls - relaxation_at_C, Xadv_at_C, Xcorr_at_C, Xdif_at_C,
!                  Yadv_at_C, Ycorr_at_C, Ydif_at_C, Zadv_at_C, Zcorr_at_C,
!                  Zdif_at_C
!
! Module calls - error_alloc, exchange_mod, num_halo, tridiag_vert
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE physpars
USE relaxation
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE nla_library, ONLY: tridiag_vert
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: num_halo

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ibcbot, ibcsur, iopt_adv, iopt_hdif, ivarid, nbcs, nbcb
INTEGER, INTENT(IN), DIMENSION(norlxzones) :: iprofrlx
REAL, INTENT(IN), DIMENSION(nobu,nz) :: psiobu
REAL, INTENT(IN), DIMENSION(nobv,nz) :: psiobv
REAL, INTENT(INOUT), DIMENSION(1-nhalo:ncloc+nhalo,&
                             & 1-nhalo:nrloc+nhalo,nz) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz) :: source
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1) :: vdifcoefatw
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcs) :: bcsur
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcb) :: bcbot

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psic*        REAL    C-node scalar to be updated                        [psic]
!*source*      REAL    Source terms in transport equations              [psic/s]
!*vdifcoefatw* REAL    Vertical diffusion coefficient at W-nodes         [m^2/s]
!*iopt_adv*    INTEGER Switch to select advection scheme (disabled if zero)
!*iopt_hdif*   INTEGER Switch to select horizontal diffusion scheme
!                      (disabled if zero)
!*psiobu*      REAL    Vertical profiles at U-open boundaries
!                      (zero gradient condition if flagged)               [psic]
!*psiobv*      REAL    Vertical profiles at V-open boundaries
!                      (zero gradient condition if flagged)               [psic]
!*iprofrlx*    INTEGER Disables/enables relaxation at open boundary zones
!*ibcsur*      INTEGER Type of surface boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!                 = 3 => Dirichlet at first C-node below the surface
!                 = 4 => Dirichlet at the surface
!*ibcbot*      INTEGER Type of bottom boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!                 = 3 => Dirichlet at first C-node above the bottom
!                 = 4 => Dirichlet at the bottom
!*nbcs*        INTEGER Last dimension of array bcsur
!*nbcb*        INTEGER Last dimension of array bcbot
!*bcsur*       REAL    Data for surface boundary condition
!              (:,:,1) => prescribed surface flux or surface value
!              (:,:,2) => transfer velocity                                [m/s]
!*bcbot*       REAL    Data for bottom boundary condition
!              (:,:,1) => prescribed surface flux or surface value
!              (:,:,2) => transfer velocity                                [m/s]
!*ivarid*      INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: time_split, vadv, vimplicit
INTEGER :: k, klo, kup, nhpsi, npcc
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhdims, nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: psic_A, psic_B, ucorr, vcorr, wcorr
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*time_split* LOGICAL .TRUE. for directional splitting in time
!*vimplicit*  LOGICAL .TRUE. for implicit integration in the vertical
!*tridcfc*    REAL    Tridiagonal matrix for implicit vertical solution
!*psic_A*     REAL    psic during first integration step            [psic]
!*psic_B*     REAL    psic during second integration step           [psic]  
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'transport_at_C_3d'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Initialise
!-------------
!
!1.1 Array shape and halo size
!-----------------------------
!

nhpsi = MAX(num_halo(iopt_adv),MIN(1,iopt_hdif))
nhdims = nhalo
klo = MERGE(1,2,ibcbot.LE.2)
kup = MERGE(nz,nz-1,ibcsur.LE.2)

!
!1.2 Time integration
!--------------------
!

vimplicit = ((iopt_adv.GT.0.AND.iopt_vadv_impl.GT.0).OR.&
           & (iopt_vdif_coef.GT.0.AND.iopt_vdif_impl.GT.0).OR.&
           & (ibcsur.EQ.4).OR.(ibcbot.EQ.4))
vadv = iopt_grid_nodim.NE.2.AND.iopt_adv.GT.0
time_split = iopt_adv.EQ.3

!
!1.3 Allocate arrays
!-------------------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (tridcfc(ncloc,nrloc,nz,4),STAT=errstat)
CALL error_alloc('tridcfc',4,(/ncloc,nrloc,nz,4/),kndrtype)
IF (time_split) THEN
   ALLOCATE (psic_A(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),&
        & STAT=errstat)
   CALL error_alloc('psic_A',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
   ALLOCATE (psic_B(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),&
        & STAT=errstat)
   CALL error_alloc('psic_B',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
ENDIF

IF (time_split) THEN
   ALLOCATE (ucorr(ncloc,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('ucorr',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (vcorr(ncloc,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('vcorr',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (wcorr(ncloc,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('wcorr',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
ENDIF

!
!1.4 Dirichlet boundary conditions
!---------------------------------
!
!---surface
IF (ibcsur.EQ.3) THEN
   WHERE (maskatc_int)
      tridcfc(:,:,nz,1) = 0.0
      tridcfc(:,:,nz,2) = 1.0
      tridcfc(:,:,nz,3) = 0.0
      tridcfc(:,:,nz,4) = bcsur(:,:,1)
   END WHERE
ELSEIF (ibcsur.EQ.4) THEN
   WHERE (maskatc_int)
      array2dc = delzatc(1:ncloc,1:nrloc,nz) + &
               & 2.0*delzatw(1:ncloc,1:nrloc,nz)
      tridcfc(:,:,nz,1) = -delzatc(1:ncloc,1:nrloc,nz)/array2dc
      tridcfc(:,:,nz,2) = 1.0
      tridcfc(:,:,nz,3) = 0.0
      tridcfc(:,:,nz,4) = 2.0*delzatw(1:ncloc,1:nrloc,nz)*bcsur(:,:,1)&
                        & /array2dc
   END WHERE
ENDIF

!---bottom
IF (ibcbot.EQ.3) THEN
   WHERE (maskatc_int)
      tridcfc(:,:,1,1) = 0.0
      tridcfc(:,:,1,2) = 1.0
      tridcfc(:,:,1,3) = 0.0
      tridcfc(:,:,1,4) = bcbot(:,:,1)
   END WHERE
ELSEIF (ibcbot.EQ.4) THEN
   WHERE (maskatc_int)
      array2dc = 2.0*delzatw(1:ncloc,1:nrloc,2) + delzatc(1:ncloc,1:nrloc,1)
      tridcfc(:,:,1,1) = 0.0
      tridcfc(:,:,1,2) = 1.0
      tridcfc(:,:,1,3) = -delzatc(1:ncloc,1:nrloc,1)/array2dc
      tridcfc(:,:,1,4) = 2.0*delzatw(1:ncloc,1:nrloc,2)*bcbot(:,:,1)/array2dc
   END WHERE
ENDIF

!
!1.5 Work space arrays
!---------------------
!

IF (.NOT.time_split) THEN
   WHERE (maskatc_int)
      array2dc = deptotatc_old(1:ncloc,1:nrloc)/deptotatc(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!1.6 Local arrays
!----------------
!

IF (time_split) THEN
   psic_A = psic; psic_B = psic
ENDIF

!
!1.7 Source terms
!----------------
!

k_170:DO k=klo,kup
   WHERE (maskatc_int)
      source(:,:,k) = delt3d*source(:,:,k)
   END WHERE 
ENDDO k_170

!
!2. Corrector terms
!------------------
!

IF (time_split) THEN
   CALL Xcorr_at_C(ucorr,klo,kup)
   CALL Ycorr_at_C(vcorr,klo,kup)
   CALL Zcorr_at_C(wcorr,klo,kup)
ENDIF

!
!3. Time integration without operator splitting
!----------------------------------------------
!

IF (.NOT.time_split) THEN

!
!3.1 Initialise
!--------------
!
   
   IF (vimplicit) THEN
      k_310: DO k=klo,kup
         WHERE (maskatc_int)
            tridcfc(:,:,k,1) = 0.0
            tridcfc(:,:,k,3) = 0.0
         END WHERE
      ENDDO k_310
   ENDIF

!
!3.2 Time derivative
!-------------------
!

   k_320: DO k=klo,kup
      WHERE (maskatc_int) 
         tridcfc(:,:,k,2) = 1.0
         tridcfc(:,:,k,4) = array2dc*psic(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_320

!
!3.3 Source term
!---------------
!

   k_330: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4) = tridcfc(:,:,k,4) + source(:,:,k)
      END WHERE
   ENDDO k_330

!
!3.4 Advection
!-------------
!
!  ---horizontal
   IF (iopt_adv.GT.0) THEN
      CALL Xadv_at_C(psic,tridcfc,1,iopt_adv,nhpsi,klo,kup,psiobu,ufvel)
      CALL Yadv_at_C(psic,tridcfc,1,iopt_adv,nhpsi,klo,kup,psiobv,vfvel)
   ENDIF
!  ---vertical
   IF (vadv) THEN
      CALL Zadv_at_C(psic,tridcfc,wfvel(1:ncloc,1:nrloc,:),1,iopt_adv,klo,kup)
   ENDIF

!
!3.5 Diffusion
!-------------
!
!  ---horizontal
   IF (iopt_hdif.GT.0) THEN
      CALL Xdif_at_C(psic,tridcfc,1,klo,kup)
      CALL Ydif_at_C(psic,tridcfc,1,klo,kup)
   ENDIF
!  ---vertical
   CALL Zdif_at_C(psic,tridcfc,vdifcoefatw,1,ibcsur,ibcbot,nbcs,nbcb,bcsur,&
                & bcbot,klo,kup)

!
!3.6 Update psic
!---------------
!
   IF (vimplicit) THEN
      CALL tridiag_vert(tridcfc,psic,nhdims,nz,1,'C  ')
   ELSE
      k_360: DO k=1,nz
         WHERE (maskatc_int)
            psic(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)/tridcfc(:,:,k,2)
         END WHERE
      ENDDO k_360
   ENDIF

   GOTO 1000

ENDIF

!
!4. Time integration with operator splitting - step A
!----------------------------------------------------
!
!4.1 X-derivative terms
!----------------------
!
!---time derivative, corrector term
k_411: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,2) = 1.0
      tridcfc(:,:,k,4) = psic(1:ncloc,1:nrloc,k)*(1.0+ucorr(:,:,k))
   END WHERE
ENDDO k_411

!---horizontal advection
IF (iopt_adv.GT.0) THEN
   CALL Xadv_at_C(psic,tridcfc,1,iopt_adv,nhpsi,klo,kup,psiobu,ufvel)
ENDIF

!---horizontal diffusion
IF (iopt_hdif.GT.0) THEN
   CALL Xdif_at_C(psic,tridcfc,1,klo,kup)
ENDIF

!---update psic_A
k_412: DO k=klo,kup
   WHERE (maskatc_int)
      psic_A(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)/tridcfc(:,:,k,2)
   END WHERE
ENDDO k_412

!
!4.2 Y-derivative terms
!----------------------
!
!---time derivative, corrector term
k_421: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,2) = 1.0
      tridcfc(:,:,k,4) = psic_A(1:ncloc,1:nrloc,k)*(1.0+vcorr(:,:,k))
   END WHERE
ENDDO k_421

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   lbounds = (/1-nhalo,1-nhalo,1/); nhexch = (/0,0,nhpsi,nhpsi/)
   CALL exchange_mod(psic_A,lbounds,nhexch,0)
ENDIF

!---horizontal advection
IF (iopt_adv.GT.0) THEN
   CALL Yadv_at_C(psic_A,tridcfc,1,iopt_adv,nhpsi,klo,kup,psiobv,vfvel)
ENDIF

!---horizontal diffusion
IF (iopt_hdif.GT.0) THEN
   CALL Ydif_at_C(psic_A,tridcfc,1,klo,kup)
ENDIF

!---update psic_A
k_422: DO k=klo,kup
   WHERE (maskatc_int)
      psic_A(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)/tridcfc(:,:,k,2)
   END WHERE
ENDDO k_422

!---exchange halo sections (if needed)
IF (iopt_MPI.EQ.1.AND.iopt_hdif_scal.GT.1) THEN
   lbounds = (/1-nhalo,1-nhalo,1/); nhexch = nhpsi
   CALL exchange_mod(psic_A,lbounds,nhexch,0,corners=.FALSE.)
ENDIF

!
!4.3 Z-derivative/source terms
!-----------------------------
!
!4.3.1 Time derivative, source and corrector terms
!-------------------------------------------------
!

IF (vimplicit) THEN
   k_4311: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,1) = 0.0
         tridcfc(:,:,k,3) = 0.0
      END WHERE
   ENDDO k_4311
ENDIF
k_4322: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,4) = psic_A(1:ncloc,1:nrloc,k)*&
                       & (1.0+(1.0-theta_vadv)*wcorr(:,:,k)) + source(:,:,k)
      tridcfc(:,:,k,2) = 1.0 - theta_vadv*wcorr(:,:,k)
   END WHERE
ENDDO k_4322

!
!4.3.2 Vertical advection
!------------------------
!

IF (vadv) THEN
   CALL Zadv_at_C(psic_A,tridcfc,wfvel(1:ncloc,1:nrloc,:),1,iopt_adv,klo,kup)
ENDIF

!
!4.3.3 Vertical diffusion
!------------------------
!

CALL Zdif_at_C(psic_A,tridcfc,vdifcoefatw,1,ibcsur,ibcbot,nbcs,nbcb,bcsur,&
               bcbot,klo,kup)

!
!4.3.4 Update psic_A
!-------------------
!

IF (vimplicit) THEN
   CALL tridiag_vert(tridcfc,psic_A,nhdims,nz,1,'C  ')
ELSE
   k_434: DO k=1,nz
      WHERE (maskatc_int)
         psic_A(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)/tridcfc(:,:,k,2)
      END WHERE
   ENDDO k_434
ENDIF

!
!5. Time integration with operator splitting - step B
!----------------------------------------------------
!
!5.1 Z-derivative/source terms
!-----------------------------
!
!5.1.1 Time derivative, source and corrector terms
!-------------------------------------------------
!

IF (vimplicit) THEN
   k_5111: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,1) = 0.0
         tridcfc(:,:,k,3) = 0.0
      END WHERE
   ENDDO k_5111
ENDIF
k_5112: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,4) = psic(1:ncloc,1:nrloc,k)*&
                       & (1.0+(1.0-theta_vadv)*wcorr(:,:,k)) + source(:,:,k)
      tridcfc(:,:,k,2) = 1.0 - theta_vadv*wcorr(:,:,k)
      
   END WHERE
ENDDO k_5112

!
!5.1.2 Vertical advection
!------------------------
!

IF (vadv) THEN
   CALL Zadv_at_C(psic,tridcfc,wfvel(1:ncloc,1:nrloc,:),1,iopt_adv,klo,kup)
ENDIF

!
!5.1.3 Vertical diffusion
!------------------------
!

CALL Zdif_at_C(psic,tridcfc,vdifcoefatw,1,ibcsur,ibcbot,nbcs,nbcb,bcsur,&
             & bcbot,klo,kup)

!
!5.1.4 Update psic_B
!-------------------
!

IF (vimplicit) THEN
   CALL tridiag_vert(tridcfc,psic_B,nhdims,nz,1,'C  ')
ELSE
   k_514: DO k=1,nz
      WHERE (maskatc_int)
         psic_B(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)
      END WHERE
   ENDDO k_514
ENDIF

!
!5.2 Y-derivative terms
!----------------------
!
!---time derivative, corrector term
k_521: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,2) = 1.0
      tridcfc(:,:,k,4) = psic_B(1:ncloc,1:nrloc,k)*(1.0+vcorr(:,:,k))
   END WHERE
ENDDO k_521

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/0,0,nhpsi,nhpsi/)
   CALL exchange_mod(psic_B,lbounds,nhexch,0)
ENDIF

!---horizontal advection
IF (iopt_adv.GT.0) THEN
   CALL Yadv_at_C(psic_B,tridcfc,1,iopt_adv,nhpsi,klo,kup,psiobv,vfvel)
ENDIF

!---horizontal diffusion
IF (iopt_hdif.GT.0) THEN
   CALL Ydif_at_C(psic_B,tridcfc,1,klo,kup)
ENDIF

!---update psic_B
k_522: DO k=klo,kup
   WHERE (maskatc_int)
      psic_B(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)/tridcfc(:,:,k,2)
   END WHERE
ENDDO k_522

!
!5.3 X-derivative terms
!----------------------
!
!---time derivative, corrector term
k_531: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,2) = 1.0
      tridcfc(:,:,k,4) = psic_B(1:ncloc,1:nrloc,k)*(1.0+ucorr(:,:,k))
   END WHERE
ENDDO k_531

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/nhpsi,nhpsi,0,0/)
   CALL exchange_mod(psic_B,lbounds,nhexch,0)
ENDIF

!---horizontal advection
IF (iopt_adv.GT.0) THEN
   CALL Xadv_at_C(psic_B,tridcfc,1,iopt_adv,nhpsi,klo,kup,psiobu,ufvel)
ENDIF

!---horizontal diffusion
IF (iopt_hdif.GT.0) THEN
   CALL Xdif_at_C(psic_B,tridcfc,1,klo,kup)
ENDIF

!---update psic_B
k_532: DO k=klo,kup
   WHERE (maskatc_int)
      psic_B(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)/tridcfc(:,:,k,2)
   END WHERE
ENDDO k_532

!
!6. Time integration with operator splitting - update psic
!---------------------------------------------------------
!

k_610: DO k=1,nz
   WHERE (maskatc_int)
      psic(1:ncloc,1:nrloc,k) = 0.5*(psic_A(1:ncloc,1:nrloc,k)+&
                                   & psic_B(1:ncloc,1:nrloc,k))
   END WHERE
ENDDO k_610

1000 CONTINUE

!
!7. Apply relaxation scheme
!--------------------------
!

IF (iopt_obc_relax.EQ.1.AND.inodesrlx(1).GT.0) THEN
   CALL relaxation_at_C(psic,psiobu,psiobv,1,iprofrlx)
ENDIF

!
!8. Deallocate arrays
!--------------------
!

DEALLOCATE (array2dc,tridcfc)
IF (time_split) DEALLOCATE (psic_A,psic_B,ucorr,vcorr,wcorr)

CALL log_timer_out(npcc,itm_trans)


RETURN

END SUBROUTINE transport_at_C_3d

!========================================================================

SUBROUTINE transport_at_C_4d1(psic,source,wsink,vdifcoefatw,novars,iopt_hadv,&
                            & iopt_vadv,iopt_hdif,psiobu,psiobv,iprofrlx,&
                            & ibcsur,ibcbot,nbcs,nbcb,bcsur,bcbot,settling,&
                            & ivarids)
!************************************************************************
!
! *transport_at_C_4d1* Solve advection-diffusion equation for a series of 3-D
!                      quantities at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Transport_Equations.f90  V2.11.2
!
! Description - each equation is separately updated
!
! Reference -
!
! Calling program -
!
! External calls - relaxation_at_C, Xadv_at_C, Xcorr_at_C, Xdif_at_C,
!                  Yadv_at_C, Ycorr_at_C, Ydif_at_C, Zadv_at_C, Zcorr_at_C,
!                  Zdif_at_C
!
! Module calls - error_alloc, exchange_mod, inquire_var, num_halo,
!                tridiag_vert, transf_vertical_fall
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE physpars
USE relaxation
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE modvars_routines, ONLY: inquire_var
USE nla_library, ONLY: tridiag_vert
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: num_halo

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: settling
INTEGER, INTENT(IN) :: ibcbot, ibcsur, iopt_hadv, iopt_hdif, iopt_vadv, nbcb, &
                     & nbcs, novars
INTEGER, INTENT(IN), DIMENSION(novars) :: ivarids
INTEGER, INTENT(IN), DIMENSION(norlxzones) :: iprofrlx
REAL, INTENT(IN), DIMENSION(nobu,nz,novars) :: psiobu
REAL, INTENT(IN), DIMENSION(nobv,nz,novars) :: psiobv
REAL, INTENT(INOUT), DIMENSION(1-nhalo:ncloc+nhalo,&
                             & 1-nhalo:nrloc+nhalo,nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,novars) :: source
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1,novars) :: vdifcoefatw
REAL, INTENT(IN), DIMENSION(0:ncloc,0:nrloc,nz+1,novars) :: wsink
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcs,novars) :: bcsur
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcb,novars) :: bcbot

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psic*        REAL    C-node scalar to be updated                       [psic]
!*source*      REAL    Source terms in transport equations              [psic/s]
!*wsink*       REAL    Vertical (positive downward) sinking velocities at
!                      W-nodes if argument settling equals .TRUE.,
!                      not used otherwise                                  [m/s]
!*vdifcoefatw* REAL    Vertical diffusion coefficients at W-nodes        [m^2/s]
!*novars*      INTEGER Number of variables
!*iopt_hadv*   INTEGER Switch to select horizontal advection scheme
!                      (disabled if zero)
!*iopt_vadv*   INTEGER Switch to select vertical advection scheme
!                      (disabled if zero, must be positive if settling is
!                      .TRUE.)
!*iopt_hdif*   INTEGER Switch to select horizontal diffusion scheme
!                      (disabled if zero)
!*psiobu*      REAL    Vertical profiles at U-open boundaries
!                      (zero gradient condition if flagged)               [psic]
!*psiobv*      REAL    Vertical profiles at V-open boundaries
!                      (zero gradient condition if flagged)               [psic]
!*iprofrlx*    INTEGER Disables/enables relaxation at open boundary zones
!*ibcsur*      INTEGER Type of surface boundary condition
!                  = 0 => Neumann (zero flux)
!                  = 1 => Neumann (prescibed flux)
!                  = 2 => Neumann (using transfer velocity)
!                  = 3 => Dirichlet at first C-node below the surface
!                  = 4 => Dirichlet at the surface
!*ibcbot*      INTEGER Type of bottom boundary condition
!                  = 0 => Neumann (zero flux)
!                  = 1 => Neumann (prescibed flux)
!                  = 2 => Neumann (using transfer velocity)
!                  = 3 => Dirichlet at first C-node above the bottom
!                  = 4 => Dirichlet at the bottom
!*nbcs*        INTEGER Last dimension of array bcsur
!*nbcb*        INTEGER Last dimension of array bcbot
!*bcsur*       REAL    Data for surface boundary condition
!              (:,:,1,:) => prescribed surface flux or surface value
!              (:,:,2,:) => transfer velocity                              [m/s]
!*bcbot*       REAL    Data for bottom boundary condition
!              (:,:,1,:) => prescribed surface flux or surface value
!              (:,:,2,:) => transfer velocity                              [m/s]
!*settling*    LOGICAL Interprets argument wsink as (downward) settling
!                      velocities if .TRUE. 
!*ivarids*     INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: time_split, vadv, vimplicit
CHARACTER (LEN=lenname) :: f90_name
INTEGER :: ivar, k, klo, kup, nhpsi, npcc
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhdims, nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc1, array2dc2, &
                                            & psiobu2d, psiobv2d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: bcbot3d, bcsur3d, psic3d, psic_A, &
              & psic_B, uadvatu, ucorr, vadvatv, vcorr, vdifatw, wadvatw, wcorr
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*time_split* LOGICAL .TRUE. for directional splitting in time
!*vimplicit*  LOGICAL .TRUE. for implicit integration in the vertical
!*tridcfc*    REAL    Tridiagonal matrix for implicit vertical solution
!*psic_A*     REAL    psic during first integration step            [psic]
!*psic_B*     REAL    psic during second integration step           [psic]  
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'transport_at_C_4d1'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!
!1.1 Array shape and halo size
!-----------------------------
!

nhpsi = MAX(num_halo(iopt_hadv),MIN(1,iopt_hdif))
nhdims = nhalo
klo = MERGE(1,2,ibcbot.LE.2)
kup = MERGE(nz,nz-1,ibcsur.LE.2)

!
!1.2 Time integration
!--------------------
!

vimplicit = ((iopt_vadv.GT.0.AND.iopt_vadv_impl.GT.0).OR.&
           & (iopt_vdif_coef.GT.0.AND.iopt_vdif_impl.GT.0).OR.&
           & (ibcsur.EQ.4).OR.(ibcbot.EQ.4))
vadv = iopt_grid_nodim.NE.2.AND.iopt_vadv.GT.0
time_split = iopt_hadv.EQ.3.OR.iopt_vadv.EQ.3

!
!1.3 Allocate arrays
!-------------------
!

ALLOCATE (array2dc1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dc2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (psiobu2d(nobu,nz),STAT=errstat)
CALL error_alloc('psiobu2d',2,(/nobu,nz/),kndrtype)
ALLOCATE (psiobv2d(nobv,nz),STAT=errstat)
CALL error_alloc('psiobv2d',2,(/nobv,nz/),kndrtype)
ALLOCATE (psic3d(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('psic3d',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
ALLOCATE (uadvatu(2-nhalo:ncloc+nhalo,nrloc,nz),STAT=errstat)
CALL error_alloc('uadvatu',3,(/ncloc+2*nhalo-1,nrloc,nz/),kndrtype)
ALLOCATE (vadvatv(ncloc,2-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('vadvatv',3,(/ncloc,nrloc+2*nhalo-1,nz/),kndrtype)
ALLOCATE (vdifatw(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('vdifatw',3,(/ncloc,nrloc,nz+1/),kndrtype)
ALLOCATE (wadvatw(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('wadvatw',3,(/ncloc,nrloc,nz+1/),kndrtype)
ALLOCATE (bcbot3d(ncloc,nrloc,nbcb),STAT=errstat)
CALL error_alloc('bcbot3d',3,(/ncloc,nrloc,nbcb/),kndrtype)
ALLOCATE (bcsur3d(ncloc,nrloc,nbcs),STAT=errstat)
CALL error_alloc('bcsur3d',3,(/ncloc,nrloc,nbcs/),kndrtype)
ALLOCATE (tridcfc(ncloc,nrloc,nz,4),STAT=errstat)
CALL error_alloc('tridcfc',4,(/ncloc,nrloc,nz,4/),kndrtype)

IF (time_split) THEN
   ALLOCATE (psic_A(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
   CALL error_alloc('psic_A',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
   ALLOCATE (psic_B(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
   CALL error_alloc('psic_B',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
   ALLOCATE (ucorr(ncloc,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('ucorr',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (vcorr(ncloc,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('vcorr',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (wcorr(ncloc,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('wcorr',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
ENDIF

!
!1.4 Work space arrays
!---------------------
!

IF (.NOT.time_split) THEN
   WHERE (maskatc_int)
      array2dc1 = deptotatc_old(1:ncloc,1:nrloc)/deptotatc(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!2. Update each variable separately
!----------------------------------
!

ivar_200: DO ivar=1,novars

   IF (loglev1.GE.pglev) THEN
      CALL inquire_var(ivarids(ivar),f90_name=f90_name,numvar=ivar)
      WRITE (iolog,'(2X,A)') REPEAT(' ',pglev-1)//TRIM(f90_name)
   ENDIF

!
!2.1 Advective velocities
!------------------------
!
!2.1.1 Initialise
!----------------
!

   IF (iopt_hadv.GT.0) THEN
      uadvatu = ufvel; vadvatv = vfvel
   ELSE
      uadvatu = 0.0; vadvatv = 0.0
   ENDIF
   wadvatw = 0.0
   
   IF (vadv) THEN

!
!2.1.2 Without settling
!----------------------
!

      IF (iopt_grid_nodim.EQ.3) THEN
         wadvatw = wfvel(1:ncloc,1:nrloc,:)
      ENDIF

!
!2.1.3 With settling
!-------------------
!

      IF (settling) THEN

!        ---without sigma-grid correction
         IF (iopt_curr_wfall.EQ.1.OR.iopt_hadv.EQ.0) THEN
            wadvatw = wadvatw - wsink(1:ncloc,1:nrloc,:,ivar)
!        ---with sigma-grid correction
         ELSEIF (iopt_curr_wfall.EQ.2) THEN
            CALL transf_vertical_fall(wsink(:,:,:,ivar),uadvatu,vadvatv,&
                                    & wadvatw,1)
         ENDIF

      ENDIF

   ENDIF

!
!2.2 Store into local arrays
!---------------------------
!

   psic3d = psic(:,:,:,ivar)
   vdifatw = vdifcoefatw(:,:,:,ivar)
   IF (time_split) THEN
      psic_A = psic3d
      psic_B = psic3d
   ENDIF
   IF (nobu.GT.0) psiobu2d = psiobu(:,:,ivar)
   IF (nobv.GT.0) psiobv2d = psiobv(:,:,ivar)
   bcbot3d = bcbot(:,:,:,ivar); bcsur3d = bcsur(:,:,:,ivar)

!
!2.3 Dirichlet boundary conditions
!---------------------------------
!
!  ---surface
   IF (ibcsur.EQ.3) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,nz,1) = 0.0
         tridcfc(:,:,nz,2) = 1.0
         tridcfc(:,:,nz,3) = 0.0
         tridcfc(:,:,nz,4) = bcsur3d(:,:,1)
      END WHERE
   ELSEIF (ibcsur.EQ.4) THEN
      WHERE (maskatc_int)
         array2dc2 = delzatc(1:ncloc,1:nrloc,nz) + &
                   & 2.0*delzatw(1:ncloc,1:nrloc,nz)
         tridcfc(:,:,nz,1) = -delzatc(1:ncloc,1:nrloc,nz)/array2dc2
         tridcfc(:,:,nz,2) = 1.0
         tridcfc(:,:,nz,3) = 0.0
         tridcfc(:,:,nz,4) = 2.0*delzatw(1:ncloc,1:nrloc,nz)*bcsur3d(:,:,1)&
                           & /array2dc2
      END WHERE
   ENDIF

!  ---bottom
   IF (ibcbot.EQ.3) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,1,1) = 0.0
         tridcfc(:,:,1,2) = 1.0
         tridcfc(:,:,1,3) = 0.0
         tridcfc(:,:,1,4) = bcbot3d(:,:,1)
      END WHERE
   ELSEIF (ibcbot.EQ.4) THEN
      WHERE (maskatc_int)
         array2dc2 = 2.0*delzatw(1:ncloc,1:nrloc,2) + &
                   & delzatc(1:ncloc,1:nrloc,1)
         tridcfc(:,:,1,1) = 0.0
         tridcfc(:,:,1,2) = 1.0
         tridcfc(:,:,1,3) = -delzatc(1:ncloc,1:nrloc,1)/array2dc2
         tridcfc(:,:,1,4) = 2.0*delzatw(1:ncloc,1:nrloc,2)*bcbot3d(:,:,1)&
                          & /array2dc2
      END WHERE
   ENDIF

!
!2.4 Source terms
!----------------
!

   k_240:DO k=klo,kup
      WHERE (maskatc_int)
         source(:,:,k,ivar) = delt3d*source(:,:,k,ivar)
      END WHERE
   ENDDO k_240
   
!
!2.5 Corrector terms
!-------------------
!

   IF (time_split.AND.ivar.EQ.1) THEN
      CALL Xcorr_at_C(ucorr,klo,kup)
      CALL Ycorr_at_C(vcorr,klo,kup)
      CALL Zcorr_at_C(wcorr,klo,kup)
   ENDIF

!
!2.6 Time integration without operator splitting
!-----------------------------------------------
!

   IF (.NOT.time_split) THEN

!
!2.6.1 Initialise
!----------------
!
   
      IF (vimplicit) THEN
         k_261: DO k=klo,kup
            WHERE (maskatc_int)
               tridcfc(:,:,k,1) = 0.0
               tridcfc(:,:,k,3) = 0.0
            END WHERE
         ENDDO k_261
      ENDIF
      
!
!2.6.2 Time derivative
!---------------------
!

      k_262: DO k=klo,kup
         WHERE (maskatc_int) 
            tridcfc(:,:,k,2) = 1.0
            tridcfc(:,:,k,4) = array2dc1*psic3d(1:ncloc,1:nrloc,k)
         END WHERE
      ENDDO k_262

!
!2.6.3 Source term
!-----------------
!

      k_263: DO k=klo,kup
         WHERE (maskatc_int)
            tridcfc(:,:,k,4) = tridcfc(:,:,k,4) + source(:,:,k,ivar)
         END WHERE
      ENDDO k_263

!
!2.6.4 Advection
!-------------
!
!     ---horizontal
      IF (iopt_hadv.GT.0) THEN
         CALL Xadv_at_C(psic3d,tridcfc,1,iopt_hadv,nhpsi,klo,kup,psiobu2d,&
                      & uadvatu)
         CALL Yadv_at_C(psic3d,tridcfc,1,iopt_hadv,nhpsi,klo,kup,psiobv2d,&
                      & vadvatv)
      ENDIF
!     ---vertical
      IF (vadv) CALL Zadv_at_C(psic3d,tridcfc,wadvatw,1,iopt_vadv,klo,kup)

!
!2.6.5 Diffusion
!---------------
!
!     ---horizontal
      IF (iopt_hdif.GT.0) THEN
         CALL Xdif_at_C(psic3d,tridcfc,1,klo,kup)
         CALL Ydif_at_C(psic3d,tridcfc,1,klo,kup)
      ENDIF
!     ---vertical
      CALL Zdif_at_C(psic3d,tridcfc,vdifatw,1,ibcsur,ibcbot,nbcs,nbcb,&
                   & bcsur3d,bcbot3d,klo,kup)

!
!2.6.6 Update psic
!-----------------
!

      IF (vimplicit) THEN
         CALL tridiag_vert(tridcfc,psic3d,nhdims,nz,1,'C  ')
      ELSE
         k_266: DO k=1,nz
            WHERE (maskatc_int)
               psic(1:ncloc,1:nrloc,k,ivar) = tridcfc(:,:,k,4)
            END WHERE
         ENDDO k_266
      ENDIF
      
      CYCLE ivar_200

   ENDIF

!
!2.7 Time integration with operator splitting - step A
!-----------------------------------------------------
!
!2.7.1 X-derivative terms
!----------------------
!
!  ---time derivative, corrector term
   k_2711: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,2) = 1.0
         tridcfc(:,:,k,4) = psic3d(1:ncloc,1:nrloc,k)*(1.0+ucorr(:,:,k))
      END WHERE
   ENDDO k_2711

!  ---horizontal advection
   IF (iopt_hadv.GT.0) THEN
      CALL Xadv_at_C(psic3d,tridcfc,1,iopt_hadv,nhpsi,klo,kup,psiobu2d,uadvatu)
   ENDIF

!  ---horizontal diffusion
   IF (iopt_hdif.GT.0) THEN
      CALL Xdif_at_C(psic3d,tridcfc,1,klo,kup)
   ENDIF

!  ---update psic_A
   k_2712: DO k=klo,kup
      WHERE (maskatc_int)
         psic_A(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)/tridcfc(:,:,k,2)
      END WHERE
   ENDDO k_2712

!
!2.7.2 Y-derivative terms
!------------------------
!
!  ---time derivative, corrector term
   k_2721: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,2) = 1.0
         tridcfc(:,:,k,4) = psic_A(1:ncloc,1:nrloc,k)*(1.0+vcorr(:,:,k))
      END WHERE
   ENDDO k_2721

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/1-nhalo,1-nhalo,1/); nhexch = (/0,0,nhpsi,nhpsi/)
      CALL exchange_mod(psic_A,lbounds,nhexch,0)
   ENDIF

!  ---horizontal advection
   IF (iopt_hadv.GT.0) THEN
      CALL Yadv_at_C(psic_A,tridcfc,1,iopt_hadv,nhpsi,klo,kup,psiobv2d,vadvatv)
   ENDIF

!  ---horizontal diffusion
   IF (iopt_hdif.GT.0) THEN
      CALL Ydif_at_C(psic_A,tridcfc,1,klo,kup)
   ENDIF

!  ---update psic_A
   k_2722: DO k=klo,kup
      WHERE (maskatc_int)
         psic_A(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)/tridcfc(:,:,k,2)
      END WHERE
   ENDDO k_2722

!  ---exchange halo sections (if needed)
   IF (iopt_MPI.EQ.1.AND.iopt_hdif_scal.GT.1) THEN
      lbounds = (/1-nhalo,1-nhalo,1/); nhexch = nhpsi
      CALL exchange_mod(psic_A,lbounds,nhexch,0,corners=.FALSE.)
   ENDIF

!
!2.7.3 Z-derivative/source terms
!-------------------------------
!
!2.7.3.1 Time derivative, source and corrector terms
!---------------------------------------------------
!

   IF (vimplicit) THEN
      k_27311: DO k=klo,kup
         WHERE (maskatc_int)
            tridcfc(:,:,k,1) = 0.0
            tridcfc(:,:,k,3) = 0.0
         END WHERE
      ENDDO k_27311
   ENDIF
   k_27312: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4) = psic_A(1:ncloc,1:nrloc,k)*(1.0+&
                          & (1.0-theta_vadv)*wcorr(:,:,k)) + source(:,:,k,ivar)
         tridcfc(:,:,k,2) = 1.0 - theta_vadv*wcorr(:,:,k)     
      END WHERE
   ENDDO k_27312

!
!2.8.3.2 Vertical advection
!--------------------------
!

   IF (vadv) CALL Zadv_at_C(psic_A,tridcfc,wadvatw,1,iopt_vadv,klo,kup)

!
!2.8.3.3 Vertical diffusion
!--------------------------
!

   CALL Zdif_at_C(psic_A,tridcfc,vdifatw,1,ibcsur,ibcbot,nbcs,nbcb,bcsur3d,&
                & bcbot3d,klo,kup)

!
!2.8.3.4 Update psic_A
!---------------------
!

   IF (vimplicit) THEN
      CALL tridiag_vert(tridcfc,psic_A,nhdims,nz,1,'C  ')
   ELSE
      k_2834: DO k=1,nz
         WHERE (maskatc_int)
            psic_A(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)/tridcfc(:,:,k,2)
         END WHERE
      ENDDO k_2834
   ENDIF

!
!2.9 Time integration with operator splitting - step B
!-----------------------------------------------------
!
!2.9.1 Z-derivative/source terms
!-------------------------------
!
!2.9.1.1 Time derivative, source and corrector terms
!---------------------------------------------------
!

   IF (vimplicit) THEN
      k_29111: DO k=klo,kup
         WHERE (maskatc_int)
            tridcfc(:,:,k,1) = 0.0
            tridcfc(:,:,k,3) = 0.0
         END WHERE
      ENDDO k_29111
   ENDIF
   k_29112: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4) = psic3d(1:ncloc,1:nrloc,k)*&
                          & (1.0+(1.0-theta_vadv)*wcorr(:,:,k)) + &
                          & source(:,:,k,ivar)
         tridcfc(:,:,k,2) = 1.0 - theta_vadv*wcorr(:,:,k)
      END WHERE
   ENDDO k_29112

!
!2.9.1.2 Vertical advection
!--------------------------
!

   IF (vadv) CALL Zadv_at_C(psic3d,tridcfc,wadvatw,1,iopt_vadv,klo,kup)

!
!2.9.1.3 Vertical diffusion
!--------------------------
!

   CALL Zdif_at_C(psic3d,tridcfc,vdifatw,1,ibcsur,ibcbot,nbcs,nbcb,bcsur3d,&
                & bcbot3d,klo,kup)

!
!2.9.1.4 Update psic_B
!---------------------
!

   IF (vimplicit) THEN
      CALL tridiag_vert(tridcfc,psic_B,nhdims,nz,1,'C  ')
   ELSE
      k_2914: DO k=1,nz
         WHERE (maskatc_int)
            psic_B(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)
         END WHERE
      ENDDO k_2914
   ENDIF

!
!2.9.2 Y-derivative terms
!------------------------
!
!  ---time derivative, corrector term
   k_2921: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,2) = 1.0
         tridcfc(:,:,k,4) = psic_B(1:ncloc,1:nrloc,k)*(1.0+vcorr(:,:,k))
      END WHERE
   ENDDO k_2921

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      nhexch = (/0,0,nhpsi,nhpsi/)
      CALL exchange_mod(psic_B,lbounds,nhexch,0)
   ENDIF

!  ---horizontal advection
   IF (iopt_hadv.GT.0) THEN
      CALL Yadv_at_C(psic_B,tridcfc,1,iopt_hadv,nhpsi,klo,kup,psiobv2d,vadvatv)
   ENDIF

!  ---horizontal diffusion
   IF (iopt_hdif.GT.0) THEN
      CALL Ydif_at_C(psic_B,tridcfc,1,klo,kup)
   ENDIF

!  ---update psic_B
   k_2922: DO k=klo,kup
      WHERE (maskatc_int)
         psic_B(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)/tridcfc(:,:,k,2)
      END WHERE
   ENDDO k_2922

!
!2.9.3 X-derivative terms
!------------------------
!
!  ---time derivative, corrector term
   k_2931: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,2) = 1.0
         tridcfc(:,:,k,4) = psic_B(1:ncloc,1:nrloc,k)*(1.0+ucorr(:,:,k))
      END WHERE
   ENDDO k_2931

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      nhexch = (/nhpsi,nhpsi,0,0/)
      CALL exchange_mod(psic_B,lbounds,nhexch,0)
   ENDIF

!  ---horizontal advection
   IF (iopt_hadv.GT.0) THEN
      CALL Xadv_at_C(psic_B,tridcfc,1,iopt_hadv,nhpsi,klo,kup,psiobu2d,uadvatu)
   ENDIF

!  ---horizontal diffusion
   IF (iopt_hdif.GT.0) THEN
      CALL Xdif_at_C(psic_B,tridcfc,1,klo,kup)
   ENDIF

!  ---update psic_B
   k_2932: DO k=klo,kup
      WHERE (maskatc_int)
         psic_B(1:ncloc,1:nrloc,k) = tridcfc(:,:,k,4)/tridcfc(:,:,k,2)
      END WHERE
   ENDDO k_2932

!
!2.10 Time integration with operator splitting - update psic
!-----------------------------------------------------------
!

   k_2100: DO k=1,nz
      WHERE (maskatc_int)
         psic(1:ncloc,1:nrloc,k,ivar) = 0.5*(psic_A(1:ncloc,1:nrloc,k)+&
                                           & psic_B(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_2100

ENDDO ivar_200

!
!3. Apply relaxation scheme
!--- ----------------------
!

IF (iopt_obc_relax.EQ.1.AND.inodesrlx(1).GT.0) THEN
   CALL relaxation_at_C(psic,psiobu2d,psiobv2d,1,iprofrlx)
ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array2dc1,array2dc2,psiobu2d,psiobv2d,psic3d,uadvatu,vadvatv,&
          & vdifatw,wadvatw,bcbot3d,bcsur3d,tridcfc)
IF (time_split) DEALLOCATE (psic_A,psic_B,ucorr,vcorr,wcorr)

CALL log_timer_out(npcc,itm_trans)


RETURN

END SUBROUTINE transport_at_C_4d1

!========================================================================

SUBROUTINE transport_at_C_4d2(psic,source,wsink,vdifcoefatw,novars,iopt_hadv,&
                            & iopt_vadv,iopt_hdif,psiobu,psiobv,iprofrlx,&
                            & ibcsur,ibcbot,nbcs,nbcb,bcsur,bcbot,settling,&
                            & ivarid)
!************************************************************************
!
! *transport_at_C_4d2* Solve advection-diffusion equation for a series of 3-D
!                      quantities at C-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Transport_Equations.f90  V2.11.2
!
! Description - equations are solved together
!
! Reference -
!
! Calling program -
!
! External calls - relaxation_at_C, Xadv_at_C, Xcorr_at_C, Xdif_at_C,
!                  Yadv_at_C, Ycorr_at_C, Ydif_at_C, Zadv_at_C, Zcorr_at_C,
!                  Zdif_at_C
!
! Module calls - error_alloc, exchange_mod, num_halo, tridiag_vert,
!                transf_vertical_fall
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE physpars
USE relaxation
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE nla_library, ONLY: tridiag_vert
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: num_halo

IMPLICIT NONE

!
!*  Arguments
!
LOGICAL, INTENT(IN) :: settling
INTEGER, INTENT(IN) :: ibcbot, ibcsur, iopt_hadv, iopt_hdif, iopt_vadv, &
                     & ivarid, nbcb, nbcs, novars
INTEGER, INTENT(IN), DIMENSION(norlxzones) :: iprofrlx
REAL, INTENT(IN), DIMENSION(nobu,nz,novars) :: psiobu
REAL, INTENT(IN), DIMENSION(nobv,nz,novars) :: psiobv
REAL, INTENT(INOUT), DIMENSION(1-nhalo:ncloc+nhalo,&
                             & 1-nhalo:nrloc+nhalo,nz,novars) :: psic
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,novars) :: source
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1,novars) :: vdifcoefatw
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1,novars) :: wsink
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcs,novars) :: bcsur
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcb,novars) :: bcbot

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*psic*        REAL    C-node quantity to be updated                      [psic]
!*source*      REAL    Source terms in transport equations              [psic/s]
!*wsink*       REAL    Vertical (positive downward) sinking velocities at
!                      W-nodes if argument settling equals .TRUE.,
!                      not used otherwise                                  [m/s]
!*vdifcoefatw* REAL    Vertical diffusion coefficient at W-nodes         [m^2/s]
!*novars*      INTEGER Number of variables
!*iopt_adv*    INTEGER Switch to select (physical) advection scheme
!                      (disabled if zero)
!*iopt_vadv*   INTEGER Switch to select vertical advection scheme
!                      (disabled if zero, must be positive if settling is
!                      .TRUE.)
!*iopt_hdif*   INTEGER Switch to select horizontal diffusion scheme
!                      (disabled if zero)
!*psiobu*      REAL    Vertical profiles at U-open boundaries
!                      (zero gradient condition if flagged)               [psic]
!*psiobv*      REAL    Vertical profiles at V-open boundaries
!                      (zero gradient condition if flagged)               [psic]
!*iprofrlx*    INTEGER Disables/enables relaxation at open boundary zones
!*ibcsur*      INTEGER Type of surface boundary condition
!                  = 0 => Neumann (zero flux)
!                  = 1 => Neumann (prescibed flux)
!                  = 2 => Neumann (using transfer velocity)
!                  = 3 => Dirichlet at first C-node below the surface
!                  = 4 => Dirichlet at the surface
!*ibcbot*      INTEGER Type of bottom boundary condition
!                  = 0 => Neumann (zero flux)
!                  = 1 => Neumann (prescibed flux)
!                  = 2 => Neumann (using transfer velocity)
!                  = 3 => Dirichlet at first C-node above the bottom
!                  = 4 => Dirichlet at the bottom
!*nbcs*        INTEGER Last dimension of array bcsur
!*nbcb*        INTEGER Last dimension of array bcbot
!*bcsur*       REAL    Data for surface boundary condition
!              (:,:,1,:) => prescribed surface flux or surface value
!              (:,:,2,:) => transfer velocity                              [m/s]
!*bcbot*       REAL    Data for bottom boundary condition
!              (:,:,1,:) => prescribed surface flux or surface value
!              (:,:,2,:) => transfer velocity                              [m/s]
!*settling*    LOGICAL Interprets argument wsink as (downward) settling
!                      velocities if .TRUE.
!*ivarid*      INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: time_split, vadv, vimplicit
INTEGER :: ivar, k, klo, kup, nhpsi, npcc
INTEGER, DIMENSION(4) :: lbounds, nhdims, nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ucorr, vcorr, wcorr
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: psic_A, psic_B, uadvatu, &
                                             & vadvatv, wadvatw
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: tridcfc

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*time_split* LOGICAL .TRUE. for directional splitting in time
!*vimplicit*  LOGICAL .TRUE. for implicit integration in the vertical
!*tridcfc*    REAL    Tridiagonal matrix for implicit vertical solution
!*psic_A*     REAL    psic during first integration step            [psic]
!*psic_B*     REAL    psic during second integration step           [psic]  
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'transport_at_C_4d2'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Initialise
!-------------
!
!1.1 Array shape and halo size
!-----------------------------
!

nhpsi = MAX(num_halo(iopt_hadv),MIN(1,iopt_hdif))
nhdims = nhalo
klo = MERGE(1,2,ibcbot.LE.2)
kup = MERGE(nz,nz-1,ibcsur.LE.2)

!
!1.2 Time integration
!--------------------
!

vimplicit = ((iopt_vadv.GT.0.AND.iopt_vadv_impl.GT.0).OR.&
           & (iopt_vdif_coef.GT.0.AND.iopt_vdif_impl.GT.0).OR.&
           & (ibcsur.EQ.4).OR.(ibcbot.EQ.4))
vadv = iopt_grid_nodim.NE.2.AND.iopt_vadv.GT.0
time_split = iopt_hadv.EQ.3.OR.iopt_vadv.EQ.3

!
!1.3 Allocate arrays
!-------------------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (tridcfc(ncloc,nrloc,nz,4,novars),STAT=errstat)
CALL error_alloc('tridcfc',5,(/ncloc,nrloc,nz,4,novars/),kndrtype)
ALLOCATE (uadvatu(2-nhalo:ncloc+nhalo,nrloc,nz,novars),STAT=errstat)
CALL error_alloc('uadvatu',4,(/ncloc+2*nhalo-1,nrloc,nz,novars/),kndrtype)
ALLOCATE (vadvatv(ncloc,2-nhalo:nrloc+nhalo,nz,novars),STAT=errstat)
CALL error_alloc('vadvatv',4,(/ncloc,nrloc+2*nhalo-1,nz,novars/),kndrtype)
ALLOCATE (wadvatw(ncloc,nrloc,nz+1,novars),STAT=errstat)
CALL error_alloc('wadvatw',4,(/ncloc,nrloc,nz+1,novars/),kndrtype)

IF (time_split) THEN
   ALLOCATE (psic_A(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz,novars),&
           & STAT=errstat)
   CALL error_alloc('psic_A',4,(/ncloc+2*nhalo,nrloc+2*nhalo,nz,novars/),&
                  & kndrtype)
   ALLOCATE (psic_B(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz,novars),&
           & STAT=errstat)
   CALL error_alloc('psic_B',4,(/ncloc+2*nhalo,nrloc+2*nhalo,nz,novars/),&
                  & kndrtype)
   ALLOCATE (ucorr(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('ucorr',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (vcorr(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('vcorr',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (wcorr(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('wcorr',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

!
!1.4 Dirichlet boundary conditions
!---------------------------------
!
!1.4.1 Surface
!-------------
!

IF (ibcsur.EQ.4) THEN
   WHERE (maskatc_int)
      array2dc = delzatc(1:ncloc,1:nrloc,nz) + 2.0*delzatw(1:ncloc,1:nrloc,nz)
   END WHERE
ENDIF

ivar_141: DO ivar=1,novars
   IF (ibcsur.EQ.3) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,nz,1,ivar) = 0.0
         tridcfc(:,:,nz,2,ivar) = 1.0
         tridcfc(:,:,nz,3,ivar) = 0.0
         tridcfc(:,:,nz,4,ivar) = bcsur(:,:,1,ivar)
      END WHERE
   ELSEIF (ibcsur.EQ.4) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,nz,1,ivar) = -delzatc(1:ncloc,1:nrloc,nz)/array2dc
         tridcfc(:,:,nz,2,ivar) = 1.0
         tridcfc(:,:,nz,3,ivar) = 0.0
         tridcfc(:,:,nz,4,ivar) = 2.0*delzatw(1:ncloc,1:nrloc,nz)*&
                                & bcsur(:,:,1,ivar)/array2dc
      END WHERE
   ENDIF
ENDDO ivar_141

!
!1.4.2 Bottom
!------------
!

IF (ibcbot.EQ.4) THEN
   WHERE (maskatc_int) 
      array2dc = 2.0*delzatw(1:ncloc,1:nrloc,2) + delzatc(1:ncloc,1:nrloc,1)
   END WHERE
ENDIF

ivar_142: DO ivar=1,novars
   IF (ibcbot.EQ.3) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,1,1,ivar) = 0.0
         tridcfc(:,:,1,2,ivar) = 1.0
         tridcfc(:,:,1,3,ivar) = 0.0
         tridcfc(:,:,1,4,ivar) = bcbot(:,:,1,ivar)
      END WHERE
   ELSEIF (ibcbot.EQ.4) THEN
      WHERE (maskatc_int)
         tridcfc(:,:,1,1,ivar) = 0.0
         tridcfc(:,:,1,2,ivar) = 1.0
         tridcfc(:,:,1,3,ivar) = -delzatc(1:ncloc,1:nrloc,1)/array2dc
         tridcfc(:,:,1,4,ivar) = 2.0*delzatw(1:ncloc,1:nrloc,2)*&
                               & bcbot(:,:,1,ivar)/array2dc
      END WHERE
   ENDIF
ENDDO ivar_142

!
!1.5 Work space arrays
!---------------------
!

IF (.NOT.time_split) THEN
   WHERE (maskatc_int)
      array2dc = deptotatc_old(1:ncloc,1:nrloc)/deptotatc(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!1.6 Local arrays
!----------------
!

IF (time_split) THEN
   psic_A = psic; psic_B = psic
ENDIF

!
!1.7 Source terms
!----------------
!

ivar_170: DO ivar=1,novars
k_170:DO k=klo,kup
   WHERE (maskatc_int)
      source(:,:,k,ivar) = delt3d*source(:,:,k,ivar)
   END WHERE
ENDDO k_170
ENDDO ivar_170

!
!1.8 Advective velocities
!------------------------
!
!1.8.1 Initialise
!----------------
!

IF (iopt_hadv.GT.0) THEN
   ivar_181: DO ivar=1,novars
      uadvatu(:,:,:,ivar) = ufvel
      vadvatv(:,:,:,ivar) = vfvel
   ENDDO ivar_181
ELSE
   uadvatu = 0.0; vadvatv = 0.0
ENDIF
wadvatw = 0.0

IF (vadv) THEN

!
!1.8.2 Without settling
!----------------------
!

   ivar_182: DO ivar=1,novars
      IF (iopt_grid_nodim.EQ.3) THEN
         wadvatw(:,:,:,ivar) = wfvel(1:ncloc,1:nrloc,:)
      ELSE
         wadvatw(:,:,:,ivar) = 0.0
      ENDIF
   ENDDO ivar_182

!
!1.8.3 With settling
!-------------------
!

   IF (settling) THEN

!     ---without sigma-grid correction
      IF (iopt_curr_wfall.EQ.1.OR.iopt_hadv.EQ.0) THEN
         ivar_183: DO ivar=1,novars
            wadvatw(:,:,:,ivar) = wadvatw(:,:,:,ivar) - &
                                & wsink(1:ncloc,1:nrloc,:,ivar)
         ENDDO ivar_183
!     ---with sigma-grid corrtection
      ELSEIF (iopt_curr_wfall.EQ.2) THEN
         CALL transf_vertical_fall(wsink,uadvatu,vadvatv,wadvatw,novars)
      ENDIF
      
   ENDIF

ENDIF

!
!1.9 Corrector terms
!--------------------
!

IF (time_split) THEN
   CALL Xcorr_at_C(ucorr,klo,kup)
   CALL Ycorr_at_C(vcorr,klo,kup)
   CALL Zcorr_at_C(wcorr,klo,kup)
ENDIF

!
!2. Time integration without operator splitting
!----------------------------------------------
!

IF (.NOT.time_split) THEN

!
!2.1 Initialise
!--------------
!
   
   IF (vimplicit) THEN
      ivar_210: DO ivar=1,novars
      k_210: DO k=klo,kup
         WHERE (maskatc_int)
            tridcfc(:,:,k,1,ivar) = 0.0
            tridcfc(:,:,k,3,ivar) = 0.0
         END WHERE
      ENDDO k_210
      ENDDO ivar_210
   ENDIF

!
!2.2 Time derivative
!-------------------
!

   ivar_220: DO ivar=1,novars
   k_220: DO k=klo,kup
      WHERE (maskatc_int) 
         tridcfc(:,:,k,2,ivar) = 1.0
         tridcfc(:,:,k,4,ivar) = array2dc*psic(1:ncloc,1:nrloc,k,ivar)
      END WHERE
   ENDDO k_220
   ENDDO ivar_220

!
!2.3 Source term
!---------------
!

   ivar_230: DO ivar=1,novars
   k_230: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,4,ivar) = tridcfc(:,:,k,4,ivar) + source(:,:,k,ivar)
      END WHERE
   ENDDO k_230
   ENDDO ivar_230

!
!2.4 Advection
!-------------
!
!  ---horizontal
   IF (iopt_hadv.GT.0) THEN
      CALL Xadv_at_C(psic,tridcfc,novars,iopt_hadv,nhpsi,klo,kup,psiobu,uadvatu)
      CALL Yadv_at_C(psic,tridcfc,novars,iopt_hadv,nhpsi,klo,kup,psiobv,vadvatv)
   ENDIF
!  ---vertical
   IF (vadv) CALL Zadv_at_C(psic,tridcfc,wadvatw,novars,iopt_vadv,klo,kup)

!
!2.5 Diffusion
!-------------
!
!  ---horizontal
   IF (iopt_hdif.GT.0) THEN
      CALL Xdif_at_C(psic,tridcfc,novars,klo,kup)
      CALL Ydif_at_C(psic,tridcfc,novars,klo,kup)
   ENDIF
!  ---vertical
   CALL Zdif_at_C(psic,tridcfc,vdifcoefatw,novars,ibcsur,ibcbot,nbcs,nbcb,&
                & bcsur,bcbot,klo,kup)

!
!2.6 Update psic
!---------------
!
   IF (vimplicit) THEN
      CALL tridiag_vert(tridcfc,psic,nhdims,nz,novars,'C  ')
   ELSE
      ivar_360: DO ivar=1,novars
      k_360: DO k=1,nz
         WHERE (maskatc_int)
            psic(1:ncloc,1:nrloc,k,ivar) = tridcfc(:,:,k,4,ivar)/&
                                         & tridcfc(:,:,k,2,ivar)
         END WHERE
      ENDDO k_360
      ENDDO ivar_360
   ENDIF

   GOTO 1000

ENDIF

!
!3. Time integration with operator splitting - step A
!----------------------------------------------------
!
!3.1 X-derivative terms
!----------------------
!
!---time derivative, corrector term
ivar_311: DO ivar=1,novars
k_311: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,2,ivar) = 1.0
      tridcfc(:,:,k,4,ivar) = psic(1:ncloc,1:nrloc,k,ivar)*(1.0+ucorr(:,:,k))
   END WHERE
ENDDO k_311
ENDDO ivar_311

!---horizontal advection
IF (iopt_hadv.GT.0) THEN
   CALL Xadv_at_C(psic,tridcfc,novars,iopt_hadv,nhpsi,klo,kup,psiobu,uadvatu)
ENDIF

!---horizontal diffusion
IF (iopt_hdif.GT.0) THEN
   CALL Xdif_at_C(psic,tridcfc,novars,klo,kup)
ENDIF

!---update psic_A
ivar_312: DO ivar=1,novars
k_312: DO k=klo,kup
   WHERE (maskatc_int)
      psic_A(1:ncloc,1:nrloc,k,ivar) = tridcfc(:,:,k,4,ivar)/&
                                     & tridcfc(:,:,k,2,ivar)
   END WHERE
ENDDO k_312
END DO ivar_312

!
!3.2 Y-derivative terms
!----------------------
!
!---time derivative, corrector term
ivar_321: DO ivar=1,novars
k_321: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,2,ivar) = 1.0
      tridcfc(:,:,k,4,ivar) = psic_A(1:ncloc,1:nrloc,k,ivar)*(1.0+vcorr(:,:,k))
   END WHERE
ENDDO k_321
ENDDO ivar_321

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   lbounds = (/1-nhalo,1-nhalo,1,1/); nhexch = (/0,0,nhpsi,nhpsi/)
   CALL exchange_mod(psic_A,lbounds,nhexch,0)
ENDIF

!---horizontal advection
IF (iopt_hadv.GT.0) THEN
   CALL Yadv_at_C(psic_A,tridcfc,novars,iopt_hadv,nhpsi,klo,kup,psiobv,vadvatv)
ENDIF

!---horizontal diffusion
IF (iopt_hdif.GT.0) THEN
   CALL Ydif_at_C(psic_A,tridcfc,novars,klo,kup)
ENDIF

!---update psic_A
ivar_322: DO ivar=1,novars
k_322: DO k=klo,kup
   WHERE (maskatc_int)
      psic_A(1:ncloc,1:nrloc,k,ivar) = tridcfc(:,:,k,4,ivar)/&
                                     & tridcfc(:,:,k,2,ivar)
   END WHERE
ENDDO k_322
END DO ivar_322

!---exchange halo sections (if needed)
IF (iopt_MPI.EQ.1.AND.iopt_hdif_scal.GT.1) THEN
   lbounds = (/1-nhalo,1-nhalo,1,1/); nhexch = nhpsi
   CALL exchange_mod(psic_A,lbounds,nhexch,0,corners=.FALSE.)
ENDIF

!
!3.3 Z-derivative/source terms
!-----------------------------
!
!3.3.1 Time derivative, source and corrector terms
!-------------------------------------------------
!

IF (vimplicit) THEN
   ivar_3311: DO ivar=1,novars
   k_3311: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,1,ivar) = 0.0
         tridcfc(:,:,k,3,ivar) = 0.0
      END WHERE
   ENDDO k_3311
   ENDDO ivar_3311
ENDIF
ivar_3312: DO ivar=1,novars
k_3312: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,4,ivar) = psic_A(1:ncloc,1:nrloc,k,ivar)*&
                            & (1.0+(1.0-theta_vadv)*wcorr(:,:,k)) + &
                            & source(:,:,k,ivar)
      tridcfc(:,:,k,2,ivar) = 1.0 - theta_vadv*wcorr(:,:,k)
   END WHERE
ENDDO k_3312
ENDDO ivar_3312

!
!3.3.2 Vertical advection
!------------------------
!

IF (vadv) CALL Zadv_at_C(psic_A,tridcfc,wadvatw,novars,iopt_vadv,klo,kup)

!
!3.3.3 Vertical diffusion
!------------------------
!

CALL Zdif_at_C(psic_A,tridcfc,vdifcoefatw,novars,ibcsur,ibcbot,nbcs,nbcb,&
             & bcsur,bcbot,klo,kup)

!
!3.3.4 Update psic_A
!-------------------
!

IF (vimplicit) THEN
   CALL tridiag_vert(tridcfc,psic_A,nhdims,nz,novars,'C  ')
ELSE
   ivar_334: DO ivar=1,novars
   k_334: DO k=1,nz
      WHERE (maskatc_int)
         psic_A(1:ncloc,1:nrloc,k,ivar) = tridcfc(:,:,k,4,ivar)/&
                                        & tridcfc(:,:,k,2,ivar)
      END WHERE
   ENDDO k_334
   ENDDO ivar_334
ENDIF

!
!4. Time integration with operator splitting - step B
!----------------------------------------------------
!
!4.1 Z-derivative/source terms
!-----------------------------
!
!4.1.1 Time derivative, source and corrector terms
!-------------------------------------------------
!

IF (vimplicit) THEN
   ivar_4111: DO ivar=1,novars
   k_4111: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfc(:,:,k,1,ivar) = 0.0
         tridcfc(:,:,k,3,ivar) = 0.0
      END WHERE
   ENDDO k_4111
   ENDDO ivar_4111
ENDIF
ivar_4112: DO ivar=1,novars
k_4112: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,4,ivar) = psic(1:ncloc,1:nrloc,k,ivar)*&
                       & (1.0+(1.0-theta_vadv)*wcorr(:,:,k)) + &
                       & source(:,:,k,ivar)
      tridcfc(:,:,k,2,ivar) = 1.0 - theta_vadv*wcorr(:,:,k)
   END WHERE
ENDDO k_4112
ENDDO ivar_4112

!
!4.1.2 Vertical advection
!------------------------
!

IF (vadv) CALL Zadv_at_C(psic,tridcfc,wadvatw,novars,iopt_vadv,klo,kup)

!
!4.1.3 Vertical diffusion
!------------------------
!

CALL Zdif_at_C(psic,tridcfc,vdifcoefatw,novars,ibcsur,ibcbot,nbcs,nbcb,&
             & bcsur,bcbot,klo,kup)

!
!4.1.4 Update psic_B
!-------------------
!

IF (vimplicit) THEN
   CALL tridiag_vert(tridcfc,psic_B,nhdims,nz,novars,'C  ')
ELSE
   ivar_414: DO ivar=1,novars
   k_414: DO k=1,nz
      WHERE (maskatc_int)
         psic_B(1:ncloc,1:nrloc,k,ivar) = tridcfc(:,:,k,4,ivar)/&
                                        & tridcfc(:,:,k,2,ivar)
      END WHERE
   ENDDO k_414
   ENDDO ivar_414
ENDIF

!
!4.2 Y-derivative terms
!----------------------
!
!---time derivative, corrector term
ivar_421: DO ivar=1,novars
k_421: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,2,ivar) = 1.0
      tridcfc(:,:,k,4,ivar) = psic_B(1:ncloc,1:nrloc,k,ivar)*(1.0+vcorr(:,:,k))
   END WHERE
ENDDO k_421
ENDDO ivar_421

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/0,0,nhpsi,nhpsi/)
   CALL exchange_mod(psic_B,lbounds,nhexch,0)
ENDIF

!---horizontal advection
IF (iopt_hadv.GT.0) THEN
   CALL Yadv_at_C(psic_B,tridcfc,novars,iopt_hadv,nhpsi,klo,kup,psiobv,vadvatv)
ENDIF

!---horizontal diffusion
IF (iopt_hdif.GT.0) THEN
   CALL Ydif_at_C(psic_B,tridcfc,novars,klo,kup)
ENDIF

!---update psic_B
ivar_422: DO ivar=1,novars
k_422: DO k=klo,kup
   WHERE (maskatc_int)
      psic_B(1:ncloc,1:nrloc,k,ivar) = tridcfc(:,:,k,4,ivar)/&
                                     & tridcfc(:,:,k,2,ivar)
   END WHERE
ENDDO k_422
ENDDO ivar_422

!
!4.3 X-derivative terms
!----------------------
!
!---time derivative, corrector term
ivar_431: DO ivar=1,novars
k_431: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfc(:,:,k,2,ivar) = 1.0
      tridcfc(:,:,k,4,ivar) = psic_B(1:ncloc,1:nrloc,k,ivar)*(1.0+ucorr(:,:,k))
   END WHERE
ENDDO k_431
ENDDO ivar_431

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/nhpsi,nhpsi,0,0/)
   CALL exchange_mod(psic_B,lbounds,nhexch,0)
ENDIF

!---horizontal advection
IF (iopt_hadv.GT.0) THEN
   CALL Xadv_at_C(psic_B,tridcfc,novars,iopt_hadv,nhpsi,klo,kup,psiobu,uadvatu)
ENDIF

!---horizontal diffusion
IF (iopt_hdif.GT.0) THEN
   CALL Xdif_at_C(psic_B,tridcfc,novars,klo,kup)
ENDIF

!---update psic_B
ivar_432: DO ivar=1,novars
k_432: DO k=klo,kup
   WHERE (maskatc_int)
      psic_B(1:ncloc,1:nrloc,k,ivar) = tridcfc(:,:,k,4,ivar)/&
                                     & tridcfc(:,:,k,2,ivar)
   END WHERE
ENDDO k_432
ENDDO ivar_432

!
!5. Time integration with operator splitting - update psic
!---------------------------------------------------------
!

ivar_510: DO ivar=1,novars
k_510: DO k=1,nz
   WHERE (maskatc_int)
      psic(1:ncloc,1:nrloc,k,ivar) = 0.5*(psic_A(1:ncloc,1:nrloc,k,ivar)+&
                                        & psic_B(1:ncloc,1:nrloc,k,ivar))
   END WHERE
ENDDO k_510
ENDDO ivar_510

1000 CONTINUE

!
!6. Apply relaxation scheme
!--------------------------
!

IF (iopt_obc_relax.EQ.1.AND.inodesrlx(1).GT.0) THEN
   CALL relaxation_at_C(psic,psiobu,psiobv,novars,iprofrlx)
ENDIF

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE (array2dc,tridcfc,uadvatu,vadvatv,wadvatw)
IF (time_split) DEALLOCATE (psic_A,psic_B,ucorr,vcorr,wcorr)

CALL log_timer_out(npcc,itm_trans)


RETURN

END SUBROUTINE transport_at_C_4d2

!========================================================================

SUBROUTINE transport_at_U_3d(source,sink,ibcsur,ibcbot,nbcs,nbcb,bcsur,bcbot)
!************************************************************************
!
! *transport_at_U_3d* Solve advection-diffusion equation for the 3-D U-current
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Transport_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - current_pred
!
! External calls - Xadv_at_U_3d, Xdif_at_U_3d, Yadv_at_U_3d, Ydif_at_U_3d,
!                  Zadv_at_U, Zdif_at_U
!
! Module calls - error_alloc, exchange_mod, num_halo, tridiag_vert, Warr_at_UW
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
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Warr_at_UW
USE error_routines, ONLY: error_alloc
USE nla_library, ONLY: tridiag_vert
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: num_halo

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ibcsur, ibcbot, nbcs, nbcb
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz) :: sink, source
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcs) :: bcsur
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcb) :: bcbot

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*source*     REAL    Source terms in transport equations             [m/s^2]
!*sink*       REAL    Sink terms in transport equations                 [1/s]
!*ibcsur*     INTEGER Type of surface boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!*ibcbot*     INTEGER Type of bottom boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!*nbcs*       INTEGER Last dimension of array bcsur
!*nbcb*       INTEGER Last dimension of array bcbot
!*bcsur*      REAL    Data for surface boundary condition
!             (:,:,1) => prescribed surface flux
!             (:,:,2) => transfer velocity                              [m/s]
!*bcbot*      REAL    Data for bottom boundary condition
!             (:,:,1) => prescribed bottom flux
!             (:,:,2) => transfer velocity                              [m/s]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fimp, fsink, time_split, vflag, vimplicit
INTEGER :: k, nhpsi, npcc
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhdims, nhexch
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advdif3d, sinkatu, sourceatu, &
             & uvel_A, uvel_A_prev, uvel_B, vdifcoefatuw, wadvatuw, wstokesatuw
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: tridcfu, tridcfu_prev, &
                                             & tridcfu_A_prev, tridcfu_B_prev

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*time_split*  LOGICAL .TRUE. for directional splitting in time
!*vimplicit*   LOGICAL .TRUE. for implicit integration in the vertical
!*maskatu*     LOGICAL Mask array at U-nodes
!*tridcfu*     REAL    Tridiagonal matrix for implicit vertical solution
!*uvel_A*      REAL    U-current during first integration step        [m/s]
!*uvel_B*      REAL    U-current during second integration step       [m/s]  
!*vdifcoefatuw*REAL    Vertical diffusion coefficient at UW-nodes   [m^2/s]
!*wadvatuw*    REAL    Advective velocity in Z-direction at UW-nodes  [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'transport_at_U_3d'
CALL log_timer_in(npcc,ivarid=iarr_uvel)

!
!1. Initialise
!-------------
!
!1.1 Array shape and halo size
!-----------------------------
!

nhpsi = MAX(num_halo(iopt_adv_3D),iopt_hdif_3D)
nhdims = nhalo

!
!1.2 Time integration
!--------------------
!

fsink = iopt_weibar.EQ.1
vimplicit = ((iopt_adv_3D.GT.0.AND.iopt_vadv_impl.GT.0).OR.&
           & (iopt_vdif_coef.GT.0.AND.iopt_vdif_impl.GT.0).OR.&
           & (ibcsur.EQ.4).OR.(ibcbot.EQ.4))
time_split = iopt_adv_3D.EQ.3
vflag = iopt_grid_nodim.EQ.3.AND.iopt_hydro_impl.EQ.0 
fimp = iopt_hydro_impl.EQ.1.AND.maxitsimp.NE.1

!
!1.3 Allocate arrays
!-------------------
!

IF (fimp.AND.itsimp.EQ.1) THEN
   IF (.NOT.time_split) THEN
      IF (ALLOCATED(tridcfu_prev)) DEALLOCATE (tridcfu_prev)
      ALLOCATE (tridcfu_prev(ncloc,nrloc,nz,4),STAT=errstat)
      CALL error_alloc('tridcfu_prev',4,(/ncloc,nrloc,nz,4/),kndrtype)
   ELSE
      IF (ALLOCATED(tridcfu_A_prev)) THEN
         DEALLOCATE (tridcfu_A_prev,tridcfu_B_prev,uvel_A_prev)
      ENDIF
      ALLOCATE (tridcfu_A_prev(ncloc,nrloc,nz,4),STAT=errstat)
      CALL error_alloc('tridcfu_A_prev',4,(/ncloc,nrloc,nz,4/),kndrtype)
      ALLOCATE (tridcfu_B_prev(ncloc,nrloc,nz,4),STAT=errstat)
      CALL error_alloc('tridcfu_B_prev',4,(/ncloc,nrloc,nz,4/),kndrtype)
      ALLOCATE (uvel_A_prev(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),&
              & STAT=errstat)
      CALL error_alloc('uvel_A_prev',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),&
                     & kndrtype)
   ENDIF
ENDIF

ALLOCATE (advdif3d(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('advdif3d',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (maskatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (sourceatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('sourceatu',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (tridcfu(ncloc,nrloc,nz,4),STAT=errstat)
CALL error_alloc('tridcfu',4,(/ncloc,nrloc,nz,4/),kndrtype)
IF (itsimp.EQ.1) THEN
   ALLOCATE (sinkatu(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('sinkatu',3,(/ncloc,nrloc,nz/),kndrtype)
   IF (iopt_vdif_coef.GT.0) THEN
      ALLOCATE (vdifcoefatuw(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('vdifcoefatuw',3,(/ncloc,nrloc,nz+1/),kndrtype)
   ENDIF
   IF (iopt_adv_3D.GT.0) THEN
      ALLOCATE (wadvatuw(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('wadvatuw',3,(/ncloc,nrloc,nz+1/),kndrtype)
      IF (iopt_waves_curr.EQ.1) THEN
         ALLOCATE (wstokesatuw(ncloc,nrloc,nz+1),STAT=errstat)
         CALL error_alloc('wstokesatuw',3,(/ncloc,nrloc,nz+1/),kndrtype)
      ENDIF
   ENDIF
ENDIF
IF (time_split) THEN
   ALLOCATE (uvel_A(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
   CALL error_alloc('uvel_A',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
   ALLOCATE (uvel_B(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
   CALL error_alloc('uvel_B',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
ENDIF

!
!1.4 Masks
!---------
!

maskatu = nodeatu(1:ncloc,1:nrloc,:).EQ.2

!
!1.5 Local arrays
!----------------
!
!---depth-integrated term
udevint = 0.0

!---temporary values of uvel
IF (time_split) THEN
   uvel_A = uvel; uvel_B = uvel
ENDIF

!
!1.6 Source/sink terms
!---------------------
!

WHERE (maskatu)
   sourceatu = delt3d*source
END WHERE
IF (itsimp.EQ.1.AND.fsink) THEN
   WHERE (maskatu)
      sinkatu = delt3d*sink
   END WHERE
ENDIF

!
!1.7 Advective currents
!----------------------
!

IF (itsimp.EQ.1.AND.iopt_adv_3D.GT.0) THEN
   CALL Warr_at_UW(wfvel(:,1:nrloc,:),wadvatuw,3,(/0,1,1/),&
                & (/ncloc,nrloc,nz+1/),1,iarr_wvel,.TRUE.)
   IF (iopt_waves_curr.EQ.1) THEN
      CALL Warr_at_UW(wstokesatw(:,1:nrloc,:),wstokesatuw,3,(/0,1,1/),&
                   & (/ncloc,nrloc,nz+1/),1,iarr_wstokesatw,.TRUE.)
      wadvatuw = wadvatuw + wstokesatuw 
   ENDIF
ENDIF

!
!1.8 Vertical diffusion coefficient
!----------------------------------
!

IF (itsimp.EQ.1.AND.iopt_vdif_coef.GT.0) THEN
   CALL Warr_at_UW(vdifcoefmom(:,1:nrloc,:),vdifcoefatuw,3,(/0,1,1/),&
                & (/ncloc,nrloc,nz+1/),1,iarr_vdifcoefmom,.TRUE.)
ENDIF

!
!2. Time integration without operator splitting
!----------------------------------------------
!

IF (.NOT.time_split) THEN

!
!2.1 Explicit scheme or first iteration
!--------------------------------------
!

   IF (itsimp.EQ.1) THEN

!
!2.1.1 Initialise
!----------------
!

      tridcfu(:,:,:,1) = 0.0
      tridcfu(:,:,:,2) = 1.0
      tridcfu(:,:,:,3) = 0.0
      tridcfu(:,:,:,4) = 0.0

!
!2.1.2 Time derivative (rhs part)
!--------------------------------
!

      WHERE (maskatu)
         tridcfu(:,:,:,4) = uvel_old(1:ncloc,1:nrloc,:)
      END WHERE

!
!2.1.3 Advection
!---------------
!

      IF (iopt_adv_3D.GT.0) THEN
         CALL Xadv_at_U_3d(uvel_old,advdif3d,nhpsi,vflag)
         WHERE (maskatu)
            tridcfu(:,:,:,4) = tridcfu(:,:,:,4) - &
                             & delt3d*advdif3d/delzatu(1:ncloc,:,:)
         END WHERE
         CALL Yadv_at_U_3d(uvel_old,advdif3d,nhpsi,vflag)
         WHERE (maskatu)
            tridcfu(:,:,:,4) = tridcfu(:,:,:,4) - &
                             & delt3d*advdif3d/delzatu(1:ncloc,:,:)
         END WHERE
         CALL Zadv_at_U(uvel_old,tridcfu,wadvatuw)
      ENDIF

!
!2.1.4 Diffusion
!---------------
!
!     ---horizontal 
      IF (iopt_hdif_3D.GT.0) THEN
         CALL Xdif_at_U_3d(uvel_old,advdif3d,vflag)
         WHERE (maskatu)
            tridcfu(:,:,:,4) = tridcfu(:,:,:,4) + &
                             & delt3d*advdif3d/delzatu(1:ncloc,:,:)
         END WHERE
         CALL Ydif_at_U_3d(uvel_old,advdif3d,vflag)
         WHERE (maskatu)
            tridcfu(:,:,:,4) = tridcfu(:,:,:,4) + &
                             & delt3d*advdif3d/delzatu(1:ncloc,:,:)
         END WHERE
      ENDIF

!     ---vertical
      IF (iopt_vdif_coef.GT.0) THEN
         CALL Zdif_at_U(uvel_old,tridcfu,vdifcoefatuw,ibcsur,ibcbot,nbcs,nbcb,&
                      & bcsur,bcbot)
      ENDIF
            
!
!2.1.5 Sink term
!---------------
!

      IF (fsink) THEN
         WHERE (maskatu)
            tridcfu(:,:,:,2) = tridcfu(:,:,:,2) + sinkatu
         END WHERE
      ENDIF
           
!
!2.1.6 Store at first iteration
!------------------------------
!

      IF (fimp) tridcfu_prev = tridcfu

   ENDIF

!
!2.2 Explicit or next iteration
!------------------------------
!
!2.2.1 Retrieve from first iteration
!-----------------------------------
!

   IF (fimp) tridcfu = tridcfu_prev
 
!
!2.2.2 Source term
!-----------------
!

   WHERE (maskatu)
      tridcfu(:,:,:,4) = tridcfu(:,:,:,4) + sourceatu
   END WHERE

!
!2.2.3 Update uvel
!-----------------
!

   IF (vimplicit) THEN
      CALL tridiag_vert(tridcfu,uvel,nhdims,nz,1,'U  ')
   ELSE
      WHERE (maskatu)
         uvel(1:ncloc,1:nrloc,:) = tridcfu(:,:,:,4)/tridcfu(:,:,:,2)
      END WHERE
      WHERE (nodeatu(1:ncloc,1:nrloc,:).EQ.1)
         uvel(1:ncloc,1:nrloc,:) = 0.0
      END WHERE
   ENDIF

   GOTO 1000

ENDIF

!
!3. Time integration with operator splitting - step A
!----------------------------------------------------
!
!3.1 Explicit scheme or first iteration
!--------------------------------------
!

IF (itsimp.EQ.1) THEN

!
!3.1.1 Initialise
!----------------
!

   tridcfu(:,:,:,1) = 0.0
   tridcfu(:,:,:,2) = 1.0
   tridcfu(:,:,:,3) = 0.0
   tridcfu(:,:,:,4) = 0.0

!   
!3.1.2 X-derivative terms
!------------------------
!
!  ---advection
   IF (iopt_adv_3D.GT.0) THEN
      CALL Xadv_at_U_3d(uvel_old,advdif3d,nhpsi,vflag)
      WHERE (maskatu)
         tridcfu(:,:,:,4) = tridcfu(:,:,:,4) - advdif3d/delzatu(1:ncloc,:,:)
      END WHERE
   ENDIF

!  ---diffusion
   IF (iopt_hdif_3D.GT.0) THEN
      CALL Xdif_at_U_3d(uvel_old,advdif3d,vflag)
      WHERE (maskatu)
         tridcfu(:,:,:,4) = tridcfu(:,:,:,4) + advdif3d/delzatu(1:ncloc,:,:)
      END WHERE
   ENDIF

!  ---update uvel_A
   WHERE (maskatu)
      uvel_A(1:ncloc,1:nrloc,:) = uvel_old(1:ncloc,1:nrloc,:) + &
                                & delt3d*tridcfu(:,:,:,4)
   END WHERE

!
!
!3.1.3 Y-derivative terms
!------------------------
!
!  ---initialise
   WHERE (maskatu) tridcfu(:,:,:,4) = 0.0

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/1-nhalo,1-nhalo,1/); nhexch = (/0,0,nhpsi,nhpsi/)
      CALL exchange_mod(uvel_A,lbounds,nhexch,iarr_uvel)
   ENDIF

!  ---horizontal advection
   IF (iopt_adv_3D.GT.0) THEN
      IF (vflag) CALL Yadv_at_U_3d(uvel_old,advdif3d,nhpsi,.TRUE.)
      CALL Yadv_at_U_3d(uvel_A,advdif3d,nhpsi,.FALSE.)
      WHERE (maskatu)
         tridcfu(:,:,:,4) = tridcfu(:,:,:,4) - advdif3d/delzatu(1:ncloc,:,:)
      END WHERE
   ENDIF

!  ---horizontal diffusion
   IF (iopt_hdif_3D.GT.0) THEN
      IF (vflag) CALL Ydif_at_U_3d(uvel_old,advdif3d,.TRUE.)
      CALL Ydif_at_U_3d(uvel_A,advdif3d,.FALSE.)
      WHERE (maskatu)
         tridcfu(:,:,:,4) = tridcfu(:,:,:,4) + advdif3d/delzatu(1:ncloc,:,:)
      END WHERE
   ENDIF
   
!  ---update uvel_A
   WHERE (maskatu)
      uvel_A(1:ncloc,1:nrloc,:) = uvel_A(1:ncloc,1:nrloc,:) + &
                                & delt3d*tridcfu(:,:,:,4)
   END WHERE

!
!3.1.4 Z-derivative terms
!------------------------
!
!  ---initialise
   WHERE (maskatu) tridcfu(:,:,:,4) = 0.0

!  ---vertical advection
   IF (iopt_adv_3D.GT.0) THEN
      CALL Zadv_at_U(uvel_A,tridcfu,wadvatuw)
   ENDIF

!  ---vertical diffusion
   IF (iopt_vdif_coef.GT.0) THEN
      CALL Zdif_at_U(uvel_A,tridcfu,vdifcoefatuw,ibcsur,ibcbot,nbcs,nbcb,&
                   & bcsur,bcbot)
   ENDIF

!  ---sink term
   IF (fsink) THEN
      WHERE (maskatu)
         tridcfu(:,:,:,2) = tridcfu(:,:,:,2) + sinkatu
      END WHERE
   ENDIF
           
!
!3.1.5 Store at first iteration
!------------------------------
!

   IF (fimp) THEN
      tridcfu_A_prev = tridcfu
      uvel_A_prev = uvel_A
   ENDIF

ENDIF

!
!3.2 Explicit or next iteration
!------------------------------
!
!3.2.1 Retrieve from first iteration
!-----------------------------------
!

IF (fimp) THEN
   tridcfu = tridcfu_A_prev
   uvel_A = uvel_A_prev
ENDIF
 
!
!3.2.2 Source term and time derivative (rhs part)
!------------------------------------------------
!

WHERE (maskatu)
   tridcfu(:,:,:,4) = tridcfu(:,:,:,4) + uvel_A(1:ncloc,1:nrloc,:) + sourceatu
END WHERE

!
!3.2.3 Update uvel_A
!-------------------
!

IF (vimplicit) THEN
   CALL tridiag_vert(tridcfu,uvel_A,nhdims,nz,1,'U  ')
ELSE
   k_323: DO k=1,nz
      WHERE (maskatu(:,:,k))
         uvel_A(1:ncloc,1:nrloc,k) = tridcfu(:,:,k,4)/tridcfu(:,:,k,2)
      END WHERE
   ENDDO k_323
ENDIF

!
!4. Time integration with operator splitting - step B
!----------------------------------------------------
!
!4.1 Explicit scheme or first iteration
!--------------------------------------
!

IF (itsimp.EQ.1) THEN

!
!4.1.1 Initialise
!----------------
!

   tridcfu(:,:,:,1) = 0.0
   tridcfu(:,:,:,2) = 1.0
   tridcfu(:,:,:,3) = 0.0
   tridcfu(:,:,:,4) = 0.0

!
!4.1.2 Z-derivative terms
!------------------------
!
!  ---vertical advection
   IF (iopt_adv_3D.GT.0) THEN
      CALL Zadv_at_U(uvel_old,tridcfu,wadvatuw)
   ENDIF

!  ---vertical diffusion
   IF (iopt_vdif_coef.GT.0) THEN
      CALL Zdif_at_U(uvel_old,tridcfu,vdifcoefatuw,ibcsur,ibcbot,nbcs,nbcb,&
                   & bcsur,bcbot)
   ENDIF

!  ---sink term
   IF (fsink) THEN
      WHERE (maskatu)
         tridcfu(:,:,:,2) = tridcfu(:,:,:,2) + sinkatu
      END WHERE
   ENDIF
           
!
!4.1.3 Store at first iteration
!------------------------------
!

   IF (fimp) tridcfu_B_prev = tridcfu

ENDIF

!
!4.2 Explicit or next iteration
!------------------------------
!
!4.2.1 Retrieve from first iteration
!-----------------------------------
!

IF (fimp) tridcfu = tridcfu_B_prev
 
!
!4.2.2 Source term
!---------------
!

WHERE (maskatu)
   tridcfu(:,:,:,4) = tridcfu(:,:,:,4) + uvel_old(1:ncloc,1:nrloc,:) + &
                    & sourceatu
END WHERE

!
!4.2.3 Update uvel_B
!-------------------
!

IF (vimplicit) THEN
   CALL tridiag_vert(tridcfu,uvel_B,nhdims,nz,1,'U  ')
ELSE
   WHERE (maskatu)
      uvel_B(1:ncloc,1:nrloc,:) = tridcfu(:,:,:,4)/tridcfu(:,:,:,2)
   END WHERE
ENDIF

!
!4.2.4 Y-derivative terms
!------------------------
!
!---initialise
WHERE (maskatu) tridcfu(:,:,:,4) = 0.0

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/0,0,nhpsi,nhpsi/)
   CALL exchange_mod(uvel_B,lbounds,nhexch,iarr_uvel)
ENDIF

!---horizontal advection
IF (iopt_adv_3D.GT.0) THEN
   CALL Yadv_at_U_3d(uvel_B,advdif3d,nhpsi,.FALSE.)
   WHERE (maskatu)
      tridcfu(:,:,:,4) = tridcfu(:,:,:,4) - advdif3d/delzatu(1:ncloc,:,:)
   END WHERE
ENDIF

!---horizontal diffusion
IF (iopt_hdif_3D.GT.0) THEN
   CALL Ydif_at_U_3d(uvel_B,advdif3d,.FALSE.)
   WHERE (maskatu)
      tridcfu(:,:,:,4) = tridcfu(:,:,:,4) + advdif3d/delzatu(1:ncloc,:,:)
   END WHERE
ENDIF

!---update uvel_B
WHERE (maskatu)
   uvel_B(1:ncloc,1:nrloc,:) = uvel_B(1:ncloc,1:nrloc,:) + &
                             & delt3d*tridcfu(:,:,:,4)
END WHERE

!
!4.2.5 X-derivative terms
!------------------------
!
!---time derivative
WHERE (maskatu) tridcfu(:,:,:,4) = 0.0

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/nhpsi,nhpsi,0,0/)
   CALL exchange_mod(uvel_B,lbounds,nhexch,iarr_uvel)
ENDIF

!---horizontal advection
IF (iopt_adv_3D.GT.0) THEN
   CALL Xadv_at_U_3d(uvel_B,advdif3d,nhpsi,.FALSE.)
   WHERE (maskatu)
      tridcfu(:,:,:,4) = tridcfu(:,:,:,4) - advdif3d/delzatu(1:ncloc,:,:)
   END WHERE
ENDIF

!---horizontal diffusion
IF (iopt_hdif_3D.GT.0) THEN
   CALL Xdif_at_U_3d(uvel_B,advdif3d,.FALSE.)
   WHERE (maskatu)
      tridcfu(:,:,:,4) = tridcfu(:,:,:,4) + advdif3d/delzatu(1:ncloc,:,:)
   END WHERE
ENDIF

!---update uvel_B
WHERE (maskatu)
      uvel_B(1:ncloc,1:nrloc,:) = uvel_B(1:ncloc,1:nrloc,:) + &
                                & delt3d*tridcfu(:,:,:,4)
END WHERE

!
!5. Time integration with operator splitting - update uvel
!---------------------------------------------------------
!

WHERE (maskatu)
   uvel(1:ncloc,1:nrloc,:) = 0.5*(uvel_A(1:ncloc,1:nrloc,:)+&
                                & uvel_B(1:ncloc,1:nrloc,:))
END WHERE

1000 CONTINUE

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (advdif3d,maskatu,sourceatu,tridcfu)
IF (itsimp.EQ.1) THEN
   DEALLOCATE (sinkatu)
   IF (iopt_vdif_coef.GT.0) DEALLOCATE (vdifcoefatuw)
   IF (iopt_adv_3D.GT.0) DEALLOCATE (wadvatuw)
   IF (iopt_adv_3D.GT.0.AND.iopt_waves_curr.EQ.1) DEALLOCATE (wstokesatuw)
ENDIF
IF (time_split) DEALLOCATE (uvel_A,uvel_B)

CALL log_timer_out(npcc,itm_trans)


RETURN

END SUBROUTINE transport_at_U_3d

!========================================================================

SUBROUTINE transport_at_U_2d(source,sink)
!************************************************************************
!
! *transport_at_U_2d* Solve advection-diffusion equation for the 2-D U-current
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Transport_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - current_2d
!
! External calls - Xadv_at_U_2d, Xdif_at_U_2d, Yadv_at_U_2d, Ydif_at_U_2d
!
! Module calls - error_alloc, exchange_mod, num_halo
!
!************************************************************************
!
USE currents
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: num_halo

IMPLICIT NONE

!
!*  Arguments
!
REAL, INTENT(IN), DIMENSION(ncloc,nrloc) :: source, sink

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*source*   REAL    Source terms in transport equations           [m^2/s^2]
!*sink*     REAL    Sink terms arising from bottom stress             [1/s]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fimp, time_split, xadvflag, xdifflag, yadvflag, ydifflag
INTEGER :: nhpsi2d, npcc
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: advdif2d, dexp, dexp_prev, &
                           & udvel_A_prev, sinkatu, sourceatu, udvel_A, udvel_B

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*time_split* LOGICAL .TRUE. for directional splitting in time
!*dexp*       REAL    Explicit terms on r.h.s.
!*udvel_A*    REAL    Integrated U-current during first integration step
!                                                                       [m^2/s]
!*udvel_B*    REAL    Integrated U-current during second integration step
!                                                                       [m^2/s]
!
!------------------------------------------------------------------------------
!

procname(pglev+1) = 'transport_at_U_2d'
CALL log_timer_in(npcc,ivarid=iarr_udvel)

!
!1. Initialise
!-------------
!
!1.1 Halo size for exchange
!--------------------------
!

nhpsi2d = MAX(num_halo(iopt_adv_2D),iopt_hdif_2D)

!
!1.2 Flags
!---------
!

time_split = iopt_adv_2D.EQ.3
fimp = iopt_hydro_impl.EQ.1.AND.maxitsimp.NE.1

xadvflag = .FALSE.; yadvflag = .FALSE.; xdifflag = .FALSE.; ydifflag = .FALSE.
IF (iopt_hydro_impl.EQ.0.AND.iopt_grid_nodim.EQ.3.AND.predstep) THEN
   IF (iopt_adv_2D.GT.0) THEN
      xadvflag = .TRUE.; yadvflag = .TRUE.
   ENDIF
   IF (iopt_hdif_2D.GT.0) THEN
      xdifflag = .TRUE.; ydifflag = .TRUE.
   ENDIF
ENDIF

!
!1.3 Allocate arrays
!-------------------
!

IF (fimp.AND.itsimp.EQ.1) THEN
   IF (.NOT.time_split) THEN
      IF (ALLOCATED(dexp_prev)) DEALLOCATE (dexp_prev)
      ALLOCATE (dexp_prev(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('dexp_prev',2,(/ncloc,nrloc/),kndrtype)
   ELSE
      IF (ALLOCATED(udvel_A_prev)) DEALLOCATE (udvel_A_prev)
      ALLOCATE (udvel_A_prev(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),&
              & STAT=errstat)
      CALL error_alloc('udvel_A_prev',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),&
                     & kndrtype)
   ENDIF
ENDIF

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (advdif2d(ncloc,nrloc),STAT=errstat)
CALL error_alloc('advdif2d',3,(/ncloc,nrloc/),kndrtype)
ALLOCATE (dexp(ncloc,nrloc),STAT=errstat)
CALL error_alloc('dexp',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sinkatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sinkatu',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sourceatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sourceatu',2,(/ncloc,nrloc/),kndrtype)
IF (time_split) THEN
   ALLOCATE (udvel_A(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
   CALL error_alloc('udvel_A',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
   ALLOCATE (udvel_B(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
   CALL error_alloc('udvel_B',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
ENDIF

!
!1.4 Masks
!---------
!

maskatu = node2du(1:ncloc,1:nrloc).EQ.2

!
!1.5 Local arrays
!----------------
!

IF (time_split) THEN
   udvel_A = udvel
   udvel_B = udvel
ENDIF

!
!1.6 Source and sink terms
!-------------------------
!

WHERE (maskatu)
   sourceatu = delt2d*source
   sinkatu = 1.0 + delt2d*sink
END WHERE

!
!2. Time integration without operator splitting
!----------------------------------------------
!

IF (.NOT.time_split) THEN

!
!2.1 Explicit scheme or first iteration
!--------------------------------------
!

   IF (itsimp.EQ.1) THEN

!
!2.1.1 Initialise
!----------------
!

      WHERE (maskatu) dexp = 0.0

!
!2.1.2 Advection
!---------------
!

      IF (iopt_adv_2D.GT.0) THEN
         CALL Xadv_at_U_2d(udvel_old,advdif2d,nhpsi2d,xadvflag)
         WHERE (maskatu) dexp = dexp - advdif2d
         CALL Yadv_at_U_2d(udvel_old,advdif2d,nhpsi2d,yadvflag)
         WHERE (maskatu) dexp = dexp - advdif2d
      ENDIF

!
!2.1.3 Diffusion
!---------------
!

      IF (iopt_hdif_2D.GT.0) THEN
         CALL Xdif_at_U_2d(udvel_old,advdif2d,xdifflag)
         WHERE (maskatu) dexp = dexp + advdif2d
         CALL Ydif_at_U_2d(udvel_old,advdif2d,ydifflag)
         WHERE (maskatu) dexp = dexp + advdif2d
      ENDIF

!
!2.1.4 Explicit terms
!--------------------
!

      IF (iopt_hydro_impl.EQ.1.OR.iopt_grid_nodim.EQ.2) THEN
         WHERE (maskatu)
            dexp = udvel_old(1:ncloc,1:nrloc) + delt2d*dexp
         END WHERE
      ELSE
         WHERE (maskatu)
            dexp = udvel_old(1:ncloc,1:nrloc) + delt2d*(dexp+udevint)
         END WHERE
      ENDIF
           
!
!2.1.5 Store at first iteration
!------------------------------
!

      IF (fimp) dexp_prev = dexp

   ENDIF

!
!2.2 Explicit or next iteration
!------------------------------
!
!2.2.1 Retrieve from first iteration
!-----------------------------------
!

   IF (fimp) dexp = dexp_prev

!
!2.2.2 Source term
!-----------------
!

   WHERE (maskatu)
      dexp = dexp + sourceatu
   END WHERE

!
!2.2.3 Update udvel
!------------------
!

   WHERE (maskatu)
      udvel(1:ncloc,1:nrloc) = dexp/sinkatu
   END WHERE
   WHERE (node2du(1:ncloc,1:nrloc).EQ.1)
      udvel(1:ncloc,1:nrloc) = 0.0
   END WHERE

   GOTO 1000

ENDIF

!
!3. Time integration with operator splitting - step A
!----------------------------------------------------
!
!3.1 Explicit scheme or first iteration
!--------------------------------------
!

IF (itsimp.EQ.1) THEN

!
!3.1.1 X-derivative terms
!------------------------
!
!  ---initialise
   WHERE (maskatu) dexp = 0.0

!  ---advection
   IF (iopt_adv_2D.GT.0) THEN
      CALL Xadv_at_U_2d(udvel_old,advdif2d,nhpsi2d,xadvflag)
      WHERE (maskatu) dexp = dexp - advdif2d
   ENDIF

!  ---diffusion
   IF (iopt_hdif_2D.GT.0) THEN
      CALL Xdif_at_U_2d(udvel_old,advdif2d,xdifflag)
      WHERE (maskatu) dexp = dexp + advdif2d
   ENDIF

!  ---update udvel_A
   WHERE (maskatu) 
      udvel_A(1:ncloc,1:nrloc) = udvel_old(1:ncloc,1:nrloc) + delt2d*dexp
   END WHERE

!
!3.1.2 Y-derivative terms
!------------------------
!
!  ---initialise
   WHERE (maskatu) dexp = 0.0

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = 1-nhalo; nhexch = (/0,0,nhpsi2d,nhpsi2d/)
      CALL exchange_mod(udvel_A,lbounds,nhexch,iarr_udvel)
   ENDIF

!  ---advection
   IF (iopt_adv_2D.GT.0) THEN
      If (yadvflag) CALL Yadv_at_U_2d(udvel_old,advdif2d,nhpsi2d,yadvflag)
      CALL Yadv_at_U_2d(udvel_A,advdif2d,nhpsi2d,yadvflag)
      WHERE (maskatu) dexp = dexp - advdif2d
   ENDIF

!  ---diffusion
   IF (iopt_hdif_2D.GT.0) THEN
      IF (ydifflag) THEN
         CALL Ydif_at_U_2d(udvel_old,advdif2d,ydifflag)
      ENDIF
      CALL Ydif_at_U_2d(udvel_A,advdif2d,ydifflag)
      WHERE (maskatu) dexp = dexp + advdif2d
   ENDIF

!  ---update udvel_A
   WHERE (maskatu) 
      udvel_A(1:ncloc,1:nrloc) = udvel_A(1:ncloc,1:nrloc) + delt2d*dexp
   END WHERE
           
!
!3.1.3 Store at first iteration
!------------------------------
!

   IF (fimp) udvel_A_prev = udvel_A

ENDIF

!
!3.2 Explicit or next iteration
!------------------------------
!
!3.2.1 Retrieve from first iteration
!-----------------------------------
!

IF (fimp) udvel_A = udvel_A_prev

!
!3.2.2 Source/sink terms
!-----------------------
!

IF (iopt_hydro_impl.EQ.1.OR.iopt_grid_nodim.EQ.2) THEN
   WHERE (maskatu)
      udvel_A(1:ncloc,1:nrloc) = (udvel_A(1:ncloc,1:nrloc)+sourceatu)/sinkatu
   END WHERE
ELSE
   WHERE (maskatu)
      udvel_A(1:ncloc,1:nrloc) = (udvel_A(1:ncloc,1:nrloc)+sourceatu+&
                                & delt2d*udevint)/sinkatu
   END WHERE
ENDIF

!
!4. Time integration with operator splitting - step B
!----------------------------------------------------
!
!4.1 Source/sink terms
!---------------------
!

IF (iopt_hydro_impl.EQ.1.OR.iopt_grid_nodim.EQ.2) THEN
   WHERE (maskatu)
      udvel_B(1:ncloc,1:nrloc) = (udvel_old(1:ncloc,1:nrloc)+sourceatu)/sinkatu
   END WHERE
ELSE
   WHERE (maskatu)
      udvel_B(1:ncloc,1:nrloc) = (udvel_old(1:ncloc,1:nrloc)+sourceatu+&
                                & delt2d*udevint)/sinkatu
   END WHERE
ENDIF

!
!4.2 Y-derivative terms
!----------------------
!
!---initialise
WHERE (maskatu) dexp = 0.0

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/0,0,nhpsi2d,nhpsi2d/)
   CALL exchange_mod(udvel_B,lbounds,nhexch,iarr_udvel)
ENDIF

!---advection
IF (iopt_adv_2D.GT.0) THEN
   CALL Yadv_at_U_2d(udvel_B,advdif2d,nhpsi2d,yadvflag)
   WHERE (maskatu) dexp = dexp - advdif2d
ENDIF

!---diffusion
IF (iopt_hdif_2D.GT.0) THEN
   CALL Ydif_at_U_2d(udvel_B,advdif2d,ydifflag)
   WHERE (maskatu) dexp = dexp + advdif2d
ENDIF

!---update udvel_B
WHERE (maskatu) 
   udvel_B(1:ncloc,1:nrloc) = udvel_B(1:ncloc,1:nrloc) + delt2d*dexp
END WHERE

!
!4.3 X-derivative terms
!----------------------
!
!---initialise
WHERE (maskatu) dexp = 0.0

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/nhpsi2d,nhpsi2d,0,0/)
   CALL exchange_mod(udvel_B,lbounds,nhexch,iarr_udvel)
ENDIF

!---advection
IF (iopt_adv_2D.GT.0) THEN
   CALL Xadv_at_U_2d(udvel_B,advdif2d,nhpsi2d,xadvflag)
   WHERE (maskatu) dexp = dexp - advdif2d
ENDIF

!---diffusion
IF (iopt_hdif_2D.GT.0) THEN
   CALL Xdif_at_U_2d(udvel_B,advdif2d,xdifflag)
   WHERE (maskatu) dexp = dexp + advdif2d
ENDIF

!---update udvel_B
WHERE (maskatu)
   udvel_B(1:ncloc,1:nrloc) = udvel_B(1:ncloc,1:nrloc) + delt2d*dexp
END WHERE

!
!5. Time integration with operator splitting - update udvel
!----------------------------------------------------------
!

WHERE (maskatu) 
   udvel(1:ncloc,1:nrloc) = 0.5*(udvel_A(1:ncloc,1:nrloc)+&
                               & udvel_B(1:ncloc,1:nrloc))
END WHERE

1000 CONTINUE

!
!6. Deallocate arrays
!--------------------
!
      
DEALLOCATE (maskatu,advdif2d,dexp,sinkatu,sourceatu)
IF (time_split) DEALLOCATE (udvel_A,udvel_B)

CALL log_timer_out(npcc,itm_trans)


RETURN

END SUBROUTINE transport_at_U_2d

!========================================================================

SUBROUTINE transport_at_V_3d(source,sink,ibcsur,ibcbot,nbcs,nbcb,bcsur,bcbot)
!************************************************************************
!
! *transport_at_V_3d* Solve advection-diffusion equation for the 3-D V-current
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Transport_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - current_pred
!
! External calls - Xadv_at_V_3d, Xdif_at_V_3d, Yadv_at_V_3d, Ydif_at_V_3d,
!                  Zadv_at_V, Zdif_at_V
!
! Module calls - error_alloc, exchange_mod, num_halo, tridiag_vert, Warr_at_VW
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
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Warr_at_VW
USE error_routines, ONLY: error_alloc
USE nla_library, ONLY: tridiag_vert
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: num_halo

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ibcsur, ibcbot, nbcs, nbcb
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz) :: sink, source
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcs) :: bcsur
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nbcb) :: bcbot

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*source*     REAL    Source terms in transport equations             [m/s^2]
!*sink*       REAL    Sink terms in transport equations                 [1/s]
!             (:,:,1) => prescribed surface flux or surface value
!             (:,:,2) => transfer velocity                              [m/s]
!*ibcsur*     INTEGER Type of surface boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!                 = 3 => Dirichlet at first C-node below the surface
!                 = 4 => Dirichlet at the surface
!*ibcbot*     INTEGER Type of bottom boundary condition
!                 = 0 => Neumann (zero flux)
!                 = 1 => Neumann (prescibed flux)
!                 = 2 => Neumann (using transfer velocity)
!                 = 3 => Dirichlet at first C-node above the bottom
!                 = 4 => Dirichlet at the bottom
!*nbcs*       INTEGER Last dimension of array bcsur
!*nbcb*       INTEGER Last dimension of array bcbot
!*bcsur*      REAL    Data for surface boundary condition
!             (:,:,1) => prescribed surface flux or surface value
!             (:,:,2) => transfer velocity                              [m/s]
!*bcbot*      REAL    Data for bottom boundary condition
!             (:,:,1) => prescribed bottom flux or bottom value
!             (:,:,2) => transfer velocity                              [m/s]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fimp, fsink, time_split, vflag, vimplicit
INTEGER :: k, nhpsi, npcc
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhdims, nhexch
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: advdif3d, sinkatv, sourceatv, &
                              & vdifcoefatvw, vvel_A, vvel_A_prev, vvel_B, &
                              & wadvatvw, wstokesatvw
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: tridcfv, tridcfv_prev, &
                                 & tridcfv_A_prev, tridcfv_B_prev

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*time_split*  LOGICAL .TRUE. for directional splitting in time
!*vimplicit*   LOGICAL .TRUE. for implicit integration in the vertical
!*maskatv*     LOGICAL Mask array at V-nodes
!*tridcfv*     REAL    Tridiagonal matrix for implicit vertical solution
!*vvel_A*      REAL    V-current during first integration step        [m/s]
!*vvel_B*      REAL    V-current during second integration step       [m/s]  
!*hdifcoefatc* REAL    Horizontal diffusion coefficient at C-nodes  [m^2/s]
!*vdifcoefatvw*REAL    Vertical diffusion coefficient at VW-nodes   [m^2/s]
!*wadvatvw*    REAL    Advective velocity in Z-direction at VW-nodes  [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'transport_at_V_3d'
CALL log_timer_in(npcc,ivarid=iarr_vvel)

!
!1. Initialise
!-------------
!
!1.1 Array shape and halo size
!--------------------------
!

nhpsi = MAX(num_halo(iopt_adv_3D),iopt_hdif_3D)
nhdims = nhalo

!
!1.2 Time integration
!--------------------
!

fsink = iopt_weibar.EQ.1
vimplicit = ((iopt_adv_3D.GT.0.AND.iopt_vadv_impl.GT.0).OR.&
           & (iopt_vdif_coef.GT.0.AND.iopt_vdif_impl.GT.0).OR.&
           & (ibcsur.EQ.4).OR.(ibcbot.EQ.4))
time_split = iopt_adv_3D.EQ.3
vflag = iopt_grid_nodim.EQ.3.AND.iopt_hydro_impl.EQ.0 
fimp = iopt_hydro_impl.EQ.1.AND.maxitsimp.NE.1

!
!1.3 Allocate arrays
!-------------------
!

IF (fimp.AND.itsimp.EQ.1) THEN
   IF (.NOT.time_split) THEN
      IF (ALLOCATED(tridcfv_prev)) DEALLOCATE (tridcfv_prev)
      ALLOCATE (tridcfv_prev(ncloc,nrloc,nz,4),STAT=errstat)
      CALL error_alloc('tridcfv_prev',4,(/ncloc,nrloc,nz,4/),kndrtype)
   ELSE
      IF (ALLOCATED(tridcfv_A_prev)) THEN
         DEALLOCATE (tridcfv_A_prev,tridcfv_B_prev,vvel_A_prev)
      ENDIF
      ALLOCATE (tridcfv_A_prev(ncloc,nrloc,nz,4),STAT=errstat)
      CALL error_alloc('tridcfv_A_prev',4,(/ncloc,nrloc,nz,4/),kndrtype)
      ALLOCATE (tridcfv_B_prev(ncloc,nrloc,nz,4),STAT=errstat)
      CALL error_alloc('tridcfv_B_prev',4,(/ncloc,nrloc,nz,4/),kndrtype)
      ALLOCATE (vvel_A_prev(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),&
              & STAT=errstat)
      CALL error_alloc('vvel_A_prev',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),&
                     & kndrtype)
   ENDIF
ENDIF

ALLOCATE (advdif3d(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('advdif3d',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (maskatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (sourceatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('sourceatv',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (tridcfv(ncloc,nrloc,nz,4),STAT=errstat)
CALL error_alloc('tridcfv',4,(/ncloc,nrloc,nz,4/),kndrtype)
IF (itsimp.EQ.1) THEN
   ALLOCATE (sinkatv(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('sinkatv',3,(/ncloc,nrloc,nz/),kndrtype)
   IF (iopt_vdif_coef.GT.0) THEN
      ALLOCATE (vdifcoefatvw(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('vdifcoefatvw',3,(/ncloc,nrloc,nz+1/),kndrtype)
   ENDIF
   IF (iopt_adv_3D.GT.0) THEN
      ALLOCATE (wadvatvw(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('wadvatvw',3,(/ncloc,nrloc,nz+1/),kndrtype)
      IF (iopt_waves_curr.EQ.1) THEN
         ALLOCATE (wstokesatvw(ncloc,nrloc,nz+1),STAT=errstat)
         CALL error_alloc('wstokesatvw',3,(/ncloc,nrloc,nz+1/),kndrtype)
      ENDIF
   ENDIF
ENDIF
IF (time_split) THEN
   ALLOCATE (vvel_A(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
   CALL error_alloc('vvel_A',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
   ALLOCATE (vvel_B(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
   CALL error_alloc('vvel_B',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
ENDIF

!
!1.4 Masks
!---------
!

maskatv = nodeatv(1:ncloc,1:nrloc,:).EQ.2

!
!1.5 Local arrays
!----------------
!
!---depth-integrated term
vdevint = 0.0

!---temporary values of vvel
IF (time_split) THEN
   vvel_A = vvel; vvel_B = vvel
ENDIF

!
!1.6 Source/sink terms
!---------------------
!

WHERE (maskatv)
   sourceatv = delt3d*source
END WHERE
IF (itsimp.EQ.1.AND.fsink) THEN
   WHERE (maskatv)
      sinkatv = delt3d*sink
   END WHERE
ENDIF

!
!1.7 Advective currents
!----------------------
!

IF (itsimp.EQ.1.AND.iopt_adv_3D.GT.0) THEN
   CALL Warr_at_VW(wfvel(1:ncloc,:,:),wadvatvw,3,(/1,0,1/),&
                & (/ncloc,nrloc,nz+1/),1,iarr_wvel,.TRUE.)
   IF (iopt_waves_curr.EQ.1) THEN
      CALL Warr_at_VW(wstokesatw(1:ncloc,:,:),wstokesatvw,3,(/1,0,1/),&
                   & (/ncloc,nrloc,nz+1/),1,iarr_wstokesatw,.TRUE.)
      wadvatvw = wadvatvw + wstokesatvw 
   ENDIF
ENDIF

!
!1.8 Vertical diffusion coefficient
!----------------------------------
!

IF (itsimp.EQ.1.AND.iopt_vdif_coef.GT.0) THEN
   CALL Warr_at_VW(vdifcoefmom(1:ncloc,:,:),vdifcoefatvw,3,(/1,0,1/),&
                & (/ncloc,nrloc,nz+1/),1,iarr_vdifcoefmom,.TRUE.)
ENDIF

!
!2. Time integration without operator splitting
!----------------------------------------------
!

IF (.NOT.time_split) THEN

!
!2.1 Explicit scheme or first iteration
!--------------------------------------
!

   IF (itsimp.EQ.1) THEN

!
!2.1.1 Initialise
!----------------
!

      tridcfv(:,:,:,1) = 0.0
      tridcfv(:,:,:,2) = 1.0
      tridcfv(:,:,:,3) = 0.0
      tridcfv(:,:,:,4) = 0.0

!
!2.1.2 Time derivative (rhs part)
!--------------------------------
!

      WHERE (maskatv)
         tridcfv(:,:,:,4) = vvel_old(1:ncloc,1:nrloc,:)
      END WHERE

!
!2.1.3 Advection
!---------------
!

      IF (iopt_adv_3D.GT.0) THEN
         CALL Xadv_at_V_3d(vvel_old,advdif3d,nhpsi,vflag)
         WHERE (maskatv)
            tridcfv(:,:,:,4) = tridcfv(:,:,:,4) - &
                             & delt3d*advdif3d/delzatv(:,1:nrloc,:)
         END WHERE
         CALL Yadv_at_V_3d(vvel_old,advdif3d,nhpsi,vflag)
         WHERE (maskatv)
            tridcfv(:,:,:,4) = tridcfv(:,:,:,4) - &
                             & delt3d*advdif3d/delzatv(:,1:nrloc,:)
         END WHERE
         CALL Zadv_at_V(vvel_old,tridcfv,wadvatvw)
      ENDIF

!
!2.1.4 Diffusion
!---------------
!
!     ---horizontal
      IF (iopt_hdif_3D.GT.0) THEN
         CALL Xdif_at_V_3d(vvel_old,advdif3d,vflag)
         WHERE (maskatv)
            tridcfv(:,:,:,4) = tridcfv(:,:,:,4) + &
                             & delt3d*advdif3d/delzatv(:,1:nrloc,:)
         END WHERE
         CALL Ydif_at_V_3d(vvel_old,advdif3d,vflag)
         WHERE (maskatv)
            tridcfv(:,:,:,4) = tridcfv(:,:,:,4) + &
                             & delt3d*advdif3d/delzatv(:,1:nrloc,:)
         END WHERE
      ENDIF

!     ---vertical
      IF (iopt_vdif_coef.GT.0) THEN
         CALL Zdif_at_V(vvel_old,tridcfv,vdifcoefatvw,ibcsur,ibcbot,nbcs,nbcb,&
                      & bcsur,bcbot)
      ENDIF
            
!
!2.1.5 Sink term
!---------------
!

      IF (fsink) THEN
         WHERE (maskatv)
            tridcfv(:,:,:,2) = tridcfv(:,:,:,2) + sinkatv
         END WHERE
      ENDIF
           
!
!2.1.6 Store at first iteration
!------------------------------
!

      IF (fimp) tridcfv_prev = tridcfv

   ENDIF

!
!2.2 Explicit or next iteration
!------------------------------
!
!2.2.1 Retrieve from first iteration
!-----------------------------------
!

   IF (fimp) tridcfv = tridcfv_prev

!
!2.2.2 Source term
!-----------------
!

   WHERE (maskatv)
      tridcfv(:,:,:,4) = tridcfv(:,:,:,4) + sourceatv
   END WHERE

!
!2.2.3 Update vvel
!-----------------
!

   IF (vimplicit) THEN
      CALL tridiag_vert(tridcfv,vvel,nhdims,nz,1,'V  ')
   ELSE
      WHERE (maskatv)
         vvel(1:ncloc,1:nrloc,:) = tridcfv(:,:,:,4)/tridcfv(:,:,:,2)
      END WHERE
      WHERE (nodeatv(1:ncloc,1:nrloc,:).EQ.1)
         vvel(1:ncloc,1:nrloc,:) = 0.0
      END WHERE
   ENDIF

   GOTO 1000

ENDIF

!
!3. Time integration with operator splitting - step A
!----------------------------------------------------
!
!3.1 Explicit scheme or first iteration
!--------------------------------------
!

IF (itsimp.EQ.1) THEN

!
!3.1.1 Initialise
!----------------
!

   tridcfv(:,:,:,1) = 0.0
   tridcfv(:,:,:,2) = 1.0
   tridcfv(:,:,:,3) = 0.0
   tridcfv(:,:,:,4) = 0.0

!
!3.1.2 X-derivative terms
!------------------------
!
!  ---advection
   IF (iopt_adv_3D.GT.0) THEN
      CALL Xadv_at_V_3d(vvel_old,advdif3d,nhpsi,vflag)
      WHERE (maskatv)
         tridcfv(:,:,:,4) = tridcfv(:,:,:,4) - advdif3d/delzatv(:,1:nrloc,:)
      END WHERE
   ENDIF

!  ---diffusion
   IF (iopt_hdif_3D.GT.0) THEN
      CALL Xdif_at_V_3d(vvel_old,advdif3d,hdifcoef3datuv,vflag)
      WHERE (maskatv)
         tridcfv(:,:,:,4) = tridcfv(:,:,:,4) + advdif3d/delzatv(:,1:nrloc,:)
      END WHERE
   ENDIF

!  ---update vvel_A
   WHERE (maskatv)
      vvel_A(1:ncloc,1:nrloc,:) = vvel_old(1:ncloc,1:nrloc,:) + &
                                & delt3d*tridcfv(:,:,:,4)
   END WHERE

!
!3.1.3 Y-derivative terms
!------------------------
!
!  ---initialise
   WHERE (maskatv) tridcfv(:,:,:,4) = 0.0

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/1-nhalo,1-nhalo,1/); nhexch = (/0,0,nhpsi,nhpsi/)
      CALL exchange_mod(vvel_A,lbounds,nhexch,iarr_vvel)
   ENDIF

!  ---horizontal advection
   IF (iopt_adv_3D.GT.0) THEN
      IF (vflag) CALL Yadv_at_V_3d(vvel_old,advdif3d,nhpsi,.TRUE.)
      CALL Yadv_at_V_3d(vvel_A,advdif3d,nhpsi,.FALSE.)
      WHERE (maskatv)
         tridcfv(:,:,:,4) = tridcfv(:,:,:,4) - advdif3d/delzatv(:,1:nrloc,:)
      END WHERE
   ENDIF

!  ---horizontal diffusion
   IF (iopt_hdif_3D.GT.0) THEN
      IF (vflag) CALL Ydif_at_V_3d(vvel_old,advdif3d,.TRUE.)
      CALL Ydif_at_V_3d(vvel_A,advdif3d,.FALSE.)
      WHERE (maskatv)
         tridcfv(:,:,:,4) = tridcfv(:,:,:,4) + advdif3d/delzatv(:,1:nrloc,:)
      END WHERE
   ENDIF

!  ---update vvel_A
   WHERE (maskatv)
      vvel_A(1:ncloc,1:nrloc,:) = vvel_A(1:ncloc,1:nrloc,:) + &
                                & delt3d*tridcfv(:,:,:,4)
   END WHERE

!
!3.1.4 Z-derivative terms
!------------------------
!
!  ---initialise
   WHERE (maskatv) tridcfv(:,:,:,4) = 0.0

!  ---vertical advection
   IF (iopt_adv_3D.GT.0) THEN
      CALL Zadv_at_V(vvel_A,tridcfv,wadvatvw)
   ENDIF

!  ---vertical diffusion
   IF (iopt_vdif_coef.GT.0) THEN
      CALL Zdif_at_V(vvel_A,tridcfv,vdifcoefatvw,ibcsur,ibcbot,nbcs,nbcb,&
                   & bcsur,bcbot)
   ENDIF

!  ---sink term
   IF (fsink) THEN
      WHERE (maskatv)
         tridcfv(:,:,:,2) = tridcfv(:,:,:,2) + sinkatv
      END WHERE
   ENDIF
           
!
!3.1.5 Store at first iteration
!------------------------------
!

   IF (fimp) THEN
      tridcfv_A_prev = tridcfv
      vvel_A_prev = vvel_A
   ENDIF

ENDIF

!
!3.2 Explicit or next iteration
!------------------------------
!
!3.2.1 Retrieve from first iteration
!-----------------------------------
!

IF (fimp) THEN
   tridcfv = tridcfv_A_prev
   vvel_A = vvel_A_prev
ENDIF
 
!
!3.2.2 Source term and time derivative (rhs part)
!------------------------------------------------
!

WHERE (maskatv)
   tridcfv(:,:,:,4) = tridcfv(:,:,:,4) + vvel_A(1:ncloc,1:nrloc,:) + sourceatv
END WHERE

!
!3.2.3 Update vvel_A
!-------------------
!

IF (vimplicit) THEN
   CALL tridiag_vert(tridcfv,vvel_A,nhdims,nz,1,'V  ')
ELSE
   k_323: DO k=1,nz
      WHERE (maskatv(:,:,k))
         vvel_A(1:ncloc,1:nrloc,k) = tridcfv(:,:,k,4)/tridcfv(:,:,k,2)
      END WHERE
   ENDDO k_323
ENDIF

!
!4. Time integration with operator splitting - step B
!----------------------------------------------------
!
!4.1 Explicit scheme or first iteration
!--------------------------------------
!

IF (itsimp.EQ.1) THEN

!
!4.1.1 Initialise
!----------------
!

   tridcfv(:,:,:,1) = 0.0
   tridcfv(:,:,:,2) = 1.0
   tridcfv(:,:,:,3) = 0.0
   tridcfv(:,:,:,4) = 0.0

!
!4.1.2 Z-derivative terms
!------------------------
!
!  ---vertical advection
   IF (iopt_adv_3D.GT.0) THEN
      CALL Zadv_at_V(vvel_old,tridcfv,wadvatvw)
   ENDIF

!  ---vertical diffusion
   IF (iopt_vdif_coef.GT.0) THEN
      CALL Zdif_at_V(vvel_old,tridcfv,vdifcoefatvw,ibcsur,ibcbot,nbcs,nbcb,&
                   & bcsur,bcbot)
   ENDIF

!  ---sink term
   IF (fsink) THEN
      WHERE (maskatv)
         tridcfv(:,:,:,2) = tridcfv(:,:,:,2) + sinkatv
      END WHERE
   ENDIF
           
!
!4.1.3 Store at first iteration
!------------------------------
!

   IF (fimp) tridcfv_B_prev = tridcfv

ENDIF

!
!4.2 Explicit or next iteration
!------------------------------
!
!4.2.1 Retrieve from first iteration
!-----------------------------------
!

IF (fimp) tridcfv = tridcfv_B_prev
 
!
!4.2.2 Source term
!-----------------
!

WHERE (maskatv)
   tridcfv(:,:,:,4) = tridcfv(:,:,:,4) + vvel_old(1:ncloc,1:nrloc,:) + sourceatv
END WHERE

!
!4.2.3 Update vvel_B
!-------------------
!

IF (vimplicit) THEN
   CALL tridiag_vert(tridcfv,vvel_B,nhdims,nz,1,'V  ')
ELSE
   WHERE (maskatv)
      vvel_B(1:ncloc,1:nrloc,:) = tridcfv(:,:,:,4)/tridcfv(:,:,:,2)
   END WHERE
ENDIF

!
!4.2.4 Y-derivative terms
!------------------------
!
!---initialise
WHERE (maskatv) tridcfv(:,:,:,4) = 0.0

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/0,0,nhpsi,nhpsi/)
   CALL exchange_mod(vvel_B,lbounds,nhexch,iarr_vvel)
ENDIF

!---horizontal advection
IF (iopt_adv_3D.GT.0) THEN
   CALL Yadv_at_V_3d(vvel_B,advdif3d,nhpsi,.FALSE.)
   WHERE (maskatv)
      tridcfv(:,:,:,4) = tridcfv(:,:,:,4) - advdif3d/delzatv(:,1:nrloc,:)
   END WHERE
ENDIF

!---horizontal diffusion
IF (iopt_hdif_3D.GT.0) THEN
   CALL Ydif_at_V_3d(vvel_B,advdif3d,.FALSE.)
   WHERE (maskatv)
      tridcfv(:,:,:,4) = tridcfv(:,:,:,4) + advdif3d/delzatv(:,1:nrloc,:)
   END WHERE
ENDIF

!---update vvel_B
WHERE (maskatv)
   vvel_B(1:ncloc,1:nrloc,:) = vvel_B(1:ncloc,1:nrloc,:) + &
                             & delt3d*tridcfv(:,:,:,4)
END WHERE

!
!4.2.5 X-derivative terms
!------------------------
!
!---time derivative
WHERE (maskatv) tridcfv(:,:,:,4) = 0.0

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/nhpsi,nhpsi,0,0/)
   CALL exchange_mod(vvel_B,lbounds,nhexch,iarr_vvel)
ENDIF

!---horizontal advection
IF (iopt_adv_3D.GT.0) THEN
   CALL Xadv_at_V_3d(vvel_B,advdif3d,nhpsi,.FALSE.)
   WHERE (maskatv)
      tridcfv(:,:,:,4) = tridcfv(:,:,:,4) - advdif3d/delzatv(:,1:nrloc,:)
   END WHERE
ENDIF

!---horizontal diffusion
IF (iopt_hdif_3D.GT.0) THEN
   CALL Xdif_at_V_3d(vvel_B,advdif3d,.FALSE.)
   WHERE (maskatv)
      tridcfv(:,:,:,4) = tridcfv(:,:,:,4) + advdif3d/delzatv(:,1:nrloc,:)
   END WHERE
ENDIF

!---update vvel_B
WHERE (maskatv)
    vvel_B(1:ncloc,1:nrloc,:) = vvel_B(1:ncloc,1:nrloc,:) + &
                              & delt3d*tridcfv(:,:,:,4)
END WHERE

!
!5. Time integration with operator splitting - update vvel
!---------------------------------------------------------
!

WHERE (maskatv)
   vvel(1:ncloc,1:nrloc,:) = 0.5*(vvel_A(1:ncloc,1:nrloc,:)+&
                                & vvel_B(1:ncloc,1:nrloc,:))
END WHERE

1000 CONTINUE

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (advdif3d,maskatv,sourceatv,tridcfv)
IF (itsimp.EQ.1) THEN
   DEALLOCATE (sinkatv)
   IF (iopt_vdif_coef.GT.0) DEALLOCATE (vdifcoefatvw)
   IF (iopt_adv_3D.GT.0) DEALLOCATE (wadvatvw)
   IF (iopt_adv_3D.GT.0.AND.iopt_waves_curr.EQ.1) DEALLOCATE (wstokesatvw)
ENDIF
IF (time_split) DEALLOCATE (vvel_A,vvel_B)

CALL log_timer_out(npcc,itm_trans)


RETURN

END SUBROUTINE transport_at_V_3d

!========================================================================

SUBROUTINE transport_at_V_2d(source,sink)
!************************************************************************
!
! *transport_at_V_2d* Solve advection-diffusion equation for the 2-D V-current
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Transport_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - current_2d
!
! External calls - Xadv_at_V_2d, Xdif_at_V_2d, Yadv_at_V_2d, Ydif_at_V_2d
!
! Module calls - error_alloc, exchange_mod, num_halo
!
!************************************************************************
!
USE currents
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: num_halo

IMPLICIT NONE

!
!*  Arguments
!
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc) :: source, sink

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*source*   REAL    Source terms in transport equations           [m^2/s^2]
!*sink*     REAL    Sink terms arising from bottom stress             [1/s]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fimp, time_split, xadvflag, xdifflag, yadvflag, ydifflag
INTEGER :: nhpsi2d, npcc
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: advdif2d, dexp, dexp_prev, sinkatv, &
                                    & sourceatv, vdvel_A, vdvel_A_prev, vdvel_B

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*time_split* LOGICAL .TRUE. for directional splitting in time
!*dexp*       REAL    Explicit terms on r.h.s.
!*vdvel_A*    REAL    Integrated U-current during first integration step
!                                                                       [m^2/s]
!*vdvel_B*    REAL    Integrated U-current during second integration step
!                                                                       [m^2/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'transport_at_V_2d'
CALL log_timer_in(npcc,ivarid=iarr_vdvel)

!
!1. Initialise
!-------------
!
!1.1 Halo size for exchange
!--------------------------
!

nhpsi2d = MAX(num_halo(iopt_adv_2D),iopt_hdif_2D)

!
!1.2 Flags
!---------
!

time_split = iopt_adv_2D.EQ.3
fimp = iopt_hydro_impl.EQ.1.AND.maxitsimp.NE.1

xadvflag = .FALSE.; yadvflag = .FALSE.; xdifflag = .FALSE.; ydifflag = .FALSE.
IF (iopt_hydro_impl.EQ.0.AND.iopt_grid_nodim.EQ.3.AND.predstep) THEN
   IF (iopt_adv_2D.GT.0) THEN
      xadvflag = .TRUE.; yadvflag = .TRUE.
   ENDIF
   IF (iopt_hdif_2D.GT.0) THEN
      xdifflag = .TRUE.; ydifflag = .TRUE.
   ENDIF
ENDIF

!
!1.3 Allocate arrays
!-------------------
!

IF (fimp.AND.itsimp.EQ.1) THEN
   IF (.NOT.time_split) THEN
      IF (ALLOCATED(dexp_prev)) DEALLOCATE (dexp_prev)
      ALLOCATE (dexp_prev(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('dexp_prev',2,(/ncloc,nrloc/),kndrtype)
   ELSE
      IF (ALLOCATED(vdvel_A_prev)) DEALLOCATE (vdvel_A_prev)
      ALLOCATE (vdvel_A_prev(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),&
              & STAT=errstat)
      CALL error_alloc('vdvel_A_prev',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),&
                     & kndrtype)
   ENDIF
ENDIF

ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (advdif2d(ncloc,nrloc),STAT=errstat)
CALL error_alloc('advdif2d',3,(/ncloc,nrloc/),kndrtype)
ALLOCATE (dexp(ncloc,nrloc),STAT=errstat)
CALL error_alloc('dexp',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sinkatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sinkatv',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sourceatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sourceatv',2,(/ncloc,nrloc/),kndrtype)
IF (time_split) THEN
   ALLOCATE (vdvel_A(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
   CALL error_alloc('vdvel_A',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
   ALLOCATE (vdvel_B(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
   CALL error_alloc('vdvel_B',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
ENDIF

!
!1.4 Masks
!---------
!

maskatv = node2dv(1:ncloc,1:nrloc).EQ.2

!
!1.5 Local arrays
!----------------
!

IF (time_split) THEN
   vdvel_A = vdvel
   vdvel_B = vdvel
ENDIF

!
!1.6 Source and sink terms
!-------------------------
!

WHERE (maskatv)
   sourceatv = delt2d*source
   sinkatv = 1.0 + delt2d*sink
END WHERE

!
!2. Time integration without operator splitting
!----------------------------------------------
!

IF (.NOT.time_split) THEN

!
!2.1 Explicit scheme or first iteration
!--------------------------------------
!

   IF (itsimp.EQ.1) THEN

!
!2.1.1 Initialise
!----------------
!

      WHERE (maskatv) dexp = 0.0

!
!2.1.2 Advection
!---------------
!

      IF (iopt_adv_2D.GT.0) THEN
         CALL Xadv_at_V_2d(vdvel_old,advdif2d,nhpsi2d,xadvflag)
         WHERE (maskatv) dexp = dexp - advdif2d
         CALL Yadv_at_V_2d(vdvel_old,advdif2d,nhpsi2d,yadvflag)
         WHERE (maskatv) dexp = dexp - advdif2d
      ENDIF

!
!2.1.3 Diffusion
!---------------
!

      IF (iopt_hdif_2D.GT.0) THEN
         CALL Xdif_at_V_2d(vdvel_old,advdif2d,xdifflag)
         WHERE (maskatv) dexp = dexp + advdif2d
         CALL Ydif_at_V_2d(vdvel_old,advdif2d,ydifflag)
         WHERE (maskatv) dexp = dexp + advdif2d
      ENDIF

!
!2.1.4 Explicit terms
!--------------------
!

      IF (iopt_hydro_impl.EQ.1.OR.iopt_grid_nodim.EQ.2) THEN
         WHERE (maskatv)
            dexp = vdvel_old(1:ncloc,1:nrloc) + delt2d*dexp
         END WHERE
      ELSE
         WHERE (maskatv)
            dexp = vdvel_old(1:ncloc,1:nrloc) + delt2d*(dexp+vdevint)
         END WHERE
      ENDIF
           
!
!2.1.5 Store at first iteration
!------------------------------
!

      IF (fimp) dexp_prev = dexp

   ENDIF

!
!2.2 Explicit or next iteration
!------------------------------
!
!2.2.1 Retrieve from first iteration
!-----------------------------------
!

   IF (fimp) dexp = dexp_prev

!
!2.2.2 Source term
!-----------------
!

   WHERE (maskatv)
      dexp = dexp + sourceatv
   END WHERE

!
!2.2.3 Update vdvel
!----------------
!

   WHERE (maskatv)
      vdvel(1:ncloc,1:nrloc) = dexp/sinkatv
   END WHERE
   WHERE (node2dv(1:ncloc,1:nrloc).EQ.1)
      vdvel(1:ncloc,1:nrloc) = 0.0
   END WHERE

   GOTO 1000

ENDIF

!
!3. Time integration with operator splitting - step A
!----------------------------------------------------
!
!
!3.1 Explicit scheme or first iteration
!--------------------------------------
!

IF (itsimp.EQ.1) THEN

!3.1 X-derivative terms
!----------------------
!
!
!  ---initialise
   WHERE (maskatv) dexp = 0.0

!  ---advection
   IF (iopt_adv_2D.GT.0) THEN
      CALL Xadv_at_V_2d(vdvel_old,advdif2d,nhpsi2d,xadvflag)
      WHERE (maskatv) dexp = dexp - advdif2d
   ENDIF

!  ---diffusion
   IF (iopt_hdif_2D.GT.0) THEN
      CALL Xdif_at_V_2d(vdvel_old,advdif2d,xdifflag)
      WHERE (maskatv) dexp = dexp + advdif2d
   ENDIF
   
!  ---update vdvel_A
   WHERE (maskatv) 
      vdvel_A(1:ncloc,1:nrloc) = vdvel_old(1:ncloc,1:nrloc) + delt2d*dexp
   END WHERE

!
!3.1.2 Y-derivative terms
!------------------------
!
!  ---initialise
   WHERE (maskatv) dexp = 0.0

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = 1-nhalo; nhexch = (/0,0,nhpsi2d,nhpsi2d/)
      CALL exchange_mod(vdvel_A,lbounds,nhexch,iarr_vdvel)
   ENDIF

!  ---advection
   IF (iopt_adv_2D.GT.0) THEN
      IF (yadvflag) CALL Yadv_at_V_2d(vdvel_old,advdif2d,nhpsi2d,yadvflag)
      CALL Yadv_at_V_2d(vdvel_A,advdif2d,nhpsi2d,yadvflag)
      WHERE (maskatv) dexp = dexp - advdif2d
   ENDIF

!  ---diffusion
   IF (iopt_hdif_2D.GT.0) THEN
      IF (ydifflag) THEN
         CALL Ydif_at_V_2d(vdvel_old,advdif2d,ydifflag)
      ENDIF
      CALL Ydif_at_V_2d(vdvel_A,advdif2d,ydifflag)
      WHERE (maskatv) dexp = dexp + advdif2d
   ENDIF

!  ---update vdvel_A
   WHERE (maskatv) 
      vdvel_A(1:ncloc,1:nrloc) = vdvel_A(1:ncloc,1:nrloc) + delt2d*dexp
   END WHERE
           
!
!3.1.3 Store at first iteration
!------------------------------
!

   IF (fimp) vdvel_A_prev = vdvel_A

ENDIF

!
!3.2 Explicit or next iteration
!------------------------------
!
!3.2.1 Retrieve from first iteration
!-----------------------------------
!

IF (fimp) vdvel_A = vdvel_A_prev

!
!3.2.2 Source/sink terms
!-----------------------
!

IF (iopt_hydro_impl.EQ.1.OR.iopt_grid_nodim.EQ.2) THEN
   WHERE (maskatv)
      vdvel_A(1:ncloc,1:nrloc) = (vdvel_A(1:ncloc,1:nrloc)+sourceatv)/sinkatv
   END WHERE
ELSE
   WHERE (maskatv)
      vdvel_A(1:ncloc,1:nrloc) = (vdvel_A(1:ncloc,1:nrloc)+sourceatv+&
                                & delt2d*vdevint)/sinkatv
   END WHERE
ENDIF

!
!4. Time integration with operator splitting - step B
!----------------------------------------------------
!
!4.1 Source/sink terms
!---------------------
!

IF (iopt_hydro_impl.EQ.1.OR.iopt_grid_nodim.EQ.2) THEN
   WHERE (maskatv)
      vdvel_B(1:ncloc,1:nrloc) = (vdvel_old(1:ncloc,1:nrloc)+sourceatv)/sinkatv
   END WHERE
ELSE
   WHERE (maskatv)
      vdvel_B(1:ncloc,1:nrloc) = (vdvel_old(1:ncloc,1:nrloc)+sourceatv+&
                                & delt2d*vdevint)/sinkatv
   END WHERE
ENDIF

!
!4.2 Y-derivative terms
!----------------------
!
!---initialise
WHERE (maskatv) dexp = 0.0

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/0,0,nhpsi2d,nhpsi2d/)
   CALL exchange_mod(vdvel_B,lbounds,nhexch,iarr_vdvel)
ENDIF

!---advection
IF (iopt_adv_2D.GT.0) THEN
   CALL Yadv_at_V_2d(vdvel_B,advdif2d,nhpsi2d,yadvflag)
   WHERE (maskatv) dexp = dexp - advdif2d
ENDIF

!---diffusion
IF (iopt_hdif_2D.GT.0) THEN
   CALL Ydif_at_V_2d(vdvel_B,advdif2d,ydifflag)
   WHERE (maskatv) dexp = dexp + advdif2d
ENDIF

!---update vdvel_B
WHERE (maskatv) 
   vdvel_B(1:ncloc,1:nrloc) = vdvel_B(1:ncloc,1:nrloc) + delt2d*dexp
END WHERE

!
!4.3 X-derivative terms
!----------------------
!
!---initialise
WHERE (maskatv) dexp = 0.0

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/nhpsi2d,nhpsi2d,0,0/)
   CALL exchange_mod(vdvel_B,lbounds,nhexch,iarr_vdvel)
ENDIF

!---advection
IF (iopt_adv_2D.GT.0) THEN
   CALL Xadv_at_V_2d(vdvel_B,advdif2d,nhpsi2d,xadvflag)
   WHERE (maskatv) dexp = dexp - advdif2d
ENDIF

!---diffusion
IF (iopt_hdif_2D.GT.0) THEN
   CALL Xdif_at_V_2d(vdvel_B,advdif2d,xdifflag)
   WHERE (maskatv) dexp = dexp + advdif2d
ENDIF

!---update vdvel_B
WHERE (maskatv)
   vdvel_B(1:ncloc,1:nrloc) = vdvel_B(1:ncloc,1:nrloc) + delt2d*dexp
END WHERE

!
!5. Time integration with operator splitting - update vdvel
!----------------------------------------------------------
!

WHERE (maskatv) 
   vdvel(1:ncloc,1:nrloc) = 0.5*(vdvel_A(1:ncloc,1:nrloc)+&
                               & vdvel_B(1:ncloc,1:nrloc))
END WHERE

1000 CONTINUE

!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatv,advdif2d,dexp,sinkatv,sourceatv)
IF (time_split) DEALLOCATE (vdvel_A,vdvel_B)

CALL log_timer_out(npcc,itm_trans)


RETURN

END SUBROUTINE transport_at_V_2d

!========================================================================

SUBROUTINE transport_at_W(psiw,source,sink,vdifcoefatw,&
                        & ibcsur,ibcbot,bcsur,bcbot,ivarid)
!************************************************************************
!
! *transport_at_W* Solve advection-diffusion equation for a quantity at W-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Transport_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - dissipation_equation, kl_equation, tke_equation
!
! External calls - Xadv_at_W, Xcorr_at_W, Xdif_at_W, Yadv_at_W, Ycorr_at_W,
!                  Ydif_at_W, Zadv_at_W, Zcorr_at_W, Zdif_at_W
!
! Module calls - error_alloc, exchange_mod, num_halo, tridiag_vert,
!                Uarr_at_UW, Varr_at_VW, Warr_at_C
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
USE physpars
USE switches
USE syspars
USE timepars
USE array_interp, ONLY: Uarr_at_UW, Varr_at_VW, Warr_at_C
USE error_routines, ONLY: error_alloc
USE nla_library, ONLY: tridiag_vert
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: num_halo

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: ibcsur, ibcbot, ivarid
REAL, INTENT(INOUT),DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1)&
                         & :: psiw
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,2:nz) :: source, sink
REAL, INTENT(IN), DIMENSION(ncloc,nrloc,nz+1) :: vdifcoefatw
REAL, INTENT(IN), DIMENSION(ncloc,nrloc) :: bcsur, bcbot

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*psiw*       REAL    W-node quantity to be updated                    [psiw]
!*source*     REAL    Source terms in transport equation             [psiw/s]
!*sink*       REAL    Sink terms in transport equation                  [1/s]
!*vdifcoefatw*REAL    Vertical diffusion coefficient at W-nodes       [m^2/s]
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
!*ivarid*     INTEGER Variable id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: time_split
INTEGER :: k, klo, kup, nhpsi, npcc
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhdims, nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: hdifcoefatuw, hdifcoefatvw, &
                                 & psiw_A, psiw_B, uadvatuw, ucorr, vadvatvw, &
                                 & vcorr, vdifcoefatc, wadvatc, wcorr
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: tridcfw

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*time_split*  LOGICAL .TRUE. for directional splitting in time
!*maskatc*     LOGICAL Mask array at C-nodes
!*maskatu*     LOGICAL Mask array at U-nodes
!*maskatv*     LOGICAL Mask array at V-nodes
!*tridcfw*     REAL    Tridiagonal matrix for implicit vertical solution
!*psiw_A*      REAL    psiw during first integration step in ADI-scheme [psiw]
!*psiw_B*      REAL    psiw during second integration step              [psiw]
!*hdifcoefatuw*REAL    Horizontal diffusion coefficient at UW-nodes    [m^2/s]
!*hdifcoefatvw*REAL    Horizontal diffusion coefficient at VW-nodes    [m^2/s]
!*uadvatuw*    REAL    Advective velocity in X-direction at UW-nodes     [m/s]
!*vadvatvw*    REAL    Advective velocity in Y-direction at VW-nodes     [m/s]
!*vdifcoefatc* REAL    Vertical diffusion coefficient at C-nodes       [m^2/s]
!*wadvatc*     REAL    Advective velocity in Z-direction at C-nodes      [m/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'transport_at_W'
CALL log_timer_in(npcc,ivarid=ivarid)

!
!1. Initialise
!-------------
!
!1.1 Array shape and halo size
!--------------------------
!

nhpsi = MAX(num_halo(iopt_adv_turb),iopt_hdif_turb)
nhdims = nhalo
klo = MERGE(2,3,ibcbot.LT.4)
kup = MERGE(nz,nz-1,ibcsur.LT.4)

!
!1.2 Time integration
!--------------------
!

time_split = iopt_adv_turb.EQ.3

!
!1.3 Allocate arrays
!-------------------
!

ALLOCATE (tridcfw(ncloc,nrloc,nz+1,4),STAT=errstat)
CALL error_alloc('tridcfw',4,(/ncloc,nrloc,nz+1,4/),kndrtype)
ALLOCATE (vdifcoefatc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('vdifcoefatc',3,(/ncloc,nrloc,nz/),kndrtype)

IF (iopt_adv_turb.GT.0) THEN
   ALLOCATE (uadvatuw(2-nhalo:ncloc+nhalo,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('uadvatuw',3,(/ncloc+2*nhalo-1,nrloc,nz-1/),kndrtype)
   ALLOCATE (vadvatvw(ncloc,2-nhalo:nrloc+nhalo,2:nz),STAT=errstat)
   CALL error_alloc('vadvatvw',3,(/ncloc,nrloc+2*nhalo-1,nz-1/),kndrtype)
   ALLOCATE (wadvatc(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('wadvatc',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

IF (iopt_hdif_turb.GT.0) THEN
   ALLOCATE (hdifcoefatuw(ncloc+1,nrloc,2:nz),STAT=errstat)
   CALL error_alloc('hdifcoefatuw',3,(/ncloc+1,nrloc,nz-1/),kndrtype)
   ALLOCATE (hdifcoefatvw(ncloc,nrloc+1,2:nz),STAT=errstat)
   CALL error_alloc('hdifcoefatvw',3,(/ncloc,nrloc+1,nz-1/),kndrtype)
ENDIF

IF (time_split) THEN
   ALLOCATE (psiw_A(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1),STAT=errstat)
   CALL error_alloc('psic_A',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz+1/),kndrtype)
   ALLOCATE (psiw_B(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1),STAT=errstat)
   CALL error_alloc('psic_B',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz+1/),kndrtype)
ENDIF

IF (time_split) THEN
   ALLOCATE (ucorr(ncloc,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('ucorr',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (vcorr(ncloc,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('vcorr',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
   ALLOCATE (wcorr(ncloc,nrloc,klo:kup),STAT=errstat)
   CALL error_alloc('wcorr',3,(/ncloc,nrloc,kup-klo+1/),kndrtype)
ENDIF

!
!1.4 Dirichlet conditions
!------------------------
!
!1.4.1 Surface
!-------------
!

IF (ibcsur.EQ.3) THEN
   WHERE (maskatc_int)
      tridcfw(:,:,nz+1,4) = bcsur
   END WHERE
ELSEIF (ibcsur.EQ.4) THEN
   WHERE (maskatc_int)
      tridcfw(:,:,nz+1,4) = psiw(1:ncloc,1:nrloc,nz+1)
      tridcfw(:,:,nz,1) = 0.0
      tridcfw(:,:,nz,2) = 1.0
      tridcfw(:,:,nz,3) = 0.0
      tridcfw(:,:,nz,4) = bcsur
   END WHERE
ELSE
   WHERE (maskatc_int)
      tridcfw(:,:,nz+1,4) = psiw(1:ncloc,1:nrloc,nz+1)
   END WHERE
ENDIF

WHERE (maskatc_int)
   tridcfw(:,:,nz+1,1) = 0.0; tridcfw(:,:,nz+1,2) = 1.0
   tridcfw(:,:,nz+1,3) = 0.0
END WHERE

!
!1.4.2 Bottom
!------------
!

IF (ibcbot.EQ.3) THEN
   WHERE (maskatc_int)
      tridcfw(:,:,1,4) = bcbot
   END WHERE
ELSEIF (ibcbot.EQ.4) THEN
   WHERE (maskatc_int)
      tridcfw(:,:,1,4) = psiw(1:ncloc,1:nrloc,1)
      tridcfw(:,:,2,1) = 0.0
      tridcfw(:,:,2,2) = 1.0
      tridcfw(:,:,2,3) = 0.0
      tridcfw(:,:,2,4) = bcbot
   END WHERE
ELSE
   WHERE (maskatc_int)
      tridcfw(:,:,1,4) = psiw(1:ncloc,1:nrloc,1)
   END WHERE
ENDIF

WHERE (maskatc_int)
   tridcfw(:,:,1,1) = 0.0; tridcfw(:,:,1,2) = 1.0; tridcfw(:,:,1,3) = 0.0
END WHERE

!
!1.5 Local arrays
!----------------
!
!---temporary values of psiw
IF (time_split) THEN
   psiw_A = psiw; psiw_B = psiw
ENDIF

!
!1.6 Source/sink terms
!---------------------
!

k_160:DO k=klo,kup
   WHERE (maskatc_int)
      source(:,:,k) = delt3d*source(:,:,k)
      sink(:,:,k) = delt3d*sink(:,:,k) 
   END WHERE
ENDDO k_160

!
!1.7 Advective velocities
!------------------------
!

IF (iopt_adv_turb.GT.0) THEN
   CALL Uarr_at_UW(uvel(2-nhalo:ncloc+nhalo,1:nrloc,:),uadvatuw,3,3,&
                & (/2-nhalo,1,1/),(/ncloc+nhalo,nrloc,nz/),1,iarr_uvel,&
                & .TRUE.)
   CALL Varr_at_VW(vvel(1:ncloc,2-nhalo:nrloc+nhalo,:),vadvatvw,3,3,&
                & (/1,2-nhalo,1/),(/ncloc,nrloc+nhalo,nz/),1,iarr_vvel,&
                & .TRUE.)
   CALL Warr_at_C(wfvel(1:ncloc,1:nrloc,:),wadvatc,(/1,1,1/),&
               & (/ncloc,nrloc,nz+1/),1,iarr_wvel,.TRUE.)
ENDIF

!
!1.8 Horizontal diffusion coefficients
!-------------------------------------
!

IF ((iopt_hdif_turb.GT.0).AND.(iopt_hdif_coef.EQ.1)) THEN
   k_190: DO k=2,nz
      WHERE (nodeatu(1:ncloc+1,1:nrloc,2:nz).GT.1) hdifcoefatuw = hdifscal_cst
      WHERE (nodeatv(1:ncloc,1:nrloc+1,2:nz).GT.1) hdifcoefatvw = hdifscal_cst
   ENDDO k_190
ELSEIF ((iopt_hdif_turb.GT.0).AND.(iopt_hdif_coef.EQ.2)) THEN
   CALL Uarr_at_UW(hdifcoef3datu,hdifcoefatuw,3,3,(/1,1,1/),&
                & (/ncloc+1,nrloc,nz/),1,iarr_hdifcoef3datu,.TRUE.)
   CALL Varr_at_VW(hdifcoef3datv,hdifcoefatvw,3,3,(/1,1,1/),&
                & (/ncloc,nrloc+1,nz/),1,iarr_hdifcoef3datv,.TRUE.)
ENDIF

!
!1.9 Vertical diffusion coefficient
!----------------------------------
!

CALL Warr_at_C(vdifcoefatw,vdifcoefatc,(/1,1,1/),(/ncloc,nrloc,nz+1/),1,0,&
            & .TRUE.) 

!
!2. Corrector terms
!------------------
!

IF (time_split) THEN
   CALL Xcorr_at_W(ucorr,uadvatuw,klo,kup)
   CALL Ycorr_at_W(vcorr,vadvatvw,klo,kup)
   CALL Zcorr_at_W(wcorr,wadvatc,klo,kup)
ENDIF

!
!3. Time integration without operator splitting
!----------------------------------------------
!

IF (.NOT.time_split) THEN

!
!3.1 Initialise
!--------------
!

   k_310: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfw(:,:,k,1) = 0.0
         tridcfw(:,:,k,3) = 0.0
      END WHERE
   ENDDO k_310

!
!3.2 Time derivative
!-------------------
!

   k_320: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfw(:,:,k,2) = 1.0
         tridcfw(:,:,k,4) = psiw(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_320

!
!3.3 Source and sink terms
!-------------------------
!

   k_330: DO k=klo,kup
      WHERE (maskatc_int)
         tridcfw(:,:,k,2) = tridcfw(:,:,k,2) + sink(:,:,k)
         tridcfw(:,:,k,4) = tridcfw(:,:,k,4) + source(:,:,k)
      END WHERE
   ENDDO k_330

!
!3.4 Advection
!-------------
!

   IF (iopt_adv_turb.GT.0) THEN
      CALL Xadv_at_W(psiw,tridcfw,uadvatuw,nhpsi,klo,kup)
      CALL Yadv_at_W(psiw,tridcfw,vadvatvw,nhpsi,klo,kup)
      CALL Zadv_at_W(psiw,tridcfw,wadvatc,klo,kup)
   ENDIF

!
!3.5 Diffusion
!-------------
!
!  ---horizontal
   IF (iopt_hdif_turb.GT.0) THEN
      CALL Xdif_at_W(psiw,tridcfw,hdifcoefatuw,klo,kup)
      CALL Ydif_at_W(psiw,tridcfw,hdifcoefatvw,klo,kup)
   ENDIF
!  ---vertical
   CALL Zdif_at_W(psiw,tridcfw,vdifcoefatc,ibcsur,ibcbot,bcsur,bcbot,klo,kup)

!
!3.6 Update psiw
!---------------
!

   CALL tridiag_vert(tridcfw,psiw,nhdims,nz+1,1,'C  ')

   GOTO 1000

ENDIF

!
!4. Time integration with operator splitting - step A
!----------------------------------------------------
!
!4.1 X-derivative terms
!----------------------
!
!---time derivative, corrector term
k_411: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfw(:,:,k,4) = psiw(1:ncloc,1:nrloc,k)*(1.0+ucorr(:,:,k))
   END WHERE
ENDDO k_411

!---horizontal advection
IF (iopt_adv_turb.GT.0) THEN
   CALL Xadv_at_W(psiw,tridcfw,uadvatuw,nhpsi,klo,kup)
ENDIF

!---horizontal diffusion
IF (iopt_hdif_turb.GT.0) THEN
   CALL Xdif_at_W(psiw,tridcfw,hdifcoefatuw,klo,kup)
ENDIF

!---update psiw_A
k_412: DO k=klo,kup
   WHERE (maskatc_int)
      psiw_A(1:ncloc,1:nrloc,k) = tridcfw(:,:,k,4)
   END WHERE
ENDDO k_412

!
!4.2 Y-derivative terms
!----------------------
!
!---time derivative, corrector term
k_421: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfw(:,:,k,4) = psiw_A(1:ncloc,1:nrloc,k)*(1.0+vcorr(:,:,k))
   END WHERE
ENDDO k_421

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   lbounds = (/1-nhalo,1-nhalo,1/); nhexch = (/0,0,nhpsi,nhpsi/)
   CALL exchange_mod(psiw_A,lbounds,nhexch,0)
ENDIF

!---horizontal advection
IF (iopt_adv_turb.GT.0) THEN
   CALL Yadv_at_W(psiw_A,tridcfw,vadvatvw,nhpsi,klo,kup)
ENDIF

!---horizontal diffusion
IF (iopt_hdif_turb.GT.0) THEN
   CALL Ydif_at_W(psiw_A,tridcfw,hdifcoefatvw,klo,kup)
ENDIF

!---update psiw_A
k_422: DO k=klo,kup
   WHERE (maskatc_int)
      psiw_A(1:ncloc,1:nrloc,k) = tridcfw(:,:,k,4)
   END WHERE
ENDDO k_422

!
!4.3 Z-derivative, source/sink terms
!-----------------------------------
!
!---time derivative, source/sink and corrector terms
k_431: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfw(:,:,k,1) = 0.0
      tridcfw(:,:,k,2) = 1.0 + sink(:,:,k) - theta_vadv*wcorr(:,:,k)
      tridcfw(:,:,k,3) = 0.0
      tridcfw(:,:,k,4) = psiw_A(1:ncloc,1:nrloc,k)*&
                       & (1.0+(1.0-theta_vadv)*wcorr(:,:,k)) + source(:,:,k)
   END WHERE
ENDDO k_431

!---vertical advection
IF (iopt_adv_turb.GT.0) THEN
   CALL Zadv_at_W(psiw_A,tridcfw,wadvatc,klo,kup)
ENDIF

!---vertical diffusion
CALL Zdif_at_W(psiw_A,tridcfw,vdifcoefatc,ibcsur,ibcbot,bcsur,bcbot,klo,kup)

!---update psiw_A
CALL tridiag_vert(tridcfw,psiw_A,nhdims,nz+1,1,'C  ')

!
!5. Time integration with operator splitting - step B
!----------------------------------------------------
!
!5.1 Z-derivative, source/sink terms
!-----------------------------------
!
!---time derivative, source/sink and corrector term
k_511: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfw(:,:,k,1) = 0.0
      tridcfw(:,:,k,2) = 1.0 + sink(:,:,k) - theta_vadv*wcorr(:,:,k)
      tridcfw(:,:,k,3) = 0.0
      tridcfw(:,:,k,4) = psiw(1:ncloc,1:nrloc,k)*&
                       & (1.0+(1.0-theta_vadv)*wcorr(:,:,k)) + source(:,:,k)
   END WHERE
ENDDO k_511

!---vertical advection
IF (iopt_adv_turb.GT.0) THEN
   CALL Zadv_at_W(psiw,tridcfw,wadvatc,klo,kup)
ENDIF

!---vertical diffusion
CALL Zdif_at_W(psiw,tridcfw,vdifcoefatc,ibcsur,ibcbot,bcsur,bcbot,klo,kup)

!---update psiw_B
CALL tridiag_vert(tridcfw,psiw_B,nhdims,nz+1,1,'C  ')

!
!5.2 Y-derivative terms
!----------------------
!
!---time derivative, corrector term
k_521: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfw(:,:,k,4) = psiw_B(1:ncloc,1:nrloc,k)*(1.0+vcorr(:,:,k))
   END WHERE
ENDDO k_521

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/0,0,nhpsi,nhpsi/)
   CALL exchange_mod(psiw_B,lbounds,nhexch,0)
ENDIF

!---horizontal advection
IF (iopt_adv_turb.GT.0) THEN
   CALL Yadv_at_W(psiw_B,tridcfw,vadvatvw,nhpsi,klo,kup)
ENDIF

!---horizontal diffusion
IF (iopt_hdif_turb.GT.0) THEN
   CALL Ydif_at_W(psiw_B,tridcfw,hdifcoefatvw,klo,kup)
ENDIF

!---update psiw_B
k_522: DO k=klo,kup
   WHERE (maskatc_int)
      psiw_B(1:ncloc,1:nrloc,k) = tridcfw(:,:,k,4)
   END WHERE
ENDDO k_522

!
!5.3 X-derivative terms
!----------------------
!
!---time derivative, corrector term
k_531: DO k=klo,kup
   WHERE (maskatc_int)
      tridcfw(:,:,k,4) = psiw_B(1:ncloc,1:nrloc,k)*(1.0+ucorr(:,:,k))
   END WHERE
ENDDO k_531

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   nhexch = (/nhpsi,nhpsi,0,0/)
   CALL exchange_mod(psiw_B,lbounds,nhexch,0)
ENDIF

!---horizontal advection
IF (iopt_adv_turb.GT.0) THEN
   CALL Xadv_at_W(psiw_B,tridcfw,uadvatuw,nhpsi,klo,kup)
ENDIF

!---horizontal diffusion
IF (iopt_hdif_turb.GT.0) THEN
   CALL Xdif_at_W(psiw_B,tridcfw,hdifcoefatuw,klo,kup)
ENDIF

!---update psiw_B
k_532: DO k=klo,kup
   WHERE (maskatc_int)
      psiw_B(1:ncloc,1:nrloc,k) = tridcfw(:,:,k,4)
   END WHERE
ENDDO k_532

!
!6. Time integration with operator splitting - update psiw
!---------------------------------------------------------
!

k_610: DO k=1,nz+1
   WHERE (maskatc_int)
      psiw(1:ncloc,1:nrloc,k) = 0.5*(psiw_A(1:ncloc,1:nrloc,k)+&
                                   & psiw_B(1:ncloc,1:nrloc,k))
   END WHERE
ENDDO k_610

1000 CONTINUE

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE (tridcfw)
IF (iopt_adv_turb.GT.0) DEALLOCATE (uadvatuw,vadvatvw,wadvatc)
IF (iopt_hdif_turb.GT.0) DEALLOCATE (hdifcoefatuw,hdifcoefatvw)
IF (iopt_vdif_coef.GT.0) DEALLOCATE (vdifcoefatc)
IF (time_split) DEALLOCATE (psiw_A,psiw_B)
IF (time_split) DEALLOCATE (ucorr,vcorr,wcorr)
   
CALL log_timer_out(npcc,itm_trans)


RETURN

END SUBROUTINE transport_at_W
