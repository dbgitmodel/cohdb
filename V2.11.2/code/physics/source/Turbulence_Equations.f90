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
! *Turbulence_Equations* Turbulence closure equations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description -
!
! Routines -  dissipation_equation, dissip_lim_conds, dissip_to_zlmix,
!             eddy_coefs_alg, eddy_coefs_tc, init_turbulence, kl_equation,
!             mixing_length, tke_equation, tke_equilibrium, zlmix_lim_conds,
!             zlmix_to_dissip
!
!************************************************************************
!


!========================================================================

SUBROUTINE dissipation_equation
!************************************************************************
!
! *dissipation_equation* Solve dissipation equation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - vertical_diff_coefs
!
! External calls - transport_at_W
!
! Module calls - error_alloc, exchange_mod
!
!************************************************************************
!
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE turbpars
USE turbulence
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: ibcbot, ibcsur, k
REAL :: vdifeps_cst
INTEGER, DIMENSION(4) :: nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, bcbot, bcsur, buoprod
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: sink, source, vdifcoefdissip

!
! Name           Type    Purpose
!------------------------------------------------------------------------------
!*buoprod*       REAL    Buoyancy production term                        [W/kg]
!*sink*          REAL    Sink terms in dissipation equation               [1/s]
!*source*        REAL    Source terms in dissipation equation          [W/kg/s]
!*vdifcoefdissip*REAL    Vertical diffusion coefficient in transport equation
!                                                                       [m^2/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'dissipation_equation'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (bcbot(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcbot',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (bcsur(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcsur',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (buoprod(ncloc,nrloc),STAT=errstat)
CALL error_alloc('buoprod',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sink(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('sink',3,(/ncloc,nrloc,nz-1/),kndrtype)
ALLOCATE (source(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('source',3,(/ncloc,nrloc,nz-1/),kndrtype)
ALLOCATE (vdifcoefdissip(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('vdifcoefdissip',3,(/ncloc,nrloc,nz+1/),kndrtype)

!
!2. Initialise
!-------------
!
!---diffusion coefficients
k_210: DO k=1,nz+1
   WHERE (.NOT.maskatc_int)
      vdifcoefdissip(:,:,k) = 0.0
   END WHERE
ENDDO k_210

IF (nt.EQ.1) THEN

!  ---surface value
   IF (zrough_sur.GT.0) THEN
      WHERE (maskatc_int)
         dissip(1:ncloc,1:nrloc,nz+1) = (sstresatc**1.5)/(ckar*zrough_sur)
      END WHERE
   ELSE
      WHERE (maskatc_int)
         dissip(1:ncloc,1:nrloc,nz+1) = 0.0
      END WHERE
   ENDIF

!  ---bottom value
   IF (zrough_bot.GT.0) THEN
      WHERE (maskatc_int)
         dissip(1:ncloc,1:nrloc,1) = (bstresatc(1:ncloc,1:nrloc)**1.5)/&
                                   & (ckar*zrough_bot)
      END WHERE
   ELSE
      WHERE (maskatc_int)
         dissip(1:ncloc,1:nrloc,1) = 0.0
      END WHERE
   ENDIF

ENDIF

!
!3. Diffusion coefficient
!------------------------
!

vdifeps_cst = vdifmom_cst*(sigma_eps-1.0)
k_310: DO k=1,nz+1
   WHERE (maskatc_int)
      vdifcoefdissip(:,:,k) = (vdifcoefmom(1:ncloc,1:nrloc,k)+vdifeps_cst)&
                            & /sigma_eps
   END WHERE
ENDDO k_310

!
!4. Boundary conditions
!----------------------
!
!4.1 Surface
!-----------
!
!---Neumann
IF (iopt_turb_dis_sbc.EQ.1) THEN
   ibcsur = 2
   WHERE (maskatc_int)
      bcsur = (0.5*eps0/ckar)*&
            & (vdifcoefdissip(:,:,nz)+vdifcoefdissip(:,:,nz+1))&
            & /(zrough_sur+0.5*delzatc(1:ncloc,1:nrloc,nz))**2
   END WHERE
!---Dirichlet
ELSEIF (iopt_turb_dis_sbc.EQ.2) THEN
   ibcsur = 4
   WHERE (maskatc_int)
      bcsur = (eps0/ckar)*(tke(1:ncloc,1:nrloc,nz)**1.5)/&
            & (zrough_sur+delzatc(1:ncloc,1:nrloc,nz))
   END WHERE
ENDIF

!
!4.2 Bottom
!----------
!
!---Neumann
IF (iopt_turb_dis_bbc.EQ.1) THEN
   ibcbot = 2
   WHERE (maskatc_int)
      bcbot = -(0.5*eps0/ckar)*&
             & (vdifcoefdissip(:,:,1)+vdifcoefdissip(:,:,2))&
             & /(zrough_bot+0.5*delzatc(1:ncloc,1:nrloc,1))**2
   END WHERE
!---Dirichlet
ELSEIF (iopt_turb_dis_bbc.EQ.2) THEN
   ibcbot = 4
   WHERE (maskatc_int)
      bcbot = (eps0/ckar)*(tke(1:ncloc,1:nrloc,2)**1.5)/&
            & (zrough_bot+delzatc(1:ncloc,1:nrloc,1))
   END WHERE
ENDIF

!
!5. Source/sink terms
!--------------------
!

WHERE (.NOT.maskatc_int) buoprod = 0.0

k_510: DO k=2,nz
   WHERE (maskatc_int)
      array2dc = dissip(1:ncloc,1:nrloc,k)/tke_old(1:ncloc,1:nrloc,k)
      source(:,:,k) = c1_eps*array2dc*(vdifcoefmom(1:ncloc,1:nrloc,k)-&
                    & vdifmom_cst)*shearfreq2(:,:,k)
      sink(:,:,k) = c2_eps*array2dc
      buoprod = c1_eps*(vdifcoefscal(:,:,k)-vdifscal_cst)*buofreq2(:,:,k)
   END WHERE
   WHERE (buoprod.LT.0.0)
      source(:,:,k) = source(:,:,k) - c32_eps*array2dc*buoprod
   END WHERE
   WHERE (buoprod.GT.0.0)
      sink(:,:,k) = sink(:,:,k) + c31_eps*buoprod/tke_old(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_510

!
!6. Solve transport equation
!---------------------------
!

CALL transport_at_W(dissip,source,sink,vdifcoefdissip,ibcsur,ibcbot,bcsur,&
                  & bcbot,iarr_dissip)

!
!7. Surface/bottom values
!------------------------
!
!---surface value
IF (zrough_sur.GT.0) THEN
   WHERE (maskatc_int)
      dissip(1:ncloc,1:nrloc,nz+1) = (sstresatc**1.5)/(ckar*zrough_sur)
   END WHERE
ELSE
   WHERE (maskatc_int)
      dissip(1:ncloc,1:nrloc,nz+1) = 0.0
   END WHERE
ENDIF

!---bottom value
IF (zrough_bot.GT.0) THEN
   WHERE (maskatc_int)
      dissip(1:ncloc,1:nrloc,1) = (bstresatc(1:ncloc,1:nrloc)**1.5)/&
                                & (ckar*zrough_bot)
   END WHERE
ELSE
   WHERE (maskatc_int)
      dissip(1:ncloc,1:nrloc,1) = 0.0
   END WHERE
ENDIF

!
!8. Apply lower limit (numerical)
!--------------------------------
!

k_810: DO k=2,nz
   WHERE (maskatc_int)
      dissip(1:ncloc,1:nrloc,k) = MAX(dissipmin,dissip(1:ncloc,1:nrloc,k))
   END WHERE
ENDDO k_810

!
!9. Exchange halo sections
!-------------------------
!

IF ((iopt_MPI.EQ.1).AND.(iopt_adv_turb.GT.0.OR.iopt_hdif_turb.GT.0)) THEN
   nhexch = nhturb
   CALL exchange_mod(dissip,(/1-nhalo,1-nhalo,1/),nhexch,iarr_dissip,&
                   & corners=.FALSE.)
ENDIF

!
!10. Deallocate arrays
!---------------------
!

DEALLOCATE (array2dc,bcbot,bcsur,buoprod,sink,source,vdifcoefdissip)

CALL log_timer_out()


RETURN

END SUBROUTINE dissipation_equation

!========================================================================

SUBROUTINE dissip_lim_conds
!************************************************************************
!
! *dissip_lim_conds* Apply limiting conditions for dissipation rate
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.0
!
! Description - 
!
! Calling program - vertical_diff_coefs
!
! Module calls -
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE turbpars
USE turbulence
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
INTEGER :: i, j, k
REAL, PARAMETER :: epstol = 1.0E-20


procname(pglev+1) = 'dissip_lim_conds'
CALL log_timer_in()

!---limiting conditions
k_110: DO k=2,nz
   j_111: DO j=1,nrloc
   i_111: DO i=1,ncloc
      IF (maskatc_int(i,j)) THEN
         tke(i,j,k) = MAX(tkelim,tke(i,j,k))
         IF (buofreq2(i,j,k).GT.epstol) THEN
            dissip(i,j,k) = MAX(tke(i,j,k)*&
                          & SQRT(buofreq2(i,j,k)/alphaN_max),dissip(i,i,k))
         ENDIF
      ENDIF
   ENDDO i_111
   ENDDO j_111
ENDDO k_110

CALL log_timer_out()


RETURN

END SUBROUTINE dissip_lim_conds

!========================================================================

SUBROUTINE dissip_to_zlmix
!************************************************************************
!
! *dissip_to_zlmix* Mixing length from dissipation rate
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.0
!
! Description - 
!
! Calling program - vertical_diff_coefs
!
! Module calls -
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE physpars
USE timepars
USE turbpars
USE turbulence
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
INTEGER :: i, j, k


procname(pglev+1) = 'dissip_to_zlmix'
CALL log_timer_in()

!---surface/bottom values
IF (nt.EQ.1) THEN
   WHERE (maskatc_int)
      zlmix(1:ncloc,1:nrloc,nz+1) = ckar*zrough_sur
      zlmix(1:ncloc,1:nrloc,1) = ckar*zrough_bot
   END WHERE
ENDIF

!---interior points
k_110: DO k=2,nz
   j_111: DO j=1,nrloc
   i_111: DO i=1,ncloc
      IF (maskatc_int(i,j)) THEN
         zlmix(i,j,k) = MAX(zlmixmin,eps0*tke(i,j,k)**1.5/dissip(i,j,k))
      ENDIF
   ENDDO i_111
   ENDDO j_111
ENDDO k_110

CALL log_timer_out()


RETURN

END SUBROUTINE dissip_to_zlmix

!========================================================================

SUBROUTINE eddy_coefs_alg
!************************************************************************
!
! *eddy_coefs_alg* Eddy coefficients from algebraic relations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - vertical_diff_coefs
!
! External calls -
!
! Module calls - error_alloc, exchange_mod, mom_strat_MA, richardson_number,
!                scal_strat_MA, Uarr_at_C, Varr_at_C
!
!************************************************************************
!
USE currents
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE turbpars
USE turbulence
USE array_interp, ONLY: Uarr_at_C, Varr_at_C
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY : exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE turbulence_routines, ONLY: mom_strat_MA, richardson_number, scal_strat_MA

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, k
REAL :: denom0, epsmin
INTEGER, DIMENSION(4) :: nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: phiprof
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: alpha, array2dc1, array2dc2, &
                                         & delta_b, denom, ricnum, vedwind

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*phiprof*    REAL    Normalised vertical profile for flow-dependent relations
!*alpha*      REAL    Current-dependence factor in flow-dependent relations
!*delta_b*    REAL    Bed layer thickness                                   [m]
!*ricnum*     REAL    Richardson number
!*vedwind*    REAL    Diffusion coefficient due to surface wind         [m^2/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'eddy_coefs_alg'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (alpha(ncloc,nrloc),STAT=errstat)
CALL error_alloc('alpha',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dc1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dc2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (delta_b(ncloc,nrloc),STAT=errstat)
CALL error_alloc('delta_b',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (denom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('denom',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (ricnum(ncloc,nrloc),STAT=errstat)
CALL error_alloc('ricnum',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (vedwind(ncloc,nrloc),STAT=errstat)
CALL error_alloc('vedwind',2,(/ncloc,nrloc/),kndrtype)

!
!2. Initialise normalised vertical profile on first call
!-------------------------------------------------------
!

IF (nt.EQ.0) THEN

   IF ((iopt_turb_alg.EQ.3).OR.(iopt_turb_alg.EQ.4).OR.&
      &(iopt_turb_alg.EQ.5)) THEN

!     ---allocate vertical profile array
      ALLOCATE (phiprof(ncloc,nrloc,nz+1),STAT=errstat)
      CALL error_alloc('phiprof',3,(/ncloc,nrloc,nz+1/),kndrtype)
      phiprof = 0.0

!     ---construct vertical profile
      denom0 = 1.0 + 0.5*delta1_ad*(r1_ad-1.0) + 0.5*delta2_ad*(r2_ad-1.0)
      i_210: DO i=1,ncloc
      j_210: DO j=1,nrloc
         IF (maskatc_int(i,j)) THEN
            k_211: DO k=1,nz+1
               IF (gscoordatw(i,j,k).LT.delta1_ad) THEN
                  phiprof(i,j,k) = (1.0-r1_ad)*gscoordatw(i,j,k)/delta1_ad+&
                                  & r1_ad
               ELSEIF (gscoordatw(i,j,k).GT.(1.0-delta2_ad)) THEN
                  phiprof(i,j,k) = r2_ad-(r2_ad-1.0)*&
                                & (1.0-gscoordatw(i,j,k))/delta2_ad
               ELSE
                  phiprof(i,j,k) = 1.0
               ENDIF
               phiprof(i,j,k) = phiprof(i,j,k)/denom0
            ENDDO k_211
         ENDIF
      ENDDO j_210
      ENDDO i_210

   ENDIF

!  ---dissipation rate at surface/bottom
   IF (zrough_sur.EQ.0) THEN
      WHERE (maskatc_int) dissip(1:ncloc,1:nrloc,nz+1) = 0.0
   ENDIF
   IF (zrough_bot.EQ.0) THEN
      WHERE (maskatc_int) dissip(1:ncloc,1:nrloc,1) = 0.0
   ENDIF

ENDIF

!
!3. Pacanowski-Philander
!-----------------------
!

IF (iopt_turb_alg.EQ.1) THEN
   epsmin = vmax_pp**(-1.0/expmom_pp)
   k_310: DO k=1,nz+1
      WHERE (maskatc_int)
         ricnum = richardson_number(k,maskatc_int)
         denom = MAX(1.0+alpha_pp*ricnum,epsmin)
         vdifcoefmom(1:ncloc,1:nrloc,k) = v0dif_pp/(denom**expmom_pp) +&
                                        & vbmom_pp
         vdifcoefscal(:,:,k) = vdifcoefmom(1:ncloc,1:nrloc,k)/denom + vbscal_pp
      END WHERE
   ENDDO k_310

!
!4. Munk-Anderson
!----------------
!

ELSEIF (iopt_turb_alg.EQ.2) THEN
   k_410: DO k=1,nz+1
      WHERE (maskatc_int)
         ricnum = richardson_number(k,maskatc_int)
         vdifcoefmom(1:ncloc,1:nrloc,k) = v0dif_ma*&
                                        & mom_strat_MA(ricnum,maskatc_int)
         vdifcoefscal(:,:,k) = v0dif_ma*scal_strat_MA(ricnum,maskatc_int)
      END WHERE
   ENDDO k_410

!
!5. Flow-dependent relations
!---------------------------
!
!5.1 Proportional to current magnitude
!-------------------------------------
!

ELSEIF (iopt_turb_alg.EQ.3) THEN
   CALL Uarr_at_C(udvel(1:ncloc+1,1:nrloc),array2dc1,1,1,(/1,1,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_udvel,.TRUE.)
   CALL Varr_at_C(vdvel(1:ncloc,1:nrloc+1),array2dc2,1,1,(/1,1,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vdvel,.TRUE.)
   WHERE (maskatc_int)
      alpha = k1_ad*SQRT(array2dc1**2+array2dc2**2)
      vedwind = lambda_ad*SQRT(sstresatc)
   END WHERE
   k_510: DO k=1,nz+1
      WHERE (maskatc_int)
         array2dc1 = alpha*phiprof(:,:,k) + vedwind
         ricnum = richardson_number(k,maskatc_int)
         vdifcoefmom(1:ncloc,1:nrloc,k) = array2dc1*&
                               & mom_strat_MA(ricnum,maskatc_int) + vdifmom_cst
         vdifcoefscal(:,:,k) = array2dc1*scal_strat_MA(ricnum,maskatc_int) + &
                               & vdifscal_cst
      END WHERE
   ENDDO k_510

!
!5.2 Proportional to squared current magnitude
!---------------------------------------------
!

ELSEIF (iopt_turb_alg.EQ.4) THEN
   CALL Uarr_at_C(umvel(1:ncloc+1,1:nrloc),array2dc1,1,1,(/1,1,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_umvel,.TRUE.)
   CALL Varr_at_C(vmvel(1:ncloc,1:nrloc+1),array2dc2,1,1,(/1,1,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vmvel,.TRUE.)
   WHERE (maskatc_int)
      alpha = (k2_ad/omega1_ad)*(array2dc1**2+array2dc2**2)
      vedwind = lambda_ad*SQRT(sstresatc)
   END WHERE
   k_520: DO k=1,nz+1
      WHERE (maskatc_int)
         array2dc1 = alpha*phiprof(:,:,k) + vedwind
         ricnum = richardson_number(k,maskatc_int)
         vdifcoefmom(1:ncloc,1:nrloc,k) = array2dc1*&
                               & mom_strat_MA(ricnum,maskatc_int) + vdifmom_cst
         vdifcoefscal(:,:,k) = array2dc1*scal_strat_MA(ricnum,maskatc_int) + &
                             & vdifscal_cst
      END WHERE
   ENDDO k_520

!
!5.3 Proportional to current magnitude and fractional layer thickness
!--------------------------------------------------------------------
!

ELSEIF (iopt_turb_alg.EQ.5) THEN
   CALL Uarr_at_C(umvel(1:ncloc+1,1:nrloc),array2dc1,1,1,(/1,1,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_umvel,.TRUE.)
   CALL Varr_at_C(vmvel(1:ncloc,1:nrloc+1),array2dc2,1,1,(/1,1,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vmvel,.TRUE.)
   WHERE (maskatc_int)
      delta_b = MIN((cnu_ad/omega1_ad)*bstresatc(1:ncloc,1:nrloc),&
                   & deptotatc(1:ncloc,1:nrloc))
      alpha = k1_ad*delta_b*SQRT(array2dc1**2+array2dc2**2)
      vedwind = lambda_ad*SQRT(sstresatc)
   END WHERE
   k_530: DO k=1,nz+1
      WHERE (maskatc_int)
         array2dc1 = alpha*phiprof(:,:,k) + vedwind
         ricnum = richardson_number(k,maskatc_int)
         vdifcoefmom(1:ncloc,1:nrloc,k) = array2dc1*&
                               & mom_strat_MA(ricnum,maskatc_int) + vdifmom_cst
         vdifcoefscal(:,:,k) = array2dc1*scal_strat_MA(ricnum,maskatc_int) + &
                             & vdifscal_cst
      END WHERE
   ENDDO k_530

!
!6. Parabolic eddy viscosity
!---------------------------
!

ELSEIF (iopt_turb_alg.EQ.6) THEN
   k_640: DO k=1,nz+1
      WHERE (maskatc_int)
         ricnum = richardson_number(k,maskatc_int)
         array2dc1 = ckar*SQRT(bstresatc(1:ncloc,1:nrloc))*&
                & gscoordatw(1:ncloc,1:nrloc,k)*&
                & (1.0-gscoordatw(1:ncloc,1:nrloc,k))*deptotatc(1:ncloc,1:nrloc)
         vdifcoefmom(1:ncloc,1:nrloc,k) = array2dc1*&
                                        & mom_strat_MA(ricnum,maskatc_int)
         vdifcoefscal(1:ncloc,1:nrloc,k) = array2dc1*&
                                         & scal_strat_MA(ricnum,maskatc_int)
      END WHERE
   ENDDO k_640
ENDIF

!
!7. Dissipation rate
!-------------------
!
!---interior points
k_710: DO k=2,nz
   WHERE (maskatc_int)
      dissip(1:ncloc,1:nrloc,k)  =  vdifcoefmom(1:ncloc,1:nrloc,k)*&
                                  & shearfreq2(1:ncloc,1:nrloc,k)
   END WHERE
   IF (iopt_dens.GT.0.OR.iopt_sed.GT.0) THEN
      WHERE (maskatc_int)
         dissip(1:ncloc,1:nrloc,k) = dissip(1:ncloc,1:nrloc,k) - &
                   & vdifcoefscal(1:ncloc,1:nrloc,k)*buofreq2(1:ncloc,1:nrloc,k)
      END WHERE
   ENDIF
ENDDO k_710

!---surface value
IF (zrough_sur.GT.0) THEN
   WHERE (maskatc_int)
      dissip(1:ncloc,1:nrloc,nz+1) = (sstresatc**1.5)/(ckar*zrough_sur)
   END WHERE
ENDIF

!---bottom value
IF (zrough_bot.GT.0) THEN
   WHERE (maskatc_int)
      dissip(1:ncloc,1:nrloc,1) = (bstresatc(1:ncloc,1:nrloc)**1.5)/&
                                & (ckar*zrough_bot)
   END WHERE
ENDIF

!
!8. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = (/1,0,1,0/)
   CALL exchange_mod(vdifcoefmom,(/0,0,1/),nhexch,iarr_vdifcoefmom)
ENDIF

!
!9. Deallocate arrays
!--------------------
!

DEALLOCATE (alpha,array2dc1,array2dc2,delta_b,denom,ricnum,vedwind)
IF ((nt+ic3d.GT.nstep.OR.cold_start).AND.ALLOCATED(phiprof)) THEN
   DEALLOCATE (phiprof)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE eddy_coefs_alg

!========================================================================

SUBROUTINE eddy_coefs_tc
!************************************************************************
!
! *eddy_coefs_tc* Eddy coefficients from turbulence closure
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.11.2
!
! Description - uses concept of stability functions
!
! Reference -
!
! Calling program - vertical_diff_coefs
!
! Module calls - error_abort, error_alloc, error_lbound_arr, exchange_mod,
!                mom_strat_MA, quadratic_root_arr, richardson_number,
!                scal_strat_MA
!
!************************************************************************
!
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE turbpars
USE turbulence
USE error_routines, ONLY: error_abort, error_alloc, error_lbound_arr
USE math_library, ONLY: quadratic_root_arr
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE turbulence_routines, ONLY: mom_strat_MA, richardson_number, scal_strat_MA

!
!* Local variables
!
INTEGER :: i, j, k, noerrs
INTEGER, DIMENSION(4) :: nhexch
REAL :: xmax
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: alphaM, alphaN, array2dc1, &
                          & array2dc2, fstabmom, fstabscal, fstabtke, ricnum
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: coef
COMPLEX, SAVE, ALLOCATABLE, DIMENSION(:,:) :: root1, root2

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*alphaM*     Shear stability parameter
!*alphaN*     Stratification stability parameter
!*fstabmom*   Stability function for momentum
!*fstabscal*  Stability function for scalars
!*fstbatke*   Stability function for t.k.e.
!*ricnum*     Richardson number
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'eddy_coefs_tc'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (alphaM(ncloc,nrloc),STAT=errstat)
CALL error_alloc('alphaM',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (alphaN(ncloc,nrloc),STAT=errstat)
CALL error_alloc('alphaN',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dc1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dc2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (fstabmom(ncloc,nrloc),STAT=errstat)
CALL error_alloc('fstabmom',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (fstabscal(ncloc,nrloc),STAT=errstat)
CALL error_alloc('fstabscal',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (fstabtke(ncloc,nrloc),STAT=errstat)
CALL error_alloc('fstabtke',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (ricnum(ncloc,nrloc),STAT=errstat)
CALL error_alloc('ricnum',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (root1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('root1',2,(/ncloc,nrloc/),kndcmplx)
ALLOCATE (root2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('root2',2,(/ncloc,nrloc/),kndcmplx)
ALLOCATE (coef(ncloc,nrloc,3),STAT=errstat)
CALL error_alloc('coef',3,(/ncloc,nrloc,nz/),kndrtype)

!
!2. Initialise
!-------------
!
!---coefficient arrays
coef = 0.0

!---limiting condition
xmax = alphaN_max-alphaN_sl*alphaM_max

!
!3. Error checking
!-----------------
!

IF (errchk) THEN
   IF (iopt_turb_param.EQ.1) THEN
      noerrs = COUNT((tke(1:ncloc,1:nrloc,2:nz).GT.0.0).AND.&
                   & (zlmix(1:ncloc,1:nrloc,2:nz).LE.0.0))
      IF (noerrs.GT.0) THEN
         i_310: DO i=1,ncloc
         j_310: DO j=1,nrloc
         k_310: DO k=1,nz
            IF ((tke(i,j,k).GT.0.0).AND.(zlmix(i,j,k).LE.0.0)) THEN
               CALL error_lbound_arr(zlmix(i,j,k),'zlmix',0.0,.FALSE.,&
                                   & 3,indx=(/i,j,k/))
            ENDIF
         ENDDO k_310
         ENDDO j_310
         ENDDO i_310
      ENDIF
   ELSEIF (iopt_turb_param.EQ.2) THEN
      noerrs = COUNT((tke(1:ncloc,1:nrloc,2:nz).GT.0.0).AND.&
                    &(dissip(1:ncloc,1:nrloc,2:nz).LE.0.0))
      IF (noerrs.GT.0) THEN
         i_320: DO i=1,ncloc
         j_320: DO j=1,nrloc
         k_320: DO k=2,nz
            IF ((tke(i,j,k).GT.0.0).AND.(dissip(i,j,k).LE.0.0)) THEN
               CALL error_lbound_arr(dissip(i,j,k),'dissip',0.0,.FALSE.,&
                                   & 3,indx=(/i,j,k/))
            ENDIF
         ENDDO k_320
         ENDDO j_320
         ENDDO i_320
      ENDIF
   ENDIF

   CALL error_abort('eddy_coefs_tc',ierrno_runval)

ENDIF

!
!4. Eddy coefficients
!--------------------
!

k_400: DO k=2,nz

!   
!4.1 Stability functions for momentum and scalars
!------------------------------------------------
!
!4.1.1 Uniform
!-------------
!

   IF (iopt_turb_stab_form.EQ.1) THEN
      WHERE (maskatc_int)
         fstabmom = f0stabmom
         fstabscal = f0stabscal
      END WHERE

!
!4.1.2 Munk-Anderson
!-------------------
!

   ELSEIF (iopt_turb_stab_form.EQ.2) THEN
      WHERE (maskatc_int)
         ricnum = richardson_number(k,maskatc_int)
         fstabmom = f0stabmom*mom_strat_MA(ricnum,maskatc_int)
         fstabscal = f0stabscal*scal_strat_MA(ricnum,maskatc_int)
      END WHERE

!
!4.1.3 Second-order closure
!--------------------------
!

   ELSEIF (iopt_turb_stab_form.EQ.3) THEN

!     ---stability parameters
      IF (iopt_turb_param.EQ.1) THEN
         WHERE (maskatc_int)
            array2dc1 = (zlmix(1:ncloc,1:nrloc,k)/eps0)**2&
                      & /tke(1:ncloc,1:nrloc,k)
         END WHERE
      ELSEIF (iopt_turb_param.EQ.2) THEN
         WHERE (maskatc_int)
            array2dc1 = (tke(1:ncloc,1:nrloc,k)/dissip(1:ncloc,1:nrloc,k))**2
         END WHERE
      ENDIF
      WHERE (maskatc_int)
         alphaN = buofreq2(:,:,k)*array2dc1
         alphaM = shearfreq2(:,:,k)*array2dc1
      END WHERE

!     ---limiting condition (unstable)
      j_4131: DO j=1,nrloc
      i_4131: DO i=1,ncloc
         IF (maskatc_int(i,j)) alphaN(i,j) = MAX(alphaN_min,alphaN(i,j))
      ENDDO i_4131
      ENDDO j_4131

!     ---limiting condition (stable)
      IF (iopt_turb_stab_lev.EQ.1) THEN
         j_4132: DO j=1,nrloc
         i_4132: DO i=1,ncloc
            IF (maskatc_int(i,j)) alphaN(i,j) = MIN(alphaN_max,alphaN(i,j))
         ENDDO i_4132
         ENDDO j_4132
      ELSEIF ((iopt_turb_iwlim.EQ.1).AND.(iopt_turb_stab_lev.EQ.2)) THEN
         j_4133: DO i=1,ncloc
         i_4133: DO j=1,nrloc
            IF (maskatc_int(i,j)) THEN
               alphaN(i,j) = MIN(alphaN_sl*alphaM(i,j)+xmax,alphaN(i,j))
            ENDIF
         ENDDO i_4133
         ENDDO j_4133
      ENDIF

!     ---quasi-equilibrium method, c22b = 0
      IF ((iopt_turb_stab_lev.EQ.1).AND.(ib22.EQ.0)) THEN
         WHERE (maskatc_int)
            array2dc1 = 1.0+cfstab1(4)*alphaN
            fstabmom = (cfstab1(1)+cfstab1(2)*alphaN)&
                     & /(array2dc1*(1.0+cfstab1(5)*alphaN))
            fstabscal = cfstab1(3)/array2dc1
         END WHERE

!     ---quasi-equilibrium method, c22b /= 0
      ELSEIF ((iopt_turb_stab_lev.EQ.1).AND.(ib22.EQ.1)) THEN
         WHERE (maskatc_int)
            coef(:,:,1) = cfstab1(5)
            coef(:,:,2) = cfstab1(3)+cfstab1(4)*alphaN
            coef(:,:,3) = (cfstab1(1)+cfstab1(2)*alphaN)*alphaN
         END WHERE
         CALL quadratic_root_arr(coef,root1,root2,(/ncloc,nrloc,1,1/),0.0,&
                              & .TRUE.,maskatc_int)
         WHERE (maskatc_int)
            fstabscal = REAL(root1)
            fstabmom = (cfstab1(6)-cfstab1(7)*alphaN*fstabscal)&
                     & /(1.0+cfstab1(8)*alphaN)
         END WHERE

!     ---non-equilibrium method
      ELSEIF (iopt_turb_stab_lev.EQ.2) THEN
         WHERE (maskatc_int)
            array2dc1 = 1.0+cfstab2(7)*alphaM+cfstab2(8)*alphaN+&
                      & cfstab2(9)*alphaM*alphaM+cfstab2(10)*alphaM*alphaN+&
                      & cfstab2(11)*alphaN*alphaN
            fstabmom = (cfstab2(1)+cfstab2(2)*alphaM+cfstab2(3)*alphaN)&
                     & /array2dc1
            fstabscal = (cfstab2(4)+cfstab2(5)*alphaM+cfstab2(6)*alphaN)&
                      & /array2dc1
         END WHERE
      ENDIF

   ENDIF

!
!4.2 Stability coefficient for t.k.e.
!------------------------------------
!

   IF (iopt_turb_stab_tke.EQ.1) THEN
      WHERE (maskatc_int) fstabtke = f0stabtke
   ELSEIF (iopt_turb_stab_tke.EQ.2) THEN
      WHERE (maskatc_int) fstabtke = fstabmom/sigma_k
   ELSEIF (iopt_turb_stab_tke.EQ.3) THEN
      IF (iopt_turb_stab_lev.EQ.1) THEN
         WHERE (maskatc_int)
            fstabtke = cfstabtke(1) + cfstabtke(2)*alphaN*fstabscal
         END WHERE
      ELSEIF (iopt_turb_stab_lev.EQ.2) THEN
         WHERE (maskatc_int)
            fstabtke = cfstabtke(1) + cfstabtke(2)*alphaM*fstabmom + &
                     & cfstabtke(3)*alphaN*fstabscal
         END WHERE
      ENDIF
   ENDIF

!
!4.3 Eddy coefficients
!---------------------
!

   IF (iopt_turb_param.EQ.1) THEN
      WHERE (maskatc_int)
         array2dc1 = SQRT(tke(1:ncloc,1:nrloc,k))*zlmix(1:ncloc,1:nrloc,k)/eps0
      END WHERE
   ELSEIF (iopt_turb_param.EQ.2) THEN
      WHERE (maskatc_int)
         array2dc1 = (tke(1:ncloc,1:nrloc,k)/dissip(1:ncloc,1:nrloc,k))*&
                    & tke(1:ncloc,1:nrloc,k)
      END WHERE
   ENDIF
   WHERE (maskatc_int)
      array2dc2 = MERGE(vdifmom_cst,kinvisc(1:ncloc,1:nrloc,k),&
                      & iopt_turb_kinvisc.EQ.0)
      vdifcoefmom(1:ncloc,1:nrloc,k) = fstabmom*array2dc1 + array2dc2
      vdifcoefscal(:,:,k) = fstabscal*array2dc1 + vdifscal_cst
      vdifcoeftke(:,:,k) = fstabtke*array2dc1 + array2dc2
   END WHERE

ENDDO k_400

!
!5. Surface and bottom values
!----------------------------
!

WHERE (maskatc_int)
!  ---surface
   array2dc1 = (ckar*zrough_sur)*SQRT(sstresatc)
   array2dc2 = MERGE(vdifmom_cst,kinvisc(1:ncloc,1:nrloc,nz+1),&
                   & iopt_turb_kinvisc.EQ.0)
   vdifcoefmom(1:ncloc,1:nrloc,nz+1) = array2dc1 + array2dc2
   vdifcoefscal(:,:,nz+1) = (f0stabscal/f0stabmom)*array2dc1 + vdifscal_cst
   vdifcoeftke(:,:,nz+1) = (f0stabtke/f0stabmom)*array2dc1 + array2dc2
!  ---bottom
   array2dc1 = (ckar*zrough_bot)*SQRT(bstresatc(1:ncloc,1:nrloc))
   array2dc2 = MERGE(vdifmom_cst,kinvisc(1:ncloc,1:nrloc,1),&
                   & iopt_turb_kinvisc.EQ.0)
   vdifcoefmom(1:ncloc,1:nrloc,1) = array2dc1 + array2dc2
   vdifcoefscal(:,:,1) = (f0stabscal/f0stabmom)*array2dc1 + vdifscal_cst
   vdifcoeftke(:,:,1) = (f0stabtke/f0stabmom)*array2dc1 + array2dc2
END WHERE

!
!6. KPP internal wave mixing scheme
!----------------------------------
!

IF (iopt_turb_iwlim.EQ.2) THEN
   k_610: DO k=2,nz
      WHERE (maskatc_int.AND.tke(1:ncloc,1:nrloc,k).LE.tkelim)
         ricnum = richardson_number(k,maskatc_int)
         array2dc1 = MERGE(MAX(0.0,(1.0-(ricnum/riccrit_iw)**2)**3),1.0,&
                                & ricnum.GT.0.0)*vdifshear_iw
         vdifcoefmom(1:ncloc,1:nrloc,k) = vdifcoefmom(1:ncloc,1:nrloc,k) +&
                                        & array2dc1 + vdifmom_iw
         vdifcoefscal(:,:,k) = vdifcoefscal(:,:,k) + array2dc1 + vdifscal_iw
      END WHERE
   ENDDO k_610
ENDIF

!
!7. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = (/1,0,1,0/)
   CALL exchange_mod(vdifcoefmom,(/0,0,1/),nhexch,iarr_vdifcoefmom)
ENDIF

!
!8. Deallocate arrays
!--------------------
!

DEALLOCATE (alphaM,alphaN,array2dc1,array2dc2,fstabmom,fstabscal,fstabtke,&
          & ricnum,root1,root2,coef)

CALL log_timer_out()


RETURN

END SUBROUTINE eddy_coefs_tc

!========================================================================

SUBROUTINE init_turbulence
!************************************************************************
!
! *init_turbulence* Initialise parameters for turbulence closures
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.7.1
!
! Description - 
!
! Reference -
!
! Calling program - default_phsics, vertical_diff_coefs
!
! Module calls - error_abort, poly_root, quadratic_root_var
!
!************************************************************************
!
USE iopars
USE physpars
USE switches
USE turbpars
USE error_routines, ONLY: error_abort
USE math_library, ONLY: poly_root, quadratic_root_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=15) :: cval
INTEGER :: ierr, npcc
REAL :: a1, a21, a22, a23, a3, b1, b21, b22, c1a, c1b, c1my, c21a, c21b, c22a,&
      & c22b, c23a, c3a, c3b, gam, q0, rbeta, rstar, x, xfac1, xfac2, y
COMPLEX :: root1, root2, z
COMPLEX (KIND=kndclong) :: z_double
REAL, DIMENSION(7) :: coef
REAL(KIND=kndrlong), DIMENSION(7) :: coef_double
REAL, PARAMETER :: alphaN_init = 4.5, a1my = 0.92, a2my = 0.74, b1my = 16.6,&
                 & b2my = 10.1, ratlim = 1.0, x13 = 1.0/3.0, x23 = 2.0/3.0

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*ierr*       INTEGER Error code number
!*c1my*       REAL    Mellor-Yamada constant C_1
!*alphaN_init*REAL    Initial value used to solve polynomial equation for
!                     alphaN_max
!*a1my*       REAL    Mellor_Yamada constant A_1
!*a2my*       REAL    Mellor_Yamada constant A_2
!*b1my*       REAL    Mellor_Yamada constant B_1
!*b2my*       REAL    Mellor_Yamada constant B_2
!*ratlim*     REAL    Critical ratio of Thorpe to Ozmidov scale
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'init_turbulence'
CALL log_timer_in(npcc)

!
!1. Model parameters
!-------------------
!

SELECT CASE(iopt_turb_stab_mod)
!  ---Mellor-Yamada
   CASE (1)
      c1my = x13 - 2.0*a1my/b1my - x13/(a1my*(b1my**x13))
      c1a = b1my/(6.0*a1my)
      c21a = 0.0; c22a = -2.0*c1my; c23a = 0.0; c3a = 0.0
      c1b = b1my/(6.0*a2my); c21b = 0.0; c22b = 0.0; c3b = 0.0
      rbeta = b2my/b1my
      ib22 = 0
!  ---Kantha_Clayson
   CASE (2)
      c1my = x13 - 2.0*a1my/b1my - x13/(a1my*(b1my**x13))
      c1a = b1my/(6.0*a1my)
      c21a = 0.0; c22a = -2.0*c1my; c23a = 0.0; c3a = 0.0
      c1b = b1my/(6.0*a2my); c21b = 0.7; c22b = 0.0; c3b = 0.2
      rbeta = b2my/b1my
      ib22 = 0
!  ---Burchard-Baumert
   CASE (3)
      c1a = 1.8; c21a = 0.6; c22a = 0.0; c23a = 0.0; c3a = 0.6
      c1b = 3.0; c21b = 0.33; c22b = 0.0; c3b = 0.333; rbeta = 0.8
      ib22 = 0
!  ---Hossain-Rodi
   CASE (4)
      c1a = 2.2; c21a = 0.55; c22a = 0.0; c23a = 0.0; c3a = 0.55
      c1b = 3.0; c21b = 0.5; c22b = 0.0; c3b = 0.5; rbeta = 0.8
      ib22 = 0
!  ---Canuto A
   CASE (5)
      c1a = 2.49; c21a = 0.777; c22a = 0.256; c23a = 0.207; c3a = 0.402
      c1b = 5.95; c21b = 0.8; c22b = 0.2; c3b = 0.333; rbeta = 0.72
      ib22 = 1
!  ---Canuto B
   CASE (6)
      c1a = 2.1; c21a = 0.803; c22a = 0.257; c23a = 0.183; c3a = 0.576
      c1b = 5.6; c21b = 0.8; c22b = 0.2; c3b = 0.333; rbeta = 0.477
      ib22 = 1
END SELECT

!---shortcut notations
a1 = 1.0/c1a
a21 = a1*(1.0-c21a)
a22 = a1*c22a
a23 = a1*c23a
a3 = a1*(1.0-c3a)
b1 = 1.0/c1b
b21 = b1*(1.0-c21b)
b22 = b1*c22b
gam = 2.0*rbeta*(1.0-c3b)
q0 = a21*a21+a23*a23+4.0*a21*a23

!
!2. Coefficients for equilibrium method
!--------------------------------------
!

cfequil(1) = x23*(q0+a23-a21)-a22-b21*b22
cfequil(2) = b1*(x23+7.0*a3/3.0+gam)
cfequil(3) = b21*b22*(x23*(a21-a23-q0)+a22)
cfequil(4) = x23*b1*(2.0*a23*(a21+a23+a3*(1.0+a21))-a22*(a21+2.0*a23+2.0*a3)+&
           & a3*b21*(1.0-a21-2.0*a23)+b22*(a21-a23+1.5*a22+a3*(2.0*a21+a23))-&
           & gam*(a21-a23-q0+1.5*a22))
cfequil(5) = b1*b1*a3*(x23*(1.0+2.0*a3)+gam)

!
!3. Coefficients for stability functions
!---------------------------------------
!
!---quasi-equilibrium, c22b = 0
IF (ib22.EQ.0) THEN
   cfstab1(1) = x23*(a21-a23-q0)+a22
   cfstab1(2) = x23*b1*(a22*(a21+2.0*a23)-2.0*a23*(a21+a23)+&
              & 2.0*a3*(a22-a23-a21*a23)-b21*a3*(1.0-a21-2.0*a23)+&
              & gam*(a21-a23+1.5*a22-q0))
   cfstab1(3) = x23*b1*(1.0-a21-2.0*a23)
   cfstab1(4) = b1*(x23*(a21+2.0*a23+2.0*a3)+gam)
   cfstab1(5) = b1*a3

!---quasi-equilibrium, c22b/=0
ELSEIF (ib22.EQ.1) THEN
   cfstab1(1) = (x23*(q0+a3*(2.0*a21+a23))+a3*b21+b21*b22)/a1/b1
   cfstab1(2) = x23*(x23*(q0*(a21+2.0*a23)+a3*(4.0*q0-3.0*a21*a23)+&
              & 2.0*a3*a3*(2.0*a21+a23))+b21*a3*(a21+2.0*a23+2.0*a3)-&
              & b22*(q0+a3*(2.0*a21+a23))+gam*(q0+a3*(2.0*a21+a23+1.5*b21)))/a1
   cfstab1(3) = (x23*(q0-a21+a23)-a22+b21*b22)/a1/b1
   cfstab1(4) = x23*(x23*(2.0*q0*(a21+2.0*a23)-2.0*a21*a21+a23*a23-&
              & 5.0*a21*a23+a3*(a23-4.0*a21+4.0*q0-3.0*a21*a23))-&
              & a22*(a21+2.0*a23+2.0*a3)+a3*b21*(a21+2.0*a23-1.0)-&
              & b22*(a23-a21-1.5*a22+2.0*q0+a3*(2.0*a21+a23))+&
              & gam*(q0-a21+a23-1.5*a22))/a1
   cfstab1(5) = x23*(x23*(a21-a23-2.0*a21*a21+a23*a23-5.0*a21*a23+&
              & q0*(a21+2.0*a23))+a22*(1.0-a21-2.0*a23)+&
              & b22*(a21-a23-q0+1.5*a22))/a1
   cfstab1(6) = x23*(a21-a23-q0)+a22
   cfstab1(7) = x23*(q0+a3*(2.0*a21+a23))+b21*a3
   cfstab1(8) = b1*a3
ENDIF

!---non-equilibrium
cfstab2(1) = x23*(a21-a23)+a22
cfstab2(2) = -b21*b22*(x23*(a21-a23)+a22)
cfstab2(3) = x23*b1*(2.0*a3*(a22-a23)-a3*b21+gam*(a21-a23+1.5*a22))
cfstab2(4) = x23*b1
cfstab2(5) = x23*b1*(2.0*a23*(a21+a23)-a22*(a21+2.0*a23)+b22*(a21-a23+1.5*a22))
cfstab2(6) = x23*a3*b1*b1
cfstab2(7) = x23*q0-b21*b22
cfstab2(8) = b1*(7.0*a3/3.0+gam)
cfstab2(9) = -x23*b21*b22*q0
cfstab2(10)= x23*b1*(a3*(2.0*a21*a23+(2.0*a21+a23)*b22-&
           & b21*(a21+2.0*a23))+gam*q0)
cfstab2(11)= a3*b1*b1*(2.0*x23*a3+gam)

!---coefficients for t.k.e
IF (iopt_turb_stab_lev.EQ.1) THEN
   cfstabtke(1) = x23*c_sk*(1.0-a21-2.0*a23)
   cfstabtke(2) = -x23*c_sk*(a21+2.0*a23+2.0*a3)
ELSEIF (iopt_turb_stab_lev.EQ.2) THEN
   cfstabtke(1) = x23*c_sk
   cfstabtke(2) = -x23*c_sk*(a21+2.0*a23)
   cfstabtke(3) = -2.0*x23*c_sk*a3
ENDIF

!
!4. Parameters for limiting conditions
!-------------------------------------
!
!4.1 Unstable case
!-----------------
!
!---quasi-equilibrium, c22b = 0
IF ((iopt_turb_stab_lev.EQ.1).AND.(ib22.EQ.0)) THEN
   alphaN_min = -1.0/(x23*(1.0+2.0*a3)+gam)/b1

!---quasi-equilibrium c22b /= 0
ELSEIF ((iopt_turb_stab_lev.EQ.1).AND.(ib22.EQ.1)) THEN
   coef(1) = cfstab1(1)*(cfstab1(1)-cfstab1(3))
   coef(2) = cfstab1(2)*(2.0*cfstab1(1)-cfstab1(3))+&
           & cfstab1(1)*(cfstab1(5)-cfstab1(4))
   coef(3) = cfstab1(2)*(cfstab1(2)-cfstab1(4)+cfstab1(5))
   CALL quadratic_root_var(coef(1:3),root1,root2,float_fill)
   IF (ABS(REAL(root1)-float_fill).GT.float_fill_eps) THEN
      alphaN_min = MAX(REAL(root1),REAL(root2))
   ENDIF

!---non-equilibrium
ELSEIF (iopt_turb_stab_lev.EQ.2) THEN
   alphaN_min = -1.0/(2.0*a3+gam)/b1
ENDIF

!
!4.2 Stable case
!---------------
!

rstar = 0.5*ratlim*ratlim/rbeta

!---(large) numerical limit
IF (iopt_turb_iwlim.NE.1) THEN
   alphaN_max = 1.0E+11

!---quasi-equilibrium, c22b = 0
ELSEIF (ib22.EQ.0) THEN
   coef(1) = -rstar
   coef(2) = 0.0
   coef(3) = -rstar*cfstab1(4)
   coef(4) = cfstab1(3)
   alphaN_max = alphaN_init
   z = CMPLX(alphaN_max)
   coef_double = coef
   z_double = z
   CALL poly_root(coef_double(1:4),z_double,3,ierr)
   z = z_double
   IF (ierr.NE.0.OR.ABS(AIMAG(z)).NE.0.0) GOTO 1001
   alphaN_max = REAL(z)**2

!---quasi-equilibrium, c22b /= 0
ELSEIF (ib22.EQ.1) THEN
   coef(1) = (rstar*cfstab1(1))**2
   coef(2) = rstar*cfstab1(1)*cfstab1(3)
   coef(3) = 2.0*rstar*rstar*cfstab1(1)*cfstab1(2)
   coef(4) = rstar*(cfstab1(1)*cfstab1(4)+cfstab1(2)*cfstab1(3))
   coef(5) = (rstar*cfstab1(2))**2 + cfstab1(1)*cfstab1(5)
   coef(6) = rstar*cfstab1(2)*cfstab1(4)
   coef(7) = cfstab1(2)*cfstab1(5)
   alphaN_max = alphaN_init
   z = CMPLX(alphaN_max)
   coef_double = coef
   z_double = z
   CALL poly_root(coef_double(1:7),z_double,6,ierr)
   z = z_double
   IF (ierr.NE.0.OR.ABS(AIMAG(z)).NE.0.0) GOTO 1001
   alphaN_max = REAL(z)**2
ENDIF

!---non-equilibrium
IF (iopt_turb_stab_lev.EQ.2) THEN
   x = alphaN_max
   coef(1) = rstar*(1.0+cfstab2(8)*x+cfstab2(11)*x*x)-&
           & (cfstab2(4)+cfstab2(6)*x)*x**1.5
   coef(2) = rstar*(cfstab2(7)+cfstab2(10)*x)-cfstab2(5)*x**1.5
   coef(3) = rstar*cfstab2(9)
   IF (ib22.EQ.0) THEN
      alphaM_max = -coef(1)/coef(2)
   ELSE
      CALL quadratic_root_var(coef(1:3),root1,root2,float_fill)
      IF (ABS(REAL(root1)-float_fill).GT.float_fill_eps) THEN
         alphaM_max = MIN(REAL(root1),REAL(root2))
      ENDIF
   ENDIF
   y = alphaM_max
   xfac1 = rstar*(cfstab2(7)+2.0*cfstab2(9)*y+cfstab2(10)*x)-cfstab2(5)*x**1.5
   xfac2 = rstar*(cfstab2(8)+cfstab2(10)*y+2.0*cfstab2(11)*x)-&
         & 0.5*SQRT(x)*(3.0*(cfstab2(4)+cfstab2(5)*y)+5.0*cfstab2(6)*x)
   alphaN_sl = -xfac1/xfac2
ENDIF

!
!5. Neutral stability coefficients
!---------------------------------
!
!---momentum
f0stabmom = x23*(a21-a23-q0)+a22
eps0 = f0stabmom**0.75

!---scalars
xfac1 = x23*b1*(x23*(a21-a23-2.0*a21*a21+a23*a23-5.0*a21*a23+&
      & q0*(a21+2.0*a23))+a22*(1.0-a21-2.0*a23)+b22*(a21-a23-q0+1.5*a22))
xfac2 = x23*(q0-a21+a23)-a22+b21*b22
f0stabscal = -xfac1/xfac2

!---t.k.e.
IF (iopt_turb_stab_tke.EQ.1) THEN
   IF (iopt_turb_param.EQ.1) THEN
      f0stabtke = SQRT(2.0)*eps0*sq_my
   ELSE
      f0stabtke = skeps
   ENDIF
ELSEIF (iopt_turb_stab_tke.EQ.2) THEN
   f0stabtke = f0stabmom/sigma_k
ELSEIF (iopt_turb_stab_tke.EQ.3) THEN
   f0stabtke = x23*c_sk*(1.0-a21-2.0*a23)
ENDIF

!
!6. Schmidt number in transport equations
!----------------------------------------
!
!---kl-equation
sigma_kl = (2.0*f0stabtke*(ckar/eps0)**2)/(1.0+e2_my-e1_my)
!---epsilon-equation
sigma_eps = ((ckar/eps0)**2)*f0stabtke/(c2_eps-c1_eps)

CALL log_timer_out(npcc,itm_init)


RETURN

1001 nerrs = 1
IF (errchk) THEN
   WRITE (cval,'(G15.7)') ABS(AIMAG(z)); cval = ADJUSTL(cval)
   WRITE (ioerr,'(A)') 'Imaginary part of alphaN_max must be zero: '&
                     & //TRIM(cval)
ENDIF
CALL error_abort('init_turbulence',ierrno_runval)

RETURN

END SUBROUTINE init_turbulence

!========================================================================

SUBROUTINE kl_equation
!************************************************************************
!
! *kl_equation* Solve kl-(q^2*l)equation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - vertical_diff_coefs
!
! External calls - transport_at_W
!
! Module calls - error_alloc, exchange_mod
!
!************************************************************************
!
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE turbpars
USE turbulence
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: k
REAL :: eps0fac, e1fac, e3fac
INTEGER, DIMENSION(4) :: nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: buoprod, tkezlbot, tkezlsur, &
                                         & wallfunc, zbot, ztop
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: sink, source, tkezl, vdifcoeftkezl

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*buoprod*      REAL    Buoyancy production term                         [W/kg]
!*fluxbot*      REAL    Bottom flux                                   [m^4/s^3]
!*fluxsur*      REAL    Surface flux                                  [m^4/s^3]
!*tkezlbot*     REAL    Bottom value of t.k.e. times l                 [J*m/kg]
!*tkezlsur*     REAL    Surface value of t.k.e. times l                [J*m/kg]
!*wallfunc*     REAL    Wall proximity function
!*zbot*         REAL    Near-bottom length scale                            [m]
!*ztop*         REAL    Near-surface length scale                           [m]
!*sink*         REAL    Sink terms in t.k.e. equation                     [1/s]
!*source*       REAL    Source terms in t.k.e. equation                [W*m/kg]
!*tkezl*        REAL    t.k.e. times l                                 [J*m/kg]
!*vdifcoeftkezl*REAL    Vertical diffusion coefficient in transport equation
!                                                                       [m^2/s]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'kl_equation'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (buoprod(ncloc,nrloc),STAT=errstat)
CALL error_alloc('buoprod',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sink(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('sink',3,(/ncloc,nrloc,nz-1/),kndrtype)
ALLOCATE (source(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('source',3,(/ncloc,nrloc,nz-1/),kndrtype)
ALLOCATE (tkezl(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1),STAT=errstat)
CALL error_alloc('tkezl',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz+1/),kndrtype)
ALLOCATE (tkezlbot(ncloc,nrloc),STAT=errstat)
CALL error_alloc('tkezlbot',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (tkezlsur(ncloc,nrloc),STAT=errstat)
CALL error_alloc('tkezlsur',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (vdifcoeftkezl(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('vdifcoeftkezl',3,(/ncloc,nrloc,nz+1/),kndrtype)
ALLOCATE (wallfunc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('waalfunc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (zbot(ncloc,nrloc),STAT=errstat)
CALL error_alloc('zbot',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (ztop(ncloc,nrloc),STAT=errstat)
CALL error_alloc('ztop',2,(/ncloc,nrloc/),kndrtype)

!
!2. Initialise
!-------------
!
!---diffusion coefficients
k_210: DO k=1,nz+1
   WHERE (.NOT.maskatc_int)
      vdifcoeftkezl(:,:,k) = 0.0
   END WHERE
ENDDO k_210

!---surface/bottom values
IF (nt.EQ.1) THEN
   WHERE (maskatc_int)
      zlmix(1:ncloc,1:nrloc,nz+1) = ckar*zrough_sur
      zlmix(1:ncloc,1:nrloc,1) = ckar*zrough_bot
   END WHERE
ENDIF

!
!3. Store old values
!-------------------
!

tkezl = tke_old*zlmix

!
!4. Boundary conditions
!----------------------
!
!4.1 Surface
!-----------
!

WHERE (maskatc_int)
   tkezlsur = ckar*tke(1:ncloc,1:nrloc,nz)*&
            & (zrough_sur+delzatc(1:ncloc,1:nrloc,nz))
END WHERE

!
!4.2 Bottom
!----------
!

WHERE (maskatc_int)
   tkezlbot = ckar*tke(1:ncloc,1:nrloc,2)*&
            & (zrough_bot+delzatc(1:ncloc,1:nrloc,1))
END WHERE

!
!5. Diffusion coefficient
!------------------------
!

k_510: DO k=1,nz+1
   WHERE (maskatc_int)
      vdifcoeftkezl(1:ncloc,1:nrloc,k) = (vdifcoeftke(1:ncloc,1:nrloc,k)-&
                                        & vdifmom_cst)/sigma_kl
   END WHERE
ENDDO k_510

!
!6. Source/sink terms
!--------------------
!

e1fac = 0.5*e1_my; e3fac = 0.5*e3_my; eps0fac = 0.5*eps0
WHERE (.NOT.maskatc_int) buoprod = 0.0

k_610: DO k=2,nz
   WHERE (maskatc_int)
      source(:,:,k) = e1fac*zlmix(1:ncloc,1:nrloc,k)*&
               & (vdifcoefmom(1:ncloc,1:nrloc,k)-vdifmom_cst)*shearfreq2(:,:,k)
      zbot = deptotatc(1:ncloc,1:nrloc)*gscoordatw(1:ncloc,1:nrloc,k)&
         & + zrough_bot
      ztop = deptotatc(1:ncloc,1:nrloc)*(1.0-gscoordatw(1:ncloc,1:nrloc,k))&
         & + zrough_sur
      wallfunc = 1.0 + e2_my*(zlmix(1:ncloc,1:nrloc,k)*(zbot+ztop)&
              & /(ckar*zbot*ztop))**2
      sink(:,:,k) = eps0fac*wallfunc*SQRT(tke_old(1:ncloc,1:nrloc,k))&
                  & /zlmix(1:ncloc,1:nrloc,k)
      buoprod = e3fac*(vdifcoefscal(:,:,k)-vdifscal_cst)*buofreq2(:,:,k)
   END WHERE
   WHERE (buoprod.LT.0.0)
      source(:,:,k) = source(:,:,k) - zlmix(1:ncloc,1:nrloc,k)*buoprod
   END WHERE
   WHERE (buoprod.GT.0.0)
      sink(:,:,k) = sink(:,:,k) + buoprod/tke_old(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_610

!
!7. Solve transport equation
!---------------------------
!

CALL transport_at_W(tkezl,source,sink,vdifcoeftkezl,4,4,tkezlsur,tkezlbot,&
                  & iarr_tkezl)

!
!8. Update mixing length
!-----------------------
!
!---interior points
k_810: DO k=2,nz
   WHERE (maskatc_int)
      zlmix(1:ncloc,1:nrloc,k) = tkezl(1:ncloc,1:nrloc,k)&
                              & /tke(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_810

!
!9. Apply lower limit (numerical)
!--------------------------------
!

k_910: DO k=1,nz+1
   WHERE (maskatc_int)
      zlmix(1:ncloc,1:nrloc,k) = MAX(zlmixmin,zlmix(1:ncloc,1:nrloc,k))
   END WHERE
ENDDO k_910

!
!10. Exchange halo sections
!--------------------------
!

IF ((iopt_MPI.EQ.1).AND.(iopt_adv_turb.GT.0.OR.iopt_hdif_turb.GT.0)) THEN
   nhexch = nhturb
   CALL exchange_mod(zlmix,(/1-nhalo,1-nhalo,1/),nhexch,iarr_zlmix,&
                   & corners=.FALSE.)
ENDIF

!
!11. Deallocate arrays
!---------------------
!

DEALLOCATE (buoprod,sink,source,tkezl,tkezlbot,tkezlsur,vdifcoeftkezl,&
          & wallfunc,zbot,ztop)

CALL log_timer_out()


RETURN

END SUBROUTINE kl_equation

!========================================================================

SUBROUTINE mixing_length
!************************************************************************
!
! *mixing_length* Algebraic mixing length formulations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - vertical_diff_coefs
!
! External calls -
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
USE turbpars
USE turbulence
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: k
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc1, array2dc2, array2dc3,&
                                         & zbot, zlasym, ztop
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dc

!
! Name     Type    Purpose
!------------------------------------------------------------------------------
!*zbot*    REAL    Near-bottom length scale                                [m]
!*zlasym*  REAL    Blackadar's asymptotic mixing length                    [m]
!*ztop*    REAL    Near-surface length scale                               [m]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'mixing_length'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (zbot(ncloc,nrloc),STAT=errstat)
CALL error_alloc('zbot',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (ztop(ncloc,nrloc),STAT=errstat)
CALL error_alloc('ztop',2,(/ncloc,nrloc/),kndrtype)
IF (iopt_turb_lmix.EQ.4) THEN
   ALLOCATE (array2dc1(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc1',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (array2dc2(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc2',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (array2dc3(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc3',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (zlasym(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('zlasym',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF ((iopt_turb_lmix.EQ.3).OR.(iopt_turb_lmix.EQ.4)) THEN
   ALLOCATE (array3dc(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('array3dc',3,(/ncloc,nrloc,nz+1/),kndrtype)
ENDIF

!
!2. Parabolic law
!----------------
!

IF (iopt_turb_lmix.EQ.1) THEN
   k_210: DO k=1,nz+1
      WHERE (maskatc_int)
         zbot = deptotatc(1:ncloc,1:nrloc)*gscoordatw(1:ncloc,1:nrloc,k)&
            & + zrough_bot
         ztop = deptotatc(1:ncloc,1:nrloc)*&
              & (1.0-gscoordatw(1:ncloc,1:nrloc,k)) + zrough_sur
         zlmix(1:ncloc,1:nrloc,k) = (ckar*zbot*ztop)/(zbot+ztop)
      END WHERE
   ENDDO k_210

!
!3. "Modified" parabolic law
!---------------------------
!

ELSEIF (iopt_turb_lmix.EQ.2) THEN
   k_310: DO k=1,nz+1
      WHERE (maskatc_int)
         zbot = deptotatc(1:ncloc,1:nrloc)*gscoordatw(1:ncloc,1:nrloc,k)&
            & + zrough_bot
         ztop = deptotatc(1:ncloc,1:nrloc)*&
              & (1.0-gscoordatw(1:ncloc,1:nrloc,k)) + zrough_sur
         zlmix(1:ncloc,1:nrloc,k) = (ckar*zbot)*SQRT(ztop/(zbot+ztop))
      END WHERE
   ENDDO k_310

!
!4. Xing formulation
!-------------------
!

ELSEIF (iopt_turb_lmix.EQ.3) THEN
   WHERE (maskatc_int)
      array3dc(:,:,1) = 1.0
      ztop = deptotatc(1:ncloc,1:nrloc) + zrough_sur
      zlmix(1:ncloc,1:nrloc,1) = (ckar*zrough_bot)*ztop/(zrough_bot+ztop)
   END WHERE
   k_410: DO k=2,nz+1
      WHERE (maskatc_int)
         array3dc(:,:,k) = array3dc(:,:,k-1)*&
              & EXP(-beta_Xing*delzatc(1:ncloc,1:nrloc,k)/&
                   & deptotatc(1:ncloc,1:nrloc))
         zbot = array3dc(:,:,k)*gscoordatw(1:ncloc,1:nrloc,k)*&
              & deptotatc(1:ncloc,1:nrloc) + zrough_bot
         ztop = deptotatc(1:ncloc,1:nrloc)*&
              & (1.0-gscoordatw(1:ncloc,1:nrloc,k)) + zrough_sur
         zlmix(1:ncloc,1:nrloc,k) = (ckar*zbot*ztop)/(zbot+ztop)
      END WHERE
   ENDDO k_410

!
!5. Using Blackadar's asymptotic length
!--------------------------------------
!

ELSEIF (iopt_turb_lmix.EQ.4) THEN

!  ---asymptotic length
   WHERE (maskatc_int)
      array2dc1 = 0.5*delzatc(1:ncloc,1:nrloc,1)*SQRT(tke(1:ncloc,1:nrloc,1))
      array2dc2 = (1.0-0.25*gscoordatw(1:ncloc,1:nrloc,2))*array2dc1
   END WHERE
   k_511: DO k=2,nz
      WHERE (maskatc_int)
         array2dc3 = delzatw(1:ncloc,1:nrloc,k)*SQRT(tke(1:ncloc,1:nrloc,k))
         array2dc1 = array2dc1 + array2dc3
         array2dc2 = array2dc2 + (1.0-gscoordatw(1:ncloc,1:nrloc,k))&
                               & *array2dc3
      END WHERE
   ENDDO k_511
   WHERE (maskatc_int)
      array2dc3 = 0.5*delzatc(1:ncloc,1:nrloc,nz)*&
                & SQRT(tke(1:ncloc,1:nrloc,nz+1))
      array2dc1 = array2dc1 + array2dc3
      array2dc2 = array2dc2 + 0.25*(1.0-gscoordatw(1:ncloc,1:nrloc,nz))*&
                & array2dc3
      zlasym = alpha_Black*deptotatc(1:ncloc,1:nrloc)*(array2dc2/array2dc1)
   END WHERE

!  ---mixing length
   k_512: DO k=1,nz+1
      WHERE (maskatc_int)
         zbot = deptotatc(1:ncloc,1:nrloc)*gscoordatw(1:ncloc,1:nrloc,k)&
            & + zrough_bot
         ztop = deptotatc(1:ncloc,1:nrloc)*&
              & (1.0-gscoordatw(1:ncloc,1:nrloc,k)) + zrough_sur
         array3dc(:,:,k) = ckar*zbot*ztop
         zlmix(1:ncloc,1:nrloc,k) = array3dc(:,:,k)&
                                  & /(zbot+ztop+array3dc(:,:,k)/zlasym)
      END WHERE
   ENDDO k_512
ENDIF

!
!6. Apply lower limit (numerical)
!--------------------------------
!

k_610: DO k=1,nz+1
   WHERE (maskatc_int)
      zlmix(1:ncloc,1:nrloc,k) = MAX(zlmixmin,zlmix(1:ncloc,1:nrloc,k))
   END WHERE
ENDDO k_610

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE (zbot,ztop)
IF (iopt_turb_lmix.EQ.4) DEALLOCATE (array2dc1,array2dc2,array2dc3,zlasym)
IF ((iopt_turb_lmix.EQ.3).OR.(iopt_turb_lmix.EQ.4)) DEALLOCATE (array3dc)

CALL log_timer_out()


RETURN

END SUBROUTINE mixing_length

!========================================================================

SUBROUTINE tke_equation
!************************************************************************
!
! *tke_equation* Solve turbulent energy equation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - vertical_diff_coefs
!
! External calls - transport_at_W
!
! Module calls - error_alloc, exchange_mod
!
!************************************************************************
!
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE turbpars
USE turbulence
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: ibcbot, ibcsur, k
REAL :: xfac
INTEGER, DIMENSION(4) :: nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: bcbot, bcsur, buoprod
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: sink, source

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*buoprod*    REAL    Buoyancy production term                           [W/kg]
!*sink*       REAL    Sink terms in t.k.e. equation                       [1/s]
!*source*     REAL    Source terms in t.k.e. equation                    [W/kg]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'tke_equation'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (bcbot(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcbot',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (bcsur(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcsur',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (buoprod(ncloc,nrloc),STAT=errstat)
CALL error_alloc('buoprod',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sink(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('sink',3,(/ncloc,nrloc,nz-1/),kndrtype)
ALLOCATE (source(ncloc,nrloc,2:nz),STAT=errstat)
CALL error_alloc('source',3,(/ncloc,nrloc,nz-1/),kndrtype)

!
!2. Initialise
!-------------
!
!---surface/bottom values on first call
xfac = SQRT(f0stabmom)
IF (nt.EQ.1) THEN
   WHERE (maskatc_int) tke(1:ncloc,1:nrloc,nz+1) = sstresatc/xfac
   WHERE (maskatc_int) tke(1:ncloc,1:nrloc,1) = bstresatc(1:ncloc,1:nrloc)/xfac
ENDIF

!---store old values
tke_old = tke

!
!3. Boundary conditions
!----------------------
!
!3.1 Surface
!-----------
!
!---Neumann
IF (iopt_turb_tke_sbc.EQ.1) THEN
   ibcsur = 1
   WHERE (maskatc_int)
      bcsur = wfltke*sstresatc**1.5
   END WHERE
!---Dirichlet
ELSEIF (iopt_turb_tke_sbc.EQ.2) THEN
   ibcsur = 3
   WHERE (maskatc_int)
      bcsur = sstresatc/xfac
   END WHERE
ENDIF

!
!3.2 Bottom
!----------
!

!---Neumann
IF (iopt_turb_tke_bbc.EQ.1) THEN
   ibcbot = 1
   WHERE (maskatc_int)
      bcbot = 0.0
   END WHERE
!---Dirichlet
ELSEIF (iopt_turb_tke_bbc.EQ.2) THEN
   ibcbot = 3
   WHERE (maskatc_int)
      bcbot = bstresatc(1:ncloc,1:nrloc)/xfac
   END WHERE
ENDIF

!
!4. Source/sink terms
!--------------------
!

WHERE (.NOT.maskatc_int) buoprod = 0.0
k_410: DO k=2,nz
   WHERE (maskatc_int)
      source(:,:,k) = (vdifcoefmom(1:ncloc,1:nrloc,k)-vdifmom_cst)*&
                     & shearfreq2(:,:,k)
      sink(:,:,k) = dissip(1:ncloc,1:nrloc,k)/tke(1:ncloc,1:nrloc,k)
      buoprod = (vdifcoefscal(:,:,k)-vdifscal_cst)*buofreq2(:,:,k)
   END WHERE
   WHERE (buoprod.LT.0.0)
      source(:,:,k) = source(:,:,k) - buoprod
   END WHERE
   WHERE (buoprod.GT.0.0)
      sink(:,:,k) = sink(:,:,k) + buoprod/tke(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_410

!
!5. Solve transport equation
!---------------------------
!

CALL transport_at_W(tke,source,sink,vdifcoeftke,ibcsur,ibcbot,bcsur,bcbot,&
                   & iarr_tke)

!
!6. Surface/bottom values
!------------------------
!
!---surface
WHERE (maskatc_int) tke(1:ncloc,1:nrloc,nz+1) = sstresatc/xfac
!---bottom
WHERE (maskatc_int) tke(1:ncloc,1:nrloc,1) = bstresatc(1:ncloc,1:nrloc)/xfac

!
!7. Apply lower limit (numerical)
!--------------------------------
!

k_610: DO k=1,nz+1
   WHERE (maskatc_int)
      tke(1:ncloc,1:nrloc,k) = MAX(tkemin,tke(1:ncloc,1:nrloc,k))
   END WHERE
ENDDO k_610

!
!8. Exchange halo sections
!-------------------------
!

IF ((iopt_MPI.EQ.1).AND.(iopt_adv_turb.GT.0.OR.iopt_hdif_turb.GT.0)) THEN
   nhexch = nhturb
   CALL exchange_mod(tke,(/1-nhalo,1-nhalo,1/),nhexch,iarr_tke,corners=.FALSE.)
ENDIF

!
!9. Deallocate arrays
!--------------------
!

DEALLOCATE (bcbot,bcsur,buoprod,sink,source)

CALL log_timer_out()


RETURN

END SUBROUTINE tke_equation

!========================================================================

SUBROUTINE tke_equilibrium
!************************************************************************
!
! *tke_equilibrium* Turbulent energy using equilibrium (level "2") method
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - vertical_diff_coefs
!
! Module calls - error_alloc, mom_strat_MA, quadratic_root_arr,
!                richardson_number, scal_strat_MA
!
!************************************************************************
!
USE fluxes
USE grid
USE gridpars
USE iopars
USE switches
USE syspars
USE turbpars
USE turbulence
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: quadratic_root_arr
USE time_routines, ONLY: log_timer_in, log_timer_out
USE turbulence_routines, ONLY: mom_strat_MA, richardson_number, scal_strat_MA

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: k
REAL :: xfac
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: ricnum
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dc, coef
COMPLEX, SAVE, ALLOCATABLE, DIMENSION(:,:) :: root1, root2


procname(pglev+1) = 'tke_equilibrium'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array3dc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('array3dc',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (coef(ncloc,nrloc,3),STAT=errstat)
CALL error_alloc('coef',3,(/ncloc,nrloc,3/),kndrtype)
ALLOCATE (ricnum(ncloc,nrloc),STAT=errstat)
CALL error_alloc('ricnum',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (root1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('root1',2,(/ncloc,nrloc/),kndcmplx)
ALLOCATE (root2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('root2',2,(/ncloc,nrloc/),kndcmplx)

!
!2. Squared turbulence time
!--------------------------
!
!2.1 Uniform value
!-----------------
!

IF (iopt_turb_stab_form.EQ.1) THEN
   k_210: DO k=2,nz
      WHERE (maskatc_int)
         array3dc(:,:,k) = f0stabmom*shearfreq2(:,:,k) - &
                         & f0stabscal*buofreq2(:,:,k)
      END WHERE
   ENDDO k_210

!
!2.2 Munk-Anderson case
!----------------------
!

ELSEIF (iopt_turb_stab_form.EQ.2) THEN
   k_220: DO k=2,nz
      WHERE (maskatc_int)
         ricnum = richardson_number(k,maskatc_int)
         array3dc(:,:,k) = f0stabmom*mom_strat_MA(ricnum,maskatc_int) - &
                         & f0stabscal*scal_strat_MA(ricnum,maskatc_int)
      END WHERE
   ENDDO k_220

!
!2.3 Second-order closure
!------------------------
!

ELSEIF (iopt_turb_stab_form.EQ.3) THEN
   coef(:,:,3) = MERGE(1.0,0.0,maskatc_int)
   k_230: DO k=2,nz
      WHERE (maskatc_int)
         coef(:,:,2) = cfequil(1)*shearfreq2(:,:,k)+cfequil(2)*buofreq2(:,:,k)
         coef(:,:,1) = cfequil(3)*shearfreq2(:,:,k)*shearfreq2(:,:,k) + &
                     & cfequil(4)*shearfreq2(:,:,k)*buofreq2(:,:,k) + &
                     & cfequil(5)*buofreq2(:,:,k)*buofreq2(:,:,k)
      END WHERE
      CALL quadratic_root_arr(coef,root1,root2,(/ncloc,nrloc,1,1/),0.0,&
                           & .TRUE.,maskatc_int)
      WHERE (maskatc_int) 
         array3dc(:,:,k) = MAX(0.0,REAL(root1),REAL(root2))
      END WHERE
   ENDDO k_230
ENDIF

!
!3. Turbulence energy
!--------------------
!
!---interior points
k_310: DO k=2,nz
   WHERE (maskatc_int)
      tke(1:ncloc,1:nrloc,k) = MAX(tkemin,array3dc(:,:,k)*&
                                & (zlmix(1:ncloc,1:nrloc,k)/eps0)**2)
   END WHERE
ENDDO k_310

!---bottom/surface
xfac = SQRT(f0stabmom)
WHERE (maskatc_int)
   tke(1:ncloc,1:nrloc,1) = bstresatc(1:ncloc,1:nrloc)/xfac
   tke(1:ncloc,1:nrloc,nz+1) = sstresatc/xfac
END WHERE

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (array3dc,coef,ricnum,root1,root2)

CALL log_timer_out()


RETURN

END SUBROUTINE tke_equilibrium

!========================================================================

SUBROUTINE zlmix_lim_conds
!************************************************************************
!
! *zlmix_lim_conds* Limiting conditions for mixing length
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.0
!
! Description - 
!
! Calling program - vertical_diff_coefs
!
! Module calls -
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE turbpars
USE turbulence
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
INTEGER :: i, j, k
REAL, PARAMETER :: epstol = 1.0E-20
REAL :: xmax


procname(pglev+1) = 'zlmix_lim_conds'
CALL log_timer_in()

!---limiting conditions
xmax = eps0*SQRT(alphaN_max)
k_110: DO k=2,nz
   j_111: DO j=1,nrloc
   i_111: DO i=1,ncloc
      IF (maskatc_int(i,j)) THEN
         tke(i,j,k) = MAX(tkelim,tke(i,j,k))
         IF (buofreq2(i,j,k).GT.epstol) THEN
            zlmix(i,j,k) = MIN(xmax*SQRT(tke(i,j,k)/buofreq2(i,j,k)),&
                             & zlmix(i,j,k))
         ENDIF
      ENDIF
   ENDDO i_111
   ENDDO j_111
ENDDO k_110

CALL log_timer_out()


RETURN

END SUBROUTINE zlmix_lim_conds

!========================================================================

SUBROUTINE zlmix_to_dissip
!************************************************************************
!
! *zlmix_to_dissip* Dissipation rate from mixing length
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Turbulence_Equations.f90  V2.0
!
! Description - 
!
! Calling program - default_phsics, vertical_diff_coefs
!
! Module calls -
!
!************************************************************************
!
USE fluxes
USE grid
USE gridpars
USE iopars
USE physpars
USE timepars
USE turbpars
USE turbulence
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
INTEGER :: i, j, k


procname(pglev+1) = 'zlmix_to_dissip'
CALL log_timer_in()


!---surface value
IF (zrough_sur.GT.0) THEN
   WHERE (maskatc_int)
      dissip(1:ncloc,1:nrloc,nz+1) = (sstresatc**1.5)/(ckar*zrough_sur)
   END WHERE
ENDIF

!---bottom value
IF (zrough_bot.GT.0) THEN
   WHERE (maskatc_int)
      dissip(1:ncloc,1:nrloc,1) = (bstresatc(1:ncloc,1:nrloc)**1.5)/&
                                & (ckar*zrough_bot)
   END WHERE
ENDIF

!---interior points
j_110: DO j=1,nrloc
i_110: DO i=1,ncloc
   IF (maskatc_int(i,j)) THEN
      k_111: DO k=2,nz
         dissip(i,j,k) = MAX(dissipmin,eps0*tke(i,j,k)**1.5/zlmix(i,j,k))
      ENDDO k_111
   ENDIF
ENDDO i_110
ENDDO j_110

CALL log_timer_out()


RETURN

END SUBROUTINE zlmix_to_dissip
