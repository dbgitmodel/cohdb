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

SUBROUTINE bottom_stress
!************************************************************************
!
! *bottom_stress* Bottom stress
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Bottom_Fluxes.f90  V2.11.1
!
! $Date: 2016-04-04 09:25:08 +0200 (Mon, 04 Apr 2016) $
!
! $Revision: 931 $
!
! Description - 
!
! Reference -
!
! Calling program - current_corr, current_corr_1d, current_2d, initialise_model
!
! Externals - weirs_bstress
!
! Module calls - error_alloc, exchange_mod, Carr_at_U, Carr_at_V, Uarr_at_C,
!                Uarr_at_V, Varr_at_C, Varr_at_U
!
!************************************************************************
!
USE currents
USE depths
USE fluxes
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE array_interp, ONLY: Carr_at_U, Carr_at_V, Uarr_at_C, Uarr_at_V, Varr_at_C, &
                      & Varr_at_U
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc
LOGICAL :: waveflag
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2, bcuratc2, &
                                         & bcuratu2, bcuratv2, bstresatc_old


procname(pglev+1) = 'bottom_stress'
CALL log_timer_in(npcc)

!
!1. Allocate and initialise
!--------------------------
!

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)
maskatu = node2du(1:ncloc,1:nrloc).GT.0

ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)
maskatv = node2dv(1:ncloc,1:nrloc).GT.0

ALLOCATE (array2d1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE (array2d2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE (bcuratc2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcuratc2',2,(/ncloc,nrloc/),kndrtype)
bcuratc2 = 0.0

ALLOCATE (bcuratu2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcuratu2',2,(/ncloc,nrloc/),kndrtype)
bcuratu2 = 0.0

ALLOCATE (bcuratv2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcuratv2',2,(/ncloc,nrloc/),kndrtype)
bcuratv2 = 0.0

IF (nt.EQ.0) THEN
   ALLOCATE (bstresatc_old(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bstresatc_old',2,(/ncloc,nrloc/),kndrtype)
   bstresatc_old = 0.0
ENDIF

waveflag = iopt_bstres_waves_bfric.GT.0.AND.&
         & (.NOT.cold_start.OR.iopt_waves_couple.EQ.0)

!
!2. Bottom currents
!------------------
!
!2.1 2-D case
!------------
!

IF (iopt_bstres_nodim.EQ.2) THEN

!  ---U-nodes
   CALL Varr_at_U(vmvel(0:ncloc,1:nrloc+1),array2d1,1,3,(/0,1,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vmvel,.TRUE.)
   WHERE (maskatu)
      bcuratu2 = umvel(1:ncloc,1:nrloc)**2+array2d1**2
      bcuratu = SQRT(bcuratu2)
   END WHERE

!  ---V-nodes
   CALL Uarr_at_V(umvel(1:ncloc+1,0:nrloc),array2d2,1,3,(/1,0,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_umvel,.TRUE.)
   WHERE (maskatv)
      bcuratv2 = array2d2**2+vmvel(1:ncloc,1:nrloc)**2
      bcuratv = SQRT(bcuratv2)
   END WHERE

!  ---C-nodes
   CALL Uarr_at_C(umvel(1:ncloc+1,1:nrloc),array2d1,1,1,(/1,1,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_umvel,.TRUE.)
   CALL Varr_at_C(vmvel(1:ncloc,1:nrloc+1),array2d2,1,1,(/1,1,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vmvel,.TRUE.)
   WHERE (maskatc_int)
      bcuratc2 = array2d1**2+array2d2**2
      bcuratc = SQRT(bcuratc2)
   END WHERE

!
!2.2 3-D case
!------------
!

ELSEIF (iopt_bstres_nodim.EQ.3) THEN

!  ---U-nodes
   CALL Varr_at_U(vvel(0:ncloc,1:nrloc+1,1),array2d1,1,3,(/0,1,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vvel,.TRUE.)
   WHERE (maskatu)
      bcuratu2 = uvel(1:ncloc,1:nrloc,1)**2+array2d1**2
      bcuratu = SQRT(bcuratu2)
   END WHERE

!  ---V-nodes
   CALL Uarr_at_V(uvel(1:ncloc+1,0:nrloc,1),array2d2,1,3,(/1,0,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_uvel,.TRUE.)
   WHERE (maskatv)
      bcuratv2 = array2d2**2+vvel(1:ncloc,1:nrloc,1)**2
      bcuratv = SQRT(bcuratv2)
   END WHERE

!  ---C-nodes
   CALL Uarr_at_C(uvel(1:ncloc+1,1:nrloc,1),array2d1,1,1,(/1,1,1/),&
               & (/ncloc+1,nrloc,1/),1,iarr_uvel,.TRUE.)
   CALL Varr_at_C(vvel(1:ncloc,1:nrloc+1,1),array2d2,1,1,(/1,1,1/),&
           & (/ncloc,nrloc+1,1/),1,iarr_vvel,.TRUE.)
   WHERE (maskatc_int)
      bcuratc2 = array2d1**2+array2d2**2
      bcuratc = SQRT(bcuratc2)
   END WHERE
   
ENDIF

!
!3. Linear bottom stress
!-----------------------
!

IF (iopt_bstres_form.EQ.1) THEN

!
!3.1. Friction velocity at initial call
!--------------------------------------
!

   IF (nt.EQ.0) THEN
      WHERE (maskatu)
         bfricatu = bdraglin
      END WHERE
      WHERE (maskatv)
         bfricatv = bdraglin
      END WHERE
   ENDIF

!
!3.2 Bottom stress
!-----------------
!
!3.2.1 2-D case
!--------------
!

   IF (iopt_bstres_nodim.EQ.2) THEN

!     ---C-nodes
      CALL Uarr_at_C(umvel(1:ncloc+1,1:nrloc),array2d1,1,1,(/1,1,nz/),&
                  & (/ncloc+1,nrloc,nz/),1,iarr_umvel,.TRUE.)
      CALL Varr_at_C(vmvel(1:ncloc,1:nrloc+1),array2d2,1,1,(/1,1,nz/),&
                  & (/ncloc,nrloc+1,nz/),1,iarr_vmvel,.TRUE.)
      WHERE (maskatc_int)
         bstresatc(1:ncloc,1:nrloc) = bdraglin*SQRT(array2d1**2+array2d2**2)
      END WHERE

!     ---U-nodes
      WHERE (maskatu)
         ubstresatu(1:ncloc,:) = bdraglin*umvel(1:ncloc,1:nrloc)
      END WHERE

!     ---V-nodes
      WHERE (maskatv)
         vbstresatv(:,1:nrloc) = bdraglin*vmvel(1:ncloc,1:nrloc)
      END WHERE

!
!3.2.2 3-D case
!--------------
!

   ELSEIF (iopt_bstres_nodim.EQ.3) THEN

!     ---C-nodes
      CALL Uarr_at_C(uvel(1:ncloc+1,1:nrloc,1),array2d1,1,1,(/1,1,1/),&
                  & (/ncloc+1,nrloc,1/),1,iarr_umvel,.TRUE.)
      CALL Varr_at_C(vvel(1:ncloc,1:nrloc+1,1),array2d2,1,1,(/1,1,1/),&
                  & (/ncloc,nrloc+1,1/),1,iarr_vmvel,.TRUE.)
      WHERE (maskatc_int)
         bstresatc(1:ncloc,1:nrloc) = bdraglin*SQRT(array2d1**2+array2d2**2)
      END WHERE

!     ---U-nodes
      WHERE (maskatu)
         ubstresatu(1:ncloc,1:nrloc) = bdraglin*uvel(1:ncloc,1:nrloc,1)
      END WHERE

!     ---V-nodes
      WHERE (maskatv)
         vbstresatv(1:ncloc,1:nrloc) = bdraglin*uvel(1:ncloc,1:nrloc,1)
      END WHERE

   ENDIF

!
!3.3 Exchange halo sections
!--------------------------
!

   IF (iopt_MPI.EQ.1) THEN
      CALL exchange_mod(ubstresatu,(/1,1/),(/0,1,0,0/),iarr_ubstresatu)
      CALL exchange_mod(vbstresatv,(/1,1/),(/0,0,0,1/),iarr_vbstresatv)
   ENDIF

   GOTO 1000

ENDIF

!
!4. Bottom drag coefficient
!--------------------------
!


IF (iopt_bstres_drag.GT.2.OR.nt.EQ.0) THEN

!
!4.1 Uniform drag coefficients
!-----------------------------
!

   IF (iopt_bstres_drag.LE.2) THEN
      CALL Carr_at_U(bdragcoefatc(:,1:nrloc),bdragcoefatu,1,1,(/0,1,nz/),&
                  & (/ncloc,nrloc,nz/),1,iarr_bdragcoefatc,.TRUE.)
      CALL Carr_at_V(bdragcoefatc(1:ncloc,:),bdragcoefatv,1,1,(/1,0,nz/),&
                  & (/ncloc,nrloc,nz/),1,iarr_bdragcoefatc,.TRUE.)
      
!
!4.2 Non-uniform drag coefficient
!--------------------------------
!

   ELSE

!         
!4.2.1 C-nodes
!-------------
!         
      
      IF (iopt_bstres_nodim.EQ.2) THEN
         WHERE (maskatc_int)
            array2d1 = MAX(deptotatc(1:ncloc,1:nrloc)/&
                        & (enap*zroughatc(1:ncloc,1:nrloc)),zbtoz0lim)
         END WHERE
      ELSEIF (iopt_bstres_nodim.EQ.3) THEN
         WHERE (maskatc_int)
            array2d1 = MAX(0.5*delzatc(1:ncloc,1:nrloc,1)/&
                         & zroughatc(1:ncloc,1:nrloc),zbtoz0lim)
         END WHERE
      ENDIF
      WHERE (maskatc_int)
         bdragcoefatc(1:ncloc,1:nrloc) = (ckar/LOG(array2d1))**2
      END WHERE

!         
!4.2.2 U-nodes
!-------------         
!

      IF (.NOT.waveflag) THEN

         IF (iopt_bstres_nodim.EQ.2) THEN
            WHERE (nodeatu(1:ncloc,1:nrloc,1).GT.1)
               array2d1 = MAX(deptotatu(1:ncloc,1:nrloc)/(enap*zroughatu),&
                            & zbtoz0lim)
            END WHERE
         ELSEIF (iopt_bstres_nodim.EQ.3) THEN
            WHERE (nodeatu(1:ncloc,1:nrloc,1).GT.1)
               array2d1 = MAX(0.5*delzatu(1:ncloc,1:nrloc,1)/zroughatu,&
                            & zbtoz0lim)
            END WHERE
         ENDIF
         WHERE (nodeatu(1:ncloc,1:nrloc,1).GT.1)
            bdragcoefatu = (ckar/LOG(array2d1))**2
         ELSEWHERE
            bdragcoefatu = 0.0
         END WHERE

!
!4.2.3 V-nodes
!-----------
!

         IF (iopt_bstres_nodim.EQ.2) THEN
            WHERE (nodeatv(1:ncloc,1:nrloc,1).GT.1)
               array2d1 = MAX(deptotatv(1:ncloc,1:nrloc)/(enap*zroughatv),&
                            & zbtoz0lim)
            END WHERE
         ELSEIF (iopt_bstres_nodim.EQ.3) THEN
            WHERE (nodeatv(1:ncloc,1:nrloc,1).GT.1)
               array2d1 = MAX(0.5*delzatv(1:ncloc,1:nrloc,1)/zroughatv,&
                            & zbtoz0lim)
            END WHERE
         ENDIF
         WHERE (nodeatv(1:ncloc,1:nrloc,1).GT.1)
            bdragcoefatv = (ckar/LOG(array2d1))**2
         ELSEWHERE
            bdragcoefatv = 0.0
         END WHERE
         
      ENDIF
      
   ENDIF
      
ENDIF

!
!5. Quadratic bottom stress
!--------------------------
!
!5.1 C-nodes
!-----------
!
!5.1.1 Without waves
!-------------------
!

IF (.NOT.waveflag) THEN

   WHERE (maskatc_int)
      bstresatc(1:ncloc,1:nrloc) = bdragcoefatc(1:ncloc,1:nrloc)*bcuratc2
   END WHERE

!
!5.1.2 With waves
!----------------
!

ELSEIF (waveflag) THEN

!  ---initialise at initial time   
   IF (nt.EQ.0) THEN
      WHERE (maskatc_int)
         bstresatc_old = bdragcoefatc(1:ncloc,1:nrloc)*bcuratc2 
      END WHERE
   ENDIF
   
!  ---wave enhanced stress
   CALL bottom_stress_waves(bstresatc_old,zroughatc,bstresatc(1:ncloc,1:nrloc),&
                          & bstresatc_max,bstresatc_wav,fwave,wavethickatc,&
                          & zaroughatc)
   
!  ---reset old value
   WHERE (maskatc_int)
      bstresatc_old = bstresatc(1:ncloc,1:nrloc)
   END WHERE
   
!  ---drag coefficient
   WHERE (bcuratc2.GT.0.0)
      bdragcoefatc(1:ncloc,1:nrloc) = bstresatc(1:ncloc,1:nrloc)/bcuratc2
   ELSEWHERE
      bdragcoefatc(1:ncloc,1:nrloc) = 0.0
   END WHERE
   
ENDIF

!
!5.1.3 Exchange halo sections
!----------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(bstresatc,(/0,0/),(/1,0,1,0/),iarr_bstresatc)
ENDIF

!
!5.2 U-nodes
!-----------
!
!---bottom stress
CALL Carr_at_U(bstresatc(:,1:nrloc),bstresatu,1,1,(/0,1,nz/),&
            & (/ncloc,nrloc,nz/),1,iarr_bstresatc,.TRUE.)

!---drag coefficient
IF (waveflag) THEN
   WHERE (bcuratu2.GT.0)
      bdragcoefatu = bstresatu/bcuratu2
   ELSEWHERE
      bdragcoefatu = 0.0
   END WHERE
ENDIF

!---bottom friction velocity
WHERE (maskatu)
   bfricatu = bdragcoefatu*bcuratu
END WHERE

!---X-component of bottom stress
IF (iopt_bstres_nodim.EQ.2) THEN
   WHERE (maskatu)
      ubstresatu(1:ncloc,:) = bfricatu*umvel(1:ncloc,1:nrloc)
   END WHERE
ELSEIF (iopt_bstres_nodim.EQ.3) THEN
   WHERE (maskatu)
      ubstresatu(1:ncloc,:) = bfricatu*uvel(1:ncloc,1:nrloc,1)
   END WHERE
ENDIF

!
!5.3 V-nodes
!-----------
!
!---bottom stress
CALL Carr_at_V(bstresatc(1:ncloc,:),bstresatv,1,1,(/1,0,nz/),&
            & (/ncloc,nrloc,nz/),1,iarr_bstresatc,.TRUE.)

!---drag coefficient
IF (waveflag) THEN
   WHERE (bcuratv2.GT.0)
      bdragcoefatv = bstresatv/bcuratv2
   ELSEWHERE
      bdragcoefatv = 0.0
   END WHERE
ENDIF
   
!---bottom friction velocity
WHERE (maskatv)
   bfricatv = bdragcoefatv*bcuratv
END WHERE

!--- Y-component of bottom stress
IF (iopt_bstres_nodim.EQ.2) THEN
   WHERE (maskatv)
      vbstresatv(:,1:nrloc) = bfricatv*vmvel(1:ncloc,1:nrloc)
   END WHERE
ELSEIF (iopt_bstres_nodim.EQ.3) THEN
   WHERE (maskatv)
      vbstresatv(:,1:nrloc) = bfricatv*vvel(1:ncloc,1:nrloc,1)
   END WHERE
ENDIF

!
!6. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(ubstresatu,(/1,1/),(/0,1,0,0/),iarr_ubstresatu)
   CALL exchange_mod(vbstresatv,(/1,1/),(/0,0,0,1/),iarr_vbstresatv)
ENDIF

!
!7. Weirs
!--------
!

IF (iopt_weibar.EQ.1) CALL weirs_bstress 

!
!8. Deallocate
!-------------
!

1000 CONTINUE
DEALLOCATE (maskatu,maskatv)
DEALLOCATE (array2d1,array2d2,bcuratc2,bcuratu2,bcuratv2)
IF (cold_start.OR.nt.EQ.nstep) DEALLOCATE (bstresatc_old)

CALL log_timer_out(npcc,itm_bconds)

RETURN

END SUBROUTINE bottom_stress

!========================================================================

SUBROUTINE bottom_stress_waves(tau_old,zrough,taumean,taumax,tauwav,fwc,wthick,&
                             & zarough)
!************************************************************************
!
! *bottom_stress_waves* calculates bottom stress and related arrays in the
!                       presence of waves using a current-wave interaction
!                       formulation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Bottom_Fluxes.f90  V2.11.1
!
! Description -
!
! Reference -
!
! Calling program - bottom_stres
!
! External calls - 
!
! Module calls -
!
!************************************************************************
!
USE currents  
USE fluxes  
USE grid
USE gridpars  
USE iopars
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: vector_mag_arr_atc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
REAL, INTENT(IN), DIMENSION(0:ncloc,0:nrloc) :: zrough
REAL, INTENT(IN), DIMENSION(ncloc,nrloc) :: tau_old
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: fwc, taumax, taumean, tauwav, &
                                           & wthick, zarough

!
! Name        Type  Purpose
!-------------------------------------------------------------------------------
!*tau_old*    REAL  Mean or maximum bottom stress under combined effects of
!                   currents and waves at old time step                  [m2/s2]
!*zrough*     REAL  Bottom roughness height                                  [m]
!*taumean*    REAL  Returned mean bottom stress under combined effects of
!                   currents and waves.                                  [m2/s2]
!*taumax*     REAL  Returned maximum bottom stress under combined effects of
!                   currents and waves.                                  [m2/s2]
!*tauwav*     REAL  Returned maximum wave bottom stress                  [m2/s2]
!*fwc*        REAL  Wave friction factor                                     [-]
!*wthick*     REAL  Thickness of the wave boundary layer                     [m]
!*zarough*    REAL  Apparent roughness height                                [m]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
INTEGER, PARAMETER :: maxiterations = 100
INTEGER :: i, iter, j, npcc
REAL ::  bstar = 2.241, epsconv = 1.0E-04, JCJcrit = 3.47, &
       & wthickpar_CJ = 0.367, wthickpar_GM = 1.0, wthickpar_MD = 0.65, &
       & wthickpar_SC = 0.65, wvel_min = 0.001, zbot_min = 0.025
REAL :: a, b, eps, JCJ, rat, ustarc2, ustarcw, ustarcw2, ustare, ustarm, &
      & ustarm2, ustarm_old, &
      & ustarp, ustarp2, ustarw, ustarw2, wthick2, xpi, xvar, zb
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: phicur, wcanglecos


procname(pglev+1) = 'bottom_stress_waves'
CALL log_timer_in(npcc)

!
!1. Cold start with wave coupling
!--------------------------------
!

IF (cold_start.AND.iopt_waves_couple.GT.0.AND.nt.EQ.0) THEN
   taumean = tau_old
   taumax = tau_old
   tauwav = 0.0
   wthick = 0.0
   zarough = zrough(1:ncloc,1:nrloc)
   GOTO 1000
ENDIF

!
!2. Current-wave angle
!---------------------
!

ALLOCATE (phicur(ncloc,nrloc),STAT=errstat)
CALL error_alloc('phicur',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (wcanglecos(ncloc,nrloc),STAT=errstat)
CALL error_alloc('wcanglecos',2,(/ncloc,nrloc/),kndrtype)

CALL vector_mag_arr_atc(uvel(1:ncloc+1,1:nrloc,1),&
                      & vvel(1:ncloc,1:nrloc+1,1),1,1,1,1,iarr_hmvelmag,&
                      & .TRUE.,vecpha=phicur)
WHERE (maskatc_int)
   wcanglecos = ABS(COS(wavedir(1:ncloc,1:nrloc)-phicur))
END WHERE

!
!3. Soulsby and Clarke
!----------------------
!

IF (iopt_bstres_waves_bfric.EQ.1) THEN

   j_310: DO j=1,nrloc
   i_310: DO i=1,ncloc   
      zb = 0.5*delzatc(i,j,1)

      IF (.NOT.maskatc_int(i,j)) THEN
         taumean(i,j) = 0.0
         taumax(i,j) = 0.0
         tauwav(i,j) = 0.0
         fwc(i,j) = 0.0
         wthick(i,j) = 0.0
         zarough(i,j) = 0.0

      ELSEIF (zb.GE.zbot_min.AND.wavevel(i,j).GE.wvel_min) THEN
         ustarc2 = bdragcoefatc(i,j)*bcuratc(i,j)*bcuratc(i,j)
         fwc(i,j) = 1.39*(wavefreq(i,j)*zrough(i,j)/wavevel(i,j))**0.52
         ustarw2 = 0.5*fwc(i,j)*wavevel(i,j)*wavevel(i,j)
         ustarw = SQRT(ustarw2)
         ustare = (ustarc2**2+ustarw2**2)**0.25
         ustarp2 = ustare*ustarw
         wthick(i,j) = wthickpar_SC*ckar*ustarw/wavefreq(i,j)
         wthick2 = wthick(i,j) + zrough(i,j)
         IF (zb.GT.wthick2) THEN
            a = LOG(wthick2/zrough(i,j))/ustare
            b = LOG(zb/wthick2)
            ustarm = (SQRT(b*b+4.0*ckar*a*bcuratc(i,j))-b)/(2.0*a)
         ELSE
            ustarm = SQRT(ckar*ustare*bcuratc(i,j)/LOG(zb/zrough(i,j)))
         ENDIF
         ustarm2 = ustarm*ustarm
         ustarcw2 = (ustarp2**2+2.0*wcanglecos(i,j)*ustarp2*ustarm2+&
                   & ustarm2**2)**0.5
         taumean(i,j) = ustarm2
         taumax(i,j) = ustarcw2
         tauwav(i,j) = ustarp2
         fwc(i,j) = 2.0*ustarp2/wavevel(i,j)**2
         rat = ustarm/ustare
         zarough(i,j) = wthick2*(zrough(i,j)/wthick2)**rat
      ELSE
         taumean(i,j) = bdragcoefatc(i,j)*bcuratc(i,j)**2
         taumax(i,j) = taumean(i,j)
         tauwav(i,j) = 0.0
         fwc(i,j) = 0.0
         wthick(i,j) = 0.0
         zarough(i,j) = zrough(i,j)

      ENDIF
      
   ENDDO i_310
   ENDDO j_310

!
!4. Malarkey and Davies
!----------------------
!

ELSEIF (iopt_bstres_waves_bfric.EQ.2) THEN

   j_410: DO j=1,nrloc
   i_410: DO i=1,ncloc   
      zb = 0.5*delzatc(i,j,1)

      IF (.NOT.maskatc_int(i,j)) THEN
         taumean(i,j) = 0.0
         taumax(i,j) = 0.0
         tauwav(i,j) = 0.0
         fwc(i,j) = 0.0
         wthick(i,j) = 0.0
         zarough(i,j) = 0.0

      ELSEIF (zb.GE.zbot_min.AND.wavevel(i,j).GE.wvel_min) THEN

         ustarc2 = bdragcoefatc(i,j)*bcuratc(i,j)*bcuratc(i,j)
         fwc(i,j) = 1.39*(wavefreq(i,j)*zrough(i,j)/wavevel(i,j))**0.52
         ustarw2 = 0.5*fwc(i,j)*wavevel(i,j)*wavevel(i,j)
         ustarw = SQRT(ustarw2)
         ustare = (ustarc2**2+2.0*wcanglecos(i,j)*ustarc2*ustarw2+&
                 & ustarw2**2)**0.25
         ustarp2 = ustare*ustarw
         wthick(i,j) = wthickpar_MD*ckar*ustarw/wavefreq(i,j)
         wthick2 = wthick(i,j) + zrough(i,j)
         IF (zb.GT.wthick2) THEN
            a = LOG(wthick2/zrough(i,j))/ustare
            b = LOG(zb/wthick2)
            ustarm = (SQRT(b*b+4.0*ckar*a*bcuratc(i,j))-b)/(2.0*a)
         ELSE
            ustarm = SQRT(ckar*ustare*bcuratc(i,j)/LOG(zb/zrough(i,j)))
         ENDIF
         ustarm2 = ustarm*ustarm
         eps = 1.0 + (ustarw2/(ustarm2+ustarw2)**2)*&
              & (wcanglecos(i,j)*ustarm2+0.25*ustarw2)
         ustarcw2 = (ustarp2**2+2.0*SQRT(eps)*wcanglecos(i,j)*ustarm2*ustarp2+&
              & eps*ustarm2**2)**0.5
         taumean(i,j) = ustarm2
         taumax(i,j) = ustarcw2
         tauwav(i,j) = ustarp2
         fwc(i,j) = 2.0*ustarp2/wavevel(i,j)**2
         rat = ustarm/ustare
         zarough(i,j) = wthick2*(zrough(i,j)/wthick2)**rat
      ELSE
         taumean(i,j) = bdragcoefatc(i,j)*bcuratc(i,j)**2
         taumax(i,j) = taumean(i,j)
         tauwav(i,j) = 0.0
         fwc(i,j) = 0.0
         wthick(i,j) = 0.0
         zarough(i,j) = zrough(i,j)

      ENDIF

   ENDDO i_410
   ENDDO j_410

!
!5. Grant and Madsen
!-------------------
!

ELSEIF (iopt_bstres_waves_bfric.EQ.3) THEN

   xpi = halfpi**2
   j_510: DO j=1,nrloc
   i_510: DO i=1,ncloc   
      zb = 0.5*delzatc(i,j,1)

      IF (.NOT.maskatc_int(i,j)) THEN
         taumean(i,j) = 0.0
         taumax(i,j) = 0.0
         tauwav(i,j) = 0.0
         fwc(i,j) = 0.0
         wthick(i,j) = 0.0
         zarough(i,j) = 0.0

      ELSEIF (zb.GE.zbot_min.AND.wavevel(i,j).GE.wvel_min) THEN
!        ---initialise
         ustarm = SQRT(tau_old(i,j))
         ustarm2 = tau_old(i,j)
         ustarcw = MAX(1.0E-08,ustarm)
         
!        ---iteration loop
         flag = .TRUE.; iter = 1
         DO WHILE (flag.AND.iter.LE.maxiterations)
            ustarm_old = ustarm
            xvar = (LOG(ckar*ustarcw/(zrough(i,j)*wavefreq(i,j)))-1.15)**2 + xpi
            ustarp2 = ckar*ustarcw*wavevel(i,j)/SQRT(xvar)
            ustarp = SQRT(ustarp2)
            ustarcw2 = (ustarp2**2+2.0*wcanglecos(i,j)*ustarp2*ustarm2+&
                      & ustarm2**2)**0.5
            ustarcw = SQRT(ustarcw2)
            wthick(i,j) =  wthickpar_GM*ckar*ustarcw/wavefreq(i,j)
            wthick2 =  wthick(i,j) + zrough(i,j)
            IF (zb.GT.wthick2) THEN
               a = LOG(wthick2/zrough(i,j))/ustarcw
               b = LOG(zb/wthick2)
               ustarm = (SQRT(b*b+4.0*ckar*a*bcuratc(i,j))-b)/(2.0*a)
            ELSE
               ustarm = SQRT(ckar*ustarcw*bcuratc(i,j)/LOG(zb/zrough(i,j)))
            ENDIF
            ustarm2 = ustarm*ustarm
            IF (ABS(ustarm-ustarm_old).GT.epsconv) THEN
               iter = iter + 1
            ELSE
               flag = .FALSE.
            ENDIF
         ENDDO
         
!        ---evaluate parameters
         taumean(i,j) = ustarm2
         taumax(i,j) = ustarcw2
         tauwav(i,j) = ustarp2
         fwc(i,j) = 2.0*ustarp2/wavevel(i,j)**2
         rat = ustarm/ustarcw
         zarough(i,j) = wthick2*(zrough(i,j)/wthick2)**rat
         
      ELSE
         taumean(i,j) = bdragcoefatc(i,j)*bcuratc(i,j)**2
         taumax(i,j) = taumean(i,j)
         tauwav(i,j) = 0.0
         fwc(i,j) = 0.0
         wthick(i,j) = 0.0
         zarough(i,j) = zrough(i,j)

      ENDIF
      
   ENDDO i_510
   ENDDO j_510
      
!
!6. Christoffersen and Johnsson
!------------------------------
!

ELSEIF (iopt_bstres_waves_bfric.EQ.4) THEN

   j_610: DO j=1,nrloc
   i_610: DO i=1,ncloc   
      zb = 0.5*delzatc(i,j,1)

      IF (.NOT.maskatc_int(i,j)) THEN
         taumean(i,j) = 0.0
         taumax(i,j) = 0.0
         tauwav(i,j) = 0.0
         fwc(i,j) = 0.0
         wthick(i,j) = 0.0
         zarough(i,j) = 0.0

      ELSEIF (zb.GE.zbot_min.AND.wavevel(i,j).GE.wvel_min) THEN
         
!
!6.2.1 Model I
!-------------
!
!        ---initialise
         ustarm = SQRT(tau_old(i,j))
         ustarm2 = ustarm*ustarm
         ustarcw = MAX(1.0E-08,ustarm)
         
!        ---iteration loop
         flag = .TRUE.; iter = 1
         DO WHILE (flag.AND.iter.LE.maxiterations)
            ustarm_old = ustarm
            ustarp2 = wavevel(i,j)*&
                    & SQRT(bstar*zrough(i,j)*wavefreq(i,j)*ustarcw)
            ustarp = SQRT(ustarp2)
            ustarcw2 = (ustarp2**2+2.0*wcanglecos(i,j)*ustarm2*ustarp2+&
                      & ustarm2**2)**0.5
            ustarcw = SQRT(ustarcw2)
            wthick(i,j) =  SQRT(bstar*zrough(i,j)*ustarcw/wavefreq(i,j))
            wthick2 =  wthick(i,j) + zrough(i,j)
            IF (zb.GT.wthick2) THEN
               a = wthick2/(zrough(i,j)*bstar*ustarcw)
               b = LOG(zb/wthick2)/ckar
               ustarm = (SQRT(b*b+4.0*a*bcuratc(i,j))-b)/(2.0*a)
            ELSE
               ustarm = SQRT(bstar*zrough(i,j)*bcuratc(i,j)*ustarcw/zb)
            ENDIF
            ustarm2 = ustarm*ustarm
            IF (ABS(ustarm-ustarm_old).GT.epsconv) THEN
               iter = iter + 1
            ELSE
               flag = .FALSE.
            ENDIF
         ENDDO
      
!        ---evaluate paramaters
         JcJ = ustarcw/(30*zrough(i,j)*wavefreq(i,j))
         taumean(i,j) = ustarm2
         taumax(i,j) = ustarcw2
         tauwav(i,j) = ustarp2
         fwc(i,j) = 2.0*ustarp2/wavevel(i,j)**2
         zarough(i,j) = wthick2*EXP(-ckar*(ustarm/ustarcw)*wthick(i,j)&
                      & /(bstar*zrough(i,j)))

!
!6.2.2 Model II
!--------------
!

         
         IF (JcJ.GT.JcJcrit) THEN
         
!           ---initialise
            ustarm = SQRT(tau_old(i,j))
            ustarm2 = ustarm*ustarm
            ustarcw = MAX(1.0E-08,ustarm)

!           ---iteration loop
            flag = .TRUE.; iter = 1
            DO WHILE (flag.AND.iter.LE.maxiterations)
               ustarm_old = ustarm
               xvar = 5.76*(LOG10(5.76*ustarcw/wavevel(i,j))-0.1164+&
                    & LOG10(wavevel(i,j)/(30*zrough(i,j)*wavefreq(i,j))))
               ustarp2 = ustarcw*wavevel(i,j)/xvar
               ustarp = SQRT(ustarp2)
               ustarcw2 = (ustarp2**2+2.0*wcanglecos(i,j)*ustarp2*ustarm2+&
                         & ustarm2**2)**0.5
               ustarcw = SQRT(ustarcw2)
               wthick(i,j) =  wthickpar_CJ*ckar*ustarcw/wavefreq(i,j)
               wthick2 =  wthick(i,j) + zrough(i,j)
               IF (zb.GT.wthick2) THEN
                  a = LOG(wthick2/zrough(i,j))/ustarcw
                  b = LOG(zb/wthick2)
                  ustarm = (SQRT(b*b+4.0*ckar*a*bcuratc(i,j))-b)/(2.0*a)
               ELSE
                  ustarm = SQRT(ckar*ustarcw*bcuratc(i,j)/LOG(zb/zrough(i,j)))
               ENDIF
               ustarm2 = ustarm*ustarm
               IF (ABS(ustarm-ustarm_old).GT.epsconv) THEN
                  iter = iter + 1
               ELSE
                  flag = .FALSE.
               ENDIF
            ENDDO
         
!           ---evaluate parameters
            taumean(i,j) = ustarm2
            taumax(i,j) = ustarcw2
            tauwav(i,j) = ustarp2
            fwc(i,j) = 2.0*ustarp2/wavevel(i,j)**2
            rat = ustarm/ustarcw
            zarough(i,j) = wthick2*(zrough(i,j)/wthick2)**rat

         ENDIF
            
      ELSE
         taumean(i,j) = bdragcoefatc(i,j)*bcuratc(i,j)**2
         taumax(i,j) = taumean(i,j)
         tauwav(i,j) = 0.0
         fwc(i,j) = 0.0
         wthick(i,j) = 0.0
         zarough(i,j) = zrough(i,j)

      ENDIF
         
   ENDDO i_610
   ENDDO j_610

ENDIF
   
!
!7. Deallocate
!-------------
!

DEALLOCATE (phicur,wcanglecos)

1000 CALL log_timer_out(npcc,itm_waves)


RETURN

END SUBROUTINE bottom_stress_waves
