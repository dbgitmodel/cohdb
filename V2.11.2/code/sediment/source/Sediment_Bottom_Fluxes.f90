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
! *Sediment_Bottom_Fluxes.f90 Near bed boundary conditions and shear stresses
!                             for the sediment transport module
!
! Author - Boudewijn Decrop, Alexander Breugem (IMDC) and Patrick Luyten
!
! Version - @(COHERENS)Sediment_Bottom_Fluxes.f90  V2.11.2
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description -
!
! Routines - critical_shear_stress, sediment_roughness
!
!************************************************************************
!

!========================================================================

SUBROUTINE bottom_stress_sed
!************************************************************************
!
! *bottom_stress_sed* Bottom stress arrays for sediment module
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Sediment_Bottom_Fluxes.f90  V2.11.1
!
! Description -
!
! Reference -
!
! Calling program - sediment_equation
!
! External calls - bottom_stress_waves, sediment_roughness, wave_boundary_layer
!
! Module calls - error_alloc, exchange_mod, Carr_at_U, Carr_at_V
!
!************************************************************************
!
USE currents  
USE depths
USE fluxes  
USE grid
USE gridpars  
USE iopars
USE physpars
USE sedarrays
USE sedids
USE sedswitches
USE switches
USE syspars
USE timepars
USE array_interp, ONLY: Carr_at_U, Carr_at_V
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d, bcuratc2, bcuratu2, &
                                         & bcuratv2, bstresatc_old


procname(pglev+1) = 'bottom_stress_sed'
CALL log_timer_in()

!
!1. Same arrays as for hydrodynamic model
!----------------------------------------
!

IF (iopt_sed_bstres.EQ.0) THEN

   bdragcoefatc_sed = bdragcoefatc
   bdragcoefatu_sed = bdragcoefatu
   bdragcoefatv_sed = bdragcoefatv
   bstresatc_sed = bstresatc
   bstresatu_sed = bstresatu
   bstresatv_sed = bstresatv
   ubstresatu_sed = ubstresatu
   vbstresatv_sed = vbstresatv
   zroughatc_sed = zroughatc
   zroughatu_sed = zroughatu
   zroughatv_sed = zroughatv
   IF (iopt_waves.GT.0) THEN
      bstresatc_wav_sed = bstresatc_wav
      fwave_sed = fwave
      wavethickatc_sed = wavethickatc
      zaroughatc_sed = zaroughatc
   ENDIF

   GOTO 1000
   
ENDIF
   
!
!2. Allocate and initialise
!--------------------------
!

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)
maskatu = node2du(1:ncloc,1:nrloc).GT.0

ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)
maskatv = node2dv(1:ncloc,1:nrloc).GT.0

ALLOCATE (array2d(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d',2,(/ncloc,nrloc/),kndrtype)
array2d = 0.0

ALLOCATE (bcuratc2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcuratc2',2,(/ncloc,nrloc/),kndrtype)
bcuratc2 = bcuratc*bcuratc

ALLOCATE (bcuratu2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcuratu2',2,(/ncloc,nrloc/),kndrtype)
bcuratu2 = bcuratu*bcuratu

ALLOCATE (bcuratv2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcuratv2',2,(/ncloc,nrloc/),kndrtype)
bcuratv2 = bcuratv*bcuratv

IF (nt.EQ.0) THEN
   ALLOCATE (bstresatc_old(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bstresatc_old',2,(/ncloc,nrloc/),kndrtype)
   bstresatc_old = 0.0
ENDIF

!
!3. Roughness length
!-------------------
!

IF (nt.EQ.0.OR.iopt_sed_rough.GT.2) CALL sediment_roughness
   
!
!4. Bottom drag coefficient
!--------------------------
!
!4.2.1 C-nodes
!-------------
!         

IF (iopt_bstres_nodim.EQ.2) THEN
   WHERE (maskatc_int)
      array2d = MAX(deptotatc(1:ncloc,1:nrloc)/&
                 & (enap*zroughatc_sed(1:ncloc,1:nrloc)),zbtoz0lim)
   END WHERE
ELSEIF (iopt_bstres_nodim.EQ.3) THEN
   WHERE (maskatc_int)
      array2d = MAX(0.5*delzatc(1:ncloc,1:nrloc,1)/zroughatc_sed,zbtoz0lim)
   END WHERE
ENDIF
WHERE (maskatc_int)
   bdragcoefatc_sed = (ckar/LOG(array2d))**2
END WHERE

!         
!4.2.2 U-nodes
!-------------         
!

IF (iopt_waves.EQ.0) THEN

   IF (iopt_bstres_nodim.EQ.2) THEN
      WHERE (nodeatu(1:ncloc,1:nrloc,1).GT.1)
         array2d = MAX(deptotatu(1:ncloc,1:nrloc)/(enap*zroughatu_sed),&
                     & zbtoz0lim)
      END WHERE
   ELSEIF (iopt_bstres_nodim.EQ.3) THEN
      WHERE (nodeatu(1:ncloc,1:nrloc,1).GT.1)
         array2d = MAX(0.5*delzatu(1:ncloc,1:nrloc,1)/zroughatu_sed,zbtoz0lim)
      END WHERE
   ENDIF
   WHERE (nodeatu(1:ncloc,1:nrloc,1).GT.1)
      bdragcoefatu_sed = (ckar/LOG(array2d))**2
   ELSEWHERE
      bdragcoefatu_sed = 0.0
   END WHERE

!
!4.2.3 V-nodes
!-----------
!

   IF (iopt_bstres_nodim.EQ.2) THEN
      WHERE (nodeatv(1:ncloc,1:nrloc,1).GT.1)
         array2d = MAX(deptotatv(1:ncloc,1:nrloc)/(enap*zroughatv_sed),&
                     & zbtoz0lim)
      END WHERE
   ELSEIF (iopt_bstres_nodim.EQ.3) THEN
      WHERE (nodeatv(1:ncloc,1:nrloc,1).GT.1)
         array2d = MAX(0.5*delzatv(1:ncloc,1:nrloc,1)/zroughatv_sed,zbtoz0lim)
      END WHERE
   ENDIF
   WHERE (nodeatv(1:ncloc,1:nrloc,1).GT.1)
      bdragcoefatv_sed = (ckar/LOG(array2d))**2
   ELSEWHERE
      bdragcoefatv_sed = 0.0
   END WHERE

ENDIF
   
!
!5. Bottom stress
!----------------
!
!5.1. C-nodes
!------------
!
!5.1.1 Without waves
!-------------------
!

IF (iopt_waves.EQ.0) THEN

   WHERE (maskatc_int)
      bstresatc_sed(1:ncloc,1:nrloc) = bdragcoefatc_sed*bcuratc2
   END WHERE

!
!5.1.2 With waves
!----------------
!

ELSE

!  ---initialise at initial time   
   IF (nt.EQ.0) THEN
      WHERE (maskatc_int)
         bstresatc_old = bdragcoefatc_sed(1:ncloc,1:nrloc)*bcuratc2 
      END WHERE
   ENDIF
         
!  ---wave enhanced stress
   CALL bottom_stress_waves(bstresatc_old,zroughatc_sed,bstresatc_mean_sed,&
                          & bstresatc_sed(1:ncloc,1:nrloc),&
                          & bstresatc_wav_sed(1:ncloc,1:nrloc),&
                          & fwave_sed(1:ncloc,1:nrloc),wavethickatc_sed,&
                          & zaroughatc_sed)
   
!  ---reset old value    
   WHERE (maskatc_int)
      bstresatc_old = bstresatc(1:ncloc,1:nrloc)
   END WHERE
   
!  ---drag coefficient
   WHERE (bcuratc2.GT.0.0)
      bdragcoefatc_sed = bstresatc_sed(1:ncloc,1:nrloc)/bcuratc2
   ELSEWHERE
      bdragcoefatc_sed = 0.0
   END WHERE
   
ENDIF

!
!5.1.3 Exchange halo sections
!----------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(bstresatc_sed,(/0,0/),(/1,0,1,0/),iarr_bstresatc_sed)
   CALL exchange_mod(bstresatc_wav_sed,(/0,0/),(/1,0,1,0/),&
                   & iarr_bstresatc_wav_sed)
   CALL exchange_mod(fwave_sed,(/0,0/),(/1,0,1,0/),iarr_fwave_sed)
ENDIF

!
!5.2 U-nodes
!-----------
!
!---bottom stress
CALL Carr_at_U(bstresatc_sed(:,1:nrloc),bstresatu_sed,1,1,(/0,1,nz/),&
            & (/ncloc,nrloc,nz/),1,iarr_bstresatc_sed,.TRUE.)

!---drag coefficient
IF (iopt_waves.GT.0) THEN
   WHERE (bcuratu2.GT.0)
      bdragcoefatu_sed = bstresatu_sed/bcuratu2
   ELSEWHERE
      bdragcoefatu_sed = 0.0
   END WHERE
ENDIF

!---X-component of bottom stress
IF (iopt_bstres_nodim.EQ.2) THEN
   WHERE (maskatu)
      ubstresatu_sed(1:ncloc,:) = bdragcoefatu*bcuratu*umvel(1:ncloc,1:nrloc)
   END WHERE
ELSEIF (iopt_bstres_nodim.EQ.3) THEN
   WHERE (maskatu)
      array2d = bdragcoefatu_sed*bcuratu
      ubstresatu_sed(1:ncloc,:) = bdragcoefatu*bcuratu*uvel(1:ncloc,1:nrloc,1)
   END WHERE
ENDIF

!
!5.3 V-nodes
!-----------
!
!---bottom stress
CALL Carr_at_V(bstresatc_sed(1:ncloc,:),bstresatv_sed,1,1,(/1,0,nz/),&
            & (/ncloc,nrloc,nz/),1,iarr_bstresatc_sed,.TRUE.)

!---drag coefficient
IF (iopt_waves.GT.0) THEN
   WHERE (bcuratv2.GT.0)
      bdragcoefatv_sed = bstresatv_sed/bcuratv2
   ELSEWHERE
      bdragcoefatv_sed = 0.0
   END WHERE
ENDIF

!--- Y-component of bottom stress
IF (iopt_bstres_nodim.EQ.2) THEN
   WHERE (maskatv)
      vbstresatv_sed(:,1:nrloc) = bdragcoefatv*bcuratv*vmvel(1:ncloc,1:nrloc)
   END WHERE
ELSEIF (iopt_bstres_nodim.EQ.3) THEN
   WHERE (maskatv)
      vbstresatv(:,1:nrloc) = bdragcoefatv*bcuratv*vvel(1:ncloc,1:nrloc,1)
   END WHERE
ENDIF

!
!5.3 Exchange halo sections
!--------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(ubstresatu_sed,(/1,1/),(/0,1,0,0/),iarr_ubstresatu_sed)
   CALL exchange_mod(vbstresatv_sed,(/1,1/),(/0,0,0,1/),iarr_vbstresatv_sed)
ENDIF

!
!6. Deallocate
!-------------
!

DEALLOCATE (maskatu,maskatv)
DEALLOCATE (array2d,bcuratc2,bcuratu2,bcuratv2)
IF (cold_start.OR.nt.EQ.nstep) DEALLOCATE (bstresatc_old)

1000 CALL log_timer_out()


RETURN

END SUBROUTINE bottom_stress_sed
  
!========================================================================

SUBROUTINE critical_shear_stress
!************************************************************************
!
! *critical_shear_stress* critical shear stress for different sediment fractions
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Bottom_Fluxes.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - sediment_equation
!
! External calls - 
!
! Module calls - error_alloc, exchange_mod, Uarr_at_C, Varr_at_C
!
!************************************************************************
!
USE density
USE diffusion
USE grid
USE gridpars
USE iopars
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE switches
USE syspars
USE array_interp, ONLY: Uarr_at_C, Varr_at_C
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: f, ff
REAL, PARAMETER :: log19 = LOG10(19.0), onethird = 1.0/3.0
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, deltas, dstar, hiding,&
                                         & ubstresatc, vbstresatc 

procname(pglev+1) = 'critical_shear_stress'
CALL log_timer_in()

!
!1. Allocate/initialise
!----------------------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
IF (iopt_sed_bstres_cr.GT.1) THEN
   ALLOCATE (deltas(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('deltas',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (dstar(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('destar',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (iopt_sed_slope.EQ.1) THEN
   ALLOCATE (ubstresatc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('ubstresatc',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (vbstresatc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('vbstresatc',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (iopt_sed_hiding.EQ.1) THEN
   ALLOCATE (hiding(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('hiding',2,(/ncloc,nrloc/),kndrtype)
ENDIF

bstres_cr= 0.0

!
!2. Critical shear stress
!------------------------
!

f_200: DO f=1,nfp

!  ---work space arrays
   IF (iopt_sed_bstres_cr.GT.1) THEN
      WHERE (maskatc_int)
         deltas = MAX(rhos(f)/dens(1:ncloc,1:nrloc,1)-1.0,0.0)
      ELSEWHERE
         deltas = 0.0
      END WHERE
   ENDIF
   IF (iopt_sed_bstres_cr.EQ.2.OR.iopt_sed_bstres_cr.EQ.3) THEN
      WHERE (maskatc_int)
         dstar = dp(f)*((deltas*gaccatc(1:ncloc,1:nrloc)&
               & /kinvisc(1:ncloc,1:nrloc,1)**2)**onethird)
      END WHERE
   ENDIF

!
!2.1 "Base" value
!----------------
!

   SELECT CASE (iopt_sed_bstres_cr)
!     ---constant
      CASE (1)
         WHERE (maskatc_int)
            bstres_cr(1:ncloc,1:nrloc,f) = bstres_cr_cst(f)
         END WHERE
!     ---Brownlie (1981)
      CASE (2)
         WHERE (maskatc_int)
            array2dc = dstar**(-0.9)
            array2dc = 0.22*array2dc + 0.06*10.0**(-7.7*array2dc)
            bstres_cr(1:ncloc,1:nrloc,f) = array2dc*deltas*&
                                      & gaccatc(1:ncloc,1:nrloc)*dp(f)
         END WHERE
!     ---Soulsby-Whitehouse
      CASE (3)
         WHERE (maskatc_int)
            array2dc = 0.3/(1.0+1.2*dstar) + 0.055*(1.0-EXP(-0.02*dstar))
            bstres_cr(1:ncloc,1:nrloc,f) = array2dc*deltas*&
                                      & gaccatc(1:ncloc,1:nrloc)*dp(f)
         END WHERE

!     ---Wu et al. (2000)
      CASE (4)
         WHERE (maskatc_int)
            bstres_cr(1:ncloc,1:nrloc,f) = 0.03*deltas*&
                                         & gaccatc(1:ncloc,1:nrloc)*dp(f)
         END WHERE
   END SELECT

!
!2.2 Bed slope effects 
!---------------------
!

   IF (iopt_sed_slope.EQ.1) THEN
      CALL Uarr_at_C(ubstresatu_sed,ubstresatc,1,1,(/1,1,nz/),&
                  & (/ncloc+1,nrloc,nz/),1,iarr_ubstresatu_sed,.TRUE.)
      CALL Varr_at_C(vbstresatv_sed,vbstresatc,1,1,(/1,1,nz/),&
                   & (/ncloc,nrloc+1,nz/),1,iarr_vbstresatv_sed,.TRUE.)
      WHERE (maskatc_int.AND.bstresatc_sed(1:ncloc,1:nrloc).GT.0.0)
         array2dc  = (bed_slope_x*ubstresatc + &
                    & bed_slope_y*vbstresatc)/bstresatc_sed(1:ncloc,1:nrloc)
         bstres_cr(1:ncloc,1:nrloc,f) = bstres_cr(1:ncloc,1:nrloc,f)*&
                                      & (1.0+array2dc)
      END WHERE
   ENDIF

!
!2.3 Hiding factor
!-----------------
!
!  ---Wu, wang & Jia (2000)
   IF (iopt_sed_hiding.EQ.1.AND.nf.GT.1) THEN
      array2dc = 0.0
      ff_231: DO ff=1,nfp
         WHERE (maskatc_int)
            array2dc = array2dc + &
                     & bed_fraction(1:ncloc,1:nrloc,1,ff)*dp(ff)/(dp(f)+dp(ff))
         END WHERE
      ENDDO ff_231
      WHERE (maskatc_int.AND.array2dc.GT.1.0E-15 .AND.array2dc.LT.1.0)
         hiding = (1.0/array2dc-1.0)**(-wu_exp)
         bstres_cr(1:ncloc,1:nrloc,f) = hiding*bstres_cr(1:ncloc,1:nrloc,f)
      END WHERE

!  ---Ashida & Michiue (1972)
   ELSEIF (iopt_sed_hiding.EQ.2) THEN
      WHERE (maskatc_int)
         array2dc = d50_bed/dp(f)
         hiding = MERGE(0.8429*array2dc,&
                 & (log19/(log19-LOG10(array2dc)))**2,array2dc.LT.0.38889)
         bstres_cr(1:ncloc,1:nrloc,f) = hiding*bstres_cr(1:ncloc,1:nrloc,f)
     END WHERE
  ENDIF

ENDDO f_200

!
!3. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(bstres_cr,(/0,0,1/),(/1,0,1,0/),iarr_bstres_cr)
ENDIF

!
!4. Deallocate
!-------------
!

DEALLOCATE (array2dc)
IF (iopt_sed_bstres_cr.GT.1) DEALLOCATE (deltas,dstar)
IF (iopt_sed_slope.EQ.1) DEALLOCATE (ubstresatc,vbstresatc)
IF (iopt_sed_hiding.EQ.1) DEALLOCATE (hiding)

CALL log_timer_out()


RETURN

END SUBROUTINE critical_shear_stress

!========================================================================

SUBROUTINE sediment_roughness
!************************************************************************
!
! *sediment_roughness* Bottom roughness length
!
! Author - Alexander Breugem (IMDC) and Patrick Luyten
!
! Version - @(COHERENS)Sediment_Bottom_Fluxes.f90  V2.11.1
!
! Description -
!
! Reference -
!
! Calling program - bottom_stres
!
! External calls -
!
! Module calls - error_alloc, exchange_mod, Carr_at_U, Carr_at_V
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
USE switches
USE array_interp, ONLY: Carr_at_U, Carr_at_V
USE paral_comms, ONLY: exchange_mod 
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'sediment_roughness'
CALL log_timer_in()

!
!1. Bottom roughness
!-------------------
!

IF (iopt_sed_rough.EQ.1) THEN
   WHERE (maskatc_int)
      zroughatc_sed(1:ncloc,1:nrloc) = zrough_sed_cst
   END WHERE
ELSEIF (iopt_sed_rough.EQ.3) THEN
   WHERE (maskatc_int)
      zroughatc_sed(1:ncloc,1:nrloc) = zrough_grain*d50_bed(1:ncloc,1:nrloc)
   END WHERE
ENDIF

!
!2. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(zroughatc_sed,(/0,0/),(/1,0,1,0/),iarr_zroughatc_sed)
ENDIF

!
!3. Roughness at velocity nodes 
!------------------------------
!

CALL Carr_at_U(zroughatc_sed(:,1:nrloc),zroughatu_sed,1,1,(/0,1,nz/),&
     & (/ncloc,nrloc,nz/),1,iarr_zroughatc_sed,.TRUE.)

CALL Carr_at_V(zroughatc_sed(1:ncloc,:),zroughatv_sed,1,1,(/1,0,nz/),&
            & (/ncloc,nrloc,nz/),1,iarr_zroughatc_sed,.TRUE.)

CALL log_timer_out()


RETURN

END SUBROUTINE sediment_roughness
