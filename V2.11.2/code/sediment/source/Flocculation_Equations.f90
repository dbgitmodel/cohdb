!
! Copyright 2013 RBINS-MUMM
!
! Licensed under the EUPL, Version 1.1 or - as soon they will be approved by
! the European Commission - subsequent versions of the EUPL (the "Licence");
! You may not use this work except in compliance with the Licence.
! You may obtain a copy fof the Licence at:
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
! *Flocculation_Equations* COHERENS flocculation model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.1
!
! $Date: 2018-01-31 13:59:45 +0100 (Wed, 31 Jan 2018) $
!
! $Revision: 1085 $
!
! Description - The model contains the equations for calculating suspended
!               sediment transport taken account of flocculation processes
!
! Routines - flocculation_model, floc_arrays, floc_settling_velocity, floc_sources
!
!************************************************************************
!

!========================================================================
  
SUBROUTINE flocculation_model
!************************************************************************
!
! *flocculation_model* Flocculation module
!
! Author - Valerie Duliere and Patrick Luyten (IMDC)
!
! Version - @(COHERENS)Flocculation_Equations.f90  V2.11.1
!
! Description -
!
! Reference -
!
! Calling program - sediment_equation
!
! External calls - beta_factor, define_profobc_spec, floc_arrays,
!                  floc_settling_velocity, floc_sources,
!                  open_boundary_conds_prof, transport_at_C_4d1,
!                  transport_at_C_4d2, update_nest_data_prof,
!                  update_profobc_data
!
! Module calls - error_alloc, exchange_mod
!
!************************************************************************
!
USE diffusion  
USE grid
USE gridpars  
USE iopars
USE nestgrids
USE physpars
USE relaxation
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out  

IMPLICIT NONE

!
!*Local variables
!
LOGICAL, SAVE :: settling
INTEGER :: f, k
INTEGER, SAVE :: maxprofs, noprofs, nofiles
INTEGER, DIMENSION(3) :: ivarids
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: iprofrlx, itypobu, itypobv, noprofsd
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: indexprof, indexvar, iprofobu, &
                                            & iprofobv
INTEGER (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:,:) :: nosecsprof
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: taurat, weight1, weight2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: cnumpobu, cnumpobv, obccnump3d, &
                                           & zeros3d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: bcbot, flocsource, profcnump


procname(pglev+1) = 'flocculation_model'
CALL log_timer_in()

IF (nt.EQ.0) THEN

   nofiles = maxdatafiles(io_sedobc,1)

!  ---variables ids
   IF (iopt_transp_full.EQ.0) ivarids = (/iarr_floc_P, iarr_floc_F,iarr_floc_T/)

!  ---allocate specifier arrays and o.b. profiles
   ALLOCATE (cnumpobu(nobu,nz,nf),STAT=errstat)
   CALL error_alloc('cnumpobu',3,(/nobu,nz,nf/),kndrtype)
   IF (nobu.GT.0) cnumpobu = real_fill
   ALLOCATE (cnumpobv(nobv,nz,nf),STAT=errstat)
   CALL error_alloc('cnumpobv',3,(/nobv,nz,nf/),kndrtype)
   IF (nobv.GT.0) cnumpobv = real_fill
   ALLOCATE (iprofrlx(norlxzones),STAT=errstat)
   CALL error_alloc('iprofrlx',1,(/norlxzones/),kndint)
   IF (norlxzones.GT.0) iprofrlx = 0
   IF (iopt_obc_sed.EQ.1) THEN
      ALLOCATE (itypobu(nobu),STAT=errstat)
      CALL error_alloc('itypobu',1,(/nobu/),kndint)
      ALLOCATE (itypobv(nobv),STAT=errstat)
      CALL error_alloc('itypobv',1,(/nobv/),kndint)
      ALLOCATE (iprofobu(nobu,nf),STAT=errstat)
      CALL error_alloc('iprofobu',2,(/nobu,nf/),kndint)
      ALLOCATE (iprofobv(nobv,nf),STAT=errstat)
      CALL error_alloc('iprofobv',2,(/nobv,nf/),kndint)
      ALLOCATE (noprofsd(2:nofiles),STAT=errstat)
      CALL error_alloc('noprofsd',1,(/nofiles-1/),kndint)
      ALLOCATE (indexprof(nf*(nobu+nobv),2:nofiles),STAT=errstat)
      CALL error_alloc('indexprof',2,(/nf*(nobu+nobv),nofiles-1/),kndint)
      ALLOCATE (indexvar(nf*(nobu+nobv),2:nofiles),STAT=errstat)
      CALL error_alloc('indexvar',2,(/nf*(nobu+nobv),nofiles-1/),kndint)
   ENDIF

!  ---define specifier arrays
   IF (iopt_obc_sed.EQ.1) THEN
      CALL define_profobc_spec(io_sedobc,itypobu,itypobv,iprofobu,iprofobv,&
                             & iprofrlx,noprofsd,indexprof,indexvar,nf,&
                             & nofiles,nobu,nobv)
   ENDIF

!  ---allocate o.b. data arrays
   noprofs = 0; maxprofs = 0
   IF (iopt_obc_sed.EQ.1) THEN
      IF (nofiles.GT.1) noprofs = MAX(MAXVAL(iprofobu),MAXVAL(iprofobv))
        ALLOCATE (obccnump3d(0:noprofs,nz,nf),STAT=errstat)
        CALL error_alloc('obccnump3d',3,(/noprofs+1,nz,nf/),kndrtype)
        obccnump3d = real_fill
   ENDIF
   IF (nofiles.GT.1) THEN
      maxprofs = MAXVAL(noprofsd)
      ALLOCATE (nosecsprof(2:nofiles,2),STAT=errstat)
      CALL error_alloc('nosecsprof',2,(/nofiles-1,2/),kndilong)
      ALLOCATE (profcnump(maxprofs,nz,2:nofiles,2),STAT=errstat)
      CALL error_alloc('profcnump',4,(/maxprofs,nz,nofiles-1,2/),kndrtype)
   ENDIF

!  ---update data at initial time
   IF (nofiles.GT.1) THEN
      CALL update_profobc_data(profcnump,obccnump3d,noprofsd,&
                             & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                             & maxprofs,noprofs,nf,nofiles,nosecsprof,io_sedobc)
   ENDIF

!  ---settling velocity
   settling = iopt_sed_vadv.GT.0
   IF (settling) CALL floc_settling_velocity

!  ---vertical diffusion coefficient
   CALL beta_factor

   f_110: DO f=1,nf
   k_110: DO k=1,nz+1
      WHERE (maskatc_int)
         vdiffcoef_sed(:,:,k,f) = beta_sed(:,:,f)*vdifcoefmom(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_110
   ENDDO f_110

!  ---write at nest locations
   IF (iopt_nests.EQ.1) THEN 
      CALL update_nest_data_prof(cnump,nf,nosednst,instsed,io_sednst)
   ENDIF

   GOTO 1000

ENDIF

!
!2. Allocate and initialise
!--------------------------
!

ALLOCATE (bcbot(ncloc,nrloc,1,nf),STAT=errstat)
CALL error_alloc('bcbot',4,(/ncloc,nrloc,1,nf/),kndrtype)
ALLOCATE (flocsource(ncloc,nrloc,nz,nf),STAT=errstat)
CALL error_alloc('flocsource',4,(/ncloc,nrloc,nz,nf/),kndrtype)
ALLOCATE (taurat(ncloc,nrloc),STAT=errstat)
CALL error_alloc('taurat',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (zeros3d(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('zeros3d',3,(/ncloc,nrloc,nf/),kndrtype)
IF (iopt_scal_depos.EQ.2) THEN
   ALLOCATE (weight1(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('weight1',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (weight2(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('weight2',2,(/ncloc,nrloc/),kndrtype)
ENDIF

bcbot = 0; flocsource = 0.0; zeros3d = 0.0


!
!3. Open boundary conditions
!---------------------------
!
!---update data
IF (nofiles.GT.1) THEN
   CALL update_profobc_data(profcnump,obccnump3d,noprofsd,&
                          & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                          & maxprofs,noprofs,nf,nofiles,nosecsprof,io_sedobc)
ENDIF

!---update o.b. profiles
IF (iopt_obc_sed.EQ.1) THEN
   CALL open_boundary_conds_prof(cnump,cnumpobu,cnumpobv,obccnump3d,&
                               & obccnumpatu,obccnumpatv,itypobu,itypobv,&
                               & iprofobu,iprofobv,noprofs,nf)
ENDIF

!
!4. Source terms
!---------------
!

CALL floc_sources(flocsource)

!
!5. Bottom fluxes
!----------------
!
!5.1 Initialise
!--------------
!

bottom_sed_ero = 0.0; bottom_sed_dep = 0.0

!
!5.2 Erosion rate
!----------------
!

WHERE (maskatc_int)
   taurat = MAX(bstresatc_sed(1:ncloc,1:nrloc)&
               /bstres_cr(1:ncloc,1:nrloc,1)-1.0,0.0)
   bottom_sed_ero(:,:,1)  = (parth_coef/massp(1))*(taurat**parth_exp)
   bcbot(:,:,1,1) = -bottom_sed_ero(:,:,1)
END WHERE

!
!5.3 Explicit part of deposition flux
!------------------------------------
!

IF (settling.AND.iopt_vadv_impl.LT.2) THEN

!  ---first order scheme
   IF (iopt_scal_depos.EQ.1) THEN
      f_531: DO f=1,nf
         bottom_sed_dep(:,:,f) = (1.0-theta_vadv)*&
                                & wfall(1:ncloc,1:nrloc,1,f)*&
                                & cnump(1:ncloc,1:nrloc,1,f)
      ENDDO f_531
!  ---second order scheme      
   ELSEIF (iopt_scal_depos.EQ.2) THEN
      WHERE (maskatc_int)
         weight1 = 0.5*delzatc(1:ncloc,1:nrloc,1)/delzatw(1:ncloc,1:nrloc,2)
         weight2 = 1.0 + weight1
      END WHERE
      f_532: DO f=1,nf
         WHERE (maskatc_int)
            bottom_sed_dep(:,:,f) = (1.0-theta_vadv)*&
                                  & wfall(1:ncloc,1:nrloc,1,f)*&
                                  & (weight2*cnump(1:ncloc,1:nrloc,1,f)-&
                                  &  weight1*cnump(1:ncloc,1:nrloc,2,f))
         END WHERE
      ENDDO f_532
   ENDIF

ENDIF

!
!6. Solve transport equation
!---------------------------
!
   
IF (iopt_transp_full.EQ.0) THEN
   CALL transport_at_C_4d1(cnump,flocsource,wfall,vdiffcoef_sed,nf,&
           & iopt_adv_scal,iopt_sed_vadv,iopt_hdif_scal,cnumpobu,cnumpobv,&
           & iprofrlx,0,1,1,1,zeros3d,bcbot,settling,ivarids)
ELSEIF (iopt_transp_full.EQ.1) THEN
   CALL transport_at_C_4d2(cnump,flocsource,wfall,vdiffcoef_sed,nf,&
           & iopt_adv_scal,iopt_sed_vadv,iopt_hdif_scal,cnumpobu,cnumpobv,&
           & iprofrlx,0,1,1,1,zeros3d,bcbot,settling,iarr_cnump)
ENDIF

!
!7. Add implicit parts to bottom flux
!------------------------------------
!


IF (iopt_vadv_impl.GT.0.AND.settling) THEN
      
!  ---first order scheme
   IF (iopt_scal_depos.EQ.1) THEN
      f_710: DO f=1,nf
         WHERE (maskatc_int)
            bottom_sed_dep(:,:,f) = bottom_sed_dep(:,:,f) + &
                                  & theta_vadv*wfall(1:ncloc,1:nrloc,1,f)*&
                                  & cnump(1:ncloc,1:nrloc,1,f)
         END WHERE
      ENDDO f_710

!  ---second order scheme
   ELSEIF (iopt_scal_depos.EQ.2) THEN
      WHERE (maskatc_int)
         weight1 = 0.5*delzatc(1:ncloc,1:nrloc,1)/delzatw(1:ncloc,1:nrloc,2)
         weight2 = 1.0 + weight1
      END WHERE
      f_720: DO f=1,nf
         WHERE (maskatc_int)
            bottom_sed_dep(:,:,f) = bottom_sed_dep(:,:,f) - &
                                  & theta_vadv*wfall(1:ncloc,1:nrloc,1,f)*&
                                  & (weight2*cnump(1:ncloc,1:nrloc,1,f)-&
                                  &  weight1*cnump(1:ncloc,1:nrloc,2,f))
         END WHERE
      ENDDO f_720
   ENDIF

ENDIF

!
!8. Macrofloc properties and sediment concentrations
!---------------------------------------------------
!

CALL floc_arrays

!
!9. Bottom sediment flux
!-----------------------
!

WHERE (maskatc_int)
   bottom_sed_flux(:,:,1) = volp(1)*(bottom_sed_ero(:,:,1)-&
                                   & bottom_sed_dep(:,:,1))
   bottom_sed_flux(:,:,2) = -floc_vol(:,:,1)*bottom_sed_dep(:,:,2)
   bottom_sed_flux(:,:,3) = -volp(1)*bottom_sed_dep(:,:,3)
END WHERE

!
!10. Settling velocity
!---------------------
!

IF (settling) CALL floc_settling_velocity

!
!11. Vertical diffusion coefficient
!----------------------------------
!
!11.1 Beta factor
!----------------
!

IF (iopt_sed_beta.EQ.3) CALL beta_factor

!
!11.2 Vertical diffusion coefficient
!-----------------------------------
!

f_1121: DO f=1,nf
k_1121: DO k=1,nz+1
   WHERE (maskatc_int)
      vdiffcoef_sed(:,:,k,f) = beta_sed(:,:,f)*vdifcoefmom(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_1121
ENDDO f_1121

!
!12. Write at nest locations
!---------------------------
!

IF (iopt_nests.EQ.1) THEN
   CALL update_nest_data_prof(cnump,nf,nosednst,instsed,io_sednst)
ENDIF

!
!13. Deallocate arrays
!---------------------
!
!---work space
DEALLOCATE (bcbot,flocsource,taurat,zeros3d)
IF (iopt_scal_depos.EQ.2) DEALLOCATE (weight1,weight2)

!---at final step
1000 CONTINUE   
IF (cold_start.OR.nt.EQ.nstep) THEN
   DEALLOCATE (cnumpobu,cnumpobv,iprofrlx)
   IF (iopt_obc_sed.EQ.1) THEN
      DEALLOCATE (itypobu,itypobv,iprofobu,iprofobv,noprofsd,indexprof,&
                & indexvar,obccnump3d)
   ENDIF
   IF (nofiles.GT.1) DEALLOCATE (nosecsprof,profcnump)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE flocculation_model

!======================================================================

SUBROUTINE floc_arrays
!************************************************************************
!
! *floc_arrays* Macrofloc properties and sediment concentrations
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Flocculation_Equations.f90  V2.11.1
!
! Description - 
!
! Calling program - flocculation_model
!
! External calls -
!
! Module calls - exchange_mod
!
!************************************************************************
!
USE grid
USE gridpars  
USE iopars
USE sedarrays
USE sedids
USE sedpars
USE switches
use timepars
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out  

IMPLICIT NONE

!
! Local variables
!
INTEGER :: k
INTEGER, DIMENSION(4) :: nhexch
REAL :: eps = 1.0E-10


procname(pglev+1) = 'floc_arrays'
CALL log_timer_in()

!
!1. Macrofloc properties
!------------------------
!

floc_ncmax = 100.0; floc_ncmin = 2.0

k_110: DO k=1,nz
   WHERE (maskatc_int)
!     ---number of microflocs bound in a macrofloc      
      floc_nc(1:ncloc,1:nrloc,k) = cnump(1:ncloc,1:nrloc,k,3) &
                   & /(cnump(1:ncloc,1:nrloc,k,2)+eps)
      floc_nc(1:ncloc,1:nrloc,k) = MIN(floc_nc(1:ncloc,1:nrloc,k),floc_ncmax)
      floc_nc(1:ncloc,1:nrloc,k) = MAX(floc_nc(1:ncloc,1:nrloc,k),floc_ncmin)
!     ---floc diameter
      floc_dia(:,:,k) = dp(1)*floc_nc(1:ncloc,1:nrloc,k)**(1.0/nfrdim)
!     ---floc volume      
      floc_vol(:,:,k) = volp(1)*floc_nc(1:ncloc,1:nrloc,k)**(3.0/nfrdim)      
   END WHERE
ENDDO k_110

!
!2. Water column sediments
!-------------------------
!

k_210: DO k=1,nz
   WHERE (maskatc_int)
      cvol(1:ncloc,1:nrloc,k,1) = cnump(1:ncloc,1:nrloc,k,1)*volp(1)
      cvol(1:ncloc,1:nrloc,k,2) = cnump(1:ncloc,1:nrloc,k,2)*floc_vol(:,:,k)
      cvol(1:ncloc,1:nrloc,k,3) = cnump(1:ncloc,1:nrloc,k,3)*volp(1)
   END WHERE
ENDDO k_210
   
!
!3. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = nhalo
   CALL exchange_mod(cnump,(/1-nhalo,1-nhalo,1,1/),nhexch,iarr_cnump)
   CALL exchange_mod(cvol,(/1-nhalo,1-nhalo,1,1/),nhexch,iarr_cvol)
   nhexch = 1
   CALL exchange_mod(floc_nc,(/0,0,1/),nhexch,iarr_floc_nc)
ENDIF

!
!4. Total sediment concentration
!-------------------------------
!

k_410: DO k=1,nz
   WHERE (nodeatc(0:ncloc+1,0:nrloc+1).GT.0)
      ctot(:,:,k) = cvol(0:ncloc+1,0:nrloc+1,k,1) + &
                  & cvol(0:ncloc+1,0:nrloc+1,k,2)
   END WHERE
ENDDO k_410

CALL log_timer_out()


RETURN

END SUBROUTINE floc_arrays
  
!========================================================================

SUBROUTINE floc_settling_velocity
!************************************************************************
!
! *floc_settling_velocity* settling velocities in the flocculation model
!
! Author - Valerie Duliere and Patrick Luyten 
!
! Version - @(COHERENS)Flocculation_Equations.f90  V2.11.1
!
! Description -
!
! Reference -
!
! Calling program - flocculation_model
!
! External calls - 
!
! Module calls - Carr_at_W, error_alloc, exchange_mod
!
!************************************************************************
!
USE density
USE depths
USE diffusion  
USE grid
USE gridpars  
USE iopars
USE modids
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE array_interp, ONLY: Carr_at_W
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out  

IMPLICIT NONE

!
! Local variables
!
INTEGER :: f, k
INTEGER, DIMENSION(4) :: lbounds, nhexch 
REAL, ALLOCATABLE, DIMENSION(:,:) ::  array2d
REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  array3dw, deltas, reynolds

   
procname(pglev+1) = 'floc_settling_velocity'
CALL log_timer_in()
   
!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2d(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array3dw(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('array3dw',3,(/ncloc,nrloc,nz+1/),kndrtype)
ALLOCATE (deltas(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('deltas',3,(/ncloc,nrloc,nz+1/),kndrtype)
ALLOCATE (reynolds(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('reynolds',3,(/ncloc,nrloc,nz+1/),kndrtype)

!
!2. Initialise arrays
!--------------------
!
!---density at W-nodes
CALL Carr_at_W(densw(1:ncloc,1:nrloc,:),array3dw(:,:,2:nz),&
             &(/1,1,1/),(/ncloc,nrloc,nz/),1,iarr_densw,.TRUE.)
WHERE (maskatc_int)
   array3dw(:,:,nz+1) = densw(1:ncloc,1:nrloc,nz)
   array3dw(:,:,1) = densw(1:ncloc,1:nrloc,1)
END WHERE

!---density factor
k_210: DO k=1,nz+1
   WHERE (maskatc_int)
      deltas(:,:,k) = MAX(rhos(1)/array3dw(:,:,k)-1.0,0.0)
   END WHERE
ENDDO k_210

!---floc diameter at W-nodes
array3dw = 0.0
CALL Carr_at_W(floc_dia,array3dw(:,:,2:nz),(/1,1,1/),(/ncloc,nrloc,nz/),1,&
              & iarr_floc_dia,.TRUE.)
WHERE (maskatc_int)
   array3dw(:,:,nz+1) = floc_dia(:,:,nz)
   array3dw(:,:,1) = floc_dia(:,:,1)
END WHERE
   
!---Reynolds number
IF (nt.EQ.0) THEN
   reynolds = 0.0
ELSE
   k_220: DO k=1,nz+1
      WHERE (maskatc_int)
         reynolds(:,:,k) = array3dw(:,:,k)*wfall(1:ncloc,1:nrloc,k,2)&
                        & /kinvisc(1:ncloc,1:nrloc,k)
      ELSEWHERE
         reynolds(:,:,k) = 0.0
      END WHERE
   ENDDO k_220
ENDIF

!
!3. Microflocs
!-------------
!

k_310: DO k=1,nz+1
   WHERE (maskatc_int)
      wfall(1:ncloc,1:nrloc,k,1) = deltas(:,:,k)*gaccatc(1:ncloc,1:nrloc)*&
                                 & dp(1)**2/(18.0*kinvisc(1:ncloc,1:nrloc,k))
   END WHERE
ENDDO k_310

!
!4. Macroflocs
!-------------
!

k_410: DO k=1,nz+1
   WHERE (maskatc_int)
      wfall(1:ncloc,1:nrloc,k,2) = wfall(1:ncloc,1:nrloc,k,1)*&
                                & (dp(1)/array3dw(:,:,k))**(1.0-nfrdim)&
                                & /(1.0+0.15*reynolds(:,:,k)**0.687)
   END WHERE
ENDDO k_410

!
!5. Third fraction
!-----------------
!

wfall(1:ncloc,1:nrloc,:,3) = wfall(1:ncloc,1:nrloc,:,2)

!
!6. Upper limit in shallow waters
!--------------------------------
!

IF (iopt_sed_ws_lim.EQ.1) THEN

   f_610: DO f=1,nf
   k_610: DO k=2,nz+1
      WHERE (maskatc_int)
         array2d = gscoordatw(1:ncloc,1:nrloc,k)*deptotatc(1:ncloc,1:nrloc)
         wfall(1:ncloc,1:nrloc,k,f) = MIN(wfall(1:ncloc,1:nrloc,k,f),&
                                        & wslimfac*array2d/delt3d)
      END WHERE
   ENDDO k_610
   ENDDO f_610

ENDIF

!
!7. Exchange halos
!-----------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = (/0,0,1,1/); nhexch = (/1,0,1,0/)
   CALL exchange_mod(wfall,lbounds,nhexch,iarr_wfall)
ENDIF

!
!8. Deallocate arrays
!--------------------
!

DEALLOCATE (array2d,array3dw,deltas,reynolds)

CALL log_timer_out()


RETURN

END SUBROUTINE floc_settling_velocity

!========================================================================

SUBROUTINE floc_sources(source)
!************************************************************************
!
! *floc_sources* sources in the transport equations of the flocculation module 
!
! Author - Valerie Duliere and Patrick Luyten 
!
! Version - @(COHERENS)Flocculation_Equations.f90  V2.11.1
!
! Description -
!
! Reference -
!
! Calling program - flocculation_model
!
! External calls - 
!
! Module calls -
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
USE sedarrays
USE sedpars
USE turbulence
USE array_interp, ONLY: Warr_at_C
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out  

  
IMPLICIT NONE

!
!*Arguments
!
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc,nz,3) :: source

!
! Name      Type Purpose
!------------------------------------------------------------------------------
!*source*   REAL Source terms in the flocculation transport equations [1/m^3/s]
!
!************************************************************************
!
!*Local variables
!
INTEGER :: k
REAL :: twothird = 2.0/3.0, fourthird = 4.0/3.0
REAL, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2, beta_FF, beta_PF, &
                                   & beta_PP
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: kinviscatc, muvisc, shear_rate


procname(pglev+1) = 'floc_sources'
CALL log_timer_in()

!
!1. Allocate
!-----------
!

ALLOCATE(array2d1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE(array2d2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE(beta_FF(ncloc,nrloc),STAT=errstat)
CALL error_alloc('beta_FF',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE(beta_PF(ncloc,nrloc),STAT=errstat)
CALL error_alloc('beta_PF',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE(beta_PP(ncloc,nrloc),STAT=errstat)
CALL error_alloc('beta_PP',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE(kinviscatc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('kinviscatc',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE(muvisc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('muvisc',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE(shear_rate(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('shear_rate',3,(/ncloc,nrloc,nz/),kndrtype)
   
!
!2. Work space arrays
!--------------------
!
!--- kinematic viscosity
CALL Warr_at_C(kinvisc(1:ncloc,1:nrloc,:),kinviscatc,(/1,1,1/),&
            & (/ncloc,nrloc,nz+1/),1,iarr_kinvisc,.TRUE.)

!---aboslute viscosity
k_210: DO k=1,nz
   WHERE (maskatc_int)
      muvisc(:,:,k) = dens(1:ncloc,1:nrloc,k)*kinviscatc(:,:,k)
   END WHERE
ENDDO k_210

!---shear rate   
WHERE (maskatc_int)
   shear_rate(:,:,1) = SQRT(dissip(1:ncloc,1:nrloc,1)/kinviscatc(:,:,1))
   shear_rate(:,:,nz) = SQRT(dissip(1:ncloc,1:nrloc,nz)/kinviscatc(:,:,nz))
END WHERE
k_220: DO k=2,nz-1
   WHERE (maskatc_int)
      shear_rate(:,:,k) = SQRT(SQRT(dissip(1:ncloc,1:nrloc,k)*&
                               & dissip(1:ncloc,1:nrloc,k+1))/kinviscatc(:,:,k))
   END WHERE
ENDDO k_220

!
!3. Aggegration terms
!--------------------
!

k_310: DO k=1,nz
   WHERE (maskatc_int)
      array2d1 = (twothird*kboltz)*(temp(1:ncloc,1:nrloc,k)+celtokelv)&
               & /muvisc(:,:,k)
      array2d2 = dp(1) + floc_dia(:,:,k)
      beta_PP = 4.0*array2d1 + fourthird*shear_rate(:,:,k)*dp(1)**3
      beta_FF = 4.0*array2d1 + fourthird*shear_rate(:,:,k)*floc_dia(:,:,k)**3
      beta_PF = array2d1*(1.0/dp(1)+1.0/floc_dia(:,:,k))*array2d2
      beta_PF = beta_PF + shear_rate(:,:,k)*array2d2**3/6.0
      beta_PF = beta_PF + 0.5*halfpi*array2d2**2*&
              & ABS(wfall(1:ncloc,1:nrloc,k,2)-wfall(1:ncloc,1:nrloc,k,1))
      array2d1 = 0.5*beta_PP*(1.0/(floc_nc(1:ncloc,1:nrloc,k)-1.0))*&
               & cnump(1:ncloc,1:nrloc,k,1)**2
      array2d2 = array2d1*floc_nc(1:ncloc,1:nrloc,k) + &
               & beta_PF*cnump(1:ncloc,1:nrloc,k,1)*&
               & cnump(1:ncloc,1:nrloc,k,2)
      source(:,:,k,1) = -agg_alpha*array2d2
      array2d2 = array2d1 - 0.5*beta_FF*cnump(1:ncloc,1:nrloc,k,2)**2
      source(:,:,k,2) = agg_alpha*array2d2
   END WHERE
ENDDO k_310

!
!4. Breakage terms
!-----------------
!

k_410: DO k=1,nz
   WHERE (maskatc_int)
      array2d1 = brk_es*shear_rate(:,:,k)*(floc_dia(:,:,k)/dp(1)-1.0)**brk_p
      array2d1 = array2d1*(muvisc(1:ncloc,1:nrloc,k)*shear_rate(:,:,k)*&
               & floc_dia(:,:,k)**2/brk_ffy)**brk_q*cnump(1:ncloc,1:nrloc,k,2)
      source(:,:,k,1) = source(:,:,k,1) + &
                      & brk_frac*floc_nc(1:ncloc,1:nrloc,k)*array2d1
      source(:,:,k,2) = source(:,:,k,2) + array2d1
      source(:,:,k,3) = -source(:,:,k,1)
   END WHERE
ENDDO k_410

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (array2d1,array2d2,beta_FF,beta_PF,beta_PP,kinviscatc,muvisc,&
          & shear_rate)

CALL log_timer_out()


RETURN

END SUBROUTINE floc_sources
