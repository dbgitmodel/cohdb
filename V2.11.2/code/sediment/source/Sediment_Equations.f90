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
! *Sediment_Equations* COHERENS sediment model
!
! Author - Alexander Breugem and Boudewijn Decrop (IMDC), Patrick Luyten
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.2
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - The model contains the equations for calculating suspended
!               sediment transport and bed load, including parametrizations of
!               important processes such as settling velocities
!
! Routines - ackerswhite_params, bartnicki_filter, bed_slope_arrays,
!            beta_factor, diff_coef_waves, equilibrium_concentration,
!            median_particle_diameter, reference_concentration,
!            sediment_advdiff, sediment_bedload, sediment_equation,
!            sediment_suspendedload, sediment_totalload,
!            sediment_settling_velocity, thetastar_engelund_hansen
!
!************************************************************************
!

!========================================================================

SUBROUTINE ackerswhite_params(dstar,maskvals,nw,mw,aw,cw)
!************************************************************************
!
! *ackerswhite_params* 
!
! Author - Boudewijn Decrop (IMDC)
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.5
!
! Description - parameters used in the  Ackers-White (1973) formula
!
! Calling program - sediment_suspendedload, sediment_totalload
!
! Module calls - 
!
!************************************************************************
!
USE gridpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
LOGICAL, INTENT(IN), DIMENSION(ncloc,nrloc) :: maskvals
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc) :: dstar
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: nw, mw, aw, cw

!
! Name            Type    Purpose
!------------------------------------------------------------------------------
!*cnode*          CHAR    Grid location 'C', 'U', 'V'
!*dstar*          REAL    Dimensionless particle diameter
!*maskvals*       LOGICAL Mask to include land points
!*nw, mw, aw, cw* REAL    Parameters for Ackers-White formula
!
!************************************************************************
!

procname(pglev+1) = 'ackerswhite_params'
CALL log_timer_in()


WHERE (maskvals)
   dstar = MIN(60.0,MAX(1.0,dstar))
   nw = 1.0 - 0.243*LOG(dstar)
   mw = 1.67 + 6.83/dstar
   aw = 0.14 + 0.23/SQRT(dstar)
   cw = 2.79*LOG(dstar) - 0.426*LOG(dstar)**2 - 7.97
   cw = EXP(cw)
ELSEWHERE
   nw = 0.0; mw = 0.0; aw = 0.0; cw = 0.0
END WHERE

CALL log_timer_out()


RETURN

END SUBROUTINE ackerswhite_params

!======================================================================

SUBROUTINE bartnicki_filter(psic,novars,ivarid)
!******************************************************************************
!
! *bartnicki_filter* Apply the bartnicki filter to eliminate negative
!                    concentrations
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - sediment_advdiff
!
! External calls - 
!
! Module calls - combine_mod, distribute_mod, error_abort, error_alloc,
!                inquire_var
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE paralpars
USE sedpars
USE syspars
USE error_routines, ONLY: error_abort, error_alloc
USE modvars_routines, ONLY: inquire_var
USE paral_comms, ONLY: combine_mod, distribute_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE
!
!* Arguments
!
INTEGER, INTENT(IN):: ivarid, novars
REAL, INTENT(INOUT), DIMENSION(1-nhalo:ncloc+nhalo,&
                             & 1-nhalo:nrloc+nhalo,nz,novars) :: psic

! Name     Type               Purpose
!------------------------------------------------------------------------------
!*psic*    REAL       C-node quantity to be updated                       [psic]
!*novars*  INTEGER    Size of variable dimension
!*ivarid*  INTEGER    Variable key id_of psic
!
!*******************************************************************************
!
!* Local variables
!
LOGICAL :: isnotfinished
CHARACTER (LEN=lenname):: f90_name
REAL :: mass, negmass, posvol 
INTEGER :: icount, ivar, k, nrnegcells
INTEGER, DIMENSION(4) :: nhdist
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: psiglb, volume, volumeloc


procname(pglev+1) = 'bartnicki_filter'
CALL log_timer_in()

!
!1. Initialise
!-------------
!
!---allocate arrays
ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype) 
ALLOCATE (psiglb(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo,nz),STAT=errstat)
CALL error_alloc('psiglb',3,(/nc+2*nhalo,nr+2*nhalo,nz/),kndrtype)
ALLOCATE (volume(nc,nr,nz),STAT=errstat)
CALL error_alloc('volume',3,(/nc,nr,nz/),kndrtype)
ALLOCATE (volumeloc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('volumeloc',3,(/ncloc,nrloc,nz/),kndrtype)

!---cell volume
array2dc = delxatc(1:ncloc,1:nrloc)*delyatc(1:ncloc,1:nrloc)
k_110: DO k=1,nz
   WHERE (maskatc_int)
      volumeloc(:,:,k) = array2dc*delzatc(1:ncloc,1:nrloc,k)
   ELSEWHERE
      volumeloc(:,:,k) = 0.0
   END WHERE
ENDDO k_110
CALL combine_mod(volume,volumeloc,(/1,1,1/),0,0.0,commall=.TRUE.)

!---variable name
CALL inquire_var(ivarid,f90_name=f90_name)

!
!2. Apply filter
!---------------
!

ivar_200: DO ivar=1,novars

!
!2.1 Concentrations on global domain
!-----------------------------------
!

   CALL combine_mod(psiglb,psic(:,:,:,ivar),(/1-nhalo,1-nhalo,1/),ivarid,0.0,&
                  & commall=.TRUE.)

!
!2.2 Initialise counter
!----------------------
!

   isnotfinished =.TRUE.; icount = 1

!
!2.3  Iteration loop
!-------------------
!

   DO WHILE (isnotfinished)

!
!2.3.1 Computation of negative masses
!----------------------------------
!

      negmass = SUM(volume*psiglb(1:nc,1:nr,:),MASK=psiglb(1:nc,1:nr,:).LT.0.0)
      posvol = SUM(volume,MASK=psiglb(1:nc,1:nr,:).GT.0.0)
      nrnegcells = COUNT(psiglb(1:nc,1:nr,:).LT.0.0)

!
!2.3.2 Adapt concentrations
!------------------------
!

      IF (nrnegcells.GT.0.AND.icount.LE.maxitbartnicki) THEN
         IF (icount.EQ.1.AND.loglev1.GT.0) THEN
            WRITE(iolog,'(A)') REPEAT(' ',pglev-1)//&
           & 'Removing negative concentrations for variable: '//TRIM(f90_name)
         ENDIF
         mass = -negmass/posvol
!        ---change location with a negative mass
         WHERE (psiglb(1:nc,1:nr,:).LT.0.0)
            psiglb(1:nc,1:nr,:) = 0.0
         END WHERE
!        ---change location with a positive mass
         WHERE (psiglb(1:nc,1:nr,:).GT.0.0)
            psiglb(1:nc,1:nr,:) = psiglb(1:nc,1:nr,:) + mass
         END WHERE
         icount = icount + 1
      ELSE
         isnotfinished = .FALSE.
         IF (master.AND.icount.EQ.maxitbartnicki.AND.loglev1.GT.0) THEN
            nerrs = 1
            IF (errchk) THEN
               WRITE (ioerr,'(A)') TRIM(f90_name)//&
                                & ': Bartnicki filter is not converging'
            ENDIF
            CALL error_abort('bartnicki_filter',ierrno_runval)
         ENDIF
      ENDIF

   ENDDO

!
!2.4 Re-distribute global array
!-------------------------------
!

   nhdist = nhalo
   CALL distribute_mod(psiglb,psic(:,:,:,ivar),(/1-nhalo,1-nhalo,1/),nhdist,&
                     & ivarid,0.0)

END DO ivar_200

!
!3. Deallocate
!-------------
!

DEALLOCATE (array2dc,psiglb,volume,volumeloc)

CALL log_timer_out()


RETURN

END SUBROUTINE bartnicki_filter

!=======================================================================

SUBROUTINE bed_slope_arrays
!***********************************************************************

! *bed_slope_arrays* bed slopes at different grid nodes
!
! Author - Alexander Breugem (IMDC) and Patrick Luyten
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.5
!
! Description -
!
! Reference -
!
! Calling program - avalanching, sediment_equation
!
! External calls - 
!
! Module calls - 
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE sedarrays
USE sedswitches
USE time_routines,ONLY: log_timer_in, log_timer_out 

IMPLICIT NONE
!
!*Local variables
!
LOGICAL, DIMENSION(ncloc,nrloc) :: mask


procname(pglev+1) = 'bed_slope_arrays'
CALL log_timer_in()

!
!1. C-nodes
!----------
!
!1.1 X-direction
!---------------
!

mask = seapoint(0:ncloc-1,1:nrloc).AND.seapoint(1:ncloc,1:nrloc).AND.&
     & seapoint(2:ncloc+1,1:nrloc)
WHERE (mask) 
       bed_slope_x = (depmeanatu(1:ncloc,1:nrloc) - &
                    & depmeanatu(2:ncloc+1,1:nrloc))/delxatc(1:ncloc,1:nrloc)
ELSEWHERE
       bed_slope_x = 0.0
END WHERE

!
!1.2 Y-direction
!---------------
!

mask = seapoint(1:ncloc,0:nrloc-1).AND.seapoint(1:ncloc,1:nrloc).AND.&
     & seapoint(1:ncloc,2:nrloc+1)
WHERE (mask) 
       bed_slope_y = (depmeanatv(1:ncloc,1:nrloc) - &
                    & depmeanatv(1:ncloc,2:nrloc+1))/delyatc(1:ncloc,1:nrloc)
ELSEWHERE
       bed_slope_y = 0.0
END WHERE

!
!2. U-nodes
!----------
!
!2.1 X-direction
!---------------
!

mask = seapoint(0:ncloc-1,1:nrloc).AND.seapoint(1:ncloc,1:nrloc)

WHERE (mask)
   bed_slope_x_atu = (depmeanatc(0:ncloc-1,1:nrloc)-&
                    & depmeanatc(1:ncloc,1:nrloc))/delxatu(1:ncloc,1:nrloc)
END WHERE

!
!2.2 Y-direction
!---------------
!

mask = seapoint(1:ncloc,0:nrloc-1).AND.seapoint(1:ncloc,1:nrloc).AND.&
     & seapoint(1:ncloc,2:nrloc+1).AND.seapoint(0:ncloc-1,0:nrloc-1).AND.&
     & seapoint(0:ncloc-1,1:nrloc).AND.seapoint(0:ncloc-1,2:nrloc+1)
WHERE (mask)
    bed_slope_y_atu = (depmeanatuv(1:ncloc,1:nrloc)-&
        & depmeanatuv(1:ncloc,2:nrloc+1))&
        & /delyatu(1:ncloc,1:nrloc)
ELSEWHERE
    bed_slope_y_atu = 0.0
END WHERE
     
!
!3. V-nodes
!----------
!
!3.1 X-direction
!---------------
!

mask = seapoint(0:ncloc-1,1:nrloc).AND.seapoint(1:ncloc,1:nrloc).AND.&
     & seapoint(2:ncloc+1,1:nrloc).AND.seapoint(0:ncloc-1,0:nrloc-1).AND.&
     & seapoint(1:ncloc,0:nrloc-1).AND.seapoint(2:ncloc+1,0:nrloc-1)
WHERE (mask)
    bed_slope_x_atv = (depmeanatuv(1:ncloc,1:nrloc)-&
        & depmeanatuv(2:ncloc+1,1:nrloc))&
        & /delxatv(1:ncloc,1:nrloc)
ELSEWHERE
    bed_slope_x_atv = 0.0
END WHERE
     
!
!3.2 Y-direction
!---------------
!

mask = seapoint(1:ncloc,0:nrloc-1).AND.seapoint(1:ncloc,1:nrloc)
WHERE (mask)
   bed_slope_y_atv = (depmeanatc(1:ncloc,0:nrloc-1)-&
                    & depmeanatc(1:ncloc,1:nrloc))/delyatv(1:ncloc,1:nrloc)
ELSEWHERE
   bed_slope_y_atv = 0.0
END WHERE

CALL log_timer_out()


RETURN

END SUBROUTINE bed_slope_arrays

!=======================================================================

SUBROUTINE beta_factor
!************************************************************************
!
! *beta_factor* Calculates the beta factor for the vertical sediment diffusion
!               coefficient
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.1
!
! Description -
!
! Calling program - sediment_advdiff
!
! Module calls -
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE sedarrays
USE sedpars
USE sedswitches
USE time_routines,ONLY: log_timer_in, log_timer_out 

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: f


procname(pglev+1) = 'beta_factor'
CALL log_timer_in()

f_100: DO f=1,nf
   SELECT CASE (iopt_sed_beta)
      CASE (1)
         WHERE (maskatc_int) beta_sed(:,:,f) = 1.0
      CASE (2)
         WHERE (maskatc_int) beta_sed(:,:,f) = beta_sed_cst
      CASE (3)
         WHERE (bstresatc_sed(1:ncloc,1:nrloc).GT.0.0)
            beta_sed(:,:,f) = 1.0+2.0*&
             & ((wfall(1:ncloc,1:nrloc,nz,f)**2)/bstresatc_sed(1:ncloc,1:nrloc))
            beta_sed(:,:,f) = MIN(beta_sed(:,:,f),beta_sed_max)
            beta_sed(:,:,f) = MAX(beta_sed(:,:,f),beta_sed_min)
         ELSEWHERE
            beta_sed(:,:,f) = 0.0
         END WHERE
   END SELECT
ENDDO f_100

CALL log_timer_out()


RETURN

END SUBROUTINE beta_factor

!========================================================================

SUBROUTINE diff_coef_waves
!************************************************************************
!
! *diff_coef_waves* Calculates the eddy diffusivity due to surface waves
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.2
!
! Description -
!
! Calling program - sediment_advdiff
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
USE sedarrays
USE sedpars
USE syspars
USE wavevars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER ::  f, k
REAL, PARAMETER :: diffmax = 0.05, period_min = 0.01
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, betaw, breakerpar, deltas,&
                                         & diffbed, difftop, vdiffcoef_wav


procname(pglev+1) = 'diff_coef_waves'
CALL log_timer_in()

!
!1. Allocate
!-----------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (betaw(ncloc,nrloc),STAT=errstat)
CALL error_alloc('betaw',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (breakerpar(ncloc,nrloc),STAT=errstat)
CALL error_alloc('breakerpar',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (deltas(ncloc,nrloc),STAT=errstat)
CALL error_alloc('deltas',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (diffbed(ncloc,nrloc),STAT=errstat)
CALL error_alloc('diffbed',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (difftop(ncloc,nrloc),STAT=errstat)
CALL error_alloc('difftop',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (vdiffcoef_wav(ncloc,nrloc),STAT=errstat)
CALL error_alloc('vdiffcoef_wav',2,(/ncloc,nrloc/),kndrtype)

!
!2. Wave parameters
!------------------
!
!---wave breaking coefficient
WHERE (maskatc_int)
   array2dc = MAX(waveheight(1:ncloc,1:nrloc)/deptotatc(1:ncloc,1:nrloc)-0.4,&
                & 0.0)
   breakerpar = 1.0 + array2dc**0.5
END WHERE

!---wave diffusion coefficient in top layer
WHERE (maskatc_int)
   difftop = MIN(0.035*breakerpar*deptotatc(1:ncloc,1:nrloc)*&
               & waveheight(1:ncloc,1:nrloc)/&
               & MAX(waveperiod(1:ncloc,1:nrloc),period_min),diffmax)
   deltas = 0.72*breakerpar*waveexcurs(1:ncloc,1:nrloc)**0.75 *&
          & (30.0*zroughatc_sed(1:ncloc,1:nrloc))**0.25
   array2dc = deltas/deptotatc(1:ncloc,1:nrloc)
END WHERE

!
!3. Current-wave diffusion coefficient
!-------------------------------------
!

f_310: DO f=1,nf
!  ---bottom layer
   WHERE (maskatc_int)
      betaw = 1.0 + 2.0*wfall(1:ncloc,1:nrloc,1,1)**2/bstresatc_wav_sed
      diffbed = 0.018*breakerpar*betaw*deltas*wavevel(1:ncloc,1:nrloc)
   END WHERE
!  ---at all layers
   k_311: DO k=1,nz+1
      WHERE (maskatc_int)
         vdiffcoef_wav = MERGE(difftop,&
            & MERGE(diffbed,diffbed+(difftop-diffbed)*&
                 & (gscoordatw(1:ncloc,1:nrloc,k)-array2dc)/(0.5-array2dc),&
                 & gscoordatw(1:ncloc,1:nrloc,k).LT.array2dc),&
                 & gscoordatw(1:ncloc,1:nrloc,k).GT.0.5)
         vdiffcoef_sed(1:ncloc,1:nrloc,k,f) = &
              & SQRT(vdiffcoef_sed(1:ncloc,1:nrloc,k,f)**2+vdiffcoef_wav**2)
     END WHERE
  ENDDO k_311
ENDDO f_310

!
!4. Deallocate
!-------------
!

DEALLOCATE (array2dc,betaw,breakerpar,deltas,diffbed,difftop,vdiffcoef_wav)

CALL log_timer_out()


RETURN

END SUBROUTINE diff_coef_waves

!=======================================================================

SUBROUTINE equilibrium_concentration
!***********************************************************************
!
! *equilibrium_concentration* equilibrium concentration and time scale
!                             (if needed)
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - sediment_advdiff
!
! External calls - sediment_suspendedload
!
! Module calls - error_alloc, gauss_quad, vector_mag_arr_atc
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE physpars
USE sedarrays
USE sedpars
USE sedswitches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: gauss_quad, vector_mag_arr_atc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
REAL, PARAMETER :: hdcurmin = 0.001, p1 = 67.342, p2 = -321.0, p3 = 614.14, &
                 & p4 = -592.75, p5 = 292.15, p6 = -62.141, p7 = 3.667, &
                 & p8 = -2.2571, p9 = 0.97978, rousemax_c = 2.5, &
                 & rousemax_d = 10.0, rousemax_t = 1.14
INTEGER :: f, n
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:) :: nodes_depth_avg, &
                                                       & weights_depth_avg
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: weights
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc1, array2dc2, array2dc3, &
                                         & conc, rouse
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: concprof, elev_depth_avg


procname(pglev+1) = 'equilibrium_concentration'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!
!1.1 At initial time
!-------------------
!   

IF (nt.EQ.0.AND.iopt_sed_nodim.EQ.2.AND.iopt_sed_bbc_eq.EQ.1) THEN
   ALLOCATE (weights_depth_avg(nrquad_sed),STAT=errstat)
   CALL error_alloc('weights_depth_avg',1,(/nrquad_sed/),kndrtype)
   weights_depth_avg = 0.0
   ALLOCATE (nodes_depth_avg(nrquad_sed),STAT=errstat)
   CALL error_alloc('nodes_depth_avg',1,(/nrquad_sed/),kndrtype)
   nodes_depth_avg = 0.0
   CALL gauss_quad(nrquad_sed,weights_depth_avg,nodes_depth_avg)
   ALLOCATE (weights(ncloc,nrloc,nrquad_sed),STAT=errstat)
   CALL error_alloc('weights',3,(/ncloc,nrloc,nrquad_sed/),kndrtype)
   weights = 0.0
   n_110: DO n=1,nrquad_sed
      weights(:,:,n) = weights_depth_avg(n)
   ENDDO n_110
ENDIF

!
!1.2 At each call
!----------------
   
   
ALLOCATE (array2dc1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (rouse(ncloc,nrloc),STAT=errstat)
CALL error_alloc('rouse',2,(/ncloc,nrloc/),kndrtype)
IF (iopt_sed_nodim.EQ.2.AND.iopt_sed_bbc_eq.EQ.1) THEN
   ALLOCATE (array2dc2(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc2',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (array2dc3(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc3',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (conc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('conc',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (concprof(ncloc,nrloc,nrquad_sed),STAT=errstat)
   CALL error_alloc('concprof',3,(/ncloc,nrloc,nrquad_sed/),kndrtype)
   concprof = 0.0
   ALLOCATE(elev_depth_avg(ncloc,nrloc,nrquad_sed),STAT=errstat)
   CALL error_alloc('elev_depth_avg',3,(/ncloc,nrloc,nrquad_sed/),kndrtype)
   elev_depth_avg = 0.0
ENDIF

!
!2. Equilibrium concentration
!----------------------------
!
!2.1 2-D case
!------------
!
   
IF (iopt_sed_nodim.EQ.2) THEN

   IF (iopt_sed_bbc_eq.EQ.1) THEN
      
!   
!2.1.1 From depth Rouse profile
!------------------------------
!

      WHERE (maskatc_int)
         array2dc1 = height_c(1:ncloc,1:nrloc,1)/deptotatc(1:ncloc,1:nrloc)
         array2dc2 = 0.5*LOG(array2dc1)
         rouse = rousemax_d
      END WHERE

      f_211: DO f=1,nf
!        ---Rouse number
         WHERE (maskatc_int)
            array2dc3 = ckar*beta_sed(:,:,f)*&
                      & SQRT(bstresatc_sed(1:ncloc,1:nrloc))
         ELSEWHERE
            array2dc3 = 0.0
         END WHERE
         WHERE (array2dc3.GT.0.0)
            rouse = MIN(wfall(1:ncloc,1:nrloc,1,f)/array2dc3,rousemax_d)
         ELSEWHERE
            rouse = rousemax_d
         END WHERE
!        ---constant term in integral
         WHERE (maskatc_int)
            array2dc3 = (1.0/array2dc1-1.0)**(-rouse)
         END WHERE
!        ---depth locations for integration
         n_2111: DO n=1,nrquad_sed
            WHERE (maskatc_int)
               elev_depth_avg(:,:,n) = array2dc2*(1.0-nodes_depth_avg(n))
            END WHERE
         ENDDO n_2111
!        ---log-transformed Rouse profile 
         n_2112: DO n=1,nrquad_sed
            WHERE (maskatc_int)
               concprof(:,:,n) = EXP(elev_depth_avg(:,:,n))*&
                              & (EXP(-elev_depth_avg(:,:,n))-1.0)**rouse
            END WHERE
         ENDDO n_2112
!        ---integration using gauss-legendre quadrature
         conc = 0.0
         n_2113: DO n=1,nrquad_sed
            WHERE (maskatc_int)
               conc = conc + weights(:,:,n)*concprof(:,:,n)
            END WHERE
         ENDDO n_2113
         WHERE (maskatc_int)
            ceq(:,:,f) = -array2dc2*array2dc3*cref(1:ncloc,1:nrloc,f)*conc
         END WHERE
      ENDDO f_211

!
!2.1.2 From suspended load
!-------------------------
!

   ELSEIF (iopt_sed_bbc_eq.GT.1) THEN

!  ---suspended load
      CALL sediment_suspendedload

!  ---current magnitude
      CALL vector_mag_arr_atc(udvel(1:ncloc+1,1:nrloc),&
                            & vdvel(1:ncloc,1:nrloc+1),1,1,1,1,iarr_hdvelmag,&
                            & .TRUE.,vecmag=array2dc1)
      WHERE (maskatc_int)
         array2dc1 = MAX(array2dc1,hdcurmin)
      END WHERE

!  ---equilibrium concentration
      f_212: DO f=1,nf
         WHERE (maskatc_int)
            ceq(:,:,f) = qsusatc(:,:,f)/array2dc1
         END WHERE
      ENDDO f_212

   ENDIF
   
!
!2.2 3-D case
!------------
!
   
ELSEIF (iopt_sed_nodim.EQ.3) THEN

   f_220: DO f = 1,nf
!     ---Rouse number
      WHERE (maskatc_int)
         array2dc1 = ckar*SQRT(bstresatc_sed(1:ncloc,1:nrloc))*beta_sed(:,:,f)
      ELSEWHERE
         array2dc1 = 0.0
      END WHERE
      WHERE (array2dc1.GT.0.0)
         rouse =  MIN(wfall(1:ncloc,1:nrloc,1,f)/array2dc1,rousemax_c)
      ELSEWHERE
         rouse = rousemax_c
      END WHERE

!     ---extrapolation from the rouse profile
      WHERE (maskatc_int)
         array2dc1  = ((1.0/gscoordatc(1:ncloc,1:nrloc,1)-1.0)*&
                    & (height_c(:,:,f)/(deptotatc(1:ncloc,1:nrloc)-&
                     & height_c(:,:,f))))**rouse
         ceq(1:ncloc,1:nrloc,f) =  array2dc1*cref(1:ncloc,1:nrloc,f)
         ceq(1:ncloc,1:nrloc,f) = MIN(ceq(1:ncloc,1:nrloc,f),cmax)
      END WHERE
   ENDDO f_220

ENDIF   
   
!
!3. Equilibrium time scale
!-------------------------
!

IF (iopt_sed_nodim.EQ.2) THEN

   f_310: DO f=1,nf
      WHERE (bstresatc_sed(1:ncloc,1:nrloc).GT.0.0)
         rouse = MIN(wfall(1:ncloc,1:nrloc,1,f)/&
                   & SQRT(bstresatc_sed(1:ncloc,1:nrloc)),rousemax_t)
      ELSEWHERE
         rouse = rousemax_t
      END WHERE
      WHERE (maskatc_int)
         t_equil(:,:,f) = p1*rouse**8 + p2*rouse**7+ p3*rouse**6 + &
                        & p4*rouse**5 + p5*rouse**4 + p6*rouse**3 + &
                        & p7*rouse**2 + p8*rouse + p9
      END WHERE
   ENDDO f_310

ENDIF

!
!4. Deallocate
!-------------
!

DEALLOCATE (array2dc1,rouse)
IF (iopt_sed_nodim.EQ.2.AND.iopt_sed_bbc_eq.EQ.1) THEN
   DEALLOCATE (array2dc2,array2dc3,conc,concprof,elev_depth_avg)
ENDIF

IF ((cold_start.OR.nt.EQ.nstep).AND.iopt_sed_nodim.EQ.2.AND.&
   & iopt_sed_bbc_eq.EQ.1) THEN
   DEALLOCATE(nodes_depth_avg,weights,weights_depth_avg)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE equilibrium_concentration

!========================================================================

SUBROUTINE median_particle_diameter(dsize,pctsize)
!***********************************************************************
!
! *median_particle_diameter* determines the xpct% particle size
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.10
!
! Description - 
!
! Reference -
!
! Calling program - initialise_sediment_arrays
!
! External calls - 
!
! Module calls - qsort_index
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE sedarrays
USE sedpars
USE sedswitches
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: qsort_index

IMPLICIT NONE

!
!*Arguments
!
REAL, INTENT(IN) :: pctsize
REAL, INTENT(OUT), DIMENSION(0:ncloc,0:nrloc) :: dsize
!
! Name      Type    Purpose
!-------------------------------------------------------------------------------
!*dsize*    REAL    % particle size                                          [m]
!*pctsize*  REAL    the % size
!
!*******************************************************************************
!
!* Local variables
!
INTEGER :: f, fc, i, j
INTEGER, DIMENSION(nf) :: sort_array
REAL :: xpct
REAL, DIMENSION(nf+1) :: cumfraction, dp_edges
REAL, DIMENSION(nf-1) :: dp_diff
REAL, DIMENSION(nf) :: dp_work, fraction_work


procname(pglev+1) = 'median_particle_diameter'
CALL log_timer_in()

!
!1. Case of one fraction only
!----------------------------
!

IF (nfp.EQ.1) THEN
   dsize = dp(1)
   GOTO 1000
ENDIF

fc = nf; xpct = 0.01*pctsize

!
!2. Without interpolation
!------------------------
!

IF (iopt_sed_median.EQ.1) THEN

   j_210: DO j=1,nrloc
   i_210: DO i=1,ncloc
      IF (maskatc_int(i,j)) THEN
         IF (SUM(bed_fraction(i,j,1,:)).EQ.0.0) THEN
            dsize(i,j) = 0.0
         ELSE
            CALL qsort_index(dp(1:nf),sort_array,1)
            dp_work = dp(sort_array)
            fraction_work = bed_fraction(i,j,1,sort_array)
            cumfraction(1) = 0.0
            f_111: DO f=2,nf+1
               cumfraction(f) = cumfraction(f-1) + fraction_work(f-1)
               IF (cumfraction(f).GE.xpct.AND.cumfraction(f-1).LE.xpct) fc = f-1
            ENDDO f_111
            dsize(i,j) = dp_work(fc)
         ENDIF
      ENDIF
   ENDDO i_210
   ENDDO j_210

!            
!3. With linear interpolation
!----------------------------
!

ELSEIF (iopt_sed_median.EQ.2) THEN

   j_310: DO j=1,nrloc
   i_310: DO i=1,ncloc
      IF (maskatc_int(i,j)) THEN
         IF (SUM(bed_fraction(i,j,1,:)).EQ.0.0) THEN
            dsize(i,j) = 0.0
         ELSE
!           ---sort data (mainly only necessary on first time)
            CALL qsort_index(dp(1:nf),sort_array,1)
            dp_work = dp(sort_array)
            fraction_work = bed_fraction(i,j,1,sort_array)
!           ---extrapolation of the edges of the particle diameter histogram
            dp_diff = dp_work(2:nf)-dp_work(1:nf-1)
            dp_edges(1) = MAX(dp_work(1)-0.5*dp_diff(1),0.0)
            dp_edges(2:nf) = dp_work(1:nf-1)+0.5*dp_diff(1:nf-1)
            dp_edges(nf+1)  = MAX(dp_work(nf)+0.5*dp_diff(nf-1),0.0)
!           ---cumulative distribution function
            cumfraction(1) = 0.0
            f_311: DO f=2,nf+1
               cumfraction(f) = cumfraction(f-1) + fraction_work(f-1)
               IF (cumfraction(f).GE.xpct.AND.cumfraction(f-1).LE.xpct) fc = f-1
            ENDDO f_311
!           ---linear interpolation
            dsize(i,j) = dp_edges(fc) + (dp_edges(fc+1)-dp_edges(fc))*&
                       & (xpct-cumfraction(fc)) &
                       & /(cumfraction(fc+1)-cumfraction(fc))
         ENDIF
      ENDIF
   ENDDO i_310
   ENDDO j_310

ENDIF

1000 CALL log_timer_out()


RETURN

END SUBROUTINE median_particle_diameter

!========================================================================

SUBROUTINE reference_concentration
!************************************************************************
!
! *reference_concentration* reference height and concentration
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Bottom_Fluxes.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - sediment_advdiff, sediment_suspendedload,
!                   sediment_totalload
!
! External calls -
!
! Module calls - error_alloc, exchange_mod
!
!************************************************************************
!
USE depths
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
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE
!
!* Local variables
!

INTEGER :: f
REAL, PARAMETER :: onethird = 1.0/3.0
REAL :: maxRV, minRV
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: densnorm, dstar, taurat


procname(pglev+1) = 'reference_concentration'
CALL log_timer_in()

!
!1. Allocate
!-----------
!

ALLOCATE (densnorm(ncloc,nrloc),STAT=errstat)
CALL error_alloc('densnorm',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (taurat(ncloc,nrloc),STAT=errstat)
CALL error_alloc('taurat',2,(/ncloc,nrloc/),kndrtype)
IF (iopt_sed_bbc_ref.EQ.2) THEN
   ALLOCATE (dstar(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('dstar',2,(/ncloc,nrloc/),kndrtype)
ENDIF

!
!2. Reference height and concentration
!-------------------------------------   
!
   
SELECT CASE (iopt_sed_bbc_ref)

!   
!2.1. Smith and McLean (1997)
!----------------------------
!

   CASE (1)

      maxRV = 0.1; minRV = 1.0E-05
      f_210: DO f=1,nf
!        ---reference height
         WHERE (maskatc_int)
            densnorm = MAX(rhos(f)/densw(1:ncloc,1:nrloc,1)-1.0,0.0)
            height_c(:,:,f) = 30.0*zroughatc_sed(1:ncloc,1:nrloc) + &
                            & 26.3*MAX((bstresatc_sed(1:ncloc,1:nrloc)-&
                                      & bstres_cr(1:ncloc,1:nrloc,f)),0.0)& 
                            & /(densnorm*gaccatc(1:ncloc,1:nrloc))
            height_c(:,:,f) = MIN(height_c(:,:,f),&
                                & maxRV*deptotatc(1:ncloc,1:nrloc))
            height_c(:,:,f) = MAX(height_c(:,:,f),&
                                & minRV*deptotatc(1:ncloc,1:nrloc))
         ELSEWHERE
            height_c(:,:,f) = 0.0
         END WHERE
         WHERE (maskatc_int)
            taurat = MAX(bstresatc_sed(1:ncloc,1:nrloc)&
                       & /bstres_cr(1:ncloc,1:nrloc,f)-1.0,0.0)
            cref(1:ncloc,1:nrloc,f) = 0.0024*cmax*(taurat/(1.0+0.0024*taurat))
            cref(1:ncloc,1:nrloc,f) = cref(1:ncloc,1:nrloc,f)*&
                                    & bed_fraction(1:ncloc,1:nrloc,1,f)
         ELSEWHERE
            cref(1:ncloc,1:nrloc,f) = 0.0
         END WHERE
      ENDDO f_210

!
!2.2 van Rijn (1984)
!-------------------
!

   CASE (2)
      maxRV = 0.1; minRV = 0.01
      WHERE (maskatc_int)
         height_c(:,:,1) = MIN(30.0*zroughatc_sed(1:ncloc,1:nrloc),&
                             & maxRV*deptotatc(1:ncloc,1:nrloc))
         height_c(:,:,1) = MAX(height_c(:,:,1),minRV*deptotatc(1:ncloc,1:nrloc))
      ELSEWHERE
         height_c(:,:,1) = 0.0
      END WHERE
 
      IF (nf.GT.1) THEN
         f_221: DO f=2,nf
            height_c(:,:,f) = height_c(:,:,1)
         ENDDO f_221
      ENDIF
      f_222: DO f = 1,nf
         WHERE (maskatc_int)
            densnorm  = MAX(rhos(f)/densw(1:ncloc,1:nrloc,1)-1.0,0.0)
            taurat = MAX((bstresatc_sed(1:ncloc,1:nrloc)&
                       & /bstres_cr(1:ncloc,1:nrloc,f)-1.0),0.0)
            dstar = dp(f)*(densnorm*gaccatc(1:ncloc,1:nrloc)&
                 & /kinvisc(1:ncloc,1:nrloc,1)**2)**onethird
            cref(1:ncloc,1:nrloc,f) = (0.015*dp(f))*taurat**1.5/dstar**0.3&
                 & /height_c(:,:,f)
            cref(1:ncloc,1:nrloc,f) = cref(1:ncloc,1:nrloc,f)*&
                                    & bed_fraction(1:ncloc,1:nrloc,1,f)
         ELSEWHERE
            cref(1:ncloc,1:nrloc,f) = 0.0   
         END WHERE
      ENDDO f_222

END SELECT

!
!3. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(cref,(/0,0,1/),(/1,0,1,0/),iarr_cref)
ENDIF

!
!4. Deallocate
!-------------
!

DEALLOCATE (densnorm,taurat)
IF (iopt_sed_bbc_ref.EQ.2) DEALLOCATE (dstar)

CALL log_timer_out()


RETURN

END SUBROUTINE reference_concentration

!========================================================================

SUBROUTINE sediment_advdiff
!************************************************************************
!
! *sediment_advdiff* Solve suspended sediment transport by advection-diffusion
!
! Author - IMDC, Boudewijn Decrop
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - sediment_equation
!
! External calls - bartnicki_filter, beta_factor,
!                  define_profobc_spec, diff_coef_waves,
!                  equilibrium_concentration, open_boundary_conds_prof,
!                  reference_concentration, sediment_settling_velocity,
!                  transport_at_C_4d1, transport_at_C_4d2,
!                  update_nest_data_prof, update_profobc_data
!
! Module calls - Carr_at_U, Carr_at_V, error_alloc, exchange_mod, mult_index 
!
!************************************************************************
!
USE currents
USE dararrays
USE depths
USE diffusion
USE grid
USE gridpars
USE iopars
USE morphpars
USE morpharrays
USE nestgrids
USE physpars
USE relaxation
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE structures
USE switches
USE syspars
USE timepars
USE array_interp, ONLY: Carr_at_U, Carr_at_V
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: update
LOGICAL, SAVE :: disch_data, settling
INTEGER :: f, i, ibcbot, ii, iiloc, iu, j, jj, jjloc, jv, k, &
         & morphsteps, nbcb
INTEGER, SAVE :: maxprofs, noprofs, nofiles
INTEGER, DIMENSION(4) :: nhexch
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: iprofrlx, itypobu, itypobv, &
                                          & ivarids, noprofsd
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: indexprof, indexvar, iprofobu, &
                                            & iprofobv
INTEGER (KIND=kndilong), SAVE, DIMENSION(2) :: nosecsdis
INTEGER (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:,:) :: nosecsprof
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: dissed, taurat, weight1, weight2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3d, disdata, obcsed3d, &
                                           & sedobu, sedobv, zeros3d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: bcbot, profsed, sedsource


procname(pglev+1) = 'sediment_advdiff'
CALL log_timer_in()

!
!1. Initialise open boundary conditions on first call
!----------------------------------------------------
!

IF (nt.EQ.0) THEN

   nofiles = maxdatafiles(io_sedobc,1)
   disch_data = iopt_dischr.EQ.1.AND.modfiles(io_dissed,1,1)%defined

!  ---variables ids
   IF (iopt_transp_full.EQ.0) THEN
      ALLOCATE (ivarids(nf),STAT=errstat)
      CALL error_alloc('ivarids',1,(/nf/),kndint)
      ivarids = iarr_cvol
   ENDIF

!  ---allocate specifier arrays and o.b. profiles
   ALLOCATE (sedobu(nobu,nz,nf),STAT=errstat)
   CALL error_alloc('sedobu',3,(/nobu,nz,nf/),kndrtype)
   IF (nobu.GT.0) sedobu = real_fill
   ALLOCATE (sedobv(nobv,nz,nf),STAT=errstat)
   CALL error_alloc('sedobv',3,(/nobv,nz,nf/),kndrtype)
   IF (nobv.GT.0) sedobv = real_fill
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
        ALLOCATE (obcsed3d(0:noprofs,nz,nf),STAT=errstat)
        CALL error_alloc('obcsed3d',3,(/noprofs+1,nz,nf/),kndrtype)
        obcsed3d = real_fill
   ENDIF
   IF (nofiles.GT.1) THEN
      maxprofs = MAXVAL(noprofsd)
      ALLOCATE (nosecsprof(2:nofiles,2),STAT=errstat)
      CALL error_alloc('nosecsprof',2,(/nofiles-1,2/),kndilong)
      ALLOCATE (profsed(maxprofs,nz,2:nofiles,2),STAT=errstat)
      CALL error_alloc('profsed',4,(/maxprofs,nz,nofiles-1,2/),kndrtype)
   ENDIF

!  ---allocate discharge data arrays
   IF (disch_data) THEN
      ALLOCATE (dissed(numdis,nf),STAT=errstat)
      CALL error_alloc('dissed',2,(/numdis,nf/),kndrtype)
      ALLOCATE (disdata(numdis,nf,2),STAT=errstat)
      CALL error_alloc('disdata',3,(/numdis,nf,2/),kndrtype)
   ENDIF

!  ---update data at initial time
   IF (nofiles.GT.1) THEN
      CALL update_profobc_data(profsed,obcsed3d,noprofsd,&
                             & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                             & maxprofs,noprofs,nf,nofiles,nosecsprof,io_sedobc)
   ENDIF

!  ---update discharge data at initial time
   IF (disch_data) THEN
      CALL update_dischr_data(io_dissed,1,disdata,dissed,numdis,nf,nosecsdis,&
                            & update)
   ENDIF

!  ---write at nest locations
   IF (iopt_nests.EQ.1) THEN 
      CALL update_nest_data_prof(cvol,nf,nosednst,instsed,io_sednst)
   ENDIF

!
!1.3 Settling velocity
!---------------------
!

   settling = iopt_sed_vadv.GT.0
   IF (settling.OR.iopt_sed_nodim.EQ.2) CALL sediment_settling_velocity

!
!1.4 Vertical diffusion coefficient
!----------------------------------
!

   CALL beta_factor

   f_140: DO f=1,nf
   k_140: DO k=1,nz+1
      WHERE (maskatc_int)
         vdiffcoef_sed(:,:,k,f) = beta_sed(:,:,f)*vdifcoefmom(1:ncloc,1:nrloc,k)
      END WHERE
      ENDDO k_140
   ENDDO f_140

!
!1.5 Reference height and concentrations
!---------------------------------------
!
 
   IF (iopt_sed_bbc.EQ.2) THEN
      CALL reference_concentration
      CALL equilibrium_concentration
   ENDIF

   GOTO 1000

ENDIF

!
!2. Allocate and initialise
!--------------------------
!

ALLOCATE (sedsource(ncloc,nrloc,nz,nf),STAT=errstat)
CALL error_alloc('sedsource',4,(/ncloc,nrloc,nz,nf/),kndrtype)
ALLOCATE (zeros3d(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('zeros3d',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (bcbot(ncloc,nrloc,2,nf),STAT=errstat)
CALL error_alloc('bcbot',4,(/ncloc,nrloc,2,nf/),kndrtype)
ALLOCATE (array3d(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('array3d',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
IF (iopt_sed_bbc.EQ.1) THEN
   ALLOCATE (taurat(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('taurat',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (iopt_scal_depos.EQ.2) THEN
   ALLOCATE (weight1(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('weight1',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (weight2(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('weight2',2,(/ncloc,nrloc/),kndrtype)
ENDIF

sedsource = 0.0; zeros3d = 0.0

!
!3. Open boundary conditions
!---------------------------
!
!---update data
IF (nofiles.GT.1) THEN
   CALL update_profobc_data(profsed,obcsed3d,noprofsd,&
                          & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                          & maxprofs,noprofs,nf,nofiles,nosecsprof,io_sedobc)
ENDIF

!---update o.b. profiles
IF (iopt_obc_sed.EQ.1) THEN
   CALL open_boundary_conds_prof(cvol,sedobu,sedobv,obcsed3d,obcsedatu,&
                               & obcsedatv,itypobu,itypobv,iprofobu,iprofobv,&
                               & noprofs,nf)
ENDIF

!---update discharge data
IF (disch_data) THEN
   CALL update_dischr_data(io_dissed,1,disdata,dissed,numdis,nf,nosecsdis,&
                         & update)
ENDIF

!
!4. Settling velocity
!--------------------
!

IF (settling.OR.iopt_sed_nodim.EQ.2) CALL sediment_settling_velocity

!
!5. Vertical diffusion coefficient
!---------------------------------
!
!5.1 Beta factor
!---------------
!

IF (iopt_sed_beta.EQ.3) CALL beta_factor

!
!5.2 Vertical diffusion coefficient
!----------------------------------
!

f_520: DO f=1,nf
k_520: DO k=1,nz+1
   WHERE (maskatc_int)
      vdiffcoef_sed(:,:,k,f) = beta_sed(:,:,f)*vdifcoefmom(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_520
ENDDO f_520

!---include wave effects
IF (iopt_sed_wave_diff.GT.0) CALL diff_coef_waves

!
!6. Source terms
!---------------
!

f_610: DO f=1,nf
k_610: DO k=1,nz
   WHERE (maskatc_int)
      sedsource(:,:,k,f) = dar_sediment(:,:,k,f)/&
                         & (garea*deptotatc(1:ncloc,1:nrloc))
   END WHERE
ENDDO k_610
ENDDO f_610

!---add discharges
IF (disch_data) CALL scalar_discharge(dissed,sedsource,nf,2) 

!
!7. Bottom fluxes
!-----------------
!
!7.1 Reference height and concentrations
!---------------------------------------
!
 
IF (iopt_sed_bbc.EQ.2) THEN
   CALL reference_concentration
   CALL equilibrium_concentration
ENDIF

!
!7.2 Initialise
!--------------
!

bottom_sed_ero = 0.0; bottom_sed_dep = 0.0

!
!7.3 2-D case
!------------
!

IF (iopt_sed_nodim.EQ.2) THEN

!
!7.3.1 Standard case
!-------------------
!
   

   IF (iopt_sed_bbc.EQ.1) THEN

      ibcbot = 3; nbcb = 2
      f_731: DO f=1,nf
         WHERE (maskatc_int)
            taurat = MAX(bstresatc_sed(1:ncloc,1:nrloc)&
                      & /bstres_cr(1:ncloc,1:nrloc,f)-1.0,0.0)
            bottom_sed_ero(:,:,f) = bed_fraction(1:ncloc,1:nrloc,1,f)*&
                                  & (parth_coef/rhos(f))*(taurat**parth_exp)
            bottom_sed_dep(:,:,f) = (1.0-theta_vdif)*&
                                  & wfall(1:ncloc,1:nrloc,1,f)*&
                                  & cvol(1:ncloc,1:nrloc,1,f)
            bcbot(:,:,1,f) = bottom_sed_dep(:,:,f) - bottom_sed_ero(:,:,f)
            bcbot(:,:,2,f) = theta_vdif*wfall(1:ncloc,1:nrloc,1,f)
         END WHERE
      ENDDO f_731
      
!
!7.3.2 Equilibrium method
!------------------------
!      

   ELSE

      ibcbot = 2; nbcb = 2
      f_732: DO f=1,nf
         WHERE (maskatc_int)
            bcbot(:,:,1,f) = ceq(:,:,f)
            bcbot(:,:,2,f) = wfall(1:ncloc,1:nrloc,1,f)/t_equil(:,:,f)
            bottom_sed_ero(:,:,f) = bcbot(:,:,2,f)*ceq(:,:,f)
            bottom_sed_dep(:,:,f) = (1.0-theta_vdif)*bcbot(:,:,2,f)*&
                                  & cvol(1:ncloc,1:nrloc,1,f)
         END WHERE
      ENDDO f_732

   ENDIF


!
!7.4 3-D case
!------------
!

ELSEIF (iopt_sed_nodim.EQ.3) THEN

!
!7.4.1 Standard case
!-------------------
!
   
   IF (iopt_sed_bbc.EQ.1) THEN
      
      ibcbot = 1; nbcb = 1
      f_741: DO f=1,nf
         WHERE (maskatc_int)
            taurat = MAX(bstresatc_sed(1:ncloc,1:nrloc)&
                      & /bstres_cr(1:ncloc,1:nrloc,f)-1.0,0.0)
            bottom_sed_ero(:,:,f)  = bed_fraction(1:ncloc,1:nrloc,1,f)*&
                                   & (parth_coef/rhos(f))*(taurat**parth_exp)
            bcbot(:,:,1,f) = -bottom_sed_ero(:,:,f)
         END WHERE
      ENDDO f_741

!
!7.4.2 Equilibrium method
!------------------------
!      

   ELSE
      
      ibcbot = 1; nbcb = 1            
      f_742: DO f=1,nf
         WHERE (maskatc_int)
            bottom_sed_ero(:,:,f) = wfall(1:ncloc,1:nrloc,1,f)*ceq(:,:,f)
            bcbot(:,:,1,f) = -bottom_sed_ero(:,:,f)
         END WHERE
      ENDDO f_742
   ENDIF

!
!7.4.3 Explicit part of deposition flux
!--------------------------------------
!

   IF (settling.AND.iopt_vadv_impl.LT.2) THEN

!     ---first order scheme
      IF (iopt_scal_depos.EQ.1) THEN
         f_7431: DO f=1,nf
            bottom_sed_dep(:,:,f) = (1.0-theta_vadv)*&
                                  & wfall(1:ncloc,1:nrloc,1,f)*&
                                  & cvol(1:ncloc,1:nrloc,1,f)
         ENDDO f_7431
!     ---second order scheme      
      ELSEIF (iopt_scal_depos.EQ.2) THEN
         WHERE (maskatc_int)
            weight1 = 0.5*delzatc(1:ncloc,1:nrloc,1)/delzatw(1:ncloc,1:nrloc,2)
            weight2 = 1.0 + weight1
         END WHERE
         f_7432: DO f=1,nf
            WHERE (maskatc_int)
               bottom_sed_dep(:,:,f) = (1.0-theta_vadv)*&
                                      & wfall(1:ncloc,1:nrloc,1,f)*&
                                      & (weight2*cvol(1:ncloc,1:nrloc,1,f)-&
                                      &  weight1*cvol(1:ncloc,1:nrloc,2,f))
            END WHERE
         ENDDO f_7432
      ENDIF

   ENDIF

ENDIF

!
!8. Solve transport equation
!---------------------------
!
   
IF (iopt_transp_full.EQ.0) THEN
   CALL transport_at_C_4d1(cvol,sedsource,wfall,vdiffcoef_sed,nf,&
           & iopt_adv_scal,iopt_sed_vadv,iopt_hdif_scal,sedobu,sedobv,&
           & iprofrlx,0,ibcbot,1,nbcb,zeros3d,&
           & bcbot(:,:,1:nbcb,:),settling,ivarids)
ELSEIF (iopt_transp_full.EQ.1) THEN
   CALL transport_at_C_4d2(cvol,sedsource,wfall,vdiffcoef_sed,nf,&
           & iopt_adv_scal,iopt_sed_vadv,iopt_hdif_scal,sedobu,sedobv,&
           & iprofrlx,0,ibcbot,1,nbcb,zeros3d,&
           & bcbot(:,:,1:nbcb,:),settling,iarr_cvol)
ENDIF

!
!9. Add implicit parts to bottom flux
!------------------------------------
!
!9.1 2-D case
!------------
!

IF (iopt_sed_nodim.EQ.2.AND.iopt_vdif_impl.GT.1) THEN

!
!9.1.1 Standard case
!-------------------
!

   IF (iopt_sed_bbc.EQ.1) THEN

      f_911: DO f=1,nf
         WHERE (maskatc_int)
            bottom_sed_dep(:,:,f) = bottom_sed_dep(:,:,f) + &
                                  & theta_vdif*wfall(1:ncloc,1:nrloc,1,f)*&
                                  & cvol(1:ncloc,1:nrloc,1,f)
         END WHERE
      ENDDO f_911

!
!9.1.2 Equilibrium method 
!------------------------
!

   ELSEIF (iopt_sed_bbc.EQ.2) THEN
      f_912: DO f=1,nf
         bottom_sed_dep(:,:,f) = bottom_sed_dep(:,:,f) + &
                               & theta_vdif*bcbot(:,:,2,f)*&
                               & cvol(1:ncloc,1:nrloc,1,f)
      ENDDO f_912
   ENDIF
      
!
!9.2 3-D case
!------------
!

ELSEIF (iopt_sed_nodim.EQ.3.AND.iopt_vadv_impl.GT.0.AND.settling) THEN
      
!  ---first order scheme
   IF (iopt_scal_depos.EQ.1) THEN
      f_9211: DO f=1,nf
         WHERE (maskatc_int)
            bottom_sed_dep(:,:,f) = bottom_sed_dep(:,:,f) + &
                                  & theta_vadv*wfall(1:ncloc,1:nrloc,1,f)*&
                                  & cvol(1:ncloc,1:nrloc,1,f)
         END WHERE
      ENDDO f_9211

!  ---second order scheme
   ELSEIF (iopt_scal_depos.EQ.2) THEN
      WHERE (maskatc_int)
         weight1 = 0.5*delzatc(1:ncloc,1:nrloc,1)/delzatw(1:ncloc,1:nrloc,2)
         weight2 = 1.0 + weight1
      END WHERE
      f_9222: DO f=1,nf
         WHERE (maskatc_int)
            bottom_sed_dep(:,:,f) = bottom_sed_dep(:,:,f) - &
                                  & theta_vadv*wfall(1:ncloc,1:nrloc,1,f)*&
                                  & (weight2*cvol(1:ncloc,1:nrloc,1,f)-&
                                  &  weight1*cvol(1:ncloc,1:nrloc,2,f))
         END WHERE
      ENDDO f_9222
   ENDIF

ENDIF

!
!10. Bottom sediment flux
!------------------------
!

f_1010: DO f=1,nf
   WHERE (maskatc_int)
      bottom_sed_flux(:,:,f) = bottom_sed_ero(:,:,f) - bottom_sed_dep(:,:,f)
   END WHERE
ENDDO f_1010

!
!11. Apply Bartnciki filter to eliminate negative concentrations
!---------------------------------------------------------------
!

IF (iopt_sed_filter.EQ.1) THEN
   CALL bartnicki_filter(cvol,nf,iarr_cvol)
   
!
!12. Exchange halo sections
!--------------------------
!

ELSEIF (iopt_MPI.EQ.1) THEN
   nhexch = nhalo
   CALL exchange_mod(cvol,(/1-nhalo,1-nhalo,1,1/),nhexch,iarr_cvol)
ENDIF

!
!13. Write at nest locations
!---------------------------
!

IF (iopt_nests.EQ.1) THEN
   CALL update_nest_data_prof(cvol,nf,nosednst,instsed,io_sednst)
ENDIF

!
!14. Total sediment concentration
!--------------------------------
!

k_1400: DO k=1,nz
   WHERE (nodeatc(0:ncloc+1,0:nrloc+1).EQ.1)
     ctot(:,:,k) = SUM(cvol(0:ncloc+1,0:nrloc+1,k,:),DIM=3)
   END WHERE
ENDDO k_1400

!
!15. Suspended and total load
!----------------------------
!
!15.1 U-nodes
!------------
!
!---suspended load
qsusatu = 0.0
f_1511: DO f=1,nf
   CALL Carr_at_U(cvol(-1:ncloc+1,1:nrloc,:,f),array3d(:,1:nrloc,:),1,3,&
               & (/-1,1,1/),(/ncloc+1,nrloc,nz/),1,iarr_cvol,.TRUE.)
   k_15111: DO k=1,nz
      WHERE (node2du(0:ncloc+1,1:nrloc).GT.1)
         qsusatu(:,:,k,f) = ufvel(0:ncloc+1,1:nrloc,k)*array3d(:,1:nrloc,k)
      END WHERE
   ENDDO k_15111
ENDDO f_1511

!---total load
qtotatu = 0.0
f_1512: DO f=1,nf
   WHERE (node2du(0:ncloc+1,1:nrloc).GT.1)
      qtotatu(:,:,f) = SUM(qsusatu(:,:,:,f)*delzatu,DIM=3)
   END WHERE
   IF (iopt_sed_mode.EQ.3) THEN
      WHERE (node2du(0:ncloc+1,1:nrloc).GT.1)
         qtotatu(:,:,f) = qtotatu(:,:,f) + qbedatu(:,:,f)
      END WHERE
   ENDIF
ENDDO f_1512

!
!15.2 V-nodes
!------------
!
!---suspended load
qsusatv = 0.0
f_1521: DO f=1,nf
   CALL Carr_at_V(cvol(1:ncloc,-1:nrloc+1,:,f),array3d(1:ncloc,:,:),1,3,&
               & (/1,-1,1/),(/ncloc,nrloc+1,nz/),1,iarr_cvol,.TRUE.)
   k_15211: DO k=1,nz
      WHERE (node2dv(1:ncloc,0:nrloc+1).GT.1)
         qsusatv(:,:,k,f) = vfvel(1:ncloc,0:nrloc+1,k)*array3d(1:ncloc,:,k)
      END WHERE
   ENDDO k_15211
ENDDO f_1521

!---total load
qtotatv = 0.0
f_1522: DO f=1,nf
   WHERE (node2dv(1:ncloc,0:nrloc+1).GT.1)
      qtotatv(:,:,f) = SUM(qsusatv(:,:,:,f)*delzatv,DIM=3)
   END WHERE
   IF (iopt_sed_mode.EQ.3) THEN
      WHERE (node2dv(1:ncloc,0:nrloc+1).GT.1)
         qtotatv(:,:,f) = qtotatv(:,:,f) + qbedatv(:,:,f)
      END WHERE
   ENDIF
ENDDO f_1522

!
!16. Upwind correction at open boundaries
!----------------------------------------
!

IF (iopt_sed_obc_flux.EQ.1) THEN

!
!16.1 U-nodes
!-----------
!

   iiloc_1610: DO iiloc=1,nobuloc
      i = iobuloc(iiloc); j = jobuloc(iiloc)
      ii = indexobu(iiloc); iu = MERGE(i+1,i-1,westobu(ii))
      qsusatu(i,j,:,:) = (delyatu(iu,j)/delyatu(i,j))*qsusatu(iu,j,:,:)
      qtotatu(i,j,:) = (delyatu(iu,j)/delyatu(i,j))*qtotatu(iu,j,:)
   ENDDO iiloc_1610
 
!
!16.2 V-nodes
!------------
!

   jjloc_1620: DO jjloc=1,nobvloc
      i = iobvloc(jjloc); j = jobvloc(jjloc)
      jj = indexobv(jjloc); jv = MERGE(j+1,j-1,soutobv(jj))
      qsusatv(i,j,:,:) = (delxatv(i,jv)/delxatv(i,j))*qsusatv(iu,j,:,:)
      qtotatv(i,j,:) = (delxatv(i,jv)/delxatv(i,j))*qtotatv(i,jv,:)
   ENDDO jjloc_1620

ENDIF

!
!17. Suspended load averaged over morphological step
!---------------------------------------------------
!

IF (iopt_morph.EQ.1) THEN
   morphsteps = icmorph/icsed
   f_1710: DO f=1,nf
      IF (icmorph.EQ.icsed) THEN
         WHERE (node2du(0:ncloc+1,1:nrloc).GT.1)
            qsusatu_int(:,:,f) = SUM(qsusatu(:,:,:,f)*delzatu,DIM=3)
         END WHERE
         WHERE (node2dv(1:ncloc,0:nrloc+1).GT.1)
            qsusatv_int(:,:,f) =  SUM(qsusatv(:,:,:,f)*delzatv,DIM=3)
         END WHERE
      ELSE
         IF (mult_index(nt,icmorph)) THEN
            qsusatu_int(:,:,f) = 0.0; qsusatv_int(:,:,f) = 0.0
         ENDIF
         WHERE (node2du(0:ncloc+1,1:nrloc).GT.1)
            qsusatu_int(:,:,f) = qsusatu_int(:,:,f) + &
                 & SUM(qsusatu(:,:,:,f)*delzatu,DIM=3)/morphsteps
         END WHERE
         WHERE (node2dv(1:ncloc,0:nrloc+1).GT.1)
            qsusatv_int(:,:,f) = qsusatv_int(:,:,f) + &
                 & SUM(qsusatv(:,:,:,f)*delzatv,DIM=3)/morphsteps
         END WHERE
      ENDIF
   ENDDO f_1710
ENDIF

!
!18. Deallocate arrays
!---------------------
!
!---work space
DEALLOCATE (array3d,bcbot,sedsource,zeros3d)
IF (iopt_sed_bbc.EQ.1) DEALLOCATE (taurat)
IF (iopt_scal_depos.EQ.2) DEALLOCATE (weight1,weight2)

!---at final step
1000 CONTINUE
IF (cold_start.OR.nt.EQ.nstep) THEN
   IF (iopt_transp_full.EQ.0) DEALLOCATE (ivarids)
   DEALLOCATE (sedobu,sedobv,iprofrlx)
   IF (iopt_obc_sed.EQ.1) THEN
      DEALLOCATE (itypobu,itypobv,iprofobu,iprofobv,noprofsd,indexprof,&
                & indexvar,obcsed3d)
   ENDIF
   IF (nofiles.GT.1) DEALLOCATE (nosecsprof,profsed)
   IF (disch_data) DEALLOCATE (dissed,disdata)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE sediment_advdiff

!========================================================================

SUBROUTINE sediment_bedload
!************************************************************************
!
! *sediment_bedload* Solve sediment bed load transport fluxes
!
! Author - Boudewijn Decrop (IMDC)
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - sediment_equation, sediment_totalload
!
! External calls -
!
! Module calls - Carr_at_U, Carr_at_V, error_alloc, exchange_mod, gauss_quad,
!                mult_index, vector_mag_arr_atu, vector_mag_arr_atv
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
USE morphpars
USE physpars
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Carr_at_U, Carr_at_V
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: gauss_quad, vector_mag_arr_atu, vector_mag_arr_atv
USE paral_comms, ONLY: exchange_mod 
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
INTEGER :: f, i, ii, iiloc, iu, j, jj, jjloc, jv, morphsteps, n
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
REAL :: onethird = 1.0/3.0
REAL :: xpi
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:) :: ndgau_tr, nodes_gauss,&
                                                       & weights_gauss
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: alpha, amplitude, array2d1, &
       & array2d2, array2d3, awave, beta, bphix, bphiy, bstresatu_wav, &
       & bstresatv_wav, curmag, delta, d50atu, d50atv, fcw, fw, ka, ks, &
       & oscwavetr, oscwavetr_per, peng, phi, phicuratu, phicuratv, phiwav, &
       & qprojx, qprojy, uc, xslope, yslope, wavevelatu, wavevelatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: deltasatu, deltasatv, dstaratu,&
       & dstaratv, excesshearatu, excesshearatv, qbedatu_per, qbedatv_per, &
       & qnormatu, qnormatv, taucritatu, taucritatv, taunormatu, taunormatv, &
       & velcurwav


procname(pglev+1) = 'sediment_bedload'
CALL log_timer_in()

!
!1.Allocate
!-----------
!

IF (nt.EQ.0.AND.iopt_waves.EQ.1.AND.iopt_sed_bedeq.GT.5) THEN
   ALLOCATE(weights_gauss(nrquad_wav),STAT=errstat)
   CALL error_alloc('weights_gauss',1,(/nrquad_wav/),kndrtype)
   weights_gauss =0.0
   ALLOCATE(nodes_gauss(nrquad_wav),STAT=errstat)
   CALL error_alloc('nodes_gauss',1,(/nrquad_wav/),kndrtype)
   nodes_gauss = 0.0
   ALLOCATE(ndgau_tr(nrquad_wav),STAT=errstat)
   CALL error_alloc('ndgau_tr',1,(/nrquad_wav/),kndrtype)      
   ndgau_tr =0.0
   CALL gauss_quad(nrquad_wav,weights_gauss,nodes_gauss)
   ndgau_tr = pi*(nodes_gauss+1.0)
ENDIF

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (array2d1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2d2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2d3(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d3',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (deltasatu(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('deltasatu',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (deltasatv(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('deltasatv',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (excesshearatu(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('excesshearatu',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (excesshearatv(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('excesshearatv',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (phicuratu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('phicuratu',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (phicuratv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('phicuratv',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (qprojx(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qprojx',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (qprojy(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qprojy',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (taucritatu(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('taucritatu',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (taucritatv(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('taucritatv',3,(/ncloc,nrloc,nf/),kndrtype)

IF (iopt_sed_bedeq.EQ.1.OR.iopt_sed_bedeq.EQ.2.OR.iopt_sed_bedeq.EQ.5) THEN
   ALLOCATE (taunormatu(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('taunormatu',3,(/ncloc,nrloc,nf/),kndrtype)
   ALLOCATE (taunormatv(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('taunormatv',3,(/ncloc,nrloc,nf/),kndrtype)
ENDIF
   
IF (iopt_sed_bedeq.EQ.3.OR.iopt_sed_bedeq.EQ.4) THEN
   ALLOCATE (qnormatu(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('qnormatu',3,(/ncloc,nrloc,nf/),kndrtype)
   ALLOCATE (qnormatv(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('qnormatv',3,(/ncloc,nrloc,nf/),kndrtype)
ENDIF

IF (iopt_sed_bedeq.EQ.2) THEN
   ALLOCATE (peng(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('peng',2,(/ncloc,nrloc/),kndrtype)
ENDIF

IF (iopt_sed_bedeq.EQ.3.OR.iopt_sed_bedeq.EQ.6.OR.iopt_sed_bedeq.EQ.7) THEN
   ALLOCATE (dstaratu(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('dstaratu',3,(/ncloc,nrloc,nf/),kndrtype)
   ALLOCATE (dstaratv(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('dstaratv',3,(/ncloc,nrloc,nf/),kndrtype)
ENDIF
   
IF (iopt_sed_bedeq.EQ.5) THEN
   ALLOCATE (bphix(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bphix',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (bphiy(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bphiy',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (bstresatu_wav(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bstresatu_wav',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (bstresatv_wav(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bstresatv_wav',2,(/ncloc,nrloc/),kndrtype)
ENDIF

IF (iopt_sed_bedeq.GE.5) THEN
   ALLOCATE (d50atu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('d50atu',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (d50atv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('d50atv',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (phi(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('phi',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (phiwav(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('phiwav',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (qbedatu_per(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('qbedatu_per',3,(/ncloc,nrloc,nf/),kndrtype)
   ALLOCATE (qbedatv_per(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('qbedatv_per',3,(/ncloc,nrloc,nf/),kndrtype)
ENDIF
   
IF (iopt_sed_bedeq.GE.6) THEN
   ALLOCATE (alpha(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('alpha',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (amplitude(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('amplitude',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (awave(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('awave',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (beta(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('beta',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (curmag(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('curmag',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (delta(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('delta',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (fcw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fcw',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (fw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fw',2,(/ncloc,nrloc/), kndrtype)
   ALLOCATE (ka(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('ka',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (ks(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('ks',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (oscwavetr(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('oscwavetr',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (oscwavetr_per(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('oscwavetr_per',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (uc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('uc',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (velcurwav(ncloc,nrloc,nrquad_wav),STAT=errstat)
   CALL error_alloc('velcurwav',3,(/ncloc,nrloc,nrquad_wav/),kndrtype)
   ALLOCATE (wavevelatu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('wavevelatu',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (wavevelatv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('wavevelatv',2,(/ncloc,nrloc/),kndrtype)
ENDIF

IF (iopt_sed_slope.EQ.1) THEN
   ALLOCATE (xslope(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('xslope',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (yslope(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('yslope',2,(/ncloc,nrloc/),kndrtype)
ENDIF

!
!2. Initialise arrays
!--------------------
!
!---masks
maskatu = nodeatu(1:ncloc,1:nrloc,1).GT.1
maskatv = nodeatv(1:ncloc,1:nrloc,1).GT.1

!---bed loads
qbedatu = 0.0; qbedatv = 0.0

!---critical and excess shear stress
f_210: DO f=1,nf
   CALL Carr_at_U(bstres_cr(0:ncloc,1:nrloc,f),taucritatu(:,:,f),1,3,(/0,1,1/),&
               & (/ncloc,nrloc,1/),1,iarr_bstres_cr,.TRUE.,&
               & mask=bed_fraction(:,1:nrloc,1,f).GT.0.0)
   CALL Carr_at_V(bstres_cr(1:ncloc,0:nrloc,f),taucritatv(:,:,f),1,3,(/1,0,1/),&
               & (/ncloc,nrloc,1/),1,iarr_bstres_cr,.TRUE.,&
               & mask=bed_fraction(1:ncloc,:,1,f).GT.0.0)
   WHERE (maskatu)
      excesshearatu(:,:,f) = MAX(bstresatu_sed-taucritatu(:,:,f),0.0)
   ELSEWHERE
      excesshearatu(:,:,f) = 0.0
   END WHERE
   WHERE (maskatv)
      excesshearatv(:,:,f) = MAX(bstresatv_sed-taucritatv(:,:,f),0.0)
   ELSEWHERE
      excesshearatv(:,:,f) = 0.0
   END WHERE
ENDDO f_210

!---normalising factors at U-nodes
CALL Carr_at_U(densw(0:ncloc,1:nrloc,1),array2d1,1,3,(/0,1,1/),&
            & (/ncloc,nrloc,1/),1,iarr_densw,.TRUE.)
CALL Carr_at_U(kinvisc(0:ncloc,1:nrloc,1),array2d2,1,3,(/0,1,1/),&
            & (/ncloc,nrloc,1/),1,iarr_kinvisc,.FALSE.)
f_220: DO f=1,nf
   WHERE (maskatu)
      deltasatu(:,:,f) = MAX(rhos(f)/array2d1-1.0,0.0)
   END WHERE
   IF (iopt_sed_bedeq.EQ.1.OR.iopt_sed_bedeq.EQ.2.OR.iopt_sed_bedeq.EQ.5) THEN
      WHERE (maskatu)
         taunormatu(:,:,f) = deltasatu(:,:,f)*gaccatu(1:ncloc,:)*dp(f)
      END WHERE
   ENDIF
   IF (iopt_sed_bedeq.EQ.3.OR.iopt_sed_bedeq.EQ.4) THEN
      WHERE (maskatu)
         qnormatu(:,:,f) = SQRT(deltasatu(:,:,f)*gaccatu(1:ncloc,:)*dp(f)**3)
      END WHERE
   ENDIF
   IF (iopt_sed_bedeq.EQ.3.OR.iopt_sed_bedeq.EQ.6.OR.iopt_sed_bedeq.EQ.7) THEN
      WHERE (maskatu)
         dstaratu(:,:,f) = dp(f)*(deltasatu(:,:,f)*gaccatu(1:ncloc,:)&
                               & /array2d2**2)**onethird
      END WHERE
   ENDIF
ENDDO f_220

!---normalising factors at V-nodes
CALL Carr_at_V(densw(1:ncloc,0:nrloc,1),array2d1,1,3,(/1,0,1/),&
            & (/ncloc,nrloc,1/),1,iarr_densw,.TRUE.)
CALL Carr_at_V(kinvisc(1:ncloc,0:nrloc,1),array2d2,1,3,(/1,0,1/),&
            & (/ncloc,nrloc,1/),1,iarr_kinvisc,.FALSE.)
f_230: DO f=1,nf
   WHERE (maskatv)
      deltasatv(:,:,f) = MAX(rhos(f)/array2d1-1.0,0.0)
   END WHERE
   IF (iopt_sed_bedeq.EQ.1.OR.iopt_sed_bedeq.EQ.2.OR.iopt_sed_bedeq.EQ.5) THEN
      WHERE (maskatv)
         taunormatv(:,:,f) = deltasatv(:,:,f)*gaccatv(:,1:nrloc)*dp(f)
      END WHERE
   ENDIF
   IF (iopt_sed_bedeq.EQ.3.OR.iopt_sed_bedeq.EQ.4) THEN
      WHERE (maskatv)
         qnormatv(:,:,f) = SQRT(deltasatv(:,:,f)*gaccatv(:,1:nrloc)*dp(f)**3)
      END WHERE
   ENDIF
   IF (iopt_sed_bedeq.EQ.3.OR.iopt_sed_bedeq.EQ.6.OR.iopt_sed_bedeq.EQ.7) THEN
      WHERE (maskatv)
         dstaratv(:,:,f) = dp(f)*(deltasatv(:,:,f)*gaccatv(:,1:nrloc)&
                               & /array2d2**2)**onethird
      END WHERE
   ENDIF
ENDDO f_230

!---median sizes at velocity nodes
IF (iopt_sed_bedeq.GE.5) THEN
   CALL Carr_at_U(d50_bed(:,1:nrloc),d50atu(1:ncloc,1:nrloc),1,3,(/0,1,1/),&
       & (/ncloc,nrloc,1/),1,iarr_d50_bed,.TRUE.,mask=d50_bed(:,1:nrloc).GT.0.0)
   CALL Carr_at_V(d50_bed(1:ncloc,:),d50atv(1:ncloc,1:nrloc),1,3,(/1,0,1/),&
       & (/ncloc,nrloc,1/),1,iarr_d50_bed,.TRUE.,mask=d50_bed(1:ncloc,:).GT.0.0)
ENDIF

!---current direction at velocity nodes
CALL vector_mag_arr_atu(uvel(1:ncloc,1:nrloc,1),vvel(0:ncloc,1:nrloc+1,1),&
                      & 1,3,1,1,iarr_hvelmag,.TRUE.,vecpha=phicuratu)
CALL vector_mag_arr_atv(uvel(1:ncloc+1,0:nrloc,1),vvel(1:ncloc,1:nrloc,1),&
                      & 1,3,1,1,iarr_hvelmag,.TRUE.,vecpha=phicuratv)

!---wave bottom stress at velocity nodes
IF (iopt_waves.GT.0.AND.iopt_sed_bedeq.EQ.5) THEN
      CALL Carr_at_U(bstresatc_wav_sed(:,1:nrloc),bstresatu_wav,1,3,(/0,1,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_bstresatc_wav_sed,.TRUE.)
      CALL Carr_at_V(bstresatc_wav_sed(1:ncloc,:),bstresatv_wav,1,3,(/1,0,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_bstresatc_wav_sed,.TRUE.)
ENDIF
      
!
!3. Bed load transport formulae
!------------------------------
!
!3.1 Meyer-Peter_Mueller (1948)
!------------------------------
!

IF (iopt_sed_bedeq.EQ.1) THEN

   f_310: DO f=1,nf
!      ---U-nodes
      WHERE (maskatu)
         qbedatu(1:ncloc,:,f) = 8.0*dp(f)*excesshearatu(:,:,f)**1.5&
                              & /taunormatu(:,:,f)
      END WHERE
!      ---V-nodes
      WHERE (maskatv)
         qbedatv(:,1:nrloc,f) = 8.0*dp(f)*excesshearatv(:,:,f)**1.5&
                              & /taunormatv(:,:,f)
      END WHERE
   ENDDO f_310

!
!3.2 Engelund-Fredsoe (1976)
!---------------------------
!

ELSEIF (iopt_sed_bedeq.EQ.2) THEN

   xpi = 0.085*pi
   f_320: DO f = 1,nf
!     ---U-nodes
      WHERE (excesshearatu(:,:,f).GT.0.0)
         peng = excesshearatu(:,:,f)/&
              & (excesshearatu(:,:,f)**4+(xpi*taunormatu(:,:,f))**4)**0.25
         qbedatu(1:ncloc,:,f) = 5.0*dp(f)*peng*&
                             & (SQRT(bstresatu_sed)-0.7*SQRT(taucritatu(:,:,f)))
      END WHERE
!     ---V-nodes
      WHERE (excesshearatv(:,:,f).GT.0.0)
         peng = excesshearatv(:,:,f)/&
              & (excesshearatv(:,:,f)**4+(xpi*taunormatv(:,:,f))**4)**0.25
         qbedatv(:,1:nrloc,f) = 5.0*dp(f)*peng*&
                             & (SQRT(bstresatv_sed)-0.7*SQRT(taucritatv(:,:,f)))
      END WHERE

   ENDDO f_320

!
!3.3 Van Rijn (1984)
!-------------------
!

ELSEIF (iopt_sed_bedeq.EQ.3) THEN


   f_330: DO f=1,nf
!     ---U-nodes
      WHERE (taucritatu(:,:,f).GT.0.0)
         qbedatu(1:ncloc,:,f) = 0.053*qnormatu(:,:,f)*&
                              & (excesshearatu(:,:,f)/taucritatu(:,:,f))**2.1/&
                              & dstaratu(:,:,f)**0.3
      END WHERE
!     ---V-nodes
      WHERE (taucritatv(:,:,f).GT.0.0)
         qbedatv(:,1:nrloc,f) = 0.053*qnormatv(:,:,f)*&
                              & (excesshearatv(:,:,f)/taucritatv(:,:,f))**2.1/&
                              & dstaratv(:,:,f)**0.3
      END WHERE
   ENDDO f_330

!
!3.4 Wu et al. (2000)
!--------------------
!

ELSEIF (iopt_sed_bedeq.EQ.4) THEN

   f_340: DO f = 1,nf
!     ---U-nodes
      WHERE (taucritatu(:,:,f).GT.0.0)
         qbedatu(1:ncloc,:,f) = 0.0053*qnormatu(:,:,f)*&
                              & (excesshearatu(:,:,f)/taucritatu(:,:,f))**2.2
      END WHERE
!     ---V-nodes
      WHERE (taucritatv(:,:,f).GT.0.0)
         qbedatv(:,1:nrloc,f) = 0.0053*qnormatv(:,:,f)*&
                              & (excesshearatv(:,:,f)/taucritatv(:,:,f))**2.2
      END WHERE
   ENDDO f_340

!
!3.5 Soulsby (1997)
!------------------
!

ELSEIF (iopt_sed_bedeq.EQ.5) THEN

!
!3.5.1 U-nodes
!-------------
!
!  ---without waves
   IF (iopt_waves.EQ.0) THEN
      f_3511: DO f=1,nf
         WHERE (maskatu)
            qbedatu(1:ncloc,:,f) = 12.0*SQRT(bstresatu_sed)*&
                    & excesshearatu(:,:,f)/(deltasatu(:,:,f)*gaccatu(1:ncloc,:))
         END WHERE
      ENDDO f_3511

!  ---with waves
   ELSE

!     --angle between currents and waves
      CALL Carr_at_U(wavedir(:,1:nrloc),phiwav,1,3,(/0,1,1/),&
           & (/ncloc,nrloc,1/),1,iarr_wavedir,.TRUE.,angle=.TRUE.)
      WHERE (maskatu)
         phi = phiwav - phicuratu
      END WHERE

!     --bed load per fraction
      f_3512: DO f=1,nf
         WHERE (taucritatu(:,:,f).GT.0.0.AND.bstresatu_sed.GT.taucritatu(:,:,f))
            array2d1 = SQRT(bstresatu_sed)*(bstresatu_sed-taucritatu(:,:,f))
            array2d2 = (0.95+0.19*COS(2.0*phi))*&
                     & SQRT(bstresatu_wav)*bstresatu_sed
            bphix = 12.0*MAX(array2d1,array2d2)
            bphiy = 2.28*bstresatu_sed*bstresatu_wav*bstresatu_wav*&
                  & SIN(2.0*phi)/(bstresatu_wav**1.5+1.5*bstresatu_sed**1.5)
         ELSEWHERE
            bphix = 0.0; bphiy = 0.0
         END WHERE
         WHERE (maskatu)
            array2d1 = taunormatu(:,:,f)/dp(f)
            qbedatu(1:ncloc,:,f) = bphix/array2d1
            qbedatu_per(:,:,f)= bphiy/array2d1
         ELSEWHERE
            qbedatu(1:ncloc,:,f) = 0.0
            qbedatu_per(:,:,f)= 0.0
         END WHERE
      ENDDO f_3512

   ENDIF

!
!3.5.2 V-nodes
!-------------
!
!  ---without waves
   IF (iopt_waves.EQ.0) THEN
      f_3521: DO f=1,nf
         WHERE (maskatv)
            qbedatv(:,1:nrloc,f) = 12.0*SQRT(bstresatv_sed)*&
                    & excesshearatv(:,:,f)/(deltasatv(:,:,f)*gaccatv(:,1:nrloc))
         END WHERE
      ENDDO f_3521

!  ---with waves
   ELSE

!     --angle between currents and waves
      CALL Carr_at_V(wavedir(1:ncloc,:),phiwav,1,3,(/1,0,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_wavedir,.TRUE.,angle=.TRUE.)
      WHERE (maskatv)
         phi = phiwav - phicuratv
      END WHERE

!     --bed load per fraction
      f_3522: DO f = 1,nf 
         WHERE (taucritatv(:,:,f).GT.0.0.AND.bstresatv_sed.GT.taucritatv(:,:,f))
            array2d1 = SQRT(bstresatv_sed)*(bstresatv_sed-taucritatv(:,:,f))
            array2d2 = (0.95+0.19*COS(2.0*phi))*&
                     & SQRT(bstresatv_wav)*bstresatv_sed
            bphix = 12.0*MAX(array2d1,array2d2)
            bphiy = 2.28*bstresatv_sed*bstresatv_wav*bstresatv_wav*&
                  & SIN(2.0*phi)/(bstresatv_wav**1.5+1.5*bstresatv_sed**1.5)
         ELSEWHERE
            bphix = 0.0; bphiy = 0.0
         END WHERE
         WHERE (maskatv)
            array2d1 = taunormatv(:,:,f)/dp(f)
            qbedatv(:,1:nrloc,f) = bphix/array2d1
            qbedatv_per(:,:,f)= bphiy/array2d1
         ELSEWHERE
            qbedatv(:,1:nrloc,f) = 0.0
            qbedatv_per(:,:,f)= 0.0
         END WHERE
      ENDDO f_3522
      
   ENDIF

!
!3.6 Van Rijn (2003/2007)
!------------------------
!

ELSEIF (iopt_sed_bedeq.EQ.6.OR.iopt_sed_bedeq.EQ.7) THEN

! If sand-mud mixtures will be considered in the future, the silt parmeter in
! the van Rijn (2007) equation can be computed. 
! In the van Rijn (2003) equation this parameter is always equal to one

!
!3.6.1 U-nodes
!-------------
!
!3.6.1.1 Without waves
!---------------------
!

   IF (iopt_waves.EQ.0) THEN
   
      f_3611: DO f=1,nf
         WHERE (taucritatu(:,:,f).GT.0.0)
            qbedatu(1:ncloc,:,f) = 0.5*dp(f)*SQRT(bstresatu_sed)*&
                 & excesshearatu(:,:,f)/taucritatu(:,:,f)/dstaratu(:,:,f)**0.3
         END WHERE
      ENDDO f_3611
      qbedatu_per = 0.0

!
!3.6.1.2 With waves
!------------------
!

   ELSE

!     ---angle between currents and waves
      CALL Carr_at_U(wavedir(:,1:nrloc),phiwav,1,3,(/0,1,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_wavedir,.TRUE.,angle=.TRUE.)
      WHERE (maskatu)
         phi = phiwav - phicuratu
      END WHERE

!     ---wave excursion and velocity at U-nodes      
      CALL Carr_at_U(waveexcurs(:,1:nrloc),awave,1,3,(/0,1,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_waveexcurs,.TRUE.)
      CALL Carr_at_U(wavevel(:,1:nrloc),wavevelatu,1,3,(/0,1,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_wavevel,.TRUE.)

!     ---wave friction factor
      CALL Carr_at_U(fwave_sed(:,1:nrloc),fw,1,3,(/0,1,1/),(/ncloc,nrloc,1/),1,&
                   & iarr_fwave,.TRUE.)

!     ---apparent grain roughness
      WHERE (maskatu)
         array2d1 = ATAN2(ABS(SIN(phi)),COS(phi))
         array2d1 = 0.8 + array2d1 - 0.3*array2d1*array2d1
      END WHERE

      CALL vector_mag_arr_atu(uvel(1:ncloc,1:nrloc,1),&
                            & vvel(0:ncloc,1:nrloc+1,1),&
                            & 1,3,1,1,iarr_hvelmag,.TRUE.,vecmag=curmag)
      ks = 1.5*d50atu
      WHERE (curmag.GT.0.0)
         ka = MAX(ks,0.01)*MIN(EXP(array2d1*wavevelatu/curmag),10.0)
      ELSEWHERE
         ka = MAX(ks,0.01)
      END WHERE

!     ---boundary layer thickness
      WHERE (awave.GT.0.0)
         delta = 0.216*awave**0.75*ks**0.25
         delta = MAX(30.0*delta,30.0*zroughatu_sed)
      ELSEWHERE
         delta = 30.0*zroughatu_sed
      END WHERE

!     ---(grain) friction factor due to currents and waves
      WHERE (maskatu)
         beta = 0.25*((LOG(15.0*deptotatu(1:ncloc,1:nrloc)/ks)-1.0)&
                    & /LOG(30.0*delta/ks))**2
         uc = curmag*LOG(30.0*delta/ka)&
                 & /(LOG(30.0*deptotatu(1:ncloc,1:nrloc)/ka-1.0))
         alpha = uc/(uc+wavevelatu)
         array2d1 =  2.0*(ckar/LOG(10.0*deptotatu(1:ncloc,1:nrloc)/d50atu))**2
         fcw = alpha*beta*array2d1 + (1.0-alpha)*fw
      END WHERE

!     ---near-bed orbital velocity over wave cycle
      n_3612: DO n=1,nrquad_wav
         WHERE (maskatu)
            array2d1 = (wavevelatu*SIN(ndgau_tr(n)))**2+&
                      & 2.0*wavevelatu*uc*SIN(ndgau_tr(n))*COS(phi)
            velcurwav(:,:,n) = SQRT(array2d1+uc**2)
         ELSEWHERE
            velcurwav(:,:,n) = 0.0
         END WHERE
      ENDDO n_3612

!     ---transport per fraction by (Gaussian) integration over wave period
      f_3613: DO f=1,nf
!        --transport at specified phases during wave cycle
         WHERE (maskatu)
            qbedatu(1:ncloc,:,f) = 0.0; qbedatu_per(:,:,f) = 0.0
         END WHERE
         n_36131: DO n=1,nrquad_wav
            WHERE (taucritatu(:,:,f).GT.0.0.AND.velcurwav(:,:,n).GT.0.0)
               array2d1 = 0.5*fcw*velcurwav(:,:,n)**2
               oscwavetr = 0.5*dp(f)*dstaratu(:,:,f)**(-0.3)*SQRT(array2d1)*&
                         & MAX((array2d1/taucritatu(:,:,f)-1.0),0.0)
               array2d1 = (uc+wavevelatu*SIN(ndgau_tr(n))*COS(phi))&
                        & /velcurwav(:,:,n)
               oscwavetr = array2d1*oscwavetr
               array2d1 = SIN(ndgau_tr(n))*wavevelatu*SIN(phi)/velcurwav(:,:,n)
               oscwavetr_per = array2d1*oscwavetr
               qbedatu(1:ncloc,:,f) = qbedatu(1:ncloc,:,f) + &
                                    & 0.5*weights_gauss(n)*oscwavetr
               qbedatu_per(:,:,f) = qbedatu_per(:,:,f) + &
                                  & 0.5*weights_gauss(n)*oscwavetr_per
            END WHERE
         ENDDO n_36131
      ENDDO f_3613

   ENDIF

!
!3.6.2 V-nodes
!-------------
!
!3.6.2.1 Without waves
!---------------------
!

   IF (iopt_waves.EQ.0) THEN
   
      f_3621: DO f=1,nf
         WHERE (taucritatv(:,:,f).GT.0.0)
            qbedatv(:,1:nrloc,f) = 0.5*dp(f)*SQRT(bstresatv_sed)*&
                 & excesshearatv(:,:,f)/taucritatv(:,:,f)/dstaratv(:,:,f)**0.3
         END WHERE
      ENDDO f_3621
      qbedatv_per = 0.0

!
!3.6.1.2 With waves
!------------------
!

   ELSE

!     ---angle between currents and waves
      CALL Carr_at_V(wavedir(1:ncloc,:),phiwav,1,3,(/1,0,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_wavedir,.TRUE.,angle=.TRUE.)
      WHERE (maskatv)
         phi = phiwav - phicuratv
      END WHERE

!     ---wave excursion and velocity at V-nodes      
      CALL Carr_at_V(waveexcurs(1:ncloc,:),awave,1,3,(/1,0,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_waveexcurs,.TRUE.)
      CALL Carr_at_V(wavevel(1:ncloc,:),wavevelatv,1,3,(/1,0,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_wavevel,.TRUE.)

!     ---wave friction factor
      CALL Carr_at_V(fwave_sed(1:ncloc,:),fw,1,3,(/1,0,1/),(/ncloc,nrloc,1/),1,&
                   & iarr_fwave,.TRUE.)
      
!     ---apparent grain roughness
      WHERE (maskatv)
         array2d1 = ATAN2(ABS(SIN(phi)),COS(phi))
         array2d1 = 0.8 + array2d1 - 0.3*array2d1*array2d1
      END WHERE

      CALL vector_mag_arr_atv(uvel(1:ncloc+1,0:nrloc,1),&
                            & vvel(1:ncloc,1:nrloc,1),&
                            & 1,3,1,1,iarr_hvelmag,.TRUE.,vecmag=curmag)
      ks = 1.5*d50atv
      WHERE (curmag.GT.0.0 )
         ka = MAX(ks,0.01)*MIN(EXP(array2d1*wavevelatv/curmag),10.0)
      ELSEWHERE
         ka = MAX(ks,0.01)
      END WHERE

!     ---boundary layer thickness
      WHERE (awave.GT.0.0)
         delta = 0.216*awave**0.75*ks**0.25
         delta = MAX(30.0*delta,30.0*zroughatv_sed)
      ELSEWHERE
         delta = 30.0*zroughatv_sed
      END WHERE

!     ---(grain) friction factor due to currents and waves
      WHERE (maskatv)
         beta = 0.25*((LOG(15.0*deptotatv(1:ncloc,1:nrloc)/ks)-1.0)&
                    & /LOG(delta/ks))**2
         uc = curmag*LOG(30.0*delta/ka)/&
                  & (LOG(30.0*deptotatv(1:ncloc,1:nrloc)/ka-1.0))
         alpha = uc/(uc+wavevelatv)
         array2d1 =  2.0*(ckar/LOG(10.0*deptotatv(1:ncloc,1:nrloc)/d50atv))**2
         fcw = alpha*beta*array2d1 + (1.0-alpha)*fw
      END WHERE

!     ---near-bed orbital velocity over wave cycle
      n_3622: DO n=1,nrquad_wav
         WHERE (maskatv)
            array2d1 = (wavevelatv*SIN(ndgau_tr(n)))**2+&
                      & 2.0*wavevelatv*uc*SIN(ndgau_tr(n))*COS(phi)
            velcurwav(:,:,n) = SQRT(array2d1+uc**2)
         ELSEWHERE
            velcurwav(:,:,n) = 0.0
         END WHERE
      ENDDO n_3622

!     ---transport per fraction by (Gaussian) integration over wave period
      f_3623: DO f=1,nf
!        --transport at specified phases during wave cycle
         WHERE (maskatv)
            qbedatv(:,1:nrloc,f) = 0.0; qbedatv_per(:,:,f) = 0.0
         END WHERE
         n_36231: DO n=1,nrquad_wav
            WHERE (taucritatv(:,:,f).GT.0.0.AND.velcurwav(:,:,n).GT.0.0)
               array2d1 = 0.5*fcw*velcurwav(:,:,n)**2
               oscwavetr = 0.5*dp(f)*dstaratv(:,:,f)**(-0.3)*SQRT(array2d1)*&
                         & MAX((array2d1/taucritatv(:,:,f)-1.0),0.0)
               array2d1 = (uc+wavevelatv*SIN(ndgau_tr(n))*COS(phi))&
                        & /velcurwav(:,:,n)
               oscwavetr = array2d1*oscwavetr
               array2d1 = SIN(ndgau_tr(n))*wavevelatv*SIN(phi)/velcurwav(:,:,n)
               oscwavetr_per = array2d1*oscwavetr
               qbedatv(:,1:nrloc,f) = qbedatv(:,1:nrloc,f) + &
                                    & 0.5*weights_gauss(n)*oscwavetr
               qbedatv_per(:,:,f) = qbedatv_per(:,:,f) + &
                                  & 0.5*weights_gauss(n)*oscwavetr_per
            END WHERE
         ENDDO n_36231
      ENDDO f_3623

   ENDIF

ENDIF

!
!4. Bed slope effects and projection to coordinate directions
!------------------------------------------------------------
!

flag = iopt_waves.GT.0.AND.iopt_sed_bedeq.GE.5.AND.iopt_sed_bedeq.LE.7

!
!4.1 U-nodes
!-----------
!

maskatu = bstresatu_sed.GT.0.0

WHERE (maskatu)
   qprojx = COS(phicuratu)
   qprojy = SIN(phicuratu)
END WHERE

!
!4.1.1 Without slope
!-------------------
!

IF (iopt_sed_slope.EQ.0) THEN

   f_411: DO f=1,nf
!     ---longitudinal component
      WHERE (maskatu)
         qbedatu(1:ncloc,:,f) = qprojx*qbedatu(1:ncloc,:,f)
      END WHERE
!     ---transverse component
      IF (flag) THEN
         WHERE (maskatu)
            qbedatu(1:ncloc,:,f) = qbedatu(1:ncloc,:,f) - &
                                 & qprojy*qbedatu_per(:,:,f)
         END WHERE
      ENDIF
   ENDDO f_411

!
!4.1.2 With slope
!----------------
!

ELSE

!  ---bed slope factors
   WHERE (maskatu)
      xslope = 1.0 - coef_bed_grad*(bed_slope_x_atu*qprojx+&
                                  & bed_slope_y_atu*qprojy)
      yslope = coef_bed_grad*(bed_slope_y_atu*qprojx-bed_slope_x_atu*qprojy)
      array2d1 = SQRT(1.0+yslope*yslope)
      array2d2 = 1.0/array2d1
      array2d3 = yslope/array2d1
   END WHERE
  
   f_412: DO f=1,nf
!     ---longitudinal component
      WHERE (maskatu)
         qbedatu(1:ncloc,:,f) = qbedatu(1:ncloc,:,f)*xslope*&
                             & (qprojx*array2d2-qprojy*array2d3)
      END WHERE
!     ---transverse component
      IF (flag) THEN
         WHERE (maskatu)
            qbedatu(1:ncloc,:,f) = qbedatu(1:ncloc,:,f)-qbedatu_per(:,:,f)*&
                                 & xslope*(qprojy*array2d2+qprojx*array2d3)
         END WHERE
      ENDIF
   ENDDO f_412
   
ENDIF

!
!4.2 V-nodes
!-----------
!

maskatv = bstresatv_sed.GT.0.0

WHERE (maskatv)
   qprojx = COS(phicuratv) 
   qprojy = SIN(phicuratv)
END WHERE

!
!4.2.1 Without slope
!-------------------
!

IF (iopt_sed_slope.EQ.0) THEN

   f_421: DO f=1,nf
!     ---longitudinal component
      WHERE (maskatv)
         qbedatv(:,1:nrloc,f) = qprojy*qbedatv(:,1:nrloc,f)
      END WHERE
!     ---transverse component
      IF (flag) THEN
         WHERE (maskatv)
            qbedatv(:,1:nrloc,f) = qbedatv(:,1:nrloc,f) + &
                                 & qprojx*qbedatv_per(:,:,f)
         END WHERE
      ENDIF
   ENDDO f_421

!
!4.2.2 With slope
!----------------
!

ELSE

!  ---bed slope factors
   WHERE (maskatv)
      xslope = 1.0 - coef_bed_grad*(bed_slope_x_atv*qprojx+&
                                  & bed_slope_y_atv*qprojy)
      yslope = coef_bed_grad*(bed_slope_y_atv*qprojx-bed_slope_x_atv*qprojy)
      array2d1 = SQRT(1.0+yslope*yslope)
      array2d2 = 1.0/array2d1
      array2d3 = yslope/array2d1
   END WHERE
  
   f_422: DO f=1,nf
!     ---longitudinal component
      WHERE (maskatv)
         qbedatv(:,1:nrloc,f) = qbedatv(:,1:nrloc,f)*xslope*&
                             & (qprojy*array2d2+qprojx*array2d3)
      END WHERE
!     ---transverse component
      IF (flag) THEN
         WHERE (maskatv)
            qbedatv(:,1:nrloc,f) = qbedatv(:,1:nrloc,f)+qbedatv_per(:,:,f)*&
                                 & xslope*(qprojx*array2d2-qprojx*array2d3)
         END WHERE
      ENDIF
   ENDDO f_422
   
ENDIF

!
!5. Fractional transport
!-----------------------
!
!5.1 U-nodes
!-----------
!

WHERE (.NOT.maskatu) array2d1 = 0.0

f_510: DO f=1,nf
   WHERE (maskatu)
      array2d1 = SIGN(1.0,qbedatu(1:ncloc,:,f))
   END WHERE
   WHERE (array2d1.GT.0.0)
      qbedatu(1:ncloc,:,f) = bed_fraction(0:ncloc-1,1:nrloc,1,f)*&
                           & qbedatu(1:ncloc,:,f)
   ELSEWHERE
      qbedatu(1:ncloc,:,f) = bed_fraction(1:ncloc,1:nrloc,1,f)*&
                           & qbedatu(1:ncloc,:,f)
   END WHERE
ENDDO f_510

!
!5.2 V-nodes
!-----------
!

WHERE (.NOT.maskatv) array2d1 = 0.0

f_520: DO f=1,nf
   WHERE (maskatv)
      array2d1 = SIGN(1.0,qbedatv(:,1:nrloc,f))
   END WHERE
   WHERE (array2d1.GT.0.0)
      qbedatv(:,1:nrloc,f) = bed_fraction(1:ncloc,0:nrloc-1,1,f)*&
                           & qbedatv(:,1:nrloc,f)
   ELSEWHERE
      qbedatv(:,1:nrloc,f) = bed_fraction(1:ncloc,1:nrloc,1,f)*&
                           & qbedatv(:,1:nrloc,f)
   END WHERE
ENDDO f_520

!
!6. Exchange halos
!-----------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = (/0,1,1/); nhexch = (/1,1,0,0/)
   CALL exchange_mod(qbedatu,lbounds,nhexch,iarr_qbedatu)
   lbounds = (/1,0,1/); nhexch = (/0,0,1,1/)
   CALL exchange_mod(qbedatv,lbounds,nhexch,iarr_qbedatv)
ENDIF

!
!7. Upwind correction at open boundaries
!---------------------------------------
!

IF (iopt_sed_obc_flux.EQ.1) THEN

!
!7.1 U-nodes
!-----------
!

   iiloc_710: DO iiloc=1,nobuloc
      i = iobuloc(iiloc); j = jobuloc(iiloc)
      ii = indexobu(iiloc); iu = MERGE(i+1,i-1,westobu(ii))
      qbedatu(i,j,:) = (delyatu(iu,j)/delyatu(i,j))*qbedatu(iu,j,:)
   ENDDO iiloc_710
 
!
!7.2 V-nodes
!-------------
!

   jjloc_720: DO jjloc=1,nobvloc
      i = iobvloc(jjloc); j = jobvloc(jjloc)
      jj = indexobv(jjloc); jv = MERGE(j+1,j-1,soutobv(jj))
      qbedatv(i,j,:) = (delxatv(i,jv)/delxatv(i,j))*qbedatv(i,jv,:)
   ENDDO jjloc_720

!
!7.3 Exchange halos
!------------------
!

   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/0,1,1/); nhexch = (/1,1,0,0/)
      CALL exchange_mod(qbedatu,lbounds,nhexch,iarr_qbedatu)
      lbounds = (/1,0,1/); nhexch = (/0,0,1,1/)
      CALL exchange_mod(qbedatv,lbounds,nhexch,iarr_qbedatv)
   ENDIF
   
ENDIF

!
!8. Bed load averaged over morphological step
!--------------------------------------------
!

IF (iopt_morph.EQ.1) THEN
   morphsteps = icmorph/icsed
   f_810: DO f=1,nf
      IF (icmorph.EQ.icsed) THEN
         qbedatu_int(:,:,f) = qbedatu(:,:,f)
         qbedatv_int(:,:,f) = qbedatv(:,:,f)
      ELSE
         IF (mult_index(nt-1,icmorph)) THEN
            qbedatu_int(:,:,f) = 0.0; qbedatv_int(:,:,f) = 0.0
         ENDIF
         WHERE (node2du(0:ncloc+1,1:nrloc).GT.1)
            qbedatu_int(:,:,f) = qbedatu_int(:,:,f) + qbedatu(:,:,f)/morphsteps
         END WHERE
         WHERE (node2dv(1:ncloc,0:nrloc+1).GT.1)
            qbedatv_int(:,:,f) = qbedatv_int(:,:,f) + qbedatv(:,:,f)/morphsteps
         END WHERE
      ENDIF
   ENDDO f_810
ENDIF

!
!9. Deallocate arrays
!--------------------
!

IF ((cold_start.OR.nt.EQ.nstep).AND.iopt_waves.EQ.1.AND.&
   & iopt_sed_bedeq.GT.5) THEN
   DEALLOCATE (ndgau_tr,nodes_gauss,weights_gauss)
ENDIF

DEALLOCATE (maskatu,maskatv)
DEALLOCATE (array2d1,array2d2,array2d3,deltasatu,deltasatv,excesshearatu,&
          & excesshearatv,phicuratu,phicuratv,qprojx,qprojy,taucritatu,&
          & taucritatv)
IF (iopt_sed_bedeq.EQ.1.OR.iopt_sed_bedeq.EQ.2.OR.iopt_sed_bedeq.EQ.5) THEN
   DEALLOCATE (taunormatu,taunormatv)
ENDIF
IF (iopt_sed_bedeq.EQ.3.OR.iopt_sed_bedeq.EQ.4) THEN
   DEALLOCATE (qnormatu,qnormatv)
ENDIF
IF (iopt_sed_bedeq.EQ.2) DEALLOCATE (peng)
IF (iopt_sed_bedeq.EQ.3.OR.iopt_sed_bedeq.EQ.6.OR.iopt_sed_bedeq.EQ.7) THEN
   DEALLOCATE (dstaratu,dstaratv)
ENDIF
IF (iopt_sed_bedeq.EQ.5) THEN
   DEALLOCATE (bphix,bphiy,bstresatu_wav,bstresatv_wav)
ENDIF
IF (iopt_sed_bedeq.GE.5) THEN
   DEALLOCATE (d50atu,d50atv,phi,phiwav,qbedatu_per,qbedatv_per)
ENDIF
IF (iopt_sed_bedeq.GE.6) THEN
   DEALLOCATE (alpha,amplitude,awave,beta,curmag,delta,fcw,fw,ka,ks,&
             & oscwavetr,oscwavetr_per,uc,velcurwav,wavevelatu,wavevelatv)
ENDIF
IF (iopt_sed_slope.EQ.1) DEALLOCATE (xslope,yslope)

CALL log_timer_out()


RETURN

END SUBROUTINE sediment_bedload

!========================================================================

SUBROUTINE sediment_equation
!************************************************************************
!
! *sediment_equation* Main sediment transport model
!
! Author - Alexander Breugem and Boudewijn Decrop (IMDC)
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.1
!
! Description -
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! External calls - bed_slope_arrays, bottom_stress_sed, critical_shear_stress, 
!                  flocculation_model, kinematic_viscosity, sediment_advdiff,
!                  sediment_bedload, sediment_totalload
!
! Module calls - 
!
!************************************************************************
!
USE iopars
USe morphswitches
USE sedpars
USE sedswitches
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc

sedstep = (nt/icsed)*icsed.EQ.nt
IF (.NOT.sedstep) RETURN

procname(pglev+1) = 'sediment_equation'
CALL log_timer_in(npcc)

!
!1. Bottom slopes
!-----------------
!

IF ((iopt_sed_slope.GT.0).OR.(iopt_morph_avalanching.GT.0)) THEN
   CALL bed_slope_arrays
ENDIF

!
!2. Kinematic viscosity
!----------------------
!

CALL kinematic_viscosity

!
!3. Bottom stress arrays
!-----------------------
!

CALL bottom_stress_sed

!
!4. Critical shear stress
!------------------------
!

CALL  critical_shear_stress

!
!5. Sediment transport
!---------------------
!
!---bed load
IF (iopt_sed_mode.EQ.1.OR.iopt_sed_mode.EQ.3) THEN
   CALL sediment_bedload
!---total load
ELSEIF (iopt_sed_mode.EQ.4) THEN
   CALL sediment_totalload
ENDIF

!---without flocculation
IF (iopt_sed.EQ.1) THEN
   IF (iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3) CALL sediment_advdiff
!---flocculation module      
ELSEIF (iopt_sed.EQ.2) THEN
   CALL flocculation_model
ENDIF

CALL log_timer_out(npcc,itm_sed)


RETURN

END SUBROUTINE sediment_equation

!========================================================================

SUBROUTINE sediment_suspendedload
!************************************************************************
!
! *sediment_suspendedload* Calculates suspended load magntide at C-nodes
!
! Author - IMDC, Alexander Breugem
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - equilibrium_concentration, sediment_equation
!
! External calls - ackerswhite_params, reference_concentration,
!                  thetastar_engelund_hansen
!
! Module calls - error_alloc, vector_mag_arr_atc
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
USE sedarrays
USE sedpars
USE sedswitches
USE switches
USE syspars
USE wavevars
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: vector_mag_arr_atc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables 
!
INTEGER :: itype, f
REAL :: onethird = 1.0/3.0, real_small = 1.0E-04, wmin = 0.03, &
      & ws_ust_max = 100.0
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2, aw, curmag, &
         & cw, dstar, fgr, ggr, hcnorm, kw, mw, nw, qnorm, theta, thetastar, &
         & wavelen, xterm1, xterm2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: deltas, taunorm


procname(pglev+1) = 'sediment_suspendedload'
CALL log_timer_in()

!
! 1.Initialise arrays
!--------------------
!
!1.1 Allocate
!------------
!

IF (iopt_sed_bbc_eq.EQ.2.OR.iopt_sed_bbc_eq.EQ.3) THEN
   ALLOCATE(theta(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('theta',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE(thetastar(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('thetastar',2,(/ncloc,nrloc/),kndrtype)
ELSEIF (iopt_sed_bbc_eq.EQ.4) THEN
   ALLOCATE(dstar(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('dstar',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE(fgr(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fgr',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE(ggr(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('ggr',2,(/ncloc,nrloc/),kndrtype)
ELSEIF (iopt_sed_bbc_eq.EQ.5) THEN
   ALLOCATE(array2d2(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE(hcnorm(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('hcnorm',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE(xterm1(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('xterm1',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE(xterm2(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('xterm2',2,(/ncloc,nrloc/),kndrtype)
   IF (iopt_waves.EQ.1) THEN
      ALLOCATE(kw(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('kw',2,(/ncloc,nrloc/),kndrtype)
      ALLOCATE(wavelen(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('wavelen',2,(/ncloc,nrloc/),kndrtype)
   ENDIF
ENDIF
IF (iopt_sed_bbc_eq.NE.5) THEN
   ALLOCATE(deltas(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('deltas',3,(/ncloc,nrloc,nf/),kndrtype)
ENDIF
IF (iopt_sed_bbc_eq.LT.5) THEN
   ALLOCATE(taunorm(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('taunorm',3,(/ncloc,nrloc,nf/),kndrtype)
ENDIF
IF (iopt_sed_bbc_eq.EQ.2.OR.iopt_sed_bbc_eq.EQ.3.OR.iopt_sed_bbc_eq.EQ.6) THEN
   ALLOCATE(qnorm(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('qnorm',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (iopt_sed_bbc_eq.EQ.4.OR.iopt_sed_bbc_eq.EQ.5) THEN
   ALLOCATE(mw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('mw',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (iopt_sed_bbc_eq.EQ.5.OR.iopt_sed_bbc_eq.EQ.6) THEN
   ALLOCATE(array2d1(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (iopt_sed_bbc_eq.GE.4.AND.iopt_sed_bbc_eq.LE.6) THEN
   ALLOCATE(curmag(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('curmag',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (iopt_sed_bbc_eq.EQ.4.OR.(iopt_sed_bbc_eq.EQ.5.AND.iopt_waves.EQ.1)) THEN
   ALLOCATE(aw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('aw',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE(cw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('cw',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE(nw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('nw',2,(/ncloc,nrloc/),kndrtype)
ENDIF

!
!1.2 Initialise
!--------------
!

f_120: DO f = 1,nf
   IF (iopt_sed_bbc_eq.NE.5) THEN
      WHERE (maskatc_int)
         deltas(:,:,f)  = MAX(rhos(f)/densw(1:ncloc,1:nrloc,1)-1.0,0.0)
      END WHERE
   ENDIF
   IF (iopt_sed_bbc_eq.LT.5) THEN
      WHERE (maskatc_int)
         taunorm(:,:,f) = deltas(:,:,f)*gaccatc(1:ncloc,1:nrloc)*dp(f)
      END WHERE
   ENDIF
 ENDDO f_120

!
!2.Transport formulae
!--------------------
!
!2.1 Engelund-Hansen
!-------------------
!

IF (iopt_sed_bbc_eq.EQ.2.OR.iopt_sed_bbc_eq.EQ.3) THEN

   itype = MERGE(1,2,iopt_sed_bbc_eq.EQ.2)
   f_210: DO f = 1,nf
     WHERE (maskatc_int)
        deltas(:,:,f)  = MAX(rhos(f)/densw(1:ncloc,1:nrloc,1)-1.0,0.0)
        theta = bstresatc_sed(1:ncloc,1:nrloc)/taunorm(:,:,f)
     END WHERE
     CALL thetastar_engelund_hansen(theta,maskatc_int,thetastar,itype)
     WHERE (maskatc_int)
        qnorm = SQRT(deltas(:,:,f)*gaccatc(1:ncloc,1:nrloc)*dp(f)**3)
        qsusatc(:,:,f) = 0.05*qnorm*thetastar**(2.5)/&
                       & bdragcoefatc_sed(1:ncloc,1:nrloc)
     END WHERE
  ENDDO f_210

!
!4.2 Ackers-White
!----------------
!

ELSEIF (iopt_sed_bbc_eq.EQ.4) THEN

   CALL vector_mag_arr_atc(uvel(1:ncloc+1,1:nrloc,1),vvel(1:ncloc,1:nrloc+1,1),&
                         & 1,1,1,1,iarr_hvelmag,.TRUE.,vecmag=curmag)

   f_420: DO f=1,nf
      WHERE (maskatc_int)
         dstar = dp(f)*(deltas(:,:,f)*gaccatc(1:ncloc,1:nrloc)&
               & /kinvisc(1:ncloc,1:nrloc,1)**2)**onethird
      END WHERE
      CALL ackerswhite_params(dstar,maskatc_int,nw,mw,aw,cw)
      WHERE (maskatc_int)
         fgr = SQRT(bstresatc_sed(1:ncloc,1:nrloc)**nw/taunorm(:,:,f))*&
          & (curmag/(2.46*LOG(10.0*deptotatc(1:ncloc,1:nrloc)/dp(f))))**(1.0-nw)
      ELSEWHERE
         fgr = 0.0
      END WHERE
      WHERE (maskatc_int.AND.fgr.GT.aw)
         ggr = cw*(fgr/aw-1.0)**mw
         qsusatc(:,:,f) =  dp(f)*curmag*ggr&
                         & /SQRT((bdragcoefatc_sed(1:ncloc,1:nrloc))**nw)
      ELSEWHERE
         qsusatc(:,:,f) = 0.0
      END WHERE
   END DO f_420

!
!4.3 Van Rijn (2003)
!-------------------
!
! Is similar to Van Rijn (1984), except that 
!    1.) Wave input is accounted for
!    2.) An explicit expression is used for turbulence damping due to sediment

ELSEIF (iopt_sed_bbc_eq.EQ.5) THEN

!
!4.3.1 Current magnitude
!-----------------------
!

   CALL vector_mag_arr_atc(uvel(1:ncloc+1,1:nrloc,1),vvel(1:ncloc,1:nrloc+1,1),&
                         & 1,1,1,1,iarr_hvelmag,.TRUE.,vecmag=curmag)

!
!4.3.2 Bed boundary condition
!----------------------------
!

   CALL reference_concentration

!
!4.3.3 Transport in wave direction (Grasmeijer & van Rijn 1998)
!--------------------------------------------------------------
!

   IF (iopt_waves.EQ.1) THEN

      WHERE (maskatc_int)
         wavelen = twopi/wavenum(1:ncloc,1:nrloc)
         array2d1 = wavelen*waveheight(1:ncloc,1:nrloc)
         array2d2 = wavevel(1:ncloc,1:nrloc)*&
             & (2.0-6.4*array2d1**0.65*array2d1**&
             & (3.4*wavelen*deptotatc(1:ncloc,1:nrloc)))
!        ---procedure to find onshore orbital velocity peak
         aw = waveperiod(1:ncloc,1:nrloc)*&
            & SQRT(gaccatc(1:ncloc,1:nrloc)/deptotatc(1:ncloc,1:nrloc))
         mw = MERGE(3.2E-3*aw**2+8E-5*aw**3,5.6E-3*aw**2+4E-5*aw**3,aw.LE.30)!L5
         aw = MERGE(-15.0+1.35*aw,-2.7+0.53*aw, aw.LE.15.0)               !L4
         nw = (0.5-mw)/(aw-1+EXP(-aw))                                    !L3
         mw = nw*aw + mw                                                  !L2
         cw = array2d2/SQRT(gaccatc(1:ncloc,1:nrloc)*deptotatc(1:ncloc,1:nrloc))
         cw = array2d2*(0.5-nw+mw*cw+nw*EXP(-aw*cw))                      !U_on
         mw = SQRT(curmag**2+cw**2)                                       !v_eff
      ELSEWHERE
         array2d2 = 0.0
      END WHERE
!     ---velocity asymmetry
      WHERE (array2d2.NE.0.0)
         nw = (cw**4-(array2d2-cw)**4)/(cw**3+(array2d2-cw)**3)
      ELSEWHERE
         nw = 0.0
      END WHERE

   ENDIF

!
!4.3.4 Transport
!---------------
!

   f_434: DO f=1,nf

      WHERE (maskatc_int)
         qsusatc(:,:,f) = deptotatc(1:ncloc,1:nrloc)*curmag*&
                        & cref(1:ncloc,1:nrloc,f)
      END WHERE

!     ---ws/u*
      WHERE (bstresatc_sed(1:ncloc,1:nrloc).GT.0.0)
         array2d1 = wfall(1:ncloc,1:nrloc,1,f)&
                  & /SQRT(bstresatc_sed(1:ncloc,1:nrloc))
      ELSEWHERE
         array2d1 = ws_ust_max
      END WHERE
  
      WHERE (maskatc_int)
!        ---beta
         array2d2 = MAX(beta_sed_min,MIN(1.0+2.0*(array2d1**2),beta_sed_max))
!        ---modified Rouse number
         array2d2 = array2d1/(array2d2*ckar) + &
                  & 2.5*array2d1**0.8*(cref(1:ncloc,1:nrloc,f)/cmax)**0.4
!        ---a/h
         hcnorm = height_c(1:ncloc,1:nrloc,f)/deptotatc(1:ncloc,1:nrloc)
!        ---prevent numerical problems for high Rouse numbers
!        ---location of maximum of F(Z') as function of a/h (fitted numerically)
         array2d2 = MIN(array2d2,1.1821*hcnorm**(-0.96038))
      END WHERE
!     ---factor F (current and waves)
      WHERE (maskatc_int.AND.(ABS(array2d2-1.2).GT.real_small))
         xterm1 = MAX((hcnorm**array2d2-hcnorm**1.2)/&
                & ((1-hcnorm)**array2d2*(1.2-array2d2)),0.0)
      ELSEWHERE
         xterm1 = 0.0
      END WHERE

      IF (iopt_waves.EQ.1) THEN
         WHERE (waveheight(1:ncloc,1:nrloc).GT.0.0)
            array2d1 = 4.0*(0.2*deptotatc(1:ncloc,1:nrloc))**0.6 *&
               & (MAX(wfall(1:ncloc,1:nrloc,1,f),wmin)*&
               & waveperiod(1:ncloc,1:nrloc)/waveheight(1:ncloc,1:nrloc))**0.8
         ELSEWHERE
            array2d1 = 0.0
         END WHERE
         WHERE (maskatc_int.AND.(ABS(array2d1-1.2).GT.real_small))
            xterm2 = MAX((hcnorm**array2d1-hcnorm**1.2)/&
                    & ((1.0-hcnorm)**array2d1*(1.2-array2d1)),0.0)
         ELSEWHERE
            xterm2 = 0.1649
         END WHERE
!        ---transport in wave direction     
         WHERE (maskatc_int)
            kw = LOG(2.0*deptotatc(1:ncloc,1:nrloc)/dp(f))
!           --v_crit
            kw = MERGE(0.19*dp(f)**0.1*kw,8.5*dp(f)**0.6*kw,dp(f).LT.5E-4)   
!           --sediment mobility number
            mw = (mw-kw)**2/taunorm(:,:,f)      
!           --suspended transport in wave direction [m2/s]
            mw = 0.0014*dp(f)*mw*nw
         END WHERE
      ELSE
         xterm2 = 0.0; mw = 0.0
      ENDIF

!     ---transport
      WHERE (maskatc_int)
         qsusatc(:,:,f) = qsusatc(1:ncloc,1:nrloc,f)*(xterm1+xterm2) + mw
      END WHERE

   ENDDO f_434

!
!4.4 Wu (2000)
!-------------
!

ELSEIF (iopt_sed_bbc_eq.EQ.6) THEN

   CALL vector_mag_arr_atc(uvel(1:ncloc+1,1:nrloc,1),vvel(1:ncloc,1:nrloc+1,1),&
                         & 1,1,1,1,iarr_hvelmag,.TRUE.,vecmag=curmag)

   f_440: DO f=1,nf
      WHERE (maskatc_int.AND.bstresatc_sed(1:ncloc,1:nrloc).GT.&
           & bstres_cr(1:ncloc,1:nrloc,f))
         array2d1 = 2.62E-05*((bstresatc_sed(1:ncloc,1:nrloc)&
                            & /bstres_cr(1:ncloc,1:nrloc,f)-1.0)*&
                            & (curmag/wfall(1:ncloc,1:nrloc,1,f)))**1.74
         qnorm = SQRT(deltas(:,:,f)*gaccatc(1:ncloc,1:nrloc)*dp(f)**3)
         qsusatc(:,:,f) = array2d1*qnorm
      ELSEWHERE
         qsusatc(:,:,f) = 0.0
      END WHERE
   ENDDO f_440

ENDIF

!
!5. Fractional transport
!-----------------------
!

f_510: DO f=1,nf
   WHERE (maskatc_int)
      qsusatc(:,:,f) = qsusatc(:,:,f)*bed_fraction(1:ncloc,1:nrloc,1,f)
   END WHERE
ENDDO f_510

!
!6. Deallocate
!-------------
!

IF (iopt_sed_bbc_eq.EQ.2.OR.iopt_sed_bbc_eq.EQ.3) THEN
   DEALLOCATE (theta,thetastar)
ELSEIF (iopt_sed_bbc_eq.EQ.4) THEN
   DEALLOCATE (dstar,fgr,ggr)
ELSEIF (iopt_sed_bbc_eq.EQ.5) THEN
   DEALLOCATE (array2d2,hcnorm,xterm1,xterm2)
   IF (iopt_waves.EQ.1) DEALLOCATE (kw,wavelen)
ENDIF
IF (iopt_sed_bbc_eq.NE.5) DEALLOCATE (deltas)
IF (iopt_sed_bbc_eq.LT.5) DEALLOCATE (taunorm)
IF (iopt_sed_bbc_EQ.EQ.2.OR.iopt_sed_bbc_eq.EQ.3.OR.&
     & iopt_sed_bbc_eq.EQ.6) THEN
   DEALLOCATE (qnorm)
ENDIF
IF (iopt_sed_bbc_eq.EQ.4.OR.iopt_sed_bbc_eq.EQ.5) DEALLOCATE (mw)
IF (iopt_sed_bbc_eq.EQ.5.OR.iopt_sed_bbc_eq.EQ.6) DEALLOCATE (array2d1)
IF (iopt_sed_bbc_eq.GE.4.AND.iopt_sed_bbc_eq.LE.6) THEN
   DEALLOCATE (curmag)
ENDIF
IF (iopt_sed_bbc_eq.EQ.4.OR.(iopt_sed_bbc_eq.EQ.5.AND.iopt_waves.EQ.1)) THEN
   DEALLOCATE (aw,cw,nw)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE sediment_suspendedload

!========================================================================

SUBROUTINE sediment_totalload
!************************************************************************
!
! *sediment_totalload* Solve sediment total load formulae
!
! Author - Boudewijn Decrop (IMDC)
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - sediment_equation
!
! External calls - ackerswhite_params, reference_concentration,
!                  sediment_bedload, sediment_settling_velocity,
!                  thetastar_engelund_hansen
!
! Module calls -  Carr_at_U, Carr_at_V, error_alloc, exchange_mod, gauss_quad,
!                 mult_index, vector_mag_arr_atu, vector_mag_arr_atv
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
USE morphpars
USE physpars
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Carr_at_U, Carr_at_V
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: gauss_quad, vector_mag_arr_atu, vector_mag_arr_atv
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: f, i, ii, iiloc, iu, j, jj, jjloc, jv, morphsteps, n
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
REAL :: onethird = 1.0/3.0, real_small = 1.0E-04, rouse_max = 5.0, wmin = 0.03
REAL(KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:) :: ndgau_tr, nodes_gauss, &
                                                      & weights_gauss
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2, array2d3, aw,&
         & ca, canorm, curmag, cw, deltaw, fc, fcw, fgr, fw, f_c, f_w, ggr, &
         & me, mw, nw, oscwavetr, phi, phicur, phiwav, qprojatu, qprojatv, &
         & qtotatu_per, qtotatv_per, theta, thetastar, ua, uc, uof, uon, vc, &
         & veff, wavevelatu, wavevelatv, wheight, wnum, wperiod, ws, zclim, &
         & z_c, z_w
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: dstaratu, dstaratv, deltasatu, &
         & deltasatv, excesshearatu, excesshearatv, qnormatu, qnormatv, &
         & taucritatu, taucritatv, taunormatu, taunormatv, velwav, velcurwav2


procname(pglev+1) = 'sediment_totalload'
CALL log_timer_in()

!
!1. Allocate
!-----------
!

IF (nt.EQ.0.AND.iopt_waves.EQ.1) THEN
   ALLOCATE(ndgau_tr(nrquad_wav),STAT=errstat)
   CALL error_alloc('ndgau_tr',1,(/nrquad_wav/),kndrlong)      
   ndgau_tr = 0.0
   ALLOCATE(nodes_gauss(nrquad_wav),STAT=errstat)
   CALL error_alloc('nodes_gauss',1,(/nrquad_wav/),kndrlong)
   nodes_gauss = 0.0
   ALLOCATE(weights_gauss(nrquad_wav),STAT=errstat)
   CALL error_alloc('weights_gauss',1,(/nrquad_wav/),kndrlong)
   weights_gauss = 0.0
   CALL gauss_quad(nrquad_wav,weights_gauss,nodes_gauss)
   ndgau_tr = pi*(nodes_gauss+1.0)
ENDIF

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (array2d1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2d2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2d3(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d3',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (aw(ncloc,nrloc),STAT=errstat)
CALL error_alloc('aw',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (curmag(ncloc,nrloc),STAT=errstat)
CALL error_alloc('curmag',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (deltasatu(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('deltasatu',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (deltasatv(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('deltasatv',3,(/ncloc,nrloc,nf/),kndrtype)
ALLOCATE (qprojatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qprojatu',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (qprojatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qprojatv',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (ws(ncloc,nrloc),STAT=errstat)
CALL error_alloc('ws',2,(/ncloc,nrloc/),kndrtype)

IF (iopt_sed_toteq.LE.2.OR.iopt_sed_toteq.EQ.5) THEN
   ALLOCATE (qnormatu(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('qnormatu',3,(/ncloc,nrloc,nf/),kndrtype)
   ALLOCATE (qnormatv(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('qnormatv',3,(/ncloc,nrloc,nf/),kndrtype)
ENDIF

IF (iopt_sed_toteq.LE.4) THEN
   ALLOCATE (taunormatu(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('taunormatu',3,(/ncloc,nrloc,nf/),kndrtype)
   ALLOCATE (taunormatv(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('taunormatv',3,(/ncloc,nrloc,nf/),kndrtype)
ENDIF

IF (iopt_sed_toteq.LE.2.OR.iopt_sed_toteq.EQ.4) THEN
   ALLOCATE (theta(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('theta',2,(/ncloc,nrloc/),kndrtype)
ENDIF

IF (iopt_sed_toteq.EQ.3.OR.iopt_sed_toteq.GT.5) THEN
   ALLOCATE (cw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('cw',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (mw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('mw',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (nw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('nw',2,(/ncloc,nrloc/),kndrtype)
ENDIF

IF (iopt_sed_toteq.EQ.4.OR.iopt_sed_toteq.GT.5) THEN
   ALLOCATE (phiwav(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('phiwav',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (wavevelatu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('wavevelatu',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (wavevelatv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('wavevelatv',2,(/ncloc,nrloc/),kndrtype)
ENDIF

IF (iopt_sed_toteq.LE.2) THEN
   ALLOCATE (thetastar(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('thetastar',2,(/ncloc,nrloc/),kndrtype)

ELSEIF (iopt_sed_toteq.EQ.3) THEN
   ALLOCATE (dstaratu(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('dstaratu',3,(/ncloc,nrloc,nf/),kndrtype)
   ALLOCATE (dstaratv(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('dstaratv',3,(/ncloc,nrloc,nf/),kndrtype)
   ALLOCATE (fgr(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fgr',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (ggr(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('ggr',2,(/ncloc,nrloc/),kndrtype)
   
ELSEIF (iopt_sed_toteq.EQ.4) THEN
   ALLOCATE (deltaw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('deltaw',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (fc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fc',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (fcw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fcw',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (fw(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fw',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (oscwavetr(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('oscwavetr',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (phi(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('phi',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (phicur(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('phicur',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (qtotatu_per(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('qtotatu_per',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (qtotatv_per(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('qtotatv_per',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (uc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('uc',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (velcurwav2(ncloc,nrloc,nrquad_wav),STAT=errstat)
   CALL error_alloc('velcurwav2',3,(/ncloc,nrloc,nrquad_wav/),kndrtype)
   ALLOCATE (velwav(ncloc,nrloc,nrquad_wav),STAT=errstat)
   CALL error_alloc('velwav',3,(/ncloc,nrloc,nrquad_wav/),kndrtype)

ELSEIF (iopt_sed_toteq.EQ.5) THEN
   ALLOCATE (excesshearatu(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('excesshearatu',3,(/ncloc,nrloc,nf/),kndrtype)
   ALLOCATE (excesshearatv(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('excesshearatv',3,(/ncloc,nrloc,nf/),kndrtype)
   ALLOCATE (taucritatu(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('taucritatu',3,(/ncloc,nrloc,nf/),kndrtype)
   ALLOCATE (taucritatv(ncloc,nrloc,nf),STAT=errstat)
   CALL error_alloc('taucritatv',3,(/ncloc,nrloc,nf/),kndrtype)

ELSEIF (iopt_sed_toteq.GT.5) THEN
   ALLOCATE (ca(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('ca',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (canorm(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('canorm',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (f_c(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('f_c',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (f_w(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('f_w',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (me(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('me',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (ua(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('ua',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (uof(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('uof',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (uon(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('uon',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (vc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('vc',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (veff(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('veff',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (wheight(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('wheight',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (wnum(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('wnum',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (wperiod(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('wperiod',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (zclim(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('zclim',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (z_c(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('z_c',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (z_w(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('z_w',2,(/ncloc,nrloc/),kndrtype)

ENDIF

!
!2. Initialise arrays
!--------------------
!
!---masks
maskatu = node2du(1:ncloc,1:nrloc).GT.1
maskatv = node2dv(1:ncloc,1:nrloc).GT.1

!---current projection factor
WHERE (bstresatu_sed.GT.0.0)
   qprojatu = ubstresatu_sed(1:ncloc,:)/bstresatu_sed
ELSEWHERE
   qprojatu = 0.0
END WHERE
WHERE (bstresatv_sed.GT.0.0)
   qprojatv = vbstresatv_sed(:,1:nrloc)/bstresatv_sed
ELSEWHERE
   qprojatv = 0.0
END WHERE

!---critical and excess shear stress
IF (iopt_sed_toteq.EQ.5) THEN
   f_210: DO f=1,nf
      CALL Carr_at_U(bstres_cr(0:ncloc,1:nrloc,f),taucritatu(:,:,f),1,3,&
                   & (/0,1,1/),(/ncloc,nrloc,1/),1,iarr_bstres_cr,.TRUE.,&
                   & mask=bed_fraction(:,1:nrloc,1,f).GT.0.0)
      CALL Carr_at_V(bstres_cr(1:ncloc,0:nrloc,f),taucritatv(:,:,f),1,3,&
                  & (/1,0,1/),(/ncloc,nrloc,1/),1,iarr_bstres_cr,.TRUE.,&
                  & mask=bed_fraction(1:ncloc,:,1,f).GT.0.0)
      WHERE (maskatu)
         excesshearatu(:,:,f) = MAX(bstresatu_sed-taucritatu(:,:,f),0.0)
      END WHERE
      WHERE (maskatv)
         excesshearatv(:,:,f) = MAX(bstresatv_sed-taucritatv(:,:,f),0.0)
      END WHERE
   ENDDO f_210
ENDIF

!---normalising factors at U-nodes
CALL Carr_at_U(densw(0:ncloc,1:nrloc,1),array2d1,1,3,(/0,1,1/),&
            & (/ncloc,nrloc,1/),1,iarr_densw,.TRUE.)
CALL Carr_at_U(kinvisc(0:ncloc,1:nrloc,1),array2d2,1,3,(/0,1,1/),&
            & (/ncloc,nrloc,1/),1,0,.FALSE.)
f_220: DO f=1,nf
   WHERE (maskatu)
      deltasatu(:,:,f) = MAX(rhos(f)/array2d1-1.0,0.0)
   END WHERE
   IF (iopt_sed_toteq.LE.4) THEN
      WHERE (maskatu)
         taunormatu(:,:,f) = deltasatu(:,:,f)*gaccatu(1:ncloc,:)*dp(f)
      END WHERE
   ENDIF
   IF (iopt_sed_toteq.LE.2.OR.iopt_sed_toteq.EQ.5) THEN
      WHERE (maskatu)
         qnormatu(:,:,f) = SQRT(deltasatu(:,:,f)*gaccatu(1:ncloc,:)*dp(f)**3)
      END WHERE
   ENDIF
   IF (iopt_sed_toteq.EQ.3) THEN
      WHERE (maskatu)
         dstaratu(:,:,f) = dp(f)*(deltasatu(:,:,f)*gaccatu(1:ncloc,:)&
                               & /array2d2**2)**onethird
      END WHERE
   ENDIF
ENDDO f_220

!---normalising factors at V-nodes
CALL Carr_at_V(densw(1:ncloc,0:nrloc,1),array2d1,1,3,(/1,0,1/),&
            & (/ncloc,nrloc,1/),1,iarr_densw,.TRUE.)
CALL Carr_at_V(kinvisc(1:ncloc,0:nrloc,1),array2d2,1,3,(/1,0,1/),&
            & (/ncloc,nrloc,1/),1,0,.FALSE.)
f_230: DO f=1,nf
   WHERE (maskatv)
      deltasatv(:,:,f) = MAX(rhos(f)/array2d1-1.0,0.0)
   END WHERE
   IF (iopt_sed_toteq.LE.4) THEN
      WHERE (maskatv)
         taunormatv(:,:,f) = deltasatv(:,:,f)*gaccatv(:,1:nrloc)*dp(f)
      END WHERE
   ENDIF
   IF (iopt_sed_toteq.LE.2.OR.iopt_sed_toteq.EQ.5) THEN
      WHERE (maskatv)
         qnormatv(:,:,f) = SQRT(deltasatv(:,:,f)*gaccatv(:,1:nrloc)*dp(f)**3)
      END WHERE
   ENDIF
   IF (iopt_sed_toteq.EQ.3) THEN
      WHERE (maskatv)
         dstaratv(:,:,f) = dp(f)*(deltasatv(:,:,f)*gaccatv(:,1:nrloc)&
                               & /array2d2**2)**onethird
      END WHERE
   ENDIF
ENDDO f_230

!
!3. Total load transport formulae
!--------------------------------
!
!3.1 Engelund-Hansen
!-------------------
!

IF (iopt_sed_toteq.LE.2) THEN

!  ---U-nodes
   f_311: DO f = 1,nf
      WHERE (maskatu)
         theta = bstresatu_sed/taunormatu(:,:,f)
      ELSEWHERE
         theta = 0.0
      END WHERE
      CALL thetastar_engelund_hansen(theta,maskatu,thetastar,iopt_sed_toteq)
      WHERE (maskatu)
         qtotatu(1:ncloc,:,f) = 0.05*qprojatu*qnormatu(:,:,f)*thetastar**2.5&
                              & /bdragcoefatu_sed
      END WHERE
   ENDDO f_311

!  ---V-nodes
   f_312: DO f = 1,nf
      WHERE (maskatv)
         theta = bstresatv_sed/taunormatv(:,:,f)
      ELSEWHERE
         theta = 0.0
      END WHERE
      CALL thetastar_engelund_hansen(theta,maskatv,thetastar,iopt_sed_toteq)
      WHERE (maskatv)
         qtotatv(:,1:nrloc,f) = 0.05*qprojatv*qnormatv(:,:,f)*thetastar**2.5&
                              & /bdragcoefatv_sed
      ELSEWHERE
         qtotatv(:,1:nrloc,f) = 0.0
      END WHERE
   ENDDO f_312

!
!3.2 Ackers-White
!----------------
!

ELSEIF (iopt_sed_toteq.EQ.3) THEN

!  ---U-nodes
   CALL vector_mag_arr_atu(umvel(1:ncloc,1:nrloc),vmvel(0:ncloc,1:nrloc+1),&
                         & 1,3,1,1,iarr_hmvelmag,.TRUE.,vecmag=curmag)
   f_321: DO f=1,nf
      CALL ackerswhite_params(dstaratu(:,:,f),maskatu,nw,mw,aw,cw)
      WHERE (maskatu)
         fgr = SQRT(bstresatu_sed**nw/taunormatu(:,:,f))*&
             & (curmag/(2.46*LOG(10.0*deptotatu(1:ncloc,1:nrloc)/dp(f))))**&
             & (1.0-nw)
      ELSEWHERE
         fgr = 0.0
      END WHERE
      WHERE (fgr.GT.aw)
         ggr = cw*(fgr/aw-1.0)**mw
         qtotatu(1:ncloc,:,f) =  qprojatu*dp(f)*curmag*ggr&
                               & /SQRT(bdragcoefatu_sed**nw)
      ELSEWHERE
         qtotatu(1:ncloc,:,f) = 0.0
      END WHERE
   ENDDO f_321

!  ---V-nodes
   CALL vector_mag_arr_atv(umvel(1:ncloc+1,0:nrloc),vmvel(1:ncloc,1:nrloc),&
                         & 1,3,1,1,iarr_hmvelmag,.TRUE.,vecmag=curmag)
   f_320: DO f=1,nf
      CALL ackerswhite_params(dstaratv(:,:,f),maskatv,nw,mw,aw,cw)
      WHERE (maskatv)
         fgr = SQRT(bstresatv_sed**nw/taunormatv(:,:,f))*&
             & (curmag/(2.46*LOG(10.0*deptotatv(1:ncloc,1:nrloc)/dp(f))))**&
             & (1.0-nw)
      ELSEWHERE
         fgr = 0.0
      END WHERE
      WHERE (fgr.GT.aw)
         ggr = cw*(fgr/aw-1.0)**mw
         qtotatv(:,1:nrloc,f) = qprojatv*dp(f)*curmag*ggr&
                              & /SQRT(bdragcoefatv_sed**nw)
      ELSEWHERE
         qtotatv(:,1:nrloc,f) = 0.0
      END WHERE
   ENDDO f_320

!
!3.3 Madsen and Grant (1976)
!---------------------------
!

ELSEIF (iopt_sed_toteq.EQ.4) THEN

!  ---settling velocity
   CALL sediment_settling_velocity

!
!3.3.1 No wave case
!------------------
!


   IF (iopt_waves.EQ.0) THEN
      
!     ---U-nodes
      f_3311: DO f=1,nf
         CALL Carr_at_U(wfall(0:ncloc,1:nrloc,nz,f),ws,1,3,(/0,1,nz/),&
                     & (/ncloc,nrloc,nz/),1,iarr_wfall,.TRUE.)
         WHERE (taunormatu(:,:,f).GT.0.0)
            qtotatu(1:ncloc,:,f) = (40.0*dp(f))*qprojatu*ws*&
                                 & (bstresatu_sed/taunormatu(:,:,f))**3
         ELSEWHERE
            qtotatu(1:ncloc,:,f) = 0.0
         END WHERE
      ENDDO f_3311

!     ---V-nodes
      f_3312: DO f=1,nf
         CALL Carr_at_V(wfall(1:ncloc,0:nrloc,nz,f),ws,1,3,(/1,0,nz/),&
                     & (/ncloc,nrloc,nz/),1,iarr_wfall,.TRUE.)
         WHERE (taunormatv(:,:,f).GT.0.0)
            qtotatv(:,1:nrloc,f) = (40.0*dp(f))*qprojatv*ws*&
                                 & (bstresatv_sed/taunormatv(:,:,f))**3
         ELSEWHERE
            qtotatv(:,1:nrloc,f) = 0.0
         END WHERE
      ENDDO f_3312
      
   ELSE

!
!3.3.2 U-nodes
!-------------
!
!     ---wave parameters and angle between waves and currents
      CALL Carr_at_U(wavevel(:,1:nrloc),wavevelatu,1,3,(/0,1,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_wavevel,.TRUE.)
      CALL Carr_at_U(waveexcurs(:,1:nrloc),aw,1,3,(/0,1,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_waveexcurs,.TRUE.)
      CALL Carr_at_U(wavedir(:,1:nrloc),phiwav,1,3,(/0,1,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_wavedir,.TRUE.,angle=.TRUE.)
      CALL vector_mag_arr_atu(umvel(1:ncloc,1:nrloc),&
                            & vmvel(0:ncloc,1:nrloc+1),1,3,1,1,iarr_hmvelmag,&
                           & .TRUE.,vecpha=phicur)
      WHERE (maskatu)
         phi = phiwav - phicur
      END WHERE

!     ---boundary thickness
      WHERE (aw.GT.0.0)
         deltaw = 0.072*aw**0.75*(30.0*zroughatu_sed)**0.25
         deltaw = MAX(30.0*deltaw,30.0*zroughatu_sed)
      ELSEWHERE
         deltaw = 30.0*zroughatu_sed
      END WHERE

!     ---current at top of boundary layer 
      WHERE (deltaw.GT.0.0)
         uc = SQRT(bstresatu_sed)*LOG(deltaw/zroughatu_sed)/ckar
         array2d1 = uc + wavevelatu
      ELSEWHERE
         uc = 0.0
         array2d1 = 0.0
      END WHERE

!     ---friction factors
      WHERE (array2d1.GT.0.0)
         array2d2 = uc/array2d1
      ELSEWHERE
         array2d2 = 0.0
      END WHERE
      WHERE (maskatu)
         fc = 0.32/(LOG(deptotatu(1:ncloc,1:nrloc)/zroughatu_sed)-1.0)**2
      ELSEWHERE
         fc = 0.0
      END WHERE
      WHERE (aw.GT.0.0)
         fw = 1.39/(aw/zroughatu_sed)**0.52
         fcw = array2d2*fc + (1.0-array2d2)*fw
      ELSEWHERE
         fcw = fc
      END WHERE
      
!     ---total current during wave period
      n_3321: DO n=1,nrquad_wav
         WHERE (maskatu)
            velwav(:,:,n) = wavevelatu*SIN(ndgau_tr(n))
            velcurwav2(:,:,n) = velwav(:,:,n)*(velwav(:,:,n)+2.0*uc*COS(phi))+&
                              & uc**2
         ELSEWHERE
            velcurwav2(:,:,n) = 0.0
         END WHERE
      ENDDO n_3321
  
!     ---total load per fraction
      f_3322: DO f=1,nf
         CALL Carr_at_U(wfall(0:ncloc,1:nrloc,nz,f),ws,1,3,(/0,1,nz/),&
                     & (/ncloc,nrloc,nz/),1,iarr_wfall,.TRUE.)
         qtotatu(1:ncloc,:,f) = 0.0; qtotatu_per = 0.0
         n_33221: DO n=1,nrquad_wav
            WHERE (velcurwav2(:,:,n).GT.0.0)
               theta = 0.5*fcw*velcurwav2(:,:,n)/taunormatu(:,:,f)
               array2d1 = 40.0*ws*dp(f)*theta**3/SQRT(velcurwav2(:,:,n))
               oscwavetr = array2d1*(uc+wavevelatu*COS(phi)*SIN(ndgau_tr(n)))
               qtotatu(1:ncloc,:,f) = qtotatu(1:ncloc,:,f) + &
                                    & 0.5*weights_gauss(n)*oscwavetr
               oscwavetr = array2d1*wavevelatu*SIN(phi)*SIN(ndgau_tr(n))
               qtotatu_per = qtotatu_per + 0.5*weights_gauss(n)*oscwavetr
            END WHERE
         ENDDO n_33221
         WHERE (maskatu)
            qtotatu(1:ncloc,:,f) = COS(phicur)*qtotatu(1:ncloc,:,f) - &
                                 & SIN(phicur)*qtotatu_per
         END WHERE
      ENDDO f_3322

!
!3.3.3 V-nodes
!-------------
!
!     ---wave parameters and angle between waves and currents
      CALL Carr_at_V(wavevel(1:ncloc,:),wavevelatv,1,3,(/1,0,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_wavevel,.TRUE.)
      CALL Carr_at_V(waveexcurs(1:ncloc,:),aw,1,3,(/1,0,1/),&
                  & (/ncloc,nrloc,1/),1,iarr_waveexcurs,.TRUE.)
      CALL Carr_at_V(wavedir(1:ncloc,:),phiwav,1,3,(/1,0,1/),&
                  & (/ncloc,nrloc,1,1/),1,iarr_wavedir,.TRUE.,angle=.TRUE.)
      CALL vector_mag_arr_atv(umvel(1:ncloc+1,0:nrloc),&
                            & vmvel(1:ncloc,1:nrloc),1,3,1,1,iarr_hmvelmag,&
                           & .TRUE.,vecpha=phicur)
      WHERE (maskatv)
         phi = phiwav - phicur
      END WHERE

!     ---boundary thickness
      WHERE (aw.GT.0.0)
         deltaw = 0.072*aw**0.75*(30.0*zroughatv_sed)**0.25
         deltaw = MAX(30.0*deltaw,30.0*zroughatv_sed)
      ELSEWHERE
         deltaw = 30.0*zroughatv_sed
      END WHERE

!     ---current at top of boundary layer 
      WHERE (deltaw.GT.0.0)
         uc = SQRT(bstresatv_sed)*LOG(deltaw/zroughatv_sed)/ckar
         array2d1 = uc + wavevelatv
      ELSEWHERE
         uc = 0.0
         array2d1 = 0.0
      END WHERE

!     ---friction factors
      WHERE (array2d1.GT.0.0)
         array2d2 = uc/array2d1
      ELSEWHERE
         array2d2 = 0.0
      END WHERE
      WHERE (maskatv)
         fc = 0.32/((LOG(deptotatv(1:ncloc,1:nrloc)/zroughatv_sed)-1.0)**2)
      ELSEWHERE
         fc = 0.0
      END WHERE
      WHERE (aw.GT.0.0)
         fw = 1.39/(aw/zroughatv_sed)**0.52
         fcw = array2d2*fc + (1.0-array2d2)*fw
      ELSEWHERE
         fcw =  fc
      END WHERE

!     ---total current during wave period
      n_3331: DO n=1,nrquad_wav
         WHERE (maskatv)
            velwav(:,:,n) = wavevelatv*SIN(ndgau_tr(n))
            velcurwav2(:,:,n) = velwav(:,:,n)*(velwav(:,:,n)+2.0*uc*COS(phi))+ &
                              & uc**2
         ELSEWHERE
            velcurwav2(:,:,n) = 0.0
         END WHERE
      ENDDO n_3331
      
!     ---total load per fraction
      f_3332: DO f=1,nf
         CALL Carr_at_V(wfall(1:ncloc,0:nrloc,nz,f),ws,1,3,(/1,0,nz/),&
                     & (/ncloc,nrloc,nz/),1,iarr_wfall,.TRUE.)
         qtotatv(:,1:nrloc,f) = 0.0; qtotatv_per = 0.0
         n_33321: DO n=1,nrquad_wav
            WHERE (velcurwav2(:,:,n).GT.0.0)
               theta = 0.5*fcw*velcurwav2(:,:,n)/taunormatv(:,:,f)
               array2d1 = 40.0*ws*dp(f)*theta**3/SQRT(velcurwav2(:,:,n))
               oscwavetr = array2d1*(uc+wavevelatv*COS(phi)*SIN(ndgau_tr(n)))
               qtotatv(:,1:nrloc,f) = qtotatv(:,1:nrloc,f) + &
                                    & 0.5*weights_gauss(n)*oscwavetr
               oscwavetr = array2d1*wavevelatv*SIN(phi)*SIN(ndgau_tr(n))
               qtotatv_per = qtotatv_per + 0.5*weights_gauss(n)*oscwavetr
            END WHERE
         ENDDO n_33321
         WHERE (maskatv)
            qtotatv(:,1:nrloc,f) = SIN(phicur)*qtotatv(:,1:nrloc,f) + &
                                 & COS(phicur)*qtotatv_per
         END WHERE
      ENDDO f_3332

   ENDIF

!
!3.4 Wu et al. (2000)
!--------------------
!

ELSEIF (iopt_sed_toteq.EQ.5) THEN

!  ---settling velocity
   CALL sediment_settling_velocity

!  ---bed load
   CALL sediment_bedload

!
!3.4.1 U-nodes
!-------------
!
!  ---current magnitude
   CALL vector_mag_arr_atu(umvel(1:ncloc,1:nrloc),vmvel(0:ncloc,1:nrloc+1),&
                         & 1,3,1,1,iarr_hmvelmag,.TRUE.,vecmag=curmag)

!   ---total load
    f_341: DO f=1,nf
       CALL Carr_at_U(wfall(0:ncloc,1:nrloc,nz,f),ws,1,3,(/0,1,nz/),&
                   & (/ncloc,nrloc,nz/),1,iarr_wfall,.TRUE.)
       WHERE (ws.GT.0.0.AND.excesshearatu(:,:,f).GT.0.0)
          array2d1 = 2.62E-05*((excesshearatu(:,:,f)/taucritatu(:,:,f))*&
                   & (curmag/ws))**1.74
       ELSEWHERE
          array2d1 = 0.0
       END WHERE
       WHERE (maskatu)
          qtotatu(1:ncloc,:,f) = qnormatu(:,:,f)*array2d1*qprojatu + &
                               & qbedatu(1:ncloc,:,f)
       END WHERE
    ENDDO f_341

!
!3.4.2 V-nodes
!-------------
!
!  ---current magnitude
   CALL vector_mag_arr_atv(umvel(1:ncloc+1,0:nrloc),vmvel(1:ncloc,1:nrloc),&
                         & 1,3,1,1,iarr_hmvelmag,.TRUE.,vecmag=curmag)

!   ---total load
    f_342: DO f=1,nf
       CALL Carr_at_V(wfall(1:ncloc,0:nrloc,nz,f),ws,1,3,(/1,0,nz/),&
                  & (/ncloc,nrloc,nz/),1,iarr_wfall,.TRUE.)
       WHERE (ws.GT.0.0.AND.taucritatv(:,:,f).GT.0.0)
          array2d1 = 2.62E-05*((excesshearatv(:,:,f)/taucritatv(:,:,f))*&
                   & (curmag/ws))**1.74
       ELSEWHERE
          array2d1 = 0.0
       END WHERE
       WHERE (maskatv)
          qtotatv(:,1:nrloc,f) = qnormatv(:,:,f)*array2d1*qprojatv + &
                               & qbedatv(:,1:nrloc,f)
       END WHERE
    ENDDO f_342

!
!3.5 Van Rijn (2003/2007)
!-------------------------
!

ELSEIF (iopt_sed_toteq.GT.5) THEN

!  ---settling velocity
   CALL sediment_settling_velocity

!  ---bed load
   CALL sediment_bedload

!  ---reference concentration
   CALL reference_concentration

!
!3.5.1 U-nodes
!-------------
!
!3.5.1.1 Initialise parameters
!-----------------------------
!
!   ---current magnitude
    CALL vector_mag_arr_atu(umvel(1:ncloc,1:nrloc),vmvel(0:ncloc,1:nrloc+1),&
                          & 1,3,1,1,iarr_hmvelmag,.TRUE.,vecmag=curmag)

!   ----wave related parameters
    IF (iopt_waves.EQ.1) THEN
!      --wave height and period
       CALL Carr_at_U(waveheight(:,1:nrloc),wheight,1,3,(/0,1,1/),&
                   & (/ncloc,nrloc,1/),1,iarr_waveheight,.TRUE.)
       CALL Carr_at_U(waveperiod(:,1:nrloc),wperiod,1,3,(/0,1,1/),&
                   & (/ncloc,nrloc,1/),1,iarr_waveperiod,.TRUE.)
!      --wave direction
       CALL Carr_at_U(wavedir(:,1:nrloc),phiwav,1,3,(/0,1,1/),&
                   & (/ncloc,nrloc,1/),1,iarr_wavedir,.TRUE.,angle=.TRUE.)
!      --wave velocity asymetry in on-shore and off-shore direction
!        (Grasmeijer en van Rijn 1998)
       CALL Carr_at_U(wavevel(:,1:nrloc),wavevelatu,1,3,(/0,1,1/),&
                   & (/ncloc,nrloc,1/),1,iarr_wavevel,.TRUE.)
       CALL Carr_at_U(wavenum(:,1:nrloc),wnum,1,3,(/0,1,1/),&
                   & (/ncloc,nrloc,1/),1,iarr_wavenum,.TRUE.)
       WHERE (maskatu)
           array2d1 = wnum*wheight/twopi
           ! u hat
           array2d1 = (2.0-6.4*array2d1**0.65*(array2d1)**&
                    & (0.541*wnum*deptotatu(1:ncloc,1:nrloc)))*wavevelatu
           !procedure to find onshore orbital velocity peak
           aw = wperiod*SQRT(gaccatu(1:ncloc,:)/deptotatu(1:ncloc,1:nrloc))
           !L5
           mw = MERGE(3.2E-3*aw**2+8e-5*aw**3,5.6E-3*aw**2+4E-5*aw**3,&
                    & aw.LE.30.0)
           !L4
           aw = MERGE(-15.0+1.35*aw,-2.7+0.53*aw,aw.LE.15.0)
           !L3
           nw = (0.5-mw)/(aw-1.0+EXP(-aw))
           !L2
           mw = nw*aw + mw
           cw = array2d1/SQRT(gaccatu(1:ncloc,:)*deptotatu(1:ncloc,1:nrloc))
           !U_on
           uon = array2d1*(0.5-nw+mw*cw+nw*EXP(-aw*cw))
           !U_of
           uof = array2d1 - uon
           !v_eff
           veff = SQRT(curmag**2+uon**2)
           !velocity asymmetry
           ua = (uon**4-uof**4)/(uon**3+uof**3)
      END WHERE

   ENDIF

!
!3.5.1.2 Total load
!------------------
!

  f_3512: DO f=1,nf

!    ---interpolate at U-nodes
     CALL Carr_at_U(wfall(0:ncloc,1:nrloc,nz,f),ws,1,3,(/0,1,nz/),&
                 & (/ncloc,nrloc,nz/),1,iarr_wfall,.TRUE.)
     CALL Carr_at_U(cref(0:ncloc,1:nrloc,f),ca,1,3,(/0,1,1/),&
                 & (/ncloc,nrloc,1/),1,iarr_cref,.TRUE.)
!    ---transport due to currents
!    --current part

     WHERE (bstresatu_sed.GT.1.0E-30)
       array2d1 = 1.0 + 2.0*(ws**2)/bstresatu_sed
       array2d1 = MAX(beta_sed_min,MIN(array2d1,beta_sed_max))
       array2d2 = ws/SQRT(bstresatu_sed)
     ELSEWHERE
        array2d1 = beta_sed_max
        array2d2 = rouse_max*ckar
     END WHERE

     IF (iopt_sed_toteq.EQ.6) THEN
        WHERE (maskatu)
           z_c = array2d2/(array2d1*ckar) +  2.5*array2d2**0.8*(ca/cmax)**0.4
        END WHERE
     ELSE
        WHERE (maskatu)
           array2d3 = (ca/cmax)**0.4
           z_c = array2d2/(array2d1*ckar*(1.0+array2d3**2-2.0*array2d3))
        END WHERE
     ENDIF

     WHERE (maskatu)
        canorm = MAX(30*zroughatu_sed,0.02)/deptotatu(1:ncloc,1:nrloc)
        canorm = MIN(canorm,0.99)
!       --prevent numerical problems for high Rouse numbers
!       --location of maximum of F(Z') as function of a/h (fitted numerically)
        zclim = 1.1821*canorm**(-0.96038)
        z_c = MIN(z_c,zclim)
     ELSEWHERE
        z_c = 1.2
     END WHERE
     WHERE (ABS(z_c-1.2).GT.real_small)
        f_c = (canorm**z_c-canorm**1.2)/((1-canorm)**z_c*(1.2-z_c))
        f_c = MAX(f_c,0.0)
     ELSEWHERE
        f_c = 0.1649
     END WHERE

!    --wave part
     IF (iopt_waves.EQ.1) THEN
        WHERE (wheight.GT.0.0)
           z_w = 4.0*(0.2*deptotatu(1:ncloc,1:nrloc))**0.6*(MAX(ws,wmin)*&
               & (wperiod/wheight))**0.8
           z_w = MIN(z_w,zclim)
        ELSEWHERE
           z_w = 1.2
        END WHERE
        WHERE (ABS(z_w-1.2).GT.real_small)
           f_w = (canorm**z_w-canorm**1.2)/((1-canorm)**z_w*(1.2-z_w))
           f_w = MAX(f_w,0.0)
        ELSEWHERE
           f_w = 0.1649
        END WHERE
     ELSE
        f_w = 0.0
     ENDIF

!    --currents+waves
     WHERE (maskatu)
        qtotatu(1:ncloc,:,f) = (f_c+f_w)*deptotatu(1:ncloc,1:nrloc)*curmag*ca
     END WHERE

!    ---transport due to waves
     IF (iopt_waves.EQ.1) THEN
        WHERE (maskatu)
           array2d1 = LOG(2.0*deptotatu(1:ncloc,1:nrloc)/dp(f))
           vc = MERGE(0.19*dp(f)**0.1*array2d1,8.5*dp(f)**0.6*array2d1,&
                    & dp(f).LT.0.0005)   
           me = (veff-vc)**2/(deltasatu(:,:,f)*gaccatu(1:ncloc,:))
           qtotatu(1:ncloc,:,f) = qtotatu(1:ncloc,:,f) + &
                                & 0.0014*me*ua*COS(phiwav)
         END WHERE
      ENDIF

!    ---project from current to X-direction
      WHERE (maskatu)
         qtotatu(1:ncloc,:,f) = qtotatu(1:ncloc,:,f)*qprojatu + &
                              & qbedatu(1:ncloc,:,f)
      END WHERE

   ENDDO f_3512

!
!3.5.2 V-nodes
!-------------
!
!3.5.2.1 Initialise parameters
!-----------------------------
!
!   ---current magnitude
    CALL vector_mag_arr_atv(umvel(1:ncloc+1,0:nrloc),vmvel(1:ncloc,1:nrloc),&
                          & 1,3,1,1,iarr_hmvelmag,.TRUE.,vecmag=curmag)

!   ----wave related parameters
    IF (iopt_waves.EQ.1) THEN
!      --wave height and period
       CALL Carr_at_V(waveheight(1:ncloc,:),wheight,1,3,(/1,0,1/),&
                   & (/ncloc,nrloc,1/),1,iarr_waveheight,.TRUE.)
       CALL Carr_at_V(waveperiod(1:ncloc,:),wperiod,1,3,(/1,0,1/),&
                   & (/ncloc,nrloc,1/),1,iarr_waveperiod,.TRUE.)
!      --wave direction
       CALL Carr_at_V(wavedir(1:ncloc,:),array2d1,1,3,(/1,0,1/),&
                   & (/ncloc,nrloc,1/),1,iarr_wavedir,.TRUE.,angle=.TRUE.)

!        (Grasmeijer en van Rijn 1998)
       CALL Carr_at_V(wavevel(1:ncloc,:),wavevelatv,1,3,(/1,0,1/),&
                   & (/ncloc,nrloc,1/),1,iarr_wavevel,.TRUE.)
       CALL Carr_at_V(wavenum(1:ncloc,:),wnum,1,3,(/1,0,1/),&
                   & (/ncloc,nrloc,1/),1,iarr_wavenum,.TRUE.)
       WHERE (maskatv)
           array2d1 = wnum*wheight/twopi
           ! u hat
           array2d1 = (2.0-6.4*array2d1**0.65*(array2d1)**&
                    & (0.541*wnum*deptotatv(1:ncloc,1:nrloc)))*wavevelatv
           !procedure to find onshore orbital velocity peak
           aw = wperiod*SQRT(gaccatv(:,1:nrloc)/deptotatv(1:ncloc,1:nrloc))
           !L5
           mw = MERGE(3.2E-3*aw**2+8.0E-5*aw**3,5.6E-3*aw**2+4E-5*aw**3,&
                    & aw.LE.30.0)
           !L4
           aw = MERGE(-15.0+1.35*aw,-2.7+0.53*aw,aw.LE.15.0)
           !L3
           nw = (0.5-mw)/(aw-1.0+EXP(-aw))
           !L2
           mw = nw*aw + mw
           cw = array2d1/SQRT(gaccatv(:,1:nrloc)*deptotatv(1:ncloc,1:nrloc))
           !U_on
           uon = array2d1*(0.5-nw+mw*cw+nw*EXP(-aw*cw))
           !U_of
           uof = array2d1 - uon
           !v_eff
           veff = SQRT(curmag**2+uon**2)
           !velocity asymmetry
           ua = (uon**4-uof**4)/(uon**3+uof**3)
      END WHERE

  ENDIF

!
!3.5.2.2 Total load
!------------------
!

  f_3522: DO f=1,nf

!    ---interpolate at V-nodes
     CALL Carr_at_V(wfall(1:ncloc,0:nrloc,nz,f),ws,1,3,(/1,0,nz/),&
                 & (/ncloc,nrloc,nz/),1,iarr_wfall,.TRUE.)
     CALL Carr_at_V(cref(1:ncloc,0:nrloc,f),ca,1,3,(/1,0,1/),&
                 & (/ncloc,nrloc,1/),1,iarr_cref,.TRUE.)

!    ---transport due to currents
!    --current part
     WHERE (bstresatv_sed.GT.1.0E-30)
       array2d1 = 1.0 + 2.0*(ws**2)/bstresatv_sed
       array2d1 = MAX(beta_sed_min,MIN(array2d1,beta_sed_max))
       array2d2 = ws/SQRT(bstresatv_sed)
     ELSEWHERE
        array2d1 = beta_sed_max
        array2d2 = rouse_max*ckar
     END WHERE
     IF (iopt_sed_toteq.EQ.6) THEN
        WHERE (maskatv)
           z_c = array2d2/(array2d1*ckar) +  2.5*array2d2**0.8*(ca/cmax)**0.4
        END WHERE
     ELSE
        WHERE (maskatv)
           array2d3 = (ca/cmax)**0.4
           z_c = array2d2/(array2d1*ckar*(1.0+array2d3**2-2.0*array2d3))
        END WHERE
     ENDIF
     WHERE (maskatv)
        canorm = MAX(30*zroughatv_sed,0.02)/deptotatv(1:ncloc,1:nrloc)
        canorm = MIN(canorm,0.99)
!       --prevent numerical problems for high Rouse numbers
!       --location of maximum of F(Z') as function of a/h (fitted numerically)
        zclim = 1.1821*canorm**(-0.96038)
        z_c = MIN(z_c,zclim)
     ELSEWHERE
        z_c = 1.2
     END WHERE
     WHERE (ABS(z_c-1.2).GT.real_small)
        f_c = (canorm**z_c-canorm**1.2)/((1-canorm)**z_c*(1.2-z_c))
        f_c = MAX(f_c,0.0)
     ELSEWHERE
        f_c = 0.1649
     END WHERE
!    --wave part
     IF (iopt_waves.EQ.1) THEN
        WHERE (wheight.GT.0.0)
           z_w = 4.0*(0.2*deptotatv(1:ncloc,1:nrloc))**0.6*(MAX(ws,wmin)*&
               & (wperiod/wheight))**0.8
           z_w = MIN(z_w,zclim)
        ELSEWHERE
           z_w = 1.2
        END WHERE
        WHERE (ABS(z_w-1.2).GT.real_small)
           f_w = (canorm**z_w-canorm**1.2)/((1-canorm)**z_w*(1.2-z_w))
           f_w = MAX(f_w,0.0)
        ELSEWHERE
           f_w = 0.1649
        END WHERE
     ELSE
        f_w = 0.0
     ENDIF
!    --currents+waves
     WHERE (maskatv)
        qtotatv(:,1:nrloc,f) = (f_c+f_w)*deptotatv(1:ncloc,1:nrloc)*curmag*ca
     END WHERE

!    ---transport due to waves
     IF (iopt_waves.EQ.1) THEN
        WHERE (maskatv)
           array2d1 = LOG(2.0*deptotatv(1:ncloc,1:nrloc)/dp(f))
           vc = MERGE(0.19*dp(f)**0.1*array2d1,8.5*dp(f)**0.6*array2d1,&
                    & dp(f).LT.0.0005)
           me = (veff-vc)**2/(deltasatv(:,:,f)*gaccatv(:,1:nrloc))
           qtotatv(:,1:nrloc,f) = qtotatv(:,1:nrloc,f) + &
                                & 0.0014*me*ua*SIN(phiwav)
        END WHERE
     ENDIF

!    ---project from current to X-direction
      WHERE (maskatv)
         qtotatv(:,1:nrloc,f) = qtotatv(:,1:nrloc,f)*qprojatv + &
                              & qbedatv(:,1:nrloc,f)
      END WHERE

   ENDDO f_3522

ENDIF

!
!4. Fractional transport
!-----------------------
!
!4.1 U-nodes
!-----------
!

WHERE (.NOT.maskatu) array2d1 = 0.0

f_410: DO f=1,nf
   WHERE (maskatu)
      array2d1 = SIGN(1.0,qtotatu(1:ncloc,:,f))
   END WHERE
   WHERE (array2d1.GT.0.0)
      qtotatu(1:ncloc,:,f) = bed_fraction(0:ncloc-1,1:nrloc,1,f)*&
                           & qtotatu(1:ncloc,:,f)
   ELSEWHERE
      qtotatu(1:ncloc,:,f) = bed_fraction(1:ncloc,1:nrloc,1,f)*&
                           & qtotatu(1:ncloc,:,f)
   END WHERE

ENDDO f_410

!
!4.2 V-nodes
!-----------
!

WHERE (.NOT.maskatv) array2d1 = 0.0

f_420: DO f=1,nf
   WHERE (maskatv)
      array2d1 = SIGN(1.0,qtotatv(:,1:nrloc,f))
   END WHERE
   WHERE (array2d1.GT.0.0)
      qtotatv(:,1:nrloc,f) = bed_fraction(1:ncloc,0:nrloc-1,1,f)*&
                           & qtotatv(:,1:nrloc,f)
   ELSEWHERE
      qtotatv(:,1:nrloc,f) = bed_fraction(1:ncloc,1:nrloc,1,f)*&
                           & qtotatv(:,1:nrloc,f)
   END WHERE

ENDDO f_420

!
!5. Exchange halos
!-----------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = (/0,1,1/); nhexch = (/1,1,0,0/)
   CALL exchange_mod(qtotatu,lbounds,nhexch,iarr_qtotatu)
   lbounds = (/1,0,1/); nhexch = (/0,0,1,1/)
   CALL exchange_mod(qtotatv,lbounds,nhexch,iarr_qtotatv)
ENDIF

!
!6. Upwind correction at open boundaries
!---------------------------------------
!

IF (iopt_sed_obc_flux.EQ.1) THEN

!
!6.1 U-nodes
!-----------
!

   iiloc_610: DO iiloc=1,nobuloc
      i = iobuloc(iiloc); j = jobuloc(iiloc)
      ii = indexobu(iiloc); iu = MERGE(i+1,i-1,westobu(ii))
      qtotatu(i,j,:) = (delyatu(iu,j)/delyatu(i,j))*qtotatu(iu,j,:)
   ENDDO iiloc_610
 
!
!6.2 V-nodes
!-------------
!

   jjloc_620: DO jjloc=1,nobvloc
      i = iobvloc(jjloc); j = jobvloc(jjloc)
      jj = indexobv(jjloc); jv = MERGE(j+1,j-1,soutobv(jj))
      qtotatv(i,j,:) = (delxatv(i,jv)/delxatv(i,j))*qtotatv(i,jv,:)
   ENDDO jjloc_620

ENDIF

!
!7. Bed load averaged over morphological step
!--------------------------------------------
!

IF (iopt_morph.EQ.1) THEN
   morphsteps = icmorph/icsed
   f_710: DO f=1,nf
      IF (icmorph.EQ.icsed) THEN
         WHERE (node2du(0:ncloc+1,1:nrloc).GT.1)
            qtotatu_int(:,:,f) = qtotatu(:,:,f)
         END WHERE
         WHERE (node2dv(1:ncloc,0:nrloc+1).GT.1)
            qtotatv_int(:,:,f) = qtotatv(:,:,f)
         END WHERE
      ELSE
         IF (mult_index(nt,icmorph)) THEN
            qtotatu_int(:,:,f) = 0.0; qtotatv_int(:,:,f) = 0.0
         ENDIF
         WHERE (node2du(0:ncloc+1,1:nrloc).GT.1)
            qtotatu_int(:,:,f) = qtotatu_int(:,:,f) + qtotatu(:,:,f)/morphsteps
         END WHERE
         WHERE (node2dv(1:ncloc,0:nrloc+1).GT.1)
            qtotatv_int(:,:,f) = qtotatv_int(:,:,f) + qtotatv(:,:,f)/morphsteps
         END WHERE
      ENDIF
   ENDDO f_710
ENDIF

!
!8. Deallocate
!-------------
!

IF ((cold_start.OR.nt.EQ.nstep).AND.iopt_waves.EQ.1) THEN
   DEALLOCATE (ndgau_tr,nodes_gauss,weights_gauss)
ENDIF

DEALLOCATE (maskatu,maskatv)
DEALLOCATE (array2d1,array2d2,array2d3,aw,curmag,deltasatu,deltasatv,&
             & qprojatu,qprojatv,ws)
IF (iopt_sed_toteq.LE.2.OR.iopt_sed_toteq.EQ.5) THEN
   DEALLOCATE (qnormatu,qnormatv)
ENDIF
IF (iopt_sed_toteq.LE.4) DEALLOCATE (taunormatu,taunormatv)
IF (iopt_sed_toteq.LE.2.OR.iopt_sed_toteq.EQ.4) DEALLOCATE (theta)
IF (iopt_sed_toteq.EQ.3.OR.iopt_sed_toteq.GT.5) DEALLOCATE (cw,mw,nw)
IF (iopt_sed_toteq.EQ.4.OR.iopt_sed_toteq.GT.5) THEN
   DEALLOCATE (phiwav,wavevelatu,wavevelatv)
ENDIF
IF (iopt_sed_toteq.LE.2) THEN
   DEALLOCATE (thetastar)
ELSEIF (iopt_sed_toteq.EQ.3) THEN
   DEALLOCATE (dstaratu,dstaratv,fgr,ggr)
ELSEIF (iopt_sed_toteq.EQ.4) THEN
   DEALLOCATE (deltaw,fc,fcw,fw,oscwavetr,phi,phicur,qtotatu_per,&
             & qtotatv_per,uc,velcurwav2,velwav)
ELSEIF (iopt_sed_toteq.EQ.5) THEN
   DEALLOCATE (excesshearatu,excesshearatv,taucritatu,taucritatv)
ELSEIF (iopt_sed_toteq.GT.5) THEN
   DEALLOCATE (ca,canorm,f_c,f_w,me,ua,uof,uon,vc,veff,wheight,wnum,wperiod,&
             & zclim,z_c,z_w)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE sediment_totalload

!========================================================================

SUBROUTINE sediment_settling_velocity
!************************************************************************
!
! *sediment_settling_velocity* settling velocity for sediments
!
! Author - Boudewijn Decrop (IMDC)
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - sediment_advdiff, sediment_totalload
!
! External calls - 
!
! Module calls - Carr_at_W, error_alloc, exchange_mod
!
!************************************************************************
!
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
USE timepars
USE turbulence
USE array_interp, ONLY: Carr_at_W
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
! Local variables
!
INTEGER :: f, k
INTEGER, DIMENSION(4) :: lbounds, nhexch
REAL :: acam, bcam, camfac1, camfac2, camfac3, mcam, xmin
REAL, PARAMETER :: half = 0.5, onethird = 1.0/3.0
REAL, ALLOCATABLE, DIMENSION(:,:) ::  array2dc, dstar, flocVL, flocVR
REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  array3dw, ctotatw
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) ::  deltas


procname(pglev+1) = 'sediment_settling_velocity'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array3dw(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('array3dw',3,(/ncloc,nrloc,nz+1/),kndrtype)
IF (iopt_sed_ws.GT.1) THEN
   ALLOCATE (deltas(ncloc,nrloc,nz+1,nf),STAT=errstat)
   CALL error_alloc('deltas',4,(/ncloc,nrloc,nz+1,nf/),kndrtype)
ENDIF
IF (iopt_sed_ws.EQ.2.OR.iopt_sed_ws.EQ.3.OR.iopt_sed_ws.EQ.5) THEN
   ALLOCATE (dstar(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('dstar',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (iopt_sed_ws_floc.GT.0) THEN
   ALLOCATE (flocVL(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('flocVL',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (flocVR(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('flocVR',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (iopt_sed_ws_floc.GT.0.OR.iopt_sed_ws_hindset.GT.0) THEN
   ALLOCATE (ctotatw(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('ctotatw',2,(/ncloc,nrloc/),kndrtype)
ENDIF

!
!2. Initialise arrays
!--------------------
!
!---density at W-nodes
IF (nz.GT.1) THEN
   CALL Carr_at_W(densw(1:ncloc,1:nrloc,:),array3dw(:,:,2:nz),&
                &(/1,1,1/),(/ncloc,nrloc,nz/),1,iarr_densw,.TRUE.)
ENDIF
WHERE (maskatc_int)
   array3dw(:,:,nz+1) = densw(1:ncloc,1:nrloc,nz)
   array3dw(:,:,1) = densw(1:ncloc,1:nrloc,1)
END WHERE

!---density factor
IF (iopt_sed_ws.GT.1) THEN
   f_210: DO f=1,nf
   k_210: DO k=1,nz+1
      WHERE (maskatc_int)
         deltas(:,:,k,f) = MAX(rhos(f)/array3dw(:,:,k)-1.0,0.0)
      END WHERE
   ENDDO k_210
   ENDDO f_210
ENDIF

!---concentration at W-nodes
IF (iopt_sed_ws_floc.GT.0.OR.iopt_sed_ws_hindset.GT.0) THEN
   IF (nz.GT.1) THEN
      CALL Carr_at_W(ctot(1:ncloc,1:nrloc,:),ctotatw(:,:,2:nz),(/1,1,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_ctot,.TRUE.)
   ENDIF
   WHERE (maskatc_int)
      ctotatw(:,:,1) = ctot(1:ncloc,1:nrloc,1)
      ctotatw(:,:,nz+1) = ctot(1:ncloc,1:nrloc,nz)
   END WHERE
ENDIF

!
!3. Single particle settling
!---------------------------
!
!3.1 Uniform velocity
!--------------------
!

IF (iopt_sed_ws.EQ.1) THEN
   f_310: DO f=1,nf
   k_310: DO k=1,nz+1
      WHERE (maskatc_int)
         wfall(1:ncloc,1:nrloc,k,f) = ws_cst(f)
      END WHERE
   ENDDO k_310
   ENDDO f_310

!
!3.2 Camenen
!-----------
!

ELSEIF (iopt_sed_ws.EQ.2.OR.iopt_sed_ws.EQ.3) THEN
   IF (iopt_sed_ws.EQ.2) THEN
!     ---sand
      acam = 24.6; bcam = 0.96; mcam = 1.53
   ELSE
!     ---mud
      acam = 26.8; bcam = 2.11; mcam = 1.19
   ENDIF
   camfac1 = 0.5*(acam/bcam)**(1.0/mcam)
   camfac2 = camfac1**2
   camfac3 = 4.0/(3.0*bcam)

   f_320: DO f=1,nf
   k_320: DO k=1,nz+1
      WHERE (maskatc_int)
         dstar = dp(f)*(deltas(:,:,k,f)*gaccatc(1:ncloc,1:nrloc)&
                     & /kinvisc(1:ncloc,1:nrloc,k)**2)**onethird
         array2dc = (SQRT(camfac2+(camfac3*dstar**3)**(1.0/mcam))-camfac1)**mcam
         wfall(1:ncloc,1:nrloc,k,f) = array2dc*kinvisc(1:ncloc,1:nrloc,k)/dp(f)
      END WHERE
   ENDDO k_320
   ENDDO f_320

!
!3.3 Stokes
!----------
!

ELSEIF (iopt_sed_ws.EQ.4) THEN
   f_330: DO f=1,nf
   k_330: DO k=1,nz+1
      WHERE (maskatc_int)
         wfall(1:ncloc,1:nrloc,k,f) = deltas(:,:,k,f)*gaccatc(1:ncloc,1:nrloc)*&
                                    & dp(f)**2/(18.0*kinvisc(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_330
   ENDDO f_330

!
!3.4 Soulsby
!-----------
!

ELSEIF (iopt_sed_ws.EQ.5) THEN
   f_340: DO f=1,nf
   k_340: DO k=1,nz+1
      WHERE (maskatc_int)
         dstar = dp(f)*(deltas(:,:,k,f)*gaccatc(1:ncloc,1:nrloc)&
                     & /kinvisc(1:ncloc,1:nrloc,k)**2)**onethird
         array2dc = SQRT(10.36**2+1.049*dstar**3)-10.36
         wfall(1:ncloc,1:nrloc,k,f) = array2dc*kinvisc(1:ncloc,1:nrloc,k)/dp(f)
      END WHERE
   ENDDO k_340
   ENDDO f_340

!
!3.5 Zhang and Xie (1993)
!------------------------
!

ELSEIF (iopt_sed_ws.EQ.6) THEN
   f_350: DO f=1,nf
   k_350: DO k=1,nz+1
      WHERE (maskatc_int)
         array2dc = (13.95/dp(f))*kinvisc(1:ncloc,1:nrloc,k)
         wfall(1:ncloc,1:nrloc,k,f) = SQRT(array2dc**2+&
             & 1.09*deltas(:,:,k,f)*gaccatc(1:ncloc,1:nrloc)*dp(f)) - array2dc
      END WHERE
   ENDDO k_350
   ENDDO f_350

ENDIF

!
!4. Effect of flocculation
!-------------------------
!

IF (iopt_sed_ws_floc.GT.0) THEN

!
!4.1 Shear rate
!--------------
!

   IF (iopt_sed_ws_floc.EQ.1.OR.iopt_sed_ws_floc.EQ.3) THEN
      k_410: DO k=2,nz
         WHERE (maskatc_int)
            array3dw(:,:,k) = SQRT(dissip(1:ncloc,1:nrloc,k)&
                                & /kinvisc(1:ncloc,1:nrloc,k))
         END WHERE
      ENDDO k_410
   ENDIF

!
!4.2 Flocculation correction factor
!----------------------------------
!

   f_420: DO f=1,nf
   k_420: DO k=2,nz

      flocVL = 1.0; flocVR = 1.0
!     ---Van Leussen (1994)
      IF (iopt_sed_ws_floc.EQ.1.OR.iopt_sed_ws_floc.EQ.3) THEN
         WHERE (maskatc_int)
            flocVL = (1.0+a_leussen*array3dw(:,:,k))/&
                   & (1.0+b_leussen*(array3dw(:,:,k)**2))
         ELSEWHERE
            flocVL = 0.0
         END WHERE
     ENDIF
!    ---Van Rijn (2007) 
     IF (iopt_sed_ws_floc.EQ.2.OR.iopt_sed_ws_floc.EQ.3) THEN
        xmin = 10**(floc_VR_min**(1.0/alpha_VR)-4.0)
        WHERE (maskatc_int)
           array2dc = 2.0*ctot(1:ncloc,1:nrloc,k)/cgel
        ELSEWHERE
           array2dc = 0.0
        END WHERE
        WHERE (array2dc.GT.xmin)
           flocVR = (4.0 + LOG10(array2dc))**alpha_VR
           flocVR = MAX(flocVR,floc_VR_min)
           flocVR = MIN(flocVR,floc_VR_max)
        END WHERE
        WHERE (.NOT.maskatc_int)
           flocVR = 0.0
        END WHERE
      ENDIF
!     ---apply
      SELECT CASE (iopt_sed_ws_floc)
         CASE(1); wfall(1:ncloc,1:nrloc,k,f) = flocVL*wfall(1:ncloc,1:nrloc,k,f)
         CASE(2); wfall(1:ncloc,1:nrloc,k,f) = flocVR*wfall(1:ncloc,1:nrloc,k,f)
         CASE(3)
            wfall(1:ncloc,1:nrloc,k,f) = flocVL*flocVR*&
                                       & wfall(1:ncloc,1:nrloc,k,f)
      END SELECT

   ENDDO k_420
   ENDDO f_420

ENDIF

!
!5. Hindered settling
!--------------------
!
!5.1 Richardson & Zaki (sand)
!----------------------------
!

IF (iopt_sed_ws_hindset.EQ.1) THEN

   k_511: DO k=1,nz+1
      WHERE (maskatc_int)
         array3dw(:,:,k) = MAX((1.0-ctotatw(:,:,k))**n_RichZaki,0.0) 
      END WHERE
   ENDDO k_511

    f_512: DO f=1,nf
    k_512: DO k=1,nz+1
       WHERE (maskatc_int)
          wfall(1:ncloc,1:nrloc,k,f) = array3dw(:,:,k)*&
                                     & wfall(1:ncloc,1:nrloc,k,f)
       END WHERE
    ENDDO k_512
    ENDDO f_512

!
!5.2 Winterwerp & van Kesteren (mud)
!-----------------------------------
!

ELSEIF (iopt_sed_ws_hindset.EQ.2) THEN

   k_521: DO k=1,nz+1
      WHERE (maskatc_int)
         array2dc = MAX(ctotatw(1:ncloc,1:nrloc,k)/cgel,0.0)
         array3dw(:,:,k) = MAX((1.0-array2dc)*(1-ctotatw(:,:,k))/&
                         & (1.0+2.5*array2dc),0.0)
      END WHERE
   ENDDO k_521

   f_522: DO f=1,nf
   k_522: DO k=1,nz+1
      WHERE (maskatc_int)
         wfall(1:ncloc,1:nrloc,k,f) = array3dw(:,:,k)*wfall(1:ncloc,1:nrloc,k,f)
      END WHERE
   ENDDO k_522
   ENDDO f_522
ENDIF

!
!6. Upper limit in shallow waters
!--------------------------------
!

IF (iopt_sed_ws_lim.EQ.1) THEN

   f_610: DO f=1,nf
   k_610: DO k=2,nz
      WHERE (maskatc_int)
         array2dc = gscoordatw(1:ncloc,1:nrloc,k)*deptotatc(1:ncloc,1:nrloc)
         wfall(1:ncloc,1:nrloc,k,f) = MIN(wfall(1:ncloc,1:nrloc,k,f),&
                                        & wslimfac*array2dc/delt3d)
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

DEALLOCATE (array2dc,array3dw)
IF (iopt_sed_ws.GT.1) DEALLOCATE (deltas)
IF (iopt_sed_ws.EQ.2.OR.iopt_sed_ws.EQ.3.OR.iopt_sed_ws.EQ.5) THEN
   DEALLOCATE (dstar)
ENDIF
IF (iopt_sed_ws_floc.GT.0) DEALLOCATE (flocVL,flocVR)
IF (iopt_sed_ws_floc.GT.0.OR.iopt_sed_ws_hindset.GT.0) DEALLOCATE (ctotatw)

CALL log_timer_out()


RETURN

END SUBROUTINE sediment_settling_velocity

!=======================================================================

SUBROUTINE thetastar_engelund_hansen(theta,mask,thetastar,itype)
!************************************************************************
!
! *thetastar_engelund_hansen* Shield parameter used in Engelund & Hansen (1967)
!                             formula
!
! Author - Boudewijn Decrop (IMDC)
!
! Version - @(COHERENS)Sediment_Equations.f90  V2.10.2
!
! Description -
!
! Calling program - sediment_suspendedload, suspended_totalload
!
! Module calls - 
!
!************************************************************************
!
USE gridpars
USE iopars
USE sedswitches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE
!
!*Arguments
!
INTEGER, INTENT(IN) :: itype
REAL, INTENT(IN), DIMENSION(ncloc,nrloc) :: theta
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: thetastar
LOGICAL, INTENT(IN), DIMENSION(ncloc,nrloc) :: mask
!
! Name      Type    Purpose
!-------------------------------------------------------------------------------
!*theta*     REAL    Shields parameter                           
!*mask*      LOGICAL Mask array to exclude dry grid nodes
!*thetastar* REAL    Modified Shields parameter for Engelund-Hansen formula
!*itype*     INTEGER Type of formulation
!             = 1 => standard formulation
!             = 2 => using Chollet-Cunge (1979)
!
!*******************************************************************************
!


procname(pglev+1) = 'thetastar_engelund_hansen'
CALL log_timer_in()

IF (itype.EQ.1) THEN
   WHERE (mask) thetastar = theta
ELSEIF (itype.EQ.2) THEN
   WHERE (theta.LT.0.06) thetastar = 0.0
   WHERE (theta.GE.0.06.AND.theta.LT.0.384)
      thetastar = SQRT(2.5*(theta-0.06))
   END WHERE
   WHERE (theta.GE.0.384.AND.theta.LE.1.08)
      thetastar = 1.065*theta**(0.176)
   END WHERE
   WHERE (theta.GT.1.08)  thetastar = theta
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE thetastar_engelund_hansen
