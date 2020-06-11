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
! *Surface_Fluxes* Meteo input and surface fluxes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Fluxes.f90  V2.11.2
!
! $Date: 2018-04-17 12:38:28 +0200 (Tue, 17 Apr 2018) $
!
!
! $Revision: 1121 $
!
! Description - 
!
! Reference -
!
! Subroutines - exchange_coefs, heat_flux, meteo_input, salinity_flux,
!               short_wave_radiation, surface_stress, vapour_pressure
!
!************************************************************************
!

!========================================================================

SUBROUTINE exchange_coefs
!************************************************************************
!
! *exchange_coefs* Exchange coefficients and atmospheric boundary layer
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Fluxes.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - meteo_input
!
! External calls - vapour_pressure
!
! Module calls -  complex_polar, error_alloc, Uarr_at_C, Varr_at_C
!
!************************************************************************
!
USE currents  
USE fluxes  
USE grid  
USE gridpars
USE iopars
USE meteo
USE modids
USE physpars
USE switches
USE syspars
USE array_interp, ONLY: Uarr_at_C, Varr_at_C
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: complex_polar
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out  

IMPLICIT NONE

!
!*Local variables
!
INTEGER, PARAMETER :: maxiterations = 10
INTEGER :: i, iter, j, l, npcc
REAL, PARAMETER :: alpha = 5.0, beta = 16.0, epstol = 1.0E-06, &
                 & ratmol = 0.62197, ratmol1 = 1.0-ratmol
REAL :: ckar2, decr, dfunc, fac, func, psih, psim, tau0, tau, tau2, tkelv, &
      & tvir, xchar, xsi, xvar
REAL, DIMENSION(3) :: xcoef, xcoef_old
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2, qspec_air, &
                                         & qspec_sur, vappres_sur, wlim
REAL, DIMENSION(5), PARAMETER :: adkon = (/0.0,0.771,0.867,1.2,0.0/),&
                               & aekon = (/0.0,0.969,1.18,1.196,1.68/),&
                               & ahkon = (/0.0,0.927,1.15,1.17,1.652/),&
                               & bdkon = (/1.08,0.0858,0.0667,0.025,0.073/),&
                               & bhkon = (/1.185,0.0546,0.01,0.0075,-0.017/),&
                               & bekon = (/1.23,0.0521,0.01,0.008,-0.016/),&
                               & cekon = (/0.0,0.0,0.0,-0.0004,0.0/),&
                               & chkon = (/0.0,0.0,0.0,-0.00045,0.0/),&
                               & pdkon = (/-0.15,1.0,1.0,1.0,1.0/),&
                               & pekon = (/-0.16,1.0,1.0,1.0,1.0/),&
                               & phkon = (/-0.157,1.0,1.0,1.0,1.0/)


procname(pglev+1) = 'exchange_coefs'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (array2d1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2d2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (qspec_air(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qspec_air',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (qspec_sur(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qspec_sur',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (vappres_sur(ncloc,nrloc),STAT=errstat)
CALL error_alloc('vappres_sur',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (wlim(ncloc,nrloc),STAT=errstat)
CALL error_alloc('wlim',2,(/ncloc,nrloc/),kndrtype)

!
!2. Atmospheric variables
!------------------------
!
!---wind speed
CALL Uarr_at_C(uvel(1:ncloc+1,1:nrloc,nz),array2d1,1,1,(/1,1,1/),&
            & (/ncloc+1,nrloc,1/),1,iarr_uvel,.TRUE.)
array2d1 = uwindatc(1:ncloc,1:nrloc) - array2d1
CALL Varr_at_C(vvel(1:ncloc,1:nrloc+1,nz),array2d2,1,1,(/1,1,1/),&
           & (/ncloc,nrloc+1,1/),1,iarr_vvel,.TRUE.)
array2d2 = vwindatc(1:ncloc,1:nrloc) - array2d2
CALL complex_polar(array2d1,array2d2,xamp=windatc)

!---air/sea differences
IF (iopt_meteo_heat.EQ.1) THEN
   tempdif = rmaskatc*(airtemp-sst)
   CALL vapour_pressure(vappres_air,'A')
   WHERE (maskatc_int)
      qspec_air = ratmol/(0.01*atmpres(1:ncloc,1:nrloc)/vappres_air-ratmol1)
   END WHERE
   CALL vapour_pressure(vappres_sur,'S')
   WHERE (maskatc_int)
      qspec_sur = ratmol/(0.01*atmpres(1:ncloc,1:nrloc)/vappres_sur-ratmol1)
      qspecdif = rmaskatc*(qspec_air-qspec_sur)
   END WHERE
ENDIF

!
!3. Surface exchange coefficients (neutral case)
!-----------------------------------------------
!
!3.1 Drag coefficient
!--------------------
!

IF (iopt_meteo_stres.EQ.1) THEN

   SELECT CASE (iopt_sflux_pars)

!     ---generic
      CASE (1)
         cds = MERGE(cdspars(1),cdspars(2)+cdspars(3)*windatc,&
                   & windatc.LT.cdspars(4))
      
!     ---Large and Pond (1981)
      CASE (2)
         cds = MERGE(0.00114,0.00049+0.000065*windatc,&
                   & windatc.LE.10.0)

!     ---Smith and Banke (1975)
      CASE (3); cds = 0.00063 + 0.000066*windatc

!     ---Geernaert et al. (1986)
      CASE (4); cds = 0.00043 + 0.000097*windatc

!     ---Kondo (1975)      
      CASE (5)
         wlim = MAX(0.3,windatc)
         j_311: DO j=1,nrloc
         i_311: DO i=1,ncloc
            IF (wlim(i,j).LT.2.2) THEN
               l = 1
            ELSEIF (wlim(i,j).LT.5.0) THEN
               l = 2
            ELSEIF (wlim(i,j).LT.8.0) THEN
               l = 3
            ELSEIF (wlim(i,j).LT.25.0) THEN
               l = 4
            ELSE
               l = 5
            ENDIF
            cds(i,j) = 0.001*(adkon(l)+bdkon(l)*(wlim(i,j)**pdkon(l)))
         ENDDO i_311
         ENDDO j_311
   
!     ---Wu (1980)      
      CASE (6)
         array2d1 = 10.0**(0.137*windatc-0.616)
         cds = 0.0012*(array2d1**0.15)

!     ---Charnock
      CASE (7)
         wlim = MAX(0.1,windatc)
         j_312: DO j=1,nrloc
         i_312: DO i=1,ncloc
            cds(i,j) = 0.001
            xchar = LOG(zref_wind*gaccatc(i,j)/ccharno)
            decr = 1.0
            DO WHILE (ABS(decr).GT.epstol)
               func = LOG(cds(i,j)) + ckar/SQRT(cds(i,j)) + &
                        & 2.0*LOG(wlim(i,j)) - xchar
               dfunc = 1.0/cds(i,j) - 0.5*ckar/cds(i,j)**1.5
               decr = func/dfunc
               cds(i,j) = ABS(cds(i,j)-decr)
            ENDDO
         ENDDO i_312
         ENDDO j_312
      
!     ---North Sea      
      CASE (8)
         j_313: DO j=1,nrloc
         i_313: DO i=1,ncloc
            IF (windatc(i,j).LE.5.0) THEN
               cds(i,j) = 0.000621
            ELSEIF (windatc(i,j).GE.19.22) THEN
               cds(i,j) = 0.002764
            ELSE
               cds(i,j) = -0.0001325 + 0.000151*windatc(i,j)
            ENDIF
         ENDDO i_313
         ENDDO j_313
       
!     ---Large and Yeager (2004)
      CASE (9)
         wlim = MAX(0.5,windatc)
         cds = 0.001*(2.7/wlim+0.142+0.0764*wlim)
      
   END SELECT

ENDIF

!
!3.2 Heat exchange coefficients
!------------------------------
!

IF (iopt_meteo_heat.EQ.1) THEN

   SELECT CASE (iopt_sflux_pars)

!     ---generic
      CASE (1)
         ces = MERGE(ces_scst,ces_ucst,tempdif.GE.0.0)
         chs = MERGE(chs_scst,chs_ucst,tempdif.GE.0.0)
      
!     ---Large and Pond (1981)
      CASE (2)
         ces = 0.00115
         chs = MERGE(0.00066,0.00113,tempdif.GE.0.0)

!     ---Smith and Banke (1975)
      CASE (3)
         ces = 0.0015; chs = ces

!     ---Geernaert et al. (1986)
      CASE (4)
         ces = MERGE(ces_scst,ces_ucst,tempdif.GE.0.0)
         chs = MERGE(chs_scst,chs_ucst,tempdif.GE.0.0)

!     ---Kondo (1975)      
      CASE (5)
         j_321: DO j=1,nrloc
         i_321: DO i=1,ncloc
            wlim(i,j) = MAX(0.3,windatc(i,j))
            IF (wlim(i,j).LT.2.2) THEN
               l = 1
            ELSEIF (wlim(i,j).LT.5.0) THEN
               l = 2
            ELSEIF (wlim(i,j).LT.8.0) THEN
               l = 3
            ELSEIF (wlim(i,j).LT.25.0) THEN
               l = 4
            ELSE
               l = 5
            ENDIF
            fac = (wlim(i,j)-8.0)**2
            ces(i,j) = 0.001*(aekon(l)+bekon(l)*(wlim(i,j)**pekon(l))&
                           & +cekon(l)*fac)
            chs(i,j) = 0.001*(ahkon(l)+bhkon(l)*(wlim(i,j)**phkon(l))&
                           & +chkon(l)*fac)
         ENDDO i_321
         ENDDO j_321

!     ---Wu (1980)
      CASE (6)
         array2d1 = 10.0**(0.137*windatc-0.616)
         ces = 0.001*(array2d1**0.11)
         chs = ces

!     ---Charnock
      CASE (7)
         ces = MERGE(ces_scst,ces_ucst,tempdif.GE.0.0)
         chs = MERGE(chs_scst,chs_ucst,tempdif.GE.0.0)
         
!     ---North Sea
      CASE (8)
         ces = MERGE(ces_scst,ces_ucst,tempdif.GE.0.0)
         chs = MERGE(chs_scst,chs_ucst,tempdif.GE.0.0)
       
!     ---Large and Yeager (2004)
      CASE (9)
         array2d1 = MERGE(0.0327,0.018,tempdif.GT.0.0)
         ces = 0.0346*SQRT(cds)
         chs = array2d1*SQRT(cds)
      
   END SELECT

ENDIF

!
!4. Surface exchange coefficients (stratified case)
!--------------------------------------------------
!

IF (iopt_sflux_strat.EQ.1) THEN

!
!4.1 Convert to new reference height
!-----------------------------------
!   

   IF (ABS(zref_wind-zref_temp).GT.0.01) THEN
      
      xvar = LOG(zref_temp/zref_wind)/ckar
      array2d1 = cds
      cds = (cds/(1.0+SQRT(cds)*xvar))**2
      array2d2 = cds/array2d1
      windatc = windatc/array2d2
      array2d2 = SQRT(array2d2)
      array2d1 = SQRT(array2d1)
      ces = ces*array2d2/(1.0+xvar*ces/array2d1)
      chs = chs*array2d2/(1.0+xvar*chs/array2d1)
      
   ENDIF

!   
!4.2 Kondo (1975)
!----------------
!

   IF (iopt_sflux_pars.EQ.5) THEN
      
      array2d1 = -rmaskatc*tempdif/wlim**2
      array2d1 = array2d1*ABS(array2d1)/(ABS(array2d1)+0.01)
      WHERE (tempdif.GE.0.0)
         array2d2 = MERGE(0.0,0.1+0.03*array2d1+0.9*EXP(4.8*array2d1),&
                        & array2d1.LE.-3.3)
         cds = cds*array2d2
         ces = ces*array2d2
         chs = chs*array2d2
      ELSEWHERE
         array2d2 = 1.0+0.47*SQRT(array2d1)
         cds = cds*array2d2
         array2d2 = 1.0+0.63*SQRT(array2d1)
         ces = ces*array2d2
         chs = chs*array2d2
      END WHERE

!
!4.3 General case
!---------------- 
!

   ELSE

      ckar2 = ckar*ckar

      j_430: DO j=1,nrloc
      i_430: DO i=1,ncloc
         wlim(i,j) = MAX(1.0,windatc(i,j))
         xcoef(1) = cds(i,j); xcoef(2) = ces(i,j); xcoef(3) = chs(i,j)
         tau0 = ckar/SQRT(xcoef(1)) + LOG(xcoef(1))
         tkelv = airtemp(i,j) + celtokelv
         tvir = tkelv*(1.0+0.61*qspec_air(i,j))
         fac = gaccatc(i,j)*ckar*zref_temp/tvir/wlim(i,j)/wlim(i,j)

         iter_431: DO iter=1,maxiterations
            
!           ---height normalised with MO length
            xvar = fac/(xcoef(1)**1.5)
            xsi = xvar*(xcoef(2)*tempdif(i,j)+0.61*tkelv*xcoef(3)*qspecdif(i,j))
!           ---flux profiles
            IF (xsi.GE.0.0) THEN
               psim = -alpha*LOG(xsi+1)
               psih = psim
            ELSE
               xvar = SQRT(1.0-beta*xsi)
               psih = 2.0*LOG(0.5+0.5*xvar)
               xvar = SQRT(xvar)
               psim = 2.0*LOG(0.5+0.5*xvar) + LOG(0.5+0.5*xvar*xvar) - &
                    & 2.0*ATAN(xvar) + halfpi
            ENDIF
!           ---work space variables            
            tau = tau0 - LOG(xcoef(1))
            tau2 = tau - psim
!           ---update values
            xcoef_old = xcoef
            xcoef(1) = 1.0/(tau/ckar-psim/ckar)**2
            xcoef(2) = ces(i,j)/(1.0-psim/tau-psih*ces(i,j)*tau2/ckar2)
            xcoef(3) = chs(i,j)/(1.0-psim/tau-psih*chs(i,j)*tau2/ckar2)
!           ---check convergence            
            IF (ALL(ABS(xcoef-xcoef_old).LE.epstol)) EXIT iter_431

         ENDDO iter_431
         cds(i,j) = xcoef(1); ces(i,j) = xcoef(2); chs(i,j) = xcoef(3)
      
      ENDDO i_430
      ENDDO j_430

   ENDIF

ENDIF

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (array2d1,array2d2,qspec_air,qspec_sur,vappres_sur,wlim)

CALL log_timer_out(npcc,itm_bconds)


RETURN

END SUBROUTINE exchange_coefs

!========================================================================

SUBROUTINE heat_flux
!************************************************************************
!
! *heat_flux* (Non-solar) surface heat fluxes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Fluxes.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model, temperature_equation
!
! External calls -
!
! Module calls - error_alloc
!
!************************************************************************
!
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc
REAL, PARAMETER :: cpair = 1005.0, emiscoef = 0.98
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: cloudcorr, sstkelv

IF (.NOT.metstep) RETURN

procname(pglev+1) = 'heat_flux'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!
!---allocate
IF (nt.EQ.0.AND.iopt_sflux_qlong.EQ.2) THEN
   IF (ALLOCATED(cloudcorr)) DEALLOCATE (cloudcorr)
   ALLOCATE (cloudcorr(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('cloudcorr',2,(/ncloc,nrloc/),kndrtype)
   cloudcorr = 0.00422*ABS(gylat(1:ncloc,1:nrloc)) + 0.5
ENDIF
ALLOCATE (sstkelv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sstkelv',2,(/ncloc,nrloc/),kndrtype)

!---convert sst to degrees Kelvin
sstkelv = sst + celtokelv

!
!1. Latent/sensible heat fluxes
!------------------------------
!
   
evaporation = -rho_air*rmaskatc*ces*windatc*qspecdif
qlatent = (2.5008*1.0E+06-2300.0*sst)*rmaskatc*evaporation
qsensible = -(rho_air*cpair)*rmaskatc*chs*windatc*tempdif

!
!2. Long wave radiative flux
!---------------------------
!
!2.1 Gill (1982)
!---------------
!
   
IF (iopt_sflux_qlong.EQ.1) THEN
   qlwave = rmaskatc*(emiscoef*stefbolz)*sstkelv**4*&
          & (0.39-0.05*SQRT(vappres_air))*(1.0-0.6*cloud_cover**2)

!   
!2.2 Clark et al. (1974)
!-----------------------
!

ELSEIF (iopt_sflux_qlong.EQ.2) THEN
   qlwave = rmaskatc*(emiscoef*stefbolz)*sstkelv**3*(sstkelv*&
        & (0.39-0.05*SQRT(vappres_air))*(1.0-cloudcorr*cloud_cover**2)-&
        & 4.0*tempdif)

!
!2.3 Bignami et al. (1995)
!-------------------------
!
   
ELSEIF (iopt_sflux_qlong.EQ.3) THEN
   qlwave = rmaskatc*stefbolz*(emiscoef*sstkelv**4-&
         & (airtemp+celtokelv)**4*(0.653+0.00535*vappres_air)*&
         & (1.0+0.1762*cloud_cover**2))
ENDIF

!
!3. Non-solar heat flux
!----------------------
!

qnonsol = qlatent + qsensible + qlwave

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (sstkelv)

CALL log_timer_out(npcc,itm_bconds)


RETURN

END SUBROUTINE heat_flux

!========================================================================

SUBROUTINE meteo_input
!************************************************************************
!
! *meteo_input* Update surface meteo input
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Fluxes.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! External calls - combine_particle_meteo_data, define_surface_input_grid,
!                  surface_exchange_MO, update_surface_data
!
! Module calls - error_alloc, error_alloc_struc, exchange_mod,
!                hrelativecoords_init, set_modfiles_atts
!
!************************************************************************
!
USE datatypes
USE density
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE modids
USE optics
USE paralpars
USE physpars
USE switches
USE syspars
USE timepars
USE datatypes_init, ONLY: hrelativecoords_init
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE grid_routines, ONLY: rotate_vec
USE modvars_routines, ONLY: set_modfiles_atts
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: ivar, nhtype, npcc
INTEGER, SAVE :: novars, n1dat, n1grd, n2dat, n2grd
INTEGER (KIND=kndilong), SAVE, DIMENSION(2) :: nosecsdat
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: angle
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: meteovals
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: meteoindat
TYPE (HRelativeCoords), SAVE, ALLOCATABLE, DIMENSION(:,:) :: meteogrid

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*nhtype*     INTEGER Type of data grid (0/1/2)
!*n1dat*      INTEGER X-dimension of data grid
!*n2dat*      INTEGER Y-dimension of data grid
!*novars*     INTEGER Number of input parameters
!*nosecsdat*  LONGINT No. of seconds between start of simulation and old/new
!                     data time
!*meteoindat* REAL    Meteo data at old and new data time
!*meteovals*  REAL    Meteo data interpolated in time and on model grid
!*meteogrid*  DERIVED Model grid relative coordinates with respect to meteo
!                     grid
!
!------------------------------------------------------------------------------
!
!1. Return if no update required
!-------------------------------
!

IF (.NOT.metstep) RETURN

procname(pglev+1) = 'meteo_input'
CALL log_timer_in(npcc)

nhtype = surfacegrids(igrd_meteo,1)%nhtype

!
!2. Initialise on first call
!---------------------------
!

IF (nt.EQ.0) THEN

!
!2.1 Surface grid
!----------------
!
!  ---grid parameters
   n1dat = surfacegrids(igrd_meteo,1)%n1dat
   n2dat = surfacegrids(igrd_meteo,1)%n2dat
   n1grd  = MERGE(ncloc,0,nhtype.GT.0.AND.nhtype.LT.4)
   n2grd  = MERGE(nrloc,0,nhtype.GT.0.AND.nhtype.LT.4)

!  ---allocate coordinate structure
   ALLOCATE (meteogrid(n1grd,n2grd),STAT=errstat)
   CALL error_alloc_struc('meteogrid',2,(/n1grd,n2grd/),'HRelativeCoords')
   CALL hrelativecoords_init(meteogrid,.FALSE.)

!  ---define surface grid
   IF (nhtype.GT.0.AND.nhtype.LT.4) THEN
      CALL define_surface_input_grid(igrd_meteo,1,meteogrid)
   ENDIF

!
!2.2 Number of meteo variables
!-----------------------------
!

   CALL set_modfiles_atts(io_metsur,1,1)
   novars = modfiles(io_metsur,1,1)%novars

!
!2.3 Allocate memory for data arrays
!-----------------------------------
!

   ALLOCATE (meteoindat(n1dat,n2dat,novars,2),STAT=errstat)
   CALL error_alloc('meteoindat',4,(/n1dat,n2dat,novars,2/),kndrtype)
   meteoindat = 0.0
   ALLOCATE (meteovals(ncloc,nrloc,novars),STAT=errstat)
   CALL error_alloc('meteovals',3,(/ncloc,nrloc,novars/),kndrtype)
   meteovals = 0.0
   ALLOCATE (angle(novars),STAT=errstat)
   CALL error_alloc('angle',1,(/novars/),kndlog)
   angle = .FALSE.

ENDIF

!
!3. Store old wind values (particle model only)
!----------------------------------------------
!

IF (iopt_part_model.GT.0) THEN
   uwindatc_old = uwindatc
   vwindatc_old = vwindatc
ENDIF

!
!4. Update meteo input
!---------------------
!

CALL update_surface_data(io_metsur,1,meteoindat,meteovals,meteogrid,&
                       & n1dat,n2dat,novars,n1grd,n2grd,nosecsdat,angle)

!
!5. Store data into surface arrays
!---------------------------------
!

ivar = 0

!---parameters for surface stress
IF (iopt_meteo_stres.EQ.1) THEN
   IF (iopt_meteo_data.EQ.1) THEN
      ivar = ivar + 1
      uwindatc(1:ncloc,1:nrloc) = meteovals(:,:,ivar)
      ivar = ivar + 1
      vwindatc(1:ncloc,1:nrloc) = meteovals(:,:,ivar)
!     --case of rotated grid
      IF (rotate_gvecs) THEN
         CALL rotate_vec(meteovals(:,:,ivar-1),meteovals(:,:,ivar),&
                       & uwindatc(1:ncloc,1:nrloc),vwindatc(1:ncloc,1:nrloc),&
                       & gangleatc(1:ncloc,1:nrloc),1)
      ENDIF
   ELSEIF (iopt_meteo_data.EQ.2) THEN
      ivar = ivar + 1
      usstresatc(1:ncloc,1:nrloc) = meteovals(:,:,ivar)
      ivar = ivar + 1
      vsstresatc(1:ncloc,1:nrloc) = meteovals(:,:,ivar)
!     --case of rotated grid
      IF (rotate_gvecs) THEN
         CALL rotate_vec(meteovals(:,:,ivar-1),meteovals(:,:,ivar),&
              & usstresatc(1:ncloc,1:nrloc),&
              & vsstresatc(1:ncloc,1:nrloc),gangleatc(1:ncloc,1:nrloc),1)
      ENDIF
   ENDIF
ENDIF

IF ((iopt_part_model.EQ.2.AND.modelid.EQ.modelidpart).OR.&
   & iopt_part_model.GT.2) THEN
   GOTO 1000
ENDIF

!---atmospheric pressure
IF (iopt_meteo_pres.EQ.1) THEN
   ivar = ivar + 1
   atmpres(1:ncloc,1:nrloc) = meteovals(:,:,ivar)
ENDIF

!---parameters for heat flux
IF (iopt_meteo_heat.EQ.1) THEN
   IF (iopt_meteo_data.EQ.1) THEN
      ivar = ivar + 1
      airtemp = meteovals(:,:,ivar)
      sst = temp(1:ncloc,1:nrloc,nz)
      ivar = ivar + 1
      relhum = meteovals(:,:,ivar)
      ivar = ivar + 1
      cloud_cover = meteovals(:,:,ivar)
   ELSEIF (iopt_meteo_data.EQ.2) THEN
      ivar = ivar + 1
      qnonsol = meteovals(:,:,ivar)
      ivar = ivar + 1
      qrad = meteovals(:,:,ivar)
   ENDIF
ENDIF
   
!---(evaporation)/precipitation
IF (iopt_meteo_precip.EQ.1) THEN
   ivar = ivar + 1
   precipitation = meteovals(:,:,ivar)
ELSEIF (iopt_meteo_precip.EQ.2) THEN
   ivar = ivar + 1
   evapminprec = meteovals(:,:,ivar)
ENDIF

!
!6. Exchange halos
!-----------------
!

IF (iopt_MPI.EQ.1) THEN

   IF (iopt_meteo_stres.EQ.1) THEN
      IF (iopt_meteo_data.EQ.1) THEN
         lbounds = 0; nhexch = 1
         CALL exchange_mod(uwindatc,lbounds,nhexch,iarr_uwindatc)
         CALL exchange_mod(vwindatc,lbounds,nhexch,iarr_vwindatc)
      ELSEIF (iopt_meteo_data.EQ.2) THEN
         lbounds = (/0,1/); nhexch = (/1,0,0,0/)
         CALL exchange_mod(usstresatc,lbounds,nhexch,iarr_usstresatc)
         lbounds = (/1,0/); nhexch = (/0,0,1,0/)
         CALL exchange_mod(vsstresatc,lbounds,nhexch,iarr_vsstresatc)
      ENDIF
   ENDIF
   IF (iopt_meteo_pres.EQ.1) THEN
      lbounds = 0; nhexch = (/1,0,1,0/)
      CALL exchange_mod(atmpres,lbounds,nhexch,iarr_atmpres)
   ENDIF

ENDIF

!
!7. Surface drag and heat exchange coefficients
!----------------------------------------------
!

IF (iopt_meteo_data.EQ.1.AND.iopt_sflux_pars.GT.0) CALL exchange_coefs

!
!8. Deallocate arrays
!--------------------
!

1000 CONTINUE

IF (cold_start.OR.&
  & nt+ABS(modfiles(io_metsur,1,1)%tlims(3)).GT.nstep) THEN
   DEALLOCATE (angle,meteogrid,meteoindat,meteovals)
ENDIF

CALL log_timer_out(npcc,itm_meteo)


RETURN

END SUBROUTINE meteo_input

!========================================================================

SUBROUTINE salinity_flux
!************************************************************************
!
! *salinity_flux* Surface salinity flux
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Fluxes.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - salinity_equation
!
! External calls -
!
! Module calls -
!
!************************************************************************
!
USE density
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE physpars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc


procname(pglev+1) = 'salinity_flux'
CALL log_timer_in(npcc)

IF (iopt_meteo_precip.EQ.1) THEN
   WHERE (maskatc_int)
      ssalflux = sal(1:ncloc,1:nrloc,nz)*(evaporation-precipitation)&
               & /(1.0-0.001*sal(1:ncloc,1:nrloc,nz))/density_ref

   END WHERE
ELSEIF (iopt_meteo_precip.EQ.2) THEN
   WHERE (maskatc_int)
      ssalflux = sal(1:ncloc,1:nrloc,nz)*evapminprec&
               & /(1.0-0.001*sal(1:ncloc,1:nrloc,nz))/density_ref
   END WHERE
ENDIF

CALL log_timer_out(npcc,itm_bconds)


RETURN

END SUBROUTINE salinity_flux

!========================================================================

SUBROUTINE short_wave_radiation
!************************************************************************
!
! *short_wave_radiation* Short wave solar radiance
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Fluxes.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - temperature_equation
!
! External calls -
!
! Module calls - day_number, error_alloc, leap_year
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE meteo
USE optics
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: day_number, leap_year, log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: l, npcc
REAL, PARAMETER :: albedo_os = 0.06, doptcld = 6.35, solar_cst = 1367.0, &
                 & tau = 0.7
REAL :: daynum, dayrad1, dayrad2, decsol, hcorr, hourday, porb, solrad
REAL, DIMENSION(0:2) :: aorb = (/1.00011,0.034221,0.000719/)
REAL, DIMENSION(2) :: borb = (/0.00128,0.000077/)
REAL, DIMENSION(0:3) :: adec = (/0.006918,-0.399912,-0.006758,-0.002697/)
REAL, DIMENSION(3) :: bdec = (/0.070257,0.000907,0.00148/)
REAL, DIMENSION(3) :: aeq = (/0.0072,-0.0528,-0.0012/),&
                    & beq = (/-0.1229,-0.1565,-0.0041/)
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: albedo_cs, altmax, array2dc, &
                                         & cloudcorr,coszen, fabsorb, hoursol

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*albedo_cs*  REAL    Sea surface albedo at under clear sky conditions
!*albedo_os*  REAL    Sea surface albedo at under overcast sky conditions
!*altmax*     REAL    Sun's altitude at noon                             [deg]
!*cloudcorr*  REAL    Correction factor for cloudiness   
!*coszen*     REAL    Cosinus of the sun's zenith angle
!*daynum*     REAL    Day number
!*decsol*     REAL    Sun's declination                                  [rad]
!*doptcld*    REAL    Dimensionless cloud optical depth
!*fabsorb*    REAL    Fraction of solar radiation absorbed at the sea surface
!*hcorr*      REAL    Difference (in hours) between real and mean sun  [hours]
!                                                                      [hours]
!*hourday*    REAl    Hour of the day (zonal time)                     [hours]
!*hoursol*    REAL    Sun's hour angle                                   [rad]
!*solar_cst*  REAL    Solar constant (anuual mean of solar radiation at top
!                     of the atmosphere                                [W/m^2]
!*solrad*     REAL    Solar constant corrected for Earth motion around
!                     the sun                                          [W/m^2]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'short_wave_radiation'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (albedo_cs(ncloc,nrloc),STAT=errstat)
CALL error_alloc('albedo_cs',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (coszen(ncloc,nrloc),STAT=errstat)
CALL error_alloc('coszen',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (fabsorb(ncloc,nrloc),STAT=errstat)
CALL error_alloc('fabsorb',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (hoursol(ncloc,nrloc),STAT=errstat)
CALL error_alloc('hoursol',2,(/ncloc,nrloc/),kndrtype)
IF (iopt_sflux_qshort.EQ.1) THEN
   ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ENDIF
IF (iopt_sflux_qshort.LE.2) THEN
   ALLOCATE (altmax(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('altmax',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (cloudcorr(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('cloudcorr',2,(/ncloc,nrloc/),kndrtype)
ENDIF

!
!2. Initialise parameters
!------------------------
!
!---day number
CALL day_number(IdateTime,daynum)
dayrad1 = twopi*(daynum-1.0)/(365.0+leap_year(IdateTime(1)))
dayrad2 = twopi*daynum/(365.0+leap_year(IdateTime(1)))

!---hour of the day
hourday = IdateTime(4) + (IdateTime(5)+IdateTime(6)/60.0)/60.0

!---solar radiation at top of atmospehere
porb = aorb(0)
l_210: DO l=1,2
   porb = porb + aorb(l)*COS(dayrad2) + borb(l)*SIN(dayrad2)
ENDDO l_210
solrad = solar_cst*porb

!---solar declination, time equation
decsol = adec(0); hcorr = 0.0
l_220: DO l=1,3
   decsol = decsol + adec(l)*COS(l*dayrad1) + bdec(l)*SIN(l*dayrad1)
   hcorr = hcorr + aeq(l)*COS(l*dayrad2) + beq(l)*SIN(l*dayrad2)
ENDDO l_220

!---hour angle of the sun
hoursol = gxlon(1:ncloc,1:nrloc) + pi*(hourday-time_zone-12.0+hcorr)/12.0

!---zenith angle
coszen = SIN(decsol)*SIN(gylat(1:ncloc,1:nrloc))+COS(decsol)*&
       & COS(gylat(1:ncloc,1:nrloc))*COS(hoursol)
coszen = MAX(0.0,coszen)

!---clear sky albedo
albedo_cs = 0.05/(1.1*coszen**1.4+0.15)

!---cloud correction factor
IF (iopt_sflux_qshort.LE.2) THEN
!  --solar altitude at noon   
   WHERE (gylat(1:ncloc,1:nrloc).GT.decsol)
      altmax = halfpi + decsol - gylat(1:ncloc,1:nrloc)
   ELSEWHERE
      altmax = halfpi - decsol + gylat(1:ncloc,1:nrloc)
   END WHERE
   altmax = radtodeg*altmax
!  --correction factor
   WHERE (maskatc_int)
      cloudcorr = 1.0-0.62*cloud_cover+0.0019*altmax
   END WHERE
ENDIF
   
!
!3. Solar irradiance
!-------------------
!
!3.1 Reed (1977)
!---------------
!

IF (iopt_sflux_qshort.EQ.1) THEN
   
!  ---attenuation factor
   WHERE (coszen.GT.0.05)
      array2dc = 0.5*tau**(1.0/coszen) + 0.455
   ELSEWHERE
      array2dc = 0.455
   END WHERE

!  ---absorbed fraction
   fabsorb = 1.0 - albedo_os*cloud_cover - albedo_cs*(1.0-cloud_cover)
   
!  ---short-wave radiation
   qrad = solrad*fabsorb*coszen*array2dc*cloudcorr

!   
!3.2 Zillman (1972)
!------------------
!   

ELSEIF (iopt_sflux_qshort.EQ.2) THEN
   
!  ---absorbed fraction
   fabsorb = 1.0 - albedo_os*cloud_cover - albedo_cs*(1.0-cloud_cover)
      
!  ---short-wave radiation
   qrad = solrad*fabsorb*cloudcorr*coszen*coszen/&
        & (1.085*coszen+0.001*vappres_air*(2.7+coszen)+0.1)

!
!3.3 Shine (1984)
!----------------
!   

ELSEIF (iopt_sflux_qshort.EQ.3) THEN
      
!  ---clear sky part)
   qrad = solrad*(1.0-albedo_cs)*(1.0-cloud_cover)*coszen*coszen/&
        & (1.2*coszen+0.001*vappres_air*(1.0+coszen)+0.0455)
         
!  ---add oversky part
   qrad = qrad + cloud_cover*(1.0-albedo_os)*SQRT(coszen)*(53.5+1274.5*coszen)/&
        & (1.0+0.139*doptcld*(1.0-0.935*albedo_os))

ENDIF
   
!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (albedo_cs,coszen,fabsorb,hoursol)
IF (iopt_sflux_qshort.EQ.1) DEALLOCATE (array2dc)
IF (iopt_sflux_qshort.LE.2) DEALLOCATE (altmax,cloudcorr)

CALL log_timer_out(npcc,itm_bconds)


RETURN

END SUBROUTINE short_wave_radiation

!========================================================================

SUBROUTINE surface_stress
!************************************************************************
!
! *surface_stress* Surface stress
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Fluxes.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - current_corr, current_corr_1d, current_2d, initialise_model
!
! External calls - surface_drag_coef
!
! Module calls - error_alloc, exchange_mod, Carr_at_U, Carr_at_V
!
!************************************************************************
!
USE currents  
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE modids
USE physpars
USE switches
USE timepars
USE array_interp, ONLY: Carr_at_U, Carr_at_V
USE paral_comms, ONLY: exchange_mod 
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc
REAL :: fac

   
IF (.NOT.metstep.OR.nt.EQ.nstep) RETURN

procname(pglev+1) = 'surface_stress'
CALL log_timer_in(npcc)

!
!1. Surface stress components
!----------------------------
!
!1.1 C-nodes
!-----------
!

IF (iopt_meteo_data.EQ.1) THEN

   fac = rho_air/density_ref
   usstresatc(1:ncloc,:) = fac*rmaskatc*cds*windatc*&
              & (uwindatc(1:ncloc,1:nrloc)-uvel(1:ncloc,1:nrloc,nz))
   vsstresatc(:,1:nrloc) = fac*rmaskatc*cds*windatc*&
              & (vwindatc(1:ncloc,1:nrloc)-vvel(1:ncloc,1:nrloc,nz))
   
!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      CALL exchange_mod(usstresatc,(/0,1/),(/1,0,0,0/),iarr_usstresatc)
      CALL exchange_mod(vsstresatc,(/1,0/),(/0,0,1,0/),iarr_vsstresatc)
   ENDIF

ENDIF

!
!1.2 U-nodes
!-----------
!

CALL Carr_at_U(usstresatc,usstresatu(1:ncloc,:),1,3,(/0,1,nz/),&
            & (/ncloc,nrloc,nz/),1,iarr_usstresatc,.TRUE.)

!
!1.3 V-nodes
!-----------
!

CALL Carr_at_V(vsstresatc,vsstresatv(:,1:nrloc),1,3,(/1,0,nz/),&
            & (/ncloc,nrloc,nz/),1,iarr_vsstresatc,.TRUE.)

!
!2. Apply reduction factor for flooding/drying scheme
!----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN
   WHERE (node2du(1:ncloc,1:nrloc).EQ.2)
      usstresatu(1:ncloc,:) = alphatu_fld*usstresatu(1:ncloc,:)
   END WHERE
   WHERE (node2dv(1:ncloc,1:nrloc).EQ.2)
      vsstresatv(:,1:nrloc) = alphatv_fld*vsstresatv(:,1:nrloc)
   END WHERE
ENDIF

!
!3. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(usstresatu,(/1,1/),(/0,1,0,0/),iarr_usstresatu)
   CALL exchange_mod(vsstresatv,(/1,1/),(/0,0,0,1/),iarr_vsstresatv)
ENDIF

!
!4. Total stress at C-nodes
!--------------------------
!

IF (iopt_meteo_data.EQ.1) THEN
   sstresatc = fac*rmaskatc*cds*windatc**2
ELSEIF (iopt_meteo_data.EQ.2) THEN
   sstresatc = rmaskatc*SQRT(usstresatc(1:ncloc,:)**2+vsstresatc(:,1:nrloc)**2)
ENDIF

CALL log_timer_out(npcc,itm_bconds)


RETURN

END SUBROUTINE surface_stress

!========================================================================

SUBROUTINE vapour_pressure(vpres,clev)
!************************************************************************
!
! *vapour_pressure* Saturated vapour pressure
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Fluxes.f90  V2.11.2
!
! Description - uses Tetens' formula
!
! Reference -
!
! Calling program - exchange_coefs
!
! External calls - 
!
! Module calls -
!
!************************************************************************
!
USE grid  
USE gridpars  
USE iopars
USE meteo
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
CHARACTER (LEN=1), INTENT(IN) :: clev
REAL, INTENT(OUT), DIMENSION(ncloc,nrloc) :: vpres 

!
! Name   Type  Purpose
!-------------------------------------------------------------------------------
!*clev*  CHAR  Result is evaluated at reference height if equal to 'A', or at
!              the sea surface if equal to 'S'
!*vpres* REAL  Returned saturation vapour pressure                        [mbar]
!
!-------------------------------------------------------------------------------
!


procname(pglev+1) = 'vapour_pressure'
CALL log_timer_in()

!
!1. At reference height
!----------------------
!

IF (clev.EQ.'A') THEN

   WHERE (airtemp.GT.0.0)
      vpres = 6.1078*rmaskatc*relhum*EXP(17.269*airtemp/(airtemp+273.29))
   ELSEWHERE
      vpres = 6.1078*rmaskatc*relhum*EXP(21.875*airtemp/(airtemp+265.49))
   END WHERE

!
!2. At sea surface
!-----------------
!   

ELSEIF (clev.EQ.'S') THEN
   
   WHERE (sst.GT.0.0)
      vpres = 6.1078*rmaskatc*EXP(17.269*sst/(sst+273.29))
   ELSEWHERE
      vpres = 6.1078*rmaskatc*EXP(21.875*sst/(sst+265.49))
   END WHERE

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE vapour_pressure
