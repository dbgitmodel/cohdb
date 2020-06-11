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
! *Density_Equations* Density (temperature and salinity) calculations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Density_Equations.f90  V2.11.2
!
! $Date: 2018-09-28 13:01:37 +0200 (Fri, 28 Sep 2018) $
!
! $Revision: 1188 $
!
! Description -
!
! Routines - convect_adjust, salinity_equation, temperature_equation,
!            heat_optics, equation_of_state, buoyancy_frequency,
!            baroclinic_gradient, hcube_deriv, hcube_fluxes
!
!************************************************************************
!

!========================================================================

SUBROUTINE salinity_equation
!************************************************************************
!
! *salinity_equation* Solve salinity equation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Density_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! External calls - define_profobc_spec, open_boundary_conds_prof,
!                  salinity_flux, scalar_discharge, transport_at_C_3d,
!                  update_dischr_data, update_profobc_data,
!                  update_nest_data_prof
!
! Module calls - error_alloc, exchange_mod
!
!************************************************************************
!
USE currents
USE density
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE modids
USE nestgrids
USE obconds
USE relaxation
USE structures
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: update
LOGICAL, SAVE :: disch_data
INTEGER :: k, npcc
INTEGER, SAVE :: maxprofs, nofiles, noprofs
INTEGER, DIMENSION(4) :: nhexch
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: iprofrlx, itypobu, itypobv, &
                                          & noprofsd, novarsnst
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: indexprof, indexvar, instvars, &
                                            & iprofobu, iprofobv
INTEGER (KIND=kndilong), SAVE, DIMENSION(2) :: nosecsdis
INTEGER (KIND=kndilong), SAVE, ALLOCATABLE, DIMENSION(:,:) :: nosecsprof
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: dissal
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: disdata
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: obcsal3d, salobu, salobv, source
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: profsal


procname(pglev+1) = 'salinity_equation'
CALL log_timer_in(npcc)

!
!1. Initialise open boundary conditions on first call
!----------------------------------------------------
!

IF (nt.EQ.0) THEN
   nofiles = maxdatafiles(io_salobc,1)
   disch_data = iopt_dischr.EQ.1.AND.modfiles(io_dissal,1,1)%defined

!  ---allocate specifier arrays and o.b. profiles
   ALLOCATE (salobu(nobu,nz,1),STAT=errstat)
   CALL error_alloc('salobu',3,(/nobu,nz,1/),kndrtype)
   IF (nobu.GT.0) salobu = float_fill
   ALLOCATE (salobv(nobv,nz,1),STAT=errstat)
   CALL error_alloc('salobv',3,(/nobv,nz,1/),kndrtype)
   IF (nobv.GT.0) salobv = float_fill
   ALLOCATE (iprofrlx(norlxzones),STAT=errstat)
   CALL error_alloc('iprofrlx',1,(/norlxzones/),kndint)
   IF (norlxzones.GT.0) iprofrlx = 0
   IF (iopt_obc_sal.EQ.1) THEN
      ALLOCATE (itypobu(nobu),STAT=errstat)
      CALL error_alloc('itypobu',1,(/nobu/),kndint)
      ALLOCATE (itypobv(nobv),STAT=errstat)
      CALL error_alloc('itypobv',1,(/nobv/),kndint)
      ALLOCATE (iprofobu(nobu,1),STAT=errstat)
      CALL error_alloc('iprofobu',2,(/nobu,1/),kndint)
      ALLOCATE (iprofobv(nobv,1),STAT=errstat)
      CALL error_alloc('iprofobv',2,(/nobv,1/),kndint)
      ALLOCATE (noprofsd(2:nofiles),STAT=errstat)
      CALL error_alloc('noprofsd',1,(/nofiles-1/),kndint)
      ALLOCATE (indexprof(nobu+nobv,2:nofiles),STAT=errstat)
      CALL error_alloc('indexprof',2,(/nobu+nobv,nofiles-1/),kndint)
      ALLOCATE (indexvar(nobu+nobv,2:nofiles),STAT=errstat)
      CALL error_alloc('indexvar',2,(/nobu+nobv,nofiles-1/),kndint)
   ENDIF
   IF (iopt_nests.EQ.1) THEN
      ALLOCATE (novarsnst(nonestsets),STAT=errstat)
      CALL error_alloc('novarsnst',1,(/nonestsets/),kndint)
      novarsnst = 1
      ALLOCATE (instvars(1,nonestsets),STAT=errstat)
      CALL error_alloc('instvars',2,(/1,nonestsets/),kndint)
      instvars = 1
   ENDIF

!  ---define specifier arrays
   IF (iopt_obc_sal.EQ.1) THEN
      CALL define_profobc_spec(io_salobc,itypobu,itypobv,iprofobu,iprofobv,&
                             & iprofrlx,noprofsd,indexprof,indexvar,1,nofiles,&
                             & nobu,nobv)
   ENDIF

!  ---allocate o.b. data arrays
   noprofs = 0; maxprofs = 0
   IF (iopt_obc_sal.EQ.1) THEN
      IF (nofiles.GT.1) noprofs = MAX(MAXVAL(iprofobu),MAXVAL(iprofobv))
      ALLOCATE (obcsal3d(0:noprofs,nz,1),STAT=errstat)
      CALL error_alloc('obcsal3d',3,(/noprofs+1,nz,1/),kndrtype)
      obcsal3d = float_fill
   ENDIF
   IF (nofiles.GT.1) THEN
      maxprofs = MAXVAL(noprofsd)
      ALLOCATE (nosecsprof(2:nofiles,2),STAT=errstat)
      CALL error_alloc('nosecsprof',2,(/nofiles-1,2/),kndilong)
      ALLOCATE (profsal(maxprofs,nz,2:nofiles,2),STAT=errstat)
      CALL error_alloc('profsal',4,(/maxprofs,nz,nofiles-1,2/),kndrtype)
   ENDIF

!  ---allocate discharge data arrays
   IF (disch_data) THEN
      ALLOCATE (dissal(numdis),STAT=errstat)
      CALL error_alloc('dissal',1,(/numdis/),kndrtype)
      ALLOCATE (disdata(numdis,2),STAT=errstat)
      CALL error_alloc('disdata',2,(/numdis,2/),kndrtype)
   ENDIF

!  ---update o.b. data at initial time
   IF (nofiles.GT.1) THEN
      CALL update_profobc_data(profsal,obcsal3d,noprofsd,&
                             & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                             & maxprofs,noprofs,1,nofiles,nosecsprof,io_salobc)
   ENDIF

!  ---update discharge data at initial time
   IF (disch_data) THEN
      CALL update_dischr_data(io_dissal,1,disdata,dissal,numdis,1,nosecsdis,&
                            & update)
   ENDIF

!  ---write at nest locations
   IF (iopt_nests.EQ.1) THEN
      CALL update_nest_data_prof(sal,1,novarsnst,instvars,io_salnst)
   ENDIF

   GOTO 1000

ENDIF

!
!2. Allocate
!-----------
!

ALLOCATE (source(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('source',3,(/ncloc,nrloc,nz/),kndrtype)

!
!3. Open boundary conditions
!---------------------------
!
!---update data
IF (nofiles.GT.1) THEN
   CALL update_profobc_data(profsal,obcsal3d,noprofsd,&
                          & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                          & maxprofs,noprofs,1,nofiles,nosecsprof,io_salobc)
ENDIF

!---update o.b. profiles
IF (iopt_obc_sal.EQ.1) THEN
   CALL open_boundary_conds_prof(sal,salobu,salobv,obcsal3d,obcsalatu,&
                               & obcsalatv,itypobu,itypobv,iprofobu,iprofobv,&
                               & noprofs,1)
ENDIF

!---update discharge data
IF (disch_data) THEN
   CALL update_dischr_data(io_dissal,1,disdata,dissal,numdis,1,nosecsdis,update)
ENDIF

!
!4. Surface fluxes
!-----------------
!

IF (iopt_sflux_precip.EQ.1.OR.iopt_sflux_precip.EQ.3) CALL salinity_flux

!
!5. Source term
!--------------
!

k_510: DO k=1,nz
   WHERE (maskatc_int) source(:,:,k) = 0.0
ENDDO k_510

!---add discharges
IF (disch_data) CALL scalar_discharge(dissal,source,1,1) 

!
!6. Solve transport equation
!---------------------------
!

CALL transport_at_C_3d(sal,source,vdifcoefscal,iopt_adv_scal,iopt_hdif_scal,&
                     & salobu,salobv,iprofrlx,1,0,1,1,ssalflux,zeros2d,iarr_sal)

!
!7. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = nhalo
   CALL exchange_mod(sal,(/1-nhalo,1-nhalo,1/),nhexch,iarr_sal,&
                   & corners=.FALSE.)
ENDIF

!
!8. Write at nest locations
!---------------------------
!

IF (iopt_nests.EQ.1) THEN
   CALL update_nest_data_prof(sal,1,novarsnst,instvars,io_salnst)
ENDIF

!
!9. Deallocate arrays
!--------------------
!
!---work space
DEALLOCATE (source)

!---at final step
1000 CONTINUE

IF (cold_start.OR.nt.EQ.nstep) THEN
   DEALLOCATE (salobu,salobv,iprofrlx)
   IF (iopt_obc_sal.EQ.1) THEN
      DEALLOCATE (itypobu,itypobv,iprofobu,iprofobv,noprofsd,indexprof,&
                & indexvar,obcsal3d)
   ENDIF
   IF (iopt_nests.EQ.1) DEALLOCATE (novarsnst,instvars)
   IF (nofiles.GT.1) DEALLOCATE (nosecsprof,profsal)
   IF (disch_data) DEALLOCATE (dissal,disdata)
ENDIF

CALL log_timer_out(npcc,itm_sal)


RETURN

END SUBROUTINE salinity_equation

!========================================================================

SUBROUTINE temperature_equation
!************************************************************************
!
! *temperature_equation* Solve temperature equation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Density_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! External calls - define_profobc_spec, define_surface_input_grid, heat_flux,
!                  heat_optics, long_wave_radiation, open_boundary_conds_prof,
!                  scalar_discharge, short_wave_radiation, transport_at_C_3d,
!                  update_dischr_data, update_nest_data_prof,
!                  update_profobc_data, update_surface_data
!
! Module calls - error_alloc, error_alloc_struc, exchange_mod,
!                hrelativecoords_init
!
!************************************************************************
!
USE currents
USE datatypes
USE density
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE modids
USE nestgrids
USE obconds
USE optics
USE physpars
USE relaxation
USE structures
USE switches
USE syspars
USE timepars
USE datatypes_init, ONLY: hrelativecoords_init
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: update
LOGICAL, SAVE :: disch_data
INTEGER :: k, npcc
INTEGER, SAVE :: ibcsur, maxprofs, nhtype, nofiles, noprofs, n1dat, n1grd, &
               & n2dat, n2grd
REAL :: xfac
INTEGER, DIMENSION(4) :: nhexch
LOGICAL, SAVE, DIMENSION(1) :: angle = .FALSE.
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: iprofrlx, itypobu, itypobv, &
                                          & noprofsd, novarsnst
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: indexprof, indexvar, instvars, &
                                            & iprofobu, iprofobv
INTEGER (KIND=kndilong), SAVE, DIMENSION(2) :: nosecsdis
INTEGER (KIND=kndilong), SAVE, ALLOCATABLE, DIMENSION(:,:) :: nosecsprof
INTEGER (KIND=kndilong), SAVE, DIMENSION(2) :: nosecssst
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: distmp
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: bcsur, disdata
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: obctemp3d, source, tempindat, &
                                           & tmpobu, tmpobv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: proftmp
TYPE (HRelativeCoords), SAVE, ALLOCATABLE, DIMENSION(:,:) :: sstgrid

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*surflux*      REAL    Surface flux                                [deg C m/s]
!*sstgrid*      DERIVED Model grid coordinates with respect to sst grid
!
!------------------------------------------------------------------------------
!

procname(pglev+1) = 'temperature_equation'
CALL log_timer_in(npcc)

!
!1. Initialise open boundary conditions on first call
!----------------------------------------------------
!

IF (nt.EQ.0) THEN

!
!1.1 Open boundary conditions
!----------------------------
!

   nofiles = maxdatafiles(io_tmpobc,1)
   disch_data = iopt_dischr.EQ.1.AND.modfiles(io_distmp,1,1)%defined

!  ---type of surface boundary conditions
   SELECT CASE (iopt_temp_sbc)
      CASE (1); ibcsur = 1
      CASE (2); ibcsur = 3
      CASE (3); ibcsur = 4
   END SELECT

!  ---allocate/initialise specifier arrays
   ALLOCATE (tmpobu(nobu,nz,1),STAT=errstat)
   CALL error_alloc('tmpobu',3,(/nobu,nz,1/),kndrtype)
   IF (nobu.GT.0) tmpobu = float_fill
   ALLOCATE (tmpobv(nobv,nz,1),STAT=errstat)
   CALL error_alloc('tmpobv',3,(/nobv,nz,1/),kndrtype)
   IF (nobv.GT.0) tmpobv = float_fill
   ALLOCATE (iprofrlx(norlxzones),STAT=errstat)
   CALL error_alloc('iprofrlx',1,(/norlxzones/),kndint)
   IF (norlxzones.GT.0) iprofrlx = 0
   IF (iopt_obc_temp.EQ.1) THEN
      ALLOCATE (itypobu(nobu),STAT=errstat)
      CALL error_alloc('itypobu',1,(/nobu/),kndint)
      ALLOCATE (itypobv(nobv),STAT=errstat)
      CALL error_alloc('itypobv',1,(/nobv/),kndint)
      ALLOCATE (iprofobu(nobu,1),STAT=errstat)
      CALL error_alloc('iprofobu',2,(/nobu,1/),kndint)
      ALLOCATE (iprofobv(nobv,1),STAT=errstat)
      CALL error_alloc('iprofobv',2,(/nobv,1/),kndint)
      ALLOCATE (noprofsd(2:nofiles),STAT=errstat)
      CALL error_alloc('noprofsd',1,(/nofiles-1/),kndint)
      ALLOCATE (indexprof(nobu+nobv,2:nofiles),STAT=errstat)
      CALL error_alloc('indexprof',2,(/nobu+nobv,nofiles-1/),kndint)
      ALLOCATE (indexvar(nobu+nobv,2:nofiles),STAT=errstat)
      CALL error_alloc('indexvar',2,(/nobu+nobv,nofiles-1/),kndint)
   ENDIF
   IF (iopt_nests.EQ.1) THEN
      ALLOCATE (novarsnst(nonestsets),STAT=errstat)
      CALL error_alloc('novarsnst',1,(/nonestsets/),kndint)
      novarsnst = 1
      ALLOCATE (instvars(1,nonestsets),STAT=errstat)
      CALL error_alloc('instvars',2,(/1,nonestsets/),kndint)
      instvars = 1
   ENDIF

!  ---define specifier arrays
   IF (iopt_obc_temp.EQ.1) THEN
      CALL define_profobc_spec(io_tmpobc,itypobu,itypobv,iprofobu,iprofobv,&
                             & iprofrlx,noprofsd,indexprof,indexvar,1,nofiles,&
                             & nobu,nobv)
   ENDIF

!  ---allocate o.b. data arrays
   noprofs = 0; maxprofs = 0
   IF (iopt_obc_temp.EQ.1) THEN
      IF (nofiles.GT.1) noprofs = MAX(MAXVAL(iprofobu),MAXVAL(iprofobv))
      ALLOCATE (obctemp3d(0:noprofs,nz,1),STAT=errstat)
      CALL error_alloc('obctemp3d',3,(/noprofs+1,nz,1/),kndrtype)
      obctemp3d = float_fill
   ENDIF
   IF (nofiles.GT.1) THEN
      maxprofs = MAXVAL(noprofsd)
      ALLOCATE (nosecsprof(2:nofiles,2),STAT=errstat)
      CALL error_alloc('nosecsprof',2,(/nofiles-1,2/),kndilong)
      ALLOCATE (proftmp(maxprofs,nz,2:nofiles,2),STAT=errstat)
      CALL error_alloc('proftmp',4,(/maxprofs,nz,nofiles-1,2/),kndrtype)
   ENDIF

!  ---allocate discharge data arrays
   IF (disch_data) THEN
      ALLOCATE (distmp(numdis),STAT=errstat)
      CALL error_alloc('distmp',1,(/numdis/),kndrtype)
      ALLOCATE (disdata(numdis,2),STAT=errstat)
      CALL error_alloc('disdata',2,(/numdis,2/),kndrtype)
   ENDIF

!  ---update o.b. data at initial time
   IF (nofiles.GT.1) THEN
      CALL update_profobc_data(proftmp,obctemp3d,noprofsd,&
                             & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                             & maxprofs,noprofs,1,nofiles,nosecsprof,io_tmpobc)
   ENDIF

!  ---update discharge data at initial time
   IF (disch_data) THEN
      CALL update_dischr_data(io_distmp,1,disdata,distmp,numdis,1,nosecsdis,&
                            & update)
   ENDIF

!  ---write at nest locations
   IF (iopt_nests.EQ.1) THEN
      CALL update_nest_data_prof(temp,1,novarsnst,instvars,io_tmpnst)
   ENDIF

!
!1.2 Surface grid
!----------------
!

   nhtype = surfacegrids(igrd_sst,1)%nhtype
   IF (iopt_temp_sbc.GT.1) THEN

!     ---grid parameters
      n1dat = surfacegrids(igrd_sst,1)%n1dat
      n2dat = surfacegrids(igrd_sst,1)%n2dat
      n1grd  = MERGE(ncloc,0,nhtype.GT.0.AND.nhtype.LT.4)
      n2grd  = MERGE(nrloc,0,nhtype.GT.0.AND.nhtype.LT.4)

!     ---allocate coordinate structure
      ALLOCATE (sstgrid(n1grd,n2grd),STAT=errstat)
      CALL error_alloc_struc('sstgrid',2,(/n1grd,n2grd/),'HRelativeCoords')
      CALL hrelativecoords_init(sstgrid,.FALSE.)

!     ---define surface grid
      IF (nhtype.GT.0.AND.nhtype.LT.4) THEN
         CALL define_surface_input_grid(igrd_meteo,1,sstgrid)
      ENDIF

   ENDIF

!
!1.3 Allocate data arrays
!------------------------
!

   ALLOCATE (bcsur(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bcsur',2,(/ncloc,nrloc/),kndrtype)
   bcsur = 0.0
   IF (iopt_temp_sbc.GT.1) THEN
      ALLOCATE (tempindat(n1dat,n2dat,2),STAT=errstat)
      CALL error_alloc('tempindat',3,(/n1dat,n2dat,2/),kndrtype)
   ENDIF

!
!1.4 Update surface data at initial time
!---------------------------------------
!
!  ---sst as surface boundary condition
   IF (iopt_temp_sbc.GT.1) THEN
      nosecssst(1) = 0
      CALL update_surface_data(io_sstsur,1,tempindat,bcsur,sstgrid,n1dat,&
                             & n2dat,1,n1grd,n2grd,nosecssst,angle)
      WHERE (ABS(bcsur-float_fill).LE.float_fill_eps)
         bcsur = temp(1:ncloc,1:nrloc,nz)
      END WHERE
      WHERE (maskatc_int)
         sst(1:ncloc,1:nrloc) = bcsur
      END WHERE
   ENDIF

!
!1.5 Heat fluxes
!---------------
!   
!  ---short wave radiation
   IF (iopt_meteo_data.EQ.1) CALL short_wave_radiation
!  ---solar radiance
   IF (iopt_temp_optic.EQ.1) CALL heat_optics
!  ---non-solar heat fluxes
   IF (iopt_temp_sbc.EQ.1.AND.iopt_meteo_data.EQ.1) THEN
      CALL heat_flux
   ENDIF

   GOTO 1000

ENDIF

!
!2. Allocate
!-----------
!

ALLOCATE (source(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('source',3,(/ncloc,nrloc,nz/),kndrtype)

!
!3. Open boundary conditions
!---------------------------
!
!---update data
IF (nofiles.GT.1) THEN
   CALL update_profobc_data(proftmp,obctemp3d,noprofsd,&
                          & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                          & maxprofs,noprofs,1,nofiles,nosecsprof,io_tmpobc)
ENDIF

!---update o.b. profiles
IF (iopt_obc_temp.EQ.1) THEN
   CALL open_boundary_conds_prof(temp,tmpobu,tmpobv,obctemp3d,obctmpatu,&
                               & obctmpatv,itypobu,itypobv,iprofobu,iprofobv,&
                               & noprofs,1)
ENDIF

!---update discharge data
IF (disch_data) THEN
   CALL update_dischr_data(io_distmp,1,disdata,distmp,numdis,1,nosecsdis,update)
ENDIF

!
!4. Source terms
!---------------
!
!---short wave radiative flux
xfac = density_ref*specheat
IF (iopt_meteo_data.EQ.1) CALL short_wave_radiation
!---solar radiance
IF (iopt_temp_optic.EQ.0) THEN
   source = 0.0
ELSE
   CALL heat_optics
   k_410: DO k=1,nz
      WHERE (maskatc_int)
         source(:,:,k) = (radiance(:,:,k+1)-radiance(:,:,k))&
                      & /(xfac*delzatc(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_410
ENDIF

!---discharges
IF (disch_data) CALL scalar_discharge(distmp,source,1,1)

!
!5. Surface boundary condition
!-----------------------------
!
!5.1 Neumann
!-----------
!

IF (iopt_temp_sbc.EQ.1) THEN
   bcsur = 0.0
   IF (iopt_meteo_data.EQ.1) CALL heat_flux
   IF (iopt_temp_optic.EQ.0) THEN
      WHERE (maskatc_int) bcsur = (qrad-qnonsol)/xfac
   ELSE
      WHERE (maskatc_int) bcsur = -qnonsol/xfac
   ENDIF

!
!5.2 Dirichlet
!-------------
!

ELSEIF (iopt_temp_sbc.GT.1) THEN
   CALL update_surface_data(io_sstsur,1,tempindat,bcsur,sstgrid,n1dat,n2dat,1,&
                          & n1grd,n2grd,nosecssst,angle)
   WHERE (ABS(bcsur-float_fill).LE.float_fill_eps)
      bcsur = temp(1:ncloc,1:nrloc,nz)
   END WHERE
   WHERE (maskatc_int)
      sst(1:ncloc,1:nrloc) = bcsur
   END WHERE
ENDIF

!
!6. Solve transport equation
!---------------------------
!

CALL transport_at_C_3d(temp,source,vdifcoefscal,iopt_adv_scal,iopt_hdif_scal,&
                     & tmpobu,tmpobv,iprofrlx,ibcsur,0,1,1,bcsur,zeros2d,&
                     & iarr_temp)

!
!7. Apply freezing point limit
!-----------------------------
!
!---imposed lower limit
IF (ABS(temp_min-float_fill).GT.float_fill_eps) THEN
   k_710: DO k=1,nz
      WHERE (maskatc_int)
         temp(1:ncloc,1:nrloc,k) = MAX(temp_min,temp(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_710

!---freezing point limit
ELSE
   k_720: DO k=1,nz
      WHERE (maskatc_int)
         temp(1:ncloc,1:nrloc,k) = MAX(-0.0575*sal(1:ncloc,1:nrloc,k),&
                                     & temp(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_720
ENDIF

!
!8. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = nhalo
   CALL exchange_mod(temp,(/1-nhalo,1-nhalo,1/),nhexch,iarr_temp,&
                   & corners=.FALSE.)
ENDIF

!
!9. Write at nest locations
!---------------------------
!

IF (iopt_nests.EQ.1) THEN
   CALL update_nest_data_prof(temp,1,novarsnst,instvars,io_tmpnst)
ENDIF

!
!10. Deallocate arrays
!---------------------
!
!---work space
DEALLOCATE (source)

!---at final step
1000 CONTINUE
IF (cold_start.OR.nt.EQ.nstep) THEN
   DEALLOCATE (tmpobu,tmpobv,iprofrlx)
   IF (iopt_obc_temp.EQ.1) THEN
      DEALLOCATE (itypobu,itypobv,iprofobu,iprofobv,noprofsd,indexprof,&
                & indexvar,obctemp3d)
   ENDIF
   IF (iopt_nests.EQ.1) DEALLOCATE (novarsnst,instvars)
   IF (nofiles.GT.1) DEALLOCATE (nosecsprof,proftmp)
   DEALLOCATE (bcsur)
   IF (iopt_temp_sbc.GT.1) DEALLOCATE (sstgrid,tempindat)
   IF (disch_data) DEALLOCATE (distmp,disdata)
ENDIF

CALL log_timer_out(npcc,itm_temp)


RETURN

END SUBROUTINE temperature_equation

!========================================================================

SUBROUTINE heat_optics
!************************************************************************
!
! *heat_optics* Optical module
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Density_Equations.f90  V2.11.2
!
! Description -  evaluates irradiance
!
! Reference -
!
! Calling program - temperature_equation
!
! External calls -
!
! Module calls - error_alloc
!
!************************************************************************
!
USE density
USE depths
USE grid
USE gridpars
USE iopars
USE optics
USE physpars
USE syspars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: k
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: rad1, rad2

!
! Name             Type    Purpose
!------------------------------------------------------------------------------
!*optattcoef1_tot* REAL Attenuation coefficient for short-wave radiation
!                       including contribution from salinity              [1/m]
!*rad1*            REAL Infrared part of irradiance                     [W/m^2]
!*rad2*            REAL Short-wave part of irradiance                   [W/m^2]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'heat_optics'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (rad1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('rad1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (rad2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('rad2',2,(/ncloc,nrloc/),kndrtype)

!
!2. Irradiance
!-------------
!
!---surface/bottom
WHERE (maskatc_int)
   rad1 = opt_frac*qrad
   rad2 = (1.0-opt_frac)*qrad
   radiance(:,:,nz+1) = qrad
   radiance(:,:,1) = 0.0
END WHERE

!---interior cells
k_210: DO k=nz,2,-1
   WHERE (maskatc_int)
      rad1 = rad1*EXP(-optattcoef1_cst*delzatc(1:ncloc,1:nrloc,k))
      rad2 = rad2*EXP(-optattcoef2*delzatc(1:ncloc,1:nrloc,k))
   END WHERE
   WHERE (maskatc_int) 
      radiance(:,:,k) = rad1 + rad2
   END WHERE
ENDDO k_210

!
!3. Deallocate arrays
!--------------------
!

DEALLOCATE (rad1,rad2)

CALL log_timer_out()


RETURN

END SUBROUTINE heat_optics

!========================================================================

SUBROUTINE equation_of_state
!************************************************************************
!
! *equation_of_state* Density and expansion coefficients
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Density_Equations.f90  V2.11.2
!
! Description - 
!
! Reference - McDougall T.J., Jackett D.R., Wright D.G. and Feistel R., 2003.
!      Accurate and computationally efficient algorithms for potential
!      temperature and density of seawater. Journal of Atmospheric and Oceanic
!      Technology, 20, 730-741.
!
! Calling program - coherens_main, initialise_model
!
! External calls - convect_adjust, equation_of_state_sed
!
! Module calls - error_alloc, mult_index
!
!************************************************************************
!
USE density
USE depths
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
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: salflag, sedimflag, tempflag
INTEGER :: k, npcc
REAL :: dp1x, dp2x, p1x, p2x, slx, s15x, tpx, tp2x
REAL, PARAMETER, DIMENSION(24) :: a = (/&
    & 999.843699,      7.3521284,      -5.45928211E-02,  3.98476704E-04, &
    & 2.96938239,     -7.23268813E-03,  2.12382341E-03,  1.04004591E-02, &
    & 1.03970529E-07,  5.1876188E-06,  -3.24041825E-08, -1.2386936E-11, &
    & 7.28606739E-03, -4.60835542E-05,  3.68390573E-07,  1.80809186E-10, &
    & 2.14691708E-03, -9.27062484E-06, -1.78343643E-10,  4.76534122E-06, &
    & 1.63410736E-09,  5.30848875E-06, -3.03175128E-16, -1.27934137E-17 /)
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, dp1, dp2, gacc, press, &
                                         & press2, p1, p2, sl, s15, tp, tp2


!
!1. Return if no update
!----------------------
!

tempflag = iopt_temp.EQ.2.OR.nt.EQ.0
salflag = iopt_sal.EQ.2.OR.nt.EQ.0
sedimflag = iopt_sed.GT.0
IF (.NOT.tempflag.AND.(.NOT.salflag).AND.(.NOT.sedimflag)) RETURN

procname(pglev+1) = 'equation_of_state'
CALL log_timer_in(npcc)

!
!2. Allocate arrays
!------------------
!


ALLOCATE (maskatc(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc+2,nrloc+2/),kndlog)

IF (iopt_dens.EQ.1) THEN
   ALLOCATE (array2dc(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('array2dc',2,(/ncloc+2,nrloc+2/),kndrtype)
ENDIF
   
IF (iopt_dens.GE.2) THEN
   ALLOCATE (dp1(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('dp1',2,(/ncloc+2,nrloc+2/),kndrtype)
   ALLOCATE (dp2(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('dp2',2,(/ncloc+2,nrloc+2/),kndrtype)
   ALLOCATE (p1(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('p1',2,(/ncloc+2,nrloc+2/),kndrtype)
   ALLOCATE (p2(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('p2',2,(/ncloc+2,nrloc+2/),kndrtype)
   ALLOCATE (sl(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('sl',2,(/ncloc+2,nrloc+2/),kndrtype)
   ALLOCATE (s15(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('s15',2,(/ncloc+2,nrloc+2/),kndrtype)
   ALLOCATE (tp(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('tp',2,(/ncloc+2,nrloc+2/),kndrtype)
   ALLOCATE (tp2(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('tp2',2,(/ncloc+2,nrloc+2/),kndrtype)
ENDIF
IF (iopt_dens.EQ.3) THEN
   ALLOCATE (gacc(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('gacc',2,(/ncloc+2,nrloc+2/),kndrtype)
   ALLOCATE (press(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('press',2,(/ncloc+2,nrloc+2/),kndrtype)
   ALLOCATE (press2(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('press2',2,(/ncloc+2,nrloc+2/),kndrtype)
ENDIF

!
!3. Initialise parameters and arrays
!-----------------------------------
!

maskatc = nodeatc(0:ncloc+1,0:nrloc+1).GT.0

!
!4. Initialise on first call
!---------------------------
!

IF (nt.EQ.0) THEN

!
!4.1 Reference values
!--------------------
!
!  ---density
   tpx = temp_ref; slx = sal_ref; tp2x = tpx*tpx; s15x = slx**1.5
   p1x = a(1)+tpx*(a(2)+tpx*(a(3)+tpx*a(4)))+slx*(a(5)+a(6)*tpx+a(7)*slx)
   p2x = 1.0+tpx*(a(13)+tpx*(a(14)+tpx*(a(15)+tpx*a(16))))&
       & +slx*(a(17)+tpx*(a(18)+tp2x*a(19)))+s15x*(a(20)+a(21)*tp2x)
   density_ref = p1x/p2x
   
!  ---thermal expansion coefficient
   dp1x = a(2)+tpx*(2.0*a(3)+3.0*a(4)*tpx)+a(6)*slx
   dp2x = a(13)+tpx*(2.0*a(14)+tpx*(3.0*a(15)+4.0*a(16)*tpx))&
        & +slx*(a(18)+3.0*a(19)*tp2x)+2.0*a(21)*s15x*tpx
   beta_temp_ref = dp2x/p2x - dp1x/p1x

!  ---salinity expansion coefficient
   dp1x = a(5)+a(6)*tpx+2.0*a(7)*slx
   dp2x = a(17)+tpx*(a(18)+a(19)*tp2x)+1.5*SQRT(slx)*(a(20)+a(21)*tp2x)
   beta_sal_ref = dp1x/p1x - dp2x/p2x

!
!4.2 Uniform density
!-------------------
!

   IF (iopt_dens.EQ.0) THEN
      k_420: DO k=1,nz
         WHERE (maskatc)
            dens(:,:,k) = density_ref
         END WHERE
      ENDDO k_420
   ENDIF

!
!4.3 Uniform expansion coefficients
!----------------------------------
!

   IF (iopt_dens.EQ.1) THEN
      k_430: DO k=1,nz
         WHERE (maskatc)
            beta_temp(:,:,k) = beta_temp_ref
            beta_sal(:,:,k) = beta_sal_ref
         END WHERE
      ENDDO k_430
   ENDIF

ENDIF

!
!5. Uniform density
!------------------
!

IF (sedimflag.AND.iopt_dens.EQ.0) THEN
   k_510: DO k=1,nz
      WHERE (maskatc)
         dens(:,:,k) = density_ref
      END WHERE
   ENDDO k_510
ENDIF

!
!6. Linear equation of state
!---------------------------
!

IF (iopt_dens.EQ.1) THEN
   k_610: DO k=1,nz
      WHERE (maskatc)
         dens(:,:,k) = 0.0
      END WHERE
      IF (iopt_temp.GT.0) THEN
         WHERE (maskatc)
            dens(:,:,k) = beta_temp_ref*(temp_ref-temp(0:ncloc+1,0:nrloc+1,k))
         END WHERE
      ENDIF
      IF (iopt_sal.GT.0) THEN
         WHERE (maskatc)
            dens(:,:,k) = dens(:,:,k) - &
                        & beta_sal_ref*(sal_ref-sal(0:ncloc+1,0:nrloc+1,k))
         END WHERE
      ENDIF
      WHERE (maskatc)
         dens(:,:,k) = density_ref*(1.0+dens(:,:,k))
      END WHERE
   ENDDO k_610

!
!7. General equation of state without pressure
!---------------------------------------------
!

ELSEIF (iopt_dens.EQ.2) THEN

   k_710: DO k=1,nz

!     ---density
      WHERE (maskatc)
         tp = temp(0:ncloc+1,0:nrloc+1,k)
         sl = sal(0:ncloc+1,0:nrloc+1,k)
         tp2 = tp*tp; s15 = sl**1.5
         p1 = a(1)+tp*(a(2)+tp*(a(3)+tp*a(4)))+sl*(a(5)+a(6)*tp+a(7)*sl)
         p2 = 1.0+tp*(a(13)+tp*(a(14)+tp*(a(15)+tp*a(16))))&
            & +sl*(a(17)+tp*(a(18)+tp2*a(19)))+s15*(a(20)+a(21)*tp2)
         dens(:,:,k) = p1/p2
      END WHERE

!     ---thermal expansion coefficient
      IF (tempflag) THEN
         WHERE (maskatc)
            dp1 = a(2)+tp*(2.0*a(3)+3.0*a(4)*tp)+a(6)*sl
            dp2 = a(13)+tp*(2.0*a(14)+tp*(3.0*a(15)+4.0*a(16)*tp))&
               & +sl*(a(18)+3.0*a(19)*tp2)+2.0*a(21)*s15*tp
            beta_temp(:,:,k) = dp2/p2 - dp1/p1
         END WHERE
      ENDIF

!     ---salinity expansion coefficient
      IF (salflag) THEN
         WHERE (maskatc)
            dp1 = a(5)+a(6)*tp+2.0*a(7)*sl
            dp2 = a(17)+tp*(a(18)+a(19)*tp2)+1.5*SQRT(sl)*(a(20)+a(21)*tp2)
            beta_sal(:,:,k) = dp1/p1 - dp2/p2
         END WHERE
      ENDIF

   ENDDO k_710

!
!8. General equation of state with pressure
!------------------------------------------
!

ELSEIF (iopt_dens.EQ.3) THEN

!  ---surface pressure
   WHERE (maskatc)
      gacc = 1.0E-04*gaccatc
      tp = temp(0:ncloc+1,0:nrloc+1,nz)
      sl = sal(0:ncloc+1,0:nrloc+1,nz)
      tp2 = tp*tp; s15 = sl**1.5
      p1 = a(1)+tp*(a(2)+tp*(a(3)+tp*a(4)))+sl*(a(5)+a(6)*tp+a(7)*sl)
      p2 = 1.0+tp*(a(13)+tp*(a(14)+tp*(a(15)+tp*a(16))))&
         & +sl*(a(17)+tp*(a(18)+tp2*a(19)))+s15*(a(20)+a(21)*tp2)
      dens(:,:,nz) = p1/p2
      press = 0.5*gacc*delzatc(:,:,nz)*dens(:,:,nz)
      press2 = press*press
   END WHERE

   k_810: DO k=nz,1,-1

!     ---pressure
      IF (k.LT.nz) THEN
         WHERE (maskatc)
            tp = temp(0:ncloc+1,0:nrloc+1,k)
            sl = sal(0:ncloc+1,0:nrloc+1,k)
            tp2 = tp*tp; s15 = sl**1.5
            press = press + gacc*dens(:,:,k+1)*delzatw(:,:,k+1)
            press2 = press*press
         END WHERE
      ENDIF

!     ---density
      WHERE (maskatc)
         p1 = a(1)+tp*(a(2)+tp*(a(3)+tp*a(4)))+sl*(a(5)+a(6)*tp+a(7)*sl)&
            & +press*(a(8)+a(9)*tp2+a(10)*sl)+press2*(a(11)+tp2*a(12))
         p2 = 1.0+tp*(a(13)+tp*(a(14)+tp*(a(15)+tp*a(16))))&
            & +sl*(a(17)+tp*(a(18)+a(19)*tp2))+s15*(a(20)+a(21)*tp2)&
            & +press*a(22)+press2*tp*(tp2*a(23)+press*a(24))  
         dens(:,:,k) = p1/p2
      END WHERE

!     ---thermal expansion coefficient
      IF (tempflag) THEN
         WHERE (maskatc)
            dp1 = a(2)+tp*(2.0*a(3)+3.0*a(4)*tp)+a(6)*sl&
               & +2.0*press*tp*(a(9)+press*a(12))
            dp2 = a(13)+tp*(2.0*a(14)+tp*(3.0*a(15)+4.0*a(16)*tp))&
               & +sl*(a(18)+3.0*a(19)*tp2)+2.0*a(21)*s15*tp&
               & +press2*(3.0*a(23)*tp2+a(24)*press) 
            beta_temp(:,:,k) = dp2/p2 - dp1/p1
         END WHERE
      ENDIF

!     ---salinity expansion coefficient
      IF (salflag) THEN
         WHERE (maskatc)
            dp1 = a(5)+a(6)*tp+2.0*a(7)*sl+a(10)*press
            dp2 = a(17)+tp*(a(18)+tp2*a(19))+1.5*SQRT(sl)*(a(20)+tp2*a(21))
            beta_sal(:,:,k) = dp1/p1 - dp2/p2
         END WHERE
      ENDIF

   ENDDO k_810

ENDIF

!
!9. Contribution of sediments
!----------------------------
!

IF (sedimflag) CALL equation_of_state_sed

!
!10. Apply convective adjustment
!------------------------------
!

IF (iopt_dens_convect.EQ.1.AND.mult_index(nt,iccvt,matchzero=.TRUE.)) THEN
   CALL convect_adjust
ENDIF

!
!11. Deallocate arrays
!---------------------
!

DEALLOCATE (maskatc)
IF (iopt_dens.EQ.1) DEALLOCATE (array2dc)
IF (iopt_dens.GE.2) DEALLOCATE (dp1,dp2,p1,p2,sl,s15,tp,tp2)
IF (iopt_dens.EQ.3) DEALLOCATE (gacc,press,press2)

CALL log_timer_out(npcc,itm_dens)


RETURN

END SUBROUTINE equation_of_state

!========================================================================

SUBROUTINE buoyancy_frequency
!************************************************************************
!
! *buoyancy_frequency* Squared buoyancy frequency
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Density_Equations.f90  V2.11.2
!
! Description - defined as vertical gradient of buoyancy
!             - uses interpolation of neighbouring cells
!
! Calling program - vertical_diff_coefs
!
! External calls - buoyancy_frequency_sed
!
! Module calls - Carr_at_W, error_alloc
!
!************************************************************************
!
USE density
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE switches
USE syspars
USE timepars
USE turbulence
USE array_interp, ONLY: Carr_at_W
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
LOGICAL :: salflag, sedimflag, tempflag
INTEGER :: k, npcc
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, wsum
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: beta, bgrad

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*wsum*     REAL    Sum of weight factors
!*bgrad*    REAL    Non-averaged buoyancy gradient                   [1/s^2]
!*beta*     REAL    Expansion coefficient at W-nodes
!
!************************************************************************
!
!1. Return if no update
!----------------------
!

tempflag = iopt_temp.EQ.2.OR.(iopt_temp.EQ.1.AND.nt.EQ.0)
salflag = iopt_sal.EQ.2.OR.(iopt_sal.EQ.1.AND.nt.EQ.0)
sedimflag = iopt_sed.GT.0
IF (.NOT.tempflag.AND.(.NOT.salflag).AND.(.NOT.sedimflag)) RETURN

procname(pglev+1) = 'buoyancy_frequency'
CALL log_timer_in(npcc)

!
!2. Allocate arrays
!------------------
!

ALLOCATE (maskatc(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('maskatc',2,(/ncloc+2,nrloc+2/),kndlog)
ALLOCATE (array2dc(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc+2,nrloc+2/),kndrtype)
ALLOCATE (beta(0:ncloc+1,0:nrloc+1,2:nz),STAT=errstat)
CALL error_alloc('beta',3,(/ncloc+2,nrloc+2,nz-1/),kndrtype)
ALLOCATE (bgrad(0:ncloc+1,0:nrloc+1,2:nz),STAT=errstat)
CALL error_alloc('bgrad',3,(/ncloc+2,nrloc+2,nz-1/),kndrtype)
ALLOCATE (wsum(ncloc,nrloc),STAT=errstat)
CALL error_alloc('wsum',2,(/ncloc,nrloc/),kndrtype)

!
!3. Initialise
!-------------
!
!---mask array
maskatc = nodeatc(0:ncloc+1,0:nrloc+1).GT.0

!---sum of weight factors
WHERE (maskatc)
   array2dc = 1.0
ELSEWHERE
   array2dc = 0.0
END WHERE
wsum = 2.0*array2dc(1:ncloc,1:nrloc) + array2dc(0:ncloc-1,1:nrloc) +&
     & array2dc(2:ncloc+1,1:nrloc) + array2dc(1:ncloc,0:nrloc-1) +&
     & array2dc(1:ncloc,2:nrloc+1)

!---initialise buoyancy gradient
bgrad = 0.0

!
!4. Buoyancy frequency
!---------------------
!
!4.1 Temperature
!---------------
!

IF (tempflag) THEN
   CALL Carr_at_W(beta_temp,beta,(/0,0,1/),(/ncloc+1,nrloc+1,nz/),1,&
                & iarr_beta_temp,.TRUE.)

   k_410: DO k=2,nz
      WHERE (maskatc)
         bgrad(:,:,k) = gaccatc*beta(:,:,k)*&
               & (temp(0:ncloc+1,0:nrloc+1,k)-temp(0:ncloc+1,0:nrloc+1,k-1))&
               & /delzatw(:,:,k)
      END WHERE
   ENDDO k_410
ENDIF

!
!4.2 Salinity
!------------
!

IF (salflag) THEN
   CALL Carr_at_W(beta_sal,beta,(/0,0,1/),(/ncloc+1,nrloc+1,nz/),1,&
                & iarr_beta_sal,.TRUE.)

   k_420: DO k=2,nz
      WHERE (maskatc)
         bgrad(:,:,k) = bgrad(:,:,k) - gaccatc*beta(:,:,k)*&
               & (sal(0:ncloc+1,0:nrloc+1,k)-sal(0:ncloc+1,0:nrloc+1,k-1))&
               & /delzatw(:,:,k)
      END WHERE
   ENDDO k_420
ENDIF

!
!4.3 Sediment
!------------
!

IF (sedimflag) CALL buoyancy_frequency_sed(bgrad)

!
!4.4 Horizontal averages
!-----------------------
!

k_430: DO k=2,nz
   WHERE (maskatc_int)
      buofreq2(:,:,k) = (2.0*bgrad(1:ncloc,1:nrloc,k) + &
                  & bgrad(0:ncloc-1,1:nrloc,k) + bgrad(2:ncloc+1,1:nrloc,k) + &
                  & bgrad(1:ncloc,0:nrloc-1,k) + bgrad(1:ncloc,2:nrloc+1,k))&
                  & /wsum
   END WHERE

ENDDO k_430

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatc,array2dc,beta,bgrad,wsum)

CALL log_timer_out(npcc,itm_dens)


RETURN

END SUBROUTINE buoyancy_frequency

!========================================================================

SUBROUTINE baroclinic_gradient
!************************************************************************
!
! *baroclinic_gradient* Baroclinic pressure gradient
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Density_Equations.f90  V2.11.2
!
! Description - 
!
! Reference - Shchepetkin A.F. and McWilliams J.C., 2003. A method for
!    computing horizontal pressure-gradient force in an oceanic model with a
!    nonaligned vertical coordinate. Journal of Geophysical Research, 108,
!    3090, doi:10.1029/2001JC001407.
!
! Calling program - coherens_main
!
! External calls - baroclinic_gradient_sed_cubic,
!                  baroclinic_gradient_sed_sigma, baroclinic_gradient_sed_z,
!                  hcube_deriv, hcube_fluxes
!
! Module calls - error_alloc, exchange_mod, Carr_at_U, Carr_at_UW, Carr_at_V,
!                Carr_at_VW, Carr_at_W, Zcoord_arr, Zcoord_var
!
!************************************************************************
!
USE currents
USE density
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE switches
USE syspars
USE timepars
USE array_interp, ONLY: Carr_at_U, Carr_at_UW, Carr_at_V, Carr_at_VW, Carr_at_W
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: Zcoord_arr, Zcoord_var
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: salflag, sedimflag, tempflag
INTEGER :: i, ii, j, jj, k, kc, l, npcc
REAL :: sigk, z
INTEGER, DIMENSION(4) :: nhexch
REAL, DIMENSION(nz+1,2) :: tsint
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: masksuratu, masksuratv
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu, maskatv
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: kint
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d, array2du, array2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: beta, dhts, dxyts, dzts, dzx, &
                           & dzy, dzz, ftsatvel, ftsatw, tsatvel, tsatw, zcoord
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: sigint


!
!1. Return if no update
!----------------------
!

tempflag = iopt_temp.EQ.2.OR.(iopt_temp.EQ.1.AND.nt.EQ.0)
salflag = iopt_sal.EQ.2.OR.(iopt_sal.EQ.1.AND.nt.EQ.0)
sedimflag = iopt_sed.GT.0
IF (.NOT.tempflag.AND.(.NOT.salflag).AND.(.NOT.sedimflag)) RETURN

procname(pglev+1) = 'baroclinic_gradient'
CALL log_timer_in(npcc)

!
!2. Allocate arrays
!------------------
!

ALLOCATE (maskatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc,nz/),kndlog)

ALLOCATE (masksuratu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('masksuratu',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (masksuratv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('masksuratv',2,(/ncloc,nrloc/),kndlog)

ALLOCATE (array2d(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2du(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc/),kndrtype)

ALLOCATE (beta(ncloc,nrloc,2:nz+1),STAT=errstat)
CALL error_alloc('beta',3,(/ncloc,nrloc,nz/),kndrtype)

IF (iopt_dens_grad.EQ.1.OR.iopt_dens_grad.EQ.3) THEN
   ALLOCATE (zcoord(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
   CALL error_alloc('zcoord',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
ENDIF

IF (iopt_dens_grad.EQ.1) THEN
   ALLOCATE (tsatw(0:ncloc,0:nrloc,2:nz),STAT=errstat)
   CALL error_alloc('tsatw',3,(/ncloc+1,nrloc+1,nz-1/),kndrtype)
   ALLOCATE (tsatvel(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('tsatvel',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

IF (iopt_dens_grad.EQ.2) THEN
   ALLOCATE (kint(ncloc,nrloc,2:nz+1,2),STAT=errstat)
   CALL error_alloc('kint',4,(/ncloc,nrloc,nz,2/),kndint)
   ALLOCATE (sigint(ncloc,nrloc,2:nz+1,2),STAT=errstat)
   CALL error_alloc('sigint',4,(/ncloc,nrloc,nz,2/),kndrtype)
   ALLOCATE (dhts(ncloc,nrloc,2:nz+1),STAT=errstat)
   CALL error_alloc('dhts',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

IF (iopt_dens_grad.EQ.3) THEN
   ALLOCATE (dzz(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
   CALL error_alloc('dzz',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
   ALLOCATE (dzx(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
   CALL error_alloc('dzx',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
   ALLOCATE (dzy(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
   CALL error_alloc('dzy',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
   ALLOCATE (dzts(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
   CALL error_alloc('dzts',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
   ALLOCATE (dxyts(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
   CALL error_alloc('dxyts',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
   ALLOCATE (ftsatw(0:ncloc,0:nrloc,nz),STAT=errstat)
   CALL error_alloc('ftsatw',3,(/ncloc+1,nrloc+1,nz/),kndrtype)
   ALLOCATE (ftsatvel(0:ncloc,0:nrloc,nz),STAT=errstat)
   CALL error_alloc('ftsatw',3,(/ncloc+1,nrloc+1,nz/),kndrtype)
ENDIF

!
!3. Initialise parameters and arrays
!-----------------------------------
!

p3dbcgradatu = 0.0; p3dbcgradatv = 0.0

maskatu = nodeatu(1:ncloc,1:nrloc,:).EQ.2
maskatv = nodeatv(1:ncloc,1:nrloc,:).EQ.2
masksuratu = nodeatu(1:ncloc,1:nrloc,nz).EQ.2
masksuratv = nodeatv(1:ncloc,1:nrloc,nz).EQ.2

WHERE (masksuratu)
   array2du = gaccatu(1:ncloc,:)/delxatu(1:ncloc,1:nrloc)
ELSEWHERE
   array2du = 0.0
END WHERE
WHERE (masksuratv)
   array2dv = gaccatv(:,1:nrloc)/delyatv(1:ncloc,1:nrloc)
ELSEWHERE
   array2dv = 0.0
END WHERE

!
!4. Sigma second order
!---------------------
!

IF (iopt_dens_grad.EQ.1) THEN

!
!4.1 Z at W-nodes
!----------------
!

   CALL Zcoord_arr(zcoord(0:ncloc,0:nrloc,:),(/0,0,1/),(/ncloc,nrloc,nz/),&
                & 'W  ',.FALSE.)

!
!4.2 Temperature
!---------------
!

   IF (tempflag) THEN

      CALL Carr_at_W(temp(0:ncloc,0:nrloc,:),tsatw,(/0,0,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_temp,.TRUE.)

!
!4.2.1 X-direction
!-----------------
!
!     ---interpolate arrays
      CALL Carr_at_U(beta_temp(0:ncloc,1:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)
      CALL Carr_at_UW(beta_temp(0:ncloc,1:nrloc,:),beta(:,:,2:nz),2,&
                   & (/0,1,1/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)
      CALL Carr_at_U(temp(0:ncloc,1:nrloc,:),tsatvel,1,2,(/0,1,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_temp,.TRUE.)

!     ---surface
      WHERE (masksuratu)
         array2d = 0.5*array2du*delzatu(1:ncloc,:,nz)*beta(:,:,nz+1)*&
                 & (temp(1:ncloc,1:nrloc,nz)-temp(0:ncloc-1,1:nrloc,nz))
         p3dbcgradatu(:,:,nz) = array2d
      END WHERE

!     ---interior points
      k_421: DO k=nz,2,-1
         WHERE (maskatu(:,:,k))
            array2d = array2d + array2du*beta(:,:,k)*&
                & ((tsatw(1:ncloc,1:nrloc,k)-tsatw(0:ncloc-1,1:nrloc,k))*&
                & delzatuw(1:ncloc,:,k)-&
                & (tsatvel(:,:,k)-tsatvel(:,:,k-1))*&
                & (zcoord(1:ncloc,1:nrloc,k)-zcoord(0:ncloc-1,1:nrloc,k)))
            p3dbcgradatu(:,:,k-1) = array2d
         END WHERE
      ENDDO k_421

!
!4.2.2 Y-direction
!-----------------
!
!     ---interpolate arrays
      CALL Carr_at_V(beta_temp(1:ncloc,0:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)
      CALL Carr_at_VW(beta_temp(1:ncloc,0:nrloc,:),beta(:,:,2:nz),2,&
                   & (/1,0,1/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)
      CALL Carr_at_V(temp(1:ncloc,0:nrloc,:),tsatvel,1,2,(/1,0,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_temp,.TRUE.)

!     ---surface
      WHERE (masksuratv)
         array2d = 0.5*array2dv*delzatv(:,1:nrloc,nz)*beta(:,:,nz+1)*&
                 & (temp(1:ncloc,1:nrloc,nz)-temp(1:ncloc,0:nrloc-1,nz))
         p3dbcgradatv(:,:,nz) = array2d
      END WHERE

!     ---interior points
      k_422: DO k=nz,2,-1
         WHERE (maskatv(:,:,k))
            array2d = array2d + array2dv*beta(:,:,k)*&
                & ((tsatw(1:ncloc,1:nrloc,k)-tsatw(1:ncloc,0:nrloc-1,k))*&
                & delzatvw(:,1:nrloc,k)-&
                & (tsatvel(:,:,k)-tsatvel(:,:,k-1))*&
                & (zcoord(1:ncloc,1:nrloc,k)-zcoord(1:ncloc,0:nrloc-1,k)))
            p3dbcgradatv(:,:,k-1) = array2d
         END WHERE
      ENDDO k_422

   ENDIF

!
!4.3 Salinity
!------------
!

   IF (salflag) THEN

      CALL Carr_at_W(sal(0:ncloc,0:nrloc,:),tsatw,(/0,0,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_sal,.TRUE.)

!
!4.3.1 X-direction
!-----------------
!
!     ---interpolate arrays
      CALL Carr_at_U(beta_sal(0:ncloc,1:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)
      CALL Carr_at_UW(beta_sal(0:ncloc,1:nrloc,:),beta(:,:,2:nz),2,&
                   & (/0,1,1/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)
      CALL Carr_at_U(sal(0:ncloc,1:nrloc,:),tsatvel,1,2,(/0,1,1/),&
                      & (/ncloc,nrloc,nz/),1,iarr_sal,.TRUE.)

!     ---surface
      WHERE (masksuratu)
         array2d = 0.5*array2du*delzatu(1:ncloc,:,nz)*beta(:,:,nz+1)*&
                 & (sal(1:ncloc,1:nrloc,nz)-sal(0:ncloc-1,1:nrloc,nz))
         p3dbcgradatu(:,:,nz) = p3dbcgradatu(:,:,nz) - array2d
      END WHERE

!     ---interior points
      k_431: DO k=nz,2,-1
         WHERE (maskatu(:,:,k))
            array2d = array2d + array2du*beta(:,:,k)*&
                & ((tsatw(1:ncloc,1:nrloc,k)-tsatw(0:ncloc-1,1:nrloc,k))*&
                & delzatuw(1:ncloc,:,k)-&
                & (tsatvel(:,:,k)-tsatvel(:,:,k-1))*&
                & (zcoord(1:ncloc,1:nrloc,k)-zcoord(0:ncloc-1,1:nrloc,k)))
            p3dbcgradatu(:,:,k-1) = p3dbcgradatu(:,:,k-1) - array2d
         END WHERE
      ENDDO k_431

!
!4.3.2 Y-direction
!-----------------
!
!     ---interpolate arrays
      CALL Carr_at_V(beta_sal(1:ncloc,0:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)
      CALL Carr_at_VW(beta_sal(1:ncloc,0:nrloc,:),beta(:,:,2:nz),2,&
                   & (/1,0,1/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)
      CALL Carr_at_V(sal(1:ncloc,0:nrloc,:),tsatvel,1,2,(/1,0,1/),&
                  & (/ncloc,nrloc,nz/),1,iarr_sal,.TRUE.)

!     ---surface
      WHERE (masksuratv)
         array2d = 0.5*array2dv*delzatv(:,1:nrloc,nz)*beta(:,:,nz+1)*&
                 & (sal(1:ncloc,1:nrloc,nz)-sal(1:ncloc,0:nrloc-1,nz))
         p3dbcgradatv(:,:,nz) = p3dbcgradatv(:,:,nz) - array2d
      END WHERE

!     ---interior points
      k_432: DO k=nz,2,-1
         WHERE (maskatv(:,:,k))
            array2d = array2d + array2dv*beta(:,:,k)*&
                & ((tsatw(1:ncloc,1:nrloc,k)-tsatw(1:ncloc,0:nrloc-1,k))*&
                & delzatvw(:,1:nrloc,k)-&
                & (tsatvel(:,:,k)-tsatvel(:,:,k-1))*&
                & (zcoord(1:ncloc,1:nrloc,k)-zcoord(1:ncloc,0:nrloc-1,k)))
            p3dbcgradatv(:,:,k-1) = p3dbcgradatv(:,:,k-1) - array2d
         END WHERE
      ENDDO k_432

   ENDIF

!
!4.4. Sediment
!-------------
!

   IF (sedimflag) THEN
      CALL baroclinic_gradient_sed_sigma(zcoord,'X')
      CALL baroclinic_gradient_sed_sigma(zcoord,'Y')
   ENDIF

!
!5. Z-level
!----------
!

ELSEIF (iopt_dens_grad.EQ.2) THEN

!
!5.1 X-direction
!----------------
!
!5.1.1 Vertical levels for interpolation
!---------------------------------------
!

   j_511: DO j=1,nrloc
   i_511: DO i=1,ncloc
      IF (node2du(i,j).EQ.2) THEN
         ii_5111: DO ii=i-1,i
            l = ii-i+2
            kc = 1
            k_51111: DO k=2,nz+1
               IF (k.EQ.nz+1) THEN
                  z = Zcoord_var(i,j,nz,'U  ',.FALSE.)
               ELSE
                  z = Zcoord_var(i,j,k,'UW ',.FALSE.)
               ENDIF
               sigk = (z+depmeanatc(ii,j))/deptotatc(ii,j)
               IF (sigk.LE.gscoordatc(ii,j,1)) THEN
                  kc = 0
                  sigint(i,j,k,l) = 0.0
               ELSEIF (sigk.GE.gscoordatc(ii,j,nz)) THEN
                  kc = nz
                  sigint(i,j,k,l) = 0.0
               ELSE
                  DO WHILE (sigk.GE.gscoordatc(ii,j,kc+1))
                     kc = kc + 1
                  END DO
                  sigint(i,j,k,l) = deptotatc(ii,j)*(sigk-gscoordatc(ii,j,kc))&
                                  & /delzatw(ii,j,kc+1)
               ENDIF
               kint(i,j,k,l) = kc
            ENDDO k_51111
         ENDDO ii_5111
        ENDIF
   ENDDO i_511
ENDDO j_511

!
!5.1.2 Temperature
!-----------------
!

   IF (tempflag) THEN

!     ---horizontal temperature difference
      j_5121: DO j=1,nrloc
      i_5121: DO i=1,ncloc
      k_5121: DO k=2,nz+1
         IF ((k.EQ.(nz+1).AND.masksuratu(i,j)).OR.nodeatuw(i,j,k).EQ.2) THEN
            ii_51211: DO ii=i-1,i
               l = ii-i+2
               kc = kint(i,j,k,l)
               IF (kc.EQ.0) THEN
                  tsint(k,l) = temp(ii,j,1)
               ELSEIF (kc.EQ.nz) THEN
                  tsint(k,l) = temp(ii,j,nz)
               ELSE
                  sigk = sigint(i,j,k,l)
                  tsint(k,l) = (1.0-sigk)*temp(ii,j,kc)+sigk*temp(ii,j,kc+1)
               ENDIF
            ENDDO ii_51211
            dhts(i,j,k) = tsint(k,2) - tsint(k,1)
         ENDIF
      ENDDO k_5121
      ENDDO i_5121
      ENDDO j_5121

!     ---interpolate arrays
      CALL Carr_at_U(beta_temp(0:ncloc,1:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)
      CALL Carr_at_UW(beta_temp(0:ncloc,1:nrloc,:),beta(:,:,2:nz),2,&
                   & (/0,1,1/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)

!     ---surface
      WHERE (masksuratu)
         array2d = 0.5*array2du*delzatu(1:ncloc,:,nz)*beta(:,:,nz+1)*&
                 & dhts(:,:,nz+1)
         p3dbcgradatu(:,:,nz) = array2d
      END WHERE
      
!     ---interior points
      k_5122: DO k=nz,2,-1
         WHERE (maskatu(:,:,k))
            array2d = array2d + array2du*beta(:,:,k)*dhts(:,:,k)*&
                    & delzatuw(1:ncloc,:,k)
            p3dbcgradatu(:,:,k-1) = array2d
         END WHERE
      ENDDO k_5122

   ENDIF

!
!5.1.3 Salinity
!--------------
!

   IF (salflag) THEN

!     ---horizontal salinity difference
      j_5131: DO j=1,nrloc
      i_5131: DO i=1,ncloc
      k_5131: DO k=2,nz+1
         IF ((k.EQ.(nz+1).AND.masksuratu(i,j)).OR.nodeatuw(i,j,k).EQ.2) THEN
            ii_51311: DO ii=i-1,i
               l = ii-i+2
               kc = kint(i,j,k,l)
               IF (kc.EQ.0) THEN
                  tsint(k,l) = sal(ii,j,1)
               ELSEIF (kc.EQ.nz) THEN
                  tsint(k,l) = sal(ii,j,nz)
               ELSE
                  sigk = sigint(i,j,k,l)
                  tsint(k,l) = (1.0-sigk)*sal(ii,j,kc)+sigk*sal(ii,j,kc+1)
               ENDIF
            ENDDO ii_51311
            dhts(i,j,k) = tsint(k,2) - tsint(k,1)
         ENDIF
      ENDDO k_5131
      ENDDO i_5131
      ENDDO j_5131

!     ---interpolate arrays
      CALL Carr_at_U(beta_sal(0:ncloc,1:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)
      CALL Carr_at_UW(beta_sal(0:ncloc,1:nrloc,:),beta(:,:,2:nz),2,&
                   & (/0,1,1/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)

!     ---surface
      WHERE (masksuratu)
         array2d = 0.5*array2du*delzatu(1:ncloc,:,nz)*beta(:,:,nz+1)*&
                 & dhts(:,:,nz+1)
         p3dbcgradatu(:,:,nz) = p3dbcgradatu(:,:,nz) - array2d
      END WHERE
      
!     ---interior points
      k_5132: DO k=nz,2,-1
         WHERE (maskatu(:,:,k))
            array2d = array2d + array2du*beta(:,:,k)*dhts(:,:,k)*&
                    & delzatuw(1:ncloc,:,k)
            p3dbcgradatu(:,:,k-1) = p3dbcgradatu(:,:,k-1) - array2d
         END WHERE
      ENDDO k_5132

   ENDIF

!
!5.1.4 Sediment
!--------------
!

   IF (sedimflag) CALL baroclinic_gradient_sed_z(sigint,kint,'X')

!
!5.2 Y-direction
!----------------
!
!5.2.1 Vertical levels for interpolation
!---------------------------------------
!

   j_521: DO j=1,nrloc
   i_521: DO i=1,ncloc
      IF (node2dv(i,j).EQ.2) THEN
         jj_5211: DO jj=j-1,j
            l = jj-j+2
            kc = 1
            k_52111: DO k=2,nz+1
               IF (k.EQ.nz+1) THEN
                  z = Zcoord_var(i,j,nz,'V  ',.FALSE.)
               ELSE
                  z = Zcoord_var(i,j,k,'VW ',.FALSE.)
               ENDIF
               sigk = (z+depmeanatc(i,jj))/deptotatc(i,jj)
               IF (sigk.LE.gscoordatc(i,jj,1)) THEN
                  kc = 0
                  sigint(i,j,k,l) = 0.0
               ELSEIF (sigk.GE.gscoordatc(i,jj,nz)) THEN
                  kc = nz
                  sigint(i,j,k,l) = 0.0
               ELSE
                  DO WHILE (sigk.GE.gscoordatc(i,jj,kc+1))
                     kc = kc + 1
                  END DO
                  sigint(i,j,k,l) = deptotatc(i,jj)*(sigk-gscoordatc(i,jj,kc))&
                                  & /delzatw(i,jj,kc+1)
               ENDIF
               kint(i,j,k,l) = kc
            ENDDO k_52111
         ENDDO jj_5211
      ENDIF
   ENDDO i_521
   ENDDO j_521

!
!5.2.2 Temperature
!-----------------
!

   IF (tempflag) THEN

!     ---horizontal temperature difference
      j_5221: DO j=1,nrloc
      i_5221: DO i=1,ncloc
      k_5221: DO k=2,nz+1
         IF ((k.EQ.(nz+1).AND.masksuratv(i,j)).OR.nodeatvw(i,j,k).EQ.2) THEN
            jj_5221: DO jj=j-1,j
               l = jj-j+2
               kc = kint(i,j,k,l)
               IF (kc.EQ.0) THEN
                  tsint(k,l) = temp(i,jj,1)
               ELSEIF (kc.EQ.nz) THEN
                  tsint(k,l) = temp(i,jj,nz)
               ELSE
                  sigk = sigint(i,j,k,l)
                  tsint(k,l) = (1.0-sigk)*temp(i,jj,kc)+sigk*temp(i,jj,kc+1)
               ENDIF
            ENDDO jj_5221
            dhts(i,j,k) = tsint(k,2) - tsint(k,1)
         ENDIF
      ENDDO k_5221
      ENDDO i_5221
      ENDDO j_5221

!     ---interpolate arrays
      CALL Carr_at_V(beta_temp(1:ncloc,0:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)
      CALL Carr_at_VW(beta_temp(1:ncloc,0:nrloc,:),beta(:,:,2:nz),2,&
                   & (/1,0,1/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)

!     ---surface
      WHERE (masksuratv)
         array2d = 0.5*array2du*delzatv(:,1:nrloc,nz)*beta(:,:,nz+1)*&
                  & dhts(:,:,nz+1)
         p3dbcgradatv(:,:,nz) = array2d
      END WHERE
      
!     ---interior points
      k_5222: DO k=nz,2,-1
         WHERE (maskatv(:,:,k))
            array2d = array2d + array2dv*beta(:,:,k)*dhts(:,:,k)*&
                    & delzatvw(:,1:nrloc,k)
            p3dbcgradatv(:,:,k-1) = array2d
         END WHERE
      ENDDO k_5222

   ENDIF

!
!5.2.3 Salinity
!--------------
!

   IF (salflag) THEN

!     ---horizontal salinity difference
      j_5231: DO j=1,nrloc
      i_5231: DO i=1,ncloc
      k_5231: DO k=2,nz+1
         IF ((k.EQ.(nz+1).AND.masksuratv(i,j)).OR.nodeatvw(i,j,k).EQ.2) THEN
            jj_52311: DO jj=j-1,j
               l = jj-j+2
               kc = kint(i,j,k,l)
               IF (kc.EQ.0) THEN
                  tsint(k,l) = sal(i,jj,1)
               ELSEIF (kc.EQ.nz) THEN
                  tsint(k,l) = sal(i,jj,nz)
               ELSE
                  sigk = sigint(i,j,k,l)
                  tsint(k,l) = (1.0-sigk)*sal(i,jj,kc)+sigk*sal(i,jj,kc+1)
               ENDIF
            ENDDO jj_52311
            dhts(i,j,k) = tsint(k,2) - tsint(k,1)
         ENDIF
      ENDDO k_5231
      ENDDO i_5231
      ENDDO j_5231

!     ---interpolate arrays
      CALL Carr_at_V(beta_sal(1:ncloc,0:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)
      CALL Carr_at_VW(beta_sal(1:ncloc,0:nrloc,:),beta(:,:,2:nz),2,&
                   & (/1,0,1/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)

!     ---surface
      WHERE (masksuratv)
         array2d = 0.5*array2dv*delzatv(:,1:nrloc,nz)*beta(:,:,nz+1)*&
                 & dhts(:,:,nz+1)
         p3dbcgradatv(:,:,nz) = p3dbcgradatv(:,:,nz) - array2d
      END WHERE
      
!     ---interior points
      k_5232: DO k=nz,2,-1
         WHERE (maskatv(:,:,k))
            array2d = array2d + array2dv*beta(:,:,k)*&
                    & dhts(:,:,k)*delzatvw(:,1:nrloc,k)
            p3dbcgradatv(:,:,k-1) = p3dbcgradatv(:,:,k-1) - array2d
         END WHERE
      ENDDO k_5232

   ENDIF

!
!5.2.4 Sediment
!--------------
!

   IF (sedimflag) CALL baroclinic_gradient_sed_z(sigint,kint,'Y')

!
!6. Cubic H
!----------
!

ELSEIF (iopt_dens_grad.EQ.3) THEN

!
!6.1 Z-coordinate and derivatives
!--------------------------------
!

   CALL Zcoord_arr(zcoord(1:ncloc,1:nrloc,:),(/1,1,1/),(/ncloc,nrloc,nz/),&
                & 'C  ',.FALSE.)
   IF (iopt_MPI.EQ.1) THEN
      nhexch = nhalo
      CALL exchange_mod(zcoord,(/1-nhalo,1-nhalo,1/),nhexch,0)
   ENDIF
   CALL hcube_deriv(zcoord,dzx,'X')
   CALL hcube_deriv(zcoord,dzy,'Y')
   CALL hcube_deriv(zcoord,dzz,'Z')

!
!6.2 Temperature
!---------------
!

   IF (tempflag) THEN

!
!6.2.1 Z-derivative and fluxes
!-----------------------------
!

      CALL hcube_deriv(temp,dzts,'Z')
      CALL hcube_fluxes(temp,zcoord,dzts,dzz,ftsatw,'Z')

!
!6.2.2 X-direction
!-----------------
!
!     ---interpolate arrays
      CALL Carr_at_U(beta_temp(0:ncloc,1:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)
      CALL Carr_at_UW(beta_temp(0:ncloc,1:nrloc,:),beta(:,:,2:nz),2,&
                   & (/0,1,1/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)

!     ---X-derivative and fluxes
      CALL hcube_deriv(temp,dxyts,'X')
      CALL hcube_fluxes(temp,zcoord,dxyts,dzx,ftsatvel,'X')

!     ---surface
      WHERE (masksuratu)
         array2d = 0.5*array2du*delzatu(1:ncloc,:,nz)*beta(:,:,nz+1)*&
                 & (temp(1:ncloc,1:nrloc,nz)-temp(0:ncloc-1,1:nrloc,nz))
         p3dbcgradatu(:,:,nz) = array2d
      END WHERE

!     ---interior points
      k_6222: DO k=nz,2,-1
         WHERE (maskatu(:,:,k))
            array2d = array2d + array2du*beta(:,:,k)*&
                  & (ftsatw(1:ncloc,1:nrloc,k)-ftsatw(0:ncloc-1,1:nrloc,k)&
                  & +ftsatvel(1:ncloc,1:nrloc,k-1)-ftsatvel(1:ncloc,1:nrloc,k))
            p3dbcgradatu(:,:,k-1) = array2d
         END WHERE
      ENDDO k_6222

!
!6.2.3 Y-direction
!-----------------
!
!     ---interpolate arrays
      CALL Carr_at_V(beta_temp(1:ncloc,0:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)
      CALL Carr_at_VW(beta_temp(1:ncloc,0:nrloc,:),beta(:,:,2:nz),2,&
                   & (/1,0,1/),(/ncloc,nrloc,nz/),1,iarr_beta_temp,.TRUE.)

!     ---X-derivative and fluxes
      CALL hcube_deriv(temp,dxyts,'Y')
      CALL hcube_fluxes(temp,zcoord,dxyts,dzy,ftsatvel,'Y')

!     ---surface
      WHERE (masksuratv)
         array2d = 0.5*array2dv*delzatv(:,1:nrloc,nz)*beta(:,:,nz+1)*&
                 & (temp(1:ncloc,1:nrloc,nz)-temp(1:ncloc,0:nrloc-1,nz))
         p3dbcgradatv(:,:,nz) = array2d
      END WHERE

!     ---interior points
      k_6232: DO k=nz,2,-1
         WHERE (maskatv(:,:,k))
            array2d = array2d + array2dv*beta(:,:,k)*&
                  & (ftsatw(1:ncloc,1:nrloc,k)-ftsatw(1:ncloc,0:nrloc-1,k)&
                  & +ftsatvel(1:ncloc,1:nrloc,k-1)-ftsatvel(1:ncloc,1:nrloc,k))
            p3dbcgradatv(:,:,k-1) = array2d
         END WHERE
      ENDDO k_6232

   ENDIF

!
!6.3 Salinity
!------------
!

   IF (salflag) THEN

!
!6.3.1 Z-derivative and fluxes
!-----------------------------
!

      CALL hcube_deriv(sal,dzts,'Z')
      CALL hcube_fluxes(sal,zcoord,dzts,dzz,ftsatw,'Z')

!
!6.3.2 X-direction
!-----------------
!
!     ---interpolate arrays
      CALL Carr_at_U(beta_sal(0:ncloc,1:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/0,1,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)
      CALL Carr_at_UW(beta_sal(0:ncloc,1:nrloc,:),beta(:,:,2:nz),2,&
                   & (/0,1,1/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)

!     ---X-derivative and fluxes
      CALL hcube_deriv(sal,dxyts,'X')
      CALL hcube_fluxes(sal,zcoord,dxyts,dzx,ftsatvel,'X')

!     ---surface
      WHERE (masksuratu)
         array2d = 0.5*array2du*delzatu(1:ncloc,:,nz)*beta(:,:,nz+1)*&
                 & (sal(1:ncloc,1:nrloc,nz)-sal(0:ncloc-1,1:nrloc,nz))
         p3dbcgradatu(:,:,nz) = p3dbcgradatu(:,:,nz) - array2d
      END WHERE

!     ---interior points
      k_6322: DO k=nz,2,-1
         WHERE (maskatu(:,:,k))
            array2d = array2d + array2du*beta(:,:,k)*&
                  & (ftsatw(1:ncloc,1:nrloc,k)-ftsatw(0:ncloc-1,1:nrloc,k)&
                  & +ftsatvel(1:ncloc,1:nrloc,k-1)-ftsatvel(1:ncloc,1:nrloc,k))
            p3dbcgradatu(:,:,k-1) = p3dbcgradatu(:,:,k-1) - array2d
         END WHERE
      ENDDO k_6322

!
!6.3.3 Y-direction
!-----------------
!
!     ---interpolate arrays
      CALL Carr_at_V(beta_sal(1:ncloc,0:nrloc,nz),beta(:,:,nz+1),1,2,&
                  & (/1,0,nz/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)
      CALL Carr_at_VW(beta_sal(1:ncloc,0:nrloc,:),beta(:,:,2:nz),2,&
                   & (/1,0,1/),(/ncloc,nrloc,nz/),1,iarr_beta_sal,.TRUE.)

!     ---X-derivative and fluxes
      CALL hcube_deriv(sal,dxyts,'Y')
      CALL hcube_fluxes(sal,zcoord,dxyts,dzy,ftsatvel,'Y')

!     ---surface
      WHERE (masksuratv)
         array2d = 0.5*array2dv*delzatv(:,1:nrloc,nz)*beta(:,:,nz+1)*&
                 & (sal(1:ncloc,1:nrloc,nz)-sal(1:ncloc,0:nrloc-1,nz))
         p3dbcgradatv(:,:,nz) = p3dbcgradatv(:,:,nz) - array2d
      END WHERE

!     ---interior points
      k_6332: DO k=nz,2,-1
         WHERE (maskatv(:,:,k))
            array2d = array2d + array2dv*beta(:,:,k)*&
                  & (ftsatw(1:ncloc,1:nrloc,k)-ftsatw(1:ncloc,0:nrloc-1,k)&
                  & +ftsatvel(1:ncloc,1:nrloc,k-1)-ftsatvel(1:ncloc,1:nrloc,k))
            p3dbcgradatv(:,:,k-1) = p3dbcgradatv(:,:,k-1) - array2d
         END WHERE
      ENDDO k_6332

   ENDIF

!
!6.3.4 Sediment
!--------------
!

   IF (sedimflag) THEN
      CALL baroclinic_gradient_sed_cubic (zcoord,dzx,dzy,dzz,'X')
      CALL baroclinic_gradient_sed_cubic (zcoord,dzx,dzy,dzz,'Y')
   ENDIF

ENDIF

!
!7. Depth-integrated gradient
!----------------------------
!

WHERE (node2du(1:ncloc,1:nrloc).EQ.2)
   p2dbcgradatu = SUM(p3dbcgradatu*delzatu(1:ncloc,:,:),DIM=3)
END WHERE

WHERE (node2dv(1:ncloc,1:nrloc).EQ.2)
   p2dbcgradatv = SUM(p3dbcgradatv*delzatv(:,1:nrloc,:),DIM=3)
END WHERE

!
!8. Apply reduction factor for flooding/drying scheme
!----------------------------------------------------
!

IF (iopt_fld_alpha.GT.0) THEN

!
!8.1 U-nodes
!-----------
!
   WHERE (node2du(1:ncloc,1:nrloc).EQ.2)
      p2dbcgradatu = alphatu_fld*p2dbcgradatu
   END WHERE
   k_810: DO k=1,nz
      WHERE (maskatu(:,:,k))
         p3dbcgradatu(:,:,k) = alphatu_fld*p3dbcgradatu(:,:,k)
      END WHERE
   ENDDO k_810

!
!8.2 V-nodes
!-----------
!
   WHERE (node2dv(1:ncloc,1:nrloc).EQ.2)
      p2dbcgradatv = alphatv_fld*p2dbcgradatv
   END WHERE
   k_820: DO k=1,nz
      WHERE (maskatv(:,:,k))
         p3dbcgradatv(:,:,k) = alphatv_fld*p3dbcgradatv(:,:,k)
      END WHERE
   ENDDO k_820

ENDIF

!
!9. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatu,maskatv,masksuratu,masksuratv,array2d,array2du,array2dv,&
          & beta)
IF (iopt_dens_grad.EQ.1.OR.iopt_dens_grad.EQ.3) DEALLOCATE (zcoord)
IF (iopt_dens_grad.EQ.1) DEALLOCATE (tsatw,tsatvel)
IF (iopt_dens_grad.EQ.2) DEALLOCATE (dhts,kint,sigint)
IF (iopt_dens_grad.EQ.3) THEN
   DEALLOCATE (dxyts,dzz,dzx,dzy,dzts,ftsatw,ftsatvel)
ENDIF

CALL log_timer_out(npcc,itm_phgrad)


RETURN

END SUBROUTINE baroclinic_gradient

!========================================================================

SUBROUTINE hcube_deriv(psi,dpsi,cdir)
!************************************************************************
!
! *hcube_deriv* Harmonic derivative of array psi
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Density_Equations.f90  V2.0
!
! Description - 
!
! Reference - Shchepetkin A.F. and McWilliams J.C., 2003. A method for
!    computing horizontal pressure-gradient force in an oceanic model with a
!    nonaligned vertical coordinate. Journal of Geophysical Research, 108,
!    3090, doi:10.1029/2001JC001407.
!
! Calling program - baroclinic_gradient
!
! Module calls - exchange_mod
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE switches
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
 
!
!*  Arguments
!
CHARACTER (LEN=1), INTENT(IN) :: cdir
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz) :: psi
REAL, INTENT(OUT), DIMENSION(0:ncloc+1,0:nrloc+1,nz) :: dpsi

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*psi*      REAL    Input array
!*dpsi*     REAL    Array of derivatives
!*cdir*     CHAR    Direction of derivative ('X','Y','Z')
!
!------------------------------------------------------------------------------
!
!*Local variables
!     
INTEGER :: i, j, k
REAL :: dprodx1, dprodx2, xfac
INTEGER, DIMENSION(4) :: nhexch
REAL, DIMENSION(nz) :: dprod1, dprod2


procname(pglev+1) = 'hcube_deriv'
CALL log_timer_in()

!
!1. Initialise
!-------------
!

xfac = 7.0/15.0
dpsi = 0.0

!
!2. X-direction
!--------------
!

IF (cdir.EQ.'X') THEN

!  ---internal points
   j_210: DO j=1,nrloc
      i_211: DO i=0,ncloc+1
         IF (ALL(nodeatc(i-1:i+1,j).GT.0)) THEN
            dprod1 = psi(i+1,j,:) - psi(i,j,:)
            dprod2 = psi(i,j,:) - psi(i-1,j,:)
            WHERE (dprod1*dprod2.GT.0.0)
               dpsi(i,j,:) = 2.0*dprod1*dprod2/(dprod1+dprod2)
            END WHERE
         ENDIF
      ENDDO i_211
      i_212: DO i=1,ncloc
         IF (nodeatc(i-1,j).EQ.0.AND.ALL(nodeatc(i:i+1,j).GT.0)) THEN
            dpsi(i,j,:) = MERGE(1.2*(psi(i+1,j,:)-psi(i,j,:))-&
                              & xfac*dpsi(i+1,j,:),&
                              & psi(i+1,j,:)-psi(i,j,:),dpsi(i+1,j,:).NE.0.0)
         ELSEIF (ALL(nodeatc(i-1:i,j).GT.0).AND.nodeatc(i+1,j).EQ.0) THEN
            dpsi(i,j,:) = MERGE(1.2*(psi(i,j,:)-psi(i-1,j,:))-&
                              & xfac*dpsi(i-1,j,:),&
                              & psi(i,j,:)-psi(i-1,j,:),dpsi(i-1,j,:).NE.0.0)
         ENDIF
      ENDDO i_212
   ENDDO j_210

!  ---halo exchange
   IF (iopt_MPI.EQ.1) THEN
      nhexch = (/1,0,0,0/)
      CALL exchange_mod(dpsi,(/0,0,1/),nhexch,0)
   ENDIF

!
!3. Y-direction
!--------------
!

ELSEIF (cdir.EQ.'Y') THEN

!  ---internal points
   i_310: DO i=1,ncloc
      j_311: DO j=0,nrloc+1
         IF (ALL(nodeatc(i,j-1:j+1).GT.0)) THEN
            dprod1 = psi(i,j+1,:) - psi(i,j,:)
            dprod2 = psi(i,j,:) - psi(i,j-1,:)
            WHERE (dprod1*dprod2.GT.0.0)
               dpsi(i,j,:) = 2.0*dprod1*dprod2/(dprod1+dprod2)
            END WHERE
         ENDIF
      ENDDO j_311
      j_312: DO j=1,nrloc
         IF (nodeatc(i,j-1).EQ.0.AND.ALL(nodeatc(i,j:j+1).GT.0)) THEN
            dpsi(i,j,:) = MERGE(1.2*(psi(i,j+1,:)-psi(i,j,:))-&
                              & xfac*dpsi(i,j+1,:),&
                              & psi(i,j+1,:)-psi(i,j,:),dpsi(i,j+1,:).NE.0.0)
         ELSEIF (ALL(nodeatc(i,j-1:j).GT.0).AND.nodeatc(i,j+1).EQ.0) THEN
            dpsi(i,j,:) = MERGE(1.2*(psi(i,j,:)-psi(i,j-1,:))-&
                              & xfac*dpsi(i,j-1,:),&
                              & psi(i,j,:)-psi(i,j-1,:),dpsi(i,j-1,:).NE.0.0)
         ENDIF
      ENDDO j_312
   ENDDO i_310

!  ---halo exchange
   IF (iopt_MPI.EQ.1) THEN
      nhexch = (/0,0,1,0/)
      CALL exchange_mod(dpsi,(/0,0,1/),nhexch,0)
   ENDIF

!
!4. Z-direction
!--------------
!

ELSEIF (cdir.EQ.'Z') THEN
   j_410: DO j=0,nrloc
   i_410: DO i=0,ncloc
      IF (nodeatc(i,j).GT.0) THEN
         k_411: DO k=2,nz-1
            dprodx1 = psi(i,j,k+1) - psi(i,j,k)
            dprodx2 = psi(i,j,k) - psi(i,j,k-1)
            IF (dprodx1*dprodx2.GT.0.0) THEN
               dpsi(i,j,k) = 2.0*dprodx1*dprodx2/(dprodx1+dprodx2)
            ENDIF
         ENDDO k_411
         dpsi(i,j,1) = 1.2*(psi(i,j,2)-psi(i,j,1))-xfac*dpsi(i,j,2)
         dpsi(i,j,nz) = 1.2*(psi(i,j,nz)-psi(i,j,nz-1))-xfac*dpsi(i,j,nz-1)
      ENDIF
   ENDDO i_410
   ENDDO j_410
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE hcube_deriv

!========================================================================

SUBROUTINE hcube_fluxes(psi,zcoord,dpsi,dzcoord,fcube,cdir)
!************************************************************************
!
! *hcube_deriv* Pseudo fluxes for baroclinic pressure gradient
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Density_Equations.f90  V2.2
!
! Description - 
!
! Reference - Shchepetkin A.F. and McWilliams J.C., 2003. A method for
!    computing horizontal pressure-gradient force in an oceanic model with a
!    nonaligned vertical coordinate. Journal of Geophysical Research, 108,
!    3090, doi:10.1029/2001JC001407.
!
! Calling program - baroclinic_gradient
!
! Module calls - error_alloc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE syspars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*  Arguments
!
CHARACTER (LEN=1), INTENT(IN) :: cdir
REAL, INTENT(IN), DIMENSION(1-nhalo:ncloc+nhalo,&
                          & 1-nhalo:nrloc+nhalo,nz) :: psi, zcoord
REAL, INTENT(IN), DIMENSION(0:ncloc+1,0:nrloc+1,nz) :: dpsi, dzcoord
REAL, INTENT(OUT), DIMENSION(0:ncloc,0:nrloc,nz) :: fcube

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*psi*      REAL    Density (temperature or salinity) array
!*zcoord*   REAL    Array of Z-coordinates
!*dpsi*     REAL    Array of psi-derivatives
!*dzcoord*  REAL    Array of Z-derivatives
!*fcube*    REAL    Pseudo fluxes
!*cdir*     CHAR    Direction of derivative ('X','Y','Z')
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2, array2d3


procname(pglev+1) = 'hcube_fluxes'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

SELECT CASE (cdir)
   CASE ('X','Y')
      ALLOCATE (array2d1(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
      ALLOCATE (array2d2(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)
      ALLOCATE (array2d3(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('array2d3',2,(/ncloc,nrloc/),kndrtype)
   CASE ('Z')
      ALLOCATE (array2d1(0:ncloc,0:nrloc),STAT=errstat)
      CALL error_alloc('array2d1',2,(/ncloc+1,nrloc+1/),kndrtype)
      ALLOCATE (array2d2(0:ncloc,0:nrloc),STAT=errstat)
      CALL error_alloc('array2d2',2,(/ncloc+1,nrloc+1/),kndrtype)
      ALLOCATE (array2d3(0:ncloc,0:nrloc),STAT=errstat)
      CALL error_alloc('array2d3',2,(/ncloc+1,nrloc+1/),kndrtype)
END SELECT

!
!2. X-direction
!--------------
!

IF (cdir.EQ.'X') THEN

   k_210: DO k=1,nz
      WHERE (nodeatu(1:ncloc,1:nrloc,k).EQ.2)
         array2d2 = zcoord(1:ncloc,1:nrloc,k)-zcoord(0:ncloc-1,1:nrloc,k)
         array2d1 = (psi(0:ncloc-1,1:nrloc,k)+psi(1:ncloc,1:nrloc,k))*array2d2
         array2d2 = array2d2 - &
               & (dzcoord(0:ncloc-1,1:nrloc,k)+dzcoord(1:ncloc,1:nrloc,k))/12.0
         array2d3 = psi(1:ncloc,1:nrloc,k)-psi(0:ncloc-1,1:nrloc,k) - &
                  & (dpsi(0:ncloc-1,1:nrloc,k)+dpsi(1:ncloc,1:nrloc,k))/12.0
         fcube(1:ncloc,1:nrloc,k) = 0.5*array2d1 - 0.1*&
          & ((dpsi(1:ncloc,1:nrloc,k)-dpsi(0:ncloc-1,1:nrloc,k))*array2d2-& 
          & (dzcoord(1:ncloc,1:nrloc,k)-dzcoord(0:ncloc-1,1:nrloc,k))*array2d3)
      END WHERE
   ENDDO k_210

!
!3. Y-direction
!--------------
!

ELSEIF (cdir.EQ.'Y') THEN

   k_310: DO k=1,nz
      WHERE (nodeatv(1:ncloc,1:nrloc,k).EQ.2)
         array2d2 = zcoord(1:ncloc,1:nrloc,k)-zcoord(1:ncloc,0:nrloc-1,k)
         array2d1 = (psi(1:ncloc,0:nrloc-1,k)+psi(1:ncloc,1:nrloc,k))*array2d2
         array2d2 = array2d2 - &
               & (dzcoord(1:ncloc,0:nrloc-1,k)+dzcoord(1:ncloc,1:nrloc,k))/12.0
         array2d3 = psi(1:ncloc,1:nrloc,k)-psi(1:ncloc,0:nrloc-1,k) - &
                  & (dpsi(1:ncloc,0:nrloc-1,k)+dpsi(1:ncloc,1:nrloc,k))/12.0
         fcube(1:ncloc,1:nrloc,k) = 0.5*array2d1 - 0.1*&
          & ((dpsi(1:ncloc,1:nrloc,k)-dpsi(1:ncloc,0:nrloc-1,k))*array2d2-& 
          & (dzcoord(1:ncloc,1:nrloc,k)-dzcoord(1:ncloc,0:nrloc-1,k))*array2d3)
      END WHERE
   ENDDO k_310

!
!4. Z-direction
!--------------
!

ELSEIF (cdir.EQ.'Z') THEN

   k_410: DO k=2,nz
      WHERE (nodeatc(0:ncloc,0:nrloc).GT.0)
         array2d2 = zcoord(0:ncloc,0:nrloc,k) - zcoord(0:ncloc,0:nrloc,k-1)
         array2d1 = (psi(0:ncloc,0:nrloc,k-1)+psi(0:ncloc,0:nrloc,k))*array2d2
         array2d2 = array2d2 - (dzcoord(0:ncloc,0:nrloc,k-1)+&
                              & dzcoord(0:ncloc,0:nrloc,k))/12.0
         array2d3 = psi(0:ncloc,0:nrloc,k)-psi(0:ncloc,0:nrloc,k-1) - &
                 & (dpsi(0:ncloc,0:nrloc,k-1)+dpsi(0:ncloc,0:nrloc,k))/12.0
         fcube(:,:,k) = 0.5*array2d1 - 0.1*&
              & ((dpsi(0:ncloc,0:nrloc,k)-dpsi(0:ncloc,0:nrloc,k-1))*array2d2-&
          & (dzcoord(0:ncloc,0:nrloc,k)-dzcoord(0:ncloc,0:nrloc,k-1))*array2d3)
      END WHERE
   ENDDO k_410

ENDIF

!
!5. Deallocate arrays
!--------------------
!

DEALLOCATE (array2d1,array2d2,array2d3)

CALL log_timer_out()


RETURN

END SUBROUTINE hcube_fluxes

!========================================================================

SUBROUTINE convect_adjust
!************************************************************************
!
! *convect_adjust* Adjust density and active tracers for convective transport
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Density_Equations.f90  V2.11.2
!
! Description -
!
! Calling program - equation_of_state
!
! External calls -
!
! Module calls - exchange_mod, mult_index
!
!************************************************************************
!
USE density  
USE grid
USE gridpars  
USE iopars
USE modids
USE switches
use timepars
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
use utility_routines, only: mult_index

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, iter, j, k, klo, kup
INTEGER, DIMENSION(4) :: nhexch
REAL :: adelz, adens, asal, atemp, delzlo
LOGICAL, DIMENSION(ncloc,nrloc) :: densflag
REAL, DIMENSION(0:nz) :: dens1d


procname(pglev+1) = 'convect_adjust'
CALL log_timer_in()

!
!1. Check for unstable profiles
!------------------------------
!

j_110: DO j=1,nrloc
i_110: DO i=1,ncloc
   densflag(i,j) = .FALSE.
   IF (maskatc_int(i,j)) THEN
      k_111: DO k=2,nz
         IF ((dens(i,j,k)-dens(i,j,k-1).GT.0.0)) densflag(i,j) = .TRUE.
      ENDDO k_111
   ENDIF
ENDDO i_110
ENDDO j_110

!
!2. Adjust unstable profiles
!---------------------------
!

dens1d = 0.0

j_200: DO j=1,nrloc
i_200: DO i=1,ncloc
   IF (densflag(i,j)) THEN
      dens1d(1:nz) = dens(i,j,:)

      iter_210: DO iter=1,2*nz
         kup = nz
         DO WHILE (kup.GT.1.AND.(dens1d(kup)-dens1d(kup-1)).LE.0.0)
            kup = kup - 1
         ENDDO
         IF (kup.EQ.1) GOTO 99
         adelz = delzatc(i,j,kup); atemp = temp(i,j,kup)
         asal = sal(i,j,kup); adens = dens1d(kup)
         k_211: DO k=kup-1,1,-1
            IF (adens.LE.dens1d(k)) THEN
               klo = k
               GOTO 98
            ENDIF
            delzlo = delzatc(i,j,k)
            adelz = adelz + delzlo
            atemp = (atemp*(adelz-delzlo)+temp(i,j,k)*delzlo)/adelz
            asal = (asal*(adelz-delzlo)+sal(i,j,k)*delzlo)/adelz
            adens = (adens*(adelz-delzlo)+dens1d(k)*delzlo)/adelz
         ENDDO k_211
         klo = 1
98       CONTINUE
         k_212: DO k=kup,klo+1,-1
            temp(i,j,k) = atemp; sal(i,j,k) = asal; dens1d(k) = adens
         ENDDO k_212
         IF (klo.EQ.1.AND.adens.GE.dens1d(klo)) THEN
            temp(i,j,klo) = atemp; sal(i,j,klo) = asal; dens1d(klo) = adens
         ENDIF
      ENDDO iter_210
99    CONTINUE
      dens(i,j,:) = dens1d(1:nz)
      
   ENDIF

ENDDO i_200
ENDDO j_200

!
!3. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = nhalo
   CALL exchange_mod(sal,(/1-nhalo,1-nhalo,1/),nhexch,iarr_sal)
   CALL exchange_mod(temp,(/1-nhalo,1-nhalo,1/),nhexch,iarr_temp)
   nhexch = 1
   CALL exchange_mod(dens,(/0,0,1/),nhexch,iarr_dens)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE convect_adjust
