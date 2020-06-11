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
! *Hydrodynamic_Equations* Hydrodynamic (currents and elevations) calculations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.11.2
!
! $Date: 2018-09-28 13:01:37 +0200 (Fri, 28 Sep 2018) $
!
! $Revision: 1188 $
!
! Description -
!
! Routines -  correct_free_surf, current_corr, current_corr_1d, current_pred,
!             current_2d, hydrodynamic_equations, physical_vertical_current,
!             shear_frequency, surface_elevation, transf_vertical_current,
!             transf_vertical_fall
!
!************************************************************************
!

!========================================================================

SUBROUTINE hydrodynamic_equations
!************************************************************************
!
! *hydrodynamic_equations* Solve momentum and continuity equations
!
! Author - Pieter Rauwoens and Patrick Luyten
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.10.2
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - bottom_stress, correct_free_surf, current_corr,
!                  current_pred, current_2d, drying_factor,
!                  open_boundary_outflow_time, physical_vertical_current,
!                  surface_stress, transf_vertical_current, update_nest_data_2d,
!                  update_data_nest_3d
!
! Module calls - exchange_mod, maxloc_vars, sum_vars
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE physpars
USE switches
USE timepars
USE wavevars
USE paral_comms, ONLY: exchange_mod
USE paral_utilities, ONLY: maxloc_vars, sum_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: cmean, cmin
CHARACTER (LEN=12), DIMENSION(3) :: cmax
INTEGER :: npcc
INTEGER, DIMENSION(3) :: lbounds, maxatu, maxatv
INTEGER, DIMENSION(4) :: nhexch
REAL :: umax, vmax


procname(pglev+1) = 'hydrodynamic_equations'
CALL log_timer_in(npcc)


!
!1. Mode Splitting
!-----------------
!

IF (iopt_hydro_impl.EQ.0) THEN
   
!  ---predictor step
   IF (predstep.AND.iopt_mode_3D.EQ.2) CALL current_pred

!  ---barotropic step
   IF (iopt_mode_2D.EQ.2) CALL current_2d

!  ---corrector step
   IF (corrstep.AND.iopt_mode_3D.EQ.2) CALL current_corr

!
!2. Implicit method
!------------------
!

ELSEIF (iopt_hydro_impl.EQ.1) THEN

!  ---initialise monitoring
   IF (nt.EQ.1.AND.monflag) THEN
      maxmgiter = 0; minmgiter = nomgiterations; meanmgiter = 0.0
   ENDIF

!  ---monitoring
   IF (monflag) WRITE (iomon,'(A)') CDateTime

!  ---store old values
   zeta_old = zeta

   itsimp_210: DO itsimp=1,maxitsimp

!     ---predictor step
      IF (iopt_mode_3D.EQ.2) THEN
         CALL current_pred
      ELSEIF (iopt_mode_2D.EQ.2) THEN
         CALL current_2d
      ENDIF

!     ---free surface corrector step
      CALL correct_free_surf

!     ---monitoring 
      IF (monflag.AND.itsimp.EQ.1) THEN
         WRITE (iomon,9001) 'O: ', itsimp, dzetaresid
         maxmgiter = MAX(mgiteration,maxmgiter)
         minmgiter = MIN(mgiteration,minmgiter)
         meanmgiter = meanmgiter + mgiteration
      ENDIF

!     ---test convergence
      IF (dzetaresid.LT.dzetaresid_conv) THEN
         IF (monflag) WRITE (iomon,*)
         EXIT itsimp_210
      ENDIF

   ENDDO itsimp_210

   noitsimp = itsimp

!  ---corrector step
   IF (iopt_mode_3D.EQ.2) CALL current_corr

!  ---filtered currents
   IF (iopt_grid_nodim.EQ.3) THEN
      ufvel = uvel(2-nhalo:ncloc+nhalo,1:nrloc,:)
      vfvel = vvel(1:ncloc,2-nhalo:nrloc+nhalo,:)
      IF (iopt_grid_vtype.EQ.3) THEN
         udfvel = udvel(1:ncloc+1,1:nrloc)
         vdfvel = udvel(1:ncloc,1:nrloc+1)
      ENDIF
   ENDIF


ENDIF

!
!3. Vertical current
!-------------------
!

IF (iopt_grid_nodim.EQ.3.AND.corrstep) THEN

!  ---transformed current
   CALL transf_vertical_current
   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/0,0,1/); nhexch = (/1,0,1,0/)
      CALL exchange_mod(wvel,lbounds,nhexch,iarr_wvel)
      IF (iopt_waves_curr.EQ.1) THEN
         CALL exchange_mod(wstokesatw,lbounds,nhexch,iarr_wvel)
      ENDIF
   ENDIF

!  ---add Stokes drift to advective velocities     
   IF (iopt_waves_curr.EQ.1) THEN
      ufvel = ufvel + ustokesatu(2-nhalo:ncloc+nhalo,1:nrloc,:)
      vfvel = vfvel + vstokesatv(1:ncloc,2-nhalo:nrloc+nhalo,:)
      wfvel = wvel + wstokesatw
   ELSE
      wfvel = wvel
   ENDIF
   
!  ---physical current
   CALL physical_vertical_current
 
ENDIF
  
!
!4. "3-D" currents for 2-D case
!------------------------------
!

IF (iopt_grid_nodim.EQ.2) THEN
   ufvel(:,:,1) = umvel(2-nhalo:ncloc+nhalo,1:nrloc)
   uvel(:,0:nrloc+1,1) = umvel
   vfvel(:,:,1) = vmvel(1:ncloc,2-nhalo:nrloc+nhalo)
   vvel(0:ncloc+1,:,1) = vmvel
   IF (iopt_waves_curr.EQ.1) THEN
      ufvel(:,:,1) = ufvel(:,:,1) + umstokesatu(2-nhalo:ncloc+nhalo,1:nrloc)
      vfvel(:,:,1) = vfvel(:,:,1) + vmstokesatv(1:ncloc,2-nhalo:nrloc+nhalo)
   ENDIF
ENDIF

!
!5. Update bottom/surface stress
!--------------------------------
!

IF (iopt_grid_nodim.NE.3.OR.corrstep) THEN
   IF (iopt_bstres_form.GT.0) CALL bottom_stress
   IF (iopt_meteo_stres.EQ.1) CALL surface_stress
ENDIF

!
!6. Flooding/drying scheme
!-------------------------
!

IF (iopt_fld.GT.0.AND.corrstep) CALL drying_factor

!
!7. Write at nest locations
!--------------------------
!

IF (iopt_nests.EQ.1) THEN
   CALL update_nest_data_2d(io_2uvnst)
   CALL update_nest_data_2d(io_2xynst)
   IF (corrstep) THEN
      CALL update_nest_data_3d(io_3uvnst)
      CALL update_nest_data_3d(io_3xynst)
   ENDIF
ENDIF

!
!8. Return times at open boundaries
!----------------------------------
!

IF (iopt_obc_th.EQ.1) CALL open_boundary_outflow_time 

!
!9. Set land mask at emptied cells
!---------------------------------
!

IF (iopt_fld.EQ.2.AND.corrstep) THEN

   WHERE (deptotatc(1:ncloc,1:nrloc).EQ.dmin_fld)
      nodeatc(1:ncloc,1:nrloc) = 0
      maskatc_int = .FALSE.
      zeta(1:ncloc,1:nrloc) = dmin_fld - depmeanatc(1:ncloc,1:nrloc)
   END WHERE

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds(1:2) = 1-nhalo; nhexch = nhalo
      CALL exchange_mod(nodeatc,lbounds(1:2),nhexch,iarr_nodeatc)
      lbounds(1:2) = 0; nhexch = 1
      CALL exchange_mod(zeta,lbounds(1:2),nhexch,iarr_zeta)
   ENDIF
   
ENDIF

!
!10. Number of active points
!---------------------------
!

IF (iopt_fld.GT.0) THEN

!  ---C-nodes
   nowetatcloc = COUNT(maskatc_int)
   CALL sum_vars(nowetatcloc,nowetatc,0,commall=.TRUE.)

!  ---U-nodes
   nowetatuloc =  COUNT(node2du(1:ncloc,1:nrloc).GT.1)
   CALL sum_vars(nowetatuloc,nowetatu,0,commall=.TRUE.)

!  ---V-nodes
   nowetatvloc =  COUNT(node2dv(1:ncloc,1:nrloc).GT.1)
   CALL sum_vars(nowetatvloc,nowetatv,0,commall=.TRUE.)

ENDIF

!
!11. Write max-values to log-file
!--------------------------------
!
!  ---2-D case
IF (iopt_grid_nodim.GT.1) THEN
   CALL maxloc_vars(ABS(umvel(1:ncloc,1:nrloc)),umax,maxatu(1:2),iarr_umvel,&
                      & mask=node2du(1:ncloc,1:nrloc).GT.1,commall=.TRUE.)
   CALL maxloc_vars(ABS(vmvel(1:ncloc,1:nrloc)),vmax,maxatv(1:2),iarr_vmvel,&
                      & mask=node2dv(1:ncloc,1:nrloc).GT.1,commall=.TRUE.)
   IF (master.AND.loglev1.GE.0) THEN
      IF (nc.GT.2.AND.ABS(umax-float_fill).GT.float_fill_eps) THEN
         WRITE (cmax(1),'(I12)') maxatu(1); cmax(1) = ADJUSTL(cmax(1))
         WRITE (cmax(2),'(I12)') maxatu(2); cmax(2) = ADJUSTL(cmax(2))
         WRITE (iolog,9002) REPEAT(' ',pglev-1)//'umvelmax = ', umax, &
                          & TRIM(cmax(1)), TRIM(cmax(2))
      ENDIF
      IF (nr.GT.2.AND.ABS(vmax-float_fill).GT.float_fill_eps) THEN
         WRITE (cmax(1),'(I12)') maxatv(1); cmax(1) = ADJUSTL(cmax(1))
         WRITE (cmax(2),'(I12)') maxatv(2); cmax(2) = ADJUSTL(cmax(2))
         WRITE (iolog,9002) REPEAT(' ',pglev-1)//'vmvelmax = ', vmax, &
                          & TRIM(cmax(1)), TRIM(cmax(2))
      ENDIF

   ENDIF

ENDIF

!---3-D case
IF (corrstep.AND.iopt_grid_nodim.NE.2) THEN
   CALL maxloc_vars(ABS(uvel(1:ncloc,1:nrloc,:)),umax,maxatu,iarr_uvel,&
                      & mask=node2du(1:ncloc,1:nrloc).GT.1)
   CALL maxloc_vars(ABS(vvel(1:ncloc,1:nrloc,:)),vmax,maxatv,iarr_vvel,&
                      & mask=node2dv(1:ncloc,1:nrloc).GT.1)
   IF (master.AND.loglev1.GE.0) THEN
      IF (nc.GT.2.AND.ABS(umax-float_fill).GT.float_fill_eps) THEN
         WRITE (cmax(1),'(I12)') maxatu(1); cmax(1) = ADJUSTL(cmax(1))
         WRITE (cmax(2),'(I12)') maxatu(2); cmax(2) = ADJUSTL(cmax(2))
         WRITE (cmax(3),'(I12)') maxatu(3); cmax(3) = ADJUSTL(cmax(3))
         WRITE (iolog,9003) REPEAT(' ',pglev-1)//'uvelmax = ', umax, &
                               & TRIM(cmax(1)), TRIM(cmax(2)), TRIM(cmax(3))
      ENDIF
      IF (nr.GT.2.AND.ABS(vmax-float_fill).GT.float_fill_eps) THEN
         WRITE (cmax(1),'(I12)') maxatv(1); cmax(1) = ADJUSTL(cmax(1))
         WRITE (cmax(2),'(I12)') maxatv(2); cmax(2) = ADJUSTL(cmax(2))
         WRITE (cmax(3),'(I12)') maxatv(3); cmax(3) = ADJUSTL(cmax(3))
         WRITE (iolog,9003) REPEAT(' ',pglev-1)//'vvelmax = ', vmax, &
                                & TRIM(cmax(1)), TRIM(cmax(2)), TRIM(cmax(3)) 
      ENDIF
   ENDIF
ENDIF

!
!12. Monitoring
!--------------
!

IF (monflag.AND.nt.EQ.nstep) THEN
   meanmgiter = meanmgiter/nstep
   WRITE (cmax(1),'(I12)') maxmgiter; cmax(1) = ADJUSTL(cmax(1))
   WRITE (cmin,'(I12)') minmgiter; cmin = ADJUSTL(cmin)
   WRITE (cmean,'(F5.1)') meanmgiter; cmean = ADJUSTL(cmean)
   WRITE (iomon,*)
   WRITE (iomon,'(A)') 'Minimum no of iterations : '//TRIM(cmin)
   WRITE (iomon,'(A)') 'Maximum no of iterations : '//TRIM(cmax(1))
   WRITE (iomon,'(A)') 'Mean no of iterations : '//TRIM(cmean)
ENDIF

CALL log_timer_out(npcc,itm_hydro)


RETURN

9001 FORMAT(A,I3,1X,G15.7)
9002 FORMAT(A,G15.7,1X,'(',A,',',A,')')
9003 FORMAT(A,G15.7,1X,'(',A,',',A,',',A,')')

END SUBROUTINE hydrodynamic_equations

!========================================================================

SUBROUTINE correct_free_surf
!************************************************************************
!
! *correct_free_surf* update free surface correction for the implicit scheme  
!
! Author - Pieter Rauwoens
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.8
!
! Description - 
!
! Reference -
!
! Calling program - hydrodynamic_equations
!
! External calls - define_2dobc_spec, multi_grid_scheme, open_boundary_conds_2d,
!                  update_2dobc_data, water_depths
!
! Module calls - error_alloc, exchange_mod, max_vars
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
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE paral_utilities, ONLY: max_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
INTEGER, SAVE :: maxdatuv, maxdatxy, nofiluv, nofilxy
INTEGER (KIND=kndilong), SAVE, ALLOCATABLE, DIMENSION(:,:) :: nosecsuv, nosecsxy
REAL :: fac
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: obdatuv, obdatxy


procname(pglev+1) = 'correct_free_surf_mg'
CALL log_timer_in(npcc)

!
!1. Initialise open boundary conditions on first call
!----------------------------------------------------
!

IF (nt.EQ.0) THEN

!
!1.1 Normal open boundary conditions
!-----------------------------------
!

   IF (iopt_obc_2D.EQ.1) THEN
!     ---define specifier arrays
      CALL define_2dobc_spec(io_2uvobc)
!     ---allocate arrays
      IF (ALLOCATED(nosecsuv)) DEALLOCATE (nosecsuv,obdatuv)
      nofiluv = maxdatafiles(io_2uvobc,1)
      ALLOCATE (nosecsuv(2:nofiluv,2),STAT=errstat)
      CALL error_alloc('nosecsuv',2,(/nofiluv-1,2/),kndilong)
      IF (nofiluv.GT.1) THEN
         nosecsuv = 0
         maxdatuv = MAXVAL(no2dobuv)
      ELSE
         maxdatuv = 0
      ENDIF
      ALLOCATE (obdatuv(maxdatuv,2,2:nofiluv,2),STAT=errstat)
      CALL error_alloc('obdatuv',4,(/maxdatuv,2,nofiluv-1,2/),kndrtype)
      IF (SIZE(obdatuv).GT.0) obdatuv = 0.0
!     ---update data at initial time
      CALL update_2dobc_data(obdatuv,no2dobuv,maxdatuv,nofiluv,nosecsuv,&
                           & io_2uvobc)
   ENDIF

!
!1.2 Tangential open boundary conditions
!---------------------------------------
!

   IF (iopt_obc_2D_tang.EQ.1) THEN
!     ---define specifier arrays
      CALL define_2dobc_spec(io_2xyobc)
!     ---allocate arrays
      IF (ALLOCATED(nosecsxy)) DEALLOCATE (nosecsxy,obdatxy)
      nofilxy = maxdatafiles(io_2xyobc,1)
      ALLOCATE (nosecsxy(2:nofilxy,2),STAT=errstat)
      CALL error_alloc('nosecsxy',2,(/nofilxy-1,2/),kndilong)
      IF (nofilxy.GT.1) THEN
         nosecsxy = 0
         maxdatxy = MAXVAL(no2dobxy)
      ELSE
         maxdatxy = 0
      ENDIF
      ALLOCATE (obdatxy(maxdatxy,2,2:nofilxy,2),STAT=errstat)
      CALL error_alloc('obdatxy',4,(/maxdatxy,2,nofilxy-1,2/),kndrtype)
      IF (SIZE(obdatxy).GT.0) obdatxy = 0.0 
!     ---update data at initial time
      CALL update_2dobc_data(obdatxy,no2dobxy,maxdatxy,nofilxy,nosecsxy,&
                           & io_2xyobc)
   ENDIF

   GOTO 1000

ENDIF

!
!2. Save values at old time 
!---------------------------
!

IF (iopt_grid_nodim.EQ.3.AND.itsimp.EQ.1) THEN
   udvel_old = udvel; vdvel_old = vdvel
ENDIF

!
!3. Update open boundary data
!----------------------------
!

IF (iopt_obc_2D.EQ.1.AND.itsimp.EQ.1) THEN
   obc2uvatu_old = obc2uvatu
   obc2uvatv_old = obc2uvatv
   CALL update_2dobc_data(obdatuv,no2dobuv,maxdatuv,nofiluv,nosecsuv,io_2uvobc)
ENDIF

!
!4. Depth integrated currents (predictor step)
!---------------------------------------------
!
!4.1 Interior points
!-------------------
!

IF (iopt_grid_nodim.EQ.3) THEN 
!   ---U-nodes
   udvel(1:ncloc,1:nrloc) = SUM(uvel(1:ncloc,1:nrloc,:)*delzatu(1:ncloc,:,:), &
                              & DIM=3,MASK=nodeatu(1:ncloc,1:nrloc,:).GT.1)
!   ---V-nodes
   vdvel(1:ncloc,1:nrloc) = SUM(vvel(1:ncloc,1:nrloc,:)*delzatv(:,1:nrloc,:),&
                              & DIM=3,MASK=nodeatv(1:ncloc,1:nrloc,:).GT.1)
ENDIF

!
!4.2 Zero currents at blocked velocity nodes
!-------------------------------------------
!

IF (iopt_weibar.EQ.1) THEN
   WHERE (node2du(1:ncloc,1:nrloc).EQ.2.AND.&
        & deptotatu(1:ncloc,1:nrloc).LE.dmin_fld)
      udvel(1:ncloc,1:nrloc) = 0.0
   END WHERE
   WHERE (node2dv(1:ncloc,1:nrloc).EQ.2.AND.&
        & deptotatv(1:ncloc,1:nrloc).LE.dmin_fld)
      vdvel(1:ncloc,1:nrloc) = 0.0
   END WHERE
ENDIF

!
!4.3 Exchange halo sections
!--------------------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = 1-nhalo; nhexch = 1
   CALL exchange_mod(udvel,lbounds,nhexch,iarr_udvel,corners=.FALSE.)
   nhexch = (/0,0,1,1/)
   CALL exchange_mod(vdvel,lbounds,nhexch,iarr_vdvel,corners=.FALSE.)
ENDIF

!
!5. Predicted currents
!----------------------
!

WHERE (node2du(0:ncloc+1,1:nrloc).GT.1)
   umpred = udvel(0:ncloc+1,1:nrloc)/deptotatu(0:ncloc+1,1:nrloc)
END WHERE
WHERE (node2dv(1:ncloc,0:nrloc+1).GT.1)
   vmpred = vdvel(1:ncloc,0:nrloc+1)/deptotatv(1:ncloc,0:nrloc+1)
END WHERE

!
!6. Apply myltigrid scheme
!-------------------------
!

CALL multi_grid_scheme

!
!7. Update surface elevation
!---------------------------
!
!7.1 Calculate residual 
!----------------------
!

CALL max_vars(ABS(dzeta(1:ncloc,1:nrloc)),dzetaresid,iarr_dzeta,&
            & mask=maskatc_int, commall = .TRUE.)
 
!
!7.2 Update surface elevation 
!----------------------------
!

WHERE (nodeatc(0:ncloc+1,0:nrloc+1).GT.0)
   zeta = zeta + dzeta
END WHERE

!
!8. Depth integrated currents
!----------------------------
!
!8.1 Internal nodes
!------------------
!
!---U-nodes
fac = theta_sur*delt3d
WHERE (node2du(1:ncloc,1:nrloc).EQ.2)
   udvel(1:ncloc,1:nrloc) = deptotatu(1:ncloc,1:nrloc)*&
                          & (umpred(1:ncloc,1:nrloc)+fac*gaccatu(1:ncloc,:)*&
                          & (dzeta(0:ncloc-1,1:nrloc)-dzeta(1:ncloc,1:nrloc))&
                          & /delxatu(1:ncloc,1:nrloc))
END WHERE

!---V-nodes
WHERE (node2dv(1:ncloc,1:nrloc).EQ.2)
   vdvel(1:ncloc,1:nrloc) = deptotatv(1:ncloc,1:nrloc)*&
                          & (vmpred(1:ncloc,1:nrloc)+fac*gaccatv(:,1:nrloc)*&
                          & (dzeta(1:ncloc,0:nrloc-1)-dzeta(1:ncloc,1:nrloc))&
                          & /delyatv(1:ncloc,1:nrloc))
END WHERE

!
!8.2 Apply open boundary conditions
!----------------------------------
!

IF (iopt_obc_2D.EQ.1) THEN

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = 1-nhalo; nhexch = 1
      CALL exchange_mod(udvel,lbounds,nhexch,iarr_udvel,corners=.FALSE.)
      nhexch = (/0,0,1,1/)
      CALL exchange_mod(vdvel,lbounds,nhexch,iarr_vdvel,corners=.FALSE.)
   ENDIF

!  ---apply conditions
   CALL open_boundary_conds_2d

ENDIF

!
!8.3 Exchange halo sections
!----------------------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = 1-nhalo; nhexch = nhalo
   CALL exchange_mod(udvel,lbounds,nhexch,iarr_udvel)
   CALL exchange_mod(vdvel,lbounds,nhexch,iarr_vdvel)
ENDIF

!
!9. Water depths
!---------------
!

CALL water_depths

!
!10. Depth mean current
!----------------------
!
!---U-nodes
WHERE (node2du(:,0:nrloc+1).GT.1)
   umvel = udvel(:,0:nrloc+1)/deptotatu
END WHERE

!---V-nodes
WHERE (node2dv(0:ncloc+1,:).GT.1)
   vmvel = vdvel(0:ncloc+1,:)/deptotatv
END WHERE

1000 CALL log_timer_out(npcc,itm_3dmode)


RETURN

END SUBROUTINE correct_free_surf

!========================================================================

SUBROUTINE current_corr
!************************************************************************
!
! *current_corr* 3-D current field at corrector step
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - hydrodynamic_equations, initialise_model
!
! External calls - current_corr_1d, define_profobc_spec, open_boundary_conds_3d,
!                  relaxation_at_U, relaxation_at_V, update_nest_data_3d,
!                  update_profobc_data
!
! Module calls - error_alloc, exchange_mod
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
USE relaxation
USE switches
USE syspars
USE timepars
USE wavevars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER ::  k, npcc
INTEGER, SAVE :: maxprofs, nofiles, noprofs
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu, maskatv
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: iprofrlx, iprofobu, iprofobv, &
                                          & itypobu, itypobv, noprofsd
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: indexprof, indexvar
INTEGER (KIND=kndilong), SAVE, ALLOCATABLE, DIMENSION(:,:) :: nosecsprof
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2du1, array2du2, array2dv1, &
                                         & array2dv2,obcvel3d, uvelobu, vvelobv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: profvel


!
!1. 1-D case
!-----------
!

IF (iopt_grid_nodim.EQ.1.AND.nt.GT.0) THEN
   CALL current_corr_1d
   RETURN
ENDIF

procname(pglev+1) = 'current_corr'
CALL log_timer_in(npcc)

!
!2. Initialise open boundary conditions on first call
!----------------------------------------------------
!

IF (nt.EQ.0) THEN

   nofiles = maxdatafiles(io_3uvobc,1)

!  ---deallocate on first call (if needed)
   IF (ALLOCATED(uvelobu)) THEN
      DEALLOCATE (uvelobu,vvelobv,iprofrlx,itypobu,itypobv,iprofobu,iprofobv,&
                & obcvel3d)
   ENDIF
   IF (ALLOCATED(noprofsd)) DEALLOCATE (noprofsd,indexprof,indexvar)
   IF (ALLOCATED(nosecsprof)) DEALLOCATE (nosecsprof,profvel)

!  ---allocate specifier arrays
   ALLOCATE (uvelobu(nobu,nz),STAT=errstat)
   CALL error_alloc('uvelobu',2,(/nobu,nz/),kndrtype)
   IF (nobu.GT.0) uvelobu = float_fill
   ALLOCATE (vvelobv(nobv,nz),STAT=errstat)
   CALL error_alloc('vvelobv',2,(/nobv,nz/),kndrtype)
   IF (nobv.GT.0) vvelobv = float_fill
   ALLOCATE (iprofrlx(norlxzones),STAT=errstat)
   CALL error_alloc('iprofrlx',1,(/norlxzones/),kndint)
   IF (norlxzones.GT.0) iprofrlx = 0
   ALLOCATE (itypobu(nobu),STAT=errstat)
   CALL error_alloc('itypobu',1,(/nobu/),kndint)
   IF (nobu.GT.0) itypobu = 0
   ALLOCATE (itypobv(nobv),STAT=errstat)
   CALL error_alloc('itypobv',1,(/nobv/),kndint)
   IF (nobv.GT.0) itypobv = 0
   ALLOCATE (iprofobu(nobu),STAT=errstat)
   CALL error_alloc('iprofobu',1,(/nobu/),kndint)
   IF (nobu.GT.0) iprofobu = 0
   ALLOCATE (iprofobv(nobv),STAT=errstat)
   CALL error_alloc('iprofobv',1,(/nobv/),kndint)
   IF (nobv.GT.0) iprofobv = 0
   IF (iopt_obc_3D.EQ.1) THEN
      ALLOCATE (noprofsd(2:nofiles),STAT=errstat)
      CALL error_alloc('noprofsd',1,(/nofiles-1/),kndint)
      ALLOCATE (indexprof(nobu+nobv,2:nofiles),STAT=errstat)
      CALL error_alloc('indexprof',2,(/nobu+nobv,nofiles-1/),kndint)
      ALLOCATE (indexvar(nobu+nobv,2:nofiles),STAT=errstat)
      CALL error_alloc('indexvar',2,(/nobu+nobv,nofiles-1/),kndint)
   ENDIF

!  ---define specifier arrays
   IF (iopt_obc_3D.EQ.1) THEN
      CALL define_profobc_spec(io_3uvobc,itypobu,itypobv,iprofobu,iprofobv,&
                             & iprofrlx,noprofsd,indexprof,indexvar,1,nofiles,&
                             & nobu,nobv)
   ENDIF

!  ---allocate o.b. data arrays
   noprofs = 0; maxprofs = 0
   IF (nofiles.GT.1) noprofs = SUM(noprofsd)
   ALLOCATE (obcvel3d(0:noprofs,nz),STAT=errstat)
   CALL error_alloc('obcvel3d',2,(/noprofs+1,nz/),kndrtype)
   obcvel3d = float_fill
   IF (nofiles.GT.1) THEN
      maxprofs = MAXVAL(noprofsd)
      ALLOCATE (nosecsprof(2:nofiles,2),STAT=errstat)
      CALL error_alloc('nosecsprof',2,(/nofiles-1,2/),kndilong)
      ALLOCATE (profvel(maxprofs,nz,2:nofiles,2),STAT=errstat)
      CALL error_alloc('profvel',4,(/maxprofs,nz,nofiles-1,2/),kndrtype)
   ENDIF

!  ---update data at initial time
   IF (nofiles.GT.1) THEN
      CALL update_profobc_data(profvel,obcvel3d,noprofsd,&
                             & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                             & maxprofs,noprofs,1,nofiles,nosecsprof,io_3uvobc)
   ENDIF

!  ---write at nest locations
   IF (iopt_nests.EQ.1) THEN
       CALL update_nest_data_3d(io_3uvnst)
       CALL update_nest_data_3d(io_3xynst)
   ENDIF

   GOTO 1000

ENDIF

!
!3. Allocate arrays
!------------------
!

ALLOCATE (maskatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (array2du1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2du1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2du2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2du2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dv1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dv1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dv2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dv2',2,(/ncloc,nrloc/),kndrtype)

!
!4. Baroclinic current
!---------------------
!
!4.1 Mask arrays
!---------------
!

maskatu = nodeatu(1:ncloc,1:nrloc,:).GT.1
maskatv = nodeatv(1:ncloc,1:nrloc,:).GT.1

!
!4.2 Substract depth-mean part
!-----------------------------
!
!---U-nodes
k_421: DO k=1,nz
   WHERE (maskatu(:,:,k))
      uvel(1:ncloc,1:nrloc,k) = uvel(1:ncloc,1:nrloc,k) - umpred(1:ncloc,:)
   END WHERE
ENDDO k_421

!---V-nodes
k_422: DO k=1,nz
   WHERE (maskatv(:,:,k))
      vvel(1:ncloc,1:nrloc,k) = vvel(1:ncloc,1:nrloc,k) - vmpred(:,1:nrloc)
   END WHERE
ENDDO k_422

!
!4.4. Add correction to baroclinic current
!-----------------------------------------
!
!4.4.1 U-nodes
!-------------
!

WHERE (node2du(1:ncloc,1:nrloc).EQ.2)
   array2du1 = deptotatu_old(1:ncloc,:)/deptotatu(1:ncloc,1:nrloc)
END WHERE
k_441: DO k=1,nz
   WHERE (nodeatu(1:ncloc,1:nrloc,k).EQ.2)
      uvel(1:ncloc,1:nrloc,k) = array2du1*uvel(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_441

!
!4.4.2  V-nodes
!--------------
!

WHERE (node2dv(1:ncloc,1:nrloc).EQ.2)
   array2dv1 = deptotatv_old(:,1:nrloc)/deptotatv(1:ncloc,1:nrloc)
END WHERE
k_442: DO k=1,nz
   WHERE (nodeatv(1:ncloc,1:nrloc,k).EQ.2)
      vvel(1:ncloc,1:nrloc,k) = array2dv1*vvel(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_442

!
!4.5 Exchange halo sections
!--------------------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = (/1-nhalo,1-nhalo,1/); nhexch = (/nhalo,nhalo,0,0/)
   CALL exchange_mod(uvel,lbounds,nhexch,iarr_uvel)
   nhexch = (/0,0,nhalo,nhalo/)
   CALL exchange_mod(vvel,lbounds,nhexch,iarr_vvel)
ENDIF

!
!5. Apply open boundary conditions
!---------------------------------
!
!---update data
IF (itsimp.EQ.1.AND.nofiles.GT.1) THEN
   CALL update_profobc_data(profvel,obcvel3d,noprofsd,&
                          & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                          & maxprofs,noprofs,1,nofiles,nosecsprof,io_3uvobc)
ENDIF

!---update baroclinic current at open boundaries
CALL open_boundary_conds_3d(uvelobu,vvelobv,obcvel3d,io_3uvobc,itypobu,itypobv,&
                          & iprofobu,iprofobv,noprofs,nobu,nobv)

!---apply relaxation scheme in open boundary zones
IF (iopt_obc_relax.EQ.1.AND.inodesrlx(2).GT.0) THEN
   CALL relaxation_at_U(uvel,uvelobu,nz,iprofrlx)
   CALL relaxation_at_V(vvel,vvelobv,nz,iprofrlx)
ENDIF

!
!6. Corrected and filtered currents 
!----------------------------------
!
!6.1 Explicit scheme
!-------------------
!

IF (iopt_hydro_impl.EQ.0) THEN

!
!6.1.1 U-nodes
!--------------
!

   WHERE (node2du(1:ncloc,1:nrloc).GT.1)
      array2du2 = udfvel(1:ncloc,:)/deptotatu(1:ncloc,1:nrloc)
   END WHERE
   k_611: DO k=1,nz
      WHERE (maskatu(:,:,k))
         ufvel(1:ncloc,:,k) = uvel(1:ncloc,1:nrloc,k) + array2du2
         uvel(1:ncloc,1:nrloc,k) = uvel(1:ncloc,1:nrloc,k) + &
                                 & umvel(1:ncloc,1:nrloc)
      END WHERE
   ENDDO k_611

!
!6.1.2 V-nodes
!-------------
!
!  ---interior points
   WHERE (node2dv(1:ncloc,1:nrloc).GT.1)
      array2dv2 = vdfvel(:,1:nrloc)/deptotatv(1:ncloc,1:nrloc)
   END WHERE
   k_612: DO k=1,nz
      WHERE (maskatv(:,:,k))
         vfvel(:,1:nrloc,k) = vvel(1:ncloc,1:nrloc,k) + array2dv2
         vvel(1:ncloc,1:nrloc,k) = vvel(1:ncloc,1:nrloc,k) + &
                                 & vmvel(1:ncloc,1:nrloc)
      END WHERE
   ENDDO k_612

!
!6.2 Implicit scheme
!-------------------
!

ELSEIF (iopt_hydro_impl.EQ.1) THEN

!
!6.2.1 U-nodes
!--------------
!
!  ---interior points
   k_621: DO k=1,nz
      WHERE (maskatu(:,:,k))
         uvel(1:ncloc,1:nrloc,k) = uvel(1:ncloc,1:nrloc,k) + &
                                 & umvel(1:ncloc,1:nrloc)
      END WHERE
   ENDDO k_621

!
!6.2.2 V-nodes
!-------------
!
!  ---interior points
   k_622: DO k=1,nz
      WHERE (maskatv(:,:,k))
         vvel(1:ncloc,1:nrloc,k) = vvel(1:ncloc,1:nrloc,k) + &
                                 & vmvel(1:ncloc,1:nrloc)
      END WHERE
   ENDDO k_622

ENDIF

!
!7. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = (/1-nhalo,1-nhalo,1/); nhexch = nh3vel
   CALL exchange_mod(uvel,lbounds,nhexch,iarr_uvel)
   CALL exchange_mod(vvel,lbounds,nhexch,iarr_vvel)
   IF (iopt_hydro_impl.EQ.0) THEN
      lbounds = (/2-nhalo,1,1/); nhexch = (/nhfvel-1,nhfvel,0,0/)
      CALL exchange_mod(ufvel,lbounds,nhexch,iarr_ufvel)
      lbounds = (/1,2-nhalo,1/); nhexch = (/0,0,nhfvel-1,nhfvel/)
      CALL exchange_mod(vfvel,lbounds,nhexch,iarr_vfvel)
   ENDIF
ENDIF

!
!8. Deallocate arrays
!--------------------
!
!---work space
DEALLOCATE (maskatu,maskatv)
DEALLOCATE (array2du1,array2du2,array2dv1,array2dv2)

1000 CALL log_timer_out(npcc,itm_3dmode)


RETURN

END SUBROUTINE current_corr

!========================================================================

SUBROUTINE current_corr_1d
!************************************************************************
!
! *current_corr_1d* Corrected current field for 1-D case
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.0
!
! Description - 
!
! Reference -
!
! Calling program - current_corr
!
! External calls -
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, ic, ii, j, jc, jj, k, npcc
LOGICAL, DIMENSION(ncloc,nrloc) :: maskatu, maskatv
REAL, DIMENSION(ncloc,nrloc) :: array2du, array2dv


procname(pglev+1) = 'current_corr_1d'
CALL log_timer_in(npcc)

!
!1. Mask arrays
!--------------
!

maskatu = node2du(1:ncloc,1:nrloc).EQ.2
maskatv = node2dv(1:ncloc,1:nrloc).EQ.2

!
!2. 3-D currents at interior points
!----------------------------------
!
!---U-nodes
WHERE (maskatu)
   array2du = deptotatu_old(1:ncloc,:)/deptotatu(1:ncloc,1:nrloc)
END WHERE
k_210: DO k=1,nz
   WHERE (maskatu)
      uvel(1:ncloc,1:nrloc,k) = array2du*uvel(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_210

!---V-nodes
WHERE (maskatv)
   array2dv = deptotatv_old(:,1:nrloc)/deptotatv(1:ncloc,1:nrloc)
END WHERE
k_220: DO k=1,nz
   WHERE (maskatv)
      vvel(1:ncloc,1:nrloc,k) = array2dv*vvel(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_220

!
!3. Open boundaries
!------------------
!
!---U-nodes
ii_310: DO ii=1,nobuloc
   i = iobuloc(ii); j = jobuloc(ii)
   ic = MERGE(i+1,i-1,westobu(ii))
   uvel(i,j,:) = uvel(ic,j,:)
ENDDO ii_310

!---V-nodes
jj_320: DO jj=1,nobvloc
   i = iobvloc(jj); j = jobvloc(jj)
   jc = MERGE(j+1,j-1,soutobv(jj))
   vvel(i,j,:) = vvel(i,jc,:)
ENDDO jj_320

!
!4. Depth-integrated and mean currents
!-------------------------------------
!
!---U-nodes
WHERE (node2du(1:ncloc,1:nrloc).GT.1)
   udvel(1:ncloc,1:nrloc) = SUM(uvel(1:ncloc,1:nrloc,:)*delzatu(1:ncloc,:,:),&
                              & DIM=3)
   umvel(1:ncloc,1:nrloc) = udvel(1:ncloc,1:nrloc)/deptotatu(1:ncloc,1:nrloc)
END WHERE

!---V-nodes
WHERE (node2dv(1:ncloc,1:nrloc).GT.1)
   vdvel(1:ncloc,1:nrloc) = SUM(vvel(1:ncloc,1:nrloc,:)*delzatv(:,1:nrloc,:),&
                              & DIM=3)
   vmvel(1:ncloc,1:nrloc) = vdvel(1:ncloc,1:nrloc)/deptotatv(1:ncloc,1:nrloc)
END WHERE

CALL log_timer_out(npcc,itm_1dmode)


RETURN

END SUBROUTINE current_corr_1d

!========================================================================

SUBROUTINE current_pred
!************************************************************************
!
! *current_pred* 3-D current field at predictor step
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - hydrodynamic_equations
!
! External calls - define_profobc_spec, momentum_discharge_3d,
!                  open_boundary_conds_3d, transport_at_U_3d, transport_at_V_3d,
!                  update_profobc_data, update_1dsur_data, weirs_loss,
!                  wave_sources_momentum, weirs_sink
!
! Module calls - error_alloc, exchange_mod, Uarr_at_V, Varr_at_U
!
!************************************************************************
!
USE currents
USE depths
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE modids
USE obconds
USE physpars
USE relaxation
USE switches
USE syspars
USE tide
USE timepars
USE array_interp, ONLY: Uarr_at_V, Varr_at_U
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL, SAVE :: info = .FALSE.
LOGICAL :: fimp
INTEGER :: i, j, k, npcc
INTEGER, SAVE :: ibcbot, maxprofs, nbcb, nofiles, noprofs
REAL :: ximp
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: iprofrlx, noprofsd
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: indexprof, indexvar
INTEGER (KIND=kndilong), SAVE, ALLOCATABLE, DIMENSION(:,:) :: nosecsprof
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask2du, mask2dv
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2du, array2du1, array2du2, &
                                         & array2dv, array2dv1, array2dv2, &
                                         & obcvel3d, udif, vdif
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: bcbotatu, bcbotatv, sinkatu, &
                                    & sinkatv, sourceatu, sourceatu_prev, &
                                    & sourceatv, sourceatv_prev
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: profvel 

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*sinkatu*      REAL    Sink (weir loss) terms at U-nodes                 [1/s]
!*sinkatu*      REAL    Sink (weir loss) terms at V-nodes                 [1/s]
!*sourceatu*    REAL    Source terms in U-current equation              [m/s^2]
!*sourceatv*    REAL    Source terms in V-current equation              [m/s^2]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'current_pred'
CALL log_timer_in(npcc)

fimp = iopt_hydro_impl.EQ.1.AND.maxitsimp.NE.1

!
!1. Initialise boundary arrays on first call
!-------------------------------------------
!

IF (nt.EQ.0) THEN

   nofiles = maxdatafiles(io_3xyobc,1)
   ibcbot = MERGE(1,2,iopt_bstres_nodim.EQ.2)
   nbcb = MERGE(1,2,iopt_bstres_nodim.EQ.2)

!  ---deallocate on first call (if needed)
   IF (ALLOCATED(itypobx)) THEN
      DEALLOCATE (itypobx,itypoby,iprofrlx,bcbotatu,bcbotatv)
   ENDIF
   IF (ALLOCATED(iprofobx)) THEN
      DEALLOCATE (iprofobx,iprofoby,vvelobx,uveloby,noprofsd,indexprof,&
                & indexvar,obcvel3d)
   ENDIF
   IF (ALLOCATED(nosecsprof)) DEALLOCATE (nosecsprof,profvel)

!  ---specifier arrays
   ALLOCATE (itypobx(nobx),STAT=errstat)
   CALL error_alloc('itypobx',1,(/nobx/),kndint)
   IF (nobx.GT.0) itypobx = 0
   ALLOCATE (itypoby(noby),STAT=errstat)
   CALL error_alloc('itypoby',1,(/noby/),kndint)
   IF (noby.GT.0) itypoby = 0
   ALLOCATE (iprofrlx(norlxzones),STAT=errstat)
   CALL error_alloc('iprofrlx',1,(/norlxzones/),kndint)
   IF (norlxzones.GT.0) iprofrlx = 0
   IF (iopt_obc_3D_tang.EQ.1) THEN
      ALLOCATE (iprofobx(nobx),STAT=errstat)
      CALL error_alloc('iprofobx',1,(/nobx/),kndint)
      IF (nobx.GT.0) iprofobx = 0
      ALLOCATE (iprofoby(noby),STAT=errstat)
      CALL error_alloc('iprofoby',1,(/noby/),kndint)
      IF (noby.GT.0) iprofoby = 0
      ALLOCATE (vvelobx(nobx,nz),STAT=errstat)
      CALL error_alloc('vvelobx',2,(/nobx,nz/),kndrtype)
      IF (nobx.GT.0) vvelobx = float_fill
      ALLOCATE (uveloby(noby,nz),STAT=errstat)
      CALL error_alloc('uveloby',2,(/noby,nz/),kndrtype)
      IF (noby.GT.0) uveloby = float_fill
      ALLOCATE (noprofsd(2:nofiles),STAT=errstat)
      CALL error_alloc('noprofsd',1,(/nofiles-1/),kndint)
      ALLOCATE (indexprof(nobx+noby,2:nofiles),STAT=errstat)
      CALL error_alloc('indexprof',2,(/nobx+noby,nofiles-1/),kndint)
      ALLOCATE (indexvar(nobx+noby,2:nofiles),STAT=errstat)
      CALL error_alloc('indexvar',2,(/nobx+noby,nofiles-1/),kndint)
   ENDIF
   
!  ---define specifier arrays
   IF (iopt_obc_3D_tang.EQ.1) THEN
      CALL define_profobc_spec(io_3xyobc,itypobx,itypoby,iprofobx,iprofoby,&
                             & iprofrlx,noprofsd,indexprof,indexvar,1,nofiles,&
                             & nobx,noby)
   ENDIF

!  ---bottom boundary conditions
   ALLOCATE (bcbotatu(ncloc,nrloc,nbcb),STAT=errstat)
   CALL error_alloc('bcbotatu',3,(/ncloc,nrloc,nbcb/),kndrtype)
   ALLOCATE (bcbotatv(ncloc,nrloc,nbcb),STAT=errstat)
   CALL error_alloc('bcbotatv',3,(/ncloc,nrloc,nbcb/),kndrtype)

!  ---allocate o.b. data arrays
   noprofs = 0; maxprofs = 0
   IF (nofiles.GT.1) noprofs = SUM(noprofsd)
   IF (iopt_obc_3D_tang.EQ.1) THEN
      ALLOCATE (obcvel3d(0:noprofs,nz),STAT=errstat)
      CALL error_alloc('obcvel3d',2,(/noprofs+1,nz/),kndrtype)
      obcvel3d = float_fill
   ENDIF
   IF (nofiles.GT.1) THEN
      maxprofs = MAXVAL(noprofsd)
      ALLOCATE (nosecsprof(2:nofiles,2),STAT=errstat)
      CALL error_alloc('nosecsprof',2,(/nofiles-1,2/),kndilong)
      ALLOCATE (profvel(maxprofs,nz,2:nofiles,2),STAT=errstat)
      CALL error_alloc('profvel',4,(/maxprofs,nz,nofiles-1,2/),kndrtype)
   ENDIF

!  ---update data at initial time
   IF (nofiles.GT.1) THEN
      CALL update_profobc_data(profvel,obcvel3d,noprofsd,&
                             & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                             & maxprofs,noprofs,1,nofiles,nosecsprof,io_3xyobc)
   ENDIF

   GOTO 1000

ENDIF

!
!2. Allocate arrays
!------------------
!

IF (itsimp.EQ.1.AND.fimp) THEN
   IF (ALLOCATED(sinkatu)) THEN
      DEALLOCATE (sinkatu,sinkatv,sourceatu_prev,sourceatv_prev)
   ENDIF
   ALLOCATE (sinkatu(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('sinkatu',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (sinkatv(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('sinkatv',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (sourceatu_prev(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('sourceatu_prev',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (sourceatv_prev(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('sourceatv_prev',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

ALLOCATE (maskatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc,nz/),kndlog)
ALLOCATE (mask2du(ncloc,nrloc),STAT=errstat)
CALL error_alloc('mask2du',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (mask2dv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('mask2dv',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (array2du(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2du1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2du1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2du2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2du2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dv1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dv1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2dv2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dv2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sourceatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('sourceatu',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (sourceatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('sourceatv',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (udif(ncloc+1,0:nrloc),STAT=errstat)
CALL error_alloc('udif',2,(/ncloc+1,nrloc+1/),kndrtype)
ALLOCATE (vdif(0:ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('vdif',2,(/ncloc+1,nrloc+1/),kndrtype)

IF (iopt_hydro_impl.EQ.0.OR.(.NOT.fimp)) THEN
   ALLOCATE (sinkatu(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('sinkatu',3,(/ncloc,nrloc,nz/),kndrtype)
   ALLOCATE (sinkatv(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('sinkatv',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

!
!3. Initialise
!-------------
!
!3.1 Mask and work space arrays
!------------------------------
!
!---U-nodes
maskatu = nodeatu(1:ncloc,1:nrloc,:).EQ.2
mask2du = node2du(1:ncloc,1:nrloc).EQ.2
array2du = 0.0

!---V-nodes
maskatv = nodeatv(1:ncloc,1:nrloc,:).EQ.2
mask2dv = node2dv(1:ncloc,1:nrloc).EQ.2
array2dv = 0.0

!
!3.2 Save values at old time
!---------------------------
!

IF (iopt_hydro_impl.EQ.0.OR.itsimp.EQ.1) THEN
   uvel_old = uvel; vvel_old = vvel
ENDIF

!
!3.3 Initialise arrays for new barotropic cycle
!----------------------------------------------
!

IF (iopt_hydro_impl.EQ.0.AND.iopt_grid_nodim.EQ.3) THEN
!  ---U-nodes
   udevint = 0.0; udfvel = 0.0
!  ---V-nodes
   vdevint = 0.0; vdfvel = 0.0
ENDIF

!
!4. Tangential open boundary conditions
!--------------------------------------
!

IF (iopt_obc_3D_tang.EQ.1) THEN
   IF (itsimp.EQ.1.AND.nofiles.GT.1) THEN
!     ---update data
      CALL update_profobc_data(profvel,obcvel3d,noprofsd,&
                             & indexprof(1:maxprofs,:),indexvar(1:maxprofs,:),&
                             & maxprofs,noprofs,1,nofiles,nosecsprof,io_3xyobc)
!     ---update baroclinic current
      CALL open_boundary_conds_3d(vvelobx,uveloby,obcvel3d,io_3xyobc,itypobx,&
                                & itypoby,iprofobx,iprofoby,noprofs,nobx,noby)
   ENDIF
ENDIF

!
!5. Source terms
!---------------
!
!5.1 Explicit scheme or first iteration
!--------------------------------------
!

IF (itsimp.EQ.1) THEN

!
!5.1.1 Initialise
!----------------
!

   sourceatu = 0.0; sinkatu = 0.0
   sourceatv = 0.0; sinkatv = 0.0   

!
!5.1.2 Atmospheric pressure gradient
!-----------------------------------
!

   IF (iopt_grid_nodim.NE.1.AND.iopt_meteo.EQ.1) THEN

!     ---U-nodes
      WHERE (mask2du)
         array2du = (atmpres(0:ncloc-1,1:nrloc)-atmpres(1:ncloc,1:nrloc))&
                  & /delxatu(1:ncloc,1:nrloc)/density_ref
      END WHERE
!     ---V-nodes
      WHERE (mask2dv)
         array2dv = (atmpres(1:ncloc,0:nrloc-1)-atmpres(1:ncloc,1:nrloc))&
                  & /delyatv(1:ncloc,1:nrloc)/density_ref
      END WHERE

   ENDIF

!
!5.1.3 Astronomical tidal force
!------------------------------
!

   IF (iopt_astro_tide.EQ.1) THEN
      WHERE (mask2du)
         array2du = array2du + fxastro(1:ncloc,:)
      END WHERE
      WHERE (mask2dv)
         array2dv = array2dv + fyastro(:,1:nrloc)
      END WHERE
   ENDIF

!
!5.1.4 Baroclinic pressure gradient
!----------------------------------
!

   IF (iopt_dens_grad.EQ.0) THEN
      k_5141: DO k=1,nz
         WHERE (maskatu(:,:,k))
            sourceatu(:,:,k) = array2du
         END WHERE
         WHERE (maskatv(:,:,k))
            sourceatv(:,:,k) = array2dv
         END WHERE
      END DO k_5141
   ELSE
      k_5142: DO k=1,nz
         WHERE (maskatu(:,:,k))
            sourceatu(:,:,k) = array2du + p3dbcgradatu(:,:,k)
         END WHERE
         WHERE (maskatv(:,:,k))
            sourceatv(:,:,k) = array2dv + p3dbcgradatv(:,:,k)
         END WHERE
      ENDDO k_5142
   ENDIF

!
!5.1.5 Coriolis terms (explicit part)
!------------------------------------
!
!  ---U-nodes
   IF (iopt_fld.EQ.0) THEN
      WHERE (mask2du)
         array2du1 = coriolatu(1:ncloc,1:nrloc)
      END WHERE
   ELSEIF (iopt_fld.GT.0) THEN
      WHERE (mask2du)
         array2du1 = alphatu_fld*coriolatu(1:ncloc,1:nrloc)
      END WHERE
   ENDIF
   
   k_5151: DO k=1,nz
      CALL Varr_at_U(vvel_old(0:ncloc,1:nrloc+1,k),array2du2,4,2,(/0,1,k/),&
                  & (/ncloc,nrloc+1,k/),1,iarr_vvel_old,info)
      WHERE (maskatu(:,:,k))
         sourceatu(:,:,k) = sourceatu(:,:,k) + array2du1*array2du2
      END WHERE
   ENDDO k_5151

!  ---V-nodes
   IF (iopt_fld.EQ.0) THEN
      WHERE (mask2dv)
         array2dv1 = coriolatv(1:ncloc,1:nrloc)
      END WHERE
   ELSEIF (iopt_fld.GT.0) THEN
      WHERE (mask2dv)
         array2dv1 = alphatv_fld*coriolatv(1:ncloc,1:nrloc)
      END WHERE
   ENDIF

   k_5152: DO k=1,nz
      CALL Uarr_at_V(uvel_old(1:ncloc+1,0:nrloc,k),array2dv2,4,2,(/1,0,k/),&
                  & (/ncloc+1,nrloc,k/),1,iarr_uvel_old,info)
      WHERE (maskatv(:,:,k))
         sourceatv(:,:,k) = sourceatv(:,:,k) - array2dv1*array2dv2
      END WHERE
   ENDDO k_5152

!
!5.1.6 Bottom stress
!-------------------
!

   IF (iopt_bstres_nodim.EQ.2) THEN
      WHERE (mask2du)
         bcbotatu(:,:,1) = ubstresatu(1:ncloc,:)
      END WHERE
      WHERE (mask2dv)
         bcbotatv(:,:,1) = vbstresatv(:,1:nrloc)
      END WHERE
   ELSEIF (iopt_bstres_nodim.EQ.3) THEN
      WHERE (mask2du)
         bcbotatu(:,:,1) = zeros2d
         bcbotatu(:,:,2) = bfricatu
      END WHERE
      WHERE (mask2dv)
         bcbotatv(:,:,1) = zeros2d
         bcbotatv(:,:,2) = bfricatv
      END WHERE
   ENDIF

!
!5.1.7 Surface waves
!-------------------
!

   IF (iopt_waves_curr.EQ.1) THEN
      CALL wave_sources_momentum(sourceatu,sourceatv,nz)
   ENDIF

!
!5.1.8 Weirs
!-----------
!

   IF (iopt_weibar.EQ.1) THEN
      CALL weirs_loss
      CALL weirs_sink(sinkatu,sinkatv,nz)
   ENDIF

!
!5.1.9 Discharges
!----------------
!

   IF (modfiles(io_discur,1,1)%defined) THEN
      CALL momentum_discharge_3d(sourceatu,sourceatv)
   ENDIF
           
!
!5.1.10 Store at first iteration
!--------------------------------
!

   IF (fimp) THEN
      sourceatu_prev = sourceatu
      sourceatv_prev = sourceatv
   ENDIF

ENDIF

!
!5.2 Explicit or next iteration
!------------------------------
!
!5.2.1 Retrieve from first iteration
!-----------------------------------
!

IF (fimp) THEN
   sourceatu = sourceatu_prev
   sourceatv = sourceatv_prev
ENDIF

!
!5.2.2 Surface slope
!-------------------
!
!5.2.2.1 1-D case
!----------------
!

IF (iopt_sur_1D.EQ.1) THEN

   CALL update_1dsur_data

   WHERE (mask2du) array2du = -gxslope
   WHERE (mask2dv) array2dv = -gyslope

!
!5.2.2.2 3-D case
!----------------
!

ELSEIF (iopt_grid_nodim.EQ.3) THEN

!  ---U-nodes
   WHERE (mask2du)
      array2du = gaccatu(1:ncloc,:)*&
               & deptotatu(1:ncloc,1:nrloc)/deptotatu_old(1:ncloc,1:nrloc)*&
               & (zeta(0:ncloc-1,1:nrloc)-zeta(1:ncloc,1:nrloc))&
               & /delxatu(1:ncloc,1:nrloc)
   END WHERE
!  ---V-nodes
   WHERE (mask2dv)
      array2dv = gaccatv(:,1:nrloc)*&
               & deptotatv(1:ncloc,1:nrloc)/deptotatv_old(1:ncloc,1:nrloc)*&
               & (zeta(1:ncloc,0:nrloc-1)-zeta(1:ncloc,1:nrloc))&
               & /delyatv(1:ncloc,1:nrloc)
   END WHERE

ENDIF

!
!5.2.2.3 Correction to prevent outflow from "dry cells"
!------------------------------------------------------
!


IF (iopt_grid_nodim.EQ.3.AND.iopt_fld.GT.0) THEN

!  ---U-nodes
   j_52231: DO j=1,nrloc
   i_52231: DO i=1,ncloc
      IF (mask2du(i,j)) THEN
         IF (ALL(deptotatc(i-1:i,j).EQ.dmin_fld)) THEN
            array2du(i,j) = 0.0
         ELSEIF (deptotatc(i,j).EQ.dmin_fld.AND.zeta(i-1,j).LT.zeta(i,j)) THEN
            array2du(i,j) = 0.0
         ELSEIF (deptotatc(i-1,j).EQ.dmin_fld.AND.zeta(i,j).LT.zeta(i-1,j)) THEN
            array2du(i,j) = 0.0
         ENDIF
      ENDIF
   ENDDO i_52231
   ENDDO j_52231

!  ---V-nodes
   j_52232: DO j=1,nrloc
   i_52232: DO i=1,ncloc
      IF (mask2dv(i,j)) THEN
         IF (ALL(deptotatc(i,j-1:j).EQ.dmin_fld)) THEN
            array2dv(i,j) = 0.0
         ELSEIF (deptotatc(i,j).EQ.dmin_fld.AND.zeta(i,j-1).LT.zeta(i,j)) THEN
            array2dv(i,j) = 0.0
         ELSEIF (deptotatc(i,j-1).EQ.dmin_fld.AND.zeta(i,j).LT.zeta(i,j-1)) THEN
            array2dv(i,j) = 0.0
         ENDIF
      ENDIF
   ENDDO i_52232
   ENDDO j_52232

ENDIF

!
!5.2.2.4 Add as source terms
!---------------------------
!

k_5214: DO k=1,nz
   WHERE (mask2du)
      sourceatu(:,:,k) = sourceatu(:,:,k) + array2du
   END WHERE
   WHERE (mask2dv)
      sourceatv(:,:,k) = sourceatv(:,:,k) + array2dv
   END WHERE
END DO k_5214

!
!6. Update predicted horizontal currents
!---------------------------------------
!

CALL transport_at_U_3d(sourceatu,sinkatu,1,ibcbot,1,nbcb,usstresatu(1:ncloc,:),&
                     & bcbotatu)
CALL transport_at_V_3d(sourceatv,sinkatv,1,ibcbot,1,nbcb,vsstresatv(:,1:nrloc),&
                     & bcbotatv)

!
!7. Coriolis terms (implicit part)
!---------------------------------
!

IF (iopt_cor_impl.GT.0) THEN

!
!7.1 Exchange halo sections
!---------------------------
!

   IF (iopt_MPI.EQ.1) THEN
      lbounds = (/1-nhalo,1-nhalo,1/); nhexch = (/0,1,1,0/)
      CALL exchange_mod(uvel,lbounds,nhexch,iarr_uvel)
      nhexch  = (/1,0,0,1/)
      CALL exchange_mod(vvel,lbounds,nhexch,iarr_vvel)
   ENDIF

!
!7.2 Add implicit corrections
!-----------------------------
!

   ximp = theta_cor*delt3d

!  ---work space arrays
   IF (iopt_fld_alpha.EQ.0) THEN
      WHERE (mask2du)
         array2du1 = ximp*coriolatu(1:ncloc,1:nrloc)
      END WHERE
      WHERE (mask2dv)
         array2dv1 = ximp*coriolatv(1:ncloc,1:nrloc)
      END WHERE
   ELSEIF (iopt_fld_alpha.GT.0) THEN
      WHERE (mask2du)
         array2du1 = ximp*alphatu_fld*coriolatu(1:ncloc,1:nrloc)
      END WHERE
      WHERE (mask2dv)
         array2dv1 = ximp*alphatv_fld*coriolatv(1:ncloc,1:nrloc)
      END WHERE
   ENDIF

!  ---update currents
   k_720: DO k=1,nz
      udif = uvel(1:ncloc+1,0:nrloc,k)-uvel_old(1:ncloc+1,0:nrloc,k)
      vdif = vvel(0:ncloc,1:nrloc+1,k)-vvel_old(0:ncloc,1:nrloc+1,k)
      CALL Varr_at_U(vdif,array2du2,4,2,(/0,1,k/),(/ncloc,nrloc+1,k/),1,0,info)
      WHERE (maskatu(:,:,k))
         uvel(1:ncloc,1:nrloc,k) = uvel(1:ncloc,1:nrloc,k) + array2du1*&
          & (array2du2-array2du1*udif(1:ncloc,1:nrloc))/(1.0+array2du1**2)
      END WHERE
      CALL Uarr_at_V(udif,array2dv2,4,2,(/1,0,k/),(/ncloc+1,nrloc,k/),1,0,info)
      WHERE (maskatv(:,:,k))
         vvel(1:ncloc,1:nrloc,k) = vvel(1:ncloc,1:nrloc,k) - array2dv1*&
          & (array2dv2+array2dv1*vdif(1:ncloc,1:nrloc))/(1.0+array2dv1**2)
      END WHERE
   ENDDO k_720

ENDIF

!
!8. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = (/1-nhalo,1-nhalo,1/); nhexch = (/1,1,0,0/)
   CALL exchange_mod(uvel,lbounds,nhexch,iarr_uvel)
   nhexch  = (/0,0,1,1/)
   CALL exchange_mod(vvel,lbounds,nhexch,iarr_vvel)
ENDIF

!
!9. Depth-integrated predicted currents
!--------------------------------------
!

IF (iopt_hydro_impl.EQ.0) THEN

   umpred = SUM(uvel(0:ncloc+1,1:nrloc,:)*delzatu(0:ncloc+1,:,:),DIM=3,&
              & MASK=nodeatu(0:ncloc+1,1:nrloc,:).GT.1)
   WHERE (node2du(0:ncloc+1,1:nrloc).GT.1)
      umpred = umpred/deptotatu(0:ncloc+1,1:nrloc)
   END WHERE
   
   vmpred = SUM(vvel(1:ncloc,0:nrloc+1,:)*delzatv(:,0:nrloc+1,:),DIM=3,&
              & MASK=nodeatv(1:ncloc,0:nrloc+1,:).GT.1)
   WHERE (node2dv(1:ncloc,0:nrloc+1).GT.1)
      vmpred = vmpred/deptotatv(1:ncloc,0:nrloc+1)
   END WHERE

ENDIF
   
!
!10. Deallocate arrays
!---------------------
!

DEALLOCATE (maskatu,mask2du,maskatv,mask2dv)
DEALLOCATE (array2du,array2du1,array2du2,array2dv,array2dv1,array2dv2)
DEALLOCATE (sourceatu,sourceatv,udif,vdif)
IF (iopt_hydro_impl.EQ.0.OR.(.NOT.fimp)) DEALLOCATE (sinkatu,sinkatv)

1000 CALL log_timer_out(npcc,itm_3dmode)


RETURN

END SUBROUTINE current_pred

!========================================================================

SUBROUTINE current_2d
!************************************************************************
!
! *current_2d* Solve 2-D mode equations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - hydrodynamic_equations, initialise_model
!
! External calls - define_2dobc_spec, momentum_discharge_2d,
!                  open_boundary_conds_2d, surface_elevation,
!                  transport_at_U_2d, transport_at_V_2d,
!                  update_2dobc_data, wave_sources_momentum, weirs_loss,
!                  weirs_sink
!
! Module calls - error_alloc, exchange_mod, Uarr_at_V, Varr_at_U
!
!************************************************************************
!
USE currents
USE depths
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE modids
USE obconds
USE physpars
USE switches
USE syspars
USE tide
USE timepars
USE array_interp, ONLY: Uarr_at_V, Varr_at_U, Carr_at_U, Carr_at_V
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: fimp
INTEGER :: i, j, npcc
INTEGER, SAVE :: maxdatuv, maxdatxy, nofiluv, nofilxy
REAL :: ximp
INTEGER, DIMENSION(2) :: lbounds
INTEGER, DIMENSION(4) :: nhexch
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu, maskatv
INTEGER (KIND=kndilong), SAVE, ALLOCATABLE, DIMENSION(:,:) :: nosecsuv, nosecsxy
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: obdatuv, obdatxy

REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2, sinkatu, &
                      & sinkatv, sinkatu_prev, sinkatv_prev, sourceatu, &
                      & sourceatu_prev, sourceatv, sourceatv_prev,udif, vdif

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*sinkatu*      REAL    Sink terms in U-current equation                  [1/s]
!*sinkatv*      REAL    Sink terms in V-current equation                  [1/s]
!*sourceatu*    REAL    Source terms in U-current equation            [m^2/s^2]
!*sourceatv*    REAL    Source terms in V-current equation            [m^2/s^2]
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'current_2d'
CALL log_timer_in(npcc)

fimp = iopt_hydro_impl.EQ.1.AND.maxitsimp.NE.1

!
!1. Initialise open boundary conditions on first call
!----------------------------------------------------
!

IF (nt.EQ.0.AND.iopt_hydro_impl.EQ.0) THEN

!
!1.1 Normal open boundary conditions
!-----------------------------------
!

   IF (iopt_obc_2D.EQ.1) THEN
!     ---define specifier arrays
      CALL define_2dobc_spec(io_2uvobc)
!     ---allocate arrays
      IF (ALLOCATED(nosecsuv)) DEALLOCATE (nosecsuv,obdatuv)
      nofiluv = maxdatafiles(io_2uvobc,1)
      ALLOCATE (nosecsuv(2:nofiluv,2),STAT=errstat)
      CALL error_alloc('nosecsuv',2,(/nofiluv-1,2/),kndilong)
      IF (nofiluv.GT.1) THEN
         nosecsuv = 0
         maxdatuv = MAXVAL(no2dobuv)
      ELSE
         maxdatuv = 0
      ENDIF
      ALLOCATE (obdatuv(maxdatuv,2,2:nofiluv,2),STAT=errstat)
      CALL error_alloc('obdatuv',4,(/maxdatuv,2,nofiluv-1,2/),kndrtype)
      IF (SIZE(obdatuv).GT.0) obdatuv = 0.0
!     ---update data at initial time
      CALL update_2dobc_data(obdatuv,no2dobuv,maxdatuv,nofiluv,nosecsuv,&
                           & io_2uvobc)
   ENDIF

!
!1.2 Tangential open boundary conditions
!---------------------------------------
!

   IF (iopt_obc_2D_tang.EQ.1) THEN
!     ---define specifier arrays
      CALL define_2dobc_spec(io_2xyobc)
!     ---allocate arrays
      IF (ALLOCATED(nosecsxy)) DEALLOCATE (nosecsxy,obdatxy)
      nofilxy = maxdatafiles(io_2xyobc,1)
      ALLOCATE (nosecsxy(2:nofilxy,2),STAT=errstat)
      CALL error_alloc('nosecsxy',2,(/nofilxy-1,2/),kndilong)
      IF (nofilxy.GT.1) THEN
         nosecsxy = 0
         maxdatxy = MAXVAL(no2dobxy)
      ELSE
         maxdatxy = 0
      ENDIF
      ALLOCATE (obdatxy(maxdatxy,2,2:nofilxy,2),STAT=errstat)
      CALL error_alloc('obdatxy',4,(/maxdatxy,2,nofilxy-1,2/),kndrtype)
      IF (SIZE(obdatxy).GT.0) obdatxy = 0.0 
!     ---update data at initial time
      CALL update_2dobc_data(obdatxy,no2dobxy,maxdatxy,nofilxy,nosecsxy,&
                           & io_2xyobc)
   ENDIF

   GOTO 1000

ENDIF

!
!2. Allocate arrays
!------------------
!

IF (itsimp.EQ.1.AND.fimp) THEN
   IF (ALLOCATED(sinkatu_prev)) THEN
      DEALLOCATE (sinkatu_prev,sinkatv_prev,sourceatu_prev,sourceatv_prev)
   ENDIF
   ALLOCATE (sinkatu_prev(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('sinkatu_prev',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (sinkatv_prev(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('sinkatv_prev',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (sourceatu_prev(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('sourceatu_prev',2,(/ncloc,nrloc/),kndrtype)
   ALLOCATE (sourceatv_prev(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('sourceatv_prev',2,(/ncloc,nrloc/),kndrtype)
ENDIF

ALLOCATE (maskatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatu',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatv',2,(/ncloc,nrloc/),kndlog)
ALLOCATE (array2d1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2d2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sinkatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sinkatu',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sinkatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sinkatv',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sourceatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sourceatu',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (sourceatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sourceatv',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (udif(ncloc+1,0:nrloc),STAT=errstat)
CALL error_alloc('udif',2,(/ncloc+1,nrloc+1/),kndrtype)
ALLOCATE (vdif(0:ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('vdif',2,(/ncloc+1,nrloc+1/),kndrtype)

!
!3. Initialise
!-------------
!
!3.1 Mask arrays
!---------------
!

maskatu = node2du(1:ncloc,1:nrloc).EQ.2
maskatv = node2dv(1:ncloc,1:nrloc).EQ.2

!
!3.2 Save values at old time 
!-----------------------------
!

IF (itsimp.EQ.1) THEN
   udvel_old = udvel; vdvel_old = vdvel
   umvel_old = umvel; vmvel_old = vmvel
ENDIF

!
!3.3 Work space arrays
!---------------------
!

udif = 0.0; vdif = 0.0

!
!3.4 Update open boundary data
!-----------------------------
!

IF (iopt_hydro_impl.EQ.0) THEN
   IF (iopt_obc_2D.EQ.1) THEN
      obc2uvatu_old = obc2uvatu
      obc2uvatv_old = obc2uvatv
      CALL update_2dobc_data(obdatuv,no2dobuv,maxdatuv,nofiluv,nosecsuv,&
                           & io_2uvobc)
   ENDIF
   IF (iopt_obc_2D_tang.EQ.1) THEN
      CALL update_2dobc_data(obdatxy,no2dobxy,maxdatxy,nofilxy,nosecsxy,&
                           & io_2xyobc)
   ENDIF
ENDIF

!
!3.5 Initialisation for 2-D case
!-------------------------------
!

IF (iopt_hydro_impl.EQ.0) THEN

!  ---2-D case
   IF (iopt_grid_nodim.EQ.2) THEN
      WHERE (node2du(1:ncloc+1,1:nrloc).GT.0)
         udfvel = udvel(1:ncloc+1,1:nrloc)/REAL(ic3d)
      END WHERE
      WHERE (node2dv(1:ncloc,1:nrloc+1).GT.0)
         vdfvel = vdvel(1:ncloc,1:nrloc+1)/REAL(ic3d)
      END WHERE

!  ---3-D case   
   ELSEIF (iopt_grid_nodim.EQ.3) THEN
      WHERE (node2du(1:ncloc+1,1:nrloc).GT.1)
         udfvel = udfvel + udvel(1:ncloc+1,1:nrloc)/REAL(ic3d)
      END WHERE
      WHERE (node2dv(1:ncloc,1:nrloc+1).GT.1)
         vdfvel = vdfvel + vdvel(1:ncloc,1:nrloc+1)/REAL(ic3d)
      END WHERE
   ENDIF

ENDIF

!
!4. Update surface elevation
!---------------------------
!

IF (iopt_hydro_impl.EQ.0) CALL surface_elevation

!
!5. Source terms
!---------------
!
!5.1 Explicit scheme or first iteration
!--------------------------------------
!

IF (itsimp.EQ.1) THEN

!
!5.1.1 Initialise
!----------------
!

   sourceatu = 0.0; sinkatu = 0.0
   sourceatv = 0.0; sinkatv = 0.0   

!
!5.1.2 Atmospheric pressure gradient
!-----------------------------------
!

   IF (iopt_meteo.EQ.1) THEN
!     ---U-nodes
      WHERE (maskatu)
         sourceatu = deptotatu(1:ncloc,1:nrloc)*&
                  & (atmpres(0:ncloc-1,1:nrloc)-atmpres(1:ncloc,1:nrloc))&
                  & /delxatu(1:ncloc,1:nrloc)/density_ref
      END WHERE
!     ---V-nodes
      WHERE (maskatv)
         sourceatv = deptotatv(1:ncloc,1:nrloc)*&
                  & (atmpres(1:ncloc,0:nrloc-1)-atmpres(1:ncloc,1:nrloc))&
                  & /delyatv(1:ncloc,1:nrloc)/density_ref
      END WHERE
   ENDIF

!
!5.1.3 Astronomical tidal force
!------------------------------
!

   IF (iopt_astro_tide.EQ.1) THEN
      WHERE (maskatu)
         sourceatu = sourceatu + deptotatu(1:ncloc,1:nrloc)*fxastro(1:ncloc,:)
      END WHERE
      WHERE (maskatv)
         sourceatv = sourceatv + deptotatv(1:ncloc,1:nrloc)*fyastro(:,1:nrloc)
      END WHERE
   ENDIF

!
!5.1.4 Baroclinic pressure gradient
!----------------------------------
!

   IF (iopt_grid_nodim.EQ.3.AND.iopt_dens_grad.GT.0) THEN
!     ---U-nodes
      WHERE (maskatu)
         sourceatu = sourceatu + p2dbcgradatu
      END WHERE
!     ---V-nodes
      WHERE (maskatv)
         sourceatv = sourceatv + p2dbcgradatv
      END WHERE
   ENDIF

!
!5.1.5 Coriolis terms (explicit part)
!------------------------------------
!
!  ---U-nodes
   CALL Varr_at_U(vdvel_old(0:ncloc,1:nrloc+1),array2d1,4,2,(/0,1,nz/),&
               & (/ncloc,nrloc+1,nz/),1,iarr_vdvel_old,.TRUE.)
   IF (iopt_fld_alpha.EQ.0) THEN
      WHERE (maskatu)
         sourceatu = sourceatu + coriolatu(1:ncloc,1:nrloc)*array2d1
      END WHERE
   ELSEIF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatu)
         sourceatu = sourceatu + alphatu_fld*coriolatu(1:ncloc,1:nrloc)*array2d1
      END WHERE
   ENDIF

!  ---V-nodes
   CALL Uarr_at_V(udvel_old(1:ncloc+1,0:nrloc),array2d1,4,2,(/1,0,nz/),&
               & (/ncloc+1,nrloc,nz/),1,iarr_udvel_old,.TRUE.)
   IF (iopt_fld_alpha.EQ.0) THEN
      WHERE (maskatv)
         sourceatv = sourceatv - coriolatv(1:ncloc,1:nrloc)*array2d1
      END WHERE
   ELSEIF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatv)
         sourceatv = sourceatv - alphatv_fld*coriolatv(1:ncloc,1:nrloc)*array2d1
      END WHERE
   ENDIF

!
!5.1.6 Stress terms
!------------------
!
!5.1.6.1 Bottom
!--------------
!

   IF (iopt_bstres_form.GT.0) THEN

!     ---U-nodes
      IF (iopt_bstres_nodim.EQ.2) THEN
         WHERE (maskatu) 
            sinkatu = bfricatu/deptotatu(1:ncloc,1:nrloc)
         END WHERE
      ELSEIF (iopt_bstres_nodim.EQ.3) THEN
         IF (iopt_fld_alpha.EQ.0) THEN
            WHERE (maskatu)
               sinkatu = bfricatu/deptotatu_old(1:ncloc,:)
               sourceatu = sourceatu + bfricatu*&
                        & (umpred(1:ncloc,:)-uvel(1:ncloc,1:nrloc,1))
            END WHERE
         ELSEIF (iopt_fld_alpha.GT.0) THEN
            WHERE (maskatu)
               sinkatu = bfricatu/deptotatu_old(1:ncloc,:)
               sourceatu = sourceatu + alphatu_fld*bfricatu*&
                        & (umpred(1:ncloc,:)-uvel(1:ncloc,1:nrloc,1))
            END WHERE
         ENDIF
      ENDIF

!     ---V-nodes
      IF (iopt_bstres_nodim.EQ.2) THEN
         WHERE (maskatv) 
            sinkatv = bfricatv/deptotatv(1:ncloc,1:nrloc)
         END WHERE
      ELSEIF (iopt_bstres_nodim.EQ.3) THEN
         IF (iopt_fld_alpha.EQ.0) THEN
            WHERE (maskatv)
               sinkatv = bfricatv/deptotatv_old(:,1:nrloc)
               sourceatv = sourceatv + bfricatv*&
                        & (vmpred(:,1:nrloc)-vvel(1:ncloc,1:nrloc,1))
            END WHERE
         ELSEIF (iopt_fld_alpha.GT.0) THEN
            WHERE (maskatv)
               sinkatv = bfricatv/deptotatv_old(:,1:nrloc)
               sourceatv = sourceatv + alphatv_fld*bfricatv*&
                        & (vmpred(:,1:nrloc)-vvel(1:ncloc,1:nrloc,1))
            END WHERE
         ENDIF
      ENDIF

   ENDIF

!
!5.1.6.2 Surface
!---------------
!
!  ---U-nodes
   WHERE (maskatu)
      sourceatu = sourceatu + usstresatu(1:ncloc,:)
   END WHERE
!  ---V-nodes
   WHERE (maskatv)
      sourceatv = sourceatv + vsstresatv(:,1:nrloc)
   END WHERE

!
!5.1.7 Surface waves
!-------------------
!

   IF (iopt_waves_curr.EQ.1) THEN
      CALL wave_sources_momentum(sourceatu,sourceatv,1)
   ENDIF

!
!5.1.8 Weirs
!-----------
!

   IF (iopt_weibar.EQ.1) THEN
      IF (.NOT.predstep.OR.iopt_grid_nodim.EQ.2) CALL weirs_loss
      CALL weirs_sink(sinkatu,sinkatv,1)
   ENDIF

!
!5.1.9 Discharges
!----------------
!

   IF (modfiles(io_discur,1,1)%defined) THEN
      CALL momentum_discharge_2d(sourceatu,sourceatv)
   ENDIF
           
!
!5.1.10 Store at first iteration
!-------------------------------
!

   IF (fimp) THEN
      sinkatu_prev = sinkatu; sinkatv_prev = sinkatv
      sourceatu_prev = sourceatu; sourceatv_prev = sourceatv
   ENDIF

ENDIF

!
!5.2 Explicit or next iteration
!------------------------------
!
!5.2.1 Retrieve from first iteration
!-----------------------------------
!

IF (fimp) THEN
   sinkatu = sinkatu_prev; sinkatv = sinkatv_prev
   sourceatu = sourceatu_prev; sourceatv = sourceatv_prev
ENDIF

!
!5.2.2 Surface slope
!-------------------
!
!5.2.2.1 Without correction
!--------------------------
!
!---U-nodes
WHERE (maskatu)
   array2d1 = deptotatu(1:ncloc,1:nrloc)*gaccatu(1:ncloc,:)*&
            & (zeta(0:ncloc-1,1:nrloc)-zeta(1:ncloc,1:nrloc)) &
            & /delxatu(1:ncloc,1:nrloc)
END WHERE

!---V-nodes
WHERE (maskatv)
   array2d2 = deptotatv(1:ncloc,1:nrloc)*gaccatv(:,1:nrloc)*&
             & (zeta(1:ncloc,0:nrloc-1)-zeta(1:ncloc,1:nrloc)) &
             & /delyatv(1:ncloc,1:nrloc)
END WHERE

!
!5.2.2.2 With correction to prevent outflow from "dry" cell
!----------------------------------------------------------
!

IF (iopt_fld.GT.0) THEN

!  ---U-nodes
   j_52231: DO j=1,nrloc
   i_52231: DO i=1,ncloc
      IF (maskatu(i,j)) THEN
         IF (ALL(deptotatc(i-1:i,j).EQ.dmin_fld)) THEN
            array2d1(i,j) = 0.0
         ELSEIF (deptotatc(i,j).EQ.dmin_fld.AND.zeta(i-1,j).LT.zeta(i,j)) THEN
            array2d1(i,j) = 0.0
         ELSEIF (deptotatc(i-1,j).EQ.dmin_fld.AND.zeta(i,j).LT.zeta(i-1,j)) THEN
            array2d1(i,j) = 0.0
         ENDIF
      ENDIF
   ENDDO i_52231
   ENDDO j_52231

!  ---V-nodes
   j_52232: DO j=1,nrloc
   i_52232: DO i=1,ncloc
      IF (maskatv(i,j)) THEN
         IF (ALL(deptotatc(i,j-1:j).EQ.dmin_fld)) THEN
            array2d2(i,j) = 0.0
         ELSEIF (deptotatc(i,j).EQ.dmin_fld.AND.zeta(i,j-1).LT.zeta(i,j)) THEN
            array2d2(i,j) = 0.0
         ELSEIF (deptotatc(i,j-1).EQ.dmin_fld.AND.zeta(i,j).LT.zeta(i,j-1)) THEN
            array2d2(i,j) = 0.0
         ENDIF
      ENDIF
   ENDDO i_52232
   ENDDO j_52232

ENDIF

!
!5.2.2.3 Add as source terms
!---------------------------
!

WHERE (maskatu)
   sourceatu = sourceatu + array2d1
END WHERE
WHERE (maskatv)
   sourceatv = sourceatv + array2d2
END WHERE

!
!6. Update 2-D currents
!----------------------
!

CALL transport_at_U_2d(sourceatu,sinkatu)
CALL transport_at_V_2d(sourceatv,sinkatv)

!
!7. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = 1-nhalo; nhexch = nhalo
   CALL exchange_mod(udvel,lbounds,nhexch,iarr_udvel)
   CALL exchange_mod(vdvel,lbounds,nhexch,iarr_vdvel)
ENDIF

!
!8. Coriolis terms (implicit part)
!---------------------------------
!

IF (iopt_cor_impl.GT.0) THEN
   ximp = theta_cor*delt2d
   udif = udvel(1:ncloc+1,0:nrloc) - udvel_old(1:ncloc+1,0:nrloc)
   vdif = vdvel(0:ncloc,1:nrloc+1) - vdvel_old(0:ncloc,1:nrloc+1)
!  ---U-nodes
   CALL Varr_at_U(vdif,array2d2,4,2,(/0,1,nz/),(/ncloc,nrloc+1,nz/),1,0,.TRUE.)
   IF (iopt_fld_alpha.EQ.0) THEN
      WHERE (maskatu)
         array2d1 = ximp*coriolatu(1:ncloc,1:nrloc)
         udvel(1:ncloc,1:nrloc) = udvel(1:ncloc,1:nrloc) + array2d1*&
                  & (array2d2-array2d1*udif(1:ncloc,1:nrloc))/(1.0+array2d1**2)
      END WHERE
   ELSEIF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatu)
         array2d1 = ximp*alphatu_fld*coriolatu(1:ncloc,1:nrloc)
         udvel(1:ncloc,1:nrloc) = udvel(1:ncloc,1:nrloc) + array2d1*&
                  & (array2d2-array2d1*udif(1:ncloc,1:nrloc))/(1.0+array2d1**2)
      END WHERE
   ENDIF
!  ---V-nodes
   CALL Uarr_at_V(udif,array2d2,4,2,(/1,0,nz/),(/ncloc+1,nrloc,nz/),1,0,.TRUE.)
   IF (iopt_fld_alpha.EQ.0) THEN
      WHERE (maskatv)
         array2d1 = ximp*coriolatv(1:ncloc,1:nrloc)
         vdvel(1:ncloc,1:nrloc) = vdvel(1:ncloc,1:nrloc) - array2d1*&
                  & (array2d2+array2d1*vdif(1:ncloc,1:nrloc))/(1.0+array2d1**2)
      END WHERE
   ELSEIF (iopt_fld_alpha.GT.0) THEN
      WHERE (maskatv)
         array2d1 = ximp*alphatv_fld*coriolatv(1:ncloc,1:nrloc)
         vdvel(1:ncloc,1:nrloc) = vdvel(1:ncloc,1:nrloc) - array2d1*&
                  & (array2d2+array2d1*vdif(1:ncloc,1:nrloc))/(1.0+array2d1**2)
      END WHERE
   ENDIF

ENDIF

!
!9. Zero currents at blocked velocity nodes
!------------------------------------------
!

IF (iopt_hydro_impl.EQ.0.AND.iopt_weibar.EQ.1) THEN
   WHERE (maskatu.AND.deptotatu(1:ncloc,1:nrloc).LE.dmin_fld)
      udvel(1:ncloc,1:nrloc) = 0.0
   END WHERE
   WHERE (maskatv.AND.deptotatv(1:ncloc,1:nrloc).LE.dmin_fld)
      vdvel(1:ncloc,1:nrloc) = 0.0
   END WHERE
ENDIF

!
!10. Open boundary conditions
!----------------------------
!

IF (iopt_hydro_impl.EQ.0.AND.iopt_obc_2D.EQ.1) THEN

!  ---exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      lbounds = 1-nhalo; nhexch = (/1,1,0,0/) 
      CALL exchange_mod(udvel,lbounds,nhexch,iarr_udvel,corners=.FALSE.)
      nhexch = (/0,0,1,1/) 
      CALL exchange_mod(vdvel,lbounds,nhexch,iarr_vdvel,corners=.FALSE.)
   ENDIF

!  ---apply conditions
   CALL open_boundary_conds_2d

ENDIF

!
!11. Exchange halo sections
!--------------------------
!

IF (iopt_hydro_impl.EQ.0.AND.iopt_MPI.EQ.1) THEN
   lbounds = 1-nhalo; nhexch = nhalo
   CALL exchange_mod(udvel,lbounds,nhexch,iarr_udvel)
   CALL exchange_mod(vdvel,lbounds,nhexch,iarr_vdvel)
ENDIF

!
!12. Depth-mean currents
!-----------------------
!

IF (iopt_hydro_impl.EQ.0) THEN
!  ---U-nodes
   WHERE (node2du(:,0:nrloc+1).GT.1)
      umvel = udvel(:,0:nrloc+1)/deptotatu
   END WHERE
!  ---V-nodes
   WHERE (node2dv(0:ncloc+1,:).GT.1)
      vmvel = vdvel(0:ncloc+1,:)/deptotatv
   END WHERE
ENDIF

!
!13. Deallocate arrays
!---------------------
!

DEALLOCATE (maskatu,maskatv)
DEALLOCATE (array2d1,array2d2,sinkatu,sinkatv,sourceatu,sourceatv,udif,vdif)

1000 CALL log_timer_out(npcc,itm_2dmode)


RETURN

END SUBROUTINE current_2d

!========================================================================

SUBROUTINE physical_vertical_current
!************************************************************************
!
! *physical_vertical_current* Update physical vertical current
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - hydrodynamic_equations
!
! External calls - 
!
! Module calls - error_alloc, Zcoord_arr
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE switches
USE syspars
USE timepars
USe wavevars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: Zcoord_arr
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: k
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, array2du, array2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dw, denom, zcoord


procname(pglev+1) = 'physical_vertical_current'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatu(ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc+1,nrloc,nz/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc+1,nz/),kndlog)
ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2du(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (array2dv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (array3dw(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('array3dw',3,(/ncloc,nrloc,nz+1/),kndrtype)
ALLOCATE (zcoord(ncloc+1,nrloc+1,nz),STAT=errstat)
CALL error_alloc('zcoord',3,(/ncloc+1,nrloc+1,nz/),kndrtype)
IF (.NOT.(dXregY.AND.dYregX)) THEN
   ALLOCATE (denom(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('denom',3,(/ncloc,nrloc,nz/),kndrtype)
ENDIF

!
!2. Initialise
!-------------
!
!---mask arrays
maskatu = nodeatu(1:ncloc+1,1:nrloc,:).GT.1
maskatv = nodeatv(1:ncloc,1:nrloc+1,:).GT.1

!---zero values at land boundaries
array2du = 0.0; array2dv = 0.0

!---denominator
IF (.NOT.(dXregY.AND.dYregX)) THEN
   k_220: DO k=1,nz
      WHERE (maskatc_int)
         denom(:,:,k) = garea*delzatc(1:ncloc,1:nrloc,k)
      END WHERE
   ENDDO k_220
ENDIF

!
!3. Time derivative
!------------------
!

CALL Zcoord_arr(zcoord(1:ncloc,1:nrloc,:),(/1,1,1/),(/ncloc,nrloc,nz/),'C  ',&
             & .FALSE.)

k_310: DO k=1,nz
   WHERE (maskatc_int)
      wphys(:,:,k) = (deptotatc(1:ncloc,1:nrloc)*zcoord(1:ncloc,1:nrloc,k)&
               & -deptotatc_old(1:ncloc,1:nrloc)*(gscoordatc(1:ncloc,1:nrloc,k)&
               & *deptotatc_old(1:ncloc,1:nrloc)-depmeanatc(1:ncloc,1:nrloc)))&
               & /(0.5*delt3d*(deptotatc(1:ncloc,1:nrloc)&
               & +deptotatc_old(1:ncloc,1:nrloc)))
   END WHERE
ENDDO k_310

!
!4. X-derivative terms
!---------------------
!

CALL Zcoord_arr(zcoord(:,1:nrloc,:),(/1,1,1/),(/ncloc+1,nrloc,nz/),'U  ',&
             & .FALSE.)

!
!4.1 Regular horizontal grid
!---------------------------
!

IF (dYregX) THEN

   k_410: DO k=1,nz

      WHERE (maskatu(:,:,k))
         array2du = delzatu(1:ncloc+1,:,k)*zcoord(:,1:nrloc,k)*&
                  & ufvel(1:ncloc+1,1:nrloc,k)
      END WHERE

      WHERE (maskatc_int)
         wphys(:,:,k) = wphys(:,:,k) + &
                      & (array2du(2:ncloc+1,:)-array2du(1:ncloc,:))&
                      & /(delxatc(1:ncloc,1:nrloc)*delzatc(1:ncloc,1:nrloc,k))
      END WHERE
   ENDDO k_410

!
!4.2 Irregular horizontal grid
!-----------------------------
!

ELSE

   k_420: DO k=1,nz

      WHERE (maskatu(:,:,k))
         array2du = delyatu(1:ncloc+1,1:nrloc)*delzatu(1:ncloc+1,:,k)*&
                  & zcoord(:,1:nrloc,k)*ufvel(1:ncloc+1,1:nrloc,k)
      END WHERE

      WHERE (maskatc_int)
         wphys(:,:,k) = wphys(:,:,k) + (array2du(2:ncloc+1,:)-&
                                      & array2du(1:ncloc,:))/denom(:,:,k)
      END WHERE

   ENDDO k_420

ENDIF

!
!5. Y-derivative terms
!---------------------
!

CALL Zcoord_arr(zcoord(1:ncloc,:,:),(/1,1,1/),(/ncloc,nrloc+1,nz/),'V  ',&
             & .FALSE.)

!
!5.1 Regular horizontal grid
!---------------------------
!

IF (dXregY) THEN

   k_510: DO k=1,nz

      WHERE (maskatv(:,:,k))
         array2dv = delzatv(:,1:nrloc+1,k)*zcoord(1:ncloc,:,k)*&
                  & vfvel(1:ncloc,1:nrloc+1,k)
      END WHERE

      WHERE (maskatc_int)
         wphys(:,:,k) = wphys(:,:,k) + &
                      & (array2dv(:,2:nrloc+1)-array2dv(:,1:nrloc))&
                      & /(delyatc(1:ncloc,1:nrloc)*delzatc(1:ncloc,1:nrloc,k))
      END WHERE

   ENDDO k_510

!
!5.2 Irregular horizontal grid
!-----------------------------
!

ELSE

   k_520: DO k=1,nz

      WHERE (maskatv(:,:,k))
         array2dv = delxatv(1:ncloc,1:nrloc+1)*delzatv(:,1:nrloc+1,k)*&
                  & zcoord(1:ncloc,:,k)*vfvel(1:ncloc,1:nrloc+1,k)
      END WHERE

      WHERE (maskatc_int)
         wphys(:,:,k) = wphys(:,:,k) + (array2dv(:,2:nrloc+1)-&
                                      & array2dv(:,1:nrloc))/denom(:,:,k)
      END WHERE

   ENDDO k_520

ENDIF

!
!6. Z-derivative term
!--------------------
!

WHERE (maskatc_int)
   array2dc = depmeanatc(1:ncloc,1:nrloc)/deptotatc(1:ncloc,1:nrloc)
   array3dw(:,:,1) = 0.0
   array3dw(:,:,nz+1) = 0.0
END WHERE
k_611: DO k=2,nz
   WHERE (maskatc_int)
      array3dw(:,:,k) = (gscoordatw(1:ncloc,1:nrloc,k)-array2dc)*&
                       & wvel(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_611
k_612: DO k=1,nz
   WHERE (maskatc_int)
      wphys(:,:,k) = wphys(:,:,k) + deptotatc(1:ncloc,1:nrloc)*&
               & (array3dw(:,:,k+1)-array3dw(:,:,k))/delzatc(1:ncloc,1:nrloc,k)
   END WHERE
ENDDO k_612

!
!7. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatu,maskatv)
DEALLOCATE (array2dc,array2du,array2dv,array3dw,zcoord)
IF (ALLOCATED(denom)) DEALLOCATE (denom)

CALL log_timer_out()


RETURN

END SUBROUTINE physical_vertical_current

!========================================================================

SUBROUTINE shear_frequency
!************************************************************************
!
! *shear_frequency* Squared shear frequency
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.11.2
!
! Description - 
!
! Calling program - vertical_diff_coefs
!
! Module calls - error_alloc, Uarr_at_C, Varr_at_C
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE syspars
USE turbulence
USE array_interp, ONLY: Uarr_at_C, Varr_at_C
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
INTEGER :: k
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: uvelatc, vvelatc

!
! Name        Type  Purpose
!------------------------------------------------------------------------------
!*uvelatc*    REAL  X-current at C-nodes                                  [m/s]
!*vvelatc*    REAL  Y-current at C-nodes                                  [m/s]
!
!************************************************************************
!


procname(pglev+1) = 'shear_frequency'
CALL log_timer_in()
  
!
!1. Allocate arrays
!------------------
!

ALLOCATE (uvelatc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('uvelatc',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (vvelatc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('vvelatc',3,(/ncloc,nrloc,nz/),kndrtype)

!
!2. Interpolate arrays
!---------------------
!

CALL Uarr_at_C(uvel(1:ncloc+1,1:nrloc,:),uvelatc,1,1,(/1,1,1/),&
                 & (/ncloc+1,nrloc,nz/),1,iarr_uvel,.TRUE.)
CALL Varr_at_C(vvel(1:ncloc,1:nrloc+1,:),vvelatc,1,1,(/1,1,1/),&
            & (/ncloc,nrloc+1,nz/),1,iarr_vvel,.TRUE.)

!
!3. Calculate shear frequency
!----------------------------
!

k_310: DO k=2,nz
   WHERE (maskatc_int)
      shearfreq2(:,:,k) = ((uvelatc(:,:,k)-uvelatc(:,:,k-1))**2 + &
                         & (vvelatc(:,:,k)-vvelatc(:,:,k-1))**2)&
                        & /(delzatw(1:ncloc,1:nrloc,k))**2
   END WHERE
ENDDO k_310

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (uvelatc,vvelatc)

CALL log_timer_out()


RETURN

END SUBROUTINE shear_frequency

!========================================================================

SUBROUTINE surface_elevation
!************************************************************************
!
! *surface_elevation* Solve 2-D continuity equation for surface elevation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - current_2d
!
! External calls - surface_discharge, water_depths
!
! Module calls - error_alloc, exchange_mod
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE meteo
USE modids
USE physpars
USE switches
USE syspars
USE timepars
USE wavevars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER, DIMENSION(4) :: nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2du, array2dv


procname(pglev+1) = 'surface_elevation'
CALL log_timer_in()

!
!1. Initialise and allocate arrays
!---------------------------------
!

IF (.NOT.dYregX) THEN
   ALLOCATE (array2du(ncloc+1,nrloc),STAT=errstat)
   CALL error_alloc('array2du',2,(/ncloc+1,nrloc/),kndrtype)
ENDIF
IF (.NOT.dXregY) THEN
   ALLOCATE (array2dv(ncloc,nrloc+1),STAT=errstat)
   CALL error_alloc('array2dv',2,(/ncloc,nrloc+1/),kndrtype)
ENDIF

!
!2. Update surface elevation 
!---------------------------
!
!2.1 Without Stokes transport
!----------------------------
!
!---X-direction
IF (dYregX) THEN
   WHERE (maskatc_int)
      zeta(1:ncloc,1:nrloc) = zeta(1:ncloc,1:nrloc) - delt2d*&
          & (udvel(2:ncloc+1,1:nrloc)-udvel(1:ncloc,1:nrloc))&
          & /delxatc(1:ncloc,1:nrloc)
   END WHERE
ELSE
   array2du = delyatu(1:ncloc+1,1:nrloc)*udvel(1:ncloc+1,1:nrloc)
   WHERE (maskatc_int)
      zeta(1:ncloc,1:nrloc) = zeta(1:ncloc,1:nrloc) - delt2d*&
          & (array2du(2:ncloc+1,:)-array2du(1:ncloc,:))/garea
   END WHERE
ENDIF

!---Y-direction
IF (dXregY) THEN
   WHERE (maskatc_int)
      zeta(1:ncloc,1:nrloc) = zeta(1:ncloc,1:nrloc) - delt2d*&
          & (vdvel(1:ncloc,2:nrloc+1)-vdvel(1:ncloc,1:nrloc))&
          & /delyatc(1:ncloc,1:nrloc)
   END WHERE
ELSE
   array2dv = delxatv(1:ncloc,1:nrloc+1)*vdvel(1:ncloc,1:nrloc+1)
   WHERE (maskatc_int)
      zeta(1:ncloc,1:nrloc) = zeta(1:ncloc,1:nrloc) - delt2d*&
          & (array2dv(:,2:nrloc+1)-array2dv(:,1:nrloc))/garea
   END WHERE
ENDIF

!
!2.2 Add Stokes transport
!------------------------
!

IF (iopt_waves_curr.EQ.1) THEN

!  ---X-direction
   IF (dYregX) THEN
      WHERE (maskatc_int)
         zeta(1:ncloc,1:nrloc) = zeta(1:ncloc,1:nrloc) - delt2d*&
             & (udstokesatu(2:ncloc+1,:)-udstokesatu(1:ncloc,:))&
             & /delxatc(1:ncloc,1:nrloc)
      END WHERE
   ELSE
      array2du = delyatu(1:ncloc+1,1:nrloc)*udstokesatu(1:ncloc+1,:)
      WHERE (maskatc_int)
         zeta(1:ncloc,1:nrloc) = zeta(1:ncloc,1:nrloc) - delt2d*&
             & (array2du(2:ncloc+1,:)-array2du(1:ncloc,:))/garea
      END WHERE
   ENDIF

!  ---Y-direction
   IF (dXregY) THEN
      WHERE (maskatc_int)
         zeta(1:ncloc,1:nrloc) = zeta(1:ncloc,1:nrloc) - delt2d*&
             & (vdstokesatv(:,2:nrloc+1)-vdstokesatv(:,1:nrloc))&
             & /delyatc(1:ncloc,1:nrloc)
      END WHERE
   ELSE
      array2dv = delxatv(1:ncloc,1:nrloc+1)*vdstokesatv(:,1:nrloc+1)
      WHERE (maskatc_int)
         zeta(1:ncloc,1:nrloc) = zeta(1:ncloc,1:nrloc) - delt2d*&
             & (array2dv(:,2:nrloc+1)-array2dv(:,1:nrloc))/garea
      END WHERE
   ENDIF

ENDIF

!
!3. Add discharges
!-----------------
!

IF (modfiles(io_disvol,1,1)%defined) CALL surface_discharge

!
!4. Add precipitation minus evaporation
!--------------------------------------
!

IF (iopt_sflux_precip.GT.1) THEN
   IF (iopt_meteo_precip.EQ.1) THEN
      zeta(1:ncloc,1:nrloc) = rmaskatc*(zeta(1:ncloc,1:nrloc)+&
                            & delt2d*(evaporation-precipitation)/density_ref)
   ELSEIF (iopt_meteo_precip.EQ.2) THEN
      zeta(1:ncloc,1:nrloc) = rmaskatc*(zeta(1:ncloc,1:nrloc)+&
                            & delt2d*evapminprec/density_ref)
   ENDIF
ENDIF
   
!
!4. Exchange halo sections
!-------------------------
!

IF (iopt_MPI.EQ.1) THEN
   nhexch = 1
   CALL exchange_mod(zeta,(/0,0/),nhexch,iarr_zeta)
ENDIF

!
!5. Update total water depths
!----------------------------
!

CALL water_depths

!
!6. Deallocate arrays
!--------------------
!

IF (ALLOCATED(array2du)) DEALLOCATE (array2du)
IF (ALLOCATED(array2dv)) DEALLOCATE (array2dv)

CALL log_timer_out()


RETURN

END SUBROUTINE surface_elevation

!========================================================================

SUBROUTINE transf_vertical_current
!************************************************************************
!
! *transf_vertical_current* Solve 3-D continuity equation for (transformed)
!                           vertical current
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - hydrodynamic_equations
!
! External calls - 
!
! Module calls - error_alloc, Uarr_at_C, Varr_at_C
!
!************************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE switches
USE syspars
USE timepars
USE wavevars
USE array_interp, ONLY: Uarr_at_C, Varr_at_C
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, k, npcc
REAL :: dsig
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: maskatu, maskatv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2dc, array2du, array2dv
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3dc


procname(pglev+1) = 'transf_vertical_current'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (maskatu(ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('maskatu',3,(/ncloc+1,nrloc,nz/),kndlog)
ALLOCATE (maskatv(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('maskatv',3,(/ncloc,nrloc+1,nz/),kndlog)
ALLOCATE (array2dc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2dc',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2du(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('array2du',2,(/ncloc+1,nrloc/),kndrtype)
ALLOCATE (array2dv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('array2dv',2,(/ncloc,nrloc+1/),kndrtype)
ALLOCATE (array3dc(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('array3dc',3,(/ncloc,nrloc,nz/),kndrtype)

!
!2. Initialise
!-------------
!
!---mask arrays
maskatu = nodeatu(1:ncloc+1,1:nrloc,:).GT.1
maskatv = nodeatv(1:ncloc,1:nrloc+1,:).GT.1

!---zero values at land boundaries
array2du = 0.0; array2dv = 0.0
k_210: DO k=1,nz
   WHERE (.NOT.maskatc_int)
      array3dc(:,:,k) = 0.0
   END WHERE
END DO k_210

!
!3. Eulerian vertical current
!----------------------------
!
!3.1 Terms in X-direction
!------------------------
!
!---regular grid
IF (dYregX) THEN
   k_311: DO k=1,nz
      WHERE (maskatu(:,:,k))
         array2du = delzatu(1:ncloc+1,:,k)*&
                 & (uvel(1:ncloc+1,1:nrloc,k)-umvel(1:ncloc+1,1:nrloc))
      END WHERE
      WHERE (maskatc_int)
         array3dc(:,:,k) = (array2du(2:ncloc+1,:)-array2du(1:ncloc,:))&
                         & /delxatc(1:ncloc,1:nrloc)
      END WHERE
   ENDDO k_311

!---irregular grid
ELSE
   k_312: DO k=1,nz
      WHERE (maskatu(:,:,k))
         array2du = delzatu(1:ncloc+1,:,k)*delyatu(1:ncloc+1,1:nrloc)*&
                 & (uvel(1:ncloc+1,1:nrloc,k)-umvel(1:ncloc+1,1:nrloc))
      END WHERE
      WHERE (maskatc_int)
         array3dc(:,:,k) = (array2du(2:ncloc+1,:)-array2du(1:ncloc,:))/garea
      END WHERE
   ENDDO k_312
ENDIF

!---non-uniform sigma grid
IF (iopt_grid_vtype.EQ.3) THEN
   CALL Uarr_at_C(udfvel,array2dc,1,1,(/1,1,nz/),(/ncloc+1,nrloc,nz/),1,&
                & iarr_udvel,.TRUE.)
   j_313: DO j=1,nrloc
   i_313: DO i=1,ncloc
      IF (ALL(node2du(i:i+1,j).GT.1)) THEN
         k_3131: DO k=1,nz
            dsig = gscoordatuw(i+1,j,k+1) - gscoordatuw(i+1,j,k) - &
                   gscoordatuw(i,j,k+1) + gscoordatuw(i,j,k)
            array3dc(i,j,k) = array3dc(i,j,k) + array2dc(i,j)*dsig/delxatc(i,j)
         ENDDO k_3131
      ENDIF
   ENDDO i_313
   ENDDO j_313
ENDIF

!
!3.2 Terms in Y-direction
!------------------------
!
!---regular grid
IF (dXregY) THEN
   k_321: DO k=1,nz
      WHERE (maskatv(:,:,k))
         array2dv = delzatv(:,1:nrloc+1,k)*&
                  & (vvel(1:ncloc,1:nrloc+1,k)-vmvel(1:ncloc,1:nrloc+1))
      END WHERE
      WHERE (maskatc_int)
         array3dc(:,:,k) = array3dc(:,:,k) + &
                         & (array2dv(:,2:nrloc+1)-array2dv(:,1:nrloc))&
                         & /delyatc(1:ncloc,1:nrloc)
      END WHERE
   ENDDO k_321

!---irregular grid
ELSE
   k_322: DO k=1,nz
      WHERE (maskatv(:,:,k))
         array2dv = delzatv(:,1:nrloc+1,k)*delxatv(1:ncloc,1:nrloc+1)*&
                  & (vvel(1:ncloc,1:nrloc+1,k)-vmvel(1:ncloc,1:nrloc+1))
      END WHERE
      WHERE (maskatc_int)
         array3dc(:,:,k) = array3dc(:,:,k) + &
                         & (array2dv(:,2:nrloc+1)-array2dv(:,1:nrloc))/garea
      END WHERE
   ENDDO k_322
ENDIF

!---non-uniform sigma grid
IF (iopt_grid_vtype.EQ.3) THEN
   CALL Varr_at_C(vdfvel,array2dc,1,1,(/1,1,nz/),(/ncloc,nrloc+1,nz/),1,&
                & iarr_vdvel,.TRUE.)
   j_323: DO j=1,nrloc
   i_323: DO i=1,ncloc
      IF (ALL(node2dv(i,j:j+1).GT.1)) THEN
         k_3231: DO k=1,nz
            dsig = gscoordatvw(i,j+1,k+1) - gscoordatvw(i,j+1,k) - &
                   gscoordatvw(i,j,k+1) + gscoordatvw(i,j,k)
            array3dc(i,j,k) = array3dc(i,j,k) + array2dc(i,j)*dsig/delyatc(i,j)
         ENDDO k_3231
      ENDIF
   ENDDO i_323
   ENDDO j_323
ENDIF

!
!3.3 Substract depth-integrated value for numerical accuracy
!----------------------------------------------------------
!

WHERE (maskatc_int)
   array2dc = SUM(array3dc,DIM=3)/REAL(nz)
END WHERE
k_330: DO k=1,nz
   WHERE (maskatc_int)
      array3dc(:,:,k) = array3dc(:,:,k) - array2dc
   END WHERE
ENDDO k_330

!
!3.4 Vertical (transformed) current
!---------------------------------
!

WHERE (maskatc_int) wvel(1:ncloc,1:nrloc,1) = 0.0
k_340: DO k=1,nz
   WHERE (maskatc_int)
      wvel(1:ncloc,1:nrloc,k+1) = wvel(1:ncloc,1:nrloc,k) - array3dc(:,:,k)
   END WHERE
ENDDO k_340
WHERE (maskatc_int) wvel(1:ncloc,1:nrloc,nz+1) = 0.0

!
!4. Vertical Stokes current
!--------------------------
!

IF (iopt_waves_curr.EQ.1) THEN

!
!4.1 Terms in X-direction
!------------------------
!
!  ---regular grid
   IF (dYregX) THEN
      k_411: DO k=1,nz
         WHERE (maskatu(:,:,k))
            array2du = delzatu(1:ncloc+1,:,k)*&
              & (ustokesatu(1:ncloc+1,1:nrloc,k)-umstokesatu(1:ncloc+1,1:nrloc))
         END WHERE
         WHERE (maskatc_int)
            array3dc(:,:,k) = (array2du(2:ncloc+1,:)-array2du(1:ncloc,:))&
                            & /delxatc(1:ncloc,1:nrloc)
         END WHERE
      ENDDO k_411

!  ---irregular grid
   ELSE
      k_412: DO k=1,nz
         WHERE (maskatu(:,:,k))
            array2du = delzatu(1:ncloc+1,:,k)*delyatu(1:ncloc+1,1:nrloc)*&
              & (ustokesatu(1:ncloc+1,1:nrloc,k)-umstokesatu(1:ncloc+1,1:nrloc))
         END WHERE
         WHERE (maskatc_int)
            array3dc(:,:,k) = (array2du(2:ncloc+1,:)-array2du(1:ncloc,:))/garea
         END WHERE
      ENDDO k_412
   ENDIF

!  ---non-uniform sigma grid
   IF (iopt_grid_vtype.EQ.3) THEN
      CALL Uarr_at_C(udstokesatu(1:ncloc+1,1:nrloc),array2dc,1,1,(/1,1,nz/),&
                  & (/ncloc+1,nrloc,nz/),1,iarr_udvel,.TRUE.)
      j_413: DO j=1,nrloc
      i_413: DO i=1,ncloc
         IF (ALL(node2du(i:i+1,j).GT.1)) THEN
            k_4131: DO k=1,nz
               dsig = gscoordatuw(i+1,j,k+1) - gscoordatuw(i+1,j,k) - &
                      gscoordatuw(i,j,k+1) + gscoordatuw(i,j,k)
               array3dc(i,j,k) = array3dc(i,j,k) + &
                               & array2dc(i,j)*dsig/delxatc(i,j)
            ENDDO k_4131
         ENDIF
      ENDDO i_413
      ENDDO j_413
   ENDIF

!
!4.2 Terms in Y-direction
!------------------------
!
!  ---regular grid
   IF (dXregY) THEN
      k_421: DO k=1,nz
         WHERE (maskatv(:,:,k))
            array2dv = delzatv(:,1:nrloc+1,k)*&
              & (vstokesatv(1:ncloc,1:nrloc+1,k)-vmstokesatv(1:ncloc,1:nrloc+1))
         END WHERE
         WHERE (maskatc_int)
            array3dc(:,:,k) = array3dc(:,:,k) + &
                            & (array2dv(:,2:nrloc+1)-array2dv(:,1:nrloc))&
                            & /delyatc(1:ncloc,1:nrloc)
         END WHERE
      ENDDO k_421

!  ---irregular grid
   ELSE
      k_422: DO k=1,nz
         WHERE (maskatv(:,:,k))
            array2dv = delzatv(:,1:nrloc+1,k)*delxatv(1:ncloc,1:nrloc+1)*&
              & (vstokesatv(1:ncloc,1:nrloc+1,k)-vmstokesatv(1:ncloc,1:nrloc+1))
         END WHERE
         WHERE (maskatc_int)
            array3dc(:,:,k) = array3dc(:,:,k) + &
                            & (array2dv(:,2:nrloc+1)-array2dv(:,1:nrloc))/garea
         END WHERE
      ENDDO k_422
   ENDIF

!  ---non-uniform sigma grid
   IF (iopt_grid_vtype.EQ.3) THEN
      CALL Varr_at_C(vdstokesatv(1:ncloc,1:nrloc+1),array2dc,1,1,(/1,1,nz/),&
                  & (/ncloc,nrloc+1,nz/),1,iarr_vdvel,.TRUE.)
      j_423: DO j=1,nrloc
      i_423: DO i=1,ncloc
         IF (ALL(node2dv(i,j:j+1).GT.1)) THEN
            k_4231: DO k=1,nz
               dsig = gscoordatvw(i,j+1,k+1) - gscoordatvw(i,j+1,k) - &
                      gscoordatvw(i,j,k+1) + gscoordatvw(i,j,k)
               array3dc(i,j,k) = array3dc(i,j,k) + &
                               & array2dc(i,j)*dsig/delyatc(i,j)
            ENDDO k_4231
         ENDIF
      ENDDO i_423
      ENDDO j_423
   ENDIF

!
!4.3 Substract depth-integrated value for numerical accuracy
!----------------------------------------------------------
!

   WHERE (maskatc_int)
      array2dc = SUM(array3dc,DIM=3)/REAL(nz)
   END WHERE
   k_430: DO k=1,nz
      WHERE (maskatc_int)
         array3dc(:,:,k) = array3dc(:,:,k) - array2dc
      END WHERE
   ENDDO k_430

!
!4.4 Vertical (transformed) current
!---------------------------------
!

   WHERE (maskatc_int) wstokesatw(1:ncloc,1:nrloc,1) = 0.0
   k_440: DO k=1,nz
      WHERE (maskatc_int)
         wstokesatw(1:ncloc,1:nrloc,k+1) = wstokesatw(1:ncloc,1:nrloc,k) - &
                                         & array3dc(:,:,k)
      END WHERE
   ENDDO k_440
   WHERE (maskatc_int) wstokesatw(1:ncloc,1:nrloc,nz+1) = 0.0

ENDIF

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (maskatu,maskatv)
DEALLOCATE (array2dc,array2du,array2dv,array3dc)

CALL log_timer_out(npcc,itm_2dmode)


RETURN

END SUBROUTINE transf_vertical_current

!========================================================================

SUBROUTINE transf_vertical_fall(wsink,utrans,vtrans,wtrans,novars)
!************************************************************************
!
! *transf_vertical_fall* Correct advective velocities for settling effects
!                        in the sigma-grid
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Hydrodynamic_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - transport_at_C_4d1, transport_at_C_4d2 
!
! External calls - 
!
! Module calls - error_alloc, exchange_mod, Warr_at_U, Warr_at_V
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE switches
USE syspars
USE array_interp, ONLY: Warr_at_U, Warr_at_V
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT(IN) :: novars
REAL, INTENT(IN), DIMENSION(0:ncloc,0:nrloc,nz+1,novars) :: wsink
REAL, INTENT(INOUT), DIMENSION(2-nhalo:ncloc+nhalo,nrloc,nz,novars):: utrans
REAL, INTENT(INOUT), DIMENSION(ncloc,2-nhalo:nrloc+nhalo,nz,novars):: vtrans
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz+1,novars) :: wtrans


! Name     Type    Purpose
!------------------------------------------------------------------------------
!*wsink*   REAL    Vertical fall velocity (positive downwards)             [m/s]
!*utrans*  REAL    Physical advective current in the X-direction (on input)
!                  plus sigma-grid correction for settling (on output)     [m/s]
!*vtrans*  REAL    Physical advective current in the Y-direction (on input)
!                  plus sigma-grid correction for settling (on output)     [m/s]
!*wtrans*  REAL    Physical advective current in the vertical direction on
!                  input, sigma-grid corrected on output                   [m/s]
!*novars*  INTEGER Size of variable dimension (number of transformed currents)
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: ivar, k
INTEGER, DIMENSION(4) :: lbounds, nhexch
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: array2d1, array2d2
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: array3d1, array3d2


procname(pglev+1) = 'transform_vertical_fall'
CALL log_timer_in()

!
!1. Allocate
!-----------
!

ALLOCATE (array2d1(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d1',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array2d2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('array2d2',2,(/ncloc,nrloc/),kndrtype)
ALLOCATE (array3d1(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('array3d1',3,(/ncloc,nrloc,nz/),kndrtype)
ALLOCATE (array3d2(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('array3d2',3,(/ncloc,nrloc,nz/),kndrtype)

!
!2. X-direction
!--------------
!

k_210: DO k=1,nz
   WHERE (maskatc_int)
      array3d1(:,:,k) = delzatu(1:ncloc,1:nrloc,k)*&
            & (gscoordatc(1:ncloc,1:nrloc,k)-gscoordatc(0:ncloc-1,1:nrloc,k))&
            & /delxatu(1:ncloc,1:nrloc)
   END WHERE
ENDDO k_210

ivar_220: DO ivar=1,novars
   CALL Warr_at_U(wsink(:,1:nrloc,:,ivar),array3d2,1,(/0,1,1/),&
               & (/ncloc,nrloc,nz+1/),1,0,.TRUE.)
   k_221: DO k=1,nz
      WHERE (nodeatu(1:ncloc,1:nrloc,k).GT.1)
         array2d1 = array3d1(:,:,k)*array3d2(:,:,k)
         utrans(1:ncloc,:,k,ivar) = utrans(1:ncloc,:,k,ivar) - array2d1
      END WHERE
   ENDDO k_221
ENDDO ivar_220

!
!3. Y-direction
!--------------
!

k_310: DO k=1,nz
   WHERE (maskatc_int)
      array3d1(:,:,k) = delzatv(1:ncloc,1:nrloc,k)*&
            & (gscoordatc(1:ncloc,1:nrloc,k)-gscoordatc(1:ncloc,0:nrloc-1,k))&
            & /delxatv(1:ncloc,1:nrloc)
   END WHERE
ENDDO k_310

ivar_320: DO ivar=1,novars
   CALL Warr_at_V(wsink(1:ncloc,:,:,ivar),array3d2,1,(/1,0,1/),&
               & (/ncloc,nrloc,nz+1/),1,0,.TRUE.)
   k_321: DO k=1,nz
      WHERE (nodeatv(1:ncloc,1:nrloc,k).GT.1)
         array2d1 = array3d1(:,:,k)*array3d2(:,:,k)
         vtrans(:,1:nrloc,k,ivar) = vtrans(:,1:nrloc,k,ivar) - array2d1
      END WHERE
   ENDDO k_321
ENDDO ivar_320

!
!4. Vertical direction 
!---------------------
!

wtrans(:,:,1,:) = -wsink(:,:,1,:)
wtrans(:,:,nz+1,:) = 0.0

k_411: DO k=2,nz
   WHERE (maskatc_int)
      array2d1  = delzatw(1:ncloc,1:nrloc,k)*&
          & (gscoordatuw(2:ncloc+1,1:nrloc,k)-gscoordatuw(1:ncloc,1:nrloc,k))&
          & /delxatc(1:ncloc,1:nrloc)
      array2d2  = delzatw(1:ncloc,1:nrloc,k)*&
          & (gscoordatvw(1:ncloc,2:nrloc+1,k)-gscoordatvw(1:ncloc,1:nrloc,k))&
          & /delyatc(1:ncloc,1:nrloc)
      array3d1(:,:,k) = 1.0 - array2d1**2 - array2d2**2
   END WHERE
ENDDO k_411

ivar_412: DO ivar=1,novars
k_412: DO k=2,nz
   WHERE (maskatc_int)
      wtrans(:,:,k,ivar) = wtrans(:,:,k,ivar) - &
                         & wsink(1:ncloc,1:nrloc,k,ivar)*array3d1(:,:,k)
   END WHERE
ENDDO k_412
ENDDO ivar_412

!
!5. Exchange halos
!-----------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = (/2-nhalo,1,1,1/); nhexch = (/nhalo-1,nhalo,0,0/)
   CALL exchange_mod(utrans,lbounds,nhexch,0)
   lbounds = (/1,2-nhalo,1,1/); nhexch = (/0,0,nhalo-1,nhalo/)
   CALL exchange_mod(vtrans,lbounds,nhexch,0)
ENDIF


!
!6. Deallocate
!-------------
!

DEALLOCATE (array2d1,array2d2,array3d1,array3d2)

CALL log_timer_out()


RETURN

END SUBROUTINE transf_vertical_fall
