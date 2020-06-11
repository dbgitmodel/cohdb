! This file is part of COHERENS. You can redistribute and/or modify it under
! the conditions of the COHERENS license. Details are found in the file
! COHERENS_License.

!************************************************************************
!
! *Usrdef_Sediment* User-defined sediment setup (example routine)
!
! Author - Alexander Breugem and Boudewijn Decrop (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.5
!
! $Date: 2013-05-28 15:26:47 +0200 (Tue, 28 May 2013) $
!
! $Revision: 574 $
!
! Description -
!
! Reference -
!
! Routines - usrdef_sed_params, usrdef_sedics, usrdef_sed_spec
!
!************************************************************************
!


!========================================================================

SUBROUTINE usrdef_sed_params
!************************************************************************
!
! *usrdef_sed_params* Define parameters for the sediment model (example routine)
!
! Author - Alexander Breugem and Boudewijn Decrop (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.5
!
! Description 
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE sedpars
USE sedswitches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'usrdef_sed_params'
CALL log_timer_in()

!
!1. Sediment switches
!--------------------
!
!1.1 General
!-----------
!
!---transport mode (1/2/3/4)
iopt_sed_mode = 2
!---grid dimension (2/3)
iopt_sed_nodim = 3
!---type of sediment (1/2)
iopt_sed_type = 2

!
!1.2 Load formula
!----------------
!
!---bed load (1/2/3/4/5/6/7)
iopt_sed_bedeq = 1
!---total load (1/2/3/4/5/6)
iopt_sed_toteq = 1

!
!1.3 Bottom stress
!-----------------
!
!---type of roughness (1/2/3)
iopt_sed_tau = 1
!---formulation for critical shear stress (0/1/2/3)
iopt_sed_taucr = 1
!---hiding factor (0/1/2)
iopt_sed_hiding = 0

!
!1.4 Settling
!------------
!
!---type of settling velocity (0/1/2/3/4/5)
iopt_sed_ws = 1
!---hindered settling (0/1/2)
iopt_sed_hindset = 0
!---flocculation (0/1/2/3)
iopt_sed_floc = 0
 
!
!1.5 Bottom boundary condition
!-----------------------------
!
!---type of near bed b.c. (0/1/2/3)
iopt_sed_bbc = 1
!---method for applying bed b.c. (1/2/3)
iopt_sed_bbc_type = 1
!---formulation for reference concentration (1/2/3/4/5)
iopt_sed_ceqeq = 1
!---form of Engelund-Hansen equation (1/2)
iopt_sed_eha = 1
!---slope effects (0/1)
iopt_sed_slope = 0

!
!1.6 Diffusion
!-------------
!
!---beta-factor (1/2/3)
iopt_sed_beta = 1
!---effect of waves (0/1)
iopt_sed_wave_diff = 0

!
!1.7 Density effects
!-------------------
!
!---disable/enable density effects (0/1)
iopt_sed_dens = 0

!
!1.8 Particle attributes
!-----------------------
!
!---method for median size (1/2)
iopt_sed_median = 1

!
!1.9 Numerical
!-------------
!
!---settling (0/1/2/3)
iopt_sed_vadv = 3
!---filtering (0/1)
iopt_sed_filter = 0

!
!2. Sediment parameters
!----------------------
!
!---number of fractions
nf = 1

!---Gauss-Legendre integration
nrquad_sed = 7; nrquad_wav = 10

!---number of iteration in Bartnicki filter
maxitbartnicki = 100

!---settling
n_RichZaki = 4.6

!---beta factor for diffusion
beta_sed_cst = 1.0
beta_sed_max = 1.5
beta_sed_min = 1.0

!---bottom boundary conditions
maxRV = 0.1
minRV = 1.0E-05,
z0_coef = 30.0
height_c_cst = 0.01

!---shear stress
parth_coef = 1.0E-08   ! [m^3/m^2/s]
parth_exp = 1.0
wu_exp = -0.6

!---flocculation
a_leussen = 0.02       ! [s]
b_leussen = 0.0024     ! [S^2]
alpha_VR = 2.19
floc_VR_max = 10.0
floc_VR_min = 1.0

!---slope factor
coef_bed_grad = 1.3

!---gelling concentration [m^3/m^3]
cgel = 0.0

!---maximum concentration in load formulae [m^3/m^3]
cmax = 0.65

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sed_params

!========================================================================

SUBROUTINE usrdef_sedics
!************************************************************************
!
! *usrdef_sedics* Define initial conditions for sediments (example routine)
!
! Author - Alexander Breugem and Boudewijn Decrop (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.5
!
! Description - initial arrays are defined on local grid including halos
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE gridpars
USE iopars
USE sedarrays
USE sedids
USE sedpars
USE syspars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_filepars, open_filepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iunit, method = ?
INTEGER, DIMENSION(4) :: lbounds, nhdist = 0
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realglb2d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realglb3d
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: realglb4d


procname(pglev+1) = 'usrdef_sedics'
CALL log_timer_in()

!
!1. "Local" method
!-----------------
!

IF (method.EQ.1) THEN

!  ---volumetric concentrations [m^3/m^3]
   cvol(1:ncloc,1:nrloc,n1:nz,1:nf) = ?

!  ---skin roughness [m]
   zroughatc_sed(1:ncloc,1:nrloc) = ?

!  ---bed fractions
   bed_fraction(1:ncloc,1:nrloc,1:nf) = ?

!
!2. "Global" method
!------------------
!

ELSEIF (method.EQ.2) THEN

!
!2.1 Allocate work space arrays
!------------------------------
!

   ALLOCATE (realglb2d(nc,nr),STAT=errstat)
   CALL error_alloc('realglb2d',2,(/nc,nr/),real_type)
   realglb2d = 0.0
   ALLOCATE (realglb3d(nc,nr,nf),STAT=errstat)
   CALL error_alloc('realglb3d',2,(/nc,nr,nf/),real_type)
   realglb3d = 0.0
   ALLOCATE (realglb4d(nc,nr,nz,nf),STAT=errstat)
   CALL error_alloc('realglb4d',2,(/nc,nr,nz,nf/),real_type)
   realglb4d = 0.0

!
!2.2 Open data file
!------------------
!
!  ---it is assumed that the name of the data file is defined in 
!     usrdef_mod_params
   CALL open_filepars(modfiles(io_inicon,ics_sed,1))
   iunit = modfiles(io_inicon,ics_sed,1)%iunit

!
!2.3 Read data
!-------------
!
!  ---volumetric concentrations
   lbounds = (/1-nhalo,1-nhalo,1,1/)
   READ (iunit,?) realglb4d
   CALL distribute_mod(realglb3d,cvol,lbounds,nhdist,iarr_cvol,0.0)

!  ---skin roughness
   lbounds(1:2) = (/0,0/)
   READ (iunit,?) realglg2d
   CALL distribute_mod(realglb2d,zroughatc_sed,lbounds(1:2),nhdist,&
                     & iarr_zroughatc_sed,0.0)

!  ---bed fractions
   lbounds(1:3) = (/0,0,1/)
   READ (iunit,?) realglg3d
   CALL distribute_mod(realglb3d,bed_fraction,lbounds(1:3),nhdist,&
                     & iarr_bed_fraction,0.0)

!
!2.4 Close data file
!--------------------
!

   CALL close_filepars(modfiles(io_inicon,ics_phys,1))

!
!2.5 Deallocate
!--------------
!

   DEALLOCATE (realglb2d,realgbl3d)

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sedics

!========================================================================

SUBROUTINE usrdef_sed_spec
!************************************************************************
!
! *usrdef_sed_spec* Define particle properties of sediment fractions
!                   (example routine)
!
! Author - Alexander Breugem and Boudewijn Decrop (IMDC)
!
! Version - @(COHERENS)Usrdef_Sediment.f90  V2.5
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE sedarrays
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'usrdef_sed_spec'
CALL log_timer_in()

!---diameter [m]
dp = ?

!---density [m^3/s]
rhos = ?

!---(uniform) critical stress [m^2/s^2]
tau_cr_cst = ?

!---(uniform) setling velocity [m/s]
ws_cst =  ?

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_sed_spec


