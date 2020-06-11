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
! *Allocate_Arrays* Allocate/deallocate memory for arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Allocate_Arrays.f90  V2.11.2
!
! $Date: 2018-08-10 16:33:34 +0200 (Fri, 10 Aug 2018) $
!
! $Revision: 1177 $
!
! Description -
!
! Routines - allocate_global_grid, allocate_mod_arrays, deallocate_mod_arrays
!
!************************************************************************
!

!========================================================================

SUBROUTINE allocate_global_grid
!************************************************************************
!
! *allocate_global_grid* Allocate global grid arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Allocate_Arrays.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
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
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'allocate_global_grid'
CALL log_timer_in()

!
!1. Grid coordinates
!-------------------
!

ALLOCATE (gxcoordglb(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('gxcoordglb',2,(/nc+2,nr+2/),kndrtype)
gxcoordglb = 0.0

ALLOCATE (gycoordglb(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('gycoordglb',2,(/nc+2,nr+2/),kndrtype)
gycoordglb = 0.0

ALLOCATE (gdelxglb(0:nc+1),STAT=errstat)
CALL error_alloc('gdelxglb',1,(/nc+2/),kndrtype)
gdelxglb = 0.0

ALLOCATE (gdelyglb(0:nr+1),STAT=errstat)
CALL error_alloc('gdelyglb',1,(/nr+2/),kndrtype)
gdelyglb = 0.0

ALLOCATE (gscoordglb(0:nc+1,0:nr+1,nz+1),STAT=errstat)
CALL error_alloc('gscoordglb',3,(/nc+2,nr+2,nz+1/),kndrtype)
gscoordglb = 0.0

ALLOCATE (gsigcoordatc(nz),STAT=errstat)
CALL error_alloc('gsigcoordatc',1,(/nz/),kndrtype)
gsigcoordatc = 0.0

ALLOCATE (gsigcoordatw(nz+1),STAT=errstat)
CALL error_alloc('gsigcoordatw',1,(/nz+1/),kndrtype)
gsigcoordatw = 0.0

!
!2. Mean water depths
!--------------------
!

ALLOCATE (depmeanglb(1-nhalo:nc+nhalo,1-nhalo:nr+nhalo),STAT=errstat)
CALL error_alloc('depmeanglb',2,(/nc+2*nhalo,nr+2*nhalo/),kndrtype)
depmeanglb = depmean_flag

ALLOCATE (seapointglb(0:nc+1,0:nr+1),STAT=errstat)
CALL error_alloc('seapointglb',2,(/nc+2,nr+2/),kndrtype)
seapointglb = .FALSE.

!
!3. Indices of open boundary locations
!-------------------------------------
!

ALLOCATE (iobu(nobu),STAT=errstat)
CALL error_alloc('iobu',1,(/nobu/),kndint)
ALLOCATE (jobu(nobu),STAT=errstat)
CALL error_alloc('jobu',1,(/nobu/),kndint)

ALLOCATE (iobv(nobv),STAT=errstat)
CALL error_alloc('iobv',1,(/nobv/),kndint)
ALLOCATE (jobv(nobv),STAT=errstat)
CALL error_alloc('jobv',1,(/nobv/),kndint)

CALL log_timer_out()


RETURN

END SUBROUTINE allocate_global_grid

!========================================================================

SUBROUTINE allocate_mod_arrays
!************************************************************************
!
! *allocate_mod_arrays* Allocate model arrays
!
! Author - Patrick Luyten
!  
! Version - @(COHERENS)Allocate_Arrays.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - error_alloc, error_alloc_struc, num_halo
!
!************************************************************************
!
USE currents
USE density
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE obconds
USE optics
USE paralpars
USE structures
USE switches
USE syspars
USE tide
USE timepars
USE turbulence
USE wavevars
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: num_halo

IMPLICIT NONE

!
!*Local variables
!

INTEGER :: nhadvd, nhadvw, nhadv2, nhadv3, nhdifd, nhdifw, nhdif2, nhdif3

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*nhadvd*     INTEGER  Halo size for horizontal advection of density    (0/1/2)
!*nhadvw*     INTEGER  Halo size for horizontal advection of turbulence (0/1/2)
!*nhadv2*     INTEGER  Halo size for horizontal advection of 2-D transport
!                                                                       (0/1/2)
!*nhadv3*     INTEGER  Halo size for horizontal advection of 3-D currents
!                                                                       (0/1/2)
!*nhdifd*     INTEGER  Halo size for horizontal diffusion of density      (0/1)
!*nhdifw*     INTEGER  Halo size for horizontal diffusion of turbulence   (0/1)
!*nhdif2*     INTEGER  Halo size for horizontal diffusion of 2-D transport
!                                                                         (0/1)
!*nhdif3*     INTEGER  Halo size for horizontal diffusion of 3-D currents (0/1)
!
!************************************************************************
!


procname(pglev+1) = 'allocate_mod_arrays'
CALL log_timer_in()

!
!1. Halo sizes
!-------------
!
!---halo sizes for advection/diffusion
nhadvd = num_halo(iopt_adv_scal)
nhadvw = num_halo(iopt_adv_turb)
nhadv2 = num_halo(iopt_adv_2D)
nhadv3 = num_halo(iopt_adv_3D)
nhdifd = iopt_hdif_scal
nhdifw = iopt_hdif_turb
nhdif2 = iopt_hdif_2D
nhdif3 = iopt_hdif_3D

!---halos sizes for global arrays
nhdens = MERGE(2,MAX(nhadvd,1),iopt_dens_grad.EQ.3)
nhscal = MAX(nhadvd,1)
nhfvel = MAX(nhadvd,nhadvw,1)
nhturb = MAX(nhadvw,nhdifw)
nh2vel = MAX(nhadv2,1)
nh3vel = MAX(nhadv3,1)

!
!2. Grid arrays
!--------------
!

ALLOCATE (coriolatu(ncloc+1,0:nrloc),STAT=errstat)
CALL error_alloc('coriolatu',2,(/ncloc+1,nrloc+1/),kndrtype)
coriolatu = 0.0

ALLOCATE (coriolatv(0:ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('coriolatv',2,(/ncloc+1,nrloc+1/),kndrtype)
coriolatv = 0.0

ALLOCATE (delxatc(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('delxatc',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
delxatc = 0.0

ALLOCATE (delxatu(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('delxatu',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
delxatu = 0.0

ALLOCATE (delxatuv(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('delxatuv',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
delxatuv = 0.0

ALLOCATE (delxatv(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('delxatv',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
delxatv = 0.0

ALLOCATE (delyatc(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('delyatc',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
delyatc = 0.0

ALLOCATE (delyatu(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('delyatu',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
delyatu = 0.0

ALLOCATE (delyatuv(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('delyatuv',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
delyatuv = 0.0

ALLOCATE (delyatv(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('delyatv',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
delyatv = 0.0

ALLOCATE (delzatc(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('delzatc',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
delzatc = 0.0

ALLOCATE (delzatu(0:ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('delzatu',3,(/ncloc+2,nrloc,nz/),kndrtype)
delzatu = 0.0

ALLOCATE (delzatuv(ncloc+1,nrloc+1,nz),STAT=errstat)
CALL error_alloc('delzatuv',3,(/ncloc+1,nrloc+1,nz/),kndrtype)
delzatuv = 0.0

ALLOCATE (delzatuw(ncloc+1,nrloc,nz+1),STAT=errstat)
CALL error_alloc('delzatuw',3,(/ncloc+1,nrloc,nz+1/),kndrtype)
IF (nz.GT.0) delzatuw = 0.0

ALLOCATE (delzatv(ncloc,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('delzatv',3,(/ncloc,nrloc+2,nz/),kndrtype)
delzatv = 0.0

ALLOCATE (delzatvw(ncloc,nrloc+1,nz+1),STAT=errstat)
CALL error_alloc('delzatvw',3,(/ncloc,nrloc+1,nz+1/),kndrtype)
IF (nz.GT.0) delzatvw = 0.0

ALLOCATE (delzatw(0:ncloc+1,0:nrloc+1,nz+1),STAT=errstat)
CALL error_alloc('delzatw',3,(/ncloc+2,nrloc+2,nz+1/),kndrtype)
IF (nz.GT.0) delzatw = 0.0

ALLOCATE (gaccatc(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('gaccatc',2,(/ncloc+2,nrloc+2/),kndrtype)
gaccatc = 0.0

ALLOCATE (gaccatu(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('gaccatu',2,(/ncloc+1,nrloc/),kndrtype)
gaccatu = 0.0

ALLOCATE (gaccatv(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('gaccatv',2,(/ncloc,nrloc+1/),kndrtype)
gaccatv = 0.0

ALLOCATE (gangleatc(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('gangleatc',2,(/ncloc+2,nrloc+2/),kndrtype)
gangleatc = 0.0

ALLOCATE (garea(ncloc,nrloc),STAT=errstat)
CALL error_alloc('garea',2,(/ncloc,nrloc/),kndrtype)
garea = 0.0

ALLOCATE (gscoordatc(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('gscoordatc',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
gscoordatc = 0.0

ALLOCATE (gscoordatu(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('gscoordatu',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
gscoordatu = 0.0

ALLOCATE (gscoordatuvw(0:ncloc+1,0:nrloc+1,nz+1),STAT=errstat)
CALL error_alloc('gscoordatuvw',3,(/ncloc+2,nrloc+2,nz+1/),kndrtype)
gscoordatuvw = 0.0

ALLOCATE (gscoordatuw(0:ncloc+1,0:nrloc+1,nz+1),STAT=errstat)
CALL error_alloc('gscoordatuw',3,(/ncloc+2,nrloc+2,nz+1/),kndrtype)
gscoordatuw = 0.0

ALLOCATE (gscoordatv(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('gscoordatv',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
gscoordatv = 0.0

ALLOCATE (gscoordatvw(0:ncloc+1,0:nrloc+1,nz+1),STAT=errstat)
CALL error_alloc('gscoordatvw',3,(/ncloc+2,nrloc+2,nz+1/),kndrtype)
gscoordatvw = 0.0

ALLOCATE (gscoordatw(0:ncloc+1,0:nrloc+1,nz+1),STAT=errstat)
CALL error_alloc('gscoordatw',3,(/ncloc+2,nrloc+2,nz+1/),kndrtype)
gscoordatw = 0.0

ALLOCATE (maskatc_int(ncloc,nrloc),STAT=errstat)
CALL error_alloc('maskatc_int',2,(/ncloc,nrloc/),kndlog)
maskatc_int = .FALSE.

ALLOCATE (nodeatc(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('nodeatc',2,(/ncloc+nhalo+1,nrloc+nhalo+1/),kndint)
nodeatc = 0

ALLOCATE (nodeatu(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('nodeatu',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndint)
nodeatu = 0

ALLOCATE (nodeatv(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('nodeatv',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndint)
nodeatv = 0

ALLOCATE (nodeatuv(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('nodeatuv',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndint)
nodeatuv = 0

ALLOCATE (nodeatuw(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1),STAT=errstat)
CALL error_alloc('nodeatuw',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz+1/),kndint)
nodeatuw = 0

ALLOCATE (nodeatvw(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1),STAT=errstat)
CALL error_alloc('nodeatvw',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz+1/),kndint)
nodeatvw = 0

ALLOCATE (node2du(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('node2du',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndint)
node2du = 0

ALLOCATE (node2duv(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('node2duv',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndint)
node2duv = 0

ALLOCATE (node2dv(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('node2dv',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndint)
node2dv = 0

ALLOCATE (gxlon(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('gxlon',2,(/ncloc+1,nrloc+1/),kndrtype)
gxlon = 0.0

ALLOCATE (gylat(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('gylat',2,(/ncloc+1,nrloc+1/),kndrtype)
gylat = 0.0

ALLOCATE (gxcoord(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('gxcoord',2,(/ncloc+2,nrloc+2/),kndrtype)
gxcoord = 0.0

ALLOCATE (gxcoordglbatc(nc,nr),STAT=errstat)
CALL error_alloc('gxcoordglbatc',2,(/nc,nr/),kndrtype)
gxcoordglbatc = 0.0

ALLOCATE (gycoord(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('gycoord',2,(/ncloc+2,nrloc+2/),kndrtype)
gycoord = 0.0

ALLOCATE (gycoordglbatc(nc,nr),STAT=errstat)
CALL error_alloc('gycoordglbatc',2,(/nc,nr/),kndrtype)
gycoordglbatc = 0.0

ALLOCATE (rmaskatc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('rmaskatc',2,(/ncloc,nrloc/),kndrtype)
rmaskatc = 0.0

ALLOCATE (seapoint(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('seapoint',2,(/ncloc+2,nrloc+2/),kndlog)
seapoint = .FALSE.

IF (iopt_fld.GT.0) THEN
   ALLOCATE (alphatc_fld(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('alphatc_fld',2,(/ncloc,nrloc/),kndrtype)
   alphatc_fld = 1.0
   ALLOCATE (alphatu_fld(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('alphatu_fld',2,(/ncloc,nrloc/),kndrtype)
   alphatu_fld = 1.0
   ALLOCATE (alphatv_fld(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('alphatv_fld',2,(/ncloc,nrloc/),kndrtype)
   alphatv_fld = 1.0
ENDIF

IF (iopt_obc_advrlx.EQ.1) THEN
   ALLOCATE (rlxobcatu(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('rlxobcatu',2,(/ncloc,nrloc/),kndrtype)
   rlxobcatu = 0.0
   ALLOCATE (rlxobcatv(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('rlxobcatv',2,(/ncloc,nrloc/),kndrtype)
   rlxobcatv = 0.0
ENDIF

!
!3. Water depths and elevations
!------------------------------
!

ALLOCATE (depmeanatc(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('depmeanatc',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
depmeanatc = 0.0

ALLOCATE (depmeanatu(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('depmeanatu',2,(/ncloc+2,nrloc+2/),kndrtype)
depmeanatu = 0.0

ALLOCATE (depmeanatuv(1-nhalo:ncloc+1+nhalo,1-nhalo:nrloc+1+nhalo),STAT=errstat)
CALL error_alloc('depmeanatuv',2,(/ncloc+1+2*nhalo,nrloc+1+2*nhalo/),kndrtype)
depmeanatuv = 0.0

ALLOCATE (depmeanatv(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('depmeanatv',2,(/ncloc+2,nrloc+2/),kndrtype)
depmeanatv = 0.0

ALLOCATE (deptotatc(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('deptotatc',2,(/ncloc+2,nrloc+2/),kndrtype)
deptotatc = 0.0

ALLOCATE (deptotatc_old(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('deptotatc_old',2,(/ncloc+2,nrloc+2/),kndrtype)
deptotatc_old = 0.0

ALLOCATE (deptotatu(1-nhalo:ncloc+nhalo,0:nrloc+1),STAT=errstat)
CALL error_alloc('deptotatu',2,(/ncloc+2*nhalo,nrloc+2/),kndrtype)
deptotatu = 0.0

ALLOCATE (deptotatu_old(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('deptotatu_old',2,(/ncloc+2,nrloc/),kndrtype)
deptotatu_old = 0.0

ALLOCATE (deptotatuv(ncloc+1,nrloc+1),STAT=errstat)
CALL error_alloc('deptotatuv',2,(/ncloc+1,nrloc+1/),kndrtype)
deptotatuv = 0.0

ALLOCATE (deptotatv(0:ncloc+1,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('deptotatv',2,(/ncloc+2,nrloc+2*nhalo/),kndrtype)
deptotatv = 0.0

ALLOCATE (deptotatv_old(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('deptotatv_old',2,(/ncloc,nrloc+2/),kndrtype)
deptotatv_old = 0.0

ALLOCATE (zeta(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('zeta',2,(/ncloc+2,nrloc+2/),kndrtype)
zeta = 0.0

IF (iopt_hydro_impl.EQ.1) THEN
   ALLOCATE (dzeta(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('dzeta',2,(/ncloc+2,nrloc+2/),kndrtype)
   dzeta = 0.0
   ALLOCATE (zeta_old(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('zeta_old',2,(/ncloc+2,nrloc+2/),kndrtype)
   zeta_old = 0.0
ENDIF

IF (iopt_fld.GT.0) THEN
   ALLOCATE (deptotatc_err(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('deptotatc_err',2,(/ncloc,nrloc/),kndrtype)
   deptotatc_err = 0.0
   ALLOCATE (deptotatc_prev(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('deptotatc_prev',2,(/ncloc,nrloc/),kndrtype)
   deptotatc_prev = 0.0
ENDIF

!
!4. Velocity arrays
!------------------
!

ALLOCATE (p2dbcgradatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('p2dbcgradatu',2,(/ncloc,nrloc/),kndrtype)
p2dbcgradatu = 0.0

ALLOCATE (p2dbcgradatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('p2dbcgradatv',2,(/ncloc,nrloc/),kndrtype)
p2dbcgradatv = 0.0

ALLOCATE (p3dbcgradatu(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('p3dbcgradatu',3,(/ncloc,nrloc,nz/),kndrtype)
p3dbcgradatu = 0.0

ALLOCATE (p3dbcgradatv(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('p3dbcgradatv',3,(/ncloc,nrloc,nz/),kndrtype)
p3dbcgradatv = 0.0

ALLOCATE (udevint(ncloc,nrloc),STAT=errstat)
CALL error_alloc('udevint',2,(/ncloc,nrloc/),kndrtype)
udevint = 0.0

ALLOCATE (udfvel(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('udfvel',2,(/ncloc+1,nrloc/),kndrtype)
udfvel = 0.0

ALLOCATE (udvel(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('udvel',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
udvel = 0.0

ALLOCATE (udvel_old(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('udvel_old',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
udvel_old = 0.0

ALLOCATE (ufvel(2-nhalo:ncloc+nhalo,nrloc,nz),STAT=errstat)
CALL error_alloc('ufvel',3,(/ncloc+2*nhalo-1,nrloc,nz/),kndrtype)
ufvel = 0.0

ALLOCATE (umpred(0:ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('umpred',2,(/ncloc+2,nrloc/),kndrtype)
umpred = 0.0

ALLOCATE (umvel(1-nhalo:ncloc+nhalo,0:nrloc+1),STAT=errstat)
CALL error_alloc('umvel',2,(/ncloc+2*nhalo,nrloc+2/),kndrtype)
umvel = 0.0

ALLOCATE (umvel_old(1-nhalo:ncloc+nhalo,0:nrloc+1),STAT=errstat)
CALL error_alloc('umvel_old',2,(/ncloc+2*nhalo,nrloc+2/),kndrtype)
umvel_old = 0.0

ALLOCATE (uvel(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('uvel',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
uvel = 0.0

ALLOCATE (uvel_old(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('uvel_old',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
uvel_old = 0.0

ALLOCATE (vdevint(ncloc,nrloc),STAT=errstat)
CALL error_alloc('vdevint',2,(/ncloc,nrloc/),kndrtype)
vdevint = 0.0

ALLOCATE (vdfvel(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('vdfvel',2,(/ncloc,nrloc+1/),kndrtype)
vdfvel = 0.0

ALLOCATE (vdvel(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('vdvel',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
vdvel = 0.0

ALLOCATE (vdvel_old(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('vdvel_old',2,(/ncloc+2*nhalo,nrloc+2*nhalo/),kndrtype)
vdvel_old = 0.0

ALLOCATE (vfvel(ncloc,2-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('vfvel',3,(/ncloc,nrloc+2*nhalo-1,nz/),kndrtype)
vfvel = 0.0

ALLOCATE (vmpred(ncloc,0:nrloc+1),STAT=errstat)
CALL error_alloc('vmpred',2,(/ncloc,nrloc+2/),kndrtype)
vmpred = 0.0

ALLOCATE (vmvel(0:ncloc+1,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('vmvel',2,(/ncloc+2,nrloc+2*nhalo/),kndrtype)
vmvel = 0.0

ALLOCATE (vmvel_old(0:ncloc+1,1-nhalo:nrloc+nhalo),STAT=errstat)
CALL error_alloc('vmvel_old',2,(/ncloc+2,nrloc+2*nhalo/),kndrtype)
vmvel_old = 0.0

ALLOCATE (wphys(ncloc,nrloc,nz),STAT=errstat)
CALL error_alloc('wphys',3,(/ncloc,nrloc,nz/),kndrtype)
wphys = 0.0

ALLOCATE (wfvel(0:ncloc,0:nrloc,nz+1),STAT=errstat)
CALL error_alloc('wfvel',3,(/ncloc+1,nrloc+1,nz+1/),kndrtype)
wfvel = 0.0

ALLOCATE (wvel(0:ncloc,0:nrloc,nz+1),STAT=errstat)
CALL error_alloc('wvel',3,(/ncloc+1,nrloc+1,nz+1/),kndrtype)
wvel = 0.0

ALLOCATE (vvel(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('vvel',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
vvel = 0.0

ALLOCATE (vvel_old(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('vvel_old',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
vvel_old = 0.0

!
!5. Density arrays
!-----------------
!

ALLOCATE (beta_sal(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('beta_sal',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
beta_sal = 0.0

ALLOCATE (beta_temp(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('beta_temp',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
beta_temp = 0.0

ALLOCATE (dens(0:ncloc+1,0:nrloc+1,nz),STAT=errstat)
CALL error_alloc('dens',3,(/ncloc+2,nrloc+2,nz/),kndrtype)
dens = 0.0

ALLOCATE (sal(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('sal',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
sal = 0.0

ALLOCATE (temp(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
CALL error_alloc('temp',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz/),kndrtype)
temp = 0.0

!
!6. Diffusion coefficients
!-------------------------
!

ALLOCATE (hdifcoef2datc(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('hdifcoef2datc',2,(/ncloc+1,nrloc+1/),kndrtype)
hdifcoef2datc = 0.0

ALLOCATE (hdifcoef2datuv(ncloc+1,nrloc+1),STAT=errstat)
CALL error_alloc('hdifcoef2datuv',2,(/ncloc+1,nrloc+1/),kndrtype)
hdifcoef2datuv = 0.0

ALLOCATE (hdifcoef3datc(0:ncloc,0:nrloc,nz),STAT=errstat)
CALL error_alloc('hdifcoef3datc',3,(/ncloc+1,nrloc+1,nz/),kndrtype)
hdifcoef3datc = 0.0

ALLOCATE (hdifcoef3datu(ncloc+1,nrloc,nz),STAT=errstat)
CALL error_alloc('hdifcoef3datu',3,(/ncloc+1,nrloc,nz/),kndrtype)
hdifcoef3datu = 0.0

ALLOCATE (hdifcoef3datuv(ncloc+1,nrloc+1,nz),STAT=errstat)
CALL error_alloc('hdifcoef3datuv',3,(/ncloc+1,nrloc+1,nz/),kndrtype)
hdifcoef3datuv = 0.0

ALLOCATE (hdifcoef3datv(ncloc,nrloc+1,nz),STAT=errstat)
CALL error_alloc('hdifcoef3datv',3,(/ncloc,nrloc+1,nz/),kndrtype)
hdifcoef3datv = 0.0

ALLOCATE (kinvisc(0:ncloc,0:nrloc,nz+1),STAT=errstat)
CALL error_alloc('kinvisc',3,(/ncloc+1,nrloc+1,nz+1/),kndrtype)
kinvisc = 0.0

ALLOCATE (vdifcoefmom(0:ncloc,0:nrloc,nz+1),STAT=errstat)
CALL error_alloc('vdifcoefmom',3,(/ncloc+1,nrloc+1,nz+1/),kndrtype)
vdifcoefmom = 0.0

ALLOCATE (vdifcoefscal(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('vdifcoefscal',3,(/ncloc,nrloc,nz+1/),kndrtype)
vdifcoefscal = 0.0

ALLOCATE (vdifcoeftke(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('vdifcoeftke',3,(/ncloc,nrloc,nz+1/),kndrtype)
vdifcoeftke = 0.0

IF (iopt_hdif_scal.GE.2) THEN

   ALLOCATE (hdifcoef3datw(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('hdifcoef3datc',3,(/ncloc,nrloc,nz+1/),kndrtype)
   hdifcoef3datw = 0.0

   ALLOCATE (xslopeatu_geo(ncloc+1,nrloc,nz),STAT=errstat)
   CALL error_alloc('xslopeatu_geo',3,(/ncloc+1,nrloc,nz/),kndrtype)
   xslopeatu_geo = 0.0

   ALLOCATE (yslopeatv_geo(ncloc,nrloc+1,nz),STAT=errstat)
   CALL error_alloc('yslopeatv_geo',3,(/ncloc,nrloc+1,nz/),kndrtype)
   yslopeatv_geo = 0.0

   ALLOCATE (vdifcoefscal_rot(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('vdifcoefscal_rot',3,(/ncloc,nrloc,nz+1/),kndrtype)
   vdifcoefscal_rot = 0.0

ENDIF

IF (iopt_hdif_scal.EQ.2) THEN

   ALLOCATE (xslopeatw_geo(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('xslopeatw_geo',3,(/ncloc,nrloc,nz+1/),kndrtype)
   xslopeatw_geo = 0.0

   ALLOCATE (yslopeatw_geo(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('yslopeatw_geo',3,(/ncloc,nrloc,nz+1/),kndrtype)
   yslopeatw_geo = 0.0

ENDIF

IF (iopt_hdif_scal.EQ.3) THEN

   ALLOCATE (xslopeatu_siso(0:ncloc+1,nrloc,nz,0:1,0:1),STAT=errstat)
   CALL error_alloc('xslopeatu_siso',5,(/ncloc+2,nrloc,nz,2,2/),kndrtype)
   xslopeatu_siso = 0.0

   ALLOCATE (xslopeatu_ziso(0:ncloc+1,nrloc,nz,0:1,0:1),STAT=errstat)
   CALL error_alloc('xslopeatu_ziso',5,(/ncloc+2,nrloc,nz,2,2/),kndrtype)
   xslopeatu_ziso = 0.0

   ALLOCATE (yslopeatv_siso(ncloc,0:nrloc+1,nz,0:1,0:1),STAT=errstat)
   CALL error_alloc('yslopeatv_siso',5,(/ncloc,nrloc+2,nz,2,2/),kndrtype)
   yslopeatv_siso = 0.0

   ALLOCATE (yslopeatv_ziso(ncloc,0:nrloc+1,nz,0:1,0:1),STAT=errstat)
   CALL error_alloc('yslopeatv_ziso',5,(/ncloc,nrloc+2,nz,2,2/),kndrtype)
   yslopeatv_ziso = 0.0

ENDIF

!
!7. Turbulence arrays
!--------------------
!

ALLOCATE (buofreq2(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('buofreq2',3,(/ncloc,nrloc,nz+1/),kndrtype)
buofreq2 = 0.0

ALLOCATE (dissip(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1),STAT=errstat)
CALL error_alloc('dissip',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz+1/),kndrtype)
dissip = 0.0

ALLOCATE (shearfreq2(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('shearfreq2',3,(/ncloc,nrloc,nz+1/),kndrtype)
shearfreq2 = 0.0

ALLOCATE (tke(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1),STAT=errstat)
CALL error_alloc('tke',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz+1/),kndrtype)
tke = 0.0

ALLOCATE (tke_old(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1),STAT=errstat)
CALL error_alloc('tke_old',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz+1/),kndrtype)
tke_old = 0.0

ALLOCATE (zlmix(1-nhalo:ncloc+nhalo,1-nhalo:nrloc+nhalo,nz+1),STAT=errstat)
CALL error_alloc('zlmix',3,(/ncloc+2*nhalo,nrloc+2*nhalo,nz+1/),kndrtype)
zlmix = 0.0

!
!8. Tidal arrays
!---------------
!

ALLOCATE (fnode_astro(nconastro),STAT=errstat)
CALL error_alloc('fnode_astro',1,(/nconastro/),kndrtype)
IF (nconastro.GT.0) fnode_astro = 1.0

ALLOCATE (fnode_obc(nconobc),STAT=errstat)
CALL error_alloc('fnode_obc',1,(/nconobc/),kndrtype)
IF (nconobc.GT.0) fnode_obc = 1.0

ALLOCATE (fxastro(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('fxastro',2,(/ncloc+1,nrloc/),kndrtype)
fxastro = 0.0

ALLOCATE (fyastro(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('fyastro',2,(/ncloc,nrloc+1/),kndrtype)
fyastro = 0.0

ALLOCATE (phase_astro(nconastro),STAT=errstat)
CALL error_alloc('phase_astro',1,(/nconastro/),kndrtype)
IF (nconastro.GT.0) phase_astro = 0.0

ALLOCATE (phase_obc(nconobc),STAT=errstat)
CALL error_alloc('phase_obc',1,(/nconobc/),kndrtype)
IF (nconobc.GT.0) phase_obc = 0.0

!
!9. Meteorological arrays
!------------------------
!

ALLOCATE (airtemp(ncloc,nrloc),STAT=errstat)
CALL error_alloc('airtemp',2,(/ncloc,nrloc/),kndrtype)
airtemp = 0.0

ALLOCATE (atmpres(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('atmpres',2,(/ncloc+1,nrloc+1/),kndrtype)
atmpres = 0.0

ALLOCATE (cloud_cover(ncloc,nrloc),STAT=errstat)
CALL error_alloc('cloud_cover',2,(/ncloc,nrloc/),kndrtype)
cloud_cover = 0.0

ALLOCATE (evapminprec(ncloc,nrloc),STAT=errstat)
CALL error_alloc('evapminprec',2,(/ncloc,nrloc/),kndrtype)
evapminprec = 0.0

ALLOCATE (evaporation(ncloc,nrloc),STAT=errstat)
CALL error_alloc('evaporation',2,(/ncloc,nrloc/),kndrtype)
evaporation = 0.0

ALLOCATE (precipitation(ncloc,nrloc),STAT=errstat)
CALL error_alloc('precipitation',2,(/ncloc,nrloc/),kndrtype)
precipitation = 0.0

ALLOCATE (qspecdif(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qspecdif',2,(/ncloc,nrloc/),kndrtype)
qspecdif = 0.0

ALLOCATE (relhum(ncloc,nrloc),STAT=errstat)
CALL error_alloc('relhum',2,(/ncloc,nrloc/),kndrtype)
relhum = 0.0

ALLOCATE (sst(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sst',2,(/ncloc,nrloc/),kndrtype)
sst = 0.0

ALLOCATE (tempdif(ncloc,nrloc),STAT=errstat)
CALL error_alloc('tempdif',2,(/ncloc,nrloc/),kndrtype)
tempdif = 0.0

ALLOCATE (uwindatc(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('uwindatc',2,(/ncloc+1,nrloc+1/),kndrtype)
uwindatc = 0.0

ALLOCATE (vappres_air(ncloc,nrloc),STAT=errstat)
CALL error_alloc('vappres_air',2,(/ncloc,nrloc/),kndrtype)
vappres_air = 0.0

ALLOCATE (vwindatc(0:ncloc+1,0:nrloc+1),STAT=errstat)
CALL error_alloc('vwindatc',2,(/ncloc+2,nrloc+2/),kndrtype)
vwindatc = 0.0

ALLOCATE (windatc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('windatc',2,(/ncloc,nrloc/),kndrtype)
windatc = 0.0

IF (iopt_part_model.GT.0) THEN

   ALLOCATE (uwindatc_old(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('uwindatc_old',2,(/ncloc+1,nrloc+1/),kndrtype)
   uwindatc_old = 0.0

   ALLOCATE (vwindatc_old(0:ncloc+1,0:nrloc+1),STAT=errstat)
   CALL error_alloc('vwindatc_old',2,(/ncloc+1,nrloc+1/),kndrtype)
   vwindatc_old = 0.0

ENDIF

!
!10. Bottom and surface fluxes
!-----------------------------
!
!10.1 Bottom fluxes
!------------------
!

ALLOCATE (bdragcoefatc(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('bdragcoefatc',2,(/ncloc+1,nrloc+1/),kndrtype)
bdragcoefatc = 0.0

ALLOCATE (bdragcoefatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bdragcoefatu',2,(/ncloc,nrloc/),kndrtype)
bdragcoefatu = 0.0

ALLOCATE (bdragcoefatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bdragcoefatv',2,(/ncloc,nrloc/),kndrtype)
bdragcoefatv = 0.0

ALLOCATE (bcuratc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcuratc',2,(/ncloc,nrloc/),kndrtype)
bcuratc = 0.0

ALLOCATE (bcuratu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcuratu',2,(/ncloc,nrloc/),kndrtype)
bcuratu = 0.0

ALLOCATE (bcuratv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bcuratv',2,(/ncloc,nrloc/),kndrtype)
bcuratv = 0.0

ALLOCATE (bfricatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bfricatu',2,(/ncloc,nrloc/),kndrtype)
bfricatu = 0.0

ALLOCATE (bfricatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bfricatv',2,(/ncloc,nrloc/),kndrtype)
bfricatv = 0.0

ALLOCATE (bstresatc(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('bstresatc',2,(/ncloc+1,nrloc+1/),kndrtype)
bstresatc = 0.0

ALLOCATE (bstresatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bstresatu',2,(/ncloc,nrloc/),kndrtype)
bstresatu = 0.0

ALLOCATE (bstresatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bstresatv',2,(/ncloc,nrloc/),kndrtype)
bstresatv = 0.0

ALLOCATE (ubstresatu(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('ubstresatu',2,(/ncloc+1,nrloc/),kndrtype)
ubstresatu = 0.0

ALLOCATE (vbstresatv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('vbstresatv',2,(/ncloc,nrloc+1/),kndrtype)
vbstresatv = 0.0

ALLOCATE (zroughatc(0:ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('zroughatc',2,(/ncloc+1,nrloc+1/),kndrtype)
zroughatc = 0.0

ALLOCATE (zroughatu(ncloc,nrloc),STAT=errstat)
CALL error_alloc('zroughatu',2,(/ncloc,nrloc/),kndrtype)
zroughatu = 0.0

ALLOCATE (zroughatv(ncloc,nrloc),STAT=errstat)
CALL error_alloc('zroughatv',2,(/ncloc,nrloc/),kndrtype)
zroughatv = 0.0

!
!10.2 Surface fluxes
!-------------------
!

ALLOCATE (cds(ncloc,nrloc),STAT=errstat)
CALL error_alloc('cds',2,(/ncloc,nrloc/),kndrtype)
cds = 0.0

ALLOCATE (ces(ncloc,nrloc),STAT=errstat)
CALL error_alloc('ces',2,(/ncloc,nrloc/),kndrtype)
ces = 0.0

ALLOCATE (chs(ncloc,nrloc),STAT=errstat)
CALL error_alloc('chs',2,(/ncloc,nrloc/),kndrtype)
chs = 0.0

ALLOCATE (qlatent(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qlatent',2,(/ncloc,nrloc/),kndrtype)
qlatent = 0.0

ALLOCATE (qlwave(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qlwave',2,(/ncloc,nrloc/),kndrtype)
qlwave = 0.0

ALLOCATE (qnonsol(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qnonsol',2,(/ncloc,nrloc/),kndrtype)
qnonsol = 0.0

ALLOCATE (qsensible(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qsensible',2,(/ncloc,nrloc/),kndrtype)
qsensible = 0.0

ALLOCATE (ssalflux(ncloc,nrloc),STAT=errstat)
CALL error_alloc('ssalflux',2,(/ncloc,nrloc/),kndrtype)
ssalflux = 0.0

ALLOCATE (sstresatc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('sstresatc',2,(/ncloc,nrloc/),kndrtype)
sstresatc = 0.0

ALLOCATE (usstresatc(0:ncloc,nrloc),STAT=errstat)
CALL error_alloc('usstresatc',2,(/ncloc+1,nrloc/),kndrtype)
usstresatc = 0.0

ALLOCATE (usstresatu(ncloc+1,nrloc),STAT=errstat)
CALL error_alloc('usstresatu',2,(/ncloc+1,nrloc/),kndrtype)
usstresatu = 0.0

ALLOCATE (vsstresatc(ncloc,0:nrloc),STAT=errstat)
CALL error_alloc('vsstresatc',2,(/ncloc,nrloc+1/),kndrtype)
vsstresatc = 0.0

ALLOCATE (vsstresatv(ncloc,nrloc+1),STAT=errstat)
CALL error_alloc('vsstresatv',2,(/ncloc,nrloc+1/),kndrtype)
vsstresatv = 0.0

ALLOCATE (zeros2d(ncloc,nrloc),STAT=errstat)
CALL error_alloc('zeros2d',2,(/ncloc,nrloc/),kndrtype)
zeros2d = 0.0

!
!10.3 Bottom fluxes with waves
!-----------------------------
!

IF (iopt_waves.GT.0) THEN

   ALLOCATE (bstresatc_max(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bstresatc_max',2,(/ncloc,nrloc/),kndrtype)
   bstresatc_max = 0.0

   ALLOCATE (bstresatc_wav(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('bstresatc_wav',2,(/ncloc,nrloc/),kndrtype)
   bstresatc_wav = 0.0

   ALLOCATE (fwave(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('fwave',2,(/ncloc,nrloc/),kndrtype)
   fwave = 0.0

   ALLOCATE (wavethickatc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('wavethickatc',2,(/ncloc,nrloc/),kndrtype)
   wavethickatc = 0.0

   ALLOCATE (zaroughatc(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('zaroughatc',2,(/ncloc,nrloc/),kndrtype)
   zaroughatc = 0.0

ENDIF

!
!11. Surface wave arrays
!-----------------------
!

IF (iopt_waves.GT.0) THEN

!  ---wave parameters
   ALLOCATE (wavedir(0:ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('wavedir',2,(/ncloc+1,nrloc+1/),kndrtype)
   wavedir = 0.0

   ALLOCATE (waveexcurs(0:ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('waveexcurs',2,(/ncloc+1,nrloc+1/),kndrtype)
   waveexcurs = 0.0

   ALLOCATE (wavefreq(0:ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('wavefreq',2,(/ncloc+1,nrloc+1/),kndrtype)
   wavefreq = 0.0

   ALLOCATE (waveheight(0:ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('waveheight',2,(/ncloc+1,nrloc+1/),kndrtype)
   waveheight = 0.0

   ALLOCATE (wavenum(0:ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('wavenum',2,(/ncloc+1,nrloc+1/),kndrtype)
   wavenum = 0.0

   ALLOCATE (waveperiod(0:ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('waveperiod',2,(/ncloc+1,nrloc+1/),kndrtype)
   waveperiod = 0.0

   ALLOCATE (wavevel(0:ncloc,0:nrloc),STAT=errstat)
   CALL error_alloc('wavevel',2,(/ncloc+1,nrloc+1/),kndrtype)
   wavevel = 0.0

   IF (iopt_waves_curr.EQ.1) THEN

!     ---Stokes drift arrays
      ALLOCATE (stokessource2du(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('stokessource2du',2,(/ncloc,nrloc/),kndrtype)

      ALLOCATE (stokessource2dv(ncloc,nrloc),STAT=errstat)
      CALL error_alloc('stokessource2dv',2,(/ncloc,nrloc/),kndrtype)

      ALLOCATE (udstokesatu(ncloc+1,nrloc),STAT=errstat)
      CALL error_alloc('udstokesatu',2,(/ncloc+2,nrloc/),kndrtype)
      udstokesatu = 0.0

      ALLOCATE (umstokesatc(1-nhalo:ncloc+nhalo,0:nrloc),STAT=errstat)
      CALL error_alloc('umstokesatc',2,(/ncloc+2*nhalo,nrloc+1/),kndrtype)
      umstokesatc = 0.0

      ALLOCATE (umstokesatu(1-nhalo:ncloc+nhalo,0:nrloc),STAT=errstat)
      CALL error_alloc('umstokesatu',2,(/ncloc+2*nhalo,nrloc+1/),kndrtype)
      umstokesatu = 0.0

      ALLOCATE (ustokesatc(1-nhalo:ncloc+nhalo,0:nrloc,nz),STAT=errstat)
      CALL error_alloc('ustokesatc',3,(/ncloc+2*nhalo,nrloc+1,nz/),kndrtype)
      ustokesatc = 0.0

      ALLOCATE (ustokesatu(1-nhalo:ncloc+nhalo,0:nrloc,nz),STAT=errstat)
      CALL error_alloc('ustokesatu',3,(/ncloc+2*nhalo,nrloc+1,nz/),kndrtype)
      ustokesatu = 0.0

      ALLOCATE (vdstokesatv(ncloc,nrloc+1),STAT=errstat)
      CALL error_alloc('vdstokesatv',2,(/ncloc,nrloc+1/),kndrtype)
      vdstokesatv = 0.0

      ALLOCATE (vmstokesatc(0:ncloc,1-nhalo:nrloc+nhalo),STAT=errstat)
      CALL error_alloc('vmstokesatc',2,(/ncloc+1,nrloc+2*nhalo/),kndrtype)
      vmstokesatc = 0.0

      ALLOCATE (vmstokesatv(0:ncloc,1-nhalo:nrloc+nhalo),STAT=errstat)
      CALL error_alloc('vmstokesatv',2,(/ncloc+1,nrloc+2*nhalo/),kndrtype)
      vmstokesatv = 0.0

      ALLOCATE (vstokesatc(0:ncloc,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
      CALL error_alloc('vstokesatc',3,(/ncloc+1,nrloc+2*nhalo,nz/),kndrtype)
      vstokesatc = 0.0

      ALLOCATE (vstokesatv(0:ncloc,1-nhalo:nrloc+nhalo,nz),STAT=errstat)
      CALL error_alloc('vstokesatv',3,(/ncloc+1,nrloc+2*nhalo,nz/),kndrtype)
      vstokesatv = 0.0

      ALLOCATE (wstokesatw(0:ncloc,0:nrloc,nz+1),STAT=errstat)
      CALL error_alloc('wstokesatw',3,(/ncloc+1,nrloc+1,nz+1/),kndrtype)
      wstokesatw = 0.0

!     ---wave induced pressure
      ALLOCATE (wavepres(0:ncloc,0:nrloc),STAT=errstat)
      CALL error_alloc('wavepres',2,(/ncloc+1,nrloc+1/),kndrtype)
      wavepres = 0.0

!     ---wave dissipation force
      ALLOCATE (ubwdissipatc(0:ncloc,nrloc,nz),STAT=errstat)
      CALL error_alloc('ubwdisipatc',3,(/ncloc+1,nrloc,nz/),kndrtype)
      ubwdissipatc = 0.0

      ALLOCATE (umbwdissipatc(0:ncloc,nrloc),STAT=errstat)
      CALL error_alloc('umbwdisipatc',2,(/ncloc+1,nrloc/),kndrtype)
      umbwdissipatc = 0.0

      ALLOCATE (umswdissipatc(0:ncloc,nrloc),STAT=errstat)
      CALL error_alloc('umswdisipatc',2,(/ncloc+1,nrloc/),kndrtype)
      umswdissipatc = 0.0

      ALLOCATE (uswdissipatc(0:ncloc,nrloc,nz),STAT=errstat)
      CALL error_alloc('uswdisipatu',3,(/ncloc+1,nrloc,nz/),kndrtype)
      uswdissipatc = 0.0

      ALLOCATE (vbwdissipatc(ncloc,0:nrloc,nz),STAT=errstat)
      CALL error_alloc('vbwdisipatc',3,(/ncloc,nrloc+1,nz/),kndrtype)
      vbwdissipatc = 0.0

      ALLOCATE (vmbwdissipatc(ncloc,0:nrloc),STAT=errstat)
      CALL error_alloc('vmbwdisipatc',2,(/ncloc,nrloc+1/),kndrtype)
      vmbwdissipatc = 0.0

      ALLOCATE (vmswdissipatc(ncloc,0:nrloc),STAT=errstat)
      CALL error_alloc('vmswdisipatc',2,(/ncloc,nrloc+1/),kndrtype)
      vmswdissipatc = 0.0

      ALLOCATE (vswdissipatc(ncloc,0:nrloc,nz),STAT=errstat)
      CALL error_alloc('vswdisipatc',3,(/ncloc,nrloc+1,nz/),kndrtype)
      vswdissipatc = 0.0

   ENDIF

ENDIF

!
!12. Optical arrays
!------------------
!

ALLOCATE (qrad(ncloc,nrloc),STAT=errstat)
CALL error_alloc('qrad',2,(/ncloc,nrloc/),kndrtype)
qrad = 0.0

ALLOCATE (optattcoef2(ncloc,nrloc),STAT=errstat)
CALL error_alloc('optattcoef2',2,(/ncloc,nrloc/),kndrtype)
optattcoef2 = 0.0

ALLOCATE (radiance(ncloc,nrloc,nz+1),STAT=errstat)
CALL error_alloc('radiance',3,(/ncloc,nrloc,nz+1/),kndrtype)
radiance = 0.0

!
!13. Open boundary arrays
!------------------------
!

IF (iopt_obc_sal.EQ.1) THEN
   ALLOCATE (obcsalatu(nobu,nz,0:2,1),STAT=errstat)
   CALL error_alloc('obcsalatu',4,(/nobu,nz,3,1/),kndrtype)
   IF (nobu.GT.0) obcsalatu = 0.0
   ALLOCATE (obcsalatv(nobv,nz,0:2,1),STAT=errstat)
   CALL error_alloc('obcsalatv',4,(/nobv,nz,3,1/),kndrtype)
   IF (nobv.GT.0) obcsalatv = 0.0
ENDIF

IF (iopt_obc_temp.EQ.1) THEN
   ALLOCATE (obctmpatu(nobu,nz,0:2,1),STAT=errstat)
   CALL error_alloc('obctmpatu',4,(/nobu,nz,3,1/),kndrtype)
   IF (nobu.GT.0) obctmpatu = 0.0
   ALLOCATE (obctmpatv(nobv,nz,0:2,1),STAT=errstat)
   CALL error_alloc('obctmpatv',4,(/nobv,nz,3,1/),kndrtype)
   IF (nobv.GT.0) obctmpatv = 0.0
ENDIF

IF (iopt_obc_2D.EQ.1) THEN
   ALLOCATE (obc2uvatu(nobu,2),STAT=errstat)
   CALL error_alloc('obc2uvatu',2,(/nobu,2/),kndrtype)
   IF (nobu.GT.0) obc2uvatu = 0.0
   ALLOCATE (obc2uvatu_old(nobu,2),STAT=errstat)
   CALL error_alloc('obc2uvatu_old',2,(/nobu,2/),kndrtype)
   IF (nobu.GT.0) obc2uvatu_old = 0.0
   ALLOCATE (obc2uvatv(nobv,2),STAT=errstat)
   CALL error_alloc('obc2uvatv',2,(/nobv,2/),kndrtype)
   IF (nobv.GT.0) obc2uvatv = 0.0
   ALLOCATE (obc2uvatv_old(nobv,2),STAT=errstat)
   CALL error_alloc('obc2uvatv_old',2,(/nobv,2/),kndrtype)
   IF (nobv.GT.0) obc2uvatv_old = 0.0
ENDIF

IF (iopt_obc_3D.EQ.1) THEN
   ALLOCATE (obc3uvatu(nobu,nz,2),STAT=errstat)
   CALL error_alloc('obc3uvatu',3,(/nobu,nz,2/),kndrtype)
   IF (nobu.GT.0) obc3uvatu = 0.0
   ALLOCATE (obc3uvatv(nobv,nz,2),STAT=errstat)
   CALL error_alloc('obc3uvatv',3,(/nobv,nz,2/),kndrtype)
   IF (nobv.GT.0) obc3uvatv = 0.0
ENDIF

IF (iopt_obc_th.EQ.1) THEN
   ALLOCATE (floutobu(nobu,nz),STAT=errstat)
   CALL error_alloc('floutobu',2,(/nobu,nz/),kndrtype)
   IF (nobu.GT.0) floutobu = 0.0
   ALLOCATE (floutobv(nobv,nz),STAT=errstat)
   CALL error_alloc('floutobv',2,(/nobv,nz/),kndrtype)
   IF (nobv.GT.0) floutobv = 0.0
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE allocate_mod_arrays

!========================================================================

SUBROUTINE deallocate_mod_arrays
!************************************************************************
!
! *deallocate_mod_arrays* Deallocate model arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Allocate_Arrays.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - simulation_end
!
!************************************************************************
!
USE currents
USE density
USE depths
USE diffusion
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE nestgrids
USE obconds
USE optics
USE paralpars
USE relaxation
USE structures
USE switches
USE tide
USE timepars
USE turbulence
USE wavevars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'deallocate_mod_arrays'
CALL log_timer_in()

!
!1. Grid arrays
!--------------
!

DEALLOCATE (maskatc_int,rmaskatc,seapoint)
DEALLOCATE (iobu,jobu,iobv,jobv)
DEALLOCATE (iobx,jobx,ioby,joby)
DEALLOCATE (iobuloc,jobuloc,iobvloc,jobvloc)
DEALLOCATE (iobxloc,jobxloc,iobyloc,jobyloc)
DEALLOCATE (westobu,soutobv,westobx,soutoby)
DEALLOCATE (nodeatc,nodeatu,nodeatuw,nodeatv,nodeatvw,nodeatuv,node2du,&
          & node2dv,node2duv)
DEALLOCATE (coriolatu,coriolatv)
DEALLOCATE (gaccatc,gaccatu,gaccatv)
DEALLOCATE (gxlon,gylat,gangleatc,garea)
DEALLOCATE (gxcoord,gxcoordglbatc,gycoord,gycoordglbatc)
DEALLOCATE (gscoordatc,gscoordatw,gscoordatu,gscoordatuvw,gscoordatv,&
          & gscoordatuw,gscoordatvw,gsigcoordatc,gsigcoordatw)
DEALLOCATE (delxatc,delxatu,delxatv,delxatuv)
DEALLOCATE (delyatc,delyatu,delyatv,delyatuv)
DEALLOCATE (delzatc,delzatu,delzatv,delzatuv,delzatw,delzatuw,delzatvw)
IF (iopt_fld.GT.0) DEALLOCATE (alphatc_fld, alphatu_fld,alphatv_fld)
IF (iopt_obc_advrlx.EQ.1) DEALLOCATE(rlxobcatu,rlxobcatv)
IF (ALLOCATED(gxcoordglb)) DEALLOCATE(gxcoordglb)
IF (ALLOCATED(gycoordglb)) DEALLOCATE(gycoordglb)
IF (ALLOCATED(gdelxglb)) DEALLOCATE(gdelxglb)
IF (ALLOCATED(gdelyglb)) DEALLOCATE(gdelyglb)
IF (ALLOCATED(gscoordglb)) DEALLOCATE(gscoordglb)
IF (ALLOCATED(depmeanglb)) DEALLOCATE(depmeanglb)
IF (ALLOCATED(seapointglb)) DEALLOCATE(seapointglb)

!
!2. Water depths
!---------------
!

DEALLOCATE (depmeanatc,depmeanatu,depmeanatv,depmeanatuv)
DEALLOCATE (deptotatc,deptotatu,deptotatv,deptotatuv)
DEALLOCATE (deptotatc_old,deptotatu_old,deptotatv_old)
DEALLOCATE (zeta)
IF (iopt_hydro_impl.EQ.1) DEALLOCATE (dzeta,zeta_old)
IF (iopt_fld.GT.0) DEALLOCATE (deptotatc_err,deptotatc_prev)

!
!3. Velocity arrays
!------------------
!

DEALLOCATE (p2dbcgradatu,p2dbcgradatv,p3dbcgradatu,p3dbcgradatv)
DEALLOCATE (ufvel,uvel,uvel_old,udfvel,udvel,udvel_old,umvel,umvel_old,udevint,&
          & umpred)
DEALLOCATE (vfvel,vvel,vvel_old,vdfvel,vdvel,vdvel_old,vmvel,vmvel_old,vdevint,&
          & vmpred)
DEALLOCATE (wfvel,wvel,wphys)

!
!4. Density arrays
!-----------------
!

DEALLOCATE (beta_sal,beta_temp,dens,sal,temp)

!
!5. Diffusion coefficients
!-------------------------
!

DEALLOCATE (hdifcoef2datc,hdifcoef2datuv)
DEALLOCATE (hdifcoef3datc,hdifcoef3datu,hdifcoef3datuv,hdifcoef3datv)
DEALLOCATE (vdifcoefmom,vdifcoefscal,vdifcoeftke)
DEALLOCATE(kinvisc)
IF (iopt_hdif_scal.GE.2) THEN
   DEALLOCATE (hdifcoef3datw,xslopeatu_geo,yslopeatv_geo,vdifcoefscal_rot)
ENDIF
IF (iopt_hdif_scal.EQ.2) DEALLOCATE (xslopeatw_geo,yslopeatw_geo)
IF (iopt_hdif_scal.EQ.3) THEN
   DEALLOCATE (xslopeatu_siso,xslopeatu_ziso,yslopeatv_siso,yslopeatv_ziso)
ENDIF

!
!6. Turbulence arrays
!--------------------
!

DEALLOCATE (dissip,tke,tke_old,zlmix,buofreq2,shearfreq2)

!
!7. Tidal and open boundary arrays
!---------------------------------
!

DEALLOCATE (fnode_obc,phase_obc)
DEALLOCATE (fnode_astro,phase_astro,fxastro,fyastro)

IF (iopt_obc_2D.EQ.1) THEN
   DEALLOCATE (ityp2dobu,ityp2dobv,ityp2dobx,ityp2doby,iloczobu,iloczobv)
   DEALLOCATE (udatobu_amp,zdatobu_amp,vdatobv_amp,zdatobv_amp)
   DEALLOCATE (udatobu_pha,zdatobu_pha,vdatobv_pha,zdatobv_pha)
   DEALLOCATE (udatobu,zdatobu,zdatobu_old,vdatobv,zdatobv,zdatobv_old)
   IF (maxdatafiles(io_2uvobc,1).GT.1) then
      DEALLOCATE (no2dobuv,index2dobuv,iobc2dtype)
   ENDIF
   DEALLOCATE (iqsecobu,jqsecobv)
ENDIF

IF (iopt_obc_2D_tang.EQ.1) THEN
   DEALLOCATE (vdatobx_amp,udatoby_amp,vdatobx_pha,udatoby_pha)
   DEALLOCATE (vdatobx,udatoby)
   IF (maxdatafiles(io_2xyobc,1).GT.1) THEN
      DEALLOCATE (no2dobxy,index2dobxy)
   ENDIF
ENDIF

IF (iopt_sur_1D.EQ.1) THEN
   DEALLOCATE (gxslope_amp,gyslope_amp,gxslope_pha,gyslope_pha,&
             & zeta_amp,zeta_pha)
ENDIF

!
!8. Relaxation zones
!-------------------
!

IF (iopt_obc_relax.GT.0) THEN
   DEALLOCATE (idirrlx,iposrlx,ityprlx,jposrlx,ncrlx,nrrlx)
   IF (inodesrlx(1).GT.0) DEALLOCATE (indexrlxatc,rlxwghtatc)
   IF (inodesrlx(2).GT.0) DEALLOCATE (indexrlxatuv,rlxwghtatuv)
ENDIF

!
!9. Nested output
!----------------
!
IF (iopt_nests.EQ.1) THEN
   DEALLOCATE (nestcoords,nohnstglbc,nohnstglbu,nohnstglbv,nohnstglbx,&
             & nohnstglby)
   DEALLOCATE (nohnstatc,nohnstatu,nohnstatv,nohnstatx,nohnstaty,novnst)
   DEALLOCATE (lbhnstatc,lbhnstatu,lbhnstatv,lbhnstatx,lbhnstaty)
   DEALLOCATE (nohnstcprocs,nohnstuvprocs,nohnstxyprocs,inst2dtype,indexnstc,&
             & indexnstuv,indexnstxy)
   DEALLOCATE (hnstctoc,hnstctou,hnstctov,hnstutou,hnstvtov,hnstvtox,hnstutoy)
   DEALLOCATE (vnstctoc,vnstutou,vnstvtov,vnstvtox,vnstutoy)
   IF (iopt_sed.GT.0) DEALLOCATE (nosednst,instsed)
   IF (iopt_biolgy.GT.0) DEALLOCATE (nobionst,instbio)
ENDIF

!
!10. Meteorological arrays
!-------------------------
!
DEALLOCATE (atmpres,uwindatc,vwindatc,windatc)
DEALLOCATE (airtemp,cloud_cover,evapminprec,evaporation,precipitation,&
          & qspecdif,relhum,sst,tempdif,vappres_air)
IF (iopt_part_model.GT.0) DEALLOCATE (uwindatc_old,vwindatc_old)

!
!11. Bottom and surface fluxes
!-----------------------------
!

DEALLOCATE (bdragcoefatc,bdragcoefatu,bdragcoefatv,bcuratc,bcuratu,bcuratv,&
          & bfricatu,bfricatv,bstresatc,bstresatu,bstresatv,ubstresatu,&
          & vbstresatv,zroughatc,zroughatu,zroughatv)
DEALLOCATE (cds,ces,chs,qlatent,qlwave,qnonsol,qsensible,ssalflux,sstresatc,&
          & usstresatc,usstresatu,vsstresatc,vsstresatv,zeros2d)
IF (iopt_waves.GT.0) THEN
   DEALLOCATE (bstresatc_max,bstresatc_wav,fwave,wavethickatc,zaroughatc)
ENDIF

!
!12. Surface wave arrays
!-----------------------
!

IF (iopt_waves.GT.0) THEN
   DEALLOCATE (wavedir,waveexcurs,wavefreq,waveheight,wavenum,waveperiod,wavevel)
   IF (iopt_waves_curr.EQ.1) THEN
      DEALLOCATE (stokessource2du,stokessource2dv)
      DEALLOCATE (udstokesatu,umstokesatc,umstokesatu,ustokesatc,ustokesatu)
      DEALLOCATE (vdstokesatv,vmstokesatc,vmstokesatv,vstokesatc,vstokesatv,&
                & wstokesatw,wavepres)
      DEALLOCATE (ubwdissipatc,umbwdissipatc,umswdissipatc,uswdissipatc)
      DEALLOCATE (vbwdissipatc,vmbwdissipatc,vmswdissipatc,vswdissipatc)
   ENDIF
ENDIF

!
!13. Optical arrays
!------------------
!

DEALLOCATE (radiance,optattcoef2,qrad)

!
!14. Open boundary arrays
!------------------------
!

IF (iopt_obc_2D.EQ.1) THEN
   DEALLOCATE (obc2uvatu,obc2uvatu_old,obc2uvatv,obc2uvatv_old)
ENDIF
IF (iopt_obc_sal.EQ.1) DEALLOCATE (obcsalatu,obcsalatv)
IF (iopt_obc_temp.EQ.1) DEALLOCATE (obctmpatu,obctmpatv)
IF (iopt_obc_3D.EQ.1) DEALLOCATE (obc3uvatu,obc3uvatv)
IF (iopt_obc_th.EQ.1) DEALLOCATE (floutobu,floutobv)

!
!15. Paralllel arrays
!--------------------
!

DEALLOCATE (ncprocs,nc1procs,nc2procs,nrprocs,nr1procs,nr2procs)
DEALLOCATE (nobuprocs,nobvprocs,nobxprocs,nobyprocs)
DEALLOCATE (nosbuprocs,nosbvprocs,nrvbuprocs,nrvbvprocs)
DEALLOCATE (comprocs,icoordprocs,jcoordprocs)
IF (ALLOCATED(iddomain)) DEALLOCATE (iddomain)
DEALLOCATE (indexobuprocs,indexobvprocs,indexobxprocs,indexobyprocs)
DEALLOCATE (indexobu,indexobv,indexobx,indexoby)
IF (iopt_weibar.EQ.1) THEN
   DEALLOCATE (indexwbaruprocs,indexwbarvprocs,nowbaruprocs,nowbarvprocs)
ENDIF
IF (iopt_dischr.EQ.1) DEALLOCATE (indexdisprocs,nodisprocs)


!
!16. Structures
!--------------
!
!---dry cells
IF (iopt_drycel.EQ.1) DEALLOCATE (idry,jdry)

!---thin dams
IF (iopt_thndam.EQ.1) THEN
   DEALLOCATE (ithinu,ithinuloc,ithinv,ithinvloc,jthinu,jthinuloc,&
             & jthinv,jthinvloc)
ENDIF

!---weirs/barriers
IF (iopt_weibar.EQ.1) THEN
   DEALLOCATE (indexwbaru,indexwbarv,iwbaru,iwbaruloc,iwbarv,iwbarvloc,jwbaru,&
             & jwbaruloc,jwbarv,jwbarvloc,oricoefu,oricoefv,oriheightu,&
             & oriheightv,orisillu,orisillv,wbarcoefu,wbarcoefv,wbarcrestu,&
             & wbarcrestv,wbarmodlu,wbarmodlv,wbarelossu,wbarelossv)
ENDIF

!---discharges
IF (iopt_dischr.EQ.1) THEN
   DEALLOCATE (disarea,disdir,disflag,disspeed,disvol,idis,indexdisloc,&
             & jdis,kdis,kdistype,xdiscoord,ydiscoord,zdiscoord)
   IF (ALLOCATED(idisloc)) DEALLOCATE (idisloc,jdisloc)
ENDIF

!
!17. User output
!---------------
!

!---time series
IF (iopt_out_tsers.EQ.1) THEN
   DEALLOCATE (tsrvars,tsrgpars)
   DEALLOCATE (tsr0d,tsr2d,tsr3d)
   DEALLOCATE (lstatstsr,tsrstatlocs,ivarstsr)
ENDIF

!---time averages
IF (iopt_out_avrgd.EQ.1) THEN
   DEALLOCATE (avrvars,avrgpars)
   DEALLOCATE (avr0d,avr2d,avr3d)
   DEALLOCATE (lstatsavr,avrstatlocs,ivarsavr)
ENDIF

!---harmonic analysis
IF (iopt_out_anal.EQ.1) THEN
   DEALLOCATE (analvars,ellvars,analgpars)
   DEALLOCATE (res0d,res2d,res3d,amp0d,amp2d,amp3d,pha0d,pha2d,pha3d,ell2d,&
             & ell3d)
   DEALLOCATE (lstatsanal,analstatlocs,ivarsanal,ivarsell,ivecell2d,ivecell3d)
   DEALLOCATE (nofreqsharm,ifreqsharm,index_anal,icanal)
   DEALLOCATE (cdate_time_ref,harm_freq_names,harm_freq)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE deallocate_mod_arrays
