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
! *Allocate_Particle_Arrays* Allocate/deallocate arrays for the particle module
!
! Author - Valetrie Duliere
!
! Version - @(COHERENS)Allocate_Particle_Arrays.f90  V2.11
!
! $Date: 2017-10-16 16:32:54 +0200 (Mon, 16 Oct 2017) $
!
! $Revision: 1058 $
!
! Description -
!
! Routines - allocate_part_arrays, deallocate_part_arrays
!
!************************************************************************

!=======================================================================

SUBROUTINE allocate_part_arrays
!************************************************************************
!
! *allocate_part_arrays* Allocate arrays for the particle module
!
! Author - Valerie Duliere
!
! Version - @(COHERENS)Allocate_Particle_Arrays.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_module
!
! Module calls - cloud_init, error_alloc, error_alloc_struc, parts_init
!
!************************************************************************
!
USE currents
USE density
USE depths
USE grid
USE gridpars
USE iopars
USE meteo
USE paralpars
USE partpars
USE partswitches
USE partvars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE parttypes_init, ONLY: cloud_init, parts_init
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: nsize


procname(pglev+1) = 'allocate_part_arrays'
CALL log_timer_in()

!
!1. Coherens model arrays
!------------------------
!

IF (iopt_part_model.GT.2.OR.modelid.EQ.modelidpart) THEN

!  ---coordinate arrays
   ALLOCATE (gxcoord(nc,nr),STAT=errstat)
   CALL error_alloc('gxcoord',2,(/nc,nr/),kndrtype)
   gxcoord = 0.0
   ALLOCATE (gycoord(nc,nr),STAT=errstat)
   CALL error_alloc('gycoordglb',2,(/nc,nr/),kndrtype)
   gycoord = 0.0
   IF (iopt_grid_nodim.NE.2) THEN
      IF (iopt_grid_vtype.LT.3) THEN
         ALLOCATE (gsigcoordatc(nz),STAT=errstat)
         CALL error_alloc('gsigcoordatc',1,(/nz/),kndrtype)
         gsigcoordatc = 0.0
         ALLOCATE (gsigcoordatw(nz+1),STAT=errstat)
         CALL error_alloc('gsigcoordatw',1,(/nz+1/),kndrtype)
         gsigcoordatw = 0.0
      ELSE
         ALLOCATE (gscoordatc(nc,nr,nz),STAT=errstat)
         CALL error_alloc('gscoordatc',3,(/nc,nr,nz/),kndrtype)
         gscoordatc = 0.0
         ALLOCATE (gscoordatu(nc,nr,nz),STAT=errstat)
         CALL error_alloc('gscoordatu',3,(/nc,nr,nz/),kndrtype)
         gscoordatu = 0.0
         ALLOCATE (gscoordatv(nc,nr,nz),STAT=errstat)
         CALL error_alloc('gscoordatv',3,(/nc,nr,nz/),kndrtype)
         gscoordatv = 0.0
         ALLOCATE (gscoordatw(nc,nr,nz+1),STAT=errstat)
         CALL error_alloc('gscoordatw',3,(/nc,nr,nz+1/),kndrtype)
         gscoordatw = 0.0
      ENDIF
   ENDIF
   
!  ---grid spacings
   ALLOCATE (delxatc(nc,nr),STAT=errstat)
   CALL error_alloc('delxatc',2,(/nc,nr/),kndrtype)
   delxatc = 0.0
   ALLOCATE (delyatc(nc,nr),STAT=errstat)
   CALL error_alloc('delyatc',2,(/nc,nr/),kndrtype)
   delyatc = 0.0

   
!  ---bathymetry
   ALLOCATE (depmeanatc(0:nc+1,0:nr+1),STAT=errstat)
   CALL error_alloc('depmeanatc',2,(/nc+2,nr+2/),kndrtype)
   depmeanatc = 0.0
   ALLOCATE (depmeanatu(0:nc+1,0:nr+1),STAT=errstat)
   CALL error_alloc('depmeanatu',2,(/nc+2,nr+2/),kndrtype)
   depmeanatu = 0.0
   ALLOCATE (depmeanatv(0:nc+1,0:nr+1),STAT=errstat)
   CALL error_alloc('depmeanatv',2,(/nc+2,nr+2/),kndrtype)
   depmeanatv = 0.0
   
!  ---total water depths
   ALLOCATE (deptotatc(0:nc+1,0:nr+1),STAT=errstat)
   CALL error_alloc('deptotatc',2,(/nc+2,nr+2/),kndrtype)
   deptotatc = 0.0
   ALLOCATE (deptotatc_old(0:nc+1,0:nr+1),STAT=errstat)
   CALL error_alloc('deptotatc_old',2,(/nc+2,nr+2/),kndrtype)
   deptotatc_old = 0.0
   ALLOCATE (deptotatu(nc,nr),STAT=errstat)
   CALL error_alloc('deptotatu',2,(/nc,nr/),kndrtype)
   deptotatu = 0.0
   ALLOCATE (deptotatv(nc,nr),STAT=errstat)
   CALL error_alloc('deptotatv',2,(/nc,nr/),kndrtype)
   deptotatv = 0.0

!  ---currents
   ALLOCATE (uvel(0:nc+1,0:nr+1,nz),STAT=errstat)
   CALL error_alloc('uvel',3,(/nc+2,nr+2,nz/),kndrtype)
   uvel = 0.0
   ALLOCATE (vvel(0:nc+1,0:nr+1,nz),STAT=errstat)
   CALL error_alloc('vvel',3,(/nc+2,nr+2,nz/),kndrtype)
   vvel = 0.0
   IF (iopt_grid_nodim.EQ.3) THEN
      ALLOCATE (wvel(0:nc+1,0:nr+1,nz+1),STAT=errstat)
      CALL error_alloc('wvel',3,(/nc+2,nr+2,nz+1/),kndrtype)
      wvel = 0.0
   ENDIF
   IF (iopt_part_hadv.EQ.2) THEN
      ALLOCATE (uvel_old(0:nc+1,0:nr+1,nz),STAT=errstat)
      CALL error_alloc('uvel_old',3,(/nc+2,nr+2,nz/),kndrtype)
      uvel_old = 0.0
      ALLOCATE (vvel_old(0:nc+1,0:nr+1,nz),STAT=errstat)
      CALL error_alloc('vvel_old',3,(/nc+2,nr+2,nz/),kndrtype)
      vvel_old = 0.0
   ENDIF
   
!  ---density
   IF (iopt_part_dens.EQ.2) THEN
      ALLOCATE (dens(nc,nr,nz),STAT=errstat)
      CALL error_alloc('dens',3,(/nc,nr,nz/),kndrtype)
      dens = 0.0
   ENDIF
   
!  ---surface winds
   IF (iopt_part_wind.EQ.1) THEN
      ALLOCATE (uwindatc(0:nc+1,0:nr+1),STAT=errstat)
      CALL error_alloc('uwindatc',2,(/nc+2,nr+2/),kndrtype)
      uwindatc = 0.0
      ALLOCATE (uwindatc_old(0:nc+1,0:nr+1),STAT=errstat)
      CALL error_alloc('uwindatc',2,(/nc+2,nr+2/),kndrtype)
      uwindatc_old = 0.0
      ALLOCATE (vwindatc(0:nc+1,0:nr+1),STAT=errstat)
      CALL error_alloc('vwindatc',2,(/nc+2,nr+2/),kndrtype)
      vwindatc = 0.0
      ALLOCATE (vwindatc_old(0:nc+1,0:nr+1),STAT=errstat)
      CALL error_alloc('vwindatc',2,(/nc+2,nr+2/),kndrtype)
      vwindatc_old = 0.0
   ENDIF

ENDIF

!---domain decomposition
IF (iopt_part_model.EQ.2.AND.(modelid.EQ.modelidpart)) THEN
   ALLOCATE(nc1procs(nprocscoh),STAT=errstat)
   CALL error_alloc('nc1procs',1,(/nprocscoh/),kndint)
   ALLOCATE(nc2procs(nprocscoh),STAT=errstat)
   CALL error_alloc('nc2procs',1,(/nprocscoh/),kndint)
   ALLOCATE(nr1procs(nprocscoh),STAT=errstat)
   CALL error_alloc('nr1procs',1,(/nprocscoh/),kndint)
   ALLOCATE(nr2procs(nprocscoh),STAT=errstat)
   CALL error_alloc('nr2procs',1,(/nprocscoh/),kndint)
   ALLOCATE(ncprocs(nprocscoh),STAT=errstat)
   CALL error_alloc('ncprocs',1,(/nprocscoh/),kndint)
   ALLOCATE(nrprocs(nprocscoh),STAT=errstat)
   CALL error_alloc('nrprocs',1,(/nprocscoh/),kndint)
ENDIF

!
!2. Particle model arrays
!------------------------
!

IF (iopt_part_model.NE.2.OR.modelid.EQ.modelidpart) THEN
   
!  ---particle attributes
   ALLOCATE (part(nopart),STAT=errstat)
   CALL error_alloc_struc('part',1,(/nopart/),'PartAtts')
   CALL parts_init(part)

!  ---particle clouds
   ALLOCATE (cloud(noclouds),STAT=errstat)
   CALL error_alloc_struc('cloud',1,(/noclouds/),'PartCloud')
   CALL cloud_init(cloud)

!  ---random numbers for particle diffusion
   nsize = MERGE(0,nopart,iopt_part_hdif.EQ.0.AND.iopt_part_vdif.EQ.0)
   ALLOCATE (ranmagpart(nsize),STAT=errstat)
   CALL error_alloc('ranmagpart',1,(/nsize/),kndrtype)
   nsize = MERGE(0,nopart,iopt_part_hdif.EQ.0)
   ALLOCATE (randirpart(nsize),STAT=errstat)
   CALL error_alloc('randirpart',1,(/nsize/),kndrtype)

!  ---mask arrays
   ALLOCATE (maskatc(0:nc+1,0:nr+1),STAT=errstat)
   CALL error_alloc('maskatc',2,(/nc+2,nr+2/),kndlog)
   maskatc = .FALSE.
   ALLOCATE (maskatu(0:nc+1,0:nr+1),STAT=errstat)
   CALL error_alloc('maskatu',2,(/nc+2,nr+2/),kndlog)
   maskatu = .FALSE.
   ALLOCATE (maskatv(0:nc+1,0:nr+1),STAT=errstat)
   CALL error_alloc('maskatv',2,(/nc+2,nr+2/),kndlog)
   maskatv = .FALSE.

!  ---diffusion
   IF (iopt_part_vdif.EQ.2) THEN
      ALLOCATE (vdifcoefpart(0:nc+1,0:nr+1,nz+1),STAT=errstat)
      CALL error_alloc('vdifcoefpart',3,(/nc+2,nr+2,nz+1/),kndrtype)
      vdifcoefpart = 0.0
   ENDIF

ENDIF

!---particulate matter
IF (iopt_part_conc.GT.1.AND.(iopt_part_model.NE.2.OR.&
   & modelid.EQ.modelidpart)) THEN
   ALLOCATE (diamconc(nolabels),STAT=errstat)
   CALL error_alloc('diamconc',1,(/nolabels/),kndrtype)
   diamconc = 0.0
   ALLOCATE (massconc(nolabels),STAT=errstat)
   CALL error_alloc('massconc',1,(/nolabels/),kndrtype)
   massconc = 0.0
   ALLOCATE (rhoconc(nolabels),STAT=errstat)
   CALL error_alloc('rhoconc',1,(/nolabels/),kndrtype)
   rhoconc = 0.0
   ALLOCATE (volumeconc(nolabels),STAT=errstat)
   CALL error_alloc('volumeconc',1,(/nolabels/),kndrtype)
   volumeconc = 0.0
ENDIF

!---volume distributions
IF (iopt_part_conc.GT.0) THEN
   ALLOCATE (ageconc(ncloc,nrloc,nz,nolabels),STAT=errstat)
   CALL error_alloc('ageconc',4,(/ncloc,nrloc,nz,nolabels/),kndrtype)
   ageconc = 0.0
   ALLOCATE (ageconc_tot(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('ageconc_tot',3,(/ncloc,nrloc,nz/),kndrtype)
   ageconc_tot = 0.0
   ALLOCATE (conc(ncloc,nrloc,nz,nolabels),STAT=errstat)
   CALL error_alloc('conc',4,(/ncloc,nrloc,nz,nolabels/),kndrtype)
   conc = 0.0
   ALLOCATE (conc_tot(ncloc,nrloc,nz),STAT=errstat)
   CALL error_alloc('conc_tot',3,(/ncloc,nrloc,nz/),kndrtype)
   conc_tot = 0.0
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE allocate_part_arrays

!========================================================================

SUBROUTINE deallocate_part_arrays
!************************************************************************
!
! *deallocate_part_arrays* Deallocate arrays for the particle module
!
! Author - Valerie Duliere
!
! Version - @(COHERENS)Allocate_Arrays.f90  V2.11
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
USE grid
USE iopars
USE meteo
USE paralpars
USE partswitches
USE partpars
USE partvars
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'deallocate_part_arrays'
CALL log_timer_in()

!
!1. Coherens model arrays
!------------------------
!

IF (iopt_part_model.GT.2.OR.modelid.EQ.modelidpart) THEN

!  ---coordinate arrays   
   DEALLOCATE (gxcoord,gycoord)
   IF (iopt_grid_nodim.NE.2) THEN
      IF (iopt_grid_vtype.LT.3) THEN
         DEALLOCATE (gsigcoordatc,gsigcoordatw)
      ELSE
         DEALLOCATE (gscoordatc,gscoordatu,gscoordatv,gscoordatw)
      ENDIF
   ENDIF
   
!  ---grid spacings
   DEALLOCATE (delxatc,delyatc)

!  ---bathymetry
   DEALLOCATE (depmeanatc,depmeanatu,depmeanatv)

!  ---total water depths
   DEALLOCATE (deptotatc,deptotatc_old,deptotatu,deptotatv)

!  ---currents
   DEALLOCATE (uvel,vvel)
   IF (iopt_grid_nodim.EQ.3) DEALLOCATE (wvel)
   IF (iopt_part_hadv.EQ.2) DEALLOCATE(uvel_old,vvel_old) 

!  ---density
   IF (iopt_part_dens.EQ.2) DEALLOCATE (dens)   

!  ---surface winds
   IF (iopt_part_wind.EQ.1) THEN
      DEALLOCATE (uwindatc,uwindatc_old,vwindatc,vwindatc_old)
   ENDIF

ENDIF

!---domain decomposition
IF (iopt_part_model.EQ.2.AND.(modelid.EQ.modelidpart)) THEN
   DEALLOCATE (nc1procs,nc2procs,nr1procs,nr2procs,ncprocs,nrprocs)
ENDIF

!
!2. Particle model arrays
!------------------------
!

IF (iopt_part_model.NE.2.OR.modelid.EQ.modelidpart) THEN

!  ---particle attributes
   DEALLOCATE (part)

!  ---particle clouds
   DEALLOCATE (cloud)

!  ---random numbers for particle diffusion
   DEALLOCATE (randirpart,ranmagpart)

!  ---masks at velocity nodes
   DEALLOCATE (maskatc,maskatu,maskatv)

!  ---diffusion
   IF (iopt_part_vdif.EQ.2) DEALLOCATE (vdifcoefpart)

ENDIF

!---particulate matter
IF (iopt_part_conc.GT.1.AND.(iopt_part_model.NE.2.OR.&
  & modelid.EQ.modelidpart)) THEN
   DEALLOCATE (diamconc,massconc,rhoconc,volumeconc)
ENDIF

!---volume distributions
IF (iopt_part_conc.GT.0) THEN
   DEALLOCATE (ageconc,ageconc_tot,conc,conc_tot)
ENDIF

!---particle output
IF (ALLOCATED(part1d)) THEN
   DEALLOCATE (part1d,outpvars,outppars)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE deallocate_part_arrays
