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
! *Particle_Forcing_Data* Forcing data for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Forcing_Data.f90  V2.11
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - 
!
! Reference -
!
! Routines - combine_particle_grid, combine_particle_phys_data,
!            define_particle_phys_data, read_particle_grid,
!            read_particle_phys_data, update_particle_phys_data,
!            write_particle_grid, write_particle_phys_data
!
!************************************************************************
!

!=======================================================================

SUBROUTINE combine_particle_grid
!************************************************************************
!
! *combine_particle_grid* Send model grid data from COHERENS to the
!                         particle module
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Forcing_Data.f90  V2.11
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - physical data are combined on the COHERENS master and received
!               by the particle master process
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
! External calls -
!
! Module calls - combine_mod, comms_recv_real, comms_send_real, error_alloc
!
!************************************************************************
!
USE depths
USE grid  
USE gridpars  
USE iopars
USE modids
USE paralpars
USE partids
USE partvars
USE physpars
USE switches
USE syspars
USE comms_MPI, ONLY: comms_recv_real, comms_send_real
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: combine_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!

INTEGER, DIMENSION(3) :: lbounds
REAL, DIMENSION(2) :: rval
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: logglb, logloc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realglb


procname(pglev+1) = 'combine_particle_grid'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE (realglb(nc,nr,nz+1),STAT=errstat)
CALL error_alloc('realglb',3,(/nc,nr,nz+1/),kndrtype)
realglb = 0.0

ALLOCATE (logglb(nc,nr),STAT=errstat)
CALL error_alloc('logglb',2,(/nc,nr/),kndlog)
logglb = .FALSE.

ALLOCATE (logloc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('logloc',2,(/ncloc,nrloc/),kndlog)
logloc = .FALSE.

lbounds = 1

!
!2. Combine data to the particle module
!--------------------------------------
!
!---coordinate arrays
CALL combine_mod(realglb(:,:,1),gxcoord(1:ncloc,1:nrloc),lbounds(1:2),&
               & iarr_gxcoord,0.0,idroot=idmasterpart,comm=comm_world_MPI,&
               & commall=.FALSE.)
IF (idlocglb.EQ.idmasterpart) gxcoord = realglb(:,:,1)
CALL combine_mod(realglb(:,:,1),gycoord(1:ncloc,1:nrloc),lbounds(1:2),&
               & iarr_gycoord,0.0,idroot=idmasterpart,comm=comm_world_MPI,&
               & commall=.FALSE.)
IF (modelid.EQ.modelidpart) gycoord = realglb(:,:,1)

!---sigma coordinates
IF (iopt_grid_nodim.EQ.3) THEN
   IF (iopt_grid_vtype.LT.3) THEN
      IF (idlocglb.EQ.idmastercoh) THEN
         CALL comms_send_real(gsigcoordatw,nz+1,idmasterpart,1,comm_world_MPI,&
                            & 1,iarr_gsigcoordatw)
      ELSEIF (modelid.EQ.modelidpart) THEN
         CALL comms_recv_real(gsigcoordatw,nz+1,idmastercoh,1,comm_world_MPI,&
                            & iarr_gsigcoordatw)
         gsigcoordatc = 0.5*(gsigcoordatw(1:nz)+gsigcoordatw(2:nz+1))
      ENDIF
   ELSE
      CALL combine_mod(realglb(:,:,1:nz),gscoordatc(1:ncloc,1:nrloc,:),&
                     & lbounds,iarr_gscoordatc,0.0,idroot=idmasterpart,&
                     & comm=comm_world_MPI,commall=.FALSE.)
      IF (modelid.EQ.modelidpart) gscoordatc = realglb(:,:,1:nz)
      CALL combine_mod(realglb(:,:,1:nz),gscoordatu(1:ncloc,1:nrloc,:),&
                     & lbounds,iarr_gscoordatu,0.0,idroot=idmasterpart,&
                     & comm=comm_world_MPI,commall=.FALSE.)
      IF (modelid.EQ.modelidpart)  gscoordatu = realglb(:,:,1:nz)
      CALL combine_mod(realglb(:,:,:nz),gscoordatv(1:ncloc,1:nrloc,:),&
                     & lbounds,iarr_gscoordatv,0.0,idroot=idmasterpart,&
                     & comm=comm_world_MPI,commall=.FALSE.)
      IF (modelid.EQ.modelidpart) gscoordatv = realglb(:,:,1:nz)
      CALL combine_mod(realglb,gscoordatw(1:ncloc,1:nrloc,:),lbounds,&
                     & iarr_gscoordatw,0.0,idroot=idmasterpart,&
                     & comm=comm_world_MPI,commall=.FALSE.)
      IF (modelid.EQ.modelidpart)  gscoordatw = realglb
   ENDIF
ENDIF

!---grid spacings
CALL combine_mod(realglb(:,:,1),delxatc(1:ncloc,1:nrloc),&
               & lbounds(1:2),iarr_delxatc,0.0,idroot=idmasterpart,&
               & comm=comm_world_MPI,commall=.FALSE.)
IF (modelid.EQ.modelidpart) delxatc = realglb(:,:,1)
CALL combine_mod(realglb(:,:,1),delyatc(1:ncloc,1:nrloc),&
               & lbounds(1:2),iarr_delyatc,0.0,idroot=idmasterpart,&
               & comm=comm_world_MPI,commall=.FALSE.)
IF (modelid.EQ.modelidpart) delyatc = realglb(:,:,1)

!---bathymetry
CALL combine_mod(realglb(:,:,1),depmeanatc(1:ncloc,1:nrloc),&
               & lbounds(1:2),iarr_depmeanatc,0.0,idroot=idmasterpart,&
               & comm=comm_world_MPI,commall=.FALSE.)
IF (modelid.EQ.modelidpart) depmeanatc(1:nc,1:nr) = realglb(:,:,1)
CALL combine_mod(realglb(:,:,1),depmeanatu(1:ncloc,1:nrloc),&
               & lbounds(1:2),iarr_depmeanatu,0.0,idroot=idmasterpart,&
               & comm=comm_world_MPI,commall=.FALSE.)
IF (modelid.EQ.modelidpart) depmeanatu(1:nc,1:nr) = realglb(:,:,1)
CALL combine_mod(realglb(:,:,1),depmeanatv(1:ncloc,1:nrloc),&
               & lbounds(1:2),iarr_depmeanatv,0.0,idroot=idmasterpart,&
               & comm=comm_world_MPI,commall=.FALSE.)
IF (modelid.EQ.modelidpart) depmeanatv(1:nc,1:nr) = realglb(:,:,1)

!---land masks
IF (modelid.EQ.modelidcoh) logloc = maskatc_int
CALL combine_mod(logglb,logloc,lbounds(1:2),iarr_maskatc,.FALSE.,&
               & idroot=idmasterpart,comm=comm_world_MPI,commall=.FALSE.)
IF (modelid.EQ.modelidpart) maskatc(1:nc,1:nr) = logglb
IF (modelid.EQ.modelidcoh) logloc = node2du(1:ncloc,1:nrloc).GT.0
CALL combine_mod(logglb,logloc,lbounds(1:2),0,.FALSE.,idroot=idmasterpart,&
               & comm=comm_world_MPI,commall=.FALSE.)
IF (modelid.EQ.modelidpart) maskatu(1:nc,1:nr) = logglb
IF (modelid.EQ.modelidcoh) logloc = node2dv(1:ncloc,1:nrloc).GT.0
CALL combine_mod(logglb,logloc,lbounds(1:2),0,.FALSE.,idroot=idmasterpart,&
               & comm=comm_world_MPI,commall=.FALSE.)
IF (modelid.EQ.modelidpart) maskatv(1:nc,1:nr) = logglb

!---model parameters
IF (idlocglb.EQ.idmaster) THEN
   rval(1) = gacc_mean; rval(2) = density_ref
   CALL comms_send_real(rval,2,idmasterpart,1,comm_world_MPI,1,0)
ELSEIF (modelid.EQ.modelidpart) THEN
   CALL comms_recv_real(rval,2,idmastercoh,1,comm_world_MPI,0)
   gacc_mean = rval(1); density_ref = rval(2)
ENDIF

!
!3. Deallocate arrays
!--------------------
!

DEALLOCATE (logglb,logloc,realglb)

CALL log_timer_out()


RETURN

END SUBROUTINE combine_particle_grid

!=======================================================================

SUBROUTINE combine_particle_phys_data
!************************************************************************
!
! *combine_particle_phys_data* Send physical forcing data from COHERENS to the
!                              particle module
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Forcing_Data.f90  V2.11
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - physical data are combined on the COHERENS master and received
!               by the particle master process
!
! Reference -
!
! Calling program - coherens_main, initialise_model, initialise_particle_model,
!                   particle_model
! External calls -
!
! Module calls - combine_mod, error_alloc
!
!************************************************************************
!
USE currents
USE density  
USE depths
USE diffusion
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE partids
USE partswitches
USE partvars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE paral_comms, ONLY: combine_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
  
IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc
INTEGER, DIMENSION(3) :: lbounds
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: logglb, logloc
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realglb, realloc


procname(pglev+1) = 'combine_particle_phys_data'
CALL log_timer_in(npcc)

!
!1. Allocate arrays
!------------------
!

ALLOCATE (realglb(nc,nr,nz+1),STAT=errstat)
CALL error_alloc('realglb',3,(/nc,nr,nz+1/),kndrtype)
realglb = 0.0

IF (iopt_part_vdif.EQ.2) THEN
   ALLOCATE (realloc(ncloc,nrloc,nz+1),STAT=errstat)
   CALL error_alloc('realloc',3,(/ncloc,nrloc,nz+1/),kndrtype)
   realloc = 0.0
ENDIF

IF (iopt_fld.EQ.2) THEN
   ALLOCATE (logglb(nc,nr),STAT=errstat)
   CALL error_alloc('logglb',2,(/nc,nr/),kndlog)
   logglb = .FALSE.
   ALLOCATE (logloc(nc,nr),STAT=errstat)
   CALL error_alloc('logloc',2,(/nc,nr/),kndlog)
   logloc = .FALSE.
ENDIF

lbounds = 1

!
!2. Combine data to the particle module
!--------------------------------------
!
!---water depths
CALL combine_mod(realglb(:,:,1),deptotatc(1:ncloc,1:nrloc),lbounds(1:2),&
               & iarr_deptotatc,0.0,idroot=idmasterpart,comm=comm_world_MPI,&
               & commall=.FALSE.)
IF (modelid.EQ.modelidpart) deptotatc(1:nc,1:nr) = realglb(:,:,1)
CALL combine_mod(realglb(:,:,1),deptotatu(1:ncloc,1:nrloc),lbounds(1:2),&
               & iarr_deptotatu,0.0,idroot=idmasterpart,comm=comm_world_MPI,&
               & commall=.FALSE.)
IF (modelid.EQ.modelidpart) deptotatu = realglb(:,:,1)
CALL combine_mod(realglb(:,:,1),deptotatv(1:ncloc,1:nrloc),lbounds(1:2),&
               & iarr_deptotatv,0.0,idroot=idmasterpart,comm=comm_world_MPI,&
               & commall=.FALSE.)
IF (modelid.EQ.modelidpart) deptotatv = realglb(:,:,1)

!---currents
CALL combine_mod(realglb(:,:,1:nz),uvel(1:ncloc,1:nrloc,:),lbounds,iarr_uvel,&
               & 0.0,idroot=idmasterpart,comm=comm_world_MPI,commall=.FALSE.)
IF (modelid.EQ.modelidpart) uvel(1:nc,1:nr,:) = realglb(:,:,1:nz)
CALL combine_mod(realglb(:,:,1:nz),vvel(1:ncloc,1:nrloc,:),lbounds,iarr_vvel,&
               & 0.0,idroot=idmasterpart,comm=comm_world_MPI,commall=.FALSE.)
IF (modelid.EQ.modelidpart) vvel(1:nc,1:nr,:) = realglb(:,:,1:nz)
IF (iopt_grid_nodim.EQ.3) THEN
   CALL combine_mod(realglb,wvel(1:ncloc,1:nrloc,:),lbounds,iarr_wvel,0.0,&
                  & idroot=idmasterpart,comm=comm_world_MPI,commall=.FALSE.)
   IF (modelid.EQ.modelidpart) wvel(1:nc,1:nr,:) = realglb
ENDIF

!---diffusion coefficient
IF (iopt_part_vdif.EQ.2) THEN
   IF (modelid.EQ.modelidcoh) realloc = vdifcoefmom(1:ncloc,1:nrloc,:)
   CALL combine_mod(realglb,realloc,lbounds,&
                  & iarr_vdifcoefmom,0.0,idroot=idmasterpart,&
                  & comm=comm_world_MPI,commall=.FALSE.)
   IF (modelid.EQ.modelidpart)  vdifcoefpart(1:nc,1:nr,:) = realglb
ENDIF

!---density
IF (iopt_part_dens.EQ.2) THEN
   CALL combine_mod(realglb(:,:,1:nz),dens(1:ncloc,1:nrloc,:),lbounds,&
                  & iarr_dens,0.0,idroot=idmasterpart,comm=comm_world_MPI,&
                  & commall=.FALSE.)
   IF (modelid.EQ.modelidpart)  dens = realglb(:,:,1:nz)
ENDIF

!---land mask
IF (iopt_fld.EQ.2) THEN
   IF (modelid.EQ.modelidcoh) logloc = maskatc_int
   CALL combine_mod(logglb,logloc,lbounds(1:2),iarr_maskatc,.FALSE.,&
                  & idroot=idmasterpart,comm=comm_world_MPI,commall=.FALSE.)
   IF (modelid.EQ.modelidpart) maskatc(1:nc,1:nr) = logglb
   IF (modelid.EQ.modelidcoh) logloc = node2du(1:ncloc,1:nrloc).GT.0
   CALL combine_mod(logglb,logloc,lbounds(1:2),0,.FALSE.,idroot=idmasterpart,&
                  & comm=comm_world_MPI,commall=.FALSE.)
   IF (modelid.EQ.modelidpart) maskatu(1:nc,1:nr) = logglb
   IF (modelid.EQ.modelidcoh) logloc = node2dv(1:ncloc,1:nrloc).GT.0
   CALL combine_mod(logglb,logloc,lbounds(1:2),0,.FALSE.,idroot=idmasterpart,&
                  & comm=comm_world_MPI,commall=.FALSE.)
   IF (modelid.EQ.modelidpart) maskatv(1:nc,1:nr) = logglb
ENDIF

!
!3. Deallocate
!-------------
!

DEALLOCATE (realglb)
IF (iopt_part_vdif.EQ.2) DEALLOCATE (realloc)
IF (iopt_fld.EQ.2) DEALLOCATE (logglb,logloc)

IF (modelid.EQ.modelidpart) THEN
   CALL log_timer_out()
ELSE
   CALL log_timer_out(npcc,itm_part)
ENDIF


RETURN

END SUBROUTINE combine_particle_phys_data

!========================================================================

SUBROUTINE define_particle_phys_data
!************************************************************************
!
! *define_particle_phys_data* Define physical arrays for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Model.f90  V2.11
!
! Description - Only used in online/serial mode (iopt_model_part=1)
!
! Calling program - initialise_model, initialise_particle_model
!
!************************************************************************
!
USE diffusion  
USE grid
USE gridpars  
USE iopars
USE partswitches
USE partvars
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc


procname(pglev+1) = 'define_particle_pĥys_data'
CALL log_timer_in(npcc)

!---set/reset mask at velocity nodes
IF (iopt_fld.EQ.2.OR.nt.EQ.0) THEN
   maskatc(1:nc,1:nr) = maskatc_int 
   maskatu(1:nc,1:nr) = node2du(1:nc,1:nr).GT.0
   maskatv(1:nc,1:nr) = node2dv(1:nc,1:nr).GT.0
ENDIF

!---vertical diffusion coefficient
IF (iopt_part_vdif.EQ.2) THEN
   vdifcoefpart(1:nc,1:nr,:) = vdifcoefmom(1:nc,1:nr,:)
ENDIF

CALL log_timer_out(npcc,itm_part)
  

RETURN

END SUBROUTINE define_particle_phys_data
   
!=======================================================================

SUBROUTINE read_particle_grid
!************************************************************************
!
! *read_particle_grid* read model grid arrays for the particle module
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Forcing_Data.f90  V2.11
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description -
!
! Reference -
!
! Calling program - initialise_particle_model
!
! External calls -
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes  
USE depths
USE grid
USE gridpars
USE iopars
USE partvars
USE physpars
USE switches
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_particle_grid'
CALL log_timer_in()

filepars = modfiles(io_pargrd,1,1)

!
!1. File header
!--------------
!

CALL open_filepars(filepars)
CALL read_glbatts_mod(filepars)
numvars = filepars%novars
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL read_varatts_mod(filepars,varatts,numvars)

!
!2. Coordinate arrays
!--------------------
!

CALL read_vars(gxcoord,filepars,1,varatts=varatts)
CALL read_vars(gycoord,filepars,2,varatts=varatts)
varid = 2

IF (iopt_grid_nodim.EQ.3) THEN
   IF (iopt_grid_vtype.LT.3) THEN
      varid = varid + 1
      CALL read_vars(gsigcoordatw,filepars,varid,varatts=varatts)
      gsigcoordatc = 0.5*(gsigcoordatw(1:nz)+gsigcoordatw(2:nz+1))
   ELSE
      varid = varid + 1
      CALL read_vars(gscoordatc,filepars,varid,varatts=varatts)
      varid = varid + 1
      CALL read_vars(gscoordatu,filepars,varid,varatts=varatts)
      varid = varid + 1
      CALL read_vars(gscoordatv,filepars,varid,varatts=varatts)
      varid = varid + 1
      CALL read_vars(gscoordatw,filepars,varid,varatts=varatts)
   ENDIF
ENDIF

!
!3. Grid spacings
!----------------
!

varid = varid + 1
CALL read_vars(delxatc,filepars,varid,varatts=varatts)
varid = varid + 1
CALL read_vars(delyatc,filepars,varid,varatts=varatts)

!
!4. Bathymetry
!-------------
!

varid = varid + 1
CALL read_vars(depmeanatc(1:nc,1:nr),filepars,varid,varatts=varatts)
varid = varid + 1
CALL read_vars(depmeanatu(1:nc,1:nr),filepars,varid,varatts=varatts)
varid = varid + 1
CALL read_vars(depmeanatv(1:nc,1:nr),filepars,varid,varatts=varatts)

!
!5. Mask array
!-------------
!

varid = varid + 1
CALL read_vars(maskatc(1:nc,1:nr),filepars,varid,varatts=varatts)
varid = varid + 1
CALL read_vars(maskatu(1:nc,1:nr),filepars,varid,varatts=varatts)
varid = varid + 1
CALL read_vars(maskatv(1:nc,1:nr),filepars,varid,varatts=varatts)

!
!6. Model parameters
!-------------------
!

varid = varid + 1
CALL read_vars(gacc_mean,filepars,varid,varatts=varatts)
varid = varid + 1
CALL read_vars(density_ref,filepars,varid,varatts=varatts)

!
!7. Finalise
!-----------
!
 
CALL close_filepars(filepars)
modfiles(io_pargrd,1,1) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE read_particle_grid

!=======================================================================

SUBROUTINE read_particle_phys_data(ciodatetime)
!************************************************************************
!
! *read_particle_phys_data* read physical data for the particle module
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Forcing_Data.f90  V2.11
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description -
!
! Reference -
!
! Calling program - update_particle_phys_data
!
! External calls -
!
! Module calls - error_alloc_struc, open_filepars, read_glbatts_mod, read_time,
!                read_varatts_mod, read_vars
!
!************************************************************************
!
USE currents  
USE datatypes
USE density
USE depths
USE grid
USE gridpars
USE iopars
USE partpars
USE partswitches
USE partvars
USE switches
USE timepars
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: open_filepars, read_glbatts_mod, read_time, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR     Date/time in data file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts

procname(pglev+1) = 'read_particle_phys_data'
CALL log_timer_in()

filepars = modfiles(io_parphs,1,1)

!
!1. File header on first call
!----------------------------
!

IF (filepars%iostat.EQ.0) THEN

!  ---open file
   CALL open_filepars(filepars)

!  ---global attributes   
   CALL read_glbatts_mod(filepars)
   filepars%timerec = MERGE(0,filepars%maxrecs+1,for_drift)
   
!  ---allocate   
   numvars = filepars%novars + 1
   IF (ALLOCATED(varatts)) DEALLOCATE (varatts)
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')

!  ---variable attributes   
   CALL read_varatts_mod(filepars,varatts,numvars)

ENDIF

!
!2. Read date/time
!-----------------
!
!---date/time
CALL read_time(ciodatetime,filepars,back=back_drift)

!---water depths
varid = 2
CALL read_vars(deptotatc(1:nc,1:nr),filepars,varid,varatts=varatts)
varid = varid + 1
CALL read_vars(deptotatu(1:nc,1:nr),filepars,varid,varatts=varatts)
varid = varid + 1
CALL read_vars(deptotatv(1:nc,1:nr),filepars,varid,varatts=varatts)

!---currents
varid = varid + 1
CALL read_vars(uvel(1:nc,1:nr,:),filepars,varid,varatts=varatts)
varid = varid + 1
CALL read_vars(vvel(1:nc,1:nr,:),filepars,varid,varatts=varatts)
IF (iopt_grid_nodim.EQ.3) THEN
   varid = varid + 1
   CALL read_vars(wvel(1:nc,1:nr,:),filepars,varid,varatts=varatts)
ENDIF

!---diffusion coefficient
IF (iopt_part_vdif.EQ.2) THEN
   varid = varid + 1
   CALL read_vars(vdifcoefpart(1:nc,1:nr,:),filepars,varid,varatts=varatts)
ENDIF

!---density
IF (iopt_part_dens.EQ.2) THEN
   varid = varid + 1
   CALL read_vars(dens(1:nc,1:nr,:),filepars,varid,varatts=varatts)
ENDIF

!---land mask
IF (iopt_fld.EQ.2) THEN
   varid = varid + 1
   CALL read_vars(maskatc(1:nc,1:nr),filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(maskatu(1:nc,1:nr),filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(maskatv(1:nc,1:nr),filepars,varid,varatts=varatts)
ENDIF

modfiles(io_parphs,1,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_particle_phys_data
   
!=======================================================================

SUBROUTINE update_particle_phys_data
!************************************************************************
!
! *update_particle_phys_data* update physical data for the particle module
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Forcing_Data.f90  V2.11
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description -
!
! Reference -
!
! Calling program - initialise_particle_model, particle_model
!
! External calls - combine_particle_phys_data, read_particle_phys_data
!
! Module calls - diff_dates, error_alloc, mult_index
!
!************************************************************************
!
USE currents
USE density
USE depths
USE grid
USE gridpars  
USE iopars
USE partpars
USE partswitches
USE partvars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=lentime) :: ciodatetime
INTEGER :: i, j
INTEGER (KIND=kndilong), SAVE :: nosecsdat_old, nosecsdat_new
REAL (KIND=kndrlong) :: ratio1, ratio2, w1, w2
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatc_old, maskatc_new
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatu_old, maskatu_new
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskatv_old, maskatv_new
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: deptotatcdat, deptotatudat, &
                                           & deptotatvdat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: densdat, uveldat, &
                                             & vdifcoefpartdat, vveldat, wveldat


procname(pglev+1) = 'update_particle_phys_data'
CALL log_timer_in()

!
!1. Alllocate on first call
!--------------------------
!

IF (nt.EQ.0.AND.iopt_part_model.GE.3.AND.icpart.GT.0.AND.icpart.NE.1) THEN   
   
   ALLOCATE (deptotatcdat(nc,nr,2),STAT=errstat)
   CALL error_alloc('deptotatcdat',3,(/nc,nr,2/),kndrtype)
   
   ALLOCATE (deptotatudat(nc,nr,2),STAT=errstat)
   CALL error_alloc('deptotatudat',3,(/nc,nr,2/),kndrtype)
   
   ALLOCATE (deptotatvdat(nc,nr,2),STAT=errstat)
   CALL error_alloc('deptotatvdat',3,(/nc,nr,2/),kndrtype)
   
   ALLOCATE (uveldat(nc,nr,nz,2),STAT=errstat)
   CALL error_alloc('uveldat',4,(/nc,nr,nz,2/),kndrtype)

   ALLOCATE (vveldat(nc,nr,nz,2),STAT=errstat)
   CALL error_alloc('vveldat',4,(/nc,nr,nz,2/),kndrtype)

   IF (iopt_grid_nodim.EQ.3) THEN
      ALLOCATE (wveldat(nc,nr,nz+1,2),STAT=errstat)
      CALL error_alloc('wveldat',4,(/nc,nr,nz+1,2/),kndrtype)
   ENDIF

   IF (iopt_part_vdif.EQ.2) THEN
      ALLOCATE (vdifcoefpartdat(nc,nr,nz+1,2),STAT=errstat)
      CALL error_alloc('vdifcoefpartdat',4,(/nc,nr,nz+1,2/),kndrtype)
   ENDIF

   IF (iopt_part_dens.EQ.2) THEN
      ALLOCATE (densdat(nc,nr,nz,2),STAT=errstat)
      CALL error_alloc('densdat',4,(/nc,nr,nz,2/),kndrtype)
   ENDIF

   IF (iopt_fld.EQ.2) THEN
      ALLOCATE (maskatc_old(nc,nr),STAT=errstat)
      CALL error_alloc('maskatc_old',2,(/nc,nr/),log_type)
      ALLOCATE (maskatc_new(nc,nr),STAT=errstat)
      CALL error_alloc('maskatc_new',2,(/nc,nr/),log_type)
      ALLOCATE (maskatu_old(nc,nr),STAT=errstat)
      CALL error_alloc('maskatu_old',2,(/nc,nr/),log_type)
      ALLOCATE (maskatu_new(nc,nr),STAT=errstat)
      CALL error_alloc('maskatu_new',2,(/nc,nr/),log_type)
      ALLOCATE (maskatv_old(nc,nr),STAT=errstat)
      CALL error_alloc('maskatv_old',2,(/nc,nr/),log_type)
      ALLOCATE (maskatv_new(nc,nr),STAT=errstat)
      CALL error_alloc('maskatv_new',2,(/nc,nr/),log_type)
   ENDIF
ENDIF
   
!
!2. Through communication
!------------------------
!

IF (iopt_part_model.EQ.2) THEN

   CALL combine_particle_phys_data

!
!3. From data file without time interpolation
!--------------------------------------------
!   

ELSEIF (iopt_part_model.GE.3.AND.&
     & ((ABS(icpart).EQ.1).OR.&
     &  (icpart.LT.0.AND.(nt.EQ.0.OR.mult_index(nt,icpart))))) THEN
   
   CALL read_particle_phys_data(ciodatetime)
   
!
!4. From data file with time interpolation
!-----------------------------------------
!   

ELSEIF (iopt_part_model.GE.3.AND.icpart.GT.0.AND.icpart.NE.1) THEN

!
!4.1 Read first record
!---------------------
!
   
   IF (nt.EQ.0) THEN
      
!     ---read data
      CALL read_particle_phys_data(ciodatetime)

!     ---date/time
      IF (for_drift) THEN
         CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat_new)
      ELSE
         CALL diff_dates(ciodatetime,CStartDateTime,0,nosecsdat_new)
      ENDIF

!     ---store data
      deptotatcdat(:,:,2) = deptotatc(1:nc,1:nr)
      deptotatudat(:,:,2) = deptotatu
      deptotatvdat(:,:,2) = deptotatv
      uveldat(:,:,:,2) = uvel(1:nc,1:nr,:)
      vveldat(:,:,:,2) = vvel(1:nc,1:nr,:)
      IF (iopt_grid_nodim.EQ.3) wveldat(:,:,:,2) = wvel(1:nc,1:nr,:)
      IF (iopt_part_vdif.EQ.2) THEN
         vdifcoefpartdat(:,:,:,2) = vdifcoefpart(1:nc,1:nr,:)
      ENDIF
      IF (iopt_part_dens.EQ.2) densdat(:,:,:,2) = dens
      IF (iopt_fld.EQ.2) THEN
         maskatc_new = maskatc(1:nc,1:nr)
         maskatu_new = maskatu(1:nc,1:nr)
         maskatv_new = maskatv(1:nc,1:nr)
      ENDIF

   ENDIF

!
!4.2 Read next record
!--------------------
!   

   IF (nt.LT.nstep.AND.(nt.EQ.0.OR.mult_index(nt,icpart))) THEN

!     ---store old data
      nosecsdat_old = nosecsdat_new
      deptotatcdat(:,:,1) = deptotatcdat(:,:,2)
      deptotatudat(:,:,1) = deptotatudat(:,:,2)
      deptotatvdat(:,:,1) = deptotatvdat(:,:,2)
      uveldat(:,:,:,1) = uveldat(:,:,:,2)
      vveldat(:,:,:,1) = vveldat(:,:,:,2)
      IF (iopt_grid_nodim.EQ.3) wveldat(:,:,:,1) = wveldat(:,:,:,2)
      IF (iopt_part_vdif.EQ.2) THEN
         vdifcoefpartdat(:,:,:,1) = vdifcoefpartdat(:,:,:,2)
      ENDIF
      IF (iopt_part_dens.EQ.2) densdat(:,:,:,1) = densdat(:,:,:,2)
      IF (iopt_fld.EQ.2) THEN
         maskatc_old = maskatc_new
         maskatu_old = maskatu_new
         maskatv_old = maskatv_new
      ENDIF
      
!     ---read data
      CALL read_particle_phys_data(ciodatetime)

!     ---date/time
      IF (for_drift) THEN
         CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat_new)
      ELSE
         CALL diff_dates(ciodatetime,CStartDateTime,0,nosecsdat_new)
      ENDIF

!     ---store new data
      deptotatcdat(:,:,2) = deptotatc(1:nc,1:nr)
      deptotatudat(:,:,2) = deptotatu
      deptotatvdat(:,:,2) = deptotatv
      uveldat(:,:,:,2) = uvel(1:nc,1:nr,:)
      vveldat(:,:,:,2) = vvel(1:nc,1:nr,:)
      IF (iopt_grid_nodim.EQ.3) wveldat(:,:,:,2) = wvel(1:nc,1:nr,:)
      IF (iopt_part_vdif.EQ.2) THEN
         vdifcoefpartdat(:,:,:,2) = vdifcoefpart(1:nc,1:nr,:)
      ENDIF
      IF (iopt_part_dens.EQ.2) densdat(:,:,:,2) = dens
      IF (iopt_fld.EQ.2) THEN
         maskatc_new = maskatc(1:nc,1:nr)
         maskatu_new = maskatu(1:nc,1:nr)
         maskatv_new = maskatv(1:nc,1:nr)
      ENDIF

   ENDIF
      
!
!4.3 Interpolate in time
!-----------------------
!
!

   ratio2 = (nosecsrun-nosecsdat_old)/REAL(nosecsdat_new-nosecsdat_old)
   ratio1 = 1.0 - ratio2

!
!4.3.1 Without drying/wetting
!----------------------------
!   

   IF (iopt_fld.LT.2) THEN

      deptotatc(1:ncloc,1:nrloc) = ratio1*deptotatcdat(:,:,1) + &
                                 & ratio2*deptotatcdat(:,:,2)
      deptotatu = ratio1*deptotatudat(:,:,1) + ratio2*deptotatudat(:,:,2)
      deptotatv = ratio1*deptotatvdat(:,:,1) + ratio2*deptotatvdat(:,:,2)
      uvel(1:nc,1:nr,:) = ratio1*uveldat(:,:,:,1) + ratio2*uveldat(:,:,:,2)
      vvel(1:nc,1:nr,:) = ratio1*vveldat(:,:,:,1) + ratio2*vveldat(:,:,:,2)
      IF (iopt_grid_nodim.EQ.3) THEN
         wvel(1:nc,1:nr,:) = ratio1*wveldat(:,:,:,1) + ratio2*wveldat(:,:,:,2)
      ENDIF
      IF (iopt_part_vdif.EQ.2) THEN
         vdifcoefpart(1:nc,1:nr,:) = ratio1*vdifcoefpartdat(:,:,:,1) + &
                      & ratio2*vdifcoefpartdat(:,:,:,2)
      ENDIF
      IF (iopt_part_dens.EQ.2) THEN
         dens = ratio1*densdat(:,:,:,1) + ratio2*densdat(:,:,:,2)
      ENDIF

!
!4.3.2 With drying/wettting
!--------------------------
!      

   ELSE

!      
!4.3.2.1 C-nodes
!---------------
!      

      j_4321: DO j=1,nr
      i_4321: DO i=1,nc
         IF (maskatc_old(i,j).AND.maskatc_new(i,j)) THEN
            maskatc(i,j) = .TRUE.
            w1 = ratio1; w2 = ratio2
         ELSEIF (.NOT.(maskatc_old(i,j).OR.maskatc_new(i,j))) THEN
            maskatc(i,j) = .FALSE.
            w1 = 0.0; w2 = 0.0
         ELSE
            maskatc(i,j) = MERGE(maskatc_old(i,j),maskatc_new(i,j),&
                                   & ratio2.LE.0.5)
            IF (ratio2.LE.0.5) THEN
               w1 = MERGE(1.0,0.0,maskatc_old(i,j))
               w2 = 1.0 - w1
            ELSE
               w2 = MERGE(1.0,0.0,maskatc_new(i,j))
               w2 = 1.0 - w1
            ENDIF
         ENDIF
         deptotatc(i,j) = w1*deptotatcdat(i,j,1) + w2*deptotatcdat(i,j,2)
         IF (iopt_grid_nodim.EQ.3) THEN
            wvel(i,j,:) = w1*wveldat(i,j,:,1) + w2*wveldat(i,j,:,2)
         ENDIF
         IF (iopt_part_vdif.EQ.2) THEN
            vdifcoefpart(i,j,:) = w1*vdifcoefpartdat(i,j,:,1) + &
                                & w2*vdifcoefpartdat(i,j,:,2)
         ENDIF
         IF (iopt_part_dens.EQ.2) THEN
            dens(i,j,:) = w1*densdat(i,j,:,1) + w2*densdat(i,j,:,2)
         ENDIF
      ENDDO i_4321
      ENDDO j_4321

!      
!4.3.2.2 U-nodes
!---------------
!      

      j_4322: DO j=1,nr
      i_4322: DO i=1,nc
         IF (maskatu_old(i,j).AND.maskatu_new(i,j)) THEN
            maskatu(i,j) = .TRUE.
            w1 = ratio1; w2 = ratio2
         ELSEIF (.NOT.(maskatu_old(i,j).OR.maskatu_new(i,j))) THEN
            maskatu(i,j) = .FALSE.
            w1 = 0.0; w2 = 0.0
         ELSE
            maskatu(i,j) = MERGE(maskatu_old(i,j),maskatu_new(i,j),&
                               & ratio2.LE.0.5)
            IF (ratio2.LE.0.5) THEN
               w1 = MERGE(1.0,0.0,maskatu_old(i,j))
               w2 = 1.0 - w1
            ELSE
               w2 = MERGE(1.0,0.0,maskatu_new(i,j))
               w2 = 1.0 - w1
            ENDIF
         ENDIF
         deptotatu(i,j) = w1*deptotatudat(i,j,1) + w2*deptotatudat(i,j,2)
         uvel(i,j,:) = w1*uveldat(i,j,:,1) + w2*uveldat(i,j,:,2)
      ENDDO i_4322
      ENDDO j_4322

!      
!4.3.2.3 V-nodes
!---------------
!      

      j_4323: DO j=1,nr
      i_4323: DO i=1,nc
         IF (maskatv_old(i,j).AND.maskatv_new(i,j)) THEN
            maskatv(i,j) = .TRUE.
            w1 = ratio1; w2 = ratio2
         ELSEIF (.NOT.(maskatv_old(i,j).OR.maskatv_new(i,j))) THEN
            maskatv(i,j) = .FALSE.
            w1 = 0.0; w2 = 0.0
         ELSE
            maskatv(i,j) = MERGE(maskatv_old(i,j),maskatv_new(i,j),&
                               & ratio2.LE.0.5)
            IF (ratio2.LE.0.5) THEN
               w1 = MERGE(1.0,0.0,maskatv_old(i,j))
               w2 = 1.0 - w1
            ELSE
               w2 = MERGE(1.0,0.0,maskatv_new(i,j))
               w2 = 1.0 - w1
            ENDIF
         ENDIF
         deptotatv(i,j) = w1*deptotatvdat(i,j,1) + w2*deptotatvdat(i,j,2)
         vvel(i,j,:) = w1*vveldat(i,j,:,1) + w2*vveldat(i,j,:,2)
      ENDDO i_4323
      ENDDO j_4323
         
   ENDIF

ENDIF

!
!4. Deallocate
!-------------
!

IF ((cold_start.OR.nt.EQ.nstep).AND.iopt_part_model.GE.3.AND.&
  & icpart.GT.0.AND.icpart.NE.1) THEN
   DEALLOCATE (deptotatcdat,deptotatudat,deptotatvdat)
   DEALLOCATE (uveldat,vveldat)
   IF (iopt_grid_nodim.EQ.3) DEALLOCATE (wveldat)
   IF (iopt_part_vdif.EQ.2) DEALLOCATE (vdifcoefpartdat)
   IF (iopt_part_dens.EQ.2) DEALLOCATE (densdat)
   IF (iopt_fld.EQ.2) THEN
      DEALLOCATE (maskatc_old,maskatc_new,maskatu_old,maskatu_new,&
                & maskatv_old,maskatv_new)
   ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE update_particle_phys_data

!=======================================================================

SUBROUTINE write_particle_grid
!************************************************************************
!
! *write_particle_grid* Write model grid arrays for the particle module
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Forcing_Data.f90  V2.11
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - combine_write_mod
!
! Module calls - close_filepars, combine_write_mod, error_alloc,
!                error_alloc_struc, open_filepars, set_modfiles_atts,
!                set_modvars_atts, write_atts_mod, write_vars
!  
!************************************************************************
!
USE datatypes  
USE depths
USE grid
USE gridpars
USE iopars
USE paralpars
USE physpars
USE switches
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE inout_paral, ONLY: combine_write_mod
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars, varid
INTEGER, DIMENSION(3) :: lbounds, ubounds
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: logloc
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'write_particle_grid'
CALL log_timer_in()

!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(io_pargrd,1,2)
filepars = modfiles(io_pargrd,1,2)
numvars = filepars%novars

!---open file
IF (master) CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(io_pargrd,1,2,varatts,numvars)

!---write
IF (master) CALL write_atts_mod(filepars,varatts,numvars)

!---allocate
ALLOCATE (logloc(ncloc,nrloc),STAT=errstat)
CALL error_alloc('logloc',2,(/ncloc,nrloc/),kndlog)

!
!2. Coordinate arrays
!--------------------
!

lbounds(1:2) = 0; ubounds(1:2) = (/nc+1,nr+1/)
CALL combine_write_mod(gxcoord,filepars,1,lbounds(1:2),ubounds(1:2),'UV ',&
                     & varatts=varatts)
CALL combine_write_mod(gycoord,filepars,2,lbounds(1:2),ubounds(1:2),'UV ',&
                     & varatts=varatts)
varid = 2

IF (iopt_grid_nodim.EQ.3) THEN
   IF (iopt_grid_vtype.LT.3) THEN
      varid = varid + 1
      CALL write_vars(gsigcoordatw,filepars,varid,varatts=varatts)
   ELSE
      varid = varid + 1
      lbounds = (/0,0,1/); ubounds = (/nc+1,nr+1,nz/)
      CALL combine_write_mod(gscoordatc,filepars,varid,lbounds,ubounds,'C  ',&
                           & varatts=varatts)
      varid = varid + 1
      CALL combine_write_mod(gscoordatu,filepars,varid,lbounds,ubounds,'U  ',&
                          & varatts=varatts)
      varid = varid + 1
      CALL combine_write_mod(gscoordatv,filepars,varid,lbounds,ubounds,'V  ',&
                           & varatts=varatts)
      varid = varid + 1
      ubounds(3) = nz + 1
      CALL combine_write_mod(gscoordatw,filepars,varid,lbounds,ubounds,'W  ',&
                           & varatts=varatts)
   ENDIF
ENDIF

!
!3. Grid spacings
!----------------
!

varid = varid + 1
lbounds(1:2) = 1-nhalo; ubounds(1:2) = (/nc+nhalo,nr+nhalo/)
CALL combine_write_mod(delxatc,filepars,varid,lbounds(1:2),ubounds(1:2),'C  ',&
                     & varatts=varatts)
varid = varid + 1
CALL combine_write_mod(delyatc,filepars,varid,lbounds(1:2),ubounds(1:2),'C  ',&
                     & varatts=varatts)

!
!4. Bathymetry
!-------------
!

varid = varid + 1
lbounds(1:2) = 1-nhalo; ubounds(1:2) = (/nc+nhalo,nr+nhalo/)
CALL combine_write_mod(depmeanatc,filepars,varid,lbounds(1:2),ubounds(1:2),&
                    & 'C  ',varatts=varatts)
varid = varid + 1
lbounds(1:2) = (/0,1/); ubounds(1:2) = (/nc+1,nr/)
CALL combine_write_mod(depmeanatu,filepars,varid,lbounds(1:2),ubounds(1:2),&
                    & 'U  ',varatts=varatts)
varid = varid + 1
lbounds(1:2) = (/1,0/); ubounds(1:2) = (/nc,nr+1/)
CALL combine_write_mod(depmeanatv,filepars,varid,lbounds(1:2),ubounds(1:2),&
                    & 'V  ',varatts=varatts)

!
!5. Mask arrays
!--------------
!

varid = varid + 1
lbounds(1:2) = 1; ubounds(1:2) = (/nc,nr/)
CALL combine_write_mod(maskatc_int,filepars,varid,lbounds(1:2),ubounds(1:2),&
                    & 'C  ',varatts=varatts)

varid = varid + 1
lbounds(1:2) = 1; ubounds(1:2) = (/nc,nr/)
logloc = node2du(1:ncloc,1:nrloc).GT.0
CALL combine_write_mod(logloc,filepars,varid,lbounds(1:2),ubounds(1:2),&
                    & 'U  ',varatts=varatts)

varid = varid + 1
lbounds(1:2) = 1; ubounds(1:2) = (/nc,nr/)
logloc = node2dv(1:ncloc,1:nrloc).GT.0
CALL combine_write_mod(logloc,filepars,varid,lbounds(1:2),ubounds(1:2),&
                    & 'V  ',varatts=varatts)

!
!6. Model parameters
!-------------------
!

varid = varid + 1
CALL write_vars(gacc_mean,filepars,varid,varatts=varatts)
varid = varid + 1
CALL write_vars(density_ref,filepars,varid,varatts=varatts)

!
!7. Finalise
!-----------
!
 
IF (master) CALL close_filepars(filepars)
modfiles(io_pargrd,1,2) = filepars
DEALLOCATE (varatts,logloc)

CALL log_timer_out()


RETURN

END SUBROUTINE write_particle_grid
   
!=======================================================================

SUBROUTINE write_particle_phys_data(ciodatetime)
!************************************************************************
!
! *write_particle_phys_data* write physical data for the particle module
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Forcing_Data.f90  V2.11
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description -
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! External calls -
!
! Module calls - close_filepars, combine_write_mod, error_alloc_struc, lim_dims,
!                loop_index, open_filepars, set_modfiles_atts, set_modvars_atts,
!                write_atts_mod, write_time
!
!************************************************************************
!
USE currents  
USE datatypes
USE density
USE depths
USE diffusion
USE grid
USE gridpars
USE iopars
USE paralpars
USE partpars
USE partswitches
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc_struc
USE inout_paral, ONLY: combine_write_mod
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, write_time
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: ciodatetime

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR     Date/time in data file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: npcc, numvars, tlimit, varid
INTEGER, DIMENSION(3) :: lbounds, tlims, ubounds
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


filepars = modfiles(io_parphs,1,2)
tlims = filepars%tlims
IF (.NOT.loop_index(tlims,nt)) RETURN

procname(pglev+1) = 'write_particle_phys_data'
CALL log_timer_in(npcc)

!
!1. File header on first call
!----------------------------
!

IF (filepars%iostat.EQ.0) THEN

!  ---file attributes
   CALL set_modfiles_atts(io_parphs,1,2)
   filepars = modfiles(io_parphs,1,2)
   numvars = filepars%novars + 1
   tlimit = lim_dims(tlims)

!  ---open file
   IF (master) CALL open_filepars(filepars)
   
!  ---variable attributes
   IF (ALLOCATED(varatts)) DEALLOCATE (varatts)
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
   CALL set_modvars_atts(io_parphs,1,2,varatts,numvars)

!  ---write attributes
   IF (master) CALL write_atts_mod(filepars,varatts,numvars,tdim=tlimit)

ENDIF

!
!2. Read date/time
!-----------------
!
!---date/time
IF (master) CALL write_time(ciodatetime,filepars,back=back_drift)

!---water depths
varid = 2
lbounds(1:2) = 0; ubounds(1:2) = (/nc+1,nr+1/)
CALL combine_write_mod(deptotatc,filepars,varid,lbounds(1:2),ubounds(1:2),&
                    & 'C  ',varatts=varatts)
varid = varid + 1
lbounds(1:2) = (/1-nhalo,0/); ubounds(1:2) = (/nc+nhalo,nr+1/)
CALL combine_write_mod(deptotatu,filepars,varid,lbounds(1:2),ubounds(1:2),&
                    & 'U  ',varatts=varatts)
varid = varid + 1
lbounds(1:2) = (/0,1-nhalo/); ubounds(1:2) = (/nc+1,nr+nhalo/)
CALL combine_write_mod(deptotatv,filepars,varid,lbounds(1:2),ubounds(1:2),&
                    & 'V  ',varatts=varatts)

!---currents
varid = varid + 1
lbounds = (/1-nhalo,1-nhalo,1/); ubounds = (/nc+nhalo,nr+nhalo,nz/)
CALL combine_write_mod(uvel,filepars,varid,lbounds,ubounds,'U  ',&
                     & varatts=varatts)
varid = varid + 1
CALL combine_write_mod(vvel,filepars,varid,lbounds,ubounds,'V  ',&
                     & varatts=varatts)
IF (iopt_grid_nodim.EQ.3) THEN
   varid = varid + 1
   lbounds = (/0,0,1/); ubounds = (/nc,nr,nz+1/)
   CALL combine_write_mod(wvel,filepars,varid,lbounds,ubounds,'W  ',&
                        & varatts=varatts)
ENDIF   

!---diffusion coefficient
IF (iopt_part_vdif.EQ.2) THEN
   varid = varid + 1
   lbounds = (/0,0,1/); ubounds = (/nc,nr,nz+1/)
   CALL combine_write_mod(vdifcoefmom,filepars,varid,lbounds,ubounds,'W  ',&
                        & varatts=varatts)
ENDIF

!---density
IF (iopt_part_dens.EQ.2) THEN
   varid = varid + 1
   lbounds = (/0,0,1/); ubounds = (/nc+1,nr+1,nz/)
   CALL combine_write_mod(dens,filepars,varid,lbounds,ubounds,'C  ',&
                        & varatts=varatts)
ENDIF

!---land mask
IF (iopt_fld.EQ.2) THEN
   varid = varid + 1
   lbounds(1:2) = 1; ubounds(1:2) = (/nc,nr/)
   CALL combine_write_mod(maskatc_int,filepars,varid,lbounds(1:2),&
                        & ubounds(1:2),'C  ',varatts=varatts)
ENDIF

modfiles(io_parphs,1,2) = filepars

!
!3. Close file
!-------------
!

IF (master.AND.nt.EQ.filepars%tlims(2)) THEN
   CALL close_filepars(filepars)
ENDIF

CALL log_timer_out(npcc,itm_part)


RETURN

END SUBROUTINE write_particle_phys_data
