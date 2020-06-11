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
! *Sediment_Initialisation* Series of routines for initialisation of the
!                           sediment model
!
! Author -  Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Initialisation.f90  V2.11.1
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
! Reference -
!
! Routines - exchange_sedics, initialise_sediment_arrays, read_dar_spec,
!            read_morphics, read_sedics, read_sed_spec
!
!************************************************************************
!


!========================================================================

SUBROUTINE exchange_sedics
!************************************************************************
!
! *exchange_sedics* Exchange initial conditions for the sediments
!
! Author - Boudewijn Decrop and  Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Initialisation.f90  V2.9
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - exchange_mod
!
!************************************************************************
!
USE gridpars
USE iopars
USE sedarrays
USE sedids
USE sedpars
USE sedswitches
USE switches
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER, DIMENSION(4) :: lbounds
INTEGER, DIMENSION(4) :: nhexch


procname(pglev+1) = 'exchange_sedics'
CALL log_timer_in()

!
!1. Transport without flocculation
!---------------------------------
!

IF (iopt_sed.EQ.1) THEN

!
!1.1 Sediment concentrations
!---------------------------
!

   lbounds = (/1-nhalo,1-nhalo,1,1/); nhexch = nhalo
   IF (iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3) THEN
      CALL exchange_mod(cvol,lbounds,nhexch,iarr_cvol)
   ENDIF

!
!1.2 Bed fractions
!-----------------
!

   IF (nf.GT.1) THEN
      lbounds = (/0,0,1,1/); nhexch = (/1,0,1,0/)
      CALL exchange_mod(bed_fraction,lbounds,nhexch,iarr_bed_fraction)
   ENDIF

!
!2.Transport with flocculation
!-----------------------------
!   

ELSEIF (iopt_sed.EQ.2) THEN
   
   lbounds = (/1-nhalo,1-nhalo,1,1/); nhexch = nhalo
   CALL exchange_mod(cnump,lbounds,nhexch,iarr_cnump)
   
ENDIF
   
CALL log_timer_out()


RETURN

END SUBROUTINE exchange_sedics

!========================================================================

SUBROUTINE initialise_sediment_arrays
!************************************************************************
!
! *initialise_sediment_arrays* Initialise arrays for the sediment model
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Sediment_Initialisation.f90  V2.11.1
!
! Description -
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - floc_arrays, median_particle_diameter, sediment roughness
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
USE sedswitches
USE switches
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: k


procname(pglev+1) = 'initialise_sediment_arrays' 
CALL log_timer_in()

!
!1. Particle properties
!----------------------
!

volp = (pi/6.0)*dp**3
massp = rhos*volp

!
!2. Median particle diameter
!---------------------------
!

CALL median_particle_diameter(d50_bed,50.0)

!---exchange halo sections
IF (iopt_MPI.EQ.1) THEN
   CALL exchange_mod(d50_bed,(/0,0/),(/1,0,1,0/),iarr_d50_bed)
ENDIF

!
!3. Bottom roughness
!-------------------
!

IF (iopt_sed_rough.GT.0) CALL sediment_roughness

!
!4. Sediment concentrations
!--------------------------
!
!4.1 Without flocculation
!------------------------
!

IF (iopt_sed.EQ.1) THEN

   k_410: DO k=1,nz
      WHERE (nodeatc(0:ncloc+1,0:nrloc+1).EQ.1)
         ctot(:,:,k) = SUM(cvol(0:ncloc+1,0:nrloc+1,k,:),DIM=3)
      END WHERE
   ENDDO k_410

!
!4.2 With flocculation
!---------------------
!   

ELSEIF (iopt_sed.EQ.2) THEN

   CALL floc_arrays

ENDIF
   
CALL log_timer_out()


RETURN

END SUBROUTINE initialise_sediment_arrays

!========================================================================

SUBROUTINE read_dar_spec
!************************************************************************
!
! *read_dar_spec* Read dredging site and ship properties in standard format
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Sediment_Initialisation.f90  V2.8
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - close_filepars, error_alloc, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE dararrays
USE darpars
USE darswitches
USE datatypes
USE iopars
USE syspars
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iloc, nolocs, numvars, n0
TYPE (FileParams) :: filepars
CHARACTER (LEN=lentime), SAVE, ALLOCATABLE, DIMENSION(:,:) :: coupletime
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: ncoordpol, nosites, notimes
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  fracs, XYpol
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) ::  XYpolsub
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_dar_spec'
CALL log_timer_in()

filepars = modfiles(io_darspc,1,1)

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
!2. Read data
!------------
!

nolocs = MERGE(nocampaigns,nodredgingsites,&
             & iopt_dar_dredging_criterion.EQ.3)
n0 = MERGE(15,0,iopt_dar_dredging_criterion.EQ.3)

!
!2.1 Allocate arrays
!-------------------
!

ALLOCATE (coupletime(nolocs,MaxTimes),STAT=errstat)
CALL error_alloc('coupletime',2,(/nolocs,MaxTimes/),kndchar,lenstr=lentime)

ALLOCATE (ncoordpol(nolocs,MaxSubzones),STAT=errstat)
CALL error_alloc('ncoordpol',2,(/nolocs,MaxSubzones/),kndint)

ALLOCATE (nosites(nolocs,MaxSites),STAT=errstat)
CALL error_alloc('nosites',2,(/nolocs,MaxSites/),kndint)

ALLOCATE (notimes(nolocs,MaxTimes),STAT=errstat)
CALL error_alloc('notimes',2,(/nolocs,MaxTimes/),kndint)

ALLOCATE (fracs(nolocs,MaxFractions,4),STAT=errstat)
CALL error_alloc('fracs',3,(/nolocs,MaxFractions,4/),kndrtype)

ALLOCATE (XYpol(nolocs,MaxCoordPol,2),STAT=errstat)
CALL error_alloc('XYpol',3,(/nolocs,MaxCoordPol,2/),kndrtype)

ALLOCATE (XYpolsub(nolocs,MaxCoordPol,MaxSubzones,2),STAT=errstat)
CALL error_alloc('XYpolsub',4,(/nolocs,MaxCoordPol,MaxSubzones,4/),kndrtype)

!
!2.2 Campaigns
!-------------
!

IF (iopt_dar_dredging_criterion.EQ.3) THEN
   CALL read_vars(Campaigns(1:nolocs)%campdepth,filepars,1,varatts=varatts)
   CALL read_vars(fracs(:,:,1),filepars,2,varatts=varatts)
   CALL read_vars(Campaigns(1:nolocs)%camplevel,filepars,3,varatts=varatts)
   CALL read_vars(Campaigns(1:nolocs)%CampaignStartDateTime,filepars,4,&
                & varatts=varatts)
   CALL read_vars(Campaigns(1:nolocs)%campvolume,filepars,5,varatts=varatts)
   CALL read_vars(coupleTime,filepars,6,varatts=varatts)
   CALL read_vars(Campaigns(1:nolocs)%dredge,filepars,7,varatts=varatts)
   CALL read_vars(Campaigns(1:nolocs)%nr_dredgingsite,filepars,8,&
                & varatts=varatts)
   CALL read_vars(Campaigns(1:nolocs)%nr_relocationsite,filepars,9,&
                & varatts=varatts)
   CALL read_vars(Campaigns(1:nolocs)%nr_ship,filepars,10,varatts=varatts)
   CALL read_vars(Campaigns(1:nolocs)%nr_zones_coupled,filepars,11,&
                & varatts=varatts)
   CALL read_vars(Campaigns(1:nolocs)%relocate,filepars,12,varatts=varatts)
   CALL read_vars(nosites,filepars,13,varatts=varatts)
   CALL read_vars(Campaigns(1:nolocs)%t_lag,filepars,14,varatts=varatts)
   CALL read_vars(Campaigns(1:nolocs)%t_spread,filepars,15,varatts=varatts)
ENDIF

iloc_220: DO iloc=1,nolocs
   Campaigns(iloc)%campfractions = fracs(iloc,:,1)
   Campaigns(iloc)%CoupleTime = coupletime(iloc,:)
   Campaigns(iloc)%relocationsites = nosites(iloc,:)
ENDDO iloc_220

!
!2.3 Dredging sites
!------------------
!

CALL read_vars(coupleTime,filepars,n0+1,varatts=varatts)
CALL read_vars(DredgingSites(1:nolocs)%critical_volume,filepars,n0+2,&
             & varatts=varatts)
CALL read_vars(DredgingSites(1:nolocs)%dredging_thickness,filepars,n0+3,&
             & varatts=varatts)
CALL read_vars(DredgingSites(1:nolocs)%monitoring_frequency,filepars,n0+4,&
             & varatts=varatts)
CALL read_vars(DredgingSites(1:nolocs)%nr_coordpol,filepars,n0+5,&
             & varatts=varatts)
CALL read_vars(notimes,filepars,n0+6,varatts=varatts)
CALL read_vars(DredgingSites(1:nolocs)%nr_ship,filepars,n0+7,varatts=varatts)
CALL read_vars(DredgingSites(1:nolocs)%nr_zones_coupled,filepars,n0+8,&
             & varatts=varatts)
CALL read_vars(DredgingSites(1:nolocs)%over_depth,filepars,n0+9,&
             & varatts=varatts)
CALL read_vars(nosites,filepars,n0+10,varatts=varatts)
CALL read_vars(DredgingSites(1:nolocs)%target_depth,filepars,n0+11,&
             & varatts=varatts)
CALL read_vars(DredgingSites(1:nolocs)%t_lag,filepars,n0+12,varatts=varatts)
CALL read_vars(DredgingSites(1:nolocs)%t_spread,filepars,n0+13,varatts=varatts)
CALL read_vars(XYpol(:,:,1),filepars,n0+14,varatts=varatts)
CALL read_vars(XYpol(:,:,2),filepars,n0+15,varatts=varatts)

iloc_230: DO iloc=1,nolocs
   DredgingSites(iloc)%CoupleTime = coupletime(iloc,:)
   DredgingSites(iloc)%nr_coupled  = notimes(iloc,:)
   DredgingSites(iloc)%relocationsites  = nosites(iloc,:)
   DredgingSites(iloc)%Xpol = XYpol(iloc,:,1)
   DredgingSites(iloc)%Ypol = XYpol(iloc,:,2)
ENDDO iloc_230

!
!2.4 Relocation sites
!--------------------
!

CALL read_vars(RelocationSites(1:nolocs)%max_volume,filepars,n0+16,&
             & varatts=varatts)
CALL read_vars(RelocationSites(1:nolocs)%min_depth,filepars,n0+17,&
             & varatts=varatts)
CALL read_vars(RelocationSites(1:nolocs)%nr_coordpol,filepars,n0+18,&
             & varatts=varatts)
CALL read_vars(XYpol(:,:,1),filepars,n0+19,varatts=varatts)
CALL read_vars(XYpol(:,:,2),filepars,n0+20,varatts=varatts)

iloc_240: DO iloc=1,nolocs
   RelocationSites(iloc)%Xpol = XYpol(iloc,:,1)
   RelocationSites(iloc)%Ypol = XYpol(iloc,:,2)
ENDDO iloc_240

!
!2.5 Ships
!---------
!

CALL read_vars(Ships(1:nolocs)%D_dredge_max,filepars,n0+21,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%D_dredge_min,filepars,n0+22,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%D_relocate_max,filepars,n0+23,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%D_relocate_min,filepars,n0+24,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%H_dredge_max,filepars,n0+25,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%H_dredge_min,filepars,n0+26,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%H_relocate_max,filepars,n0+27,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%H_relocate_min,filepars,n0+28,varatts=varatts)
CALL read_vars(fracs(:,:,1),filepars,n0+29,varatts=varatts)
CALL read_vars(fracs(:,:,2),filepars,n0+30,varatts=varatts)
CALL read_vars(fracs(:,:,3),filepars,31,varatts=varatts)
CALL read_vars(fracs(:,:,4),filepars,n0+32,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%R_circle,filepars,n0+33,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%t_fill,filepars,n0+34,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%t_idle,filepars,n0+35,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%t_relocate,filepars,n0+36,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%t_sail_dredge,filepars,n0+37,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%t_sail_relocate,filepars,n0+38,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%t_wait,filepars,n0+39,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%U_ship,filepars,n0+40,varatts=varatts)
CALL read_vars(Ships(1:nolocs)%V_ship,filepars,n0+41,varatts=varatts)

iloc_250: DO iloc=1,nolocs
   Ships(iloc)%p_overflow  = fracs(iloc,:,1)
   Ships(iloc)%p_passive  = fracs(iloc,:,2)
   Ships(iloc)%p_relocation  = fracs(iloc,:,3)
   Ships(iloc)%p_suctionhead  = fracs(iloc,:,4)
ENDDO iloc_250

!
!2.6 Sub-zone dredging
!---------------------
!

CALL read_vars(ncoordpol,filepars,n0+42,varatts)
CALL read_vars(SubzonesDredging(1:nolocs)%nr_subzones,filepars,n0+43,&
             & varatts=varatts)
CALL read_vars(XYpolsub(:,:,:,1),filepars,n0+44,varatts=varatts)
CALL read_vars(XYpolsub(:,:,:,2),filepars,n0+45,varatts=varatts)

iloc_260: DO iloc=1,nolocs
   SubzonesDredging(iloc)%nr_coordpol = ncoordpol(iloc,:)
   SubzonesDredging(iloc)%Xpol = XYpolsub(iloc,:,:,1)
   SubzonesDredging(iloc)%Ypol = XYpolsub(iloc,:,:,2)
ENDDO iloc_260

!
!2.7  Sub-zone relocation
!------------------------
!

CALL read_vars(ncoordpol,filepars,n0+46,varatts=varatts)
CALL read_vars(SubzonesRelocation(1:nolocs)%nr_subzones,filepars,n0+47,&
             & varatts=varatts)
CALL read_vars(XYpolsub(:,:,:,1),filepars,n0+48,varatts=varatts)
CALL read_vars(XYpolsub(:,:,:,2),filepars,n0+49,varatts=varatts)

iloc_270: DO iloc=1,nolocs
   SubzonesRelocation(iloc)%nr_coordpol = ncoordpol(iloc,:)
   SubzonesRelocation(iloc)%Xpol = XYpolsub(iloc,:,:,1)
   SubzonesRelocation(iloc)%Ypol = XYpolsub(iloc,:,:,2)
ENDDO iloc_270

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_darspc,1,1) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE read_dar_spec

!========================================================================

SUBROUTINE read_morphics
!************************************************************************
!
! *read_morphics* Read initial conditions for the morphological model
!
! Author - Boudewijn Decrop and  Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Initialisation.f90  V2.10
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_distribute_mod, read_glbatts_mod, read_time, 
!                read_varatts_mod, varatts_init
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE morpharrays
USE morphswitches
USE sedpars
USE switches
USE syspars
USE timepars
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_alloc_struc
USE inout_paral, ONLY: read_distribute_mod
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_time, read_varatts_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: fcomm
CHARACTER (LEN=lentime) :: ciodatetime
INTEGER :: ivar, numvars, varid
INTEGER, DIMENSION(3) :: lbounds, ubounds
TYPE(FileParams) :: filepars
INTEGER, DIMENSION(4) :: nhdims
INTEGER, DIMENSION(nf) :: vecids
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_morphics' 
CALL log_timer_in()

filepars = modfiles(io_inicon,ics_morph,1)

!
!1. File header
!--------------
!

CALL open_filepars(filepars)
CALL read_glbatts_mod(filepars)
numvars = filepars%novars + 1
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL varatts_init(varatts)
CALL read_varatts_mod(filepars,varatts,numvars)

!
!2. Read sediment data
!---------------------
!

fcomm = .FALSE.

DO WHILE (.NOT.fcomm)

!
!2.1 Date/time
!-------------
!

   CALL read_time(ciodatetime,filepars)
   fcomm = ciodatetime.EQ.CStartDateTime
   varid = 2
   nhdims = 0

!
!2.2 Morphological arrays
!------------------------
!
!  ---bed layer thickness
   IF (iopt_morph_fixed_layer.EQ.1) THEN
      lbounds = 1; ubounds = (/nc,nr,nb/)
      CALL read_distribute_mod(bed_layer_thickness,filepars,varid,lbounds,&
           & ubounds,nhdims,'C  ',fcomm,varatts=varatts)
      varid = varid + 1
   ENDIF

!  ---bed update
   IF (iopt_morph_time_int.GT.1) THEN
      lbounds = 1; ubounds = (/nc,nr,nf/)
      vecids = (/(ivar,ivar=varid,varid+nf-1)/)     
      CALL read_distribute_mod(bed_update_dep_old1,filepars,0,lbounds,ubounds,&
                             & nhdims,'C  ',fcomm,varatts=varatts,vecids=vecids)
      varid = varid + nf
      vecids = (/(ivar,ivar=varid,varid+nf-1)/)     
      CALL read_distribute_mod(bed_update_ero_old1,filepars,0,lbounds,ubounds,&
                             & nhdims,'C  ',fcomm,varatts=varatts,vecids=vecids)
      varid = varid + nf
      IF (iopt_morph_time_int.EQ.3) THEN
         vecids = (/(ivar,ivar=varid,varid+nf-1)/)     
         CALL read_distribute_mod(bed_update_dep_old2,filepars,0,lbounds,&
                                & ubounds,nhdims,'C  ',fcomm,varatts=varatts,&
                                & vecids=vecids)
         varid = varid + nf
         vecids = (/(ivar,ivar=varid,varid+nf-1)/)
         CALL read_distribute_mod(bed_update_ero_old2,filepars,0,lbounds,&
                                & ubounds,nhdims,'C  ',fcomm,varatts=varatts,&
                                & vecids=vecids)
         varid = varid + nf
      ENDIF
   ENDIF
   lbounds(1:2) = 1; ubounds(1:2) = (/nc,nr/)
   CALL read_distribute_mod(bed_update_int,filepars,varid,lbounds,ubounds,&
                          & nhdims,'C  ',fcomm,varatts=varatts)

ENDDO

!
!3. Close
!--------
!

CALL close_filepars(filepars)

DEALLOCATE (varatts)

modfiles(io_inicon,ics_sed,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_morphics

!========================================================================

SUBROUTINE read_sedics
!************************************************************************
!
! *read_sedics* Read initial conditions for the sediment transport model
!
! Author - Boudewijn Decrop and  Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Initialisation.f90  V2.11.1
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_distribute_mod, read_glbatts_mod, read_time, read_vars,
!                read_varatts_mod, varatts_init
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE sedarrays
USE sedpars
USE sedswitches
USE switches
USE syspars
USE timepars
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_alloc_struc
USE inout_paral, ONLY: read_distribute_mod
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                       & read_time, read_vars, read_varatts_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: fcomm
CHARACTER (LEN=lentime) :: ciodatetime
INTEGER :: ivar, numvars, varid
TYPE(FileParams) :: filepars
INTEGER, DIMENSION(4) :: lbounds, ubounds
INTEGER, DIMENSION(4) :: nhdims
INTEGER, DIMENSION(nf) :: vecids
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_sedics' 
CALL log_timer_in()

filepars = modfiles(io_inicon,ics_sed,1)

!
!1. File header
!--------------
!

CALL open_filepars(filepars)
CALL read_glbatts_mod(filepars)
numvars = filepars%novars + 1
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL varatts_init(varatts)
CALL read_varatts_mod(filepars,varatts,numvars)

!
!2. Read sediment data
!---------------------
!

fcomm = .FALSE.

DO WHILE (.NOT.fcomm)

!
!2.1 Date/time
!-------------
!

   CALL read_time(ciodatetime,filepars)
   fcomm = ciodatetime.EQ.CStartDateTime
   varid = 2

!
!2.2 Sediment concentrations
!---------------------------
!

   lbounds = (/1-nhalo,1-nhalo,1,1/); ubounds = (/nc+nhalo,nr+nhalo,nz,nf/)
   nhdims = nhalo
!  ---sediment concentrations (per fraction)
   IF (iopt_sed.EQ.1) THEN
      IF (iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3) THEN
         vecids = (/(ivar,ivar=varid,varid+nf-1)/)
         CALL read_distribute_mod(cvol,filepars,0,lbounds,ubounds,nhdims,&
                               & 'C  ',fcomm,varatts=varatts,vecids=vecids)
         varid = varid + nf
      ENDIF
   ELSEIF (iopt_sed.EQ.2) THEN
      vecids = (/(ivar,ivar=varid,varid+nf-1)/)
      CALL read_distribute_mod(cnump,filepars,0,lbounds,ubounds,nhdims,&
                            & 'C  ',fcomm,varatts=varatts,vecids=vecids)
      varid = varid + nf
   ENDIF
   
!
!2.3 Bed fractions
!-----------------
!

   IF (nfp.GT.1) THEN
      lbounds = (/0,0,1,1/); ubounds = (/nc,nr,nb,nfp/)
      nhdims = (/1,0,1,0/)
      vecids = (/(ivar,ivar=varid,varid+nfp-1)/)
      CALL read_distribute_mod(bed_fraction,filepars,0,lbounds,ubounds,&
                               nhdims,'C  ',fcomm,varatts=varatts,vecids=vecids)
      varid = varid + nfp
   ENDIF

!
!2.4 Skin roughness
!------------------
!   

   IF (iopt_sed_rough.EQ.2) THEN
      lbounds(1:2) = 0; ubounds(1:2) = (/nc,nr/)
      nhdims = (/1,0,1,0/)
      CALL read_distribute_mod(zroughatc_sed,filepars,varid,&
                             & lbounds(1:2),ubounds(1:2),nhdims,'C  ',fcomm,&
                             & varatts=varatts)
      varid = varid + 1
   ENDIF
      
!
!2.5 Open boundary arrays
!------------------------
!

   IF (iopt_obc_sed.EQ.1) THEN
      IF (nobu.GT.0) THEN
         vecids = (/(ivar,ivar=varid,varid+nf-1)/)
         IF (iopt_sed.EQ.1) THEN
            CALL read_vars(obcsedatu,filepars,0,varatts=varatts,vecids=vecids)
         ELSEIF (iopt_sed.EQ.2) THEN
            CALL read_vars(obccnumpatu,filepars,0,varatts=varatts,vecids=vecids)
         ENDIF
         varid = varid + nf
      ENDIF
      IF (nobv.GT.0) THEN
         vecids = (/(ivar,ivar=varid,varid+nf-1)/)
         IF (iopt_sed.EQ.1) THEN
            CALL read_vars(obcsedatv,filepars,0,varatts=varatts,vecids=vecids)
         ELSEIF (iopt_sed.EQ.2) THEN
            CALL read_vars(obccnumpatv,filepars,0,varatts=varatts,vecids=vecids)
         ENDIF
         varid = varid + nf
      ENDIF
   ENDIF

ENDDO

!
!3. Close
!--------
!

CALL close_filepars(filepars)
DEALLOCATE (varatts)

modfiles(io_inicon,ics_sed,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_sedics

!========================================================================

SUBROUTINE read_sed_spec
!************************************************************************
!
! *read_sed_spec* Read particle properties in standard format
!
! Author - Boudewijn Decrop and  Alexander Breugem (IMDC),
!          Patrick Luyten
!
! Version - @(COHERENS)Sediment_Initialisation.f90  V2.9
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE sedarrays
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_sed_spec'
CALL log_timer_in()

filepars = modfiles(io_sedspc,1,1)

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
!2. Read data
!------------
!

CALL read_vars(dp,filepars,1,varatts=varatts)
CALL read_vars(rhos,filepars,2,varatts=varatts)
CALL read_vars(bstres_cr_cst,filepars,3,varatts=varatts)
CALL read_vars(ws_cst,filepars,4,varatts=varatts)

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_sedspc,1,1) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE read_sed_spec
