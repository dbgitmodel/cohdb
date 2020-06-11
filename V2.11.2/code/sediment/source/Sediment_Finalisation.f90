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
! *Sediment_Finalisation* Finalise sediment model
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Finalisation.f90  V2.11.1
!
! $Date: 2018-05-23 10:09:17 +0200 (Wed, 23 May 2018) $
!
! $Revision: 1138 $
!
! Description -
!
! Reference -
!
! Routines - write_dar_spec, write_morphics, write_sedics, write_sed_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE write_dar_spec
!************************************************************************
!
! *write_dar_spec* Write dredging/relocation attributes in standard format
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Finalisation.f90  V2.8
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - close_filepars, error_alloc, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_vars
!
!************************************************************************
!
USE dararrays
USE darpars
USE darswitches
USE iopars
USE paralpars
USE syspars
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, &
                        & write_atts_mod, write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
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


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_dar_spec'
CALL log_timer_in()

!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(io_darspc,1,2)
filepars = modfiles(io_darspc,1,2)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(io_darspc,1,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2. Write data
!-------------
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

   iloc_220: DO iloc=1,nolocs
      fracs(iloc,:,1) = Campaigns(iloc)%campfractions
      coupletime(iloc,:) =  Campaigns(iloc)%CoupleTime
      nosites(iloc,:) =  Campaigns(iloc)%relocationsites
   ENDDO iloc_220

   CALL write_vars(Campaigns(1:nolocs)%campdepth,filepars,1,varatts=varatts)
   CALL write_vars(fracs(:,:,1),filepars,2,varatts=varatts)
   CALL write_vars(Campaigns(1:nolocs)%camplevel,filepars,3,varatts=varatts)
   CALL write_vars(Campaigns(1:nolocs)%CampaignStartDateTime,filepars,4,&
                 & varatts=varatts)
   CALL write_vars(Campaigns(1:nolocs)%campvolume,filepars,5,varatts=varatts)
   CALL write_vars(coupletime,filepars,6,varatts=varatts)
   CALL write_vars(Campaigns(1:nolocs)%dredge,filepars,7,varatts=varatts)
   CALL write_vars(Campaigns(1:nolocs)%nr_dredgingsite,filepars,8,&
                 & varatts=varatts)
   CALL write_vars(Campaigns(1:nolocs)%nr_relocationsite,filepars,9,&
                 & varatts=varatts)
   CALL write_vars(Campaigns(1:nolocs)%nr_ship,filepars,10,varatts=varatts)
   CALL write_vars(Campaigns(1:nolocs)%nr_zones_coupled,filepars,11,&
                 & varatts=varatts)
   CALL write_vars(Campaigns(1:nolocs)%relocate,filepars,12,varatts=varatts)
   CALL write_vars(nosites,filepars,13,varatts=varatts)
   CALL write_vars(Campaigns(1:nolocs)%t_lag,filepars,14,varatts=varatts)
   CALL write_vars(Campaigns(1:nolocs)%t_spread,filepars,15,varatts=varatts)

ENDIF

!
!2.3 Dredging sites
!------------------
!

iloc_230: DO iloc=1,nolocs
   coupletime(iloc,:) =  DredgingSites(iloc)%CoupleTime
   notimes(iloc,:) =  DredgingSites(iloc)%nr_coupled
   nosites(iloc,:) =  DredgingSites(iloc)%relocationsites
   XYpol(iloc,:,1)  = DredgingSites(iloc)%Xpol
   XYpol(iloc,:,2)  = DredgingSites(iloc)%Ypol
ENDDO iloc_230

CALL write_vars(coupleTime,filepars,n0+1,varatts=varatts)
CALL write_vars(DredgingSites(1:nolocs)%critical_volume,filepars,n0+2,&
              & varatts=varatts)
CALL write_vars(DredgingSites(1:nolocs)%dredging_thickness,filepars,n0+3,&
              & varatts=varatts)
CALL write_vars(DredgingSites(1:nolocs)%monitoring_frequency,filepars,n0+4,&
              & varatts=varatts)
CALL write_vars(DredgingSites(1:nolocs)%nr_coordpol,filepars,n0+5,&
              & varatts=varatts)
CALL write_vars(notimes,filepars,n0+6,varatts=varatts)
CALL write_vars(DredgingSites(1:nolocs)%nr_ship,filepars,n0+7,varatts=varatts)
CALL write_vars(DredgingSites(1:nolocs)%nr_zones_coupled,filepars,n0+8,&
              & varatts=varatts)
CALL write_vars(DredgingSites(1:nolocs)%over_depth,filepars,n0+9,&
              & varatts=varatts)
CALL write_vars(nosites,filepars,n0+10,varatts=varatts)
CALL write_vars(DredgingSites(1:nolocs)%target_depth,filepars,n0+11,&
              & varatts=varatts)
CALL write_vars(DredgingSites(1:nolocs)%t_lag,filepars,n0+12,varatts=varatts)
CALL write_vars(DredgingSites(1:nolocs)%t_spread,filepars,n0+13,varatts=varatts)
CALL write_vars(XYpol(:,:,1),filepars,n0+14,varatts=varatts)
CALL write_vars(XYpol(:,:,2),filepars,n0+15,varatts=varatts)

!
!2.4 Relocation sites
!--------------------
!

iloc_240: DO iloc=1,nolocs
   XYpol(iloc,:,1)  = RelocationSites(iloc)%Xpol
   XYpol(iloc,:,2)  = RelocationSites(iloc)%Ypol
ENDDO iloc_240

CALL write_vars(RelocationSites(1:nolocs)%max_volume,filepars,n0+16,&
              & varatts=varatts)
CALL write_vars(RelocationSites(1:nolocs)%min_depth,filepars,n0+17,&
              & varatts=varatts)
CALL write_vars(RelocationSites(1:nolocs)%nr_coordpol,filepars,n0+18,&
              & varatts=varatts)
CALL write_vars(XYpol(iloc,:,1),filepars,n0+19,varatts=varatts)
CALL write_vars(XYpol(iloc,:,2),filepars,n0+20,varatts=varatts)

!
!2.5 Ships
!---------
!

iloc_250: DO iloc=1,nolocs
   fracs(iloc,:,1) = Ships(iloc)%p_overflow
   fracs(iloc,:,2) = Ships(iloc)%p_passive
   fracs(iloc,:,3) = Ships(iloc)%p_relocation
   fracs(iloc,:,4) = Ships(iloc)%p_suctionhead
ENDDO iloc_250

CALL write_vars(Ships(1:nolocs)%D_dredge_max,filepars,n0+21,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%D_dredge_min,filepars,n0+22,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%D_relocate_max,filepars,n0+23,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%D_relocate_min,filepars,n0+24,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%H_dredge_max,filepars,n0+25,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%H_dredge_min,filepars,n0+26,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%H_relocate_max,filepars,n0+27,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%H_relocate_min,filepars,n0+28,varatts=varatts)
CALL write_vars(fracs(:,:,1),filepars,n0+29,varatts=varatts)
CALL write_vars(fracs(:,:,2),filepars,n0+30,varatts=varatts)
CALL write_vars(fracs(:,:,3),filepars,n0+31,varatts=varatts)
CALL write_vars(fracs(:,:,4),filepars,n0+32,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%R_circle,filepars,n0+33,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%t_fill,filepars,n0+34,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%t_idle,filepars,n0+35,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%t_relocate,filepars,n0+36,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%t_sail_dredge,filepars,n0+37,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%t_sail_relocate,filepars,n0+38,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%t_wait,filepars,n0+39,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%U_ship,filepars,n0+40,varatts=varatts)
CALL write_vars(Ships(1:nolocs)%V_ship,filepars,n0+41,varatts=varatts)

!
!2.6 Sub-zone dredging
!---------------------
!

iloc_260: DO iloc=1,nolocs
   ncoordpol(iloc,:) = SubzonesDredging(iloc)%nr_coordpol
   XYpolsub(iloc,:,:,1) = SubzonesDredging(iloc)%Xpol
   XYpolsub(iloc,:,:,2) = SubzonesDredging(iloc)%Ypol
ENDDO iloc_260

CALL write_vars(ncoordpol,filepars,n0+42,varatts)
CALL write_vars(SubzonesDredging(1:nolocs)%nr_subzones,filepars,n0+43,&
              & varatts=varatts)
CALL write_vars(XYpolsub(:,:,:,1),filepars,n0+44,varatts=varatts)
CALL write_vars(XYpolsub(:,:,:,2),filepars,n0+45,varatts=varatts)

!
!2.7  Sub-zone relocation
!------------------------
!

iloc_270: DO iloc=1,nolocs
   ncoordpol(iloc,:) = SubzonesRelocation(iloc)%nr_coordpol
   XYpolsub(iloc,:,:,1) = SubzonesRelocation(iloc)%Xpol
   XYpolsub(iloc,:,:,2) = SubzonesRelocation(iloc)%Ypol
ENDDO iloc_270

CALL write_vars(ncoordpol,filepars,n0+46,varatts=varatts)
CALL write_vars(SubzonesRelocation(1:nolocs)%nr_subzones,filepars,n0+47,&
              & varatts=varatts)
CALL write_vars(XYpolsub(:,:,:,1),filepars,n0+48,varatts=varatts)
CALL write_vars(XYpolsub(:,:,:,2),filepars,n0+49,varatts=varatts)

!
!2.8 Deallocate
!--------------
!

DEALLOCATE (coupletime,ncoordpol,nosites,notimes,fracs,XYpol,XYpolsub)

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_darspc,1,2) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE write_dar_spec

!========================================================================

SUBROUTINE write_morphics
!************************************************************************
!
! *write_morphics* Write initial conditions for the morphological model
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Finalisation.f90  V2.11.1
!
! Description - 
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! Module calls - close_filepars, combine_write_mod, error_alloc_struc,
!                open_filepars, set_modfiles_atts, set_modvars_atts,
!                write_atts_mod, write_time
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE morpharrays
USE morphswitches
USE paralpars
USE sedpars
USE timepars
USE error_routines, ONLY: error_alloc_struc
USE inout_paral, ONLY: combine_write_mod
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_time
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
INTEGER :: ivar, numvars, varid
INTEGER, SAVE :: itrestart
TYPE(FileParams) :: filepars
INTEGER, DIMENSION(3) :: lbounds, ubounds
INTEGER, DIMENSION(nf) :: vecids
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


IF (norestarts.EQ.0.OR.modfiles(io_fincon,ics_morph,2)%status.EQ.'0') RETURN
IF (nt.EQ.0) itrestart = 0
IF (itrestart.GT.norestarts) RETURN

procname(pglev+1) = 'write_morphics' 
flag = .FALSE.
filepars = modfiles(io_fincon,ics_morph,2)

!
!1. Write file header on first call
!----------------------------------
!

IF (nt.EQ.ntrestart(1)) THEN

   CALL log_timer_in()
   flag = .TRUE.

!  ---file attributes
   CALL set_modfiles_atts(io_fincon,ics_morph,2)
   filepars = modfiles(io_fincon,ics_morph,2)
   numvars = filepars%novars + 1

!  ---open file
   IF (master) CALL open_filepars(filepars)

!  ---variable attributes
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
   CALL set_modvars_atts(io_fincon,ics_morph,2,varatts,numvars)

!  ---write attributes
   IF (master) CALL write_atts_mod(filepars,varatts,numvars)

   itrestart = 1

ENDIF

!
!2. Write morphological data
!---------------------------
!

IF (itrestart.GT.0.AND.itrestart.LE.norestarts) THEN

   IF (nt.EQ.ntrestart(itrestart)) THEN
      IF (itrestart.GT.1) THEN
         CALL log_timer_in()
         flag = .TRUE.
      ENDIF

!
!2.1 Date/time
!-------------
!

      CALL write_time(CDateTime,filepars)
      varid = 2

!
!2.2 Morphological arrays
!------------------------
!
!     ---bed layer thickness
      IF (iopt_morph_fixed_layer.EQ.1) THEN
         lbounds = 1; ubounds = (/nc,nr,nb/)
         CALL combine_write_mod(bed_layer_thickness,filepars,varid,lbounds,&
              & ubounds,'C  ',varatts=varatts)
         varid = varid + 1
      ENDIF
!     ---bed update
      IF (iopt_morph_time_int.GT.1) THEN
         lbounds = 1; ubounds = (/nc,nr,nf/)
         vecids = (/(ivar,ivar=varid,varid+nf-1)/)
         CALL combine_write_mod(bed_update_dep_old1,filepars,0,lbounds,&
                              & ubounds,'C  ',varatts=varatts,vecids=vecids)
         varid = varid + nf
         vecids = (/(ivar,ivar=varid,varid+nf-1)/)
         CALL combine_write_mod(bed_update_ero_old1,filepars,0,lbounds,&
              & ubounds,'C  ',varatts=varatts,vecids=vecids)
         varid = varid + nf
         IF (iopt_morph_time_int.EQ.3) THEN

            vecids = (/(ivar,ivar=varid,varid+nf-1)/)
            CALL combine_write_mod(bed_update_dep_old2,filepars,0,lbounds,&
                                 & ubounds,'C  ',varatts=varatts,vecids=vecids)
            varid = varid + nf
            vecids = (/(ivar,ivar=varid,varid+nf-1)/)
            CALL combine_write_mod(bed_update_ero_old2,filepars,0,lbounds,&
                 & ubounds,'C  ',varatts=varatts,vecids=vecids)
            varid = varid + nf
         ENDIF
      ENDIF
      lbounds(1:2) = 1; ubounds(1:2) = (/nc,nr/)
      CALL combine_write_mod(bed_update_int,filepars,varid,lbounds(1:2),&
                           & ubounds(1:2),'C  ',varatts=varatts)

!
!2.3 Next output index
!---------------------
!

      itrestart = itrestart + 1

   ENDIF

ENDIF

!
!3. Finalise
!-----------
!

IF ((cold_start.AND.ntrestart(1).EQ.0).OR.&
  & (itrestart.GT.norestarts)) THEN
   IF (master) CALL close_filepars(filepars)
   DEALLOCATE (varatts)
ENDIF

modfiles(io_fincon,ics_morph,2) = filepars

IF (flag) CALL log_timer_out() 


RETURN

END SUBROUTINE write_morphics

!========================================================================

SUBROUTINE write_sedics
!************************************************************************
!
! *write_sedics* Write initial conditions for the sediment model
!
! Author - Boudewijn Decrop and Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Finalisation.f90  V2.11.1
!
! Description - 
!
! Reference -
!
! Calling program - coherens_main, initialise_model
!
! Module calls - close_filepars, combine_write_mod, combine_write_stats_glb,
!                error_alloc_struc, open_filepars, set_modfiles_atts,
!                set_modvars_atts, write_atts_mod, write_time
!
!************************************************************************
!
USE datatypes
USE grid
USE gridpars
USE iopars
USE paralpars
USE sedarrays
USE sedpars
USE sedswitches
USE switches
USE timepars
USE error_routines, ONLY: error_alloc_struc
USE inout_paral, ONLY: combine_write_mod, combine_write_stats_glb
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_time
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
INTEGER :: ivar, numvars, varid
INTEGER, SAVE :: itrestart
TYPE(FileParams) :: filepars
INTEGER, DIMENSION(4) :: lbounds, ubounds
INTEGER, DIMENSION(nf) :: vecids
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


IF (norestarts.EQ.0.OR.modfiles(io_fincon,ics_sed,2)%status.EQ.'0') RETURN
IF (nt.EQ.0) itrestart = 0
IF (itrestart.GT.norestarts) RETURN

procname(pglev+1) = 'write_sedics' 
flag = .FALSE.
filepars = modfiles(io_fincon,ics_sed,2)

!
!1. Write file header on first call
!----------------------------------
!

IF (nt.EQ.ntrestart(1)) THEN

   CALL log_timer_in()
   flag = .TRUE.

!  ---file attributes
   CALL set_modfiles_atts(io_fincon,ics_sed,2)
   filepars = modfiles(io_fincon,ics_sed,2)
   numvars = filepars%novars + 1

!  ---open file
   IF (master) CALL open_filepars(filepars)

!  ---variable attributes
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
   CALL set_modvars_atts(io_fincon,ics_sed,2,varatts,numvars)

!  ---write attributes
   IF (master) CALL write_atts_mod(filepars,varatts,numvars)

   itrestart = 1

ENDIF

!
!2. Write sediment data
!----------------------
!

IF (itrestart.GT.0.AND.itrestart.LE.norestarts) THEN

   IF (nt.EQ.ntrestart(itrestart)) THEN
      IF (itrestart.GT.1) THEN
         CALL log_timer_in()
         flag = .TRUE.
      ENDIF

!
!2.1 Date/time
!-------------
!

      CALL write_time(CDateTime,filepars)
      varid = 2

!
!2.2  Sediment concentrations
!----------------------------
!

      IF (iopt_sed.EQ.1) THEN
         IF (iopt_sed_mode.EQ.2.OR.iopt_sed_mode.EQ.3) THEN
            lbounds = (/1-nhalo,1-nhalo,1,1/)
            ubounds = (/nc+nhalo,nr+nhalo,nz,nf/)
            vecids = (/(ivar,ivar=varid,varid+nf-1)/)
            CALL combine_write_mod(cvol,filepars,0,lbounds,ubounds,'C  ',&
                                 & varatts=varatts,vecids=vecids)
            varid = varid + nf
         ENDIF
      ELSEIF (iopt_sed.EQ.2) THEN
         lbounds = (/1-nhalo,1-nhalo,1,1/)
         ubounds = (/nc+nhalo,nr+nhalo,nz,nf/)
         vecids = (/(ivar,ivar=varid,varid+nf-1)/)
         CALL combine_write_mod(cnump,filepars,0,lbounds,ubounds,'C  ',&
                              & varatts=varatts,vecids=vecids)
         varid = varid + nf
      ENDIF
   
!
!2.3 Bed fractions
!-----------------
!

      IF (nfp.GT.1) THEN
         lbounds = (/0,0,1,1/); ubounds = (/nc,nr,nb,nfp/)
         vecids = (/(ivar,ivar=varid,varid+nfp-1)/)
         CALL combine_write_mod(bed_fraction,filepars,0,lbounds,ubounds,&
              & 'C  ',varatts=varatts,vecids=vecids)
         varid = varid + nfp
      ENDIF

!
!2.4 Skin roughness
!------------------
!   

      IF (iopt_sed_rough.EQ.2) THEN
         lbounds(1:2) = 0; ubounds(1:2) = (/nc,nr/)
         CALL combine_write_mod(zroughatc_sed,filepars,varid,lbounds(1:2),&
              & ubounds(1:2),'C  ',varatts=varatts)
         varid = varid + 1         
      ENDIF
      
!
!2.5 Open boundary arrays
!------------------------
!

      IF (iopt_obc_sed.EQ.1) THEN
         vecids(1:nf) = (/(ivar,ivar=varid,varid+nf-1)/)
         IF (nobu.GT.0) THEN
            IF (iopt_sed.EQ.1) THEN
               CALL combine_write_stats_glb(obcsedatu,filepars,0,nobu,&
                                          & nobuprocs,indexobuprocs,&
                                          & varatts=varatts,vecids=vecids)
            ELSEIF (iopt_sed.EQ.2) THEN
               CALL combine_write_stats_glb(obccnumpatu,filepars,0,nobu,&
                                          & nobuprocs,indexobuprocs,&
                                          & varatts=varatts,vecids=vecids)
            ENDIF
            varid = varid + nf
         ENDIF
         IF (nobv.GT.0) THEN
            IF (iopt_sed.EQ.1) THEN
               CALL combine_write_stats_glb(obcsedatv,filepars,0,nobv,&
                                          & nobvprocs,indexobvprocs,&
                                          & varatts=varatts,vecids=vecids)
            ELSEIF (iopt_sed.EQ.2) THEN
               CALL combine_write_stats_glb(obccnumpatv,filepars,0,nobv,&
                                          & nobvprocs,indexobvprocs,&
                                          & varatts=varatts,vecids=vecids)
            ENDIF
            varid = varid + nf
         ENDIF
      ENDIF

!
!2.5 Next output index
!---------------------
!

      itrestart = itrestart + 1

   ENDIF

ENDIF

!
!3. Finalise
!-----------
!

IF ((cold_start.AND.ntrestart(1).EQ.0).OR.&
  & (itrestart.GT.norestarts)) THEN
   IF (master) CALL close_filepars(filepars)
   DEALLOCATE (varatts)
ENDIF

modfiles(io_fincon,ics_sed,2) = filepars

IF (flag) CALL log_timer_out() 


RETURN

END SUBROUTINE write_sedics

!========================================================================

SUBROUTINE write_sed_spec
!************************************************************************
!
! *write_sed_spec* Write particle properties in standard format
!
! Author - Boudewijn Decrop and  Alexander Breugem (IMDC),
!          Patrick Luyten (MUMM)
!
! Version - @(COHERENS)Sediment_Finalisation.f90  V2.9
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE sedarrays
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_sed_spec'
CALL log_timer_in()

!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(io_sedspc,1,2)
filepars = modfiles(io_sedspc,1,2)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(io_sedspc,1,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2. Write data
!-------------
!

CALL write_vars(dp,filepars,1,varatts=varatts)
CALL write_vars(rhos,filepars,2,varatts=varatts)
CALL write_vars(bstres_cr_cst,filepars,3,varatts=varatts)
CALL write_vars(ws_cst,filepars,4,varatts=varatts)

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_sedspc,1,2) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE write_sed_spec
