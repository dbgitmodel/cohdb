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
! limitations under the Licence
!******************************************************************************
!
! *Dredging_Equations* 
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.11.2
!
! $Date: 2016-04-04 09:25:08 +0200 (Mon, 04 Apr 2016) $
!
! $Revision: 931 $
!
! Description - Dredging and relocation routines
!
! Routines - calculate_dredging_volume, check_eligibility,
!      create_dredging_track, dar_dredging, dar_equation, dar_init,
!      dar_masks, dar_relocation, dar_sorting_dredge, dar_sorting_relocate,
!      determine_campaign_params, determine_dar_times,
!      determine_plume_location, equal_dredging, get_spatial_counter,
!      priority_dredging, relocate, set_dar_flags, spatial_coupling,
!      track_dredging, Zcoord_dredge, Zcoord_relocate
!
!******************************************************************************
!

SUBROUTINE calculate_dredging_volume(n) 
!******************************************************************************
!
! *calculate_dredging_volume* calculates the volume to be dredged (if needed)
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - dar_equation
!
! External calls - 
!
! Module calls - collect_log_expr, sum2_vars
!
!******************************************************************************
!
USE dararrays
USE depths
USE grid
USE gridpars
USE iopars
USE paral_comms, ONLY: collect_log_expr
USE paral_utilities, ONLY: sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name     Type   Purpose
!------------------------------------------------------------------------------
!*n*       INTEGER Counter indicating the campaign or dredging site number
!
!******************************************************************************
!
!*Local variables
!
LOGICAL :: flagloc, flag
INTEGER, DIMENSION(4) :: nhdims
LOGICAL, DIMENSION(ncloc,nrloc) :: mask_crit
REAL, DIMENSION(ncloc,nrloc) :: array2d


procname(pglev+1) = 'calculate_dredging_volume'
CALL log_timer_in()

!
!1. Check whether a dredging campaign should start
!-------------------------------------------------
!

mask_crit = mask_dredging(:,:,n).AND. &
          & depmeanatc(1:ncloc,1:nrloc).LT.DredgingSites(n)%target_depth

flagloc = ANY(mask_crit)
flag = collect_log_expr(flagloc,1)
IF (.NOT.flag) GOTO 1000

!
!2. Total dredged volume
!-----------------------
!

IF (flag) THEN
   WHERE (mask_crit)
      array2d = (DredgingSites(n)%target_depth+&
               & DredgingSites(n)%over_depth-depmeanatc(1:ncloc,1:nrloc))*garea
   END WHERE
   CALL sum2_vars(array2d,V_dredge_total(n),nhdims,'C  ',0,commall=.TRUE.,&
                & mask=mask_crit)
   campaign(n) = .TRUE.
ENDIF

1000 CALL log_timer_out()


RETURN

END SUBROUTINE calculate_dredging_volume

!==============================================================================

SUBROUTINE check_eligibility(siteno)
!******************************************************************************
!
! *check_eligibility* determines whether a relocation site is eligible
!                     for relocation according to the relocation criterium
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - spatial_coupling
!
! External calls - 
!
! Module calls - collect_vars, error_alloc
!
!******************************************************************************
!
USE dararrays
USE darswitches
USE depths
USE gridpars
USE iopars
USE paral_comms, ONLY: collect_log_expr
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

INTEGER, INTENT(IN) :: siteno

! Name      Type    Purpose
!------------------------------------------------------------------------------
!*siteno*   INTEGER Number of the relocation site
!
!*******************************************************************************
!
!*Local variables
!
LOGICAL :: flag, flagloc


procname(pglev+1) = 'check_eligibility'
CALL log_timer_in()


!---check minimum depth criterion
IF (iopt_dar_relocation_criterion.EQ.1) THEN
   flagloc = ANY(mask_relocation(1:ncloc,1:nrloc,siteno).AND.&
             & depmeanatc(1:ncloc,1:nrloc).LT.RelocationSites(siteno)%min_depth)
   flag = collect_log_expr(flagloc,1)
   IF (flag) RelocationSites(siteno)%eligibility = .FALSE.

!---check maximum volume criterion
ELSEIF (V_relocate_total(siteno).GT.RelocationSites(siteno)%max_volume) THEN
   RelocationSites(siteno)%eligibility = .FALSE.
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE check_eligibility

!==============================================================================

SUBROUTINE create_dredging_track(n)
!******************************************************************************
!
! *create_dredging_track* define the grid cells along the dredging track and
!                         the time the ship stays in each cell
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - dar_init
!
! External calls - 
!
! Module calls - data_to_model_hcoords, distance_pts, error_alloc,
!                hrelativecoords_init, local_proc
!
!******************************************************************************
!
USE dararrays
USE darpars
USE darswitches
USE datatypes
USE grid
USE gridpars
USE iopars
USE timepars
USE syspars
USE datatypes_init, ONLY: hrelativecoords_init
USE error_routines, ONLY: error_alloc, error_alloc_struc
USE grid_interp, ONLY: data_to_model_hcoords
USE grid_routines, ONLY: distance_pts, local_proc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER Number of campaign or dredging site
!
!******************************************************************************
!
!*Local variables
!
INTEGER :: i, j, l, ll, nopts, noptsnew, ntstart, ntx, shipno, siteno
INTEGER, DIMENSION(1) :: iloc, ilocnew, jloc, jlocnew
REAL, DIMENSION(1) :: x, y
REAL :: dist, track_delx, track_dely, track_dist, Uship 
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask_track
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ipos, iposnew, jpos, jposnew
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: xpos, xposnew, ypos, yposnew
TYPE (HRelativeCoords), SAVE, ALLOCATABLE, DIMENSION(:,:) :: trackcoords



procname(pglev+1) = 'create_dredging_track'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ntstart = OperationTimes(n)%dredging_times(1)

ALLOCATE(mask_track(ncloc,nrloc),STAT=errstat)
CALL error_alloc('mask_track',2,(/ncloc,nrloc/),kndlog)
mask_track = .FALSE.

ALLOCATE(ipos(maxtrackpoints),STAT=errstat)
CALL error_alloc('ipos',1,(/maxtrackpoints/),kndint)
ipos = 0

ALLOCATE(iposnew(maxtrackpoints),STAT=errstat)
CALL error_alloc('iposnew',1,(/maxtrackpoints/),kndint)
iposnew = 0

ALLOCATE(jpos(maxtrackpoints),STAT=errstat)
CALL error_alloc('jpos',1,(/maxtrackpoints/),kndint)
jpos = 0

ALLOCATE(jposnew(maxtrackpoints),STAT=errstat)
CALL error_alloc('jposnew',1,(/maxtrackpoints/),kndint)
jposnew = 0

ALLOCATE(xpos(maxtrackpoints),STAT=errstat)
CALL error_alloc('xpos',1,(/maxtrackpoints/),kndrtype)
xpos = 0.0

ALLOCATE(ypos(maxtrackpoints),STAT=errstat)
CALL error_alloc('ypos',1,(/maxtrackpoints/),kndrtype)
ypos = 0.0

ALLOCATE(xposnew(maxtrackpoints),STAT=errstat)
CALL error_alloc('xposnew',1,(/maxtrackpoints/),kndrtype)
xpos = 0.0

ALLOCATE(yposnew(maxtrackpoints),STAT=errstat)
CALL error_alloc('yposnew',1,(/maxtrackpoints/),kndrtype)
yposnew = 0.0

ALLOCATE (trackcoords(maxtrackpoints,1),STAT=errstat)
CALL error_alloc_struc('trackcoords',2,(/maxtrackpoints,1/),&
                     & 'HRelativeCoords')

!
!2. Initialise parameters
!------------------------
!
!---relative coordinates
CALL hrelativecoords_init(trackcoords,.TRUE.)

!---ship and site number
shipno = MERGE(Campaigns(n)%nr_ship,DredgingSites(n)%nr_ship,&
             & iopt_dar_dredging_criterion.EQ.3)
siteno = MERGE(Campaigns(n)%nr_dredgingsite,n,&
             & iopt_dar_dredging_criterion.EQ.3)

!
!3. Create the dredging track
!----------------------------
!
!3.1 Cell locations
!------------------
!

nopts = notrackpoints(siteno)

CALL data_to_model_hcoords(trackcoords,'UV ',nopts,1,xtrack(1:nopts,n),&
                         & ytrack(1:nopts,n))

!
!3.2 Eliminate dry locations
!---------------------------
!

noptsnew = 0
l_320: DO l=1,nopts
   i = trackcoords(l,1)%icoord; j = trackcoords(l,1)%jcoord
   IF (i.NE.int_fill) THEN
      noptsnew = noptsnew + 1
      iposnew(noptsnew) = i; jposnew(noptsnew) = j
      xposnew(noptsnew) = xtrack(l,n); yposnew(noptsnew) = ytrack(l,n)
   ENDIF
ENDDO l_320

!
!3.3 Account for possible double dredging at turning points of the track
!-----------------------------------------------------------------------
!

itrack(1,n) = iposnew(1); jtrack(1,n) = jposnew(1)
xpos(1) = xposnew(1); yposnew(1) = ypos(1)

ll = 1
l_330: DO l=2,nopts-1
   ll = ll + 1
   itrack(ll,n) = iposnew(l); jtrack(ll,n) = jposnew(l)
   xpos(l) = xposnew(ll); ypos(l) = yposnew(ll)
   IF (iposnew(l+1).EQ.iposnew(l-1).AND.jposnew(l+1).EQ.jposnew(l-1)) THEN
      ll = ll + 1
      itrack(ll,n) = iposnew(l); jtrack(ll,n) = jposnew(l)
      xpos(ll) = xposnew(l); ypos(ll) = yposnew(l)
   ENDIF
ENDDO l_330
ll = ll + 1
itrack(ll,n) = iposnew(nopts); jtrack(ll,n) = jposnew(nopts)
xpos(ll) = xposnew(nopts); yposnew(ll) = ypos(nopts)

!
!3.4 Cell locations and time spent in cells along dredging track
!---------------------------------------------------------------
!
!---initialise track counters
track_points_counter(n) = 1; nt_track(1,n) = 1
!---ship's velocity
Uship = Ships(shipno)%U_ship
!---initial position of the ship
l = 1
iloc(1) = ipos(1); jloc(1) = jpos(1)
x(1) = xpos(1); y(1) = ypos(1)
!---initialise distance with next tracking position and projection
!   onto X- and Y-axes
track_dist = distance_pts(xpos(1),xpos(2),ypos(1),ypos(2),2,3)
track_delx = SIGN(1.0,xpos(2)-xpos(1))*&
           & distance_pts(xpos(1),xpos(2),ypos(1),ypos(1),2,3)
track_dely = SIGN(1.0,ypos(2)-ypos(1))*&
           & distance_pts(xpos(1),xpos(1),ypos(1),ypos(2),2,3)

ntx_340: DO ntx=ntstart,nstep/ic3d
!  ---ship's position at next time
   dist = Uship*delt3d
!  ---ship reached next tracking position
   DO WHILE (dist.GT.track_dist)
      l = l + 1
!     --ship reached last tracking position
      IF (l.GT.nopts) THEN
         EXIT ntx_340
!     --update distance to next track location
      ELSE
         iloc(1) = ipos(l); jloc(1) = jpos(l)
         x(1) = xpos(l); y(1) = ypos(l)
         dist = dist - track_dist
         track_dist = distance_pts(xpos(l),xpos(l+1),ypos(l),ypos(l+1),2,3)
         track_delx = SIGN(1.0,xpos(l+1)-xpos(l))*&
                    & distance_pts(xpos(l),xpos(l+1),ypos(l),ypos(l),2,3)
         track_dely = SIGN(1.0,ypos(l+1)-ypos(l))*&
                    & distance_pts(xpos(l),xpos(l),ypos(l),ypos(l+1),2,3)
      ENDIF
   ENDDO
!  ---ship has not reached next tracking position
   IF (dist.LT.track_dist) THEN
!     --new ship's position
      x(1) = x(1) + dist*track_delx/track_dist
      y(1) = y(1) + dist*track_dely/track_dist
      CALL data_to_model_hcoords(trackcoords(1:1,:),'UV ',1,1,x,y)
      ilocnew(1) = trackcoords(1,1)%icoord
      jlocnew(1) = trackcoords(1,1)%jcoord
!     --ship remains in the same cell
      IF (ilocnew(1).EQ.iloc(1).AND.jlocnew(1).EQ.jloc(1)) THEN
         nt_track(track_points_counter(n),n) = &
       & nt_track(track_points_counter(n),n) + 1
!     --ship sails to new cell
      ELSE
         track_points_counter(n) = track_points_counter(n) + 1
         nt_track(track_points_counter(n),n) = 1
         itrack(track_points_counter(n),n) = ilocnew(1)
         jtrack(track_points_counter(n),n) = jlocnew(1)
      ENDIF
   ENDIF
ENDDO ntx_340

!
!3.5 Dredging mask
!-----------------
!

l_350: DO l=1,track_points_counter(n)
   i = itrack(l,n); j = jtrack(l,n)
   iloc(1) = i - nc1loc + 1; jloc(1) = j - nr1loc + 1
   IF (local_proc(iloc(1),jloc(1))) THEN
      mask_track(iloc(1),jloc(1)) = maskatc_int(iloc(1),jloc(1))
   ENDIF
ENDDO l_350
mask_dredging(:,:,siteno) = mask_dredging(:,:,siteno).AND.mask_track

!
!4. Deallocate arrays
!--------------------
!

DEALLOCATE (ipos,iposnew,jpos,jposnew,mask_track,xpos,xposnew,ypos,yposnew)
DEALLOCATE (trackcoords)

CALL log_timer_out()


RETURN

END SUBROUTINE create_dredging_track

!==============================================================================

SUBROUTINE dar_dredging(n)
!******************************************************************************
!
! *dar_dredging* selcts type of method for dredging operation
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - dar_equation
!
! External calls - equal_dredging, priority_dredging, track_dredging 
!
! Module calls - 
!
!*****************************************************************************
!
USE darswitches
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER Number of campaign or dredging site
!
!******************************************************************************
!


procname(pglev+1) = 'dar_dredging'
CALL log_timer_in()

SELECT CASE (iopt_dar_time_dredge)
   CASE (1); CALL equal_dredging(n)
   CASE (2); CALL priority_dredging(n)
   CASE (3); CALL track_dredging(n)
END SELECT

CALL log_timer_out()


RETURN

END SUBROUTINE dar_dredging

!==============================================================================

SUBROUTINE dar_equation
!******************************************************************************
!
! *dar_equation* Mian program unit for dredging and relocation
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - calculate_dredging_volume, dar_relocation, depth_at_nodes,
!                  determine_campaign_params, determine_dar_times,
!                  set_dar_flags, water_depths
!
! Module calls - exchange_mod
!
!******************************************************************************
!
USE dararrays
USE darpars
USE darswitches
USE depths
USE gridpars
USE iopars
USE modids
USE paralpars
USE switches
USE timepars
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: cval
INTEGER :: n, nolocs
INTEGER, DIMENSION(4) :: nhdims


procname(pglev+1) = 'dar_equation'
CALL log_timer_in()

nolocs = MERGE(nocampaigns,nodredgingsites,iopt_dar_dredging_criterion.EQ.3)

n_100: DO n=1,nolocs
   
 !*************************
 ! Determine if dredging or relocation takes place at the current timestep nt
 !   a) criterium-based dredging
 !     --> If current timestep corrsponds to a monitoring step: check the
 !         dredging criterion
 !     --> When the dredging criterion is met: start a dredging campaign
 !     --> During a dredging campaign: check the dredging/relocation flags
 !         each time step
 !   b) user-defined dredging:
 !     --> If current timestep corrsponds to the start of a dredging campaign
 !         (calculated at initialisation). start a dredging campaign
 !     --> During a dredging campaign: check the dredging/relocation flags each
 !         time step
   !**************************
   
   IF (.NOT.campaign(n)) THEN

!     ---criterion-based dredging and relocation
      IF (iopt_dar_dredging_criterion.NE.3) THEN
         IF (MOD(nt,DredgingSites(n)%monitoring_frequency).EQ.0) THEN
            CALL calculate_dredging_volume(n)
            IF (campaign(n)) THEN
               CALL determine_campaign_params(n)
               CALL determine_dar_times(n)
            ENDIF
         ENDIF

 !    ---user-defined dredging/relocation
      ELSEIF (iopt_dar_dredging_criterion.EQ.3.AND.nt.EQ.campaign_start(n).AND.&
           & ANY(mask_dredging(:,:,Campaigns(n)%nr_dredgingsite))) THEN
         campaign(n) = .TRUE.
      ENDIF

   ENDIF
   
!   ---set operational flags for dredging and/or relocation
   IF (campaign(n)) CALL set_dar_flags(n)
   
 !*******************************
 !.Carry out the actual dredging and relocation activities
 ! a) if dredging/relocation flag is set to true -> dredging/relocation is
 !    carried out this timestep
 ! b) if dredging/relocation is instantaneous, dredging/relocation only take
 !     one timestep.
 ! The flags are set to .FALSE. here instead of in set_dar_flags. After
 ! relocation, the campaign is finished and set to .FALSE.
 !*******************************

!  ---perform dredging   
   IF (dredging(n)) THEN
      CALL dar_dredging(n)
      IF (iopt_dar_duration.EQ.1) THEN
         dredging(n) = .FALSE.
         dredge_counter(n) = 1
      ENDIF
   ENDIF

!  ---perform relocation   
   IF (relocation(n)) THEN
      CALL dar_relocation(n)
      IF (iopt_dar_duration.EQ.1) THEN
         relocation(n) = .FALSE.
         campaign(n) = .FALSE.
      ENDIF
   ENDIF

!  --exchange halo sections
   IF (iopt_MPI.EQ.1) THEN
      nhdims = nhalo
      CALL exchange_mod(depmeanatc,(/1-nhalo,1-nhalo/),nhdims,iarr_depmeanatc,&
                      & corners=.FALSE.)
   ENDIF

!  ---update bathymetry at nodes and total water depths 
   CALL depth_at_nodes
   CALL water_depths

ENDDO n_100

!---monitoring
IF (nt.EQ.nstep) THEN
   IF (master.AND.loglev1.GT.0) THEN
      WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//&
                         &'Superfluous volumes of relocated sediment:'
      n_200: DO n=1,norelocationsites
         WRITE (cval,'(I12)') n
         WRITE (iolog,'(A,G12.7,A)') REPEAT(' ',pglev-1)//&
                                 & 'Relocation site' //TRIM(cval)//': ', &
                                 & V_relocate_extra(n),' m^3'
      ENDDO n_200
   ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE dar_equation

!==============================================================================

SUBROUTINE dar_init
!******************************************************************************
!
! *dar_init* initialise dredging and relocation parameters and arrays
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!  
! Description -
!
! Reference -
!
! Calling program - dar_equation
!
! External calls - create_dredging_track, dar_masks, determine_campaign_params,
!                  determine_dar_times
!
! Module calls - num_time_steps
!
!******************************************************************************
!
USE dararrays
USE darpars
USE darswitches
USE iopars
USE sedpars
USE syspars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out, &
                       & num_time_steps

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=lentime) :: CoupleTime
INTEGER :: m, n, nolocs, nstep_campaign, nstep_coupling, shipno, siteno


procname(pglev+1) = 'dar_init'
CALL log_timer_in()

nolocs = MERGE(nocampaigns,nodredgingsites,iopt_dar_dredging_criterion.EQ.3)

!
!1. Define grid masks for dredging and relocation (sub)zones
!-----------------------------------------------------------
!

CALL dar_masks

!
!2. Dredging/relocation parameters
!---------------------------------
!

n_210: DO n=1,nolocs

!  ---dredging/relocation activators   
   dredging(n) = .FALSE.; relocation(n) = .FALSE.
   campaign(n) = .FALSE.; eligibility(n) = .TRUE.

!  ---ship and site number
   shipno = MERGE(Campaigns(n)%nr_ship,DredgingSites(n)%nr_ship,&
                & iopt_dar_dredging_criterion.EQ.3)
   siteno = MERGE(Campaigns(n)%nr_dredgingsite,n,&
                & iopt_dar_dredging_criterion.EQ.3)

!  ---time spread dredging
   IF (iopt_dar_duration.EQ.2)THEN
      nt_spread(n) = MERGE(Campaigns(n)%t_spread,DredgingSites(n)%t_spread,&
                         & iopt_dar_dredging_criterion.EQ.3)
   ENDIF
   
!  ---time lagged dredging
   IF (iopt_dar_coupling_time.EQ.2) THEN
      nt_lag(n) = MERGE(Campaigns(n)%t_lag,DredgingSites(n)%t_lag,&
                      & iopt_dar_dredging_criterion.EQ.3)

!  ---full dredging cycle
!  ---PLU (times ic3d to get 2-D time index)
   ELSEIF (iopt_dar_coupling_time.EQ.3) THEN
      nt_wait(n) = ic3d*CEILING(Ships(shipno)%t_wait/delt3d)
      nt_fill(n) = ic3d*CEILING(Ships(shipno)%t_fill/delt3d)
      nt_sail_dredge(n) = ic3d*CEILING(Ships(shipno)%t_sail_dredge/delt3d)
      nt_sail_relocate(n) = ic3d*CEILING(Ships(shipno)%t_sail_relocate/delt3d)
      nt_idle(n) = ic3d*CEILING(Ships(shipno)%t_idle/delt3d)
      nt_relocate(n) = ic3d*CEILING(Ships(shipno)%t_relocate/delt3d)
      nt_cycle(n) = nt_fill(n) + nt_sail_dredge(n) + nt_relocate(n) + &
                  & nt_sail_relocate(n) + nt_idle(n)
   ENDIF

!  ---plume parameters
   IF (iopt_dar_effect.EQ.2.OR.iopt_dar_effect.EQ.3) THEN
      p_suctionhead(n,:) = Ships(shipno)%p_suctionhead(1:nf)
      p_overflow(n,:) = Ships(shipno)%p_overflow(1:nf)
      p_relocation(n,:) = Ships(shipno)%p_relocation(1:nf)
      p_passive(n,:) = Ships(shipno)%p_passive(1:nf)
      D_dredge_max(n) = Ships(shipno)%D_dredge_max
      D_dredge_min(n) = Ships(shipno)%D_dredge_min
      D_relocate_max(n) = Ships(shipno)%D_relocate_max
      D_relocate_min(n) = Ships(shipno)%D_relocate_min
      H_dredge_max(n) = Ships(shipno)%H_dredge_max
      H_dredge_min(n) = Ships(shipno)%H_dredge_min
      H_relocate_max(n) = Ships(shipno)%H_relocate_max
      H_relocate_min(n) = Ships(shipno)%H_relocate_min
   ENDIF

!  ---user-defined dredging/relocation
   IF (iopt_dar_dredging_criterion.EQ.3) THEN
!     --start date of the campaign
      CALL num_time_steps(CStartDateTime,Campaigns(n)%CampaignStartDateTime,&
                        & delt3d,nstep_campaign)
      campaign_start(n) = nstep_campaign
!     ---relevant parameters for the campaign
      CALL determine_campaign_params(n)
!     --dredging and relocation times for the campaign
      CALL determine_dar_times(n)
   ENDIF

!  ---track dredging
   IF (iopt_dar_time_dredge.EQ.3) THEN
      CALL create_dredging_track(n,siteno,shipno)
      V_ship_rest(n)= Ships(shipno)%V_ship
   ENDIF

!  ---coupling times
   IF (iopt_dar_coupling_space.EQ.3) THEN
      m_211: DO m=1,nocouplingtimes(n)
         IF (iopt_dar_dredging_criterion.EQ.3) THEN
            CoupleTime = Campaigns(n)%CoupleTime(m)
         ELSE
            CoupleTime = DredgingSites(n)%CoupleTime(m)
         ENDIF
         CALL num_time_steps(CStartDateTime,CoupleTime,delt3d,nstep_coupling)
         OperationTimes(n)%coupling_times(m) = nstep_coupling
      ENDDO m_211
   ENDIF

!  ---initialise counters
   couple_counter(n) = 1; dredge_counter(n) = 1; relocate_counter(n) = 1
   cycle_counter(n) = 1; priority_counter(n) = 1; track_cell_counter(n) = 1

ENDDO n_210

CALL log_timer_out()

RETURN

END SUBROUTINE dar_init

!==============================================================================

SUBROUTINE dar_masks
!******************************************************************************
!
! *dar_masks* Mask arrays for dredging/relocation (sub)zones
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - dar_initialisation
!
! Module calls - polygon_mask
!
!******************************************************************************
!
USE dararrays
USE darpars
USE darswitches
USE grid
USE gridpars
USE iopars
USE math_library, ONLY: polygon_mask
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables

INTEGER :: m, n 

procname(pglev+1) = 'dar_masks'
CALL log_timer_in()

!
!1. Dredging sites
!-----------------
!

n_110: DO n=1,nodredgingsites
   CALL polygon_mask(&
       & DredgingSites(n)%Xpol(1:DredgingSites(n)%nr_coordpol),&
       & DredgingSites(n)%Ypol(1:DredgingSites(n)%nr_coordpol),&
       & DredgingSites(n)%nr_coordpol,&
       & gxcoord(1:ncloc,1:nrloc),gycoord(1:ncloc,1:nrloc),&
       & ncloc,nrloc,mask_dredging(:,:,n),.TRUE.)
   mask_dredging(:,:,n) = maskatc_int.AND.mask_dredging(:,:,n)
ENDDO n_110

!
!2. Relocation sites
!-------------------
!

n_210: DO n=1,norelocationsites
   CALL polygon_mask(&
       & RelocationSites(n)%Xpol(1:RelocationSites(n)%nr_coordpol),&
       & RelocationSites(n)%Ypol(1:RelocationSites(n)%nr_coordpol),&
       & RelocationSites(n)%nr_coordpol,&
       & gxcoord(1:ncloc,1:nrloc),gycoord(1:ncloc,1:nrloc),&
       & ncloc,nrloc,mask_relocation(:,:,n),.TRUE.)
   mask_relocation(:,:,n) = maskatc_int.AND.mask_relocation(:,:,n)
ENDDO n_210

!
!3. Dredging subzones
!--------------------
!

IF (iopt_dar_time_dredge.EQ.2) THEN
   
   n_310: DO n=1,nodredgingsites
   m_310: DO m=1,SubzonesDredging(n)%nr_subzones
      CALL polygon_mask(&
          & SubzonesDredging(n)%Xpol(1:SubzonesDredging(n)%nr_coordpol(m),m),&
          & SubzonesDredging(n)%Ypol(1:SubzonesDredging(n)%nr_coordpol(m),m),&
          & SubzonesDredging(n)%nr_coordpol(m),&
          & gxcoord(1:ncloc,1:nrloc),gycoord(1:ncloc,1:nrloc),&
          & ncloc,nrloc,mask_subzones_dredging(:,:,n,m),.TRUE.)
      mask_subzones_dredging(:,:,n,m) = &
          & mask_dredging(:,:,n).AND.mask_subzones_dredging(:,:,n,m)
   ENDDO m_310
   ENDDO n_310

ENDIF

!
!4. Relocation subzones
!----------------------
!

IF (iopt_dar_time_relocate.EQ.3) THEN

   n_410: DO n=1,norelocationsites
   m_410: DO m=1,SubzonesRelocation(n)%nr_subzones
      CALL polygon_mask(&
        & SubzonesRelocation(n)%Xpol(1:SubzonesRelocation(n)%nr_coordpol(m),m),&
        & SubzonesRelocation(n)%Ypol(1:SubzonesRelocation(m)%nr_coordpol(m),m),&
        & SubzonesRelocation(n)%nr_coordpol(m),&
        & gxcoord(1:ncloc,1:nrloc),gycoord(1:ncloc,1:nrloc),&
        & ncloc,nrloc,mask_subzones_relocation(:,:,n,m),.TRUE.)
      mask_subzones_relocation(:,:,n,m) = &
        & mask_relocation(:,:,n).AND.mask_subzones_relocation(:,:,n,m)
   ENDDO m_410
   ENDDO n_410

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE dar_masks

!==============================================================================

SUBROUTINE dar_relocation(n)
!******************************************************************************
!
! *dar_relocation* relocate sediment
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - dar_equation
!
! External calls - 
!
! Module calls - copy_vars, error_alloc, maxloc_vars, polygon_mask, rng_close,
!                rng_init, rng_open, rng_standard_uniform
!
!******************************************************************************
!
USE dararrays
USE darpars
USE darswitches
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE math_library, ONLY: polygon_mask
USE paral_comms, ONLY: copy_vars
USE paral_utilities, ONLY: maxloc_vars
USE rng_library, ONLY: rng_close, rng_init, rng_open, rng_standard_uniform
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER Number of campaign or dredging site
!
!******************************************************************************
!
!*Local variables
!
INTEGER :: ic, iglb, iran, jglb, l, ngen, shipno, siteno
INTEGER, DIMENSION(2) :: maxpos
REAL :: depmax, originx, originy, xran
REAL, DIMENSION(nocirclepoints) :: theta
REAL, DIMENSION(nocirclepoints+1) :: X_circle, Y_circle
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: mask_circle


procname(pglev+1) = 'dar_relocation'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

IF (iopt_dar_time_relocate.EQ.2) THEN
   ALLOCATE(mask_circle(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('mask_circle',2,(/ncloc,nrloc/),kndlog)
ENDIF

!
!2. Initialise parameters
!------------------------
!

IF (iopt_dar_dredging_criterion.EQ.3) THEN
   shipno = Campaigns(n)%nr_ship
   siteno = Campaigns(n)%relocationsites(spatial_counter(n))
ELSE
   shipno = DredgingSites(n)%nr_ship
   siteno = DredgingSites(n)%relocationsites(spatial_counter(n))
 ENDIF

!
!3. Define relocation area
!-------------------------
!
!3.1 Equal spreading
!-------------------
!

IF (iopt_dar_time_relocate.EQ.1) THEN
   CALL relocate(siteno,mask_relocation(:,:,siteno))

!
!3.2 Deepest point
!-----------------
!

ELSEIF (iopt_dar_time_relocate.EQ.2) THEN

!  ---location of deepest point
   CALL maxloc_vars(depmeanatc(1:ncloc,1:nrloc),depmax,maxpos,iarr_depmeanatc,&
                  & commall=.TRUE.,mask=maskatc_int)
   iglb = maxpos(1); jglb = maxpos(2)

!  ---define points on the circle wirth given radius
   theta = (/(l,l=1,nocirclepoints)/)
   theta = theta*twopi/nocirclepoints
   originx = 0.25*(gxcoordglb(iglb,jglb)+gxcoordglb(iglb+1,jglb)+&
                 & gxcoordglb(iglb,jglb+1)+gxcoordglb(iglb+1,jglb+1))
   originy = 0.25*(gycoordglb(iglb,jglb)+gycoordglb(iglb+1,jglb)+&
                 & gycoordglb(iglb,jglb+1)+gycoordglb(iglb+1,jglb+1))
   ic_320: DO ic = 1,nocirclepoints
      X_circle(ic) = originx + Ships(shipno)%R_circle*COS(theta(ic))
      Y_circle(ic) = originy + Ships(shipno)%R_circle*SIN(theta(ic))
   ENDDO ic_320
   X_circle(nocirclepoints+1) = X_circle(1)
   Y_circle(nocirclepoints+1) = Y_circle(1)

!  ---define mask for points inside the circle
   CALL polygon_mask(X_circle,Y_circle,nocirclepoints+1,&
                   & gxcoord(1:ncloc,1:nrloc),gycoord(1:ncloc,1:nrloc),&
                   & ncloc,nrloc,mask_circle,info=.TRUE.)
   mask_circle = mask_circle.AND.mask_relocation(1:ncloc,1:nrloc,siteno)

!  ---relocate   
   CALL relocate(siteno,mask_circle)

! 
!3.3 Random order (subzones)
!---------------------------
!
  
ELSEIF (iopt_dar_time_relocate.EQ.3) THEN
   IF (master) THEN
      CALL rng_init
      CALL rng_open(ngen)
      CALL rng_standard_uniform(xran,ngen)
      iran = FLOOR(SubzonesRelocation(siteno)%nr_subzones*xran) + 1
      CALL rng_close(ngen)
   ENDIF
   IF (iopt_MPI.EQ.1) CALL copy_vars(iran,0)
   CALL relocate(siteno,mask_subzones_relocation(:,:,siteno,iran))
ENDIF

!
!4. Deallocate arrays
!--------------------
!

IF (iopt_dar_time_relocate.EQ.2) THEN
   DEALLOCATE (mask_circle)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE dar_relocation

!==============================================================================

SUBROUTINE dar_sorting_dredge(mask,erosion)
!******************************************************************************
!
! *dar_sorting_dredge* adjust layer thicknesses and sediment fractions 
!                      in the bed due to dredging
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - dar_equation
!
! External calls - 
!
! Module calls - error_alloc
!
!******************************************************************************
!
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE sedarrays
USE sedpars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
LOGICAL, DIMENSION(1:ncloc,1:nrloc), INTENT(IN) :: mask
REAL, DIMENSION(1:ncloc,1:nrloc), INTENT(IN) :: erosion

! Name      Type    Purpose
!-------------------------------------------------------------------------------
!*mask*     LOGICAL Mask array for dredging
!*erosion*  REAL    Dredging depth                                           [m]
!
!*******************************************************************************
!
!*Local variables
!
INTEGER :: i, j, k, removed_layers
REAL :: bed_layer_depth, bed_layer_thickness_rest, thickness_rest
REAL, DIMENSION(nb) :: bed_layer_thick_old
REAL, DIMENSION(nb,nf) :: bed_frac_old


procname(pglev+1) = 'dar_sorting_dredge'
CALL log_timer_in()

j_110: DO j=1,nrloc
i_110: DO i=1,ncloc

   IF (mask(i,j)) THEN

!
!1. Initialise
!-------------
!

      bed_layer_thick_old = bed_layer_thickness(i,j,:)
      bed_frac_old = bed_fraction(i,j,:,:)
      removed_layers = 0

!
!2. Number of removed layers
!---------------------------
!
      bed_layer_thickness_rest = bed_layer_thickness(i,j,1) - erosion(i,j)
      thickness_rest = ABS(bed_layer_thickness_rest)
      IF (bed_layer_thickness_rest.LT.0.0) THEN
         removed_layers = 1
         bed_layer_depth = bed_layer_thick_old(1)
      ENDIF

      k_210: DO k=2,nb
         bed_layer_depth = bed_layer_depth + bed_layer_thick_old(k)
         IF (bed_layer_thickness_rest.LT.0.AND.&
           & thickness_rest.GE.(bed_layer_depth-bed_layer_thickness(i,j,1))) &
       & THEN
            removed_layers = removed_layers + 1
         ENDIF
      ENDDO k_210

!
!3. Bed layer thickness
!----------------------
!

      bed_layer_thickness(i,j,1) = bed_layer_thick_old(1+removed_layers) - &
                                 & thickness_rest
      
      k_310: DO k= 2,nb-removed_layers
         bed_layer_thickness(i,j,k) = bed_layer_thick_old(k+removed_layers)
      ENDDO k_310
      
      k_320: DO k=nb-removed_layers+1,nb
         bed_layer_thickness(i,j,k) = 0.0
      ENDDO k_320

!
!4. Bed fractions
!----------------
!
!     PLU include upper layer !
      k_410: DO k=1,nb-removed_layers
         bed_fraction(i,j,k,:) = bed_frac_old(k+removed_layers,:)
      ENDDO k_410
         
      k_420: DO k=nb-removed_layers+1,nb
         bed_fraction(i,j,k,:) = 0
      ENDDO k_420

   ENDIF

ENDDO i_110
ENDDO j_110

CALL log_timer_out()


RETURN

END SUBROUTINE dar_sorting_dredge

!===============================================================================

SUBROUTINE dar_sorting_relocate(mask,fractions,update)
!*******************************************************************************
!
! *dar_sorting_relocate* adjusts the fractions and layer thicknesses
!                                     in the bed due to relocation
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - dar_equation
!
! External calls - 
!
! Module calls -
!
!******************************************************************************
!
USE gridpars
USE iopars
USE morpharrays
USE morphpars
USE sedarrays
USE sedpars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
LOGICAL, DIMENSION(ncloc,nrloc), INTENT(IN) :: mask
REAL, INTENT(IN) :: update
REAL, DIMENSION(nf), INTENT(IN) :: fractions

! Name       Type     Purpose
!------------------------------------------------------------------------------
!*mask*      LOGICAL .TRUE. for cells satisfying the relocation criterion
!*fractions* REAL    Fraction distribution of the sediment that is deposited
!*update*    REAL    Bed update due to deposition of relocated sediment 
!******************************************************************************
!
!*Local variables
!
INTEGER :: f, i, j, k
LOGICAL :: similarity
REAL :: bedupdate
REAL, DIMENSION(nb) :: bed_layer_thick_old
REAL, DIMENSION(nb,nf) :: bed_frac_old


procname(pglev+1) = 'dar_sorting_relocate'
CALL log_timer_in()

bedupdate = SUM(fractions)*update

j_000: DO j=1,nrloc
i_000: DO i=1,ncloc
   IF (mask(i,j)) THEN
      bed_layer_thick_old = bed_layer_thickness(i,j,:)
      bed_frac_old = bed_fraction(i,j,:,:)

!
!1. Check similarity
!--------------------
!

      IF (bed_layer_thickness(i,j,1).GT.0.0) THEN
         similarity  = .TRUE.
         f_110: DO f=1,nf
            IF (fractions(f).LT.(bed_fraction(i,j,1,f) - similarity_range).OR.&
              & fractions(f).GT.(bed_fraction(i,j,1,f) + similarity_range)) THEN
               similarity = .FALSE.
            ENDIF
         ENDDO f_110
      ELSE
         similarity = .FALSE.
      ENDIF

!
!2. First layer
!--------------
!

      IF (similarity) THEN
         bed_layer_thickness(i,j,1) = bed_layer_thick_old(1) + bedupdate
         bed_fraction(i,j,1,:) = &
              & (bed_frac_old(1,:)*bed_layer_thick_old(1)+fractions(f)*update)/&
              & (bed_layer_thick_old(1)+bedupdate)
      ELSE
         bed_layer_thickness(i,j,1) = bedupdate
         bed_fraction(i,j,1,:) = fractions(f)
      ENDIF

!
!3. Shifting layers 2 to nb-1
!----------------------------
!

      IF (.NOT.similarity) THEN
         k_310: DO k=2,nb-1
            bed_layer_thickness(i,j,k) = bed_layer_thick_old(k-1)
            bed_fraction(i,j,k,:) = bed_frac_old(k-1,:)
         ENDDO k_310

!
!4. Merging of bottom layers
!---------------------------
!

         IF (bed_layer_thick_old(nb).GT.0.0.AND.&
           & bed_layer_thick_old(nb-1).GT.0.0) THEN
            bed_layer_thickness(i,j,nb) = bed_layer_thick_old(nb) + &
                                        & bed_layer_thick_old(nb-1)
            bed_fraction(i,j,nb,:) = &
                           & (bed_frac_old(nb,:)*bed_layer_thick_old(nb)+&
                           & bed_frac_old(nb-1,:)*bed_layer_thick_old(nb-1))/&
                           & (bed_layer_thick_old(nb)+bed_layer_thick_old(nb-1))
         ELSEIF (bed_layer_thick_old(nb).EQ.0.0.AND.&
                 bed_layer_thick_old(nb-1).GT.0.0) THEN
            bed_layer_thickness(i,j,nb) = bed_layer_thick_old(nb-1) 
            bed_fraction(i,j,nb,:) = bed_frac_old(nb-1,:)
         ELSE
            bed_layer_thickness(i,j,nb) = 0.0
            bed_fraction(i,j,nb,:) = 0.0
         ENDIF
      ENDIF

   ENDIF

ENDDO i_000
ENDDO j_000

CALL log_timer_out()


RETURN

END SUBROUTINE dar_sorting_relocate

!==============================================================================

SUBROUTINE determine_campaign_params(n)
!******************************************************************************
!
! *determine_campaign_params initialise campaign parameters
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - dar_equation
!
! External calls - error_alloc
!
! Module calls - sum2_vars
!
!******************************************************************************
!
USE dararrays
USE darpars
USE darswitches
USE depths
USE grid
USE gridpars
USE iopars
USE error_routines, ONLY: error_alloc
USE paral_utilities, ONLY: sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER Number of campaign or dredging site
!
!******************************************************************************
!
!*Local variables
!
INTEGER :: m, nolocs, shipno, siteno
INTEGER, DIMENSION(4) :: nhdims
REAL :: V_subzone_total
REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: V_subzone


procname(pglev+1) = 'determine_campaign_params'
CALL log_timer_in()

!
!1. Iniitalise parameters and arrays
!-----------------------------------
!
!---initialise parameters
nhdims = 0
nolocs = MERGE(nodredgingsites,nocampaigns,iopt_dar_dredging_criterion.EQ.3)
shipno = MERGE(Campaigns(n)%nr_ship,DredgingSites(n)%nr_ship,&
             & iopt_dar_dredging_criterion.EQ.3)
siteno = MERGE(Campaigns(n)%nr_dredgingsite,n,&
             & iopt_dar_dredging_criterion.EQ.3)

!---allocate arrays
IF (iopt_dar_time_dredge.EQ.2) THEN
   ALLOCATE(V_subzone(ncloc,nrloc,SubzonesDredging(siteno)%nr_subzones),&
          & STAT=errstat)
   CALL error_alloc('V_subzone',1,&
               & (/ncloc,nrloc,SubzonesDredging(siteno)%nr_subzones/),kndrtype)
ENDIF

V_subzone = 0.0

!
!2. Amount of dredged volume
!---------------------------
!
!2.1 User-defined dredging
!-------------------------
!

IF (iopt_dar_dredging_criterion.EQ.3) THEN

   SELECT CASE (iopt_dar_user_volume)

!     ---campaign depth
      CASE (1)
         WHERE (mask_dredging(:,:,siteno))
            dredging_depth(:,:,n) = Campaigns(n)%campdepth
            V_dredge(:,:,n) = dredging_depth(:,:,n)*garea
         END WHERE
         CALL sum2_vars(V_dredge(:,:,n),V_dredge_total(n),nhdims,'C  ',0,&
                      & commall=.TRUE.,mask=mask_dredging(:,:,siteno))

!     ---campaign volume
      CASE (2)
         WHERE (mask_dredging(:,:,siteno))
            V_dredge(:,:,n) = Campaigns(n)%campvolume*garea/gareatot
         END WHERE
         CALL sum2_vars(V_dredge(:,:,n),V_dredge_total(n),nhdims,'C  ',0,&
                      & commall=.TRUE.,mask=mask_dredging(:,:,n))
         WHERE (mask_dredging(:,:,siteno))
            dredging_depth(:,:,n) = V_dredge_total(n)/gareatot
         END WHERE

!     ---campaign level    
      CASE (3)
         WHERE (mask_dredging(:,:,n).AND. &
              & depmeanatc(1:ncloc,1:nrloc).LT.Campaigns(n)%camplevel)
            dredging_depth(:,:,n) = Campaigns(n)%camplevel - &
                                  & depmeanatc(1:ncloc,1:nrloc)
            V_dredge(:,:,n) = dredging_depth(:,:,n)*garea
         END WHERE
         CALL sum2_vars(V_dredge(:,:,n),V_dredge_total(n),nhdims,'C  ',0,&
                      & commall=.TRUE.,mask=mask_dredging(:,:,n))
   END SELECT

!
!2.2. Criterium based dredging
!-----------------------------
!

ELSE

   WHERE (mask_dredging(:,:,n))
      dredging_depth(:,:,n) = DredgingSites(n)%target_depth + &
                            & DredgingSites(n)%over_depth - &
                            & depmeanatc(1:ncloc,1:nrloc)
      V_dredge(:,:,n) = dredging_depth(:,:,n)*garea
   END WHERE
ENDIF

!
!3. Full dredging cycle
!----------------------
!

IF (iopt_dar_coupling_time.EQ.3) THEN 

   nr_cycles(n) = CEILING(V_dredge_total(n)/Ships(shipno)%V_ship)
   V_lastcycle(n) = V_dredge_total(n)-(nr_cycles(n)-1)*Ships(shipno)%V_ship
   nt_dredge_lastcycle(n) = CEILING(V_lastcycle(n)*nt_fill(n)/&
                                  & Ships(shipno)%V_ship)
   nt_relocate_lastcycle(n) = CEILING(V_lastcycle(n)*nt_relocate(n)/&
                                    & Ships(shipno)%V_ship)
   nt_dredge(n) = (nr_cycles(n)-1)*nt_fill(n) + nt_dredge_lastcycle(n)
   IF (iopt_dar_time_dredge.EQ.2) THEN 
      nt_passed(n) = 0
      priority_counter(n) = 1
   ENDIF

ENDIF

!
!4. Priority dredging
!--------------------
!

IF (iopt_dar_time_dredge.EQ.2) THEN

   m_410: DO m=1,SubzonesDredging(siteno)%nr_subzones
      WHERE (mask_subzones_dredging(:,:,siteno,m))
         V_subzone(:,:,m) = V_dredge(:,:,n)
      END WHERE
      CALL sum2_vars(V_subzone(:,:,m),V_subzone_total,nhdims,'C  ',0,&
               & commall=.TRUE.,mask=mask_subzones_dredging(:,:,siteno,m))
         nt_dredge_subzone(n,m) = &
               & CEILING(100*V_subzone_total/V_dredge_total(m))*nt_dredge(n)/100
      ENDDO m_410

!  ---avoid rounding errors
   nt_dredge_subzone(n,SubzonesDredging(siteno)%nr_subzones) = nt_dredge(n) - &
 & SUM(nt_dredge_subzone(n,1:SubzonesDredging(siteno)%nr_subzones-1))

ENDIF

!
!5. Track dredging
!-----------------
!

IF (iopt_dar_time_dredge.EQ.3) THEN 
   V_dredge_rest(:,:,n) = V_dredge(:,:,n)
   WHERE (mask_dredging(:,:,siteno))
      V_dredge_max(:,:,n) = DredgingSites(siteno)%dredging_thickness*garea
   END WHERE
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE determine_campaign_params

!==============================================================================

SUBROUTINE determine_dar_times(n)
!******************************************************************************
!
! *determine_dar_times* determines timesteps at which dredging and relocation
!                       takes place
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - check_dar_criterium
!
! External calls - 
!
! Module calls - 
!
!******************************************************************************
!
USE dararrays
USE darswitches
USE grid
USE gridpars
USE iopars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name      Type  Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER counter indicating the campaign or dredging site
!
!******************************************************************************
!
!*Local variables
!
INTEGER :: m


procname(pglev+1) = 'determine_dar_times'
CALL log_timer_in()

!
!1. User-defined dredging
!------------------------
!

IF (iopt_dar_time_dredge.EQ.3) THEN
   OperationTimes(n)%dredging_times(1) = campaign_start(n) + nt_wait(n) + &
                                       & nt_idle(n)
   GOTO 1000
ENDIF

!
!2. Duration of dredging/relocation campaign 
!-------------------------------------------
!
!2.1 Full dredging cycle
!-----------------------
!

IF (iopt_dar_duration.EQ.3.AND.nr_cycles(n).GT.0) THEN

   OperationTimes(n)%dredging_times(1) = nt + nt_wait(n) + nt_idle(n)
   OperationTimes(n)%dredging_times(2) = &
 & OperationTimes(n)%dredging_times(1) + nt_fill(n)
   OperationTimes(n)%relocation_times(1) = &
 & OperationTimes(n)%dredging_times(2) + nt_sail_relocate(n)
   OperationTimes(n)%relocation_times(2) = &
 & OperationTimes(n)%relocation_times(1) + nt_relocate(n)
   IF (nr_cycles(n).EQ.1) THEN
      OperationTimes(n)%campaign_end = OperationTimes(n)%relocation_times(2)
   ELSE
      m_210: DO m=3,2*(nr_cycles(n)-1),2
         OperationTimes(n)%dredging_times(m) = &
       & OperationTimes(n)%relocation_times(m-1) + &
       & nt_sail_dredge(n) + nt_idle(n)
         OperationTimes(n)%dredging_times(m+1) = &
       & OperationTimes(n)%dredging_times(m) + nt_fill(n)
         OperationTimes(n)%relocation_times(m) =  &
       & OperationTimes(n)%dredging_times(m+1) + nt_sail_relocate(n)
         OperationTimes(n)%relocation_times(m+1) = &
       & OperationTimes(n)%relocation_times(m) + nt_relocate(n)
      ENDDO m_210
         ! Final cycle
      OperationTimes(n)%dredging_times(2*nr_cycles(n)-1) = &
    & OperationTimes(n)%dredging_times(2*nr_cycles(n)-3) + nt_cycle(n)
      OperationTimes(n)%dredging_times(2*nr_cycles(n)) = &
    & OperationTimes(n)%dredging_times(2*nr_cycles(n)-1) + &
    & nt_dredge_lastcycle(n)
      OperationTimes(n)%relocation_times(2*nr_cycles(n)-1) =  &
    & OperationTimes(n)%dredging_times(2*nr_cycles(n)) + nt_sail_relocate(n)
      OperationTimes(n)%relocation_times(2*nr_cycles(n)) = &
    & OperationTimes(n)%relocation_times(2*nr_cycles(n)-1) + &
    & nt_relocate_lastcycle(n)
      OperationTimes(n)%campaign_end = &
    & OperationTimes(n)%relocation_times(2*nr_cycles(n))
   ENDIF

!
!2.2 Instantaneous
!-----------------
!

ELSEIF (iopt_dar_duration.EQ.1) THEN 
   OperationTimes(n)%dredging_times(1) = nt
!  ---simultaneous
   IF (iopt_dar_coupling_time.EQ.1) THEN 
      OperationTimes(n)%relocation_times(1) = nt
      OperationTimes(n)%campaign_end = nt
!  ---time lagged
   ELSEIF (iopt_dar_coupling_time.EQ.2) THEN
      OperationTimes(n)%relocation_times(1) = nt + nt_lag(n)
      OperationTimes(n)%campaign_end = nt + nt_lag(n)
   ENDIF

!
!2.3 Spread over time
!--------------------
!

ELSEIF (iopt_dar_duration.EQ.2) THEN
        
!  ---time spread
   OperationTimes(n)%dredging_times (1) = nt
   OperationTimes(n)%dredging_times (2) = nt + nt_spread(n)
!  ---simultaneous
   IF (iopt_dar_coupling_time.EQ.1) THEN 
      OperationTimes(n)%relocation_times = OperationTimes(n)%dredging_times
      OperationTimes(n)%campaign_end = OperationTimes(n)%relocation_times(2)
!  ---time lagged
   ELSEIF (iopt_dar_coupling_time.EQ.2) THEN 
      OperationTimes(n)%relocation_times = &
    & OperationTimes(n)%dredging_times + nt_lag (n)
      OperationTimes(n)%campaign_end = OperationTimes(n)%relocation_times(2)
   ENDIF
ENDIF

1000 CALL log_timer_out()


RETURN

END SUBROUTINE determine_dar_times

!==============================================================================

SUBROUTINE determine_plume_location(i,j,kbot,ktop,zbot,ztop,frac,idir)
!******************************************************************************
!
! *determine_plume_location* define the lower and upper vertical levels
!                            where dredged/relocated material enters the water
!                            column
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description - returns also for the fractional of input sediment for each
!               vertical level
!
! Reference -
!
! Calling program - equal_dredging_priority_dredging, track_dredging, 
!                   relocate
!
!******************************************************************************
!
USE depths
USE gridpars
USE grid

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, idir, j
INTEGER, INTENT(OUT) :: kbot, ktop
REAL, INTENT(IN) :: zbot, ztop
REAL, INTENT(OUT), DIMENSION(nz) :: frac
!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*i*        INTEGER X-index of dredging/relocation location 
!*j*        INTEGER Y-index of dredging/relocation location
!*kbot*     INTEGER Lower vertical level for sediment input
!*ktop*     INTEGER Upper vertical level for sediment input
!*zbot*     REAL    Height above the bed or distance below the surface of the
!                   lowest input level
!*ztop*     REAL    Height above the bed or distance below the surface of the
!                   highest input level
!*frac*     REAL    Fraction of the input dredged/relocated sediment at each
!                   vertical level
!*idir*     INTEGER Defines the vertical direction from which zbot ztop are
!                   taken
!              = 1 => distance from the sea bed
!              = 2 => distance from the sea surface 
!
!*******************************************************************************
!
!*Local variables
!
INTEGER :: k
REAL ::  zlo, zup, z1, z2


!---return if input vertical section is outside the water column
IF ((idir.EQ.1.AND.zbot.GE.deptotatc(i,j)).OR.&
  & (idir.EQ.2.AND.ztop.GE.deptotatc(i,j))) THEN
   frac = 0.0; kbot = 0; ktop = 0
ENDIF

!---convert to height above the bed
z1 = MERGE(zbot,deptotatc(i,j)-zbot,idir.EQ.1)
z2 = MERGE(ztop,deptotatc(i,j)-ztop,idir.EQ.2)

!---lower and upper levels
kbot = 1; zlo = delzatc(i,j,1)
DO WHILE (zlo.GT.z1)
   kbot = kbot + 1
   zlo = zlo + delzatc(i,j,kbot)
ENDDO
ktop = kbot; zup = zlo
DO WHILE (zup.GT.z2.AND.ktop.LE.nz)
   ktop = ktop + 1
   IF (ktop.LE.nz) zup = zup + delzatc(i,j,ktop)
ENDDO
ktop = MIN(ktop,nz)

!---fractions
frac = 0.0
IF (kbot.EQ.ktop) THEN
   frac(kbot) = 1.0
ELSE
   k_210 : DO k=1,nz
      IF (k.LT.kbot.OR.k.GT.ktop) THEN
         frac(k) = 0.0
      ELSEIF (k.EQ.kbot) THEN
         frac(k) = zlo-zbot
      ELSEIF (k.EQ.ktop) THEN
         frac(k) = ztop-zup+delzatc(i,j,ktop)
      ELSE
         frac(k) = delzatc(i,j,k)
      ENDIF
   ENDDO k_210
   frac = frac/(z2-z1)
ENDIF


RETURN

END SUBROUTINE determine_plume_location

!==============================================================================

SUBROUTINE equal_dredging(n)
!******************************************************************************
!
! *equal_dredging* perform dredging using equal spreading in time
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - dar_dredging
!
! External calls - dar_sorting_dredge, determine_plume_location
!
! Module calls - sum2_vars
!
!******************************************************************************
!
USE dararrays
USE darswitches
USE depths
USE gridpars
USE iopars
USE morpharrays
USE morphswitches
USE sedarrays
USE sedpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE paral_utilities, ONLY: sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER Campaign or dredging site number
!
!******************************************************************************
!
!*Local variables
!
INTEGER :: i, j, f, kbot, ktop, shipno, siteno, volume_ship_old
INTEGER, DIMENSION(4) :: nhdims 
REAL :: volume_red
REAL, DIMENSION(nz) :: frac
REAL, DIMENSION(nf) :: fraction_ship_old, fraction_volume
REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: bed_fraction_red_sum, hdredge
REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: bed_fraction_red


procname(pglev+1) = 'equal_dredging'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE(bed_fraction_red(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('bed_fraction_red',3,(/ncloc,nrloc,nf/),kndrtype)
bed_fraction_red = 0.0

ALLOCATE(bed_fraction_red_sum(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bed_fraction_red_sum',2,(/ncloc,nrloc/),kndrtype)
bed_fraction_red_sum = 0.0

ALLOCATE(hdredge(ncloc,nrloc),STAT=errstat)
CALL error_alloc('hdredge',2,(/ncloc,nrloc/),kndrtype)
hdredge = 0.0

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---dredging site annd ship number 
IF (iopt_dar_dredging_criterion.EQ.3) THEN
   shipno = Campaigns(n)%nr_ship
   siteno = Campaigns(n)%nr_dredgingsite
ELSE
   shipno = DredgingSites(n)%nr_ship
   siteno = n
ENDIF

!---work space
fraction_ship_old = fraction_ship(shipno,:)
volume_ship_old = volume_ship(shipno)
SELECT CASE (iopt_dar_duration)
   CASE (1)
      WHERE (mask_dredging(:,:,siteno))
         hdredge = dredging_depth(:,:,n)
      END WHERE
   CASE (2)
      WHERE (mask_dredging(:,:,siteno))
         hdredge = dredging_depth(:,:,n)/nt_spread(n)
      END WHERE
   CASE (3)
      WHERE (mask_dredging(:,:,siteno))
         hdredge = dredging_depth(:,:,n)/nt_dredge(n)
      END WHERE
END SELECT

!
!3. Update volume and fraction distribution on the ship
!------------------------------------------------------
!
!3.1 Calculation of stored volume
!--------------------------------
!

f_310: DO f=1,nf
   WHERE (mask_dredging(:,:,siteno))
      bed_fraction_red(:,:,f) = (1-p_suctionhead(n,f))*(1-p_overflow(n,f))*&
                              & bed_fraction(:,:,1,f)
   END WHERE
ENDDO f_310
   
WHERE (mask_dredging(:,:,siteno))
   bed_fraction_red_sum = SUM(bed_fraction_red,DIM=3)
END WHERE
   
CALL sum2_vars(V_dredge(:,:,n)*bed_fraction_red_sum,volume_red,nhdims,'C  ',&
             & 0,commall=.TRUE.,mask=mask_dredging(:,:,siteno))

volume_ship(shipno) = volume_ship(shipno) + volume_red

!
!3.2 Fraction distribution
!-------------------------
!

f_320: DO f=1,nf
   CALL sum2_vars(V_dredge(:,:,n)*bed_fraction_red(:,:,f),fraction_volume(f),&
                & nhdims,'C  ',0,commall=.TRUE.,mask=mask_dredging(:,:,siteno))
   fraction_ship(shipno,f) = (fraction_ship_old(f)*volume_ship_old + &
                            & fraction_volume(f))/volume_ship(shipno)
ENDDO f_320

!
!4. Morphology
!-------------
!

IF (iopt_dar_effect.EQ.1.OR.iopt_dar_effect.EQ.3) THEN

   WHERE (mask_dredging(:,:,siteno))
      depmeanatc(1:ncloc,1:nrloc) = depmeanatc(1:ncloc,1:nrloc) + hdredge
   END WHERE

   IF (iopt_morph.EQ.1) THEN
      IF(nb.GT.1) THEN
         CALL dar_sorting_dredge(mask_dredging(:,:,n),hdredge)
      ELSEIF (iopt_morph_fixed_layer.EQ.1) THEN
         bed_layer_thickness(:,:,1) = bed_layer_thickness(:,:,1) - hdredge
      ENDIF
   ENDIF

ENDIF

!
!5. Sediment plumes
!------------------
!

IF (iopt_dar_effect.EQ.2 .OR. iopt_dar_effect.EQ.3) THEN

   j_510: DO j = 1,nrloc
   i_510: DO i = 1,ncloc
      IF (mask_dredging(i,j,n)) THEN

!        ---source term at suctionhead
         CALL determine_plume_location(i,j,kbot,ktop,H_dredge_min(n),&
                                     & H_dredge_max(n),frac,1)
         IF (kbot.GT.0) THEN
            f_511: DO f=1,nf
               dar_sediment(i,j,kbot:ktop,f) = dar_sediment(i,j,kbot:ktop,f) + &
                                   & frac(kbot:ktop)*p_suctionhead(n,f)*&
                                   & bed_fraction(i,j,1,f)*hdredge(i,j)/delt3d
            ENDDO f_511
         ENDIF

!        ---overflow source term
         CALL determine_plume_location(i,j,kbot,ktop,D_dredge_max(n),&
                                     & H_dredge_min(n),frac,2)
         IF (kbot.GT.0) THEN
            f_512: DO f=1,nf
               dar_sediment(i,j,kbot:ktop,f) = dar_sediment(i,j,kbot:ktop,f) + &
                    & frac(kbot:ktop)*p_overflow(n,f)*(1.0-p_suctionhead(n,f))*&
                    & bed_fraction(i,j,1,f)*hdredge(i,j)/delt3d
            ENDDO f_512
         ENDIF

      ENDIF

   ENDDO i_510
   ENDDO j_510
      
ENDIF

!
!6. Deallocate arrays
!----- --------------
!

DEALLOCATE (bed_fraction_red,bed_fraction_red_sum,hdredge)

CALL log_timer_out()


RETURN

END SUBROUTINE equal_dredging

!==============================================================================

SUBROUTINE get_spatial_counter(n)
!******************************************************************************
!
! *get_spatial_counter* obtain the spatial counter which sets the relocation
!                       site in case of a user-defined spatial coupling through
!                       time series
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - spatial_coupling
!
! External calls - 
!
! Module calls - 
!
!******************************************************************************
!

USE dararrays
USE darswitches
USE grid
USE gridpars
USE iopars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE
!
!*Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name      Type  Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER counter indicating the campaign or dredging site
!
!******************************************************************************
!

procname(pglev+1) = 'get_spatial_counter'
CALL log_timer_in()


IF (nt.GT.OperationTimes(n)%coupling_times(couple_counter(n))) THEN
   couple_counter(n) = couple_counter(n) + 1
   DO WHILE (nt.GT.OperationTimes(n)%coupling_times(couple_counter(n)).AND.&
           & couple_counter(n).LE.nocouplingtimes(n))
      couple_counter(n) = couple_counter(n) + 1
   ENDDO
ENDIF
couple_counter(n) = couple_counter(n) - 1

IF (iopt_dar_dredging_criterion.EQ.3) THEN
   spatial_counter(n) = DredgingSites(&
                   & Campaigns(n)%nr_dredgingsite)%nr_coupled(couple_counter(n))
ELSE
   spatial_counter(n) = DredgingSites(n)%nr_coupled(couple_counter(n))
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE get_spatial_counter

!===============================================================================

SUBROUTINE priority_dredging(n)
!*******************************************************************************
!
! *priority_dredging* predetermined priority order on dredging site
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - dar_dredging
!
! External calls - 
!
! Module calls - error_alloc, sum2_vars
!
!******************************************************************************
!
USE dararrays
USE darswitches
USE depths
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphswitches
USE sedarrays
USE sedpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE paral_utilities, ONLY: sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER Campaign or dredging site number
!
!******************************************************************************
!
!*Local variables
!
INTEGER :: i, j, f, kbot, ktop, shipno, siteno, subzone, volume_ship_old
INTEGER, DIMENSION(4) :: nhdims
REAL :: volume_red 
REAL, DIMENSION(nz) :: frac
REAL, DIMENSION(nf) :: fraction_ship_old, fraction_volume
REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: bed_fraction_red_sum, hdredge
REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: bed_fraction_red


procname(pglev+1) = 'priority_dredging'
CALL log_timer_in()

!
!1. Allocate arrays
!------------------
!

ALLOCATE(bed_fraction_red_sum(ncloc,nrloc),STAT=errstat)
CALL error_alloc('bed_fraction_red_sum',2,(/ncloc,nrloc/),kndrtype)
bed_fraction_red_sum = 0.0

ALLOCATE(bed_fraction_red(ncloc,nrloc,nf),STAT=errstat)
CALL error_alloc('bed_fraction_red',3,(/ncloc,nrloc,nf/),kndrtype)
bed_fraction_red = 0.0

ALLOCATE(hdredge(ncloc,nrloc),STAT=errstat)
CALL error_alloc('hdredge',2,(/ncloc,nrloc/),kndrtype)
hdredge = 0.0

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---dredging site annd ship number 
IF (iopt_dar_dredging_criterion.EQ.3) THEN
   shipno = Campaigns(n)%nr_ship
   siteno = Campaigns(n)%nr_dredgingsite
ELSE
   shipno = DredgingSites(n)%nr_ship
   siteno = n
ENDIF
subzone = DredgingSites(siteno)%priority_order(priority_counter(n))

!---work space
fraction_ship_old = fraction_ship(shipno,:)
volume_ship_old = volume_ship(shipno)
WHERE (mask_subzones_dredging(:,:,siteno,subzone))
   hdredge = dredging_depth(:,:,n)/nt_dredge_subzone(n,subzone)
END WHERE

!
!3. Update volume and fraction distribution on the ship
!------------------------------------------------------
!
!3.1 Calculation of stored volume
!--------------------------------
!

f_310: DO f=1,nf
   WHERE (mask_subzones_dredging(:,:,siteno,subzone))
      bed_fraction_red(:,:,f) = (1-p_suctionhead(n,f))*(1-p_overflow(n,f))*&
                              & bed_fraction(:,:,1,f)
   END WHERE
ENDDO f_310
   
WHERE (mask_subzones_dredging(:,:,siteno,subzone))
   bed_fraction_red_sum = SUM(bed_fraction_red,DIM=3)
END WHERE
   
CALL sum2_vars(hdredge*garea*bed_fraction_red_sum,volume_red,nhdims,&
             & 'C  ',0,commall=.TRUE.,&
             & mask=mask_subzones_dredging(:,:,siteno,subzone))

volume_ship(shipno) = volume_ship(shipno) + volume_red

!
!3.2 Fraction distribution
!-------------------------
!

f_320: DO f=1,nf
   CALL sum2_vars(hdredge*garea*bed_fraction_red(:,:,f),&
                & fraction_volume(f),nhdims,'C  ',0,commall=.TRUE.,&
                & mask=mask_subzones_dredging(:,:,siteno,subzone))
   fraction_ship(shipno,f) = (fraction_ship_old(f)*volume_ship_old + &
                            & fraction_volume(f))/volume_ship(shipno)
ENDDO f_320

!
!4. Morphology
!-------------
!

IF (iopt_dar_effect.EQ.1 .OR. iopt_dar_effect.EQ.3) THEN

   WHERE (mask_subzones_dredging(:,:,siteno,subzone))
      depmeanatc(1:ncloc,1:nrloc) = depmeanatc(1:ncloc,1:nrloc) + hdredge
   END WHERE

   IF (iopt_morph.EQ.1) THEN
      IF(nb.GT.1) THEN
         CALL dar_sorting_dredge(mask_subzones_dredging(:,:,n,subzone),hdredge)
      ELSEIF (iopt_morph_fixed_layer.EQ.1) THEN
         bed_layer_thickness(:,:,1) = bed_layer_thickness(:,:,1) - hdredge
      ENDIF
   ENDIF
   
   nt_passed(n) = nt_passed(n) + 1
   IF (nt_passed(n).EQ.nt_dredge_subzone(n,subzone)) THEN
      priority_counter(n) = priority_counter(n) + 1
      nt_passed(n) = 0
   ENDIF

ENDIF

!
!5. Sediment plumes
!------------------
!

IF (iopt_dar_effect.EQ.2 .OR. iopt_dar_effect.EQ.3) THEN

   j_510: DO j=1,nrloc
   i_510: DO i=1,ncloc
      IF (mask_subzones_dredging(i,j,siteno,subzone)) THEN

!        --- source term at suctionhead
         CALL determine_plume_location(i,j,kbot,ktop,H_dredge_min(n),&
                                     & H_dredge_max(n),frac,1)
         IF (kbot.GT.0) THEN
            f_511: DO f=1,nf
               dar_sediment(i,j,kbot:ktop,f) = dar_sediment(i,j,kbot:ktop,f) + &
                                   & frac(kbot:ktop)*p_suctionhead(n,f)*&
                                   & bed_fraction(i,j,1,f)*hdredge(i,j)/delt3d
            ENDDO f_511
         ENDIF

!        --- overloss source term
         CALL determine_plume_location(i,j,kbot,ktop,D_dredge_max(n),&
                                     & H_dredge_min(n),frac,2)
         IF (kbot.GT.0) THEN
            f_512: DO f=1,nf
               dar_sediment(i,j,kbot:ktop,f) = dar_sediment(i,j,kbot:ktop,f) + &
                                   & frac(kbot:ktop)*p_suctionhead(n,f)*&
                                   & bed_fraction(i,j,1,f)*hdredge(i,j)/delt3d
            ENDDO f_512
         ENDIF

      ENDIF

   ENDDO i_510
   ENDDO j_510
      
ENDIF
   
!
!6. Deallocate arrays
!--------------------
!

DEALLOCATE (bed_fraction_red,bed_fraction_red_sum,hdredge)

CALL log_timer_out()

RETURN

END SUBROUTINE priority_dredging

!==============================================================================

SUBROUTINE relocate(n,mask)
!******************************************************************************
!
! *relocate* relocation algorithm
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.10
!
! Description -
!
! Reference -
!
! Calling program - dar_relocation
!
! External calls - 
!
! Module calls - sum_vars
!
!******************************************************************************
!
USE dararrays
USE darswitches
USE depths
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphswitches
USE sedarrays
USE sedpars
USE switches
USE timepars
USE paral_utilities, ONLY: sum2_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Arguments
!
INTEGER, INTENT (IN) :: n
LOGICAL, DIMENSION(ncloc,nrloc), INTENT(IN) :: mask

! Name      Type  Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER Campaign or dredging site number
!*mask*     LOGICAL .TRUE. for cells satisfying the relocation criterion
!
!******************************************************************************
!
!*Local variables
!
INTEGER :: i, j, f, kbot, ktop, shipno, siteno, time_relocate
INTEGER, DIMENSION(4) :: nhdims
REAL :: deposition, gareatot_rel, update, volume_ship_rel
REAL, DIMENSION(nz) :: frac
REAL, DIMENSION(nf) :: fraction_rel


procname(pglev+1) = 'relocate'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

IF (iopt_dar_dredging_criterion.EQ.3) THEN
   shipno = Campaigns(n)%nr_ship
   siteno = Campaigns(n)%relocationsites(spatial_counter(n))
ELSE
   shipno = DredgingSites(n)%nr_ship
   siteno = DredgingSites(n)%relocationsites(spatial_counter(n))
 ENDIF

nhdims = 0

!
!2. Ship volume and fractions
!----------------------------
!

IF (.NOT.Campaigns(n)%dredge) THEN
   IF (iopt_dar_coupling_time.EQ.3) THEN
      IF( cycle_counter(n).LT.nr_cycles(n)) THEN
         volume_ship(shipno) = Ships(shipno)%V_ship
      ELSE
         volume_ship(shipno) = V_lastcycle(n)
      ENDIF
   ELSE
      volume_ship(shipno) = V_dredge_total(n)
   ENDIF
   fraction_ship(shipno,:) = Campaigns(n)%campfractions(1:nf)
ENDIF

!
!3. Relocated volume at current time step
!----------------------------------------
!

SELECT CASE (iopt_dar_duration)
!  ---instantaneous
   CASE (1); volume_ship_rel = volume_ship(shipno)
!  ---spread over time
   CASE (2); volume_ship_rel = volume_ship(shipno)/nt_spread(n)
!  ---spread over cycles
   CASE (3)
      IF (cycle_counter(n).EQ.nr_cycles(n)) THEN
         time_relocate = nt_relocate_lastcycle(n)
      ELSE
         time_relocate = nt_relocate(n)
      ENDIF
      volume_ship_rel = volume_ship(shipno)/time_relocate
END SELECT

!
!4. Morhpology
!------------
!

IF (iopt_dar_effect.EQ.1.OR.iopt_dar_effect.EQ.3) THEN

   fraction_rel = (1-p_relocation(n,:))*(1-p_passive(n,:))*&
                & fraction_ship(shipno,:)
   CALL sum2_vars(garea,gareatot_rel,nhdims,'C  ',0,commall=.TRUE.,mask=mask)
!  ---update bathymetry
   deposition = SUM(fraction_rel)*volume_ship_rel/gareatot_rel
   WHERE (mask)
      depmeanatc(1:ncloc,1:nrloc) = depmeanatc(1:ncloc,1:nrloc) - deposition
   END WHERE
!  ---amount of relocated sediment
   V_relocate_total(siteno) = V_relocate_total(siteno) + &
                            & SUM(fraction_rel)*volume_ship_rel
!  ---amount of non-relocated sediment
   IF (iopt_dar_relocation_criterion.EQ.1 .AND. &
        & (.NOT.RelocationSites(siteno)%eligibility)) THEN
      V_relocate_extra(siteno) = V_relocate_extra(siteno) + volume_ship(shipno)
   ELSEIF (iopt_dar_relocation_criterion.EQ.2) THEN
      V_relocate_extra(siteno) = MAX(V_relocate_total(siteno)-&
                                   & RelocationSites(siteno)%max_volume,0.0)
   ENDIF
 ! ---vertical sorting
   IF (iopt_morph.EQ.1) THEN
      IF (nb.GT.1) THEN
         update = volume_ship_rel/gareatot
         CALL dar_sorting_relocate(mask,fraction_rel,update)
      ELSEIF (iopt_morph_fixed_layer.EQ.1) THEN
         WHERE(mask)
            bed_layer_thickness(:,:,1) = bed_layer_thickness(:,:,1) + deposition
         END WHERE
      ENDIF
   ENDIF

ENDIF

!
!5. Sediment plumes
!------------------
!

IF (iopt_dar_effect.EQ.2.OR.iopt_dar_effect.EQ.3) THEN
      
   j_510: DO j=1,nrloc
   i_510: DO i=1,ncloc
      IF (mask(i,j)) THEN

!        ---passive plume source term
         CALL determine_plume_location(i,j,kbot,ktop,H_relocate_min(n),&
                                     & H_relocate_max(n),frac,1)
         IF (kbot.GT.0) THEN
            f_511: DO f=1,nf
               dar_sediment(i,j,kbot:ktop,f) = dar_sediment(i,j,kbot:ktop,f) + &
                       & p_relocation(n,f)*frac(kbot:ktop)*&
                       & fraction_ship(shipno,f)*volume_ship_rel/delt3d
            ENDDO f_511
         ENDIF

!        ---overflow source term
         CALL determine_plume_location(i,j,kbot,ktop,D_relocate_min(n),&
                                     & D_relocate_max(n),frac,1)
         IF (kbot.GT.0) THEN
            f_512: DO f=1,nf
               dar_sediment(i,j,kbot:ktop,f) = dar_sediment(i,j,kbot:ktop,f) + &
                      & p_passive(n,f)*(1.0-p_relocation(n,f))*frac(kbot:ktop)*&
                      & fraction_ship(shipno,f)*volume_ship_rel/delt3d
            ENDDO f_512
         ENDIF

      ENDIF

   ENDDO i_510
   ENDDO j_510
      
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE relocate

!==============================================================================

SUBROUTINE set_dar_flags(n)
!******************************************************************************
!
! *set_dar_flags* decides about dredging/relocation activities at the current
!                 time step
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description - determine whether dredging and/or relocation activities 
!               start/end on the current timestep
!             - define which relocation site is used
!
! Reference -
!
! Calling program - dar_equation
!
! External calls - spatial_coupling
!
! Module calls - 
!
!******************************************************************************
!
USE dararrays
USE darswitches
USE iopars
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE
!
!*Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name      Type  Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER counter indicating the campaign or dredging site
!
!******************************************************************************
!
!*Local variables
!
LOGICAL :: dredge, relocate
INTEGER :: shipno, norelsites


procname(pglev+1) = 'set_dar_flags'
CALL log_timer_in()

!
!1. Check if dredging and/or relocation functionality is enabled
!---------------------------------------------------------------
!

IF (iopt_dar_dredging_criterion.EQ.3) THEN
   dredge = Campaigns(n)%dredge
   relocate = Campaigns(n)%relocate
   shipno = Campaigns(n)%nr_ship
ELSE
   dredge = .TRUE.
   relocate = .TRUE.
   shipno = DredgingSites(n)%nr_ship
ENDIF

!
!2. Start dredging
!-----------------
!

IF (dredge.AND.(.NOT.dredging(n)).AND.&
  & OperationTimes(n)%dredging_times(dredge_counter(n)).EQ.nt) THEN
   dredging(n) = .TRUE.
   dredge_counter(n) = dredge_counter(n) + 1
   volume_ship(shipno) = 0.0
   fraction_ship(shipno,:) = 0.0
!  ---simultaneous dredging and relocation
   IF (relocate.AND.iopt_dar_coupling_time.EQ.1.AND.&
     & OperationTimes(n)%relocation_times(relocate_counter(n)).EQ.nt) THEN
      relocation(n) = .TRUE.
      relocate_counter(n) = relocate_counter(n) + 1
   ENDIF
ENDIF

!
!3. Start relocation
!-------------------
!

IF (relocate.AND.(.NOT.relocation(n)).AND.&
  & OperationTimes(n)%relocation_times(relocate_counter(n)).EQ.nt) THEN
   relocation(n) = .TRUE.
   relocate_counter(n) = relocate_counter(n) + 1
   IF (iopt_dar_duration.EQ.3.AND.eligibility(n)) CALL spatial_coupling(n)
ENDIF

!
!4. End dredging
!---------------
!

IF (dredge.AND.dredging(n).AND.&
  & OperationTimes(n)%dredging_times(dredge_counter(n)).EQ.nt) THEN
   dredging(n) = .FALSE.
   dredge_counter(n) = dredge_counter(n) + 1
ENDIF

!
!5. End relocation
!-----------------
!

IF (relocate.AND.relocation(n).AND.&
  & OperationTimes(n)%relocation_times(relocate_counter(n)).EQ.nt) THEN
   relocation(n) = .FALSE.
   relocate_counter(n) = relocate_counter(n) + 1
   cycle_counter(n) = cycle_counter(n) + 1
!  ---switching to a new relocation site in case of cyclical shifting
!     (spatial coupling)
   IF (iopt_dar_coupling_space.EQ.1) THEN
      IF (iopt_dar_dredging_criterion.EQ.3) THEN
         norelsites = Campaigns(n)%nr_zones_coupled
      ELSE
         norelsites = DredgingSites(n)%nr_zones_coupled
      ENDIF
!     --cyclical shifting --> automatically shifts to a new site
      IF (spatial_counter(n).EQ.norelsites) THEN
         spatial_counter(n) = 1
      ELSE
         spatial_counter(n) = spatial_counter(n)+1
      ENDIF
   ENDIF
ENDIF

!
!6. End of the campaign
!----------------------
!

IF (OperationTimes(n)%campaign_end.EQ.nt) THEN
   campaign(n) = .FALSE.; couple_counter(n) = 1; dredge_counter(n) = 1
   relocate_counter(n) = 1; cycle_counter(n) = 1
   priority_counter(n) = 1; track_cell_counter(n) = 1
ENDIF

CALL log_timer_out()

RETURN

END SUBROUTINE set_dar_flags

!==============================================================================

SUBROUTINE spatial_coupling(n)
!******************************************************************************
!
! *spatial_coupling* identifies the designated relocation site based
!                    on one of three spatial coupling criteria
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! The relocation site the ships sails to is chosen from a list of possible
! relocation sites. The spatial counter acts as an index to this list and
! depends on the switch iopt_dar_coupling_space.
!  = 1 --> relocation site changes every cycle (spatial_counter ++)
!  = 2 --> relocation site only changes when entering the site is prohibited
!          (determined by criterion)
!  = 3 --> relocation site determined from timeseries
!
! If iopt_dar_coupling_space equals 3, the spatial coupling is governed by a
! time series, in which certain timesteps are linked to certain relocation
! sites. So in order to know which relocation site to sail to, coherens needs
! to link the current timestep to the time series. The spatial counter is
! defined by iopt_dar_coupling_space
!  = 1 --> Spatial counter is transferred from last timestep and updated at the
!          end of the routine
!  = 2 --> Spatial counter is transferred from last timestep and not updated.
!  = 3 --> Spatial counter is obtained from timeseries and not updated.
!
! If iopt_dar_coupling_space not equals 3, check if the site is eligible for
! relocation. If not, shift to the next site in line and check for eligibility.
! If false, continue shifting until an eligible site is found or until all
! sites are found not eligible. In that case, stick with the current site and
! label all the relocation sites as not eligible (i.e. the site never shifts
! anymore)
!
! Reference -
!
! Calling program - set_dar_flags
!
! External calls - check_eligibility, get_spatial_counter 
!
! Module calls -
!
!*******************************************************************************
!
USE dararrays
USE darswitches
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name      Type  Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER counter indicating the campaign or dredging site
!
!******************************************************************************
!
!*Local variables
!
INTEGER :: norelsites, siteno, site_counter


procname(pglev+1) = 'spatial_coupling'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

IF (iopt_dar_dredging_criterion.EQ.3) THEN
   norelsites = Campaigns(n)%nr_zones_coupled
ELSE
   norelsites = DredgingSites(n)%nr_zones_coupled
ENDIF

siteno = MERGE(Campaigns(n)%relocationsites(spatial_counter(n)),&
             & DredgingSites(n)%relocationsites(spatial_counter(n)),&
             & iopt_dar_dredging_criterion.EQ.3)

site_counter = 0

!
!2. Set spatial counter
!----------------------
!
!2.1 From time series
!--------------------
!

IF (iopt_dar_coupling_space.EQ.3) THEN
   CALL get_spatial_counter(n)

!
!2.2 Explore for eligible relocation site
!----------------------------------------
!

ELSE
   CALL check_eligibility(siteno)

   IF (.NOT.RelocationSites(siteno)%eligibility)THEN
      
      DO WHILE (.NOT.RelocationSites(siteno)%eligibility.AND.eligibility(n))
            IF (spatial_counter(n).EQ.norelsites) THEN
               spatial_counter(n) = 1
            ELSE
               spatial_counter(n) = spatial_counter(n) + 1
            ENDIF
            site_counter = site_counter + 1
            IF (site_counter.EQ.norelsites) eligibility(n) = .FALSE.
            siteno = MERGE(&
                 & Campaigns(n)%relocationsites(spatial_counter(n)),&
                 & DredgingSites(n)%relocationsites(spatial_counter(n)),&
                 & iopt_dar_dredging_criterion.EQ.3)
            CALL check_eligibility(siteno)
         ENDDO

      ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE spatial_coupling

!==============================================================================

SUBROUTINE track_dredging(n)
!******************************************************************************
!
! *track_dredging* dredging along a dredging track
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - dar_dredging
!
! External calls - 
!
! Module calls - copy_vars, local_proc
!
!******************************************************************************
!
USE dararrays
USE darswitches
USE depths
USE grid
USE gridpars
USE iopars
USE morpharrays
USE morphswitches
USE paralpars
USE sedarrays
USE sedpars
USE switches
USE timepars
USE paral_comms, ONLY: copy_vars
USE grid_routines, ONLY: local_proc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n

!
! Name      Type  Purpose
!------------------------------------------------------------------------------
!*n*        INTEGER counter indicating the campaign or dredging site
!
!******************************************************************************
!
!*Local variables
!
INTEGER :: f, i, iglb, j, jglb, kbot, ktop, shipno, siteno 
REAL :: bed_fraction_red_sum, hdredge, V_dredge_pass, volume_red, &
      & volume_ship_old
REAL, DIMENSION(nz) :: frac
REAL, DIMENSION(nf) :: bed_fraction_red, fraction_ship_old 


procname(pglev+1) = 'track_dredging'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

IF (iopt_dar_dredging_criterion.EQ.3) THEN
   siteno = Campaigns(n)%nr_dredgingsite
   shipno = Campaigns(n)%nr_ship
ELSE
   siteno = n
   shipno = DredgingSites(n)%nr_ship
ENDIF

!
!2. Perform dredging
!-------------------
!

iglb = itrack(track_cell_counter(n),siteno)
jglb = jtrack(track_cell_counter(n),siteno)
i = iglb - nc1loc + 1; j = jglb - nr1loc + 1

IF (local_proc(i,j)) THEN

!
!2.1 Amount of dredging
!----------------------
!

   IF (V_dredge_rest(i,j,n).GT.V_dredge_max(i,j,n).AND. &
     & V_ship_rest(n).GT.V_dredge_max(i,j,n)) THEN
      V_dredge_pass = V_dredge_max(i,j,n)
   ELSEIF(V_dredge_rest(i,j,n).GT.V_dredge_max(i,j,n).AND. &
        & V_dredge_max(i,j,n).GT.V_ship_rest(n)) THEN
      V_dredge_pass = V_ship_rest(n)
   ELSEIF(V_dredge_rest(i,j,n).LT.V_dredge_max(i,j,n).AND. &
        & V_dredge_rest(i,j,n).LT.V_ship_rest(n)) THEN
      V_dredge_pass = V_dredge_rest(i,j,n)
   ELSE
      V_dredge_pass = V_ship_rest(n)
   ENDIF
   hdredge = V_dredge_pass/(garea(i,j)*nt_track(track_cell_counter(n),n))

!   
!2.2 Morphology
!--------------
!

   IF (iopt_dar_effect.EQ.1.OR.iopt_dar_effect.EQ.3) THEN
      depmeanatc(i,j) = depmeanatc(i,j) + hdredge
   ENDIF
   IF (iopt_morph.EQ.1) THEN
      IF (nb.GT.1) THEN
         CALL dar_sorting_dredge(mask_dredging(i,j,n),hdredge)
      ELSEIF (iopt_morph_fixed_layer.EQ.1) THEN
         bed_layer_thickness(:,:,1) = bed_layer_thickness(:,:,1) - hdredge
      ENDIF
   ENDIF

!
!2.3 Sediment plumes
!-------------------
!

   IF (iopt_dar_effect.EQ.2.OR.iopt_dar_effect.EQ.3) THEN

!     ---source term at suctionhead
      CALL determine_plume_location(i,j,kbot,ktop,H_dredge_min(n),&
                                  & H_dredge_max(n),frac,1)
      IF (kbot.GT.0) THEN
         f_231: DO f=1,nf
            dar_sediment(i,j,kbot:ktop,f) = dar_sediment(i,j,kbot:ktop,f) + &
                    & frac(kbot:ktop)*p_overflow(n,f)*(1.0-p_suctionhead(n,f))*&
                    & bed_fraction(i,j,1,f)*hdredge/delt3d
         ENDDO f_231
      ENDIF

!     ---overflow source term
      CALL determine_plume_location(i,j,kbot,ktop,D_dredge_max(n),&
                                  & H_dredge_min(n),frac,2)
      IF (kbot.GT.0) THEN
         f_232: DO f=1,nf
            dar_sediment(i,j,kbot:ktop,f) = dar_sediment(i,j,kbot:ktop,f) + &
                    & frac(kbot:ktop)*p_overflow(n,f)*(1.0-p_suctionhead(n,f))*&
                    & bed_fraction(i,j,1,f)*hdredge/delt3d
         ENDDO f_232
      ENDIF
   ENDIF

!   
!2.4 Stored sediment volume and fraction distribution on the ship
!----------------------------------------------------------------
!

   IF (V_dredge_pass.GT.0.0) THEN
      fraction_ship_old = fraction_ship(shipno,:)
      volume_ship_old = volume_ship(shipno)
      
!     ---stored volume
      bed_fraction_red = (1-p_suctionhead(n,:))*(1-p_overflow(n,:))*&
                       & bed_fraction(i,j,1,:)
      bed_fraction_red_sum = SUM(bed_fraction_red)
      volume_red = bed_fraction_red_sum*V_dredge_pass
      volume_ship(shipno) = volume_ship(shipno) + volume_red
      
!     ---fraction distribution
      fraction_ship(shipno,:) = (fraction_ship_old(:)*volume_ship_old+&
              & bed_fraction_red*volume_red)/volume_ship(shipno)
   ENDIF

!   ---adjusting volume of the ship and the site volume matrix
   V_dredge_rest(i,j,n) = V_dredge_rest(i,j,n) - V_dredge_pass
   V_ship_rest(n) = V_ship_rest(n) - V_dredge_pass

ENDIF

!
!3. Copy to other processes
!--------------------------
!

IF (iopt_MPI.EQ.1) CALL copy_vars(V_ship_rest(n),0,idroot=idloc)

!
!4. Advancing the ship along the track
!-------------------------------------
!

IF (track_nt_counter(n).EQ.nt_track(track_cell_counter(n),n)) THEN
   track_nt_counter(n) = 1
   IF (track_cell_counter(n).NE.track_points_counter(n)) THEN
!     ---ship has not completed the track   
      track_cell_counter(n) = track_cell_counter(n) + 1
!     ---ship has reached the end of the track -> go back to beginning of track
   ELSE
      track_cell_counter(n) = 1
   ENDIF
ELSE
   track_nt_counter(n) =  track_nt_counter(n) + 1
ENDIF

!
!5. Updating operation times
!----------------------------
!

IF (V_ship_rest(n).EQ.0.0) THEN
!  --- end of dredging + start and end of relocation
   OperationTimes(n)%dredging_times(dredge_counter(n)) = nt + 1
   OperationTimes(n)%relocation_times(dredge_counter(n)-1) =  &
 & OperationTimes(n)%dredging_times(dredge_counter(n)) + nt_sail_relocate(n)
   OperationTimes(n)%relocation_times(dredge_counter(n)) = &
 & OperationTimes(n)%relocation_times(dredge_counter(n)-1) + nt_relocate(n)
!  ---no sediment left, the campaign ends -> set end of campaign
   IF (.NOT.ANY(V_dredge_rest(:,:,n).GT.1E-08.AND.&
     & mask_dredging(:,:,siteno))) THEN
     OperationTimes(n)%campaign_end = &
   & OperationTimes(n)%relocation_times(dredge_counter(n))
  ELSE
!    ---still sediment left -> reset ship volume and start of new dredging cycle
     V_ship_rest(n)= Ships(shipno)%V_ship
     OperationTimes(n)%dredging_times(dredge_counter(n)+1) = &
   & OperationTimes(n)%relocation_times(dredge_counter(n)) + &
   & nt_sail_dredge(n) + nt_idle(n)
   ENDIF
!   ---ship is not full, but no sediment left to dredge -> set end of dredging, 
!   ---start and end of relocation and end of campaign
ELSEIF (.NOT.ANY(V_dredge_rest(:,:,n).GT.1E-08.AND.&
      & mask_dredging(:,:,siteno))) THEN
   OperationTimes(n)%dredging_times(dredge_counter(n)) = nt + 1
   OperationTimes(n)%relocation_times(dredge_counter(n)-1) =  &
 & OperationTimes(n)%dredging_times(dredge_counter(n)) + nt_sail_relocate(n)
   OperationTimes(n)%relocation_times(dredge_counter(n)) = &
 & OperationTimes(n)%relocation_times(dredge_counter(n)-1) + nt_relocate(n)
   OperationTimes(n)%campaign_end = &
 & OperationTimes(n)%relocation_times(dredge_counter(n))
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE track_dredging

!===============================================================================

SUBROUTINE Zcoord_dredge(zcoord,lbounds,ubounds)
!*******************************************************************************
!
! *Zcoord_dredge* Array of Z-coordinates
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description - returns zero at dry cells or nodes
!
! Calling program - time_series_grid
!
!*******************************************************************************
!

USE darpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE
!
!*  Arguments
!
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1), &
     & lbounds(2):ubounds(2),lbounds(3):ubounds(3)) :: zcoord

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*zcoord*    REAL    Array of Z-coordinates of the dredging arrays

!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, l1, l2, l3, u1, u2, u3

procname(pglev+1) = 'Zcoord_dredge'
CALL log_timer_in()

!
!1. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!2. Calculation
!---------------
!

k_210: DO k = 2,nodredgingsites
   zcoord(l1:u1,l2:u2,k) = k
ENDDO k_210

CALL log_timer_out()


RETURN

END SUBROUTINE Zcoord_dredge

!==============================================================================

SUBROUTINE Zcoord_relocate(zcoord,lbounds,ubounds)
!******************************************************************************
!
! *Zcoord_relocate* Array of Z-coordinates
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Dredging_Equations.f90  V2.8
!
! Description - returns zero at dry cells or nodes
!
! Calling program - time_series_grid
!
!******************************************************************************
!

USE darpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE
!
!*  Arguments
!
INTEGER, INTENT(IN), DIMENSION(3) :: lbounds, ubounds
REAL, INTENT(OUT), DIMENSION(lbounds(1):ubounds(1), &
     & lbounds(2):ubounds(2),lbounds(3):ubounds(3)) :: zcoord

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*zcoord*    REAL    Array of Z-coordinates of the relocation arrays

!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: k, l1, l2, l3, u1, u2, u3

procname(pglev+1) = 'Zcoord_relocate'
CALL log_timer_in()

!
!1. Array bounds
!---------------
!

l1 = lbounds(1); l2 = lbounds(2); l3 = lbounds(3)
u1 = ubounds(1); u2 = ubounds(2); u3 = ubounds(3)

!
!2. Calculation
!---------------
!

k_210: DO k = 2,norelocationsites
   zcoord(l1:u1,l2:u2,k) = k
ENDDO k_210

CALL log_timer_out()


RETURN

END SUBROUTINE Zcoord_relocate
