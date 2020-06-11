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

MODULE dartypes_init
!************************************************************************
!
! *dartypes_init* Initialise derived type scalars and arrays for the
!                 dredging/relocation module
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)dartypes_init.f90  V2.8
!
! $Date: 2015-01-16 16:39:38 +0100 (Fri, 16 Jan 2015) $
!
! $Revision: 796 $
!
! Description - 
!
! Routines - campaigns_init, dredgingsites_init, error_alloc_dar_struc,
!            operationtimes_init, relocationsites_init, ships_init,
!            subzonesdredging_init, subzonesrelocation_init
!
!************************************************************************
!

CONTAINS

!========================================================================

SUBROUTINE campaigns_init(Campaigns)
!************************************************************************
!
! *campaigns_init* Initialise derived type variable of type 'CampaignParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)dartypes_init.f90  V2.8
!
! Description - argument of rank 1
!
!************************************************************************
!
USE dartypes

!
!*Arguments
!
TYPE (CampaignParams), INTENT(OUT), DIMENSION(:) :: Campaigns

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*Campaigns*  DERIVED Attributes of dredging and relocation campaigns
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize

nsize = SIZE(Campaigns)
IF (nsize.EQ.0) RETURN

n_110: DO n=1,nsize
   Campaigns(n)%dredge= .TRUE.
   Campaigns(n)%relocate = .TRUE.
   Campaigns(n)%nr_relocationsite = 1
   Campaigns(n)%nr_zones_coupled = 1
   Campaigns(n)%nr_dredgingsite = 1
   Campaigns(n)%nr_ship = 1
   Campaigns(n)%t_spread = 1
   Campaigns(n)%t_lag = 1
   Campaigns(n)%campvolume = 0.0
   Campaigns(n)%camplevel = 0.0
   Campaigns(n)%campdepth = 1.0
   Campaigns(n)%CampaignStartDateTime = ''
   Campaigns(n)%relocationsites = 1
   Campaigns(n)%CoupleTime = ''
   Campaigns(n)%campfractions = 0.0
ENDDO n_110


RETURN

END SUBROUTINE campaigns_init

!========================================================================

SUBROUTINE dredgingsites_init(DredgingSites)
!************************************************************************
!
! *dredgingsites_init* Initialise derived type variable of type 'DredgingParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)dartypes_init.f90  V2.8
!
! Description - argument of rank 1
!
!************************************************************************
!
USE dartypes

!
!*Arguments
!
TYPE (DredgingParams), INTENT(OUT), DIMENSION(:) :: DredgingSites

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*DredgingSites*  DERIVED Attributes of dredging sites
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize

nsize = SIZE(DredgingSites)
IF (nsize.EQ.0) RETURN

n_110: DO n=1,nsize
   DredgingSites(n)%nr_coordpol = 1
   DredgingSites(n)%nr_ship= 1
   DredgingSites(n)%nr_zones_coupled = 1
   DredgingSites(n)%t_spread = 1
   DredgingSites(n)%t_lag = 1
   DredgingSites(n)%target_depth = 0.0
   DredgingSites(n)%over_depth = 0.0
   DredgingSites(n)%critical_volume = 0.0
   DredgingSites(n)%dredging_thickness = 0.0
   DredgingSites(n)%monitoring_frequency = 1
   DredgingSites(n)%relocationsites = 1
   DredgingSites(n)%priority_order = 1
   DredgingSites(n)%Xpol = 1.0
   DredgingSites(n)%Ypol = 1.0
   DredgingSites(n)%nr_coupled = 1
   DredgingSites(n)%CoupleTime = ''
ENDDO n_110


RETURN

END SUBROUTINE dredgingsites_init

!========================================================================

SUBROUTINE error_alloc_dar_struc(arrname,ndims,nshape,structype,abort)
!************************************************************************
!
! *error_alloc_dar_struc* Displays error message when a dynamic
!                         dredging/relocation array of derived type cannot be
!                         allocated
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)error_routines.F90  V2.8
!
! Description -
!
! Module calls - error_abort, error_arg_var
!
!************************************************************************
!
USE darpars
USE iopars
USE syspars
USE error_routines, ONLY: nerrs, error_abort, error_arg_var

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: abort
CHARACTER (LEN=*), INTENT(IN) :: arrname, structype
INTEGER, INTENT(IN) :: ndims
INTEGER, INTENT(IN), DIMENSION(ndims) :: nshape

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*arrname*   CHAR    Array name
!*ndims*     INTEGER Array rank
!*nshape*    INTEGER Array shape
!*structype* CHAR    Type of derived structure
!*abort*     LOGICAL Abort program if PRESENT and .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!

CHARACTER (LEN=12) :: cerrstat
CHARACTER (LEN=22) :: csize
INTEGER :: n
INTEGER (KIND=kndilong) :: nsize
CHARACTER (LEN=12), DIMENSION(ndims) :: cshape


IF (errstat.NE.0) THEN

   SELECT CASE (TRIM(structype))
      CASE ('Campaigns')
         nsize = kndlog*2 + kndchar*lentime*(1+MaxTimes) + &
               & kndint*(6+MaxSites) +  kndreal*(3+MaxFractions)
      CASE ('DredgingParams') 
         nsize = kndchar*lentime*MaxTimes + &
               & kndint*(6+MaxTimes+MaxSites+MaxSubzones) + &
               & kndreal*(4+2*MaxCoordPol)
      CASE ('OperationTimesParams')
         nsize = kndint*(1+3*MaxTimes)
      CASE ('RelocationParams')
         nsize = kndlog + kndint + kndreal*(2+2*MaxCoordPol)
      CASE ('ShipParams')
         nsize = kndint*6 + kndreal*(11+MaxFractions*4)
      CASE ('SubzonesDredgingParams')
         nsize = kndint*(1+MaxSubzones) + kndreal*2*MaxCoordPol*MaxSubzones
      CASE ('SubzonesRelocationParams')
         nsize = kndint*(1+MaxSubzones) + kndreal*2*MaxCoordPol*MaxSubzones
      CASE DEFAULT
         procname(pglev+1) = 'error_alloc_struc'
         CALL error_arg_var(structype,'structype')
   END SELECT
   nsize = nsize*PRODUCT(nshape)

   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (ioerr,'(2A)') 'Unable to allocate derived type array: ', arrname
      n_110: DO n=1,ndims
         WRITE (cshape(n),'(I12)') nshape(n)
         cshape(n) = ADJUSTL(cshape(n))
      ENDDO n_110
      WRITE (ioerr,'(8A)') 'Shape:', (' '//TRIM(cshape(n)),n=1,ndims)
      WRITE (csize,'(I12)') nsize
      csize = ADJUSTL(csize)
      WRITE (ioerr,'(3A)') 'Size: ', TRIM(csize), ' bytes'
      WRITE (cerrstat,'(I12)') errstat
      cerrstat = ADJUSTL(cerrstat)
      WRITE (ioerr,'(A)') 'Exit code '//TRIM(cerrstat)
   ENDIF
   
   IF (.NOT.PRESENT(abort)) THEN
      CALL error_abort(procname(pglev),ierrno_alloc)
   ELSEIF (abort) THEN
      CALL error_abort(procname(pglev),ierrno_alloc)
   ENDIF

ENDIF


RETURN

END SUBROUTINE error_alloc_dar_struc

!========================================================================

SUBROUTINE operationtimes_init(OperationTimes)
!************************************************************************
!
! *operationtimes_init*  Initialise derived type variable of type
!                       'OperationTimesParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)dartypes_init.f90  V2.8
!
! Description - argument of rank 1
!
!************************************************************************
!
USE dartypes

!
!*Arguments
!
TYPE (OperationTimesParams), INTENT(OUT), DIMENSION(:) :: OperationTimes

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*OperationTimes*  DERIVED Attributes of dar timesteps
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize

nsize = SIZE(OperationTimes)
IF (nsize.EQ.0) RETURN

n_110: DO n=1,nsize
   OperationTimes(n)%campaign_end = 0
   OperationTimes(n)%coupling_times= 0
   OperationTimes(n)%dredging_times= 0
   OperationTimes(n)%relocation_times = 0
ENDDO n_110


RETURN

END SUBROUTINE operationtimes_init

!========================================================================

SUBROUTINE relocationsites_init(RelocationSites)
!************************************************************************
!
! *relocationsites_init* Initialise derived type variable of type
!                       'RelocationParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)dartypes_init.f90  V2.8
!
! Description - argument of rank 1
!
!************************************************************************
!
USE dartypes

!
!*Arguments
!
TYPE (RelocationParams), INTENT(OUT), DIMENSION(:) :: RelocationSites

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*RelocationSites*  DERIVED Attributes of relocation sites
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize

nsize = SIZE(RelocationSites)
IF (nsize.EQ.0) RETURN

n_110: DO n=1,nsize

   RelocationSites(n)%nr_coordpol = 1
   RelocationSites(n)%eligibility = .TRUE.
   RelocationSites(n)%min_depth = 1.0
   RelocationSites(n)%max_volume = 1.0
   RelocationSites(n)%Xpol = 1.0
   RelocationSites(n)%Ypol = 1.0

ENDDO n_110


RETURN

END SUBROUTINE relocationsites_init

!========================================================================

SUBROUTINE ships_init(Ships)
!************************************************************************
!
! *ships_init* Initialise derived type variable of type 'ShipParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)dartypes_init.f90  V2.8
!
! Description - argument of rank 1
!
!************************************************************************
!
USE dartypes

!
!*Arguments
!
TYPE (ShipParams), INTENT(OUT), DIMENSION(:) :: Ships

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*Ships*  DERIVED Attributes of dredging ships
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize

nsize = SIZE(Ships)
IF (nsize.EQ.0) RETURN

n_110: DO n=1,nsize

   Ships(n)%U_ship=1.0
   Ships(n)%V_ship= 1.0
   Ships(n)%t_fill= 1
   Ships(n)%t_idle= 1
   Ships(n)%t_relocate= 1
   Ships(n)%t_sail_dredge= 1
   Ships(n)%t_sail_relocate= 1
   Ships(n)%t_wait= 1
   Ships(n)%D_dredge_max= 0.0
   Ships(n)%D_dredge_min= 0.0
   Ships(n)%D_relocate_max = 0.0
   Ships(n)%D_relocate_min = 0.0
   Ships(n)%H_dredge_max = 0.0
   Ships(n)%H_dredge_min = 0.0
   Ships(n)%H_relocate_max = 0.0
   Ships(n)%H_relocate_min = 0.0
   Ships(n)%R_circle= 1.0
   Ships(n)%p_suctionhead= 0.0
   Ships(n)%p_overflow= 0.0
   Ships(n)%p_relocation= 0.0
   Ships(n)%p_passive= 0.0

ENDDO n_110


RETURN

END SUBROUTINE ships_init

!========================================================================

SUBROUTINE subzonesdredging_init(SubzonesDredging)
!************************************************************************
!
! *subzonesdredging_init* Initialise derived type variable of type
!                        'SubzonesDredgingParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)dartypes_init.f90  V2.8
!
! Description - argument of rank 1
!
!************************************************************************
!
USE dartypes

!
!*Arguments
!
TYPE (SubzonesDredgingParams), INTENT(OUT), DIMENSION(:) :: SubzonesDredging

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*SubzonesDredging*  DERIVED Attributes of dredging subzones
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize

nsize = SIZE(SubzonesDredging)
IF (nsize.EQ.0) RETURN

n_110: DO n=1,nsize

   SubzonesDredging(n)%nr_subzones= 1
   SubzonesDredging(n)%nr_coordpol = 1
   SubzonesDredging(n)%Xpol = 1.0
   SubzonesDredging(n)%Ypol = 1.0

ENDDO n_110


RETURN

END SUBROUTINE subzonesdredging_init

!========================================================================

SUBROUTINE subzonesrelocation_init(SubzonesRelocation)
!************************************************************************
!
! *subzonesrelocation_init* Initialise derived type variable of type
!                          'SubzonesRelocationParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)dartypes_init.f90  V2.8
!
! Description - argument of rank 1
!
!************************************************************************
!
USE dartypes

!
!*Arguments
!
TYPE (SubzonesRelocationParams), INTENT(OUT), DIMENSION(:) :: SubzonesRelocation

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*SubzonesRelocation*  DERIVED Attributes of relocation subzones

!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize

nsize = SIZE(SubzonesRelocation)
IF (nsize.EQ.0) RETURN

n_110: DO n=1,nsize

   SubzonesRelocation(n)%nr_subzones= 1
   SubzonesRelocation(n)%nr_coordpol= 1
   SubzonesRelocation(n)%Xpol = 1.0
   SubzonesRelocation(n)%Ypol = 1.0

ENDDO n_110


RETURN

END SUBROUTINE subzonesrelocation_init


END MODULE dartypes_init
