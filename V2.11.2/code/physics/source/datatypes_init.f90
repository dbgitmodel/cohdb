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

MODULE datatypes_init
!************************************************************************
!
! *datatypes_init* Initialise derived type scalars and arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.11.2
!
! $Date: 2018-05-31 16:41:05 +0200 (Thu, 31 May 2018) $
!
! $Revision: 1141 $
!
! Description - 
!
! Generic routines - filepars_init, hrelativecoords_init, varatts_init,
!                    vrelativecoords_init
!
! Routines - exchangecomms_init, gridpars_init, outgpars_init, statslocs_init
!
!************************************************************************
!
USE datatypes

INTERFACE filepars_init
   MODULE PROCEDURE filepars_init_0d, filepars_init_1d, filepars_init_2d, &
                  & filepars_init_3d
END INTERFACE

INTERFACE hrelativecoords_init
   MODULE PROCEDURE hrelativecoords_init_0d, hrelativecoords_init_1d, &
                  & hrelativecoords_init_2d
END INTERFACE

INTERFACE varatts_init
   MODULE PROCEDURE varatts_init_0d, varatts_init_1d, varatts_init_2d
END INTERFACE

INTERFACE vrelativecoords_init
   MODULE PROCEDURE vrelativecoords_init_0d, vrelativecoords_init_1d, &
                  & vrelativecoords_init_2d, vrelativecoords_init_3d, &
                  & vrelativecoords_init_4d
END INTERFACE

CONTAINS

!========================================================================

SUBROUTINE exchangecomms_init(comms)
!************************************************************************
!
! *exchangecomms_init* Type 'ExchComms'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.0
!
! Description - argument of rank 1
!
!************************************************************************
!
!*Arguments
!F
TYPE (ExchComms), INTENT(OUT), DIMENSION(:) :: comms

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*comms*     DERIVED Parameters for halo communications
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize


nsize = SIZE(comms)
IF (nsize.EQ.0) RETURN

n_110: DO n=1,nsize
   comms(n)%sfirst = .FALSE.
   comms(n)%iddest = 0
   comms(n)%idsrce = 0
   comms(n)%tag = 0
   comms(n)%i1rcv = 0
   comms(n)%i2rcv = 0
   comms(n)%j1rcv = 0
   comms(n)%j2rcv = 0
   comms(n)%i1snd = 0
   comms(n)%i2snd = 0
   comms(n)%j1snd = 0
   comms(n)%j1snd = 0
ENDDO n_110


RETURN

END SUBROUTINE exchangecomms_init

!========================================================================

SUBROUTINE filepars_init_0d(filepars)
!************************************************************************
!
! *filepars_init_0d* Type 'FileParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.11
!
! Description - argument of rank 0
!
!************************************************************************
!
USE iopars

!
!*Arguments
!
TYPE (FileParams), INTENT(OUT) :: filepars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED File attributes
!
!------------------------------------------------------------------------------
!


filepars%defined = .FALSE.
filepars%fill = .FALSE.
filepars%info = .FALSE.
filepars%opened = .FALSE.
filepars%time_regular = .FALSE.
filepars%floattype = ''
filepars%form = ''
filepars%status = '0'
filepars%filename = ''
filepars%pathname = ''
filepars%title = ''
filepars%comment = ''
filepars%endfile = 0
filepars%iunit = int_fill
filepars%maxrecs = 0
filepars%nocoords = 0
filepars%nodim = 0
filepars%novars = 0
filepars%iostat = 0
filepars%timeid = 0
filepars%timerec = 0
filepars%tskips = 0
filepars%varid = 0
filepars%zvarid = 0
filepars%nstepout = 0
filepars%tlims = 0


RETURN

END SUBROUTINE filepars_init_0d

!========================================================================

SUBROUTINE filepars_init_1d(filepars)
!************************************************************************
!
! *filepars_init_1d* Type 'FileParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.11
!
! Description - argument of rank 1
!
!************************************************************************
!
USE iopars

!
!*Arguments
!
TYPE (FileParams), INTENT(OUT), DIMENSION(:) :: filepars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED File attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize


nsize = SIZE(filepars)
IF (nsize.EQ.0) RETURN

n_110: DO n=1,nsize
   filepars(n)%defined = .FALSE.
   filepars(n)%fill = .FALSE.
   filepars(n)%info = .FALSE.
   filepars(n)%opened = .FALSE.
   filepars(n)%packing = .FALSE.
   filepars(n)%time_regular = .FALSE.
   filepars(n)%floattype = '' 
   filepars(n)%form = ''
   filepars(n)%status = '0'
   filepars(n)%filename = ''
   filepars(n)%pathname = ''
   filepars(n)%title = ''
   filepars(n)%comment = ''
   filepars(n)%endfile = 0
   filepars(n)%iunit = int_fill
   filepars(n)%maxrecs = 0
   filepars(n)%nocoords = 0
   filepars(n)%nodim = 0
   filepars(n)%novars = 0
   filepars(n)%iostat = 0
   filepars(n)%timeid = 0
   filepars(n)%timerec = 0
   filepars(n)%tskips = 0
   filepars(n)%varid = 0
   filepars(n)%zvarid = 0
   filepars(n)%nstepout = 0
   filepars(n)%tlims = 0
ENDDO n_110


RETURN

END SUBROUTINE filepars_init_1d

!========================================================================

SUBROUTINE filepars_init_2d(filepars)
!************************************************************************
!
! *filepars_init_2d* Type 'FileParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.11
!
! Description - argument of rank 2
!
!************************************************************************
!
USE iopars

!
!*Arguments
!
TYPE (FileParams), INTENT(OUT), DIMENSION(:,:) :: filepars
!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED File attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n1, n2
INTEGER, DIMENSION(2) :: ndims


IF (SIZE(filepars).EQ.0) RETURN
ndims = SHAPE(filepars)

n2_110: DO n2=1,ndims(2)
n1_110: DO n1=1,ndims(1)
   filepars(n1,n2)%defined = .FALSE.
   filepars(n1,n2)%fill = .FALSE.
   filepars(n1,n2)%info = .FALSE.
   filepars(n1,n2)%opened = .FALSE.
   filepars(n1,n2)%packing = .FALSE.
   filepars(n1,n2)%time_regular = .FALSE.
   filepars(n1,n2)%floattype = '' 
   filepars(n1,n2)%form = ''
   filepars(n1,n2)%status = '0'
   filepars(n1,n2)%filename = ''
   filepars(n1,n2)%pathname = ''
   filepars(n1,n2)%title = ''
   filepars(n1,n2)%comment = ''
   filepars(n1,n2)%endfile = 0
   filepars(n1,n2)%iunit = int_fill
   filepars(n1,n2)%maxrecs = 0
   filepars(n1,n2)%nocoords = 0
   filepars(n1,n2)%nodim = 0
   filepars(n1,n2)%novars = 0
   filepars(n1,n2)%iostat = 0
   filepars(n1,n2)%timeid = 0
   filepars(n1,n2)%timerec = 0
   filepars(n1,n2)%tskips = 0
   filepars(n1,n2)%varid = 0
   filepars(n1,n2)%nstepout = 0
   filepars(n1,n2)%zvarid = 0
   filepars(n1,n2)%tlims = 0
ENDDO n1_110
ENDDO n2_110


RETURN

END SUBROUTINE filepars_init_2d

!========================================================================

SUBROUTINE filepars_init_3d(filepars)
!************************************************************************
!
! *filepars_init_3d* Type 'FileParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.11
!
! Description - argument of rank 3
!
!************************************************************************
!  
USE iopars
  
!
!*Arguments
!
TYPE (FileParams), INTENT(OUT), DIMENSION(:,:,:) :: filepars
!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED File attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n1, n2, n3
INTEGER, DIMENSION(3) :: ndims


IF (SIZE(filepars).EQ.0) RETURN
ndims = SHAPE(filepars)

n3_110: DO n3=1,ndims(3)
n2_110: DO n2=1,ndims(2)
n1_110: DO n1=1,ndims(1)
   filepars(n1,n2,n3)%defined = .FALSE.
   filepars(n1,n2,n3)%info = .FALSE.
   filepars(n1,n2,n3)%opened = .FALSE.
   filepars(n1,n2,n3)%packing = .FALSE.
   filepars(n1,n2,n3)%time_regular = .FALSE.
   filepars(n1,n2,n3)%floattype = '' 
   filepars(n1,n2,n3)%form = ''
   filepars(n1,n2,n3)%status = '0'
   filepars(n1,n2,n3)%filename = ''
   filepars(n1,n2,n3)%pathname = ''
   filepars(n1,n2,n3)%title = ''
   filepars(n1,n2,n3)%comment = ''
   filepars(n1,n2,n3)%endfile = 0
   filepars(n1,n2,n3)%iunit = int_fill
   filepars(n1,n2,n3)%maxrecs = 0
   filepars(n1,n2,n3)%nocoords = 0
   filepars(n1,n2,n3)%nodim = 0
   filepars(n1,n2,n3)%novars = 0
   filepars(n1,n2,n3)%iostat = 0
   filepars(n1,n2,n3)%timeid = 0
   filepars(n1,n2,n3)%timerec = 0
   filepars(n1,n2,n3)%varid = 0
   filepars(n1,n2,n3)%nstepout = 0
   filepars(n1,n2,n3)%zvarid = 0
   filepars(n1,n2,n3)%tlims = 0
ENDDO n1_110
ENDDO n2_110
ENDDO n3_110


RETURN

END SUBROUTINE filepars_init_3d

!========================================================================

SUBROUTINE gridpars_init(gridpars)
!************************************************************************
!
! *gridpars_init* Type 'GridParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.11.2
!
! Description - argument of rank 1
!
!************************************************************************
!
!*Arguments
!
TYPE (GridParams), INTENT(OUT), DIMENSION(:,:) :: gridpars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*gridpars*  DERIVED Attributes of surface grids
!
!------------------------------------------------------------------------------
!


IF (SIZE(gridpars).EQ.0) RETURN

gridpars%rotated = .FALSE.
gridpars%extrapol = .FALSE.
gridpars%land = .TRUE.
gridpars%datflag = 0
gridpars%nhtype = 0
gridpars%n1dat = 1
gridpars%n2dat = 1
gridpars%surcoords = 1
gridpars%delxdat = 0.0
gridpars%delydat = 0.0
gridpars%gridangle = 0.0
gridpars%longpole = 0.0
gridpars%x0dat = 0.0
gridpars%y0dat = 0.0
gridpars%x0rot = 0.0
gridpars%y0rot = 0.0


RETURN

END SUBROUTINE gridpars_init

!========================================================================

SUBROUTINE hrelativecoords_init_0d(hgrid,flag_undef)
!************************************************************************
!
! *hrelativecoords_init_0d* Type 'HRelativeCoords'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.9
!
! Description - initialise to undefined values if flag_undef is .TRUE.,
!               to zero otherwise
!             - argument of rank 0
!
!************************************************************************
!
USE iopars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: flag_undef
TYPE (HRelativeCoords), INTENT(OUT) :: hgrid

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*hgrid*      DERIVED Horizontal relative coordinates
!*flag_undef* LOGICAL Flag values if .TRUE.
!
!------------------------------------------------------------------------------
!


IF (flag_undef) THEN
   hgrid%icoord = int_fill
   hgrid%jcoord = int_fill
   hgrid%weights = float_fill
ELSE
   hgrid%icoord = 0
   hgrid%jcoord = 0
   hgrid%weights = 0.0
ENDIF


RETURN

END SUBROUTINE hrelativecoords_init_0d

!========================================================================

SUBROUTINE hrelativecoords_init_1d(hgrid,flag_undef)
!************************************************************************
!
! *hrelativecoords_init_1d* Type 'HRelativeCoords'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.9
!
! Description - initialise to undefined values if flag_undef is .TRUE.,
!               to zero otherwise
!             - argument of rank 1
!
!************************************************************************
!
USE iopars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: flag_undef
TYPE (HRelativeCoords), INTENT(OUT), DIMENSION(:) :: hgrid

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*hgrid*      DERIVED Horizontal relative coordinates
!*flag_undef* LOGICAL Flag values if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize

nsize = SIZE(hgrid)
IF (nsize.EQ.0) RETURN

IF (flag_undef) THEN
   hgrid%icoord = int_fill
   hgrid%jcoord = int_fill
   n_110: DO n=1,nsize
      hgrid(n)%weights = float_fill
   ENDDO n_110
ELSE
   hgrid%icoord = 0
   hgrid%jcoord = 0
   n_120: DO n=1,nsize
      hgrid(n)%weights= 0.0
   ENDDO n_120
ENDIF


RETURN

END SUBROUTINE hrelativecoords_init_1d

!========================================================================

SUBROUTINE hrelativecoords_init_2d(hgrid,flag_undef)
!************************************************************************
!
! *hrelativecoords_init_2d* Type 'HRelativeCoords'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.9
!
! Description - initialise to undefined values if flag_undef is .TRUE.,
!               to zero otherwise
!             - argument of rank 2
!
!************************************************************************
!
USE iopars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: flag_undef
TYPE (HRelativeCoords), INTENT(OUT), DIMENSION(:,:) :: hgrid

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*hgrid*      DERIVED Horizontal relative coordinates
!*flag_undef* LOGICAL Flag values if .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n1, n2
INTEGER, DIMENSION(2) :: ndims


IF (SIZE(hgrid).EQ.0) RETURN
ndims = SHAPE(hgrid)

IF (flag_undef) THEN
   hgrid%icoord = int_fill
   hgrid%jcoord = int_fill
   n2_110: DO n2=1,ndims(2)
   n1_110: DO n1=1,ndims(1)
      hgrid(n1,n2)%weights = float_fill
   ENDDO n1_110
   ENDDO n2_110
ELSE
   hgrid%icoord = 0
   hgrid%jcoord = 0
   n2_120: DO n2=1,ndims(2)
   n1_120: DO n1=1,ndims(1)
      hgrid(n1,n2)%weights = 0.0
   ENDDO n1_120
   ENDDO n2_120
ENDIF


RETURN

END SUBROUTINE hrelativecoords_init_2d

!========================================================================

SUBROUTINE outgpars_init(outgpars)
!************************************************************************
!
! *outgpars_init* Type 'OutGridParams'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.10.2
!
! Description - argument of rank 1
!
!************************************************************************
!
!*Arguments
!
TYPE (OutGridParams), INTENT(OUT), DIMENSION(:) :: outgpars

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*outgpars*   DERIVED Space/time locations for user-defined output
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: l, n, nsize


nsize = SIZE(outgpars)
IF (nsize.EQ.0) RETURN

n_110: DO n=1,nsize
   outgpars(n)%dredgevars = .FALSE.
   outgpars(n)%relocatevars = .FALSE.
   outgpars(n)%gridded = .FALSE.
   outgpars(n)%packing = .FALSE.
   outgpars(n)%time_grid = .FALSE.
   outgpars(n)%enddate = ''
   outgpars(n)%refdate = ''
   outgpars(n)%startdate = ''
   outgpars(n)%nodim = 0
   outgpars(n)%nostats = 0
   outgpars(n)%nowetout = 0
   outgpars(n)%nrout = 0
   outgpars(n)%nstepout = 0
   outgpars(n)%nzout = 0
   outgpars(n)%time_format = 1
   outgpars(n)%vcoord = 0
   l_111: DO l=1,3
      outgpars(n)%tlims(l) = 0
      outgpars(n)%xlims(l) = 0
      outgpars(n)%ylims(l) = 0
      outgpars(n)%zlims(l) = 0
   ENDDO l_111
   outgpars(n)%deltout = 0.0
ENDDO n_110


RETURN

END SUBROUTINE outgpars_init

!========================================================================

SUBROUTINE statlocs_init(statlocs)
!************************************************************************
!
! *statlocs_init* Type 'StationLocs'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.0
!
! Description - argument of rank 1
!
!************************************************************************
!
!*Arguments
!
TYPE (StationLocs), INTENT(OUT), DIMENSION(:) :: statlocs

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*outlocs*  DERIVED Station locations for irregular user-defined output
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER :: n, nostats


nostats  = SIZE(statlocs)
IF (nostats.EQ.0) RETURN

n_110: DO n=1,nostats
   statlocs(n)%ipos = 0
   statlocs(n)%jpos = 0
   statlocs(n)%name = ''
ENDDO n_110


RETURN

END SUBROUTINE statlocs_init

!========================================================================

SUBROUTINE varatts_init_0d(varatts)
!************************************************************************
!
! *varatts_init_0d* Type 'VariableAtts'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.7.1
!
! Description - argument of rank 0
!
!************************************************************************
!
USE iopars
USE switches

!*Arguments
!
TYPE (VariableAtts), INTENT(OUT) :: varatts

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*varatts*  DERIVED Variable attributes
!
!------------------------------------------------------------------------------
!


varatts%fill = .FALSE.
varatts%axis = ''
varatts%f90_name = ''
varatts%comment = ''
varatts%coordinates = ''
varatts%long_name = ''
varatts%standard_name = ''
varatts%vector_name = ''
varatts%units = ''
varatts%node = 'C'
varatts%data_type = 0
varatts%ivarid = 0
varatts%klev = 0
varatts%nrank = -1
varatts%numvar = -1
varatts%oopt = oopt_null
varatts%global_dims = 0
varatts%halo_dims = 0
varatts%local_dims = 0
varatts%dimids = 0
varatts%dep = 0.0
varatts%fill_value = double_fill


RETURN

END SUBROUTINE varatts_init_0d

!========================================================================

SUBROUTINE varatts_init_1d(varatts)
!************************************************************************
!
! *varatts_init_1d* Type 'VariableAtts'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.7.1
!
! Description - argument of rank 1
!
!************************************************************************
!
USE iopars
USE switches

!
!*Arguments
!
TYPE (VariableAtts), INTENT(OUT), DIMENSION(:) :: varatts

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*varatts*  DERIVED Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize


nsize = SIZE(varatts)
IF (nsize.EQ.0) RETURN

n_110: DO n=1,nsize
   varatts(n)%fill = .FALSE.
   varatts(n)%axis = ''
   varatts(n)%f90_name = ''
   varatts(n)%comment = ''
   varatts(n)%coordinates = ''
   varatts(n)%long_name = ''
   varatts(n)%standard_name = ''
   varatts(n)%vector_name = ''
   varatts(n)%units = ''
   varatts(n)%node = 'C'
   varatts(n)%data_type = 0
   varatts(n)%ivarid = 0
   varatts(n)%klev = 0
   varatts(n)%nrank = -1
   varatts(n)%numvar = -1
   varatts(n)%oopt = oopt_null
   varatts(n)%global_dims = 0
   varatts(n)%halo_dims = 0
   varatts(n)%local_dims = 0
   varatts(n)%dimids = 0
   varatts(n)%dep = 0.0
   varatts(n)%fill_value = double_fill
ENDDO n_110


RETURN

END SUBROUTINE varatts_init_1d

!========================================================================

SUBROUTINE varatts_init_2d(varatts)
!************************************************************************
!
! *varatts_init_2d* Type 'VariableAtts'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.7.1
!
! Description - argument of rank 2
!
!************************************************************************
!
USE iopars
USE switches

!
!*Arguments
!
TYPE (VariableAtts), INTENT(OUT), DIMENSION(:,:) :: varatts

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*varatts*  DERIVED Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n1, n2
INTEGER, DIMENSION(2) :: ndims


IF (SIZE(varatts).EQ.0) RETURN
ndims = SHAPE(varatts)

n2_110: DO n2=1,ndims(2)
n1_110: DO n1=1,ndims(1)
   varatts(n1,n2)%fill = .FALSE.
   varatts(n1,n2)%axis = ''
   varatts(n1,n2)%f90_name = ''
   varatts(n1,n2)%comment = ''
   varatts(n1,n2)%coordinates = ''
   varatts(n1,n2)%long_name = ''
   varatts(n1,n2)%standard_name = ''
   varatts(n1,n2)%vector_name = ''
   varatts(n1,n2)%units = ''
   varatts(n1,n2)%node = 'C'
   varatts(n1,n2)%data_type = 0
   varatts(n1,n2)%ivarid = 0
   varatts(n1,n2)%klev = 0
   varatts(n1,n2)%nrank = -1
   varatts(n1,n2)%numvar = -1
   varatts(n1,n2)%oopt = oopt_null
   varatts(n1,n2)%global_dims = 0
   varatts(n1,n2)%halo_dims = 0
   varatts(n1,n2)%local_dims = 0
   varatts(n1,n2)%dimids = 0
   varatts(n1,n2)%dep = 0.0
   varatts(n1,n2)%fill_value = double_fill
ENDDO n1_110
ENDDO n2_110


RETURN

END SUBROUTINE varatts_init_2d

!========================================================================

SUBROUTINE vrelativecoords_init_0d(vgrid,flag_undef)
!************************************************************************
!
! *vrelativecoords_init_0d* Type 'VRelativeCoords'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.9
!
! Description - initialise to undefined values if flag_undef is .TRUE.,
!               to zero otherwise
!             - argument of rank 0
!
!************************************************************************
!
USE iopars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: flag_undef
TYPE (VRelativeCoords), INTENT(OUT) :: vgrid

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*vgrid*      DERIVED Vertical relative coordinates
!*flag_undef* LOGICAL Flag values if .TRUE.
!
!------------------------------------------------------------------------------
!


IF (flag_undef) THEN
   vgrid%kcoord = int_fill
   vgrid%weights = float_fill
ELSE
   vgrid%kcoord = 0
   vgrid%weights = 0.0
ENDIF


RETURN

END SUBROUTINE vrelativecoords_init_0d

!========================================================================

SUBROUTINE vrelativecoords_init_1d(vgrid,flag_undef)
!************************************************************************
!
! *vrelativecoords_init_1d* Type 'VRelativeCoords'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.9
!
! Description - initialise to undefined values if flag_undef is .TRUE.,
!               to zero otherwise
!             - argument of rank 1
!
!************************************************************************
!
USE iopars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: flag_undef
TYPE (VRelativeCoords), INTENT(OUT), DIMENSION(:) :: vgrid

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*vgrid*      DERIVED Vertical relative coordinates
!*flag_undef* LOGICAL Flag values if .TRUE.
!
!------------------------------------------------------------------------------
!
!
!*Local variables
!
INTEGER :: n, nsize


nsize = SIZE(vgrid)
IF (nsize.EQ.0) RETURN

IF (flag_undef) THEN
   vgrid%kcoord = int_fill
   n_110: DO n=1,nsize
      vgrid%weights(n) = float_fill
   ENDDO n_110
ELSE
   vgrid%kcoord = 0
   n_120: DO n=1,nsize
      vgrid%weights(n) = 0.0
   ENDDO n_120
ENDIF


RETURN

END SUBROUTINE vrelativecoords_init_1d

!========================================================================

SUBROUTINE vrelativecoords_init_2d(vgrid,flag_undef)
!************************************************************************
!
! *vrelativecoords_init_2d* Type 'VRelativeCoords'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.9
!
! Description - initialise to undefined values if flag_undef is .TRUE.,
!               to zero otherwise
!             - argument of rank 2
!
!************************************************************************
!
USE iopars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: flag_undef
TYPE (VRelativeCoords), INTENT(OUT), DIMENSION(:,:) :: vgrid

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*vgrid*      DERIVED Vertical relative coordinates
!*flag_undef* LOGICAL Flag values if .TRUE.
!
!------------------------------------------------------------------------------
!


IF (SIZE(vgrid).EQ.0) RETURN

IF (flag_undef) THEN
   vgrid%kcoord = int_fill
   vgrid%weights(1) = float_fill
   vgrid%weights(2) = float_fill
ELSE
   vgrid%kcoord = 0
   vgrid%weights(1) = 0.0
   vgrid%weights(2) = 0.0
ENDIF


RETURN

END SUBROUTINE vrelativecoords_init_2d

!========================================================================

SUBROUTINE vrelativecoords_init_3d(vgrid,flag_undef)
!************************************************************************
!
! *vrelativecoords_init_3d* Type 'VRelativeCoords'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.9
!
! Description - initialise to undefined values if flag_undef is .TRUE.,
!               to zero otherwise
!             - argument of rank 3
!
!************************************************************************
!
USE iopars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: flag_undef
TYPE (VRelativeCoords), INTENT(OUT), DIMENSION(:,:,:) :: vgrid

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*vgrid*      DERIVED Vertical relative coordinates
!*flag_undef* LOGICAL Flag values if .TRUE.
!
!------------------------------------------------------------------------------
!


IF (SIZE(vgrid).EQ.0) RETURN 

IF (flag_undef) THEN
   vgrid%kcoord = int_fill
   vgrid%weights(1) = float_fill
   vgrid%weights(2) = float_fill
ELSE
   vgrid%kcoord = 0
   vgrid%weights(1) = 0.0
   vgrid%weights(2) = 0.0
ENDIF


RETURN

END SUBROUTINE vrelativecoords_init_3d

!========================================================================

SUBROUTINE vrelativecoords_init_4d(vgrid,flag_undef)
!************************************************************************
!
! *vrelativecoords_init_4d* Type 'VRelativeCoords'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)datatypes_init.f90  V2.9
!
! Description - initialise to undefined values if flag_undef is .TRUE.,
!               to zero otherwise
!             - argument of rank 4
!
!************************************************************************
!
USE iopars

!
!*Arguments
!
LOGICAL, INTENT(IN) :: flag_undef
TYPE (VRelativeCoords), INTENT(OUT), DIMENSION(:,:,:,:) :: vgrid

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*vgrid*      DERIVED Vertical relative coordinates
!*flag_undef* LOGICAL Flag values if .TRUE.
!
!------------------------------------------------------------------------------
!


IF (SIZE(vgrid).EQ.0) RETURN 

IF (flag_undef) THEN
   vgrid%kcoord = int_fill
   vgrid%weights(1) = float_fill
   vgrid%weights(2) = float_fill
ELSE
   vgrid%kcoord = 0
   vgrid%weights(1) = 0.0
   vgrid%weights(2) = 0.0
ENDIF


RETURN

END SUBROUTINE vrelativecoords_init_4d


END MODULE datatypes_init
