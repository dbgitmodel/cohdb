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

MODULE obconds
!************************************************************************
!
! *obconds* Arrays for 2-D open boundary conditions and 1-D surface forcing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)obconds.f90  V2.7
!
! $Date: 2017-08-25 15:58:53 +0200 (Fri, 25 Aug 2017) $
!
! $Revision: 1046 $
!
! Description - 
!
! Modifications:  Addition new arrays for the new boundary conditions
!                 types, each new array corresponds to the new boundary 
!                 condition type 
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!---2-D open boundary conditions
INTEGER, ALLOCATABLE, DIMENSION(:) :: ityp2dobu, ityp2dobv, ityp2dobx, ityp2doby
INTEGER, ALLOCATABLE, DIMENSION(:) :: iloczobu, iloczobv
INTEGER, ALLOCATABLE, DIMENSION(:) :: iobc2dtype, no2dobuv, no2dobxy
INTEGER, ALLOCATABLE, DIMENSION(:) :: iprofobx, iprofoby, itypobx, itypoby
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: index2dobuv, index2dobxy
REAL, ALLOCATABLE, DIMENSION(:) :: zdatobu_old, zdatobv_old
REAL, ALLOCATABLE, DIMENSION(:,:) :: udatobu, udatoby, vdatobv, vdatobx, &
                                   & zdatobu, zdatobv
REAL, ALLOCATABLE, DIMENSION(:,:) :: uveloby, vvelobx
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: obcsalatu, obcsalatv
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: obctmpatu, obctmpatv
REAL, ALLOCATABLE, DIMENSION(:,:) :: obc2uvatu, obc2uvatu_old, &
                                   & obc2uvatv, obc2uvatv_old
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: obc3uvatu, obc3uvatv
REAL, ALLOCATABLE, DIMENSION(:,:) :: &
                         & udatobu_amp, udatobu_pha, udatoby_amp, udatoby_pha, &
                         & vdatobv_amp, vdatobv_pha, vdatobx_amp, vdatobx_pha, &
                         & zdatobu_amp, zdatobu_pha, zdatobv_amp, zdatobv_pha
!---discharge along sections
INTEGER :: nqsecobu, nqsecobv
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iqsecobu, jqsecobv

!---1-D surface forcing
INTEGER :: isur1dtype
REAL :: gxslope = 0.0, gyslope = 0.0
REAL, ALLOCATABLE, DIMENSION(:) :: gxslope_amp, gyslope_amp, zeta_amp
REAL, ALLOCATABLE, DIMENSION(:) :: gxslope_pha, gyslope_pha, zeta_pha

!---Thatcher-Harleman open boundary condition
REAL, DIMENSION(MaxVLevels) :: return_time = 0.0
REAL, ALLOCATABLE, DIMENSION(:,:) :: floutobu, floutobv


SAVE

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*floutobu*     REAL    Time since last ouflow at U-open boundaries (used by
!                       the Thatcher-Harleman open boundary condition)       [s]
!*floutobv*     REAL    Time since last ouflow at V-open boundaries (used by
!                       the Thatcher-Harleman open boundary condition)       [s]
!*gxslope*      REAL    X-component of pressure gradient                 [m/s^2]
!*gxslope_amp*  REAL    Amplitudes of X-component of pressure gradient   [m/s^2]
!*gxslope_pha*  REAL    Phases of X-component of pressure gradient       [m/s^2]
!*gyslope*      REAL    Y-component of pressure gradient                 [m/s^2]
!*gyslope_amp*  REAL    Amplitudes of Y-component of pressure gradient   [m/s^2]
!*gyslope_pha*  REAL    Phases of Y-component of pressure gradient       [m/s^2]
!*iloczobu*     INTEGER Position of elevation points used for o.b.c. at U-nodes
!                       (0/2)
!                 =  0 => not specified
!                 =  1 => at U-open boundary
!                 =  2 => at external C-node
!*iloczobv*     INTEGER Position of elevation points used for o.b.c. at V-nodes
!                      (0/2)
!*index2dobuv*  INTEGER Index maps to store 2-D input into full array
!                       (indices for U-nodes are between 1 and nobu,
!                       indices for V-nodes are between nobu+1 and nobu+nobv)
!*index2dobxy* INTEGER Index maps to store 2-D input into full array
!                       (indices for X-nodes are between 1 and nobx,
!                       indices for Y-nodes are between nobx+1 and nobx+noby)
!*iobc2dtype*   INTEGER Type of data in 2-D data files
!                 =  1  => transports and elevations
!                 =  2  => elevations
!                 =  3  => transports
!*iqsecobu*     INTEGER Start/end o.b. index of each discharge sections along
!                       U-o.b.
!*isur1dtype*   INTEGER Type of data in 1-D surface data files
!                 =  1  => slopes and elevation
!                 =  2  => elevation
!                 =  3  => slopes
!*iprofobx*     INTEGER Profile number for 3-D tangential o.b.c. at X-nodes
!*iprofoby*     INTEGER Profile number for 3-D tangential o.b.c. at Y-nodes
!*itypobx*      INTEGER Type of tangential o.b.c. at X-nodes for 3-D currents
!                 =  0  => flux is set to its nearest interior value (default)
!                 =  1  => zero gradient condition
!                 =  2  => using external data
!*itypoby*      INTEGER Type of tangential o.b.c. at Y-nodes for 3-D currents
!                       (0/2)
!*ityp2dobu*    INTEGER Type of 2-D o.b.c. at U-nodes (0/17)
!                 =  0  => clamped
!                 =  1  => zero surface slope and U from local solution
!                 =  2  => zero gradient for U, elevation specified at internal
!                         node
!                 =  3  => specified elevation and U from local solution
!                 =  4  => specified transport
!                 =  5  => radiation condition using shallow water speed
!                 =  6  => Orlanski
!                 =  7  => Camerlengo and O'Brien
!                 =  8  => Flather condition with specified elevation and U 
!                 =  9  => Flather condition with specified elevation and U from
!                         local solution
!                 = 10  => Roed and Smedstad
!                 = 11  => characteristic method using specified elevation and U
!                 = 12  => characteristic method using specified elevation and U
!                          from local solution
!                 = 13  => characteristic method using zero gradient condition
!                 = 14  => specified depth mean current
!                 = 15  => specified discharge (at local position)
!                 = 16  => discharge specified along o.b. sections
!                 = 17  => specified surface slope
!*ityp2dobv*    INTEGER Type of 2-D o.b.c. at V-nodes (0/17)
!*ityp2dobx*    INTEGER Type of 2-D o.b.c. at X-nodes (0/2)
!                 =  0  => flux is set to its nearest interior value (default)
!                 =  1  => zero gradient condition
!                 =  2  => using external data
!*ityp2doby*    INTEGER Type of 2-D o.b.c. at Y-nodes (0/2)
!*jqsecobv*     INTEGER Start/end o.b. index of each discharge sections along
!                       V-o.b.
!*no2dobuv*     INTEGER Number of input data at (U,V)-o.b.for each data file
!*no2dobxy*     INTEGER Number of input data at (X,Y)-o.b.for each data file
!*nqsecobu*     INTEGER Number of discharge sections along U-open boundaries
!*nqsecobv*     INTEGER Number of discharge sectionsat along V-open boundaries
!*obcsalatu*    REAL    Storage array for salinity radiation o.b.c. at U-nodes
!*obcsalatv*    REAL    Storage array for salinity radiation o.b.c. at V-nodes
!*obctmpatu*    REAL    Storage array for temperature radiation o.b.c. at
!                       U-nodes
!*obctmpatv*    REAL    Storage array for temperature radiation o.b.c. at
!                       V-nodes
!*obc2uvatu*    REAL    Storage array for 2-D Orlanski and characteristic
!                       o.b.c. at U-nodes
!*obc2uvatu_old*REAL    Value of obc2uvatu at the first (outer) iteration
!                       with the implicit scheme
!*obc2uvatv*    REAL    Storage array for 2-D Orlanski and characteristic
!                       o.b.c. at V-nodes
!*obc2uvatv_old*REAL    Value of obc2uvatv at the first (outer) iteration
!                       with the implicit scheme
!*obc3uvatu*    REAL    Storage array for 3-D radiation o.b.c. at U-nodes
!*obc3uvatv*    REAL    Storage array for 3-D rdiation o.b.c. at V-nodes
!*return_time*  REAL    Vertical profile of return times as used by the
!                       Thatcher-Harleman open boundary condition            [s]
!*udatobu*      REAL    X-component of transport/depth-mean current/discharge
!                       at U-o.b. (residual or full)           [m^2/s,m/s,m^3/s]
!*udatobu_amp*  REAL    Amplitudes of X-transport at U-nodes             [m^2/s]
!*udatobu_pha*  REAL    Phases of X-transport at U-nodes                   [rad]
!*udatoby*      REAL    X-component of transport at Y-o.b. (residual or full)
!                                                                        [m^2/s]
!*udatoby_amp*  REAL    Amplitudes of X-transport at Y-nodes             [m^2/s]
!*udatoby_pha*  REAL    Phases of X-transport at Y-nodes                   [rad]
!*uveloby*      REAL    Value of uvel at external U-nodes, used as 3-D
!                       tangential boundary condition at Y-nodes           [m/s]
!*vdatobv*      REAL    Y-component of transport/depth-mean current/discharge
!                       at V-o.b. (residual or full)           [m^2/s,m/s,m^3/s]
!*vdatobv_amp*  REAL    Amplitudes of Y-transport at V-nodes             [m^2/s]
!*vdatobv_pha*  REAL    Phases of Y-transport at V-nodes                   [rad]
!*vdatobx*      REAL    Y-component of transport at X-o.b. (residual or full)
!
!*vdatobx_amp*  REAL    Amplitudes of Y-transport at X-nodes             [m^2/s]
!*vdatobx_pha*  REAL    Phases of Y-transport at X-nodes                   [rad]
!*vvelobx*    REAL      Value of vvel at external V-nodes, used as 3-D 
!                       tangential boundary condition at X-nodes           [m/s]
!*zdatobu*      REAL    Surface elevation/slope at U-open boundaries
!                       (residual or full)                                 [m,-]
!*zdatobu_amp*  REAL    Amplitudes of elevation at U-open boundaries         [m]
!*zdatobu_pha*  REAL    Phases of elevation at U-open boundaries           [rad]
!*zdatobv*      REAL    Surface elevation/slope at V-open boundaries
!                       (residual or full)                                 [m,-]
!*zdatobv_amp*   REAL    Amplitudes of elevation at V-open boundaries        [m]
!*zdatobv_pha*   REAL    Phases of elevation at V-open boundaries          [rad]
!*zeta_amp*     REAL    Amplitudes of surface elevation                      [m]
!*zeta_pha*     REAL    Phases of surface elevation                          [m]
!
!************************************************************************
!

END MODULE obconds
