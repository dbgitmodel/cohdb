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

MODULE structures
!
!*****************************************************************
!
! *structures* Parameters and arrays for structures and discharges model unit
!
! Author - ANTEA
!
! Version - @(COHERENS)structures.f90  V2.7
!
! $Date: 2017-02-13 14:59:12 +0100 (Mon, 13 Feb 2017) $
!
! $Revision: 1004 $
!
! Description -
!
!*****************************************************************
!
USE syspars

IMPLICIT NONE

!---dry cells
INTEGER :: numdry
INTEGER, ALLOCATABLE, DIMENSION(:) :: idry, jdry

!---thin dams
INTEGER :: numthinu, numthinuloc, numthinv, numthinvloc
INTEGER, ALLOCATABLE, DIMENSION(:) :: ithinu, ithinv, jthinu, jthinv
INTEGER, ALLOCATABLE, DIMENSION(:) :: ithinuloc, ithinvloc, jthinuloc, jthinvloc

!---weirs/barriers/orifices
INTEGER :: numwbaru, numwbaruloc, numwbarv, numwbarvloc
INTEGER, ALLOCATABLE, DIMENSION(:) :: iwbaru, iwbarv, jwbaru, jwbarv
INTEGER, ALLOCATABLE, DIMENSION(:) :: indexwbaru, indexwbarv, iwbaruloc, &
                                    & iwbarvloc, jwbaruloc, jwbarvloc
REAL :: wbarrlxu = 1.0, wbarrlxv = 1.0
REAL, ALLOCATABLE, DIMENSION(:) :: oricoefu, oricoefv, oriheightu, oriheightv, &
                                 & orisillu, orisillv
REAL, ALLOCATABLE, DIMENSION(:) :: wbarcoefu, wbarcoefv, wbarcrestu, &
                                 & wbarcrestv, wbarmodlu, wbarmodlv
REAL, ALLOCATABLE, DIMENSION(:) :: wbarelossu, wbarelossv

!---discharges
INTEGER :: numdis, numdisloc, numdisloc_ext
LOGICAL, ALLOCATABLE, DIMENSION(:) :: disflag
INTEGER, ALLOCATABLE, DIMENSION(:) :: idis, idisloc, indexdisloc, jdis, &
                                    & jdisloc, kdis, kdistype
REAL, ALLOCATABLE, DIMENSION(:) :: disarea, disdir, disspeed, disvol, &
                                 & xdiscoord, ydiscoord, zdiscoord

!---parallel processing
INTEGER, ALLOCATABLE, DIMENSION(:) :: nodisprocs, nowbaruprocs, nowbarvprocs
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: indexdisprocs, indexwbaruprocs, &
                                      & indexwbarvprocs 

SAVE

!
! Name            Type    Purpose
!------------------------------------------------------------------------------
!---dry cells
!*idry*           INTEGER X-indices of dry cells
!*jdry*           INTEGER Y-indices of dry cells
!*numdry*         INTEGER (Global) number of dry cells
!---thin dams
!*ithinu*         INTEGER Global X-indices of thin dams at U-nodes
!*ithinuloc*      INTEGER Local X-indices of thin dams at U-nodes
!*ithinv*         INTEGER Global X-indices of thin dams at V-nodes
!*ithinvloc*      INTEGER Local X-indices of thin dams at V-nodes
!*jthinu*         INTEGER Global Y-indices of thin dams at U-nodes
!*jthinuloc*      INTEGER Local Y-indices of thin dams at U-nodes
!*jthinv*         INTEGER Global Y-indices of thin dams at V-nodes
!*jthinvloc*      INTEGER Local Y-indices of thin dams at V-nodes
!*numthinu*       INTEGER Global number of thin dams at U-nodes
!*numthinuloc*    INTEGER Local number of thin dams at U-nodes
!*numthinv*       INTEGER Global number of thin dams at V-nodes
!*numthinvloc*    INTEGER Local number of thin dams at V-nodes
!---weirs/barriers/orifices
!*indexwbaru*     INTEGER Global indices of the local U-node weirs/barriers
!*indexwbaruprocs*INTEGER Global indices per process of the local U-node
!                         weir/barrier points
!*indexwbarv*     INTEGER Global indices of the local V-node weirs/barriers
!*indexwbarvprocs*INTEGER Global indices per process of the local V-node
!                         weir/barrier points
!*iwbaru*         INTEGER Global X-indices of weirs/barriers at U-nodes
!*iwbaruloc*      INTEGER Local X-indices of weirs/barriers at U-nodes
!*iwbarv*         INTEGER Global X-indices of weirs/barriers at V-nodes
!*iwbarvloc*      INTEGER Local X-indices of weirs/barriers at V-nodes
!*jwbaru*         INTEGER Global Y-indices of weirs/barriers at U-nodes
!*jwbaruloc*      INTEGER Local Y-indices of weirs/barriers at U-nodes
!*jwbarv*         INTEGER Global Y-indices of weirs/barriers at V-nodes
!*jwbarvloc*      INTEGER Local Y-indices of weirs/barriers at V-nodes
!*nowbaruprocs*   INTEGER Array with number of U-node weirs/barriers per process
!*nowbarvprocs*   INTEGER Array with number of V-node weirs/barriers per process
!*numwbaru*       INTEGER Global number of weirs/barriers at U-nodes
!*numwbaruloc*    INTEGER Local number of weirs/barriers at U-nodes
!*numwbarv*       INTEGER Global number of weirs/barriers at V-nodes
!*numwbarvloc*    INTEGER Local number of weirs/barriers at V-nodes
!*oricoefu*       REAL    Discharge coefficient for orifices at U-nodes
!                                                                   [m^1/2 s^-1]
!*oricoefv*       REAL    Discharge coefficient for orifices at V-nodes
!                                                                   [m^1/2 s^-1]
!*oriheightu*     REAL    Orifice width at U-nodes                           [m]
!*oriheightv*     REAL    Orifice width at V-nodes                           [m]
!*orisillu*       REAL    Orifice height at U-nodes                          [m]
!*orisillv*       REAL    Orifice height at V-nodes                          [m]
!*wbarcoefu*      REAL    Discharge coefficient for weirs/barriers at U-nodes
!                                                                   [m^1/2 s^-1]
!*wbarcoefv*      REAL    Discharge coefficient for weirs/barriers at V-nodes
!                                                                   [m^1/2 s^-1]
!*wbarcrestu*     REAL    Height of weir crest at U-nodes                    [m]
!*wbarcrestv*     REAL    Height of weir crest at V-nodes                    [m]
!*wbarelossu*     REAL    Energy loss sink term at U-node weirs/barriers   [1/s]
!*wbarelossv*     REAL    Energy loss sink term at V-node weirs/barriers   [1/s]
!*wbarmodlu*      REAL    Modular limit at U-node weirs/barriers
!*wbarmodlv*      REAL    Modular limit at V-node weirs/barriers
!*wbarrlxu*       REAL    Time relaxation coefficient at U-node weirs/barriers
!*wbarrlxv*       REAL    Time relaxation coefficient at V-node weirs/barriers
!---discharges
!*disarea*        REAL    Area over wwhich the discharge takes place       [m^2]
!*disdir*         REAL    Discharge direction                              [rad]
!*disflag*        LOGICAL .TRUE./.FALSE. if a discharge locations is located
!                         on a wet/dry cell
!*disspeed*       REAL    Discharge speed defined as the volume discharge
!                         rate divided by the cell area                    [m^2]
!*disvol*         REAL    Volume discharge rate [m^3/s]
!*idis*           INTEGER Global X-indices of discharge locations
!*idisloc*        INTEGER Local X-indices of discharge locations
!*indexdisloc*    INTEGER Global indices of the local discharge points
!*indexdisprocs*  INTEGER Global indices per process of the local discharge
!                         points
!*jdis*           INTEGER Global Y-indices of discharge locations
!*jdisloc*        INTEGER Local Y-indices of discharge locations
!*kdis*           INTEGER Vertical grid indices of the discharges locations.
!                         If zero, discharges are taken as homogeneous over
!                         the vertical
!*kdistype*       INTEGER Selects type of vertical location of the discharge
!                  = 0 =>  Uniformly distributed over the vertical
!                  = 1 =>  At a givven cell location
!                  = 2 =>  At a fixed distance from the sea bed
!                  = 3 =>  At a fixed distance from the sea surface
!*nodisprocs*     INTEGER Array with number of discharge locations per process 
!*numdis*         INTEGER Global number of discharge locations
!*numdisloc*      INTEGER Local number of discharge locations
!*numdisloc_ext*  INTEGER Local number of discharge locations including points
!                         within the first column of the western and first row
!                         of the southern halo
!*xdiscoord*      REAL    X-coordinates of discharge locations
!*ydiscoord*      REAL    Y-coordinates of discharge locations
!*zdiscoord*      REAL    Vertical coordinates (distance from sea bed or
!                         surface) of discharge locations
!
!************************************************************************
!


END MODULE structures
