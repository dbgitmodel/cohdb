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

MODULE morpharrays
!************************************************************************
!
! *morpharrays* Morphology arrays
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)morpharrays.f90  V2.10
!
! $Date: 2016-02-18 12:05:20 +0100 (Thu, 18 Feb 2016) $
!
! $Revision: 908 $
!
! Description - 
!
!************************************************************************
!
!
 
IMPLICIT NONE

INTEGER, ALLOCATABLE, DIMENSION(:) :: tidalstep
REAL, ALLOCATABLE, DIMENSION(:,:) :: active_layer_thickness, bed_update, &
                                   & bed_update_int, sed_avail_tot

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: bed_layer_thickness, &
     & bed_layer_thickness_old, bed_porosity, bed_update_dep, &
     & bed_update_dep_old1, bed_update_dep_old2, bed_update_ero, &
     & bed_update_ero_old1, bed_update_ero_old2, dar_morph, depth_guess, &
     & sed_avail, umvel_guess, vmvel_guess, zeta_guess

REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: bed_fraction_old

!
! Name                   Type  Purpose
!------------------------------------------------------------------------------
!
!*active_layer_thickness* REAL    The thickness of the layer used in the
!                                 vertical sorting model                     [m]
!*bed_fraction_old*       REAL    Bed fractions at previous time step        [-]
!*bed_layer_thickness*    REAL    Bed layer thicknesses                      [m]
!*bed_layer_thickness_old*REAL    Bed layer thicknesses at previous time
!                                 step                                       [m]
!*bed_porosity*           REAL    Fraction of sea water (porosity) within the
!                                 sea bed                                    [-]
!*bed_update*             REAL    The change in the bed level elevation between
!                                 two morphological time steps               [m]
!*bed_update_dep*         REAL    Change in bed level elevation between two
!                                 morphological time steps due to source
!                                 ("deposition") terms [m]                   [m]
!*bed_update_dep_old1*    REAL    bed_update_dep at previous time step       [m]
!*bed_update_dep_old2*    REAL    bed_update_dep at second previous time
!                                 step                                       [m]
!*bed_update_ero*         REAL    Change in bed level elevation between two
!                                 morphological time steps due to sink
!                                 ("erosion") terms                          [m]
!*bed_update_ero_old1*    REAL    bed_update_ero at previous time step       [m]
!*bed_update_ero_old2*    REAL    bed_update_ero at second previous time
!                                 step                                       [m]
!*bed_update_int*         REAL    The change in the bed level elevation
!                                 accumulated over time                      [m]
!*dar_morph*              REAL    Morphology source/sinks due to dredging and
!                                 relocation                             [m^3/s]
!*depth_guess*            REAL    Total water depth stored during the
!                                 hydrodynamical phase of a tidal averaging
!                                 cycle                                      [m]
!*sed_avail*              REAL    Amount of available sediment per fraction in
!                                 the bed layers(s)                          [m]
!*sed_avail_tot*          REAL    Total amount of available sediment in the
!                                 bed layers(s)                              [m]
!*tidalstep*              INTEGER Time (expressed as number of 2-D time steps)
!                                 between two subsequent updates for
!                                 storing/retrieving hydrodynamic data (used in
!                                 case of tidal acceleration)
!*umvelguess*             REAL    X-component of the depth-mean current stored
!                                 during the hydrodynamical phase of a tidal
!                                 averaging cycle                          [m/s]
!*vmvelguess*             REAL    Y-component of the depth-mean current stored
!                                 during the hydrodynamical phase of a tidal
!                                 averaging cycle                          [m/s]
!*zeta_guess*             REAL    Surface elevation stored during the
!                                 hydrodynamical phase of a tidal averaging
!                                 cycle                                      [m]
!
!************************************************************************

SAVE

END MODULE morpharrays
