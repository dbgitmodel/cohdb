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

MODULE sedids
!************************************************************************
!
! *sedids* Sediment variable key ids
!
! Author - Alexander Breugem and Patrick Luyten (IMDC)
!
! Version - @(COHERENS)sedids.f90  V2.11.1
!
! $Date: 2018-05-23 10:09:17 +0200 (Wed, 23 May 2018) $
!
! $Revision: 1138 $
!
! Description -
!
! Reference -
!
!************************************************************************
!

USE syspars

IMPLICIT NONE


!---start key id number
INTEGER, PARAMETER, PRIVATE :: n0 = MaxModArids

!
!1. Sediment particle attributes
!-------------------------------
!
!---particle attributes
INTEGER, PARAMETER :: &
& iarr_bstres_cr_cst = n0+1, iarr_dp = n0+2, iarr_massp = n0+3, &
& iarr_rhos = n0+4, iarr_volp = n0+5, iarr_ws_cst = n0+6

!---sediment concentration
INTEGER, PARAMETER :: &
& iarr_ceq = n0+7, iarr_cref = n0+8, iarr_ctot = n0+9, iarr_cvol = n0+10, &
& iarr_dar_sediment = n0+11

!---flocculation module
INTEGER, PARAMETER :: &
& iarr_cnump = n0+12, iarr_floc_dens = n0+13, iarr_floc_dia = n0+14, &
& iarr_floc_F = n0+15, iarr_floc_nc = n0+16, iarr_floc_P = n0+17, &
& iarr_floc_T = n0+18, iarr_floc_vol = n0+19

!---sediment loads
INTEGER, PARAMETER :: &
& iarr_qbedatc = n0+20, iarr_qbedatu = n0+21, iarr_qbedatu_int = n0+22, &
& iarr_qbedatv = n0+23, iarr_qbedatv_int = n0+24, iarr_qsusatc = n0+25, &
& iarr_qsusatu = n0+26, iarr_qsusatu_int = n0+27, iarr_qsusatv = n0+28, &
& iarr_qsusatv_int = n0+29, iarr_qtotatc = n0+30, iarr_qtotatu = n0+31, &
& iarr_qtotatu_int = n0+32, iarr_qtotatv = n0+33, iarr_qtotatv_int = n0+34

!---bottom fluxes
INTEGER, PARAMETER :: &
& iarr_bed_fraction = n0+35, iarr_bottom_sed_dep = n0+36, &
& iarr_bottom_sed_ero = n0+37, iarr_bottom_sed_flux = n0+38, &
& iarr_height_c = n0+39, iarr_t_equil = n0+40

!---bottom stress related arrays
INTEGER, PARAMETER :: &
& iarr_bdragcoefatc_sed = n0+41, iarr_bfricvel_sed = n0+42, &
& iarr_bfricvel_mean_sed = n0+43, iarr_bfricvel_wav_sed = n0+44, &
& iarr_bstres_cr = n0+45, iarr_bstresatc_sed = n0+46, &
& iarr_bstresatc_mean_sed = n0+47, iarr_bstresatc_wav_sed = n0+48, &
& iarr_bstresatu_sed = n0+49, iarr_bstresatv_sed = n0+50, &
& iarr_d50_bed = n0+51, iarr_fwave_sed = n0+52, iarr_ubstresatc_sed = n0+53, &
& iarr_ubstresatu_sed = n0+54, iarr_vbstresatc_sed = n0+55, &
& iarr_vbstresatv_sed = n0+56, iarr_wavethickatc_sed = n0+57, &
& iarr_zaroughatc_sed = n0+58, iarr_zroughatc_sed = n0+59, &
& iarr_zroughatu_sed = n0+60, iarr_zroughatv_sed = n0+61

!---diffusion coefficients
INTEGER, PARAMETER :: &
& iarr_beta_sed = n0+62, iarr_vdiffcoef_sed = n0+63

!---bed slopes
INTEGER, PARAMETER :: &
& iarr_bed_slope_x = n0+64, iarr_bed_slope_x_atu = n0+65, &
& iarr_bed_slope_x_atv = n0+66, iarr_bed_slope_y = n0+67, &
& iarr_bed_slope_y_atu = n0+68, iarr_bed_slope_y_atv = n0+69

!---density
INTEGER, PARAMETER :: &
& iarr_beta_state_sed = n0+70, iarr_densw = n0+71

!---fall velocity
INTEGER, PARAMETER :: &
& iarr_wfall = n0+72

!---miscellaneous
INTEGER, PARAMETER :: &
& iarr_tidalstep = n0+73

!---open boundary conditions
INTEGER, PARAMETER :: &
& iarr_obccnumpatu = n0+74, iarr_obccnumpatv = n0+75, iarr_obcsedatu = n0+76, &
& iarr_obcsedatv = n0+77

!
!2. Morphology
!-------------
!

INTEGER, PARAMETER :: &
& iarr_active_layer_thickness = n0+78, iarr_bed_layer_thickness = n0+79, &
& iarr_bed_porosity = n0+80, iarr_bed_update = n0+81, &
& iarr_bed_update_dep = n0+82, iarr_bed_update_dep_old1 = n0+83, &
& iarr_bed_update_dep_old2 = n0+84, iarr_bed_update_ero = n0+85, &
& iarr_bed_update_ero_old1 = n0+86, iarr_bed_update_ero_old2 = n0+87, &
& iarr_bed_update_int = n0+88, iarr_dar_morph = n0+89, &
& iarr_depth_guess = n0+90, iarr_sed_avail = n0+91, &
& iarr_sed_avail_tot = n0+92, iarr_sediment_balance = n0+93, &
& iarr_sediment_obc = n0+94, iarr_sediment_vol = n0+95, &
& iarr_umvel_guess = n0+96, iarr_vmvel_guess = n0+97, iarr_zeta_guess = n0+98

!
!3. Dredging and relocation
!--------------------------
!
!---campaigns
INTEGER, PARAMETER :: &
& iarr_camp_campdepth = n0+99, iarr_camp_campfractions = n0+100, &
& iarr_camp_camplevel = n0+101, iarr_camp_CampaignStartDateTime = n0+102, &
& iarr_camp_campvolume = n0+103, iarr_camp_CoupleTime = n0+104, &
& iarr_camp_dredge = n0+105, iarr_camp_nr_dredgingsite = n0+106, &
& iarr_camp_nr_relocationsite = n0+107, iarr_camp_nr_ship = n0+108, &
& iarr_camp_nr_zones_coupled = n0+109, iarr_camp_relocate = n0+110, &
& iarr_camp_relocationsites = n0+111, iarr_camp_t_lag = n0+112, &
& iarr_camp_t_spread = n0+113

!---dredging sites
INTEGER, PARAMETER :: &
& iarr_dredge_CoupleTime = n0+114, iarr_dredge_critical_volume = n0+115, &
& iarr_dredge_dredging_thickness = n0+116, &
& iarr_dredge_monitoring_frequency = n0+117, iarr_dredge_nr_coordpol = n0+118, &
& iarr_dredge_nr_coupled = n0+119, iarr_dredge_nr_ship = n0+120, &
& iarr_dredge_nr_zones_coupled = n0+121, iarr_dredge_over_depth = n0+122, &
& iarr_dredge_relocationsites = n0+123, iarr_dredge_target_depth = n0+124, &
& iarr_dredge_t_lag = n0+125, iarr_dredge_t_spread = n0+126, &
& iarr_dredge_Xpol = n0+127, iarr_dredge_Ypol = n0+128

!---relocation sites
INTEGER, PARAMETER :: &
& iarr_reloc_max_volume = n0+129, iarr_reloc_min_depth = n0+130, &
& iarr_reloc_nr_coordpol = n0+131, iarr_reloc_Xpol = n0+132, &
& iarr_reloc_Ypol = n0+133

!---ships
INTEGER, PARAMETER :: &
& iarr_ship_D_dredge_max = n0+134, iarr_ship_D_dredge_min = n0+135, &
& iarr_ship_D_relocate_max = n0+136, iarr_ship_D_relocate_min = n0+137, &
& iarr_ship_H_dredge_max = n0+138, iarr_ship_H_dredge_min = n0+139, &
& iarr_ship_H_relocate_max = n0+140, iarr_ship_H_relocate_min = n0+141, &
& iarr_ship_p_overflow = n0+142, iarr_ship_p_passive = n0+143, &
& iarr_ship_p_relocation = n0+144, iarr_ship_p_suctionhead = n0+145, &
& iarr_ship_R_circle = n0+146, iarr_ship_t_fill = n0+147, &
& iarr_ship_t_idle = n0+148, iarr_ship_t_relocate = n0+149, &
& iarr_ship_t_sail_dredge = n0+150, iarr_ship_t_sail_relocate = n0+151, &
& iarr_ship_t_wait = n0+152, iarr_ship_U_ship = n0+153, &
& iarr_ship_V_ship = n0+154

!--- sub-zone dredging
INTEGER, PARAMETER :: &
& iarr_subdred_nr_coordpol = n0+155, iarr_subdred_nr_subzones = n0+156, &
& iarr_subdred_Xpol = n0+157, iarr_subdred_Ypol = n0+158

!---sub-zone relocation
INTEGER, PARAMETER :: &
& iarr_subrel_nr_coordpol = n0+159, iarr_subrel_nr_subzones = n0+160, &
& iarr_subrel_Xpol = n0+161, iarr_subrel_Ypol = n0+162

!---other arrays
INTEGER, PARAMETER :: &
& iarr_dredging_depth = n0+163, iarr_D_dredge_max = n0+164, &
& iarr_D_dredge_min = n0+165, iarr_D_relocate_max = n0+166, &
& iarr_D_relocate_min = n0+167, iarr_fraction_ship = n0+168, &
& iarr_H_dredge_max = n0+169, iarr_H_dredge_min = n0+170, &
& iarr_H_relocate_max = n0+171, iarr_H_relocate_min = n0+172, &
& iarr_itrack = n0+173, iarr_jtrack = n0+174, iarr_mask_dredging = n0+175, &
& iarr_mask_relocation = n0+176, iarr_mask_subzones_dredging = n0+177, &
& iarr_mask_subzones_relocation = n0+178, iarr_nocouplingtimes = n0+179, &
& iarr_notrackpoints = n0+180, iarr_p_overflow = n0+181, &
& iarr_p_passive = n0+182, iarr_p_relocation = n0+183, &
& iarr_p_suctionhead = n0+184, iarr_track_points_counter = n0+185, &
& iarr_volume_ship = n0+186, iarr_V_dredge = n0+187, &
& iarr_V_dredge_max = n0+188, iarr_V_dredge_rest = n0+189, &
& iarr_V_dredge_total = n0+190, iarr_V_lastcycle = n0+191, &
& iarr_V_relocate_extra = n0+192, iarr_V_relocate_total = n0+193, &
& iarr_V_ship_rest = n0+194, iarr_xtrack = n0+195, iarr_ytrack = n0+196


END MODULE sedids
