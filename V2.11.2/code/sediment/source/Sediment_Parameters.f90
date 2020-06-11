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
! *Sediment_Parameters* Routines for reading and writing a sediment CIF
!
! Author -  Alexander Breugem, Kevin Delecluyse and Patrick Luyten
!
! Version - @(COHERENS)Sediment_Parameters.f90  V2.11.1
!
! $Date: 2018-05-23 10:09:17 +0200 (Wed, 23 May 2018) $
!
! $Revision: 1138 $
!
! Description -
!
! Routines - assign_cif_vars_dar, assign_cif_vars_morph, assign_cif_vars_sed,
!            write_cif_vars_dar, write_cif_vars_morph, write_cif_vars_sed
!
!************************************************************************
!

!============================================================================

SUBROUTINE assign_cif_vars_dar(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_dar* convert the string data from an input line in the CIF
!                       block with dredging/relocation setup parameters to the
!                       appropriate numeric or non-numeric format 
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Sediment_Parameters.f90  V2.8
!
! Description - 
!
! Reference -
!
! Calling program - read_cif_params
!
! External calls -
!
! Module calls - conv_from_chars
!
!************************************************************************
!
USE iopars
USE darpars
USE darswitches
USE syspars
USE cif_routines, ONLY: conv_from_chars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lenname), INTENT(IN), OPTIONAL :: cname
CHARACTER (LEN=lencifvar), INTENT(IN), DIMENSION(MaxCIFvars) :: cvals
INTEGER, INTENT(IN) :: numvars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cname*     CHAR    Name of the variable(s)
!*cvals*     CHAR    Data values on input data line
!*numvars*   INTEGER Number of data variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: lvar


lvar = 1

SELECT CASE (TRIM(cname))

!
!1. Switches
!-----------
!

CASE('IOPT_DAR_COUPLING_SPACE')
   CALL conv_from_chars(cvals(1),iopt_dar_coupling_space,lvar)
CASE('IOPT_DAR_COUPLING_TIME')
   CALL conv_from_chars(cvals(1),iopt_dar_coupling_time,lvar)
CASE('IOPT_DAR_DREDGING_CRITERION')
   CALL conv_from_chars(cvals(1),iopt_dar_dredging_criterion,lvar)
CASE('IOPT_DAR_DURATION')
   CALL conv_from_chars(cvals(1),iopt_dar_duration,lvar)
CASE('IOPT_DAR_EFFECT')
   CALL conv_from_chars(cvals(1),iopt_dar_effect,lvar)
CASE('IOPT_DAR_RELOCATION_CRITERION')
   CALL conv_from_chars(cvals(1),iopt_dar_relocation_criterion,lvar)
CASE('IOPT_DAR_SCENARIO')
   CALL conv_from_chars(cvals(1),iopt_dar_scenario,lvar)
CASE('IOPT_DAR_TIME_DREDGE')
   CALL conv_from_chars(cvals(1),iopt_dar_time_dredge,lvar)
CASE('IOPT_DAR_TIME_RELOCATE')
   CALL conv_from_chars(cvals(1),iopt_dar_time_relocate,lvar)
CASE('IOPT_DAR_USER_VOLUME')
   CALL conv_from_chars(cvals(1),iopt_dar_user_volume,lvar)

!2.2 Model parameters
!--------------------
!

CASE('MAXTRACKPOINTS')
   CALL conv_from_chars(cvals(1),maxtrackpoints,lvar)
CASE('NOCAMPAIGNS')
   CALL conv_from_chars(cvals(1),nocampaigns,lvar)
CASE('NOCIRCLEPOINTS')
   CALL conv_from_chars(cvals(1),nocirclepoints,lvar)
CASE('NODREDGINGSITES')
CALL conv_from_chars(cvals(1),nodredgingsites,lvar)
CASE('NORELOCATIONSITES')
CALL conv_from_chars(cvals(1),norelocationsites,lvar)

END SELECT


RETURN

END SUBROUTINE assign_cif_vars_dar

!========================================================================

SUBROUTINE assign_cif_vars_morph(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_morph* convert the string data from an input line in the CIF
!                         block with parameters for the morphological module to
!                         the appropriate numeric or non-numeric format  
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Sediment_Parameters.f90  V2.10
!
! Description - 
!
! Reference -
!
! Calling program - read_cif_mod_params
!
! External calls -
!
! Module calls - conv_from_chars
!
!************************************************************************
!

USE iopars
USE switches
USE syspars
USE morphpars
USE morphswitches
USE cif_routines, ONLY: conv_from_chars

IMPLICIT NONE
!
!*Arguments
!
CHARACTER (LEN=lenname), INTENT(IN), OPTIONAL :: cname
CHARACTER (LEN=lencifvar), INTENT(IN), DIMENSION(MaxCIFvars) :: cvals
INTEGER, INTENT(IN) :: numvars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cname*     CHAR    Name of the variable(s)
!*cvals*     CHAR    Data values on input data line
!*numvars*   INTEGER Number of data variables

!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: lvar


lvar = 1
SELECT CASE (TRIM(cname))

!
!1. Switches
!------------
!

CASE('IOPT_MORPH_ACTIVE_LAYER')
   CALL conv_from_chars(cvals(1),iopt_morph_active_layer,lvar)
CASE('IOPT_MORPH_AVALANCHING')
   CALL conv_from_chars(cvals(1),iopt_morph_avalanching,lvar)
   CASE('IOPT_MORPH_CORR')
   CALL conv_from_chars(cvals(1),iopt_morph_corr,lvar)
CASE('IOPT_MORPH_FIXED_LAYER')
   CALL conv_from_chars(cvals(1),iopt_morph_fixed_layer,lvar)
CASE('IOPT_MORPH_HDIFF')
   CALL conv_from_chars(cvals(1),iopt_morph_hdiff,lvar)
CASE('IOPT_MORPH_TIDAL_SCHEME')
   CALL conv_from_chars(cvals(1),iopt_morph_tidal_scheme,lvar)
CASE('IOPT_MORPH_TIME_INT')
   CALL conv_from_chars(cvals(1),iopt_morph_time_int,lvar)
CASE('IOPT_MORPH_VERT_FLUXES')
   CALL conv_from_chars(cvals(1),iopt_morph_vert_fluxes,lvar)

!
!2. Parameters
!-------------
!

CASE('ACCELSTEP')
   CALL conv_from_chars(cvals(1),accelstep,lvar)
CASE('ALPHA_JAMESON')
   CALL conv_from_chars(cvals(1),alpha_jameson,lvar)
CASE('AVERAGE_BEDFORM_HEIGHT')
   CALL conv_from_chars(cvals(1),average_bedform_height,lvar)
CASE('AVERAGE_BEDFORM_LENGTH')
   CALL conv_from_chars(cvals(1),average_bedform_length,lvar)
CASE('BED_POROSITY_CST')
   CALL conv_from_chars(cvals(1),bed_porosity_cst,lvar)
CASE('DUNE_CELERITY'); CALL conv_from_chars(cvals(1),dune_celerity,lvar)
CASE('DYN_ANGLE_DRY'); CALL conv_from_chars(cvals(1),dyn_angle_dry,lvar)
CASE('DYN_ANGLE_WET'); CALL conv_from_chars(cvals(1),dyn_angle_wet,lvar)
CASE('ICAVAL'); CALL conv_from_chars(cvals(1),icaval,lvar)
CASE('ICMORPH'); CALL conv_from_chars(cvals(1),icmorph,lvar)
CASE('ICSEDBAL'); CALL conv_from_chars(cvals(1),icsedbal,lvar)
CASE('K1_HARRIS'); CALL conv_from_chars(cvals(1),k1_harris,lvar)
CASE('K2_HARRIS'); CALL conv_from_chars(cvals(1),k2_harris,lvar)
CASE('K2_JAMESON'); CALL conv_from_chars(cvals(1),k2_jameson,lvar)
CASE('K4_JAMESON'); CALL conv_from_chars(cvals(1),k4_jameson,lvar)
CASE('MAX_ITERATIONS_AVAL')
   CALL conv_from_chars(cvals(1),max_iterations_aval,lvar)
CASE('MORPH_FACTOR'); CALL conv_from_chars(cvals(1),morph_factor,lvar)
CASE('MORPH_STEPS')
   CALL conv_from_chars(cvals(1),morph_steps,lvar)
CASE ('NSTEP_HYDRO'); CALL conv_from_chars(cvals(1),nstep_hydro,lvar)
CASE('NUMBER_TIDAL_STEPS')
   CALL conv_from_chars(cvals(1),number_tidal_steps,lvar)
CASE('SIMILARITY_RANGE'); CALL conv_from_chars(cvals(1),similarity_range,lvar)
CASE('STAT_ANGLE_DRY'); CALL conv_from_chars(cvals(1),stat_angle_dry,lvar)
CASE('STAT_ANGLE_WET'); CALL conv_from_chars(cvals(1),stat_angle_wet,lvar)
CASE('TROUGH_PROBABILITY')
   CALL conv_from_chars(cvals(1),trough_probability,lvar)


END SELECT


RETURN

END SUBROUTINE assign_cif_vars_morph

!============================================================================

SUBROUTINE assign_cif_vars_sed(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_sed* convert the string data from an input line in the CIF
!                       block with parameters for the sediment module to
!                       the appropriate numeric or non-numeric format    
!
! Author - Alexander Breugem (IMDC)
!
! Version - @(COHERENS)Sediment_Parameters.f90  V2.11.1
!
! Description - 
!
! Reference -
!
! Calling program - read_cif_params
!
! External calls -
!
! Module calls - conv_from_chars
!
!************************************************************************
!
USE iopars
USE sedpars
USE sedswitches
USE syspars
USE cif_routines, ONLY: conv_from_chars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lenname), INTENT(IN), OPTIONAL :: cname
CHARACTER (LEN=lencifvar), INTENT(IN), DIMENSION(MaxCIFvars) :: cvals
INTEGER, INTENT(IN) :: numvars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cname*     CHAR    Name of the variable(s)
!*cvals*     CHAR    Data values on input data line
!*numvars*   INTEGER Number of data variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: lvar

lvar = 1
SELECT CASE (TRIM(cname))

!
!1. Switches
!-----------
!

CASE('IOPT_SED_BBC'); CALL conv_from_chars(cvals(1),iopt_sed_bbc,lvar)
CASE('IOPT_SED_BBC_EQ')
   CALL conv_from_chars(cvals(1),iopt_sed_bbc_eq,lvar)
CASE('IOPT_SED_BBC_REF')
   CALL conv_from_chars(cvals(1),iopt_sed_bbc_ref,lvar)
CASE('IOPT_SED_BEDEQ')
   CALL conv_from_chars(cvals(1),iopt_sed_bedeq,lvar)
CASE('IOPT_SED_BETA'); CALL conv_from_chars(cvals(1),iopt_sed_beta,lvar)
CASE('IOPT_SED_BSTRES'); CALL conv_from_chars(cvals(1),iopt_sed_bstres,lvar)
CASE('IOPT_SED_BSTRES_CR')
   CALL conv_from_chars(cvals(1),iopt_sed_bstres_cr,lvar)
CASE('IOPT_SED_DENS_GRAD')
   CALL conv_from_chars(cvals(1),iopt_sed_dens_grad,lvar)
CASE('IOPT_SED_FILTER'); CALL conv_from_chars(cvals(1),iopt_sed_filter,lvar)
CASE('IOPT_SED_HIDING'); CALL conv_from_chars(cvals(1),iopt_sed_hiding,lvar)
CASE('IOPT_SED_WS_FLOC'); CALL conv_from_chars(cvals(1),iopt_sed_ws_floc,lvar)
CASE('IOPT_SED_WS_HINDSET')
   CALL conv_from_chars(cvals(1),iopt_sed_ws_hindset,lvar)
CASE('IOPT_SED_MEDIAN')
   CALL conv_from_chars(cvals(1),iopt_sed_median,lvar)
CASE('IOPT_SED_MODE'); CALL conv_from_chars(cvals(1),iopt_sed_mode,lvar)
CASE('IOPT_SED_NODIM'); CALL conv_from_chars(cvals(1),iopt_sed_nodim,lvar)
CASE('IOPT_SED_OBC_FLUX'); CALL conv_from_chars(cvals(1),iopt_sed_obc_flux,lvar)
CASE('IOPT_SED_ROUGH'); CALL conv_from_chars(cvals(1),iopt_sed_rough,lvar)
CASE('IOPT_SED_SLOPE'); CALL conv_from_chars(cvals(1),iopt_sed_slope,lvar)
CASE('IOPT_SED_SOURCE_IMPL')
   CALL conv_from_chars(cvals(1),iopt_sed_source_impl,lvar)
CASE('IOPT_SED_TOTEQ'); CALL conv_from_chars(cvals(1),iopt_sed_toteq,lvar)
CASE('IOPT_SED_TYPE'); CALL conv_from_chars(cvals(1),iopt_sed_type,lvar)
CASE('IOPT_SED_VADV'); CALL conv_from_chars(cvals(1),iopt_sed_vadv,lvar)
   CALL conv_from_chars(cvals(1),iopt_sed_vadv,lvar)
CASE('IOPT_SED_WAVE_DIFF')
   CALL conv_from_chars(cvals(1),iopt_sed_wave_diff,lvar)
CASE('IOPT_SED_WS'); CALL conv_from_chars(cvals(1),iopt_sed_ws,lvar)
CASE('IOPT_SED_WS_LIM'); CALL conv_from_chars(cvals(1),iopt_sed_ws_lim,lvar)

!
!2. Parameters
!-------------
!
!2.1 Sediment model
!------------------
!
   
CASE('ICSED'); CALL conv_from_chars(cvals(1),icsed,lvar)
CASE('MAXITBARTNICKI')
   CALL conv_from_chars(cvals(1),maxitbartnicki,lvar)
CASE('NB'); CALL conv_from_chars(cvals(1),nb,lvar)
CASE('NF'); CALL conv_from_chars(cvals(1),nf,lvar)
CASE('NRQUAD_SED'); CALL conv_from_chars(cvals(1),nrquad_sed,lvar)
CASE('NRQUAD_WAV'); CALL conv_from_chars(cvals(1),nrquad_wav,lvar)
CASE('ALPHA_VR'); CALL conv_from_chars(cvals(1),alpha_VR,lvar)
CASE('A_LEUSSEN'); CALL conv_from_chars(cvals(1),a_leussen,lvar)
CASE('BETA_SED_CST'); CALL conv_from_chars(cvals(1),beta_sed_cst,lvar)
CASE('B_LEUSSEN'); CALL conv_from_chars(cvals(1),b_leussen,lvar)
CASE('CGEL'); CALL conv_from_chars(cvals(1),cgel,lvar)
CASE('CMAX'); CALL conv_from_chars(cvals(1),cmax,lvar)
CASE('COEF_BED_GRAD'); CALL conv_from_chars(cvals(1),coef_bed_grad,lvar)
CASE('FLOC_VR_MAX'); CALL conv_from_chars(cvals(1),floc_VR_max,lvar)
CASE('FLOC_VR_MIN'); CALL conv_from_chars(cvals(1),floc_VR_min,lvar)
CASE('HEIGHT_C_CST'); CALL conv_from_chars(cvals(1),height_c_cst,lvar)
CASE('N_RICHZAKI'); CALL conv_from_chars(cvals(1),n_RichZaki,lvar)
CASE('PARTH_COEF'); CALL conv_from_chars(cvals(1),parth_coef,lvar)
CASE('PARTH_EXP'); CALL conv_from_chars(cvals(1),parth_exp,lvar)
CASE('ZROUGH_GRAIN'); CALL conv_from_chars(cvals(1),zrough_grain,lvar)
CASE('ZROUGH_SED_CST'); CALL conv_from_chars(cvals(1),zrough_sed_cst,lvar)
CASE('Z0_COEF'); CALL conv_from_chars(cvals(1),z0_coef,lvar)
CASE('WSLIMFAC'); CALL conv_from_chars(cvals(1),wslimfac,lvar)
CASE('WU_EXP'); CALL conv_from_chars(cvals(1),wu_exp,lvar)

!   
!2.2 Flocculation model
!----------------------
!
   
CASE('AGG_ALPHA'); CALL conv_from_chars(cvals(1),agg_alpha,lvar)
CASE('BRK_ES'); CALL conv_from_chars(cvals(1),brk_es,lvar)
CASE('BRK_FFY'); CALL conv_from_chars(cvals(1),brk_ffy,lvar)
CASE('BRK_FRAC'); CALL conv_from_chars(cvals(1),brk_frac,lvar)
CASE('BRK_P'); CALL conv_from_chars(cvals(1),brk_p,lvar)
CASE('BRK_Q'); CALL conv_from_chars(cvals(1),brk_q,lvar)
CASE('FLOC_NCMAX'); CALL conv_from_chars(cvals(1),floc_ncmax,lvar)
CASE('FLOC_NCMIN'); CALL conv_from_chars(cvals(1),floc_ncmin,lvar)
CASE('NFRDIM'); CALL conv_from_chars(cvals(1),nfrdim,lvar)
   
END SELECT


RETURN

END SUBROUTINE assign_cif_vars_sed

!========================================================================

SUBROUTINE write_cif_vars_dar
!************************************************************************
!
! *write_cif_vars_dar* write the CIF block with parameters for the dredging and
!                      relocation model
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Sediment_Parameters.f90  V2.11.1
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - conv_to_chars, write_cif_line
!
!************************************************************************
!
USE darpars
USE darswitches
USE iopars
USE paralpars
USE cif_routines, ONLY: conv_to_chars, write_cif_line
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: iunit


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_cif_vars_dar' 
CALL log_timer_in()

!
!1. Write block header
!--------------------
!

iunit = ciffile%iunit
WRITE (iunit,'(A)') cifend//TRIM(cif_block_names(icif_dar))

!
!2. Switches
!-----------
!

CALL conv_to_chars(cvals(1),iopt_dar_coupling_space)
CALL write_cif_line(cvals(1:1),'IOPT_DAR_COUPLING_SPACE')

CALL conv_to_chars(cvals(1),iopt_dar_coupling_time)
CALL write_cif_line(cvals(1:1),'IOPT_DAR_COUPLING_TIME')

CALL conv_to_chars(cvals(1),iopt_dar_dredging_criterion)
CALL write_cif_line(cvals(1:1),'IOPT_DAR_DREDGING_CRITERION')

CALL conv_to_chars(cvals(1),iopt_dar_duration)
CALL write_cif_line(cvals(1:1),'IOPT_DAR_DURATION')

CALL conv_to_chars(cvals(1),iopt_dar_effect)
CALL write_cif_line(cvals(1:1),'IOPT_DAR_EFFECT')

CALL conv_to_chars(cvals(1),iopt_dar_relocation_criterion)
CALL write_cif_line(cvals(1:1),'IOPT_DAR_RELOCATION_CRITERION')

CALL conv_to_chars(cvals(1),iopt_dar_scenario)
CALL write_cif_line(cvals(1:1),'IOPT_DAR_SCENARIO')

CALL conv_to_chars(cvals(1),iopt_dar_time_dredge)
CALL write_cif_line(cvals(1:1),'IOPT_DAR_TIME_DREDGE')

CALL conv_to_chars(cvals(1),iopt_dar_time_relocate)
CALL write_cif_line(cvals(1:1),'IOPT_DAR_TIME_RELOCATE')

CALL conv_to_chars(cvals(1),iopt_dar_user_volume)
CALL write_cif_line(cvals(1:1),'IOPT_DAR_USER_VOLUME')

!
!3. Parameters
!-------------
!

CALL conv_to_chars(cvals(1),maxtrackpoints)
CALL write_cif_line(cvals(1:1),'MAXTRACKPOINTS')

CALL conv_to_chars(cvals(1),nocampaigns)
CALL write_cif_line(cvals(1:1),'NOCAMPAIGNS')

CALL conv_to_chars(cvals(1),nocirclepoints)
CALL write_cif_line(cvals(1:1),'NOCIRCLEPOINTS')

CALL conv_to_chars(cvals(1),nodredgingsites)
CALL write_cif_line(cvals(1:1),'NODREDGINGSITES')

CALL conv_to_chars(cvals(1),norelocationsites)
CALL write_cif_line(cvals(1:1),'NORELOCATIONSITES')

!
!4. Write end of block
!----------------------
!

WRITE (iunit,'(A)') cifend
ciflinenum = ciflinenum + 1

CALL log_timer_out()


RETURN

END SUBROUTINE write_cif_vars_dar

!========================================================================

SUBROUTINE write_cif_vars_morph
!************************************************************************
!
! *write_cif_vars_morph* write the CIF block with parameters for the dredging
!                        and morphological model
!
! Author - Kevin Delecluyse (IMDC)
!
! Version - @(COHERENS)Sediment_Parameters.f90  V2.11.1
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - close_filepars, conv_to_chars, open_filepars, write_cif_line
!
!************************************************************************
!
USE iopars
USE morphpars
USE morphswitches
USE paralpars
USE cif_routines, ONLY: conv_to_chars, write_cif_line
USE inout_routines, ONLY: close_filepars, open_filepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: iunit


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_cif_vars_morph' 
CALL log_timer_in()

!
!1. Write block header
!--------------------
!

iunit = ciffile%iunit
WRITE (iunit,'(A)') cifend//TRIM(cif_block_names(icif_morph))

!
!2. Switches
!-----------
!

CALL conv_to_chars(cvals(1),iopt_morph_active_layer)
CALL write_cif_line(cvals(1:1),'IOPT_MORPH_ACTIVE_LAYER')

CALL conv_to_chars(cvals(1),iopt_morph_avalanching)
CALL write_cif_line(cvals(1:1),'IOPT_MORPH_AVALANCHING')

CALL conv_to_chars(cvals(1),iopt_morph_corr)
CALL write_cif_line(cvals(1:1),'IOPT_MORPH_CORR')

CALL conv_to_chars(cvals(1),iopt_morph_fixed_layer)
CALL write_cif_line(cvals(1:1),'IOPT_MORPH_FIXED_LAYER')

CALL conv_to_chars(cvals(1),iopt_morph_hdiff)
CALL write_cif_line(cvals(1:1),'IOPT_MORPH_HDIFF')

CALL conv_to_chars(cvals(1),iopt_morph_tidal_scheme)
CALL write_cif_line(cvals(1:1),'IOPT_MORPH_TIDAL_SCHEME')

CALL conv_to_chars(cvals(1),iopt_morph_time_int)
CALL write_cif_line(cvals(1:1),'IOPT_MORPH_TIME_INT')

CALL conv_to_chars(cvals(1),iopt_morph_vert_fluxes)
CALL write_cif_line(cvals(1:1),'IOPT_MORPH_VERT_FLUXES')

!
!3. Parameters
!-------------
!

CALL conv_to_chars(cvals(1),accelstep)
CALL write_cif_line(cvals(1:1),'ACCELSTEP')

CALL conv_to_chars(cvals(1),alpha_jameson)
CALL write_cif_line(cvals(1:1),'ALPHA_JAMESON')

CALL conv_to_chars(cvals(1),average_bedform_height)
CALL write_cif_line(cvals(1:1),'AVERAGE_BEDFORM_HEIGHT')

CALL conv_to_chars(cvals(1),average_bedform_length)
CALL write_cif_line(cvals(1:1),'AVERAGE_BEDFORM_LENGTH')

CALL conv_to_chars(cvals(1),bed_porosity_cst)
CALL write_cif_line(cvals(1:1),'BED_POROSITY_CST')

CALL conv_to_chars(cvals(1),dune_celerity)
CALL write_cif_line(cvals(1:1),'DUNE_CELERITY')

CALL conv_to_chars(cvals(1),dyn_angle_dry)
CALL write_cif_line(cvals(1:1),'DYN_ANGLE_DRY')

CALL conv_to_chars(cvals(1),dyn_angle_wet)
CALL write_cif_line(cvals(1:1),'DYN_ANGLE_WET')

CALL conv_to_chars(cvals(1),icaval)
CALL write_cif_line(cvals(1:1),'ICAVAL')

CALL conv_to_chars(cvals(1),icmorph)
CALL write_cif_line(cvals(1:1),'ICMORPH')

CALL conv_to_chars(cvals(1),icsedbal)
CALL write_cif_line(cvals(1:1),'ICSEDBAL')

CALL conv_to_chars(cvals(1),k1_harris)
CALL write_cif_line(cvals(1:1),'K1_HARRIS')

CALL conv_to_chars(cvals(1),k2_harris)
CALL write_cif_line(cvals(1:1),'K2_HARRIS')

CALL conv_to_chars(cvals(1),k2_jameson)
CALL write_cif_line(cvals(1:1),'K2_JAMESON')

CALL conv_to_chars(cvals(1),k4_jameson)
CALL write_cif_line(cvals(1:1),'K4_JAMESON')

CALL conv_to_chars(cvals(1),max_iterations_aval)
CALL write_cif_line(cvals(1:1),'MAX_ITERATIONS_AVAL')

CALL conv_to_chars(cvals(1),morph_factor)
CALL write_cif_line(cvals(1:1),'MORPH_FACTOR')

CALL conv_to_chars(cvals(1),morph_steps)
CALL write_cif_line(cvals(1:1),'MORPH_STEPS')

CALL conv_to_chars(cvals(1),nstep_hydro)
CALL write_cif_line(cvals(1:1),'NSTEP_HYDRO')

CALL conv_to_chars(cvals(1),number_tidal_steps)
CALL write_cif_line(cvals(1:1),'NUMBER_TIDAL_STEPS')

CALL conv_to_chars(cvals(1),similarity_range)
CALL write_cif_line(cvals(1:1),'SIMILARITY_RANGE')

CALL conv_to_chars(cvals(1),stat_angle_dry)
CALL write_cif_line(cvals(1:1),'STAT_ANGLE_DRY')

CALL conv_to_chars(cvals(1),stat_angle_wet)
CALL write_cif_line(cvals(1:1),'STAT_ANGLE_WET')

CALL conv_to_chars(cvals(1),trough_probability)
CALL write_cif_line(cvals(1:1),'TROUGH_PROBABILITY')

!
!4. Write end of block
!----------------------
!

WRITE (iunit,'(A)') cifend
ciflinenum = ciflinenum + 1

CALL log_timer_out()


RETURN

END SUBROUTINE write_cif_vars_morph

!========================================================================

SUBROUTINE write_cif_vars_sed
!************************************************************************
!
! *write_cif_vars_morph* write the CIF block with parameters for the sediment
!                        model
!
! Author - IMDC, Alexander Breugem
!
! Version - @(COHERENS)Sediment_Parameters.f90  V2.11.1
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - conv_to_chars, write_cif_line
!
!************************************************************************
!
USE iopars
USE paralpars
USE sedpars
USE sedswitches
USE cif_routines, ONLY: conv_to_chars, write_cif_line
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: iunit


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_cif_vars_sed' 
CALL log_timer_in()

!
!1. Write block header
!--------------------
!

iunit = ciffile%iunit
WRITE (iunit,'(A)') cifend//TRIM(cif_block_names(icif_sed))

!
!2. Switches
!-----------
!

CALL conv_to_chars(cvals(1),iopt_sed_bbc)
CALL write_cif_line(cvals(1:1),'IOPT_SED_BBC')

CALL conv_to_chars(cvals(1),iopt_sed_bbc_eq)
CALL write_cif_line(cvals(1:1),'IOPT_SED_BBC_EQ')

CALL conv_to_chars(cvals(1),iopt_sed_bbc_ref)
CALL write_cif_line(cvals(1:1),'IOPT_SED_BBC_REF')

CALL conv_to_chars(cvals(1),iopt_sed_bedeq)
CALL write_cif_line(cvals(1:1),'IOPT_SED_BEDEQ')

CALL conv_to_chars(cvals(1),iopt_sed_beta)
CALL write_cif_line(cvals(1:1),'IOPT_SED_BETA')

CALL conv_to_chars(cvals(1),iopt_sed_bstres)
CALL write_cif_line(cvals(1:1),'IOPT_SED_BSTRES')

CALL conv_to_chars(cvals(1),iopt_sed_bstres_cr)
CALL write_cif_line(cvals(1:1),'IOPT_SED_BSTRES_CR')

CALL conv_to_chars(cvals(1),iopt_sed_dens_grad)
CALL write_cif_line(cvals(1:1),'IOPT_SED_DENS_GRAD')

CALL conv_to_chars(cvals(1),iopt_sed_filter)
CALL write_cif_line(cvals(1:1),'IOPT_SED_FILTER')

CALL conv_to_chars(cvals(1),iopt_sed_hiding)
CALL write_cif_line(cvals(1:1),'IOPT_SED_HIDING')

CALL conv_to_chars(cvals(1),iopt_sed_ws_floc)
CALL write_cif_line(cvals(1:1),'IOPT_SED_WS_FLOC')

CALL conv_to_chars(cvals(1),iopt_sed_ws_hindset)
CALL write_cif_line(cvals(1:1),'IOPT_SED_WS_HINDSET')

CALL conv_to_chars(cvals(1),iopt_sed_median)
CALL write_cif_line(cvals(1:1),'IOPT_SED_MEDIAN')

CALL conv_to_chars(cvals(1),iopt_sed_mode)
CALL write_cif_line(cvals(1:1),'IOPT_SED_MODE')

CALL conv_to_chars(cvals(1),iopt_sed_nodim)
CALL write_cif_line(cvals(1:1),'IOPT_SED_NODIM')

CALL conv_to_chars(cvals(1),iopt_sed_obc_flux)
CALL write_cif_line(cvals(1:1),'IOPT_SED_OBC_FLUX')

CALL conv_to_chars(cvals(1),iopt_sed_rough)
CALL write_cif_line(cvals(1:1),'IOPT_SED_ROUGH')

CALL conv_to_chars(cvals(1),iopt_sed_slope)
CALL write_cif_line(cvals(1:1),'IOPT_SED_SLOPE')

CALL conv_to_chars(cvals(1),iopt_sed_source_impl)
CALL write_cif_line(cvals(1:1),'IOPT_SED_SOURCE_IMPL')

CALL conv_to_chars(cvals(1),iopt_sed_toteq)
CALL write_cif_line(cvals(1:1),'IOPT_SED_TOTEQ')

CALL conv_to_chars(cvals(1),iopt_sed_type)
CALL write_cif_line(cvals(1:1),'IOPT_SED_TYPE')

CALL conv_to_chars(cvals(1),iopt_sed_vadv)
CALL write_cif_line(cvals(1:1),'IOPT_SED_VADV')

CALL conv_to_chars(cvals(1),iopt_sed_wave_diff)
CALL write_cif_line(cvals(1:1),'IOPT_SED_WAVE_DIFF')

CALL conv_to_chars(cvals(1),iopt_sed_ws)
CALL write_cif_line(cvals(1:1),'IOPT_SED_WS')

CALL conv_to_chars(cvals(1),iopt_sed_ws_lim)
CALL write_cif_line(cvals(1:1),'IOPT_SED_WS_LIM')

!
!3. Parameters
!-------------
!
!3.1 Sediment model
!------------------
!

CALL conv_to_chars(cvals(1),icsed)
CALL write_cif_line(cvals(1:1),'ICSED')

CALL conv_to_chars(cvals(1),maxitbartnicki)
CALL write_cif_line(cvals(1:1),'MAXITBARTNICKI')

CALL conv_to_chars(cvals(1),nb)
CALL write_cif_line(cvals(1:1),'NB')

CALL conv_to_chars(cvals(1),nf)
CALL write_cif_line(cvals(1:1),'NF')

CALL conv_to_chars(cvals(1),nrquad_sed)
CALL write_cif_line(cvals(1:1),'NRQUAD_SED')

CALL conv_to_chars(cvals(1),nrquad_wav)
CALL write_cif_line(cvals(1:1),'NRQUAD_WAV')

CALL conv_to_chars(cvals(1),alpha_VR)
CALL write_cif_line(cvals(1:1),'ALPHA_VR')

CALL conv_to_chars(cvals(1),a_leussen)
CALL write_cif_line(cvals(1:1),'A_LEUSSEN')

CALL conv_to_chars(cvals(1),beta_sed_cst)
CALL write_cif_line(cvals(1:1),'BETA_SED_CST')

CALL conv_to_chars(cvals(1),b_leussen)
CALL write_cif_line(cvals(1:1),'B_LEUSSEN')

CALL conv_to_chars(cvals(1),cgel)
CALL write_cif_line(cvals(1:1),'CGEL')

CALL conv_to_chars(cvals(1),cmax)
CALL write_cif_line(cvals(1:1),'CMAX')

CALL conv_to_chars(cvals(1),coef_bed_grad)
CALL write_cif_line(cvals(1:1),'COEF_BED_GRAD')

CALL conv_to_chars(cvals(1),floc_VR_max)
CALL write_cif_line(cvals(1:1),'FLOC_VR_MAX')

CALL conv_to_chars(cvals(1),floc_VR_min)
CALL write_cif_line(cvals(1:1),'FLOC_VR_MIN')

CALL conv_to_chars(cvals(1),height_c_cst)
CALL write_cif_line(cvals(1:1),'HEIGHT_C_CST')

CALL conv_to_chars(cvals(1),n_RichZaki)
CALL write_cif_line(cvals(1:1),'N_RICHZAKI')

CALL conv_to_chars(cvals(1),parth_coef)
CALL write_cif_line(cvals(1:1),'PARTH_COEF')

CALL conv_to_chars(cvals(1),parth_exp)
CALL write_cif_line(cvals(1:1),'PARTH_EXP')

CALL conv_to_chars(cvals(1),zrough_grain)
CALL write_cif_line(cvals(1:1),'ZROUGH_GRAIN')

CALL conv_to_chars(cvals(1),zrough_sed_cst)
CALL write_cif_line(cvals(1:1),'ZROUGH_SED_CST')

CALL conv_to_chars(cvals(1),z0_coef)
CALL write_cif_line(cvals(1:1),'Z0_COEF')

CALL conv_to_chars(cvals(1),wslimfac)
CALL write_cif_line(cvals(1:1),'WSLIMFAC')

CALL conv_to_chars(cvals(1),wu_exp)
CALL write_cif_line(cvals(1:1),'WU_EXP')

!   
!3.2 Flocculation model
!----------------------
!

CALL conv_to_chars(cvals(1),agg_alpha)
CALL write_cif_line(cvals(1:1),'AGG_ALPHA')

CALL conv_to_chars(cvals(1),brk_es)
CALL write_cif_line(cvals(1:1),'BRK_ES')

CALL conv_to_chars(cvals(1),brk_ffy)
CALL write_cif_line(cvals(1:1),'BRK_FFY')

CALL conv_to_chars(cvals(1),brk_frac)
CALL write_cif_line(cvals(1:1),'BRK_FRAC')

CALL conv_to_chars(cvals(1),brk_p)
CALL write_cif_line(cvals(1:1),'BRK_P')

CALL conv_to_chars(cvals(1),brk_q)
CALL write_cif_line(cvals(1:1),'BRK_Q')

CALL conv_to_chars(cvals(1),floc_ncmax)
CALL write_cif_line(cvals(1:1),'FLOC_NCMAX')

CALL conv_to_chars(cvals(1),floc_ncmin)
CALL write_cif_line(cvals(1:1),'FLOC_NCMIN')

CALL conv_to_chars(cvals(1),nfrdim)
CALL write_cif_line(cvals(1:1),'NFRDIM')

!
!4. Write end of block
!----------------------
!

WRITE (iunit,'(A)') cifend
ciflinenum = ciflinenum + 1

CALL log_timer_out()


RETURN

END SUBROUTINE write_cif_vars_sed
