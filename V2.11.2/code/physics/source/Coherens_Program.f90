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

PROGRAM coherens_main

!************************************************************************
!
! *coherens_main* Coherens main program
!
! Author - Patrick Luyten
!
! Version - 4 Nov 2008  @(COHERENS)Coherens_Program.f90  V2.11.2
!
! $Date: 2018-09-28 13:01:37 +0200 (Fri, 28 Sep 2018) $
!
! $Revision: 1188 $
!
! Description - 
!
! Reference -
!
! External calls - astronomical_tides, baroclinic_gradient, biological_model,
!   coherens_end, coherens_start, combine_particle_phys_data, dar_equation,
!   define_particle_phys_data, define_phsics, discharges, equation_of_state,
!   harmonic_analysis, horizontal_diff_coeffs, hydrodynamic_equations,
!   initialise_model, mask_function, meteo_input, morphology_accel,
!   morphology_equation, particle_concentrations, particle_drift,
!   particle_model, particle_trajects, salinity_equation, sediment_equation,
!   simulation_end, simulation_start, store_depths_old, temperature_equation,
!   time_averages, time_series, usrdef_output, usrdef_part_output,
!   vertical_diff_coefs, wave_input, wave_model, weirs_mask, write_bioics,
!   write_morphics, write_partics, write_particle_phys_data, write_phsics,
!   write_sedics
!
! Module calls - comms_barrier, error_abort, loop_index, monitor_files,
!                update_time, write_error_message
!
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE syspars
USE timepars
USE comms_MPI, ONLY: comms_barrier
USE error_routines, ONLY: error_abort, write_error_message
USE inout_routines, ONLY: monitor_files
USE time_routines, ONLY: update_time
USE utility_routines, ONLY: loop_index

IMPLICIT NONE


pglev = 1
procname(1) = 'coherens_main'

!
!1. Start MPI (if required)
!--------------------------
!

CALL coherens_start

!
!2. Title of next simulation or exit
!-----------------------------------
!

1000 CALL simulation_start
IF (.NOT.next_simul) GOTO 2000 

!
!2.1 Start coupled model(s)
!--------------------------
!
!---particle model
IF ((iopt_part_model.EQ.2.AND.modelid.EQ.modelidpart).OR.&
   & iopt_part_model.GE.3) THEN
   CALL particle_model
   GOTO 3000
ENDIF
!---wave model  
IF (iopt_waves_model.GT.0.AND.modelid.EQ.modelidwav) THEN
   CALL wave_model
   GOTO 3000
ENDIF

!
!3. Initialise simulation
!-------------------------
!

CALL initialise_model

!---terminate simulation in case of cold start
IF (cold_start) GOTO 3000

!
!4. Tidal acceleration method
!----------------------------
!

IF (iopt_tidal_accel.EQ.1) THEN 
   CALL morphology_accel
   GOTO 3000
ENDIF

!
!4. Time loop
!------------
!

nt_400: DO nt=1,nstep

!
!4.1 Restart log file
!--------------------
!

   CALL monitor_files

!
!4.2 Date/time
!-------------
!

   CALL update_time

!
!4.3 Define new initial conditions
!---------------------------------
!

   IF (physinit) CALL define_phsics
 
!
!
!4.4 Meteo input
!---------------
!

   IF (iopt_meteo.EQ.1) CALL meteo_input

!
!4.5 Surface waves
!-----------------
!

   IF (iopt_waves.EQ.1) CALL wave_input

!
!4.6 Discharge locations
!-----------------------
!

   IF (iopt_dischr.EQ.1) CALL discharges

!
!4.7 Predictor step
!------------------
!

   IF (predstep) THEN
!     ---store old water depths
      CALL store_depths_old
!     ---dynamic mask
      IF (iopt_fld.EQ.2) CALL mask_function
      IF (iopt_weibar.EQ.1) CALL weirs_mask
      IF (iopt_mode_3D.EQ.2) THEN
!        ---vertical diffusion coefficients
         IF (iopt_vdif_coef.GT.1) CALL vertical_diff_coefs
!        ---baroclinic pressure gradient
         IF (iopt_dens_grad.GT.0) CALL baroclinic_gradient
      ENDIF
      IF (iopt_vdif_rot.EQ.1) CALL vertical_diff_rot
!     ---horizontal diffusion coefficients
      IF (iopt_hdif_coef.GT.0) CALL horizontal_diff_coefs
   ENDIF

!
!4.8 Astronomical tide
!---------------------
!
   
   IF (iopt_astro_tide.EQ.1) CALL astronomical_tides
   
!
!4.9 Currents and elevations
!---------------------------
!

   IF (iopt_curr.GT.0) CALL hydrodynamic_equations

!
!4.10 Corrector step
!------------------
!
!  ---physics
   IF (corrstep) THEN
      IF (iopt_temp.EQ.2) CALL temperature_equation
      IF (iopt_sal.EQ.2) CALL salinity_equation
      IF (iopt_dar.GT.0) CALL dar_equation
      IF (iopt_biolgy.GT.0) CALL biological_model
   ENDIF

!  ---sediments   
   IF (iopt_sed.GT.0) CALL sediment_equation
   IF (iopt_morph.GT.0) CALL morphology_equation
   
!  ---density and expansion coefficients
   IF (corrstep.AND.iopt_dens.GT.0.OR.iopt_sed.GT.0) THEN
      CALL equation_of_state
   ENDIF
   
!  ---particle model   
   IF (corrstep) THEN
      IF (iopt_part_write.EQ.1) CALL write_particle_phys_data(CDateTime)
      IF (iopt_part_model.EQ.1) THEN
         CALL define_particle_phys_data
      ELSEIF (iopt_part_model.EQ.2) THEN
         CALL combine_particle_phys_data
      ENDIF
      IF (iopt_part_model.EQ.1) CALL particle_drift
      IF (iopt_part_model.GT.0) THEN
         CALL particle_concentrations
      ENDIF
   ENDIF
   
!
!4.11 Model output
!-----------------
!
!  ---user-defined
   IF (iopt_out_tsers.EQ.1) CALL time_series
   IF (iopt_out_avrgd.EQ.1.AND.corrstep) CALL time_averages
   IF (iopt_out_anal.EQ.1.AND.corrstep) CALL harmonic_analysis
   IF (iopt_part_model.EQ.1) CALL particle_trajects
   CALL usrdef_output
   IF (iopt_part_model.EQ.1) CALL usrdef_part_output

!  ---restart files
   CALL write_phsics
   IF (iopt_sed.GT.0) CALL write_sedics
   IF (iopt_morph.GT.0) CALL write_morphics
   IF (iopt_biolgy.GT.0) CALL write_bioics
   IF (iopt_part_model.EQ.1) CALL write_partics

ENDDO nt_400

!
!5. Synchronise processes for new simulation
!-------------------------------------------
!

3000 IF (parallel_set) CALL comms_barrier(comm_world_MPI)

!
!6. Finalise simulation
!----------------------
!

CALL simulation_end


GOTO 1000

!
!7. Finalise program
!-------------------
!

2000 CALL coherens_end


STOP 'Main program terminated'


END PROGRAM coherens_main
