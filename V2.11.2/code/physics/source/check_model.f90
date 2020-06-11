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

MODULE check_model
!************************************************************************
!
! *check_model* Check model parameters and arrays for errors
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.11.2
!
! $Date: 2018-09-28 13:01:37 +0200 (Fri, 28 Sep 2018) $
!
! $Revision: 1188 $
!
! Description - 
!
! Reference -
!
! Generic routines - check_out_filepars
!
! Routines - check_grid_arrays, check_hrel_coords, check_initial_conditions,
!            check_init_params, check_mod_filepars, check_mod_params,
!            check_out_gpars, check_partition, check_struct_locs,
!            check_statlocs, check_variables
!
!************************************************************************
!
USE iopars

IMPLICIT NONE


INTERFACE check_out_filepars
   MODULE PROCEDURE check_out_filepars_1d, check_out_filepars_2d
END INTERFACE

CONTAINS

!========================================================================

SUBROUTINE check_grid_arrays
!************************************************************************
!
! *check_grid_arrays* Check grid arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.9
!
! Description - 
!
! Reference -
!
! Module calls - error_abort, error_alloc, error_lbound_arr, error_limits_arr
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE paralpars
USE physpars
USE switches
USE syspars
USE error_routines, ONLY: nerrs, error_abort, error_alloc, error_lbound_arr, &
                        & error_limits_arr
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: ci, cii, cj, cjj
INTEGER :: i, ii, j, jj, k, noerrs


IF (.NOT.master) RETURN

procname(pglev+1) = 'check_grid_arrays'
CALL log_timer_in()

!
!1. Longitude and latitude
!-------------------------
!

IF (iopt_grid_sph.EQ.1) THEN

!  ---X-coordinates
   noerrs = COUNT(gxcoordglb(1:nc,1:nr).LT.-180.0.OR.&
                & gxcoordglb(1:nc,1:nr).GT.180.0)
   IF (noerrs.GT.0) THEN
      i_110: DO i=1,nc
      j_110: DO j=1,nr
         CALL error_limits_arr(gxcoordglb(i,j),'gxcoordglb',-180.0,180.0,2,&
                             & indx=(/i,j/))
      ENDDO j_110
      ENDDO i_110
   ENDIF

!  ---Y-coordinates
   noerrs = COUNT(gycoordglb(1:nc,1:nr).LT.-90.0.OR.&
                & gycoordglb(1:nc,1:nr).GT.90.0)
   IF (noerrs.GT.0) THEN
      i_120: DO i=1,nc
      j_120: DO j=1,nr
         CALL error_limits_arr(gycoordglb(i,j),'gycoordglb',-90.0,90.0,2,&
                             & indx=(/i,j/))
      ENDDO j_120
      ENDDO i_120
   ENDIF

ENDIF

!
!2. Sigma-coordinates
!--------------------
!

i_210: DO i=1,nc-1
j_210: DO j=1,nr-1
   IF (seapointglb(i,j)) THEN
      k_211: DO k=2,nz+1
         CALL error_lbound_arr(gscoordglb(i,j,k),'gscoordglb',&
                             & gscoordglb(i,j,k-1),.FALSE.,3,&
                             & indx=(/i,j,k/))
      ENDDO k_211
   ENDIF
ENDDO j_210
ENDDO i_210

!
!3. Open boundary locations
!--------------------------
!
!3.1 Check if locations are inside the domain
!--------------------------------------------
! 
!---U-nodes
ii_311: DO ii=1,nobu
   i = iobu(ii); j = jobu(ii)
   CALL error_limits_arr(i,'iobu',1,nc,1,indx=(/ii/))
   CALL error_limits_arr(j,'jobu',1,nr,1,indx=(/ii/))
ENDDO ii_311

!---V-nodes
jj_312: DO jj=1,nobv
   i = iobv(jj); j = jobv(jj)
   CALL error_limits_arr(i,'iobv',1,nc,1,indx=(/jj/))
   CALL error_limits_arr(j,'jobv',1,nr,1,indx=(/jj/))
ENDDO jj_312

!
!3.2 Check if locations are at coastal locations
!-----------------------------------------------
!
!---U-nodes
ii_321: DO ii=1,nobu
   i = iobu(ii); j = jobu(ii)
   IF ((seapointglb(i-1,j).AND.seapointglb(i,j)).OR.&
     & .NOT.(seapointglb(i-1,j).OR.seapointglb(i,j))) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cii,'(I12)') ii; cii = ADJUSTL(cii)
         WRITE (ci,'(I12)') i; ci = ADJUSTL(ci)
         WRITE (cj,'(I12)') j; cj = ADJUSTL(cj)
         WRITE (ioerr,'(A)') 'Invalid location of U-open boundary point '&
                       & //TRIM(cii)//' at: '//'('//TRIM(ci)//','//TRIM(cj)//')'
         WRITE (ioerr,'(A)') 'Open boundaries must be located between a'//&
                           & ' dry and a wet cell'
      ENDIF
   ENDIF
ENDDO ii_321

!---V-nodes
jj_322: DO jj=1,nobv
   i = iobv(jj); j = jobv(jj)
   IF ((seapointglb(i,j-1).AND.seapointglb(i,j)).OR.&
     & .NOT.(seapointglb(i,j-1).OR.seapointglb(i,j))) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cjj,'(I12)') jj; cjj = ADJUSTL(cjj)
         WRITE (ci,'(I12)') i; ci = ADJUSTL(ci)
         WRITE (cj,'(I12)') j; cj = ADJUSTL(cj)
         WRITE (ioerr,'(A)') 'Invalid location of V-open boundary point '&
                       & //TRIM(cjj)//' at: '//'('//TRIM(ci)//','//TRIM(cj)//')'
         WRITE (ioerr,'(A)') 'Open boundaries must be located between a'//&
                           & ' dry and a wet cell'
      ENDIF
   ENDIF
ENDDO jj_322

CALL error_abort('check_grid_arrays',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_grid_arrays

!========================================================================

SUBROUTINE check_init_params
!************************************************************************
!
! *check_init_params* Check initial parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.11
!
! Description - 
!
! Reference -
!
! Module calls - error_value_var, error_limits_var, warning_reset_var
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE error_routines, ONLY: error_value_var, error_limits_var, warning_reset_var
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
CHARACTER (LEN=12) :: cval


procname(pglev+1)  ='check_init_params'
CALL log_timer_in()

!---MCT can only be used with MPI 
IF (iopt_MPI.EQ.0.AND.iopt_MCT.EQ.1) THEN
   CALL warning_reset_var(iopt_MCT,'iopt_MCT',0)
ENDIF

!---MPI and MCT must activated if wave coupling is enabled
IF (iopt_waves_model.GT.0) THEN
   IF (iopt_MPI.EQ.0) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (ioerr,'(A)') 'MPI must be enabled if wave '//&
                           & 'coupling is activated'
      ENDIF
   ENDIF
   IF (iopt_MCT.EQ.0) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (ioerr,'(A)') 'The MCT library must be enabled '//&
                           & 'if wave coupling is activated'
      ENDIF
   ENDIF
ENDIF

!---switch for activating particle model
CALL error_limits_var(iopt_part_model,'iopt_part_model',0,4)

!---MPI must be enabled if the particle model is activated and COHERENS
!  runs in parallel
IF (iopt_part_model.EQ.2) THEN
   IF (iopt_MPI.EQ.0) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (ioerr,'(A)') 'MPI must be enabled if the particle model '//&
                           & 'switch equals 2'
      ENDIF
   ENDIF
ENDIF

!---no parallel or wave coupling if particle model runs offline
IF (iopt_part_model.GT.2) THEN
   IF (iopt_MPI.EQ.1) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (ioerr,'(A)') 'MPI must be disabled if the particle model '//&
                           & 'runs offline'
      ENDIF
   ENDIF
   IF (iopt_waves_model.EQ.1) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (ioerr,'(A)') 'Wave coupling must be disabled if the '//&
                           & ' particle model runs offline'
      ENDIF
   ENDIF
ENDIF

!---sum of processes per model must be equal to the number in MPI_comm_world
IF (nprocscoh+nprocspart+nprocswav.NE.npworld) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(I12)') nprocscoh; cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Number of processes attributed to flow model: '&
                      & //TRIM(cval)
      WRITE (cval,'(I12)') nprocspart; cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Number of processes attributed to particle model: '&
                      & //TRIM(cval)
      WRITE (cval,'(I12)') nprocswav; cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Number of processes attributed to wave model: '//&
                      & TRIM(cval)
      WRITE (cval,'(I12)') npworld; cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Total number of processes must be equal to: '//&
                      & TRIM(cval)
   ENDIF
ENDIF

!---MCT activated for the wave coupling
IF (master) THEN
   CALL error_limits_var(iopt_waves_model,'iopt_waves_model',0,1)
   IF (iopt_waves_model.GT.0) THEN
      CALL error_value_var(iopt_MCT,'iopt_MCT',1)
      IF (iopt_MCT.EQ.0) THEN
         nerrs = nerrs + 1
         IF (errchk.AND.nerrs.LE.maxerrors) THEN
            WRITE (ioerr,'(A)') 'Coupling with an external wave model is '//&
                              & 'only allowed if the MCT library is activated'
         ENDIF
      ENDIF
   ENDIF
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE check_init_params

!========================================================================

SUBROUTINE check_hrel_coords(hcoords,nhdat,varname)
!************************************************************************
!
! *check_hrel_coords* Check the components of a vector of type 'HRelativeCoords'
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.9
!
! Description - 
!
! Reference -
!
! Module calls - error_abort, error_limits_arr_struc, num_proc
!
!************************************************************************
!
USE datatypes
USE paralpars  
USE syspars
USE error_routines, ONLY: error_abort, error_limits_arr_struc
USE grid_routines, ONLY: num_proc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: varname
INTEGER, INTENT(IN) :: nhdat
TYPE (HRelativeCoords), INTENT(IN), DIMENSION(nhdat) :: hcoords

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*varname*  CHAR    Name of coordinate array
!*nhdat*    INTEGER Number of horozontal data points
!*hcoords*  DERIVED Horizontal coordinate array structure
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: ci, cj, cl
INTEGER :: i, iproc, j, l


IF (.NOT.master.OR.nhdat.EQ.0) RETURN

procname(pglev+1)  ='check_hrel_coords_1d'
CALL log_timer_in(logname=TRIM(procname(pglev+1))//': '//TRIM(varname))

!
!1. Index coordinates
!--------------------
!

l_110: DO l=1,nhdat
   i = hcoords(l)%icoord; j = hcoords(l)%jcoord
   iproc = num_proc(i,j)
   IF (iproc.LE.0) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cl,'(I12)') l; cl = ADJUSTL(cl)
         WRITE (ci,'(I12)') i; ci = ADJUSTL(ci)
         WRITE (cj,'(I12)') j; cj = ADJUSTL(cj)
         WRITE (ioerr,'(A)') 'Invalid value for element '//TRIM(cl)//&
                           & ' of coordinate array '//TRIM(varname)
         WRITE (ioerr,'(A)') 'Data point ('//TRIM(ci)//','//TRIM(cj)//&
                           & ') outside computational domain'
      ENDIF
   ENDIF
ENDDO l_110
CALL error_abort('check_hrel_coords',ierrno_alloc)

!
!2. Relative x- and y-coordinates
!--------------------------------
!

l_210: DO l=1,nhdat
   CALL error_limits_arr_struc(hcoords(l)%weights(1,1),'hcoords',&
                            & 'weights(1,1)',0.0,1.0,1,(/l/))
   CALL error_limits_arr_struc(hcoords(l)%weights(2,1),'hcoords',&
                            & 'weights(2,1)',0.0,1.0,1,(/l/))
   CALL error_limits_arr_struc(hcoords(l)%weights(1,2),'hcoords',&
                            & 'weights(1,2)',0.0,1.0,1,(/l/))
   CALL error_limits_arr_struc(hcoords(l)%weights(2,2),'hcoords',&
                            & 'weights(2,2)',0.0,1.0,1,(/l/))
ENDDO l_210

CALL error_abort('check_hrel_coords',ierrno_alloc)

CALL log_timer_out()


RETURN

END SUBROUTINE check_hrel_coords

!========================================================================

SUBROUTINE check_initial_conditions
!************************************************************************
!
! *check_initial_conditions* Check initial conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.10.2
!
! Description - 
!
! Reference -
!
! Module calls - error_abort, error_lbound_arr, max_vars, min_vars,
!                warning_ubound_var_real
!
!************************************************************************
!
USE depths
USE fluxes
USE grid
USE gridpars
USE modids
USE paralpars
USE switches
USE timepars
USE error_routines, ONLY: error_abort, error_lbound_arr, &
                        & warning_ubound_var_real
USE paral_utilities, ONLY: max_vars, min_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=15) :: cval
INTEGER :: i, iglb, j, jglb, noerrs
REAL :: depmax, dtmax, gmax


procname(pglev+1) = 'check_initial_conditions'
CALL log_timer_in()

!
!1. Bottom drag coefficient or roughness length
!----------------------------------------------
!

IF (iopt_bstres_form.EQ.2) THEN

!  ---bottom drag coefficient
   IF (iopt_bstres_drag.LE.2) THEN
      noerrs = COUNT(bdragcoefatc(1:ncloc,1:nrloc).LE.0.0.AND.&
                   & seapoint(1:ncloc,1:nrloc))
      IF (noerrs.GT.0) THEN
         i_410: DO i=1,ncloc
         j_410: DO j=1,nrloc
            IF (seapoint(i,j)) THEN
               iglb = nc1loc+i-1; jglb = nr1loc+j-1
               CALL error_lbound_arr(bdragcoefatc(i,j),'bdragcoefatc',0.0,&
                                  & .FALSE.,2,indx=(/iglb,jglb/))
            ENDIF
         ENDDO j_410
         ENDDO i_410
      ENDIF

!  ---roughness length
   ELSEIF (iopt_sed.EQ.0) THEN
      noerrs = COUNT(zroughatc(1:ncloc,1:nrloc).LE.0.0.AND.&
                   & seapoint(1:ncloc,1:nrloc))
      IF (noerrs.GT.0) THEN
         i_420: DO i=1,ncloc
         j_420: DO j=1,nrloc
            IF (seapoint(i,j)) THEN
               iglb = nc1loc+i-1; jglb = nr1loc+j-1
               CALL error_lbound_arr(zroughatc(i,j),'zroughatc',0.0,.FALSE.,2,&
                                   & indx=(/iglb,jglb/))
            ENDIF
         ENDDO j_420
         ENDDO i_420
      ENDIF
   ENDIF

ENDIF

CALL error_abort('check_initial_conditions',ierrno_inival)

!
!2. CFL-limit for 2-D time step
!------------------------------
!

CALL max_vars(gaccatc(1:ncloc,1:nrloc),gmax,iarr_gaccatc,mask=maskatc_int,&
            & commall=.TRUE.)
CALL max_vars(depmeanatc(1:ncloc,1:nrloc),depmax,iarr_deptotatc,&
            & mask=maskatc_int,commall=.TRUE.)
dtmax = 0.5*MIN(delxatc_min,delyatc_min)/SQRT(gmax*depmax)
IF (master.AND.loglev1.GT.0) THEN
   WRITE (cval,'(G15.7)') dtmax; cval = ADJUSTL(cval)
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//' dtmax = '//cval
   WRITE (cval,'(G15.7)') MIN(delxatc_min,delyatc_min); cval = ADJUSTL(cval)
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//' delhmin = '//cval
   WRITE (cval,'(G15.7)') MAX(delxatc_max,delyatc_max); cval = ADJUSTL(cval)
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//' delhmax = '//cval
ENDIF
IF (iopt_mode_2D.EQ.2.AND.iopt_hydro_impl.EQ.0) THEN
   CALL warning_ubound_var_real(delt2d,'delt2d',dtmax,.TRUE.)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE check_initial_conditions

!========================================================================

SUBROUTINE check_mod_filepars
!************************************************************************
!
! *check_mod_filepars* Check parameters of model files
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.7.1
!
! Description -
!
! Reference -
!
! Module calls - check_time_limits_arr_struc, error_abort,
!                error_limits_arr_struc, error_limits_var,
!                error_vals_arr_struc_char, error_vals_arr_struc_int,
!                error_value_arr_struc, warning_reset_arr_struc
!
!************************************************************************
!
USE nestgrids
USE paralpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: check_time_limits_arr_struc, error_abort, &
                        & error_limits_arr_struc, error_limits_var, &
                        & error_vals_arr_struc_char, error_vals_arr_struc_int, &
                        & error_value_arr_struc, warning_reset_arr_struc
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
CHARACTER (LEN=12) :: cdesc, cfil
INTEGER :: idesc, ifil, iotype, iset


IF (.NOT.master) RETURN

procname(pglev+1) = 'check_mod_filepars'
CALL log_timer_in()

!
!1. Status
!---------
!

IF (norestarts.EQ.0) THEN
   CALL warning_reset_arr_struc(modfiles(io_fincon,ics_phys,2)%status,&
                             & 'modfiles','%status','0',3,&
                             & (/io_inicon,ics_phys,2/))

!ABR: Bugfix  with respect to the status
   IF (iopt_sed.GT.0) THEN
      CALL warning_reset_arr_struc(modfiles(io_inicon,ics_sed,2)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_inicon,ics_sed,2/))
   ENDIF

   IF (iopt_biolgy.GT.0) THEN
      CALL warning_reset_arr_struc(modfiles(io_fincon,ics_bio,2)%status,&
                                & 'modfiles','%status','0',3,&
                                & (/io_inicon,ics_bio,2/))
   ENDIF

ENDIF

IF (iopt_meteo.EQ.1.AND.surfacegrids(igrd_meteo,1)%nhtype.GT.1.AND.&
     & surfacegrids(igrd_meteo,1)%nhtype.LT.4) THEN
   IF (surfacegrids(igrd_meteo,1)%surcoords.EQ.1) THEN
      CALL error_vals_arr_struc_char(modfiles(io_metabs,1,1)%status,'modfiles',&
                                 & '%status','"R" "N"',3,(/io_metabs,1,1/))
   ELSEIF (surfacegrids(igrd_meteo,1)%surcoords.EQ.2) THEN
      CALL error_vals_arr_struc_char(modfiles(io_metrel,1,1)%status,'modfiles',&
                                 & '%status','"R" "N"',3,(/io_metabs,1,1/))
   ENDIF
ENDIF

IF (iopt_waves.GT.0) THEN
   CALL error_vals_arr_struc_int(surfacegrids(igrd_waves,1)%nhtype,&
                              & 'surfacegrid','%nhtype',2,3,(/0,4/),(/io_wavsur,1,1/))
   IF (surfacegrids(igrd_waves,1)%surcoords.EQ.1) THEN
      CALL error_vals_arr_struc_char(modfiles(io_wavabs,1,1)%status,'modfiles',&
                                 & '%status','"R" "N"',3,(/io_wavabs,1,1/))
   ELSEIF (surfacegrids(igrd_meteo,1)%surcoords.EQ.2) THEN
      CALL error_vals_arr_struc_char(modfiles(io_wavrel,1,1)%status,'modfiles',&
                                 & '%status','"R" "N"',3,(/io_wavabs,1,1/))
   ENDIF
ENDIF

IF (iopt_temp_sbc.GT.1.AND.surfacegrids(igrd_sst,1)%nhtype.GT.1.AND.&
     & surfacegrids(igrd_sst,1)%nhtype.LT.4) THEN
   IF (surfacegrids(igrd_sst,1)%surcoords.EQ.1) THEN
      CALL error_vals_arr_struc_char(modfiles(io_sstabs,1,1)%status,'modfiles',&
                                 & '%status','"R" "N"',3,(/io_sstabs,1,1/))
   ELSEIF (surfacegrids(igrd_sst,1)%surcoords.EQ.2) THEN
      CALL error_vals_arr_struc_char(modfiles(io_sstrel,1,1)%status,'modfiles',&
                                 & '%status','"R" "N"',3,(/io_sstabs,1,1/))
   ENDIF
ENDIF

IF (iopt_obc_relax.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_rlxobc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_rlxobc,1,1/))
ENDIF

IF (iopt_drycel.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_drycel,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_drycel,1,1/))
ENDIF

IF (iopt_thndam.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_thndam,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_thndam,1,1/))
ENDIF

IF (iopt_weibar.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_weibar,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_weibar,1,1/))
ENDIF

IF (iopt_dischr.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_disspc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_weibar,1,1/))
   CALL error_vals_arr_struc_char(modfiles(io_disloc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_weibar,1,1/))
ENDIF

IF (iopt_nests.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_nstspc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_nstspc,1,1/))
ENDIF

IF (iopt_meteo.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_metsur,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_metsur,1,1/))
ENDIF

IF (iopt_waves.GT.0) THEN
   CALL error_vals_arr_struc_char(modfiles(io_wavsur,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_wavsur,1,1/))
ENDIF

IF (iopt_temp_sbc.GT.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_sstsur,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_sstsur,1,1/))
ENDIF

IF (iopt_sur_1D.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_1uvsur,1,1)%status,&
                      & 'modfiles','%status','"R" "N"',3,(/io_1uvsur,1,1/))
ENDIF

IF (iopt_obc_2D.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_2uvobc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_2uvobc,1,1/))
ENDIF

IF (iopt_obc_2D_tang.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_2xyobc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_2xyobc,1,1/))
ENDIF

IF (iopt_obc_3D.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_3uvobc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_3uvobc,1,1/))
ENDIF

IF (iopt_obc_3D_tang.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_3xyobc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_3xyobc,1,1/))
ENDIF

IF (iopt_obc_sal.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_salobc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_salobc,1,1/))
ENDIF

IF (iopt_obc_temp.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_tmpobc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_tmpobc,1,1/))
ENDIF

IF (iopt_obc_sed.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_sedobc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_sedobc,1,1/))
ENDIF
IF (iopt_obc_bio.EQ.1) THEN
   CALL error_vals_arr_struc_char(modfiles(io_bioobc,1,1)%status,'modfiles',&
                          & '%status','"R" "N"',3,(/io_bioobc,1,1/))
ENDIF

idesc_120: DO idesc=1,MaxIOTypes
ifil_120: DO ifil=1,MaxIOFiles
   CALL error_vals_arr_struc_char(modfiles(idesc,ifil,1)%status,'modfiles',&
                          & '%status','"0" "R" "N"',3,(/ifil,idesc,1/))
   CALL error_vals_arr_struc_char(modfiles(idesc,ifil,2)%status,'modfiles',&
                          & '%status','"0" "W"',3,(/ifil,idesc,2/))
ENDDO ifil_120
ENDDO idesc_120

!
!2. Form
!-------
!

iotype_210: DO iotype=1,2
ifil_210: DO ifil=1,MaxIOFiles
idesc_210: DO idesc=1,MaxIOTypes
   IF (modfiles(idesc,ifil,iotype)%status.NE.'0') THEN
      CALL error_vals_arr_struc_char(modfiles(idesc,ifil,iotype)%form,&
                             & 'modfiles','%form','"A" "U" "N"',3,&
                             & (/idesc,ifil,iotype/))
   ENDIF
ENDDO idesc_210
ENDDO ifil_210
ENDDO iotype_210

!
!3. Filename
!-----------
!

ifil_310: DO ifil=1,MaxIOFiles
idesc_310: DO idesc=1,MaxIOTypes
   IF (idesc.NE.io_inicon.AND.modfiles(idesc,ifil,1)%defined.AND.&
                            & modfiles(idesc,ifil,2)%defined) THEN
      IF (TRIM(modfiles(idesc,ifil,1)%filename).EQ.&
        & TRIM(modfiles(idesc,ifil,2)%filename)) THEN
         nerrs = nerrs + 1
         IF (errchk.AND.nerrs.LE.maxerrors) THEN
            WRITE (cdesc,'(I12)') idesc; cdesc = ADJUSTL(cdesc)
            WRITE (cfil,'(I12)') ifil; cfil = ADJUSTL(cfil)
            WRITE (ioerr,'(A)') &
               & 'Invalid value for character component filename and element ('&
               & //TRIM(cdesc)//','//TRIM(cfil)//') of array modfiles: '&
               & //TRIM(modfiles(idesc,ifil,2)%filename)
            WRITE (ioerr,'(A)') &
               & 'Input and output forcing files must have different names'
         ENDIF
      ENDIF
   ENDIF
ENDDO idesc_310
ENDDO ifil_310

!
!4. Time limits
!--------------
!

idesc_410: DO idesc=1,MaxIOTypes
   ifil_411: DO ifil=2,MaxIOFiles
      SELECT CASE (idesc)
         CASE (io_1uvsur,io_2uvobc,io_2xyobc,io_3uvobc,io_3xyobc,io_salobc,&
             & io_tmpobc,io_sedobc,io_bioobc)
            IF (ifil.GT.1.AND.modfiles(idesc,ifil,1)%status.NE.'0') THEN
               CALL check_time_limits_arr_struc(modfiles(idesc,ifil,1)%tlims,&
                                   & 'modfiles','tlims',0,0,3,(/idesc,ifil,1/))
            ENDIF
      END SELECT
   ENDDO ifil_411

   ifil_412: DO ifil=1,MaxIOFiles
      SELECT CASE (idesc)
         CASE (io_metsur,io_sstsur,io_wavsur,io_biosur,io_disloc,io_discur,&
             & io_dissal,io_distmp)
            IF (modfiles(idesc,ifil,1)%status.NE.'0') THEN
               CALL check_time_limits_arr_struc(modfiles(idesc,ifil,1)%tlims,&
                                    & 'modfiles','tlims',0,0,3,(/idesc,ifil,1/))
               CALL error_value_arr_struc(modfiles(idesc,ifil,1)%tlims(1),&
                                   & 'modfiles','tlims(1)',0,3,(/idesc,ifil,1/))
            ENDIF
      END SELECT
   ENDDO ifil_412
   
   IF (iopt_nests.EQ.1) THEN
      iset_413:DO iset=1,nonestsets
         SELECT CASE (idesc)
            CASE (io_2uvnst,io_2xynst)
               IF (modfiles(idesc,iset,2)%status.NE.'0') THEN
                  CALL check_time_limits_arr_struc(&
                                     & modfiles(idesc,iset,1)%tlims,&
                                    & 'modfiles','tlims',0,0,3,&
                                    & (/idesc,iset,2/))
               ENDIF
            CASE (io_3uvnst,io_3xynst,io_salnst,io_tmpnst,io_sednst,io_bionst)
               IF (modfiles(idesc,iset,2)%status.NE.'0') THEN
                  CALL check_time_limits_arr_struc(&
                                     & modfiles(idesc,iset,1)%tlims,&
                                    & 'modfiles','tlims',0,0,3,&
                                    & (/idesc,iset,2/),flag3d=.TRUE.)
               ENDIF
         END SELECT
      ENDDO iset_413
   ENDIF
ENDDO idesc_410

!
!5. Time skips
!-------------
!

idesc_510: DO idesc=1,MaxIOTypes
ifil_510: DO ifil=1,MaxIOFiles
   IF (modfiles(idesc,ifil,2)%defined) THEN
      CALL error_limits_arr_struc(modfiles(idesc,ifil,2)%tskips,'modfiles',&
                               & '%tskips',1,9,3,(/idesc,ifil,2/))
   ENDIF
ENDDO ifil_510
ENDDO idesc_510

!
!6. Counter for log file
!-----------------------
!

IF (iopt_part_model.LE.2.AND.loglev1.GT.0) THEN
   CALL error_limits_var(runlog_count,'runlog_count',1,nstep)
ENDIF

CALL error_abort('check_mod_filepars',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_mod_filepars

!========================================================================

SUBROUTINE check_mod_params
!************************************************************************
!
! *check_mod_params* Check model parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Module calls - add_secs_to_date, check_filepars, check_mult_counter,
!                error_abort, error_lbound_var, error_limits_arr,
!                error_limits_arr_struc, error_limits_var, error_vals_var_int,
!                error_value_var, mult_index
!
!************************************************************************
!
USE gridpars
USE paralpars
USE physpars
USE structures
USE switches
USE syspars
USE tide
USE timepars
USE error_routines, ONLY: nerrs, check_mult_counter, error_abort, &
                        & error_lbound_var, error_limits_arr, &
                        & error_limits_arr_struc, error_limits_var, &
                        & error_vals_var_int, error_value_var
USE time_routines, ONLY: add_secs_to_date,  log_timer_in, log_timer_out
USE utility_routines, ONLY: mult_index

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: cnum
CHARACTER (LEN=15) :: cval
CHARACTER (LEN=lentime) :: cdatetimex
INTEGER :: icon, icwav, igrd, l


IF (.NOT.master) RETURN

procname(pglev+1) = 'check_mod_params'
CALL log_timer_in()

!
!1. Grid dimensions
!------------------
!

CALL error_lbound_var(nc,'nc',0,.FALSE.)
CALL error_lbound_var(nr,'nr',0,.FALSE.)
CALL error_lbound_var(nz,'nz',0,.FALSE.)

!
!2. Switches
!-----------
!
!2.1 Ranges
!----------
!
!---grid
CALL error_limits_var(iopt_grid_htype,'iopt_grid_htype',1,3)
CALL error_limits_var(iopt_grid_nodim,'iopt_grid_nodim',1,3)
CALL error_limits_var(iopt_grid_sph,'iopt_grid_sph',0,1)
CALL error_limits_var(iopt_grid_vtype,'iopt_grid_vtype',1,3)
IF (iopt_grid_vtype.EQ.2) THEN
   CALL error_vals_var_int(iopt_grid_vtype_transf,'iopt_grid_vtype_transf',&
                        &  4,(/0,11,12,13/))
ELSEIF (iopt_grid_vtype.EQ.3) THEN
   CALL error_vals_var_int(iopt_grid_vtype_transf,'iopt_grid_vtype_transf',&
                         & 2,(/0,21/))
ENDIF

!---interpolation
CALL error_limits_var(iopt_arrint_depths,'iopt_arrkndint',1,2)
CALL error_limits_var(iopt_arrint_hreg,'iopt_arrint_hreg',0,1)
CALL error_limits_var(iopt_arrint_vreg,'iopt_arrint_vreg',0,1)
CALL error_limits_var(iopt_arrint_3D,'iopt_arrint_3D',0,1)

!---hydrodynamics
CALL error_limits_var(iopt_curr,'iopt_curr',0,2)
!CALL error_limits_var(iopt_curr_wfall,'iopt_curr_wfall',1,2)
CALL error_value_var(iopt_curr_wfall,'iopt_curr_wfall',1)

!---density
CALL error_limits_var(iopt_dens,'iopt_dens',0,3)
CALL error_limits_var(iopt_dens_convect,'iopt_dens_convect',0,1)
CALL error_limits_var(iopt_dens_grad,'iopt_dens_grad',0,3)
CALL error_limits_var(iopt_sal,'iopt_sal',0,2)
CALL error_limits_var(iopt_temp,'iopt_temp',0,2)
CALL error_limits_var(iopt_temp_optic,'iopt_temp_optic',0,1)
CALL error_limits_var(iopt_temp_sbc,'iopt_temp_sbc',1,3)

!---biology
CALL error_limits_var(iopt_biolgy,'iopt_biolgy',0,1)

!---sediments
CALL error_limits_var(iopt_dar,'iopt_dar',0,1)
CALL error_limits_var(iopt_morph,'iopt_morph',0,1)
CALL error_limits_var(iopt_sed,'iopt_sed',0,2)
CALL error_limits_var(iopt_tidal_accel,'iopt_tidal_accel',0,1)

!---particle model
CALL error_limits_var(iopt_part_write,'iopt_part_write',0,1)

!---bottom stress
CALL error_limits_var(iopt_bstres_drag,'iopt_bstres_drag',0,4)
CALL error_limits_var(iopt_bstres_form,'iopt_bstres_form',0,2)
CALL error_limits_var(iopt_bstres_nodim,'iopt_bstres_nodim',2,3)
CALL error_limits_var(iopt_bstres_waves_bfric,'iopt_bstres_waves_bfric',0,4)

!---transport
CALL error_limits_var(iopt_transp_full,'iopt_transp_full',0,1)

!---advection
CALL error_limits_var(iopt_adv_scal,'iopt_adv_scal',0,3)
CALL error_limits_var(iopt_adv_turb,'iopt_adv_turb',0,3)
CALL error_limits_var(iopt_adv_tvd,'iopt_adv_tvd',1,2)
CALL error_limits_var(iopt_adv_2D,'iopt_adv_2D',0,3)
CALL error_limits_var(iopt_adv_3D,'iopt_adv_3D',0,3)

!---diffusion
CALL error_limits_var(iopt_hdif_coef,'iopt_hdif_coef',0,2)
CALL error_limits_var(iopt_hdif_lim,'iopt_hdif_lim',0,3)
CALL error_limits_var(iopt_hdif_scal,'iopt_hdif_scal',0,3)
CALL error_limits_var(iopt_hdif_turb,'iopt_hdif_turb',0,1)
CALL error_limits_var(iopt_hdif_2D,'iopt_hdif_2D',0,1)
CALL error_limits_var(iopt_hdif_3D,'iopt_hdif_3D',0,1)
CALL error_limits_var(iopt_kinvisc,'iopt_kinvisc',0,1)
CALL error_limits_var(iopt_vdif_coef,'iopt_vdif_coef',0,3)
CALL error_limits_var(iopt_vdif_rot,'iopt_vdif_rot',0,1)

!---turbulence
CALL error_limits_var(iopt_turb_alg,'iopt_turb_alg',1,6)
CALL error_limits_var(iopt_turb_dis_bbc,'iopt_turb_bbc',1,2)
CALL error_limits_var(iopt_turb_dis_sbc,'iopt_turb_sbc',1,2)
CALL error_limits_var(iopt_turb_iwlim,'iopt_turb_iwlim',0,2)
CALL error_limits_var(iopt_turb_kinvisc,'iopt_turb_kinvisc',0,1)
CALL error_limits_var(iopt_turb_lmix,'iopt_turb_lmix',1,4)
CALL error_limits_var(iopt_turb_ntrans,'iopt_turb_ntrans',0,2)
CALL error_limits_var(iopt_turb_param,'iopt_turb_param',1,2)
CALL error_limits_var(iopt_turb_stab_form,'iopt_turb_stab_form',1,3)
CALL error_limits_var(iopt_turb_stab_lev,'iopt_turb_stab_lev',1,2)
CALL error_limits_var(iopt_turb_stab_mod,'iopt_turb_stab_mod',1,6)
CALL error_limits_var(iopt_turb_stab_tke,'iopt_turb_stab_tke',1,3)
CALL error_limits_var(iopt_turb_tke_bbc,'iopt_turb_bbc',1,2)
CALL error_limits_var(iopt_turb_tke_sbc,'iopt_turb_sbc',1,2)

!---drying/wetting
CALL error_limits_var(iopt_fld,'iopt_fld',0,2)
CALL error_limits_var(iopt_fld_alpha,'iopt_fld_alpha',0,2)

!---structures
CALL error_limits_var(iopt_drycel,'iopt_drycel',0,1)
CALL error_limits_var(iopt_thndam,'iopt_thndam',0,1)
CALL error_limits_var(iopt_weibar,'iopt_weibar',0,1)

!---discharges
CALL error_limits_var(iopt_dischr,'iopt_dischr',0,1)
CALL error_limits_var(iopt_dischr_land,'iopt_dischr_land',1,2)

!---time-integration
CALL error_limits_var(iopt_cor_impl,'iopt_cor_impl',0,2)
CALL error_limits_var(iopt_hydro_impl,'iopt_hydro_impl',0,1)
CALL error_limits_var(iopt_scal_depos,'iopt_scal_depos',0,2)
CALL error_limits_var(iopt_vadv_impl,'iopt_vadv_impl',0,2)
CALL error_limits_var(iopt_vdif_impl,'iopt_vdif_impl',0,2)

!---open boundary conditions
CALL error_limits_var(iopt_obc_advflux,'iopt_obc_advflux',1,2)
CALL error_limits_var(iopt_obc_advrlx,'iopt_obc_advrlx',0,1)
CALL error_limits_var(iopt_obc_bio,'iopt_obc_bio',0,1)
CALL error_limits_var(iopt_obc_invbar,'iopt_obc_invbar',0,1)
CALL error_limits_var(iopt_obc_relax,'iopt_obc_relax',0,1)
CALL error_limits_var(iopt_obc_sal,'iopt_obc_sal',0,1)
CALL error_limits_var(iopt_obc_sed,'iopt_obc_sed',0,1)
CALL error_limits_var(iopt_obc_temp,'iopt_obc_temp',0,1)
CALL error_limits_var(iopt_obc_th,'iopt_obc_th',0,1)
CALL error_limits_var(iopt_obc_2D,'iopt_obc_2D',0,1)
CALL error_limits_var(iopt_obc_2D_tang,'iopt_obc_2D_tang',0,1)
CALL error_limits_var(iopt_obc_3D,'iopt_obc_3D',0,1)
CALL error_limits_var(iopt_obc_3D_tang,'iopt_obc_3D_tang',0,1)

!---tides
CALL error_limits_var(iopt_astro_anal,'iopt_astro_anal',0,1)
CALL error_limits_var(iopt_astro_pars,'iopt_astro_pars',0,2)
CALL error_limits_var(iopt_astro_tide,'iopt_astro_tide',0,1)

!---1-D applications
CALL error_limits_var(iopt_sur_1D,'iopt_sur_1D',0,1)
   
!---meteo
CALL error_limits_var(iopt_meteo,'iopt_meteo',0,1)
CALL error_limits_var(iopt_meteo_data,'iopt_meteo_data',0,2)
CALL error_limits_var(iopt_meteo_heat,'iopt_meteo_heat',0,1)
CALL error_limits_var(iopt_meteo_precip,'iopt_meteo_precip',0,2)
CALL error_limits_var(iopt_meteo_pres,'iopt_meteo_pres',0,1)
CALL error_limits_var(iopt_meteo_stres,'iopt_meteo_stres',0,1)

!---surface waves
CALL error_limits_var(iopt_waves,'iopt_waves',0,1)
CALL error_limits_var(iopt_waves_couple,'iopt_waves_couple',0,3)
CALL error_limits_var(iopt_waves_curr,'iopt_waves_curr',0,1)
CALL error_limits_var(iopt_waves_dissip,'iopt_waves_dissip',0,1)
CALL error_limits_var(iopt_waves_extrapol,'iopt_waves_extrapol',0,1)
CALL error_limits_var(iopt_waves_form,'iopt_waves_form',1,2)
CALL error_limits_var(iopt_waves_pres,'iopt_waves_pres',0,1)

!---surface fluxes
CALL error_limits_var(iopt_sflux_pars,'iopt_sflux_pars',0,9)
CALL error_limits_var(iopt_sflux_precip,'iopt_sflux_precip',0,3)
CALL error_limits_var(iopt_sflux_qlong,'iopt_sflux_qlong',1,3)
CALL error_limits_var(iopt_sflux_qshort,'iopt_sflux_qshort',1,3)
CALL error_limits_var(iopt_sflux_strat,'iopt_sflux_strat',0,1)

!---surface waves
CALL error_limits_var(iopt_waves,'iopt_waves',0,1)

!---nesting
CALL error_limits_var(iopt_nests,'iopt_nests',0,1)

!---MPI
CALL error_limits_var(iopt_MPI_abort,'iopt_MPI_abort',0,1)
CALL error_limits_var(iopt_MPI_comm_all,'iopt_MPI_comm_all',1,2)
CALL error_limits_var(iopt_MPI_comm_coll,'iopt_MPI_comm_coll',0,1)
CALL error_vals_var_int(iopt_MPI_comm_exch,'iopt_MPI_comm_exch',3,(/1,2,5/))
CALL error_limits_var(iopt_MPI_comm_full,'iopt_MPI_comm_full',0,1)
CALL error_limits_var(iopt_MPI_comm_gath,'iopt_MPI_comm_gath',1,2)
CALL error_limits_var(iopt_MPI_comm_scat,'iopt_MPI_comm_scat',1,2)
CALL error_limits_var(iopt_MPI_partit,'iopt_MPI_partit',1,2)
CALL error_limits_var(iopt_MPI_sync,'iopt_MPI_sync',0,1)

!---output
CALL error_limits_var(iopt_out_anal,'iopt_out_anal',0,1)
CALL error_limits_var(iopt_out_avrgd,'iopt_out_avrgd',0,1)
CALL error_limits_var(iopt_out_tsers,'iopt_out_tsers',0,1)

!---multigrid
CALL error_limits_var(iopt_mg_cycle,'iopt_mg_cycle',1,2)
CALL error_limits_var(iopt_mg_prolong,'iopt_mg_prolong',1,2)
CALL error_limits_var(iopt_mg_smoother,'iopt_mg_smoother',1,2)

!---netCDF
CALL error_limits_var(iopt_CDF_abort,'iopt_CDF_abort',0,1)
CALL error_limits_var(iopt_CDF_fill,'iopt_CDF_fill',0,1)
CALL error_limits_var(iopt_CDF_format,'iopt_CDF_format',1,2)
CALL error_limits_var(iopt_CDF_shared,'iopt_CDF_shared',0,1)
CALL error_limits_var(iopt_CDF_sync,'iopt_CDF_sync',0,1)
CALL error_limits_var(iopt_CDF_tlim,'iopt_CDF_tlim',1,2)

!---random generator
CALL error_limits_var(iopt_rng_seed,'iopt_rng_seed',1,2)

!---verification procedure
!CALL error_limits_var(iopt_verif,'iopt_verif',0,1)

!
!2.2 Incompatibilities
!---------------------
!
!---meteo
IF (iopt_temp.LT.2.AND.iopt_meteo_precip.EQ.2) THEN
   nerrs = nerrs + 1
   IF (errchk) THEN
      WRITE (ioerr,'(A)') 'Invalid value for switches iopt_meteo_precip or '//&
                        & 'iopt_temp'
      WRITE (ioerr,'(A)') 'iopt_meteo_precip = 2'
      WRITE (ioerr,'(A,I1)') 'iopt_temp = ', iopt_temp
   ENDIF
ENDIF

!---astronomical tide
IF (iopt_astro_tide.EQ.1) CALL error_lbound_var(nconastro,'nconastro',0,.FALSE.)

!---bottom stress
IF (iopt_bstres_form.EQ.2.AND.iopt_bstres_drag.EQ.0) THEN
   CALL error_lbound_var(iopt_bstres_drag,'iopt_bstres_drag',0,.FALSE.)
ENDIF
IF (iopt_grid_nodim.EQ.2.AND.iopt_bstres_form.EQ.2.AND.&
  & iopt_bstres_nodim.EQ.3) THEN
   CALL error_value_var(iopt_bstres_nodim,'iopt_bstres_nodim',2)
ENDIF

!---open boundary conditions
IF (iopt_obc_2D.EQ.0) THEN
   CALL error_value_var(iopt_obc_2D_tang,'iopt_obc_2D_tang',0)
ENDIF
IF (iopt_obc_2D.EQ.0.OR.iopt_obc_3D.EQ.0) THEN
   CALL error_value_var(iopt_obc_3D_tang,'iopt_obc_3D_tang',0)
ENDIF

!---structures
IF (iopt_weibar.EQ.1.AND.iopt_bstres_form.LT.2) THEN
   CALL error_value_var(iopt_bstres_form,'iopt_bstres_form',2)
ENDIF

!
!3. Date/time parameters
!-----------------------
!

!---non-zero 3D time step
CALL error_lbound_var(ic3d,'ic3d',0,.FALSE.)

!---integer number of external time steps within simulation period
IF (iopt_part_model.LE.2) THEN
   CALL add_secs_to_date(CStartDateTime,cdatetimex,nstep,delt2d)
   IF (cdatetimex.NE.CEndDateTime) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval,'(G15.7)') delt2d
         cval = ADJUSTL(cval)
         WRITE (ioerr,'(A)') 'Invalid external time step: '//TRIM(cval)
         WRITE (ioerr,'(A)') 'Simulation period must contain an integer '//&
                          &  'number of external time steps'
      ENDIF
   ENDIF
ENDIF

!---integer number of external time steps within simulation period
IF (.NOT.mult_index(nstep,ic3d)) THEN
   nerrs = nerrs + 1
   IF (errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cval,'(G15.7)') delt2d*ic3d
      cval = ADJUSTL(cval)
      WRITE (ioerr,'(A)') 'Invalid internal time step: '//TRIM(cval)
      WRITE (ioerr,'(A)') 'Simulation period must contain an integer number'//&
                        & ' of internal time steps'
   ENDIF
ENDIF

!---counter for convective adjustment
IF (iopt_dens_convect.EQ.1) THEN
   CALL check_mult_counter(iccvt,'iccvt',ic3d,matchzero=.TRUE.)
ENDIF
   
!---time zone
CALL error_limits_var(time_zone,'time_zone',-12.0,12.0)

!---wave current interaction
IF (iopt_waves.GT.0) THEN
   IF (modfiles(io_wavsur,1,1)%status.NE.'0') THEN
      icwav = modfiles(io_wavsur,1,1)%tlims(3)
      IF (.NOT.mult_index(icwav,ic3d)) THEN
         nerrs = nerrs + 1
         IF (errchk.AND.nerrs.LE.maxerrors) THEN
            WRITE (cval,'(G15.7)') delt2d*icwav
            cval = ADJUSTL(cval)
            WRITE (ioerr,'(A)') 'Invalid wave current coupling-time step: '//&
                              & TRIM(cval)
            WRITE (ioerr,'(A)') 'Must be a multiple of the 3-D time step'
         ENDIF
      ENDIF
   ENDIF
ENDIF

!
!4. I/O file properties
!----------------------
!

CALL check_mod_filepars

!
!5. Model parameters
!-------------------
!
!---number of open boundaries
CALL error_lbound_var(nosbu,'nosbu',0,.TRUE.)
CALL error_lbound_var(nosbv,'nosbv',0,.TRUE.)
CALL error_lbound_var(nrvbu,'nrvbu',0,.TRUE.)
CALL error_lbound_var(nrvbv,'nrvbv',0,.TRUE.)

!---number of tidal frequencies
CALL error_limits_var(nconobc,'nconobc',0,MaxConstituents)
CALL error_limits_var(nconastro,'nconastro',0,MaxAstroTides)

!---grid
CALL error_limits_var(dlat_ref,'dlat_ref',-90.0,90.0)
CALL error_limits_var(dlon_ref,'dlon_ref',-180.0,180.0)
CALL error_limits_var(dlon_ref_anal,'dlon_ref_anal',-180.0,180.0)
CALL error_limits_var(dlon_ref_obc,'dlon_ref_obc',-180.0,180.0)

!---reference values
CALL error_lbound_var(atmpres_ref,'atmpres_ref',0.0,.FALSE.)
CALL error_lbound_var(sal_ref,'sal_ref',0.0,.TRUE.)
CALL error_lbound_var(temp_ref,'temp_ref',0.0,.TRUE.)

!---diffusion coefficients
IF (iopt_hdif_coef.EQ.1) THEN
   CALL error_lbound_var(hdifmom_cst,'hdifmom_cst',0.0,.TRUE.)
   CALL error_lbound_var(hdifscal_cst,'hdifscal_cst',0.0,.TRUE.)
ENDIF
IF (iopt_hdif_coef.EQ.2) THEN
   CALL error_lbound_var(smag_coef_mom,'smag_coef_mom',0.0,.FALSE.)
   CALL error_lbound_var(smag_coef_scal,'smag_coef_scal',0.0,.FALSE.)
ENDIF
IF (iopt_vdif_coef.EQ.1) THEN
   CALL error_lbound_var(vdifmom_cst,'vdifmom_cst',0.0,.TRUE.)
   CALL error_lbound_var(vdifscal_cst,'vdifscal_cst',0.0,.TRUE.)
ELSEIF (iopt_vdif_coef.GT.0) THEN
   CALL error_lbound_var(vdifmom_cst,'vdifmom_cst',0.0,.TRUE.)
   CALL error_lbound_var(vdifscal_cst,'vdifscal_cst',0.0,.TRUE.)
ENDIF
IF (iopt_hdif_scal.EQ.3) THEN
   CALL error_lbound_var(rho_crit_iso,'rho_crit_iso',0.0,.FALSE.)
   CALL error_lbound_var(skewdiff_cst,'skewdiff_cst',0.0,.TRUE.)
   CALL error_lbound_var(slopemax_iso,'slopemax_iso',0.0,.FALSE.)

ENDIF
CALL error_lbound_var(kinvisc_cst,'kinvisc_cst',0.0,.FALSE.)

!---bottom fluxes
IF (iopt_bstres_form.EQ.1) THEN
   CALL error_lbound_var(bdraglin,'bdraglin',0.0,.TRUE.)
ELSEIF (iopt_bstres_form.EQ.2) THEN
   CALL error_lbound_var(zbtoz0lim,'zbtoz0lim',1.0,.TRUE.)
ENDIF

!---optical parameters
IF (iopt_temp_optic.EQ.1) THEN
   CALL error_lbound_var(optattcoef1_cst,'optattcoef1_cst',0.0,.TRUE.)
   CALL error_lbound_var(optattcoef2_cst,'optattcoef2_cst',0.0,.TRUE.)
ENDIF

!---thickness of the WBL
IF (iopt_waves.EQ.1.AND.iopt_bstres_waves_bfric.EQ.0) THEN
   CALL error_lbound_var(wavethick_cst,'wavethick_cst',0.0,.FALSE.)
ENDIF

!---number of user output files
IF (iopt_out_tsers.EQ.1) THEN
   CALL error_lbound_var(nosetstsr,'nosetstsr',1,.TRUE.)
ENDIF
IF (iopt_out_avrgd.EQ.1) THEN
   CALL error_lbound_var(nosetsavr,'nosetsavr',1,.TRUE.)
ENDIF
IF (iopt_out_anal.EQ.1) THEN
   CALL error_lbound_var(nosetsanal,'nosetsanal',1,.TRUE.)
ENDIF

!---mask criteria
IF (iopt_fld.EQ.2) THEN
   l_511: DO l=1,nofldmasks
      CALL error_limits_arr(fld_mask(l),'fld_mask',0,1,1,indx=(/l/))
   ENDDO l_511
ENDIF

!---tidal indices
IF (nconobc.GT.0) THEN
   icon_512: DO icon=1,nconobc
      CALL error_limits_arr(index_obc(icon),'index_obc',0,MaxConstituents,&
                          & 1,indx=(/icon/))
   ENDDO icon_512
ENDIF
IF (nconastro.GT.0) THEN
   icon_513: DO icon=1,nconastro
      CALL error_limits_arr(index_astro(icon),'index_astro',0,MaxAstroTides,&
                          & 1,indx=(/icon/))
   ENDDO icon_513
ENDIF

!---structure parameters
IF (iopt_drycel.EQ.1) CALL error_lbound_var(numdry,'numdry',0,.FALSE.)
IF (iopt_thndam.EQ.1) THEN
   IF (numthinu+numthinv.LE.0) THEN
      nerrs = nerrs + 1
      IF (errchk) THEN
         WRITE (ioerr,'(A)') 'Number of thin dams must be greater than zero'
         WRITE (cnum,'(I12)') numthinu 
         WRITE (ioerr,'(A)') 'numthinu = '//TRIM(cnum)
         WRITE (cnum,'(I12)') numthinv 
         WRITE (ioerr,'(A)') 'numthinv = '//TRIM(cnum)
      ENDIF
   ENDIF
ENDIF
IF (iopt_weibar.EQ.1) THEN
   IF (numwbaru+numwbarv.LE.0) THEN
      nerrs = nerrs + 1
      IF (errchk) THEN
         WRITE (ioerr,'(A)') 'Number of weirs must be greater than zero'
         WRITE (cnum,'(I12)') numwbaru 
         WRITE (ioerr,'(A)') 'numwbaru = '//TRIM(cnum)
         WRITE (cnum,'(I12)') numwbarv 
         WRITE (ioerr,'(A)') 'numwbarv = '//TRIM(cnum)
      ENDIF
   ENDIF
   CALL error_limits_var(wbarrlxu,'wbarrlxu',0.0,1.0)
   CALL error_limits_var(wbarrlxv,'wbarrlxv',0.0,1.0)
ENDIF

!---discharges
IF (iopt_dischr.EQ.1) CALL error_lbound_var(numdis,'numdis',0,.FALSE.)

!
!6. Surface grids
!----------------
!

IF (iopt_grid_sph.EQ.1) THEN

   igrd_610: DO igrd=1,MaxGridTypes
      CALL error_limits_arr_struc(surfacegrids(igrd,1)%x0dat,'surfacegrids',&
                               & '%x0dat',-180.0,180.0,2,(/igrd,1/))
      CALL error_limits_arr_struc(surfacegrids(igrd,1)%y0dat,'surfacegrids',&
                               & '%y0dat',-90.0,90.0,2,(/igrd,1/))
   ENDDO igrd_610
   IF (surfacegrids(igrd_model,1)%rotated) THEN
      CALL error_limits_arr_struc(surfacegrids(igrd_model,1)%y0rot,&
                               & 'surfacegrids','%y0rot',-90.0,90.0,2,&
                               & (/igrd_model,1/))
   ENDIF

ENDIF

IF (surfacegrids(igrd_model,1)%rotated) THEN
   CALL error_limits_arr_struc(surfacegrids(igrd_model,1)%gridangle,&
                            & 'surfacegrids','%gridangle',0.0,180.0,2,&
                            & (/igrd_model,1/))
ENDIF


CALL error_abort('check_mod_params',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_mod_params

!========================================================================

SUBROUTINE check_out_filepars_1d(filepars,arrname,maxvars)
!************************************************************************
!
! *check_out_filepars_1d* Check parameters of user output files
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.7.1
!
! Description - 
!
! Reference -
!
! Module calls - check_time_limits_arr_struc, error_abort,
!                error_limits_arr_struc, error_vals_arr_struc_char
!
!************************************************************************
!
USE datatypes
USE paralpars  
USE error_routines, ONLY: check_time_limits_arr_struc, error_abort, &
                        & error_limits_arr_struc, error_vals_arr_struc_char
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!* Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN), OPTIONAL :: maxvars
TYPE (FileParams), INTENT(INOUT), DIMENSION(:) :: filepars

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*filepars* DERIVED File parameters
!*arrname*  CHAR    Name of derived type array
!*maxvars*  INTEGER Maximum allowed number of variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, ndim


ndim = SIZE(filepars)
IF (.NOT.master.OR.ndim.EQ.0) RETURN

procname(pglev+1) = 'check_out_filepars_1d'
CALL log_timer_in()

n_110: DO n=1,ndim

!  ---status component
   CALL error_vals_arr_struc_char(filepars(n)%status,arrname,'%status',&
                          & '"0" "W"',1,(/n/))

!  ---form component
   CALL error_vals_arr_struc_char(filepars(n)%form,arrname,'%form',&
                               & '"A" "U" "N"',1,(/n/))

!  ---time limits
   CALL check_time_limits_arr_struc(filepars(n)%tlims,arrname,'%tlims',0,0,1,&
                                 & (/n/))

!  ---number of variables
   IF (PRESENT(maxvars)) THEN
      CALL error_limits_arr_struc(filepars(n)%novars,arrname,'%novars',&
                                & 1,maxvars,1,(/n/))
   ENDIF

ENDDO n_110

CALL error_abort('check_out_filepars_1d',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_out_filepars_1d

!========================================================================

SUBROUTINE check_out_filepars_2d(filepars,arrname,maxvars)
!************************************************************************
!
! *check_out_filepars_2d* Check parameters of user output files
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.7.1
!
! Description - 
!
! Reference -
!
! Module calls - check_time_limits_arr_struc, error_abort,
!                error_limits_arr_struc, error_vals_arr_struc_char
!
!************************************************************************
!
USE datatypes
USE paralpars  
USE error_routines, ONLY: check_time_limits_arr_struc, error_abort, &
                        & error_limits_arr_struc, error_vals_arr_struc_char
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!* Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
INTEGER, INTENT(IN), OPTIONAL :: maxvars
TYPE (FileParams), INTENT(INOUT), DIMENSION(:,:) :: filepars

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*filepars* DERIVED File parameters
!*arrname*  CHAR    Name of derived type array
!*maxvars*  INTEGER Maximum allowed number of variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n1, n2
INTEGER, DIMENSION(2) :: ndims


IF (.NOT.master.OR.SIZE(filepars).EQ.0) RETURN

procname(pglev+1) = 'check_out_filepars_2d'
CALL log_timer_in()
 
ndims = SHAPE(filepars)

n2_110: DO n2=1,ndims(2)
n1_110: DO n1=1,ndims(1)

!  ---status component
   CALL error_vals_arr_struc_char(filepars(n1,n2)%status,arrname,'%status',&
                          & '"0" "W" "P"',2,(/n1,n2/))

!  ---form component
   CALL error_vals_arr_struc_char(filepars(n1,n2)%form,arrname,'%form',&
                          & '"A" "U" "N"',2,(/n1,n2/))

!  ---time limits
   CALL check_time_limits_arr_struc(filepars(n1,n2)%tlims,arrname,'tlims',&
                                  & 0,0,2,(/n1,n2/))

!  ---number of variables
   IF (PRESENT(maxvars)) THEN
      CALL error_limits_arr_struc(filepars(n1,n2)%novars,arrname,'%novars',&
                                & 1,maxvars,2,(/n1,n2/))
   ENDIF

ENDDO n1_110
ENDDO n2_110

CALL error_abort('check_out_filepars_2d',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_out_filepars_2d

!========================================================================

SUBROUTINE check_out_gpars(outgpars,arrname,maxstats)
!************************************************************************
!
! *check_out_gpars* Check parameters of user output grids
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.10.2
!
! Description - 
!
! Reference -
!
! Module calls - check_space_limits_arr_struc, check_time_limits_struc,
!                error_abort, error_limits_arr_struc, error_vals_arr_struc_char
!
!************************************************************************
!
USE datatypes
USE gridpars
USE syspars
USE paralpars
USE timepars
USE error_routines, ONLY: check_space_limits_arr_struc, &
                        & check_time_limits_arr_struc, error_abort, &
                        & error_limits_arr_struc, error_vals_arr_struc_char
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
CHARACTER (LEN=*) :: arrname
INTEGER, INTENT(IN) :: maxstats
TYPE (OutGridParams), INTENT(IN), DIMENSION(:) :: outgpars

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*outgpars*  DERIVED  Output grid parameters
!*maxstats*  INTEGER  Maximum number of stations
!*arrname*   CHAR     Name of derived type array
!
!------------------------------------------------------------------------------
!
!* Local variables
!
LOGICAL :: gridded
INTEGER :: iset, nodim, nosets
INTEGER, DIMENSION(3) :: slims, tlims


nosets = SIZE(outgpars)
IF (.NOT.master.OR.nosets.EQ.0) RETURN

procname(pglev+1) = 'check_out_gpars'
CALL log_timer_in()

!
!1. Status
!---------
!

iset_110: DO iset=1,nosets
   CALL error_vals_arr_struc_char(outgpars(iset)%status,arrname,'%status',&
                          & '"W" "P"',1,(/iset/))
ENDDO iset_110

!
!2. Output times
!---------------
!

iset_210: DO iset=1,nosets
   tlims = outgpars(iset)%tlims
   CALL check_time_limits_arr_struc(tlims,arrname,'%tlims',0,nstep,1,(/iset/))
ENDDO iset_210

!
!3. Vertical coordinate
!----------------------
!

iset_310: DO iset=1,nosets
   CALL error_limits_arr_struc(outgpars(iset)%vcoord,arrname,'%vcoord',1,2,1,&
                            & (/iset/))
ENDDO iset_310

!
!4. Number of stations
!---------------------
!

iset_410: DO iset=1,nosets
   IF (.NOT.outgpars(iset)%gridded) THEN
      CALL error_limits_arr_struc(outgpars(iset)%nostats,arrname,'%nostats',1,&
                                & maxstats,1,(/iset/))
   ENDIF
ENDDO iset_410

!
!5. Output locations
!-------------------
!

iset_510: DO iset=1,nosets
   nodim = outgpars(iset)%nodim
   gridded = outgpars(iset)%gridded

!  ---X-direction
   IF (gridded.AND.nodim.GT.0) THEN
      slims = outgpars(iset)%xlims
      CALL check_space_limits_arr_struc(slims,arrname,'%xlims',nc-1,1,(/iset/))
   ENDIF

!  ---Y-direction
   IF (gridded.AND.nodim.GT.0) THEN
      slims = outgpars(iset)%ylims
      CALL check_space_limits_arr_struc(slims,arrname,'%ylims',nr-1,1,(/iset/))
   ENDIF

!  ---Z-direction
!   IF (nodim.EQ.3) THEN
!      slims = outgpars(iset)%zlims
!      CALL check_space_limits_arr_struc(slims,arrname,'%zlims',nz,1,(/iset/))
!   ENDIF

ENDDO iset_510

!
!6. Time format
!--------------
!

iset_610: DO iset=1,nosets
   CALL error_limits_arr_struc(outgpars(iset)%time_format,arrname,&
                            & '%time_format',0,8,1,(/iset/))
ENDDO iset_610

CALL error_abort('check_out_gpars',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_out_gpars

!========================================================================

SUBROUTINE check_partition
!************************************************************************
!
! *check_partition* Check domain decomposition
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.5.1
!
! Description - checks if all sea points are inside a process domain
!             - checks if all open boundaries are inside a process domain
!
! Reference -
!
! Module calls - error_abort, error_alloc, num_proc
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE paralpars
USE physpars
USE syspars
USE error_routines, ONLY: error_abort, error_alloc
USE grid_routines, ONLY: num_proc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=12) :: ci, cj
INTEGER :: i, ii, iproc, j, jj


IF (.NOT.master) RETURN

procname(pglev+1) = 'check_partition'
CALL log_timer_in()


!
!1. Open sea grid points
!------------------------
!

i_110: DO i=1,nc-1
j_110: DO j=1,nr-1
   IF (seapointglb(i,j)) THEN
      iproc = num_proc(i,j)
      IF (iproc.LE.0) THEN
         nerrs = nerrs + 1
         IF (errchk.AND.nerrs.LE.maxerrors) THEN
            WRITE (ci,'(I12)') i; ci = ADJUSTL(ci)
            WRITE (cj,'(I12)') j; cj = ADJUSTL(cj)
            WRITE (ioerr,'(A)') 'Open sea grid point ('//TRIM(ci)//','&
                             & //TRIM(cj)//') outside computational domain'
         ENDIF
      ENDIF
   ENDIF
ENDDO j_110
ENDDO i_110

!
!2. Open boundary locations
!--------------------------
!
!2.1 U-boundaries
!----------------
!

ii_210: DO ii=1,nobu
   i = iobu(ii); j = jobu(ii)
   iproc = num_proc(i,j)
   IF (iproc.LE.0) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (ci,'(I12)') i; ci = ADJUSTL(ci)
         WRITE (cj,'(I12)') j; cj = ADJUSTL(cj)
         WRITE (ioerr,'(A)') 'U-open boundary grid point ('//TRIM(ci)//','&
                          & //TRIM(cj)//') outside computational domain'
      ENDIF
   ENDIF
ENDDO ii_210

!
!2.2 V-boundaries
!----------------
!

jj_220: DO jj=1,nobv
   i = iobv(jj); j = jobv(jj)
   iproc = num_proc(i,j)
   IF (iproc.LE.0) THEN
      nerrs = nerrs + 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (ci,'(I12)') i; ci = ADJUSTL(ci)
         WRITE (cj,'(I12)') j; cj = ADJUSTL(cj)
         WRITE (ioerr,'(A)') 'V-open boundary grid point ('//TRIM(ci)//','&
                          & //TRIM(cj)//') outside computational domain'
      ENDIF
   ENDIF
ENDDO jj_220

!
!3. Exit in case of error
!------------------------
!

IF (nerrs.GT.0) THEN
   WRITE (ioerr,'(A)') 'Invalid domain decomposition'
   CALL error_abort('check_partition',ierrno_inival)
ENDIF

CALL log_timer_out()


RETURN
   
END SUBROUTINE check_partition

!========================================================================

SUBROUTINE check_statlocs(outlocs,arrname)
!************************************************************************
!
! *check_statlocs* Check station locations
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.0
!
! Description - 
!
! Reference -
!
! Module calls - error_abort, error_limits_arr_struc
!
!************************************************************************
!
USE datatypes
USE gridpars
USE paralpars
USE error_routines, ONLY: error_abort, error_limits_arr_struc
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
TYPE (StationLocs), INTENT(IN), DIMENSION(:) :: outlocs

!
!*Local variables
!
INTEGER :: istat, nostats

nostats = SIZE(outlocs)
IF (.NOT.master.OR.nostats.EQ.0) RETURN
 
procname(pglev+1) = 'check_statlocs'
CALL log_timer_in()

istat_100: DO istat=1,nostats
   CALL error_limits_arr_struc(outlocs(istat)%ipos,arrname,'%ipos',1,nc-1,1,&
                            & (/istat/))
   CALL error_limits_arr_struc(outlocs(istat)%jpos,arrname,'%jpos',1,nr-1,1,&
                            & (/istat/))
ENDDO istat_100

CALL error_abort('check_statlocs',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_statlocs

!========================================================================

SUBROUTINE check_struct_locs
!************************************************************************
!
! *check_struct_locs* Check strucure locations
!
! Author - ANTEA
!
! Version - @(COHERENS)check_model.f90  V2.6
!
! Description - 
!
! Reference -
!
! Module calls - error_abort, error_alloc, error_limits_arr
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE paralpars
USE physpars
USE structures
USE switches
USE syspars
USE error_routines, ONLY: nerrs, error_abort, error_alloc, error_limits_arr
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
CHARACTER (LEN=12) :: ci, cii, cj, cjj
INTEGER :: ii, i, j, jj


IF (.NOT.master) RETURN

procname(pglev+1) = 'check_struct_locs'
CALL log_timer_in()


!
!1. Thin dams
!------------
!

IF (iopt_thndam.EQ.1) THEN

!
!1.1 Check whether locations are inside the domain
!-------------------------------------------------
!
!  ---U-nodes
   ii_111: DO ii=1,numthinu
      i = ithinu(ii); j = jthinu(ii)
      CALL error_limits_arr(i,'ithinu',1,nc-1,1,indx=(/ii/))
      CALL error_limits_arr(j,'jthinu',1,nr-1,1,indx=(/ii/))
   ENDDO ii_111

!  ---V-nodes
   jj_112: DO jj=1,numthinv
      i = ithinv(jj); j = jthinv(jj)
      CALL error_limits_arr(i,'ithinv',1,nc-1,1,indx=(/ii/))
      CALL error_limits_arr(j,'jthinv',1,nr-1,1,indx=(/ii/))
   ENDDO jj_112

!
!1.2 Check whether surrounding cells are wet
!-------------------------------------------
!
!  ---U-nodes
   ii_121: DO ii=1,numthinu
      i = ithinu(ii); j = jthinu(ii)
      IF (.NOT.(seapointglb(i-1,j).AND.seapointglb(i,j))) THEN
         nerrs = nerrs + 1
         IF (errchk.AND.nerrs.LE.maxerrors) THEN
            WRITE (cii,'(I12)') ii; cii = ADJUSTL(cii)
            WRITE (ci,'(I12)') i; ci = ADJUSTL(ci)
            WRITE (cj,'(I12)') j; cj = ADJUSTL(cj)
            WRITE (ioerr,'(A)') 'Invalid location of U-node thin dam point '&
                       & //TRIM(cii)//' at: '//'('//TRIM(ci)//','//TRIM(cj)//')'
            WRITE (ioerr,'(A)') 'Thin dams must be surrounded by wet cells'
         ENDIF
      ENDIF
   ENDDO ii_121

!  ---V-nodes
   jj_122: DO jj=1,numthinv
      i = ithinv(jj); j = jthinv(jj)
      IF (.NOT.(seapointglb(i,j-1).AND.seapointglb(i,j))) THEN
         nerrs = nerrs + 1
         IF (errchk.AND.nerrs.LE.maxerrors) THEN
            WRITE (cjj,'(I12)') jj; cjj = ADJUSTL(cjj)
            WRITE (ci,'(I12)') i; ci = ADJUSTL(ci)
            WRITE (cj,'(I12)') j; cj = ADJUSTL(cj)
            WRITE (ioerr,'(A)') 'Invalid location of V-node thin dam point '&
                       & //TRIM(cjj)//' at: '//'('//TRIM(ci)//','//TRIM(cj)//')'
            WRITE (ioerr,'(A)') 'Thin dams must be surrounded by wet cells'
         ENDIF
      ENDIF
   ENDDO jj_122

ENDIF

!
!2. Weirs/barriers
!-----------------
!

IF (iopt_weibar.EQ.1) THEN

!
!3.1 Check whether locations are inside the domain
!-------------------------------------------------
!
!  ---U-nodes
   ii_211: DO ii=1,numwbaru
      i = iwbaru(ii); j = jwbaru(ii)
      CALL error_limits_arr(i,'iwbaru',1,nc-1,1,indx=(/ii/))
      CALL error_limits_arr(j,'jwbaru',1,nr-1,1,indx=(/ii/))
   ENDDO ii_211

!  ---V-nodes
   jj_212: DO jj=1,numwbarv
      i = iwbarv(jj); j = jwbarv(jj)
      CALL error_limits_arr(i,'iwbarv',1,nc-1,1,indx=(/ii/))
      CALL error_limits_arr(j,'jwbarv',1,nr-1,1,indx=(/ii/))
   ENDDO jj_212

!
!2.2 Check whether surrounding cells are wet
!-------------------------------------------
!
!  ---U-nodes
   ii_221: DO ii=1,numwbaru
      i = iwbaru(ii); j = jwbaru(ii)
      IF (.NOT.(seapointglb(i-1,j).AND.seapointglb(i,j))) THEN
         nerrs = nerrs + 1
         IF (errchk.AND.nerrs.LE.maxerrors) THEN
            WRITE (cii,'(I12)') ii; cii = ADJUSTL(cii)
            WRITE (ci,'(I12)') i; ci = ADJUSTL(ci)
            WRITE (cj,'(I12)') j; cj = ADJUSTL(cj)
            WRITE (ioerr,'(A)') 'Invalid location of U-node weir point '&
                       & //TRIM(cii)//' at: '//'('//TRIM(ci)//','//TRIM(cj)//')'
            WRITE (ioerr,'(A)') 'Weirs and barriers must be surrounded '//&
                              & 'by wet cells'
         ENDIF
      ENDIF
   ENDDO ii_221

!  ---V-nodes
   jj_222: DO jj=1,numwbarv
      i = iwbarv(jj); j = jwbarv(jj)
      IF (.NOT.(seapointglb(i,j-1).AND.seapointglb(i,j))) THEN
         nerrs = nerrs + 1
         IF (errchk.AND.nerrs.LE.maxerrors) THEN
            WRITE (cjj,'(I12)') jj; cjj = ADJUSTL(cjj)
            WRITE (ci,'(I12)') i; ci = ADJUSTL(ci)
            WRITE (cj,'(I12)') j; cj = ADJUSTL(cj)
            WRITE (ioerr,'(A)') 'Invalid location of V-node weir point '&
                       & //TRIM(cjj)//' at: '//'('//TRIM(ci)//','//TRIM(cj)//')'
            WRITE (ioerr,'(A)') 'Weirs and barriers must be surrounded '//&
                              & 'by wet cells'
         ENDIF
      ENDIF
   ENDDO jj_222

ENDIF

CALL error_abort('check_struct_locs',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_struct_locs

!========================================================================

SUBROUTINE check_variables(varatts,arrname)
!************************************************************************
!
! *check_variables* Check output variables
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)check_model.f90  V2.11
!
! Description - 
!
! Reference -
!
! Module calls - error_abort, error_lbound_arr_struc, error_limits_arr_struc,
!                error_vals_arr_struc_int, inquire_var
!                
!
!************************************************************************
!
USE datatypes
USE gridpars
USE paralpars
USE switches
USE error_routines, ONLY: error_abort, error_lbound_arr_struc, &
                        & error_limits_arr_struc, error_vals_arr_struc_int
USE modvars_routines, ONLY: inquire_var
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: arrname
TYPE (VariableAtts), INTENT(IN), DIMENSION(:) :: varatts

!
!*Local variables
!
INTEGER :: ivar, klev, nmax, nodim, novars, numvar, oopt
INTEGER, DIMENSION(5) :: ndims
REAL :: dep


novars = SIZE(varatts)
IF (.NOT.master.OR.novars.EQ.0) RETURN

procname(pglev+1) = 'check_variables'
CALL log_timer_in()

ivar_100: DO ivar=1,novars

!  ---variable rank
   nodim = varatts(ivar)%nrank
   CALL error_vals_arr_struc_int(nodim,arrname,'%nrank',3,1,(/0,2,3/),(/ivar/))

!  ---applied operator
   oopt = varatts(ivar)%oopt
   CALL error_vals_arr_struc_int(oopt,arrname,'%oopt',7,1,&
      & (/oopt_null,oopt_dep,oopt_int,oopt_klev,oopt_min,oopt_max,oopt_mean/),&
      & (/ivar/))
   IF (oopt.EQ.oopt_dep) THEN
      dep = varatts(ivar)%dep
      CALL error_lbound_arr_struc(dep,arrname,'%dep',0.0,.FALSE.,1,(/ivar/))
   ENDIF
   IF (oopt.EQ.oopt_klev) THEN
      klev = varatts(ivar)%klev
      CALL error_limits_arr_struc(klev,arrname,'%klev',1,nz,1,(/ivar/))
   ENDIF

!  ---variable number
   numvar = varatts(ivar)%numvar
   IF (numvar.GT.0.AND.varatts(ivar)%ivarid.GT.0) THEN
      CALL inquire_var(varatts(ivar)%ivarid,nrank=nodim)
      CALL inquire_var(varatts(ivar)%ivarid,global_dims=ndims)
      nmax = ndims(nodim)
      CALL error_limits_arr_struc(numvar,arrname,'%numvar',1,nmax,1,(/ivar/))
   ENDIF
     
ENDDO ivar_100

CALL error_abort('check_variables',ierrno_inival)

CALL log_timer_out()


RETURN

END SUBROUTINE check_variables


END MODULE check_model
