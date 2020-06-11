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

SUBROUTINE harmonic_analysis
!************************************************************************
!
! *harmonic_analysis* Perform harmonic analysis
!
! Author - Patrick Luyten
!
! Version - 4 Nov 2008  @(COHERENS)Harmonic_Analysis.f90  V2.11
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description - 
!
! Reference -
!
! Calling program - coherens_main
!
! External calls - read_cif_params
!
! Internal calls - ellips_params, harmonic_analysis_data,
!                  harmonic_analysis_grid, harmonic_analysis_init,
!                  harmonic_analysis_reset, harmonic_analysis_update
!
! Module calls - close_filepars, error_alloc
!
!************************************************************************
!
USE datatypes
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE physpars
USE switches
USE syspars
USE tide
USE timepars
USE check_model, ONLY: check_out_filepars, check_out_gpars, check_statlocs, &
                     & check_variables
USE datatypes_init, ONLY: filepars_init, outgpars_init, statlocs_init, &
                        & varatts_init
USE default_model, ONLY: default_ellvars, default_out_files, &
                       & default_out_gpars
USE error_routines, ONLY: check_mult_counter, error_abort, error_alloc, &
                        & error_alloc_struc, error_diff_vals_arrlist, &
                        & error_diff_vals_varlist, error_limits_arr, &
                        & error_vals_arr_int, warning_reset_arr_struc
USE grid_routines, ONLY: global_mask, local_proc, mask_array, Zcoord_var
USE inout_paral, ONLY: combine_write_stats_loc, combine_write_submod
USE inout_routines, ONLY: close_filepars, inout_atts_out, write_time, &
                        & write_vars
USE math_library, ONLY: complex_polar
USE model_output, ONLY: define_out0d_vals, define_out2d_vals, define_out3d_vals
USE nla_library, ONLY: cholesky_decomp, cholesky_solve
USE paral_comms, ONLY: combine_mod
USE paral_utilities, ONLY: sum_vars
USE reset_model, ONLY: reset_out_files, reset_out_gpars, reset_out_stats, &
                     & reset_out_vars
USE time_routines, ONLY: add_secs_to_date, add_secs_to_phase, check_date, &
                       & date_to_year, diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index, mult_index

IMPLICIT none

!
!* Local variables
!
LOGICAL, SAVE :: info = .FALSE., last_call
LOGICAL :: astro_ref, first, flag, fmask, gridded, newvals, oldvals, packing, second
CHARACTER (LEN=12) :: cset
CHARACTER (LEN=lentime) :: ctime, refdate
INTEGER :: i, ierr, ifreq, ii, iifreq, iivar, imax, iproc, iset, istat, itype, ivar, &
         & i0, i1, i2, i3, j, jj, jmax, j1, j2, j3, k, kk, l, lbufa, lbufa1, &
         & lbufa2, lbufb, lbufb1, lbufb2, lbuf0da, lbuf0db, lbuf2da, lbuf2db, &
         & lbuf3da, lbuf3db, lloc, lnew, lold, l0data, l0datb, l1, l2, &
         & l2data, l2datb, l3data, l3datb, m, maxstatsloc, maxvars, mm, n, &
         & nanalvars0d, nanalvars2d, nanalvars3d, ncout, nhlocout, nn, &
         & nocoords0d, nocoords2d, nocoords3d, nodat, nodim, nofreqs, nosecs2, &
         & nostats, novars, nowetlocout, nowetout, npcc, nrank, nrout, nzout, &
         & n1, n2, n3, time_format, vcoord
INTEGER, SAVE ::  maxstats, novals0da_new, novals0da_old, novals0db_new, &
                & novals0db_old, novals2da_new, novals2da_old, novals2db_new, &
                & novals2db_old, novals3da_new, novals3da_old, novals3db_new, &
                & novals3db_old, novars0d, novars2d, novars3d
INTEGER (KIND=kndilong) :: nosecs
REAL, PARAMETER :: epstol = 1.0E-20 
REAL :: amin, aplus, cfmin, cfplus, deltanal, fill_value, phi, rmin, rmini, &
      & rminr, rplus, rplusi, rplusr, sig, theta, xarg
REAL (KIND=kndrlong) :: rtime
TYPE (FileParams) :: filepars

INTEGER, DIMENSION(3) :: tlims, xlims, xlimsglb, ylims, ylimsglb, zlims
INTEGER, DIMENSION(4) :: ndims
INTEGER, DIMENSION(0:nofreqsanal) :: indexarr
INTEGER, DIMENSION(novarsanal) :: vecids
INTEGER, DIMENSION(2,2,nprocs) :: lprocs
REAL, DIMENSION(7) :: outell
REAL, DIMENSION(0:nofreqsanal) :: fcos
REAL, DIMENSION(nofreqsanal) :: fnode_anal, fsin, phase_anal, phase_ref
REAL, DIMENSION(0:nofreqsanal) :: xvec
REAL, DIMENSION(novarsanal) :: outdat, outsub

LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskglb, maskvals
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: icount, lmask, nanaltot, nclocout, &
                               & ndat0da, ndat0db, ndat2da, ndat2db, ndat3da, &
                               & ndat3db, nostatsloc, nrlocout, ntout
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:) :: istatpos, ivarsanal0d, &
                               & ivarsanal2d, ivarsanal3d, ivarsell2d, &
                               & ivarsell3d, jstatpos, nostatsprocs
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: limloc, lstatsprocs
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: limprocs
REAL, ALLOCATABLE, DIMENSION(:) :: buffold
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: buff0da, buff0db, buff2da, buff2db, &
                               & buff3da, buff3db, gsigout, out1damp, &
                               & out1dpha, out1dres
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: depout, ell2dsub, out2damp, &
                               & out2dpha, out2dres, xout, yout
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: acoef, bcoef, ell3dsub, gsout, &
                               & out3damp, out3dpha, out3dres, zout
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: ell4dsub, out4damp, out4dpha, &
                                             & out4dres
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: outvars

!
! Name           Type     Purpose
!------------------------------------------------------------------------------
!*last_call*     LOGICAL .TRUE. if the routine has been called for the last
!                        time
!*icount*        INTEGER Counter increased after each update time
!*istatpos*      INTEGER X-index of local output station
!*ivarsanal0d*   INTEGER Array of indices (per file) of 0-D output variables
!*ivarsanal2d*   INTEGER Array of indices (per file) of 2-D output variables
!*ivarsanal3d*   INTEGER Array of indices (per file) of 3-D output variables
!*ivarsell2d*    INTEGER Array of indices (per file) of 2-D elliptic parameters
!*ivarsell3d*    INTEGER Array of indices (per file) of 3-D elliptic parameters
!*jstatpos*      INTEGER Y-index of local output station
!*limloc*        INTEGER Start/End/Increment values of local (regular) output
!                        grid in X- and Y-directions
!*limprocs*      INTEGER Start/End indices of local array section within global
!                        array (regular grid)
!*lstatsprocs*   INTEGER Output labels of local stations for each process
!                        domain
!*nanaltot*      INTEGER Number of time steps within one analysis period (+1)
!*nclocout*      INTEGER X-dimension of local (regular) output grid
!*ndat0da*       INTEGER Total number (per file) of 0-D data in 'a'-buffer
!*ndat0db*       INTEGER Total number (per file) of 0-D data in 'b'-buffer
!*ndat2da*       INTEGER Total number (per file) of 2-D data in 'a'-buffer
!*ndat2db*       INTEGER Total number (per file) of 2-D data in 'b'-buffer
!*ndat3da*       INTEGER Total number (per file) of 3-D data in 'a'-buffer
!*ndat3db*       INTEGER Total number (per file) of 3-D data in 'b'-buffer
!*nostatsloc*    INTEGER Number of local output stations (per file)
!*nostatsprocs*  INTEGER Number of local output stations for each process
!                        domain
!*novals0da_new* INTEGER Total number (over all files) of data in 0-D
!                        'a'-buffer at new time
!*novals0da_old* INTEGER Total number (over all files) of data in 0-D
!                        'a'-buffer at old time
!*novals0db_new* INTEGER Total number (over all files) of data in 0-D
!                        'b'-buffer at new time
!*novals0db_old* INTEGER Total number (over all files) of data in 0-D
!                        'b'-buffer at old time
!*novals2da_new* INTEGER Total number (over all files) of data in 2-D
!                        'a'-buffer at new time
!*novals2da_old* INTEGER Total number (over all files) of data in 2-D
!                        'a'-buffer at old time
!*novals2db_new* INTEGER Total number (over all files) of data in 2-D
!                        'b'-buffer at new time
!*novals2db_old* INTEGER Total number (over all files) of data in 2-D
!                        'b'-buffer at old time
!*novals3da_new* INTEGER Total number (over all files) of data in 3-D
!                        'a'-buffer at new time
!*novals3da_old* INTEGER Total number (over all files) of data in 3-D
!                        'a'-buffer at old time
!*novals3db_new* INTEGER Total number (over all files) of data in 3-D
!                        'b'-buffer at new time
!*novals3db_old* INTEGER Total number (over all files) of data in 3-D
!                        'b'-buffer at old time
!*novars0d*      INTEGER Total number of 0-D output variables
!*novars2d*      INTEGER Total number of 2-D output variables
!*novars3d*      INTEGER Total number of 3-D output variables
!*nrlocout*      INTEGER Y-dimension of local (regular) output grid
!*ntout*         INTEGER Output time counter
!*acoef*         REAL    Matrix of linear system to be solved
!*bcoef*         REAL    R.h.s. of linear system to be solved
!*buff0da*       REAL    Buffer for temporary storage of 0-D 'a'-data
!*buff0db*       REAL    Buffer for temporary storage of 0-D 'b'-data
!*buff2da*       REAL    Buffer for temporary storage of 2-D 'a'-data
!*buff2db*       REAL    Buffer for temporary storage of 2-D 'b'-data
!*buff3da*       REAL    Buffer for temporary storage of 3-D 'a'-data
!*buff3db*       REAL    Buffer for temporary storage of 3-D 'b'-data
!*depout*        REAL    Output depth array                                 [m]
!*fnode_anal*    REAL    Nodel factors for analysed frequencies
!*phase_anal*    REAL    Astronomical argument with respect to central time
!                                                                      [radian]
!*phase_ref*     REAL    Astronomical argument with respect to reference date
!                                                                      [radian]
!*xout*          REAL    Output X-coordinates                [degree_east or m]
!*yout*          REAL    Output Y-coordinates               [degree_north or m]
!
!------------------------------------------------------------------------------
!

IF (nt.EQ.0) last_call=.FALSE.
IF (last_call) RETURN

procname(pglev+1) = 'harmonic_analysis'
CALL log_timer_in(npcc)

!
!1. Initialise on first call
!---------------------------
!

IF (nt.EQ.0) CALL harmonic_analysis_init

!
!2. Write output grid
!--------------------
!

CALL harmonic_analysis_grid

!
!3. Create buffers for data storage
!----------------------------------
!
!3.1 Buffer sizes
!----------------
!

ndat0da = 0; ndat2da = 0; ndat3da = 0
ndat0db = 0; ndat2db = 0; ndat3db = 0

iset_310: DO iset=1,nosetsanal
   tlims = analgpars(iset)%tlims
   IF ((nt.GE.tlims(1)).AND.(nt.LE.tlims(2))) THEN
      nofreqs = nofreqsharm(iset)
      nanalvars0d = res0d(iset)%novars
      nanalvars2d = res2d(iset)%novars
      nanalvars3d = res3d(iset)%novars
      nhlocout = nclocout(iset)*nrlocout(iset)
      nzout = analgpars(iset)%nzout
      ndat0da(iset) = ndat0da(iset) + nanalvars0d*(nofreqs+1)
      ndat0db(iset) = ndat0db(iset) + nanalvars0d*nofreqs
      IF (analgpars(iset)%gridded) THEN
         ndat2da(iset) = ndat2da(iset) + nhlocout*nanalvars2d*(nofreqs+1)
         ndat2db(iset) = ndat2db(iset) + nhlocout*nanalvars2d*nofreqs
         ndat3da(iset) = ndat3da(iset) + nhlocout*nzout*nanalvars3d*(nofreqs+1)
         ndat3db(iset) = ndat3db(iset) + nhlocout*nzout*nanalvars3d*nofreqs
      ELSE
         ndat2da(iset) = ndat2da(iset) + nostatsloc(iset)*nanalvars2d*&
                                       & (nofreqs+1)
         ndat2db(iset) = ndat2db(iset) + nostatsloc(iset)*nanalvars2d*nofreqs
         ndat3da(iset) = ndat3da(iset) + nostatsloc(iset)*nzout*nanalvars3d*&
                                       & (nofreqs+1)
         ndat3db(iset) = ndat3db(iset) + nostatsloc(iset)*nzout*nanalvars3d*&
                                       & nofreqs
      ENDIF
   ENDIF
ENDDO iset_310

!
!3.2 Refresh buffers (if necessary)
!----------------------------------
!
!3.2.1 0-D buffer
!----------------
!
!---'a'-buffer
novals0da_new = SUM(ndat0da)
IF ((novals0da_old.EQ.0).AND.(novals0da_new.GT.0)) THEN
   ALLOCATE (buff0da(novals0da_new),STAT=errstat)
   CALL error_alloc('buff0da',1,(/novals0da_new/),kndrtype)
ELSEIF ((novals0da_old.GT.0).AND.(novals0da_new.EQ.0)) THEN
   DEALLOCATE (buff0da)
ELSEIF (novals0da_old.NE.novals0da_new) THEN
   ALLOCATE (buffold(novals0da_old),STAT=errstat)
   CALL error_alloc('buffold',1,(/novals0da_old/),kndrtype)
   buffold = buff0da
   DEALLOCATE (buff0da)
   ALLOCATE (buff0da(novals0da_new),STAT=errstat)
   CALL error_alloc('buff0da',1,(/novals0da_new/),kndrtype)

   lold = 0; lnew = 0
   iset_3211: DO iset=1,nosetsanal
      tlims = analgpars(iset)%tlims
      IF (res0d(iset)%novars.GT.0) THEN
         nodat = ndat0da(iset)
         oldvals = ((nt-ic3d).GE.tlims(1)).AND.((nt-ic3d).LE.tlims(2))
         newvals = (nt.GE.tlims(1)).AND.(nt.LE.tlims(2))
         IF (oldvals.AND.newvals) THEN
            buff0da(lnew+1:lnew+nodat) = buffold(lold+1:lold+nodat)
            lold = lold + nodat; lnew = lnew + nodat 
         ELSEIF (oldvals.AND.(.NOT.newvals)) THEN
            lold = lold + nodat
         ELSEIF ((.NOT.oldvals).AND.newvals) THEN
            lnew = lnew + nodat
         ENDIF
      ENDIF
   ENDDO iset_3211

   DEALLOCATE (buffold)

ENDIF

novals0da_old = novals0da_new

!---'b'-buffer
novals0db_new = SUM(ndat0db)
IF ((novals0db_old.EQ.0).AND.(novals0db_new.GT.0)) THEN
   ALLOCATE (buff0db(novals0db_new),STAT=errstat)
   CALL error_alloc('buff0db',1,(/novals0db_new/),kndrtype)
ELSEIF ((novals0db_old.GT.0).AND.(novals0db_new.EQ.0)) THEN
   DEALLOCATE (buff0db)
ELSEIF (novals0db_old.NE.novals0db_new) THEN
   ALLOCATE (buffold(novals0db_old),STAT=errstat)
   CALL error_alloc('buffold',1,(/novals0db_old/),kndrtype)
   buffold = buff0db
   DEALLOCATE (buff0db)
   ALLOCATE (buff0db(novals0db_new),STAT=errstat)
   CALL error_alloc('buff0db',1,(/novals0db_new/),kndrtype)

   lold = 0; lnew = 0
   iset_3212: DO iset=1,nosetsanal
      tlims = analgpars(iset)%tlims
      IF (res0d(iset)%novars.GT.0) THEN
         nodat = ndat0db(iset)
         oldvals = ((nt-ic3d).GE.tlims(1)).AND.((nt-ic3d).LE.tlims(2))
         newvals = (nt.GE.tlims(1)).AND.(nt.LE.tlims(2))
         IF (oldvals.AND.newvals) THEN
            buff0db(lnew+1:lnew+nodat) = buffold(lold+1:lold+nodat)
            lold = lold + nodat; lnew = lnew + nodat 
         ELSEIF (oldvals.AND.(.NOT.newvals)) THEN
            lold = lold + nodat
         ELSEIF ((.NOT.oldvals).AND.newvals) THEN
            lnew = lnew + nodat
         ENDIF
      ENDIF
   ENDDO iset_3212

   DEALLOCATE (buffold)

ENDIF

novals0db_old = novals0db_new

!
!3.2.2 2-D buffer
!----------------
!
!---'a'-buffer
novals2da_new = SUM(ndat2da)
IF ((novals2da_old.EQ.0).AND.(novals2da_new.GT.0)) THEN
   ALLOCATE (buff2da(novals2da_new),STAT=errstat)
   CALL error_alloc('buff2da',1,(/novals2da_new/),kndrtype)
ELSEIF ((novals2da_old.GT.0).AND.(novals2da_new.EQ.0)) THEN
   DEALLOCATE (buff2da)
ELSEIF (novals2da_old.NE.novals2da_new) THEN
   ALLOCATE (buffold(novals2da_old),STAT=errstat)
   CALL error_alloc('buffold',1,(/novals2da_old/),kndrtype)
   buffold = buff2da
   DEALLOCATE (buff2da)
   ALLOCATE (buff2da(novals2da_new),STAT=errstat)
   CALL error_alloc('buff2da',1,(/novals2da_new/),kndrtype)

   lold = 0; lnew = 0
   iset_3221: DO iset=1,nosetsanal
      tlims = analgpars(iset)%tlims
      IF (res2d(iset)%novars.GT.0) THEN
         nodat = ndat2da(iset)
         oldvals = ((nt-ic3d).GE.tlims(1)).AND.((nt-ic3d).LE.tlims(2))
         newvals = (nt.GE.tlims(1)).AND.(nt.LE.tlims(2))
         IF (oldvals.AND.newvals) THEN
            buff2da(lnew+1:lnew+nodat) = buffold(lold+1:lold+nodat)
            lold = lold + nodat; lnew = lnew + nodat 
         ELSEIF (oldvals.AND.(.NOT.newvals)) THEN
            lold = lold + nodat
         ELSEIF ((.NOT.oldvals).AND.newvals) THEN
            lnew = lnew + nodat
         ENDIF
      ENDIF
   ENDDO iset_3221

   DEALLOCATE (buffold)

ENDIF

novals2da_old = novals2da_new

!---'b'-buffer
novals2db_new = SUM(ndat2db)
IF ((novals2db_old.EQ.0).AND.(novals2db_new.GT.0)) THEN
   ALLOCATE (buff2db(novals2db_new),STAT=errstat)
   CALL error_alloc('buff2db',1,(/novals2db_new/),kndrtype)
ELSEIF ((novals2db_old.GT.0).AND.(novals2db_new.EQ.0)) THEN
   DEALLOCATE (buff2db)
ELSEIF (novals2db_old.NE.novals2db_new) THEN
   ALLOCATE (buffold(novals2db_old),STAT=errstat)
   CALL error_alloc('buffold',1,(/novals2db_old/),kndrtype)
   buffold = buff2db
   DEALLOCATE (buff2db)
   ALLOCATE (buff2db(novals2db_new),STAT=errstat)
   CALL error_alloc('buff2db',1,(/novals2db_new/),kndrtype)

   lold = 0; lnew = 0
   iset_3222: DO iset=1,nosetsanal
      tlims = analgpars(iset)%tlims
      IF (res2d(iset)%novars.GT.0) THEN
         nodat = ndat2db(iset)
         oldvals = ((nt-ic3d).GE.tlims(1)).AND.((nt-ic3d).LE.tlims(2))
         newvals = (nt.GE.tlims(1)).AND.(nt.LE.tlims(2))
         IF (oldvals.AND.newvals) THEN
            buff2db(lnew+1:lnew+nodat) = buffold(lold+1:lold+nodat)
            lold = lold + nodat; lnew = lnew + nodat 
         ELSEIF (oldvals.AND.(.NOT.newvals)) THEN
            lold = lold + nodat
         ELSEIF ((.NOT.oldvals).AND.newvals) THEN
            lnew = lnew + nodat
         ENDIF
      ENDIF
   ENDDO iset_3222

   DEALLOCATE (buffold)

ENDIF

novals2db_old = novals2db_new

!
!3.2.3 3-D buffer
!----------------
!
!---'a'-buffer
novals3da_new = SUM(ndat3da)
IF ((novals3da_old.EQ.0).AND.(novals3da_new.GT.0)) THEN
   ALLOCATE (buff3da(novals3da_new),STAT=errstat)
   CALL error_alloc('buff3da',1,(/novals3da_new/),kndrtype)
ELSEIF ((novals3da_old.GT.0).AND.(novals3da_new.EQ.0)) THEN
   DEALLOCATE (buff3da)
ELSEIF (novals3da_old.NE.novals3da_new) THEN
   ALLOCATE (buffold(novals3da_old),STAT=errstat)
   CALL error_alloc('buffold',1,(/novals3da_old/),kndrtype)
   buffold = buff3da
   DEALLOCATE (buff3da)
   ALLOCATE (buff3da(novals3da_new),STAT=errstat)
   CALL error_alloc('buff3da',1,(/novals3da_new/),kndrtype)

   lold = 0; lnew = 0
   iset_3231: DO iset=1,nosetsanal
      tlims = analgpars(iset)%tlims
      IF (res3d(iset)%novars.GT.0) THEN
         nodat = ndat3da(iset)
         oldvals = ((nt-ic3d).GE.tlims(1)).AND.((nt-ic3d).LE.tlims(2))
         newvals = (nt.GE.tlims(1)).AND.(nt.LE.tlims(2))
         IF (oldvals.AND.newvals) THEN
            buff3da(lnew+1:lnew+nodat) = buffold(lold+1:lold+nodat)
            lold = lold + nodat; lnew = lnew + nodat 
         ELSEIF (oldvals.AND.(.NOT.newvals)) THEN
            lold = lold + nodat
         ELSEIF ((.NOT.oldvals).AND.newvals) THEN
            lnew = lnew + nodat
         ENDIF
      ENDIF
   ENDDO iset_3231

   DEALLOCATE (buffold)

ENDIF

novals3da_old = novals3da_new

!---'b'-buffer
novals3db_new = SUM(ndat3db)
IF ((novals3db_old.EQ.0).AND.(novals3db_new.GT.0)) THEN
   ALLOCATE (buff3db(novals3db_new),STAT=errstat)
   CALL error_alloc('buff3db',1,(/novals3db_new/),kndrtype)
ELSEIF ((novals3db_old.GT.0).AND.(novals3db_new.EQ.0)) THEN
   DEALLOCATE (buff3db)
ELSEIF (novals3db_old.NE.novals3db_new) THEN
   ALLOCATE (buffold(novals3db_old),STAT=errstat)
   CALL error_alloc('buffold',1,(/novals3db_old/),kndrtype)
   buffold = buff3db
   DEALLOCATE (buff3db)
   ALLOCATE (buff3db(novals3db_new),STAT=errstat)
   CALL error_alloc('buff3db',1,(/novals3db_new/),kndrtype)

   lold = 0; lnew = 0
   iset_3232: DO iset=1,nosetsanal
      tlims = analgpars(iset)%tlims
      IF (res3d(iset)%novars.GT.0) THEN
         nodat = ndat3db(iset)
         oldvals = ((nt-ic3d).GE.tlims(1)).AND.((nt-ic3d).LE.tlims(2))
         newvals = (nt.GE.tlims(1)).AND.(nt.LE.tlims(2))
         IF (oldvals.AND.newvals) THEN
            buff3db(lnew+1:lnew+nodat) = buffold(lold+1:lold+nodat)
            lold = lold + nodat; lnew = lnew + nodat 
         ELSEIF (oldvals.AND.(.NOT.newvals)) THEN
            lold = lold + nodat
         ELSEIF ((.NOT.oldvals).AND.newvals) THEN
            lnew = lnew + nodat
         ENDIF
      ENDIF
   ENDDO iset_3232

   DEALLOCATE (buffold)

ENDIF

novals3db_old = novals3db_new

!
!4. Harmonic analysis and output
!-------------------------------
!
!4.1 Initialise buffers
!----------------------
!

first = .TRUE.; second = .FALSE.
CALL harmonic_analysis_reset

!
!4.2 Update buffers
!------------------
!

CALL harmonic_analysis_update

!
!4.3 Write harmonic output
!-------------------------
!

CALL harmonic_analysis_data

!
!4.4 Reinitialise/update buffers
!-------------------------------
!

first = .FALSE.; second = .TRUE.
CALL harmonic_analysis_reset
CALL harmonic_analysis_update

!
!5. Finalise after last output
!-----------------------------
!
!5.1 Close data files
!--------------------
!

last_call = .TRUE.
iset_510: DO iset=1,nosetsanal
   IF (.NOT.cold_start.AND.nt.LT.analgpars(iset)%tlims(2)) last_call = .FALSE.
   IF (master.AND.(cold_start.OR.nt.EQ.analgpars(iset)%tlims(2))) THEN
      nofreqs = nofreqsharm(iset)

!     ---0-D case
      IF (res0d(iset)%defined) CALL close_filepars(res0d(iset))
      ifreq_511: DO ifreq=1,nofreqs
         IF (amp0d(iset,ifreq)%defined) CALL close_filepars(amp0d(iset,ifreq))
         IF (pha0d(iset,ifreq)%defined) CALL close_filepars(pha0d(iset,ifreq))
      ENDDO ifreq_511

!     ---2-D case
      IF (res2d(iset)%defined) CALL close_filepars(res2d(iset))
      ifreq_512: DO ifreq=1,nofreqs
         IF (amp2d(iset,ifreq)%defined) CALL close_filepars(amp2d(iset,ifreq))
         IF (pha2d(iset,ifreq)%defined) CALL close_filepars(pha2d(iset,ifreq))
         IF (ell2d(iset,ifreq)%defined) CALL close_filepars(ell2d(iset,ifreq))
      ENDDO ifreq_512

!     ---3-D case
      IF (res3d(iset)%defined) CALL close_filepars(res3d(iset))
      ifreq_513: DO ifreq=1,nofreqs
         IF (amp3d(iset,ifreq)%defined) CALL close_filepars(amp3d(iset,ifreq))
         IF (pha3d(iset,ifreq)%defined) CALL close_filepars(pha3d(iset,ifreq))
         IF (ell3d(iset,ifreq)%defined) CALL close_filepars(ell3d(iset,ifreq))
      ENDDO ifreq_513

   ENDIF
ENDDO iset_510

!
!5.2 Deallocate
!--------------
!

IF (last_call) THEN

   DEALLOCATE (ivarsanal0d,ivarsanal2d,ivarsanal3d,ivarsell2d,ivarsell3d,&
             & outvars)
   DEALLOCATE (nclocout,nrlocout,limprocs,limloc)
   DEALLOCATE (nostatsprocs,nostatsloc,lstatsprocs,istatpos,jstatpos)
   DEALLOCATE (ndat0da,ndat2da,ndat3da,ndat0db,ndat2db,ndat3db)
   DEALLOCATE (icount,nanaltot,ntout)
   DEALLOCATE (acoef,bcoef)

   IF (ALLOCATED(buff0da)) DEALLOCATE (buff0da)
   IF (ALLOCATED(buff2da)) DEALLOCATE (buff2da)
   IF (ALLOCATED(buff3da)) DEALLOCATE (buff3da)
   IF (ALLOCATED(buff0db)) DEALLOCATE (buff0db)
   IF (ALLOCATED(buff2db)) DEALLOCATE (buff2db)
   IF (ALLOCATED(buff3db)) DEALLOCATE (buff3db)

ENDIF

CALL log_timer_out(npcc,itm_user)


RETURN

CONTAINS

!========================================================================

SUBROUTINE harmonic_analysis_init
!************************************************************************
!
! *harmonic_analysis_init* Initialise harmonic analysis and output parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Harmonic_Analysis.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - harmonic_analysis
!
! External calls - read_cif_params, usrdef_anal_freqs, usrdef_anal_params
!
! Module calls - check_date, check_out_filepars, check_out_gpars,
!   check_statlocs, check_variables, cholesky_decomp, default_ellvars,
!   default_out_files, default_out_gpars, error_abort, error_alloc,
!   error_alloc_struc, error_diff_vals_arrlist, error_diff_vals_varlist,
!   error_limits_arr, error_vals_arr_int, filepars_init, inout_atts_out,
!   lim_dims, local_proc, mask_array, reset_out_files, outgpars_init,
!   reset_out_gpars, reset_out_stats, reset_out_vars, statlocs_init, sum_vars,
!   varatts_init, warning_reset_arr_struc
!
!************************************************************************
!


procname(pglev+1) = 'harmonic_analysis_init'
CALL log_timer_in()

!
!1. Harmonic parameters
!----------------------
!
!1.1 Allocate and initialise
!---------------------------
!
!---frequencies
ALLOCATE (nofreqsharm(nosetsanal),STAT=errstat)
CALL error_alloc('nofreqsharm',1,(/nosetsanal/),kndint)
nofreqsharm = 0
ALLOCATE (ifreqsharm(nosetsanal,nofreqsanal),STAT=errstat)
CALL error_alloc('ifreqsharm',2,(/nosetsanal,nofreqsanal/),kndint)
ifreqsharm = 0

!---parameters for analysis
ALLOCATE (index_anal(nofreqsanal),STAT=errstat)
CALL error_alloc('index_anal',1,(/nofreqsanal/),kndint)
index_anal = 0
ALLOCATE (icanal(nosetsanal),STAT=errstat)
CALL error_alloc('icanal',1,(/nosetsanal/),kndint)
icanal = ic3d
ALLOCATE (cdate_time_ref(nosetsanal),STAT=errstat)
CALL error_alloc('cdate_time_ref',1,(/nosetsanal/),kndchar,lenstr=lentime)
cdate_time_ref = cdatetime_undef
ALLOCATE (harm_freq_names(nofreqsanal),STAT=errstat)
CALL error_alloc('harm_freq_names',1,(/nofreqsanal/),kndchar,lenstr=lenname)
harm_freq_names = ''
ALLOCATE (harm_freq(nofreqsanal),STAT=errstat)
CALL error_alloc('harm_freq',1,(/nofreqsanal/),kndrtype)
harm_freq = 0.0

!
!1.2 Define
!----------
!

IF (ciffile%status.EQ.'R') THEN
   CALL read_cif_params(icif_harm)
ELSE
   CALL usrdef_anal_freqs
ENDIF

!
!1.3 Reset
!---------
!

ifreq_130: DO ifreq=1,nofreqsanal
   CALL error_limits_arr(index_anal(ifreq),'index_anal',0,MaxConstituents,1,&
                       & indx=(/ifreq/))
   CALL error_diff_vals_varlist(index_anal,'index_anal',.TRUE.)
   IF (index_anal(ifreq).GT.0) THEN
      iifreq = index_anal(ifreq)
      harm_freq(ifreq) = tidal_spectrum(iifreq)
      harm_freq_names(ifreq) = tidal_freq_names(iifreq)
   ENDIF
ENDDO ifreq_130

!
!1.4 Check
!---------
!
!1.4.1 Frequencies
!-----------------
!

iset_141: DO iset=1,nosetsanal
   nofreqs = nofreqsharm(iset)
   CALL error_limits_arr(nofreqs,'nofreqsharm',1,nofreqsanal,1,(/iset/))
   ifreq_1411: DO ifreq=1,nofreqs
      CALL error_limits_arr(ifreqsharm(iset,ifreq),'ifreqsharm',1,&
                          & nofreqsanal,2,(/iset,ifreq/))
   ENDDO ifreq_1411
   CALL error_diff_vals_arrlist(ifreqsharm(iset,1:nofreqs),'ifreqsharm',2,2,&
                             & (/iset,1/))
ENDDO iset_141

!
!1.4.2 Reference dates
!---------------------
!

iset_142: DO iset=1,nosetsanal
   IF (cdate_time_ref(iset).NE.cdatetime_undef) THEN
      WRITE (cset,'(I12)') iset; cset = ADJUSTL(cset)
      CALL check_date(cdate_time_ref(iset),&
                   & 'cdate_time_ref'//'('//TRIM(cset)//')')
   ENDIF
ENDDO iset_142

!
!2. Parameters for harmonic output
!---------------------------------
!
!2.1 Allocate and initialise
!---------------------------
!
!---variables
ALLOCATE (analvars(novarsanal),STAT=errstat)
CALL error_alloc_struc('analvars',1,(/novarsanal/),'VariableAtts')
CALL varatts_init(analvars)
ALLOCATE (ellvars(14),STAT=errstat)
CALL error_alloc_struc('ellvars',1,(/14/),'VariableAtts')
CALL varatts_init(ellvars)

!---file parameters
ALLOCATE (res0d(nosetsanal),STAT=errstat)
CALL error_alloc_struc('res0d',1,(/nosetsanal/),'FileParams')
CALL filepars_init(res0d)
ALLOCATE (res2d(nosetsanal),STAT=errstat)
CALL error_alloc_struc('res2d',1,(/nosetsanal/),'FileParams')
CALL filepars_init(res2d)
ALLOCATE (res3d(nosetsanal),STAT=errstat)
CALL error_alloc_struc('res3d',1,(/nosetsanal/),'FileParams')
CALL filepars_init(res3d)
ALLOCATE (amp0d(nosetsanal,nofreqsanal),STAT=errstat)
CALL error_alloc_struc('amp0d',2,(/nosetsanal,nofreqsanal/),'FileParams')
CALL filepars_init(amp0d)
ALLOCATE (amp2d(nosetsanal,nofreqsanal),STAT=errstat)
CALL error_alloc_struc('amp2d',2,(/nosetsanal,nofreqsanal/),'FileParams')
CALL filepars_init(amp2d)
ALLOCATE (amp3d(nosetsanal,nofreqsanal),STAT=errstat)
CALL error_alloc_struc('amp3d',2,(/nosetsanal,nofreqsanal/),'FileParams')
CALL filepars_init(amp3d)
ALLOCATE (pha0d(nosetsanal,nofreqsanal),STAT=errstat)
CALL error_alloc_struc('pha0d',2,(/nosetsanal,nofreqsanal/),'FileParams')
CALL filepars_init(pha0d)
ALLOCATE (pha2d(nosetsanal,nofreqsanal),STAT=errstat)
CALL error_alloc_struc('pha2d',2,(/nosetsanal,nofreqsanal/),'FileParams')
CALL filepars_init(pha2d)
ALLOCATE (pha3d(nosetsanal,nofreqsanal),STAT=errstat)
CALL error_alloc_struc('pha3d',2,(/nosetsanal,nofreqsanal/),'FileParams')
CALL filepars_init(pha3d)
ALLOCATE (ell2d(nosetsanal,nofreqsanal),STAT=errstat)
CALL error_alloc_struc('ell2d',2,(/nosetsanal,nofreqsanal/),'FileParams')
CALL filepars_init(ell2d)
ALLOCATE (ell3d(nosetsanal,nofreqsanal),STAT=errstat)
CALL error_alloc_struc('ell3d',2,(/nosetsanal,nofreqsanal/),'FileParams')
CALL filepars_init(ell3d)

!---output specifiers
ALLOCATE (analgpars(nosetsanal),STAT=errstat)
CALL error_alloc_struc('analgpars',1,(/nosetsanal/),'OutFileParams')
CALL outgpars_init(analgpars)

!---output stations
ALLOCATE (lstatsanal(nosetsanal,nostatsanal),STAT=errstat)
CALL error_alloc('lstatsanal',2,(/nosetsanal,nostatsanal/),kndint)
lstatsanal = 0
ALLOCATE (analstatlocs(nostatsanal),STAT=errstat)
CALL error_alloc_struc('analstatloc',1,(/nostatsanal/),'StationLocs')
CALL statlocs_init(analstatlocs)

!---variable indices
ALLOCATE (ivarsanal(nosetsanal,novarsanal),STAT=errstat)
CALL error_alloc('ivarsanal',2,(/nosetsanal,novarsanal/),kndint)
ivarsanal = 0

!---elliptic parameters
ALLOCATE (ivarsell(nosetsanal,14),STAT=errstat)
CALL error_alloc('ivarsell',2,(/nosetsanal,14/),kndint)
ivarsell = 0
ALLOCATE (ivecell2d(nosetsanal,2),STAT=errstat)
CALL error_alloc('ivecell2d',2,(/nosetsanal,2/),kndint)
ivecell2d = 0
ALLOCATE (ivecell3d(nosetsanal,2),STAT=errstat)
CALL error_alloc('ivecell3d',2,(/nosetsanal,2/),kndint)
ivecell3d = 0

!
!2.2 Defaults
!------------
!
!---variables and vectors
CALL default_ellvars(ellvars)

!---file parameters
CALL default_out_files(res0d,'HR')
CALL default_out_files(res2d,'HR')
CALL default_out_files(res3d,'HR')
CALL default_out_files(amp0d,'HA')
CALL default_out_files(amp2d,'HA')
CALL default_out_files(amp3d,'HA')
CALL default_out_files(pha0d,'HP')
CALL default_out_files(pha2d,'HP')
CALL default_out_files(pha3d,'HP')
CALL default_out_files(ell2d,'HE')
CALL default_out_files(ell3d,'HE')

!---output parameters
CALL default_out_gpars(analgpars)

!
!2.3 Define
!---------
!

IF (ciffile%status.EQ.'R') THEN
   CALL read_cif_params(icif_anal)
ELSE
   CALL usrdef_anal_params
ENDIF

!
!2.4 Reset
!---------
!
!---file parameters
CALL reset_out_files(res0d,'res0d',analgpars,.FALSE.)
CALL reset_out_files(res2d,'res2d',analgpars,.FALSE.)
CALL reset_out_files(res3d,'res3d',analgpars,.FALSE.)
CALL reset_out_files(amp0d,'amp0d',analgpars,.FALSE.)
CALL reset_out_files(amp2d,'amp2d',analgpars,.FALSE.)
CALL reset_out_files(amp3d,'amp3d',analgpars,.FALSE.)
CALL reset_out_files(pha0d,'pha0d',analgpars,.FALSE.)
CALL reset_out_files(pha2d,'pha2d',analgpars,.FALSE.)
CALL reset_out_files(pha3d,'pha3d',analgpars,.FALSE.)
CALL reset_out_files(ell2d,'ell2d',analgpars,.FALSE.)
CALL reset_out_files(ell3d,'ell3d',analgpars,.FALSE.)

!---variable atrributes
CALL reset_out_vars(analvars)

!---output grid dimension
iset_241: DO iset=1,nosetsanal
   nodim = 0
   ivar_2411: DO ivar=1,novarsanal
      iivar = ivarsanal(iset,ivar)
      IF (iivar.GT.0) THEN
         IF (analvars(iivar)%nrank.EQ.3) THEN
            nodim = 3; EXIT ivar_2411
         ENDIF
      ENDIF
   ENDDO ivar_2411
   IF (nodim.EQ.0) THEN
      ivar_2412: DO ivar=1,novarsanal
         iivar = ivarsanal(iset,ivar)
         IF (iivar.GT.0) THEN
            IF (analvars(iivar)%nrank.EQ.2) THEN
               nodim = 2; EXIT ivar_2412
            ENDIF
         ENDIF
      ENDDO ivar_2412
   ENDIF
   analgpars(iset)%nodim = nodim
   IF (iopt_grid_nodim.EQ.2.AND.analgpars(iset)%nodim.EQ.3) THEN
      CALL warning_reset_arr_struc(analgpars(iset)%nodim,'analgpars','%nodim',&
                                 & 2,1,(/iset/))
   ENDIF
ENDDO iset_241

!---output grid parameters
CALL reset_out_gpars(analgpars,'analgpars',.FALSE.,.FALSE.)

!---station attributes
CALL reset_out_stats(analstatlocs)

!---file parameters
res0d%nodim = 0; res2d%nodim = 2; res3d%nodim = 3
amp0d%nodim = 0; amp2d%nodim = 2; amp3d%nodim = 3
pha0d%nodim = 0; pha2d%nodim = 2; pha3d%nodim = 3
ell2d%nodim = 2; ell3d%nodim = 3
iset_242: DO iset=1,nosetsanal
   IF (analgpars(iset)%nodim.LT.3) THEN
      CALL warning_reset_arr_struc(res3d(iset)%defined,'res3d','%defined',&
                                & .FALSE.,1,(/iset/))
      ifreq_2421: DO ifreq=1,nofreqsanal
         CALL warning_reset_arr_struc(amp3d(iset,ifreq)%defined,'amp3d',&
                                   & '%defined',.FALSE.,2,(/iset,ifreq/))
         CALL warning_reset_arr_struc(pha3d(iset,ifreq)%defined,'pha3d',&
                                   & '%defined',.FALSE.,2,(/iset,ifreq/))
         CALL warning_reset_arr_struc(ell3d(iset,ifreq)%defined,'ell3d',&
                                   & '%defined',.FALSE.,2,(/iset,ifreq/))
      ENDDO ifreq_2421
   ENDIF
   IF (analgpars(iset)%nodim.EQ.0) THEN
      CALL warning_reset_arr_struc(res2d(iset)%defined,'res2d','%defined',&
                                & .FALSE.,1,(/iset/))
      ifreq_2422: DO ifreq=1,nofreqsanal
         CALL warning_reset_arr_struc(amp2d(iset,ifreq)%defined,'amp2d',&
                                   & '%defined',.FALSE.,2,(/iset,ifreq/))
         CALL warning_reset_arr_struc(pha2d(iset,ifreq)%defined,'pha2d',&
                                   & '%defined',.FALSE.,2,(/iset,ifreq/))
         CALL warning_reset_arr_struc(ell2d(iset,ifreq)%defined,'ell2d',&
                                   & '%defined',.FALSE.,2,(/iset,ifreq/))
      ENDDO ifreq_2422
   ENDIF
ENDDO iset_242

!
!2.5 Check
!---------
!
!---time step for harmonic analysis
iset_251: DO iset=1,nosetsanal
   CALL check_mult_counter(icanal(iset),'icanal',ic3d)
ENDDO iset_251

!---output variables and vectors
CALL check_variables(analvars,'analvars')

!---file parameters
CALL check_out_filepars(res0d,'res0d')
CALL check_out_filepars(res2d,'res2d')
CALL check_out_filepars(res3d,'res3d')
CALL check_out_filepars(amp0d,'amp0d')
CALL check_out_filepars(amp2d,'amp2d')
CALL check_out_filepars(amp3d,'amp3d')
CALL check_out_filepars(pha0d,'pha0d')
CALL check_out_filepars(pha2d,'pha2d')
CALL check_out_filepars(pha3d,'pha3d')

!---output parameters
CALL check_out_gpars(analgpars,'analgpars',nostatsanal)

!---output stations
iset_252: DO iset=1,nosetsanal
   istat_2521: DO istat=1,analgpars(iset)%nostats
      CALL error_limits_arr(lstatsanal(iset,istat),'lstatsanal',1,&
                          & nostatsanal,2,(/iset,istat/))
   ENDDO istat_2521
ENDDO iset_252
CALL check_statlocs(analstatlocs,'analstatlocs')

!---analysed variable indices
iset_253: DO iset=1,nosetsanal
   ivar_2531: DO ivar=1,novarsanal
      CALL error_limits_arr(ivarsanal(iset,ivar),'ivarsanal',0,novarsanal,2,&
                         & (/iset,ivar/))
   ENDDO ivar_2531
   CALL error_diff_vals_arrlist(ivarsanal(iset,:),'ivarsanal',2,2,(/iset,1/),&
                             & .TRUE.)
ENDDO iset_253

!---elliptic variable indices
iset_254: DO iset=1,nosetsanal
   ivar_2541: DO ivar=1,14
      CALL error_limits_arr(ivarsell(iset,ivar),'ivarsell',0,14,2,&
                         & (/iset,ivar/))
   ENDDO ivar_2541
   CALL error_diff_vals_arrlist(ivarsell(iset,:),'ivarsell',2,2,(/iset,1/),&
                             & .TRUE.)
ENDDO iset_254

!
!3. Indices of analysed variables
!--------------------------------
!
!3.1 Number of variables per dimension
!-------------------------------------
!

novars0d = COUNT(analvars%nrank.EQ.0)
novars2d = COUNT(analvars%nrank.EQ.2)
novars3d = COUNT(analvars%nrank.EQ.3)

!
!3.2 Allocate
!------------
!

ALLOCATE (ivarsanal0d(nosetsanal,novars0d),STAT=errstat)
CALL error_alloc('ivarsanal0d',2,(/nosetsanal,novars0d/),kndint)
IF (novars0d.GT.0) ivarsanal0d = 0
ALLOCATE (ivarsanal2d(nosetsanal,novars2d),STAT=errstat)
CALL error_alloc('ivarsanal2d',2,(/nosetsanal,novars2d/),kndint)
IF (novars2d.GT.0) ivarsanal2d = 0
ALLOCATE (ivarsanal3d(nosetsanal,novars3d),STAT=errstat)
CALL error_alloc('ivarsanal3d',2,(/nosetsanal,novars3d/),kndint)
IF (novars3d.GT.0) ivarsanal3d = 0
maxvars = MAX(novars0d,novars2d,novars3d,7)
ALLOCATE (outvars(maxvars),STAT=errstat)
CALL error_alloc_struc('outvars',1,(/maxvars/),'VariableAtts')

!
!3.3 Define
!----------
!

iset_330: DO iset=1,nosetsanal

   i0 = 0; i2 = 0; i3 = 0
   ivar_331: DO ivar=1,novarsanal
      l = ivarsanal(iset,ivar)
      IF (l.GT.0) THEN
         nrank = analvars(l)%nrank
         SELECT CASE (nrank)
         CASE (0)
            i0 = i0 + 1
            ivarsanal0d(iset,i0) = l
         CASE (2)
            i2 = i2 + 1
            ivarsanal2d(iset,i2) = l
         CASE (3)
            i3 = i3 + 1
            ivarsanal3d(iset,i3) = l
         END SELECT
      ENDIF
   ENDDO ivar_331
   res0d(iset)%novars = i0; res2d(iset)%novars = i2
   res3d(iset)%novars = i3; amp0d(iset,:)%novars = i0
   amp2d(iset,:)%novars = i2; amp3d(iset,:)%novars = i3
   pha0d(iset,:)%novars = i0; pha2d(iset,:)%novars = i2
   pha3d(iset,:)%novars = i3
   ivarsanal2d(iset,:) = ivarsanal2d(iset,:) - novars0d
   ivarsanal3d(iset,:) = ivarsanal3d(iset,:) - novars0d - novars2d

ENDDO iset_330

!
!3.4 Write output only if number of variables is greater than zero
!-----------------------------------------------------------------
!

iset_340: DO iset=1,nosetsanal
   IF (res0d(iset)%novars.EQ.0) THEN
      CALL warning_reset_arr_struc(res0d(iset)%defined,'res0d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
   IF (res2d(iset)%novars.EQ.0) THEN
      CALL warning_reset_arr_struc(res2d(iset)%defined,'res2d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
   IF (res3d(iset)%novars.EQ.0) THEN
      CALL warning_reset_arr_struc(res3d(iset)%defined,'res3d','%defined',&
                                & .FALSE.,1,(/iset/))
   ENDIF
   ifreq_341: DO ifreq=1,nofreqsanal
      IF (amp0d(iset,ifreq)%novars.EQ.0) THEN
         CALL warning_reset_arr_struc(amp0d(iset,ifreq)%defined,'amp0d',&
                                   & '%defined',.FALSE.,2,(/iset,ifreq/))
      ENDIF
      IF (amp2d(iset,ifreq)%novars.EQ.0) THEN
         CALL warning_reset_arr_struc(amp2d(iset,ifreq)%defined,'amp2d',&
                                  & '%defined',.FALSE.,2,(/iset,ifreq/))
      ENDIF
      IF (amp3d(iset,ifreq)%novars.EQ.0) THEN
         CALL warning_reset_arr_struc(amp3d(iset,ifreq)%defined,'amp3d',&
                                   & '%defined',.FALSE.,2,(/iset,ifreq/))
      ENDIF
      IF (pha0d(iset,ifreq)%novars.EQ.0) THEN
         CALL warning_reset_arr_struc(pha0d(iset,ifreq)%defined,'pha0d',&
                                   & '%defined',.FALSE.,2,(/iset,ifreq/))
      ENDIF
      IF (pha2d(iset,ifreq)%novars.EQ.0) THEN
         CALL warning_reset_arr_struc(pha2d(iset,ifreq)%defined,'pha2d',&
                                   & '%defined',.FALSE.,2,(/iset,ifreq/))
      ENDIF
      IF (pha3d(iset,ifreq)%novars.EQ.0) THEN
         CALL warning_reset_arr_struc(pha3d(iset,ifreq)%defined,'pha3d',&
                                   & '%defined',.FALSE.,2,(/iset,ifreq/))
      ENDIF
   ENDDO ifreq_341

ENDDO iset_340

!
!4. Indices of elliptic variables
!--------------------------------
!
!4.1 Allocate
!------------
!

ALLOCATE (ivarsell2d(nosetsanal,7),STAT=errstat)
CALL error_alloc('ivarsell2d',2,(/nosetsanal,7/),kndint)
ivarsell2d = 0
ALLOCATE (ivarsell3d(nosetsanal,7),STAT=errstat)
CALL error_alloc('ivarsell3d',2,(/nosetsanal,7/),kndint)
ivarsell3d = 0

!
!4.2 Define
!----------
!

iset_420: DO iset=1,nosetsanal
   i2 = 0; i3 = 0
   ivar_421: DO ivar=1,14
      l = ivarsell(iset,ivar)
      IF (l.GT.0) THEN
         IF (l.LE.7) THEN
            i2 = i2 + 1
            ivarsell2d(iset,i2) = l
         ELSEIF (l.LE.14) THEN
            i3 = i3 + 1
            ivarsell3d(iset,i3) = l-7
         ENDIF
      ENDIF
   ENDDO ivar_421
   ell2d(iset,:)%novars = MERGE(i2,0,ell2d(iset,:)%defined)
   ell3d(iset,:)%novars = MERGE(i3,0,ell3d(iset,:)%defined)
ENDDO iset_420

!
!4.3 Write output only if number of variables is greater than zero
!-----------------------------------------------------------------
!

iset_430: DO iset=1,nosetsanal
ifreq_440: DO ifreq=1,nofreqsanal
   IF (ell2d(iset,ifreq)%novars.EQ.0) THEN
      CALL warning_reset_arr_struc(ell2d(iset,ifreq)%defined,'ell2d',&
                                & '%defined',.FALSE.,2,(/iset,ifreq/))
   ENDIF
   IF (ell3d(iset,ifreq)%novars.EQ.0) THEN
      CALL warning_reset_arr_struc(ell3d(iset,ifreq)%defined,'ell3d',&
                                & '%defined',.FALSE.,2,(/iset,ifreq/))
   ENDIF
ENDDO ifreq_440
ENDDO iset_430

!
!4.4 Check vector components for ellipses
!----------------------------------------
!

iset_440: DO iset=1,nosetsanal

!  ---2-D case
   IF (ANY(ell2d(iset,:)%defined)) THEN
      CALL error_diff_vals_arrlist(ivecell2d(iset,:),'ivecell2d',2,2,&
                                 & (/iset,1/),.FALSE.)
      novars = res2d(iset)%novars
      l_441: DO l=1,2
         CALL error_vals_arr_int(ivecell2d(iset,l),'ivecell2d',&
                               & novars,2,ivarsanal2d(iset,1:novars),(/iset,l/))
      ENDDO l_441
   ENDIF

!    ---3-D case
   IF (ANY(ell3d(iset,:)%defined)) THEN
      CALL error_diff_vals_arrlist(ivecell3d(iset,:),'ivecell3d',2,2,&
                                & (/iset,1/),.FALSE.)
      novars = res3d(iset)%novars
      l_442: DO l=1,2
         CALL error_vals_arr_int(ivecell3d(iset,l),'ivecell3d',&
                               & novars,2,ivarsanal3d(iset,1:novars),(/iset,l/))
      ENDDO l_442
   ENDIF

ENDDO iset_440

!
!5. Space output specifiers
!--------------------------
!
!5.1 Allocate and initialise
!---------------------------
!
!---regular output
ALLOCATE (nclocout(nosetsanal),STAT=errstat)
CALL error_alloc('nclocout',1,(/nosetsanal/),kndint)
nclocout = 0
ALLOCATE (nrlocout(nosetsanal),STAT=errstat)
CALL error_alloc('nrlocout',1,(/nosetsanal/),kndint)
nrlocout = 0
ALLOCATE (limprocs(2,2,nprocs,nosetsanal),STAT=errstat)
CALL error_alloc('limprocs',4,(/2,2,nprocs,nosetsanal/),kndint)
limprocs = 0
ALLOCATE (limloc(3,2,nosetsanal),STAT=errstat)
CALL error_alloc('limloc',3,(/3,2,nosetsanal/),kndint)
limloc = 0

!---station output
ALLOCATE (nostatsprocs(nprocs,nosetsanal),STAT=errstat)
CALL error_alloc('nostatsprocs',2,(/nprocs,nosetsanal/),kndint)
nostatsprocs = 0
ALLOCATE (nostatsloc(nosetsanal),STAT=errstat)
CALL error_alloc('nostatsloc',1,(/nosetsanal/),kndint)
nostatsloc = 0

!
!5.2 Regular output
!------------------
!

iset_520: DO iset=1,nosetsanal

   IF (analgpars(iset)%nodim.GT.0.AND.analgpars(iset)%gridded) THEN

      analgpars(iset)%ncout = lim_dims(analgpars(iset)%xlims)
      analgpars(iset)%nrout = lim_dims(analgpars(iset)%ylims)
      analgpars(iset)%nzout = lim_dims(analgpars(iset)%zlims)

      iproc_521: DO iproc=1,nprocs

!        ---X-direction
         xlims = analgpars(iset)%xlims
         nodat = 0; limprocs(1,2,iproc,iset) = 0; ii = 0
         IF (idloc.EQ.idprocs(iproc)) nclocout(iset) = 0
         i_5211: DO i=xlims(1),xlims(2),xlims(3)
            ii = ii + 1
            IF (i.GE.nc1procs(iproc).AND.i.LE.nc2procs(iproc)) THEN
               limprocs(1,2,iproc,iset) = ii
               nodat = nodat + 1
               IF (idloc.EQ.idprocs(iproc)) THEN
                  nclocout(iset) = nclocout(iset) + 1
                  i2 = i
               ENDIF
            ENDIF
         ENDDO i_5211
         limprocs(1,1,iproc,iset) = limprocs(1,2,iproc,iset) - nodat + 1
         IF (idloc.EQ.idprocs(iproc)) THEN
            IF (nclocout(iset).GT.0) THEN
               i2 = i2 - nc1loc + 1
               i1 = i2 + (1-nclocout(iset))*xlims(3)
               i3 = xlims(3)
            ELSE
               i1 = 1; i2 = 0; i3 = 1
            ENDIF
            limloc(:,1,iset) = (/i1,i2,i3/)
         ENDIF
            
!        ---Y-direction
         ylims = analgpars(iset)%ylims
         nodat = 0; limprocs(2,2,iproc,iset) = 0; jj = 0
         IF (idloc.EQ.idprocs(iproc)) nrlocout(iset) = 0
         j_5212: DO j=ylims(1),ylims(2),ylims(3)
            jj = jj + 1
            IF (j.GE.nr1procs(iproc).AND.j.LE.nr2procs(iproc)) THEN
               limprocs(2,2,iproc,iset) = jj
               nodat = nodat + 1
               IF (idloc.EQ.idprocs(iproc)) THEN
                  nrlocout(iset) = nrlocout(iset) + 1
                  j2 = j
               ENDIF
            ENDIF
         ENDDO j_5212
         limprocs(2,1,iproc,iset) = limprocs(2,2,iproc,iset) - nodat + 1
         IF (idloc.EQ.idprocs(iproc)) THEN
            IF (nrlocout(iset).GT.0) THEN
               j2 = j2 - nr1loc + 1
               j1 = j2 + (1-nrlocout(iset))*ylims(3)
               j3 = ylims(3)
            ELSE
               j1 = 1; j2 = 0; j3 = 1
            ENDIF
            limloc(:,2,iset) = (/j1,j2,j3/)
         ENDIF

      ENDDO iproc_521

   ENDIF

ENDDO iset_520
            
!
!5.3 Station output
!------------------
!
!5.3.1 Number of local stations
!----------------------------
!

iset_531: DO iset=1,nosetsanal
   IF (analgpars(iset)%nodim.GT.0.AND.(.NOT.analgpars(iset)%gridded)) THEN
      istat_5311: DO istat=1,analgpars(iset)%nostats
         l = lstatsanal(iset,istat)
         i = analstatlocs(l)%ipos; j = analstatlocs(l)%jpos
         iproc_53111: DO iproc=1,nprocs
            IF (local_proc(i,j,iproc=iproc)) THEN
               nostatsprocs(iproc,iset) = nostatsprocs(iproc,iset) + 1
               IF (idloc.EQ.idprocs(iproc)) THEN
                  nostatsloc(iset) = nostatsloc(iset) + 1
               ENDIF
            ENDIF
         ENDDO iproc_53111
      ENDDO istat_5311
   ENDIF
ENDDO iset_531

!
!5.3.2 Allocate
!--------------
!

maxstats = MAXVAL(nostatsprocs)
ALLOCATE (lstatsprocs(maxstats,nprocs,nosetsanal),STAT=errstat)
CALL error_alloc('lstatsprocs',3,(/maxstats,nprocs,nosetsanal/),kndint)
IF (maxstats.GT.0) lstatsprocs = 0

maxstatsloc = MAXVAL(nostatsloc)
ALLOCATE (istatpos(maxstatsloc,nosetsanal),STAT=errstat)
CALL error_alloc('istatpos',2,(/maxstatsloc,nosetsanal/),kndint)
IF (maxstatsloc.GT.0) istatpos = 0
ALLOCATE (jstatpos(maxstatsloc,nosetsanal),STAT=errstat)
CALL error_alloc('jstatpos',2,(/maxstatsloc,nosetsanal/),kndint)
IF (maxstatsloc.GT.0) jstatpos = 0

!
!5.3.3 Define
!------------
!

iproc_533: DO iproc=1,nprocs
iset_533: DO iset=1,nosetsanal
   IF (analgpars(iset)%nodim.GT.0.AND.(.NOT.analgpars(iset)%gridded)) THEN
      lloc = 0
      istat_5331: DO istat=1,analgpars(iset)%nostats
         l = lstatsanal(iset,istat)
         i = analstatlocs(l)%ipos; j = analstatlocs(l)%jpos
         IF (local_proc(i,j,iproc=iproc)) THEN
            lloc = lloc + 1
            lstatsprocs(lloc,iproc,iset) = l
            IF (idloc.EQ.idprocs(iproc)) THEN
               istatpos(lloc,iset) = i - nc1loc + 1
               jstatpos(lloc,iset) = j - nr1loc + 1
            ENDIF
         ENDIF
      ENDDO istat_5331
   ENDIF
ENDDO iset_533
ENDDO iproc_533

!
!5.4 Number of wet points
!------------------------
!

iset_540: DO iset=1,nosetsanal
   IF (analgpars(iset)%packing) THEN
      imax = lim_dims(limloc(:,1,iset))
      jmax = lim_dims(limloc(:,2,iset))
      ALLOCATE (maskvals(imax,jmax),STAT=errstat)
      CALL error_alloc('maskvals',2,(/imax,jmax/),kndlog)
      maskvals = mask_array(imax,jmax,limloc(:,:,iset),'C')
      nowetlocout = COUNT(maskvals)
      CALL sum_vars(nowetlocout,nowetout,0,commall=.TRUE.)
      analgpars(iset)%nowetout = nowetout   
      DEALLOCATE (maskvals)
   ELSE
      analgpars(iset)%nowetout = 0
   ENDIF
ENDDO iset_540

!
!5.5 Number of data
!------------------
!

ALLOCATE (ndat0da(nosetsanal),STAT=errstat)
CALL error_alloc('ndat0da',1,(/nosetsanal/),kndint)
ndat0da = 0
ALLOCATE (ndat2da(nosetsanal),STAT=errstat)
CALL error_alloc('ndat2da',1,(/nosetsanal/),kndint)
ndat2da = 0
ALLOCATE (ndat3da(nosetsanal),STAT=errstat)
CALL error_alloc('ndat3da',1,(/nosetsanal/),kndint)
ndat3da = 0
ALLOCATE (ndat0db(nosetsanal),STAT=errstat)
CALL error_alloc('ndat0db',1,(/nosetsanal/),kndint)
ndat0db = 0
ALLOCATE (ndat2db(nosetsanal),STAT=errstat)
CALL error_alloc('ndat2db',1,(/nosetsanal/),kndint)
ndat2db = 0
ALLOCATE (ndat3db(nosetsanal),STAT=errstat)
CALL error_alloc('ndat3db',1,(/nosetsanal/),kndint)
ndat3db = 0

!
!6. Time specifiers
!------------------
!
!6.1 Allocate and initialise
!---------------------------
!

ALLOCATE (icount(nosetsanal),STAT=errstat)
CALL error_alloc('icount',1,(/nosetsanal/),kndint)
icount = 0
ALLOCATE (nanaltot(nosetsanal),STAT=errstat)
CALL error_alloc('nanaltot',1,(/nosetsanal/),kndint)
ALLOCATE (ntout(nosetsanal),STAT=errstat)
CALL error_alloc('ntout',1,(/nosetsanal/),kndint)
ntout = 0

!
!6.2 Initialise
!--------------
!

iset_620: DO iset=1,nosetsanal
   nanaltot(iset) = analgpars(iset)%tlims(3)/icanal(iset) + 1
ENDDO iset_620

!
!7. Work space arrays for analysis
!---------------------------------
!
!7.1 Allocate and initialise
!---------------------------
!

ALLOCATE (acoef(0:nofreqsanal,0:nofreqsanal,nosetsanal),STAT=errstat)
CALL error_alloc('acoef',3,(/nofreqsanal+1,nofreqsanal+1,nosetsanal/),&
               & kndrtype)
acoef = 0.0
ALLOCATE (bcoef(nofreqsanal,nofreqsanal,nosetsanal),STAT=errstat)
CALL error_alloc('bcoef',3,(/nofreqsanal,nofreqsanal,nosetsanal/),kndrtype)
bcoef = 0.0

!
!7.2 Define
!----------
!

iset_720: DO iset=1,nosetsanal
   deltanal = icanal(iset)*delt2d
   nofreqs = nofreqsharm(iset)

!  ---matrices for linear equations
   acoef(0,0,iset) = nanaltot(iset)
   m_721: DO m=1,nofreqs
      mm = ifreqsharm(iset,m)
      IF (MOD(harm_freq(mm)*deltanal,twopi).NE.0.0) THEN
         acoef(0,m,iset) = SIN(0.5*nanaltot(iset)*harm_freq(mm)*deltanal)/&
                         & SIN(0.5*harm_freq(mm)*deltanal)
      ELSE
         acoef(0,m,iset) = nanaltot(iset)
      ENDIF
      acoef(m,0,iset) = acoef(0,m,iset)
      n_7211: DO n=1,m
         nn = ifreqsharm(iset,n)
         sig = harm_freq(mm) + harm_freq(nn)
         IF (MOD(sig*deltanal,twopi).NE.0.0) THEN
            cfplus = SIN(0.5*nanaltot(iset)*sig*deltanal)/SIN(0.5*sig*deltanal)
         ELSE
            cfplus = nanaltot(iset)
         ENDIF
         sig = harm_freq(mm) - harm_freq(nn)
         IF (MOD(sig*deltanal,twopi).NE.0.0) THEN
            cfmin = SIN(0.5*nanaltot(iset)*sig*deltanal)/SIN(0.5*sig*deltanal)
         ELSE
            cfmin = nanaltot(iset)
         ENDIF
         acoef(m,n,iset) = 0.5*(cfplus+cfmin)
         bcoef(m,n,iset) = 0.5*(cfmin-cfplus)
         IF (n.LT.m) THEN
            acoef(n,m,iset) = acoef(m,n,iset)
            bcoef(n,m,iset) = bcoef(m,n,iset)
         ENDIF
      ENDDO n_7211
   ENDDO m_721

!  ---LU-decomposition
   CALL cholesky_decomp(acoef(:,:,iset),nofreqs+1,nofreqsanal+1,'acoef',ierr,&
                     & .TRUE.)
   CALL cholesky_decomp(bcoef(:,:,iset),nofreqs,nofreqsanal,'bcoef',ierr,&
                     & .TRUE.)
   IF (ierr.NE.0) CALL error_abort('harmonic_analysis_init',ierrno_runval)

ENDDO iset_720

!
!7.3 Number of output data for buffering
!---------------------------------------
!

novals0da_old = 0; novals2da_old = 0; novals3da_old = 0
novals0db_old = 0; novals2db_old = 0; novals3db_old = 0

!
!8. Abort if necessary
!----------------------
!

CALL error_abort('harmonic_analysis_init',ierrno_inival)

!
!9. Write file headers and initial arrays
!----------------------------------------
!

CALL inout_atts_out('HR')
CALL inout_atts_out('HA')
CALL inout_atts_out('HP')
CALL inout_atts_out('HE')

CALL log_timer_out()


RETURN

END SUBROUTINE harmonic_analysis_init

!========================================================================

SUBROUTINE harmonic_analysis_grid
!************************************************************************
!
! *harmonic_analysis_grid* Write output grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Harmonic_Analysis.f90  V2.8
!
! Description -
!
! Reference -
!
! Calling program - harmonic_analysis
!
! External calls -
!
! Module calls - add_secs_to_date, combine_mod, date_to_year, diff_dates,
  !                error_alloc, global_mask, loop_index, write_time, Zcoord_var
!
!************************************************************************
!
!1. Log info
!-----------
!

flag = .FALSE.
iset_110: DO iset=1,nosetsanal
   IF (loop_index(analgpars(iset)%tlims,nt)) THEN
      flag = .TRUE.
   ENDIF
ENDDO iset_110
IF (flag) THEN
   procname(pglev+1) = 'harmonic_analysis_grid'
   CALL log_timer_in()
ELSE
   RETURN
ENDIF

!
!2. First output time
!--------------------
!

iset_200: DO iset=1,nosetsanal

   nodim = analgpars(iset)%nodim
   IF (nodim.GT.0.AND.nt.EQ.analgpars(iset)%tlims(1)) THEN

      gridded = analgpars(iset)%gridded
      packing = analgpars(iset)%packing
      ncout = analgpars(iset)%ncout
      nrout = analgpars(iset)%nrout
      nzout = analgpars(iset)%nzout
      nostats = analgpars(iset)%nostats
      nowetout = analgpars(iset)%nowetout
      vcoord = analgpars(iset)%vcoord
      fill_value = MERGE(double_fill,0.0_kndrlong,packing)
      xlimsglb = analgpars(iset)%xlims
      ylimsglb = analgpars(iset)%ylims
      xlims = limloc(:,1,iset)
      ylims = limloc(:,2,iset)
      zlims = analgpars(iset)%zlims
      lprocs = limprocs(:,:,:,iset)
      fill_value = double_fill
      
!
!2.1 Allocate arrays
!-------------------
!
!2.1.1 Land mask arrays
!----------------------
!      

      IF (packing) THEN
         ALLOCATE(maskglb(nc,nr),STAT=errstat)
         CALL error_alloc('maskglb',2,(/nc,nr/),kndlog)
         ALLOCATE(lmask(nowetout),STAT=errstat)
         CALL error_alloc('lmask',1,(/nowetout/),kndint)
      ENDIF
 
!
!2.1.2 Gridded
!-------------
!

      IF (gridded) THEN
         
!        ---horizontal grid
         ALLOCATE (xout(ncout,nrout),STAT=errstat)
         CALL error_alloc('xout',2,(/ncout,nrout/),kndrtype)
         ALLOCATE (yout(ncout,nrout),STAT=errstat)
         CALL error_alloc('yout',2,(/ncout,nrout/),kndrtype)
         
!        ---bathymetry
         ALLOCATE (depout(nclocout(iset),nrlocout(iset)),STAT=errstat)
         CALL error_alloc('depout',2,(/nclocout(iset),nrlocout(iset)/),kndrtype)
         IF (SIZE(depout).GT.0) depout = fill_value

!        ---vertical grid
         IF (nodim.EQ.3) THEN
            IF (vcoord.EQ.1) THEN
               ALLOCATE (zout(nclocout(iset),nrlocout(iset),nzout),STAT=errstat)
               CALL error_alloc('zout',3,&
                             & (/nclocout(iset),nrlocout(iset),nzout/),kndrtype)
               IF (SIZE(zout).GT.0) zout = fill_value
            ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.LE.2) THEN
               ALLOCATE (gsigout(nzout),STAT=errstat)
               CALL error_alloc('gsigout',1,(/nzout/),kndrtype)
            ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
               ALLOCATE (gsout(nclocout(iset),nrlocout(iset),nzout),&
                       & STAT=errstat)
               CALL error_alloc('gsout',3,&
                             & (/nclocout(iset),nrlocout(iset),nzout/),kndrtype)
               IF (SIZE(gsout).GT.0) gsout = fill_value
            ENDIF
         ENDIF
 
!
!2.1.3 Non-gridded
!-----------------      
!
         
      ELSE

!        ---horizontal grid         
         ALLOCATE (xout(nostats,1),STAT=errstat)
         CALL error_alloc('xout',2,(/nostats,1/),kndrtype)
         ALLOCATE (yout(nostats,1),STAT=errstat)
         CALL error_alloc('yout',2,(/nostats,1/),kndrtype)

!        ---bathymetry         
         ALLOCATE (depout(nostatsloc(iset),1),STAT=errstat)
         CALL error_alloc('depout',2,(/nostatsloc(iset),1/),kndrtype)

!        ---vertical grid
         IF (nodim.EQ.3) THEN
            IF (vcoord.EQ.1) THEN
               ALLOCATE (zout(nostatsloc(iset),1,nzout),STAT=errstat)
               CALL error_alloc('zout',3,(/nostatsloc(iset),1,nzout/),kndrtype)
            ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.LE.2) THEN
               ALLOCATE (gsigout(nzout),STAT=errstat)
               CALL error_alloc('gsigout',1,(/nzout/),kndrtype)
            ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
               ALLOCATE (gsout(nostatsloc(iset),1,nzout),STAT=errstat)
               CALL error_alloc('gsout',3,(/nostatsloc(iset),1,nzout/),kndrtype)
            ENDIF
         ENDIF

      ENDIF

!
!2.2 Output grid
!---------------
!
!2.2.1 Index vecor of sea points
!-------------------------------         
!

      IF (packing) THEN
         CALL global_mask(maskglb,'C  ',1)            
         l = 0
         j_221: DO j=ylimsglb(1),ylimsglb(2),ylimsglb(3)
         i_221: DO i=xlimsglb(1),xlimsglb(2),xlimsglb(3)
            IF (maskglb(i,j)) THEN 
               l = l + 1
               lmask(l) = (j-1)*ncout + i
            ENDIF
         ENDDO i_221
         ENDDO j_221
      ENDIF
         
!         
!2.2.2 Gridded
!-------------
!      

      IF (gridded) THEN

!        ---horizontal grid         
         j_2221: DO j=ylimsglb(1),ylimsglb(2),ylimsglb(3)
         i_2221: DO i=xlimsglb(1),xlimsglb(2),xlimsglb(3)
            ii = (i-xlimsglb(1))/xlimsglb(3) + 1
            jj = (j-ylimsglb(1))/ylimsglb(3) + 1
            xout(ii,jj) = gxcoordglbatc(i,j)
            yout(ii,jj) = gycoordglbatc(i,j)
         ENDDO i_2221
         ENDDO j_2221
            
!        ---bathymetry and vertical grid         
         j_2222: DO j=ylims(1),ylims(2),ylims(3)
         i_2222: DO i=xlims(1),xlims(2),xlims(3)
            ii = (i-xlims(1))/xlims(3) + 1
            jj = (j-ylims(1))/ylims(3) + 1      
            depout(ii,jj) = MERGE(depmeanatc(i,j),fill_value,maskatc_int(i,j))
            IF (nodim.EQ.3) THEN
               IF (maskatc_int(i,j)) THEN
                  k_22221: DO k=zlims(1),zlims(2),zlims(3)
                     kk = (k-zlims(1))/zlims(3) + 1
                     IF (vcoord.EQ.1) THEN
                        zout(ii,jj,kk) = Zcoord_var(i,j,k,'C  ',.TRUE.)
                     ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
                        gsout(ii,jj,kk) = gscoordatc(i,j,k) - 1.0
                     ENDIF
                  ENDDO k_22221
               ENDIF
            ENDIF
         ENDDO i_2222
         ENDDO j_2222

!
!2.2.3 Non-gridded
!-----------------
!         

      ELSE

!        ---horizontal grid
         l_2231: DO l=1,nostats
            i = analstatlocs(l)%ipos; j = analstatlocs(l)%jpos
            xout(l,1) = gxcoordglbatc(i,j); yout(l,1) = gycoordglbatc(i,j)
         ENDDO l_2231
         
!        ---bathymetry and vertical grid 
         l_2232: DO l=1,nostatsloc(iset)
            i = istatpos(l,iset); j = jstatpos(l,iset)
            depout(l,1) = depmeanatc(i,j)
            IF (nodim.EQ.3) THEN
               k_22321: DO k=zlims(1),zlims(2),zlims(3)
                  kk = (k-zlims(1))/zlims(3) + 1
                  IF (vcoord.EQ.1) THEN
                     zout(l,1,kk) = Zcoord_var(i,j,k,'C  ',.TRUE.)
                  ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
                     gsout(l,1,kk) = gscoordatc(i,j,k) - 1.0
                  ENDIF
               ENDDO k_22321
            ENDIF
         ENDDO l_2232

      ENDIF

!
!2.2.4 Sigma coordinate
!----------------------
!      

      IF (nodim.EQ.3.AND.vcoord.EQ.2.AND.iopt_grid_vtype.LE.2) THEN
         k_224: DO k=zlims(1),zlims(2),zlims(3)
            kk = (k-zlims(1))/zlims(3) + 1
            gsigout(kk)= gsigcoordatc(k) - 1.0
         ENDDO k_224
      ENDIF
      
!
!2.3 Write
!---------
!

      l_230: DO l=2,3
      itype_231: DO itype=1,4
         nofreqs = MERGE(1,nofreqsharm(iset),itype.EQ.1)
         ifreq_2311: DO ifreq=1,nofreqs
            SELECT CASE (l)
               CASE (2)
                  SELECT CASE (itype)
                     CASE (1); filepars = res2d(iset)
                     CASE (2); filepars = amp2d(iset,ifreq)
                     CASE (3); filepars = pha2d(iset,ifreq)
                     CASE (4); filepars = ell2d(iset,ifreq)
                  END SELECT
               CASE (3)
                  SELECT CASE (itype)
                     CASE (1); filepars = res3d(iset)
                     CASE (2); filepars = amp3d(iset,ifreq)
                     CASE (3); filepars = pha3d(iset,ifreq)
                     CASE (4); filepars = ell3d(iset,ifreq)
                  END SELECT
            END SELECT
            IF (.NOT.filepars%defined) CYCLE ifreq_2311
            ivar = 1
            
!
!2.3.1 Gridded
!--------------            
!

            IF (gridded) THEN
!              ---horizontal grid
               IF (master) THEN
                  IF (iopt_grid_htype.LE.2) THEN
                     CALL write_vars(xout(:,1),filepars,ivar+1)
                     CALL write_vars(yout(1,:),filepars,ivar+2)
                  ELSE
                     CALL write_vars(xout,filepars,ivar+1)
                     CALL write_vars(yout,filepars,ivar+2)
                  ENDIF
               ENDIF
               ivar = ivar + 2
!              ---index vector of sea points               
               IF (packing) THEN
                  IF (master) CALL write_vars(lmask,filepars,ivar+1)
                  ivar = ivar + 1
               ENDIF
!              ---bathymetry
               CALL combine_write_submod(depout,filepars,ivar+1,&
                                   & (/ncout,nrout/),lprocs,nowet=nowetout)
               ivar = ivar + 1
!              ---vertical grid
               IF (filepars%nodim.EQ.3) THEN
                  IF (vcoord.EQ.1) THEN
                     CALL combine_write_submod(zout,filepars,ivar+1,&
                                  & (/ncout,nrout,nzout/),lprocs,nowet=nowetout)
                     ivar = ivar + 1
                  ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
                     CALL combine_write_submod(gsout,filepars,ivar+1,&
                                  & (/ncout,nrout,nzout/),lprocs,nowet=nowetout)
                     ivar = ivar + 1
                  ENDIF
               ENDIF
            
!
!2.3.2 Non-gridded
!-----------------            
!
               
            ELSE
!              ---horizontal grid
               IF (master) THEN
                  CALL write_vars(xout(:,1),filepars,ivar+1)
                  CALL write_vars(yout(:,1),filepars,ivar+2)
                  ivar = ivar + 2
               ENDIF
!              ---bathymetry
               CALL combine_write_stats_loc(depout(:,1),filepars,&
                                   & ivar+1,maxstats,nostats,&
                                   & nostatsprocs(:,iset),lstatsprocs(:,:,iset))
               ivar = ivar + 1
!              ---vertical grid
               IF (filepars%nodim.EQ.3) THEN               
                  IF (vcoord.EQ.1) THEN
                     CALL combine_write_stats_loc(zout(:,1,:),filepars,&
                                & ivar+1,maxstats,nostats,nostatsprocs(:,iset),&
                                & lstatsprocs(:,:,iset),reduced=.FALSE.)
                     ivar = ivar + 1
                  ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
                     CALL combine_write_stats_loc(gsout(:,1,:),filepars,&
                                & ivar+1,maxstats,nostats,nostatsprocs(:,iset),&
                                & lstatsprocs(:,:,iset),reduced=.FALSE.)
                     ivar = ivar + 1
                  ENDIF
               ENDIF
                  
            ENDIF

!
!2.3.3 Sigma coordinate
!----------------------
!
            
            IF (filepars%nodim.EQ.3.AND.vcoord.EQ.2.AND.&
                 & iopt_grid_vtype.EQ.1) THEN
               IF (master) CALL write_vars(gsigout,filepars,ivar+1)
               ivar = ivar + 1
            ENDIF

         ENDDO ifreq_2311
      ENDDO itype_231
         
      ENDDO l_230

      
!
!2.4 Deallocate
!--------------
!

      DEALLOCATE (xout,yout,depout)
      IF (packing) DEALLOCATE (maskglb,lmask)
      IF (nodim.EQ.3) THEN
         IF (vcoord.EQ.1) THEN
            DEALLOCATE (zout)
         ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.LE.2) THEN
            DEALLOCATE (gsigout)
         ELSEIF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.3) THEN
            DEALLOCATE (gsout)
         ENDIF
      ENDIF
   
   ENDIF
      
ENDDO iset_200

!
!3. Time data
!------------
!

iset_300: DO iset=1,nosetsanal

   IF (loop_index(analgpars(iset)%tlims,nt).AND.&
     & nt.GT.analgpars(iset)%tlims(1)) THEN
      time_format = analgpars(iset)%time_format
      refdate = analgpars(iset)%refdate
      CALL add_secs_to_date(CDateTime,ctime,-analgpars(iset)%tlims(3),&
                          & 0.5*delt2d)

!
!3.1 Time coordinate
!-------------------
!

      IF (time_format.EQ.8) THEN
         CALL date_to_year(ctime,rtime)
      ELSEIF (time_format.GT.0) THEN
         CALL diff_dates(refdate,ctime,time_format,rtime=rtime)
      ENDIF

!
!3.2 Write
!---------
!

      l_320: DO l=0,3
         IF (l.EQ.1) CYCLE l_320
         itype_321: DO itype=1,4         
            nofreqs = MERGE(1,nofreqsharm(iset),itype.EQ.1)
            ifreq_3211: DO ifreq=1,nofreqs
               SELECT CASE (l)
                  CASE (0)
                     SELECT CASE (itype)
                        CASE (1); filepars = res0d(iset)
                        CASE (2); filepars = amp0d(iset,ifreq)
                        CASE (3); filepars = pha0d(iset,ifreq)
                     END SELECT
                     flag = MERGE(filepars%defined,.FALSE.,itype.LT.4)
                  CASE (2)
                     SELECT CASE (itype)
                        CASE (1); filepars = res2d(iset)
                        CASE (2); filepars = amp2d(iset,ifreq)
                        CASE (3); filepars = pha2d(iset,ifreq)
                        CASE (4); filepars = ell2d(iset,ifreq)
                     END SELECT
                     flag = filepars%defined
                  CASE (3)
                     SELECT CASE (itype)
                        CASE (1); filepars = res3d(iset)
                        CASE (2); filepars = amp3d(iset,ifreq)
                        CASE (3); filepars = pha3d(iset,ifreq)
                        CASE (4); filepars = ell3d(iset,ifreq)
                     END SELECT
                     flag = filepars%defined
               END SELECT
               IF (.NOT.flag.OR.(l.EQ.0.AND.itype.EQ.4)) CYCLE ifreq_3211

               IF (time_format.EQ.0) THEN
                  CALL write_time(ctime,filepars)
               ELSE
                  CALL write_time(rtime,filepars)
               ENDIF
!              ---store changes in file structure
               SELECT CASE (l)
                  CASE (0)
                     SELECT CASE (itype)
                     CASE (1); res0d(iset) = filepars
                     CASE (2); amp0d(iset,ifreq) = filepars
                     CASE (3); pha0d(iset,ifreq) = filepars
                     END SELECT
                  CASE (2)
                     SELECT CASE (itype)
                     CASE (1); res2d(iset) = filepars
                     CASE (2); amp2d(iset,ifreq) = filepars
                     CASE (3); pha2d(iset,ifreq) = filepars
                     CASE (4); ell2d(iset,ifreq) = filepars
                     END SELECT
                  CASE (3)
                     SELECT CASE (itype)
                     CASE (1); res3d(iset) = filepars
                     CASE (2); amp3d(iset,ifreq) = filepars
                     CASE (3); pha3d(iset,ifreq) = filepars
                     CASE (4); ell3d(iset,ifreq) = filepars
                     END SELECT
               END SELECT
            ENDDO ifreq_3211
         ENDDO itype_321
      ENDDO l_320
   ENDIF
   
ENDDO iset_300

CALL log_timer_out()


RETURN

END SUBROUTINE harmonic_analysis_grid

!========================================================================

SUBROUTINE harmonic_analysis_reset
!************************************************************************
!
! *harmonic_analysis_reset* Initialise buffers with harmonic data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Harmonic_Analysis.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - harmonic_analysis
!
! Module calls - loop_index
!
!************************************************************************
!


procname(pglev+1) = 'harmonic_analysis_reset'
CALL log_timer_in()

lbuf0da = 0; lbuf2da = 0; lbuf3da = 0
lbuf0db = 0; lbuf2db = 0; lbuf3db = 0

iset_1000: DO iset=1,nosetsanal
   tlims = analgpars(iset)%tlims
   l0data = ndat0da(iset); l2data = ndat2da(iset); l3data = ndat3da(iset)
   l0datb = ndat0db(iset); l2datb = ndat2db(iset); l3datb = ndat3db(iset)
   IF ((first.AND.(nt.EQ.tlims(1))).OR.&
     & (second.AND.loop_index(tlims,nt).AND.&
     & nt.GT.tlims(1).AND.nt.LT.tlims(2))) THEN
      IF (l0data.GT.0) buff0da(lbuf0da+1:lbuf0da+l0data) = 0.0
      IF (l2data.GT.0) buff2da(lbuf2da+1:lbuf2da+l2data) = 0.0
      IF (l3data.GT.0) buff3da(lbuf3da+1:lbuf3da+l3data) = 0.0
      IF (l0datb.GT.0) buff0db(lbuf0db+1:lbuf0db+l0datb) = 0.0
      IF (l2datb.GT.0) buff2db(lbuf2db+1:lbuf2db+l2datb) = 0.0
      IF (l3datb.GT.0) buff3db(lbuf3db+1:lbuf3db+l3datb) = 0.0
      icount(iset) = -tlims(3)/(2*icanal(iset))
   ENDIF
   lbuf0da = lbuf0da + l0data
   lbuf2da = lbuf2da + l2data
   lbuf3da = lbuf3da + l3data
   lbuf0db = lbuf0db + l0datb
   lbuf2db = lbuf2db + l2datb
   lbuf3db = lbuf3db + l3datb
ENDDO iset_1000

CALL log_timer_out()


RETURN

END SUBROUTINE harmonic_analysis_reset

!========================================================================

SUBROUTINE harmonic_analysis_update
!************************************************************************
!
! *harmonic_analysis_update* Update buffers with harmonic data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Harmonic_Analysis.f90  V2.9
!
! Description -
!
! Reference -
!
! Calling program - harmonic_analysis
!
! External calls - usrdef_anal0d_vals, usrdef_anal2d_vals, usrdef_anal3d_vals
!
! Module calls - define_out0d_vals, define_out2d_vals, define_out3d_vals,
!                loop_index
!
!************************************************************************
!


procname(pglev+1) = 'harmonic_analysis_update'
CALL log_timer_in()

lbuf0da = 0; lbuf2da = 0; lbuf3da = 0
lbuf0db = 0; lbuf2db = 0; lbuf3db = 0

fcos(0) = 1.0

iset_1000: DO iset=1,nosetsanal
   deltanal = icanal(iset)*delt2d
   tlims = analgpars(iset)%tlims
   IF ((first.AND.nt.GE.tlims(1).AND.nt.LE.tlims(2).AND.&
     & mult_index(nt-tlims(1),icanal(iset))).OR.&
     & (second.AND.loop_index(tlims,nt).AND.&
     & nt.GT.tlims(1).AND.nt.LT.tlims(2))) THEN

!
!1. Initialise
!-------------
!
!     ---resolution
      xlims = limloc(:,1,iset)
      ylims = limloc(:,2,iset)
      zlims = analgpars(iset)%zlims

!     ---number of variables
      nanalvars0d = res0d(iset)%novars
      nanalvars2d = res2d(iset)%novars
      nanalvars3d = res3d(iset)%novars

!     --- frequency arrays
      nofreqs = nofreqsharm(iset)
      n_110: DO n=1,nofreqs
         ifreq = ifreqsharm(iset,n)
         xarg = MOD(harm_freq(ifreq)*icount(iset)*deltanal,twopi)
         fcos(n) = COS(xarg)
         fsin(n) = SIN(xarg)
      ENDDO n_110

!
!2. 0-D case
!-----------
!

      IF (nanalvars0d.GT.0) THEN
         outvars(1:nanalvars0d) = analvars(ivarsanal0d(iset,1:nanalvars0d)) 
         CALL define_out0d_vals(outsub(1:nanalvars0d),nanalvars0d,&
                              & outvars=outvars(1:nanalvars0d))
         IF (ANY(outvars(1:nanalvars0d)%ivarid.EQ.0)) THEN
            CALL usrdef_anal0d_vals(outdat(1:novars0d),novars0d)
            WHERE (outvars(1:nanalvars0d)%ivarid.EQ.0)
               outsub(1:nanalvars0d) = outdat(ivarsanal0d(iset,1:nanalvars0d))
            END WHERE
         ENDIF
         n_210: DO n=0,nofreqs
            buff0da(lbuf0da+1:lbuf0da+nanalvars0d) = &
               & buff0da(lbuf0da+1:lbuf0da+nanalvars0d) + &
               & outsub(1:nanalvars0d)*fcos(n)
            lbuf0da = lbuf0da + nanalvars0d
         ENDDO n_210
         n_220: DO n=1,nofreqs
            buff0db(lbuf0db+1:lbuf0db+nanalvars0d) = &
               & buff0db(lbuf0db+1:lbuf0db+nanalvars0d) + &
               & outsub(1:nanalvars0d)*fsin(n)
            lbuf0db = lbuf0db + nanalvars0d
         ENDDO n_220

      ENDIF

!
!3. 2-D case
!-----------
!

      IF (nanalvars2d.GT.0) THEN
         outvars(1:nanalvars2d) = &
                           & analvars(novars0d+ivarsanal2d(iset,1:nanalvars2d))
         flag = ANY(outvars(1:nanalvars2d)%ivarid.EQ.0)

!
!3.1 Regular data
!----------------
!

         IF (analgpars(iset)%gridded) THEN
            i_310: DO i=xlims(1),xlims(2),xlims(3)
            j_310: DO j=ylims(1),ylims(2),ylims(3)
               CALL define_out2d_vals(outsub(1:nanalvars2d),i,j,nanalvars2d,&
                                    & outvars=outvars(1:nanalvars2d))
               IF (flag) THEN
                  CALL usrdef_anal2d_vals(outdat(1:novars2d),i,j,novars2d)
                  WHERE (outvars(1:nanalvars2d)%ivarid.EQ.0)
                     outsub(1:nanalvars2d) = &
                                       & outdat(ivarsanal2d(iset,1:nanalvars2d))
                  END WHERE
               ENDIF
               n_311: DO n=0,nofreqs
                  buff2da(lbuf2da+1:lbuf2da+nanalvars2d) = &
                         & buff2da(lbuf2da+1:lbuf2da+nanalvars2d) + &
                         & outsub(1:nanalvars2d)*fcos(n)
                  lbuf2da = lbuf2da + nanalvars2d
               ENDDO n_311
               n_312: DO n=1,nofreqs
                  buff2db(lbuf2db+1:lbuf2db+nanalvars2d) = &
                         & buff2db(lbuf2db+1:lbuf2db+nanalvars2d) + &
                         & outsub(1:nanalvars2d)*fsin(n)
                  lbuf2db = lbuf2db + nanalvars2d
               ENDDO n_312
            ENDDO j_310
            ENDDO i_310
      
!
!3.2 Irregular data
!------------------
!

         ELSE
            l_320: DO l=1,nostatsloc(iset)
               i = istatpos(l,iset); j = jstatpos(l,iset)
               CALL define_out2d_vals(outsub(1:nanalvars2d),i,j,nanalvars2d,&
                                    & outvars=outvars(1:nanalvars2d))
               IF (flag) THEN
                  CALL usrdef_anal2d_vals(outdat(1:novars2d),i,j,novars2d)
                  WHERE (outvars(1:nanalvars2d)%ivarid.EQ.0)
                     outsub(1:nanalvars2d) = &
                                       & outdat(ivarsanal2d(iset,1:nanalvars2d))
                  END WHERE
               ENDIF
               n_321: DO n=0,nofreqs
                  buff2da(lbuf2da+1:lbuf2da+nanalvars2d) = &
                         & buff2da(lbuf2da+1:lbuf2da+nanalvars2d) + &
                         & outsub(1:nanalvars2d)*fcos(n)
                  lbuf2da = lbuf2da + nanalvars2d
               ENDDO n_321
               n_322: DO n=1,nofreqs
                  buff2db(lbuf2db+1:lbuf2db+nanalvars2d) = &
                         & buff2db(lbuf2db+1:lbuf2db+nanalvars2d) + &
                         & outsub(1:nanalvars2d)*fsin(n)
                  lbuf2db = lbuf2db + nanalvars2d
               ENDDO n_322
            ENDDO l_320
         ENDIF

      ENDIF

!
!4. 3-D case
!-----------
!

      IF (nanalvars3d.GT.0) THEN
         outvars(1:nanalvars3d) = &
                   & analvars(novars0d+novars2d+ivarsanal3d(iset,1:nanalvars3d))

         flag = ANY(outvars(1:nanalvars3d)%ivarid.EQ.0)

!
!4.1 Regular data
!----------------
!

         IF (analgpars(iset)%gridded) THEN
            i_410: DO i=xlims(1),xlims(2),xlims(3)
            j_410: DO j=ylims(1),ylims(2),ylims(3)
            k_410: DO k=zlims(1),zlims(2),zlims(3)
               CALL define_out3d_vals(outsub(1:nanalvars3d),i,j,k,nanalvars3d,&
                                    & outvars=outvars(1:nanalvars3d))
               IF (flag) THEN
                  CALL usrdef_anal3d_vals(outdat(1:novars3d),i,j,k,novars3d)
                  WHERE (outvars(1:nanalvars3d)%ivarid.EQ.0)
                     outsub(1:nanalvars3d) = &
                                       & outdat(ivarsanal3d(iset,1:nanalvars3d))
                  END WHERE
               ENDIF
               n_4111: DO n=0,nofreqs
                  buff3da(lbuf3da+1:lbuf3da+nanalvars3d) = &
                          & buff3da(lbuf3da+1:lbuf3da+nanalvars3d) + &
                          & outsub(1:nanalvars3d)*fcos(n)
                  lbuf3da = lbuf3da + nanalvars3d
               ENDDO n_4111
               n_4112: DO n=1,nofreqs
                  buff3db(lbuf3db+1:lbuf3db+nanalvars3d) = &
                          & buff3db(lbuf3db+1:lbuf3db+nanalvars3d) + &
                          & outsub(1:nanalvars3d)*fsin(n)
                  lbuf3db = lbuf3db + nanalvars3d
               ENDDO n_4112

            ENDDO k_410
            ENDDO j_410
            ENDDO i_410
      
!
!4.2 Irregular data
!------------------
!

         ELSE
            l_420: DO l=1,nostatsloc(iset)
               i = istatpos(l,iset); j = jstatpos(l,iset)
               k_421: DO k=zlims(1),zlims(2),zlims(3)
                  CALL define_out3d_vals(outsub(1:nanalvars3d),i,j,k,&
                                   & nanalvars3d,outvars=outvars(1:nanalvars3d))
                  IF (flag) THEN
                     CALL usrdef_anal3d_vals(outdat(1:novars3d),i,j,k,novars3d)
                     WHERE (outvars(1:nanalvars3d)%ivarid.EQ.0)
                        outsub(1:nanalvars3d) = &
                                       & outdat(ivarsanal3d(iset,1:nanalvars3d))
                     END WHERE
                  ENDIF
                  n_4211: DO n=0,nofreqs
                     buff3da(lbuf3da+1:lbuf3da+nanalvars3d) = &
                            & buff3da(lbuf3da+1:lbuf3da+nanalvars3d) + &
                            & outsub(1:nanalvars3d)*fcos(n)
                     lbuf3da = lbuf3da + nanalvars3d
                  ENDDO n_4211
                  n_4212: DO n=1,nofreqs
                     buff3db(lbuf3db+1:lbuf3db+nanalvars3d) = &
                            & buff3db(lbuf3db+1:lbuf3db+nanalvars3d) + &
                            & outsub(1:nanalvars3d)*fsin(n)
                     lbuf3db = lbuf3db + nanalvars3d
                  ENDDO n_4212
               ENDDO k_421
            ENDDO l_420
         ENDIF

      ENDIF

!     ---update counter
      icount(iset) = icount(iset) + 1

   ELSE

      lbuf0da = lbuf0da + ndat0da(iset)
      lbuf2da = lbuf2da + ndat2da(iset)
      lbuf3da = lbuf3da + ndat3da(iset)
      lbuf0db = lbuf0db + ndat0db(iset)
      lbuf2db = lbuf2db + ndat2db(iset)
      lbuf3db = lbuf3db + ndat3db(iset)

   ENDIF

ENDDO iset_1000

CALL log_timer_out()


RETURN

END SUBROUTINE harmonic_analysis_update

!========================================================================

SUBROUTINE harmonic_analysis_data
!************************************************************************
!
! *harmonic_analysis_data* Evaluate and write harmonic data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Harmonic_Analysis.f90  V2.10.1
!
! Description -
!
! Reference -
!
! Calling program - harmonic_analysis
!
! External calls - astro_params
!
! Internal calls - ellips_params
!
! Module calls - add_secs_to_date, add_secs_to_phase, cholesky_solve,
!                combine_write_stats_loc, combine_write_submod,
!                complex_polar, diff_dates, error_alloc, loop_index,
!                write_vars
!
!************************************************************************
!
!1. Log info
!-----------
!

flag = .FALSE.
iset_110: DO iset=1,nosetsanal
   tlims = analgpars(iset)%tlims
   IF (loop_index(tlims,nt).AND.(nt.GT.tlims(1))) THEN
      flag = .TRUE.
   ENDIF
ENDDO iset_110
IF (flag) THEN
   procname(pglev+1) = 'harmonic_analysis_data'
   CALL log_timer_in()
ELSE
   RETURN
ENDIF

!
!2. Harmonic analysis and output
!-------------------------------
!

lbuf0da = 0; lbuf2da = 0; lbuf3da = 0
lbuf0db = 0; lbuf2db = 0; lbuf3db = 0

iset_2000: DO iset=1,nosetsanal

   tlims = analgpars(iset)%tlims
   IF (loop_index(tlims,nt).AND.(nt.GT.tlims(1))) THEN

!
!2.1 Initialise
!--------------
!
!2.1.1 Initialise parameters
!---------------------------
!

      gridded = analgpars(iset)%gridded
      packing = analgpars(iset)%packing
      ncout = analgpars(iset)%ncout
      nrout = analgpars(iset)%nrout
      nzout = analgpars(iset)%nzout
      nowetout = analgpars(iset)%nowetout
      nostats = analgpars(iset)%nostats
      nofreqs = nofreqsharm(iset)
      lprocs = limprocs(:,:,:,iset)
      ntout(iset) = ntout(iset) + 1
      fill_value = double_fill

      nanalvars0d = res0d(iset)%novars
      nanalvars2d = res2d(iset)%novars
      nanalvars3d = res3d(iset)%novars

      nocoords0d = MAX(res0d(iset)%nocoords,MAXVAL(amp0d(iset,:)%nocoords),&
                 & MAXVAL(pha0d(iset,:)%nocoords))
      nocoords2d = MAX(res2d(iset)%nocoords,MAXVAL(amp2d(iset,:)%nocoords),&
                & MAXVAL(pha2d(iset,:)%nocoords),MAXVAL(ell2d(iset,:)%nocoords))
      nocoords3d = MAX(res3d(iset)%nocoords,MAXVAL(amp3d(iset,:)%nocoords),&
                & MAXVAL(pha3d(iset,:)%nocoords),MAXVAL(ell3d(iset,:)%nocoords))
 
!
!2.1.2 Create mask array
!-----------------------
!

      IF (gridded) THEN
         ALLOCATE (maskvals(nclocout(iset),nrlocout(iset)),STAT=errstat)
         CALL error_alloc('maskvals',2,(/nclocout(iset),nrlocout(iset)/),&
                        & kndlog)
         IF (SIZE(maskvals).GT.0) THEN
            fmask = .TRUE.
            maskvals = mask_array(nclocout(iset),nrlocout(iset),&
                                & limloc(:,:,iset),'C')
         ELSE
            fmask = .FALSE.
         ENDIF
      ENDIF

!
!2.1.3 Reference phases and nodal factors
!----------------------------------------
!

      IF (cdate_time_ref(iset).EQ.cdatetime_undef) THEN
         IF (iopt_astro_anal.EQ.0) THEN
!           ---with respect to central time
            phase_ref = 0.0
         ELSEIF (iopt_astro_anal.EQ.1) THEN
!           ---with respect to astronomical phase at central time
            CALL add_secs_to_date(CDateTime,ctime,-tlims(3)/2,delt2d)
            CALL astro_params(index_anal,fnode_anal,phase_anal,nofreqsanal,&
                            & dlon_ref_anal,ctime)
         ENDIF
      ELSE
!        ---with respect to reference date
         CALL diff_dates(cdate_time_ref(iset),CDateTime,0,nosecs=nosecs)
         nosecs = tlims(3)/2 - nosecs
         nosecs2 = nosecs
         ifreq_213: DO ifreq=1,nofreqsanal
            CALL add_secs_to_phase(0.0,phase_ref(ifreq),nosecs2,1.0,&
                                 & harm_freq(ifreq))
         ENDDO ifreq_213
      ENDIF

!
!2.2 0-D output
!--------------
!

      IF (master.AND.nanalvars0d.GT.0) THEN
         outvars(1:nanalvars0d) = analvars(ivarsanal0d(iset,1:nanalvars0d)) 
         vecids(1:nanalvars0d) = (/(nocoords0d+ivar,ivar=1,nanalvars0d)/)

!
!2.2.1 Check whether phases are referenced with respect astronomical argument
!----------------------------------------------------------------------------
!

         indexarr(1:nofreqs) = ifreqsharm(iset,1:nofreqs)
         astro_ref = iopt_astro_anal.EQ.1.AND.&
                   & cdate_time_ref(iset).EQ.cdatetime_undef.AND.&
                   & ALL(index_anal(indexarr(1:nofreqs)).GT.0)

!
!2.2.2 Coefficients of harmonic expansion
!----------------------------------------
!

         ivar_222: DO ivar=1,nanalvars0d
            indexarr(0:nofreqs) = (/(lbuf0da+ivar+n*nanalvars0d,&
                                     & n=0,nofreqs)/)
            xvec(0:nofreqs) = buff0da(indexarr(0:nofreqs))
            CALL cholesky_solve(acoef(:,:,iset),xvec,nofreqs+1,&
                              & nofreqsanal+1,'acoef',info)
            buff0da(indexarr(0:nofreqs)) = xvec(0:nofreqs)
            indexarr(1:nofreqs) = (/(lbuf0db+ivar+(n-1)*nanalvars0d,&
                                     & n=1,nofreqs)/)
            xvec(1:nofreqs) = buff0db(indexarr(1:nofreqs))
            CALL cholesky_solve(bcoef(:,:,iset),xvec(1:nofreqsanal),&
                              & nofreqs,nofreqsanal,'bcoef',info)
            buff0db(indexarr(1:nofreqs)) = xvec(1:nofreqs)
         ENDDO ivar_222

!
!2.2.3 Residuals
!---------------
!

         IF (res0d(iset)%defined) THEN

!           ---allocate
            ALLOCATE (out1dres(nanalvars0d),STAT=errstat)
            CALL error_alloc('out1dres',1,(/nanalvars0d/),kndrtype)

!           ---evaluate
            out1dres = buff0da(lbuf0da+1:lbuf0da+nanalvars0d)

!           ---write
            CALL write_vars(out1dres,res0d(iset),0,&
                          & varatts=outvars(1:nanalvars0d),&
                          & vecids=vecids(1:nanalvars0d))

!           ---deallocate
            DEALLOCATE (out1dres)

         ENDIF

!
!2.2.4 Amplitudes and phases
!---------------------------
!

         n_224: DO n=1,nofreqs

            IF (amp0d(iset,n)%defined.OR.pha0d(iset,n)%defined) THEN

!              ---allocate
               ALLOCATE (out1damp(nanalvars0d),STAT=errstat)
               CALL error_alloc('out1damp',1,(/nanalvars0d/),kndrtype)
               ALLOCATE (out1dpha(nanalvars0d),STAT=errstat)
               CALL error_alloc('out1dpha',1,(/nanalvars0d/),kndrtype)

!              ---evaluate
               lbufa = lbuf0da + n*nanalvars0d
               lbufb = lbuf0db + (n-1)*nanalvars0d
               CALL complex_polar(buff0da(lbufa+1:lbufa+nanalvars0d),&
                                & buff0db(lbufb+1:lbufb+nanalvars0d),&
                                & xamp=out1damp,xpha=out1dpha) 

!              ---remove nodal factor
               ifreq = ifreqsharm(iset,n)
               IF (astro_ref) out1damp = out1damp/fnode_anal(ifreq)

!              ---add reference phase
               IF (astro_ref) THEN
                  out1dpha = out1dpha + phase_anal(ifreq)
               ELSE
                  out1dpha = out1dpha + phase_ref(ifreq)
               ENDIF
               out1dpha = MOD(out1dpha,twopi)
               IF (DegreesOut) out1dpha = out1dpha*radtodeg

!              ---write
               IF (amp0d(iset,n)%defined) THEN
                  CALL write_vars(out1damp,amp0d(iset,n),1,&
                                & varatts=outvars(1:nanalvars0d),&
                                & vecids=vecids(1:nanalvars0d))
               ENDIF
               IF (pha0d(iset,n)%defined) THEN
                  CALL write_vars(out1dpha,pha0d(iset,n),1,&
                                & varatts=outvars(1:nanalvars0d),&
                                & vecids=vecids(1:nanalvars0d))
               ENDIF

!              ---deallocate
               DEALLOCATE (out1damp,out1dpha)

            ENDIF
         ENDDO n_224
      ENDIF

!
!2.3 2-D output
!--------------
!

      IF (nanalvars2d.GT.0) THEN
         outvars(1:nanalvars2d) = &
                           & analvars(novars0d+ivarsanal2d(iset,1:nanalvars2d))
         vecids(1:nanalvars2d) = (/(nocoords2d+ivar,ivar=1,nanalvars2d)/)

!
!2.3.1 Check whether phases are referenced with respect astronomical argument
!----------------------------------------------------------------------------
!

         astro_ref = iopt_astro_anal.EQ.1.AND.&
                   & cdate_time_ref(iset).EQ.cdatetime_undef.AND.&
                   & ALL(index_anal(ifreqsharm(iset,1:nofreqs)).GT.0)

!
!2.3.2 Output dimensions
!-----------------------
!

         IF (gridded) THEN
            ndims = (/nclocout(iset),nrlocout(iset),1,nanalvars2d/)
         ELSE
            ndims = (/nostatsloc(iset),1,1,nanalvars2d/)
         ENDIF

!
!2.3.3 Coefficients of harmonic expansion
!----------------------------------------
!

         n1_233: DO n1=1,ndims(1)
         n2_233: DO n2=1,ndims(2)
         ivar_233: DO ivar=1,nanalvars2d
            indexarr(0:nofreqs) = (/(lbuf2da+ivar+nanalvars2d*(n+(nofreqs+1)*&
                                  & (n2-1+ndims(2)*(n1-1))),n=0,nofreqs)/)
            xvec(0:nofreqs) = buff2da(indexarr(0:nofreqs))
            CALL cholesky_solve(acoef(:,:,iset),xvec,nofreqs+1,&
                              & nofreqsanal+1,'acoef',info)
            buff2da(indexarr(0:nofreqs)) = xvec(0:nofreqs)
            indexarr(1:nofreqs) = (/(lbuf2db+ivar+nanalvars2d*(n-1+nofreqs*&
                                  & (n2-1+ndims(2)*(n1-1))),n=1,nofreqs)/)
            xvec(1:nofreqs) = buff2db(indexarr(1:nofreqs))
            CALL cholesky_solve(bcoef(:,:,iset),xvec(1:nofreqsanal),&
                              & nofreqs,nofreqsanal,'bcoef',info)
            buff2db(indexarr(1:nofreqs)) = xvec(1:nofreqs)
         ENDDO ivar_233
         ENDDO n2_233
         ENDDO n1_233

!
!2.3.4 Residuals
!---------------
!

         IF (res2d(iset)%defined) THEN

!
!2.3.4.1 Gridded
!---------------
!

            IF (gridded) THEN

!              ---allocate
               ALLOCATE(out3dres(ndims(1),ndims(2),nanalvars2d),STAT=errstat)
               CALL error_alloc('out3dres',3,(/ndims(1),ndims(2),nanalvars2d/),&
                              & kndrtype)

!              ---evaluate
               n1_23411: DO n1=1,ndims(1)
               n2_23411: DO n2=1,ndims(2)
                  lbufa = lbuf2da+nanalvars2d*(nofreqs+1)*(n2-1+ndims(2)*(n1-1))
                  out3dres(n1,n2,:) = buff2da(lbufa+1:lbufa+nanalvars2d)
               ENDDO n2_23411
               ENDDO n1_23411

!              ---insert fill values
               If (res2d(iset)%fill.AND.fmask) THEN
                  ivar_23412: DO ivar=1,nanalvars2d
                     WHERE (.NOT.maskvals)
                        out3dres(:,:,ivar) = fill_value
                     END WHERE
                  ENDDO ivar_23412
               ENDIF

!              ---write
               CALL combine_write_submod(out3dres,res2d(iset),0,&
                                      & (/ncout,nrout,nanalvars2d/),lprocs,&
                                      & nowet=nowetout,&
                                      & varatts=outvars(1:nanalvars2d),&
                                      & vecids=vecids(1:nanalvars2d))

!              ---deallocate
               DEALLOCATE (out3dres)

!
!2.3.4.2 Non-gridded
!-------------------
!

            ELSE

!              ---allocate
               ALLOCATE(out2dres(ndims(1),nanalvars2d),STAT=errstat)
               CALL error_alloc('out2dres',2,(/ndims(1),nanalvars2d/),&
                               & kndrtype)

!              ---evaluate
               n1_2342: DO n1=1,ndims(1)
                  lbufa = lbuf2da+nanalvars2d*(nofreqs+1)*(n1-1)
                  out2dres(n1,:) = buff2da(lbufa+1:lbufa+nanalvars2d)
               ENDDO n1_2342

!              ---write
               CALL combine_write_stats_loc(out2dres,res2d(iset),0,maxstats,&
                                      & nostats,nostatsprocs(:,iset),&
                                      & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                      & varatts=outvars(1:nanalvars2d),&
                                      & vecids=vecids(1:nanalvars2d))

!              ---deallocate
               DEALLOCATE (out2dres)

            ENDIF
         ENDIF

!
!2.3.5 Amplitudes and phases
!---------------------------
!

         n_235: DO n=1,nofreqs

            IF (amp2d(iset,n)%defined.OR.pha2d(iset,n)%defined) THEN

!
!2.3.5.1 Gridded
!---------------
!

               IF (gridded) THEN

!                 ---allocate
                  ALLOCATE (out3damp(ndims(1),ndims(2),nanalvars2d),&
                          & STAT=errstat)
                  CALL error_alloc('out3damp',3,&
                                & (/ndims(1),ndims(2),nanalvars2d/),kndrtype)
                  ALLOCATE (out3dpha(ndims(1),ndims(2),nanalvars2d),&
                          & STAT=errstat)
                  CALL error_alloc('out3dpha',3,&
                                & (/ndims(1),ndims(2),nanalvars2d/),kndrtype)

!                 ---evaluate
                  n1_23511: DO n1=1,ndims(1)
                  n2_23511: DO n2=1,ndims(2)
                     lbufa = lbuf2da+nanalvars2d*(n+(nofreqs+1)*(n2-1+ndims(2)*&
                                              & (n1-1)))
                     lbufb = lbuf2db+nanalvars2d*(n-1+nofreqs*(n2-1+ndims(2)*&
                                              & (n1-1)))
                     CALL complex_polar(buff2da(lbufa+1:lbufa+nanalvars2d),&
                                      & buff2db(lbufb+1:lbufb+nanalvars2d),&
                                      & xamp=out3damp(n1,n2,:),&
                                      & xpha=out3dpha(n1,n2,:))
                  ENDDO n2_23511
                  ENDDO n1_23511

!                 ---remove nodal factor
                  ifreq = ifreqsharm(iset,n)
                  IF (astro_ref) out3damp = out3damp/fnode_anal(ifreq)

!                 ---add reference phase
                  IF (astro_ref) THEN
                     out3dpha = out3dpha + phase_anal(ifreq)
                  ELSE
                     out3dpha = out3dpha + phase_ref(ifreq)
                  ENDIF
                  out3dpha = MOD(out3dpha,twopi)
                  IF (DegreesOut) out3dpha = out3dpha*radtodeg

!                 ---insert fill values
                  IF (amp2d(iset,n)%fill.AND.fmask) THEN
                     ivar_23512: DO ivar=1,nanalvars2d
                        WHERE (.NOT.maskvals)
                           out3damp(:,:,ivar) = fill_value
                        END WHERE
                     ENDDO ivar_23512
                  ENDIF
                  IF (pha2d(iset,n)%fill.AND.fmask) THEN
                     ivar_23513: DO ivar=1,nanalvars2d
                        WHERE (.NOT.maskvals)
                           out3dpha(:,:,ivar) = fill_value
                        END WHERE
                     ENDDO ivar_23513
                  ENDIF

!                 ---write
                  IF (amp2d(iset,n)%defined) THEN
                     CALL combine_write_submod(out3damp,amp2d(iset,n),0,&
                                            & (/ncout,nrout,nanalvars2d/),&
                                            & lprocs,nowet=nowetout,&
                                            & varatts=outvars(1:nanalvars2d),&
                                            & vecids=vecids(1:nanalvars2d))
                  ENDIF
                  IF (pha2d(iset,n)%defined) THEN
                     CALL combine_write_submod(out3dpha,pha2d(iset,n),0,&
                                            & (/ncout,nrout,nanalvars2d/),&
                                            & lprocs,nowet=nowetout,&
                                            & varatts=outvars(1:nanalvars2d),&
                                            & vecids=vecids(1:nanalvars2d))
                  ENDIF

!                  ---deallocate
                  DEALLOCATE (out3damp,out3dpha)

!
!2.3.5.2 Non-gridded
!-------------------
!

               ELSE

!                 ---allocate
                  ALLOCATE (out2damp(ndims(1),nanalvars2d),STAT=errstat)
                  CALL error_alloc('out2damp',2,(/ndims(1),nanalvars2d/),&
                                 & kndrtype)
                  ALLOCATE (out2dpha(ndims(1),nanalvars2d),STAT=errstat)
                  CALL error_alloc('out2dpha',2,(/ndims(1),nanalvars2d/),&
                                 & kndrtype)

!                 ---evaluate
                  n1_2352: DO n1=1,ndims(1)
                     lbufa = lbuf2da+nanalvars2d*(n+(nofreqs+1)*(n1-1))
                     lbufb = lbuf2db+nanalvars2d*(n-1+nofreqs*(n1-1))
                     CALL complex_polar(buff2da(lbufa+1:lbufa+nanalvars2d),&
                                      & buff2db(lbufb+1:lbufb+nanalvars2d),&
                                      & xamp=out2damp(n1,:),&
                                      & xpha=out2dpha(n1,:))
                  ENDDO n1_2352

!                 ---remove nodal factor
                  ifreq = ifreqsharm(iset,n)
                  IF (astro_ref) out2damp = out2damp/fnode_anal(ifreq)

!                 ---add reference phase
                  IF (astro_ref) THEN
                     out2dpha = out2dpha + phase_anal(ifreq)
                  ELSE
                     out2dpha = out2dpha + phase_ref(ifreq)
                  ENDIF
                  out2dpha = MOD(out2dpha,twopi)
                  IF (DegreesOut) out2dpha = out2dpha*radtodeg

!                 ---write
                  IF (amp2d(iset,n)%defined) THEN
                     CALL combine_write_stats_loc(out2damp,amp2d(iset,n),0,&
                                       & maxstats,nostats,nostatsprocs(:,iset),&
                                       & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                       & varatts=outvars(1:nanalvars2d),&
                                       & vecids=vecids(1:nanalvars2d))
                  ENDIF
                  IF (pha2d(iset,n)%defined) THEN
                     CALL combine_write_stats_loc(out2dpha,pha2d(iset,n),0,&
                                       & maxstats,nostats,nostatsprocs(:,iset),&
                                       & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                       & varatts=outvars(1:nanalvars2d),&
                                       & vecids=vecids(1:nanalvars2d))
                  ENDIF

!                 ---deallocate
                  DEALLOCATE (out2damp,out2dpha)
         
               ENDIF
            ENDIF
         ENDDO n_235
         
!
!2.3.6 Elliptic parameters
!-------------------------
!
!

         n_236: DO n=1,nofreqs
            IF (ell2d(iset,n)%defined) THEN

               l1 = ivecell2d(iset,1); l2 = ivecell2d(iset,2)
               novars = ell2d(iset,n)%novars
               outvars(1:novars) = ellvars(ivarsell2d(iset,1:novars))
               vecids(1:novars) = (/(nocoords2d+ivar,ivar=1,novars)/)

!
!2.3.6.1 Gridded
!----------------
!

               IF (gridded) THEN

!                  ---allocate
                  ALLOCATE (ell3dsub(ndims(1),ndims(2),novars),STAT=errstat)
                  CALL error_alloc('ell3dsub',3,(/ndims(1),ndims(2),novars/),&
                                 & kndrtype)

!                 ---evaluate
                  n1_23611: DO n1=1,ndims(1)
                  n2_23611: DO n2=1,ndims(2)
                     lbufa1 = lbuf2da+l1+nanalvars2d*(n+(nofreqs+1)*&
                                                  & (n2-1+ndims(2)*(n1-1)))
                     lbufa2 = lbufa1+l2-l1 
                     lbufb1 = lbuf2db+l1+nanalvars2d*(n-1+nofreqs*&
                                                  & (n2-1+ndims(2)*(n1-1)))
                     lbufb2 = lbufb1+l2-l1
                     outell = ellips_params(buff2da(lbufa1),buff2da(lbufa2),&
                                          & buff2db(lbufb1),buff2db(lbufb2))
                     ell3dsub(n1,n2,1:novars) = &
                   & outell(ivarsell2d(iset,1:novars))
                  ENDDO n2_23611
                  ENDDO n1_23611

!                 ---insert fill values
                  IF (ell2d(iset,n)%fill.AND.fmask) THEN
                     ivar_23612: DO ivar=1,novars
                        WHERE (.NOT.maskvals)
                           ell3dsub(:,:,ivar) = fill_value
                        END WHERE
                     ENDDO ivar_23612
                  ENDIF

!                 ---write
                  CALL combine_write_submod(ell3dsub,ell2d(iset,n),0,&
                                         & (/ncout,nrout,novars/),lprocs,&
                                         & nowet=nowetout,&
                                         & varatts=outvars(1:novars),&
                                         & vecids=vecids(1:novars))

!                 ---deallocate
                  DEALLOCATE (ell3dsub)

!
!2.3.5.2 Non-gridded
!-------------------
!

               ELSE

!                 ---allocate
                  ALLOCATE (ell2dsub(ndims(1),novars),STAT=errstat)
                  CALL error_alloc('ell2dsub',2,(/ndims(1),novars/),&
                                 & kndrtype)

!                 ---evaluate
                  n1_2361: DO n1=1,ndims(1)
                     lbufa1 = lbuf2da+l1+nanalvars2d*(n+(nofreqs+1)*(n1-1))
                     lbufa2 = lbufa1+l2-l1 
                     lbufb1 = lbuf2db+l1+nanalvars2d*(n-1+nofreqs*(n1-1))
                     lbufb2 = lbufb1+l2-l1
                     outell = ellips_params(buff2da(lbufa1),buff2da(lbufa2),&
                                          & buff2db(lbufb1),buff2db(lbufb2))
                     ell2dsub(n1,1:novars) = outell(ivarsell2d(iset,1:novars))
                  ENDDO n1_2361

!                 ---write
                  CALL combine_write_stats_loc(ell2dsub,ell2d(iset,n),0,&
                                       & maxstats,nostats,nostatsprocs(:,iset),&
                                       & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                       & varatts=outvars(1:novars),&
                                       & vecids=vecids(1:novars))

!                 ---deallocate
                  DEALLOCATE (ell2dsub)

               ENDIF
            ENDIF
         ENDDO n_236
      ENDIF

!
!2.4 3-D output
!--------------
!

      IF (nanalvars3d.GT.0) THEN
         outvars(1:nanalvars3d) = &
                   & analvars(novars0d+novars2d+ivarsanal3d(iset,1:nanalvars3d))
         vecids(1:nanalvars3d) = (/(nocoords3d+ivar,ivar=1,nanalvars3d)/)

!
!2.4.1 Check whether phases are referenced with respect astronomical argument
!----------------------------------------------------------------------------
!

         astro_ref = iopt_astro_anal.EQ.1.AND.&
                   & cdate_time_ref(iset).EQ.cdatetime_undef.AND.&
                   & ALL(index_anal(ifreqsharm(iset,1:nofreqs)).GT.0)

!
!2.4.2 Output dimensions
!-----------------------
!

         IF (gridded) THEN
            ndims = (/nclocout(iset),nrlocout(iset),nzout,nanalvars3d/)
         ELSE
            ndims = (/nostatsloc(iset),1,nzout,nanalvars3d/)
         ENDIF

!
!2.4.3 Coefficients of harmonic expansion
!----------------------------------------
!

         n1_243: DO n1=1,ndims(1)
         n2_243: DO n2=1,ndims(2)
         n3_243: DO n3=1,ndims(3)
         ivar_243: DO ivar=1,nanalvars3d
            indexarr(0:nofreqs) = (/(lbuf3da+ivar+nanalvars3d*(n+(nofreqs+1)*&
                                  & (n3-1+ndims(3)*(n2-1+ndims(2)*(n1-1)))),&
                                  & n=0,nofreqs)/)
            xvec(0:nofreqs) = buff3da(indexarr(0:nofreqs))
            CALL cholesky_solve(acoef(:,:,iset),xvec,nofreqs+1,&
                              & nofreqsanal+1,'acoef',info)
            buff3da(indexarr(0:nofreqs)) = xvec(0:nofreqs)
            indexarr(1:nofreqs) = (/(lbuf3db+ivar+nanalvars3d*(n-1+nofreqs*&
                                  & (n3-1+ndims(3)*(n2-1+ndims(2)*(n1-1)))),&
                                  & n=1,nofreqs)/)
            xvec(1:nofreqs) = buff3db(indexarr(1:nofreqs))
            CALL cholesky_solve(bcoef(:,:,iset),xvec(1:nofreqsanal),&
                              & nofreqs,nofreqsanal,'bcoef',info)
            buff3db(indexarr(1:nofreqs)) = xvec(1:nofreqs)
         ENDDO ivar_243
         ENDDO n3_243
         ENDDO n2_243
         ENDDO n1_243

!
!2.4.4 Residuals
!---------------
!

         IF (res3d(iset)%defined) THEN

!
!2.4.4.1 Gridded
!---------------
!

            IF (gridded) THEN

!              ---allocate
               ALLOCATE(out4dres(ndims(1),ndims(2),ndims(3),nanalvars3d),&
                      & STAT=errstat)
               CALL error_alloc('out4dres',4,ndims,kndrtype)

!              ---evaluate
               n1_24411: DO n1=1,ndims(1)
               n2_24411: DO n2=1,ndims(2)
               n3_24411: DO n3=1,ndims(3)
                  lbufa = lbuf3da+nanalvars3d*(nofreqs+1)*(n3-1+ndims(3)*&
                                                      & (n2-1+ndims(2)*(n1-1)))
                  out4dres(n1,n2,n3,:) = buff3da(lbufa+1:lbufa+nanalvars3d)
               ENDDO n3_24411
               ENDDO n2_24411
               ENDDO n1_24411

!              ---inseert fill values
               IF (res3d(iset)%fill.AND.fmask) THEN
                  ivar_24412: DO ivar=1,nanalvars3d
                  n3_24412: do n3=1,ndims(3)
                     WHERE (.NOT.maskvals)
                        out4dres(:,:,n3,ivar) = fill_value
                     END WHERE
                  ENDDO n3_24412
                  ENDDO ivar_24412
               ENDIF

!              ---write
               CALL combine_write_submod(out4dres,res3d(iset),0,&
                                      & (/ncout,nrout,nzout,nanalvars3d/),&
                                      & lprocs,nowet=nowetout,&
                                      & varatts=outvars(1:nanalvars3d),&
                                      & vecids=vecids(1:nanalvars3d))

!              ---deallocate
               DEALLOCATE (out4dres) 

!
!2.4.4.2 Non-gridded
!-------------------
!

            ELSE

!              ---allocate
               ALLOCATE(out3dres(ndims(1),ndims(3),nanalvars3d),STAT=errstat)
               CALL error_alloc('out3dres',3,(/ndims(1),ndims(3),nanalvars3d/),&
                               & kndrtype)

!              ---evaluate
               n1_24421: DO n1=1,ndims(1)
               n3_24421: DO n3=1,ndims(3)
                  lbufa = lbuf3da+nanalvars3d*(nofreqs+1)*(n3-1+ndims(3)*(n1-1))
                  out3dres(n1,n3,:) = buff3da(lbufa+1:lbufa+nanalvars3d)
               ENDDO n3_24421
               ENDDO n1_24421

!              ---write
               CALL combine_write_stats_loc(out3dres,res3d(iset),0,maxstats,&
                                      & nostats,nostatsprocs(:,iset),&
                                      & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                      & varatts=outvars(1:nanalvars3d),&
                                      & vecids=vecids(1:nanalvars3d))

!              ---deallocate
               DEALLOCATE (out3dres) 

            ENDIF
         ENDIF

!
!2.4.5 Amplitudes and phases
!---------------------------
!

         n_245: DO n=1,nofreqs

            IF (amp3d(iset,n)%defined.OR.pha3d(iset,n)%defined) THEN

!
!2.4.5.1 Gridded
!---------------
!

               IF (gridded) THEN

!                 ---allocate
                  ALLOCATE (out4damp(ndims(1),ndims(2),ndims(3),nanalvars3d),&
                          & STAT=errstat)
                  CALL error_alloc('out4damp',4,ndims,kndrtype)
                  ALLOCATE (out4dpha(ndims(1),ndims(2),ndims(3),nanalvars3d),&
                          & STAT=errstat)
                  CALL error_alloc('out4dpha',4,ndims,kndrtype)

!                 ---evaluate
                  n1_24511: DO n1=1,ndims(1)
                  n2_24511: DO n2=1,ndims(2)
                  n3_24511: DO n3=1,ndims(3)
                     lbufa = lbuf3da+nanalvars3d*(n+(nofreqs+1)*(n3-1+ndims(3)*&
                                              & (n2-1+ndims(2)*(n1-1))))
                     lbufb = lbuf3db+nanalvars3d*(n-1+nofreqs*(n3-1+ndims(3)*&
                                              & (n2-1+ndims(2)*(n1-1))))
                     CALL complex_polar(buff3da(lbufa+1:lbufa+nanalvars3d),&
                                      & buff3db(lbufb+1:lbufb+nanalvars3d),&
                                      & xamp=out4damp(n1,n2,n3,:),&
                                      & xpha=out4dpha(n1,n2,n3,:))
                  ENDDO n3_24511
                  ENDDO n2_24511
                  ENDDO n1_24511

!                 ---remove nodal factor
                  ifreq = ifreqsharm(iset,n)
                  IF (astro_ref) out4damp = out4damp/fnode_anal(ifreq)

!                 ---add reference phase
                  IF (astro_ref) THEN
                     out4dpha = out4dpha + phase_anal(ifreq)
                  ELSE
                     out4dpha = out4dpha + phase_ref(ifreq)
                  ENDIF
                  out4dpha = MOD(out4dpha,twopi)
                  IF (DegreesOut) out4dpha = out4dpha*radtodeg

!                 ---insert fill values
                  IF (amp3d(iset,n)%fill.AND.fmask) THEN
                     ivar_24512: DO ivar=1,nanalvars3d
                     n3_24512: DO n3=1,ndims(3)
                        WHERE (.NOT.maskvals)
                           out4damp(:,:,n3,ivar) = fill_value
                        END WHERE
                     ENDDO n3_24512
                     ENDDO ivar_24512
                  ENDIF
                  IF (pha3d(iset,n)%fill.AND.fmask) THEN
                     ivar_24513: DO ivar=1,nanalvars3d
                     n3_24513: DO n3=1,ndims(3)
                        WHERE (.NOT.maskvals)
                           out4dpha(:,:,n3,ivar) = fill_value
                        END WHERE
                     ENDDO n3_24513
                     ENDDO ivar_24513
                  ENDIF
!                 ---write
                  IF (amp3d(iset,n)%defined) THEN
                     CALL combine_write_submod(out4damp,amp3d(iset,n),0,&
                                           & (/ncout,nrout,nzout,nanalvars3d/),&
                                           & lprocs,nowet=nowetout,&
                                           & varatts=outvars(1:nanalvars3d),&
                                           & vecids=vecids(1:nanalvars3d))
                  ENDIF
                  IF (pha3d(iset,n)%defined) THEN
                     CALL combine_write_submod(out4dpha,pha3d(iset,n),0,&
                                           & (/ncout,nrout,nzout,nanalvars3d/),&
                                           & lprocs,nowet=nowetout,&
                                           & varatts=outvars(1:nanalvars3d),&
                                           & vecids=vecids(1:nanalvars3d))
                  ENDIF

!                 ---deallocate
                  DEALLOCATE (out4damp,out4dpha)

!
!2.4.5.2 Non-gridded
!-------------------
!

               ELSE

!                 ---allocate
                  ALLOCATE (out3damp(ndims(1),ndims(3),nanalvars3d),&
                          & STAT=errstat)
                  CALL error_alloc('out3damp',3,&
                                & (/ndims(1),ndims(3),nanalvars3d/),kndrtype)
                  ALLOCATE (out3dpha(ndims(1),ndims(3),nanalvars3d),&
                          & STAT=errstat)
                  CALL error_alloc('out3dpha',3,&
                                & (/ndims(1),ndims(3),nanalvars3d/),kndrtype)

!                 ---evaluate
                  n1_2452: DO n1=1,ndims(1)
                  n3_2452: DO n3=1,ndims(3)
                     lbufa = lbuf3da+nanalvars3d*(n+(nofreqs+1)*(n3-1+ndims(3)*&
                                              & (n1-1)))
                     lbufb = lbuf3db+nanalvars3d*(n-1+nofreqs*(n3-1+ndims(3)*&
                                              & (n1-1)))
                     CALL complex_polar(buff3da(lbufa+1:lbufa+nanalvars3d),&
                                      & buff3db(lbufb+1:lbufb+nanalvars3d),&
                                      & xamp=out3damp(n1,n3,:),&
                                      & xpha=out3dpha(n1,n3,:))
                  ENDDO n3_2452
                  ENDDO n1_2452

!                 ---remove nodal factor
                  ifreq = ifreqsharm(iset,n)
                  IF (astro_ref) out3damp = out3damp/fnode_anal(ifreq)

!                 ---add reference phase
                  IF (astro_ref) THEN
                     out3dpha = out3dpha + phase_anal(ifreq)
                  ELSE
                     out3dpha = out3dpha + phase_ref(ifreq)
                  ENDIF
                  out3dpha = MOD(out3dpha,twopi)
                  IF (DegreesOut) out3dpha = out3dpha*radtodeg

!                 ---write
                  IF (amp3d(iset,n)%defined) THEN
                     CALL combine_write_stats_loc(out3damp,amp3d(iset,n),0,&
                                       & maxstats,nostats,nostatsprocs(:,iset),&
                                       & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                       & varatts=outvars(1:nanalvars3d),&
                                       & vecids=vecids(1:nanalvars3d))
                  ENDIF
                  IF (pha3d(iset,n)%defined) THEN
                     CALL combine_write_stats_loc(out3dpha,pha3d(iset,n),0,&
                                       & maxstats,nostats,nostatsprocs(:,iset),&
                                       & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                       & varatts=outvars(1:nanalvars3d),&
                                       & vecids=vecids(1:nanalvars3d))
                  ENDIF

!                 ---deallocate
                  DEALLOCATE (out3damp,out3dpha)
               
               ENDIF
            ENDIF
         ENDDO n_245

!
!2.4.6 Elliptic parameters
!-------------------------
!
!

         n_246: DO n=1,nofreqs
            IF (ell3d(iset,n)%defined) THEN
               l1 = ivecell3d(iset,1); l2 = ivecell3d(iset,2)
               novars = ell3d(iset,n)%novars
               outvars(1:novars) = ellvars(ivarsell3d(iset,1:novars))
               vecids(1:novars) = (/(nocoords3d+ivar,ivar=1,novars)/)

!
!2.4.6.1 Gridded
!---------------
!

               IF (gridded) THEN

!                 ---allocate
                  ALLOCATE (ell4dsub(ndims(1),ndims(2),ndims(3),novars),&
                       & STAT=errstat)
                  CALL error_alloc('ell4dsub',4,&
                                & (/ndims(1),ndims(2),ndims(3),novars/),&
                                & kndrtype)

!                 ---evaluate
                  n1_24611: DO n1=1,ndims(1)
                  n2_24611: DO n2=1,ndims(2)
                  n3_24611: DO n3=1,ndims(3)
                     lbufa1 = lbuf3da+l1+nanalvars3d*(n+(nofreqs+1)*&
                                                  & (n3-1+ndims(3)*&
                                                  & (n2-1+ndims(2)*(n1-1))))
                     lbufa2 = lbufa1+l2-l1 
                     lbufb1 = lbuf3db+l1+nanalvars3d*(n-1+nofreqs*&
                                                  & (n3-1+ndims(3)*&
                                                  & (n2-1+ndims(2)*(n1-1))))
                     lbufb2 = lbufb1+l2-l1
                     outell = ellips_params(buff3da(lbufa1),buff3da(lbufa2),&
                                          & buff3db(lbufb1),buff3db(lbufb2))
                     ell4dsub(n1,n2,n3,1:novars) = &
                   & outell(ivarsell3d(iset,1:novars))
                  ENDDO n3_24611
                  ENDDO n2_24611
                  ENDDO n1_24611

!                 ---insert fill values
                  IF (ell3d(iset,n)%fill.AND.fmask) THEN
                     ivar_24612: DO ivar=1,novars
                     n3_24612: DO n3=1,ndims(3)
                        WHERE (.NOT.maskvals)
                           ell4dsub(:,:,n3,ivar) = fill_value
                        END WHERE
                     ENDDO n3_24612
                     ENDDO ivar_24612
                  ENDIF

!                 ---write
                  CALL combine_write_submod(ell4dsub,ell3d(iset,n),0,&
                                        & (/ncout,nrout,nzout,novars/),lprocs,&
                                        & nowet=nowetout,&
                                        & varatts=outvars(1:novars),&
                                        & vecids=vecids(1:novars))

!                 ---deallocate
                  DEALLOCATE (ell4dsub)

!
!2.4.6.2 Non-gridded
!-------------------
!
                  
               ELSE

!                 ---allocate
                  ALLOCATE (ell3dsub(ndims(1),ndims(3),novars),&
                       & STAT=errstat)
                  CALL error_alloc('ell3dsub',3,&
                                & (/ndims(1),ndims(3),novars/),kndrtype)

!                 ---evaluate
                  n1_2462: DO n1=1,ndims(1)
                  n3_2462: DO n3=1,ndims(3)
                     lbufa1 = lbuf3da+l1+nanalvars3d*(n+(nofreqs+1)*&
                                                  & (n3-1+ndims(3)*(n1-1)))
                     lbufa2 = lbufa1+l2-l1 
                     lbufb1 = lbuf3db+l1+nanalvars3d*(n-1+nofreqs*&
                                                  & (n3-1+ndims(3)*(n1-1)))
                     lbufb2 = lbufb1+l2-l1
                     outell = ellips_params(buff3da(lbufa1),buff3da(lbufa2),&
                                          & buff3db(lbufb1),buff3db(lbufb2))
                     ell3dsub(n1,n3,1:novars) = &
                   & outell(ivarsell3d(iset,1:novars))
                  ENDDO n3_2462
                  ENDDO n1_2462

!                 ---write
                  CALL combine_write_stats_loc(ell3dsub,ell3d(iset,n),0,&
                                       & maxstats,nostats,nostatsprocs(:,iset),&
                                       & lstatsprocs(:,:,iset),reduced=.FALSE.,&
                                       & varatts=outvars(1:novars),vecids=vecids(1:novars))

!                 ---deallocate
                  DEALLOCATE (ell3dsub)

               ENDIF
            ENDIF
         ENDDO n_246

      ENDIF

!     ---deallocate
      IF (ALLOCATED(maskvals)) DEALLOCATE (maskvals)

   ENDIF

   lbuf0da = lbuf0da + ndat0da(iset); lbuf0db = lbuf0db + ndat0db(iset)
   lbuf2da = lbuf2da + ndat2da(iset); lbuf2db = lbuf2db + ndat2db(iset)
   lbuf3da = lbuf3da + ndat3da(iset); lbuf3db = lbuf3db + ndat3db(iset)

ENDDO iset_2000

CALL log_timer_out()


RETURN

END SUBROUTINE harmonic_analysis_data

!========================================================================

FUNCTION ellips_params(uacoef,vacoef,ubcoef,vbcoef)
!************************************************************************
!
! *ellips_params* Tidal ellipse parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Harmonic_Analysis.f90  V2.0
!
! Description - Returns parameters as a 7-element vector
!                1: major axis
!                2: minor axis
!                3: ellipticity
!                4: inclination
!                5: elliptic angle
!                6: amplitude of cyclonic current
!                7: amplitude of anticyclonic current
!
! Reference -
!
! Calling program - harmonic_analysis_data
!
! Module calls - complex_polar
!
!************************************************************************
!
!* Arguments
!
REAL, INTENT(IN) :: uacoef, ubcoef, vacoef, vbcoef
REAL, DIMENSION(7) :: ellips_params

!
! Name     Type  Purpose
!-----------------------------------------------------------------------------
!*uacoef*  REAL  Coefficient of cos-term in harmonic expansion for u-current
!*ubcoef*  REAL  Coefficient of sin-term in harmonic expansion for u-current
!*vacoef*  REAL  Coefficient of cos-term in harmonic expansion for v-current
!*vbcoef*  REAL  Coefficient of sin-term in harmonic expansion for v-current
!
!------------------------------------------------------------------------------
!


rplusr = 0.5*(uacoef+vbcoef)
rplusi = 0.5*(vacoef-ubcoef)
rminr = 0.5*(uacoef-vbcoef)
rmini = -0.5*(vacoef+ubcoef)
CALL complex_polar(rplusr,rplusi,xamp=rplus,xpha=aplus)
CALL complex_polar(rminr,rmini,xamp=rmin,xpha=amin)
theta = 0.5*(aplus-amin)
IF (theta.LT.0.0) theta = theta + pi
phi = pi - 0.5*(aplus+amin)
IF (phi.LT.0.0) phi = phi + pi
IF (phi.EQ.pi) phi = 0.0
ellips_params(1) = rplus+rmin
ellips_params(2) = ABS(rplus-rmin)
IF ((rplus+rmin).GT.epstol) THEN
   ellips_params(3) = (rplus-rmin)/(rplus+rmin)
ELSE
   ellips_params(3) = 0.0
ENDIF
ellips_params(4) = theta
ellips_params(5) = phi
ellips_params(6) = rplus
ellips_params(7) = rmin
IF (DeGreesOut) THEN
   ellips_params(4) = ellips_params(4)*radtodeg
   ellips_params(5) = ellips_params(5)*radtodeg
ENDIF


RETURN

END FUNCTION ellips_params


END SUBROUTINE harmonic_analysis
