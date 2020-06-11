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
! *Surface_Boundary_Data_1D* Ensemble of routines for 1-D surface forcing
!                            (slopes and elevation) conditions and data 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Boundary_Data_1D.f90  V2.11
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - 
!
! Reference -
!
! Routines - define_1dsur_data, define_1dsur_spec, read_1dsur_data,
!            read_1dsur_spec, update_1dsur_data, write_1dsur_data,
!            write_1dsur_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE define_1dsur_data(ciodatetime,data1d,novars)
!************************************************************************
!
! *define_1dsur_data* Obtain data for 1-D forcing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Boundary_Data_1D.f90  V2.1.0
!
! Description - 
!
! Reference -
!
! Calling program - update_1dsur_data
!
! External calls - read_1dsur_data, usrdef_1dsur_data, write_1dsur_data
!
! Module calls - error_file
!
!************************************************************************
!
USE iopars
USE syspars
USE error_routines, ONLY: error_file
USE time_routines, ONLY: log_timer_in, log_timer_out, suspend_proc

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: novars
REAL, INTENT(INOUT), DIMENSION(novars) :: data1d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR    Date/time in data file
!*data1d*      REAL    Surface data
!*novars*      INTEGER Number of surface data
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: define_data
INTEGER :: endfile, nosecs


procname(pglev+1) = 'define_1dsur_data'
CALL log_timer_in()

!
!1. Initialise
!-------------
!

endfile = modfiles(io_1uvsur,2,1)%endfile
nosecs = 0
data1d = 0.0

!
!2. Open file on first call
!--------------------------
!

DO WHILE (modfiles(io_1uvsur,2,1)%iostat.EQ.0)

!  ---standard
   IF (modfiles(io_1uvsur,2,1)%status.EQ.'R') THEN
      CALL read_1dsur_data(ciodatetime,data1d,novars)
!  ---user-defined
   ELSEIF (modfiles(io_1uvsur,2,1)%status.EQ.'N') THEN
      CALL usrdef_1dsur_data(ciodatetime,data1d,novars)
   ENDIF

!  ---suspend/abort if needed
   IF (modfiles(io_1uvsur,2,1)%iostat.EQ.-1) THEN
      SELECT CASE (endfile)
         CASE (0:1)
            CALL error_file(ierrno_fopen,filepars=modfiles(io_1uvsur,2,1))
         CASE (2)
            CALL suspend_proc(nowaitsecs)
            nosecs = nosecs + nowaitsecs
            IF (nosecs.GT.maxwaitsecs) THEN
               CALL error_file(ierrno_fopen,filepars=modfiles(io_1uvsur,2,1))
            ENDIF
      END SELECT
   ENDIF
         
ENDDO

!
!3. Obtain data
!--------------
!

define_data = .TRUE.

DO WHILE (define_data)

!  ---read
   IF (modfiles(io_1uvsur,2,1)%status.EQ.'R') THEN
      CALL read_1dsur_data(ciodatetime,data1d,novars)
!  ---user-defined
   ELSEIF (modfiles(io_1uvsur,2,1)%status.EQ.'N') THEN
      CALL usrdef_1dsur_data(ciodatetime,data1d,novars)
   ENDIF
   IF (modfiles(io_1uvsur,2,1)%iostat.EQ.3) ciodatetime = cdatetime_undef

!  ---suspend/abort if needed
   IF (modfiles(io_1uvsur,2,1)%iostat.EQ.3) THEN
      SELECT CASE (endfile)
         CASE (0)
            CALL error_file(ierrno_fend,filepars=modfiles(io_1uvsur,2,1))
         CASE (1)
            define_data = .FALSE.
         CASE (2)
            CALL suspend_proc(nowaitsecs)
            nosecs = nosecs + nowaitsecs
            IF (nosecs.GT.maxwaitsecs) THEN
               CALL error_file(ierrno_fend,filepars=modfiles(io_1uvsur,2,1))
            ENDIF
      END SELECT
   ELSE
      define_data = .FALSE.
   ENDIF

ENDDO
 
!
!4. Write data
!-------------
!

IF (modfiles(io_1uvsur,2,2)%defined.AND.&
  & modfiles(io_1uvsur,2,1)%iostat.LT.3) THEN
   CALL write_1dsur_data(ciodatetime,data1d,novars)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE define_1dsur_data

!========================================================================

SUBROUTINE define_1dsur_spec
!************************************************************************
!
! *define_1dsur_spec* Define specifier arrays for 1-D forcing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Boundary_Data_1D.f90  V2.1.0
!
! Description - 
!
! Reference -
!
! Calling program - pressure_gradient_1d
!
! External calls - read_1dsur_spec, usrdef_1dsur_spec, write_1dsur_spec
!
! Module calls - error_abort, error_alloc, error_limits_var
!
!************************************************************************
!
USE iopars
USE obconds
USE paralpars
USE switches
USE syspars
USE tide
USE error_routines, ONLY: error_abort, error_alloc, error_limits_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE


procname(pglev+1) = 'define_1dsur_spec'
CALL log_timer_in()

!
!1. Allocate specifier arrays
!----------------------------
!
!---pressure gradient
ALLOCATE (gxslope_amp(nconobc),STAT=errstat)
CALL error_alloc('gxslope_amp',1,(/nconobc/),kndrtype)
IF (nconobc.GT.0) gxslope_amp = 0.0
ALLOCATE (gyslope_amp(nconobc),STAT=errstat)
CALL error_alloc('gyslope_amp',1,(/nconobc/),kndrtype)
IF (nconobc.GT.0) gyslope_amp = 0.0
ALLOCATE (gxslope_pha(nconobc),STAT=errstat)
CALL error_alloc('gxslope_pha',1,(/nconobc/),kndrtype)
IF (nconobc.GT.0) gxslope_pha = 0.0
ALLOCATE (gyslope_pha(nconobc),STAT=errstat)
CALL error_alloc('gyslope_pha',1,(/nconobc/),kndrtype)
IF (nconobc.GT.0) gyslope_pha = 0.0

!---surface elevations
ALLOCATE (zeta_amp(nconobc),STAT=errstat)
CALL error_alloc('zeta_amp',1,(/nconobc/),kndrtype)
IF (nconobc.GT.0) zeta_amp = 0.0
ALLOCATE (zeta_pha(nconobc),STAT=errstat)
CALL error_alloc('zeta_pha',1,(/nconobc/),kndrtype)
IF (nconobc.GT.0) zeta_pha = 0.0

!---type of data
isur1dtype = 0

!
!2. Obtain specifier arrays
!--------------------------
!
!---read
IF (modfiles(io_1uvsur,1,1)%status.EQ.'R') THEN
   CALL read_1dsur_spec
!---user-defined
ELSEIF (modfiles(io_1uvsur,1,1)%status.EQ.'N') THEN
   CALL usrdef_1dsur_spec
ENDIF

!---write
IF (modfiles(io_1uvsur,1,2)%defined) CALL write_1dsur_spec

!
!3. Check type of data variables
!-------------------------------
!

IF (master.AND.modfiles(io_1uvsur,2,1)%defined) THEN
   CALL error_limits_var(isur1dtype,'isur1dtype',1,3)
ENDIF

CALL error_abort('define_1dsur_spec',ierrno_input)

CALL log_timer_out()


RETURN

END SUBROUTINE define_1dsur_spec

!========================================================================

SUBROUTINE read_1dsur_data(ciodatetime,data1d,novars)
!************************************************************************
!
! *read_1dsur_data* Read 1-D forcing data in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Boundary_Data_1D.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - define_1dsur_data
!
! Module calls - error_abort, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_time, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE syspars
USE error_routines, ONLY: error_abort, error_alloc_struc
USE inout_routines, ONLY: open_filepars, read_glbatts_mod, read_time, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: novars
REAL, INTENT(INOUT), DIMENSION(novars) :: data1d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR    Date/time in data file
!*data1d*      REAL    Surface data
!*novars*      INTEGER Number of surface data
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cval
INTEGER :: ivar, numvars
TYPE (FileParams) :: filepars
INTEGER, DIMENSION(novars) :: vecids
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_1dsur_data'
CALL log_timer_in()

filepars = modfiles(io_1uvsur,2,1)
numvars = novars + 1

!
!1. File header on first call
!----------------------------
!

IF (filepars%iostat.EQ.0) THEN

!  ---open file
   filepars = modfiles(io_1uvsur,2,1)
   CALL open_filepars(filepars)
   IF (filepars%endfile.EQ.2.AND.filepars%iostat.EQ.-1) GOTO 1000

!  ---global attributes
   CALL read_glbatts_mod(filepars)
   IF (filepars%novars.NE.novars) THEN
      nerrs = 1
      IF (master.AND.errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval,'(I12)') filepars%novars; cval = ADJUSTL(cval)
         WRITE (ioerr,'(A)') 'Wrong number of variables in file '&
                            & //TRIM(filepars%filename)//': '//TRIM(cval)
         WRITE (cval,'(I12)') novars; cval = ADJUSTL(cval)
         WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cval)
      ENDIF
   ENDIF

!  ---allocate
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
   
!  ---variable attributes
   CALL read_varatts_mod(filepars,varatts,numvars)

   CALL error_abort('read_1dsur_data',ierrno_input)

!  ---deallocate   
   DEALLOCATE (varatts)

   GOTO 1000
   
ENDIF

!
!2. Read data
!------------
!
!2.1 Date/time
!-------------
!

CALL read_time(ciodatetime,filepars)

!
!2.2 Data
!--------
!

IF (filepars%iostat.LT.3) THEN
   vecids = (/(ivar,ivar=2,numvars)/)
   CALL read_vars(data1d,filepars,0,vecids=vecids)
ENDIF

1000 modfiles(io_1uvsur,2,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_1dsur_data

!========================================================================

SUBROUTINE read_1dsur_spec
!************************************************************************
!
! *real_1dsur_spec* Read specifier arrays for 1-D forcing in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Boundary_Data_1D.f90  V2.1.0
!
! Description - 
!
! Reference -
!
! Calling program - define_1dsur_spec
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars,
!
!************************************************************************
!
USE datatypes
USE iopars
USE obconds
USE switches
USE tide
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_1dsur_spec'
CALL log_timer_in()

filepars = modfiles(io_1uvsur,1,1)

!
!1. File header
!--------------
!

CALL open_filepars(filepars)
CALL read_glbatts_mod(filepars)
numvars = filepars%novars
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL read_varatts_mod(filepars,varatts,numvars)

!
!2. Read data
!------------
!

!---type of data
varid = 1
CALL read_vars(isur1dtype,filepars,varid,varatts=varatts)

!---amplitudes and phases of X-component of pressure gradient
IF (isur1dtype.EQ.1.OR.isur1dtype.EQ.3) THEN
   varid = varid + 1
   CALL read_vars(gxslope_amp,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(gxslope_pha,filepars,varid,varatts=varatts)
ENDIF

!---amplitudes and phases of Y-component of pressure gradient
IF (isur1dtype.EQ.1.OR.isur1dtype.EQ.3) THEN
   varid = varid + 1
   CALL read_vars(gyslope_amp,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(gyslope_pha,filepars,varid,varatts=varatts)
ENDIF

!---surface elevation
IF (isur1dtype.EQ.1.OR.isur1dtype.EQ.2) THEN
   varid = varid + 1
   CALL read_vars(zeta_amp,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL read_vars(zeta_pha,filepars,varid,varatts=varatts)
ENDIF

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_1uvsur,1,1) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE read_1dsur_spec

!========================================================================

SUBROUTINE update_1dsur_data
!************************************************************************
!
! *update_1dsur_data* Update surface slopes and elevation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Boundary_Data_1D.f90  V2.8
!
! Description - 
!
! Reference -
!
! Calling program - pressure_gradient_1d
!
! External calls - astro_params, define_1dsur_spec, define_1dsur_data,
!                  water_depths
!
! Module calls - diff_dates, error_file, lim_dims, loop_index, mult_index
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE obconds
USE physpars
USE switches
USE syspars
USE tide
USE timepars
USE error_routines, ONLY: error_alloc, error_file
USE time_routines, ONLY: diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index, mult_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL, SAVE :: hargxslope, hargyslope, harzeta
LOGICAL :: no_time_series
CHARACTER (LEN=lentime) :: ciodatetime
INTEGER :: icon, iicon, i1, i2, i3
INTEGER, SAVE :: noread, novars
REAL (KIND=kndrlong) :: ratio1, ratio2, zetax
REAL, SAVE :: gxslopedat, gyslopedat, zetadat
INTEGER (KIND=kndilong), SAVE, DIMENSION(2) :: nosecsdat
INTEGER, DIMENSION(3) :: tlims
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: surdata
REAL (KIND=kndrlong), DIMENSION(nconobc) :: phase
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:) :: phase_obc_d


procname(pglev+1) = 'update_1dsur_data'
CALL log_timer_in()

!
!1. Initialise on first call
!---------------------------
!

IF (nt.EQ.0) THEN

!
!1.1 Define o.b.c. specifiers
!----------------------------
!

   CALL define_1dsur_spec

!
!1.2 Allocate
!------------
!

   SELECT CASE (isur1dtype)
      CASE (1); novars = 3
      CASE (2); novars = 1
      CASE (3); novars = 2
   END SELECT
   ALLOCATE (surdata(novars,2),STAT=errstat)
   CALL error_alloc('surdata',2,(/novars,2/),kndrtype)
   ALLOCATE (phase_obc_d(nconobc),STAT=errstat)
   CALL error_alloc('phase_obc_d',1,(/nconobc/),kndrlong)

!
!1.3 Initialise
!--------------
!

   nosecsdat = 0
   surdata = 0.0

   IF (nconobc.GT.0) THEN
      hargxslope = ANY(gxslope_amp.NE.0.0)
      hargyslope = ANY(gyslope_amp.NE.0.0)
      harzeta = ANY(zeta_amp.NE.0.0)
      phase_obc_d = phase_obc
   ELSE
      hargxslope = .FALSE.
      hargyslope = .FALSE.
      harzeta = .FALSE.
      gxslopedat = 0.0
      gyslopedat = 0.0
      zetadat = 0
   ENDIF

ENDIF

!
!2. Update data
!--------------
!

IF (modfiles(io_1uvsur,2,1)%defined) THEN

   tlims = modfiles(io_1uvsur,2,1)%tlims
   no_time_series = lim_dims(ABS(tlims)).EQ.1
   IF ((nt.GT.tlims(1).AND.no_time_series).OR.&
     & (.NOT.loop_index(tlims,nt)).OR.&
     & modfiles(io_1uvsur,2,1)%iostat.EQ.3) GOTO 300

!
!2.1 Set parameters
!------------------
!

   SELECT CASE (isur1dtype)
      CASE (1); i1 = 1; i2 = 2; i3 = 3
      CASE (2); i1 = 0; i2 = 0; i3 = 1
      CASE (3); i1 = 1; i2 = 2; i3 = 0
   END SELECT

!
!2.2 Time-independent data or regular data without time interpolation
!--------------------------------------------------------------------
!

   IF (no_time_series.OR.&
    & (modfiles(io_1uvsur,2,1)%time_regular.AND.tlims(3).LT.0)) THEN
      IF (nt.EQ.tlims(1)) THEN
         ciodatetime = cdatetime_undef
         DO WHILE (ciodatetime.NE.CDateTime)
            CALL define_1dsur_data(ciodatetime,surdata(:,1),novars)
            IF (modfiles(io_1uvsur,2,1)%iostat.EQ.3) THEN
               CALL error_file(ierrno_fend,filepars=modfiles(io_1uvsur,2,1))
            ENDIF
         ENDDO
      ELSE
         CALL define_1dsur_data(ciodatetime,surdata(:,1),novars)
         IF (modfiles(io_1uvsur,2,1)%iostat.EQ.3) GOTO 300
         IF (ciodatetime.NE.CDateTime) GOTO 1001
      ENDIF
      
!
!2.3 Irregular data or with time interpolation
!---------------------------------------------
!

   ELSE

!
!2.3.1 Read data on first call
!-----------------------------
!

      IF (nt.EQ.tlims(1)) THEN
         CALL define_1dsur_data(ciodatetime,surdata(:,2),novars)
         IF (modfiles(io_1uvsur,2,1)%iostat.EQ.3) THEN
            CALL error_file(ierrno_fend,filepars=modfiles(io_1uvsur,2,1))
         ENDIF
         CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat(2))
         IF (nosecsdat(2).GT.nosecsrun) GOTO 1002
         noread = 1
      ENDIF

!
!2.3.2 Update data
!-----------------
!

      IF (nt.LT.nstep) THEN
         DO WHILE (nosecsdat(2).LE.nosecsrun)
!           ---store old data
            surdata(:,1) = surdata(:,2)
            nosecsdat(1) = nosecsdat(2)
!           ---read new data
            CALL define_1dsur_data(ciodatetime,surdata(:,2),novars)
            noread = noread + 1
            IF (noread.EQ.2.AND.modfiles(io_1uvsur,2,1)%iostat.EQ.3) THEN
               CALL error_file(ierrno_fend,filepars=modfiles(io_1uvsur,2,1))
            ENDIF
            CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat(2))
         ENDDO
      ENDIF

   ENDIF

!
!2.4 Evaluate data at current time
!---------------------------------
!
!  ---without time interpolation
   IF (no_time_series.OR.tlims(3).LT.0) THEN
      IF (i1.GT.0) gxslopedat = surdata(i1,1)
      IF (i2.GT.0) gyslopedat = surdata(i2,1)
      IF (i3.GT.0) zetadat = surdata(i3,1)

!  ---with time interpolation
   ELSE
      ratio2 = MERGE(0.0,(nosecsrun-nosecsdat(1))&
                  & /REAL(nosecsdat(2)-nosecsdat(1)),tlims(3).LT.0)
      ratio1 = 1.0 - ratio2
      gxslopedat = ratio1*surdata(i1,1) + ratio2*surdata(i1,2)
      gyslopedat = ratio1*surdata(i2,1) + ratio2*surdata(i2,2)
      zetadat = ratio1*surdata(i3,1) + ratio2*surdata(i3,2)
   ENDIF

ENDIF

!
!3. Update nodal corrections and phases
!--------------------------------------
!

300 CONTINUE

IF (iopt_astro_pars.GT.0.AND.(nt.EQ.0.OR.mult_index(nt,icnodal))) THEN
   CALL astro_params(index_obc,fnode_obc,phase_obc,nconobc,dlon_ref_obc,&
                   & CDateTime)
   phase_obc_d = phase_obc
ELSEIF (nt.GT.0) THEN
   icon_310: DO icon=1,nconobc
      iicon = index_obc(icon)
      phase_obc_d(icon) = MOD(phase_obc_d(icon)+&
                            & delt2d*tidal_spectrum(iicon),twopi_d)
   ENDDO icon_310
   phase_obc = phase_obc_d
ENDIF

!
!4. Surface forcing
!------------------
!
!---X-component of pressure gradient
IF (hargxslope) THEN
   phase = MOD(phase_obc_d-gxslope_pha,twopi_d)
   gxslope = gxslopedat +  SUM(fnode_obc*gxslope_amp*COS(phase))
ELSE
   gxslope = gxslopedat
ENDIF

!---Y-component of pressure gradient
IF (hargyslope) THEN
   phase = MOD(phase_obc_d-gyslope_pha,twopi_d)
   gyslope = gyslopedat +  SUM(fnode_obc*gyslope_amp*COS(phase))
ELSE
   gyslope = gyslopedat
ENDIF

!---surface elevation
IF (harzeta) THEN
   phase = MOD(phase_obc_d-zeta_pha,twopi_d)
   zetax = zetadat +  SUM(fnode_obc*zeta_amp*COS(phase))
ELSE
   zetax = zetadat
ENDIF

!
!5. Update elevations and depth arrays
!-------------------------------------
!

WHERE (maskatc_int)
   zeta(1:ncloc,1:nrloc) = zetax
   deptotatc_old(1:ncloc,1:nrloc) = deptotatc(1:ncloc,1:nrloc)
END WHERE

WHERE (node2du(1:ncloc,1:nrloc).GT.1)
   deptotatu_old(1:ncloc,:) = deptotatu(1:ncloc,1:nrloc)
END WHERE

WHERE (node2dv(1:ncloc,1:nrloc).GT.1)
   deptotatv_old(:,1:nrloc) = deptotatv(1:ncloc,1:nrloc)
END WHERE

CALL water_depths

!
!6. Deallocate arrays
!--------------------
!

IF (cold_start.OR.(iopt_tidal_accel.EQ.0.AND.nt.EQ.nstep).OR.&
 & (iopt_tidal_accel.EQ.1.AND.nt.EQ.hydro_end)) THEN
   DEALLOCATE (surdata,phase_obc_d)
ENDIF

CALL log_timer_out()


RETURN

!---exit code in case of error
1001 nerrs = 1
IF (errchk.AND.nerrs.LE.maxerrors) THEN
   WRITE (ioerr,'(A)') 'Read invalid date/time: '//ciodatetime
   WRITE (ioerr,'(A)') 'Should be: '//CDateTime
ENDIF
CALL error_file(ierrno_read,filepars=modfiles(io_1uvsur,2,1))
1002 nerrs = 1
IF (errchk.AND.nerrs.LE.maxerrors) THEN
   WRITE (ioerr,'(A)') 'Read invalid date/time: '//ciodatetime
   WRITE (ioerr,'(A)') 'Should be no later than: '//CDateTime
ENDIF
CALL error_file(ierrno_read,filepars=modfiles(io_1uvsur,2,1))

END SUBROUTINE update_1dsur_data

!========================================================================

SUBROUTINE write_1dsur_data(ciodatetime,data1d,novars)
!************************************************************************
!
! *write_1dsur_data* Write 1-D forcing data in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Boundary_Data_1D.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - define_1dsur_data
!
! Module calls - error_alloc_struc, open_filepars, output_flag,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_time, write_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE syspars
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: open_filepars, output_flag, write_atts_mod, &
                        & write_time, write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: ciodatetime
INTEGER, INTENT(IN) :: novars
REAL, INTENT(IN), DIMENSION(novars) :: data1d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR    Date/time in data file
!*data1d*      REAL    Surface data
!*novars*      INTEGER Number of surface data
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
INTEGER :: ivar, numvars
TYPE(FileParams) :: filepars
INTEGER, DIMENSION(novars) :: vecids
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_1dsur_data'
CALL log_timer_in()

!
!1. Disable writing if needed
!----------------------------
!

filepars = modfiles(io_1uvsur,2,2)
flag = output_flag(ciodatetime,filepars%tskips)
IF (.NOT.flag) GOTO 1000
numvars = novars + 1

!
!2. Write header on first call
!-----------------------------
!

IF (filepars%iostat.EQ.0) THEN

!  ---allocate
   IF (ALLOCATED(varatts)) DEALLOCATE(varatts)
   ALLOCATE (varatts(novars+1),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/novars+1/),'VariableAtts')

!  ---file attributes
   CALL set_modfiles_atts(io_1uvsur,2,2)
   filepars = modfiles(io_1uvsur,2,2)

!  ---open file
   CALL open_filepars(filepars)

!  ---variable attributes
   CALL set_modvars_atts(io_1uvsur,2,2,varatts,numvars)

!  ---write attributes
   CALL write_atts_mod(filepars,varatts,numvars)

ENDIF

!
!3. Write data
!-------------
!
!---date/time
CALL write_time(ciodatetime,filepars)

!---data
vecids = (/(ivar,ivar=2,numvars)/)
!CALL write_vars(data1d,filepars,0,varatts=varatts(2:numvars),vecids=vecids)
CALL write_vars(data1d,filepars,0,vecids=vecids)

modfiles(io_1uvsur,2,2) = filepars

1000 CALL log_timer_out()


RETURN

END SUBROUTINE write_1dsur_data

!========================================================================

SUBROUTINE write_1dsur_spec
!************************************************************************
!
! *write_1dsur_spec* Write specifier arrays for 1-D forcing in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Boundary_Data_1D.f90  V2.0
!
! Description - 
!
! Reference -
!
! Calling program - define_1dsur_spec
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE obconds
USE paralpars
USE switches
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_1dsur_spec'
CALL log_timer_in()

!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(io_1uvsur,1,2)
filepars = modfiles(io_1uvsur,1,2)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(io_1uvsur,1,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2. Write data
!-------------
!
!---type of data
varid = 1
CALL write_vars(isur1dtype,filepars,varid,varatts)

!---amplitudes and phases of X-component of pressure gradient
IF (isur1dtype.EQ.1.OR.isur1dtype.EQ.3) THEN
   varid = varid + 1
   CALL write_vars(gxslope_amp,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(gxslope_pha,filepars,varid,varatts=varatts)
ENDIF

!---amplitudes and phases of Y-component of pressure gradient
IF (isur1dtype.EQ.1.OR.isur1dtype.EQ.3) THEN
   varid = varid + 1
   CALL write_vars(gyslope_amp,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(gyslope_pha,filepars,varid,varatts=varatts)
ENDIF

!---surface elevation
IF (isur1dtype.EQ.1.OR.isur1dtype.EQ.2) THEN
   varid = varid + 1
   CALL write_vars(zeta_amp,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(zeta_pha,filepars,varid,varatts=varatts)
ENDIF

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_1uvsur,1,2) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE write_1dsur_spec
