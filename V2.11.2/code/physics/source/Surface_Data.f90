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
! *Surface_Data* Define surface grids and surface data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Data.f90  V2.11.1
!
! $Date: 2018-05-04 16:55:52 +0200 (Fri, 04 May 2018) $
!
! $Revision: 1128 $
!
! Description - 
!
! Reference -
!
! Subroutines - define_surface_data, read_surface_data, update_surface_data,
!               write_surface_data
!
!************************************************************************
!

!========================================================================

SUBROUTINE define_surface_data(iddesc,ifil,ciodatetime,surdata,n1dat,n2dat,&
                             & novars)
!************************************************************************
!
! *define_surface_data* Obtain new set of global surface data 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Data.f90  V2.1.0
!
! Description - 
!
! Reference -
!
! Calling program - checking_surface_data, update_surface_data
!
! External calls - read_surface_data, usrdef_surface_data, write_surface_data
!
! Module calls - error_file, suspend_proc
!
!************************************************************************
!
USE iopars
USE syspars
USE timepars
USE error_routines, ONLY: error_file
USE time_routines, ONLY: log_timer_in, log_timer_out, suspend_proc

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, n1dat, n2dat, novars
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat,novars) :: surdata

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*surdata*     REAL     Data array
!*n1dat*       INTEGER  Number of data points in X-direction
!*n2dat*       INTEGER  Number of data points in Y-direction
!*novars*      INTEGER  Number of surface data
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: define_data
INTEGER :: endfile, nosecs


procname(pglev+1) = 'define_surface_data'
CALL log_timer_in()

!
!1. Initialise
!-------------
!

endfile = modfiles(iddesc,ifil,1)%endfile
nosecs = 0
surdata = 0.0

!
!2. Open file on first call
!--------------------------
!

DO WHILE (modfiles(iddesc,ifil,1)%iostat.EQ.0)

!  ---standard
   IF (modfiles(iddesc,ifil,1)%status.EQ.'R') THEN
      CALL read_surface_data(iddesc,ifil,ciodatetime,surdata,n1dat,n2dat,novars)
!  ---user-defined
   ELSEIF (modfiles(iddesc,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_surface_data(iddesc,ifil,ciodatetime,surdata,n1dat,n2dat,&
                             & novars)
   ENDIF

!  ---suspend/abort if needed
   IF (modfiles(iddesc,ifil,1)%iostat.EQ.-1) THEN
      SELECT CASE (endfile)
         CASE (0:1)
            CALL error_file(ierrno_fopen,filepars=modfiles(iddesc,ifil,1))
         CASE (2)
            CALL suspend_proc(nowaitsecs)
            nosecs = nosecs + nowaitsecs
            IF (nosecs.GT.maxwaitsecs) THEN
               CALL error_file(ierrno_fopen,filepars=modfiles(iddesc,ifil,1))
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
   IF (modfiles(iddesc,ifil,1)%status.EQ.'R') THEN
      CALL read_surface_data(iddesc,ifil,ciodatetime,surdata,n1dat,n2dat,novars)
!  ---user-defined
   ELSEIF (modfiles(iddesc,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_surface_data(iddesc,ifil,ciodatetime,surdata,n1dat,n2dat,&
                                & novars)
   ENDIF
   IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) ciodatetime = cdatetime_undef

!  ---suspend/abort if needed
   IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
      SELECT CASE (endfile)
         CASE (0)
            CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
         CASE (1)
            define_data = .FALSE.
         CASE (2)
            CALL suspend_proc(nowaitsecs)
            nosecs = nosecs + nowaitsecs
            IF (nosecs.GT.maxwaitsecs) THEN
               CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
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

IF (modfiles(iddesc,ifil,2)%defined.AND.&
  & modfiles(iddesc,ifil,1)%iostat.LT.3) THEN
   CALL write_surface_data(iddesc,ifil,2,ciodatetime,surdata,n1dat,n2dat,novars)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE define_surface_data

!========================================================================

SUBROUTINE read_surface_data(iddesc,ifil,ciodatetime,surdata,n1dat,n2dat,&
                           & novars)
!************************************************************************
!
! *read_surface_data* Read surface data on data grid in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Data.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - define_surface_data
!
! External calls -
!
! Module calls - error_abort, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_time, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE switches
USE syspars
USE error_routines, ONLY: nerrs, error_abort, error_alloc_struc
USE inout_routines, ONLY: open_filepars, read_glbatts_mod, read_time, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, n1dat, n2dat, novars
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat,novars) :: surdata

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*surdata*     REAL     Data array
!*n1dat*       INTEGER  X-dimension of data array
!*n2dat*       INTEGER  Y-dimension of data array
!*novars*      INTEGER  Number of data parameters
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12), DIMENSION(2) :: cval
INTEGER :: idgrd, ivar, l, numvars
INTEGER, DIMENSION(2) :: shape
INTEGER, DIMENSION(novars) :: vecids
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_surface_data'
CALL log_timer_in()

filepars = modfiles(iddesc,ifil,1)
numvars = novars + 1

!
!1. Surface grid id
!------------------
!

SELECT CASE (iddesc)
   CASE (io_metsur); idgrd = igrd_meteo
   CASE (io_sstsur); idgrd = igrd_sst
   CASE (io_wavsur); idgrd = igrd_waves
   CASE (io_biosur); idgrd = igrd_bio
END SELECT

!
!2. File header on first call
!----------------------------
!

IF (filepars%iostat.LE.0) THEN

!  ---open file
   filepars = modfiles(iddesc,ifil,1)
   CALL open_filepars(filepars)
   IF (filepars%endfile.EQ.2.AND.filepars%iostat.EQ.-1) GOTO 1000

!  ---global attributes
   CALL read_glbatts_mod(filepars)
   IF (filepars%novars.NE.novars) THEN
      nerrs = 1
      IF (errchk.AND.nerrs.LE.maxerrors) THEN
         WRITE (cval(1),'(I12)') filepars%novars; cval(1) = ADJUSTL(cval(1))
         WRITE (ioerr,'(A)') 'Wrong number of variables in file '&
                            & //TRIM(filepars%filename)//': '//TRIM(cval(1))
         WRITE (cval(1),'(I12)') novars; cval(1) = ADJUSTL(cval(1))
         WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cval(1))
      ENDIF
      CALL error_abort('read_surface_data',ierrno_input)
   ENDIF
   filepars%timerec = MERGE(0,filepars%maxrecs+1,iopt_part_model.LT.4)
   
!  ---allocate
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')

!  ---variable atrributes
   CALL read_varatts_mod(filepars,varatts,numvars)
   ivar_210: DO ivar=2,numvars
      IF (varatts(ivar)%global_dims(1).NE.n1dat.OR.&
        & varatts(ivar)%global_dims(2).NE.n2dat) THEN
         nerrs = 1
         IF (errchk.AND.nerrs.LE.maxerrors) THEN
            WRITE (ioerr,'(A)') 'File: '//TRIM(filepars%filename)
            l_211: DO l=1,2
               WRITE (cval(l),'(I12)') varatts(ivar)%global_dims(l)
               cval(l) = ADJUSTL(cval(l))
            ENDDO l_211
            WRITE (ioerr,'(A)') 'Wrong shape for variable '&
                              & //TRIM(varatts(ivar)%f90_name)//': '&
                              & //TRIM(cval(1))//','//TRIM(cval(2))
            shape = (/n1dat,n2dat/)
            l_212: DO l=1,2
               WRITE (cval(l),'(I12)') shape(l); cval(l) = ADJUSTL(cval(l))
            ENDDO l_212
            WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cval(1))//','&
                                                  & //TRIM(cval(2))
         ENDIF
      ENDIF
   ENDDO ivar_210

   CALL error_abort('read_surface_data',ierrno_input)

!  ---deallocate   
   DEALLOCATE (varatts)

   GOTO 1000

ENDIF

!
!2. Read date/time
!-----------------
!
!---date/time
CALL read_time(ciodatetime,filepars)

!
!3. Read surface data
!--------------------
!

IF (filepars%iostat.LT.3) THEN
   vecids = (/(ivar,ivar=2,numvars)/)
   IF (surfacegrids(idgrd,ifil)%nhtype.GT.0.OR.filepars%form.EQ.'N') THEN
      CALL read_vars(surdata,filepars,0,vecids=vecids)
   ELSE
      CALL read_vars(surdata(1,1,:),filepars,0,vecids=vecids)
   ENDIF
ENDIF

1000 modfiles(iddesc,ifil,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_surface_data

!========================================================================

SUBROUTINE update_surface_data(iddesc,ifil,surindat,survals,surfgrid,n1dat,&
                             & n2dat,novars,n1grd,n2grd,nosecsdat,angle)
!************************************************************************
!
! *update_surface_data* Update surface input defined on data grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Data.f90  V2.11.1
!
! Description - 
!
! Reference -
!
! Calling program - meteo_input, temperature_equation, wave_input
!
! External calls - define_surface_data
!
! Module calls - diff_dates, distribute_mod, error_file,
!                intpol_data_to_model_2d, lim_dims, loop_index
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE syspars
USE timepars
USE error_routines, ONLY: nerrs, error_file
USE grid_interp, ONLY: intpol_data_to_model_2d
USE paral_comms, ONLY: distribute_mod
USE time_routines, ONLY: diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index

IMPLICIT NONE

!
!*Arguments
!
LOGICAL, INTENT(IN), DIMENSION(novars) :: angle
INTEGER, INTENT(IN) :: iddesc, ifil, novars, n1dat, n1grd, n2dat, n2grd
INTEGER (KIND=kndilong), INTENT(INOUT), DIMENSION(2) :: nosecsdat
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat,novars,2) :: surindat
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,novars) :: survals
TYPE (HRelativeCoords), INTENT(IN), DIMENSION(n1grd,n2grd) :: surfgrid

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*     INTEGER  Data file id
!*ifil*       INTEGER  No. of data file
!*surindat*   REAL    Surface data at old and new data time
!*survals*    REAL    Surface data interpolated in time and on model grid
!*surfgrid*   DERIVED Model grid coordinates with respect to surface data grid
!*n1dat*      INTEGER X-dimension of input data
!*n2dat*      INTEGER Y-dimension of input data
!*novars*     INTEGER Number of data parameters
!*n1grd*      INTEGER X-dimension of data grid
!*n2grd*      INTEGER Y-dimension of data grid
!*nosecsdat*  LONGINT No. of seconds between start of simulation and old/new
!                     data time
!*angle*      LOGICAL .TRUE. if a data variable represnets an angle,
!                     .FALSE. otherwise 
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: no_time_series
CHARACTER (LEN=lentime) :: ciodatetime
INTEGER :: datflag, idgrd, ivar
REAL (KIND=kndrlong) :: ratio1, ratio2
INTEGER, DIMENSION(3) :: tlims
INTEGER, DIMENSION(4) :: nhdist
INTEGER, SAVE, DIMENSION(MaxIOTypes,MaxIOFiles) :: noread
REAL, DIMENSION(n1dat,n2dat,novars) :: surdat

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR    Input date/time
!*idgrd*       INTEGER Surface grid id
!
!------------------------------------------------------------------------------
!
!1. Initialise
!-------------
!
!---return if no update
tlims = modfiles(iddesc,ifil,1)%tlims
no_time_series = lim_dims(ABS(tlims)).EQ.1
IF ((nt.GT.tlims(1).AND.no_time_series).OR.(.NOT.loop_index(tlims,nt)).OR.&
   & modfiles(iddesc,ifil,1)%iostat.EQ.3) RETURN

procname(pglev+1) = 'update_surface_data'
CALL log_timer_in()

!---surface grid id
SELECT CASE (iddesc)
   CASE (io_metsur); idgrd = igrd_meteo
   CASE (io_sstsur); idgrd = igrd_sst
   CASE (io_wavsur); idgrd = igrd_waves
   CASE (io_biosur); idgrd = igrd_bio
END SELECT

!---data flagging   
datflag = surfacegrids(idgrd,1)%datflag

!
!2. Read/update data in time
!---------------------------
!
!2.1 Time-independent data or regular data without time interpolation
!--------------------------------------------------------------------
!

IF (no_time_series.OR.&
 & (modfiles(iddesc,ifil,1)%time_regular.AND.tlims(3).LT.0.AND.&
  & nt.LT.nstep)) THEN
   IF (nt.EQ.tlims(1)) THEN
      ciodatetime = cdatetime_undef
      DO WHILE (ciodatetime.NE.CDateTime)
         CALL define_surface_data(iddesc,ifil,ciodatetime,surindat(:,:,:,1),&
                                & n1dat,n2dat,novars)
         IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
            CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
         ENDIF
      ENDDO
   ELSE
      CALL define_surface_data(iddesc,ifil,ciodatetime,surindat(:,:,:,1),&
                             & n1dat,n2dat,novars)
      IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) GOTO 1000
      IF (ciodatetime.NE.CDateTime) GOTO 1001
   ENDIF
      
!
!2.2 Irregular data or with time interpolation
!---------------------------------------------
!

ELSE

!
!2.2.1 Read data on first call
!-----------------------------
!

   IF (nt.EQ.tlims(1)) THEN
      CALL define_surface_data(iddesc,ifil,ciodatetime,surindat(:,:,:,2),&
                             & n1dat,n2dat,novars)
      IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
         CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
      ENDIF
      CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat(2))
      IF (nosecsdat(2).GT.nosecsrun) GOTO 1002
      noread(iddesc,ifil) = 1
   ENDIF

!
!2.2.2 Update data
!-----------------
!

   IF (nt.LT.nstep) THEN
      DO WHILE (nosecsdat(2).LE.nosecsrun)
!        ---store old data
         surindat(:,:,:,1) = surindat(:,:,:,2)      
         nosecsdat(1) = nosecsdat(2)
!        ---read new data
         CALL define_surface_data(iddesc,ifil,ciodatetime,surindat(:,:,:,2),&
                                & n1dat,n2dat,novars)
         noread(iddesc,ifil) = noread(iddesc,ifil) + 1
         IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
            IF (modfiles(iddesc,ifil,1)%endfile.EQ.1) THEN
               GOTO 1000
            ELSEIF (noread(iddesc,ifil).EQ.2) THEN
               CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
            ENDIF
         ENDIF
         CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat(2))
      ENDDO
   ENDIF

ENDIF

!
!3. Evaluate data at current time
!--------------------------------
!
!---without time interpolation
IF (no_time_series.OR.tlims(3).LT.0) THEN
   surdat = MERGE(surindat(:,:,:,1),float_fill,&
                & ABS(surindat(:,:,:,1)-float_fill).GT.float_fill_eps)

!---with time interpolation
ELSE
   ratio2 = (nosecsrun-nosecsdat(1))/REAL(nosecsdat(2)-nosecsdat(1))
   ratio1 = 1.0-ratio2
   ivar_310: DO ivar=1,novars
      IF (.NOT.angle(ivar)) THEN
         WHERE ((ABS(surindat(:,:,ivar,1)-float_fill).GT.float_fill_eps).AND.&
              & (ABS(surindat(:,:,ivar,2)-float_fill).GT.float_fill_eps))
            surdat(:,:,ivar) = ratio1*surindat(:,:,ivar,1) + &
                             & ratio2*surindat(:,:,ivar,2)
         ELSEWHERE
            surdat(:,:,ivar) = float_fill
         END WHERE
      ELSE
         WHERE ((ABS(surindat(:,:,ivar,1)-float_fill).GT.float_fill_eps).AND.&
              & (ABS(surindat(:,:,ivar,2)-float_fill).GT.float_fill_eps))
            surdat(:,:,ivar) = ratio1*surindat(:,:,ivar,1) + &
                             & ratio2*surindat(:,:,ivar,2) + &
                             & pi*(1.0-SIGN(1.0,pi-&
                             & ABS(surindat(:,:,ivar,2)-surindat(:,:,ivar,1)))) 
         ELSEWHERE
            surdat(:,:,ivar) = float_fill
         END WHERE
      ENDIF
   ENDDO ivar_310
ENDIF

!
!4. Interpolate in space
!-----------------------
!

SELECT CASE (surfacegrids(idgrd,ifil)%nhtype)
   CASE (0)
      ivar_410: DO ivar=1,novars
         survals(:,:,ivar) = surdat(1,1,ivar)
      ENDDO ivar_410
   CASE (1:3)
      CALL intpol_data_to_model_2d(surdat,survals,n1dat,n2dat,novars,surfgrid,&
                                 & datflag,.FALSE.,inflag=float_fill)
   CASE (4)
      nhdist  = 0
      CALL distribute_mod(surdat,survals,(/1,1,1/),nhdist,0,0.0,0)
END SELECT

1000 CALL log_timer_out()


RETURN

!---exit code in case of error
1001 nerrs = 1
IF (errchk.AND.nerrs.LE.maxerrors) THEN
   WRITE (ioerr,'(A)') 'Read invalid date/time: '//ciodatetime
   WRITE (ioerr,'(A)') 'Should be: '//CDateTime
ENDIF
CALL error_file(ierrno_read,filepars=modfiles(iddesc,ifil,1))
1002 nerrs = 1
IF (errchk.AND.nerrs.LE.maxerrors) THEN
   WRITE (ioerr,'(A)') 'Read invalid date/time: '//ciodatetime
   WRITE (ioerr,'(A)') 'Should be no later than: '//CDateTime
ENDIF
CALL error_file(ierrno_read,filepars=modfiles(iddesc,ifil,1))

END SUBROUTINE update_surface_data

!========================================================================

SUBROUTINE write_surface_data(iddesc,ifil,iotype,ciodatetime,surdata,n1dat,&
                            & n2dat,novars)
!************************************************************************
!
! *write_surface_data* Write surface data in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Surface_Data.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - define_surface_data
!
! External calls -
!
! Module calls - error_alloc_struc, later, noearlier, open_filepars,
!                output_flag, set_modfiles_atts, set_modvars_atts,
!                write_atts_mod, write_time, write_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE syspars
USE error_routines, ONLY: error_abort, error_alloc_struc
USE inout_routines, ONLY: open_filepars, output_flag, write_atts_mod, &
                        & write_time, write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, iotype, n1dat, n2dat, novars
REAL, INTENT(IN), DIMENSION(n1dat,n2dat,novars) :: surdata

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*iotype*      INTEGER  Type of surface grid
!                 = 1 => input
!                 = 2 => output
!*ciodatetime* CHAR     Date/time in data file
!*surdata*     REAL     Data array
!*n1dat*       INTEGER  X-dimension of data array
!*n2dat*       INTEGER  Y-dimension of data array
!*novars*      INTEGER  Number of data parameters
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
INTEGER :: ivar, numvars
TYPE (FileParams) :: filepars
INTEGER, DIMENSION(novars) :: vecids
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_surface_data'
CALL log_timer_in()

!
!1. Disable writing if needed
!----------------------------
!

filepars = modfiles(iddesc,ifil,2)
flag = output_flag(ciodatetime,filepars%tskips)
IF (.NOT.flag) GOTO 1000
numvars = novars + 1

!
!3. Write file header on first call
!----------------------------------
!

IF (filepars%iostat.EQ.0) THEN

!  ---file attributes
   CALL set_modfiles_atts(iddesc,ifil,2)
   filepars = modfiles(iddesc,ifil,2)

!  ---open file
   CALL open_filepars(filepars)

!  ---allocate
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')

!  ---variable attributes
   CALL set_modvars_atts(iddesc,ifil,iotype,varatts,numvars)

!  ---write attributes
   CALL write_atts_mod(filepars,varatts,numvars)

!  ---deallocate
   DEALLOCATE (varatts)
   
ENDIF

!
!4. Write data
!-------------
!
!---date/time
CALL write_time(ciodatetime,filepars)

!---surface data
vecids = (/(ivar,ivar=2,numvars)/)
IF (n1dat*n2dat.GT.1.OR.filepars%form.EQ.'N') THEN
   CALL write_vars(surdata,filepars,0,vecids=vecids)
ELSE
   CALL write_vars(surdata(1,1,:),filepars,0,vecids=vecids)
ENDIF

modfiles(iddesc,ifil,2) = filepars

1000 CALL log_timer_out()


RETURN

END SUBROUTINE write_surface_data
