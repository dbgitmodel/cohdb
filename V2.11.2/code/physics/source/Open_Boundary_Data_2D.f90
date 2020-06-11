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
! *Open_Boundary_Data_2D* Ensemble of routines for 2-D normal and tangential
!                         open boundary conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_2D.f90  V2.11
!
! $Date: 2017-11-07 14:38:09 +0100 (Tue, 07 Nov 2017) $
!
! $Revision: 1060 $
!
! Description - 
!
! Reference -
!
! Subroutines - define_2dobc_data, define_2dobc_spec, read_2dobc_data,
!               read_2dobc_spec, update_2dobc_data, write_2dobc_data,
!               write_2dobc_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE define_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
!************************************************************************
!
! *define_2dobc_data* Obtain open boundary data for 2-D normal or
!                     tangential open boundary conditions 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_2D.f90  V2.7
!
! Description - 
!
! Reference -
!
! Calling program - update_2dobc_data
!
! External calls - read_2dobc_data, usrdef_2dobc_data, write_2dobc_data
!
! Module calls - error_file, suspend_proc
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE obconds
USE syspars
USE timepars
USE error_routines, ONLY: error_file
USE time_routines, ONLY: log_timer_in, log_timer_out, suspend_proc

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
REAL, INTENT(INOUT), DIMENSION(nodat,novars) :: data2d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER   Data file id (io_2uvobc or io_2xyobc)
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*data2d*      REAL    Input data
!*nodat*       INTEGER Number of input data
!*novars*      INTEGER Number of input parameters
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: define_data
INTEGER :: endfile, l, ll, nosecs


procname(pglev+1) = 'define_2dobc_data'
CALL log_timer_in()

!
!1. Initialise
!-------------
!

endfile = modfiles(iddesc,ifil,1)%endfile
nosecs = 0

!
!2. Open file on first call
!--------------------------
!

DO WHILE (modfiles(iddesc,ifil,1)%iostat.EQ.0)

!  ---standard
   IF (modfiles(iddesc,ifil,1)%status.EQ.'R') THEN
      CALL read_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
!  ---user-defined
   ELSEIF (modfiles(iddesc,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
   ENDIF

!   ---suspend/abort if needed
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


!3. Obtain data
!--------------
!

define_data = .TRUE.

DO WHILE (define_data)

!  ---read
   IF (modfiles(iddesc,ifil,1)%status.EQ.'R') THEN
      CALL read_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
!  ---user-defined
   ELSEIF (modfiles(iddesc,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
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
   CALL write_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
ENDIF

!
!5. Outflow conditions at river boundaries
!-----------------------------------------
!

IF (iddesc.EQ.io_2uvobc.AND.modfiles(io_2uvobc,ifil,1)%iostat.LT.3.AND.&
  & iobc2dtype(ifil).NE.2) THEN
   l_510: DO l=1,nodat
      ll = index2dobuv(l,ifil)
      IF (ll.LE.nobu) THEN
         IF (ll.GT.nosbu.AND.(.NOT.westobu(ll))) data2d(l,1) = -data2d(l,1)
      ELSE
         ll = ll - nobu
         IF (ll.GT.nosbv.AND.(.NOT.soutobv(ll))) data2d(l,1) = -data2d(l,1)
      ENDIF
   ENDDO l_510
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE define_2dobc_data

!========================================================================

SUBROUTINE define_2dobc_spec(iddesc)
!************************************************************************
!
! *define_2dobc_spec* Define specifier arrays for 2-D normal or tangential
!                     open boundary conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_2D.f90  V2.7
!
! Description - 
!
! Reference -
!
! Calling program - current_2d
!
! External calls - read_2dobc_spec, usrdef_2dobc_spec, write_2dobc_spec
!
! Module calls - error_abort, error_alloc, error_limits_arr
!
!************************************************************************
!
USE gridpars
USE iopars
USE obconds
USE syspars
USE tide
USE error_routines, ONLY: error_abort, error_alloc, error_limits_arr
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name      Type      Purpose
!-------------------------------------------------------------------------
!*iddesc*   INTEGER   Data file id (io_2uvobc or io_2xyobc)
!
!-------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: idat, ifil, ii, jj, nofiles


procname(pglev+1) = 'define_2dobc_spec'
CALL log_timer_in()


nofiles = maxdatafiles(iddesc,1)

!
!1. Allocate specifier arrays
!----------------------------
!

SELECT CASE (iddesc)

!
!1.1 (U,V) open boundaries
!-------------------------
!

   CASE (io_2uvobc)

!     ---type of conditions at open boundaries
      ALLOCATE (ityp2dobu(nobu),STAT=errstat)
      CALL error_alloc('ityp2dobu',1,(/nobu/),kndint)
      IF (nobu.GT.0) ityp2dobu = 0
      ALLOCATE (ityp2dobv(nobv),STAT=errstat)
      CALL error_alloc('ityp2dobv',1,(/nobv/),kndint)
      IF (nobv.GT.0) ityp2dobv = 0
      ALLOCATE (iloczobu(nobu),STAT=errstat)
      CALL error_alloc('iloczobu',1,(/nobu/),kndint)
      IF (nobu.GT.0) iloczobu = 1
      ALLOCATE (iloczobv(nobv),STAT=errstat)
      CALL error_alloc('iloczobv',1,(/nobv/),kndint)
      IF (nobv.GT.0) iloczobv = 1
      ALLOCATE (ityp2dobx(nobx),STAT=errstat)
      CALL error_alloc('ityp2dobx',1,(/nobx/),kndint)
      IF (nobx.GT.0) ityp2dobx = 0
      ALLOCATE (ityp2doby(noby),STAT=errstat)
      CALL error_alloc('ityp2doby',1,(/noby/),kndint)
      IF (noby.GT.0) ityp2doby = 0

!     ---discharge along open boundary sections
      ALLOCATE (iqsecobu(nqsecobu,2),STAT=errstat)
      CALL error_alloc('iqsecobu',2,(/nqsecobu,2/),kndrtype)
      IF (nqsecobu.GT.0) iqsecobu = 0.0
      ALLOCATE (jqsecobv(nqsecobv,2),STAT=errstat)
      CALL error_alloc('jqsecobv',2,(/nqsecobv,2/),kndrtype)
      IF (nqsecobv.GT.0) jqsecobv = 0.0

!     ---number of data and index maps per input file
      IF (nofiles.GT.1) THEN
         ALLOCATE (no2dobuv(2:nofiles),STAT=errstat)
         CALL error_alloc('no2dobuv',1,(/nofiles-1/),kndint)
         no2dobuv = 0
         ALLOCATE (index2dobuv(nobu+nobv,2:nofiles),STAT=errstat)
         CALL error_alloc('index2dobuv',2,(/nobu+nobv,nofiles-1/),kndint)
         IF (SIZE(index2dobuv).GT.0) index2dobuv = 0
         ALLOCATE (iobc2dtype(2:nofiles),STAT=errstat)
         CALL error_alloc('iobc2dtype',1,(/nofiles-1/),kndint)
         iobc2dtype = 0
      ENDIF

!     ---amplitudes
      ALLOCATE (udatobu_amp(nobu,nconobc),STAT=errstat)
      CALL error_alloc('udatobu_amp',2,(/nobu,nconobc/),kndrtype)
      IF (SIZE(udatobu_amp).GT.0) udatobu_amp = 0.0
      ALLOCATE (zdatobu_amp(nobu,nconobc),STAT=errstat)
      CALL error_alloc('zdatobu_amp',2,(/nobu,nconobc/),kndrtype)
      IF (SIZE(zdatobu_amp).GT.0) zdatobu_amp = 0.0
      ALLOCATE (vdatobv_amp(nobv,nconobc),STAT=errstat)
      CALL error_alloc('vdatobv_amp',2,(/nobv,nconobc/),kndrtype)
      IF (SIZE(vdatobv_amp).GT.0) vdatobv_amp = 0.0
      ALLOCATE (zdatobv_amp(nobv,nconobc),STAT=errstat)
      CALL error_alloc('zdatobv_amp',2,(/nobv,nconobc/),kndrtype)
      IF (SIZE(zdatobv_amp).GT.0) zdatobv_amp = 0.0

!     ---phases
      ALLOCATE (udatobu_pha(nobu,nconobc),STAT=errstat)
      CALL error_alloc('udatobu_pha',2,(/nobu,nconobc/),kndrtype)
      IF (SIZE(udatobu_pha).GT.0) udatobu_pha = 0.0
      ALLOCATE (zdatobu_pha(nobu,nconobc),STAT=errstat)
      CALL error_alloc('zdatobu_pha',2,(/nobu,nconobc/),kndrtype)
      IF (SIZE(zdatobu_pha).GT.0) zdatobu_pha = 0.0
      ALLOCATE (vdatobv_pha(nobv,nconobc),STAT=errstat)
      CALL error_alloc('vdatobv_pha',2,(/nobv,nconobc/),kndrtype)
      IF (SIZE(vdatobv_pha).GT.0) vdatobv_pha = 0.0
      ALLOCATE (zdatobv_pha(nobv,nconobc),STAT=errstat)
      CALL error_alloc('zdatobv_pha',2,(/nobv,nconobc/),kndrtype)
      IF (SIZE(zdatobv_pha).GT.0) zdatobv_pha = 0.0

!     ---transport and elevation at open boundaries
      ALLOCATE (udatobu(nobu,2),STAT=errstat)
      CALL error_alloc('udatobu',2,(/nobu,2/),kndrtype)
      IF (nobu.GT.0) udatobu = 0.0
      ALLOCATE (zdatobu(nobu,2),STAT=errstat)
      CALL error_alloc('zdatobu',2,(/nobu,2/),kndrtype)
      IF (nobu.GT.0) zdatobu = 0.0
      ALLOCATE (zdatobu_old(nobu),STAT=errstat)
      CALL error_alloc('zdatobu_old',1,(/nobu/),kndrtype)
      IF (nobu.GT.0) zdatobu_old = 0.0
      ALLOCATE (vdatobv(nobv,2),STAT=errstat)
      CALL error_alloc('vdatobv',2,(/nobv,2/),kndrtype)
      IF (nobv.GT.0) vdatobv = 0.0
      ALLOCATE (zdatobv(nobv,2),STAT=errstat)
      CALL error_alloc('zdatobv',2,(/nobv,2/),kndrtype)
      IF (nobv.GT.0) zdatobv = 0.0
      ALLOCATE (zdatobv_old(nobv),STAT=errstat)
      CALL error_alloc('zdatobv_old',1,(/nobv/),kndrtype)
      IF (nobu.GT.0) zdatobv_old = 0.0

!
!1.2 Tangential open boundary data
!------- -------------------------
!

   CASE (io_2xyobc)

!     ---number of data and index maps per input file
      IF (nofiles.GT.1) THEN
         ALLOCATE (no2dobxy(2:nofiles),STAT=errstat)
         CALL error_alloc('no2dobxy',1,(/nofiles-1/),kndint)
         no2dobxy = 0
         ALLOCATE (index2dobxy(nobx+noby,2:nofiles),STAT=errstat)
         CALL error_alloc('index2dobxy',2,(/nobx+noby,nofiles-1/),kndint)
         IF (SIZE(index2dobxy).GT.0) index2dobxy = 0
      ENDIF

!     ---amplitudes
      ALLOCATE (vdatobx_amp(nobx,nconobc),STAT=errstat)
      CALL error_alloc('vdatobx_amp',2,(/nobx,nconobc/),kndrtype)
      IF (SIZE(vdatobx_amp).GT.0) vdatobx_amp = 0.0
      ALLOCATE (udatoby_amp(noby,nconobc),STAT=errstat)
      CALL error_alloc('udatoby_amp',2,(/noby,nconobc/),kndrtype)      
      IF (SIZE(udatoby_amp).GT.0) udatoby_amp = 0.0
 
!     ---phases
      ALLOCATE (vdatobx_pha(nobx,nconobc),STAT=errstat)
      CALL error_alloc('vdatobx_pha',2,(/nobx,nconobc/),kndrtype)
      IF (SIZE(vdatobx_pha).GT.0) vdatobx_pha = 0.0
      ALLOCATE (udatoby_pha(noby,nconobc),STAT=errstat)
      CALL error_alloc('udatoby_pha',2,(/noby,nconobc/),kndrtype)
      IF (SIZE(udatoby_pha).GT.0) udatoby_pha = 0.0

!     ---transport and elevation at open boundaries
      ALLOCATE (vdatobx(nobx,2),STAT=errstat)
      CALL error_alloc('vdatobx',2,(/nobx,2/),kndrtype)
      IF (nobx.GT.0) vdatobx = 0.0
      ALLOCATE (udatoby(noby,2),STAT=errstat)
      CALL error_alloc('udatoby',2,(/noby,2/),kndrtype)
      IF (noby.GT.0) udatoby = 0.0

END SELECT

!
!2. Obtain specifier arrays
!--------------------------
!
!---read
IF (modfiles(iddesc,1,1)%status.EQ.'R') THEN
   CALL read_2dobc_spec(iddesc)
!---user-defined
ELSEIF (modfiles(iddesc,1,1)%status.EQ.'N') THEN
   CALL usrdef_2dobc_spec(iddesc)
ENDIF

!
!3. Default values
!-----------------
!

IF (nofiles.EQ.2)THEN

   SELECT CASE (iddesc)

!     ---(U,V)-open boundaries
      CASE (io_2uvobc)
         IF (ALL(index2dobuv.EQ.0)) THEN
            index2dobuv(1:no2dobuv(2),2) = (/(idat,idat=1,no2dobuv(2))/)
         ENDIF

!     ---(X,Y)-open boundaries
      CASE (io_2xyobc)
         IF (ALL(index2dobxy.EQ.0)) THEN
            index2dobxy(1:no2dobxy(2),2) = (/(idat,idat=1,no2dobxy(2))/)
         ENDIF

   END SELECT

ENDIF

!
!4. Write specifier arrays
!-------------------------
!

IF (modfiles(iddesc,1,2)%defined) CALL write_2dobc_spec(iddesc)

!
!5. Check specifier arrays
!-------------------------
!
!5.1 Type of open boundary conditions
!------------------------------------
!

SELECT CASE (iddesc)

   CASE (io_2uvobc)
!   ---U-nodes
      ii_511: DO ii=1,nobu
         CALL error_limits_arr(ityp2dobu(ii),'ityp2dobu',0,17,1,indx=(/ii/))
         CALL error_limits_arr(iloczobu(ii),'iloczobu',0,2,1,indx=(/ii/))
      ENDDO ii_511
!    ---V-nodes
      jj_512: DO jj=1,nobv
         CALL error_limits_arr(ityp2dobv(jj),'ityp2dobv',0,17,1,indx=(/jj/))
         CALL error_limits_arr(iloczobv(jj),'iloczobv',0,2,1,indx=(/jj/))
      ENDDO jj_512

   CASE (io_2xyobc)
!     ---X-nodes
      ii_513: DO ii=1,nobx
         CALL error_limits_arr(ityp2dobx(ii),'ityp2dobx',0,2,1,indx=(/ii/))
      ENDDO ii_513
!     ---Y-nodes
      jj_514: DO jj=1,noby
         CALL error_limits_arr(ityp2doby(jj),'ityp2doby',0,2,1,indx=(/jj/))
      ENDDO jj_514

END SELECT

!
!5.2 Data file specifiers
!------------------------
!

SELECT CASE (iddesc)

!
!5.2.1 (U,V)-open boundaries
!---------------------------
!

   CASE (io_2uvobc)
      ifil_521: DO ifil=2,nofiles
         IF (modfiles(io_2uvobc,ifil,1)%defined) THEN
            CALL error_limits_arr(no2dobuv(ifil),'no2dobuv',1,nobu+nobv,&
                                & 1,(/ifil/))
            CALL error_limits_arr(iobc2dtype(ifil),'iobc2dtype',1,3,1,(/ifil/))
            idat_5211: DO idat=1,no2dobuv(ifil)
               CALL error_limits_arr(index2dobuv(idat,ifil),'index2dobuv',1,&
                                   & nobu+nobv,2,(/idat,ifil/))
            ENDDO idat_5211
         ENDIF
      ENDDO ifil_521

!
!5.2.2 (X,Y)-open boundaries
!---------------------------
!

   CASE (io_2xyobc)

      ifil_522: DO ifil=2,nofiles
         IF (modfiles(io_2xyobc,ifil,1)%defined) THEN
            CALL error_limits_arr(no2dobxy(ifil),'no2dobcxy',1,nobx+noby,&
                                & 1,(/ifil/))
            idat_5221: DO idat=1,no2dobxy(ifil)
               CALL error_limits_arr(index2dobxy(idat,ifil),'index2dobxy',1,&
                                   & nobx+noby,2,(/idat,ifil/))
            ENDDO idat_5221
         ENDIF
      ENDDO ifil_522

END SELECT

CALL error_abort('define_2dobc_spec',ierrno_input)

CALL log_timer_out()
   

RETURN

END SUBROUTINE define_2dobc_spec

!========================================================================

SUBROUTINE read_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
!************************************************************************
!
! *read_2dobc_data* Read data for 2-D normal or tangential open boundary
!                   conditions in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_2D.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - define_2dobc_data
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
CHARACTER (LEN=lentime), INTENT(OUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
REAL, INTENT(INOUT), DIMENSION(nodat,novars) :: data2d

!
! Name         Type Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER Data file id (io_2uvobc or io_2xyobc)
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*data2d*      REAL    Input data
!*nodat*       INTEGER Number of input data
!*novars*      INTEGER Number of input parameters
!*iddesc*      INTEGER Data file id
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


procname(pglev+1) = 'read_2dobc_data'
CALL log_timer_in()

filepars = modfiles(iddesc,ifil,1)
numvars = novars + 1

!
!1. File header on first call
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
         WRITE (cval,'(I12)') filepars%novars; cval = ADJUSTL(cval)
         WRITE (ioerr,'(A)') 'Wrong number of variables in file '&
                            & //TRIM(filepars%filename)//': '//TRIM(cval)
         WRITE (cval,'(I12)') novars; cval = ADJUSTL(cval)
         WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cval)
      ENDIF
      CALL error_abort('read_2dobc_data',ierrno_input)
   ENDIF

!  ---allocate
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')

!  ---variable atrributes
   CALL read_varatts_mod(filepars,varatts,numvars)
   ivar_110: DO ivar=2,novars+1
      IF (varatts(ivar)%global_dims(1).NE.nodat) THEN
         nerrs = 1
         IF (master.AND.errchk.AND.nerrs.LE.maxerrors) THEN
            WRITE (ioerr,'(A)') 'File: '//TRIM(filepars%filename)
            WRITE (cval,'(I12)') varatts(ivar)%global_dims(1)
            cval = ADJUSTL(cval)
            WRITE (ioerr,'(A)') 'Wrong data size for variable '&
                              & //TRIM(varatts(ivar)%f90_name)//': '&
                              & //TRIM(cval)
            WRITE (cval,'(I12)') nodat; cval = ADJUSTL(cval)
            WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cval)
         ENDIF
         CALL error_abort('read_2dobc_data',ierrno_input)
      ENDIF
   ENDDO ivar_110

!  ---deallocate
   DEALLOCATE (varatts)
   
   GOTO 1000   

ENDIF

!
!2. Read data
!------------
!
!---date/time
CALL read_time(ciodatetime,filepars)

!---data
IF (filepars%iostat.LT.3) THEN
   vecids = (/(ivar,ivar=2,numvars)/)
   CALL read_vars(data2d,filepars,0,vecids=vecids)
ENDIF

1000 modfiles(iddesc,ifil,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_2dobc_data

!========================================================================

SUBROUTINE read_2dobc_spec(iddesc)
!************************************************************************
!
! *read_2dobc_spec* Read specifier arrays for 2-D normal or tangential open
!                   boundary conditions in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_2D.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - define_2dobc_spec
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE obconds
USE tide
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name      Type      Purpose
!-------------------------------------------------------------------------
!*iddesc*   INTEGER   Data file id (io_2uvobc or io_2xyobc)
!
!-------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_2dobc_spec'
CALL log_timer_in()

filepars = modfiles(iddesc,1,1)

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
!2. Read specifier arrays
!------------------------
!

SELECT CASE (iddesc)

!
!2.1 Normal (U,V)-open boundaries
!--------------------------------
!

   CASE (io_2uvobc)

!
!2.1.1 Type of conditions
!------------------------
!
      
      varid = 0
!     ---U-nodes
      IF (nobu.GT.0) THEN
         varid = varid + 1
         CALL read_vars(ityp2dobu,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL read_vars(iloczobu,filepars,varid,varatts=varatts)
      ENDIF
!     ---V-nodes
      IF (nobv.GT.0) THEN
         varid = varid + 1
         CALL read_vars(ityp2dobv,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL read_vars(iloczobv,filepars,varid,varatts=varatts)
      ENDIF

!
!2.1.2 Discharge along open boundary sections
!--------------------------------------------
!

      IF (nqsecobu.GT.0) THEN
         varid = varid + 1
         CALL read_vars(iqsecobu,filepars,varid,varatts=varatts)
      ENDIF
      IF (nqsecobv.GT.0) THEN
         varid = varid + 1
         CALL read_vars(jqsecobv,filepars,varid,varatts=varatts)
      ENDIF

!
!2.1.3 Number of data and index maps per input file
!--------------------------------------------------
!

      IF (maxdatafiles(iddesc,1).GT.1) THEN
         varid = varid + 1
         CALL read_vars(no2dobuv,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL read_vars(iobc2dtype,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL read_vars(index2dobuv,filepars,varid,varatts=varatts)
      ENDIF
      
!
!2.1.4 Amplitudes and phases
!---------------------------
!
!     ---U-nodes
      IF (nobu.GT.0.AND.nconobc.GT.0) THEN
         varid = varid + 1
         CALL read_vars(udatobu_amp,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL read_vars(zdatobu_amp,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL read_vars(udatobu_pha,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL read_vars(zdatobu_pha,filepars,varid,varatts=varatts)
      ENDIF

!     ---V-nodes
      IF (nobv.GT.0.AND.nconobc.GT.0) THEN
         varid = varid + 1
         CALL read_vars(vdatobv_amp,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL read_vars(zdatobv_amp,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL read_vars(vdatobv_pha,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL read_vars(zdatobv_pha,filepars,varid,varatts=varatts)
      ENDIF

!
!2.2 Tangential (X,Y)-open boundaries
!------------------------------------
!

   CASE (io_2xyobc)

!
!2.2.1 Type of conditions
!------------------------
!

     varid = 0
!     ---X-nodes
      IF (nobx.GT.0) THEN
         varid = varid + 1
         CALL read_vars(ityp2dobx,filepars,varid,varatts=varatts)
      ENDIF
!     ---Y-nodes
      IF (noby.GT.0) THEN
         varid = varid + 1
         CALL read_vars(ityp2doby,filepars,varid,varatts=varatts)
     ENDIF

!
!2.2.2 Number of data and index maps per input file
!--------------------------------------------------
!


     IF (maxdatafiles(iddesc,1).GT.1) THEN
        varid = varid + 1
        CALL read_vars(no2dobxy,filepars,varid,varatts=varatts)
        varid = varid + 1
        CALL read_vars(index2dobxy,filepars,varid,varatts=varatts)
     ENDIF

!
!2.2.3 Amplitudes and phases
!---------------------------
!
!    ---X-nodes
     IF (nobx.GT.0.AND.nconobc.GT.0) THEN
        varid = varid + 1
        CALL read_vars(vdatobx_amp,filepars,varid,varatts=varatts)
        varid = varid + 1
        CALL read_vars(vdatobx_pha,filepars,varid,varatts=varatts)
     ENDIF

!    ---Y-nodes
     IF (noby.GT.0.AND.nconobc.GT.0) THEN
        varid = varid + 1
        CALL read_vars(udatoby_amp,filepars,varid,varatts=varatts)
        varid = varid + 1
        CALL read_vars(udatoby_pha,filepars,varid,varatts=varatts)
     ENDIF

END SELECT

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(iddesc,1,1) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE read_2dobc_spec

!========================================================================

SUBROUTINE update_2dobc_data(datavals2d,no2dobc,maxdat,nofiles,nosecsdat,iddesc)
!************************************************************************
!
! *update_2dobc_data* Update data for 2-D normal or tangential open boundary
!                     conditions
!
! Author - Patrick Luyten
!
! Description - 
!
! Version - @(COHERENS)Open_Boundary_Data_2D.f90  V2.7
!
! Reference -
!
! Calling program - current_2d, initialise_model
!
! External calls - astro_params, define_2dobc_data
!
! Module calls - combine_stats_glb, diff_dates, error_abort, error_alloc,
!                error_file, lim_dims, loop_index, mult_index
!
!************************************************************************
!
USE depths
USE fluxes
USE grid
USE gridpars
USE iopars
USE meteo
USE obconds
USE physpars
USE switches
USE syspars
USE tide
USE timepars
USE error_routines, ONLY: nerrs, error_abort, error_alloc, error_file
USE paral_comms, ONLY: combine_stats_glb
USE time_routines, ONLY: diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index, mult_index

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, maxdat, nofiles
INTEGER, INTENT(IN) , DIMENSION(2:nofiles) :: no2dobc
INTEGER (KIND=kndilong), INTENT(INOUT), DIMENSION(2:nofiles,2) :: nosecsdat
REAL, INTENT(INOUT), DIMENSION(maxdat,2,2:nofiles,2) :: datavals2d

!
! Name        Type      Purpose
!-------------------------------------------------------------------------
!*datavals2d* REAL      Data from input file
!          1st dimension => data index number
!          2nd dimension => number of variables (1 or 2)
!          3rd dimension => file number
!          4th dimension => time level (old,new)
!*no2dobc*    INTEGER   Number of o.b. data values per variable and file
!*maxdat*     INTEGER   Maximum number of data values in any data file and each
!                       variable
!*nofiles*    INTEGER   Number of data files (plus 1)
!*nosecsdat*  LONGINT   Number of seconds since start of simulation and old/new
!                       data time
!*iddesc*     INTEGER   Data file id (io_2uvobc or io_2xyobc)
!
!-------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: no_time_series
CHARACTER (LEN=12) :: cpos, ctyp
CHARACTER (LEN=lentime) :: ciodatetime
INTEGER :: i, ic, icon, idat, idesctype, ifil, ii, iicon, iiloc, isec, itype, &
         & i1, i2, j, jc, jj, jjloc, l, ll, nodat, noread, novars, npcc
REAL (KIND=kndrlong) :: alpha, ratio1, ratio2
INTEGER, DIMENSION(3) :: tlims
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: defudatobu, defudatoby, defvdatobv,&
                                         &  defvdatobx, defzdatobu, defzdatobv 
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: harudatobu, harudatoby, harvdatobv,&
                                         &  harvdatobx, harzdatobu, harzdatobv
REAL, DIMENSION(nqsecobu) :: uwsum
REAL, DIMENSION(nqsecobv) :: vwsum
REAL, DIMENSION(nobu,nqsecobu) :: uweights
REAL, DIMENSION(nobv,nqsecobv) :: vweights
REAL (KIND=kndrlong), DIMENSION(nconobc) :: phase
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:) :: phase_obc_d
 
!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*defudatobu* LOGICAL .TRUE. if current/transport/discharge data are externally
!                     specified at U-o.b.
!*defudatoby* LOGICAL .TRUE. if X-transport data are externally specified at
!                      Y-o.b.
!*defvdatobv* LOGICAL .TRUE. if current/transport/discharge data are externally
!                     specified at V-o.b.
!*defvdatobx* LOGICAL .TRUE. if Y-transport data are externally specified at
!                     X-o.b.
!*defzdatobu* LOGICAL .TRUE. if elevation data are externally specified at
!                     U-o.b.
!*defzdatobv* LOGICAL .TRUE. if elevation data are externally specified at
!                     V-o.b.
!*harudatobu* LOGICAL .TRUE. if current/transport/discharge data are specified
!                     in harmonic form at U-o.b.
!*harudatoby* LOGICAL .TRUE. if transport data are specified in harmonic form
!                     at Y-o.b.
!*harvdatobv* LOGICAL .TRUE. if current/transport/discharge data are specified
!                     in harmonic form at V-o.b.
!*harvdatobx* LOGICAL .TRUE. if transport data are specified in harmonic form
!                     at X-o.b.
!*harzdatobu* LOGICAL .TRUE. if elevation data are specified in harmonic form
!                     at U-o.b.
!*harzdatobv* LOGICAL .TRUE. if elevation data are specified in harmonic form
!                     at V-o.b.
!
!------------------------------------------------------------------------------
!


procname(pglev+1) = 'update_2dobc_data'
CALL log_timer_in(npcc)

idesctype = MERGE(1,2,iddesc.EQ.io_2uvobc)

!
!1. Initialise on first call
!---------------------------
!

IF (nt.EQ.0) THEN

!
!1.1 Allocate
!------------
!
!1.1.1 Phases
!------------
!

   IF (idesctype.EQ.1) THEN
      ALLOCATE(phase_obc_d(nconobc),STAT=errstat)
      CALL error_alloc('phase_obc_d',1,(/nconobc/),kndrlong)
      IF (nconobc.GT.0) phase_obc_d = phase_obc

!
!1.1.2 Arrays at (U,V)-open boundaries
!-------------------------------------
!
!     ---arrays defining o.b. points where external data are defined
      ALLOCATE (defudatobu(nobu),STAT=errstat)
      CALL error_alloc('defudatobu',1,(/nobu/),kndlog)
      IF (nobu.GT.0) defudatobu = .FALSE.
      ALLOCATE (defvdatobv(nobv),STAT=errstat)
      CALL error_alloc('defvdatobv',1,(/nobv/),kndlog)
      IF (nobv.GT.0) defvdatobv = .FALSE.
      ALLOCATE (defzdatobu(nobu),STAT=errstat)
      CALL error_alloc('defzdatobu',1,(/nobu/),kndlog)
      IF (nobu.GT.0) defzdatobu = .FALSE.
      ALLOCATE (defzdatobv(nobv),STAT=errstat)
      CALL error_alloc('defzdatobv',1,(/nobv/),kndlog)
      IF (nobv.GT.0) defzdatobv = .FALSE.

!     ---arrays defining o.b. points where harmonic external data are defined
      ALLOCATE (harudatobu(nobu),STAT=errstat)
      CALL error_alloc('harudatobu',1,(/nobu/),kndlog)
      IF (nobu.GT.0) harudatobu = .FALSE.
      ALLOCATE (harvdatobv(nobv),STAT=errstat)
      CALL error_alloc('harvdatobv',1,(/nobv/),kndlog)
      IF (nobv.GT.0) harvdatobv = .FALSE.
      ALLOCATE (harzdatobu(nobu),STAT=errstat)
      CALL error_alloc('harzdatobu',1,(/nobu/),kndlog)
      IF (nobu.GT.0) harzdatobu = .FALSE.
      ALLOCATE (harzdatobv(nobv),STAT=errstat)
      CALL error_alloc('harzdatobv',1,(/nobv/),kndlog)
      IF (nobv.GT.0) harzdatobv = .FALSE.

!
!1.1.3 Arrays at (X,Y)-open boundaries
!-------------------------------------
!

   ELSEIF (idesctype.EQ.2) THEN

!     ---arrays defining o.b. points where external data are defined
      ALLOCATE (defudatoby(noby),STAT=errstat)
      CALL error_alloc('defudatoby',1,(/noby/),kndlog)
      IF (noby.GT.0) defudatoby = .FALSE.
      ALLOCATE (defvdatobx(nobx),STAT=errstat)
      CALL error_alloc('defvdatobx',1,(/nobx/),kndlog)
      IF (nobx.GT.0) defvdatobx = .FALSE.

!     ---arrays defining o.b. points where harmonic external data are defined
      ALLOCATE (harudatoby(noby),STAT=errstat)
      CALL error_alloc('harudatoby',1,(/noby/),kndlog)
      IF (noby.GT.0) harudatoby = .FALSE.
      ALLOCATE (harvdatobx(nobx),STAT=errstat)
      CALL error_alloc('harvdatobx',1,(/nobx/),kndlog)
      IF (nobx.GT.0) harvdatobx = .FALSE.

   ENDIF

!
!1.2 Initialise
!--------------
!
!1.2.1 Arrays at (U,V)-open boundaries
!-------------------------------------
!

   IF (idesctype.EQ.1) THEN

      IF (nobu.GT.0) THEN
         defudatobu = .FALSE.; defzdatobu = .FALSE.
         harudatobu = .FALSE.; harzdatobu = .FALSE.
      ENDIF
      IF (nobv.GT.0) THEN
         defvdatobv = .FALSE.; defzdatobv = .FALSE.
         harvdatobv = .FALSE.; harzdatobv = .FALSE.
      ENDIF

!
!1.2.2 Arrays at (X,Y)-open boundaries
!-------------------------------------
!

   ELSEIF (idesctype.EQ.2) THEN

      IF (noby.GT.0) THEN
         defudatoby = .FALSE.; harudatoby = .FALSE.
      ENDIF
      IF (nobx.GT.0) THEN
         defvdatobx = .FALSE.; harvdatobx = .FALSE.
      ENDIF

   ENDIF

!
!1.3 Check whether external/harmonic data are in harmonic form
!-------------------------------------------------------------
!

   IF (nconobc.GT.0) THEN

!
!1.3.1 (U,V)-open boundaries
!---------------------------
!

      IF (idesctype.EQ.1) THEN

         ii_1311: DO ii=1,nobu
            IF (ANY(udatobu_amp(ii,:).NE.0.0)) THEN
               defudatobu(ii) = .TRUE.; harudatobu(ii) = .TRUE.
            ENDIF
            IF (ANY(zdatobu_amp(ii,:).NE.0.0)) THEN
               defzdatobu(ii) = .TRUE.; harzdatobu(ii) = .TRUE.
            ENDIF
         ENDDO ii_1311

         jj_1312: DO jj=1,nobv
            IF (ANY(vdatobv_amp(jj,:).NE.0.0)) THEN
               defvdatobv(jj) = .TRUE.; harvdatobv(jj) = .TRUE.
            ENDIF
            IF (ANY(zdatobv_amp(jj,:).NE.0.0)) THEN
               defzdatobv(jj) = .TRUE.; harzdatobv(jj) = .TRUE.
            ENDIF
         ENDDO jj_1312

!
!1.3.2 (X,Y)-open boundaries
!---------------------------
!

      ELSEIF (idesctype.EQ.2) THEN

         ii_1321: DO ii=1,nobx
            IF (ityp2dobx(ii).EQ.2.AND.ANY(vdatobx_amp(ii,:).NE.0.0)) THEN
               defvdatobx(ii) = .TRUE.; harvdatobx(ii) = .TRUE.
            ENDIF
         ENDDO ii_1321

         jj_1322: DO jj=1,noby
            IF (ityp2doby(jj).EQ.2.AND.ANY(udatoby_amp(jj,:).NE.0.0)) THEN
               defudatoby(jj) = .TRUE.; harudatoby(jj) = .TRUE.
            ENDIF
         ENDDO jj_1322

      ENDIF

   ENDIF

!
!1.4 Check if external data are provided
!---------------------------------------
!
!1.4.1 (U,V)-open boundaries
!---------------------------
!

   IF (idesctype.EQ.1) THEN

      ifil_141: DO ifil=2,nofiles
         IF (modfiles(iddesc,ifil,1)%defined) THEN
            nodat = no2dobc(ifil); itype = iobc2dtype(ifil)
            idat_1411: DO idat=1,nodat
               l = index2dobuv(idat,ifil)
               IF (l.GT.0.AND.l.LE.nobu) THEN
                  ii = l
                  SELECT CASE (itype)
                     CASE (1)
                        defudatobu(ii) = .TRUE.; defzdatobu(ii) = .TRUE.
                     CASE (2)
                        defzdatobu(ii) = .TRUE.
                     CASE (3)
                        defudatobu(ii) = .TRUE.
                  END SELECT
               ELSEIF (l.GT.0.AND.l.LE.(nobu+nobv)) THEN
                  jj = l - nobu
                  SELECT CASE (itype)
                     CASE (1)
                        defvdatobv(jj) = .TRUE.; defzdatobv(jj) = .TRUE.
                     CASE (2)
                        defzdatobv(jj) = .TRUE.
                     CASE (3)
                        defvdatobv(jj) = .TRUE.
                  END SELECT
               ENDIF
            ENDDO idat_1411
         ENDIF
      ENDDO ifil_141

!
!1.4.2 (X,Y)-open boundaries
!---------------------------
!

   ELSEIF (idesctype.EQ.2) THEN

      ifil_142: DO ifil=2,nofiles
         IF (modfiles(iddesc,ifil,1)%defined) THEN
            nodat = no2dobc(ifil)
            idat_1421: DO idat=1,nodat
               l = index2dobxy(idat,ifil)
               IF (l.GT.0.AND.l.LE.nobx) THEN
                  ii = l
                  defvdatobx(ii) = .TRUE.
               ELSEIF (l.GT.0.AND.l.LE.(nobx+noby)) THEN
                  jj = l - nobx
                  defudatoby(jj) = .TRUE.
               ENDIF
            ENDDO idat_1421
         ENDIF
      ENDDO ifil_142

   ENDIF

!
!1.5 Check whether data are available if needed
!----------------------------------------------
!

   IF (errchk) THEN

!
!1.5.1 (U,V)-open boundaries
!---------------------------
!

      IF (idesctype.EQ.1) THEN

!
!1.5.1.1 U-nodes
!---------------
!

         ii_1511: DO ii=1,nobu

!           ---currents (depth integrated, depth mean, discharge)
            SELECT CASE (ityp2dobu(ii))
               CASE (4,8,11,14:16)
                  IF (.NOT.defudatobu(ii)) THEN
                     nerrs = nerrs + 1
                     IF (nerrs.LE.maxerrors) THEN
                        WRITE (cpos,'(I12)') ii
                        WRITE (ctyp,'(I12)') ityp2dobu(ii)
                        cpos = ADJUSTL(cpos); ctyp = ADJUSTL(ctyp)
                        WRITE (ioerr,'(A)') 'No current data defined at&
                                          & U-open boundary point: '//TRIM(cpos)
                        WRITE (ioerr,'(A)') 'Type of open boundary condition: '&
                                          & //TRIM(ctyp)
                     ENDIF
                  ENDIF
            END SELECT

!           ---elevation data
            SELECT CASE (ityp2dobu(ii))
               CASE (2:3,8:9,11:12,17)
                  IF (.NOT.defzdatobu(ii)) THEN
                     nerrs = nerrs + 1
                     IF (nerrs.LE.maxerrors) THEN
                        WRITE (cpos,'(I12)') ii
                        WRITE (ctyp,'(I12)') ityp2dobu(ii)
                        cpos = ADJUSTL(cpos); ctyp = ADJUSTL(ctyp)
                        WRITE (ioerr,'(A)') 'No elevation data defined at&
                                          & U-open boundary point: '//TRIM(cpos)
                        WRITE (ioerr,'(A)') 'Type of open boundary condition: '&
                                          & //TRIM(ctyp)
                     ENDIF
                  ENDIF
                  
            END SELECT

         ENDDO ii_1511

!
!1.5.1.2 V-nodes
!---------------
!

         jj_1512: DO jj=1,nobv

!           ---currents (depth integrated, depth mean, discharge)
            SELECT CASE (ityp2dobv(jj))
               CASE (4,8,11,14:16)
                  IF (.NOT.defvdatobv(jj)) THEN
                     nerrs = nerrs + 1
                     IF (nerrs.LE.maxerrors) THEN
                        WRITE (cpos,'(I12)') jj
                        WRITE (ctyp,'(I12)') ityp2dobv(jj)
                        cpos = ADJUSTL(cpos); ctyp = ADJUSTL(ctyp)
                        WRITE (ioerr,'(A)') 'No current data defined at&
                                          & V-open boundary point: '//TRIM(cpos)
                        WRITE (ioerr,'(A)') 'Type of open boundary condition: '&
                                          & //TRIM(ctyp)
                     ENDIF
                  ENDIF
            END SELECT

!           ---elevation data
            SELECT CASE (ityp2dobv(jj))
               CASE (2:3,8:9,11:12,17)
                  IF (.NOT.defzdatobv(jj)) THEN
                     nerrs = nerrs + 1
                     IF (nerrs.LE.maxerrors) THEN
                        WRITE (cpos,'(I12)') jj
                        WRITE (ctyp,'(I12)') ityp2dobv(jj)
                        cpos = ADJUSTL(cpos); ctyp = ADJUSTL(ctyp)
                        WRITE (ioerr,'(A)') 'No elevation data defined at&
                                          & V-open boundary point: '//TRIM(cpos)
                        WRITE (ioerr,'(A)') 'Type of open boundary condition: '&
                                          & //TRIM(ctyp)
                     ENDIF
                  ENDIF
            END SELECT

         ENDDO jj_1512

!
!1.5.2 (X,Y)-open boundaries
!---------------------------
!

      ELSEIF (idesctype.EQ.2) THEN

!        ---X-nodes
         ii_1521: DO ii=1,nobx
            IF (ityp2dobx(ii).EQ.2.AND.(.NOT.defvdatobx(ii))) THEN
               nerrs = nerrs + 1
               IF (nerrs.LE.maxerrors) THEN
                  WRITE (cpos,'(I12)') ii; WRITE (ctyp,'(I12)') ityp2dobx(ii)
                  cpos = ADJUSTL(cpos); ctyp = ADJUSTL(ctyp)
                  WRITE (iowarn,'(A)') 'No tangential current data defined at&
                                      & X-open boundary point: '//TRIM(cpos)
                  WRITE (iowarn,'(A)') 'Type of open boundary condition: '&
                                     & //TRIM(ctyp)
               ENDIF
            ENDIF
         ENDDO ii_1521

!        ---Y-nodes
         jj_1522: DO jj=1,noby
            IF (ityp2doby(jj).EQ.2.AND.(.NOT.defudatoby(jj))) THEN
               nerrs = nerrs + 1
               IF (nerrs.LE.maxerrors) THEN
                  WRITE (cpos,'(I12)') jj; WRITE (ctyp,'(I12)') ityp2doby(jj)
                  cpos = ADJUSTL(cpos); ctyp = ADJUSTL(ctyp)
                  WRITE (ioerr,'(A)') 'No current data defined at&
                                     & Y-open boundary point: '//TRIM(cpos)
                  WRITE (ioerr,'(A)') 'Type of open boundary condition: '&
                                    & //TRIM(ctyp)
               ENDIF
            ENDIF
         ENDDO jj_1522

      ENDIF

   ENDIF

   CALL error_abort('update_2dobc_data',ierrno_inival)

ENDIF

!
!2. Update data
!--------------
!

ifil_200: DO ifil=2,nofiles

   IF (modfiles(iddesc,ifil,1)%defined) THEN

      tlims = modfiles(iddesc,ifil,1)%tlims
      no_time_series = lim_dims(ABS(tlims)).EQ.1

      IF ((nt.GT.tlims(1).AND.no_time_series).OR.&
        & (.NOT.loop_index(tlims,nt)).OR.&
         & modfiles(iddesc,ifil,1)%iostat.EQ.3) CYCLE ifil_200

!
!2.1 Set parameters
!------------------
!

      nodat = no2dobc(ifil)
      IF (idesctype.EQ.1) THEN
         novars = MERGE(1,2,iobc2dtype(ifil).GT.1)
         SELECT CASE (iobc2dtype(ifil))
            CASE (1); i1 = 1; i2 = 2
            CASE (2); i1 = 0; i2 = 1
            CASE (3); i1 = 1; i2 = 0
         END SELECT
      ELSE
         novars = 1
         i1 = 1; i2 = 0
      ENDIF

!
!2.2 Time-independent data or regular data without time interpolation
!--------------------------------------------------------------------
!

      IF (no_time_series.OR.&
       & (modfiles(iddesc,ifil,1)%time_regular.AND.tlims(3).LT.0)) THEN
 
         IF (nt.EQ.tlims(1)) THEN
            ciodatetime = cdatetime_undef
            DO WHILE (ciodatetime.NE.CDateTime)
               CALL define_2dobc_data(iddesc,ifil,ciodatetime,&
                                    & datavals2d(1:nodat,1:novars,ifil,1),&
                                    & nodat,novars)
               IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
                  CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
               ENDIF
            ENDDO
         ELSE
            CALL define_2dobc_data(iddesc,ifil,ciodatetime,&
                                 & datavals2d(1:nodat,1:novars,ifil,1),&
                                 & nodat,novars)
            IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) CYCLE ifil_200
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
            CALL define_2dobc_data(iddesc,ifil,ciodatetime,&
                                 & datavals2d(1:nodat,1:novars,ifil,2),&
                                 & nodat,novars)
            IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
               CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
            ENDIF
            CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat(ifil,2))
            IF (nosecsdat(ifil,2).GT.nosecsrun) GOTO 1002
            noread = 1
         ENDIF

!
!2.3.2 Update data
!-----------------
!

         IF (nt.LT.nstep) THEN
            DO WHILE (nosecsdat(ifil,2).LE.nosecsrun)
!              ---store old data
               datavals2d(1:nodat,1:novars,ifil,1) = &
             & datavals2d(1:nodat,1:novars,ifil,2)
               nosecsdat(ifil,1) = nosecsdat(ifil,2)
!              ---read new data
               CALL define_2dobc_data(iddesc,ifil,ciodatetime,&
                                   & datavals2d(1:nodat,1:novars,ifil,2),&
                                   & nodat,novars)
               noread = noread + 1
               IF (noread.EQ.2.AND.&
                 & modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
                  CALL error_file(ierrno_fend,&
                                & filepars=modfiles(iddesc,ifil,1))
               ENDIF
               CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat(ifil,2))
            ENDDO
         ENDIF

      ENDIF

!
!2.4 Evaluate data at current time
!---------------------------------
!
!2.4.1 Normal conditions
!-----------------------
!
!     ---without time interpolation
      IF (no_time_series.OR.tlims(3).LT.0) THEN
         l_241: DO l=1,nodat
            IF (idesctype.EQ.1) THEN
               ll = index2dobuv(l,ifil)
               IF (ll.LE.nobu) THEN
                  IF (i1.GT.0) udatobu(ll,1) = datavals2d(l,i1,ifil,1)
                  IF (i2.GT.0) zdatobu(ll,1) = datavals2d(l,i2,ifil,1)
               ELSE
                  ll = ll - nobu
                  IF (i1.GT.0) vdatobv(ll,1) = datavals2d(l,i1,ifil,1)
                  IF (i2.GT.0) zdatobv(ll,1) = datavals2d(l,i2,ifil,1)

               ENDIF
            ELSE
               ll = index2dobxy(l,ifil)
               IF (ll.LE.nobx) THEN
                  vdatobx(ll,1) = datavals2d(l,1,ifil,1)
               ELSE
                  ll = ll - nobx
                  udatoby(ll,1) = datavals2d(l,1,ifil,1)
               ENDIF
            ENDIF
         ENDDO l_241

!     ---with time interpolation
      ELSE
         ratio2 = MERGE(0.0,(nosecsrun-nosecsdat(ifil,1))&
                     & /REAL(nosecsdat(ifil,2)-nosecsdat(ifil,1)),tlims(3).LT.0)
         ratio1 = 1.0 - ratio2
         l_242: DO l=1,nodat
            IF (idesctype.EQ.1) THEN
               ll = index2dobuv(l,ifil)
               IF (ll.LE.nobu) THEN
                  IF (i1.GT.0) THEN
                       udatobu(ll,1) = ratio1*datavals2d(l,i1,ifil,1) + &
                                     & ratio2*datavals2d(l,i1,ifil,2)
                  ENDIF
                  IF (i2.GT.0) THEN
                       zdatobu(ll,1) = ratio1*datavals2d(l,i2,ifil,1) + &
                                     & ratio2*datavals2d(l,i2,ifil,2)
                  ENDIF
               ELSE
                  ll = ll - nobu
                  IF (i1.GT.0) THEN
                       vdatobv(ll,1) = ratio1*datavals2d(l,i1,ifil,1) + &
                                     & ratio2*datavals2d(l,i1,ifil,2)
                  ENDIF
                  IF (i2.GT.0) THEN
                       zdatobv(ll,1) = ratio1*datavals2d(l,i2,ifil,1) + &
                                     & ratio2*datavals2d(l,i2,ifil,2)
                  ENDIF
               ENDIF
            ELSE
               ll = index2dobxy(l,ifil)
               IF (ll.LE.nobx) THEN
                  vdatobx(ll,1) =  ratio1*datavals2d(l,1,ifil,1) + &
                                 & ratio2*datavals2d(l,1,ifil,2)
               ELSE
                  ll = ll - nobx
                  udatoby(ll,1) = ratio1*datavals2d(l,1,ifil,1) + &
                                & ratio2*datavals2d(l,1,ifil,2)
               ENDIF
            ENDIF
         ENDDO l_242
      ENDIF

   ENDIF

ENDDO ifil_200

!
!3. Update nodal corrections and phases
!--------------------------------------
!

IF (idesctype.EQ.1.AND.nconobc.GT.0) THEN
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
ENDIF

!
!4. Prescribed transports and elevations
!---------------------------------------
!
!4.1 Normal conditions
!---------------------
!

IF (idesctype.EQ.1) THEN

!
!4.1.1 U-nodes
!-------------
!

   zdatobu_old = zdatobu(:,2)

   iiloc_411: DO iiloc=1,nobuloc_ext
      i = iobuloc(iiloc); j = jobuloc(iiloc)
      ii = indexobu(iiloc)
      IF (westobu(ii).AND.i.EQ.ncloc+1) CYCLE iiloc_411
      IF (defudatobu(ii)) THEN
         IF (harudatobu(ii)) THEN 
            phase = MOD(phase_obc_d-udatobu_pha(ii,:),twopi_d)
            udatobu(ii,2) = udatobu(ii,1) + &
                          & SUM(fnode_obc*udatobu_amp(ii,:)*COS(phase))
         ELSE
            udatobu(ii,2) = udatobu(ii,1)
         ENDIF
      ENDIF
      IF (defzdatobu(ii)) THEN
         IF (harzdatobu(ii)) THEN 
            phase = MOD(phase_obc_d-zdatobu_pha(ii,:),twopi_d)
            zdatobu(ii,2) = zdatobu(ii,1) + &
                      & SUM(fnode_obc*zdatobu_amp(ii,:)*COS(phase))
         ELSE
            zdatobu(ii,2) = zdatobu(ii,1)
         ENDIF
         IF (iopt_obc_invbar.EQ.1.AND.ityp2dobu(ii).NE.17) THEN
            ic = MERGE(i,i-1,westobu(ii))
            zdatobu(ii,2) = zdatobu(ii,2) + (atmpres_ref-atmpres(ic,j))&
                                         & /(density_ref*gaccatu(i,j))
         ENDIF
      ENDIF
   ENDDO iiloc_411

!
!4.1.2 V-nodes
!-------------
!

   zdatobv_old = zdatobv(:,2)

   jjloc_412: DO jjloc=1,nobvloc_ext
      i = iobvloc(jjloc); j = jobvloc(jjloc)
      jj = indexobv(jjloc)
      IF (soutobv(jj).AND.j.EQ.nrloc+1) CYCLE jjloc_412
      IF (defvdatobv(jj)) THEN
         IF (harvdatobv(jj)) THEN 
            phase = MOD(phase_obc_d-vdatobv_pha(jj,:),twopi_d)
            vdatobv(jj,2) = vdatobv(jj,1) + &
                          & SUM(fnode_obc*vdatobv_amp(jj,:)*COS(phase))
         ELSE
            vdatobv(jj,2) = vdatobv(jj,1)
         ENDIF
      ENDIF
      IF (defzdatobv(jj)) THEN
         IF (harzdatobv(jj)) THEN 
            phase = MOD(phase_obc_d-zdatobv_pha(jj,:),twopi_d)
            zdatobv(jj,2) = zdatobv(jj,1) + &
                          & SUM(fnode_obc*zdatobv_amp(jj,:)*COS(phase))
         ELSE
            zdatobv(jj,2) = zdatobv(jj,1)
         ENDIF
         IF (iopt_obc_invbar.EQ.1.AND.ityp2dobv(jj).NE.17) THEN
            jc = MERGE(j,j-1,soutobv(jj))
            zdatobv(jj,2) = zdatobv(jj,2) + (atmpres_ref-atmpres(i,jc))&
                                         & /(density_ref*gaccatv(i,j))
         ENDIF
      ENDIF
   ENDDO jjloc_412

!
!4.2 Tangential conditions
!-------------------------
!

ELSE

!
!4.2.1 X-nodes
!-------------
!

   iiloc_421: DO iiloc=1,nobxloc
      ii = indexobx(iiloc)
      IF (defvdatobx(ii)) THEN
         IF (harvdatobx(ii)) THEN 
            phase = MOD(phase_obc_d-vdatobx_pha(ii,:),twopi_d)
            vdatobx(ii,2) = vdatobx(ii,1) + &
                          & SUM(fnode_obc*vdatobx_amp(ii,:)*COS(phase))
         ELSE
            vdatobx(ii,2) = vdatobx(ii,1)
         ENDIF
      ENDIF
   ENDDO iiloc_421

!
!4.2.2 Y-nodes
!-------------
!

   jjloc_422: DO jjloc=1,nobyloc
      jj = indexoby(jjloc)
      IF (defudatoby(jj)) THEN
         IF (harudatoby(jj)) THEN 
            phase = MOD(phase_obc_d-udatoby_pha(jj,:),twopi_d)
            udatoby(jj,2) = udatoby(jj,1) + &
                          & SUM(fnode_obc*udatoby_amp(jj,:)*COS(phase))
         ELSE
            udatoby(jj,2) = udatoby(jj,1)
         ENDIF
      ENDIF
   ENDDO jjloc_422

ENDIF

!
!5. Discharge conditions along sections
!--------------------------------------
!
!5.1 U-nodes
!-----------
!

IF (nqsecobu.GT.0) THEN

!  ---weight factors
   uweights = 0.0
   iiloc_511: DO iiloc=1,nobuloc
      ii = indexobu(iiloc)
      IF (ityp2dobu(ii).EQ.16) THEN
         i = iobuloc(iiloc); j = jobuloc(iiloc)
         isec_5111: DO isec=1,nqsecobu
            IF (ii.GE.iqsecobu(isec,1).AND.ii.LE.iqsecobu(isec,2)) THEN
               uweights(ii,isec) = SQRT(gaccatu(i,j)/bdragcoefatu(i,j))* &
                                      & delyatu(i,j)*deptotatu(i,j)**1.5
            ENDIF
         ENDDO isec_5111
      ENDIF
   ENDDO iiloc_511

!  ---combine on all processes
   IF (iopt_MPI.EQ.1) THEN
      CALL combine_stats_glb(uweights,nobu,nobuprocs,indexobuprocs,0,&
                           & commall=.TRUE.)
   ENDIF

!  ---weights sum for each section
   isec_512: DO isec=1,nqsecobu
      uwsum(isec) = SUM(uweights(:,isec),MASK=ityp2dobu.EQ.16)
   ENDDO isec_512

!  ---discharge value
   iiloc_513: DO iiloc=1,nobuloc
      ii = indexobu(iiloc)
      IF (ityp2dobu(ii).EQ.16) THEN
         isec_5131: DO isec=1,nqsecobu
            IF (ii.GE.iqsecobu(isec,1).AND.ii.LE.iqsecobu(isec,2)) THEN
               udatobu(ii,2) = (uweights(ii,isec)/uwsum(isec))*udatobu(ii,2)
            ENDIF
         ENDDO isec_5131
      ENDIF
   ENDDO iiloc_513

ENDIF

!
!5.2 V-nodes
!-----------
!

IF (nqsecobv.GT.0) THEN

!  ---weight factors
   vweights = 0.0
   jjloc_521: DO jjloc=1,nobvloc
      jj = indexobv(jjloc)
      IF (ityp2dobv(jj).EQ.16) THEN
         i = iobvloc(jjloc); j = jobvloc(jjloc)
         isec_5211: DO isec=1,nqsecobv
            IF (jj.GE.jqsecobv(isec,1).AND.ii.LE.jqsecobv(isec,2)) THEN
               vweights(jj,isec) = SQRT(gaccatv(i,j)/bdragcoefatv(i,j))* &
                                      & delxatv(i,j)*deptotatv(i,j)**1.5
            ENDIF
         ENDDO isec_5211
      ENDIF
   ENDDO jjloc_521

!  ---combine on all processes
   IF (iopt_MPI.EQ.1) THEN
      CALL combine_stats_glb(vweights,nobv,nobvprocs,indexobvprocs,0,&
                           & commall=.TRUE.)
   ENDIF

!  ---weights sum for each section
   isec_522: DO isec=1,nqsecobv
      vwsum(isec) = SUM(vweights(:,isec),MASK=ityp2dobv.EQ.16)
   ENDDO isec_522

!  ---discharge value
   jjloc_523: DO jjloc=1,nobvloc
      jj = indexobv(jjloc)
      IF (ityp2dobv(jj).EQ.16) THEN
         isec_5231: DO isec=1,nqsecobv
            IF (jj.GE.jqsecobv(isec,1).AND.jj.LE.jqsecobv(isec,2)) THEN
               vdatobv(jj,2) = (vweights(jj,isec)/vwsum(isec))*vdatobv(jj,2)
            ENDIF
         ENDDO isec_5231
      ENDIF
   ENDDO jjloc_523

ENDIF

!
!6. Apply relaxation condition
!-----------------------------
!

IF (ntobcrlx.GT.0.AND.nt.LT.ntobcrlx) THEN
   alpha = nt/REAL(ntobcrlx)

!
!6.1 Normal conditions
!---------------------
!

   IF (idesctype.EQ.1) THEN

!     ---U-nodes
      iiloc_611: DO iiloc=1,nobuloc
         ii = indexobu(iiloc)
         IF (defudatobu(ii)) udatobu(ii,2) = alpha*udatobu(ii,2)
         IF (defzdatobu(ii)) zdatobu(ii,2) = alpha*zdatobu(ii,2)
      ENDDO iiloc_611

!     ---V-nodes
      jjloc_612: DO jjloc=1,nobvloc
         jj = indexobv(jjloc)
         IF (defvdatobv(jj)) vdatobv(jj,2) = alpha*vdatobv(jj,2)
         IF (defzdatobv(jj)) zdatobv(jj,2) = alpha*zdatobv(jj,2)
      ENDDO jjloc_612

!
!6.2 Tangential conditions
!-------------------------
!

   ELSE

!  ---X-nodes
      iiloc_621: DO iiloc=1,nobxloc
         ii = indexobx(iiloc)
         IF (defvdatobx(ii)) vdatobx(ii,2) = alpha*vdatobx(ii,2)
      ENDDO iiloc_621

!  ---Y-nodes
      jjloc_622: DO jjloc=1,nobyloc
         jj = indexoby(jjloc)
         IF (defudatoby(jj)) udatoby(jj,2) = alpha*udatoby(jj,2)
      ENDDO jjloc_622

   ENDIF

ENDIF

!
!7. Deallocate arrays
!--------------------
!

IF (cold_start.OR.nt.GT.0.AND.nt.EQ.hydro_end) THEN
   IF ((idesctype.EQ.1.AND.iopt_obc_2d_tang.EQ.0).OR.idesctype.EQ.2) THEN
      DEALLOCATE (phase_obc_d)
   ENDIF
   IF (idesctype.EQ.1) THEN
      DEALLOCATE (defudatobu,defvdatobv,defzdatobu,defzdatobv)
      DEALLOCATE (harudatobu,harvdatobv,harzdatobu,harzdatobv)
   ELSE
      DEALLOCATE (defudatoby,defvdatobx,harudatoby,harvdatobx)
   ENDIF
ENDIF

CALL log_timer_out(npcc,itm_bconds)


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

END SUBROUTINE update_2dobc_data

!========================================================================

SUBROUTINE write_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
!************************************************************************
!
! *write_2dobc_data* Write data for 2-D normal or tangential open boundary
!                    conditions in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_2D.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - define_2dobc_data
!
! Module calls - error_alloc_struc, later, open_filepars,
!                output_flag, set_modfiles_atts, set_modvars_atts,
!                write_atts_mod, write_time, write_vars
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
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
REAL, INTENT(IN), DIMENSION(nodat,novars) :: data2d

!
! Name         Type Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER   Data file id (io_2uvobc or io_2xyobc)
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*data2d*      REAL    Output data
!*nodat*       INTEGER Number of output data
!*novars*      INTEGER Number of output parameters
!*iddesc*      INTEGER Data file id
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

procname(pglev+1) = 'write_2dobc_data'
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
!2. Write header on first call
!-----------------------------
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
   CALL set_modvars_atts(iddesc,ifil,2,varatts,numvars)

!  ---write attributes
   CALL write_atts_mod(filepars,varatts,numvars)

!  ---deallocate
   DEALLOCATE (varatts)

ENDIF

!
!3. Write data
!-------------
!
!---date/time
CALL write_time(ciodatetime,filepars)

!---data
vecids = (/(ivar,ivar=2,numvars)/)
CALL write_vars(data2d,filepars,0,vecids=vecids)

modfiles(iddesc,ifil,2) = filepars

1000 CALL log_timer_out()


RETURN

END SUBROUTINE write_2dobc_data

!========================================================================

SUBROUTINE write_2dobc_spec(iddesc)
!************************************************************************
!
! *write_2dobc_spec* Write specifier arrays for 2-D normal or tangential open
!                    boundary conditions in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_2D.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - define_2dobc_spec
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_vars
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE obconds
USE paralpars
USE switches
USE tide
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name      Type      Purpose
!-------------------------------------------------------------------------
!*iddesc*   INTEGER   Data file id
!
!-------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_2dobc_spec'
CALL log_timer_in()

!
!1. Write file header
!--------------------
!
!---file attributes
CALL set_modfiles_atts(iddesc,1,2)
filepars = modfiles(iddesc,1,2)
numvars = filepars%novars


!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL set_modvars_atts(iddesc,1,2,varatts,numvars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2. Write specifier arrays
!-------------------------
!

SELECT CASE (iddesc)

!
!2.1 Normal (U,V)-open boundaries
!--------------------------------
!

   CASE (io_2uvobc)

!
!2.1.1 Type of conditions
!------------------------
!

      varid = 0
!     ---U-nodes
      IF (nobu.GT.0) THEN
         varid = varid + 1
         CALL write_vars(ityp2dobu,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(iloczobu,filepars,varid,varatts=varatts)
      ENDIF

!     ---V-nodes
      IF (nobv.GT.0) THEN
         varid = varid + 1
         CALL write_vars(ityp2dobv,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(iloczobv,filepars,varid,varatts=varatts)
      ENDIF

!
!2.1.2 Discharge along open boundary sections
!--------------------------------------------
!

      IF (nqsecobu.GT.0) THEN
         varid = varid + 1
         CALL write_vars(iqsecobu,filepars,varid,varatts=varatts)
      ENDIF
      IF (nqsecobv.GT.0) THEN
         varid = varid + 1
         CALL write_vars(jqsecobv,filepars,varid,varatts=varatts)
      ENDIF

!
!2.1.3 Number of data and index maps per input file
!--------------------------------------------------
!

      IF (maxdatafiles(iddesc,1).GT.1) THEN
         varid = varid + 1
         CALL write_vars(no2dobuv,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(iobc2dtype,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(index2dobuv,filepars,varid,varatts=varatts)
      ENDIF
      
!
!2.1.4 Amplitudes and phases
!---------------------------
!
!     ---U-nodes
      IF (nobu.GT.0.AND.nconobc.GT.0) THEN
         varid = varid + 1
         CALL write_vars(udatobu_amp,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(zdatobu_amp,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(udatobu_pha,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(zdatobu_pha,filepars,varid,varatts=varatts)
      ENDIF

!     ---V-nodes
      IF (nobv.GT.0.AND.nconobc.GT.0) THEN
         varid = varid + 1
         CALL write_vars(vdatobv_amp,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(zdatobv_amp,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(vdatobv_pha,filepars,varid,varatts=varatts)
         varid = varid + 1
         CALL write_vars(zdatobv_pha,filepars,varid,varatts=varatts)
      ENDIF

!
!2.2 Tangential (X,Y)-open boundaries
!------------------------------------
!

   CASE (io_2xyobc)

!
!2.2.1 Type of conditions
!------------------------
!

      varid = 0
!     ---X-nodes
      IF (nobx.GT.0) THEN
         varid = varid + 1
         CALL write_vars(ityp2dobx,filepars,varid,varatts=varatts)
      ENDIF

!     ---Y-nodes
      IF (noby.GT.0) THEN
         varid = varid + 1
         CALL write_vars(ityp2doby,filepars,varid,varatts=varatts)
      ENDIF

!
!2.2.2 Number of data and index maps per input file
!--------------------------------------------------
!

     IF (maxdatafiles(iddesc,1).GT.1) THEN
        varid = varid + 1
        CALL write_vars(no2dobxy,filepars,varid,varatts=varatts)
        varid = varid + 1
        CALL write_vars(index2dobxy,filepars,varid,varatts=varatts)
     ENDIF

!
!2.2.3 Amplitudes and phases
!---------------------------
!
!    ---X-nodes
     IF (nobx.GT.0.AND.nconobc.GT.0) THEN
        varid = varid + 1
        CALL write_vars(vdatobx_amp,filepars,varid,varatts=varatts)
        varid = varid + 1
        CALL write_vars(vdatobx_pha,filepars,varid,varatts=varatts)
     ENDIF

!    ---Y-nodes
     IF (noby.GT.0.AND.nconobc.GT.0) THEN
        varid = varid + 1
        CALL write_vars(udatoby_amp,filepars,varid,varatts=varatts)
        varid = varid + 1
        CALL write_vars(udatoby_pha,filepars,varid,varatts=varatts)
     ENDIF

END SELECT

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(iddesc,1,2) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE write_2dobc_spec
