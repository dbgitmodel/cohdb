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
! *Open_Boundary_Data_Prof* Define open boundary conditions and data for
!                           scalar quantities at C-nodes or 3-D currents
!                           at U- or V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_Prof.f90  V2.10.1
!
! $Date: 2017-06-12 09:17:33 +0200 (Mon, 12 Jun 2017) $
!
! $Revision: 1033 $
!
! Description - 
!
! Reference -
!
! Subroutines - define_profobc_data, define_profobc_spec, read_profobc_data,
!               read_profobc_spec, update_profobc_data, write_profobc_data,
!               write_profobc_spec
!
!************************************************************************
!

!========================================================================

SUBROUTINE define_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
!************************************************************************
!
! *define_profobc_data* Obtain open boundary profiles
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_Prof.f90  V2.1.2
!
! Description - 
!
! Reference -
!
! Calling program - update_profobc_data
!
! External calls - read_profobc_data, usrdef_profobc_data, write_profobc_data
!
! Module calls - error_file, suspend_proc
!
!************************************************************************
!
USE gridpars
USE iopars
USE syspars
USE error_routines, ONLY: error_file
USE time_routines, ONLY: log_timer_in, log_timer_out, suspend_proc

IMPLICIT NONE

!
!*  Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, numprofs
REAL, INTENT(INOUT), DIMENSION(numprofs,nz) :: psiprofdat

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*psiprofdat*  REAL     Profile arrays
!*numprofs*    INTEGER  Number of profiles in data file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: define_data
INTEGER :: endfile, nosecs


procname(pglev+1) = 'define_profobc_data'
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
      CALL read_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
!  ---user-defined
   ELSEIF (modfiles(iddesc,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
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

!
!3. Obtain data
!--------------
!

define_data = .TRUE.

DO WHILE (define_data)

!  ---read
   IF (modfiles(iddesc,ifil,1)%status.EQ.'R') THEN
      CALL read_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
!  ---user-defined
   ELSEIF (modfiles(iddesc,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
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
   CALL write_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE define_profobc_data

!========================================================================

SUBROUTINE define_profobc_spec(iddesc,itypobux,itypobvy,iprofobux,iprofobvy,&
                             & iprofrlx,noprofsd,indexprof,indexvar,novars,&
                             & nofiles,nobux,nobvy)
!************************************************************************
!
! *define_profobc_spec* Define specifier arrays for profile open boundary
!                       conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_Prof.f90  V2.7
!
! Description - 
!
! Reference -
!
! Calling program - current_corr, salinity_equation, temperature_equation
!
! External calls - read_profobc_spec, usrdef_profobc_spec, write_profobc_spec
!
! Module calls - error_abort, error_limits_arr, error_value_var
!
!************************************************************************
!
USE iopars
USE paralpars
USE relaxation
USE switches
USE error_routines, ONLY: error_abort, error_limits_arr, error_value_var
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, nobux, nobvy, nofiles, novars
INTEGER, INTENT(INOUT), DIMENSION(2:nofiles) :: noprofsd
INTEGER, INTENT(INOUT), DIMENSION(nobux) :: itypobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy) :: itypobvy
INTEGER, INTENT(INOUT), DIMENSION(nobux,novars) :: iprofobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy,novars) :: iprofobvy
INTEGER, INTENT(INOUT), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexprof
INTEGER, INTENT(INOUT), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexvar
INTEGER, INTENT(INOUT), DIMENSION(norlxzones) :: iprofrlx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id
!*itypobux*  INTEGER Type of U- or X-open boundary condition
!*itypobvy*  INTEGER Type of V- or Y-open boundary condition
!*iprofobux* INTEGER Profile numbers at U- or X-open boundaries
!*iprofobvy* INTEGER Profile numbers at V- or Y-open boundaries
!*iprofrlx*  INTEGER Disables/enables relaxation at open boundary zones
!*noprofsd*  INTEGER Number of profiles per data file
!*indexprof* INTEGER Mapping array of the profile numbers in the data files to
!                    the profile numbers assigned to the open boundaries. The
!                    physical size of the first dimension equals the number of
!                    profiles in a data file.
!*indexvar*  INTEGER Defines the variable number of the profiles in a data
!                    file. The physical size of the first dimension equals the
!                    number of profiles in a data file.
!*novars*    INTEGER Total number of variables
!*nofiles*   INTEGER Number of data files (+1)
!*nobux*     INTEGER Number of nodes at U- or X-open boundaries
!*nobvy*     INTEGER Number of nodes at V- or Y-open boundaries
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cprof, cpvars, cvar
INTEGER :: ifil, ii, iprof, iprofdat, ivar, jj, maxtypes, ncount, noprofs, &
         & noprofsivar, npcc


procname(pglev+1) = 'define_profobc_spec'
CALL log_timer_in(npcc)

!
!1. Initialise arrays
!--------------------
!

IF (nobux.GT.0) THEN
   itypobux = 0; iprofobux = 0
ENDIF
IF (nobvy.GT.0) THEN
   itypobvy = 0; iprofobvy = 0
ENDIF

IF (norlxzones.GT.0) iprofrlx = 0

IF (nofiles.GT.1) THEN
   noprofsd = 0
   indexprof = 0
   indexvar = 0
ENDIF

!
!2. Obtain specifier arrays
!--------------------------
!
!---read
IF (modfiles(iddesc,1,1)%status.EQ.'R') THEN
   CALL read_profobc_spec(iddesc,itypobux,itypobvy,iprofobux,iprofobvy,&
                        & iprofrlx,noprofsd,indexprof,indexvar,novars,nofiles,&
                        & nobux,nobvy)
!---user-defined
ELSEIF (modfiles(iddesc,1,1)%status.EQ.'N') THEN
   CALL usrdef_profobc_spec(iddesc,itypobux,itypobvy,iprofobux,iprofobvy,&
                          & iprofrlx,noprofsd,indexprof,indexvar,novars,&
                          & nofiles,nobux,nobvy)
ENDIF

!
!3. Reset to defaults
!--------------------
!

IF (novars.EQ.1) indexvar = 1


!
!4. Check specifier arrays
!-------------------------
!

IF (master) THEN

   IF (nofiles.GT.1) noprofs = SUM(noprofsd)

   SELECT CASE (iddesc)
      CASE (io_3uvobc); maxtypes = 5
      CASE (io_3xyobc); maxtypes = 2
      CASE DEFAULT; maxtypes = 3
   END SELECT

!  ---U- or X-nodes
   IF (nofiles.GT.1) THEN
      ii_411: DO ii=1,nobux
         CALL error_limits_arr(itypobux(ii),'itypobux',0,maxtypes,1,indx=(/ii/))
         ivar_4111: DO ivar=1,novars
            CALL error_limits_arr(iprofobux(ii,ivar),'iprofobux',0,noprofs,2,&
                                & indx=(/ii,ivar/))
         ENDDO ivar_4111
      ENDDO ii_411
   ENDIF

!  ---V- or Y-nodes
   IF (nofiles.GT.1) THEN
      jj_412: DO jj=1,nobvy
         CALL error_limits_arr(itypobvy(jj),'itypobvy',0,maxtypes,1,indx=(/jj/))
         ivar_4121: DO ivar=1,novars
            CALL error_limits_arr(iprofobvy(jj,ivar),'iprofobvy',0,noprofs,2,&
                                & indx=(/jj,ivar/))
         ENDDO ivar_4121
      ENDDO jj_412
   ENDIF

!  ---discharge conditions
   IF (iddesc.EQ.io_3uvobc) THEN
      ii_413: DO ii=1,nobux
         IF (itypobux(ii).EQ.15.OR.itypobux(ii).EQ.16) THEN
            CALL error_value_var(itypobux(ii),'itypobu',5)
         ENDIF
      ENDDO ii_413
      jj_414: DO jj=1,nobvy
         IF (itypobvy(jj).EQ.15.OR.itypobvy(jj).EQ.16) THEN
            CALL error_value_var(itypobvy(jj),'itypobv',5)
         ENDIF
      ENDDO jj_414
   ENDIF

!  ---no missing profiles at open boundaries
   ivar_414: DO ivar=1,novars
      noprofsivar = MAX(MAXVAL(iprofobux(:,ivar)),MAXVAL(iprofobvy(:,ivar)))
      iprof_4141: DO iprof=1,noprofsivar
         IF (ALL(iprofobux(:,ivar).NE.iprof).AND.&
           & ALL(iprofobvy(:,ivar).NE.iprof)) THEN
            nerrs = nerrs + 1
            IF (master.AND.errchk.AND.nerrs.LE.maxerrors) THEN
               WRITE (cprof,'(I12)') iprof; cprof = ADJUSTL(cprof)
               WRITE (cvar,'(I12)') ivar; cvar = ADJUSTL(cvar)
               WRITE (cpvars,'(I12)') noprofsivar; cpvars = ADJUSTL(cpvars) 
               WRITE (ioerr,'(A)') 'Profile number '//TRIM(cprof)//&
                                 & ' for variable '//TRIM(cvar)//&
                                 & ' is not defined at any open boundary'
               WRITE (ioerr,'(A)') 'All numbers between 1 and '//&
                                & TRIM(cpvars)//' must be defined at least once'
            ENDIF
         ENDIF
      ENDDO iprof_4141
   ENDDO ivar_414

 !  ---variable number in a data profile
   ifil_415: DO ifil=2,nofiles
   iprofdat_415: DO iprofdat=1,noprofsd(ifil)
      CALL error_limits_arr(indexvar(iprofdat,ifil),'indexvar',1,novars,2,&
                          & indx=(/iprofdat,ifil/))
   ENDDO iprofdat_415
   ENDDO ifil_415

!  ---data profile number must match a number defined at an open boundary
   ivar_416: DO ivar=1,novars
      noprofsivar = MAX(MAXVAL(iprofobux(:,ivar)),MAXVAL(iprofobvy(:,ivar)))
      IF (noprofsivar.GT.0) THEN
         ifil_4161: DO ifil=2,nofiles
         iprofdat_4161: DO iprofdat=1,noprofsd(ifil)
            IF (indexvar(iprofdat,ifil).EQ.ivar) THEN
               CALL error_limits_arr(indexprof(iprofdat,ifil),'indexprof',1,&
                                   & noprofsivar,2,indx=(/iprofdat,ifil/))
            ENDIF
         ENDDO iprofdat_4161
         ENDDO ifil_4161
      ENDIF
   ENDDO ivar_416

!  ---all profiles must be defined in only one data file
   ivar_417: DO ivar=1,novars
      noprofsivar = MAX(MAXVAL(iprofobux(:,ivar)),MAXVAL(iprofobvy(:,ivar)))
      iprof_4171: DO iprof=1,noprofsivar
         ncount = 0
         ifil_41711: DO ifil=2,nofiles
         iprofdat_41711: DO iprofdat=1,noprofsd(ifil)
            IF (indexprof(iprofdat,ifil).EQ.iprof.AND.&
              & indexvar(iprofdat,ifil).EQ.ivar) THEN
               ncount = ncount + 1
            ENDIF
         ENDDO iprofdat_41711
         ENDDO ifil_41711

         IF (ncount.EQ.0) THEN
            nerrs = nerrs + 1
            IF (errchk.AND.nerrs.LE.maxerrors) THEN
               WRITE (cprof,'(I12)') iprof; cprof = ADJUSTL(cprof)
               WRITE (cvar,'(I12)') ivar; cvar = ADJUSTL(cvar)
               WRITE (ioerr,'(A)') 'Profile number '//TRIM(cprof)//&
                                 & ' for variable '//TRIM(cvar)//&
                                 & ' not included in any data file'
            ENDIF
         ELSEIF (ncount.GT.1) THEN
            nerrs = nerrs + 1
            IF (errchk.AND.nerrs.LE.maxerrors) THEN
               WRITE (cprof,'(I12)') iprof; cprof = ADJUSTL(cprof)
               WRITE (cvar,'(I12)') ivar; cvar = ADJUSTL(cvar)
               WRITE (ioerr,'(A)') 'Profile number '//TRIM(cprof)//&
                                 & ' for variable '//TRIM(cvar)//&
                               & ' included more than once in the data file(s)' 
            ENDIF
         ENDIF
         
      ENDDO iprof_4171
   ENDDO ivar_417

ENDIF

CALL error_abort('define_profobc_spec',ierrno_input)

!
!5. Write
!--------

IF (modfiles(iddesc,1,2)%defined) THEN
   CALL write_profobc_spec(iddesc,itypobux,itypobvy,iprofobux,iprofobvy,&
                         & iprofrlx,noprofsd,indexprof,indexvar,novars,nofiles,&
                         & nobux,nobvy)
ENDIF

CALL log_timer_out(npcc,itm_bconds)


RETURN

END SUBROUTINE define_profobc_spec

!========================================================================

SUBROUTINE read_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
!************************************************************************
!
! *read_profobc_data* Read open boundary profiles in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_Prof.f90  V2.1.2
!
! Description -
!
! Reference -
!
! Calling program - define_profobc_data
!
! Module calls - error_abort, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_time, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE paralpars
USE syspars
USE error_routines, ONLY: nerrs, error_abort, error_alloc_struc
USE inout_routines, ONLY: open_filepars, read_glbatts_mod, read_time, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*  Arguments
!
CHARACTER (LEN=lentime), INTENT(OUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, numprofs
REAL, INTENT(OUT), DIMENSION(numprofs,nz) :: psiprofdat

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Input date/time
!*psiprofdat*  REAL     Profile arrays
!*numprofs*    INTEGER  Number of profiles in data file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12), DIMENSION(2) :: cval
INTEGER :: l
TYPE (FileParams) :: filepars
INTEGER, DIMENSION(2) :: shape
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_profobc_data'
CALL log_timer_in()

filepars = modfiles(iddesc,ifil,1)

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

!  ---allocate
   ALLOCATE (varatts(2),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/2/),'VariableAtts')

!  ---variable atrributes
   CALL read_varatts_mod(filepars,varatts,2)
   IF (varatts(2)%global_dims(1).NE.numprofs.OR.&
     & varatts(2)%global_dims(2).NE.nz) THEN
      nerrs = 1
      IF (master.AND.errchk) THEN
         WRITE (ioerr,'(A)') 'File: '//TRIM(filepars%filename)
         l_111: DO l=1,2
            WRITE (cval(l),'(I12)') varatts(2)%global_dims(l)
            cval(l) = ADJUSTL(cval(l))
         ENDDO l_111
         WRITE (ioerr,'(A)') 'Wrong shape for variable '&
                           & //TRIM(varatts(2)%f90_name)//': '&
                           & //TRIM(cval(1))//','//TRIM(cval(2))
         shape = (/numprofs,nz/)
         l_112: DO l=1,2
            WRITE (cval(l),'(I12)') shape(l); cval(l) = ADJUSTL(cval(l))
         ENDDO l_112
         WRITE (ioerr,'(A)') 'Must be equal to: '//TRIM(cval(1))//','&
                                               & //TRIM(cval(2))
      ENDIF
      CALL error_abort('read_profobc_data',ierrno_input)
   ENDIF

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

!---profiles
IF (filepars%iostat.LT.3) CALL read_vars(psiprofdat,filepars,2)

1000 modfiles(iddesc,ifil,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_profobc_data

!========================================================================

SUBROUTINE read_profobc_spec(iddesc,itypobux,itypobvy,iprofobux,iprofobvy,&
                           & iprofrlx,noprofsd,indexprof,indexvar,novars,&
                           & nofiles,nobux,nobvy)
!************************************************************************
!
! *read_profobc_spec* Read specifier arrays used for profile open boundary
!                     conditions in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_Prof.f90  V2.10.1
!
! Description -
!
! Reference -
!
! Calling program - define_profobc_spec
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE relaxation
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, nobux, nobvy, nofiles, novars
INTEGER, INTENT(INOUT), DIMENSION(2:nofiles) :: noprofsd
INTEGER, INTENT(OUT), DIMENSION(nobux) :: itypobux
INTEGER, INTENT(OUT), DIMENSION(nobvy) :: itypobvy
INTEGER, INTENT(INOUT), DIMENSION(nobux,novars) :: iprofobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy,novars) :: iprofobvy
INTEGER, INTENT(INOUT), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexprof
INTEGER, INTENT(INOUT), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexvar
INTEGER, INTENT(INOUT), DIMENSION(norlxzones) :: iprofrlx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id
!*itypobux*  INTEGER Type of U- or X-open boundary condition
!*itypobvy*  INTEGER Type of V- or Y-open boundary condition
!*iprofobux* INTEGER Profile numbers at U- or X-open boundaries
!*iprofobvy* INTEGER Profile numbers at V- or Y-open boundaries
!*iprofrlx*  INTEGER Disables/enables relaxation at open boundary zones
!*noprofsd*  INTEGER Number of profiles per data file
!*indexprof* INTEGER Mapping array of the profile numbers in the data files to
!                    the profile numbers assigned to the open boundaries. The
!                    physical size of the first dimension equals the number of
!                    profiles in a data file.
!*indexvar*  INTEGER Defines the variable number of the profiles in a data
!                    file. The physical size of the first dimension equals the
!                    number of profiles in a data file.
!*novars*    INTEGER Total number of variables
!*nofiles*   INTEGER Number of data files
!*nobux*     INTEGER Number of U- or X-open boundary nodes
!*nobvy*     INTEGER Number of V- or Y-open boundary nodes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ifil, numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_profobc_spec'
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
!---U- or X-nodes
varid = 0
IF (nobux.GT.0) THEN 
   varid = varid + 1
   CALL read_vars(itypobux,filepars,varid,varatts)
   varid = varid + 1
   CALL read_vars(iprofobux,filepars,varid,varatts)
ENDIF

!---V- or Y-nodes
IF (nobvy.GT.0) THEN 
   varid = varid + 1
   CALL read_vars(itypobvy,filepars,varid,varatts)
   varid = varid + 1
   CALL read_vars(iprofobvy,filepars,varid,varatts)
ENDIF

!---number of profiles
IF (nofiles.GT.1) THEN
   varid = varid + 1
   CALL read_vars(noprofsd,filepars,varid,varatts)
ENDIF

!---profile and variable number mapping arrays
ifil_211: DO ifil=2,nofiles
   varid = varid + 1
   CALL read_vars(indexprof(1:noprofsd(ifil),ifil),filepars,varid,&
                & varatts=varatts)
   varid = varid + 1
   CALL read_vars(indexvar(1:noprofsd(ifil),ifil),filepars,varid,&
                & varatts=varatts)
ENDDO ifil_211

!---relaxation zones
IF (norlxzones.GT.0) THEN
   varid = varid + 1
   CALL read_vars(iprofrlx,filepars,varid,varatts=varatts)
ENDIF

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(iddesc,1,1) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE read_profobc_spec

!========================================================================

SUBROUTINE update_profobc_data(profdata,obcdata,noprofsd,indexprof,indexvar,&
                             & maxprofs,noprofs,novars,nofiles,nosecsdat,iddesc)
!************************************************************************
!
! *update_profobc_data* Update profiles at open boundaries
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_Prof.f90  V2.7.1
!
! Description - 
!
! Reference -
!
! Calling program - current_corr, salinity_equation, temperature_equation
!
! External calls - define_profobc_data
!
! Module calls - diff_dates, error_file, lim_dims, loop_index
!
!************************************************************************
!
USE gridpars
USE iopars
USE syspars
USE timepars
USE error_routines, ONLY: nerrs, error_file
USE time_routines, ONLY: diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index

IMPLICIT NONE

!
!*  Arguments
!
INTEGER, INTENT(IN) :: iddesc, maxprofs, nofiles, noprofs, novars
INTEGER, INTENT(IN) , DIMENSION(2:nofiles) :: noprofsd
INTEGER, INTENT(IN), DIMENSION(maxprofs,2:nofiles) :: indexprof
INTEGER, INTENT(IN), DIMENSION(maxprofs,2:nofiles) :: indexvar
INTEGER (KIND=kndilong), INTENT(INOUT), DIMENSION(2:nofiles,2) :: nosecsdat
REAL, INTENT(INOUT), DIMENSION(maxprofs,nz,2:nofiles,2) :: profdata
REAL, INTENT(INOUT), DIMENSION(0:noprofs,nz,novars) :: obcdata

!
! Name        Type    Purpose
!------------------------------------------------------------------------------
!*profdata*   REAL    Profile arrays from data file
!          1st dimension => profile number 
!          2nd dimension => vertical level
!          3rd dimension => file number
!          4th dimension => time level (old,new)
!*obcdata*    REAL    Profile arrays at current time
!*noprofsd*   INTEGER Number of profiles per data file
!*indexprof*  INTEGER Mapping array of the profile numbers in the data files to
!                     the profile numbers assigned to the open boundaries. The
!                     physical size of the first dimension equals the number of
!                     profiles in a data file.
!*indexvar*   INTEGER Defines the variable number of the profiles in a data
!                     file. The physical size of the first dimension equals the
!                     number of profiles in a data file.
!*maxprofs*   INTEGER Maximum number of profiles in any data file
!*noprofs*    INTEGER Total number of profiles
!*novars*     INTEGER Total number of variables
!*nofiles*    INTEGER Number of data files (plus 1)
!*nosecsdat*  LONGINT Number of seconds since start of simulation and old/new
!                     data time
!*iddesc*     INTEGER Data file id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: no_time_series
CHARACTER (LEN=12) :: cfil
CHARACTER (LEN=lentime) :: ciodatetime
INTEGER :: ifil, iprof, ivar, l, npcc, numprofs
REAL (KIND=kndrlong) :: ratio1, ratio2
INTEGER, DIMENSION(3) :: tlims
INTEGER, DIMENSION(MaxIOFiles) :: noread


!
!1. Update data
!--------------
!

ifil_100: DO ifil=2,nofiles

   IF (modfiles(iddesc,ifil,1)%defined) THEN

      tlims = modfiles(iddesc,ifil,1)%tlims
      no_time_series = lim_dims(ABS(tlims)).EQ.1
      IF ((nt.GT.tlims(1).AND.no_time_series).OR.&
        & (.NOT.loop_index(tlims,nt)).OR.&
         & modfiles(iddesc,ifil,1)%iostat.EQ.3) CYCLE ifil_100
      WRITE (cfil,'(I12)') ifil; cfil = ADJUSTL(cfil)
      procname(pglev+1) = 'update_profobc_data'
      CALL log_timer_in(npcc=npcc,logname=TRIM(procname(pglev+1))//': '//cfil)

      numprofs = noprofsd(ifil)

!
!1.1 Time-independent data or regular data without time interpolation
!--------------------------------------------------------------------
!

      IF (no_time_series.OR.&
       & (modfiles(iddesc,ifil,1)%time_regular.AND.tlims(3).LT.0.AND.&
       &  nt.LT.nstep)) THEN

         IF (nt.EQ.tlims(1)) THEN
            ciodatetime = cdatetime_undef
            DO WHILE (ciodatetime.NE.CDateTime)
               CALL define_profobc_data(iddesc,ifil,ciodatetime,&
                                      & profdata(1:numprofs,:,ifil,1),numprofs)
               IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
                  CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
               ENDIF
            ENDDO
         ELSE
            CALL define_profobc_data(iddesc,ifil,ciodatetime,&
                                   & profdata(1:numprofs,:,ifil,1),numprofs)
            IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) CYCLE ifil_100
            IF (ciodatetime.NE.CDateTime) GOTO 1001
         ENDIF
      
!
!1.2 Irregular data or with time interpolation
!---------------------------------------------
!

      ELSE

!
!1.2.1 Read data on first call
!-----------------------------
!

         IF (nt.EQ.tlims(1)) THEN
            CALL define_profobc_data(iddesc,ifil,ciodatetime,&
                                   & profdata(1:numprofs,:,ifil,2),numprofs)
            IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
               CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
            ENDIF
            CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat(ifil,2))
            IF (nosecsdat(ifil,2).GT.nosecsrun) GOTO 1002
            noread(ifil) = 1
         ENDIF

!
!1.2.2 Update data
!-----------------
!

         IF (nt.LT.nstep) THEN
            DO WHILE (nosecsdat(ifil,2).LE.nosecsrun)
!              ---store old data
               profdata(1:numprofs,:,ifil,1) = profdata(1:numprofs,:,ifil,2)
               nosecsdat(ifil,1) = nosecsdat(ifil,2)
!              ---read new data
               CALL define_profobc_data(iddesc,ifil,ciodatetime,&
                                      & profdata(1:numprofs,:,ifil,2),numprofs)
               noread(ifil) = noread(ifil) + 1
               IF (noread(ifil).EQ.2.AND.&
                 & modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
                  CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
               ENDIF
               CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat(ifil,2))
            ENDDO
         ENDIF

      ENDIF

!
!1.3 Evaluate data at current time
!---------------------------------
!
!     ---without time interpolation
      IF (no_time_series.OR.tlims(3).LT.0) THEN
         l_131: DO l=1,numprofs
            iprof = indexprof(l,ifil)
            ivar = indexvar(l,ifil)
            obcdata(iprof,:,ivar) = MERGE(profdata(l,:,ifil,1),real_fill,&
                         & ABS(profdata(l,:,ifil,1)-real_fill).GT.real_fill_eps)
         ENDDO l_131

!     ---with time interpolation
      ELSE
         ratio2 = MERGE(0.0,(nosecsrun-nosecsdat(ifil,1))&
                    & /REAL(nosecsdat(ifil,2)-nosecsdat(ifil,1)),tlims(3).LT.0)
         ratio1 = 1.0 - ratio2
         l_132: DO l=1,numprofs
            iprof = indexprof(l,ifil)
            ivar = indexvar(l,ifil)
            WHERE ((ABS(profdata(l,:,ifil,1)-real_fill).GT.real_fill_eps).AND.&
                 & (ABS(profdata(l,:,ifil,2)-real_fill).GT.real_fill_eps))
               obcdata(iprof,:,ivar) = ratio1*profdata(l,:,ifil,1) + &
                                     & ratio2*profdata(l,:,ifil,2)
            ELSEWHERE
               obcdata(iprof,:,ivar) = real_fill
            END WHERE
         ENDDO l_132
      ENDIF

      CALL log_timer_out(npcc,itm_bconds)

   ENDIF

ENDDO ifil_100


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

END SUBROUTINE update_profobc_data

!========================================================================

SUBROUTINE write_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
!************************************************************************
!
! *write_profobc_data* Write open boundary profiles in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_Prof.f90  V2.X
!
! Description -
!
! Reference -
!
! Calling program - define_profobc_data
!
! Module calls - error_alloc_struc, later, noearlier, open_filepars,
!                output_flag, set_modfiles_atts, set_modvars_atts,
!                write_atts_mod, write_time, write_vars
!
!************************************************************************
!
USE datatypes
USE gridpars
USE iopars
USE paralpars
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: open_filepars, output_flag, write_atts_mod, &
                        & write_time, write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines

IMPLICIT NONE

!
!*  Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, numprofs
REAL, INTENT(IN), DIMENSION(numprofs,nz) :: psiprofdat

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*psiprofdat*  REAL     Profile arrays
!*numprofs*    INTEGER  Number of profiles in data file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_profobc_data'
CALL log_timer_in()

!
!1. Disable writing if needed
!----------------------------
!

filepars = modfiles(iddesc,ifil,2)
flag = output_flag(ciodatetime,filepars%tskips)

!
!2. Write file header on first call
!----------------------------------
!

IF (filepars%iostat.EQ.0) THEN
   
!  ---file attributes
   CALL set_modfiles_atts(iddesc,ifil,2)
   filepars = modfiles(iddesc,ifil,2)

!  ---open file
   CALL open_filepars(filepars)

!  ---allocate
   ALLOCATE (varatts(2),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/2/),'VariableAtts')

!  ---variable attributes
   CALL set_modvars_atts(iddesc,ifil,2,varatts,2,numprofs=numprofs)

!  ---write attributes
   CALL write_atts_mod(filepars,varatts,2)

!  ---deallocate
   DEALLOCATE (varatts)

ENDIF

!
!3. Write data
!-------------
!
!---date/time
CALL write_time(ciodatetime,filepars)

!---profiles
CALL write_vars(psiprofdat,filepars,2)

modfiles(iddesc,ifil,2) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE write_profobc_data

!========================================================================

SUBROUTINE write_profobc_spec(iddesc,itypobux,itypobvy,iprofobux,iprofobvy,&
                            & iprofrlx,noprofsd,indexprof,indexvar,novars,&
                            & nofiles,nobux,nobvy)
!************************************************************************
!
! *write_profobc_spec* Write specifier arrays used for profile open boundary
!                      conditions in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Open_Boundary_Data_Prof.f90  V2.1.2
!
! Description -
!
! Reference -
!
! Calling program - define_profobc_spec
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_vars
!
!************************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE relaxation
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, nobux, nobvy, nofiles, novars
INTEGER, INTENT(IN), DIMENSION(2:nofiles) :: noprofsd
INTEGER, INTENT(OUT), DIMENSION(nobux) :: itypobux
INTEGER, INTENT(OUT), DIMENSION(nobvy) :: itypobvy
INTEGER, INTENT(INOUT), DIMENSION(nobux,novars) :: iprofobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy,novars) :: iprofobvy
INTEGER, INTENT(IN), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexprof
INTEGER, INTENT(IN), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexvar
INTEGER, INTENT(IN), DIMENSION(norlxzones) :: iprofrlx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id
!*itypobux*  INTEGER Type of U- or X-open boundary condition
!*itypobvy*  INTEGER Type of V- or Y-open boundary condition
!*iprofobux* INTEGER Profile numbers at U- or X-open boundaries
!*iprofobvy* INTEGER Profile numbers at V- or Y-open boundaries
!*iprofrlx*  INTEGER Disables/enables relaxation at 4 boundary zones
!*noprofsd*  INTEGER Number of profiles per data file
!*indexprof*  INTEGER Mapping array of the profile numbers in the data files to
!                    the profile numbers assigned to the open boundaries. The
!                    physical size of the first dimension equals the number of
!                    profiles in a data file.
!*indexvar*   INTEGER Defines the variable number of the profiles in a data
!                    file. The physical size of the first dimension equals the
!                    number of profiles in a data file.
!*novars*    INTEGER Total number of variables
!*nofiles*   INTEGER Number of data files
!*nobux*     INTEGER Number of nodes at U- or X-open boundaries
!*nobvy*     INTEGER Number of nodes at V- or Y-open boundaries
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ifil, numvars, varid
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_profobc_spec'
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
CALL set_modvars_atts(iddesc,1,2,varatts,numvars,noprofsd=noprofsd,&
                    & novars=novars)

!---write
CALL write_atts_mod(filepars,varatts,numvars)

!
!2. Write data
!-------------
!
!---U- or X-nodes
varid = 0
IF (nobux.GT.0) THEN
   varid = varid + 1
   CALL write_vars(itypobux,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(iprofobux,filepars,varid,varatts=varatts)
ENDIF

!---V- or Y-nodes
IF (nobvy.GT.0) THEN
   varid = varid + 1
   CALL write_vars(itypobvy,filepars,varid,varatts=varatts)
   varid = varid + 1
   CALL write_vars(iprofobvy,filepars,varid,varatts=varatts)
ENDIF

!---number of profiles
IF (nofiles.GT.1) THEN
   varid = varid + 1
   CALL write_vars(noprofsd,filepars,varid,varatts=varatts)
ENDIF

!---profile and variable indices
ifil_211: DO ifil=2,nofiles
   varid = varid + 1
   CALL write_vars(indexprof(1:noprofsd(ifil),ifil),filepars,varid,&
                 & varatts=varatts)
   varid = varid + 1
   CALL write_vars(indexvar(1:noprofsd(ifil),ifil),filepars,varid,&
                          & varatts=varatts)
ENDDO ifil_211

!---relaxation zones
IF (norlxzones.GT.0) THEN
   varid = varid + 1
   CALL write_vars(iprofrlx,filepars,varid,varatts=varatts)
ENDIF

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(iddesc,1,2) = filepars
DEALLOCATE (varatts)

CALL log_timer_out()


RETURN

END SUBROUTINE write_profobc_spec
