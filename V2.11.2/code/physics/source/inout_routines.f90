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

MODULE inout_routines
!************************************************************************
!
! *inout_routines* Routines for standard I/O operations
!
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description -
!
! Generic routines - read_submod, read_time, read_vars, write_time, write_vars
!
! Routines - close_file, close_filepars, inout_atts_out, inquire_modfil_atts,
!            inquire_modfil_dims, inquire_modfil_vars, monitor_files,
!            open_file, open_filepars, read_glbatts_mod, read_metadata_line,
!            read_output_metadata, read_station_names, read_varatts_mod,
!            write_atts_mod, write_metadata_line, write_info_mod
! Functions - get_unit, inquire_unit, output_flag
!
!************************************************************************
!
USE iopars
USE syspars
USE time_routines

IMPLICIT NONE

INTERFACE read_submod
   MODULE PROCEDURE read_submod_real_2d, read_submod_real_3d, &
                  & read_submod_real_4d
END INTERFACE

INTERFACE read_time
   MODULE PROCEDURE read_time_char, read_time_real
END INTERFACE

INTERFACE read_vars
   MODULE PROCEDURE read_vars_char_0d, read_vars_char_1d, read_vars_char_2d, &
                  & read_vars_int_0d, read_vars_int_1d, read_vars_int_2d, &
                  & read_vars_int_3d, read_vars_int_4d, read_vars_log_0d, &
                  & read_vars_log_1d, read_vars_log_2d, read_vars_real_0d, &
                  & read_vars_real_1d, read_vars_real_2d, read_vars_real_3d, &
                  & read_vars_real_4d, read_vars_double_0d, &
                  & read_vars_double_1d, read_vars_double_2d, &
                  & read_vars_double_3d, read_vars_double_4d
END INTERFACE

INTERFACE write_submod
   MODULE PROCEDURE write_submod_real_2d, write_submod_real_3d, &
                  & write_submod_real_4d
END INTERFACE

INTERFACE write_time
   MODULE PROCEDURE write_time_char, write_time_real
END INTERFACE

INTERFACE write_vars
   MODULE PROCEDURE write_vars_char_0d, write_vars_char_1d, write_vars_char_2d,&
                  & write_vars_int_0d, write_vars_int_1d, write_vars_int_2d, &
                  & write_vars_int_3d, write_vars_int_4d, write_vars_log_0d, &
                  & write_vars_log_1d, write_vars_log_2d, write_vars_real_0d, &
                  & write_vars_real_1d, write_vars_real_2d, &
                  & write_vars_real_3d, write_vars_real_4d, &
                  & write_vars_double_0d, write_vars_double_1d, &
                  & write_vars_double_2d, write_vars_double_3d, &
                  & write_vars_double_4d
END INTERFACE

CONTAINS


!========================================================================

SUBROUTINE close_file(iounit,ioform,filename,fildel)
!************************************************************************
!
! *close_file* Close a file connection
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - file is deleted if fildel is present and .TRUE.
!             - file deletion is  
!
! Module calls - cf90_close, error_abort, error_file, inquire_file, inquire_unit
!
!************************************************************************
!
USE cf90_routines, ONLY: cf90_close
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: fildel
CHARACTER (LEN=lenform), INTENT(IN) :: ioform
CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: filename
INTEGER, INTENT(INOUT) :: iounit

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iounit*    INTEGER File unit
!*ioform*    CHAR    File format ('A', 'U', or 'N')
!*filename*  CHAR    File name
!*fildel*    LOGICAL File is deleted if PRESENT and .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: file_delete, flag
CHARACTER (LEN=12) :: cunit


pglev = pglev + 1
procname(pglev) = 'close_file'

!
!1. Write info to log file
!-------------------------
!

WRITE (cunit,'(I12)') iounit; cunit = ADJUSTL(cunit)
IF (loglev1.GT.0) THEN
   IF (PRESENT(filename)) THEN
      WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Close file '//TRIM(filename)//&
                        & ' on unit '//TRIM(cunit)//' ('//ioform//')'
   ELSE
      WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Close file on unit '//&
                        & TRIM(cunit)//' ('//ioform//')'
   ENDIF
ENDIF

!
!2. Optional argument
!--------------------
!

IF (PRESENT(fildel)) THEN
   file_delete = fildel
ELSE
   file_delete = .FALSE.
ENDIF

!
!3. Check whether file unit is opened
!------------------------------------
!

flag = inquire_unit(iounit,ioform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,iounit=iounit)

!
!4. Close file
!-------------
!

IF (ioform.NE.'N') THEN
   IF (file_delete) THEN
      CLOSE (iounit,STATUS='DELETE',ERR=1000)
   ELSE
      CLOSE (iounit,ERR=1000)
   ENDIF
ELSE
   CALL cf90_close(iounit)
   IF (errstat.NE.noerr_NF90) GOTO 1000
ENDIF

!
!5. Update number of opened files
!--------------------------------
!

nopenf = nopenf - 1

!
!6. Reset unit number
!--------------------
!

iounit = int_fill

pglev = pglev - 1


RETURN

!
!6. Process close error
!----------------------
!

1000 nerrs = 1
IF (errchk) THEN
   IF (PRESENT(filename)) THEN
      WRITE (ioerr,'(A)') 'Unable to close file '//TRIM(filename)//&
                        & ' on unit '//cunit
   ELSE
      WRITE (ioerr,'(A)') 'Unable to close file on unit '//cunit
   ENDIF
ENDIF
CALL error_abort(procname(pglev),ierrno_fclose)

END SUBROUTINE close_file

!========================================================================

SUBROUTINE close_filepars(filepars,fildel)
!************************************************************************
!
! *close_filepars* Close a file connection
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - file is deleted if fildel is present and .TRUE.
!             - file deletion is not allowed for netCDF files
!             - uses file atrributes stored in 'filepars'
!
! Module calls - cf90_close, cf90_enddef, cf90_put_att, cf90_redef,
!                error_abort, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_close, cf90_enddef, cf90_put_att, cf90_redef
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: fildel
TYPE(FileParams), INTENT(INOUT) :: filepars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED File attributes
!*fildel*    LOGICAL File is deleted if PRESENT and .TRUE.
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: file_delete, flag
CHARACTER (LEN=lenform) :: ioform
CHARACTER (LEN=12) :: cunit
CHARACTER (LEN=leniofile) :: filename 
INTEGER :: iounit


procname(pglev+1) = 'close_filepars'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

filename = filepars%filename
ioform = filepars%form
iounit = filepars%iunit

!
!2. Write info to log file
!-------------------------
!

WRITE (cunit,'(I12)') iounit; cunit = ADJUSTL(cunit)
IF (loglev1.GT.0) THEN
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Close file '//TRIM(filename)//&
                     & ' on unit '//TRIM(cunit)//' ('//ioform//')'
ENDIF

!
!3. Optional argument
!--------------------
!

IF (PRESENT(fildel)) THEN
   file_delete = fildel
ELSE
   file_delete = .FALSE.
ENDIF

!
!4. Check whether file unit is opened
!------------------------------------
!

flag = inquire_unit(iounit,ioform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!5. Close file
!-------------
!

IF (ioform.NE.'N') THEN
   IF (file_delete) THEN
      CLOSE (iounit,STATUS='DELETE',ERR=1000)
   ELSE
      CLOSE (iounit,ERR=1000)
   ENDIF
ELSE
   CALL cf90_close(iounit)
   IF (errstat.NE.noerr_NF90) GOTO 1000
ENDIF

!
!6. Update/reset parameters
!--------------------------
!

nopenf = nopenf - 1
filepars%iunit = int_fill
filepars%iostat = 0
filepars%timerec = 0

CALL log_timer_out()


RETURN

!
!7. Process close error
!----------------------
!

1000 nerrs = 1
IF (errchk) THEN
   WRITE (ioerr,'(A)') 'Unable to close file '//TRIM(filepars%filename)//&
                     & ' on unit '//cunit
ENDIF
CALL error_abort(procname(pglev),ierrno_fclose)

END SUBROUTINE close_filepars

!========================================================================

FUNCTION get_unit()
!************************************************************************
!
! *get_unit* Return next available unit number for file connection
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.0
!
! Description -
!
! Module calls - error_abort
!
!************************************************************************
!
USE error_routines, ONLY: error_abort

!
!*Arguments
!
INTEGER :: get_unit

!
!*Local variables
!
LOGICAL :: flag, fopen
CHARACTER (LEN=12) :: cunext
INTEGER :: iunext

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*flag*    LOGICAL  Error flag (.TRUE. if unit is not wothin allowed range)
!*fopen*   LOGICAL  .TRUE. if unit number is already connected to a file
!*iunext*  INTEGER  Next available unit number
!
!------------------------------------------------------------------------------
!

iunext = 1000
fopen =.TRUE.

DO WHILE (fopen)
   
   INQUIRE (UNIT=iunext,EXIST=flag,OPENED=fopen)
   IF (.NOT.flag) THEN
      nerrs = 1
      IF (errchk) THEN
         WRITE (cunext,'(I12)') iunext; cunext = ADJUSTL(cunext)
         WRITE (ioerr,'(A)') 'Invalid value for logical unit number: '&
                           & //TRIM(cunext)
         WRITE (ioerr,'(A)') 'Value is not in the range allowed by the&
                            & processor'
      ENDIF
      CALL error_abort('get_unit',ierrno_fopen)
   ENDIF
   
   IF (fopen.OR.iunext.EQ.5.OR.iunext.EQ.6) THEN
      iunext = iunext + 1
      fopen  = .TRUE.
   ELSE
      get_unit = iunext
   ENDIF
   
ENDDO


RETURN

END FUNCTION get_unit

!========================================================================

SUBROUTINE inout_atts_out(file_type)
!************************************************************************
!
! *inout_atts_out* Write the global, coordinate and variable attributes of a
!                  user output file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - write info if required
!
! Module calls - cf90_def_dim, cf90_def_var, cf90_enddef, cf90_put_att,
!                cf90_put_att_chars, cf90_put_var_chars, cf90_put_var,
!                close_file, close_filepars, conv_to_chars,
!                conv_to_chars_modvars, convert_date, error_alloc,
!                error_alloc_struc, error_file, inquire_var, open_file,
!                open_filepars, varatts_init
!
!************************************************************************
!
USE datatypes
USE modids
USE paralpars
USE physpars
USE switches
USE timepars
USE cf90_routines, ONLY: cf90_def_dim, cf90_def_var, cf90_enddef, cf90_put_att,&
                      &  cf90_put_att_chars, cf90_put_var_chars, cf90_put_var
USE cif_routines, ONLY: conv_to_chars, conv_to_chars_modvars
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_alloc, error_alloc_struc, error_file
USE modvars_routines, ONLY: file_suffix, inquire_var

!
!*Arguments
!
CHARACTER (LEN=2), INTENT(IN) :: file_type

!
! Name       Type     Purpose
!------------------------------------------------------------------------------
!*file_type* CHAR     Type of output file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: info, gridded, packing, time_grid
CHARACTER (LEN=1) :: cdim, floattype
CHARACTER (LEN=3) :: suffix, xname, yname
CHARACTER (LEN=12) :: cfreq, cset
CHARACTER (LEN=lenform) :: dataform, outform
CHARACTER (LEN=lendesc) :: comment, coordinates, formula_terms
CHARACTER (LEN=lenunit) :: tunits
CHARACTER (LEN=leniofile) :: infofile
CHARACTER (LEN=lencifline) :: cline
CHARACTER (LEN=lendesc), DIMENSION(15) :: cvals
CHARACTER (LEN=lenname), SAVE, ALLOCATABLE, DIMENSION(:) :: dimnames
INTEGER :: b_id, dimname, dimspt, dimstat, dimt, dimtlen, dimx, dimy, dimz, &
         & hcrit_id, ifreq, iifreq, iivar, ifil, iglb, iodat, ioinfo, iset, &
         & istat, iunit, ivar, l, ncout, nfvars, nocoords, nodim, nodims, &
         & noglbatts, nosets, nostats, novars, nowetout, npcc, nrank, nrout, &
         & nstepout, numfreqs, numstats, nvals, nzout, statid, theta_id, &
         & time_format, varid, vcoord
REAL :: deltout
REAL (KIND=kndrlong) :: fill_value
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: outdims
INTEGER, ALLOCATABLE, DIMENSION(:) :: indstats, indvars
TYPE (FileParams) :: filepars
TYPE (OutGridParams) :: outgpars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: coordatts, filevars, varatts
TYPE (StationLocs), ALLOCATABLE, DIMENSION(:) :: outlocs, outstatlocs


procname(pglev+1) = 'inout_atts_out'
CALL log_timer_in(npcc)

!
!1. Initialise parameters and arrays
!-----------------------------------
!
!1.1 Parameters
!--------------
!

SELECT CASE (file_type)
   CASE ('TS')
      nosets = nosetstsr
      nfvars = novarstsr
      numstats = nostatstsr
      numfreqs = 1
   CASE ('TA')
      nosets = nosetsavr
      nfvars = novarsavr
      numstats = nostatsavr
      numfreqs = 1
   CASE ('HR')
      nosets = nosetsanal
      nfvars = novarsanal
      numstats = nostatsanal
      numfreqs = 1
   CASE ('HA','HP')
      nosets = nosetsanal
      nfvars = novarsanal
      numstats = nostatsanal
      numfreqs = nofreqsanal
   CASE ('HE')
      nosets = nosetsanal
      nfvars = 14
      numstats = nostatsanal
      numfreqs = nofreqsanal
END SELECT

!
!1.2 Allocate
!------------
!

ALLOCATE (indvars(nfvars),STAT=errstat)
CALL error_alloc('indvars',1,(/nfvars/),kndint)
ALLOCATE (indstats(numstats),STAT=errstat)
CALL error_alloc('indstats',1,(/numstats/),kndint)
ALLOCATE (outstatlocs(numstats),STAT=errstat)
CALL error_alloc_struc('outstatlocs',1,(/numstats/),'StationLocs')
ALLOCATE (outlocs(numstats),STAT=errstat)
CALL error_alloc_struc('outlocs',1,(/numstats/),'StationLocs')
ALLOCATE (varatts(nfvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/nfvars/),'VariableAtts')
CALL varatts_init(varatts)

!
!1.3 Variable and station attributes
!-----------------------------------
!

SELECT CASE (file_type)
   CASE ('TS')
      varatts = tsrvars
   CASE ('TA')
      varatts = avrvars
   CASE ('HR','HA','HP')
      varatts = analvars
   CASE ('HE')
      varatts = ellvars
END SELECT

!
!2. Define and write metadata
!----------------------------
!

iset_200: DO iset=1,nosets
ifreq_200: DO ifreq=1,numfreqs
nodim_200: DO nodim=0,3
   IF (nodim.EQ.1.OR.(nodim.EQ.0.AND.file_type.EQ.'HE')) CYCLE nodim_200

!
!2.1 Define parameters
!---------------------
!

   SELECT CASE (file_type)
      CASE ('TS')
         IF (nodim.EQ.0) filepars = tsr0d(iset)
         IF (nodim.EQ.2) filepars = tsr2d(iset)
         IF (nodim.EQ.3) filepars = tsr3d(iset)
         outgpars = tsrgpars(iset)
         indvars = ivarstsr(iset,:)
         IF (numstats.GT.0) THEN
            indstats = lstatstsr(iset,:)
            outstatlocs = tsrstatlocs
         ENDIF
      CASE ('TA')
         IF (nodim.EQ.0) filepars = avr0d(iset)
         IF (nodim.EQ.2) filepars = avr2d(iset)
         IF (nodim.EQ.3) filepars = avr3d(iset)
         outgpars = avrgpars(iset)
         indvars = ivarsavr(iset,:)
         IF (numstats.GT.0) THEN
            indstats = lstatsavr(iset,:)
            outstatlocs = avrstatlocs
         ENDIF
      CASE ('HR')
         IF (nodim.EQ.0) filepars = res0d(iset)
         IF (nodim.EQ.2) filepars = res2d(iset)
         IF (nodim.EQ.3) filepars = res3d(iset)
         outgpars = analgpars(iset)
         indvars = ivarsanal(iset,:)
         IF (numstats.GT.0) THEN
            indstats = lstatsanal(iset,:)
            outstatlocs = analstatlocs
         ENDIF
      CASE ('HA')
         IF (nodim.EQ.0) filepars = amp0d(iset,ifreq)
         IF (nodim.EQ.2) filepars = amp2d(iset,ifreq)
         IF (nodim.EQ.3) filepars = amp3d(iset,ifreq)
         outgpars = analgpars(iset)
         indvars = ivarsanal(iset,:)
         IF (numstats.GT.0) THEN
            indstats = lstatsanal(iset,:)
            outstatlocs = analstatlocs
         ENDIF
      CASE ('HP')
         IF (nodim.EQ.0) filepars = pha0d(iset,ifreq)
         IF (nodim.EQ.2) filepars = pha2d(iset,ifreq)
         IF (nodim.EQ.3) filepars = pha3d(iset,ifreq)
         outgpars = analgpars(iset)
         indvars = ivarsanal(iset,:)
         IF (numstats.GT.0) THEN
            indstats = lstatsanal(iset,:)
            outstatlocs = analstatlocs
         ENDIF
      CASE ('HE')
         IF (nodim.EQ.2) filepars = ell2d(iset,ifreq)
         IF (nodim.EQ.3) filepars = ell3d(iset,ifreq)
         outgpars = analgpars(iset)
         indvars = ivarsell(iset,:)
         IF (numstats.GT.0) THEN
            indstats = lstatsanal(iset,:)
            outstatlocs = analstatlocs
         ENDIF
   END SELECT

!
!2.2 Return if output is disabled
!--------------------------------
!

   IF (.NOT.(master.AND.filepars%defined)) CYCLE nodim_200

!
!2.3 Output parameters
!---------------------
!
!  ---initialise
   info = filepars%info
   novars = filepars%novars
   outform = filepars%form
   floattype = filepars%floattype
   packing = outgpars%packing
   filepars%packing = outgpars%packing
   gridded = outgpars%gridded
   vcoord = outgpars%vcoord
   time_grid = outgpars%time_grid.AND.nodim.EQ.3
   time_format = outgpars%time_format
   IF (floattype.EQ.'S') THEN
      fill_value = real_fill
   ELSE
      fill_value = double_fill
   ENDIF
   
!  ---dimensions
   ncout = MERGE(0,outgpars%ncout,nodim.EQ.0)
   nrout = MERGE(0,outgpars%nrout,nodim.EQ.0)
   nzout = MERGE(0,outgpars%nzout,nodim.EQ.0.OR.nodim.EQ.2)
   nowetout = MERGE(0,outgpars%nowetout,nodim.EQ.0)
   nostats = MERGE(0,outgpars%nostats,nodim.EQ.0)

!  ---global attributes
   glbatts(2)%value = TRIM(filepars%title)

!  ---station attributes
   IF (.NOT.gridded) THEN
      istat_230: DO istat=1,nostats
         outlocs(istat) = outstatlocs(indstats(istat))
      ENDDO istat_230
   ENDIF
   
!  ---number of time steps
   nstepout = outgpars%nstepout
   
!  ---time step
   deltout = outgpars%deltout

!  ---coordinates attribute
   IF (.NOT.gridded.OR.iopt_grid_htype.EQ.3) THEN
      IF (iopt_grid_sph.EQ.0) THEN
         IF (vcoord.EQ.1) THEN
            coordinates = 'x y z time'
         ELSE
            coordinates = 'x y lev time'
         ENDIF
      ELSE
         IF (vcoord.EQ.1) THEN
            coordinates = 'lon lat z time'
         ELSE
            coordinates = 'lon lat lev time'
         ENDIF
      ENDIF
   ELSE
      coordinates = ''
   ENDIF

!  ---time unit in CF format
   SELECT CASE (time_format)
      CASE (0); tunits = 'Calendar date and time'
      CASE (1); tunits = 'seconds since'
      CASE (2); tunits = 'minutes since'
      CASE (3); tunits = 'hours since'
      CASE (4); tunits = 'days since'
      CASE (5); tunits = 'weeks since'
      CASE (6); tunits = 'months since'
      CASE (7,8); tunits = 'years since'
   END SELECT
   IF (time_format.GT.0) tunits = TRIM(tunits)//' '//outgpars%refdate(1:19)

!
!2.4 File name
!-------------
!

   IF (TRIM(filepars%filename).EQ.'') THEN
      suffix = file_suffix(outform)
      WRITE (cset,'(I12)') iset; cset = ADJUSTL(cset)
      WRITE (cdim,'(I1)') nodim
      WRITE (cfreq,'(I12)') ifreq; cfreq = ADJUSTL(cfreq)
      SELECT CASE (file_type)
         CASE ('TS')
            filepars%filename = TRIM(outtitle)//'_'//TRIM(cset)//'.'//&
                             & 'tsout'//cdim//'.'//TRIM(suffix)
         CASE ('TA')
            filepars%filename = TRIM(outtitle)//'_'//TRIM(cset)//'.'//&
                             & 'avrgd'//cdim//'.'//TRIM(suffix)
         CASE ('HR')
            filepars%filename = TRIM(outtitle)//'_'//TRIM(cset)//'.'//&
                             & 'resid'//cdim//'.'//TRIM(suffix)
         CASE ('HA')
            filepars%filename = TRIM(outtitle)//'_'//TRIM(cset)//'.'//&
                              & TRIM(cfreq)//'amplt'//cdim//'.'//TRIM(suffix)
         CASE ('HP')
            filepars%filename = TRIM(outtitle)//'_'//TRIM(cset)//'.'//&
                              & TRIM(cfreq)//'phase'//cdim//'.'//TRIM(suffix)
         CASE ('HE') 
            filepars%filename = TRIM(outtitle)//'_'//TRIM(cset)//'.'//&
                              & TRIM(cfreq)//'ellip'//cdim//'.'//TRIM(suffix)
      END SELECT
      IF (TRIM(filepars%pathname).NE.'') THEN
         filepars%filename = TRIM(filepars%pathname)//TRIM(filepars%filename)
      ENDIF
   ENDIF
   IF (info) infofile = TRIM(filepars%filename)//'.inf'

!
!2.5 Open files
!--------------
!

   CALL open_filepars(filepars)
   iodat = filepars%iunit
   IF (info) CALL open_file(ioinfo,infofile,'OUT','A')

!
!2.6 Define netcdf dimensions
!----------------------------
!

   IF (outform.EQ.'N') THEN
      IF (time_format.EQ.0) THEN
         CALL cf90_def_dim(iodat,'tlendim',lentime,dimtlen)
      ENDIF
      IF (iopt_CDF_tlim.EQ.1) THEN
         CALL cf90_def_dim(iodat,'time',nstepout,dimt)
      ELSE
         CALL cf90_def_dim(iodat,'time',unlimited_NF90,dimt)
      ENDIF
      IF (nodim.GT.0) THEN
         IF (gridded) THEN
            IF (iopt_grid_htype.LE.2) THEN
               IF (iopt_grid_sph.EQ.0) THEN
                  CALL cf90_def_dim(iodat,'x',ncout,dimx)
                  CALL cf90_def_dim(iodat,'y',nrout,dimy)
               ELSE
                  CALL cf90_def_dim(iodat,'lon',ncout,dimx)
                  CALL cf90_def_dim(iodat,'lat',nrout,dimy)
               ENDIF
            ELSE
               CALL cf90_def_dim(iodat,'nc',ncout,dimx)
               CALL cf90_def_dim(iodat,'nr',nrout,dimy)
            ENDIF
            IF (packing) THEN
               CALL cf90_def_dim(iodat,'seapoint',nowetout,dimspt)       
            ENDIF
         ELSE
            CALL cf90_def_dim(iodat,'station',nostats,dimstat)
            CALL cf90_def_dim(iodat,'lenname',lendesc,dimname)
         ENDIF
         IF (nodim.EQ.3) THEN
            IF (vcoord.EQ.1) THEN
               CALL cf90_def_dim(iodat,'z',nzout,dimz)
            ELSE
               CALL cf90_def_dim(iodat,'lev',nzout,dimz)
            ENDIF
         ENDIF
      ENDIF
   ELSE
      dimt = 0; dimtlen = 0; dimx = 0; dimy = 0; dimspt = 0
      dimstat = 0; dimz = 0; dimname = 0
   ENDIF

!
!2.7 Coordinate variables
!------------------------
!
!2.7.1 Number of coordinates
!---------------------------
!

   nocoords = 1
   IF (nodim.GT.0) THEN
      nocoords = nocoords + 3
      IF (packing) nocoords = nocoords + 1
      IF (nodim.EQ.3) THEN
         IF (vcoord.EQ.1) THEN
            nocoords = nocoords + 1
         ELSEIF (vcoord.EQ.2) THEN
            nocoords = nocoords + 1
            IF (time_grid) nocoords = nocoords + 1
            IF (iopt_grid_vtype_transf.EQ.21) nocoords = nocoords + 3
         ENDIF
      ENDIF
      IF (.NOT.gridded) nocoords = nocoords + 1
   ENDIF

!
!2.7.2 Initialise
!----------------
!

   ALLOCATE (coordatts(nocoords),STAT=errstat)
   CALL error_alloc_struc('coordatts',1,(/nocoords/),'VariableAtts')
   CALL varatts_init(coordatts)

!
!2.7.3 Coordinate attributes
!---------------------------
!
!  ---data type
   IF (outform.EQ.'N') THEN
      coordatts(1)%data_type = MERGE(char_NF90,double_NF90,time_format.EQ.0)
      coordatts(2:nocoords)%data_type = MERGE(real_NF90,double_NF90,&
                                            & floattype.EQ.'S')
   ELSE
      coordatts(1)%data_type = MERGE(char_type,rlong_type,time_format.EQ.0)
      coordatts(2:nocoords)%data_type = MERGE(real_type,rlong_type,&
                                            & floattype.EQ.'S')
   ENDIF

!  ---fill value
   coordatts%fill = .FALSE.
   coordatts%fill_value = fill_value

!  ---number of coordinates
   filepars%nocoords = nocoords

!  ---f90_name, standard_name, long_name, units
   CALL inquire_var(iarr_time,f90_name=coordatts(1)%f90_name,&
                  & standard_name=coordatts(1)%standard_name,&
                  & long_name=coordatts(1)%long_name)
   coordatts(1)%units = TRIM(tunits)
   ivar = 1
   IF (nodim.GT.0) THEN
      CALL inquire_var(iarr_xout,f90_name=coordatts(ivar+1)%f90_name,&
                     & standard_name=coordatts(ivar+1)%standard_name,&
                     & long_name=coordatts(ivar+1)%long_name,&
                     & units=coordatts(ivar+1)%units)
      CALL inquire_var(iarr_yout,f90_name=coordatts(ivar+2)%f90_name,&
                     & standard_name=coordatts(ivar+2)%standard_name,&
                     & long_name=coordatts(ivar+2)%long_name,&
                     & units=coordatts(ivar+2)%units)
      ivar = ivar + 2
      IF (packing) THEN
         coordatts(ivar+1)%data_type = MERGE(int_NF90,int_type,&
                                           & outform.EQ.'N')
         coordatts(ivar+1)%f90_name = 'seapoint'
         coordatts(ivar+1)%standard_name = 'sea_point_indices'
         coordatts(ivar+1)%long_name = 'Sea point grid indices'
         coordatts(ivar+1)%units = '_'
         ivar = ivar + 1
      ENDIF
      CALL inquire_var(iarr_depout,&
                     & f90_name=coordatts(ivar+1)%f90_name,&
                     & standard_name=coordatts(ivar+1)%standard_name,&
                     & long_name=coordatts(ivar+1)%long_name,&
                     & units=coordatts(ivar+1)%units)
      coordatts(ivar+1)%fill = filepars%fill
      ivar = ivar + 1
      IF (nodim.EQ.3.AND.vcoord.EQ.1) THEN
         IF (.NOT.time_grid) THEN
            CALL inquire_var(iarr_zcmean,&
                     & f90_name=coordatts(ivar+1)%f90_name,&
                     & standard_name=coordatts(ivar+1)%standard_name,&
                     & long_name=coordatts(ivar+1)%long_name,&
                     & units=coordatts(ivar+1)%units)
         ELSE
            CALL inquire_var(iarr_zcoord,&
                     & f90_name=coordatts(ivar+1)%f90_name,&
                     & standard_name=coordatts(ivar+1)%standard_name,&
                     & long_name=coordatts(ivar+1)%long_name,&
                     & units=coordatts(ivar+1)%units)
         ENDIF
         ivar = ivar + 1
      ELSEIF (nodim.EQ.3.AND.vcoord.EQ.2) THEN
         CALL inquire_var(iarr_levout,f90_name=coordatts(ivar+1)%f90_name,&
                     & standard_name=coordatts(ivar+1)%standard_name,&
                     & long_name=coordatts(ivar+1)%long_name,&
                     & units=coordatts(ivar+1)%units)
         ivar = ivar + 1
         IF (time_grid) THEN
            CALL inquire_var(iarr_zeta,&
                     & standard_name=coordatts(ivar+1)%standard_name,&
                     & long_name=coordatts(ivar+1)%long_name,&
                     & units=coordatts(ivar+1)%units)
            coordatts(ivar+1)%fill = filepars%fill
            coordatts(ivar+1)%f90_name = 'eta'
            ivar = ivar + 1
         ENDIF
         IF (iopt_grid_vtype_transf.EQ.21) THEN
            coordatts(ivar+1:ivar+3)%f90_name = &
                                    & (/'theta_SH','b_SH    ','hcrit_SH'/)
            coordatts(ivar+1:ivar+3)%standard_name = &
                                    & 'parameters_for_ocean_s-coordinate'
            coordatts(ivar+1:ivar+3)%long_name = &
                                    & 'Parameters for ocean s-coordinate'
            coordatts(ivar+1:ivar+3)%units = (/' ',' ','m'/)
            ivar = ivar + 3
         ENDIF
      ENDIF
      IF (.NOT.gridded) THEN
         coordatts(ivar+1)%f90_name = 'station_name'
         coordatts(ivar+1)%long_name = 'Station name'
         coordatts(ivar+1)%data_type = MERGE(char_NF90,char_type,&
                                           & outform.EQ.'N')
         ivar = ivar + 1
      ENDIF
   ENDIF

!  ---rank, dimids, other atrributes
   IF (time_format.GT.0) THEN
      coordatts(1)%nrank = 1
      coordatts(1)%global_dims(1) = nstepout
      coordatts(1)%dimids(1) = dimt
   ELSE
      coordatts(1)%nrank = 2
      coordatts(1)%global_dims(1:2) = (/lentime,nstepout/)
      coordatts(1)%dimids(1:2) = (/dimtlen,dimt/)
   ENDIF
   ivar = 1
   IF (nodim.GT.0) THEN
!     ---gridded
      IF (gridded) THEN
!        --horizontal coordinates
         IF (iopt_grid_htype.LE.2) THEN
            coordatts(ivar+1:ivar+2)%nrank = 1
            coordatts(ivar+1)%global_dims(1) = ncout
            coordatts(ivar+1)%dimids(1) = dimx
            coordatts(ivar+2)%global_dims(1) = nrout
            coordatts(ivar+2)%dimids(1) = dimy
         ELSE
            coordatts(ivar+1:ivar+2)%nrank = 2
            coordatts(ivar+1)%global_dims(1:2) = (/ncout,nrout/)
            coordatts(ivar+1)%dimids(1:2) = (/dimx,dimy/)
            coordatts(ivar+1)%axis = 'X'
            coordatts(ivar+2)%global_dims(1:2) = (/ncout,nrout/)
            coordatts(ivar+2)%dimids(1:2) = (/dimx,dimy/)
            coordatts(ivar+2)%axis = 'Y'
         ENDIF
         ivar = ivar + 2
!        --land mask + masked bathymetry
         IF (packing) THEN
            coordatts(ivar+1:ivar+2)%nrank = 1
            coordatts(ivar+1:ivar+2)%global_dims(1) = nowetout
            coordatts(ivar+1:ivar+2)%dimids(1) = dimspt
            ivar = ivar + 2
!        --bathymetry without land mask
         ELSE
            coordatts(ivar+1)%nrank = 2
            coordatts(ivar+1)%global_dims(1:2) = (/ncout,nrout/)
            coordatts(ivar+1)%dimids(1:2) = (/dimx,dimy/)
            ivar = ivar + 1
         ENDIF

!     ---non-gridded
      ELSE
!        --horizontal coordinates and bathymetry
         coordatts(ivar+1:ivar+3)%nrank = 1
         coordatts(ivar+1:ivar+3)%global_dims(1) = nostats
         coordatts(ivar+1:ivar+3)%dimids(1) = dimstat
         ivar = ivar + 3
      ENDIF

!     ---vertical grid
      IF (nodim.EQ.3) THEN
         IF (vcoord.EQ.1) THEN
!           --vertical coordinate 
            IF (gridded) THEN
               IF (.NOT.packing) THEN
                  coordatts(ivar+1)%nrank = MERGE(4,3,time_grid)
                  coordatts(ivar+1)%global_dims(1:3) = (/ncout,nrout,nzout/)
                  coordatts(ivar+1)%dimids(1:3) = (/dimx,dimy,dimz/)
               ELSE
                  coordatts(ivar+1)%nrank = MERGE(3,2,time_grid)
                  coordatts(ivar+1)%global_dims(1:2) = (/nowetout,nzout/)
                  coordatts(ivar+1)%dimids(1:2) = (/dimspt,dimz/)
               ENDIF
            ELSE
               coordatts(ivar+1)%nrank = MERGE(3,2,time_grid)
               coordatts(ivar+1)%global_dims(1:2) = (/nostats,nzout/)
               coordatts(ivar+1)%dimids(1:2) = (/dimstat,dimz/)
            ENDIF
            IF (time_grid) THEN
               nrank = coordatts(ivar+1)%nrank
               coordatts(ivar+1)%global_dims(nrank) = nstepout
               coordatts(ivar+1)%dimids(nrank) = dimt
            ENDIF
            ivar = ivar + 1
         ELSE
!           --sigma levels            
            IF (iopt_grid_vtype.LE.2) THEN
               coordatts(ivar+1)%nrank = 1
               coordatts(ivar+1)%global_dims(1) = nzout
               coordatts(ivar+1)%dimids(1) = dimz
               IF (time_grid) THEN
                  formula_terms = 'sigma: lev eta: zeta depth: depout'
               ELSE
                  formula_terms = 'sigma: lev depth: depout'
               ENDIF
            ELSEIF (iopt_grid_vtype_transf.EQ.21) THEN
               coordatts(ivar+1)%nrank = 3
               coordatts(ivar+1)%global_dims(1:3) = (/ncout,nrout,nzout/)
               coordatts(ivar+1)%dimids(1:3) = (/dimx,dimy,dimz/)
               IF (time_grid) THEN
                  formula_terms = 's: lev eta: zeta depth: depout '//& 
                                & 'a: theta_SH b: b_SH depth_c : hcrit_SH'
               ELSE
                  formula_terms = 's: lev depth: depout '//& 
                                & 'a: theta_SH b: b_SH depth_c : hcrit_SH'
               ENDIF
            ENDIF
            ivar = ivar + 1
            IF (time_grid) THEN
               IF (gridded) THEN
                  IF (.NOT.packing) THEN
                     coordatts(ivar+1)%nrank = 3
                     coordatts(ivar+1)%global_dims(1:3) = &
                                                     & (/ncout,nrout,nstepout/)
                     coordatts(ivar+1)%dimids(1:3) = (/dimx,dimy,dimt/)
                  ELSE
                     coordatts(ivar+1)%nrank = 2
                     coordatts(ivar+1)%global_dims(1:2) = (/nowetout,nstepout/)
                     coordatts(ivar+1)%dimids(1:2) = (/dimspt,dimt/)
                  ENDIF
               ELSE
                  coordatts(ivar+1)%nrank = 2
                  coordatts(ivar+1)%global_dims(1:2) = (/nostats,nstepout/)
                  coordatts(ivar+1)%dimids(1:2) = (/dimstat,dimt/)
               ENDIF
               ivar = ivar + 1
            ENDIF
            IF (iopt_grid_vtype_transf.EQ.21) THEN
               coordatts(ivar+1:ivar+3)%nrank = 0
               ivar = ivar + 3
            ENDIF
         ENDIF
      ENDIF
      IF (.NOT.gridded) THEN
         coordatts(ivar+1)%nrank = 2
         coordatts(ivar+1)%global_dims(1:2) = (/lendesc,nostats/)
         coordatts(ivar+1)%dimids(1:2) = (/dimname,dimstat/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!
!2.7.4 Define coordinate variables (netCDF)
!------------------------------------------
!

   IF (outform.EQ.'N') THEN
      IF (time_format.EQ.0) THEN
         CALL cf90_def_var(iodat,coordatts(1)%f90_name,char_NF90,&
                        & (/dimtlen,dimt/),varid)
      ELSE
         CALL cf90_def_var(iodat,coordatts(1)%f90_name,double_NF90,(/dimt/),&
                         & varid)
      ENDIF
      ivar_274: DO ivar=2,nocoords
         nrank = coordatts(ivar)%nrank
         CALL cf90_def_var(iodat,coordatts(ivar)%f90_name,&
                         & coordatts(ivar)%data_type,&
                         & coordatts(ivar)%dimids(1:nrank),varid)
      ENDDO ivar_274
   ENDIF

!
!2.7.5 Variable ids and coordinate attributes
!--------------------------------------------
!

   ivar_275: DO ivar=1,nocoords
      SELECT CASE (TRIM(coordatts(ivar)%f90_name))
         CASE ('time'); filepars%timeid = ivar
         CASE ('eta','z')
            filepars%zvarid = ivar
            IF (.NOT.gridded.OR.iopt_grid_htype.EQ.3) THEN
               coordatts(ivar)%coordinates = TRIM(coordinates)
            ENDIF
         CASE ('depout')
            IF (.NOT.gridded.OR.iopt_grid_htype.EQ.3) THEN
               coordatts(ivar)%coordinates = TRIM(coordinates)
            ENDIF
         CASE ('theta_SH')
            theta_id = ivar
         CASE ('b_SH')
            b_id = ivar
         CASE ('hcrit_SH')
            hcrit_id = ivar
         CASE ('station_name'); statid = ivar
      END SELECT
   ENDDO ivar_275

!
!2.7.6 Store output dimensions
!-----------------------------
!
!  ---number of dimensions
   nodims = 1
   IF (nodim.GT.0) THEN
      IF (gridded) THEN
         nodims = nodims + 2
         IF (packing) nodims = nodims + 1
      ELSE
         nodims = nodims + 1
      ENDIF
      IF (nodim.EQ.3) nodims = nodims + 1
   ENDIF

!  ---allocate
   ALLOCATE (outdims(nodims),STAT=errstat)
   CALL error_alloc('outdims',1,(/nodims/),kndint)
   ALLOCATE (dimnames(nodims),STAT=errstat)
   CALL error_alloc('dimnames',1,(/nodims/),kndchar,lenstr=lenname)

!  ---values and names
   outdims(1) = nstepout
   dimnames(1) = 'time'
   ivar = 1
   IF (nodim.GT.0) THEN
      IF (gridded) THEN
         IF (iopt_grid_htype.LE.2) THEN
            xname = MERGE ('x  ','lon',iopt_grid_sph.EQ.0)
            yname = MERGE ('y  ','lat',iopt_grid_sph.EQ.0)
         ELSE
            xname = 'nc'; yname = 'nr'
         ENDIF
         outdims(ivar+1:ivar+2) = (/ncout,nrout/)
         dimnames(ivar+1:ivar+2) = (/xname,yname/)
         ivar = ivar + 2
         IF (packing) THEN
            outdims(ivar+1) = nowetout
            dimnames(ivar+1) = 'seapoint'
            ivar = ivar + 1
         ENDIF
      ELSE
         outdims(ivar+1) = nostats
         dimnames(ivar+1) = 'station'
         ivar = ivar + 1
      ENDIF
      IF (nodim.EQ.3) THEN
         outdims(ivar+1) = nzout
         dimnames(ivar+1) = MERGE('z  ', 'lev',vcoord.EQ.1)
         ivar = ivar + 1
      ENDIF
   ENDIF
   
!
!2.8 Output variables
!--------------------
!
!2.8.1 Initialise
!----------------
!

   ALLOCATE (filevars(novars),STAT=errstat)
   CALL error_alloc_struc('filevars',1,(/novars/),'VariableAtts')
   CALL varatts_init(filevars)
   iivar = 0
   ivar_281: DO ivar=1,nfvars
      l = indvars(ivar)
      IF (l.GT.0) THEN
         IF (varatts(l)%nrank.EQ.nodim) THEN
            iivar = iivar + 1
            filevars(iivar) = varatts(l)
         ENDIF
      ENDIF
   ENDDO ivar_281

!
!2.8.2  Variable attributes
!--------------------------
!
   
   ivar_282: DO ivar=1,novars

!     ---rank
      SELECT CASE (nodim)
         CASE (0)
            filevars(ivar)%nrank = 1
         CASE (2)
            filevars(ivar)%nrank = MERGE (3,2,.NOT.packing.AND.gridded)
         CASE (3)
            filevars(ivar)%nrank = MERGE (4,3,.NOT.packing.AND.gridded)
      END SELECT

!     ---data type
      IF (outform.EQ.'N') THEN
         filevars%data_type = MERGE(real_NF90,double_NF90,floattype.EQ.'S')
      ELSE
         filevars%data_type = MERGE(real_type,rlong_type,floattype.EQ.'S')
      ENDIF
      
!     ---shape, dimension IDs
      SELECT CASE (nodim)
         CASE (0)
            filevars(ivar)%global_dims(1) = nstepout
            filevars(ivar)%dimids(1) = dimt 
         CASE (2)
            IF (gridded) THEN
               IF (.NOT.packing) THEN
                  filevars(ivar)%global_dims(1:3) = (/ncout,nrout,nstepout/)
                  filevars(ivar)%dimids(1:3) = (/dimx,dimy,dimt/)
               ELSE
                  filevars(ivar)%global_dims(1:2) = (/nowetout,nstepout/)
                  filevars(ivar)%dimids(1:2) = (/dimspt,dimt/)
               ENDIF
            ELSE
               filevars(ivar)%global_dims(1:2) = (/nostats,nstepout/)
               filevars(ivar)%dimids(1:2) = (/dimstat,dimt/)
            ENDIF
         CASE (3)
            IF (gridded) THEN
               IF (.NOT.packing) THEN
                  filevars(ivar)%global_dims(1:4) = &
                                             & (/ncout,nrout,nzout,nstepout/)
                  filevars(ivar)%dimids(1:4) = (/dimx,dimy,dimz,dimt/)
               ELSE
                  filevars(ivar)%global_dims(1:3) = (/nowetout,nzout,nstepout/)
                  filevars(ivar)%dimids(1:3) = (/dimspt,dimz,dimt/)
               ENDIF
            ELSE
               filevars(ivar)%global_dims(1:3) = (/nostats,nzout,nstepout/)
               filevars(ivar)%dimids(1:3) = (/dimstat,dimz,dimt/)
            ENDIF
      END SELECT

!     ---standard_name, long_name, units
      comment = 'standard name constructed following the netCDF CF-1.6 '//&
              & 'guidelines'
      SELECT CASE (file_type)
         CASE ('TS','TA'); comment = ''
         CASE ('HR')
            filevars(ivar)%standard_name = 'residual_'//&
                                      & TRIM(filevars(ivar)%standard_name) 
            filevars(ivar)%long_name = 'Residual '//&
                                     & TRIM(filevars(ivar)%long_name)
            IF (TRIM(filevars(ivar)%vector_name).NE.'') THEN
               filevars(ivar)%vector_name = 'Residual '//&
                                        & TRIM(filevars(ivar)%vector_name)
            ENDIF
         CASE ('HA')
            iifreq = ifreqsharm(iset,ifreq)
            filevars(ivar)%standard_name = 'amplitude_of_'//&
	                              & TRIM(harm_freq_names(iifreq))//'_'//&
                                      & TRIM(filevars(ivar)%standard_name)
            filevars(ivar)%long_name = 'Amplitude of '//&
                                     & TRIM(harm_freq_names(iifreq))//'-'//&
                                     & TRIM(filevars(ivar)%long_name)
            IF (TRIM(filevars(ivar)%vector_name).NE.'') THEN
               filevars(ivar)%vector_name = 'Amplitude of '//&
                                        & TRIM(filevars(ivar)%vector_name)
            ENDIF
         CASE ('HP')
            iifreq = ifreqsharm(iset,ifreq)
            filevars(ivar)%standard_name = 'phase _of_'//&
	                              & TRIM(harm_freq_names(iifreq))//'_'//&
                                      & TRIM(filevars(ivar)%standard_name)
            filevars(ivar)%long_name = 'Phase of '//&
                                     & TRIM(harm_freq_names(iifreq))//'-'//&
                                     & TRIM(filevars(ivar)%long_name)
            filevars(ivar)%units = MERGE('degrees','radian ',DeGreesout)
            IF (TRIM(filevars(ivar)%vector_name).NE.'') THEN
               filevars(ivar)%vector_name = 'Phase of '//&
                                        & TRIM(filevars(ivar)%vector_name)
            ENDIF
         CASE ('HE')
            iifreq = ifreqsharm(iset,ifreq)
            filevars(ivar)%standard_name = TRIM(harm_freq_names(iifreq))//'_'//&
                                      & TRIM(filevars(ivar)%standard_name)
            IF (TRIM(filevars(ivar)%vector_name).NE.'') THEN
               filevars(ivar)%vector_name = &
                                      & TRIM(harm_freq_names(iifreq))//'-'//&
                                      & TRIM(filevars(ivar)%vector_name)
            ENDIF
            filevars(ivar)%long_name = TRIM(harm_freq_names(iifreq))//'-'//&
                                     & TRIM(filevars(ivar)%long_name)
            IF (TRIM(filevars(ivar)%vector_name).NE.'') THEN
               filevars(ivar)%vector_name = &
                                      & TRIM(harm_freq_names(iifreq))//'-'//&
                                      & TRIM(filevars(ivar)%vector_name)
            ENDIF
      END SELECT

!     ---comment
      filevars(ivar)%comment = comment

!     ---coordinates
      IF (.NOT.gridded.OR.iopt_grid_htype.EQ.3) THEN
         filevars(ivar)%coordinates = 'lon lat lev time'
      ENDIF

!     ---fill value
      filevars%fill = filepars%fill
      filevars%fill_value = fill_value

   ENDDO ivar_282

!
!2.8.3 Define output variables (netCDF)
!-------------------------------------
!

   IF (outform.EQ.'N') THEN
      ivar_283: DO ivar=1,novars
         nrank = filevars(ivar)%nrank
         CALL cf90_def_var(iodat,filevars(ivar)%f90_name,&
                         & filevars(ivar)%data_type,&
                         & filevars(ivar)%dimids(1:nrank),varid)
      ENDDO ivar_283
   ENDIF

!
!2.9 Write metadata
!------------------
!
!2.9.1 ASCII/binary/info
!-----------------------
!

   ifil_291: DO ifil=1,2

      IF (ifil.EQ.1.AND.outform.EQ.'N') CYCLE ifil_291
      IF (ifil.EQ.2.AND.(.NOT.info)) EXIT ifil_291
      iunit = MERGE(iodat,ioinfo,ifil.EQ.1)
      dataform = MERGE(outform,'A',ifil.EQ.1)

!     ---global attributes
      noglbatts = COUNT(LEN_TRIM(glbatts%value).GT.0)
      CALL conv_to_chars(cvals(1),noglbatts)
      CALL write_metadata_line(iunit,cvals(1:1),'nGlbatts',dataform)
      iglb_2911: DO iglb=1,numglbatts
         IF (TRIM(glbatts(iglb)%value).NE.''.AND.&
           & TRIM(glbatts(iglb)%value).NE.'netcdf') THEN
            cvals(1) = TRIM(glbatts(iglb)%value)
            CALL write_metadata_line(iunit,cvals(1:1),&
                                   & TRIM(glbatts(iglb)%name),dataform)
         ENDIF
      ENDDO iglb_2911
      IF (.NOT.gridded) THEN
         IF (nodim.EQ.2) THEN
            cvals(1) = 'timeSeries'
         ELSEIF (nodim.EQ.3) THEN
            cvals(1) = 'timeSeriesProfile'
         ENDIF
         CALL write_metadata_line(iunit,cvals(1:1),'featureType',dataform)
      ENDIF
      
!     ---output info
      CALL conv_to_chars(cvals(1),nodim)
      CALL write_metadata_line(iunit,cvals(1:1),'nodim',dataform)
      CALL conv_to_chars(cvals(1),nocoords)
      CALL write_metadata_line(iunit,cvals(1:1),'nCoordinates',dataform)
      CALL conv_to_chars(cvals(1),novars)
      CALL write_metadata_line(iunit,cvals(1:1),'nVariables',dataform)
      cvals(1) = floattype
      CALL write_metadata_line(iunit,cvals(1:1),'Floattype',dataform)
      IF (vcoord.EQ.2.AND.nodim.EQ.3) THEN
         cvals(1) = TRIM(formula_terms)
         CALL write_metadata_line(iunit,cvals(1:1),'formula_terms',dataform)
      ENDIF
      CALL conv_to_chars(cvals(1),nodims)
      CALL write_metadata_line(iunit,cvals(1:1),'nDimensions',dataform)
      CALL write_metadata_line(iunit,dimnames(1:nodims),'Dimnames',dataform)
      CALL conv_to_chars(cvals(1:nodims),outdims(1:nodims))
      CALL write_metadata_line(iunit,cvals(1:nodims),'Dimensions',dataform)
      cvals(1) = outgpars%startdate
      CALL write_metadata_line(iunit,cvals(1:1),'StartDate',dataform)
      CALL conv_to_chars(cvals(1),deltout)
      CALL write_metadata_line(iunit,cvals(1:1),'time_step',dataform)
      cvals(1) = TRIM(coordatts(1)%units)
      CALL write_metadata_line(iunit,cvals(1:1),'time_format',dataform)
      
!     ---coordinate attributes
      ivar_2912: DO ivar=1,nocoords
         CALL conv_to_chars(cvals(1),ivar)
         nvals = 11 + coordatts(ivar)%nrank
         CALL conv_to_chars_modvars(cvals(2:nvals),coordatts(ivar))
         CALL write_metadata_line(iunit,cvals(1:nvals),'coordatts',dataform)
      ENDDO ivar_2912

!     ---variable attributes
      ivar_2913: DO ivar=1,novars
         CALL conv_to_chars(cvals(1),nocoords+ivar)
         nvals = 11 + filevars(ivar)%nrank
         CALL conv_to_chars_modvars(cvals(2:nvals),filevars(ivar))
         CALL write_metadata_line(iunit,cvals(1:nvals),'varatts',dataform)
      ENDDO ivar_2913
      
!     ---SH-parameters
      IF (vcoord.EQ.2.AND.iopt_grid_vtype.EQ.21) THEN
         CALL conv_to_chars(cvals(1:3),(/theta_SH,b_SH,hcrit_SH/))
         CALL write_metadata_line(iunit,cvals(1:3),'SH_parameters',dataform)
      ENDIF

!     ---station names
      IF (.NOT.gridded) THEN
         istat_2914: DO istat=1,nostats
            CALL conv_to_chars(cvals(1),istat)
            cvals(2) = outlocs(istat)%name
            CALL write_metadata_line(iunit,cvals(1:2),'StationNames',dataform)
         ENDDO istat_2914
      ENDIF

!     ---end of header
      cline = '#'
      IF (dataform.EQ.'A') THEN
         WRITE (iunit,'(A)') TRIM(cline)
      ELSEIF (dataform.EQ.'U') THEN
         WRITE (iunit) cline
      ENDIF

   ENDDO ifil_291

!  ---close
   IF (info) CALL close_file(ioinfo,'A',TRIM(infofile))

!
!2.9.2 Netcdf
!------------
!

   IF (outform.EQ.'N') THEN

!     ---global attributes
      iglb_2921: DO iglb=1,numglbatts
         IF (TRIM(glbatts(iglb)%value).NE.'') THEN
            l = LEN_TRIM(glbatts(iglb)%value)
            CALL cf90_put_att_chars(iodat,global_NF90,&
                                  & TRIM(glbatts(iglb)%name),l,&
                                  & TRIM(glbatts(iglb)%value))
         ENDIF
      ENDDO iglb_2921
      IF (.NOT.gridded) THEN
         IF (nodim.EQ.2) THEN
            CALL cf90_put_att_chars(iodat,global_NF90,'featureType',10,&
                                 & 'timeSeries')
         ELSEIF (nodim.EQ.3) THEN
            CALL cf90_put_att_chars(iodat,global_NF90,'featureType',17,&
                                 & 'timeSeriesProfile')
         ENDIF
      ENDIF

!     ---coordinate attributes
      varid = 0
      ivar_2922: DO ivar=1,nocoords
         varid = varid + 1
         CALL cf90_put_att_chars(iodat,varid,'standard_name',lendesc,&
                               & coordatts(ivar)%standard_name)
         CALL cf90_put_att_chars(iodat,varid,'long_name',lendesc,&
                               & coordatts(ivar)%long_name)
         IF (TRIM(coordatts(ivar)%units).NE.'') THEN
            CALL cf90_put_att_chars(iodat,varid,'units',lenunit,&
                                  & coordatts(ivar)%units)
         ENDIF
         IF (coordatts(ivar)%fill) THEN
            CALL cf90_put_att(iodat,varid,'_FillValue',fill_value)
            CALL cf90_put_att(iodat,varid,'missing_value',fill_value)
         ENDIF
         IF (TRIM(coordatts(ivar)%comment).NE.'') THEN
            CALL cf90_put_att_chars(iodat,varid,'comment',lendesc,&
                                  & coordatts(ivar)%comment)
         ENDIF
         IF (TRIM(coordatts(ivar)%axis).NE.'') THEN
            CALL cf90_put_att_chars(iodat,varid,'axis',1,&
                                  & coordatts(ivar)%axis)
         ENDIF
!        --other
         SELECT CASE (TRIM(coordatts(ivar)%f90_name))
            CASE ('time')
               CALL cf90_put_att_chars(iodat,varid,'calendar',9,'gregorian')
            CASE ('lev')
               CALL cf90_put_att_chars(iodat,varid,'positive',2,'up')
               l = LEN_TRIM(formula_terms)
               CALL cf90_put_att_chars(iodat,varid,'formula_terms',l,&
                                     & TRIM(formula_terms))
            CASE ('z')
               CALL cf90_put_att_chars(iodat,varid,'positive',2,'up')
            CASE ('seapoint')
               IF (iopt_grid_sph.EQ.0) THEN
                  CALL cf90_put_att_chars(iodat,varid,'compress',3,'x y')
               ELSE
                  CALL cf90_put_att_chars(iodat,varid,'compress',7,'lon lat')
               ENDIF
            CASE ('station_name')
               CALL cf90_put_att_chars(iodat,varid,'cf_role',13,&
                                    & 'timeseries_id')
         END SELECT
         IF (TRIM(coordatts(ivar)%coordinates).NE.'') THEN
            CALL cf90_put_att_chars(iodat,varid,'coordinates',lendesc,&
                                  & coordinates)
         ENDIF
      ENDDO ivar_2922

!     ---variable attributes
      ivar_2923: DO ivar=1,novars
         varid = varid + 1
         CALL cf90_put_att_chars(iodat,varid,'standard_name',lendesc,&
                               & filevars(ivar)%standard_name)
         CALL cf90_put_att_chars(iodat,varid,'long_name',lendesc,&
                               & filevars(ivar)%long_name)
         IF (TRIM(filevars(ivar)%units).NE.'') THEN
            CALL cf90_put_att_chars(iodat,varid,'units',lenunit,&
                                  & filevars(ivar)%units)
         ENDIF
         IF (TRIM(filevars(ivar)%vector_name).NE.'') THEN
            CALL cf90_put_att_chars(iodat,varid,'vector_name',lendesc,&
                                  & filevars(ivar)%vector_name)
         ENDIF
         IF (filepars%fill) THEN
            CALL cf90_put_att(iodat,varid,'_FillValue',fill_value)
            CALL cf90_put_att(iodat,varid,'missing_value',fill_value)
         ENDIF
         IF (TRIM(filevars(ivar)%comment).NE.'') THEN
            CALL cf90_put_att_chars(iodat,varid,'comment',lendesc,&
                                  & filevars(ivar)%comment)
         ENDIF
         IF (TRIM(filevars(ivar)%coordinates).NE.'') THEN
            CALL cf90_put_att_chars(iodat,varid,'coordinates',lendesc,&
                                  & coordinates)
         ENDIF
      ENDDO ivar_2923

!     ---leave define mode
      CALL cf90_enddef(iodat)

!     ---SH-parameters
      IF (iopt_grid_vtype.EQ.21) THEN
         CALL cf90_put_var(iodat,theta_id,theta_SH,0)
         CALL cf90_put_var(iodat,b_id,b_SH,0)
         CALL cf90_put_var(iodat,hcrit_id,hcrit_SH,0)
      ENDIF

!      --station names
      IF (.NOT.gridded) THEN
         CALL cf90_put_var_chars(iodat,statid,outlocs(1:nostats)%name,0)
      ENDIF
      
   ENDIF

!
!2.10 Abort if needed
!--------------------
!

   IF (nerrs.GT.0) CALL error_file(ierrno_write,filepars=filepars)

!
!2.11 Deallocate
!---------------
!

   DEALLOCATE (coordatts,dimnames,filevars,outdims)

!
!2.12 Store
!----------
!

   SELECT CASE (file_type)
      CASE ('TS')
         IF (nodim.EQ.0) tsr0d(iset) = filepars
         IF (nodim.EQ.2) tsr2d(iset) = filepars
         IF (nodim.EQ.3) tsr3d(iset) = filepars
     CASE ('TA')
         IF (nodim.EQ.0) avr0d(iset) = filepars
         IF (nodim.EQ.2) avr2d(iset) = filepars
         IF (nodim.EQ.3) avr3d(iset) = filepars
     CASE ('HR')
         IF (nodim.EQ.0) res0d(iset) = filepars
         IF (nodim.EQ.2) res2d(iset) = filepars
         IF (nodim.EQ.3) res3d(iset) = filepars
     CASE ('HA')
         IF (nodim.EQ.0) amp0d(iset,ifreq) = filepars
         IF (nodim.EQ.2) amp2d(iset,ifreq) = filepars
         IF (nodim.EQ.3) amp3d(iset,ifreq) = filepars
     CASE ('HP')
         IF (nodim.EQ.0) pha0d(iset,ifreq) = filepars
         IF (nodim.EQ.2) pha2d(iset,ifreq) = filepars
         IF (nodim.EQ.3) pha3d(iset,ifreq) = filepars
     CASE ('HE')
         IF (nodim.EQ.2) ell2d(iset,ifreq) = filepars
         IF (nodim.EQ.3) ell3d(iset,ifreq) = filepars
   END SELECT

ENDDO nodim_200
ENDDO ifreq_200
ENDDO iset_200

!
!3. Deallocate
!-------------
!

DEALLOCATE (indvars,indstats,outstatlocs,outlocs,varatts)

CALL log_timer_out(npcc,itm_output)


RETURN

END SUBROUTINE inout_atts_out

!========================================================================

SUBROUTINE inquire_modfil_atts(iounit,ioform,varid,nrank,data_type,global_dims,&
                             & f90_name,standard_name,comment,long_name,units,&
                             & vector_name,fill,fill_value,reset)
!************************************************************************
!
! *inquire_modfil_atts* Returns information about names and attributes of
!                       coordinate and data variables from a file in
!                       CF-standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.10.1
!
! Description -
!
! Module calls - cf90_get_att, cf90_get_att_chars, cf90_get_var_chars,
!                cf90_inquire_attribute, cf90_inquire_dimension,
!                cf90_inquire_variable, cf90_inq_attname, cf90_inq_varid,
!                conv_from_chars, read_metadata_line
!
!************************************************************************
!
USE cf90_routines, ONLY: cf90_get_att, cf90_get_att_chars, cf90_get_var_chars, &
                       & cf90_inquire_attribute, cf90_inquire_dimension, &
                       & cf90_inquire_variable, cf90_inq_attname, &
                       & cf90_inq_varid 
USE cif_routines, ONLY: conv_from_chars

!
!*Arguments
!

LOGICAL, INTENT(IN), OPTIONAL :: reset
CHARACTER (LEN=lenform), INTENT(IN) :: ioform
INTEGER, INTENT(IN) :: iounit, varid
INTEGER, INTENT(OUT), OPTIONAL :: data_type, nrank
LOGICAL, INTENT(OUT), OPTIONAL :: fill
INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(5) :: global_dims
REAL (KIND=kndrlong), INTENT(OUT), OPTIONAL :: fill_value
CHARACTER (LEN=lenname), INTENT(OUT), OPTIONAL :: f90_name
CHARACTER (LEN=lendesc), INTENT(OUT), OPTIONAL :: comment, long_name, &
                                                & standard_name, vector_name
CHARACTER (LEN=lenunit), INTENT(OUT), OPTIONAL :: units

!
! Name           Type    Purpose
!------------------------------------------------------------------------------
!*iounit*        INTEGER File unit number
!*ioform*        CHAR    File format
!*varid*         INTEGER Variable id
!*nrank*         INTEGER Variable array rank
!*data_type*     INTEGER Data type
!*global_dims*   INTEGER (Global) variable array shape 
!*f90_name*      CHAR    FORTRAN name
!*standard_name* CHAR    netCDF CF compliant standard name
!*comment*       CHAR    Miscellaneous information about the data or methods
!*long_name*     CHAR    Long descriptive name
!*units*         CHAR    Variable unit
!*vector_name*   CHAR    Associated vector name
!*fill*          LOGICAL .TRUE. if a fill value is defined
!*fill_value*    DOUBLE  Fill value
!*reset*         LOGICAL Rewinds the file if PRESENT and .TRUE.
!                        (ASCII and binary format only)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=lendesc), DIMENSION(17) :: cvals
CHARACTER (LEN=lenname) :: cname, fname, name
CHARACTER (LEN=lendesc) :: commt, lname, sname, vname
CHARACTER (LEN=lenunit) :: unit
LOGICAL :: comflag, fillx, vecflag
INTEGER :: datatype, idim, istat, ivar, l, linenum, n, natts, nodim, numvars, &
         & xtype
INTEGER, DIMENSION(5) :: dimids, ngdims
REAL :: fillvalue_s
REAL (KIND=kndrlong) :: fillvalue


procname(pglev+1) = 'inquire_modfil_atts'
CALL log_timer_in()

!
!1. ASCII or binary format
!-------------------------
!

IF (ioform.EQ.'A'.OR.ioform.EQ.'U') THEN

   REWIND (iounit)

   istat = 1; linenum = 0
   DO WHILE (istat.EQ.1)
      CALL read_metadata_line(iounit,ioform,cvals,numvars,cname,linenum,istat)
      SELECT CASE (TRIM(cname))
         CASE ('coordatts','varatts') 
            CALL conv_from_chars(cvals(1),ivar,1)
            IF (ivar.EQ.varid) THEN
               CALL conv_from_chars(cvals(2),datatype,2)
               CALL conv_from_chars(cvals(3),nodim,3)
               l_110: DO l=1,nodim
                  CALL conv_from_chars(cvals(l+3),ngdims(l),l+3)
               ENDDO l_110
               l = nodim + 3
               CALL conv_from_chars(cvals(l+1),fname,l+1)
               CALL conv_from_chars(cvals(l+2),sname,l+2)
               CALL conv_from_chars(cvals(l+3),commt,l+3)
               CALL conv_from_chars(cvals(l+4),lname,l+4)
               CALL conv_from_chars(cvals(l+5),unit,l+5)
               CALL conv_from_chars(cvals(l+6),vname,l+6)
               CALL conv_from_chars(cvals(l+7),fillx,l+7)
               IF (fillx) THEN
                  CALL conv_from_chars(cvals(l+8),fillvalue,l+8)
               ELSE
                  fillvalue = 0.0
               ENDIF
            ENDIF
      END SELECT

   ENDDO

   IF (PRESENT(reset)) THEN
      IF (reset) REWIND (iounit)
   ENDIF

!
!2. NetCDF format
!----------------
!

ELSEIF (ioform.EQ.'N') THEN

   fillx = .FALSE.; vecflag = .FALSE.
   CALL cf90_inquire_variable(iounit,varid,name=fname,xtype=xtype,ndims=nodim,&
                            & natts=natts)
   CALL cf90_inquire_variable(iounit,varid,dimids=dimids(1:nodim))
   IF (xtype.EQ.char_NF90) THEN
      datatype = char_type
   ELSEIF (xtype.EQ.int_NF90) THEN
         CALL cf90_get_att_chars(iounit,varid,'units',lenunit,lname)
         datatype = MERGE(log_type,int_type,TRIM(units).EQ.'log')
   ELSEIF (xtype.EQ.real_NF90) THEN
      datatype = real_type
   ELSEIF (xtype.EQ.double_NF90) THEN
      datatype = rlong_type
   ENDIF
   IF (PRESENT(global_dims)) THEN
      idim_211: DO idim=1,nodim
         CALL cf90_inquire_dimension(iounit,dimids(idim),len=ngdims(idim))
      ENDDO idim_211
   ENDIF
   IF (PRESENT(f90_name)) f90_name = fname
   n_210: DO n=1,natts
      CALL cf90_inq_attname(iounit,varid,n,name=name)
      SELECT CASE (TRIM(name))
         CASE ('standard_name')
            CALL cf90_get_att_chars(iounit,varid,'standard_name',lendesc,sname)
         CASE ('comment')
            CALL cf90_get_att_chars(iounit,varid,'comment',lendesc,commt)
            comflag = .TRUE.
         CASE ('long_name')
            CALL cf90_get_att_chars(iounit,varid,'long_name',lendesc,lname)
         CASE ('units')
            CALL cf90_get_att_chars(iounit,varid,'units',lenunit,unit)
         CASE ('vector_name')
            CALL cf90_get_att_chars(iounit,varid,'vector_name',lendesc,vname)
            vecflag = .TRUE.
         CASE ('_FillValue')
            fillx = .TRUE.
            CALL cf90_inquire_attribute(iounit,varid,'_FillValue',xtype=xtype)
            IF (xtype.EQ.real_NF90) THEN
               CALL cf90_get_att(iounit,varid,'_FillValue',fillvalue_s)
               fillvalue = fillvalue_s
            ELSE
               CALL cf90_get_att(iounit,varid,'_FillValue',fillvalue)
            ENDIF
      END SELECT
   ENDDO n_210

   IF (.NOT.comflag) commt = ''
   IF (.NOT.vecflag) vname = ''
   IF (.NOT.fillx) fillvalue = 0.0

ENDIF

!
!3. Store
!--------
!

IF (PRESENT(nrank)) nrank = nodim
IF (PRESENT(data_type)) data_type = datatype
IF (PRESENT(global_dims)) THEN
   global_dims = 0
   global_dims(1:nodim) = ngdims(1:nodim)
ENDIF
IF (PRESENT(f90_name)) f90_name = fname
IF (PRESENT(standard_name)) standard_name = sname
IF (PRESENT(comment)) comment = commt
IF (PRESENT(long_name)) long_name = lname
IF (PRESENT(units)) units = unit
IF (PRESENT(vector_name)) vector_name = vname
IF (PRESENT(fill)) fill = fillx 
IF (PRESENT(fill_value)) fill_value = fillvalue

CALL log_timer_out()


RETURN

END SUBROUTINE inquire_modfil_atts

!========================================================================

SUBROUTINE inquire_modfil_dims(iounit,ioform,ndimensions,ncoordinates,&
                             & nvariables,nglbatts,dimensions,dimnames,&
                             & maxrecs,floattype,packing,vcoord,ftype,timeid,&
                             & reset)
!************************************************************************
!
! *inquire_modfil_dims* Returns information about dimensions and number of
!                       variables, coordinates, global attributes from a file
!                       in CF-standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls -  cf90_get_att, cf90_get_att_chars, cf90_inquire,
!                 cf90_inquire_dimension, cf90_inquire_variable,
!                 cf90_inq_attname, conv_from_chars, error_alloc,
!                 read_metadata_line
!
!************************************************************************
!
USE cf90_routines, ONLY: cf90_get_att, cf90_get_att_chars, cf90_inquire, &
                       & cf90_inquire_dimension, cf90_inquire_variable, &
                       & cf90_inq_attname
USE cif_routines, ONLY: conv_from_chars
USE error_routines, ONLY: error_alloc

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: reset
CHARACTER (LEN=lenform), INTENT(IN) :: ioform
CHARACTER (LEN=lendesc), OPTIONAL :: ftype
INTEGER, INTENT(IN) :: iounit
LOGICAL, INTENT(OUT), OPTIONAL :: packing
CHARACTER (LEN=1), INTENT(OUT), OPTIONAL :: floattype
INTEGER, INTENT(OUT), OPTIONAL ::  maxrecs, nglbatts, ncoordinates, &
                                 & ndimensions, nvariables, timeid, vcoord
INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(:) :: dimensions
CHARACTER (LEN=lenname), INTENT(OUT), OPTIONAL, DIMENSION(:) :: dimnames

!
! Name            Type     Purpose
!------------------------------------------------------------------------------
!*iounit*         INTEGER File unit (or NerCDF ID)
!*ioform*         CHAR    File format
!*ndimensions*    INTEGER Returned number of dimensions
!*ncoordinates*   INTEGER Returned number of coordinate variables
!*nvariables*     INTEGER Returned number of non-coordinate variables
!*nglbatts*       INTEGER Returned number of global attributes
!*dimensions*     INTEGER Returned vector with dimension values
!*dimnames*       INTEGER Returned vector with dimension names
!*maxrecs*        INTEGER Returned number of time records (if any)
!*floattype*      CHAR    Returned data type for real variables ('S' or 'D')
!*packing*        LOGICAL Returns .TRUE. if data are in compressed ("packed")
!                         format
!*vcoord*         INTEGER Returns type of vertical coordinate
!                       = 1 => Z-coordinate
!                       = 2 => sigma-coordinate
!*reset*          LOGICAL Rewinds the file if PRESENT and .TRUE.
!                         (ASCII and binary format only)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packingx
CHARACTER (LEN=1) :: floattypex
CHARACTER (LEN=lenname) :: cname, f90_name
CHARACTER (LEN=lendesc) :: attname, ftypex
CHARACTER (LEN=lendesc), DIMENSION(15) :: cvals
INTEGER :: idim, ipack, istat, ivar, linenum, n, nocoords, nodim, nodims, &
         & noglbatts, norecs, novars, numvars, timeidx, xtype, vcoordx
CHARACTER (LEN=lenname), SAVE, ALLOCATABLE, DIMENSION(:) :: namedims
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: outdims


procname(pglev+1) = 'inquire_modfil_dims'
CALL log_timer_in()

nodims = 0; norecs = 0; packingx = .FALSE.; vcoordx = 0; ftypex = ' '

!
!1. ASCII or binary format
!-------------------------
!

IF (ioform.EQ.'A'.OR.ioform.EQ.'U') THEN

   REWIND (iounit)

   istat = 1; linenum = 0
   DO WHILE (istat.EQ.1)

      CALL read_metadata_line(iounit,ioform,cvals,numvars,cname,linenum,istat)

      SELECT CASE (TRIM(cname))

!        ---number of dimensions
         CASE ('nDimensions')
            CALL conv_from_chars(cvals(1),nodims,1)
!           --allocate
            ALLOCATE (namedims(nodims),STAT=errstat)
            CALL error_alloc('namedims',1,(/nodims/),kndchar,lenstr=lenname)
            ALLOCATE (outdims(nodims),STAT=errstat)
            CALL error_alloc('outdims',1,(/nodims/),kndint)
            
!        ---number of coordinates
         CASE ('nCoordinates')
            CALL conv_from_chars(cvals(1),nocoords,1)

!        ---number of variables            
         CASE ('nVariables')
            CALL conv_from_chars(cvals(1),novars,1)
            
!        ---number of global attributes
         CASE ('nGlbatts')
            CALL conv_from_chars(cvals(1),noglbatts,1)
            
!        ---dimension values and names            
         CASE ('Dimensions')
            idim_110: DO idim=1,nodims
               CALL conv_from_chars(cvals(idim),outdims(idim),idim)
            ENDDO idim_110
         CASE ('Dimnames')
            idim_120: DO idim=1,nodims
               CALL conv_from_chars(cvals(idim),namedims(idim),idim)
               IF (TRIM(namedims(idim)).EQ.'z') THEN
                  vcoordx = 1
               ELSEIF (TRIM(namedims(idim)).EQ.'lev') THEN
                  vcoordx = 2
               ENDIF
            ENDDO idim_120
!           --number of time records
            idim_130: DO idim=1,nodims
               IF (TRIM(namedims(idim)).EQ.'time') norecs = outdims(idim)
            ENDDO idim_130
            
!        ---real data type             
         CASE ('Floattype')
            CALL conv_from_chars(cvals(1),floattypex,1)

!        ---packing            
         CASE ('packing')
            CALL conv_from_chars(cvals(1),packingx,1)

!        ---featureType
         CASE ('featureType')
            CALL conv_from_chars(cvals(1),ftypex,1)

      END SELECT

   ENDDO

!  ---id of time coordinate         
   timeidx = 1   

!  ---rewind if needed         
   IF (PRESENT(reset)) THEN
      IF (reset) REWIND (iounit)
   ENDIF

!
!2. NetCDF format
!----------------
!

ELSEIF (ioform.EQ.'N') THEN

!  ---number of dimensions
   CALL cf90_inquire(iounit,ndimensions=nodims,nattributes=noglbatts)
   nodim = MERGE(0,2,nodims.EQ.1)

!  ---allocate
   ALLOCATE (namedims(nodims),STAT=errstat)
   CALL error_alloc('namedims',1,(/nodims/),kndchar,lenstr=lenname)
   ALLOCATE (outdims(nodims),STAT=errstat)
   CALL error_alloc('outdims',1,(/nodims/),kndint)

!  ---dimension values and names
   idim_210: DO idim=1,nodims
      CALL cf90_inquire_dimension(iounit,idim,name=namedims(idim),&
                                & len=outdims(idim))
      IF (TRIM(namedims(idim)).EQ.'z') THEN
         nodim = 3; vcoordx = 1
      ELSEIF (TRIM(namedims(idim)).EQ.'lev') THEN
         nodim = 3; vcoordx = 2
      ELSEIF (TRIM(namedims(idim)).EQ.'T'.OR.TRIM(namedims(idim)).EQ.'t'.OR.&
            & TRIM(namedims(idim)).EQ.'time') THEN
         norecs = outdims(idim)
      ENDIF
   ENDDO idim_210

!  ---packing (forcing files only)
   n_220: DO n=1,noglbatts
      CALL cf90_inq_attname(iounit,global_NF90,n,attname)
      IF (TRIM(attname).EQ.'packing') THEN
         CALL cf90_get_att(iounit,global_NF90,'packing',ipack)
         packingx = ipack.EQ.1
      ENDIF
      IF (TRIM(attname).EQ.'featureType') THEN
         CALL cf90_get_att_chars(iounit,global_NF90,'featureType',lendesc,&
                               & ftypex)
      ENDIF
   ENDDO n_220

!  ---number of coordinates/variables/global attributes
   CALL cf90_inquire(iounit,nvariables=numvars)
   novars = 0; nocoords = 0
   ivar_230: DO ivar=1,numvars
      CALL cf90_inquire_variable(iounit,ivar,name=f90_name,xtype=xtype)
      IF (TRIM(f90_name).NE.'time') THEN
         IF (xtype.EQ.real_NF90) floattypex = 'S'
         IF (xtype.EQ.double_NF90) floattypex = 'D'
      ELSE
         timeidx = ivar
      ENDIF
      SELECT CASE (TRIM(f90_name))
         CASE ('time'); nocoords = nocoords + 1
         CASE ('x','y'); nocoords = nocoords + 1
         CASE ('lon','lat'); nocoords = nocoords + 1
         CASE ('lev','zout'); nocoords = nocoords + 1
         CASE ('depout'); nocoords = nocoords + 1
         CASE ('eta')
            IF (nodim.EQ.3) THEN
               nocoords = nocoords + 1
            ELSE
               novars = novars + 1
            ENDIF
         CASE ('seapoint','station_name')
            nocoords = nocoords + 1
!           packing (user output only)            
            packingx = .TRUE.
         CASE ('theta_SH','b_SH','hcrit_SH'); nocoords = nocoords + 1
         CASE ('trajectory'); nocoords = nocoords + 1
         CASE DEFAULT; novars = novars + 1
      END SELECT
   ENDDO ivar_230
   
ENDIF

!
!3. Store
!--------
!

IF (PRESENT(ndimensions)) ndimensions = nodims
IF (PRESENT(ncoordinates)) ncoordinates = nocoords
IF (PRESENT(nvariables)) nvariables = novars
IF (PRESENT(nglbatts)) nglbatts = noglbatts
IF (PRESENT(dimensions)) dimensions = outdims
IF (PRESENT(dimnames)) dimnames = namedims
IF (PRESENT(maxrecs)) maxrecs = norecs
IF (PRESENT(floattype)) floattype = floattypex
IF (PRESENT(packing)) packing = packingx
IF (PRESENT(vcoord)) vcoord = vcoordx
IF (PRESENT(timeid)) timeid = timeidx
IF (PRESENT(ftype)) ftype = ftypex

!
!4. Deallocate
!-------------
!

IF (ALLOCATED(namedims)) DEALLOCATE (namedims)
IF (ALLOCATED(outdims)) DEALLOCATE (outdims)


CALL log_timer_out()

RETURN

END SUBROUTINE inquire_modfil_dims

!========================================================================

SUBROUTINE inquire_modfil_vars(iounit,ioform,station_names,SH_params,reset)
!************************************************************************
!
! *inquire_modfil_vars* Returns the values of some auxiliary coordinate
!                       variables from a file in CF-standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.7.1
!
! Description -
!
! Module calls - cf90_get_var, cf90_get_var_chars, conv_from_chars,
!                read_metadata_line
!
!************************************************************************
!
USE cf90_routines, ONLY: cf90_get_var, cf90_get_var_chars, cf90_inq_varid 
USE cif_routines, ONLY: conv_from_chars

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: reset
CHARACTER (LEN=lenform), INTENT(IN) :: ioform
INTEGER, INTENT(IN) :: iounit
CHARACTER (LEN=lendesc), INTENT(OUT), OPTIONAL, DIMENSION(:) :: station_names
REAL, INTENT(OUT), OPTIONAL, DIMENSION(3) :: SH_params 

!
! Name           Type    Purpose
!------------------------------------------------------------------------------
!*iounit*        INTEGER File unit number
!*ioform*        CHAR    File format
!*station_names* CHAR    Names of output stations
!*SH_params*     REL     Parameters for s-coordinate formulation
!*reset*         LOGICAL Rewinds the file if PRESENT and .TRUE.
!                        (ASCII and binary format only)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=lendesc), DIMENSION(14) :: cvals
CHARACTER (LEN=lenname) :: cname
INTEGER :: istat, linenum, n, numvars, varid


procname(pglev+1) = 'inquire_modfil_vars'
CALL log_timer_in()

!
!1. ASCII or binary format
!-------------------------
!

IF (ioform.EQ.'A'.OR.ioform.EQ.'U') THEN

   REWIND(iounit)

   istat = 1; linenum = 0
   DO WHILE (istat.EQ.1)

      CALL read_metadata_line(iounit,ioform,cvals,numvars,cname,linenum,istat)
      SELECT CASE (TRIM(cname))
         CASE ('StationNames')
            IF (PRESENT(station_names)) THEN
               CALL conv_from_chars(cvals(1),n,1)
               CALL conv_from_chars(cvals(2),station_names(n),1)
            ENDIF
         CASE ('SH_parameters')
            IF (PRESENT(SH_params)) THEN
               CALL conv_from_chars(cvals(1),SH_params(1),1)
               CALL conv_from_chars(cvals(1),SH_params(2),2)
               CALL conv_from_chars(cvals(1),SH_params(3),3)
            ENDIF
      END SELECT

   ENDDO

   IF (PRESENT(reset)) THEN
      IF (reset) REWIND(iounit)
   ENDIF

!
!2. NetCDF format
!----------------
!

ELSEIF (ioform.EQ.'N') THEN

   IF (PRESENT(station_names)) THEN
      CALL cf90_inq_varid(iounit,'station_name',varid)
      CALL cf90_get_var_chars(iounit,varid,station_names,0)
   ENDIF
   IF (PRESENT(SH_params)) THEN
      CALL cf90_inq_varid(iounit,'theta_SH',varid)
      CALL cf90_get_var(iounit,varid,SH_params(1),0)
      CALL cf90_inq_varid(iounit,'b_SH',varid)
      CALL cf90_get_var(iounit,varid,SH_params(2),0)
      CALL cf90_inq_varid(iounit,'hcrit_SH',varid)
      CALL cf90_get_var(iounit,varid,SH_params(3),0)
   ENDIF

ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE inquire_modfil_vars

!========================================================================

FUNCTION inquire_unit(iounit,ioform)
!************************************************************************
!
! *inquire_unit* Inquires whether file unit 'iunit' is connected to a file 
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_inquire
!
!************************************************************************
!
USE cf90_routines, ONLY: cf90_inquire

!
!*Arguments
!
CHARACTER (LEN=lenform) :: ioform
INTEGER, INTENT(IN) :: iounit
LOGICAL :: inquire_unit

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*iunit*   INTEGER  File unit
!*ioform*  CHAR    File format ('A', 'U' or 'N')
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag


IF (ioform.NE.'N') THEN
   INQUIRE (UNIT=iounit,OPENED=flag)
   inquire_unit = flag
ELSE
   CALL cf90_inquire(iounit,abort=.FALSE.)
   inquire_unit = errstat.EQ.noerr_NF90
ENDIF


RETURN

END FUNCTION inquire_unit

!========================================================================

SUBROUTINE monitor_files
!************************************************************************
!
! *monitor_files* Set parameters and open/close log, error and warning files
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.7.1
!
! Description - 
!
! Module calls - close_file, comms_comm_rank, mult_index, open_file
!
!************************************************************************
!
USE paralpars
USE timepars
USE utility_routines, ONLY: mult_index

!
!*Local variables
!
CHARACTER (LEN=12) :: cdloc


!
!1. Initialise parameters
!------------------------
! 

IF (parallel_set) THEN
   WRITE (cdloc,'(I12)') idlocglb; cdloc = ADJUSTL(cdloc)
ELSE
   cdloc = ' '
ENDIF

!
!2. Initialise at setup
!----------------------
!

IF (nt.EQ.0) THEN

!  ---set parameters
   nopenf = 0
   timer = levtimer.GT.0
   
!  ---initial log file
   loglev1 = levprocs_ini(iprocglb)
   loglev2 = MERGE(loglev1,0,exitlog)
   IF (loglev1.GT.0) THEN
      inilog_file = TRIM(inilog_file)
      IF (parallel_set) inilog_file = TRIM(inilog_file)//TRIM(cdloc)
      CALL open_file(iolog,inilog_file,'OUT','A')
   ENDIF

!  ---error file
   errchk = levprocs_err(iprocglb).GT.0
   IF (errchk) THEN
      errlog_file = TRIM(errlog_file)
      IF (parallel_set) errlog_file = TRIM(errlog_file)//TRIM(cdloc)
      CALL open_file(ioerr,errlog_file,'OUT','A')
   ENDIF

!  ---warning file
   warnflag = warning.AND.mastermod
   IF (warnflag) THEN
      CALL open_file(iowarn,warlog_file,'OUT','A')
   ENDIF

!
!3. Reset at run time
!--------------------
!

ELSEIF (nt.GT.0) THEN

   IF (nt.EQ.1) THEN
!     ---error checking
      errchk = levprocs_err(iprocglb).GT.1
!     ---close file
      IF (loglev1.GT.0) CALL close_file(iolog,'A')
!     ---set parameters
      loglev1 = levprocs_run(iprocglb)
      loglev2 = MERGE(loglev1,0,exitlog)
   ENDIF

!  ---close/open log file
   IF (loglev1.GT.0.AND.&
    & (mult_index(nt-1,runlog_count).OR.nt.EQ.1)) THEN
      IF (nt.GT.1) CALL close_file(iolog,'A')
      runlog_file = TRIM(runlog_file)
      IF (parallel_set) runlog_file = TRIM(runlog_file)//TRIM(cdloc)
      CALL open_file(iolog,runlog_file,'OUT','A')
   ENDIF

ENDIF


RETURN

END SUBROUTINE monitor_files

!========================================================================

SUBROUTINE open_file(iounit,filename,iotype,ioform)
!************************************************************************
!
! *open_file* Open a file connection
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.7.1
!
! Description -
!
! Module calls - cf90_create, cf90_open, cf90_set_fill, error_abort,
!                error_vals_var_char, get_unit
!
!************************************************************************
!
USE switches
USE cf90_routines, ONLY: cf90_create, cf90_open, cf90_set_fill
USE error_routines, ONLY: error_abort, error_vals_var_char

!
!*Arguments
!
CHARACTER (LEN=lenform), INTENT(IN) :: ioform
CHARACTER (LEN=*), INTENT(IN) :: filename, iotype
INTEGER, INTENT(INOUT) :: iounit

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iounit*    INTEGER File unit
!*filename*  CHAR    Name of file
!*iotype*    CHAR    File access ('IN', 'OUT', 'INOUT')
!*ioform*    CHAR    File format ('A', 'U' or 'N')
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=12) :: cunit
CHARACTER (LEN=20) :: accessspec, actionspec, formspec, positionspec, &
                    & statusspec
INTEGER :: cmode, fillmode, iact, old_mode


pglev = pglev + 1
procname(pglev) = 'open_file'

!
!1. Open specifiers
!------------------
!
!1.1 ACTION, STATUS and POSITION specifiers
!------------------------------------------
!

SELECT CASE (iotype)
CASE ('IN')
   iact = 1; actionspec = 'READ'; statusspec='OLD'; positionspec = 'ASIS'
CASE ('OUT')
   iact = 2; actionspec = 'WRITE'; statusspec='REPLACE'; positionspec = 'ASIS'
CASE ('INOUT')
   iact = 3; actionspec = 'READWRITE'; statusspec='UNKNOWN'
   positionspec = 'APPEND'
CASE DEFAULT
   CALL error_vals_var_char(TRIM(iotype),'iotype','"IN" "OUT" "INOUT"')
END SELECT

!
!1.2 ACCESS and FORM specifiers
!------------------------------
!

SELECT CASE (ioform)
CASE ('A')
   accessspec = 'SEQUENTIAL'; formspec = 'FORMATTED'
CASE ('U')
   accessspec = 'SEQUENTIAL'; formspec = 'UNFORMATTED'
CASE ('N')
CASE DEFAULT
   CALL error_vals_var_char(TRIM(ioform),'ioform','"A" "U" "N"')
END SELECT

!
!2. Error checking
!-----------------
!
!---inquire whether file exists
IF (iact.NE.2) THEN
   IF (TRIM(filename).EQ.'') THEN
      flag = .FALSE.
   ELSE
      INQUIRE (FILE=filename,EXIST=flag)
   ENDIF
   IF (.NOT.flag) THEN
      nerrs = 1
      IF (errchk) THEN
         WRITE (ioerr,'(A)') 'Unable to open non-existing file: '&
                           & //TRIM(filename)
      ENDIF
      CALL error_abort('open_file',ierrno_fopen)
   ENDIF
ENDIF

!---check whether file is not already connected
INQUIRE (FILE=filename,OPENED=flag)
IF (flag) THEN
   INQUIRE (FILE=filename,NUMBER=iounit)
   WRITE (cunit,'(I12)') iounit; cunit = ADJUSTL(cunit)
   nerrs = 1
   IF (errchk) THEN
      WRITE (ioerr,'(A)') 'Unable to open  file: '//TRIM(filename)
      WRITE (ioerr,'(A)') 'File is already connected to unit: '//TRIM(cunit)
   ENDIF
   CALL error_abort('open_file',ierrno_fopen)
ENDIF

!---number of opened files
nopenf = nopenf + 1

!
!3. Unit number
!--------------
!

IF (ioform.NE.'N') THEN
   iounit = get_unit()
ENDIF
WRITE (cunit,'(I12)') iounit; cunit = ADJUSTL(cunit)
     
!
!4. Open file connection
!-----------------------
!

SELECT CASE (ioform)
CASE ('A','U')
   OPEN (UNIT=iounit,FILE=filename,ACTION=actionspec,STATUS=statusspec,&
       & ACCESS=accessspec,FORM=formspec,POSITION=positionspec,ERR=1000)
CASE ('N')
   IF (iact.EQ.1) THEN
      cmode = MERGE(share_NF90,nowrite_NF90,iopt_CDF_shared.EQ.1)
      CALL cf90_open(filename,cmode,iounit)
   ELSEIF (iact.EQ.2) THEN
      cmode = MERGE(share_NF90,clobber_NF90,iopt_CDF_shared.EQ.1)
      cmode = MERGE(cmode,offset_64bit_NF90,iopt_CDF_format.EQ.1)
      CALL cf90_create(filename,cmode,iounit)
   ELSEIF (iact.EQ.3) THEN
      cmode = IOR(write_NF90,share_NF90)
      CALL cf90_open(filename,cmode,iounit)
   ENDIF
   IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

!
!5. Fill mode
!------------
!

IF (ioform.EQ.'N'.AND.iotype.NE.'IN') THEN
   fillmode = MERGE(nofill_NF90,fill_NF90,iopt_CDF_fill.EQ.0)
   CALL cf90_set_fill(iounit,fillmode,old_mode)
ENDIF

!
!6. Write log info
!-----------------
!

IF (loglev1.GT.0) THEN
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Open file '//TRIM(filename)//&
                     & ' on unit '//TRIM(cunit)//' type '//iotype//' &
                     & ('//ioform//')'
ENDIF

pglev = pglev - 1


RETURN

!
!7. Process open error
!---------------------
!

1000 nerrs = 1
IF (errchk) THEN
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Unable to open file '//&
                     & TRIM(filename)//' on unit '//TRIM(cunit)//' type '//&
                     & iotype//' ('//ioform//')'
ENDIF
CALL error_abort(procname(pglev),ierrno_fopen)

END SUBROUTINE open_file

!========================================================================

SUBROUTINE open_filepars(filepars,iotype)
!************************************************************************
!
! *open_filepars* Open a file conneCtion
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.7.1
!
! Description - uses file attributes stored in 'filepars'
!
! Module calls - cf90_create, cf90_open, cf90_set_fill, error_abort,
!                error_vals_var_char, get_unit
!
!************************************************************************
!
USE datatypes
USE switches
USE cf90_routines, ONLY: cf90_create, cf90_open, cf90_set_fill
USE error_routines, ONLY: error_abort, error_vals_var_char

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: iotype
TYPE(FileParams), INTENT(INOUT) :: filepars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*filepars*  DERIVED File attributes
!*iotype*    CHAR    File access ('IN', 'OUT', 'INOUT')
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=1) :: status
CHARACTER (LEN=lenform) :: ioform
CHARACTER (LEN=12) :: cunit
CHARACTER (LEN=20) :: accessspec, actionspec, formspec, positionspec, statusspec
CHARACTER (LEN=leniofile) :: filename
INTEGER :: cmode, fillmode, iact, iounit, old_mode
CHARACTER (LEN=6), DIMENSION(3) :: cstatus = (/'IN   ','OUT  ','INOUT'/)


procname(pglev+1) = 'open_filepars'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

filename = filepars%filename
ioform = filepars%form
status = filepars%status

!
!2. Open specifiers
!------------------
!
!2.1 ACTION, STATUS and POSITION specifiers
!------------------------------------------
!

IF (PRESENT(iotype)) THEN
   SELECT CASE (iotype)
      CASE ('IN')
         iact = 1
         actionspec = 'READ'; statusspec='OLD'; positionspec = 'ASIS'
      CASE ('OUT')
         iact = 2
         actionspec = 'WRITE'; statusspec='REPLACE'; positionspec = 'ASIS'
      CASE ('INOUT')
         iact = 3
         actionspec = 'READWRITE'; statusspec='UNKNOWN'; positionspec = 'APPEND'
      CASE DEFAULT
         CALL error_vals_var_char(TRIM(iotype),'iotype','"IN" "OUT" "INOUT"')
   END SELECT
ELSE
   SELECT CASE (status)
      CASE ('R','N')
         iact = 1
         actionspec = 'READ'; statusspec='OLD'; positionspec = 'ASIS'
      CASE ('W')
         iact = 2
         actionspec = 'WRITE'; statusspec='REPLACE'; positionspec = 'ASIS'
      CASE ('P')
         iact = 3
         actionspec = 'READWRITE'; statusspec='UNKNOWN'; positionspec = 'APPEND'
      CASE ('T')
         iact = 3
         actionspec = 'READWRITE'; statusspec='UNKNOWN'; positionspec = 'ASIS'
      CASE ('0')
         nerrs = 1
         IF (errchk) THEN
            WRITE (ioerr,'(A)') 'Unable to open non-existing file: '&
                              & //TRIM(filename)
         ENDIF
         CALL error_abort('open_filepars',ierrno_fopen)
   END SELECT
ENDIF

!
!2.2 ACCESS and FORM specifiers
!------------------------------
!

SELECT CASE (ioform)
   CASE ('A')
      accessspec = 'SEQUENTIAL'; formspec = 'FORMATTED'
   CASE ('U')
      accessspec = 'SEQUENTIAL'; formspec = 'UNFORMATTED'
   CASE ('N')
   CASE DEFAULT
      CALL error_vals_var_char(TRIM(ioform),'ioform','"A" "U" "N"')
END SELECT

!
!3. Error checking
!-----------------
!
!---inquire whether file exists
IF (iact.NE.2.AND.ioform.NE.'N') THEN
   IF (TRIM(filename).EQ.'') THEN
      flag = .FALSE.
   ELSE
      INQUIRE (FILE=filename,EXIST=flag)
   ENDIF
   IF (.NOT.flag) THEN
      IF (filepars%endfile.EQ.2) THEN
         filepars%iostat = -1
         GOTO 1001
      ELSE
         nerrs = 1
         IF (errchk) THEN
            WRITE (ioerr,'(A)') 'Unable to open non-existing file: '&
                              & //TRIM(filename)
         ENDIF
      ENDIF
      CALL error_abort('open_filepars',ierrno_fopen)
   ENDIF
ENDIF

!---check whether file is not already connected
IF (ioform.NE.'N') THEN
   INQUIRE (FILE=filename,OPENED=flag)
   IF (flag) THEN
      INQUIRE (FILE=filename,NUMBER=iounit)
      WRITE (cunit,'(I12)') iounit; cunit = ADJUSTL(cunit)
      nerrs = 1
      IF (errchk) THEN
         WRITE (ioerr,'(A)') 'Unable to open  file: '//TRIM(filename)
         WRITE (ioerr,'(A)') 'File is already connected to unit: '//TRIM(cunit)
      ENDIF
      CALL error_abort('open_file',ierrno_fopen)
   ENDIF
ENDIF

!---number of opened files
nopenf = nopenf + 1

!
!4. Unit number
!--------------
!

IF (ioform.NE.'N') THEN
   iounit = get_unit()
ENDIF
     
!
!5. Open file connection
!-----------------------
!

SELECT CASE (ioform)
CASE ('A','U')
   OPEN (UNIT=iounit,FILE=filename,ACTION=actionspec,STATUS=statusspec,&
       & ACCESS=accessspec,FORM=formspec,POSITION=positionspec,ERR=1000)
   nerrs = MERGE(1,0,errstat.NE.0)
CASE ('N')
   IF (iact.EQ.1) THEN
      cmode = MERGE(share_NF90,nowrite_NF90,iopt_CDF_shared.EQ.1)
      CALL cf90_open(filename,cmode,iounit)
   ELSEIF (iact.EQ.2) THEN
      cmode = MERGE(share_NF90,clobber_NF90,iopt_CDF_shared.EQ.1)
      cmode = MERGE(cmode,offset_64bit_NF90,iopt_CDF_format.EQ.1)
      CALL cf90_create(filename,cmode,iounit)
   ELSEIF (iact.EQ.3) THEN 
      cmode = IOR(write_NF90,share_NF90)
      CALL cf90_open(filename,cmode,iounit)
   ENDIF
   IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

!
!6. Fill mode
!------------
!

IF (ioform.EQ.'N'.AND.iact.EQ.2) THEN
   fillmode = MERGE(nofill_NF90,fill_NF90,iopt_CDF_fill.EQ.0)
   CALL cf90_set_fill(iounit,fillmode,old_mode)
ENDIF

!
!7. Write log info
!-----------------
!

WRITE (cunit,'(I12)') iounit; cunit = ADJUSTL(cunit)
IF (loglev1.GT.0) THEN
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Open file '//TRIM(filename)//&
                     & ' on unit '//TRIM(cunit)//' type '//TRIM(cstatus(iact))&
                     & //' ('//ioform//')'
ENDIF

!
!8. Update parameters in 'filepars'
!----------------------------------
!

filepars%iunit = iounit
filepars%iostat = 1
filepars%timerec = 0

1001 CALL log_timer_out()


RETURN

!
!9. Process open error
!---------------------
!

1000 nerrs = 1
WRITE (cunit,'(I12)') iounit; cunit = ADJUSTL(cunit)
IF (errchk) THEN
   WRITE (iolog,'(A)') REPEAT(' ',pglev-1)//'Unable to open file '//&
                     & TRIM(filename)//' on unit '//TRIM(cunit)//' type '//&
                     & TRIM(cstatus(iact))//' ('//ioform//')'
ENDIF
CALL error_abort(procname(pglev),ierrno_fopen)

END SUBROUTINE open_filepars

!========================================================================

FUNCTION output_flag(ciodatetime,tskips)
!************************************************************************
!
! *output_flag* Decode whether output is written for a given date/time  
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.0
!
! Description - result depends on value of the tskips (1-9)
!
! Module calls - error_abort
!
!************************************************************************
!
USE timepars

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: ciodatetime
INTEGER, INTENT(IN) :: tskips
LOGICAL :: output_flag

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR    Date/time within data file
!*tskips*      INTEGER Selects criterion for output permission
!
!------------------------------------------------------------------------------
!

SELECT CASE (tskips)
   CASE (1); output_flag = .TRUE.
   CASE (2); output_flag = ciodatetime.noearlier.CStartDateTime
   CASE (3); output_flag = ciodatetime.later.CStartDateTime
   CASE (4); output_flag = ciodatetime.nolater.CEndDatetime
   CASE (5)
      output_flag = (ciodatetime.noearlier.CStartDateTime).AND.&
                  & (ciodatetime.nolater.CEndDatetime)
   CASE (6)
      output_flag = (ciodatetime.later.CStartDateTime).AND.&
                  & (ciodatetime.nolater.CEndDatetime)
   CASE (7); output_flag = ciodatetime.earlier.CEndDatetime
   CASE (8)
      output_flag = (ciodatetime.noearlier.CStartDateTime).AND.&
                  & (ciodatetime.earlier.CEndDatetime)
   CASE (9)
      output_flag = (ciodatetime.later.CStartDateTime).AND.&
                  & (ciodatetime.earlier.CEndDatetime)
END SELECT


RETURN

END FUNCTION output_flag

!========================================================================

SUBROUTINE read_glbatts_mod(filepars)
!************************************************************************
!
! *read_glbatts_mod* Read the global atrributes of a model file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - error_file, inquire_modfil_dims, inquire_unit
!
!************************************************************************
!
USE datatypes
USE error_routines, ONLY: error_file

!
!*Arguments
!
TYPE (FileParams), INTENT(INOUT) :: filepars

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*filepars* DERIVED  File attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: inform
INTEGER :: iunit, npcc


procname(pglev+1) = 'read_glbatts_mod'
CALL log_timer_in(npcc,&
                & logname=TRIM(procname(pglev+1))//': '//filepars%filename)

!
!1. File parameters
!------------------
!

iunit = filepars%iunit
inform = filepars%form

!
!2. Check whether file unit is opened
!------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!3. Global attributes
!--------------------
!

CALL inquire_modfil_dims(iunit,inform,ncoordinates=filepars%nocoords,&
                       & nvariables=filepars%novars,&
                       & maxrecs=filepars%maxrecs,&
                       & floattype=filepars%floattype,&
                       & packing=filepars%packing,&
                       & timeid=filepars%timeid)

CALL log_timer_out(npcc,itm_input)


RETURN

END SUBROUTINE read_glbatts_mod

!========================================================================

SUBROUTINE read_metadata_line(iunit,ioform,cvals,numvars,cname,linenum,istat)
!************************************************************************
!
! *read_metadata_line* Parse a series of variables, separated by a delimiter as
!                      character data on an output metadata line in a standard
!                      output file 
!
! Author - SORESMA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.7.1
!
! Description - 
!
! Module calls -
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=lenform), INTENT(IN) :: ioform 
CHARACTER (LEN=lenname), INTENT(OUT), OPTIONAL :: cname
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(:) :: cvals
INTEGER, INTENT(IN) :: iunit
INTEGER, INTENT(OUT) :: istat, numvars
INTEGER, INTENT(INOUT) :: linenum

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iunit*     INTEGER File unit number
!*ioform*    CHAR    File format
!*cvals*     CHAR    Data values in string format
!*numvars*   INTEGER Number of variables on the input line
!*cname*     CHAR    Name of the variable(s)
!*linenum*   INTEGER Line number
!*istat*     INTEGER Read status (2 at end of read block, 1 otherwise)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=lencifline) :: cline, xline
INTEGER :: iend, ieq, isep, ival, l


!
!1. Initialise
!-------------
!

cvals = ''; numvars = 0

!
!2. Read input line
!------------------
!

IF (ioform.EQ.'A') THEN
   READ (iunit,'(A)') cline
ELSE
   READ (iunit) cline
ENDIF
linenum = linenum + 1

IF (cline(1:1).EQ.'#') THEN
   istat = 2
ELSE
   istat = 1
ENDIF

!
!3. Removing comments from input line
!------------------------------------
!

cline = ADJUSTL(cline)
iend = INDEX(cline,'!')
IF (iend.EQ.1) THEN
   cline = ''; cvals = ''; numvars = 0
   cname = ''
   RETURN
ENDIF
IF (iend.GT.0) cline = ADJUSTL(cline(1:iend-1))

!
!4. Data variable name
!---------------------
!

ieq = INDEX(cline,':')
cname = TRIM(ADJUSTL(cline(1:ieq-1)))

!
!5. Separate data
!----------------
!
!---first value
ival = 1
isep = INDEX(cline,',')
IF (isep.EQ.ieq+1) THEN
   cvals(1) = ''
ELSEIF (isep.EQ.0) THEN
   IF (ieq.EQ.LEN(cline)) THEN
      cvals(1) = ''
   ELSE
      cvals(1) = cline(ieq+1:LEN(cline))
   ENDIF
ELSE
   cvals(1) = cline(ieq+1:isep-1)
ENDIF
IF (isep.EQ.LEN(cline)) THEN
   xline = ''
ELSE
   xline = cline(isep+1:LEN(cline))
ENDIF

!---next value
DO WHILE (isep.GT.0)
   ival = ival + 1
   l = INDEX(xline,',')
   isep = MERGE(0,l,xline.EQ.'')
   IF (isep.EQ.0) THEN
      cvals(ival) = TRIM(xline)
   ELSEIF (isep.EQ.1) THEN
      cvals(ival) = ''
   ELSE
      cvals(ival) = xline(1:isep-1)
   ENDIF
   IF (isep.EQ.LEN(xline)) THEN
      xline = ''
   ELSE
      xline = xline(isep+1:LEN(xline))
   ENDIF
ENDDO

!---number of data values
numvars = ival

!---removing trailing blanks
cvals = ADJUSTL(cvals)


RETURN

END SUBROUTINE read_metadata_line

!========================================================================

SUBROUTINE read_output_metadata(filepars)
!************************************************************************
!
! *read_output_metadata* read file header and coordinate arrays without storing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.7.1
!
! Description - 
!
! Module calls - error_alloc, error_alloc_struct, inquire_modfil_atts,
!                inquire_modfil_dims, read_vars, varatts_init
!
!************************************************************************
!
USE datatypes
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_alloc, error_alloc_struc

!
!*Arguments
!
TYPE (FileParams), INTENT(INOUT) :: filepars

!
!*Local variables
!
LOGICAL :: gridded, hrect, packing
CHARACTER (LEN=lenform) :: outform
INTEGER :: i, icoord, idim, iunit, j, l, n, ncout, nocoords, nodim, nodims, &
         & nrout, nostats, nowetout, nzout, varid
LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: maskout
CHARACTER (LEN=lenname), SAVE, ALLOCATABLE, DIMENSION(:) :: dimnames
INTEGER, DIMENSION(3) :: ndims
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: lmaskout, outdims
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: gsigout
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: depout, xout, yout
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: coordvars


procname(pglev+1) = 'read_output_metadata'
CALL log_timer_in()

!
!1. Number of variables, dimensions
!----------------------------------
!

iunit = filepars%iunit
outform = filepars%form
IF (outform.EQ.'N') GOTO 1001
CALL inquire_modfil_dims(iunit,outform,ndimensions=nodims,&
                       & ncoordinates=nocoords,floattype=filepars%floattype)

!
!2. Allocate
!-----------
!

ALLOCATE (outdims(nodims),STAT=errstat)
CALL error_alloc('outdims',1,(/nodims/),kndint)
outdims = 0

ALLOCATE (dimnames(nodims),STAT=errstat)
CALL error_alloc('dimnames',1,(/nodims/),kndchar)
dimnames = ''

!
!3. Array dimensions
!-------------------
!

CALL inquire_modfil_dims(iunit,outform,dimensions=outdims,dimnames=dimnames)

ncout = 0; nrout = 0; nowetout = 0; nostats = 0; nzout = 0
gridded = .TRUE.; packing = .FALSE. 

idim_310: DO idim=1,nodims
   SELECT CASE (TRIM(dimnames(idim)))
      CASE ('x','lon')
         ncout = outdims(idim)
         hrect = .TRUE.
      CASE ('y','lat')
         nrout = outdims(idim)
      CASE ('nc')
         ncout = outdims(idim)
         hrect = .FALSE.
      CASE ('nr')
         nrout = outdims(idim)
      CASE ('seapoint')
         nowetout = outdims(idim)
         packing = .TRUE.
      CASE ('station')
         nostats = outdims(idim)
         gridded = .FALSE.
      CASE ('lev')
         nzout = outdims(idim)
   END SELECT
ENDDO idim_310

!
!4. File dimension
!-----------------
!

IF (nodims.EQ.1) THEN
   nodim = 0
   DEALLOCATE (dimnames,outdims)
   GOTO 1001
ELSE
   nodim = MERGE(2,3,nzout.EQ.0)
ENDIF

!
!5. Coordinate attributes
!------------------------
!

ALLOCATE (coordvars(nocoords),STAT=errstat)
CALL error_alloc_struc('coordvars',1,(/nocoords/),'VariableAtts')
CALL varatts_init(coordvars)

icoord_510: DO icoord=1,nocoords
   CALL inquire_modfil_atts(iunit,outform,icoord,&
                          & nrank=coordvars(icoord)%nrank,&
                          & data_type=coordvars(icoord)%data_type,&        
                          & global_dims=coordvars(icoord)%global_dims,&
                          & f90_name=coordvars(icoord)%f90_name,&
                          & standard_name=coordvars(icoord)%standard_name,&
                          & long_name=coordvars(icoord)%long_name,&
                          & units=coordvars(icoord)%units,&
                          & fill=coordvars(icoord)%fill,&
                          & fill_value=coordvars(icoord)%fill_value)
ENDDO icoord_510

!
!6. Allocate coordinate arrays
!-----------------------------
!

IF (gridded) THEN
   ndims = (/ncout,nrout,nzout/)
ELSE
   ndims = (/nostats,1,nzout/)
ENDIF
ALLOCATE (xout(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('xout',2,ndims(1:2),kndrtype)
ALLOCATE (yout(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('yout',2,ndims(1:2),kndrtype)
ALLOCATE (gsigout(ndims(3)),STAT=errstat)
CALL error_alloc('gsigout',1,ndims(3),kndrtype)
IF (nzout.GT.0) gsigout = 0.0
ALLOCATE (depout(ndims(1),ndims(2)),STAT=errstat)
CALL error_alloc('depout',2,ndims(1:2),kndrtype)
IF (packing) THEN
   ALLOCATE (lmaskout(nowetout),STAT=errstat)
   CALL error_alloc('lmaskout',1,(/nowetout/),kndint)
   ALLOCATE (maskout(ndims(1),ndims(2)),STAT=errstat)
   CALL error_alloc('maskout',2,ndims(1:2),kndlog)
   maskout = .FALSE.
ENDIF

!
!7. Skip coordinate arrays
!-------------------------
!

varid = 1

IF (gridded) THEN

!  ---horizontal coordinates (regular grid)
   IF (hrect) THEN
      varid = varid + 1
      CALL read_vars(xout(:,1),filepars,varid,coordvars)
      varid = varid + 1
      CALL read_vars(yout(1,:),filepars,varid,coordvars)
!  ---horizontal coordinates (curvilinear grid)
   ELSE
      varid = varid + 1
      CALL read_vars(xout,filepars,varid,coordvars)
      varid = varid + 1
      CALL read_vars(yout,filepars,varid,coordvars)
   ENDIF
!  ---land mask array
   IF (packing) THEN
      varid = varid + 1
      CALL read_vars(lmaskout,filepars,varid,coordvars)
      l_710: DO l=1,nowetout
         n = lmaskout(l)
         j = (n-1)/ncout + 1
         i = n - (j-1)*ncout
         maskout(i,j) = .TRUE.
      ENDDO l_710
      varid = varid + 1
      CALL read_submod(depout,filepars,varid,varatts=coordvars,maskvals=maskout)
   ENDIF
!  ---bathymetry
   IF (.NOT.packing) THEN
      varid = varid + 1
      CALL read_vars(depout,filepars,varid,coordvars)
   ENDIF

ELSE

!  ---horizontal coordinates and bathymetry (irregular grid)
   varid = varid + 1
   CALL read_vars(xout(:,1),filepars,varid,coordvars)
   varid = varid + 1
   CALL read_vars(yout(:,1),filepars,varid,coordvars)
   varid = varid + 1
   CALL read_vars(depout(:,1),filepars,varid,coordvars)
ENDIF

!---sigma/s-coordinates
IF (nzout.GT.0) THEN
   varid = varid + 1
   CALL read_vars(gsigout,filepars,varid,coordvars)
ENDIF

!
!8. Deallocate
!-------------
!

DEALLOCATE (coordvars,dimnames,outdims)
DEALLOCATE(depout,gsigout,xout,yout)
IF (packing) DEALLOCATE (lmaskout,maskout)

1001 CALL log_timer_out()

END SUBROUTINE read_output_metadata

!========================================================================

SUBROUTINE read_station_names(filepars,station_names,nostats)
!************************************************************************
!
! *read_station_names* Read the names of stations in user output file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.7.1
!
! Description -
!
! Module calls - cf90_get_att_chars, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_att_chars
USE error_routines, ONLY: nerrs, error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: nostats
TYPE (FileParams), INTENT(INOUT) :: filepars
CHARACTER (LEN=lenname), INTENT(OUT), DIMENSION(nostats) :: station_names

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*filepars*      DERIVED  Data file attributes
!*station_names* CHAR     Station names
!*nostats*       INTEGER  Number of stations 
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: inform
INTEGER :: istat, iunit, l, npcc


procname(pglev+1) = 'read_station_names'
CALL log_timer_in(npcc,&
                & logname=TRIM(procname(pglev+1))//': '//filepars%filename)

!
!1. File parameters
!------------------
!

iunit = filepars%iunit
inform = filepars%form

!
!2 Check whether file unit is opened
!-----------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!3. Read station names
!-------------------------
!

SELECT CASE (inform)
   CASE ('A')
      istat_310: DO istat=1,nostats
         READ (iunit,'(I5,1X,A)') l, station_names(istat)
      ENDDO istat_310
   CASE ('U')
      istat_320: DO istat=1,nostats
         READ (iunit) l, station_names(istat)
      ENDDO istat_320
   CASE ('N') 
      CALL cf90_get_att_chars(iunit,global_NF90,'station_names',lenname,&
                            & station_names)
END SELECT

IF (nerrs.GT.0) CALL error_file(ierrno_read,filepars=filepars)

CALL log_timer_out(npcc,itm_input)


RETURN

END SUBROUTINE read_station_names

!========================================================================

SUBROUTINE read_submod_real_2d(realsub,filepars,varid,varatts,maskvals)
!************************************************************************
!
! *read_submod_real_2d* Read a global real sub-model 2-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.10.1
!
! Description - only called in serial mode
!
! Module calls - error_alloc, read_vars
!
!************************************************************************
!
USE datatypes
USE error_routines, ONLY: error_alloc

!
!*  Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE(FileParams), INTENT(IN) :: filepars
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: maskvals
REAL, INTENT(INOUT), DIMENSION(:,:) :: realsub
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realsub*   REAL     Sub-model array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*varatts*   DERIVED  Variable attributes
!*maskvals*  LOGICAL  Array of mask values used in case a land mask is applied
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: ivarid, nowet
REAL :: fill_value
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: realdat


IF (SIZE(realsub).EQ.0) RETURN

procname(pglev+1) = 'read_submod_real_2d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Initialise paramaters and arrays
!------------------------------------
!
!---fill parameters
IF (filepars%fill) THEN
   IF (PRESENT(varatts)) THEN
      fill_value = varatts(varid)%fill_value
   ELSE
      fill_value = double_fill
   ENDIF
ELSE
   fill_value = 0.0
ENDIF

!---packing
packing = filepars%packing

!---allocate array for packing
IF (packing) THEN
   nowet = COUNT(maskvals)
   ALLOCATE (realdat(nowet),STAT=errstat)
   CALL error_alloc('realdat',1,(/nowet/),kndrtype)
ENDIF

!
!2. Read data
!------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realdat,filepars,varid,varatts)
   ELSE
      CALL read_vars(realdat,filepars,varid)
   ENDIF
ELSE
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realsub,filepars,varid,varatts)
   ELSE
      CALL read_vars(realsub,filepars,varid)
   ENDIF
ENDIF

!
!3. Unpack data using land mask
!------------------------------
!

IF (packing) THEN
   realsub = UNPACK(realdat,MASK=maskvals,FIELD=fill_value)
ENDIF

!
!4. Deallocate array
!-------------------
!

IF (packing) DEALLOCATE (realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE read_submod_real_2d

!========================================================================

SUBROUTINE read_submod_real_3d(realsub,filepars,varid,varatts,vecids,maskvals)
!************************************************************************
!
! *read_submod_real_3d* Read a global real sub-model 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.10.1
!
! Description -
!
! Module calls - error_alloc, read_vars
!
!************************************************************************
!
USE datatypes
USE error_routines, ONLY: error_alloc

!
!*  Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE(FileParams), INTENT(IN) :: filepars
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: maskvals
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL, INTENT(INOUT), DIMENSION(:,:,:) :: realsub
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realsub*   REAL     Sub-model array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*   DERIVED  Variable attributes
!*vecids*    INTEGER  Arrays of variable IDs if varid = 0
!*maskvals*  LOGICAL  Array of mask values used in case a land mask is applied
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: ivar, ivarid, k, nowet, n3dim
REAL :: fill_value
INTEGER, DIMENSION(SIZE(realsub,DIM=3)) :: idvars
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realdat


IF (SIZE(realsub).EQ.0) RETURN

procname(pglev+1) = 'read_submod_real_3d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

n3dim = SIZE(realsub,DIM=3)
IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,n3dim)/)
ENDIF

!
!2. Initialise paramaters and arrays
!------------------------------------
!
!---fill parameters
IF (filepars%fill) THEN
   IF (PRESENT(varatts)) THEN
      fill_value = varatts(varid)%fill_value
   ELSE
      fill_value = double_fill
   ENDIF
ELSE
   fill_value = 0.0
ENDIF

!---packing
packing = filepars%packing

!---allocate array for packing
IF (packing) THEN
   nowet = COUNT(maskvals)
   ALLOCATE (realdat(nowet,n3dim),STAT=errstat)
   CALL error_alloc('realdat',2,(/nowet,n3dim/),kndrtype)
ENDIF

!
!3. Read data
!------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realdat,filepars,varid,varatts,idvars)
   ELSE
      CALL read_vars(realdat,filepars,varid,vecids=idvars)
   ENDIF
ELSE
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realsub,filepars,varid,varatts,idvars)
   ELSE
      CALL read_vars(realsub,filepars,varid,vecids=idvars)
   ENDIF
ENDIF

!
!4. Unpack data using land mask
!------------------------------
!

IF (packing) THEN
   k_410: DO k=1,n3dim
      realsub(:,:,k) = UNPACK(realdat(:,k),MASK=maskvals,FIELD=fill_value)
   ENDDO k_410
ENDIF

!
!5. Deallocate array
!-------------------
!

IF (packing) DEALLOCATE (realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE read_submod_real_3d

!========================================================================

SUBROUTINE read_submod_real_4d(realsub,filepars,varid,varatts,vecids,maskvals)
!************************************************************************
!
! *read_submod_real_4d* Read a global real sub-model 4-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.10.1
!
! Description -
!
! Module calls - error_alloc, read_vars
!
!************************************************************************
!
USE datatypes
USE error_routines, ONLY: error_alloc

!
!*  Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE(FileParams), INTENT(IN) :: filepars
LOGICAL, INTENT(IN), OPTIONAL, DIMENSION(:,:) :: maskvals
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL, INTENT(INOUT), DIMENSION(:,:,:,:) :: realsub
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realsub*   REAL     Sub-model array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*varatts*   DERIVED  Variable attributes
!*vecids*    INTEGER  Arrays of variable IDs if varid = 0
!*maskvals*  LOGICAL  Array of mask values used in case a land mask is applied
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: packing
INTEGER :: ivar, ivarid, k, l, nowet, n3dim, n4dim
REAL :: fill_value
INTEGER, DIMENSION(SIZE(realsub,DIM=4)) :: idvars
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realdat


IF (SIZE(realsub).EQ.0) RETURN

procname(pglev+1) = 'read_submod_real_4d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Optional arguments
!---------------------
!

n3dim = SIZE(realsub,DIM=3); n4dim = SIZE(realsub,DIM=4)
IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,n4dim)/)
ENDIF

!
!2. Initialise paramaters and arrays
!------------------------------------
!
!---fill parameters
IF (filepars%fill) THEN
   IF (PRESENT(varatts)) THEN
      fill_value = varatts(varid)%fill_value
   ELSE
      fill_value = double_fill
   ENDIF
ELSE
   fill_value = 0.0
ENDIF

!---packing
packing = filepars%packing

!---allocate array for packing
IF (packing) THEN
   nowet = COUNT(maskvals)
   ALLOCATE (realdat(nowet,n3dim,n4dim),STAT=errstat)
   CALL error_alloc('realdat',3,(/nowet,n3dim,n4dim/),kndrtype)
ENDIF

!
!3. Read data
!------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realdat,filepars,varid,varatts,idvars)
   ELSE
      CALL read_vars(realdat,filepars,varid,vecids=idvars)
   ENDIF
ELSE
   IF (PRESENT(varatts)) THEN
      CALL read_vars(realsub,filepars,varid,varatts,idvars)
   ELSE
      CALL read_vars(realsub,filepars,varid,vecids=idvars)
   ENDIF
ENDIF

!
!4. Unpack data using land mask
!------------------------------
!

IF (packing) THEN
   l_410: DO l=1,n4dim
   k_410: DO k=1,n3dim
      realsub(:,:,k,l) = UNPACK(realdat(:,k,l),MASK=maskvals,FIELD=fill_value)
   ENDDO k_410
   ENDDO l_410
ENDIF

!
!5. Deallocate array
!-------------------
!

IF (packing) DEALLOCATE (realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE read_submod_real_4d

!========================================================================

SUBROUTINE read_time_char(ciodatetime,filepars,timerec,back)
!************************************************************************
!
! *read_time_char* Read a date/time in string format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var_chars, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var_chars
USE error_routines, ONLY: error_file

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: back
CHARACTER (LEN=lentime), INTENT(OUT)  :: ciodatetime
INTEGER, INTENT(IN), OPTIONAL :: timerec
TYPE(FileParams), INTENT(INOUT) :: filepars

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR    Date/time
!*filepars*    DERIVED File attributes
!*timerec*     INTEGER Time record number
!*back*        LOGICAL If present and .TRUE., reading is performed backwards
!                      (i.e. starting with the last time record)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, maxrecs, npcc, nsign, rstat, timerecx


procname(pglev+1) = 'read_time_char'
CALL log_timer_in(npcc,&
                & logname=TRIM(procname(pglev+1))//': '//filepars%filename)

!
!1. Check and initialise
!-----------------------
!
!---parameters for reading
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
maxrecs = filepars%maxrecs
IF (PRESENT(timerec)) THEN
   timerecx = timerec
ELSE
   IF (PRESENT(back)) THEN
      nsign = MERGE (-1,1,back)
   ELSE
      nsign = 1
   ENDIF
   timerecx = filepars%timerec + nsign
   filepars%timerec = timerecx
ENDIF

!---check whether file unit is opened
flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!

SELECT CASE (inform)
   CASE ('A')
      READ (iunit,'(A)',IOSTAT=errstat) ciodatetime
      rstat = MERGE(1,3,errstat.GE.0)
      errstat = MAX(0,errstat)
      IF (errstat.NE.0) GOTO 1000
   CASE ('U')
      READ (iunit,IOSTAT=errstat) ciodatetime
      rstat = MERGE(1,3,errstat.GE.0)
      errstat = MAX(0,errstat)
      IF (errstat.NE.0) GOTO 1000
   CASE ('N')
      IF (timerecx.LE.0.OR.timerecx.GT.maxrecs) THEN
         rstat = 2
      ELSE
         rstat = 1
         CALL cf90_get_var_chars(iunit,filepars%timeid,ciodatetime,timerecx)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      ENDIF
END SELECT

!
!3. Reset file parameters
!------------------------
!

filepars%iostat = rstat

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read time'
CALL error_file(ierrno_read,filepars=filepars)

END SUBROUTINE read_time_char

!========================================================================

SUBROUTINE read_time_real(realtime,filepars,timerec,back)
!************************************************************************
!
! *read_time_real* Read a real date/time
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_file

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: back
INTEGER, INTENT(IN), OPTIONAL :: timerec
REAL (KIND=kndrlong), INTENT(OUT)  :: realtime
TYPE(FileParams), INTENT(INOUT) :: filepars

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*realtime*    REAL    Date/time
!*filepars*    DERIVED File attributes
!*timerec*     INTEGER Time record number
!*back*        LOGICAL If present and .TRUE., reading is performed backwards
!                      (i.e. starting with the last time record)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, maxrecs, npcc, nsign, timerecx, rstat


procname(pglev+1) = 'read_time_real'
CALL log_timer_in(npcc,&
                & logname=TRIM(procname(pglev+1))//': '//filepars%filename)

!
!1. Check and initialise
!-----------------------
!
!---parameters for reading
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
maxrecs = filepars%maxrecs
IF (PRESENT(timerec)) THEN
   timerecx = timerec
ELSE
   IF (PRESENT(back)) THEN
      nsign = MERGE (-1,1,back)
   ELSE
      nsign = 1
   ENDIF
   timerecx = filepars%timerec + nsign
   filepars%timerec = timerecx
ENDIF

!---check whether file unit is opened
flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!

SELECT CASE (inform)
   CASE ('A')
      READ (iunit,DoubleFormat,IOSTAT=errstat) realtime
      rstat = MERGE(1,3,errstat.GE.0)
      errstat = MAX(0,errstat)
      IF (errstat.NE.0) GOTO 1000
   CASE ('U')
      READ (iunit,IOSTAT=errstat) realtime
      rstat = MERGE(1,3,errstat.GE.0)
      errstat = MAX(0,errstat)
      IF (errstat.NE.0) GOTO 1000
   CASE ('N')
      IF (timerecx.LE.0.OR.timerecx.GT.maxrecs) THEN
         rstat = 2
      ELSE
         rstat = 1
         CALL cf90_get_var(iunit,filepars%timeid,realtime,timerecx)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      ENDIF
END SELECT

!
!3. Reset file parameters
!------------------------
!

filepars%iostat = rstat

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read time'
CALL error_file(ierrno_read,filepars=filepars)

END SUBROUTINE read_time_real

!========================================================================

SUBROUTINE read_varatts_mod(filepars,varatts,numvars)
!************************************************************************
!
! *read_varatts_mod* Read the variable atrributes of a model file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - error_file, inquire_modfil_atts, inquire_unit, write_info_mod
!
!************************************************************************
!
USE datatypes
USE error_routines, ONLY: error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: numvars
TYPE (FileParams), INTENT(INOUT) :: filepars
TYPE (VariableAtts), INTENT(OUT), DIMENSION(numvars) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*filepars* DERIVED  File attributes
!*varatts*  DERIVED  Variable attributes
!*numvars*  INTEGER  Number of variables in file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: inform
INTEGER :: icoord, iunit, ivar, nocoords, npcc


procname(pglev+1) = 'read_varatts_mod'
CALL log_timer_in(npcc,&
                & logname=TRIM(procname(pglev+1))//': '//filepars%filename)

!
!1. File parameters
!------------------
!

iunit = filepars%iunit
inform = filepars%form
nocoords = filepars%nocoords

!
!2. Check whether file unit is opened
!------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!3. Read variable attributes
!---------------------------
!

IF (nocoords.EQ.1) filepars%timeid = 1

!---coordinate attributes
icoord_310: DO icoord=1,nocoords
   CALL inquire_modfil_atts(iunit,inform,icoord,&
                          & nrank=varatts(icoord)%nrank,&
                          & data_type=varatts(icoord)%data_type,&
                          & global_dims=varatts(icoord)%global_dims,&
                          & f90_name=varatts(icoord)%f90_name,&
                          & long_name=varatts(icoord)%long_name,&
                          & units=varatts(icoord)%units,&
                          & fill=varatts(icoord)%fill,&
                          & fill_value=varatts(icoord)%fill_value)
ENDDO icoord_310

!---coordinate attributes
ivar_320: DO ivar=nocoords+1,numvars
   CALL inquire_modfil_atts(iunit,inform,ivar,&
                          & nrank=varatts(ivar)%nrank,&
                          & data_type=varatts(ivar)%data_type,&
                          & global_dims=varatts(ivar)%global_dims,&
                          & f90_name=varatts(ivar)%f90_name,&
                          & long_name=varatts(ivar)%long_name,&
                          & units=varatts(ivar)%units,&
                          & fill=varatts(ivar)%fill,&
                          & fill_value=varatts(ivar)%fill_value)
ENDDO ivar_320

!
!4. Write information file
!-------------------------
!

IF (filepars%info) THEN
   CALL write_info_mod(filepars,varatts,numvars)
ENDIF

IF (nerrs.GT.0) CALL error_file(ierrno_read,filepars=filepars)

CALL log_timer_out(npcc,itm_input)


RETURN

END SUBROUTINE read_varatts_mod

!========================================================================

SUBROUTINE read_vars_char_0d(chardat,filepars,varid,varatts)
!************************************************************************
!
! *read_vars_char_0d* Read a scalar character scalar variable in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var_chars, error_abort, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var_chars
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
CHARACTER (LEN=*), INTENT(OUT) :: chardat
TYPE (FileParams), INTENT(IN) :: filepars
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*chardat*  CHAR     Character input variable
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  variable id
!*varatts*  DERIVED  Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=1) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, novars, npcc, timerec


procname(pglev+1) = 'read_vars_char_0d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!
!1.1 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.2 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!

SELECT CASE (inform)
CASE ('A')
   READ (iunit,*)
   READ (iunit,'(A)',IOSTAT=errstat,ERR=1000) chardat
CASE ('U')
   READ (iunit,IOSTAT=errstat,ERR=1000) chardat
CASE ('N')
   CALL cf90_get_var_chars(iunit,varid,chardat,timerec)
   IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
nerrs = 1
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read character variable'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_char_0d

!========================================================================

SUBROUTINE read_vars_char_1d(chardat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_char_1d* Read a 1-D character array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var_chars, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var_chars
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(:) :: chardat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*chardat*  CHAR     Character input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: allrec, flag
CHARACTER (LEN=1) :: inform
CHARACTER (LEN=12) :: clen, csize
CHARACTER (LEN=29) :: cfmt
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(chardat)) :: idvars


IF (SIZE(chardat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_char_1d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(chardat)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---type of input record
allrec = inform.NE.'N'.OR.varid.NE.0

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (allrec) THEN

   SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      WRITE (csize,'(I12)') nosize; csize = ADJUSTL(csize)
      WRITE (clen,'(I12)') LEN(chardat(1)); clen = ADJUSTL(clen)
      cfmt = '('//TRIM(csize)//'A'//TRIM(clen)//')'
      READ (iunit,TRIM(cfmt),IOSTAT=errstat,ERR=1000) chardat
   CASE ('U')
      READ (iunit,IOSTAT=errstat,ERR=1000) chardat
   CASE ('N')
      CALL cf90_get_var_chars(iunit,varid,chardat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      CALL cf90_get_var_chars(iunit,idvars(l),chardat(l),timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read character array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_char_1d

!========================================================================

SUBROUTINE read_vars_char_2d(chardat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_char_2d* Read a 2-D character array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var_chars, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var_chars
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(:,:) :: chardat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*chardat*  CHAR     Character input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=1) :: inform
CHARACTER (LEN=12) :: clen, csize
CHARACTER (LEN=29) :: cfmt
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(chardat,DIM=2)) :: idvars


IF (SIZE(chardat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_char_2d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(chardat,DIM=2)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      WRITE (csize,'(I12)') nosize; csize = ADJUSTL(csize)
      WRITE (clen,'(I12)') LEN(chardat(1,1)); clen = ADJUSTL(clen)
      cfmt = '('//TRIM(csize)//'A'//TRIM(clen)//')'
      READ (iunit,TRIM(cfmt),IOSTAT=errstat,ERR=1000) chardat
   CASE ('U')
      READ (iunit,IOSTAT=errstat,ERR=1000) chardat
   CASE ('N')
      CALL cf90_get_var_chars(iunit,varid,chardat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension   
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      SELECT CASE (inform)
      CASE ('A')
         READ (iunit,*)
         WRITE (csize,'(I12)') nosize; csize = ADJUSTL(csize)
         WRITE (clen,'(I12)') LEN(chardat(1,1)); clen = ADJUSTL(clen)
         cfmt = '('//TRIM(csize)//'A'//TRIM(clen)//')'
         READ (iunit,TRIM(cfmt),IOSTAT=errstat,ERR=1000) chardat(:,l)
      CASE ('U')
         READ (iunit,IOSTAT=errstat,ERR=1000) chardat(:,l)
      CASE ('N')
         CALL cf90_get_var_chars(iunit,idvars(l),chardat(:,l),timerec)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      END SELECT
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read character array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_char_2d

!========================================================================

SUBROUTINE read_vars_int_0d(intdat,filepars,varid,varatts)
!************************************************************************
!
! *read_vars_int_0d* Read a scalar integer scalar variable in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(OUT) :: intdat
TYPE (FileParams), INTENT(IN) :: filepars
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER  Integer input variable
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  variable id
!*varatts*  DERIVED  Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, novars, npcc, timerec


procname(pglev+1) = 'read_vars_int_0d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!
!1.1 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.2 Check whether file unit is opened
!-------------------------------------
!
flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!

SELECT CASE (inform)
CASE ('A')
   READ (iunit,*)
   READ (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat
CASE ('U')
   READ (iunit,IOSTAT=errstat,ERR=1000) intdat
CASE ('N')
   CALL cf90_get_var(iunit,varid,intdat,timerec)
   IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
nerrs = 1
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read integer variable'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_int_0d

!========================================================================

SUBROUTINE read_vars_int_1d(intdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_int_1d* Read a 1-D integer array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(OUT), DIMENSION(:) :: intdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER  Integer input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: allrec, flag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(intdat)) :: idvars


IF (SIZE(intdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_int_1d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(intdat)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---type of input record
allrec = inform.NE.'N'.OR.varid.NE.0

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (allrec) THEN

   SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      READ (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat
   CASE ('U')
      READ (iunit,IOSTAT=errstat,ERR=1000) intdat
   CASE ('N')
      CALL cf90_get_var(iunit,varid,intdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      CALL cf90_get_var(iunit,idvars(l),intdat(l),timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read integer array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_int_1d

!========================================================================

SUBROUTINE read_vars_int_2d(intdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_int_2d* Read a 2-D integer array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(OUT), DIMENSION(:,:) :: intdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER  Integer input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(intdat,DIM=2)) :: idvars


IF (SIZE(intdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_int_2d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(intdat,DIM=2)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!

iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      READ (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat
   CASE ('U')
      READ (iunit,IOSTAT=errstat,ERR=1000) intdat
   CASE ('N')
      CALL cf90_get_var(iunit,varid,intdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      SELECT CASE (inform)
      CASE ('A')
         READ (iunit,*)
         READ (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat(:,l)
      CASE ('U')
         READ (iunit,IOSTAT=errstat,ERR=1000) intdat(:,l)
      CASE ('N')
         CALL cf90_get_var(iunit,idvars(l),intdat(:,l),timerec)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      END SELECT
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read integer array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_int_2d

!========================================================================

SUBROUTINE read_vars_int_3d(intdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_int_3d* Read a 3-D integer array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(OUT), DIMENSION(:,:,:) :: intdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER  Integer input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(intdat,DIM=3)) :: idvars


IF (SIZE(intdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_int_3d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(intdat,DIM=3)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!

iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      READ (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat
   CASE ('U')
      READ (iunit,IOSTAT=errstat,ERR=1000) intdat
   CASE ('N')
      CALL cf90_get_var(iunit,varid,intdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      SELECT CASE (inform)
      CASE ('A')
         READ (iunit,*)
         READ (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat(:,:,l)
      CASE ('U')
         READ (iunit,IOSTAT=errstat,ERR=1000) intdat(:,:,l)
      CASE ('N')
         CALL cf90_get_var(iunit,idvars(l),intdat(:,:,l),timerec)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      END SELECT
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read integer array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_int_3d

!========================================================================

SUBROUTINE read_vars_int_4d(intdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_int_4d* Read a 4-D integer array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(OUT), DIMENSION(:,:,:,:) :: intdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER  Integer input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(intdat,DIM=4)) :: idvars


IF (SIZE(intdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_int_4d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(intdat,DIM=4)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      READ (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat
   CASE ('U')
      READ (iunit,IOSTAT=errstat,ERR=1000) intdat
   CASE ('N')
      CALL cf90_get_var(iunit,varid,intdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      SELECT CASE (inform)
      CASE ('A')
         READ (iunit,*)
         READ (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat(:,:,:,l)
      CASE ('U')
         READ (iunit,IOSTAT=errstat,ERR=1000) intdat(:,:,:,l)
      CASE ('N')
         CALL cf90_get_var(iunit,idvars(l),intdat(:,:,:,l),timerec)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      END SELECT
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read integer array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_int_4d

!========================================================================

SUBROUTINE read_vars_log_0d(logdat,filepars,varid,varatts)
!************************************************************************
!
! *read_vars_log_0d* Read a logical scalar variable in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
LOGICAL, INTENT(OUT) :: logdat
TYPE (FileParams), INTENT(IN) :: filepars
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*logdat*  LOGICAL     Logical input variable
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  variable id
!*varatts*  DERIVED  Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=1) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, novars, npcc, timerec

procname(pglev+1) = 'read_vars_log_0d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!
!1.1 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.2 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!

SELECT CASE (inform)
CASE ('A')
   READ (iunit,*)
   READ (iunit,LogicalFormat,IOSTAT=errstat,ERR=1000) logdat
CASE ('U')
   READ (iunit,IOSTAT=errstat,ERR=1000) logdat
CASE ('N')
   CALL cf90_get_var(iunit,varid,logdat,timerec)
   IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read logical variable'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_log_0d

!========================================================================

SUBROUTINE read_vars_log_1d(logdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_log_1d* Read a 1-D logical array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
LOGICAL, INTENT(OUT), DIMENSION(:) :: logdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*logdat*  LOGICAL     Logical input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: allrec, flag
CHARACTER (LEN=1) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(logdat)) :: idvars


IF (SIZE(logdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_log_1d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(logdat)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---type of input record
allrec = inform.NE.'N'.OR.varid.NE.0

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (allrec) THEN

   SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      READ (iunit,LogicalFormat,IOSTAT=errstat,ERR=1000) logdat
   CASE ('U')
      READ (iunit,IOSTAT=errstat,ERR=1000) logdat
   CASE ('N')
      CALL cf90_get_var(iunit,varid,logdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      CALL cf90_get_var(iunit,idvars(l),logdat(l),timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   ENDDO l_220
ENDIF

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read logical array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_log_1d

!========================================================================

SUBROUTINE read_vars_log_2d(logdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_log_2d* Read a 2-D logical array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
LOGICAL, INTENT(OUT), DIMENSION(:,:) :: logdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*logdat*   LOGICAL  Logical input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(logdat,DIM=2)) :: idvars


IF (SIZE(logdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_log_2d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(logdat,DIM=2)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!

iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      READ (iunit,LogicalFormat,IOSTAT=errstat,ERR=1000) logdat
   CASE ('U')
      READ (iunit,IOSTAT=errstat,ERR=1000) logdat
   CASE ('N')
      CALL cf90_get_var(iunit,varid,logdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      SELECT CASE (inform)
      CASE ('A')
         READ (iunit,*)
         READ (iunit,LogicalFormat,IOSTAT=errstat,ERR=1000) logdat(:,l)
      CASE ('U')
         READ (iunit,IOSTAT=errstat,ERR=1000) logdat(:,l)
      CASE ('N')
         CALL cf90_get_var(iunit,idvars(l),logdat(:,l),timerec)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      END SELECT
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read logical array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_log_2d

!========================================================================

SUBROUTINE read_vars_real_0d(realdat,filepars,varid,varatts)
!************************************************************************
!
! *read_vars_real_0d* Read a single precision real scalar variable in standard
!                     format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
REAL (KIND=kndreal), INTENT(OUT) :: realdat
TYPE (FileParams), INTENT(IN) :: filepars
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real input variable
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  variable id
!*varatts*  DERIVED  Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, inflag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, novars, npcc, timerec
REAL (KIND=kndrlong) :: realtmp


procname(pglev+1) = 'read_vars_real_0d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!
!1.1 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   inflag = varatts(varid)%data_type.EQ.real_type
ELSE
   inflag = filepars%floattype.EQ.'S'
ENDIF

!
!1.2 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!

SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      IF (inflag) THEN
         READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('U')
      IF (inflag) THEN
         READ (iunit,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         READ (iunit,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('N')
      IF (inflag) THEN
         CALL cf90_get_var(iunit,varid,realdat,timerec)
      ELSE
         CALL cf90_get_var(iunit,varid,realtmp,timerec)
      ENDIF
      IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

IF (.NOT.inflag) realdat = realtmp 

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read real variable'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_real_0d

!========================================================================

SUBROUTINE read_vars_real_1d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_real_1d* Read a 1-D single precision real array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndreal), INTENT(OUT), DIMENSION(:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: allrec, flag, inflag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(realdat)) :: idvars
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:) :: realtmp


IF (SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_real_1d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   inflag = varatts(varid)%data_type.EQ.real_type
ELSE
   inflag = filepars%floattype.EQ.'S'
ENDIF

!---type of input record
allrec = inform.NE.'N'.OR.varid.NE.0

!
!1.4 Allocate
!------------
!

IF (.NOT.inflag) THEN
   ALLOCATE (realtmp(nosize),STAT=errstat)
   CALL error_alloc('realtmp',1,(/nosize/),rlong_type)
ENDIF

!
!1.5 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (allrec) THEN

   SELECT CASE (inform)
      CASE ('A')
         READ (iunit,*)
         IF (inflag) THEN
            READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('U')
         IF (inflag) THEN
            READ (iunit,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('N')
         IF (inflag) THEN
            CALL cf90_get_var(iunit,varid,realdat,timerec)
         ELSE
            CALL cf90_get_var(iunit,varid,realtmp,timerec)
         ENDIF
         CALL cf90_get_var(iunit,varid,realdat,timerec)
         IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT
   IF (.NOT.inflag) realdat = realtmp 

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (inflag) THEN
      l_221: DO l=1,nosize
         CALL cf90_get_var(iunit,idvars(l),realdat(l),timerec)
      ENDDO l_221
   ELSE
      l_222: DO l=1,nosize
         CALL cf90_get_var(iunit,idvars(l),realtmp(l),timerec)
      END DO l_222
   ENDIF
   IF (errstat.NE.noerr_NF90) GOTO 1000
   
   IF (.NOT.inflag) realdat = realtmp 

ENDIF

!
!3. Deallocate
!-------------
!

IF (.NOT.inflag) DEALLOCATE (realtmp)


CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read real array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_real_1d
 
!========================================================================

SUBROUTINE read_vars_real_2d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_real_2d* Read a 2-D single precision real array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var, error_shape

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndreal), INTENT(OUT), DIMENSION(:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, inflag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(2) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=2)) :: idvars
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:,:) :: realtmp


IF (SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_real_2d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=2)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   inflag = varatts(varid)%data_type.EQ.real_type
ELSE
   inflag = filepars%floattype.EQ.'S'
ENDIF

!
!1.4 Allocate
!------------
!

IF (.NOT.inflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2)),STAT=errstat)
   CALL error_alloc('realtmp',2,nshape,rlong_type)
ENDIF

!
!1.5 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (inform)
      CASE ('A')
         READ (iunit,*)
         IF (inflag) THEN
            READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('U')
         IF (inflag) THEN
            READ (iunit,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('N')
         IF (inflag) THEN
            CALL cf90_get_var(iunit,varid,realdat,timerec)
         ELSE
            CALL cf90_get_var(iunit,varid,realtmp,timerec)
         ENDIF
         IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT
   
!
!2.2 Multiple records specified by last dimension   
!------------------------------------------------
!

ELSE

   IF (inflag) THEN

      l_221: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat(:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realdat(:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realdat(:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp(:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realtmp(:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222



   ENDIF
      
ENDIF

!
!3. Convert to single precision
!------------------------------
!

IF (.NOT.inflag) realdat = realtmp

!
!4. Deallocate
!-------------
!

IF (.NOT.inflag) DEALLOCATE (realtmp)

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read real array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_real_2d

!========================================================================

SUBROUTINE read_vars_real_3d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_real_3d* Read a 3-D single precision real array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndreal), INTENT(OUT), DIMENSION(:,:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, inflag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(3) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=3)) :: idvars
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realtmp


IF (SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_real_3d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=3)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   inflag = varatts(varid)%data_type.EQ.real_type
ELSE
   inflag = filepars%floattype.EQ.'S'
ENDIF

!
!1.4 Allocate
!------------
!

IF (.NOT.inflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2),nshape(3)),STAT=errstat)
   CALL error_alloc('realtmp',3,nshape,rlong_type)
ENDIF

!
!1.5 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (inform)
      CASE ('A')
         READ (iunit,*)
         IF (inflag) THEN
            READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('U')
         IF (inflag) THEN
            READ (iunit,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('N')
         IF (inflag) THEN
            CALL cf90_get_var(iunit,varid,realdat,timerec)
         ELSE
            CALL cf90_get_var(iunit,varid,realtmp,timerec)
         ENDIF
         IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (inflag) THEN
      
      l_221: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat(:,:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realdat(:,:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realdat(:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp(:,:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realtmp(:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222

   ENDIF

ENDIF

!
!3. Convert to single precision
!------------------------------
!

IF (.NOT.inflag) realdat = realtmp

!
!4. Deallocate
!-------------
!

IF (.NOT.inflag) DEALLOCATE (realtmp)

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read real array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_real_3d

!========================================================================

SUBROUTINE read_vars_real_4d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_real_4d* Read a 4-D single precision real array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndreal), INTENT(OUT), DIMENSION(:,:,:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, inflag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(4) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=4)) :: idvars
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: realtmp


IF (SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_real_4d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=4)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   inflag = varatts(varid)%data_type.EQ.real_type
ELSE
   inflag = filepars%floattype.EQ.'S'
ENDIF

!
!1.4 Allocate
!------------
!

IF (.NOT.inflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2),nshape(3),nshape(4)),STAT=errstat)
   CALL error_alloc('realtmp',4,nshape,rlong_type)
ENDIF

!
!1.5 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (inform)
      CASE ('A')
         READ (iunit,*)
         IF (inflag) THEN
            READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('U')
         IF (inflag) THEN
            READ (iunit,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('N')
         IF (inflag) THEN
            CALL cf90_get_var(iunit,varid,realdat,timerec)
         ELSE
            CALL cf90_get_var(iunit,varid,realtmp,timerec)
         ENDIF
         IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (inflag) THEN

      l_221: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat(:,:,:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realdat(:,:,:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realdat(:,:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) &
                                                              & realtmp(:,:,:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,:,:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realtmp(:,:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222

   ENDIF

ENDIF

!
!3. Convert to single precision
!------------------------------
!

IF (.NOT.inflag) realdat = realtmp

!
!4. Deallocate
!-------------
!

IF (.NOT.inflag) DEALLOCATE (realtmp)


CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read real array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_real_4d

!========================================================================

SUBROUTINE read_vars_double_0d(realdat,filepars,varid,varatts)
!************************************************************************
!
! *read_vars_double_0d* Read a double precision real scalar variable in
!                       standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
REAL (KIND=kndrlong), INTENT(OUT) :: realdat
TYPE (FileParams), INTENT(IN) :: filepars
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  DOUBLE   Real input variable
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  variable id
!*varatts*  DERIVED  Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, inflag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, novars, npcc, timerec
REAL (KIND=kndreal) :: realtmp


procname(pglev+1) = 'read_vars_double_0d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!
!1.1 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data data are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   inflag = varatts(varid)%data_type.EQ.rlong_type
ELSE
   inflag = filepars%floattype.EQ.'D'
ENDIF

!
!1.2 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!

SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      IF (inflag) THEN
         READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('U')
      IF (inflag) THEN
         READ (iunit,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         READ (iunit,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('N')
      IF (inflag) THEN
         CALL cf90_get_var(iunit,varid,realdat,timerec)
      ELSE
         CALL cf90_get_var(iunit,varid,realtmp,timerec)
      ENDIF
      IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

IF (.NOT.inflag) realdat = realtmp 

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read real variable'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_double_0d

!========================================================================

SUBROUTINE read_vars_double_1d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_double_1d* Read a 1-D double precision real array in standard
!                       format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndrlong), INTENT(OUT), DIMENSION(:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: allrec, flag, inflag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(realdat)) :: idvars
REAL (KIND=kndreal), SAVE, ALLOCATABLE, DIMENSION(:) :: realtmp


IF (SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_double_1d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   inflag = varatts(varid)%data_type.EQ.rlong_type
ELSE
   inflag = filepars%floattype.EQ.'D'
ENDIF

!---type of input record
allrec = inform.NE.'N'.OR.varid.NE.0

!
!1.4 Allocate
!------------
!

IF (.NOT.inflag) THEN
   ALLOCATE (realtmp(nosize),STAT=errstat)
   CALL error_alloc('realtmp',1,(/nosize/),real_type)
ENDIF

!
!1.5 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (allrec) THEN

   SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      IF (inflag) THEN
         READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('U')
      IF (inflag) THEN
         READ (iunit,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         READ (iunit,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('N')
      IF (inflag) THEN
         CALL cf90_get_var(iunit,varid,realdat,timerec)
      ELSE
         CALL cf90_get_var(iunit,varid,realtmp,timerec)
      ENDIF
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT
   IF (.NOT.inflag) realdat = realtmp 

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (inflag) THEN
      l_221: DO l=1,nosize
         CALL cf90_get_var(iunit,idvars(l),realdat(l),timerec)
      ENDDO l_221
   ELSE
      l_222: DO l=1,nosize
         CALL cf90_get_var(iunit,idvars(l),realtmp(l),timerec)
      END DO l_222
   ENDIF
   IF (errstat.NE.noerr_NF90) GOTO 1000
   
   IF (.NOT.inflag) realdat = realtmp 

ENDIF

!
!3. Deallocate
!-------------
!

IF (.NOT.inflag) DEALLOCATE (realtmp)


CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read real array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_double_1d

!========================================================================

SUBROUTINE read_vars_double_2d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_double_2d* Read a 2-D double precision real array in standard
!                       format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndrlong), INTENT(OUT), DIMENSION(:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, inflag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(2) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=2)) :: idvars
REAL (KIND=kndreal), SAVE, ALLOCATABLE, DIMENSION(:,:) :: realtmp


IF (SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_double_2d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=2)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   inflag = varatts(varid)%data_type.EQ.rlong_type
ELSE
   inflag = filepars%floattype.EQ.'D'
ENDIF

!
!1.4 Allocate
!------------
!

IF (.NOT.inflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2)),STAT=errstat)
   CALL error_alloc('realtmp',2,nshape,rlong_type)
ENDIF

!
!1.5 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (inform)
      CASE ('A')
         READ (iunit,*)
         IF (inflag) THEN
            READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('U')
         IF (inflag) THEN
            READ (iunit,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('N')
         IF (inflag) THEN
            CALL cf90_get_var(iunit,varid,realdat,timerec)
         ELSE
            CALL cf90_get_var(iunit,varid,realtmp,timerec)
         ENDIF
         IF (errstat.NE.noerr_NF90) GOTO 1000
      END SELECT

   IF (.NOT.inflag) realdat = realtmp
   
!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (inflag) THEN

      l_221: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat(:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realdat(:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realdat(:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp(:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realtmp(:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222

   ENDIF
      
ENDIF
 
!
!3. Convert to single precision
!------------------------------
!

IF (.NOT.inflag) realdat = realtmp

!
!4. Deallocate
!-------------
!

IF (.NOT.inflag) DEALLOCATE (realtmp)

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_double_2d

!========================================================================

SUBROUTINE read_vars_double_3d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_double_3d* Read a 3-D double precision real array in standard!
!                       format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndrlong), INTENT(OUT), DIMENSION(:,:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, inflag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(3) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=3)) :: idvars
REAL (KIND=kndreal), SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realtmp

IF (SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_double_3d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=3)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   inflag = varatts(varid)%data_type.EQ.rlong_type
ELSE
   inflag = filepars%floattype.EQ.'D'
ENDIF

!
!1.4 Allocate
!------------
!

IF (.NOT.inflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2),nshape(3)),STAT=errstat)
   CALL error_alloc('realtmp',3,nshape,real_type)
ENDIF

!
!1.5 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (inform)
      CASE ('A')
         READ (iunit,*)
         IF (inflag) THEN
            READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('U')
         IF (inflag) THEN
            READ (iunit,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('N')
         IF (inflag) THEN
            CALL cf90_get_var(iunit,varid,realdat,timerec)
         ELSE
            CALL cf90_get_var(iunit,varid,realtmp,timerec)
         ENDIF
         IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (inflag) THEN
      
      l_221: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat(:,:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realdat(:,:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realdat(:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp(:,:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realtmp(:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222

   ENDIF

ENDIF
 
!
!3. Convert to single precision
!------------------------------
!

IF (.NOT.inflag) realdat = realtmp

!
!4. Deallocate
!-------------
!

IF (.NOT.inflag) DEALLOCATE (realtmp)

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read real array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_double_3d

!========================================================================

SUBROUTINE read_vars_double_4d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *read_vars_double_4d* Read a 4-D double precision real array in standard
!                       format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_get_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_get_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndrlong), INTENT(OUT), DIMENSION(:,:,:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real input array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of input array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, inflag
CHARACTER (LEN=lenform) :: inform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: ierr, iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(4) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=4)) :: idvars
REAL (KIND=kndreal), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: realtmp

IF (SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'read_vars_double_4d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=4)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for reading
!--------------------------
!
!---file attributes
iunit = filepars%iunit
inform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   inflag = varatts(varid)%data_type.EQ.rlong_type
ELSE
   inflag = filepars%floattype.EQ.'D'
ENDIF

!
!1.4 Allocate
!------------
!

IF (.NOT.inflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2),nshape(3),nshape(4)),STAT=errstat)
   CALL error_alloc('realtmp',4,nshape,real_type)
ENDIF

!
!1.5 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,inform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Read data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (inform)
   CASE ('A')
      READ (iunit,*)
      IF (inflag) THEN
            READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('U')
         IF (inflag) THEN
            READ (iunit,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            READ (iunit,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('N')
         IF (inflag) THEN
            CALL cf90_get_var(iunit,varid,realdat,timerec)
         ELSE
            CALL cf90_get_var(iunit,varid,realtmp,timerec)
         ENDIF
         IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (inflag) THEN
      
      l_221: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) &
                                                     & realdat(:,:,:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realdat(:,:,:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realdat(:,:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (inform)
            CASE ('A')
               READ (iunit,*)
               READ (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp(:,:,:,l)
            CASE ('U')
               READ (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,:,:,l)
            CASE ('N')
               CALL cf90_get_var(iunit,idvars(l),realtmp(:,:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222

   ENDIF

ENDIF
 
!
!3. Convert to single precision
!------------------------------
!

IF (.NOT.inflag) realdat = realtmp

!
!4. Deallocate
!-------------
!

IF (.NOT.inflag) DEALLOCATE (realtmp)

CALL log_timer_out(npcc,itm_input)


RETURN

!---process read error
1000 CONTINUE
IF (inform.NE.'N') THEN
   ierr = MERGE(ierrno_read,ierrno_fend,errstat.GT.0)
ELSE
   ierr = ierrno_read
ENDIF
IF (errchk) WRITE (ioerr,'(A)') 'Unable to read real array'
CALL error_file(ierr,filepars=filepars)

END SUBROUTINE read_vars_double_4d

!========================================================================

SUBROUTINE write_atts_mod(filepars,varatts,numvars,tdim)
!************************************************************************
!
! *write_atts_mod* Write the atrributes (global and variable) of a model file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_def_dim, cf90_def_var, cf90_enddef, cf90_put_att,
!                cf90_put_att_chars, conv_to_chars, conv_to_chars_modvars,
!                error_file, inquire_unit, reset_mod_vars, write_metadata_line
!
!************************************************************************
!
USE datatypes
USE cf90_routines, ONLY: cf90_def_dim, cf90_def_var, cf90_enddef, &
                       & cf90_put_att, cf90_put_att_chars
USE cif_routines, ONLY: conv_to_chars, conv_to_chars_modvars
USE error_routines, ONLY: error_file
USE reset_model, ONLY: reset_mod_vars

!
!*Arguments
!
INTEGER, INTENT(IN) :: numvars
INTEGER, INTENT(IN), OPTIONAL :: tdim
TYPE (FileParams), INTENT(INOUT) :: filepars
TYPE (VariableAtts), INTENT(INOUT), DIMENSION(numvars) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*filepars* DERIVED  File attributes
!*varatts*  DERIVED  Variable attributes
!*numvars*  INTEGER  Number of (including coordinate) variables in file
!*tdim*     INTEGER  Dimension of the time variable (unlimited if not present)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, packing
CHARACTER (LEN=1) :: floattype
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=4) :: cdim
CHARACTER (LEN=lendesc) :: comment, long_name, standard_name
CHARACTER (LEN=lenunit) :: units
CHARACTER (LEN=lencifline) :: cline
INTEGER :: icoord, iglb, ipack, iunit, ivar, l, n, nocoords, nodims, &
         & noglbatts, novars, npcc, nrank, nstepout, nvals, varid
CHARACTER (LEN=1), DIMENSION(4) :: cdir = (/'X','Y','Z','V'/)
CHARACTER (LEN=lendesc), DIMENSION(14) :: cvals
INTEGER, DIMENSION(2) :: dimstime, dimtimeid
INTEGER, DIMENSION(5) :: dimids


procname(pglev+1) = 'write_atts_mod'
CALL log_timer_in(npcc,logname=TRIM(procname(pglev+1))//': '//filepars%filename)

!
!1. File parameters
!------------------
!

iunit = filepars%iunit
outform = filepars%form
nocoords = filepars%nocoords
novars = filepars%novars
nstepout = filepars%nstepout
floattype = filepars%floattype
packing = filepars%packing
ipack = MERGE(1,0,packing)

!
!2. Check whether file unit is opened
!------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!3. Global attributes
!--------------------
!

glbatts(2)%value = TRIM(filepars%title)

!
!4. Write global attributes
!--------------------------
!

SELECT CASE (outform)
   
!  ---ASCII, binary
   CASE ('A','U')
      noglbatts = COUNT(LEN_TRIM(glbatts%value).GT.0)
      CALL conv_to_chars(cvals(1),noglbatts)
      CALL write_metadata_line(iunit,cvals(1:1),'nGlbAtts',outform)
      iglb_411: DO iglb=1,numglbatts
         IF (TRIM(glbatts(iglb)%value).NE.''.AND.&
           & TRIM(glbatts(iglb)%value).NE.'netcdf') THEN
            cvals(1) = TRIM(glbatts(iglb)%value)
            CALL write_metadata_line(iunit,cvals(1:1),TRIM(glbatts(iglb)%name),&
                                   & outform)
         ENDIF
      ENDDO iglb_411
      CALL conv_to_chars(cvals(1),nocoords)
      CALL write_metadata_line(iunit,cvals(1:1),'nCoordinates',outform)
      CALL conv_to_chars(cvals(1),novars)
      CALL write_metadata_line(iunit,cvals(1:1),'nVariables',outform)
      cvals(1) = floattype
      CALL write_metadata_line(iunit,cvals(1:1),'Floattype',outform)
      CALL conv_to_chars(cvals(1),packing)
      CALL write_metadata_line(iunit,cvals(1:1),'packing',outform)

!  ---netCDF
   CASE ('N')
      iglb_412: DO iglb=1,numglbatts
         IF (TRIM(glbatts(iglb)%value).NE.'') THEN
            l = LEN_TRIM(glbatts(iglb)%value)
            CALL cf90_put_att_chars(iunit,global_NF90,TRIM(glbatts(iglb)%name),&
                                  & l,glbatts(iglb)%value)
         ENDIF
      ENDDO iglb_412
      CALL cf90_put_att(iunit,global_NF90,'packing',ipack)
      
END SELECT

!
!5. Define variable attributes
!-----------------------------
!

ivar_510: DO ivar=1,numvars
   CALL reset_mod_vars(varatts(ivar))
   varatts(ivar)%comment = ''
ENDDO ivar_510

!
!6. Write variable attributes
!----------------------------
!

SELECT CASE (outform)

   CASE ('A','U')

      IF (nocoords.EQ.1) filepars%timeid = 1
      
      icoord_610: DO icoord=1,nocoords
         CALL conv_to_chars(cvals(1),icoord)
         nvals = 11 + varatts(icoord)%nrank
         CALL conv_to_chars_modvars(cvals(2:nvals),varatts(icoord))
         CALL write_metadata_line(iunit,cvals(1:nvals),'coordatts',outform)
      ENDDO icoord_610

      ivar_610: DO ivar=nocoords+1,numvars
         CALL conv_to_chars(cvals(1),ivar)
         nvals = 11 + varatts(ivar)%nrank
         CALL conv_to_chars_modvars(cvals(2:nvals),varatts(ivar))
         CALL write_metadata_line(iunit,cvals(1:nvals),'varatts',outform)
      ENDDO ivar_610

   CASE ('N')

!     ---define time dimension
      IF (nocoords.EQ.1) THEN
         IF (PRESENT(tdim)) THEN
            dimstime = (/lentime,tdim/)
         ELSE
            dimstime = (/lentime,unlimited_NF90/)
         ENDIF
         CALL cf90_def_dim(iunit,'tlendim',dimstime(1),dimtimeid(1))
         CALL cf90_def_dim(iunit,'time',dimstime(2),dimtimeid(2))
      ENDIF

!     ---define time coordinate
      IF (nocoords.EQ.1) THEN
         CALL cf90_def_var(iunit,varatts(1)%f90_name,char_NF90,dimtimeid,&
                         & filepars%timeid)
         standard_name = varatts(1)%standard_name
         CALL cf90_put_att_chars(iunit,1,'standard_name',lendesc,standard_name)
         long_name = varatts(1)%long_name
         CALL cf90_put_att_chars(iunit,1,'long_name',lendesc,long_name)
         units = varatts(1)%units
         CALL cf90_put_att_chars(iunit,1,'units',lenunit,units)
         comment = varatts(1)%comment
         IF (TRIM(comment).NE.'') THEN
            CALL cf90_put_att_chars(iunit,1,'comment',lendesc,comment)
         ENDIF
      ENDIF

!     ---define variables
      ivar_630: DO ivar=nocoords+1,numvars
         nrank = varatts(ivar)%nrank
         IF (nocoords.EQ.1) THEN
            nodims = nrank + 1
            dimids(nodims) = dimtimeid(2)
         ELSE
            nodims = nrank
         ENDIF
         n_631: DO n=1,nrank
            WRITE (cdim,'(A,I3.3)') cdir(n), ivar
            CALL cf90_def_dim(iunit,cdim,varatts(ivar)%global_dims(n),&
                            & dimids(n))
         ENDDO n_631
         IF (varatts(ivar)%data_type.EQ.int_type.OR.&
           & varatts(ivar)%data_type.EQ.log_type) THEN
            CALL cf90_def_var(iunit,varatts(ivar)%f90_name,int_NF90,&
                            & dimids(1:nodims),varid)
         ELSEIF (varatts(ivar)%data_type.EQ.real_type) THEN
            CALL cf90_def_var(iunit,varatts(ivar)%f90_name,real_NF90,&
                            & dimids(1:nodims),varid)
         ELSEIF (varatts(ivar)%data_type.EQ.rlong_type) THEN
            CALL cf90_def_var(iunit,varatts(ivar)%f90_name,double_NF90,&
                            & dimids(1:nodims),varid)
         ENDIF
         standard_name = varatts(ivar)%standard_name
         CALL cf90_put_att_chars(iunit,varid,'standard_name',lendesc,&
                               & standard_name)
         long_name = varatts(ivar)%long_name
         CALL cf90_put_att_chars(iunit,varid,'long_name',lendesc,long_name)
         units = varatts(ivar)%units
         CALL cf90_put_att_chars(iunit,varid,'units',lenunit,units)
         IF (varatts(ivar)%fill) THEN
            IF (floattype.EQ.'S') THEN
               CALL cf90_put_att(iunit,varid,'_FillValue',real_fill)
            ELSE
               CALL cf90_put_att(iunit,varid,'_FillValue',double_fill)
            ENDIF
         ENDIF
      ENDDO ivar_630

END SELECT

!
!7. End of header
!----------------
!

cline = '#'
IF (outform.EQ.'A') THEN
   WRITE (iunit,'(A)') TRIM(cline)
ELSEIF (outform.EQ.'U') THEN
   WRITE (iunit) cline
ENDIF

!
!8. Initialise time record number
!--------------------------------
!

filepars%timerec = 0

!
!9. Leave define mode
!--------------------
!

IF (outform.EQ.'N') CALL cf90_enddef(iunit)

!
!10. Write information file
!--------------------------
!

IF (filepars%info) CALL write_info_mod(filepars,varatts,numvars)

IF (nerrs.GT.0) CALL error_file(ierrno_write,filepars=filepars)

CALL log_timer_out(npcc,itm_output)


RETURN

END SUBROUTINE write_atts_mod

!========================================================================

SUBROUTINE write_info_mod(filepars,varatts,numvars)
!************************************************************************
!
! *write_info_mod* Write a model information file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.7.1
!
! Description -
!
! Module calls - close_file, open_file, write_metadata_line
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cif_routines, ONLY: conv_to_chars

!
!*Arguments
!
INTEGER, INTENT(IN) :: numvars
TYPE (FileParams), INTENT(IN) :: filepars
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(numvars) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*filepars* DERIVED  File attributes
!*varatts*  DERIVED  Variable attributes
!*numvars*  INTEGER  Number of (including coordinate) variables in file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=leniofile) :: infofile
CHARACTER (LEN=lendesc), DIMENSION(1) :: cvals
INTEGER :: data_type, iglb, iunit, ivar, ivarid, nocoords, noglbatts, novars, &
         & nrank


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_info_mod'
CALL log_timer_in()

!
!1. Open file
!------------
!

infofile = TRIM(filepars%filename)//'.inf'
CALL open_file(iunit,infofile,'OUT','A')

!
!2. Write info
!-------------
!
!---global attributes
noglbatts = COUNT(glbatts%value.NE.'')
CALL conv_to_chars(cvals(1),noglbatts)
CALL write_metadata_line(iunit,cvals(1:1),'nGlobalAtts','A')
iglb_210: DO iglb=1,numglbatts
   IF (TRIM(glbatts(iglb)%value).NE.'') THEN
      cvals(1) = TRIM(glbatts(iglb)%value)
      CALL write_metadata_line(iunit,cvals(1:1),&
                           & TRIM(glbatts(iglb)%name),'A')
   ENDIF
ENDDO iglb_210

nocoords = filepars%nocoords
CALL conv_to_chars(cvals(1),nocoords)
CALL write_metadata_line(iunit,cvals(1:1),'nCoordinates','A')
novars = filepars%novars
CALL conv_to_chars(cvals(1),novars)
CALL write_metadata_line(iunit,cvals(1:1),'nVariables','A')

!---coordinate and variable attributes
ivar_220: DO ivar=1,numvars
   ivarid = varatts(ivar)%ivarid
   data_type = varatts(ivar)%data_type
   nrank = varatts(ivar)%nrank
   WRITE (iunit,*) ivarid, data_type, nrank, &
                 & varatts(ivar)%global_dims(1:nrank)
   WRITE (iunit,'(A)') TRIM(varatts(ivar)%f90_name)
   WRITE (iunit,'(A)') TRIM(varatts(ivar)%long_name)
   IF (LEN_TRIM(varatts(ivar)%units).EQ.0) THEN
      WRITE (iunit,'(A)') '_'
   ELSE
      WRITE (iunit,'(A)') TRIM(varatts(ivar)%units)
   ENDIF
ENDDO ivar_220

!
!3. Close
!--------
!

CALL close_file(iunit,'A',filename=infofile)

CALL log_timer_out()


RETURN

END SUBROUTINE write_info_mod

!========================================================================

SUBROUTINE write_metadata_line(iunit,cvals,cname,form)
!************************************************************************
!
! *write_metadata_line* write header info in output files
!
! Author - SORESMA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.7.1
!
! Description - 
!
! Module calls -
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: cname
CHARACTER (LEN=*), INTENT(INOUT), DIMENSION(:) :: cvals
CHARACTER (len=lenform), INTENT(IN), OPTIONAL :: form
INTEGER, INTENT(IN) :: iunit

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER CIF file key id
!*cvals*     CHAR    Data values in string format
!*cname*     CHAR    Name of the variable(s)
!*form*      CHAR    Format of the output file ('A' or 'U')
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=lencifline) :: cline
INTEGER :: n, nosize


nosize = SIZE(cvals)

IF (nosize.GT.1) THEN
   cline = TRIM(cvals(1))
   n_110: DO n=2,nosize
      cline = TRIM(cline)//','//TRIM(cvals(n))
   ENDDO n_110
ELSE
   cline = cvals(1)
ENDIF

IF (form.EQ.'A') THEN
   WRITE (iunit,'(A)') TRIM(cname)//' : '//TRIM(cline)
ELSE
   cline = TRIM(cname)//' : '//TRIM(cline)
   WRITE (iunit) cline
ENDIF


RETURN

END SUBROUTINE write_metadata_line

!========================================================================

SUBROUTINE write_submod_real_2d(realsub,filepars,varid,nowet,varatts)
!************************************************************************
!
! *write_submod_real_2d* Write a global sub-model real 2-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.11
!
! Description -
!
! Module calls - error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE paralpars  
USE syspars
USE error_routines, ONLY: error_alloc

!
!*  Arguments
!
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: nowet
TYPE(FileParams), INTENT(INOUT) :: filepars
REAL, INTENT(IN), DIMENSION(:,:) :: realsub
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realsub*   REAL     Output sub-model 2-D array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*nowet*     INTEGER  Number of no-wet horizontal data points in case the data
!                     are witten in "packed" (vector) format
!*varatts*   DERIVED  Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fill, packing
INTEGER :: ivarid
REAL :: fill_value, fill_value_eps
REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: realdat


IF (SIZE(realsub).EQ.0) RETURN
procname(pglev+1) = 'write_submod_real_2d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

!
!1. Initialise parameters and arrays
!-----------------------------------
!
!---fill parameters
fill = filepars%fill
IF (fill) THEN
   fill_value = double_fill
   fill_value_eps = 0.00001*ABS(fill_value)
ELSE
   fill_value = 0.0
ENDIF

!---packing
packing = filepars%packing.AND.mastermod

!---data array for packing
IF (packing) THEN
   ALLOCATE (realdat(nowet),STAT=errstat)
   CALL error_alloc('realdat',1,(/nowet/),kndrtype)
ENDIF

!
!2. Pack data using land mask
!----------------------------
!

IF (packing) THEN
   realdat = PACK(realsub,MASK=ABS(realsub-fill_value).GT.fill_value_eps)
ENDIF

!
!3. Write data
!-------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realdat,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(realdat,filepars,varid)
   ENDIF
ELSE
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realsub,filepars,varid,varatts=varatts)
   ELSE
      CALL write_vars(realsub,filepars,varid)
   ENDIF
ENDIF

!
!4. Deallocate array
!-------------------
!

IF (packing) DEALLOCATE (realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE write_submod_real_2d

!========================================================================

SUBROUTINE write_submod_real_3d(realsub,filepars,varid,nowet,varatts,vecids)
!************************************************************************
!
! *write_submod_real_3d* Write a global sub-model real 3-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.11
!
! Description -
!
! Module calls - error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE paralpars  
USE syspars
USE error_routines, ONLY: error_alloc

!
!*  Arguments
!
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: nowet
TYPE(FileParams), INTENT(INOUT) :: filepars
REAL, INTENT(IN), DIMENSION(:,:,:) :: realsub
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realsub*   REAL     Output sub-model 3-D array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*nowet*     INTEGER  Number of no-wet horizontal data points in case the data
!                     are witten in "packed" (vector) format
!*varatts*   DERIVED  Variable attributes
!*vecids*    INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fill, packing
INTEGER :: ivar, ivarid, k, n3dim
REAL :: fill_value, fill_value_eps
INTEGER, DIMENSION(SIZE(realsub,DIM=3)) :: idvars 
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: realdat


IF (SIZE(realsub).EQ.0) RETURN
procname(pglev+1) = 'write_submod_real_3d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

n3dim = SIZE(realsub,DIM=3)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSE
   idvars = (/(ivar,ivar=1,n3dim)/)
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---fill parameters
fill = filepars%fill
IF (fill) THEN
   fill_value = double_fill
   fill_value_eps = 0.00001*ABS(fill_value)
ELSE
   fill_value = 0.0
ENDIF

!---packing
packing = filepars%packing.AND.mastermod

!---data array for packing
IF (packing) THEN
   ALLOCATE (realdat(nowet,n3dim),STAT=errstat)
   CALL error_alloc('realdat',2,(/nowet,n3dim/),kndrtype)
ENDIF

!
!3. Pack data using land mask
!----------------------------
!

IF (packing) THEN
   k_310: DO k=1,n3dim
      realdat(:,k) = PACK(realsub(:,:,k),&
                   & MASK=ABS(realsub(:,:,k)-fill_value).GT.fill_value_eps)
   ENDDO k_310
ENDIF

!
!4. Write data
!-------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realdat,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realdat,filepars,varid,vecids=idvars)
   ENDIF
ELSE
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realsub,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realsub,filepars,varid,vecids=idvars)
   ENDIF
ENDIF

!
!5. Deallocate array
!-------------------
!

IF (packing) DEALLOCATE (realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE write_submod_real_3d

!========================================================================

SUBROUTINE write_submod_real_4d(realsub,filepars,varid,nowet,varatts,vecids)
!************************************************************************
!
! *write_submod_real_4d* Write a global sub-model real 4-D array
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_paral.f90  V2.11
!
! Description -
!
! Module calls - error_alloc, write_vars
!
!************************************************************************
!
USE datatypes
USE paralpars
USE syspars
USE error_routines, ONLY: error_alloc

!
!*  Arguments
!
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: nowet
TYPE(FileParams), INTENT(INOUT) :: filepars
REAL, INTENT(IN), DIMENSION(:,:,:,:) :: realsub
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name       Type      Purpose
!------------------------------------------------------------------------------
!*realsub*   REAL     Output sub-model 4-D array
!*filepars*  DERIVED  File attributes
!*varid*     INTEGER  If >0, variable id of input array.
!                     If =0, last dimension is a variable dimension
!*nowet*     INTEGER  Number of no-wet horizontal data points in case the data
!                     are witten in "packed" (vector) format
!*varatts*   DERIVED  Variable attributes
!*vecids*    INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: fill, packing
INTEGER :: ivar, ivarid, k, l, n3dim, n4dim
REAL :: fill_value, fill_value_eps
INTEGER, DIMENSION(SIZE(realsub,DIM=3)) :: idvars 
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realdat


IF (SIZE(realsub).EQ.0) RETURN
procname(pglev+1) = 'write_submod_real_3d'
IF (varid.GT.0.AND.PRESENT(varatts)) THEN
   ivarid = varatts(varid)%ivarid
ELSE
   ivarid = 0
ENDIF
CALL log_timer_in(ivarid=ivarid)

n3dim = SIZE(realsub,DIM=3)
n4dim = SIZE(realsub,DIM=4)

!
!1. Optional arguments
!---------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSE
   idvars = (/(ivar,ivar=1,n3dim)/)
ENDIF

!
!2. Initialise parameters and arrays
!-----------------------------------
!
!---fill parameters
fill = filepars%fill
IF (fill) THEN
   fill_value = double_fill
   fill_value_eps = 0.00001*ABS(fill_value)
ELSE
   fill_value = 0.0
ENDIF

!---packing
packing = filepars%packing.AND.mastermod

!---data array for packing
IF (packing) THEN
   ALLOCATE (realdat(nowet,n3dim,n4dim),STAT=errstat)
   CALL error_alloc('realdat',3,(/nowet,n3dim,n4dim/),kndrtype)
ENDIF

!
!3. Pack data using land mask
!----------------------------
!

IF (packing) THEN
   l_410: DO l=1,n4dim
   k_310: DO k=1,n3dim
      realdat(:,k,l) = PACK(realsub(:,:,k,l),&
                     & MASK=ABS(realsub(:,:,k,l)-fill_value).GT.fill_value_eps)
   ENDDO k_310
   ENDDO l_410
ENDIF

!
!4. Write data
!-------------
!

IF (packing) THEN
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realdat,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realdat,filepars,varid,vecids=idvars)
   ENDIF
ELSE
   IF (PRESENT(varatts)) THEN
      CALL write_vars(realsub,filepars,varid,varatts=varatts,vecids=idvars)
   ELSE
      CALL write_vars(realsub,filepars,varid,vecids=idvars)
   ENDIF
ENDIF

!
!5. Deallocate array
!-------------------
!

IF (packing) DEALLOCATE (realdat)

CALL log_timer_out()


RETURN

END SUBROUTINE write_submod_real_4d

!========================================================================

SUBROUTINE write_time_char(ciodatetime,filepars,timerec,back)
!************************************************************************
!
! *write_time_char* Write a date/time in string format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var_chars, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var_chars
USE error_routines, ONLY: error_file

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: back
INTEGER, INTENT(IN), OPTIONAL :: timerec
CHARACTER (LEN=lentime), INTENT(IN)  :: ciodatetime
TYPE(FileParams), INTENT(INOUT) :: filepars

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR    Date/time
!*filepars*    DERIVED File attributes
!*timerec*     INTEGER Time record number
!*back*        LOGICAL If present and .TRUE., writing is performed backwards
!                      (i.e. starting with the last time record)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, npcc, nsign, timerecx


IF (.NOT.mastermod) RETURN

procname(pglev+1) = 'write_time_char'
CALL log_timer_in(npcc,&
                & logname=TRIM(procname(pglev+1))//': '//filepars%filename)

!
!1. Check and initialise
!-----------------------
!
!---parameters for reading
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
IF (PRESENT(timerec)) THEN
   timerecx = timerec
ELSE
   IF (PRESENT(back)) THEN
      nsign = MERGE (-1,1,back)
   ELSE
      nsign = 1
   ENDIF
   timerecx = filepars%timerec + nsign
   filepars%timerec = timerecx
ENDIF

!---check whether file unit is opened
flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!

SELECT CASE (outform)
CASE ('A')
   WRITE (iunit,'(A)',IOSTAT=errstat,ERR=1000) ciodatetime
CASE ('U')
   WRITE (iunit,IOSTAT=errstat,ERR=1000) ciodatetime
CASE ('N')
   CALL cf90_put_var_chars(iunit,filepars%timeid,ciodatetime,timerecx)
   IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write time'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_time_char

!========================================================================

SUBROUTINE write_time_real(realtime,filepars,timerec,back)
!************************************************************************
!
! *write_time_real* Write a real date/time
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_file

!
!*Arguments
!
LOGICAL, INTENT(IN), OPTIONAL :: back
INTEGER, INTENT(IN), OPTIONAL :: timerec
REAL( KIND=kndrlong), INTENT(IN)  :: realtime
TYPE(FileParams), INTENT(INOUT) :: filepars

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*realtime*    REAL    Date/time
!*filepars*    DERIVED File attributes
!*timerec*     INTEGER Time record number
!*back*        LOGICAL If present and .TRUE., writing is performed backwards
!                      (i.e. starting with the last time record)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, npcc, nsign, timerecx


IF (.NOT.mastermod) RETURN

procname(pglev+1) = 'write_time_real'
CALL log_timer_in(npcc,&
                & logname=TRIM(procname(pglev+1))//': '//filepars%filename)

!
!1. Check and initialise
!-----------------------
!
!---parameters for reading
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
IF (PRESENT(timerec)) THEN
   timerecx = timerec
ELSE
   IF (PRESENT(back)) THEN
      nsign = MERGE (-1,1,back)
   ELSE
      nsign = 1
   ENDIF
   timerecx = filepars%timerec + nsign
   filepars%timerec = timerecx
ENDIF

!---check whether file unit is opened
flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!

SELECT CASE (outform)
CASE ('A')
   WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtime
CASE ('U')
   WRITE (iunit,IOSTAT=errstat,ERR=1000) realtime
CASE ('N')
   CALL cf90_put_var(iunit,filepars%timeid,realtime,timerecx)
   IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write time'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_time_real

!========================================================================

SUBROUTINE write_vars_char_0d(chardat,filepars,varid,varatts)
!************************************************************************
!
! *write_vars_char_0d* Write a character scalar variable in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var_chars, error_abort, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var_chars
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
CHARACTER (LEN=*), INTENT(IN) :: chardat
TYPE (FileParams), INTENT(IN) :: filepars
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*chardat*  CHAR     Character output variable
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  variable id
!*varatts*  DERIVED  Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, varflag
CHARACTER (LEN=1) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, novars, npcc, timerec


IF (.NOT.mastermod) RETURN

procname(pglev+1) = 'write_vars_char_0d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!
!1.1 Parameters for writing
!--------------------------
!

varflag = PRESENT(varatts)
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.2 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!-------------
!

SELECT CASE (outform)
CASE ('A')
   IF (varflag) THEN
      WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
   ELSE
      WRITE (iunit,*)
   ENDIF
   WRITE (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) chardat
CASE ('U')
   WRITE (iunit,IOSTAT=errstat,ERR=1000) chardat
CASE ('N')
   CALL cf90_put_var_chars(iunit,varid,chardat,timerec)
   IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write character variable'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_char_0d

!========================================================================

SUBROUTINE write_vars_char_1d(chardat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_char_1d* Write a 1-D character array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var_chars, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var_chars
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
CHARACTER (LEN=*), INTENT(IN), DIMENSION(:) :: chardat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*chardat*  CHAR     Character output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: allrec, flag, varflag
CHARACTER (LEN=1) :: outform
CHARACTER (LEN=12) :: clen, csize
CHARACTER (LEN=29) :: cfmt
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(chardat)) :: idvars


IF (.NOT.mastermod.OR.SIZE(chardat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_char_1d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(chardat)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---type of input record
allrec = outform.NE.'N'.OR.varid.NE.0

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!-------------
!
!2.1 As one record
!-----------------
!

IF (allrec) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag.AND.varid.GT.0) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      WRITE (csize,'(I12)') nosize; csize = ADJUSTL(csize)
      WRITE (clen,'(I12)') LEN(chardat(1)); clen = ADJUSTL(clen)
      cfmt = '('//TRIM(csize)//'A'//TRIM(clen)//')'
      WRITE (iunit,TRIM(cfmt),IOSTAT=errstat,ERR=1000) chardat
   CASE ('U')
      WRITE (iunit,IOSTAT=errstat,ERR=1000) chardat
   CASE ('N')
      CALL cf90_put_var_chars(iunit,varid,chardat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      CALL cf90_put_var_chars(iunit,idvars(l),chardat(l),timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write character array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_char_1d

!========================================================================

SUBROUTINE write_vars_char_2d(chardat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_char_2d* Write a 1-D character array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var_chars, error_abort, error_file, error_dim_arr,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var_chars
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
CHARACTER (LEN=*), INTENT(IN), DIMENSION(:,:) :: chardat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*chardat*  CHAR     Character output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, varflag
CHARACTER (LEN=1) :: outform
CHARACTER (LEN=12) :: clen, csize
CHARACTER (LEN=29) :: cfmt
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(chardat)) :: idvars


IF (.NOT.mastermod.OR.SIZE(chardat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_char_2d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(chardat,DIM=2)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!-------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag.AND.varid.GT.0) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      WRITE (csize,'(I12)') nosize; csize = ADJUSTL(csize)
      WRITE (clen,'(I12)') LEN(chardat(1,1)); clen = ADJUSTL(clen)
      cfmt = '('//TRIM(csize)//'A'//TRIM(clen)//')'
      WRITE (iunit,TRIM(cfmt),IOSTAT=errstat,ERR=1000) chardat
   CASE ('U')
      WRITE (iunit,IOSTAT=errstat,ERR=1000) chardat
   CASE ('N')
      CALL cf90_put_var_chars(iunit,varid,chardat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      SELECT CASE (outform)
      CASE ('A')
         IF (varflag) THEN
            WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
         ELSE
            WRITE (iunit)
         ENDIF
         WRITE (csize,'(I12)') nosize; csize = ADJUSTL(csize)
         WRITE (clen,'(I12)') LEN(chardat(1,1)); clen = ADJUSTL(clen)
         cfmt = '('//TRIM(csize)//'A'//TRIM(clen)//')'
         WRITE (iunit,TRIM(cfmt),IOSTAT=errstat,ERR=1000) chardat(:,l)
      CASE ('U')
         WRITE (iunit,IOSTAT=errstat,ERR=1000) chardat(:,l)
      CASE ('N')
         CALL cf90_put_var_chars(iunit,idvars(l),chardat(:,l),timerec)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      END SELECT
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write character array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_char_2d

!========================================================================

SUBROUTINE write_vars_int_0d(intdat,filepars,varid,varatts)
!************************************************************************
!
! *write_vars_int_0d* Write an integer scalar variable in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN) :: intdat
TYPE (FileParams), INTENT(IN) :: filepars
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER  Integer output variable
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  variable id
!*varatts*  DERIVED  Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, novars, npcc, timerec


IF (.NOT.mastermod) RETURN

procname(pglev+1) = 'write_vars_int_0d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!
!1.1 Parameters for writing
!--------------------------
!

varflag = PRESENT(varatts)
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.2 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!-------------
!

SELECT CASE (outform)
CASE ('A')
   IF (varflag) THEN
      WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
   ELSE
      WRITE (iunit,*)
   ENDIF
   WRITE (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat
CASE ('U')
   WRITE (iunit,IOSTAT=errstat,ERR=1000) intdat
CASE ('N')
   CALL cf90_put_var(iunit,varid,intdat,timerec)
   IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write integer variable'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_int_0d

!========================================================================

SUBROUTINE write_vars_int_1d(intdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_int_1d* Write a 1-D integer array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var, error_shape

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(IN), DIMENSION(:) :: intdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER  Integer output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: allrec, flag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(intdat)) :: idvars


IF (.NOT.mastermod.OR.SIZE(intdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_int_1d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(intdat)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---type of input record
allrec = outform.NE.'N'.OR.varid.NE.0

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!-------------
!
!2.1 As one record
!-----------------
!

IF (allrec) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag.AND.varid.GT.0) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      WRITE (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat
   CASE ('U')
      WRITE (iunit,IOSTAT=errstat,ERR=1000) intdat
   CASE ('N')
      CALL cf90_put_var(iunit,varid,intdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      CALL cf90_put_var(iunit,idvars(l),intdat(l),timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write integer array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_int_1d

!========================================================================

SUBROUTINE write_vars_int_2d(intdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_int_2d* Write a 2-D integer array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_file, error_dim_arr,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(IN), DIMENSION(:,:) :: intdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER  Integer output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(intdat,DIM=2)) :: idvars


IF (.NOT.mastermod.OR.SIZE(intdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_int_2d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(intdat,DIM=2)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!

iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      WRITE (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat
   CASE ('U')
      WRITE (iunit,IOSTAT=errstat,ERR=1000) intdat
   CASE ('N')
      CALL cf90_put_var(iunit,varid,intdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      SELECT CASE (outform)
      CASE ('A')
         IF (varflag) THEN
            WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
         ELSE
            WRITE (iunit)
         ENDIF
         WRITE (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat(:,l)
      CASE ('U')
         WRITE (iunit,IOSTAT=errstat,ERR=1000) intdat(:,l)
      CASE ('N')
         CALL cf90_put_var(iunit,idvars(l),intdat(:,l),timerec)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      END SELECT
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write integer array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_int_2d

!========================================================================

SUBROUTINE write_vars_int_3d(intdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_int_3d* Write a 3-D integer array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(IN), DIMENSION(:,:,:) :: intdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER  Integer output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(intdat,DIM=3)) :: idvars


IF (.NOT.mastermod.OR.SIZE(intdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_int_3d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(intdat,DIM=3)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!

iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      WRITE (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat
   CASE ('U')
      WRITE (iunit,IOSTAT=errstat,ERR=1000) intdat
   CASE ('N')
      CALL cf90_put_var(iunit,varid,intdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      SELECT CASE (outform)
      CASE ('A')
         IF (varflag) THEN
            WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
         ELSE
            WRITE (iunit,*)
         ENDIF
         WRITE (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat(:,:,l)
      CASE ('U')
         WRITE (iunit,IOSTAT=errstat,ERR=1000) intdat(:,:,l)
      CASE ('N')
         CALL cf90_put_var(iunit,idvars(l),intdat(:,:,l),timerec)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      END SELECT
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write integer array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_int_3d

!========================================================================

SUBROUTINE write_vars_int_4d(intdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_int_4d* Write a 4-D integer array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
INTEGER, INTENT(IN), DIMENSION(:,:,:,:) :: intdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*intdat*   INTEGER  Integer output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(intdat,DIM=4)) :: idvars


IF (.NOT.mastermod.OR.SIZE(intdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_int_4d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(intdat,DIM=4)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!

iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      WRITE (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat
   CASE ('U')
      WRITE (iunit,IOSTAT=errstat,ERR=1000) intdat
   CASE ('N')
      CALL cf90_put_var(iunit,varid,intdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      SELECT CASE (outform)
      CASE ('A')
         IF (varflag) THEN
            WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
         ELSE
            WRITE (iunit)
         ENDIF
         WRITE (iunit,IntegerFormat,IOSTAT=errstat,ERR=1000) intdat(:,:,:,l)
      CASE ('U')
         WRITE (iunit,IOSTAT=errstat,ERR=1000) intdat(:,:,:,l)
      CASE ('N')
         CALL cf90_put_var(iunit,idvars(l),intdat(:,:,:,l),timerec)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      END SELECT
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write integer array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_int_4d

!========================================================================

SUBROUTINE write_vars_log_0d(logdat,filepars,varid,varatts)
!************************************************************************
!
! *write_vars_log_0d* Write a logical variable in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
LOGICAL, INTENT(IN) :: logdat
TYPE (FileParams), INTENT(IN) :: filepars
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*logdat*  LOGICAL     Logical output variable
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  variable id
!*varatts*  DERIVED  Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, varflag
CHARACTER (LEN=1) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, novars, npcc, timerec


IF (.NOT.mastermod) RETURN

procname(pglev+1) = 'write_vars_log_0d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!
!1.1 Parameters for writing
!--------------------------
!

varflag = PRESENT(varatts)
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.2 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!-------------
!

SELECT CASE (outform)
CASE ('A')
   IF (varflag) THEN
      WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
   ELSE
      WRITE (iunit,*)
   ENDIF
   WRITE (iunit,LogicalFormat,IOSTAT=errstat,ERR=1000) logdat
CASE ('U')
   WRITE (iunit,IOSTAT=errstat,ERR=1000) logdat
CASE ('N')
   CALL cf90_put_var(iunit,varid,logdat,timerec)
   IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write logical variable'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_log_0d

!========================================================================

SUBROUTINE write_vars_log_1d(logdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_log_1d* Write a 1-D logical array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_dim_arr, error_file,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var, error_shape

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
LOGICAL, INTENT(IN), DIMENSION(:) :: logdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*logdat*  LOGICAL     Logical output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: allrec, flag, varflag
CHARACTER (LEN=1) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars,timerec
INTEGER, DIMENSION(SIZE(logdat)) :: idvars


IF (.NOT.mastermod.OR.SIZE(logdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_log_1d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF


!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(logdat)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---type of input record
allrec = outform.NE.'N'.OR.varid.NE.0

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!-------------
!
!2.1 As one record
!-----------------
!

IF (allrec) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag.AND.varid.GT.0) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      WRITE (iunit,LogicalFormat,IOSTAT=errstat,ERR=1000) logdat
   CASE ('U')
      WRITE (iunit,IOSTAT=errstat,ERR=1000) logdat
   CASE ('N')
      CALL cf90_put_var(iunit,varid,logdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      CALL cf90_put_var(iunit,idvars(l),logdat(l),timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   ENDDO l_220

ENDIF


CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write logical array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_log_1d

!========================================================================

SUBROUTINE write_vars_log_2d(logdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_log_2d* Write a 2-D logical array in standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_file, error_dim_arr,
!                error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
LOGICAL, INTENT(IN), DIMENSION(:,:) :: logdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*logdat*   LOGICAL  Logical output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(logdat,DIM=2)) :: idvars


IF (.NOT.mastermod.OR.SIZE(logdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_log_2d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(logdat,DIM=2)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!

iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!
!1.4 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      WRITE (iunit,LogicalFormat,IOSTAT=errstat,ERR=1000) logdat
   CASE ('U')
      WRITE (iunit,IOSTAT=errstat,ERR=1000) logdat
   CASE ('N')
      CALL cf90_put_var(iunit,varid,logdat,timerec)
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   l_220: DO l=1,nosize
      SELECT CASE (outform)
      CASE ('A')
         IF (varflag) THEN
            WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
         ELSE
            WRITE (iunit)
         ENDIF
         WRITE (iunit,LogicalFormat,IOSTAT=errstat,ERR=1000) logdat(:,l)
      CASE ('U')
         WRITE (iunit,IOSTAT=errstat,ERR=1000) logdat(:,l)
      CASE ('N')
         CALL cf90_put_var(iunit,idvars(l),logdat(:,l),timerec)
         IF (errstat.NE.noerr_NF90) GOTO 1000
      END SELECT
   ENDDO l_220

ENDIF

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write logical array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_log_2d

!========================================================================

SUBROUTINE write_vars_real_0d(realdat,filepars,varid,varatts)
!************************************************************************
!
! *write_vars_real_0d* Write a single precision real scalar variable in
!                      standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_inquire, cf90_put_var, error_abort, error_file,
!                inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
REAL (KIND=kndreal), INTENT(IN) :: realdat
TYPE (FileParams), INTENT(IN) :: filepars
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real output variable
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  variable id
!*varatts*  DERIVED  Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, outflag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, novars, npcc, timerec
REAL (KIND=kndrlong) :: realtmp


IF (.NOT.mastermod) RETURN

procname(pglev+1) = 'write_vars_real_0d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!
!1.1 Parameters for writing
!--------------------------
!
!---file attributes
varflag = PRESENT(varatts)
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether output and model data data are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   outflag = varatts(varid)%data_type.EQ.real_type
ELSE
   outflag = filepars%floattype.EQ.'S'
ENDIF

!---store to new data type
IF (.NOT.outflag) realtmp = realdat

!
!1.2 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!-------------
!

SELECT CASE (outform)
   CASE ('A')
      IF (varflag) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      IF (outflag) THEN
         WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('U')
      IF (outflag) THEN
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('N')
      IF (outflag) THEN
         CALL cf90_put_var(iunit,varid,realdat,timerec)
      ELSE
         CALL cf90_put_var(iunit,varid,realtmp,timerec)
      ENDIF
      IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write real variable'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_real_0d

!========================================================================

SUBROUTINE write_vars_real_1d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_real_1d* Write a 1-D single precision real array in standard
!                      format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndreal), INTENT(IN), DIMENSION(:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: allrec, flag, outflag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(realdat)) :: idvars
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:) :: realtmp


IF (.NOT.mastermod.OR.SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_real_1d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether output and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   outflag = varatts(varid)%data_type.EQ.real_type
ELSE
   outflag = filepars%floattype.EQ.'S'
ENDIF

!---type of input record
allrec = outform.NE.'N'.OR.varid.NE.0

!
!1.4 Allocate and store to new data type
!---------------------------------------
!

IF (.NOT.outflag) THEN
   ALLOCATE (realtmp(nosize),STAT=errstat)
   CALL error_alloc('realtmp',1,(/nosize/),rlong_type)
   realtmp = realdat
ENDIF

!
!1.5 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!-------------
!
!2.1 As one record
!-----------------
!

IF (allrec) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag.AND.varid.GT.0) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      IF (outflag) THEN
         WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('U')
      IF (outflag) THEN
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('N')
      IF (outflag) THEN
         CALL cf90_put_var(iunit,varid,realdat,timerec)
      ELSE
         CALL cf90_put_var(iunit,varid,realtmp,timerec)
      ENDIF
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (outflag) THEN
      l_221: DO l=1,nosize
         CALL cf90_put_var(iunit,idvars(l),realdat(l),timerec)
      ENDDO l_221
   ELSE
      l_222: DO l=1,nosize
         CALL cf90_put_var(iunit,idvars(l),realtmp(l),timerec)
      END DO l_222
   ENDIF
   IF (errstat.NE.noerr_NF90) GOTO 1000

ENDIF

!
!3. Deallocate
!-------------
!

IF (.NOT.outflag) DEALLOCATE (realtmp)

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write real array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_real_1d

!========================================================================

SUBROUTINE write_vars_real_2d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_real_2d* Write a 2-D single precision real array in standard
!                      format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndreal), INTENT(IN), DIMENSION(:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, outflag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(2) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=2)) :: idvars
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:,:) :: realtmp


IF (.NOT.mastermod.OR.SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_real_2d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=2)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   outflag = varatts(varid)%data_type.EQ.real_type
ELSE
   outflag = filepars%floattype.EQ.'S'
ENDIF

!
!1.4 Allocate
!------------
!

IF (.NOT.outflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2)),STAT=errstat)
   CALL error_alloc('realtmp',2,nshape,rlong_type)
ENDIF

!
!1.5 Convert to double precision
!-------------------------------
!

IF (.NOT.outflag) realtmp = realdat

!
!1.6 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (outform)
      CASE ('A')
         IF (varflag) THEN
            WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
         ELSE
            WRITE (iunit,*)
         ENDIF
         IF (outflag) THEN
            WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('U')
         IF (outflag) THEN
            WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('N')
         IF (outflag) THEN
            CALL cf90_put_var(iunit,varid,realdat,timerec)
         ELSE
            CALL cf90_put_var(iunit,varid,realtmp,timerec)
         ENDIF
         IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (outflag) THEN

      l_221: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat(:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat(:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realdat(:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp(:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realtmp(:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222

   ENDIF

ENDIF

!
!3. Deallocate
!-------------
!

IF (.NOT.outflag) DEALLOCATE (realtmp)

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write real array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_real_2d

!========================================================================

SUBROUTINE write_vars_real_3d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_real_3d* Write a 3-D single precision real array in standard
!                      format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndreal), INTENT(IN), DIMENSION(:,:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, outflag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(3) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=3)) :: idvars
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realtmp


IF (.NOT.mastermod.OR.SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_real_3d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=3)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   outflag = varatts(varid)%data_type.EQ.real_type
ELSE
   outflag = filepars%floattype.EQ.'S'
ENDIF

!
!1.4 Allocate
!------------
!

IF (.NOT.outflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2),nshape(3)),STAT=errstat)
   CALL error_alloc('realtmp',3,nshape,rlong_type)
ENDIF

!
!1.5 Convert to double precision
!-------------------------------
!

IF (.NOT.outflag) realtmp = realdat

!
!1.6 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      IF (outflag) THEN
         WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('U')
      IF (outflag) THEN
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('N')
      IF (outflag) THEN
         CALL cf90_put_var(iunit,varid,realdat,timerec)
      ELSE
         CALL cf90_put_var(iunit,varid,realtmp,timerec)
      ENDIF
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (outflag) THEN

      l_221: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat(:,:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat(:,:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realdat(:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp(:,:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realtmp(:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222

   ENDIF

ENDIF

!
!3. Deallocate
!-------------
!

IF (.NOT.outflag) DEALLOCATE (realtmp)

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write real array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_real_3d

!========================================================================

SUBROUTINE write_vars_real_4d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_real_4d* Write a 4-D single precision real array in standard
!                      format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var, error_shape

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndreal), INTENT(IN), DIMENSION(:,:,:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, outflag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(4) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=4)) :: idvars
REAL (KIND=kndrlong), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: realtmp


IF (.NOT.mastermod.OR.SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_real_4d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=4)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   outflag = varatts(varid)%data_type.EQ.real_type
ELSE
   outflag = filepars%floattype.EQ.'S'
ENDIF

!
!1.4 Allocate
!------------
!

IF (.NOT.outflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2),nshape(3),nshape(4)),STAT=errstat)
   CALL error_alloc('realtmp',4,nshape,rlong_type)
ENDIF

!
!1.5 Convert to double precision
!-------------------------------
!

IF (.NOT.outflag) realtmp = realdat

!
!1.6 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      IF (outflag) THEN
         WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('U')
      IF (outflag) THEN
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('N')
      IF (outflag) THEN
         CALL cf90_put_var(iunit,varid,realdat,timerec)
      ELSE
         CALL cf90_put_var(iunit,varid,realtmp,timerec)
      ENDIF
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension 
!------------------------------------------------
!

ELSE

   IF (outflag) THEN

      l_221: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realdat(:,:,:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat(:,:,:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realdat(:,:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) &
                                                              & realtmp(:,:,:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,:,:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realtmp(:,:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222

   ENDIF

ENDIF

!
!3. Deallocate
!-------------
!

IF (.NOT.outflag) DEALLOCATE (realtmp)

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write real array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_real_4d

!========================================================================

SUBROUTINE write_vars_double_0d(realdat,filepars,varid,varatts)
!************************************************************************
!
! *write_vars_double_0d* Write a double precision real scalar variable in
!                        standard format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_file, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_file

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
REAL (KIND=kndrlong), INTENT(IN) :: realdat
TYPE (FileParams), INTENT(IN) :: filepars
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real output variable
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  variable id
!*varatts*  DERIVED  Variable attributes
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, outflag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, novars, npcc, timerec
REAL (KIND=kndreal) :: realtmp


IF (.NOT.mastermod) RETURN

procname(pglev+1) = 'write_vars_double_0d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!
!1.1 Parameters for writing
!--------------------------
!
!---file attributes
varflag = PRESENT(varatts)
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether output and model data data are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   outflag = varatts(varid)%data_type.EQ.rlong_type
ELSE
   outflag = filepars%floattype.EQ.'D'
ENDIF

!---store to new data type
IF (.NOT.outflag) realtmp = realdat

!
!1.2 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!-------------
!

SELECT CASE (outform)
   CASE ('A')
      IF (varflag) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      IF (outflag) THEN
         WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('U')
      IF (outflag) THEN
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('N')
      IF (outflag) THEN
         CALL cf90_put_var(iunit,varid,realdat,timerec)
      ELSE
         CALL cf90_put_var(iunit,varid,realtmp,timerec)
      ENDIF
      IF (errstat.NE.noerr_NF90) GOTO 1000
END SELECT

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write real variable'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_double_0d

!========================================================================

SUBROUTINE write_vars_double_1d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_double_1d* Write a 1-D double precision real array in standard
!                        format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: allrec, flag, outflag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(SIZE(realdat)) :: idvars
REAL (KIND=kndreal), SAVE, ALLOCATABLE, DIMENSION(:) :: realtmp


IF (.NOT.mastermod.OR.SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_double_1d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether output and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   outflag = varatts(varid)%data_type.EQ.rlong_type
ELSE
   outflag = filepars%floattype.EQ.'D'
ENDIF

!---type of input record
allrec = outform.NE.'N'.OR.varid.NE.0

!
!1.4 Allocate and store to new data type
!---------------------------------------
!

IF (.NOT.outflag) THEN
   ALLOCATE (realtmp(nosize),STAT=errstat)
   CALL error_alloc('realtmp',1,(/nosize/),real_type)
   realtmp = realdat
ENDIF

!
!1.5 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!-------------
!
!2.1 As one record
!-----------------
!

IF (allrec) THEN

   SELECT CASE (outform)
   CASE ('A')
      IF (varflag.AND.varid.GT.0) THEN
         WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
      ELSE
         WRITE (iunit,*)
      ENDIF
      IF (outflag) THEN
         WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('U')
      IF (outflag) THEN
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat
      ELSE
         WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp
      ENDIF
   CASE ('N')
      IF (outflag) THEN
         CALL cf90_put_var(iunit,varid,realdat,timerec)
      ELSE
         CALL cf90_put_var(iunit,varid,realtmp,timerec)
      ENDIF
      IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (outflag) THEN
      l_221: DO l=1,nosize
         CALL cf90_put_var(iunit,idvars(l),realdat(l),timerec)
      ENDDO l_221
   ELSE
      l_222: DO l=1,nosize
         CALL cf90_put_var(iunit,idvars(l),realtmp(l),timerec)
      END DO l_222
   ENDIF
   IF (errstat.NE.noerr_NF90) GOTO 1000
   
ENDIF

CALL log_timer_out(npcc,itm_output)

!
!3. Deallocate
!-------------
!

IF (.NOT.outflag) DEALLOCATE (realtmp)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write real array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_double_1d

!========================================================================

SUBROUTINE write_vars_double_2d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_double_2d* Write a 2-D double precision real array in standard
!                        format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, outflag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(2) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=2)) :: idvars
REAL (KIND=kndreal), SAVE, ALLOCATABLE, DIMENSION(:,:) :: realtmp


IF (.NOT.mastermod.OR.SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_double_2d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=2)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   outflag = varatts(varid)%data_type.EQ.rlong_type
ELSE
   outflag = filepars%floattype.EQ.'D'
ENDIF

!
!1.4 Allocate and store to new data type
!---------------------------------------
!

IF (.NOT.outflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2)),STAT=errstat)
   CALL error_alloc('realtmp',2,nshape,real_type)
ENDIF

!
!1.5 Convert to single precision
!-------------------------------
!

IF (.NOT.outflag) realtmp = realdat

!
!1.6 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (outform)
      CASE ('A')
         IF (varflag) THEN
            WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
         ELSE
            WRITE (iunit,*)
         ENDIF
         IF (outflag) THEN
            WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('U')
         IF (outflag) THEN
            WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('N')
         IF (outflag) THEN
            CALL cf90_put_var(iunit,varid,realdat,timerec)
         ELSE
            CALL cf90_put_var(iunit,varid,realtmp,timerec)
         ENDIF
         IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (outflag) THEN

      l_221: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat(:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat(:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realdat(:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp(:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realtmp(:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222

   ENDIF

ENDIF

!
!3. Deallocate
!-------------
!

IF (.NOT.outflag) DEALLOCATE (realtmp)

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write real array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_double_2d

!========================================================================

SUBROUTINE write_vars_double_3d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_double_3d* Write a 3-D real double precision array in standard
!                        format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(:,:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, outflag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(3) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=3)) :: idvars
REAL (KIND=kndreal), SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: realtmp


IF (.NOT.mastermod.OR.SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_double_3d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=3)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   outflag = varatts(varid)%data_type.EQ.rlong_type
ELSE
   outflag = filepars%floattype.EQ.'D'
ENDIF

!
!1.4 Allocate and store to new data type
!---------------------------------------
!

IF (.NOT.outflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2),nshape(3)),STAT=errstat)
   CALL error_alloc('realtmp',3,nshape,real_type)
ENDIF

!
!1.5 Convert to single precision
!-------------------------------
!

IF (.NOT.outflag) realtmp = realdat

!
!1.6 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (outform)
      CASE ('A')
         IF (varflag) THEN
            WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
         ELSE
            WRITE (iunit,*)
         ENDIF
         IF (outflag) THEN
            WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('U')
         IF (outflag) THEN
            WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('N')
         IF (outflag) THEN
            CALL cf90_put_var(iunit,varid,realdat,timerec)
         ELSE
            CALL cf90_put_var(iunit,varid,realtmp,timerec)
         ENDIF
         IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (outflag) THEN

      l_221: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat(:,:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat(:,:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realdat(:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp(:,:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realtmp(:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222

   ENDIF

ENDIF

!
!3. Deallocate
!-------------
!

IF (.NOT.outflag) DEALLOCATE (realtmp)


CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write real array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_double_3d

!========================================================================

SUBROUTINE write_vars_double_4d(realdat,filepars,varid,varatts,vecids)
!************************************************************************
!
! *write_vars_real_4d* Write a 4-D real double precision array in standard
!                      format
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - cf90_put_var, error_abort, error_alloc, error_dim_arr,
!                error_file, error_limits_arr, error_limits_var, inquire_unit
!
!************************************************************************
!
USE datatypes
USE paralpars
USE cf90_routines, ONLY: cf90_put_var
USE error_routines, ONLY: error_abort, error_alloc, error_dim_arr, error_file, &
                        & error_limits_arr, error_limits_var

!
!*Arguments
!
INTEGER, INTENT(IN) :: varid
TYPE (FileParams), INTENT(IN) :: filepars
INTEGER, INTENT(IN), OPTIONAL, DIMENSION(:) :: vecids
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(:,:,:,:) :: realdat
TYPE (VariableAtts), INTENT(IN), OPTIONAL, DIMENSION(:) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*realdat*  REAL     Real output array
!*filepars* DERIVED  File attributes
!*varid*    INTEGER  If >0, variable id of output array.
!                    If =0, last dimension is a variable dimension
!*varatts*  DERIVED  Variable attributes
!*vecids*   INTEGER  Arrays of variable IDs if varid = 0
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: flag, outflag, varflag
CHARACTER (LEN=lenform) :: outform
CHARACTER (LEN=leniofile) :: filename
INTEGER :: iunit, ivar, l, nosize, novars, npcc, numvars, timerec
INTEGER, DIMENSION(4) :: nshape
INTEGER, DIMENSION(SIZE(realdat,DIM=4)) :: idvars
REAL (KIND=kndreal), SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) :: realtmp


IF (.NOT.mastermod.OR.SIZE(realdat).EQ.0) RETURN

procname(pglev+1) = 'write_vars_double_4d'

IF (.NOT.PRESENT(varatts).OR.varid.EQ.0) THEN
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename))
ELSE
   CALL log_timer_in(npcc,&
        & logname=TRIM(procname(pglev+1))//': '//TRIM(filepars%filename)//&
        & ':'//TRIM(varatts(varid)%f90_name))
ENDIF

!
!1. Check and initialise
!-----------------------
!

nosize = SIZE(realdat,DIM=4)
varflag = PRESENT(varatts)

!
!1.1 Check variables ids
!-----------------------
!

IF (errchk) THEN
   numvars = filepars%nocoords + filepars%novars   
   IF (varid.GT.0) THEN
      CALL error_limits_var(varid,'varid',1,numvars)
   ELSEIF (PRESENT(vecids)) THEN
      CALL error_dim_arr(SIZE(vecids),'vecids',nosize)
      ivar_110: DO ivar=1,nosize
         CALL error_limits_arr(vecids(ivar),'vecids',1,numvars,1,(/ivar/))
      ENDDO ivar_110
   ELSE
      CALL error_limits_var(nosize,'nosize',1,numvars)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_arg)
ENDIF

!
!1.2 Initialise variable ids
!---------------------------
!

IF (PRESENT(vecids)) THEN
   idvars = vecids
ELSEIF (varid.EQ.0) THEN
   idvars = (/(ivar,ivar=1,nosize)/)
ENDIF

!
!1.3 Parameters for writing
!--------------------------
!
!---file attributes
iunit = filepars%iunit
outform = filepars%form
filename = filepars%filename
novars = filepars%novars
timerec = filepars%timerec

!---check whether input and model data type are the same
IF (PRESENT(varatts).AND.varid.GT.0) THEN
   outflag = varatts(varid)%data_type.EQ.rlong_type
ELSE
   outflag = filepars%floattype.EQ.'D'
ENDIF

!
!1.4 Allocate and store to new data type
!---------------------------------------
!

IF (.NOT.outflag) THEN
   nshape = SHAPE(realdat)
   ALLOCATE (realtmp(nshape(1),nshape(2),nshape(3),nshape(4)),STAT=errstat)
   CALL error_alloc('realtmp',4,nshape,real_type)
ENDIF

!
!1.5 Convert to single precision
!-------------------------------
!

IF (.NOT.outflag) realtmp = realdat

!
!1.6 Check whether file unit is opened
!-------------------------------------
!

flag = inquire_unit(iunit,outform)
IF (.NOT.flag) CALL error_file(ierrno_fopen,filepars=filepars)

!
!2. Write data
!------------
!
!2.1 As one record
!-----------------
!

IF (varid.GT.0) THEN

   SELECT CASE (outform)
      CASE ('A')
         IF (varflag) THEN
            WRITE (iunit,'(A)') TRIM(varatts(varid)%f90_name)
         ELSE
            WRITE (iunit,*)
         ENDIF
         IF (outflag) THEN
            WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('U')
         IF (outflag) THEN
            WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat
         ELSE
            WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp
         ENDIF
      CASE ('N')
         IF (outflag) THEN
            CALL cf90_put_var(iunit,varid,realdat,timerec)
         ELSE
            CALL cf90_put_var(iunit,varid,realtmp,timerec)
         ENDIF
         IF (errstat.NE.noerr_NF90) GOTO 1000
   END SELECT

!
!2.2 Multiple records specified by last dimension
!------------------------------------------------
!

ELSE

   IF (outflag) THEN

      l_221: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,DoubleFormat,IOSTAT=errstat,ERR=1000) &
                                                              & realdat(:,:,:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realdat(:,:,:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realdat(:,:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_221

   ELSE

      l_222: DO l=1,nosize
         SELECT CASE (outform)
            CASE ('A')
               IF (varflag) THEN
                  WRITE (iunit,'(A)') TRIM(varatts(l)%f90_name)
               ELSE
                  WRITE (iunit,*)
               ENDIF
               WRITE (iunit,RealFormat,IOSTAT=errstat,ERR=1000) realtmp(:,:,:,l)
            CASE ('U')
               WRITE (iunit,IOSTAT=errstat,ERR=1000) realtmp(:,:,:,l)
            CASE ('N')
               CALL cf90_put_var(iunit,idvars(l),realtmp(:,:,:,l),timerec)
               IF (errstat.NE.noerr_NF90) GOTO 1000
         END SELECT
      ENDDO l_222

   ENDIF

ENDIF

!
!3. Deallocate
!-------------
!

IF (.NOT.outflag) DEALLOCATE (realtmp)

CALL log_timer_out(npcc,itm_output)


RETURN

!---process write error
1000 CONTINUE
IF (errchk) WRITE (ioerr,'(A)') 'Unable to write real array'
CALL error_file(ierrno_write,filepars=filepars)

END SUBROUTINE write_vars_double_4d


END MODULE inout_routines
