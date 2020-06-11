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

MODULE cif_routines
!************************************************************************
!
! *cif_routines* Routines for reading and converting data from a central
!                input file
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11.2
!
! $Date: 2018-05-31 16:41:05 +0200 (Thu, 31 May 2018) $
!
! $Revision: 1141 $
!
! Description -
!
! Generic routines - conv_from_chars, conv_to_chars
!
! Subroutines - check_cif_lbound_novars, conv_from_chars_gridpars,
!               conv_from_chars_infiles, conv_from_chars_modvars,
!               conv_from_chars_outfiles, conv_from_chars_outvars,
!               conv_from_chars_statlocs, conv_from_chars_surfgrd,
!               conv_to_chars_gridpars, conv_to_chars_infiles,
!               conv_to_chars_modvars, conv_to_chars_outfiles,
!               conv_to_chars_outvars, conv_to_chars_statlocs,
!               conv_to_chars_surfgrd, error_cif_var, read_cif_line,
!               write_cif_line
!
!************************************************************************
!

INTERFACE conv_from_chars
   MODULE PROCEDURE conv_from_chars_chars, conv_from_chars_int, &
        & conv_from_chars_log, conv_from_chars_real, &
        & conv_from_chars_double
END INTERFACE

INTERFACE conv_to_chars
   MODULE PROCEDURE conv_to_chars_chars_0d, conv_to_chars_int_0d, &
                  & conv_to_chars_int_1d, conv_to_chars_log_0d, &
                  & conv_to_chars_log_1d, conv_to_chars_real_0d, &
                  & conv_to_chars_real_1d, conv_to_chars_double_0d, &
                  & conv_to_chars_double_1d
END INTERFACE


CONTAINS

!========================================================================

SUBROUTINE check_cif_lbound_novars(novars,novarsmin)
!************************************************************************
!
! *check_novars_lread* Check whether the number of parameters on the input line
!                      of a CIF file is not lower than a minimum limit
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - displays error message if needed
!
! Module calls - error_abort
!
!************************************************************************
!
USE iopars
USE paralpars  
USE error_routines, ONLY: nerrs, error_abort

!
!*Arguments
!
INTEGER, INTENT(IN) :: novars, novarsmin

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*novars*    INTEGER Number of parameters on input line
!*novarsmax* INTEGER Minimum number of input parameters required on input line
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: clnum, cvars, cvarsmin


IF (novars.LT.novarsmin) THEN
   nerrs = nerrs + 1
   IF (master.AND.errchk.AND.nerrs.LE.maxerrors) THEN
      WRITE (cvars,'(I12)') novars; cvars = ADJUSTL(cvars)
      WRITE (cvarsmin,'(I12)') novarsmin; cvarsmin = ADJUSTL(cvarsmin)
      WRITE (clnum,'(I12)') ciflinenum; clnum = ADJUSTL(clnum)
      WRITE (ioerr,'(A)') 'Wrong number of input parameters on line '//&
                 & TRIM(clnum)//' of file '&
                 & //TRIM(ciffile%filename)//': '//TRIM(cvars)
      WRITE (ioerr,'(A)') 'Should be no lower than: '//TRIM(cvarsmin)
   ENDIF
   CALL error_abort(procname(pglev),ierrno_inival)
ENDIF


RETURN

END SUBROUTINE check_cif_lbound_novars

!========================================================================

SUBROUTINE conv_from_chars_chars(string,cdat,lvar)
!************************************************************************
!
! *conv_from_chars_chars* Assign the value of 'string' to 'cvar' if 'string' is
!                         non-empty
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
! Module calls - error_cif_var
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: string
CHARACTER (LEN=*), INTENT(OUT) :: cdat
INTEGER, INTENT(IN) :: lvar

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Input string
!*cdat*     CHAR    Output string variable
!*lvar*     INTEGER Position of variable on input line
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: l1, l2


l1 = LEN_TRIM(string); l2 = LEN(cdat)

IF (l1.GT.0) THEN
   IF (l1.LE.l2) THEN
      cdat = TRIM(string)
   ELSE
      CALL error_cif_var(lvar)
   ENDIF
ENDIF


RETURN

END SUBROUTINE conv_from_chars_chars

!========================================================================

SUBROUTINE conv_from_chars_int(string,idat,lvar)
!************************************************************************
!
! *conv_from_chars_int* Convert a numerical string into an integer variable
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
! Module calls - error_cif_var
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: lvar
INTEGER, INTENT(OUT) :: idat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Input string
!*idat*     INTEGER Output integer variable
!*lvar*     INTEGER Position of variable on input line
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ierr


IF (LEN_TRIM(string).GT.0) THEN
   READ (string,*,IOSTAT=ierr) idat
   IF (ierr.NE.0) CALL error_cif_var(lvar)
ENDIF 


RETURN

END SUBROUTINE conv_from_chars_int

!========================================================================

SUBROUTINE conv_from_chars_log(string,ldat,lvar)
!************************************************************************
!
! *conv_from_chars_log* Convert a numerical string into a logical variable
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
! Module calls - error_cif_var
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: lvar
LOGICAL, INTENT(OUT) :: ldat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Input string
!*ldat*     LOGICAL Output logical variable
!*lvar*     INTEGER Position of variable on input line
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ierr


IF (LEN_TRIM(string).GT.0) THEN
   READ (string,*,IOSTAT=ierr) ldat
   IF (ierr.NE.0) CALL error_cif_var(lvar)
ENDIF 


RETURN

END SUBROUTINE conv_from_chars_log

!========================================================================

SUBROUTINE conv_from_chars_real(string,rdat,lvar)
!************************************************************************
!
! *conv_from_chars_real* Convert a numerical string into a single precision
!                        real variable
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
! Module calls - error_cif_var
!
!************************************************************************
!  
USE syspars

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: lvar
REAL (KIND=kndreal), INTENT(OUT) :: rdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Input string
!*rdat*     REAL    Output real variable
!*lvar*     INTEGER Position of variable on input line
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ierr


IF (LEN_TRIM(string).GT.0) THEN
   READ (string,*,IOSTAT=ierr) rdat
   IF (ierr.NE.0) CALL error_cif_var(lvar)
ENDIF


RETURN

END SUBROUTINE conv_from_chars_real

!========================================================================

SUBROUTINE conv_from_chars_double(string,rdat,lvar)
!************************************************************************
!
! *conv_from_chars_double* Convert a numerical string into a double precision
!                          real variable
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
! Module calls - error_cif_var
!
!************************************************************************
!
USE syspars 

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: lvar
REAL (KIND=kndrlong), INTENT(OUT) :: rdat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Input string
!*rdat*     REAL    Output real variable
!*lvar*     INTEGER Position of variable on input line
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: ierr


IF (LEN_TRIM(string).GT.0) THEN
   READ (string,*,IOSTAT=ierr) rdat
   IF (ierr.NE.0) CALL error_cif_var(lvar)
ENDIF


RETURN

END SUBROUTINE conv_from_chars_double

!========================================================================

SUBROUTINE conv_from_chars_gridpars(cvals,gridpars,lvar)
!************************************************************************
!
! *conv_from_chars_gridpars* Obtain output grid attributes from a array of 
!                            character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_from_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), DIMENSION(19) :: cvals
INTEGER, INTENT(INOUT) :: lvar
TYPE (OutGridParams), INTENT(INOUT) :: gridpars

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*gridpars* DERIVED  Returned output grid attributes
!*cvals*    CHAR     Data values on input data line     
!*lvar*     INTEGER  Position of variable on input line
!
!------------------------------------------------------------------------------
!


lvar = lvar + 1
CALL conv_from_chars(cvals(1),gridpars%gridded,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(2),gridpars%packing,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(3),gridpars%time_grid,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(4),gridpars%refdate,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(5),gridpars%nostats,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(6),gridpars%time_format,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(7),gridpars%vcoord,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(8),gridpars%xlims(1),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(9),gridpars%xlims(2),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(10),gridpars%xlims(3),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(11),gridpars%ylims(1),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(12),gridpars%ylims(2),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(13),gridpars%ylims(3),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(14),gridpars%zlims(1),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(15),gridpars%zlims(2),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(16),gridpars%zlims(3),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(17),gridpars%tlims(1),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(18),gridpars%tlims(2),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(19),gridpars%tlims(3),lvar)


RETURN

END SUBROUTINE conv_from_chars_gridpars

!========================================================================

SUBROUTINE conv_from_chars_infiles(cvals,filepars,lvar)
!************************************************************************
!
! *conv_from_chars_infiles* Obtain forcing file attributes from a array of 
!                           character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_from_chars
!
!***********************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), DIMENSION(12) :: cvals
INTEGER, INTENT(INOUT) :: lvar
TYPE (FileParams), INTENT(INOUT) :: filepars

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*filepars* DERIVED  Returned file attributes
!*cvals*    CHAR     Data values on input data line     
!*lvar*     INTEGER  Position of variable on input line
!
!------------------------------------------------------------------------------
!


lvar = lvar + 1
CALL conv_from_chars(cvals(1),filepars%form,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(2),filepars%status,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(3),filepars%filename,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(4),filepars%floattype,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(5),filepars%tlims(1),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(6),filepars%tlims(2),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(7),filepars%tlims(3),lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(8),filepars%packing,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(9),filepars%endfile,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(10),filepars%info,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(11),filepars%time_regular,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(12),filepars%pathname,lvar)


RETURN

END SUBROUTINE conv_from_chars_infiles

!========================================================================

SUBROUTINE conv_from_chars_modvars(cvals,outvars,lvar)
!************************************************************************
!
! *conv_from_chars_modvars* Obtain attributes of model output variables from
!                           an array of  character strings
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_from_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), DIMENSION(:) :: cvals
INTEGER, INTENT(INOUT) :: lvar
TYPE (VariableAtts), INTENT(INOUT) :: outvars

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*outvars*  DERIVED  Returned variable attributes
!*cvals*    CHAR     Data values on input data line     
!*lvar*     INTEGER  Position of variable on input line
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: l


lvar = lvar + 1
CALL conv_from_chars(cvals(1),outvars%data_type,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(2),outvars%nrank,lvar)
ivar_110: DO l=1,outvars%nrank
   lvar = lvar + 1
   CALL conv_from_chars(cvals(l+2),outvars%global_dims(l),lvar)
ENDDO ivar_110
l = outvars%nrank + 2
lvar = lvar + 1
CALL conv_from_chars(cvals(l+1),outvars%f90_name,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(l+2),outvars%standard_name,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(l+3),outvars%comment,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(l+4),outvars%long_name,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(l+5),outvars%units,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(l+6),outvars%vector_name,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(l+7),outvars%fill,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(l+8),outvars%fill_value,lvar)


RETURN

END SUBROUTINE conv_from_chars_modvars

!========================================================================

SUBROUTINE conv_from_chars_outfiles(cvals,filepars,lvar)
!************************************************************************
!
! *conv_from_chars_outfiles* Obtain user output file attributes from a array of 
!                            character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_from_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), DIMENSION(6) :: cvals
INTEGER, INTENT(INOUT) :: lvar
TYPE (FileParams), INTENT(INOUT) :: filepars

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*filepars* DERIVED  Returned file attributes
!*cvals*    CHAR     Data values on input data line     
!*lvar*     INTEGER  Position of variable on input line
!
!------------------------------------------------------------------------------
!


lvar = lvar + 1
CALL conv_from_chars(cvals(1),filepars%defined,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(2),filepars%form,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(3),filepars%filename,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(4),filepars%floattype,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(5),filepars%info,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(6),filepars%pathname,lvar)


RETURN

END SUBROUTINE conv_from_chars_outfiles

!========================================================================

SUBROUTINE conv_from_chars_outvars(cvals,varatts,lvar)
!************************************************************************
!
! *conv_from_chars_outvars* Obtain variable attributes from a array of 
!                           character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_from_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), DIMENSION(11) :: cvals
INTEGER, INTENT(INOUT) :: lvar
TYPE (VariableAtts), INTENT(INOUT) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*varatts*  DERIVED  Returned variable attributes
!*cvals*    CHAR     Data values on input data line     
!*lvar*     INTEGER  Position of variable on input line
!
!------------------------------------------------------------------------------
!


lvar = lvar + 1
CALL conv_from_chars(cvals(1),varatts%ivarid,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(2),varatts%nrank,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(3),varatts%numvar,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(4),varatts%oopt,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(5),varatts%klev,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(6),varatts%dep,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(7),varatts%node,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(8),varatts%f90_name,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(9),varatts%long_name,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(10),varatts%units,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(11),varatts%vector_name,lvar)


RETURN

END SUBROUTINE conv_from_chars_outvars

!========================================================================

SUBROUTINE conv_from_chars_statlocs(cvals,statlocs,lvar)
!************************************************************************
!
! *conv_from_chars_statlocs* Obtain output station attributes from a array of 
!                            character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_from_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), DIMENSION(3) :: cvals
INTEGER, INTENT(INOUT) :: lvar
TYPE (StationLocs), INTENT(INOUT) :: statlocs

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*statlocs* DERIVED  Returned output station attributes
!*cvals*    CHAR     Data values on input data line     
!*lvar*     INTEGER  Position of variable on input line
!
!------------------------------------------------------------------------------
!


lvar = lvar + 1
CALL conv_from_chars(cvals(1),statlocs%ipos,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(2),statlocs%jpos,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(3),statlocs%name,lvar)


RETURN

END SUBROUTINE conv_from_chars_statlocs

!========================================================================

SUBROUTINE conv_from_chars_surfgrd(cvals,surfgrd,lvar)
!************************************************************************
!
! *conv_from_chars_surfgrd* Obtain surface grid attributes from a array of 
!                           character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11.2
!
! Description -
!
! Module calls - conv_from_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN), DIMENSION(14) :: cvals
INTEGER, INTENT(INOUT) :: lvar
TYPE (GridParams), INTENT(INOUT) :: surfgrd

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*surfgrd*  DERIVED  Returned surface grid attributes
!*cvals*    CHAR     Data values on input data line     
!*lvar*     INTEGER  Position of variable on input line
!
!------------------------------------------------------------------------------
!


lvar = lvar + 1
CALL conv_from_chars(cvals(1),surfgrd%nhtype,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(2),surfgrd%n1dat,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(3),surfgrd%n2dat,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(4),surfgrd%surcoords,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(5),surfgrd%delxdat,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(6),surfgrd%delydat,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(7),surfgrd%x0dat,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(8),surfgrd%y0dat,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(9),surfgrd%rotated,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(10),surfgrd%gridangle,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(11),surfgrd%y0rot,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(12),surfgrd%extrapol,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(13),surfgrd%land,lvar)
lvar = lvar + 1
CALL conv_from_chars(cvals(14),surfgrd%datflag,lvar)


RETURN

END SUBROUTINE conv_from_chars_surfgrd

!========================================================================

SUBROUTINE conv_to_chars_chars_0d(string,cvar)
!************************************************************************
!
! *conv_to_chars_chars_0d* Assign the string 'cdat' to the string 'string'
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT) :: string
CHARACTER (LEN=*), INTENT(IN) :: cvar

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Output string
!*cvar*     CHAR    Character variable
!
!------------------------------------------------------------------------------
!


IF (TRIM(cvar).NE.'') THEN
   string = cvar
ELSE
   string = ''
ENDIF


RETURN

END SUBROUTINE conv_to_chars_chars_0d

!========================================================================

SUBROUTINE conv_to_chars_int_0d(string,ivar)
!************************************************************************
!
! *conv_to_chars_int_0d* Convert an integer scalar variable into a string
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT) :: string
INTEGER, INTENT(IN) :: ivar

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Output string
!*ivar*     INTEGER Input integer variable
!
!------------------------------------------------------------------------------
!


string = ''
WRITE(string(1:12),'(I12)') ivar
string = ADJUSTL(string)


RETURN

END SUBROUTINE conv_to_chars_int_0d

!========================================================================

SUBROUTINE conv_to_chars_int_1d(string,ivar)
!************************************************************************
!
! *conv_to_chars_int_1d* Convert an integer vector into a string array
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(:) :: string
INTEGER, INTENT(IN), DIMENSION(:) :: ivar

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Output string array
!*ivar*     INTEGER Input integer vector
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize


nsize = SIZE(ivar)

n_110: DO n=1,nsize
   string(n) = ''
   WRITE (string(n)(1:12),'(I12)') ivar(n)
   string(n) = ADJUSTL(string(n))
ENDDO n_110

RETURN

END SUBROUTINE conv_to_chars_int_1d

!========================================================================

SUBROUTINE conv_to_chars_log_0d(string,lvar)
!************************************************************************
!
! *conv_to_chars_log_0d* Convert a logical scalar variable into a string
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT) :: string
LOGICAL, INTENT(IN) :: lvar

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Output string
!*lvar*     LOGICAL Input logical variable
!
!------------------------------------------------------------------------------
!


string = ''
WRITE(string(1:12),'(L1)') lvar
string = ADJUSTL(string)


RETURN

END SUBROUTINE conv_to_chars_log_0d

!========================================================================

SUBROUTINE conv_to_chars_log_1d(string,lvar)
!************************************************************************
!
! *conv_to_chars_log_1d* Convert a logical vector into a string array
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
!************************************************************************
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(:) :: string
LOGICAL, INTENT(IN), DIMENSION(:) :: lvar

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Output string array
!*lvar*     LOGICAL Input logical vector
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: n, nsize


nsize = SIZE(lvar)

n_110: DO n=1,nsize
   string(n) = ''
   WRITE (string(n)(1:12),'(L1)') lvar(n)
   string(n) = ADJUSTL(string(n))
ENDDO n_110

RETURN

END SUBROUTINE conv_to_chars_log_1d

!========================================================================

SUBROUTINE conv_to_chars_real_0d(string,rvar)
!************************************************************************
!
! *conv_to_chars_real_0d* Convert a single precision real scalar variable into
!                         a string
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
!************************************************************************
!
USE syspars
  
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT) :: string
REAL (KIND=kndreal), INTENT(IN) :: rvar

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Output string
!*rvar*     REAL    Input real variable
!
!------------------------------------------------------------------------------
!
!
!* Local variables
!
INTEGER :: l, lposE


string = ''
WRITE(string(1:15),'(G15.7)') rvar

lposE = INDEX(string,'E')
l = MERGE(LEN_TRIM(string),lposE-1,lposE.EQ.0)
DO WHILE(string(l:l).EQ.'0')
   string(l:l) = ''
   l = l-1
ENDDO
IF (lposE.GT.0) string = string(1:l)//string(lposE:LEN_TRIM(string))

string = ADJUSTL(string)


RETURN

END SUBROUTINE conv_to_chars_real_0d

!========================================================================

SUBROUTINE conv_to_chars_real_1d(string,rvar)
!************************************************************************
!
! *conv_to_chars_real_1d* Convert a single precision real vector into a string
!                         array
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
!************************************************************************
!
USE syspars  

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(:) :: string
REAL (KIND=kndreal), INTENT(IN), DIMENSION(:) :: rvar

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Output string array
!*rvar*     REAL    Input real vector
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: l, lposE, n, nsize


nsize = SIZE(rvar)

n_110: DO n=1,nsize

   string(n) = ''
   WRITE (string(n)(1:15),'(G15.7)') rvar(n)

   lposE = INDEX(string(n),'E')
   l = MERGE(LEN_TRIM(string(n)),lposE-1,lposE.EQ.0)
   DO WHILE(string(n)(l:l).EQ.'0')
      string(n)(l:l) = ''
      l = l-1
   ENDDO
   IF (lposE.GT.0) THEN
      string(n) = string(n)(1:l)//string(n)(lposE:LEN_TRIM(string(n)))
   ENDIF
   string(n) = ADJUSTL(string(n))

ENDDO n_110


RETURN

END SUBROUTINE conv_to_chars_real_1d

!========================================================================

SUBROUTINE conv_to_chars_double_0d(string,rvar)
!************************************************************************
!
! *conv_to_chars_double_0d* Convert a double precision real scalar variable
!                           into a string
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
!************************************************************************
!
USE syspars

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT) :: string
REAL (KIND=kndrlong), INTENT(IN) :: rvar

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Output string
!*rvar*     REAL    Input real variable
!
!------------------------------------------------------------------------------
!
!
!* Local variables
!
INTEGER :: l, lposE


string = ''
WRITE(string(1:22),'(G22.14)') rvar

lposE = INDEX(string,'E')
l = MERGE(LEN_TRIM(string),lposE-1,lposE.EQ.0)
DO WHILE(string(l:l).EQ.'0')
   string(l:l) = ''
   l = l-1
ENDDO
IF (lposE.GT.0) string = string(1:l)//string(lposE:LEN_TRIM(string))

string = ADJUSTL(string)


RETURN

END SUBROUTINE conv_to_chars_double_0d

!========================================================================

SUBROUTINE conv_to_chars_double_1d(string,rvar)
!************************************************************************
!
! *conv_to_chars_double_1d* Convert a double precision real vector into a
!                           string array
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
!************************************************************************
!
USE syspars
  
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(:) :: string
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(:) :: rvar

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*string*   CHAR    Output string array
!*rvar*     REAL    Input real vector
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: l, lposE, n, nsize


nsize = SIZE(rvar)

n_110: DO n=1,nsize

   string(n) = ''
   WRITE (string(n)(1:22),'(G22.14)') rvar(n)

   lposE = INDEX(string(n),'E')
   l = MERGE(LEN_TRIM(string(n)),lposE-1,lposE.EQ.0)
   DO WHILE(string(n)(l:l).EQ.'0')
      string(n)(l:l) = ''
      l = l-1
   ENDDO
   IF (lposE.GT.0) THEN
      string(n) = string(n)(1:l)//string(n)(lposE:LEN_TRIM(string(n)))
   ENDIF
   string(n) = ADJUSTL(string(n))

ENDDO n_110


RETURN

END SUBROUTINE conv_to_chars_double_1d

!========================================================================

SUBROUTINE conv_to_chars_gridpars(cvals,gridpars)
!************************************************************************
!
! *conv_to_chars_gridpars* Write output grid attributes to an array of 
!                          character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_to_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(19) :: cvals
TYPE (OutGridParams), INTENT(IN) :: gridpars

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*gridpars* DERIVED  Output grid attributes
!*cvals*    CHAR     Returned data string array
!
!------------------------------------------------------------------------------
!


CALL conv_to_chars(cvals(1),gridpars%gridded)
CALL conv_to_chars(cvals(2),gridpars%packing)
CALL conv_to_chars(cvals(3),gridpars%time_grid)
cvals(4) = gridpars%refdate
CALL conv_to_chars(cvals(5),gridpars%nostats)
CALL conv_to_chars(cvals(6),gridpars%time_format)
CALL conv_to_chars(cvals(7),gridpars%vcoord)
CALL conv_to_chars(cvals(8),gridpars%xlims(1))
CALL conv_to_chars(cvals(9),gridpars%xlims(2))
CALL conv_to_chars(cvals(10),gridpars%xlims(3))
CALL conv_to_chars(cvals(11),gridpars%ylims(1))
CALL conv_to_chars(cvals(12),gridpars%ylims(2))
CALL conv_to_chars(cvals(13),gridpars%ylims(3))
CALL conv_to_chars(cvals(14),gridpars%zlims(1))
CALL conv_to_chars(cvals(15),gridpars%zlims(2))
CALL conv_to_chars(cvals(16),gridpars%zlims(3))
CALL conv_to_chars(cvals(17),gridpars%tlims(1))
CALL conv_to_chars(cvals(18),gridpars%tlims(2))
CALL conv_to_chars(cvals(19),gridpars%tlims(3))


RETURN

END SUBROUTINE conv_to_chars_gridpars

!========================================================================

SUBROUTINE conv_to_chars_infiles(cvals,filepars)
!************************************************************************
!
! *conv_to_chars_infiles* Write forcing file attributes to an array of 
!                         character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_to_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(12) :: cvals
TYPE (FileParams), INTENT(IN) :: filepars

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*filepars* DERIVED  File attributes
!*cvals*    CHAR     Returned data string array
!
!------------------------------------------------------------------------------
!


cvals(1) = filepars%form
cvals(2) = filepars%status
cvals(3) = filepars%filename
cvals(4) = filepars%floattype
CALL conv_to_chars(cvals(5),filepars%tlims(1))
CALL conv_to_chars(cvals(6),filepars%tlims(2))
CALL conv_to_chars(cvals(7),filepars%tlims(3))
CALL conv_to_chars(cvals(8),filepars%packing)
CALL conv_to_chars(cvals(9),filepars%endfile)
CALL conv_to_chars(cvals(10),filepars%info)
CALL conv_to_chars(cvals(11),filepars%time_regular)
cvals(12) = filepars%pathname


RETURN

END SUBROUTINE conv_to_chars_infiles

!========================================================================

SUBROUTINE conv_to_chars_modvars(cvals,outvars)
!************************************************************************
!
! *conv_to_chars_modvars* Write attributes of model output variables to an array
!                         of character strings
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_to_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(:) :: cvals
TYPE (VariableAtts), INTENT(IN) :: outvars

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*varatts*  DERIVED  Variable attributes
!*cvals*    CHAR     Returned data string array
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: l


CALL conv_to_chars(cvals(1),outvars%data_type)
CALL conv_to_chars(cvals(2),outvars%nrank)
ivar_110: DO l=1,outvars%nrank
   CALL conv_to_chars(cvals(l+2),outvars%global_dims(l))
ENDDO ivar_110
l = outvars%nrank + 2
cvals(l+1) = outvars%f90_name
cvals(l+2) = outvars%standard_name
cvals(l+3) = outvars%comment
cvals(l+4) = outvars%long_name
cvals(l+5) = outvars%units
cvals(l+6) = outvars%vector_name
CALL conv_to_chars(cvals(l+7),outvars%fill)
IF (outvars%fill) THEN
   CALL conv_to_chars(cvals(l+8),outvars%fill_value)
ELSE
   cvals(l+8) = ''
ENDIF


RETURN

END SUBROUTINE conv_to_chars_modvars

!========================================================================

SUBROUTINE conv_to_chars_outfiles(cvals,filepars)
!************************************************************************
!
! *conv_to_chars_outfiles* Write user output file attributes to an array of 
!                          character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_to_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(6) :: cvals
TYPE (FileParams), INTENT(IN) :: filepars

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*filepars* DERIVED  File attributes
!*cvals*    CHAR     Returned data string array
!
!------------------------------------------------------------------------------
!


CALL conv_to_chars(cvals(1),filepars%defined)
cvals(2) = filepars%form 
cvals(3) = filepars%filename
cvals(4) = filepars%floattype
CALL conv_to_chars(cvals(5),filepars%info)
cvals(6) = filepars%pathname 


RETURN

END SUBROUTINE conv_to_chars_outfiles

!========================================================================

SUBROUTINE conv_to_chars_outvars(cvals,varatts)
!************************************************************************
!
! *conv_to_chars_outvars* Write variable attributes to an array of 
!                         character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_to_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(11) :: cvals
TYPE (VariableAtts), INTENT(IN) :: varatts

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*varatts*  DERIVED  Variable attributes
!*cvals*    CHAR     Returned data string array
!
!------------------------------------------------------------------------------
!


CALL conv_to_chars(cvals(1),varatts%ivarid)
CALL conv_to_chars(cvals(2),varatts%nrank)
CALL conv_to_chars(cvals(3),varatts%numvar)
CALL conv_to_chars(cvals(4),varatts%oopt)
CALL conv_to_chars(cvals(5),varatts%klev)
CALL conv_to_chars(cvals(6),varatts%dep)
cvals(7) = varatts%node
cvals(8) = varatts%f90_name
cvals(9) = varatts%long_name
cvals(10) = varatts%units
cvals(11) = varatts%vector_name


RETURN

END SUBROUTINE conv_to_chars_outvars

!========================================================================

SUBROUTINE conv_to_chars_statlocs(cvals,statlocs)
!************************************************************************
!
! *conv_to_chars_statlocs* Write output station attributes to a array of 
!                          character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description -
!
! Module calls - conv_to_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(3) :: cvals
TYPE (StationLocs), INTENT(IN) :: statlocs

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*statlocs* DERIVED  Output station attributes
!*cvals*    CHAR     Returned data string array
!
!------------------------------------------------------------------------------
!


CALL conv_to_chars(cvals(1),statlocs%ipos)
CALL conv_to_chars(cvals(2),statlocs%jpos)
cvals(3) = statlocs%name


RETURN

END SUBROUTINE conv_to_chars_statlocs

!========================================================================

SUBROUTINE conv_to_chars_surfgrd(cvals,surfgrd)
!************************************************************************
!
! *conv_to_chars_surfgrd* Write surface grid attributes to an array of 
!                         character strings
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11.2
!
! Description -
!
! Module calls - conv_to_chars
!
!************************************************************************
!
USE datatypes

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(OUT), DIMENSION(14) :: cvals
TYPE (GridParams), INTENT(IN) :: surfgrd

!
! Name      Type     Purpose
!------------------------------------------------------------------------------
!*surfgrd*  DERIVED  Surface grid attributes
!*cvals*    CHAR     Returned data string array
!
!------------------------------------------------------------------------------
!


CALL conv_to_chars(cvals(1),surfgrd%nhtype)
CALL conv_to_chars(cvals(2),surfgrd%n1dat)
CALL conv_to_chars(cvals(3),surfgrd%n2dat)
CALL conv_to_chars(cvals(4),surfgrd%surcoords)
CALL conv_to_chars(cvals(5),surfgrd%delxdat)
CALL conv_to_chars(cvals(6),surfgrd%delydat)
CALL conv_to_chars(cvals(7),surfgrd%x0dat)
CALL conv_to_chars(cvals(8),surfgrd%y0dat)
CALL conv_to_chars(cvals(9),surfgrd%rotated)
CALL conv_to_chars(cvals(10),surfgrd%gridangle)
CALL conv_to_chars(cvals(11),surfgrd%y0rot)
CALL conv_to_chars(cvals(12),surfgrd%extrapol)
CALL conv_to_chars(cvals(13),surfgrd%land)
CALL conv_to_chars(cvals(14),surfgrd%datflag)

RETURN

END SUBROUTINE conv_to_chars_surfgrd

!========================================================================

SUBROUTINE error_cif_var(lvar)
!************************************************************************
!
! *error_cif_var* write an error message when a variable in the CIF cannot
!                 be read
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
! Module calls - error_abort
!
!************************************************************************
!
USE iopars
USE paralpars  
USE error_routines, ONLY: nerrs, error_abort

!
!*Arguments
!
INTEGER, INTENT(IN) :: lvar

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*lvar*      INTEGER Position of the variable on the CIF input line
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cnum


nerrs = nerrs + 1
IF (master.AND.errchk.AND.nerrs.LE.maxerrors) THEN
   WRITE (cnum,'(I12)') lvar; cnum = ADJUSTL(cnum)
   WRITE (ioerr,'(A)') 'Unable to read parameter '//TRIM(cnum)
   CALL error_abort(procname(pglev),ierrno_inival)
ENDIF


RETURN

END SUBROUTINE error_cif_var

!========================================================================

SUBROUTINE read_cif_line(cline,cvals,numvars,cname,filename)
!************************************************************************
!
! *read_cif_line* Parse a series of variables, separated by a delimiter as
!                 character data on a data (CIF) input line
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
! Module calls - error_abort, upper_case
!
!************************************************************************
!
USE iopars
USE paralpars
USE syspars
USE error_routines, ONLY: nerrs, error_abort
USE utility_routines, ONLY: upper_case

!
!*Arguments
!
CHARACTER (LEN=lencifline), INTENT(INOUT) :: cline
CHARACTER (LEN=lenname), INTENT(INOUT), OPTIONAL :: cname
CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: filename
CHARACTER (LEN=lencifvar), INTENT(OUT), DIMENSION(:) :: cvals
INTEGER, INTENT(OUT) :: numvars

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cline*     CHAR    Input line string
!*cvals*     CHAR    Data values in string format
!*numvars*   INTEGER Number of data variables
!*cname*     CHAR    Name of the variable(s)
!*filename*  CHAR    Name of the file if PRESENT. Nmae of the CIF file
!                    otherwise. 
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cnum
CHARACTER (LEN=lencifline) :: xline
INTEGER :: iend, ieq, isep, ival, l


!
!1. Initialise
!-------------
!

cvals = ''; numvars = 0

!
!2. Removing comments from input line
!------------------------------------
!

cline = ADJUSTL(cline)
iend = INDEX(cline,cifcom)
IF (iend.EQ.1) THEN
   cline = ''; cvals = ''; numvars = 0
   IF (PRESENT(cname)) cname = ''
   RETURN
ENDIF
IF (iend.GT.0) cline = ADJUSTL(cline(1:iend-1))

!
!3. Data variable name
!---------------------
!

IF (PRESENT(cname)) THEN
   ieq = INDEX(cline,'=')
   IF (ieq.GT.1) THEN
      cname = TRIM(ADJUSTL(cline(1:ieq-1)))
      CALL upper_case(cname)
   ELSE
      GOTO 1001
   ENDIF
ELSE
   ieq = 0
ENDIF

!
!4. Separate data
!----------------
!
!---first value
ival = 1
isep = INDEX(cline,cifsep)
IF (isep.GT.0.AND.isep.LT.ieq) THEN
   GOTO 1002
ELSEIF (isep.EQ.0) THEN
   IF (ieq.EQ.LEN(cline)) THEN
      cvals(1) = ''
   ELSE
      cvals(1) = cline(ieq+1:LEN(cline))
   ENDIF
ELSE
   IF (isep.EQ.ieq+1) THEN
      cvals(1) = ''
   ELSE
      cvals(1) = cline(ieq+1:isep-1)
   ENDIF
   IF (isep.EQ.LEN(cline)) THEN
      xline = ''
   ELSE
      xline = cline(isep+1:LEN(cline))
   ENDIF
ENDIF

!---next values
DO WHILE (isep.GT.0)
   ival = ival + 1
   l = INDEX(xline,cifsep)
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

!---exit code in case of error
1001 nerrs = nerrs + 1
IF (master.AND.errchk.AND.nerrs.LE.maxerrors) THEN

   WRITE (cnum,'(I12)') ciflinenum; cnum = ADJUSTL(cnum)
   IF (PRESENT(filename)) THEN
      WRITE (ioerr,'(A)') 'Error occurred on line '//TRIM(cnum)//' of file '&
                        & //TRIM(filename)
   ELSE
      WRITE (ioerr,'(A)') 'Error occurred on line '//TRIM(cnum)//' of file '&
                        & //TRIM(ciffile%filename)
   ENDIF
   WRITE (ioerr,'(A)') 'Variable name is not defined'
ENDIF
CALL error_abort(procname(pglev),ierrno_inival)
1002 nerrs = nerrs + 1
IF (master.AND.errchk.AND.nerrs.LE.maxerrors) THEN
   WRITE (cnum,'(I12)') ciflinenum; cnum = ADJUSTL(cnum)
   IF (PRESENT(filename)) THEN
      WRITE (ioerr,'(A)') 'Error occurred on line '//TRIM(cnum)//' of file '&
                        & //TRIM(filename)
      WRITE (ioerr,'(A)') 'Error occurred on line '//TRIM(cnum)//' of file '&
                        & //TRIM(ciffile%filename)
   ENDIF
   WRITE (ioerr,'(A)') 'Data separator appears before equal sign'
ENDIF
CALL error_abort(procname(pglev),ierrno_inival)

END SUBROUTINE read_cif_line

!========================================================================

SUBROUTINE write_cif_line(cvals,cname)
!************************************************************************
!
! *write_cif_line* write setup parameter(s) to a CIF file 
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)inout_routines.f90  V2.11
!
! Description - 
!
! Module calls -
!
!************************************************************************
!
USE iopars
USE syspars
!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: cname
CHARACTER (LEN=*), INTENT(INOUT), DIMENSION(:) :: cvals

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*cvals*     CHAR    Data values in string format
!*cname*     CHAR    Name of the variable(s)
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
   cline = TRIM(cvals(1))
ENDIF

WRITE (ciffile%iunit,'(A)') TRIM(cname)//' = '//TRIM(cline)

ciflinenum = ciflinenum + 1


RETURN

END SUBROUTINE write_cif_line


END MODULE cif_routines
