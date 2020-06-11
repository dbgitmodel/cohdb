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
! *Sediment_Parameters* Routines for reading and writing a particle CIF
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Parameters.f90  V2.11
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description -
!
! Routines - assign_cif_vars_part, assign_cif_vars_tspart, write_cif_vars_part,
!            write_cif_vars_tspart
!
!************************************************************************
!

!============================================================================

SUBROUTINE assign_cif_vars_part(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_part* convert the string data from the particle model CIF
!                        input line to the appropriate numeric or non-numeric
!                        format 
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Parameters.f90  V2.11
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
USE partpars
USE partswitches
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
INTEGER :: iddesc, lvar


iddesc = icif_mod
lvar = 1

SELECT CASE (TRIM(cname))

!
!1. Switches
!-----------
!

CASE('IOPT_PART_CLOUD')
   CALL conv_from_chars(cvals(1),iopt_part_cloud,lvar)
CASE('IOPT_PART_CONC')
   CALL conv_from_chars(cvals(1),iopt_part_conc,lvar)
CASE('IOPT_PART_DENS')
   CALL conv_from_chars(cvals(1),iopt_part_dens,lvar)
CASE('IOPT_PART_HADV')
   CALL conv_from_chars(cvals(1),iopt_part_hadv,lvar)
CASE('IOPT_PART_HDIF')
   CALL conv_from_chars(cvals(1),iopt_part_hdif,lvar)
CASE('IOPT_PART_LEEWAY')
   CALL conv_from_chars(cvals(1),iopt_part_leeway,lvar)
CASE('IOPT_PART_OUT')
   CALL conv_from_chars(cvals(1),iopt_part_out,lvar)
CASE('IOPT_PART_VADV')
   CALL conv_from_chars(cvals(1),iopt_part_vadv,lvar)
CASE('IOPT_PART_VDIF')
   CALL conv_from_chars(cvals(1),iopt_part_vdif,lvar)
CASE('IOPT_PART_WIND')
   CALL conv_from_chars(cvals(1),iopt_part_wind,lvar)   

!
!2. Parameters
!-------------
!

CASE('REFDATETIME_PART')
   CALL conv_from_chars(cvals(1),RefDateTime_part,lvar)
CASE('CDRIFT_FAC')
   CALL conv_from_chars(cvals(1),cdrift_fac,lvar)
CASE('ICPART'); CALL conv_from_chars(cvals(1),icpart,lvar)
CASE('NOCLOUDS'); CALL conv_from_chars(cvals(1),noclouds,lvar)
CASE('NOLABELS'); CALL conv_from_chars(cvals(1),nolabels,lvar)
CASE('NOPART'); CALL conv_from_chars(cvals(1),nopart,lvar)
CASE('NOSETSPART'); CALL conv_from_chars(cvals(1),nosetspart,lvar)
CASE('PTIME_UNIT')
   CALL conv_from_chars(cvals(1),ptime_unit,lvar)   
CASE('WDRIFT_ANGLE1')
   CALL conv_from_chars(cvals(1),wdrift_angle1,lvar)
CASE('WDRIFT_ANGLE2')
      CALL conv_from_chars(cvals(1),wdrift_angle2,lvar)
CASE('WDRIFT_SLOPE')
   CALL conv_from_chars(cvals(1),wdrift_slope,lvar)
CASE('XDIF_PART_CST')
   CALL conv_from_chars(cvals(1),xdifpart_cst,lvar)
CASE('YDIF_PART_CST')
   CALL conv_from_chars(cvals(1),ydifpart_cst,lvar)
CASE('ZDIF_PART_CST')
   CALL conv_from_chars(cvals(1),zdifpart_cst,lvar)
   
END SELECT


RETURN

END SUBROUTINE assign_cif_vars_part

!========================================================================

SUBROUTINE assign_cif_vars_tspart(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_tspart* convert string data from an input line in the CIF
!                          block with parameters for particle trajectory output
!                          to the appropriate numeric or non-numeric format
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Parameters.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - read_cif_params
!
! External calls -
!
! Module calls - check_cif_lbound_novars, conv_from_chars,
!                conv_from_chars_outfiles, error_array_index
!
!************************************************************************
!
USE partpars
USE partvars
USE syspars
USE cif_routines, ONLY: check_cif_lbound_novars, conv_from_chars, &
                      & conv_from_chars_outfiles
USE error_routines, ONLY: error_array_index

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lenname), INTENT(IN) :: cname
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
INTEGER :: iset, lvar


lvar = 1

SELECT CASE (TRIM(cname))

!---file attributes
CASE('PART1D')
   CALL check_cif_lbound_novars(numvars,7)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'part1d',1,nosetspart,1)
   CALL conv_from_chars_outfiles(cvals(2:7),part1d(iset),lvar)

!---output attributes
CASE('OUTPPARS')
   CALL check_cif_lbound_novars(numvars,10)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'outppars',1,nosetspart,1)
   lvar = 2
   CALL conv_from_chars(cvals(lvar),outppars(iset)%aging,lvar)
   lvar = 3
   CALL conv_from_chars(cvals(lvar),outppars(iset)%refdate,lvar)
   lvar = 4
   CALL conv_from_chars(cvals(lvar),outppars(iset)%ktype,lvar)
   lvar = 5
   CALL conv_from_chars(cvals(lvar),outppars(iset)%label,lvar)
   lvar = 6
   CALL conv_from_chars(cvals(lvar),outppars(iset)%nodim,lvar)
   lvar = 7
   CALL conv_from_chars(cvals(lvar),outppars(iset)%ntype,lvar)
   lvar = 8
   CALL conv_from_chars(cvals(lvar),outppars(iset)%tlims(1),lvar)
   lvar = 9
   CALL conv_from_chars(cvals(lvar),outppars(iset)%tlims(2),lvar)
   lvar = 10
   CALL conv_from_chars(cvals(lvar),outppars(iset)%tlims(3),lvar)

END SELECT


RETURN

END SUBROUTINE assign_cif_vars_tspart

!========================================================================

SUBROUTINE write_cif_vars_part
!************************************************************************
!
! *write_cif_vars_part* write particle model parameters to a CIF
!
! Author - Valerie Duliere
!
! Version - @(COHERENS)Particle_Parameters.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model, initialise_particle_model
!
! External calls -
!
! Module calls - conv_to_chars, write_cif_line
!
!************************************************************************
!
USE iopars
USE paralpars
USE partpars
USE partswitches
USE syspars
USE cif_routines, ONLY: conv_to_chars, write_cif_line
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: iunit


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_cif_vars_part' 
CALL log_timer_in()

!
!1. Write block header
!---------------------
!

iunit = ciffile%iunit
WRITE (iunit,'(A)') cifend//TRIM(cif_block_names(icif_part))

!
!2. Switches
!-----------
!

CALL conv_to_chars(cvals(1),iopt_part_cloud)
CALL write_cif_line(cvals(1:1),'IOPT_PART_CLOUD')
CALL conv_to_chars(cvals(1),iopt_part_conc)
CALL write_cif_line(cvals(1:1),'IOPT_PART_CONC')
CALL conv_to_chars(cvals(1),iopt_part_dens)
CALL write_cif_line(cvals(1:1),'IOPT_PART_DENS')
CALL conv_to_chars(cvals(1),iopt_part_hadv)
CALL write_cif_line(cvals(1:1),'IOPT_PART_HADV')
CALL conv_to_chars(cvals(1),iopt_part_hdif)
CALL write_cif_line(cvals(1:1),'IOPT_PART_HDIF')
CALL conv_to_chars(cvals(1),iopt_part_leeway)
CALL write_cif_line(cvals(1:1),'IOPT_PART_LEEWAY')
CALL conv_to_chars(cvals(1),iopt_part_out)
CALL write_cif_line(cvals(1:1),'IOPT_PART_OUT')
CALL conv_to_chars(cvals(1),iopt_part_vadv)
CALL write_cif_line(cvals(1:1),'IOPT_PART_VADV')
CALL conv_to_chars(cvals(1),iopt_part_vdif)
CALL write_cif_line(cvals(1:1),'IOPT_PART_VDIF')
CALL conv_to_chars(cvals(1),iopt_part_wind)
CALL write_cif_line(cvals(1:1),'IOPT_PART_WIND')

!
!3. Parameters
!-------------
!

CALL conv_to_chars(cvals(1),RefDateTime_part)
CALL write_cif_line(cvals(1:1),'REFDATETIME_PART')
CALL conv_to_chars(cvals(1),cdrift_fac)
CALL write_cif_line(cvals(1:1),'CDRIFT_FAC')
CALL conv_to_chars(cvals(1),icpart)
CALL write_cif_line(cvals(1:1),'ICPART')
CALL conv_to_chars(cvals(1),noclouds)
CALL write_cif_line(cvals(1:1),'NOCLOUDS')
CALL conv_to_chars(cvals(1),nolabels)
CALL write_cif_line(cvals(1:1),'NOLABELS')
CALL conv_to_chars(cvals(1),nopart)
CALL write_cif_line(cvals(1:1),'NOPART')
CALL conv_to_chars(cvals(1),nosetspart)
CALL write_cif_line(cvals(1:1),'NOSETSPART')
CALL conv_to_chars(cvals(1),ptime_unit)
CALL write_cif_line(cvals(1:1),'PTIME_UNIT')
CALL conv_to_chars(cvals(1),wdrift_angle1)
CALL write_cif_line(cvals(1:1),'WDRIFT_ANGLE1')
CALL conv_to_chars(cvals(1),wdrift_angle2)
CALL write_cif_line(cvals(1:1),'WDRIFT_ANGLE2')
CALL conv_to_chars(cvals(1),wdrift_slope)
CALL write_cif_line(cvals(1:1),'WDRIFT_SLOPE')
CALL conv_to_chars(cvals(1),xdifpart_cst)
CALL write_cif_line(cvals(1:1),'XDIF_PART_CST')
CALL conv_to_chars(cvals(1),ydifpart_cst)
CALL write_cif_line(cvals(1:1),'YDIF_PART_CST')
CALL conv_to_chars(cvals(1),zdifpart_cst)
CALL write_cif_line(cvals(1:1),'ZDIF_PART_CST')

!
!4. Write end of block
!----------------------
!

WRITE (ciffile%iunit,'(A)') cifend
ciflinenum = ciflinenum + 1

CALL log_timer_out()


RETURN

END SUBROUTINE write_cif_vars_part

!========================================================================

SUBROUTINE write_cif_vars_tspart
!************************************************************************
!
! *write_cif_vars_tspart* Write the CIF block with parameters for particle
!                         trajectory output
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Parameters.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - conv_to_chars, conv_to_chars_outfiles, write_cif_line
!
!************************************************************************
!
USE iopars
USE paralpars
USE partpars
USE partvars
USE syspars
USE cif_routines, ONLY: conv_to_chars, conv_to_chars_outfiles, write_cif_line
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: iset, iunit


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_cif_vars_tspart'
CALL log_timer_in()

!---write block header
iunit = ciffile%iunit
WRITE (iunit,'(A)') cifend//TRIM(cif_block_names(icif_tspart))

!---file attributes
iset_110: DO iset=1,nosetspart
   IF (part1d(iset)%defined) THEN
      CALL conv_to_chars(cvals(1),iset)
      CALL conv_to_chars_outfiles(cvals(2:7),part1d(iset))
      CALL write_cif_line(cvals(1:7),'PART1D')
   ENDIF
ENDDO iset_110

!---output grid (space/time) attributes
iset_120: DO iset=1,nosetspart
   CALL conv_to_chars(cvals(1),iset)
   CALL conv_to_chars(cvals(2),outppars(iset)%aging)
   CALL conv_to_chars(cvals(3),outppars(iset)%refdate)
   CALL conv_to_chars(cvals(4),outppars(iset)%ktype)
   CALL conv_to_chars(cvals(5),outppars(iset)%label)
   CALL conv_to_chars(cvals(6),outppars(iset)%nodim)
   CALL conv_to_chars(cvals(7),outppars(iset)%ntype)
   CALL conv_to_chars(cvals(8),outppars(iset)%tlims(1))
   CALL conv_to_chars(cvals(9),outppars(iset)%tlims(2))
   CALL conv_to_chars(cvals(10),outppars(iset)%tlims(3))
   CALL write_cif_line(cvals(1:10),'OUTPPARS')
ENDDO iset_120

WRITE (iunit,'(A)') cifend
ciflinenum = ciflinenum + 1

CALL log_timer_out()

RETURN

END SUBROUTINE write_cif_vars_tspart
