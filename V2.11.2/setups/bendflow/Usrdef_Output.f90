! This file is part of COHERENS. You can redistribute and/or modify it under
! the conditions of the COHERENS license. Details are found in the file
! COHERENS_License.

SUBROUTINE usrdef_output
!************************************************************************
!
! *usrdef_output* User-formatted output
!
! Author - IMDC
!
! Version - @(COHERENS)Usrdef_Output.f90  V2.5
!
! $Date: 2017-08-21 13:49:46 +0200 (Mon, 21 Aug 2017) $
!
! $Revision: 1044 $
!
! Description - test case bendflow
!
! Reference -
!
! Calling program - coherens_main
!
!************************************************************************
!
USE iopars
USE timepars

IMPLICIT NONE


!
!1. Reset forcing attributes for CIF file
!----------------------------------------
!

IF (nt.EQ.0.AND.ciffile%status.EQ.'W') THEN
!  ---model grid
   modfiles(io_modgrd,1,1)%status = 'R'
   modfiles(io_modgrd,1,1)%form = 'U'
   modfiles(io_modgrd,1,1)%filename = modfiles(io_modgrd,1,2)%filename
   modfiles(io_modgrd,1,2)%status = '0'
!  ---initial conditions
   modfiles(io_inicon,ics_phys,1)%status = 'R'
   modfiles(io_inicon,ics_phys,1)%form = 'U'
   modfiles(io_inicon,ics_phys,1)%filename = modfiles(io_fincon,1,2)%filename
   modfiles(io_fincon,ics_phys,2)%status = '0'
!  ---open boundary conditions (2-D)
   modfiles(io_2uvobc,1:3,1)%status = 'R'
   modfiles(io_2uvobc,1:3,1)%form = 'U'
   modfiles(io_2uvobc,1:3,1)%filename = modfiles(io_2uvobc,1:3,2)%filename
   modfiles(io_2uvobc,1:3,2)%status = '0'
!  ---sediment particle attributes
   modfiles(io_sedspc,1,1)%status = 'R'
   modfiles(io_sedspc,1,1)%form = 'U'
   modfiles(io_sedspc,1,1)%filename = modfiles(io_sedspc,1,2)%filename
   modfiles(io_sedspc,1,2)%status = '0'
ENDIF


RETURN

END SUBROUTINE usrdef_output







