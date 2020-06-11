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
! *Model_Parameters* Define model parameters through CIF
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11.2
!
! $Date: 2018-07-13 10:07:41 +0200 (Fri, 13 Jul 2018) $
!
! $Revision: 1163 $
!
! Description -
!
! Routines - assign_cif_vars_anal, assign_cif_vars_avrgd, assign_cif_vars_harm,
!            assign_cif_vars_mod, assign_cif_vars_mon, assign_cif_vars_tsout,
!            read_cif_params, write_cif_vars_anal, write_cif_vars_avrgd,
!            write_cif_vars_harm, write_cif_vars_mod, write_cif_vars_mon,
!            write_cif_vars_tsout
!
!************************************************************************
!

!========================================================================

SUBROUTINE assign_cif_vars_anal(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_anal* convert string data from an input line in the CIF
!                        block with parameters for harmonic output to the
!                        appropriate numeric or non-numeric format   
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11
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
!                conv_from_chars_gridpars, conv_from_chars_outfiles,
!                conv_from_chars_statlocs, conv_from_chars_outvars,
!                error_array_index
!
!************************************************************************
!
USE iopars
USE syspars
USE cif_routines, ONLY: check_cif_lbound_novars, conv_from_chars, &
                      & conv_from_chars_gridpars, conv_from_chars_outfiles, &
                      & conv_from_chars_statlocs, conv_from_chars_outvars
USE error_routines, ONLY: error_array_index

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
INTEGER :: ifreq, iset, istat, ivar, lvar, nostats


lvar = 1

SELECT CASE (TRIM(cname))

!
!1. Variable attributes
!----------------------
!

CASE('ANALVARS')
   CALL check_cif_lbound_novars(numvars,12)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),ivar,lvar)
   CALL error_array_index(ivar,'analvars',1,novarsanal,1)
   CALL conv_from_chars_outvars(cvals(2:12),analvars(ivar),lvar)

!
!2.Variable indices
!--------------------
!

CASE('IVARSANAL')
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'ivarsanal',1,nosetsanal,1)
   lvar_210: DO lvar=2,numvars
      CALL conv_from_chars(cvals(lvar),ivarsanal(iset,lvar-1),lvar)
   ENDDO lvar_210

!
!3. Elliptic variables
!---------------------
!
!---variable numbers
CASE('IVARSELL')
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'ivarsell',1,nosetsanal,1)
   lvar_310: DO lvar=2,numvars
      CALL conv_from_chars(cvals(lvar),ivarsell(iset,lvar-1),lvar)
   ENDDO lvar_310

!---components of elliptic vector (2-D)
CASE('IVECELL2D')
   CALL check_cif_lbound_novars(numvars,3)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'ivecell2d',1,nosetsanal,1)
   lvar_320: DO lvar=2,3
      CALL conv_from_chars(cvals(lvar),ivecell2d(iset,lvar-1),lvar)
   ENDDO lvar_320

!---components of elliptic vector (3-D)
CASE('IVECELL3D')
   CALL check_cif_lbound_novars(numvars,3)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'ivecell3d',1,nosetsanal,1)
   lvar_330: DO lvar=2,3
      CALL conv_from_chars(cvals(lvar),ivecell3d(iset,lvar-1),lvar)
   ENDDO lvar_330

!
!4. File attributes
!------------------
!

CASE('RES0D')
   CALL check_cif_lbound_novars(numvars,7)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'res0d',1,nosetsanal,1)
   CALL conv_from_chars_outfiles(cvals(2:7),res0d(iset),lvar)

CASE('RES2D')
   CALL check_cif_lbound_novars(numvars,7)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'res2d',1,nosetsanal,1)
   CALL conv_from_chars_outfiles(cvals(2:7),res2d(iset),lvar)

CASE('RES3D')
   CALL check_cif_lbound_novars(numvars,7)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'res3d',1,nosetsanal,1)
   CALL conv_from_chars_outfiles(cvals(2:7),res3d(iset),lvar)

CASE('AMP2D')
   CALL check_cif_lbound_novars(numvars,8)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'amp2d',1,nosetsanal,1)
   lvar = 2
   CALL conv_from_chars(cvals(lvar),ifreq,lvar)
   CALL error_array_index(ifreq,'amp2d',1,nofreqsharm(iset),2)
   CALL conv_from_chars_outfiles(cvals(3:8),amp2d(iset,ifreq),lvar)

CASE('AMP3D')
   CALL check_cif_lbound_novars(numvars,8)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'amp3d',1,nosetsanal,1)
   lvar = 2
   CALL conv_from_chars(cvals(lvar),ifreq,lvar)
   CALL error_array_index(ifreq,'amp3d',1,nofreqsharm(iset),2)
   CALL conv_from_chars_outfiles(cvals(3:8),amp3d(iset,ifreq),lvar)

CASE('PHA2D')
   CALL check_cif_lbound_novars(numvars,8)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'pha2d',1,nosetsanal,1)
   lvar = 2
   CALL conv_from_chars(cvals(lvar),ifreq,lvar)
   CALL error_array_index(ifreq,'pha2d',1,nofreqsharm(iset),2)
   CALL conv_from_chars_outfiles(cvals(3:8),pha2d(iset,ifreq),lvar)

CASE('PHA3D')
   CALL check_cif_lbound_novars(numvars,8)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'pha3d',1,nosetsanal,1)
   lvar = 2
   CALL conv_from_chars(cvals(lvar),ifreq,lvar)
   CALL error_array_index(ifreq,'pha3d',1,nofreqsharm(iset),2)
   CALL conv_from_chars_outfiles(cvals(3:8),pha3d(iset,ifreq),lvar)

CASE('ELL2D')
   CALL check_cif_lbound_novars(numvars,8)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'ell2d',1,nosetsanal,1)
   lvar = 2
   CALL conv_from_chars(cvals(lvar),ifreq,lvar)
   CALL error_array_index(ifreq,'ell2d',1,nofreqsharm(iset),2)
   CALL conv_from_chars_outfiles(cvals(3:8),ell2d(iset,ifreq),lvar)

CASE('ELL3D')
   CALL check_cif_lbound_novars(numvars,8)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'ell3d',1,nosetsanal,1)
   lvar = 2
   CALL conv_from_chars(cvals(lvar),ifreq,lvar)
   CALL error_array_index(ifreq,'ell3d',1,nofreqsharm(iset),2)
   CALL conv_from_chars_outfiles(cvals(3:8),ell3d(iset,ifreq),lvar)

!
!5. Output grid (space/time) attributes
!--------------------------------------
!

CASE('ANALGPARS')
   CALL check_cif_lbound_novars(numvars,20)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'analgpars',1,nosetsanal,1)
   CALL conv_from_chars_gridpars(cvals(2:20),analgpars(iset),lvar)

!
!6. Station attributes
!---------------------
!

CASE('ANALSTATLOCS')
   CALL check_cif_lbound_novars(numvars,4)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),istat,lvar)
   CALL error_array_index(istat,'analstatlocs',1,nostatsanal,1)
   CALL conv_from_chars_statlocs(cvals(2:4),analstatlocs(istat),lvar)

!
!7. Station labels
!-----------------
!

CASE('LSTATSANAL')
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   nostats = analgpars(iset)%nostats
   CALL error_array_index(iset,'lstatsanal',1,nosetsanal,1)
   CALL check_cif_lbound_novars(numvars,nostats+1)
   lvar_710: DO lvar=2,nostats+1
      CALL conv_from_chars(cvals(lvar),lstatsanal(iset,lvar-1),lvar)
   ENDDO lvar_710

END SELECT


RETURN

END SUBROUTINE assign_cif_vars_anal

!========================================================================

SUBROUTINE assign_cif_vars_avrgd(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_avrgd* convert string data from an input line in the CIF
!                         block with parameters for time averaging to the
!                         appropriate numeric or non-numeric format   
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11
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
!                conv_from_chars_gridpars, conv_from_chars_outfiles,
!                conv_from_chars_statlocs, conv_from_chars_outvars,
!                error_array_index
!
!************************************************************************
!
USE iopars
USE syspars
USE cif_routines, ONLY: check_cif_lbound_novars, conv_from_chars, &
                      & conv_from_chars_gridpars, conv_from_chars_outfiles, &
                      & conv_from_chars_statlocs, conv_from_chars_outvars
USE error_routines, ONLY: error_array_index

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
INTEGER :: iset, istat, ivar, lvar, nostats


lvar = 1

SELECT CASE (TRIM(cname))

!---variable attributes
CASE('AVRVARS')
   CALL check_cif_lbound_novars(numvars,12)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),ivar,lvar)
   CALL error_array_index(ivar,'avrvars',1,novarsavr,1)
   CALL conv_from_chars_outvars(cvals(2:12),avrvars(ivar),lvar)

!---variable indices
CASE('IVARSAVR')
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'ivarstsr',1,nosetstsr,1)
   ivar_110: DO lvar=2,numvars
      CALL conv_from_chars(cvals(lvar),ivarsavr(iset,lvar-1),lvar)
   ENDDO ivar_110

!---file attributes
CASE('AVR0D')
   CALL check_cif_lbound_novars(numvars,7)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'avr0d',1,nosetsavr,1)
   CALL conv_from_chars_outfiles(cvals(2:7),avr0d(iset),lvar)

CASE('AVR2D')
   CALL check_cif_lbound_novars(numvars,7)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'avr2d',1,nosetsavr,1)
   CALL conv_from_chars_outfiles(cvals(2:7),avr2d(iset),lvar)

CASE('AVR3D')
   CALL check_cif_lbound_novars(numvars,7)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'avr3d',1,nosetsavr,1)
   CALL conv_from_chars_outfiles(cvals(2:7),avr3d(iset),lvar)

!---output grid (space/time) attributes
CASE('AVRGPARS')
   CALL check_cif_lbound_novars(numvars,20)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'avrgpars',1,nosetsavr,1)
   CALL conv_from_chars_gridpars(cvals(2:20),avrgpars(iset),lvar)

!---station attributes
CASE('AVRSTATLOCS')
   CALL check_cif_lbound_novars(numvars,4)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),istat,lvar)
   CALL error_array_index(istat,'avrstatlocs',1,nostatsavr,1)
   CALL conv_from_chars_statlocs(cvals(2:4),avrstatlocs(istat),lvar)

!--station labels
CASE('LSTATSAVR')
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   nostats = avrgpars(iset)%nostats
   CALL error_array_index(iset,'lstatsavr',1,nosetsavr,1)
   CALL check_cif_lbound_novars(numvars,nostats+1)
   lvar_120: DO lvar=2,nostats+1
      CALL conv_from_chars(cvals(lvar),lstatsavr(iset,lvar-1),lvar)
   ENDDO lvar_120

END SELECT


RETURN

END SUBROUTINE assign_cif_vars_avrgd

!========================================================================

SUBROUTINE assign_cif_vars_harm(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_harm* convert string data from an input line in the CIF
!                        block with parameters for harmonic analysis to the
!                        appropriate numeric or non-numeric format   
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - read_cif_params
!
! External calls -
!
! Module calls - check_cif_lbound_novars, conv_from_chars, error_array_index
!
!************************************************************************
!
USE iopars
USE syspars
USE cif_routines, ONLY: check_cif_lbound_novars, conv_from_chars
USE error_routines, ONLY: error_array_index

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

!---harmonic tidal index ID
CASE('INDEX_ANAL')
   CALL check_cif_lbound_novars(numvars,nofreqsanal)
   lvar_110: DO lvar=1,nofreqsanal
      CALL conv_from_chars(cvals(lvar),index_anal(lvar),lvar)
   ENDDO lvar_110

!---frequency names
CASE('HARM_FREQ_NAMES')
   CALL check_cif_lbound_novars(numvars,nofreqsanal)
   lvar_120: DO lvar=1,nofreqsanal
      CALL conv_from_chars(cvals(lvar),harm_freq_names(lvar),lvar)
   ENDDO lvar_120

!---frequency values
CASE('HARM_FREQ')
   CALL check_cif_lbound_novars(numvars,nofreqsanal)
   lvar_130: DO lvar=1,nofreqsanal
      CALL conv_from_chars(cvals(lvar),harm_freq(lvar),lvar)
   ENDDO lvar_130

!---reference times
CASE('CDATE_TIME_REF')
   CALL check_cif_lbound_novars(numvars,nosetsanal)
   lvar_140: DO lvar=1,nosetsanal
      CALL conv_from_chars(cvals(lvar),cdate_time_ref(lvar),lvar)
   ENDDO lvar_140

!---number of requencies per set
CASE('NOFREQSHARM')
   CALL check_cif_lbound_novars(numvars,nosetsanal)
   lvar_150: DO lvar=1,nosetsanal
      CALL conv_from_chars(cvals(lvar),nofreqsharm(lvar),lvar)
   ENDDO lvar_150

!---frequency indices
CASE('IFREQSHARM')
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'ifreqsharm',1,nosetsanal,1)
   CALL check_cif_lbound_novars(numvars,nofreqsharm(iset)+1)
   lvar_160: DO lvar=2,nofreqsharm(iset)+1
      CALL conv_from_chars(cvals(lvar),ifreqsharm(iset,lvar-1),lvar)
   ENDDO lvar_160

!---time step for harmonic analysis
CASE('ICANAL')
   CALL check_cif_lbound_novars(numvars,nosetsanal)
   lvar_170: DO lvar=1,nosetsanal
      CALL conv_from_chars(cvals(lvar),icanal(lvar),lvar)
   ENDDO lvar_170

END SELECT


RETURN

END SUBROUTINE assign_cif_vars_harm

!========================================================================

SUBROUTINE assign_cif_vars_mod(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_mod* convert string data from an input line in the CIF
!                       block with model setup parameters to the
!                       appropriate numeric or non-numeric format   
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11.2
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
!                conv_from_chars_infiles, conv_from_chars_surfgrd,
!                error_array_index
!
!************************************************************************
!
USE gridpars
USE iopars
USE nestgrids
USE obconds
USE paralpars
USE physpars
USE relaxation
USE structures
USE switches
USE syspars
USE tide
USE timepars
USE turbpars
USE cif_routines, ONLY: check_cif_lbound_novars, conv_from_chars, &
                      & conv_from_chars_infiles, conv_from_chars_surfgrd
USE error_routines, ONLY: error_array_index

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
CHARACTER (LEN=6) :: fdesc
INTEGER :: idesc, ifil, iotype, lvar, n


lvar = 1

SELECT CASE (TRIM(cname))

!
!1. Switches
!------------
!

CASE('IOPT_ADV_SCAL'); CALL conv_from_chars(cvals(1),iopt_adv_scal,lvar)
CASE('IOPT_ADV_TURB'); CALL conv_from_chars(cvals(1),iopt_adv_turb,lvar)
CASE('IOPT_ADV_TVD'); CALL conv_from_chars(cvals(1),iopt_adv_tvd,lvar)
CASE('IOPT_ADV_2D'); CALL conv_from_chars(cvals(1),iopt_adv_2D,lvar)
CASE('IOPT_ADV_3D'); CALL conv_from_chars(cvals(1),iopt_adv_3D,lvar)
CASE('IOPT_ARRINT_DEPTHS');
   CALL conv_from_chars(cvals(1),iopt_arrint_hreg,lvar)
CASE('IOPT_ARRINT_HREG');
   CALL conv_from_chars(cvals(1),iopt_arrint_hreg,lvar)
CASE('IOPT_ARRINT_VREG')
   CALL conv_from_chars(cvals(1),iopt_arrint_vreg,lvar)
CASE('IOPT_ARRINT_3D')
   CALL conv_from_chars(cvals(1),iopt_arrint_3D,lvar)
CASE('IOPT_ASTRO_ANAL')
   CALL conv_from_chars(cvals(1),iopt_astro_anal,lvar)
CASE('IOPT_ASTRO_PARS')
   CALL conv_from_chars(cvals(1),iopt_astro_pars,lvar)
CASE('IOPT_ASTRO_TIDE')
   CALL conv_from_chars(cvals(1),iopt_astro_tide,lvar)
CASE('IOPT_BIOLGY'); CALL conv_from_chars(cvals(1),iopt_biolgy,lvar)
CASE('IOPT_BSTRES_DRAG')
   CALL conv_from_chars(cvals(1),iopt_bstres_drag,lvar)
CASE('IOPT_BSTRES_FORM')
   CALL conv_from_chars(cvals(1),iopt_bstres_form,lvar)
CASE('IOPT_BSTRES_NODIM')
   CALL conv_from_chars(cvals(1),iopt_bstres_nodim,lvar)
CASE('IOPT_BSTRES_WAVES_BFRIC')
   CALL conv_from_chars(cvals(1),iopt_bstres_waves_bfric,lvar)
CASE('IOPT_CDF_ABORT')
   CALL conv_from_chars(cvals(1),iopt_CDF_abort,lvar)
CASE('IOPT_CDF_FILL'); CALL conv_from_chars(cvals(1),iopt_CDF_fill,lvar)
CASE('IOPT_CDF_FORMAT')
   CALL conv_from_chars(cvals(1),iopt_CDF_format,lvar)
CASE('IOPT_CDF_SHARED')
   CALL conv_from_chars(cvals(1),iopt_CDF_shared,lvar)
CASE('IOPT_CDF_SYNC'); CALL conv_from_chars(cvals(1),iopt_CDF_sync,lvar)
CASE('IOPT_CDF_TLIM'); CALL conv_from_chars(cvals(1),iopt_CDF_tlim,lvar)
CASE('IOPT_COR_IMPL')
   CALL conv_from_chars(cvals(1),iopt_cor_impl,lvar)
CASE('IOPT_CURR'); CALL conv_from_chars(cvals(1),iopt_curr,lvar)
CASE('IOPT_CURR_WFALL')
   CALL conv_from_chars(cvals(1),iopt_curr_wfall,lvar)
CASE('IOPT_DAR'); CALL conv_from_chars(cvals(1),iopt_dar,lvar)
CASE('IOPT_DENS'); CALL conv_from_chars(cvals(1),iopt_dens,lvar)
CASE('IOPT_DENS_CONVECT')
   CALL conv_from_chars(cvals(1),iopt_dens_convect,lvar)
CASE('IOPT_DENS_GRAD')
   CALL conv_from_chars(cvals(1),iopt_dens_grad,lvar)
CASE('IOPT_DISCHR'); CALL conv_from_chars(cvals(1),iopt_dischr,lvar)
CASE('IOPT_DISCHR_LAND')
   CALL conv_from_chars(cvals(1),iopt_dischr_land,lvar)
CASE('IOPT_DRYCEL'); CALL conv_from_chars(cvals(1),iopt_drycel,lvar)
CASE('IOPT_FLD'); CALL conv_from_chars(cvals(1),iopt_fld,lvar)
CASE('IOPT_FLD_ALPHA'); CALL conv_from_chars(cvals(1),iopt_fld_alpha,lvar)
CASE('IOPT_GRID_HTYPE')
   CALL conv_from_chars(cvals(1),iopt_grid_htype,lvar)
CASE('IOPT_GRID_NODIM')
   CALL conv_from_chars(cvals(1),iopt_grid_nodim,lvar)
CASE('IOPT_GRID_SPH')
   CALL conv_from_chars(cvals(1),iopt_grid_sph,lvar)
CASE('IOPT_GRID_VTYPE')
   CALL conv_from_chars(cvals(1),iopt_grid_vtype,lvar)
CASE('IOPT_GRID_VTYPE_TRANSF')
   CALL conv_from_chars(cvals(1),iopt_grid_vtype_transf,lvar)
CASE('IOPT_HDIF_COEF')
   CALL conv_from_chars(cvals(1),iopt_hdif_coef,lvar)
CASE('IOPT_HDIF_LIM')
   CALL conv_from_chars(cvals(1),iopt_hdif_lim,lvar)
CASE('IOPT_HDIF_SCAL')
   CALL conv_from_chars(cvals(1),iopt_hdif_scal,lvar)
CASE('IOPT_HDIF_TURB')
   CALL conv_from_chars(cvals(1),iopt_hdif_turb,lvar)
CASE('IOPT_HDIF_2D'); CALL conv_from_chars(cvals(1),iopt_hdif_2D,lvar)
CASE('IOPT_HDIF_3D'); CALL conv_from_chars(cvals(1),iopt_hdif_3D,lvar)
CASE('IOPT_HYDRO_IMPL')
   CALL conv_from_chars(cvals(1),iopt_hydro_impl,lvar)
CASE('IOPT_KINVISC'); CALL conv_from_chars(cvals(1),iopt_kinvisc,lvar)
CASE('IOPT_METEO'); CALL conv_from_chars(cvals(1),iopt_meteo,lvar)
CASE('IOPT_METEO_DATA')
   CALL conv_from_chars(cvals(1),iopt_meteo_data,lvar)
CASE('IOPT_METEO_HEAT')
   CALL conv_from_chars(cvals(1),iopt_meteo_heat,lvar)
   CASE('IOPT_METEO_PRECIP')
   CALL conv_from_chars(cvals(1),iopt_meteo_precip,lvar)
CASE('IOPT_METEO_PRES')
   CALL conv_from_chars(cvals(1),iopt_meteo_pres,lvar)
CASE('IOPT_METEO_STRES')
   CALL conv_from_chars(cvals(1),iopt_meteo_stres,lvar)
CASE('IOPT_MG_CYCLE')
   CALL conv_from_chars(cvals(1),iopt_mg_cycle,lvar)
CASE('IOPT_MG_PROLONG')
   CALL conv_from_chars(cvals(1),iopt_mg_prolong,lvar)
CASE('IOPT_MG_SMOOTHER')
   CALL conv_from_chars(cvals(1),iopt_mg_smoother,lvar)
CASE('IOPT_MORPH') 
   CALL conv_from_chars(cvals(1),iopt_morph,lvar)
CASE('IOPT_MPI_ABORT')
   CALL conv_from_chars(cvals(1),iopt_MPI_abort,lvar)
CASE('IOPT_MPI_COMM_ALL')
   CALL conv_from_chars(cvals(1),iopt_MPI_comm_all,lvar)
CASE('IOPT_MPI_COMM_COLL')
   CALL conv_from_chars(cvals(1),iopt_MPI_comm_coll,lvar)
CASE('IOPT_MPI_COMM_EXCH')
   CALL conv_from_chars(cvals(1),iopt_MPI_comm_exch,lvar)
CASE('IOPT_MPI_COMM_FULL')
   CALL conv_from_chars(cvals(1),iopt_MPI_comm_full,lvar)
CASE('IOPT_MPI_COMM_GATH')
   CALL conv_from_chars(cvals(1),iopt_MPI_comm_gath,lvar)
CASE('IOPT_MPI_COMM_SCAT');
   CALL conv_from_chars(cvals(1),iopt_MPI_comm_scat,lvar)
CASE('IOPT_MPI_PARTIT')
   CALL conv_from_chars(cvals(1),iopt_MPI_partit,lvar)
CASE('IOPT_MPI_SYNC'); CALL conv_from_chars(cvals(1),iopt_MPI_sync,lvar)
CASE('IOPT_NESTS'); CALL conv_from_chars(cvals(1),iopt_nests,lvar)
CASE('IOPT_OBC_ADVFLUX')
   CALL conv_from_chars(cvals(1),iopt_obc_advflux,lvar)
CASE('IOPT_OBC_ADVRLX')
   CALL conv_from_chars(cvals(1),iopt_obc_advrlx,lvar)
CASE('IOPT_OBC_BIO'); CALL conv_from_chars(cvals(1),iopt_obc_bio,lvar)
CASE('IOPT_OBC_INVBAR')
   CALL conv_from_chars(cvals(1),iopt_obc_invbar,lvar)
CASE('IOPT_OBC_RELAX')
   CALL conv_from_chars(cvals(1),iopt_obc_relax,lvar)
CASE('IOPT_OBC_SAL'); CALL conv_from_chars(cvals(1),iopt_obc_sal,lvar)
CASE('IOPT_OBC_SED'); CALL conv_from_chars(cvals(1),iopt_obc_sed,lvar)
CASE('IOPT_OBC_TEMP'); CALL conv_from_chars(cvals(1),iopt_obc_temp,lvar)
CASE('IOPT_OBC_TH'); CALL conv_from_chars(cvals(1),iopt_obc_th,lvar)
CASE('IOPT_OBC_2D'); CALL conv_from_chars(cvals(1),iopt_obc_2D,lvar)
CASE('IOPT_OBC_2D_TANG')
   CALL conv_from_chars(cvals(1),iopt_obc_2D_tang,lvar)
CASE('IOPT_OBC_3D'); CALL conv_from_chars(cvals(1),iopt_obc_3D,lvar)
CASE('IOPT_OBC_3D_TANG')
   CALL conv_from_chars(cvals(1),iopt_obc_3D_tang,lvar)
CASE('IOPT_OUT_ANAL'); CALL conv_from_chars(cvals(1),iopt_out_anal,lvar)
CASE('IOPT_OUT_AVRGD')
   CALL conv_from_chars(cvals(1),iopt_out_avrgd,lvar)
CASE('IOPT_OUT_TSERS')
   CALL conv_from_chars(cvals(1),iopt_out_tsers,lvar)
CASE('IOPT_PART_WRITE'); CALL conv_from_chars(cvals(1),iopt_part_write,lvar)
CASE('IOPT_SAL'); CALL conv_from_chars(cvals(1),iopt_sal,lvar)
CASE('IOPT_RNG_SEED'); CALL conv_from_chars(cvals(1),iopt_rng_seed,lvar)
CASE('IOPT_SCAL_DEPOS')
   CALL conv_from_chars(cvals(1),iopt_scal_depos,lvar)
CASE('IOPT_SED'); CALL conv_from_chars(cvals(1),iopt_sed,lvar)
CASE('IOPT_SFLUX_PARS')
   CALL conv_from_chars(cvals(1),iopt_sflux_pars,lvar)
CASE('IOPT_SFLUX_PRECIP')
   CALL conv_from_chars(cvals(1),iopt_sflux_precip,lvar)
CASE('IOPT_SFLUX_QLONG')
   CALL conv_from_chars(cvals(1),iopt_sflux_qlong,lvar)
CASE('IOPT_SFLUX_QSHORT')
   CALL conv_from_chars(cvals(1),iopt_sflux_qshort,lvar)
CASE('IOPT_SFLUX_STRAT')
   CALL conv_from_chars(cvals(1),iopt_sflux_strat,lvar)
CASE('IOPT_SUR_1D'); CALL conv_from_chars(cvals(1),iopt_sur_1D,lvar)
CASE('IOPT_TEMP'); CALL conv_from_chars(cvals(1),iopt_temp,lvar)
CASE('IOPT_TEMP_OPTIC')
   CALL conv_from_chars(cvals(1),iopt_temp_optic,lvar)
CASE('IOPT_TEMP_SBC'); CALL conv_from_chars(cvals(1),iopt_temp_sbc,lvar)
CASE('IOPT_THNDAM'); CALL conv_from_chars(cvals(1),iopt_thndam,lvar)
CASE('IOPT_TIDAL_ACCEL') 
   CALL conv_from_chars(cvals(1),iopt_tidal_accel,lvar)
CASE('IOPT_TRANSP_FULL')
   CALL conv_from_chars(cvals(1),iopt_transp_full,lvar)
CASE('IOPT_TURB_ALG'); CALL conv_from_chars(cvals(1),iopt_turb_alg,lvar)
CASE('IOPT_TURB_DIS_BBC')
   CALL conv_from_chars(cvals(1),iopt_turb_dis_bbc,lvar)
CASE('IOPT_TURB_DIS_SBC')
   CALL conv_from_chars(cvals(1),iopt_turb_dis_sbc,lvar)
CASE('IOPT_TURB_IWLIM')
   CALL conv_from_chars(cvals(1),iopt_turb_iwlim,lvar)
CASE('IOPT_TURB_KINVISC')
   CALL conv_from_chars(cvals(1),iopt_turb_kinvisc,lvar)
CASE('IOPT_TURB_LMIX')
   CALL conv_from_chars(cvals(1),iopt_turb_lmix,lvar)
CASE('IOPT_TURB_NTRANS')
   CALL conv_from_chars(cvals(1),iopt_turb_ntrans,lvar)
CASE('IOPT_TURB_PARAM')
   CALL conv_from_chars(cvals(1),iopt_turb_param,lvar)
CASE('IOPT_TURB_STAB_FORM');
   CALL conv_from_chars(cvals(1),iopt_turb_stab_form,lvar)
CASE('IOPT_TURB_STAB_LEV')
   CALL conv_from_chars(cvals(1),iopt_turb_stab_lev,lvar)
CASE('IOPT_TURB_STAB_MOD')
   CALL conv_from_chars(cvals(1),iopt_turb_stab_mod,lvar)
CASE('IOPT_TURB_STAB_TKE')
   CALL conv_from_chars(cvals(1),iopt_turb_stab_tke,lvar)
CASE('IOPT_TURB_TKE_BBC')
   CALL conv_from_chars(cvals(1),iopt_turb_tke_bbc,lvar)
CASE('IOPT_TURB_TKE_SBC')
   CALL conv_from_chars(cvals(1),iopt_turb_tke_sbc,lvar)
CASE('IOPT_VADV_IMPL')
   CALL conv_from_chars(cvals(1),iopt_vadv_impl,lvar)
CASE('IOPT_VDIF_COEF')
   CALL conv_from_chars(cvals(1),iopt_vdif_coef,lvar)
CASE('IOPT_VDIF_IMPL')
   CALL conv_from_chars(cvals(1),iopt_vdif_impl,lvar)
CASE('IOPT_VDIF_ROT')
   CALL conv_from_chars(cvals(1),iopt_vdif_rot,lvar)
CASE('IOPT_WAVES'); CALL conv_from_chars(cvals(1),iopt_waves,lvar)
CASE('IOPT_WAVES_COUPLE')
   CALL conv_from_chars(cvals(1),iopt_waves_couple,lvar)
CASE('IOPT_WAVES_CURR')
   CALL conv_from_chars(cvals(1),iopt_waves_curr,lvar)
CASE('IOPT_WAVES_DISSIP')
   CALL conv_from_chars(cvals(1),iopt_waves_dissip,lvar)
CASE('IOPT_WAVES_EXTRAPOL')
   CALL conv_from_chars(cvals(1),iopt_waves_extrapol,lvar)
CASE('IOPT_WAVES_FORM')
   CALL conv_from_chars(cvals(1),iopt_waves_form,lvar)
CASE('IOPT_WAVES_PRES')
   CALL conv_from_chars(cvals(1),iopt_waves_pres,lvar)
CASE('IOPT_WEIBAR'); CALL conv_from_chars(cvals(1),iopt_weibar,lvar)

!
!2. Model grid
!-------------
!

CASE('NC'); CALL conv_from_chars(cvals(1),nc,lvar)
CASE('NR'); CALL conv_from_chars(cvals(1),nr,lvar)
CASE('NZ'); CALL conv_from_chars(cvals(1),nz,lvar)
CASE('NOSBU'); CALL conv_from_chars(cvals(1),nosbu,lvar)
CASE('NOSBV'); CALL conv_from_chars(cvals(1),nosbv,lvar)
CASE('NRVBU'); CALL conv_from_chars(cvals(1),nrvbu,lvar)
CASE('NRVBV'); CALL conv_from_chars(cvals(1),nrvbv,lvar)

!
!3. Dynamic masks
!----------------
!

CASE('FLD_MASK')
   lvar_310: DO lvar=1,nofldmasks
      CALL conv_from_chars(cvals(lvar),fld_mask(lvar),lvar)
   ENDDO lvar_310

!
!4. Number of nested sub-grids and relaxation zones
!--------------------------------------------------
!

CASE('NONESTSETS'); CALL conv_from_chars(cvals(1),nonestsets,lvar)
CASE('NORLXZONES'); CALL conv_from_chars(cvals(1),norlxzones,lvar)

!
!5. Parameters for parallelisation
!---------------------------------
!

CASE('NPROCSX'); CALL conv_from_chars(cvals(1),nprocsx,lvar)
CASE('NPROCSY'); CALL conv_from_chars(cvals(1),nprocsy,lvar)

!
!6. Date/time parameters
!-----------------------
!

CASE('CSTARTDATETIME')
   CALL conv_from_chars(cvals(1),CStartDateTime,lvar)
CASE('CENDDATETIME'); CALL conv_from_chars(cvals(1),CEndDateTime,lvar)
CASE('DELT2D'); CALL conv_from_chars(cvals(1),delt2d,lvar)
CASE('IC3D'); CALL conv_from_chars(cvals(1),ic3d,lvar)
CASE('ICCVT'); CALL conv_from_chars(cvals(1),iccvt,lvar)
CASE('ICNODAL'); CALL conv_from_chars(cvals(1),icnodal,lvar)
CASE('TIME_ZONE'); CALL conv_from_chars(cvals(1),time_zone,lvar)
CASE('NTOBCRLX'); CALL conv_from_chars(cvals(1),ntobcrlx,lvar)
CASE('RETURN_TIME')
   CALL check_cif_lbound_novars(numvars,nz)
   lvar_610: DO lvar=1,nz
      CALL conv_from_chars(cvals(lvar),return_time(lvar),lvar)
   ENDDO lvar_610

!
!7. Model constants
!------------------
!
!7.1 Physical parameters
!-----------------------
!

CASE('ATMPRES_REF'); CALL conv_from_chars(cvals(1),atmpres_ref,lvar)
CASE('BDRAGCOEF_CST'); CALL conv_from_chars(cvals(1),bdragcoef_cst,lvar)
CASE('BDRAGLIN'); CALL conv_from_chars(cvals(1),bdraglin,lvar)
CASE('B_SH'); CALL conv_from_chars(cvals(1),b_SH,lvar)
CASE('CCHARNO'); CALL conv_from_chars(cvals(1),ccharno,lvar)
CASE('CDSPARS')
   CALL check_cif_lbound_novars(numvars,4)
   lvar_710: DO lvar=1,4
      CALL conv_from_chars(cvals(lvar),cdspars(lvar),lvar)
   ENDDO lvar_710
CASE('CES_SCST'); CALL conv_from_chars(cvals(1),ces_scst,lvar)
CASE('CES_UCST'); CALL conv_from_chars(cvals(1),ces_ucst,lvar)
CASE('CHS_SCST'); CALL conv_from_chars(cvals(1),chs_scst,lvar)
CASE('CHS_UCST'); CALL conv_from_chars(cvals(1),chs_ucst,lvar)
CASE('CGRAVRATIO'); CALL conv_from_chars(cvals(1),cgravratio,lvar)
CASE('CKAR'); CALL conv_from_chars(cvals(1),ckar,lvar)
CASE('DCRIT_FLD'); CALL conv_from_chars(cvals(1),dcrit_fld,lvar)
CASE('DEPMEAN_CST'); CALL conv_from_chars(cvals(1),depmean_cst,lvar)
CASE('DEPMEAN_FLAG'); CALL conv_from_chars(cvals(1),depmean_flag,lvar)
CASE('DISTRLX_OBC'); CALL conv_from_chars(cvals(1),distrlx_obc,lvar)
CASE('DLAT_REF'); CALL conv_from_chars(cvals(1),dlat_ref,lvar)
CASE('DLON_REF'); CALL conv_from_chars(cvals(1),dlon_ref,lvar)
CASE('DLON_REF_ANAL'); CALL conv_from_chars(cvals(1),dlon_ref_anal,lvar)
CASE('DLON_REF_OBC'); CALL conv_from_chars(cvals(1),dlon_ref_obc,lvar)
CASE('DL_BB'); CALL conv_from_chars(cvals(1),dl_BB,lvar)
CASE('DMIN_FLD'); CALL conv_from_chars(cvals(1),dmin_fld,lvar)
CASE('DTHD_FLD'); CALL conv_from_chars(cvals(1),dthd_fld,lvar)
CASE('DU_BB'); CALL conv_from_chars(cvals(1),du_BB,lvar)
CASE('DZETARESID_CONV')
   CALL conv_from_chars(cvals(1),dzetaresid_conv,lvar)
CASE('GACC_REF'); CALL conv_from_chars(cvals(1),gacc_ref,lvar)
CASE('HDIFMOM_CST'); CALL conv_from_chars(cvals(1),hdifmom_cst,lvar)
CASE('HDIFSCAL_CST'); CALL conv_from_chars(cvals(1),hdifscal_cst,lvar)
CASE('KINVISC_CST'); CALL conv_from_chars(cvals(1),kinvisc_cst,lvar)
CASE('MAXITSIMP'); CALL conv_from_chars(cvals(1),maxitsimp,lvar)
CASE('MG_TOL'); CALL conv_from_chars(cvals(1),mg_tol,lvar)
CASE('NDWEXP'); CALL conv_from_chars(cvals(1),ndwexp,lvar)
CASE('NOMGLEVELS'); CALL conv_from_chars(cvals(1),nomglevels,lvar)
CASE('NOMGITERATIONS')
   CALL conv_from_chars(cvals(1),nomgiterations,lvar)
CASE('NOPRESWEEPS'); CALL conv_from_chars(cvals(1),nopresweeps,lvar)
CASE('NOPOSTSWEEPS'); CALL conv_from_chars(cvals(1),nopostsweeps,lvar)
CASE('NOSMOOTHSTEPS'); CALL conv_from_chars(cvals(1),nosmoothsteps,lvar)
CASE('OPTATTCOEF1_CST')
   CALL conv_from_chars(cvals(1),optattcoef1_cst,lvar)
CASE('OPTATTCOEF2_CST')
   CALL conv_from_chars(cvals(1),optattcoef2_cst,lvar)
CASE('OPT_FRAC'); CALL conv_from_chars(cvals(1),opt_frac,lvar)
CASE('RHO_AIR'); CALL conv_from_chars(cvals(1),rho_air,lvar)
CASE('SAL_REF'); CALL conv_from_chars(cvals(1),sal_ref,lvar)
CASE('SIGSTAR_DJ'); CALL conv_from_chars(cvals(1),sigstar_DJ,lvar)
CASE('SIG0_DJ'); CALL conv_from_chars(cvals(1),sig0_DJ,lvar)
CASE('RHO_CRIT_ISO'); CALL conv_from_chars(cvals(1),rho_crit_iso,lvar)   
CASE('SKEWDIFF_CST'); CALL conv_from_chars(cvals(1),skewdiff_cst,lvar)
CASE('SLOPEMAX_ISO'); CALL conv_from_chars(cvals(1),slopemax_iso,lvar)
CASE('SMAG_COEF_MOM'); CALL conv_from_chars(cvals(1),smag_coef_mom,lvar)
CASE('SMAG_COEF_SCAL')
   CALL conv_from_chars(cvals(1),smag_coef_scal,lvar)
CASE('SPECHEAT'); CALL conv_from_chars(cvals(1),specheat,lvar)
CASE('SST_REF'); CALL conv_from_chars(cvals(1),sst_ref,lvar)
CASE('TEMP_MIN'); CALL conv_from_chars(cvals(1),temp_min,lvar)
CASE('TEMP_REF'); CALL conv_from_chars(cvals(1),temp_ref,lvar)
CASE('THETA_COR'); CALL conv_from_chars(cvals(1),theta_cor,lvar)
CASE('THETA_SH'); CALL conv_from_chars(cvals(1),theta_SH,lvar)
CASE('THETA_SUR'); CALL conv_from_chars(cvals(1),theta_sur,lvar)
CASE('THETA_VADV'); CALL conv_from_chars(cvals(1),theta_vadv,lvar)
CASE('THETA_VDIF'); CALL conv_from_chars(cvals(1),theta_vdif,lvar)
CASE('UR_MG'); CALL conv_from_chars(cvals(1),ur_mg,lvar)
CASE('UR_SMOOTH'); CALL conv_from_chars(cvals(1),ur_smooth,lvar)
CASE('VDIFMOM_CST'); CALL conv_from_chars(cvals(1),vdifmom_cst,lvar)
CASE('VDIFSCAL_CST'); CALL conv_from_chars(cvals(1),vdifscal_cst,lvar)
CASE('WAVE_PENETRATION_BED')
  CALL conv_from_chars(cvals(1),wave_penetration_bed,lvar)
CASE('WAVE_PENETRATION_SURF')
   CALL conv_from_chars(cvals(1),wave_penetration_surf,lvar)
CASE('WAVETHICK_CST')
  CALL conv_from_chars(cvals(1),wavethick_cst,lvar)
CASE('ZREF_TEMP'); CALL conv_from_chars(cvals(1),zref_temp,lvar)
CASE('ZREF_WIND'); CALL conv_from_chars(cvals(1),zref_wind,lvar)
CASE('ZBTOZ0LIM'); CALL conv_from_chars(cvals(1),zbtoz0lim,lvar)
CASE('ZROUGH_CST'); CALL conv_from_chars(cvals(1),zrough_cst,lvar)

!
!7.2 Tidal constituents
!----------------------
!
!---open boundaries
CASE('NCONOBC'); CALL conv_from_chars(cvals(1),nconobc,lvar)
CASE('INDEX_OBC')
   CALL check_cif_lbound_novars(numvars,nconobc)
   lvar_721: DO lvar=1,nconobc
      CALL conv_from_chars(cvals(lvar),index_obc(lvar),lvar)
   ENDDO lvar_721
CASE('NQSECOBU'); CALL conv_from_chars(cvals(1),nqsecobu,lvar)
CASE('NQSECOBV'); CALL conv_from_chars(cvals(1),nqsecobv,lvar)

!---tidal force
CASE('NCONASTRO'); CALL conv_from_chars(cvals(1),nconastro,lvar)
CASE('INDEX_ASTRO')
   CALL check_cif_lbound_novars(numvars,nconastro)
   lvar_722: DO lvar=1,nconastro
      CALL conv_from_chars(cvals(lvar),index_astro(lvar),lvar)
   ENDDO lvar_722

!
!7.3 Structure module parameters
!-------------------------------
!

CASE('NUMDRY'); CALL conv_from_chars(cvals(1),numdry,lvar)
CASE('NUMTHINU'); CALL conv_from_chars(cvals(1),numthinu,lvar)
CASE('NUMTHINV'); CALL conv_from_chars(cvals(1),numthinv,lvar)
CASE('NUMWBARU'); CALL conv_from_chars(cvals(1),numwbaru,lvar)
CASE('NUMWBARV'); CALL conv_from_chars(cvals(1),numwbarv,lvar)
CASE('WBARRLXU'); CALL conv_from_chars(cvals(1),wbarrlxu,lvar)
CASE('WBARRLXV'); CALL conv_from_chars(cvals(1),wbarrlxv,lvar)
CASE('NUMDIS'); CALL conv_from_chars(cvals(1),numdis,lvar)

!
!7.4 Turbulence parameters
!-------------------------
!

CASE('ALPHA_BLACK'); CALL conv_from_chars(cvals(1),alpha_black,lvar)
CASE('ALPHA_MA'); CALL conv_from_chars(cvals(1),alpha_ma,lvar)
CASE('ALPHA_PP'); CALL conv_from_chars(cvals(1),alpha_pp,lvar)
CASE('BETA_MA'); CALL conv_from_chars(cvals(1),beta_ma,lvar)
CASE('BETA_XING'); CALL conv_from_chars(cvals(1),beta_xing,lvar)
CASE('CNU_AD'); CALL conv_from_chars(cvals(1),cnu_ad,lvar)
CASE('C_SK'); CALL conv_from_chars(cvals(1),c_sk,lvar)
CASE('C1_EPS'); CALL conv_from_chars(cvals(1),c1_eps,lvar)
CASE('C2_EPS'); CALL conv_from_chars(cvals(1),c2_eps,lvar)
CASE('C31_EPS'); CALL conv_from_chars(cvals(1),c31_eps,lvar)
CASE('C32_EPS'); CALL conv_from_chars(cvals(1),c32_eps,lvar)
CASE('DELTA1_AD'); CALL conv_from_chars(cvals(1),delta1_ad,lvar)
CASE('DELTA2_AD'); CALL conv_from_chars(cvals(1),delta2_ad,lvar)
CASE('DISSIPMIN'); CALL conv_from_chars(cvals(1),dissipmin,lvar)
CASE('EXPMOM_MA'); CALL conv_from_chars(cvals(1),expmom_ma,lvar)
CASE('EXPMOM_PP'); CALL conv_from_chars(cvals(1),expmom_pp,lvar)
CASE('EXPSCAL_MA'); CALL conv_from_chars(cvals(1),expscal_ma,lvar)
CASE('E1_MY'); CALL conv_from_chars(cvals(1),e1_my,lvar)
CASE('E2_MY'); CALL conv_from_chars(cvals(1),e2_my,lvar)
CASE('E3_MY'); CALL conv_from_chars(cvals(1),e3_my,lvar)
CASE('K1_AD'); CALL conv_from_chars(cvals(1),k1_ad,lvar)
CASE('K2_AD'); CALL conv_from_chars(cvals(1),k2_ad,lvar)
CASE('LAMBDA_AD'); CALL conv_from_chars(cvals(1),lambda_ad,lvar)
CASE('OMEGA1_AD'); CALL conv_from_chars(cvals(1),omega1_ad,lvar)
CASE('RICCRIT_IW'); CALL conv_from_chars(cvals(1),riccrit_iw,lvar)
CASE('R1_AD'); CALL conv_from_chars(cvals(1),r1_ad,lvar)
CASE('R2_AD'); CALL conv_from_chars(cvals(1),r2_ad,lvar)
CASE('SIGMA_K'); CALL conv_from_chars(cvals(1),sigma_k,lvar)
CASE('SKEPS'); CALL conv_from_chars(cvals(1),skeps,lvar)
CASE('SQ_MY'); CALL conv_from_chars(cvals(1),sq_my,lvar)
CASE('TKELIM'); CALL conv_from_chars(cvals(1),tkelim,lvar)
CASE('TKEMIN'); CALL conv_from_chars(cvals(1),tkemin,lvar)
CASE('VBMOM_PP'); CALL conv_from_chars(cvals(1),vbmom_pp,lvar)
CASE('VBSCAL_PP'); CALL conv_from_chars(cvals(1),vbscal_pp,lvar)
CASE('VDIFMOM_IW'); CALL conv_from_chars(cvals(1),vdifmom_iw,lvar)
CASE('VDIFSCAL_IW'); CALL conv_from_chars(cvals(1),vdifscal_iw,lvar)
CASE('VDIFSHEAR_IW'); CALL conv_from_chars(cvals(1),vdifshear_iw,lvar)
CASE('VMAXMOM_MA'); CALL conv_from_chars(cvals(1),vmaxmom_ma,lvar)
CASE('VMAXSCAL_MA'); CALL conv_from_chars(cvals(1),vmaxscal_ma,lvar)
CASE('VMAX_PP'); CALL conv_from_chars(cvals(1),vmax_pp,lvar)
CASE('V0DIF_MA'); CALL conv_from_chars(cvals(1),v0dif_ma,lvar)
CASE('V0DIF_PP'); CALL conv_from_chars(cvals(1),v0dif_pp,lvar)
CASE('WFLTKE'); CALL conv_from_chars(cvals(1),wfltke,lvar)
CASE('ZLMIXMIN'); CALL conv_from_chars(cvals(1),zlmixmin,lvar)
CASE('ZROUGH_BOT'); CALL conv_from_chars(cvals(1),zrough_bot,lvar)
CASE('ZROUGH_SUR'); CALL conv_from_chars(cvals(1),zrough_sur,lvar)

!
!8. I/O specifiers
!----------------
!
!---restart times
CASE('NORESTARTS'); CALL conv_from_chars(cvals(1),norestarts,lvar)
CASE('NTRESTART')
   CALL check_cif_lbound_novars(numvars,norestarts)
   lvar_810: DO lvar=1,norestarts
      CALL conv_from_chars(cvals(lvar),ntrestart(lvar),lvar)
   ENDDO lvar_810

!---input/output titles
CASE('INTITLE'); CALL conv_from_chars(cvals(1),intitle,lvar)
CASE('OUTTITLE'); CALL conv_from_chars(cvals(1),outtitle,lvar)

!---CF global attributes
CASE('INSTITUTION_CF');
   CALL conv_from_chars(cvals(1),institution_CF,lvar)
CASE('COMMENT_CF'); CALL conv_from_chars(cvals(1),comment_CF,lvar)
CASE('REFERENCES_CF'); CALL conv_from_chars(cvals(1),references_CF,lvar)

!---I/O parameters
CASE('MAXWAITSECS'); CALL conv_from_chars(cvals(1),maxwaitsecs,lvar)
CASE('NOWAITSECS'); CALL conv_from_chars(cvals(1),nowaitsecs,lvar)
CASE('NRECUNIT'); CALL conv_from_chars(cvals(1),nrecunit,lvar)

!---user-output parameters
CASE('NOSETSTSR'); CALL conv_from_chars(cvals(1),nosetstsr,lvar)
CASE('NOSTATSTSR'); CALL conv_from_chars(cvals(1),nostatstsr,lvar)
CASE('NOVARSTSR'); CALL conv_from_chars(cvals(1),novarstsr,lvar)
CASE('NOSETSAVR'); CALL conv_from_chars(cvals(1),nosetsavr,lvar)
CASE('NOSTATSAVR'); CALL conv_from_chars(cvals(1),nostatsavr,lvar)
CASE('NOVARSAVR'); CALL conv_from_chars(cvals(1),novarsavr,lvar)
CASE('NOSETSANAL'); CALL conv_from_chars(cvals(1),nosetsanal,lvar)
CASE('NOFREQSANAL'); CALL conv_from_chars(cvals(1),nofreqsanal,lvar)
CASE('NOSTATSANAL'); CALL conv_from_chars(cvals(1),nostatsanal,lvar)
CASE('NOVARSANAL'); CALL conv_from_chars(cvals(1),novarsanal,lvar)

!
!9. Attributes of forcing files
!------------------------------
!

CASE('MODFILES')
   CALL check_cif_lbound_novars(numvars,15)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),fdesc,lvar)
   n_910: DO n=1,MaxIOTypes
      IF (fdesc.EQ.modfiles_desc(n)) idesc = n
   ENDDO n_910
   CALL error_array_index(idesc,'modfiles',1,MaxIOTypes,1)
   lvar = 2
   CALL conv_from_chars(cvals(lvar),ifil,lvar)
   CALL error_array_index(ifil,'modfiles',1,MaxIOFiles,2)
   lvar = 3
   CALL conv_from_chars(cvals(lvar),iotype,lvar)
   CALL error_array_index(iotype,'modfiles',1,2,3)
   CALL conv_from_chars_infiles(cvals(4:15),modfiles(idesc,ifil,iotype),&
                              & lvar)

!
!10. Data grid attributes
!------------------------
!

CASE('SURFACEGRIDS')
   CALL check_cif_lbound_novars(numvars,14)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),idesc,lvar)
   CALL error_array_index(idesc,'modfiles',1,MaxGridTypes,1)
   CALL conv_from_chars_surfgrd(cvals(2:15),surfacegrids(idesc,1),lvar)

END SELECT

RETURN


END SUBROUTINE assign_cif_vars_mod

!========================================================================

SUBROUTINE assign_cif_vars_mon(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_mon* convert string data from an input line in the CIF
!                       block with general/monitoring parameters to the
!                       appropriate numeric or non-numeric format 
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11
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
USE paralpars
USE switches
USE syspars
USE cif_routines, ONLY: conv_from_chars

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
INTEGER :: lvar


lvar = 1

SELECT CASE (TRIM(cname))

!---cold/warm start
CASE('COLD_START'); CALL conv_from_chars(cvals(1),cold_start,lvar)

!---log files
CASE('LEVPROCS_INI')
   lvar_110: DO lvar=1,numvars
      CALL conv_from_chars(cvals(lvar),levprocs_ini(lvar),lvar)
   ENDDO lvar_110
CASE('LEVPROCS_RUN')
   lvar_120: DO lvar=1,numvars
      CALL conv_from_chars(cvals(lvar),levprocs_run(lvar),lvar)
   ENDDO lvar_120
CASE('INILOG_FILE'); CALL conv_from_chars(cvals(1),inilog_file,lvar)
CASE('RUNLOG_FILE'); CALL conv_from_chars(cvals(1),runlog_file,lvar)
CASE('RUNLOG_COUNT'); CALL conv_from_chars(cvals(1),runlog_count,lvar)

!---error files
CASE('MAXERRORS'); CALL conv_from_chars(cvals(1),maxerrors,lvar)
CASE('LEVPROCS_ERR')
   lvar_130: DO lvar=1,numvars
      CALL conv_from_chars(cvals(lvar),levprocs_err(lvar),lvar)
   ENDDO lvar_130
CASE('ERRLOG_FILE'); CALL conv_from_chars(cvals(1),errlog_file,lvar)

!---warning file
CASE('WARNING'); CALL conv_from_chars(cvals(1),warning,lvar)
CASE('WARLOG_FILE'); CALL conv_from_chars(cvals(1),warlog_file,lvar)

!---monitoring
CASE ('MONLOG'); CALL conv_from_chars(cvals(1),monlog,lvar)
CASE('MONLOG_FILE'); CALL conv_from_chars(cvals(1),monlog_file,lvar)
CASE ('SEDLOG'); CALL conv_from_chars(cvals(1),sedlog,lvar)
CASE ('SEDLOG_FILE'); CALL conv_from_chars(cvals(1),sedlog_file,lvar)

!---timer report
CASE('LEVTIMER'); CALL conv_from_chars(cvals(1),levtimer,lvar)
CASE('TIMING_FILE'); CALL conv_from_chars(cvals(1),timing_file,lvar)
CASE('TIMER_FORMAT'); CALL conv_from_chars(cvals(1),timer_format,lvar)

!---number of processes for each model component
CASE('NPROCSCOH'); CALL conv_from_chars(cvals(1),nprocscoh,lvar)
CASE('NPROCSWAV'); CALL conv_from_chars(cvals(1),nprocswav,lvar)
   
!---switches for model coupling
CASE('IOPT_PART_MODEL')
   CALL conv_from_chars(cvals(1),iopt_part_model,lvar)
CASE('IOPT_WAVES_MODEL')
   CALL conv_from_chars(cvals(1),iopt_waves_model,lvar)

END SELECT


RETURN

END SUBROUTINE assign_cif_vars_mon

!========================================================================

SUBROUTINE assign_cif_vars_tsout(cname,cvals,numvars)
!************************************************************************
!
! *assign_cif_vars_tsout* convert string data from an input line in the CIF
!                         block with time series parameters to the appropriate
!                         numeric or non-numeric format   
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11
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
!                conv_from_chars_gridpars, conv_from_chars_outfiles,
!                conv_from_chars_statlocs, conv_from_chars_outvars,
!                error_array_index
!
!************************************************************************
!
USE iopars
USE syspars
USE cif_routines, ONLY: check_cif_lbound_novars, conv_from_chars, &
                      & conv_from_chars_gridpars, conv_from_chars_outfiles, &
                      & conv_from_chars_statlocs, conv_from_chars_outvars
USE error_routines, ONLY: error_array_index

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
INTEGER :: iset, istat, ivar, lvar, nostats


lvar = 1

SELECT CASE (TRIM(cname))

!---variable attributes
CASE('TSRVARS')
   CALL check_cif_lbound_novars(numvars,12)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),ivar,lvar)
   CALL error_array_index(ivar,'tsrvars',1,novarstsr,1)
   CALL conv_from_chars_outvars(cvals(2:12),tsrvars(ivar),lvar)

!---variable indices
CASE('IVARSTSR')
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'ivarstsr',1,nosetstsr,1)
   lvar_110: DO lvar=2,numvars
      CALL conv_from_chars(cvals(lvar),ivarstsr(iset,lvar-1),lvar)
   ENDDO lvar_110

!---file attributes
CASE('TSR0D')
   CALL check_cif_lbound_novars(numvars,7)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'tsr0d',1,nosetstsr,1)
   CALL conv_from_chars_outfiles(cvals(2:7),tsr0d(iset),lvar)

CASE('TSR2D')
   CALL check_cif_lbound_novars(numvars,7)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'tsr2d',1,nosetstsr,1)
   CALL conv_from_chars_outfiles(cvals(2:7),tsr2d(iset),lvar)

CASE('TSR3D')
   CALL check_cif_lbound_novars(numvars,7)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'tsr3d',1,nosetstsr,1)
   CALL conv_from_chars_outfiles(cvals(2:7),tsr3d(iset),lvar)

!---output grid (space/time) attributes
CASE('TSRGPARS')
   CALL check_cif_lbound_novars(numvars,20)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   CALL error_array_index(iset,'tsrgpars',1,nosetstsr,1)
   CALL conv_from_chars_gridpars(cvals(2:20),tsrgpars(iset),lvar)

!---station attributes
CASE('TSRSTATLOCS')
   CALL check_cif_lbound_novars(numvars,4)
   lvar = 1
   CALL conv_from_chars(cvals(lvar),istat,lvar)
   CALL error_array_index(istat,'tsrstatlocs',1,nostatstsr,1)
   CALL conv_from_chars_statlocs(cvals(2:4),tsrstatlocs(istat),lvar)

!---station labels
CASE('LSTATSTSR')
   lvar = 1
   CALL conv_from_chars(cvals(lvar),iset,lvar)
   nostats = tsrgpars(iset)%nostats
   CALL error_array_index(iset,'lstatstsr',1,nosetstsr,1)
   CALL check_cif_lbound_novars(numvars,nostats+1)
   lvar_120: DO lvar=2,nostats+1
      CALL conv_from_chars(cvals(lvar),lstatstsr(iset,lvar-1),lvar)
   ENDDO lvar_120

END SELECT


RETURN

END SUBROUTINE assign_cif_vars_tsout

!========================================================================

SUBROUTINE read_cif_params(iddesc)
!************************************************************************
!
! *read_cif_params* Read model setup data from the CIF
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - harmonic_analysis, initialise_model,
!                   simulation_start, time_averages_init,
!                   time_series
!
! External calls - assign_cif_vars_anal, assign_cif_vars_avrgd,
!                  assign_cif_vars_bio, assign_cif_vars_dar,
!                  assign_cif_vars_harm, assign_vars_mod,
!                  assign_cif_vars_mon, assign_vars_morph,
!                  assign_cif_vars_part, assign_vars_sed,
!                  assign_cif_vars_tsout, assign_vars_tspart
!
! Module calls - error_file, read_cif_line
!
!************************************************************************
!
USE iopars
USE syspars
USE cif_routines, ONLY: read_cif_line
USE error_routines, ONLY: error_file 
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER CIF file key id
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: end_block, search
CHARACTER (LEN=lencifline) :: cline
CHARACTER (LEN=lenname) :: cname
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: iunit, numvars


procname(pglev+1) = 'read_cif_params'
CALL log_timer_in()

!
!1. Initialise parameters
!------------------------
!

iunit = ciffile%iunit
REWIND(iunit)

end_block = .FALSE.
cline = ''

!
!2. Search block header name in CIF file
!---------------------------------------
! 

search = .FALSE.
DO WHILE (.NOT.search)
   READ (iunit,'(A)',END=1001) cline
   ciflinenum = ciflinenum + 1
   IF (LEN_TRIM(cline).GT.0.AND.cline(1:1).EQ.cifend.AND.&
        & TRIM(ADJUSTL(cline(2:3))).EQ.TRIM(cif_block_names(iddesc)(1:2))) THEN
      search = .TRUE.
   ENDIF
ENDDO

!
!3. Read with assignment
!-----------------------
!

DO WHILE (.NOT.end_block)
   READ (iunit,'(A)',END=1002) cline
   ciflinenum = ciflinenum + 1
   IF (LEN_TRIM(cline).EQ.0) CYCLE
   cline = ADJUSTL(cline)
   IF (cline(1:1).EQ.cifend) THEN
      end_block = .TRUE.
   ELSE
      CALL read_cif_line(cline,cvals,numvars,cname=cname)
      SELECT CASE (iddesc)
         CASE (icif_mon)
            CALL assign_cif_vars_mon(cname,cvals,numvars)
         CASE (icif_mod)
            CALL assign_cif_vars_mod(cname,cvals,numvars)
         CASE (icif_bio)
            CALL assign_cif_vars_bio(cname,cvals,numvars)
         CASE (icif_sed)
            CALL assign_cif_vars_sed(cname,cvals,numvars)
         CASE (icif_morph) 
            CALL assign_cif_vars_morph(cname,cvals,numvars)
         CASE (icif_dar)
            CALL assign_cif_vars_dar(cname,cvals,numvars)
         CASE (icif_part)
            CALL assign_cif_vars_part(cname,cvals,numvars)
         CASE (icif_tsout)
            CALL assign_cif_vars_tsout(cname,cvals,numvars)
         CASE (icif_avrgd)
            CALL assign_cif_vars_avrgd(cname,cvals,numvars)
         CASE (icif_harm)
            CALL assign_cif_vars_harm(cname,cvals,numvars)
         CASE (icif_anal)
            CALL assign_cif_vars_anal(cname,cvals,numvars)
         CASE (icif_tspart)
            CALL assign_cif_vars_tspart(cname,cvals,numvars)
      END SELECT
   ENDIF
ENDDO

CALL log_timer_out()


RETURN

1001 nerrs = nerrs + 1
IF (errchk.AND.nerrs.LE.maxerrors) THEN
   WRITE (ioerr,'(A)') 'Cannot find CIF block '//TRIM(cif_block_names(iddesc))
   CALL error_file(ierrno_fend,filepars=ciffile)
ENDIF

1002 nerrs = nerrs + 1
CALL error_file(ierrno_fend,filepars=ciffile)

END SUBROUTINE read_cif_params

!========================================================================

SUBROUTINE write_cif_vars_anal
!************************************************************************
!
! *write_cif_vars_anal* Write the CIF block with parameters for harmonic output
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - conv_to_chars, conv_to_chars_gridpars, conv_to_chars_outfiles,
!                conv_to_chars_outvars, conv_to_chars_statlocs, write_cif_line
!
!************************************************************************
!
USE iopars
USE paralpars  
USE syspars
USE cif_routines, ONLY: conv_to_chars, conv_to_chars_gridpars, &
                      & conv_to_chars_outfiles, conv_to_chars_outvars, &
                      & conv_to_chars_statlocs, write_cif_line
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: ifreq, iset, istat, iunit, ivar, lvar, nostats


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_cif_vars_anal'
CALL log_timer_in()

!
!1. Write block header
!---------------------
!

iunit = ciffile%iunit
WRITE (iunit,'(A)') cifend//TRIM(cif_block_names(icif_anal))

!
!2. Variable attributes
!----------------------
!

ivar_210: DO ivar=1,novarsanal
   CALL conv_to_chars(cvals(1),ivar)
   CALL conv_to_chars_outvars(cvals(2:12),analvars(ivar))
   CALL write_cif_line(cvals(1:12),'ANALVARS')
ENDDO ivar_210

!
!3. Variable indices
!-------------------
!

iset_310: DO iset=1,nosetsanal
   lvar = 1
   CALL conv_to_chars(cvals(1),iset)
   ivar_311: DO ivar=1,novarsanal
      IF (ivarsanal(iset,ivar).GT.0) THEN
         lvar = lvar + 1
         CALL conv_to_chars(cvals(lvar),ivarsanal(iset,ivar))
      ENDIF
   ENDDO ivar_311
   CALL write_cif_line(cvals(1:lvar),'IVARSANAL')
ENDDO iset_310

!
!4. Elliptic variables
!---------------------
!
!---variable numbers
iset_410: DO iset=1,nosetsanal
   lvar = 1
   CALL conv_to_chars(cvals(1),iset)
   ivar_411: DO ivar=1,14
      IF (ivarsell(iset,ivar).GT.0) THEN
         lvar = lvar + 1
         CALL conv_to_chars(cvals(lvar),ivarsell(iset,ivar))
      ENDIF
   ENDDO ivar_411
   CALL write_cif_line(cvals(1:lvar),'IVARSELL')
ENDDO iset_410

!---components of elliptic vector (2-D)
iset_420: DO iset=1,nosetsanal
   CALL conv_to_chars(cvals(1),iset)
   CALL conv_to_chars(cvals(2:3),ivecell2d(iset,:))
   CALL write_cif_line(cvals(1:3),'IVECELL2D')
ENDDO iset_420

!---components of elliptic vector (3-D)
iset_430: DO iset=1,nosetsanal
   CALL conv_to_chars(cvals(1),iset)
   CALL conv_to_chars(cvals(2:3),ivecell3d(iset,:))
   CALL write_cif_line(cvals(1:3),'IVECELL3D')
ENDDO iset_430

!
!5. File attributes
!------------------
!

iset_510: DO iset=1,nosetsanal
   CALL conv_to_chars(cvals(1),iset)
   IF (res0d(iset)%defined) THEN
      CALL conv_to_chars_outfiles(cvals(2:7),res0d(iset))
      CALL write_cif_line(cvals(1:7),'RES0D')
   ENDIF
   IF (res2d(iset)%defined) THEN
      CALL conv_to_chars_outfiles(cvals(2:7),res2d(iset))
      CALL write_cif_line(cvals(1:7),'RES2D')
   ENDIF
   IF (res3d(iset)%defined) THEN
      CALL conv_to_chars_outfiles(cvals(2:7),res3d(iset))
      CALL write_cif_line(cvals(1:7),'RES3D')
   ENDIF
   
   ifreq_511: DO ifreq=1,nofreqsharm(iset)
      CALL conv_to_chars(cvals(2),ifreq)
      IF (amp2d(iset,ifreq)%defined) THEN
         CALL conv_to_chars_outfiles(cvals(3:8),amp2d(iset,ifreq))
         CALL write_cif_line(cvals(1:8),'AMP2D')
      ENDIF
      IF (amp3d(iset,ifreq)%defined) THEN
         CALL conv_to_chars_outfiles(cvals(3:8),amp3d(iset,ifreq))
         CALL write_cif_line(cvals(1:8),'AMP3D')
      ENDIF
      IF (pha2d(iset,ifreq)%defined) THEN
         CALL conv_to_chars_outfiles(cvals(3:8),pha2d(iset,ifreq))
         CALL write_cif_line(cvals(1:8),'PHA2D')
      ENDIF
      IF (pha3d(iset,ifreq)%defined) THEN
         CALL conv_to_chars_outfiles(cvals(3:8),pha3d(iset,ifreq))
         CALL write_cif_line(cvals(1:8),'PHA3D')
      ENDIF
      IF (ell2d(iset,ifreq)%defined) THEN
         CALL conv_to_chars_outfiles(cvals(3:8),ell2d(iset,ifreq))
         CALL write_cif_line(cvals(1:8),'ELL2D')
      ENDIF
      IF (ell3d(iset,ifreq)%defined) THEN
         CALL conv_to_chars_outfiles(cvals(3:8),ell3d(iset,ifreq))
         CALL write_cif_line(cvals(1:8),'ELL3D')
      ENDIF
   ENDDO ifreq_511
      
ENDDO iset_510

!
!6.Output grid (space/time) attributes
!---------------------------------------
!

iset_610: DO iset=1,nosetsanal
   CALL conv_to_chars(cvals(1),iset)
   CALL conv_to_chars_gridpars(cvals(2:20),analgpars(iset))
   CALL write_cif_line(cvals(1:20),'ANALGPARS')
ENDDO iset_610

!
!7. Station attributes
!---------------------
!

istat_710: DO istat=1,nostatsanal
   CALL conv_to_chars(cvals(1),istat)
   CALL conv_to_chars_statlocs(cvals(2:4),analstatlocs(istat))
   CALL write_cif_line(cvals(1:4),'ANALSTATLOCS')
ENDDO istat_710

!
!8. Station labels
!-----------------
!

iset_810: DO iset=1,nosetsanal
   CALL conv_to_chars(cvals(1),iset)
   nostats = analgpars(iset)%nostats
   IF (nostats.GT.0) THEN
      CALL conv_to_chars(cvals(2:nostats+1),lstatsanal(iset,1:nostats))
      CALL write_cif_line(cvals(1:nostats+1),'LSTATSANAL')
   ENDIF
ENDDO iset_810

WRITE (iunit,'(A)') cifend
ciflinenum = ciflinenum + 1

CALL log_timer_out()


RETURN

END SUBROUTINE write_cif_vars_anal

!========================================================================

SUBROUTINE write_cif_vars_avrgd
!************************************************************************
!
! *write_cif_vars_avrgd* Write the CIF block with parameters for time averaging
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - conv_to_chars, conv_to_chars_gridpars, conv_to_chars_outfiles,
!                conv_to_chars_outvars, conv_to_chars_statlocs, write_cif_line
!
!************************************************************************
!
USE iopars
USE paralpars
USE syspars
USE cif_routines, ONLY: conv_to_chars, conv_to_chars_gridpars, &
                      & conv_to_chars_outfiles, conv_to_chars_outvars, &
                      & conv_to_chars_statlocs, write_cif_line
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: iset, istat, iunit, ivar, lvar, nostats


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_cif_vars_avrgd'
CALL log_timer_in()

!---write block header
iunit = ciffile%iunit
WRITE (iunit,'(A)') cifend//TRIM(cif_block_names(icif_avrgd))

!---variable attributes
ivar_110: DO ivar=1,novarsavr
   CALL conv_to_chars(cvals(1),ivar)
   CALL conv_to_chars_outvars(cvals(2:12),avrvars(ivar))
   CALL write_cif_line(cvals(1:12),'AVRVARS')
ENDDO ivar_110

!---variable indices
iset_120: DO iset=1,nosetsavr
   lvar = 1
   CALL conv_to_chars(cvals(1),iset)
   ivar_121: DO ivar=1,novarsavr
      IF (ivarsavr(iset,ivar).GT.0) THEN
         lvar = lvar + 1
         CALL conv_to_chars(cvals(lvar),ivarsavr(iset,ivar))
      ENDIF
   ENDDO ivar_121
   CALL write_cif_line(cvals(1:lvar),'IVARSAVR')
ENDDO iset_120

!---file attributes
iset_130: DO iset=1,nosetsavr
   IF (avr0d(iset)%defined) THEN
      CALL conv_to_chars(cvals(1),iset)
      CALL conv_to_chars_outfiles(cvals(2:7),avr0d(iset))
      CALL write_cif_line(cvals(1:7),'AVR0D')
   ENDIF
   IF (avr2d(iset)%defined) THEN
      CALL conv_to_chars(cvals(1),iset)
      CALL conv_to_chars_outfiles(cvals(2:7),avr2d(iset))
      CALL write_cif_line(cvals(1:7),'AVR2D')
   ENDIF
   IF (avr3d(iset)%defined) THEN
      CALL conv_to_chars(cvals(1),iset)
      CALL conv_to_chars_outfiles(cvals(2:7),avr3d(iset))
      CALL write_cif_line(cvals(1:7),'AVR3D')
   ENDIF
ENDDO iset_130

!---output grid (space/time) attributes
iset_140: DO iset=1,nosetsavr
   CALL conv_to_chars(cvals(1),iset)
   CALL conv_to_chars_gridpars(cvals(2:20),avrgpars(iset))
   CALL write_cif_line(cvals(1:20),'AVRGPARS')
ENDDO iset_140

!---station attributes
istat_150: DO istat=1,nostatsavr
   CALL conv_to_chars(cvals(1),istat)
   CALL conv_to_chars_statlocs(cvals(2:4),avrstatlocs(istat))
   CALL write_cif_line(cvals(1:4),'AVRSTATLOCS')
ENDDO istat_150

!---station labels
iset_160: DO iset=1,nosetsavr
   CALL conv_to_chars(cvals(1),iset)
   nostats = avrgpars(iset)%nostats
   IF (nostats.GT.0) THEN
      CALL conv_to_chars(cvals(2:nostats+1),lstatsavr(iset,1:nostats))
      CALL write_cif_line(cvals(1:nostats+1),'LSTATSAVR')
   ENDIF
ENDDO iset_160

WRITE (iunit,'(A)') cifend
ciflinenum = ciflinenum + 1

CALL log_timer_out()


RETURN

END SUBROUTINE write_cif_vars_avrgd

!========================================================================

SUBROUTINE write_cif_vars_harm
!************************************************************************
!
! *write_cif_vars_harm* Write the CIF block with parameters for harmonic analysis
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - conv_to_chars, write_cif_line
!
!************************************************************************
!
USE iopars
USE paralpars
USE syspars
USE cif_routines, ONLY: conv_to_chars, write_cif_line
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: iset, iunit, nofreqs


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_cif_vars_harm'
CALL log_timer_in()

!---write block header
iunit = ciffile%iunit
WRITE (iunit,'(A)') cifend//TRIM(cif_block_names(icif_harm))

!---harmonic tidal index ids
CALL conv_to_chars(cvals(1:nofreqsanal),index_anal)
CALL write_cif_line(cvals(1:nofreqsanal),'INDEX_ANAL')

!---frequency names
IF (ANY(index_anal.EQ.0)) THEN
   cvals(1:nofreqsanal) = harm_freq_names
   CALL write_cif_line(cvals(1:nofreqsanal),'HARM_FREQ_NAMES')
ENDIF

!---frequency values
IF (ANY(index_anal.EQ.0)) THEN
   CALL conv_to_chars(cvals(1:nofreqsanal),harm_freq)
   CALL write_cif_line(cvals(1:nofreqsanal),'HARM_FREQ')
ENDIF

!---reference times
IF (ANY(cdate_time_ref.NE.cdatetime_undef)) THEN
   cvals(1:nosetsanal) = cdate_time_ref
   CALL write_cif_line(cvals(1:nosetsanal),'CDATE_TIME_REF')
ENDIF

!---number of requencies per set
CALL conv_to_chars(cvals(1:nosetsanal),nofreqsharm)
CALL write_cif_line(cvals(1:nosetsanal),'NOFREQSHARM')

!---frequency indices
iset_110: DO iset=1,nosetsanal
   CALL conv_to_chars(cvals(1),iset)
   nofreqs = nofreqsharm(iset)
   CALL conv_to_chars(cvals(2:nofreqs+1),ifreqsharm(iset,1:nofreqs))
   CALL write_cif_line(cvals(1:nofreqs+1),'IFREQSHARM')
ENDDO iset_110

!---time step for harmonic analysis
CALL conv_to_chars(cvals(1:nosetsanal),icanal)
CALL write_cif_line(cvals(1:nosetsanal),'ICANAL')

WRITE (iunit,'(A)') cifend
ciflinenum = ciflinenum + 1

CALL log_timer_out()


RETURN

END SUBROUTINE write_cif_vars_harm

!========================================================================

SUBROUTINE write_cif_vars_mod
!************************************************************************
!
! *write_cif_vars_mod* Write the CIF block with model setup parameters
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11.2
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - conv_to_chars, conv_to_chars_infiles, conv_to_chars_surfgrd,
!                write_cif_line,
!
!************************************************************************
!
USE gridpars
USE iopars
USE nestgrids
USE obconds
USE paralpars
USE physpars
USE relaxation
USE structures
USE switches
USE syspars
USE tide
USE timepars
USE turbpars
USE cif_routines, ONLY: conv_to_chars, conv_to_chars_infiles, &
                      & conv_to_chars_surfgrd, write_cif_line
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: idesc, ifil, iotype, iunit


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_cif_vars_mod'
CALL log_timer_in()

!
!1. Write block header
!--------------------
!

iunit = ciffile%iunit
WRITE (iunit,'(A)') cifend//TRIM(cif_block_names(icif_mod))

!
!2. Switches
!-----------
!

CALL conv_to_chars(cvals(1),iopt_adv_scal)
CALL write_cif_line(cvals(1:1),'IOPT_ADV_SCAL')
CALL conv_to_chars(cvals(1),iopt_adv_turb)
CALL write_cif_line(cvals(1:1),'IOPT_ADV_TURB')
CALL conv_to_chars(cvals(1),iopt_adv_tvd)
CALL write_cif_line(cvals(1:1),'IOPT_ADV_TVD')
CALL conv_to_chars(cvals(1),iopt_adv_2D)
CALL write_cif_line(cvals(1:1),'IOPT_ADV_2D')
CALL conv_to_chars(cvals(1),iopt_adv_3D)
CALL write_cif_line(cvals(1:1),'IOPT_ADV_3D')
CALL conv_to_chars(cvals(1),iopt_arrint_depths)
CALL write_cif_line(cvals(1:1),'IOPT_ARRINT_DEPTHS')
CALL conv_to_chars(cvals(1),iopt_arrint_hreg)
CALL write_cif_line(cvals(1:1),'IOPT_ARRINT_HREG')
CALL conv_to_chars(cvals(1),iopt_arrint_vreg)
CALL write_cif_line(cvals(1:1),'IOPT_ARRINT_VREG')
CALL conv_to_chars(cvals(1),iopt_arrint_3D)
CALL write_cif_line(cvals(1:1),'IOPT_ARRINT_3D')
CALL conv_to_chars(cvals(1),iopt_astro_anal)
CALL write_cif_line(cvals(1:1),'IOPT_ASTRO_ANAL')
CALL conv_to_chars(cvals(1),iopt_astro_pars)
CALL write_cif_line(cvals(1:1),'IOPT_ASTRO_PARS')
CALL conv_to_chars(cvals(1),iopt_astro_tide)
CALL write_cif_line(cvals(1:1),'IOPT_ASTRO_TIDE')
CALL conv_to_chars(cvals(1),iopt_biolgy)
CALL write_cif_line(cvals(1:1),'IOPT_BIOLGY')
CALL conv_to_chars(cvals(1),iopt_bstres_drag)
CALL write_cif_line(cvals(1:1),'IOPT_BSTRES_DRAG')
CALL conv_to_chars(cvals(1),iopt_bstres_form)
CALL write_cif_line(cvals(1:1),'IOPT_BSTRES_FORM')
CALL conv_to_chars(cvals(1),iopt_bstres_nodim)
CALL write_cif_line(cvals(1:1),'IOPT_BSTRES_NODIM')
CALL conv_to_chars(cvals(1),iopt_bstres_waves_bfric)
CALL write_cif_line(cvals(1:1),'IOPT_BSTRES_WAVES_BFRIC')
CALL conv_to_chars(cvals(1),iopt_CDF_abort)
CALL write_cif_line(cvals(1:1),'IOPT_CDF_ABORT')
CALL conv_to_chars(cvals(1),iopt_CDF_fill)
CALL write_cif_line(cvals(1:1),'IOPT_CDF_FILL')
CALL conv_to_chars(cvals(1),iopt_CDF_format)
CALL write_cif_line(cvals(1:1),'IOPT_CDF_FORMAT')
CALL conv_to_chars(cvals(1),iopt_CDF_shared)
CALL write_cif_line(cvals(1:1),'IOPT_CDF_SHARED')
CALL conv_to_chars(cvals(1),iopt_CDF_sync)
CALL write_cif_line(cvals(1:1),'IOPT_CDF_SYNC')
CALL conv_to_chars(cvals(1),iopt_CDF_tlim)
CALL write_cif_line(cvals(1:1),'IOPT_CDF_TLIM')
CALL conv_to_chars(cvals(1),iopt_cor_impl)
CALL write_cif_line(cvals(1:1),'IOPT_COR_IMPL')
CALL conv_to_chars(cvals(1),iopt_curr)
CALL write_cif_line(cvals(1:1),'IOPT_CURR')
CALL conv_to_chars(cvals(1),iopt_curr_wfall)
CALL write_cif_line(cvals(1:1),'IOPT_CURR_WFALL')
CALL conv_to_chars(cvals(1),iopt_dar)
CALL write_cif_line(cvals(1:1),'IOPT_DAR')
CALL conv_to_chars(cvals(1),iopt_dens)
CALL write_cif_line(cvals(1:1),'IOPT_DENS')
CALL conv_to_chars(cvals(1),iopt_dens_convect)
CALL write_cif_line(cvals(1:1),'IOPT_DENS_CONVECT')
CALL conv_to_chars(cvals(1),iopt_dens_grad)
CALL write_cif_line(cvals(1:1),'IOPT_DENS_GRAD')
CALL conv_to_chars(cvals(1),iopt_dischr)
CALL write_cif_line(cvals(1:1),'IOPT_DISCHR')
CALL conv_to_chars(cvals(1),iopt_dischr_land)
CALL write_cif_line(cvals(1:1),'IOPT_DISCHR_LAND')
CALL conv_to_chars(cvals(1),iopt_drycel)
CALL write_cif_line(cvals(1:1),'IOPT_DRYCEL')
CALL conv_to_chars(cvals(1),iopt_fld)
CALL write_cif_line(cvals(1:1),'IOPT_FLD')
CALL conv_to_chars(cvals(1),iopt_fld_alpha)
CALL write_cif_line(cvals(1:1),'IOPT_FLD_ALPHA')
CALL conv_to_chars(cvals(1),iopt_grid_htype)
CALL write_cif_line(cvals(1:1),'IOPT_GRID_HTYPE')
CALL conv_to_chars(cvals(1),iopt_grid_nodim)
CALL write_cif_line(cvals(1:1),'IOPT_GRID_NODIM')
CALL conv_to_chars(cvals(1),iopt_grid_sph)
CALL write_cif_line(cvals(1:1),'IOPT_GRID_SPH')
CALL conv_to_chars(cvals(1),iopt_grid_vtype)
CALL write_cif_line(cvals(1:1),'IOPT_GRID_VTYPE')
CALL conv_to_chars(cvals(1),iopt_grid_vtype_transf)
CALL write_cif_line(cvals(1:1),'IOPT_GRID_VTYPE_TRANSF')
CALL conv_to_chars(cvals(1),iopt_hdif_coef)
CALL write_cif_line(cvals(1:1),'IOPT_HDIF_COEF')
CALL conv_to_chars(cvals(1),iopt_hdif_lim)
CALL write_cif_line(cvals(1:1),'IOPT_HDIF_LIM')
CALL conv_to_chars(cvals(1),iopt_hdif_scal)
CALL write_cif_line(cvals(1:1),'IOPT_HDIF_SCAL')
CALL conv_to_chars(cvals(1),iopt_hdif_turb)
CALL write_cif_line(cvals(1:1),'IOPT_HDIF_TURB')
CALL conv_to_chars(cvals(1),iopt_hdif_2D)
CALL write_cif_line(cvals(1:1),'IOPT_HDIF_2D')
CALL conv_to_chars(cvals(1),iopt_hdif_3D)
CALL write_cif_line(cvals(1:1),'IOPT_HDIF_3D')
CALL conv_to_chars(cvals(1),iopt_hydro_impl)
CALL write_cif_line(cvals(1:1),'IOPT_HYDRO_IMPL')
CALL conv_to_chars(cvals(1),iopt_kinvisc)
CALL write_cif_line(cvals(1:1),'IOPT_KINVISC')
CALL conv_to_chars(cvals(1),iopt_meteo)
CALL write_cif_line(cvals(1:1),'IOPT_METEO')
CALL conv_to_chars(cvals(1),iopt_meteo_data)
CALL write_cif_line(cvals(1:1),'IOPT_METEO_DATA')
CALL conv_to_chars(cvals(1),iopt_meteo_heat)
CALL write_cif_line(cvals(1:1),'IOPT_METEO_HEAT')
CALL conv_to_chars(cvals(1),iopt_meteo_precip)
CALL write_cif_line(cvals(1:1),'IOPT_METEO_PRECIP')
CALL conv_to_chars(cvals(1),iopt_meteo_pres)
CALL write_cif_line(cvals(1:1),'IOPT_METEO_PRES')
CALL conv_to_chars(cvals(1),iopt_meteo_stres)
CALL write_cif_line(cvals(1:1),'IOPT_METEO_STRES')
CALL conv_to_chars(cvals(1),iopt_mg_cycle)
CALL write_cif_line(cvals(1:1),'IOPT_MG_CYCLE')
CALL conv_to_chars(cvals(1),iopt_mg_prolong)
CALL write_cif_line(cvals(1:1),'IOPT_MG_PROLONG')
CALL conv_to_chars(cvals(1),iopt_mg_smoother)
CALL write_cif_line(cvals(1:1),'IOPT_MG_SMOOTHER')
CALL conv_to_chars(cvals(1),iopt_morph)
CALL write_cif_line(cvals(1:1),'IOPT_MORPH')
CALL conv_to_chars(cvals(1),iopt_MPI_abort)
CALL write_cif_line(cvals(1:1),'IOPT_MPI_ABORT')
CALL conv_to_chars(cvals(1),iopt_MPI_comm_all)
CALL write_cif_line(cvals(1:1),'IOPT_MPI_COMM_ALL')
CALL conv_to_chars(cvals(1),iopt_MPI_comm_coll)
CALL write_cif_line(cvals(1:1),'IOPT_MPI_COMM_COLL')
CALL conv_to_chars(cvals(1),iopt_MPI_comm_exch)
CALL write_cif_line(cvals(1:1),'IOPT_MPI_COMM_EXCH')
CALL conv_to_chars(cvals(1),iopt_MPI_comm_full)
CALL write_cif_line(cvals(1:1),'IOPT_MPI_COMM_FULL')
CALL conv_to_chars(cvals(1),iopt_MPI_comm_gath)
CALL write_cif_line(cvals(1:1),'IOPT_MPI_COMM_GATH')
CALL conv_to_chars(cvals(1),iopt_MPI_comm_scat)
CALL write_cif_line(cvals(1:1),'IOPT_MPI_COMM_SCAT')
CALL conv_to_chars(cvals(1),iopt_MPI_partit)
CALL write_cif_line(cvals(1:1),'IOPT_MPI_PARTIT')
CALL conv_to_chars(cvals(1),iopt_MPI_sync)
CALL write_cif_line(cvals(1:1),'IOPT_MPI_SYNC')
CALL conv_to_chars(cvals(1),iopt_nests)
CALL write_cif_line(cvals(1:1),'IOPT_NESTS')
CALL conv_to_chars(cvals(1),iopt_obc_advflux)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_ADVFLUX')
CALL conv_to_chars(cvals(1),iopt_obc_advrlx)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_ADVRLX')
CALL conv_to_chars(cvals(1),iopt_obc_bio)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_BIO')
CALL conv_to_chars(cvals(1),iopt_obc_invbar)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_INVBAR')
CALL conv_to_chars(cvals(1),iopt_obc_relax)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_RELAX')
CALL conv_to_chars(cvals(1),iopt_obc_sal)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_SAL')
CALL conv_to_chars(cvals(1),iopt_obc_sed)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_SED')
CALL conv_to_chars(cvals(1),iopt_obc_temp)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_TEMP')
CALL conv_to_chars(cvals(1),iopt_obc_th)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_TH')
CALL conv_to_chars(cvals(1),iopt_obc_2D)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_2D')
CALL conv_to_chars(cvals(1),iopt_obc_2D_tang)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_2D_TANG')
CALL conv_to_chars(cvals(1),iopt_obc_3D)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_3D')
CALL conv_to_chars(cvals(1),iopt_obc_3D_tang)
CALL write_cif_line(cvals(1:1),'IOPT_OBC_3D_TANG')
CALL conv_to_chars(cvals(1),iopt_out_anal)
CALL write_cif_line(cvals(1:1),'IOPT_OUT_ANAL')
CALL conv_to_chars(cvals(1),iopt_out_avrgd)
CALL write_cif_line(cvals(1:1),'IOPT_OUT_AVRGD')
CALL conv_to_chars(cvals(1),iopt_out_tsers)
CALL write_cif_line(cvals(1:1),'IOPT_OUT_TSERS')
CALL conv_to_chars(cvals(1),iopt_part_write)
CALL write_cif_line(cvals(1:1),'IOPT_PART_WRITE')
CALL conv_to_chars(cvals(1),iopt_rng_seed)
CALL write_cif_line(cvals(1:1),'IOPT_RNG_SEED')
CALL conv_to_chars(cvals(1),iopt_sal)
CALL write_cif_line(cvals(1:1),'IOPT_SAL')
CALL conv_to_chars(cvals(1),iopt_scal_depos)
CALL write_cif_line(cvals(1:1),'IOPT_SCAL_DEPOS')
CALL conv_to_chars(cvals(1),iopt_sed)
CALL write_cif_line(cvals(1:1),'IOPT_SED')
CALL conv_to_chars(cvals(1),iopt_sflux_pars)
CALL write_cif_line(cvals(1:1),'IOPT_SFLUX_PARS')
CALL conv_to_chars(cvals(1),iopt_sflux_precip)
CALL write_cif_line(cvals(1:1),'IOPT_SFLUX_PRECIP')
CALL conv_to_chars(cvals(1),iopt_sflux_qlong)
CALL write_cif_line(cvals(1:1),'IOPT_SFLUX_QLONG')
CALL conv_to_chars(cvals(1),iopt_sflux_qshort)
CALL write_cif_line(cvals(1:1),'IOPT_SFLUX_QSHORT')
CALL conv_to_chars(cvals(1),iopt_sflux_strat)
CALL write_cif_line(cvals(1:1),'IOPT_SFLUX_STRAT')
CALL conv_to_chars(cvals(1),iopt_sur_1D)
CALL write_cif_line(cvals(1:1),'IOPT_SUR_1D')
CALL conv_to_chars(cvals(1),iopt_temp)
CALL write_cif_line(cvals(1:1),'IOPT_TEMP')
CALL conv_to_chars(cvals(1),iopt_temp_optic)
CALL write_cif_line(cvals(1:1),'IOPT_TEMP_OPTIC')
CALL conv_to_chars(cvals(1),iopt_temp_sbc)
CALL write_cif_line(cvals(1:1),'IOPT_TEMP_SBC')
CALL conv_to_chars(cvals(1),iopt_thndam)
CALL write_cif_line(cvals(1:1),'IOPT_THNDAM')
CALL conv_to_chars(cvals(1),iopt_tidal_accel)
CALL write_cif_line(cvals(1:1),'IOPT_TIDAL_ACCEL')
CALL conv_to_chars(cvals(1),iopt_transp_full)
CALL write_cif_line(cvals(1:1),'IOPT_TRANSP_FULL')
CALL conv_to_chars(cvals(1),iopt_turb_alg)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_ALG')
CALL conv_to_chars(cvals(1),iopt_turb_dis_bbc)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_DIS_BBC')
CALL conv_to_chars(cvals(1),iopt_turb_dis_sbc)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_DIS_SBC')
CALL conv_to_chars(cvals(1),iopt_turb_iwlim)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_IWLIM')
CALL conv_to_chars(cvals(1),iopt_turb_kinvisc)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_KINVISC')
CALL conv_to_chars(cvals(1),iopt_turb_lmix)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_LMIX')
CALL conv_to_chars(cvals(1),iopt_turb_ntrans)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_NTRANS')
CALL conv_to_chars(cvals(1),iopt_turb_param)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_PARAM')
CALL conv_to_chars(cvals(1),iopt_turb_stab_form)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_STAB_FORM')
CALL conv_to_chars(cvals(1),iopt_turb_stab_lev)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_STAB_LEV')
CALL conv_to_chars(cvals(1),iopt_turb_stab_mod)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_STAB_MOD')
CALL conv_to_chars(cvals(1),iopt_turb_stab_tke)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_STAB_TKE')
CALL conv_to_chars(cvals(1),iopt_turb_tke_bbc)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_TKE_BBC')
CALL conv_to_chars(cvals(1),iopt_turb_tke_sbc)
CALL write_cif_line(cvals(1:1),'IOPT_TURB_TKE_SBC')
CALL conv_to_chars(cvals(1),iopt_vadv_impl)
CALL write_cif_line(cvals(1:1),'IOPT_VADV_IMPL')
CALL conv_to_chars(cvals(1),iopt_vdif_coef)
CALL write_cif_line(cvals(1:1),'IOPT_VDIF_COEF')
CALL conv_to_chars(cvals(1),iopt_vdif_impl)
CALL write_cif_line(cvals(1:1),'IOPT_VDIF_IMPL')
CALL conv_to_chars(cvals(1),iopt_vdif_rot)
CALL write_cif_line(cvals(1:1),'IOPT_VDIF_ROT')
CALL conv_to_chars(cvals(1),iopt_waves)
CALL write_cif_line(cvals(1:1),'IOPT_WAVES')
CALL conv_to_chars(cvals(1),iopt_waves_couple)
CALL write_cif_line(cvals(1:1),'IOPT_WAVES_COUPLE')
CALL conv_to_chars(cvals(1),iopt_waves_curr)
CALL write_cif_line(cvals(1:1),'IOPT_WAVES_CURR')
CALL conv_to_chars(cvals(1),iopt_waves_dissip)
CALL write_cif_line(cvals(1:1),'IOPT_WAVES_DISSIP')
CALL conv_to_chars(cvals(1),iopt_waves_extrapol)
CALL write_cif_line(cvals(1:1),'IOPT_WAVES_EXTRAPOL')
CALL conv_to_chars(cvals(1),iopt_waves_form)
CALL write_cif_line(cvals(1:1),'IOPT_WAVES_FORM')
CALL conv_to_chars(cvals(1),iopt_waves_pres)
CALL write_cif_line(cvals(1:1),'IOPT_WAVES_PRES')
CALL conv_to_chars(cvals(1),iopt_weibar)
CALL write_cif_line(cvals(1:1),'IOPT_WEIBAR')

!
!3. Model grid
!-------------
!

CALL conv_to_chars(cvals(1),nc)
CALL write_cif_line(cvals(1:1),'NC')
CALL conv_to_chars(cvals(1),nr)
CALL write_cif_line(cvals(1:1),'NR')
CALL conv_to_chars(cvals(1),nz)
CALL write_cif_line(cvals(1:1),'NZ')
CALL conv_to_chars(cvals(1),nosbu)
CALL write_cif_line(cvals(1:1),'NOSBU')
CALL conv_to_chars(cvals(1),nosbv)
CALL write_cif_line(cvals(1:1),'NOSBV')
CALL conv_to_chars(cvals(1),nrvbu)
CALL write_cif_line(cvals(1:1),'NRVBU')
CALL conv_to_chars(cvals(1),nrvbv)
CALL write_cif_line(cvals(1:1),'NRVBV')

!
!4. Dynamic masks
!----------------
!

CALL conv_to_chars(cvals(1:nofldmasks),fld_mask(1:nofldmasks))
CALL write_cif_line(cvals(1:nofldmasks),'FLD_MASK')

!
!5. Number of nested sub-grids and relaxation zones
!--------------------------------------------------
!

CALL conv_to_chars(cvals(1),nonestsets)
CALL write_cif_line(cvals(1:1),'NONESTSETS')
CALL conv_to_chars(cvals(1),norlxzones)
CALL write_cif_line(cvals(1:1),'NORLXZONES')

!
!6. Parameters for parallelisation
!---------------------------------
!

CALL conv_to_chars(cvals(1),nprocsx)
CALL write_cif_line(cvals(1:1),'NPROCSX')
CALL conv_to_chars(cvals(1),nprocsy)
CALL write_cif_line(cvals(1:1),'NPROCSY')

!
!7. Date/time parameters
!-----------------------
!

cvals(1) = CStartDateTime
CALL write_cif_line(cvals(1:1),'CSTARTDATETIME')
cvals(1) = CEndDateTime
CALL write_cif_line(cvals(1:1),'CENDDATETIME')
CALL conv_to_chars(cvals(1),delt2d)
CALL write_cif_line(cvals(1:1),'DELT2D')
CALL conv_to_chars(cvals(1),ic3d)
CALL write_cif_line(cvals(1:1),'IC3D')
CALL conv_to_chars(cvals(1),iccvt)
CALL write_cif_line(cvals(1:1),'ICCVT')
CALL conv_to_chars(cvals(1),icnodal)
CALL write_cif_line(cvals(1:1),'ICNODAL')
CALL conv_to_chars(cvals(1),time_zone)
CALL write_cif_line(cvals(1:1),'TIME_ZONE')
CALL conv_to_chars(cvals(1),ntobcrlx)
CALL write_cif_line(cvals(1:1),'NTOBCRLX')
CALL conv_to_chars(cvals(1:nz),return_time(1:nz))
CALL write_cif_line(cvals(1:nz),'RETURN_TIME')

!
!8. Model constants
!------------------
!
!8.1 Physical parameters
!-----------------------
!

CALL conv_to_chars(cvals(1),atmpres_ref)
CALL write_cif_line(cvals(1:1),'ATMPRES_REF')
CALL conv_to_chars(cvals(1),bdragcoef_cst)
CALL write_cif_line(cvals(1:1),'BDRAGCOEF_CST')
CALL conv_to_chars(cvals(1),bdraglin)
CALL write_cif_line(cvals(1:1),'BDRAGLIN')
CALL conv_to_chars(cvals(1),b_SH)
CALL write_cif_line(cvals(1:1),'B_SH')
CALL conv_to_chars(cvals(1),ccharno)
CALL write_cif_line(cvals(1:1),'CCHARNO')
CALL conv_to_chars(cvals(1:4),cdspars)
CALL write_cif_line(cvals(1:4),'CDSPARS')
CALL conv_to_chars(cvals(1),ces_scst)
CALL write_cif_line(cvals(1:1),'CES_SCST')
CALL conv_to_chars(cvals(1),ces_ucst)
CALL write_cif_line(cvals(1:1),'CES_UCST')
CALL conv_to_chars(cvals(1),chs_scst)
CALL write_cif_line(cvals(1:1),'CHS_SCST')
CALL conv_to_chars(cvals(1),chs_ucst)
CALL write_cif_line(cvals(1:1),'CHS_UCST')
CALL conv_to_chars(cvals(1),cgravratio)
CALL write_cif_line(cvals(1:1),'CGRAVRATIO')
CALL conv_to_chars(cvals(1),ckar)
CALL write_cif_line(cvals(1:1),'CKAR')
CALL conv_to_chars(cvals(1),dcrit_fld)
CALL write_cif_line(cvals(1:1),'DCRIT_FLD')
CALL conv_to_chars(cvals(1),depmean_cst)
CALL write_cif_line(cvals(1:1),'DEPMEAN_CST')
CALL conv_to_chars(cvals(1),depmean_flag)
CALL write_cif_line(cvals(1:1),'DEPMEAN_FLAG')
CALL conv_to_chars(cvals(1),distrlx_obc)
CALL write_cif_line(cvals(1:1),'DISTRLX_OBC')
CALL conv_to_chars(cvals(1),dlat_ref)
CALL write_cif_line(cvals(1:1),'DLAT_REF')
CALL conv_to_chars(cvals(1),dlon_ref)
CALL write_cif_line(cvals(1:1),'DLON_REF')
CALL conv_to_chars(cvals(1),dlon_ref_anal)
CALL write_cif_line(cvals(1:1),'DLON_REF_ANAL')
CALL conv_to_chars(cvals(1),dlon_ref_obc)
CALL write_cif_line(cvals(1:1),'DLON_REF_OBC')
CALL conv_to_chars(cvals(1),dl_BB)
CALL write_cif_line(cvals(1:1),'DL_BB')
CALL conv_to_chars(cvals(1),dmin_fld)
CALL write_cif_line(cvals(1:1),'DMIN_FLD')
CALL conv_to_chars(cvals(1),dthd_fld)
CALL write_cif_line(cvals(1:1),'DTHD_FLD')
CALL conv_to_chars(cvals(1),du_BB)
CALL write_cif_line(cvals(1:1),'DU_BB')
CALL conv_to_chars(cvals(1),dzetaresid_conv)
CALL write_cif_line(cvals(1:1),'DZETARESID_CONV')
CALL conv_to_chars(cvals(1),gacc_ref)
CALL write_cif_line(cvals(1:1),'GACC_REF')
CALL conv_to_chars(cvals(1),hdifmom_cst)
CALL write_cif_line(cvals(1:1),'HDIFMOM_CST')
CALL conv_to_chars(cvals(1),hdifscal_cst)
CALL write_cif_line(cvals(1:1),'HDIFSCAL_CST')
CALL conv_to_chars(cvals(1),kinvisc_cst)
CALL write_cif_line(cvals(1:1),'KINVISC_CST')
CALL conv_to_chars(cvals(1),maxitsimp)
CALL write_cif_line(cvals(1:1),'MAXITSIMP')
CALL conv_to_chars(cvals(1),mg_tol)
CALL write_cif_line(cvals(1:1),'MG_TOL')
CALL conv_to_chars(cvals(1),ndwexp)
CALL write_cif_line(cvals(1:1),'NDWEXP')
CALL conv_to_chars(cvals(1),nomglevels)
CALL write_cif_line(cvals(1:1),'NOMGLEVELS')
CALL conv_to_chars(cvals(1),nomgiterations)
CALL write_cif_line(cvals(1:1),'NOMGITERATIONS')
CALL conv_to_chars(cvals(1),nopresweeps)
CALL write_cif_line(cvals(1:1),'NOPRESWEEPS')
CALL conv_to_chars(cvals(1),nopostsweeps)
CALL write_cif_line(cvals(1:1),'NOPOSTSWEEPS')
CALL conv_to_chars(cvals(1),nosmoothsteps)
CALL write_cif_line(cvals(1:1),'NOSMOOTHSTEPS')
CALL conv_to_chars(cvals(1),optattcoef1_cst)
CALL write_cif_line(cvals(1:1),'OPTATTCOEF1_CST')
CALL conv_to_chars(cvals(1),optattcoef2_cst)
CALL write_cif_line(cvals(1:1),'OPTATTCOEF2_CST')
CALL conv_to_chars(cvals(1),opt_frac)
CALL write_cif_line(cvals(1:1),'OPT_FRAC')
CALL conv_to_chars(cvals(1),rho_air)
CALL write_cif_line(cvals(1:1),'RHO_AIR')
CALL conv_to_chars(cvals(1),sal_ref)
CALL write_cif_line(cvals(1:1),'SAL_REF')
CALL conv_to_chars(cvals(1),sigstar_DJ)
CALL write_cif_line(cvals(1:1),'SIGSTAR_DJ')
CALL conv_to_chars(cvals(1),sig0_DJ)
CALL write_cif_line(cvals(1:1),'SIG0_DJ')
CALL conv_to_chars(cvals(1),skewdiff_cst)
CALL write_cif_line(cvals(1:1),'SKEWDIFF_CST')
CALL conv_to_chars(cvals(1),rho_crit_iso)
CALL write_cif_line(cvals(1:1),'RHO_CRIT_ISO')
CALL conv_to_chars(cvals(1),slopemax_iso)
CALL write_cif_line(cvals(1:1),'slopemax_iso')
CALL conv_to_chars(cvals(1),smag_coef_mom)
CALL write_cif_line(cvals(1:1),'SMAG_COEF_MOM')
CALL conv_to_chars(cvals(1),smag_coef_scal)
CALL write_cif_line(cvals(1:1),'SMAG_COEF_SCAL')
CALL conv_to_chars(cvals(1),specheat)
CALL write_cif_line(cvals(1:1),'SPECHEAT')
CALL conv_to_chars(cvals(1),sst_ref)
CALL write_cif_line(cvals(1:1),'SST_REF')
CALL conv_to_chars(cvals(1),temp_min)
CALL write_cif_line(cvals(1:1),'TEMP_MIN')
CALL conv_to_chars(cvals(1),temp_ref)
CALL write_cif_line(cvals(1:1),'TEMP_REF')
CALL conv_to_chars(cvals(1),theta_cor)
CALL write_cif_line(cvals(1:1),'THETA_COR')
CALL conv_to_chars(cvals(1),theta_SH)
CALL write_cif_line(cvals(1:1),'THETA_SH')
CALL conv_to_chars(cvals(1),theta_sur)
CALL write_cif_line(cvals(1:1),'THETA_SUR')
CALL conv_to_chars(cvals(1),theta_vadv)
CALL write_cif_line(cvals(1:1),'THETA_VADV')
CALL conv_to_chars(cvals(1),theta_vdif)
CALL write_cif_line(cvals(1:1),'THETA_VDIF')
CALL conv_to_chars(cvals(1),ur_mg)
CALL write_cif_line(cvals(1:1),'UR_MG')
CALL conv_to_chars(cvals(1),ur_smooth)
CALL write_cif_line(cvals(1:1),'UR_SMOOTH')
CALL conv_to_chars(cvals(1),vdifmom_cst)
CALL write_cif_line(cvals(1:1),'VDIFMOM_CST')
CALL conv_to_chars(cvals(1),vdifscal_cst)
CALL write_cif_line(cvals(1:1),'VDIFSCAL_CST')
CALL conv_to_chars(cvals(1),wave_penetration_bed)
CALL write_cif_line(cvals(1:1),'WAVE_PENETRATION_BED')
CALL conv_to_chars(cvals(1),wave_penetration_surf)
CALL write_cif_line(cvals(1:1),'WAVE_PENETRATION_SURF')
CALL conv_to_chars(cvals(1),wavethick_cst)
CALL write_cif_line(cvals(1:1),'WAVETHICK_CST')
CALL conv_to_chars(cvals(1),zbtoz0lim)
CALL write_cif_line(cvals(1:1),'ZBTOZ0LIM')
CALL conv_to_chars(cvals(1),zref_temp)
CALL write_cif_line(cvals(1:1),'ZREF_TEMP')
CALL conv_to_chars(cvals(1),zref_wind)
CALL write_cif_line(cvals(1:1),'ZREF_WIND')
CALL conv_to_chars(cvals(1),zrough_cst)
CALL write_cif_line(cvals(1:1),'ZROUGH_CST')

!
!8.2 Tidal constituents
!----------------------
!
!---open boundaries
CALL conv_to_chars(cvals(1),nconobc)
CALL write_cif_line(cvals(1:1),'NCONOBC')
IF (nconobc.GT.0) THEN
   CALL conv_to_chars(cvals(1:nconobc),index_obc(1:nconobc))
   CALL write_cif_line(cvals(1:nconobc),'INDEX_OBC')
ENDIF
CALL conv_to_chars(cvals(1),nqsecobu)
CALL write_cif_line(cvals(1:1),'NQSECOBU')
CALL conv_to_chars(cvals(1),nqsecobv)
CALL write_cif_line(cvals(1:1),'NQSECOBV')

!---tidal force
CALL conv_to_chars(cvals(1),nconastro)
CALL write_cif_line(cvals(1:1),'NCONASTRO')
IF (nconastro.GT.0) THEN
   CALL conv_to_chars(cvals(1:nconastro),index_astro(1:nconastro))
   CALL write_cif_line(cvals(1:nconastro),'INDEX_ASTRO')
ENDIF

!
!8..3 Structure module parameters
!--------------------------------
!
CALL conv_to_chars(cvals(1),numdry)
CALL write_cif_line(cvals(1:1),'NUMDRY')

CALL conv_to_chars(cvals(1),numthinu)
CALL write_cif_line(cvals(1:1),'NUMTHINU')

CALL conv_to_chars(cvals(1),numthinv)
CALL write_cif_line(cvals(1:1),'NUMTHINV')

CALL conv_to_chars(cvals(1),numwbaru)
CALL write_cif_line(cvals(1:1),'NUMWBARU')

CALL conv_to_chars(cvals(1),numwbarv)
CALL write_cif_line(cvals(1:1),'NUMWBARV')

CALL conv_to_chars(cvals(1),wbarrlxu)
CALL write_cif_line(cvals(1:1),'WBARRLXU')

CALL conv_to_chars(cvals(1),wbarrlxv)
CALL write_cif_line(cvals(1:1),'WBARRLXV')

CALL conv_to_chars(cvals(1),numdis)
CALL write_cif_line(cvals(1:1),'NUMDIS')

!
!8..4 Turbulence parameters
!--------------------------
!

CALL conv_to_chars(cvals(1),alpha_black)
CALL write_cif_line(cvals(1:1),'ALPHA_BLACK')
CALL conv_to_chars(cvals(1),alpha_ma)
CALL write_cif_line(cvals(1:1),'ALPHA_MA')
CALL conv_to_chars(cvals(1),alpha_pp)
CALL write_cif_line(cvals(1:1),'ALPHA_PP')
CALL conv_to_chars(cvals(1),beta_ma)
CALL write_cif_line(cvals(1:1),'BETA_MA')
CALL conv_to_chars(cvals(1),beta_xing)
CALL write_cif_line(cvals(1:1),'BETA_XING')
CALL conv_to_chars(cvals(1),cnu_ad)
CALL write_cif_line(cvals(1:1),'CNU_AD')
CALL conv_to_chars(cvals(1),c_sk)
CALL write_cif_line(cvals(1:1),'C_SK')
CALL conv_to_chars(cvals(1),c1_eps)
CALL write_cif_line(cvals(1:1),'C1_EPS')
CALL conv_to_chars(cvals(1),c2_eps)
CALL write_cif_line(cvals(1:1),'C2_EPS')
CALL conv_to_chars(cvals(1),c31_eps)
CALL write_cif_line(cvals(1:1),'C31_EPS')
CALL conv_to_chars(cvals(1),c32_eps)
CALL write_cif_line(cvals(1:1),'C32_EPS')
CALL conv_to_chars(cvals(1),delta1_ad)
CALL write_cif_line(cvals(1:1),'DELTA1_AD')
CALL conv_to_chars(cvals(1),delta2_ad)
CALL write_cif_line(cvals(1:1),'DELTA2_AD')
CALL conv_to_chars(cvals(1),dissipmin)
CALL write_cif_line(cvals(1:1),'DISSIPMIN')
CALL conv_to_chars(cvals(1),expmom_ma)
CALL write_cif_line(cvals(1:1),'EXPMOM_MA')
CALL conv_to_chars(cvals(1),expmom_pp)
CALL write_cif_line(cvals(1:1),'EXPMOM_PP')
CALL conv_to_chars(cvals(1),expscal_ma)
CALL write_cif_line(cvals(1:1),'EXPSCAL_MA')
CALL conv_to_chars(cvals(1),e1_my)
CALL write_cif_line(cvals(1:1),'E1_MY')
CALL conv_to_chars(cvals(1),e2_my)
CALL write_cif_line(cvals(1:1),'E2_MY')
CALL conv_to_chars(cvals(1),e3_my)
CALL write_cif_line(cvals(1:1),'E3_MY')
CALL conv_to_chars(cvals(1),k1_ad)
CALL write_cif_line(cvals(1:1),'K1_AD')
CALL conv_to_chars(cvals(1),k2_ad)
CALL write_cif_line(cvals(1:1),'K2_AD')
CALL conv_to_chars(cvals(1),lambda_ad)
CALL write_cif_line(cvals(1:1),'LAMBDA_AD')
CALL conv_to_chars(cvals(1),omega1_ad)
CALL write_cif_line(cvals(1:1),'OMEGA1_AD')
CALL conv_to_chars(cvals(1),riccrit_iw)
CALL write_cif_line(cvals(1:1),'RICCRIT_IW')
CALL conv_to_chars(cvals(1),r1_ad)
CALL write_cif_line(cvals(1:1),'R1_AD')
CALL conv_to_chars(cvals(1),r2_ad)
CALL write_cif_line(cvals(1:1),'R2_AD')
CALL conv_to_chars(cvals(1),sigma_k)
CALL write_cif_line(cvals(1:1),'SIGMA_K')
CALL conv_to_chars(cvals(1),skeps)
CALL write_cif_line(cvals(1:1),'SKEPS')
CALL conv_to_chars(cvals(1),sq_my)
CALL write_cif_line(cvals(1:1),'SQ_MY')
CALL conv_to_chars(cvals(1),tkelim)
CALL write_cif_line(cvals(1:1),'TKELIM')
CALL conv_to_chars(cvals(1),tkemin)
CALL write_cif_line(cvals(1:1),'TKEMIN')
CALL conv_to_chars(cvals(1),vbmom_pp)
CALL write_cif_line(cvals(1:1),'VBMOM_PP')
CALL conv_to_chars(cvals(1),vbscal_pp)
CALL write_cif_line(cvals(1:1),'VBSCAL_PP')
CALL conv_to_chars(cvals(1),vdifmom_iw)
CALL write_cif_line(cvals(1:1),'VDIFMOM_IW')
CALL conv_to_chars(cvals(1),vdifscal_iw)
CALL write_cif_line(cvals(1:1),'VDIFSCAL_IW')
CALL conv_to_chars(cvals(1),vdifshear_iw)
CALL write_cif_line(cvals(1:1),'VDIFSHEAR_IW')
CALL conv_to_chars(cvals(1),vmaxmom_ma)
CALL write_cif_line(cvals(1:1),'VMAXMOM_MA')
CALL conv_to_chars(cvals(1),vmaxscal_ma)
CALL write_cif_line(cvals(1:1),'VMAXSCAL_MA')
CALL conv_to_chars(cvals(1),vmax_pp)
CALL write_cif_line(cvals(1:1),'VMAX_PP')
CALL conv_to_chars(cvals(1),v0dif_ma)
CALL write_cif_line(cvals(1:1),'V0DIF_MA')
CALL conv_to_chars(cvals(1),v0dif_pp)
CALL write_cif_line(cvals(1:1),'V0DIF_PP')
CALL conv_to_chars(cvals(1),wfltke)
CALL write_cif_line(cvals(1:1),'WFLTKE')
CALL conv_to_chars(cvals(1),zlmixmin)
CALL write_cif_line(cvals(1:1),'ZLMIXMIN')
CALL conv_to_chars(cvals(1),zrough_bot)
CALL write_cif_line(cvals(1:1),'ZROUGH_BOT')
CALL conv_to_chars(cvals(1),zrough_sur)
CALL write_cif_line(cvals(1:1),'ZROUGH_SUR')

!
!9. I/O specifiers
!-----------------
!
!---restart times
CALL conv_to_chars(cvals(1),norestarts)
CALL write_cif_line(cvals(1:1),'NORESTARTS')
IF (norestarts.GT.0) THEN
   CALL conv_to_chars(cvals(1:norestarts),ntrestart(1:norestarts))
   CALL write_cif_line(cvals(1:norestarts),'NTRESTART')
ENDIF

!---input/output titles
cvals(1) = TRIM(intitle)
CALL write_cif_line(cvals(1:1),'INTITLE')
cvals(1) = TRIM(outtitle)
CALL write_cif_line(cvals(1:1),'OUTTITLE')

!---CF global atrributes
cvals(1) = TRIM(institution_CF)
CALL write_cif_line(cvals(1:1),'INSTITUTION_CF')
cvals(1) = TRIM(comment_CF)
CALL write_cif_line(cvals(1:1),'COMMENT_CF')
cvals(1) = TRIM(references_CF)
CALL write_cif_line(cvals(1:1),'REFERENCES_CF')

!---I/O parameters
CALL conv_to_chars(cvals(1),maxwaitsecs)
CALL write_cif_line(cvals(1:1),'MAXWAITSECS')
CALL conv_to_chars(cvals(1),nowaitsecs)
CALL write_cif_line(cvals(1:1),'NOWAITSECS')
CALL conv_to_chars(cvals(1),nrecunit)
CALL write_cif_line(cvals(1:1),'NRECUNIT')

!---user-output parameters
CALL conv_to_chars(cvals(1),nosetstsr)
CALL write_cif_line(cvals(1:1),'NOSETSTSR')
CALL conv_to_chars(cvals(1),nostatstsr)
CALL write_cif_line(cvals(1:1),'NOSTATSTSR')
CALL conv_to_chars(cvals(1),novarstsr)
CALL write_cif_line(cvals(1:1),'NOVARSTSR')
CALL conv_to_chars(cvals(1),nosetsavr)
CALL write_cif_line(cvals(1:1),'NOSETSAVR')
CALL conv_to_chars(cvals(1),nostatsavr)
CALL write_cif_line(cvals(1:1),'NOSTATSAVR')
CALL conv_to_chars(cvals(1),novarsavr)
CALL write_cif_line(cvals(1:1),'NOVARSAVR')
CALL conv_to_chars(cvals(1),nosetsanal)
CALL write_cif_line(cvals(1:1),'NOSETSANAL')
CALL conv_to_chars(cvals(1),nofreqsanal)
CALL write_cif_line(cvals(1:1),'NOFREQSANAL')
CALL conv_to_chars(cvals(1),nostatsanal)
CALL write_cif_line(cvals(1:1),'NOSTATSANAL')
CALL conv_to_chars(cvals(1),novarsanal)
CALL write_cif_line(cvals(1:1),'NOVARSANAL')

!
!10. Attributes of forcing files
!-------------------------------
!

iotype_1010: DO iotype=1,2
ifil_1010: DO ifil=1,MaxIOFiles
idesc_1010: DO idesc=1,MaxIOTypes
   IF (modfiles(idesc,ifil,iotype)%status.NE.'0') THEN
      cvals(1) = modfiles_desc(idesc)
      CALL conv_to_chars(cvals(2),ifil)
      CALL conv_to_chars(cvals(3),iotype)
      CALL conv_to_chars_infiles(cvals(4:15),modfiles(idesc,ifil,iotype))
      CALL write_cif_line(cvals(1:15),'MODFILES')
   ENDIF
ENDDO idesc_1010
ENDDO ifil_1010
ENDDO iotype_1010

!
!11. Data grid attributes
!------------------------
!

idesc_1110: DO idesc=1,MaxGridTypes
   IF (surfacegrids(idesc,1)%nhtype.GT.0) THEN
      CALL conv_to_chars(cvals(1),idesc)
      CALL conv_to_chars_surfgrd(cvals(2:15),surfacegrids(idesc,1))
      CALL write_cif_line(cvals(1:14),'SURFACEGRIDS')
   ENDIF
ENDDO idesc_1110

WRITE (iunit,'(A)') cifend
ciflinenum = ciflinenum + 1

CALL log_timer_out()


RETURN

END SUBROUTINE write_cif_vars_mod

!========================================================================

SUBROUTINE write_cif_vars_mon
!************************************************************************
!
! *write_cif_vars_mon* Write the CIF block with general/monitoring parameters
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - conv_to_chars, write_cif_line
!
!************************************************************************
!
USE iopars
USE paralpars
USE switches
USE syspars
USE cif_routines, ONLY: conv_to_chars, write_cif_line
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: iunit, l


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_cif_vars_mon'
CALL log_timer_in()

!---write block header
iunit = ciffile%iunit
WRITE (iunit,'(A)') cifend//TRIM(cif_block_names(icif_mon))

!---cold/warm start
CALL conv_to_chars(cvals(1),.FALSE.)
CALL write_cif_line(cvals(1:1),'COLD_START')

!---log files
CALL conv_to_chars(cvals(1:npworld),levprocs_ini(1:npworld))
CALL write_cif_line(cvals(1:npworld),'LEVPROCS_INI')
CALL conv_to_chars(cvals(1:npworld),levprocs_run(1:npworld))
CALL write_cif_line(cvals(1:npworld),'LEVPROCS_RUN')
IF (parallel_set) THEN
   l = LEN_TRIM(inilog_file)
   cvals(1) = inilog_file(1:l-1)
ELSE
   cvals(1) = inilog_file
ENDIF
CALL write_cif_line(cvals(1:1),'INILOG_FILE')
cvals(1) = runlog_file
CALL write_cif_line(cvals(1:1),'RUNLOG_FILE')
CALL conv_to_chars(cvals(1),runlog_count)
CALL write_cif_line(cvals(1:1),'RUNLOG_COUNT')

!---error files
CALL conv_to_chars(cvals(1),maxerrors)
CALL write_cif_line(cvals(1:1),'MAX_ERRORS')
CALL conv_to_chars(cvals(1:npworld),levprocs_err(1:npworld))
CALL write_cif_line(cvals(1:npworld),'LEVPROCS_ERR')
IF (parallel_set) THEN
   l = LEN_TRIM(errlog_file)
   cvals(1) = errlog_file(1:l-1)
ELSE
   cvals(1) = errlog_file
ENDIF
CALL write_cif_line(cvals(1:1),'ERRLOG_FILE')

!---warning file
CALL conv_to_chars(cvals(1),warning)
CALL write_cif_line(cvals(1:1),'WARNING')
cvals(1) = warlog_file
CALL write_cif_line(cvals(1:1),'WARLOG_FILE')

!---monitoring
CALL conv_to_chars(cvals(1),monlog)
CALL write_cif_line(cvals(1:1),'MONLOG')
cvals(1) = monlog_file
CALL write_cif_line(cvals(1:1),'MONLOG_FILE')
CALL conv_to_chars(cvals(1),sedlog)
CALL write_cif_line(cvals(1:1),'SEDLOG')
cvals(1) = sedlog_file
CALL write_cif_line(cvals(1:1),'SEDLOG_FILE')

!---timer report
CALL conv_to_chars(cvals(1),levtimer)
CALL write_cif_line(cvals(1:1),'LEVTIMER')
cvals(1) = timing_file
CALL write_cif_line(cvals(1:1),'TIMING_FILE')
CALL conv_to_chars(cvals(1),timer_format)
CALL write_cif_line(cvals(1:1),'TIMER_FORMAT')

!---number of processes for each model component
CALL conv_to_chars(cvals(1),nprocscoh)
CALL write_cif_line(cvals(1:1),'NPROCSCOH')
CALL conv_to_chars(cvals(1),nprocswav)
CALL write_cif_line(cvals(1:1),'NPROCSWAV')

!---switches for model coupling
CALL conv_to_chars(cvals(1),iopt_part_model)
CALL write_cif_line(cvals(1:1),'IOPT_PART_MODEL')
CALL conv_to_chars(cvals(1),iopt_waves_model)
CALL write_cif_line(cvals(1:1),'IOPT_WAVES_MODEL')

WRITE (iunit,'(A)') cifend
ciflinenum = ciflinenum + 1

CALL log_timer_out()

RETURN

END SUBROUTINE write_cif_vars_mon

!========================================================================

SUBROUTINE write_cif_vars_tsout
!************************************************************************
!
! *write_cif_vars_tsout* Write the CIF block with time series parameters
!
! Author - ANTEA and Patrick Luyten
!
! Version - @(COHERENS)Model_Parameters.f90  V2.11
!
! Description - 
!
! Reference -
!
! Calling program - initialise_model
!
! External calls -
!
! Module calls - conv_to_chars, conv_to_chars_gridpars, conv_to_chars_outfiles,
!                conv_to_chars_outvars, conv_to_chars_statlocs, write_cif_line
!
!************************************************************************
!
USE iopars
USE paralpars  
USE syspars
USE cif_routines, ONLY: conv_to_chars, conv_to_chars_gridpars, &
                      & conv_to_chars_outfiles, conv_to_chars_outvars, &
                      & conv_to_chars_statlocs, write_cif_line
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
CHARACTER (LEN=lencifvar), DIMENSION(MaxCIFVars) :: cvals
INTEGER :: iset, istat, iunit, ivar, lvar, nostats


IF (.NOT.master) RETURN

procname(pglev+1) = 'write_cif_vars_tsout'
CALL log_timer_in()

!---write block header
iunit = ciffile%iunit
WRITE (iunit,'(A)') cifend//TRIM(cif_block_names(icif_tsout))

!---variable attributes
ivar_110: DO ivar=1,novarstsr
   CALL conv_to_chars(cvals(1),ivar)
   CALL conv_to_chars_outvars(cvals(2:12),tsrvars(ivar))
   CALL write_cif_line(cvals(1:12),'TSRVARS')
ENDDO ivar_110

!---variable indices
iset_120: DO iset=1,nosetstsr
   lvar = 1
   CALL conv_to_chars(cvals(1),iset)
   ivar_121: DO ivar=1,novarstsr
      IF (ivarstsr(iset,ivar).GT.0) THEN
         lvar = lvar + 1
         CALL conv_to_chars(cvals(lvar),ivarstsr(iset,ivar))
      ENDIF
   ENDDO ivar_121
   CALL write_cif_line(cvals(1:lvar),'IVARSTSR')
ENDDO iset_120

!---file attributes
iset_130: DO iset=1,nosetstsr
   IF (tsr0d(iset)%defined) THEN
      CALL conv_to_chars(cvals(1),iset)
      CALL conv_to_chars_outfiles(cvals(2:7),tsr0d(iset))
      CALL write_cif_line(cvals(1:7),'TSR0D')
   ENDIF
   IF (tsr2d(iset)%defined) THEN
      CALL conv_to_chars(cvals(1),iset)
      CALL conv_to_chars_outfiles(cvals(2:7),tsr2d(iset))
      CALL write_cif_line(cvals(1:7),'TSR2D')
   ENDIF
   IF (tsr3d(iset)%defined) THEN
      CALL conv_to_chars(cvals(1),iset)
      CALL conv_to_chars_outfiles(cvals(2:7),tsr3d(iset))
      CALL write_cif_line(cvals(1:7),'TSR3D')
   ENDIF
ENDDO iset_130

!---output grid (space/time) attributes
iset_140: DO iset=1,nosetstsr
   CALL conv_to_chars(cvals(1),iset)
   CALL conv_to_chars_gridpars(cvals(2:20),tsrgpars(iset))
   CALL write_cif_line(cvals(1:20),'TSRGPARS')
ENDDO iset_140

!---station attributes
istat_150: DO istat=1,nostatstsr
   CALL conv_to_chars(cvals(1),istat)
   CALL conv_to_chars_statlocs(cvals(2:4),tsrstatlocs(istat))
   CALL write_cif_line(cvals(1:4),'TSRSTATLOCS')
ENDDO istat_150

!---station labels
iset_160: DO iset=1,nosetstsr
   CALL conv_to_chars(cvals(1),iset)
   nostats = tsrgpars(iset)%nostats
   IF (nostats.GT.0) THEN 
      CALL conv_to_chars(cvals(2:nostats+1),lstatstsr(iset,1:nostats))
      CALL write_cif_line(cvals(1:nostats+1),'LSTATSTSR')
   ENDIF
ENDDO iset_160

WRITE (iunit,'(A)') cifend
ciflinenum = ciflinenum + 1

CALL log_timer_out()

RETURN

END SUBROUTINE write_cif_vars_tsout
