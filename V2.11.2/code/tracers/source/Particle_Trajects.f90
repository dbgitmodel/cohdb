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

SUBROUTINE particle_trajects
!************************************************************************
!
! *particle_tsout* Particle trajectory output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Particle_Trajects.f90  V2.11
!
! $Date: 2017-08-25 15:58:53 +0200 (Fri, 25 Aug 2017) $
!
! $Revision: 1046 $
!
! Description - particle trajectories output
!
! Reference -
!
! Calling program - coherens_main, particle_model
!
! External calls - particle_trajects_init
!
! Module calls - close_filepars, Cint_at_p, diff_dates, error_alloc, loop_index,
!                write_time, write_vars
!
!************************************************************************
!
USE depths  
USE iopars
USE paralpars
USE partpars  
USE partswitches
USE partvars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: close_filepars, write_time, write_vars
USE particle_routines, ONLY: Cint_at_p
USE time_routines, ONLY: diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: loop_index

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: aging, drift, outflag
CHARACTER (LEN=lentime) :: refdate
INTEGER :: i, iset, ivar, j, l, label, nocoords, nodim, nopartout, novars, &
         & npcc, ntype, p, pp
INTEGER, DIMENSION(4) :: vecids
REAL :: depmeanp, deptotp, fill_value, x, y
REAL (KIND=kndrlong) :: rtime
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: outdat


IF (.NOT.(iopt_part_model.EQ.1.OR.&
       & (iopt_part_model.EQ.2.AND.modelid.EQ.modelidpart).OR.&
        & iopt_part_model.GT.2).OR.iopt_part_out.EQ.0) RETURN

procname(pglev+1) = 'particle_trajects'
CALL log_timer_in(npcc)

!
!1. Initialisation on first call
!-------------------------------
!    

IF (nt.EQ.0) CALL particle_trajects_init

!
!2. Define/write trajectory data
!-------------------------------
!

iset_200: DO iset=1,nosetspart
   IF (part1d(iset)%defined.AND.loop_index(outppars(iset)%tlims,nt)) THEN

!
!2.1 Initialise parameters
!-------------------------
!
      
      refdate = outppars(iset)%refdate
      nodim = outppars(iset)%nodim
      label = outppars(iset)%label
      ntype = outppars(iset)%ntype
      nopartout = outppars(iset)%nopartout
      aging = outppars(iset)%aging
      fill_value = double_fill
      novars = part1d(iset)%novars
      nocoords = part1d(iset)%nocoords

!
!2.2 Output time
!--------------
!      
!     ---define      
      CALL diff_dates(refdate,CDateTime,ptime_unit,rtime=rtime)

!     ---write
      CALL write_time(rtime,part1d(iset))

!
!2.3 Allocate      
!------------
!

      ALLOCATE (outdat(nopartout,novars),STAT=errstat)
      CALL error_alloc('outdat',2,(/nopartout,novars/),kndrtype)
      outdat = fill_value
      
!
!2.4 Store
!---------
!

      pp = 0
      p_240: DO p=1,nopart
         drift = part(p)%drift_state.EQ.2.OR.part(p)%drift_state.EQ.3
         IF (nodim.EQ.0.OR.(nodim.EQ.2.AND.part(p)%state.EQ.1).OR.&
          & (nodim.EQ.3.AND.part(p)%state.EQ.2)) THEN
            
            SELECT CASE (ntype)
               CASE (0); outflag = .TRUE.
               CASE (1); outflag = part(p)%label.EQ.label
               CASE (2); outflag = part(p)%center
               CASE (3); outflag = part(p)%label.EQ.label.AND.part(p)%center
            END SELECT

            IF (outflag) THEN
               pp = pp + 1
               IF (drift) THEN
                  outdat(pp,1) = part(p)%xpos
                  outdat(pp,2) = part(p)%ypos
               ENDIF
               l = 2
               IF (nodim.NE.2) THEN
                  l = l + 1
                  IF (drift) THEN
                     IF (outppars(iset)%ktype.EQ.0) THEN
                        outdat(pp,l) = part(p)%zpos
                     ELSE
                        i = part(p)%icoord; j = part(p)%jcoord
                        x = part(p)%xcoord; y = part(p)%ycoord
                        depmeanp = Cint_at_p(depmeanatc(i-1:i+1,j-1:j+1),i,j,x,y)
                        IF (outppars(iset)%ktype.EQ.1) THEN
                           outdat(pp,l) = part(p)%zpos + depmeanp
                        ELSE
                           deptotp = Cint_at_p(deptotatc(i-1:i+1,j-1:j+1),i,j,x,y)
                           outdat(pp,l) = deptotp - depmeanp - part(p)%zpos
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
               IF (aging) THEN
                  l = l + 1
                  IF (drift) outdat(pp,l) = part(p)%age
               ENDIF
            ENDIF
            
         ENDIF
      ENDDO p_240

!      
!2.5 Write data
!--------------
!      
      vecids(1:novars) = (/(ivar,ivar=nocoords+1,nocoords+novars)/)
      CALL write_vars(outdat,part1d(iset),0,varatts=outpvars(iset,1:novars),&
                    & vecids=vecids(1:novars))

!
!2.6 Deallocate
!--------------  
!

      DEALLOCATE (outdat)
      
   ENDIF

!
!2.7 Close files after last output
!---------------------------------
!

   IF (part1d(iset)%defined.AND.&
    & (nt.EQ.outppars(iset)%tlims(2).OR.cold_start)) THEN
      CALL close_filepars(part1d(iset))
   ENDIF
   
ENDDO iset_200

IF (modelid.EQ.modelidpart) THEN
   CALL log_timer_out()
ELSE
   CALL log_timer_out(npcc,itm_part)
ENDIF


RETURN

END SUBROUTINE particle_trajects

!========================================================================

SUBROUTINE particle_trajects_init
!************************************************************************
!
! *particle_trajects_init* Initialise parameters for particle trajectory
!                          output
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Trajects.f90  V2.11
!
! Description -
!
! Reference -
!
! Calling program - particle_trajects
!
! External calls - read_cif_params, usrdef_pout_params
!
! Module calls - check_out_filepars, check_out_ppars, cf90_def_dim,
!                cf90_def_var, cf90_enddef, cf90_put_att, cf90_put_att_chars,
!                close_file, conv_to_chars, conv_to_chars_modvars,
!                default_out_files, default_out_ppars, error_alloc,
!                error_alloc_struc, error_file, file_suffix, filepars_init,
!                inquire_var, open_file, open_filepars, outppart_init,
!                reset_out_files_part, reset_out_ppars, varatts_init,
!                warning_reset_arr_struc, write_meta_data_line, write_vars,
!                UVarr_at_C
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
USE partids
USE partpars
USE partswitches
USE partvars
USE switches
USE syspars
USE array_interp, ONLY: UVarr_at_C
USE cf90_routines, ONLY: cf90_def_dim, cf90_def_var, cf90_enddef, &
                       & cf90_put_att, cf90_put_att_chars
USE check_model, ONLY: check_out_filepars
USE check_particles, ONLY: check_out_ppars
USE cif_routines, ONLY: conv_to_chars, conv_to_chars_modvars
USE datatypes_init, ONLY: filepars_init, varatts_init
USE default_model, ONLY: default_out_files
USE default_particles, ONLY: default_out_ppars
USE error_routines, ONLY: error_alloc, error_alloc_struc, error_file, &
                        & warning_reset_arr_struc
USE inout_routines, ONLY: close_file, open_file, open_filepars, &
                        & write_metadata_line, write_vars
USE modvars_routines, ONLY: file_suffix, inquire_var
USE parttypes_init, ONLY: outppart_init
USE reset_particles, ONLY: reset_out_files_part, reset_out_ppars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: aging, info
CHARACTER (LEN=1) :: dataform, floattype, outform
CHARACTER (LEN=3) :: suffix
CHARACTER (LEN=12) :: cset
CHARACTER (LEN=lenunit) :: tunits
CHARACTER (LEN=leniofile) :: infofile
CHARACTER (LEN=lencifline) :: cline
CHARACTER (LEN=lendesc), DIMENSION(12) :: cvals
CHARACTER (LEN=lenname), SAVE, ALLOCATABLE, DIMENSION(:) :: dimnames
INTEGER :: dimname, dimp, dimt, dimx, dimy, ifil, iglb, iodat, ioinfo, iset, &
         & iunit, ivar, l, label, nocoords, nodim, nodims, noglbatts, &
         & nopartout, novars, nrank, nstepout, ntype, nvals, varid
INTEGER, DIMENSION(4) :: ivarid
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: outdims
REAL :: fill_value, deltout
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: depout, xout, yout 
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: coordatts


IF (iopt_part_out.EQ.0) RETURN

procname(pglev+1) = 'particle_trajects_init'
CALL log_timer_in()

!
!1. Parameters for particle output
!---------------------------------
!
!1.1 Allocate
!------------
!
!---variable attributes
ALLOCATE (outpvars(nosetspart,4),STAT=errstat)
CALL error_alloc_struc('outpvars',2,(/nosetspart,4/),'VariableAtts')
CALL varatts_init(outpvars)

!---file attributes
ALLOCATE (part1d(nosetspart),STAT=errstat)
CALL error_alloc_struc('part1d',1,(/nosetspart/),'FileParams')
CALL filepars_init(part1d)

!---output parameters
ALLOCATE (outppars(nosetspart),STAT=errstat)
CALL error_alloc_struc('outppars',1,(/nosetspart/),'OutPartParams')
CALL outppart_init(outppars)

!
!1.2 Defaults
!------------
!
!---file attributes
CALL default_out_files(part1d,'PT')

!---output parameters
CALL default_out_ppars

!
!1.3 Define
!----------
!

IF (ciffile%status.EQ.'R') THEN
   CALL read_cif_params(icif_tspart)
ELSE
   CALL usrdef_pout_params
ENDIF

!
!1.4 Reset
!---------
!
!---file parameters
CALL reset_out_files_part(part1d)

!---output parameters
CALL reset_out_ppars(.TRUE.)

!
!1.5 Check
!---------
!
!---file attributes
CALL check_out_filepars(part1d,'part1d')

!---output parameters
CALL check_out_ppars

!
!1.6 Return if needed
!--------------------
!

IF (iopt_part_model.EQ.2.AND.(modelid.EQ.modelidcoh)) RETURN

!
!1.7 Number and attributes of output variables
!---------------------------------------------
!

iset_170: DO iset=1,nosetspart

!  ---number of output variables   
   novars = 2
   ivarid(1:2) = (/iarr_part_xpos,iarr_part_ypos/)
   IF (outppars(iset)%nodim.EQ.0.OR.&
    & (outppars(iset)%nodim.EQ.3.AND.ANY(part%state.EQ.2))) THEN
      novars = novars + 1
      ivarid(novars) = iarr_part_zpos
   ENDIF
   IF (outppars(iset)%aging) THEN
      novars = novars + 1
      ivarid(novars) = iarr_part_age
   ENDIF
   part1d(iset)%novars = novars

!  ---variable attributes   
   ivar_171: DO ivar=1,novars
      CALL inquire_var(ivarid(ivar),f90_name=outpvars(iset,ivar)%f90_name,&
                     & standard_name=outpvars(iset,ivar)%standard_name,&
                     & long_name=outpvars(iset,ivar)%long_name,&
                     & units=outpvars(iset,ivar)%units)
      outpvars(iset,ivar)%nrank = 1
   ENDDO ivar_171
   
ENDDO iset_170

!
!1.8 Number of output particles
!------------------------------
!

iset_180: DO iset=1,nosetspart

   IF (part1d(iset)%defined) THEN

      nodim = outppars(iset)%nodim
      label = outppars(iset)%label
      ntype = outppars(iset)%ntype
      aging = outppars(iset)%aging

      IF (nodim.EQ.0) THEN
         SELECT CASE (ntype)
            CASE (0); nopartout = nopart
            CASE (1); nopartout = COUNT(part%label.EQ.label)
            CASE (2); nopartout = COUNT(part%center)
            CASE (3); nopartout = COUNT(part%label.EQ.label.AND.part%center)
         END SELECT
      ELSEIF (nodim.EQ.2) THEN
         SELECT CASE (ntype)
            CASE (0); nopartout = COUNT(part%state.EQ.1)
            CASE (1); nopartout = COUNT(part%state.EQ.1.AND.part%label.EQ.label)
            CASE (2); nopartout = COUNT(part%state.EQ.1.AND.part%center)
            CASE (3); nopartout = COUNT(part%state.EQ.1.AND.&
                                      & part%label.EQ.label.AND.part%center)
         END SELECT
      ELSEIF (nodim.EQ.3) THEN
         SELECT CASE (ntype)
            CASE (0); nopartout = COUNT(part%state.EQ.2)
            CASE (1); nopartout = COUNT(part%state.EQ.2.AND.part%label.EQ.label)
            CASE (2); nopartout = COUNT(part%state.EQ.2.AND.part%center)
            CASE (3); nopartout = COUNT(part%state.EQ.2.AND.&
                                      & part%label.EQ.label.AND.part%center)
         END SELECT
      ENDIF

      outppars(iset)%nopartout = nopartout

!     ---write output only if number of particles is greater than zero
      IF (nopartout.EQ.0) THEN
         CALL warning_reset_arr_struc(part1d(iset)%defined,'part1d','%defined',&
                                   & .FALSE.,1,(/iset/))
      ENDIF
      
   ENDIF

ENDDO iset_180

!
!2. Open output files and write metadata
!----------------------------------------
!

iset_200: DO iset=1,nosetspart
   filepars = part1d(iset)
   IF (.NOT.filepars%defined) CYCLE iset_200
   
!
!2.1 Output parameters
!---------------------
!   
!  ---initialise
   info = filepars%info
   novars = filepars%novars
   outform = filepars%form
   floattype = filepars%floattype
   fill_value = double_fill
   nocoords = 4

!  ---global attributes
   glbatts(2)%value = TRIM(filepars%title)

!  ---number of output particles
   nopartout = outppars(iset)%nopartout 

!  ---number of time steps
   nstepout = outppars(iset)%nstepout

!  ---time step
   deltout = outppars(iset)%deltout

!  ---time unit in CF format
   SELECT CASE (ptime_unit)
      CASE (1); tunits = 'seconds since'
      CASE (2); tunits = 'minutes since'
      CASE (3); tunits = 'hours since'
      CASE (4); tunits = 'days since'
   END SELECT
   tunits = TRIM(tunits)//' '//outppars(iset)%refdate(1:19)

!  ---file name
   IF (TRIM(filepars%filename).EQ.'') THEN
      suffix = file_suffix(outform)
      WRITE (cset,'(I12)') iset; cset = ADJUSTL(cset)
      filepars%filename = TRIM(outtitle)//'_'//TRIM(cset)//'.'//'tspar1.'//&
                        & TRIM(suffix)
      IF (TRIM(filepars%pathname).NE.'') THEN
         filepars%filename = TRIM(filepars%pathname)//TRIM(filepars%filename)
      ENDIF
   ENDIF
   IF (info) infofile = TRIM(filepars%filename)//'.inf'

!
!2.2 Open files
!--------------
!

   CALL open_filepars(filepars)
   iodat = filepars%iunit
   IF (info) CALL open_file(ioinfo,infofile,'OUT','A')

!
!2.3 Define netCDF dimensions
!----------------------------
!   

   IF (outform.EQ.'N') THEN
      IF (iopt_CDF_tlim.EQ.1) THEN
         CALL cf90_def_dim(iodat,'time',nstepout,dimt)
      ELSE
         CALL cf90_def_dim(iodat,'time',unlimited_NF90,dimt)
      ENDIF
      IF (iopt_grid_htype.LE.2) THEN
         IF (iopt_grid_sph.EQ.0) THEN
            CALL cf90_def_dim(iodat,'x',nc-1,dimx)
            CALL cf90_def_dim(iodat,'y',nr-1,dimy)
         ELSE
            CALL cf90_def_dim(iodat,'lon',nc-1,dimx)
            CALL cf90_def_dim(iodat,'lat',nr-1,dimy)
         ENDIF
      ELSE
         CALL cf90_def_dim(iodat,'nc',nc-1,dimx)
         CALL cf90_def_dim(iodat,'nr',nr-1,dimy)
      ENDIF
      CALL cf90_def_dim(iodat,'lenname',13,dimname)
      CALL cf90_def_dim(iodat,'trajectory',nopartout,dimp)
   ELSE
      dimt = 0
   ENDIF

!   
!2.4 Coordinate variables   
!------------------------
!
!2.4.1 Initialise
!----------------
!

   ALLOCATE (coordatts(nocoords),STAT=errstat)
   CALL error_alloc_struc('coordatts',1,(/nocoords/),'VariableAtts')
   CALL varatts_init(coordatts)

!
!2.4.2 Attributes
!----------------
!
!  ---data type
   IF (outform.EQ.'N') THEN
      coordatts(1)%data_type = double_NF90
      coordatts(2:nocoords)%data_type = MERGE(real_NF90,double_NF90,&
                                            & floattype.EQ.'S')
   ELSE
      coordatts(1)%data_type = rlong_type
      coordatts(2:nocoords)%data_type = MERGE(real_type,rlong_type,&
                                            & floattype.EQ.'S')
   ENDIF
   
!  ---fill value
   coordatts(1:3)%fill = .FALSE.
   coordatts(4)%fill = filepars%fill
   coordatts(4)%fill_value = fill_value
   
!  ---f90_name, standard_name, long_name, units
   CALL inquire_var(iarr_time,f90_name=coordatts(1)%f90_name,&
                  & standard_name=coordatts(1)%standard_name,&
                  & long_name=coordatts(1)%long_name)
   coordatts(1)%units = TRIM(tunits)
   CALL inquire_var(iarr_xout,f90_name=coordatts(2)%f90_name,&
                  & standard_name=coordatts(2)%standard_name,&
                  & long_name=coordatts(2)%long_name,&
                  & units=coordatts(2)%units)
   CALL inquire_var(iarr_yout,f90_name=coordatts(3)%f90_name,&
                  & standard_name=coordatts(3)%standard_name,&
                  & long_name=coordatts(3)%long_name,&
                  & units=coordatts(3)%units)
   CALL inquire_var(iarr_depout,&
                  & f90_name=coordatts(4)%f90_name,&
                  & standard_name=coordatts(4)%standard_name,&
                  & long_name=coordatts(4)%long_name,&
                  & units=coordatts(4)%units)

!  ---rank, dimids, other atrributes
   coordatts(1)%nrank = 1
   coordatts(1)%global_dims(1) = nstepout
   coordatts(1)%dimids(1) = dimt
   IF (iopt_grid_htype.LE.2) THEN
      coordatts(2:3)%nrank = 1
      coordatts(2)%global_dims(1) = nc-1
      coordatts(2)%dimids(1) = dimx
      coordatts(3)%global_dims(1) = nr-1
      coordatts(3)%dimids(1) = dimy
   ELSE
      coordatts(2:3)%nrank = 2
      coordatts(2)%global_dims(1:2) = (/nc-1,nr-1/)
      coordatts(2)%dimids(1:2) = (/dimx,dimy/)
      coordatts(2)%axis = 'X'
      coordatts(3)%global_dims(1:2) = (/nc-1,nr-1/)
      coordatts(3)%dimids(1:2) = (/dimx,dimy/)
      coordatts(3)%axis = 'Y'
   ENDIF
   coordatts(4)%nrank = 2
   coordatts(4)%global_dims(1:2) = (/nc-1,nr-1/)
   coordatts(4)%dimids(1:2) = (/dimx,dimy/)

!
!2.4.3 Define coordinate variables (netCDF)
!------------------------------------------
!

   IF (outform.EQ.'N') THEN
      CALL cf90_def_var(iodat,coordatts(1)%f90_name,double_NF90,(/dimt/),&
                      & filepars%timeid)
      ivar_243: DO ivar=2,nocoords
         nrank = coordatts(ivar)%nrank
         CALL cf90_def_var(iodat,coordatts(ivar)%f90_name,&
                         & coordatts(ivar)%data_type,&
                         & coordatts(ivar)%dimids(1:nrank),varid)
      ENDDO ivar_243
      CALL cf90_def_var(iodat,'trajectory',char_NF90,(/dimname/),varid)
   ENDIF

!  ---number of coordinates
   filepars%nocoords = MERGE(nocoords+1,nocoords,filepars%form.EQ.'N')

!
!2.4.4 Store output dimensions
!-----------------------------
!
!  ---number of dimensions
   nodims = 4

!  ---allocate
   ALLOCATE (outdims(nodims),STAT=errstat)
   CALL error_alloc('outdims',1,(/nodims/),kndint)
   ALLOCATE (dimnames(nodims),STAT=errstat)
   CALL error_alloc('dimnames',1,(/nodims/),kndchar,lenstr=lenname)

!  ---values and names
   outdims(1) = nstepout
   dimnames(1) = 'time'
   IF (iopt_grid_htype.LE.2) THEN
      dimnames(2) = MERGE ('x  ','lon',iopt_grid_sph.EQ.0)
      dimnames(3) = MERGE ('y  ','lat',iopt_grid_sph.EQ.0)
   ELSE
      dimnames(2) = 'nc'; dimnames(3) = 'nr'
   ENDIF
   outdims(2:3) = (/nc-1,nr-1/)
   dimnames(4) = 'trajectory'
   outdims(4) = nopartout
   
!
!2.5 Output variables   
!--------------------
!
!  ---attributes
   IF (outform.EQ.'N') THEN
      outpvars(iset,:)%data_type = MERGE(real_NF90,double_NF90,floattype.EQ.'S')
   ELSE
      outpvars(iset,:)%data_type = MERGE(real_type,rlong_type,floattype.EQ.'S')
   ENDIF
   outpvars(iset,:)%global_dims(1) = nstepout
   outpvars(iset,:)%dimids(1) = dimp
   outpvars(iset,:)%dimids(2) = dimt
   outpvars(iset,:)%fill = filepars%fill
   outpvars(iset,:)%fill_value = fill_value

!  ---define (netCDF)   
   IF (outform.EQ.'N') THEN
      ivar_250: DO ivar=1,novars
         CALL cf90_def_var(iodat,outpvars(iset,ivar)%f90_name,&
                         & outpvars(iset,ivar)%data_type,&
                         & outpvars(iset,ivar)%dimids(1:2),varid)
      ENDDO ivar_250
   ENDIF

!
!2.6 Write metadata
!------------------
!   
!2.6.1 ASCII/binary/info
!-----------------------
!

   ifil_261: DO ifil=1,2

      IF (ifil.EQ.1.AND.outform.EQ.'N') CYCLE ifil_261
      IF (ifil.EQ.2.AND.(.NOT.info)) EXIT ifil_261
      iunit = MERGE(iodat,ioinfo,ifil.EQ.1)
      dataform = MERGE(outform,'A',ifil.EQ.1) 
   
!     ---global attributes
      noglbatts = COUNT(LEN_TRIM(glbatts%value).GT.0)
      CALL conv_to_chars(cvals(1),noglbatts)
      CALL write_metadata_line(iunit,cvals(1:1),'nGlbatts',dataform)
      iglb_2611: DO iglb=1,numglbatts
         IF (TRIM(glbatts(iglb)%value).NE.''.AND.&
           & TRIM(glbatts(iglb)%value).NE.'netcdf') THEN
            cvals(1) = TRIM(glbatts(iglb)%value)
            CALL write_metadata_line(iunit,cvals(1:1),&
                                   & TRIM(glbatts(iglb)%name),dataform)
         ENDIF
      ENDDO iglb_2611
      cvals(1) = 'trajectory'
      CALL write_metadata_line(iunit,cvals(1:1),'featureType',dataform)

!     ---output info
      CALL conv_to_chars(cvals(1),1)
      CALL write_metadata_line(iunit,cvals(1:1),'nodim',dataform)
      CALL conv_to_chars(cvals(1),nocoords)
      CALL write_metadata_line(iunit,cvals(1:1),'nCoordinates',dataform)
      CALL conv_to_chars(cvals(1),novars)
      CALL write_metadata_line(iunit,cvals(1:1),'nVariables',dataform)
      cvals(1) = floattype
      CALL write_metadata_line(iunit,cvals(1:1),'Floattype',dataform)
      CALL conv_to_chars(cvals(1),nodims)
      CALL write_metadata_line(iunit,cvals(1:1),'nDimensions',dataform)
      CALL write_metadata_line(iunit,cvals(1:nodims),'Dimnames',dataform)
      CALL conv_to_chars(cvals(1:nodims),outdims(1:nodims))
      CALL write_metadata_line(iunit,cvals(1:nodims),'Dimensions',dataform)
      CALL write_metadata_line(iunit,cvals(1:1),'Dimensions',dataform)
      cvals(1) = outppars(iset)%startdate
      CALL write_metadata_line(iunit,cvals(1:1),'StartDate',dataform)
      CALL conv_to_chars(cvals(1),deltout)
      CALL write_metadata_line(iunit,cvals(1:1),'time_step',dataform)
      cvals(1) = TRIM(coordatts(1)%units)
      CALL write_metadata_line(iunit,cvals(1:1),'time_format',dataform)

!     ---coordinate attributes
      ivar_2622: DO ivar=1,nocoords
         CALL conv_to_chars(cvals(1),ivar)
         nvals = 11 + coordatts(ivar)%nrank
         CALL conv_to_chars_modvars(cvals(2:nvals),coordatts(ivar))
         CALL write_metadata_line(iunit,cvals(1:nvals),'coordatts',dataform)
      ENDDO ivar_2622

!     ---variable attributes
      ivar_2612: DO ivar=1,novars
         CALL conv_to_chars(cvals(1),nocoords+ivar)
         nvals = 11 + outpvars(iset,ivar)%nrank
         CALL conv_to_chars_modvars(cvals(2:nvals),outpvars(iset,ivar))
         CALL write_metadata_line(iunit,cvals(1:nvals),'varatts',dataform)
      ENDDO ivar_2612

!     ---end of header
      cline = '#'
      IF (dataform.EQ.'A') THEN
         WRITE (iunit,'(A)') TRIM(cline)
      ELSEIF (dataform.EQ.'U') THEN
         WRITE (iunit) cline
      ENDIF

   ENDDO ifil_261

!  ---close
   IF (info) CALL close_file(ioinfo,'A',TRIM(infofile))

!
!2.6.2 Netcdf
!------------
!

   IF (outform.EQ.'N') THEN

!     ---global attributes
      iglb_621: DO iglb=1,numglbatts
         IF (TRIM(glbatts(iglb)%value).NE.'') THEN
            l = LEN_TRIM(glbatts(iglb)%value)
            CALL cf90_put_att_chars(iodat,global_NF90,&
                                  & TRIM(glbatts(iglb)%name),l,&
                                  & TRIM(glbatts(iglb)%value))
         ENDIF
      ENDDO iglb_621
      CALL cf90_put_att_chars(iodat,global_NF90,'featureType',10,'trajectory')

!     ---coordinate attributes
      ivar_622: DO ivar=1,nocoords
         CALL cf90_put_att_chars(iodat,ivar,'standard_name',lendesc,&
                               & coordatts(ivar)%standard_name)
         CALL cf90_put_att_chars(iodat,ivar,'long_name',lendesc,&
                               & coordatts(ivar)%long_name)
         IF (TRIM(coordatts(ivar)%units).NE.'') THEN
            CALL cf90_put_att_chars(iodat,ivar,'units',lenunit,&
                                  & coordatts(ivar)%units)
         ENDIF
         IF (coordatts(ivar)%fill) THEN
            CALL cf90_put_att(iodat,ivar,'_FillValue',&
                            & coordatts(ivar)%fill_value)
            CALL cf90_put_att(iodat,ivar,'missing_value',&
                            & coordatts(ivar)%fill_value)
         ENDIF
         IF (TRIM(coordatts(ivar)%axis).NE.'') THEN
            CALL cf90_put_att_chars(iodat,ivar,'axis',1,&
                                  & coordatts(ivar)%axis)
         ENDIF
!        --other
         IF (ivar.EQ.1) THEN
            CALL cf90_put_att_chars(iodat,1,'calendar',9,'gregorian')
         ENDIF
      ENDDO ivar_622
      CALL cf90_put_att_chars(iodat,nocoords+1,'cf_role',13,'trajectory_id')

!     ---variable attributes
      ivar_262: DO ivar=1,novars
         varid = filepars%nocoords + ivar
         CALL cf90_put_att_chars(iodat,varid,'standard_name',lendesc,&
                               & outpvars(iset,ivar)%standard_name)
         CALL cf90_put_att_chars(iodat,varid,'long_name',lendesc,&
                               & outpvars(iset,ivar)%long_name)
         CALL cf90_put_att_chars(iodat,varid,'units',lenunit,&
                               & outpvars(iset,ivar)%units)
         IF (filepars%fill) THEN
            IF (floattype.EQ.'S') THEN
               CALL cf90_put_att(iodat,varid,'_FillValue',real_fill)
               CALL cf90_put_att(iodat,varid,'missing_value',real_fill)
            ELSE
               CALL cf90_put_att(iodat,varid,'_FillValue',double_fill)
               CALL cf90_put_att(iodat,varid,'missing_value',double_fill)
            ENDIF
         ENDIF
      ENDDO ivar_262

!     ---leave define mode
      CALL cf90_enddef(iodat)
      
   ENDIF

!
!2.7 Abort if needed
!--------------------
!

   IF (nerrs.GT.0) CALL error_file(ierrno_write,filepars=filepars)

!
!2.8 Deallocate
!--------------
!   

   DEALLOCATE (coordatts,dimnames,outdims)
   
!
!2.9 Store
!---------   
!

   part1d(iset) = filepars

ENDDO iset_200

!
!3. Write horizontal coordinate arrays
!-------------------------------------
!
!3.1 Allocate
!------------
!

ALLOCATE (xout(nc-1,nr-1),STAT=errstat)
CALL error_alloc('xout',2,(/nc-1,nr-1/),kndrtype)
ALLOCATE (yout(nc-1,nr-1),STAT=errstat)
CALL error_alloc('yout',2,(/nc-1,nr-1/),kndrtype)
ALLOCATE (depout(nc-1,nr-1),STAT=errstat)
CALL error_alloc('depout',2,(/nc-1,nr-1/),kndrtype)
IF (SIZE(depout).GT.0) depout = fill_value

!
!3.2 Output grid
!---------------
!

CALL UVarr_at_C(gxcoord(1:nc,1:nr),xout,0,0,(/1,1,nz/),&
             & (/nc,nr,nz/),1,iarr_gxcoord,.TRUE.)
CALL UVarr_at_C(gycoord(1:nc,1:nr),yout,0,0,(/1,1,nz/),&
             & (/nc,nr,nz/),1,iarr_gycoord,.TRUE.)
depout = MERGE(depmeanatc(1:nc-1,1:nr-1),fill_value,maskatc(1:nc-1,1:nr-1))

!
!3.3 Write
!---------
!

iset_330: DO iset=1,nosetspart
   filepars = part1d(iset)
   IF (filepars%defined) THEN
      IF (iopt_grid_htype.LE.2) THEN
         CALL write_vars(xout(:,1),filepars,2)
         CALL write_vars(yout(1,:),filepars,3)
      ELSE
         CALL write_vars(xout,filepars,2)
         CALL write_vars(yout,filepars,3)
      ENDIF
      CALL write_vars(depout,filepars,4)
   ENDIF
ENDDO iset_330

!
!3.4 Deallocate
!--------------
!

DEALLOCATE (xout,yout,depout)

CALL log_timer_out()


RETURN

END SUBROUTINE particle_trajects_init
