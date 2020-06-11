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
! *Particle_Cloud* Discharge module for the particle model
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Cloud.f90  V2.11
!
! $Date: 2018-06-14 15:25:04 +0200 (Thu, 14 Jun 2018) $
!
! $Revision: 1155 $
!
! Description - particle clouds are generated at discharge locations
!
! Reference -
!
! Routines - define_parcld_data, particle_cloud, read_parcld_data,
!            write_parcld_data
!
!************************************************************************
!

!========================================================================

SUBROUTINE define_parcld_data(ifil,ciodatetime,discoords)
!************************************************************************
!
! *define_parcld_data* Obtain locations for particle cloud releases
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Cloud.f90 V2.11
!
! Description -
!
! Reference -
!
! Calling program -
!
! External calls - read_parcld_data, usrdef_parcld_data, write_parcld_data
!
! Module calls - error_file, suspend_proc
!
!************************************************************************
!
USE iopars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: ifil
REAL (KIND=kndrlong), INTENT(INOUT), DIMENSION(3) :: discoords

!
! Name         Type    Purpose
!----------------------------------------------------------------------------
!*ifil*        INTEGER Cloud number
!*ciodatetime* CHAR    Date/time in data file
!*discoords*   REAL    Location of cloud center
!
!----------------------------------------------------------------------------
!


procname(pglev+1) = 'define_cloud_data'
CALL log_timer_in()

discoords = 0.0


!
!1. Initialise
!-------------
!

discoords = 0.0

!
!2. Open file on first call
!--------------------------
!


IF (modfiles(io_parcld,ifil,1)%iostat.EQ.0) THEN

!  ---standard
   IF (modfiles(io_parcld,ifil,1)%status.EQ.'R') THEN
      CALL read_parcld_data(ifil,ciodatetime,discoords)
!  ---user-defined
   ELSEIF (modfiles(io_parcld,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_parcld_data(ifil,ciodatetime,discoords)
   ENDIF

ENDIF   

!
!3. Read data
!------------
!
!---standard
IF (modfiles(io_parcld,ifil,1)%status.EQ.'R') THEN
   CALL read_parcld_data(ifil,ciodatetime,discoords)
!---user-defined
ELSEIF (modfiles(io_parcld,ifil,1)%status.EQ.'N') THEN
   CALL usrdef_parcld_data(ifil,ciodatetime,discoords)
ENDIF
 
!
!2. Write data
!-------------
!

IF (modfiles(io_parcld,ifil,2)%status.EQ.'W') THEN
   CALL write_parcld_data(ifil,ciodatetime,discoords)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE define_parcld_data

!========================================================================

SUBROUTINE particle_cloud
!************************************************************************
!
! *particle_cloud* Generate initial conditions for the particle model in the
!                  form of cloud discharge
!
! Author -  Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Cloud.f90  V2.11
!
! $Date: 2018-06-14 15:25:04 +0200 (Thu, 14 Jun 2018) $
!
! $Revision: 1155 $
!
! Description - cloud discharge are enabled when iopt_part_cloud = 1
!
! Reference -
!
! Calling program - particle_model
!
! External calls - define_parcld_data
!
! Module calls - Cint_at_p, diff_dates, earlier, error_abort,
!                hrel_coords_curv_abs, hrel_coords_curv_rel,
!                hrel_coords_rect_nonunif_abs, hrel_coords_rect_nonunif_rel,
!                hrel_coords_rect_unif_abs, hrel_coords_rect_unif_rel, later,
!                rng_close, rng_open, rng_uniform_var
!
!************************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE partpars
USE partswitches
USE partvars
USE physpars
USE switches
USE syspars
USE timepars
USE error_routines, ONLY: error_abort
USE grid_interp, ONLY: hrel_coords_curv_abs, hrel_coords_curv_rel, &
                     & hrel_coords_rect_nonunif_abs, &
                     & hrel_coords_rect_nonunif_rel, &
                     & hrel_coords_rect_unif_abs, hrel_coords_rect_unif_rel
USE particle_routines, ONLY: Cint_at_p
USE rng_library, ONLY: rng_close, rng_open, rng_uniform_var  
USE time_routines
  
IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flag
CHARACTER (LEN=lentime) :: ciodatetime
INTEGER :: i, icld, icoord, irel, istat, j, jcoord, k, lbuf1, lbuf2, &
         & nopartcld, noreleases, ntstart, numgen, p, state
REAL :: cosangle, delxdat, delydat, depmeancen, depmeanp, deptotcen, deptotp, &
      & displx, disply, gsigcen, phiran, rran, sellip, sinangle, spos, &
      & thetaran, x, xcoord, xellip, xfac, x0dat, y, ycoord, yellip, yfac, &
      & y0dat, zellip
REAL (KIND=kndrlong) :: tstart, xcen, xpos, ycen, ypos, zpos, zcen, ztmp
REAL (KIND=kndrlong), DIMENSION(3) :: discoords


procname(pglev+1) = 'particle_cloud'
CALL log_timer_in()

!
!1. Initialise
!-------------
!
!----random generator
CALL rng_open(numgen)

!---work space
x0dat = surfacegrids(igrd_model,1)%x0dat
y0dat = surfacegrids(igrd_model,1)%y0dat
delxdat = surfacegrids(igrd_model,1)%delxdat
delydat = surfacegrids(igrd_model,1)%delydat

!
!2. Define particle clouds
!-------------------------
!

icld_200: DO icld=1,noclouds

!
!2.1 Initialise parameters
!-------------------------
!

   nopartcld = cloud(icld)%nopart
   noreleases = cloud(icld)%noreleases
   state = MERGE(1,cloud(icld)%state,iopt_grid_nodim.EQ.2)
   xellip = 0.5*cloud(icld)%width
   yellip = 0.5*cloud(icld)%length
   zellip = 0.5*cloud(icld)%thick
   cosangle = COS(degtorad*cloud(icld)%orientation)
   sinangle = SIN(degtorad*cloud(icld)%orientation)
   IF (icld.EQ.1) THEN
      lbuf1 = 1; lbuf2 = cloud(1)%nopart
   ENDIF
  
!
!2.2 Generate particles at each release
!--------------------------------------
!   

   irel_220: DO irel=1,noreleases
      part(lbuf1:lbuf2)%label = cloud(icld)%label
      part(lbuf1:lbuf2)%state = state

!
!2.2.1 Define release
!--------------------
!      

      flag = .FALSE.
      CALL define_parcld_data(icld,ciodatetime,discoords)

!
!2.2.2 Check whether release is within simulated time interval and at sea
!------------------------------------------------------------------------
!

      IF (for_drift) THEN
         flag = ciodatetime.earlier.CStartDateTime
      ELSE
         flag = ciodatetime.later.CStartDateTime
      ENDIF
      IF (flag) THEN
         part(lbuf1:lbuf2)%drift_state = 4
         GOTO 999
      ENDIF
      xcen = discoords(1); ycen = discoords(2); ztmp = discoords(3)
      IF (iopt_grid_sph.EQ.0) THEN
         xfac = 1.0; yfac = 1.0
      ELSE
         xfac = degtorad/(Rearth*COS(degtorad*ycen))
         yfac = degtorad/Rearth
      ENDIF
      SELECT CASE (iopt_grid_htype)
         CASE (1)
            CALL hrel_coords_rect_unif_abs(xcen,ycen,x0dat,y0dat,delxdat,&
                                         & delydat,(/1,1/),(/nc,nr/),istat,i,j,&
                                         & xcoord=x,ycoord=y)
         CASE (2)
            CALL hrel_coords_rect_nonunif_abs(xcen,ycen,gxcoord(1:nc,1),&
                                            & gycoord(1,1:nr),nc,nr,(/1,1/),&
                                            & (/nc,nr/),istat,i,j,xcoord=x,&
                                            & ycoord=y)
         CASE (3)
            CALL hrel_coords_curv_abs(xcen,ycen,gxcoord(1:nc,1:nr),&
                                    & gycoord(1:nc,1:nr),nc,nr,(/1,1/),&
                                    & (/nc,nr/),istat,2,i,j,&
                                    & xcoord=x,ycoord=y)
      END SELECT
      IF (istat.EQ.0.AND.(.NOT.maskatc(i,j))) THEN
         part(lbuf1:lbuf2)%drift_state = 3
      ELSEIF (istat.EQ.1) THEN
         part(lbuf1:lbuf2)%drift_state = 4
         GOTO 999
      ELSE
         part(lbuf1:lbuf2)%drift_state = 1
      ENDIF

!
!2.2.3 Release time
!------------------      
!

      IF (for_drift) THEN
         CALL diff_dates(CStartDateTime,ciodatetime,1,rtime=tstart)
         ntstart = NINT(tstart/delt2d)
      ELSE
         CALL diff_dates(ciodatetime,CStartDateTime,1,rtime=tstart)
         ntstart = NINT(tstart/delt2d)
      ENDIF
      CALL diff_dates(RefDateTime_part,ciodatetime,ptime_unit,rtime=tstart)
      part(lbuf1:lbuf2)%ntstart = ntstart
      part(lbuf1:lbuf2)%tstart = tstart

!
!2.2.4 Verical location of cloud centre
!--------------------------------------      
!

      IF (state.EQ.2) THEN
         deptotcen = Cint_at_p(deptotatc(i-1:i+1,j-1:j+1),i,j,x,y)
         depmeancen = Cint_at_p(depmeanatc(i-1:i+1,j-1:j+1),i,j,x,y)
         SELECT CASE (cloud(icld)%kdistype)
            CASE (1)
               k = NINT(ztmp)
               IF (iopt_grid_vtype.LT.3) THEN
                  gsigcen = gsigcoordatc(k)
               ELSE
                  gsigcen = gscoordatc(i,j,k)
               ENDIF
               zcen = gsigcen*deptotatc(i,j) - depmeanatc(i,j)
            CASE (2)
               ztmp = MAX(0.0,MIN(ztmp,deptotcen))
               gsigcen = ztmp/deptotcen
               zcen = ztmp - depmeancen
            CASE (3)
               ztmp = deptotcen - ztmp
               ztmp = MAX(0.0,MIN(ztmp,deptotcen))
               gsigcen = ztmp/deptotcen
               zcen = ztmp - depmeancen
         END SELECT
         sellip = zellip/deptotcen
      ENDIF

!      
!2.2.5 Define cloud release particles
!------------------------------------
!

      p_225: DO p=lbuf1,lbuf2

!        ---first particle release at cloud centre         
         IF (p.EQ.lbuf1) THEN
            part(p)%icoord = i
            part(p)%jcoord = j
            part(p)%xcoord = x
            part(p)%ycoord = y
            part(p)%xpos = xcen
            part(p)%ypos = ycen
            part(p)%center = .TRUE.
            IF (state.EQ.2) part(p)%zpos = zcen

!        ---cloud releases            
         ELSE
            part(p)%center  = .FALSE.
!           --horizontal locations       
            flag = .FALSE.
            DO WHILE (.NOT.flag) 
               CALL rng_uniform_var(rran,numgen,0.0,1.0)
               CALL rng_uniform_var(phiran,numgen,0.0,twopi)
               CALL rng_uniform_var(thetaran,numgen,0.0,twopi)
               xpos = rran*xellip*COS(phiran)*SIN(thetaran)
               ypos = rran*yellip*SIN(phiran)*SIN(thetaran)
               displx = xpos*cosangle + ypos*sinangle
               disply = ypos*cosangle - xpos*sinangle
               xpos = xcen + displx
               ypos = ycen + disply
               icoord = i; jcoord = j; xcoord = x; ycoord = y
               SELECT CASE (iopt_grid_htype)
                  CASE (1)
                     CALL hrel_coords_rect_unif_rel(xcen,ycen,displx,disply,&
                                    & delxdat,delydat,(/1,1/),&
                                    & (/nc,nr/),istat,icoord,jcoord,&
                                    & xcoord=xcoord,ycoord=ycoord)
                  CASE (2)
                     CALL hrel_coords_rect_nonunif_rel(xcen,ycen,displx,disply,&
                                    & gxcoord(1:nc,1),gycoord(1,1:nr),nc,nr,&
                                    & (/1,1/),(/nc,nr/),istat,icoord,jcoord,&
                                    & xcoord=xcoord,ycoord=ycoord)
                  CASE (3)
                     CALL hrel_coords_curv_rel(xcen,ycen,displx,disply,&
                                    & gxcoord(1:nc,1:nr),gycoord(1:nc,1:nr),&
                                    & nc,nr,(/1,1/),(/nc,nr/),nsize_sb,istat,2,&
                                    & icoord,jcoord,xcoord=xcoord,ycoord=ycoord)
                  END SELECT
                  IF (istat.EQ.0) THEN
                     flag = maskatc(icoord,jcoord)
                  ENDIF
            ENDDO
!           --vertical locations            
            IF (state.EQ.2) THEN
               spos = rran*sellip*COS(thetaran)
               IF (icoord.NE.i.OR.jcoord.NE.j) THEN
                  deptotp = Cint_at_p(deptotatc(i-1:i+1,j-1:j+1),i,j,x,y)
                  depmeanp = Cint_at_p(depmeanatc(i-1:i+1,j-1:j+1),i,j,x,y)
               ELSE
                  deptotp = deptotcen; depmeanp = depmeancen
               ENDIF
               spos = gsigcen + spos
               spos = MAX(0.0,MIN(spos,1.0))
               zpos = spos*deptotp - depmeanp
            ENDIF
!           --store            
            part(p)%icoord = icoord
            part(p)%jcoord = jcoord
            part(p)%xcoord = xcoord
            part(p)%ycoord = ycoord
            part(p)%xpos = xpos
            part(p)%ypos = ypos
            IF (state.EQ.2) part(p)%zpos = zpos

         ENDIF

      ENDDO p_225
         
!     ---update particle index buffer
999   CONTINUE
      lbuf1 = lbuf2 + 1
      lbuf2 = lbuf2 + nopartcld
   ENDDO irel_220

   CALL error_abort('particle_cloud',ierrno_inival)
   
ENDDO icld_200

!
!3. Finalise
!-----------
!
!---close random generator
CALL rng_close(numgen)

CALL log_timer_out()


RETURN

END SUBROUTINE particle_cloud

!=================================================================

SUBROUTINE read_parcld_data(ifil,ciodatetime,disdata)
!*****************************************************************
!
! *read_parcld_data* Read discharge locations for particle cloud releases in
!                    standard format
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Cloud.f90 V2.11
!
! Description -
!
! Reference -
!
! Calling program - define_parcld_data
!
! Module calls - error_alloc_struc, open_filepars, read_glbatts_mod, read_time,
!                read_varatts_mod, read_vars
!
!*****************************************************************
!
USE datatypes
USE iopars
USE syspars
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: open_filepars, read_glbatts_mod, &
                        & read_time, read_varatts_mod, read_vars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(OUT) :: ciodatetime
INTEGER, INTENT(IN) :: ifil
REAL (KIND=kndrlong), INTENT(INOUT), DIMENSION(3) :: disdata

!
! Name         Type Purpose
!------------------------------------------------------------------------------
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*disdata*     REAL    Input data
!
!------------------------------------------------------------------------------
!
!* Local variables
!
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'read_parcld_data'
CALL log_timer_in()

filepars = modfiles(io_parcld,ifil,1)

!
!1. File Header on first call
!----------------------------
!

IF (filepars%iostat.LE.0) THEN

!  ---open file
   filepars = modfiles(io_parcld,ifil,1)
   CALL open_filepars(filepars)
   
!  ---global attributes
   CALL read_glbatts_mod(filepars)
 
!  ---allocate
   ALLOCATE (varatts(2),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/2/),'VariableAtts')

!  ---variable attributes
   CALL read_varatts_mod(filepars,varatts,2)

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

IF (filepars%iostat.LT.3) CALL read_vars(disdata,filepars,2)

1000 modfiles(io_parcld,ifil,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_parcld_data

!=================================================================

SUBROUTINE write_parcld_data(ifil,ciodatetime,disdata)
!*****************************************************************
!
! *write_parcld_data* Write discharge discharge locations for particle cloud
!                     model in standard format
!
! Author - Valerie Duliere and Patrick Luyten
!
! Version - @(COHERENS)Particle_Cloud.f90 V2.11
!
! Description -
!
! Reference -
!
! Calling program - define_parcld_data
!
! Module calls - error_alloc_struc, open_filepars, output_flag, 
!                set_modfiles_atts, set_modvars_atts, write_atts_mod,
!                write_time, write_vars
!
!*****************************************************************
!
USE datatypes
USE iopars
USE syspars
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: open_filepars, output_flag, write_atts_mod, &
                        & write_time, write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: ciodatetime
INTEGER, INTENT(IN) :: ifil
REAL (KIND=kndrlong), INTENT(IN), DIMENSION(3) :: disdata

!
! Name         Type Purpose
!------------------------------------------------------------------------------
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*disdata*     REAL    Input data
!
!------------------------------------------------------------------------------
!
!* Local variables
!
LOGICAL :: flag
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:) :: varatts


procname(pglev+1) = 'write_parcld_data'
CALL log_timer_in()

!
!1. Disable writing if needed
!----------------------------
!

filepars = modfiles(io_parcld,ifil,2)
flag = output_flag(ciodatetime,filepars%tskips)
IF (.NOT.flag) GOTO 1000

!
!2. Write header on first call
!-----------------------------
!

IF (filepars%iostat.EQ.0) THEN

!  ---file attributes
   CALL set_modfiles_atts(io_parcld,ifil,2)
   filepars = modfiles(io_parcld,ifil,2)

!  ---open file
   CALL open_filepars(filepars)
   
!  ---allocate
   ALLOCATE (varatts(2),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/4/),'VariableAtts')

!  ---variable attributes
   CALL set_modvars_atts(io_parcld,ifil,2,varatts,2)

!  ---write
   CALL write_atts_mod(filepars,varatts,2)

!  ---deallocate
   DEALLOCATE (varatts)
   
ENDIF
 
!
!3. Write data of discharges
!---------------------------
!
!---date/time
CALL write_time(ciodatetime,filepars)

!---data
CALL write_vars(disdata,filepars,2)

modfiles(io_parcld,ifil,2) = filepars

1000 CALL log_timer_out()


RETURN

END SUBROUTINE write_parcld_data
