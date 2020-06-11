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
! *Usrdef_Time_Series* Define time series output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.11.2
!
! $Date: 2018-07-23 16:55:25 +0200 (Mon, 23 Jul 2018) $
!
! $Revision: 1171 $
!
! Description - test case isotest
!
! Reference -
!
! Routines - usrdef_tsr_params, usrdef_tsr0d_vals, usrdef_tsr2d_vals,
!            usrdef_tsr3d_vals
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_tsr_params
!************************************************************************
!
! *usrdef_tsr_params* Specifiers for time series output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90 V2.11.2
!
! Description - test case isotest
!
! Reference -
!
! Calling program - time_series_init
!
!************************************************************************
!
USE iopars
USE modids
USE switches
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: zflag
INTEGER :: l, ntmod


procname(pglev+1) = 'usrdef_tsr_params'
CALL log_timer_in()

zflag = runtitle(8:8).GT.'C'

!
!1. Output variables
!-------------------
!
!---fortran name
tsrvars(4)%f90_name = 'xdifflux_temp'
tsrvars(5)%f90_name = 'xdifflux_sal'
tsrvars(6)%f90_name = 'xdifflux_dens'
IF (zflag) THEN
   tsrvars(7)%f90_name = 'zdifflux_temp'
   tsrvars(8)%f90_name = 'zdifflux_sal'
   tsrvars(9)%f90_name = 'zdifflux_dens'
ENDIF

!---standard name
tsrvars(4:novarstsr)%standard_name = ''

!---comment
tsrvars(4:novarstsr)%comment = ''

!---long name
tsrvars(4)%long_name = 'X-component of temperature diffusive flux'
tsrvars(5)%long_name = 'X-component of salinity diffusive flux'
tsrvars(6)%long_name = 'X-component of density diffusive flux'
IF (zflag) THEN
   tsrvars(7)%long_name = 'Vertical component of temperature diffusive flux'
   tsrvars(8)%long_name = 'Vertical component of salinity diffusive flux'
   tsrvars(9)%long_name = 'Vertical component of density diffusive flux'
ENDIF

!---units
tsrvars(4)%units = 'degC m s-1'
tsrvars(5)%units = 'PSU m s-1'
tsrvars(6)%units = 'kg m-2 s-1'
IF (zflag) THEN
   tsrvars(7)%units = 'degC m s-1'
   tsrvars(8)%units = 'PSU m s-1'
   tsrvars(9)%units = 'kg m-2 s-1'
ENDIF

!---varids and ranks
tsrvars(1:6)%ivarid = (/iarr_temp,iarr_sal,iarr_dens,0,0,0/)
IF (zflag) tsrvars(7:9)%ivarid = (/0,0,0/)
tsrvars%nrank = 3

!
!2. Variable indices
!-------------------
!

ivarstsr(1,1:novarstsr) = (/(l,l=1,novarstsr)/)

!
!3. File parameters
!------------------
!

tsr3d(1)%defined = .TRUE.

!
!4. Output grid
!--------------
!

tsrgpars(1)%tlims = (/0,nstep,8640/)
tsrgpars(1)%time_format = 4

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_tsr_params

!========================================================================

SUBROUTINE usrdef_tsr0d_vals(out0ddat,n0vars)
!************************************************************************
!
! *usrdef_tsr0d_vals* 0-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - time_series
!
! Module calls - energy_0d, max_vars, vector_mag_arr_atc
!
!************************************************************************
!

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: n0vars
REAL, INTENT(OUT), DIMENSION(n0vars) :: out0ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out0ddat* REAL    output data
!*n0vars*   INTEGER number of 0-D variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_tsr0d_vals

!========================================================================

SUBROUTINE usrdef_tsr2d_vals(out2ddat,i,j,n2vars)
!************************************************************************
!
! *usrdef_tsr2d_vals* 2-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.1.2
!
! Description -
!
! Reference -
!
! Calling program - time_series
!
!************************************************************************
!

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, n2vars
REAL, INTENT(OUT), DIMENSION(n2vars) :: out2ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out2ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*n2vars*   INTEGER number of 2-D variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_tsr2d_vals

!========================================================================

SUBROUTINE usrdef_tsr3d_vals(out3ddat,i,j,k,n3vars)
!************************************************************************
!
! *usrdef_tsr3d_vals* 3-D output data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Time_Series.f90  V2.11.2
!
! Description -
!
! Reference -
!
! Calling program - time_series
!
!************************************************************************
!
USE density
USE diffusion
USE grid
USE gridpars
USE iopars
USE physpars
USE switches

  
IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: i, j, k, n3vars
REAL, INTENT(OUT), DIMENSION(n3vars) :: out3ddat

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*out3ddat* REAL    output data
!*i*        INTEGER X-index of output location
!*j*        INTEGER Y-index of output location
!*k*        INTEGER Vertical index of output location
!*n3vars*   INTEGER number of 3-D variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL zflag
INTEGER :: iip, ip, kkp, kp, l
REAL :: hdifu, hdifw, slopefac
REAL, DIMENSION(2) :: xdiff, zdiff


zflag = runtitle(8:8).GT.'C'
hdifu = hdifscal_fac*hdifcoef3datu(i,j,k)
IF (iopt_hdif_scal.EQ.2) hdifw = hdifscal_fac*hdifcoef3datw(i,j,k)

out3ddat(4:6) = 0.0
IF (zflag) out3ddat(7:9) = 0.0

!
!1. Diffusive flux in X-direction
!--------------------------------
!

IF (node2du(i,j).EQ.2) THEN
   
!
!1.1  Laplacian diffusion
!------------------------
!   

   IF (iopt_hdif_scal.EQ.1) THEN
      
      out3ddat(4) = -hdifcoef3datu(i,j,k)*(temp(i,j,k)-temp(i-1,j,k))&
                  & /delxatu(i,j)
      out3ddat(5) = -hdifcoef3datu(i,j,k)*(sal(i,j,k)-sal(i-1,j,k))/delxatu(i,j)
      out3ddat(6) = -beta_temp(i,j,k)*out3ddat(4)+beta_sal(i,j,k)*out3ddat(5)

!      
!1.2 Isolevel diffusion
!----------------------
!
      
   ELSEIF (iopt_hdif_scal.EQ.2) THEN

!     ---temperature
      IF (k.EQ.1) THEN
         zdiff(1) = 0.5*(temp(i-1,j,2)+temp(i,j,2)-&
                       & temp(i-1,j,1)-temp(i,j,1))/delzatuw(i,j,2)
      ELSE
         zdiff(1) = 0.5*(temp(i-1,j,k)+temp(i,j,k)-&
                       & temp(i-1,j,k-1)-temp(i,j,k-1))/delzatuw(i,j,k)
      ENDIF
      IF (k.LT.nz) THEN
         zdiff(2) = 0.5*(temp(i-1,j,k+1)+temp(i,j,k+1)-&
                       & temp(i-1,j,k)-temp(i,j,k))/delzatuw(i,j,k+1)
      ELSE
         zdiff(2) = 0.5*(temp(i-1,j,nz)+temp(i,j,nz)-&
                       & temp(i-1,j,nz-1)-temp(i,j,nz-1))/delzatuw(i,j,nz)
      ENDIF
      out3ddat(4) = -hdifcoef3datu(i,j,k)*(temp(i,j,k)-temp(i-1,j,k))&
                  & /delxatu(i,j)
      out3ddat(4) = out3ddat(4) + xslopeatu_geo(i,j,k)*0.5*(zdiff(1)+zdiff(2))

!     ---salinity
      IF (k.EQ.1) THEN
         zdiff(1) = 0.5*(sal(i-1,j,2)+sal(i,j,2)-&
                       & sal(i-1,j,1)-sal(i,j,1))/delzatuw(i,j,2)
      ELSE
         zdiff(1) = 0.5*(sal(i-1,j,k)+sal(i,j,k)-&
                       & sal(i-1,j,k-1)-sal(i,j,k-1))/delzatuw(i,j,k)
      ENDIF
      IF (k.LT.nz) THEN
         zdiff(2) = 0.5*(sal(i-1,j,k+1)+sal(i,j,k+1)-&
                       & sal(i-1,j,k)-sal(i,j,k))/delzatuw(i,j,k+1)
      ELSE
         zdiff(2) = 0.5*(sal(i-1,j,nz)+sal(i,j,nz)-&
                       & sal(i-1,j,nz-1)-sal(i,j,nz-1))/delzatuw(i,j,nz)
      ENDIF
      out3ddat(5) = -hdifcoef3datu(i,j,k)*(sal(i,j,k)-sal(i-1,j,k))&
                  & /delxatu(i,j)
      out3ddat(5) = out3ddat(5) + xslopeatu_geo(i,j,k)*0.5*(zdiff(1)+zdiff(2))

!     ---density
      out3ddat(6) = -beta_temp(i,j,k)*out3ddat(4)+beta_sal(i,j,k)*out3ddat(5)

!      
!1.3 Isoneutral diffusion
!------------------------
!

   ELSEIF (iopt_hdif_scal.EQ.3) THEN
      
      kp_131: DO kp=0,1
      ip_131: DO ip=0,1

         iip = i+ip-1
         kkp = k+1-kp
            
         IF (kkp.GT.1.AND.kkp.LE.nz) THEN
            slopefac = 0.25*(hdifcoef3datu(i,j,k)*&
                     & xslopeatu_siso(iip,j,k,1-ip,kp))/delzatw(iip,j,kkp)
            out3ddat(4) = out3ddat(4) + &
                        & slopefac*(temp(iip,j,kkp)-temp(iip,j,kkp-1))
            out3ddat(5) = out3ddat(5) + &
                        & slopefac*(sal(iip,j,kkp)-sal(iip,j,kkp-1))
         ELSEIF (kkp.EQ.nz+1) THEN
            slopefac = 0.25*(hdifcoef3datu(i,j,nz-1)*&
                     & xslopeatu_siso(iip,j,nz-1,1-ip,kp))/delzatw(iip,j,nz)
            out3ddat(4) = out3ddat(4) + &
                        & slopefac*(temp(iip,j,nz)-temp(iip,j,nz-1))
            out3ddat(5) = out3ddat(5) + slopefac*(sal(iip,j,nz)-sal(iip,j,nz-1))
         ELSEIF (kkp.EQ.1) THEN
            slopefac = 0.25*(hdifcoef3datu(i,j,2)*&
                     & xslopeatu_siso(iip,j,2,1-ip,kp))/delzatw(iip,j,2)
            out3ddat(4) = out3ddat(4) + slopefac*(temp(iip,j,2)-temp(iip,j,1))
            out3ddat(5) = out3ddat(5) + slopefac*(sal(iip,j,2)-sal(iip,j,1))
         ENDIF

      ENDDO ip_131
      ENDDO kp_131

!     ---temperature
      out3ddat(4) =  out3ddat(4) - hdifcoef3datu(i,j,k)*&
                   & (temp(i,j,k)-temp(i-1,j,k))/delxatu(i,j)
!     ---salinity
      out3ddat(5) =  out3ddat(5) - hdifcoef3datu(i,j,k)*&
                   & (sal(i,j,k)-sal(i-1,j,k))/delxatu(i,j)
!     ---density      
      out3ddat(6) = -beta_temp(i,j,k)*out3ddat(4)+beta_sal(i,j,k)*out3ddat(5)
   ENDIF

ENDIF

!
!2. Diffusive flux in Z-direction
!--------------------------------
!

IF (runtitle(8:8).GT.'C') THEN

!      
!2.1 Isolevel diffusion
!----------------------
!
      
   IF (iopt_hdif_scal.EQ.2) THEN

      IF (k.GT.1) THEN
         
!        ---temperature
         xdiff = 0.0
         l_211: DO l=1,2
            IF (i+l.GT.2) THEN
               xdiff(l) = 0.5*(temp(i+l-1,j,k-1)+temp(i+l-1,j,k)-&
                             & temp(i+l-2,j,k-1)-temp(i+l-2,j,k))&
                             & /delxatu(i+l-1,j)
            ENDIF
         ENDDO l_211
         out3ddat(7) = -vdifcoefscal_rot(i,j,k)*(temp(i,j,k)-temp(i,j,k-1))&
                     & /delzatw(i,j,k)
         out3ddat(7) = out3ddat(7) +  &
                     & hdifcoef3datw(i,j,k)*xslopeatw_geo(i,j,k)*&
                     & 0.5*(xdiff(1)+xdiff(2))

!        ---salinity
         xdiff = 0.0
         l_212: DO l=1,2
            IF (i+l.GT.2) THEN
               xdiff(l) = 0.5*(sal(i+l-1,j,k-1)+sal(i+l-1,j,k)-&
                             & sal(i+l-2,j,k-1)-sal(i+l-2,j,k))&
                             & /delxatu(i+l-1,j)
            ENDIF
         ENDDO l_212
         out3ddat(8) = -vdifcoefscal_rot(i,j,k)*(sal(i,j,k)-sal(i,j,k-1))&
                     & /delzatw(i,j,k)
         out3ddat(8) = out3ddat(8) +  &
                     & hdifcoef3datw(i,j,k)*xslopeatw_geo(i,j,k)*&
                     & 0.5*(xdiff(1)+xdiff(2))
      ENDIF
         
!     ---density
      out3ddat(9) = -beta_temp(i,j,k)*out3ddat(7)+beta_sal(i,j,k)*out3ddat(8)
      
!      
!2.2 Isoneutral diffusion
!------------------------
!

   ELSEIF (iopt_hdif_scal.EQ.3) THEN

      IF (k.GT.1) THEN
         
         kp_220: DO kp=0,1
         ip_220: DO ip=0,1

            iip = i+1-ip
            kkp = k+kp-1

            slopefac = -hdifcoef3datu(iip,j,k)*xslopeatu_siso(i,j,kkp,1-ip,kp)
            out3ddat(7) = out3ddat(7) + slopefac*&
                        & (temp(iip,j,kkp)-temp(iip-1,j,kkp))/delxatu(iip,j)
            out3ddat(8) = out3ddat(8) + slopefac*&
                        & (sal(iip,j,kkp)-sal(iip-1,j,kkp))/delxatu(iip,j)

         ENDDO ip_220
         ENDDO kp_220

!        ---temperature
         out3ddat(7) = out3ddat(7) - vdifcoefscal_rot(i,j,k)*&
                     & (temp(i,j,k)-temp(i,j,k-1))/delzatw(i,1,k)
!        ---salinity
         out3ddat(8) = out3ddat(8) - vdifcoefscal_rot(i,j,k)*&
                     & (sal(i,j,k)-sal(i,j,k-1))/delzatw(i,1,k)
         !        ---density
         out3ddat(9) = -beta_temp(i,j,k)*out3ddat(7)+beta_sal(i,j,k)*out3ddat(8)

      ENDIF
         
   ENDIF

ENDIF
   
RETURN

END SUBROUTINE usrdef_tsr3d_vals
