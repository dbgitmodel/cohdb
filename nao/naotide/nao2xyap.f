c$main
      program nao2xyap
c
c Convert NAO cmp format to standard Longitude, Latitude,
c Amplitude, Phase format
c
c Code : Koji Matsumoto
c Date : 1999.06.08
c
      implicit double precision (a-h,o-z)
c
      parameter (mmax   = 721, nmax   = 541)
      parameter (iunt05 = 5)
c
      dimension amp(mmax,nmax), phs(mmax,nmax)
c     
      call rdmap(iunt05, amp   , phs   , dx    , dy    ,
     +           mmax  , nmax  , xmin  , ymax  , mend  ,
     +           nend  , undef                          )
c
      do n = 1,nend
         y = ymax - dfloat(n-1)*dy
         do m = 1,mend
            x = xmin + dfloat(m-1)*dx
c
            if (amp(m,n).lt.undef) then
               write(6,'(4f9.4)')x,y,amp(m,n),phs(m,n)
            endif
c
         enddo
      enddo
c
 19   stop
      end
c
c$rdcmp
c----------------------------------------------------------------
      subroutine rdcmp(iu    , amp   , phs   , mend  , nend  ,
     +                 mmax  , nmax  , fmt   , iamp  , iphs  ,
     +                 aunit , punit  )
c----------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      parameter (kc = 10)
c
      dimension amp(mmax,nmax), phs(mmax,nmax)
      dimension iamp(mmax)    , iphs(mmax)
      character fmt*6
c
c -----< Reading loop >-----
c
      kend = mend/kc
      krem = mod(mend,kc)
c
      do n = 1,nend
c
         do k = 1,kend
c
            m1 = (k-1)*kc + 1
            m2 = m1 + kc - 1
            read(iu,fmt) (iamp(m),m=m1,m2)
c
         enddo
c
         if (krem.ne.0) then
c
            m1 = kend*kc + 1
            m2 = kend*kc + krem
            read(iu,fmt) (iamp(m),m=m1,m2)
c
         endif
c
         do k = 1,kend
c
            m1 = (k-1)*kc + 1
            m2 = m1 + kc - 1
            read(iu,fmt) (iphs(m),m=m1,m2)
c
         enddo
c
         if (krem.ne.0) then
c
            m1 = kend*kc + 1
            m2 = kend*kc + krem
            read(iu,fmt) (iphs(m),m=m1,m2)
c
         endif
c
         do m = 1,mend
c
            amp(m,n) = dfloat(iamp(m))*aunit  ! in centimeters
            phs(m,n) = dfloat(iphs(m))*punit  ! in degrees
c
         enddo
c
      enddo
c
      close(iu)
c
      return
      end
c
c$rdhead
c----------------------------------------------------------------
      subroutine rdhead(iu    , name  , wave  , date  , xmin  ,
     +                  xmax  , ymin  , ymax  , dx    , dy    ,
     +                  mend  , nend  , ideff , fmt   , aunit ,
     +                  punit                                  )
c----------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      character fmt*6, wave*3, name*50, date*50, buf*50
c
      read(iu,'(13x,a50)') name
      read(iu,'(13x,a3 )') wave
      read(iu,'(13x,f5.3,19x,f4.2)') aunit,punit
      read(iu,'(13x,a50)') date
      read(iu,'(7x,f7.2,3(9x,f7.2))') xmin,xmax,ymin,ymax
      read(iu,'(12x,i2,14x,i2,2(9x,i7))') idx,idy,mend,nend
      read(iu,'(16x,i6,11x,a6)') ideff, fmt
c
      if (idx.eq.50) then
         dx = 0.5d0
         dy = 0.5d0
      else
         dx = 1.d0/dfloat(idx)
         dy = 1.d0/dfloat(idy)
      endif
c
      return
      end
c
c$rdmap
c------------------------------------------------------------
      subroutine rdmap(iu    , amp   , phs   , dx    , dy    ,
     +                 mmax  , nmax  , xmin  , ymax  , mend  ,
     +                 nend  , undef                          )
c------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      dimension amp(mmax,nmax), phs(mmax,nmax)
      dimension iwa(mmax)     , iwp(mmax)
c
      character fmt*6, wave*3, name*50, date*50
c
      call rdhead(iu    , name  , wave  , date  , xmin  ,
     +            xmax  , ymin  , ymax  , dx    , dy    ,
     +            mend  , nend  , ideff , fmt   , aunit ,
     +            punit                                  )
c
c      write(6,6002) name
c 6002 format('   Model name   = ',a50)
c      write(6,6003) date
c 6003 format('   Created date = ',a50)
c      write(6,101) xmin, xmax, ymin, ymax
c      write(6,102) dx, dy, mend, nend
c      write(6,103) aunit, punit
c 101  format(4x,'xmin = ',f6.2,', xmax = ',f6.2,', ymin = ',f6.2,
c     +          ', ymax = ',f6.2)
c 102  format(4x,'dx   = ',f6.3,', dy =   ',f6.3,', mend = ',i4,
c     +          '  , nend = ',i4)
c 103  format(4x,'amplitude unit = ',f6.3,', phase unit =  ',f6.3)
c
      undef = dfloat(ideff)*aunit
c
      call rdcmp(iu    , amp   , phs   , mend  , nend  ,
     +           mmax  , nmax  , fmt   , iwa   , iwp   ,
     +           aunit , punit  )
c
      close (iu)
c
      return
c
 98   print*,'Error : Can not open file.'
      stop
c
 99   return
      end
c
c ----------------------< End of program >----------------------
c \(^_^)/ \(^o^)/ \(^_^)/ \(^o^)/ \(^_^)/ \(^o^)/ \(^_^)/ \(^o^)/
