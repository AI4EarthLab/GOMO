#include "common.h"
#include "macro.h"

subroutine internal_t(ff,f,fb,wfsurf,fsurf,nbc,frad,fclim,fbe,fbw,fbn,fbs)
  use openarray
  use variables
  use config
  implicit none
  integer, intent(in) :: nbc
  integer :: k
  real*8 rr(5),ad1(5),ad2(5)
  type(array), intent(inout):: ff,f,fb,wfsurf,fsurf,frad,fclim,fbe,fbw,fbn,fbs
  type(array) :: dh
  integer :: ierr, i

  rr =(/0.58d0,0.62d0,0.67d0,0.77d0,0.78d0/)
  ad1=(/0.35d0,0.60d0,1.d0,  1.5d0, 1.4d0/)
  ad2=(/0.23d0,20.d0 ,17.d0, 14.d0, 7.9d0/)

  if(nadv==1) then
     call tic("internal_t_ff")
     ff=( fb*dhb_3d &
          - dti2*(DXF(axbdt_3d *AXB(f)*u &
          -AXB(aam)*AXB(h_3d)*tprni*DXB(fb)*dum_3d) + &
          DYF(aybdt_3d *AYB(f)*v-AYB(aam)*AYB(h_3d)*tprni &
          *DYB(fb)*dvm_3d)-DZF(AZB(f)*w)))/dhf_3d
     call toc("internal_t_ff")     
  else
     ff=mat_zeros
     call advt2(ff,f,fb,fclim)
  endif
 
  rad=mat_zeros
  dh=h+etf       ! ; dh_3d=h_3d+etf_3d

  !call set(sub(a,':',':',r(1,kbm2)) , -dti2*( sub(kh,':',':',r(2,kbm1))+umol ) )
  !a=a/(dz_3d*dzz_3d*dhf_3d*dhf_3d)              !a(kbm1:kb)=0

  do i = 1, kbm2
     A(i) = -dti2 * (KH(i+1) + umol) / (dz1(i) * dzz1(i) * DHF_3D(i)**2)
  end do
  
  A(kbm1)= 0.d0
  A(kb)= 0.d0

  !call set(sub(c,':',':',r(2,kbm1)) , sub(dzz_3d,':',':',r(1,kbm2)))
  do i = 2, kbm1
     C(i) = dzz1(i-1)
  enddo
  
  !c=-dti2*(kh+umol)/(dz_3d*c*dhf_3d*dhf_3d)
  do i = 2, kbm1
     C(i) = -dti2*(KH(i) + umol) / (dz1(i) * C(i) * DHF_3D(i)**2)
  end do
  
  C(1)=0.d0
  C(kb)=0.d0

  if(nbc==2 .or. nbc==4) then
     do k=1,kb 
        RAD(k)=FRAD
     enddo
     rad= rad*(rr(ntp)*exp(z_3d*dhf_3d/ad1(ntp)) &
          +(1.d0-rr(ntp))*exp(z_3d*dhf_3d/ad2(ntp)))
     call set(sub(rad,':',':',kb) , 0.d0)
  end if

  if(nbc==1) then
!     call set( EE(1), A(1)/( A(1)-1.d0 ))
!     call set( GG(1), ( dti2*wfsurf/( DZ(1)*dh )-FF(1) )/(A(1)-1.d0) )
     EE(1)=A(1)/( A(1)-1.d0 )
     GG(1)=( dti2*WFSURF_2D/( DZ(1)*DH_2D )-FF(1) )/(A(1)-1.d0)
  elseif(nbc==2) then
!     call set( EE(1), A(1)/( A(1)-1.d0 ))
!     call set( GG(1), ( dti2*(wfsurf+RAD(1)-RAD(2))/( DZ(1)*dh )-FF(1) )/(A(1)-1.d0) )
      EE(1)= A(1)/( A(1)-1.d0 )
      GG(1)= (dti2*(WFSURF_2D+RAD(1)-RAD(2))/( DZ(1)*DH_2D )-FF(1) )/(A(1)-1.d0) 
  elseif(nbc==3 .or. nbc==4) then
!     call set(EE(1), 0.d0)
!     call set(GG(1), fsurf)
     EE(1)= 0.d0
     GG(1)= FSURF_2D
  endif
  
  do k=2,kbm1
!     call set(GG(k), 1.d0 / (A(k)+C(k)*(1.d0-EE(k-1)) -1.d0))
!     call set(EE(k), A(k) * GG(k))
!     call set(GG(k), (C(k) * GG(k-1) - FF(k) &
!          +dti2*(RAD(k)-RAD(k+1))/(DZ(k)*dh)) * GG(k))
     GG(k)= 1.d0 / (A(k)+C(k)*(1.d0-EE(k-1)) -1.d0)
     EE(k)= A(k) * GG(k)
     GG(k)= (C(k) * GG(k-1) - FF(k) &
          +dti2*(RAD(k)-RAD(k+1))/(DZ(k)*DH_2D)) * GG(k)
  enddo

!  call disp(ee,'internal_t ee=')
!  call disp(gg,'internal_t gg=')


  !call disp_info(f, 'f_'//trim(i2s(__LINE__)))
  
  do k=kbm1,1,-1
  !   call set(FF(k),EE(k)*FF(k+1)+GG(k))
     FF(k)=EE(k)*FF(k+1)+GG(k)
  enddo

!  call disp(ff, 'ff = after internal_t')

  call tic("bcond4")
  call bcond4(ff,f,fb,fbe,fbw,fbn,fbs)
  call toc("bcond4")
  
!  call disp(ff,'ff = after bcond4')
  !call disp_info(f, 'f_'//trim(i2s(__LINE__)))

  call tic("smoth_update1")
  call smoth_update1(ff,f,fb)
  call toc("smoth_update1")
  
!  call disp(ff,'ff')
!  call disp(f,'f')
!  call disp(fb,'fb')
!  print*,"end of internal_t"

  call destroy(dh, ierr)
!  call destroy(dh_3d, ierr)
  
end subroutine
