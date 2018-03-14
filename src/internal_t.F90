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
  type(array) :: dh, tmp, tmp1
  integer :: ierr, i

  rr =(/0.58d0,0.62d0,0.67d0,0.77d0,0.78d0/)
  ad1=(/0.35d0,0.60d0,1.d0,  1.5d0, 1.4d0/)
  ad2=(/0.23d0,20.d0 ,17.d0, 14.d0, 7.9d0/)

  if(nadv==1) then
    call tic("internal_t_ff")

    tmp = AZB(f)*w
    tmp1 = f*w
    call set(sub(tmp, ':',':',1), sub(tmp1, ':',':',1))
    call set(sub(tmp, ':',':',kb), 0.d0)
    
    ff=( dhb*fb &
      - dti2*(DXF(AXB(dt)*AXB(f)*u &
      -AXB(aam)*AXB(h)*tprni*DXB(fb)*dum) + &
      DYF(AYB(dt)*AYB(f)*v-AYB(aam)*AYB(h)*tprni &
      *DYB(fb)*dvm)-DZF(tmp)))/dhf
    call set(sub(ff, ':',':',kb), 0.d0)
    call toc("internal_t_ff")     
  else
    ff=mat_zeros
    call advt2(ff,f,fb,fclim)
  endif
 
  rad=mat_zeros
  dh=h+etf       

  call set(sub(a,':',':',[1,kbm2]), -dti2*(sub(kh,':',':',[2,kbm1])+umol))
  a = a/(dz*dzz*dhf*dhf) 
  call set(A(kbm1), 0.d0)
  call set(A(kb), 0.d0)

  do i = 2, kbm1
    call set(C(i), dzz1(i-1))
  enddo
  c=-dti2*(kh+umol)/(dz*c*dhf*dhf)
  call set(C(1), 0.d0)
  call set(C(kb), 0.d0)

  if(nbc==2 .or. nbc==4) then
    rad = mat_ones*frad*(rr(ntp)*exp(z*dhf/ad1(ntp)) &
        +(1.d0-rr(ntp))*exp(z*dhf/ad2(ntp)))
    call set(sub(rad,':',':',kb) , 0.d0)
  end if

  if(nbc==1) then
    call set(EE(1), A(1)/(A(1)-1.d0))
    call set(GG(1), (dti2*wfsurf/(dz1(1)*dh)-FF(1))/(A(1)-1.d0) )
  elseif(nbc==2) then
    call set( EE(1), A(1)/( A(1)-1.d0 ))
    call set( GG(1), ( dti2*(wfsurf+RAD(1)-RAD(2))/( dz1(1)*dh )-FF(1) )/(A(1)-1.d0) )
  elseif(nbc==3 .or. nbc==4) then
    call set(EE(1), 0.d0)
    call set(GG(1), fsurf)
  endif
  
  do k=2,kbm1
    call set(GG(k), 1.d0 / (A(k)+C(k)*(1.d0-EE(k-1)) -1.d0))
    call set(EE(k), A(k) * GG(k))
    call set(GG(k), (C(k) * GG(k-1) - FF(k) &
         +dti2*(RAD(k)-RAD(k+1))/(dh*dz1(k))) * GG(k))
  enddo

  call set(FF(kbm1), (C(kbm1)*GG(kbm2)-FF(kbm1)+dti2*(RAD(kbm1)-RAD(kb)) &
          /(dh*dz1(kbm1))) / (C(kbm1)*(1.d0-EE(kbm2))-1.d0))

  do k=kbm2,1,-1
    call set(FF(k),EE(k)*FF(k+1)+GG(k))
  enddo

  call tic("bcond4")
  call bcond4(ff,f,fb,fbe,fbw,fbn,fbs)
  call toc("bcond4")
  
  call tic("smoth_update1")
  call smoth_update1(ff,f,fb)
  call toc("smoth_update1")
  
end subroutine
