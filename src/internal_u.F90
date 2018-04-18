#include "common.h"
subroutine internal_u()
  use openarray
  use config
  use variables
  implicit none
  type(array) :: dh, tmpa, bondu, tmparray, tmp
  integer :: ierr, k, i, pos(3), imax, jmax
  real(kind=8) :: tmpmax

  dh=AXB(dhf)

  tmp = AXB(w) * AZB(u)
  call set(sub(tmp, ':',':',kb), 0.d0)

  uf=(AXB(dhb)*ub-dti2*(advx + &
       drhox - AXB(cor*dt*AYF(v) )+  &
       grav* AXB(dt)*(DXB(egf+egb)+ &
       DXB(e_atmos)*2.d0)*0.5d0-  &
       DZF(tmp))) /dh
  
  call set(sub(dh,1,':',':'),1.d0)
  call set(sub(dh,':',1,':'),1.d0)

  c=AXB(km)
  call set(sub(c,1,':',':'),0.d0)
  tmpa=mat_zeros
  call set(sub(tmpa,':',':',[1,kbm2]),-dti2*( sub(c,':',':',[2,kbm1])+ umol))
  a=tmpa/(dz*dzz*dh*dh)
  call set(A(kbm1), 0.d0)
  call set(A(kb), 0.d0)

  !bondu = AXB(w) * AZB(u)
  !call set(sub(uf,im,':',':') , sub(bondu,im,':',':')) 

  c=-dti2*(c+umol)/(dz*tmpdzz*dh*dh)
  
  call set(C(1) ,0)
  call set(C(kb) ,0)
  
  call set(EE(1),A(1)/(A(1)-1.d0) )
  call set(GG(1),(dti2*wusurf/(dz1(1)*dh)-UF(1))/(A(1)-1.d0))

  do k=2,kbm2
    call set(GG(k), 1.d0 / (A(k) + C(k) * (1.d0 - EE(k-1)) -1.d0))
    call set(EE(k), A(k) * GG(k))
    call set(GG(k), (C(k) * GG(k-1) - UF(k))* GG(k))
  enddo
 
  tps=AXB(cbc) * sqrt( UB(kbm1)*UB(kbm1) + AXB( AYF( VB(kbm1)))*AXB( AYF( VB(kbm1))))

  
 call set(UF(kbm1) , (C(kbm1)*GG(kbm2)-UF(kbm1))/ &
      (tps*dti2 /(-dz1(kbm1)*dh)-&
      1.d0-C(kbm1)*( EE(kbm2)-1.d0))) 
  
  do k=kbm2,1,-1
     call set(UF(k) , EE(k)*UF(k+1)+GG(k))
  end do

  uf=uf*dum;
  wubot=-tps*UF(kbm1)

  call bcond3_u()

end subroutine
