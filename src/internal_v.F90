#include "macro.h"
subroutine internal_v()
  use openarray
  use config
  use variables
  implicit none
  type(array) dh_3d, tmpa,tmpdzz_3d
  integer :: ierr, k, i

  dh=AYB(dhf)

  vf =(AYB(dhb)*vb-dti2*(advy+ &
       drhoy+AYB(cor*dt*AXF(u))+     &
       grav * AYB(dt)*(DYB(egf+egb)+ &
       DYB(e_atmos)*2.d0)/2.d0 &
       -DZF(AYB(w) * AZB(v))))/dh
  
  call set(VF(kb), 0.d0)

  call set(sub(dh,1,':',':'),1.d0)
  call set(sub(dh,':',1,':'),1.d0)

  c=AYB(km)
  tmpa=mat_zeros
  call set(sub(tmpa,':',':',[1,kbm2]),-dti2*( sub(c,':',':',[2,kbm1])+ umol))
  a=tmpa/(dz*dzz*dh*dh) 
  call set(A(kbm1), 0.d0)
  call set(A(kb), 0.d0)
  
  tmpdzz=dzz
  call set(sub(tmpdzz,':',':',[2,kbm1]),sub(dzz,':',':',[1,kbm2]))
  c=-dti2*(c+umol)/(dz*tmpdzz*dh*dh) 
  call set(C(1), 0.d0)
  call set(C(kb), 0.d0)

  call set( EE(1), A(1)/( A(1)-1.d0 ) )
  call set( GG(1), (dti2*wvsurf/( DZ(1)*DH(1))-VF(1))/(A(1)-1.d0)  )

  do k=2,kbm2
    call set(GG(k), 1.d0 / (A(k) + C(k) * (1.d0 - EE(k-1)) -1.d0))
    call set(EE(k), A(k) * GG(k))
    call set(GG(k), ( C(k) * GG(k-1) - VF(k) ) * GG(k)) 
  enddo

  tps=AYB(cbc) * sqrt( AYB( AXF( UB(kbm1) ) )**2 + VB(kbm1)**2 )
  call set(sub(tps,1, ':',':'), 0.d0)
  call set(sub(tps,im,':',':'), 0.d0)

  call set(VF(kbm1) , ( C(kbm1)* GG(kbm2)-VF(kbm1) )/ &
       (tps*dti2 /(-DZ(kbm1)*sub(dh,':',':',kbm1))-1.d0-C(kbm1)*( EE(kbm2)-1.d0 ) ) ) 
 
  do k=kbm2,1,-1
    call set(VF(k) , EE(k)*VF(k+1)+GG(k) )
  enddo

  vf=vf*dvm
  wvbot=-tps * VF(kbm1)

  call bcond3_v()

  print*,"vf after bcond"
  call disp(vf,'vf = ')
  
end subroutine
