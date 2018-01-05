#include "common.h"

subroutine adjust_uv()
  use openarray
  use variables
  use config
  use openarray
  implicit none
  integer::k,ierr 

  ! sum_udz=sum(u*dz,3)
  ! sum_vdz=sum(v*dz,3)
  ! ave_ut=0.5d0*(utb+utf)/AXB(dt)
  ! ave_vt=0.5d0*(vtb+vtf)/AYB(dt)

  ! do k=1,kbm1
  !    U(k)=U(k)-SUM_UDZ1+AVE_UT
  !    V(k)=V(k)-SUM_VDZ1+AVE_VT 
  ! enddo
  
  ! U(kb)=0.d0
  ! V(kb)=0.d0
  
  call set(u, u-sum(u*dz, 3) + (utb+utf)/(2.d0 * AXB(dt)))
  call set(sub(u,':',':',kb) , 0.d0)
  !axf_u=AXF(u)

  call set(v, v-sum(v*dz, 3) + (vtb+vtf)/(2.d0 * AYB(dt)))
  call set(sub(v,':',':',kb) , 0.d0)
  !ayf_v=AYF(v)

end subroutine
