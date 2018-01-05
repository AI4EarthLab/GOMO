#include "macro.h"

subroutine adjust_ufvf()
  use openarray
  use config
  use variables
  implicit none
  integer::k,ierr
  type(array)::fluc_u,fluc_v,sum_udz,sum_vdz

  fluc_u =uf+ub-2.d0*u
  sum_udz=sum(fluc_u * dz, 3)

  fluc_v =vf+vb-2.d0*v
  sum_vdz=sum(fluc_v * dz, 3)

  do k=1, kbm1
    U(k)=U(k)+0.5d0*smoth*(FLUC_U(k)-SUM_UDZ1)
    V(k)=V(k)+0.5d0*smoth*(FLUC_V(k)-SUM_VDZ1)
  enddo

!  u=u+0.5d0*smoth*(uf+ub-2.d0*u-rep(sum((uf+ub-2.d0*u) * dz, 3),1,1,kb))
!  call set(sub(u,':',':',kb),0.d0)
  U(kb)=0.d0
  ub = u
  u  = uf
  axf_u=AXF(u)

!  v=v+ 0.5d0*smoth*(vf+vb-2.d0*v-rep(sum((vf+vb-2.d0*v) * dz, 3),1,1,kb))
!  call set(sub(v,':',':',kb),0.d0)
  V(kb)=0.d0
  vb = v
  v  = vf
  ayf_v=AYF(v)

end subroutine
