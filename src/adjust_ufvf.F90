#include "common.h"
subroutine adjust_ufvf()
  use openarray
  use config
  use variables
  implicit none
  integer::k,ierr
  type(array)::fluc_u,fluc_v,sum_udz,sum_vdz

  u=u+0.5d0*smoth*(uf+ub-2.d0*u-sum((uf+ub-2.d0*u) * dz, 3))
  call set(sub(u,':',':',kb),0.d0)

  v=v+ 0.5d0*smoth*(vf+vb-2.d0*v-sum((vf+vb-2.d0*v) * dz, 3))
  call set(sub(v,':',':',kb),0.d0)

  ub = u
  u  = uf

  vb = v
  v  = vf

end subroutine
