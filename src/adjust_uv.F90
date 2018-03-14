#include "common.h"

subroutine adjust_uv()
  use openarray
  use variables
  use config
  use openarray
  implicit none

  u = u-sum(u*dz, 3) + (utb+utf)/(2.d0 * AXB(dt))
  call set(sub(u,':',':',kb) , 0.d0)

  v = v-sum(v*dz, 3) + (vtb+vtf)/(2.d0 * AYB(dt))
  call set(sub(v,':',':',kb) , 0.d0)

end subroutine
