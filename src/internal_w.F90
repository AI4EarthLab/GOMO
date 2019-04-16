#include "common.h"
subroutine internal_w()
  use openarray
  use config
  use variables
  implicit none
  type(array) :: del_w
 
  del_w = shift(csum(dz*(DXF(AXB(dt)*u)+DYF(AYB(dt) * v) + &
       (etf-etb)/dti2), 3), 0, 0, 1)

  w = (0.5d0*(vfluxb+vfluxf) + del_w) * fsm

  call set(sub(w, 1, ':',':'), 0.d0)
  call set(sub(w, im,':',':'), 0.d0)
  call set(sub(w, ':', ':', kb), 0.d0)

end subroutine
