#include "macro.h"

subroutine internal_w()
  use openarray
  use config
  use variables
  implicit none
  type(array) :: del_w,vflux
  integer :: ierr,k
 
  del_w = shift(csum(dz*(DXF(AXB(dt)*u)+DYF(AYB(dt) * v) + &
       (etf-etb)/dti2), 3), 0, 0, 1)
  
  !vflux=0.5d0*(vfluxb+vfluxf)

  w = (0.5d0*(vfluxb+vfluxf) + del_w) * fsm
  
  ! call toc(2)
  ! do k=1,kbm1
  !   W(k)=(DEL_W(k)+VFLUX)*FSM_3D(k)
  ! enddo

  call set(sub(w, 1, ':',':'), 0.d0)
  call set(sub(w, im,':',':'), 0.d0)
  call set(sub(w, ':', ':', kb), 0.d0)

  !W(kb)=0.d0
end subroutine
