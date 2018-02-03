#include "common.h"

subroutine advct()
  use openarray
  use config  
  use variables
  implicit none
  integer:: ierr, i
  type(array) :: tmp
  
  curv = (AYF(v) * DXB(AXF(dy))  &
       - AXF(u) * DYB(AYF(dx))) / dx / dy
  ! the equation has two fusion kernels, split to two in order to have same
  ! answer with matlab

  !advx  = DXB(AXF(AXB(dt) * u) * AXF(u) - dt*aam*2.d0*DXF(ub)) &
  !     + DYF(AXB(AYB(dt) * v) * AYB(u) - AYB(AXB(dt)) &
  !     * AYB(AXB(aam))*(DYB(ub) + DXB(vb))) &
  !     - AXB(curv * dt * AYF(v));

  ! call open_debug()

  tmp  = AXF(AXB(dt) * u) * AXF(u) - dt*aam*2.d0*DXF(ub)
  ! call close_debug()

  call set(sub(tmp,1,':',':'), 0.d0)
  call set(sub(tmp,im,':',':'), 0.d0)

  advx  = DXB(tmp) &
       + DYF(AXB(AYB(dt) * v) * AYB(u) - AYB(AXB(dt)) &
       * AYB(AXB(aam))*(DYB(ub) + DXB(vb))) &
       - AXB(curv * dt * AYF(v));

  call set(sub(advx,1,':',':'), 0.d0)
  call set(sub(advx,':',1,':'), 0.d0)
  call set(sub(advx,':',jm,':'), 0.d0)

  !call disp(advx, "advx = ")
  

  ! same as advx, advy should split into two equations

  !advy  = DXF(AYB(AXB(dt) * u) * AXB(v) &
  !     - AYB(AXB(dt))*AYB(AXB(aam))*(DYB(ub) + DXB(vb))) &
  !     + DYB(AYF(AYB(dt) * v) * AYF(v) - dt * aam * 2.d0 *DYF(vb)) &
  !     + AYB(curv * dt * AXF(u));

  ! call open_debug()
  tmp = AYF(AYB(dt) * v) * AYF(v) - dt*aam*2.d0*DYF(vb)
!  call close_debug()
  call set(sub(tmp,':',1,':'), 0.d0)
  call set(sub(tmp,':',jm,':'), 0.d0)
  
  advy  = DXF(AYB(AXB(dt) * u) * AXB(v) &
       - AYB(AXB(dt))*AYB(AXB(aam))*(DYB(ub) + DXB(vb))) &
       + DYB(tmp) &
       + AYB(curv * dt * AXF(u));

  call set(sub(advy,':',1,':'), 0.d0)
  call set(sub(advy,1,':',':'), 0.d0)
  call set(sub(advy,im,':',':'), 0.d0)
  
  !call disp(advy, "advy = ")

end subroutine
