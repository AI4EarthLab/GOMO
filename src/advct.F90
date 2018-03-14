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

  advx  = DXB(AXF(AXB(dt) * u) * AXF(u) - dt*aam*2.d0*DXF(ub)) &
       + DYF(AXB(AYB(dt) * v) * AYB(u) - AYB(AXB(dt)) &
       * AYB(AXB(aam))*(DYB(ub) + DXB(vb))) &
       - AXB(curv * dt * AYF(v));

  advy  = DXF(AYB(AXB(dt) * u) * AXB(v) &
       - AYB(AXB(dt))*AYB(AXB(aam))*(DYB(ub) + DXB(vb))) &
       + DYB(AYF(AYB(dt) * v) * AYF(v) - dt * aam * 2.d0 *DYF(vb)) &
       + AYB(curv * dt * AXF(u));

end subroutine
