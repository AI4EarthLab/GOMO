#include "macro.h"

subroutine internal_update()
  use openarray
  use config
  use variables
  implicit none
  integer::k

  egb = egf
  etb = et
  et  = etf
  dt  = h + et

  ! do k=1,kb
  !  DT_3D(k) =DT_2D
  !  EGB_3D(k)=EGB_2D
  !  ETB_3D(k)=ETB_2D
  ! enddo

  dhb = etb + h
  ! axbdt=AXB(dt);        aybdt=AYB(dt)
  utb = utf
  vtb = vtf
  vfluxb = vfluxf
  
end subroutine internal_update
