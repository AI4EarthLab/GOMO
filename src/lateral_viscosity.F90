#include "common.h"

subroutine lateral_viscosity()
  use openarray
  use config
  use variables
  implicit none
  type(node) :: aam_expr
  
  call tic("advct")
  call advct()
  call toc("advct")

  if(npg == 1) then
     call tic("baropg")
     call baropg()
     call toc("baropg")
  elseif(npg == 2) then
     call baropg_mcc()
  else
     print*, "Error: invalid value for npg"
     stop
  endif
  call tic("aam")
  

  ! pos of u and v has changed before, need check OpenArray
  call grid_bind(u, 2)
  call grid_bind(v, 1)

  aam=horcon * dx * dy &
       * sqrt(DXF(u)**2 + DYF(v)**2 &
       + 0.5*(DYB(AYF(AXF(u))) + DXB(AXF(AYF(v))))**2);

  call toc("aam")
  
  call set(sub(aam, ':',1,':'),  aam_init)
  call set(sub(aam, ':',jm,':'), aam_init)
  call set(sub(aam, 1,':',':'),  aam_init)
  call set(sub(aam, im,':',':'), aam_init)
  call set(sub(aam, ':',':',kb), aam_init)
  
  !call disp(aam, "aam = ")

end subroutine
