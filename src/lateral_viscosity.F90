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

  aam=horcon * dx * dy &
       * sqrt(DXF(u)*DXF(u) + DYF(v)*DYF(v) &
       + 0.5*(DYB(AYF(AXF(u))) + DXB(AXF(AYF(v))))*(DYB(AYF(AXF(u))) +   DXB(AXF(AYF(v)))));

  call toc("aam")
  
  call set(sub(aam, ':',1,':'),  aam_init)
  call set(sub(aam, ':',jm,':'), aam_init)
  call set(sub(aam, 1,':',':'),  aam_init)
  call set(sub(aam, im,':',':'), aam_init)
  call set(sub(aam, ':',':',kb), aam_init)

end subroutine
