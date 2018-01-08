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
  print*, "xxxxxxxxx"
  call tic("aam")

  aam=horcon * dx * dy &
       * sqrt(DXF(u)**2 + DYF(v)**2 &
       + 0.5*(DYB(AYF(AXF(u))) + DXB(AXF(AYF(v))))**2);

  ! do i = 1, kb
  !    AAM(i) = horcon * DX_2D * DY_2D * AAM(i)
  ! end do
  call toc("aam")
  
  ! aam_expr=horcon * dx_3d * dy_3d &
  !      *sqrt(DXF(u)**2 + DYF(v)**2 &
  !      +0.5*(DYB(AYF(axf_u)) + DXB(AXF(ayf_v)) )**2);
  ! call gen_hash(aam_expr)
  ! call write_graph(aam_expr, file="aam_expr.dot")
  
  call set(sub(aam, ':',1,':'),  aam_init)
  call set(sub(aam, ':',jm,':'), aam_init)
  call set(sub(aam, 1,':',':'),  aam_init)
  call set(sub(aam, im,':',':'), aam_init)
  call set(sub(aam, ':',':',kb), aam_init)

end subroutine
