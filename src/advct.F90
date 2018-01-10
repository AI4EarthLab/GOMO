#include "common.h"

subroutine advct()
  use openarray
  use config  
  use variables
  implicit none
  integer:: ierr, i
  type(array) :: tmp
  
  ! curv= (ayf_v * DXB(AXF(dy_3d)) &
  !      - axf_u * DYB(AYF(dx_3d))) / (dx_3d * dy_3d);

  !! Warning, this expression will stuck...
  ! tmp = (AYF(v) * DXB(AXF(dy))  &
  !      - AXF(u) * DYB(AYF(dx))) / (dx * dy)
  
  ! call disp(tmp, 'tmp = ')


  curv = (AYF(v) * DXB(AXF(dy))  &
       - AXF(u) * DYB(AYF(dx))) / dx / dy

  ! do i = 1, kb
  !    ! CURV(i) = (AYF_V(i) * DXB_AXF_DY - AXF_U(i) * DYB_AYF_DX) &
  !    !      / (DX_2D * DY_2D)

  !    CURV_X(i) = (AYF_V(i) * DXB_AXF_DY - AXF_U(i) * DYB_AYF_DX) &
  !         * DT_3D(1) * AYF_V(1)/ (DX_2D * DY_2D)

  !    CURV_Y(i) = (AYF_V(i) * DXB_AXF_DY - AXF_U(i) * DYB_AYF_DX) &
  !         * DT_3D(1) * AXF_U(1)/ (DX_2D * DY_2D)

  ! end do

  ! call grid_bind(curv_x, gs, ayf_v%grid_pos)
  ! call grid_bind(curv_y, gs, axf_u%grid_pos)
  
  ! tmp1 = dt_3d*aam*2.d0
  ! tmp2 = axbdt_3d * u
  ! tmp3 = aybdt_3d * v
  ! tmp4 = AYB(axbdt_3d ) * AYB(AXB(aam))*(DYB(ub) + DXB(vb))

  ! advx  = DXB(AXF( tmp2 ) * axf_u - tmp1*DXF(ub)) &
  !       + DYF(AXB( tmp3 ) * AYB(u) - tmp4 ) &
  !       - AXB(curv * dt_3d * ayf_v );

  !dt = make_psudo3d(dt)

  advx  = DXB(AXF(AXB(dt) * u) * AXF(u) - dt*aam*2.d0*DXF(ub)) &
       + DYF(AXB(AYB(dt) * v) * AYB(u) - AYB(AXB(dt)) &
       * AYB(AXB(aam))*(DYB(ub) + DXB(vb))) &
       - AXB(curv * dt * AYF(v));
  
  call set(sub(advx,1,':',':'), 0.d0)

  ! advy  = DXF(AYB( tmp2 ) * AXB(v) - tmp4 ) &
  !       + DYB(AYF( tmp3 ) * ayf_v - tmp1*DYF(vb)) &
  !       + AYB(curv * dt_3d * axf_u );

  advy  = DXF(AYB(AXB(dt) * u) * AXB(v) &
       - AYB(AXB(dt))*AYB(AXB(aam))*(DYB(ub) + DXB(vb))) &
       + DYB(AYF(AYB(dt) * v) * AYF(v) - dt * aam * 2.d0 *DYF(vb)) &
       + AYB(curv * dt * AXF(u));
  
  call set(sub(advy,':',1,':'), 0.d0)

end subroutine
