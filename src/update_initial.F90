#include "common.h"

subroutine update_initial()
  use openarray
  use variables
  use config
  implicit none
  integer::k

  ua = uab;
  va = vab;
  el = elb;

  d  = h + el

  et = etb;
  etf= et;

  dt = h + et;

  ! do k=1,kb
  !   DT_3D(k) =DT_2D
  !   ETB_3D(k)=ETB_2D
  !   D_3D(k)  =D_2D
  !   E_ATMOS_3D(k)=E_ATMOS_2D
  ! enddo

  ! axbdt_3d=AXB(dt_3d)
  ! aybdt_3d=AYB(dt_3d)

  ! axbdt = AXB(dt)
  ! aybdt = AYB(dt)

  l = dt * 0.1d0

  l = make_psudo3d(l)

  call grid_bind(l, 7)

  q2b = small
  q2 = q2b
  q2lb = l * q2b
  q2l = q2lb
  kh = l * sqrt(q2b)
  km = kh
  kq = kh
  aam = aam_init

  call set(sub(w, ':',':',1), vfluxf)

  t = tb
  s = sb
  u = ub
  !axf_u=AXF(u)
  v = vb
  !ayf_v=AYF(v)

!  call dens(rho, t, s)

!  if(npg == 1) then
     ![drhox,drhoy]   = baropg(rho, rmean, dt_3d, ramp);
!     call baropg()
!  elseif(npg == 2) then
     ![drhox,drhoy]   = baropg_mcc(rho,rmean,d_3d,dt_3d,ramp);
!     call baropg_mcc()     
!  else
!     print*, "Error: invalid value for npg";
!     stop
!  end if

  ! drx2d=sum(drhox.*dz_3d, 3);
  ! dry2d=sum(drhoy.*dz_3d, 3);
  ! drx2d = sum(drhox * dz_3d, 3)
  ! dry2d = sum(drhoy * dz_3d, 3);
  !call disp(drx2d, 'dxr2d = ')
end subroutine
