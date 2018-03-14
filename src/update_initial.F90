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
  et = etb;
  etf = et;

  d = h + el;
  dt = h + et;

  l = dt * 0.1d0

  call grid_bind(l, 7)

  q2b = small
  q2lb = l * q2b
  kh = l * sqrt(q2b)
  km = kh
  kq = kh
  aam = aam_init

  call set(sub(w, ':',':',1), vfluxf)

  q2 = q2b
  q2l = q2lb
  t = tb
  s = sb
  u = ub
  v = vb

  call dens(rho, t, s)

  if(npg == 1) then
    call baropg()
  elseif(npg == 2) then
    call baropg_mcc()     
  else
     print*, "Error: invalid value for npg";
     stop
  end if

  drx2d = sum(drhox * dz, 3)
  dry2d = sum(drhoy * dz, 3)

end subroutine
