#include "common.h"

!function [cbc] = bottom_friction(kappa,zz,h,z0b,cbcmin,cbcmax)
subroutine bottom_friction()
  use openarray
  use variables
  use config
  implicit none
  real(kind=8) :: zz_kbm1
  !global kbm1;
  zz_kbm1 = sub(zz, 1,1,kbm1)
  !cbc=(kappa./log((1.0+zz(kbm1))*h/z0b)).^2;
  cbc=(kappa/log((1.0+zz_kbm1)*h/z0b))**2;

  ! cbc=max(cbcmin, cbc);
  ! cbc=min(cbcmax, cbc);

  call set(cbc, cbcmin, cbc < cbcmin)
  call set(cbc, cbcmax, cbc > cbcmax)
end subroutine
