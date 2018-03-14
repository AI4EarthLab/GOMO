#include "common.h"

subroutine bottom_friction()
  use openarray
  use variables
  use config
  implicit none
  real(kind=8) :: zz_kbm1
  
  !zz_kbm1 = sub(zz, 1,1,kbm1)
  zz_kbm1 = zz1(kbm1)
  cbc=(kappa/log((1.0+zz_kbm1)*h/z0b))*(kappa/log((1.0+zz_kbm1)*h/z0b));

  call set(cbc, cbcmin, cbc < cbcmin)
  call set(cbc, cbcmax, cbc > cbcmax)

end subroutine
