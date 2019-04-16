#include "common.h"
subroutine internal_update()
  use openarray
  use config
  use variables
  implicit none

  egb = egf
  etb = et
  et = etf
  dt = h + et

  dhb = etb + h
  utb = utf
  vtb = vtf
  vfluxb = vfluxf
  
end subroutine internal_update
