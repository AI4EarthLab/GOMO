#include "common.h"

subroutine bcond2_ua()
  use openarray
  use variables
  use config
  implicit none
  type(array):: tmpua,tmpel,tmp
  integer :: ierr
  
  tmpua=mat_zeros_im_jm_1
  tmpel=mat_zeros_im_jm_1
  tmp  =mat_zeros_im_jm_1

  call set(sub(tmpua, imm1,':'), uabe)
  call set(sub(tmpua, 2   ,':'), uabw)
  call set(sub(tmpel, imm1,':'), ele)    
  call set(sub(tmpel, 2   ,':'), elw)
 
  tmp = ramp * (tmpua + rfe * sqrt(grav/h) * (el - tmpel))
  call set(sub(uaf, im, [2,jmm1]), sub(tmp, imm1, [2,jmm1]))

  tmp = ramp * (tmpua - rfw * sqrt(grav/h) * (el - tmpel))
  call set(sub(uaf, 2, [2,jmm1]), sub(tmp, 2, [2,jmm1]))
  call set(sub(uaf, 1, [2,jmm1]), sub(tmp, 2, [2,jmm1]))
  call set(sub(uaf, [1,im], jm), 0.d0)
  call set(sub(uaf, [1,im], 1 ), 0.d0)
  uaf = uaf * dum

end subroutine bcond2_ua
