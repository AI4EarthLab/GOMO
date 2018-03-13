#include "common.h"

subroutine bcond2_va()
  use openarray
  use variables
  use config
  implicit none
  type(array):: tmpva,tmpel,tmp
  integer :: ierr
  
  tmpva=mat_zeros_im_jm_1
  tmpel=mat_zeros_im_jm_1
  tmp=  mat_zeros_im_jm_1
    
  call set(sub(tmpva,':', jmm1), vabn)
  call set(sub(tmpva,':', 2), vabs)
  call set(sub(tmpel,':', jmm1), eln)    
  call set(sub(tmpel,':', 2), els)

  tmp = ramp * (tmpva + rfn * sqrt(grav/h) * (el - tmpel))
  call set(sub(vaf,[2,imm1], jm), sub(tmp,[2,imm1], jmm1))

  tmp = ramp * (tmpva - rfs * sqrt(grav/h) * (el - tmpel))
  call set(sub(vaf,[2,imm1], 2), sub(tmp,[2,imm1], 2))
  call set(sub(vaf,[2,imm1], 1), sub(tmp,[2,imm1], 2))  
  call set(sub(vaf, im, [1,jm]), 0.d0)
  call set(sub(vaf, 1,  [1,jm]), 0.d0)
  
  vaf = vaf * dvm

end subroutine 
