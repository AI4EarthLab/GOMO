#include "common.h"
subroutine bcond3_u()
  use openarray
  use variables
  use config
  implicit none
  type(array):: tmpu,tmph,tmp

  tmph=h 
  call set(sub(tmph,  2,':',':'), sub(tmph,    1,':',':'))
  tmpu=u
  call set(sub(tmpu, im,':',':'), sub(tmpu, imm1,':',':'))
  call set(sub(tmpu, 2 ,':',':'), sub(tmpu, 3   ,':',':'))  
 
  tmp = sqrt(tmph/hmax)*AYF(AYB(tmpu)) + (1.d0-sqrt(tmph/hmax))*AYF(AYB(u))
  call set(sub(uf, im,[2,jmm1],[1,kbm1]), sub(tmp, im,[2,jmm1],[1,kbm1]) )
  call set(sub(uf, 2 ,[2,jmm1],[1,kbm1]), sub(tmp, 2 ,[2,jmm1],[1,kbm1]) )
  call set(sub(uf, 1 ,[2,jmm1],[1,kbm1]), sub(tmp, 2 ,[2,jmm1],[1,kbm1]) )
  call set(sub(uf, [1,im], jm ,[1,kbm1]), 0.d0 )
  call set(sub(uf, [1,im], 1  ,[1,kbm1]), 0.d0 )

  uf = uf * dum
end subroutine
