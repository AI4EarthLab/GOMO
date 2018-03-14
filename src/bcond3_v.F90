#include "common.h"

subroutine bcond3_v()
  use openarray
  use variables
  use config
  
  implicit none
  type(array):: tmph,tmp,tmpv

  tmph=h 
  call set(sub(tmph,  2,':',':'), sub(tmph,    1,':',':'))
  tmpv=v
  call set(sub(tmpv,':', jm,':'), sub(tmpv,':', jmm1,':'))
  call set(sub(tmpv,':', 2 ,':'), sub(tmpv,':', 3   ,':'))  

  tmp =sqrt(tmph/hmax)* AXF(AXB(tmpv)) + (1.d0 -sqrt(tmph/hmax))* AXF(AXB(v))
  call set(sub(vf, [2,imm1],jm,[1,kbm1]), sub(tmp, [2,imm1],jm,[1,kbm1]))
  call set(sub(vf, [2,imm1],2 ,[1,kbm1]), sub(tmp, [2,imm1],2 ,[1,kbm1]))
  call set(sub(vf, [2,imm1],1 ,[1,kbm1]), sub(tmp, [2,imm1],2 ,[1,kbm1]))  
  call set(sub(vf, im ,[1,jm], [1,kbm1]), 0.d0 )
  call set(sub(vf, 1  ,[1,jm], [1,kbm1]), 0.d0 )

  vf = vf * dvm

end subroutine
