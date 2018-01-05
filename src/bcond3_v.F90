#include "common.h"

subroutine bcond3_v()
  use openarray
  use variables
  use config
  
  implicit none
  type(array):: tmph,tmp,tmpv
  integer :: ierr
 !   tmph = h; tmph(:, 2) = tmph(:, 1); tmph = repmat(tmph, 1, 1, kb);
 !   tmpv = v; tmpv(:, jm, :) = tmpv(:, jmm1, :); tmpv(:, 2, :) = tmpv(:, 3, :);
  tmph=h_3d 
  call set(sub(tmph,  2,':',':'), sub(tmph,    1,':',':'))
  tmpv=v
  call set(sub(tmpv,':', jm,':'), sub(tmpv,':', jmm1,':'))
  call set(sub(tmpv,':', 2 ,':'), sub(tmpv,':', 3   ,':'))  
 
!    tmp = sqrt(tmph/hmax) .* AXF(AXB(tmpv)) + (1.0 - sqrt(tmph/hmax)) .* AXF(AXB(v));
!    vf(2:imm1, jm, 1:kbm1) = tmp(2:imm1, jm, 1:kbm1);
!    vf(2:imm1, 2, 1:kbm1) = tmp(2:imm1, 2, 1:kbm1);
!    vf(2:imm1, 1, 1:kbm1) = tmp(2:imm1, 2, 1:kbm1);
!    vf(im, 2:jmm1, 1:kbm1) = 0.e0;
!    vf(1, 2:jmm1, 1:kbm1) = 0.e0;
!    vf = vf .* dvm_3d;  


  tmp =sqrt(tmph/hmax)* AXF(AXB(tmpv)) + (1.d0 -sqrt(tmph/hmax))* AXF(AXB(v))
  call set(sub(vf, [2,imm1],jm,r(1,kbm1)), sub(tmp, [2,imm1],jm,r(1,kbm1)))
  call set(sub(vf, [2,imm1],2 ,r(1,kbm1)), sub(tmp, [2,imm1],2 ,r(1,kbm1)))
  call set(sub(vf, [2,imm1],1 ,r(1,kbm1)), sub(tmp, [2,imm1],2 ,r(1,kbm1)))  
  call set(sub(vf, im ,r(2,jmm1), r(1,kbm1)), 0.d0 )
  call set(sub(vf, 1  ,r(2,jmm1), r(1,kbm1)), 0.d0 )

  vf = vf * dvm_3d

  call destroy(tmpv,ierr); call destroy(tmph,ierr)
  call destroy(tmp,ierr)

end subroutine
