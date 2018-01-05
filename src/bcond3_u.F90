#include "common.h"
subroutine bcond3_u()
  use openarray
  use variables
  use config
  implicit none
  type(array):: tmpu,tmph,tmp
  integer :: ierr
!  tmph = h; tmph(2, :) = tmph(1, :); tmph = repmat(tmph, 1, 1, kb);
!  tmpu = u; tmpu(im, :, :) = tmpu(imm1, :, :); tmpu(2, :, :) = tmpu(3, :, :);
  tmph=h_3d 
  call set(sub(tmph,  2,':',':'), sub(tmph,    1,':',':'))
  tmpu=u
  call set(sub(tmpu, im,':',':'), sub(tmpu, imm1,':',':'))
  call set(sub(tmpu, 2 ,':',':'), sub(tmpu, 3   ,':',':'))  
 
!  tmp = sqrt(tmph/hmax) .* AYF(AYB(tmpu)) + (1.0 - sqrt(tmph/hmax)) .* AYF(AYB(u));
!  uf(im, 2:jmm1, 1:kbm1) = tmp(im, 2:jmm1, 1:kbm1);
!    uf(2, 2:jmm1, 1:kbm1) = tmp(2, 2:jmm1, 1:kbm1);
!  uf(1, 2:jmm1, 1:kbm1) = tmp(2, 2:jmm1, 1:kbm1);
!    uf(2:imm1, jm, 1:kbm1) = 0.e0;
!    uf(2:imm1, 1, 1:kbm1) = 0.e0;
!    uf = uf .* dum_3d;    

  tmp = sqrt(tmph/hmax)*AYF(AYB(tmpu)) + (1.d0-sqrt(tmph/hmax))*AYF(AYB(u))
  call set(sub(uf, im,r(2,jmm1),r(1,kbm1)), sub(tmp, im,r(2,jmm1),r(1,kbm1)) )
  call set(sub(uf, 2 ,r(2,jmm1),r(1,kbm1)), sub(tmp, 2 ,r(2,jmm1),r(1,kbm1)) )
  call set(sub(uf, 1 ,r(2,jmm1),r(1,kbm1)), sub(tmp, 2 ,r(2,jmm1),r(1,kbm1)) )
  call set(sub(uf, r(2,imm1), jm ,r(1,kbm1)), 0.d0 )
  call set(sub(uf, r(2,imm1), 1  ,r(1,kbm1)), 0.d0 )

  uf = uf * dum_3d

  call destroy(tmpu,ierr); call destroy(tmph,ierr)
  call destroy(tmp,ierr)

end subroutine
