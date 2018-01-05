#include "common.h"
subroutine bcond6()
  use openarray
  use variables
  use config
  implicit none
  type(array):: tmpu,tmpv,tmp1,ewflag,snflag
  real(kind=8) :: ewflag_im, ewflag_1, snflag_jm, snflag_1
  integer :: ierr

  ! tmpu = u;               tmpu(1, :, :) = tmpu(2, :, :);
  ! tmpv = v;               tmpv(:, 1, :) = tmpv(:, 2, :);
  ! ewflag=any(fsm,2);      snflag=any(fsm,1);
  tmpu=u;             tmpv=v
  call set(sub(tmpu,  1   ,':',':'), sub(tmpu,  2 ,':',':'))
  call set(sub(tmpv,  ':',1   ,':'), sub(tmpv,  ':',2 ,':'))
  ewflag=sum(fsm,2);  snflag=sum(fsm,1)

  !------------------EAST   (uf stands for q2f; vf stands for q2lf)
  if(ewflag_im  /= 0) then
     !    tmp1 = q2 - dti * ((.5e0 * (u - abs(u))) .* (small - q2) ./ AXB(dx_3d) + (.5e0 * (u + abs(u))) .* DXB(q2));
     !    uf(im, :, :) = tmp1(im, :, :);
     tmp1 = q2-dti*((0.5d0*(u-abs(u)))*(small - q2)/AXB(dx_3d)&
          +(0.5d0*(u + abs(u)))* DXB(q2))
     call set(sub(q2f, im,':',':' ), sub(tmp1,im,':',':' ))

     !tmp1 = q2l - dti * ((.5e0 * (u - abs(u)))
     !     .* (small - q2l) ./ AXB(dx_3d) 
     !     + (.5e0 * (u + abs(u))) .* DXB(q2l));
     !vf(im, :, :) = tmp1(im, :, :);
     tmp1 = q2l-dti*((.5d0*(u-abs(u)))*(small &
          - q2l)/AXB(dx_3d)+(0.5d0*(u + abs(u)))* DXB(q2l))
     call set(sub(q2lf,im,':',':' ), sub(tmp1,  im,':',':' ))
  endif
  !------------------WEST 
  if(ewflag_1  /= 0) then
     !    tmp1 = q2 - dti * ((.5e0 * (tmpu + abs(tmpu))) .* (q2 - small) ./ AXF(dx_3d) + (.5e0 * (tmpu - abs(tmpu))) .* DXF(q2));
     !  uf(1, :, :) = tmp1(1, :, :);
     tmp1 = q2 -dti*((0.5d0*(tmpu &
          + abs(tmpu)))*(q2 - small)/AXF(dx_3d) &
          +(0.5d0*(tmpu - abs(tmpu)))* DXF(q2) )
     call set(sub(q2f,1  ,':',':'), sub(tmp1,  1,':',':'))
     !    tmp1 = q2l - dti * ((.5e0 * (tmpu + abs(tmpu)))
     !     .* (q2l - small) ./ AXF(dx_3d)
     !     + (.5e0 * (tmpu - abs(tmpu))) .* DXF(q2l));
     ! vf(1, :, :) = tmp1(1, :, :);    
     tmp1 = q2l-dti*((0.5d0*(tmpu+abs(tmpu) )) &
          *(q2l-small)/AXF(dx_3d)+(0.5d0*(tmpu-abs(tmpu)))* DXF(q2l))
     call set(sub(q2lf,1,':',':' ), sub(tmp1,  1,':',':' ))    
  endif
  !------------------NORTH
  if(snflag_jm  /= 0) then
     !    tmp1 = q2 - dti * ((.5e0 * (v - abs(v))) .* (small - q2) ./ AYB(dy_3d) + (.5e0 * (v + abs(v))) .* DYB(q2));
     !  uf(:, jm, :) = tmp1(:, jm, :);
     tmp1=q2-dti*((0.5d0 * (v - abs(v)))* (small &
          - q2)/ AYB(dy_3d) + (0.5d0 * (v + abs(v)))* DYB(q2))
     call set(sub(q2f,':', jm ,':'), sub(tmp1,':', jm   ,':'))     

     !    tmp1 = q2l - dti * ((.5e0 * (v - abs(v))) .* (small - q2l) ./ AYB(dy_3d) + (.5e0 * (v + abs(v))) .* DYB(q2l));
     !  vf(:, jm, :) = tmp1(:, jm, :);
     tmp1=q2l-dti * ((0.5d0*(v - abs(v) ) )* (small - q2l)/ AYB(dy_3d) +(0.5d0*(v + abs(v) ))* DYB(q2l))
     call set(sub(q2lf,':', jm ,':'), sub(tmp1,':', jm   ,':') ) 
  endif
  !------------------SOUTH 
  if(snflag_1 /= 0) then
     !  tmp1 = q2 - dti * ((.5e0 * (tmpv + abs(tmpv))) .* (q2 - small) ./ AYF(dy_3d) + (.5e0 * (tmpv - abs(tmpv))) .* DYF(q2));
     !  uf(:, 1, :) = tmp1(:, 1, :);
     tmp1=q2-dti*((0.5d0*(tmpv + abs(tmpv))) &
          * (q2 - small)/ AYF(dy_3d) + (0.5d0*(tmpv - abs(tmpv) ))* DYF(q2))
     call set(sub(q2f,':', 1 ,':'), sub(tmp1,':', 1 ,':') )   

     !  tmp1 = q2l - dti * ((.5e0 * (tmpv + abs(tmpv))) .* (q2l - small) ./ AYF(dy_3d) + (.5e0 * (tmpv - abs(tmpv))) .* DYF(q2l));
     !  vf(:, 1, :) = tmp1(:, 1, :);  
     tmp1 =q2l-dti* ((0.5d0*(tmpv + abs(tmpv))) &
          *(q2l - small)/AYF(dy_3d) &
          +(0.5d0*(tmpv - abs(tmpv)))* DYF(q2l))
     call set(sub(q2lf,':', 1 ,':'), sub(tmp1,':', 1 ,':') )   
  endif

  !  uf = uf .* fsm_3d + 1.e-10;
  !  vf = vf .* fsm_3d + 1.e-10;
  q2f  = q2f  *fsm_3d + 1.d-10
  q2lf = q2lf *fsm_3d + 1.d-10

  !tmpu,tmpv,tmp1,ewflag,snflag
  call destroy(tmpu,ierr);  call destroy(tmpv,ierr); call destroy(tmp1,ierr)
  call destroy(ewflag,ierr);call destroy(snflag,ierr)
end subroutine bcond6
