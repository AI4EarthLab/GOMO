#include "common.h"
subroutine bcond6()
  use openarray
  use variables
  use config
  implicit none
  type(array):: tmpu,tmpv,tmp1,ewflag,snflag
  real(kind=8) :: ewflag_im, ewflag_1, snflag_jm, snflag_1
  integer :: ierr

  tmpu=u;             tmpv=v
  call set(sub(tmpu,  1   ,':',':'), sub(tmpu,  2 ,':',':'))
  call set(sub(tmpv,  ':',1   ,':'), sub(tmpv,  ':',2 ,':'))
  ! ewflag=sum(fsm,2);  snflag=sum(fsm,1)

  ! ewflag_im = sub(ewflag, im,1,1)
  ! ewflag_1  = sub(ewflag, 1,1,1)
  ! snflag_jm = sub(snflag, 1,jm,1)
  ! snflag_1  = sub(snflag, 1,1,1)

  ! call disp(ewflag, "ewflag = ")
  ! call disp(snflag, "snflag = ")

  !------------------EAST   (uf stands for q2f; vf stands for q2lf)
  !if(ewflag_im  /= 0) then
    tmp1 = q2-dti*((0.5d0*(u-abs(u)))*(small - q2)/AXB(dx)&
         +(0.5d0*(u + abs(u)))* DXB(q2))
    call set(sub(q2f, im,':',':' ), sub(tmp1,im,':',':' ))

    tmp1 = q2l-dti*((.5d0*(u-abs(u)))*(small &
         - q2l)/AXB(dx)+(0.5d0*(u + abs(u)))* DXB(q2l))
    call set(sub(q2lf,im,':',':' ), sub(tmp1,  im,':',':' ))
  !endif

  !------------------WEST 
  !if(ewflag_1  /= 0) then
    tmp1 = q2 -dti*((0.5d0*(tmpu &
         + abs(tmpu)))*(q2 - small)/AXF(dx) &
         +(0.5d0*(tmpu - abs(tmpu)))* DXF(q2) )
    call set(sub(q2f,1  ,':',':'), sub(tmp1,  1,':',':'))

    tmp1 = q2l-dti*((0.5d0*(tmpu+abs(tmpu) )) &
         *(q2l-small)/AXF(dx)+(0.5d0*(tmpu-abs(tmpu)))* DXF(q2l))
    call set(sub(q2lf,1,':',':' ), sub(tmp1,  1,':',':' ))    
  !endif

  !------------------NORTH
  !if(snflag_jm  /= 0) then
    tmp1=q2-dti*((0.5d0 * (v - abs(v)))* (small &
         - q2)/ AYB(dy) + (0.5d0 * (v + abs(v)))* DYB(q2))
    call set(sub(q2f,':', jm ,':'), sub(tmp1,':', jm   ,':'))     
    
    tmp1=q2l-dti * ((0.5d0*(v - abs(v) ) )* (small - q2l)/ AYB(dy) +(0.5d0*(v + abs(v) ))* DYB(q2l))
    call set(sub(q2lf,':', jm ,':'), sub(tmp1,':', jm   ,':') ) 
  !endif

  !------------------SOUTH 
  !if(snflag_1 /= 0) then
    tmp1=q2-dti*((0.5d0*(tmpv + abs(tmpv))) &
         * (q2 - small)/ AYF(dy) + (0.5d0*(tmpv - abs(tmpv) ))* DYF(q2))
    call set(sub(q2f,':', 1 ,':'), sub(tmp1,':', 1 ,':') )   

    tmp1 =q2l-dti* ((0.5d0*(tmpv + abs(tmpv))) &
         *(q2l - small)/AYF(dy) &
         +(0.5d0*(tmpv - abs(tmpv)))* DYF(q2l))
    call set(sub(q2lf,':', 1 ,':'), sub(tmp1,':', 1 ,':') )   
  !endif

  q2f  = q2f  *fsm + 1.d-10
  q2lf = q2lf *fsm + 1.d-10

end subroutine bcond6
