subroutine bcond4(ff,f,fb,fbe,fbw,fbn,fbs)
  use openarray
  use config
  !  use namelist_constant
  use variables
  implicit none
  integer::ierr
  type(array):: tmptb,tmpu,tmpv,tmp1,tmp2,ewflag,snflag
  type(array), intent(inout):: ff,f,fb,fbe,fbw,fbn,fbs
  real*8::eflag, wflag, sflag, nflag

  tmptb = mat_zeros  
  
  call set(sub(tmptb,  im,':',':'), fbe)
  call set(sub(tmptb,  1 ,':',':'), fbw)
  call set(sub(tmptb,  ':',jm,':'), fbn)
  call set(sub(tmptb,  ':',1 ,':'), fbs)

  tmpu=u; tmpv=v
  call set(sub(tmpu,  imm1,':',':'), sub(tmpu,  im,':',':'))
  call set(sub(tmpu,  1   ,':',':'), sub(tmpu,  2 ,':',':'))
  call set(sub(tmpv,  ':',jmm1,':'), sub(tmpv,  ':',jm,':'))
  call set(sub(tmpv,  ':',1   ,':'), sub(tmpv,  ':',2 ,':'))

  ! ewflag=sum(fsm,2);    snflag=sum(fsm,1)

  ! wflag = sub(ewflag, 1)
  ! eflag = sub(ewflag, im)
  ! sflag = sub(snflag, 1,1,1)  
  ! nflag = sub(snflag, 1,jm,1)

  ! !------------------EAST
  ! if(eflag .gt. 0.d0) then
     tmp1 = f-dti*((0.5d0*(u-abs(u)))* &
          (tmptb - f)/AXB(dx) + &
          (0.5d0*(u + abs(u)))* DXB(f))
    
     call set(sub(ff, im,':',[1,kbm1]), &
          sub(tmp1, im,':',[1,kbm1]))
     
     call set(sub(tmp1,imm1,':',':'), &
          sub(tmp1,  im,':',':'))

     tmp2 = tmp1+0.5d0*dti*(tmpu + abs(tmpu))/tmpu *&
          AZF(w)/dt*DZF(AZB(f))*dz/AZB(dzz)
     
     call set(sub(ff,im,':',[2,kbm2]), &
          sub(tmp2,  imm1,':',[2,kbm2]))
  ! endif

  !print*, __FILE__, __LINE__
  
  !------------------WEST 
  ! if(wflag .gt. 0.d0 ) then
     tmp1 = f -dti*((.5d0*(tmpu + abs(tmpu)))* &
          (f - tmptb)/AXF(dx)+ &
          (.5d0*(tmpu - abs(tmpu)))* DXF(f)*dx/(AXF(dx)))
     
     call set(sub(ff,  1  ,':',[1,kbm1]), &
          sub(tmp1,  1,':',[1,kbm1]))
     
     call set(sub(tmp1, 2,':',':'), &
          sub(tmp1, 1,':',':'))
     
     tmp2 =tmp1+0.5d0*dti *(u - abs(u))/u * &
          AZF(w)/ dt* DZF(AZB(f))*dz/AZB(dzz)
     
     call set(sub(ff, 1,':',[2,kbm2]), &
          sub(tmp2, 2,':',[2,kbm2]))
  ! endif

  !print*, __FILE__, __LINE__
  
  !------------------NORTH
  ! if(nflag .gt. 0.d0) then
     tmp1 = f - dti*((0.5d0*(v-abs(v)))* &
          (tmptb - f)/AYB(dy) + &
          (0.5d0*(v+abs(v)))*DYB(f))
     
     call set(sub(ff,':', jm,[1,kbm1]), &
          sub(tmp1,':', jm,[1,kbm1]))
     
     call set(sub(tmp1,':', jmm1 ,':'), &
          sub(tmp1,':', jm,':'))
     
     tmp2 = tmp1 + 0.5d0*dti*(tmpv + abs(tmpv))/tmpv * &
          AZF(w)/dt *DZF(AZB(f))*dz/AZB(dzz)
     
     !    ff(:, jm, 2:kbm2) = tmp2(:, jmm1, 2:kbm2);
     call set(sub(ff,':',jm,[2,kbm2] ), &
          sub(tmp2,  ':',jmm1,[2,kbm2]))
  ! endif

  !print*, __FILE__, __LINE__
  
  !------------------SOUTH 
  !if(sflag .gt. 0.d0) then
     tmp1 = f-dti*((0.5d0*(tmpv + abs(tmpv))) * &
          (f - tmptb)/AYF(dy)+ &
          (0.5d0 * (tmpv - abs(tmpv)))*DYF(f))
     
     call set(sub(ff,':', 1,[1,kbm1]), &
          sub(tmp1,':', 1,[1,kbm1]))
     
     call set(sub(tmp1,':', 2 ,':'), &
          sub(tmp1,':', 1 ,':'))
     
     tmp2 = tmp1+0.5d0*dti*(v - abs(v))/v * &
          AZF(w) / dt* DZF(AZB(f))*dz/AZB(dzz)
     
     call set(sub(ff, ':',1,[2,kbm2] ), &
          sub(tmp2, ':',2,[2,kbm2]))
  !endif

  call set(sub(ff, ':',1,[2,kbm2]), 0.d0)
  call set(sub(ff, ':',jm,[2,kbm2]), 0.d0)

  !print*, __FILE__, __LINE__
  
  ff = ff * fsm

end subroutine bcond4
