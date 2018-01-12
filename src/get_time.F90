
#include "common.h"

!  function [time,ramp] = get_time(iint,time0,lramp,dti,cor)
subroutine get_time(iint)
  use openarray
  use variables
  use config
  implicit none
  real(kind=8) :: cor_im2_jm2
  integer :: iint
  
  time=dti*iint*1.0/86400.0+time0;
  
  if(lramp) then
     cor_im2_jm2 = sub(cor, int(im/2),int(jm/2), 1)
     ramp = time/((2.0*pi)/abs(cor_im2_jm2)/86400.0);

     if(ramp > 1.0) then
        ramp=1.0;
     endif
  else
     ramp=1.0;
  end if
end subroutine
