#include "common.h"
subroutine advt2(ff,f,fb,fclim)
  use variables
  use config
  use openarray  
  implicit none
  integer:: itera, ierr
  type(array) :: ff,f,fb,fclim
  type(array) :: xmassflux,ymassflux,eta
  type(array) :: zwflux,fbmem,tmp1,tmpim

  xmassflux = AXB(dt)*u;
  ymassflux = AYB(dt)*v;
  eta = etb;
  zwflux = w;
  fbmem = fb;

  call set(sub(fb,':',':',kb) , sub(fb,':',':',kbm1)) 
  tmp1 = sub(ff, 1,':',':')
  tmpim = sub(ff,im,':',':')

  do itera=1,nitera

     xflux=(xmassflux>0.d0)*xmassflux*shift(fbmem,1,1) &
          + (xmassflux<0.d0)*xmassflux*fbmem

     yflux=(ymassflux>0.d0)*ymassflux*shift(fbmem,1,2) &
          + (ymassflux<0.d0)*ymassflux*fbmem

     zflux=(zwflux<0.d0)*zwflux*shift(fbmem,1,3) &
          + (zwflux>0.d0)*zwflux*fbmem

     call set(sub(zflux,':',':',1), 0.d0)

     if(itera==1) &
          call set(sub(zflux,':',':',1), &
          sub(w,':',':',1)*sub(f,':',':',1)) 

     call set(sub(zflux,':',':',kb), 0.d0)            

     ff=(fbmem*(h+eta)-dti2*(DXF(xflux) &
          + DYF(yflux) - DZF(zflux)))/(h+etf)

     call set(sub(ff,1,':',':'), tmp1);
     call set(sub(ff,im,':',':'), tmpim);

     call smol_adif(xmassflux,ymassflux,zwflux,ff)

     eta=etf;
     fbmem=ff;
  enddo

  fb=fb-fclim;
  xmassflux=AXB(aam);
  ymassflux=AYB(aam);
  xflux=-xmassflux*AXB(h)*tprni*DXB(fb)*dum
  yflux=-ymassflux*AYB(h)*tprni*DYB(fb)*dvm
  ff=ff-dti2* (DXF(xflux)+ DYF(yflux))/(h+etf)

  call set(sub(ff,':',':',kb),0.d0)
  call set(sub(ff,im,':',':'),tmpim)

end subroutine
