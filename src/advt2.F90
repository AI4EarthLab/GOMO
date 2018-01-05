#include "common.h"
subroutine advt2(ff,f,fb,fclim)
  use variables
  use config
  use openarray  
  implicit none
  integer:: itera, ierr
  type(array) :: ff,f,fb,fclim
  type(array) :: xmassflux,ymassflux,eta_3d
  type(array) :: zwflux,fbmem,tmp1,tmpim

  xmassflux = axbdt_3d *u;
  ymassflux = aybdt_3d *v;
  eta_3d = etb_3d;
  zwflux = w;
  fbmem = fb;

  call set(sub(fb,':',':',kb) , sub(fb,':',':',kbm1)) 
  call set(tmp1 ,  sub(ff, 1,':',':'))
  call set(tmpim,  sub(ff,im,':',':'))

  do itera=1,nitera

     xflux=(xmassflux>0.e0)*xmassflux*shift(fbmem,1,1) &
          + (xmassflux<0.e0)*xmassflux* fbmem

     yflux=(ymassflux>0.e0)*ymassflux*shift(fbmem,1,2) &
          + (ymassflux<0.e0)*ymassflux* fbmem

     zflux=(zwflux   <0.e0)*zwflux   *shift(fbmem,1,3) &
          + (zwflux   >0.e0)*zwflux   * fbmem

     call set(sub(zflux,':',':',1), 0.e0)

     if(itera==1) &
          call set(sub(zflux,':',':',1), &
          sub(w,':',':',1)*sub(f,':',':',1)) 

     call set(sub(zflux,':',':',kb), 0.e0)            

     ff=(fbmem*(h_3d+eta_3d)-dti2*(DXF(xflux) &
          + DYF(yflux) - DZF(zflux)))/(h_3d+etf_3d)

     call set(sub(ff,1,':',':'), tmp1);
     call set(sub(ff,im,':',':'), tmpim);

     call smol_adif(xmassflux,ymassflux,zwflux,ff)

     !eta=etf;
     fbmem=ff;
  enddo

  fb=fb-fclim;
  xmassflux=AXB(aam);
  ymassflux=AYB(aam);
  xflux=-xmassflux*AXB(h_3d)*tprni*DXB(fb)*dum_3d
  yflux=-ymassflux*AYB(h_3d)*tprni*DYB(fb)*dvm_3d
  ff=ff-dti2* (DXF(xflux)+ DYF(yflux))/(h_3d+etf_3d)

  call set(sub(ff,':',':',kb),0.e0)
  call set(sub(ff,im,':',':'),tmpim)

  call destroy(xmassflux, ierr); call destroy(ymassflux,ierr);
  call destroy(eta_3d, ierr);    call destroy(zwflux,ierr);
  call destroy(fbmem, ierr);     call destroy(tmp1,ierr);
  call destroy(tmpim, ierr);     
end subroutine
