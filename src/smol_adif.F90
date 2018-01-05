subroutine smol_adif(xmassflux,ymassflux,zwflux,ff)
  use openarray
  use config
  use variables
  implicit none
  real*8:: value_min,epsilon
  type(array) :: ff
  type(array) :: xmassflux,ymassflux,zwflux
  type(array) :: flagno,udx,u2dt,vdy,v2dt,wdz,w2dt
  integer :: ierr
  
  !	value_min=1.0e-9;   epsilon=1.0e-14;
  !	ff=ff*fsm_3d;      dt_3d=repmat(dt,1,1,kb);
  value_min = 1.0e-9;
  epsilon   = 1.0e-14;
  
  ff = ff * fsm_3d;
  xmassflux = xmassflux * AXB(dy_3d)
  ymassflux = ymassflux * AYB(dx_3d)

  !%     Recalculate mass fluxes with antidiffusion velocity:
  !      xmassflux(double(or(ff<value_min, shift(ff,1,1)<value_min)))=0.e0;
  !      flagno= and(ff>=value_min, shift(ff,1,1)>=value_min);
  !      udx=abs(xmassflux).*flagno;
  !      u2dt=DIVISION(dti2.*xmassflux.^2 .*flagno,(aru_3d.*AXB(dt_3d)));
  !      xmassflux=(udx-u2dt).*(DXB(ff).*AXB(dx_3d)./(2*AXB(ff)+epsilon)).*sw;
  !      xmassflux(double(lt(abs(udx),abs(u2dt))))=0.e0;
  !xmassflux((ff<value_min) .or. (shift(ff,1,1)<value_min))=0.e0
  call set(xmassflux, 0.0, &
       filter=(ff<value_min) .or. (shift(ff,1,1)<value_min))
  
  flagno = ((ff>=value_min) .and. (shift(ff,1,1)>=value_min))
  
  udx = abs(xmassflux) * flagno
  
  u2dt = dti2 * xmassflux * xmassflux * flagno/(aru_3d*axbdt_3d )
  
  xmassflux=(udx-u2dt)*DXB(ff)*AXB(dx_3d)/(2.0*AXB(ff)+epsilon)*sw
  !xmassflux(abs(udx)<abs(u2dt))=0.e0 
  call set(xmassflux, 0.0, filter=abs(udx)<abs(u2dt))
  
  !ymassflux((ff<value_min).or.(shift(ff,1,2)<value_min))=0.e0
  call set(ymassflux, 0.0, &
       filter=(ff<value_min).or.(shift(ff,1,2)<value_min))
  
  flagno = (ff>=value_min) .and. (shift(ff,1,2)>=value_min)
  
  vdy = abs(ymassflux) * flagno
  
  v2dt=dti2*ymassflux*ymassflux/(arv_3d* aybdt_3d ) *flagno
  
  ymassflux=(vdy-v2dt)*DYB(ff)*AYB(dy_3d)/(2*AYB(ff)+epsilon)*sw
  
  !ymassflux(abs(vdy)<abs(v2dt))=0.e0
  call set(ymassflux, 0.0, abs(vdy)<abs(v2dt))
  
  !zwflux((ff<value_min) .or. (shift(ff,1,3)<value_min))=0.e0
  call set(zwflux, 0.0, &
       filter=(ff<value_min) .or. (shift(ff,1,3)<value_min))
  
  flagno=((ff>=value_min) .and. (shift(ff,1,3)>=value_min))
  
  wdz  = abs(zwflux) * flagno
  w2dt = dti2*zwflux*zwflux* flagno/(circshift(dzz_3d,1,3)*dt_3d) !change shift parameters.
  zwflux = -(wdz-w2dt)*DZB(ff)* AZB(dz_3d)/(2.0*AZB(ff)+epsilon)*sw

  !zwflux(abs(wdz)<abs(w2dt))=0.e0
  call set(zwflux, 0.0, filter=abs(wdz) < abs(w2dt))
  
  xmassflux=xmassflux/AXB(dy_3d)
  ymassflux=ymassflux/AYB(dx_3d)

  !flagno,udx,u2dt,vdy,v2dt,wdz,w2dt
  call destroy(flagno, ierr); call destroy(udx,ierr);
  call destroy(u2dt, ierr);   call destroy(vdy,ierr);
  call destroy(v2dt, ierr);   call destroy(wdz,ierr);
  call destroy(w2dt, ierr);
end subroutine smol_adif
