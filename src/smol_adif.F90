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
  
  value_min = 1.0d-9;
  epsilon   = 1.0d-14;
  
  ff = ff * fsm;
  xmassflux = xmassflux * AXB(dy)
  ymassflux = ymassflux * AYB(dx)

  !%     Recalculate mass fluxes with antidiffusion velocity:
  call set(xmassflux, 0.d0, &
       mask=(ff<value_min) .or. (shift(ff,1,1)<value_min))
  
  flagno = ((ff>=value_min) .and. (shift(ff,1,1)>=value_min))
  
  udx = abs(xmassflux) * flagno
  
  u2dt = dti2 * xmassflux * xmassflux * flagno/(aru*AXB(dt))
  
  xmassflux=(udx-u2dt)*DXB(ff)*AXB(dx)/(2.d0*AXB(ff)+epsilon)*sw

  call set(xmassflux, 0.d0, mask=abs(udx)<abs(u2dt))
  
  call set(ymassflux, 0.d0, &
       mask=(ff<value_min).or.(shift(ff,1,2)<value_min))
  
  flagno = (ff>=value_min) .and. (shift(ff,1,2)>=value_min)
  
  vdy = abs(ymassflux) * flagno
  
  v2dt=dti2*ymassflux*ymassflux/(arv* AYB(dt) ) *flagno
  
  ymassflux=(vdy-v2dt)*DYB(ff)*AYB(dy)/(2*AYB(ff)+epsilon)*sw
  
  call set(ymassflux, 0.d0, abs(vdy)<abs(v2dt))
  
  call set(zwflux, 0.d0, &
       mask=(ff<value_min) .or. (shift(ff,1,3)<value_min))
  
  flagno=((ff>=value_min) .and. (shift(ff,1,3)>=value_min))
  
  wdz  = abs(zwflux) * flagno
  w2dt = dti2*zwflux*zwflux* flagno/(circshift(dzz,1,3)*dt) !change shift parameters.
  zwflux = -(wdz-w2dt)*DZB(ff)* AZB(dz)/(2.d0*AZB(ff)+epsilon)*sw

  call set(zwflux, 0.d0, mask=abs(wdz) < abs(w2dt))
  
  xmassflux=xmassflux/AXB(dy)
  ymassflux=ymassflux/AYB(dx)

end subroutine smol_adif
