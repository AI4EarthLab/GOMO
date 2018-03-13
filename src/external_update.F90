subroutine external_update(iext, vamax, imax, jmax)
  use openarray
  use config
  use variables
  implicit none
  integer :: iext, pos(3)
  integer, intent(out) :: imax, jmax
  real(kind=8), intent(out) :: vamax
 
  if(iext==(isplit-2)) then
     etf=0.25*smoth*elf 
  else if(iext==(isplit-1)) then
     etf=etf+0.5*(1.0-0.5*smoth)*elf
  else if(iext==isplit) then
     etf=(etf+0.5*elf)*fsm
  endif


  vamax = abs_max(vaf)

  call set(pos, abs_max_at(vaf))
  imax = pos(1)
  jmax = pos(2)
  
  if(vamax<=vmaxl) then            
     uab=ua+0.5*smoth*(uab-2.0*ua+uaf)
     ua=uaf
     vab=va+0.5*smoth*(vab-2.0*va+vaf)
     va=vaf
     elb=el+0.5*smoth*(elb-2.0*el+elf)
     el=elf
     d=h+el
     if(iext/=isplit) then
        egf=egf+el*ispi
        utf=utf+2.0* ua * AXB(d) * isp2i
        vtf=vtf+2.0* va * AYB(d) * isp2i
     endif
  endif

end subroutine
