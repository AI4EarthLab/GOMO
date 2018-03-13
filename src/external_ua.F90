subroutine external_ua(iext)
  use openarray
  use variables
  use config
  implicit none
  integer :: iext

  if(mod(iext,ispadv)==0) then
    advua = DXB( AXF(AXB(d)*ua)* AXF(ua)-2.0*d* aam2d* DXF(uab) ) &
         +DYF(AXB(AYB(d)*va)*AYB(ua)-AYB(AXB(d)) &
         *AXB(AYB(aam2d))*(DYB(uab)+DXB(vab))) 
  endif

  uaf=(AXB(h+elb)*uab - 2.0*dte*(adx2d+advua &
       - AXB(cor*d*AYF(va))       &
       + grav*AXB(d)*( (1.0-2.0*alpha)*DXB(el) &
       + alpha*(DXB(elb)+DXB(elf))+DXB(e_atmos) )  &
       + drx2d+(wusurf-wubot)) ) / AXB(h+elf)

  call bcond2_ua()

end subroutine external_ua
