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
    print *, "in external_ua"
    call disp(advua, "advua = ")
  endif

  ! call disp(uaf, "before uaf = ")
  ! call disp(h, "h = ")
  ! call disp(elb, "elb = ")
  ! call disp(uab, "uab = ")
  ! print *, "dte = ", dte
  ! call disp(adx2d, "adx2d = ")
  ! call disp(advua, "advua = ")
  ! call disp(cor, "cor = ")
  ! call disp(d, "d = ")
  ! call disp(va, "va = ")
  ! print *, grav, "grav = "
  ! print *, alpha, "alpha = "
  ! call disp(el, "el = ")
  ! call disp(elb, "elb = ")
  ! call disp(elf, "elf = ")
  ! call disp(e_atmos, "e_atmos = ")
  ! call disp(drx2d, "drx2d = ")
  ! call disp(wusurf, "wusurf = ")
  ! call disp(wubot, "wubot = ")

  uaf=(AXB(h+elb)*uab - 2.0*dte*(adx2d+advua &
       - AXB(cor*d*AYF(va))       &
       + grav*AXB(d)*( (1.0-2.0*alpha)*DXB(el) &
       + alpha*(DXB(elb)+DXB(elf))+DXB(e_atmos) )  &
       + drx2d+(wusurf-wubot)) ) / AXB(h+elf)

  call bcond2_ua()
  call disp(uaf, "uaf = ")

end subroutine external_ua
