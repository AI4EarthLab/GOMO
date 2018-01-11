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

  ! call disp(h, "-----h-----")
  ! call disp(elb, "-----elb-----")
  ! call disp(uab, "-----uab-----")
  ! call disp(adx2d, "-----adx2d-----")
  ! call disp(advua, "-----advua-----")
  ! call disp(cor, "-----cor-----")
  ! call disp(d, "-----d-----")
  ! call disp(va, "-----va-----")
  ! call disp(el, "-----el-----")
  ! call disp(elf, "-----elf-----")
  ! call disp(e_atmos, "-----e_atmos-----")
  ! call disp(drx2d, "-----drx2d-----")
  ! call disp(wusurf, "-----wusurf-----")
  ! call disp(wubot, "-----wubot-----")
  ! 
  ! call disp(uaf, "-----uaf-----")

  uaf=(AXB(h+elb)*uab -2.0* dte* (adx2d + advua &
       - AXB(cor* d* AYF(va))       &
       + grav*AXB(d)*( (1.0-2.0*alpha)*DXB(el) &
       +alpha*(DXB(elb)+DXB(elf))+DXB(e_atmos) )  &
       + drx2d+(wusurf-wubot) ) )/ AXB(h+elf)

  ! call disp(uaf, "-----uaf-----")
  ! uaf=-2.0* dte* (adx2d + advua - AXB(cor* d* AYF(va))       &
  !      + grav*AXB(d)*( (1.0-2.0*alpha)*DXB(el) &
  !      +alpha*(DXB(elb)+DXB(elf))+DXB(e_atmos) )  &
  !      + drx2d) 
  
  !call disp(h+elf, 'h+elf = ')
  ! print *, "in external_ua"
  call bcond2_ua()

end subroutine external_ua
