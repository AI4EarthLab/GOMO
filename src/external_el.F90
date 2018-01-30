subroutine external_el()
  use openarray
  use variables
  use config
  implicit none

!    elf= elb-dte2.*((DXF(AXB(d).*ua)+DYF(AYB(d).*va))-vfluxf);
!    elf = bcond1(elf);
    elf= elb-2.0*dte*((DXF(AXB(d)*ua)+DYF(AYB(d)*va))-vfluxf)
    
    call bcond1()
    !call disp(elf, "elf = ")

end subroutine external_el
