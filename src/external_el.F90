subroutine external_el()
  use openarray
  use variables
  use config
  implicit none

    elf= elb-2.0*dte*((DXF(AXB(d)*ua)+DYF(AYB(d)*va))-vfluxf)
    call bcond1()

end subroutine external_el
