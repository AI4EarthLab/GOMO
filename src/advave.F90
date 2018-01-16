
subroutine advave()
  use openarray
  use variables
  use config
  implicit none
  integer:: ierr

  tps = AYB(AXB(d)) * AXB(AYB(aam2d)) *  (DYB(uab)  + DXB(vab));

  advua = DXB(AXF( AXB(d) * ua ) * AXF(ua) &
       - 2.d0 * d * aam2d * DXF(uab))  &
       + DYF(AXB(AYB(d) * va ) * AYB(ua)-tps);
  call set(sub(advua,1,':',':'), 0.d0)

  advva = DXF(AYB( AXB(d) * ua ) * AXB(va)-tps) &
       +DYB(AYF( AYB(d) * va ) * AYF(va) &
       - 2.d0 * d * aam2d * DYF(vab));
  call set(sub(advva,':',1,':'), 0.d0)

  !call disp(advua, "advua = ");
  !call disp(advva, "advva = ");

end subroutine

