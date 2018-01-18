
subroutine advave()
  use openarray
  use variables
  use config
  implicit none
  integer:: ierr

  type(array) tmp

  print *, "in advave()"

  tps = AYB(AXB(d)) * AXB(AYB(aam2d)) *  (DYB(uab)  + DXB(vab));
  call disp(tps, "tps = ")

  ! split into two equations in order to keep same as matlab
  ! advua = DXB(AXF( AXB(d) * ua ) * AXF(ua) &
  !      - 2.d0 * d * aam2d * DXF(uab))  &
  !      + DYF(AXB(AYB(d) * va ) * AYB(ua)-tps);
  
  tmp = AXF( AXB(d) * ua ) * AXF(ua) &
       - 2.d0 * d * aam2d * DXF(uab);
  call set(sub(tmp, 1,':',':'),0.d0)
  call set(sub(tmp, im,':',':'),0.d0)
  call disp(tmp, "tmp_advua = ")
  advua = DXB(tmp)  &
       + DYF(AXB(AYB(d) * va ) * AYB(ua)-tps);
  call set(sub(advua, 1,':',':'), 0.d0)
  call set(sub(advua, ':',1,':'), 0.d0)
  call set(sub(advua, ':',jm,':'), 0.d0)

  ! split into two equations in order to keep same as matlab
  ! advva = DXF(AYB( AXB(d) * ua ) * AXB(va)-tps) &
  !      +DYB(AYF( AYB(d) * va ) * AYF(va) &
  !      - 2.d0 * d * aam2d * DYF(vab));
  
  tmp = AYF( AYB(d) * va ) * AYF(va) &
       - 2.d0 * d * aam2d * DYF(vab);
  call set(sub(tmp, ':',1,':'), 0.d0)
  call set(sub(tmp, ':',jm,':'), 0.d0)
  call disp(tmp, "tmp_advva = ")
  advva = DXF(AYB( AXB(d) * ua ) * AXB(va)-tps) &
       +DYB(tmp);
  call set(sub(advva,':',1,':'), 0.d0)
  call set(sub(advva, 1,':',':'),0.d0)
  call set(sub(advva, im,':',':'),0.d0)

  call disp(advua, "advua = ");
  call disp(advva, "advva = ");

end subroutine

