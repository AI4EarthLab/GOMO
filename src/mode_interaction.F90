
subroutine mode_interaction()
  use openarray
  use config
  use variables
  implicit none

  adx2d = sum(advx*dz,  3);
  ady2d = sum(advy*dz,  3);
  drx2d = sum(drhox*dz, 3);
  dry2d = sum(drhoy*dz, 3);
  aam2d = sum(aam*dz,   3);
  
  !call disp(aam*dz, 'dz*aam = ')
  ! call disp(advx, 'advx = ')
  !stop
  ! call disp(advy, 'advy = ')

  ![advua,advva]   = advave(aam2d,uab,vab,ua,va,d);
  call advave()

  adx2d=adx2d-advua;
  ady2d=ady2d-advva;

  egf=el * ispi;
  utf=ua * 2.0 * AXB(d) * isp2i;
  vtf=va * 2.0 * AYB(d) * isp2i;
end subroutine


