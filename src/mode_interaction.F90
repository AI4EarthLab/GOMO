
subroutine mode_interaction()
  use openarray
  use config
  use variables
  implicit none

  print *, "in mode_interaction"
  adx2d = sum(advx*dz,  3);
  ady2d = sum(advy*dz,  3);
  drx2d = sum(drhox*dz, 3);
  dry2d = sum(drhoy*dz, 3);
  aam2d = sum(aam*dz,   3);

  !call disp(adx2d, "adx2d = ")
  !call disp(ady2d, "ady2d = ")
  !call disp(drx2d, "drx2d = ")
  !call disp(dry2d, "dry2d = ")
  !call disp(aam2d, "aam2d = ")
  
  call advave()

  adx2d=adx2d-advua;
  ady2d=ady2d-advva;

  !call disp(adx2d, "adx2d = ")
  !call disp(ady2d, "ady2d = ")

  egf=el * ispi;
  utf=ua * 2.0 * AXB(d) * isp2i;
  vtf=va * 2.0 * AYB(d) * isp2i;

  !call disp(utf, "utf = ")
  !call disp(vtf, "vtf = ")
end subroutine


