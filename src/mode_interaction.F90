
subroutine mode_interaction()
  use openarray
  use config
  use variables
  implicit none

  ! print *, "start0 of mode_interaction"
  ! call disp(advx, "-------advx------")
  ! call disp(drhox, "-------drhox------")
  ! call disp(dz, "-------dz------")
  adx2d = sum(advx*dz,  3);
  ady2d = sum(advy*dz,  3);
  drx2d = sum(drhox*dz, 3);
  dry2d = sum(drhoy*dz, 3);
  aam2d = sum(aam*dz,   3);
  ! print *, "start2 of mode_interaction"
  
  !call disp(aam*dz, 'dz*aam = ')
  ! call disp(advx, 'advx = ')
  !stop
  ! call disp(advy, 'advy = ')

  ![advua,advva]   = advave(aam2d,uab,vab,ua,va,d);
  call advave()
  ! print *, "start3 of mode_interaction"
  ! call disp(adx2d, "-------adx2d---------")
  ! call disp(ady2d, "-------ady2d---------")
  ! call disp(advua, "-------advua---------")
  ! call disp(advva, "-------advva---------")

  adx2d=adx2d-advua;
  ady2d=ady2d-advva;

  egf=el * ispi;
  utf=ua * 2.0 * AXB(d) * isp2i;
  vtf=va * 2.0 * AYB(d) * isp2i;
  ! print *, "start4 of mode_interaction"
end subroutine


